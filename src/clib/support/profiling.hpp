//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Lightweight, compile-gated profiling harness for the chemistry/cooling
/// solver.
///
/// This module exists to gather *baseline* performance data ahead of the GPU
/// port. It answers two questions that the existing `cxx_omp_example` benchmark
/// cannot, because that benchmark only times whole public-API calls:
///
///   1. Where does time go *inside* `solve_rate_cool` (per kernel)?
///   2. How do cells distribute across solver paths (Gauss-Seidel vs.
///      Newton-Raphson) and how many subcycles do they take? This is the key
///      thread-divergence risk for the GPU port.
///
/// EVERYTHING here is gated behind the `GRACKLE_PROFILE` macro. When it is not
/// defined, every macro below expands to a no-op and there is zero runtime or
/// memory overhead -- so the instrumentation can stay in the source tree
/// permanently and be switched on only when measuring (e.g. compile with
/// `-DGRACKLE_PROFILE`).
///
/// Usage (in instrumented code):
/// @code
///   GRACKLE_PROF_SCOPE(cool1d_multi_g);          // time the enclosing scope
///   GRACKLE_PROF_COUNT(cells_gs, n_gs_cells);    // bump a counter
///   GRACKLE_PROF_HIST_SUBCYCLE(iter);            // record an i-slice's #subcycles
///   GRACKLE_PROF_REPORT();                       // dump the table (once, at end)
/// @endcode
///
//===----------------------------------------------------------------------===//

#ifndef GRACKLE_SUPPORT_PROFILING_HPP
#define GRACKLE_SUPPORT_PROFILING_HPP

#ifndef GRACKLE_PROFILE

// ---------------------------------------------------------------------------
// Disabled build: every hook is a no-op. No state, no headers, no overhead.
// ---------------------------------------------------------------------------
#define GRACKLE_PROF_SCOPE(name)        ((void)0)
#define GRACKLE_PROF_COUNT(name, n)     ((void)0)
#define GRACKLE_PROF_HIST_SUBCYCLE(n)   ((void)0)
#define GRACKLE_PROF_REPORT()           ((void)0)
#define GRACKLE_PROF_REPORT_NOW()       ((void)0)
#define GRACKLE_PROF_RESET()            ((void)0)

#else  // GRACKLE_PROFILE

#include <cstdint>
#include <cstdio>
#include <cstdlib>   // std::getenv
#include <cstring>   // std::strcmp

#ifdef _OPENMP
#include <omp.h>
#else
#include <chrono>
#endif

namespace grackle::impl::prof {

// ---------------------------------------------------------------------------
// The set of timed kernels and counters is declared via X-macro lists, so
// adding a new slot is a single line and the enum / name-table / storage all
// stay in sync automatically.
// ---------------------------------------------------------------------------

#define GRACKLE_PROF_TIMER_LIST                                                \
  X(extended_gas_props)                                                        \
  X(cool1d_multi_g)                                                            \
  X(lookup_cool_rates1d)                                                       \
  X(rate_timestep)                                                             \
  X(set_subcycle_dt)                                                           \
  X(step_rate_gauss_seidel)                                                    \
  X(step_rate_newton_raphson)                                                  \
  /* sub-timers INSIDE step_rate_newton_raphson (per Newton iteration): */     \
  X(nr_residual)         /* the single derivatives() eval (residual F)    */   \
  X(nr_jacobian)         /* finite-difference Jacobian: ~nsp derivs evals */   \
  X(nr_gaussj)           /* the dense nsp x nsp Gauss-Jordan linear solve */   \
  X(make_consistent)                                                           \
  X(solve_rate_cool_total)

#define GRACKLE_PROF_COUNTER_LIST                                              \
  X(cells_gs)            /* cell-subcycles solved with Gauss-Seidel        */  \
  X(cells_nr)            /* cell-subcycles solved with Newton-Raphson      */  \
  X(cells_nr_coevolve)   /* subset of NR that co-evolves internal energy   */  \
  X(newton_iterations)   /* total Newton iterations across all NR cells    */  \
  X(islices)             /* number of i-slices (j,k pairs) processed       */  \
  X(subcycles)           /* total subcycle iterations across all i-slices  */  \
  X(islices_maxed_out)   /* i-slices that hit max_iterations               */

enum TimerId {
#define X(name) name,
  GRACKLE_PROF_TIMER_LIST
#undef X
  N_TIMERS
};

enum CounterId {
#define X(name) name,
  GRACKLE_PROF_COUNTER_LIST
#undef X
  N_COUNTERS
};

inline const char* timer_name(int id) {
  static const char* names[] = {
#define X(name) #name,
    GRACKLE_PROF_TIMER_LIST
#undef X
  };
  return names[id];
}

inline const char* counter_name(int id) {
  static const char* names[] = {
#define X(name) #name,
    GRACKLE_PROF_COUNTER_LIST
#undef X
  };
  return names[id];
}

// ---------------------------------------------------------------------------
// Storage: per-thread accumulators, padded to a cache line to avoid false
// sharing. Reduced across threads only at report time.
// ---------------------------------------------------------------------------

inline constexpr int    MAX_THREADS = 1024;
inline constexpr int    N_HIST_BINS = 24;  // log2 bins for subcycle counts

struct alignas(64) ThreadTimers {
  double   seconds[N_TIMERS] = {};
  uint64_t calls[N_TIMERS]   = {};
};

struct alignas(64) ThreadCounters {
  uint64_t counter[N_COUNTERS] = {};
  uint64_t subcycle_hist[N_HIST_BINS] = {};  // bin = floor(log2(iters))
};

// Defined inline (C++17) so the header needs no .cpp file.
inline ThreadTimers   g_timers[MAX_THREADS];
inline ThreadCounters g_counters[MAX_THREADS];

inline int thread_id() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

inline double now_seconds() {
#ifdef _OPENMP
  return omp_get_wtime();
#else
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
#endif
}

// ---------------------------------------------------------------------------
// Hot-path primitives (all trivially inlinable)
// ---------------------------------------------------------------------------

class ScopedTimer {
 public:
  explicit ScopedTimer(TimerId id)
      : id_(id), tid_(thread_id()), start_(now_seconds()) {}
  ~ScopedTimer() {
    ThreadTimers& s = g_timers[tid_];
    s.seconds[id_] += now_seconds() - start_;
    s.calls[id_]   += 1;
  }
  ScopedTimer(const ScopedTimer&)            = delete;
  ScopedTimer& operator=(const ScopedTimer&) = delete;

 private:
  TimerId id_;
  int     tid_;
  double  start_;
};

inline void add_counter(CounterId id, uint64_t n) {
  g_counters[thread_id()].counter[id] += n;
}

inline void record_subcycle_hist(int iters) {
  if (iters < 1) iters = 1;
  int bin = 0;
  // floor(log2(iters)) without <cmath>
  while ((iters >> (bin + 1)) != 0 && bin + 1 < N_HIST_BINS) ++bin;
  g_counters[thread_id()].subcycle_hist[bin] += 1;
}

inline void reset() {
  for (int t = 0; t < MAX_THREADS; ++t) {
    g_timers[t]   = ThreadTimers{};
    g_counters[t] = ThreadCounters{};
  }
}

// ---------------------------------------------------------------------------
// Reporting (called once, outside any parallel region)
// ---------------------------------------------------------------------------

inline void report() {
  // Pick destination: GRACKLE_PROFILE_OUT=<path> appends to a file, else stderr.
  std::FILE* out = stderr;
  bool close_out = false;
  if (const char* path = std::getenv("GRACKLE_PROFILE_OUT")) {
    if (path[0] != '\0' && std::strcmp(path, "-") != 0) {
      if (std::FILE* f = std::fopen(path, "a")) { out = f; close_out = true; }
    }
  }

  // Reduce across threads.
  double   tsec[N_TIMERS]   = {};
  uint64_t tcalls[N_TIMERS] = {};
  uint64_t cnt[N_COUNTERS]  = {};
  uint64_t hist[N_HIST_BINS] = {};
  int      used_threads = 0;
  for (int t = 0; t < MAX_THREADS; ++t) {
    bool touched = false;
    for (int i = 0; i < N_TIMERS; ++i) {
      tsec[i]   += g_timers[t].seconds[i];
      tcalls[i] += g_timers[t].calls[i];
      if (g_timers[t].calls[i]) touched = true;
    }
    for (int i = 0; i < N_COUNTERS; ++i) {
      cnt[i] += g_counters[t].counter[i];
      if (g_counters[t].counter[i]) touched = true;
    }
    for (int b = 0; b < N_HIST_BINS; ++b) hist[b] += g_counters[t].subcycle_hist[b];
    if (touched) ++used_threads;
  }

  const double total = tsec[solve_rate_cool_total];
  const double denom = (total > 0.0) ? total : 1.0;

  std::fprintf(out,
      "\n================ Grackle solve_rate_cool profile ================\n");
  std::fprintf(out, "threads that did work: %d\n\n", used_threads);

  // --- Per-kernel timing (wall-clock summed over threads) ---
  std::fprintf(out,
      "%-28s %14s %12s %10s %14s\n",
      "kernel", "thread-s", "calls", "%total", "us/call");
  std::fprintf(out, "%s\n", "----------------------------------------------"
                            "--------------------------------");
  for (int i = 0; i < N_TIMERS; ++i) {
    if (i == solve_rate_cool_total) continue;
    double us_per_call =
        tcalls[i] ? (tsec[i] / (double)tcalls[i]) * 1e6 : 0.0;
    std::fprintf(out, "%-28s %14.4f %12llu %9.1f%% %14.2f\n",
                 timer_name(i), tsec[i],
                 (unsigned long long)tcalls[i],
                 100.0 * tsec[i] / denom, us_per_call);
  }
  std::fprintf(out, "%-28s %14.4f\n", "solve_rate_cool_total (wall)", total);
  std::fprintf(out,
      "  (note: per-kernel sums are over all threads; total is wall-clock,\n"
      "   so kernel%% can exceed 100%% under OpenMP -- compare kernels to each\n"
      "   other, not to the wall total.)\n\n");

  // --- Iteration-count summary (for CPU/GPU iteration-divergence comparison) ---
  std::fprintf(out,
      "ITERATION COUNTS (for GPU/CPU convergence comparison):\n"
      "  Total subcycle iterations: %llu\n",
      (unsigned long long)cnt[subcycles]);
  // Sanity check: for kernels inside the subcycle loop (cool1d_multi_g,
  // lookup_cool_rates1d, etc.), their 'calls' should equal total subcycles.
  // Verify the most representative one:
  if (tcalls[cool1d_multi_g] != cnt[subcycles]) {
    std::fprintf(out,
        "  WARNING: cool1d_multi_g calls (%llu) != subcycles counter (%llu)\n",
        (unsigned long long)tcalls[cool1d_multi_g],
        (unsigned long long)cnt[subcycles]);
  } else {
    std::fprintf(out,
        "  --> The 'calls' column for subcycle-loop kernels (cool1d_multi_g,\n"
        "      lookup_cool_rates1d, rate_timestep, step_rate_*) equals the\n"
        "      subcycle count. Compare this across CPU/GPU runs to detect\n"
        "      iteration-count divergence due to numerical differences.\n");
  }
  std::fprintf(out, "\n");

  // --- Cell-path distribution ---
  std::fprintf(out, "cell-path counters:\n");
  for (int i = 0; i < N_COUNTERS; ++i) {
    std::fprintf(out, "  %-22s %16llu\n", counter_name(i),
                 (unsigned long long)cnt[i]);
  }
  uint64_t path_total = cnt[cells_gs] + cnt[cells_nr];
  if (path_total) {
    std::fprintf(out,
        "  --> Gauss-Seidel: %.2f%%   Newton-Raphson: %.2f%% "
        "(of which %.2f%% co-evolve energy)\n",
        100.0 * (double)cnt[cells_gs] / (double)path_total,
        100.0 * (double)cnt[cells_nr] / (double)path_total,
        cnt[cells_nr] ? 100.0 * (double)cnt[cells_nr_coevolve] /
                            (double)cnt[cells_nr]
                      : 0.0);
  }
  if (cnt[islices]) {
    std::fprintf(out, "  --> mean subcycles / i-slice: %.2f\n",
                 (double)cnt[subcycles] / (double)cnt[islices]);
  }

  // --- Newton-Raphson internal breakdown ---
  // Where does the per-NR-cell cost go: residual eval, finite-difference
  // Jacobian (~nsp derivative evals), or the dense linear solve?
  double nr_total = tsec[nr_residual] + tsec[nr_jacobian] + tsec[nr_gaussj];
  if (nr_total > 0.0) {
    std::fprintf(out, "\nNewton-Raphson internal breakdown:\n");
    std::fprintf(out, "  %-14s %12.4f s  %6.1f%%\n", "residual",
                 tsec[nr_residual], 100.0 * tsec[nr_residual] / nr_total);
    std::fprintf(out, "  %-14s %12.4f s  %6.1f%%   (finite-diff, ~nsp evals)\n",
                 "jacobian", tsec[nr_jacobian],
                 100.0 * tsec[nr_jacobian] / nr_total);
    std::fprintf(out, "  %-14s %12.4f s  %6.1f%%   (gaussj dense solve)\n",
                 "linear solve", tsec[nr_gaussj],
                 100.0 * tsec[nr_gaussj] / nr_total);
    if (cnt[newton_iterations]) {
      std::fprintf(out, "  total Newton iterations: %llu\n",
                   (unsigned long long)cnt[newton_iterations]);
      if (cnt[cells_nr]) {
        std::fprintf(out, "  --> mean Newton iters / NR cell-subcycle: %.2f\n",
                     (double)cnt[newton_iterations] / (double)cnt[cells_nr]);
      }
    }
  }

  // --- Subcycle-count histogram (drives GPU divergence analysis) ---
  std::fprintf(out, "\nsubcycle-count histogram (per i-slice):\n");
  // (the accurate mean is the "mean subcycles / i-slice" printed above, which
  // is subcycles/islices; here we only summarize the busiest bin, since a
  // histogram-derived mean using bin edges would be biased.)
  int max_bin_with_data = -1;
  for (int b = 0; b < N_HIST_BINS; ++b) {
    if (hist[b] > 0) max_bin_with_data = b;
  }
  if (max_bin_with_data >= 0) {
    std::fprintf(out, "  busiest i-slice falls in [%d .. %d] subcycles\n",
                 1 << max_bin_with_data,
                 (1 << (max_bin_with_data + 1)) - 1);
  }

  for (int b = 0; b < N_HIST_BINS; ++b) {
    if (!hist[b]) continue;
    int lo = 1 << b;
    int hi = (1 << (b + 1)) - 1;
    std::fprintf(out, "  [%6d .. %6d] %12llu\n", lo, hi,
                 (unsigned long long)hist[b]);
  }
  std::fprintf(out,
      "================================================================\n\n");

  if (close_out) std::fclose(out);
}

/// Register report() to run once at program exit.
///
/// This is what GRACKLE_PROF_REPORT() actually calls. Reporting at exit (rather
/// than at the end of each solve_rate_cool invocation) means:
///   * the accumulators sum over *all* invocations (whole-run totals), and
///   * the per-call `solve_rate_cool_total` ScopedTimer has already destructed
///     and recorded its time before we read it.
/// The function-local static guarantees the std::atexit registration happens
/// exactly once, thread-safely (C++11 static-init rules).
inline void report_atexit() {
  static int once = (std::atexit(report), 0);
  (void)once;
}

}  // namespace grackle::impl::prof

// ---------------------------------------------------------------------------
// Public macros (enabled build)
// ---------------------------------------------------------------------------
#define GRACKLE_PROF_CONCAT_(a, b) a##b
#define GRACKLE_PROF_CONCAT(a, b) GRACKLE_PROF_CONCAT_(a, b)

#define GRACKLE_PROF_SCOPE(name)                                               \
  ::grackle::impl::prof::ScopedTimer GRACKLE_PROF_CONCAT(_prof_scope_,         \
                                                         __LINE__)(            \
      ::grackle::impl::prof::TimerId::name)

#define GRACKLE_PROF_COUNT(name, n)                                           \
  ::grackle::impl::prof::add_counter(::grackle::impl::prof::CounterId::name,   \
                                     (uint64_t)(n))

#define GRACKLE_PROF_HIST_SUBCYCLE(n)                                         \
  ::grackle::impl::prof::record_subcycle_hist((int)(n))

// Registers a single report() dump at program exit (safe to call repeatedly;
// only the first call registers). Use GRACKLE_PROF_REPORT_NOW() to force an
// immediate dump instead.
#define GRACKLE_PROF_REPORT() ::grackle::impl::prof::report_atexit()
#define GRACKLE_PROF_REPORT_NOW() ::grackle::impl::prof::report()
#define GRACKLE_PROF_RESET() ::grackle::impl::prof::reset()

#endif  // GRACKLE_PROFILE

#endif  // GRACKLE_SUPPORT_PROFILING_HPP

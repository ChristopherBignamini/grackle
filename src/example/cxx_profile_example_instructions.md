# Profiling `solve_rate_cool` with `cxx_profile_example`

A driver + in-library timing harness for measuring where the chemistry/cooling
solver spends its time. The harness (`src/clib/support/profiling.hpp`) is a
compile-time no-op unless the library is built with `-DGRACKLE_PROFILE=ON`.

## Build

Profiling is a **library** build flag — it must be set when the library
compiles, not just the example. Use a fresh build dir.

```bash
cmake -B build \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \   # optimized + debug symbols
  -DGRACKLE_USE_OPENMP=ON \             # for multithreaded runs
  -DGRACKLE_PROFILE=ON                  # <-- enables the timing harness
cmake --build build --target cxx_profile_example -j
```

Sanity check: if the profile table prints at program exit, it worked. If you
only see `done (N iterations).` and no table, the library was built **without**
`-DGRACKLE_PROFILE=ON`.

## Run

Run **from `build/examples`** so the example's relative data-file path
(`../../input/CloudyData_UVB=HM2012.h5`) resolves.

```bash
cd build/examples
OMP_SCHEDULE=dynamic ./cxx_profile_example -s 2 -P 5.5 -n 32 -t 8 -a 5 -d 10
```

The per-kernel timing + cell-path + subcycle-histogram table is printed at exit.

## Flags

| flag | meaning | default |
|------|---------|---------|
| `-n` | cells per dimension (grid is n^3) | 16 |
| `-t` | OpenMP threads | max |
| `-a` | iterations (repeats the solve, amortizes init) | 10 |
| `-s` | solver: 1=auto hybrid, 2=force Gauss-Seidel, 3=force Newton-Raphson | 1 |
| `-i` | subcycle cap per cell (bounds runaway stiff cells) | 1000 |
| `-d` | timestep in years | 1.0 |
| `-p` / `-P` | log10 min/max hydrogen number density (cm^-3) | -3 / 8.5 |
| `-y` / `-Y` | log10 min/max temperature (K) | 2 / 6 |
| `-m` | multi_metals (1 is unsupported by this driver) | 0 |

The grid maps each axis to a parameter sweep: x -> density, y -> temperature,
z -> metallicity. A single run therefore covers a broad slice of parameter
space at once.

## Runtime environment

| var | effect |
|-----|--------|
| `GRACKLE_PROFILE_OUT=<file>` | append the profile table to a file instead of stderr |
| `OMP_NUM_THREADS` | thread count (or use `-t`) |
| `OMP_SCHEDULE=dynamic` | dynamic loop scheduling — important on heterogeneous CPUs and for the uneven per-cell cost |

## Useful recipes

```bash
# Full maximal config, realistic solver mix (both GS and NR):
./cxx_profile_example -n 16 -t 8 -a 5

# Gauss-Seidel path only, in the density regime where GS actually runs
# (forcing GS onto dense cells produces NaNs):
./cxx_profile_example -s 2 -P 5.5 -n 32 -t 8 -a 5 -d 10

# Dust in the regime where grains survive (below ~2000 K); above the grain
# sublimation temperature the dust-temperature solve is forced into bisection:
./cxx_profile_example -s 2 -P 5.5 -Y 3.3 -n 32 -t 8 -a 5 -d 10

# Single-thread run for clean per-kernel attribution (no allocator/bandwidth
# contention inflating thread-seconds):
./cxx_profile_example -s 2 -P 5.5 -n 32 -t 1 -a 5 -d 10
```

## Reading the table

- **Per-kernel table**: compare kernels to *each other*. `%total` can exceed
  100% because it sums thread-seconds over all threads against a wall-clock
  total; that is expected under OpenMP.
- **`cells_gs` / `cells_nr`**: how cells split across the two solvers.
- **`islices_maxed_out`**: i-slices that hit the `-i` subcycle cap (did not
  converge). A large fraction means the run is cap-limited, not converged —
  lower `-d` or raise `-i`.
- **subcycle histogram**: per-i-slice iteration spread — the key GPU
  thread-divergence signal.

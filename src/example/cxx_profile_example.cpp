/***********************************************************************
/
/ Profiling driver for the GPU-porting / optimization effort.
/
/ Exercises the *most complete* Grackle configuration (primordial_chemistry=4
/ + metal_chemistry + full dust-species network) over a 3D grid whose cells
/ span a wide density AND temperature range. The wide density range is the key
/ ingredient: it drives cells across the Gauss-Seidel <-> Newton-Raphson
/ threshold (ddom ~ 1e8), so BOTH chemistry solvers -- and the full set of
/ per-cell kernels -- are exercised in a single run.
/
/ This driver produces no output of its own beyond a short banner. The useful
/ output is the per-kernel timing + cell-path + subcycle-histogram table emitted
/ at program exit by the GRACKLE_PROFILE harness (support/profiling.hpp). Build
/ the library with -DGRACKLE_PROFILE=ON (CMake: -DGRACKLE_PROFILE=ON) to enable
/ it; otherwise this is just a plain (silent) correctness run.
/
/ Distributed under the terms of the Enzo Public Licence.
/ The full license is in the file LICENSE, distributed with this software.
************************************************************************/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <grackle.h>

#define mh     1.67262171e-24
#define kboltz 1.3806504e-16

// ---------------------------------------------------------------------------
// Parameter-space sweep, mapped onto the three grid axes so that a single grid
// covers a broad slice of the problem space:
//   x-axis (i) -> hydrogen number density (log-spaced)  [drives GS vs NR]
//   y-axis (j) -> gas temperature         (log-spaced)  [drives table lookups]
//   z-axis (k) -> metallicity in solar units (log-spaced) [drives metal path]
// ---------------------------------------------------------------------------
// Density-axis defaults; override at runtime with -p / -P (log10 nH min/max).
//
// With density_units == mh the internal `dom` factor is ~1, so ddom ~ nH/0.76.
// Reference points on this axis:
//   nH ~ 7.6e5 (log 5.88) -> ddom = 1e6 : NR threshold WHEN metals are present
//   nH ~ 7.6e7 (log 7.88) -> ddom = 1e8 : NR threshold in general
//
// Choosing the ceiling is a trade-off:
//  - it must clear the NR threshold to exercise Newton-Raphson at all, but the
//    further above it you go the stiffer the cell. The per-subcycle step is
//    stiffness-limited, so the subcycle count needed to integrate `dt` grows
//    ~linearly with nH. The default (3e8, ~0.6 decade above the threshold)
//    gives a real NR band while still converging within the subcycle cap at the
//    small default dt. Raise it (with -P) for more NR stress -- and lower dt to
//    keep it converging.
//  - to profile the *Gauss-Seidel* path in the regime it actually runs in
//    production, use `-s 2 -P 5.5`. Forcing GS onto denser cells makes it
//    diverge into NaNs: those cells are precisely what NR exists for.
static const double DEFAULT_LOG_NH_MIN = -3.0;  // cm^-3  (very diffuse)
static const double DEFAULT_LOG_NH_MAX =  8.5;  // cm^-3  (~0.6 dec above NR thresh)
// Temperature axis (override with -y / -Y). Note: above the grain sublimation
// temperature (~1500-2000 K, log10 ~ 3.2-3.3) the dust-temperature solve is
// forced into its (expensive) bisection branch. Cap at -Y 3.3 to profile dust
// in the regime where grains physically survive.
static const double DEFAULT_LOG_T_MIN =  2.0;   // K      (100 K, molecular regime)
static const double DEFAULT_LOG_T_MAX =  6.0;   // K      (1e6 K, coll.-ionization)
static const double LOG_ZSOL_MIN = -4.0;   // Z/Zsun
static const double LOG_ZSOL_MAX =  0.0;   // Z/Zsun

static double logspace(double lo, double hi, int idx, int n) {
  double frac = (n > 1) ? (double)idx / (double)(n - 1) : 0.0;
  return std::pow(10.0, lo + (hi - lo) * frac);
}

int main(int argc, char* argv[]) {

  // ---- command-line options ----
#ifdef _OPENMP
  int NThread = omp_get_max_threads();
#else
  int NThread = 1;
#endif
  int NIter = 10, NCell1D = 16, c;
  int SolverMethod = 1;      // 1=auto (density split), 2=force GS, 3=force NR
  int MultiMetals  = 0;      // 0=single metal source, 1=per-pathway (see below)
  int MaxIter      = 1000;   // per-cell subcycle cap (bounds worst-case work)
  double dt_years  = 1.0;    // timestep in years. The subcycle count for a stiff
                             // cell scales with the time span integrated, so a
                             // short dt is what lets high-density (NR) cells
                             // finish within MaxIter. Raise it to stress the
                             // solver (and expect the stiffest cells to max out).
  double log_nh_min = DEFAULT_LOG_NH_MIN;  // -p : log10 of min hydrogen density
  double log_nh_max = DEFAULT_LOG_NH_MAX;  // -P : log10 of max hydrogen density
  double log_t_min  = DEFAULT_LOG_T_MIN;   // -y : log10 of min temperature (K)
  double log_t_max  = DEFAULT_LOG_T_MAX;   // -Y : log10 of max temperature (K)
  while ((c = getopt(argc, argv, "ht:a:n:s:m:i:d:p:P:y:Y:")) != -1) {
    switch (c) {
      case 't': NThread      = atoi(optarg); break;
      case 'a': NIter        = atoi(optarg); break;
      case 'n': NCell1D      = atoi(optarg); break;
      case 's': SolverMethod = atoi(optarg); break;
      case 'm': MultiMetals  = atoi(optarg); break;
      case 'i': MaxIter      = atoi(optarg); break;
      case 'd': dt_years     = atof(optarg); break;
      case 'p': log_nh_min   = atof(optarg); break;
      case 'P': log_nh_max   = atof(optarg); break;
      case 'y': log_t_min    = atof(optarg); break;
      case 'Y': log_t_max    = atof(optarg); break;
      case 'h':
      case '?':
      default:
        fprintf(stderr,
                "usage: %s [-t nthreads] [-a niters] [-n ncells_per_dim]\n"
                "          [-s solver] [-m multi_metals] [-i max_iter] [-d dt_yr]\n"
                "          [-p log10_nH_min] [-P log10_nH_max]\n"
                "          [-y log10_T_min] [-Y log10_T_max]\n",
                argv[0]);
        fprintf(stderr,
                "  Runs the full primordial_chemistry=4 + metal + dust config\n"
                "  over an n^3 grid spanning density/temperature/metallicity.\n"
                "  Build the library with -DGRACKLE_PROFILE=ON to get the\n"
                "  per-kernel profiling table (dumped at exit).\n\n"
                "  -s solver : 1=auto density split [default], 2=force\n"
                "              Gauss-Seidel everywhere, 3=force Newton-Raphson\n"
                "              everywhere. Use 2 vs 3 for A/B solver profiling.\n"
                "  -m 0|1    : multi_metals (0=single metal source [default]).\n"
                "              1 is currently unsupported by this driver (see\n"
                "              the error message it prints).\n"
                "  -i max_iter : subcycle cap per cell [1000]. Bounds worst-case\n"
                "              work; the stiff high-density cells otherwise run\n"
                "              to the library default of 10000 and appear hung.\n"
                "  -d dt_yr  : timestep in years [1e3]. Larger => more subcycles.\n"
                "  -p, -P    : log10 of the min/max hydrogen number density\n"
                "              [%g .. %g] cm^-3. Reference points: log10(nH)=5.88\n"
                "              is the NR threshold with metals, 7.88 without.\n"
                "              To profile the Gauss-Seidel path in the regime it\n"
                "              actually runs, use -s 2 -P 5.5 (forcing GS onto\n"
                "              denser cells makes it produce NaNs).\n"
                "  -y, -Y    : log10 of the min/max temperature [%g .. %g] K.\n"
                "              Above the grain sublimation temp (log10 ~ 3.3) the\n"
                "              dust-temperature solve is forced into bisection;\n"
                "              use -Y 3.3 to profile dust where grains survive.\n\n"
                "  Tip: start small to confirm it completes, e.g. -n 4 -a 1.\n",
                DEFAULT_LOG_NH_MIN, DEFAULT_LOG_NH_MAX,
                DEFAULT_LOG_T_MIN, DEFAULT_LOG_T_MAX);
        exit(1);
    }
  }
  if (NThread < 1 || NIter < 1 || NCell1D < 1 || MaxIter < 1 || dt_years <= 0.0) {
    fprintf(stderr, "ERROR: -t, -a, -n, -i must be >= 1 and -d > 0\n");
    exit(EXIT_FAILURE);
  }
  if (SolverMethod < 1 || SolverMethod > 3) {
    fprintf(stderr, "ERROR: -s must be 1 (auto), 2 (force GS), or 3 (force NR)\n");
    exit(EXIT_FAILURE);
  }
  if (log_nh_min >= log_nh_max) {
    fprintf(stderr, "ERROR: -p (%g) must be < -P (%g)\n", log_nh_min, log_nh_max);
    exit(EXIT_FAILURE);
  }
  if (log_t_min >= log_t_max) {
    fprintf(stderr, "ERROR: -y (%g) must be < -Y (%g)\n", log_t_min, log_t_max);
    exit(EXIT_FAILURE);
  }
  // Forcing Gauss-Seidel onto cells the hybrid would route to Newton-Raphson
  // makes GS diverge (NaN species densities). Warn rather than silently produce
  // garbage / crash.
  if (SolverMethod == 2 && log_nh_max > 5.88) {
    fprintf(stderr,
            "WARNING: -s 2 (force Gauss-Seidel) with -P %g exceeds the NR\n"
            "  threshold (log10 nH ~ 5.88 with metals). GS is not designed for\n"
            "  those cells and may produce NaNs. Consider -P 5.5\n",
            log_nh_max);
  }
  if (MultiMetals != 0) {
    // multi_metals=1 makes Grackle read per-injection-pathway metal densities
    // from my_fields->inject_pathway_metal_density[]. The number of pathways is
    // determined internally at init time and is not exposed through the public
    // API, so this driver cannot allocate/populate those fields correctly.
    // Leaving them unset would segfault inside solve_rate_cool. Fail loudly.
    fprintf(stderr,
            "ERROR: -m 1 (multi_metals) is not supported by this driver.\n"
            "  It requires my_fields->inject_pathway_metal_density[] to hold\n"
            "  n_pathways valid pointers, but n_pathways is only known inside\n"
            "  the library after initialization. Run with -m 0.\n");
    exit(EXIT_FAILURE);
  }

  if (gr_check_consistency() != GR_SUCCESS) {
    fprintf(stderr, "Error in gr_check_consistency.\n");
    return EXIT_FAILURE;
  }
  // keep verbose off: at high density the stiff cells emit per-i-slice
  // "MULTI_COOL iter" warnings, and we only care about the profiling table
  grackle_verbose = 0;

  // ---- units: non-cosmological, z=0 (mirrors cxx_omp_example) ----
  code_units my_units;
  my_units.comoving_coordinates = 0;
  my_units.density_units        = mh;      // 1 code density unit ~ 1 H atom/cm^3
  my_units.length_units         = 1.0;
  my_units.time_units           = 1.0e12;
  my_units.a_units              = 1.0;
  my_units.a_value              = 1.0;
  set_velocity_units(&my_units);

  // ---- chemistry parameters: the maximal (pc=4 + metal + dust) config ----
  // Based on the C++ test preset `primchem4_dustspecies3` (known-good with the
  // HM2012 data file). multi_metals is left at 0 to avoid the inject-pathway
  // field plumbing; flip it to 1 (and allocate inject_pathway_metal_density[])
  // as a follow-up if you want that path profiled too.
  chemistry_data* my_grackle_data = new chemistry_data;
  if (set_default_chemistry_parameters(my_grackle_data) == 0) {
    fprintf(stderr, "Error in set_default_chemistry_parameters.\n");
    return EXIT_FAILURE;
  }

  grackle_data->use_grackle            = 1;
  grackle_data->with_radiative_cooling = 1;
  grackle_data->primordial_chemistry   = 4;   // full new-chem network + D
  grackle_data->metal_cooling          = 1;
  grackle_data->metal_chemistry        = 1;   // molecular/atomic metal network
  grackle_data->dust_chemistry         = 1;
  grackle_data->dust_species           = 3;   // full grain-species set
  grackle_data->use_dust_density_field = 1;
  grackle_data->multi_metals           = MultiMetals;   // 0 (enforced above)
  grackle_data->solver_method          = SolverMethod;  // 1=auto, 2=GS, 3=NR
  grackle_data->max_iterations         = MaxIter;       // bound worst-case work
  grackle_data->UVbackground           = 1;
  grackle_data->use_isrf_field         = 1;
  grackle_data->grackle_data_file      = "../../input/CloudyData_UVB=HM2012.h5";
#ifdef _OPENMP
  grackle_data->omp_nthreads           = NThread;
#endif

  fprintf(stdout, "[progress] initializing chemistry data ...\n");
  fflush(stdout);
  if (initialize_chemistry_data(&my_units) == 0) {
    fprintf(stderr, "Error in initialize_chemistry_data.\n");
    return EXIT_FAILURE;
  }
  fprintf(stdout, "[progress] chemistry data initialized\n");
  fflush(stdout);

  // ---- grid setup ----
  const int N3 = NCell1D * NCell1D * NCell1D;
  grackle_field_data my_fields;
  gr_initialize_field_data(&my_fields);   // NULLs every field pointer
  my_fields.grid_rank      = 3;
  my_fields.grid_dimension = new int[3];
  my_fields.grid_start     = new int[3];
  my_fields.grid_end       = new int[3];
  my_fields.grid_dx        = 0.0;
  for (int d = 0; d < 3; d++) {
    my_fields.grid_dimension[d] = NCell1D;
    my_fields.grid_start[d]     = 0;
    my_fields.grid_end[d]       = NCell1D - 1;
  }

  // ---- allocate fields ----
  auto A = [&]() { return new gr_float[N3]; };

  // core / bulk-carrying fields
  my_fields.density         = A();
  my_fields.internal_energy = A();
  my_fields.x_velocity      = A();
  my_fields.y_velocity      = A();
  my_fields.z_velocity      = A();
  my_fields.metal_density   = A();
  my_fields.dust_density    = A();
  my_fields.isrf_habing     = A();
  my_fields.HI_density      = A();
  my_fields.HeI_density     = A();
  my_fields.DI_density      = A();

  // species initialized as a tiny fraction of the *gas* density
  std::vector<gr_float*> tiny_of_gas;
  auto Ag = [&]() { gr_float* p = A(); tiny_of_gas.push_back(p); return p; };
  my_fields.e_density     = Ag();
  my_fields.HII_density   = Ag();
  my_fields.HeII_density  = Ag();
  my_fields.HeIII_density = Ag();
  my_fields.HM_density    = Ag();
  my_fields.H2I_density   = Ag();
  my_fields.H2II_density  = Ag();
  my_fields.DII_density   = Ag();
  my_fields.HDI_density   = Ag();
  my_fields.DM_density    = Ag();
  my_fields.HDII_density  = Ag();
  my_fields.HeHII_density = Ag();

  // species initialized as a tiny fraction of the *metal* density
  std::vector<gr_float*> tiny_of_metal;
  auto Am = [&]() { gr_float* p = A(); tiny_of_metal.push_back(p); return p; };
  my_fields.CI_density    = Am();
  my_fields.CII_density   = Am();
  my_fields.CO_density    = Am();
  my_fields.CO2_density   = Am();
  my_fields.OI_density    = Am();
  my_fields.OH_density    = Am();
  my_fields.H2O_density   = Am();
  my_fields.O2_density    = Am();
  my_fields.SiI_density   = Am();
  my_fields.SiOI_density  = Am();
  my_fields.SiO2I_density = Am();
  my_fields.CH_density    = Am();
  my_fields.CH2_density   = Am();
  my_fields.COII_density  = Am();
  my_fields.OII_density   = Am();
  my_fields.OHII_density  = Am();
  my_fields.H2OII_density = Am();
  my_fields.H3OII_density = Am();
  my_fields.O2II_density  = Am();
  my_fields.Mg_density    = Am();
  my_fields.Al_density    = Am();
  my_fields.S_density     = Am();
  my_fields.Fe_density    = Am();

  // dust grain species initialized as a tiny fraction of the *dust* density
  std::vector<gr_float*> tiny_of_dust;
  auto Ad = [&]() { gr_float* p = A(); tiny_of_dust.push_back(p); return p; };
  my_fields.MgSiO3_dust_density  = Ad();
  my_fields.AC_dust_density      = Ad();
  my_fields.SiM_dust_density     = Ad();
  my_fields.FeM_dust_density     = Ad();
  my_fields.Mg2SiO4_dust_density = Ad();
  my_fields.Fe3O4_dust_density   = Ad();
  my_fields.SiO2_dust_density    = Ad();
  my_fields.MgO_dust_density     = Ad();
  my_fields.FeS_dust_density     = Ad();
  my_fields.Al2O3_dust_density   = Ad();
  my_fields.ref_org_dust_density = Ad();
  my_fields.vol_org_dust_density = Ad();
  my_fields.H2O_ice_dust_density = Ad();

  // ---- initialize field values ----
  const double temperature_units = get_temperature_units(&my_units);
  const double H_frac  = grackle_data->HydrogenFractionByMass;
  const double Zsol    = grackle_data->SolarMetalFractionByMass;
  const double d2g_sol = grackle_data->local_dust_to_gas_ratio;
  const gr_float tiny  = 1.0e-20;

  for (int k = 0; k < NCell1D; k++) {
    double Zsolar = logspace(LOG_ZSOL_MIN, LOG_ZSOL_MAX, k, NCell1D); // Z/Zsun
    for (int j = 0; j < NCell1D; j++) {
      double T = logspace(log_t_min, log_t_max, j, NCell1D);          // K
      for (int i = 0; i < NCell1D; i++) {
        double nH = logspace(log_nh_min, log_nh_max, i, NCell1D);     // cm^-3

        int idx = i + NCell1D * (j + NCell1D * k);

        // mass density in code units (density_units == mh, H_frac of mass is H)
        gr_float rho = (gr_float)(nH * mh / my_units.density_units / H_frac);
        my_fields.density[idx]         = rho;
        my_fields.x_velocity[idx]      = 0.0;
        my_fields.y_velocity[idx]      = 0.0;
        my_fields.z_velocity[idx]      = 0.0;
        // internal energy for target T (mu~1 approximation; Grackle recomputes T)
        my_fields.internal_energy[idx] = (gr_float)(T / temperature_units);
        my_fields.isrf_habing[idx]     = 1.0;

        gr_float metal = (gr_float)(Zsolar * Zsol * rho);
        gr_float dust  = (gr_float)(Zsolar * d2g_sol * rho);
        my_fields.metal_density[idx]   = metal;
        my_fields.dust_density[idx]    = dust;

        // Seed every evolved species with a tiny floor first...
        for (gr_float* p : tiny_of_gas)   p[idx] = tiny * rho;
        for (gr_float* p : tiny_of_metal) p[idx] = tiny * metal;
        for (gr_float* p : tiny_of_dust)  p[idx] = tiny * dust;

        // ...then override H/He with a TEMPERATURE-APPROPRIATE ionization
        // state. Starting the gas neutral at high T (as a naive init does) puts
        // it far from equilibrium: the collisional-ionization/cooling rate then
        // blows up (edot -> +/-1e50 / NaN) and the subcycle never converges.
        // fion : neutral->singly-ionized transition (~1e4 K)
        // fion2: HeII->HeIII transition (~1e5 K)
        // The electron density matches make_consistent's charge balance
        // (de = HII + HeII/4 + HeIII/2).
        double logT  = std::log10(T);
        double fion  = 1.0 / (1.0 + std::pow(10.0, -(logT - 4.2) / 0.3));
        double fion2 = 1.0 / (1.0 + std::pow(10.0, -(logT - 5.0) / 0.3));
        double He_mass = (1.0 - H_frac) * rho;   // metals are trace mass

        gr_float HII   = (gr_float)(fion * H_frac * rho);
        gr_float HeII  = (gr_float)((fion - fion2) * He_mass);
        gr_float HeIII = (gr_float)(fion2 * He_mass);
        my_fields.HI_density[idx]    = (gr_float)((1.0 - fion) * H_frac * rho);
        my_fields.HII_density[idx]   = HII;
        my_fields.HeI_density[idx]   = (gr_float)((1.0 - fion) * He_mass);
        my_fields.HeII_density[idx]  = HeII;
        my_fields.HeIII_density[idx] = HeIII;
        my_fields.e_density[idx]     =
            HII + HeII / (gr_float)4.0 + HeIII / (gr_float)2.0;
        my_fields.DI_density[idx]    = (gr_float)(2.0 * 3.4e-5 * rho);
      }
    }
  }

  // ---- run ----
  const double dt = dt_years * 3.15e7 / my_units.time_units;  // years -> code
#ifdef _OPENMP
  omp_set_num_threads(NThread);
#endif
  const char* solver_desc = (SolverMethod == 2) ? "force Gauss-Seidel"
                          : (SolverMethod == 3) ? "force Newton-Raphson"
                          : "auto (density split)";
  fprintf(stdout,
          "profile run: %d^3 = %d cells, %d threads, %d iterations\n"
          "  solver   : method=%d (%s)\n"
          "  dt       : %g yr    max subcycles/cell : %d\n"
          "  density  : nH   in [1e%g .. 1e%g] cm^-3  (spans GS/NR threshold)\n"
          "  temp     : T    in [1e%g .. 1e%g] K\n"
          "  metals   : Z    in [1e%g .. 1e%g] Zsun\n"
          "  (profiling table, if enabled, is printed at exit)\n",
          NCell1D, N3, NThread, NIter, SolverMethod, solver_desc,
          dt_years, MaxIter,
          log_nh_min, log_nh_max, log_t_min, log_t_max,
          LOG_ZSOL_MIN, LOG_ZSOL_MAX);
  fprintf(stdout, "[progress] fields initialized, starting solve loop\n");
  fflush(stdout);

  for (int t = 0; t < NIter; t++) {
#ifdef _OPENMP
    double t0 = omp_get_wtime();
#endif
    if (solve_chemistry(&my_units, &my_fields, dt) == 0) {
      fprintf(stderr, "Error in solve_chemistry (iteration %d).\n", t);
      return EXIT_FAILURE;
    }
#ifdef _OPENMP
    fprintf(stdout, "[progress] iteration %d/%d done (%.3f s)\n",
            t + 1, NIter, omp_get_wtime() - t0);
#else
    fprintf(stdout, "[progress] iteration %d/%d done\n", t + 1, NIter);
#endif
    fflush(stdout);
  }

  fprintf(stdout, "done (%d iterations).\n", NIter);

  // ---- cleanup ----
  delete[] my_fields.density;
  delete[] my_fields.internal_energy;
  delete[] my_fields.x_velocity;
  delete[] my_fields.y_velocity;
  delete[] my_fields.z_velocity;      
  delete[] my_fields.metal_density;
  delete[] my_fields.dust_density;    
  delete[] my_fields.isrf_habing;
  delete[] my_fields.HI_density;      
  delete[] my_fields.HeI_density;
  delete[] my_fields.DI_density;
  for (gr_float* p : tiny_of_gas)   delete[] p;
  for (gr_float* p : tiny_of_metal) delete[] p;
  for (gr_float* p : tiny_of_dust)  delete[] p;
  delete[] my_fields.grid_dimension;
  delete[] my_fields.grid_start;
  delete[] my_fields.grid_end;

  return 0;
}

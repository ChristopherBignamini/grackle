#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdint>  // per tipi interi a 64 bit
#include <omp.h>    // per OpenMP, se necessario

const double DKIND = 1.0;  // Usato per coerenza con il codice Fortran

void cool1D_cloudy_old_tables_g(
    const std::vector<std::vector<std::vector<double>>>& d,
    const std::vector<std::vector<std::vector<double>>>& de,
    const std::vector<double>& rhoH,
    const std::vector<double>& metallicity,
    int in, int jn, int kn, int is, int ie, int j, int k,
    const std::vector<double>& logtem,
    std::vector<double>& edot,
    double comp2, int ispecies, double dom, double zr,
    int icmbTfloor, int iClHeat,
    double clEleFra, int64_t clGridRank, const std::vector<int64_t>& clGridDim,
    const std::vector<double>& clPar1, const std::vector<double>& clPar2,
    const std::vector<double>& clPar3, const std::vector<double>& clPar4,
    const std::vector<double>& clPar5,
    int64_t clDataSize, const std::vector<double>& clCooling,
    const std::vector<double>& clHeating,
    const std::vector<bool>& itmask)
{
    // Variabili locali
    double inv_log10 = 1.0 / std::log(10.0);
    double log10_tCMB = std::log10(comp2);

    // Calcola le pendenze dei parametri
    std::vector<double> dclPar(clGridRank);
    dclPar[0] = (clPar1[clGridDim[0] - 1] - clPar1[0]) / static_cast<double>(clGridDim[0] - 1);
    if (clGridRank > 1) {
        dclPar[1] = (clPar2[clGridDim[1] - 1] - clPar2[0]) / static_cast<double>(clGridDim[1] - 1);
    }
    if (clGridRank > 2) {
        dclPar[2] = (clPar3[clGridDim[2] - 1] - clPar3[0]) / static_cast<double>(clGridDim[2] - 1);
    }
    if (clGridRank > 3) {
        dclPar[3] = (clPar4[clGridDim[3] - 1] - clPar4[0]) / static_cast<double>(clGridDim[3] - 1);
    }
    if (clGridRank > 4) {
        dclPar[4] = (clPar5[clGridDim[4] - 1] - clPar5[0]) / static_cast<double>(clGridDim[4] - 1);
    }

    // Vettori locali per la slice
    std::vector<double> log_Z(in);
    std::vector<double> e_frac(in);
    std::vector<double> log_e_frac(in);
    std::vector<double> cl_e_frac(in);
    std::vector<double> fh(in);
    std::vector<double> log_n_h(in);
    std::vector<double> log_cool(in);
    std::vector<double> log_cool_cmb(in);
    std::vector<double> log_heat(in);
    std::vector<double> edot_met(in);
    std::vector<double> log10tem(in);

    for (int i = is; i <= ie; ++i) {
        if (itmask[i]) {
            log10tem[i] = logtem[i] * inv_log10;

            // Calcola la frazione di idrogeno
            fh[i] = rhoH[i] / d[i][j][k];

            // Calcola il logaritmo della densità dell'idrogeno in unità proprie
            if (clGridRank > 1) {
                log_n_h[i] = std::log10(rhoH[i] * dom);
            }

            // Calcola il logaritmo della metallicità
            if (clGridRank > 2) {
                log_Z[i] = std::log10(metallicity[i]);
            }

            // Calcola la frazione di elettroni
            if (clGridRank > 3) {
                e_frac[i] = 2.0 * de[i][j][k] / (d[i][j][k] * (1.0 + fh[i]));
                // Assicura che la frazione di elettroni non superi 1
                log_e_frac[i] = std::min(std::log10(e_frac[i]), 0.0);

                // Calcola gli elettroni extra contribuiti dai metalli
                cl_e_frac[i] = e_frac[i] * (1.0 + (2.0 * clEleFra * metallicity[i] * fh[i]) / (1.0 + fh[i]));
            }

            // Interpolazione per ottenere il raffreddamento/riscaldamento
            if (clGridRank == 1) {
                interpolate_1D_g(log10tem[i], clGridDim, clPar1, dclPar[0], clDataSize, clCooling, log_cool[i]);
                edot_met[i] = -std::pow(10.0, log_cool[i]);

                // Ignora il termine del CMB se T >> T_CMB
                if ((icmbTfloor == 1) && ((log10tem[i] - log10_tCMB) < 2.0)) {
                    interpolate_1D_g(log10_tCMB, clGridDim, clPar1, dclPar[0], clDataSize, clCooling, log_cool_cmb[i]);
                    edot_met[i] += std::pow(10.0, log_cool_cmb[i]);
                }

                if (iClHeat == 1) {
                    interpolate_1D_g(log10tem[i], clGridDim, clPar1, dclPar[0], clDataSize, clHeating, log_heat[i]);
                    edot_met[i] += std::pow(10.0, log_heat[i]);
                }
            } else if (clGridRank == 2) {
                interpolate_2D_g(log_n_h[i], log10tem[i], clGridDim, clPar1, dclPar[0], clPar2, dclPar[1], clDataSize, clCooling, log_cool[i]);
                edot_met[i] = -std::pow(10.0, log_cool[i]);

                // Ignora il termine del CMB se T >> T_CMB
                if ((icmbTfloor == 1) && ((log10tem[i] - log10_tCMB) < 2.0)) {
                    interpolate_2D_g(log_n_h[i], log10_tCMB, clGridDim, clPar1, dclPar[0], clPar2, dclPar[1], clDataSize, clCooling, log_cool_cmb[i]);
                    edot_met[i] += std::pow(10.0, log_cool_cmb[i]);
                }

                if (iClHeat == 1) {
                    interpolate_2D_g(log_n_h[i], log10tem[i], clGridDim, clPar1, dclPar[0], clPar2, dclPar[1], clDataSize, clHeating, log_heat[i]);
                    edot_met[i] += std::pow(10.0, log_heat[i]);
                }
            } else if (clGridRank == 3) {
                interpolate_3D_g(log_n_h[i], log_Z[i], log10tem[i], clGridDim, clPar1, dclPar[0], clPar2, dclPar[1], clPar3, dclPar[2], clDataSize, clCooling, log_cool[i]);
                edot_met[i] = -std::pow(10.0, log_cool[i]);

                // Ignora il termine del CMB se T >> T_CMB
                if ((icmbTfloor == 1) && ((log10tem[i] - log10_tCMB) < 2.0)) {
                    interpolate_3D_g(log_n_h[i], log_Z[i], log10_tCMB, clGridDim, clPar1, dclPar[0], clPar2, dclPar[1], clPar3, dclPar[2], clDataSize, clCooling, log_cool_cmb[i]);
                    edot_met[i] += std::pow(10.0, log_cool_cmb[i]);
                }

                if (iClHeat == 1) {
                    interpolate_3D_g(log_n_h[i], log_Z[i], log10tem[i], clGridDim, clPar1, dclPar[0], clPar2, dclPar[1], clPar3, dclPar[2], clDataSize, clHeating, log_heat[i]);
                    edot_met[i] += std::pow(10.0, log_heat[i]);
                }
            } else if (clGridRank == 4) {
                interpolate_4D_g(log_n_h[i], log_Z[i], log_e_frac[i], log10tem[i], clGridDim, clPar1, dclPar[0], clPar2, dclPar[1], clPar3, dclPar[2], clPar4, dclPar[3], clDataSize, clCooling, log_cool[i]);
                edot_met[i] = -std::pow(10.0, log_cool[i]);

                // Ignora il termine del CMB se T >> T_CMB
                if ((icmbTfloor == 1) && ((log10tem[i] - log10_tCMB) < 2.0)) {
                    interpolate_4D_g(log_n_h[i], log_Z[i], log_e_frac[i], log10_tCMB, clGridDim, clPar1, dclPar[0], clPar2, dclPar[1], clPar3, dclPar[2], clPar4, dclPar[3], clDataSize, clCooling, log_cool_cmb[i]);
                    edot_met[i] += std::pow(10.0, log_cool_cmb[i]);
                }

                if (iClHeat == 1) {
                    interpolate_4D_g(log_n_h[i], log_Z[i], log_e_frac[i], log10tem[i], clGridDim, clPar1, dclPar[0], clPar2, dclPar[1], clPar3, dclPar[2], clPar4, dclPar[3], clDataSize, clHeating, log_heat[i]);
                    edot_met[i] += std::pow(10.0, log_heat[i]);
                }

                // Scala il contributo con la frazione di elettroni
                edot_met[i] *= cl_e_frac[i];
            } else {
                // Caso per clGridRank == 5
                interpolate_5D_g(log_n_h[i], log_Z[i], log_e_frac[i], zr, log10tem[i], clGridDim, clPar1, dclPar[0], clPar2, dclPar[1], clPar3, dclPar[2], clPar4, dclPar[3], clPar5, dclPar[4], clDataSize, clCooling, log_cool[i]);
                edot_met[i] = -std::pow(10.0, log_cool[i]);

                // Ignora il termine del CMB se T >> T_CMB
                if ((icmbTfloor == 1) && ((log10tem[i] - log10_tCMB) < 2.0)) {
                    interpolate_5D_g(log_n_h[i], log_Z[i], log_e_frac[i], zr, log10_tCMB, clGridDim, clPar1, dclPar[0], clPar2, dclPar[1], clPar3, dclPar[2], clPar4, dclPar[3], clPar5, dclPar[4], clDataSize, clCooling, log_cool_cmb[i]);
                    edot_met[i] += std::pow(10.0, log_cool_cmb[i]);
                }

                if (iClHeat == 1) {
                    interpolate_5D_g(log_n_h[i], log_Z[i], log_e_frac[i], zr, log10tem[i], clGridDim, clPar1, dclPar[0], clPar2, dclPar[1], clPar3, dclPar[2], clPar4, dclPar[3], clPar5, dclPar[4], clDataSize, clHeating, log_heat[i]);
                    edot_met[i] += std::pow(10.0, log_heat[i]);
                }

                // Scala il contributo con la frazione di elettroni
                edot_met[i] *= cl_e_frac[i];
            }

            // Aggiorna edot con il contributo del raffreddamento/riscaldamento dei metalli
            edot[i] += edot_met[i] * rhoH[i] * d[i][j][k];
        }
    }
}

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cstdint>  // per tipi interi a 64 bit

// Definizione di costanti
const double DKIND = 1.0;  // Usato come moltiplicatore per double, in Fortran può essere usato per specificare la precisione

void interpolate_1D_g(double x, const std::vector<int64_t>& clGridDim,
                      const std::vector<double>& clPar1, double dclPar1,
                      int64_t clDataSize, const std::vector<double>& clData,
                      double& result);

void interpolate_2D_g(double x, double y, const std::vector<int64_t>& clGridDim,
                      const std::vector<double>& clPar1, double dclPar1,
                      const std::vector<double>& clPar2, double dclPar2,
                      int64_t clDataSize, const std::vector<double>& clData,
                      double& result);

void interpolate_3Dz_g(double x, double z, double y,
                       const std::vector<int64_t>& clGridDim,
                       const std::vector<double>& clPar1, double dclPar1,
                       const std::vector<double>& clPar2, int64_t zindex,
                       const std::vector<double>& clPar3, double dclPar3,
                       int64_t clDataSize, const std::vector<double>& clData,
                       bool end_int, double& result);

void cool1d_cloudy_g(
    const std::vector<std::vector<std::vector<double>>>& d,
    const std::vector<double>& rhoH,
    const std::vector<double>& metallicity,
    int in, int jn, int kn, int is, int ie, int j, int k,
    const std::vector<double>& logtem,
    std::vector<double>& edot,
    double comp2, double dom, double zr,
    int icmbTfloor, int iClHeat, int iZscale,
    int64_t clGridRank,
    const std::vector<int64_t>& clGridDim,
    const std::vector<double>& clPar1,
    const std::vector<double>& clPar2,
    const std::vector<double>& clPar3,
    int64_t clDataSize,
    const std::vector<double>& clCooling,
    const std::vector<double>& clHeating,
    const std::vector<bool>& itmask)
{
    // Variabili locali
    bool end_int = false;
    int get_heat = iClHeat;

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

    // Vettori locali per la slice
    std::vector<double> log_n_h(in);
    std::vector<double> log_cool(in), log_cool_cmb(in), log_heat(in);
    std::vector<double> edot_met(in), log10tem(in);

    for (int i = is; i <= ie; ++i) {
        if (itmask[i]) {
            log10tem[i] = logtem[i] * inv_log10;

            // Calcola log(n_H) in unità "proper"
            log_n_h[i] = std::log10(rhoH[i] * dom);

            int64_t zindex = 0, zmidpt = 0, zhighpt = 0;

            // Calcola l'indice per la dimensione del redshift
            if (clGridRank > 2) {
                if (zr <= clPar2[0]) {
                    zindex = 0;
                } else if (zr >= clPar2[clGridDim[1] - 2]) {
                    zindex = clGridDim[1] - 1;
                    end_int = true;
                    get_heat = 0;
                } else if (zr >= clPar2[clGridDim[1] - 3]) {
                    zindex = clGridDim[1] - 3;
                } else {
                    zindex = 0;
                    zhighpt = clGridDim[1] - 3;
                    while ((zhighpt - zindex) > 1) {
                        zmidpt = (zhighpt + zindex) / 2;
                        if (zr >= clPar2[zmidpt]) {
                            zindex = zmidpt;
                        } else {
                            zhighpt = zmidpt;
                        }
                    }
                }
            }

            // Interpolazione per ottenere il raffreddamento/riscaldamento
            if (clGridRank == 1) {
                interpolate_1D_g(log10tem[i], clGridDim, clPar1, dclPar[0],
                                 clDataSize, clCooling, log_cool[i]);
                edot_met[i] = -std::pow(10.0, log_cool[i]);

                // Ignora il termine del CMB se T >> T_CMB
                if ((icmbTfloor == 1) && ((log10tem[i] - log10_tCMB) < 2.0)) {
                    interpolate_1D_g(log10_tCMB, clGridDim, clPar1, dclPar[0],
                                     clDataSize, clCooling, log_cool_cmb[i]);
                    edot_met[i] += std::pow(10.0, log_cool_cmb[i]);
                }

                if (get_heat == 1) {
                    interpolate_1D_g(log10tem[i], clGridDim, clPar1, dclPar[0],
                                     clDataSize, clHeating, log_heat[i]);
                    edot_met[i] += std::pow(10.0, log_heat[i]);
                }
            } else if (clGridRank == 2) {
                interpolate_2D_g(log_n_h[i], log10tem[i], clGridDim,
                                 clPar1, dclPar[0], clPar2, dclPar[1],
                                 clDataSize, clCooling, log_cool[i]);
                edot_met[i] = -std::pow(10.0, log_cool[i]);

                // Ignora il termine del CMB se T >> T_CMB
                if ((icmbTfloor == 1) && ((log10tem[i] - log10_tCMB) < 2.0)) {
                    interpolate_2D_g(log_n_h[i], log10_tCMB, clGridDim,
                                     clPar1, dclPar[0], clPar2, dclPar[1],
                                     clDataSize, clCooling, log_cool_cmb[i]);
                    edot_met[i] += std::pow(10.0, log_cool_cmb[i]);
                }

                if (get_heat == 1) {
                    interpolate_2D_g(log_n_h[i], log10tem[i], clGridDim,
                                     clPar1, dclPar[0], clPar2, dclPar[1],
                                     clDataSize, clHeating, log_heat[i]);
                    edot_met[i] += std::pow(10.0, log_heat[i]);
                }
            } else if (clGridRank == 3) {
                interpolate_3Dz_g(log_n_h[i], zr, log10tem[i],
                                  clGridDim, clPar1, dclPar[0],
                                  clPar2, zindex, clPar3, dclPar[2],
                                  clDataSize, clCooling, end_int, log_cool[i]);
                edot_met[i] = -std::pow(10.0, log_cool[i]);

                // Ignora il termine del CMB se T >> T_CMB
                if ((icmbTfloor == 1) && ((log10tem[i] - log10_tCMB) < 2.0)) {
                    interpolate_3Dz_g(log_n_h[i], zr, log10_tCMB,
                                      clGridDim, clPar1, dclPar[0],
                                      clPar2, zindex, clPar3, dclPar[2],
                                      clDataSize, clCooling, end_int, log_cool_cmb[i]);
                    edot_met[i] += std::pow(10.0, log_cool_cmb[i]);
                }

                if (get_heat == 1) {
                    interpolate_3Dz_g(log_n_h[i], zr, log10tem[i],
                                      clGridDim, clPar1, dclPar[0],
                                      clPar2, zindex, clPar3, dclPar[2],
                                      clDataSize, clHeating, end_int, log_heat[i]);
                    edot_met[i] += std::pow(10.0, log_heat[i]);
                }
            } else {
                #pragma omp critical
                {
                    std::cerr << "Maximum cooling data grid rank is 3!" << std::endl;
                }
                return;
            }

            // Scala il raffreddamento per la metallicità
            if (iZscale == 1) {
                edot_met[i] *= metallicity[i];
            }

            edot[i] += (edot_met[i] * rhoH[i] * rhoH[i]);
        }
    }
}

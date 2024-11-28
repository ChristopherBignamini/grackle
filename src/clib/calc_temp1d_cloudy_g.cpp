#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cstdint>  // per tipi interi a 64 bit
#include <omp.h>    // per OpenMP, se necessario

const double mu_metal = 16.0;  // Peso molecolare medio approssimativo dei metalli
const int ti_max = 20;         // Numero massimo di iterazioni
const double DKIND = 1.0;      // Costante per la precisione (in Fortran)

void calc_temp1d_cloudy_g(
    const std::vector<std::vector<std::vector<double>>>& d,
    const std::vector<std::vector<std::vector<double>>>& metal,
    const std::vector<std::vector<std::vector<double>>>& e,
    const std::vector<double>& rhoH,
    int in, int jn, int kn, int is, int ie, int j, int k,
    std::vector<double>& tgas,
    std::vector<double>& mmw,
    double dom, double zr,
    double temstart, double temend,
    double gamma, double utem, int imetal,
    int64_t clGridRank,
    const std::vector<int64_t>& clGridDim,
    const std::vector<double>& clPar1,
    const std::vector<double>& clPar2,
    const std::vector<double>& clPar3,
    int64_t clDataSize,
    const std::vector<double>& clMMW,
    std::vector<bool>& itmask)
{
    // Variabili locali
    bool end_int = false;
    double inv_log10 = 1.0 / std::log(10.0);

    // Calcola le pendenze dei parametri
    std::vector<double> dclPar(clGridRank);
    dclPar[0] = (clPar1[clGridDim[0] - 1] - clPar1[0]) / static_cast<double>(clGridDim[0] - 1);
    if (clGridRank > 1) {
        dclPar[1] = (clPar2[clGridDim[1] - 1] - clPar2[0]) / static_cast<double>(clGridDim[1] - 1);
    }
    if (clGridRank > 2) {
        dclPar[2] = (clPar3[clGridDim[2] - 1] - clPar3[0]) / static_cast<double>(clGridDim[2] - 1);
    }

    // Calcola l'indice per la dimensione del redshift
    int64_t zindex = 0;
    if (clGridRank > 2) {
        if (zr <= clPar2[0]) {
            zindex = 0;
        } else if (zr >= clPar2[clGridDim[1] - 2]) {
            zindex = clGridDim[1] - 1;
            end_int = true;
        } else if (zr >= clPar2[clGridDim[1] - 3]) {
            zindex = clGridDim[1] - 3;
        } else {
            zindex = 0;
            int64_t zhighpt = clGridDim[1] - 3;
            while ((zhighpt - zindex) > 1) {
                int64_t zmidpt = (zhighpt + zindex) / 2;
                if (zr >= clPar2[zmidpt]) {
                    zindex = zmidpt;
                } else {
                    zhighpt = zmidpt;
                }
            }
        }
    }

    // Vettori locali per la slice
    std::vector<double> log_n_h(in);
    std::vector<double> log10tem(in);
    std::vector<double> logtem(in);

    // Calcolo di log(n_H) per ogni elemento
    for (int i = is; i <= ie; ++i) {
        if (itmask[i]) {
            log_n_h[i] = std::log10(rhoH[i] * dom);
        }
    }

    // Iterazione principale per calcolare tgas e mmw
    for (int i = is; i <= ie; ++i) {
        if (itmask[i]) {
            double munew = 1.0;
            double muold = 0.0;

            for (int ti = 0; ti < ti_max; ++ti) {
                muold = munew;

                // Calcolo della temperatura
                tgas[i] = std::max((gamma - 1.0) * e[i][j][k] * munew * utem, temstart);
                logtem[i] = std::log(tgas[i]);
                log10tem[i] = logtem[i] * inv_log10;

                // Interpolazione per ottenere munew
                if (clGridRank == 1) {
                    interpolate_1D_g(log10tem[i], clGridDim, clPar1, dclPar[0],
                                     clDataSize, clMMW, munew);
                } else if (clGridRank == 2) {
                    interpolate_2D_g(log_n_h[i], log10tem[i], clGridDim,
                                     clPar1, dclPar[0], clPar2, dclPar[1],
                                     clDataSize, clMMW, munew);
                } else if (clGridRank == 3) {
                    interpolate_3Dz_g(log_n_h[i], zr, log10tem[i],
                                      clGridDim, clPar1, dclPar[0],
                                      clPar2, zindex, clPar3, dclPar[2],
                                      clDataSize, clMMW, end_int, munew);
                } else {
                    #pragma omp critical
                    {
                        std::cerr << "Maximum mmw data grid rank is 3!" << std::endl;
                    }
                    return;
                }

                // Aggiornamento di munew e tgas
                munew = 0.5 * (munew + muold);
                tgas[i] = tgas[i] * munew / muold;

                // Controllo della convergenza
                if (std::abs((munew / muold) - 1.0) <= 1.e-2) {
                    muold = munew;

                    // Aggiunta dell'effetto dei metalli se necessario
                    if (imetal == 1) {
                        munew = d[i][j][k] /
                                ((d[i][j][k] - metal[i][j][k]) / munew +
                                 metal[i][j][k] / mu_metal);
                        tgas[i] = tgas[i] * munew / muold;
                    }

                    mmw[i] = munew;
                    break;  // Esce dal ciclo di iterazione
                }
            }

            // Se non convergente dopo ti_max iterazioni
            if (std::abs((munew / muold) - 1.0) > 1.e-2) {
                mmw[i] = munew;
                #pragma omp critical
                {
                    std::cerr << "Mean molecular weight not converged! "
                              << "munew: " << munew << ", muold: " << muold
                              << ", difference: " << std::abs((munew / muold) - 1.0) << std::endl;
                }
            }
        }
    }
}

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdint>  // per tipi interi a 64 bit
#include <omp.h>    // per OpenMP, se necessario

// Definizione delle costanti fisiche
const double mass_h = 1.6735575e-24;  // Massa dell'idrogeno in grammi (CGS)
const double mh = mass_h;

// Altre costanti
const double DKIND = 1.0;  // Usato per coerenza con il codice Fortran

void calc_temp_cloudy_g(
    std::vector<std::vector<std::vector<double>>>& d,
    std::vector<std::vector<std::vector<double>>>& e,
    std::vector<std::vector<std::vector<double>>>& metal,
    std::vector<std::vector<std::vector<double>>>& temperature,
    int in, int jn, int kn, int iexpand, int imetal,
    int is, int js, int ks, int ie, int je, int ke,
    double aye, double temstart, double temend,
    double utem, double uxyz, double uaye, double urho, double utim,
    double gamma, double fh,
    int64_t priGridRank, const std::vector<int64_t>& priGridDim,
    const std::vector<double>& priPar1, const std::vector<double>& priPar2,
    const std::vector<double>& priPar3, int64_t priDataSize,
    const std::vector<double>& priMMW)
{
    // Variabili locali
    double dom, zr;
    double dbase1, tbase1, xbase1;

    // Vettori temporanei per una slice
    std::vector<double> tgas(in);
    std::vector<double> rhoH(in);
    std::vector<double> mmw(in);
    std::vector<bool> itmask(in, true);

    // Calcolo delle unità
    dom = urho * std::pow(aye, 3) / mh;
    tbase1 = utim;
    xbase1 = uxyz / (aye * uaye);
    dbase1 = urho * std::pow(aye * uaye, 3);
    zr = 1.0 / (aye * uaye) - 1.0;

    // Converti le densità da comoventi a proprie, se necessario
    if (iexpand == 1) {
        scale_fields_table_g(d, metal, is, ie, js, je, ks, ke, in, jn, kn, imetal, std::pow(aye, -3));
    }

    // Loop su k e j
    int dj = je - js + 1;
    int dk = ke - ks + 1;

    #pragma omp parallel for schedule(runtime) private(tgas, rhoH, mmw, itmask)
    for (int t = 0; t < dk * dj; ++t) {
        int k = t / dj + ks;
        int j = t % dj + js;

        // Inizializza la maschera di iterazione e calcola rhoH
        for (int i = is; i <= ie; ++i) {
            itmask[i] = true;
            rhoH[i] = fh * d[i][j][k];
        }

        // Calcola la temperatura e il peso molecolare medio
        calc_temp1d_cloudy_g(d, metal, e, rhoH, in, jn, kn, is, ie, j, k,
                             tgas, mmw, dom, zr,
                             temstart, temend,
                             gamma, utem, imetal,
                             priGridRank, priGridDim,
                             priPar1, priPar2, priPar3,
                             priDataSize, priMMW,
                             itmask);

        // Copia i valori calcolati nell'array temperature
        for (int i = is; i <= ie; ++i) {
            temperature[i][j][k] = tgas[i];
        }
    }

    // Converti le densità da proprie a comoventi, se necessario
    if (iexpand == 1) {
        scale_fields_table_g(d, metal, is, ie, js, je, ks, ke, in, jn, kn, imetal, std::pow(aye, 3));
    }
}

void scale_fields_table_g(
    std::vector<std::vector<std::vector<double>>>& d,
    std::vector<std::vector<std::vector<double>>>& metal,
    int is, int ie, int js, int je, int ks, int ke,
    int in, int jn, int kn, int imetal, double factor)
{
    for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
            for (int i = is; i <= ie; ++i) {
                d[i][j][k] *= factor;
            }
        }
    }

    if (imetal == 1) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    metal[i][j][k] *= factor;
                }
            }
        }
    }
}

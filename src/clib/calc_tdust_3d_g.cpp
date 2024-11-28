#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <omp.h>  // Per OpenMP

// Definizione delle costanti fisiche
const double mass_h = 1.6735575e-24;  // Massa dell'atomo di idrogeno in grammi
const double mh = mass_h;

void calc_tdust_3d_g(
    const std::vector<std::vector<std::vector<double>>>& d,
    const std::vector<std::vector<std::vector<double>>>& de,
    const std::vector<std::vector<std::vector<double>>>& HI,
    const std::vector<std::vector<std::vector<double>>>& HII,
    const std::vector<std::vector<std::vector<double>>>& HeI,
    const std::vector<std::vector<std::vector<double>>>& HeII,
    const std::vector<std::vector<std::vector<double>>>& HeIII,
    const std::vector<std::vector<std::vector<double>>>& HM,
    const std::vector<std::vector<std::vector<double>>>& H2I,
    const std::vector<std::vector<std::vector<double>>>& H2II,
    int in, int jn, int kn,
    int nratec, int iexpand,
    int ispecies, int idim,
    int is, int js, int ks,
    int ie, int je, int ke,
    double aye, double temstart, double temend,
    double fgr, const std::vector<double>& gasgra,
    double gamma_isrfa, double isrf,
    double utem, double uxyz, double uaye,
    double urho, double utim,
    std::vector<std::vector<std::vector<double>>>& gas_temp,
    std::vector<std::vector<std::vector<double>>>& dust_temp,
    int iisrffield,
    const std::vector<std::vector<std::vector<double>>>& isrf_habing) {

    // Definizione delle costanti
    const double mh = mass_h;

    // Variabili locali
    double trad, zr, logtem0, logtem9, dlogtem;
    double coolunit, dbase1, tbase1, xbase1;

    // Calcolo dei valori logaritmici per le tabelle
    logtem0 = std::log(temstart);
    logtem9 = std::log(temend);
    dlogtem = (std::log(temend) - std::log(temstart)) / static_cast<double>(nratec - 1);

    // Impostazione delle unità
    tbase1 = utim;
    xbase1 = uxyz / (aye * uaye);
    dbase1 = urho * std::pow((aye * uaye), 3);
    coolunit = (std::pow(uaye, 5) * std::pow(xbase1, 2) * std::pow(mh, 2)) / (std::pow(tbase1, 3) * dbase1);
    zr = 1.0 / (aye * uaye) - 1.0;

    // Impostazione della temperatura del CMB
    trad = 2.73 * (1.0 + zr);

    // Dimensioni locali per i loop
    int dk = ke - ks + 1;
    int dj = je - js + 1;

    // Variabili per il parallelismo
    #pragma omp parallel for schedule(runtime) collapse(2)
    for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
            // Vettori locali per il calcolo in una slice
            std::vector<double> indixe(in);
            std::vector<double> t1(in), t2(in), logtem(in), tdef(in);
            std::vector<double> tgas(in), tdust(in), nh(in), gasgr_local(in);
            std::vector<double> myisrf(in);
            std::vector<bool> itmask(in, true);

            for (int i = is; i <= ie; ++i) {
                // Calcolo del campo di radiazione interstellare
                if (iisrffield > 0) {
                    myisrf[i] = isrf_habing[i][j][k];
                } else {
                    myisrf[i] = isrf;
                }

                // Calcolo della densità numerica dell'idrogeno
                nh[i] = HI[i][j][k] + HII[i][j][k];
                if (ispecies > 1) {
                    nh[i] += H2I[i][j][k] + H2II[i][j][k];
                }

                // Non abbiamo convertito a densità "proper", quindi usiamo urho
                nh[i] = nh[i] * urho / mh;

                // Calcolo del logaritmo della temperatura del gas
                tgas[i] = gas_temp[i][j][k];
                logtem[i] = std::log(tgas[i]);
                logtem[i] = std::max(logtem[i], logtem0);
                logtem[i] = std::min(logtem[i], logtem9);

                // Calcolo dell'indice nella tabella e interpolazione lineare
                indixe[i] = std::min(static_cast<double>(nratec - 1),
                                     std::max(1.0, std::floor((logtem[i] - logtem0) / dlogtem) + 1.0));
                t1[i] = logtem0 + (indixe[i] - 1) * dlogtem;
                t2[i] = logtem0 + indixe[i] * dlogtem;
                tdef[i] = (logtem[i] - t1[i]) / (t2[i] - t1[i]);

                // Interpolazione e conversione in unità cgs
                int idx = static_cast<int>(indixe[i]);
                gasgr_local[i] = gasgra[idx] + tdef[i] * (gasgra[idx + 1] - gasgra[idx]);
                gasgr_local[i] = gasgr_local[i] * fgr * coolunit / mh;
            }

            // Calcolo della temperatura della polvere in una slice
            calc_tdust_1d_g(tdust, tgas, nh, gasgr_local, gamma_isrfa, myisrf, itmask, trad, in, is, ie, j, k);

            // Copia dei valori calcolati nella griglia 3D
            for (int i = is; i <= ie; ++i) {
                dust_temp[i][j][k] = tdust[i];
            }
        }
    }
}

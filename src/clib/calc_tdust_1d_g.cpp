#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

const double t_subl = 1.5e3;  // grain sublimation temperature
const double radf = 4.0 * sigma_sb;
const double kgr1 = 4.0e-4;
const double tol = 1.e-5;
const double bi_tol = 1.e-3;
const double minpert = 1.e-10;
const int itmax = 50;
const int bi_itmax = 30;

void calc_kappa_gr_g(const std::vector<double>& tdust, std::vector<double>& kgr,
                     const std::vector<bool>& itmask, int in, int is, int ie, double t_subl);

void calc_gr_balance_g(const std::vector<double>& tdust, const std::vector<double>& tgas,
                       const std::vector<double>& kgr, double trad4, const std::vector<double>& gasgr,
                       const std::vector<double>& gamma_isrf, const std::vector<double>& nh,
                       const std::vector<bool>& itmask, std::vector<double>& sol, int in, int is, int ie);

void calc_tdust_1d_g(std::vector<double>& tdust, const std::vector<double>& tgas,
                     const std::vector<double>& nh, const std::vector<double>& gasgr,
                     double gamma_isrfa, std::vector<double>& isrf, const std::vector<bool>& itmask,
                     double trad, int in, int is, int ie, int j, int k) {
    double pert_i = 1.e-3;
    trad = std::max(1.0, trad);
    double trad4 = std::pow(trad, 4);

    int c_done = 0;
    int nm_done = 0;
    int c_total = ie - is + 1;

    std::vector<double> gamma_isrf(in);
    for (int i = is; i <= ie; ++i) {
        if (itmask[i]) {
            gamma_isrf[i] = isrf[i] * gamma_isrfa;
        }
    }

    std::vector<bool> nm_itmask = itmask;
    std::vector<bool> bi_itmask = itmask;
    std::vector<double> tdustnow(in);
    std::vector<double> pert(in, pert_i);

    for (int i = is; i <= ie; ++i) {
        if (nm_itmask[i]) {
            if (trad >= tgas[i]) {
                tdustnow[i] = trad;
                nm_itmask[i] = false;
                bi_itmask[i] = false;
                c_done++;
                nm_done++;
            } else if (tgas[i] > t_subl) {
                nm_itmask[i] = false;
                nm_done++;
            } else {
                tdustnow[i] = std::max(trad, std::pow(gamma_isrf[i] / (radf * kgr1), 0.17));
            }
        } else {
            c_done++;
            nm_done++;
        }
    }

    std::vector<double> tdplus(in), kgr(in), kgrplus(in), sol(in), solplus(in);
    std::vector<double> slope(in), tdustold(in);
    for (int iter = 1; iter <= itmax; ++iter) {
        for (int i = is; i <= ie; ++i) {
            if (nm_itmask[i]) {
                tdplus[i] = std::max(1.e-3, (1.0 + pert[i]) * tdustnow[i]);
            }
        }

        calc_kappa_gr_g(tdustnow, kgr, nm_itmask, in, is, ie, t_subl);
        calc_kappa_gr_g(tdplus, kgrplus, nm_itmask, in, is, ie, t_subl);

        calc_gr_balance_g(tdustnow, tgas, kgr, trad4, gasgr, gamma_isrf, nh, nm_itmask, sol, in, is, ie);
        calc_gr_balance_g(tdplus, tgas, kgrplus, trad4, gasgr, gamma_isrf, nh, nm_itmask, solplus, in, is, ie);

        for (int i = is; i <= ie; ++i) {
            if (nm_itmask[i]) {
                slope[i] = (solplus[i] - sol[i]) / (pert[i] * tdustnow[i]);
                tdustold[i] = tdustnow[i];
                tdustnow[i] = tdustnow[i] - (sol[i] / slope[i]);

                pert[i] = std::max(std::min(pert[i],
                    0.5 * std::abs(tdustnow[i] - tdustold[i]) / tdustnow[i]), minpert);

                if (tdustnow[i] < trad) {
                    nm_itmask[i] = false;
                    nm_done++;
                } else if (std::abs(sol[i] / solplus[i]) < tol) {
                    nm_itmask[i] = false;
                    c_done++;
                    bi_itmask[i] = false;
                    nm_done++;
                }
            }
        }

        if (c_done >= c_total) break;
        if (nm_done >= c_total) break;
    }

    if (c_done < c_total) {
        std::vector<double> bi_t_mid(in), bi_t_high(in);
        for (int i = is; i <= ie; ++i) {
            if (bi_itmask[i]) {
                tdustnow[i] = trad;
                bi_t_high[i] = tgas[i];
            }
        }

        for (int iter = 1; iter <= bi_itmax; ++iter) {
            for (int i = is; i <= ie; ++i) {
                if (bi_itmask[i]) {
                    bi_t_mid[i] = 0.5 * (tdustnow[i] + bi_t_high[i]);
                    if (iter == 1) {
                        bi_t_mid[i] = std::min(bi_t_mid[i], t_subl);
                    }
                }
            }

            calc_kappa_gr_g(bi_t_mid, kgr, bi_itmask, in, is, ie, t_subl);
            calc_gr_balance_g(bi_t_mid, tgas, kgr, trad4, gasgr, gamma_isrf, nh, bi_itmask, sol, in, is, ie);

            for (int i = is; i <= ie; ++i) {
                if (bi_itmask[i]) {
                    if (sol[i] > 0.0) {
                        tdustnow[i] = bi_t_mid[i];
                    } else {
                        bi_t_high[i] = bi_t_mid[i];
                    }

                    if ((std::abs(bi_t_high[i] - tdustnow[i]) / tdustnow[i]) <= bi_tol) {
                        bi_itmask[i] = false;
                        c_done++;
                    }

                    if (c_done >= c_total) break;
                }
            }
            if (iter > itmax) {
                std::cerr << "CALC_TDUST_1D_G failed using both methods for " << (c_total - c_done) << " cells.\n";
            }
        }
    }

    for (int i = is; i <= ie; ++i) {
        if (itmask[i]) {
            if (tdustnow[i] < 0.0) {
                std::cerr << "CALC_TDUST_1D_G Newton method - T_dust < 0: i = " << i << " j = " << j
                          << " k = " << k << " nh = " << nh[i] << " t_gas = " << tgas[i]
                          << " t_rad = " << trad << " t_dust = " << tdustnow[i] << "\n";
            }
            tdust[i] = tdustnow[i];
        }
    }
}

void calc_kappa_gr_g(const std::vector<double>& tdust, std::vector<double>& kgr,
                     const std::vector<bool>& itmask, int in, int is, int ie, double t_subl) {
    const double kgr1 = 4.0e-4;
    const double kgr200 = 16.0;

    for (int i = is; i <= ie; ++i) {
        if (itmask[i]) {
            if (tdust[i] < 200.0) {
                kgr[i] = kgr1 * tdust[i] * tdust[i];
            } else if (tdust[i] < t_subl) {
                kgr[i] = kgr200;
            } else {
                kgr[i] = std::max(std::numeric_limits<double>::min(),
                                  kgr200 * std::pow(tdust[i] / 1.5e3, -12));
            }
        }
    }
}

void calc_gr_balance_g(const std::vector<double>& tdust, const std::vector<double>& tgas,
                       const std::vector<double>& kgr, double trad4, const std::vector<double>& gasgr,
                       const std::vector<double>& gamma_isrf, const std::vector<double>& nh,
                       const std::vector<bool>& itmask, std::vector<double>& sol, int in, int is, int ie) {
    const double radf = 4.0 * sigma_sb;

    for (int i = is; i <= ie; ++i) {
        if (itmask[i]) {
            sol[i] = gamma_isrf[i] + radf * kgr[i] * (trad4 - std::pow(tdust[i], 4))
                     + (gasgr[i] * nh[i] * (tgas[i] - tdust[i]));
        }
    }
}

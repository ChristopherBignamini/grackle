#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdint>  // Per tipi interi a 64 bit
#include <omp.h>    // Per OpenMP, se necessario

// Definizione delle costanti fisiche
const double mass_h = 1.6735575e-24;  // Massa dell'idrogeno in grammi (CGS)
const double mh = mass_h;
const double mu_metal = 16.0;  // Peso molecolare medio dei metalli (approssimativo)
const int ti_max = 20;         // Numero massimo di iterazioni

// Costanti aggiuntive
const double DKIND = 1.0;  // Per coerenza con il codice Fortran
const double DIKIND = 1.0; // Per gli indici interi
const double tiny = 1e-20; // Un valore molto piccolo per evitare divisioni per zero
const double tiny8 = 1e-8;

// Macro per gestire le opzioni di compilazione del Fortran
#define USE_GLOVER_ABEL2008
#define OPTICAL_DEPTH_FUDGE

void cool1d_multi_g(
    // Array 3D di densità, energia interna, velocità, ecc.
    const std::vector<std::vector<std::vector<double>>>& d,
    const std::vector<std::vector<std::vector<double>>>& e,
    const std::vector<std::vector<std::vector<double>>>& u,
    const std::vector<std::vector<std::vector<double>>>& v,
    const std::vector<std::vector<std::vector<double>>>& w,
    const std::vector<std::vector<std::vector<double>>>& de,
    const std::vector<std::vector<std::vector<double>>>& HI,
    const std::vector<std::vector<std::vector<double>>>& HII,
    const std::vector<std::vector<std::vector<double>>>& HeI,
    const std::vector<std::vector<std::vector<double>>>& HeII,
    const std::vector<std::vector<std::vector<double>>>& HeIII,
    int in, int jn, int kn, int nratec,
    int iexpand, int ispecies, int imetal, int imcool,
    int idust, int idustall, int idustfield, int idustrec,
    int idim, int is, int ie, int j, int k, int ih2co, int ipiht, int iter, int igammah,
    double aye, double temstart, double temend, double z_solar, double fgr,
    double utem, double uxyz, double uaye, double urho, double utim,
    double gamma_val, double fh,
    // Array di parametri e tabelle di raffreddamento
    const std::vector<double>& ceHIa, const std::vector<double>& ceHeIa, const std::vector<double>& ceHeIIa,
    const std::vector<double>& ciHIa, const std::vector<double>& ciHeIa,
    const std::vector<double>& ciHeISa, const std::vector<double>& ciHeIIa,
    const std::vector<double>& reHIIa, const std::vector<double>& reHeII1a,
    const std::vector<double>& reHeII2a, const std::vector<double>& reHeIIIa, const std::vector<double>& brema,
    double compa, double gammaha,
    double isrf, const std::vector<double>& regra, const std::vector<double>& gamma_isrfa,
    double comp_xraya, double comp_temp,
    double piHI, double piHeI, double piHeII,
    // Altri array di specie chimiche
    const std::vector<std::vector<std::vector<double>>>& HM,
    const std::vector<std::vector<std::vector<double>>>& H2I,
    const std::vector<std::vector<std::vector<double>>>& H2II,
    const std::vector<std::vector<std::vector<double>>>& DI,
    const std::vector<std::vector<std::vector<double>>>& DII,
    const std::vector<std::vector<std::vector<double>>>& HDI,
    const std::vector<std::vector<std::vector<double>>>& metal,
    const std::vector<std::vector<std::vector<double>>>& dust,
    // Array per le tabelle di raffreddamento H2
    const std::vector<double>& hyd01ka, const std::vector<double>& h2k01a, const std::vector<double>& vibha,
    const std::vector<double>& rotha, const std::vector<double>& rotla,
    const std::vector<double>& gpldla, const std::vector<double>& gphdla,
    const std::vector<double>& hdltea, const std::vector<double>& hdlowa,
    const std::vector<double>& gaHIa, const std::vector<double>& gaH2a, const std::vector<double>& gaHea, const std::vector<double>& gaHpa, const std::vector<double>& gaela,
    const std::vector<double>& h2ltea, const std::vector<double>& gasgra,
    // Parametri per il calcolo del raffreddamento
    const std::vector<double>& ciecoa,
    // Parametri per le griglie Cloudy
    int64_t priGridRank, const std::vector<int64_t>& priGridDim,
    const std::vector<double>& priPar1, const std::vector<double>& priPar2, const std::vector<double>& priPar3, const std::vector<double>& priPar4, const std::vector<double>& priPar5,
    int64_t priDataSize, const std::vector<double>& priCooling, const std::vector<double>& priHeating, const std::vector<double>& priMMW,
    int64_t metGridRank, const std::vector<int64_t>& metGridDim,
    const std::vector<double>& metPar1, const std::vector<double>& metPar2, const std::vector<double>& metPar3, const std::vector<double>& metPar4, const std::vector<double>& metPar5,
    int64_t metDataSize, const std::vector<double>& metCooling, const std::vector<double>& metHeating, int clnew,
    int icmbTfloor, int iClHeat, double clEleFra,
    // Altri parametri
    int iVheat, int iMheat, const std::vector<std::vector<std::vector<double>>>& Vheat, const std::vector<std::vector<std::vector<double>>>& Mheat,
    int iTfloor, double Tfloor_scalar, const std::vector<std::vector<std::vector<double>>>& Tfloor,
    int iisrffield, const std::vector<std::vector<std::vector<double>>>& isrf_habing,
    int iradshield, double avgsighi, double avgsighei,
    double avgsigheii,
    double k24, double k26, int iradtrans, const std::vector<std::vector<std::vector<double>>>& photogamma,
    int ih2optical, int iciecool,
    int is, int ie, int j, int k, int iter,
    std::vector<double>& tgas, std::vector<double>& tgasold, std::vector<double>& edot,
    std::vector<bool>& itmask)
void cool1d_multi_g(
    // Tutti i parametri come definiti sopra
)
{
    // Variabili locali
    double dom, qq, vibl, logtem0, logtem9, dlogtem, zr;
    double comp1_local, comp2_local;
    double hdlte1, hdlow1, gamma2, x, fudge, fH2;
    double gphdl1, dom_inv, tau, ciefudge;
    double coolunit, dbase1, tbase1, xbase1;
    double nH2, nother, nSSh, nratio, nSSh_he, nratio_he;
    double fSShHI, fSShHeI, pe_eps, pe_X, grbeta;
    bool anydust, interp;
    int iradfield, mycmbTfloor;
    int i;
    int64_t ind;
    double nratio_h2;

    // Vettori locali per le slice
    int size = ie - is + 1;
    std::vector<double> t1(size), t2(size), logtem(size), tdef(size), p2d(size);
    std::vector<double> mmw(size), tdust(size), rhoH(size);
    std::vector<double> mynh(size), metallicity(size), dust2gas(size);
    std::vector<double> myde(size), gammaha_eff(size), gasgr_tdust(size), regr(size), myisrf(size);
    std::vector<double> ceHI(size), ceHeI(size), ceHeII(size), ciHI(size), ciHeI(size), ciHeIS(size), ciHeII(size);
    std::vector<double> reHII(size), reHeII1(size), reHeII2(size), reHeIII(size), brem(size), cieco(size);
    std::vector<double> hyd01k(size), h2k01(size), vibh(size), roth(size), rotl(size);
    std::vector<double> gpldl(size), gphdl(size), hdlte(size), hdlow(size);
    std::vector<double> gaHI(size), gaH2(size), gaHe(size), gaHp(size), gael(size);
    std::vector<double> h2lte(size), galdl(size);
    std::vector<int64_t> indixe(size);

    // Impostazione dei flag
    anydust = (idust > 0) || (idustall > 0) || (idustrec > 0);
    interp = (ispecies > 0) || (idustall > 0);

    // Impostazione dei valori di logaritmo
    logtem0 = std::log(temstart);
    logtem9 = std::log(temend);
    dlogtem = (std::log(temend) - std::log(temstart)) / static_cast<double>(nratec - 1);

    // Impostazione delle unità
    dom = urho * std::pow(aye, 3) / mh;
    dom_inv = 1.0 / dom;
    tbase1 = utim;
    xbase1 = uxyz / (aye * uaye);
    dbase1 = urho * std::pow(aye * uaye, 3);
    coolunit = (std::pow(uaye, 5) * std::pow(xbase1, 2) * std::pow(mh, 2)) / (std::pow(tbase1, 3) * dbase1);
    zr = 1.0 / (aye * uaye) - 1.0;
    fudge = 1.0;
    iradfield = -1;

    // Coefficienti di raffreddamento Compton
    comp1_local = compa * std::pow(1.0 + zr, 4);
    comp2_local = 2.73 * (1.0 + zr);

    // Inizializzazione di edot
    for (i = is; i <= ie; ++i) {
        if (itmask[i]) {
            edot[i] = 0.0;
        }
    }

    // Calcolo della pressione
    for (i = is; i <= ie; ++i) {
        if (itmask[i]) {
            p2d[i - is] = (gamma_val - 1.0) * d[i][j][k] * e[i][j][k];
        }
    }

    // Calcolo della temperatura
    if (ispecies == 0) {
        // Calcolo di rhoH
        if (imetal == 1) {
            for (i = is; i <= ie; ++i) {
                if (itmask[i]) {
                    rhoH[i - is] = fh * (d[i][j][k] - metal[i][j][k]);
                }
            }
        } else {
            for (i = is; i <= ie; ++i) {
                if (itmask[i]) {
                    rhoH[i - is] = fh * d[i][j][k];
                }
            }
        }

        // Chiamata alla funzione calc_temp1d_cloudy_g (da implementare)
        // Dovrai implementare questa funzione in base alle tue esigenze
        // Ad esempio:
        // calc_temp1d_cloudy_g(...);
    } else {
        // Calcolo del peso molecolare medio
        for (i = is; i <= ie; ++i) {
            if (itmask[i]) {
                mmw[i - is] = (HeI[i][j][k] + HeII[i][j][k] + HeIII[i][j][k]) / 4.0 +
                              HI[i][j][k] + HII[i][j][k] + de[i][j][k];
                rhoH[i - is] = HI[i][j][k] + HII[i][j][k];
                myde[i - is] = de[i][j][k];
            }
        }

        // Include H2
        if (ispecies > 1) {
            for (i = is; i <= ie; ++i) {
                if (itmask[i]) {
                    mmw[i - is] += HM[i][j][k] + (H2I[i][j][k] + H2II[i][j][k]) / 2.0;
                    rhoH[i - is] += H2I[i][j][k] + H2II[i][j][k];
                }
            }
        }

        // Include metalli
        if (imetal == 1) {
            for (i = is; i <= ie; ++i) {
                if (itmask[i]) {
                    mmw[i - is] += metal[i][j][k] / mu_metal;
                }
            }
        }

        // Calcolo della temperatura
        for (i = is; i <= ie; ++i) {
            if (itmask[i]) {
                tgas[i] = std::max(p2d[i - is] * utem / mmw[i - is], temstart);
                mmw[i - is] = d[i][j][k] / mmw[i - is];
            }
        }

        // Correzione per gamma da H2
        if (ispecies > 1) {
            for (i = is; i <= ie; ++i) {
                if (itmask[i]) {
                    nH2 = 0.5 * (H2I[i][j][k] + H2II[i][j][k]);
                    nother = (HeI[i][j][k] + HeII[i][j][k] + HeIII[i][j][k]) / 4.0 +
                             HI[i][j][k] + HII[i][j][k] + de[i][j][k];
                    nratio_h2 = nH2 / nother;
                    if (nratio_h2 > 1.0e-3) {
                        x = 6100.0 / tgas[i];
                        if (x > 10.0) {
                            gamma2 = 0.5 * 5.0;
                        } else {
                            gamma2 = 0.5 * (5.0 + 2.0 * x * x * std::exp(x) / std::pow(std::exp(x) - 1.0, 2));
                        }
                    } else {
                        gamma2 = 2.5;
                    }
                    gamma2 = 1.0 + (nH2 + nother) / (nH2 * gamma2 + nother / (gamma_val - 1.0));
                    tgas[i] = tgas[i] * (gamma2 - 1.0) / (gamma_val - 1.0);
                }
            }
        }
    }

    // Controllo della temperatura minima
    if (iTfloor == 1) {
        for (i = is; i <= ie; ++i) {
            if (itmask[i]) {
                if (tgas[i] <= Tfloor_scalar) {
                    edot[i] = tiny;
                    itmask[i] = false;
                }
            }
        }
    } else if (iTfloor == 2) {
        for (i = is; i <= ie; ++i) {
            if (itmask[i]) {
                if (tgas[i] <= Tfloor[i][j][k]) {
                    edot[i] = tiny;
                    itmask[i] = false;
                }
            }
        }
    }

    // Calcolo della metallicità e densità di idrogeno
    if (imetal == 1) {
        for (i = is; i <= ie; ++i) {
            if (itmask[i]) {
                metallicity[i - is] = metal[i][j][k] / d[i][j][k] / z_solar;
            }
        }
    }

    for (i = is; i <= ie; ++i) {
        if (itmask[i]) {
            mynh[i - is] = rhoH[i - is] * dom;
        }
    }

    // Inizializzazione di tgasold
    if (iter == 1) {
        for (i = is; i <= ie; ++i) {
            if (itmask[i]) {
                tgasold[i] = tgas[i];
            }
        }
    }

    // Calcolo del logaritmo della temperatura
    for (i = is; i <= ie; ++i) {
        if (itmask[i]) {
            logtem[i - is] = std::log(0.5 * (tgas[i] + tgasold[i]));
            logtem[i - is] = std::max(logtem[i - is], logtem0);
            logtem[i - is] = std::min(logtem[i - is], logtem9);
        }
    }

    // Calcolo degli indici di interpolazione
    if (interp) {
        for (i = is; i <= ie; ++i) {
            if (itmask[i]) {
                ind = static_cast<int64_t>((logtem[i - is] - logtem0) / dlogtem) + 1;
                indixe[i - is] = std::min(static_cast<int64_t>(nratec - 1), std::max(static_cast<int64_t>(1), ind));
                t1[i - is] = logtem0 + (indixe[i - is] - 1) * dlogtem;
                t2[i - is] = logtem0 + indixe[i - is] * dlogtem;
                tdef[i - is] = (logtem[i - is] - t1[i - is]) / (t2[i - is] - t1[i - is]);
            }
        }
    }

    // --- Raffreddamento per le 6 specie ---
    if (ispecies > 0) {
        // Lookup dei valori di raffreddamento e interpolazione
        for (i = is; i <= ie; ++i) {
            if (itmask[i]) {
                ind = indixe[i - is] - 1; // Fortran indexing starts from 1
                ceHI[i - is] = ceHIa[ind] + tdef[i - is] * (ceHIa[ind + 1] - ceHIa[ind]);
                ceHeI[i - is] = ceHeIa[ind] + tdef[i - is] * (ceHeIa[ind + 1] - ceHeIa[ind]);
                ceHeII[i - is] = ceHeIIa[ind] + tdef[i - is] * (ceHeIIa[ind + 1] - ceHeIIa[ind]);
                ciHI[i - is] = ciHIa[ind] + tdef[i - is] * (ciHIa[ind + 1] - ciHIa[ind]);
                ciHeI[i - is] = ciHeIa[ind] + tdef[i - is] * (ciHeIa[ind + 1] - ciHeIa[ind]);
                ciHeIS[i - is] = ciHeISa[ind] + tdef[i - is] * (ciHeISa[ind + 1] - ciHeISa[ind]);
                ciHeII[i - is] = ciHeIIa[ind] + tdef[i - is] * (ciHeIIa[ind + 1] - ciHeIIa[ind]);
                reHII[i - is] = reHIIa[ind] + tdef[i - is] * (reHIIa[ind + 1] - reHIIa[ind]);
                reHeII1[i - is] = reHeII1a[ind] + tdef[i - is] * (reHeII1a[ind + 1] - reHeII1a[ind]);
                reHeII2[i - is] = reHeII2a[ind] + tdef[i - is] * (reHeII2a[ind + 1] - reHeII2a[ind]);
                reHeIII[i - is] = reHeIIIa[ind] + tdef[i - is] * (reHeIIIa[ind + 1] - reHeIIIa[ind]);
                brem[i - is] = brema[ind] + tdef[i - is] * (brema[ind + 1] - brema[ind]);
            }
        }

        // Calcolo della funzione di raffreddamento
        for (i = is; i <= ie; ++i) {
            if (itmask[i]) {
                double de_local = de[i][j][k];
                double dom_local = dom;
                edot[i] += (
                    // Eccitazioni collisioni
                    -ceHI[i - is] * HI[i][j][k] * de_local
                    -ceHeI[i - is] * HeII[i][j][k] * de_local * de_local * dom_local / 4.0
                    -ceHeII[i - is] * HeII[i][j][k] * de_local / 4.0
                    // Ionizzazioni collisioni
                    -ciHI[i - is] * HI[i][j][k] * de_local
                    -ciHeI[i - is] * HeI[i][j][k] * de_local / 4.0
                    -ciHeII[i - is] * HeII[i][j][k] * de_local / 4.0
                    -ciHeIS[i - is] * HeII[i][j][k] * de_local * de_local * dom_local / 4.0
                    // Ricombinazioni
                    -reHII[i - is] * HII[i][j][k] * de_local
                    -reHeII1[i - is] * HeII[i][j][k] * de_local / 4.0
                    -reHeII2[i - is] * HeII[i][j][k] * de_local / 4.0
                    -reHeIII[i - is] * HeIII[i][j][k] * de_local / 4.0
                    // Bremsstrahlung
                    -brem[i - is] * (HII[i][j][k] + HeII[i][j][k] / 4.0 + HeIII[i][j][k]) * de_local
                );
            }
        }
    }

    // --- Raffreddamento da H2 ---
    if (ispecies > 1) {
        // Uso delle tabelle di Glover & Abel 2008
        for (i = is; i <= ie; ++i) {
            if (itmask[i]) {
                ind = indixe[i - is] - 1; // Indice base 0
                gaHI[i - is] = gaHIa[ind] + tdef[i - is] * (gaHIa[ind + 1] - gaHIa[ind]);
                gaH2[i - is] = gaH2a[ind] + tdef[i - is] * (gaH2a[ind + 1] - gaH2a[ind]);
                gaHe[i - is] = gaHea[ind] + tdef[i - is] * (gaHea[ind + 1] - gaHea[ind]);
                gaHp[i - is] = gaHpa[ind] + tdef[i - is] * (gaHpa[ind + 1] - gaHpa[ind]);
                gael[i - is] = gaela[ind] + tdef[i - is] * (gaela[ind + 1] - gaela[ind]);
                h2lte[i - is] = h2ltea[ind] + tdef[i - is] * (h2ltea[ind + 1] - h2ltea[ind]);
                cieco[i - is] = ciecoa[ind] + tdef[i - is] * (ciecoa[ind + 1] - ciecoa[ind]);
            }
        }

        for (i = is; i <= ie; ++i) {
            if (itmask[i]) {
                double H2I_local = H2I[i][j][k];
                double HI_local = HI[i][j][k];
                double HII_local = HII[i][j][k];
                double HeI_local = HeI[i][j][k];
                double de_local = de[i][j][k];
                double dom_local = dom;

                galdl[i - is] = gaHI[i - is] * HI_local +
                                gaH2[i - is] * H2I_local / 2.0 +
                                gaHe[i - is] * HeI_local / 4.0 +
                                gaHp[i - is] * HII_local +
                                gael[i - is] * de_local;

                double h2lte_local = h2lte[i - is];
                gphdl1 = h2lte_local / dom_local;
                double denominator = 1.0 + gphdl1 / galdl[i - is];

                // Correzione per profondità ottica
                if (ih2optical == 1) {
                    fudge = std::pow(0.76 * d[i][j][k] * dom_local / 8e9, -0.45);
                    fudge = std::min(fudge, 1.0);
                } else {
                    fudge = 1.0;
                }

                edot[i] -= ih2co * fudge * H2I_local * h2lte_local / denominator / (2.0 * dom_local);
            }
        }

        // Raffreddamento CIE
        if (iciecool == 1) {
            for (i = is; i <= ie; ++i) {
                if (itmask[i]) {
                    if (d[i][j][k] * dom > 1e10) {
                        tau = std::pow((d[i][j][k] / 2e16) * dom, 2.8);
                        tau = std::max(tau, 1e-5);
                        ciefudge = (1.0 - std::exp(-tau)) / tau;
                        tau = std::pow((d[i][j][k] / 2e18) * dom, 8.0);
                        tau = std::max(tau, 1e-5);
                        ciefudge *= (1.0 - std::exp(-tau)) / tau;
                        ciefudge = std::min(ciefudge, 1.0);
                        edot[i] = ciefudge * (edot[i] - H2I[i][j][k] * d[i][j][k] * cieco[i - is]);
                    }
                }
            }
        }
    }

    // --- Raffreddamento da HD ---
    if (ispecies > 2) {
        for (i = is; i <= ie; ++i) {
            if (itmask[i]) {
                ind = indixe[i - is] - 1;
                if (tgas[i] > comp2_local) {
                    hdlte[i - is] = hdltea[ind] + tdef[i - is] * (hdltea[ind + 1] - hdltea[ind]);
                    hdlow[i - is] = hdlowa[ind] + tdef[i - is] * (hdlowa[ind + 1] - hdlowa[ind]);
                } else {
                    hdlte[i - is] = tiny;
                    hdlow[i - is] = tiny;
                }
            }
        }

        for (i = is; i <= ie; ++i) {
            if (itmask[i]) {
                double HDI_local = HDI[i][j][k];
                double HI_local = HI[i][j][k];
                double dom_local = dom;
                hdlte1 = hdlte[i - is] / (HI_local * dom_local);
                hdlow1 = std::max(hdlow[i - is], tiny);
                edot[i] -= HDI_local * (hdlte[i - is] / (1.0 + hdlte1 / hdlow1)) / (3.0 * dom_local);
            }
        }
    }

    // --- Calcolo del rapporto polvere-gas ---
    if (anydust || (igammah > 0)) {
        if (idustfield > 0) {
            for (i = is; i <= ie; ++i) {
                if (itmask[i]) {
                    dust2gas[i - is] = dust[i][j][k] / d[i][j][k];
                }
            }
        } else {
            for (i = is; i <= ie; ++i) {
                if (itmask[i]) {
                    dust2gas[i - is] = fgr * metallicity[i - is];
                }
            }
        }
    }

    // --- Calcolo del campo di radiazione interstellare ---
    if (anydust || (igammah > 1)) {
        if (iisrffield > 0) {
            for (i = is; i <= ie; ++i) {
                if (itmask[i]) {
                    myisrf[i - is] = isrf_habing[i][j][k];
                }
            }
        } else {
            for (i = is; i <= ie; ++i) {
                if (itmask[i]) {
                    myisrf[i - is] = isrf;
                }
            }
        }
    }

    // --- Trasferimento di calore gas-grano ---
    if (anydust) {
        // Lookup dei tassi di trasferimento calore gas-grano
        for (i = is; i <= ie; ++i) {
            if (itmask[i]) {
                ind = indixe[i - is] - 1;
                double gasgr_value = gasgra[ind] + tdef[i - is] * (gasgra[ind + 1] - gasgra[ind]);
                gasgr_tdust[i - is] = fgr * gasgr_value * coolunit / mh;
            }
        }

        Calcolo della temperatura della polvere
         Dovrai implementare la funzione calc_tdust_1d_g in base alle tue esigenze
         Ad esempio:
         calc_tdust_1d_g(...);

        // Calcolo del tasso di raffreddamento
        for (i = is; i <= ie; ++i) {
            if (itmask[i]) {
                edot[i] -= (gasgr_tdust[i - is] * (tgas[i] - tdust[i - is]) *
                            dust2gas[i - is] * rhoH[i - is] * rhoH[i - is]);
            }
        }
    }

    // --- Riscaldamento radiativo esterno ---
    if (ispecies > 0) {
        if (iradshield == 0) {
            // Nessun shielding
            for (i = is; i <= ie; ++i) {
                if (itmask[i]) {
                    edot[i] += ipiht * (
                        piHI * HI[i][j][k] +
                        piHeI * HeI[i][j][k] * 0.25 +
                        piHeII * HeII[i][j][k] * 0.25
                    ) / dom;
                }
            }
        } else if (iradshield == 1) {
            // Shielding approssimato
            for (i = is; i <= ie; ++i) {
                if (itmask[i]) {
                    if (k24 < tiny8) {
                        fSShHI = 1.0;
                    } else {
                        nSSh = 6.73e-3 *
                               std::pow(avgsighi / 2.49e-18, -2.0 / 3.0) *
                               std::pow(tgas[i] / 1.0e4, 0.17) *
                               std::pow(k24 / (tbase1 * 1.0e-12), 2.0 / 3.0);
                        nratio = (HI[i][j][k] + HII[i][j][k]) * dom / nSSh;
                        fSShHI = 0.98 * std::pow(1.0 + std::pow(nratio, 1.64), -2.28) +
                                 0.02 * std::pow(1.0 + nratio, -0.84);
                    }

                    edot[i] += ipiht * (
                        piHI * HI[i][j][k] * fSShHI +
                        piHeI * HeI[i][j][k] * 0.25 +
                        piHeII * HeII[i][j][k] * 0.25
                    ) / dom;
                }
            }
        }
        // Altri casi di iradshield possono essere aggiunti qui
    }

    // --- Riscaldamento fotoelettrico ---
    if (igammah > 0) {
        if (igammah == 1) {
            for (i = is; i <= ie; ++i) {
                if (itmask[i]) {
                    if (tgas[i] > 2e4) {
                        gammaha_eff[i - is] = 0.0;
                    } else {
                        gammaha_eff[i - is] = gammaha;
                    }
                }
            }
        } else if (igammah == 2) {
            // Implementazione per igammah == 2
        } else if (igammah == 3) {
            // Implementazione per igammah == 3
        }

        // Calcolo del riscaldamento
        for (i = is; i <= ie; ++i) {
            if (itmask[i]) {
                edot[i] += gammaha_eff[i - is] * rhoH[i - is] *
                           dom_inv * dust2gas[i - is] / fgr;
            }
        }
    }

    // --- Aggiornamento di tgasold ---
    for (i = is; i <= ie; ++i) {
        if (itmask[i]) {
            tgasold[i] = tgas[i];
        }
    }
}

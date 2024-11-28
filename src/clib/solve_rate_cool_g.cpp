// Inclusione delle librerie necessarie
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <omp.h> // Per la parallelizzazione con OpenMP

// Definizione delle costanti fisiche
const double mass_h = /* valore appropriato per la massa dell'idrogeno in cgs */;
const double pi_val = M_PI;
const double kboltz = /* Costante di Boltzmann in cgs */;
const double GravConst = /* Costante gravitazionale in cgs */;

// Definizione dei tipi per compatibilità con il codice Fortran
typedef double R_PREC;

// Definizione delle costanti per la tolleranza
#ifdef GRACKLE_FLOAT_4
const R_PREC tolerance = 1.0e-05;
#endif

#ifdef GRACKLE_FLOAT_8
const R_PREC tolerance = 1.0e-10;
#endif

// Definizione della funzione solve_rate_cool_g
void solve_rate_cool_g(
    // Argomenti generali
    int icool, int in, int jn, int kn, int is, int js, int ks, int ie, int je, int ke, int nratec,
    int iexpand, int ih2co, int ipiht, int ispecies, int imetal, int idim,
    int &ierr, int imcool, int idust, int idustall, int idustfield, int idustrec,
    int igammah, int ih2optical, int iciecool, int ithreebody,
    int ndratec, int clnew, int iVheat, int iMheat, int iTfloor,
    int iH2shield, int iradshield,
    int iradtrans, int iradcoupled, int iradstep, int irt_honly,
    int iisrffield, int iH2shieldcustom, int itmax, int exititmax,
    double dx, double dt, double aye, double temstart, double temend, double gamma_val,
    double utim, double uxyz, double uaye, double urho, double utem, double fh, double dtoh, double z_solar,
    double fgr, double dtemstart, double dtemend, double clEleFra, double Tfloor_scalar,
    // Campi di densità, energia e velocità (array 3D)
    std::vector<std::vector<std::vector<R_PREC>>> &de, std::vector<std::vector<std::vector<R_PREC>>> &HI,
    std::vector<std::vector<std::vector<R_PREC>>> &HII, std::vector<std::vector<std::vector<R_PREC>>> &HeI,
    std::vector<std::vector<std::vector<R_PREC>>> &HeII, std::vector<std::vector<std::vector<R_PREC>>> &HeIII,
    std::vector<std::vector<std::vector<R_PREC>>> &HM, std::vector<std::vector<std::vector<R_PREC>>> &H2I,
    std::vector<std::vector<std::vector<R_PREC>>> &H2II, std::vector<std::vector<std::vector<R_PREC>>> &DI,
    std::vector<std::vector<std::vector<R_PREC>>> &DII, std::vector<std::vector<std::vector<R_PREC>>> &HDI,
    std::vector<std::vector<std::vector<R_PREC>>> &d, std::vector<std::vector<std::vector<R_PREC>>> &e,
    std::vector<std::vector<std::vector<R_PREC>>> &u, std::vector<std::vector<std::vector<R_PREC>>> &v,
    std::vector<std::vector<std::vector<R_PREC>>> &w, std::vector<std::vector<std::vector<R_PREC>>> &metal,
    std::vector<std::vector<std::vector<R_PREC>>> &dust,
    std::vector<std::vector<std::vector<R_PREC>>> &Vheat, std::vector<std::vector<std::vector<R_PREC>>> &Mheat,
    std::vector<std::vector<std::vector<R_PREC>>> &Tfloor,
    // Campi per il trasferimento radiativo
    std::vector<std::vector<std::vector<R_PREC>>> &kphHI, std::vector<std::vector<std::vector<R_PREC>>> &kphHeI,
    std::vector<std::vector<std::vector<R_PREC>>> &kphHeII, std::vector<std::vector<std::vector<R_PREC>>> &kdissH2I,
    std::vector<std::vector<std::vector<R_PREC>>> &photogamma,
    // Campo di lunghezza di auto-schermatura H2
    std::vector<std::vector<std::vector<R_PREC>>> &xH2shield,
    // Campo di radiazione interstellare per il riscaldamento della polvere
    std::vector<std::vector<std::vector<R_PREC>>> &isrf_habing,
    // Fattore di schermatura H2 personalizzato
    std::vector<std::vector<std::vector<R_PREC>>> &f_shield_custom,
    // Tabelle di raffreddamento (tassi di raffreddamento in funzione della temperatura)
    std::vector<double> &hyd01ka, std::vector<double> &h2k01a, std::vector<double> &vibha,
    std::vector<double> &rotha, std::vector<double> &rotla, std::vector<double> &gpldla,
    std::vector<double> &gphdla, std::vector<double> &hdltea, std::vector<double> &hdlowa,
    std::vector<double> &gaHIa, std::vector<double> &gaH2a, std::vector<double> &gaHea,
    std::vector<double> &gaHpa, std::vector<double> &gaela, std::vector<double> &h2ltea,
    std::vector<double> &gasgra, std::vector<double> &ciecoa,
    std::vector<double> &ceHIa, std::vector<double> &ceHeIa, std::vector<double> &ceHeIIa,
    std::vector<double> &ciHIa, std::vector<double> &ciHeIa, std::vector<double> &ciHeISa,
    std::vector<double> &ciHeIIa, std::vector<double> &reHIIa, std::vector<double> &reHeII1a,
    std::vector<double> &reHeII2a, std::vector<double> &reHeIIIa, std::vector<double> &brema,
    double compa, double piHI, double piHeI, double piHeII, double comp_xraya, double comp_temp,
    double gammaha, double isrf, std::vector<double> &regra, double gamma_isrfa,
    double avgsighi, double avgsighei, double avgsigheii,
    // Tabelle di chimica (tassi in funzione della temperatura)
    std::vector<double> &k1a, std::vector<double> &k2a, std::vector<double> &k3a, std::vector<double> &k4a,
    std::vector<double> &k5a, std::vector<double> &k6a, std::vector<double> &k7a, std::vector<double> &k8a,
    std::vector<double> &k9a, std::vector<double> &k10a, std::vector<double> &k11a, std::vector<double> &k12a,
    std::vector<double> &k13a, std::vector<double> &k14a, std::vector<double> &k15a, std::vector<double> &k16a,
    std::vector<double> &k17a, std::vector<double> &k18a, std::vector<double> &k19a, std::vector<double> &k22a,
    std::vector<double> &k50a, std::vector<double> &k51a, std::vector<double> &k52a, std::vector<double> &k53a,
    std::vector<double> &k54a, std::vector<double> &k55a, std::vector<double> &k56a,
    std::vector<double> &k57a, std::vector<double> &k58a,
    std::vector<std::vector<double>> &k13dda, std::vector<std::vector<double>> &h2dusta,
    std::vector<double> &ncrna, std::vector<double> &ncrd1a, std::vector<double> &ncrd2a,
    double k24, double k25, double k26, double k27, double k28, double k29, double k30, double k31,
    // Dati di raffreddamento Cloudy
    int icmbTfloor, int iClHeat,
    int64_t priGridRank, int64_t priDataSize,
    int64_t metGridRank, int64_t metDataSize,
    std::vector<int64_t> &priGridDim, std::vector<int64_t> &metGridDim,
    std::vector<double> &priPar1, std::vector<double> &priPar2, std::vector<double> &priPar3,
    std::vector<double> &priPar4, std::vector<double> &priPar5,
    std::vector<double> &metPar1, std::vector<double> &metPar2, std::vector<double> &metPar3,
    std::vector<double> &metPar4, std::vector<double> &metPar5,
    std::vector<double> &priCooling, std::vector<double> &priHeating, std::vector<double> &priMMW,
    std::vector<double> &metCooling, std::vector<double> &metHeating
) {
    // Impostazione dell'indicatore di errore
    ierr = 1;

    // Impostazione del flag per le opzioni legate alla polvere
    bool anydust = (idust > 0) || (idustall > 0);

    // Impostazione delle unità
    double dom = urho * pow(aye, 3) / mass_h;
    double tbase1 = utim;
    double xbase1 = uxyz / (aye * uaye); // uxyz è [x]*a = [x]*[a]*a'
    double dbase1 = urho * pow(aye * uaye, 3); // urho è [dens]/a^3 = [dens]/([a]*a')^3
    double coolunit = (pow(uaye, 5) * pow(xbase1, 2) * pow(mass_h, 2)) / (pow(tbase1, 3) * dbase1);
    double uvel = (uxyz / aye) / utim;
    double chunit = (1.60218e-12) / (2.0 * uvel * uvel * mass_h); // 1 eV per H2 formato

    double dx_cgs = dx * xbase1;
    double c_ljeans = sqrt((gamma_val * pi_val * kboltz) / (GravConst * mass_h * dbase1));

    double dlogtem = (log(temend) - log(temstart)) / static_cast<double>(nratec - 1);

    // Converti le densità da comoventi a proprie
    if (iexpand == 1) {
        scale_fields_g(
            d, de, HI, HII, HeI, HeII, HeIII,
            HM, H2I, H2II, DI, DII, HDI, metal, dust,
            is, ie, js, je, ks, ke,
            in, jn, kn, ispecies, imetal, idustfield,
            pow(aye, -3)
        );
    }

    // Applica il ceiling alle specie
    ceiling_species_g(
        d, de, HI, HII, HeI, HeII, HeIII,
        HM, H2I, H2II, DI, DII, HDI, metal,
        is, ie, js, je, ks, ke,
        in, jn, kn, ispecies, imetal
    );

    // Loop sulle zone, elaborando un'intera colonna i in una volta
    int dk = ke - ks + 1;
    int dj = je - js + 1;

    // Parallelizza i loop su k e j con OpenMP
    // Flat j e k loops per migliorare il parallelismo
    #pragma omp parallel for schedule(runtime) \
    private(\
        i, j, k, iter, \
        ttmin, energy, comp1, comp2, \
        heq1, heq2, eqk221, eqk222, eqk131, eqk132, \
        eqt1, eqt2, eqtdef, dheq, heq, \
        indixe, \
        t1, t2, logtem, tdef, \
        dtit, ttot, p2d, tgas, tgasold, \
        tdust, metallicity, dust2gas, rhoH, mmw, \
        mynh, myde, gammaha_eff, gasgr_tdust, regr, olddtit, \
        HIp, HIIp, HeIp, HeIIp, HeIIIp, \
        HMp, H2Ip, H2IIp, \
        dep, dedot, HIdot, dedot_prev, \
        DIp, DIIp, HDIp, HIdot_prev, \
        k24shield, k25shield, k26shield, \
        k28shield, k29shield, k30shield, \
        k31shield, \
        k1, k2, k3, k4, k5, \
        k6, k7, k8, k9, k10, \
        k11, k12, k13, k14, k15, \
        k16, k17, k18, k19, k22, \
        k50, k51, k52, k53, k54, \
        k55, k56, k57, k58, k13dd, h2dust, \
        ncrn, ncrd1, ncrd2, \
        ceHI, ceHeI, ceHeII, \
        ciHI, ciHeI, ciHeIS, ciHeII, \
        reHII, reHeII1, reHeII2, reHeIII, \
        brem, edot, \
        hyd01k, h2k01, vibh, roth, rotl, \
        gpldl, gphdl, hdlte, hdlow, cieco, \
        itmask \
    ) shared(ierr)
    for (int t = 0; t < dk * dj; ++t) {
        k = t / dj + ks + 1;
        j = t % dj + js + 1;

        // Inizializza la maschera di iterazione a true per tutte le celle
        std::vector<bool> itmask(ie - is + 2, true);

        // Se stiamo usando la radiazione accoppiata con step intermedi
        if (iradtrans == 1) {
            if (iradcoupled == 1 && iradstep == 1) {
                for (int i = is + 1; i <= ie + 1; ++i) {
                    if (kphHI[i][j][k] > 0) {
                        itmask[i] = true;
                    } else {
                        itmask[i] = false;
                    }
                }
            }
            // Solver di rate normale, ma non conteggiare due volte le celle con radiazione
            if (iradcoupled == 1 && iradstep == 0) {
                for (int i = is + 1; i <= ie + 1; ++i) {
                    if (kphHI[i][j][k] > 0) {
                        itmask[i] = false;
                    } else {
                        itmask[i] = true;
                    }
                }
            }
        }

        // Inizializza il tempo trascorso a zero per ogni cella
        std::vector<double> ttot(ie - is + 2, 0.0);

        // ------------------ Loop sulle subiterazioni ----------------
        for (iter = 1; iter <= itmax; ++iter) {
            // Inizializza dtit
            std::vector<double> dtit(ie - is + 2, std::numeric_limits<double>::max());
            
            // Calcola il rate di raffreddamento, tgas, tdust e metallicity per questa riga
            // (Implementa la funzione cool1d_multi_g in base al codice originale)
            //            call cool1d_multi_g(
            //     &                d, e, u, v, w, de, HI, HII, HeI, HeII, HeIII,
            //     &                in, jn, kn, nratec,
            //     &                iexpand, ispecies, imetal, imcool,
            //     &                idust, idustall, idustfield, idustrec,
            //     &                idim, is, ie, j, k, ih2co, ipiht, iter, igammah,
            //     &                aye, temstart, temend, z_solar, fgr,
            //     &                utem, uxyz, uaye, urho, utim,
            //     &                gamma, fh,
            //     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa,
            //     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a,
            //     &                reHeII2a, reHeIIIa, brema, compa, gammaha,
            //     &                isrf, regra, gamma_isrfa, comp_xraya, comp_temp,
            //     &                piHI, piHeI, piHeII, comp1, comp2,
            //     &                HM, H2I, H2II, DI, DII, HDI, metal, dust,
            //     &                hyd01ka, h2k01a, vibha, rotha, rotla,
            //     &                hyd01k, h2k01, vibh, roth, rotl,
            //     &                gpldla, gphdla, gpldl, gphdl,
            //     &                hdltea, hdlowa, hdlte, hdlow,
            //     &                gaHIa, gaH2a, gaHea, gaHpa, gaela,
            //     &                h2ltea, gasgra,
            //     &                ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeIS, ciHeII,
            //     &                reHII, reHeII1, reHeII2, reHeIII, brem,
            //     &                indixe, t1, t2, logtem, tdef, edot,
            //     &                tgas, tgasold, mmw, p2d, tdust, metallicity,
            //     &                dust2gas, rhoH, mynh, myde,
            //     &                gammaha_eff, gasgr_tdust, regr,
            //     &                iradshield, avgsighi, avgsighei, avgsigheii,
            //     &                k24, k26,
            //     &                iradtrans, photogamma,
            //     &                ih2optical, iciecool, ciecoa, cieco,
            //     &                icmbTfloor, iClHeat, clEleFra,
            //     &                priGridRank, priGridDim,
            //     &                priPar1, priPar2, priPar3, priPar4, priPar5,
            //     &                priDataSize, priCooling, priHeating, priMMW,
            //     &                metGridRank, metGridDim,
            //     &                metPar1, metPar2, metPar3, metPar4, metPar5,
            //     &                metDataSize, metCooling, metHeating, clnew,
            //     &                iVheat, iMheat, Vheat, Mheat,
            //     &                iTfloor, Tfloor_scalar, Tfloor,
            //     &                iisrffield, isrf_habing, itmask)
            
            if (ispecies > 0) {
                // Look-up dei rate in funzione della temperatura per un set 1D di zone
                // (Implementa la funzione lookup_cool_rates1d_g in base al codice originale)
                //                call lookup_cool_rates1d_g(temstart, temend, nratec, j, k,
                //         &               is, ie, ithreebody,
                //         &               in, jn, kn, ispecies, anydust,
                //         &               iH2shield, iradshield,
                //         &               tgas, mmw, d, HI, HII, HeI, HeII, HeIII,
                //         &               HM, H2I, H2II, DI, DII, HDI,
                //         &               tdust, dust2gas,
                //         &               k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
                //         &               k11a, k12a, k13a, k13dda, k14a, k15a, k16a,
                //         &               k17a, k18a, k19a, k22a,
                //         &               k50a, k51a, k52a, k53a, k54a, k55a, k56a,
                //         &               k57a, k58a, ndratec, dtemstart, dtemend, h2dusta,
                //         &               ncrna, ncrd1a, ncrd2a,
                //         &               avgsighi, avgsighei, avgsigheii, piHI, piHeI,
                //         &               k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
                //         &               k11, k12, k13, k14, k15, k16, k17, k18,
                //         &               k19, k22, k24, k25, k26, k28, k29, k30, k31,
                //         &               k50, k51, k52, k53, k54, k55, k56, k57,
                //         &               k58, k13dd, k24shield, k25shield, k26shield,
                //         &               k28shield, k29shield, k30shield,
                //         &               k31shield, h2dust, ncrn, ncrd1, ncrd2,
                //         &               t1, t2, tdef, logtem, indixe,
                //         &               dom, coolunit, tbase1, xbase1, dx_cgs, c_ljeans,
                //         &               iradtrans, kdissH2I, xH2shield, iH2shieldcustom,
                //         &               f_shield_custom,  itmask)
                // Calcola dedot e HIdot, i rate di cambiamento di de e HI
                // (Implementa la funzione rate_timestep_g in base al codice originale)
                //                call rate_timestep_g(
                //         &                     dedot, HIdot, ispecies, anydust,
                //         &                     de, HI, HII, HeI, HeII, HeIII, d,
                //         &                     HM, H2I, H2II,
                //         &                     in, jn, kn, is, ie, j, k,
                //         &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
                //         &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
                //         &                     k24, k25, k26, k27, k28, k29, k30,
                //         &                     k50, k51, k52, k53, k54, k55, k56, k57, k58,
                //         &                     h2dust, ncrn, ncrd1, ncrd2, rhoH,
                //         &                     k24shield, k25shield, k26shield,
                //         &                     k28shield, k29shield, k30shield, k31shield,
                //         &                     iradtrans, irt_honly,
                //         &                     kphHI, kphHeI, kphHeII,
                //         &                     itmask, edot, chunit, dom)
                
                // Trova il timestep che mantiene i cambiamenti chimici relativi sotto il 10%
                for (int i = is + 1; i <= ie + 1; ++i) {
                    if (itmask[i]) {
                        if (std::abs(dedot[i]) < std::numeric_limits<double>::epsilon())
                            dedot[i] = std::min(tiny, de[i][j][k]);
                        if (std::abs(HIdot[i]) < std::numeric_limits<double>::epsilon())
                            HIdot[i] = std::min(tiny, HI[i][j][k]);
                        
                        if (std::min(std::abs(k1[i] * de[i][j][k] * HI[i][j][k]),
                                     std::abs(k2[i] * HII[i][j][k] * de[i][j][k])) /
                            std::max(std::abs(dedot[i]), std::abs(HIdot[i])) >
                            1.0e6) {
                            dedot[i] = tiny;
                            HIdot[i] = tiny;
                        }
                        
                        if (iter > 50) {
                            dedot[i] = std::min(std::abs(dedot[i]), std::abs(dedot_prev[i]));
                            HIdot[i] = std::min(std::abs(HIdot[i]), std::abs(HIdot_prev[i]));
                        }
                        
                        double olddtit = dtit[i];
                        dtit[i] = std::min({std::abs(0.1 * de[i][j][k] / dedot[i]),
                            std::abs(0.1 * HI[i][j][k] / HIdot[i]),
                            dt - ttot[i], 0.5 * dt});
                        
                        // Ulteriori calcoli per specifiche condizioni (da implementare)
                        
                        
                        // Check the condition as per the Fortran code
                        if (d[i][j][k] * dom > 1e8 && edot[i] > 0.0 && ispecies > 1) {
                            // Equilibrium value for H is calculated as:
                            // H = (-1.0 / (4.0 * k22)) * (k13 - sqrt(8 * k13 * k22 * rho + k13^2))
                            
                            // Compute eqt2
                            double eqt2 = std::min(std::log(tgas[i]) + 0.1 * dlogtem, t2[i]);
                            
                            // Compute eqtdef
                            double eqtdef = (eqt2 - t1[i]) / (t2[i] - t1[i]);
                            
                            // Compute eqk222
                            double eqk222 = k22a[indixe[i]] +
                            (k22a[indixe[i] + 1] - k22a[indixe[i]]) * eqtdef;
                            
                            // Compute eqk132
                            double eqk132 = k13a[indixe[i]] +
                            (k13a[indixe[i] + 1] - k13a[indixe[i]]) * eqtdef;
                            
                            // Compute heq2
                            double heq2 = (-1.0 / (4.0 * eqk222)) *
                            (eqk132 - sqrt(8.0 * eqk132 * eqk222 * fh * d[i][j][k] +
                                           eqk132 * eqk132));
                            
                            // Compute eqt1
                            double eqt1 = std::max(std::log(tgas[i]) - 0.1 * dlogtem, t1[i]);
                            
                            // Compute eqtdef again
                            eqtdef = (eqt1 - t1[i]) / (t2[i] - t1[i]);
                            
                            // Compute eqk221
                            double eqk221 = k22a[indixe[i]] +
                            (k22a[indixe[i] + 1] - k22a[indixe[i]]) * eqtdef;
                            
                            // Compute eqk131
                            double eqk131 = k13a[indixe[i]] +
                            (k13a[indixe[i] + 1] - k13a[indixe[i]]) * eqtdef;
                            
                            // Compute heq1
                            double heq1 = (-1.0 / (4.0 * eqk221)) *
                            (eqk131 - sqrt(8.0 * eqk131 * eqk221 * fh * d[i][j][k] +
                                           eqk131 * eqk131));
                            
                            // Compute dheq
                            double dheq = (std::abs(heq2 - heq1) / (std::exp(eqt2) - std::exp(eqt1))) *
                            (tgas[i] / p2d[i]) * edot[i];
                            
                            // Compute heq
                            double heq = (-1.0 / (4.0 * k22[i])) *
                            (k13[i] - sqrt(8.0 * k13[i] * k22[i] * fh * d[i][j][k] +
                                           k13[i] * k13[i]));
                            
                            // Optional debugging output
                            if (d[i][j][k] * dom > 1e18 && i == 4) {
#pragma omp critical
                                {
                                    std::cout << HI[i][j][k] / heq << " " << edot[i] << " " << tgas[i] << std::endl;
                                }
                            }
                            
                            // Update dtit[i]
                            dtit[i] = std::min(dtit[i], 0.1 * heq / dheq);
                        }
                        
                        // Adjust dtit[i] if iter > 10
                        if (iter > 10) {
                            dtit[i] = std::min(olddtit * 1.5, dtit[i]);
                        }
                        
                        // Debugging output (if enabled)
#ifndef DONT_WRITE_COOLING_DEBUG
#ifdef WRITE_COOLING_DEBUG
                        if ((dtit[i] / dt < 1.0e-2) && (iter > 800) && (std::abs((dt - ttot[i]) / dt) > 1.0e-3)) {
#pragma omp critical
                            {
                                // Output debugging information
                                std::ofstream debug_file("cooling_debug.log", std::ios::app);
                                if (debug_file.is_open()) {
                                    debug_file << std::scientific;
                                    debug_file.precision(11);
                                    // Write iteration info
                                    debug_file << iter << " " << i << " " << j << " " << k << " "
                                    << dtit[i] << " " << ttot[i] << " " << dt << " "
                                    << de[i][j][k] << " " << dedot[i] << " "
                                    << HI[i][j][k] << " " << HIdot[i] << " "
                                    << tgas[i] << " " << dedot_prev[i] << " " << HIdot_prev[i] << "\n";
                                    // Write species densities
                                    debug_file << HI[i][j][k] << " " << HII[i][j][k] << " "
                                    << HeI[i][j][k] << " " << HeII[i][j][k] << " " << HeIII[i][j][k] << " "
                                    << HM[i][j][k] << " " << H2I[i][j][k] << " "
                                    << H2II[i][j][k] << " " << de[i][j][k] << "\n";
                                    // Write rate calculations
                                    debug_file
                                    << -k1[i] * de[i][j][k] * HI[i][j][k] << " "
                                    << -k7[i] * de[i][j][k] * HI[i][j][k] << " "
                                    << -k8[i] * HM[i][j][k] * HI[i][j][k] << " "
                                    << -k9[i] * HII[i][j][k] * HI[i][j][k] << " "
                                    << -k10[i] * H2II[i][j][k] * HI[i][j][k] / 2.0 << " "
                                    << -2.0 * k22[i] * HI[i][j][k] * HI[i][j][k] * HI[i][j][k] << " "
                                    << k2[i] * HII[i][j][k] * de[i][j][k] << " "
                                    << 2.0 * k13[i] * HI[i][j][k] * H2I[i][j][k] / 2.0 << " "
                                    << k11[i] * HII[i][j][k] * H2I[i][j][k] / 2.0 << " "
                                    << 2.0 * k12[i] * de[i][j][k] * H2I[i][j][k] / 2.0 << " "
                                    << k14[i] * HM[i][j][k] * de[i][j][k] << " "
                                    << k15[i] * HM[i][j][k] * HI[i][j][k] << " "
                                    << 2.0 * k16[i] * HM[i][j][k] * HII[i][j][k] << " "
                                    << 2.0 * k18[i] * H2II[i][j][k] * de[i][j][k] / 2.0 << " "
                                    << k19[i] * H2II[i][j][k] * HM[i][j][k] / 2.0 << " "
                                    << -k57[i] * HI[i][j][k] * HI[i][j][k] << " "
                                    << -k58[i] * HI[i][j][k] * HeI[i][j][k] / 4.0 << "\n";
                                    debug_file.close();
                                }
                            }
                        }
#endif
#endif // DONT_WRITE_COOLING_DEBUG
                        
                        
                    } else {
                        dtit[i] = dt;
                    }
                }
            }
            
            // Calcola il massimo timestep per il raffreddamento/riscaldamento
            for (int i = is + 1; i <= ie + 1; ++i) {
                if (itmask[i]) {
                    double energy = std::max(p2d[i] / (gamma_val - 1.0), tiny);
                    if (tgas[i] <= 1.01 * temstart && edot[i] < 0.0)
                        edot[i] = tiny;
                    if (std::abs(edot[i]) < tiny)
                        edot[i] = tiny;
                    
                    dtit[i] = std::min({std::abs(0.1 * energy / edot[i]),
                        dt - ttot[i], dtit[i]});
                    
                    // Debugging e controllo degli errori (da implementare)
                    if (std::isnan(dtit[i])) {
#pragma omp critical
                        {
                            std::cout << "HUGE dtit :: " << energy << " " << edot[i] << " " << dtit[i] << " "
                            << dt << " " << ttot[i] << " " << std::abs(0.1 * energy / edot[i]) << " "
                            << std::abs(0.1 * energy / edot[i]) << std::endl;
                        }
                    }
#ifdef WRITE_COOLING_DEBUG
                    // If the timestep is too small, then output some debugging info
                    if (((dtit[i] / dt < 1.0e-2) && (iter > 800) || (iter > itmax - 100)) &&
                        (std::abs((dt - ttot[i]) / dt) > 1.0e-3)) {
                        
#pragma omp critical
                        {
                            // Open a file for writing debug information or use standard output
                            // For this example, we'll write to a file named "cooling_debug.log"
                            std::ofstream debug_file("cooling_debug.log", std::ios::app);
                            if (debug_file.is_open()) {
                                // Set formatting to match the Fortran format
                                debug_file << std::scientific << std::setprecision(3);
                                
                                // Output the data with formatting similar to the Fortran format
                                debug_file << std::setw(4) << i << " "
                                << std::setw(4) << j << " "
                                << std::setw(4) << k << " "
                                << std::setw(4) << iter << " ";
                                debug_file << std::setw(14) << e[i][j][k] << " "
                                << std::setw(14) << edot[i] << " "
                                << std::setw(14) << tgas[i] << " "
                                << std::setw(14) << energy << " "
                                << std::setw(14) << de[i][j][k] << " "
                                << std::setw(14) << ttot[i] << " "
                                << std::setw(14) << d[i][j][k] << " "
                                << std::setw(14) << e[i][j][k] << " "
                                << std::setw(14) << dtit[i] << "\n";
                                
                                debug_file.close();
                            } else {
                                // If the file cannot be opened, output to standard error
                                std::cerr << "Error opening cooling_debug.log for writing.\n";
                            }
                        } // End of critical section
                    }
#endif /* WRITE_COOLING_DEBUG */
                    
                }
            }
            
            // Aggiorna l'energia totale e del gas
            if (icool == 1) {
                for (int i = is + 1; i <= ie + 1; ++i) {
                    if (itmask[i]) {
                        e[i][j][k] += static_cast<R_PREC>(edot[i] / d[i][j][k] * dtit[i]);
#ifdef WRITE_COOLING_DEBUG
                        // Check if e[i][j][k] is NaN (Not a Number)
                        if (std::isnan(e[i][j][k])) {
#pragma omp critical
                            {
                                // Output edot[i], d[i][j][k], dtit[i]
                                // For example, write to a debug file or standard error
                                std::ofstream debug_file("cooling_debug.log", std::ios::app);
                                if (debug_file.is_open()) {
                                    debug_file << std::scientific << std::setprecision(5);
                                    debug_file << edot[i] << " " << d[i][j][k] << " " << dtit[i] << "\n";
                                    debug_file.close();
                                } else {
                                    // If the file cannot be opened, output to standard error
                                    std::cerr << "edot: " << edot[i] << ", d: " << d[i][j][k]
                                    << ", dtit: " << dtit[i] << std::endl;
                                }
                            }
                        }
#endif /* WRITE_COOLING_DEBUG */
                    }
                }
            }
            
            if (ispecies > 0) {
                // Risolvi le equazioni dei rate con un sweep di Gauss-Seidel
                // (Implementa la funzione step_rate_g in base al codice originale)
                
                //                call step_rate_g(de, HI, HII, HeI, HeII, HeIII, d,
                //         &                     HM, H2I, H2II, DI, DII, HDI, dtit,
                //         &                     in, jn, kn, is, ie, j, k,
                //         &                     ispecies, anydust,
                //         &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
                //         &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
                //         &                     k24, k25, k26, k27, k28, k29, k30,
                //         &                     k50, k51, k52, k53, k54, k55, k56, k57, k58,
                //         &                     h2dust, rhoH,
                //         &                     k24shield, k25shield, k26shield,
                //         &                     k28shield, k29shield, k30shield, k31shield,
                //         &                     HIp, HIIp, HeIp, HeIIp, HeIIIp, dep,
                //         &                     HMp, H2Ip, H2IIp, DIp, DIIp, HDIp,
                //         &                     dedot_prev, HIdot_prev,
                //         &                     iradtrans, irt_honly,
                //         &                     kphHI, kphHeI, kphHeII,
                //         &                     itmask)
                
            }
            
            // Aggiungi il timestep al tempo trascorso per ogni cella e trova il minimo timestep
            ttmin = std::numeric_limits<double>::max();
            for (int i = is + 1; i <= ie + 1; ++i) {
                ttot[i] = std::min(ttot[i] + dtit[i], dt);
                if (std::abs(dt - ttot[i]) < tolerance * dt)
                    itmask[i] = false;
                if (ttot[i] < ttmin)
                    ttmin = ttot[i];
            }
            
            // Se tutte le celle sono completate, esci dal loop
            if (std::abs(dt - ttmin) < tolerance * dt)
                break;
        } // Fine del loop sulle subiterazioni
            // Controllo se il conteggio delle iterazioni supera il massimo
            if (iter > itmax) {
                #pragma omp critical
                {
                    std::cerr << "MULTI_COOL iter > " << itmax << " at j,k = " << j << "," << k << std::endl;
                    std::cerr << "FATAL error (2) in MULTI_COOL" << std::endl;
                    std::cerr << "dt = " << dt << " ttmin = " << ttmin << std::endl;
                    if (exititmax == 1) {
                        ierr = 0;
                    }
                }
                break;
            }

            if (iter > itmax / 2) {
                #pragma omp critical
                {
                    std::cout << "MULTI_COOL iter,j,k = " << iter << "," << j << "," << k << std::endl;
                }
            }


    } // Fine del loop parallelo su t

    // Se è stato prodotto un errore, ritorna ora
    if (ierr == 0) {
        return;
    }

    // Converti le densità da proprie a comoventi
    if (iexpand == 1) {
        scale_fields_g(
            d, de, HI, HII, HeI, HeII, HeIII,
            HM, H2I, H2II, DI, DII, HDI, metal, dust,
            is, ie, js, je, ks, ke,
            in, jn, kn, ispecies, imetal, idustfield,
            pow(aye, 3)
        );
    }

    if (ispecies > 0) {
        // Correggi le specie per garantire la consistenza
        make_consistent_g(
            de, HI, HII, HeI, HeII, HeIII,
            HM, H2I, H2II, DI, DII, HDI, metal,
            d, is, ie, js, je, ks, ke,
            in, jn, kn, ispecies, imetal, fh, dtoh
        );
    }

    return;
}

#include <vector>

// Define R_PREC as double for consistency with Fortran's real*8
typedef double R_PREC;

void scale_fields_g(
    // Density fields (3D arrays)
    std::vector<std::vector<std::vector<R_PREC>>> &d,
    std::vector<std::vector<std::vector<R_PREC>>> &de,
    std::vector<std::vector<std::vector<R_PREC>>> &HI,
    std::vector<std::vector<std::vector<R_PREC>>> &HII,
    std::vector<std::vector<std::vector<R_PREC>>> &HeI,
    std::vector<std::vector<std::vector<R_PREC>>> &HeII,
    std::vector<std::vector<std::vector<R_PREC>>> &HeIII,
    std::vector<std::vector<std::vector<R_PREC>>> &HM,
    std::vector<std::vector<std::vector<R_PREC>>> &H2I,
    std::vector<std::vector<std::vector<R_PREC>>> &H2II,
    std::vector<std::vector<std::vector<R_PREC>>> &DI,
    std::vector<std::vector<std::vector<R_PREC>>> &DII,
    std::vector<std::vector<std::vector<R_PREC>>> &HDI,
    std::vector<std::vector<std::vector<R_PREC>>> &metal,
    std::vector<std::vector<std::vector<R_PREC>>> &dust,
    // Indices
    int is, int ie, int js, int je, int ks, int ke,
    int in, int jn, int kn,
    // Flags
    int ispecies, int imetal, int idustfield,
    // Scaling factor
    double factor)
{
    // Multiply all fields by factor (1/a^3 or a^3)
    for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
            for (int i = is; i <= ie; ++i) {
                d[i][j][k] *= factor;
            }
        }
    }

    if (ispecies > 0) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    de[i][j][k]    *= factor;
                    HI[i][j][k]    *= factor;
                    HII[i][j][k]   *= factor;
                    HeI[i][j][k]   *= factor;
                    HeII[i][j][k]  *= factor;
                    HeIII[i][j][k] *= factor;
                }
            }
        }
    }
    if (ispecies > 1) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    HM[i][j][k]   *= factor;
                    H2I[i][j][k]  *= factor;
                    H2II[i][j][k] *= factor;
                }
            }
        }
    }
    if (ispecies > 2) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    DI[i][j][k]  *= factor;
                    DII[i][j][k] *= factor;
                    HDI[i][j][k] *= factor;
                }
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
    if (idustfield == 1) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    dust[i][j][k] *= factor;
                }
            }
        }
    }
}

#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>

// Define R_PREC as double for consistency with Fortran's real*8
typedef double R_PREC;

// Define the tiny threshold value (adjust as appropriate)
const R_PREC tiny = 1e-20; // Or use the appropriate value from your program
const R_PREC RKIND = 1.0;  // For scaling tiny in HeIII

void ceiling_species_g(
    // Density fields (3D arrays)
    std::vector<std::vector<std::vector<R_PREC>>> &d,
    std::vector<std::vector<std::vector<R_PREC>>> &de,
    std::vector<std::vector<std::vector<R_PREC>>> &HI,
    std::vector<std::vector<std::vector<R_PREC>>> &HII,
    std::vector<std::vector<std::vector<R_PREC>>> &HeI,
    std::vector<std::vector<std::vector<R_PREC>>> &HeII,
    std::vector<std::vector<std::vector<R_PREC>>> &HeIII,
    std::vector<std::vector<std::vector<R_PREC>>> &HM,
    std::vector<std::vector<std::vector<R_PREC>>> &H2I,
    std::vector<std::vector<std::vector<R_PREC>>> &H2II,
    std::vector<std::vector<std::vector<R_PREC>>> &DI,
    std::vector<std::vector<std::vector<R_PREC>>> &DII,
    std::vector<std::vector<std::vector<R_PREC>>> &HDI,
    std::vector<std::vector<std::vector<R_PREC>>> &metal,
    // Indices
    int is, int ie, int js, int je, int ks, int ke,
    int in, int jn, int kn,
    // Flags
    int ispecies, int imetal)
{
    // Loop over the specified indices and apply the ceiling
    if (ispecies > 0) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    de[i][j][k]    = std::max(de[i][j][k], tiny);
                    HI[i][j][k]    = std::max(HI[i][j][k], tiny);
                    HII[i][j][k]   = std::max(HII[i][j][k], tiny);
                    HeI[i][j][k]   = std::max(HeI[i][j][k], tiny);
                    HeII[i][j][k]  = std::max(HeII[i][j][k], tiny);
                    HeIII[i][j][k] = std::max(HeIII[i][j][k], 1e-5 * RKIND * tiny);
                }
            }
        }
    }
    if (ispecies > 1) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    HM[i][j][k]   = std::max(HM[i][j][k], tiny);
                    H2I[i][j][k]  = std::max(H2I[i][j][k], tiny);
                    H2II[i][j][k] = std::max(H2II[i][j][k], tiny);
                }
            }
        }
    }
    if (ispecies > 2) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    DI[i][j][k]  = std::max(DI[i][j][k], tiny);
                    DII[i][j][k] = std::max(DII[i][j][k], tiny);
                    HDI[i][j][k] = std::max(HDI[i][j][k], tiny);
                }
            }
        }
    }
    if (imetal == 1) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                for (int i = is; i <= ie; ++i) {
                    metal[i][j][k] = std::max(metal[i][j][k], tiny);
                    if (metal[i][j][k] > d[i][j][k]) {
                        std::cerr << "WARNING: metal density exceeds total density!\n";
                        std::cerr << "i, j, k, metal, density = "
                                  << i << ", " << j << ", " << k << ", "
                                  << metal[i][j][k] << ", " << d[i][j][k] << "\n";
                    }
                }
            }
        }
    }
}

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

// Define types for consistency with Fortran's real*8 and integer*8
typedef double real8;
typedef int64_t int8;

// Constants (adjust these based on your application's requirements)
const real8 tiny8 = 1e-20;
const real8 DKIND = 1.0; // For casting purposes

void lookup_cool_rates1d_g(
    // Inputs
    double temstart, double temend, int nratec, int j, int k,
    int is, int ie, int ithreebody, int in, int jn, int kn,
    int ispecies, bool anydust, int iH2shield, int iradshield,
    std::vector<double> &tgas1d, std::vector<double> &mmw,
    std::vector<std::vector<std::vector<real8>>> &d,
    std::vector<std::vector<std::vector<real8>>> &HI,
    std::vector<std::vector<std::vector<real8>>> &HII,
    std::vector<std::vector<std::vector<real8>>> &HeI,
    std::vector<std::vector<std::vector<real8>>> &HeII,
    std::vector<std::vector<std::vector<real8>>> &HeIII,
    std::vector<std::vector<std::vector<real8>>> &HM,
    std::vector<std::vector<std::vector<real8>>> &H2I,
    std::vector<std::vector<std::vector<real8>>> &H2II,
    std::vector<std::vector<std::vector<real8>>> &DI,
    std::vector<std::vector<std::vector<real8>>> &DII,
    std::vector<std::vector<std::vector<real8>>> &HDI,
    std::vector<double> &tdust, std::vector<double> &dust2gas,
    // Chemistry rates as a function of temperature
    std::vector<double> &k1a, std::vector<double> &k2a, std::vector<double> &k3a,
    std::vector<double> &k4a, std::vector<double> &k5a, std::vector<double> &k6a,
    std::vector<double> &k7a, std::vector<double> &k8a, std::vector<double> &k9a,
    std::vector<double> &k10a, std::vector<double> &k11a, std::vector<double> &k12a,
    std::vector<double> &k13a, std::vector<std::vector<double>> &k13dda,
    std::vector<double> &k14a, std::vector<double> &k15a, std::vector<double> &k16a,
    std::vector<double> &k17a, std::vector<double> &k18a, std::vector<double> &k19a,
    std::vector<double> &k22a,
    std::vector<double> &k50a, std::vector<double> &k51a, std::vector<double> &k52a,
    std::vector<double> &k53a, std::vector<double> &k54a, std::vector<double> &k55a,
    std::vector<double> &k56a, std::vector<double> &k57a, std::vector<double> &k58a,
    int ndratec, double dtemstart, double dtemend,
    std::vector<std::vector<double>> &h2dusta,
    std::vector<double> &ncrna, std::vector<double> &ncrd1a, std::vector<double> &ncrd2a,
    double avgsighi, double avgsighei, double avgsigheii, double piHI, double piHeI,
    // Output rate values
    std::vector<double> &k1, std::vector<double> &k2, std::vector<double> &k3,
    std::vector<double> &k4, std::vector<double> &k5, std::vector<double> &k6,
    std::vector<double> &k7, std::vector<double> &k8, std::vector<double> &k9,
    std::vector<double> &k10, std::vector<double> &k11, std::vector<double> &k12,
    std::vector<double> &k13, std::vector<double> &k14, std::vector<double> &k15,
    std::vector<double> &k16, std::vector<double> &k17, std::vector<double> &k18,
    std::vector<double> &k19, std::vector<double> &k22,
    double k24, double k25, double k26, double k28, double k29, double k30, double k31,
    std::vector<double> &k50, std::vector<double> &k51, std::vector<double> &k52,
    std::vector<double> &k53, std::vector<double> &k54, std::vector<double> &k55,
    std::vector<double> &k56, std::vector<double> &k57, std::vector<double> &k58,
    std::vector<std::vector<double>> &k13dd, std::vector<double> &h2dust,
    std::vector<double> &ncrn, std::vector<double> &ncrd1, std::vector<double> &ncrd2,
    std::vector<double> &k24shield, std::vector<double> &k25shield,
    std::vector<double> &k26shield, std::vector<double> &k28shield,
    std::vector<double> &k29shield, std::vector<double> &k30shield,
    std::vector<double> &k31shield,
    // 1D temporaries
    std::vector<double> &t1, std::vector<double> &t2, std::vector<double> &tdef,
    std::vector<double> &logtem, std::vector<int64_t> &indixe,
    // Additional parameters
    double dom, double coolunit, double tbase1, double xbase1, double dx_cgs,
    double c_ljeans,
    int iradtrans,
    std::vector<std::vector<std::vector<real8>>> &kdissH2I,
    std::vector<std::vector<std::vector<real8>>> &xH2shield,
    int iH2shieldcustom,
    std::vector<std::vector<std::vector<real8>>> &f_shield_custom,
    std::vector<bool> &itmask
)
{
    // Constants
    const double everg = 1.602176634e-12; // Conversion factor from eV to erg
    const double e24 = 13.6;   // Ionization energy of hydrogen in eV
    const double e26 = 24.6;   // Ionization energy of helium in eV
    const double kboltz = 1.380649e-16; // Boltzmann constant in erg/K
    const double mass_h = 1.6726219e-24; // Mass of hydrogen atom in grams

    // Logarithmic values for the lookup tables
    double logtem0 = std::log(temstart);
    double logtem9 = std::log(temend);
    double dlogtem = (std::log(temend) - std::log(temstart)) / static_cast<double>(nratec - 1);

    // Loop over 'i' from 'is+1' to 'ie+1'
    for (int i = is + 1; i <= ie + 1; ++i) {
        if (itmask[i]) {
            // Compute log of temperature
            logtem[i] = std::log(tgas1d[i]);
            logtem[i] = std::max(logtem[i], logtem0);
            logtem[i] = std::min(logtem[i], logtem9);

            // Find index into table and precompute interpolation values
            indixe[i] = std::min(static_cast<int64_t>(nratec - 1),
                                 std::max(static_cast<int64_t>(1),
                                          static_cast<int64_t>((logtem[i] - logtem0) / dlogtem) + 1));
            t1[i] = logtem0 + (indixe[i] - 1) * dlogtem;
            t2[i] = logtem0 + indixe[i] * dlogtem;
            tdef[i] = (logtem[i] - t1[i]) / (t2[i] - t1[i]);

            // Linear table lookup (in log temperature)
            k1[i] = k1a[indixe[i] - 1] + (k1a[indixe[i]] - k1a[indixe[i] - 1]) * tdef[i];
            k2[i] = k2a[indixe[i] - 1] + (k2a[indixe[i]] - k2a[indixe[i] - 1]) * tdef[i];
            k3[i] = k3a[indixe[i] - 1] + (k3a[indixe[i]] - k3a[indixe[i] - 1]) * tdef[i];
            k4[i] = k4a[indixe[i] - 1] + (k4a[indixe[i]] - k4a[indixe[i] - 1]) * tdef[i];
            k5[i] = k5a[indixe[i] - 1] + (k5a[indixe[i]] - k5a[indixe[i] - 1]) * tdef[i];
            k6[i] = k6a[indixe[i] - 1] + (k6a[indixe[i]] - k6a[indixe[i] - 1]) * tdef[i];
            k57[i] = k57a[indixe[i] - 1] + (k57a[indixe[i]] - k57a[indixe[i] - 1]) * tdef[i];
            k58[i] = k58a[indixe[i] - 1] + (k58a[indixe[i]] - k58a[indixe[i] - 1]) * tdef[i];
        }
    }

    // Look-up for 9-species model
    if (ispecies > 1) {
        for (int i = is + 1; i <= ie + 1; ++i) {
            if (itmask[i]) {
                k7[i] = k7a[indixe[i] - 1] + (k7a[indixe[i]] - k7a[indixe[i] - 1]) * tdef[i];
                k8[i] = k8a[indixe[i] - 1] + (k8a[indixe[i]] - k8a[indixe[i] - 1]) * tdef[i];
                k9[i] = k9a[indixe[i] - 1] + (k9a[indixe[i]] - k9a[indixe[i] - 1]) * tdef[i];
                k10[i] = k10a[indixe[i] - 1] + (k10a[indixe[i]] - k10a[indixe[i] - 1]) * tdef[i];
                k11[i] = k11a[indixe[i] - 1] + (k11a[indixe[i]] - k11a[indixe[i] - 1]) * tdef[i];
                k12[i] = k12a[indixe[i] - 1] + (k12a[indixe[i]] - k12a[indixe[i] - 1]) * tdef[i];
                k13[i] = k13a[indixe[i] - 1] + (k13a[indixe[i]] - k13a[indixe[i] - 1]) * tdef[i];
                k14[i] = k14a[indixe[i] - 1] + (k14a[indixe[i]] - k14a[indixe[i] - 1]) * tdef[i];
                k15[i] = k15a[indixe[i] - 1] + (k15a[indixe[i]] - k15a[indixe[i] - 1]) * tdef[i];
                k16[i] = k16a[indixe[i] - 1] + (k16a[indixe[i]] - k16a[indixe[i] - 1]) * tdef[i];
                k17[i] = k17a[indixe[i] - 1] + (k17a[indixe[i]] - k17a[indixe[i] - 1]) * tdef[i];
                k18[i] = k18a[indixe[i] - 1] + (k18a[indixe[i]] - k18a[indixe[i] - 1]) * tdef[i];
                k19[i] = k19a[indixe[i] - 1] + (k19a[indixe[i]] - k19a[indixe[i] - 1]) * tdef[i];
                k22[i] = k22a[indixe[i] - 1] + (k22a[indixe[i]] - k22a[indixe[i] - 1]) * tdef[i];

                // H2 formation heating terms
                ncrn[i] = ncrna[indixe[i] - 1] + (ncrna[indixe[i]] - ncrna[indixe[i] - 1]) * tdef[i];
                ncrd1[i] = ncrd1a[indixe[i] - 1] + (ncrd1a[indixe[i]] - ncrd1a[indixe[i] - 1]) * tdef[i];
                ncrd2[i] = ncrd2a[indixe[i] - 1] + (ncrd2a[indixe[i]] - ncrd2a[indixe[i] - 1]) * tdef[i];
            }
        }

        // Loop over k13dd array
        for (int n1 = 0; n1 < 14; ++n1) {
            for (int i = is + 1; i <= ie + 1; ++i) {
                if (itmask[i]) {
                    k13dd[i][n1] = k13dda[indixe[i] - 1][n1] +
                                   (k13dda[indixe[i]][n1] - k13dda[indixe[i] - 1][n1]) * tdef[i];
                }
            }
        }
    }

    // Look-up for 12-species model
    if (ispecies > 2) {
        for (int i = is + 1; i <= ie + 1; ++i) {
            if (itmask[i]) {
                k50[i] = k50a[indixe[i] - 1] + (k50a[indixe[i]] - k50a[indixe[i] - 1]) * tdef[i];
                k51[i] = k51a[indixe[i] - 1] + (k51a[indixe[i]] - k51a[indixe[i] - 1]) * tdef[i];
                k52[i] = k52a[indixe[i] - 1] + (k52a[indixe[i]] - k52a[indixe[i] - 1]) * tdef[i];
                k53[i] = k53a[indixe[i] - 1] + (k53a[indixe[i]] - k53a[indixe[i] - 1]) * tdef[i];
                k54[i] = k54a[indixe[i] - 1] + (k54a[indixe[i]] - k54a[indixe[i] - 1]) * tdef[i];
                k55[i] = k55a[indixe[i] - 1] + (k55a[indixe[i]] - k55a[indixe[i] - 1]) * tdef[i];
                k56[i] = k56a[indixe[i] - 1] + (k56a[indixe[i]] - k56a[indixe[i] - 1]) * tdef[i];
            }
        }
    }

    // Look-up for H2 formation on dust
    if (anydust) {
        double d_logtem0 = std::log(dtemstart);
        double d_logtem9 = std::log(dtemend);
        double d_dlogtem = (std::log(dtemend) - std::log(dtemstart)) / static_cast<double>(ndratec - 1);

        std::vector<double> d_logtem(ie + 2, 0.0);
        std::vector<int64_t> d_indixe(ie + 2, 0);
        std::vector<double> d_t1(ie + 2, 0.0);
        std::vector<double> d_t2(ie + 2, 0.0);
        std::vector<double> d_tdef(ie + 2, 0.0);
        std::vector<double> dusti1(ie + 2, 0.0);
        std::vector<double> dusti2(ie + 2, 0.0);

        for (int i = is + 1; i <= ie + 1; ++i) {
            if (itmask[i]) {
                if (tdust[i] > dtemend) {
                    h2dust[i] = tiny8;
                } else {
                    // Get log dust temperature
                    d_logtem[i] = std::log(tdust[i]);
                    d_logtem[i] = std::max(d_logtem[i], d_logtem0);
                    d_logtem[i] = std::min(d_logtem[i], d_logtem9);

                    // Find index into table and precompute interpolation values
                    d_indixe[i] = std::min(static_cast<int64_t>(ndratec - 1),
                                           std::max(static_cast<int64_t>(1),
                                                    static_cast<int64_t>((d_logtem[i] - d_logtem0) / d_dlogtem) + 1));
                    d_t1[i] = d_logtem0 + (d_indixe[i] - 1) * d_dlogtem;
                    d_t2[i] = d_logtem0 + d_indixe[i] * d_dlogtem;
                    d_tdef[i] = (d_logtem[i] - d_t1[i]) / (d_t2[i] - d_t1[i]);

                    // Get rate from 2D interpolation
                    dusti1[i] = h2dusta[indixe[i] - 1][d_indixe[i] - 1] +
                                (h2dusta[indixe[i]][d_indixe[i] - 1] - h2dusta[indixe[i] - 1][d_indixe[i] - 1]) * tdef[i];
                    dusti2[i] = h2dusta[indixe[i] - 1][d_indixe[i]] +
                                (h2dusta[indixe[i]][d_indixe[i]] - h2dusta[indixe[i] - 1][d_indixe[i]]) * tdef[i];
                    h2dust[i] = dusti1[i] + (dusti2[i] - dusti1[i]) * d_tdef[i];

                    // Multiply by dust to gas ratio
                    h2dust[i] *= dust2gas[i];
                }
            }
        }
    }

    // Initialize shielding factors
    for (int i = is + 1; i <= ie + 1; ++i) {
        if (itmask[i]) {
            k24shield[i] = k24;
            k25shield[i] = k25;
            k26shield[i] = k26;
            k28shield[i] = k28;
            k29shield[i] = k29;
            k30shield[i] = k30;
        }
    }

    // Inside the lookup_cool_rates1d_g function, after previous calculations

    if (ispecies > 1) {
        if (iradtrans == 0) {
            for (int i = is + 1; i <= ie + 1; ++i) {
                if (itmask[i]) {
                    k31shield[i] = k31;
                }
            }
        } else {
            for (int i = is + 1; i <= ie + 1; ++i) {
                if (itmask[i]) {
                    k31shield[i] = k31 + kdissH2I[i][j][k];
                }
            }
        }
    }


    if (iH2shield > 0) {
        for (int i = is + 1; i <= ie + 1; ++i) {
            if (itmask[i]) {
                double l_H2shield = 0.0;

                if (iH2shield == 1) {
                    // Calculate a Sobolev-like length assuming a 3D grid
                    std::vector<double> divrhoa(6);
                    double divrho = tiny;
                    int ip1 = std::min(i + 1, in - 1);??? che cazzo sono queste variabili???
                    int im1 = std::max(i - 1, 0);
                    int jp1 = std::min(j + 1, jn - 1);
                    int jm1 = std::max(j - 1, 0);
                    int kp1 = std::min(k + 1, kn - 1);
                    int km1 = std::max(k - 1, 0);

                    divrhoa[0] = d[ip1][j][k] - d[i][j][k];
                    divrhoa[1] = d[im1][j][k] - d[i][j][k];
                    divrhoa[2] = d[i][jp1][k] - d[i][j][k];
                    divrhoa[3] = d[i][jm1][k] - d[i][j][k];
                    divrhoa[4] = d[i][j][kp1] - d[i][j][k];
                    divrhoa[5] = d[i][j][km1] - d[i][j][k];

                    // Exclude directions with (drho/ds > 0)
                    for (int n1 = 0; n1 < 6; ++n1) {
                        if (divrhoa[n1] < 0.0) {
                            divrho += divrhoa[n1];
                        }
                    }

                    // Calculate Sobolev-like length
                    l_H2shield = std::min(dx_cgs * d[i][j][k] / std::abs(divrho), xbase1);
                } else if (iH2shield == 2) {
                    // User-supplied length-scale field
                    l_H2shield = xH2shield[i][j][k] * xbase1;
                } else if (iH2shield == 3) {
                    // Jeans Length
                    l_H2shield = c_ljeans * std::sqrt(tgas1d[i] / (d[i][j][k] * mmw[i]));
                } else {
                    l_H2shield = 0.0;
                }

                // Calculate N_H2
                double N_H2 = dom * H2I[i][j][k] * l_H2shield;

                // Self-shielding following Wolcott-Green & Haiman (2019)
                double tgas_touse = std::clamp(tgas1d[i], 1.0e2, 8.0e3);
                double ngas_touse = std::min(d[i][j][k] * dom / mmw[i], 1.0e7);

                double aWG2019 = (0.8711 * std::log10(tgas_touse) - 1.928) *
                                 std::exp(-0.2856 * std::log10(ngas_touse)) +
                                 (-0.9639 * std::log10(tgas_touse) + 3.892);

                double x = 2.0e-15 * N_H2;
                double b_doppler = 1.0e-5 * std::sqrt(2.0 * kboltz * tgas1d[i] / mass_h);

                double f_shield = 0.965 / std::pow(1.0 + x / b_doppler, aWG2019) +
                                  0.035 * std::exp(-8.5e-4 * std::sqrt(1.0 + x)) / std::sqrt(1.0 + x);

                // Ensure f_shield <= 1
                f_shield = std::min(f_shield, 1.0);

                // Update k31shield[i]
                k31shield[i] *= f_shield;
            }
        }
    }
    

    if (iH2shieldcustom > 0) {
        for (int i = is + 1; i <= ie + 1; ++i) {
            if (itmask[i]) {
                k31shield[i] *= f_shield_custom[i][j][k];
            }
        }
    }
    

    if (iradshield > 0) {
        // Compute shielding factors
        for (int i = is + 1; i <= ie + 1; ++i) {
            if (itmask[i]) {

                // Compute shielding factor for H
                double nSSh = 6.73e-3 *
                    std::pow(avgsighi / 2.49e-18, -2.0 / 3.0) *
                    std::pow(tgas1d[i] / 1.0e4, 0.17) *
                    std::pow(k24 / tbase1 / 1.0e-12, 2.0 / 3.0);

                // Compute the total Hydrogen number density
                double nratio = HI[i][j][k] + HII[i][j][k];
                if (ispecies > 1) {
                    nratio += HM[i][j][k] + H2I[i][j][k] + H2II[i][j][k];

                    if (ispecies > 2) {
                        nratio += 0.5 * (DI[i][j][k] + DII[i][j][k]) +
                                  (2.0 / 3.0) * HDI[i][j][k];
                    }
                }

                nratio = nratio * dom / nSSh;

                f_shield_H[i] = (0.98 *
                    std::pow(1.0 + std::pow(nratio, 1.64), -2.28) +
                    0.02 * std::pow(1.0 + nratio, -0.84));

                // Compute shielding factor for He

                nSSh = 6.73e-3 *
                    std::pow(avgsighei / 2.49e-18, -2.0 / 3.0) *
                    std::pow(tgas1d[i] / 1.0e4, 0.17) *
                    std::pow(k26 / tbase1 / 1.0e-12, 2.0 / 3.0);

                nratio = 0.25 *
                    (HeI[i][j][k] + HeII[i][j][k] + HeIII[i][j][k]) * dom / nSSh;

                f_shield_He[i] = (0.98 *
                    std::pow(1.0 + std::pow(nratio, 1.64), -2.28) +
                    0.02 * std::pow(1.0 + nratio, -0.84));

            }
        }
    }
    
    // Inside the lookup_cool_rates1d_g function, after previous calculations

    if (iradshield == 1) {
        // Approximate self-shielding using equations from Rahmati et al. 2013
        // Shield HI while leaving HeI and HeII optically thin
        // Attenuate radiation rates for direct H2 ionization using the same scaling (rate k29)
        for (int i = is + 1; i <= ie + 1; ++i) {
            if (itmask[i]) {
                if (k24 < tiny8) {
                    k24shield[i] = 0.0;
                } else {
                    k24shield[i] *= f_shield_H[i];
                }

                if (k29 < tiny8) {
                    k29shield[i] = 0.0;
                } else {
                    k29shield[i] *= f_shield_H[i];
                }

                k25shield[i] = k25;
                k26shield[i] = k26;
            }
        }
    } else if (iradshield == 2) {
        // Better self-shielding in HI and approximate self-shielding in HeI and HeII
        // Attenuate radiation rates for direct H2 ionization (k29) and H2+ dissociation (k28 and k30)
        for (int i = is + 1; i <= ie + 1; ++i) {
            if (itmask[i]) {
                if (k24 < tiny8) {
                    k24shield[i] = 0.0;
                } else {
                    k24shield[i] *= f_shield_H[i];
                }

                if (k29 < tiny8) {
                    k29shield[i] = 0.0;
                } else {
                    k29shield[i] *= f_shield_H[i];
                }

                if (k26 < tiny8) {
                    k26shield[i] = 0.0;
                } else {
                    k26shield[i] *= f_shield_He[i];
                }

                if (k28 < tiny8) {
                    k28shield[i] = 0.0;
                } else {
                    k28shield[i] *= f_shield_He[i];
                }

                if (k30 < tiny8) {
                    k30shield[i] = 0.0;
                } else {
                    k30shield[i] *= f_shield_He[i];
                }

                k25shield[i] = k25;
            }
        }
    } else if (iradshield == 3) {
        // Shielding in HI and HeI, ignoring HeII heating entirely
        for (int i = is + 1; i <= ie + 1; ++i) {
            if (itmask[i]) {
                if (k24 < tiny8) {
                    k24shield[i] = 0.0;
                } else {
                    k24shield[i] *= f_shield_H[i];
                }

                if (k29 < tiny8) {
                    k29shield[i] = 0.0;
                } else {
                    k29shield[i] *= f_shield_H[i];
                }

                if (k26 < tiny8) {
                    k26shield[i] = 0.0;
                } else {
                    k26shield[i] *= f_shield_He[i];
                }

                if (k28 < tiny8) {
                    k28shield[i] = 0.0;
                } else {
                    k28shield[i] *= f_shield_He[i];
                }

                if (k30 < tiny8) {
                    k30shield[i] = 0.0;
                } else {
                    k30shield[i] *= f_shield_He[i];
                }

                k25shield[i] = 0.0;
            }
        }
    }
    

    // If secondary ionization is considered
    #ifdef SECONDARY_IONIZATION_NOT_YET_IMPLEMENTED
        // If using a high-energy radiation field, account for effects of secondary electrons
        // (Shull & Steenberg 1985)
        for (int i = is + 1; i <= ie + 1; ++i) {
            if (itmask[i]) {
                double x = std::max(HII[i][j][k] / (HI[i][j][k] + HII[i][j][k]), 1.0e-4);
                double factor = 0.3908 * std::pow(1.0 - std::pow(x, 0.4092), 1.7592);
                k24shield[i] += factor * (piHI + 0.08 * piHeI) / (e24 * everg) * coolunit * tbase1;
                factor = 0.0554 * std::pow(1.0 - std::pow(x, 0.4614), 1.6660);
                k26shield[i] += factor * (piHI / 0.08 + piHeI) / (e26 * everg) * coolunit * tbase1;
            }
        }
    #endif

    // If using the density-dependent H2 dissociation rate
    #ifdef USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE
        if (ispecies > 1 && ithreebody == 0) {
            for (int i = is + 1; i <= ie + 1; ++i) {
                if (itmask[i]) {
                    double nh = std::min(HI[i][j][k] * dom, 1.0e9);
                    k13[i] = tiny8;
                    if (tgas1d[i] >= 500.0 && tgas1d[i] < 1.0e6) {
                        // Direct collisional dissociation
                        double k13_CID = k13dd[i][0] - k13dd[i][1] /
                            (1.0 + std::pow(nh / k13dd[i][4], k13dd[i][6]))
                            + k13dd[i][2] - k13dd[i][3] /
                            (1.0 + std::pow(nh / k13dd[i][5], k13dd[i][6]));
                        k13_CID = std::max(std::pow(10.0, k13_CID), tiny8);
                        // Dissociative tunneling
                        double k13_DT = k13dd[i][7] - k13dd[i][8] /
                            (1.0 + std::pow(nh / k13dd[i][11], k13dd[i][13]))
                            + k13dd[i][9] - k13dd[i][10] /
                            (1.0 + std::pow(nh / k13dd[i][12], k13dd[i][13]));
                        k13_DT = std::max(std::pow(10.0, k13_DT), tiny8);
                        // Total rate
                        k13[i] = k13_DT + k13_CID;
                    }
                }
            }
        }
    #endif
    
}

#include <vector>
#include <cmath>
#include <algorithm>

// Define types for consistency with Fortran's real*8
typedef double real8;
typedef int64_t int8;

void rate_timestep_g(
    // Outputs
    std::vector<double>& dedot, // Electron density rate-of-change
    std::vector<double>& HIdot, // HI density rate-of-change

    // Inputs
    int ispecies, // Species model identifier
    bool anydust, // Flag indicating the presence of dust

    // Density fields
    const std::vector<std::vector<std::vector<double>>>& de,   // Electron density
    const std::vector<std::vector<std::vector<double>>>& HI,   // Neutral hydrogen density
    const std::vector<std::vector<std::vector<double>>>& HII,  // Ionized hydrogen density
    const std::vector<std::vector<std::vector<double>>>& HeI,  // Neutral helium density
    const std::vector<std::vector<std::vector<double>>>& HeII, // Singly ionized helium density
    const std::vector<std::vector<std::vector<double>>>& HeIII,// Doubly ionized helium density
    const std::vector<std::vector<std::vector<double>>>& d,    // Total density

    const std::vector<std::vector<std::vector<double>>>& HM,   // Negative hydrogen ion density
    const std::vector<std::vector<std::vector<double>>>& H2I,  // Molecular hydrogen density
    const std::vector<std::vector<std::vector<double>>>& H2II, // Ionized molecular hydrogen density

    // Grid parameters
    int in, int jn, int kn, // Grid dimensions
    int is, int ie,         // Start and end indices in i-direction
    int j, int k,           // Fixed indices in j and k-directions

    // Rate coefficients (arrays of size 'in')
    const std::vector<double>& k1, const std::vector<double>& k2,
    const std::vector<double>& k3, const std::vector<double>& k4,
    const std::vector<double>& k5, const std::vector<double>& k6,
    const std::vector<double>& k7, const std::vector<double>& k8,
    const std::vector<double>& k9, const std::vector<double>& k10,
    const std::vector<double>& k11, const std::vector<double>& k12,
    const std::vector<double>& k13, const std::vector<double>& k14,
    const std::vector<double>& k15, const std::vector<double>& k16,
    const std::vector<double>& k17, const std::vector<double>& k18,
    const std::vector<double>& k19, const std::vector<double>& k22,

    double k24, double k25, double k26, double k27, double k28, double k29, double k30,

    const std::vector<double>& k50, const std::vector<double>& k51,
    const std::vector<double>& k52, const std::vector<double>& k53,
    const std::vector<double>& k54, const std::vector<double>& k55,
    const std::vector<double>& k56, const std::vector<double>& k57,
    const std::vector<double>& k58,

    const std::vector<double>& h2dust, const std::vector<double>& ncrn,
    const std::vector<double>& ncrd1, const std::vector<double>& ncrd2,
    const std::vector<double>& rhoH,

    const std::vector<double>& k24shield, const std::vector<double>& k25shield,
    const std::vector<double>& k26shield, const std::vector<double>& k28shield,
    const std::vector<double>& k29shield, const std::vector<double>& k30shield,
    const std::vector<double>& k31shield,

    int iradtrans, int irt_honly,

    // Radiative Transfer Fields
    const std::vector<std::vector<std::vector<double>>>& kphHI,
    const std::vector<std::vector<std::vector<double>>>& kphHeI,
    const std::vector<std::vector<std::vector<double>>>& kphHeII,

    // Mask and other parameters
    const std::vector<bool>& itmask, // Mask indicating active cells
    std::vector<double>& edot,       // Energy rate-of-change
    double chunit, double dom        // Unit conversion factors
)
{
    // Constants
    const double tiny = 1e-20;

    // Local variables
    std::vector<double> h2heatfac(in, 0.0);
    std::vector<double> H2delta(in, 0.0);
    double H2dmag, atten, tau;

    if (ispecies == 1) {
        for (int i = is + 1; i <= ie + 1; ++i) {
            if (itmask[i]) {
                // Compute the electron density rate-of-change
                dedot[i] =
                    + k1[i] * HI[i][j][k] * de[i][j][k]
                    + k3[i] * HeI[i][j][k] * de[i][j][k] / 4.0
                    + k5[i] * HeII[i][j][k] * de[i][j][k] / 4.0
                    - k2[i] * HII[i][j][k] * de[i][j][k]
                    - k4[i] * HeII[i][j][k] * de[i][j][k] / 4.0
                    - k6[i] * HeIII[i][j][k] * de[i][j][k] / 4.0
                    + k57[i] * HI[i][j][k] * HI[i][j][k]
                    + k58[i] * HI[i][j][k] * HeI[i][j][k] / 4.0
                    + (k24shield[i] * HI[i][j][k]
                    + k25shield[i] * HeII[i][j][k] / 4.0
                    + k26shield[i] * HeI[i][j][k] / 4.0);

                // Compute the HI density rate-of-change
                HIdot[i] =
                    - k1[i] * HI[i][j][k] * de[i][j][k]
                    + k2[i] * HII[i][j][k] * de[i][j][k]
                    - k57[i] * HI[i][j][k] * HI[i][j][k]
                    - k58[i] * HI[i][j][k] * HeI[i][j][k] / 4.0
                    - k24shield[i] * HI[i][j][k];
            }
        }
    } else {
        // Include molecular hydrogen rates for HIdot
        for (int i = is + 1; i <= ie + 1; ++i) {
            if (itmask[i]) {
                HIdot[i] =
                    - k1[i] * de[i][j][k] * HI[i][j][k]
                    - k7[i] * de[i][j][k] * HI[i][j][k]
                    - k8[i] * HM[i][j][k] * HI[i][j][k]
                    - k9[i] * HII[i][j][k] * HI[i][j][k]
                    - k10[i] * H2II[i][j][k] * HI[i][j][k] / 2.0
                    - 2.0 * k22[i] * std::pow(HI[i][j][k], 3)
                    + k2[i] * HII[i][j][k] * de[i][j][k]
                    + 2.0 * k13[i] * HI[i][j][k] * H2I[i][j][k] / 2.0
                    + k11[i] * HII[i][j][k] * H2I[i][j][k] / 2.0
                    + 2.0 * k12[i] * de[i][j][k] * H2I[i][j][k] / 2.0
                    + k14[i] * HM[i][j][k] * de[i][j][k]
                    + k15[i] * HM[i][j][k] * HI[i][j][k]
                    + 2.0 * k16[i] * HM[i][j][k] * HII[i][j][k]
                    + 2.0 * k18[i] * H2II[i][j][k] * de[i][j][k] / 2.0
                    + k19[i] * H2II[i][j][k] * HM[i][j][k] / 2.0
                    - k57[i] * HI[i][j][k] * HI[i][j][k]
                    - k58[i] * HI[i][j][k] * HeI[i][j][k] / 4.0
                    - k24shield[i] * HI[i][j][k]
                    + 2.0 * k31shield[i] * H2I[i][j][k] / 2.0;

                // Add H2 formation on dust grains
                if (anydust) {
                    HIdot[i] -= 2.0 * h2dust[i] * rhoH[i];
                }

                // Compute the electron density rate-of-change
                dedot[i] =
                    + k1[i] * HI[i][j][k] * de[i][j][k]
                    + k3[i] * HeI[i][j][k] * de[i][j][k] / 4.0
                    + k5[i] * HeII[i][j][k] * de[i][j][k] / 4.0
                    + k8[i] * HM[i][j][k] * HI[i][j][k]
                    + k15[i] * HM[i][j][k] * HI[i][j][k]
                    + k17[i] * HM[i][j][k] * HII[i][j][k]
                    + k14[i] * HM[i][j][k] * de[i][j][k]
                    - k2[i] * HII[i][j][k] * de[i][j][k]
                    - k4[i] * HeII[i][j][k] * de[i][j][k] / 4.0
                    - k6[i] * HeIII[i][j][k] * de[i][j][k] / 4.0
                    - k7[i] * HI[i][j][k] * de[i][j][k]
                    - k18[i] * H2II[i][j][k] * de[i][j][k] / 2.0
                    + k57[i] * HI[i][j][k] * HI[i][j][k]
                    + k58[i] * HI[i][j][k] * HeI[i][j][k] / 4.0
                    + (k24shield[i] * HI[i][j][k]
                    + k25shield[i] * HeII[i][j][k] / 4.0
                    + k26shield[i] * HeI[i][j][k] / 4.0);

                // H2 formation heating
                h2heatfac[i] = std::pow(1.0 + (ncrn[i] / (dom *
                    (HI[i][j][k] * ncrd1[i] +
                    H2I[i][j][k] * 0.5 * ncrd2[i]))), -1.0);

                H2delta[i] =
                    HI[i][j][k] *
                    (4.48 * k22[i] * std::pow(HI[i][j][k], 2.0)
                     - 4.48 * k13[i] * H2I[i][j][k] / 2.0);

                // Apply heating factor if formation dominates
                if (H2delta[i] > 0.0) {
                    H2delta[i] *= h2heatfac[i];
                }

                if (anydust) {
                    H2delta[i] += h2dust[i] * HI[i][j][k] * rhoH[i] *
                        (0.2 + 4.2 * h2heatfac[i]);
                }

                atten = 1.0;
                edot[i] += chunit * H2delta[i] * atten;
            }
        }
    }

    // Add photo-ionization rates if needed
    if (iradtrans == 1) {
        if (irt_honly == 0) {
            for (int i = is + 1; i <= ie + 1; ++i) {
                if (itmask[i]) {
                    HIdot[i] -= kphHI[i][j][k] * HI[i][j][k];
                    dedot[i] += kphHI[i][j][k] * HI[i][j][k]
                        + kphHeI[i][j][k] * HeI[i][j][k] / 4.0
                        + kphHeII[i][j][k] * HeII[i][j][k] / 4.0;
                }
            }
        } else {
            for (int i = is + 1; i <= ie + 1; ++i) {
                if (itmask[i]) {
                    HIdot[i] -= kphHI[i][j][k] * HI[i][j][k];
                    dedot[i] += kphHI[i][j][k] * HI[i][j][k];
                }
            }
        }
    }
}

#include <vector>
#include <cmath>
#include <algorithm>

// Function to advance rate equations by one sub-cycle using backward-Euler method
void step_rate_g(
    // Density fields (inputs and outputs)
    std::vector<std::vector<std::vector<double>>>& de,
    std::vector<std::vector<std::vector<double>>>& HI,
    std::vector<std::vector<std::vector<double>>>& HII,
    std::vector<std::vector<std::vector<double>>>& HeI,
    std::vector<std::vector<std::vector<double>>>& HeII,
    std::vector<std::vector<std::vector<double>>>& HeIII,
    const std::vector<std::vector<std::vector<double>>>& d,

    std::vector<std::vector<std::vector<double>>>& HM,
    std::vector<std::vector<std::vector<double>>>& H2I,
    std::vector<std::vector<std::vector<double>>>& H2II,
    std::vector<std::vector<std::vector<double>>>& DI,
    std::vector<std::vector<std::vector<double>>>& DII,
    std::vector<std::vector<std::vector<double>>>& HDI,

    // Time step
    const std::vector<double>& dtit,

    // Grid parameters
    int in, int jn, int kn, int is, int ie, int j, int k,

    // Species model and dust presence
    int ispecies, bool anydust,

    // Rate coefficients (arrays of size 'in')
    const std::vector<double>& k1, const std::vector<double>& k2,
    const std::vector<double>& k3, const std::vector<double>& k4,
    const std::vector<double>& k5, const std::vector<double>& k6,
    const std::vector<double>& k7, const std::vector<double>& k8,
    const std::vector<double>& k9, const std::vector<double>& k10,
    const std::vector<double>& k11, const std::vector<double>& k12,
    const std::vector<double>& k13, const std::vector<double>& k14,
    const std::vector<double>& k15, const std::vector<double>& k16,
    const std::vector<double>& k17, const std::vector<double>& k18,
    const std::vector<double>& k19, const std::vector<double>& k22,

    double k24, double k25, double k26, double k27, double k28, double k29, double k30,

    const std::vector<double>& k50, const std::vector<double>& k51,
    const std::vector<double>& k52, const std::vector<double>& k53,
    const std::vector<double>& k54, const std::vector<double>& k55,
    const std::vector<double>& k56, const std::vector<double>& k57,
    const std::vector<double>& k58,

    const std::vector<double>& h2dust, const std::vector<double>& rhoH,

    const std::vector<double>& k24shield, const std::vector<double>& k25shield,
    const std::vector<double>& k26shield, const std::vector<double>& k28shield,
    const std::vector<double>& k29shield, const std::vector<double>& k30shield,
    const std::vector<double>& k31shield,

    // Temporary variables (passed in and updated)
    std::vector<double>& HIp, std::vector<double>& HIIp,
    std::vector<double>& HeIp, std::vector<double>& HeIIp,
    std::vector<double>& HeIIIp, std::vector<double>& dep,

    std::vector<double>& HMp, std::vector<double>& H2Ip, std::vector<double>& H2IIp,
    std::vector<double>& DIp, std::vector<double>& DIIp, std::vector<double>& HDIp,

    std::vector<double>& dedot_prev, std::vector<double>& HIdot_prev,

    // Radiative transfer parameters
    int iradtrans, int irt_honly,
    const std::vector<std::vector<std::vector<double>>>& kphHI,
    const std::vector<std::vector<std::vector<double>>>& kphHeI,
    const std::vector<std::vector<std::vector<double>>>& kphHeII,

    // Mask
    const std::vector<bool>& itmask
) {
    // Constants
    const double tiny = 1e-20;
    const double tiny8 = 1e-20;
    const double RKIND = 1.0;
    const double DKIND = 1.0;

    // Local variables
    double scoef, acoef;

    // A) The 6-species integrator
    if (ispecies == 1) {
        for (int i = is + 1; i <= ie + 1; ++i) {
            if (itmask[i]) {
                // 1) HI
                scoef = k2[i] * HII[i][j][k] * de[i][j][k];
                acoef = k1[i] * de[i][j][k]
                    + k57[i] * HI[i][j][k]
                    + k58[i] * HeI[i][j][k] / 4.0
                    + k24shield[i];
                if (iradtrans == 1) acoef += kphHI[i][j][k];
                HIp[i] = (scoef * dtit[i] + HI[i][j][k]) / (1.0 + acoef * dtit[i]);

                // 2) HII
                scoef = k1[i] * HIp[i] * de[i][j][k]
                    + k57[i] * HIp[i] * HIp[i]
                    + k58[i] * HIp[i] * HeI[i][j][k] / 4.0
                    + k24shield[i] * HIp[i];
                if (iradtrans == 1)
                    scoef += kphHI[i][j][k] * HIp[i];
                acoef = k2[i] * de[i][j][k];
                HIIp[i] = (scoef * dtit[i] + HII[i][j][k]) / (1.0 + acoef * dtit[i]);

                // 3) Electron density
                scoef = 0.0
                    + k57[i] * HIp[i] * HIp[i]
                    + k58[i] * HIp[i] * HeI[i][j][k] / 4.0
                    + k24shield[i] * HI[i][j][k]
                    + k25shield[i] * HeII[i][j][k] / 4.0
                    + k26shield[i] * HeI[i][j][k] / 4.0;

                if (iradtrans == 1 && irt_honly == 0)
                    scoef += kphHI[i][j][k] * HI[i][j][k]
                        + kphHeI[i][j][k] * HeI[i][j][k] / 4.0
                        + kphHeII[i][j][k] * HeII[i][j][k] / 4.0;
                else if (iradtrans == 1 && irt_honly == 1)
                    scoef += kphHI[i][j][k] * HI[i][j][k];

                acoef = - (k1[i] * HI[i][j][k] - k2[i] * HII[i][j][k]
                    + k3[i] * HeI[i][j][k] / 4.0 - k6[i] * HeIII[i][j][k] / 4.0
                    + k5[i] * HeII[i][j][k] / 4.0 - k4[i] * HeII[i][j][k] / 4.0);

                dep[i] = (scoef * dtit[i] + de[i][j][k]) / (1.0 + acoef * dtit[i]);
            }
        }
    }

    // B) Do helium chemistry in any case (for all ispecies values)
    for (int i = is + 1; i <= ie + 1; ++i) {
        if (itmask[i]) {
            // 4) HeI
            scoef = k4[i] * HeII[i][j][k] * de[i][j][k];
            acoef = k3[i] * de[i][j][k] + k26shield[i];
            if (iradtrans == 1 && irt_honly == 0)
                acoef += kphHeI[i][j][k];
            HeIp[i] = (scoef * dtit[i] + HeI[i][j][k]) / (1.0 + acoef * dtit[i]);

            // 5) HeII
            scoef = k3[i] * HeIp[i] * de[i][j][k]
                + k6[i] * HeIII[i][j][k] * de[i][j][k]
                + k26shield[i] * HeIp[i];
            if (iradtrans == 1 && irt_honly == 0)
                scoef += kphHeI[i][j][k] * HeIp[i];

            acoef = k4[i] * de[i][j][k] + k5[i] * de[i][j][k] + k25shield[i];
            if (iradtrans == 1 && irt_honly == 0)
                acoef += kphHeII[i][j][k];

            HeIIp[i] = (scoef * dtit[i] + HeII[i][j][k]) / (1.0 + acoef * dtit[i]);

            // 6) HeIII
            scoef = k5[i] * HeIIp[i] * de[i][j][k] + k25shield[i] * HeIIp[i];
            if (iradtrans == 1 && irt_honly == 0)
                scoef += kphHeII[i][j][k] * HeIIp[i];
            acoef = k6[i] * de[i][j][k];
            HeIIIp[i] = (scoef * dtit[i] + HeIII[i][j][k]) / (1.0 + acoef * dtit[i]);
        }
    }

    // C) Now do extra species for molecular hydrogen
    if (ispecies > 1) {
        // First, do HI/HII with molecular hydrogen terms
        for (int i = is + 1; i <= ie + 1; ++i) {
            if (itmask[i]) {
                // 1) HI
                scoef = k2[i] * HII[i][j][k] * de[i][j][k]
                    + 2.0 * k13[i] * HI[i][j][k] * H2I[i][j][k] / 2.0
                    + k11[i] * HII[i][j][k] * H2I[i][j][k] / 2.0
                    + 2.0 * k12[i] * de[i][j][k] * H2I[i][j][k] / 2.0
                    + k14[i] * HM[i][j][k] * de[i][j][k]
                    + k15[i] * HM[i][j][k] * HI[i][j][k]
                    + 2.0 * k16[i] * HM[i][j][k] * HII[i][j][k]
                    + 2.0 * k18[i] * H2II[i][j][k] * de[i][j][k] / 2.0
                    + k19[i] * H2II[i][j][k] * HM[i][j][k] / 2.0
                    + 2.0 * k31shield[i] * H2I[i][j][k] / 2.0;

                acoef = k1[i] * de[i][j][k]
                    + k7[i] * de[i][j][k]
                    + k8[i] * HM[i][j][k]
                    + k9[i] * HII[i][j][k]
                    + k10[i] * H2II[i][j][k] / 2.0
                    + 2.0 * k22[i] * std::pow(HI[i][j][k], 2)
                    + k57[i] * HI[i][j][k]
                    + k58[i] * HeI[i][j][k] / 4.0
                    + k24shield[i];

                if (iradtrans == 1) acoef += kphHI[i][j][k];
                if (anydust)
                    acoef += 2.0 * h2dust[i] * rhoH[i];

                HIp[i] = (scoef * dtit[i] + HI[i][j][k]) / (1.0 + acoef * dtit[i]);

                // 2) HII
                scoef = k1[i] * HI[i][j][k] * de[i][j][k]
                    + k10[i] * H2II[i][j][k] * HI[i][j][k] / 2.0
                    + k57[i] * HI[i][j][k] * HI[i][j][k]
                    + k58[i] * HI[i][j][k] * HeI[i][j][k] / 4.0
                    + k24shield[i] * HI[i][j][k];
                if (iradtrans == 1)
                    scoef += kphHI[i][j][k] * HI[i][j][k];

                acoef = k2[i] * de[i][j][k]
                    + k9[i] * HI[i][j][k]
                    + k11[i] * H2I[i][j][k] / 2.0
                    + k16[i] * HM[i][j][k]
                    + k17[i] * HM[i][j][k];

                HIIp[i] = (scoef * dtit[i] + HII[i][j][k]) / (1.0 + acoef * dtit[i]);

                // 3) Electrons
                scoef = k8[i] * HM[i][j][k] * HI[i][j][k]
                    + k15[i] * HM[i][j][k] * HI[i][j][k]
                    + k17[i] * HM[i][j][k] * HII[i][j][k]
                    + k57[i] * HI[i][j][k] * HI[i][j][k]
                    + k58[i] * HI[i][j][k] * HeI[i][j][k] / 4.0
                    + k24shield[i] * HIp[i]
                    + k25shield[i] * HeIIp[i] / 4.0
                    + k26shield[i] * HeIp[i] / 4.0;

                if (iradtrans == 1 && irt_honly == 0)
                    scoef += kphHI[i][j][k] * HIp[i]
                        + kphHeI[i][j][k] * HeIp[i] / 4.0
                        + kphHeII[i][j][k] * HeIIp[i] / 4.0;
                else if (iradtrans == 1 && irt_honly == 1)
                    scoef += kphHI[i][j][k] * HIp[i];

                acoef = - (k1[i] * HI[i][j][k] - k2[i] * HII[i][j][k]
                    + k3[i] * HeI[i][j][k] / 4.0 - k6[i] * HeIII[i][j][k] / 4.0
                    + k5[i] * HeII[i][j][k] / 4.0 - k4[i] * HeII[i][j][k] / 4.0
                    + k14[i] * HM[i][j][k]
                    - k7[i] * HI[i][j][k]
                    - k18[i] * H2II[i][j][k] / 2.0);

                dep[i] = (scoef * dtit[i] + de[i][j][k]) / (1.0 + acoef * dtit[i]);

                // 7) H2I
                scoef = 2.0 * (k8[i] * HM[i][j][k] * HI[i][j][k]
                    + k10[i] * H2II[i][j][k] * HI[i][j][k] / 2.0
                    + k19[i] * H2II[i][j][k] * HM[i][j][k] / 2.0
                    + k22[i] * HI[i][j][k] * std::pow(HI[i][j][k], 2.0));

                acoef = k13[i] * HI[i][j][k] + k11[i] * HII[i][j][k]
                    + k12[i] * de[i][j][k]
                    + k29shield[i] + k31shield[i];

                if (anydust)
                    scoef += 2.0 * h2dust[i] * HI[i][j][k] * rhoH[i];

                H2Ip[i] = (scoef * dtit[i] + H2I[i][j][k]) / (1.0 + acoef * dtit[i]);

                // 8) HM
                scoef = k7[i] * HI[i][j][k] * de[i][j][k];
                acoef = (k8[i] + k15[i]) * HI[i][j][k]
                    + (k16[i] + k17[i]) * HII[i][j][k]
                    + k14[i] * de[i][j][k]
                    + k19[i] * H2II[i][j][k] / 2.0
                    + k27;
                HMp[i] = (scoef * dtit[i] + HM[i][j][k]) / (1.0 + acoef * dtit[i]);

                // 9) H2II
                H2IIp[i] = 2.0 * (k9[i] * HIp[i] * HIIp[i]
                    + k11[i] * H2Ip[i] / 2.0 * HIIp[i]
                    + k17[i] * HMp[i] * HIIp[i]
                    + k29shield[i] * H2Ip[i])
                    / (k10[i] * HIp[i] + k18[i] * dep[i]
                        + k19[i] * HMp[i]
                        + (k28shield[i] + k30shield[i]));
            }
        }
    }

    // D) Now do extra species for molecular HD
    if (ispecies > 2) {
        for (int i = is + 1; i <= ie + 1; ++i) {
            if (itmask[i]) {
                // 1) DI
                scoef = k2[i] * DII[i][j][k] * de[i][j][k]
                    + k51[i] * DII[i][j][k] * HI[i][j][k]
                    + 2.0 * k55[i] * HDI[i][j][k] * HI[i][j][k] / 3.0;

                acoef = k1[i] * de[i][j][k]
                    + k50[i] * HII[i][j][k]
                    + k54[i] * H2I[i][j][k] / 2.0
                    + k56[i] * HM[i][j][k]
                    + k24shield[i];
                if (iradtrans == 1) acoef += kphHI[i][j][k];

                DIp[i] = (scoef * dtit[i] + DI[i][j][k]) / (1.0 + acoef * dtit[i]);

                // 2) DII
                scoef = k1[i] * DI[i][j][k] * de[i][j][k]
                    + k50[i] * HII[i][j][k] * DI[i][j][k]
                    + 2.0 * k53[i] * HII[i][j][k] * HDI[i][j][k] / 3.0
                    + k24shield[i] * DI[i][j][k];
                if (iradtrans == 1)
                    scoef += kphHI[i][j][k] * DI[i][j][k];

                acoef = k2[i] * de[i][j][k]
                    + k51[i] * HI[i][j][k]
                    + k52[i] * H2I[i][j][k] / 2.0;

                DIIp[i] = (scoef * dtit[i] + DII[i][j][k]) / (1.0 + acoef * dtit[i]);

                // 3) HDI
                scoef = 3.0 * (k52[i] * DII[i][j][k] * H2I[i][j][k] / 4.0
                    + k54[i] * DI[i][j][k] * H2I[i][j][k] / 4.0
                    + 2.0 * k56[i] * DI[i][j][k] * HM[i][j][k] / 2.0);

                acoef = k53[i] * HII[i][j][k]
                    + k55[i] * HI[i][j][k];

                HDIp[i] = (scoef * dtit[i] + HDI[i][j][k]) / (1.0 + acoef * dtit[i]);
            }
        }
    }

    // E) Set densities from 1D temps to 3D fields
    for (int i = is + 1; i <= ie + 1; ++i) {
        if (itmask[i]) {
            HIdot_prev[i] = std::abs(HI[i][j][k] - HIp[i]) / std::max(dtit[i], tiny8);
            HI[i][j][k] = std::max(HIp[i], tiny);
            HII[i][j][k] = std::max(HIIp[i], tiny);
            HeI[i][j][k] = std::max(HeIp[i], tiny);
            HeII[i][j][k] = std::max(HeIIp[i], tiny);
            HeIII[i][j][k] = std::max(HeIIIp[i], 1e-5 * tiny);

            // Use charge conservation to determine electron fraction
            dedot_prev[i] = de[i][j][k];
            de[i][j][k] = HII[i][j][k] + HeII[i][j][k] / 4.0 + HeIII[i][j][k] / 2.0;
            if (ispecies > 1)
                de[i][j][k] = de[i][j][k] - HM[i][j][k] + H2II[i][j][k] / 2.0;
            dedot_prev[i] = std::abs(de[i][j][k] - dedot_prev[i]) / std::max(dtit[i], tiny8);

            if (ispecies > 1) {
                HM[i][j][k] = std::max(HMp[i], tiny);
                H2I[i][j][k] = std::max(H2Ip[i], tiny);
                H2II[i][j][k] = std::max(H2IIp[i], tiny);
            }

            if (ispecies > 2) {
                DI[i][j][k] = std::max(DIp[i], tiny);
                DII[i][j][k] = std::max(DIIp[i], tiny);
                HDI[i][j][k] = std::max(HDIp[i], tiny);
            }
        }
    }
}

#include <vector>
#include <cmath>
#include <algorithm>

// Function to correct species abundances to ensure conservation
void make_consistent_g(
    // Density fields (inputs and outputs)
    std::vector<std::vector<std::vector<double>>>& de,
    std::vector<std::vector<std::vector<double>>>& HI,
    std::vector<std::vector<std::vector<double>>>& HII,
    std::vector<std::vector<std::vector<double>>>& HeI,
    std::vector<std::vector<std::vector<double>>>& HeII,
    std::vector<std::vector<std::vector<double>>>& HeIII,
    std::vector<std::vector<std::vector<double>>>& HM,
    std::vector<std::vector<std::vector<double>>>& H2I,
    std::vector<std::vector<std::vector<double>>>& H2II,
    std::vector<std::vector<std::vector<double>>>& DI,
    std::vector<std::vector<std::vector<double>>>& DII,
    std::vector<std::vector<std::vector<double>>>& HDI,
    std::vector<std::vector<std::vector<double>>>& metal,
    const std::vector<std::vector<std::vector<double>>>& d,
    // Indices
    int is, int ie, int js, int je, int ks, int ke,
    // Grid dimensions
    int in, int jn, int kn,
    // Species and metal flags
    int ispecies, int imetal,
    // Fractions
    double fh, double dtoh)
{
    // Local variables
    double totalD;
    std::vector<double> totalH(in, 0.0);
    std::vector<double> totalHe(in, 0.0);
    std::vector<double> metalfree(in, 0.0);

    // Loop over all zones
    for (int k = ks + 1; k <= ke + 1; ++k) {
        for (int j = js + 1; j <= je + 1; ++j) {

            // Compute total densities of H and He
            // (ensure non-negativity)

            // Compute metal-free density
            if (imetal == 1) {
                for (int i = is + 1; i <= ie + 1; ++i) {
                    metalfree[i] = d[i][j][k] - metal[i][j][k];
                }
            } else {
                for (int i = is + 1; i <= ie + 1; ++i) {
                    metalfree[i] = d[i][j][k];
                }
            }

            // Ensure species densities are non-negative and compute total H and He densities
            for (int i = is + 1; i <= ie + 1; ++i) {
                HI[i][j][k] = std::abs(HI[i][j][k]);
                HII[i][j][k] = std::abs(HII[i][j][k]);
                HeI[i][j][k] = std::abs(HeI[i][j][k]);
                HeII[i][j][k] = std::abs(HeII[i][j][k]);
                HeIII[i][j][k] = std::abs(HeIII[i][j][k]);

                totalH[i] = HI[i][j][k] + HII[i][j][k];
                totalHe[i] = HeI[i][j][k] + HeII[i][j][k] + HeIII[i][j][k];
            }

            // Include molecular hydrogen if ispecies > 1
            if (ispecies > 1) {
                for (int i = is + 1; i <= ie + 1; ++i) {
                    HM[i][j][k] = std::abs(HM[i][j][k]);
                    H2II[i][j][k] = std::abs(H2II[i][j][k]);
                    H2I[i][j][k] = std::abs(H2I[i][j][k]);

                    totalH[i] += HM[i][j][k] + H2I[i][j][k] + H2II[i][j][k];
                }
            }

            // Correct densities by keeping fractions the same
            for (int i = is + 1; i <= ie + 1; ++i) {
                double correctH = fh * metalfree[i] / totalH[i];
                HI[i][j][k] *= correctH;
                HII[i][j][k] *= correctH;

                double correctHe = (1.0 - fh) * metalfree[i] / totalHe[i];
                HeI[i][j][k] *= correctHe;
                HeII[i][j][k] *= correctHe;
                HeIII[i][j][k] *= correctHe;

                // Correct molecular hydrogen-related fractions
                if (ispecies > 1) {
                    HM[i][j][k] *= correctH;
                    H2II[i][j][k] *= correctH;
                    H2I[i][j][k] *= correctH;
                }
            }

            // Do the same thing for deuterium (ignore HD) Assumes dtoh is small
            if (ispecies > 2) {
                for (int i = is + 1; i <= ie + 1; ++i) {
                    DI[i][j][k] = std::abs(DI[i][j][k]);
                    DII[i][j][k] = std::abs(DII[i][j][k]);
                    HDI[i][j][k] = std::abs(HDI[i][j][k]);

                    totalD = DI[i][j][k] + DII[i][j][k] + (2.0 / 3.0) * HDI[i][j][k];

                    double correctD = fh * dtoh * metalfree[i] / totalD;
                    DI[i][j][k] *= correctD;
                    DII[i][j][k] *= correctD;
                    HDI[i][j][k] *= correctD;
                }
            }

            // Set the electron density
            for (int i = is + 1; i <= ie + 1; ++i) {
                de[i][j][k] = HII[i][j][k] + HeII[i][j][k] / 4.0 + HeIII[i][j][k] / 2.0;
                if (ispecies > 1) {
                    de[i][j][k] = de[i][j][k] - HM[i][j][k] + H2II[i][j][k] / 2.0;
                }
            }

        } // end loop over j
    } // end loop over k
}

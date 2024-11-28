// Necessary includes
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdint> // For 64-bit integers
#ifdef _OPENMP
#include <omp.h>
#endif

// Function declaration
void cool_multi_time_g(
    // Arrays
    std::vector<std::vector<std::vector<double>>>& d,
    std::vector<std::vector<std::vector<double>>>& e,
    std::vector<std::vector<std::vector<double>>>& u,
    std::vector<std::vector<std::vector<double>>>& v,
    std::vector<std::vector<std::vector<double>>>& w,
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
    std::vector<std::vector<std::vector<double>>>& dust,
    std::vector<std::vector<std::vector<double>>>& Vheat,
    std::vector<std::vector<std::vector<double>>>& Mheat,
    std::vector<std::vector<std::vector<double>>>& Tfloor,
    std::vector<std::vector<std::vector<double>>>& isrf_habing,
    std::vector<std::vector<std::vector<double>>>& cooltime,
    std::vector<std::vector<std::vector<double>>>& photogamma,
    // Scalars
    int in, int jn, int kn, int nratec, int iexpand,
    int ispecies, int imetal, int imcool, int idust,
    int idustall, int idustfield, int idustrec, int idim,
    int is, int js, int ks, int ie, int je, int ke, int ih2co, int ipiht, int igammah,
    double aye, double temstart, double temend,
    double utem, double uxyz, double uaye, double urho, double utim,
    double gamma_val, double fh, double z_solar, double fgr, double clEleFra, double Tfloor_scalar,
    // Cooling arrays
    const std::vector<double>& ceHIa, const std::vector<double>& ceHeIa, const std::vector<double>& ceHeIIa,
    const std::vector<double>& ciHIa, const std::vector<double>& ciHeIa, const std::vector<double>& ciHeISa,
    const std::vector<double>& ciHeIIa, const std::vector<double>& reHIIa, const std::vector<double>& reHeII1a,
    const std::vector<double>& reHeII2a, const std::vector<double>& reHeIIIa, const std::vector<double>& brema,
    double compa, double piHI, double piHeI, double piHeII, double comp_xraya, double comp_temp,
    double gammaha, double isrf, const std::vector<double>& regra, const std::vector<double>& gamma_isrfa,
    double avgsighi, double avgsighei, double avgsigheii,
    double k24, double k26,
    // Cloudy cooling data
    int icmbTfloor, int iClHeat,
    int64_t priGridRank, const std::vector<int64_t>& priGridDim,
    const std::vector<double>& priPar1, const std::vector<double>& priPar2, const std::vector<double>& priPar3,
    const std::vector<double>& priPar4, const std::vector<double>& priPar5,
    int64_t priDataSize, const std::vector<double>& priCooling, const std::vector<double>& priHeating,
    const std::vector<double>& priMMW,
    int64_t metGridRank, const std::vector<int64_t>& metGridDim,
    const std::vector<double>& metPar1, const std::vector<double>& metPar2, const std::vector<double>& metPar3,
    const std::vector<double>& metPar4, const std::vector<double>& metPar5,
    int64_t metDataSize, const std::vector<double>& metCooling, const std::vector<double>& metHeating, int clnew,
    // Additional parameters
    int iVheat, int iMheat, int iTfloor, int iradtrans, int iradshield, int iisrffield,
    // H2 cooling arrays
    const std::vector<double>& hyd01ka, const std::vector<double>& h2k01a, const std::vector<double>& vibha,
    const std::vector<double>& rotha, const std::vector<double>& rotla,
    const std::vector<double>& gpldla, const std::vector<double>& gphdla,
    const std::vector<double>& hdltea, const std::vector<double>& hdlowa,
    const std::vector<double>& gaHIa, const std::vector<double>& gaH2a, const std::vector<double>& gaHea,
    const std::vector<double>& gaHpa, const std::vector<double>& gaela,
    const std::vector<double>& h2ltea, const std::vector<double>& gasgra,
    int ih2optical, int iciecool, const std::vector<double>& ciecoa
)
{
    // Local variables
    int i, j, k;
    int t, dj, dk;
    double comp1, comp2, energy;

    // Slice locals
    std::vector<int64_t> indixe(in);
    std::vector<double> t1(in), t2(in), logtem(in), tdef(in), p2d(in);
    std::vector<double> tgas(in), tgasold(in), mmw(in);
    std::vector<double> tdust(in), metallicity(in), dust2gas(in), rhoH(in);
    std::vector<double> mynh(in), myde(in), gammaha_eff(in);
    std::vector<double> gasgr_tdust(in), regr(in), edot(in);

    // Cooling/heating slice locals
    std::vector<double> ceHI(in), ceHeI(in), ceHeII(in);
    std::vector<double> ciHI(in), ciHeI(in), ciHeIS(in), ciHeII(in);
    std::vector<double> reHII(in), reHeII1(in), reHeII2(in), reHeIII(in);
    std::vector<double> brem(in), cieco(in);
    std::vector<double> hyd01k(in), h2k01(in), vibh(in), roth(in), rotl(in);
    std::vector<double> gpldl(in), gphdl(in), hdlte(in), hdlow(in);

    // Iteration mask
    std::vector<bool> itmask(in, true);

    // Convert densities from comoving to proper units if iexpand == 1
    if (iexpand == 1) {
        dk = ke - ks + 1;
        dj = je - js + 1;

        #pragma omp parallel for schedule(runtime) private(i, j, k) if (_OPENMP)
        for (t = 0; t < dk * dj; ++t) {
            k = t / dj + ks;
            j = t % dj + js;

            double aye_cubed = pow(aye, 3);

            for (i = is; i <= ie; ++i) {
                d[i][j][k] /= aye_cubed;
            }

            if (ispecies > 0) {
                for (i = is; i <= ie; ++i) {
                    de[i][j][k] /= aye_cubed;
                    HI[i][j][k] /= aye_cubed;
                    HII[i][j][k] /= aye_cubed;
                    HeI[i][j][k] /= aye_cubed;
                    HeII[i][j][k] /= aye_cubed;
                    HeIII[i][j][k] /= aye_cubed;
                }
            }

            if (ispecies > 1) {
                for (i = is; i <= ie; ++i) {
                    HM[i][j][k] /= aye_cubed;
                    H2I[i][j][k] /= aye_cubed;
                    H2II[i][j][k] /= aye_cubed;
                }
            }

            if (ispecies > 2) {
                for (i = is; i <= ie; ++i) {
                    DI[i][j][k] /= aye_cubed;
                    DII[i][j][k] /= aye_cubed;
                    HDI[i][j][k] /= aye_cubed;
                }
            }

            if (imetal == 1) {
                for (i = is; i <= ie; ++i) {
                    metal[i][j][k] /= aye_cubed;
                }
            }

            if (idustfield == 1) {
                for (i = is; i <= ie; ++i) {
                    dust[i][j][k] /= aye_cubed;
                }
            }
        }
    }

    // Loop over slices in the k-direction
    dk = ke - ks + 1;
    dj = je - js + 1;

    #pragma omp parallel for schedule(runtime) private(i, j, k, t, energy) if (_OPENMP)
    for (t = 0; t < dk * dj; ++t) {
        k = t / dj + ks;
        j = t % dj + js;

        // Reset the iteration mask
        std::fill(itmask.begin(), itmask.end(), true);

        // Compute the cooling rate using cool1d_multi_g (function must be implemented)
        cool1d_multi_g(
            d, e, u, v, w, de, HI, HII, HeI, HeII, HeIII,
            in, jn, kn, nratec,
            iexpand, ispecies, imetal, imcool,
            idust, idustall, idustfield, idustrec,
            idim, is, ie, j, k, ih2co, ipiht, 1, igammah,
            aye, temstart, temend, z_solar, fgr,
            utem, uxyz, uaye, urho, utim,
            gamma_val, fh,
            ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa,
            ciHeISa, ciHeIIa, reHIIa, reHeII1a,
            reHeII2a, reHeIIIa, brema, compa, gammaha,
            isrf, regra, gamma_isrfa, comp_xraya, comp_temp,
            piHI, piHeI, piHeII, comp1, comp2,
            HM, H2I, H2II, DI, DII, HDI, metal, dust,
            hyd01ka, h2k01a, vibha, rotha, rotla,
            hyd01k, h2k01, vibh, roth, rotl,
            gpldla, gphdla, gpldl, gphdl,
            hdltea, hdlowa, hdlte, hdlow,
            gaHIa, gaH2a, gaHea, gaHpa, gaela,
            h2ltea, gasgra,
            ciecoa, cieco,
            priGridRank, priGridDim,
            priPar1, priPar2, priPar3, priPar4, priPar5,
            priDataSize, priCooling, priHeating, priMMW,
            metGridRank, metGridDim,
            metPar1, metPar2, metPar3, metPar4, metPar5,
            metDataSize, metCooling, metHeating, clnew,
            icmbTfloor, iClHeat, clEleFra,
            iVheat, iMheat, Vheat, Mheat,
            iTfloor, Tfloor_scalar, Tfloor,
            iisrffield, isrf_habing,
            iradshield, avgsighi, avgsighei, avgsigheii,
            k24, k26, iradtrans, photogamma,
            ih2optical, iciecool,
            is, ie, j, k, 1,
            tgas, tgasold, edot,
            itmask
        );

        // Compute the cooling time on the slice
        for (i = is; i <= ie; ++i) {
            energy = std::max(p2d[i], tiny); // Ensure energy is not zero or negative
            cooltime[i][j][k] = energy / edot[i];
        }
    }

    // Convert densities back to comoving units if iexpand == 1
    if (iexpand == 1) {
        #pragma omp parallel for schedule(runtime) private(i, j, k) if (_OPENMP)
        for (t = 0; t < dk * dj; ++t) {
            k = t / dj + ks;
            j = t % dj + js;

            double aye_cubed = pow(aye, 3);

            for (i = is; i <= ie; ++i) {
                d[i][j][k] *= aye_cubed;
            }

            if (ispecies > 0) {
                for (i = is; i <= ie; ++i) {
                    de[i][j][k] *= aye_cubed;
                    HI[i][j][k] *= aye_cubed;
                    HII[i][j][k] *= aye_cubed;
                    HeI[i][j][k] *= aye_cubed;
                    HeII[i][j][k] *= aye_cubed;
                    HeIII[i][j][k] *= aye_cubed;
                }
            }

            if (ispecies > 1) {
                for (i = is; i <= ie; ++i) {
                    HM[i][j][k] *= aye_cubed;
                    H2I[i][j][k] *= aye_cubed;
                    H2II[i][j][k] *= aye_cubed;
                }
            }

            if (ispecies > 2) {
                for (i = is; i <= ie; ++i) {
                    DI[i][j][k] *= aye_cubed;
                    DII[i][j][k] *= aye_cubed;
                    HDI[i][j][k] *= aye_cubed;
                }
            }

            if (imetal == 1) {
                for (i = is; i <= ie; ++i) {
                    metal[i][j][k] *= aye_cubed;
                }
            }

            if (idustfield == 1) {
                for (i = is; i <= ie; ++i) {
                    dust[i][j][k] *= aye_cubed;
                }
            }
        }
    }
}

//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares signature of make_consistent_g
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// make_consistent_g function from FORTRAN to C++

#include <cstdio>
#include <vector>

#include "grackle.h"
#include "fortran_func_decls.h"
#include "utils-cpp.hpp"

#include "make_consistent.hpp"

namespace grackle::impl {


void make_consistent(
  const int* imetal, const double* dom, chemistry_data* my_chemistry,
  chemistry_data_storage* my_rates, grackle_field_data* my_fields
)
{
  // -------------------------------------------------------------------


  // Arguments

  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  grackle::impl::View<gr_float***> de(my_fields->e_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HI(my_fields->HI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HII(my_fields->HII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeI(my_fields->HeI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeII(my_fields->HeII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeIII(my_fields->HeIII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> d(my_fields->density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal(my_fields->metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HM(my_fields->HM_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2I(my_fields->H2I_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2II(my_fields->H2II_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> DI(my_fields->DI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> DII(my_fields->DII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HDI(my_fields->HDI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  grackle::impl::View<gr_float***> DM(my_fields->DM_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HDII(my_fields->HDII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeHII(my_fields->HeHII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CI(my_fields->CI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CII(my_fields->CII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CO(my_fields->CO_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CO2(my_fields->CO2_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OI(my_fields->OI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OH(my_fields->OH_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2O(my_fields->H2O_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> O2(my_fields->O2_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiI(my_fields->SiI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiOI(my_fields->SiOI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiO2I(my_fields->SiO2I_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CH(my_fields->CH_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CH2(my_fields->CH2_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> COII(my_fields->COII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OII(my_fields->OII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OHII(my_fields->OHII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2OII(my_fields->H2OII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H3OII(my_fields->H3OII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> O2II(my_fields->O2II_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Mg(my_fields->Mg_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Al(my_fields->Al_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> S(my_fields->S_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Fe(my_fields->Fe_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiM(my_fields->SiM_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> FeM(my_fields->FeM_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Mg2SiO4(my_fields->Mg2SiO4_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> MgSiO3(my_fields->MgSiO3_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Fe3O4(my_fields->Fe3O4_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> AC(my_fields->AC_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiO2D(my_fields->SiO2_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> MgO(my_fields->MgO_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> FeS(my_fields->FeS_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Al2O3(my_fields->Al2O3_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> reforg(my_fields->ref_org_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> volorg(my_fields->vol_org_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2Oice(my_fields->H2O_ice_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  // -- removed line (previously just declared arg types) -- 
  grackle::impl::View<gr_float***> metal_loc(my_fields->local_ISM_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_C13(my_fields->ccsn13_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_C20(my_fields->ccsn20_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_C25(my_fields->ccsn25_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_C30(my_fields->ccsn30_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_F13(my_fields->fsn13_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_F15(my_fields->fsn15_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_F50(my_fields->fsn50_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_F80(my_fields->fsn80_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_P170(my_fields->pisn170_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_P200(my_fields->pisn200_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_Y19(my_fields->y19_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 

  // locals

  int i, j, k;
  double totalD;
  std::vector<double> totalH(my_fields->grid_dimension[0]);
  std::vector<double> totalHe(my_fields->grid_dimension[0]);
  std::vector<double> metalfree(my_fields->grid_dimension[0]);
  gr_float correctH, correctHe, correctD;
  double totalZ;
  double totalC, totalO, totalMg, totalAl, totalSi, totalS, totalFe;
  double totalCg, totalOg, totalMgg, totalAlg, totalSig, totalSg, totalFeg;
  double totalCd, totalOd, totalMgd, totalAld, totalSid, totalSd, totalFed;
  gr_float correctC, correctO, correctMg, correctAl, correctSi, correctS, correctFe;
  gr_float correctCg, correctOg, correctMgg, correctAlg, correctSig, correctSg, correctFeg;
  gr_float correctCd, correctOd, correctMgd, correctAld, correctSid, correctSd, correctFed;
  gr_float correctZ;
  int iSN, nSN, iSN0;
  std::vector<int> SN_i(my_rates->SN0_N);
  std::vector<gr_float> SN_metal_data_(my_fields->grid_dimension[0] * my_rates->SN0_N);
  grackle::impl::View<gr_float**> SN_metal(SN_metal_data_.data(), my_fields->grid_dimension[0], my_rates->SN0_N);
  std::vector<double> Ct(my_fields->grid_dimension[0]);
  std::vector<double> Ot(my_fields->grid_dimension[0]);
  std::vector<double> Mgt(my_fields->grid_dimension[0]);
  std::vector<double> Alt(my_fields->grid_dimension[0]);
  std::vector<double> Sit(my_fields->grid_dimension[0]);
  std::vector<double> St(my_fields->grid_dimension[0]);
  std::vector<double> Fet(my_fields->grid_dimension[0]);
  std::vector<double> Cg(my_fields->grid_dimension[0]);
  std::vector<double> Og(my_fields->grid_dimension[0]);
  std::vector<double> Mgg(my_fields->grid_dimension[0]);
  std::vector<double> Alg(my_fields->grid_dimension[0]);
  std::vector<double> Sig(my_fields->grid_dimension[0]);
  std::vector<double> Sg(my_fields->grid_dimension[0]);
  std::vector<double> Feg(my_fields->grid_dimension[0]);
  std::vector<double> Cd(my_fields->grid_dimension[0]);
  std::vector<double> Od(my_fields->grid_dimension[0]);
  std::vector<double> Mgd(my_fields->grid_dimension[0]);
  std::vector<double> Ald(my_fields->grid_dimension[0]);
  std::vector<double> Sid(my_fields->grid_dimension[0]);
  std::vector<double> Sd(my_fields->grid_dimension[0]);
  std::vector<double> Fed(my_fields->grid_dimension[0]);

  // Loop over all zones

  for (k = my_fields->grid_start[2] + 1; k<=(my_fields->grid_end[2] + 1); k++) {
    for (j = my_fields->grid_start[1] + 1; j<=(my_fields->grid_end[1] + 1); j++) {

      // Compute total densities of H and He
      //     (ensure non-negativity)

      if ((*imetal) == 1)  {
        for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
          metalfree[i-1] = d(i-1,j-1,k-1) - metal(i-1,j-1,k-1);
        }
      } else {
        for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
          metalfree[i-1] = d(i-1,j-1,k-1);
        }
      }

      for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
        HI   (i-1,j-1,k-1) = std::fabs(HI   (i-1,j-1,k-1));
        HII  (i-1,j-1,k-1) = std::fabs(HII  (i-1,j-1,k-1));
        HeI  (i-1,j-1,k-1) = std::fabs(HeI  (i-1,j-1,k-1));
        HeII (i-1,j-1,k-1) = std::fabs(HeII (i-1,j-1,k-1));
        HeIII(i-1,j-1,k-1) = std::fabs(HeIII(i-1,j-1,k-1));
        totalH[i-1] = HI(i-1,j-1,k-1) + HII(i-1,j-1,k-1);
        totalHe[i-1] = HeI(i-1,j-1,k-1) + HeII(i-1,j-1,k-1) + HeIII(i-1,j-1,k-1);
      }

      // include molecular hydrogen

      if (my_chemistry->primordial_chemistry > 1)  {
        for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
          HM   (i-1,j-1,k-1) = std::fabs(HM   (i-1,j-1,k-1));
          H2II (i-1,j-1,k-1) = std::fabs(H2II (i-1,j-1,k-1));
          H2I  (i-1,j-1,k-1) = std::fabs(H2I  (i-1,j-1,k-1));
          totalH[i-1] = totalH[i-1] + HM(i-1,j-1,k-1) + H2I(i-1,j-1,k-1) + H2II(i-1,j-1,k-1);
        }
      }

      if(my_chemistry->primordial_chemistry > 2)  {
        for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
          HDI(i-1,j-1,k-1) = std::fabs(HDI(i-1,j-1,k-1));
          totalH [i-1] = totalH [i-1]
               + 1./3.*HDI(i-1,j-1,k-1);
        }
      }
      // ! GC202005

      if(my_chemistry->primordial_chemistry > 3)  {
        for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
          HDII (i-1,j-1,k-1) = std::fabs(HDII (i-1,j-1,k-1));
          HeHII(i-1,j-1,k-1) = std::fabs(HeHII(i-1,j-1,k-1));
          totalH [i-1] = totalH [i-1]
               + 1./3.*HDII (i-1,j-1,k-1)
               + 1./5.*HeHII(i-1,j-1,k-1);
          totalHe[i-1] = totalHe[i-1]
               + 4./5.*HeHII(i-1,j-1,k-1);
        }
      }

      // Iteration mask for metal-rich cells

      // do i = is+1, ie + 1
      //    itmask_metal(i) = .false.
      // enddo
      // if (imetal .eq. 1) then
      //     do i = is+1, ie + 1
      //        if (metal(i,j,k) .gt. 1.e-9_DKIND * d(i,j,k)) then
      //           itmask_metal(i) = .true.
      //        endif
      //     enddo
      // endif

      if(my_chemistry->metal_chemistry > 0)  {
        if(my_chemistry->multi_metals == 0)  {
          iSN0 = my_chemistry->metal_abundances + 1;
          for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
            Ct[i-1] = my_rates->SN0_XC [iSN0-1] * metal(i-1,j-1,k-1);
            Ot[i-1] = my_rates->SN0_XO [iSN0-1] * metal(i-1,j-1,k-1);
            Mgt[i-1] = my_rates->SN0_XMg[iSN0-1] * metal(i-1,j-1,k-1);
            Alt[i-1] = my_rates->SN0_XAl[iSN0-1] * metal(i-1,j-1,k-1);
            Sit[i-1] = my_rates->SN0_XSi[iSN0-1] * metal(i-1,j-1,k-1);
            St[i-1] = my_rates->SN0_XS [iSN0-1] * metal(i-1,j-1,k-1);
            Fet[i-1] = my_rates->SN0_XFe[iSN0-1] * metal(i-1,j-1,k-1);
         
            Cg[i-1] = my_rates->SN0_fC [iSN0-1] * metal(i-1,j-1,k-1);
            Og[i-1] = my_rates->SN0_fO [iSN0-1] * metal(i-1,j-1,k-1);
            Mgg[i-1] = my_rates->SN0_fMg[iSN0-1] * metal(i-1,j-1,k-1);
            Alg[i-1] = my_rates->SN0_fAl[iSN0-1] * metal(i-1,j-1,k-1);
            Sig[i-1] = my_rates->SN0_fSi[iSN0-1] * metal(i-1,j-1,k-1);
            Sg[i-1] = my_rates->SN0_fS [iSN0-1] * metal(i-1,j-1,k-1);
            Feg[i-1] = my_rates->SN0_fFe[iSN0-1] * metal(i-1,j-1,k-1);
          }
         
        } else {

          //        do i = is+1, ie+1
          //           totalZ = metal_loc(i,j,k)
          // &           + metal_C13(i,j,k) + metal_C20(i,j,k)
          // &           + metal_C25(i,j,k) + metal_C30(i,j,k)
          // &           + metal_F13(i,j,k) + metal_F15(i,j,k)
          // &           + metal_F50(i,j,k) + metal_F80(i,j,k)
          // &           + metal_P170(i,j,k)+ metal_P200(i,j,k)
          // &           + metal_Y19(i,j,k)
          //           correctZ = metal(i,j,k) / totalZ
          //           metal_loc(i,j,k) = metal_loc(i,j,k) * correctZ
          //           metal_C13(i,j,k) = metal_C13(i,j,k) * correctZ
          //           metal_C20(i,j,k) = metal_C20(i,j,k) * correctZ
          //           metal_C25(i,j,k) = metal_C25(i,j,k) * correctZ
          //           metal_C30(i,j,k) = metal_C30(i,j,k) * correctZ
          //           metal_F13(i,j,k) = metal_F13(i,j,k) * correctZ
          //           metal_F15(i,j,k) = metal_F15(i,j,k) * correctZ
          //           metal_F50(i,j,k) = metal_F50(i,j,k) * correctZ
          //           metal_F80(i,j,k) = metal_F80(i,j,k) * correctZ
          //           metal_P170(i,j,k)= metal_P170(i,j,k)* correctZ
          //           metal_P200(i,j,k)= metal_P200(i,j,k)* correctZ
          //           metal_Y19(i,j,k) = metal_Y19(i,j,k) * correctZ
          //        enddo

          nSN = 12;
          SN_i[ 1-1] = 1;
  //_// PORT:             SN_metal(:, 1) = metal_loc(:,j,k)
          SN_i[ 2-1] = 2;
  //_// PORT:             SN_metal(:, 2) = metal_C13(:,j,k)
          SN_i[ 3-1] = 3;
  //_// PORT:             SN_metal(:, 3) = metal_C20(:,j,k)
          SN_i[ 4-1] = 4;
  //_// PORT:             SN_metal(:, 4) = metal_C25(:,j,k)
          SN_i[ 5-1] = 5;
  //_// PORT:             SN_metal(:, 5) = metal_C30(:,j,k)
          SN_i[ 6-1] = 6;
  //_// PORT:             SN_metal(:, 6) = metal_F13(:,j,k)
          SN_i[ 7-1] = 7;
  //_// PORT:             SN_metal(:, 7) = metal_F15(:,j,k)
          SN_i[ 8-1] = 8;
  //_// PORT:             SN_metal(:, 8) = metal_F50(:,j,k)
          SN_i[ 9-1] = 9;
  //_// PORT:             SN_metal(:, 9) = metal_F80(:,j,k)
          SN_i[10-1] =10;
  //_// PORT:             SN_metal(:,10) = metal_P170(:,j,k)
          SN_i[11-1] =11;
  //_// PORT:             SN_metal(:,11) = metal_P200(:,j,k)
          SN_i[12-1] =12;
  //_// PORT:             SN_metal(:,12) = metal_Y19(:,j,k)
         
          for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
            Ct[i-1] = 0.;
            Cg[i-1] = 0.;
            Ot[i-1] = 0.;
            Og[i-1] = 0.;
            Mgt[i-1] = 0.;
            Mgg[i-1] = 0.;
            Alt[i-1] = 0.;
            Alg[i-1] = 0.;
            Sit[i-1] = 0.;
            Sig[i-1] = 0.;
            St[i-1] = 0.;
            Sg[i-1] = 0.;
            Fet[i-1] = 0.;
            Feg[i-1] = 0.;
            for (iSN = 1; iSN<=(nSN); iSN++) {
              iSN0 = SN_i[iSN-1];

              Ct[i-1] =  Ct[i-1] + my_rates->SN0_XC [iSN0-1] * SN_metal(i-1,iSN-1);
              Ot[i-1] =  Ot[i-1] + my_rates->SN0_XO [iSN0-1] * SN_metal(i-1,iSN-1);
              Mgt[i-1] = Mgt[i-1] + my_rates->SN0_XMg[iSN0-1] * SN_metal(i-1,iSN-1);
              Alt[i-1] = Alt[i-1] + my_rates->SN0_XAl[iSN0-1] * SN_metal(i-1,iSN-1);
              Sit[i-1] = Sit[i-1] + my_rates->SN0_XSi[iSN0-1] * SN_metal(i-1,iSN-1);
              St[i-1] =  St[i-1] + my_rates->SN0_XS [iSN0-1] * SN_metal(i-1,iSN-1);
              Fet[i-1] = Fet[i-1] + my_rates->SN0_XFe[iSN0-1] * SN_metal(i-1,iSN-1);

              Cg[i-1] =  Cg[i-1] + my_rates->SN0_fC [iSN0-1] * SN_metal(i-1,iSN-1);
              Og[i-1] =  Og[i-1] + my_rates->SN0_fO [iSN0-1] * SN_metal(i-1,iSN-1);
              Mgg[i-1] = Mgg[i-1] + my_rates->SN0_fMg[iSN0-1] * SN_metal(i-1,iSN-1);
              Alg[i-1] = Alg[i-1] + my_rates->SN0_fAl[iSN0-1] * SN_metal(i-1,iSN-1);
              Sig[i-1] = Sig[i-1] + my_rates->SN0_fSi[iSN0-1] * SN_metal(i-1,iSN-1);
              Sg[i-1] =  Sg[i-1] + my_rates->SN0_fS [iSN0-1] * SN_metal(i-1,iSN-1);
              Feg[i-1] = Feg[i-1] + my_rates->SN0_fFe[iSN0-1] * SN_metal(i-1,iSN-1);
            }
          }
            
        }
            
        for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
          Cd[i-1] =  Ct[i-1] -  Cg[i-1];
          Od[i-1] =  Ot[i-1] -  Og[i-1];
          Mgd[i-1] = Mgt[i-1] - Mgg[i-1];
          Ald[i-1] = Alt[i-1] - Alg[i-1];
          Sid[i-1] = Sit[i-1] - Sig[i-1];
          Sd[i-1] =  St[i-1] -  Sg[i-1];
          Fed[i-1] = Fet[i-1] - Feg[i-1];
        }

        for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
          // if (itmask_metal(i)) then
          OH   (i-1,j-1,k-1) = std::fabs(OH   (i-1,j-1,k-1));
          H2O  (i-1,j-1,k-1) = std::fabs(H2O  (i-1,j-1,k-1));
          CH   (i-1,j-1,k-1) = std::fabs(CH   (i-1,j-1,k-1));
          CH2  (i-1,j-1,k-1) = std::fabs(CH2  (i-1,j-1,k-1));
          OHII (i-1,j-1,k-1) = std::fabs(OHII (i-1,j-1,k-1));
          H2OII(i-1,j-1,k-1) = std::fabs(H2OII(i-1,j-1,k-1));
          H3OII(i-1,j-1,k-1) = std::fabs(H3OII(i-1,j-1,k-1));
          totalH[i-1] = totalH[i-1]
            + OH   (i-1,j-1,k-1)/17.
            + H2O  (i-1,j-1,k-1)/18.*2.
            + CH   (i-1,j-1,k-1)/13.
            + CH2  (i-1,j-1,k-1)/14.*2.
            + OHII (i-1,j-1,k-1)/17.
            + H2OII(i-1,j-1,k-1)/18.*2.
            + H3OII(i-1,j-1,k-1)/19.*3.;
          // endif
        }
      }

      // Correct densities by keeping fractions the same

      for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
        correctH = (gr_float)(my_chemistry->HydrogenFractionByMass*metalfree[i-1]/totalH[i-1] ); // TODO: which one is correct?
        // correctH = (gr_float)(my_chemistry->HydrogenFractionByMass*(1. - my_chemistry->DeuteriumToHydrogenRatio)*metalfree[i-1]/totalH[i-1]
        //              );
        // // ! GC202005
        // //- !       correctH = real(fh*metalfree(i)/totalH(i), RKIND)
        HI(i-1,j-1,k-1)  = HI(i-1,j-1,k-1)*correctH;
        HII(i-1,j-1,k-1) = HII(i-1,j-1,k-1)*correctH;

        correctHe = (gr_float)((1. - my_chemistry->HydrogenFractionByMass)*
             metalfree[i-1]/totalHe[i-1] );
        HeI(i-1,j-1,k-1)   = HeI(i-1,j-1,k-1)*correctHe;
        HeII(i-1,j-1,k-1)  = HeII(i-1,j-1,k-1)*correctHe;
        HeIII(i-1,j-1,k-1) = HeIII(i-1,j-1,k-1)*correctHe;

        // Correct molecular hydrogen-related fractions

        if (my_chemistry->primordial_chemistry > 1)  {
          HM   (i-1,j-1,k-1) = HM(i-1,j-1,k-1)*correctH;
          H2II (i-1,j-1,k-1) = H2II(i-1,j-1,k-1)*correctH;
          H2I  (i-1,j-1,k-1) = H2I(i-1,j-1,k-1)*correctH;
        }
        if(my_chemistry->primordial_chemistry > 3)  {
          // !          HDII (i,j,k) = HDII (i,j,k)*correctH
          HeHII(i-1,j-1,k-1) = HeHII(i-1,j-1,k-1)*correctHe;
        }
      }

      // Do the same thing for deuterium (ignore HD) Assumes dtoh is small

      if (my_chemistry->primordial_chemistry > 2)  {
        for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
          DI  (i-1,j-1,k-1) = std::fabs(DI  (i-1,j-1,k-1));
          DII (i-1,j-1,k-1) = std::fabs(DII (i-1,j-1,k-1));
          HDI (i-1,j-1,k-1) = std::fabs(HDI (i-1,j-1,k-1));
          totalD = DI(i-1,j-1,k-1) + DII(i-1,j-1,k-1) +
               2./3.*HDI(i-1,j-1,k-1);
          if(my_chemistry->primordial_chemistry > 3)  {
            DM   (i-1,j-1,k-1) = std::fabs(DM   (i-1,j-1,k-1));
            HDII (i-1,j-1,k-1) = std::fabs(HDII (i-1,j-1,k-1));
            totalD = totalD + DM(i-1,j-1,k-1) +
              2./3.*HDII(i-1,j-1,k-1);
          }
          correctD = (gr_float)(my_chemistry->HydrogenFractionByMass*my_chemistry->DeuteriumToHydrogenRatio*metalfree[i-1]/totalD );
          DI  (i-1,j-1,k-1) = DI (i-1,j-1,k-1)*correctD;
          DII (i-1,j-1,k-1) = DII(i-1,j-1,k-1)*correctD;
          HDI (i-1,j-1,k-1) = HDI(i-1,j-1,k-1)*correctD;
          if(my_chemistry->primordial_chemistry > 3)  {
            DM   (i-1,j-1,k-1) = DM   (i-1,j-1,k-1)*correctD;
            HDII (i-1,j-1,k-1) = HDII (i-1,j-1,k-1)*correctD;
          }
        }
      }

      // Do the same thing for metal species

      if (my_chemistry->metal_chemistry == 1)  {
        for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
          // if (itmask_metal(i)) then
          CI(i-1,j-1,k-1)      = std::fabs(CI(i-1,j-1,k-1)     );
          CII(i-1,j-1,k-1)     = std::fabs(CII(i-1,j-1,k-1)    );
          CO(i-1,j-1,k-1)      = std::fabs(CO(i-1,j-1,k-1)     );
          CO2(i-1,j-1,k-1)     = std::fabs(CO2(i-1,j-1,k-1)    );
          OI(i-1,j-1,k-1)      = std::fabs(OI(i-1,j-1,k-1)     );
          OH(i-1,j-1,k-1)      = std::fabs(OH(i-1,j-1,k-1)     );
          H2O(i-1,j-1,k-1)     = std::fabs(H2O(i-1,j-1,k-1)    );
          O2(i-1,j-1,k-1)      = std::fabs(O2(i-1,j-1,k-1)     );
          SiI(i-1,j-1,k-1)     = std::fabs(SiI(i-1,j-1,k-1)    );
          SiOI(i-1,j-1,k-1)    = std::fabs(SiOI(i-1,j-1,k-1)   );
          SiO2I(i-1,j-1,k-1)   = std::fabs(SiO2I(i-1,j-1,k-1)  );
          CH(i-1,j-1,k-1)      = std::fabs(CH(i-1,j-1,k-1)     );
          CH2(i-1,j-1,k-1)     = std::fabs(CH2(i-1,j-1,k-1)    );
          COII(i-1,j-1,k-1)    = std::fabs(COII(i-1,j-1,k-1)   );
          OII(i-1,j-1,k-1)     = std::fabs(OII(i-1,j-1,k-1)    );
          OHII(i-1,j-1,k-1)    = std::fabs(OHII(i-1,j-1,k-1)   );
          H2OII(i-1,j-1,k-1)   = std::fabs(H2OII(i-1,j-1,k-1)  );
          H3OII(i-1,j-1,k-1)   = std::fabs(H3OII(i-1,j-1,k-1)  );
          O2II(i-1,j-1,k-1)    = std::fabs(O2II(i-1,j-1,k-1)   );
          if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
            if (my_chemistry->dust_species > 0)  {
              Mg(i-1,j-1,k-1)      = std::fabs(Mg(i-1,j-1,k-1)     );
            }
            if (my_chemistry->dust_species > 1)  {
              Al(i-1,j-1,k-1)      = std::fabs(Al(i-1,j-1,k-1)     );
              S(i-1,j-1,k-1)       = std::fabs(S(i-1,j-1,k-1)      );
              Fe(i-1,j-1,k-1)      = std::fabs(Fe(i-1,j-1,k-1)     );
            }
          }
          // endif
        }
      }

      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
          // if (itmask_metal(i)) then
          if (my_chemistry->dust_species > 0)  {
            MgSiO3(i-1,j-1,k-1)  = std::fabs(MgSiO3(i-1,j-1,k-1) );
            AC(i-1,j-1,k-1)      = std::fabs(AC(i-1,j-1,k-1)     );
          }
          if (my_chemistry->dust_species > 1)  {
            SiM(i-1,j-1,k-1)     = std::fabs(SiM(i-1,j-1,k-1)    );
            FeM(i-1,j-1,k-1)     = std::fabs(FeM(i-1,j-1,k-1)    );
            Mg2SiO4(i-1,j-1,k-1) = std::fabs(Mg2SiO4(i-1,j-1,k-1));
            Fe3O4(i-1,j-1,k-1)   = std::fabs(Fe3O4(i-1,j-1,k-1)  );
            SiO2D(i-1,j-1,k-1)   = std::fabs(SiO2D(i-1,j-1,k-1)  );
            MgO(i-1,j-1,k-1)     = std::fabs(MgO(i-1,j-1,k-1)    );
            FeS(i-1,j-1,k-1)     = std::fabs(FeS(i-1,j-1,k-1)    );
            Al2O3(i-1,j-1,k-1)   = std::fabs(Al2O3(i-1,j-1,k-1)  );
          }
          if (my_chemistry->dust_species > 2)  {
            reforg(i-1,j-1,k-1)  = std::fabs(reforg(i-1,j-1,k-1)  );
            volorg(i-1,j-1,k-1)  = std::fabs(volorg(i-1,j-1,k-1)  );
            H2Oice(i-1,j-1,k-1)  = std::fabs(H2Oice(i-1,j-1,k-1)  );
          }
          // endif
        }
      }

      if (my_chemistry->metal_chemistry == 1)  {
        for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
          // if (itmask_metal(i)) then
          //- if (d(i,j,k)*dom .lt. 1.e-2_DKIND) then
          // !       if (d(i,j,k)*dom .lt.
          // !   &    min(1.e6_DKIND/(metal(i,j,k)/d(i,j,k)/0.02d-4)**2
          // !   &       ,1.e6_DKIND)) then
          if ( ( ((*imetal) == 0)
            &&  (d(i-1,j-1,k-1)*(*dom) < 1.e8) )
           ||  ( ((*imetal) == 1)
            &&  ( ( (metal(i-1,j-1,k-1) <= 1.e-9 * d(i-1,j-1,k-1))
                &&  (d(i-1,j-1,k-1)*(*dom) < 1.e8) )
               ||  ( (metal(i-1,j-1,k-1) > 1.e-9 * d(i-1,j-1,k-1))
                &&  (d(i-1,j-1,k-1)*(*dom) < 1.e6) ) ) ) )  {

            totalOg = 16./28.*   CO(i-1,j-1,k-1)
                    + 32./44.*  CO2(i-1,j-1,k-1)
                    +                        OI(i-1,j-1,k-1)
                    + 16./17.*   OH(i-1,j-1,k-1)
                    + 16./18.*  H2O(i-1,j-1,k-1)
                    +                        O2(i-1,j-1,k-1)
                    + 16./44.* SiOI(i-1,j-1,k-1)
                    + 32./60.*SiO2I(i-1,j-1,k-1)
                    + 16./28.* COII(i-1,j-1,k-1)
                    +                       OII(i-1,j-1,k-1)
                    + 16./17.* OHII(i-1,j-1,k-1)
                    + 16./18.*H2OII(i-1,j-1,k-1)
                    + 16./19.*H3OII(i-1,j-1,k-1)
                    +                      O2II(i-1,j-1,k-1);
            correctOg = (gr_float)(Og[i-1]/totalOg );
            CO(i-1,j-1,k-1) =    CO(i-1,j-1,k-1)*correctOg;
            CO2(i-1,j-1,k-1) =   CO2(i-1,j-1,k-1)*correctOg;
            OI(i-1,j-1,k-1) =    OI(i-1,j-1,k-1)*correctOg;
            OH(i-1,j-1,k-1) =    OH(i-1,j-1,k-1)*correctOg;
            H2O(i-1,j-1,k-1) =   H2O(i-1,j-1,k-1)*correctOg;
            O2(i-1,j-1,k-1) =    O2(i-1,j-1,k-1)*correctOg;
            SiOI(i-1,j-1,k-1) =  SiOI(i-1,j-1,k-1)*correctOg;
            SiO2I(i-1,j-1,k-1) = SiO2I(i-1,j-1,k-1)*correctOg;
            COII(i-1,j-1,k-1) =  COII(i-1,j-1,k-1)*correctOg;
            OII(i-1,j-1,k-1) =   OII(i-1,j-1,k-1)*correctOg;
            OHII(i-1,j-1,k-1) =  OHII(i-1,j-1,k-1)*correctOg;
            H2OII(i-1,j-1,k-1) = H2OII(i-1,j-1,k-1)*correctOg;
            H3OII(i-1,j-1,k-1) = H3OII(i-1,j-1,k-1)*correctOg;
            O2II(i-1,j-1,k-1) =  O2II(i-1,j-1,k-1)*correctOg;
            if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
              if (my_chemistry->dust_species > 0)  {
                totalOd = 48./100.* MgSiO3(i-1,j-1,k-1);
              }
              if (my_chemistry->dust_species > 1)  {
                totalOd = totalOd
                        +  64./140.*Mg2SiO4(i-1,j-1,k-1)
                        + 64./232.*  Fe3O4(i-1,j-1,k-1)
                        + 32./ 60.*  SiO2D(i-1,j-1,k-1)
                        + 16./ 40.*    MgO(i-1,j-1,k-1)
                        + 48./102.*  Al2O3(i-1,j-1,k-1);
              }
              if (my_chemistry->dust_species > 2)  {
                totalOd = totalOd
                        +  8./22.68*reforg(i-1,j-1,k-1)
                        + 16./32.  *volorg(i-1,j-1,k-1)
                        + 16./18.  *H2Oice(i-1,j-1,k-1);
              }
              correctOd = (gr_float)(Od[i-1]/totalOd );
              if (my_chemistry->dust_species > 0)  {
                MgSiO3(i-1,j-1,k-1) =  MgSiO3(i-1,j-1,k-1)*correctOd;
              }
              if (my_chemistry->dust_species > 1)  {
                Mg2SiO4(i-1,j-1,k-1) = Mg2SiO4(i-1,j-1,k-1)*correctOd;
                Fe3O4(i-1,j-1,k-1) =   Fe3O4(i-1,j-1,k-1)*correctOd;
                SiO2D(i-1,j-1,k-1) =   SiO2D(i-1,j-1,k-1)*correctOd;
                MgO(i-1,j-1,k-1) =     MgO(i-1,j-1,k-1)*correctOd;
                Al2O3(i-1,j-1,k-1) =   Al2O3(i-1,j-1,k-1)*correctOd;
              }
              if (my_chemistry->dust_species > 2)  {
                reforg(i-1,j-1,k-1) =  reforg(i-1,j-1,k-1)*correctOd;
                volorg(i-1,j-1,k-1) =  volorg(i-1,j-1,k-1)*correctOd;
                H2Oice(i-1,j-1,k-1) =  H2Oice(i-1,j-1,k-1)*correctOd;
              }
            }
  
            totalCg =                       CI(i-1,j-1,k-1)
                    +                      CII(i-1,j-1,k-1)
                    + 12./28.*  CO(i-1,j-1,k-1)
                    + 12./44.* CO2(i-1,j-1,k-1)
                    + 12./13.*  CH(i-1,j-1,k-1)
                    + 12./14.* CH2(i-1,j-1,k-1)
                    + 12./28.*COII(i-1,j-1,k-1);
            correctCg = (gr_float)(Cg[i-1]/totalCg );
            CI(i-1,j-1,k-1) =   CI(i-1,j-1,k-1)*correctCg;
            CII(i-1,j-1,k-1) =  CII(i-1,j-1,k-1)*correctCg;
            CO(i-1,j-1,k-1) =   CO(i-1,j-1,k-1)*correctCg;
            CO2(i-1,j-1,k-1) =  CO2(i-1,j-1,k-1)*correctCg;
            CH(i-1,j-1,k-1) =   CH(i-1,j-1,k-1)*correctCg;
            CH2(i-1,j-1,k-1) =  CH2(i-1,j-1,k-1)*correctCg;
            COII(i-1,j-1,k-1) = COII(i-1,j-1,k-1)*correctCg;
            if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
              if (my_chemistry->dust_species > 0)  {
                totalCd =                           AC(i-1,j-1,k-1);
              }
              if (my_chemistry->dust_species > 2)  {
                totalCd = totalCd
                        + 12./22.68*reforg(i-1,j-1,k-1)
                        + 12./32.  *volorg(i-1,j-1,k-1);
              }
              correctCd = (gr_float)(Cd[i-1]/totalCd );
              if (my_chemistry->dust_species > 0)  {
                AC(i-1,j-1,k-1) =     AC(i-1,j-1,k-1)*correctCd;
              }
              if (my_chemistry->dust_species > 2)  {
                reforg(i-1,j-1,k-1) = reforg(i-1,j-1,k-1)*correctCd;
                volorg(i-1,j-1,k-1) = volorg(i-1,j-1,k-1)*correctCd;
              }
            }
            
            totalSig =                        SiI(i-1,j-1,k-1)
                     + 28./ 44.* SiOI(i-1,j-1,k-1)
                     + 28./ 60.*SiO2I(i-1,j-1,k-1);
            correctSig = (gr_float)(Sig[i-1]/totalSig );
            SiI(i-1,j-1,k-1) =   SiI(i-1,j-1,k-1)*correctSig;
            SiOI(i-1,j-1,k-1) =  SiOI(i-1,j-1,k-1)*correctSig;
            SiO2I(i-1,j-1,k-1) = SiO2I(i-1,j-1,k-1)*correctSig;
            if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
              if (my_chemistry->dust_species > 0)  {
                totalSid = 28./100.* MgSiO3(i-1,j-1,k-1);
              }
              if (my_chemistry->dust_species > 1)  {
                totalSid = totalSid
                         +                          SiM(i-1,j-1,k-1)
                         + 28./140.*Mg2SiO4(i-1,j-1,k-1)
                         + 28./ 60.*  SiO2D(i-1,j-1,k-1);
              }
              correctSid = (gr_float)(Sid[i-1]/totalSid );
              if (my_chemistry->dust_species > 0)  {
                MgSiO3(i-1,j-1,k-1) =  MgSiO3(i-1,j-1,k-1)*correctSid;
              }
              if (my_chemistry->dust_species > 1)  {
                SiM(i-1,j-1,k-1) =     SiM(i-1,j-1,k-1)*correctSid;
                Mg2SiO4(i-1,j-1,k-1) = Mg2SiO4(i-1,j-1,k-1)*correctSid;
                SiO2D(i-1,j-1,k-1) =   SiO2D(i-1,j-1,k-1)*correctSid;
              }
            }

            if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
              if (my_chemistry->dust_species > 1)  {
                totalFeg = Fe(i-1,j-1,k-1);
                correctFeg = (gr_float)(Feg[i-1]/totalFeg );
                Fe(i-1,j-1,k-1) =    Fe(i-1,j-1,k-1)*correctFeg;

                totalFed =                        FeM(i-1,j-1,k-1)
                         +168./232.*Fe3O4(i-1,j-1,k-1)
                         + 56./ 88.*  FeS(i-1,j-1,k-1);
                correctFed = (gr_float)(Fed[i-1]/totalFed );
                FeM(i-1,j-1,k-1) =   FeM(i-1,j-1,k-1)*correctFed;
                Fe3O4(i-1,j-1,k-1) = Fe3O4(i-1,j-1,k-1)*correctFed;
                FeS(i-1,j-1,k-1) =   FeS(i-1,j-1,k-1)*correctFed;
              }

              if (my_chemistry->dust_species > 0)  {
                totalMgg =                   Mg(i-1,j-1,k-1);
                correctMgg = (gr_float)( Mgg[i-1]/totalMgg );
                Mg(i-1,j-1,k-1)      = Mg(i-1,j-1,k-1)     *correctMgg;
                totalMgd = 24./100.* MgSiO3(i-1,j-1,k-1);
                if (my_chemistry->dust_species > 1)  {
                  totalMgd = totalMgd
                           + 48./140.*Mg2SiO4(i-1,j-1,k-1)
                           + 24./ 40.*    MgO(i-1,j-1,k-1);
                }
                correctMgd = (gr_float)( Mgd[i-1]/totalMgd );
                MgSiO3(i-1,j-1,k-1)  = MgSiO3(i-1,j-1,k-1) *correctMgd;
                if (my_chemistry->dust_species > 1)  {
                  Mg2SiO4(i-1,j-1,k-1) = Mg2SiO4(i-1,j-1,k-1)*correctMgd;
                  MgO(i-1,j-1,k-1)     = MgO(i-1,j-1,k-1)    *correctMgd;
                }
              }

              if (my_chemistry->dust_species > 1)  {
                S(i-1,j-1,k-1) =                         Sg[i-1];
                FeS(i-1,j-1,k-1) = 88. / 32. * Sd[i-1];

                Al(i-1,j-1,k-1) =                         Alg[i-1];
                Al2O3(i-1,j-1,k-1) =  102./54. * Ald[i-1];
              }
            }

          } else {

            totalO  = 16./28.*   CO(i-1,j-1,k-1)
                    + 32./44.*  CO2(i-1,j-1,k-1)
                    +                        OI(i-1,j-1,k-1)
                    + 16./17.*   OH(i-1,j-1,k-1)
                    + 16./18.*  H2O(i-1,j-1,k-1)
                    +                        O2(i-1,j-1,k-1)
                    + 16./44.* SiOI(i-1,j-1,k-1)
                    + 32./60.*SiO2I(i-1,j-1,k-1)
                    + 16./28.* COII(i-1,j-1,k-1)
                    +                       OII(i-1,j-1,k-1)
                    + 16./17.* OHII(i-1,j-1,k-1)
                    + 16./18.*H2OII(i-1,j-1,k-1)
                    + 16./19.*H3OII(i-1,j-1,k-1)
                    +                      O2II(i-1,j-1,k-1);
            if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
              if (my_chemistry->dust_species > 0)  {
                totalO  = totalO
                        + 48./100.* MgSiO3(i-1,j-1,k-1);
              }
              if (my_chemistry->dust_species > 1)  {
                totalO  = totalO
                        + 64./140.*Mg2SiO4(i-1,j-1,k-1)
                        + 64./232.*  Fe3O4(i-1,j-1,k-1)
                        + 32./ 60.*  SiO2D(i-1,j-1,k-1)
                        + 16./ 40.*    MgO(i-1,j-1,k-1)
                        + 48./102.*  Al2O3(i-1,j-1,k-1);
              }
              if (my_chemistry->dust_species > 2)  {
                totalO  = totalO
                        +  8./22.68*reforg(i-1,j-1,k-1)
                        + 16./32.  *volorg(i-1,j-1,k-1)
                        + 16./18.  *H2Oice(i-1,j-1,k-1);
              }
            }
            if ( ( my_chemistry->grain_growth == 0 )  &&  ( my_chemistry->dust_sublimation == 0) )  {
              correctO  = (gr_float)(Og[i-1]/ totalO );
              CO(i-1,j-1,k-1)     =    CO(i-1,j-1,k-1)*correctO;
              CO2(i-1,j-1,k-1)     =   CO2(i-1,j-1,k-1)*correctO;
              OI(i-1,j-1,k-1)     =    OI(i-1,j-1,k-1)*correctO;
              OH(i-1,j-1,k-1)     =    OH(i-1,j-1,k-1)*correctO;
              H2O(i-1,j-1,k-1)     =   H2O(i-1,j-1,k-1)*correctO;
              O2(i-1,j-1,k-1)     =    O2(i-1,j-1,k-1)*correctO;
              SiOI(i-1,j-1,k-1)     =  SiOI(i-1,j-1,k-1)*correctO;
              SiO2I(i-1,j-1,k-1)     = SiO2I(i-1,j-1,k-1)*correctO;
              COII(i-1,j-1,k-1)     =  COII(i-1,j-1,k-1)*correctO;
              OII(i-1,j-1,k-1)     =   OII(i-1,j-1,k-1)*correctO;
              OHII(i-1,j-1,k-1)     =  OHII(i-1,j-1,k-1)*correctO;
              H2OII(i-1,j-1,k-1)     = H2OII(i-1,j-1,k-1)*correctO;
              H3OII(i-1,j-1,k-1)     = H3OII(i-1,j-1,k-1)*correctO;
              O2II(i-1,j-1,k-1)     =  O2II(i-1,j-1,k-1)*correctO;
            } else {
              correctO  = (gr_float)(Ot[i-1]/ totalO );
              CO(i-1,j-1,k-1)     =    CO(i-1,j-1,k-1)*correctO;
              CO2(i-1,j-1,k-1)     =   CO2(i-1,j-1,k-1)*correctO;
              OI(i-1,j-1,k-1)     =    OI(i-1,j-1,k-1)*correctO;
              OH(i-1,j-1,k-1)     =    OH(i-1,j-1,k-1)*correctO;
              H2O(i-1,j-1,k-1)     =   H2O(i-1,j-1,k-1)*correctO;
              O2(i-1,j-1,k-1)     =    O2(i-1,j-1,k-1)*correctO;
              SiOI(i-1,j-1,k-1)     =  SiOI(i-1,j-1,k-1)*correctO;
              SiO2I(i-1,j-1,k-1)     = SiO2I(i-1,j-1,k-1)*correctO;
              COII(i-1,j-1,k-1)     =  COII(i-1,j-1,k-1)*correctO;
              OII(i-1,j-1,k-1)     =   OII(i-1,j-1,k-1)*correctO;
              OHII(i-1,j-1,k-1)     =  OHII(i-1,j-1,k-1)*correctO;
              H2OII(i-1,j-1,k-1)     = H2OII(i-1,j-1,k-1)*correctO;
              H3OII(i-1,j-1,k-1)     = H3OII(i-1,j-1,k-1)*correctO;
              O2II(i-1,j-1,k-1)     =  O2II(i-1,j-1,k-1)*correctO;
              if (my_chemistry->dust_species > 0)  {
                MgSiO3(i-1,j-1,k-1) =  MgSiO3(i-1,j-1,k-1)*correctO;
              }
              if (my_chemistry->dust_species > 1)  {
                Mg2SiO4(i-1,j-1,k-1) = Mg2SiO4(i-1,j-1,k-1)*correctO;
                Fe3O4(i-1,j-1,k-1) =   Fe3O4(i-1,j-1,k-1)*correctO;
                SiO2D(i-1,j-1,k-1) =   SiO2D(i-1,j-1,k-1)*correctO;
                MgO(i-1,j-1,k-1) =     MgO(i-1,j-1,k-1)*correctO;
                Al2O3(i-1,j-1,k-1) =   Al2O3(i-1,j-1,k-1)*correctO;
              }
              if (my_chemistry->dust_species > 2)  {
                reforg(i-1,j-1,k-1) =  reforg(i-1,j-1,k-1)*correctO;
                volorg(i-1,j-1,k-1) =  volorg(i-1,j-1,k-1)*correctO;
                H2Oice(i-1,j-1,k-1) =  H2Oice(i-1,j-1,k-1)*correctO;
              }
            }
  
            totalC  =                           CI(i-1,j-1,k-1)
                    +                          CII(i-1,j-1,k-1)
                    + 12./28.    *  CO(i-1,j-1,k-1)
                    + 12./44.    * CO2(i-1,j-1,k-1)
                    + 12./13.    *  CH(i-1,j-1,k-1)
                    + 12./14.    * CH2(i-1,j-1,k-1)
                    + 12./28.    *COII(i-1,j-1,k-1);
            if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1 ) )  {
              if (my_chemistry->dust_species > 0)  {
                totalC  = totalC
                        +                           AC(i-1,j-1,k-1);
              }
              if (my_chemistry->dust_species > 2)  {
                totalC  = totalC
                        + 12./22.68*reforg(i-1,j-1,k-1)
                        + 12./32.  *volorg(i-1,j-1,k-1);
              }
            }
            if ( ( my_chemistry->grain_growth == 0 )  &&  ( my_chemistry->dust_sublimation == 0 ) )  {
              correctC = (gr_float)(Cg[i-1]/ totalC );
              CI(i-1,j-1,k-1) =   CI(i-1,j-1,k-1)*correctC;
              CII(i-1,j-1,k-1) =  CII(i-1,j-1,k-1)*correctC;
              CO(i-1,j-1,k-1) =   CO(i-1,j-1,k-1)*correctC;
              CO2(i-1,j-1,k-1) =  CO2(i-1,j-1,k-1)*correctC;
              CH(i-1,j-1,k-1) =   CH(i-1,j-1,k-1)*correctC;
              CH2(i-1,j-1,k-1) =  CH2(i-1,j-1,k-1)*correctC;
              COII(i-1,j-1,k-1) = COII(i-1,j-1,k-1)*correctC;
            } else {
              correctC  = (gr_float)(  Ct[i-1]/ totalC );
              CI(i-1,j-1,k-1) =     CI(i-1,j-1,k-1)*correctC;
              CII(i-1,j-1,k-1) =    CII(i-1,j-1,k-1)*correctC;
              CO(i-1,j-1,k-1) =     CO(i-1,j-1,k-1)*correctC;
              CO2(i-1,j-1,k-1) =    CO2(i-1,j-1,k-1)*correctC;
              CH(i-1,j-1,k-1) =     CH(i-1,j-1,k-1)*correctC;
              CH2(i-1,j-1,k-1) =    CH2(i-1,j-1,k-1)*correctC;
              COII(i-1,j-1,k-1) =   COII(i-1,j-1,k-1)*correctC;
              if (my_chemistry->dust_species > 0)  {
                AC(i-1,j-1,k-1) =     AC(i-1,j-1,k-1)*correctC;
              }
              if (my_chemistry->dust_species > 2)  {
                reforg(i-1,j-1,k-1) = reforg(i-1,j-1,k-1)*correctC;
                volorg(i-1,j-1,k-1) = volorg(i-1,j-1,k-1)*correctC;
              }
            }
            
            totalSi =                       SiI(i-1,j-1,k-1)
                    + 28./ 44.*SiOI(i-1,j-1,k-1)
                    + 28./ 60.*SiO2I(i-1,j-1,k-1);
            if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1 ) )  {
              if (my_chemistry->dust_species > 0)  {
                totalSi = totalSi
                        + 28./100.* MgSiO3(i-1,j-1,k-1);
              }
              if (my_chemistry->dust_species > 1)  {
                totalSi = totalSi
                        +                          SiM(i-1,j-1,k-1)
                        + 28./140.*Mg2SiO4(i-1,j-1,k-1)
                        + 28./ 60.*  SiO2D(i-1,j-1,k-1);
              }
            }
            if ( ( my_chemistry->grain_growth == 0 )  &&  ( my_chemistry->dust_sublimation == 0 ) )  {
              correctSi = (gr_float)(Sig[i-1]/totalSi );
              SiI(i-1,j-1,k-1) =   SiI(i-1,j-1,k-1)*correctSi;
              SiOI(i-1,j-1,k-1) =  SiOI(i-1,j-1,k-1)*correctSi;
              SiO2I(i-1,j-1,k-1) = SiO2I(i-1,j-1,k-1)*correctSi;
            } else {
              correctSi = (gr_float)(Sit[i-1]/totalSi );
              SiI(i-1,j-1,k-1) =   SiI(i-1,j-1,k-1)*correctSi;
              SiOI(i-1,j-1,k-1) =  SiOI(i-1,j-1,k-1)*correctSi;
              SiO2I(i-1,j-1,k-1) = SiO2I(i-1,j-1,k-1)*correctSi;
              if (my_chemistry->dust_species > 0)  {
                MgSiO3(i-1,j-1,k-1) =  MgSiO3(i-1,j-1,k-1)*correctSi;
              }
              if (my_chemistry->dust_species > 1)  {
                SiM(i-1,j-1,k-1) =     SiM(i-1,j-1,k-1)*correctSi;
                Mg2SiO4(i-1,j-1,k-1) = Mg2SiO4(i-1,j-1,k-1)*correctSi;
                SiO2D(i-1,j-1,k-1) =   SiO2D(i-1,j-1,k-1)*correctSi;
              }
            }

            if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1 ) )  {
              if (my_chemistry->dust_species > 1)  {
                totalFe =                         Fe(i-1,j-1,k-1)
                        +                        FeM(i-1,j-1,k-1)
                        +168./232.*Fe3O4(i-1,j-1,k-1)
                        + 56./ 88.*  FeS(i-1,j-1,k-1);
                correctFe = (gr_float)( Fet[i-1]/totalFe );
                Fe(i-1,j-1,k-1)    =    Fe(i-1,j-1,k-1)*correctFe;
                FeM(i-1,j-1,k-1)   =   FeM(i-1,j-1,k-1)*correctFe;
                Fe3O4(i-1,j-1,k-1) = Fe3O4(i-1,j-1,k-1)*correctFe;
                FeS(i-1,j-1,k-1)   =   FeS(i-1,j-1,k-1)*correctFe;
              }

              if (my_chemistry->dust_species > 0)  {
                totalMg =                           Mg(i-1,j-1,k-1)
                        + 24./100.* MgSiO3(i-1,j-1,k-1);
                if (my_chemistry->dust_species > 1)  {
                  totalMg = totalMg
                          + 48./140.*Mg2SiO4(i-1,j-1,k-1)
                          + 24./ 40.*    MgO(i-1,j-1,k-1);
                }
                correctMg = (gr_float)( Mgt[i-1]/totalMg );
                Mg(i-1,j-1,k-1) =      Mg(i-1,j-1,k-1)*correctMg;
                MgSiO3(i-1,j-1,k-1) =  MgSiO3(i-1,j-1,k-1)*correctMg;
                if (my_chemistry->dust_species > 1)  {
                  Mg2SiO4(i-1,j-1,k-1) = Mg2SiO4(i-1,j-1,k-1)*correctMg;
                  MgO(i-1,j-1,k-1) =     MgO(i-1,j-1,k-1)*correctMg;
                }
              }

              if (my_chemistry->dust_species > 1)  {
                totalS  =                        S(i-1,j-1,k-1)
                        + 32./ 88.*FeS(i-1,j-1,k-1);
                correctS  = (gr_float)(  St[i-1]/totalS  );
                S(i-1,j-1,k-1) =   S(i-1,j-1,k-1)*correctS;
                FeS(i-1,j-1,k-1) = FeS(i-1,j-1,k-1)*correctS;

                totalAl =                         Al(i-1,j-1,k-1)
                        + 54./102.*Al2O3(i-1,j-1,k-1);
                correctAl = (gr_float)( Alt[i-1]/totalAl );
                Al(i-1,j-1,k-1) =    Al(i-1,j-1,k-1)*correctAl;
                Al2O3(i-1,j-1,k-1) = Al2O3(i-1,j-1,k-1)*correctAl;
              }
            }

          }

          //    CI(i,j,k)      = max(CI(i,j,k), tiny)
          //    CII(i,j,k)     = max(CII(i,j,k), tiny)
          //    CO(i,j,k)      = max(CO(i,j,k), tiny)
          //    CO2(i,j,k)     = max(CO2(i,j,k), tiny)
          //    OI(i,j,k)      = max(OI(i,j,k), tiny)
          //    OH(i,j,k)      = max(OH(i,j,k), tiny)
          //    H2O(i,j,k)     = max(H2O(i,j,k), tiny)
          //    O2(i,j,k)      = max(O2(i,j,k), tiny)
          //    SiI(i,j,k)     = max(SiI(i,j,k), tiny)
          //    SiOI(i,j,k)    = max(SiOI(i,j,k), tiny)
          //    SiO2I(i,j,k)   = max(SiO2I(i,j,k), tiny)
          //    CH(i,j,k)      = max(CH(i,j,k), tiny)
          //    CH2(i,j,k)     = max(CH2(i,j,k), tiny)
          //    COII(i,j,k)    = max(COII(i,j,k), tiny)
          //    OII(i,j,k)     = max(OII(i,j,k), tiny)
          //    OHII(i,j,k)    = max(OHII(i,j,k), tiny)
          //    H2OII(i,j,k)   = max(H2OII(i,j,k), tiny)
          //    H3OII(i,j,k)   = max(H3OII(i,j,k), tiny)
          //    O2II(i,j,k)    = max(O2II(i,j,k), tiny)
          // if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1 ) ) then
          // if (idspecies .gt. 0) then
          //    Mg(i,j,k)      = max(Mg(i,j,k), tiny)
          //    MgSiO3(i,j,k)  = max(MgSiO3(i,j,k), tiny)
          //    AC(i,j,k)      = max(AC(i,j,k), tiny)
          // endif
          // if (idspecies .gt. 1) then
          //    Al(i,j,k)      = max(Al(i,j,k), tiny)
          //    S(i,j,k)       = max(S(i,j,k), tiny)
          //    Fe(i,j,k)      = max(Fe(i,j,k), tiny)
          //    SiM(i,j,k)     = max(SiM(i,j,k), tiny)
          //    FeM(i,j,k)     = max(FeM(i,j,k), tiny)
          //    Mg2SiO4(i,j,k) = max(Mg2SiO4(i,j,k), tiny)
          //    Fe3O4(i,j,k)   = max(Fe3O4(i,j,k), tiny)
          //    SiO2D(i,j,k)   = max(SiO2D(i,j,k), tiny)
          //    MgO(i,j,k)     = max(MgO(i,j,k), tiny)
          //    FeS(i,j,k)     = max(FeS(i,j,k), tiny)
          //    Al2O3(i,j,k)   = max(Al2O3(i,j,k), tiny)
          // endif
          // if (idspecies .gt. 2) then
          //    reforg(i,j,k)  = max(reforg(i,j,k), tiny)
          //    volorg(i,j,k)  = max(volorg(i,j,k), tiny)
          //    H2Oice(i,j,k)  = max(H2Oice(i,j,k), tiny)
          // endif
          // endif

          // endif
        }
      }

      //      if ( (idustfield .gt. 0) .and. (idspecies .gt. 0) ) then
      //         do i = is+1, ie+1
      // !          if ( itmask_metal(i) ) then
      //            if (idspecies .gt. 0) then
      //               dust(i,j,k) = MgSiO3  (i,j,k)
      //     &                     + AC      (i,j,k)
      //            endif
      //            if (idspecies .gt. 1) then
      //               dust(i,j,k) = dust(i,j,k)
      //     &                     + SiM     (i,j,k)
      //     &                     + FeM     (i,j,k)
      //     &                     + Mg2SiO4 (i,j,k)
      //     &                     + Fe3O4   (i,j,k)
      //     &                     + SiO2D   (i,j,k)
      //     &                     + MgO     (i,j,k)
      //     &                     + FeS     (i,j,k)
      //     &                     + Al2O3   (i,j,k)
      //            endif
      //            if (idspecies .gt. 2) then
      //               dust(i,j,k) = dust(i,j,k)
      //     &                     + reforg  (i,j,k)
      //     &                     + volorg  (i,j,k)
      //     &                     + H2Oice  (i,j,k)
      //            endif
      // !          endif
      //         enddo
      //      endif

      // Set the electron density

      for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
        de (i-1,j-1,k-1) = HII(i-1,j-1,k-1) + HeII(i-1,j-1,k-1)/(gr_float)(4.) +
             HeIII(i-1,j-1,k-1)/(gr_float)(2.);
        if (my_chemistry->primordial_chemistry > 1) { de(i-1,j-1,k-1) = de(i-1,j-1,k-1)
             - HM(i-1,j-1,k-1) + H2II(i-1,j-1,k-1)/(gr_float)(2.); }
        if (my_chemistry->primordial_chemistry > 3) { de(i-1,j-1,k-1) = de(i-1,j-1,k-1)
             - DM   (i-1,j-1,k-1)/(gr_float)(2.)
             + HDII (i-1,j-1,k-1)/(gr_float)(3.)
             + HeHII(i-1,j-1,k-1)/(gr_float)(5.); }
        if (my_chemistry->metal_chemistry == 1)  {
          // if (itmask_metal(i)) then
          de(i-1,j-1,k-1) = de(i-1,j-1,k-1)
            + CII  (i-1,j-1,k-1)/(gr_float)(12.)
            + COII (i-1,j-1,k-1)/(gr_float)(28.)
            + OII  (i-1,j-1,k-1)/(gr_float)(16.)
            + OHII (i-1,j-1,k-1)/(gr_float)(17.)
            + H2OII(i-1,j-1,k-1)/(gr_float)(18.)
            + H3OII(i-1,j-1,k-1)/(gr_float)(19.)
            + O2II (i-1,j-1,k-1)/(gr_float)(32.);
          // endif
        }
      }

    }
  }

  return;
}

} // namespace grackle::impl

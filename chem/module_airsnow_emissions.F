! 2018/12/06, Louis Marelle
!
! Purpose:
! Compute Br2 emissions from the snowpack (surface snow) and Br2+NaCl emissions
! from blowing snow, and add emissions to chem
!
! References:
! - Louis Marelle, Jennie L. Thomas, Shaddy Ahmed, Katie Tuite, Jochen Stutz,
! Aurélien Dommergue, William R. Simpson, Markus M.  Frey, Foteini Baladima.
! Implementation and impacts of surface and blowing snow sources of Arctic
! bromine activation within WRF-Chem 4.1.1. Submitted to the Journal of Advances
! in Modeling Earth Systems, 2021
! - Toyota, K., McConnell, J. C., Lupu, A., Neary, L., McLinden, C. A., Richter,
! A., Kwok, R., Semeniuk, K., Kaminski, J. W., Gong, S.-L., Jarosz, J.,
! Chipperfield, M. P., and Sioris, C. E.: Analysis of reactive bromine
! production and ozone depletion in the Arctic boundary layer using 3-D
! simulations with GEM-AQ: inference from synoptic-scale patterns, Atmos. Chem.
! Phys., 11, 3949-3979, https://doi.org/10.5194/acp-11-3949-2011, 2011.
! and
! - Yang, X., Pyle, J. A., and Cox, R. A. ( 2008), Sea salt aerosol production
! and bromine release: Role of snow on sea ice, Geophys. Res. Lett., 35, L16815,
! doi:10.1029/2008GL034536.
! - Déry, S.J. & Yau, M.K. Simulation Of Blowing Snow In The Canadian Arctic
! Using A Double-Moment Model, Boundary-Layer Meteorology (2001) 99: 297.
! https://doi.org/10.1023/A:1018965008049
! -  Huang, J. and Jaeglé, L.: Wintertime enhancements of sea salt aerosol in
! polar regions consistent with a sea ice source from blowing snow, Atmos. Chem.
! Phys., 17, 3699–3712, https://doi.org/10.5194/acp-17-3699-2017, 2017.
!
! Remarks:
! - For now, this is only coupled to the SAPRC99_MOSAIC_8BIN_VBS2_AQ_HALOGENS_KPP
! gas-phase/aerosol mechanism
! - At the very least, the blowing snow code assumes that MOSAIC is used with
! the standard 8 size bins.
! - In order to use the blowing snow code with a different aerosol module, naersizes,
! mosaic_massfrac_nacl and depletion_factor need to be updated to correspond to
! the aerosol sizes in the new module (+ all of the other instances using isize,
! itype, or iphase must be updated).
! - In order to use a new snow salinity for blowing snow, the size distribution
! of the emitted seasalt, mosaic_massfrac_nacl, also has to change consistently
!
!TODO add a check for lake ice
!TODO Ideally, blowing snow should be in a separate module
!TODO If someone tries to use this option with another mechanism, include a call
! to WRF_FATAL in chemics_init.F.

MODULE module_airsnow_emissions

CONTAINS

!----------------------------------------------------------------------
  SUBROUTINE airsnow_emissions( id, dtstep, dz8w, moist, t_phy, p_phy, u_phy, &
         v_phy, rho_phy, chem, u10, v10, ust, T2, snowh, tsk,          &
         xice, lakemask, meltedsnow, coszen,                           &
         drydepflx_o3, drydepflx_hobr, drydepflx_brno3,                &
         drydepflx_hbr, drydepflx_br2, emiss_br2,                      &
         ids,ide, jds,jde, kds,kde,                                    &
         ims,ime, jms,jme, kms,kme,                                    &
         its,ite, jts,jte, kts,kte,                                    &
         snow_emissions_opt, blowing_snow_opt                          ) 

  USE module_configure
  USE module_state_description, ONLY: param_first_scalar
  USE module_data_mosaic_asect, ONLY: dcen_sect, numptr_aer, lptr_cl_aer, &
          lptr_na_aer, nsize_aer, ai_phase
  USE module_data_mosaic_constants, ONLY:  pi
  USE module_model_constants, ONLY: mwdry, ep_2, r_v
  USE module_mp_morr_two_moment, ONLY: POLYSVP
  USE module_wrf_error ! For wrf_error_fatal

  IMPLICIT NONE


  !-------- Subroutine arguments  --------
  INTEGER, INTENT(IN) :: id,                                      &
                         ids,ide, jds,jde, kds,kde,               &
                         ims,ime, jms,jme, kms,kme,               &
                         its,ite, jts,jte, kts,kte
  ! Logical flag for snow emissions and blowing snow options
  LOGICAL, INTENT(IN) :: snow_emissions_opt, blowing_snow_opt
  ! Dynamical time step (s)
  REAL, INTENT(IN) :: dtstep
  ! Vertical layer thickness (m), air temperature (K), air pressure (Pa),
  ! x-direction windspeed (m/s), y-direction windspeed (m/s), air density
  ! (km/m3)
  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),                     &
        INTENT(IN) :: dz8w, t_phy, p_phy, u_phy, v_phy, rho_phy
  ! Moist scalars, used to obtain qv, the water vapor mixing ratio (kg/kg)
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),        &
        INTENT(IN ) :: moist
  ! 10m-x-direction-windspeed (m/s), 10m-y-direction-windspeed (m/s), friction
  ! velocity (m/s), 2m-temperature (K), Snow height (m), skin temperature (K),
  ! ice fraction (0-1), lake mask (0-1), cosine of solar zenith angle (dimensionless)
  REAL, DIMENSION(ims:ime, jms:jme)               ,               &
        INTENT(IN) :: u10, v10, ust, t2, snowh, tsk, xice, lakemask, coszen
  ! Past snow melt flag (dimensionless flag)
  REAL, DIMENSION(ims:ime, jms:jme)               ,               &
        INTENT(INOUT) :: meltedsnow
  ! Dry deposition fluxes of O3, HOBr, BrNO3, and HBr (mole/m2/s), calculated in
  ! dry_dep_driver
  REAL,  DIMENSION(ims:ime, jms:jme),                             &
        INTENT(IN) :: drydepflx_o3, drydepflx_hobr,               &
        drydepflx_brno3, drydepflx_hbr, drydepflx_br2
  ! Accumulated br2 emission flux diagnosis (mole/m2)
  REAL, DIMENSION(ims:ime, jms:jme)               ,               &
        INTENT(INOUT) :: emiss_br2
  ! Trace species mixing ratios (gases in ppm, aerosols in microg/kg)
  REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem),           &
        INTENT(INOUT) :: chem


  !-------- Parameters --------
  ! Temperature threshold for Br2 emissions & recycling (K), sea ice cover
  ! threshold for surface snow and blowing snow emissions (area fraction)
  REAL, PARAMETER :: tskcrit = 273.15, seaice_crit = 0.75

  !---- Parameters for the blowing snow model
  INTEGER, PARAMETER :: naersizes = 8
  ! Br/NaCl mass ratio in sea salt (kg/kg, Yang et al., 2008)
  REAL, PARAMETER :: br_nacl_massratio = 0.00233
  ! Average Antarctic snow salinity (psu)
  ! The size distribution of the emitted seasalt, mosaic_massfrac_nacl, also
  ! depends on snow_salinity. Do not update one without updating the other.
  ! ! Massom et al. (2001), Antarctic upper bound
  ! REAL :: snow_salinity = 8.5
  ! Huang and Jaeglé, 2017, Arctic
  REAL, PARAMETER :: snow_salinity = 0.1
  ! Conversion factor (from kg/m2/day to kg/m2/s)
  REAL, PARAMETER :: kconv_day_s = 1.1574E-5
  ! Snow age (hours), Huang and Jaeglé (2017)
  REAL, PARAMETER :: snow_age = 72
  ! Reference 10m-wind speed (m/s)
  REAL, PARAMETER :: w10_0 = 6.975
  ! Parameters for qsnow_subl_prime parameterization, from Dery and Yau (2001)
  REAL, PARAMETER :: a0 = 3.78407E-1, a1 = -8.64089E-2, a2 = -1.60570E-2, &
          a3 = 7.25516E-4, a4 = -1.25650E-1, a5 = 2.48430E-2, &
          a6 = -9.56871E-4, a7 = 1.24600E-2, a8 = 1.56862E-3, &
          a9 = -2.93002E-4
  ! Conversion factor for the ksi parameter (unitless)
  REAL, PARAMETER :: convert_ksi = -1.0E12
  ! Density of ice (kg/m3) and sea salt (kg/m3)
  REAL, PARAMETER :: rho_ice = 914.0, rho_nacl = 2160.
  ! Molar weight of bromine (g/mole)
  REAL, PARAMETER :: mw_br = 79.904
  ! Na and Cl molar weights in g/mol
  REAL, PARAMETER :: MW_NA=22.990, MW_CL=35.453


  !-------- Local variables --------
  CHARACTER (LEN=1000) :: msg
  ! Loop indices i j k
  INTEGER :: i, j, k
  ! Conversion factor from emissions to ppm, conversion factor from emissions to µg/kg
  REAL :: conv, conv2
  ! Switches for enabling surface snow or blowing snow emissions
  LOGICAL :: do_snow_emissions, do_blowing_snow
  ! Br2 and BrCl emissions (mole/m2/s)
  REAL :: br2_emi = 0.0, brcl_emi = 0.0

  !---- Local variables for the blowing snow model
  ! Br emissions (kg/m2/s), NaCl emissions (kg/m2/s), Snow sublimation flux (kg/m2/s), Normalized
  ! column-integrated snow sublimation flux (kg/m2/day), Snow mixing ratio in
  ! saltation layer (kg/kg), Snow mixing ratio at reference wind speed (kg/kg)
  REAL :: br_emi, nacl_emi, qsnow_subl, qsnow_subl_prime, qsnow_saltation, qsnow_0
  ! 10m-windspeed (m/s), critical 10m-windspeed (m/s), temperature-moisture
  ! parameter from Dery and Yau 2001 (-1.0E-12 m2/s-1), ice saturation mixing
  ! ratio (kg/kg), conductivity term in the sublimation process (m s kg−1),
  ! diffusion term in the sublimation process (m s kg−1), saturation water vapor
  ! pressure wrt ice (Pa)
  REAL :: w10, w10_crit, ksi, qvapor_sat_ice, fconduction, fdiffusion, eis
  ! Aerosol number concentration change (#/m2/s), Conversion factor from
  ! number emissions to #/kg
  REAL :: aer_num_change, conv3
  ! Molecular diffusion coefficient of water vapor (m2 s-1), latent heat of
  ! sublimation (J kg-1), coefficient of thermal conductivity of dry air (J m-1
  ! s-1 K-1),
  REAL :: d_v, l_s, k_d
  ! MOSAIC aerosol type, size, phase indexes, and "pointers" to number, Na and
  ! Cl in chem
  INTEGER :: itype, isize, iphase, l_num, l_na, l_cl
  ! Adjustment factor for the emission flux (unitless)
  REAL ::  a_prime
  ! Mass proportion of total NaCl emissions going into each MOSAIC size bin,
  ! calculated offline from the size+salinity integral in Yang et al.
  ! (2008), eq. 8 TODO include a check on naersizes
  REAL :: mosaic_massfrac_nacl(naersizes)
  ! Bromine depletion factor in sea salt (unitless), from Sander et al. [2003],
  ! interpolated to the MOSAIC center size bin
  !LMarelle - this is disabled for now, see below
  REAL :: depletion_factor(naersizes)


  !---- Data declaration
  ! Calculated based on Yang et al. (2018), eq. 8, for the MOSAIC size bins
  ! ! Based on snow salinity from Massom et al. (2001)
  ! mosaic_massfrac_nacl = (/ 0.2189E-09, 0.6873E-08, 0.2122E-06, 0.6306E-05, &
  !                           0.1742E-03, 0.4182E-02, 0.7729E-01, 0.9183E+00 /)
  ! Based on 0.1 psu uniform snow salinity and N=5 aerosol fractioning, Huang and Jaeglé (2017)
  mosaic_massfrac_nacl = (/ 0.2711E-01, 0.8600E-01, 0.2181E+00, 0.3580E+00, &
                            0.2629E+00, 0.4721E-01, 0.7154E-03, 0.9098E-07 /)
  ! depletion_factor = (/ -3.8200, -3.8200, -3.6856, -2.2739, &
  !                        0.0965, 0.3379, 0.2332, 0.1393  /)
  ! With realistic Arctic snow salinities, most of the sea salt is
  ! emitted in submicron sizes, where depletion factors are negative.
  ! For this reason it does not seem to make sense to use size resolved
  ! depletion factors to calculate instantaneous Br emissions.
  ! Replace depletion factors by the upper bound of 0.38 in Yang et al., (2008).
  !TODO - use explicit br aerosols instead, and br2 emissions by heterogeneous
  ! reactions (the only thing already included is heterogeneous recycling on
  ! aerosols)
  depletion_factor = (/ 0.38, 0.38, 0.38, 0.38, &
                        0.38, 0.38, 0.38, 0.38  /)
  do_snow_emissions = snow_emissions_opt
  do_blowing_snow = blowing_snow_opt
  itype = 1
  iphase = ai_phase


  !-------- Add emissions from snow, recycling and blowing snow --------
  ! From Toyota et al. (2011) and Yang et al. (2008)

  DO j = jts,jte
    DO i = its,ite
      br2_emi = 0.0
      brcl_emi = 0.0
      nacl_emi = 0.0
      aer_num_change = 0.0

      a_prime = 1.0 / (1.038 + 0.03758 * snow_age - 0.00014349 * snow_age ** 2.0 &
                + 1.911315E-7 * snow_age ** 3.0)

      !-------- Halogen emissions from surface snow and recycling over snow
      ! Only valid over snow or over sea ice (assuming that ice is always snow
      ! covered below tskcrit
      ! The snowh threshold needs to be consistent with the one used in the dry
      ! deposition scheme to trigger snow-covered deposition.
      IF (((snowh(i, j) .GT. 0.01 .AND. xice(i, j) .LE. 0.01) &
           .OR. xice(i, j) .GT. seaice_crit) &
          .AND. meltedsnow(i,j) .LE. 0.5) THEN
        ! If snow melts, turn off surface emissions there for the rest of the
        ! simulation 
        !TODO this needs to be reset in the fall, e.g. when sea-ice reforms in
        !october and later or with new snow fall !TODO this should ideally drift
        !with sea ice !TODO this should ideally read in observed sea ice state
        !(snow state, snow salinity etc.)
        IF(tsk(i, j) .GT. tskcrit) THEN
          meltedsnow(i,j) = 1.
        ENDIF

        ! Br2 emissions occur only when the surface temperature is below TSKcrit
        IF (tsk(i, j) .LE. tskcrit) THEN

          !----- Br2 "trigger" release

          IF (do_snow_emissions) THEN
            ! O3 deposition over first-year sea ice is converted to Br2 with a
            ! molar yield of 0.075 for SZA <= 85°
            IF (coszen(i,j) .GE. 0.087 .AND. xice(i,j) .GT. seaice_crit &
                .AND. lakemask(i,j) .LT. 0.5) THEN
              br2_emi = br2_emi + drydepflx_o3(i,j) * 0.075 * xice(i,j)
            ELSEIF (coszen(i,j) .LT. 0.087 .AND. xice(i,j) .GT. seaice_crit &
                    .AND. lakemask(i,j) .LT. 0.5) THEN
              br2_emi = br2_emi + drydepflx_o3(i,j) * 0.001 * xice(i,j)

            END IF ! coszen >= 0.087 <=> SZA <= 85°, and over ice
          END IF ! do_snow_emissions

          !---- Br2 "recycling" release

          ! Over sea ice
          IF (xice(i,j) .GT. seaice_crit .AND. lakemask(i,j) .LT. 0.5) THEN
            ! Over FY ice, HOBr and BrONO2 deposition gets converted to Br2
            br2_emi = br2_emi + (drydepflx_hobr(i,j) + drydepflx_brno3(i,j)) * xice(i,j)

          ! ! This is disabled for now, as there is no strong evidence that MYI
          ! and FYI act differently (Peterson et al., 2019)
          ! ELSEIF (xice(i,j) .GT. 0.15 .AND. seaice_age(i,j) >= 1.5) THEN
          !   !Over MY ice, 2 cases
          !   !   1 - if HOBr+BrONO2 deposition <= HBr deposition, it gets converted to
          !   !   Br2
          !   !   2 - if HOBr+BrONO2 deposition > HBr deposition, an amount of [HBr]
          !   !   gets converted to Br2 and the excess part gets converted to BrCl
          !   !   (or 50% Br2 if there is no chlorine in the mechanism)
          !   IF ((drydepflx_hobr(i,j) + drydepflx_brno3(i,j)) .LE. drydepflx_hbr(i,j)) THEN
          !     br2_emi = br2_emi + (drydepflx_hobr(i,j) + drydepflx_brno3(i,j))
          !   ELSE
          !     br2_emi = br2_emi + drydepflx_hbr(i,j)
          !     brcl_emi = brcl_emi + (drydepflx_hobr(i,j) + drydepflx_brno3(i,j) - drydepflx_hbr(i,j))
          !   END IF

          END IF

          ! Over continental snow
          ! Make sure that this is not ocean/seaice. WRF considers sea-ice as
          ! land, so to check for continental snow we use, xice < 0.01 and
          ! snowh>0.01 (there is no snow height over the open ocean)
          !TODO use xland + seaice check, would be more legible
          IF (xice(i,j) .LE. 0.01 .AND. snowh(i, j) .GT. 0.01) THEN
            IF ((drydepflx_hobr(i,j) + drydepflx_brno3(i,j)) .LE. drydepflx_hbr(i,j)) THEN
              br2_emi = br2_emi + (drydepflx_hobr(i,j) + drydepflx_brno3(i,j))
            ELSE
              br2_emi = br2_emi + drydepflx_hbr(i,j)
            END IF
          END IF
        END IF ! tsk <= tskcrit
      END IF ! (snowh(i, j) .GT. 0.01 .OR. xice(i, j) .GT. seaice_crit), snow or sea ice cell

      !-------- Halogen emissions from blowing snow

      IF (do_blowing_snow) THEN
        ! Only valid over sea-ice and for snow that has not experienced recent melt
        IF (xice(i, j) .GT. seaice_crit .AND. meltedsnow(i,j) .LT. 0.5 &
          .AND. lakemask(i,j) .LT. 0.5) THEN
          ! Calculate the critical 10m wind speed for blowing snow, from Yang et
          ! al. (2008)
          w10_crit = w10_0 + 0.0033 * (t2(i,j) - 273.15 + 27.27) ** 2.0

          ! Wind speed in m s-1
          w10 = sqrt( u10(i,j) ** 2.0 + v10(i,j) ** 2.0 )
          ! We donṫ trust the u10,v10 values, is model layers are very
          ! thin near surface (this recommandation is from the GOCART modules)
          if ( dz8w(i,kts,j) .lt. 12.) &
            w10 = sqrt( u_phy(i,kts,j) * u_phy(i,kts,j) + v_phy(i,kts,j) * v_phy(i,kts,j) )

          IF (w10 .gt. w10_crit) THEN
            ! Calculate the snow sublimation flux, qsnow_subl, based on Yang et
            ! al. (2008) and Dery and Yau (2001), with help from Markus Frey for
            ! the fconduction and fdiffusion terms
            l_s = 2.83458E6 - (t2(i,j) - 273.15) * (340. + 10.46 * (t2(i,j) - 273.15))
            k_d = 0.023807 + 7.1128E-5 * (t2(i,j) - 273.15)
            d_v = 2.11E-5 * (t2(i,j) / 273.15) ** 1.94 * (1013.25 / (p_phy(i, kts, j) * 1.0E-2))
            ! qvapor_sat_ice from the morrison microphysics module
            eis = min(0.99 * p_phy(i,kts,j), POLYSVP(t_phy(i,kts,j), 1))
            qvapor_sat_ice = ep_2 * eis / (p_phy(i,kts,j) - eis)
            fconduction = (l_s / (r_v * t2(i,j)) - 1.0) * l_s / (k_d * t2(i,j))
            fdiffusion = (r_v * t2(i,j)) / (d_v * eis)
            ! From Dery and Yau (2001) eq. 11
            ksi = convert_ksi * (moist(i, kts, j, p_qv) / qvapor_sat_ice - 1.0) &
                  / (2.0 * rho_ice * (fconduction + fdiffusion))
            ! From Dery and Yau (2001) eq. 14
            qsnow_subl_prime = a0 + a1 * ksi + a2 * ksi ** 2.0 + a3 * ksi ** 3.0 &
              + a4 * w10 + a5 * ksi * w10 + a6 * ksi ** 2.0 * w10 + a7 * w10 ** 2.0 &
              + a8 * ksi * w10 ** 2 + a9 * w10 ** 3.0
            ! From Yang et al., (2008)
            qsnow_saltation = 0.385 * (1. - w10_crit / w10) ** 2.59 / ust(i,j)
            qsnow_0 = 0.385 * (1. - w10_0 / w10) ** 2.59 / ust(i,j)
            qsnow_subl = kconv_day_s * a_prime * qsnow_subl_prime * qsnow_saltation / qsnow_0

            ! Sea salt emission flux (kg/m2/s), from Yang et al., (2008),
            ! scaled by the ice fraction.
            nacl_emi = qsnow_subl * snow_salinity / 1000. * xice(i, j)
            IF (nacl_emi < 0.) THEN
            !  write(msg,*) 'module_airsnow_emissions: nacl_emi < 0., nacl_emi=',nacl_emi,' at i,j=',i,j
            !  CALL wrf_error_fatal( trim(msg) )
              nacl_emi = 0.
            ENDIF

            !---- Br/Br2 emissions from Yang et al., (2008)
            DO isize = 1,nsize_aer(itype)
              ! (kg/m2/s)
              br_emi = nacl_emi * mosaic_massfrac_nacl(isize) * br_nacl_massratio * depletion_factor(isize)
              ! Convert from kg(br)/m2/s to mole(br2)/m2/s
              br2_emi = br2_emi + br_emi * 1.0E3 / mw_br / 2.0
            END DO ! isize

          END IF ! w10 .gt. w10_crit
        END IF ! xice(i, j) .GT. 0.15 .AND. tsk(i,j) .LT. tskcrit
      END IF ! do_blowing_snow


      !-------- Add Br2 emissions from airsnow and blowing snow, and sea salt emissions from blowing snow  --------
      ! Conversion factor from mole/m2/s to ppm
      conv = 1.0E3 * mwdry / rho_phy(i,kts,j) * dtstep / (dz8w(i,kts,j))
      ! Conversion from kg/m2/s to µg/kg
      conv2 = 1.0E9 / (rho_phy(i,kts,j) * dz8w(i,kts,j)) * dtstep
      ! Conversion from #/m2/s to #/kg
      conv3 = dtstep / (rho_phy(i,kts,j) * (dz8w(i,kts,j)))

      DO isize = 1, nsize_aer(itype)
        !---- Na emissions
        l_na = lptr_na_aer(isize,itype,iphase)
        IF (l_na >= param_first_scalar) THEN
          chem(i,kts,j,l_na) = chem(i,kts,j,l_na) + nacl_emi * MW_NA/(MW_NA+MW_CL) * mosaic_massfrac_nacl(isize) * conv2
        END IF
        !---- Cl emissions
        l_cl = lptr_cl_aer(isize,itype,iphase)
        IF (l_cl >= param_first_scalar) THEN
          chem(i,kts,j,l_cl) = chem(i,kts,j,l_cl) + nacl_emi * MW_CL/(MW_NA+MW_CL) * mosaic_massfrac_nacl(isize) * conv2
        END IF
        !---- Aerosol number "emissions" from the mass change
        l_num = numptr_aer(isize,itype,iphase)
        IF (l_num >= param_first_scalar) THEN
          ! Calculate the change in aerosol number concentrations due to
          ! the sea salt emissions (#/m2/s). dcen_sect is in cm in MOSAIC
          aer_num_change = nacl_emi / rho_nacl * 6.0/pi * 1.0/(dcen_sect(isize,itype)*1.0E-2)**3.0
          ! #/m2/s -> #/kg
          chem(i,kts,j,l_num) = chem(i,kts,j,l_num) + aer_num_change * conv3
        END IF
        !---- Br2 emissions
        IF (chem(i,kts,j,p_br2)  + br2_emi * conv <= 0) THEN
          br2_emi = - chem(i,kts,j,p_br2) / conv 
          chem(i,kts,j,p_br2) = 0.0
        ELSE
          chem(i,kts,j,p_br2) = chem(i,kts,j,p_br2) + br2_emi * conv
        END IF
        ! !---- BrCl emissions
        ! IF (chem(i,kts,j,p_brcl)  + brcl_emi * conv <= 0) THEN
        !   brcl_emi = - chem(i,kts,j,p_brcl) / conv
        !   chem(i,kts,j,p_brcl) = 0.0
        ! ELSE
        !   chem(i,kts,j,p_brcl) = chem(i,kts,j,p_brcl) + brcl_emi * conv
        ! END IF

        !---- Output diagnosis, accumulated br2 emissions in mole/m2
        emiss_br2(i,j) = emiss_br2(i,j) + br2_emi * dtstep

      END DO !isize

    END DO ! i loop
  END DO ! j loop

END SUBROUTINE airsnow_emissions

END MODULE module_airsnow_emissions

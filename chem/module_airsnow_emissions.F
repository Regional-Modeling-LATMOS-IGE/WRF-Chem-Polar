! 2018/12/06, Louis Marelle
!
! Purpose:
! Compute Br2 emissions from the snowpack and add them to chem(:, :, :, p_br2)
!
! References:
! - Toyota, K., McConnell, J. C., Lupu, A., Neary, L., McLinden, C. A., Richter,
! A., Kwok, R., Semeniuk, K., Kaminski, J. W., Gong, S.-L., Jarosz, J.,
! Chipperfield, M. P., and Sioris, C. E.: Analysis of reactive bromine
! production and ozone depletion in the Arctic boundary layer using 3-D
! simulations with GEM-AQ: inference from synoptic-scale patterns, Atmos. Chem.
! Phys., 11, 3949-3979, https://doi.org/10.5194/acp-11-3949-2011, 2011.
! and
! - Falk, S. and Sinnhuber, B.-M.: Polar boundary layer bromine explosion and
! ozone depletion events in the chemistry–climate model EMAC v2.52:
! implementation and evaluation of AirSnow algorithm, Geosci. Model Dev., 11,
! 1115-1131, https://doi.org/10.5194/gmd-11-1115-2018, 2018.
! - Yang, X., Pyle, J. A., and Cox, R. A. ( 2008), Sea salt aerosol production
! and bromine release: Role of snow on sea ice, Geophys. Res. Lett., 35, L16815,
! doi:10.1029/2008GL034536.
! - Déry, S.J. & Yau, M.K. Simulation Of Blowing Snow In The Canadian Arctic
! Using A Double-Moment Model, Boundary-Layer Meteorology (2001) 99: 297.
! https://doi.org/10.1023/A:1018965008049
!TODO Massom et al. (2001) + Sander et al 2003 references
!TODO Blowing snow should be in a separate module, and should use snow salinity
!data for the Arctic instead of the Antarctic data used here
!
! Remarks:
! For now, this is only coupled to the SAPRC99_MOSAIC_8BIN_VBS2_AQ_HALOGENS_KPP
! gas-phase/aerosol mechanism
!
!TODO check blowing snow results vs. Yang2008, fig. 2
!
!TODO scale by snowc, xice??
!
!TODO add checks for numerical issues, eg divisions
!
!TODO, write down assumptions, e.g. it's a MOSAIC module with the standard 8
!size bins, maybe compare the bin sizes to what I assumed offline
!
!TODO define all variables and units in declaration section
! Check all units !
!
!TODO change the critical temperature to 0K but include a memory to this,
!everytime the temperature rises above 0K or if liquid rain falls bromine
!emissions and recycling are turned off, until new snow falls again and
!temperatures are below 0K.
!TODO
! - If someone tries to use this option with another mechanism include a call
! to WRF_FATAL in chemics_init.F.
! - Write a caveat for lake ice, which should not emit Br2.

MODULE module_airsnow_emissions

CONTAINS

!----------------------------------------------------------------------
  SUBROUTINE airsnow_emissions( id, dtstep, dz8w, moist, t_phy, p_phy, u_phy, &
         v_phy, rho_phy, chem, u10, v10, ust, T2, snowh, tsk,   &
         xice, coszen, drydepflx_o3, drydepflx_hobr, drydepflx_brno3,  &
         drydepflx_hbr, drydepflx_br2, emiss_br2,                      &
         ids,ide, jds,jde, kds,kde,                                    &
         ims,ime, jms,jme, kms,kme,                                    &
         its,ite, jts,jte, kts,kte                                     )

  USE module_configure
  USE module_state_description, ONLY: param_first_scalar
  USE module_data_mosaic_asect, ONLY: dcen_sect, numptr_aer, lptr_cl_aer, &
          lptr_na_aer, nsize_aer, ai_phase
  USE module_data_mosaic_constants, ONLY:  pi
  USE module_model_constants, ONLY: mwdry, ep_2, r_v
  USE module_mp_morr_two_moment, ONLY: POLYSVP
  USE module_wrf_error ! For wrf_error_fatal

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: id,                                      &
                         ids,ide, jds,jde, kds,kde,               &
                         ims,ime, jms,jme, kms,kme,               &
                         its,ite, jts,jte, kts,kte

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
  ! ice fraction (0-1), cosine of solar zenith angle (dimensionless)
  REAL, DIMENSION(ims:ime, jms:jme)               ,               &
        INTENT(IN) :: u10, v10, ust, t2, snowh, tsk, xice, coszen

  ! Dry deposition fluxes of O3, HOBr, BrNO3, and HBr (mole/m2/s), calculated in
  ! dry_dep_driver
  REAL,  DIMENSION(ims:ime, jms:jme),                             &
        INTENT(IN) :: drydepflx_o3, drydepflx_hobr,               &
        drydepflx_brno3, drydepflx_hbr, drydepflx_br2

  ! Accumulated br2 emission flux (mole/m2)
  REAL, DIMENSION(ims:ime, jms:jme)               ,               &
        INTENT(INOUT) :: emiss_br2

  ! Trace species mixing ratios (gases in ppm, aerosols in microg/m3)
  REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem),           &
        INTENT(INOUT) :: chem

  !-------- Local variables --------
  CHARACTER (LEN=1000) :: msg
  ! Loop indices i j k, ice landuse category
  INTEGER :: i, j, k, isice_temp = 24
  ! Br2 emissions (mole/m2/s), Conversion factor from emissions to ppm, from
  ! emissions to µg/kg, critical temperature for Br2
  ! emissions & recycling (K)
  REAL :: br2_emi = 0.0, conv, conv2, tskcrit = 273.15

  !---- Local variables for the blowing snow model
  INTEGER, PARAMETER :: naersizes = 8
  ! Br emissions (kg/m2/s), NaCL emissions (kg/m2/s), Snow sublimation flux (kg/m2/s), Normalized
  ! column-integrated snow sublimation flux (kg/m2/day), Snow mixing ratio in
  ! saltation layer (kg/kg), Snow mixing ratio at reference wind speed (kg/kg)
  REAL :: br_emi, nacl_emi, qsnow_subl, qsnow_subl_prime, qsnow_saltation, qsnow_0
  ! 10m-windspeed (m/s), critical 10m-windspeed (m/s), temperature-moisture
  ! parameter from Dery and Yau 2001 (-1E-12 m2/s-1), ice saturation mixing
  ! ratio (kg/kg), conductivity term in the sublimation process (m s kg−1),
  ! diffusion term in the sublimation process (m s kg−1), saturation water vapor
  ! pressure wrt ice (Pa)
  REAL :: w10, w10_crit, ksi, qvapor_sat_ice, fconduction, fdiffusion, eis
  ! Aerosol number concentration change (#/m2/s), Conversion factor from
  ! number emissions to #/kg
  REAL :: aer_num_change, conv3
  ! Molecular diffusion coefficient of H2O vapor (m2 s-1), latent heat of
  ! sublimation (J kg-1), coefficient of thermal conductivity of dry air (J m-1
  ! s-1 K-1),
  REAL :: d_v, l_s, k_d
  ! MOSAIC aerosol type, size, phase indexes, and "pointers" to number, Na and
  ! Cl in chem
  INTEGER :: itype, isize, iphase, l_num, l_na, l_cl
  ! Br/NaCl mass ratio in sea salt (kg/kg)
  REAL :: br_nacl_massratio = 0.00233
  ! Average Antarctic snow salinity (psu) - from Massom et al. (2001)
  REAL :: snow_salinity = 8.5
  ! Conversion factor (from kg/m2/day to kg/m2/s)
  REAL :: kconv_day_s = 1.1574E-5
  ! Adjustment factor for the emission flux, based on snow age (unitless),
  ! snow age (hours)
  ! TODO add a snow age variable, starts at 3 days, resets to 0 when there is
  ! more than XX mm of snow per day (5 mm ?)
  ! TODO update salinity and snow age to follow Frey
  REAL ::  a_prime, snow_age = 3 * 24
  ! Reference 10m-wind speed (m/s)
  REAL :: w10_0 = 6.975
  ! Parameters for the qsnow_subl_prime parameterization, from Dery and
  ! Yau (2001)
  REAL :: a0 = 3.78407E-1, a1 = -8.64089E-2, a2 = -1.60570E-2, &
          a3 = 7.25516E-4, a4 = -1.25650E-1, a5 = 2.48430E-2, &
          a6 = -9.56871E-4, a7 = 1.24600E-2, a8 = 1.56862E-3, &
          a9 = -2.93002E-4
  ! Conversion factor for the ksi parameter (unitless)
  REAL :: convert_ksi = -1.0E12
  ! Density of ice (kg/m3) and sea salt (kg/m3)
  REAL :: rho_ice = 914.0, rho_nacl = 2160.
  ! Molar weight of bromine (g/mole)
  REAL :: mw_br2 = 159.81
  ! Mass proportion of total NaCl emissions going into each MOSAIC size bin,
  ! calculated offline by integrating the size distribution in Yang et al.
  ! (2008) TODO check naersizes
  REAL :: mosaic_massfrac_nacl(naersizes)
  ! Bromine depletion factor in sea salt (unitless), from Sander et al. [2003],
  ! interpolated to the MOSAIC center size bin
  REAL :: depletion_factor(naersizes)

  LOGICAL :: do_snow_emissions, do_blowing_snow

  !---- Data declaration
  mosaic_massfrac_nacl = (/ 0.2189E-09, 0.6873E-08, 0.2122E-06, 0.6306E-05, &
                            0.1742E-03, 0.4182E-02, 0.7729E-01, 0.9183E+00 /)
  depletion_factor = (/ -3.8200, -3.8200, -3.6856, -2.2739, &
                         0.0965, 0.3379, 0.2332, 0.1393  /)
  do_snow_emissions = .true.
  do_blowing_snow = .true.

  !-------- Add emissions from snow, recycling and blowing snow --------
  ! From Toyoyota et al. (2011) and Yang et al. (2008)

  DO j = jts,jte
    DO i = its,ite
      br2_emi = 0.0
      nacl_emi = 0.0
      aer_num_change = 0.0

      a_prime = 1.0 / (1.038 + 0.03758 * snow_age - 0.00014349 * snow_age ** 2.0 &
                + 1.911315E-7 * snow_age ** 3.0)

      !-------- Halogen emissions from snow and recycling over snow
      ! Only valid over snow or over sea ice (assuming that ice is always snow
      ! covered below tskcrit, this may not be true for land ice in the
      ! Antarctic)
      ! The snowh threshold needs to be consistent with the one used in the dry
      ! deposition scheme to trigger snow-covered conditions. The xice threshold
      ! used in the whole routine also needs to be consistent with the one used
      ! to set ice-covered conditions in dry deposition (that is, the threshold
      ! used to set ivgtyp = is_ice)
      ! xice(i, j) .GT. 0.9 to avoid cases with high deposition because of
      ! fractional sea ice
      IF ((snowh(i, j) .GT. 0.01 .AND. xice(i, j) .LE. 0.01) .OR. xice(i, j) .GT. 0.9) THEN

        ! Br2 emissions occur only when the surface temperature is below Tcrit
        IF (tsk(i, j) .LE. tskcrit) THEN

          ! Br2 deposition is reemitted directly over snow and ice (meaning it
          ! is not deposited)
          br2_emi = br2_emi + drydepflx_br2(i,j)

          !----- Br2 "trigger" release

          IF (do_snow_emissions) THEN
            ! O3 deposition over first-year sea ice is converted to Br2 with a
            ! molar yield of 0.075 for SZA <= 85°
            !TODO add ice age, for now assume that all ice is first year ice
            !TODO tuned up the factors by x10
            IF (coszen(i,j) .GE. 0.087 .AND. xice(i,j) .GT. 0.9) THEN
              br2_emi = br2_emi + drydepflx_o3(i,j) * 0.075 * 10
            ELSEIF (coszen(i,j) .LT. 0.087 .AND. xice(i,j) .GT. 0.9 ) THEN
              br2_emi = br2_emi + drydepflx_o3(i,j) * 0.001 * 0.01 * 10

            END IF ! coszen >= 0.087 <=> SZA <= 85°, and over ice
          END IF ! do_snow_emissions

          !---- Br2 "recycling" release

          ! Over FY ice, HOBr and BrONO2 deposition gets converted to Br2
          !TODO add ice age, for now assume that all ice is first year ice
          IF (xice(i,j) .GT. 0.9) THEN
            br2_emi = br2_emi + (drydepflx_hobr(i,j) + drydepflx_brno3(i,j))
          END IF

          !TODO Over MY ice, 2 cases. Need to add ice age to do this.
          !   1 - if HOBr+BrONO2 deposition <= HBr deposition, it gets converted to
          !   Br2
          !   2 - if HOBr+BrONO2 deposition > HBr deposition, an amount of [HBr]
          !   gets converted to Br2 and the excess part gets converted to BrCl
          !   (or 50% Br2 if there is no chlorine in the mechanism)

          ! Over continental snow, the smaller of both fluxes gets converted to
          ! Br2
          !TODO Make sure that this is not ocean/seaice. Technically, xice <
          !0.01 should be enough to make sure that it is not a sea ice cell, and
          !there is no snow height over the open ocean, but it would be safer to
          !also use xland
          IF (xice(i,j) .LE. 0.01 .AND. snowh(i, j) .GT. 0.01) THEN
            IF ((drydepflx_hobr(i,j) + drydepflx_brno3(i,j)) .LE. drydepflx_hbr(i,j)) THEN
              br2_emi = br2_emi + (drydepflx_hobr(i,j) + drydepflx_brno3(i,j))
            ELSE
              br2_emi = br2_emi + drydepflx_hbr(i,j)
            END IF
          END IF
        END IF ! tsk <= tskcrit
      END IF ! (snowh(i, j) .GT. 0.01 .OR. xice(i, j) .GT. 0.9), snow or sea ice cell

      !-------- Halogen emissions from blowing snow

      IF (do_blowing_snow) THEN
        ! Only over sea-ice, and below freezing
        IF (xice(i, j) .GT. 0.01 .AND. tsk(i,j) .LT. tskcrit) THEN
          ! Calculate the critical 10m wind speed for blowing snow, from Yang et
          ! al. (2008)
          !TODO It would be better to use T10
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
            !TODO p_phy should be p10m or psfc, t2 should be 10m
            !TODO cross-check the l_s, k_d and d_v parameterizations and output
            !values
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

            ! Sea salt emission flux (kg/m2/s), from Yang et al., (2008), using
            ! the mean Antarctic snow salinity from Massom et al. (2001)
            nacl_emi = qsnow_subl * snow_salinity / 1000.
            IF (nacl_emi < 0.) THEN
            !TODO understand why this is sometimes < 0
            !  write(msg,*) 'module_airsnow_emissions: nacl_emi < 0., nacl_emi=',nacl_emi,' at i,j=',i,j
            !  CALL wrf_error_fatal( trim(msg) )
              nacl_emi = 0.
            ENDIF

            !---- Br/Br2 emissions from Yang et al., (2008)
            ! (kg/m2/s)
            DO isize = 1,nsize_aer(itype)
              br_emi = nacl_emi * mosaic_massfrac_nacl(isize) * br_nacl_massratio * depletion_factor(isize)
              ! Divide by 2 to conserve bromine, and scale by the ice fraction,
              ! convert from kg/m2/s to mole/m2/s
              br2_emi = br2_emi + br_emi / 2. * xice(i,j) * 1.0E3 / mw_br2
            END DO ! isize

          END IF ! w10 .gt. w10_crit
        END IF ! xice(i, j) .GT. 0.01 .AND. tsk(i,j) .LT. tskcrit
      END IF ! do_blowing_snow


      !-------- Add Br2 emissions from airsnow and blowing snow, and sea salt emissions from blowing snow  --------
      !TODO Add p1st checks and variables everywhere in the module
      itype = 1
      iphase = ai_phase

      ! Conversion factor from mole/m2/s to ppm
      conv = 1.0E3 * mwdry / rho_phy(i,kts,j) * dtstep / (dz8w(i,kts,j))
      ! Conversion from kg/m2/s to µg/m3
      conv2 = 1.0E9 / (dz8w(i,kts,j)) * dtstep
      ! Conversion from #/m2/s to #/kg
      conv3 = dtstep / (rho_phy(i,kts,j) * (dz8w(i,kts,j)))

      DO isize = 1, nsize_aer(itype)
        !---- Na emissions
        l_na = lptr_na_aer(isize,itype,iphase)
        IF (l_na >= param_first_scalar) THEN
          chem(i,kts,j,l_na) = chem(i,kts,j,l_na) + nacl_emi / 2.0 * mosaic_massfrac_nacl(isize) * conv2
        END IF
        !---- Cl emissions
        l_cl = lptr_cl_aer(isize,itype,iphase)
        IF (l_cl >= param_first_scalar) THEN
          chem(i,kts,j,l_cl) = chem(i,kts,j,l_cl) + nacl_emi / 2.0 * mosaic_massfrac_nacl(isize) * conv2
        END IF
        !---- Aerosol number "emissions" from the mass change
        l_num = numptr_aer(isize,itype,iphase)
        IF (l_num >= param_first_scalar) THEN
          ! Calculate the change in aerosol number concentrations due to
          ! the sea salt emissions (#/m2/s). dcen_sect is in cm in MOSAIC
          aer_num_change = nacl_emi / rho_nacl * 6.0/pi * 1.0/(dcen_sect(isize,itype)*1E-2)**3.0
          ! #/m2/s -> #/kg
          chem(i,kts,j,l_num) = chem(i,kts,j,l_num) + aer_num_change * conv3
        END IF
        ! Add Br2 emissions
        IF (chem(i,kts,j,p_br2)  + br2_emi * conv <= 0) THEN
          br2_emi = - chem(i,kts,j,p_br2) / conv
          chem(i,kts,j,p_br2) = 0.0
        ELSE
          chem(i,kts,j,p_br2) = chem(i,kts,j,p_br2) + br2_emi * conv
        END IF

        !---- Output diagnosis, accumulated br2 emissions in mole/m2
        ! TODO separate recycling, blowing snow, snow trigger, br2 reemissions
        ! TODO remove br2 reemissions and replace by Br2 deposition set to 0
        ! over snow/ice ?
        emiss_br2(i,j) = emiss_br2(i,j) + br2_emi * dtstep

      END DO !isize

    END DO ! i loop
  END DO ! j loop

END SUBROUTINE airsnow_emissions

END MODULE module_airsnow_emissions

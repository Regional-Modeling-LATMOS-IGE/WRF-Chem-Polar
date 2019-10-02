! 2018/12/06, Louis Marelle
!
! Purpose:
! Compute Br2 emissions from the snowpack and add them to chem(:, :, :, p_br2)
!
! References:
! Toyota, K., McConnell, J. C., Lupu, A., Neary, L., McLinden, C. A., Richter,
! A., Kwok, R., Semeniuk, K., Kaminski, J. W., Gong, S.-L., Jarosz, J.,
! Chipperfield, M. P., and Sioris, C. E.: Analysis of reactive bromine
! production and ozone depletion in the Arctic boundary layer using 3-D
! simulations with GEM-AQ: inference from synoptic-scale patterns, Atmos. Chem.
! Phys., 11, 3949-3979, https://doi.org/10.5194/acp-11-3949-2011, 2011.
! and
! Falk, S. and Sinnhuber, B.-M.: Polar boundary layer bromine explosion and
! ozone depletion events in the chemistry–climate model EMAC v2.52:
! implementation and evaluation of AirSnow algorithm, Geosci. Model Dev., 11,
! 1115-1131, https://doi.org/10.5194/gmd-11-1115-2018, 2018.
!
! Remarks:
! For now, this is only coupled to the SAPRC99_MOSAIC_8BIN_VBS2_AQ_HALOGENS_KPP
! gas-phase/aerosol mechanism 
!
!TODO change the critical temperature to OK but include a memory to this,
!everytime the temperature rises above 0K or if liquid rain falls bromine
!emissions and recycling are turned off, until new snow falls again and
!temperatures are below 0K.
!TODO 
! - If someone tries to use this option with another mechanism include a call
! to WRF_FATAL in chemics_init.F.
! TODO This is hardcoded to USGS for now read in IVGTYP and MMINLU, abort if
! MMINLU is not USGS. Also, determine if land type is sea ice by reading xice .gt. 0
! and ivgtyp .eq. isice, determine if land snow by xice .eq. 0 and ivgtyp .ne.
! isocean, also treat land ice by xice .neq. 0 and ivgtyp .ne. isocean. Write a
! caveat for lake ice, which should not emit Br2.

MODULE module_airsnow_emissions

CONTAINS

!----------------------------------------------------------------------
  SUBROUTINE  airsnow_emissions( id, dtstep, dz8w, rho_phy,         &
         chem, snowh, snowc, tsk, xice, coszen, drydepflx_o3,       &
         drydepflx_hobr, drydepflx_brno3, drydepflx_hbr,            &
         drydepflx_br2, emiss_br2,                                  &
         ids,ide, jds,jde, kds,kde,                                 &
         ims,ime, jms,jme, kms,kme,                                 &
         its,ite, jts,jte, kts,kte                                  )

  USE module_configure
  USE module_state_description
  USE module_model_constants, only: mwdry

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: id,                                      &
                         ids,ide, jds,jde, kds,kde,               &
                         ims,ime, jms,jme, kms,kme,               &
                         its,ite, jts,jte, kts,kte

  ! Dynamical time step (s)
  REAL, INTENT(IN) :: dtstep

  ! Vertical layer thickness (m), air density (km/m3)
  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),                     &
        INTENT(IN) :: dz8w, rho_phy

  ! Snow height (m), snow cover (0-1), skin temperature (K), ice fraction (0-1),
  ! cosine of solar zenith angle (dimensionless)
  !TODO: snow cover is not used for now, remove ?
  REAL, DIMENSION(ims:ime, jms:jme)               ,               &
        INTENT(IN) :: snowh, snowc, tsk, xice, coszen

  ! Dry deposition fluxes of O3, HOBr, BrNO3, and HBr (mol/m2/s), calculated in
  ! dry_dep_driver
  REAL,  DIMENSION(ims:ime, jms:jme),                             &
        INTENT(IN) :: drydepflx_o3, drydepflx_hobr,               &
        drydepflx_brno3, drydepflx_hbr, drydepflx_br2

  ! Accumulated br2 emission flux (mol/m2)
  REAL, DIMENSION(ims:ime, jms:jme)               ,               &
        INTENT(INOUT) :: emiss_br2

  ! Trace species mixing ratios (gases in ppm, aerosols in microg/m3)
  REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem),           &
        INTENT(INOUT) :: chem

  ! Local variables
  ! Loop indices i j k, ice landuse category
  INTEGER :: i, j, k, isice_temp = 24
  ! Br2 emissions (moles/m2/s), Conversion factor from emissions to ppm,
  ! critical temperature for Br2 emissions & recycling (K)
  REAL :: br2_emi = 0.0, conv, tskcrit = 273.15


!-------- Add emissions --------
! From Toyoyota et al. (2011)

  DO j = jts,jte  
    DO i = its,ite 

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
      IF (snowh(i, j) .GT. 0.01 .OR. xice(i, j) .GT. 0.9) THEN
        br2_emi = 0.0

        ! Br2 emissions occur only when the surface temperature is below Tcrit (-15°?). 
        !TODO Toyota et al. (2011) say "surface temperature", is this TSK or T2 ? 
        IF (tsk(i, j) .LE. tskcrit) THEN

          ! Br2 deposition is reemitted directly over snow and ice (meaning it
          ! is not deposited)
          br2_emi = br2_emi + drydepflx_br2(i,j)

          !----- Br2 "trigger" release

          ! O3 deposition over first-year sea ice is converted to Br2 with a
          ! molar yield of 0.075 for SZA <= 85°
          !TODO add ice age, for now assume that all ice is first year ice
          !TODO tuned up the factors by x10
          IF (coszen(i,j) .GE. 0.087 .AND. xice(i,j) .GT. 0.9) THEN
            br2_emi = br2_emi + drydepflx_o3(i,j) * 0.075 * 10
          ELSEIF (coszen(i,j) .LT. 0.087 .AND. xice(i,j) .GT. 0.9 ) THEN
            br2_emi = br2_emi + drydepflx_o3(i,j) * 0.001 * 0.01 * 10
          END IF ! coszen >= 0.087 <=> SZA <= 85°, and over ice

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
          !TODO make sure that this is not ocean/seaice. Technically, xice < 0.9
          !should be enough to make sure that it is not a sea ice cell, and
          !there is no snow height over the open ocean, but it would be safer to
          !also use xland
          IF (xice(i,j) .LE. 0.9 .AND. snowh(i, j) .GT. 0.01) THEN
            IF ((drydepflx_hobr(i,j) + drydepflx_brno3(i,j)) .LE. drydepflx_hbr(i,j)) THEN
              br2_emi = br2_emi + (drydepflx_hobr(i,j) + drydepflx_brno3(i,j))
            ELSE
              br2_emi = br2_emi + drydepflx_hbr(i,j)
            END IF
          END IF
          
          ! Conversion factor from mol/m2/s to ppm
          conv = 1.0E3 * mwdry / rho_phy(i,kts,j) * dtstep / (dz8w(i,kts,j))

          ! Add emissions to chem at the surface (k == kts)
          ! (i.e., emissions tendencies over time dtstep are applied 
          ! to the surface trace gas concentrations)
          chem(i, kts, j, p_br2) = chem(i, kts, j, p_br2) + br2_emi * conv

          ! Output diagnosis, accumulated br2 emissions in mol/m2
          emiss_br2(i, j) = emiss_br2(i, j) + br2_emi * dtstep

        END IF ! tsk <= tskcrit
      END IF ! (snowh(i, j) .GT. 0.01 .OR. xice(i, j) .GT. 0.9), snow or sea ice cell

    END DO ! i loop
  END DO ! j loop

END SUBROUTINE airsnow_emissions

END MODULE module_airsnow_emissions

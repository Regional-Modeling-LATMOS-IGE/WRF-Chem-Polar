! 2025/03 Louis Marelle, Rémy lapere
!
! Purpose:
!  Compute DMS emissions
!
! Remarks:
!  - This requires that DMS_OCEAN is present in the wrflowinp file and that
!    the DMS tracer is defined in chem for this mechanism
!
! References:
! ==============================================================================
! Air-sea exchange functions:
!  (1) Nightingale et al., 2000: https://doi.org/10.1029/1999GB900091
! Schmidt number for DMS
!  (1) Saltzman et al., 1993: https://doi.org/10.1029/93JC01858
! Sea ice correction factor
!  (1) Loose et al., 2009: doi:10.1029/2008GL036318
! Implementation and evaluation:
!  (1) Marelle et al., 2017: doi:10.5194/gmd-10-3661-2017
!  (2) Lapere et al., 2025: in prep.

MODULE module_dms_emis

  IMPLICIT NONE

  CONTAINS

!------------------------------------------------------------------------------

  SUBROUTINE dms_emis( id, dtstep, dz8w, dms_opt,                   &
         rho_phy, chem,                            &
         u10, v10, xland, tsk, xice, lakemask, dms_ocean,           &
         ids,ide, jds,jde, kds,kde,                                 &
         ims,ime, jms,jme, kms,kme,                                 &
         its,ite, jts,jte, kts,kte                                  )

    ! Purpose:
    !  Compute DMS emissions

    ! Arguments:
    ! Input:
    !  id = domain ID number
    !  dtstep = model time step (s)
    !  dz8w = 3D vertical level depth (m)
    !  dms_opt = dms emission option from config_flags
    !  rho_phy = 3D air density (kg/m3)
    !  u10 = windpeed at 10m in x (positive west-east) direction (m/s)
    !  v10 = windpeed at 10m in y (positive south-north) direction (m/s)
    !  xland = landmask flag (1 land, 2 water)
    !  tsk = surface skin temperature (K)
    !  xice = sea ice or lake ice fraction (0-1)
    !  lakemask = lake mask (1 for lakes, 0 for rest)
    !  dms_ocean = oceanic concentration of DMS in micro mol/m3
    !  ids,ims,its etc. = domain, memory, and tile start and end indices
    ! Input/output
    !  chem = advected chemical species (gases in ppm, aerosols in microg/kg)

    USE module_state_description

    IMPLICIT NONE

    !---- Arguments
    ! Input
    INTEGER, INTENT(IN) :: dms_opt
    INTEGER, INTENT(IN) :: id,                                      &
                           ids,ide, jds,jde, kds,kde,               &
                           ims,ime, jms,jme, kms,kme,               &
                           its,ite, jts,jte, kts,kte
    REAL, INTENT(IN) :: dtstep
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: dz8w, &
                                                              rho_phy
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: u10, v10, xland, &
                                          tsk, xice, lakemask, dms_ocean
    ! Input/output
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), &
        INTENT(INOUT) :: chem

    !---- Local variables
    ! DMS emission flux (mol km-2 hr-1)
    REAL :: dms_emiss
    ! Surface oceanic dms content in mol/m^3
    REAL :: dms_ocean_sfc
    ! 10-m wind (m/s), sea surface temperature (celsius)
    REAL :: w10, sst
    ! Schmidt number for DMS
    REAL :: sc_dms
    ! Sea-air exchange coefficient for DMS (m/hour)
    REAL :: k_dms
    ! Conversion factor from mol km-2 hr-1 to ppmv
    REAL :: conv
    ! Sea ice fraction
    REAL :: icefrac
    ! Indices
    INTEGER :: i, j

    !---- Parameters
    ! Schmidt number for co2 at 20°C = 600 (reference value)
    REAL, PARAMETER :: sc_co2=600.

    !-------- Compute DMS emissions --------

    DO i=its,ite
      DO j=jts,jte

        ! Only perform DMS emissions over open water (xland>1.5),
        ! Sea ice is considered land, so include fractional sea ice points
        ! explicitly. We cannot use xice(i,j) .LT. 0.0 on its own because this
        ! is true for all land points.
        icefrac = xice(i,j)
        IF ( (xland(i,j) .GT. 1.5 &
            .OR. (icefrac .GT. 1.E-9 .AND. icefrac .LT. 1.0)) ) THEN
!Could exclude lakes            .AND. lakemask(i,j) .LT. 0.5) THEN

          !-------- Initialize DMS emission parameters --------
          dms_emiss = 0.0
          ! DMS surface concentration in mol/m3
          dms_ocean_sfc = dms_ocean(i,j)*1.E-6

          IF(dms_ocean_sfc > 0.) THEN
            ! 10-m wind speed (m/s)
            w10 = SQRT(u10(i, j)**2.0 + v10(i, j)**2.0)
            ! Conversion factor from mol km-2 hr-1 to ppmv
            conv = 4.828e-4 / rho_phy(i,kts,j) * dtstep / (dz8w(i,kts,j) * 60.)

            !---- DMS emission flux from Nightingale et al. (2000)
            IF(dms_opt == DMSNIGHTINGALE) THEN
              ! SST - Schmidt number parameterization valid 5 to 30°C
              sst = ( tsk(i,j) - 273.15 )
              sst = MAX(MIN(sst, 30.0), 5.0)
              ! Schmidt number for DMS (Saltzman et al., 1993)
              sc_dms = 2674.0 - 147.12 * sst + 3.726 * sst ** 2.0 &
                          - 0.038 * sst ** 3
              ! Sea - air exchange coefficient from Nightingale et al. 2000 (m/hour)
              ! k_dms = ( 0.222*w10**2 + 0.333*w10) * (sc_dms/sc_co2)**(-0.5) / 3600. / 100.
              k_dms = (0.00222*w10**2 + 0.00333*w10) * (sc_dms/sc_co2)**(-0.5)
              ! DMS emission flux in mol m-2 hour-1
              dms_emiss = k_dms * dms_ocean_sfc
              ! Convert to  WRF-Chem units (mol km-2 hour-1)
              dms_emiss = dms_emiss * 1.0E6

            ELSE
              !---- Other dms_opt cases are not supported
              CALL wrf_error_fatal ('dms_emis: dms_opt not supported')

            ENDIF ! dms_opt

            ! Account for fractional sea ice (Loose et al., 2009)
            IF(icefrac > 0.01) THEN
              dms_emiss = dms_emiss*(1.0-icefrac)**0.4
            ENDIF
            ! Add emissions to chem, at the surface
            IF(p_dms > param_first_scalar) THEN
              chem(i,kts,j,p_dms) = chem(i,kts,j,p_dms) + dms_emiss*conv
            ENDIF

          ENDIF ! dms_ocean_sfc
        END IF ! xland
      END DO ! j
    END DO ! i

  END SUBROUTINE dms_emis

END MODULE module_dms_emis

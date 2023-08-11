! Louis Marelle, Lefteris Ioannidis, 2023/06/12
!
! Purpose:
!  Compute size-resolved sea spray aerosol emissions (mass and number) for the
!  sectional 4-bin or 8-bin MOSAIC aerosols, or for the GOCART 4 sea salt bins.
!
! Remarks:
!  - For now this only does sea-spray NaCl and sea-spray sulfate emissions
!  - TODO include a namelist option for enabling ss_so4
!  - TODO Include ssorganics as size-resolved fraction of total sea spray
!    (ssa_oa_frac, ssa_nacl_frac)

! References:
! ==============================================================================
! Source functions:
! (1) Grythe et al., 2014: doi:10.5194/acp-14-1277-2014
! (2) Monahan et al., 1986: doi:10.1007/978-94-009-4668-2_16
!TODO replace with 1997 (3) Gong, 2003: doi:10.1029:2003GB002079
! (5) Ioannidis et al., 2023: TODO doi
!TODO O'Dowd 1997
! Correction factors:
! (6) Jaegle et al., 2011: doi:10.5194/acp-11-3137-2011
! (7) Salisbury et al., 2013: doi:10.1002/2013JC008797
! SSSO4 fraction:
! (8) Calhoun et al., 1991, doi: https://doi.org/10.1029/91GL02304
!

MODULE module_seaspray_emis

  IMPLICIT NONE

  !---- Global module-wide variables
  ! These are initialized in seaspray_emis_init

  ! Allocatable size-resolved aerosol arrays (dim is size bin number and depends
  ! on aerosol mechanism, unit for vol_emiss_dist is um3/m2/s,
  ! bin_diameters_cen in um)
  REAL, DIMENSION(:), ALLOCATABLE :: vol_emiss_dist, &
    bin_diameters_cen

  ! Number of sea spray aerosol bins for chem_opt in use
  INTEGER :: naersizes

  ! chem_opt aerosol model family
  LOGICAL :: aer_is_gocart
  LOGICAL :: aer_is_mosaic4bin
  LOGICAL :: aer_is_mosaic8bin

  CONTAINS

!------------------------------------------------------------------------------

 SUBROUTINE seaspray_emis( id, dtstep, dz8w, chem_opt, seas_opt,     &
         rho_phy, chem, e_nacl_ocean,                               &
         u10, v10, xland, tsk, xice, lakemask,                      &
         ids,ide, jds,jde, kds,kde,                                 &
         ims,ime, jms,jme, kms,kme,                                 &
         its,ite, jts,jte, kts,kte                                  )

    ! Purpose:
    !  Compute size-resolved sea spray aerosol number and mass emissions

    ! Arguments:
    ! Input:
    !  id = domain ID number
    !  dtstep = model time step (s)
    !  dz8w = 3D vertical level depth (m)
    !  chem_opt = chemical mechanism + aerosol mechanism option from config_flags
    !  seas_opt = sea spray emission option from config_flags
    !  rho_phy = 3D air density (kg/m3)
    !  u10 = windpeed at 10m in x (positive west-east) direction (m/s)
    !  v10 = windpeed at 10m in y (positive south-north) direction (m/s)
    !  xland = landmask flag (1 land, 2 water)
    !  tsk = surface skin temperature (K)
    !  xice = sea ice or lake ice fraction (0-1)
    !  lakemask = lake mask (1 for lakes, 0 for rest)
    !  ids,ims,its etc. = domain, memory, and tile start and end indices
    ! Input/output
    !  chem = advected chemical species (gases in ppm, aerosols in microg/kg)
    !  e_nacl_ocean = accumulated NaCl emission flux from sea spray (microg/m2) !TODO should be mg/m2

    !TODO include all explicit ONLY: associations for module_state_description
    ! chem_opt, seas_opt, param_first_scalar, num_chem, p_seas_1, p_seas_2, p_seas_3, p_seas_4
    USE module_state_description
    USE module_configure, ONLY : grid_config_rec_type
    !TODO Is it possible to only include this for mosaic runs?
    USE module_data_mosaic_asect, ONLY: numptr_aer, lptr_cl_aer, &
          lptr_na_aer, lptr_so4_aer, nsize_aer

    IMPLICIT NONE

    !---- Dummy arguments
    ! Input
    INTEGER, INTENT(IN) :: chem_opt, seas_opt
    INTEGER, INTENT(IN) :: id,                                      &
                           ids,ide, jds,jde, kds,kde,               &
                           ims,ime, jms,jme, kms,kme,               &
                           its,ite, jts,jte, kts,kte
    REAL, INTENT(IN) :: dtstep
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: dz8w, &
                                                              rho_phy
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: u10, v10, xland, &
                                                     tsk, xice, lakemask
    ! Input/output
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), &
        INTENT(INOUT) :: chem
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: e_nacl_ocean

    !---- Local variables
    ! Sea spray emission fluxes, microm3/m2/s, #/m2/s
    REAL :: seaspray_vol_emiss, seaspray_num_emiss
    ! Sea salt, sulfate, organic emission fluxes, microg/m2/s
    REAL :: nacl_mass_emiss, so4_mass_emiss, org_mass_emiss
    ! 10-m wind (m/s), sea surface temperature (celsius)
    REAL :: w10, sst
    ! Open ocean fraction in grid cell (0-1)
    REAL :: open_ocean_frac
    ! Whitecap fraction in grid cell (0-1)
    REAL :: whitecap_frac
    ! SST scaling factor (unitless), SST constrained to valid range (celsius)
    REAL :: sst_dependence, sst_valid
    ! aerosol bin center diameter (micro m)
    REAL :: dcen
    ! Conversion factor from 1/m2/s to 1/kg_dry_air
    REAL :: conv
    ! Conversion factor from sea spray volume to nacl mass (ug(nacl)/um3(seaspray))
    REAL :: nacl_mass_to_ss_vol_ratio
    ! Conversion factor from sea spray volume to sulfate mass (ug(so4)/um3(seaspray))
    REAL :: so4_mass_to_ss_vol_ratio
    ! Sea-salt sulfate / nacl mass ratio (kg/kg)
    ! Calhoun et al. (1991), https://doi.org/10.1029/91GL02304, converted from kg(so4)/kg(na) to kg(so4)/kg(nacl)
    !TODO this should probably be size-resolved
    REAL :: so4_nacl_ratio 
    ! Indices
    INTEGER :: i, j, isize, iphase
    INTEGER :: iaerbin
    INTEGER :: l_na, l_cl, l_so4, l_num, l_seasalt

    !---- Parameters
    !TODO Read from model constants instead
    REAL, PARAMETER :: PI = 3.14159
    !TODO should be in init and/or data module if you want to avoid hardcoded
    !values
    ! REAL, PARAMETER :: NACL_AERDEN = 2.17 * 1.0E6 * 1.0E-12 ! g/cm3 -> microg/cm3 -> microg/um3
    ! REAL, PARAMETER :: SO4_AERDEN = 1.8 * 1.0E6 * 1.0E-12 ! g/cm3 -> microg/cm3 -> microg/um3, MOSAIC value
    REAL, PARAMETER :: NACL_AERDEN = 2.17E-6 ! microg/um3
    REAL, PARAMETER :: SO4_AERDEN = 1.8E-6 ! microg/um3, MOSAIC value
    REAL, PARAMETER :: MW_NA = 22.99, MW_CL = 35.45, MW_SO4 = 96.06
    ! Logicals
    LOGICAL, PARAMETER :: DO_SEAS_SO4 = .true.

    !-------- Initialize the mass fractions of so4 and nacl in sea spray --------
    ! For convenience, we integrate the so4/nacl/oa fractioning into the
    ! conversion factor from sea spray aerosol volume to sea spray aerosol
    ! mass.
    !TODO the following is no longer valid if marine OA emissions are added
    !TODO should be in init module
    IF(DO_SEAS_SO4) THEN
      so4_nacl_ratio = 0.252 * MW_NA / (MW_NA + MW_CL)
      nacl_mass_to_ss_vol_ratio = 1. / (so4_nacl_ratio/(1.-so4_nacl_ratio)/SO4_AERDEN + 1./NACL_AERDEN)
      so4_mass_to_ss_vol_ratio = 1. / ((1.-so4_nacl_ratio)/so4_nacl_ratio/NACL_AERDEN + 1./SO4_AERDEN)
    ELSE
      so4_nacl_ratio = 0.0
      nacl_mass_to_ss_vol_ratio = NACL_AERDEN
      so4_mass_to_ss_vol_ratio = 0.0
    !TODO add case for OA emissions, and oa_mass_to_ss_vol_ratio. Careful that
    !this will also change the nacl and so4 ratios, who will need to be
    !redefined above
    ENDIF

    !-------- Compute SSA emissions --------
    iphase = 1 ! No emissions to cloud-phase !TODO aiphase, check
    DO i=its,ite
      DO j=jts,jte

        ! Only perform sea spray emissions over open sea water (xland>1.5,
        ! lakemask(i,j)<0.5)
        ! Sea ice is considered land, so include fractional sea ice points
        ! explicitly. We cannot use xice(i,j) .LT. 1.0 on its own because this
        ! is true for all land points.
        IF ( (xland(i,j) .GT. 1.5 &
            .OR. (xice(i,j) .GT. 0.01 .AND. xice(i,j) .LT. 1.0)) &
            .AND. lakemask(i,j) .LT. 0.5) THEN

          !-------- Initialize sea spray emission parameters --------
          nacl_mass_emiss = 0.0
          so4_mass_emiss = 0.0
          whitecap_frac = 0.0
          sst_dependence = 1.0
          ! 10-m wind speed (m/s)
          w10 = SQRT(u10(i, j)**2.0 + v10(i, j)**2.0)
          ! Conversion factor from 1/m2/s to 1/kg_dry_air
          conv = 1.0/(rho_phy(i,kts,j) * dz8w(i,kts,j)) * dtstep
          !---- Open ocean fraction
          ! When fractional sea ice is present, scale the emission flux by the
          ! open ocean fraction
          open_ocean_frac = 1.0 - xice(i,j)
          !---- Whitecap fraction, used later in sea spray source functions
          IF(seas_opt == SEASIOANNIDIS) THEN
            ! Whitecap fraction from Salisbury et al. (2013), recommended by
            ! Ioannidis et al. (2023)
            whitecap_frac = 4.60E-5 * w10**2.26
          ELSE
            ! Whitecap fraction from Monahan et al. (1986)
            whitecap_frac = 3.84E-6 * w10**3.41
          ENDIF ! seas_opt
          !---- SST correction factors
          IF(seas_opt == SEASMONAHAN .OR. seas_opt == SEASIOANNIDIS) THEN
            ! SST correction from Jaegle et al. (2011), valid between -2 C and
            ! 30 C, recommended by Ioannidis et al. (2023)
            sst_valid = MIN(MAX(tsk(i,j)-273.15, -2.0), 30.0)
            sst_dependence = 0.3 + 0.1*sst_valid - 0.0076*sst_valid**2.0 + 0.00021*sst_valid**3.0
          ENDIF ! seas_opt


          !-------- Compute sea-spray aerosol mass/num emissions for each aer size bin --------
          DO iaerbin = 1,naersizes
            !---- Compute sea-spray aerosol  mass and number emission fluxes
            ! Calculate the volume emission flux seaspray_vol_emiss(iaerbin)
            ! (um3/m2/s), from size-resolved volume emission scaling factors
            ! vol_emiss_dist(iaerbin) (unit is source function dependent,
            ! um3/m2/s for most source functions)
            seaspray_vol_emiss = open_ocean_frac * sst_dependence * whitecap_frac &
                          * vol_emiss_dist(iaerbin)
            ! Do not allow negative emissions
            seaspray_vol_emiss = MAX(seaspray_vol_emiss, 0.0)
            !---- Convert sea spray aerosol volume emission flux (microm3/m2/s)
            !  to speciated aerosol mass emission flux (microg/m2/s)
            nacl_mass_emiss = seaspray_vol_emiss * nacl_mass_to_ss_vol_ratio
            IF(DO_SEAS_SO4) THEN
              so4_mass_emiss = seaspray_vol_emiss * so4_mass_to_ss_vol_ratio
            ENDIF
            !---- Convert sea spray aerosol volume emission flux to number
            ! emission flux (#/m2/s), assuming aerosol volume is the volume at
            ! bin log-center diameter dcen
            dcen = bin_diameters_cen(iaerbin)
            seaspray_num_emiss = seaspray_vol_emiss / ((PI/6.0)*dcen**3.0)
            seaspray_num_emiss = MAX(seaspray_num_emiss, 0.0)
            !---- Accumulate sea-spray nacl emissions from open ocean, output diagnosis
            e_nacl_ocean(i, j) = e_nacl_ocean(i, j) + nacl_mass_emiss * dtstep


            !-------- Emit sea spray emissions (add them to chem(i,k,j,l)) --------
            !---- Emit MOSAIC aerosols to chem
            IF(aer_is_mosaic4bin .OR. aer_is_mosaic8bin) THEN
              ! MOSAIC sodium aerosol
              l_na = lptr_na_aer(iaerbin,1,iphase)
              IF (l_na >= param_first_scalar) THEN
                chem(i,kts,j,l_na) = chem(i,kts,j,l_na) + nacl_mass_emiss * MW_NA/(MW_NA+MW_CL) * conv
              END IF
              ! MOSAIC chlorine aerosol
              l_cl = lptr_cl_aer(iaerbin,1,iphase)
              IF (l_cl >= param_first_scalar) THEN
                chem(i,kts,j,l_cl) = chem(i,kts,j,l_cl) + nacl_mass_emiss * MW_CL/(MW_NA+MW_CL) * conv
              END IF
              ! MOSAIC sulfate aerosol
              IF(DO_SEAS_SO4) THEN
                l_so4 = lptr_so4_aer(iaerbin,1,iphase)
                IF (l_cl >= param_first_scalar) THEN
                  chem(i,kts,j,l_so4) = chem(i,kts,j,l_so4) + so4_mass_emiss * conv
                END IF
              ENDIF
              ! MOSAIC aerosol number
              l_num = numptr_aer(iaerbin,1,iphase)
              IF (l_num >= param_first_scalar) THEN
                chem(i,kts,j,l_num) = chem(i,kts,j,l_num) + seaspray_num_emiss * conv
              END IF
            !---- Emit GOCART aerosols to chem
            ELSEIF(aer_is_gocart) THEN
              ! GOCART seasalt aerosol
              IF(iaerbin == 1 ) THEN
                l_seasalt = p_seas_1
              ELSEIF(iaerbin == 2 ) THEN
                l_seasalt = p_seas_2
              ELSEIF(iaerbin == 3 ) THEN
                l_seasalt = p_seas_3
              ELSEIF(iaerbin == 4 ) THEN
                l_seasalt = p_seas_4
              ENDIF
              IF (l_seasalt >= param_first_scalar) THEN
                chem(i,kts,j,l_seasalt) = chem(i,kts,j,l_seasalt) + nacl_mass_emiss * conv
              ENDIF
              !TODO include sulfate emissions for GOCART, only one size class,
              !     either move out of iaerbin loop or do for iaerbin = 1
            !---- Other chem_opt cases are not supported
            ELSE
              CALL wrf_error_fatal ('seaspray_emis: chem_opt not supported')
           ENDIF ! aer_is_


          END DO ! iaerbin

        END IF ! xland
      END DO ! j
    END DO ! i

  END SUBROUTINE seaspray_emis

!------------------------------------------------------------------------------

  SUBROUTINE seaspray_emis_init(chem_opt, seas_opt)

    ! Purpose:
    !  Initialize the parameters for the sea spray aerosol source functions used
    !  in seaspray_emis, which depend on chem_opt and seas_opt
    !
    ! Remarks:
    ! - vol_emiss_dist was precalculated by integrating the source functions
    !   size distributions from the original publications over the WRF-Chem
    !   aerosol size ranges. For this reason if the aerosol bin limits are changed
    !   vol_emiss_dist needs to be recalculated.
    ! - TODO Center diameters and naersizes should be compared to the ones used at
    !   runtime (for MOSAIC, dcen_sect(iaerbin,itype) and nsize_aer(itype) from
    !   module_data_mosaic_asect), and model should crash if there is a mismatch
    !   - TODO If comparing to MOSAIC nsizes and dcen_sect, careful the MOSAIC top
    !     level initializations should have been done before.
    ! - TODO Allocatables here are never deallocated
    ! - TODO some of this could be in a separate data module

    !TODO only use with only:, specify all variables here
    USE module_state_description ! For chem_opt and seas_opt option values

    IMPLICIT NONE

    ! Dummy arguments:
    ! Input
    INTEGER, INTENT(IN) :: chem_opt, seas_opt

    !---- Initialize chem mechanism families IDs
    !TODO Could be automated by checking the presence of mechanism-specific pointers
    !TODO this needs to be updated if a new mechanism is defined
    aer_is_gocart = chem_opt == MOZCART_KPP .OR. &
                    chem_opt == T1_MOZCART_KPP .OR. &
                    chem_opt == GOCART_SIMPLE .OR. &
                    chem_opt == GOCARTRADM2 .OR. &
                    chem_opt == GOCARTRACM_KPP

    aer_is_mosaic4bin = chem_opt == CBMZ_MOSAIC_4BIN .OR. &
                        chem_opt == CBMZ_MOSAIC_DMS_4BIN .OR. &
                        chem_opt == CBMZ_MOSAIC_4BIN_AQ .OR. &
                        chem_opt == CBMZ_MOSAIC_DMS_4BIN_AQ .OR. &
                        chem_opt == MOZART_MOSAIC_4BIN_KPP .OR. &
                        chem_opt == MOZART_MOSAIC_4BIN_AQ_KPP .OR. &
                        chem_opt == CRI_MOSAIC_4BIN_AQ_KPP .OR. &
                        chem_opt == SAPRC99_MOSAIC_4BIN_VBS2_KPP

    aer_is_mosaic8bin = chem_opt == CBMZ_MOSAIC_8BIN .OR. &
                        chem_opt == CBMZ_MOSAIC_DMS_8BIN .OR. &
                        chem_opt == CBMZ_MOSAIC_8BIN_AQ .OR. &
                        chem_opt == CBMZ_MOSAIC_DMS_8BIN_AQ .OR. &
                        chem_opt == CRI_MOSAIC_8BIN_AQ_KPP .OR. &
                        chem_opt == SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP .OR. &
                        chem_opt == SAPRC99_MOSAIC_8BIN_VBS2_KPP

    !-------- Initialize sea spray aerosol bin number for chem_opt --------
    ! This could be recovered from the aerosol scheme data directly, but
    ! vol_emiss_dist is already chem_opt-dependent and depends on naersizes so
    ! there is no harm in redefining it here just for this module
    IF(aer_is_gocart) THEN
      ! GOCART with 4 sea salt bins cases
      naersizes = 4
    ELSEIF(aer_is_mosaic8bin) THEN
      ! MOSAIC-8bin cases
      naersizes = 8
    ELSEIF(aer_is_mosaic4bin) THEN
      ! MOSAIC-4bin cases
      naersizes = 4
    ELSE
      CALL wrf_error_fatal ('seaspray_emis_init: chem_opt not supported')
    ENDIF

    !-------- Initialize size-resolved aerosol arrays for SSA emissions --------
    ! Initialize vol_emiss_dist (microm3/m2/s) and bin_diameters_cen
    ! (micrometer) for each family of aerosol options selected by
    ! chem_opt

    ! Allocate the size-resolved aerosol arrays
    IF (.NOT. ALLOCATED(vol_emiss_dist)) THEN
      ALLOCATE(vol_emiss_dist(naersizes))
    ENDIF
    IF (.NOT. ALLOCATED(bin_diameters_cen)) THEN
      ALLOCATE(bin_diameters_cen(naersizes))
    ENDIF

    !-- Initialize SSA source functions for GOCART aerosols (4 sea salt bins)
    IF(aer_is_gocart) THEN
      ! GOCART sea salt bin center diameters (micrometer)
      ! These should be log-averaged but I follow the GOCART definition instead
      bin_diameters_cen = (/ 0.60, 2.00, 6.5, 15.0 /)
      ! Size-resolved volume emission scaling factors for GOCART
      IF(seas_opt == SEASMONAHAN) THEN
        ! Monahan et al. (1986) SSA source function
        vol_emiss_dist = (/ 367221., 3534527., 2057189., 3732253. /)
      ELSEIF(seas_opt == SEASIOANNIDIS .OR. seas_opt == SEASGONG) THEN
        ! Gong et al. (1997) SSA source function
        ! vol_emiss_dist = (/ 553577., 1701853., 124828., 40680. /)
        ! Same as Monahan above 100 nm
        vol_emiss_dist = (/ 367221., 3534527., 2057189., 3732253. /)
      ELSE
        CALL wrf_error_fatal ('seaspray_emis_init: seas_opt not supported')
      ENDIF

    !-- Initialize SSA source functions for MOSAIC-8bin aerosols
    ELSEIF(aer_is_mosaic8bin) THEN
      ! MOSAIC-8bin bin center diameters (micrometer)
      bin_diameters_cen = (/ 0.0552, 0.1105, 0.2210, 0.4419, &
                             0.8839, 1.7678, 3.5355, 7.0711 /)
      ! Size-resolved volume emission scaling factors for MOSAIC-8bin
      IF(seas_opt == SEASMONAHAN) THEN
        ! Monahan et al. (1986) SSA source function
        vol_emiss_dist = (/ 0., 10673., 29885., 72482., &
                            786274., 2695160., 977044., 1406417. /)
      ELSEIF(seas_opt == SEASIOANNIDIS .OR. seas_opt == SEASGONG) THEN
        ! Gong et al. (1997) SSA source function
        ! vol_emiss_dist = (/ 2.4470606E+02, 8.0060064E+03, 6.1196694E+04, 1.6120765E+05, &
        !                     7.7027059E+05, 1.2101714E+06, 1.3779869E+05, 5.1837796E+04 /)
        vol_emiss_dist = (/ 5037., 14721., 29885., 72482., &
                            786274., 2695160., 977044., 1406417. /)
      ELSE
        CALL wrf_error_fatal ('seaspray_emis_init: seas_opt not supported')
      ENDIF

    !-- Initialize SSA source functions for MOSAIC-4bin aerosols
    ELSEIF(aer_is_mosaic4bin) THEN
      ! MOSAIC-4bin bin center diameters (micrometer)
      bin_diameters_cen = (/ 0.0781, 0.3125, 1.2500, 5.0000 /)
      ! Size-resolved volume emission scaling factors for MOSAIC-4bin
      IF(seas_opt == SEASMONAHAN) THEN
        ! Monahan et al. (1986) SSA source function
        vol_emiss_dist = (/ 10673., 102367., 3481434., 2383461. /)
      ELSEIF(seas_opt == SEASIOANNIDIS .OR. seas_opt == SEASGONG) THEN
        ! Gong et al. (1997) + O'Dowd et al. (1997) SSA source function
        ! vol_emiss_dist = (/ 8251., 222404., 1980442., 189636. /)
        vol_emiss_dist = (/ 19757., 102367., 3481434., 2383461./)
      ELSE
        CALL wrf_error_fatal ('seaspray_emis_init: seas_opt not supported')
      ENDIF

    ELSE ! aer_is_
      CALL wrf_error_fatal ('seaspray_emis_init: chem_opt not supported')

    ENDIF ! aer_is_

  END SUBROUTINE seaspray_emis_init

END MODULE module_seaspray_emis


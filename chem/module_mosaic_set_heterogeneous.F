    MODULE module_mosaic_set_heterogeneous

! 2019/01/16, Louis Marelle, Jennie L. Thomas
!
! Calculates heterogeneous reaction rates for the halogen mechanism coupled with
! MOSAIC
!
!TODO check the resulting rates
!TODO include a control on Cl- content and on pH: species should react only on
!  acidic aerosols containing some Cl-

    CONTAINS

!-----------------------------------------------------------------------

    SUBROUTINE mosaic_set_heterogeneous(                               &
        chem, t_phy, p_phy, alt,                                       &
        ph_het_hobr, ph_het_brono2, ph_het_hcl, ph_het_hbr,            &
        ph_het_n2o5,                                        &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte                                      )

    ! USE associations
    ! Get MOSAIC data (aerosol density by chem composition, molecular weights,
    ! bin diameters, bin volumes, etc.). dens_* is in (g/cm3) d*_sect etc. are
    ! in (cm), vol_* in cm3
    ! Also get MOSAIC pointers (e.g. numptr_aer) to aerosol number, aerosol
    ! water and aerosol components in chem
    ! Also get some mechanism constants (pi, r8)
    USE module_data_mosaic_asect, ONLY:  &
        dens_water_aer, mw_water_aer, dens_aer, & 
        dlo_sect, dcen_sect, dhi_sect, volumlo_sect, &
        volumhi_sect, ai_phase, ncomp_aer, nsize_aer, &
        numptr_aer, waterptr_aer, massptr_aer
    USE module_configure, ONLY:  num_chem
    USE module_data_mosaic_other, ONLY: pi 
    USE module_data_mosaic_kind,  ONLY: r8

    IMPLICIT NONE

    ! Subroutine arguments
    INTEGER, INTENT(IN) ::   &
        ids, ide, jds, jde, kds, kde,   &
        ims, ime, jms, jme, kms, kme,   &
        its, ite, jts, jte, kts, kte
    REAL, INTENT(IN),   &
        DIMENSION( ims:ime, kms:kme, jms:jme ) :: &
        t_phy, p_phy, alt
    REAL, INTENT(IN),   &
        DIMENSION( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
        chem
    REAL, INTENT(OUT), &
        DIMENSION( ims:ime, kms:kme, jms:jme ) :: ph_het_hobr, ph_het_brono2, &
        ph_het_hcl, ph_het_hbr, ph_het_n2o5

    ! Local variables
    INTEGER it, jt, kt
    INTEGER iphase, itype, ibin, ll, l1, ispecies
    REAL(r8) :: temp, press, airdens
    REAL(r8) :: dumfact, dummass, rnum
    REAL(r8) :: wetdp, wetdens, wetmass, wetvol, &
         drydp, drydens, drymass, dryvol, dryarea, wetarea
    REAL(r8) :: freepath, molar_weight, diffusion_constant
    REAL(r8) :: mean_speed, knudsen, j_corr, upper, lower
    INTEGER,  PARAMETER :: nspecies = 5
    REAL(r8), DIMENSION(nspecies) :: k_het

    ! Parameters and data declarations
    REAL(r8), PARAMETER :: densdefault = 2.0 ! g/cm3
    REAL(r8), PARAMETER :: smallmassaa = 1.0e-19 ! g/kg air
    REAL(r8), PARAMETER :: smallmassbb = 1.0e-29 ! g/kg air
    REAL(r8), PARAMETER :: piover6 = pi/6.0
    REAL(r8), PARAMETER :: onethird = 1.0/3.0
    ! i_<name> are the indices (in the data structures specnames, gamma_factors,
    ! molar_weights) of each species that are interacting with aerosols
    ! If new species are added, new indices should be added here, nspecies
    ! should be raised  and all the data structures below should be updated
    INTEGER,  PARAMETER :: i_hobr = 1
    INTEGER,  PARAMETER :: i_brono2 = 2
    INTEGER,  PARAMETER :: i_hcl = 3
    INTEGER,  PARAMETER :: i_hbr = 4
    INTEGER,  PARAMETER :: i_n2o5 = 5
    ! CHARACTER*6, DIMENSION(nspecies) :: specnames
    REAL(r8), DIMENSION(nspecies) :: gamma_factors
    REAL(r8), DIMENSION(nspecies) :: molar_weights

    ! Data. 
    ! Specnames is not used for now, but needs to be updated for reference when
    ! adding new species
    ! data specnames &
    !   / &
    !   'HOBr  ', 'BrONO2', 'HCl   ', 'HBr   ', 'N2O5  ' &
    !   /
    data gamma_factors &
      / &
      0.1_r8, 0.14_r8, 0.02_r8, 0.02_r8, 0.044_r8 &
      /
    data molar_weights &
      / &
      97._r8, 142._r8, 36.5_r8, 81._r8, 108._r8 &
      /

! ********************** Description of variables ****************************
! *** Subroutine arguments ***
!    t_phy - 3D air temperature (K)
!    p_phy - 3D air pressure (Pa)
!    alt - 3D inverse air density (m3/kg)
!    chem - Advected chemical/aer species (ppm or µg/kg air or #/kg air)
!    ph_het_hobr, ph_het_brono2, ph_het_hcl, ph_het_hbr, ...
!        ph_het_n2o5 - 3D heterogeneous reaction rate arrays for
!        halogen mechanims (min-1)
!
! *** Local variables ***
!    it, jt, kt, iphase, itype, ibin, ll - Loop indices for wrf i dimension, wrf
!        j dimension, wrf k dimension, aerosol phase, aerosol type, aerosol
!        size, aerosol chemical composition
!    temp - Air temperature (K)
!    press - Air pressure (Pa)
!    airdens - Air density (kg/cm3)
!    rnum - Aerosol number (# kg-1 air)
!    wetdp, wetdens, wetmass, wetvol, wetarea - Aerosol wet diameter (cm),
!        density (g cm-3), mass (g kg-1 air), volume (cm-3 kg-1 air), area
!        (cm2/cm3 air)
!    drydp, drydens, drymass, dryvol, dryarea - Same for dry aerosol
!    freepath - Molecular free path in air (cm)
!    molar_weights, molar_weight - Molar weight (g/mol)
!    diffusion_constant - Diffusion constant (m2/s)
!    mean_speed - Mean molecular speed (cm/s)
!    knudsen - Knudsen number ( )
!    j_corr - Fuchs and Sutugin, 1971, correction factor for k_het ( )
!    k_het - Heterogeneous reaction rate (s-1)
!    gamma_factors - probability that a species uptaken by an aerosol undergoes
!        reaction ( )
!    specnames - List of names of halogen mechanism species interacting with
!        aerosols. Not used in the code for now but included for reference.
! ****************************************************************************

!   loop over i and j
    DO jt = jts, jte
    DO it = its, ite
    DO kt = kts, kte

!   Air temperature (K), pressure (Pa), density (kg/cm3)
    temp = t_phy(it,kt,jt)
    press = p_phy(it,kt,jt)
    airdens = 1._r8 / alt(it,kt,jt) * 1.E-6_r8

!   Get aerosol wet density (g/cm3) and wet diameter (cm)
    ! Only works for the interstitial phase
    iphase = ai_phase
    ! For now, only works for itype = 1, ntypes_aer = 1 (one aerosol type)
    itype = 1

    ! Initalize reaction constant k_het to 0 s-1
    k_het(:) = 0._r8

    DO ibin = 1, nsize_aer(itype)
      dryvol = 0.0_r8
      drymass = 0.0_r8
  
      ! Sum over aerosol chemical compounds to calculate total aerosol dry
      ! volume in cm3/kg air. E.g., for bin 1 drymass = bc_a01+so4_a01+... etc
      ! and dryvol = vol_bc_a01+...
      DO ll = 1, ncomp_aer(itype) 
        l1 = massptr_aer(ll,ibin,itype,iphase)
        dummass = chem(it,kt,jt,l1) * 1.0E-6_r8 ! chem is in ug/kg, convert to g/kg air
        drymass = drymass + dummass
        dryvol = dryvol + dummass/dens_aer(ll,itype) ! cm3/kg air
      END DO
  
      ! Add water mass to calculate wet volume in cm3/kg air
      l1 = waterptr_aer(ibin,itype)
      dummass = chem(it,kt,jt,l1) * 1.0E-6_r8 ! chem is in ug/kg, convert to g/kg air
      wetmass = drymass + dummass
      wetvol = dryvol + dummass/dens_water_aer ! cm3/kg air, dens is in g/cm3
  
      ! Retrieve aerosol number concentration
      l1 = numptr_aer(ibin,itype,iphase)
      rnum = chem(it,kt,jt,l1) ! #/kg
  
      ! Calculate the dry and wet diameters drydp and wetdp in cm
      IF (drymass .LE. smallmassbb) THEN
        drydp = dcen_sect(ibin,itype)
        drydens = densdefault
        wetdp = drydp
        wetdens = drydens
  
      ELSE ! (drymass > smallmassbb)
        IF (drymass .LE. smallmassaa) THEN 
          wetmass = drymass
          wetvol = dryvol
        ENDIF
        drydens = drymass/dryvol ! g/cm3
        wetdens = wetmass/wetvol ! g/cm3

        IF (rnum .ge. dryvol/volumlo_sect(ibin,itype)) THEN 
        ! This is capping the dp diameters to [dlo - dhi],
        ! This routine is assuming internal mixing in the bin
          drydp = dlo_sect(ibin,itype)
        ELSE IF (rnum .le. dryvol/volumhi_sect(ibin,itype)) THEN 
          drydp = dhi_sect(ibin,itype) 
        ELSE
          drydp = (dryvol/(piover6*rnum))**onethird
        ENDIF

        IF(abs(wetvol).gt.(1000._r8*abs(dryvol))) THEN
          dumfact=10.0_r8
        ELSE
          dumfact=abs(wetvol/dryvol)**onethird
          dumfact=max(1.0_r8,min(dumfact,10.0_r8))
        ENDIF
        wetdp = drydp*dumfact

      ENDIF !if (drymass .le. smallmassbb)

      ! Calculate mean dry and wet particle surface areas in
      ! cm2(aerosol)/cm3(air)
      dryarea = pi * rnum * airdens * drydp ** 2
      wetarea = pi * rnum * airdens * wetdp ** 2
  
      ! Calculate the Fuchs and Sutugin (1971) expression from Table 11.1 
      ! Seinfeld and Pandis, 1998 
      !TODO define everything
      DO ispecies = 1,nspecies

        molar_weight = molar_weights(ispecies)

        ! Calculate the mean molecular speed (cm/s)
        CALL calc_mean_speed(temp, molar_weight, mean_speed)

        ! Calculate the diffusion constant (m2/s) for species ispecies
        CALL calc_diffusion_constant(temp, press, molar_weight, diffusion_constant)

        ! Calculate the  molecular free path in air (cm)
        CALL calc_freepath(temp, molar_weight, diffusion_constant, freepath)

        ! Calculate the Knudsen number
        knudsen = freepath / wetdp

        ! Calculate the Fuchs and Sutugin (1871) expression j_corr from Table
        ! 11.1, Seinfeld and Pandis, 1998
        ! Only calculate j_corr for species with a reaction probability
        ! gamma_factors(ispecies) greater than zero
        IF (gamma_factors(ispecies) > 0._r8) THEN
            upper = 0.75_r8 * gamma_factors(ispecies) * (1._r8 + knudsen)
            lower = knudsen ** 2 + knudsen + 0.283_r8 * knudsen * gamma_factors(ispecies) &
                    + 0.75_r8 * gamma_factors(ispecies)
            j_corr = upper / lower
        ELSE
            j_corr = 0._r8
        ENDIF

        ! Calculate the effective rate constant k_het (s-1) for surface reactions
        ! Sum k_het over bins to calculate the total resulting effective rate constant
        !TODO check wetarea vs. dryarea with Jennie, check k_het sum over bins
        k_het(ispecies) = k_het(ispecies) + 0.25_r8 * mean_speed * gamma_factors(ispecies) &
                           * wetarea * j_corr
      END DO !ispecies
  
    END DO !ibin, size

    ! Update the 3D reaction constants and convert from s-1 to min-1
    ph_het_hobr(it, kt, jt) = k_het(i_hobr) * 60._r8
    ph_het_brono2(it, kt, jt) = k_het(i_brono2) * 60._r8
    ph_het_hcl(it, kt, jt) = k_het(i_hcl) * 60._r8
    ph_het_hbr(it, kt, jt) = k_het(i_hbr) * 60._r8
    ph_het_n2o5(it, kt, jt) = k_het(i_n2o5) * 60._r8

    END DO !kt
    END DO !it
    END DO !jt

    END SUBROUTINE mosaic_set_heterogeneous

!------------------------------------------------------------------------

    SUBROUTINE calc_diffusion_constant(temp, press, molar_weight, diffusion_constant)

!   Calculates the diffusion constant of the species of molar weight
!   molar_weight
!
!   input parameters
!    temp - Temperature (K)
!    press - Pressure (Pa)
!    molar_weight - Molar weight (g/mol)
!   output parameters
!    diffusion_constant - Diffusion constant (m2/s)
!

    USE module_data_mosaic_kind,  ONLY: r8

    IMPLICIT NONE

    ! Subroutine arguments
    REAL(r8), INTENT(IN) :: temp, press, molar_weight
    REAL(r8), INTENT(OUT) :: diffusion_constant
    ! Local variables
    REAL(r8) :: mw_kg
    ! Parameters
    REAL(r8), PARAMETER :: rgas = 8.3144598_r8 !TODO use the one from model constants (J K-1 mol-1, or kg*m2/s2/K/mol)

    ! Molar weight in kg/mol
    mw_kg =  molar_weight / 1000._r8
    
    ! diffusion_constant calculated as D = mean_free_path in air*molecular_velocity/3
    ! mean_free_path in air calculated according to Seinfeld and Pandis, 1998 edition
    ! equation 8.6
    ! molecular_velocity calculated according to Modeling chemistry in and
    ! above snow at Summit, J.L. Thomas
    diffusion_constant = 2. * 1.8_r8 * 10._r8 ** (-5) * rgas * temp / &
        (3. * press * SQRT(0.029_r8) * SQRT(mw_kg))

    END SUBROUTINE calc_diffusion_constant

!------------------------------------------------------------------------

    SUBROUTINE calc_freepath(temp, molar_weight, diffusion_constant, freepath)

!   Calculates the molecular free path in air, in cm
!
!   Input parameters
!    temp - Temperature (K)
!    diffusion_constant - Diffusion constant (m2/s)
!    molar_weight - Molar weight (g/mol)
!   Output parameters
!    freepath - Molecular freepath in air (cm)
!    

    USE module_data_mosaic_other, ONLY: pi 
    USE module_data_mosaic_kind,  ONLY: r8

    IMPLICIT NONE

    ! Subroutine arguments
    REAL(r8), INTENT(IN) :: temp, molar_weight, diffusion_constant
    REAL(r8), INTENT(OUT) :: freepath
    ! Local variables
    REAL mw_kg, c_bar
    ! Parameters
    REAL(r8), PARAMETER :: rgas = 8.3144598_r8 !TODO use the one from model constants (J K-1 mol-1, or kg*m2/s2/K/mol)
    REAL(r8), PARAMETER :: m_to_cm = 100._r8
    
    ! Molar weight in kg/mol
    mw_kg =  molar_weight / 1000._r8

     
    ! Mean speed of molecules in m/s
    c_bar = SQRT(8._r8 * rgas * temp / (pi * mw_kg))

    ! Calculate mean free path (cm) based on Fuchs and Sutugin (1971), equation
    ! from Seinfeld and Pandis, 1998, Table 11.1
    freepath = 3._r8 * diffusion_constant / c_bar * m_to_cm

    END SUBROUTINE calc_freepath

!------------------------------------------------------------------------

    SUBROUTINE calc_mean_speed(temp, molar_weight, mean_speed)

!   Calculates the mean speed of molecules in cm/s
!
!   Input parameters
!    temp - Temperature (K)
!    molar_weight - Molar weight (g/mol)
!   Output parameters
!    mean_speed - Mean molecular speed (cm/s)

    USE module_data_mosaic_other, ONLY: pi 
    USE module_data_mosaic_kind,  ONLY: r8

    IMPLICIT NONE

    ! Subroutine arguments
    REAL(r8), INTENT(IN) :: temp, molar_weight
    REAL(r8), INTENT(OUT) :: mean_speed
    ! Local variables
    REAL(r8) :: mw_kg
    ! Parameters
    REAL(r8), PARAMETER :: rgas = 8.3144598_r8 !TODO use the one from model constants (J K-1 mol-1, or kg*m2/s2/K/mol)
    REAL(r8), PARAMETER :: m_to_cm = 100._r8 

    ! Molar weight in kg/mol
    mw_kg =  molar_weight / 1000._r8
         
    ! Mean speed of molecules (cm/s), equation PP from Pits & Pits (book) in Box
    ! 5.2
    mean_speed = SQRT(8._r8 * rgas * temp / (pi * mw_kg)) * m_to_cm
    
    END SUBROUTINE calc_mean_speed

!------------------------------------------------------------------------

    END MODULE module_mosaic_set_heterogeneous


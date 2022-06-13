    MODULE module_mosaic_set_heterogeneous

! 2019/01/16, Louis Marelle, Jennie L. Thomas
!
! Calculates heterogeneous reaction rates for the halogen mechanism coupled with
! MOSAIC
!
! References:
! - Louis Marelle, Jennie L. Thomas, Shaddy Ahmed, Katie Tuite, Jochen Stutz,
! Aurélien Dommergue, William R. Simpson, Markus M.  Frey, Foteini Baladima.
! Implementation and impacts of surface and blowing snow sources of Arctic
! bromine activation within WRF-Chem 4.1.1. Submitted to the Journal of Advances
! in Modeling Earth Systems.
! - H. M. Amos, D. J. Jacob, C. D. Holmes, et al., Gas-particle partitioning of 
! atmospheric Hg(II) and its effect on global mercury deposition, Atmospheric
! Chemistry and Physics, 2012.
!
! Revision history:
! Ahmed S. 2021/11/26: Include mercury heterogeneous uptake of HgII(g) species to 
!                      aerosols and liquid clouds.
! Ahmed S. 2022/05/30: Add HgII(g) gas-particle partitioning based on the empirical 
!                      equilibrium of Amos et al. (2012).
    CONTAINS

!-----------------------------------------------------------------------

    SUBROUTINE mosaic_set_heterogeneous(                               &
        chem, t_phy, p_phy, alt, dtstepc,                              &
        ph_het_hobr, ph_het_brono2, ph_het_hcl, ph_het_hbr,            &
        ph_het_n2o5, ph_het_hocl, ph_het_clono2, ph_het_oh,            &
        ph_het_hgbroh, ph_het_hgbrno2, ph_het_hgbrho2, ph_het_hgbrbro, &
        ph_het_hgbrclo, ph_het_hgbr2, ph_het_hgbrcl, ph_het_hgcloh,    &
        ph_het_hgclno2, ph_het_hgclho2, ph_het_hgclbro, ph_het_hgclclo,&
        ph_het_hgohoh, ph_het_hgohno2, ph_het_hgohho2, ph_het_hgohbro, &
        ph_het_hgohclo, ph_het_hgx,                                    &
        aer_ph_01, aer_ph_02, aer_ph_03, aer_ph_04, &
        aer_ph_05, aer_ph_06, aer_ph_07, aer_ph_08, &
        qcloud, qndrop, cldfra, pm2_5_water,        &
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
    USE module_configure, ONLY:  num_chem, p_hgbroh, p_hgbrno2, p_hgbrho2, &
            p_hgbrbro, p_hgbrclo, p_hgbr2, p_hgbrcl, p_hgcloh, p_hgclno2, &
            p_hgclho2, p_hgclbro, p_hgclclo, p_hgohoh, p_hgohno2, p_hgohho2, &
            p_hgohbro, p_hgohclo, p_hgx, p_hgp
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
    REAL, INTENT(INOUT),   &
        DIMENSION( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
        chem
    REAL, INTENT(IN) :: dtstepc
    REAL, INTENT(OUT), &
        DIMENSION( ims:ime, kms:kme, jms:jme ) :: ph_het_hobr, ph_het_brono2, &
        ph_het_hcl, ph_het_hbr, ph_het_n2o5, ph_het_hocl, ph_het_clono2, &
        ph_het_oh, ph_het_hgbroh, ph_het_hgbrno2, ph_het_hgbrho2, &
        ph_het_hgbrbro, ph_het_hgbrclo, ph_het_hgbr2, ph_het_hgbrcl, &
        ph_het_hgcloh, ph_het_hgclno2, ph_het_hgclho2, ph_het_hgclbro, &
        ph_het_hgclclo, ph_het_hgohoh, ph_het_hgohno2, ph_het_hgohho2, &
        ph_het_hgohbro, ph_het_hgohclo, ph_het_hgx
    REAL, INTENT(OUT), &
        DIMENSION( ims:ime, kms:kme, jms:jme ) :: aer_ph_01, aer_ph_02, aer_ph_03, aer_ph_04, &
        aer_ph_05, aer_ph_06, aer_ph_07, aer_ph_08
    REAL, INTENT(IN), &
        DIMENSION( ims:ime, kms:kme, jms:jme ) :: pm2_5_water
    REAL, INTENT(IN),   &
        DIMENSION( ims:ime, kms:kme, jms:jme) :: &
        qcloud, qndrop, cldfra

    ! Local variables
    INTEGER it, jt, kt
    INTEGER iphase, itype, ibin, ll, l1, ispecies
    REAL(r8) :: temp, press, airdens
    REAL(r8) :: dumfact, dummass, rnum
    REAL(r8) :: wetdp, wetdens, wetmass, wetvol, &
         drydp, drydens, drymass, dryvol, dryarea, wetarea
    REAL(r8) :: freepath, molar_weight, diffusion_constant
    REAL(r8) :: mean_speed, knudsen, j_corr, upper, lower, aerosol_bin_ph
    REAL(r8) :: cld_num, drop_num, cld_fra, cld_area
    REAL(r8) :: droplet_mass, droplet_volume, reff_droplet
    REAL(r8) :: HgGasTot, HgAer, HgGasAerTot, HgGasTot_eq, dHgGasTot_eq, &
                HgX_conc, new_HgX_conc, new_HgAer_conc
    REAL(r8) :: pm25, K_part, gas_frac, conv_ppm_molcm3
    INTEGER,  PARAMETER :: nspecies = 26
    REAL(r8), DIMENSION(nspecies) :: k_het
    REAL(r8), DIMENSION(nspecies) :: k_het_cld
    REAL(r8) :: k_het_hgx

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
    INTEGER,  PARAMETER :: i_hocl = 6
    INTEGER,  PARAMETER :: i_clono2 = 7
    INTEGER,  PARAMETER :: i_oh = 8
    INTEGER,  PARAMETER :: i_hgbroh = 9
    INTEGER,  PARAMETER :: i_hgbrno2 = 10
    INTEGER,  PARAMETER :: i_hgbrho2 = 11
    INTEGER,  PARAMETER :: i_hgbrbro = 12
    INTEGER,  PARAMETER :: i_hgbrclo = 13
    INTEGER,  PARAMETER :: i_hgbr2 = 14
    INTEGER,  PARAMETER :: i_hgbrcl = 15
    INTEGER,  PARAMETER :: i_hgcloh = 16
    INTEGER,  PARAMETER :: i_hgclno2 = 17
    INTEGER,  PARAMETER :: i_hgclho2 = 18
    INTEGER,  PARAMETER :: i_hgclbro = 19
    INTEGER,  PARAMETER :: i_hgclclo = 20
    INTEGER,  PARAMETER :: i_hgohoh = 21
    INTEGER,  PARAMETER :: i_hgohno2 = 22
    INTEGER,  PARAMETER :: i_hgohho2 = 23
    INTEGER,  PARAMETER :: i_hgohbro = 24
    INTEGER,  PARAMETER :: i_hgohclo = 25
    INTEGER,  PARAMETER :: i_hgx = 26

    ! CHARACTER*6, DIMENSION(nspecies) :: specnames
    REAL(r8), DIMENSION(nspecies) :: gamma_factors
    REAL(r8), DIMENSION(nspecies) :: molar_weights

    ! Data. 
    ! Specnames is not used for now, but needs to be updated for reference when
    ! adding new species
    ! data specnames &
    !   / &
    !   'HOBr   ', 'BrONO2 ', 'HCl    ', 'HBr    ', 'N2O5   ', &
    !   'HOCl   ', 'ClONO2 ', 'OH     ', 'HgBrOH ', 'HgBrNO2', &
    !   'HgBrHO2', 'HgBrBrO', 'HgBrClO', 'HgBr2  ', 'HgBrCl ', &
    !   'HgClOH ', 'HgClNO2', 'HgClHO2', 'HgClBrO', 'HgClClO', &
    !   'HgOHOH ', 'HgOHNO2', 'HgOHHO2', 'HgOHBrO', 'HgOHClO', &
    !   'HgX (treated as HgCl2)' /
    data gamma_factors &
      / &
      0.1_r8, 0.14_r8, 0.02_r8, 0.02_r8, 0.044_r8, &
      4.0E-4_r8, 0.11_r8, 0.2_r8, 0.1_r8, 0.1_r8,  &
      0.1_r8, 0.1_r8, 0.1_r8, 0.1_r8, 0.1_r8, &
      0.1_r8, 0.1_r8, 0.1_r8, 0.1_r8, 0.1_r8, &
      0.1_r8, 0.1_r8, 0.1_r8, 0.1_r8, 0.1_r8, &
      0.1_r8 &
      /
    data molar_weights &
      / &
      97._r8, 142._r8, 36.5_r8, 81._r8, 108._r8, &
      52.5_r8, 97.5_r8, 17._r8, 297.5_r8, 326.5_r8, &
      313.5_r8, 376.5_r8, 332._r8, 360.5_r8, 316._r8, &
      253._r8, 282._r8, 269._r8, 332._r8, 287.5_r8, &
      234.5_r8, 263.5_r8, 250.5_r8, 313.5_r8, 269._r8, &
      271.5_r8 &
      /

! ********************** Description of variables ****************************
! *** Subroutine arguments ***
!    t_phy - 3D air temperature (K)
!    p_phy - 3D air pressure (Pa)
!    alt - 3D inverse air density (m3/kg)
!    chem - Advected chemical/aer species (ppm or µg/kg air or #/kg air)
!    dtstepc - Chemical timestep (s)
!    ph_het_hobr, ph_het_brono2, ph_het_hcl, ph_het_hbr, ...
!        ph_het_n2o5 - 3D heterogeneous reaction rate arrays for
!        halogen mechanims (min-1)
!    pm2_5_water - PM2.5 wet concentration (ug m-3)
!    qcloud - Cloud water mixing ratio (kg kg-1)
!    qndrop - Droplet number mixing ratio (# kg-1)
!    cldfra - Cloud fraction (0 - 1)
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
!    cld_num, drop_num, cld_fra - Cloud water mixing ratio (kg kg-1), droplet
!    number mixing ratio (# kg-1), cloud fraction (0 - 1)
!    cld_area - Effective cloud surface area (cm2/cm3 air)
!    droplet_mass, droplet_volume, reff_droplet - Mass of a single cloud droplet
!    (g), volume of a single cloud droplet (cm3), effective mean radius of a
!    single cloud droplet (cm)
!    HgGasTot -  Total HgII gas-phase concentration (molec cm-3)
!    HgAer - Aerosol-phase Hg concentration (HgP, molec cm-3 )
!    HgGasAerTot - Sum of total HgII gas and particle phases (molec cm-3)
!    HgGasTot_eq - Amount of HgII gas in equilibrium with the aerosol phase, following
!                  gas-particle partitioning function of Amos et al. (2012) (molec cm-3)
!    dHgGasTot_eq - Amount of HgII gas formed after volatilization from the aerosol-phase
!                   to the unspeciated HgII gas-phase species HgX (molec cm-3)
!    HgX_conc - Concentration of the unspeciated HgII gas-phase species formed
!               via volatilization of HgII on aerosolsosols (molec cm-3)
!    new_HgX_conc, new_HgAer_conc - Updated concentrations of HgX and HgP following volatilization (molec cm-3)
!    pm25 - PM2.5 wet concentration (ug m-3)
!    K_part - HgII gas-particle partitioning equilibrium function following Amos et al. (2012)
!    gas_frac - Amount of HgII in the gas-phase (unitless)
!    conv_ppm_molcm3 - Conversion factor from ppm to molec cm-3
!    k_het - Heterogeneous reaction rate (s-1)
!    k_het_cld - Heterogeneous reaction rate to clouds (s-1)
!    k_het_hgx - Mass transfer rate coefficient (s-1) for volatilization of HgII(g) from HgP
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

!   Get cloud water mixing ratio (kg kg-1), droplet number mixing ratio (# kg-1) 
!   and cloud fraction (unitless) 
    cld_num = qcloud(it,kt,jt)
    drop_num = qndrop(it,kt,jt)
    cld_fra = cldfra(it,kt,jt)

!   Get aerosol wet density (g/cm3) and wet diameter (cm)
    ! Only works for the interstitial phase
    iphase = ai_phase
    ! For now, only works for itype = 1, ntypes_aer = 1 (one aerosol type)
    itype = 1

    ! Initalize reaction constant k_het, k_het_cld, and k_het_hgx to 0 s-1
    k_het(:) = 0._r8
    k_het_cld(:) = 0._r8
    k_het_hgx = 0.0_r8

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
      DO ispecies = 1,nspecies

        molar_weight = molar_weights(ispecies)

        ! Calculate the mean molecular speed (cm/s)
        CALL calc_mean_speed(temp, molar_weight, mean_speed)

        ! Calculate the diffusion constant (m2/s) for species ispecies
        CALL calc_diffusion_constant(temp, press, molar_weight, diffusion_constant)

        ! Calculate the  molecular free path in air (cm)
        CALL calc_freepath(temp, molar_weight, diffusion_constant, freepath)
        
        ! Calculate heterogeneous uptake for halogens + mercury species
        IF (ispecies .le. 8) THEN   ! Calculate uptake for halogen species

          ! Calculate the Knudsen number
          knudsen = freepath * 2.0_r8 / wetdp

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

          ! Aerosol pH, calculated in module_data_mosaic_therm
          IF(ibin == 1) THEN
            aerosol_bin_ph = aer_ph_01(it, kt, jt)
          ELSEIF(ibin == 2) THEN
            aerosol_bin_ph = aer_ph_02(it, kt, jt)
          ELSEIF(ibin == 3) THEN
            aerosol_bin_ph = aer_ph_03(it, kt, jt)
          ELSEIF(ibin == 4) THEN
            aerosol_bin_ph = aer_ph_04(it, kt, jt)
          ELSEIF(ibin == 5) THEN
            aerosol_bin_ph = aer_ph_05(it, kt, jt)
          ELSEIF(ibin == 6) THEN
            aerosol_bin_ph = aer_ph_06(it, kt, jt)
          ELSEIF(ibin == 7) THEN
            aerosol_bin_ph = aer_ph_07(it, kt, jt)
          ELSEIF(ibin == 8) THEN
            aerosol_bin_ph = aer_ph_08(it, kt, jt)
          ENDIF
          ! Turn off the het reaction if the aerosol is not sufficiently acidic
          ! (using a ph == 5 threshold)
          IF(aerosol_bin_ph > 5._r8 .OR. aerosol_bin_ph <= 0._r8) THEN
            j_corr = 0._r8
          ENDIF

          ! Calculate the effective rate constant k_het (s-1) for surface reactions
          ! Sum k_het over bins to calculate the total resulting effective rate constant
          k_het(ispecies) = k_het(ispecies) + 0.25_r8 * mean_speed * gamma_factors(ispecies) &
                             * wetarea * j_corr
        ELSE    ! Calculate uptake for mercury species
          ! Uptake to aerosols 
          ! Calculate k_het for uptake to aerosols
          k_het(ispecies) = k_het(ispecies) + wetarea * &
          (((0.5_r8*wetdp)/diffusion_constant)+(4.0_r8/(mean_speed*gamma_factors(ispecies))))**-1.0_r8

          ! Uptake to clouds
          ! Do this once to make sure it is not calculated for each aerosol size bin.
          IF (ibin == 1) THEN
            ! Check cloud fraction in each grid cell, if there are no clouds return khet as 0 for that grid cell
            IF (cld_fra == 0._r8) THEN
              k_het_cld = 0._r8
            ELSE
              ! Calculate the liquid cloud droplet effective mean radius (cm) and surface area (cm2/cm3 air)
              CALL calc_cloud_properties(cld_num, drop_num, airdens, reff_droplet, cld_area)

              k_het_cld(ispecies) = k_het_cld(ispecies) + cld_area * &
                        ((reff_droplet/diffusion_constant)+(4.0_r8/(mean_speed*gamma_factors(ispecies))))**-1.0_r8
            END IF
          ELSE ! Set to zero for other aerosol bins so its not counted multiple times
            k_het_cld = 0._r8
          END IF
          ! Calculate total heterogeneous uptake of mercury species
          k_het(ispecies) = k_het(ispecies) + k_het_cld(ispecies)

        END IF ! If Hg or halogen species           

      END DO !ispecies
  
      ! Volatilization of particle-bound mercury (HgP) to unspeciated gas-phase HgII (HgX(g))
      ! Following the empirical gas-particle equilibrium of Amos et al. (2012)
      ! HgX(g) is treated as HgCl2.
      ! Do volatilization for fine-mode aerosols (<2.5 um)
      IF (ibin .LE. 6) THEN
        ! Calculate the mean molecular speed (cm/s)
        CALL calc_mean_speed(temp, molar_weights(i_hgx), mean_speed)

        ! Calculate the diffusion constant (m2/s) for species ispecies
        CALL calc_diffusion_constant(temp, press, molar_weights(i_hgx), diffusion_constant)

        ! Calculate the  molecular free path in air (cm)
        CALL calc_freepath(temp, molar_weights(i_hgx), diffusion_constant, freepath)

        ! Calculate mercury mass transfer from aerosol to gas-phase
        ! Calculate k_het for HgX and sum over aerosol bins 1-6
        k_het_hgx = k_het_hgx + wetarea * &
        (((0.5_r8*wetdp)/diffusion_constant)+(4.0_r8/(mean_speed*gamma_factors(i_hgx))))**-1.0_r8
      END IF ! Volatilization of fine aerosols

    END DO !ibin, size
    
    ! Update concentrations of HgX and HgP following volatilization
    ! Initialize total gas-phase HgII concentration
    HgGasTot = 0.0_r8

    ! Get total gas-phase concentration of HgII species
    HgGasTot = HgGasTot + chem(it,kt,jt,p_hgbroh) + chem(it,kt,jt,p_hgbrno2) &
                        + chem(it,kt,jt,p_hgbrho2) + chem(it,kt,jt,p_hgbrbro) &
                        + chem(it,kt,jt,p_hgbrclo) + chem(it,kt,jt,p_hgbr2) &
                        + chem(it,kt,jt,p_hgbrcl) + chem(it,kt,jt,p_hgcloh) &
                        + chem(it,kt,jt,p_hgclno2) + chem(it,kt,jt,p_hgclho2) &
                        + chem(it,kt,jt,p_hgclbro) + chem(it,kt,jt,p_hgclclo) &
                        + chem(it,kt,jt,p_hgohoh)  + chem(it,kt,jt,p_hgohno2) &
                        + chem(it,kt,jt,p_hgohho2) + chem(it,kt,jt,p_hgohbro) &
                        + chem(it,kt,jt,p_hgohclo) + chem(it,kt,jt,p_hgx)

    ! Get the conversion from ppmv to molec cm-3
    CALL ppmv_to_molcm3(temp, press, conv_ppm_molcm3)

    ! Convert chem units from ppmv to molec cm-3
    HgGasTot = HgGasTot * conv_ppm_molcm3

    ! Get PM2.5 concentration (in ug m-3)
    pm25 = pm2_5_water(it,kt,jt)

    ! Calculate gas-particle partitioning coefficient (m-3/ug) based on Amos et al. (2012)
    K_part = 10.0_r8 ** ((2500.0_r8/temp) - 10.0_r8)

    ! Calculate the HgII gas fraction
    gas_frac = 1.0_r8 / (1.0_r8 + (K_part * pm25))

    ! Get aerosol-phase concentration (in molec cm-3)
    HgAer = chem(it,kt,jt,p_hgp) * conv_ppm_molcm3

    ! Sum the HgII gas-phase and aerosol concentrations 
    HgGasAerTot = HgGasTot + HgAer

    ! Calculate the concentration of HgII gas at equilibrium
    HgGasTot_eq = HgGasAerTot * gas_frac

    ! Calculate the change in concentration of HgII at equilibrium after volatilization
    ! This is the amount of HgII gas lost from particles and forming HgX
    dHgGasTot_eq = HgGasTot_eq * (1.0_r8 - EXP(- k_het_hgx * dtstepc))

    ! Make sure amount of outgas is less than amount of HgP present
    dHgGasTot_eq = MIN( dHgGasTot_eq, HgAer )

    ! Add this change of volatilized mercury to the HgX species
    HgX_conc = chem(it,kt,jt,p_hgx)*conv_ppm_molcm3
    new_HgX_conc = HgX_conc + dHgGasTot_eq
    ! Convert back from molecules/cm3 to ppmv
    chem(it,kt,jt,p_hgx) = new_HgX_conc / conv_ppm_molcm3

    ! Update Hg particle concentration
    new_HgAer_conc = HgAer - dHgGasTot_eq
    ! Convert back from molecules/cm3 to ppmv
    chem(it,kt,jt,p_hgp) = new_HgAer_conc / conv_ppm_molcm3


    ! Update the 3D reaction constants and convert from s-1 to min-1
    !LMarelle TODO: this is passed to KPP using the existing code for passing
    !photolysis rates, which explains the ph_ naming. In the future, k_het
    !should be passed directly to KPP using the include files in the coupler
    ph_het_hobr(it, kt, jt) = k_het(i_hobr) * 60._r8
    ph_het_brono2(it, kt, jt) = k_het(i_brono2) * 60._r8
    ph_het_hcl(it, kt, jt) = k_het(i_hcl) * 60._r8
    ph_het_hbr(it, kt, jt) = k_het(i_hbr) * 60._r8
    ph_het_n2o5(it, kt, jt) = k_het(i_n2o5) * 60._r8
    ph_het_hocl(it, kt, jt) = k_het(i_hocl) * 60._r8
    ph_het_clono2(it, kt, jt) = k_het(i_clono2) * 60._r8
    ph_het_oh(it, kt, jt) = k_het(i_oh) * 60._r8
    ph_het_hgbroh(it, kt, jt) = k_het(i_hgbroh) * 60._r8
    ph_het_hgbrno2(it, kt, jt) = k_het(i_hgbrno2) * 60._r8
    ph_het_hgbrho2(it, kt, jt) = k_het(i_hgbrho2) * 60._r8
    ph_het_hgbrbro(it, kt, jt) = k_het(i_hgbrbro) * 60._r8
    ph_het_hgbrclo(it, kt, jt) = k_het(i_hgbrclo) * 60._r8
    ph_het_hgbr2(it, kt, jt) = k_het(i_hgbr2) * 60._r8
    ph_het_hgbrcl(it, kt, jt) = k_het(i_hgbrcl) * 60._r8
    ph_het_hgcloh(it, kt, jt) = k_het(i_hgcloh) * 60._r8
    ph_het_hgclno2(it, kt, jt) = k_het(i_hgclno2) * 60._r8
    ph_het_hgclho2(it, kt, jt) = k_het(i_hgclho2) * 60._r8
    ph_het_hgclbro(it, kt, jt) = k_het(i_hgclbro) * 60._r8
    ph_het_hgclclo(it, kt, jt) = k_het(i_hgclclo) * 60._r8
    ph_het_hgohoh(it, kt, jt) = k_het(i_hgohoh) * 60._r8
    ph_het_hgohno2(it, kt, jt) = k_het(i_hgohno2) * 60._r8
    ph_het_hgohho2(it, kt, jt) = k_het(i_hgohho2) * 60._r8
    ph_het_hgohbro(it, kt, jt) = k_het(i_hgohbro) * 60._r8
    ph_het_hgohclo(it, kt, jt) = k_het(i_hgohclo) * 60._r8
    ph_het_hgx(it, kt, jt) = k_het(i_hgx) * 60._r8

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
    !TODO use the value from model_constants
    REAL(r8), PARAMETER :: rgas = 8.3144598_r8 !(J K-1 mol-1, or kg*m2/s2/K/mol)

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
    !TODO use the value from model_constants
    REAL(r8), PARAMETER :: rgas = 8.3144598_r8 !(J K-1 mol-1, or kg*m2/s2/K/mol)
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
    !TODO use the value from model_constants
    REAL(r8), PARAMETER :: rgas = 8.3144598_r8 !(J K-1 mol-1, or kg*m2/s2/K/mol)
    REAL(r8), PARAMETER :: m_to_cm = 100._r8 

    ! Molar weight in kg/mol
    mw_kg =  molar_weight / 1000._r8
         
    ! Mean speed of molecules (cm/s), equation PP from Pits & Pits (book) in Box
    ! 5.2
    mean_speed = SQRT(8._r8 * rgas * temp / (pi * mw_kg)) * m_to_cm
    
    END SUBROUTINE calc_mean_speed

!------------------------------------------------------------------------

    SUBROUTINE calc_cloud_properties(cld_num, drop_num, airdens, reff_droplet, cld_area)

!   Calculates the mean radius of liquid cloud droplet (in cm) and mean cloud surface area (in cm2/cm3 air)
!
!   Input parameters
!    cld_num - Cloud water mixing ratio (kg kg-1)
!    drop_num - droplet number mixing ratio (# kg-1)
!    airdens - Air density (kg/cm3)
!   Output parameters
!    reff_droplet - Effective mean radius of a single cloud droplet (cm)
!    cld_area - Effective cloud surface area (cm2/cm3 air)

    USE module_data_mosaic_other, ONLY: pi
    USE module_data_mosaic_kind,  ONLY: r8

    IMPLICIT NONE

    ! Subroutine arguments
    REAL(r8), INTENT(IN) :: cld_num, drop_num, airdens
    REAL(r8), INTENT(OUT) :: reff_droplet, cld_area
    ! Local variables
    REAL(r8) :: droplet_mass
    REAL(r8) :: droplet_volume
    ! Parameters
    REAL(r8), PARAMETER :: kg_to_g = 1000._r8
    REAL(r8), PARAMETER :: water_density = 1.0 ! g/cm3

    ! Calculate the mass of a single cloud droplet (in g)
    droplet_mass = (cld_num/drop_num)*kg_to_g

    ! Convert mass to volume (in cm3) by using a water density of 1g/cm3
    droplet_volume = droplet_mass / water_density

    ! Calculate the effective mean radius (in cm), assuming a spherical
    ! droplet
    reff_droplet = (3./(4.*pi) * droplet_volume) ** (1./3.)

    ! Calculate effective mean cloud surface area (cm2/cm3 air)
    cld_area = 4. * pi * drop_num * airdens * reff_droplet ** 2

    END SUBROUTINE calc_cloud_properties

!------------------------------------------------------------------------

    SUBROUTINE ppmv_to_molcm3(temp, press, conv_ppm_molcm3)

!   Calculate conversion factor from ppmv to molecules/cm3.
!
!   Input parameters
!    temp - Temperature (K)
!    press - Pressure (Pa)
!   Output parameters
!    conv_ppm_molcm3 - Conversion factor from ppmv to molecules/cm3 

    USE module_data_mosaic_kind,  ONLY: r8

    IMPLICIT NONE

    ! Subroutine arguments
    REAL(r8), INTENT(IN) :: temp, press
    REAL(r8), INTENT(OUT) :: conv_ppm_molcm3
    ! Local variables
    REAL(r8) :: air_vol
    ! Parameters
    !TODO use the value from model_constants
    REAL(r8), PARAMETER :: rgas = 8.3144598E6_r8 !(cm3 Pa K-1 mol-1)
    REAL(r8), PARAMETER :: avogadro = 6.022E23_r8 ! (molecules mol-1)
    REAL(r8), PARAMETER :: ppm_vmr = 1.0E6_r8 ! Parts per million to VMR

    ! Volume of air (in cm3)
    air_vol =  (rgas * temp) / press

    ! Conversion factor from ppmv to molecules/cm3
    conv_ppm_molcm3 = avogadro/(air_vol * ppm_vmr)

    END SUBROUTINE ppmv_to_molcm3

!------------------------------------------------------------------------

    END MODULE module_mosaic_set_heterogeneous


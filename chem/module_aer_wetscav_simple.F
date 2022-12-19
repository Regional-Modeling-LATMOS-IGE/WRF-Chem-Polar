! 2022/12/09 Louis Marelle for CICERO
!
! Purpose:
! Perform wet scavenging of aerosols by grid scale rainout/washout
!
! Remarks:
! Washout and rainout from Luo et al. (2019, 2020)
! TODO - mosaic aerosol should be deposited assuming internal mixing (using
!        the volume-averaged activation efficiency in each bin)
!
! References:
! ==============================================================================
! (1) Luo, G., Yu, F., and Moch, J. M.: Further improvement of wet process
!     treatments in GEOS-Chem v12.6.0: impact on global distributions of
!     aerosols and aerosol precursors, Geosci. Model Dev., 13, 2879–2903,
!     https://doi.org/10.5194/gmd-13-2879-2020, 2020.
! (2) Luo, G., Yu, F., and Schwab, J.: Revised treatment of wet scavenging
!     processes dramatically improves GEOS-Chem 12.0.0 simulations of surface
!     nitric acid, nitrate, and ammonium over the United States, Geosci. Model
!     Dev., 12, 3439–3447, https://doi.org/10.5194/gmd-12-3439-2019, 2019.
!
!TODO use reasonable thresholds instead of 0.0 in all the module
!TODO threshold on chem, do not use wet scav if aerosol conc is below pg/kg or something similar
!
!TODO remove unused arguments, e.g. dx, dy, config_flags, rho_phy, p_phy
!TODO maybe use double precision
!TODO qc,qi,qr,qs,qg should use b4mp arrays
!

MODULE module_aer_wetscav_simple

   ! USE module_wrf_error, ONLY : wrf_debug

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: aer_wetscav_simple, aer_wetscav_simple_init

   !---- Global module-wide variables
   ! chem_opt aerosol model family
   LOGICAL :: aer_is_mosaic

   !-- Allocatable size-resolved aerosol arrays
   ! Number of aerosol size bins/modes for chem_opt in use
   INTEGER :: naersizes

   !-- Allocatable species-resolved aerosol arrays
   ! Dim is (naerspecies,naersizes), number of aerosol species, number of
   ! aerosol sizes, and depends on aerosol mechanism. These are initialized in
   ! aer_wetscav_simple_init
   !
   ! Number of wet-scavenged aerosol species for chem_opt (dim of
   ! aer_pointers)
   INTEGER :: naerspecies
   ! List of aerosol species experiencing rainout/washout
   INTEGER, DIMENSION(:, :), ALLOCATABLE :: aer_pointers
   ! Rainout activation efficiency of species (0-1)
   REAL, DIMENSION(:, :), ALLOCATABLE :: aq_aer_ratio, ice_aer_ratio
   ! List of fine/coarse aerosol species in mechanism for washout (r > or <= 1 micrometer)
   LOGICAL, DIMENSION(:, :), ALLOCATABLE :: aer_is_fine, aer_is_coarse
   ! Deposition flux of species during timestep (ug/m2)
   REAL, DIMENSION(:, :), ALLOCATABLE :: deposition_flux

CONTAINS

!----------------------------------------------------------------------------

   SUBROUTINE aer_wetscav_simple( id, dtstep, config_flags,                &
                              t_phy, p_phy, dz8w,                          &
                              chem, rho_phy, cldfra, rainprod, evapprod,   &
                              precr, preci, precs, precg, dx, dy,          &
                              qv, qc, qi, qs, qr,                          &
                              ids,ide, jds,jde, kds,kde,                   &
                              ims,ime, jms,jme, kms,kme,                   &
                              its,ite, jts,jte, kts,kte                    )

   ! Purpose:
   ! Calculates and removes wet scavenging tendencies for aerosols
   !
   ! Arguments:
   ! Input:
   !   id = domain ID
   !   dtstep = model time step (s)
   !   config_flags = namelist configuration flags
   !   t_phy = 3D physical temperature (K)
   !   p_phy = 3D pressure (Pa)
   !   dz8w = 3D vertical level depth (m)
   !   rho_phy = 3D air density (kg/m3)
   !   cldfra = 3D cloud fraction (unitless, 0-1)
   !   rainprod = 3D rain production rate (kg/kg dry air/s)
   !   evapprod = 3D rain evaporation rate (kg/kg dry air/s)
   !   precr = 3D rain precipitation rate through grid cell (kg/m2/s)
   !   preci = 3D ice precipitation rate through grid cell (kg/m2/s)
   !   precs = 3D snow precipitation rate through grid cell (kg/m2/s)
   !   precg = 3D graupel precipitation rate through grid cell (kg/m2/s)
   !   dx, dy = model resolution in x and y (km)
   !   qv = specific humidity = water vapor mixing ratio (kg/kg dry air)
   !   qc = cloud water mixing ratio (kg/kg dry air)
   !   qi = cloud ice mixing ratio (kg/kg dry air)
   !   qs = snow water mixing ratio (kg/kg dry air)
   !   qr = rain water mixing ratio (kg/kg dry air)
   !   ids,ims,its etc. = domain, memory, and tile start and end indices
   ! Input/output:
   !   chem = advected chemical species (gases in ppm, aerosols in ug/kg)

   ! USE associations:
   !TODO Use ONLY:
   USE module_configure
   USE module_state_description
   ! p_seas_1, p_seas_2, p_seas_3, p_seas_4, p_bc2, p_bc1, p_oc2, p_dust_1,
   ! p_dust_2, p_dust_3, p_dust_4, p_sulf from module_state_description

   IMPLICIT NONE

   !---- Arguments
   TYPE(grid_config_rec_type), INTENT(IN) :: config_flags
   INTEGER, INTENT(IN) :: ids,ide, jds,jde, kds,kde,    &
                          ims,ime, jms,jme, kms,kme,    &
                          its,ite, jts,jte, kts,kte,    &
                          id
   REAL, INTENT(IN) :: dx, dy
   REAL, INTENT(IN) :: dtstep
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme)         ,        &
       INTENT(IN) :: t_phy, p_phy, rho_phy, dz8w
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
       INTENT(IN) :: QV, QI, QC, QS, QR
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
       INTENT(IN) :: rainprod, evapprod
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
       INTENT(IN) :: precr, preci, precs, precg
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme)         ,        &
       INTENT(IN) :: cldfra
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem),          &
       INTENT(INOUT) :: chem

   !---- Local variables and parameters
   ! x,y,z, indices; aerosol size index; aerosol species index; chemical species pointer
   INTEGER :: i, j, k, isize, ispec, p_chem
   ! Cell fractions where rainout/washout is applied (See Luo et al., 2020)
   REAL :: f, fr
   ! Vertical dimension of cldfra for cloud fraction calculation
   INTEGER :: dimz
   ! Aerosol conc. in cloud phase, aerosol conc. lost to rainout, aerosols
   ! conc. resuspended from precip (ug/kg)
   REAL :: cloudborne_aer, rainout_aer, resusp_aer
   ! Condensed water mixing ratio (kg/kg dry air), precipitation mixing ratio
   ! (kg/kg dry air), liquid fraction of condensed water (0-1)
   REAL :: condensed_water, precip_water, liquid_fraction
   ! Rainout activation efficiency of species (0-1)
   REAL :: aq_aer_ratio_ispec, ice_aer_ratio_ispec
   ! Total 3D precipitation trough cell for rain/snow/ice (kg/m2/hr), rain
   ! production and rain evaporation rate at i,k,j (kg/kg/s)
   REAL :: ppr, pps, ppi, rainrate, evaprate
   ! Temperature (K)
   REAL :: tk
   ! Fraction of aerosol species washed out/rained out/resuspended (0-1)
   REAL :: washfrac, rainfrac, resusp_frac
   ! First-order rainout rate (s-1)
   REAL :: k_rain
   ! Value of chem(i,k,j,p) before deposition step (ug/kg)
   REAL :: chem_spec_before
   !TODO MOSAIC aerosol stuff
   INTEGER :: itype,iphase,ispec_mosaic
   REAL :: aer_dryvol_before_scav,aer_dryvol_after_scav ! mm3/kg dry air
   REAL, DIMENSION(num_chem) :: chem_species

   itype = 1

   !-------- Perform wet scavenging for aerosol species -------
   DO j=jts,jte
   DO i=its,ite
     fr = 0.0
     deposition_flux(:, :) = 0.0

     DO k=kte-1,kts,-1 ! Do not change, needs to loop in descending order
       !---- Aerosol wet scavenging - initialize parameters for rainout and washout schemes
       tk = t_phy(i, k, j)
       rainrate = rainprod(i,k,j)
       evaprate = evapprod(i,k,j)
       ! Total precipitation through grid cell in kg/m2/hr
       ppr = precr(i,k,j) * 3600.
       pps = precs(i,k,j) * 3600.
       ppi = preci(i,k,j) * 3600.
       ppr = MAX(ppr, 0.)
       pps = MAX(pps, 0.)
       ppi = MAX(ppi, 0.)
       ! Add evaporated precip water back to precip_water because evapprod was
       ! substracted just before in previous mp call
       condensed_water = qc(i,k,j) + qi(i,k,j) + qr(i,k,j) + qs(i,k,j) + evaprate*dtstep
       precip_water = qr(i,k,j) + qs(i,k,j) + evaprate*dtstep
       IF(condensed_water .GT. 1.0E-12) THEN
         liquid_fraction = MIN(1.0, MAX(0.0, (qc(i,k,j)+qr(i,k,j)) / condensed_water))
       ELSE
         liquid_fraction = 1.0
       ENDIF
       ! Calculate f, the precipitating fraction of the grid cell as defined
       ! in Luo et al. (2020)
       ! Also calculate rainfrac, the fraction of rained out aerosols
       IF(cldfra(i, k, j) .GT. 0.001 .AND. condensed_water .GT. 1.0E-12 &
          .AND. rainrate .GT. 0.0) THEN
         CALL rainfrac_aerosol(rainfrac, f, rainrate, dtstep, condensed_water)
       ELSE
         rainfrac = 0.0
         f = 0.0
       END IF
       ! The grid cell fraction experiencing washout, fr, is the max between f
       ! at k and fr at k+1 (this is why we loop over k in descending order)
       fr = MIN(1.0, MAX(fr, f))

       !---- Aerosol wet scavenging - loop over aerosol sizes and perform rainout + washout
       !TODO re-check that this size loop does not cause issues of carried over
       !values that should be reinitialized at each loop step, e.g. deposition
       !flux
       DO isize = 1,naersizes
         ! For MOSAIC aerosols, calculate total aerosol dry volume mixing
         ! ratio before wet scavenging.
         ! This is necessary in order to calculate number changes and aerosol
         ! changes at the end of the isize loop
         IF(aer_is_mosaic) THEN
           chem_species(:) = chem(i,k,j,:)
           CALL calc_mosaic_aer_dryvol(aer_dryvol_before_scav, chem_species, isize, itype)
         ENDIF ! aer_is_mosaic

         !-- Aerosol wet scavenging for isize - perform rainout+washout for all
         ! species of size isize and species ispec
         DO ispec = 1,naerspecies
           p_chem = aer_pointers(ispec,isize)

           ! Skip wet scavenging if aerosol (ispec,isize) is undefined
           !TODO this is needed because some species are not defined for all size bins
           IF(p_chem .GE. PARAM_FIRST_SCALAR) THEN

             !---- Rainout ("in-cloud scavenging") for aerosol (ispec,isize)
             IF(f .GT. 0 .AND. cldfra(i,k,j) .GT. 0.001) THEN

               !-- Aerosol rainout (ispec, isize) - Calculate the mass concentration of cloudborne aerosols
               !TODO Can take this out of chem loop, the only thing that changes here is the aq_aer_ratio_ispec and
               !chem(i,k,j,p_chem) that could be applied in the end to active_cloud_water_frac, active_cloud_ice_frac to get
               !cloudborne_aer (TODO could also be put into subroutine I think)
               ! REAL :: clfrac ! cldfra at i,k,j
               ! REAL :: active_cloud_water_frac, active_cloud_ice_frac
               ! In init:
               ! clfrac = cldfra(i,k,j)
               !
               ! Before the isize,ispec loop, careful to include cldfra checks
               ! active_cloud_water_frac = 0.0
               ! active_cloud_ice_frac = 0.0
               ! CALL get_active_cloud_frac(active_cloud_water_frac, active_cloud_ice_frac,
               !                            clfrac, liquid_fraction)
               !TODO make sure this does not sum to more than 100%
               ! cloudborne_aer = chem(i,k,j,p_chem) * aq_aer_ratio_ispec * active_cloud_water_frac &
               !                  + chem(i,k,j,p_chem) * ice_aer_ratio_ispec * active_cloud_ice_frac
               cloudborne_aer = 0.0
               ! Liquid-phase cloudborne aerosols
               aq_aer_ratio_ispec = aq_aer_ratio(ispec,isize)
               IF(tk .GT. 237.0 .AND. aq_aer_ratio_ispec .GT. 0.0) THEN
                 cloudborne_aer = chem(i,k,j,p_chem) * aq_aer_ratio_ispec &
                                  * cldfra(i,k,j)
                 ! In mixed-phase clouds, correct by the liquid water fraction
                 IF(tk .LT. 258.0 .AND. liquid_fraction .GT. 0.0) THEN
                   cloudborne_aer = MIN(cloudborne_aer,  &
                                    chem(i,k,j,p_chem) * aq_aer_ratio_ispec &
                                    * cldfra(i,k,j) * liquid_fraction)
                 ENDIF
               ENDIF ! aq_aer_ratio_ispec .GT. 0.0 .AND. tk .GT. 237.0
               ! Ice-phase cloudborne aerosols
               ice_aer_ratio_ispec = ice_aer_ratio(ispec,isize)
               IF(tk .LT. 258.0 .AND. ice_aer_ratio_ispec .GT. 0.0) THEN
                 ! Reduce ice-phase scavenging efficiency at high temperatures
                 IF(tk .GT. 237.0) THEN
                   ice_aer_ratio_ispec = ice_aer_ratio_ispec &
                                      * EXP(0.46*(273.16-tk)-11.6)/153.5
                 ENDIF
                 cloudborne_aer = chem(i,k,j,p_chem) * ice_aer_ratio_ispec * cldfra(i,k,j)
                 ! Correct by the ice fraction (1.0-liquid_fraction) at high
                 ! temperatures
                 IF(tk .GT. 237.0) THEN
                   cloudborne_aer = MIN(cloudborne_aer, chem(i,k,j,p_chem) &
                                    * ice_aer_ratio_ispec * cldfra(i,k,j) &
                                    * (1.0 - liquid_fraction))
                 ENDIF
               ENDIF ! ice_aer_ratio_ispec .GT. 0.0

               !-- Appply rainout to chem
               IF(cloudborne_aer .GT. 0.0) THEN
                 rainout_aer = rainfrac * cloudborne_aer
                 chem_spec_before = chem(i,k,j,p_chem)
                 chem(i,k,j,p_chem) = MAX(0.0, chem(i,k,j,p_chem) - rainout_aer)
                 deposition_flux(ispec,isize) = deposition_flux(ispec,isize) + (chem_spec_before - chem(i,k,j,p_chem)) &
                                       * rho_phy(i,k,j) * dz8w(i,k,j)
               ENDIF
             ENDIF ! f .GT. 0 .AND. cldfra(i, k, j) .GT. 0.001


             !---- Perform washout ("below-cloud scavenging") for species (ispec,isize)
             IF(ppr+pps+ppi .GT. 0.0 .AND. fr .GT. 0.001) THEN
               !-- Calculate washfrac, the fraction of washed out aerosols
               ! Fine aerosol washout (D<=2.5 micrometer)
               IF(aer_is_fine(ispec,isize)) THEN
                  CALL washfrac_fine_aerosol(washfrac, ppr, pps, ppi, tk, fr, dtstep)
               ! Coarse aerosol washout (D>2.5 micrometer)
               ELSEIF(aer_is_coarse(ispec,isize)) THEN
                  CALL washfrac_coarse_aerosol(washfrac, ppr, pps, ppi, tk, fr, dtstep)
               ! Any other aerosol is not washed out
               ELSE
                 washfrac = 0.0
               END IF
               !-- Apply washout to chem
               IF(washfrac .GT. 0.0) THEN
                 chem_spec_before = chem(i,k,j,p_chem)
                 chem(i,k,j,p_chem) = MAX(0.0, (1.0-washfrac) * chem(i,k,j,p_chem))
                 deposition_flux(ispec,isize) = deposition_flux(ispec,isize) + (chem_spec_before - chem(i,k,j,p_chem)) &
                                       * rho_phy(i,k,j) * dz8w(i,k,j)
               ENDIF
             END IF ! ppr+pps+ppi .GT. 0 .AND. fr .GT. 0.01


             !---- Perform resuspension from evaporating precip. for species ispec,isize
             ! We apply resuspension to the total deposition flux from the column
             ! above, so this assumes that the timescale of precipitation to level
             ! k is smaller than dtstep
             IF(precip_water .GT. 0. .AND. evaprate .GT. 0.) THEN
               ! Scale by 0.5 (TODO Luo ref), to take into account that
               ! evaporation of hydrometeors can be incomplete
               resusp_frac = 0.5 * evaprate * dtstep / precip_water
               resusp_frac = MAX(0.0, MIN(1.0, resusp_frac))
               resusp_aer = (resusp_frac * deposition_flux(ispec,isize) / (rho_phy(i,k,j) * dz8w(i,k,j)))
               chem(i,k,j,p_chem) = MAX(0.0, chem(i,k,j,p_chem) + resusp_aer)
               deposition_flux(ispec,isize) = deposition_flux(ispec,isize) * (1.0 - resusp_frac)
             END IF

          END IF ! p_chem .GE. PARAM_FIRST_SCALAR

         END DO ! ispec

         !TODO recalculate aer dry volume for isize, calculate deposited
         ! ratio, apply to num and water for isize
         !-- Aerosol wet scavenging - calculate MOSAIC aerosol parameters for bin isize after wet
         ! scav
         IF(aer_is_mosaic) THEN
           chem_species(:) = chem(i,k,j,:)
           CALL calc_mosaic_aer_dryvol(aer_dryvol_after_scav, chem_species, isize, itype)
           !TODO need to check that num is still within expected bounds at the
           !end (something went wrong if not) - check total_vol_bin_min and
           !total_vol_bin_max = num * 4./3. * pi * (dhi/2.)**3. (same with dlo,
           !careful with underflow with **3.0, convert units for everything if
           !possible) TODO or check mosaic settling for this
           ! WRITE(*,*) (aer_dryvol_after_scav-aer_dryvol_before_scav)/aer_dryvol_before_scav
           ! TODO not sure what to do if resusp brings volume from ~0 to >0
           ! TODO handle resuspension, maybe need to move around water in
           ! deposition_flux (aer_water_dep_flux, aer_dep_flux) but will have
           ! issues with cw species (remove them later if this is used)
         ENDIF ! aer_is_mosaic

       END DO ! isize
     END DO ! k
   END DO ! i
   END DO ! j

   END SUBROUTINE aer_wetscav_simple

!----------------------------------------------------------------------------

   SUBROUTINE rainfrac_aerosol(rainfrac, f, rainrate, dtstep, condensed_water)

     ! Purpose:
     ! Calculates rainfrac, the fraction of cloudborne aerosols lost to rainout.
     !
     ! Arguments:
     ! Input:
     !   rainrate - precipitation formation rate (kg/kg/s)
     !   dtstep - time step (s)
     !   condensed_water - condensed water mixing ratio (kg/kg)
     ! Output:
     !   rainfrac - fraction of aerosols lost to rainout (unitless)
     !   f - rainout area fraction (unitless)
     !
     ! Remarks:
     ! From GEOS-Chem, Luo et al (2020)

     IMPLICIT NONE

     ! Arguments
     REAL, INTENT(IN) :: rainrate
     REAL, INTENT(IN) :: dtstep
     REAL, INTENT(INOUT) :: condensed_water
     REAL, INTENT(OUT) :: rainfrac
     REAL, INTENT(OUT) :: f
     ! Local variables and parameters
     REAL :: k_rain ! first-order rainout rate (s-1)
     REAL,  PARAMETER :: k_0 = 1.0E-4 ! s-1

     ! Do not allow more than 100% of the cloud water to precipitate
     condensed_water = MAX(rainrate*dtstep, condensed_water)
     k_rain = k_0 + rainrate / condensed_water
     ! Area of grid box undergoing precipitation
     f = rainrate / (k_rain * condensed_water)
     f = MIN(f, 1.0)
     ! Fraction of aerosols lost to rainout
     rainfrac = f * (1.0 - EXP(-k_rain*dtstep))

   END SUBROUTINE rainfrac_aerosol

!----------------------------------------------------------------------------

   SUBROUTINE washfrac_fine_aerosol(washfrac, ppr, pps, ppi, tk, fr, dtstep)

     ! Purpose:
     ! Calculates washfrac, the fraction of fine (r_eff < 1 µm) aerosols
     ! lost to washout.
     !
     ! Arguments:
     ! Input:
     !   ppr = rain precipitation through cell (kg/m2/hr = mm/hr)
     !   pps = snow precipitation through cell (kg/m2/hr = mm/hr)
     !   ppi = ice precipitation through cell (kg/m2/hr = mm/hr)
     !   tk = temperature (K)
     !   fr = fraction of cell experiencing washout
     !   dtstep - time step (s)
     ! Output:
     !   washfrac - fraction of aerosols lost to washout (unitless)
     !
     ! Comments:
     ! From GEOS-Chem, Luo et al (2020), using GC12 parameters

     IMPLICIT NONE

     ! Arguments
     REAL, INTENT(IN) :: ppr
     REAL, INTENT(IN) :: pps
     REAL, INTENT(IN) :: ppi
     REAL, INTENT(IN) :: tk
     REAL, INTENT(IN) :: fr
     REAL, INTENT(IN) :: dtstep
     REAL, INTENT(OUT) :: washfrac

     ! Use the parameteters from the dominant precipitation type
     washfrac = 0.0
     IF(ppr >= (pps+ppi)) THEN
     ! Mostly rain
       washfrac = fr * (1.0 - EXP(-4.3E-6 * ((ppr+pps+ppi)/fr)**0.61 * dtstep))
     ELSEIF(pps >= (ppr+ppi)) THEN
     ! Mostly snow
       washfrac = fr * (1.0 - EXP(-8.8E-6 * ((ppr+pps+ppi)/fr)**0.96 * dtstep))
     ELSEIF(ppi >= (ppr+pps)) THEN
     ! Mostly ice
       washfrac = fr * (1.0 - EXP(-8.8E-6/5.0 * ((ppr+pps+ppi)/fr)**0.96 * dtstep))
     ELSE
     ! If mixed-phase precip with no clear main phase, default to ice
       washfrac = fr * (1.0 - EXP(-8.8E-6/5.0 * ((ppr+pps+ppi)/fr)**0.96 * dtstep))
     END IF

   END SUBROUTINE washfrac_fine_aerosol

!----------------------------------------------------------------------------

   SUBROUTINE washfrac_coarse_aerosol(washfrac, ppr, pps, ppi, tk, fr, dtstep)

     ! Purpose:
     ! Calculates washfrac, the fraction of coarse (r_eff > 1 µm) aerosols lost
     ! to washout.
     !
     ! Arguments:
     ! Input:
     !   ppr = rain precipitation through cell (kg/m2/hr = mm/hr)
     !   pps = snow precipitation through cell (kg/m2/hr = mm/hr)
     !   ppi = ice precipitation through cell (kg/m2/hr = mm/hr)
     !   tk = temperature (K)
     !   fr = fraction of cell experiencing washout
     !   dtstep - time step (s)
     ! Output:
     !   washfrac - fraction of aerosols lost to washout (unitless)
     !
     ! Comments:
     ! From GEOS-Chem, Luo et al (2020), using GC12 parameters

     IMPLICIT NONE

     ! Arguments
     REAL, INTENT(IN) :: ppr
     REAL, INTENT(IN) :: pps
     REAL, INTENT(IN) :: ppi
     REAL, INTENT(IN) :: tk
     REAL, INTENT(IN) :: fr
     REAL, INTENT(IN) :: dtstep
     REAL, INTENT(OUT) :: washfrac

     ! Use the parameteters from the dominant precipitation type
     washfrac = 0.0
     IF(ppr >= (pps+ppi)) THEN
     ! Mostly rain
       washfrac = fr * (1.0 - EXP(-2.6E-4 * ((ppr+pps+ppi)/fr)**0.79 * dtstep))
     ELSEIF(pps >= (ppr+ppi)) THEN
     ! Mostly snow
       washfrac = fr * (1.0 - EXP(-4.2E-4 * ((ppr+pps+ppi)/fr)**0.96 * dtstep))
     ELSEIF(ppi >= (ppr+pps)) THEN
     ! Mostly ice
       washfrac = fr * (1.0 - EXP(-4.2E-4/5.0 * ((ppr+pps+ppi)/fr)**0.96 * dtstep))
     ELSE
     ! If mixed-phase precip with no clear main phase, default to ice
       washfrac = fr * (1.0 - EXP(-4.2E-4/5.0 * ((ppr+pps+ppi)/fr)**0.96 * dtstep))
     END IF

   END SUBROUTINE washfrac_coarse_aerosol

!----------------------------------------------------------------------------

   SUBROUTINE  calc_mosaic_aer_dryvol(aer_dryvol, chem, isize, itype)
     ! Purpose
     !  Calculate the total aerosol dry volume for MOSAIC aerosol bin isize
     !
     !TODO arguments, remarks
     !TODO add depend.chem link to module_data_mosaic_asect
     
     USE module_state_description, ONLY: PARAM_FIRST_SCALAR, num_chem
     USE module_data_mosaic_asect, ONLY: dens_aer, ncomp_aer, massptr_aer

     IMPLICIT NONE

     ! Arguments
     !  Input  TODO
     INTEGER, INTENT(IN) :: isize,itype
     REAL, DIMENSION(num_chem), INTENT(IN) :: chem_species
     REAL, INTENT(OUT) :: aer_dryvol

     ! Local variables TODO
     INTEGER :: iphase,p_chem

     iphase = 1
     aer_dryvol = 0.0
     DO ispec_mosaic = 1, ncomp_aer(itype)
       p_chem = massptr_aer(ispec_mosaic,isize,itype,iphase)
       IF(p_chem .GE. PARAM_FIRST_SCALAR) THEN
         ! ug/kg / mg/cm3 -> 1.0E-3 cm3/kg = mm3/kg dry air
         ! We don't convert to SI units because values would be too small
         aer_dryvol = aer_dryvol + chem_species(p_chem) &
                     / (dens_aer(ispec_mosaic,itype) * 1.0E3)
       ENDIF
     ENDDO

   END SUBROUTINE calc_mosaic_aer_dryvol

!----------------------------------------------------------------------------

   SUBROUTINE aer_wetscav_simple_init(chem_opt)

     ! Purpose
     !  Initialize aer_wetscav_simple parameters, which depend on chem_opt
     !
     ! Arguments
     !  chem_opt TODO
     !TODO arguments, remarks
     !
     ! Remarks
     !TODO arguments, remarks

     USE module_state_description ! For chem_opt option values

     IMPLICIT NONE

     ! Arguments
     !  Input
     !TODO
     INTEGER, INTENT(IN) :: chem_opt

     ! Local variables
     ! Aeorosol size loop index
     INTEGER :: isize

    !---- Initialize chem_opt aerosol ID
    aer_is_mosaic = chem_opt == CBMZ_MOSAIC_4BIN .OR. &
                    chem_opt == CBMZ_MOSAIC_DMS_4BIN .OR. &
                    chem_opt == CBMZ_MOSAIC_4BIN_AQ .OR. &
                    chem_opt == CBMZ_MOSAIC_DMS_4BIN_AQ .OR. &
                    chem_opt == MOZART_MOSAIC_4BIN_KPP .OR. &
                    chem_opt == MOZART_MOSAIC_4BIN_AQ_KPP .OR. &
                    chem_opt == CRI_MOSAIC_4BIN_AQ_KPP .OR. &
                    chem_opt == SAPRC99_MOSAIC_4BIN_VBS2_KPP .OR. &
                    chem_opt == CBMZ_MOSAIC_8BIN .OR. &
                    chem_opt == CBMZ_MOSAIC_DMS_8BIN .OR. &
                    chem_opt == CBMZ_MOSAIC_8BIN_AQ .OR. &
                    chem_opt == CBMZ_MOSAIC_DMS_8BIN_AQ .OR. &
                    chem_opt == CRI_MOSAIC_8BIN_AQ_KPP .OR. &
                    chem_opt == SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP .OR. &
                    chem_opt == SAPRC99_MOSAIC_8BIN_VBS2_KPP

    !---- Initialize aerosol tracer number for chem_opt
    ! This could be recovered from the aerosol scheme data directly
    IF(chem_opt == MOZCART_KPP) THEN
      ! MOZCART
      naerspecies = 12
      naersizes = 1
    ELSEIF(chem_opt == MOZART_MOSAIC_4BIN_AQ_KPP) THEN
      ! MOZART-MOSAIC-4bin-AQ
      naerspecies = 8
      naersizes = 4
    ELSE
      CALL wrf_error_fatal ('wetscav_simple_init: chem_opt not supported')
    ENDIF

    !---- Allocate aerosol mappings and parameters for wet scavenging
    !-- Species-resolved
    IF (.NOT. ALLOCATED(aer_pointers)) THEN
      ALLOCATE(aer_pointers(naerspecies, naersizes))
    ENDIF
    IF (.NOT. ALLOCATED(aq_aer_ratio)) THEN
      ALLOCATE(aq_aer_ratio(naerspecies, naersizes))
    ENDIF
    IF (.NOT. ALLOCATED(ice_aer_ratio)) THEN
      ALLOCATE(ice_aer_ratio(naerspecies, naersizes))
    ENDIF
    IF (.NOT. ALLOCATED(aer_is_fine)) THEN
      ALLOCATE(aer_is_fine(naerspecies, naersizes))
    ENDIF
    IF (.NOT. ALLOCATED(aer_is_coarse)) THEN
      ALLOCATE(aer_is_coarse(naerspecies, naersizes))
    ENDIF
    IF (.NOT. ALLOCATED(deposition_flux)) THEN
      ALLOCATE(deposition_flux(naerspecies, naersizes))
    ENDIF
    !-- Size-resolved
    !TODO maybe replace is_coarse, is_fine by spec center diameter, and only
    !calculate number change if p_num_* > p1st and bin_cen_diam == spec_cen_diam(ispec)
    !
    !TODO aer_densities
    !TODO aer_num_pointers, aer_water_pointers, aer_hysw_pointers + allocate +
    ! define in init

    !---- Initialize aerosol mappings and parameters for wet scavenging
    deposition_flux(:, :) = 0.0

    IF(chem_opt == MOZCART_KPP) THEN
      !-- MOZCART
      ! bin_diameters_cen = (/ 0.0 /) ! Placeholder since this is unused for
      !                                ! single moment schemes
      ! Species that experience wet scavenging (sulfate deposition needs to
      ! be disabled in mozcart_wetscav)
      aer_pointers(:, 1) = (/ p_sulf, p_seas_1, p_seas_2, p_seas_3, &
                        p_seas_4, p_bc2, p_oc2, p_dust_1, &
                        p_dust_2, p_dust_3, p_dust_4, p_bc1 /)
      ! Fraction of activated aer_pointers species (Luo et al. 2020)
      aq_aer_ratio(:, 1) = (/ 1., 1., 1., 1., 1., 0.5, 0.5, 0., 0., 0., 0., 0. /)
      ice_aer_ratio(:, 1) = (/ 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 0.5 /)
      ! List of washed out fine aer_pointers species (D<=2.5 micrometer)
      aer_is_fine(:, 1) = (/ .true., .true., .true., .false., &
                       .false., .true., .true. , .true., &
                       .false., .false., .false., .false. /)
      ! List of washed out coarse aer_pointers species (D>2.5 micrometer)
      aer_is_coarse(:, 1) = (/ .false., .false., .false., .true., &
                       .true., .false., .false. , .false., &
                       .true., .true., .true., .false. /)

    ELSEIF(chem_opt == MOZART_MOSAIC_4BIN_AQ_KPP) THEN
      ! so4_a01,no3_a01,asoaX_a01,asoa1_a01,asoa2_a01,asoa3_a01,asoa4_a01,bsoaX_a01,bsoa1_a01,bsoa2_a01,bsoa3_a01,bsoa4_a01,glysoa_r1_a01,glysoa_r2_a01,glysoa_sfc_a01,glysoa_nh4_a01,glysoa_oh_a01,cl_a01,co3_a01,nh4_a01,na_a01,ca_a01,oin_a01,oc_a01,bc_a01,hysw_a01,water_a01,num_a01
      ! so4_a02,no3_a02,asoaX_a02,asoa1_a02,asoa2_a02,asoa3_a02,asoa4_a02,bsoaX_a02,bsoa1_a02,bsoa2_a02,bsoa3_a02,bsoa4_a02,glysoa_r1_a02,glysoa_r2_a02,glysoa_sfc_a02,glysoa_nh4_a02,glysoa_oh_a02,cl_a02,co3_a02,nh4_a02,na_a02,ca_a02,oin_a02,oc_a02,bc_a02,hysw_a02,water_a02,num_a02
      ! so4_a03,no3_a03,asoaX_a03,asoa1_a03,asoa2_a03,asoa3_a03,asoa4_a03,bsoaX_a03,bsoa1_a03,bsoa2_a03,bsoa3_a03,bsoa4_a03,glysoa_r1_a03,glysoa_r2_a03,glysoa_sfc_a03,glysoa_nh4_a03,glysoa_oh_a03,cl_a03,co3_a03,nh4_a03,na_a03,ca_a03,oin_a03,oc_a03,bc_a03,hysw_a03,water_a03,num_a03
      ! so4_a04,no3_a04,asoaX_a04,asoa1_a04,asoa2_a04,asoa3_a04,asoa4_a04,bsoaX_a04,bsoa1_a04,bsoa2_a04,bsoa3_a04,bsoa4_a04,glysoa_r1_a04,glysoa_r2_a04,glysoa_sfc_a04,glysoa_nh4_a04,glysoa_oh_a04,cl_a04,co3_a04,nh4_a04,na_a04,ca_a04,oin_a04,oc_a04,bc_a04,hysw_a04,water_a04,num_a04
      ! so4_cw01,no3_cw01,asoaX_cw01,asoa1_cw01,asoa2_cw01,asoa3_cw01,asoa4_cw01,bsoaX_cw01,bsoa1_cw01,bsoa2_cw01,bsoa3_cw01,bsoa4_cw01,glysoa_r1_cw01,glysoa_r2_cw01,glysoa_sfc_cw01,glysoa_nh4_cw01,glysoa_oh_cw01,cl_cw01,co3_cw01,nh4_cw01,na_cw01,ca_cw01,oin_cw01,oc_cw01,bc_cw01,num_cw01
      ! so4_cw02,no3_cw02,asoaX_cw02,asoa1_cw02,asoa2_cw02,asoa3_cw02,asoa4_cw02,bsoaX_cw02,bsoa1_cw02,bsoa2_cw02,bsoa3_cw02,bsoa4_cw02,glysoa_r1_cw02,glysoa_r2_cw02,glysoa_sfc_cw02,glysoa_nh4_cw02,glysoa_oh_cw02,cl_cw02,co3_cw02,nh4_cw02,na_cw02,ca_cw02,oin_cw02,oc_cw02,bc_cw02,num_cw02
      ! so4_cw03,no3_cw03,asoaX_cw03,asoa1_cw03,asoa2_cw03,asoa3_cw03,asoa4_cw03,bsoaX_cw03,bsoa1_cw03,bsoa2_cw03,bsoa3_cw03,bsoa4_cw03,glysoa_r1_cw03,glysoa_r2_cw03,glysoa_sfc_cw03,glysoa_nh4_cw03,glysoa_oh_cw03,cl_cw03,co3_cw03,nh4_cw03,na_cw03,ca_cw03,oin_cw03,oc_cw03,bc_cw03,num_cw03
      ! so4_cw04,no3_cw04,asoaX_cw04,asoa1_cw04,asoa2_cw04,asoa3_cw04,asoa4_cw04,bsoaX_cw04,bsoa1_cw04,bsoa2_cw04,bsoa3_cw04,bsoa4_cw04,glysoa_r1_cw04,glysoa_r2_cw04,glysoa_sfc_cw04,glysoa_nh4_cw04,glysoa_oh_cw04,cl_cw04,co3_cw04,nh4_cw04,na_cw04,ca_cw04,oin_cw04,oc_cw04,bc_cw04,num_cw04
      !-- MOZART-MOSAIC-4bin-AQ
      ! MOSAIC-4bin bin center diameters (micrometer)
      ! bin_diameters_cen = (/ 0.0781, 0.3125, 1.2500, 5.0000 /)
      ! Species that experience wet scavenging
      !TODO include all species
      !TODO could be automated, only set main species common to all mechs, then set parameters automatically based on size and kappa
      !TODO in this case important to print to rsl the final table with wrf_message
      aer_pointers(:, 1) = (/ p_so4_a01, p_no3_a01, p_nh4_a01, p_na_a01, &
                              p_cl_a01, p_oc_a01, p_bc_a01, p_oin_a01 /)
      aer_pointers(:, 2) = (/ p_so4_a02, p_no3_a02, p_nh4_a02, p_na_a02, &
                              p_cl_a02, p_oc_a02, p_bc_a02, p_oin_a02 /)
      aer_pointers(:, 3) = (/ p_so4_a03, p_no3_a03, p_nh4_a03, p_na_a03, &
                              p_cl_a03, p_oc_a03, p_bc_a03, p_oin_a03 /)
      aer_pointers(:, 4) = (/ p_so4_a04, p_no3_a04, p_nh4_a04, p_na_a04, &
                              p_cl_a04, p_oc_a04, p_bc_a04, p_oin_a04 /)
      ! Fraction of activated aer_pointers species (Luo et al. 2020)
      !TODO define accross cross sections? all s04,nh4,no3,na,cl species
      ! should be 1., all OA & BC species should be 0.5, all oin species should
      ! be 0? TODO leave this guide to mapping on in the module init remarks before
      ! commit
      DO isize = 1,naersizes
        aq_aer_ratio(:, isize) = (/ 1., 1., 1., 1., &
                                    1., 0.25, 0.25, 0./)
        ice_aer_ratio(:, isize) = (/ 0., 0., 0., 0., &
                                     0., 0.25, 0.25, 1. /)
      ENDDO
      ! List of washed out fine aer_pointers species (D<=2.5 micrometer)
      !TODO could be automated based on bin_diameters_cen
      aer_is_fine(:, 1) = .true.
      aer_is_fine(:, 2) = .true.
      aer_is_fine(:, 3) = .true.
      aer_is_fine(:, 4) = .false.
      ! List of washed out coarse aer_pointers species (D>2.5 micrometer)
      aer_is_coarse(:, 1) = .false.
      aer_is_coarse(:, 2) = .false.
      aer_is_coarse(:, 3) = .false.
      aer_is_coarse(:, 4) = .true.

    ELSE
      CALL wrf_error_fatal ('wetscav_simple_init: chem_opt not supported')

    ENDIF ! chem_opt

   END SUBROUTINE aer_wetscav_simple_init

!----------------------------------------------------------------------------

END MODULE module_aer_wetscav_simple


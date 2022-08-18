! 2021/09/22 Louis Marelle for CICERO
!
! Purpose:
! Perform wet scavenging of GOCART aerosols by grid scale rainout/washout
!
! Remarks:
! Washout and rainout from Luo et al. (2019, 2020)
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
!TODO also perform resuspension when rain evaporates (assume
!precipitation to ground timescale is < dtstep)
!TODO get wet removal fluxes out for validation
!TODO use reasonable thresholds instead of 0.0 in all the module
!TODO remove unused arguments, e.g. dx, dy, config_flags, rho_phy, p_phy
!TODO maybe use double precision
!

MODULE module_gocart_wetscav

   ! USE module_wrf_error, ONLY : wrf_debug

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: wetscav_gocart

CONTAINS

!------------------------------------------------------------------------------

   SUBROUTINE wetscav_gocart( id, dtstep, config_flags,                    &
                              t_phy, p_phy,                                &
                              chem, rho_phy, cldfra, rainprod, evapprod,   &
                              precr, preci, precs, precg, dx, dy,          &
                              qv, qc, qi, qs, qr,                          &
                              ids,ide, jds,jde, kds,kde,                   &
                              ims,ime, jms,jme, kms,kme,                   &
                              its,ite, jts,jte, kts,kte                    )

   ! Purpose:
   ! Calculates and removes wet scavenging tendencies for GOCART aerosols
   !
   ! Dummy arguments:
   ! Input:
   !   id = domain ID
   !   dtstep = model time step (s)
   !   config_flags = namelist configuration flags
   !   t_phy = 3D physical temperature (K)
   !   p_phy = 3D pressure (Pa)
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
   !   chem = advected chemical species (gases in ppm, aerosols in µg/kg)

   ! USE associations:
   !TODO Use ONLY:
   USE module_configure
   USE module_state_description
   ! p_seas_1, p_seas_2, p_seas_3, p_seas_4, p_bc2, p_bc1, p_oc2, p_dust_1,
   ! p_dust_2, p_dust_3, p_dust_4, p_sulf from module_state_description
   ! USE module_data_gocart_seas, only: reff_seas
   ! USE module_data_gocart_dust, only: reff_dust

   !---- Dummy arguments
   TYPE(grid_config_rec_type), INTENT(IN) :: config_flags
   INTEGER, INTENT(IN) :: ids,ide, jds,jde, kds,kde,    &
                          ims,ime, jms,jme, kms,kme,    &
                          its,ite, jts,jte, kts,kte,    &
                          id
   REAL, INTENT(IN) :: dx, dy
   REAL, INTENT(IN) :: dtstep
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme)         ,        &
       INTENT(IN) :: t_phy, p_phy, rho_phy
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
   ! x,y,z, and chemical species index, chemical species pointer
   INTEGER :: i, j, k, nv, p_chem
   ! Cell fractions where rainout/washout is applied (See Luo et al., 2020)
   REAL :: f, fr
   ! Vertical dimension of cldfra for cloud fraction calculation
   INTEGER :: dimz
   ! Aerosol conc. in cloud phase, aerosols conc. lost to rainout(µg/kg)
   REAL :: cloudborne_aer, rainout_aer
   ! Condensed water mixing ratio (kg/kg dry air), liquid fraction of condensed
   ! water (0-1)
   REAL :: condensed_water, liquid_fraction
   ! Rainout activation efficiency of species (0-1)
   REAL :: aq_aer_ratio_nv, ice_aer_ratio_nv
   ! Total 3D precipitation trough cell for rain/snow/ice (kg/m2/hr), rainprod at i,k,j (kg/kg/s)
   REAL :: ppr, pps, ppi, rainrate
   ! Temperature (K)
   REAL :: tk
   ! Fraction of aerosol species washed out/rained out (0-1)
   REAL :: washfrac, rainfrac
   ! First-order rainout rate (s-1)
   REAL :: k_rain
   ! Number of wet-scavenged aerosol species in gocart_pointers
   INTEGER, PARAMETER :: NAERSPECIES = 12
   ! List of gocart species experiencing rainout/washout
   INTEGER :: gocart_pointers(NAERSPECIES)
   ! Rainout activation efficiency of species (0-1)
   REAL :: aq_aer_ratio(NAERSPECIES), ice_aer_ratio(NAERSPECIES)

   !---- Data declaration
   ! GOCART species that experience wet scavenging (sulfate deposition needs to be
   ! disabled in mozcart_wetscav)
   gocart_pointers = (/ p_sulf, p_seas_1, p_seas_2, p_seas_3, &
                        p_seas_4, p_bc2, p_oc2, p_dust_1, &
                        p_dust_2, p_dust_3, p_dust_4, p_bc1 /)
   ! Fraction of activated aerosols, values from Luo et al. (2020)
   aq_aer_ratio = (/ 1., 1., 1., 1., 1., 0.5, 0.5, 0., 0., 0., 0., 0. /)
   ice_aer_ratio = (/ 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 0.5 /)


   !-------- Perform wet scavenging for GOCART species -------
   DO j=jts,jte
     DO i=its,ite
       fr = 0.0

       DO k=kte-1,kts,-1 ! Do not change, needs to loop in descending order
         !---- Calculate parameters for rainout and washout schemes
         tk = t_phy(i, k, j)
         condensed_water = QC(i, k, j) + QI(i, k, j) + QR(i, k, j) + QS(i, k, j)
         IF(condensed_water .GT. 1.0E-12) THEN
           liquid_fraction = MIN(1.0, MAX(0.0, (QC(i, k, j)+QR(i, k, j)) / condensed_water))
         ELSE
           liquid_fraction = 1.0
         ENDIF
         rainrate = MAX(rainprod(i, k, j) - evapprod(i, k, j), 0.0)
         ! Total precipitation through grid cell in kg/m2/hr
         ppr = precr(i, k, j) * 3600.
         pps = precs(i, k, j) * 3600.
         ppi = preci(i, k, j) * 3600.
         ppr = MAX(ppr, 0.)
         pps = MAX(pps, 0.)
         ppi = MAX(ppi, 0.)
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

         DO nv = 1,NAERSPECIES
           p_chem = gocart_pointers(nv)

           !---- Perform rainout ("in-cloud scavenging")
           IF(f .GT. 0 .AND. cldfra(i, k, j) .GT. 0.001) THEN

             !-- Calculate the mass concentration of cloudborne aerosols
             cloudborne_aer = 0.0
             ! Liquid-phase cloudborne aerosols
             aq_aer_ratio_nv = aq_aer_ratio(nv)
             IF(tk .GT. 237.0 .AND. aq_aer_ratio_nv .GT. 0.0) THEN
               cloudborne_aer = chem(i ,k, j, p_chem) * aq_aer_ratio_nv &
                                * cldfra(i, k, j)
               ! In mixed-phase clouds, correct by the liquid water fraction
               IF(tk .LT. 258.0 .AND. liquid_fraction .GT. 0.0) THEN
                 cloudborne_aer = MIN(cloudborne_aer,  &
                                  chem(i ,k, j, p_chem) * aq_aer_ratio_nv &
                                  * cldfra(i, k, j) * liquid_fraction)
               ENDIF
             ENDIF ! aq_aer_ratio_nv .GT. 0.0 .AND. tk .GT. 237.0
             ! Ice-phase cloudborne aerosols
             ice_aer_ratio_nv = ice_aer_ratio(nv)
             IF(tk .LT. 258.0 .AND. ice_aer_ratio_nv .GT. 0.0) THEN
               ! Reduce ice-phase scavenging efficiency at high temperatures
               IF(tk .GT. 237.0) THEN
                 ice_aer_ratio_nv = ice_aer_ratio_nv &
                                    * EXP(0.46*(273.16-tk)-11.6)/153.5
               ENDIF
               cloudborne_aer = chem(i ,k, j, p_chem) * ice_aer_ratio_nv * cldfra(i, k, j)
               ! Correct by the ice fraction (1.0-liquid_fraction) at high
               ! temperatures
               IF(tk .GT. 237.0) THEN
                 cloudborne_aer = MIN(cloudborne_aer, chem(i ,k, j, p_chem) &
                                  * ice_aer_ratio_nv * cldfra(i, k, j) &
                                  * (1.0 - liquid_fraction))
               ENDIF
             ENDIF ! ice_aer_ratio_nv .GT. 0.0

             !-- Appply rainout to chem
             IF(cloudborne_aer .GT. 0.0) THEN
               rainout_aer = rainfrac * cloudborne_aer
               chem(i ,k, j, p_chem) = MAX(0.0, chem(i ,k, j, p_chem) - rainout_aer)
             ENDIF
           ENDIF ! f .GT. 0 .AND. cldfra(i, k, j) .GT. 0.001


           !---- Perform washout ("below-cloud scavenging")
           IF(ppr+pps+ppi .GT. 0.0 .AND. fr .GT. 0.001) THEN
             !-- Calculate washfrac, the fraction of washed out aerosols
             ! Fine aerosol washout (r < 1 µm)
             IF(p_chem == p_bc2 .OR. p_chem == p_oc2 .OR. p_chem == p_seas_1 &
                .OR. p_chem == p_seas_2 .OR. p_chem == p_sulf &
                .OR. p_chem == p_dust_1) THEN
                CALL washfrac_fine_aerosol(washfrac, ppr, pps, ppi, tk, fr, dtstep)
             ! Coarse aerosol washout (r > 1 µm)
             ELSEIF(p_chem == p_seas_3 .OR. p_chem == p_seas_4 &
                    .OR. p_chem == p_dust_2 .OR. p_chem == p_dust_3 &
                    .OR. p_chem == p_dust_4) THEN
                CALL washfrac_coarse_aerosol(washfrac, ppr, pps, ppi, tk, fr, dtstep)
             ! Any other aerosol is not washed out
             ELSE
               washfrac = 0.0
             END IF
             !-- Apply washout to chem
             IF(washfrac .GT. 0.0) THEN
               chem(i ,k, j, p_chem) = MAX(0.0, (1.0 - washfrac)*chem(i ,k, j, p_chem))
             ENDIF
           END IF ! ppr+pps+ppi .GT. 0 .AND. fr .GT. 0.01

         END DO ! nv
       END DO ! k
     END DO ! i
   END DO ! j

   END SUBROUTINE wetscav_gocart

!------------------------------------------------------------------------------

   SUBROUTINE rainfrac_aerosol(rainfrac, f, rainrate, dtstep, condensed_water)

     ! Purpose:
     ! Calculates rainfrac, the fraction of cloudborne aerosols lost to rainout.
     !
     ! Dummy arguments:
     ! Input:
     !   rainrate - preciptation formation rate (kg/kg/s)
     !   dtstep - time step (s)
     !   condensed_water - condensed water mixing ratio (kg/kg)
     ! Output:
     !   rainfrac - fraction of aerosols lost to rainout (unitless)
     !   f - rainout area fraction (unitless)
     !
     ! Remarks:
     ! From GEOS-Chem, Luo et al (2020)

     IMPLICIT NONE

     ! Dummy arguments
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

!------------------------------------------------------------------------------

   SUBROUTINE washfrac_fine_aerosol(washfrac, ppr, pps, ppi, tk, fr, dtstep)

     ! Purpose:
     ! Calculates washfrac, the fraction of fine (r_eff < 1 µm) aerosols
     ! lost to washout.
     !
     ! Dummy arguments:
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

     ! Dummy arguments
     REAL, INTENT(IN) :: ppr
     REAL, INTENT(IN) :: pps
     REAL, INTENT(IN) :: ppi
     REAL, INTENT(IN) :: tk
     REAL, INTENT(IN) :: fr
     REAL, INTENT(IN) :: dtstep
     REAL, INTENT(OUT) :: washfrac

     ! Use the parameteters from the dominant precipitation type
     washfrac = 0.0
     !TODO could be rewritten as a function, changing only parameters
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

!------------------------------------------------------------------------------

   SUBROUTINE washfrac_coarse_aerosol(washfrac, ppr, pps, ppi, tk, fr, dtstep)

     ! Purpose:
     ! Calculates washfrac, the fraction of coarse (r_eff > 1 µm) aerosols lost
     ! to washout.
     !
     ! Dummy arguments:
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

     ! Dummy arguments
     REAL, INTENT(IN) :: ppr
     REAL, INTENT(IN) :: pps
     REAL, INTENT(IN) :: ppi
     REAL, INTENT(IN) :: tk
     REAL, INTENT(IN) :: fr
     REAL, INTENT(IN) :: dtstep
     REAL, INTENT(OUT) :: washfrac

     ! Use the parameteters from the dominant precipitation type
     washfrac = 0.0
     !TODO could be rewritten as a function, changing only parameters
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

!------------------------------------------------------------------------------

END MODULE module_gocart_wetscav


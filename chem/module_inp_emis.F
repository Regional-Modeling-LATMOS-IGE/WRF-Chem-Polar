! Anderson Da Silva, Louis Marelle, 2026
!
! Purpose:
!  Prognostic ice nucleating particle (INP) scheme for WRF-Chem-Polar. INP number
!  concentrations are created from the model aerosol ONLY at the moment of
!  emission (species-specific ice-active parameterization), then carried as
!  independent advected chem tracers and removed during transport by dry
!  deposition, wet deposition, and in-cloud activation. Intended for comparison
!  with in-situ immersion-freezing observations (e.g. MOSAiC, MOCCHA, ARTofMelt,
!  Tara Polar Station) and, later, for coupling to the microphysics ice nucleation.
!
!  Design ("emit-once"): turning aerosol into INP at emission (not continuously
!  during transport) makes the INP genuinely independent of the bulk aerosol
!  afterwards. The temperature-selective activation sink then permanently depletes
!  them -- cloud-processed air stays INP-poor -- with NO temperature-history
!  variable and NO artificial regeneration. Fresh emission is the only INP source.
!
!  Storage: DIFFERENTIAL activation-temperature bins (# kg-1, immersion mode), per
!  species (dust, marine sea-spray, biological PBAP). The cumulative spectrum
!  n_INP(T) is the partial sum of the bins from bin 1 to the requested temperature.
!    bin 1: freezing temperature      T >= -10 C
!    bin 2:                     -10 >  T >= -15 C
!    bin 3:                     -15 >  T >= -20 C
!    bin 4:                     -20 >  T >= -25 C
!  Cumulative TOTAL n_INP(T) (dust+marine+bio) is written to the diagnostic
!  inp_tot_m10/m15/m20/m25 fields each step (for output and microphysics coupling).
!  Convert to per-litre for obs comparison externally: n[#/L] = n[#/kg]*rho*1e-3.
!
! Processes:
!   - Emission       (inp_emit_bins)     : aerosol emission increment -> differential
!                                          INP bins, added to the tracers at the
!                                          surface. Called from the emission code
!                                          (GOCART dust, sea salt; PBAP number).
!   - Advection / nesting / boundaries   : handled by the host as for any chem tracer.
!   - Dry / wet deposition               : the INP tracers are removed with their
!                                          parent aerosol (wired in the deposition code).
!   - In-cloud activation (inp_prog_driver): bin k is depleted wherever a cell is in
!                                          cloud and its temperature has reached the
!                                          bin cold edge. Idempotent, needs no memory.
!
! ============================================================================
! !!! COEFFICIENT VERIFICATION REQUIRED BEFORE PRODUCTION RUNS !!!
!   VERIFIED here (published coefficients):
!     - Niemand et al. (2012)   desert-dust immersion n_s
!     - DeMott et al. (2015)    dust, number-based (D15)
!     - McCluskey et al. (2018) sea-spray-aerosol immersion n_s (M18)
!   PLACEHOLDERS -- functional form coded, COEFFICIENTS MUST BE SET FROM THE PAPERS
!   (search "TODO VERIFY"):
!     - Tobo et al. (2019)      Icelandic/high-latitude dust immersion n_s
!     - Hummel et al. (2018)    bacteria / fungal-spore frozen fractions
!   Default dust parameterization is Niemand-2012 (inp_dust_param=2).
! ============================================================================
!
! References:
!  (1) Da Silva, A., Marelle, L., Lapere, R., and Raut, J.-C. Modelling Ice
!      Nucleating Particles from High-Latitude Sources to Reproduce Arctic In
!      Situ Concentrations. In preparation, 2026
!  (2) Niemand et al. (2012)  - https://doi.org/10.1175/JAS-D-11-0249.1
!  (3) DeMott et al. (2015)   - https://doi.org/10.5194/acp-15-393-2015
!  (4) Tobo et al. (2019)     - https://doi.org/10.1038/s41561-019-0314-x
!  (5) McCluskey et al. (2018)- https://doi.org/10.1175/JAS-D-17-0155.1
!  (6) Ullrich et al. (2017)  - https://doi.org/10.1175/JAS-D-16-0074.1
!  (7) Hummel et al. (2018)   - https://doi.org/10.5194/acp-18-15437-2018

MODULE module_inp_emis

  USE module_model_constants, ONLY: PICONST

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: inp_prog_driver, inp_emit_dust, inp_emit_marine, inp_emit_bio

  !---- (C) PHYSICS TUNING KNOB: dust parameterization selector (inp_dust_param)
  INTEGER, PARAMETER :: INP_DUST_TOBO2019    = 1
  INTEGER, PARAMETER :: INP_DUST_NIEMAND2012 = 2
  INTEGER, PARAMETER :: INP_DUST_DEMOTT2015  = 3
  INTEGER, PARAMETER :: INP_DUST_ULLRICH2017 = 4

  !---- (B) PARTICLE DEFINITION: activation-temperature bins (deg C). INP_BIN_TC(k)
  !     is the COLD edge of differential bin k; bin k holds INP whose singular
  !     freezing temperature lies between INP_BIN_TC(k-1) (exclusive, warm) and
  !     INP_BIN_TC(k) (inclusive, cold), with INP_BIN_TC(0) taken as 0 C. Must
  !     match the bins declared in Registry/registry.chem (inp_*1..inp_*4).
  INTEGER, PARAMETER :: NUM_INP_BINS = 4
  REAL, PARAMETER :: INP_BIN_TC(NUM_INP_BINS) = (/ -10.0, -15.0, -20.0, -25.0 /)

  !---- (B) PARTICLE DEFINITION: GOCART size bins (dry log-centre radii, m) and
  !     densities (kg m-3), matching the NWFA/NIFA construction in the microphysics
  !     driver so the INP conversion assumes the same particle sizes as the model.
  REAL, PARAMETER :: rcen_dust(5) = (/ 0.3162E-6, 1.342E-6, 2.324E-6, 4.243E-6, 7.746E-6 /)
  REAL, PARAMETER :: rcen_seas(4) = (/ 0.2236E-6, 0.8660E-6, 2.739E-6, 7.071E-6 /)
  REAL, PARAMETER :: DENS_DUST = 2650.0   ! TODO VERIFY vs module_data_gocart_dust den_dust(:)
  REAL, PARAMETER :: DENS_SEAS = 2200.0   ! module_data_gocart_seas den_seas

  !---- (C) PHYSICS TUNING KNOBS: in-cloud activation e-folding timescale (s) and
  !     the cloud-condensate mixing ratio (kg kg-1) that defines "in cloud".
  REAL, PARAMETER :: TAU_ACT      = 300.0
  REAL, PARAMETER :: QC_CLOUD_MIN = 1.0E-6

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE inp_prog_driver( chem, t_phy, qcloud, qice, dt,             &
                  inp_tot_m10, inp_tot_m15, inp_tot_m20, inp_tot_m25,    &
                  ids,ide, jds,jde, kds,kde,                             &
                  ims,ime, jms,jme, kms,kme,                             &
                  its,ite, jts,jte, kts,kte                              )

    ! In-cloud activation sink on the differential INP bins, plus the cumulative
    ! TOTAL diagnostics. Called every chem step (after emission and transport).
    !
    ! Activation (singular hypothesis): a cell "in cloud" (qcloud+qice > QC_CLOUD_MIN)
    ! at temperature Tc removes the INP whose freezing temperature is warmer than or
    ! equal to Tc -- i.e. every bin whose cold edge INP_BIN_TC(b) is >= Tc. Removal
    ! is a first-order loss (e-folding TAU_ACT). This is idempotent: a bin already at
    ! zero stays at zero, and a bin whose cold edge has not been reached is untouched,
    ! so no coldest-temperature-experienced memory variable is needed.

    USE module_state_description, ONLY: num_chem, param_first_scalar,       &
         p_inp_dust1, p_inp_dust2, p_inp_dust3, p_inp_dust4,               &
         p_inp_mar1,  p_inp_mar2,  p_inp_mar3,  p_inp_mar4,                &
         p_inp_bio1,  p_inp_bio2,  p_inp_bio3,  p_inp_bio4

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ids,ide, jds,jde, kds,kde,                      &
                           ims,ime, jms,jme, kms,kme,                      &
                           its,ite, jts,jte, kts,kte
    REAL, INTENT(IN) :: dt
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), INTENT(INOUT) :: chem
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN)  :: t_phy, qcloud, qice
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT) ::           &
         inp_tot_m10, inp_tot_m15, inp_tot_m20, inp_tot_m25

    INTEGER :: i, j, k, b
    INTEGER :: pdu(NUM_INP_BINS), pma(NUM_INP_BINS), pbi(NUM_INP_BINS)
    REAL    :: tcC, f_act, cdu(NUM_INP_BINS), cma(NUM_INP_BINS), cbi(NUM_INP_BINS)

    ! INP tracers are only defined for the GOCART chem packages. If absent, their
    ! p-indices fall below param_first_scalar and there is nothing to do.
    IF (p_inp_dust1 < param_first_scalar) RETURN

    pdu = (/ p_inp_dust1, p_inp_dust2, p_inp_dust3, p_inp_dust4 /)
    pma = (/ p_inp_mar1,  p_inp_mar2,  p_inp_mar3,  p_inp_mar4  /)
    pbi = (/ p_inp_bio1,  p_inp_bio2,  p_inp_bio3,  p_inp_bio4  /)
    f_act = 1.0 - EXP(-dt/TAU_ACT)

    DO j = jts, jte
    DO k = kts, kte
    DO i = its, ite

      ! In-cloud activation: deplete the bins whose cold edge has been reached.
      IF (qcloud(i,k,j) + qice(i,k,j) > QC_CLOUD_MIN) THEN
        tcC = t_phy(i,k,j) - 273.15
        DO b = 1, NUM_INP_BINS
          IF (tcC <= INP_BIN_TC(b)) THEN
            chem(i,k,j,pdu(b)) = chem(i,k,j,pdu(b)) * (1.0 - f_act)
            chem(i,k,j,pma(b)) = chem(i,k,j,pma(b)) * (1.0 - f_act)
            chem(i,k,j,pbi(b)) = chem(i,k,j,pbi(b)) * (1.0 - f_act)
          END IF
        END DO
      END IF

      ! Cumulative TOTAL n_INP(T) = partial sums over bins, all species.
      ! Clamp >= 0 to absorb the tiny negatives advection leaves behind.
      DO b = 1, NUM_INP_BINS
        cdu(b) = MAX(chem(i,k,j,pdu(b)), 0.0)
        cma(b) = MAX(chem(i,k,j,pma(b)), 0.0)
        cbi(b) = MAX(chem(i,k,j,pbi(b)), 0.0)
      END DO
      inp_tot_m10(i,k,j) = cdu(1)+cma(1)+cbi(1)
      inp_tot_m15(i,k,j) = inp_tot_m10(i,k,j) + cdu(2)+cma(2)+cbi(2)
      inp_tot_m20(i,k,j) = inp_tot_m15(i,k,j) + cdu(3)+cma(3)+cbi(3)
      inp_tot_m25(i,k,j) = inp_tot_m20(i,k,j) + cdu(4)+cma(4)+cbi(4)

    END DO
    END DO
    END DO

  END SUBROUTINE inp_prog_driver

!------------------------------------------------------------------------------
  SUBROUTINE inp_emit_dust( edust_ugkg, rho, inp_dust_param, du_diff )
    ! Emit-once: this step's emitted dust mass increment (ug kg-1, per GOCART bin)
    ! -> differential dust INP bins (# kg-1) to ADD to the dust INP tracers.
    !   cumulative n_INP(T) = dust surface area * n_s(T)  [Niemand/Tobo/Ullrich]
    ! NOTE: the DeMott-2015 option is number-based and NON-LINEAR in dust number,
    ! so it is not strictly additive over emission events; it is applied here to
    ! the per-step increment as an approximation. The n_s options (default
    ! Niemand-2012) are linear in surface area and the natural choice for emit-once.
    IMPLICIT NONE
    REAL,    INTENT(IN)  :: edust_ugkg(5), rho
    INTEGER, INTENT(IN)  :: inp_dust_param
    REAL,    INTENT(OUT) :: du_diff(NUM_INP_BINS)
    INTEGER :: n, b
    REAL :: num, sa, sadust_tot, ndust_tot, n05_cm3, tc, cdu(NUM_INP_BINS)

    sadust_tot = 0.0 ; ndust_tot = 0.0
    DO n = 1, 5
      CALL gocart_mass_to_num_sa(edust_ugkg(n), rcen_dust(n), DENS_DUST, num, sa)
      sadust_tot = sadust_tot + sa ; ndust_tot = ndust_tot + num
    END DO
    n05_cm3 = ndust_tot * rho * 1.0E-6   ! emitted dust number > 0.5 um diam, cm-3 (D15)

    DO b = 1, NUM_INP_BINS
      tc = INP_BIN_TC(b)
      IF (inp_dust_param == INP_DUST_DEMOTT2015) THEN
        cdu(b) = ninp_demott2015_dust(n05_cm3, tc+273.15) / MAX(rho*1.0E-3, 1.0E-12)
      ELSE
        cdu(b) = sadust_tot * ns_dust(tc, inp_dust_param)
      END IF
    END DO
    CALL cum_to_diff(cdu, du_diff)
  END SUBROUTINE inp_emit_dust

  SUBROUTINE inp_emit_marine( eseas_ugkg, mar_diff )
    ! Emit-once: this step's emitted sea-spray mass increment (ug kg-1, per GOCART
    ! bin) -> differential marine INP bins (# kg-1). Cumulative n_INP(T) = sea-spray
    ! surface area * n_s(T) (McCluskey-2018 M18).
    IMPLICIT NONE
    REAL, INTENT(IN)  :: eseas_ugkg(4)
    REAL, INTENT(OUT) :: mar_diff(NUM_INP_BINS)
    INTEGER :: n, b
    REAL :: num, sa, saseas_tot, cma(NUM_INP_BINS)

    saseas_tot = 0.0
    DO n = 1, 4
      CALL gocart_mass_to_num_sa(eseas_ugkg(n), rcen_seas(n), DENS_SEAS, num, sa)
      saseas_tot = saseas_tot + sa
    END DO
    DO b = 1, NUM_INP_BINS
      cma(b) = saseas_tot * ns_marine(INP_BIN_TC(b))
    END DO
    CALL cum_to_diff(cma, mar_diff)
  END SUBROUTINE inp_emit_marine

  SUBROUTINE inp_emit_bio( ebact_num, espore_num, bio_diff )
    ! Emit-once: this step's emitted PBAP number increment (# kg-1) -> differential
    ! biological INP bins (# kg-1). Cumulative n_INP(T) = N_bact*ff_bact(T)
    ! + N_spore*ff_spore(T) (per-particle frozen fractions).
    IMPLICIT NONE
    REAL, INTENT(IN)  :: ebact_num, espore_num
    REAL, INTENT(OUT) :: bio_diff(NUM_INP_BINS)
    INTEGER :: b
    REAL :: cbi(NUM_INP_BINS)

    DO b = 1, NUM_INP_BINS
      cbi(b) = ebact_num*ff_bact(INP_BIN_TC(b)) + espore_num*ff_spore(INP_BIN_TC(b))
    END DO
    CALL cum_to_diff(cbi, bio_diff)
  END SUBROUTINE inp_emit_bio

  SUBROUTINE cum_to_diff( cum, diff )
    ! Cumulative n_INP at the bin cold edges -> non-negative differential bins.
    IMPLICIT NONE
    REAL, INTENT(IN)  :: cum(NUM_INP_BINS)
    REAL, INTENT(OUT) :: diff(NUM_INP_BINS)
    INTEGER :: b
    diff(1) = MAX(cum(1), 0.0)
    DO b = 2, NUM_INP_BINS
      diff(b) = MAX(cum(b) - cum(b-1), 0.0)
    END DO
  END SUBROUTINE cum_to_diff

!------------------------------------------------------------------------------
  SUBROUTINE gocart_mass_to_num_sa(mass_ugkg, rcen, dens, num, sa)
    ! Convert a bulk/bin aerosol mass mixing ratio (ug kg-1) to number (# kg-1) and
    ! surface area (m2 kg-1), assuming all mass sits at the bin centre radius rcen (m)
    ! with density dens (kg m-3). Mirrors the GOCART NWFA/NIFA construction.
    IMPLICIT NONE
    REAL, INTENT(IN)  :: mass_ugkg, rcen, dens
    REAL, INTENT(OUT) :: num, sa
    REAL :: mass_per_particle_ug

    ! mass per particle (ug) = 1e9 * (4/3) pi r^3 rho
    mass_per_particle_ug = 1.0E9 * (4.0/3.0)*PICONST * dens * rcen**3
    num = MAX(mass_ugkg, 0.0) / MAX(mass_per_particle_ug, 1.0E-30)
    sa  = num * 4.0*PICONST * rcen**2
  END SUBROUTINE gocart_mass_to_num_sa

!==============================================================================
! Layer B physics: (C) PHYSICS TUNING KNOBS -- ice-active surface-site densities
! n_s(T) [m-2] and biological frozen fractions. Tc is the temperature in deg C.
!==============================================================================

  REAL FUNCTION ns_dust(tc, param)
    IMPLICIT NONE
    REAL,    INTENT(IN) :: tc
    INTEGER, INTENT(IN) :: param
    REAL :: tcl

    SELECT CASE (param)

    CASE (INP_DUST_NIEMAND2012)
      ! VERIFIED: Niemand et al. (2012), desert-dust immersion.
      ! n_s = exp(-0.517*Tc + 8.934) m-2, valid -36 <= Tc <= -12 C.
      IF (tc > -12.0) THEN
        ns_dust = 0.0
      ELSE
        tcl = MAX(tc, -36.0)
        ns_dust = EXP(-0.517*tcl + 8.934)
      END IF

    CASE (INP_DUST_TOBO2019)
      ! TODO VERIFY: Tobo et al. (2019) Icelandic/high-latitude dust immersion.
      ! Functional form exp(a*Tc + b); SET a,b FROM THE PAPER. Placeholder below
      ! reuses the Niemand slope with a +1 decade offset (NOT the paper values).
      IF (tc > -7.0) THEN
        ns_dust = 0.0
      ELSE
        tcl = MAX(tc, -35.0)
        ns_dust = EXP(-0.517*tcl + 8.934 + 2.3026)   !!! TODO VERIFY coefficients
      END IF

    CASE (INP_DUST_ULLRICH2017)
      ! TODO VERIFY: Ullrich et al. (2017) dust immersion n_s polynomial.
      ns_dust = 0.0                                   !!! TODO VERIFY: implement

    CASE DEFAULT
      ns_dust = 0.0
    END SELECT
  END FUNCTION ns_dust

  REAL FUNCTION ns_marine(tc)
    ! VERIFIED: McCluskey et al. (2018) marine sea-spray-aerosol immersion n_s
    ! (the "M18" parameterization). The published marine INP fit is
    !   n_INP[m-3] = exp(-0.545*Tc + 1.0125) * SSA[m2 m-3],
    ! i.e. an ice-active surface-site density
    !   n_s(Tc) = exp(-0.545*Tc + 1.0125)  [m-2],  Tc in deg C.
    ! Derived from pristine North Atlantic air; evaluated against Southern Ocean
    ! INP down to ~ -25 C. Represents the regularly-occurring marine INP class
    ! only (heat-labile organic INP not included). Applied over its marine
    ! mixed-phase range: warm bound -5 C (warmer INP is the biological tail),
    ! cold extrapolation clamped at -30 C.
    IMPLICIT NONE
    REAL, INTENT(IN) :: tc
    REAL :: tcl
    IF (tc > -5.0) THEN
      ns_marine = 0.0
    ELSE
      tcl = MAX(tc, -30.0)
      ns_marine = EXP(-0.545*tcl + 1.0125)            ! McCluskey et al. (2018), M18
    END IF
  END FUNCTION ns_marine

  REAL FUNCTION ff_bact(tc)
    ! TODO VERIFY: Hummel et al. (2018) bacterial frozen fraction f(T) in [0,1].
    ! Bacteria carry the warm tail (active up to ~ -5 C for ice-active strains).
    IMPLICIT NONE
    REAL, INTENT(IN) :: tc
    ff_bact = 0.0                                     !!! TODO VERIFY: implement
    ! Placeholder smooth ramp (NOT paper values), so the field is non-trivial:
    IF (tc <= -2.0) ff_bact = MIN(MAX((-2.0 - tc)/18.0, 0.0), 1.0)
  END FUNCTION ff_bact

  REAL FUNCTION ff_spore(tc)
    ! TODO VERIFY: Hummel et al. (2018) fungal-spore frozen fraction f(T) in [0,1].
    ! Spores activate colder than bacteria (~ -15 to -25 C).
    IMPLICIT NONE
    REAL, INTENT(IN) :: tc
    ff_spore = 0.0                                    !!! TODO VERIFY: implement
    IF (tc <= -12.0) ff_spore = MIN(MAX((-12.0 - tc)/15.0, 0.0), 1.0)
  END FUNCTION ff_spore

  REAL FUNCTION ninp_demott2015_dust(n05_cm3, tk)
    ! VERIFIED form: DeMott et al. (2015) D15, dust.
    !   n_INP[std L-1] = cf * (n05)^(a*(273.16-Tk)+b) * exp(c*(273.16-Tk)+d)
    !   dust: cf=3, a=0, b=1.25, c=0.46, d=-11.6 ; n05 in cm-3
    IMPLICIT NONE
    REAL, INTENT(IN) :: n05_cm3, tk
    REAL, PARAMETER :: cf=3.0, a=0.0, b=1.25, c=0.46, d=-11.6
    REAL :: dt
    IF (n05_cm3 <= 0.0 .OR. tk >= 273.16) THEN
      ninp_demott2015_dust = 0.0
      RETURN
    END IF
    dt = 273.16 - tk
    ninp_demott2015_dust = cf * (n05_cm3**(a*dt + b)) * EXP(c*dt + d)
  END FUNCTION ninp_demott2015_dust

END MODULE module_inp_emis

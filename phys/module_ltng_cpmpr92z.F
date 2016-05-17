! WRF:MODEL_LAYER:PHYSICS
!
! Lightning flash rate prediction based on cloud-top height. Implemented
! for models using convective parameterization. Assume use of sub-grid LNB.
!
! Price, C., and D. Rind (1992), A Simple Lightning Parameterization for Calculating
!   Global Lightning Distributions, J. Geophys. Res., 97(D9), 9919-9933, doi:10.1029/92JD00719.
!
! Wong, J., M. Barth, and D. Noone (2012), Evaluating a Lightning Parameterization
!   at Resolutions with Partially-Resolved Convection, GMDD, in preparation.
!
! Contact: J. Wong <johnwong@ucar.edu>
!
!**********************************************************************
!
!lm Louis Marelle 2016/5/17: Add lightning parameterization for KF-CuP cumuli

 MODULE module_ltng_cpmpr92z
 CONTAINS

 SUBROUTINE ltng_cpmpr92z ( &
                          ! Frequently used prognostics
                            dx, dy, xland, ht, z, t,              &
                          ! Scheme specific prognostics
                            kLNB,                                 &
                          ! Scheme specific namelist inputs
                            cldtop_adjustment,                    &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                          ! Mandatory output for all quantitative schemes
                            total_flashrate                       &
                          )
!-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z, t

! Scheme specific prognostics
 INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: kLNB     ! model LNB from cu_physics

! Scheme specific namelist inputs
 REAL,    INTENT(IN   )    ::       cldtop_adjustment

! Order dependent args for domain, mem, and tile dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: total_flashrate

! Local variables
 REAL :: dA              ! grid area dx*dy in km2
 REAL :: zkm             ! AGL z in km

 REAL, PARAMETER:: baseArea=1296. ! base-case area, dx = 36 km

 INTEGER :: i,k,j

!-----------------------------------------------------------------

 dA = dx*dy/1E6

 total_flashrate( ips:ipe,jps:jpe ) = 0.

! Compute AGL heights in km
 jloop: DO j=jps,jpe
    iloop: DO i=ips,ipe
        IF ( t(i,kLNB(i,j),j) .lt. 273.15 .and. &
            kLNB(i,j) .ge. kps .and. kLNB(i,j) .le. kpe ) THEN              ! Cloud top must freeze
            zkm = ( z(i,kLNB(i,j),j) - ht(i,j) )/1E3 + cldtop_adjustment    ! Compute AGL height in km
            IF ( zkm .gt. 0. ) THEN                                         ! Must be above ground
              IF ( xland(i,j) .lt. 1.5 ) THEN
                total_flashrate(i,j) = 3.44E-5 * (zkm**4.9) /60.            ! Continental equation
              ELSE
                total_flashrate(i,j) = 6.57E-6 * (zkm**4.9) /60.            ! Marine equation (Michalon 99)
              ENDIF
            ENDIF
        ENDIF
    ENDDO iloop
 ENDDO jloop

! Scale by areal ratio
 total_flashrate(ips:ipe,jps:jpe) = total_flashrate(ips:ipe,jps:jpe) * dA/baseArea

 END SUBROUTINE ltng_cpmpr92z
 
!lm Add a subroutine for KFCUP
 SUBROUTINE ltng_cpmpr92z_cup ( &
                          ! Frequently used prognostics
                            dx, dy, xland, ht, z, t,              &
                          ! Scheme specific namelist inputs
                            cldtop_adjustment,                    &
                          ! Additional inputs for lightning param with KF-CuP cumuli
                            msft, &
                            cldfra_cup, &
                            htop, &
                            shall, &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                          ! Mandatory output for all quantitative schemes
                            total_flashrate                       &
                          )
!-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z, t

! Scheme specific namelist inputs
 REAL,    INTENT(IN   )    ::       cldtop_adjustment

! Additional inputs for lightning param with KF-CuP cumuli
REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ),           INTENT(IN   ) :: cldfra_cup ! KF-CuP cloud fraction 
REAL,    DIMENSION( ims:ime,          jms:jme ),           INTENT(IN   ) :: htop ! cloud top index
REAL,    DIMENSION( ims:ime,          jms:jme ),           INTENT(IN   ) :: shall ! flag indicating shallow/deep/no convection
REAL,    DIMENSION( ims:ime,          jms:jme ),           INTENT(IN   ) :: msft ! map factors for grid area scaling

! Order dependent args for domain, mem, and tile dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: total_flashrate

! Local variables
 REAL :: dA              ! grid area dx*dy in km2
 REAL :: zkm             ! AGL z in km
! Scaling factors for lightning rates (param is grid dependent)
 REAL :: cldfra_scaling_f    !average cloud fraction between cloud bottom and cloud top
 REAL :: mapfac_scaling_f    !squared map factor, should be mapfac_x * mapfax_y
                             !but the param is a rough approximation anyway 
 INTEGER :: kbtm, ktop  ! cloud bottom index, top index

 REAL, PARAMETER:: baseArea=1296. ! base-case area, dx = 36 km

 INTEGER :: i,k,j

!-----------------------------------------------------------------

!Initialize scaling factors 
 cldfra_scaling_f=1. !cloud fraction = 1
 mapfac_scaling_f=1. !no grid distortion
!Initialize other local variables
 ktop = -1
 kbtm = -1
 zkm = 0.

 total_flashrate( ips:ipe,jps:jpe ) = 0.

!Grid area in km, assuming no projection distortion for now
 dA = dx*dy/1E6

! Compute flash rates from AGL cloud top height 
 jloop: DO j=jps,jpe
    iloop: DO i=ips,ipe
        IF (shall(i,j) .le. 0.5) THEN ! deep convection
          ! Cloud top level ktop directly from KF-CuP
          ktop=nint(htop(i,j))
          ! Look for cloud bottom level kbtm (used here to compute the average cloud fraction)
          k = kps
          DO WHILE( cldfra_cup(i,k,j) .le. 0.01 .and. k .le. ktop )
            k = k+1
          ENDDO
          kbtm = k
          IF ( (ktop - kbtm) > 0 ) THEN !cloud depth is > 0
            cldfra_scaling_f=SUM(cldfra_cup(i,kbtm:ktop,j))/REAL(ktop-kbtm+1) ! "cloud-wide" average cloud fraction
            IF ( t(i,ktop,j) .lt. 273.15 .and. &
              ktop .ge. kps .and. ktop .le. kpe ) THEN              ! Cloud top must freeze for lightning
              zkm = ( z(i,ktop,j) - ht(i,j) )/1E3 + cldtop_adjustment    ! Compute AGL height in km
              IF ( zkm .gt. 0. ) THEN                                         ! Must be above ground
                IF ( xland(i,j) .lt. 1.5 ) THEN
                  total_flashrate(i,j) = 3.44E-5 * (zkm**4.9) /60.            ! Continental equation
                ELSE
                  total_flashrate(i,j) = 6.57E-6 * (zkm**4.9) /60.            ! Marine equation (Michalon 99)
                ENDIF
! Scale flashrate by cloud fraction (the parameterization above assumes that the cloud fills the whole
! grid)
              total_flashrate(i,j) = total_flashrate(i,j) * cldfra_scaling_f 
! Scale flashrate by areal ratio (the parameterization above is cell-area dependent). Take into account map projection distortions
              mapfac_scaling_f = 1.0 / (msft(i,j) * msft(i,j))
              total_flashrate(i,j) = total_flashrate(i,j) * dA/baseArea * mapfac_scaling_f  
              ENDIF ! zkm .gt. 0.
            ENDIF ! t(i,ktop,j) .lt. 273.15
          ENDIF ! (ktop - kbtm) > 0
        ENDIF ! shall(i,j) .le. 0.5
    ENDDO iloop
 ENDDO jloop

 END SUBROUTINE ltng_cpmpr92z_cup
!lm

 END MODULE module_ltng_cpmpr92z

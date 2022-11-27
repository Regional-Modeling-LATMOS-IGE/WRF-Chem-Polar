! 2018/10/04 Louis Marelle
! Compute DMS emissions and add them to chem(:,:,:,p_dms)
! I developed this for CBMZ_MOSAIC with DMS and SAPRC_MOSAIC with DMS reactions
! added to KPP, but this should work for any mechanism where DMS is defined in
! chem. However, if someone tries to use this option with another mechanism I
! included a call to WRF_FATAL in chemics_init.F.
! This routines requires that e_dms_oc, the oceanic dms content, is
! present in emis_ant (read from the wrfchemi files). This is checked in
! chem_init.F
!
! This routine needs dms oceanic content as an input, read in wrfchemi files
! (emis_ant) for now. This oceanic content can be overwritten by another
! value here if needed (dms_ocean_sfc), see below.
!
! This module is based on module_cbmz_addemiss.F
!
! This code was developped for the following publication, also showing tests in
! the Arctic:
! Marelle, L., J.-C. Raut, K. S. Law, L. K. Berg, J. D. Fast, R. C. Easter, M.
! Shrivastava, and J. L. Thomas (2017), Improvements to the WRF-Chem 3.5.1 model
! for quasi- hemispheric simulations of aerosols and ozone in the Arctic,
! Geoscientific Model Devel- opment, 10(10), 3661–3677,
! doi:10.5194/gmd-10-3661-2017.

MODULE module_nightingale_dmsemis

CONTAINS

!----------------------------------------------------------------------
   subroutine nightingale_dmsemis( id, dtstep, dz8w, config_flags,              &
               u_phy, v_phy, rho_phy, chem, emis_ant, alt,                &
               u10, v10, xland, tsk,                                      &
               ids,ide, jds,jde, kds,kde,                                 &
               ims,ime, jms,jme, kms,kme,                                 &
               its,ite, jts,jte, kts,kte                                  )
!
! calculates and adds emissions to dms in chem
! (i.e., emissions tendencies over time dtstep are applied
! to the trace gas concentrations)
!

  USE module_configure
  USE module_state_description
  USE module_data_radm2

  IMPLICIT NONE

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: id,                                      &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte

   REAL, INTENT(IN   ) ::    dtstep

! trace species mixing ratios (gases=ppm)
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::   chem
!
! emissions arrays (v.1: ppm m/min; v.2: mole km^-2 hr^-1)
!
!   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                          &
   REAL, DIMENSION( ims:ime, kms:config_flags%kemit, jms:jme,num_emis_ant),&
         INTENT(IN ) ::                                                    &
                         emis_ant
! layer thickness (m)
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   dz8w, rho_phy,alt, u_phy, v_phy

   REAL,  DIMENSION( ims:ime , jms:jme )         ,               &
          INTENT(IN   ) ::   u10, v10, xland, tsk

! local variables
   integer :: i,j,k
   real :: conv
! dms_ocean_sfc surface oceanic dms content in mol/m^3, e.g.  Lana et. al 2011
! dms_ocean_sfc is overwritten by the value in the wrfchemi emission files
! sc_co2 schmidt number for co2 at 20°C = 600 (reference value)
   real :: w10, dms_emi, k_dms, dms_ocean_sfc = 0.0, sc_co2 = 600.0, sc_dms, &
                sst_ij_cels

!
! add emissions
!

   do j=jts,jte
      do i=its,ite

         conv = 4.828e-4 / rho_phy(i,kts,j) * dtstep / (dz8w(i,kts,j) * 60.)

         ! Nightingale et al., 2000 parameterization using a Schmidt number
         ! for DMS depending on the actual SST
             ! only over open water (xland>1.5)
             if (xland(i,j) .gt. 1.5) then

                 ! Wind speed in m s-1
                 w10 = sqrt( u10(i,j) ** 2.0 + v10(i,j) ** 2.0 )

                 ! we donṫ trust the u10,v10 values, is model layers are very
                 ! thin near surface (this recommandation is from the gocart modules)
                 if ( dz8w(i,kts,j) .lt. 12.) &
                   w10 = sqrt( u_phy(i,kts,j) * u_phy(i,kts,j) + v_phy(i,kts,j) * v_phy(i,kts,j) )

                 ! Convert sst (=tskin at sea surface) to celsius. The Schmidt
                 ! number parameterization is only valid between 5 and 30 °C,
                 ! so constrain the SST between these two values
                 sst_ij_cels = ( tsk(i,j) - 273.15 )
                 if (sst_ij_cels < 5.0) sst_ij_cels = 5.0
                 if (sst_ij_cels > 30.0) sst_ij_cels = 30.0

                 ! Schmidt number for DMS, using the sst in celsius, from
                 ! Saltzman et al., 1993
                 sc_dms = 2674.0 - 147.12 * sst_ij_cels + 3.726 * sst_ij_cels ** 2.0 &
                             - 0.038 * sst_ij_cels ** 3

                 ! Sea - air exchange coefficient from Nightingale et al. 2000,
                 ! and unit change from cm hr-1 to m s-1
                 k_dms = ( 0.222 * w10**2 + 0.333 * w10) * (sc_dms / sc_co2)**(-0.5) / 3600. / 100.

                 ! This can be commented in order to use a single
                 ! (0D) hardcoded value of dms_ocean_sfc if needed
                 dms_ocean_sfc=emis_ant(i,1,j,p_e_dms_oc)

                 ! DMS emission flux in SI units (mol m-2 s-1)
                 dms_emi = k_dms * dms_ocean_sfc
                 ! Change to  WRF-Chem units (mol km-2 hr-1)
                 dms_emi = dms_emi * 1.0E6 * 3600.0

                 ! Add emissions to chem, at the surface (k==kts)
                 chem(i,kts,j,p_dms)  =  chem(i,kts,j,p_dms)               &
                             +dms_emi*conv

             end if !xland .gt 1.5, over open water
      end do
   end do

END subroutine nightingale_dmsemis

END MODULE module_nightingale_dmsemis

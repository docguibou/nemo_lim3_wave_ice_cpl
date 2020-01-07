MODULE sbc_ice
!!======================================================================
!!                 ***  MODULE  sbc_ice  ***
!! Surface module - LIM-3: parameters & variables defined in memory
!!======================================================================
!! History :  3.0  ! 2006-08  (G. Madec)  Surface module
!!            3.2  ! 2009-06  (S. Masson) merge with ice_oce
!!            3.3.1! 2011-01  (A. R. Porter, STFC Daresbury) dynamical allocation
!!            3.4  ! 2011-11  (C. Harris) CICE added as an option
!!----------------------------------------------------------------------

!!----------------------------------------------------------------------
!!   'key_lim2' or 'key_lim3' :             LIM-2 or LIM-3 sea-ice model
!!----------------------------------------------------------------------
   USE par_oce          ! ocean parameters
   USE sbc_oce          ! surface boundary condition: ocean

   USE ice              ! LIM-3 parameters



   USE lib_mpp          ! MPP library
   USE in_out_manager   ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC sbc_ice_alloc ! called in iceini(_2).F90



   LOGICAL         , PUBLIC, PARAMETER ::   lk_lim2    = .FALSE.  !: no LIM-2
   LOGICAL         , PUBLIC, PARAMETER ::   lk_lim3    = .TRUE.   !: LIM-3 ice model
   LOGICAL         , PUBLIC, PARAMETER ::   lk_cice    = .FALSE.  !: no CICE
   LOGICAL         , PUBLIC            ::   ln_wavice  = .FALSE.  !WAVE ICE MODELLING GB
   LOGICAL         , PUBLIC            ::   cpl_taiwx  = .FALSE.  !WAVE ICE MODELLING GB
   LOGICAL         , PUBLIC            ::   cpl_taiwy  = .FALSE.  !WAVE ICE MODELLING GB
   LOGICAL         , PUBLIC            ::   cpl_icef   = .FALSE.  !WAVE ICE MODELLING GB
   CHARACTER(len=1), PUBLIC, PARAMETER ::   cp_ice_msh = 'C'      !: 'C'-grid ice-velocity



   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qns_ice        !: non solar heat flux over ice                  [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qsr_ice        !: solar heat flux over ice                      [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qla_ice        !: latent flux over ice                          [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   dqla_ice       !: latent sensibility over ice                 [W/m2/K]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   dqns_ice       !: non solar heat flux over ice (LW+SEN+LA)    [W/m2/K]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tn_ice         !: ice surface temperature                          [K]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   alb_ice        !: ice albedo                                       [-]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   utau_ice       !: atmos-ice u-stress. VP: I-pt ; EVP: U,V-pts   [N/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   vtau_ice       !: atmos-ice v-stress. VP: I-pt ; EVP: U,V-pts   [N/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   fr1_i0         !: Solar surface transmission parameter, thick ice  [-]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   fr2_i0         !: Solar surface transmission parameter, thin ice   [-]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   emp_ice        !: sublimation - precip over sea ice          [kg/m2/s]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   topmelt            !: category topmelt
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   botmelt            !: category botmelt


   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   evap_ice       !: sublimation                              [kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   devap_ice      !: sublimation sensitivity                [kg/m2/s/K]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qns_oce        !: non solar heat flux over ocean              [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qsr_oce        !: non solar heat flux over ocean              [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qemp_oce       !: heat flux of precip and evap over ocean     [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qemp_ice       !: heat flux of precip and evap over ice       [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qevap_ice      !: heat flux of evap over ice                  [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qprec_ice      !: enthalpy of precip over ice                 [J/m3]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   emp_oce        !: evap - precip over ocean                 [kg/m2/s]
!GB ice wave coupling
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   TAIWX        !: radiative stress from waves to ice       [m2/s2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   TAIWY        !: radiative stress from waves to ice       [m2/s2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ICEF_r         !: Ice floe diameter from wave model        [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ICEF_s         !: Ice floe diameter sent o wave model        [m]


   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wndm_ice       !: wind speed module at T-point                 [m/s]



   




   REAL(wp), PUBLIC, SAVE ::   cldf_ice = 0.81    !: cloud fraction over sea ice, summer CLIO value   [-]

!!----------------------------------------------------------------------
!! NEMO/OPA 4.0 , NEMO Consortium (2011)
!! $Id: sbc_ice.F90 6416 2016-04-01 12:22:17Z clem $
!! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sbc_ice_alloc()
!!----------------------------------------------------------------------
!!                     ***  FUNCTION sbc_ice_alloc  ***
!!----------------------------------------------------------------------
      INTEGER :: ierr(5)
!!----------------------------------------------------------------------
      ierr(:) = 0


      ALLOCATE( qns_ice (jpi,jpj,jpl) , qsr_ice (jpi,jpj,jpl) ,     &
         &      qla_ice (jpi,jpj,jpl) , dqla_ice(jpi,jpj,jpl) ,     &
         &      dqns_ice(jpi,jpj,jpl) , tn_ice  (jpi,jpj,jpl) , alb_ice (jpi,jpj,jpl) ,   &
         &      utau_ice(jpi,jpj)     , vtau_ice(jpi,jpj)     , wndm_ice(jpi,jpj)     ,   &
         &      fr1_i0  (jpi,jpj)     , fr2_i0  (jpi,jpj)     , TAIWX(jpi,jpj)        ,   &
         &      TAIWY(jpi,jpj)        , ICEF_s(jpi,jpj)       ,ICEF_r(jpi,jpj)        ,   &


         &      evap_ice(jpi,jpj,jpl) , devap_ice(jpi,jpj,jpl) , qprec_ice(jpi,jpj) ,   &
         &      qemp_ice(jpi,jpj)     , qevap_ice(jpi,jpj,jpl) , qemp_oce (jpi,jpj) ,   &
         &      qns_oce (jpi,jpj)     , qsr_oce  (jpi,jpj)     , emp_oce (jpi,jpj)  ,   &

         &      emp_ice(jpi,jpj)      ,  STAT= ierr(1) )



!


      sbc_ice_alloc = MAXVAL( ierr )
      IF( lk_mpp            )   CALL mpp_sum ( sbc_ice_alloc )
      IF( sbc_ice_alloc > 0 )   CALL ctl_warn('sbc_ice_alloc: allocation of arrays failed')
   END FUNCTION sbc_ice_alloc



!!======================================================================
END MODULE sbc_ice

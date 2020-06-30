MODULE usrdef_sbc
   !!======================================================================
   !!                       ***  MODULE  usrdef_sbc  ***
   !! 
   !!                     ===  SAS_BIPER configuration  ===
   !!
   !! User defined :   surface forcing of a user configuration
   !!======================================================================
   !! History :  4.0   ! 2016-03  (S. Flavoni, G. Madec)  user defined interface
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_sbc    : user defined surface bounday conditions in SAS_BIPER case
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE sbc_ice         ! Surface boundary condition: ice fields
   USE phycst          ! physical constants
   USE ice, ONLY       : pfrld, a_i_b
   USE limthd_dh       ! for CALL lim_thd_snwblow
   !
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined) 
   USE wrk_nemo

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usrdef_sbc_oce      ! routine called by sbcmod.F90 for sbc ocean
   PUBLIC   usrdef_sbc_ice_tau  ! routine called by sbcice_lim.F90 for ice dynamics
   PUBLIC   usrdef_sbc_ice_flx  ! routine called by sbcice_lim.F90 for ice thermo

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2016)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usrdef_sbc_oce( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE usr_def_sbc  ***
      !!              
      !! ** Purpose :   provide at each time-step the surface boundary
      !!              condition, i.e. the momentum, heat and freshwater fluxes.
      !!
      !! ** Method  :   all 0 fields, for SAS_BIPER case
      !!                CAUTION : never mask the surface stress field !
      !!
      !! ** Action  : - set to ZERO all the ocean surface boundary condition, i.e.   
      !!                   utau, vtau, taum, wndm, qns, qsr, emp, sfx
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !!---------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         !
         IF(lwp)   WRITE(numout,*)' usrdef_sbc_oce : SAS_BIPER case: NO surface forcing'
         ! --- oce variables --- !
         utau(:,:) = 0._wp
         vtau(:,:) = 0._wp
         taum(:,:) = 0._wp
         wndm(:,:) = 0._wp
         !
         emp (:,:) = 0._wp
         sfx (:,:) = 0._wp
         qns (:,:) = 0._wp
         qsr (:,:) = 0._wp
         !
      ENDIF
      !
   END SUBROUTINE usrdef_sbc_oce

   SUBROUTINE usrdef_sbc_ice_tau( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE usrdef_sbc_ice_tau  ***
      !!
      !! ** Purpose :   provide the surface boundary (momentum) condition over sea-ice
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !!---------------------------------------------------------------------
      IF( kt==nit000 .AND. lwp)   WRITE(numout,*)' usrdef_sbc_ice : SAS_BIPER case: NO stress forcing'
      !
      utau_ice(:,:) = 0._wp
      vtau_ice(:,:) = 0._wp
      !
   END SUBROUTINE usrdef_sbc_ice_tau

   SUBROUTINE usrdef_sbc_ice_flx( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE usrdef_sbc_ice_flx  ***
      !!
      !! ** Purpose :   provide the surface boundary (flux) condition over sea-ice
      !!---------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), POINTER ::   zsnw       ! snw distribution after wind blowing
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !!---------------------------------------------------------------------
      CALL wrk_alloc( jpi,jpj, zsnw )
      !
      IF( kt==nit000 .AND. lwp)   WRITE(numout,*)' usrdef_sbc_ice : SAS_BIPER case: NO flux forcing'
      !
      ! ocean variables (renaming)
      emp_oce (:,:)   = 0._wp   ! uniform value for freshwater budget (E-P)
      qsr_oce (:,:)   = 0._wp   ! uniform value for     solar radiation
      qns_oce (:,:)   = 0._wp   ! uniform value for non-solar radiation
      !sst_m(:,:)      = 1._wp
      ! ice variables
      alb_ice (:,:,:) = 0.7_wp  ! useless
      qsr_ice (:,:,:) = 1._wp   ! uniform value for     solar radiation
      qns_ice (:,:,:) = 1.0_wp   ! uniform value for non-solar radiation
      sprecip (:,:)   = 0._wp   ! uniform value for snow precip
      evap_ice(:,:,:) = 0._wp   ! uniform value for sublimation

      ! ice fields deduced from above
      zsnw(:,:) = 1._wp
      !!CALL lim_thd_snwblow( pfrld, zsnw )  ! snow distribution over ice after wind blowing 
      emp_ice  (:,:)   = SUM( a_i_b(:,:,:) * evap_ice(:,:,:), dim=3 ) - sprecip(:,:) * zsnw(:,:)
      emp_oce  (:,:)   = emp_oce(:,:) - sprecip(:,:) * (1._wp - zsnw(:,:) )
      qevap_ice(:,:,:) =   0._wp
      qprec_ice(:,:)   =   rhosn * ( sst_m(:,:) * cpic - lfus ) * tmask(:,:,1) !  in J/m3
      qemp_oce (:,:)   = - emp_oce(:,:) * sst_m(:,:) * rcp
      qemp_ice (:,:)   =   sprecip(:,:) * zsnw * ( sst_m(:,:) * cpic - lfus ) * tmask(:,:,1) ! solid precip (only)

      ! total fluxes
      emp_tot (:,:) = emp_ice  + emp_oce
      qns_tot (:,:) = pfrld(:,:) * qns_oce(:,:) + SUM( a_i_b(:,:,:) * qns_ice(:,:,:), dim=3 ) + qemp_ice(:,:) + qemp_oce(:,:)
      qsr_tot (:,:) = pfrld(:,:) * qsr_oce(:,:) + SUM( a_i_b(:,:,:) * qsr_ice(:,:,:), dim=3 )

      !--------------------------------------------------------------------
      ! FRACTIONs of net shortwave radiation which is not absorbed in the
      ! thin surface layer and penetrates inside the ice cover
      ! ( Maykut and Untersteiner, 1971 ; Ebert and Curry, 1993 )
      fr1_i0(:,:) = ( 0.18 * ( 1.0 - cldf_ice ) + 0.35 * cldf_ice )
      fr2_i0(:,:) = ( 0.82 * ( 1.0 - cldf_ice ) + 0.65 * cldf_ice )

      CALL wrk_dealloc( jpi,jpj, zsnw )

   END SUBROUTINE usrdef_sbc_ice_flx


   !!======================================================================
END MODULE usrdef_sbc

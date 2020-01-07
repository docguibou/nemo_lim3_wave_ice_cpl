MODULE usrdef_hgr
!!======================================================================
!!                       ***  MODULE  usrdef_hgr  ***
!!
!!                     ===  SAS_BIPER configuration  ===
!!
!! User defined :   mesh and Coriolis parameter of a user configuration
!!======================================================================
!! History :  NEMO  ! 2016-08  (S. Flavoni, G. Madec)  Original code
!!----------------------------------------------------------------------

!!----------------------------------------------------------------------
!!   usr_def_hgr    : initialize the horizontal mesh for SAS_BIPER configuration
!!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE par_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE usrdef_nam, ONLY: rn_dx, rn_dy   ! horizontal resolution in meters
!
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE timing          ! Timing
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_hgr   ! called by domhgr.F90

!!----------------------------------------------------------------------
!! NEMO/OPA 4.0 , NEMO Consortium (2016)
!! $Id$
!! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_hgr( plamt , plamu , plamv  , plamf  ,   &   ! geographic position (required)
      &                    pphit , pphiu , pphiv  , pphif  ,   &   !
      &                    kff   , pff_f , pff_t  ,            &   ! Coriolis parameter  (if domain not on the sphere)
      &                    pe1t  , pe1u  , pe1v   , pe1f   ,   &   ! scale factors       (required)
      &                    pe2t  , pe2u  , pe2v   , pe2f   ,   &   !
      &                    ke1e2u_v      , pe1e2u , pe1e2v     )   ! u- & v-surfaces (if gridsize reduction is used in strait(s))
!!----------------------------------------------------------------------
!!                  ***  ROUTINE usr_def_hgr  ***
!!
!! ** Purpose :   user defined mesh and Coriolis parameter
!!
!! ** Method  :   set all intent(out) argument to a proper value
!!                SAS_BIPER configuration : uniform grid spacing (rn_dx)
!!                without Coriolis force (f=0)
!!
!! ** Action  : - define longitude & latitude of t-, u-, v- and f-points (in degrees)
!!              - define coriolis parameter at f-point if the domain in not on the sphere (on beta-plane)
!!              - define i- & j-scale factors at t-, u-, v- and f-points (in meters)
!!              - define u- & v-surfaces (if gridsize reduction is used in some straits) (in m2)
!!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   plamt, plamu, plamv, plamf   ! longitude outputs                     [degrees]
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pphit, pphiu, pphiv, pphif   ! latitude outputs                      [degrees]
      INTEGER                 , INTENT(out) ::   kff                          ! =1 Coriolis parameter computed here, =0 otherwise
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pff_f, pff_t                 ! Coriolis factor at f-point                [1/s]
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe1t, pe1u, pe1v, pe1f       ! i-scale factors                             [m]
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe2t, pe2u, pe2v, pe2f       ! j-scale factors                             [m]
      INTEGER                 , INTENT(out) ::   ke1e2u_v                     ! =1 u- & v-surfaces computed here, =0 otherwise
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe1e2u, pe1e2v               ! u- & v-surfaces (if reduction in strait)   [m2]
!
      INTEGER  ::   ji, jj   ! dummy loop indices
      REAL(wp) ::   zphi0, zlam0
      REAL(wp) ::   zti, zui, zvi, zfi   ! local scalars
      REAL(wp) ::   ztj, zuj, zvj, zfj   !   -      -
!!-------------------------------------------------------------------------------
!
      IF( nn_timing == 1 )  CALL timing_start('usr_def_hgr')
!
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_hgr : SAS_BIPER configuration bassin'
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '          f-plane with irregular grid-spacing (+- 10%)'
      IF(lwp) WRITE(numout,*) '          the max is given by rn_dx and rn_dy' 

!                          ==========
      zlam0 = 0._wp
      zphi0 = 0._wp
         


      DO jj = 1, jpj
         DO ji = 1, jpi
            zti = FLOAT( ji - 1 + nimpp - 1 )         ;   ztj = FLOAT( jj - 1 + njmpp - 1 )
            zui = FLOAT( ji - 1 + nimpp - 1 ) + 0.5   ;   zuj = FLOAT( jj - 1 + njmpp - 1 )
            zvi = FLOAT( ji - 1 + nimpp - 1 )         ;   zvj = FLOAT( jj - 1 + njmpp - 1 ) + 0.5
            zfi = FLOAT( ji - 1 + nimpp - 1 ) + 0.5   ;   zfj = FLOAT( jj - 1 + njmpp - 1 ) + 0.5

            plamt(ji,jj) = zlam0 + rn_dx * 1.e-5 * zti
            plamu(ji,jj) = zlam0 + rn_dx * 1.e-5 * zui
            plamv(ji,jj) = zlam0 + rn_dx * 1.e-5 * zvi
            plamf(ji,jj) = zlam0 + rn_dx * 1.e-5 * zfi
   
            pphit(ji,jj) = zphi0 + rn_dy * 1.e-5 * ztj
            pphiu(ji,jj) = zphi0 + rn_dy * 1.e-5 * zuj
            pphiv(ji,jj) = zphi0 + rn_dy * 1.e-5 * zvj
            pphif(ji,jj) = zphi0 + rn_dy * 1.e-5 * zfj
         END DO
      END DO
         
! Horizontal scale factors (in meters)
!                              ======
!! ==> EITHER 1) variable scale factors
         DO jj = 1, jpj
            DO ji = 1, jpi
!!pe1t(ji,jj) = rn_dx * EXP( -0.8/REAL(jpiglo**2) * (mi0(ji)-REAL(jpiglo+1)*0.5)**2 )  ! gaussian shape
!!pe2t(ji,jj) = rn_dy * EXP( -0.8/REAL(jpjglo**2) * (mj0(jj)-REAL(jpjglo+1)*0.5)**2 )  ! gaussian shape
               pe1t(ji,jj) = rn_dx * ( 1. -0.1 * ABS(REAL(mi0(ji))-REAL(jpiglo+1)*0.5) / (1.-REAL(jpiglo+1)*0.5) ) ! linear shape
               pe2t(ji,jj) = rn_dy * ( 1. -0.1 * ABS(REAL(mj0(jj))-REAL(jpjglo+1)*0.5) / (1.-REAL(jpjglo+1)*0.5) ) ! linear shape
            END DO
         END DO

!! ==> OR 2) constant scale factors
!!         pe1t(:,:) = rn_dx
!!         pe2t(:,:) = rn_dy
         
         pe1u(:,:) = pe1t(:,:)      ;      pe2u(:,:) = pe2t(:,:)
         pe1v(:,:) = pe1t(:,:)      ;      pe2v(:,:) = pe2t(:,:)
         pe1f(:,:) = pe1t(:,:)      ;      pe2f(:,:) = pe2t(:,:)

!                             ! NO reduction of grid size in some straits
      ke1e2u_v = 0                  !    ==>> u_ & v_surfaces will be computed in dom_ghr routine
      pe1e2u(:,:) = 0._wp           !    CAUTION: set to zero to avoid error with some compilers that
      pe1e2v(:,:) = 0._wp           !             require an initialization of INTENT(out) arguments
!
!
!                       !==  Coriolis parameter  ==!
      kff = 1                       !  indicate not to compute Coriolis parameter afterward
!
      pff_f(:,:) = 0._wp            ! here No earth rotation: f=0
      pff_t(:,:) = 0._wp
!
      IF( nn_timing == 1 )  CALL timing_stop('usr_def_hgr')
!
      END SUBROUTINE usr_def_hgr

!!======================================================================
END MODULE usrdef_hgr

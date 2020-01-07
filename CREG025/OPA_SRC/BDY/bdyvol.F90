MODULE bdyvol
   !!======================================================================
   !!                       ***  MODULE  bdyvol  ***
   !! Ocean dynamic :  Volume constraint when unstructured boundary 
   !!                  and filtered free surface are used
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!             -   !  2006-01  (J. Chanut) Bug correction
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!            3.6  !  2016     (J. Chanut) Generalize for all free surface options
   !!----------------------------------------------------------------------
#if defined key_bdy 
   !!----------------------------------------------------------------------
   !!   'key_bdy'      unstructured open boundary conditions
   !!----------------------------------------------------------------------
   USE timing          ! Timing
   USE oce             ! ocean dynamics and tracers 
   USE sbcisf          ! ice shelf
   USE dom_oce         ! ocean space and time domain 
   USE phycst          ! physical constants
   USE bdy_oce         ! ocean open boundary conditions
   USE lib_mpp         ! for mppsum
   USE in_out_manager  ! I/O manager
   USE sbc_oce         ! ocean surface boundary conditions
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   PUBLIC bdy_vol        ! routine called by dynspg_flt.h90

   REAL(wp), PRIVATE, SAVE ::   cflxemp   ! Total water flux at the surface

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: bdyvol.F90 5628 2015-07-22 20:26:35Z mathiot $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_vol( kt, jit, pua, pva, phu, phv )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE bdyvol  ***
      !!
      !! ** Purpose :   This routine is called in dynspg_ts (split explicit free surface) 
      !!       or in bdydyn (explicit/filtered free surface) to control the volume of the system. 
      !!      A correction velocity is calculated to correct the total transport through 
      !!      the unstructured OBC. 
      !!
      !! ** Method  :   The correction velocity (zubtpecor here) is defined calculating
      !!      the total transport through all open boundaries minus
      !!      the cumulate E-P flux (z_cflxemp) divided by the total lateral 
      !!      surface (bdysurftot) of the unstructured boundary. 
      !!         zubtpecor = [trans_bdy - z_cflxemp ]*(1./bdysurftot)
      !!      with z_cflxemp => sum of (Evaporation minus Precipitation)
      !!                       over all the domain in m3/s at each time step.
      !!      z_cflxemp < 0 when precipitation dominate
      !!      z_cflxemp > 0 when evaporation dominate
      !!
      !!      There are 2 options (user's desiderata): 
      !!         1/ The volume changes according to E-P, this is the default
      !!            option. In this case the cumulate E-P flux are setting to
      !!            zero (z_cflxemp=0) to calculate the correction velocity. So
      !!            it will only balance the flux through open boundaries.
      !!            (set nn_volctl to 0 in tne namelist for this option)
      !!         2/ The volume is constant even with E-P flux. In this case
      !!            the correction velocity must balance both the flux 
      !!            through open boundaries and the ones through the free
      !!            surface. 
      !!            (set nn_volctl to 1 in tne namelist for this option)
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   		! ocean baroclinic time-step
      INTEGER, INTENT( in ), OPTIONAL ::   jit  ! ocean barotropic time-step
      REAL(wp), DIMENSION(:,:), INTENT(inout) :: pua , pva  ! barotropic velocities
      REAL(wp), DIMENSION(:,:), INTENT(in   ) :: phu , phv  ! Total depths
      !!
      LOGICAL  ::   lk_2d, lk_fstp
      INTEGER  ::   ji, jj, jk, jb, jgrd
      INTEGER  ::   ib_bdy, ii, ij
      REAL(wp) ::   zubtpecor, ztranst
      TYPE(OBC_INDEX), POINTER :: idx
      !!---------------------------------------------------------------------------

      IF (.NOT.ln_vol) RETURN

      IF( nn_timing == 1 ) CALL timing_start('bdy_vol')

      IF( PRESENT(jit) ) THEN ; lk_2d = .TRUE.; ELSE ; lk_2d = .FALSE. ; ENDIF
      lk_fstp =.FALSE.
      IF (( lk_2d.AND.(kt==nit000).AND.(jit==1) ).OR.(.NOT.lk_2d.AND.(kt==nit000) )) lk_fstp = .TRUE.

      IF( lk_fstp ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*)'bdy_vol : Correction of velocities along unstructured OBC'
         IF(lwp) WRITE(numout,*)'~~~~~~~'
      END IF 

      ! Calculate the cumulate surface Flux z_cflxemp (m3/s) over all the domain
      ! ------------------------------------------------------------------------
      IF( ( nn_volctl==1 ).AND.( (lk_2d.AND.jit==1).OR.(.NOT.lk_2d) ) ) THEN ! First external time step only in 2d case
         !
         IF ( lk_2d.AND.ln_bt_fw ) THEN
            cflxemp = SUM ( ( emp(:,:)-rnf(:,:)+fwfisf(:,:) ) * bdytmask(:,:) * tmask_i(:,:) * e12t(:,:) ) / rau0
         ELSE
            cflxemp = SUM ( ( emp(:,:) + emp_b(:,:) - rnf(:,:) - rnf_b(:,:) & 
                      &    + fwfisf(:,:) + fwfisf_b(:,:) ) * bdytmask(:,:) * tmask_i(:,:) * e12t(:,:) ) * 0.5_wp / rau0
         ENDIF
         IF( lk_mpp )   CALL mpp_sum( cflxemp )     ! sum over the global domain
         !
      ENDIF

      ! Transport through the unstructured open boundary
      ! ------------------------------------------------
      zubtpecor = 0._wp

      DO ib_bdy = 1, nb_bdy
         idx => idx_bdy(ib_bdy)
!CT add volume correction {
!CT Compute the transport through all BDYs, remove U component for safety reason only
!CT since the flagu parameter should be zero for both northern and southern BDYs
!CT      jgrd = 2                               ! cumulate u component contribution first 
!CT      DO jb = 1, idx%nblenrim(jgrd)
!CT         ii = idx%nbi(jb,jgrd)
!CT         ij = idx%nbj(jb,jgrd)
!CT         zubtpecor = zubtpecor + idx%flagu(jb,jgrd) * pua(ii,ij)    & 
!CT                   &             * e2u(ii,ij) * phu(ii,ij) * tmask_i(ii,ij) * tmask_i(ii+1,ij)
!CT      END DO
         jgrd = 3                               ! then add v component contribution
         DO jb = 1, idx%nblenrim(jgrd)
            ii = idx%nbi(jb,jgrd)
            ij = idx%nbj(jb,jgrd)
            zubtpecor = zubtpecor + idx%flagv(jb,jgrd) * pva(ii,ij)    &
                      &             * e1v(ii,ij) * phv(ii,ij) * tmask_i(ii,ij) * tmask_i(ii,ij+1) 
         END DO
      END DO
!CT add volume correction }

      IF( lk_mpp )   CALL mpp_sum( zubtpecor )   ! sum over the global domain

      ! The normal velocity correction
      ! ------------------------------
      IF( nn_volctl==1 ) THEN   ;   zubtpecor = ( zubtpecor - cflxemp) / bdysurftot 
      ELSE                      ;   zubtpecor =   zubtpecor            / bdysurftot
      END IF

      ! Correction of the barotropic velocity on the unstructured boundary to respect the mass flux conservation
      ! --------------------------------------------------------------------------------------------------------
      ztranst = 0._wp

!CT add volume correction {
!CT   DO ib_bdy = 1, nb_bdy
!CT Apply the volume correction to the southern BDY only
      DO ib_bdy = 1, 1
         idx => idx_bdy(ib_bdy)
!CT      jgrd = 2                               ! correct u component
!CT      DO jb = 1, idx%nblenrim(jgrd)
!CT         ii = idx%nbi(jb,jgrd)
!CT         ij = idx%nbj(jb,jgrd)
!CT         pua(ii,ij) = pua(ii,ij) - idx%flagu(jb,jgrd) * zubtpecor * umask(ii,ij,1)
!CT         ztranst = ztranst + idx%flagu(jb,jgrd) * pua(ii,ij) * e2u(ii,ij) * phu(ii,ij) &
!CT                 &           * tmask_i(ii,ij) * tmask_i(ii+1,ij) 
!CT      END DO
         jgrd = 3                               ! correct v component
         DO jb = 1, idx%nblenrim(jgrd)
            ii = idx%nbi(jb,jgrd)
            ij = idx%nbj(jb,jgrd)
            pva(ii,ij) = pva(ii,ij) - idx%flagv(jb,jgrd) * zubtpecor * vmask(ii,ij,1)
            ztranst = ztranst + idx%flagv(jb,jgrd) * pva(ii,ij) * e1v(ii,ij) * phv(ii,ij) &
                       &           * tmask_i(ii,ij) * tmask_i(ii,ij+1)            
         END DO
!CT add volume correction }

         CALL lbc_bdy_lnk( pua, 'U', -1., ib_bdy ) 
         CALL lbc_bdy_lnk( pva, 'V', -1., ib_bdy)   ! Boundary points should be updated
      END DO

      IF( lk_mpp )   CALL mpp_sum( ztranst )   ! sum over the global domain
 
      ! Check the cumulated transport through unstructured OBC once barotropic velocities corrected
      ! ------------------------------------------------------
      IF (( lwp .AND. MOD( kt, nwrite ) == 0).AND.(lk_2d.AND.(jit==1))) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*)'bdy_vol : time step :', kt
         IF(lwp) WRITE(numout,*)'~~~~~~~ '
         IF(lwp) WRITE(numout,*)'          cumulate flux EMP             =', cflxemp   , '(m3/s)'
         IF(lwp) WRITE(numout,*)'          total lateral surface of OBC  =', bdysurftot, '(m2)'
         IF(lwp) WRITE(numout,*)'          correction velocity zubtpecor =', zubtpecor , '(m/s)'
         IF(lwp) WRITE(numout,*)'          cumulated transport ztranst   =', ztranst   , '(m3/s)'
      END IF 
      !
      IF( nn_timing == 1 ) CALL timing_stop('bdy_vol')
      !

   END SUBROUTINE bdy_vol

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                   NO Unstruct Open Boundary Conditions
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_vol( kt )        ! Empty routine
      INTEGER, INTENT( in ) ::   kt   		! ocean baroclinic time-step
      WRITE(*,*) 'bdy_vol: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_vol
#endif

   !!======================================================================
END MODULE bdyvol

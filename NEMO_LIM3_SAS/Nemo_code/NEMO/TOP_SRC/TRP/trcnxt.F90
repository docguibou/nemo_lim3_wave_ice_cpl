MODULE trcnxt
   !!======================================================================
   !!                       ***  MODULE  trcnxt  ***
   !! Ocean passive tracers:  time stepping on passives tracers
   !!======================================================================
   !! History :  7.0  !  1991-11  (G. Madec)  Original code
   !!                 !  1993-03  (M. Guyon)  symetrical conditions
   !!                 !  1995-02  (M. Levy)   passive tracers
   !!                 !  1996-02  (G. Madec & M. Imbard)  opa release 8.0
   !!            8.0  !  1996-04  (A. Weaver)  Euler forward step
   !!            8.2  !  1999-02  (G. Madec, N. Grima)  semi-implicit pressure grad.
   !!  NEMO      1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!                 !  2002-08  (G. Madec)  F90: Free form and module
   !!                 !  2002-11  (C. Talandier, A-M Treguier) Open boundaries
   !!                 !  2004-03  (C. Ethe) passive tracers
   !!                 !  2007-02  (C. Deltel) Diagnose ML trends for passive tracers
   !!            2.0  !  2006-02  (L. Debreu, C. Mazauric) Agrif implementation
   !!            3.0  !  2008-06  (G. Madec)  time stepping always done in trazdf
   !!            3.1  !  2009-02  (G. Madec, R. Benshila)  re-introduce the vvl option
   !!            3.3  !  2010-06  (C. Ethe, G. Madec) Merge TRA-TRC
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_nxt     : time stepping on passive tracers
   !!----------------------------------------------------------------------
   USE oce_trc         ! ocean dynamics and tracers variables
   USE trc             ! ocean passive tracers variables
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl_trc      ! Print control for debbuging
   USE trd_oce
   USE trdtra
   USE tranxt
   USE bdy_oce   , ONLY: ln_bdy
   USE trcbdy          ! BDY open boundaries
# if defined key_agrif
   USE agrif_top_interp
# endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nxt          ! routine called by step.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnxt.F90 7753 2017-03-03 11:46:59Z mocavero $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_nxt( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trcnxt  ***
      !!
      !! ** Purpose :   Compute the passive tracers fields at the 
      !!      next time-step from their temporal trends and swap the fields.
      !! 
      !! ** Method  :   Apply lateral boundary conditions on (ua,va) through 
      !!      call to lbc_lnk routine
      !!   default:
      !!      arrays swap
      !!         (trn) = (tra) ; (tra) = (0,0)
      !!         (trb) = (trn) 
      !!
      !!   For Arakawa or TVD Scheme : 
      !!      A Asselin time filter applied on now tracers (trn) to avoid
      !!      the divergence of two consecutive time-steps and tr arrays
      !!      to prepare the next time_step:
      !!         (trb) = (trn) + atfp [ (trb) + (tra) - 2 (trn) ]
      !!         (trn) = (tra) ; (tra) = (0,0)
      !!
      !!
      !! ** Action  : - update trb, trn
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step index
      !
      INTEGER  ::   jk, jn   ! dummy loop indices
      REAL(wp) ::   zfact            ! temporary scalar
      CHARACTER (len=22) :: charout
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::  ztrdt 
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_nxt')
      !
      IF( kt == nittrc000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_nxt : time stepping on passive tracers'
      ENDIF
      !
#if defined key_agrif
      CALL Agrif_trc                   ! AGRIF zoom boundaries
#endif
      DO jn = 1, jptra                 ! Update after tracer on domain lateral boundaries
         CALL lbc_lnk( tra(:,:,:,jn), 'T', 1. )   
      END DO

      IF( ln_bdy )  CALL trc_bdy( kt )

      IF( l_trdtrc )  THEN             ! trends: store now fields before the Asselin filter application
         CALL wrk_alloc( jpi, jpj, jpk, jptra, ztrdt )
         ztrdt(:,:,:,:)  = trn(:,:,:,:)
      ENDIF
      !                                ! Leap-Frog + Asselin filter time stepping
      IF( neuler == 0 .AND. kt == nittrc000 ) THEN    ! Euler time-stepping at first time-step (only swap)
         DO jn = 1, jptra
            DO jk = 1, jpkm1
               trn(:,:,jk,jn) = tra(:,:,jk,jn)
            END DO
         END DO
      ELSE                                            ! Asselin filter + swap
         IF( ln_linssh ) THEN   ;   CALL tra_nxt_fix( kt, nittrc000,         'TRC', trb, trn, tra, jptra )  !     linear ssh
         ELSE                   ;   CALL tra_nxt_vvl( kt, nittrc000, rdttrc, 'TRC', trb, trn, tra,      &
           &                                                                   sbc_trc, sbc_trc_b, jptra )  ! non-linear ssh
         ENDIF
         !
         DO jn = 1, jptra
            CALL lbc_lnk( trb(:,:,:,jn), 'T', 1._wp ) 
            CALL lbc_lnk( trn(:,:,:,jn), 'T', 1._wp )
            CALL lbc_lnk( tra(:,:,:,jn), 'T', 1._wp )
         END DO
      ENDIF
      !
      IF( l_trdtrc ) THEN              ! trends: send Asselin filter trends to trdtra manager for further diagnostics
         DO jn = 1, jptra
            DO jk = 1, jpkm1
               zfact = 1._wp / r2dttrc  
               ztrdt(:,:,jk,jn) = ( trb(:,:,jk,jn) - ztrdt(:,:,jk,jn) ) * zfact 
               CALL trd_tra( kt, 'TRC', jn, jptra_atf, ztrdt )
            END DO
         END DO
         CALL wrk_dealloc( jpi, jpj, jpk, jptra, ztrdt ) 
      END IF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('nxt')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=trn, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_nxt')
      !
   END SUBROUTINE trc_nxt

#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nxt( kt )  
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_nxt: You should not have seen this print! error?', kt
   END SUBROUTINE trc_nxt
#endif
   !!======================================================================
END MODULE trcnxt

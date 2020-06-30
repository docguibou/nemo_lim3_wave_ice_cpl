#define SPONGE_TOP

MODULE agrif_top_sponge
   !!======================================================================
   !!                ***  MODULE agrif_top_sponge  ***
   !! AGRIF :   define in memory AGRIF variables for sea-ice
   !!======================================================================
   !! History :  2.0  ! 2006-08  (R. Benshila, L. Debreu)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   Agrif_Sponge_trc : 
   !!   interptrn_sponge :  
   !!----------------------------------------------------------------------
#if defined key_agrif && defined key_top
   USE par_oce
   USE par_trc
   USE oce
   USE trc
   USE dom_oce
   USE agrif_oce
   USE agrif_opa_sponge
   !
   USE in_out_manager
   USE lib_mpp
   USE wrk_nemo  

   IMPLICIT NONE
   PRIVATE

   PUBLIC Agrif_Sponge_trc, interptrn_sponge

   !!----------------------------------------------------------------------
   !! NEMO/NST 3.7 , NEMO Consortium (2015)
   !! $Id: agrif_top_sponge.F90 6140 2015-12-21 11:35:23Z timgraham $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE Agrif_Sponge_trc
      !!----------------------------------------------------------------------
      !!                   *** ROUTINE Agrif_Sponge_Trc ***
      !!----------------------------------------------------------------------
      REAL(wp) ::   timecoeff
      !!----------------------------------------------------------------------
      !
#if defined SPONGE_TOP
      timecoeff = REAL( Agrif_NbStepint(), wp ) / Agrif_rhot()
      CALL Agrif_sponge
      Agrif_SpecialValue    = 0._wp
      Agrif_UseSpecialValue = .TRUE.
      tabspongedone_trn     = .FALSE.
      CALL Agrif_Bc_Variable( trn_sponge_id, calledweight=timecoeff, procname=interptrn_sponge )
      Agrif_UseSpecialValue = .FALSE.
#endif
      !
   END SUBROUTINE Agrif_Sponge_Trc


   SUBROUTINE interptrn_sponge( tabres, i1, i2, j1, j2, k1, k2, n1, n2, before )
      !!----------------------------------------------------------------------
      !!                   *** ROUTINE interptrn_sponge ***
      !!----------------------------------------------------------------------
      INTEGER                                     , INTENT(in   ) ::   i1, i2, j1, j2, k1, k2, n1, n2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2,n1:n2), INTENT(inout) ::   tabres
      LOGICAL                                     , INTENT(in   ) ::   before
      !
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      REAL(wp) ::   zabe1, zabe2
      REAL(wp), DIMENSION(i1:i2,j1:j2)             ::   ztu, ztv
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2,n1:n2) ::   trbdiff
      !!----------------------------------------------------------------------
      !
      IF( before ) THEN
         tabres(i1:i2,j1:j2,k1:k2,n1:n2) = trn(i1:i2,j1:j2,k1:k2,n1:n2)
      ELSE      
!!gm line below use of :,:  versus i1:i2,j1:j2  ....   strange, not wrong.    ===>> to be corrected
         trbdiff(:,:,:,:) = trb(i1:i2,j1:j2,:,:) - tabres(:,:,:,:)      
         DO jn = 1, jptra
            DO jk = 1, jpkm1
               DO jj = j1,j2-1
                  DO ji = i1,i2-1
                     zabe1 = fsaht_spu(ji,jj) * e2_e1u(ji,jj) * e3u_n(ji,jj,jk) * umask(ji,jj,jk)
                     zabe2 = fsaht_spv(ji,jj) * e1_e2v(ji,jj) * e3v_n(ji,jj,jk) * vmask(ji,jj,jk)
                     ztu(ji,jj) = zabe1 * ( trbdiff(ji+1,jj  ,jk,jn) - trbdiff(ji,jj,jk,jn) )
                     ztv(ji,jj) = zabe2 * ( trbdiff(ji  ,jj+1,jk,jn) - trbdiff(ji,jj,jk,jn) )
                  END DO
               END DO
               !
               DO jj = j1+1,j2-1
                  DO ji = i1+1,i2-1
                     IF( .NOT. tabspongedone_trn(ji,jj) ) THEN 
                        tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + (  ztu(ji,jj) - ztu(ji-1,jj  )     &
                           &                                   + ztv(ji,jj) - ztv(ji  ,jj-1)  )  &
                           &                                * r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk)
                     ENDIF
                  END DO
               END DO
            END DO
            !
         END DO
         !
         tabspongedone_trn(i1+1:i2-1,j1+1:j2-1) = .TRUE.
      ENDIF
      !                 
   END SUBROUTINE interptrn_sponge

#else

CONTAINS
   SUBROUTINE agrif_top_sponge_empty
      WRITE(*,*)  'agrif_top_sponge : You should not have seen this print! error?'
   END SUBROUTINE agrif_top_sponge_empty
#endif

   !!======================================================================
END MODULE agrif_top_sponge

#define TWO_WAY
#undef DECAL_FEEDBACK

MODULE agrif_top_update
   !!======================================================================
   !!                ***  MODULE agrif_top_update  ***
   !! AGRIF :   
   !!----------------------------------------------------------------------
   !! History :  
   !!----------------------------------------------------------------------

#if defined key_agrif && defined key_top
   USE par_oce
   USE oce
   USE par_trc
   USE trc
   USE dom_oce
   USE agrif_oce
   USE wrk_nemo  

   IMPLICIT NONE
   PRIVATE

   PUBLIC Agrif_Update_Trc

   INTEGER, PUBLIC ::   nbcline_trc = 0   !: ???

   !!----------------------------------------------------------------------
   !! NEMO/NST 3.7 , NEMO Consortium (2015)
   !! $Id: agrif_top_update.F90 6140 2015-12-21 11:35:23Z timgraham $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE Agrif_Update_Trc( kt )
      !!----------------------------------------------------------------------
      !!                   *** ROUTINE Agrif_Update_Trc ***
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt
      !!----------------------------------------------------------------------
      ! 
      IF((Agrif_NbStepint() .NE. (Agrif_irhot()-1)).AND.(kt /= 0)) RETURN
#if defined TWO_WAY   
      Agrif_UseSpecialValueInUpdate = .TRUE.
      Agrif_SpecialValueFineGrid    = 0._wp
      ! 
      IF( MOD(nbcline_trc,nbclineupdate) == 0 ) THEN
# if ! defined DECAL_FEEDBACK
         CALL Agrif_Update_Variable(trn_id, procname=updateTRC )
# else
         CALL Agrif_Update_Variable(trn_id, locupdate=(/1,0/),procname=updateTRC )
# endif
      ELSE
# if ! defined DECAL_FEEDBACK
         CALL Agrif_Update_Variable( trn_id, locupdate=(/0,2/), procname=updateTRC )
# else
         CALL Agrif_Update_Variable( trn_id, locupdate=(/1,2/), procname=updateTRC )
# endif
      ENDIF
      !
      Agrif_UseSpecialValueInUpdate = .FALSE.
      nbcline_trc = nbcline_trc + 1
#endif
      !
   END SUBROUTINE Agrif_Update_Trc


   SUBROUTINE updateTRC( ptab, i1, i2, j1, j2, k1, k2, n1, n2, before )
      !!----------------------------------------------------------------------
      !!                      *** ROUTINE updateT ***
      !!----------------------------------------------------------------------
      INTEGER                                    , INTENT(in   ) ::   i1, i2, j1, j2, k1, k2, n1, n2
      REAL(wp),DIMENSION(i1:i2,j1:j2,k1:k2,n1:n2), INTENT(inout) ::   ptab
      LOGICAL                                    , INTENT(in   ) ::   before
      !!
      INTEGER ::   ji, jj, jk, jn
      !!----------------------------------------------------------------------
      !
      IF( before ) THEN
         ptab(i1:i2,j1:j2,k1:k2,n1:n2) = trn(i1:i2,j1:j2,k1:k2,n1:n2)
      ELSE
         IF( .NOT.(lk_agrif_fstep.AND.(neuler==0)) ) THEN
            ! Add asselin part
            DO jn = n1,n2
               DO jk = k1, k2
                  DO jj = j1, j2
                     DO ji = i1, i2
                        IF( ptab(ji,jj,jk,jn) /= 0._wp ) THEN
                           trb(ji,jj,jk,jn) = trb(ji,jj,jk,jn)             & 
                              &             + atfp * ( ptab(ji,jj,jk,jn)   &
                                 &                    - trn(ji,jj,jk,jn) ) * tmask(ji,jj,jk)
                        ENDIF
                     END DO
                  END DO
               END DO
            END DO
         ENDIF
         DO jn = n1, n2
            DO jk = k1, k2
               DO jj = j1, j2
                  DO ji = i1, i2
                     IF( ptab(ji,jj,jk,jn) /= 0._wp ) THEN 
                        trn(ji,jj,jk,jn) = ptab(ji,jj,jk,jn) * tmask(ji,jj,jk)
                     END IF
                  END DO
               END DO
            END DO
         END DO
      ENDIF
      ! 
   END SUBROUTINE updateTRC

#else
CONTAINS
   SUBROUTINE agrif_top_update_empty
      !!---------------------------------------------
      !!   *** ROUTINE agrif_Top_update_empty ***
      !!---------------------------------------------
      WRITE(*,*)  'agrif_top_update : You should not have seen this print! error?'
   END SUBROUTINE agrif_top_update_empty
#endif

   !!======================================================================
END MODULE agrif_top_update

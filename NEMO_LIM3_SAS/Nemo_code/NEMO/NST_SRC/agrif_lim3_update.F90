#define TWO_WAY

MODULE agrif_lim3_update
   !!=====================================================================================
   !!                       ***  MODULE agrif_lim3_update ***
   !! Nesting module :  update surface ocean boundary condition over ice from a child grid
   !! Sea-Ice model  :  LIM 3.6 Sea ice model time-stepping
   !!=====================================================================================
   !! History :  2.0   !  04-2008  (F. Dupont)  initial version
   !!            3.4   !  08-2012  (R. Benshila, C. Herbaut) update and EVP
   !!            3.6   !  05-2016  (C. Rousset)  Add LIM3 compatibility
   !!----------------------------------------------------------------------
#if defined key_agrif && defined key_lim3
   !!----------------------------------------------------------------------
   !!   'key_lim3'  :                                 LIM 3.6 sea-ice model
   !!   'key_agrif' :                                 AGRIF library 
   !!----------------------------------------------------------------------
   !!   agrif_update_lim3  : update sea-ice on boundaries or total
   !!                        child domain for velocities and ice properties
   !!   update_tra_ice     : sea-ice properties
   !!   update_u_ice       : zonal      ice velocity
   !!   update_v_ice       : meridional ice velocity
   !!----------------------------------------------------------------------
   USE dom_oce
   USE sbc_oce
   USE agrif_oce
   USE ice
   USE agrif_ice 

   IMPLICIT NONE
   PRIVATE

   PUBLIC agrif_update_lim3

   !!----------------------------------------------------------------------
   !! NEMO/NST 3.6 , LOCEAN-IPSL (2016)
   !! $Id: agrif_lim3_update.F90 6204 2016-01-04 13:47:06Z cetlod $
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE agrif_update_lim3( kt )
      !!----------------------------------------------------------------------
      !!                     *** ROUTINE agrif_update_lim3 ***
      !! ** Method  :   Call the hydrostaticupdate pressure at the boundary or the entire domain 
      !!
      !! ** Action : - Update (u_ice,v_ice) and ice tracers
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt
      !!
      !!----------------------------------------------------------------------
      !
      !! clem: I think the update should take place each time the ocean sees the surface forcings
      !!       (but maybe I am wrong and we should update every rhot time steps) 
      IF( ( MOD( (kt-nit000)/nn_fsbc + 1, Agrif_irhot() * Agrif_Parent(nn_fsbc) / nn_fsbc ) /=0 ) .AND. (kt /= 0) ) RETURN ! do not update if nb of child time steps differ from time refinement
                                                                                                                           ! i.e. update only at the parent time step
      Agrif_UseSpecialValueInUpdate = .TRUE.
      Agrif_SpecialValueFineGrid = -9999.
# if defined TWO_WAY
      IF( MOD(nbcline,nbclineupdate) == 0) THEN ! update the whole basin at each nbclineupdate (=nn_cln_update) baroclinic parent time steps
                                                ! nbcline is incremented (+1) at the end of each parent time step from 0 (1st time step)
         CALL Agrif_Update_Variable( tra_ice_id , procname = update_tra_ice  )
         CALL Agrif_Update_Variable( u_ice_id   , procname = update_u_ice    )
         CALL Agrif_Update_Variable( v_ice_id   , procname = update_v_ice    )
      ELSE                                      ! update only the boundaries defined par locupdate
         CALL Agrif_Update_Variable( tra_ice_id , locupdate=(/0,2/), procname = update_tra_ice  )
         CALL Agrif_Update_Variable( u_ice_id   , locupdate=(/0,1/), procname = update_u_ice    )
         CALL Agrif_Update_Variable( v_ice_id   , locupdate=(/0,1/), procname = update_v_ice    )
      ENDIF
# endif
      Agrif_UseSpecialValueInUpdate = .FALSE.
      !
   END SUBROUTINE agrif_update_lim3


   !!------------------
   !! Local subroutines
   !!------------------
   SUBROUTINE update_tra_ice( ptab, i1, i2, j1, j2, k1, k2, before )
      !!-----------------------------------------------------------------------
      !!                        *** ROUTINE update_tra_ice ***
      !! ** Method  : Compute the mass properties on the fine grid and recover
      !!              the properties per mass on the coarse grid
      !!-----------------------------------------------------------------------
      INTEGER , INTENT(in) :: i1, i2, j1, j2, k1, k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: ptab
      LOGICAL , INTENT(in) :: before
      !!
      INTEGER  :: jk, jl, jm
      !!-----------------------------------------------------------------------
      ! it is ok not to multiply by e1*e2 since we conserve tracers here (same as in the ocean).
      IF( before ) THEN
         jm = 1
         DO jl = 1, jpl
            ptab(:,:,jm) = a_i  (i1:i2,j1:j2,jl) ; jm = jm + 1
            ptab(:,:,jm) = v_i  (i1:i2,j1:j2,jl) ; jm = jm + 1
            ptab(:,:,jm) = v_s  (i1:i2,j1:j2,jl) ; jm = jm + 1
            ptab(:,:,jm) = smv_i(i1:i2,j1:j2,jl) ; jm = jm + 1
            ptab(:,:,jm) = oa_i (i1:i2,j1:j2,jl) ; jm = jm + 1
            DO jk = 1, nlay_s
               ptab(:,:,jm) = e_s(i1:i2,j1:j2,jk,jl) ; jm = jm + 1
            ENDDO
            DO jk = 1, nlay_i
               ptab(:,:,jm) = e_i(i1:i2,j1:j2,jk,jl) ; jm = jm + 1
            ENDDO
         ENDDO

         DO jk = k1, k2
            WHERE( tmask(i1:i2,j1:j2,1) == 0. )  ptab(:,:,jk) = -9999.
         ENDDO
                  
      ELSE
         jm = 1
         DO jl = 1, jpl
            a_i  (i1:i2,j1:j2,jl) = ptab(:,:,jm) * tmask(i1:i2,j1:j2,1) ; jm = jm + 1
            v_i  (i1:i2,j1:j2,jl) = ptab(:,:,jm) * tmask(i1:i2,j1:j2,1) ; jm = jm + 1
            v_s  (i1:i2,j1:j2,jl) = ptab(:,:,jm) * tmask(i1:i2,j1:j2,1) ; jm = jm + 1
            smv_i(i1:i2,j1:j2,jl) = ptab(:,:,jm) * tmask(i1:i2,j1:j2,1) ; jm = jm + 1
            oa_i (i1:i2,j1:j2,jl) = ptab(:,:,jm) * tmask(i1:i2,j1:j2,1) ; jm = jm + 1
            DO jk = 1, nlay_s
               e_s(i1:i2,j1:j2,jk,jl) = ptab(:,:,jm) * tmask(i1:i2,j1:j2,1) ; jm = jm + 1
            ENDDO
            DO jk = 1, nlay_i
               e_i(i1:i2,j1:j2,jk,jl) = ptab(:,:,jm) * tmask(i1:i2,j1:j2,1) ; jm = jm + 1
            ENDDO
         ENDDO

         ! integrated values
         vt_i (i1:i2,j1:j2) = SUM( v_i(i1:i2,j1:j2,:), dim=3 )
         vt_s (i1:i2,j1:j2) = SUM( v_s(i1:i2,j1:j2,:), dim=3 )
         at_i (i1:i2,j1:j2) = SUM( a_i(i1:i2,j1:j2,:), dim=3 )
         et_s(i1:i2,j1:j2)  = SUM( SUM( e_s(i1:i2,j1:j2,:,:), dim=4 ), dim=3 )
         et_i(i1:i2,j1:j2)  = SUM( SUM( e_i(i1:i2,j1:j2,:,:), dim=4 ), dim=3 )
         
      ENDIF
      !
   END SUBROUTINE update_tra_ice


   SUBROUTINE update_u_ice( ptab, i1, i2, j1, j2, before )
      !!-----------------------------------------------------------------------
      !!                        *** ROUTINE update_u_ice ***
      !! ** Method  : Update the fluxes and recover the properties (C-grid)
      !!-----------------------------------------------------------------------
      INTEGER , INTENT(in) :: i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: ptab
      LOGICAL , INTENT(in) :: before
      !!
      REAL(wp) :: zrhoy
      !!-----------------------------------------------------------------------
      !
      IF( before ) THEN
         zrhoy = Agrif_Rhoy()
         ptab(:,:) = e2u(i1:i2,j1:j2) * u_ice(i1:i2,j1:j2) * zrhoy
         WHERE( umask(i1:i2,j1:j2,1) == 0. )  ptab(:,:) = -9999.
      ELSE
         u_ice(i1:i2,j1:j2) = ptab(:,:) / e2u(i1:i2,j1:j2) * umask(i1:i2,j1:j2,1)
      ENDIF
      ! 
   END SUBROUTINE update_u_ice


   SUBROUTINE update_v_ice( ptab, i1, i2, j1, j2, before )
      !!-----------------------------------------------------------------------
      !!                    *** ROUTINE update_v_ice ***
      !! ** Method  : Update the fluxes and recover the properties (C-grid)
      !!-----------------------------------------------------------------------
      INTEGER , INTENT(in) :: i1,i2,j1,j2
      REAL(wp), DIMENSION(i1:i2,j1:j2),  INTENT(inout) :: ptab
      LOGICAL , INTENT(in) :: before
      !!
      REAL(wp) :: zrhox
      !!-----------------------------------------------------------------------
      !
      IF( before ) THEN
         zrhox = Agrif_Rhox()
         ptab(:,:) = e1v(i1:i2,j1:j2) * v_ice(i1:i2,j1:j2) * zrhox
         WHERE( vmask(i1:i2,j1:j2,1) == 0. )  ptab(:,:) = -9999.
      ELSE
         v_ice(i1:i2,j1:j2) = ptab(:,:) / e1v(i1:i2,j1:j2) * vmask(i1:i2,j1:j2,1)
      ENDIF
      !
   END SUBROUTINE update_v_ice

#else
CONTAINS
   SUBROUTINE agrif_lim3_update_empty
      !!---------------------------------------------
      !!   *** ROUTINE agrif_lim3_update_empty ***
      !!---------------------------------------------
      WRITE(*,*)  'agrif_lim3_update : You should not have seen this print! error?'
   END SUBROUTINE agrif_lim3_update_empty
#endif
END MODULE agrif_lim3_update

MODULE step
!!======================================================================
!!                       ***  MODULE step  ***
!! Time-stepping    : manager of the ocean, tracer and ice time stepping
!!                    version for standalone surface scheme
!!======================================================================
!! History :  OPA  !  1991-03  (G. Madec)  Original code
!!             .   !    .
!!             .   !    .
!!   NEMO     3.5  !  2012-03  (S. Alderson)
!!----------------------------------------------------------------------

!!----------------------------------------------------------------------
!!   stp             : OPA system time-stepping
!!----------------------------------------------------------------------
   USE oce              ! ocean dynamics and tracers variables
   USE dom_oce          ! ocean space and time domain variables
   USE daymod           ! calendar                         (day     routine)
   USE sbc_oce          ! surface boundary condition: fields
   USE sbcmod           ! surface boundary condition       (sbc     routine)
   USE sbcrnf           ! surface boundary condition: runoff variables
   USE sbccpl           ! surface boundary condition: coupled interface
   USE eosbn2           ! equation of state                (eos_bn2 routine)
   USE diawri           ! Standard run outputs             (dia_wri routine)
   USE bdy_oce   , ONLY: ln_bdy
   USE bdydta           ! clem: mandatory for LIM3
   USE stpctl           ! time stepping control            (stp_ctl routine)
!
   USE in_out_manager   ! I/O manager
   USE prtctl           ! Print control                    (prt_ctl routine)
   USE iom              !
   USE lbclnk           !
   USE timing           ! Timing

   USE xios



   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   stp   ! called by nemogcm.F90

!!----------------------------------------------------------------------
!! NEMO/OPA 3.3 , NEMO Consortium (2010)
!! $Id: step.F90 7761 2017-03-06 17:58:35Z clem $
!! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE stp( kstp )
      INTEGER, INTENT(in) ::   kstp   ! ocean time-step index

!!----------------------------------------------------------------------
!!                     ***  ROUTINE stp  ***
!!
!! ** Purpose : - Time stepping of SBC (surface boundary)
!!
!! ** Method  : -1- Update forcings and data
!!              -2- Outputs and diagnostics
!!----------------------------------------------------------------------
      INTEGER ::   indic    ! error indicator if < 0
!! ---------------------------------------------------------------------


                             indic = 0                    ! although indic is not changed in stp_ctl
! need to keep the same interface
      IF( kstp == nit000 )   CALL iom_init( cxios_context ) ! iom_put initialization (must be done after nemo_init for AGRIF+XIOS+OASIS)
      IF( kstp /= nit000 )   CALL day( kstp )             ! Calendar (day was already called at nit000 in day_init)
                             CALL iom_setkt( kstp - nit000 + 1, cxios_context )   ! tell iom we are at time step kstp

! ==> clem: open boundaries is mandatory for LIM3 because ice BDY is not decoupled from
!           the environment of ocean BDY. Therefore bdy is called in both OPA and SAS modules.
!           From SAS: ocean bdy data are wrong  (but we do not care) and ice bdy data are OK.
!           This is not clean and should be changed in the future.
      IF( ln_bdy     )       CALL bdy_dta ( kstp, time_offset=+1 )   ! update dynamic & tracer data at open boundaries
! ==>
                             CALL sbc    ( kstp )         ! Sea Boundary Condition (including sea-ice)

                             CALL dia_wri( kstp )         ! ocean model: outputs


                             
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! Control
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                             CALL stp_ctl( kstp, indic )
      IF( indic < 0  )  THEN
                             CALL ctl_stop( 'step: indic < 0' )
                             CALL dia_wri_state( 'output.abort', kstp )
      ENDIF
      IF( kstp == nit000   ) CALL iom_close( numror )     ! close input  ocean restart file (clem: not sure...)
      
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! Coupled mode
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( lk_oasis    )  CALL sbc_cpl_snd( kstp )     ! coupled mode : field exchanges if OASIS-coupled ice


      IF( kstp == nitend .OR. indic < 0 ) THEN 
                             CALL iom_context_finalize( cxios_context ) ! needed for XIOS+AGRIF
      ENDIF

!
      IF( nn_timing == 1 .AND.  kstp == nit000  )   CALL timing_reset
!
   END SUBROUTINE stp

!!======================================================================
END MODULE step

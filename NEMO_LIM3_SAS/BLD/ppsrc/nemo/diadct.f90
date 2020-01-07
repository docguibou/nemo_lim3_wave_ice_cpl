MODULE diadct
!!======================================================================
!!                       ***  MODULE  diadct  ***
!! Ocean diagnostics: Compute the transport trough a sec.
!!======================================================================
!! History :  OPA  ! 02/1999 (Y Drillet)  original code
!!                 ! 10/2001 (Y Drillet, R Bourdalle Badie)
!!   NEMO     1.0  ! 10/2005 (M Laborie) F90
!!            3.0  ! 04/2007 (G Garric) Ice sections
!!             -   ! 04/2007 (C Bricaud) test on sec%nb_point, initialisation of ztransp1,ztransp2,...
!!            3.4  ! 09/2011 (C Bricaud)
!!----------------------------------------------------------------------

!!----------------------------------------------------------------------
!!   Default option :                                       Dummy module
!!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_diadct = .FALSE.    !: diamht flag
   PUBLIC 
!! $Id: diadct.F90 7646 2017-02-06 09:25:03Z timgraham $
CONTAINS

   SUBROUTINE dia_dct_init          ! Dummy routine
      WRITE(*,*) 'dia_dct_init: You should not have seen this print! error?', kt
   END SUBROUTINE dia_dct_init

   SUBROUTINE dia_dct( kt )         ! Dummy routine
      INTEGER, INTENT( in ) :: kt   ! ocean time-step index
      WRITE(*,*) 'dia_dct: You should not have seen this print! error?', kt
   END SUBROUTINE dia_dct


!!======================================================================
END MODULE diadct

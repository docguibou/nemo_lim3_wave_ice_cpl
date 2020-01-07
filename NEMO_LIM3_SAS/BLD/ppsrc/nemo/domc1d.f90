MODULE domc1d
!!======================================================================
!!                     ***  MODULE  domc1d  ***
!! Ocean Domain : 1D column position from lat/lon namelist specification
!!======================================================================
!! History :  3.5  !  2013-04  (D. Calvert)  Original code
!!----------------------------------------------------------------------

!!----------------------------------------------------------------------
!!   Default option                                  NO 1D Configuration
!!----------------------------------------------------------------------
CONTAINS  
   SUBROUTINE dom_c1d( plat, plon )     ! Empty routine
      REAL :: plat, plon
      WRITE(*,*) 'dom_c1d: You should not have seen this print! error?',plat,plon
   END SUBROUTINE dom_c1d


!!======================================================================
END MODULE domc1d

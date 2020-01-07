MODULE limthd_da
!!======================================================================
!!                       ***  MODULE limthd_da ***
!! LIM-3 sea-ice :  computation of lateral melting in the ice
!!======================================================================
!! History :   4.0   ! 2016-03 (C. Rousset) original code
!!---------------------------------------------------------------------

!!----------------------------------------------------------------------
!!   'key_lim3'                                      LIM-3 sea-ice model
!!----------------------------------------------------------------------
!!   lim_thd_da   : sea ice lateral melting
!!----------------------------------------------------------------------
   USE par_oce                ! ocean parameters
   USE phycst                 ! physical constants (ocean directory)
   USE sbc_oce, ONLY: sst_m   ! Surface boundary condition: ocean fields
   USE ice                    ! LIM variables
   USE sbc_ice , ONLY: ICEF_s ! LIM variables
   USE lib_mpp                ! MPP library
   USE wrk_nemo               ! work arrays
   USE lib_fortran            ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_thd_da        ! called by limthd module

!!----------------------------------------------------------------------
!! NEMO/LIM3 4.0 , UCL - NEMO Consortium (2011)
!! $Id: limthd_da.F90 5123 2015-03-04 16:06:03Z clem $
!! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lim_thd_da
!!-------------------------------------------------------------------
!!                ***  ROUTINE lim_thd_da  ***
!!
!! ** Purpose :   computes sea ice lateral melting
!!
!! ** Method  :   dA/dt = - P * W   [s-1]
!!                   W = melting velocity [m.s-1]
!!                   P = perimeter of ice-ocean lateral interface normalized by grid cell area [m.m-2]
!!
!!                   W = m1 * (Tw -Tf)**m2                    --- originally from Josberger 1979 ---
!!                      (Tw - Tf) = elevation of water temp above freezing
!!                      m1 and m2 = (1.6e-6 , 1.36) best fit from field experiment near the coast of Prince Patrick Island (Perovich 1983) => static ice
!!                      m1 and m2 = (3.0e-6 , 1.36) best fit from MIZEX 84 experiment (Maykut and Perovich 1987) => moving ice
!!
!!                   P = N * pi * D                           --- from Rothrock and Thorndike 1984 ---
!!                      D = mean floe caliper diameter
!!                      N = number of floes = ice area / floe area(average) = A / (Cs * D**2)
!!                         A = ice concentration
!!                         Cs = deviation from a square (square:Cs=1 ; circle:Cs=pi/4 ; floe:Cs=0.66)
!!
!!                   D = Dmin * ( Astar / (Astar-A) )**beta   --- from Lupkes et al., 2012 (eq. 26-27) ---
!!
!!                      Astar = 1 / ( 1 - (Dmin/Dmax)**(1/beta) )
!!                      Dmin = minimum floe diameter (recommended to be 8m +- 20%)
!!                      Dmax = maximum floe diameter (recommended to be 300m, but it does not impact melting much except for Dmax<100m)
!!                      beta = 1.0 +-20% (recommended value)
!!                           = 0.3 best fit for western Fram Strait and Antarctica
!!                           = 1.4 best fit for eastern Fram Strait
!!
!! ** Tunable parameters  :   We propose to tune the lateral melting via 2 parameters
!!                               Dmin [6-10m]   => 6  vs 8m = +40% melting at the peak (A~0.5)
!!                                                 10 vs 8m = -20% melting
!!                               beta [0.8-1.2] => decrease = more melt and melt peaks toward higher concentration
!!                                                                  (A~0.5 for beta=1 ; A~0.8 for beta=0.2)
!!                                                 0.3 = best fit for western Fram Strait and Antarctica
!!                                                 1.4 = best fit for eastern Fram Strait
!!
!! ** Note   :   Former and more simple formulations for floe diameters can be found in Mai (1995),
!!               Birnbaum and Lupkes (2002), Lupkes and Birnbaum (2005). They are reviewed in Lupkes et al 2012
!!               A simpler implementation for CICE can be found in Bitz et al (2001) and Tsamados et al (2015)
!!
!! ** References
!!    Bitz, C. M., Holland, M. M., Weaver, A. J., & Eby, M. (2001).
!!              Simulating the ice‐thickness distribution in a coupled climate model.
!!              Journal of Geophysical Research: Oceans, 106(C2), 2441-2463.
!!    Josberger, E. G. (1979).
!!              Laminar and turbulent boundary layers adjacent to melting vertical ice walls in salt water
!!              (No. SCIENTIFIC-16). WASHINGTON UNIV SEATTLE DEPT OF ATMOSPHERIC SCIENCES.
!!    Lüpkes, C., Gryanik, V. M., Hartmann, J., & Andreas, E. L. (2012).
!!              A parametrization, based on sea ice morphology, of the neutral atmospheric drag coefficients
!!              for weather prediction and climate models.
!!              Journal of Geophysical Research: Atmospheres, 117(D13).
!!    Maykut, G. A., & Perovich, D. K. (1987).
!!              The role of shortwave radiation in the summer decay of a sea ice cover.
!!              Journal of Geophysical Research: Oceans, 92(C7), 7032-7044.
!!    Perovich, D. K. (1983).
!!              On the summer decay of a sea ice cover. (Doctoral dissertation, University of Washington).
!!    Rothrock, D. A., & Thorndike, A. S. (1984).
!!              Measuring the sea ice floe size distribution.
!!              Journal of Geophysical Research: Oceans, 89(C4), 6477-6486.
!!    Tsamados, M., Feltham, D., Petty, A., Schroeder, D., & Flocco, D. (2015).
!!              Processes controlling surface, bottom and lateral melt of Arctic sea ice in a state of the art sea ice model.
!!              Phil. Trans. R. Soc. A, 373(2052), 20140167.
!!---------------------------------------------------------------------
      INTEGER             ::   ji, jj, jl, jf      ! dummy loop indices
      REAL(wp)            ::   zastar, zdfloe, zperi, zwlat, zda, zgamma, zfrag, zperi_fl, zlmelt, zamelt
      REAL(wp), PARAMETER ::   zdmax = 1000._wp
      REAL(wp), PARAMETER ::   zcs   = 0.66_wp
      REAL(wp), PARAMETER ::   zm1   = 3.e-6_wp
      REAL(wp), PARAMETER ::   zm2   = 1.36_wp
!
      REAL(wp), POINTER, DIMENSION(:) ::     zfl_dl, zdfl_dl
      REAL(wp), POINTER, DIMENSION(:,:) ::   zda_tot, zda_fl
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zafl_i_old
!!---------------------------------------------------------------------
      CALL wrk_alloc( jpi,jpj, zda_tot, zda_fl )
      CALL wrk_alloc( jpf+1,  zfl_dl)
      CALL wrk_alloc( jpf, zdfl_dl )
      CALL wrk_alloc( jpi,jpj,jpf, zafl_i_old )
!------------------------------------------------------------!
! --- Calculate reduction of total sea ice concentration --- !
!------------------------------------------------------------!
      zastar = 1._wp / ( 1._wp - (rn_dmin / zdmax)**(1._wp/rn_beta) )
      zfrag=0.9                       ! fragility (see Toyota et al. 2011), may be passed as a namelist parameter
      zgamma= 2 + log(zfrag)/log(2.)
      
! Be sure that there is initially no difference
!afl_i(:,:,:)  = afl_i(:,:,:) * SUM(a_i(:,:,:),dim=3)/MAX(SUM(afl_i(:,:,:),dim=3),epsi20)
            
      at_i(:,:) = SUM( a_i(:,:,:), dim=3 )

      zda_fl(:,:) = 0.
      zafl_i_old(:,:,:)=afl_i(:,:,:)

      DO jj = 1, jpj
         DO ji = 1, jpi     
! Melt speed rate [m/s]
           zwlat = zm1 * ( MAX( 0._wp, sst_m(ji,jj) - ( t_bo(ji,jj) - rt0 ) ) )**zm2       
           rswitch=MAX( 0._wp , SIGN( 1._wp,at_i(ji,jj))-epsi10 )
                 
           IF (nn_lamp ==1) THEN
              zwlat = 2._wp*zwlat  ! If Horvat/Roach only
              IF (ABS(SUM(zafl_i_old(ji,jj,:))-at_i(ji,jj) ).GT.epsi10 ) THEN
! WRITE(*,*) "WARNING : SUM FSD CONC different from SUM ITD TO BEGIN WITH"
! WRITE(*,*) "SUM FSD CONC :",SUM(zafl_i_old(ji,jj,:)) ,"  ,  SUM ITD CONC :",SUM(a_i(ji,jj,:))
! WRITE(*,*) "DIFF =", SUM(zafl_i_old(ji,jj,:))-SUM(a_i(ji,jj,:))
! WRITE(*,*) "Coordinates :",ji,jj
                 DO jf = 1, jpf 
                    zafl_i_old(ji,jj,jf)=zafl_i_old(ji,jj,jf)*at_i(ji,jj)*rswitch/MAX( SUM(zafl_i_old(ji,jj,:)), epsi20 )
                 END DO
              END IF
              zda_fl(ji,jj)=0.
              DO jf = 1, jpf 
                  
!! Mean floe caliper diameter [m]
!zdfloe = rswitch * di_mean(jf)
! Choice 1:
!! Mean perimeter of the floe = N*pi*D = (A/cs*D^2)*pi*D [m.m-2]
!zperi_fl = afl_i(ji,jj,jf) * rpi / ( zcs * zdfloe )
!! sea ice concentration decrease by floe size cat.
!zda_fl(ji,jj,jf)      =  - MIN( zwlat * zperi_fl *rdt_ice, afl_i(ji,jj,jf) )
!! Choice 2:
!! Mean perimeter of the floe = N*pi*D = (A/cs*D^2)*pi*D [m.m-2]
!zperi_fl=rpi / ( zcs * zdfloe )
!! sea ice concentration decrease by floe size cat.
!zda_fl(ji,jj,jf)   =  - MIN( afl_i(ji,jj,jf) * (1._wp-EXP(- rdt_ice * zwlat * zperi_fl) ) , afl_i(ji,jj,jf) )
                 
!! corresponding sea ice floe diameter decrease (It's wrong actually I think)
!zdl_fl(jf)          =    MIN(zwlat* rpi /( zcs*2._wp) * rdt_ice, di_max(jf)-di_max(jf-1))

! Roach et al. (2018)
                 zda_fl(ji,jj)   =zda_fl(ji,jj)  -  rdt_ice* 2._wp *zwlat *zafl_i_old(ji,jj,jf)/di_mean(jf)
! End Roach et al. (2018)

               ENDDO

               zda_fl(ji,jj)    = zda_fl(ji,jj)  - zafl_i_old(ji,jj,1)/(di_max(2)-di_max(1)) *rdt_ice*zwlat

! sea ice concentration decrease
               zda_tot(ji,jj) = MIN(0._wp,- MIN( -zda_fl(ji,jj), at_i(ji,jj) ))
           ELSE
! Mean floe caliper diameter [m]
               zdfloe = rn_dmin * ( zastar / ( zastar - at_i(ji,jj) ) )**rn_beta
! Mean perimeter of the floe = N*pi*D = (A/cs*D^2)*pi*D [m.m-2]
               zperi = at_i(ji,jj) * rpi / ( zcs * zdfloe )

! sea ice concentration decrease
               zda_tot(ji,jj) = - MIN( zwlat * zperi * rdt_ice, at_i(ji,jj) )
           ENDIF

         END DO
      END DO
     
!---------------------------------------------------------------------------------------------!
! --- Distribute reduction among ice categories and calculate associated ice-ocean fluxes --- !
!---------------------------------------------------------------------------------------------!
      DO jl = jpl, 1, -1
         DO jj = 1, jpj
            DO ji = 1, jpi
               
! decrease of concentration for the category jl
!    1st option: each category contributes to melting in proportion to its concentration
               rswitch = MAX( 0._wp , SIGN( 1._wp, at_i(ji,jj) - epsi10 ) )
               zda     = rswitch * zda_tot(ji,jj) * a_i(ji,jj,jl) / MAX( at_i(ji,jj), epsi10 )

!    2d option: melting of the upper cat first
!!zda = MAX( zda_tot(ji,jj), - a_i(ji,jj,jl) )
!!zda_tot(ji,jj) = zda_tot(ji,jj) + zda
               
! Contribution to salt flux
               sfx_lam(ji,jj) = sfx_lam(ji,jj) - rhoic *  ht_i(ji,jj,jl) * zda * sm_i(ji,jj,jl) * r1_rdtice
               
! Contribution to heat flux into the ocean [W.m-2], <0
               hfx_thd(ji,jj) = hfx_thd(ji,jj) + zda * r1_rdtice * ( ht_i(ji,jj,jl) * SUM( e_i(ji,jj,:,jl) ) * r1_nlay_i  &
                  &                                                + ht_s(ji,jj,jl) *      e_s(ji,jj,1,jl)   * r1_nlay_s )
               
! Contribution to mass flux
               wfx_lam(ji,jj) =  wfx_lam(ji,jj) - zda * r1_rdtice * ( rhoic * ht_i(ji,jj,jl) + rhosn * ht_s(ji,jj,jl) )
! IF (dmaxfl_i(ji,jj).GT.0) WRITE(*,*) 'wfx',wfx_lam(ji,jj),'dmax',dmaxfl_i(ji,jj) ,'lat melt', 'i',ji,'j',jj,'zda',zda,r1_rdtice,ht_i(ji,jj,jl)
! new concentration
               a_i(ji,jj,jl) = a_i(ji,jj,jl) + zda
            END DO
         END DO
      END DO
!---------------------------------------------------------------------------------------------!
! --- Distribute reduction among ice floes categories                                     --- !
!---------------------------------------------------------------------------------------------!
      IF (nn_lamp ==1) THEN ! SHOULD BE OK FOR GENERAL CASE NORMALLY BUT...
! Computation needed for the delta term (only fsd in this case)
        
         DO jj = 1, jpj
            DO ji = 1, jpi
               zwlat =  zm1 * ( MAX( 0._wp, sst_m(ji,jj) - ( t_bo(ji,jj) - rt0 ) ) )**zm2       
               zwlat = 2*zwlat
               zfl_dl(:)=0.

               DO jf = 2, jpf
                  zfl_dl(jf)=zwlat*zafl_i_old(ji,jj,jf)/(di_max(jf) -di_max(jf-1) )
               END DO

               zdfl_dl(:)=0.
               DO jf = 1, jpf
                  zdfl_dl(jf) = zfl_dl(jf+1)-zfl_dl(jf)
               END DO

               IF ( ABS(SUM(zdfl_dl(:))).GT.epsi10) THEN
                  WRITE(*,*) "SUM delta not equal to 0", SUM(zdfl_dl(:))
                  WRITE(*,*) "zwlat :"  ,zwlat
                  WRITE(*,*) "zfl_dl :" ,zfl_dl(:)
                  WRITE(*,*) "zdfl_dl :",zdfl_dl(:)
               ENDIF

               DO jf = 1, jpf
! decrease of concentration for the category jf
                  rswitch = MAX( 0._wp , SIGN( 1._wp, at_i(ji,jj) - epsi10 ) )
                  afl_i(ji,jj,jf) =   rswitch * ( zafl_i_old(ji,jj,jf) - rdt_ice * (&
                                      - zdfl_dl(jf)                            &
                                      + zwlat*zafl_i_old(ji,jj,jf)* 2._wp      &
                                      /di_mean(jf)                         )) 
               END DO
               afl_i(ji,jj,1) = afl_i(ji,jj,1) - zafl_i_old(ji,jj,1)/ (di_max(2)-di_max(1) )*zwlat*rdt_ice
            END DO
         END DO 

      END IF    


! total concentration
! Raw conservation test
      at_i(:,:) = SUM( a_i(:,:,:), dim=3 )
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF (at_i(ji,jj)>0) THEN
               IF (at_i(ji,jj)-SUM(a_i(ji,jj,:))>epsi06) WRITE(*,*) "NOT CONSERVED LAT. MELT THICK"
            ENDIF
         ENDDO
      ENDDO
      IF (nn_lamp ==1) THEN ! SHOULD BE OK FOR GENERAL CASE NORMALLY BUT...
        DO jj = 1, jpj
          DO ji = 1, jpi
            IF (at_i(ji,jj)>epsi06) THEN
               IF ( (ABS(at_i(ji,jj)-SUM( afl_i(ji,jj,:)) )>epsi06).AND.(ABS(zda_tot(ji,jj)-zda_fl(ji,jj))>epsi10) )  THEN
                  WRITE(*,*)" NOT CONSERVED FSD LAT. MELT -------------------------------------"
                  WRITE(*,*) "diff conc_tot - sum_fsd_conc =",at_i(ji,jj)-SUM( afl_i(ji,jj,:)),"\n"
                  WRITE(*,*) "conc tot =",at_i(ji,jj)," ; sum_concfsd=",SUM( afl_i(ji,jj,:))
                  WRITE(*,*) "old_conc_tot_fsd :", SUM(zafl_i_old(ji,jj,:))
                  WRITE(*,*) "Coordinates :",ji,jj
                  WRITE(*,*) "old fsd :", zafl_i_old(ji,jj,:)
                  WRITE(*,*) "zdafl :", zda_fl(ji,jj)
                  zwlat =  2*zm1 * ( MAX( 0._wp, sst_m(ji,jj) - ( t_bo(ji,jj) - rt0 ) ) )**zm2       
                  WRITE(*,*) "zdafl in eq redis:", - rdt_ice*zwlat* 2._wp* SUM( zafl_i_old(ji,jj,1:jpf)/di_mean(1:jpf) )
                  WRITE(*,*) "lower cat area loss :", - zafl_i_old(ji,jj,1)/ (di_max(2)-di_max(1) )*zwlat*rdt_ice
                  WRITE(*,*) "zdafl_final : ",  - rdt_ice*zwlat*  ( zafl_i_old(ji,jj,1)/(di_max(2)-di_max(1)) + 2._wp*SUM( zafl_i_old(ji,jj,1:jpf)/di_mean(1:jpf) ) ) 
                  WRITE(*,*) "zda_tot =", zda_tot(ji,jj)
               END IF
            END IF
          END DO
        END DO 
      END IF     
! --- ensure that ht_i = 0 where a_i = 0 ---
      WHERE( a_i == 0._wp ) ht_i = 0._wp
!
      CALL wrk_dealloc( jpi,jpj, zda_tot,zda_fl)
      CALL wrk_dealloc( jpf+1,  zfl_dl )
      CALL wrk_dealloc( jpf, zdfl_dl )
      CALL wrk_dealloc( jpi,jpj,jpf, zafl_i_old )
!
   END SUBROUTINE lim_thd_da
   

!!======================================================================
END MODULE limthd_da

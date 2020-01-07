MODULE limistate
   !!======================================================================
   !!                     ***  MODULE  limistate  ***
   !!              Initialisation of diagnostics ice variables
   !!======================================================================
   !! History :  2.0  ! 2004-01 (C. Ethe, G. Madec)  Original code
   !!            4.0  ! 2011-02 (G. Madec) dynamical allocation
   !!             -   ! 2014    (C. Rousset) add N/S initializations
   !!----------------------------------------------------------------------
#if defined key_lim3
   !!----------------------------------------------------------------------
   !!   'key_lim3' :                                    LIM3 sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_istate      :  Initialisation of diagnostics ice variables
   !!   lim_istate_init :  initialization of ice state and namelist read
   !!----------------------------------------------------------------------
   USE phycst           ! physical constant
   USE oce              ! dynamics and tracers variables
   USE dom_oce          ! ocean domain
   USE sbc_oce          ! Surface boundary condition: ocean fields
   USE sbc_ice          ! Surface boundary condition: ice fields
   USE eosbn2           ! equation of state
   USE ice              ! sea-ice variables
   USE par_oce          ! ocean parameters
   USE dom_ice          ! sea-ice domain
   USE limvar           ! lim_var_salprof
   USE in_out_manager   ! I/O manager
   USE lib_mpp          ! MPP library
   USE lib_fortran      ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  
   USE wrk_nemo         ! work arrays
!{ CREG modifications 
   USE fldread          ! read input fields
   USE iom
   USE prtctl           ! Print control
!}

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_istate      ! routine called by lim_init.F90

   !                          !!** init namelist (namiceini) **
   REAL(wp) ::   rn_thres_sst   ! threshold water temperature for initial sea ice
   REAL(wp) ::   rn_hts_ini_n   ! initial snow thickness in the north
   REAL(wp) ::   rn_hts_ini_s   ! initial snow thickness in the south
   REAL(wp) ::   rn_hti_ini_n   ! initial ice thickness in the north
   REAL(wp) ::   rn_hti_ini_s   ! initial ice thickness in the south
   REAL(wp) ::   rn_ati_ini_n   ! initial leads area in the north
   REAL(wp) ::   rn_ati_ini_s   ! initial leads area in the south
   REAL(wp) ::   rn_smi_ini_n   ! initial salinity 
   REAL(wp) ::   rn_smi_ini_s   ! initial salinity
   REAL(wp) ::   rn_tmi_ini_n   ! initial temperature
   REAL(wp) ::   rn_tmi_ini_s   ! initial temperature

!{ CREG modifications 
   INTEGER , PARAMETER ::   jpfldi    = 7           ! maximum number of files to read
   INTEGER , PARAMETER ::   jp_hicif = 1           ! index of thick (m)    at T-point
   INTEGER , PARAMETER ::   jp_hsnif = 2           ! index of thick (m)    at T-point
   INTEGER , PARAMETER ::   jp_frld  = 3           ! index of ice fraction (%) at T-point
   INTEGER , PARAMETER ::   jp_sist  = 4           ! index of ice surface temp (K)    at T-point
   INTEGER , PARAMETER ::   jp_tbif1 = 5           ! index of ice temp lev1 (K) at T-point
   INTEGER , PARAMETER ::   jp_tbif2 = 6           ! index of ice temp lev2 (K) at T-point
   INTEGER , PARAMETER ::   jp_tbif3 = 7           ! index of ice temp lev3 (K) at T-point
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   si    ! structure of input fields (file informations, fields read)

   REAL(wp),DIMENSION(:,:)  ,ALLOCATABLE :: hicif_ini,hsnif_ini,frld_ini,sist_ini
   REAL(wp),DIMENSION(:,:,:),ALLOCATABLE :: tbif_ini
   REAL(wp), POINTER, DIMENSION(:,:)   :: zswitch    ! ice indicator
!}

   LOGICAL  ::  ln_iceini    ! initialization or not
!{ CREG modifications 
   LOGICAL  ::  ln_limini_file   ! Ice initialization state from 2D netcdf file
!}
   !!----------------------------------------------------------------------
   !!   LIM 3.0,  UCL-LOCEAN-IPSL (2008)
   !! $Id: limistate.F90 6469 2016-04-13 15:15:13Z clem $
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lim_istate
      !!-------------------------------------------------------------------
      !!                    ***  ROUTINE lim_istate  ***
      !!
      !! ** Purpose :   defined the sea-ice initial state
      !!
      !! ** Method  :   
      !!                This routine will put some ice where ocean
      !!                is at the freezing point, then fill in ice 
      !!                state variables using prescribed initial 
      !!                values in the namelist            
      !!
      !! ** Steps   :   
      !!                1) Read namelist
      !!                2) Basal temperature; ice and hemisphere masks
      !!                3) Fill in the ice thickness distribution using gaussian
      !!                4) Fill in space-dependent arrays for state variables
      !!                5) Diagnostic arrays
      !!                6) Lateral boundary conditions
      !!
      !! ** Notes   : o_i, t_su, t_s, t_i, s_i must be filled everywhere, even
      !!              where there is no ice (clem: I do not know why, is it mandatory?) 
      !!
      !! History :
      !!   2.0  !  01-04  (C. Ethe, G. Madec)  Original code
      !!   3.0  !  2007   (M. Vancoppenolle)   Rewrite for ice cats
      !!   4.0  !  09-11  (M. Vancoppenolle)   Enhanced version for ice cats
      !!--------------------------------------------------------------------

      !! * Local variables
      INTEGER    :: ji, jj, jk, jl, ja, jf         ! dummy loop indices
      REAL(wp)   :: ztmelts, zdh
      INTEGER    :: i_hemis, i_fill, jl0  
      REAL(wp)   :: ztest_1, ztest_2, ztest_3, ztest_4, ztests, zsigma, zarg, zA, zV, zA_cons, zV_cons, zconv
      REAL(wp), POINTER, DIMENSION(:)     :: zht_i_ini, zat_i_ini, zvt_i_ini, zht_s_ini, zsm_i_ini, ztm_i_ini
      REAL(wp), POINTER, DIMENSION(:,:) :: zafl_i_ini        !data by cattegories to fill
      REAL(wp), POINTER, DIMENSION(:,:)   :: zh_i_ini, za_i_ini, zv_i_ini
      !REAL(wp), POINTER, DIMENSION(:,:)   :: zswitch    ! ice indicator
      INTEGER,  POINTER, DIMENSION(:,:)   :: zhemis   ! hemispheric index
      !--------------------------------------------------------------------

      CALL wrk_alloc( jpi, jpj, zswitch )
      CALL wrk_alloc( jpi, jpj, zhemis )
      CALL wrk_alloc( jpf,   2, zafl_i_ini )
      CALL wrk_alloc( jpl,   2, zh_i_ini,  za_i_ini,  zv_i_ini )
      CALL wrk_alloc(   2,      zht_i_ini, zat_i_ini, zvt_i_ini, zht_s_ini, zsm_i_ini, ztm_i_ini )

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'lim_istate : Ice initialization '
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '

      !--------------------------------------------------------------------
      ! 1) Read namelist
      !--------------------------------------------------------------------

      CALL lim_istate_init     !  reading the initials parameters of the ice

      ! surface temperature
      DO jl = 1, jpl ! loop over categories
         t_su  (:,:,jl) = rt0 * tmask(:,:,1)
         tn_ice(:,:,jl) = rt0 * tmask(:,:,1)
      END DO

      ! basal temperature (considered at freezing point)
      CALL eos_fzp( sss_m(:,:), t_bo(:,:) )
      t_bo(:,:) = ( t_bo(:,:) + rt0 ) * tmask(:,:,1) 

      IF( ln_iceini ) THEN

         !--------------------------------------------------------------------
         ! 2) Basal temperature, ice mask and hemispheric index
         !--------------------------------------------------------------------

         DO jj = 1, jpj                                       ! ice if sst <= t-freez + ttest
            DO ji = 1, jpi
               IF( ( sst_m(ji,jj)  - ( t_bo(ji,jj) - rt0 ) ) * tmask(ji,jj,1) >= rn_thres_sst ) THEN 
                  zswitch(ji,jj) = 0._wp * tmask(ji,jj,1)    ! no ice
               ELSE                                                                                   
                  zswitch(ji,jj) = 1._wp * tmask(ji,jj,1)    !    ice
               ENDIF
            END DO
         END DO


         ! Hemispheric index
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( fcor(ji,jj) >= 0._wp ) THEN    
                  zhemis(ji,jj) = 1 ! Northern hemisphere
               ELSE
                  zhemis(ji,jj) = 2 ! Southern hemisphere
               ENDIF
            END DO
         END DO

         !--------------------------------------------------------------------
         ! 3) Initialization of sea ice state variables
         !--------------------------------------------------------------------
!{ CREG modifications 
         IF( ln_limini_file )THEN

            CALL limini_file
            IF(ln_ctl) THEN         ! print mean trends (used for debugging)
               DO jl = 1, jpl
                  CALL prt_ctl_info(' ')
                  CALL prt_ctl_info(' - Category : ', ivar1=jl)
                  CALL prt_ctl(tab2d_1=a_i   (:,:,jl)   , clinfo1= ' lim_istate  : a_i      : ')
                  CALL prt_ctl(tab2d_1=o_i   (:,:,jl)   , clinfo1= ' lim_istate  : o_i      : ')
                  CALL prt_ctl(tab2d_1=ht_i  (:,:,jl)   , clinfo1= ' lim_istate  : ht_i     : ')
                  CALL prt_ctl(tab2d_1=ht_s  (:,:,jl)   , clinfo1= ' lim_istate  : ht_s     : ')
                  CALL prt_ctl(tab2d_1=v_i   (:,:,jl)   , clinfo1= ' lim_istate  : v_i      : ')
                  CALL prt_ctl(tab2d_1=v_s   (:,:,jl)   , clinfo1= ' lim_istate  : v_s      : ')
                  CALL prt_ctl(tab2d_1=t_su  (:,:,jl)   , clinfo1= ' lim_istate  : t_su     : ')
                  CALL prt_ctl(tab2d_1=sm_i  (:,:,jl)   , clinfo1= ' lim_istate  : sm_i     : ')
                  CALL prt_ctl(tab2d_1=smv_i (:,:,jl)   , clinfo1= ' lim_istate  : smv_i    : ')
                  CALL prt_ctl(tab2d_1=oa_i  (:,:,jl)   , clinfo1= ' lim_istate  : oa_i     : ')
                  DO ja = 1, nlay_i
                     CALL prt_ctl_info(' ')
                     CALL prt_ctl_info(' - Layer : ', ivar1=ja)
                     CALL prt_ctl_info('   ~~~~~~~')
                     CALL prt_ctl(tab2d_1=t_i(:,:,ja,jl) , clinfo1= ' lim_istate  : t_i      : ')
                     CALL prt_ctl(tab2d_1=e_i(:,:,ja,jl) , clinfo1= ' lim_istate  : e_i      : ')
                     CALL prt_ctl(tab2d_1=s_i(:,:,ja,jl) , clinfo1= ' lim_istate  : s_i      : ')
                     CALL prt_ctl(tab2d_1=t_s(:,:,ja,jl) , clinfo1= ' lim_istate  : t_s      : ')
                     CALL prt_ctl(tab2d_1=e_s(:,:,ja,jl) , clinfo1= ' lim_istate  : e_s      : ')
                  END DO
               END DO
               CALL prt_ctl(tab2d_1=zswitch              , clinfo1=' zswitch      - : ', mask1=tmask, ovlap=1 )
               CALL prt_ctl(tab2d_1=t_bo                 , clinfo1=' t_bo         - : ', mask1=tmask, ovlap=1 )
            ENDIF

         ELSE
!}

            !-----------------------------
            ! 3.1) Hemisphere-dependent arrays
            !-----------------------------
            ! assign initial thickness, concentration, snow depth and salinity to an hemisphere-dependent array
            zht_i_ini(1) = rn_hti_ini_n ; zht_i_ini(2) = rn_hti_ini_s  ! ice thickness
            zht_s_ini(1) = rn_hts_ini_n ; zht_s_ini(2) = rn_hts_ini_s  ! snow depth
            zat_i_ini(1) = rn_ati_ini_n ; zat_i_ini(2) = rn_ati_ini_s  ! ice concentration
            zsm_i_ini(1) = rn_smi_ini_n ; zsm_i_ini(2) = rn_smi_ini_s  ! bulk ice salinity
            ztm_i_ini(1) = rn_tmi_ini_n ; ztm_i_ini(2) = rn_tmi_ini_s  ! temperature (ice and snow)

            zvt_i_ini(:) = zht_i_ini(:) * zat_i_ini(:)   ! ice volume

            !---------------------------------------------------------------------
            ! 3.2) Distribute ice concentration and thickness into the categories
            !---------------------------------------------------------------------
            ! a gaussian distribution for ice concentration is used
            ! then we check whether the distribution fullfills
            ! volume and area conservation, positivity and ice categories bounds
            DO i_hemis = 1, 2

            ztest_1 = 0 ; ztest_2 = 0 ; ztest_3 = 0 ; ztest_4 = 0

            ! note for the great nemo engineers: 
            ! only very few of the WRITE statements are necessary for the reference version
            ! they were one day useful, but now i personally doubt of their
            ! potential for bringing anything useful

            DO i_fill = jpl, 1, -1
               IF ( ( ztest_1 + ztest_2 + ztest_3 + ztest_4 ) .NE. 4 ) THEN
                  !----------------------------
                  ! fill the i_fill categories
                  !----------------------------
                  ! *** 1 category to fill
                  IF ( i_fill .EQ. 1 ) THEN
                     zh_i_ini(1,i_hemis)       = zht_i_ini(i_hemis)
                     za_i_ini(1,i_hemis)       = zat_i_ini(i_hemis)
                     zafl_i_ini(jpf,i_hemis)   = zat_i_ini(i_hemis) ! Only the last category is filled since floes are supposed to be unbroken. (may be changed ?)
                     zh_i_ini(2:jpl,i_hemis)   = 0._wp
                     za_i_ini(2:jpl,i_hemis)   = 0._wp
                     zafl_i_ini(1:jpf-1,i_hemis)=0._wp 
                  ELSE

                     ! *** >1 categores to fill
                     !--- Ice thicknesses in the i_fill - 1 first categories
                     DO jl = 1, i_fill - 1
                        zh_i_ini(jl,i_hemis) = hi_mean(jl)
                     END DO
                     
                     !--- jl0: most likely index where cc will be maximum
                     DO jl = 1, jpl
                        IF ( ( zht_i_ini(i_hemis) >  hi_max(jl-1) ) .AND. &
                           & ( zht_i_ini(i_hemis) <= hi_max(jl)   ) ) THEN
                           jl0 = jl
                        ENDIF
                     END DO
                     jl0 = MIN(jl0, i_fill)
                     
                     !--- Concentrations
                     za_i_ini(jl0,i_hemis)      = zat_i_ini(i_hemis) / SQRT(REAL(jpl))
                     DO jl = 1, i_fill - 1
                        IF ( jl .NE. jl0 ) THEN
                           zsigma               = 0.5 * zht_i_ini(i_hemis)
                           zarg                 = ( zh_i_ini(jl,i_hemis) - zht_i_ini(i_hemis) ) / zsigma
                           za_i_ini(jl,i_hemis) = za_i_ini(jl0,i_hemis) * EXP(-zarg**2)
                        ENDIF
                     END DO
                     
                     zA = 0. ! sum of the areas in the jpl categories 
                     DO jl = 1, i_fill - 1
                       zA = zA + za_i_ini(jl,i_hemis)
                     END DO
                     za_i_ini(i_fill,i_hemis)   = zat_i_ini(i_hemis) - zA ! ice conc in the last category
                     IF ( i_fill .LT. jpl ) za_i_ini(i_fill+1:jpl, i_hemis) = 0._wp
               
                     !--- Ice thickness in the last category
                     zV = 0. ! sum of the volumes of the N-1 categories
                     DO jl = 1, i_fill - 1
                        zV = zV + za_i_ini(jl,i_hemis)*zh_i_ini(jl,i_hemis)
                     END DO
                     zh_i_ini(i_fill,i_hemis) = ( zvt_i_ini(i_hemis) - zV ) / za_i_ini(i_fill,i_hemis) 
                     IF ( i_fill .LT. jpl ) zh_i_ini(i_fill+1:jpl, i_hemis) = 0._wp

                     !--- volumes
                     zv_i_ini(:,i_hemis) = za_i_ini(:,i_hemis) * zh_i_ini(:,i_hemis)
                     IF ( i_fill .LT. jpl ) zv_i_ini(i_fill+1:jpl, i_hemis) = 0._wp

                  ENDIF ! i_fill

                  !---------------------
                  ! Compatibility tests
                  !---------------------
                  ! Test 1: area conservation
                  zA_cons = SUM(za_i_ini(:,i_hemis)) ; zconv = ABS(zat_i_ini(i_hemis) - zA_cons )
                  IF ( zconv .LT. 1.0e-6 ) THEN
                     ztest_1 = 1
                  ELSE 
                    ! this write is useful
                    IF(lwp)  WRITE(numout,*) ' * TEST1 AREA NOT CONSERVED *** zA_cons = ', zA_cons,' zat_i_ini = ',zat_i_ini(i_hemis) 
                     ztest_1 = 0
                  ENDIF

                  ! Test 2: volume conservation
                  zV_cons = SUM(zv_i_ini(:,i_hemis))
                  zconv = ABS(zvt_i_ini(i_hemis) - zV_cons)

                  IF ( zconv .LT. 1.0e-6 ) THEN
                     ztest_2 = 1
                  ELSE
                    ! this write is useful
                    IF(lwp)  WRITE(numout,*) ' * TEST2 VOLUME NOT CONSERVED *** zV_cons = ', zV_cons, &
                                  ' zvt_i_ini = ', zvt_i_ini(i_hemis)
                     ztest_2 = 0
                  ENDIF

                  ! Test 3: thickness of the last category is in-bounds ?
                  IF ( zh_i_ini(i_fill, i_hemis) > hi_max(i_fill-1) ) THEN
                     ztest_3 = 1
                  ELSE
                     ! this write is useful
                     IF(lwp) WRITE(numout,*) ' * TEST 3 THICKNESS OF THE LAST CATEGORY OUT OF BOUNDS *** zh_i_ini(i_fill,i_hemis) = ', &
                     zh_i_ini(i_fill,i_hemis), ' hi_max(jpl-1) = ', hi_max(i_fill-1)
                     ztest_3 = 0
                  ENDIF

                  ! Test 4: positivity of ice concentrations
                  ztest_4 = 1
                  DO jl = 1, jpl
                     IF ( za_i_ini(jl,i_hemis) .LT. 0._wp ) THEN 
                        ! this write is useful
                        IF(lwp) WRITE(numout,*) ' * TEST 4 POSITIVITY NOT OK FOR CAT ', jl, ' WITH A = ', za_i_ini(jl,i_hemis)
                        ztest_4 = 0
                     ENDIF
                  END DO

               ENDIF ! ztest_1 + ztest_2 + ztest_3 + ztest_4
 
               ztests = ztest_1 + ztest_2 + ztest_3 + ztest_4

            END DO ! i_fill

            IF(lwp) THEN 
               WRITE(numout,*) ' ztests : ', ztests
               IF ( ztests .NE. 4 ) THEN
                  WRITE(numout,*)
                  WRITE(numout,*) ' !!!! ALERT                  !!! '
                  WRITE(numout,*) ' !!!! Something is wrong in the LIM3 initialization procedure '
                  WRITE(numout,*)
                  WRITE(numout,*) ' *** ztests is not equal to 4 '
                  WRITE(numout,*) ' *** ztest_i (i=1,4) = ', ztest_1, ztest_2, ztest_3, ztest_4
                  WRITE(numout,*) ' zat_i_ini : ', zat_i_ini(i_hemis)
                  WRITE(numout,*) ' zht_i_ini : ', zht_i_ini(i_hemis)
               ENDIF ! ztests .NE. 4
            ENDIF
            
            END DO ! i_hemis

            !---------------------------------------------------------------------
            ! 3.3) Space-dependent arrays for ice state variables
            !---------------------------------------------------------------------

            ! Ice concentration, thickness and volume, ice salinity, ice age, surface temperature
            DO jl = 1, jpl ! loop over categories
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     a_i(ji,jj,jl)   = zswitch(ji,jj) * za_i_ini (jl,zhemis(ji,jj))  ! concentration
                     ht_i(ji,jj,jl)  = zswitch(ji,jj) * zh_i_ini(jl,zhemis(ji,jj))   ! ice thickness
                     ht_s(ji,jj,jl)  = ht_i(ji,jj,jl) * ( zht_s_ini( zhemis(ji,jj) ) / zht_i_ini( zhemis(ji,jj) ) )  ! snow depth
                     sm_i(ji,jj,jl)  = zswitch(ji,jj) * zsm_i_ini(zhemis(ji,jj))     ! salinity
                     o_i(ji,jj,jl)   = zswitch(ji,jj) * 1._wp                        ! age (1 day)
                     t_su(ji,jj,jl)  = zswitch(ji,jj) * ztm_i_ini(zhemis(ji,jj)) + ( 1._wp - zswitch(ji,jj) ) * rt0 ! surf temp

                     ! This case below should not be used if (ht_s/ht_i) is ok in namelist
                     ! In case snow load is in excess that would lead to transformation from snow to ice
                     ! Then, transfer the snow excess into the ice (different from limthd_dh)
                     zdh = MAX( 0._wp, ( rhosn * ht_s(ji,jj,jl) + ( rhoic - rau0 ) * ht_i(ji,jj,jl) ) * r1_rau0 ) 
                     ! recompute ht_i, ht_s avoiding out of bounds values
                     ht_i(ji,jj,jl) = MIN( hi_max(jl), ht_i(ji,jj,jl) + zdh )
                     ht_s(ji,jj,jl) = MAX( 0._wp, ht_s(ji,jj,jl) - zdh * rhoic * r1_rhosn )

                     ! ice volume, salt content, age content
                     v_i(ji,jj,jl)   = ht_i(ji,jj,jl) * a_i(ji,jj,jl)              ! ice volume
                     v_s(ji,jj,jl)   = ht_s(ji,jj,jl) * a_i(ji,jj,jl)              ! snow volume
                     smv_i(ji,jj,jl) = MIN( sm_i(ji,jj,jl) , sss_m(ji,jj) ) * v_i(ji,jj,jl) ! salt content
                     oa_i(ji,jj,jl)  = o_i(ji,jj,jl) * a_i(ji,jj,jl)               ! age content
                  END DO
               END DO
            END DO
         DO jf = 1, jpf ! loop over ice floes categories (GB)
            DO jj = 1, jpj
               DO ji = 1, jpi
                 afl_i(ji,jj,jf)   = zswitch(ji,jj) * zafl_i_ini(jf,zhemis(ji,jj))
               END DO
            END DO
         END DO

            ! for constant salinity in time
            IF( nn_icesal == 1 .OR. nn_icesal == 3 )  THEN
               CALL lim_var_salprof
               smv_i = sm_i * v_i
            ENDIF

            ! Snow temperature and heat content
            DO jk = 1, nlay_s
               DO jl = 1, jpl ! loop over categories
                  DO jj = 1, jpj
                     DO ji = 1, jpi
                         t_s(ji,jj,jk,jl) = zswitch(ji,jj) * ztm_i_ini(zhemis(ji,jj)) + ( 1._wp - zswitch(ji,jj) ) * rt0
                         ! Snow energy of melting
                         e_s(ji,jj,jk,jl) = zswitch(ji,jj) * rhosn * ( cpic * ( rt0 - t_s(ji,jj,jk,jl) ) + lfus )

                         ! Mutliply by volume, and divide by number of layers to get heat content in J/m2
                         e_s(ji,jj,jk,jl) = e_s(ji,jj,jk,jl) * v_s(ji,jj,jl) * r1_nlay_s
                     END DO
                  END DO
               END DO
            END DO

            ! Ice salinity, temperature and heat content
            DO jk = 1, nlay_i
               DO jl = 1, jpl ! loop over categories
                  DO jj = 1, jpj
                     DO ji = 1, jpi
                         t_i(ji,jj,jk,jl) = zswitch(ji,jj) * ztm_i_ini(zhemis(ji,jj)) + ( 1._wp - zswitch(ji,jj) ) * rt0 
                         s_i(ji,jj,jk,jl) = zswitch(ji,jj) * zsm_i_ini(zhemis(ji,jj)) !+ ( 1._wp - zswitch(ji,jj) ) * rn_simin
                         ztmelts          = - tmut * s_i(ji,jj,jk,jl) + rt0 !Melting temperature in K

                         ! heat content per unit volume
                         e_i(ji,jj,jk,jl) = zswitch(ji,jj) * rhoic * (   cpic    * ( ztmelts - t_i(ji,jj,jk,jl) ) &
                            +   lfus    * ( 1._wp - (ztmelts-rt0) / MIN((t_i(ji,jj,jk,jl)-rt0),-epsi20) ) &
                            -   rcp     * ( ztmelts - rt0 ) )

                         ! Mutliply by ice volume, and divide by number of layers to get heat content in J/m2
                         e_i(ji,jj,jk,jl) = e_i(ji,jj,jk,jl) * v_i(ji,jj,jl) * r1_nlay_i
                     END DO
                  END DO
               END DO
            END DO

            tn_ice (:,:,:) = t_su (:,:,:)

!{ CREG modifications 
            ENDIF !ln_limini_file
!}
      ELSE 
         ! if ln_iceini=false
         a_i  (:,:,:) = 0._wp
         v_i  (:,:,:) = 0._wp
         v_s  (:,:,:) = 0._wp
         smv_i(:,:,:) = 0._wp
         oa_i (:,:,:) = 0._wp
         ht_i (:,:,:) = 0._wp
         ht_s (:,:,:) = 0._wp
         sm_i (:,:,:) = 0._wp
         o_i  (:,:,:) = 0._wp
         afl_i(:,:,:) = 0._wp

         e_i(:,:,:,:) = 0._wp
         e_s(:,:,:,:) = 0._wp

         DO jl = 1, jpl
            DO jk = 1, nlay_i
               t_i(:,:,jk,jl) = rt0 * tmask(:,:,1)
            END DO
            DO jk = 1, nlay_s
               t_s(:,:,jk,jl) = rt0 * tmask(:,:,1)
            END DO
         END DO
      
      ENDIF ! ln_iceini
      
      at_i (:,:) = 0.0_wp
      DO jl = 1, jpl
         at_i (:,:) = at_i (:,:) + a_i (:,:,jl)
      END DO
            DO jj = 1, jpj
               DO ji = 1, jpi
               IF (ABS(SUM(afl_i(ji,jj,:))-at_i(ji,jj))>1e-6) WRITE(*,*) "LIMISTATE CONS. PAS BON", ji, jj, SUM(afl_i(ji,jj,:)), at_i(ji,jj)
               END DO
            END DO
      !
      !--------------------------------------------------------------------
      ! 4) Global ice variables for output diagnostics                    | 
      !--------------------------------------------------------------------
      u_ice (:,:)     = 0._wp
      v_ice (:,:)     = 0._wp
      stress1_i(:,:)  = 0._wp
      stress2_i(:,:)  = 0._wp
      stress12_i(:,:) = 0._wp

      !--------------------------------------------------------------------
      ! 5) Moments for advection
      !--------------------------------------------------------------------

      sxopw (:,:) = 0._wp 
      syopw (:,:) = 0._wp 
      sxxopw(:,:) = 0._wp 
      syyopw(:,:) = 0._wp 
      sxyopw(:,:) = 0._wp

      sxice (:,:,:)  = 0._wp   ;   sxsn (:,:,:)  = 0._wp   ;   sxa  (:,:,:)  = 0._wp
      syice (:,:,:)  = 0._wp   ;   sysn (:,:,:)  = 0._wp   ;   sya  (:,:,:)  = 0._wp
      sxxice(:,:,:)  = 0._wp   ;   sxxsn(:,:,:)  = 0._wp   ;   sxxa (:,:,:)  = 0._wp
      syyice(:,:,:)  = 0._wp   ;   syysn(:,:,:)  = 0._wp   ;   syya (:,:,:)  = 0._wp
      sxyice(:,:,:)  = 0._wp   ;   sxysn(:,:,:)  = 0._wp   ;   sxya (:,:,:)  = 0._wp

      sxc0  (:,:,:)  = 0._wp   ;   sxe  (:,:,:,:)= 0._wp   
      syc0  (:,:,:)  = 0._wp   ;   sye  (:,:,:,:)= 0._wp   
      sxxc0 (:,:,:)  = 0._wp   ;   sxxe (:,:,:,:)= 0._wp   
      syyc0 (:,:,:)  = 0._wp   ;   syye (:,:,:,:)= 0._wp   
      sxyc0 (:,:,:)  = 0._wp   ;   sxye (:,:,:,:)= 0._wp   

      sxsal  (:,:,:)  = 0._wp
      sysal  (:,:,:)  = 0._wp
      sxxsal (:,:,:)  = 0._wp
      syysal (:,:,:)  = 0._wp
      sxysal (:,:,:)  = 0._wp

      sxage  (:,:,:)  = 0._wp
      syage  (:,:,:)  = 0._wp
      sxxage (:,:,:)  = 0._wp
      syyage (:,:,:)  = 0._wp
      sxyage (:,:,:)  = 0._wp

      sxfl  (:,:,:)  = 0._wp
      syfl  (:,:,:)  = 0._wp
      sxxfl (:,:,:)  = 0._wp
      syyfl (:,:,:)  = 0._wp
      sxyfl (:,:,:)  = 0._wp


      CALL wrk_dealloc( jpi, jpj, zswitch )
      CALL wrk_dealloc( jpi, jpj, zhemis )
      CALL wrk_dealloc( jpl,   2, zh_i_ini,  za_i_ini,  zv_i_ini )
      CALL wrk_dealloc( jpf,   2,   zafl_i_ini)
      CALL wrk_dealloc(   2,      zht_i_ini, zat_i_ini, zvt_i_ini, zht_s_ini, zsm_i_ini, ztm_i_ini )

   END SUBROUTINE lim_istate

   SUBROUTINE lim_istate_init
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE lim_istate_init  ***
      !!        
      !! ** Purpose : Definition of initial state of the ice 
      !!
      !! ** Method : Read the namiceini namelist and check the parameter 
      !!       values called at the first timestep (nit000)
      !!
      !! ** input : 
      !!        Namelist namiceini
      !!
      !! history :
      !!  8.5  ! 03-08 (C. Ethe) original code 
      !!  8.5  ! 07-11 (M. Vancoppenolle) rewritten initialization
      !!-----------------------------------------------------------------------------
!{ CREG modifications 
      NAMELIST/namiceini/ ln_iceini, ln_limini_file, rn_thres_sst, rn_hts_ini_n, rn_hts_ini_s, rn_hti_ini_n, rn_hti_ini_s,  &
         &                rn_ati_ini_n, rn_ati_ini_s, rn_smi_ini_n, rn_smi_ini_s, rn_tmi_ini_n, rn_tmi_ini_s,               &
         &                sn_hicif, sn_hsnif, sn_frld, sn_sist,                                                             &
         &                sn_tbif1, sn_tbif2, sn_tbif3, cn_dir
      INTEGER :: ios ,ierr,inum_ice  ! Local integer output status for namelist read
      INTEGER :: ji,jj
      INTEGER :: ifpr, ierror
      !
      CHARACTER(len=100) ::  cn_dir          ! Root directory for location of ice files
      TYPE(FLD_N)                    ::   sn_hicif, sn_hsnif, sn_frld, sn_sist
      TYPE(FLD_N)                    ::   sn_tbif1, sn_tbif2, sn_tbif3
      TYPE(FLD_N), DIMENSION(jpfldi) ::   slf_i                 ! array of namelist informations on the fields to read
!}
      !!-----------------------------------------------------------------------------
      !
      REWIND( numnam_ice_ref )              ! Namelist namiceini in reference namelist : Ice initial state
      READ  ( numnam_ice_ref, namiceini, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namiceini in reference namelist', lwp )

      REWIND( numnam_ice_cfg )              ! Namelist namiceini in configuration namelist : Ice initial state
      READ  ( numnam_ice_cfg, namiceini, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namiceini in configuration namelist', lwp )
      IF(lwm) WRITE ( numoni, namiceini )

!{ CREG modifications 
      slf_i(jp_hicif) = sn_hicif  ;  slf_i(jp_hsnif) = sn_hsnif
      slf_i(jp_frld)  = sn_frld   ;  slf_i(jp_sist)  = sn_sist
      slf_i(jp_tbif1) = sn_tbif1  ;  slf_i(jp_tbif2) = sn_tbif2  ; slf_i(jp_tbif3) = sn_tbif3
!}

      ! Define the initial parameters
      ! -------------------------

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'lim_istate_init : ice parameters inititialisation '
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   initialization with ice (T) or not (F)       ln_iceini     = ', ln_iceini
!{ CREG modifications 
         WRITE(numout,*) '   initialization with ice (T) or not (F)     ln_limini_file  = ', ln_limini_file
!}
         WRITE(numout,*) '   threshold water temp. for initial sea-ice    rn_thres_sst  = ', rn_thres_sst
         WRITE(numout,*) '   initial snow thickness in the north          rn_hts_ini_n  = ', rn_hts_ini_n
         WRITE(numout,*) '   initial snow thickness in the south          rn_hts_ini_s  = ', rn_hts_ini_s 
         WRITE(numout,*) '   initial ice thickness  in the north          rn_hti_ini_n  = ', rn_hti_ini_n
         WRITE(numout,*) '   initial ice thickness  in the south          rn_hti_ini_s  = ', rn_hti_ini_s
         WRITE(numout,*) '   initial ice concentr.  in the north          rn_ati_ini_n  = ', rn_ati_ini_n
         WRITE(numout,*) '   initial ice concentr.  in the north          rn_ati_ini_s  = ', rn_ati_ini_s
         WRITE(numout,*) '   initial  ice salinity  in the north          rn_smi_ini_n  = ', rn_smi_ini_n
         WRITE(numout,*) '   initial  ice salinity  in the south          rn_smi_ini_s  = ', rn_smi_ini_s
         WRITE(numout,*) '   initial  ice/snw temp  in the north          rn_tmi_ini_n  = ', rn_tmi_ini_n
         WRITE(numout,*) '   initial  ice/snw temp  in the south          rn_tmi_ini_s  = ', rn_tmi_ini_s
      ENDIF

!{ CREG modifications 
      IF( ln_limini_file ) THEN                      ! Ice initialization using input file
         !
         ierr = alloc_lim_istate_init()
         !
!         CALL iom_open( 'Ice_initialization.nc', inum_ice )
!         !
!         IF( inum_ice > 0 ) THEN
!            IF(lwp) WRITE(numout,*)
!            IF(lwp) WRITE(numout,*) '                  ice state initialization with : Ice_initialization.nc'
!
!            CALL iom_get( inum_ice, jpdom_data, 'hicif', hicif_ini )
!            CALL iom_get( inum_ice, jpdom_data, 'hsnif', hsnif_ini )
!            CALL iom_get( inum_ice, jpdom_data, 'frld' , frld_ini  )
!            CALL iom_get( inum_ice, jpdom_data, 'ts'   , sist_ini  )
!            CALL iom_get( inum_ice, jpdom_unknown, 'tbif', tbif_ini(1:nlci,1:nlcj,:),   &
!                 &        kstart = (/ mig(1),mjg(1),1 /), kcount = (/ nlci,nlcj,3 /) )
!            ! put some values in the extra-halo...

         ! set si structure
         ALLOCATE( si(jpfldi), STAT=ierror )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'Ice_ini in limistate: unable to allocate si structure' )   ;   RETURN
         ENDIF

         DO ifpr= 1, jpfldi
            ALLOCATE( si(ifpr)%fnow(jpi,jpj,1) )
            ALLOCATE( si(ifpr)%fdta(jpi,jpj,1,2) )
         END DO

         ! fill si with slf_i and control print
         CALL fld_fill( si, slf_i, cn_dir, 'lim_istate', 'lim istate ini', 'numnam_ice' )

         CALL fld_read( nit000, 1, si )                ! input fields provided at the current time-step

         hicif_ini(:,:)  = si(jp_hicif)%fnow(:,:,1)
         hsnif_ini(:,:)  = si(jp_hsnif)%fnow(:,:,1)
         frld_ini(:,:)   = si(jp_frld)%fnow(:,:,1)
         sist_ini(:,:)   = si(jp_sist)%fnow(:,:,1)
         tbif_ini(:,:,1) = si(jp_tbif1)%fnow(:,:,1)
         tbif_ini(:,:,2) = si(jp_tbif2)%fnow(:,:,1)
         tbif_ini(:,:,3) = si(jp_tbif3)%fnow(:,:,1)

         DO jj = nlcj+1, jpj   ;   tbif_ini(1:nlci,jj,:) = tbif_ini(1:nlci,nlej,:)   ;   END DO
         DO ji = nlci+1, jpi   ;   tbif_ini(ji    ,: ,:) = tbif_ini(nlei  ,:   ,:)   ;   END DO

!            CALL iom_close( inum_ice)
!            !
!         ENDIF
      ENDIF
!}

   END SUBROUTINE lim_istate_init

!{ CREG modifications 
   SUBROUTINE limini_file
      !!-----------------------------------------------------------------------------
      !!
      !!
      !!
      !!
      !!-----------------------------------------------------------------------------
      INTEGER    :: jl,ji,jj,jk, jf
      INTEGER    :: jl0
      INTEGER    :: i_fill,jit,jjt
      REAL(wp)   :: ztest_1, ztest_2, ztest_3, ztest_4, ztests, zsigma, zarg, zA, zV, zA_cons, zV_cons, zconv,zH
      REAL(wp)   :: eps=1.e-6
      REAL(wp)   :: zmin,zmax
      !rbb REAL(wp)   :: epsi20,ztmelts,zdh
      REAL(wp)   ::ztmelts,zdh

      REAL(wp), POINTER, DIMENSION(:,:)   :: zhm_i_ini, zat_i_ini, zvt_i_ini, zhm_s_ini, zsm_i_ini
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zv_i_ini
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zht_i_ini,za_i_ini
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zafl_i_ini
      REAL(wp), POINTER, DIMENSION(:,:)   :: zidto    ! ice indicator
       !-----------------------------------------------------------------------------
      IF(lwp)WRITE(numout,*)"limistate: read file : "

      CALL wrk_alloc(jpl,jpi,jpj, zv_i_ini)
      CALL wrk_alloc(    jpi,jpj, zhm_i_ini, zat_i_ini, zvt_i_ini, zhm_s_ini, zsm_i_ini )
      CALL wrk_alloc(    jpl,jpi,jpj,zht_i_ini,za_i_ini)
      CALL wrk_alloc(    jpi,jpj,jpf,zafl_i_ini)
      CALL wrk_alloc(    jpi,jpj,zidto )

      zhm_i_ini(:,:) = hicif_ini(:,:)  ! ice thickness
      zat_i_ini(:,:) = 1._wp - frld_ini(:,:)   ! ice concentration
      zvt_i_ini(:,:) = zhm_i_ini(:,:) * zat_i_ini(:,:)   ! ice volume
      zhm_s_ini(:,:) = hsnif_ini(:,:)  ! snow depth

      zht_i_ini(:,:,:) = 0._wp
      za_i_ini(:,:,:) = 0._wp
      zv_i_ini(:,:,:) = 0._wp
!{ CREG modifications this is required to avoid crash when changing the MPI decomposition
      zsm_i_ini(:,:) = 0._wp
!{ CREG modifications 

      zat_i_ini(:,:) = MIN( zat_i_ini(:,:) , 1.0_wp )
      zafl_i_ini(:,:,:) = 0._wp

      DO ji = 1, jpi
      DO jj = 1, jpj

         IF( zat_i_ini(ji,jj) .GT. 0._wp .AND. zhm_i_ini(ji,jj) .GT. 0._wp )THEN


            IF( gphit(ji,jj) .GE. 0._wp )THEN ; zsm_i_ini(ji,jj) = rn_smi_ini_n
            ELSE                              ; zsm_i_ini(ji,jj) = rn_smi_ini_s
            ENDIF

            jl0 = 1
            DO jl = 2, jpl
               IF ( ( zhm_i_ini(ji,jj) .GT. hi_max(jl-1) ) .AND. &
                  (   zhm_i_ini(ji,jj) .LE. hi_max(jl)   )       ) THEN
               jl0 = jl
               ENDIF
            END DO

            IF( jl0==1 )THEN

               zht_i_ini(1,ji,jj)       = zhm_i_ini(ji,jj)
               za_i_ini(1,ji,jj)        = zat_i_ini(ji,jj)
               zht_i_ini(2:jpl,ji,jj)   = 0._wp
               za_i_ini(2:jpl,ji,jj)    = 0._wp

            ELSE ! jl0 ne 1
               ztest_1 = 0 ; ztest_2 = 0 ; ztest_3 = 0 ; ztest_4 = 0

               DO i_fill = jpl, 1, -1
                  IF( ( ztest_1 + ztest_2 + ztest_3 + ztest_4 ) .NE. 4 ) THEN

                     !----------------------------
                     ! fill the i_fill categories
                     !----------------------------
                     ! *** 1 category to fill
                     IF( i_fill .EQ. 1 ) THEN
                        zht_i_ini(1,ji,jj)       = zhm_i_ini(ji,jj)
                        za_i_ini(1,ji,jj)        = zat_i_ini(ji,jj)
                        zht_i_ini(2:jpl,ji,jj)   = 0._wp
                        za_i_ini(2:jpl,ji,jj)    = 0._wp
                     ELSE

                        ! *** >1 categores to fill
                        !--- Ice thicknesses in the i_fill - 1 first categories
                        DO jl = 1, i_fill - 1
                           zht_i_ini(jl,ji,jj)    = 0.5 * ( hi_max(jl) + hi_max(jl-1) )
                        END DO

                        !--- jl0: most likely index where cc will be maximum
                        DO jl = 1, jpl
                        IF ( ( zhm_i_ini(ji,jj) .GT. hi_max(jl-1) ) .AND. &
                              ( zhm_i_ini(ji,jj) .LE. hi_max(jl)   ) ) THEN
                            jl0 = jl
                        ENDIF
                        END DO
                        jl0 = MIN(jl0, i_fill)

                        !--- Concentrations
                        za_i_ini(jl0,ji,jj)      = zat_i_ini(ji,jj) / SQRT(REAL(jpl))
                        DO jl = 1, i_fill - 1
                        IF ( jl .NE. jl0 ) THEN
                             zsigma               = 0.5 * zhm_i_ini(ji,jj)
                             zarg                 = ( zht_i_ini(jl,ji,jj) - zhm_i_ini(ji,jj) ) / zsigma
                             za_i_ini(jl,ji,jj) = za_i_ini(jl0,ji,jj) * EXP(-zarg**2)
                        ENDIF
                        END DO

                        zA = 0. ! sum of the areas in the jpl categories
                        DO jl = 1, i_fill - 1
                           zA = zA + za_i_ini(jl,ji,jj)
                        END DO
                        za_i_ini(i_fill,ji,jj)   = zat_i_ini(ji,jj) - zA ! ice conc in the last category
                        IF ( i_fill .LT. jpl ) za_i_ini(i_fill+1:jpl, ji,jj) = 0._wp

                        !--- Ice thickness in the last category
                        zV = 0. ! sum of the volumes of the N-1 categories
                        DO jl = 1, i_fill - 1
                           zV = zV + za_i_ini(jl,ji,jj)*zht_i_ini(jl,ji,jj)
                        END DO
                        zht_i_ini(i_fill,ji,jj) = ( zvt_i_ini(ji,jj) - zV ) /za_i_ini(i_fill,ji,jj)
                        IF ( i_fill .LT. jpl ) zht_i_ini(i_fill+1:jpl, ji,jj) = 0._wp

                        !--- volumes
                        zv_i_ini(:,ji,jj) = za_i_ini(:,ji,jj) * zht_i_ini(:,ji,jj)
                        IF ( i_fill .LT. jpl ) zv_i_ini(i_fill+1:jpl, ji,jj) = 0._wp

                     ENDIF ! i_fill

                     !---------------------
                     ! Compatibility tests
                     !---------------------
                     ! Test 1: area conservation
                     zA_cons = SUM(za_i_ini(:,ji,jj)) ; zconv = ABS(zat_i_ini(ji,jj) - zA_cons )
                     IF ( zconv .LT. 1.0e-6 ) THEN
                        ztest_1 = 1
                     ELSE
                        WRITE(*,*) ' * TEST1 AREA NOT CONSERVED *** zA_cons = ', zA_cons,' zat_i_ini = ',zat_i_ini(ji,jj)
                      ! this write is useful
                      !WRITE(numout,*) ' * TEST1 AREA NOT CONSERVED *** zA_cons = ', zA_cons,' zat_i_ini = ',zat_i_ini(ji,jj)
                      !WRITE(numout,*) 'ji,jj,narea ',ji,jj,narea
                      !WRITE(numout,*) ' zat_i_ini : ', zat_i_ini(ji,jj)
                      !WRITE(numout,*) ' zhm_i_ini : ', zhm_i_ini(ji,jj)
                      !WRITE(numout,*) ' zht_i_ini(:,jij,jj) ',zht_i_ini(:,ji,jj)
                      !WRITE(numout,*) ' za_i_ini(:,jij,jj) ',za_i_ini(:,ji,jj)
                      !WRITE(numout,*) ' hi_max ',hi_max
                      !WRITE(numout,*) ' jl0 = ',jl0
                      !WRITE(numout,*) ' vol = ',zvt_i_ini(ji,jj),SUM(zv_i_ini(:,ji,jj))
                      ztest_1 = 0
                     ENDIF

                     ! Test 2: volume conservation
                     zV_cons = SUM(zv_i_ini(:,ji,jj))
                     zconv = ABS(zvt_i_ini(ji,jj) - zV_cons)

                     IF ( zconv .LT. 1.0e-6 ) THEN
                        ztest_2 = 1
                     ELSE
                        ! this write is useful
                        !WRITE(numout,*) ' * TEST2 VOLUME NOT CONSERVED *** zV_cons = ', zV_cons, &
                        !    ' zvt_i_ini = ', zvt_i_ini(ji,jj)
                        !WRITE(numout,*) 'ji,jj,narea ',ji,jj,narea
                        !WRITE(numout,*) ' zat_i_ini : ', zat_i_ini(ji,jj)
                        !WRITE(numout,*) ' zhm_i_ini : ', zhm_i_ini(ji,jj)
                        !WRITE(numout,*) ' zht_i_ini(:,jij,jj) ',zht_i_ini(:,ji,jj)
                        !WRITE(numout,*) ' za_i_ini(:,jij,jj) ',za_i_ini(:,ji,jj)
                        !WRITE(numout,*) ' hi_max ',hi_max
                        !WRITE(numout,*) ' jl0 = ',jl0
                        ztest_2 = 0
                     ENDIF

                     ! Test 3: thickness of the last category is in-bounds ? 
                     IF ( zht_i_ini(i_fill, ji,jj) .GT. hi_max(i_fill-1) ) THEN
                     ztest_3 = 1
                     ELSE
                     ! this write is useful
                     !WRITE(numout,*) ' * TEST 3 THICKNESS OF THE LAST CATEGORY OUT OF BOUNDS *** zht_i_ini(i_fill,ji,jj) = ', &
                     !zht_i_ini(i_fill,ji,jj), ' hi_max(jpl-1) = ', hi_max(i_fill-1)
                     !WRITE(numout,*) 'ji,jj,narea ',ji,jj,narea
                     !WRITE(numout,*) ' zat_i_ini : ', zat_i_ini(ji,jj)
                     !WRITE(numout,*) ' zhm_i_ini : ', zhm_i_ini(ji,jj)
                     !WRITE(numout,*) ' zht_i_ini(:,jij,jj) ',zht_i_ini(:,ji,jj)
                     !WRITE(numout,*) ' za_i_ini(:,jij,jj) ',za_i_ini(:,ji,jj)
                     !WRITE(numout,*) ' hi_max ',hi_max
                     !WRITE(numout,*) ' jl0 = ',jl0
                     ztest_3 = 0
                     ENDIF

                     ! Test 4: positivity of ice concentrations
                     ztest_4 = 1
                     DO jl = 1, jpl
                     IF ( za_i_ini(jl,ji,jj) .LT. 0._wp ) THEN
                        ! this write is useful
                        !WRITE(numout,*) ' * TEST 4 POSITIVITY NOT OK FOR CAT ', jl, 'WITH A = ', za_i_ini(jl,ji,jj)
                        !WRITE(numout,*) 'ji,jj,narea ',ji,jj,narea
                        !WRITE(numout,*) ' zat_i_ini : ', zat_i_ini(ji,jj)
                        !WRITE(numout,*) ' zhm_i_ini : ', zhm_i_ini(ji,jj)
                        !WRITE(numout,*) ' zht_i_ini(:,jij,jj) ',zht_i_ini(:,ji,jj)
                        !WRITE(numout,*) ' za_i_ini(:,jij,jj) ',za_i_ini(:,ji,jj)
                        !WRITE(numout,*) ' hi_max ',hi_max
                        !WRITE(numout,*) ' jl0 = ',jl0
                        !WRITE(numout,*)
                        ztest_4 = 0
                     ENDIF
                     END DO

                  ENDIF ! ztest_1 + ztest_2 + ztest_3 + ztest_4

                  ztests = ztest_1 + ztest_2 + ztest_3 + ztest_4

               END DO ! i_fill

               !WRITE(numout,*) ' ztests : ', ztests
               !IF ( ztests .NE. 4 ) THEN
               !WRITE(numout,*)
               !WRITE(numout,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
               !WRITE(numout,*) ' !!!! RED ALERT                  !!! '
               !WRITE(numout,*) ' !!!! BIIIIP BIIIP BIIIIP BIIIIP !!!'
               !WRITE(numout,*) ' !!!! Something is wrong in the LIM3 initialization procedure '
               !WRITE(numout,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
               !WRITE(numout,*) 'ji,jj,narea ',ji,jj,narea
               !WRITE(numout,*) ' *** ztests is not equal to 4 '
               !WRITE(numout,*) ' *** ztest_i (i=1,4) = ', ztest_1, ztest_2,ztest_3,ztest_4
               !WRITE(numout,*) ' zat_i_ini : ', zat_i_ini(ji,jj)
               !WRITE(numout,*) ' zhm_i_ini : ', zhm_i_ini(ji,jj)
               !WRITE(numout,*) ' zht_i_ini(:,jij,jj) ',zht_i_ini(:,ji,jj)
               !WRITE(numout,*) ' za_i_ini(:,jij,jj) ',za_i_ini(:,ji,jj)
               !WRITE(numout,*) ' hi_max ',hi_max
               !ENDIF ! ztests .NE. 4

            ENDIF  !  jl0 ne 1
         zafl_i_ini(ji,jj,jpf) = zat_i_ini(ji,jj) ! Only the last category is filled since floes are supposed to be unbroken.
         zA_cons = SUM(zafl_i_ini(ji,jj,:)) ; zconv = ABS(zat_i_ini(ji,jj) - zA_cons )
         IF ( zconv .GT. 1.0e-6 ) THEN
            WRITE(*,*) ' * TEST1 AREA FLOE NOT CONSERVED *** zA_cons = ', zA_cons,' zat_i_ini = ',zat_i_ini(ji,jj) 
         ENDIF
 
         ENDIF  !  zat_i_ini ne 0
      END DO ! jj
      END DO ! ji


      !---------------------------------------------------------------------
      ! 3.3) Space-dependent arrays for ice state variables
      !---------------------------------------------------------------------

      ! Ice concentration, thickness and volume, ice salinity, ice age, surface
      ! temperature
      DO jl = 1, jpl ! loop over categories
         DO jj = 1, jpj
            DO ji = 1, jpi
               a_i(ji,jj,jl)   = zswitch(ji,jj) * za_i_ini (jl,ji,jj)  ! concentration
               ht_i(ji,jj,jl)  = zswitch(ji,jj) * zht_i_ini(jl,ji,jj)   !ice thickness

               IF( zhm_i_ini( ji,jj ) .GT. 0_wp )THEN ; ht_s(ji,jj,jl)  = ht_i(ji,jj,jl) * ( zhm_s_ini( ji,jj ) / zhm_i_ini( ji,jj ) )
               ELSE                                   ; ht_s(ji,jj,jl)  = 0._wp
               ENDIF
               sm_i(ji,jj,jl)  = zswitch(ji,jj) * zsm_i_ini(ji,jj) !+ (1._wp - zswitch(ji,jj) ) * rn_simin ! salinity
               o_i(ji,jj,jl)   = zswitch(ji,jj) * 1._wp + ( 1._wp -zswitch(ji,jj) ) ! age
               t_su(ji,jj,jl)  = sist_ini(ji,jj)

               ! This case below should not be used if (ht_s/ht_i) is ok in
               ! namelist
               ! In case snow load is in excess that would lead to
               ! transformation from snow to ice
               ! Then, transfer the snow excess into the ice (different from
               ! limthd_dh)
               zdh = MAX( 0._wp, ( rhosn * ht_s(ji,jj,jl) + ( rhoic - rau0 ) *ht_i(ji,jj,jl) ) * r1_rau0 )
               ! recompute ht_i, ht_s avoiding out of bounds values
               ht_i(ji,jj,jl) = MIN( hi_max(jl), ht_i(ji,jj,jl) + zdh )
               ht_s(ji,jj,jl) = MAX( 0._wp, ht_s(ji,jj,jl) - zdh * rhoic *r1_rhosn )

               ! ice volume, salt content, age content
               v_i(ji,jj,jl)   = ht_i(ji,jj,jl) * a_i(ji,jj,jl)              !ice volume
               v_s(ji,jj,jl)   = ht_s(ji,jj,jl) * a_i(ji,jj,jl)              !snow volume
               smv_i(ji,jj,jl) = MIN( sm_i(ji,jj,jl) , sss_m(ji,jj) ) *v_i(ji,jj,jl) ! salt content
               oa_i(ji,jj,jl)  = o_i(ji,jj,jl) * a_i(ji,jj,jl)               !age content
            END DO ! ji
         END DO ! jj
      END DO ! jl

      DO jf = 1, jpf ! loop over ice floes categories (GB)
         DO jj = 1, jpj
            DO ji = 1, jpi
              IF (ABS(SUM(a_i(ji,jj,:))-SUM(za_i_ini (:,ji,jj))*zswitch(ji,jj))>1e-6) WRITE(*,*) "GROS PROBLEME 2",SUM(a_i(ji,jj,:)),SUM(za_i_ini (:,ji,jj)),zswitch(ji,jj)
              afl_i(ji,jj,jf)   = zswitch(ji,jj) * zafl_i_ini(ji,jj,jf)
            END DO ! ji
         END DO ! jj
      END DO ! jf

         DO jj = 1, jpj
            DO ji = 1, jpi
           !   IF (SUM(zafl_i_ini(ji,jj,:))>0) WRITE(*,*) "AIE AIE AIE", SUM(zafl_i_ini(ji,jj,:)),SUM(afl_i(ji,jj,:)),zswitch(ji,jj),SUM(a_i(ji,jj,:))
            IF (ABS(SUM(afl_i(ji,jj,:))-SUM(a_i(ji,jj,:)))>1e-6) WRITE(*,*) "AIE AIE AIE", SUM(zafl_i_ini(ji,jj,:)),SUM(afl_i(ji,jj,:)),zswitch(ji,jj),SUM(a_i(ji,jj,:))
            END DO
         END DO 
      !cbr
      DO jk = 1, nlay_s
         DO  jl = 1, jpl ! loop over categories
            !rbb t_s(:,:,1,jl) =  tbif_ini(:,:,1)
            t_s(:,:,1,jl) =  tbif_ini(:,:,1)*zswitch(:,:)+ ( 1._wp - zswitch(:,:) ) * rt0
         END DO ! jl
      END DO ! jk

      ! Snow temperature and heat content
      DO jk = 1, nlay_s
         DO jl = 1, jpl ! loop over categories
            DO jj = 1, jpj
               DO ji = 1, jpi
!cbr???                   t_s(ji,jj,jk,jl) = zswitch(ji,jj) * ztm_i_ini(ji,jj) + ( 1._wp - zswitch(ji,jj) ) * rt0
                   ! Snow energy of melting
                   e_s(ji,jj,jk,jl) = zswitch(ji,jj) * rhosn * ( cpic * ( rt0 - t_s(ji,jj,jk,jl) ) + lfus )

                   ! Mutliply by volume, and divide by number of layers to get
                   ! heat content in J/m2
                   e_s(ji,jj,jk,jl) = e_s(ji,jj,jk,jl) * v_s(ji,jj,jl) *r1_nlay_s
               END DO ! ji
            END DO ! jj
         END DO ! jl
      END DO ! jk

      ! Ice salinity, temperature and heat content
      DO  jk = 1, nlay_i
         DO jl = 1, jpl ! loop over categories
            DO jj = 1, jpj
               DO ji = 1, jpi
!cbr???                   t_i(ji,jj,jk,jl) = zswitch(ji,jj) * ztm_i_ini(ji,jj) + ( 1._wp - zswitch(ji,jj) ) * rt0
                   t_i(ji,jj,jk,jl) =  tbif_ini(ji,jj,2)*zswitch(ji,jj)+ ( 1._wp - zswitch(ji,jj) ) * rt0
                   s_i(ji,jj,jk,jl) = zswitch(ji,jj) * zsm_i_ini(ji,jj) !+ ( 1._wp - zswitch(ji,jj) ) * rn_simin
                   ztmelts          = - tmut * s_i(ji,jj,jk,jl) + rt0           !Melting temperature in K

                   ! heat content per unit volume
                   e_i(ji,jj,jk,jl) = zswitch(ji,jj) * rhoic * (   cpic    * ( ztmelts - t_i(ji,jj,jk,jl) ) &
                      +   lfus    * ( 1._wp - (ztmelts-rt0) /MIN((t_i(ji,jj,jk,jl)-rt0),-epsi20) ) &
                      -   rcp     * ( ztmelts - rt0 ) )

                   ! Mutliply by ice volume, and divide by number of layers to
                   ! get heat content in J/m2
                   e_i(ji,jj,jk,jl) = e_i(ji,jj,jk,jl) * v_i(ji,jj,jl) * r1_nlay_i
               END DO ! ji
            END DO ! jj
         END DO ! jl
      END DO ! jk

      !cbr tmp CALL wrk_dealloc(jpl,jpi,jpj, zht_i_ini, za_i_ini, zv_i_ini)
      CALL wrk_dealloc(jpl,jpi,jpj, zv_i_ini)
      CALL wrk_dealloc(    jpl,jpi,jpj,zht_i_ini,za_i_ini)
      CALL wrk_dealloc(    jpi,jpj, zhm_i_ini, zat_i_ini, zvt_i_ini, zhm_s_ini,zsm_i_ini )
      CALL wrk_dealloc(    jpi,jpj,zidto )

  END SUBROUTINE limini_file
!}


  INTEGER FUNCTION alloc_lim_istate_init()
      !!-----------------------------------------------------------------------------
      !!
      !!
      !!
      !!
      !!-----------------------------------------------------------------------------
      INTEGER :: ierr(1)
      !!-----------------------------------------------------------------------------
      ALLOCATE( hicif_ini(jpi,jpj) , hsnif_ini(jpi,jpj) , frld_ini(jpi,jpj) , sist_ini(jpi,jpj) , tbif_ini(jpi,jpj,3) , Stat=ierr(1) )
      alloc_lim_istate_init = MAXVAL(ierr)
      IF( alloc_lim_istate_init /= 0 )   CALL ctl_warn( 'lim_istate_init: failed to allocate arrays')

   END FUNCTION alloc_lim_istate_init
#else
   !!----------------------------------------------------------------------
   !!   Default option :         Empty module          NO LIM sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_istate          ! Empty routine
   END SUBROUTINE lim_istate
#endif

   !!======================================================================
END MODULE limistate

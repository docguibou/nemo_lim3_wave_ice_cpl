MODULE nemogcm
!!======================================================================
!!                       ***  MODULE nemogcm   ***
!! StandAlone Surface module : surface fluxes + sea-ice + iceberg floats
!!======================================================================
!! History :  3.6  ! 2011-11  (S. Alderson, G. Madec) original code
!!             -   ! 2013-06  (I. Epicoco, S. Mocavero, CMCC) nemo_northcomms: setup avoiding MPI communication
!!             -   ! 2014-12  (G. Madec) remove KPP scheme and cross-land advection (cla)
!!            4.0  ! 2016-10  (G. Madec, S. Flavoni)  domain configuration / user defined interface
!!----------------------------------------------------------------------

!!----------------------------------------------------------------------
!!   nemo_gcm      : solve ocean dynamics, tracer, biogeochemistry and/or sea-ice
!!   nemo_init     : initialization of the NEMO system
!!   nemo_ctl      : initialisation of the contol print
!!   nemo_closefile: close remaining open files
!!   nemo_alloc    : dynamical allocation
!!   nemo_partition: calculate MPP domain decomposition
!!   factorise     : calculate the factors of the no. of MPI processes
!!----------------------------------------------------------------------
   USE step_oce       ! module used in the ocean time stepping module
   USE sbc_oce        ! surface boundary condition: ocean
   USE phycst         ! physical constant                  (par_cst routine)
   USE domain         ! domain initialization   (dom_init & dom_cfg routines)
   USE usrdef_nam     ! user defined configuration
   USE daymod         ! calendar
   USE step           ! NEMO time-stepping                 (stp     routine)
   USE cpl_oasis3     !
   USE sbcssm         !
   USE icbini         ! handle bergs, initialisation
   USE icbstp         ! handle bergs, calving, themodynamics and transport
   USE bdyini         ! open boundary cond. setting       (bdy_init routine). clem: mandatory for LIM3
   USE bdydta         ! open boundary cond. setting   (bdy_dta_init routine). clem: mandatory for LIM3
!
   USE lib_mpp        ! distributed memory computing
   USE mppini         ! shared/distributed memory setting (mpp_init routine)
   USE lbcnfd   , ONLY: isendto, nsndto, nfsloop, nfeloop ! Setup of north fold exchanges
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)

   USE xios           ! xIOserver


   IMPLICIT NONE
   PRIVATE

   PUBLIC   nemo_gcm    ! called by model.F90
   PUBLIC   nemo_init   ! needed by AGRIF

   CHARACTER(lc) ::   cform_aaa="( /, 'AAAAAAAA', / ) "     ! flag for output listing

!!----------------------------------------------------------------------
!! NEMO/OPA 4.0 , NEMO Consortium (2016)
!! $Id: nemogcm.F90 7761 2017-03-06 17:58:35Z clem $
!! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE nemo_gcm
!!----------------------------------------------------------------------
!!                     ***  ROUTINE nemo_gcm  ***
!!
!! ** Purpose :   NEMO solves the primitive equations on an orthogonal
!!              curvilinear mesh on the sphere.
!!
!! ** Method  : - model general initialization
!!              - launch the time-stepping (stp routine)
!!              - finalize the run by closing files and communications
!!
!! References : Madec, Delecluse, Imbard, and Levy, 1997:  internal report, IPSL.
!!              Madec, 2008, internal report, IPSL.
!!----------------------------------------------------------------------
      INTEGER ::   istp   ! time step index
!!----------------------------------------------------------------------
!

!
!                            !-----------------------!
      CALL nemo_init               !==  Initialisations  ==!
!                            !-----------------------!

! check that all process are still there... If some process have an error,
! they will never enter in step and other processes will wait until the end of the cpu time!
      IF( lk_mpp )   CALL mpp_max( nstop )

      IF(lwp) WRITE(numout,cform_aaa)   ! Flag AAAAAAA

!                            !-----------------------!
!                            !==   time stepping   ==!
!                            !-----------------------!
      istp = nit000


      DO WHILE ( istp <= nitend .AND. nstop == 0 )

         IF ( .NOT. ln_diurnal_only ) THEN
            CALL stp( istp )                 ! standard time stepping
         ELSE
            CALL stp_diurnal( istp )        ! time step only the diurnal SST
         ENDIF 

         istp = istp + 1
         IF( lk_mpp )   CALL mpp_max( nstop )
         END DO
!
      IF( ln_icebergs )   CALL icb_end( nitend )

!                            !------------------------!
!                            !==  finalize the run  ==!
!                            !------------------------!
      IF(lwp) WRITE(numout,cform_aaa)        ! Flag AAAAAAA
!
      IF( nstop /= 0 .AND. lwp ) THEN        ! error print
         WRITE(numout,cform_err)
         WRITE(numout,*) nstop, ' error have been found'
      ENDIF
!

      IF( nn_timing == 1 )   CALL timing_finalize
!
      CALL nemo_closefile
!

      CALL xios_finalize                     ! end mpp communications with xios
      IF( lk_oasis )   CALL cpl_finalize     ! end coupling and mpp communications with OASIS

!
   END SUBROUTINE nemo_gcm


   SUBROUTINE nemo_init
!!----------------------------------------------------------------------
!!                     ***  ROUTINE nemo_init  ***
!!
!! ** Purpose :   initialization of the NEMO GCM
!!----------------------------------------------------------------------
      INTEGER  ::   ji            ! dummy loop indices
      INTEGER  ::   ilocal_comm   ! local integer
      INTEGER  ::   ios, inum     !   -      -
      CHARACTER(len=120), DIMENSION(30) ::   cltxt, cltxt2, clnam
      CHARACTER(len=80)                 ::   clname
!
      NAMELIST/namctl/ ln_ctl   , nn_print, nn_ictls, nn_ictle,   &
         &             nn_isplt , nn_jsplt, nn_jctls, nn_jctle,   &
         &             nn_timing, nn_diacfl
      NAMELIST/namcfg/ ln_read_cfg, cn_domcfg, ln_write_cfg, cn_domcfg_out, ln_use_jattr
!!----------------------------------------------------------------------
!
      cltxt  = ''
      cltxt2 = ''
      clnam  = ''  
      cxios_context = 'nemo'
!
!                             ! Open reference namelist and configuration namelist files
      IF( lk_oasis ) THEN 
         CALL ctl_opn( numnam_ref, 'namelist_sas_ref', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE. )
         CALL ctl_opn( numnam_cfg, 'namelist_sas_cfg', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE. )
         cxios_context = 'sas'
         clname = 'output.namelist_sas.dyn'
      ELSE
         CALL ctl_opn( numnam_ref, 'namelist_ref', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE. )
         CALL ctl_opn( numnam_cfg, 'namelist_cfg', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE. )
         cxios_context = 'nemo'
         clname = 'output.namelist.dyn'
   ENDIF
!
      REWIND( numnam_ref )              ! Namelist namctl in reference namelist : Control prints
      READ  ( numnam_ref, namctl, IOSTAT = ios, ERR = 901 )
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namctl in reference namelist', .TRUE. )
!
      REWIND( numnam_cfg )              ! Namelist namctl in confguration namelist
      READ  ( numnam_cfg, namctl, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namctl in configuration namelist', .TRUE. )
!
      REWIND( numnam_ref )              ! Namelist namcfg in reference namelist : Control prints
      READ  ( numnam_ref, namcfg, IOSTAT = ios, ERR = 903 )
903   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namcfg in reference namelist', .TRUE. )

      REWIND( numnam_cfg )              ! Namelist namcfg in confguration namelist : Control prints & Benchmark
      READ  ( numnam_cfg, namcfg, IOSTAT = ios, ERR = 904 )
904   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namcfg in configuration namelist', .TRUE. )   

!                             !--------------------------!
!                             !  Set global domain size  !   (control print return in cltxt2)
!                             !--------------------------!
      IF( ln_read_cfg ) THEN              ! Read sizes in domain configuration file
         CALL domain_cfg ( cltxt2,        cn_cfg, nn_cfg, jpiglo, jpjglo, jpkglo, jperio )
!
      ELSE                                ! user-defined namelist
         CALL usr_def_nam( cltxt2, clnam, cn_cfg, nn_cfg, jpiglo, jpjglo, jpkglo, jperio )
      ENDIF
!
      jpk = jpkglo
!

!
!                             !--------------------------------------------!
!                             !  set communicator & select the local node  !
!                             !  NB: mynode also opens output.namelist.dyn !
!                             !      on unit number numond on first proc   !
!                             !--------------------------------------------!

      IF( Agrif_Root() ) THEN
         IF( lk_oasis ) THEN
            CALL cpl_init( "sas", ilocal_comm )                          ! nemo local communicator given by oasis
            CALL xios_initialize( "not used",local_comm=ilocal_comm )    ! send nemo communicator to xios
         ELSE
            CALL xios_initialize( "for_xios_mpi_id",return_comm=ilocal_comm )    ! nemo local communicator given by xios
         ENDIF
      ENDIF
      narea = mynode ( cltxt, clname, numnam_ref, numnam_cfg, numond , nstop, ilocal_comm )  ! Nodes selection


      narea = narea + 1                                     ! mynode return the rank of proc (0 --> jpnij -1 )

      lwm = (narea == 1)                                    ! control of output namelists
      lwp = (narea == 1) .OR. ln_ctl                        ! control of all listing output print

      IF(lwm) THEN               ! write merged namelists from earlier to output namelist
!                       ! now that the file has been opened in call to mynode.
!                       ! NB: nammpp has already been written in mynode (if lk_mpp_mpi)
         WRITE( numond, namctl )
         WRITE( numond, namcfg )
         IF( .NOT.ln_read_cfg ) THEN
            DO ji = 1, SIZE(clnam)
               IF( TRIM(clnam(ji)) /= '' )   WRITE(numond, * ) clnam(ji)     ! namusr_def print
            END DO
         ENDIF
      ENDIF

! If dimensions of processor grid weren't specified in the namelist file
! then we calculate them here now that we have our communicator size
      IF( jpni < 1 .OR. jpnj < 1 ) THEN

         IF( Agrif_Root() )   CALL nemo_partition( mppsize )

      ENDIF

      IF( Agrif_Root() ) THEN       ! AGRIF mother: specific setting from jpni and jpnj

         jpi = ( jpiglo     -2*jpreci + (jpni-1) ) / jpni + 2*jpreci    ! first  dim.
         jpj = ( jpjglo     -2*jprecj + (jpnj-1) ) / jpnj + 2*jprecj    ! second dim.

      ENDIF


      jpim1 = jpi-1                                            ! inner domain indices
      jpjm1 = jpj-1                                            !   "           "
      jpkm1 = MAX( 1, jpk-1 )                                  !   "           "
      jpij  = jpi*jpj                                          !  jpi x j

      IF(lwp) THEN                            ! open listing units
!
         IF( lk_oasis ) THEN
            CALL ctl_opn( numout,   'sas.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
         ELSE
            CALL ctl_opn( numout, 'ocean.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
         ENDIF
!
         WRITE(numout,*)
         WRITE(numout,*) '   CNRS - NERC - Met OFFICE - MERCATOR-ocean - INGV - CMCC'
         WRITE(numout,*) '                       NEMO team'
         WRITE(numout,*) '            Ocean General Circulation Model'
         WRITE(numout,*) '                  version 3.7  (2016) '
         WRITE(numout,*) '             StandAlone Surface version (SAS) '
         WRITE(numout,*)
         WRITE(numout,*)
         DO ji = 1, SIZE(cltxt)
            IF( TRIM(cltxt (ji)) /= '' )   WRITE(numout,*) cltxt(ji)    ! control print of mynode
         END DO
         WRITE(numout,*)
         WRITE(numout,*)
         DO ji = 1, SIZE(cltxt2)
            IF( TRIM(cltxt2(ji)) /= '' )   WRITE(numout,*) cltxt2(ji)   ! control print of domain size
         END DO
!
         WRITE(numout,cform_aaa)                                        ! Flag AAAAAAA
!
      ENDIF

! Now we know the dimensions of the grid and numout has been set: we can allocate arrays
      CALL nemo_alloc()
!                             !-------------------------------!
!                             !  NEMO general initialization  !
!                             !-------------------------------!

      CALL nemo_ctl                          ! Control prints

!                                      ! Domain decomposition
      IF( jpni*jpnj == jpnij ) THEN   ;   CALL mpp_init      ! standard cutting out
      ELSE                            ;   CALL mpp_init2     ! eliminate land processors
      ENDIF
!
      IF( nn_timing == 1 )  CALL timing_init
!
!                                     ! General initialization
                            CALL phy_cst    ! Physical constants
                            CALL eos_init   ! Equation of state
                            CALL dom_init   ! Domain

     IF( ln_nnogather )     CALL nemo_northcomms   ! Initialise the northfold neighbour lists (must be done after the masks are defined)

      IF( ln_ctl      )     CALL prt_ctl_init   ! Print control
                            CALL day_init   ! model calendar (using both namelist and restart infos)

                            CALL sbc_init   ! Forcings : surface module

! ==> clem: open boundaries init. is mandatory for LIM3 because ice BDY is not decoupled from
!           the environment of ocean BDY. Therefore bdy is called in both OPA and SAS modules.
!           This is not clean and should be changed in the future.
                            CALL bdy_init
! ==>
                            CALL icb_init( rdt, nit000)   ! initialise icebergs instance
      
      IF(lwp) WRITE(numout,*) 'Euler time step switch is ', neuler
!
   END SUBROUTINE nemo_init


   SUBROUTINE nemo_ctl
!!----------------------------------------------------------------------
!!                     ***  ROUTINE nemo_ctl  ***
!!
!! ** Purpose :   control print setting
!!
!! ** Method  : - print namctl information and check some consistencies
!!----------------------------------------------------------------------
!
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'nemo_ctl: Control prints'
         WRITE(numout,*) '~~~~~~~ '
         WRITE(numout,*) '   Namelist namctl'
         WRITE(numout,*) '      run control (for debugging)     ln_ctl     = ', ln_ctl
         WRITE(numout,*) '      level of print                  nn_print   = ', nn_print
         WRITE(numout,*) '      Start i indice for SUM control  nn_ictls   = ', nn_ictls
         WRITE(numout,*) '      End i indice for SUM control    nn_ictle   = ', nn_ictle
         WRITE(numout,*) '      Start j indice for SUM control  nn_jctls   = ', nn_jctls
         WRITE(numout,*) '      End j indice for SUM control    nn_jctle   = ', nn_jctle
         WRITE(numout,*) '      number of proc. following i     nn_isplt   = ', nn_isplt
         WRITE(numout,*) '      number of proc. following j     nn_jsplt   = ', nn_jsplt
         WRITE(numout,*) '      timing activated    (0/1)       nn_timing  = ', nn_timing
      ENDIF
!
      nprint    = nn_print          ! convert DOCTOR namelist names into OLD names
      nictls    = nn_ictls
      nictle    = nn_ictle
      njctls    = nn_jctls
      njctle    = nn_jctle
      isplt     = nn_isplt
      jsplt     = nn_jsplt

      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'namcfg  : configuration initialization through namelist read'
         WRITE(numout,*) '~~~~~~~ '
         WRITE(numout,*) '   Namelist namcfg'
         WRITE(numout,*) '      read domain configuration files               ln_read_cfg      = ', ln_read_cfg
         WRITE(numout,*) '         filename to be read                           cn_domcfg     = ', TRIM(cn_domcfg)
         WRITE(numout,*) '      write  configuration definition files         ln_write_cfg     = ', ln_write_cfg
         WRITE(numout,*) '         filename to be written                        cn_domcfg_out = ', TRIM(cn_domcfg_out)
         WRITE(numout,*) '      use file attribute if exists as i/p j-start   ln_use_jattr     = ', ln_use_jattr
      ENDIF
!                             ! Parameter control
!
      IF( ln_ctl ) THEN                 ! sub-domain area indices for the control prints
         IF( lk_mpp .AND. jpnij > 1 ) THEN
            isplt = jpni   ;   jsplt = jpnj   ;   ijsplt = jpni*jpnj   ! the domain is forced to the real split domain
         ELSE
            IF( isplt == 1 .AND. jsplt == 1  ) THEN
               CALL ctl_warn( ' - isplt & jsplt are equal to 1',   &
                  &           ' - the print control will be done over the whole domain' )
            ENDIF
            ijsplt = isplt * jsplt            ! total number of processors ijsplt
         ENDIF
         IF(lwp) WRITE(numout,*)'          - The total number of processors over which the'
         IF(lwp) WRITE(numout,*)'            print control will be done is ijsplt : ', ijsplt
!
!                              ! indices used for the SUM control
         IF( nictls+nictle+njctls+njctle == 0 )   THEN    ! print control done over the default area
            lsp_area = .FALSE.
         ELSE                                             ! print control done over a specific  area
            lsp_area = .TRUE.
            IF( nictls < 1 .OR. nictls > jpiglo )   THEN
               CALL ctl_warn( '          - nictls must be 1<=nictls>=jpiglo, it is forced to 1' )
               nictls = 1
            ENDIF
            IF( nictle < 1 .OR. nictle > jpiglo )   THEN
               CALL ctl_warn( '          - nictle must be 1<=nictle>=jpiglo, it is forced to jpiglo' )
               nictle = jpiglo
            ENDIF
            IF( njctls < 1 .OR. njctls > jpjglo )   THEN
               CALL ctl_warn( '          - njctls must be 1<=njctls>=jpjglo, it is forced to 1' )
               njctls = 1
            ENDIF
            IF( njctle < 1 .OR. njctle > jpjglo )   THEN
               CALL ctl_warn( '          - njctle must be 1<=njctle>=jpjglo, it is forced to jpjglo' )
               njctle = jpjglo
            ENDIF
         ENDIF
      ENDIF
!
      IF( 1_wp /= SIGN(1._wp,-0._wp)  )   CALL ctl_stop( 'nemo_ctl: The intrinsec SIGN function follows ',  &
         &                                               'f2003 standard. '                              ,  &
         &                                               'Compile with key_nosignedzero enabled' )
!
   END SUBROUTINE nemo_ctl


   SUBROUTINE nemo_closefile
!!----------------------------------------------------------------------
!!                     ***  ROUTINE nemo_closefile  ***
!!
!! ** Purpose :   Close the files
!!----------------------------------------------------------------------
!
      IF( lk_mpp )   CALL mppsync
!
      CALL iom_close                                 ! close all input/output files managed by iom_*
!
      IF( numstp          /= -1 )   CLOSE( numstp      )   ! time-step file
      IF( numnam_ref      /= -1 )   CLOSE( numnam_ref      )   ! oce reference namelist
      IF( numnam_cfg      /= -1 )   CLOSE( numnam_cfg      )   ! oce configuration namelist
      IF( lwm.AND.numond  /= -1 )   CLOSE( numond          )   ! oce output namelist
      IF( numnam_ice_ref  /= -1 )   CLOSE( numnam_ice_ref  )   ! ice reference namelist
      IF( numnam_ice_cfg  /= -1 )   CLOSE( numnam_ice_cfg  )   ! ice configuration namelist
      IF( lwm.AND.numoni  /= -1 )   CLOSE( numoni          )   ! ice output namelist
      IF( numevo_ice      /= -1 )   CLOSE( numevo_ice  )   ! ice variables (temp. evolution)
      IF( numout          /=  6 )   CLOSE( numout      )   ! standard model output file
!
      numout = 6                                     ! redefine numout in case it is used after this point...
!
   END SUBROUTINE nemo_closefile


   SUBROUTINE nemo_alloc
!!----------------------------------------------------------------------
!!                     ***  ROUTINE nemo_alloc  ***
!!
!! ** Purpose :   Allocate all the dynamic arrays of the OPA modules
!!
!! ** Method  :
!!----------------------------------------------------------------------
      USE diawri    , ONLY: dia_wri_alloc
      USE dom_oce   , ONLY: dom_oce_alloc
      USE bdy_oce   , ONLY: ln_bdy, bdy_oce_alloc
      USE oce         ! clem: mandatory for LIM3 because needed for bdy arrays
!
      INTEGER :: ierr
!!----------------------------------------------------------------------
!
      ierr =        dia_wri_alloc   ()
      ierr = ierr + dom_oce_alloc   ()          ! ocean domain
      ierr = ierr + oce_alloc       ()          ! (tsn...) needed for agrif and/or lim3 and bdy
      ierr = ierr + bdy_oce_alloc   ()          ! bdy masks (incl. initialization)
!
      IF( lk_mpp    )   CALL mpp_sum( ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'nemo_alloc : unable to allocate standard ocean arrays' )
!
   END SUBROUTINE nemo_alloc


   SUBROUTINE nemo_partition( num_pes )
!!----------------------------------------------------------------------
!!                 ***  ROUTINE nemo_partition  ***
!!
!! ** Purpose :
!!
!! ** Method  :
!!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   num_pes   ! The number of MPI processes we have
!
      INTEGER, PARAMETER :: nfactmax = 20
      INTEGER :: nfact ! The no. of factors returned
      INTEGER :: ierr  ! Error flag
      INTEGER :: ji
      INTEGER :: idiff, mindiff, imin ! For choosing pair of factors that are closest in value
      INTEGER, DIMENSION(nfactmax) :: ifact ! Array of factors
!!----------------------------------------------------------------------
!
      ierr = 0
!
      CALL factorise( ifact, nfactmax, nfact, num_pes, ierr )
!
      IF( nfact <= 1 ) THEN
         WRITE (numout, *) 'WARNING: factorisation of number of PEs failed'
         WRITE (numout, *) '       : using grid of ',num_pes,' x 1'
         jpnj = 1
         jpni = num_pes
      ELSE
! Search through factors for the pair that are closest in value
         mindiff = 1000000
         imin    = 1
         DO ji = 1, nfact-1, 2
            idiff = ABS( ifact(ji) - ifact(ji+1) )
            IF( idiff < mindiff ) THEN
               mindiff = idiff
               imin = ji
            ENDIF
         END DO
         jpnj = ifact(imin)
         jpni = ifact(imin + 1)
      ENDIF
!
      jpnij = jpni*jpnj
!
   END SUBROUTINE nemo_partition


   SUBROUTINE factorise( kfax, kmaxfax, knfax, kn, kerr )
!!----------------------------------------------------------------------
!!                     ***  ROUTINE factorise  ***
!!
!! ** Purpose :   return the prime factors of n.
!!                knfax factors are returned in array kfax which is of
!!                maximum dimension kmaxfax.
!! ** Method  :
!!----------------------------------------------------------------------
      INTEGER                    , INTENT(in   ) ::   kn, kmaxfax
      INTEGER                    , INTENT(  out) ::   kerr, knfax
      INTEGER, DIMENSION(kmaxfax), INTENT(  out) ::   kfax
!
      INTEGER :: ifac, jl, inu
      INTEGER, PARAMETER :: ntest = 14
      INTEGER, DIMENSION(ntest) ::   ilfax
!!----------------------------------------------------------------------
!
! lfax contains the set of allowed factors.
      ilfax(:) = (/(2**jl,jl=ntest,1,-1)/)
!
! Clear the error flag and initialise output vars
      kerr  = 0
      kfax  = 1
      knfax = 0
!
! Find the factors of n.
      IF( kn == 1 )   GOTO 20

! nu holds the unfactorised part of the number.
! knfax holds the number of factors found.
! l points to the allowed factor list.
! ifac holds the current factor.
!
      inu   = kn
      knfax = 0
!
      DO jl = ntest, 1, -1
!
         ifac = ilfax(jl)
         IF( ifac > inu )   CYCLE

! Test whether the factor will divide.

         IF( MOD(inu,ifac) == 0 ) THEN
!
            knfax = knfax + 1            ! Add the factor to the list
            IF( knfax > kmaxfax ) THEN
               kerr = 6
               write (*,*) 'FACTOR: insufficient space in factor array ', knfax
               return
            ENDIF
            kfax(knfax) = ifac
! Store the other factor that goes with this one
            knfax = knfax + 1
            kfax(knfax) = inu / ifac
!WRITE (*,*) 'ARPDBG, factors ',knfax-1,' & ',knfax,' are ', kfax(knfax-1),' and ',kfax(knfax)
         ENDIF
!
      END DO
!
   20 CONTINUE      ! Label 20 is the exit point from the factor search loop.
!
   END SUBROUTINE factorise



   SUBROUTINE nemo_northcomms
!!----------------------------------------------------------------------
!!                     ***  ROUTINE  nemo_northcomms  ***
!! ** Purpose :   Setup for north fold exchanges with explicit
!!                point-to-point messaging
!!
!! ** Method :   Initialization of the northern neighbours lists.
!!----------------------------------------------------------------------
!!    1.0  ! 2011-10  (A. C. Coward, NOCS & J. Donners, PRACE)
!!    2.0  ! 2013-06 Setup avoiding MPI communication (I. Epicoco, S. Mocavero, CMCC)
!!----------------------------------------------------------------------
      INTEGER  ::   sxM, dxM, sxT, dxT, jn
      INTEGER  ::   njmppmax
!!----------------------------------------------------------------------
!
      njmppmax = MAXVAL( njmppt )
!
!initializes the north-fold communication variables
      isendto(:) = 0
      nsndto     = 0
!
!if I am a process in the north
      IF ( njmpp == njmppmax ) THEN
!sxM is the first point (in the global domain) needed to compute the
!north-fold for the current process
          sxM = jpiglo - nimppt(narea) - nlcit(narea) + 1
!dxM is the last point (in the global domain) needed to compute the
!north-fold for the current process
          dxM = jpiglo - nimppt(narea) + 2

!loop over the other north-fold processes to find the processes
!managing the points belonging to the sxT-dxT range
  
          DO jn = 1, jpni
!sxT is the first point (in the global domain) of the jn
!process
                sxT = nfiimpp(jn, jpnj)
!dxT is the last point (in the global domain) of the jn
!process
                dxT = nfiimpp(jn, jpnj) + nfilcit(jn, jpnj) - 1
                IF ((sxM .gt. sxT) .AND. (sxM .lt. dxT)) THEN
                   nsndto = nsndto + 1
                     isendto(nsndto) = jn
                ELSEIF ((sxM .le. sxT) .AND. (dxM .ge. dxT)) THEN
                   nsndto = nsndto + 1
                   isendto(nsndto) = jn
                ELSEIF ((dxM .lt. dxT) .AND. (sxT .lt. dxM)) THEN
                   nsndto = nsndto + 1
                   isendto(nsndto) = jn
                ENDIF
          END DO
          nfsloop = 1
          nfeloop = nlci
          DO jn = 2,jpni-1
           IF(nfipproc(jn,jpnj) .eq. (narea - 1)) THEN
              IF (nfipproc(jn - 1 ,jpnj) .eq. -1) THEN
                 nfsloop = nldi
              ENDIF
              IF (nfipproc(jn + 1,jpnj) .eq. -1) THEN
                 nfeloop = nlei
              ENDIF
           ENDIF
        END DO

      ENDIF
      l_north_nogather = .TRUE.
   END SUBROUTINE nemo_northcomms



!!======================================================================
END MODULE nemogcm

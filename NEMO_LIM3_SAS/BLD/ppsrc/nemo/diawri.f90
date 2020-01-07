MODULE diawri
!!======================================================================
!!                     ***  MODULE  diawri  ***
!! Ocean diagnostics :  write ocean output files
!!=====================================================================
!! History :  OPA  ! 1991-03  (M.-A. Foujols)  Original code
!!            4.0  ! 1991-11  (G. Madec)
!!                 ! 1992-06  (M. Imbard)  correction restart file
!!                 ! 1992-07  (M. Imbard)  split into diawri and rstwri
!!                 ! 1993-03  (M. Imbard)  suppress writibm
!!                 ! 1998-01  (C. Levy)  NETCDF format using ioipsl INTERFACE
!!                 ! 1999-02  (E. Guilyardi)  name of netCDF files + variables
!!            8.2  ! 2000-06  (M. Imbard)  Original code (diabort.F)
!!   NEMO     1.0  ! 2002-06  (A.Bozec, E. Durand)  Original code (diainit.F)
!!             -   ! 2002-09  (G. Madec)  F90: Free form and module
!!             -   ! 2002-12  (G. Madec)  merge of diabort and diainit, F90
!!                 ! 2005-11  (V. Garnier) Surface pressure gradient organization
!!            3.2  ! 2008-11  (B. Lemaire) creation from old diawri
!!----------------------------------------------------------------------

!!----------------------------------------------------------------------
!!   dia_wri       : create the standart output files
!!   dia_wri_state : create an output NetCDF file for a single instantaeous ocean state and forcing fields
!!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE zdf_oce         ! ocean vertical physics
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE sbc_ice         ! Surface boundary condition: ice fields
   USE sbcssr          ! restoring term toward SST/SSS climatology
   USE phycst          ! physical constants
   USE zdfmxl          ! mixed layer
   USE dianam          ! build name of file (routine)
   USE zdfddm          ! vertical  physics: double diffusion
   USE diahth          ! thermocline diagnostics
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! I/O manager
   USE iom
   USE ioipsl

   USE limwri

   USE lib_mpp         ! MPP library
   USE timing          ! preformance summary
   USE wrk_nemo        ! working array

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_wri                 ! routines called by step.F90
   PUBLIC   dia_wri_state
   PUBLIC   dia_wri_alloc           ! Called by nemogcm module

   INTEGER ::   nid_T, nz_T, nh_T, ndim_T, ndim_hT   ! grid_T file
   INTEGER ::   nid_U, nz_U, nh_U, ndim_U, ndim_hU   ! grid_U file
   INTEGER ::   nid_V, nz_V, nh_V, ndim_V, ndim_hV   ! grid_V file
   INTEGER ::   ndex(1)                              ! ???
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_hT, ndex_hU, ndex_hV

!! * Substitutions
!!----------------------------------------------------------------------
!!                   ***  vectopt_loop_substitute  ***
!!----------------------------------------------------------------------
!! ** purpose :   substitute the inner loop start/end indices with CPP macro
!!                allow unrolling of do-loop (useful with vector processors)
!!----------------------------------------------------------------------
!!----------------------------------------------------------------------
!! NEMO/OPA 3.7 , NEMO Consortium (2014)
!! $Id: vectopt_loop_substitute.h90 4990 2014-12-15 16:42:49Z timgraham $
!! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------




!!----------------------------------------------------------------------
!! NEMO/OPA 3.3 , NEMO Consortium (2010)
!! $Id: diawri.F90 7761 2017-03-06 17:58:35Z clem $
!! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION dia_wri_alloc()
!!----------------------------------------------------------------------
      INTEGER :: ierr
!!----------------------------------------------------------------------
!
      ALLOCATE( ndex_hT(jpi*jpj), ndex_hU(jpi*jpj), ndex_hV(jpi*jpj), STAT=dia_wri_alloc )
      IF( lk_mpp )   CALL mpp_sum( dia_wri_alloc )
!
  END FUNCTION dia_wri_alloc

!!----------------------------------------------------------------------
!!   Default option                                   NetCDF output file
!!----------------------------------------------------------------------

!!----------------------------------------------------------------------
!!   'key_iomput'                                        use IOM library
!!----------------------------------------------------------------------

   SUBROUTINE dia_wri( kt )
!!---------------------------------------------------------------------
!!                  ***  ROUTINE dia_wri  ***
!!
!! ** Purpose :   Standard output of opa: dynamics and tracer fields
!!      NETCDF format is used by default
!!      Standalone surface scheme
!!
!! ** Method  :  use iom_put
!!----------------------------------------------------------------------
!!
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
!!----------------------------------------------------------------------
!
! Output the initial state and forcings
      IF( ninist == 1 ) THEN
         CALL dia_wri_state( 'output.init', kt )
         ninist = 0
      ENDIF
!
   END SUBROUTINE dia_wri



   SUBROUTINE dia_wri_state( cdfile_name, kt )
!!---------------------------------------------------------------------
!!                 ***  ROUTINE dia_wri_state  ***
!!
!! ** Purpose :   create a NetCDF file named cdfile_name which contains
!!      the instantaneous ocean state and forcing fields.
!!        Used to find errors in the initial state or save the last
!!      ocean state in case of abnormal end of a simulation
!!
!! ** Method  :   NetCDF files using ioipsl
!!      File 'output.init.nc'  is created if ninist = 1 (namelist)
!!      File 'output.abort.nc' is created in case of abnormal job end
!!----------------------------------------------------------------------
      CHARACTER (len=* ), INTENT( in ) ::   cdfile_name      ! name of the file created
      INTEGER           , INTENT( in ) ::   kt               ! ocean time-step index
!!
      CHARACTER (len=32) :: clname
      CHARACTER (len=40) :: clop
      INTEGER  ::   id_i , nz_i, nh_i       
      INTEGER, DIMENSION(1) ::   idex             ! local workspace
      REAL(wp) ::   zsto, zout, zmax, zjulian
!!----------------------------------------------------------------------
!
      IF( nn_timing == 1 )   CALL timing_start('dia_wri_state')

! 0. Initialisation
! -----------------

! Define name, frequency of output and means
      clname = cdfile_name
      IF( .NOT. Agrif_Root() ) clname = TRIM(Agrif_CFixed())//'_'//TRIM(clname)
      zsto = rdt
      clop = "inst(x)"           ! no use of the mask value (require less cpu time)
      zout = rdt
      zmax = ( nitend - nit000 + 1 ) * rdt

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'dia_wri_state : single instantaneous ocean state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~   and forcing fields file created '
      IF(lwp) WRITE(numout,*) '                and named :', clname, '.nc'


! 1. Define NETCDF files and fields at beginning of first time step
! -----------------------------------------------------------------

! Compute julian date from starting date of the run
      CALL ymds2ju( nyear, nmonth, nday, rdt, zjulian )         ! time axis
      zjulian = zjulian - adatrj   !   set calendar origin to the beginning of the experiment
      CALL histbeg( clname, jpi, glamt, jpj, gphit,   &
          1, jpi, 1, jpj, nit000-1, zjulian, rdt, nh_i, id_i, domain_id=nidom, snc4chunks=snc4set ) ! Horizontal grid : glamt and gphit
      CALL histvert( id_i, "deptht", "Vertical T levels",   &    ! Vertical grid : gdept
          "m", jpk, gdept_1d, nz_i, "down")

! Declare all the output fields as NetCDF variables

      CALL histdef( id_i, "sowaflup", "Net Upward Water Flux" , "Kg/m2/S",   &   ! net freshwater
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( id_i, "sohefldo", "Net Downward Heat Flux", "W/m2"   ,   &   ! net heat flux
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( id_i, "soshfldo", "Shortwave Radiation"   , "W/m2"   ,   &   ! solar flux
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( id_i, "soicecov", "Ice fraction"          , "[0,1]"  ,   &   ! fr_i
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( id_i, "sozotaux", "Zonal Wind Stress"     , "N/m2"   ,   &   ! i-wind stress
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( id_i, "sometauy", "Meridional Wind Stress", "N/m2"   ,   &   ! j-wind stress
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )


      CALL lim_wri_state( kt, id_i, nh_i )


! 2. Start writing data
! ---------------------
! idex(1) est utilise ssi l'avant dernier argument est diffferent de
! la taille du tableau en sortie. Dans ce cas , l'avant dernier argument
! donne le nombre d'elements, et idex la liste des indices a sortir
      idex(1) = 1   ! init to avoid compil warning

! Write all fields on T grid
      CALL histwrite( id_i, "sowaflup", kt, emp              , jpi*jpj    , idex )    ! freshwater budget
      CALL histwrite( id_i, "sohefldo", kt, qsr + qns        , jpi*jpj    , idex )    ! total heat flux
      CALL histwrite( id_i, "soshfldo", kt, qsr              , jpi*jpj    , idex )    ! solar heat flux
      CALL histwrite( id_i, "soicecov", kt, fr_i             , jpi*jpj    , idex )    ! ice fraction
      CALL histwrite( id_i, "sozotaux", kt, utau             , jpi*jpj    , idex )    ! i-wind stress
      CALL histwrite( id_i, "sometauy", kt, vtau             , jpi*jpj    , idex )    ! j-wind stress

! 3. Close the file
! -----------------
      CALL histclo( id_i )

       
      IF( nn_timing == 1 )   CALL timing_stop('dia_wri_state')
!

   END SUBROUTINE dia_wri_state
!!======================================================================
END MODULE diawri

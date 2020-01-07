MODULE diaar5
   !!======================================================================
   !!                       ***  MODULE  diaar5  ***
   !! AR5 diagnostics
   !!======================================================================
   !! History :  3.2  !  2009-11  (S. Masson)  Original code
   !!            3.3  !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase + merge TRC-TRA
   !!----------------------------------------------------------------------
#if defined key_diaar5   || defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_diaar5'  :                           activate ar5 diagnotics
   !!----------------------------------------------------------------------
   !!   dia_ar5       : AR5 diagnostics
   !!   dia_ar5_init  : initialisation of AR5 diagnostics
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers 
   USE dom_oce        ! ocean space and time domain
   USE eosbn2         ! equation of state                (eos_bn2 routine)
   USE lib_mpp        ! distribued memory computing library
   USE iom            ! I/O manager library
   USE timing         ! preformance summary
   USE wrk_nemo       ! working arrays
   USE fldread        ! type FLD_N
   USE phycst         ! physical constant
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_ar5        ! routine called in step.F90 module
   PUBLIC   dia_ar5_init   ! routine called in opa.F90 module
   PUBLIC   dia_ar5_alloc  ! routine called in nemogcm.F90 module

   LOGICAL, PUBLIC, PARAMETER :: lk_diaar5 = .TRUE.   ! coupled flag

   REAL(wp)                         ::   vol0         ! ocean volume (interior domain)
   REAL(wp)                         ::   area_tot     ! total ocean surface (interior domain)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:  ) ::   area         ! cell surface (interior domain)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:  ) ::   thick0       ! ocean thickness (interior domain)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sn0          ! initial salinity
      
   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: diaar5.F90 6664 2016-06-03 13:42:46Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   FUNCTION dia_ar5_alloc()
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_ar5_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: dia_ar5_alloc
      !!----------------------------------------------------------------------
      !
      ALLOCATE( area(jpi,jpj), thick0(jpi,jpj) , sn0(jpi,jpj,jpk) , STAT=dia_ar5_alloc )
      !
      IF( lk_mpp             )   CALL mpp_sum ( dia_ar5_alloc )
      IF( dia_ar5_alloc /= 0 )   CALL ctl_warn('dia_ar5_alloc: failed to allocate arrays')
      !
   END FUNCTION dia_ar5_alloc


   SUBROUTINE dia_ar5( kt )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_ar5  ***
      !!
      !! ** Purpose :   compute and output some AR5 diagnostics
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk                      ! dummy loop arguments
      REAL(wp) ::   zvolssh, zvol, zssh_steric, zztmp, zarho, ztemp, zsal, zmass
      !
      REAL(wp), POINTER, DIMENSION(:,:)     :: zarea_ssh , zbotpres       ! 2D workspace 
      REAL(wp), POINTER, DIMENSION(:,:,:)   :: zrhd , zrhop               ! 3D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:,:) :: ztsn                       ! 4D workspace
      !!--------------------------------------------------------------------
      IF( nn_timing == 1 )   CALL timing_start('dia_ar5')
 
      CALL wrk_alloc( jpi , jpj              , zarea_ssh , zbotpres )
      CALL wrk_alloc( jpi , jpj , jpk        , zrhd      , zrhop    )
      CALL wrk_alloc( jpi , jpj , jpk , jpts , ztsn                 )

      zarea_ssh(:,:) = area(:,:) * sshn(:,:)

      !                                         ! total volume of liquid seawater
      zvolssh = SUM( zarea_ssh(:,:) ) 
      IF( lk_mpp )   CALL mpp_sum( zvolssh )
      zvol = vol0 + zvolssh
      
      CALL iom_put( 'voltot', zvol               )
      CALL iom_put( 'sshtot', zvolssh / area_tot )

      !                     
!CT      ztsn(:,:,:,jp_tem) = tsn(:,:,:,jp_tem)                    ! thermosteric ssh
!CT      ztsn(:,:,:,jp_sal) = sn0(:,:,:)
!CT      CALL eos( ztsn, zrhd, fsdept_n(:,:,:) )                       ! now in situ density using initial salinity
!CT      !
!CT      zbotpres(:,:) = 0._wp                        ! no atmospheric surface pressure, levitating sea-ice
!CT      DO jk = 1, jpkm1
!CT         zbotpres(:,:) = zbotpres(:,:) + fse3t(:,:,jk) * zrhd(:,:,jk)
!CT      END DO
!CT      IF( .NOT.lk_vvl ) THEN
!CT         IF ( ln_isfcav ) THEN
!CT            DO ji=1,jpi
!CT               DO jj=1,jpj
!CT                  zbotpres(ji,jj) = zbotpres(ji,jj) + sshn(ji,jj) * zrhd(ji,jj,mikt(ji,jj)) + riceload(ji,jj)
!CT               END DO
!CT            END DO
!CT         ELSE
!CT            zbotpres(:,:) = zbotpres(:,:) + sshn(:,:) * zrhd(:,:,1)
!CT         END IF
!CT      END IF
!CT      !                                         
!CT      zarho = SUM( area(:,:) * zbotpres(:,:) ) 
!CT      IF( lk_mpp )   CALL mpp_sum( zarho )
!CT      zssh_steric = - zarho / area_tot
!CT      CALL iom_put( 'sshthster', zssh_steric )
      
      !                                         ! steric sea surface height
      CALL eos( tsn, zrhd, zrhop, fsdept_n(:,:,:) )                 ! now in situ and potential density
      zrhop(:,:,jpk) = 0._wp
      CALL iom_put( 'rhop', zrhop )
      !
      zbotpres(:,:) = 0._wp                        ! no atmospheric surface pressure, levitating sea-ice
      DO jk = 1, jpkm1
         zbotpres(:,:) = zbotpres(:,:) + fse3t(:,:,jk) * zrhd(:,:,jk)
      END DO
      IF( .NOT.lk_vvl ) THEN
         IF ( ln_isfcav ) THEN
            DO ji=1,jpi
               DO jj=1,jpj
                  zbotpres(ji,jj) = zbotpres(ji,jj) + sshn(ji,jj) * zrhd(ji,jj,mikt(ji,jj)) + riceload(ji,jj)
               END DO
            END DO
         ELSE
            zbotpres(:,:) = zbotpres(:,:) + sshn(:,:) * zrhd(:,:,1)
         END IF
      END IF
      !    
!CT      zarho = SUM( area(:,:) * zbotpres(:,:) ) 
!CT      IF( lk_mpp )   CALL mpp_sum( zarho )
!CT      zssh_steric = - zarho / area_tot
!CT      CALL iom_put( 'sshsteric', zssh_steric )
      
      !                                         ! ocean bottom pressure
      zztmp = rau0 * grav * 1.e-4_wp               ! recover pressure from pressure anomaly and cover to dbar = 1.e4 Pa
      zbotpres(:,:) = zztmp * ( zbotpres(:,:) + sshn(:,:) + thick0(:,:) )
      CALL iom_put( 'botpres', zbotpres )

      !                                         ! Mean density anomalie, temperature and salinity
!CT      ztemp = 0._wp
!CT      zsal  = 0._wp
!CT      DO jk = 1, jpkm1
!CT         DO jj = 1, jpj
!CT            DO ji = 1, jpi
!CT               zztmp = area(ji,jj) * fse3t(ji,jj,jk)
!CT               ztemp = ztemp + zztmp * tsn(ji,jj,jk,jp_tem)
!CT               zsal  = zsal  + zztmp * tsn(ji,jj,jk,jp_sal)
!CT            END DO
!CT         END DO
!CT      END DO
!CT      IF( .NOT.lk_vvl ) THEN
!CT         IF ( ln_isfcav ) THEN
!CT            DO ji=1,jpi
!CT               DO jj=1,jpj
!CT                  ztemp = ztemp + zarea_ssh(ji,jj) * tsn(ji,jj,mikt(ji,jj),jp_tem) 
!CT                  zsal  = zsal  + zarea_ssh(ji,jj) * tsn(ji,jj,mikt(ji,jj),jp_sal) 
!CT               END DO
!CT            END DO
!CT         ELSE
!CT            ztemp = ztemp + SUM( zarea_ssh(:,:) * tsn(:,:,1,jp_tem) )
!CT            zsal  = zsal  + SUM( zarea_ssh(:,:) * tsn(:,:,1,jp_sal) )
!CT         END IF
!CT      ENDIF
!CT      IF( lk_mpp ) THEN  
!CT         CALL mpp_sum( ztemp )
!CT         CALL mpp_sum( zsal  )
!CT      END IF
!CT      !
!CT      zmass = rau0 * ( zarho + zvol )                 ! total mass of liquid seawater
!CT      ztemp = ztemp / zvol                            ! potential temperature in liquid seawater
!CT      zsal  = zsal  / zvol                            ! Salinity of liquid seawater
!CT      !
!CT      CALL iom_put( 'masstot', zmass )
!CT      CALL iom_put( 'temptot', ztemp )
!CT      CALL iom_put( 'saltot' , zsal  )
      !
      CALL wrk_dealloc( jpi , jpj              , zarea_ssh , zbotpres )
      CALL wrk_dealloc( jpi , jpj , jpk        , zrhd      , zrhop    )
      CALL wrk_dealloc( jpi , jpj , jpk , jpts , ztsn                 )
      !
      IF( nn_timing == 1 )   CALL timing_stop('dia_ar5')
      !
   END SUBROUTINE dia_ar5


   SUBROUTINE dia_ar5_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dia_ar5_init  ***
      !!                   
      !! ** Purpose :   initialization for AR5 diagnostic computation
      !!----------------------------------------------------------------------
      INTEGER  ::   inum
      INTEGER  ::   ik
      INTEGER  ::   ji, jj, jk  ! dummy loop indices
      REAL(wp) ::   zztmp  
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::   zsaldta   ! Jan/Dec levitus salinity
      !
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('dia_ar5_init')
      !
!CT      CALL wrk_alloc( jpi, jpj, jpk, 2, zsaldta )
      !                                      ! allocate dia_ar5 arrays
      IF( dia_ar5_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'dia_ar5_init : unable to allocate arrays' )

!CT CREG To exclude the BDY in the volume control calculation
!CT CREG use the bmask which is set to zero along the BDY in bdyini.F90 {
      area(:,:) = e1t(:,:) * e2t(:,:) * tmask_i(:,:) * bmask(:,:)
!CT }

      area_tot = SUM( area(:,:) )   ;   IF( lk_mpp )   CALL mpp_sum( area_tot )

      vol0        = 0._wp
      thick0(:,:) = 0._wp
      DO jk = 1, jpkm1
         vol0        = vol0        + SUM( area (:,:) * tmask(:,:,jk) * e3t_0(:,:,jk) )
         thick0(:,:) = thick0(:,:) +    tmask_i(:,:) * tmask(:,:,jk) * e3t_0(:,:,jk)
      END DO
      IF( lk_mpp )   CALL mpp_sum( vol0 )

!CT      CALL iom_open ( 'sali_ref_clim_monthly', inum )
!CT      CALL iom_get  ( inum, jpdom_data, 'vosaline' , zsaldta(:,:,:,1), 1  )
!CT      CALL iom_get  ( inum, jpdom_data, 'vosaline' , zsaldta(:,:,:,2), 12 )
!CT      CALL iom_close( inum )
!CT
!CT      sn0(:,:,:) = 0.5_wp * ( zsaldta(:,:,:,1) + zsaldta(:,:,:,2) )        
!CT      sn0(:,:,:) = sn0(:,:,:) * tmask(:,:,:)
!CT      IF( ln_zps ) THEN               ! z-coord. partial steps
!CT         DO jj = 1, jpj               ! interpolation of salinity at the last ocean level (i.e. the partial step)
!CT            DO ji = 1, jpi
!CT               ik = mbkt(ji,jj)
!CT               IF( ik > 1 ) THEN
!CT                  zztmp = ( gdept_1d(ik) - gdept_0(ji,jj,ik) ) / ( gdept_1d(ik) - gdept_1d(ik-1) )
!CT                  sn0(ji,jj,ik) = ( 1._wp - zztmp ) * sn0(ji,jj,ik) + zztmp * sn0(ji,jj,ik-1)
!CT               ENDIF
!CT            END DO
!CT         END DO
!CT      ENDIF
!CT      !
!CT      CALL wrk_dealloc( jpi, jpj, jpk, 2, zsaldta )
      !
      IF( nn_timing == 1 )   CALL timing_stop('dia_ar5_init')
      !
   END SUBROUTINE dia_ar5_init

#else
   !!----------------------------------------------------------------------
   !!   Default option :                                         NO diaar5
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER :: lk_diaar5 = .FALSE.   ! coupled flag
CONTAINS
   SUBROUTINE dia_ar5_init    ! Dummy routine
   END SUBROUTINE dia_ar5_init
   SUBROUTINE dia_ar5( kt )   ! Empty routine
      INTEGER ::   kt
      WRITE(*,*) 'dia_ar5: You should not have seen this print! error?', kt
   END SUBROUTINE dia_ar5
#endif

   !!======================================================================
END MODULE diaar5

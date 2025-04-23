module module_add_emiss_burn
!RAR: significantly modified for the new BB emissions
CONTAINS
  subroutine add_emis_burn(id,dtstep,ktau,dz8w,rho_phy,rel_hum,      &
                           chem,julday,gmt,xlat,xlong,               &
                           luf_igbp,nlcat,                           &
                           lu_fire1,peak_hr,                          &
                           time_int,ebu,hwp,                         &   ! RAR
                           r_q,fhist,aod3d,                          &
                           poa_vol_params,                           &
 !                          nwfa,nifa,                                &
                           scale_factor_NO2,scale_factor_NO,         &
                           new_NO2_EF,new_NO_EF,                     &
                           rainc,rainnc, swdown,                     &
                           qv,                                       &
                           config_flags,                             &
                           ids,ide, jds,jde, kds,kde,                &
                           ims,ime, jms,jme, kms,kme,                &
                           its,ite, jts,jte, kts,kte                 )
! xlat, xlong, 

   USE module_configure, only: grid_config_rec_type
   USE module_state_description
   IMPLICIT NONE

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: ktau, id, julday, nlcat,        &
                                  ids,ide, jds,jde, kds,kde,      &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                                   chem
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT) :: qv

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_ebu ),                 &
         INTENT(IN) ::                                   ebu
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN) :: hwp
   REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN) ::  xlat,xlong, rainc,rainnc,swdown, lu_fire1, peak_hr
   REAL, DIMENSION(ims:ime,jms:jme), INTENT(OUT)    ::  r_q    ! RAR:
   REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT)  ::  fhist  ! RAR:
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(OUT) ::  aod3d    ! RAR:

   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: dz8w,rho_phy,rel_hum
   REAL, DIMENSION(ims:ime,1:nlcat,jms:jme), INTENT(IN) :: luf_igbp

   REAL, DIMENSION (5), INTENT(IN)   ::      poa_vol_params         ! RAR:/ [JLS]

!   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
!          OPTIONAL, INTENT(INOUT   ) ::                 nwfa,nifa                   ! RAR:

    REAL, INTENT(IN   ) ::    dtstep,gmt
    REAL(KIND=8), INTENT(IN) ::  time_int       ! RAR: time in seconds since start of simulation
!
! Lindsey added intents
!

  REAL, DIMENSION(ims:ime, jms:jme), INTENT (OUT) :: scale_factor_NO2, scale_factor_NO,new_NO2_EF,new_NO_EF

!
!
    integer :: i,j,k,n,m
    real :: conv_rho, conv, ext2, dm_smoke, daero_num_wfa, daero_num_ifa !, lu_sum1_5, lu_sum12_14
    real :: alpha_i, alpha_j,FS_RH,FL_RH,aerext,F_SS
    !real :: ebumax
!    CHARACTER (LEN=80) :: message

    INTEGER, PARAMETER :: kfire_max=35    ! max vertical level for BB plume rise
    ! Diameters and standard deviations for emissions
    ! the diameters are the volume (mass) geometric mean diameters, following MADE_SORGAM
    REAL, PARAMETER :: dgvem_i= 0.08E-6 !0.03E-6 ! [ m ]
    REAL, PARAMETER :: sgem_i = 1.8     !1.7

    ! *** Accumulation mode:
    REAL, PARAMETER :: dgvem_j= 0.3E-6 ! [ m ]
    REAL, PARAMETER :: sgem_j = 2.0

    ! *** Coarse mode
    REAL, PARAMETER :: dgvem_c= 6.0E-6 ! [ m ] 
    REAL, PARAMETER :: sgem_c=  2.2
    REAL, PARAMETER :: pic= 3.14159

    ! RAR: factors for getting number emissions rate from mass emissions rate following made_sorgam
    REAL, PARAMETER :: fact_numn= 1.e-9*6.0/pic*exp(4.5*log(sgem_i)**2)/dgvem_i**3       ! Aitken mode
    REAL, PARAMETER :: fact_numa= 1.e-9*6.0/pic*exp(4.5*log(sgem_j)**2)/dgvem_j**3       ! accumulation mode
    REAL, PARAMETER :: fact_numc= 1.e-9*6.0/pic*exp(4.5*log(sgem_c)**2)/dgvem_c**3       ! coarse mode

    REAL, PARAMETER :: dens_oc_aer=1.4e3, dens_ec_aer=1.7e3  ! kg/m3
!    REAL, PARAMETER :: rinti=2.1813936e-8, cx=2.184936* 3600, timeq_max=3600.*24.   ! constants for the diurnal cycle calculations
    REAL, PARAMETER :: ax1=531., cx1=7800.  ! For cropland, urban and small fires
!    REAL, PARAMETER :: rinti=2.1813936e-8, ax2=3200., const2=100., coef2=10.6712963e-4, cx=2.184936* 3600, timeq_max=3600.*24.
    REAL, PARAMETER :: rinti=2.1813936e-8, ax2=3400., const2=130., coef2=10.6712963e-4, cx2=7200., timeq_max=3600.*24.    ! New parameters
    REAL, PARAMETER :: sc_me= 4.0, ab_me=0.5     ! m2/g, scattering and absorption efficiency for smoke


    REAL, PARAMETER :: pirs    =3.14159265358979324
    REAL, PARAMETER :: f6dpi   = 6.0/pirs
    REAL, PARAMETER :: f6dpim9 = 1.0E-9*f6dpi
    REAL, PARAMETER :: rhoorg  = 1.4E3 
    REAL, PARAMETER :: rhoanth = 2.2E3
    REAL, PARAMETER :: rhono3=1.8E3
    REAL, PARAMETER :: rhoso4=1.8E3
    REAL, PARAMETER :: orgfac  = f6dpim9/rhoorg
    REAL, PARAMETER :: anthfac = f6dpim9/rhoanth
    REAL, PARAMETER :: no3fac=f6dpim9/rhono3
    REAL, PARAMETER :: so4fac=f6dpim9/rhoso4
    REAL, PARAMETER :: rayleigh_scat_param = 8.
    REAL, PARAMETER :: sm_ammsulf_fact = 2.2 
    REAL, PARAMETER :: lg_ammsulf_fact = 4.8
    REAL, PARAMETER :: sm_ammnitr_fact = 2.4
    REAL, PARAMETER :: lg_ammnitr_fact = 5.1
    REAL, PARAMETER :: sm_org_fact     = 2.8
    REAL, PARAMETER :: lg_org_fact     = 6.1
    REAL, PARAMETER :: ec_fact         = 10.0
    REAL, PARAMETER :: no2_fact        = 0.33
    REAL, PARAMETER :: dcycle_adjust   = 0.75 ! The computed diurnal cycle gives weights that sum to ~32 over a 24 hr day so we reduce by 25% and apply to gases only !!!  JLS
!   Parameters used for the wfa and ifa in mp physics per Trude E. (NCAR)
!   Water friendly: radius: 0.04 micron, standard deviation: 1.8, kappa (for hygroscopic growth): 0.2, real index of refraction: 1.53, imaginary index of refraction: 1e-7
!   Ice friendly: radius: 0.4 micron, standard deviation: 1.8, kappa : 0.04, real index of refraction: 1.56, imaginary index of refraction: 3e-3

    !    real, parameter :: cx        =  2.184936 * 3600., rinti     =  2.1813936e-8    , ax        =  2000.6038
    ! bx_bburn  =  20.041288 * 3600.,  RAR: this depends on the vegetation class, location (local time) etc.
    REAL :: timeq, dt1,dt2,dtm         ! For BB emis. diurnal cycle calculation

    INTEGER, SAVE ::  icall

!
! Lindsey Local
!

   ! L. Anderson start declaring variables
   ! Declare constants
    REAL, PARAMETER :: NOx_as_NO_to_NO2_CONSTANT = 0.77
    REAL, PARAMETER :: NOx_as_NO_to_NO_CONSTANT = 0.50 ! /scratch1/BMC/wrf-chem/Jordan/others/forLindsey/convert_AeM_to_RADM_WRF_FIM.f90
    REAL, PARAMETER :: scale_NO2CO_HWP_slope = 4.
    REAL, PARAMETER :: min_NO2CO_ER = 0.0005
   ! Add other constants as needed

   ! Declare variables
   REAL :: min_lat, max_lat, min_lon, max_lon

   ! Notes about choices:
   ! western US: 31.211480,-125.419922,49.531336,-104.001786
   ! http://bboxfinder.com/#31.211480,-125.419922,49.531336,-104.001786
   ! Northern CA, OR, WA, Montana 38.869423,-125.456829,49.287419,-104.004478
   ! http://bboxfinder.com/#38.869423,-125.456829,49.287419,-104.004478
   
   ! Declare loop indices
   INTEGER :: loop_index

   ! Declare arrays
   REAL :: NO2CO_ER, estimated_CO_EF, estimated_NOx_as_NO_EF, estimated_NO2_EF,lufrac ! hwp previously defined (last two intent out)
   REAL :: estimated_NO_EF, NOxCO_ER 
   REAL, DIMENSION(1:nlcat) :: CO_EF, NOx_as_NO_EF

   ! Set the bounding box values (replace with your actual values)
   min_lat = 31.0
   max_lat = 49.5
   min_lon = -125.4 ! I think lats and lons are currently -180 to 180
   max_lon = -104.0
   CO_EF(1:nlcat)        = (/112.68, 112.68, 121.35, 112.68, 90.929, 69.177, 69.177, 69.177, 69.177, 69.177, 0.0, 69.177, 0.0, 69.177, 0.0, 69.177, 0.0, 121.35, 95.264, 0.0, 0.0/)
   NOX_as_NO_EF(1:nlcat) = (/1.315,  1.315,  1.177,  1.315,  1.9,    2.486,  2.486,  2.486,  2.486,  2.486,  0.0, 2.486,  0.0, 2.486,  0.0, 2.486,  0.0, 1.177,  1.832,  0.0, 0.0/)
   ! zeros at 11, 13, 15, 17

    timeq= gmt*3600. + REAL(time_int,4)
    timeq= mod(timeq,timeq_max)

          do j = jts, jte
          do i = its, ite
          scale_factor_NO2(i,j) = 1.0
          scale_factor_NO(i,j) = 1.0
          enddo
          enddo

! Main loops to add BB emissions
    do j=jts,jte
       do i=its,ite
          if( luf_igbp(i,17,j)>0.99 .OR. ebu(i,1,j,p_ebu_no) < 1.e-6) cycle       ! no BB emissions or water pixels

          ! RAR: the decrease in the BB emissions after >18 hrs of forecast, the decrease occurs at night. The decrease occurs at night.
          IF (time_int>64800. .AND. swdown(i,j)<.1 .AND. fhist(i,j)>.75 ) THEN
              fhist(i,j)= 0.75
          ENDIF

          IF (time_int>129600. .AND. swdown(i,j)<.1 .AND. fhist(i,j)>.5 ) THEN      ! After 36 hr forecast 
              fhist(i,j)= 0.5
          ENDIF

          IF ( (rainc(i,j) + rainnc(i,j))>=10. .AND. fhist(i,j)>.3 ) THEN    ! If it rains more than 1cm, then the BB emissions are reduced
              fhist(i,j)= 0.3
          ENDIF

! RAR: Grasslands (29% of ther western HRRR CONUS domain) probably also need to be added below, check this later
! RAR: In the HRRR CONUS domain (western part) crop 11%, 2% cropland/natural vegetation and 0.4% urban of pixels
!.OR. lu_index(i,j)==14)  then ! Croplands(12), Urban and Built-Up(13), cropland/natural vegetation (14) mosaic in MODI-RUC vegetation classes
! Peak hours for the fire activity depending on the latitude
!                   if (xlong(i,j)<-130.) then  max_ti= 24.041288* 3600.    ! peak at 24 UTC, fires in  Alaska
!                   elseif (xlong(i,j)<-100.) then   max_ti= 22.041288* 3600.    ! peak at 22 UTC, fires in the western US
!                   elseif (xlong(i,j)<-70.) then   ! peak at 20 UTC, fires in the eastern US,   max_ti= 20.041288* 3600.
!                   else   max_ti= 18.041288* 3600.
!                   endif

           IF ( lu_fire1(i,j)>0.9 ) then    !Ag, urban fires, bare land etc.
           !   these fires will have exponentially decreasing diurnal cycle, these fires decrease 55% in 2 hours, end in 5 hours
               r_q(i,j) = rinti* ax1 * exp(- (time_int**2)/(cx1**2) )
           ELSE
               ! RAR: Gaussian profile for wildfires
               dt1= abs(timeq - peak_hr(i,j))
               dt2= timeq_max - peak_hr(i,j) + timeq   ! peak hour is always <86400.
               dtm= MIN(dt1,dt2)
               r_q(i,j) = rinti*( ax2 * exp(- dtm**2/(2.*cx2**2) ) + const2 - coef2*timeq )
           ENDIF

           r_q(i,j) = fhist(i,j)* max(0.,r_q(i,j)*timeq_max)

           IF (swdown(i,j)<.1) THEN
               r_q(i,j)= MIN(0.5,r_q(i,j))   ! lower BB emissions at night
           ENDIF

           IF (.NOT. config_flags%bb_dcycle) THEN
               r_q(i,j)= 1. + 1./3. !fhist(i,j)    ! no diurnal cycle
           END IF
! INSERT LINDSEY EDITS HERE

!             if ( ebu(i,1,j,p_ebu_in_no2) .gt. 0.0 .and. hwp(i,j) .eq. hwp(i,j) .and. hwp(i,j) .gt. 0.0 &
!                 .and. config_flags%mod_no2_hwp .eq. 1 .and. xlat(i,j) > min_lat .and. xlat(i,j) < max_lat .and. &
!                  xlong(i,j) > min_lon .and. xlong(i,j) < max_lon) then ! if inside bounding box for region of interest:
             if ( ebu(i,1,j,p_ebu_in_no2) .gt. 0.0 .and. hwp(i,j) .gt. 0.0 &
                 .and. config_flags%mod_no2_hwp .eq. 1 .and. xlat(i,j) > min_lat .and. xlat(i,j) < max_lat .and. &
                  xlong(i,j) > min_lon .and. xlong(i,j) < max_lon) then ! if inside bounding box for region of interest:
                 
!                 NO2CO_ER = scale_NO2CO_HWP_slope*2.234e-6 * hwp(i,j) + 0.000127 + min_NO2CO_ER
                 NO2CO_ER = scale_NO2CO_HWP_slope*2.234e-6 * max(hwp(i,j),5.0) + 0.000127 + min_NO2CO_ER
                 NOxCO_ER = 1.2/1.0 * NO2CO_ER
                 if (hwp(i,j) .ne. hwp(i,j) ) then
                    write(*,*) 'JLS, hwp is NaN at i,j = ',xlat(i,j), xlong(i,j)
                 endif
                 estimated_CO_EF = 0.0
                 estimated_NOx_as_NO_EF = 0.0
                 do loop_index=1,19 ! only to 19, last two lutypes have emis rate = 0. 
                 ! calculate EFs based on MODIS land types
                     lufrac = luf_igbp(i,loop_index,j)
                     if (lufrac .gt. 0.0) then
                     estimated_CO_EF = estimated_CO_EF +  luf_igbp(i,loop_index,j) * CO_EF(loop_index)                
                     estimated_NOx_as_NO_EF = estimated_NOx_as_NO_EF +  luf_igbp(i,loop_index,j) * NOx_as_NO_EF(loop_index) 
                     endif
                 end do
                                                 
                 estimated_NO2_EF = NOx_as_NO_to_NO2_CONSTANT * estimated_NOx_as_NO_EF
                 estimated_NO_EF = NOx_as_NO_to_NO_CONSTANT * estimated_NOx_as_NO_EF ! LA 11292023
                 ! NOx_as_NO_to_NO2_CONSTANT = the scaling factor used to estimate NO2 based on the NOx_as_NO_ER
                 new_NO2_EF(i,j) = NOxCO_ER * 0.5 * 1.6424 * estimated_CO_EF ! (MW(NO2)/(MW(CO)) = 1.6424 
                 new_NO_EF(i,j) = NOxCO_ER * 0.5 * 1.0714 * estimated_CO_EF ! (MW(NO)/(MW(CO)) = 1.0714 ! LA 11292023
                 ! new_NO2_EF(i,j) = NO2CO_ER(i,j) * 1.0714 * estimated_CO_EF(i,j) ! (MW(NO)/(MW(CO)) = 1.0714
                 ! new_NO2_EF = NO2CO_ER * (MW(NO since NOx as NO)/(MW(CO)) * CO_EF
                 ! new_NO2_EF = NO2CO_ER * 1.6424 * estimated_CO_EF ! (MW(NO2)/(MW(CO)) = 1.6424 
                 if ( estimated_NO2_EF .gt. 0.0 ) then !.and. estimated_NO2_EF .eq. estimated_NO2_EF .and. &
!                      new_NO2_EF(i,j) .eq. new_NO2_EF(i,j) ) then
                      scale_factor_NO2(i,j) = min(max(new_NO2_EF(i,j)/estimated_NO2_EF,0.0),10.0)
                 else
                    write(*,*) 'JLS, i,j = ',i,j
                    write(*,*) 'JLS, xlat,xlong = ',xlat(i,j), xlong(i,j)
                    write(*,*) 'JLS, estimated_NO2_EF(i,j) = ', estimated_NO2_EF
                    write(*,*) 'JLS, new_NO2_EF(i,j) = ', new_NO2_EF(i,j)
                    do loop_index =1,19
                        write(*,*) 'JLS, index, luf_igbp(i,index,j) = ', luf_igbp(i,loop_index,j)
                    enddo
                    scale_factor_NO2(i,j) = 1.0
                 endif
                 if ( estimated_NO_EF .gt. 0.0 ) then !.and. estimated_NO_EF .eq. estimated_NO_EF .and. &
!                      new_NO_EF(i,j) .eq. new_NO_EF(i,j) ) then
                    scale_factor_NO(i,j) = min(max(new_NO_EF(i,j)/estimated_NO_EF,0.0),10.0)                ! LA 11292023 end
                 else
!                    write(*,*) 'JLS, estimated_NO_EF(i,j) = ', estimated_NO_EF
!                    write(*,*) 'JLS, new_NO_EF(i,j) = ', new_NO_EF(i,j)
                    scale_factor_NO(i,j) = 1.0
                 endif
             endif ! Lindsey scale factor

           do k=kts,kfire_max
! Aero conversion       
              conv= r_q(i,j)*dtstep/(rho_phy(i,k,j) *  dz8w(i,k,j))
! Gas conversion
              conv_rho=dcycle_adjust * r_q(i,j)*4.828e-4/rho_phy(i,k,j)*dtstep/60./dz8w(i,k,j)
! Add gases
              chem(i,k,j,p_no)     = chem(i,k,j,p_no)     + scale_factor_NO(i,j)*conv_rho*ebu(i,k,j,p_ebu_no)
              chem(i,k,j,p_co)     = chem(i,k,j,p_co)     + conv_rho*ebu(i,k,j,p_ebu_co)
              chem(i,k,j,p_no2)    = chem(i,k,j,p_no2)    + scale_factor_NO2(i,j)*conv_rho*ebu(i,k,j,p_ebu_no2) !
              chem(i,k,j,p_no)     = chem(i,k,j,p_no)     + conv_rho*ebu(i,k,j,p_ebu_hono)!
              chem(i,k,j,p_so2)    = chem(i,k,j,p_so2)    + conv_rho*ebu(i,k,j,p_ebu_so2)
              chem(i,k,j,p_nh3)    = chem(i,k,j,p_nh3)    + conv_rho*ebu(i,k,j,p_ebu_nh3)
              chem(i,k,j,p_ch3cho) = chem(i,k,j,p_ch3cho) + conv_rho*ebu(i,k,j,p_ebu_ch3cho)
              chem(i,k,j,p_hcho)   = chem(i,k,j,p_hcho)   + conv_rho*ebu(i,k,j,p_ebu_hcho)
              chem(i,k,j,p_isop)   = chem(i,k,j,p_isop)   + conv_rho*ebu(i,k,j,p_ebu_isop)
              chem(i,k,j,p_c10h16) = chem(i,k,j,p_c10h16) + conv_rho*ebu(i,k,j,p_ebu_c10h16)
              chem(i,k,j,p_mgly)   = chem(i,k,j,p_mgly)   + conv_rho*ebu(i,k,j,p_ebu_mgly)
              chem(i,k,j,p_c2h4)   = chem(i,k,j,p_c2h4)   + conv_rho*ebu(i,k,j,p_ebu_c2h4)
              chem(i,k,j,p_par)    = chem(i,k,j,p_par)    + conv_rho*ebu(i,k,j,p_ebu_par)
              chem(i,k,j,p_ole)    = chem(i,k,j,p_ole)    + conv_rho*ebu(i,k,j,p_ebu_ole)
              chem(i,k,j,p_soaalk) = chem(i,k,j,p_soaalk) + conv_rho*ebu(i,k,j,p_ebu_soaalk)
              chem(i,k,j,p_ivoc)   = chem(i,k,j,p_ivoc)   + conv_rho*ebu(i,k,j,p_ebu_ivoc)
              chem(i,k,j,p_tol)    = chem(i,k,j,p_tol)    + conv_rho*ebu(i,k,j,p_ebu_tol)
! Add Aerosols
              chem(i,k,j,p_eci)   = chem(i,k,j,p_eci)  + conv*ebu(i,k,j,p_ebu_eci)
              chem(i,k,j,p_ecj)   = chem(i,k,j,p_ecj)  + conv*ebu(i,k,j,p_ebu_ecj) 
              chem(i,k,j,p_p25i)  = chem(i,k,j,p_p25i) + conv*ebu(i,k,j,p_ebu_pm25i)
              chem(i,k,j,p_p25j)  = chem(i,k,j,p_p25j) + conv*ebu(i,k,j,p_ebu_pm25j)
              chem(i,k,j,p_no3ai)  = chem(i,k,j,p_no3ai) + conv*ebu(i,k,j,p_ebu_no3i)
              chem(i,k,j,p_no3aj)  = chem(i,k,j,p_no3aj) + conv*ebu(i,k,j,p_ebu_no3j)
              chem(i,k,j,p_so4ai)  = chem(i,k,j,p_so4ai) + conv*ebu(i,k,j,p_ebu_so4i)
              chem(i,k,j,p_so4aj)  = chem(i,k,j,p_so4aj) + conv*ebu(i,k,j,p_ebu_so4j)

              chem(i,k,j,p_poa0j) = chem(i,k,j,p_poa0j) + poa_vol_params(1)*conv*ebu(i,k,j,p_ebu_orgj)
              chem(i,k,j,p_poa0i) = chem(i,k,j,p_poa0i) + poa_vol_params(1)*conv*ebu(i,k,j,p_ebu_orgi)
              chem(i,k,j,p_poa1j) = chem(i,k,j,p_poa1j) + poa_vol_params(2)*conv*ebu(i,k,j,p_ebu_orgj)
              chem(i,k,j,p_poa1i) = chem(i,k,j,p_poa1i) + poa_vol_params(2)*conv*ebu(i,k,j,p_ebu_orgi)
              chem(i,k,j,p_poa2j) = chem(i,k,j,p_poa2j) + poa_vol_params(3)*conv*ebu(i,k,j,p_ebu_orgj)
              chem(i,k,j,p_poa2i) = chem(i,k,j,p_poa2i) + poa_vol_params(3)*conv*ebu(i,k,j,p_ebu_orgi)
              chem(i,k,j,p_poa3j) = chem(i,k,j,p_poa3j) + poa_vol_params(4)*conv*ebu(i,k,j,p_ebu_orgj)
              chem(i,k,j,p_poa3i) = chem(i,k,j,p_poa3i) + poa_vol_params(4)*conv*ebu(i,k,j,p_ebu_orgi)
              IF ( config_flags%add_firewater==2 ) THEN
              qv(i,k,j)   = qv(i,k,j)   + 1.e-9 * conv*ebu(i,k,j,p_ebu_h2o)
              ENDIF

! Add to antha since we don't have a coarse mode bb tracer
!              chem(i,k,j,p_antha) = chem(i,k,j,p_antha) + conv*ebu(i,k,j,p_ebu_pm10)  
! Number density
              chem(i,k,j,p_nu0)   = chem(i,k,j,p_nu0)        +                     &
                                    conv*fact_numn*(                               &
                     anthfac*( ebu(i,k,j,p_ebu_pm25i) +  ebu(i,k,j,p_ebu_eci) ) +       &
                     so4fac*ebu(i,k,j,p_ebu_so4i) + no3fac*ebu(i,k,j,p_ebu_no3i) +    &
                     orgfac*ebu(i,k,j,p_ebu_orgi) )
              chem(i,k,j,p_ac0) = chem(i,k,j,p_ac0) +                   &
                                    conv*fact_numa*(              &
                     anthfac*( ebu(i,k,j,p_ebu_pm25j) +  ebu(i,k,j,p_ebu_ecj) ) +       &
                     so4fac*ebu(i,k,j,p_ebu_so4j) + no3fac*ebu(i,k,j,p_ebu_no3j) +     &
                     orgfac*ebu(i,k,j,p_ebu_orgj) )

!              chem(i,k,j,p_corn) = chem(i,k,j,p_corn) + conv*fact_numc*anthfac*ebu(i,k,j,p_ebu_pm10)


             if (icall<5000 .AND. config_flags%debug_chem) then
               if ( k==kts ) then
                 WRITE(6,*) 'add_emiss_burn: ktau,gmt,dtstep,time_int ',ktau,gmt,dtstep,time_int
                 WRITE(*,*) 'add_emiss_burn: i,j,xlat(i,j),xlong(i,j) ',i,j,xlat(i,j),xlong(i,j)
                 WRITE(*,*) 'add_emiss_burn: luf_igbp(i,:,j) ',luf_igbp(i,:,j)
                 WRITE(*,*) 'add_emiss_burn: lu_fire1(i,j) ',lu_fire1(i,j)
                 WRITE(6,*) 'add_emiss_burn: timeq,peak_hr(i,j),fhist(i,j),r_q(i,j) ',timeq,peak_hr(i,j),fhist(i,j),r_q(i,j)
                 WRITE(*,*) 'add_emiss_burn: rainc(i,j),rainnc(i,j) ', rainc(i,j),rainnc(i,j)
                 icall= icall+1
               endif
               if ( k==kts .OR. k==kfire_max ) then
                 WRITE(6,*) 'add_emiss_burn: i,j,k ',i,j,k
                 WRITE(6,*) 'add_emiss_burn: rho_phy(i,k,j),dz8w(i,k,j),conv_rho,conv ',rho_phy(i,k,j),dz8w(i,k,j),conv_rho,conv
                 WRITE(6,*) 'add_emiss_burn: ebu(i,k,j,p_ebu_no),dm_no ', ebu(i,k,j,p_no),conv_rho*ebu(i,k,j,p_ebu_no)
               endif
             endif


              enddo
            enddo
          enddo

          ext2= sc_me + ab_me   
          do j=jts,jte
           do k=kts,kte
            do i=its,ite
 ! Fraction of nh4 as amm-sulf and amm-nitr for i/j modes
             alpha_i = 2. * chem(i,k,j,p_so4ai) / &
                       (  chem(i,k,j,p_no3ai) +  2. *  chem(i,k,j,p_so4ai) )
             alpha_j = 2. * chem(i,k,j,p_so4aj) / &
                       (  chem(i,k,j,p_no3aj) +  2. *  chem(i,k,j,p_so4aj) )
 ! Relative humidity growth factors ( modeled as exponentials instead of reading  table
             IF ( rel_hum(i,k,j) <= 0.36 ) THEN
                FS_RH = 1.
                FL_RH = 1.
             ELSEIF ( rel_hum(i,k,j) > 0.95 ) THEN
                FS_RH = 9.34
                FL_RH = 5.57
             ELSE
                FS_RH = 1.719 + 4.635E-05 * exp(12.51*rel_hum(i,k,j))
                FL_RH = 1.481 + 4.422E-04 * exp(9.484*rel_hum(i,k,j))
             ENDIF 
             IF ( rel_hum(i,k,j) <= 0.46 ) THEN
                F_SS = 1.
             ELSEIF ( rel_hum(i,k,j) > 0.95 ) THEN
                F_SS = 7.3492
             ELSE
                F_SS = 2.5177 + 0.0005 * exp(9.5155*rel_hum(i,k,j))
             ENDIF 
 ! IMPROVE extinction algorithm: Pitchford et al., JAWMA, 2007 ( Mm^-1 )
             aerext = (    sm_ammsulf_fact  * FS_RH * ( chem(i,k,j,p_so4ai) + &
                                      alpha_i * chem(i,k,j,p_nh4ai) )        + &
                           lg_ammsulf_fact  * FL_RH * ( chem(i,k,j,p_so4aj) + &
                                      alpha_j * chem(i,k,j,p_nh4ai) )        + &
                           sm_ammnitr_fact  * FS_RH * ( chem(i,k,j,p_no3ai) + &
                                      (1. - alpha_i) * chem(i,k,j,p_nh4ai) ) + &
                           lg_ammnitr_fact  * FL_RH * ( chem(i,k,j,p_no3aj) + &
                                      (1. - alpha_j) * chem(i,k,j,p_nh4aj) ) + &
                           sm_org_fact * ( chem(i,k,j,p_poa0i) +  &
                              chem(i,k,j,p_poa1i)  +  &
                              chem(i,k,j,p_poa2i)  +  &
                              chem(i,k,j,p_poa3i)  +  &
                              chem(i,k,j,p_asoa0i) +  &
                              chem(i,k,j,p_asoa1i) +  &
                              chem(i,k,j,p_asoa2i) +  &
                              chem(i,k,j,p_asoa3i) +  &
                              chem(i,k,j,p_bsoa1i) +  &
                              chem(i,k,j,p_bsoa2i) +  &
                              chem(i,k,j,p_bsoa3i) ) +  &
                           lg_org_fact  * ( chem(i,k,j,p_poa0j) +  &
                              chem(i,k,j,p_poa1j)  +  &
                              chem(i,k,j,p_poa2j)  +  &
                              chem(i,k,j,p_poa3j)  +  &
                              chem(i,k,j,p_asoa0j) +  &
                              chem(i,k,j,p_asoa1j) +  &
                              chem(i,k,j,p_asoa2j) +  &
                              chem(i,k,j,p_asoa3j) +  &
                              chem(i,k,j,p_bsoa1j) +  &
                              chem(i,k,j,p_bsoa2j) +  &
                              chem(i,k,j,p_bsoa3j) )   +  &
                           ec_fact  * ( chem(i,k,j,p_eci) + chem(i,k,j,p_ecj) ) + &
                           0.6 * (chem(i,k,j,p_soila) + chem(i,k,j,p_antha) ) +   &
                           1.7 * F_SS * (chem(i,k,j,p_clai) + chem(i,k,j,p_naai) + &
                                         chem(i,k,j,p_claj) + chem(i,k,j,p_naaj) + &
                                         chem(i,k,j,p_seas) ) + &
                           1.0 * (chem(i,k,j,p_p25i) + chem(i,k,j,p_p25j)) ) * rho_phy(i,k,j) + &
                           no2_fact * 1.e3 * chem(i,k,j,p_no2) + rayleigh_scat_param
                          
! Convert to ( m ) and multiply by height of layer
             aod3d(i,k,j)= 1.e-6 * aerext * dz8w(i,k,j)
            enddo
           enddo
          enddo

     IF ( icall<2000 .AND. config_flags%debug_chem ) then
         WRITE(*,*) 'add_emis_burn: i,j,k,ext2 ',i,j,k,ext2
         WRITE(*,*) 'add_emis_burn: rel_hum(its,kts,jts),rel_hum(ite,kfire_max,jte) ',rel_hum(its,kts,jts),rel_hum(ite,kfire_max,jte)
         WRITE(*,*) 'add_emis_burn: aod3d(its,kts,jts),aod3d(ite,kfire_max,jte) ',aod3d(its,kts,jts),aod3d(ite,kfire_max,jte)
     END IF



!     CASE DEFAULT
!       call wrf_debug(15,'nothing done with burn emissions for chem array')
!    END SELECT emiss_select

    END subroutine add_emis_burn

END module module_add_emiss_burn

 module module_plumerise1
!- Implementing the fire radiative power (FRP) methodology for biomass burning
!- emissions and convective energy estimation.
!- Saulo Freitas, Gabriel Pereira (INPE/UFJS, Brazil)
!- Ravan Ahmadov, Georg Grell (NOAA, USA)
!- The flag "plumerise_flag" defines the method:
!-    =1 => original method
!-    =2 => FRP based
!-------------------------------------------------------------------------
!
! use module_zero_plumegen_coms
!  integer, parameter :: nveg_agreg      = 4
! integer, parameter :: tropical_forest = 1
! integer, parameter :: boreal_forest   = 2
! integer, parameter :: savannah        = 3

! integer, parameter :: grassland       = 4
!  real, dimension(nveg_agreg) :: firesize,mean_fct
! character(len=20), parameter :: veg_name(nveg_agreg) = (/ &
!                              'Tropical-Forest', &
!                              'Boreal-Forest  ', &
!                              'Savanna        ', &
!                              'Grassland      ' /)
! character(len=20), parameter :: spc_suf(nveg_agreg) = (/ &
!                              'agtf' , &  ! trop forest
!                              'agef' , &  ! extratrop forest
!                              'agsv' , &  ! savanna
!                              'aggr'   /) ! grassland

CONTAINS
subroutine plumerise_driver (id,                                                            &
                             flam_frac,ebu,ebu_in,hwp,                                      &
                             grasslandfrac,                                                 &
                             config_flags,                                                  &
                             t_phy,q_vap,                                                   &       ! RAR: moist is replaced with q_vap
                             rho_phy,vvel,u_phy,v_phy,p_phy,                                &
                             z_at_w,z,ktau,                                                 &       ! scale_fire_emiss is part of config_flags
                             plume_fre, k_min, k_max,                                       &       ! RAR:
                             ids,ide, jds,jde, kds,kde,                                     &
                             ims,ime, jms,jme, kms,kme,                                     &
                             its,ite, jts,jte, kts,kte                                      )

  USE module_configure
  USE module_model_constants
  USE module_state_description
  USE module_chem_plumerise_scalar
  USE module_zero_plumegen_coms
!  USE module_smoke_plumerise
  IMPLICIT NONE
! integer, parameter :: nveg_agreg      = 4
! integer, parameter :: nveg_agreg      = 4
! integer, parameter :: tropical_forest = 1
! integer, parameter :: boreal_forest   = 2
! integer, parameter :: savannah        = 3

!-fre
!   INTEGER, PARAMETER ::  iflam_frac  =1

   REAL, DIMENSION( ims:ime, jms:jme, 4 ), INTENT(IN ) :: plume_fre          ! RAR: FRP etc. array

   TYPE(grid_config_rec_type),  INTENT(IN )    :: config_flags

   INTEGER,      INTENT(IN   ) :: id,ktau,                                 &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
!   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
!         INTENT(IN ) ::                                   moist
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_ebu ), INTENT(INOUT ) ::  ebu

   REAL, DIMENSION( ims:ime, 1 , jms:jme, num_ebu_in ), INTENT(IN )  :: ebu_in ! Edit this
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT ) :: flam_frac
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN  ) :: hwp

!   REAL, DIMENSION( ims:ime, 1, jms:jme, num_ebu_in ),                 &
!         INTENT(IN ) ::                                   ebu_in
!   REAL, DIMENSION( ims:ime, jms:jme ),                 &
!         INTENT(IN ) ::                                                &
!           mean_fct_agtf,mean_fct_agef,&
!           mean_fct_agsv,mean_fct_aggr,firesize_agtf,firesize_agef,       &
!           firesize_agsv,firesize_aggr

   REAL, DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   t_phy,z,z_at_w,vvel,u_phy,v_phy,rho_phy,p_phy,q_vap                     ! RAR
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN ) :: grasslandfrac
  ! REAL, INTENT(IN   ) ::      dtstep

! Local variables...
      INTEGER :: nv, i, j, k,  kp1, kp2
      INTEGER, DIMENSION(ims:ime, jms:jme), INTENT (OUT) :: k_min, k_max      ! Min and max ver. levels for BB injection spread
      !real, dimension (num_ebu) :: eburn_in
      !real, dimension (kte,num_ebu) :: eburn_out
      real, dimension (kte) :: u_in ,v_in ,w_in ,theta_in ,pi_in, rho_phyin ,qv_in ,zmid, z_lev
      REAL :: dz_plume

      !INTEGER, PARAMETER :: kfire_max=30  
      INTEGER, SAVE :: icall
! real, dimension(nveg_agreg) :: firesize,mean_fct
!      real :: sum, ffirs, ratio
!     real,save,dimension(its:ite,jts:jte) :: ffirs
!      nspecies=num_ebu
!     write(0,*)'plumerise'

! RAR:
!       if (config_flags%biomass_burn_opt == BIOMASSB_REDHC) then
          do j=jts,jte
             do i=its,ite
                ebu(i,kts,j,p_ebu_so2)    = ebu_in(i,1,j,p_ebu_in_so2)
                ebu(i,kts,j,p_ebu_no)     = ebu_in(i,1,j,p_ebu_in_no)
                ebu(i,kts,j,p_ebu_co)     = ebu_in(i,1,j,p_ebu_in_co)
                ebu(i,kts,j,p_ebu_no2)    = ebu_in(i,1,j,p_ebu_in_no2)
                ebu(i,kts,j,p_ebu_hono)   = ebu_in(i,1,j,p_ebu_in_hono)
                ebu(i,kts,j,p_ebu_nh3)    = ebu_in(i,1,j,p_ebu_in_nh3)
                ebu(i,kts,j,p_ebu_ch3cho) = ebu_in(i,1,j,p_ebu_in_ch3cho)
                ebu(i,kts,j,p_ebu_hcho)   = ebu_in(i,1,j,p_ebu_in_hcho)   
                ebu(i,kts,j,p_ebu_isop)   = ebu_in(i,1,j,p_ebu_in_isop)
                ebu(i,kts,j,p_ebu_c10h16) = ebu_in(i,1,j,p_ebu_in_c10h16)
                ebu(i,kts,j,p_ebu_mgly)   = ebu_in(i,1,j,p_ebu_in_mgly)
                ebu(i,kts,j,p_ebu_c2h4)   = ebu_in(i,1,j,p_ebu_in_c2h4)
                ebu(i,kts,j,p_ebu_par)    = ebu_in(i,1,j,p_ebu_in_par)
                ebu(i,kts,j,p_ebu_ole)    = ebu_in(i,1,j,p_ebu_in_ole)
                ebu(i,kts,j,p_ebu_tol)    = ebu_in(i,1,j,p_ebu_in_tol)
                ebu(i,kts,j,p_ebu_soaalk) = ebu_in(i,1,j,p_ebu_in_soaalk)
                ebu(i,kts,j,p_ebu_ivoc)   = ebu_in(i,1,j,p_ebu_in_ivoc)
                ebu(i,kts,j,p_ebu_pm25i)  = ebu_in(i,1,j,p_ebu_in_pm25i)
                ebu(i,kts,j,p_ebu_pm25j)  = ebu_in(i,1,j,p_ebu_in_pm25j)
                ebu(i,kts,j,p_ebu_so4i)   = ebu_in(i,1,j,p_ebu_in_so4i)
                ebu(i,kts,j,p_ebu_so4j)   = ebu_in(i,1,j,p_ebu_in_so4j)
                ebu(i,kts,j,p_ebu_no3i)   = ebu_in(i,1,j,p_ebu_in_no3i)
                ebu(i,kts,j,p_ebu_no3j)   = ebu_in(i,1,j,p_ebu_in_no3j)
                ebu(i,kts,j,p_ebu_orgi)   = ebu_in(i,1,j,p_ebu_in_orgi)
                ebu(i,kts,j,p_ebu_orgj)   = ebu_in(i,1,j,p_ebu_in_orgj)
                ebu(i,kts,j,p_ebu_eci)    = ebu_in(i,1,j,p_ebu_in_eci)
                ebu(i,kts,j,p_ebu_ecj)    = ebu_in(i,1,j,p_ebu_in_ecj)
                ebu(i,kts,j,p_ebu_h2o)    = ebu_in(i,1,j,p_ebu_in_h2o)
                ebu(i,kts,j,p_ebu_pm10)   = ebu_in(i,1,j,p_ebu_in_pm10)
            enddo
          enddo
   
        IF (icall<2000 .AND. config_flags%debug_chem) then
           WRITE(*,*) 'module_plumerise1: its,ite,jts,jte ', its,ite,jts,jte
           WRITE(*,*) 'module_plumerise1: ims,ime,jms,jme ', ims,ime,jms,jme
           WRITE(*,*) 'SIZE OF EBU_IN IS: ',SIZE(EBU_IN)
           WRITE(*,*) 'SIZE OF EBU IS:',SIZE(EBU)
           WRITE(*,*) 'module_plumerise1: p_ebu_no,p_ebu_co,p_ebu_oc,p_ebu_pm25,num_ebu: ', p_ebu_no,p_ebu_co,p_ebu_oc,p_ebu_pm25,num_ebu
           WRITE(*,*) 'module_plumerise1: p_ebu_in_no,p_ebu_in_co,p_ebu_in_oc,p_ebu_in_pm25,num_ebu_in: ', p_ebu_in_no,p_ebu_in_co,p_ebu_in_oc,p_ebu_in_pm25,num_ebu_in
           WRITE(*,*) 'module_plumerise1: maxval(ebu(:,kts,:,p_no)) ', maxval(ebu(:,kts,:,p_ebu_no))
           WRITE(*,*) 'module_plumerise1: maxval(ebu(:,kts,:,p_co)) ', maxval(ebu(:,kts,:,p_ebu_co)) 
           WRITE(*,*) 'module_plumerise1: maxval(ebu(:,kts,:,p_oc)) ', maxval(ebu(:,kts,:,p_ebu_oc))
           WRITE(*,*) 'module_plumerise1: maxval(ebu(:,kts,:,p_bc)) ', maxval(ebu(:,kts,:,p_ebu_bc))
           WRITE(*,*) 'module_plumerise1: maxval(ebu(:,kts,:,p_pm25)) ', maxval(ebu(:,kts,:,p_ebu_pm25))
           icall=icall+1
         END IF
!      endif

! RAR: setting to zero the ebu emissions at the levels k>1, this is necessary when the plumerise is called, so the emissions at k>1 are updated
       do nv=2,num_ebu
          do j=jts,jte
            do k=kts+1,kte
               do i=its,ite
                 ebu(i,k,j,nv)=0.
               enddo
            enddo
          enddo
       enddo
     
! For now the flammable fraction is constant, based on the namelist. The next
! step to use LU index and meteorology to parameterize it
!    IF (ktau==2) THEN
       do j=jts,jte
        do i=its,ite
           flam_frac(i,j)= config_flags%flam_part;
        enddo
       enddo
 !   ENDIF


! RAR: new FRP based approach
check_pl:  IF (config_flags%plumerise_flag == 2 ) THEN    ! if the namelist option is set for plumerise 
       do j=jts,jte
          do i=its,ite
               k_min(i,j)=0
               k_max(i,j)=0
            
            if (plume_fre(i,j,1) < 1.e+7 .or. grasslandfrac(i,j) > 0.70 ) cycle  ! FRP is too small, <0.1 MW -- JLS modified to 10 MW
               !do nv=2,num_ebu
               !   eburn_in(nv)= ebb_smoke(i,j)
               !enddo

               do k=kts,kte
                  u_in(k)=  u_phy(i,k,j)
                  v_in(k)=  v_phy(i,k,j)
                  w_in(k)=  vvel(i,k,j)
                  qv_in(k)= q_vap(i,k,j)    ! RAR: moist(i,k,j,p_qv)
                  pi_in(k)= cp*(p_phy(i,k,j)/p1000mb)**rcp
                  zmid(k)=  z(i,k,j)-z_at_w(i,kts,j)
                  z_lev(k)= z_at_w(i,k,j)-z_at_w(i,kts,j)
                  rho_phyin(k)= rho_phy(i,k,j)
                  theta_in(k)= t_phy(i,k,j)/pi_in(k)*cp
               enddo

             IF (icall<2000 .AND. config_flags%debug_chem) then
               WRITE(*,*) 'module_plumerise1: i,j ',i,j
               WRITE(*,*) 'module_plumerise1: plume_fre(i,j,:) ',plume_fre(i,j,:)
               WRITE(*,*) 'module_plumerise1: ebu(i,kts,j,p_ebu_no) ',ebu(i,kts,j,p_ebu_no)
               WRITE(*,*) 'module_plumerise1: u_in(10),v_in(10),w_in(kte),qv_in(10),pi_in(10) ',u_in(10),v_in(10),w_in(kte),qv_in(10),pi_in(10)
               WRITE(*,*) 'module_plumerise1: zmid(kte),z_lev(kte),rho_phyin(kte),theta_in(kte) ',zmid(kte),z_lev(kte),rho_phyin(kte),theta_in(kte)
             END IF

! RAR: the plume rise calculation step:
!               CALL plumerise(kte,1,1,1,1,1,1,         &
!                              !firesize,mean_fct,                    & 
                              !num_ebu, eburn_in, eburn_out,         &
!                              u_in, v_in, w_in, theta_in ,pi_in,    &
!                              rho_phyin, qv_in, zmid, z_lev,        &
!                              plume_fre(i,j,:), k_min(i,j),         & 
!                              k_max(i,j), ktau, config_flags%debug_chem )
! RAR: the plume rise calculation step:
               CALL plumerise(kte,1,1,1,1,1,1,         &
                              u_in, v_in, w_in, theta_in ,pi_in,    &
                              rho_phyin, qv_in, zmid, z_lev,        &
                              plume_fre(i,j,:), k_min(i,j),         &
                              k_max(i,j), ktau, config_flags%debug_chem, &
                              config_flags%scale_factor_hfx_fire )

!               call plumerise(kte,1,1,1,1,1,1,                      &
!                              firesize,mean_fct,                    &
!                              num_ebu,eburn_in,eburn_out,           &
!                              u_in ,v_in ,w_in ,theta_in ,pi_in,    &
!                              rho_phyin ,qv_in ,zmid, z_lev,        &
!                              config_flags%plumerise_flag,          &
!                              plume_fre(i,j,:), k_min(i,j), k_max(i,j) )

               kp1= k_min(i,j)
               kp2= k_max(i,j)   
               dz_plume= z_at_w(i,kp2,j) - z_at_w(i,kp1,j)
! RAR: add here
               do nv=2,num_ebu
                  do k=kp1,kp2-1     
                     ebu(i,k,j,nv)= flam_frac(i,j)* ebu_in(i,1,j,nv)* (z_at_w(i,k+1,j)-z_at_w(i,k,j))/dz_plume
                  enddo
                  ebu(i,kts,j,nv)= (1.-flam_frac(i,j))* ebu_in(i,1,j,nv)
               enddo

               IF (icall<2000 .AND. config_flags%debug_chem) then
                   WRITE(*,*) 'module_plumerise1: kts,flam_frac(i,j) ',kts,flam_frac(i,j) 
                   WRITE(*,*) 'module_plumerise1: k_min(i,j), k_max(i,j) ',k_min(i,j), k_max(i,j)  
               END IF
            enddo
          enddo

        ENDIF check_pl

end subroutine plumerise_driver

END module module_plumerise1

module module_gocart_dust_fengsha
!
!  This module developed by Barry Baker (NOAA ARL)
!  For serious questions contact barry.baker@noaa.gov
!
!  07/16/2019 - Adapted for NUOPC/GOCART, R. Montuoro
!  02/01/2020 - Adapted for FV3/CCPP, Haiqin Li

  use module_dust_data
  USE module_data_sorgam, ONLY: factnuma,factnumc,factnumn,soilfac
  USE module_configure
  USE dust_fengsha
  USE module_state_description

  implicit none

!  private

!  public :: gocart_dust_fengsha_driver

contains

  subroutine gocart_dust_fengsha_driver(dt,config_flags, &
       chem,rho_phy,smois,p8w,ssm,                       &
       isltyp,vegfra,snowh,xland,area,emis_dust,       &
       ust,znt,clay,sand,rdrag,uthr,                     &
       num_emis_dust,num_moist,num_chem,num_soil_layers, &
       ids,ide, jds,jde, kds,kde,                        &
       ims,ime, jms,jme, kms,kme,                        &
       its,ite, jts,jte, kts,kte)
    IMPLICIT NONE

    TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags
    INTEGER,      INTENT(IN   ) ::                       &
         ids,ide, jds,jde, kds,kde,                      &
         ims,ime, jms,jme, kms,kme,                      &
         its,ite, jts,jte, kts,kte,                      &
         num_emis_dust,num_moist,num_chem,num_soil_layers
    INTEGER,DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: isltyp
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ), INTENT(INOUT) :: chem
    REAL, DIMENSION( ims:ime, 1, jms:jme,num_emis_dust),OPTIONAL, INTENT(INOUT) :: emis_dust
    REAL, DIMENSION( ims:ime, num_soil_layers, jms:jme ), INTENT(IN) :: smois
    REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: ssm
    REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: vegfra,     &
                                                        snowh,      &
                                                        xland,      &
                                                        ust,        &
                                                        znt,        &
                                                        clay,       &
                                                        sand,       &
                                                        rdrag,      &
                                                        uthr
    REAL, DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) :: &
         p8w,             &
         rho_phy
    REAL, INTENT(IN) :: dt,area

    ! Local variables
    integer, parameter :: nsmois_lay=2 ! number of soil moisture layers to average over
    integer :: nmx,i,j,k,n
    integer :: dust_calcdrag
    integer :: ilwi
    real(8) :: erodtot
    REAL(8) :: smois_avg 
    REAL(8) :: drylimit
    real(8), dimension(ndust) :: tc,bems
    real(8) :: airden,airmas,ustar
    real(8) :: dxy
    real(8), dimension (3) :: massfrac
    real(8) :: conver,converi
    real(8) :: R, g

    real ::   dstfracna,dstfrac
    real ::   dust_alpha,dust_gamma

    dstfracna = 0.0236
    dstfrac   = 1.0 - dstfracna

    dust_alpha = config_flags%dust_alpha
    dust_gamma = config_flags%dust_gamma
    g = 9.81

    dust_calcdrag = 1 ! hard-coded to always use input - could/should make namelist

    ! threshold values
    conver=1.e-9
    converi=1.e9

    ! Number of dust bins

    nmx=ndust

    k=kts
    do j=jts,jte
       do i=its,ite

          ! Don't do dust over water!!!

          ilwi=0
          if(xland(i,j).lt.1.5)then
             ilwi=1

             ! Total concentration at lowest model level. This is still hardcoded for 5 bins.

             !    if(config_flags%chem_opt == 2 .or. config_flags%chem_opt == 11 ) then
             !       tc(:)=1.e-16*conver
             !    else
             tc(1)=0. !chem(i,kts,j,p_dust_1)*conver
             tc(2)=0. !chem(i,kts,j,p_dust_2)*conver
             tc(3)=0. !chem(i,kts,j,p_dust_3)*conver
             tc(4)=0. !chem(i,kts,j,p_dust_4)*conver
             tc(5)=0. !chem(i,kts,j,p_dust_5)*conver
             !    endif

             ! Air mass and density at lowest model level.

             airmas=-(p8w(i,kts+1,j)-p8w(i,kts,j))*area/g
             airden=rho_phy(i,kts,j)
             ustar=ust(i,j)
             dxy=area

             ! Mass fractions of clay, silt, and sand.
             massfrac(1)=clay(i,j)
             massfrac(2)=1-(clay(i,j)+sand(i,j))
             massfrac(3)=sand(i,j)


             ! Total erodibility.
             
             erodtot = ssm(i,j) ! SUM(erod(i,j,:))
             
             ! Don't allow roughness lengths greater than 20 cm to be lofted.
             ! This kludge accounts for land use types like urban areas and
             ! forests which would otherwise show up as high dust emitters.
             ! This is a placeholder for a more widely accepted kludge
             ! factor in the literature, which reduces lofting for rough areas.
             ! Forthcoming...

             IF (znt(i,j) .gt. 0.2) then
                ilwi=0
             endif

             ! limit where there is lots of vegetation
             if (vegfra(i,j) .gt. 17.) then
                ilwi = 0
             endif

             ! limit where there is snow on the ground
             if (snowh(i,j) .gt. 0) then
                ilwi = 0
             endif

             ! Do not allow areas with bedrock, lava, or land-ice to loft

             IF (isltyp(i,j) .eq. 15 .or. isltyp(i,j) .eq. 16. .or. &
                  isltyp(i,j) .eq. 18) then
                ilwi=0
             ENDIF
             IF (isltyp(i,j) .eq. 0)then
                ilwi=0
             endif
             if(ilwi == 0 ) cycle

             ! Calculate gravimetric soil moisture and drylimit.
             ! JLS - take the average of teh first nsmois_layers
!             gravsm=0.
!             do n=1,nsmois_lay 
!                gravsm=gravsm+100.*smois(i,n,j)/((1.-maxsmc(isltyp(i,j)))*(2.65*(1.-clay(i,j))+2.50*clay(i,j)))
!             enddo
!             gravsm=gravsm/nsmois_lay
             smois_avg=0.
             do n=1,nsmois_lay
                smois_avg=smois_avg+smois(i,n,j)
             enddo
             smois_avg=smois_avg/nsmois_lay
             
             drylimit=14.0*clay(i,j)*clay(i,j)+17.0*clay(i,j)

             ! get drag partition
             ! FENGSHA uses the drag partition correction of MacKinnon et al 2004
             !     doi:10.1016/j.geomorph.2004.03.009
             if (dust_calcdrag .ne. 1) then
                call fengsha_drag(znt(i,j),R)
             else
                ! use the precalculated version derived from ASCAT; Prigent et al. (2012,2015)
                ! doi:10.1109/TGRS.2014.2338913 & doi:10.5194/amt-5-2703-2012
                ! pick only valid values
                if (rdrag(i,j) > 0.) then
                  R = real(rdrag(i,j), kind=8)
                else
                  cycle
                endif
             endif

             ! Call dust emission routine.
             call source_dust( nmx, dt, tc, ustar, massfrac, & 
                  erodtot, dxy, smois_avg, airden, airmas, &
                  bems, g, drylimit, dust_alpha, dust_gamma, R, uthr(i,j))

            if ( config_flags%dust_ijk_opt == 1 ) then

           chem(i,kts,j,p_naai) = chem(i,kts,j,p_naai) &
                         +.25*(tc(1)+.286*tc(2))*converi * dstfracna
!             if (config_flags%do_isorropia ) then
!                  mgai(i,kts,j) = mgai(i,kts,j)        &
!                         +.25*(tc(1)+.286*tc(2))*converi * dstfracmg
!                  caai(i,kts,j) = caai(i,kts,j)        &
!                         +.25*(tc(1)+.286*tc(2))*converi * dstfracca
!                   kai(i,kts,j) = kai(i,kts,j)         &
!                         +.25*(tc(1)+.286*tc(2))*converi * dstfrack
!             endif
           chem(i,kts,j,p_p25i)=chem(i,kts,j,p_p25i) &
                         +.25*(tc(1)+.286*tc(2))*converi  * dstfrac

           chem(i,kts,j,p_naaj)=chem(i,kts,j,p_naaj) &
                         +.75*(tc(1)+.286*tc(2))*converi * dstfracna
!             if ( config_flags%do_isorropia ) then
!                       mgaj(i,kts,j) = mgaj(i,kts,j) &
!                         +.75*(tc(1)+.286*tc(2))*converi * dstfracmg
!                       caaj(i,kts,j) = caaj(i,kts,j) &
!                         +.75*(tc(1)+.286*tc(2))*converi * dstfracca
!                       kaj(i,kts,j) = kaj(i,kts,j) &
!                         +.75*(tc(1)+.286*tc(2))*converi * dstfrack
!             endif
           chem(i,kts,j,p_p25j)=chem(i,kts,j,p_p25j) &
                         +.75*(tc(1)+.286*tc(2))*converi  * dstfrac
           chem(i,kts,j,p_soila)=chem(i,kts,j,p_soila) &
                         +(.714*tc(2)+tc(3)+tc(4))*converi
           chem(i,kts,j,p_nu0)  =  chem(i,kts,j,p_nu0) &
                         + .25*(tc(1)+.286*tc(2))*converi*factnumn*soilfac
           chem(i,kts,j,p_ac0)  =  chem(i,kts,j,p_ac0) &
                         + .75*(tc(1)+.286*tc(2))*converi*factnuma*soilfac
           chem(i,kts,j,p_corn) =  chem(i,kts,j,p_corn) &
                         + (.714*tc(2)+tc(3)+tc(4))*converi*factnumc*soilfac
          elseif ( config_flags%dust_ijk_opt == 2 ) then

           chem(i,kts,j,p_p25j) =chem(i,kts,j,p_p25j) +  & 
                                 0.07*sum(tc(1:5))*converi*dstfrac 
           chem(i,kts,j,p_naaj) =chem(i,kts,j,p_naaj) + &
                                 0.07*sum(tc(1:5))*converi*dstfracna
!                if (config_flags%do_isorropia ) then
!                  mgaj(i,kts,j) =mgaj(i,kts,j)        +
!                  0.07*sum(tc(1:5))*converi*dstfracmg
!                  caaj(i,kts,j) =caaj(i,kts,j)        +
!                  0.07*sum(tc(1:5))*converi*dstfracca
!                  kaj(i,kts,j) =kaj(i,kts,j)         +
!                  0.07*sum(tc(1:5))*converi*dstfrack
!                endif
           chem(i,kts,j,p_soila)=chem(i,kts,j,p_soila) + &
                                  0.93*1.02*sum(tc(1:5))*converi
           chem(i,kts,j,p_ac0) =  chem(i,kts,j,p_ac0) +  &
                                  0.07*sum(tc(1:5))*converi*factnuma* soilfac
           chem(i,kts,j,p_corn) =  chem(i,kts,j,p_corn) + &
                                   0.93*1.02*sum(tc(1:5))*converi*factnumc*soilfac
          endif

             ! For output diagnostics

             emis_dust(i,1,j,p_edust1)=bems(1)
             emis_dust(i,1,j,p_edust2)=bems(2)
             emis_dust(i,1,j,p_edust3)=bems(3)
             emis_dust(i,1,j,p_edust4)=bems(4)
             emis_dust(i,1,j,p_edust5)=bems(5)
          endif
       enddo
    enddo
    !

  end subroutine gocart_dust_fengsha_driver


  SUBROUTINE source_dust( nmx, dt1, tc, ustar, massfrac, &
       erod, dxy, smois, airden, airmas, bems, g0, drylimit, alpha,  &
       gamma, R, uthres)
    
    ! ****************************************************************************
    ! *  Evaluate the source of each dust particles size bin by soil emission
    ! *
    ! *  Input:
    ! *         EROD      Fraction of erodible grid cell                (-)
    ! *         GRAVSM    Gravimetric soil moisture                     (g/g)
    ! *         DRYLIMIT  Upper GRAVSM limit for air-dry soil           (g/g)
    ! *         ALPHA     Constant to fudge the total emission of dust  (1/m)
    ! *         GAMMA     Tuning constant for erodibility               (-)
    ! *         DXY       Surface of each grid cell                     (m2)
    ! *         AIRMAS    Mass of air for each grid box                 (kg)
    ! *         AIRDEN    Density of air for each grid box              (kg/m3)
    ! *         USTAR     Friction velocity                             (m/s)
    ! *         DT1       Time step                                     (s)
    ! *         NMX       Number of dust bins                           (-)
    ! *         SMX       Number of saltation bins                      (-)
    ! *         IMX       Number of I points                            (-)
    ! *         JMX       Number of J points                            (-)
    ! *         LMX       Number of L points                            (-)
    ! *         R         Drag Partition                                (-)
    ! *         UTHRES    FENGSHA Dry Threshold Velocities              (m/s)
    ! *
    ! *  Data:
    ! *         MASSFRAC  Fraction of mass in each of 3 soil classes    (-)
    ! *         SPOINT    Pointer to 3 soil classes                     (-)
    ! *         DEN_DUST  Dust density                                  (kg/m3)
    ! *         DEN_SALT  Saltation particle density                    (kg/m3)
    ! *         REFF_SALT Reference saltation particle diameter         (m)
    ! *         REFF_DUST Reference dust particle diameter              (m)
    ! *         LO_DUST   Lower diameter limits for dust bins           (m)
    ! *         UP_DUST   Upper diameter limits for dust bins           (m)
    ! *         FRAC_SALT Soil class mass fraction for saltation bins   (-)
    ! *
    ! *  Parameters:
    ! *         CMB       Constant of proportionality                   (-)
    ! *         MMD_DUST  Mass median diameter of dust                  (m)
    ! *         GSD_DUST  Geometric standard deviation of dust          (-)
    ! *         LAMBDA    Side crack propagation length                 (m)
    ! *         CV        Normalization constant                        (-)
    ! *         G0        Gravitational acceleration                    (m/s2)
    ! *         G         Gravitational acceleration in cgs             (cm/s2)
    ! *
    ! *  Working:
    ! *         U_TS0     "Dry" threshold friction velocity             (m/s)
    ! *         U_TS      Moisture-adjusted threshold friction velocity (m/s)
    ! *         RHOA      Density of air in cgs                         (g/cm3)
    ! *         DEN       Dust density in cgs                           (g/cm3)
    ! *         DIAM      Dust diameter in cgs                          (cm)
    ! *         DMASS     Saltation mass distribution                   (-)
    ! *         DSURFACE  Saltation surface area per unit mass          (m2/kg)
    ! *         DS_REL    Saltation surface area distribution           (-)
    ! *         SALT      Saltation flux                                (kg/m/s)
    ! *         DLNDP     Dust bin width                                (-)
    ! *         EMIT      Total vertical mass flux                      (kg/m2/s)
    ! *         EMIT_VOL  Total vertical volume flux                    (m/s)
    ! *         DSRC      Mass of emitted dust               (kg/timestep/cell)
    ! *
    ! *  Output:
    ! *         TC        Total concentration of dust        (kg/kg/timestep/cell)
    ! *         BEMS      Source of each dust type           (kg/timestep/cell)
    ! *
    ! ****************************************************************************
    implicit none

    INTEGER,            INTENT(IN)    :: nmx
    REAL, INTENT(IN)    :: dt1
    REAL(8), INTENT(INOUT) :: tc(nmx)
    REAL(8), INTENT(IN)    :: ustar
    REAL(8), INTENT(IN)    :: massfrac(3)
    REAL(8), INTENT(IN)    :: erod
    REAL(8), INTENT(IN)    :: dxy
    REAL(8), INTENT(IN)    :: smois 
    REAL(8), INTENT(IN)    :: airden
    REAL(8), INTENT(IN)    :: airmas
    REAL(8), INTENT(OUT)   :: bems(nmx)
    REAL(8), INTENT(IN)    :: g0
    REAL(8), INTENT(IN)    :: drylimit
    !
    !  REAL, PARAMETER :: alpha=1.8E-8  ! (m^-1)
    REAL, INTENT(IN)    :: alpha
    ! Experimental optional exponential tuning constant for erodibility.
    ! 0 < gamma < 1 -> more relative impact by low erodibility regions.
    REAL, INTENT(IN)    :: gamma
    REAL(8), INTENT(IN)    :: R
    REAL, INTENT(IN)    :: uthres

    REAL(8)    :: den(nmx), diam(nmx)
    REAL(8)    :: dvol(nmx), distr_dust(nmx), dlndp(nmx)
    REAL(8)    :: dsurface(nmx), ds_rel(nmx)
    REAL(8)    :: u_ts0, u_ts, dsrc, dmass, dvol_tot
    REAL(8)    :: salt,emit, emit_vol, stotal
    REAL(8)    :: g
    real(8)    :: H 
    INTEGER   :: i, j, n

    ! Sandblasting mass efficiency, beta.
    ! Beta maxes out for clay fractions above 0.2 = betamax.

    REAL(8), PARAMETER :: betamax=5.25E-4
    REAL(8) :: beta
    integer :: styp

    ! Constant of proportionality from Marticorena et al, 1997 (unitless)
    ! Arguably more ~consistent~ fudge than alpha, which has many walnuts
    ! sprinkled throughout the literature. - GC

    REAL(8), PARAMETER :: cmb=1.0
    ! REAL, PARAMETER :: cmb=2.61   ! from White,1979

    ! Parameters used in Kok distribution function. Advise not to play with
    ! these without the expressed written consent of someone who knows what
    ! they're doing. - GC

    REAL(8), PARAMETER :: mmd_dust=3.4D-6  ! median mass diameter (m)
    REAL(8), PARAMETER :: gsd_dust=3.0     ! geom. std deviation
    REAL(8), PARAMETER :: lambda=12.0D-6   ! crack propagation length (m)
    REAL(8), PARAMETER :: cv=12.62D-6      ! normalization constant
    REAL(8), PARAMETER :: kvhmax=2.0e-4
    REAL(8), PARAMETER :: RHOSOIL=2650.

    ! Calculate total dust emission due to saltation of sand sized particles.
    ! Begin by calculating DRY threshold friction velocity (u_ts0).  Next adjust
    ! u_ts0 for moisture to get threshold friction velocity (u_ts). Then
    ! calculate saltation flux (salt) where ustar has exceeded u_ts.  Finally,
    ! calculate total dust emission (tot_emit), taking into account erodibility.

    ! Set DRY threshold friction velocity to input value
    u_ts0 = uthres
    g = g0  ! *1.0E2      ! JLS - Is the correct?
    emit=0.0

    call DustEmissionFENGSHA(smois,massfrac(1),massfrac(3), massfrac(2), &
                                erod, R, airden, ustar, uthres, alpha, gamma,kvhmax, &
                                g, RHOSOIL, salt)

    emit = emit + salt
     ! FENGSHA uses the 13 category soil type from the USDA
     ! call calc_fengsha_styp(massfrac(1),massfrac(3),massfrac(2),styp)
     ! Fengsha uses threshold velocities based on dale gilletes data
     ! call fengsha_utst(styp,uthres,u_ts0)

     ! Friction velocity threshold correction function based on physical
     ! properties related to moisture tension. Soil moisture greater than
     ! dry limit serves to increase threshold friction velocity (making
     ! it more difficult to loft dust). When soil moisture has not reached
     ! dry limit, treat as dry

    dvol_tot=0.
    DO n=1,nmx
       dlndp(n)=LOG(up_dust(n)/lo_dust(n))
       dvol(n)=(2.0*reff_dust(n)/cv)*(1.+ERF(LOG(2.0*reff_dust(n)/mmd_dust)/(SQRT(2.)*LOG(gsd_dust))))*&
            EXP(-(2.0*reff_dust(n)/lambda)**3.0)*dlndp(n)
       dvol_tot=dvol_tot+dvol(n)
    END DO
    DO n=1,nmx
       distr_dust(n)=dvol(n)/dvol_tot
    END DO

    DO n=1,nmx
       ! Calculate total mass emitted
       dsrc = emit*distr_dust(n)*dxy*dt1  ! (kg)
       IF (dsrc < 0.0) dsrc = 0.0

       ! Update dust mixing ratio at first model level.
       tc(n) = tc(n) + dsrc / airmas ! (kg/kg)
       bems(n) = 1.e+9*dsrc/(dxy*dt1) ! diagnostic (ug/m2/s) !lzhang
    END DO
 
  END SUBROUTINE SOURCE_DUST  

  subroutine fengsha_utst(styp,uth, ut)
    implicit none

    integer,                            intent(in)  :: styp
    real(8), dimension(fengsha_maxstypes), intent(in)  :: uth
    real(8),                 intent(out) :: ut
    ut = uth(styp)
!     real (8) :: uth(13) = &
!          (/ 0.08,   & ! Sand          - 1
!          0.20,    & ! Loamy Sand      - 2
!          0.30,    & ! Sandy Loam      - 3
!          0.30,    & ! Silt Loam       - 4
!          0.35,    & ! Silt            - 5
!          0.60,    & ! Loam            - 6
!          0.30,    & ! Sandy Clay Loam - 7
!          0.35,    & ! Silty Clay Loam - 8
!          0.45,    & ! Clay Loam       - 9
!          0.45,    & ! Sandy Clay      - 10
!          0.45,    & ! Silty Clay      - 11
!          0.60,    & ! Clay            - 12
!          9.999 /)   ! Other           - 13
    return
  end subroutine fengsha_utst


  subroutine fengsha_drag(z0,R)
    implicit none

    real, intent(in) :: z0
    real(8), intent(out) :: R
    real(8), parameter :: z0s = 1.0e-04 !Surface roughness for ideal bare surface [m]
    ! ------------------------------------------------------------------------
    ! Function: Calculates the MacKinnon et al. 2004 Drag Partition Correction
    !
    !   R = 1.0 - log(z0 / z0s) / log( 0.7 * (12255./z0s) ** 0.8)
    !
    !--------------------------------------------------------------------------
    ! Drag partition correction. See MacKinnon et al. (2004),
    !     doi:10.1016/j.geomorph.2004.03.009
    R = 1.0 - log(z0 / z0s) / log( 0.7 * (12255./z0s) ** 0.8)

    ! Drag partition correction. See Marticorena et al. (1997),
    !     doi:10.1029/96JD02964
    !R = 1.0 - log(z0 / z0s) / log( 0.7 * (10./z0s) ** 0.8)

    return
  end subroutine fengsha_drag

  subroutine DustEmissionFENGSHA(slc, clay, sand, silt,  &
                                  ssm, rdrag, airdens, ustar, uthrs, alpha, gamma, &
                                  kvhmax, grav, rhop, emissions)

    ! !USES:
    implicit NONE

! !INPUT PARAMETERS:
    REAL(8), intent(in) :: slc      ! liquid water content of soil layer, volumetric fraction [1]
    REAL(8), intent(in) :: clay     ! fractional clay content [1]
    REAL(8), intent(in) :: sand     ! fractional sand content [1]
    REAL(8), intent(in) :: silt     ! fractional silt content [1]
    REAL(8), intent(in) :: ssm      ! erosion map [1]
    REAL(8), intent(in) :: rdrag    ! drag partition [1/m]
    REAL(8), intent(in) :: airdens  ! air density at lowest level [kg/m^3]
    REAL(8), intent(in) :: ustar    ! friction velocity [m/sec]
    REAL, intent(in) :: uthrs    ! threshold velocity [m/2]
    REAL, intent(in) :: alpha    ! scaling factor [1]
    REAL, intent(in) :: gamma    ! scaling factor [1]
    REAL(8), intent(in) :: kvhmax   ! max. vertical to horizontal mass flux ratio [1]
    REAL(8), intent(in) :: grav     ! gravity [m/sec^2]
    REAL(8), intent(in) :: rhop     ! soil class density [kg/m^3]

    ! !OUTPUT PARAMETERS:
    REAL(8), intent(inout) :: emissions ! binned surface emissions [kg/(m^2 sec)]

    ! !DESCRIPTION: Compute dust emissions using NOAA/ARL FENGSHA model
    !
    ! !REVISION HISTORY:
    !
    ! 22Feb2020 B.Baker/NOAA    - Original implementation
    ! 29Mar2021 R.Montuoro/NOAA - Refactored for process library
    ! 09Aug2022 B.Baker/NOAA    - Adapted for CCPP-Physics

    ! !Local Variables
    real(8)                  :: alpha_grav
    real(8)                  :: h
    real(8)                  :: kvh
    real(8)                  :: q
    real(8)                  :: rustar
    real(8)                  :: total_emissions
    real(8)                  :: u_sum, u_thresh

!EOP
!-------------------------------------------------------------------------
!  Begin

!  Initialize emissions
!  --------------------
   emissions = 0.

!  Prepare scaling factor
!  ----------------------
   alpha_grav = alpha / grav

   ! Compute vertical-to-horizontal mass flux ratio
   ! ----------------------------------------------
   kvh = DustFluxV2HRatioMB95(clay, kvhmax)

   ! Compute total emissions
   ! -----------------------
   emissions = alpha_grav * (ssm ** gamma) * airdens * kvh

   !  Compute threshold wind friction velocity using drag partition
   !  -------------------------------------------------------------
   rustar = rdrag * ustar

   !  Now compute size-dependent total emission flux
   !  ----------------------------------------------
   ! Fecan moisture correction
   ! -------------------------
   h = moistureCorrectionFecan(slc, sand, clay, rhop)

   ! Adjust threshold
   ! ----------------
   u_thresh = uthrs * h

   u_sum = rustar + u_thresh

   ! Compute Horizontal Saltation Flux according to Eq (9) in Webb et al. (2020)
   ! ---------------------------------------------------------------------------
   q = max(0., rustar - u_thresh) * u_sum * u_sum

   ! Distribute emissions to bins and convert to mass flux (kg s-1)
   ! --------------------------------------------------------------
   emissions = emissions * q


 end subroutine DustEmissionFENGSHA
!-----------------------------------------------------------------
  real function soilMoistureConvertVol2Grav(vsoil, sandfrac, rhop)

! !USES:
    implicit NONE

! !INPUT PARAMETERS:
    REAL(8), intent(in) :: vsoil       ! volumetric soil moisture fraction [1]
    REAL(8), intent(in) :: sandfrac    ! fractional sand content [1]
    REAL(8), intent(in) :: rhop        ! dry dust density [kg m-3]

! !DESCRIPTION: Convert soil moisture fraction from volumetric to gravimetric.
!
! !REVISION HISTORY:
!
!  02Apr2020, B.Baker/NOAA    - Original implementation
!  01Apr2020, R.Montuoro/NOAA - Adapted for GOCART process library

!  !Local Variables
    real :: vsat

!  !CONSTANTS:
    REAL(8), parameter :: rhow = 1000.    ! density of water [kg m-3]

!EOP
!-------------------------------------------------------------------------
!  Begin...

!  Saturated volumetric water content (sand-dependent) ! [m3 m-3]
    vsat = 0.489 - 0.00126 * ( 100. * sandfrac )

!  Gravimetric soil content
    soilMoistureConvertVol2Grav = vsoil * rhow / (rhop * (1. - vsat))

  end function soilMoistureConvertVol2Grav
!----------------------------------------------------------------
  real function moistureCorrectionFecan(slc, sand, clay, rhop)

! !USES:
    implicit NONE

! !INPUT PARAMETERS:
    REAL(8), intent(in) :: slc     ! liquid water content of top soil layer, volumetric fraction [1]
    REAL(8), intent(in) :: sand    ! fractional sand content [1]
    REAL(8), intent(in) :: clay    ! fractional clay content [1]
    REAL(8), intent(in) :: rhop    ! dry dust density [kg m-3]

! !DESCRIPTION: Compute correction factor to account for Fecal soil moisture
!
! !REVISION HISTORY:
!
!  02Apr2020, B.Baker/NOAA    - Original implementation
!  01Apr2020, R.Montuoro/NOAA - Adapted for GOCART process library

!  !Local Variables
    real :: grvsoilm
    real :: drylimit

!EOP
!---------------------------------------------------------------
!  Begin...


!  Convert soil moisture from volumetric to gravimetric
    grvsoilm = soilMoistureConvertVol2Grav(slc, sand, rhop )

!  Compute fecan dry limit
    drylimit = clay * (14.0 * clay + 17.0)

!  Compute soil moisture correction
    moistureCorrectionFecan = sqrt(1.0 + 1.21 * max(0., grvsoilm - drylimit)**0.68)

  end function moistureCorrectionFecan
!---------------------------------------------------------------
  real function DustFluxV2HRatioMB95(clay, kvhmax)

! !USES:
    implicit NONE

! !INPUT PARAMETERS:
    REAL(8), intent(in) :: clay      ! fractional clay content [1]
    REAL(8), intent(in) :: kvhmax    ! maximum flux ratio [1]

!  !CONSTANTS:
    REAL(8), parameter :: clay_thresh = 0.2    ! clay fraction above which the maximum flux ratio is returned

! !DESCRIPTION: Computes the vertical-to-horizontal dust flux ratio according to
!               B.Marticorena, G.Bergametti, J.Geophys.Res., 100(D8), 164!
!               doi:10.1029/95JD00690
!
! !REVISION HISTORY:
!
! 22Feb2020 B.Baker/NOAA    - Original implementation
! 01Apr2021 R.Montuoro/NOAA - Adapted for GOCART process library
!
!EOP
!-------------------------------------------------------------------------
!  Begin...

    if (clay > clay_thresh) then
       DustFluxV2HRatioMB95 = kvhmax
    else
       DustFluxV2HRatioMB95 = 10.0**(13.4*clay-6.0)
    end if

  end function DustFluxV2HRatioMB95

end module module_gocart_dust_fengsha

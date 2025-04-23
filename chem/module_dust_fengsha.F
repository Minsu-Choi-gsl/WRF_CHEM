!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      module dust_fengsha

!-----------------------------------------------------------------------
! Description:
! * Calculate windblown dust emissions based on the FENGSHA scheme.       

! Revision History:
! 16 Dec 10 J.Young: Adapting Daniel Tong's work on windblown dust
! 11 Nov 11 J.Young: generalizing land use/cover
! 30 Sep 13 J.Young: corrected diag file units description; added snow cover adjustment;
! 19 Nov 16 Y. Tang: revised the codes to be used outside CMAQ.
! 28 Nov 16 D. Tong: made soil type mapping to be LSM specific; 
!	             modified threshold velocities based on a recent reanalysis.
!-----------------------------------------------------------------------

! Number of soil types: For both WRF and MM5-PX met models, there are 16 types;
! the first 12 soil types are used and the rest lumped into Other.
      integer, parameter :: nsltyp = 13

! Variables for the windblown dust diagnostic file:
      logical, save      :: dustem_diag   ! flag for dustemis diagnostic file
      integer, parameter :: fndust_diag = 42 ! number of fixed diagnostic output vars
      integer, save      :: ndust_diag       ! number of diagnostic output vars
												 
! Module shared variables:									 
!      real, allocatable, save :: agland( :,: )   ! agriculture land fraction
!      real, allocatable, save :: wmax  ( :,: )   ! max adsorb water percent
!      real, allocatable, save :: kvh   ( :,: )   ! ratio of vertical flux / horizontal (k factor)
!      real, allocatable, save :: sd_ep ( :,: )   ! soil->dust erodiblity potential
!      real, allocatable, save :: tfb   ( :,: )   ! transport fraction in and below canopy

!      integer, save :: sdate, stime  ! scenario start date & time

! Modified values for 12 NAM soil types; Silt values taken from Fig.1 of Gillette et
! al., 1980. Values of Other are set to zero so no clay is available for fine particle
! emissions.
      real :: clayc( nsltyp ) = & ! the amount of clay in each soil type [Kg/Kg]
     	     (/ 0.03,	&! Sand
     		0.06,	&! Loamy Sand
     		0.10,	&! Sandy Loam
     		0.13,	&! Silt Loam
     		0.06,	&! Silt
     		0.18,	&! Loam
     		0.27,	&! Sandy Clay Loam
     		0.34,	&! Silty Clay Loam
     		0.34,	&! Clay Loam
     		0.42,	&! Sandy Clay
     		0.47,	&! Silty Clay
     		0.58,	&! Clay
     		0.00 /)  ! Other

      real :: siltc( nsltyp ) = & ! the amount of silt in each soil type [Kg/Kg]
     	     (/ 0.05,	&! Sand
     		0.12,	&! Loamy Sand
     		0.32,	&! Sandy Loam
     		0.70,	&! Silt Loam
     		0.84,	&! Silt
     		0.05,	&! Loam
     		0.15,	&! Sandy Clay Loam
     		0.56,	&! Silty Clay Loam
     		0.34,	&! Clay Loam
     		0.06,	&! Sandy Clay
     		0.47,	&! Silty Clay
     		0.20,	&! Clay
     		0.00 /)  ! Other

      real :: sandc( nsltyp ) = & ! the amount of sand in each soil type [Kg/Kg]
     	     (/ 0.92,	&! Sand
     		0.82,	&! Loamy Sand
     		0.58,	&! Sandy Loam
     		0.17,	&! Silt Loam
     		0.10,	&! Silt
     		0.43,	&! Loam
     		0.58,	&! Sandy Clay Loam
     		0.10,	&! Silty Clay Loam
     		0.32,	&! Clay Loam
     		0.52,	&! Sandy Clay
     		0.06,	&! Silty Clay
     		0.22,	&! Clay
     		0.00 /)  ! Other

!      real :: eropot( 3 ) =   &  ! erodible potential of soil components
!     	     (/ 0.08,	&! clay
!     		1.00,	&! silt
!     		0.12 /)  ! sand

      public dust_emis_fengsha
      CONTAINS


         subroutine dust_emis_fengsha(ii,jj,lsm,  ncols, nrows, n_dlcat, uland, lay1hgt, &
	 ustar, z0, wspd10, soimt,poro, sltyp, snocov, lwmask, rn, rc, ladut, dmsk,      &
         dust_em,eropot,massfrac)

! Arguments:

         integer, intent (in ) :: ii,jj
         integer, intent (in ) :: lsm   ! Land Surface model: 1- NOAH; 2- Pleim-Xiu (PX)	 
	 integer, intent (in ) :: ncols ! I grid dimension
	 integer, intent (in ) :: nrows ! J grid dimension
         integer, intent (in ) :: n_dlcat ! number of erodible landuse types, 4 for MODIS_NOAH landuse
	                                  ! 6 Closed Shrublands , 7 Open Shrublands, 
		&			  ! 16 Barren or Sparsely Vegetated, 20 Barren Tundra
        REAL*8, intent (in), DIMENSION (3) :: massfrac
        real,  intent (in) :: uland(ncols,nrows,4), & ! % of 4 land categories : water, vegetation, urban, Forest
	        ladut(ncols,nrows,n_dlcat),  & !  % of certain landuses that contain desert
	        lay1hgt(ncols,nrows),  & ! thickness of the lowest model layer (m)
		ustar(ncols,nrows),    & ! friction velocity (m/s)
		z0(ncols,nrows),       & ! surface roughtness length (m)
		soimt(ncols,nrows),    & ! gravimetric soil moisture (kg/kg)
		poro(ncols,nrows),     & ! 
	        lwmask(ncols,nrows),   & ! land-water mask: 0 for water
		snocov(ncols,nrows),   & ! snocov indicator [1=yes, 0=no]
		wspd10(ncols,nrows),   & ! 10-m wind speed (m/s)
		rn(ncols,nrows),       & ! non-convective precipitation rate (cm/hr)
		rc(ncols,nrows),       & ! convective precipitation rate (cm/hr)
		dmsk(n_dlcat)            ! fraction mask of desert in certain landuses
	 integer,  intent (in) ::   sltyp(ncols,nrows)       ! soil type in 13 categories
	 real, intent (out) :: dust_em(ncols, nrows)  ! total dust emissions [g/m**3/s]
        real, intent (in)  :: eropot( 3 )
! Includes:

! External Functions:

! Parameters:
         real zi                           ! Height of Internal Boundary Layer [m]
         real, parameter :: sc  = 32.0     ! SENSIT-based scaling parameter
         real, parameter :: z0s = 1.0e-05  ! Surface roughness for dust [m],
         real, parameter :: amen = 32.0    ! values adjusted based on Dale's MNORT measurements
         real, parameter :: radg = 122.87  ! rho_air/grav [g s^2/m^4] 122.87->120.8
         real, parameter :: fdsa = 1.0     ! factor of dust source area 1/.1/.01
         real, parameter :: gmad = 0.78    ! Factor to convert geometric to aerodynamic diameter
         real, parameter :: cfac = amen * radg * fdsa * gmad  ! combined factor for qam
                               ! = 3066.8352,250.13874,95.8386 w/ amen=32.0,2.61,1.0 resp.

         real, parameter :: gpkg = 1.0e03     ! g/kg

         character( 16 ) :: pname = 'GET_DUST_EMIS'
         character( 16 ) :: vname
         character( 96 ) :: xmsg
         integer status
         integer c, r, j, m, n, v
         real soilgrv(ncols, nrows)

!         integer, save   :: wstep = 0      ! local write counter
!         integer         :: mdate, mtime   ! diagnostic file write date&time

                       ! automatic arrays
!            allocate ( utstar( ncols,nrows,n_dlcat ), &
!     			utem( ncols,nrows,n_dlcat ), &
!     			 qam( ncols,nrows,n_dlcat ), &
!     			elus( ncols,nrows,n_dlcat ), &
!			wmax (ncols, nrows), tfb (ncols, nrows), &
!			kvh (ncols, nrows), sd_ep (ncols, nrows), stat = status )
         real       :: fruf  ( ncols,nrows )
         real       :: fmoit ( ncols,nrows )      ! factor of soil moisture on u*t
         real       :: ustr  ( ncols,nrows )      ! friction velocity [m/s]
         real       :: tfa   ( ncols,nrows )      ! transport fraction above canopy
         real       :: wrbuf ( ncols,nrows )      ! diagnositc write buffer
         real, dimension ( ncols,nrows,n_dlcat ) :: utstar,utem,qam,elus
         real, dimension ( ncols,nrows ) :: wmax,kvh,tfb,sd_ep
         real, dimension ( n_dlcat ) :: qad
!         real, allocatable, save :: utstar( :,:,: ) ! threshold U* [m/s]
!         real, allocatable, save :: utem  ( :,:,: ) ! U term = U*(U**2-Ut**2) [(m/s)**3]
!         real, allocatable, save :: qam   ( :,:,: ) ! emis for landuse type [g/m**2/s]
!         real, allocatable, save :: elus  ( :,:,: ) ! Erodible landuse percent (0~100)
!         real, allocatable, save :: qad   ( : )     ! temp var for qam
!         real, allocatable, save :: wmax  (:,:)     ! max adsorb water percent
!         real, allocatable, save :: tfb   ( :,: )   ! transport fraction in and below canopy
!         real, allocatable, save :: kvh   ( :,: )   ! ratio of vertical flux / horizontal (k factor)
!         real, allocatable, save :: sd_ep ( :,: )   ! soil->dust erodiblity potential
 
         real       :: sumdfr           ! sum var for desert fraction
         real       :: rlay1hgt         ! reciprocal of layer-1 height [1/m]
         real       :: m3j              ! 3rd moment accumulation (J) mode emis rates [m3/m3/s]
         real       :: m3k              ! 3rd moment coarse mode (K) emis rates [m3/m3/s]

         character( 16 ), save :: rc_name, rn_name    ! new names: RC -> RCA, RN -> RNA
         logical, save :: firstime = .true.

! Soil moisture limit: 13 types and 3 variables, which are:
!    1 - saturation moisture limit, (gravimetric units assumed, Kg/Kg)
             !!!! This is probably volumetric-- NOAH
!    2 - fill capacity, and   <- not used
!    3 - wilting point        <- not used
! Modified values compatiable with both MM5 & NAM.
! Silt values are based on NAM documentation on soil types.
! Other includes all types higher than 12. The values of Other, serving as
! placeholders, are randomly chosen.  Values of Other, however, have no effect
! on dust emissions as the threshold velocity of Other will be high.
! Since only soilml( nsltyp,1 ) is used, set the following:
         real :: soilml1( nsltyp ) = &
     		(/ 0.395,    & ! Sand
     		   0.410,    & ! Loamy Sand
     		   0.435,    & ! Sandy Loam
     		   0.485,    & ! Silt Loam
     		   0.476,    & ! Silt
     		   0.451,    & ! Loam
     		   0.420,    & ! Sandy Clay Loam
     		   0.477,    & ! Silty Clay Loam
     		   0.476,    & ! Clay Loam
     		   0.426,    & ! Sandy Clay
     		   0.482,    & ! Silty Clay
     		   0.482,    & ! Clay
     		   0.482 /)    ! Other
! Threshold Friction Velocity [m/s] for 13 soil types and 3 land use types, which are:
!    1 -shrubland,
!    2 -shrubgrass land, and
!    3 -barren land
! Local (surface) Threshold Friction Velocity from field measurements
!   Ref. Gillette et al., JGR, 1980 for desert soils (J=2).
!   Ref. Gillette, JGR, 1988 for Loose Agr. Soils (J=3).
! The following soil types were not measured for desert land (we chose to use
! agr. data): Sandy Clay Loam, Clay Loam, Sandy Clay, and Silty Clay.
! Modified values compatible with both MM5 & NAM.
! There is no measurement of this value for Silt. The values for Silt are chosen from
! Silty Loam since the soil composition is close.
! Other includes all types higher than 12. The values of Other are too high to allow
! any dust emission.
         real :: uthfv( nsltyp,3 ) = reshape ( &
     		(/ 0.80,   0.42,   0.23,    & ! Sand
     		   1.00,   0.51,   0.24,    & ! Loamy Sand
     		   1.40,   0.66,   0.16,    & ! Sandy Loam
     		   1.70,   0.34,   0.51,    & ! Silt Loam
     		   1.70,   0.34,   0.51,    & ! Silt
     		   1.70,   0.49,   0.40,    & ! Loam
     		   1.70,   0.78,   0.58,    & ! Sandy Clay Loam
     		   1.70,   0.33,   0.32,    & ! Silty Clay Loam
     		   1.70,   0.71,   0.63,    & ! Clay Loam
     		   1.70,   0.71,   0.58,    & ! Sandy Clay
     		   1.70,   0.56,   0.19,    & ! Silty Clay
     		   1.70,   0.78,   0.29,    & ! Clay
     		   9.99,   9.99,   9.99 /), & ! Other
     		(/ nsltyp,3 /), order = (/ 2,1 /) ) ! fill columns first


         real :: agmsk = 0.75   ! agriculture (cropland) fraction


!-----------------------------------------------------------------------

            clayc(1)=massfrac(1)
            siltc(1)=massfrac(2)
            sandc(1)=massfrac(3)
            utstar=0.
            utem=0.
            qam=0.
            qad=0.
            elus=0.
            wmax=0.
            kvh=0.
            tfb=0.
            sd_ep=0.
!         if ( firstime ) then
            firstime = .false.

!            allocate ( utstar( ncols,nrows,n_dlcat ), &
!     			utem( ncols,nrows,n_dlcat ), &
!     			 qam( ncols,nrows,n_dlcat ), &
!     			elus( ncols,nrows,n_dlcat ), &
!			wmax (ncols, nrows), tfb (ncols, nrows), &
!			kvh (ncols, nrows), sd_ep (ncols, nrows), stat = status )
!            if ( status .ne. 0 ) then
!               print*, '*** UTSTAR, UTEM, QAM, ELUS or WMAX memory allocation failed'
!               stop
!            end if
!            allocate ( qad( n_dlcat ), stat = status )
!            if ( status .ne. 0 ) then
!               print*,'*** QAD memory allocation failed'
!               stop
!            end if
            
	    call tfbelow(ncols, nrows, uland, tfb)
	    
         do r = 1, nrows
            do c = 1, ncols
               
!     Get soil type index (LWMASK = 0 for water)
               if ( lwmask( c,r ) .gt. 0.0 ) then
                  j = sltyp( c,r )
!     Map land use types between NOAH and PX land surface model
!     Check FENGSHA documentation if other LSMs will be used.
                  if ( j .gt. 4 .and. lsm .eq. 2) j = j + 1 ! PX combines "silt" with "silt loam"
                  if ( j .gt. 12 ) j=13 ! all other land use types treated as "other" 
!     Calculate maximum amount of the adsorbed water
!     w` = 0.0014(%clay)**2 + 0.17(%clay) - w` in %
!     Fecan et al. 1999 Annales Geophys.,17,144-157.
                  wmax( c,r ) = ( 14.0 * clayc( 1 ) + 17.0 ) * clayc( 1 ) ! [%]
                  
!     log(F/G) = 0.134(%clayX100)-6
!     F/G (vertical to horizontal flux) ratio [1/cm] is a function of clay content 0~20%
!     Marticorena and Bergametti, 1995
                  if ( clayc( 1 ) .le. 0.2 ) then
                     kvh( c,r ) = 10.0 ** ( 13.4 * clayc( 1 ) - 6.0 )
                  else
                     KVH(C, R) = 2.E-4
!kvh( c,r ) = 4.786e-4 ! = 10.0 ** ( 13.4 * 0.2 - 6.0 )
                  end if

!     erodibility potential of soil component
                  sd_ep( c,r ) = clayc( 1 ) * eropot( 1 ) &
                      + siltc( 1 ) * eropot( 2 )  &
                      + sandc( 1 ) * eropot( 3 )
                  sd_ep( c,r ) = sum(eropot)
!               if(ii.eq.25.and.jj.eq.22)write(0,*)'sd_ep ',sd_ep,sum(eropot)
                  
               end if           ! over land
               
            end do              ! c
         end do                 ! r

!         end if

! Convert layer-average U* into surface local U*
!    -> U*t has to be surface values from measurements.
!               if(ii.eq.25.and.jj.eq.22)write(0,*)'sd_ep2 ',lwmask,sd_ep,sum(eropot)
!               if(ii.eq.25.and.jj.eq.22)write(0,*)'sd_ep2 ',eropot,massfrac
!               if(ii.eq.25.and.jj.eq.22)write(0,*)'sd_ep2 ',clayc(1)*eropot(1),siltc(1)*eropot(2),sandc(1)*eropot(3)
!               if(ii.eq.25.and.jj.eq.22)write(0,*)'sd_ep2 ',clayc(1),siltc(1),sandc(1)

         do r = 1, nrows
            do c = 1, ncols
               zi = sc * z0( c,r )
               ustr( c,r ) = ustar( c,r ) * log( zi / z0( c,r ) ) / log( zi / z0s )
            end do
         end do

! calculate transport factor above the canopy
!          WRITE(6,*) 'zlzltfa', tfa
         call tfabove ( ncols, nrows, wspd10, ustar, uland, tfa )

!-----------------------------------------------------------------------
!  formula of Draxler & Gillette (2001) Atmos. Environ.
!      F   =  K A (r/g) U* ( U*^2 - Ut*^2 )
!  where
!      F   = emission flux                         [g/m**2-s]
!      K   = constant 2.0E-04                      [1/m]
!      A   = 0~3.5  mean = 2.8
!      r/g = ratio of air density to gravity=120.8 [g-s**2/m**4]
!      U*  = friction velocity                     [m/s]
!      Ut* = threshold friction velocity           [m/s]
!-----------------------------------------------------------------------


! Initialize windblown dust diagnostics output buffer
         if ( dustem_diag .and. wstep .eq. 0 ) then
            dustbf = 0.0   ! array assignment
         end if

         do r = 1, nrows
         do c = 1, ncols
            dust_em( c,r ) = 0.0
            fmoit( c,r ) = 0.0   ! for diagnostic output visualization
            fruf ( c,r ) = 0.0
            do m = 1, n_dlcat
               utstar( c,r,m ) = 0.0   ! for diagnostic output visualization
               qam   ( c,r,m ) = 0.0
               elus  ( c,r,m ) = 0.0
               utem  ( c,r,m ) = 0.0
            end do
            qad = 0.0   ! array assignment
            rlay1hgt = 1./lay1hgt ( c,r )

            if ( snocov( c,r ) .gt. 0.0 ) cycle   ! grid cell snow covered
! Dust possiblity if not over water and rain < 1/100 in. (1 in. = 2.540 cm)
            if ( ( lwmask( c,r ) .gt. 0.0 ) .and. &
                 ( rn( c,r ) + rc( c,r ) ) .le. 0.0254 ) then
               j = sltyp( c,r )
               if ( j .gt. 13 ) j = 13  

!lzhang        convert the soil moisture to gravimetric soil moisture:
                soilgrv(c,r)=100.*soimt( c,r)/(poro(c,r)*(2.65*(1.-clayc(1))+2.5*clayc(1)))
!                if(ii.eq.25.and.jj.eq.22)write(0,*)'soilgrv,soimt = ',soilgrv(c,r),soimt( c,r)
!               if(ii.eq.25.and.jj.eq.22)write(0,*)'sd_ep ',sd_ep,sum(eropot)

! soil moisture effect on U*t
               !if ( soimt( c,r ) .le. 0.01 * wmax( c,r ) ) then   ! wmax in [%]
            if ( soilgrv( c,r ) .le. wmax( c,r ) ) then   ! wmax in[%]lzhang
                  fmoit( c,r ) = 1.0
               else if ( soilgrv( c,r ) .le. soilml1( j ) ) then  
                  fmoit( c,r ) = sqrt( 1.0 + 1.21 &
                             ! * ( 100.0 * soimt( c,r ) - wmax( c,r ) ) ** 0.68 )
                              * ( soilgrv( c,r ) - wmax( c,r ) ) ** 0.68 ) !lzhang
               end if

! roughness effect on effect on U*t
               fruf( c,r ) = 1.0   ! roughness factor = 1.0 (for now)
               if ( soimt( c,r ) .le. soilml1( j ) ) then
                  
#ifdef verbose_wbdust
                  dryhit = dryhit + 1
#endif

! Adjust U*t based on soil moisture ---
                  do m = 1, n_dlcat     ! desert type landuse category
                     elus( c,r,m ) = ladut( c,r,m ) * dmsk( m )  ! desert land [%]
                  end do

                  do m = 1, n_dlcat   ! desert type & crop landuse category

!                     if ( elus( c,r,m ) .gt. 100.0 .or. elus( c,r,m ) .lt. 0.0 ) then
!                        print*,'invalid elus ',elus( c,r,m ), c, r, m
!			stop
!                     end if
!                     WRITE(6,*) 'zlzlelus', elus( c,r,: )

                     if ( elus( c,r,m ) .gt. 0.0 ) then
                        utstar( c,r,m ) = uthfv( j,3 ) * fmoit( c,r ) * fruf( c,r )
!                     WRITE(6,*) 'zlzlutsar',ustr( c,r ), utstar( c,r,m )
                        if ( ustr( c,r ) .gt. utstar( c,r,m ) ) then
!                     WRITE(6,*) 'zlzlhi'
                           utem( c,r,m ) =  ustr( c,r ) * ( ustr( c,r ) * ustr( c,r ) &
                                        -   utstar( c,r,m ) * utstar( c,r,m ) )
                           qam( c,r,m ) = cfac * kvh( c,r ) * utem( c,r,m ) *  sd_ep( c,r )
                              
                           qad( m ) = qad( m ) + qam( c,r,m ) * rlay1hgt &
                                   * ( elus( c,r,m ) * 0.01 )          
                        end if
                     end if
!                     WRITE(6,*) 'zlzlutem', utem
!                     WRITE(6,*) 'zlzlqad', qad
!                     WRITE(6,*) 'zlzlqam', qam
!                     if ( elus( c,r,m ) .eq. 0.0 .and. qad( m ) .ne. 0.0 ) then
!                        print*,'ERROR: erodible land use = 0, but emissions .ne. 0'
!                        stop
!                     end if

!                    dust_em( c,r ) = dust_em( c,r ) + qam( c,r,m )
                     dust_em( c,r ) = dust_em( c,r ) + qad( m )
!            if(ii.eq.25.and.jj.eq.22)then
!                  write(0,*)'in fengsha',dust_em( c,r ),qam( c,r,m),qad( m )
!            endif

                  end do   ! m landuse

! Dust removal by surrounding vegetation <-???
! Adjust dust emissions for transport factors

                  dust_em( c,r ) = dust_em( c,r ) * tfa( c,r ) * tfb( c,r )
!            if(ii.eq.25.and.jj.eq.22)then
!                  write(0,*)'in fengsha2',dust_em( c,r ),tfa( c,r ),tfb( c,r )
!            endif

               end if   ! if soil moisture
            end if   ! if rain & land

         end do   ! c
         end do   ! r


         end subroutine dust_emis_fengsha



      subroutine tfabove (ncols, nrows, wspd10, ustar, uland, tfa )

!-----------------------------------------------------------------------
! Description:
!   Calculate transport fraction considering removal by flow above canopy
!
! Subroutines and Functions Called: None
!
! Revison History:
!  Shan He at RTP 2003
!  Jun 2009 D. Tong
!  Jan 2011 J. Young: mods for inline wind-blown dust module
!  Jun 2011 J. Young: add reshape to vd initialization for pgf90 compiler
!  Nov 2016 Y. Tang : change to callable subroutine besides CMAQ
!-----------------------------------------------------------------------

      implicit none
 
      integer, intent (in) :: ncols, nrows
      real, intent (in) :: wspd10(ncols, nrows), ustar(ncols, nrows), uland(ncols,nrows,4)
      real, intent( out ) :: tfa  ( ncols, nrows )   ! above canopy transport factor

! Parameters:

! Local Variables:

      integer c, r, i   ! loop indicies
      integer indx      ! windspeed index

! Based on Slinn 1982, fugitive dust in coarse mode, Vd average from
! PM2~PM10; 3 wind velocities: 10, 5, 1 m/s; 4 canopy characteristics:
! water and lamda = 1.0, 3.5, 5.0

      real, parameter :: vd( 4,3 ) = reshape ( &  ! deposition velocity [m/s]
     		 (/ 0.0051,   & ! natural water  Vd(1,1), WindSpeed = 10 m/s
     		    0.0152,   & ! lamda = 2.0	 Vd(2,1), WindSpeed = 10 m/s
     		    0.0268,   & ! lamda = 3.5	 Vd(3,1), WindSpeed = 10 m/s 
     		    0.0382,   & ! lamda = 5.0	 Vd(4,1), WindSpeed = 10 m/s
     		    0.0018,   & ! natural water  Vd(1,2), WindSpeed =  5 m/s
     		    0.0056,   & ! lamda = 2.0	 Vd(2,2), WindSpeed =  5 m/s
     		    0.0099,   & ! lamda = 3.5	 Vd(3,2), WindSpeed =  5 m/s
     		    0.0141,   & ! lamda = 5.0	 Vd(4,2), WindSpeed =  5 m/s
     		    0.0018,   & ! natural water  Vd(1,3), WindSpeed =  1 m/s
     		    0.0020,   & ! lamda = 2.0	 Vd(2,3), WindSpeed =  1 m/s
     		    0.0021,   & ! lamda = 3.5	 Vd(3,3), WindSpeed =  1 m/s
     		    0.0030 /),& ! lamda = 5.0	 Vd(4,3), WindSpeed =  1 m/s
     		(/ 4,3 /), order = (/ 1,2 /) )

      character( 16 ) :: pname = 'tfabove'
 
      real   :: ul( 4 )
      real   :: ku              ! k=0.08u*
      real   :: omeg            ! transport factor variable

! ----------------------------------------------------------------------
      do r = 1, nrows
      do c = 1, ncols
         tfa( c,r ) = 0.0
         if ( wspd10( c,r ) .le. 2.0 ) then        ! surface wind [0-2] range
            indx = 3
         else if ( wspd10( c,r ) .lt. 7.0 ) then   ! surface wind (2-7) range
            indx = 2
         else                                     ! surface wind [7-^) range
            indx = 1
         end if
         ul = uland( c,r,: )   ! array asignment
         ku = 0.08 * ustar( c,r )
         do i = 1, 4                              ! landuse loop
            omeg = 0.0
            if ( ul( i ) .gt. 0.0 ) then          ! non-zero landuse
               omeg = ku / ( vd( i,indx ) + ku )  ! ratio dqup/dqrd    
               tfa( c,r ) = tfa( c,r ) + 0.01 * omeg * ul( i )
            end if
         end do
      end do
      end do

      return
      end subroutine tfabove

      subroutine tfbelow ( ncols, nrows, uland, tfb )

!-----------------------------------------------------------------------
! Description:
!   Extracts selected landuse types from BELD01 and BELD03 and merges
!   the selections into a dust-related landuse array (ULAND).
!   Optionally, reads 3 gridded crop calendar file and calculates an
!   erodible agriculture land fraction.
!   Applies a predetermined removal fraction in and below canopy to
!   ULAND and determines a transport factor (TFB) for this regime.
 
! Subroutines and Functions Called:
!      OPEN3, INTERPX, XTRACT3, M3EXIT

! Revison History:
!  Jun 2009 D. Tong
!  Jan 2011 J. Young: mods for inline wind-blown dust module
!  Apr 2011 S.Roselle: replaced I/O API include files with UTILIO_DEFN
!  Jul 2011 J.Young: pad vnmld for uniform string lengths
!  Nov 2016 Y. Tang: change to callable subroutine used outside of CMAQ
!-----------------------------------------------------------------------

      integer, intent( in )  :: ncols, nrows
      real,    intent( in )  :: uland (ncols,nrows,4)
      real,    intent( out ) :: tfb(ncols, nrows ) ! dust transport factor from
                                           ! flow into canopy
 
      character(  16 ) :: pname = 'TFBELOW'
      character( 128 ) :: xmsg = ' '

      real              :: ufc  ! total weighted removal fraction for each cell

      integer i, c, r

      real   :: fr( 4 ) = & ! removal fraction in and below canopy
     	      (/ 0.0,	  & ! water
     		 0.10,    & ! grass&crop
     		 0.50,    & ! building
     		 0.95 /)    ! forest

      do r = 1, nrows
      do c = 1, ncols

         ufc = 0.0
         do i = 1, 4
            ufc = ufc + fr( i ) * uland( c,r,i )
         end do
         tfb( c,r ) = 1.0 - ufc * 0.01

!         if ( tfb( c,r ) .gt. 1.0 .or. tfb( c,r ) .lt. 0.0 ) then
!            write( 6,'( a, e10.3, a, 2i4 )' ) '*** Transport factor error: ',  &
!                 tfb( c,r ), ' At col, row: ', c, r 
!            stop
!         end if

      end do   ! ncols
      end do   ! nrows

      return
      end subroutine tfbelow

      end module dust_fengsha

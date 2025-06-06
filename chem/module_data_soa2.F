MODULE module_data_soa2
! RAR, 05/19/2013: Modified to handle semi-volatile POA
!
! This module is based on module_data_sorgam.F, it has been updated to use
! for the new SOA scheme - SOA_VBS

      IMPLICIT NONE

 !     LOGICAL :: diagn_opt
      
      INTEGER, PARAMETER :: blksize=1
      INTEGER, PARAMETER :: numcells=1  ! actual number of cells in arrays ( default is 1 in box model)

!mh    ldesn is a number of deposition species, true number of deposited species may be larger since there
!mh    are species which are deposited with the same rate
!      INTEGER, PARAMETER :: ldesn=13

! These are OCV species for POA, ASOA and BSOA, total=15
!      INTEGER, PARAMETER :: lcv_cvpoa=4, lcv_cvopoa=3, lcv_cvasoa=4, lcv_cvbsoa=4
!      PARAMETER (lsv_ocv= lcv_cvpoa + lcv_cvopoa + lcv_cvasoa + lcv_cvbsoa)   ! total number of OCVs
      INTEGER, PARAMETER :: lcv_cvpoa=4, lcv_cvasoa=4, lcv_cvbsoa=4

      INTEGER, PARAMETER :: lsv_ocv= lcv_cvpoa+lcv_cvasoa+lcv_cvbsoa    ! 12 - total number of OCVs

     ! we have 2 type of SOA - anthropogenic and biogenic
      !sm number of SAPRC species variables in CMU lumped partitioning table
      !sm 1=ALK4(hc5), 2=ALK5(hc8), 3=IVOC,  4=OLE1(ol2), 5=OLE2(oli), 6=ARO1(tol)
      !sm 7=AOR2(xyl), 8=ISOP(iso), 9=SESQ(?),10=TERP(alp)
      INTEGER, PARAMETER :: nsaprc=7  ! number of SOA precursor classes, SAPRC99

! RAR: for this option we have 32 advect aerosol spfp25ies (laerdvc)
! Non-advected species are: H2OA, NU3, AC3, sulf, hno3, nh3 and OCVs
!      INTEGER laerdvc, lnonaerdvc, l1ae, laero, imodes, aspec
!      PARAMETER (laerdvc=32, lnonaerdvc=7 + lsv_ocv)
!      PARAMETER (nspcsda= laerdvc + lnonaerdvc)   ! nspcsda=51 cblk species

!      INTEGER l1ae, laero, imodes, aspec, nspcsda
!      PARAMETER (laero=3,imodes=4,aspec=1)

! 32 - advected, 12 - non advected OCV species, 7 - other non-advected aerosol species
! Advected aerosol species: Inorganic species (6), POA species (8), ASOA species (8), BSOA species (6),
!                           Primary EC, P25 (4), Number concentrations (2); Total= 32
!
! Non-advected aerosol species: Water in aeorosl (2), 3rd mode (2), 3 species for inorganic part (vsulf, nitate, nh3);
!                               Total=7
      INTEGER, PARAMETER :: nspcsda=34+7+lsv_ocv+4+4+1+15+2   ! nspcsda= 34 + 7 + 12 + 4 + 4 + 1 = 62 +15 n2o5 = 77  cblk species

      INTEGER, PARAMETER :: Nbin_voc_ox= 4, Nbin_all= 5

      REAL, DIMENSION (Nbin_all) :: poa_vol_frac, vol_zero, vol_svoc, vol_ivoc, &
                                    vol_high, vol_wood, vol_test1, vol_test2

! Volatility bins:     0.1    1    10   100   1000
! Below are the fractions of POA emissions at different volatility bins
! The purpose of these fractions is to take into account different volatility options
! for POA emissions (poa_vol_opt in namelist.input):
! 0) sum is 1.0, non-volatile POA
     DATA vol_zero / 1.00, 0.00, 0.00, 0.00, 0.00 /

! 1) sum is 1.0, semi-volatile POA,  From Robinson et al. (2007)
     DATA vol_svoc / 0.09, 0.09, 0.14, 0.18, 0.50 /

! 2) sum is 2.5, semi-volatile POA and IVOC emissions = x1.5 POA,  from Robinson et al. (2007)
     DATA vol_ivoc / 0.09, 0.09, 0.14, 0.18, 2.00 /

! 3) sum is 7.5, like vol_ivoc, but SVOCs= x3 POA, from Dzepina et al. (2009)
     DATA vol_high / 0.27, 0.27, 0.42, 0.54, 6.00 /

! 4) sum is 1.0, for wood smoke emisisons, from Grieshop et al. (2009)
     DATA vol_wood / 0.00, 0.10, 0.14, 0.33, 0.43 /

! 5) sum is 1.0, partitioned equally between non-volatile and bin-1; For testing purposes only!
     DATA vol_test1 / 0.50, 0.50, 0.00, 0.00, 0.00 /

! 6) sum is 1.0, low volatility test case; For testing purposes only!
     DATA vol_test2 / 0.18, 0.18, 0.28, 0.36, 0.00 /
!-------------------------------------------------------------------------------

!     LAERDVC  number of advected aerosol dynamic parameters for a given
!     component species
!ia     L1AE        advected parameters+non-advected parameters
!ia     LAERO       number of aerosol component species
!ia     imodes      number of aerosol modes
!ia     ASPEC       number of gas phase comp. that are added dynamically
!ia		    currently only sulfate (=1)
!bs
! updated ldrog numbers for the new SOA mechanism
      INTEGER, PARAMETER :: ldroga=4    ! anthropogenic: ALK4,ALK5,IVOC,OLE1,OLE2,ARO1,ARO2
      INTEGER, PARAMETER :: ldrogb=3    ! biogenic: ISOP,SESQ,TERP
      INTEGER, PARAMETER :: ldrogr=1    ! for branching ratio
      INTEGER, PARAMETER :: ldrog0=1    ! for CVASOA0 production
      INTEGER, PARAMETER :: ldrog_soa2=ldroga+ldrogb+ldrogr+ldrog0        !=11, I've renamed this parameter to separate from "ldrog" for MADE/SORGAM

!!!WARNING: The current code assumes that ldrog_vbs is equal to ldrog_vbs in module_soa_vbs, otherwise the chem_driver etc. has to be modified!!!

! *** set up indices for array  CBLK
! Aerosol species in chem (from registry):
! so4aj,so4ai,nh4aj,nh4ai,no3aj,no3ai,poa0j,poa0i,poa1j,poa1i,poa2j,poa2i,poa3j,poa3i,asoa1j,asoa1i,asoa2j,asoa2i,asoa3j,asoa3i,
! bsoa1j,bsoa1i,bsoa2j,bsoa2i,bsoa3j,bsoa3i,ecj,eci,p25j,p25i,nu0,ac0

! Inorganic species:
      INTEGER, PARAMETER ::  vso4aj= 1
      INTEGER, PARAMETER ::  vso4ai= 2

      INTEGER, PARAMETER ::  vnh4aj= 3
      INTEGER, PARAMETER ::  vnh4ai= 4

      INTEGER, PARAMETER ::  vno3aj= 5
      INTEGER, PARAMETER ::  vno3ai= 6

      INTEGER, PARAMETER ::  vnaaj= 7
      INTEGER, PARAMETER ::  vnaai= 8
      INTEGER, PARAMETER ::  vclaj= 9
      INTEGER, PARAMETER ::  vclai= 10


! indices for accumulation and aitken modes of semi-volatile POA:
      INTEGER, PARAMETER ::  vpoa0i= 11
      INTEGER, PARAMETER ::  vpoa0j= 12

      INTEGER, PARAMETER ::  vpoa1i= 13
      INTEGER, PARAMETER ::  vpoa1j= 14

      INTEGER, PARAMETER ::  vpoa2i= 15
      INTEGER, PARAMETER ::  vpoa2j= 16

      INTEGER, PARAMETER ::  vpoa3i= 17
      INTEGER, PARAMETER ::  vpoa3j= 18

! indices for accumulation and aitken modes of anthropogenic SOA
      INTEGER, PARAMETER ::  vasoa0i= 19
      INTEGER, PARAMETER ::  vasoa0j= 20

      INTEGER, PARAMETER ::  vasoa1i= 21
      INTEGER, PARAMETER ::  vasoa1j= 22

      INTEGER, PARAMETER ::  vasoa2i= 23
      INTEGER, PARAMETER ::  vasoa2j= 24

      INTEGER, PARAMETER ::  vasoa3i= 25
      INTEGER, PARAMETER ::  vasoa3j= 26

! indices for accumulation and aitken modes of biogenic SOA
      INTEGER, PARAMETER ::  vbsoa1j= 27
      INTEGER, PARAMETER ::  vbsoa1i= 28

      INTEGER, PARAMETER ::  vbsoa2j= 29
      INTEGER, PARAMETER ::  vbsoa2i= 30

      INTEGER, PARAMETER ::  vbsoa3j= 31
      INTEGER, PARAMETER ::  vbsoa3i= 32

! Elemental carbon
      INTEGER, PARAMETER ::  vecj   = 33
      INTEGER, PARAMETER ::  veci   = 34

      INTEGER, PARAMETER :: vpols  = 35

      INTEGER, PARAMETER ::  vp25aj = 36
      INTEGER, PARAMETER ::  vp25ai = 37
!-------------------------------------------------------------------------------

      INTEGER, PARAMETER :: vantha = 38
! Seasalt
      INTEGER, PARAMETER :: vseas  = 39
      
      INTEGER, PARAMETER :: vpolp  = 40
! soil
      INTEGER, PARAMETER :: vsoila = 41       

! Aitken and accumulation mode numbers
      INTEGER, PARAMETER ::  vnu0= 42
      INTEGER, PARAMETER ::  vac0= 43
      INTEGER, PARAMETER ::  vcorn= 44
! aerosol water concentration, not advected
      INTEGER, PARAMETER ::  vh2oaj= 45
      INTEGER, PARAMETER ::  vh2oai= 46

! index for Aitken mode 3'rd moment, not advected
      INTEGER, PARAMETER ::  vnu3  = 47

! index for Accumulation mode 3'rd moment, not advected
      INTEGER, PARAMETER ::  vac3 = 48
     
! -- coarse
     INTEGER, PARAMETER :: vcor3 = 49

! index for sulfuric acid vapor concentration
      INTEGER, PARAMETER ::  vsulf= 50

! index for nitric acid vapor concentration
      INTEGER, PARAMETER ::  vhno3= 51

! index for ammonia gas concentration
      INTEGER, PARAMETER ::  vnh3 = 52

! OCVs from POA:
      INTEGER, PARAMETER ::  vcvpoa1= 53
      INTEGER, PARAMETER ::  vcvpoa2= 54
      INTEGER, PARAMETER ::  vcvpoa3= 55
      INTEGER, PARAMETER ::  vcvpoa4= 56

! Anthropogenic OCVs:
      INTEGER, PARAMETER ::  vcvasoa1= 57
      INTEGER, PARAMETER ::  vcvasoa2= 58
      INTEGER, PARAMETER ::  vcvasoa3= 59
      INTEGER, PARAMETER ::  vcvasoa4= 60
! Biogenic OCVs:
      INTEGER, PARAMETER ::  vcvbsoa1= 61
      INTEGER, PARAMETER ::  vcvbsoa2= 62
      INTEGER, PARAMETER ::  vcvbsoa3= 63
      INTEGER, PARAMETER ::  vcvbsoa4= 64
! For N2O5 HET / ISORROPIA
      INTEGER, PARAMETER ::  vn2o5=    65
      INTEGER, PARAMETER ::  vclno2=   66
      INTEGER, PARAMETER ::  vgamn2o5= 67
      INTEGER, PARAMETER ::  vcn2o5=   68
      INTEGER, PARAMETER ::  vkn2o5=   69
      INTEGER, PARAMETER ::  vyclno2=  70
      INTEGER, PARAMETER ::  vsnu=     71
      INTEGER, PARAMETER ::  vsac=     72
      INTEGER, PARAMETER ::  vhcl=     73
      !turned off Minsu choi
!      INTEGER, PARAMETER ::  vcaaj=    74
!      INTEGER, PARAMETER ::  vcaai=    75
!      INTEGER, PARAMETER ::  vkaj=     76
!      INTEGER, PARAMETER ::  vkai=     77
!      INTEGER, PARAMETER ::  vmgaj=    78
!      INTEGER, PARAMETER ::  vmgai=    79

!-----------------------------------------------------------------------------
! RAR: I left deposition velocities here as before, to keep VDVG_2 intact
! All the species in this section are for deposition
! *** set up species dimension and indices for sedimentation
!     velocity array VSED

! number of sedimentation velocities         
      INTEGER, PARAMETER :: naspcssed=6

! index for Aitken mode number                  
      INTEGER, PARAMETER :: vsnnuc=1

! index for Accumulation mode number            
      INTEGER, PARAMETER :: vsnacc=2

! index for coarse mode number                  
      INTEGER, PARAMETER :: vsncor=3

! index for Aitken mode mass                     
      INTEGER, PARAMETER :: vsmnuc=4

! index for accumulation mode mass               
      INTEGER, PARAMETER :: vsmacc=5

! index for coarse mass                         
      INTEGER, PARAMETER :: vsmcor=6

! *** set up species dimension and indices for deposition
!     velocity array VDEP

! number of deposition velocities            
      INTEGER, PARAMETER :: naspcsdep=7

! index for Aitken mode number                  
      INTEGER, PARAMETER :: vdnnuc=1

! index for accumulation mode number            
      INTEGER, PARAMETER :: vdnacc=2

! index for coarse mode number                  
      INTEGER, PARAMETER :: vdncor=3

! index for Aitken mode mass                    
      INTEGER, PARAMETER :: vdmnuc=4

! index for accumulation mode                   
      INTEGER, PARAMETER :: vdmacc=5

! index for fine mode mass (Aitken + accumulation)
      INTEGER, PARAMETER :: vdmfine=6

! index for coarse mode mass                    
      INTEGER, PARAMETER :: vdmcor=7
!------------------------------------------------------------------------------

    ! SOA precursors + OH, O3, NO3
    ! anthropogenic
    INTEGER, PARAMETER :: psoaalk=1 
    INTEGER, PARAMETER :: pole=2
    INTEGER, PARAMETER :: paro=3
    INTEGER, PARAMETER :: pivoc=4
!    INTEGER, PARAMETER :: palk5=2
!    INTEGER, PARAMETER :: pivoc=3
!    INTEGER, PARAMETER :: pole1=4
!    INTEGER, PARAMETER :: pole2=5
!    INTEGER, PARAMETER :: paro1=6
!    INTEGER, PARAMETER :: paro2=7

    ! biogenic
    INTEGER, PARAMETER :: pisop=5
    INTEGER, PARAMETER :: pterp=6
    INTEGER, PARAMETER :: psesq=7

    ! for branching
    INTEGER, PARAMETER :: pbrch=8

    ! for CVASOA0
    INTEGER, PARAMETER :: pocv0=9

    ! new indices for OCVs
    INTEGER, PARAMETER :: ipoa1=1
    INTEGER, PARAMETER :: ipoa2=2
    INTEGER, PARAMETER :: ipoa3=3
    INTEGER, PARAMETER :: ipoa4=4

    INTEGER, PARAMETER :: iasoa1=5
    INTEGER, PARAMETER :: iasoa2=6
    INTEGER, PARAMETER :: iasoa3=7
    INTEGER, PARAMETER :: iasoa4=8

    INTEGER, PARAMETER :: ibsoa1=9
    INTEGER, PARAMETER :: ibsoa2=10
    INTEGER, PARAMETER :: ibsoa3=11
    INTEGER, PARAMETER :: ibsoa4=12
!-----------------------------------------------
!  CONTAINS: Fundamental constants for air quality modeling

!  DEPENDENT UPON:  none

!  REVISION HISTORY:

!    Adapted 6/92 by CJC from ROM's PI.EXT.

!    Revised 3/1/93 John McHenry to include constants needed by
!    LCM aqueous chemistry
!    Revised 9/93 by John McHenry to include additional constants
!    needed for FMEM clouds and aqueous chemistry

!    Revised 3/4/96 by Dr. Francis S. Binkowski to reflect current
!    Models3 view that MKS units should be used wherever possible,
!    and that sources be documentated. Some variables have been added
!    names changed, and values revised.

!    Revised 3/7/96 to have universal gas constant input and compute
!    gas constant is chemical form. TWOPI is now calculated rather than

!    Revised 3/13/96 to group declarations and parameter statements.

!    Revised 9/13/96 to include more physical constants.
!    Revised 12/24/96 eliminate silly EPSILON, AMISS

!    Revised 1/06/97 to eliminate most derived constants
!    10/12/11- Modified to use with soa_vbs, by Ravan Ahmadov

! FSB REFERENCES:
!      CRC76,        CRC Handbook of Chemistry and Physics (76th Ed),
!                     CRC Press, 1995
!      Hobbs, P.V.   Basic Physical Chemistry for the Atmospheric Scien
!                     Cambridge Univ. Press, 206 pp, 1995.
!      Snyder, J.P., Map Projections-A Working Manual, U.S. Geological
!                     Paper 1395 U.S.GPO, Washington, DC, 1987.
!      Stull, R. B., An Introduction to Bounday Layer Meteorology, Klu
!                     Dordrecht, 1988

! Geometric Constants:
      REAL*8, PARAMETER :: pirs=3.14159265358979324

! Fundamental Constants: ( Source: CRC76, pp 1-1 to 1-6)

! Avogadro's Constant [ 1/mol ]
      REAL, PARAMETER :: avo=6.0221367E23

! universal gas constant [ J/mol-K ]
      REAL, PARAMETER :: rgasuniv=8.314510

! standard atmosphere  [ Pa ]
      REAL, PARAMETER :: stdatmpa=101325.0

! Standard Temperature [ K ]
      REAL, PARAMETER :: stdtemp=273.15

! Stefan-Boltzmann [ W/(m**2 K**4) ]
!      REAL, PARAMETER :: stfblz=5.67051E-8

! mean gravitational acceleration [ m/sec**2 ]
      REAL, PARAMETER :: grav=9.80622

! FSB Non MKS quantities:
! Molar volume at STP [ L/mol ] Non MKS units
!      REAL, PARAMETER :: molvol=22.41410

! molecular weights [g/mol] for SO4, HNO3, NH3, NH4
      REAL, PARAMETER :: mwso4=96.0576, mwhno3=63.01287, &
                         mwnh3=17.03061, mwnh4=18.03858, mwno3=62.0049

! molecular weight for HSO4
!      REAL mwhso4
!      PARAMETER (mwhso4=mwso4+1.0080)

! Atmospheric Constants:      78.06%  N2, 21% O2 and 0.943% A on a mole
      REAL, PARAMETER :: mwair=28.9628          ! mean molecular weight for dry air [ g/mol ]

! dry-air gas constant [ J / kg-K ]
      REAL, PARAMETER :: rdgas=1.0E3*rgasuniv/mwair

!  3*PI
      REAL, PARAMETER :: threepi=3.0*pirs

!  6/PI
      REAL, PARAMETER :: f6dpi=6.0/pirs

!  1.0e9 * 6/PIRS
      REAL, PARAMETER :: f6dpi9=1.0E9*f6dpi

! 1.0e-9 * 6/PIRS
      REAL, PARAMETER :: f6dpim9=1.0E-9*f6dpi

!  SQRT( PI )
      REAL, PARAMETER :: sqrtpi=1.7724539

!  SQRT( 2 )
      REAL, PARAMETER :: sqrt2=1.4142135623731

!  ln( sqrt( 2 ) )
      REAL, PARAMETER :: lgsqt2=0.34657359027997

!  1/ln( sqrt( 2 ) )
      REAL, PARAMETER :: dlgsqt2=1.0/lgsqt2

!  1/3 and 2/3
      REAL, PARAMETER :: one3=1.0/3.0, two3=2.0/3.0

! *** physical constants:

! Boltzmann's Constant [ J / K ]
      REAL, PARAMETER :: boltz=rgasuniv/avo

! *** component densities [ kg/m**3 ]:
!  bulk density of aerosol sulfate
      REAL, PARAMETER :: rhoso4=1.8E3

!  bulk density of aerosol ammonium
      REAL, PARAMETER :: rhonh4=1.8E3

! bulk density of aerosol nitrate
      REAL, PARAMETER :: rhono3=1.8E3

!  bulk density of aerosol water
      REAL, PARAMETER :: rhoh2o=1.0E3

! bulk density for aerosol organics
      REAL, PARAMETER :: rhoorg=1.4E3        ! RAR: Check this later

! OM:OC ratio for POA emissions
!      REAL, PARAMETER :: rat_om_oc=1.25     ! RAR: not needed, NEI emissions come in ug/m3

! bulk density for aerosol soil dust
      REAL, PARAMETER :: rhosoil=2.6E3

! bulk density for marine aerosol
      REAL, PARAMETER :: rhoseas=2.2E3

! bulk density for anthropogenic aerosol
      REAL, PARAMETER :: rhoanth=2.2E3

! bulk density for biogenic aerosol, pollen primary    !ts
      REAL, PARAMETER :: rhopolp=1.5E3                       !ts

! bulk density for biogenic aerosol, pollen secondary  !ts
      REAL, PARAMETER :: rhopols=1.2E3                       !ts

! bulk density of aerosol sodium
      REAL, PARAMETER :: rhona=2.2E3

! bulk density of aerosol chloride
      REAL rhocl
      PARAMETER (rhocl=2.2E3)

! bulk density of aerosol casium. equals to dust.
                REAL rhoca
                PARAMETER (rhoca=2.6E3)

! bulk density of aerosol k. equals to dust.
                REAL rhok
                PARAMETER (rhok=2.6E3)

! bulk density of aerosol mg. equals to dust.
                REAL rhomg
                PARAMETER (rhomg=2.6E3)

! *** Factors for converting aerosol mass concentration [ ug m**-3] to
!      	  to 3rd moment concentration [ m**3 m^-3]

      REAL, PARAMETER :: so4fac=f6dpim9/rhoso4

      REAL, PARAMETER :: nh4fac=f6dpim9/rhonh4

      REAL, PARAMETER :: h2ofac=f6dpim9/rhoh2o

      REAL, PARAMETER :: no3fac=f6dpim9/rhono3

      REAL, PARAMETER :: orgfac=f6dpim9/rhoorg

      REAL soilfac
      PARAMETER (soilfac=f6dpim9/rhosoil)
!
      REAL seasfac
      PARAMETER (seasfac=f6dpim9/rhoseas)

      REAL polpfac                         !ts
      PARAMETER (polpfac=f6dpim9/rhopolp)  !ts

      REAL polsfac                         !ts
      PARAMETER (polsfac=f6dpim9/rhopols)  !ts

!
      REAL, PARAMETER :: anthfac=f6dpim9/rhoanth

      REAL nafac
      PARAMETER (nafac=f6dpim9/rhona)

      REAL clfac
      PARAMETER (clfac=f6dpim9/rhocl)

                REAL cafac
                PARAMETER (cafac=f6dpim9/rhoca)

                REAL kfac
                PARAMETER (kfac=f6dpim9/rhok)

                REAL mgfac
                PARAMETER (mgfac=f6dpim9/rhomg)

!  starting standard surface pressure [ Pa ]  
      REAL, PARAMETER :: pss0=101325.0

!  starting standard surface temperature [ K ]
      REAL, PARAMETER :: tss0=288.15

!  initial sigma-G for nucleimode                 
      REAL, PARAMETER :: sginin=1.70

!  initial sigma-G for accumulation mode          
      REAL, PARAMETER :: sginia=2.00

! initial sigma-G for coarse mode               
      REAL, PARAMETER :: sginic=2.5

!  initial mean diameter for nuclei mode [ m ]    
      REAL, PARAMETER :: dginin=0.01E-6

!  initial mean diameter for accumulation mode [ m ]
      REAL, PARAMETER :: dginia=0.07E-6

! initial mean diameter for coarse mode [ m ]  
      REAL, PARAMETER :: dginic=1.0E-6

!     LOGICAL diagnostics
! *** Scalar variables for fixed standard deviations.

! Flag for writing diagnostics to file       
! nuclei mode exp( log^2( sigmag )/8 )  
      REAL en1
! accumulation mode exp( log^2( sigmag )
      REAL ea1

      REAL ec1
! coarse mode exp( log^2( sigmag )/8 )  
! nuclei        **4                    
      REAL esn04
! accumulation                         
      REAL esa04

      REAL esc04
! coarse                               
! nuclei        **5                    
      REAL esn05

      REAL esa05
! accumulation                         
! nuclei        **8                    
      REAL esn08
! accumulation                         
      REAL esa08

      REAL esc08
! coarse                               
! nuclei        **9                    
      REAL esn09

      REAL esa09
! accumulation                         
! nuclei        **12                   
      REAL esn12
! accumulation                         
      REAL esa12

      REAL esc12
! coarse mode                          
! nuclei        **16                   
      REAL esn16
! accumulation                         
      REAL esa16

      REAL esc16
! coarse                               
! nuclei        **20                   
      REAL esn20
! accumulation                         
      REAL esa20

      REAL esc20
! coarse                               
! nuclei        **25                   
      REAL esn25

      REAL esa25
! accumulation                         
! nuclei        **24                   
      REAL esn24
! accumulation                         
      REAL esa24

      REAL esc24
! coarse                               
! nuclei        **28                   
      REAL esn28
! accumulation                         
      REAL esa28

      REAL esc28
! coarse                               
! nuclei        **32                   
      REAL esn32
! accumulation                         
      REAL esa32

      REAL esc32
! coarese                              
! nuclei        **36                   
      REAL esn36
! accumulation                         
      REAL esa36

      REAL esc36
! coarse                               
! nuclei        **49                   
      REAL esn49

      REAL esa49
! accumulation                         
! nuclei        **52                   
      REAL esn52

      REAL esa52
! accumulation                         
! nuclei        **64                   
      REAL esn64
! accumulation                         
      REAL esa64

      REAL esc64
! coarse                               

      REAL esn100
! nuclei        **100                  
! nuclei        **(-20)                
      REAL esnm20
! accumulation                         
      REAL esam20

      REAL escm20
! coarse                               
! nuclei        **(-32)                
      REAL esnm32
! accumulation                         
      REAL esam32

      REAL escm32
! coarse                               
! log(sginin)                           
      REAL xxlsgn
! log(sginia)                           
      REAL xxlsga

      REAL xxlsgc
! log(sginic )                          
! log(sginin ) ** 2                           
      REAL l2sginin
! log(sginia ) ** 2                           
      REAL l2sginia

      REAL l2sginic

! *** set up COMMON blocks for esg's:

! *** SET NUCLEATION FLAG:
      INTEGER, PARAMETER :: inucl=2      ! INUCL = 0, Kerminen & Wexler Mechanism
                                         ! INUCL = 1, Youngblood and Kreidenweis mech
                                         ! INUCL = 2, Kulmala et al. mechanism

! *** Set flag for sedimentation velocities:
      LOGICAL, PARAMETER :: icoarse=.FALSE. ! *** END AERO_INTERNAL.EXT
! *** Diameters and standard deviations for emissions
!     the diameters are the volume (mass) geometric mean diameters

! *** Aitken mode:
! special factor to compute mass transfer           
      REAL, PARAMETER :: dgvem_i=0.03E-6 ! [ m ]
      REAL, PARAMETER :: sgem_i=1.7

! *** Accumulation mode:
      REAL, PARAMETER :: dgvem_j=0.3E-6 ! [ m ]
      REAL, PARAMETER :: sgem_j=2.0

! *** Coarse mode
      REAL, PARAMETER :: dgvem_c=6.0E-6 ! [ m ] <<< Corrected 11/19/97
      REAL, PARAMETER :: sgem_c=2.2

! *** factors for getting number emissions rate from mass emissions rate
      REAL factnumn       ! Aitken mode
      REAL factnuma       ! accumulation mode
      REAL factnumc       ! coarse mode

      REAL, PARAMETER :: facatkn_min=0.04
      REAL facacc_min
      PARAMETER (facacc_min=1.0-facatkn_min)
      REAL xxm3
      REAL, PARAMETER ::  conmin= 1.E-16
      REAL, PARAMETER ::  epsilc= 1.E-16
! [ ug/m**3 ] ! changed 1/6/98 
      REAL*8    nummin_i                ! factor to set minimum for Aitken mode number
      REAL*8    nummin_j                ! factor to set minimum for accumulation mode nu
      REAL*8    nummin_c                ! factor to set minimum for coarse mode number

!bs      REAL ALPHSULF ! Accommodation coefficient for sulfuric acid
!bs      PARAMETER ( ALPHSULF = 0.05 ) ! my be set to one in future
!bs
!bs      REAL DIFFSULF ! molecular diffusivity for sulfuric acid [ m**2
!bs      PARAMETER( DIFFSULF = 0.08E-4 ) ! may be changed in future
!bs
!bs * 23/03/99 updates of ALPHSULF and DIFFSULF adopted fro new code fro
!bs * DIFFSULF is calculated from Reid, Prausnitz, and Poling, The prope
!bs * of gases and liquids, 4th edition, McGraw-Hill, 1987, pp 587-588.
!bs * Equation (11-4.4) was used.
!bs * The value is at T = 273.16 K and P = 1.01325E05 Pa
!bs * Temperature dependence is included for DIFFSULF via DIFFCORR (see
!bs
! Accommodation coefficient for sulfuric
      REAL, PARAMETER ::  alphsulf=1.0

! molecular weight for sulfuric acid [ kg/mole ] MKS
!      REAL, PARAMETER ::  mwh2so4=98.07354E-3

!cia corrected error 24/11/97
! molecular diffusivity for sulfuric acid [ m**2 /se
      REAL, PARAMETER ::  diffsulf=9.362223E-06

!bs Accomodation coefficient for organic
      REAL, PARAMETER ::  alphaorg=1.0                                    !bs Kleeman et al. '99 propose alpha
!bs Bowman et al. '97 uses alpha = 1.
!bs mean molecular weight of organics [k
!      REAL mworg
!      PARAMETER (mworg=175.0E-03)
!      REAL, PARAMETER ::  mworg=250.0E-03       ! Molecular weight for POA from Murphy et al., 2009

!bs * DIFFORG is calculated from the same formula as DIFFSULF.
!bs * An average elemental composition of C=8, O=3, N=1, H=17 is asuumed
!bs * to calculate DIFFORG at T = 273.16K and  P = 1.01325E05 Pa.
!bs * Temperature dependence is included below.

!bs molecular diffusivity for organics [
      REAL, PARAMETER ::  difforg=5.151174E-06

! *** CCONC is the factor for near-continuum condensation.
! ccofm * sqrt( ta )                    
      REAL, PARAMETER :: cconc=2.0*pirs*diffsulf

!bs * factor for NC condensation for organics
! [ m**2 / sec ]       
      REAL, PARAMETER :: cconc_org=2.0*pirs*difforg

!bs analogue to CCOFM but for organics  
      REAL ccofm_org

! FSB  CCOFM is  the accommodation coefficient
!      times the mean molecular velocity for h2so4 without the temperatu
!      after some algebra

!bs CCOFM_ORG * sqrt(TA)                
! set to a value below                  
      REAL ccofm
! minimum aerosol sulfate concentration
      REAL, PARAMETER :: aeroconcmin=0.0001

!*******************************************************************
!*                                                                 *
!*  start parameters and variables for aerosol-cloud interactions  *
!*                                                                 *
!*******************************************************************
!
!   maxd_atype = maximum allowable number of aerosol types
!   maxd_asize = maximum allowable number of aerosol size bins
!   maxd_acomp = maximum allowable number of chemical components
!	in each aerosol size bin
!   maxd_aphase = maximum allowable number of aerosol phases (gas, cloud, ice, rain, ...)
!
!   ntype_aer = number of aerosol types
!   nsize_aer(t) = number of aerosol size bins for aerosol type t. each bin w/ same set of components
!   nphase_aer = number of aerosol phases
!
!   msectional - if positive, moving-center sectional code is utilized,
!	and each mode is actually a section.
!   maerosolincw - if positive, both unactivated/interstitial and activated
!	aerosol species are simulated.  if zero/negative, only the
!	unactivated are simulated.
!
!   ncomp_aer(t) = number of chemical components for aerosol type t
!   ncomp_aer_nontracer(t) = number of "non-tracer" chemical components while in gchm code
!   mastercompptr_aer(c,t) = mastercomp type/i.d. for chemical component c
!	(1=sulfate, others to be defined) and aerosol type t.
!   massptr_aer(c,s,t,p) = gchm r-array index for the mixing ratio
!	(moles-x/mole-air) for chemical component c in size bin s for type t and phase p
!
!   waterptr_aer(s,t) = mixing ratio (moles-water/mole-air) for water
!       associated with aerosol size bin s and type t
!   hygroptr_aer(s,t) = gchm r-array index for the bulk hygroscopicity of the size bin and type
!   numptr_aer(s,t,p) = gchm r-array index for the number mixing ratio
!	(particles/mole-air) for aerosol size bin s, type t, and phase p
!       If zero or negative, then number is not being simulated.
!
!   mprognum_aer(s,t,p) - if positive, number mixing-ratio for size s, type t,
!       and phase p will be prognosed.  Otherwise, no.
!
!   ntot_mastercomp_aer = number of aerosol chemical components defined
!   dens_mastercomp_aer(mc) = dry density (g/cm^3) of aerosol master chemical component type c
!   mw_mastercomp_aer(mc) = molecular weight of aerosol master chemical component type mc
!   name_mastercomp_aer(mc) = name of aerosol master chemical component type mc
!   mc=mastercompptr_aer(c,t)
!   dens_aer(c,t) = dry density (g/cm^3) of aerosol chemical component type c and type t
!   mw_aer(c,t) = molecular weight of aerosol chemical component type c and type t
!   name_aer(c,t) = name of aerosol chemical component type c and type t
!
!   lptr_so4_aer(s,t,p) = gchm r-array index for the
!	mixing ratio for sulfate associated with aerosol size bin s, type t, and phase p
!   (similar for msa, oc, bc, nacl, dust)
!
!-----------------------------------------------------------------------
!
!   volumcen_sect(s,t)= volume (cm^3) at center of section m
!   volumlo_sect(s,t) = volume (cm^3) at lower boundary of section m
!   volumhi_sect(s,t) = volume (cm^3) at upper boundary of section m
!
!   dlo_sect(s,t) = diameter (cm) at lower boundary of section m
!   dhi_sect(s,t) = diameter (cm) at upper boundary of section m
!   dcen_sect(s,t) = volume arithmetic-mean diameter (cm) of section m
!	(corresponds to volumcen_sect == 0.5*(volumlo_sect + volumhi_sect)
!
!-----------------------------------------------------------------------
!   nov-04 sg ! replaced amode with aer and expanded aerosol dimension to include type and phase

	integer, parameter :: maxd_atype  = 2
	integer, parameter :: maxd_asize = 2
	integer, parameter :: maxd_acomp = 23
	integer, parameter :: maxd_aphase = 2
	integer, save :: ai_phase ! interstitial phase of aerosol
	integer, save :: cw_phase ! cloud water phase of aerosol
	integer, save :: ci_phase ! cloud ice  phase of aerosol
	integer, save :: cr_phase ! rain  phase of aerosol
	integer, save :: cs_phase ! snow  phase of aerosol
	integer, save :: cg_phase ! graupel phase of aerosol

	integer, save :: ntype_aer = 0 ! number of types
	integer, save :: ntot_mastercomp_aer = 0 ! number of master components
	integer, save :: nphase_aer = 0 ! number of phases

	integer, save ::   &
      	  msectional, maerosolincw,   &
      	  nsize_aer( maxd_atype ),   & ! number of size bins
      	  ncomp_aer( maxd_atype ),   & ! number of chemical components
      	  ncomp_aer_nontracer( maxd_atype ),   &
          mastercompptr_aer(maxd_acomp, maxd_atype), &   !  mastercomp index
      	  massptr_aer( maxd_acomp, maxd_asize, maxd_atype, maxd_aphase ), & ! index for mixing ratio
      	  waterptr_aer( maxd_asize, maxd_atype ), & ! index for aerosol water
      	  hygroptr_aer( maxd_asize, maxd_atype ), & ! index for aerosol hygroscopicity
      	  numptr_aer( maxd_asize, maxd_atype, maxd_aphase ), & ! index for the number mixing ratio
          mprognum_aer(maxd_asize,maxd_atype,maxd_aphase)

	real, save ::   &
          dens_aer( maxd_acomp, maxd_atype ),   &
          dens_mastercomp_aer( maxd_acomp ),   &
      	  mw_mastercomp_aer( maxd_acomp ), &
      	  mw_aer( maxd_acomp, maxd_atype ),  &
      	  hygro_mastercomp_aer( maxd_acomp ), &
      	  hygro_aer( maxd_acomp, maxd_atype )
	character*10, save ::   &
      	  name_mastercomp_aer( maxd_acomp ), &
      	  name_aer( maxd_acomp, maxd_atype )

	real, save ::   &
          volumcen_sect( maxd_asize, maxd_atype ),   &
          volumlo_sect( maxd_asize, maxd_atype ),   &
          volumhi_sect( maxd_asize, maxd_atype ),   &
          dcen_sect( maxd_asize, maxd_atype ),   &
          dlo_sect( maxd_asize, maxd_atype ),   &
          dhi_sect( maxd_asize, maxd_atype ),   &
	  sigmag_aer(maxd_asize, maxd_atype)

	integer, save ::                     &
      	  lptr_so4_aer(maxd_asize,maxd_atype,maxd_aphase),        &
      	  lptr_nh4_aer(maxd_asize,maxd_atype,maxd_aphase),        &
      	  lptr_no3_aer(maxd_asize,maxd_atype,maxd_aphase),        &

! These don't need anyway, no cloud chemistry for the SOA2 option
!      	  lptr_orgpa_aer(maxd_asize,maxd_atype,maxd_aphase),      &
          lptr_poa0_aer(maxd_asize,maxd_atype,maxd_aphase),    &
          lptr_poa1_aer(maxd_asize,maxd_atype,maxd_aphase),    &
          lptr_poa2_aer(maxd_asize,maxd_atype,maxd_aphase),    &
          lptr_poa3_aer(maxd_asize,maxd_atype,maxd_aphase),     &
          lptr_asoa0_aer(maxd_asize,maxd_atype,maxd_aphase),    &
       	  lptr_asoa1_aer(maxd_asize,maxd_atype,maxd_aphase),    &
      	  lptr_asoa2_aer(maxd_asize,maxd_atype,maxd_aphase),    &
      	  lptr_asoa3_aer(maxd_asize,maxd_atype,maxd_aphase),     &
!      	  lptr_asoa4_aer(maxd_asize,maxd_atype,maxd_aphase),     &
      	  lptr_bsoa1_aer(maxd_asize,maxd_atype,maxd_aphase),     &
      	  lptr_bsoa2_aer(maxd_asize,maxd_atype,maxd_aphase),     &
      	  lptr_bsoa3_aer(maxd_asize,maxd_atype,maxd_aphase),     &
!      	  lptr_bsoa4_aer(maxd_asize,maxd_atype,maxd_aphase),     &

      	  lptr_ec_aer(maxd_asize,maxd_atype,maxd_aphase),         &
      	  lptr_p25_aer(maxd_asize,maxd_atype,maxd_aphase),        &
          lptr_anth_aer(maxd_asize,maxd_atype,maxd_aphase),       &
     	  lptr_cl_aer(maxd_asize,maxd_atype,maxd_aphase),         &
      	  lptr_na_aer(maxd_asize,maxd_atype,maxd_aphase),         &
      	  lptr_seas_aer(maxd_asize,maxd_atype,maxd_aphase),       &
     	  lptr_soil_aer(maxd_asize,maxd_atype,maxd_aphase),       &
          lptr_pols_aer(maxd_asize,maxd_atype,maxd_aphase),       &
          lptr_polp_aer(maxd_asize,maxd_atype,maxd_aphase)

	logical, save ::  do_cloudchem_aer(maxd_asize,maxd_atype)

!   molecular weights (g/mol)
	real, parameter :: mw_so4_aer   = 96.066
	real, parameter :: mw_no3_aer   = 62.007
	real, parameter :: mw_nh4_aer   = 18.042
	real, parameter :: mw_oc_aer    = 250.0
	real, parameter :: mw_ec_aer    = 1.0
	real, parameter :: mw_oin_aer   = 1.0
	real, parameter :: mw_dust_aer  = 100.087
	real, parameter :: mw_seas_aer  = 58.440
	real, parameter :: mw_cl_aer    = 35.450
	real, parameter :: mw_na_aer    = 22.990
	real, parameter :: mw_water_aer = 18.016
        real, parameter :: mw_polp_aer  = 80.440 !ts will be changed: arbitrary
        real, parameter :: mw_pols_aer  = 44.440 !ts will be changed: arbitrary 


!   dry densities (g/cm3)
	real, parameter :: dens_so4_aer  = 1.80   ! = rhoso4
	real, parameter :: dens_no3_aer  = 1.80   ! = rhono3
	real, parameter :: dens_nh4_aer  = 1.80   ! = rhonh4
	real, parameter :: dens_oc_aer   = 1.40   ! = rhoorg ! (Ng et al., 2007) RAR
	real, parameter :: dens_ec_aer   = 1.70
	real, parameter :: dens_dust_aer = 2.60  ! = rhosoil
	real, parameter :: dens_oin_aer  = 2.20  ! = rhoanth
	real, parameter :: dens_seas_aer = 2.20  ! = rhoseas
	real, parameter :: dens_cl_aer   = 2.20
	real, parameter :: dens_na_aer   = 2.20
        real, parameter :: dens_polp_aer = 1.50  ! ts= rhopolp
        real, parameter :: dens_pols_aer = 1.20  ! ts= rhopols

!   water density (g/cm3)
	real, parameter :: dens_water_aer  = 1.0

!   hygroscopicity (dimensionless)
	real, parameter :: hygro_so4_aer  = 0.5
	real, parameter :: hygro_no3_aer  = 0.5
	real, parameter :: hygro_nh4_aer  = 0.5
	real, parameter :: hygro_oc_aer   = 0.14
	real, parameter :: hygro_ec_aer   = 1.e-6
	real, parameter :: hygro_oin_aer  = 0.14
	real, parameter :: hygro_dust_aer = 0.1
	real, parameter :: hygro_seas_aer = 1.16
	real, parameter :: hygro_cl_aer   = 1.16
	real, parameter :: hygro_na_aer   = 1.16
        real, parameter :: hygro_polp_aer = 0.25 !same like organics increased a bit for primary
        real, parameter :: hygro_pols_aer = 0.14 !same like organics increased a bit for primary

! table lookup of aerosol impaction/interception scavenging rates
	real dlndg_nimptblgrow
	integer nimptblgrow_mind, nimptblgrow_maxd
	parameter (nimptblgrow_mind=-14, nimptblgrow_maxd=24)
     	real scavimptblnum(4, nimptblgrow_mind:nimptblgrow_maxd, maxd_asize, maxd_atype), &
     	     scavimptblvol(4, nimptblgrow_mind:nimptblgrow_maxd, maxd_asize, maxd_atype)

!SAM 10/08 Gaussian quadrature constants for SOA_VBS deposition numerical integration
      INTEGER, PARAMETER ::  NGAUSdv= 7   ! Number of Gaussian Quadrature Points - constants defined in aerosols_sorgam_init
      REAL Y_GQ(NGAUSdv), WGAUS(NGAUSdv)

!*****************************************************************
!*                                                               *
!*  end parameters and variables for aerosol-cloud interactions  *
!*                                                               *
!*****************************************************************


END Module module_data_soa2

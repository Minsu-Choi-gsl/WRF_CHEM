! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Parameter Module File
! 
! Generated by KPP-2.1 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : redhc_made_soa_vbs_aq_Parameters.f90
! Time                 : Wed Jan 29 03:54:24 2025
! Working directory    : /scratch1/BMC/wrf-chem/minsusee/WRFV4.6/chem/KPP/mechanisms/redhc_made_soa_vbs_aq
! Equation file        : redhc_made_soa_vbs_aq.kpp
! Output root filename : redhc_made_soa_vbs_aq
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE redhc_made_soa_vbs_aq_Parameters

  USE redhc_made_soa_vbs_aq_Precision
  PUBLIC
  SAVE


! NSPEC - Number of chemical species
  INTEGER, PARAMETER :: NSPEC = 50 
! NVAR - Number of Variable species
  INTEGER, PARAMETER :: NVAR = 48 
! NVARACT - Number of Active species
  INTEGER, PARAMETER :: NVARACT = 45 
! NFIX - Number of Fixed species
  INTEGER, PARAMETER :: NFIX = 2 
! NREACT - Number of reactions
  INTEGER, PARAMETER :: NREACT = 96 
! NVARST - Starting of variables in conc. vect.
  INTEGER, PARAMETER :: NVARST = 1 
! NFIXST - Starting of fixed in conc. vect.
  INTEGER, PARAMETER :: NFIXST = 49 
! NONZERO - Number of nonzero entries in Jacobian
  INTEGER, PARAMETER :: NONZERO = 326 
! LU_NONZERO - Number of nonzero entries in LU factoriz. of Jacobian
  INTEGER, PARAMETER :: LU_NONZERO = 343 
! CNVAR - (NVAR+1) Number of elements in compressed row format
  INTEGER, PARAMETER :: CNVAR = 49 
! NLOOKAT - Number of species to look at
  INTEGER, PARAMETER :: NLOOKAT = 0 
! NMONITOR - Number of species to monitor
  INTEGER, PARAMETER :: NMONITOR = 0 
! NMASS - Number of atoms to check mass balance
  INTEGER, PARAMETER :: NMASS = 1 
! PI - Value of pi
  REAL(kind=dp), PARAMETER :: PI = 3.14159265358979 

! Index declaration for variable species in C and VAR
!   VAR(ind_spc) = C(ind_spc)

  INTEGER, PARAMETER :: ind_IVOC = 1 
  INTEGER, PARAMETER :: ind_SOAALK = 2 
  INTEGER, PARAMETER :: ind_CVBSOA1 = 3 
  INTEGER, PARAMETER :: ind_CVBSOA2 = 4 
  INTEGER, PARAMETER :: ind_CVBSOA3 = 5 
  INTEGER, PARAMETER :: ind_CVBSOA4 = 6 
  INTEGER, PARAMETER :: ind_CVASOA0 = 7 
  INTEGER, PARAMETER :: ind_CVASOA1 = 8 
  INTEGER, PARAMETER :: ind_CVASOA2 = 9 
  INTEGER, PARAMETER :: ind_CVASOA3 = 10 
  INTEGER, PARAMETER :: ind_CVASOA4 = 11 
  INTEGER, PARAMETER :: ind_CVPOA4 = 12 
  INTEGER, PARAMETER :: ind_CVPOA3 = 13 
  INTEGER, PARAMETER :: ind_CVPOA2 = 14 
  INTEGER, PARAMETER :: ind_CVPOA1 = 15 
  INTEGER, PARAMETER :: ind_SULF = 16 
  INTEGER, PARAMETER :: ind_SO2 = 17 
  INTEGER, PARAMETER :: ind_NH3 = 18 
  INTEGER, PARAMETER :: ind_H2O2 = 19 
  INTEGER, PARAMETER :: ind_RXPAR_nt = 20 
  INTEGER, PARAMETER :: ind_N2O5 = 21 
  INTEGER, PARAMETER :: ind_TOLUENE = 22 
  INTEGER, PARAMETER :: ind_ROR_nt = 23 
  INTEGER, PARAMETER :: ind_ROOH = 24 
  INTEGER, PARAMETER :: ind_HO2NO2 = 25 
  INTEGER, PARAMETER :: ind_C2H4 = 26 
  INTEGER, PARAMETER :: ind_PAN = 27 
  INTEGER, PARAMETER :: ind_HNO3 = 28 
  INTEGER, PARAMETER :: ind_CO = 29 
  INTEGER, PARAMETER :: ind_CH3COCHO = 30 
  INTEGER, PARAMETER :: ind_MBO = 31 
  INTEGER, PARAMETER :: ind_SESQ = 32 
  INTEGER, PARAMETER :: ind_PAR = 33 
  INTEGER, PARAMETER :: ind_C10H16 = 34 
  INTEGER, PARAMETER :: ind_ONIT = 35 
  INTEGER, PARAMETER :: ind_ISOP = 36 
  INTEGER, PARAMETER :: ind_OLE = 37 
  INTEGER, PARAMETER :: ind_CH3CHO = 38 
  INTEGER, PARAMETER :: ind_XO2_nt = 39 
  INTEGER, PARAMETER :: ind_XO2N_nt = 40 
  INTEGER, PARAMETER :: ind_CH2O = 41 
  INTEGER, PARAMETER :: ind_C2O3_nt = 42 
  INTEGER, PARAMETER :: ind_OH_nt = 43 
  INTEGER, PARAMETER :: ind_O3 = 44 
  INTEGER, PARAMETER :: ind_NO = 45 
  INTEGER, PARAMETER :: ind_NO2 = 46 
  INTEGER, PARAMETER :: ind_NO3 = 47 
  INTEGER, PARAMETER :: ind_HO2_nt = 48 

! Index declaration for fixed species in C
!   C(ind_spc)

  INTEGER, PARAMETER :: ind_H2O = 49 
  INTEGER, PARAMETER :: ind_M = 50 

! Index declaration for fixed species in FIX
!    FIX(indf_spc) = C(ind_spc) = C(NVAR+indf_spc)

  INTEGER, PARAMETER :: indf_H2O = 1 
  INTEGER, PARAMETER :: indf_M = 2 

END MODULE redhc_made_soa_vbs_aq_Parameters


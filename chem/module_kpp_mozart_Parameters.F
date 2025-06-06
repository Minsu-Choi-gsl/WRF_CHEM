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
! File                 : mozart_Parameters.f90
! Time                 : Wed Jan 29 03:52:03 2025
! Working directory    : /scratch1/BMC/wrf-chem/minsusee/WRFV4.6/chem/KPP/mechanisms/mozart
! Equation file        : mozart.kpp
! Output root filename : mozart
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE mozart_Parameters

  USE mozart_Precision
  PUBLIC
  SAVE


! NSPEC - Number of chemical species
  INTEGER, PARAMETER :: NSPEC = 83 
! NVAR - Number of Variable species
  INTEGER, PARAMETER :: NVAR = 81 
! NVARACT - Number of Active species
  INTEGER, PARAMETER :: NVARACT = 80 
! NFIX - Number of Fixed species
  INTEGER, PARAMETER :: NFIX = 2 
! NREACT - Number of reactions
  INTEGER, PARAMETER :: NREACT = 196 
! NVARST - Starting of variables in conc. vect.
  INTEGER, PARAMETER :: NVARST = 1 
! NFIXST - Starting of fixed in conc. vect.
  INTEGER, PARAMETER :: NFIXST = 82 
! NONZERO - Number of nonzero entries in Jacobian
  INTEGER, PARAMETER :: NONZERO = 741 
! LU_NONZERO - Number of nonzero entries in LU factoriz. of Jacobian
  INTEGER, PARAMETER :: LU_NONZERO = 838 
! CNVAR - (NVAR+1) Number of elements in compressed row format
  INTEGER, PARAMETER :: CNVAR = 82 
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

  INTEGER, PARAMETER :: ind_SO4 = 1 
  INTEGER, PARAMETER :: ind_NH3 = 2 
  INTEGER, PARAMETER :: ind_N2O = 3 
  INTEGER, PARAMETER :: ind_C2H6 = 4 
  INTEGER, PARAMETER :: ind_C3H8 = 5 
  INTEGER, PARAMETER :: ind_BIGENE = 6 
  INTEGER, PARAMETER :: ind_BIGALK = 7 
  INTEGER, PARAMETER :: ind_SO2 = 8 
  INTEGER, PARAMETER :: ind_H2O2 = 9 
  INTEGER, PARAMETER :: ind_TOLUENE = 10 
  INTEGER, PARAMETER :: ind_CRESOL = 11 
  INTEGER, PARAMETER :: ind_EO = 12 
  INTEGER, PARAMETER :: ind_N2O5 = 13 
  INTEGER, PARAMETER :: ind_XOH = 14 
  INTEGER, PARAMETER :: ind_XOOH = 15 
  INTEGER, PARAMETER :: ind_DMS = 16 
  INTEGER, PARAMETER :: ind_C2H5OH = 17 
  INTEGER, PARAMETER :: ind_CH3OOH = 18 
  INTEGER, PARAMETER :: ind_H2 = 19 
  INTEGER, PARAMETER :: ind_MACROOH = 20 
  INTEGER, PARAMETER :: ind_C2H5OOH = 21 
  INTEGER, PARAMETER :: ind_C3H7OOH = 22 
  INTEGER, PARAMETER :: ind_ROOH = 23 
  INTEGER, PARAMETER :: ind_ENEO2 = 24 
  INTEGER, PARAMETER :: ind_MEKOOH = 25 
  INTEGER, PARAMETER :: ind_HO2NO2 = 26 
  INTEGER, PARAMETER :: ind_C2H4 = 27 
  INTEGER, PARAMETER :: ind_HYDRALD = 28 
  INTEGER, PARAMETER :: ind_ONIT = 29 
  INTEGER, PARAMETER :: ind_EO2 = 30 
  INTEGER, PARAMETER :: ind_TOLOOH = 31 
  INTEGER, PARAMETER :: ind_CH3COOOH = 32 
  INTEGER, PARAMETER :: ind_MEK = 33 
  INTEGER, PARAMETER :: ind_CH3COOH = 34 
  INTEGER, PARAMETER :: ind_POOH = 35 
  INTEGER, PARAMETER :: ind_PAN = 36 
  INTEGER, PARAMETER :: ind_GLYOXAL = 37 
  INTEGER, PARAMETER :: ind_TERPOOH = 38 
  INTEGER, PARAMETER :: ind_HNO3 = 39 
  INTEGER, PARAMETER :: ind_CH4 = 40 
  INTEGER, PARAMETER :: ind_CH3OH = 41 
  INTEGER, PARAMETER :: ind_ISOPOOH = 42 
  INTEGER, PARAMETER :: ind_MPAN = 43 
  INTEGER, PARAMETER :: ind_TOLO2 = 44 
  INTEGER, PARAMETER :: ind_BIGALD = 45 
  INTEGER, PARAMETER :: ind_C10H16 = 46 
  INTEGER, PARAMETER :: ind_O1D_CB4 = 47 
  INTEGER, PARAMETER :: ind_ALKOOH = 48 
  INTEGER, PARAMETER :: ind_MEKO2 = 49 
  INTEGER, PARAMETER :: ind_PO2 = 50 
  INTEGER, PARAMETER :: ind_CH3COCH3 = 51 
  INTEGER, PARAMETER :: ind_O = 52 
  INTEGER, PARAMETER :: ind_GLYALD = 53 
  INTEGER, PARAMETER :: ind_ALKO2 = 54 
  INTEGER, PARAMETER :: ind_C2H5O2 = 55 
  INTEGER, PARAMETER :: ind_ISOPNO3 = 56 
  INTEGER, PARAMETER :: ind_C3H7O2 = 57 
  INTEGER, PARAMETER :: ind_CO = 58 
  INTEGER, PARAMETER :: ind_ISOP = 59 
  INTEGER, PARAMETER :: ind_ONITR = 60 
  INTEGER, PARAMETER :: ind_C3H6 = 61 
  INTEGER, PARAMETER :: ind_CH3CHO = 62 
  INTEGER, PARAMETER :: ind_RO2 = 63 
  INTEGER, PARAMETER :: ind_TERPO2 = 64 
  INTEGER, PARAMETER :: ind_HYAC = 65 
  INTEGER, PARAMETER :: ind_MACR = 66 
  INTEGER, PARAMETER :: ind_CH3COCHO = 67 
  INTEGER, PARAMETER :: ind_XO2 = 68 
  INTEGER, PARAMETER :: ind_MVK = 69 
  INTEGER, PARAMETER :: ind_MACRO2 = 70 
  INTEGER, PARAMETER :: ind_ISOPO2 = 71 
  INTEGER, PARAMETER :: ind_CH2O = 72 
  INTEGER, PARAMETER :: ind_NO2 = 73 
  INTEGER, PARAMETER :: ind_MCO3 = 74 
  INTEGER, PARAMETER :: ind_OH = 75 
  INTEGER, PARAMETER :: ind_NO = 76 
  INTEGER, PARAMETER :: ind_CH3O2 = 77 
  INTEGER, PARAMETER :: ind_CH3CO3 = 78 
  INTEGER, PARAMETER :: ind_HO2 = 79 
  INTEGER, PARAMETER :: ind_NO3 = 80 
  INTEGER, PARAMETER :: ind_O3 = 81 

! Index declaration for fixed species in C
!   C(ind_spc)

  INTEGER, PARAMETER :: ind_H2O = 82 
  INTEGER, PARAMETER :: ind_M = 83 

! Index declaration for fixed species in FIX
!    FIX(indf_spc) = C(ind_spc) = C(NVAR+indf_spc)

  INTEGER, PARAMETER :: indf_H2O = 1 
  INTEGER, PARAMETER :: indf_M = 2 

END MODULE mozart_Parameters


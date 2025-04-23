MODULE module_made_soa2
! RAR, 05/09/2013: Modified to handle semi-volatile POA; Also, all coarse particles, CLA and NAA are removed.
! There aren't needed for the simulations focused on SOA.
!
! 10/12/2011: This module is a modified version of the "module_aerosols_sorgam.F". The sorgam subroutine
! has been replaced by a new SOA scheme based on the Volatiliry Basis Set (VBS) approach, recent smog chamber yields
! and multi-generational VOC oxidation mechanism (aging) for SOA formation. The SOA_VBS code has been
! developed by Ravan Ahmadov (ravan.ahmadov@noaa.gov) and Stuart McKeen (Stuart.A.McKeen@noaa.gov) at NOAA/ESRL/CSD.
! This module has been coupled to the modified version of RACM_ESRL_KPP gas chemistry mechanism. Major modifications to the gas
! gas chemistry are inclusion of Sesquiterpenes and separation of MBO from OLI.
! Unlike MOSAIC_VBS this option is for modal approach - MADE aerosol scheme
!
! Some references for the SOA_VBS scheme:
! 1) Ahmadov R., McKeen S.A., Robinson A.L., Bahreini R., Middlebrook A., deGouw J., Meagher J., Hsie E.-Y.,
! Edgerton E., Shaw S., Trainer M. (2012), A volatility basis set model for summertime secondary organic aerosols
! over the eastern U.S. in 2006. J. Geophys. Res.,117, D06301, doi:10.1029/2011JD016831.
! 2) Murphy, B. N. and S. N. Pandis (2009). "Simulating the Formation of Semivolatile Primary and Secondary Organic Aerosol
! in a Regional Chemical Transport Model." Environmental Science & Technology 43(13): 4722-4728.
! 3) Donahue, N. M., A. L. Robinson, et al. (2006). "Coupled partitioning, dilution, and chemical aging of semivolatile
! organics." Environmental Science & Technology 40(8): 2635-2643.
!
! A reference for the MADE aerosol parameterization:
! Ackermann, I.J., H. Hass, M. Memmesheimer, A. Ebel, F.S. Binkowski, and U. Shankar (1998),
! Modal aerosol dynamics model for Europe: Development and first applications, Atmos. Environ., 32(17), 2981-2999.
!
!!WARNING! This aerosol option does not support cloud phase aerosol!
!
!!WARNING! The deposition of organic condensable vapours (cvasoa* and cvbsoa*) are highly uncertain due to lack of observations.
! Currently this process is parameterized using modeled dry deposition velocities of HNO3 (multiplied by "depo_fact" for OCVs).
! Paper by Ahmadov et al. (2012) desribes this approach. The default value for "depo_fact" in WRF-CHEM is 0.25.
! A user can set a different value for "depo_fact" in namelist.input.
!
!!WARNING! Another uncertainty is wet removal of OCVs! This is neglected in the current version of the WRF-CHEM code.
!
USE module_state_description
USE module_data_soa2

IMPLICIT NONE
#define cw_species_are_in_registry

CONTAINS

   SUBROUTINE  soa2_driver ( id,ktau,dtstep,t_phy,                 &
                             alt,p_phy,chem,dz8w,rh,z,z_at_w,      &
                             gamn2o5,cn2o5,kn2o5,yclno2,snu,sac,   &
                             h2oaj,h2oai,nu3,ac3,cor3,             &
                             asulf,ahno3,anh3, &
                             vcsulf_old,kemit,                     &
                             vdrog3,brch_ratio,drog,               &
!                             mgai,mgaj,caai,caaj,kai,kaj,          &
                             ids,ide, jds,jde, kds,kde,            &
                             ims,ime, jms,jme, kms,kme,            &
                             its,ite, jts,jte, kts,kte,            &
                             diagn_opt,do_isorropia,do_n2o5_hetchem   )

!   USE module_configure, only: grid_config_rec_type
!   TYPE (grid_config_rec_type), INTENT (in) :: config_flags
LOGICAL, INTENT(IN ) :: diagn_opt,do_isorropia
INTEGER, INTENT(IN ) :: do_n2o5_hetchem

   INTEGER, INTENT(IN )  ::  ids,ide, jds,jde, kds,kde, &
                             ims,ime, jms,jme, kms,kme, &
                             its,ite, jts,jte, kts,kte, &
                             kemit,   id, ktau

!   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),        &
!         INTENT(IN ) ::                                      moist

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),         &
         INTENT(INOUT ) ::                                   chem
!
! following are aerosol arrays that are not advected
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
         INTENT(INOUT ) ::  h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3, &
                            gamn2o5,cn2o5,kn2o5,yclno2,snu,sac !,        &
!                            mgai,mgaj,caai,caaj,kai,kaj           

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
         INTENT(INOUT ) ::  brch_ratio

! Accumulated dROGs (oxidized VOC) for option 109
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_drog ),          &
         INTENT(INOUT ) ::  drog

   REAL, DIMENSION(ims:ime,kms:kme,jms:jme,ldrog_soa2),                &       ! ldrog_vbs=10
         INTENT(IN   ) :: vdrog3

   REAL, DIMENSION (ims:ime,kms:kme,jms:jme), INTENT(IN ) ::       &
         t_phy, alt, p_phy, dz8w, rh, z, z_at_w                                !rh - fractional relative humidity

   REAL, DIMENSION(ims:ime,kms:kme-0,jms:jme), INTENT(IN ) ::   vcsulf_old
   REAL, INTENT(IN) ::   dtstep

   REAL drog_in(ldrog_soa2)    ! anthropogenic AND biogenic organic aerosol precursors [ug m**-3 s**-1]

!      REAL condvap_in(lspcv) ! condensable vapors [ug m**-3]
   REAL convfac,convfac2

 !  number of species (gas + aerosol)
!      INTEGER nspcsda
!      PARAMETER (nspcsda=l1ae) !bs      ! =54 is the size of CBLK
! (internal aerosol dynamics)
!bs # of anth. cond. vapors in SOA_VBS
!      INTEGER nacv
!      PARAMETER (nacv=lcv_asoa) !bs # of anth. cond. vapors in CTM
!bs total # of cond. vapors in SOA_VBS
      INTEGER ncv
!      PARAMETER (ncv=lsv_ocv) ! RAR
!bs total # of cond. vapors in CTM
      REAL cblk(blksize,nspcsda) ! main array of variables
                                   ! particles [ug/m^3/s]
      REAL soilrat_in
                    ! emission rate of soil derived coars
                    ! input HNO3 to CBLK [ug/m^3]
      REAL nitrate_in     ! input NH3 to CBLK  [ug/m^3]
      REAL nh3_in         ! input SO4 vapor    [ug/m^3]
      REAL vsulf_in


      REAL so4rat_in      ! input SO4 formation[ug/m^3/sec]
      REAL epm25i,epm25j 
      REAL epmcoarse(blksize)
                        ! Emission rate of i-mode EC [ug m**-3 s**-1]
      REAL eeci_in      ! Emission rate of j-mode EC [ug m**-3 s**-1]
      REAL eecj_in      ! Emission rate of j-mode org. aerosol [ug m**-
      REAL eorgi_in

      REAL eorgj_in  ! Emission rate of j-mode org. aerosol [ug m**-
      REAL pres      ! pressure in cb
      REAL temp      ! temperature in K
 !     REAL relhum   ! rel. humidity (0,1)
      REAL brrto 

! JLS for het
      REAL hcl_in 

      REAL mwn2o5
      PARAMETER (mwn2o5=108.009)

      REAL mwhcl
      PARAMETER (mwhcl=36.46100)
      
      REAL mwclno2
      PARAMETER (mwclno2=81.458)
 

      REAL :: p(kts:kte),t(kts:kte),rh0(kts:kte)

      INTEGER :: i,j,k,l
! convert advected aerosol variables to ug/m3 from mixing ratio
! they will be converted back at the end of this driver
INTEGER, SAVE :: icall
!
   do l=p_so4aj,num_chem
      do j=jts,jte
         do k=kts,kte
            do i=its,ite
               chem(i,k,j,l)= max(epsilc,chem(i,k,j,l)/alt(i,k,j))  ! to get aerosol concentrations ug/m3, epsilc = 1.E-16
            enddo
         enddo
      enddo
   enddo

   ! Use RH from phys/??? 
      do 100 j=jts,jte
         do 100 i=its,ite
             do k=kts,kte

             ! added here
               t(k) = t_phy(i,k,j)
               p(k) = .001*p_phy(i,k,j)
               rh0(k) = rh(i,k,j)               ! rh0 mush be between 0.1 and 0.95

!               cblk=0.  ! ???

               convfac =  p(k)/rgasuniv/t(k)*1000.    ! = rho/mw
               so4rat_in= ( chem(i,k,j,p_sulf) - vcsulf_old(i,k,j) )/dtstep*convfac*mwso4
               soilrat_in = 0.
               vsulf_in =   max(epsilc,chem(i,k,j,p_sulf)*convfac*mwso4)
               nitrate_in = max(epsilc,chem(i,k,j,p_hno3)*convfac*mwhno3)
               nh3_in =     max(epsilc,chem(i,k,j,p_nh3)*convfac*mwnh3)
        
           ! JLS for het 
!               hcl_in = max(epsilc,chem(i,k,j,p_hcl)*convfac*mwhcl) 
               cblk(1,vn2o5) = max(epsilc,chem(i,k,j,p_n2o5)*convfac*mwn2o5)
!               cblk(1,vclno2) =max(epsilc,chem(i,k,j,p_clno2)*convfac*mwclno2)

! * organic aerosol precursors DeltaROG and SOA production (10 species)
              ! drog_in(palk4) = vdrog3(i,k,j,palk4)
              ! drog_in(palk5) = vdrog3(i,k,j,palk5)
               drog_in(psoaalk) = vdrog3(i,k,j,psoaalk)
               drog_in(pole) = vdrog3(i,k,j,pole)
               drog_in(paro) = vdrog3(i,k,j,paro)
               drog_in(pivoc) = vdrog3(i,k,j,pivoc)
!               drog_in(pole1) = vdrog3(i,k,j,pole1)
!               drog_in(pole2) = vdrog3(i,k,j,pole2)
!               drog_in(paro1) = vdrog3(i,k,j,paro1)
!               drog_in(paro2) = vdrog3(i,k,j,paro2)
               drog_in(pisop) = vdrog3(i,k,j,pisop)
               drog_in(pterp) = vdrog3(i,k,j,pterp)
               drog_in(psesq) = vdrog3(i,k,j,psesq)
               drog_in(pbrch) = vdrog3(i,k,j,pbrch)
               drog_in(pocv0) = vdrog3(i,k,j,pocv0)


               drog(i,k,j,p_drog_soaalk)= drog(i,k,j,p_drog_soaalk) + drog_in(psoaalk)
               drog(i,k,j,p_drog_ole)= drog(i,k,j,p_drog_ole) + drog_in(pole)
               drog(i,k,j,p_drog_aro)= drog(i,k,j,p_drog_aro) + drog_in(paro)
!               drog(i,k,j,p_drog_alk4)= drog(i,k,j,p_drog_alk4) + drog_in(palk4)
!               drog(i,k,j,p_drog_alk5)= drog(i,k,j,p_drog_alk5) + drog_in(palk5)
               drog(i,k,j,p_drog_ivoc)= drog(i,k,j,p_drog_ivoc) + drog_in(pivoc)
!               drog(i,k,j,p_drog_ole1)= drog(i,k,j,p_drog_ole1) + drog_in(pole1)
!               drog(i,k,j,p_drog_ole2)= drog(i,k,j,p_drog_ole2) + drog_in(pole2)
!               drog(i,k,j,p_drog_aro1)= drog(i,k,j,p_drog_aro1) + drog_in(paro1)
!               drog(i,k,j,p_drog_aro2)= drog(i,k,j,p_drog_aro2) + drog_in(paro2)
               drog(i,k,j,p_drog_isop)= drog(i,k,j,p_drog_isop) + drog_in(pisop)
               drog(i,k,j,p_drog_terp)= drog(i,k,j,p_drog_terp) + drog_in(pterp)
               drog(i,k,j,p_drog_sesq)= drog(i,k,j,p_drog_sesq) + drog_in(psesq)

              if (diagn_opt .AND. icall<200 .AND. k<3) then
                  WRITE(6,*) 'soa2_driver: i,k,j ',i,k,j
                  WRITE(6,*) 'soa2_driver: mwso4,mwhno3,mwnh3 ',mwso4,mwhno3,mwnh3
                  WRITE(6,*) 'soa2_driver: num_drog,p_drog_soaalk ', num_drog,p_drog_soaalk
                  WRITE(6,*) 'soa2_driver: t(k),p(k),rh0(k) ',t(k),p(k),rh0(k)
                  WRITE(6,*) 'soa2_driver: drog_in', drog_in
                  WRITE(6,*) 'soa2_driver: drog(i,k,j,:)', drog(i,k,j,:)
              end if

! CBLK initialization starts here:
! Inorganic species (6)
        cblk(1,vso4aj) =  chem(i,k,j,p_so4aj)
        cblk(1,vso4ai) =  chem(i,k,j,p_so4ai)
        cblk(1,vno3aj) =  chem(i,k,j,p_no3aj)
        cblk(1,vno3ai) =  chem(i,k,j,p_no3ai)
        cblk(1,vnh4aj) =  chem(i,k,j,p_nh4aj)
        cblk(1,vnh4ai) =  chem(i,k,j,p_nh4ai)

! POA species (8)
        cblk(1,vpoa0j) =  chem(i,k,j,p_poa0j)
        cblk(1,vpoa0i) =  chem(i,k,j,p_poa0i)
        cblk(1,vpoa1j) =  chem(i,k,j,p_poa1j)
        cblk(1,vpoa1i) =  chem(i,k,j,p_poa1i)
        cblk(1,vpoa2j) =  chem(i,k,j,p_poa2j)
        cblk(1,vpoa2i) =  chem(i,k,j,p_poa2i)
        cblk(1,vpoa3j) =  chem(i,k,j,p_poa3j)
        cblk(1,vpoa3i) =  chem(i,k,j,p_poa3i)

! ASOA species (8)
! ASOA0* species are secondary, but non-volatile
        cblk(1,vasoa0j) =  chem(i,k,j,p_asoa0j)
        cblk(1,vasoa0i) =  chem(i,k,j,p_asoa0i)

! ASOA1_3* species are secondary and semi-volatile
        cblk(1,vasoa1j) =  chem(i,k,j,p_asoa1j)
        cblk(1,vasoa1i) =  chem(i,k,j,p_asoa1i)
        cblk(1,vasoa2j) =  chem(i,k,j,p_asoa2j)
        cblk(1,vasoa2i) =  chem(i,k,j,p_asoa2i)
        cblk(1,vasoa3j) =  chem(i,k,j,p_asoa3j)
        cblk(1,vasoa3i) =  chem(i,k,j,p_asoa3i)

! BSOA species (6)
        cblk(1,vbsoa1j) =  chem(i,k,j,p_bsoa1j)
        cblk(1,vbsoa1i) =  chem(i,k,j,p_bsoa1i)
        cblk(1,vbsoa2j) =  chem(i,k,j,p_bsoa2j)
        cblk(1,vbsoa2i) =  chem(i,k,j,p_bsoa2i)
        cblk(1,vbsoa3j) =  chem(i,k,j,p_bsoa3j)
        cblk(1,vbsoa3i) =  chem(i,k,j,p_bsoa3i)

! Primary EC, P25 (4)
        cblk(1,vecj   ) =  chem(i,k,j,p_ecj)
        cblk(1,veci   ) =  chem(i,k,j,p_eci)
        cblk(1,vp25aj ) =  chem(i,k,j,p_p25j)
        cblk(1,vp25ai ) =  chem(i,k,j,p_p25i)

        cblk(1,VANTHA   ) =   chem(i,k,j,p_antha)
        cblk(1,VSEAS    ) =   chem(i,k,j,p_seas)
        cblk(1,VSOILA   ) =   chem(i,k,j,p_soila)

! Number concentrations (2)
        cblk(1,vnu0  ) =   max(1.e7,chem(i,k,j,p_nu0))
        cblk(1,vac0  ) =   max(1.e7,chem(i,k,j,p_ac0))
        cblk(1,vcorn ) =   chem(i,k,j,p_corn)
! Non-advected part:
! Water in aerosol, (2)
        cblk(1,vh2oaj) =   max(epsilc,h2oaj(i,k,j))
        cblk(1,vh2oai) =   max(epsilc,h2oai(i,k,j))

! 3rd mode
        cblk(1,vnu3  ) =   max(epsilc,nu3(i,k,j))
        cblk(1,vac3  ) =   max(epsilc,ac3(i,k,j))
        cblk(1,VCOR3 ) =   max(epsilc,cor3(i,k,j))

! JLS N2O5
        cblk(1,vgamn2o5) = max(epsilc,gamn2o5(i,k,j))
        cblk(1,vcn2o5)   = max(epsilc,cn2o5(i,k,j))
        cblk(1,vkn2o5)   = max(epsilc,kn2o5(i,k,j))
        cblk(1,vyclno2)  = max(epsilc,yclno2(i,k,j))
        cblk(1,vsnu)     = max(epsilc,snu(i,k,j))
        cblk(1,vsac)     = max(epsilc,sac(i,k,j))

! Primary OCVs, not oxidized (4)
        cblk(1,vcvpoa1) =  chem(i,k,j,p_cvpoa1)
        cblk(1,vcvpoa2) =  chem(i,k,j,p_cvpoa2)
        cblk(1,vcvpoa3) =  chem(i,k,j,p_cvpoa3)
        cblk(1,vcvpoa4) =  chem(i,k,j,p_cvpoa4)

! Secondary, oxidized from CVPOA (3)
!        cblk(1,vcvopoa1)  =   chem(i,k,j,p_cvopoa1)
!        cblk(1,vcvopoa2)  =   chem(i,k,j,p_cvopoa2)
!        cblk(1,vcvopoa3)  =   chem(i,k,j,p_cvopoa3)

! Secondary, oxidized from VOCs (8)
        cblk(1,vcvasoa1)  =   chem(i,k,j,p_cvasoa1)
        cblk(1,vcvasoa2)  =   chem(i,k,j,p_cvasoa2)
        cblk(1,vcvasoa3)  =   chem(i,k,j,p_cvasoa3)
        cblk(1,vcvasoa4)  =   chem(i,k,j,p_cvasoa4)

        cblk(1,vcvbsoa1)  =   chem(i,k,j,p_cvbsoa1)
        cblk(1,vcvbsoa2)  =   chem(i,k,j,p_cvbsoa2)
        cblk(1,vcvbsoa3)  =   chem(i,k,j,p_cvbsoa3)
        cblk(1,vcvbsoa4)  =   chem(i,k,j,p_cvbsoa4)

!rs. nitrate, nh3, sulf (3)
        cblk(1,vsulf)  =   vsulf_in
        cblk(1,vhno3)  =   nitrate_in
        cblk(1,vnh3)   =   nh3_in
! JLS - n2o5 het
!        cblk(1,vhcl)   =   hcl_in          
!      if ( do_isorropia ) then
!        cblk(1,vcaaj) = caaj(i,k,j) ! chem(i,k,j,p_caaj)
!        cblk(1,vcaai) = caai(i,k,j) ! chem(i,k,j,p_caai)
!        cblk(1,vkaj)  = kaj(i,k,j)  ! chem(i,k,j,p_kaj)
!        cblk(1,vkai)  = kai(i,k,j)  ! chem(i,k,j,p_kai)
!        cblk(1,vmgaj) = mgaj(i,k,j) ! chem(i,k,j,p_mgaj)
!        cblk(1,vmgai) = mgai(i,k,j) ! chem(i,k,j,p_mgai)
!      endif

         cblk(1,VNAAJ   ) = chem(i,k,j,p_naaj)
         cblk(1,VNAAI   ) = chem(i,k,j,p_naai)
         cblk(1,VCLAJ   ) = chem(i,k,j,p_claj)
         cblk(1,VCLAI   ) = chem(i,k,j,p_clai)

        cblk(1,VPOLP    ) =   chem(i,k,j,p_polp) !ts for primary pollen
        cblk(1,VPOLS    ) =   chem(i,k,j,p_pols) !ts for secondary pollen

        if (diagn_opt .AND. icall<200 .AND. k<3) then
            WRITE(6,*) 'soa2_driver: chem(i,k,j,p_poa0j), chem(i,k,j,p_cvpoa1)',chem(i,k,j,p_poa0j), chem(i,k,j,p_cvpoa1)
            WRITE(6,*) 'soa2_driver: nspcsda, cblk', nspcsda, cblk
            WRITE(6,*) 'soa2_driver: vsulf_in,nitrate_in,nh3_in', vsulf_in,nitrate_in,nh3_in
            WRITE(6,*) 'soa2_driver: vsulf,vhno3,vnh3', vsulf,vhno3,vnh3
            WRITE(6,*) 'soa2_driver: cblk(1,vecj),cblk(1,veci) ',cblk(1,vecj),cblk(1,veci)
        end if
! Set emissions to zero      !! RAR: Do we need these here???
         epmcoarse(1)     = 0.
         epm25i        = 0.
         epm25j        = 0.
         eeci_in       = 0.
         eecj_in       = 0.
         eorgi_in      = 0.
         eorgj_in      = 0.


! End of the cblk initialization
!-------------------------------------------------------------------------------

! the following operation updates cblk, which includes the vapors and SOA species; condvap_in is removed
      CALL rpmmod3(k,dtstep,10.*p(k),t(k),rh0(k),              &
                   nitrate_in,nh3_in,vsulf_in,so4rat_in,       &
                   drog_in,                                    &
                   eeci_in,eecj_in,eorgi_in,eorgj_in,          &
                   epm25i,epm25j,epmcoarse,soilrat_in,         &
                   cblk,i,j,k,brrto,diagn_opt,do_isorropia,do_n2o5_hetchem )

! calculation of brch_ratio
        brch_ratio(i,k,j)= brrto

        chem(i,k,j,p_so4aj) = cblk(1,vso4aj )
        chem(i,k,j,p_so4ai) = cblk(1,vso4ai )
        chem(i,k,j,p_nh4aj) = cblk(1,vnh4aj )
        chem(i,k,j,p_nh4ai) = cblk(1,vnh4ai )
        chem(i,k,j,p_no3aj) = cblk(1,vno3aj )
        chem(i,k,j,p_no3ai) = cblk(1,vno3ai )

        chem(i,k,j,p_naaj)  = cblk(1,VNAAJ   )
        chem(i,k,j,p_naai)  = cblk(1,VNAAI   )
        chem(i,k,j,p_claj) = cblk(1,VCLAJ   )
        chem(i,k,j,p_clai) = cblk(1,VCLAI   )


        chem(i,k,j,p_poa0j) =   cblk(1,vpoa0j)
        chem(i,k,j,p_poa0i) =   cblk(1,vpoa0i)
        chem(i,k,j,p_poa1j) =   cblk(1,vpoa1j)
        chem(i,k,j,p_poa1i) =   cblk(1,vpoa1i)
        chem(i,k,j,p_poa2j) =   cblk(1,vpoa2j)
        chem(i,k,j,p_poa2i) =   cblk(1,vpoa2i)
        chem(i,k,j,p_poa3j) =   cblk(1,vpoa3j)
        chem(i,k,j,p_poa3i) =   cblk(1,vpoa3i)

        chem(i,k,j,p_asoa0j) =   cblk(1,vasoa0j)
        chem(i,k,j,p_asoa0i) =   cblk(1,vasoa0i)
        chem(i,k,j,p_asoa1j) =   cblk(1,vasoa1j)
        chem(i,k,j,p_asoa1i) =   cblk(1,vasoa1i)
        chem(i,k,j,p_asoa2j) =   cblk(1,vasoa2j)
        chem(i,k,j,p_asoa2i) =   cblk(1,vasoa2i)
        chem(i,k,j,p_asoa3j) =   cblk(1,vasoa3j)
        chem(i,k,j,p_asoa3i) =   cblk(1,vasoa3i)

        chem(i,k,j,p_bsoa1j) =   cblk(1,vbsoa1j)
        chem(i,k,j,p_bsoa1i) =   cblk(1,vbsoa1i)
        chem(i,k,j,p_bsoa2j) =   cblk(1,vbsoa2j)
        chem(i,k,j,p_bsoa2i) =   cblk(1,vbsoa2i)
        chem(i,k,j,p_bsoa3j) =   cblk(1,vbsoa3j)
        chem(i,k,j,p_bsoa3i) =   cblk(1,vbsoa3i)

        chem(i,k,j,p_ecj)    = cblk(1,vecj  )
        chem(i,k,j,p_eci)    = cblk(1,veci  )
        chem(i,k,j,p_p25j)   = cblk(1,vp25aj)
        chem(i,k,j,p_p25i)   = cblk(1,vp25ai)

        chem(i,k,j,p_antha)  = cblk(1,VANTHA   )
        chem(i,k,j,p_seas)   = cblk(1,VSEAS    )
        chem(i,k,j,p_soila)  = cblk(1,VSOILA   )

        chem(i,k,j,p_polp)   = cblk(1,VPOLP    ) !ts
        chem(i,k,j,p_pols)   = cblk(1,VPOLS    ) !ts

        chem(i,k,j,p_nu0)    = max(1.e7,cblk(1,vnu0))
        chem(i,k,j,p_ac0)    = max(1.e7,cblk(1,vac0))
        chem(i,k,j,p_corn) = cblk(1,VCORN    )

        h2oaj(i,k,j) = cblk(1,vh2oaj)
        h2oai(i,k,j) = cblk(1,vh2oai)
        nu3(i,k,j)   = cblk(1,vnu3  )
        ac3(i,k,j)   = cblk(1,vac3  )
        cor3(i,k,j)  = cblk(1,vcor3 )

        chem(i,k,j,p_cvpoa1)= cblk(1,vcvpoa1)
        chem(i,k,j,p_cvpoa2)= cblk(1,vcvpoa2)
        chem(i,k,j,p_cvpoa3)= cblk(1,vcvpoa3)
        chem(i,k,j,p_cvpoa4)= cblk(1,vcvpoa4)

!        chem(i,k,j,p_cvopoa1)= cblk(1,vcvopoa1)
!        chem(i,k,j,p_cvopoa2)= cblk(1,vcvopoa2)
!        chem(i,k,j,p_cvopoa3)= cblk(1,vcvopoa3)

        chem(i,k,j,p_cvasoa1)= cblk(1,vcvasoa1)
        chem(i,k,j,p_cvasoa2)= cblk(1,vcvasoa2)
        chem(i,k,j,p_cvasoa3)= cblk(1,vcvasoa3)
        chem(i,k,j,p_cvasoa4)= cblk(1,vcvasoa4)

        chem(i,k,j,p_cvbsoa1)= cblk(1,vcvbsoa1)
        chem(i,k,j,p_cvbsoa2)= cblk(1,vcvbsoa2)
        chem(i,k,j,p_cvbsoa3)= cblk(1,vcvbsoa3)
        chem(i,k,j,p_cvbsoa4)= cblk(1,vcvbsoa4)

        chem(i,k,j,p_sulf)= max( epsilc,cblk(1,vsulf)/convfac/mwso4 )
        chem(i,k,j,p_hno3)= max( epsilc,cblk(1,vhno3)/convfac/mwhno3 )
        chem(i,k,j,p_nh3) = max( epsilc,cblk(1,vnh3)/convfac/mwnh3 )
! ISORROPIA only
!      if ( do_isorropia ) then
!        caaj(i,k,j) = cblk(1,vcaaj)
!        caai(i,k,j) = cblk(1,vcaai)
!        kaj(i,k,j)  = cblk(1,vkaj)
!        kai(i,k,j)  = cblk(1,vkai)
!        mgaj(i,k,j) = cblk(1,vmgaj)
!        mgai(i,k,j) = cblk(1,vmgai)
!      endif
! JLS - n2o5 het
        gamn2o5(i,k,j)= cblk(1,vgamn2o5)
        cn2o5(i,k,j)  = cblk(1,vcn2o5)
        kn2o5(i,k,j)  = cblk(1,vkn2o5)
        yclno2(i,k,j) = cblk(1,vyclno2)
        snu(i,k,j)    = cblk(1,vsnu)
        sac(i,k,j)    = cblk(1,vsac)
     
      
      !  chem(i,k,j,p_hcl) = max(epsilc,cblk(1,vhcl)/CONVFAC/MWHCL)
        chem(i,k,j,p_n2o5) = max(epsilc,cblk(1,vn2o5)/CONVFAC/MWN2O5)
      !  chem(i,k,j,p_clno2) = max(epsilc,cblk(1,vclno2)/CONVFAC/MWCLNO2)

! RAR: debugging
  if (diagn_opt .AND. icall<200 .AND. k<3) then
     WRITE(6,*) 'soa2_driver, after rpmmod3: cblk', cblk
     WRITE(6,*) 'soa2_driver, after rpmmod3: chem(i,k,j,p_poa0j), chem(i,k,j,p_cvpoa1)',chem(i,k,j,p_poa0j), chem(i,k,j,p_cvpoa1)
     WRITE(6,*) 'soa2_driver, after rpmmod3: chem(i,k,j,p_ecj), chem(i,k,j,p_eci) ', chem(i,k,j,p_ecj), chem(i,k,j,p_eci)
     icall=icall+1
  end if
      enddo          ! k-loop
100  continue ! i,j-loop ends

! convert aerosol variables back to mixing ratio from ug/m3 to ug/kg(dry air)
  do l=p_so4aj,num_chem
     do j=jts,jte
        do k=kts,kte
           do i=its,ite
              chem(i,k,j,l)= max(epsilc,chem(i,k,j,l)*alt(i,k,j))
           enddo
        enddo
     enddo
  enddo

END SUBROUTINE soa2_driver
!-------------------------------------------------------------------------------

SUBROUTINE sum_pm_soa2 ( alt, chem, h2oaj, h2oai, pm1,               &
                         pm2_5_dry, pm2_5_water, pm2_5_dry_ec,pm10,  &
!                         mgai,mgaj,caai,caaj,kai,kaj,                &
                         do_isorropia,                               &
                         ids,ide, jds,jde, kds,kde,                  &
                         ims,ime, jms,jme, kms,kme,                  &
                         its,ite, jts,jte, kts,kte             )

   INTEGER, INTENT(IN   ) ::     ids,ide, jds,jde, kds,kde,       &
                                 ims,ime, jms,jme, kms,kme,       &
                                 its,ite, jts,jte, kts,kte

   LOGICAL, INTENT(IN   ) :: do_isorropia

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),             &
         INTENT(IN ) :: chem

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
         INTENT(IN ) :: alt,h2oaj,h2oai !,mgai,mgaj,caai,caaj,kai,kaj

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
         INTENT(OUT) :: pm2_5_dry,pm2_5_water,pm2_5_dry_ec,pm10,pm1 

   INTEGER :: i,ii,j,jj,k,n
!
! sum up pm2_5 and pm10 output
!
      pm2_5_dry(its:ite, kts:kte, jts:jte)    = 0.
      pm2_5_water(its:ite, kts:kte, jts:jte)  = 0.
      pm2_5_dry_ec(its:ite, kts:kte, jts:jte) = 0.
      do j=jts,jte
         jj=min(jde-1,j)
      do k=kts,kte
      do i=its,ite
         ii=min(ide-1,i)
         do n=p_so4aj,p_p25i
            pm2_5_dry(i,k,j) = pm2_5_dry(i,k,j) + chem(ii,k,jj,n)
         enddo
    
         if ( do_isorropia ) then
            pm2_5_dry(i,k,j) = pm2_5_dry(i,k,j) !+ &
!                               mgai(ii,k,jj) + mgaj(ii,k,jj) + \
!                               caai(ii,k,jj) + caaj(ii,k,jj) + \
!                                kai(ii,k,jj) +  kaj(ii,k,jj)
         endif
    
         pm2_5_dry_ec(i,k,j) = pm2_5_dry_ec(i,k,j) + chem(ii,k,jj,p_ecj) + chem(ii,k,jj,p_eci)
         pm2_5_water(i,k,j)  = pm2_5_water(i,k,j)  + h2oaj(i,k,j) + h2oai(i,k,j)

         !Convert the units from mixing ratio to concentration (ug m^-3)
         pm2_5_dry(i,k,j)    = pm2_5_dry(i,k,j) / alt(ii,k,jj)
         pm2_5_dry_ec(i,k,j) = pm2_5_dry_ec(i,k,j) / alt(ii,k,jj)
         pm2_5_water(i,k,j)  = pm2_5_water(i,k,j) / alt(ii,k,jj)
      enddo
      enddo
      enddo

! Hardcode pm1
      pm1(its:ite, kts:kte, jts:jte) = 0.
      do j=jts,jte
         jj=min(jde-1,j)
      do k=kts,kte
      do i=its,ite
         ii=min(ide-1,i)
         pm1(i,k,j) = chem(ii,k,jj,p_so4ai)  + chem(ii,k,jj,p_nh4ai)  &
                    + chem(ii,k,jj,p_no3ai)  + chem(ii,k,jj,p_naai)   &
                    + chem(ii,k,jj,p_clai)   + chem(ii,k,jj,p_poa0i)  &
                    + chem(ii,k,jj,p_poa1i)  + chem(ii,k,jj,p_poa2i)  &
                    + chem(ii,k,jj,p_poa3i)  + chem(ii,k,jj,p_asoa0i) &
                    + chem(ii,k,jj,p_asoa1i) + chem(ii,k,jj,p_asoa2i) &
                    + chem(ii,k,jj,p_asoa3i) + chem(ii,k,jj,p_bsoa1i) &
                    + chem(ii,k,jj,p_bsoa2i) + chem(ii,k,jj,p_bsoa3i) &
                    + chem(ii,k,jj,p_eci)    + chem(ii,k,jj,p_p25i)
      enddo
      enddo
      enddo
! RAR: isn't used in this option
      do j=jts,jte
         jj=min(jde-1,j)
         do k=kts,kte
            do i=its,ite
               ii=min(ide-1,i)
               pm10(i,k,j) = pm2_5_dry(i,k,j)                       &
                           + ( chem(ii,k,jj,p_antha)               &
                           + chem(ii,k,jj,p_soila)                 &
                           + chem(ii,k,jj,p_seas)                  &
                           + chem(ii,k,jj,p_polp) ) / alt(ii,k,jj)
            enddo
         enddo
      enddo
    END SUBROUTINE sum_pm_soa2
! ///////////////////////////////////////////////////

SUBROUTINE     soa2_depdriver (id,config_flags,ktau,dtstep,                    &
                               ust,t_phy,rmol,znt,pbl,                         &
                               alt,p_phy,chem,dz8w,rh,z,z_at_w,                &
                               h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3,      &
!                               mgai,mgaj,caai,caaj,kai,kaj,                    &
                               aer_res,vgsa,                                   &
                               numaer,                                         &
                               ids,ide, jds,jde, kds,kde,                      &
                               ims,ime, jms,jme, kms,kme,                      &
                               its,ite, jts,jte, kts,kte                       )

   USE module_configure,only:  grid_config_rec_type
   TYPE (grid_config_rec_type) , INTENT (IN) :: config_flags

   INTEGER, INTENT(IN   )    ::       numaer,                    &
                                      ids,ide, jds,jde, kds,kde, &
                                      ims,ime, jms,jme, kms,kme, &
                                      its,ite, jts,jte, kts,kte, &
                                      id,ktau

!   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),        &
!         INTENT(IN ) ::                                   moist
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),         &
         INTENT(INOUT ) ::                                   chem
!
! following are aerosol arrays that are not advected
!
   REAL, DIMENSION( its:ite, jts:jte, numaer ),                       &
         INTENT(INOUT ) ::  vgsa

   REAL, DIMENSION( its:ite, jts:jte ),                       &
         INTENT(INOUT ) ::           aer_res

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                        &
         INTENT(INOUT ) ::   h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3 !, &
!                             mgai,mgaj,caai,caaj,kai,kaj

   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,    &
          INTENT(IN   ) :: t_phy,alt,p_phy,dz8w,rh,z,z_at_w

   REAL,  DIMENSION( ims:ime ,  jms:jme )                  ,    &
          INTENT(IN   ) ::                     ust,rmol, pbl, znt
   REAL,  INTENT(IN   ) ::                                 dtstep
                                                                                               
   REAL convfac,convfac2
!...blksize set to one in column model ciarev02

!      INTEGER, PARAMETER   :: blksize=1

!...number of aerosol species
!  number of species (gas + aerosol)
!      INTEGER nspcsda
!      PARAMETER (nspcsda=l1ae) !bs
! (internal aerosol dynamics)
!bs # of anth. cond. vapors in SOA_VBS
!      INTEGER nacv
!      PARAMETER (nacv=lcva) !bs # of anth. cond. vapors in CTM
!bs total # of cond. vapors in SOA_VBS
!     INTEGER, PARAMETER :: ncv=lsv_ocv   ! number of bins=8
!bs total # of cond. vapors in CTM
      REAL cblk(blksize,nspcsda) ! main array of variables
                                   ! particles [ug/m^3/s]
      REAL soilrat_in
                    ! emission rate of soil derived coars
                    ! input HNO3 to CBLK [ug/m^3]
      REAL nitrate_in
                    ! input NH3 to CBLK  [ug/m^3]
      REAL nh3_in
                    ! input SO4 vapor    [ug/m^3]
      REAL vsulf_in

      REAL so4rat_in
                    ! input SO4 formation [ug/m^3/sec]
                    ! pressure in cb
      REAL pres
                    ! temperature in K
      REAL temp
                    !bs
      REAL relhum
                    ! rel. humidity (0,1)   
      REAL ::  p(kts:kte),t(kts:kte),rh0(kts:kte)

      REAL ra(blksize )             ! aerodynamic resistance [ s m**-1 ]
      REAL ustar( blksize )         ! surface friction velocity [ m s**-1 ]
      REAL wstar( blksize )         ! convective velocity scale [ m s**-1 ]
      REAL pblh( blksize )          ! PBL height (m)
      REAL zntt( blksize )          ! Surface roughness length (m)
      REAL rmolm( blksize )         ! Inverse of Monin-Obukhov length (1/m)

      REAL blkprs(blksize)         ! pressure in cb
      REAL blkta(blksize)          ! temperature in K
      REAL blkdens(blksize)        ! Air density in kg/m3
!
! *** OUTPUT:
! *** atmospheric properties
      REAL xlm( blksize )           ! atmospheric mean free path [ m ]
      REAL amu( blksize )           ! atmospheric dynamic viscosity [ kg/m s ]
      
! *** followng is for future version       
      REAL vsed( blksize, naspcssed) ! sedimentation velocity [ m s**-1 ]
      REAL vdep( blksize, naspcsdep) ! deposition velocity [ m s**-1 ]

! *** modal diameters: [ m ]
      REAL dgnuc( blksize )         ! nuclei mode geometric mean diameter  [ m ]
      REAL dgacc( blksize )         ! accumulation geometric mean diameter [ m ]
      REAL dgcor( blksize )         ! coarse mode geometric mean diameter  [ m ]

! *** aerosol properties:
! *** Modal mass concentrations [ ug m**3 ]
      REAL pmassn( blksize )        ! mass concentration in aitken mode
      REAL pmassa( blksize )        ! mass concentration in accumulation mode
      REAL pmassc( blksize )        ! mass concentration in coarse mode

! *** average modal particle densities  [ kg/m**3 ]
      REAL pdensn( blksize )        ! average particle density in nuclei mode
      REAL pdensa( blksize )        ! average particle density in accumulation mode
      REAL pdensc( blksize )        ! average particle density in coarse mode

! *** average modal knudsen numbers
      REAL knnuc ( blksize )        ! nuclei mode  Knudsen number
      REAL knacc ( blksize )        ! accumulation Knudsen number
      REAL kncor ( blksize )        ! coarse mode  Knudsen number
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INTEGER :: i,j,k,l

INTEGER, SAVE :: icall

      do l=1,numaer
       do i=its,ite
        do j=jts,jte
           vgsa(i,j,l)=0.
        enddo
       enddo
      enddo
      vdep=0.

      do 100 j=jts,jte
         do 100 i=its,ite
!            cblk=epsilc
            do k=kts,kte
               t(k) = t_phy(i,k,j)
               p(k) = .001*p_phy(i,k,j)
               rh0(k) = rh(i,k,j)
            end do

            k=kts
               convfac = p(k)/rgasuniv/t(k)*1000.
               nitrate_in = chem(i,k,j,p_hno3)*convfac*mwhno3
               nh3_in =     chem(i,k,j,p_nh3)*convfac*mwnh3
               vsulf_in =   chem(i,k,j,p_sulf)*convfac*mwso4
               
!rs. nitrate, nh3, sulf
      blkprs(blksize)= 1.e3*p(k)                ! pressure in pa
      blkta(blksize) = t(k)         ! temperature in k
      ustar(blksize) = max(1.e-1,ust(i,j))
      wstar(blksize) = 0.
      pblh(blksize) = pbl(i,j)
      zntt(blksize) = znt(i,j)
      rmolm(blksize)= rmol(i,j)
      convfac2= 1./alt(i,k,j)    ! density of dry air
      blkdens(blksize)= convfac2
      cblk(1,vsulf) = max(epsilc,vsulf_in)
      cblk(1,vhno3) = max(epsilc,nitrate_in)
      cblk(1,vnh3)  = max(epsilc,nh3_in)

      if (p_naai >= param_first_scalar) &
         cblk(1,VNAAI ) =   max(epsilc,chem(i,k,j,p_naai)*convfac2)
      if (p_naaj >= param_first_scalar) &
         cblk(1,VNAAJ ) =   max(epsilc,chem(i,k,j,p_naaj)*convfac2)
      if (p_clai >= param_first_scalar) &
         cblk(1,VCLAI ) =   max(epsilc,chem(i,k,j,p_clai)*convfac2)
      if (p_claj >= param_first_scalar) &
         cblk(1,VCLAJ ) =   max(epsilc,chem(i,k,j,p_claj)*convfac2)

!liqy
!      if ( config_flags%do_isorropia ) then
!         cblk(1,VCAAI ) =   max(epsilc,caai(i,k,j)*convfac2)
!         cblk(1,VCAAJ ) =   max(epsilc,caaj(i,k,j)*convfac2)
!         cblk(1,VKAI )  =   max(epsilc,kai(i,k,j)*convfac2)
!         cblk(1,VKAJ )  =   max(epsilc,kaj(i,k,j)*convfac2)
!         cblk(1,VMGAI ) =   max(epsilc,mgai(i,k,j)*convfac2)
!         cblk(1,VMGAJ ) =   max(epsilc,mgaj(i,k,j)*convfac2)
!      endif

! Inorganic:
      cblk(1,vso4aj) =   max(epsilc,chem(i,k,j,p_so4aj)*convfac2)
      cblk(1,vso4ai) =   max(epsilc,chem(i,k,j,p_so4ai)*convfac2)
      cblk(1,vnh4aj) =   max(epsilc,chem(i,k,j,p_nh4aj)*convfac2)
      cblk(1,vnh4ai) =   max(epsilc,chem(i,k,j,p_nh4ai)*convfac2)
      cblk(1,vno3aj) =   max(epsilc,chem(i,k,j,p_no3aj)*convfac2)
      cblk(1,vno3ai) =   max(epsilc,chem(i,k,j,p_no3ai)*convfac2)

! POA:
      cblk(1,vpoa0j) =   max(epsilc,chem(i,k,j,p_poa0j)*convfac2)    ! ug/kg-air to ug/m3
      cblk(1,vpoa0i) =   max(epsilc,chem(i,k,j,p_poa0i)*convfac2)
      cblk(1,vpoa1j) =   max(epsilc,chem(i,k,j,p_poa1j)*convfac2)
      cblk(1,vpoa1i) =   max(epsilc,chem(i,k,j,p_poa1i)*convfac2)
      cblk(1,vpoa2j) =   max(epsilc,chem(i,k,j,p_poa2j)*convfac2)
      cblk(1,vpoa2i) =   max(epsilc,chem(i,k,j,p_poa2i)*convfac2)
      cblk(1,vpoa3j) =   max(epsilc,chem(i,k,j,p_poa3j)*convfac2)
      cblk(1,vpoa3i) =   max(epsilc,chem(i,k,j,p_poa3i)*convfac2)

! ASOA
      cblk(1,vasoa0j) =  max(epsilc,chem(i,k,j,p_asoa0j)*convfac2)  ! ug/kg-air to ug/m3
      cblk(1,vasoa0i) =  max(epsilc,chem(i,k,j,p_asoa0i)*convfac2)
      cblk(1,vasoa1j) =  max(epsilc,chem(i,k,j,p_asoa1j)*convfac2)
      cblk(1,vasoa1i) =  max(epsilc,chem(i,k,j,p_asoa1i)*convfac2)
      cblk(1,vasoa2j) =  max(epsilc,chem(i,k,j,p_asoa2j)*convfac2)
      cblk(1,vasoa2i) =  max(epsilc,chem(i,k,j,p_asoa2i)*convfac2)
      cblk(1,vasoa3j) =  max(epsilc,chem(i,k,j,p_asoa3j)*convfac2)
      cblk(1,vasoa3i) =  max(epsilc,chem(i,k,j,p_asoa3i)*convfac2)

! BSOA
      cblk(1,vbsoa1j) =  max(epsilc,chem(i,k,j,p_bsoa1j)*convfac2)
      cblk(1,vbsoa1i) =  max(epsilc,chem(i,k,j,p_bsoa1i)*convfac2)
      cblk(1,vbsoa2j) =  max(epsilc,chem(i,k,j,p_bsoa2j)*convfac2)
      cblk(1,vbsoa2i) =  max(epsilc,chem(i,k,j,p_bsoa2i)*convfac2)
      cblk(1,vbsoa3j) =  max(epsilc,chem(i,k,j,p_bsoa3j)*convfac2)
      cblk(1,vbsoa3i) =  max(epsilc,chem(i,k,j,p_bsoa3i)*convfac2)

! EC and P25
      cblk(1,vecj  )  =  max(epsilc,chem(i,k,j,p_ecj)*convfac2)
      cblk(1,veci  )  =  max(epsilc,chem(i,k,j,p_eci)*convfac2)
      cblk(1,vp25aj)  =  max(epsilc,chem(i,k,j,p_p25j)*convfac2)
      cblk(1,vp25ai)  =  max(epsilc,chem(i,k,j,p_p25i)*convfac2)

      cblk(1,VANTHA   ) =   max(epsilc,chem(i,k,j,p_antha)*convfac2)
      cblk(1,VSEAS    ) =   max(epsilc,chem(i,k,j,p_seas)*convfac2)
      cblk(1,VSOILA   ) =   max(epsilc,chem(i,k,j,p_soila)*convfac2)

      cblk(1,VPOLP    ) =   max(epsilc,chem(i,k,j,p_polp)*convfac2) !ts
      cblk(1,VPOLS    ) =   max(epsilc,chem(i,k,j,p_pols)*convfac2) !ts
      
      cblk(1,vnu0) =   max(epsilc,chem(i,k,j,p_nu0)*convfac2)
      cblk(1,vac0) =   max(epsilc,chem(i,k,j,p_ac0)*convfac2)

      cblk(1,vcorn    ) =   max(epsilc,chem(i,k,j,p_corn)*convfac2)
      cblk(1,vh2oaj) =   h2oaj(i,k,j)
      cblk(1,vh2oai) =   h2oai(i,k,j)
      cblk(1,vnu3  ) =   nu3(i,k,j)
      cblk(1,vac3  ) =   ac3(i,k,j)
      cblk(1,VCOR3    ) =   cor3(i,k,j)

! RAR: Debugging:
IF (config_flags%debug_soa .AND. icall<200) then
    WRITE(6,*) 'soa2_depdriver: icall,p_ac0,vac0, vac3 ',icall,p_so4aj,p_ac0, vac3
    WRITE(6,*) 'soa2_depdriver: convfac, convfac2 ',convfac, convfac2
    WRITE(6,*) 'soa2_depdriver: t(1), p(1), rh0(1), pblh ',t(1), p(1), rh0(1)
    WRITE(6,*) 'soa2_depdriver: cblk(1,vecj) ',cblk(1,vecj)
END IF

CALL modpar ( cblk,blkta,blkprs,      &
              pmassn, pmassa, pmassc, &
              pdensn, pdensa, pdensc, &
              xlm, amu,               &
              dgnuc, dgacc, dgcor,    &
              knnuc, knacc, kncor,    &
              config_flags%debug_soa,  &
              config_flags%do_isorropia )

! RAR: Only this option must be used!
        if (config_flags%aer_drydep_opt == 1) then
! for aerosol dry deposition, no CBLK in VDVG_2
             CALL VDVG_2( k,blkta, blkdens, aer_res(i,j), ustar, pblh, &
                          zntt, rmolm, amu, dgnuc, dgacc, dgcor, xlm,  &
                          knnuc, knacc, kncor,                         &
                          pdensn, pdensa, pdensc,                      &
                          vsed, vdep, config_flags%debug_soa )
        else
             CALL wrf_error_fatal ( 'RAR: Do not use another aer_drydep_opt!' )
        endif

        VGSA(i, j, vso4aj ) =  VDEP(1, vdmacc )
        VGSA(i, j, vso4ai ) =  VDEP(1, vdmnuc )
        VGSA(i, j, vnh4aj ) =  VGSA(i, j, vso4aj )
        VGSA(i, j, vnh4ai ) =  VGSA(i, j, vso4ai )
        VGSA(i, j, vno3aj ) =  VGSA(i, j, vso4aj )
        VGSA(i, j, vno3ai ) =  VGSA(i, j, vso4ai )

        if (p_naai >= param_first_scalar) VGSA(i, j, VNAAI )  =  VGSA(i, j, VSO4AI )
        if (p_naaj >= param_first_scalar) VGSA(i, j, VNAAJ )  =  VGSA(i, j, VSO4AJ )
        if (p_clai >= param_first_scalar) VGSA(i, j, VCLAI )  =  VGSA(i, j, VSO4AI )
        if (p_claj >= param_first_scalar) VGSA(i, j, VCLAJ )  =  VGSA(i, j, VSO4AJ )
       
!        if ( config_flags%do_isorropia ) then
!        VGSA(i, j, VCAAI )  =  VGSA(i,j,VSO4AI)
!        VGSA(i, j, VCAAJ )  = VGSA(i,j,VSO4AJ)
!        VGSA(i, j, VKAI )   =  VGSA(i, j,VSO4AI)
!        VGSA(i, j, VKAJ )   =  VGSA(i, j,VSO4AJ)
!        VGSA(i, j, VMGAI )  =  VGSA(i,j,VSO4AI)
!        VGSA(i, j, VMGAJ )  =  VGSA(i,j,VSO4AJ)
!        endif

        VGSA(i, j, vpoa0j ) =  VGSA(i, j, vso4aj )
        VGSA(i, j, vpoa0i ) =  VGSA(i, j, vso4ai )
        VGSA(i, j, vpoa1j ) =  VGSA(i, j, vso4aj )
        VGSA(i, j, vpoa1i ) =  VGSA(i, j, vso4ai )
        VGSA(i, j, vpoa2j ) =  VGSA(i, j, vso4aj )
        VGSA(i, j, vpoa2i ) =  VGSA(i, j, vso4ai )
        VGSA(i, j, vpoa3j ) =  VGSA(i, j, vso4aj )
        VGSA(i, j, vpoa3i ) =  VGSA(i, j, vso4ai )

        VGSA(i, j, vasoa0j ) =  VGSA(i, j, vso4aj )
        VGSA(i, j, vasoa0i ) =  VGSA(i, j, vso4ai )
        VGSA(i, j, vasoa1j ) =  VGSA(i, j, vso4aj )
        VGSA(i, j, vasoa1i ) =  VGSA(i, j, vso4ai )
        VGSA(i, j, vasoa2j ) =  VGSA(i, j, vso4aj )
        VGSA(i, j, vasoa2i ) =  VGSA(i, j, vso4ai )
        VGSA(i, j, vasoa3j ) =  VGSA(i, j, vso4aj )
        VGSA(i, j, vasoa3i ) =  VGSA(i, j, vso4ai )
!        VGSA(i, j, vasoa4J ) =  VGSA(i, j, vso4aj )
!        VGSA(i, j, vasoa4I ) =  VGSA(i, j, vso4ai )

        VGSA(i, j, vbsoa1j ) =  VGSA(i, j, vso4aj )
        VGSA(i, j, vbsoa1i ) =  VGSA(i, j, vso4ai )
        VGSA(i, j, vbsoa2j ) =  VGSA(i, j, vso4aj )
        VGSA(i, j, vbsoa2i ) =  VGSA(i, j, vso4ai )
        VGSA(i, j, vbsoa3j ) =  VGSA(i, j, vso4aj )
        VGSA(i, j, vbsoa3i ) =  VGSA(i, j, vso4ai )
!        VGSA(i, j, VBSOA4J ) =  VGSA(i, j, VSO4AJ )
!        VGSA(i, j, VBSOA4I ) =  VGSA(i, j, VSO4AI )
        !----------------------------------------------------------------------

!        VGSA(i, j, VORGPAJ )  =  VGSA(i, j, VSO4AJ )
!        VGSA(i, j, VORGPAI )  =  VGSA(i, j, VSO4AI )
        VGSA(i, j, vecj    )  =  VGSA(i, j, vso4aj )
        VGSA(i, j, veci    )  =  VGSA(i, j, vso4ai )
        VGSA(i, j, vp25aj  )  =  VGSA(i, j, vso4aj )
        VGSA(i, j, vp25ai  )  =  VGSA(i, j, vso4ai )


        VGSA(i, j, VANTHA  )  =  VDEP(1, VDMCOR )
        VGSA(i, j, VSEAS   )  =  VGSA(i, j, VANTHA )
        VGSA(i, j, VSOILA  )  =  VGSA(i, j, VANTHA )
        VGSA(i, j, vnu0    )  =  VDEP(1, vdnnuc )
        VGSA(i, j, vac0    )  =  VDEP(1, vdnacc )
        VGSA(i, j, VCORN   )  =  VDEP(1, VDNCOR )

        VGSA(i, j, VPOLP )  =  VGSA(i, j, VANTHA ) !ts similar to that of seas #DOUBT 
        VGSA(i, j, VPOLS )  =  VGSA(i, j, vso4aj  ) !ts similar to that od seas #DOUBT
!     enddo         ! k-loop

! RAR: Debugging:
IF (config_flags%debug_soa .AND. icall<200) then
    WRITE(6,*) 'soa2_depdriver: i,j,vpoa0j,vac0',i,j,vpoa0j,vac0
    WRITE(6,*) 'soa2_depdriver: vgsa(i,j,:),vdep ',vgsa(i,j,:),vdep
    WRITE(6,*) 'soa2_depdriver: vgsa(i,j,vecj) ',vgsa(i,j,vecj)
    icall= icall+1
END IF

100  continue      ! i,j-loop

END SUBROUTINE soa2_depdriver
! ///////////////////////////////////////////////////

    SUBROUTINE actcof(cat,an,gama,molnu,phimult)
! DESCRIPTION:
!  This subroutine computes the activity coefficients of (2NH4+,SO4--),
!  (NH4+,NO3-),(2H+,SO4--),(H+,NO3-),AND (H+,HSO4-) in aqueous
!  multicomponent solution, using Bromley's model and Pitzer's method.

! REFERENCES:
!   Bromley, L.A. (1973) Thermodynamic properties of strong electrolytes
!     in aqueous solutions.  AIChE J. 19, 313-320.

!   Chan, C.K. R.C. Flagen, & J.H.  Seinfeld (1992) Water Activities of
!     NH4NO3 / (NH4)2SO4 solutions, Atmos. Environ. (26A): 1661-1673.

!   Clegg, S.L. & P. Brimblecombe (1988) Equilibrium partial pressures
!     of strong acids over saline solutions - I HNO3,
!     Atmos. Environ. (22): 91-100

!   Clegg, S.L. & P. Brimblecombe (1990) Equilibrium partial pressures
!     and mean activity and osmotic coefficients of 0-100% nitric acid
!     as a function of temperature,   J. Phys. Chem (94): 5369 - 5380

!   Pilinis, C. and J.H. Seinfeld (1987) Continued development of a
!     general equilibrium model for inorganic multicomponent atmospheric
!     aerosols.  Atmos. Environ. 21(11), 2453-2466.

! ARGUMENT DESCRIPTION:
!     CAT(1) : conc. of H+    (moles/kg)
!     CAT(2) : conc. of NH4+  (moles/kg)
!     AN(1)  : conc. of SO4-- (moles/kg)
!     AN(2)  : conc. of NO3-  (moles/kg)
!     AN(3)  : conc. of HSO4- (moles/kg)
!     GAMA(2,1)    : mean molal ionic activity coeff for (2NH4+,SO4--)
!     GAMA(2,2)    :                                     (NH4+,NO3-)
!     GAMA(2,3)    :                                     (NH4+. HSO4-)
!     GAMA(1,1)    :                                     (2H+,SO4--)
!     GAMA(1,2)    :                                     (H+,NO3-)
!     GAMA(1,3)    :                                     (H+,HSO4-)
!     MOLNU   : the total number of moles of all ions.
!     PHIMULT : the multicomponent paractical osmotic coefficient.

! REVISION HISTORY:
!      Who       When        Detailed description of changes
!   ---------   --------  -------------------------------------------
!   S.Roselle   7/26/89   Copied parts of routine BROMLY, and began this
!                         new routine using a method described by Pilini
!                         and Seinfeld 1987, Atmos. Envirn. 21 pp2453-24
!   S.Roselle   7/30/97   Modified for use in Models-3
!   F.Binkowski 8/7/97    Modified coefficients BETA0, BETA1, CGAMA

!-----------------------------------------------------------------------
!...........INCLUDES and their descriptions
!      INCLUDE SUBST_XSTAT     ! M3EXIT status codes
!....................................................................

! Normal, successful completion           
      INTEGER, PARAMETER :: xstat0=0
! File I/O error                          
      INTEGER, PARAMETER :: xstat1=1
! Execution error                         
      INTEGER, PARAMETER :: xstat2=2
! Special  error                          
      INTEGER, PARAMETER :: xstat3=3
      CHARACTER*120 xmsg
! number of cations
      INTEGER, PARAMETER :: ncat=2
! number of anions
      INTEGER, PARAMETER :: nan=3

!...........ARGUMENTS and their descriptions
! tot # moles of all ions       
      REAL molnu
! multicomponent paractical osmo
      REAL phimult
      REAL cat(ncat) ! cation conc in moles/kg (input
      REAL an(nan) ! anion conc in moles/kg (input)
      REAL gama(ncat,nan) 
!...........SCRATCH LOCAL VARIABLES and their descriptions:
! mean molal ionic activity coef
      CHARACTER*16  pname ! driver program name
      SAVE pname

! anion indX                    
      INTEGER ian
      INTEGER icat
! cation indX                   

      REAL fgama
! ionic strength                
      REAL i
      REAL r
      REAL s
      REAL ta
      REAL tb
      REAL tc
      REAL texpv
      REAL trm
! 2*ionic strength              
      REAL twoi
! 2*sqrt of ionic strength      
      REAL twosri
      REAL zbar
      REAL zbar2
      REAL zot1
! square root of ionic strength 
      REAL sri
      REAL f2(ncat)
      REAL f1(nan)
      REAL zp(ncat) ! absolute value of charges of c
      REAL zm(nan) ! absolute value of charges of a
      REAL bgama(ncat,nan)
      REAL x(ncat,nan)
      REAL m(ncat,nan) ! molality of each electrolyte  
      REAL lgama0(ncat,nan) ! binary activity coefficients  
      REAL y(nan,ncat)
      REAL beta0(ncat,nan) ! binary activity coefficient pa
      REAL beta1(ncat,nan) ! binary activity coefficient pa
      REAL cgama(ncat,nan) ! binary activity coefficient pa
      REAL v1(ncat,nan) ! number of cations in electroly
      REAL v2(ncat,nan) 
! number of anions in electrolyt
      DATA zp/1.0, 1.0/
      DATA zm/2.0, 1.0, 1.0/
      DATA xmsg/' '/
      DATA pname/'ACTCOF'/

! *** Sources for the coefficients BETA0, BETA1, CGAMA:
! *** (1,1);(1,3)  - Clegg & Brimblecombe (1988)
! *** (2,3)        - Pilinis & Seinfeld (1987), cgama different
! *** (1,2)        - Clegg & Brimblecombe (1990)
! *** (2,1);(2,2)  - Chan, Flagen & Seinfeld (1992)

! *** now set the basic constants, BETA0, BETA1, CGAMA

  DATA beta0(1,1)/2.98E-2/, beta1(1,1)/0.0/, cgama(1,1)/4.38E-2    /        ! 2H+SO4
  DATA beta0(1,2)/1.2556E-1/, beta1(1,2)/2.8778E-1/, cgama(1,2)/ -5.59E-3 / ! HNO3
  DATA beta0(1,3)/2.0651E-1/, beta1(1,3)/5.556E-1/, cgama(1,3)/0.0       /  ! H+HSO4
  DATA beta0(2,1)/4.6465E-2/, beta1(2,1)/ -0.54196/,cgama(2,1)/ -1.2683E-3/ ! (NH4)2
  DATA beta0(2,2)/ -7.26224E-3/, beta1(2,2)/ -1.168858/,cgama(2,2)/3.51217E-5/ ! NH4NO3
  DATA beta0(2,3)/4.494E-2/, beta1(2,3)/2.3594E-1/, cgama(2,3)/ -2.962E-3 /
! NH4HSO
      DATA v1(1,1), v2(1,1)/2.0, 1.0/  ! 2H+SO4-
      DATA v1(2,1), v2(2,1)/2.0, 1.0/  ! (NH4)2SO4
      DATA v1(1,2), v2(1,2)/1.0, 1.0/  ! HNO3
      DATA v1(2,2), v2(2,2)/1.0, 1.0/  ! NH4NO3
      DATA v1(1,3), v2(1,3)/1.0, 1.0/  ! H+HSO4-
      DATA v1(2,3), v2(2,3)/1.0, 1.0/
!-----------------------------------------------------------------------
!  begin body of subroutine ACTCOF

!...compute ionic strength
! NH4HSO4                  
      i = 0.0
      DO icat = 1, ncat
        i = i + cat(icat)*zp(icat)*zp(icat)
      END DO

      DO ian = 1, nan
        i = i + an(ian)*zm(ian)*zm(ian)
      END DO

      i = 0.5*i
!...check for problems in the ionic strength
      IF (i==0.0) THEN
        DO ian = 1, nan
          DO icat = 1, ncat
            gama(icat,ian) = 0.0
          END DO
        END DO

!       xmsg = 'Ionic strength is zero...returning zero activities'
!       WRITE (6,*) xmsg
        RETURN

      ELSE IF (i<0.0) THEN
!        xmsg = 'Ionic strength below zero...negative concentrations'
!        CALL wrf_error_fatal ( xmsg )

! RAR: included from Jan Kazil
        xmsg = 'WARNING: Ionic strength below zero (= negative ion concentrations) - setting ion concentrations to zero.'
        call wrf_message(xmsg)
        DO ian = 1, nan
          DO icat = 1, ncat
            gama(icat,ian) = 0.0
          END DO
        END DO
        RETURN

      END IF

!...compute some essential expressions
      sri = sqrt(i)
      twosri = 2.0*sri
      twoi = 2.0*i
      texpv = 1.0 - exp(-twosri)*(1.0+twosri-twoi)
      r = 1.0 + 0.75*i
      s = 1.0 + 1.5*i
      zot1 = 0.511*sri/(1.0+sri)

!...Compute binary activity coeffs
      fgama = -0.392*((sri/(1.0+1.2*sri)+(2.0/1.2)*alog(1.0+1.2*sri)))
      DO icat = 1, ncat
        DO ian = 1, nan

          bgama(icat,ian) = 2.0*beta0(icat,ian) + ( 2.0*beta1(icat,ian)/(4.0*i) )*texpv

!...compute the molality of each electrolyte for given ionic strength

          m(icat,ian) = (cat(icat)**v1(icat,ian)*an(ian)**v2(icat,ian))** &
            (1.0/(v1(icat,ian)+v2(icat,ian)))

!...calculate the binary activity coefficients

          lgama0(icat,ian) = (zp(icat)*zm(ian)*fgama+m(icat,ian)*(2.0*v1(icat, &
            ian)*v2(icat,ian)/(v1(icat,ian)+v2(icat,ian))*bgama(icat, &
            ian))+m(icat,ian)*m(icat,ian)*(2.0*(v1(icat,ian)* &
            v2(icat,ian))**1.5/(v1(icat,ian)+v2(icat,ian))*cgama(icat, &
            ian)))/2.302585093

        END DO
      END DO

!...prepare variables for computing the multicomponent activity coeffs

      DO ian = 1, nan
        DO icat = 1, ncat
          zbar = (zp(icat)+zm(ian))*0.5
          zbar2 = zbar*zbar
          y(ian,icat) = zbar2*an(ian)/i
          x(icat,ian) = zbar2*cat(icat)/i
        END DO
      END DO

      DO ian = 1, nan
        f1(ian) = 0.0
        DO icat = 1, ncat
          f1(ian) = f1(ian) + x(icat,ian)*lgama0(icat,ian) + &
            zot1*zp(icat)*zm(ian)*x(icat,ian)
        END DO
      END DO

      DO icat = 1, ncat
        f2(icat) = 0.0
        DO ian = 1, nan
          f2(icat) = f2(icat) + y(ian,icat)*lgama0(icat,ian) + &
            zot1*zp(icat)*zm(ian)*y(ian,icat)
        END DO
      END DO

!...now calculate the multicomponent activity coefficients

      DO ian = 1, nan
        DO icat = 1, ncat

          ta = -zot1*zp(icat)*zm(ian)
          tb = zp(icat)*zm(ian)/(zp(icat)+zm(ian))
          tc = (f2(icat)/zp(icat)+f1(ian)/zm(ian))
          trm = ta + tb*tc

          IF (trm>30.0) THEN
            gama(icat,ian) = 1.0E+30
!           xmsg = 'Multicomponent activity coefficient is extremely large'
!           WRITE (6,*) xmsg
          ELSE
            gama(icat,ian) = 10.0**trm
          END IF

        END DO
      END DO

      RETURN
!ia*********************************************************************
    END SUBROUTINE actcof

!ia
!ia     AEROSOL DYNAMICS DRIVER ROUTINE					*
!ia     based on MODELS3 formulation by FZB
!ia     Modified by IA in November 97
!ia
!ia     Revision history
!ia     When    WHO     WHAT
!ia     ----    ----    ----
!ia     ????    FZB     BEGIN
!ia     05/97   IA      Adapted for use in CTM2-S
!ia     11/97   IA      Modified for new model version
!ia                     see comments under iarev02
!ia
!ia     Called BY:      RPMMOD3
!ia
!ia     Calls to:       EQL3, MODPAR, COAGRATE, NUCLCOND, AEROSTEP
!ia                     GETVSED
!ia
!ia*********************************************************************

SUBROUTINE aeroproc(layer,cblk,dt,blkta,blkprs,blkdens,blkrh,so4rat,          &
                orgpoa1rat,orgpoa2rat,orgpoa3rat,  &
                organt0rat,organt1rat,organt2rat,organt3rat,  &
                orgbio1rat,orgbio2rat,orgbio3rat,  &

            !    orgbio1rat,orgbio2rat,orgbio3rat,orgbio4rat,drog,ldrog_vbs,ncv,nacv,epm25i, &
            !    epm25j,eorgi,eorgj,eeci,eecj,epmcoarse,esoil,eseas,xlm,amu,dgnuc, &

                drog_in,epm25i,epm25j,eorgi,eorgj,eeci,eecj,epmcoarse,esoil,eseas,epolp,epols,xlm,amu,    &
                dgnuc,dgacc,dgcor,pmassn,pmassa,pmassc,pdensn,pdensa,pdensc,                  &
                knnuc,knacc,kncor,fconcn,fconca,fconcn_org,fconca_org,dmdt,dndt,cgrn3,cgra3,  &
                urn00,ura00,brna01,c30,deltaso4a,igrid,jgrid,kgrid,brrto, diagn_opt,          &
                do_isorropia,do_n2o5_hetchem)

!USE module_configure, only: grid_config_rec_type
!TYPE (grid_config_rec_type), INTENT (in) :: config_flags

!     IMPLICIT NONE
! dimension of arrays             
!
! number of species in CBLK       
!      INTEGER nspcsda
! actual number of cells in arrays
!      INTEGER numcells

LOGICAL :: diagn_opt, do_isorropia
INTEGER :: do_n2o5_hetchem

! number of k-level               
      INTEGER layer
! of organic aerosol precursor  
      REAL cblk(blksize,nspcsda) ! main array of variables (INPUT a

      REAL dt
! *** Meteorological information:

! synchronization time  [s]       
      REAL blkta(blksize) ! Air temperature [ K ]                  
      REAL blkprs(blksize) ! Air pressure in [ Pa ]                 
      REAL blkdens(blksize) ! Air density  [ kg/ m**3 ]              
      REAL blkrh(blksize) 
! *** Chemical production rates: [ ug / m**3 s ]

! Fractional relative humidity           
      REAL so4rat(blksize) 
! sulfate gas-phase production rate
! total # of cond. vapors & SOA species 
!      INTEGER ncv
!      INTEGER nacv
!bs * organic condensable vapor production rate
! # of anthrop. cond. vapors & SOA speci
      REAL drog_in(ldrog_soa2) !bs

! POA growth rates
      REAL orgpoa1rat(blksize), orgpoa2rat(blksize), orgpoa3rat(blksize)

! *** anthropogenic organic aerosol mass production rates from aromatics
! Delta ROG conc. [ppm]              
      REAL organt0rat(blksize), organt1rat(blksize), organt2rat(blksize), organt3rat(blksize)

! *** biogenic organic aerosol production rates
      REAL orgbio1rat(blksize), orgbio2rat(blksize), orgbio3rat(blksize)

! *** Primary emissions rates: [ ug / m**3 s ]
! *** emissions rates for unidentified PM2.5 mass
      REAL epm25i, epm25j
! *** emissions rates for primary organic aerosol
      REAL eorgi, eorgj
! *** emissions rates for elemental carbon
      REAL eeci, eecj

! *** emissions rates for coarse mode particles
! Accumululaton mode                    
      REAL esoil(blksize) ! soil derived coarse aerosols
      REAL eseas(blksize) ! marine coarse aerosols
      REAL epmcoarse(blksize)
! ts added for pollens biological particles
      REAL epolp(blksize) !ts biological coarse aerosols   
      REAL epols(blksize) !ts biological accum# aerosols  
! *** OUTPUT:
! *** atmospheric properties
! anthropogenic coarse aerosols
      REAL xlm(blksize) ! atmospheric mean free path [ m ]  
      REAL amu(blksize) 
! *** modal diameters: [ m ]

! atmospheric dynamic viscosity [ kg
      REAL dgnuc(blksize) ! nuclei mode geometric mean diamete
      REAL dgacc(blksize) ! accumulation geometric mean diamet
      REAL dgcor(blksize) 

! *** aerosol properties:
! *** Modal mass concentrations [ ug m**3 ]
! coarse mode geometric mean diamete
      REAL pmassn(blksize) ! mass concentration in Aitken mode 
      REAL pmassa(blksize) ! mass concentration in accumulation
      REAL pmassc(blksize) 
! *** average modal particle densities  [ kg/m**3 ]

! mass concentration in coarse mode 
      REAL pdensn(blksize) ! average particle density in nuclei
      REAL pdensa(blksize) ! average particle density in accumu
      REAL pdensc(blksize) 
! *** average modal Knudsen numbers

! average particle density in coarse
      REAL knnuc(blksize) ! nuclei mode  Knudsen number       
      REAL knacc(blksize) ! accumulation Knudsen number       
      REAL kncor(blksize) 
! ***  modal condensation factors ( see comments in NUCLCOND )

! coarse mode  Knudsen number       
      REAL fconcn(blksize)
      REAL fconca(blksize)
!bs
      REAL fconcn_org(blksize)
      REAL fconca_org(blksize)

! *** Rates for secondary particle formation:
! *** production of new mass concentration [ ug/m**3 s ]
      REAL dmdt(blksize) !                                 by particle formation

! *** production of new number concentration [ number/m**3 s ]
! rate of production of new mass concen
      REAL dndt(blksize) !                                 by particle formation

! *** growth rate for third moment by condensation of precursor
!      vapor on existing particles [ 3rd mom/m**3 s ]

! rate of producton of new particle num
      REAL cgrn3(blksize) !  Aitken mode                          
      REAL cgra3(blksize) 
! *** Rates for coaglulation: [ m**3/s ]

! *** Unimodal Rates:

!  Accumulation mode                    
      REAL urn00(blksize) ! Aitken mode 0th moment self-coagulation ra
      REAL ura00(blksize) 

! *** Bimodal Rates:  Aitken mode with accumulation mode ( d( Aitken mod

! accumulation mode 0th moment self-coagulat
      REAL brna01(blksize) 
! *** 3rd moment intermodal transfer rate replaces coagulation rate ( FS
! rate for 0th moment                     
      REAL c30(blksize)                                                        ! by intermodal c
      REAL brrto

! *** other processes

! intermodal 3rd moment transfer r
      REAL deltaso4a(blksize) !                                  sulfate aerosol by condensation   [ u

!      INTEGER NN, VV ! loop indICES
! increment of concentration added to

! ////////////////////// Begin code ///////////////////////////////////
! concentration lower limit
      CHARACTER*16 pname
      PARAMETER (pname=' AEROPROC       ')

      INTEGER, PARAMETER :: unit=20, isorop=0
      INTEGER igrid,jgrid,kgrid

!        if(blkta(1).ge.233.15.and.blkrh(1).ge.0.1 .and. do_isorropia)then
!            CALL eql3(blksize,nspcsda,numcells,cblk,blkta,blkrh)
!        else if (blkta(1).ge.233.15.and.blkrh(1).ge.0.1 .and. (.not. do_isorropia) )then
!            CALL eql4(cblk,blkta,blkrh,diagn_opt)
!        endif
      CALL eql4(cblk,blkta,blkrh,diagn_opt)
      if ( do_n2o5_hetchem .eq. 1 ) then
      CALL n2o5het(blksize,nspcsda,numcells,dt,cblk,blkta,blkrh,dgnuc,dgacc,dgcor,igrid,jgrid,kgrid)
      endif
! *** get size distribution information:

      CALL modpar(cblk,blkta,blkprs,pmassn,pmassa, &
             pmassc,pdensn,pdensa,pdensc,xlm,amu,dgnuc,dgacc,dgcor, &
             knnuc,knacc,kncor, diagn_opt, do_isorropia )

! *** Calculate coagulation rates for fine particles:

      CALL coagrate(cblk,blkta,pdensn,pdensa,amu, &
              dgnuc,dgacc,knnuc,knacc,urn00,ura00,brna01,c30)

! *** get condensation and particle formation (nucleation) rates:
      CALL nuclcond(cblk,dt,layer,blkta,blkprs,blkrh,so4rat, &
             orgpoa1rat,orgpoa2rat,orgpoa3rat, &
             organt0rat,organt1rat,organt2rat,organt3rat, &
             orgbio1rat,orgbio2rat,orgbio3rat, &
             drog_in,dgnuc,dgacc, &
             fconcn,fconca,fconcn_org,fconca_org,dmdt,dndt,deltaso4a, &
             cgrn3,cgra3,brrto,diagn_opt, &
             igrid,jgrid,kgrid)

! *** advance forward in time DT seconds:
      CALL aerostep(layer,cblk,dt,so4rat, &
             orgpoa1rat,orgpoa2rat,orgpoa3rat, &
             organt0rat,organt1rat,organt2rat,organt3rat, &
             orgbio1rat,orgbio2rat,orgbio3rat, &
             epm25i,epm25j,eorgi,eorgj,eeci,eecj,epmcoarse,esoil,eseas,epolp,epols, &
             dgnuc,dgacc,fconcn,fconca,fconcn_org,fconca_org,pmassn,pmassa,pmassc, &
             dmdt,dndt,deltaso4a,urn00,ura00,brna01,c30,cgrn3,cgra3,igrid,jgrid,kgrid,diagn_opt,do_isorropia)

! *** get new distribution information:
      CALL modpar(cblk,blkta,blkprs,pmassn,pmassa, &
             pmassc,pdensn,pdensa,pdensc,xlm,amu,dgnuc,dgacc,dgcor,knnuc,knacc,kncor,diagn_opt, do_isorropia)

      RETURN
    END SUBROUTINE aeroproc
!//////////////////////////////////////////////////////////////////
   
   SUBROUTINE n2o5het(blksize,nspcsda,numcells,dt,cblk,blkta,blkrh,dgnuc,dgacc,dgcor,igrid,jgrid,kgrid)

! dimension of arrays
      INTEGER blksize
! actual number of cells in arrays
      INTEGER numcells
! nmber of species in CBLK        
      INTEGER nspcsda
      REAL cblk(blksize,nspcsda)
      REAL dt
! *** Meteorological information in blocked arays:
      REAL blkta(blksize) ! Air temperature [ K ]                   
      REAL blkrh(blksize) ! Fractional relative humidity            
      REAL dgnuc(blksize) ! nuclei mode geometric mean diamete
      REAL dgacc(blksize) ! accumulation geometric mean diamet
      REAL dgcor(blksize)
!        
      Integer igrid,jgrid,kgrid

      INTEGER lcell ! loop counter                                   
! air temperature                             
      REAL temp
!relative humidity.
      REAL rh
!aerosol number density
      REAL nnu
      REAL nac
!aerosol mean diameter
      REAL dnu!nuclei
      REAL dac !accumulation
!aerosol surface area density
      REAL snu
      REAL sac
!uptake of n2o5 on aerosols
      REAL gamn2o5
!n2o5 molecular speed 
      REAL cn2o5
!reaction rate constants of N2O5 hydrolysis 
      REAL kn2o5
!yield of clno2   
      REAL yclno2


      REAL ah2o
      REAL acl
      REAL ano3
      REAL gn2o5
      
      REAL mwh2o
      PARAMETER (mwh2o = 18.015)
      REAL mwcl
      PARAMETER (mwcl = 35.453)
      REAL mwno3
      PARAMETER (mwno3 = 62.004)
      REAL mwn2o5
      PARAMETER (mwn2o5 = 108.009)
      REAL mwclno2
      PARAMETER (mwclno2 = 81.458)
      REAL deln2o5
      REAL pclno2
      REAL pno3

      REAL fraci,fracj,fracij
      REAL rgasuniv
      PARAMETER (rgasuniv = 8.314510)
      REAL pirs
      PARAMETER (pirs = 3.14)
      INTEGER xxx
      PARAMETER (xxx = 1)

      real vaer
!==================================================          

        DO lcell = 1, numcells


        temp = blkta(lcell)
        rh = blkrh(lcell)
        nnu = cblk(lcell,vnu0)          !#/m3-dry air
        nac = cblk(lcell,vac0)
        dnu = dgnuc(lcell)              !m
        dac = dgacc(lcell)
        vaer = (pirs/6.0) * (cblk(lcell,vnu3) + cblk(lcell,vac3))
!aerosol volume in i and j mode.
!=================================================      
!convert the unit from ug/m3 to mol/L (in aerosol solution)
        ah2o = ( cblk(lcell,vh2oaj) + cblk(lcell,vh2oai) ) * 1.0E-9 / (mwh2o*vaer)

!convert the unit from ug/m3 to mol/L (in aerosol solution)
        acl  = ( cblk(lcell,vclaj) + cblk(lcell,vclai) ) * 1.0E-9/(mwcl*vaer)
        ano3 = ( cblk(lcell,vno3aj) + cblk(lcell,vno3ai) ) * 1.0E-9/(mwno3*vaer)

! convert the unit from ug/m3 to mol/L in air atmosphere.
        gn2o5 = cblk(lcell,vn2o5) * 1.0E-9 /mwn2o5

        cblk(lcell,vgamn2o5) = 3.2E-8 * ( 1.15E6 - 1.15E6 * exp(-1.3E-1* ah2o )) * &
                              ( 1 - (1/((6E-2*ah2o/ano3)+1+(29*acl/ano3))))

        cblk(lcell,vsnu) = nnu*dnu*dnu*esn16*pirs
        cblk(lcell,vsac) = nac*dac*dac*esa16*pirs

        cblk(lcell,vcn2o5) = SQRT( 8.0 * rgasuniv * temp * 1000 / ( pirs* mwn2o5) )
        cblk(lcell,vkn2o5) = cblk(lcell,vcn2o5) * ( cblk(lcell,vsnu) & 
                             +cblk(lcell,vsac) ) * cblk(lcell,vgamn2o5) / 4
        deln2o5 = gn2o5-gn2o5*exp(-1*cblk(lcell,vkn2o5)*dt)      !mole/L in atmosphere

        cblk(lcell,vyclno2)= 1/(1+ah2o/(483*acl))

        pclno2=deln2o5*cblk(lcell,vyclno2)   !mol/L in atmosphere

        if (acl*vaer .lt. pclno2) then
                pclno2=abs(acl*vaer-epsilc*epsilc)
                cblk(lcell,vyclno2)=pclno2/deln2o5
        end if


        pno3 = deln2o5 * ( 2 - cblk(lcell,vyclno2) ) !mole/L in atmosphere

        cblk(lcell,vclno2) = cblk(lcell,vclno2) + pclno2*mwclno2*1.0E9

        fraci = cblk(lcell,vso4ai)/(cblk(lcell,vso4aj)+cblk(lcell,vso4ai))
        fracj = 1.0 - fraci

        cblk(lcell,vclaj)=max(epsilc*epsilc,(acl*vaer-pclno2))*mwcl*1.0E9*fracj
        cblk(lcell,vclai)=max(epsilc*epsilc,(acl*vaer-pclno2))*mwcl*1.0E9*fraci

        cblk(lcell,vn2o5) = cblk(lcell,vn2o5)*exp(-1*cblk(lcell,vkn2o5)*dt)
        cblk(lcell,vno3ai) = (ano3*vaer + pno3) * mwno3 * 1.0E9* fraci
        cblk(lcell,vno3aj) = (ano3*vaer + pno3) * mwno3 * 1.0E9* fracj

        END DO

        END SUBROUTINE n2o5het


! *** Time stepping code advances the aerosol moments one timestep;
    SUBROUTINE aerostep(layer,cblk,dt,so4rat, &
                 orgpoa1rat,orgpoa2rat,orgpoa3rat, &
                 organt0rat,organt1rat,organt2rat,organt3rat, &
                 orgbio1rat,orgbio2rat,orgbio3rat, &
                 epm25i,epm25j,eorgi,eorgj,eeci,eecj,epmcoarse,esoil,eseas,epolp,epols, &
                 dgnuc,dgacc,fconcn,fconca,fconcn_org,fconca_org,pmassn,pmassa,pmassc, &
                 dmdt,dndt,deltaso4a,urn00,ura00,brna01,c30,cgrn3,cgra3,igrid,jgrid,kgrid, diagn_opt, do_isorropia)

! ***  DESCRIPTION: Integrate the Number and Mass equations
!                   for each mode over the time interval DT.
!      PRECONDITIONS:
!       AEROSTEP() must follow calls to all other dynamics routines.

! ***   Revision history:
!       Adapted 3/95 by UAS and CJC from EAM2's code.
!       Revised 7/29/96 by FSB to use block structure
!       Revised 11/15/96 by FSB dropped flow-through and cast
!                           number solver into Riccati equation form.
!       Revised 8/8/97 by FSB to have mass in Aitken and accumulation mode
!                        each predicted rather than total mass and
!                        Aitken mode mass. Also used a local approximation
!                        the error function. Also added coarse mode.
!       Revised 9/18/97 by FSB to fix mass transfer from Aitken to
!                       accumulation mode by coagulation
!       Revised 10/27/97 by FSB to modify code to use primay emissions
!                        and to correct 3rd moment updates.
!                        Also added coarse mode.
!       Revised 11/4/97 by FSB to fix error in other anthropogenic PM2.5
!       Revised  11/5/97 by FSB to fix error in MSTRNSFR
!       Revised  11/6/97 FSB to correct the expression for FACTRANS to
!                        remove the 6/pi coefficient. UAS found this.
!       Revised 12/15/97 by FSB to change equations for mass concentratin
!                        to a chemical production form with analytic
!                        solutions for the Aitken mode and to remove
!                        time stepping of the 3rd moments. The mass concentration
!                        in the accumulation mode is updated with a forward
!                        Eulerian step.
!       Revised 1/6/98   by FSB Lowered minimum concentration for
!                        sulfate aerosol to 0.1 [ ng / m**3 ].
!       Revised 1/12/98  C30 replaces BRNA31 as a variable. C30 represents
!                        intermodal transfer rate of 3rd moment in place
!                        of 3rd moment coagulation rate.
!       Revised 5/5/98   added new renaming criterion based on diameters
!       Added   3/23/98  by BS condensational groth factors for organics

!**********************************************************************
!     IMPLICIT NONE

! *** ARGUMENTS:

! dimension of arrays             
!      INTEGER blksize
! actual number of cells in arrays
!      INTEGER numcells
! nmber of species in CBLK        
!      INTEGER nspcsda

LOGICAL :: diagn_opt,do_isorropia

! model layer                     
      INTEGER layer
      REAL cblk(blksize,nspcsda) ! main array of variables          
      REAL dt
! *** Chemical production rates: [ ug / m**3 s ]

! time step [sec]                  
      REAL so4rat(blksize)  ! sulfate gas-phase production rate

      REAL orgpoa1rat(blksize), orgpoa2rat(blksize), orgpoa3rat(blksize)

! anthropogenic organic aerosol mass production rates
      REAL organt0rat(blksize), organt1rat(blksize), organt2rat(blksize), organt3rat(blksize)

! biogenic organic aerosol production rates
      REAL orgbio1rat(blksize), orgbio2rat(blksize), orgbio3rat(blksize)

! *** Primary emissions rates: [ ug / m**3 s ]
! *** emissions rates for unidentified PM2.5 mass
      REAL epm25i, epm25j
! *** emissions rates for primary organic aerosol
      REAL eorgi, eorgj
! *** emissions rates for elemental carbon
      REAL eeci, eecj

! *** emissions rates for coarse mode particles
! Accumululaton mode                    
      REAL esoil(blksize) ! soil derived coarse aerosols
      REAL eseas(blksize) ! marine coarse aerosols
      REAL epmcoarse(blksize)
! ts added for pollens biological particles
      REAL epolp(blksize) !ts biological coarse aerosols   
      REAL epols(blksize) !ts biological accum# aerosols  
! anthropogenic coarse aerosols         
      REAL dgnuc(blksize) ! nuclei mode mean diameter [ m ]
      REAL dgacc(blksize) 
! accumulation                          
      REAL fconcn(blksize)                                 ! Aitken mode  [ 1 / s ]
! reciprocal condensation rate          
      REAL fconca(blksize)                                 ! acclumulation mode [ 1 / s ]
! reciprocal condensation rate          
      REAL fconcn_org(blksize)                                 ! Aitken mode  [ 1 / s ]
! reciprocal condensation rate for organ
      REAL fconca_org(blksize)                                 ! acclumulation mode [ 1 / s ]
! reciprocal condensation rate for organ
      REAL dmdt(blksize)                                 ! by particle formation [ ug/m**3 /s ]
! rate of production of new mass concent
      REAL dndt(blksize)                                 ! by particle formation [ number/m**3 /s
! rate of producton of new particle numb
      REAL deltaso4a(blksize)                                 ! sulfate aerosol by condensation [ ug/m
! increment of concentration added to   
      REAL urn00(blksize) ! Aitken intramodal coagulation rate    
      REAL ura00(blksize) ! Accumulation mode intramodal coagulati
      REAL brna01(blksize) ! bimodal coagulation rate for number   
      REAL c30(blksize)       							! by intermodal coagulation
! intermodal 3rd moment transfer rate by
      REAL cgrn3(blksize) ! growth rate for 3rd moment for Aitken 
      REAL cgra3(blksize) 
! *** Modal mass concentrations [ ug m**3 ]

! growth rate for 3rd moment for Accumul
      REAL pmassn(blksize) ! mass concentration in Aitken mode 
      REAL pmassa(blksize) ! mass concentration in accumulation
      REAL pmassc(blksize) 

! *** Local Variables
! mass concentration in coarse mode 
      INTEGER l, lcell, spc
! ** following scratch variables are used for solvers

! *** variables needed for modal dynamics solvers:
! Loop indices                   
      REAL*8 a, b, c
      REAL*8 m1, m2, y0, y
      REAL*8 dhat, p, pexpdt, expdt
      REAL*8 loss, prod, pol, lossinv
! mass intermodal transfer by coagulation           
      REAL mstrnsfr

      REAL factrans

! *** CODE additions for renaming
      REAL getaf2
      REAL aaa, xnum, xm3, fnum, fm3, phnum, phm3 ! Defined below
      REAL erf, erfc ! Error and complementary error function
      REAL xx
! dummy argument for ERF and ERFC          
! a numerical value for a minimum concentration
      INTEGER igrid,jgrid,kgrid

      INTEGER, SAVE :: icall
! *** This value is smaller than any reported tropospheric concentration

! *** Statement function given for error function. Source is
!     Meng, Z., and J.H.Seinfeld (1994) On the source of the submicromet
!      droplet mode of urban and regional aerosols. Aerosol Sci. and Tec
!      20:253-265. They cite Reasearch & Education Asociation (REA), (19
!      Handbook of Mathematical, Scientific, and Engineering Formulas,
!      Tables, Functions, Graphs, Transforms: REA, Piscataway, NJ. p. 49

      erf(xx) = sqrt(1.0-exp(-4.0*xx*xx/pirs))
      erfc(xx) = 1.0 - erf(xx)
!     ::::::::::::::::::::::::::::::::::::::::

! ///// begin code
! *** set up time-step integration

      DO l = 1, numcells

! *** code to move number forward by one time step.
! *** solves the Ricatti equation:
!     dY/dt = C - A * Y ** 2 - B * Y
!     Coded 11/21/96 by Dr. Francis S. Binkowski

! *** Aitken mode:
! *** coefficients
        a = urn00(l)
        b = brna01(l)*cblk(l,vac0)
        c = dndt(l) + factnumn*( anthfac*(epm25i + eeci) + orgfac*eorgi )

! includes primary emissions 
        y0 = cblk(l,vnu0) 
! ***  trap on C = 0

! initial condition                           
        IF (c>0.0D0) THEN

          dhat = sqrt(b*b+4.0D0*a*c)

          m1 = 2.0D0*a*c/(b+dhat)

          m2 = -0.5D0*(b+dhat)

          p = -(m1-a*y0)/(m2-a*y0)

          pexpdt = p*exp(-dhat*dt)

          y = (m1+m2*pexpdt)/(a*(1.0D0+pexpdt)) 
! solution                       
        ELSE

! *** rearrange solution for NUMERICAL stability
!     note If B << A * Y0, the following form, although
!     seemingly awkward gives the correct answer.

          expdt = exp(-b*dt)
          IF (expdt<1.0D0) THEN
            y = b*y0*expdt/(b+a*y0*(1.0D0-expdt))
          ELSE
            y = y0
          END IF

        END IF

        cblk(l,vnu0) = max(nummin_i,y) 

! *** now do accumulation mode number

! *** coefficients

! update                     
        a = ura00(l)
        b = 0.0D0        ! NOTE B = 0.0
        c = factnuma*( anthfac*(epm25j+eecj) + orgfac*eorgj + polsfac*epols(l))
! includes primary emissi
        y0 = cblk(l,vac0) 
! *** this equation requires special handling, because C can be zero.
!     if this happens, the form of the equation is different:

! initial condition                           
!       print *,vac0,y0,c,nummin_j,a
        IF (c>0.0D0) THEN

          dhat = sqrt(4.0D0*a*c)

          m1 = 2.0D0*a*c/dhat

          m2 = -0.5D0*dhat

          p = -(m1-a*y0)/(m2-a*y0)

!       print *,p,-dhat,dt,-dhat*dt
!       print *,exp(-dhat*dt)
          pexpdt = p*exp(-dhat*dt)

          y = (m1+m2*pexpdt)/(a*(1.0D0+pexpdt)) 
! solution                       
        ELSE

          y = y0/(1.0D0+dt*a*y0) 
!       print *,dhat,y0,dt,a
          y = y0/(1.+dt*a*y0) 
!       print *,y
! correct solution to equation
        END IF

        cblk(l,vac0) = max(nummin_j,y) 
! *** now do coarse mode number neglecting coagulation
! update                     
!       print *,soilfac,seasfac,anthfac,esoil(l),eseas(l),epmcoarse(l)
        prod = soilfac*esoil(l) + seasfac*eseas(l) + anthfac*epmcoarse(l) + polpfac*epolp(l)

!       print *,cblk(l,vcorn),factnumc,prod
        cblk(l,vcorn) = cblk(l,vcorn) + factnumc*prod*dt

! *** Prepare to advance modal mass concentration one time step.

! *** Set up production and and intermodal transfer terms terms:
!       print *,cgrn3(l),epm25i(l),eeci(l),orgfac,eorgi(l)
      cgrn3(l) = cgrn3(l) + anthfac*(epm25i + eeci) + orgfac*eorgi

! includes growth from pri
      cgra3(l) = cgra3(l) + c30(l) + anthfac*(epm25j + eecj) + orgfac*eorgj ! + &  ! and transfer of 3rd moment
                            !polsfac*epols(l)
                                             ! intermodal coagulation

! *** set up transfer coefficients for coagulation between Aitken and ac

! *** set up special factors for mass transfer from the Aitken to accumulation
!     intermodal coagulation. The mass transfer rate is proportional to
!     transfer rate, C30. The proportionality factor is p/6 times the the
!     density. The average particle density for a species is the species
!     divided by the particle volume concentration, pi/6 times the 3rd m
!     The p/6 coefficients cancel.

! includes growth from prim
!       print *,'loss',vnu3,c30(l),cblk(l,vnu3)
        loss = c30(l)/cblk(l,vnu3) 

! Normalized coagulation transfer r
        factrans = loss*dt                            ! yields an estimate of the amount of mass t
     ! the Aitken to the accumulation mode in the

! Multiplying this factor by the species con
!       print *,'factrans = ',factrans,loss
        expdt = exp(-factrans)                               ! decay term is common to all Aitken mode
!       print *,'factrans = ',factrans,loss,expdt
! variable name is re-used here. This expo
        lossinv = 1.0/loss
! *** now advance mass concentrations one time step.

! ***  update sulfuric acid vapor concentration by removing mass concent
!      condensed sulfate and newly produced particles.
! *** The method follows Youngblood and Kreidenweis, Further Development
!     of a Bimodal Aerosol Dynamics Model, Colorado State University Dep
!     Atmospheric Science Paper Number 550, April,1994, pp 85-89.
! set up for multiplication rather than divi
        cblk(l,vsulf) = max(conmin, cblk(l,vsulf) - (deltaso4a(l)+dmdt(l)*dt) )

! *** Solve Aitken-mode equations of form: dc/dt = P - L*c
! *** Solution is:     c(t0 + dt) = p/L + ( c(0) - P/L ) * exp(-L*dt)
! RAR: Debugging:
if (diagn_opt .AND. icall<200 .AND. kgrid<3) then
    WRITE(6,*) 'aerostep: icall, c ', icall, c
    WRITE(6,*) 'aerostep: epm25i,epm25j,eorgi,eorgj,eeci,eecj ',epm25i,epm25j,eorgi,eorgj,eeci,eecj
    WRITE(6,*) 'aerostep: deltaso4a,factrans,lossinv',deltaso4a,factrans,lossinv
    WRITE(6,*) 'aerostep: fconcn,fconca',fconcn,fconca
    WRITE(6,*) 'aerostep: fconcn_org,fconca_org',fconcn_org,fconca_org
    WRITE(6,*) 'aerostep: orgpoa1rat,orgpoa2rat,orgpoa3rat',orgpoa1rat,orgpoa2rat,orgpoa3rat
    WRITE(6,*) 'aerostep: organt0rat,organt1rat,organt2rat,organt3rat ',organt0rat,organt1rat,organt2rat,organt3rat
    WRITE(6,*) 'aerostep: orgbio1rat,orgbio2rat,orgbio3rat',orgbio1rat,orgbio2rat,orgbio3rat
    icall=icall+1
end if

! RAR: for NH4, NO3 it's done in eql4

! *** sulfate:
        mstrnsfr = cblk(l,vso4ai)*factrans
        prod = deltaso4a(l)*fconcn(l)/dt + dmdt(l) ! Condensed mass +
        pol = prod*lossinv

        cblk(l,vso4ai) = pol + (cblk(l,vso4ai)-pol)*expdt
        cblk(l,vso4ai) = max(aeroconcmin,cblk(l,vso4ai))
        cblk(l,vso4aj) = cblk(l,vso4aj) + deltaso4a(l)*fconca(l) + mstrnsfr

! RAR: Debugging:
if (diagn_opt .AND. icall<200) then
    WRITE(6,*) 'aerostep, 1st: cblk(l,vpoa0i)+cblk(l,vpoa0j)', cblk(l,vpoa0i)+cblk(l,vpoa0j)
end if

! *** Non-volatile POA at bin-0, no condensational growth possible for this species
        mstrnsfr = cblk(l,vpoa0i)*factrans
        prod = eorgi
        pol = prod*lossinv

        cblk(l,vpoa0i) = pol + (cblk(l,vpoa0i)-pol)*expdt
        cblk(l,vpoa0i) = max(conmin,cblk(l,vpoa0i))
        cblk(l,vpoa0j) = cblk(l,vpoa0j) + eorgj*dt + mstrnsfr

! RAR: Debugging:
if (diagn_opt .AND. icall<200) then
    WRITE(6,*) 'aerostep, 2nd: cblk(l,vpoa0i)+cblk(l,vpoa0j)', cblk(l,vpoa0i)+cblk(l,vpoa0j)
end if

! RAR: Debugging:
if (diagn_opt .AND. icall<200) then
    WRITE(6,*) 'aerostep, 1st: cblk(l,vpoa1i)+cblk(l,vpoa1j)', cblk(l,vpoa1i)+cblk(l,vpoa1j)
end if

! Semi-volatile POA:
        mstrnsfr = cblk(l,vpoa1i)*factrans
        prod = orgpoa1rat(l)*fconcn_org(l)  ! fconcn_org from NUCLCOND
        pol = prod*lossinv

        cblk(l,vpoa1i) = pol + (cblk(l,vpoa1i)-pol)*expdt
        cblk(l,vpoa1i) = max(conmin,cblk(l,vpoa1i))
        cblk(l,vpoa1j) = cblk(l,vpoa1j) + orgpoa1rat(l)*fconca_org(l)*dt + mstrnsfr

! RAR: Debugging:
if (diagn_opt .AND. icall<200) then
    WRITE(6,*) 'aerostep, 2nd: cblk(l,vpoa1i)+cblk(l,vpoa1j) ', cblk(l,vpoa1i)+cblk(l,vpoa1j)
    WRITE(6,*) 'aerostep: mstrnsfr,pol', mstrnsfr,pol
end if

        mstrnsfr = cblk(l,vpoa2i)*factrans
        prod = orgpoa2rat(l)*fconcn_org(l)  ! fconcn_org from NUCLCOND
        pol = prod*lossinv

        cblk(l,vpoa2i) = pol + (cblk(l,vpoa2i)-pol)*expdt
        cblk(l,vpoa2i) = max(conmin,cblk(l,vpoa2i))
        cblk(l,vpoa2j) = cblk(l,vpoa2j) + orgpoa2rat(l)*fconca_org(l)*dt + mstrnsfr
        !
        mstrnsfr = cblk(l,vpoa3i)*factrans
        prod = orgpoa3rat(l)*fconcn_org(l)  ! fconcn_org from NUCLCOND
        pol = prod*lossinv

        cblk(l,vpoa3i) = pol + (cblk(l,vpoa1i)-pol)*expdt
        cblk(l,vpoa3i) = max(conmin,cblk(l,vpoa3i))
        cblk(l,vpoa3j) = cblk(l,vpoa3j) + orgpoa3rat(l)*fconca_org(l)*dt + mstrnsfr
!----------------------------------------------------------------------------------------

! RAR: Debugging:
if (diagn_opt .AND. icall<200) then
    WRITE(6,*) 'aerostep, 1st: cblk(l,vasoa0i)+cblk(l,vasoa0j) ', cblk(l,vasoa0i)+cblk(l,vasoa0j)
    WRITE(6,*) 'aerostep, 1st: cblk(l,vasoa1i)+cblk(l,vasoa1j) ', cblk(l,vasoa1i)+cblk(l,vasoa1j)
end if

! ASOA: bin-0
        mstrnsfr = cblk(l,vasoa0i)*factrans
        prod = organt0rat(l)*fconcn_org(l)
        pol = prod*lossinv

        cblk(l,vasoa0i) = pol + (cblk(l,vasoa0i)-pol)*expdt
        cblk(l,vasoa0i) = max(conmin,cblk(l,vasoa0i))
        cblk(l,vasoa0j) = cblk(l,vasoa0j) + organt0rat(l)*fconca_org(l)*dt + mstrnsfr
        !!!!!!!!!!!!!

! RAR: Debugging:
if (diagn_opt .AND. icall<200) then
    WRITE(6,*) 'aerostep, 2nd: cblk(l,vasoa0i)+cblk(l,vasoa0j)', cblk(l,vasoa0i)+cblk(l,vasoa0j)
    WRITE(6,*) 'aerostep: mstrnsfr', mstrnsfr
end if
        
! bin-1:
        mstrnsfr = cblk(l,vasoa1i)*factrans
        prod = organt1rat(l)*fconcn_org(l)
        pol = prod*lossinv

        cblk(l,vasoa1i) = pol + (cblk(l,vasoa1i)-pol)*expdt
        cblk(l,vasoa1i) = max(conmin,cblk(l,vasoa1i))
        cblk(l,vasoa1j) = cblk(l,vasoa1j) + organt1rat(l)*fconca_org(l)*dt + mstrnsfr
        !!!!!!!!!!!!!

! RAR: Debugging:
if (diagn_opt .AND. icall<200) then
    WRITE(6,*) 'aerostep, 2nd: cblk(l,vasoa1i)+cblk(l,vasoa1j)', cblk(l,vasoa1i)+cblk(l,vasoa1j)
    WRITE(6,*) 'aerostep: mstrnsfr', mstrnsfr
end if
        mstrnsfr = cblk(l,vasoa2i)*factrans
        prod = organt2rat(l)*fconcn_org(l)
        pol = prod*lossinv

        cblk(l,vasoa2i) = pol + (cblk(l,vasoa2i)-pol)*expdt
        cblk(l,vasoa2i) = max(conmin,cblk(l,vasoa2i))
        cblk(l,vasoa2j) = cblk(l,vasoa2j) + organt2rat(l)*fconca_org(l)*dt + mstrnsfr
        !!!!!!!!!!!!!

        mstrnsfr = cblk(l,vasoa3i)*factrans
        prod = organt3rat(l)*fconcn_org(l)
        pol = prod*lossinv

        cblk(l,vasoa3i) = pol + (cblk(l,vasoa3i)-pol)*expdt
        cblk(l,vasoa3i) = max(conmin,cblk(l,vasoa3i))
        cblk(l,vasoa3j) = cblk(l,vasoa3j) + organt3rat(l)*fconca_org(l)*dt + mstrnsfr
!----------------------------------------------------------------------------------------
! BSOA:
        mstrnsfr = cblk(l,vbsoa1i)*factrans
        prod = orgbio1rat(l)*fconcn_org(l)
        pol = prod*lossinv

        cblk(l,vbsoa1i) = pol + (cblk(l,vbsoa1i)-pol)*expdt
        cblk(l,vbsoa1i) = max(conmin,cblk(l,vbsoa1i))
        cblk(l,vbsoa1j) = cblk(l,vbsoa1j) + orgbio1rat(l)*fconca_org(l)*dt + mstrnsfr
        !!!!!!!!!!!!!

        mstrnsfr = cblk(l,vbsoa2i)*factrans
        prod = orgbio2rat(l)*fconcn_org(l)
        pol = prod*lossinv

        cblk(l,vbsoa2i) = pol + (cblk(l,vbsoa2i)-pol)*expdt
        cblk(l,vbsoa2i) = max(conmin,cblk(l,vbsoa2i))
        cblk(l,vbsoa2j) = cblk(l,vbsoa2j) + orgbio2rat(l)*fconca_org(l)*dt + mstrnsfr
        !!!!!!!!!!!!!

        mstrnsfr = cblk(l,vbsoa3i)*factrans
        prod = orgbio3rat(l)*fconcn_org(l)
        pol = prod*lossinv

        cblk(l,vbsoa3i) = pol + (cblk(l,vbsoa3i)-pol)*expdt
        cblk(l,vbsoa3i) = max(conmin,cblk(l,vbsoa3i))
        cblk(l,vbsoa3j) = cblk(l,vbsoa3j) + orgbio3rat(l)*fconca_org(l)*dt + mstrnsfr
!----------------------------------------------------------------------------------------



! *** other anthropogenic PM2.5
        mstrnsfr = cblk(l,vp25ai)*factrans
        prod = epm25i
        pol = prod*lossinv

        cblk(l,vp25ai) = pol + (cblk(l,vp25ai)-pol)*expdt
        cblk(l,vp25ai) = max(conmin,cblk(l,vp25ai))
        cblk(l,vp25aj) = cblk(l,vp25aj) + epm25j*dt + mstrnsfr

! ***  elemental carbon
        mstrnsfr = cblk(l,veci)*factrans
        prod = eeci
        pol = prod*lossinv

        cblk(l,veci) = pol + (cblk(l,veci)-pol)*expdt
        cblk(l,veci) = max(conmin,cblk(l,veci))
        cblk(l,vecj) = cblk(l,vecj) + eecj*dt + mstrnsfr

if (diagn_opt .AND. icall<200) then
    WRITE(6,*) 'aerostep: cblk(l,vecj),cblk(l,veci) ', cblk(l,vecj),cblk(l,veci)
    WRITE(6,*) 'aerostep: mstrnsfr, prod ', mstrnsfr, prod
end if

! *** coarse mode
! *** soil dust
        cblk(l,vsoila) = cblk(l,vsoila) + esoil(l)*dt
        cblk(l,vsoila) = max(conmin,cblk(l,vsoila))

! *** sea salt
        cblk(l,vseas) = cblk(l,vseas) + eseas(l)*dt
        cblk(l,vseas) = max(conmin,cblk(l,vseas))

! *** anthropogenic PM10 coarse fraction
        cblk(l,vantha) = cblk(l,vantha) + epmcoarse(l)*dt
        cblk(l,vantha) = max(conmin,cblk(l,vantha))
!ts *** secondary pollen
!        mstrnsfr = cblk(l,vpols)*factrans   !ts
!        prod = epols(l)                      !ts
!        pol = prod*lossinv                    !ts
!ts *** secondary pollen
!        cblk(l,vpols) = pol + (cblk(l,vpols)-pol)*expdt       !ts
!        cblk(l,vpols) = max(conmin,cblk(l,vpols))               !ts
!        cblk(l,vpols) = cblk(l,vpols) + epols(l)*dt + mstrnsfr  !ts
!ts *** primary pollen
!        cblk(l,vpolp) = cblk(l,vpolp) + epolp(l)*dt !ts
!        cblk(l,vpolp) = max(conmin,cblk(l,vpolp))  !ts

!! -- THIS IS FROM MODULE_AEROSOLS_SORGAM -- why the difference????
!ts *** primary pollen
        cblk(l,vpolp) = cblk(l,vpolp) + epolp(l)*dt !ts
        cblk(l,vpolp) = max(conmin,cblk(l,vpolp))  !ts
!ts *** secondary pollen
        cblk(l,vpols) = cblk(l,vpols) + epols(l)*dt !ts
        cblk(l,vpols) = max(conmin,cblk(l,vpols))  !ts
      END DO    ! numcells

! *** Check for mode merging,if Aitken mode is growing faster than j-mod
!     then merge modes by renaming.

! *** use Binkowski-Kreidenweis paradigm, now including emissions

! end of time-step loop for total mass                 
      DO lcell = 1, numcells

!       IF( CGRN3(LCELL) .GT. CGRA3(LCELL) .AND.
!    &      CBLK(LCELL,VNU0) .GT. CBLK(LCELL,VAC0) ) THEN ! check if mer
        IF (cgrn3(lcell)>cgra3(lcell) .OR. dgnuc(lcell)>.03E-6 .AND. cblk( &
            lcell,vnu0)>cblk(lcell,vac0)) THEN

! check if mer
          aaa = getaf(cblk(lcell,vnu0),cblk(lcell,vac0),dgnuc(lcell), &
            dgacc(lcell),xxlsgn,xxlsga,sqrt2)

! *** AAA is the value of ln( dd / DGNUC ) / ( SQRT2 * XXLSGN ), where
!        dd is the diameter at which the Aitken-mode and accumulation-mo
!        distributions intersect (overap).

          xnum = max(aaa,xxm3)                                    ! this means that no more than one ha
                                   ! total Aitken mode number may be tra per call.

! do not let XNUM become negative bec
          xm3 = xnum - xxm3
! set up for 3rd moment and mass tran
          IF (xm3>0.0) THEN
! do mode merging if  overlap is corr
            phnum = 0.5*(1.0+erf(xnum))
            phm3 = 0.5*(1.0+erf(xm3))
            fnum = 0.5*erfc(xnum)
            fm3 = 0.5*erfc(xm3)

if (diagn_opt .AND. icall<200) then
    WRITE(6,*) 'aerostep: fm3 ',fm3
end if
            
!     In the Aitken mode:

! *** FNUM and FM3 are the fractions of the number and 3rd moment
!     distributions with  diameters greater than dd respectively.

! *** PHNUM and PHM3 are the fractions of the number and 3rd moment
!     distributions with diameters less than dd.

! *** rename the  Aitken mode particle number as accumulation mode
!     particle number

    cblk(lcell,vac0) = cblk(lcell,vac0) + fnum*cblk(lcell,vnu0)

! *** adjust the Aitken mode number

    cblk(lcell,vnu0) = phnum*cblk(lcell,vnu0)

! *** Rename mass from Aitken mode to acumulation mode. The mass transfer
!     to the accumulation mode is proportional to the amount of 3rd moment
!     transferred, therefore FM3 is used for mass transfer.

    cblk(lcell,vso4aj) = cblk(lcell,vso4aj) + cblk(lcell,vso4ai)*fm3

    cblk(lcell,vnh4aj) = cblk(lcell,vnh4aj) + cblk(lcell,vnh4ai)*fm3

    cblk(lcell,vno3aj) = cblk(lcell,vno3aj) + cblk(lcell,vno3ai)*fm3
    cblk(lcell,vnaaj) = cblk(lcell,vnaaj) + cblk(lcell,vnaai)*fm3
    cblk(lcell,vclaj) = cblk(lcell,vclaj) + cblk(lcell,vclai)*fm3
!  if ( do_isorropia ) then
!    cblk(lcell,vcaaj) = cblk(lcell,vcaaj) + cblk(lcell,vcaai)*fm3
!    cblk(lcell,vkaj)  = cblk(lcell,vkaj)  + cblk(lcell,vkai)*fm3
!    cblk(lcell,vmgaj) = cblk(lcell,vmgaj) + cblk(lcell,vmgai)*fm3
!  endif

! POA:
    cblk(lcell,vpoa0j) = cblk(lcell,vpoa0j) + cblk(lcell,vpoa0i)*fm3

    cblk(lcell,vpoa1j) = cblk(lcell,vpoa1j) + cblk(lcell,vpoa1i)*fm3

    cblk(lcell,vpoa2j) = cblk(lcell,vpoa2j) + cblk(lcell,vpoa2i)*fm3

    cblk(lcell,vpoa3j) = cblk(lcell,vpoa3j) + cblk(lcell,vpoa3i)*fm3

! RAR: Debugging:
if (diagn_opt .AND. icall<200) then
    WRITE(6,*) 'aerostep: fm3,vpoa0j,vpoa1j,vpoa2j,vpoa3j ', fm3,vpoa0j,vpoa1j,vpoa2j,vpoa3j
    WRITE(6,*) 'aerostep: cblk(lcell,vpoa3j),cblk(lcell,vpoa3i) ',cblk(lcell,vpoa3j),cblk(lcell,vpoa3i)
end if

! ASOA:
    cblk(lcell,vasoa0j) = cblk(lcell,vasoa0j) + cblk(lcell,vasoa0i)*fm3

    cblk(lcell,vasoa1j) = cblk(lcell,vasoa1j) + cblk(lcell,vasoa1i)*fm3

    cblk(lcell,vasoa2j) = cblk(lcell,vasoa2j) + cblk(lcell,vasoa2i)*fm3

    cblk(lcell,vasoa3j) = cblk(lcell,vasoa3j) + cblk(lcell,vasoa3i)*fm3

!    cblk(lcell,vasoa4j) = cblk(lcell,vasoa4j) + cblk(lcell,vasoa4i)*fm3

! BSOA:
    cblk(lcell,vbsoa1j) = cblk(lcell,vbsoa1j) + cblk(lcell,vbsoa1i)*fm3

    cblk(lcell,vbsoa2j) = cblk(lcell,vbsoa2j) + cblk(lcell,vbsoa2i)*fm3

    cblk(lcell,vbsoa3j) = cblk(lcell,vbsoa3j) + cblk(lcell,vbsoa3i)*fm3

!    cblk(lcell,vbsoa4j) = cblk(lcell,vbsoa4j) + cblk(lcell,vbsoa4i)*fm3

!    cblk(lcell,vorgpaj) = cblk(lcell,vorgpaj) + cblk(lcell,vorgpai)*fm3

    cblk(lcell,vp25aj) = cblk(lcell,vp25aj) + cblk(lcell,vp25ai)*fm3

    cblk(lcell,vecj) = cblk(lcell,vecj) + cblk(lcell,veci)*fm3

    ! RAR: Debugging:
if (diagn_opt .AND. icall<200) then
    WRITE(6,*) 'aerostep: fm3,vecj,veci ', fm3,vecj,veci
    WRITE(6,*) 'aerostep: cblk(lcell,vecj),cblk(lcell,veci) ',cblk(lcell,vecj),cblk(lcell,veci)
end if
! *** update Aitken mode for mass loss to accumulation mode
          cblk(lcell,vso4ai) = cblk(lcell,vso4ai)*phm3

          cblk(lcell,vnh4ai) = cblk(lcell,vnh4ai)*phm3

          cblk(lcell,vno3ai) = cblk(lcell,vno3ai)*phm3

          cblk(lcell,vclai) = cblk(lcell,vclai)*phm3
          cblk(lcell,vnaai) = cblk(lcell,vnaai)*phm3
!       if ( do_isorropia ) then
!          cblk(lcell,vcaai) = cblk(lcell,vcaai)*phm3
!          cblk(lcell,vkai)  = cblk(lcell,vkai)*phm3
!          cblk(lcell,vmgai) = cblk(lcell,vmgai)*phm3
!       endif
! POA
          cblk(lcell,vpoa0i) = cblk(lcell,vpoa0i)*phm3

          cblk(lcell,vpoa1i) = cblk(lcell,vpoa1i)*phm3

          cblk(lcell,vpoa2i) = cblk(lcell,vpoa2i)*phm3

          cblk(lcell,vpoa3i) = cblk(lcell,vpoa3i)*phm3
! ASOA
          cblk(lcell,vasoa0i) = cblk(lcell,vasoa0i)*phm3

          cblk(lcell,vasoa1i) = cblk(lcell,vasoa1i)*phm3

          cblk(lcell,vasoa2i) = cblk(lcell,vasoa2i)*phm3

          cblk(lcell,vasoa3i) = cblk(lcell,vasoa3i)*phm3

!          cblk(lcell,vasoa4i) = cblk(lcell,vasoa4i)*phm3
! BSOA
          cblk(lcell,vbsoa1i) = cblk(lcell,vbsoa1i)*phm3

          cblk(lcell,vbsoa2i) = cblk(lcell,vbsoa2i)*phm3

          cblk(lcell,vbsoa3i) = cblk(lcell,vbsoa3i)*phm3

!          cblk(lcell,vbsoa4i) = cblk(lcell,vbsoa4i)*phm3

!          cblk(lcell,vorgpai) = cblk(lcell,vorgpai)*phm3

          cblk(lcell,vp25ai) = cblk(lcell,vp25ai)*phm3

          cblk(lcell,veci) = cblk(lcell,veci)*phm3

    END IF      ! end check on whether modal overlap is OK

   END IF  ! end check on necessity for merging

END DO

!     set min value for all concentrations
! loop for merging
      DO spc = 1, nspcsda
        DO lcell = 1, numcells
          cblk(lcell,spc) = max(cblk(lcell,spc),conmin)
        END DO
      END DO
!---------------------------------------------------------------------------------

RETURN
END SUBROUTINE aerostep
!#######################################################################

SUBROUTINE awater(irhx,mso4,mnh4,mno3,wh2o)
! NOTE!!! wh2o is returned in micrograms / cubic meter
!         mso4,mnh4,mno3 are in microMOLES / cubic meter

!  This  version uses polynomials rather than tables, and uses empirical
! polynomials for the mass fraction of solute (mfs) as a function of wat
!   where:

!            mfs = ms / ( ms + mw)
!             ms is the mass of solute
!             mw is the mass of water.

!  Define y = mw/ ms

!  then  mfs = 1 / (1 + y)

!    y can then be obtained from the values of mfs as

!             y = (1 - mfs) / mfs


!     the aerosol is assumed to be in a metastable state if the rh is
!     is below the rh of deliquescence, but above the rh of crystallizat

!     ZSR interpolation is used for sulfates with x ( the molar ratio of
!     ammonium to sulfate in eh range 0 <= x <= 2, by sections.
!     section 1: 0 <= x < 1
!     section 2: 1 <= x < 1.5
!     section 3: 1.5 <= x < 2.0
!     section 4: 2 <= x
!     In sections 1 through 3, only the sulfates can affect the amount o
!     on the particles.
!     In section 4, we have fully neutralized sulfate, and extra ammoniu
!     allows more nitrate to be present. Thus, the ammount of water is c
!     using ZSR for ammonium sulfate and ammonium nitrate. Crystallizati
!     assumed to occur in sections 2,3,and 4. See detailed discussion be

! definitions:
!     mso4, mnh4, and mno3 are the number of micromoles/(cubic meter of
!      for sulfate, ammonium, and nitrate respectively
!     irhx is the relative humidity (%)
!     wh2o is the returned water amount in micrograms / cubic meter of a
!     x is the molar ratio of ammonium to sulfate
!     y0,y1,y1.5, y2 are the water contents in mass of water/mass of sol
!     for pure aqueous solutions with x equal 1, 1.5, and 2 respectively
!     y3 is the value of the mass ratio of water to solute for
!     a pure ammonium nitrate  solution.

!coded by Dr. Francis S. Binkowski, 4/8/96.

!     IMPLICIT NONE
      INTEGER irhx, irh
      REAL mso4, mnh4, mno3, mw2, mwano3
      REAL tso4, tnh4, tno3, wh2o, x
      REAL aw, awc
!     REAL poly4, poly6
      REAL mfs0, mfs1, mfs15, mfs2
      REAL c0(4), c1(4), c15(4), c2(4)
      REAL y, y0, y1, y15, y2, y3, y40, y140, y1540, yc
      REAL kso4(6), kno3(6), mfsso4, mfsno3
!      REAL mwso4, mwnh4, mwno3, mw2, mwano3

! *** molecular weights:
      PARAMETER ( mw2=mwso4+2.0*mwnh4, mwano3=mwno3+mwnh4 )

!     The polynomials use data for aw as a function of mfs from Tang and
!     Munkelwitz, JGR 99: 18801-18808, 1994.
!     The polynomials were fit to Tang's values of water activity as a
!     function of mfs.

! *** coefficients of polynomials fit to Tang and Munkelwitz data
!     now give mfs as a function of water activity.

      DATA c1 /0.9995178, -0.7952896, 0.99683673, -1.143874/
      DATA c15 /1.697092, -4.045936, 5.833688, -3.463783/
      DATA c2 /2.085067, -6.024139, 8.967967, -5.002934/

! *** the following coefficients are a fit to the data in Table 1 of
!     Nair & Vohra, J. Aerosol Sci., 6: 265-271, 1975
!      data c0/0.8258941, -1.899205, 3.296905, -2.214749 /
! *** New data fit to data from
!       Nair and Vohra J. Aerosol Sci., 6: 265-271, 1975
!       Giaque et al. J.Am. Chem. Soc., 82: 62-70, 1960
!       Zeleznik J. Phys. Chem. Ref. Data, 20: 157-1200
      DATA c0 /0.798079, -1.574367, 2.536686, -1.735297/

! *** polynomials for ammonium nitrate and ammonium sulfate are from:
!     Chan et al.1992, Atmospheric Environment (26A): 1661-1673.

      DATA kno3 /0.2906, 6.83665, -26.9093, 46.6983, -38.803, 11.8837/
      DATA kso4 /2.27515, -11.147, 36.3369, -64.2134, 56.8341, -20.0953/

! *** check range of per cent relative humidity
      irh = irhx
      irh = max(1,irh)
      irh = min(irh,100)
      aw = float(irh)/ & ! water activity = fractional relative h
        100.0
      tso4 = max(mso4,0.0)
      tnh4 = max(mnh4,0.0)
      tno3 = max(mno3,0.0)
      x = 0.0
! *** if there is non-zero sulfate calculate the molar ratio
      IF (tso4>0.0) THEN
        x = tnh4/tso4
      ELSE
! *** otherwise check for non-zero nitrate and ammonium
        IF (tno3>0.0 .AND. tnh4>0.0) x = 10.0
      END IF

! *** begin screen on x for calculating wh2o
      IF (x<1.0) THEN

        mfs0 = poly4(c0,aw)
        mfs1 = poly4(c1,aw)
        y0 = (1.0-mfs0)/mfs0
        y1 = (1.0-mfs1)/mfs1
        y = (1.0-x)*y0 + x*y1

      ELSE IF (x<1.5) THEN

        IF (irh>=40) THEN
          mfs1 = poly4(c1,aw)
          mfs15 = poly4(c15,aw)
          y1 = (1.0-mfs1)/mfs1
          y15 = (1.0-mfs15)/mfs15
          y = 2.0*(y1*(1.5-x)+y15*(x-1.0))
        ELSE
! *** set up for crystalization

! *** Crystallization is done as follows:
!      For 1.5 <= x, crystallization is assumed to occur at rh = 0.4
!      For x <= 1.0, crystallization is assumed to occur at an rh < 0.01
!      and since the code does not allow ar rh < 0.01, crystallization
!      is assumed not to occur in this range.
!      For 1.0 <= x <= 1.5 the crystallization curve is a straignt line
!      from a value of y15 at rh = 0.4 to a value of zero at y1. From
!      point B to point A in the diagram.
!      The algorithm does a double interpolation to calculate the amount
!      water.

!        y1(0.40)               y15(0.40)
!         +                     + Point B

!         +--------------------+
!       x=1                   x=1.5
!      Point A

          awc = 0.80*(x-1.0) ! rh along the crystallization curve.
          y = 0.0
          IF (aw>=awc) & ! interpolate using crystalization 
              THEN
            mfs1 = poly4(c1,0.40)
            mfs15 = poly4(c15,0.40)
            y140 = (1.0-mfs1)/mfs1
            y1540 = (1.0-mfs15)/mfs15
            y40 = 2.0*(y140*(1.5-x)+y1540*(x-1.0))
            yc = 2.0*y1540*(x-1.0) ! y along crystallization cur
            y = y40 - (y40-yc)*(0.40-aw)/(0.40-awc)
! end of checking for aw                             
          END IF

        END IF
! end of checking on irh                               
      ELSE IF (x<1.9999) THEN

        y = 0.0
        IF (irh>=40) THEN
          mfs15 = poly4(c15,aw)
          mfs2 = poly4(c2,aw)
          y15 = (1.0-mfs15)/mfs15
          y2 = (1.0-mfs2)/mfs2
          y = 2.0*(y15*(2.0-x)+y2*(x-1.5))

        END IF

! end of check for crystallization

      ELSE
! regime where ammonium sulfate and ammonium nitrate are in solution.

! *** following cf&s for both ammonium sulfate and ammonium nitrate
! *** check for crystallization here. their data indicate a 40% value
!     is appropriate.
! 1.9999 < x                                                 
        y2 = 0.0
        y3 = 0.0
        IF (irh>=40) THEN
          mfsso4 = poly6(kso4,aw)
          mfsno3 = poly6(kno3,aw)
          y2 = (1.0-mfsso4)/mfsso4
          y3 = (1.0-mfsno3)/mfsno3

        END IF

      END IF
! *** now set up output of wh2o

!      wh2o units are micrograms (liquid water) / cubic meter of air

! end of checking on x                                    
      IF (x<1.9999) THEN

        wh2o = y*(tso4*mwso4+mwnh4*tnh4)

      ELSE

! *** this is the case that all the sulfate is ammonium sulfate
!     and the excess ammonium forms ammonum nitrate

        wh2o = y2*tso4*mw2 + y3*tno3*mwano3

      END IF

      RETURN
    END SUBROUTINE awater
!//////////////////////////////////////////////////////////////////////

    SUBROUTINE coagrate(cblk,blkta,pdensn,pdensa,amu, &
                        dgnuc,dgacc,knnuc,knacc,urn00,ura00,brna01,c30)
!***********************************************************************
!**    DESCRIPTION:  calculates aerosol coagulation rates for unimodal
!       and bimodal coagulation using E. Whitby 1990's prescription.

!.......   Rates for coagulation:
!.......   Unimodal Rates:
!.......   URN00:  nuclei       mode 0th moment self-coagulation rate
!.......   URA00:  accumulation mode 0th moment self-coagulation rate

!.......   Bimodal Rates:  (only 1st order coeffs appear)
!.......   NA-- nuclei  with accumulation coagulation rates,
!.......   AN-- accumulation with nuclei coagulation rates
!.......   BRNA01:  rate for 0th moment ( d(nuclei mode 0) / dt  term)
!.......   BRNA31:           3rd        ( d(nuclei mode 3) / dt  term)
!**    Revision history:
!       prototype 1/95 by Uma and Carlie
!       Revised   8/95 by US for calculation of density from stmt func
!                 and collect met variable stmt funcs in one include fil
!      REVISED 7/25/96 by FSB to use block structure
!      REVISED 9/13/96 BY FSB for Uma's FIXEDBOTH case only.
!      REVISED 11/08/96 BY FSB the Whitby Shankar convention on signs
!                              changed. All coagulation coefficients
!                              returned with positive signs. Their
!                              linearization is also abandoned.
!                              Fixed values are used for the corrections
!                              to the free-molecular coagulation integra
!                              The code forces the harmonic means to be
!                              evaluated in 64 bit arithmetic on 32 bit
!     REVISED 11/14/96 BY FSB  Internal units are now MKS, moment / unit

!      REVISED 1/12/98 by FSB   C30 replaces BRNA31 as an array. This wa
!                              because BRNA31 can become zero on a works
!                              because of limited precision. With the ch
!                              aerostep to omit update of the 3rd moment
!                              C30 is the only variable now needed.
!                              the logic using ONE88 to force REAL*8 ari
!                              has been removed and all intermediates ar
!                              REAL*8.
!     IMPLICIT NONE

! dimension of arrays             
!      INTEGER blksize
! actual number of cells in arrays
!      INTEGER numcells
!      INTEGER nspcsda

! nmber of species in CBLK        
      REAL cblk(blksize,nspcsda) ! main array of variables         
      REAL blkta(blksize) ! Air temperature [ K ]           
      REAL pdensn(blksize) ! average particel density in Aitk
      REAL pdensa(blksize) ! average particel density in accu
      REAL amu(blksize) ! atmospheric dynamic viscosity [ 
      REAL dgnuc(blksize) ! Aitken mode mean diameter [ m ] 
      REAL dgacc(blksize) ! accumulation mode mean diameter 
      REAL knnuc(blksize) ! Aitken mode Knudsen number      
      REAL knacc(blksize) 
! *** output:

! accumulation mode Knudsen number
      REAL urn00(blksize) ! intramodal coagulation rate (Ait
      REAL ura00(blksize) 
! intramodal coagulation rate (acc
      REAL brna01(blksize) ! intermodal coagulaton rate (numb
      REAL c30(blksize)                                                               ! by inter

! *** Local variables:
! intermodal 3rd moment transfer r
      REAL*8 kncnuc, kncacc  ! coeffs for unimodal NC coag rate
      REAL*8 kfmnuc, kfmacc  ! coeffs for unimodal FM coag rate
      REAL*8 knc, kfm        ! coeffs for bimodal NC, FM coag rate
      REAL*8 bencnn, bencna  ! NC 0th moment coag rate (both modes)
      REAL*8 bencm3n         ! NC 3rd moment coag rate (nuc mode)
      REAL*8 befmnn, befmna  ! FM 0th moment coag rate (both modes)
      REAL*8 befm3n          ! FM 3rd moment coag rate (nuc mode)
      REAL*8 betann, betana  ! composite coag rates, mom 0 (both mode
      REAL*8 brna31          ! intermodal coagulation rate for 3rd mode
      REAL*8 s1              ! scratch subexpression
      REAL*8 t1, t2          ! scratch subexpressions
      REAL*8 t16, t26        ! T1**6, T2**6
      REAL*8 rat, rin        ! ratio of acc to nuc size and its inver
      REAL*8 rsqt,rsq4       ! sqrt( rat ), rsqt**4
      REAL*8 rsqti,rsqi3     ! sqrt( 1/rat ), sqrt( 1/rat**3 )
      REAL*8 dgn3            ! dgnuc**3
      REAL*8 dga3            ! in 64 bit arithmetic
! dgacc**3

      INTEGER lcell
! *** Fixed values for correctionss to coagulation
!      integrals for free-molecular case.
! loop counter                                      
      REAL*8, PARAMETER :: bm0=0.8D0
      REAL*8, PARAMETER :: bm0i=0.9D0
      REAL*8, PARAMETER :: bm3i=0.9D0
      REAL*8, PARAMETER :: a=1.246D0   ! approx Cunningham corr. factor

!   begin body of subroutine  COAGRATE
!...........   Main computational grid-traversal loops
!...........   for computing coagulation rates.

! *** Both modes have fixed std devs.
      DO lcell = 1,numcells
! *** moment independent factors

!  loop on LCELL               
        s1 = two3*boltz*blkta(lcell)/amu(lcell)

! For unimodal coagualtion:
        kncnuc = s1
        kncacc = s1

        kfmnuc = sqrt(3.0*boltz*blkta(lcell)/pdensn(lcell))
        kfmacc = sqrt(3.0*boltz*blkta(lcell)/pdensa(lcell))

! For bimodal coagulation:
        knc = s1
        kfm = sqrt(6.0*boltz*blkta(lcell)/(pdensn(lcell)+pdensa(lcell)))

!...........   Begin unimodal coagulation rate calculations:
!...........   Near-continuum regime.
        dgn3 = dgnuc(lcell)**3
        dga3 = dgacc(lcell)**3

        t1 = sqrt(dgnuc(lcell))
        t2 = sqrt(dgacc(lcell))
        t16 = & ! = T1**6                               
          dgn3
        t26 = & 
          dga3
!.......   Note rationalization of fractions and subsequent cancellation
!.......   from the formulation in  Whitby et al. (1990)

! = T2**6                               
        bencnn = kncnuc*(1.0+esn08+a*knnuc(lcell)*(esn04+esn20))
        bencna = kncacc*(1.0+esa08+a*knacc(lcell)*(esa04+esa20))

!...........   Free molecular regime. Uses fixed value for correction
!               factor BM0

        befmnn = kfmnuc*t1*(en1+esn25+2.0*esn05)*bm0
        befmna = kfmacc*t2*(ea1+esa25+2.0*esa05)*bm0

!...........   Calculate half the harmonic mean between unimodal rates
!...........   free molecular and near-continuum regimes

! FSB       64 bit evaluation

        betann = bencnn*befmnn/(bencnn+befmnn)
        betana = bencna*befmna/(bencna+befmna)

        urn00(lcell) = betann
        ura00(lcell) = betana

! *** End of unimodal coagulation calculations.

!...........   Begin bimodal coagulation rate calculations:

        rat = dgacc(lcell)/dgnuc(lcell)
        rin = 1.0D0/rat
        rsqt = sqrt(rat)
        rsq4 = rat**2

        rsqti = 1.0D0/rsqt
        rsqi3 = rin*rsqti

!...........   Near-continuum coeffs:
!...........   0th moment nuc mode bimodal coag coefficient

        bencnn = knc*(2.0+a*knnuc(lcell)*(esn04+rat*esn16*esa04)+a*knacc(lcell &
          )*(esa04+rin*esa16*esn04)+(rat+rin)*esn04*esa04)

!...........   3rd moment nuc mode bimodal coag coefficient

        bencm3n = knc*dgn3*(2.0*esn36+a*knnuc(lcell)*(esn16+rat*esn04*esa04)+a &
          *knacc(lcell)*(esn36*esa04+rin*esn64*esa16)+rat*esn16*esa04+ &
          rin*esn64*esa04)

!...........   Free molecular regime coefficients:
!...........   Uses fixed value for correction
!               factor BM0I, BM3I
!...........   0th moment nuc mode coeff

        befmnn = kfm*bm0i*t1*(en1+rsqt*ea1+2.0*rat*en1*esa04+rsq4*esn09*esa16+ &
          rsqi3*esn16*esa09+2.0*rsqti*esn04*ea1)

!...........   3rd moment nuc mode coeff

        befm3n = kfm*bm3i*t1*t16*(esn49+rsqt*esn36*ea1+2.0*rat*esn25*esa04+ &
          rsq4*esn09*esa16+rsqi3*esn100*esa09+2.0*rsqti*esn64*ea1)

!...........   Calculate half the harmonic mean between bimodal rates
!...........   free molecular and near-continuum regimes
! FSB       Force 64 bit evaluation

        brna01(lcell) = bencnn*befmnn/(bencnn+befmnn)
        brna31 = bencm3n*befm3n/(bencm3n+befm3n)              ! BRNA31 now is a scala

        c30(lcell) = brna31*cblk(lcell,vac0)*cblk(lcell,vnu0)

      END DO
! end of main lop over cells                            
      RETURN
END SUBROUTINE coagrate
!------------------------------------------------------------------

! subroutine  to find the roots of a cubic equation / 3rd order polynomi
! formulae can be found in numer. recip.  on page 145
!   kiran  developed  this version on 25/4/1990
!   dr. francis binkowski modified the routine on 6/24/91, 8/7/97
! ***
!234567
! coagrate                                     
    SUBROUTINE cubic(a2,a1,a0,nr,crutes)
!     IMPLICIT NONE
      INTEGER nr
      REAL*8 a2, a1, a0
      REAL crutes(3)
      REAL*8 qq, rr, a2sq, theta, sqrt3, one3rd
      REAL*8 dum1, dum2, part1, part2, part3, rrsq, phi, yy1, yy2, yy3
      REAL*8 costh, sinth
      DATA sqrt3/1.732050808/, one3rd/0.333333333/
!bs
      REAL*8 onebs
      PARAMETER (onebs=1.0)
!bs
      a2sq = a2*a2
      qq = (a2sq-3.*a1)/9.
      rr = (a2*(2.*a2sq-9.*a1)+27.*a0)/54.
! CASE 1 THREE REAL ROOTS or  CASE 2 ONLY ONE REAL ROOT
      dum1 = qq*qq*qq
      rrsq = rr*rr
      dum2 = dum1 - rrsq
      IF (dum2>=0.) THEN
! NOW WE HAVE THREE REAL ROOTS
        phi = sqrt(dum1)
        IF (abs(phi)<1.E-20) THEN
          print *, ' cubic phi small, phi = ',phi
          crutes(1) = 0.0
          crutes(2) = 0.0
          crutes(3) = 0.0
          nr = 0
          CALL wrf_error_fatal ( 'PHI < CRITICAL VALUE')
        END IF
        theta = acos(rr/phi)/3.0
        costh = cos(theta)
        sinth = sin(theta)
! *** use trig identities to simplify the expressions
! *** binkowski's modification
        part1 = sqrt(qq)
        yy1 = part1*costh
        yy2 = yy1 - a2/3.0
        yy3 = sqrt3*part1*sinth
        crutes(3) = -2.0*yy1 - a2/3.0
        crutes(2) = yy2 + yy3
        crutes(1) = yy2 - yy3
! *** SET NEGATIVE ROOTS TO A LARGE POSITIVE VALUE
        IF (crutes(1)<0.0) crutes(1) = 1.0E9
        IF (crutes(2)<0.0) crutes(2) = 1.0E9
        IF (crutes(3)<0.0) crutes(3) = 1.0E9
! *** put smallest positive root in crutes(1)
        crutes(1) = min(crutes(1),crutes(2),crutes(3))
        nr = 3
!     NOW HERE WE HAVE ONLY ONE REAL ROOT
      ELSE
! dum IS NEGATIVE                                           
        part1 = sqrt(rrsq-dum1)
        part2 = abs(rr)
        part3 = (part1+part2)**one3rd
        crutes(1) = -sign(onebs,rr)*(part3+(qq/part3)) - a2/3.
!bs     &        -sign(1.0,rr) * ( part3 + (qq/part3) ) - a2/3.
        crutes(2) = 0.
        crutes(3) = 0.
!IAREV02...ADDITIONAL CHECK on NEGATIVE ROOTS
! *** SET NEGATIVE ROOTS TO A LARGE POSITIVE VALUE
!     if(crutes(1) .lt. 0.0) crutes(1) = 1.0e9
        nr = 1
      END IF
      RETURN
    END SUBROUTINE cubic
!///////////////////////////////////////////////////////////////////////
! RAR: The eql3 is removed, since it isn't used here
!///////////////////////////////////////////////////////////////////////

!///////////////////////////////////////////////////////////////////////

!liqy
!    Calculate the aerosol chemical speciation and water content.
! cubic

SUBROUTINE eql4(cblk,blkta,blkrh,diagn_opt)
!***********************************************************************
!**    DESCRIPTION:
!       Calculates the distribution of ammonia/ammonium, nitric acid/nitrate,
!       and water between the gas and aerosol phases as the total sulfate,
!       ammonia, and nitrate concentrations, relative humidity and
!       temperature change.  The evolution of the aerosol mass concentration
!       due to the change in aerosol chemical composition is calculated.
!**    REVISION HISTORY:
!       prototype 1/95 by Uma and Carlie
!       Revised   8/95 by US to calculate air density in stmt func
!                 and collect met variable stmt funcs in one include fil
!       Revised 7/26/96 by FSB to use block concept.
!       Revise 12/1896 to do do i-mode calculation.
!**********************************************************************
!     IMPLICIT NONE

! dimension of arrays             
!      INTEGER blksize
! actual number of cells in arrays
!      INTEGER numcells
! nmber of species in CBLK        
!      INTEGER nspcsda

LOGICAL :: diagn_opt

      REAL cblk(blksize,nspcsda) 
! *** Meteorological information in blocked arays:

! main array of variables         
      REAL blkta(blksize) ! Air temperature [ K ]                   
      REAL blkrh(blksize) 
! Fractional relative humidity

   INTEGER, SAVE :: icall
   INTEGER lcell

! loop counter
! air temperature                             
      REAL temp
!iamodels3
      REAL rh
! relative humidity                           
      REAL so4, no3, nh3, nh4, hno3
      REAL aso4, ano3, ah2o, anh4, gnh3, gno3
! Fraction of dry sulfate mass in i-mode         
      REAL fraci
!.......................................................................
      REAL fracj

      REAL ana,acl
! Fraction of dry sulfate mass in j-mode         
      DO lcell = 1,numcells
! *** Fetch temperature, fractional relative humidity, and
!     air density

!  loop on cells                    
        temp = blkta(lcell)
        rh = blkrh(lcell)

! *** the following is an interim procedure. Assume the i-mode has the
!     same relative mass concentrations as the total mass. Use SO4 as
!     the surrogate. The results of this should be the same as those
!     from the original RPM.

! *** do total aerosol
        so4 = cblk(lcell,vso4aj) + cblk(lcell,vso4ai)

        no3 = cblk(lcell,vno3aj) + cblk(lcell,vno3ai)

        hno3= cblk(lcell,vhno3)

        nh4 = cblk(lcell,vnh4aj) + cblk(lcell,vnh4ai)

        nh3 = cblk(lcell,vnh3)

! RAR: Debugging:
IF (diagn_opt .AND. icall<200) then
    WRITE(6,*) 'eql4: rh,temp', rh,temp
    WRITE(6,*) 'eql4: so4,hno3,no3,nh3,nh4', so4,hno3,no3,nh3,nh4
    WRITE(6,*) 'eql4: aso4,ano3,ah2o,anh4', so4,hno3,no3,nh3,nh4
    WRITE(6,*) 'eql4: gnh3,gno3', gnh3,gno3
ENDIF

 CALL rpmares_old(so4,hno3,no3,nh3,nh4,rh,temp,aso4,ano3,ah2o,anh4,gnh3,gno3,diagn_opt)

! RAR: Debugging:
IF (diagn_opt .AND. icall<200) then
    WRITE(6,*) 'eql4, after rpmares_old: rh,temp ', rh,temp
    WRITE(6,*) 'eql4, after rpmares_old: so4,hno3,no3,nh3,nh4 ', so4,hno3,no3,nh3,nh4
    WRITE(6,*) 'eql4, after rpmares_old: aso4,ano3,ah2o,anh4 ', so4,hno3,no3,nh3,nh4
    WRITE(6,*) 'eql4, after rpmares_old: gnh3,gno3 ', gnh3,gno3
    icall=icall+1
ENDIF

! *** get modal fraction
        fraci = cblk(lcell,vso4ai)/(cblk(lcell,vso4aj)+cblk(lcell,vso4ai))
        fracj = 1.0 - fraci

! *** update do i-mode
        cblk(lcell,vh2oai) = fraci*ah2o
        cblk(lcell,vnh4ai) = fraci*anh4
        cblk(lcell,vno3ai) = fraci*ano3
!        cblk(lcell,vnaai) = fraci*ana
!        cblk(lcell,vclai) = fraci*acl

! *** update accumulation mode:
        cblk(lcell,vh2oaj) = fracj*ah2o
        cblk(lcell,vnh4aj) = fracj*anh4
        cblk(lcell,vno3aj) = fracj*ano3
!        cblk(lcell,vnaaj) = fracj*ana
!        cblk(lcell,vclaj) = fracj*acl


! *** update gas / vapor phase
        cblk(lcell,vnh3) = gnh3
        cblk(lcell,vhno3)= gno3

      END DO
!  end loop on cells                     
      RETURN

    END SUBROUTINE eql4
!------------------------------------------------------------------------------

    REAL FUNCTION getaf(ni,nj,dgni,dgnj,xlsgi,xlsgj,sqrt2)
! *** set up new processor for renaming of particles from i to j modes
!     IMPLICIT NONE
      REAL aa, bb, cc, disc, qq, alfa, l, yji
      REAL ni, nj, dgni, dgnj, xlsgi, xlsgj, sqrt2

      alfa = xlsgi/xlsgj
      yji = log(dgnj/dgni)/(sqrt2*xlsgi)
      aa = 1.0 - alfa*alfa
      l = log(alfa*nj/ni)
      bb = 2.0*yji*alfa*alfa
      cc = l - yji*yji*alfa*alfa
      disc = bb*bb - 4.0*aa*cc
      IF (disc<0.0) THEN
        getaf = - & ! error in intersection                     
          5.0
        RETURN
      END IF
      qq = -0.5*(bb+sign(1.0,bb)*sqrt(disc))
      getaf = cc/qq
      RETURN
! *** subroutine to implement Kulmala, Laaksonen, Pirjola
    END FUNCTION getaf
!     Parameterization for sulfuric acid/water
!     nucleation rates, J. Geophys. Research (103), pp 8301-8307,
!     April 20, 1998.

!ia rev01 27.04.99 changes made to calculation of MDOT see RBiV p.2f
!ia rev02 27.04.99 security check on MDOT > SO4RAT

!ia      subroutine klpnuc( Temp, RH, H2SO4,NDOT, MDOT, M2DOT)
! getaf                                                     
    SUBROUTINE klpnuc(temp,rh,h2so4,ndot1,mdot1,so4rat)
!     IMPLICIT NONE

! *** Input:

! ambient temperature [ K ]                            
      REAL temp
! fractional relative humidity                         
      REAL rh
! sulfuric acid concentration [ ug / m**3 ]            
      REAL h2so4

      REAL so4rat
! *** Output:

!sulfuric acid production rate [ ug / ( m**3 s )]     
! particle number production rate [ # / ( m**3 s )]   
      REAL ndot1
! particle mass production rate [ ug / ( m**3 s )]    
      REAL mdot1
                 ! [ m**2 / ( m**3 s )]
      REAL m2dot

! *** Internal:

! *** NOTE, all units are cgs internally.
! particle second moment production rate               

      REAL ra
! fractional relative acidity                           
! sulfuric acid vaper concentration [ cm ** -3 ]        
      REAL nav
! water vapor concentration   [ cm ** -3 ]              
      REAL nwv
! equilibrium sulfuric acid vapor conc. [ cm ** -3 ]    
      REAL nav0
                ! to produce a nucleation rate of 1 [ cm ** -3  s ** -1
      REAL nac
! critical sulfuric acid vapor concentration [ cm ** -3 
! mole fractio of the critical nucleus                  
      REAL xal
      REAL nsulf, & ! see usage                                    
        delta
      REAL*8 & ! factor to calculate Jnuc                             
        chi
      REAL*8 & 
        jnuc
! nucleation rate [ cm ** -3  s ** -1 ]               
      REAL tt, & ! dummy variables for statement functions              
        rr
      REAL pi
      PARAMETER (pi=3.14159265)

      REAL pid6
      PARAMETER (pid6=pi/6.0)

! avogadro's constant [ 1/mol ]                   
      REAL avo
      PARAMETER (avo=6.0221367E23)

! universal gas constant [ j/mol-k ]         
!      REAL rgasuniv
!      PARAMETER (rgasuniv=8.314510)

! 1 atmosphere in pascals                               
      REAL, PARAMETER :: atm=1013.25E+02

! formula weight for h2so4 [ g mole **-1 ]
      REAL, PARAMETER :: mwh2so4=98.07948

! diameter of a 3.5 nm particle in cm                  
      REAL, PARAMETER :: d35=3.5E-07

      REAL d35sq
      PARAMETER (d35sq=d35*d35)
! volume of a 3.5 nm particle in cm**3                 
      REAL v35
      PARAMETER (v35=pid6*d35*d35sq)
!ia rev01

      REAL mp
! ***  conversion factors:
! mass of sulfate in a 3.5 nm particle               
                     ! number per cubic cm.
      REAL ugm3_ncm3
! micrograms per cubic meter to                    
      PARAMETER (ugm3_ncm3=(avo/mwh2so4)*1.0E-12)
!ia rev01
! molecules to micrograms                          
      REAL nc_ug
      PARAMETER (nc_ug=(1.0E6)*mwh2so4/avo)

! *** statement functions **************

      REAL pdens, & 
        rho_p
! particle density [ g / cm**3]                 
      REAL ad0, ad1, ad2, & 
        ad3
! coefficients for density expression    
      PARAMETER (ad0=1.738984,ad1=-1.882301,ad2=2.951849,ad3=-1.810427) 
! *** Nair and Vohra, Growth of aqueous sulphuric acid droplets
!     as a function of relative humidity,
!     J. Aerosol Science, 6, pp 265-271, 1975.

!ia rev01

! fit to Nair & Vohra data                  
                ! the mass of sulfate in a 3.5 nm particle
      REAL mp35
! arithmetic statement function to compute              
      REAL a0, a1, a2, & ! coefficients for cubic in mp35                 
        a3
      PARAMETER (a0=1.961385E2,a1=-5.564447E2,a2=8.828801E2,a3=-5.231409E2)

      REAL ph2so4, &                         ! for h2so4 and h2o vapor pressures [ Pa ]
        ph2o

! arithmetic statement functions                
      pdens(rr) = ad0 + rr*(ad1+rr*(ad2+rr*ad3))

      ph2o(tt) = exp(77.34491296-7235.4246512/tt-8.2*log(tt)+tt*5.7113E-03)

      ph2so4(tt) = exp(27.78492066-10156.0/tt)

! *** both ph2o and ph2so4 are  as in Kulmala et al.  paper

!ia rev01

! *** function for the mass of sulfate in   a 3.5 nm sphere
! *** obtained from a fit to the number of sulfate monomers in
!     a 3.5 nm particle. Uses data from Nair & Vohra
      mp35(rr) = nc_ug*(a0+rr*(a1+rr*(a2+rr*a3)))

! *** begin code:

!     The 1.0e-6 factor in the following converts from MKS to cgs units

! *** get water vapor concentration [ molecles / cm **3 ]

      nwv = rh*ph2o(temp)/(rgasuniv*temp)*avo*1.0E-6

! *** calculate the equilibrium h2so4 vapor concentration.

! *** use Kulmala corrections:

! ***
      nav0 = ph2so4(temp)/(rgasuniv*temp)*avo*1.0E-6

! *** convert sulfuric acid vapor concentration from micrograms
!     per cubic meter to molecules per cubic centimeter.

      nav = ugm3_ncm3*h2so4

! *** calculate critical concentration of sulfuric acid vapor

      nac = exp(-14.5125+0.1335*temp-10.5462*rh+1958.4*rh/temp)

! *** calculate relative acidity

      ra = nav/nav0

! *** calculate temperature correction

      delta = 1.0 + (temp-273.15)/273.14

! *** calculate molar fraction

      xal = 1.2233 - 0.0154*ra/(ra+rh) + 0.0102*log(nav) - 0.0415*log(nwv) + &
        0.0016*temp

! *** calculate Nsulf
      nsulf = log(nav/nac)

! *** calculate particle produtcion rate [ # / cm**3 ]

      chi = 25.1289*nsulf - 4890.8*nsulf/temp - 1743.3/temp - &
        2.2479*delta*nsulf*rh + 7643.4*xal/temp - 1.9712*xal*delta/rh

      jnuc = exp(chi) 
! [ # / cm**3 ]                                   
      ndot1 = (1.0E06)*jnuc

! *** calculate particle density

      rho_p = pdens(rh)

! *** get the mass of sulfate in a 3.5 nm particle

      mp = mp35(rh)                      ! in a 3.5 nm particle at ambient RH

! *** calculate mass production rate [ ug / m**3]
!     assume that the particles are 3.5 nm in diameter.

!     MDOT1 =  (1.0E12) * rho_p * v35 * Jnuc

!ia rev01

! number of micrograms of sulfate                  
      mdot1 = mp*ndot1

!ia rev02

      IF (mdot1>so4rat) THEN

        mdot1 = & 
          so4rat
! limit nucleated mass by available ma
        ndot1 = mdot1/ & 
          mp
! adjust DNDT to this                 
      END IF

      IF (mdot1==0.) ndot1 = 0.

! *** calculate M2 production rate [ m**2 / (m**3 s)]

      m2dot = 1.0E-04*d35sq*ndot1

      RETURN

END SUBROUTINE klpnuc
!------------------------------------------------------------------------------

 SUBROUTINE modpar( cblk,blkta,blkprs,     &
                    pmassn,pmassa,pmassc,  &
                    pdensn,pdensa,pdensc,  &
                    xlm,amu,               &
                    dgnuc,dgacc,dgcor,     &
                    knnuc,knacc,kncor,     &
                    diagn_opt,             & 
                    do_isorropia           )

!**    DESCRIPTION:
!       Calculates modal parameters and derived variables,
!       log-squared of std deviation, mode mean size, Knudsen number)
!       based on current values of moments for the modes.
! FSB   Now calculates the 3rd moment, mass, and density in all 3 modes.
!**
!**    Revision history:
!       Adapted 3/95 by US and CJC from EAM2's MODPAR and INIT3
!       Revised  7/23/96 by FSB to use COMMON blocks and small blocks
!        instead of large 3-d arrays, and to assume a fixed std.
!       Revised 12/06/96 by FSB to include coarse mode
!       Revised 1/10/97 by FSB to have arrays passed in call vector
!**********************************************************************
! dimension of arrays
!      INTEGER blksize
! actual number of cells in arrays
!      INTEGER numcells
!      INTEGER nspcsda
LOGICAL :: diagn_opt

! nmber of species in CBLK        
      REAL cblk(blksize,nspcsda) ! main array of variables          
      REAL blkta(blksize) ! Air temperature [ K ]            
      REAL blkprs(blksize) 
! *** output:

! Air pressure in [ Pa ]           
! concentration lower limit [ ug/m*]
! lowest particle diameter ( m )   
      REAL, PARAMETER :: dgmin=1.0E-09

! lowest particle density ( Kg/m**3 )
      REAL, PARAMETER :: densmin=1.0E03

      REAL pmassn(blksize) ! mass concentration in nuclei mode 
      REAL pmassa(blksize) ! mass concentration in accumulation
      REAL pmassc(blksize) ! mass concentration in coarse mode 
      REAL pdensn(blksize) ! average particel density in Aitken
      REAL pdensa(blksize) ! average particel density in accumu
      REAL pdensc(blksize) ! average particel density in coarse
      REAL xlm(blksize)    ! atmospheric mean free path [ m]
      REAL amu(blksize)    ! atmospheric dynamic viscosity [kg]
      REAL dgnuc(blksize)  ! Aitken mode mean diameter [ m ]
      REAL dgacc(blksize)  ! accumulation
      REAL dgcor(blksize)  ! coarse mode
      REAL knnuc(blksize)  ! Aitken mode Knudsen number
      REAL knacc(blksize)  ! accumulation
      REAL kncor(blksize) 

! coarse mode
      INTEGER lcell
   
      LOGICAL do_isorropia

      INTEGER, SAVE :: icall

! *** set up  aerosol  3rd moment, mass, density

! loop counter                            
      DO lcell = 1, numcells

! *** Aitken-mode (16 species)
         cblk(lcell,vnu3)= &
          so4fac*cblk(lcell,vso4ai)    + nh4fac*cblk(lcell,vnh4ai) + &
          h2ofac*cblk(lcell,vh2oai)    + no3fac*cblk(lcell,vno3ai) + &
          nafac*cblk(lcell,vnaai)+  clfac*cblk(lcell,vclai)+   &
          orgfac*( cblk(lcell,vpoa0i)  + cblk(lcell,vpoa1i)  + cblk(lcell,vpoa2i) + cblk(lcell,vpoa3i) ) + &
          orgfac*( cblk(lcell,vasoa0i) + cblk(lcell,vasoa1i) + cblk(lcell,vasoa2i) + cblk(lcell,vasoa3i) ) + &
          orgfac*( cblk(lcell,vbsoa1i) + cblk(lcell,vbsoa2i) + cblk(lcell,vbsoa3i) ) + &
          anthfac*cblk(lcell,vp25ai)   + anthfac*cblk(lcell,veci)                 ! density for P25 and EC are assumed the same
!     if ( do_isorropia ) then
!         cblk(lcell,vnu3) = cblk(lcell,vnu3) + & 
!          mgfac*cblk(lcell,vmgai) + cafac*cblk(lcell,vcaai) + kfac*cblk(lcell,vkai)
!     endif 

! *** Accumulation-mode
        cblk(lcell,vac3) = &
          so4fac*cblk(lcell,vso4aj)    + nh4fac*cblk(lcell,vnh4aj) + &
          h2ofac*cblk(lcell,vh2oaj)    + no3fac*cblk(lcell,vno3aj) + &
          nafac*cblk(lcell,vnaaj)+  clfac*cblk(lcell,vclaj)+   &
          orgfac*( cblk(lcell,vpoa0j)  + cblk(lcell,vpoa1j)  + cblk(lcell,vpoa2j) + cblk(lcell,vpoa3j) ) + &
          orgfac*( cblk(lcell,vasoa0j) + cblk(lcell,vasoa1j) + cblk(lcell,vasoa2j) + cblk(lcell,vasoa3j) ) + &
          orgfac*( cblk(lcell,vbsoa1j) + cblk(lcell,vbsoa2j) + cblk(lcell,vbsoa3j) ) + &
          polsfac*( cblk(lcell,vpols) ) + &
          anthfac*cblk(lcell,vp25aj)   + anthfac*cblk(lcell,vecj)                 ! density for P25 and EC are assumed the same
!     if ( do_isorropia ) then
!        cblk(lcell,vac3) = cblk(lcell,vac3) + &
!          mgfac*cblk(lcell,vmgaj) + cafac*cblk(lcell,vcaaj) + kfac*cblk(lcell,vkaj)
!     endif

! RAR: Debugging:
IF (diagn_opt .AND. icall<100) then
    WRITE(6,*) 'modpar: so4fac,nh4fac,h2ofac,no3fac,orgfac,anthfac ', so4fac, nh4fac, h2ofac, no3fac, orgfac, anthfac
    WRITE(6,*) 'modpar: cblk(lcell,vecj),cblk(lcell,veci) ', cblk(lcell,vecj),cblk(lcell,veci)
    WRITE(6,*) 'modpar: cblk(lcell,vnu3),cblk(lcell,vac3) ', cblk(lcell,vnu3),cblk(lcell,vac3)
ENDIF
!    WRITE(6,*) 'modpar: cblk(lcell,vpols),cblk(lcell,vpolp) ', cblk(lcell,vpols),cblk(lcell,vpolp)

! *** coarse mode
        cblk(lcell,vcor3) = max(conmin,(soilfac*cblk(lcell, & ! ghan rely on conmin applied to mass, not moment
          vsoila)+seasfac*cblk(lcell,vseas)+anthfac*cblk(lcell,vantha)+polpfac*cblk(lcell,vpolp)))
!        cblk(lcell,vcor3) = soilfac*cblk(lcell, &
!          vsoila)+seasfac*cblk(lcell,vseas)+anthfac*cblk(lcell,vantha)

! *** now get particle mass and density
! *** Aitken-mode:
        pmassn(lcell) = max( conmin, &                                                              ! conmin = 1.E-16
         ( cblk(lcell,vso4ai)  + cblk(lcell,vnh4ai)  + cblk(lcell,vh2oai)  + cblk(lcell,vno3ai) +  &
           cblk(lcell,vpoa0i)  + cblk(lcell,vpoa1i)  + cblk(lcell,vpoa2i)  + cblk(lcell,vpoa3i) + &
           cblk(lcell,vasoa0i) + cblk(lcell,vasoa1i) + cblk(lcell,vasoa2i) + cblk(lcell,vasoa3i) + &
           cblk(lcell,vbsoa1i) + cblk(lcell,vbsoa2i) + cblk(lcell,vbsoa3i) + &
           cblk(lcell,vp25ai)  + cblk(lcell,veci)    + cblk(lcell,vclai)   + cblk(lcell,vnaai) ) )
!    if ( do_isorropia ) then
!       pmassn(lcell) = max( conmin, ( pmassn(lcell) +  & 
!                                        cblk(lcell,vmgai) + cblk(lcell,vcaai) + cblk(lcell,vkai) ) )
!    endif
       

! *** Accumulation-mode:
        pmassa(lcell) = max( conmin, &
         ( cblk(lcell,vso4aj)  + cblk(lcell,vnh4aj)  + cblk(lcell,vh2oaj)  + cblk(lcell,vno3aj) +  &
           cblk(lcell,vpoa0j)  + cblk(lcell,vpoa1j)  + cblk(lcell,vpoa2j)  + cblk(lcell,vpoa3j) + &
           cblk(lcell,vasoa0j) + cblk(lcell,vasoa1j) + cblk(lcell,vasoa2j) + cblk(lcell,vasoa3j) + &
           cblk(lcell,vbsoa1j) + cblk(lcell,vbsoa2j) + cblk(lcell,vbsoa3j) + &
           cblk(lcell,vp25aj)  + cblk(lcell,vecj)    + cblk(lcell,vclaj)   + cblk(lcell,vnaaj) ))
!           cblk(lcell,vpols) ) ) ! Why isn't this added here?
!    if ( do_isorropia ) then
!       pmassa(lcell) = max( conmin, ( pmassa(lcell) + & 
!                                        cblk(lcell,vmgaj) + cblk(lcell,vcaaj) + cblk(lcell,vkaj) ) )
!    endif

! RAR: Debugging:
IF (diagn_opt .AND. icall<100) then
    WRITE(6,*) 'modpar: pmassn(lcell), pmassa(lcell) ', pmassn(lcell), pmassa(lcell)
    WRITE(6,*) 'modpar: cblk(lcell,vecj), cblk(lcell,veci) ', cblk(lcell,vecj), cblk(lcell,veci)
ENDIF
!    WRITE(6,*) 'modpar: cblk(lcell,vpols), cblk(lcell,vpolp) ', cblk(lcell,vpols),cblk(lcell,vpolp)
! *** coarse mode:
        pmassc(lcell) = max(conmin,cblk(lcell,vsoila)+cblk(lcell,vseas)+cblk( &
          lcell,vantha) + cblk(lcell,vpolp))
!        pmassc(lcell)= conmin       ! RAR: just a placeholder

      END DO
! *** now get particle density, mean free path, and dynamic viscosity

! aerosol  3rd moment and  mass
      DO lcell = 1,numcells

! Density and mean free path              ! *** density in [ kg m**-3 ]
        pdensn(lcell) = max(densmin,(f6dpim9*pmassn(lcell)/cblk(lcell,vnu3)))
        pdensa(lcell) = max(densmin,(f6dpim9*pmassa(lcell)/cblk(lcell,vac3)))
        pdensc(lcell) = max(densmin,(f6dpim9*pmassc(lcell)/cblk(lcell,vcor3)))  ! RAR:

! *** Calculate mean free path [ m ]:
        xlm(lcell) = 6.6328E-8*pss0*blkta(lcell)/(tss0*blkprs(lcell))

! *** 6.6328E-8 is the sea level values given in Table I.2.8
! *** on page 10 of U.S. Standard Atmosphere 1962

! *** 	Calculate dynamic viscosity [ kg m**-1 s**-1 ]:

! *** U.S. Standard Atmosphere 1962 page 14 expression
!     for dynamic viscosity is:
!     dynamic viscosity =  beta * T * sqrt(T) / ( T + S)
!     where beta = 1.458e-6 [ kg sec^-1 K**-0.5 ], s = 110.4 [ K ].

      amu(lcell) = 1.458E-6*blkta(lcell)*sqrt(blkta(lcell))/ &
          (blkta(lcell)+110.4)

      END DO

! RAR: Debugging:
IF (diagn_opt .AND. icall<100) then
    WRITE(6,*) 'modpar: blksize,nspcsda,xlm,amu ', blksize,nspcsda,xlm,amu
ENDIF

!...............   Standard deviation fixed in both modes, so
!...............   diagnose diameter from 3rd moment and number concentration

!  density and mean free path 
      DO lcell = 1,numcells

! calculate diameters             
        dgnuc(lcell) = max(dgmin,(cblk(lcell,vnu3)/(cblk(lcell,vnu0)*esn36))**one3)

        dgacc(lcell) = max(dgmin,(cblk(lcell,vac3)/(cblk(lcell,vac0)*esa36))**one3)

!        dgcor(lcell) = dgmin    ! RAR:
        dgcor(lcell) = max(dgmin,(cblk(lcell,vcor3)/(cblk(lcell,vcorn)*esc36))**one3)

! RAR: Debugging:
IF (diagn_opt .AND. icall<200) then
    WRITE(6,*) 'modpar: pdensn,pdensa,pdensc ', pdensn,pdensa,pdensc
    WRITE(6,*) 'modpar: dgnuc,dgacc,dgcor ', dgnuc,dgacc,dgcor
!    WRITE(6,*) 'modpar: cw_phase ', cw_phase
    icall= icall+1
ENDIF

! when running with cloudborne aerosol, apply some very mild bounding
! to avoid unrealistic dg values
!      if (cw_phase > 0) then
!        dgnuc(lcell) = max( dgnuc(lcell), dginin*0.2  )  !  > 0.002 um
!        dgnuc(lcell) = min( dgnuc(lcell), dginin*10.0 )  !  < 0.10  um
!        dgacc(lcell) = max( dgacc(lcell), dginia*0.2  )  !  > 0.014 um
!        dgacc(lcell) = min( dgacc(lcell), dginia*10.0 )  !  < 0.7 um
!        dgcor(lcell) = max( dgcor(lcell), dginic*0.2  )  !  > 0.2 um
!        dgcor(lcell) = min( dgcor(lcell), dginic*10.0 )  ! < 10.0 um
!      end if

      END DO
! end loop on diameters                              
      DO lcell = 1,numcells
! Calculate Knudsen numbers           
        knnuc(lcell) = 2.0*xlm(lcell)/dgnuc(lcell)

        knacc(lcell) = 2.0*xlm(lcell)/dgacc(lcell)

        kncor(lcell) = 2.0*xlm(lcell)/dgcor(lcell)

      END DO

! end loop for  Knudsen numbers                       
      RETURN

END SUBROUTINE modpar
!------------------------------------------------------------------------------

SUBROUTINE nuclcond(cblk,dt,layer,blkta,blkprs,blkrh,so4rat, &
             orgpoa1rat,orgpoa2rat,orgpoa3rat, &
             organt0rat,organt1rat,organt2rat,organt3rat, &
             orgbio1rat,orgbio2rat,orgbio3rat, &
             drog_in,dgnuc,dgacc,    &
             fconcn,fconca,fconcn_org,fconca_org,dmdt,dndt,deltaso4a, &
             cgrn3,cgra3,brrto,diagn_opt, &
             igrid,jgrid,kgrid)

!***********************************************************************
!**    DESCRIPTION:  calculates aerosol nucleation and condensational
!**    growth rates using Binkowski and Shankar (1995) method.

! *** In this version, the method od RPM is followed where
!     the diffusivity, the average molecular velocity, and
!     the accomodation coefficient for sulfuric acid are used for
!     the organics. This is for consistency.
!       Future versions will use the correct values.  FSB 12/12/96
!**
!**    Revision history:
!       prototype 1/95 by Uma and Carlie
!       Corrected 7/95 by Uma for condensation of mass not nucleated
!       and mass conservation check
!       Revised   8/95 by US to calculate air density in stmt function
!                 and collect met variable stmt funcs in one include fil
!       Revised 7/25/96 by FSB to use block structure.
!       Revised 9/17/96 by FSB to use Y&K or K&W Nucleation mechanism
!       Revised 11/15/96 by FSB to use MKS,  and mom m^-3 units.
!       Revised 1/13/97 by FSB to pass arrays and simplify code.
!       Added   23/03/99 by BS growth factors for organics
!**********************************************************************
!     IMPLICIT NONE

!     Includes:
! *** arguments

! *** input;
!USE module_configure, only: grid_config_rec_type
!TYPE (grid_config_rec_type), INTENT (in) :: config_flags

LOGICAL :: diagn_opt

! dimension of arrays
!      INTEGER blksize
      INTEGER layer
! number of species in CBLK       
!      INTEGER nspcsda
! actual number of cells in arrays
!      INTEGER numcells
! # of organic aerosol precursor
      REAL cblk(blksize,nspcsda) ! main array of variables         
! model time step in  SECONDS     
      REAL dt
      REAL blkta(blksize) ! Air temperature [ K ]           
      REAL blkprs(blksize) ! Air pressure in [ Pa ]          
      REAL blkrh(blksize) ! Fractional relative humidity    
      REAL so4rat(blksize) ! rate [  ug/m**3 /s ]
      REAL brrto
!bs
! sulfate gas-phase production    
! total # of cond. vapors & SOA spe
!      INTEGER ncv
!      INTEGER nacv
!bs * anthropogenic organic condensable vapor production rate
! # of anthrop. cond. vapors & SOA 
      REAL drog_in(ldrog_soa2) !bs
! Delta ROG conc. [ppm]             

! POA growth rates
      REAL orgpoa1rat(blksize), orgpoa2rat(blksize), orgpoa3rat(blksize)

! anthropogenic vapor production rates
      REAL organt0rat(blksize), organt1rat(blksize), organt2rat(blksize), organt3rat(blksize)

! biogenic vapor production rates
      REAL orgbio1rat(blksize), orgbio2rat(blksize), orgbio3rat(blksize)

      REAL dgnuc(blksize)
      REAL dgacc(blksize) 

! coarse mode                           
      REAL fconcn(blksize)                                 ! Aitken mode  [ 1 / s ]
! reciprocal condensation rate          
      REAL fconca(blksize)                                 ! acclumulation mode [ 1 / s ]
! reciprocal condensation rate          
      REAL fconcn_org(blksize)                                 ! Aitken mode  [ 1 / s ]
! reciprocal condensation rate          
      REAL fconca_org(blksize)                                 ! acclumulation mode [ 1 / s ]
! reciprocal condensation rate          
      REAL dmdt(blksize)                                 ! by particle formation [ ug/m**3 /s ]
! rate of production of new mass concent
      REAL dndt(blksize)                                 ! concentration by particle formation [#
! rate of producton of new particle numb
      REAL deltaso4a(blksize)                                 ! sulfate aerosol by condensation [ ug/m
! increment of concentration added to   
      REAL cgrn3(blksize)                                 ! Aitken mode [ 3rd mom/m **3 s ]
! growth rate for 3rd moment for        
      REAL cgra3(blksize)                                 ! Accumulation mode   

!...........    SCRATCH local variables and their descriptions:

! growth rate for 3rd moment for        

      INTEGER lcell

! LOOP INDEX                                     
! conv rate so2 --> so4 [mom-3/g/s]     
      REAL chemrat
! conv rate for organics [mom-3/g/s]    
      REAL chemrat_org
      REAL am1n, & ! 1st mom density (nuc, acc modes) [mom_
        am1a
      REAL am2n, & ! 2nd mom density (nuc, acc modes) [mom_
        am2a
      REAL gnc3n, & ! near-cont fns (nuc, acc) for mom-3 den
        gnc3a
      REAL gfm3n, & ! free-mol  fns (nuc, acc) for mom-3 den
        gfm3a
! total reciprocal condensation rate    
      REAL fconc

      REAL td
! d * tinf (cgs)                        
      REAL*8 & ! Cnstant to force 64 bit evaluation of 
        one88
      PARAMETER (one88=1.0D0)
!  *** variables to set up sulfate and organic condensation rates

! sulfuric acid vapor at current time step            
      REAL vapor1
!                                    chemistry and emissions
      REAL vapor2
! Sulfuric acid vapor prior to addition from          
!bs
      REAL deltavap
!bs * start update
!bs
! change to vapor at previous time step 
      REAL diffcorr

!bs *
      REAL csqt_org
!bs * end update
!bs
      REAL csqt

   INTEGER igrid,jgrid,kgrid

   INTEGER, SAVE :: icall
!.......................................................................
!   begin body of subroutine  NUCLCOND


!...........   Main computational grid-traversal loop nest
!...........   for computing condensation and nucleation:

      DO lcell = 1,numcells
! *** First moment:

!  1st loop over NUMCELLS                  
        am1n = cblk(lcell,vnu0)*dgnuc(lcell)*esn04
        am1a = cblk(lcell,vac0)*dgacc(lcell)*esa04

!..............   near-continuum factors [ 1 / sec ]
!bs
!bs * adopted from code of FSB
!bs * correction to DIFFSULF and DIFFORG for temperature and pressure
!bs
        diffcorr = (pss0/blkprs(lcell))*(blkta(lcell)/273.16)**1.
!bs
        gnc3n = cconc*am1n*diffcorr
        gnc3a = cconc*am1a*diffcorr

! *** Second moment:

        am2n = cblk(lcell,vnu0)*dgnuc(lcell)*dgnuc(lcell)*esn16
        am2a = cblk(lcell,vac0)*dgacc(lcell)*dgacc(lcell)*esa16

        csqt = ccofm*sqrt(blkta(lcell))
!...............   free molecular factors [ 1 / sec ]

! put in temperature fac
        gfm3n = csqt*am2n
        gfm3a = csqt*am2a

! *** Condensation factors in [ s**-1] for h2so4
! *** In the future, separate factors for condensing organics will
!      be included. In this version, the h2so4 values are used.

!...............   Twice the harmonic mean of fm, nc functions:
! *** Force 64 bit evaluation:

        fconcn(lcell) = one88*gnc3n*gfm3n/(gnc3n+gfm3n)         ! one88=1.0
        fconca(lcell) = one88*gnc3a*gfm3a/(gnc3a+gfm3a)
        fconc = fconcn(lcell) + fconca(lcell)

! RAR: Debugging:
IF (diagn_opt .AND. icall<200 .AND. kgrid<3) then
    WRITE(6,*) 'nuclcond: gnc3n,gfm3n ',gnc3n,gfm3n
    WRITE(6,*) 'nuclcond: fconcn,fconca ',fconcn,fconca
end if

! *** NOTE: FCONCN and FCONCA will be redefined below <<<<<<
!bs * start modifications for organcis
        gnc3n = cconc_org*am1n*diffcorr         ! cconc_org - molecular diffusivity for organics
        gnc3a = cconc_org*am1a*diffcorr

        csqt_org = ccofm_org*sqrt(blkta(lcell))
        gfm3n = csqt_org*am2n
        gfm3a = csqt_org*am2a

        fconcn_org(lcell) = one88*gnc3n*gfm3n/(gnc3n+gfm3n)
        fconca_org(lcell) = one88*gnc3a*gfm3a/(gnc3a+gfm3a)

! RAR: Debugging:
IF (diagn_opt .AND. icall<200 .AND. kgrid<3) then
    WRITE(6,*) 'nuclcond: cconc_org,ccofm_org ',cconc_org,ccofm_org
    WRITE(6,*) 'nuclcond: gnc3n,gfm3n,gnc3a,gfm3a ',gnc3n,gfm3n,gnc3a,gfm3a
    WRITE(6,*) 'nuclcond: fconcn_org,fconca_org ',fconcn_org,fconca_org
end if

!bs * end modifications for organics
! *** calculate the total change to sulfuric acid vapor from production and condensation

        vapor1 = cblk(lcell,vsulf) ! current sulfuric acid vapor
        vapor2 = cblk(lcell,vsulf) - so4rat(lcell)* dt
! vapor at prev
        vapor2 = max(0.0,vapor2)
        deltavap = max(0.0,(so4rat(lcell)/fconc-vapor2)*(1.0-exp(-fconc*dt)))

! *** Calculate increment in total sufate aerosol mass concentration

! *** This follows the method of Youngblood & Kreidenweis.!bs
!bs        DELTASO4A( LCELL) = MAX( 0.0, SO4RAT(LCELL) * DT - DELTAVAP)
!bs
!bs * allow DELTASO4A to be negative, but the change must not be larger
!bs * than the amount of vapor available.
!bs
        deltaso4a(lcell) = max(-1.*cblk(lcell,vsulf), so4rat(lcell)*dt-deltavap)

! *** zero out growth coefficients
        cgrn3(lcell) = 0.0
        cgra3(lcell) = 0.0

      END DO

! *** Select method of nucleation
! End 1st loop over NUMCELLS
      IF (inucl==1) THEN

! *** Do Youngblood & Kreidenweis  Nucleation

!         CALL BCSUINTF(DT,SO4RAT,FCONCN,FCONCA,BLKTA,BLKRH,
!     &        DNDT,DMDT,NUMCELLS,blksize,
!     &        VAPOR1)

      ELSE IF (inucl==0) THEN

! *** Do Kerminen & Wexler Nucleation

!         CALL nuclKW(DT,SO4RAT,FCONCN,FCONCA,BLKTA,BLKRH,
!     &        DNDT,DMDT,NUMCELLS,blksize)

      ELSE IF (inucl==2) THEN   ! That is a default choise!

!bs ** Do Kulmala et al. Nucleation
!       if(dndt(1).lt.-10.)print *,'before klpnuc',blkta(1),blkrh(1),vapor1,dndt(1),dmdt(1),so4rat(1)

        if(blkta(1).ge.233.15.and.blkrh(1).ge.0.1)then
           CALL klpnuc(blkta(1),blkrh(1),vapor1,dndt(1),dmdt(1),so4rat(1))
        else
           dndt(1)=0.
           dmdt(1)=0.
        endif

!       CALL klpnuc(blkta(1),blkrh(1),vapor1,dndt(1),dmdt(1),so4rat(1))
!       if(dndt(1).lt.-10.)print *,'after klpnuc',dndt(1),dmdt(1)
        IF (dndt(1)==0.) dmdt(1) = 0.
        IF (dmdt(1)==0.) dndt(1) = 0.

        dmdt(1) = max(dmdt(1),1e-30) ! Bug fix from Kai Wang

      END IF

! RAR: Debugging:
IF (diagn_opt .AND. icall<200 .AND. kgrid<3) then
    WRITE(6,*) 'nuclcond: inucl,deltaso4a,deltavap', inucl,deltaso4a,deltavap
ENDIF

!bs
!bs * Secondary organic aerosol module (SOA_VBS)
!bs
! end of selection of nucleation method

!      CALL soa2(layer,blkta,blkprs,organt1rat,organt2rat,organt3rat, &
!        organt4rat,orgbio1rat,orgbio2rat,orgbio3rat,orgbio4rat,drog,ldrog_vbs,ncv, &
!        nacv,cblk,blksize,nspcsda,numcells,dt,igrid,jgrid,kgrid,brrto )

 CALL soa2_vbs( blkta,blkprs,&
                orgpoa1rat,orgpoa2rat,orgpoa3rat, &
                organt0rat,organt1rat,organt2rat,organt3rat, &
                orgbio1rat,orgbio2rat,orgbio3rat, &
                drog_in,cblk,dt,brrto, diagn_opt, &
                igrid,jgrid,kgrid )

      DO lcell = 1, numcells

! *** redefine FCONCN & FCONCA to be the nondimensional fractionaL
!     condensation factors

        td = 1.0/(fconcn(lcell)+fconca(lcell))
        fconcn(lcell) = td*fconcn(lcell)
        fconca(lcell) = td*fconca(lcell)
!bs
        td = 1.0/(fconcn_org(lcell)+fconca_org(lcell))
        fconcn_org(lcell) = td*fconcn_org(lcell)
        fconca_org(lcell) = td*fconca_org(lcell)
!bs
      END DO

! ***  Begin second loop over cells

      DO lcell = 1,numcells
! *** note CHEMRAT includes  species other than sulfate.

! 3rd loop on NUMCELLS                     
        chemrat = so4fac*so4rat(lcell) ! [mom3 m**-3 s-

        chemrat_org= orgfac*( &
          orgpoa1rat(lcell) + orgpoa2rat(lcell) + orgpoa3rat(lcell) + &
          organt0rat(lcell) + organt1rat(lcell) + organt2rat(lcell) + organt3rat(lcell) + &
          orgbio1rat(lcell) + orgbio2rat(lcell) + orgbio3rat(lcell) )

!        chemrat_org = orgfac*(organt1rat(lcell)+organt2rat(lcell)+organt3rat( &
!          lcell)+organt4rat(lcell)+orgbio1rat(lcell)+orgbio2rat(lcell)+ &
!          orgbio3rat(lcell)+orgbio4rat(lcell))

! *** Calculate the production rates for new particle
! [mom3 m**-3 s-
        cgrn3(lcell) = so4fac*dmdt(lcell) 
! Rate of increase of 3rd
        chemrat = chemrat - cgrn3(lcell)                                            !bs 3rd moment production fro

!bs Remove the rate of new pa
        chemrat = max(chemrat,0.0) 
! *** Now calculate the rate of condensation on existing particles.

! Prevent CHEMRAT from being negativ
        cgrn3(lcell) = cgrn3(lcell) + chemrat*fconcn(lcell) + chemrat_org*fconcn_org(lcell)
        cgra3(lcell) =                chemrat*fconca(lcell) + chemrat_org*fconca_org(lcell)

! RAR: Debugging:
IF (diagn_opt .AND. icall<200) then
    WRITE(6,*) 'nuclcond: cgrn3,cgra3,chemrat,chemrat_org',cgrn3,cgra3,chemrat,chemrat_org
    icall=icall+1
end if

      END DO
!  end 2nd loop over NUMCELLS           
      RETURN

    END SUBROUTINE nuclcond
!------------------------------------------------------------------------------

! nuclcond                              
REAL FUNCTION poly4(a,x)
  REAL a(4), x

  poly4 = a(1) + x*(a(2)+x*(a(3)+x*(a(4))))
  RETURN
END FUNCTION poly4
REAL FUNCTION poly6(a,x)
  REAL a(6), x

  poly6 = a(1) + x*(a(2)+x*(a(3)+x*(a(4)+x*(a(5)+x*(a(6))))))
  RETURN
END FUNCTION poly6
!-----------------------------------------------------------------------

SUBROUTINE rpmares_old(so4,hno3,no3,nh3,nh4,rh,temp,aso4,ano3,ah2o,anh4,gnh3,gno3,diagn_opt)
! Description:

!   ARES calculates the chemical composition of a sulfate/nitrate/
!   ammonium/water aerosol based on equilibrium thermodynamics.

!   This code considers two regimes depending upon the molar ratio
!   of ammonium to sulfate.

!   For values of this ratio less than 2,the code solves a cubic for
!   hydrogen ion molality, HPLUS,  and if enough ammonium and liquid
!   water are present calculates the dissolved nitric acid. For molal
!   ionic strengths greater than 50, nitrate is assumed not to be present

!   For values of the molar ratio of 2 or greater, all sulfate is assumed
!   to be ammonium sulfate and a calculation is made for the presence of
!   ammonium nitrate.

!   The Pitzer multicomponent approach is used in subroutine ACTCOF to
!   obtain the activity coefficients. Abandoned -7/30/97 FSB

!   The Bromley method of calculating the activity coefficients is used in this version

!   The calculation of liquid water is done in subroutine water. Details for both calculations are given
!   in the respective subroutines.

!   Based upon MARS due to
!   P. Saxena, A.B. Hudischewskyj, C. Seigneur, and J.H. Seinfeld,
!   Atmos. Environ., vol. 20, Number 7, Pages 1471-1483, 1986.

!   and SCAPE due to
!   Kim, Seinfeld, and Saxeena, Aerosol Ceience and Technology,
!   Vol 19, number 2, pages 157-181 and pages 182-198, 1993.

! NOTE: All concentrations supplied to this subroutine are TOTAL
!       over gas and aerosol phases

! Parameters:

!  SO4   : Total sulfate in MICROGRAMS/M**3 as sulfate (IN)
!  HNO3  : Nitric Acid in MICROGRAMS/M**3 as nitric acid (IN)
!  NO3   : Total nitrate in MICROGRAMS/M**3 as nitric acid (IN)
!  NH3   : Total ammonia in MICROGRAMS/M**3 as ammonia (IN)
!  NH4   : Ammonium in MICROGRAMS/M**3 as ammonium (IN)
!  RH    : Fractional relative humidity (IN)
!  TEMP  : Temperature in Kelvin (IN)
!  GNO3  : Gas phase nitric acid in MICROGRAMS/M**3 (OUT)
!  GNH3  : Gas phase ammonia in MICROGRAMS/M**3 (OUT)
!  ASO4  : Aerosol phase sulfate in MICROGRAMS/M**3 (OUT)
!  ANO3  : Aerosol phase nitrate in MICROGRAMS/M**3 (OUT)
!  ANH4  : Aerosol phase ammonium in MICROGRAMS/M**3 (OUT)
!  AH2O  : Aerosol phase water in MICROGRAMS/M**3 (OUT)
!  NITR  : Number of iterations for obtaining activity coefficients  (OU
!  NR    : Number of real roots to the cubic in the low ammonia case (OU

! Revision History:
!      Who       When        Detailed description of changes
!   ---------   --------  -------------------------------------------
!   S.Roselle   11/10/87  Received the first version of the MARS code
!   S.Roselle   12/30/87  Restructured code
!   S.Roselle   2/12/88   Made correction to compute liquid-phase
!                         concentration of H2O2.
!   S.Roselle   5/26/88   Made correction as advised by SAI, for
!                         computing H+ concentration.
!   S.Roselle   3/1/89    Modified to operate with EM2
!   S.Roselle   5/19/89   Changed the maximum ionic strength from
!                         100 to 20, for numerical stability.
!   F.Binkowski 3/3/91    Incorporate new method for ammonia rich case
!                         using equations for nitrate budget.
!   F.Binkowski 6/18/91   New ammonia poor case which
!                         omits letovicite.
!   F.Binkowski 7/25/91   Rearranged entire code, restructured
!                         ammonia poor case.
!   F.Binkowski 9/9/91    Reconciled all cases of ASO4 to be output
!                         as SO4--
!   F.Binkowski 12/6/91   Changed the ammonia defficient case so that
!                         there is only neutralized sulfate (ammonium
!                         sulfate) and sulfuric acid.
!   F.Binkowski 3/5/92    Set RH bound on AWAS to 37 % to be in agreemen
!                          with the Cohen et al. (1987)  maximum molalit
!                          of 36.2 in Table III.( J. Phys Chem (91) page
!                          4569, and Table IV p 4587.)
!   F.Binkowski 3/9/92    Redid logic for ammonia defficient case to rem
!                         possibility for denomenator becoming zero;
!                         this involved solving for HPLUS first.
!                         Note that for a relative humidity
!                          less than 50%, the model assumes that there i
!                          aerosol nitrate.
!   F.Binkowski 4/17/95   Code renamed  ARES (AeRosol Equilibrium System
!                          Redid logic as follows
!                         1. Water algorithm now follows Spann & Richard
!                         2. Pitzer Multicomponent method used
!                         3. Multicomponent practical osmotic coefficien
!                            use to close iterations.
!                         4. The model now assumes that for a water
!                            mass fraction WFRAC less than 50% there is
!                            no aerosol nitrate.
!   F.Binkowski 7/20/95   Changed how nitrate is calculated in ammonia p
!                         case, and changed the WFRAC criterion to 40%.
!                         For ammonium to sulfate ratio less than 1.0
!                         all ammonium is aerosol and no nitrate aerosol
!                         exists.
!   F.Binkowski 7/21/95   Changed ammonia-ammonium in ammonia poor case
!                         allow gas-phase ammonia to exist.
!   F.Binkowski 7/26/95   Changed equilibrium constants to values from
!                         Kim et al. (1993)
!   F.Binkowski 6/27/96   Changed to new water format
!   F.Binkowski 7/30/97   Changed to Bromley method for multicomponent
!                         activity coefficients. The binary activity coe
!                         are the same as the previous version
!   F.Binkowski 8/1/97    Chenged minimum sulfate from 0.0 to 1.0e-6 i.e
!                         1 picogram per cubic meter

! molecular weight for NH3
!      REAL mwnh3
!      PARAMETER (mwnh3=17.03061)

! molecular weight for NH4           
!      REAL mwnh4
!      PARAMETER (mwnh4=18.03858)

! molecular weight for Chloride
!      REAL mwcl
!      PARAMETER (mwcl=35.453)

! molecular weight for AIR           
!      REAL mwair
!      PARAMETER (mwair=28.964)

! molecular weight for Organic Species
!      REAL mworg
!      PARAMETER (mworg=16.0)

LOGICAL :: diagn_opt

! molecular weight for Letovicite    
      REAL mwlct
      PARAMETER (mwlct=3.0*mwnh4+2.0*mwso4+1.0080)

! molecular weight for Ammonium Sulfa
      REAL mwas
      PARAMETER (mwas=2.0*mwnh4+mwso4)

! molecular weight for Ammonium Bisul
      REAL mwabs
      PARAMETER (mwabs=mwnh4+mwso4+1.0080)

! molecular weight for NaCl          
      REAL mwnacl
      PARAMETER (mwnacl=58.44277)

! molecular weight for Chloride      
      REAL mwcl
      PARAMETER (mwcl=35.453)


!...........ARGUMENTS and their descriptions

!iamodels3
      REAL so4
! Total sulfate in micrograms / m**3 
! Total nitric acid in micrograms / m
      REAL hno3
! Total nitrate in micrograms / m**3 
      REAL no3
! Total ammonia in micrograms / m**3 
      REAL nh3
! Total ammonium in micrograms / m**3
      REAL nh4
! Fractional relative humidity       
      REAL rh
! Temperature in Kelvin              
      REAL temp
! Aerosol sulfate in micrograms / m**
      REAL aso4
! Aerosol nitrate in micrograms / m**
      REAL ano3
! Aerosol liquid water content water 
      REAL ah2o
! Aerosol ammonium in micrograms / m*
      REAL anh4
! Gas-phase nitric acid in micrograms
      REAL gno3
      REAL gnh3
!...........SCRATCH LOCAL VARIABLES and their descriptions:

! Gas-phase ammonia in micrograms / m
! Index set to percent relative humid
      INTEGER irh
! Number of iterations for activity c
      INTEGER nitr
! Loop index for iterations          
      INTEGER nnn
      INTEGER nr
! Number of roots to cubic equation f
      REAL*8 & ! Coefficients and roots of        
        a0
      REAL*8 & ! Coefficients and roots of        
        a1
      REAL*8 & ! Coefficients and roots of        
        a2
! Coefficients and discriminant for q
      REAL aa
! internal variables ( high ammonia c
      REAL bal
! Coefficients and discriminant for q
      REAL bb
! Variables used for ammonia solubili
      REAL bhat
! Coefficients and discriminant for q
      REAL cc
! Factor for conversion of units     
      REAL convt
! Coefficients and discriminant for q
      REAL dd
! Coefficients and discriminant for q
      REAL disc
! Relative error used for convergence
      REAL eror
!  Free ammonia concentration , that 
      REAL fnh3
! Activity Coefficient for (NH4+, HSO
      REAL gamaab
! Activity coefficient for (NH4+, NO3
      REAL gamaan
! Variables used for ammonia solubili
      REAL gamahat
! Activity coefficient for (H+ ,NO3-)
      REAL gamana
! Activity coefficient for (2H+, SO4-
      REAL gamas1
! Activity coefficient for (H+, HSO4-
      REAL gamas2
! used for convergence of iteration  
      REAL gamold
! internal variables ( high ammonia c
      REAL gasqd
! Hydrogen ion (low ammonia case) (mo
      REAL hplus
! Equilibrium constant for ammoniua t
      REAL k1a
! Equilibrium constant for sulfate-bi
      REAL k2sa
! Dissociation constant for ammonium 
      REAL k3
! Equilibrium constant for ammonium n
      REAL kan
! Variables used for ammonia solubili
      REAL khat
! Equilibrium constant for nitric aci
      REAL kna
! Henry's Law Constant for ammonia   
      REAL kph
! Equilibrium constant for water diss
      REAL kw
! Internal variable using KAN        
      REAL kw2
! Nitrate (high ammonia case) (moles 
      REAL man
! Sulfate (high ammonia case) (moles 
      REAL mas
! Bisulfate (low ammonia case) (moles
      REAL mhso4
! Nitrate (low ammonia case) (moles /
      REAL mna
! Ammonium (moles / kg water)        
      REAL mnh4
! Total number of moles of all ions  
      REAL molnu
! Sulfate (low ammonia case) (moles /
      REAL mso4
! Practical osmotic coefficient      
      REAL phibar
! Previous value of practical osmotic
      REAL phiold
! Molar ratio of ammonium to sulfate 
      REAL ratio
! Internal variable using K2SA       
      REAL rk2sa
! Internal variables using KNA       
      REAL rkna
! Internal variables using KNA       
      REAL rknwet
      REAL rr1
      REAL rr2
! Ionic strength                     
      REAL stion
! Internal variables for temperature 
      REAL t1
! Internal variables for temperature 
      REAL t2
! Internal variables of convenience (
      REAL t21
! Internal variables of convenience (
      REAL t221
! Internal variables for temperature 
      REAL t3
! Internal variables for temperature 
      REAL t4
! Internal variables for temperature 
      REAL t6
! Total ammonia and ammonium in micro
      REAL tnh4
! Total nitrate in micromoles / meter
      REAL tno3
! Tolerances for convergence test    
      REAL toler1
! Tolerances for convergence test    
      REAL toler2
! Total sulfate in micromoles / meter
      REAL tso4
! 2.0 * TSO4  (high ammonia case) (mo
      REAL twoso4
! Water mass fraction                
      REAL wfrac
                                   ! micrograms / meter **3 on output
      REAL wh2o
                                   ! internally it is 10 ** (-6) kg (wat
                                   ! the conversion factor (1000 g = 1 k
                                   ! for AH2O output
! Aerosol liquid water content (inter
! internal variables ( high ammonia c
      REAL wsqd
! Nitrate aerosol concentration in mi
      REAL xno3
! Variable used in quadratic solution
      REAL xxq
! Ammonium aerosol concentration in m
      REAL ynh4
! Water variable saved in case ionic 
      REAL zh2o

      REAL zso4
! Total sulfate molality - mso4 + mhs
      REAL cat(2) ! Array for cations (1, H+); (2, NH4+
      REAL an(3) ! Array for anions (1, SO4--); (2, NO
      REAL crutes(3) ! Coefficients and roots of          
      REAL gams(2,3) ! Array of activity coefficients     
! Minimum value of sulfate laerosol c
      REAL minso4
      PARAMETER (minso4=1.0E-6/mwso4)
      REAL floor
      PARAMETER (floor=1.0E-30)

      INTEGER, SAVE :: icall
!-----------------------------------------------------------------------
!  begin body of subroutine RPMARES

!...convert into micromoles/m**3
!iamodels3 merge NH3/NH4 , HNO3,NO3 here
! minimum concentration              
      tso4 = max(0.0,so4/mwso4)
      tno3 = max(0.0,(no3/mwno3+hno3/mwhno3))
      tnh4 = max(0.0,(nh3/mwnh3+nh4/mwnh4))

!...now set humidity index IRH as a percent

      irh = nint(100.0*rh)

!...Check for valid IRH

      irh = max(1,irh)
      irh = min(99,irh)

!...Specify the equilibrium constants at  correct
!...  temperature.  Also change units from ATM to MICROMOLE/M**3 (for KA
!...  KPH, and K3 )
!...  Values from Kim et al. (1993) except as noted.

      convt = 1.0/(0.082*temp)
      t6 = 0.082E-9*temp
      t1 = 298.0/temp
      t2 = alog(t1)
      t3 = t1 - 1.0
      t4 = 1.0 + t2 - t1
      kna = 2.511E+06*exp(29.17*t3+16.83*t4)*t6
      k1a = 1.805E-05*exp(-1.50*t3+26.92*t4)
      k2sa = 1.015E-02*exp(8.85*t3+25.14*t4)
      kw = 1.010E-14*exp(-22.52*t3+26.92*t4)
      kph = 57.639*exp(13.79*t3-5.39*t4)*t6
!cc      K3   =  5.746E-17 * EXP( -74.38 * T3 + 6.12  * T4 ) * T6 * T6
      khat = kph*k1a/kw
      kan = kna*khat

!...Compute temperature dependent equilibrium constant for NH4NO3
!...  ( from Mozurkewich, 1993)
      k3 = exp(118.87-24084.0/temp-6.025*alog(temp))

!...Convert to (micromoles/m**3) **2
      k3 = k3*convt*convt
      wh2o = 0.0
      stion = 0.0
      ah2o = 0.0
      mas = 0.0
      man = 0.0
      hplus = 0.0
      toler1 = 0.00001
      toler2 = 0.001
      nitr = 0
      nr = 0
      ratio = 0.0
      gamaan = 1.0
      gamold = 1.0

!...set the ratio according to the amount of sulfate and nitrate
      IF (tso4>minso4) THEN
        ratio = tnh4/tso4

!...If there is no sulfate and no nitrate, there can be no ammonium
!...  under the current paradigm. Organics are ignored in this version.

      ELSE

        IF (tno3==0.0) THEN

! *** If there is very little sulfate and no nitrate set concentrations
!      to a very small value and return.
          aso4 = max(floor,aso4)
          ano3 = max(floor,ano3)
          wh2o = 0.0
          ah2o = 0.0
          gnh3 = max(floor,gnh3)
          gno3 = max(floor,gno3)
          RETURN
        END IF

!...For the case of no sulfate and nonzero nitrate, set ratio to 5
!...  to send the code to the high ammonia case

        ratio = 5.0
      END IF

!....................................
!......... High Ammonia Case ........
!....................................

      IF (ratio>2.0) THEN

        gamaan = 0.1

!...Set up twice the sulfate for future use.

        twoso4 = 2.0*tso4
        xno3 = 0.0
        ynh4 = twoso4

!...Treat different regimes of relative humidity

!...ZSR relationship is used to set water levels. Units are
!...  10**(-6) kg water/ (cubic meter of air)
!...  start with ammomium sulfate solution without nitrate

! RAR: Debugging:
IF (diagn_opt .AND. icall<200) then
    WRITE(6,*) 'rpmares_old: irh,tso4,ynh4,tno3,ah2o', irh,tso4,ynh4,tno3,ah2o
    icall= icall+1
ENDIF

        CALL awater(irh,tso4,ynh4,tno3,ah2o) !**** note TNO3
        wh2o = 1.0E-3*ah2o
        aso4 = tso4*mwso4
        ano3 = 0.0
        anh4 = ynh4*mwnh4
        wfrac = ah2o/(aso4+anh4+ah2o)
!cc        IF ( WFRAC .EQ. 0.0 )  RETURN   ! No water
        IF (wfrac<0.2) THEN

!... dry  ammonium sulfate and ammonium nitrate
!...  compute free ammonia

          fnh3 = tnh4 - twoso4
          cc = tno3*fnh3 - k3

!...check for not enough to support aerosol

          IF (cc<=0.0) THEN
            xno3 = 0.0
          ELSE
            aa = 1.0
            bb = -(tno3+fnh3)
            disc = bb*bb - 4.0*cc

!...Check for complex roots of the quadratic
!...  set nitrate to zero and RETURN if complex roots are found

          IF (disc<0.0) THEN
            xno3 = 0.0
            ah2o = 1000.0*wh2o
            ynh4 = twoso4
            gno3 = tno3*mwhno3
            gnh3 = (tnh4-ynh4)*mwnh3
            aso4 = tso4*mwso4
            ano3 = 0.0
            anh4 = ynh4*mwnh4
            RETURN
          END IF

!...to get here, BB .lt. 0.0, CC .gt. 0.0 always

            dd = sqrt(disc)
            xxq = -0.5*(bb+sign(1.0,bb)*dd)

!...Since both roots are positive, select smaller root.

            xno3 = min(xxq/aa,cc/xxq)

          END IF
          ah2o = 1000.0*wh2o
          ynh4 = 2.0*tso4 + xno3
          gno3 = (tno3-xno3)*mwhno3
          gnh3 = (tnh4-ynh4)*mwnh3
          aso4 = tso4*mwso4
          ano3 = xno3*mwno3
          anh4 = ynh4*mwnh4
          RETURN

        END IF

!...liquid phase containing completely neutralized sulfate and
!...  some nitrate.  Solve for composition and quantity.

        mas = tso4/wh2o
        man = 0.0
        xno3 = 0.0
        ynh4 = twoso4
        phiold = 1.0

!...Start loop for iteration

!...The assumption here is that all sulfate is ammonium sulfate,
!...  and is supersaturated at lower relative humidities.

        DO nnn = 1, 150
          nitr = nnn
          gasqd = gamaan*gamaan
          wsqd = wh2o*wh2o
          kw2 = kan*wsqd/gasqd
          aa = 1.0 - kw2
          bb = twoso4 + kw2*(tno3+tnh4-twoso4)
          cc = -kw2*tno3*(tnh4-twoso4)

!...This is a quadratic for XNO3 [MICROMOLES / M**3] of nitrate in solut

          disc = bb*bb - 4.0*aa*cc

!...Check for complex roots, if so set nitrate to zero and RETURN

          IF (disc<0.0) THEN
            xno3 = 0.0
            ah2o = 1000.0*wh2o
            ynh4 = twoso4
            gno3 = tno3*mwhno3
            gnh3 = (tnh4-ynh4)*mwnh3
            aso4 = tso4*mwso4
            ano3 = 0.0
            anh4 = ynh4*mwnh4
!cc            WRITE( 10, * ) ' COMPLEX ROOTS '
            RETURN
          END IF

          dd = sqrt(disc)
          xxq = -0.5*(bb+sign(1.0,bb)*dd)
          rr1 = xxq/aa
          rr2 = cc/xxq

! RAR: Included from Jan:
!...Check for two non-positive roots, if so set nitrate to zero and RETURN
          IF (rr1 <= 0.0 .AND. rr2 <= 0.0) THEN
            xno3 = 0.0
            ah2o = 1000.0*wh2o
            ynh4 = twoso4
            gno3 = tno3*mwhno3
            gnh3 = (tnh4-ynh4)*mwnh3
            aso4 = tso4*mwso4
            ano3 = 0.0
            anh4 = ynh4*mwnh4
!            WRITE(*,*) 'TWO NON-POSITIVE ROOTS!!! '
            RETURN
          END IF

!...choose minimum positve root

          IF ((rr1*rr2)<0.0) THEN
            xno3 = max(rr1,rr2)
          ELSE
            xno3 = min(rr1,rr2)
          END IF
          xno3 = min(xno3,tno3)

!...This version assumes no solid sulfate forms (supersaturated )
!...  Now update water

          CALL awater(irh,tso4,ynh4,xno3,ah2o)

!...ZSR relationship is used to set water levels. Units are
!...  10**(-6) kg water/ (cubic meter of air)
!...  The conversion from micromoles to moles is done by the units of WH

          wh2o = 1.0E-3*ah2o

!...Ionic balance determines the ammonium in solution.

          man = xno3/wh2o
          mas = tso4/wh2o
          mnh4 = 2.0*mas + man
          ynh4 = mnh4*wh2o

!...MAS, MAN and MNH4 are the aqueous concentrations of sulfate, nitrate
!...  and ammonium in molal units (moles/(kg water) ).

          stion = 3.0*mas + man
          cat(1) = 0.0
          cat(2) = mnh4
          an(1) = mas
          an(2) = man
          an(3) = 0.0
          CALL actcof(cat,an,gams,molnu,phibar)
          gamaan = gams(2,2)

!...Use GAMAAN for convergence control

          eror = abs(gamold-gamaan)/gamold
          gamold = gamaan

!...Check to see if we have a solution

          IF (eror<=toler1) THEN
!cc            WRITE( 11, * ) RH, STION, GAMS( 1, 1 ),GAMS( 1, 2 ), GAMS
!cc     &      GAMS( 2, 1 ), GAMS( 2, 2 ), GAMS( 2, 3 ), PHIBAR

            aso4 = tso4*mwso4
            ano3 = xno3*mwno3
            anh4 = ynh4*mwnh4
            gno3 = (tno3-xno3)*mwhno3
            gnh3 = (tnh4-ynh4)*mwnh3
            ah2o = 1000.0*wh2o
            RETURN
          END IF

        END DO

!...If after NITR iterations no solution is found, then:

        aso4 = tso4*mwso4
        ano3 = 0.0
        ynh4 = twoso4
        anh4 = ynh4*mwnh4
        CALL awater(irh,tso4,ynh4,xno3,ah2o)
        gno3 = tno3*mwhno3
        gnh3 = (tnh4-ynh4)*mwnh3
        RETURN

      ELSE
!......................................
!......... Low Ammonia Case ...........
!......................................

!...coded by Dr. Francis S. Binkowski 12/8/91.(4/26/95)

!...All cases covered by this logic
        wh2o = 0.0
        CALL awater(irh,tso4,tnh4,tno3,ah2o)
        wh2o = 1.0E-3*ah2o
        zh2o = ah2o
!...convert 10**(-6) kg water/(cubic meter of air) to micrograms of wate
!...  per cubic meter of air (1000 g = 1 kg)

        aso4 = tso4*mwso4
        anh4 = tnh4*mwnh4
        ano3 = 0.0
        gno3 = tno3*mwhno3
        gnh3 = 0.0

!...Check for zero water.
        IF (wh2o==0.0) RETURN
        zso4 = tso4/wh2o

!...ZSO4 is the molality of total sulfate i.e. MSO4 + MHSO4

!cc         IF ( ZSO4 .GT. 11.0 ) THEN

!...do not solve for aerosol nitrate for total sulfate molality
!...  greater than 11.0 because the model parameters break down
!...  greater than  9.0 because the model parameters break down

        IF (zso4>9.0) & ! 18 June 97                        
            THEN
          RETURN
        END IF

!...First solve with activity coeffs of 1.0, then iterate.
        phiold = 1.0
        gamana = 1.0
        gamas1 = 1.0
        gamas2 = 1.0
        gamaab = 1.0
        gamold = 1.0

!...All ammonia is considered to be aerosol ammonium.
        mnh4 = tnh4/wh2o

!...MNH4 is the molality of ammonium ion.
        ynh4 = tnh4

!...loop for iteration
        DO nnn = 1, 150
          nitr = nnn

!...set up equilibrium constants including activities
!...  solve the system for hplus first then sulfate & nitrate
!          print*,'gamas,gamana',gamas1,gamas2,gamana
          rk2sa = k2sa*gamas2*gamas2/(gamas1*gamas1*gamas1)
          rkna = kna/(gamana*gamana)
          rknwet = rkna*wh2o
          t21 = zso4 - mnh4
          t221 = zso4 + t21

!...set up coefficients for cubic

          a2 = rk2sa + rknwet - t21
          a1 = rk2sa*rknwet - t21*(rk2sa+rknwet) - rk2sa*zso4 - rkna*tno3
          a0 = -(t21*rk2sa*rknwet+rk2sa*rknwet*zso4+rk2sa*rkna*tno3)

          CALL cubic(a2,a1,a0,nr,crutes)

!...Code assumes the smallest positive root is in CRUTES(1)

          hplus = crutes(1)
          bal = hplus**3 + a2*hplus**2 + a1*hplus + a0
          mso4 = rk2sa*zso4/(hplus+rk2sa) ! molality of sulfat
          mhso4 = zso4 - & ! molality of bisulf
            mso4
          mna = rkna*tno3/(hplus+rknwet) ! molality of nitrat
          mna = max(0.0,mna)
          mna = min(mna,tno3/wh2o)
          xno3 = mna*wh2o
          ano3 = mna*wh2o*mwno3
          gno3 = (tno3-xno3)*mwhno3

!...Calculate ionic strength
          stion = 0.5*(hplus+mna+mnh4+mhso4+4.0*mso4)

!...Update water
          CALL awater(irh,tso4,ynh4,xno3,ah2o)

!...Convert 10**(-6) kg water/(cubic meter of air) to micrograms of wate
!...  per cubic meter of air (1000 g = 1 kg)

          wh2o = 1.0E-3*ah2o
          cat(1) = hplus
          cat(2) = mnh4
          an(1) = mso4
          an(2) = mna
          an(3) = mhso4
!          print*,'actcof',cat(1),cat(2),an(1),an(2),an(3),gams,molnu,phibar
          CALL actcof(cat,an,gams,molnu,phibar)

          gamana = gams(1,2)
          gamas1 = gams(1,1)
          gamas2 = gams(1,3)
          gamaan = gams(2,2)

          gamahat = (gamas2*gamas2/(gamaab*gamaab))
          bhat = khat*gamahat
!cc          EROR = ABS ( ( PHIOLD - PHIBAR ) / PHIOLD )
!cc          PHIOLD = PHIBAR
          eror = abs(gamold-gamahat)/gamold
          gamold = gamahat

!...write out molalities and activity coefficient
!...  and return with good solution

          IF (eror<=toler2) THEN
!cc            WRITE(12,*) RH, STION,HPLUS,ZSO4,MSO4,MHSO4,MNH4,MNA
!cc            WRITE(11,*) RH, STION, GAMS(1,1),GAMS(1,2),GAMS(1,3),
!cc     &                  GAMS(2,1),GAMS(2,2),GAMS(2,3), PHIBAR
            RETURN
          END IF

        END DO

!...after NITR iterations, failure to solve the system, no ANO3

        gno3 = tno3*mwhno3
        ano3 = 0.0
        CALL awater(irh,tso4,tnh4,tno3,ah2o)
        RETURN

      END IF
! ratio .gt. 2.0
END SUBROUTINE rpmares_old

!ia*********************************************************
!ia                                                        *
!ia BEGIN OF AEROSOL ROUTINE				   *
!ia							   *
!ia*********************************************************

!***********************************************************************
!   	BEGIN OF AEROSOL CALCULATIONS
!***********************************************************************
!ia  									*
!ia	MAIN AEROSOL DYNAMICS ROUTINE					*
!ia	based on MODELS3 formulation by FZB				*
!ia	Modified by IA in May 97					*
!ia     THIS PROGRAMME IS THE LINK BETWEEN GAS PHASE AND AEROSOL PHASE
!ia     CALCULATIONS IN THE COLUMN MODEL. IT CONVERTS ALL DATA AND
!ia     VARIABLES BETWEEN THE TWO PARTS AND DRIVES THE AEROSOL
!ia     CALCULATIONS.
!ia     INPUT DATA REQUIRED FOR AEROSOL DYNAMICS ARE SET UP HERE FOR
!ia     ONE GRID CELL!!!!
!ia     and passed to dynamics calcs. subroutines.
!ia									*
!ia	Revision history						*
!ia	When	WHO	WHAT						*
!ia	----	----	----						*
!ia	????	FZB	BEGIN						*
!ia	05/97	IA	Adapted for use in CTM2-S			*
!ia			Modified renaming/bug fixing			*
!ia     11/97   IA      Modified for new model version
!ia                     see comments under iarev02
!ia     03/98   IA      corrected error on pressure units
!ia									*
!ia	Called BY:	CHEM						*
!ia									*
!ia	Calls to:	OUTPUT1,AEROPRC					*
!ia									*
!ia*********************************************************************

! end RPMares
! convapr_in is removed, it wasn't used indeed
    SUBROUTINE rpmmod3(layer,dtsec,pres,temp,relhum, &
                       nitrate_in,nh3_in,vsulf_in,so4rat_in, &
                       drog_in,                              &
                       eeci_in,eecj_in,eorgi_in,eorgj_in,    &
                       epm25i,epm25j,epmcoarse,soilrat_in,   & 
                       cblk,igrid,jgrid,kgrid,brrto,         &
                       diagn_opt,do_isorropia,do_n2o5_hetchem )

!     Includes:
!iarev02       INCLUDE  AEROINCL.EXT 
! block size, set to 1 in column model  ciarev0
!      INTEGER blksize
!ia   			  kept to 1 in current version of column model
! actual number of cells in arrays ( default is

LOGICAL :: diagn_opt, do_isorropia
INTEGER :: do_n2o5_hetchem

      INTEGER layer
! number of layer (default is 1 in

! index for cell in blocked array (default is 1 in
      INTEGER, PARAMETER :: ncell=1
! *** inputs
! Input temperature [ K ]                      
      REAL temp
! Input relative humidity  [ fraction ]        
      REAL relhum
! Input pressure [ hPa ]                       
      REAL pres
! Input number for Aitken mode [ m**-3 ]       
      REAL numnuc_in
! Input number for accumulation mode [ m**-3 ] 
      REAL numacc_in
! Input number for coarse mode  [ m**-3 ]      
      REAL numcor_in
                         ! sulfuric acid [ ug m**-3 ]
      REAL vsulf_in
! total sulfate vapor as sulfuric acid as      
                         ! sulfuric acid [ ug m**-3 ]
      REAL asulf_in
! total sulfate aerosol as sulfuric acid as    
! i-mode sulfate input as sulfuric acid [ ug m*
      REAL asulfi_in
! ammonia gas [  ug m**-3 ]                    
      REAL nh3_in
! input value of nitric acid vapor [ ug m**-3 ]
      REAL nitrate_in
! Production rate of sulfuric acid   [ ug m**-3
      REAL so4rat_in
                         ! aerosol [ ug m**-3 s**-1 ]
      REAL soilrat_in
! Production rate of soil derived coarse       
! Emission rate of i,j-mode EC [ug m**-3 s**-1]
      REAL eeci_in, eecj_in

! Emission rate of i,j-mode org. aerosol [ug m**-
      REAL eorgi_in, eorgj_in

! Emission rate of j-mode org. aerosol [ug m**-
! total # of cond. vapors & SOA species 
!      INTEGER ncv
! # of anthrop. cond. vapors & SOA speci
!      INTEGER nacv
! # of organic aerosol precursor        
      REAL drog_in(ldrog_soa2)                                 ! organic aerosol precursor [ppm]
! Input delta ROG concentration of      
!      REAL condvap_in(ncv) ! cond. vapor input [ug m^-3]
!      REAL drog(blksize,ldrog_vbs)                                 ! organic aerosol precursor [ppm]
      REAL brrto        ! Branching ratio

! *** Primary emissions rates: [ ug / m**3 s ]

! *** emissions rates for unidentified PM2.5 mass
! Delta ROG concentration of            
      REAL epm25i  ! Aitken mode
      REAL epm25j
! *** emissions rates for primary organic aerosol
! Accumululaton mode                  
      REAL eorgi   ! Aitken mode
      REAL eorgj
! *** emissions rates for elemental carbon
! Accumululaton mode                   
      REAL eeci    ! Aitken mode
      REAL eecj
! *** Primary emissions rates [ ug m**-3 s -1 ] :

! Accumululaton mode                    
      REAL epm25(blksize) ! emissions rate for PM2.5 mass           
      REAL esoil(blksize) ! emissions rate for soil derived coarse a
      REAL eseas(blksize) ! emissions rate for marine coarse aerosol
      REAL epmcoarse(blksize)
!ts added pollen
      REAL epolp(blksize) ! emissions rate for biological coarse aerosol
      REAL epols(blksize) ! ts emissions rate for biological coarse aerosol
! emissions rate for anthropogenic coarse 

      REAL dtsec
! time step [ s ], PASSED FROM MAIN COLUMN MODE

      REAL newm3
      REAL totaersulf
! total aerosol sulfate                   
! loop index for time steps                     
      INTEGER numsteps
      REAL step

! *** arrays for aerosol model codes:
!      INTEGER nspcsda

! number of species in CBLK ciarev02           
      REAL cblk(blksize,nspcsda)    ! RAR: the unit for cblk aeorosol species is ug/m3

! *** Meteorological information in blocked arays:

! *** Thermodynamic variables:

! main array of variables            
      REAL blkta(blksize)    ! Air temperature [ K ]
      REAL blkprs(blksize)   ! Air pressure in [ Pa ]
      REAL blkdens(blksize)  ! Air density  [ kg m^-3 ]
      REAL blkrh(blksize)    ! Relative humidity

! *** Chemical production rates [ ug m**-3 s -1 ] :

! Fractional relative humidity              
      REAL so4rat(blksize)                                 ! rate [ug/m^3/s]
! sulfuric acid vapor-phase production  
                                       ! production rate from aromatics [ ug /

! POA: production rates [ ug /m3 s ]
      REAL orgpoa1rat(blksize), orgpoa2rat(blksize), orgpoa3rat(blksize)

! production rate for ASOA [ ug /m3 s ]
      REAL organt0rat(blksize), organt1rat(blksize), organt2rat(blksize), organt3rat(blksize)

! production rate for ASOA [ ug /m3 s ]
      REAL orgbio1rat(blksize), orgbio2rat(blksize), orgbio3rat(blksize)

! biogenic organic aerosol production   
      REAL xlm(blksize) ! atmospheric mean free path [ m ]  
      REAL amu(blksize) 
! *** aerosol properties:

! *** modal diameters:
! atmospheric dynamic viscosity [ kg
      REAL dgnuc(blksize) ! nuclei mode geometric mean diamete
      REAL dgacc(blksize) ! accumulation geometric mean diamet
      REAL dgcor(blksize) 

! *** Modal mass concentrations [ ug m**3 ]
! coarse mode geometric mean diamete
      REAL pmassn(blksize) ! mass concentration in Aitken mode 
      REAL pmassa(blksize) ! mass concentration in accumulation
      REAL pmassc(blksize) 
! *** average modal particle densities  [ kg/m**3 ]

! mass concentration in coarse mode 
      REAL pdensn(blksize) ! average particle density in nuclei
      REAL pdensa(blksize) ! average particle density in accumu
      REAL pdensc(blksize) 
! *** average modal Knudsen numbers

! average particle density in coarse
      REAL knnuc(blksize) ! nuclei mode  Knudsen number       
      REAL knacc(blksize) ! accumulation Knudsen number       
      REAL kncor(blksize) 
! *** reciprocal modal condensation rates for sulfuric acid [ 1/s ]

! coarse mode  Knudsen number       
      REAL fconcn(blksize) 
! reciprocal condensation rate Aitke
      REAL fconca(blksize) !bs
! reciprocal condensation rate acclu
      REAL fconcn_org(blksize)
      REAL fconca_org(blksize)

! *** production of new mass concentration [ ug/m**3 s ]
      REAL dmdt(blksize) !                                 by particle formation

! *** production of new number concentration [ number/m**3 s ]

! rate of production of new mass concen
      REAL dndt(blksize) !                                 by particle formation
! *** growth rate for third moment by condensation of precursor
!      vapor on existing particles [ 3rd mom/m**3 s ]

! rate of producton of new particle num
      REAL cgrn3(blksize)       !  Aitken mode
      REAL cgra3(blksize) 
! *** Rates for coaglulation: [ m**3/s ]

! *** Unimodal Rates:

!  Accumulation mode                    
      REAL urn00(blksize) ! Aitken mode 0th moment self-coagulation ra
      REAL ura00(blksize) 

! *** Bimodal Rates:  Aitken mode with accumulation mode ( Aitken mode)
! accumulation mode 0th moment self-coagulat
      REAL brna01(blksize) ! rate for 0th moment                     
      REAL brna31(blksize) 
! *** other processes

! rate for 3rd moment                     
      REAL deltaso4a(blksize) !                                  sulfate aerosol by condensation   [ u

! *** housekeeping variables:
! increment of concentration added to
      INTEGER, PARAMETER :: unit=30
      CHARACTER*16 pname
      PARAMETER (pname=' BOX            ')
      INTEGER isp,igrid,jgrid,kgrid

      INTEGER, SAVE :: icall

! loop index for species.                             
      INTEGER ii !, iimap(8)
!      DATA iimap/1, 2, 18, 19, 21, 22, 23, 24/

      step = dtsec    ! set time step

      blkta(blksize) = temp     ! T in Kelvin

      blkprs(blksize)= pres*100. ! P in  Pa (pres is given in

      blkrh(blksize) = relhum ! fractional RH

      blkdens(blksize) = blkprs(blksize)/(rdgas*blkta(blksize))   !rs      CBLK(blksize,VSULF) = vsulf_in

      so4rat(blksize) = so4rat_in

!...INITIALISE EMISSION RATES

!     epm25i(blksize) = & ! unidentified PM2.5 mass                  
!       0.0
!     epm25j(blksize) = & 
!       0.0
! unidentified PM2.5 m
      eorgi= eorgi_in
      eorgj= eorgj_in

      eeci= eeci_in
      eecj= eecj_in

      epm25 = 0.0       !currently from input file ACTIONIA

      esoil(blksize) = & ! ACTIONIA
        soilrat_in
      eseas(blksize) = & !currently from input file ACTIONIA
        0.0
!     epmcoarse(blksize) = & !currently from input file ACTIONIA    
!       0.0
!ts
      epolp(blksize) = & !#DOUBT currently from input file ACTIONIA        
        0.0
      epols(blksize) = & !#DOUBT currently from input file ACTIONIA        
        0.0
!ts pollen
      dgnuc(blksize) = dginin
      dgacc(blksize) = dginia
      dgcor(blksize) = dginic
      newm3 = 0.0

! *** Set up initial total 3rd moment factors
      totaersulf = 0.0

! *** Call aerosol routines
      CALL aeroproc(layer,cblk,step,blkta,blkprs,blkdens,blkrh,so4rat, &
             orgpoa1rat,orgpoa2rat,orgpoa3rat, &
             organt0rat,organt1rat,organt2rat,organt3rat, &
             orgbio1rat,orgbio2rat,orgbio3rat, &
             drog_in,epm25i,epm25j,eorgi,eorgj,eeci,eecj, &
             epmcoarse,esoil,eseas,epolp,epols,xlm,amu,    &
             dgnuc,dgacc,dgcor,pmassn,pmassa,pmassc,pdensn,pdensa,pdensc,      &
             knnuc,knacc,kncor,fconcn,fconca,fconcn_org,fconca_org,dmdt,dndt,cgrn3,cgra3, &
             urn00,ura00,brna01,brna31,deltaso4a,igrid,jgrid,kgrid,brrto,diagn_opt,       &
             do_isorropia, do_n2o5_hetchem)

! RAR: Debugging:
IF (diagn_opt .AND. icall<100) then
    WRITE(6,*) 'rpmmod3: blkta,blkprs,blkdens,blkrh', blkta,blkprs,blkdens,blkrh
    icall=icall+1
ENDIF

      RETURN
END SUBROUTINE rpmmod3
!---------------------------------------------------------------------------

SUBROUTINE soa2_vbs( blkta,blkprs,                     &
                     orgpoa1rat,orgpoa2rat,orgpoa3rat, &
                     organt0rat,organt1rat,organt2rat,organt3rat, &
                     orgbio1rat,orgbio2rat,orgbio3rat, &
                     drog_in,cblk,dt,brrto,diagn_opt, &
                     igrid,jgrid,kgrid )

!bs  SOA_VBS calculates the formation and partitioning of secondary  !
!bs  organic aerosol based on (pseudo-)ideal solution thermodynamics.  !
!bs                                                                    !
!sam The original Schell model (JGR, vol. 106, D22, 28275-28293, 2001) !
!sam is modified drastically to incorporate the SOA vapor-pressure     !
!sam basis set approach developed by Carnegie Mellon folks.            !
!sam Recommended changes according to Allen Robinson, 9/15/09          !
!sam The treatment is done very similar to Lane et al., Atmos. Envrn., !
!sam vol 42, 7439-7451, 2008.                                          !
!sam Four basis vapor-pressures for anthropogenic and 4 basis vp's     !
!sam for biogenic SOA are used.  The SAPRC-99 yield information for    !
!sam low and high NOx conditions (Lane, Donahue and Pandis, ES&T,      !
!sam vol. 42, 6022-6027, 2008) are mapped to RADM2/RACM species.       !
!sam                                                                   !
!sam Basis vapor pressures (@ 300K)
!    Anthro (0.1 ug/m3) -    asoa0
!sam Anthro (1 ug/m3)   -    asoa1   Biogenic (1 ug/m3)   -    bsoa1   !
!sam Anthro (10 ug/m3)  -    asoa2   Biogenic (10 ug/m3)  -    bsoa2   !
!sam Anthro (100 ug/m3) -    asoa3   Biogenic (100 ug/m3) -    bsoa3   !
!bs                                                                    !
!bs  This code considers two cases:                                    !
!bs   i) initil absorbing mass is existend in the aerosol phase        !
!bs  ii) a threshold has to be exeeded before partitioning (even below !
!bs      saturation) will take place.                                  !
!bs                                                                    !
!bs  The temperature dependence of the saturation concentrations are   !
!bs  calculated using the Clausius-Clapeyron equation.                 !
!bs                                                                    !
!bs  Input files:   None                                               !
!bs  Output files:  None                                               !
! Initially developed by Schell
! 052011   McKeen/Ahmadov   Subroutine development            !
! RAR May, 2013: This subroutine has been substantially modified to     !
! handle semi-volatile POA and reduced bins

      USE module_configure, only: grid_config_rec_type

      LOGICAL :: diagn_opt

      REAL, INTENT(OUT) :: brrto   ! branching ratio for NOx conditions

      ! lsv_ocv=12
      REAL ctot(lsv_ocv)  !bs total conc. of cond. vapor aerosol +
      REAL cgas(lsv_ocv)  !bs gas phase concentration of cond. vapors
      REAL caer(lsv_ocv)  !bs aerosol phase concentration of cond.
      REAL aold(lsv_ocv)  !bs saved CAER for rate determination

      REAL cblk(blksize,nspcsda)              ! main array of variables
      REAL blkta(blksize), blkprs(blksize)    ! Air temperature [ K ] and air pressure in [ Pa ]

! POA growth rates:
      REAL orgpoa1rat(blksize), orgpoa2rat(blksize), orgpoa3rat(blksize)
! ASOA:
      REAL organt0rat(blksize), organt1rat(blksize), organt2rat(blksize), organt3rat(blksize)
! BSOA:
      REAL orgbio1rat(blksize), orgbio2rat(blksize), orgbio3rat(blksize)

! Delta ROG conc. [ppm]:
      REAL drog_in(ldrog_soa2)      !blksize=1, ldrog_vbs=9+1, the last ldrog_vbs is actually is the branching ratio

!sam reference temperature T0 = 300 K, a change from original 298K
      REAL, PARAMETER :: tnull=300.

!     REAL, PARAMETER :: CABSMIN=.00001   ! Minimum amount of absorbing material - needed in iteration method
!      REAL, PARAMETER :: CABSMIN= 0.001?  ! Higher CABSMIN=0.01, check this later, it isn't used here

      ! we have 2 type of SOA - anthropogenic and biogenic
      !sm number of SAPRC species variables in CMU lumped partitioning table
      !sm 1=ALK4(hc5),2=ALK5(hc8),3=OLE1(ol2),4=OLE2(oli),5=ARO1(tol)
      !sm 6=AOR2(xyl),7=ISOP(iso),8=SESQ(?),9=TERP(alp)

      INTEGER n, l, ll, bn, cls      !bs loop indices
      INTEGER, PARAMETER :: lcell=1
      INTEGER igrid,jgrid,kgrid

      REAL dt                 ! model time step in  SECONDS
      REAL convfac            !bs conversion factor ppm --> ug/m^3
      REAL ttinv              !bs difference of inverse temperatures
      REAL minit              !bs initial organic absorbing mass [ug/m^3]
      REAL mtot               !bs total organic mass [ug/m^3]
      REAL prod(lsv_ocv)      ! production of OCVs from VOC oxidation, lsv_ocv=12
      REAL f(ldrog_soa2-2)    ! scaling factor for prod, for testing purposes only

! mass-based stoichometric yield for product i and csti is the effective saturation
      REAL alphlowN (Nbin_voc_ox,nsaprc)    ! SOA yields for low NOx conditions, Nbin_voc_ox=4, nsaprc=9
      REAL alphhiN  (Nbin_voc_ox,nsaprc)    ! SOA yields for high NOx conditions
      REAL alphai   (Nbin_voc_ox,nsaprc)    ! SOA yields calculed for a given condition

      REAL mw_voc(nsaprc)         ! molecular weight of the SOA precusors

      REAL PnGtotal,dum,dum2,fmtot,fmtot2    ! Real constants used in Newton iteration

      REAL mw_ocv(lsv_ocv), dhvap(lsv_ocv)   ! RAR molecular weight of OCVs [g/mol], lsv_ocv=12
      REAL csat00(lsv_ocv), csat(lsv_ocv)    ! Saturation concentrations

      INTEGER, SAVE :: icall

! RAR: For consistency all the species are assumed to have 4 bins, however partitioning occurs only for first 3 bins!

      ! Enthalpy of vaporization (J/mol), Murphy and Pandis, 2009; Tsimpidi et al., 2010;
      ! 4- CVPOA, 4- CVASOA, 4- CVBSOA
      DATA csat00 /1.0, 1.0E01, 1.0E02, 1.0E03, 1.0, 1.0E01, 1.0E02, 1.0E03, 1.0, 1.0E01, 1.0E02, 1.0E03/      ! Volatility bins at 300K, for which partitioning takes place
      DATA dhvap /100.0E03, 94.0E03, 88.0E03, 82.0E03, 4*30.0E03, 4*30.0E03/     ! J/mol
      DATA mw_ocv /4*250., 4*150., 4*180./     !RAR: Molecular weights [g/mol] of OCVs as in Murphy and Pandis, 2009

! SAM:  from Murphy et al. 2009; SOA yields at 298K
! Vol. bins: 1       10      100     1000
!      DATA alphlowN /   &
!             0.0000, 0.0750, 0.0000, 0.0000,   & ! 1) ALK4
!             0.0000, 0.3000, 0.0000, 0.0000,   & ! 2) ALK5
!             0.0890, 0.5500, 0.2000, 0.0000,   & ! 3) HC16 yield for IVOCs
!             0.0045, 0.0090, 0.0600, 0.2250,   & ! 4) OLE1
!             0.0225, 0.0435, 0.1290, 0.3750,   & ! 5) OLE2
!             0.0750, 0.2250, 0.3750, 0.5250,   & ! 6) ARO1
!             0.0750, 0.3000, 0.3750, 0.5250,   & ! 7) ARO2
!             ! biogenic
!             0.0090, 0.0300, 0.0150, 0.0000,   & ! 8) ISOP
!             0.0750, 0.1500, 0.7500, 0.9000,   & ! 9) SESQ
!             0.1073, 0.0918, 0.3587, 0.6075    / ! 10)TERP
!
!      DATA alphhiN /    &
!             0.0000, 0.0375, 0.0000, 0.0000,   & ! 1) ALK4
!             0.0000, 0.1500, 0.0000, 0.0000,   & ! 2) ALK5
!             0.0890, 0.5500, 0.2000, 0.0000,   & ! 3) HC16 yield for IVOCs, no NOx dependence
!             0.0008, 0.0045, 0.0375, 0.1500,   & ! 4) OLE1
!             0.0030, 0.0255, 0.0825, 0.2700,   & ! 5) OLE2
!             0.0030, 0.1650, 0.3000, 0.4350,   & ! 6) ARO1
!             0.0015, 0.1950, 0.3000, 0.4350,   & ! 7) ARO2
!             ! biogenic
!             0.0003, 0.0225, 0.0150, 0.0000,   & ! 8) ISOP
!             0.0750, 0.1500, 0.7500, 0.9000,   & ! 9) SESQ
!             0.0120, 0.1215, 0.2010, 0.5070    / ! 10 TERP
!
!      DATA mw_voc / 73.23,   &    ! 1) ALK4        ! 10 VOC species
!                    106.97,  &    ! 2) ALK5
!                    137.00,  &    ! 3) IVOC, provided in the EPA NEI-2011 inventory for unidentified, for hexadecane is 226 g/mol
!                    61.68,   &    ! 4) OLE1
!                    79.05,   &    ! 5) OLE2
!                    100.47,  &    ! 6) ARO1
!                    113.93,  &    ! 7) ARO2
!                    68.12,   &    ! 8) ISOP
!                    204.0,   &    ! 9) SESQ
!                    136.24   /    ! 10 TERP

      DATA alphlowN /   &
             0.0000, 0.1875, 0.0000, 0.0000,   & ! 1) SOAALK
             0.0135, 0.02625, 0.0945, 0.300,   & ! 2) OLE
             0.0750, 0.2625, 0.3750, 0.5250,   & ! 3) ARO
             0.0890, 0.5500, 0.2000, 0.0000,   & ! 4) HC16 yield for IVOCs
             ! biogenic
             0.0090, 0.0300, 0.0150, 0.0000,   & ! 5) ISOP
             0.0750, 0.1500, 0.7500, 0.9000,   & ! 6) SESQ
             0.1073, 0.0918, 0.3587, 0.6075    / ! 7)TERP

      DATA alphhiN /    &
             0.0000, 0.09375, 0.0000, 0.0000,  & ! 1) SOAALK
             0.0019, 0.0150, 0.0600, 0.2100,   & ! 3) OLE
             0.00225, 0.180, 0.3000, 0.4350,   & ! 4) ARO
             0.0890, 0.5500, 0.2000, 0.0000,   & ! 2) HC16 yield for IVOCs, no NOx dependence             
! biogenic
             0.0003, 0.0225, 0.0150, 0.0000,   & ! 5) ISOP
             0.0750, 0.1500, 0.7500, 0.9000,   & ! 6) SESQ
             0.0120, 0.1215, 0.2010, 0.5070    / ! 70 TERP

      DATA mw_voc / 90.1,   &     ! 1) SOAALK        ! 7 VOC species
                    70.365,   &    ! 2) OLE
                    107.2,  &    ! 3) ARO
                    137.00,  &    ! 4) IVOC, provided in the EPA NEI-2011 inventory for unidentified, for hexadecane is 226 g/mol                    
                    68.12,   &    ! 5) ISOP
                    204.0,   &    ! 6) SESQ
                    136.24   /    ! 7 TERP

    ! scaling factors, for testing purposes, check TOL and ISO only
    ! 05/23/2011: for testing all are zero!
     f(psoaalk) = 1.
     f(pole) = 1.
     f(paro) = 1.
!    f(palk4) = 1.
!    f(palk5) = 1.
    f(pivoc) = 1.
!    f(pole1) = 1.
!    f(pole2) = 1.
!    f(paro1) = 1.
!    f(paro2) = 1.
    f(pisop) = 1.
    f(psesq) = 1.
    f(pterp) = 1.

      DO ll= 1,ldrog_soa2 - 2
          drog_in(ll)= f(ll)*drog_in(ll)
      END DO

      ! calculation of the yields using the branching ratio
      brrto= drog_in(pbrch)    ! temporary variable for the branching ratio
      DO bn=1,Nbin_voc_ox      ! Nbin
         DO cls=1,nsaprc       ! SOA precursor classes=10
            alphai(bn,cls)= mw_voc(cls)*( alphhiN(bn,cls)*brrto + alphlowN(bn,cls)*(1.-brrto) )
         ENDDO
      ENDDO

      ttinv= 1./tnull - 1./blkta(lcell)
      convfac= blkprs(lcell)/(rgasuniv*blkta(lcell))

      minit= cblk(lcell,vpoa0j) + cblk(lcell,vpoa0i)
      minit= AMAX1(minit,aeroconcmin)   ! aeroconcmin=0.0001

! RAR: Debugging:
IF (diagn_opt .AND. icall<100 .AND. kgrid<3 ) then
    WRITE(6,*) 'soa2_vbs: ipoa1,ipoa4,iasoa1,iasoa4,ibsoa1,ibsoa4,lsv_ocv ',ipoa1,ipoa4,iasoa1,iasoa4,ibsoa1,ibsoa4,lsv_ocv
    WRITE(6,*) 'soa2_vbs: vcvpoa1,vcvasoa1,vcvbsoa4 ',vcvpoa1,vcvasoa1,vcvbsoa4
    WRITE(6,*) 'soa2_vbs: alphai(1,:)', alphai(1,:)
!    WRITE(6,*) 'soa2_vbs: alphhiN(8,2),alphlowN(10,3)', alphhiN(8,2),alphlowN(10,3)
    WRITE(6,*) 'soa2_vbs: mw_voc(1)', mw_voc(1)
    WRITE(6,*) 'soa2_vbs: nsaprc,Nbin_voc_ox,lcell ,', nsaprc,Nbin_voc_ox,lcell
    WRITE(6,*) 'soa2_vbs: blkta,blkprs ', blkta,blkprs
    WRITE(6,*) 'soa2_vbs: ttinv,convfac,dt', ttinv,convfac,dt
END IF

      ! cblk for gases comes in ppmv, we get the density in ug/kg_air
      ! by multiplying it by (convfac=rho_air/mw_air)x mw_ocv
      cgas(ipoa1)  = cblk(lcell,vcvpoa1)*convfac* mw_ocv(ipoa1)
      cgas(ipoa2)  = cblk(lcell,vcvpoa2)*convfac* mw_ocv(ipoa2)
      cgas(ipoa3)  = cblk(lcell,vcvpoa3)*convfac* mw_ocv(ipoa3)
      cgas(ipoa4)  = cblk(lcell,vcvpoa4)*convfac* mw_ocv(ipoa4)

      cgas(iasoa1) = cblk(lcell,vcvasoa1)*convfac* mw_ocv(iasoa1)
      cgas(iasoa2) = cblk(lcell,vcvasoa2)*convfac* mw_ocv(iasoa2)
      cgas(iasoa3) = cblk(lcell,vcvasoa3)*convfac* mw_ocv(iasoa3)
      cgas(iasoa4) = cblk(lcell,vcvasoa4)*convfac* mw_ocv(iasoa4)

      cgas(ibsoa1) = cblk(lcell,vcvbsoa1)*convfac* mw_ocv(ibsoa1)
      cgas(ibsoa2) = cblk(lcell,vcvbsoa2)*convfac* mw_ocv(ibsoa2)
      cgas(ibsoa3) = cblk(lcell,vcvbsoa3)*convfac* mw_ocv(ibsoa3)
      cgas(ibsoa4) = cblk(lcell,vcvbsoa4)*convfac* mw_ocv(ibsoa4)

      ! cblk for aerosols come in density (ug/kg_air), converted already in soa2_driver
      caer(ipoa1)  = cblk(lcell,vpoa1j) + cblk(lcell,vpoa1i)
      caer(ipoa2)  = cblk(lcell,vpoa2j) + cblk(lcell,vpoa2i)
      caer(ipoa3)  = cblk(lcell,vpoa3j) + cblk(lcell,vpoa3i)
      caer(ipoa4)  = 0.0

      caer(iasoa1) = cblk(lcell,vasoa1j) + cblk(lcell,vasoa1i)
      caer(iasoa2) = cblk(lcell,vasoa2j) + cblk(lcell,vasoa2i)
      caer(iasoa3) = cblk(lcell,vasoa3j) + cblk(lcell,vasoa3i)
      caer(iasoa4) = 0.0

      caer(ibsoa1) = cblk(lcell,vbsoa1j) + cblk(lcell,vbsoa1i)
      caer(ibsoa2) = cblk(lcell,vbsoa2j) + cblk(lcell,vbsoa2i)
      caer(ibsoa3) = cblk(lcell,vbsoa3j) + cblk(lcell,vbsoa3i)
      caer(ibsoa4) = 0.0

      ! Production of SOA by oxidation of VOCs
      ! There are 6 classes of the precursors for anthropogenic SOA
      prod(ipoa1:ipoa4)= 0.0;

!      prod(iasoa1) = alphai(1,1)*drog_in(palk4) + alphai(1,2)*drog_in(palk5) + &
!                                                + alphai(1,3)*drog_in(pivoc) + &    ! RAR:
!                     alphai(1,4)*drog_in(pole1) + alphai(1,5)*drog_in(pole2) + &
!                     alphai(1,6)*drog_in(paro1) + alphai(1,7)*drog_in(paro2)
!
!      prod(iasoa2) = alphai(2,1)*drog_in(palk4) + alphai(2,2)*drog_in(palk5) + &
!                                                + alphai(2,3)*drog_in(pivoc) + &
!                     alphai(2,4)*drog_in(pole1) + alphai(2,5)*drog_in(pole2) + &
!                     alphai(2,6)*drog_in(paro1) + alphai(2,7)*drog_in(paro2)
!
!      prod(iasoa3) = alphai(3,1)*drog_in(palk4) + alphai(3,2)*drog_in(palk5) + &
!                                                + alphai(3,3)*drog_in(pivoc) + &
!                     alphai(3,4)*drog_in(pole1) + alphai(3,5)*drog_in(pole2) + &
!                     alphai(3,6)*drog_in(paro1) + alphai(3,7)*drog_in(paro2)
!
!      prod(iasoa4) = alphai(4,1)*drog_in(palk4) + alphai(4,2)*drog_in(palk5) + &
!                                                + alphai(4,3)*drog_in(pivoc) + &
!                     alphai(4,4)*drog_in(pole1) + alphai(4,5)*drog_in(pole2) + &
!                     alphai(4,6)*drog_in(paro1) + alphai(4,7)*drog_in(paro2)

      prod(iasoa1) = alphai(1,1)*drog_in(psoaalk) + &
                     alphai(1,2)*drog_in(pole)    + &
                     alphai(1,3)*drog_in(paro)    + &
                     alphai(1,4)*drog_in(pivoc)
      prod(iasoa2) = alphai(2,1)*drog_in(psoaalk) + &
                     alphai(2,2)*drog_in(pole)    + &
                     alphai(2,3)*drog_in(paro)    + &
                     alphai(2,4)*drog_in(pivoc)
      prod(iasoa3) = alphai(3,1)*drog_in(psoaalk) + &
                     alphai(3,2)*drog_in(pole)    + &
                     alphai(3,3)*drog_in(paro)    + &
                     alphai(3,4)*drog_in(pivoc)
      prod(iasoa4) = alphai(4,1)*drog_in(psoaalk) + &
                     alphai(4,2)*drog_in(pole)    + &
                     alphai(4,3)*drog_in(paro)    + &
                     alphai(4,4)*drog_in(pivoc)

! There are 3 classes of the precursors for biogenic SOA;
      prod(ibsoa1) = alphai(1,5)*drog_in(pisop) + alphai(1,6)*drog_in(psesq) + &
                     alphai(1,7)*drog_in(pterp)

      prod(ibsoa2) = alphai(2,5)*drog_in(pisop) + alphai(2,6)*drog_in(psesq) + &
                     alphai(2,7)*drog_in(pterp)

      prod(ibsoa3) = alphai(3,5)*drog_in(pisop) + alphai(3,6)*drog_in(psesq) + &
                     alphai(3,7)*drog_in(pterp)

      prod(ibsoa4) = alphai(4,5)*drog_in(pisop) + alphai(4,6)*drog_in(psesq) + &
                     alphai(4,7)*drog_in(pterp)

! RAR: Debugging:
IF (diagn_opt .AND. icall<200 .AND. kgrid<3 ) then
    WRITE(6,*) 'soa2_vbs: cgas ',cgas
    WRITE(6,*) 'soa2_vbs: caer ',caer
    WRITE(6,*) 'soa2_vbs: prod ',prod
    WRITE(6,*) 'soa2_vbs: drog_in ',drog_in
END IF

!bs * calculate actual production from gasphase reactions [ug/m^3]
!bs * calculate vapor pressure of pure compound as a liquid using the Clausius-Clapeyron equation and the actual saturation concentration.
!bs * calculate the threshold for partitioning if no initial mass is present to partition into.

PnGtotal=0.  ! track total SOA+OCVs over bins for limits on Newton Iteration of total SOA mass
      DO ll=1,lsv_ocv                          ! we've total 3*4 bins for SOA production for VOC oxidation, no alpha is needed here
          prod(ll)= convfac*prod(ll)           ! get in density units (ug/kg_air) from ppmv, (convfac=rho_air/mu_air)
          ctot(ll)= prod(ll) + cgas(ll) + caer(ll)
          aold(ll)= caer(ll)
          csat(ll)= csat00(ll)* tnull/ blkta(lcell)* exp( dhvap(ll)/rgasuniv*ttinv )
          PnGtotal= PnGtotal + ctot(ll)
      END DO

! RAR: Debugging:
IF (diagn_opt .AND. icall<200 .AND. kgrid<3 ) then
    WRITE(6,*) 'soa2_vbs: mw_ocv',mw_ocv
    WRITE(6,*) 'soa2_vbs: csat ',csat
    WRITE(6,*) 'soa2_vbs: prod ',prod
    WRITE(6,*) 'soa2_vbs: ctot ',ctot
    WRITE(6,*) 'soa2_vbs: aold ',aold
END iF

!SAM 9/8/09 - Include absorbing SOA material within aerosols in calculation !
! mtot is initial guess to SOA mass (aerosol plus extra absorbing mass (minit))
      mtot= 0.0
      DO ll=1,lsv_ocv
         mtot= mtot + AMIN1(1.,ctot(ll)/csat(ll))*ctot(ll)
      ENDDO
      mtot= mtot + minit

! RAR: Debugging:
IF (diagn_opt .AND. icall<200 .AND. kgrid<3 ) then
    WRITE(6,*) 'soa2_vbs: mtot,minit ',mtot,minit
END IF
      
!SAM: Find total SOA mass from newton iteration, needs only 5 iterations for exact solution
loop_newt: DO l=1,5   ! Fixed Newton iteration number
              fmtot=0.
              fmtot2=0.
              DO ll=1,lsv_ocv
                  dum= ctot(ll)/(1.+csat(ll)/mtot)
                  fmtot= fmtot+dum
                  fmtot2= fmtot2+dum**2
              ENDDO
              fmtot= fmtot+minit   ! forecast total soa mass
              dum=   mtot-fmtot
              dum2=  ((fmtot-minit)/mtot)-mtot*fmtot2
              mtot=  mtot-dum/(1.-dum2)
              mtot=  AMAX1(mtot,minit)                           ! Limit MTOT to min possible in case of instability
              mtot=  AMIN1(mtot,PnGtotal+minit)                  ! Limit MTOT to max possible in case of instability
           ENDDO  loop_newt   ! LL iteration number loop

! RAR: Debugging:
IF (diagn_opt .AND. icall<200 .AND. kgrid<3) then
    WRITE(6,*) 'soa2_vbs: mtot after partitioning ',mtot
    WRITE(6,*) 'soa2_vbs: fmtot,fmtot2,dum,dum2 ',fmtot,fmtot2,dum,dum2
END IF

! Have total mass MTOT, get aerosol mass from semi-ideal partitioning equation
      DO ll=1,lsv_ocv
          caer(ll)= ctot(ll)* mtot/( mtot+csat(ll) )
          caer(ll)= AMAX1(caer(ll),conmin)
          cgas(ll)= ctot(ll) - caer(ll)
      ENDDO

! These species are in gas phase only, I do partitioning anyway just to keep consistency of the code
      cgas(ipoa4) = ctot(ipoa4)
      cgas(iasoa4)= ctot(iasoa4)
      cgas(ibsoa4)= ctot(ibsoa4)

! RAR: Debugging:
IF (diagn_opt .AND. icall<200 .AND. kgrid<3) then
    WRITE(6,*) 'soa2_vbs: cgas after part. ',cgas
    WRITE(6,*) 'soa2_vbs: caer after part. ',caer
END IF

! assigning values to CBLK array (gases), convert to ppm since it goes to chem
      cblk(lcell,vcvpoa1) = max(cgas(ipoa1),conmin)/convfac/mw_ocv(ipoa1)
      cblk(lcell,vcvpoa2) = max(cgas(ipoa2),conmin)/convfac/mw_ocv(ipoa2)
      cblk(lcell,vcvpoa3) = max(cgas(ipoa3),conmin)/convfac/mw_ocv(ipoa3)
      cblk(lcell,vcvpoa4) = max(cgas(ipoa4),conmin)/convfac/mw_ocv(ipoa4)

      cblk(lcell,vcvasoa1) = max(cgas(iasoa1),conmin)/convfac/mw_ocv(iasoa1)
      cblk(lcell,vcvasoa2) = max(cgas(iasoa2),conmin)/convfac/mw_ocv(iasoa2)
      cblk(lcell,vcvasoa3) = max(cgas(iasoa3),conmin)/convfac/mw_ocv(iasoa3)
      cblk(lcell,vcvasoa4) = max(cgas(iasoa4),conmin)/convfac/mw_ocv(iasoa4)

      cblk(lcell,vcvbsoa1) = max(cgas(ibsoa1),conmin)/convfac/mw_ocv(ibsoa1)
      cblk(lcell,vcvbsoa2) = max(cgas(ibsoa2),conmin)/convfac/mw_ocv(ibsoa2)
      cblk(lcell,vcvbsoa3) = max(cgas(ibsoa3),conmin)/convfac/mw_ocv(ibsoa3)
      cblk(lcell,vcvbsoa4) = max(cgas(ibsoa4),conmin)/convfac/mw_ocv(ibsoa4)

! RAR: Debugging:
IF (diagn_opt .AND. icall<200 .AND. kgrid<3) then
    WRITE(6,*) 'soa2_vbs: cblk ',cblk
END IF

! Condensational growth rates for POA, ASOA and BSOA species
      orgpoa1rat(lcell) = ( caer(ipoa1)-aold(ipoa1) )/dt
      orgpoa2rat(lcell) = ( caer(ipoa2)-aold(ipoa2) )/dt
      orgpoa3rat(lcell) = ( caer(ipoa3)-aold(ipoa3) )/dt

      organt1rat(lcell) = ( caer(iasoa1)-aold(iasoa1) )/dt    ! Semi-volatile ASOA in bins 1,2,3
      organt2rat(lcell) = ( caer(iasoa2)-aold(iasoa2) )/dt
      organt3rat(lcell) = ( caer(iasoa3)-aold(iasoa3) )/dt

      orgbio1rat(lcell) = ( caer(ibsoa1)-aold(ibsoa1) )/dt
      orgbio2rat(lcell) = ( caer(ibsoa2)-aold(ibsoa2) )/dt
      orgbio3rat(lcell) = ( caer(ibsoa3)-aold(ibsoa3) )/dt

! To convert from ppmv to ug/m3, divide by dt to get the production rate      
      organt0rat(lcell) = drog_in(pocv0)* convfac* mw_ocv(iasoa1) /dt    ! Non-volatile ASOA production in bin-0

! RAR: Debugging:
IF (diagn_opt .AND. icall<200 .AND. kgrid<3) then
    WRITE(6,*) 'soa2_vbs: orgpoa1rat,orgpoa2rat,orgpoa3rat ',orgpoa1rat,orgpoa2rat,orgpoa3rat
    WRITE(6,*) 'soa2_vbs: organt0rat,organt1rat,organt2rat,organt3rat ',organt0rat,organt1rat,organt2rat,organt3rat
    WRITE(6,*) 'soa2_vbs: orgbio1rat,orgbio2rat,orgbio3rat ',orgbio1rat,orgbio2rat,orgbio3rat

    icall=icall+1
END IF

  RETURN
END SUBROUTINE soa2_vbs
!------------------------------------------------------------------------------
!
! *** this routine calculates the dry deposition and sedimentation
!     velocities for the three modes. 
!     coded 1/23/97 by Dr. Francis S. Binkowski. Follows 
!     FSB's original method, i.e. uses Jon Pleim's expression for deposition
!     velocity but includes Marv Wesely's wstar contribution. 
!ia eliminated Stokes term for coarse mode deposition calcs.,
!ia see comments below
! RAR: the VDVG subroutine is removed, since it isn't used
!
SUBROUTINE VDVG_2( layer, blkta, blkdens,               &        ! layer -> k
                   ra, ustar, pblh, zntt, rmolm,  amu,  &
                   dgnuc, dgacc, dgcor, xlm,            &        ! geometric mean diameters
                   knnuc, knacc,kncor,                  &        ! knudsen numbers
                   pdensn, pdensa, pdensc,              &        ! average particle density
                   vsed, vdep, diagn_opt )                       ! sedimentation and deposition velocities

! *** calculate size-averaged particle dry deposition and 
!     size-averaged sedimentation velocities.
!     IMPLICIT NONE

LOGICAL :: diagn_opt
!      INTEGER blksize                 ! dimension of arrays
!      INTEGER NSPCSDA                 ! number of species in CBLK
!      INTEGER NUMCELLS                ! actual number of cells in arrays
      INTEGER LAYER                   ! number of layer
      INTEGER, PARAMETER :: iprnt = 0

!     REAL CBLK( blksize, NSPCSDA ) ! main array of variables
      real blkta( blksize )         ! air temperature [ k ]
      real blkdens(blksize)         ! air density  [ kg m^-3 ]
      real ra(blksize )             ! aerodynamic resistance [ s m**-1 ]
      real ustar( blksize )         ! surface friction velocity [ m s**-1 ]
      real pblh( blksize )          ! pbl height (m)
      real zntt( blksize )          ! surface roughness length (m)
      real rmolm( blksize )         ! inverse of monin-obukhov length (1/m)
      real amu( blksize )           ! atmospheric dynamic viscosity [ kg m**-1 s**-1 ]
      real xlm( blksize )           ! mean free path of dry air [ m ]
      real dgnuc( blksize )         ! nuclei mode mean diameter [ m ]
      real dgacc( blksize )         ! accumulation  
      real dgcor( blksize )         ! coarse mode
      real knnuc( blksize )         ! nuclei mode knudsen number 
      real knacc( blksize )         ! accumulation  
      real kncor( blksize )         ! coarse mode
      real pdensn( blksize )        ! average particle density in nuclei mode [ kg / m**3 ]
      real pdensa( blksize )        ! average particle density in accumulation mode [ kg / m**3 ]
      real pdensc( blksize )        ! average particle density in coarse mode [ kg / m**3 ]

! *** deposition and sedimentation velocities

      REAL VDEP( blksize, naspcsdep) ! sedimentation velocity [ m s**-1 ]
      REAL VSED( blksize, naspcssed) ! deposition  velocity [ m s**-1 ]

      INTEGER LCELL,N
      REAL DCONST1, DCONST2, DCONST3, DCONST3N, DCONST3A,DCONST3C
      REAL UTSCALE,CZH   ! scratch functions of USTAR and WSTAR.
      REAL NU            ! kinematic viscosity [ m**2 s**-1 ]
      REAL BHAT
      PARAMETER( BHAT =  1.246 ) ! Constant from Binkowski-Shankar approx to Cunningham slip correction.
      REAL COLCTR_BIGD,COLCTR_SMALD
      PARAMETER ( COLCTR_BIGD=2.E-3,COLCTR_SMALD=20.E-6 )  ! Collector diameters in Stokes number and Interception Efficiency (Needleleaf Forest)
      REAL SUM0, SUM3, DQ, KNQ, CUNQ, VSEDQ, SCQ, STQ, RSURFQ, vdplim
      REAL Eff_dif, Eff_imp, Eff_int, RBcor
      INTEGER ISTOPvd0,IdoWesCor
      PARAMETER (ISTOPvd0 = 0)  ! ISTOPvd0 = 1 means dont call VDVG, particle dep. velocities are set = 0; ISTOPvd0 = 0 means do depvel calcs.

      INTEGER, SAVE :: icall

      ! no Wesley deposition, otherwise EC is too low
      PARAMETER (IdoWesCor = 0) ! IdoWesCor = 1 means do Wesley (85) convective correction to PM dry dep velocities; 0 means don't do correction
      IF (ISTOPvd0.EQ.1)THEN
      RETURN
      ENDIF
! *** check layer value. 

         IF ( LAYER .EQ. 1 ) THEN ! calculate diffusitities and sedimentation velocities
                 
         DO LCELL = 1, NUMCELLS
            DCONST1 = BOLTZ * BLKTA(LCELL) /                                         &
                    ( THREEPI * AMU(LCELL) )
            DCONST2 = GRAV / ( 18.0 * AMU(LCELL) )
            DCONST3 =  USTAR(LCELL)/(9.*AMU(LCELL)*COLCTR_BIGD)
 
! *** now calculate the deposition velocities at layer 1

         NU = AMU(LCELL) / BLKDENS(LCELL) 

         UTSCALE =  1.
        IF (IdoWesCor.EQ.1)THEN
! Wesley (1985) Monin-Obukov dependence for convective conditions (SAM 10/08)
           IF(RMOLM(LCELL).LT.0.)THEN
                CZH = -1.*PBLH(LCELL)*RMOLM(LCELL)
                IF(CZH.GT.30.0)THEN
                  UTSCALE=0.45*CZH**0.6667
                ELSE
                  UTSCALE=1.+(-300.*RMOLM(LCELL))**0.6667
                ENDIF
           ENDIF
        ENDIF   ! end of (IdoWesCor.EQ.1) test

        UTSCALE = USTAR(LCELL)*UTSCALE


! RAR: Debugging:
IF (diagn_opt .AND. icall<100) then
    WRITE(6,*) 'VDVG_2: grav,ngausdv,xxlsga,ustar,utscale ', grav,ngausdv,xxlsga,ustar,utscale
    WRITE(6,*) 'VDVG_2: dconst2,pdensa,dgacc,amu ', dconst2,pdensa(lcell),dgacc(lcell),amu(lcell)
    icall=icall+1
end if

! *** nuclei mode
        SUM0=0.
        SUM3=0.
        DO N=1,NGAUSdv
         DQ=DGNUC(LCELL)*EXP(Y_GQ(N)*sqrt2*xxlsgn)  ! Diameter (m) at quadrature point
            KNQ=2.*XLM(LCELL)/DQ                    ! Knudsen number at quadrature point
            CUNQ=1.+KNQ*(1.257+.4*exp(-1.1/KNQ))    ! Cunningham correction factor; Pruppacher and Klett (1980) Eq (12-16)
            VSEDQ=PDENSN(LCELL)*DCONST2*CUNQ*DQ*DQ  ! Gravitational sedimentation velocity m/s
            SCQ=NU*DQ/DCONST1/CUNQ                  ! Schmidt number, Brownian diffusion parameter - Same as Binkowski and Shankar
            Eff_dif=SCQ**(-TWO3)                    ! Efficiency term for diffusion - Same as Binkowski and Shankar
            STQ=DCONST3*PDENSN(LCELL)*DQ**2         ! Stokes number, Peters and Eiden (1992)
            Eff_imp=(STQ/(0.8+STQ))**2              ! Efficiency term for impaction - Peters and Eiden (1992)
    !       Eff_int=0.3*DQ/(COLCTR_SMALD+DQ) ! Slinn (1982) Interception term, 0.3 prefac insures .2 cm/s at .2 micron diam.
            Eff_int=(0.00116+0.0061*ZNTT(LCELL))*DQ/1.414E-7 ! McKeen(2008) Intercptn trm, val of .00421 @ ustr=0.475, diam=.1414 micrn, stable, needleleaf evergreen
            RBcor=1.                                ! Rebound correction factor
            vdplim=UTSCALE*(Eff_dif+Eff_imp+Eff_int)*RBcor
    !       vdplim=.002*UTSCALE
            vdplim=min(vdplim,.02)
            RSURFQ=RA(LCELL)+1./vdplim
    !       RSURFQ=RA(LCELL)+1./(UTSCALE*(Eff_dif+Eff_imp+Eff_int)*RBcor) ! Total surface resistence
    !
!   limit this here to be consisten with the gocart routine, which bases this on Walcek et al. 1986
!
    !       RSURFQ=max(RSURFQ,50.)
            SUM0=SUM0+WGAUS(N)*(VSEDQ + 1./RSURFQ)  ! Quadrature sum for 0 moment
            SUM3=SUM3+WGAUS(N)*(VSEDQ + 1./RSURFQ)*DQ**3  ! Quadrature sum for 3rd moment
            ENDDO
            VDEP(LCELL, VDNNUC) = SUM0/sqrtpi  ! normalize 0 moment vdep quadrature sum to sqrt(pi) (and number =1 per unit volume)
            VDEP(LCELL, VDMNUC) = SUM3/(sqrtpi*EXP((1.5*sqrt2*xxlsgn)**2)*DGNUC(LCELL)**3) !normalize 3 moment quad. sum to sqrt(pi) and 3rd moment analytic sum

! *** accumulation mode

            SUM0=0.
            SUM3=0.
            DO N=1,NGAUSdv
            DQ=DGACC(LCELL)*EXP(Y_GQ(N)*sqrt2*xxlsga)  ! Diameter (m) at quadrature point
            KNQ=2.*XLM(LCELL)/DQ  ! Knudsen number at quadrature point
            CUNQ=1.+KNQ*(1.257+.4*exp(-1.1/KNQ))  ! Cunningham correction factor; Pruppacher and Klett (1980) Eq (12-16)
            VSEDQ=PDENSA(LCELL)*DCONST2*CUNQ*DQ*DQ  ! Gravitational sedimentation velocity m/s
            SCQ=NU*DQ/DCONST1/CUNQ  ! Schmidt number, Brownian diffusion parameter - Same as Binkowski and Shankar
            Eff_dif=SCQ**(-TWO3)    ! Efficiency term for diffusion - Same as Binkowski and Shankar
            STQ=DCONST3*PDENSA(LCELL)*DQ**2  ! Stokes number, Peters and Eiden (1992)
            Eff_imp=(STQ/(0.8+STQ))**2   ! Efficiency term for impaction - Peters and Eiden (1992)
    !       Eff_int=0.3*DQ/(COLCTR_SMALD+DQ) ! Slinn (1982) Interception term, 0.3 prefac insures .2 cm/s at .2 micron diam.
            Eff_int=(0.00116+0.0061*ZNTT(LCELL))*DQ/1.414E-7 ! McKeen(2008) Intercptn term, val of .00421 @ ustr=0.475, diam=.1414 micrn, stable, needleleaf evergreen
            RBcor=1. ! Rebound correction factor
            vdplim=UTSCALE*(Eff_dif+Eff_imp+Eff_int)*RBcor
            vdplim=min(vdplim,.02)
            RSURFQ=RA(LCELL)+1./vdplim
!       RSURFQ=RA(LCELL)+1./(UTSCALE*(Eff_dif+Eff_imp+Eff_int)*RBcor) ! Total surface resistence
!
!   limit this here to be consisten with the gocart routine, which bases this on Walcek et al. 1986
!
!       RSURFQ=max(RSURFQ,50.)
        SUM0=SUM0+WGAUS(N)*(VSEDQ + 1./RSURFQ)  ! Quadrature sum for 0 moment
        SUM3=SUM3+WGAUS(N)*(VSEDQ + 1./RSURFQ)*DQ**3  ! Quadrature sum for 3rd moment

        ENDDO
        VDEP(LCELL, VDNACC) = SUM0/sqrtpi  ! normalize 0 moment vdep quadrature sum to sqrt(pi) (and number =1 per unit volume)
        VDEP(LCELL, VDMACC) = SUM3/(sqrtpi*EXP((1.5*sqrt2*xxlsga)**2)*DGACC(LCELL)**3) !normalize 3 moment quad. sum to sqrt(pi) and 3rd moment analytic sum
        
! *** coarse mode 
        
        SUM0=0.
        SUM3=0.
        DO N=1,NGAUSdv
           DQ=DGCOR(LCELL)*EXP(Y_GQ(N)*sqrt2*xxlsgc)  ! Diameter (m) at quadrature point
           KNQ=2.*XLM(LCELL)/DQ  ! Knudsen number at quadrature point
           CUNQ=1.+KNQ*(1.257+.4*exp(-1.1/KNQ))  ! Cunningham correction factor; Pruppacher and Klett (1980) Eq (12-16)
           VSEDQ=PDENSC(LCELL)*DCONST2*CUNQ*DQ*DQ  ! Gravitational sedimentation velocity m/s
           SCQ=NU*DQ/DCONST1/CUNQ  ! Schmidt number, Brownian diffusion parameter - Same as Binkowski and Shankar
           Eff_dif=SCQ**(-TWO3)    ! Efficiency term for diffusion - Same as Binkowski and Shankar
           STQ=DCONST3*PDENSC(LCELL)*DQ**2  ! Stokes number, Peters and Eiden (1992)
           Eff_imp=(STQ/(0.8+STQ))**2   ! Efficiency term for impaction - Peters and Eiden (1992)
!          Eff_int=0.3*DQ/(COLCTR_SMALD+DQ) ! Slinn (1982) Interception term, 0.3 prefac insures .2 cm/s at .2 micron diam.
           Eff_int=(0.00116+0.0061*ZNTT(LCELL))*DQ/1.414E-7 ! McKeen(2008) Interception term, val of .00421 @ ustr=0.475, diam=.1414 micrn, stable, needleleaf evergreen
           EFF_int=min(1.,EFF_int)
           RBcor=exp(-2.0*(STQ**0.5)) ! Rebound correction factor used in Slinn (1982)
           vdplim=UTSCALE*(Eff_dif+Eff_imp+Eff_int)*RBcor
           vdplim=min(vdplim,.02)
           vdplim=max(vdplim,1e-35) !wig: add check since occasionally a lg particle causes overflow of rsurfq
           RSURFQ=RA(LCELL)+1./vdplim
!       RSURFQ=RA(LCELL)+1./(UTSCALE*(Eff_dif+Eff_imp+Eff_int)*RBcor) ! Total surface resistence
!
!   limit this here to be consisten with the gocart routine, which bases this on Walcek et al. 1986
!
!       RSURFQ=max(RSURFQ,50.)
           sum0=sum0+wgaus(n)*(vsedq + 1./rsurfq)  ! quadrature sum for 0 moment
           sum3=sum3+wgaus(n)*(vsedq + 1./rsurfq)*dq**3  ! quadrature sum for 3rd moment
        ENDDO
            vdep(lcell, vdncor) = sum0/sqrtpi  ! normalize 0 moment vdep quadrature sum to sqrt(pi) (and number =1 per unit volume)
            vdep(lcell, vdmcor) = sum3/(sqrtpi*exp((1.5*sqrt2*xxlsgc)**2)*dgcor(lcell)**3) !normalize 3 moment quad. sum to sqrt(pi) and 3rd moment analytic sum
        END DO
             
        ENDIF  ! ENDOF LAYER = 1 test
        
! *** Calculate gravitational sedimentation velocities for all layers - as in Binkowski and Shankar (1995)

         DO lcell = 1, numcells
         
            dconst2 = grav / ( 18.0 * amu(lcell) )
            dconst3n = dconst2 * pdensn(lcell) * dgnuc( lcell )**2
            dconst3a = dconst2 * pdensa(lcell) * dgacc( lcell )**2
            dconst3c = dconst2 * pdensc(lcell) * dgcor( lcell )**2
               
! *** nucleation mode number and mass sedimentation velociticies
            vsed( lcell, vsnnuc) = dconst3n                         &
               * ( esn16 + bhat * knnuc( lcell ) * esn04 )
            vsed( lcell, vsmnuc) = dconst3n                         &
               * (esn64 + bhat * knnuc( lcell ) * esn28 )
        
! *** accumulation mode number and mass sedimentation velociticies
            vsed( lcell, vsnacc) = dconst3a                          &
              * ( esa16 + bhat * knacc( lcell ) * esa04 )
            vsed( lcell, vsmacc) = dconst3a                          &
              * ( esa64 + bhat * knacc( lcell ) * esa28 )

! *** coarse mode number and mass sedimentation velociticies
            vsed( lcell, vsncor) = dconst3c                          &
              * ( esc16 + bhat * kncor( lcell ) * esc04 )
            vsed( lcell, vsmcor) = dconst3c                          &
              * ( esc64 + bhat * kncor( lcell ) * esc28 )
         END DO
END SUBROUTINE VDVG_2
!------------------------------------------------------------------------------

SUBROUTINE aerosols_soa2_init(chem,convfac,z_at_w,pm1,               &
                   pm2_5_dry,pm2_5_water,pm2_5_dry_ec,pm10,          &
!                   mgai,mgaj,caai,caaj,kai,kaj,                      &
                   chem_in_opt,aer_ic_opt,                           &
                   poa_vol_choice, poa_vol_params,                   &
                   init_crustal,do_isorropia,                        &
                   ids,ide, jds,jde, kds,kde,                        &
                   ims,ime, jms,jme, kms,kme,                        &
                   its,ite, jts,jte, kts,kte,diagn_opt               )

   USE module_configure, only: grid_config_rec_type

   IMPLICIT NONE
   INTEGER,      INTENT(IN   ) ::  chem_in_opt, aer_ic_opt, poa_vol_choice
   INTEGER,      INTENT(IN   ) ::  ids,ide, jds,jde, kds,kde,    &
                                   ims,ime, jms,jme, kms,kme,    &
                                   its,ite, jts,jte, kts,kte

 !  LOGICAL, INTENT(OUT) :: is_aerosol(num_chem)

   REAL, DIMENSION (5), INTENT(INOUT) :: poa_vol_params

   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme, num_chem ) ,     &
          INTENT(INOUT  ) ::  chem
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ) ,               &
          INTENT(INOUT  ) ::  pm2_5_dry,pm2_5_water,pm2_5_dry_ec,pm1,pm10 !, &
!                              mgai,mgaj,caai,caaj,kai,kaj
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ) ,               &
          INTENT(IN     ) ::  convfac
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ) ,               &
          INTENT(IN  ) ::   z_at_w
   LOGICAL, INTENT(IN ) :: init_crustal,do_isorropia
!   TYPE (grid_config_rec_type) , INTENT (in) :: config_flags

   integer i,j,k,l,ii,jj,kk,itype,ai_phase,isize
   real tempfac,zz
!  real,dimension(its:ite,kts:kte,jts:jte) :: convfac
   REAL splitfac
                        !between gas and aerosol phase
   REAL so4vaptoaer
!factor for splitting initial conc. of SO4
!3rd moment i-mode [3rd moment/m^3]
   REAL m3nuc
!3rd MOMENT j-mode [3rd moment/m^3]
   REAL m3acc
!    REAL ESN36
   REAL m3cor
   DATA splitfac/.98/
   DATA so4vaptoaer/.999/

! molecular weights for sulfuric acid and oganics [ kg/mol ] MKS
REAL, PARAMETER ::  mworg=250.0E-03, mwh2so4=98.07354E-3

INTEGER, SAVE :: icall

LOGICAL :: diagn_opt

        ai_phase = 1
        lptr_so4_aer(    :,:,:) = 1
        lptr_nh4_aer(    :,:,:) = 1
        lptr_no3_aer(    :,:,:) = 1
        lptr_na_aer(     :,:,:) = 1
        lptr_cl_aer(     :,:,:) = 1
        lptr_asoa0_aer(:,:,:) = 1
        lptr_asoa1_aer(:,:,:) = 1
        lptr_asoa2_aer(:,:,:) = 1
        lptr_asoa3_aer( :,:,:) = 1
!        lptr_asoa4_aer( :,:,:) = 1
        lptr_bsoa1_aer( :,:,:) = 1
        lptr_bsoa2_aer( :,:,:) = 1
        lptr_bsoa3_aer( :,:,:) = 1
!        lptr_bsoa4_aer( :,:,:) = 1
!        lptr_orgpa_aer(  :,:,:) = 1
        lptr_poa0_aer( :,:,:)   = 1
        lptr_poa1_aer( :,:,:)   = 1
        lptr_poa2_aer( :,:,:)   = 1
        lptr_poa3_aer( :,:,:)   = 1
        lptr_ec_aer(     :,:,:) = 1
        lptr_p25_aer(    :,:,:) = 1
        lptr_polp_aer(   :,:,:) = 1
        lptr_pols_aer(   :,:,:) = 1
        lptr_anth_aer(   :,:,:) = 1
        lptr_seas_aer(   :,:,:) = 1
        lptr_soil_aer(   :,:,:) = 1
        numptr_aer(      :,:,:) = 1



! Aitken mode
        itype = 1
        isize = 1
        ncomp_aer(itype) = 19
        numptr_aer(      isize,itype,ai_phase) = p_nu0
        lptr_so4_aer(    isize,itype,ai_phase) = p_so4ai
        lptr_nh4_aer(    isize,itype,ai_phase) = p_nh4ai
        lptr_no3_aer(    isize,itype,ai_phase) = p_no3ai
        lptr_na_aer(     isize,itype,ai_phase) = p_naai
        lptr_cl_aer(     isize,itype,ai_phase) = p_clai
        lptr_asoa0_aer(isize,itype,ai_phase)   = p_asoa0i
        lptr_asoa1_aer(isize,itype,ai_phase)   = p_asoa1i
        lptr_asoa2_aer(isize,itype,ai_phase)   = p_asoa2i
        lptr_asoa3_aer( isize,itype,ai_phase)  = p_asoa3i
!        lptr_asoa4_aer( isize,itype,ai_phase)  = p_asoa4i
        lptr_bsoa1_aer( isize,itype,ai_phase)  = p_bsoa1i
        lptr_bsoa2_aer( isize,itype,ai_phase)  = p_bsoa2i
        lptr_bsoa3_aer( isize,itype,ai_phase)  = p_bsoa3i
!        lptr_bsoa4_aer( isize,itype,ai_phase)  = p_bsoa4i
!        lptr_orgpa_aer(  isize,itype,ai_phase) = p_orgpai
        lptr_poa0_aer(   isize,itype,ai_phase) = p_poa0i
        lptr_poa1_aer(   isize,itype,ai_phase) = p_poa1i
        lptr_poa2_aer(   isize,itype,ai_phase) = p_poa2i
        lptr_poa3_aer(   isize,itype,ai_phase) = p_poa3i
        lptr_ec_aer(     isize,itype,ai_phase) = p_eci
        lptr_p25_aer(    isize,itype,ai_phase) = p_p25i
! Accumulation mode
        itype = 1
        isize = 2
        ncomp_aer(itype) = 19
        numptr_aer(      isize,itype,ai_phase) = p_ac0
        lptr_so4_aer(    isize,itype,ai_phase) = p_so4aj
        lptr_nh4_aer(    isize,itype,ai_phase) = p_nh4aj
        lptr_no3_aer(    isize,itype,ai_phase) = p_no3aj
        lptr_na_aer(     isize,itype,ai_phase) = p_naaj
        lptr_cl_aer(     isize,itype,ai_phase) = p_claj
        lptr_asoa0_aer(isize,itype,ai_phase)   = p_asoa0j
        lptr_asoa1_aer(isize,itype,ai_phase)   = p_asoa1j
        lptr_asoa2_aer(isize,itype,ai_phase)   = p_asoa2j
        lptr_asoa3_aer( isize,itype,ai_phase)  = p_asoa3j
!        lptr_asoa4_aer( isize,itype,ai_phase)  = p_asoa4j
        lptr_bsoa1_aer( isize,itype,ai_phase)  = p_bsoa1j
        lptr_bsoa2_aer( isize,itype,ai_phase)  = p_bsoa2j
        lptr_bsoa3_aer( isize,itype,ai_phase)  = p_bsoa3j
!        lptr_bsoa4_aer( isize,itype,ai_phase)  = p_bsoa4j
!        lptr_orgpa_aer(  isize,itype,ai_phase) = p_orgpaj
        lptr_poa0_aer(   isize,itype,ai_phase) = p_poa0j
        lptr_poa1_aer(   isize,itype,ai_phase) = p_poa1j
        lptr_poa2_aer(   isize,itype,ai_phase) = p_poa2j
        lptr_poa3_aer(   isize,itype,ai_phase) = p_poa3j
        lptr_ec_aer(     isize,itype,ai_phase) = p_ecj
        lptr_p25_aer(    isize,itype,ai_phase) = p_p25j
        lptr_pols_aer(   isize,itype,ai_phase) = p_pols
! coarse mode
        itype = 2
        isize = 1
        ncomp_aer(itype) = 4
        numptr_aer(      isize,itype,ai_phase) = p_corn
        lptr_anth_aer(   isize,itype,ai_phase) = p_antha
        lptr_seas_aer(   isize,itype,ai_phase) = p_seas
        lptr_soil_aer(   isize,itype,ai_phase) = p_soila
        lptr_polp_aer(   isize,itype,ai_phase) = p_polp
! *** Compute these once and they will all  be saved in COMMON
   xxlsgn = log(sginin)
   xxlsga = log(sginia)
   xxlsgc = log(sginic)

   l2sginin = xxlsgn**2
   l2sginia = xxlsga**2
   l2sginic = xxlsgc**2

   en1 = exp(0.125*l2sginin)
   ea1 = exp(0.125*l2sginia)
   ec1 = exp(0.125*l2sginic)

   esn04 = en1**4
   esa04 = ea1**4
   esc04 = ec1**4

   esn05 = esn04*en1
   esa05 = esa04*ea1

   esn08 = esn04*esn04
   esa08 = esa04*esa04
   esc08 = esc04*esc04

   esn09 = esn04*esn05
   esa09 = esa04*esa05

   esn12 = esn04*esn04*esn04
   esa12 = esa04*esa04*esa04
   esc12 = esc04*esc04*esc04

   esn16 = esn08*esn08
   esa16 = esa08*esa08
   esc16 = esc08*esc08

   esn20 = esn16*esn04
   esa20 = esa16*esa04
   esc20 = esc16*esc04

   esn24 = esn12*esn12
   esa24 = esa12*esa12
   esc24 = esc12*esc12

   esn25 = esn16*esn09
   esa25 = esa16*esa09

   esn28 = esn20*esn08
   esa28 = esa20*esa08
   esc28 = esc20*esc08


   esn32 = esn16*esn16
   esa32 = esa16*esa16
   esc32 = esc16*esc16

   esn36 = esn16*esn20
   esa36 = esa16*esa20
   esc36 = esc16*esc20

   esn49 = esn25*esn20*esn04
   esa49 = esa25*esa20*esa04

   esn52 = esn16*esn36
   esa52 = esa16*esa36

   esn64 = esn32*esn32
   esa64 = esa32*esa32
   esc64 = esc32*esc32

   esn100 = esn36*esn64

   esnm20 = 1.0/esn20
   esam20 = 1.0/esa20
   escm20 = 1.0/esc20

   esnm32 = 1.0/esn32
   esam32 = 1.0/esa32
   escm32 = 1.0/esc32

   xxm3 = 3.0*xxlsgn/ sqrt2
! factor used in error function cal
        nummin_i = facatkn_min*so4fac*aeroconcmin/(dginin**3*esn36)

        nummin_j = facacc_min*so4fac*aeroconcmin/(dginia**3*esa36)

        nummin_c = anthfac*aeroconcmin/(dginic**3*esc36)

! *** Note, DGVEM_I, DGVEM_J, DGVEM_C are for the mass (volume)
!     size distribution , then
!        vol = (p/6) * density * num * (dgemv_xx**3) *
!                            exp(- 4.5 * log( sgem_xx)**2 ) )
!        note minus sign!!

        factnumn = exp(4.5*log(sgem_i)**2)/dgvem_i**3
        factnuma = exp(4.5*log(sgem_j)**2)/dgvem_j**3
        factnumc = exp(4.5*log(sgem_c)**2)/dgvem_c**3

        ccofm = alphsulf*sqrt(pirs*rgasuniv/(2.0*mwh2so4))
        ccofm_org = alphaorg*sqrt(pirs*rgasuniv/(2.0*mworg))

! RAR: Debugging:
IF (diagn_opt .AND. icall<200) then
    WRITE(6,*) 'aerosols_soa2_init: ccofm, ccofm_org', ccofm, ccofm_org
end if

!        mwso4=96.03
! initialize pointers used by aerosol-cloud-interaction routines
!        call aerosols_soa_vbs_init_aercld_ptrs( &
!           num_chem, is_aerosol, config_flags )
        pm1(its:ite, kts:kte-1, jts:jte)    = 0.
        pm10(its:ite, kts:kte-1, jts:jte)    = 0.
        pm2_5_dry(its:ite, kts:kte-1, jts:jte)    = 0.
        pm2_5_water(its:ite, kts:kte-1, jts:jte)  = 0.
        pm2_5_dry_ec(its:ite, kts:kte-1, jts:jte) = 0.

!        if ( init_crustal .and. do_isorropia ) then
!        do j=jts,jte
!        do k=kts,kte-1
!        do i=its,ite
!          caaj(i,k,j)  = 10.E-05
!          caai(i,k,j)  = 10.E-05
!          kaj(i,k,j)   = 10.E-05
!          kai(i,k,j)   = 10.E-05
!          mgaj(i,k,j)  = 10.E-05
!          mgai(i,k,j)  = 10.E-05
!        enddo
!        enddo
!        enddo
!        endif

!SAM 10/08 Add in Gaussian quadrature points and weights - Use 7 points = NGAUSdv

        Y_GQ(1)=-2.651961356835233
        WGAUS(1)=0.0009717812450995
        Y_GQ(2)=-1.673551628767471
        WGAUS(2)=0.05451558281913
        Y_GQ(3)=-0.816287882858965
        WGAUS(3)=0.4256072526101
        Y_GQ(4)=-0.0
        WGAUS(4)=0.8102646175568
        Y_GQ(5)=0.816287882858965
        WGAUS(5)=WGAUS(3)
        Y_GQ(6)=1.673551628767471
        WGAUS(6)=WGAUS(2)
        Y_GQ(7)=2.651961356835233
        WGAUS(7)=WGAUS(1)

! RAR: initialization
        sel_pars: SELECT CASE(poa_vol_choice)          ! vol_zero  / 1.00, 0.00, 0.00, 0.00, 0.00 /  ! sum is 1.0
        CASE (0)        ! Zero volatility              ! vol_svoc  / 0.09, 0.09, 0.14, 0.18, 0.50 /  ! sum is 1.0
          poa_vol_params= vol_zero                     ! vol_ivoc  / 0.09, 0.09, 0.14, 0.18, 2.00 /  ! sum is 2.5
                                                       ! vol_high  / 0.27, 0.27, 0.42, 0.54, 6.00 /  ! sum is 7.5
        CASE (1)        ! Base case,                   ! vol_wood  / 0.00, 0.10, 0.14, 0.33, 0.43 /  ! sum is 1.0
          poa_vol_params= vol_svoc                     ! vol_test1 / 0.50, 0.50, 0.00, 0.00, 0.00 /  ! sum is 1.0
                                                       ! vol_test2 / 0.18, 0.18, 0.28, 0.36, 0.00 /  ! sum is 1.0
        CASE (2)        ! IVOCs, =1.5xPOA
          poa_vol_params= vol_ivoc

        CASE (3)        ! High volatility, =7.5xPOA
          poa_vol_params= vol_high

        CASE (4)        ! Wood burning, =1.0
          poa_vol_params= vol_wood

        CASE (5)        ! Test case, evenly split between bins 0 and 1
          poa_vol_params= vol_test1

        CASE (6)        ! Test case, low volatility
          poa_vol_params= vol_test2

        CASE DEFAULT
          CALL wrf_error_fatal ("NO VOLATILITY TABLE IS INCLUDED FOR THIS POA_EMISS_VOL OPTION!")
        END SELECT sel_pars



!  IF USING OLD SIMULATION, DO NOT REINITIALIZE!
        if(chem_in_opt == 1 ) return

        do l=p_so4aj,num_chem
           chem(ims:ime,kms:kme,jms:jme,l)=epsilc
        enddo
        chem(ims:ime,kms:kme,jms:jme,p_nu0)=1.e8
        chem(ims:ime,kms:kme,jms:jme,p_ac0)=1.e8
        do j=jts,jte
           jj=min(jde-1,j)
        do k=kts,kte-1
           kk=min(kde-1,k)
        do i=its,ite
           ii=min(ide-1,i)

!Option for alternate ic's
        if( aer_ic_opt == AER_IC_DEFAULT ) then
          chem(i,k,j,p_so4aj)= chem(ii,kk,jj,p_sulf)*convfac(ii,kk,jj)*mwso4*splitfac*so4vaptoaer
          chem(i,k,j,p_so4ai)= chem(ii,kk,jj,p_sulf)*convfac(ii,kk,jj)*mwso4*(1.-splitfac)*so4vaptoaer
          chem(i,k,j,p_sulf) = chem(ii,kk,jj,p_sulf)*(1.-so4vaptoaer)
          chem(i,k,j,p_nh4aj) = 10.E-05
          chem(i,k,j,p_nh4ai) = 10.E-05
          chem(i,k,j,p_no3aj) = 10.E-05
          chem(i,k,j,p_no3ai) = 10.E-05
          chem(i,k,j,p_naaj)  = 10.E-05
          chem(i,k,j,p_naai)  = 10.E-05
          chem(i,k,j,p_claj)  = 10.E-05
          chem(i,k,j,p_clai)  = 10.E-05

          chem(i,k,j,p_p25i)  = 10.E-05
          chem(i,k,j,p_p25j)  = 10.E-05
          chem(i,k,j,p_seas)  = 10.E-05
          chem(i,k,j,p_antha) = 10.E-05
          chem(i,k,j,p_soila) = 10.E-05
          chem(i,k,j,p_eci)    = 10.E-05
          chem(i,k,j,p_ecj)    = 10.E-05   
 
!          if (do_isorropia ) then
!          caaj(i,k,j)  = 10.E-05
!          caai(i,k,j)  = 10.E-05
!          kaj(i,k,j)   = 10.E-05
!          kai(i,k,j)   = 10.E-05
!          mgaj(i,k,j)  = 10.E-05
!          mgai(i,k,j)  = 10.E-05
!          endif
!        elseif( aer_ic_opt == AER_IC_PNNL ) then
!           zz = (z_at_w(ii,k,jj)+z_at_w(ii,k+1,jj))*0.5
!           call soa_vbs_init_aer_ic_pnnl(   &
!                chem, zz, i,k,j, ims,ime,jms,jme,kms,kme )
        else
           call wrf_error_fatal("aerosols_soa_vbs_init: unable to parse aer_ic_opt" )
        end if

      enddo
      enddo
      enddo

        do j=jts,jte
        do k=kts,kte
        do i=its,ite
        !... i-mode                     RAR: in future use different anthfac for P25 and EC
              m3nuc= so4fac*chem(i,k,j,p_so4ai) + nh4fac*chem(i,k,j,p_nh4ai) + no3fac*chem(i,k,j,p_no3ai) + &
                 orgfac*( chem(i,k,j,p_poa0i)  + chem(i,k,j,p_poa1i)  + chem(i,k,j,p_poa2i) + chem(i,k,j,p_poa3i) ) + &
                 orgfac*( chem(i,k,j,p_asoa0i) + chem(i,k,j,p_asoa1i) + chem(i,k,j,p_asoa2i) + chem(i,k,j,p_asoa3i) ) + &
                 orgfac*( chem(i,k,j,p_bsoa1i) + chem(i,k,j,p_bsoa2i) + chem(i,k,j,p_bsoa3i) ) + &
                 anthfac*chem(i,k,j,p_p25i)    + anthfac*chem(i,k,j,p_eci) + &
                 nafac*chem(i,k,j,p_naai)      + clfac*chem(i,k,j,p_clai)
!         if ( do_isorropia ) then
!             m3nuc=m3nuc + mgfac*mgai(i,k,j) + cafac*caai(i,k,j) + kfac*kai(i,k,j)
!         endif
                 
        !... j-mode
              m3acc= so4fac*chem(i,k,j,p_so4aj) + nh4fac*chem(i,k,j,p_nh4aj) + no3fac*chem(i,k,j,p_no3aj) + &
                 orgfac*( chem(i,k,j,p_poa0j)  + chem(i,k,j,p_poa1j)  + chem(i,k,j,p_poa2j)  + chem(i,k,j,p_poa3j) ) + &
                 orgfac*( chem(i,k,j,p_asoa0j) + chem(i,k,j,p_asoa1j) + chem(i,k,j,p_asoa2j) + chem(i,k,j,p_asoa3j) ) + &
                 orgfac*( chem(i,k,j,p_bsoa1j) + chem(i,k,j,p_bsoa2j) + chem(i,k,j,p_bsoa3j) ) + &
                 anthfac*chem(i,k,j,p_p25j)    + anthfac*chem(i,k,j,p_ecj) + &
                 nafac*chem(i,k,j,p_naaj)      + clfac*chem(i,k,j,p_claj) + polsfac*chem(i,k,j,p_pols)
!          if ( do_isorropia ) then
!             m3acc= m3acc + mgfac*mgaj(i,k,j) + cafac*caaj(i,k,j) + kfac*kaj(i,k,j)
!          endif

        ! RAR: Debugging:
        IF (diagn_opt .AND. icall<200) then
            WRITE(6,*) 'aerosols_soa2_init: orgfac,anthfac,so4fac,nh4fac,no3fac ', orgfac,anthfac,so4fac,nh4fac,no3fac
            WRITE(6,*) 'aerosols_soa2_init: m3nuc,m3acc ', m3nuc,m3acc
            icall=icall+1
        ENDIF
        !...c-mode
              m3cor = soilfac*chem(i,k,j,p_soila) + seasfac*chem(i,k,j,p_seas) + &
                anthfac*chem(i,k,j,p_antha)+polpfac*chem(i,k,j,p_polp)

        !...NOW CALCULATE INITIAL NUMBER CONCENTRATION
              chem(i,k,j,p_nu0) = m3nuc/((dginin**3)*esn36)
              chem(i,k,j,p_ac0) = m3acc/((dginia**3)*esa36)
              chem(i,k,j,p_corn) = m3cor/((dginic**3)*esc36)

        enddo
        enddo
        enddo
    return
    END SUBROUTINE aerosols_soa2_init
!------------------------------------------------------------------------

!SUBROUTINE soa_vbs_addemiss( id, dtstep, u10, v10, alt, dz8w, xland, chem,                      &
!                             ebu,                                                               &
!                             slai,ust,smois,ivgtyp,isltyp,                                      &
!                             emis_ant,dust_emiss_active,                                        &
!                             seasalt_emiss_active,kemit,biom,num_soil_layers,emissopt,          &
!                             dust_opt,ktau,p8w,u_phy,v_phy,rho_phy,g,dx,erod,                   &

SUBROUTINE soa2_addemiss( id,seasalt_emiss_active, dust_opt, dtstep,  &
                          alt, dz8w, chem,scalar,                                  &
                          mgai,mgaj,caai,caaj,kai,kaj,                             &
                          chem_is_wfa,chem_is_ifa,                                 &
                          use_initial_IFA_as_dust,                                 &
                          use_initial_WFA_as_ss,                                   &
                          do_isorropia,                                            &
                          emit_crustal,                                            &
                          u10,v10,xland,xlat,xlong,                                &  ! For seasalt
                          ktau,num_soil_layers,u_phy,v_phy,rho_phy,smois,          & ! For Dust
                          p8w,erod,ivgtyp,isltyp,dx,g,                             & ! For Dust
                          ebu,emis_ant, kemit, emissopt,                           &
                          biom,poa_vol_params,                                     &
                          pollen_emiss_active,pollen_opt,                          &
                          polp_emis,                                               &
                          rainc, rainnc,rain_old,relhum,                           &
                          swdown,swdnbc,                                           &
                          ids,ide, jds,jde, kds,kde,                               &
                          ims,ime, jms,jme, kms,kme,                               &
                          its,ite, jts,jte, kts,kte, diagn_opt                     )
!
! Routine to apply aerosol emissions for MADE/SOA_VBS...
! William.Gustafson@pnl.gov; 3-May-2007
! Modified by steven.peckham@noaa.gov;   8-Jan-2008
! Modified by R.Ahmadov, 05/23/2013
!
  USE module_state_description, only:  num_chem,p_qnifa,p_qnwfa,num_scalar

  LOGICAL :: diagn_opt
 
  INTEGER, INTENT(IN   )   ::    kemit, emissopt, biom,                 &
                                 ids,ide, jds,jde, kds,kde,             &
                                 ims,ime, jms,jme, kms,kme,             &
                                 its,ite, jts,jte, kts,kte,             &
                                 seasalt_emiss_active,dust_opt,         &
                                 ktau,num_soil_layers,id
  INTEGER, INTENT(IN )      ::   pollen_emiss_active,pollen_opt

  REAL, INTENT(IN   ) ::    dtstep

  REAL,  DIMENSION (5), INTENT(IN)   ::      poa_vol_params

! trace species mixing ratios (aerosol mass = ug/kg-air; number = #/kg-air)
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),               &
       INTENT(INOUT ) ::   chem

  REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_scalar ),               &
       INTENT(INOUT ) ::   scalar

  REAL, DIMENSION( ims:ime, kms:kme, jms:jme),               &
       INTENT(INOUT ) ::   mgai,mgaj,caai,caaj,kai,kaj
!
! aerosol emissions arrays ((ug/m3)*m/s)
!
  REAL, DIMENSION( ims:ime, kms:kemit, jms:jme,num_emis_ant ),         &
        INTENT(IN ) ::    emis_ant

! biomass burning aerosol emissions arrays ((ug/m3)*m/s)
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme,num_ebu ),              &
         INTENT(IN    ) ::    ebu

! 1/(dry air density) and layer thickness (m)
  REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                      &
         INTENT(IN ) ::    alt, dz8w,relhum
   
  REAL, DIMENSION( ims:ime, jms:jme), INTENT(OUT ) :: polp_emis

  ! add for gocart dust
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                        &
        INTENT(IN    ) :: p8w,u_phy,v_phy,rho_phy
  REAL, INTENT(IN    ) :: dx, g
  REAL, DIMENSION( ims:ime, jms:jme, 3 ),                              &
         INTENT(IN    ) :: erod

  REAL,  DIMENSION( ims:ime , jms:jme ),                                &
       INTENT(IN   ) ::                                                 &
       u10, v10, xland,xlat,xlong,rainc,rainnc,rain_old,swdown,swdnbc
  INTEGER,  DIMENSION( ims:ime , jms:jme ),                             &
       INTENT(IN   ) ::   ivgtyp, isltyp
  REAL, DIMENSION( ims:ime, num_soil_layers, jms:jme ),    &
       INTENT(INOUT) ::   smois
  LOGICAL, INTENT(IN) :: chem_is_wfa, chem_is_ifa,use_initial_IFA_as_dust, &
                         use_initial_WFA_as_ss,do_isorropia,emit_crustal

! Local variables...
  real, dimension(its:ite,kts:kte,jts:jte) :: factor

  integer :: i,j

  INTEGER, SAVE :: icall

    REAL, PARAMETER :: dens_nh4so4           = 1.8e+3
    REAL, PARAMETER :: dens_ss               = 2.2e+3
    REAL, PARAMETER :: dens_dust             = 2.6e+3
    REAL, PARAMETER :: dens_oc               = 1.4e+3
    REAL, PARAMETER :: dens_bc               = 1.4e+3
    REAL, PARAMETER :: dens_pols             = 1425.
    REAL, PARAMETER :: frac_seas_wfa         = 0.025
    REAL, PARAMETER :: frac_bc_wfa           = 0.1

    REAL, PARAMETER :: mean_diameter1        = 4e-8
    REAL, PARAMETER :: mean_diameter_ss      = 4e-6
    REAL, PARAMETER :: mean_diameter_soa     = 4e-8
    REAL, PARAMETER :: mean_diameter_poa     = 4e-8
    REAL, PARAMETER :: mean_diameter_nh4so4  = 4e-8
    REAL, PARAMETER :: sigma1                = 1.8
    REAL, PARAMETER :: sigma_ss              = 1.8
    REAL, PARAMETER :: sigma_soa             = 1.8
    REAL, PARAMETER :: sigma_poa             = 1.8
    REAL, PARAMETER :: sigma_nh4so4          = 1.8

    REAL, PARAMETER :: mean_diameter_dust    = 4e-7
    REAL, PARAMETER :: sigma_dust            = 1.8
    REAL, PARAMETER :: frac_dust_ifa         = 0.025
    
    REAL, PARAMETER :: pirs = 3.14159265358979324

    REAL, PARAMETER :: kappa_oc          = 0.2
    REAL, PARAMETER :: kappa_bc          = 0.2
    REAL, PARAMETER :: kappa_bsoa        = 0.2   ! 0.44, Müller et al., Scientific Reports, 2017
    REAL, PARAMETER :: kappa_asoa        = 0.2   ! ?
    REAL, PARAMETER :: kappa_poa         = 0.2
    REAL, PARAMETER :: kappa_nh4so4      = 0.2   ! 0.6 ?
    REAL, PARAMETER :: kappa_ss          = 0.2
    REAL, PARAMETER :: kappa_dust        = 0.04

    REAL, PARAMETER :: fact_ifa = 1.e-9*6.0/pirs*exp(4.5*log(sigma_dust)**2)/mean_diameter_dust**3
    REAL, PARAMETER :: fact_wfa = 1.e-9*6.0/pirs*exp(4.5*log(sigma1)**2)/mean_diameter1**3
    REAL, PARAMETER :: fact_wfa_ss = 1.e-9*6.0/pirs*exp(4.5*log(sigma1)**2)/mean_diameter_ss**3
!
! Get the emissions unit conversion factor including the time step.
! Changes emissions from [ug/m3 m/s] to [ug/kg_dryair/timestep]
!
  factor(its:ite,kts:kte,jts:jte) = alt(its:ite,kts:kte,jts:jte)*dtstep/ &
                  dz8w(its:ite,kts:kte,jts:jte)
!
! Increment the aerosol numbers...
!    if(emissopt  .lt. 5 )then
    if ( emissopt == 3 .or. emissopt == 23 ) then      ! For my runs emiss_opt=3
!
! Aitken mode first...
! RAR: For nu0 and ac0 I use only non-volatile part of the e_org* emissions!

! RAR: Debugging:
IF (diagn_opt .AND. icall<10) then
    WRITE(6,*) 'soa2_addemiss: its,ite,jts,jte ', its,ite,jts,jte
    WRITE(6,*) 'soa2_addemiss: factor(its,1,jts) ', factor(its,1,jts)
    WRITE(6,*) 'soa2_addemiss: anthfac,orgfac ', anthfac,orgfac
    WRITE(6,*) 'soa2_addemiss: poa_vol_params ', poa_vol_params
    WRITE(6,*) 'soa2_addemiss: factnumn,factnuma ', factnumn,factnuma
end if

! Nucleation mode:
  chem(its:ite,kts:kemit,jts:jte,p_nu0) =                         &
       chem(its:ite,kts:kemit,jts:jte,p_nu0) +                    &
       factor(its:ite,kts:kemit,jts:jte)*factnumn*(               &
       anthfac*( emis_ant(its:ite,kts:kemit,jts:jte,p_e_pm25i) +  &
                 emis_ant(its:ite,kts:kemit,jts:jte,p_e_eci) ) +  &
         orgfac* emis_ant(its:ite,kts:kemit,jts:jte,p_e_orgi)  +  &
         no3fac* emis_ant(its:ite,kts:kemit,jts:jte,p_e_no3i)  +  &
         so4fac* emis_ant(its:ite,kts:kemit,jts:jte,p_e_so4i)     ) 

! Accumulation mode:
  chem(its:ite,kts:kemit,jts:jte,p_ac0) =                         &
       chem(its:ite,kts:kemit,jts:jte,p_ac0) +                    &
       factor(its:ite,kts:kemit,jts:jte)*factnuma*(               &
       anthfac*( emis_ant(its:ite,kts:kemit,jts:jte,p_e_pm25j) +  &
                 emis_ant(its:ite,kts:kemit,jts:jte,p_e_ecj) ) +  &
       orgfac*   emis_ant(its:ite,kts:kemit,jts:jte,p_e_orgj)  +  &
       no3fac*   emis_ant(its:ite,kts:kemit,jts:jte,p_e_no3j) +  &
       so4fac*   emis_ant(its:ite,kts:kemit,jts:jte,p_e_so4j)     )

! RAR: Debugging:
if (diagn_opt .AND. icall<100) then
    WRITE(6,*) 'soa2_addemiss: its,jts,ite,jte ',its,jts,ite,jte
    WRITE(6,*) 'soa2_addemiss: chem(its,1,jts,p_nu0)',chem(its,1,jts,p_nu0)
    WRITE(6,*) 'soa2_addemiss: chem(its,1,jts,p_ac0)',chem(its,1,jts,p_ac0)
    icall=icall+1
end if

! And now the coarse mode...

  chem(its:ite,kts:kemit,jts:jte,p_corn) =                       &
       chem(its:ite,kts:kemit,jts:jte,p_corn) +                  &
       factor(its:ite,kts:kemit,jts:jte)*factnumc*anthfac*                           &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_pm_10)
!
! Increment the aerosol masses...
!
  chem(its:ite,kts:kemit,jts:jte,p_antha) =                      &
       chem(its:ite,kts:kemit,jts:jte,p_antha) +                 &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_pm_10)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_p25j) =                       &
       chem(its:ite,kts:kemit,jts:jte,p_p25j) +                  &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_pm25j)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_p25i) =                       &
       chem(its:ite,kts:kemit,jts:jte,p_p25i) +                  &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_pm25i)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_ecj) =                        &
       chem(its:ite,kts:kemit,jts:jte,p_ecj) +                   &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_ecj)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_eci) =                        &
       chem(its:ite,kts:kemit,jts:jte,p_eci) +                   &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_eci)*factor(its:ite,kts:kemit,jts:jte)

!  chem(its:ite,kts:kemit,jts:jte,p_naaj) =                        &
!       chem(its:ite,kts:kemit,jts:jte,p_naaj) +                   &
!       emis_ant(its:ite,kts:kemit,jts:jte,p_e_naaj)*factor(its:ite,kts:kemit,jts:jte)
!  chem(its:ite,kts:kemit,jts:jte,p_naai) =                        &
!       chem(its:ite,kts:kemit,jts:jte,p_naai) +                   &
!       emis_ant(its:ite,kts:kemit,jts:jte,p_e_naai)*factor(its:ite,kts:kemit,jts:jte)

!  chem(its:ite,kts:kemit,jts:jte,p_orgpaj) =                     &
!       chem(its:ite,kts:kemit,jts:jte,p_orgpaj) +                &
!       emis_ant(its:ite,kts:kemit,jts:jte,p_e_orgj)*factor(its:ite,kts:kemit,jts:jte)
!
!  chem(its:ite,kts:kemit,jts:jte,p_orgpai) =                     &
!       chem(its:ite,kts:kemit,jts:jte,p_orgpai) +                &
!       emis_ant(its:ite,kts:kemit,jts:jte,p_e_orgi)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_so4aj) =                      &
       chem(its:ite,kts:kemit,jts:jte,p_so4aj) +                 &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_so4j)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_so4ai) =                      &
       chem(its:ite,kts:kemit,jts:jte,p_so4ai) +                 &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_so4i)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_no3aj) =                      &
       chem(its:ite,kts:kemit,jts:jte,p_no3aj) +                 &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_no3j)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_no3ai) =                      &
       chem(its:ite,kts:kemit,jts:jte,p_no3ai) +                 &
       emis_ant(its:ite,kts:kemit,jts:jte,p_e_no3i)*factor(its:ite,kts:kemit,jts:jte)

!  chem(its:ite,kts:kemit,jts:jte,p_orgpaj) =                     &
!       chem(its:ite,kts:kemit,jts:jte,p_orgpaj) +                &
!       emis_ant(its:ite,kts:kemit,jts:jte,p_e_orgj)*factor(its:ite,kts:kemit,jts:jte)
!
!  chem(its:ite,kts:kemit,jts:jte,p_orgpai) =                     &
!       chem(its:ite,kts:kemit,jts:jte,p_orgpai) +                &
!       emis_ant(its:ite,kts:kemit,jts:jte,p_e_orgi)*factor(its:ite,kts:kemit,jts:jte)
!
!-------------------------------------------------------------------------------------------
! Adding emissions to the semi-volatile POA species:
!! We assume that all the POA emissions are in ug/m3

! Volatility bin 0:     (0.1 ug/m3 at 298K)
  chem(its:ite,kts:kemit,jts:jte,p_poa0j) =                     &
    chem(its:ite,kts:kemit,jts:jte,p_poa0j) +                &
    poa_vol_params(1)* emis_ant(its:ite,kts:kemit,jts:jte,p_e_orgj)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_poa0i) =                     &
    chem(its:ite,kts:kemit,jts:jte,p_poa0i) +                &
    poa_vol_params(1)* emis_ant(its:ite,kts:kemit,jts:jte,p_e_orgi)*factor(its:ite,kts:kemit,jts:jte)

! Volatility bin 1:     (1.0 ug/m3 at 298K)
  chem(its:ite,kts:kemit,jts:jte,p_poa1j) =                     &
    chem(its:ite,kts:kemit,jts:jte,p_poa1j) +                &
    poa_vol_params(2)* emis_ant(its:ite,kts:kemit,jts:jte,p_e_orgj)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_poa1i) =                     &
    chem(its:ite,kts:kemit,jts:jte,p_poa1i) +                &
    poa_vol_params(2)* emis_ant(its:ite,kts:kemit,jts:jte,p_e_orgi)*factor(its:ite,kts:kemit,jts:jte)

! Volatility bin 2:     (10.0 ug/m3 at 298K)
  chem(its:ite,kts:kemit,jts:jte,p_poa2j) =                     &
    chem(its:ite,kts:kemit,jts:jte,p_poa2j) +                &
    poa_vol_params(3)* emis_ant(its:ite,kts:kemit,jts:jte,p_e_orgj)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_poa2i) =                     &
    chem(its:ite,kts:kemit,jts:jte,p_poa2i) +                &
    poa_vol_params(3)* emis_ant(its:ite,kts:kemit,jts:jte,p_e_orgi)*factor(its:ite,kts:kemit,jts:jte)

! Volatility bin 3:     (100.0 ug/m3 at 298K)
  chem(its:ite,kts:kemit,jts:jte,p_poa3j) =                     &
    chem(its:ite,kts:kemit,jts:jte,p_poa3j) +                &
    poa_vol_params(4)* emis_ant(its:ite,kts:kemit,jts:jte,p_e_orgj)*factor(its:ite,kts:kemit,jts:jte)

  chem(its:ite,kts:kemit,jts:jte,p_poa3i) =                     &
    chem(its:ite,kts:kemit,jts:jte,p_poa3i) +                &
    poa_vol_params(4)* emis_ant(its:ite,kts:kemit,jts:jte,p_e_orgi)*factor(its:ite,kts:kemit,jts:jte)

! Volatility bin 4:     (1000.0 ug/m3 at 298K), this is gas species only!
! Thia bin contains gas phase OCVs only, therefore the emissions are added in add_anthropogenics
!-------------------------------------------------------------------------------------------

!  if(biom >= 1 )then
!!
!! Aitken mode first...
!
!  chem(its:ite,kts:kte,jts:jte,p_nu0) =                        &
!       chem(its:ite,kts:kte,jts:jte,p_nu0) +                   &
!       factor(its:ite,kts:kte,jts:jte)*factnumn*(              &
!       anthfac*( .25*ebu(its:ite,kts:kte,jts:jte,p_ebu_pm25) +       &
!              .25*ebu(its:ite,kts:kte,jts:jte,p_ebu_bc) ) +          &
!       orgfac*.25*ebu(its:ite,kts:kte,jts:jte,p_ebu_oc) )
!
!! Accumulation mode next...
!
!  chem(its:ite,kts:kte,jts:jte,p_ac0) =                        &
!       chem(its:ite,kts:kte,jts:jte,p_ac0) +                   &
!       factor(its:ite,kts:kte,jts:jte)*factnuma*(              &
!       anthfac*(.75*ebu(its:ite,kts:kte,jts:jte,p_ebu_pm25) +        &
!      .75*ebu(its:ite,kts:kte,jts:jte,p_ebu_bc) ) +                  &
!       orgfac*.75*ebu(its:ite,kts:kte,jts:jte,p_ebu_oc) )
!! coarse
!  chem(its:ite,kts:kte,jts:jte,p_corn) =                     &
!       chem(its:ite,kts:kte,jts:jte,p_corn) +                  &
!       factor(its:ite,kts:kte,jts:jte)*factnumc*anthfac*       &
!       ebu(its:ite,kts:kte,jts:jte,p_ebu_pm10)
!
!!
!! Increment the aerosol masses...
!!
!
!  chem(its:ite,kts:kte,jts:jte,p_ecj) =                        &
!       chem(its:ite,kts:kte,jts:jte,p_ecj) +                   &
!       .75*ebu(its:ite,kts:kte,jts:jte,p_ebu_bc)*factor(its:ite,kts:kte,jts:jte)
!
!  chem(its:ite,kts:kte,jts:jte,p_eci) =                        &
!       chem(its:ite,kts:kte,jts:jte,p_eci) +                   &
!       .25*ebu(its:ite,kts:kte,jts:jte,p_ebu_bc)*factor(its:ite,kts:kte,jts:jte)
!
!  chem(its:ite,kts:kte,jts:jte,p_orgpaj) =                     &
!       chem(its:ite,kts:kte,jts:jte,p_orgpaj) +                &
!       .75*ebu(its:ite,kts:kte,jts:jte,p_ebu_oc)*factor(its:ite,kts:kte,jts:jte)
!
!  chem(its:ite,kts:kte,jts:jte,p_orgpai) =                     &
!       chem(its:ite,kts:kte,jts:jte,p_orgpai) +                &
!       .25*ebu(its:ite,kts:kte,jts:jte,p_ebu_oc)*factor(its:ite,kts:kte,jts:jte)
!
!  chem(its:ite,kts:kte,jts:jte,p_antha) =                      &
!       chem(its:ite,kts:kte,jts:jte,p_antha) +                 &
!       ebu(its:ite,kts:kte,jts:jte,p_ebu_pm10)*factor(its:ite,kts:kte,jts:jte)
!
!  chem(its:ite,kts:kte,jts:jte,p_p25j) =                       &
!       chem(its:ite,kts:kte,jts:jte,p_p25j) +                  &
!       .75*ebu(its:ite,kts:kte,jts:jte,p_ebu_pm25)*factor(its:ite,kts:kte,jts:jte)
!
!  chem(its:ite,kts:kte,jts:jte,p_p25i) =                       &
!       chem(its:ite,kts:kte,jts:jte,p_p25i) +                  &
!       .25*ebu(its:ite,kts:kte,jts:jte,p_ebu_pm25)*factor(its:ite,kts:kte,jts:jte)
!
!   endif !end biomass burning


 else
     CALL wrf_debug(15,'emissions not read in, needs opt 3 or 23')
!    CALL wrf_error_fatal ( 'emiss_opt in namelist must be =3 or 23')
 end if

!
! Get the sea salt emissions...
!
  if( seasalt_emiss_active == 1 ) then
     call soa2_seasalt_emiss(                                    &
          dtstep, u10, v10, alt, dz8w, xland, xlat,xlong,chem,   &
          emit_crustal,do_isorropia,mgai,mgaj,caai,caaj,kai,kaj, &
          ids,ide, jds,jde, kds,kde,                             &
          ims,ime, jms,jme, kms,kme,                             &
          its,ite, jts,jte, kts,kte                              )
  end if
  if ( seasalt_emiss_active == 2) then
     call soa2_seasalt_emiss_feuntes(                                        &
               dtstep, u10, v10, alt, dz8w, xland, chem, &    
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                        )   
  endif 
  if( dust_opt == 99 ) then
      call soa2_dust_gocartemis(                                  &
           ktau,dtstep,num_soil_layers,alt,u_phy,                 &
           v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,        &
           ivgtyp,isltyp,xland,dx,g,                              &
           ids,ide, jds,jde, kds,kde,                            &
           ims,ime, jms,jme, kms,kme,                             &
           its,ite, jts,jte, kts,kte                              )
  end if
  if (pollen_opt == 1) then
       call wrf_message("Pollen emissions are added in SORGAM:MADE by ts Nov 2020, @  UMICH")
       call sorgam_vbs_pollen_emiss(                                      &
            id, dtstep, u10, v10, alt, dz8w, xland, chem,   &
            rainc, rainnc, relhum,                                            & !ts
            ids,ide, jds,jde, kds,kde,                                    &!& ! domain dims
            ims,ime, jms,jme, kms,kme,                                    &!& ! memory dims
            its,ite, jts,jte, kts,kte, pollen_emiss_active                )
  end if

  if (pollen_opt == 2) then
       call wrf_message("Pollen emissions with surface level rupture are added in SORGAM:MADE by ts Nov 2020, @  UMICH")
       call sorgam_pollen_surface_emiss(                                         &
            id, dtstep, u10, v10, alt, dz8w, xland, chem,   &
            rainc, rainnc,rain_old,relhum,ivgtyp,swdown,swdnbc,                                            & !ts
            ids,ide, jds,jde, kds,kde,                                    &!& ! domain dims
            ims,ime, jms,jme, kms,kme,                                    &!& ! memory dims
            its,ite, jts,jte, kts,kte, pollen_emiss_active,polp_emis                                     )
  end if

  if (pollen_opt == 3) then
       call wrf_message("Pollen emissions with in-atmosphere rupture are added in SORGAM:MADE by ts Nov 2020, @  UMICH")
       call sorgam_pollen_atmosphere_emiss(                                         &
            id, dtstep, u10, v10, alt, dz8w, xland, chem,   &
            rainc, rainnc,rain_old,relhum,ivgtyp,swdown,swdnbc,                                            & !ts
            ids,ide, jds,jde, kds,kde,                                    &!& ! domain dims
            ims,ime, jms,jme, kms,kme,                                    &!& ! memory dims
            its,ite, jts,jte, kts,kte, pollen_emiss_active                                     )
  end if

    IF ( ktau == 1 .and. use_initial_WFA_as_ss ) THEN
       do j = jts,jte
       do i = its,ite
       if ( xland(i,j) > 1.5 ) then
          chem(i, kms:kme, j, p_seas) =                   &
             scalar(i, kms:kme, j, p_qnwfa)  /            &
             frac_seas_wfa / fact_wfa * dens_ss
          chem(i, kms:kme, j, p_clai) =                   &
             scalar(i, kms:kme, j, p_qnwfa) /             &
             frac_seas_wfa / fact_wfa * dens_ss    
          chem(i, kms:kme, j, p_claj) =                   &
             scalar(i, kms:kme, j, p_qnwfa) /             &
             frac_seas_wfa / fact_wfa * dens_ss
          chem(i, kms:kme, j, p_naai) =                   &
             scalar(i, kms:kme, j, p_qnwfa) /             &
             frac_seas_wfa / fact_wfa * dens_ss
          chem(i, kms:kme, j, p_naaj) =                   &
             scalar(i, kms:kme, j, p_qnwfa) /             &
             frac_seas_wfa / fact_wfa * dens_ss
        endif
        enddo
        enddo
    ENDIF
    if ( chem_is_wfa ) then
  ! Update WFA/IFA - convert chem arrays to number conc.
    scalar(its:ite, kms:kme, jms:jme, p_qnwfa) = &
        ( ((chem(its:ite, kms:kme, jms:jme, p_nh4ai)               + &
           chem(its:ite, kms:kme, jms:jme, p_so4ai))/dens_nh4so4 * kappa_nh4so4)  + &
           chem(its:ite, kms:kme, jms:jme, p_pols)/dens_pols * kappa_oc + &
         ((chem(its:ite, kms:kme, jms:jme, p_poa0i)               + &
           chem(its:ite, kms:kme, jms:jme, p_poa1i)               + &
           chem(its:ite, kms:kme, jms:jme, p_poa2i)               + &
           chem(its:ite, kms:kme, jms:jme, p_poa3i)               + &
           chem(its:ite, kms:kme, jms:jme, p_asoa0i)              + &
           chem(its:ite, kms:kme, jms:jme, p_asoa1i)              + &
           chem(its:ite, kms:kme, jms:jme, p_asoa2i)              + &
           chem(its:ite, kms:kme, jms:jme, p_asoa3i)              + &
           chem(its:ite, kms:kme, jms:jme, p_bsoa1i)              + &
           chem(its:ite, kms:kme, jms:jme, p_bsoa2i)              + &
           chem(its:ite, kms:kme, jms:jme, p_bsoa3i))/dens_oc * kappa_oc)) * fact_wfa   + &
         ((chem(its:ite, kms:kme, jms:jme, p_clai)                + &
           chem(its:ite, kms:kme, jms:jme, p_naai))/dens_ss * kappa_ss) * fact_wfa_ss
     endif


! First time step, use the IFA distribution for the dust distribution
    IF ( ktau == 1 .and. use_initial_IFA_as_dust ) THEN

    chem(its:ite, kms:kme, jms:jme, p_soila) =                   &
           scalar(its:ite, kms:kme, jms:jme, p_qnifa)  /         &
           frac_dust_ifa / fact_ifa * dens_dust  
    chem(its:ite, kms:kme, jms:jme, p_p25i) =                   &
           scalar(its:ite, kms:kme, jms:jme, p_qnifa)  /         &
           fact_ifa * dens_dust                 
    chem(its:ite, kms:kme, jms:jme, p_p25j) =                   &
           scalar(its:ite, kms:kme, jms:jme, p_qnifa)  /         &
           fact_ifa * dens_dust

! Afterwards, update scalar based on dust
    ENDIF
    IF ( chem_is_ifa ) THEN

    scalar(its:ite, kms:kme, jms:jme, p_qnifa) = &
           chem(its:ite, kms:kme, jms:jme, p_p25i)/dens_dust  *  fact_ifa * kappa_dust

    ENDIF
          

END SUBROUTINE soa2_addemiss



!*---------------------------------------------------------------------
!*---------------------------------------------------------------------
!* ts----                begin pollen 
!*--------------------------------------------------------------------
!* 
   subroutine sorgam_vbs_pollen_emiss(                                     &
               id, dtstep, u10, v10, alt, dz8w, xland, chem, &
               rainc, rainnc, relhum,                                          & !ts
               ids,ide, jds,jde, kds,kde,                                  &!& ! domain dims
               ims,ime, jms,jme, kms,kme,                                  &!& ! memory dims
               its,ite, jts,jte, kts,kte, pollen_emiss_active                                   )
   USE module_state_description, only:  num_chem, param_first_scalar
   USE module_data_mosaic_asect
   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: id,                                      &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte

   REAL, INTENT(IN   ) ::    dtstep
   INTEGER, INTENT(IN ) ::    pollen_emiss_active
   REAL,  DIMENSION( ims:ime , jms:jme ),                                  &!ts
          INTENT(IN   ) ::   u10, v10, xland
   REAL,  DIMENSION( ims:ime , jms:jme ),                                  &!ts
          INTENT(IN   ) ::   rainc, rainnc
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::   chem
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   alt, dz8w, relhum
! local variables
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme) :: rh
        integer i, j, k, l, l_oc,ibin                                      !ts
        integer iphase, itype
        integer p1st
        real factaa,factbb,factcc
        real fr
        real fh
        real fw

        real :: ppemfact_numb_j, ppemfact_numb_c
        real :: ppemfact_mass_j, ppemfact_mass_c 
       
        real rupt
        real pr_low                                                          !ts ??????????
        real pr_high                                                         !ts ??????????
        real rh_low                                                          !ts ??????????
        real rh_high                                                         !ts ??????????
        real pollencount
        p1st = PARAM_FIRST_SCALAR

        !   for now just do itype=1
        itype = 1
        iphase = ai_phase
    !    rho_pp = 1200.0     ! [=] g/cm3 density of primary pollens (i=1,4) based on Wozniak and Steiner,   2017
    !    rho_spp = 1425.0   ! [=] g/cm3 density of sub-pollen particles (i=5) 
               ! based on Wozniak and Steiner, 2017
        pr_low  = 0.0  ! [mm/day] lowest rainfall at 10 m
        pr_high = 0.5  ! [mm/day] highest rainfall at 10 m
        rh_low  = 50.0 ! [%] lowest relative humidity at 10 m
        rh_high = 80.0 ! [%] highest relative humidity at 10 m
     !   frupt   = 0.70 ! [~] the fraction of moisture-exposed pollen grains that rupture
     !   nsgp    = 1000.0! [SPPs/grain] the number of SPPs produced per grain of pollen
    !    reffpollen = 20.0 ! effective dimaters ( and not radius!)  in micrometer
        pollencount = 5000.0
        rh = 100.*relhum
        k = kts
        do j = jts, jte
           do i = its, ite
              if( xland(i,j) > 1.5 ) cycle
              rupt = 0.0
              if (75 <  rh(i,1,j) < 100)  then
            !  if ( rh(i,1,j) * 1.0 > 70.0)  then
                 rupt = 1.0
              end if

              if ( rainnc(i,j) < 0.5 ) then
                  fr = 1.0
              else if ( rainnc(i,j) > 1.0) then
                  fr = 0.0
              else
                  fr = (1.0 - rainnc(i,j)) / (1.0 - 0.5)
              end if

              fw = 1.5 * ( 1.0 -  exp( -u10(i,j)/ 5.0 ))

              if ( (rh(i,1,j) * 1.0) < rh_low ) then
                  fh = 1.0
              else if ( (rh(i,1,j) * 1.0) > rh_high ) then
                  fh = 0.1
              else
                  fh = (rh_high - (rh(i,1,j) * 1.0 )) / (rh_high - rh_low)
              end if


              ppemfact_numb_c=pollencount * fh*fw*fr
              ppemfact_mass_c=ppemfact_numb_c* 3.1416/6.*1.425*(20000.66)**3*1.e-15


              factaa = (dtstep/dz8w(i,k,j))*alt(i,k,j)
              factcc = factaa * 1.0
              chem(i,k,j,p_polp) = chem(i,k,j,p_polp) +   &
                                   factcc * ppemfact_mass_c
              chem(i,k,j,p_corn)  = chem(i,k,j,p_corn) +   &
                                   factcc * ppemfact_numb_c

              factbb = factaa * 0.0
              chem(i,k,j,p_pols) = chem(i,k,j,p_pols) +   &
                                   factbb * ppemfact_mass_j
              chem(i,k,j,p_ac0)  = chem(i,k,j,p_ac0) +   &
                                   factbb * ppemfact_numb_j

           end do ! i = its, ite
        end do! j = jts, jte

      return

   END subroutine sorgam_vbs_pollen_emiss
!   Add the emissions into the chem array...

!* ts---- begin pollen_opt=2---this subroutine includes primary pollen and humidity induced surface rupture 
!*--------------------------------------------------------------------
!* 
   subroutine sorgam_pollen_surface_emiss(                                     &
               id, dtstep, u10, v10, alt, dz8w, xland, chem, &
               rainc, rainnc,rain_old,relhum,ivgtyp,swdown,swdnbc,                                    & !ts
               ids,ide, jds,jde, kds,kde,                                  &!& ! domain dims
               ims,ime, jms,jme, kms,kme,                                  &!& ! memory dims
               its,ite, jts,jte, kts,kte, pollen_emiss_active,polp_emis                                   )
   USE module_state_description, only:  num_chem, param_first_scalar
   USE module_data_mosaic_asect
   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: id,                                      &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte

   REAL, INTENT(IN   ) ::    dtstep
   INTEGER, INTENT(IN ) ::    pollen_emiss_active
   REAL,  DIMENSION( ims:ime , jms:jme ),                                  &!ts
          INTENT(IN   ) ::   u10, v10, xland,swdown,swdnbc
   REAL,  DIMENSION( ims:ime , jms:jme ),                                  &!ts
          INTENT(IN   ) ::   rainc, rainnc,rain_old
   INTEGER,  DIMENSION( ims:ime , jms:jme ),                               &
          INTENT(IN   ) ::   ivgtyp

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::   chem
   REAL, DIMENSION( ims:ime, jms:jme), INTENT(OUT ) :: polp_emis
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   alt, dz8w, relhum
   REAL, DIMENSION( ims:ime , kms:kme , jms:jme ) ::  rh !ts
! local variables
        integer i, j, k, l, l_oc,ibin                                      !ts
        integer iphase, itype
        integer p1st
        real factaa,factbb,factcc,factmm,factnn
        real fr
        real fh
        real fw
       ! real flux_n, flux_m, flux_m1, flux_m2, flux_m3, flux_m4
       ! real flux_n1, flux_n2, flux_n3, flux_n4, emispol_rupt
       ! real :: flux_n_bin( maxd_asize, maxd_atype )
       ! real :: flux_m_bin( maxd_asize, maxd_atype )
       ! real :: pol_flux(ime, kme, jme)

        real :: ppemfact_numb_j, ppemfact_numb_c
        real :: ppemfact_mass_j, ppemfact_mass_c

        real :: rain_diff
        real :: relh
        real rupt
        real pr_low                                                          !ts 
        real pr_high                                                         !ts 
        real rh_low                                                          !ts 
        real rh_high                                                         !ts
        real pollencount
        integer ierr
  !      real :: emispot(111,111)
        real, dimension(953,834) :: emispot
       p1st = PARAM_FIRST_SCALAR

        !   for now just do itype=1
        itype = 1
        iphase = ai_phase
    !    rho_pp = 1200.0     ! [=] kg/m3 density of primary pollens (i=1,4) based on Wozniak and Steiner,   2017
    !    rho_spp = 1425.0   ! [=] kg/m3 density of sub-pollen particles (i=5) 
               ! based on Wozniak and Steiner, 2017
        relh = 100 !can be removed
        pr_low  = 0.0  ! [mm/day] lowest rainfall at 10 m
        pr_high = 0.5  ! [mm/day] highest rainfall at 10 m
        rh_low  = 50.0 ! [%] lowest relative humidity at 10 m
        rh_high = 70.0 ! [%] highest relative humidity at 10 m
     !   frupt   = 0.70 ! [~] the fraction of moisture-exposed pollen grains that rupture
     !   nsgp    = 1000.0! [SPPs/grain] the number of SPPs produced per grain of pollen
    !    reffpollen = 20.0 ! effective dimaters ( and not radius!)  in micrometer

        rh = 100.*relhum        

        open (unit=2, file = 'pollendata.txt', status = 'old',iostat=ierr) 

        if( ierr /= 0 ) then
        call wrf_error_fatal( 'o3_rei: Failed to open pollendata.txt' )
        endif
      !  READ (2, *)  ((emispot(i, j), i = its,ite), j = jts,jte)
     !   allocate(emispot(111,111))
           do j = 1,834
!             do k = 1,44
              do i =1,953
                 read(2,*) emispot(i, j)
              end do
           end do
!           do j = 1,111
!              do i =1,111
!                 write(*,*) emispot(i,j)
!              end do
!           end do
!           do i = 1,111
!               write(*,*)  i
!           end do


        do j = jts, jte
        k = kts
      !  do k = kts, kte  
        do i = its, ite
              if( xland(i,j) > 1.5 ) cycle
             ! relh(i,j) = relhum(i,j) * 100
              rupt = 0.0
              rain_diff = ((rainc(i,j) + rainnc(i,j)) - rain_old(i,j))
              if ( rain_diff < 0.5 ) then
                  fr = 1.0
              else if ( rain_diff > 1.0) then
                  fr = 0.0
              else
                  fr = (1.0 - rain_diff) / (1.0 - 0.5)
              end if


              fw = 1.5 * ( 1.0 -  exp( -abs(u10(i,j))/ 5.0 ))


              if ( (relhum(i,kts,j) * 100.0) < rh_low ) then
                  fh = 1.0
              else if ( (relhum(i,kts,j) * 100.0) > rh_high ) then
                  fh = 0.1
              else
                  fh = (rh_high - (relhum(i,kts,j) * 100.0 )) / (rh_high - rh_low)
              end if
                   pollencount = 0.0
             ! end if
               if (swdnbc(i,j) < 25.0) then
                   emispot(i,j)=0.0
               end if
                 ppemfact_numb_c=emispot(i,j) * fh*fw*fr/3600. !(60.*60./72.*11.)   !28800
                 ppemfact_mass_c=ppemfact_numb_c* 3.1416/6.*1200.0*(20.0)**3*1.e-9

                 ppemfact_numb_j=0 * fh*fw*fr
                 ppemfact_mass_j=ppemfact_numb_j* 3.1416/6.*1425.0*(.150)**3*1.e-9

              factaa = (dtstep/dz8w(i,kts,j))*alt(i,kts,j)
              factcc = factaa * 1.0
              
              chem(i,kts,j,p_polp) = chem(i,kts,j,p_polp) +   &
                                   factcc * ppemfact_mass_c !10000 !ppemfact_mass_c
              chem(i,kts,j,p_corn)  = chem(i,kts,j,p_corn) +   &
                                   factcc * ppemfact_numb_c !decreased by 10000 times  100000 !ppemfact_numb_c
              polp_emis(i,j) = factcc * ppemfact_mass_c
        end do
        end do
        do j = jts, jte
        k = kts
        !   do k = kts,kte  
        do i = its, ite
               if (( relhum(i,k,j) * 100 ) < 80 ) then
                  chem(i,k,j,p_pols) = chem(i,k,j,p_pols) +   &
                                   chem(i,k,j,p_polp)*0.0*1425.0/1200.0*((0.150)**3)/((20.0)**3)*1000 !10 * fw * fh * fr !ppemfact_mass_j
                  chem(i,k,j,p_ac0)  = chem(i,k,j,p_ac0) +   &
                                  chem(i,k,j,p_corn)*0.0*1000 !100 *    
               else if (80 < ( relhum(i,k,j) * 100 ) < 100 ) then
                  chem(i,k,j,p_pols) = chem(i,k,j,p_pols) +   &
                                   chem(i,k,j,p_polp)*0.7*1425.0/1200.0*((0.150)**3)/((20.0)**3)*1000 !10 * fw * fh * fr !ppemfact_mass_j
                  chem(i,k,j,p_ac0)  = chem(i,k,j,p_ac0) +   &
                                        ((chem(i,k,j,p_polp)*0.7*6.*10.**9.)/(3.14*20.**3.*1200.))*1000.
!                                  chem(i,k,j,p_corn)*0.7*1000 
                  chem(i,k,j,p_polp) = chem(i,k,j,p_polp) -   &
                                   chem(i,k,j,p_polp)*0.7 
                  chem(i,k,j,p_corn)  = chem(i,k,j,p_corn) -   &
                                        ((chem(i,k,j,p_polp)*0.7*6.*10.**9.)/(3.14*20.**3.*1200.))
!                                  chem(i,k,j,p_corn)*0.7
                  chem(i,k,j,p_corn) = max(1.e-9,chem(i,k,j,p_corn))
                  chem(i,k,j,p_polp) = max(1.e-12,chem(i,k,j,p_polp)) 
                end if
        end do
        end do
      !  end do
      close(2)
      return

   END subroutine sorgam_pollen_surface_emiss


!* ts----                begin pollen_opt=3---this subroutine includes primary pollen and humidity induced surface + in atmosphere rupture
!*--------------------------------------------------------------------
!* 
   subroutine sorgam_pollen_atmosphere_emiss(                                     &
               id, dtstep, u10, v10, alt, dz8w, xland, chem, &
               rainc, rainnc,rain_old,relhum,ivgtyp,swdown,swdnbc,                                    & !ts
               ids,ide, jds,jde, kds,kde,                                  &!& ! domain dims
               ims,ime, jms,jme, kms,kme,                                  &!& ! memory dims
               its,ite, jts,jte, kts,kte, pollen_emiss_active                                   )
   USE module_state_description, only:  num_chem, param_first_scalar
   USE module_data_mosaic_asect
   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: id,                                      &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte

   REAL, INTENT(IN   ) ::    dtstep
   INTEGER, INTENT(IN ) ::    pollen_emiss_active
   REAL,  DIMENSION( ims:ime , jms:jme ),                                  &!ts
          INTENT(IN   ) ::   u10, v10, xland,swdown,swdnbc
   REAL,  DIMENSION( ims:ime , jms:jme ),                                  &!ts
          INTENT(IN   ) ::   rainc, rainnc,rain_old
   INTEGER,  DIMENSION( ims:ime , jms:jme ),                               &
          INTENT(IN   ) ::   ivgtyp

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::   chem
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   alt, dz8w, relhum
REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ) :: rh !ts
! local variables
        integer i, j, k, l, l_oc,ibin                                      !ts
        integer iphase, itype
        integer p1st
        real factaa,factbb,factcc,factmm,factnn
        real fr
        real fh
        real fw
        real :: ppemfact_numb_j, ppemfact_numb_c, wind
        real :: ppemfact_mass_j, ppemfact_mass_c
        real :: rain_diff
        real rupt
        real pr_low                                                          !ts
        real pr_high                                                         !ts 
        real rh_low                                                          !ts 
        real rh_high                                                         !ts 
        real pollencount
        integer ierr
  !      real :: emispot(111,111)
        real, dimension(953,834) :: emispot
       p1st = PARAM_FIRST_SCALAR

        itype = 1
        iphase = ai_phase
    !    rho_pp = 1200.0     ! [=] kg/m3 density of primary pollens (i=1,4) based on Wozniak and Steiner,   2017
    !    rho_spp = 1425.0   ! [=] kg/m3 density of sub-pollen particles (i=5) 
               ! based on Wozniak and Steiner, 2017
        pr_low  = 0.0  ! [mm/day] lowest rainfall at 10 m
        pr_high = 0.5  ! [mm/day] highest rainfall at 10 m
        rh_low  = 50.0 ! [%] lowest relative humidity at 10 m
        rh_high = 70.0 ! [%] highest relative humidity at 10 m
     !   frupt   = 0.70 ! [~] the fraction of moisture-exposed pollen grains that rupture
     !   nsgp    = 1000.0! [SPPs/grain] the number of SPPs produced per grain of pollen
    !    reffpollen = 20.0 ! effective dimaters ( and not radius!)  in micrometer
        rh = 100.*relhum
          !    if( IVGTYP(i,j) == 10 ) then
          !         pollencount = 1000.0
          !    end if

        open (unit=2, file = 'pollendata.txt', status = 'old',iostat=ierr)
        if( ierr /= 0 ) then
        call wrf_error_fatal( 'o3_rei: Failed to open pollendata.txt' )
        endif
           do j = 1,834
!             do k = 1,44
              do i =1,953
                 read(2,*) emispot(i, j)
              end do
           end do
  !         do j = 1,111
  !            do i =1,111
  !               write(*,*) emispot(i,j)
  !            end do
  !         end do
  !         do i = 1,111
  !             write(*,*)  i
  !         end do
        do j = jts, jte
         ! do k = kte,kts,-1
        !   do k = kts,kte  
        do i = its, ite

              if( xland(i,j) > 1.5 ) cycle
             ! relh(i,j) = relhum(i,j) * 100
              rupt = 0.0
              rain_diff = ((rainc(i,j) + rainnc(i,j)) - rain_old(i,j))
              if ( rain_diff < 0.5 ) then
                  fr = 1.0
              else if ( rain_diff > 1.0) then
                  fr = 0.0
              else
                  fr = (1.0 - rain_diff) / (1.0 - 0.5)
              end if

              wind = sqrt( (u10(i,j)**2 + v10(i,j)**2 ))
              fw = 1.5 * ( 1.0 -  exp( -1.*wind/ 5.0) )

              if ( (relhum(i,kts,j) * 100.0) < rh_low ) then
                  fh = 1.0
              else if ( (relhum(i,kts,j) * 100.0) > rh_high ) then
                  fh = 0.1
              else
                  fh = (rh_high - (relhum(i,kts,j) * 100.0 )) / (rh_high - rh_low)
              end if
                   pollencount = 0.0
             ! end if
               if (swdnbc(i,j) < 25.0) then
                   emispot(i,j)=0.0
               end if

                 ppemfact_numb_c=emispot(i,j) * fh*fw*fr/3600. ! (60.*60./72.*11.)     !28800.
                 ppemfact_mass_c=ppemfact_numb_c* 3.1416/6.*1200.0*(20.0)**3.*1.e-9

                 ppemfact_numb_j=0 * fh*fw*fr
                 ppemfact_mass_j=ppemfact_numb_j* 3.1416/6.*1425.0*(.150)**3.*1.e-9


              factaa = (dtstep/dz8w(i,kts,j))*alt(i,kts,j)
              factcc = factaa * 1.0
              chem(i,kts,j,p_polp) = chem(i,kts,j,p_polp) +   &
                                   factcc * ppemfact_mass_c !10000 !ppemfact_mass_c
              chem(i,kts,j,p_corn)  = chem(i,kts,j,p_corn) +   &
                                   factcc * ppemfact_numb_c !decreased by 10000 times  100000 !ppemfact_numb_c

        
        end do
        end do
        do j = jts, jte
        !k = kts
        do k = kts,kte  
        do i = its, ite
!               if (( relhum(i,k,j) * 100. ) < 80. ) then
!                  chem(i,k,j,p_pols) = chem(i,k,j,p_pols) +   &
!                                   chem(i,k,j,p_polp)*0.0*1425.0/1200.0*((0.150)**3.)/((20.0)**3.)*1000.
!                  chem(i,k,j,p_ac0)  = chem(i,k,j,p_ac0) +   &
!                                  chem(i,k,j,p_corn)*0.0*1000.     
                   ! No Rupture
               if ( relhum(i,k,j) > 0.8 ) then
                  chem(i,k,j,p_pols) = chem(i,k,j,p_pols) +   &
                                       chem(i,k,j,p_polp)*0.7*1425.0/1200.0*((0.150)**3.)/((20.0)**3.)*1000.
                  chem(i,k,j,p_ac0)  = chem(i,k,j,p_ac0) +   &
                                       ((chem(i,k,j,p_polp)*0.7*6.*10.**9.)/(3.14*20.**3.*1200.))*1000.
                  chem(i,k,j,p_polp) = chem(i,k,j,p_polp) -   &
                                   chem(i,k,j,p_polp)*0.7 
                  chem(i,k,j,p_corn)  = chem(i,k,j,p_corn) -   &
                                       ((chem(i,k,j,p_polp)*0.7*6.*10.**9.)/(3.14*20.**3.*1200.))
                 ! Make sure we did not zero out number conc or polp
                  chem(i,k,j,p_corn) = max(1.e-12,chem(i,k,j,p_corn))
                  chem(i,k,j,p_polp) = max(1.e-16,chem(i,k,j,p_polp))
               end if
        end do
        end do
        end do
      close(2)
      return

   END subroutine sorgam_pollen_atmosphere_emiss

!*---------------------------------------------------------------------
!*---------------------------------------------------------------------
!* ts----                end pollen 
!*--------------------------------------------------------------------
!* 

!------------------------------------------------------------------------
!------------------------------------------------------------------------
SUBROUTINE soa2_seasalt_emiss(                                          &
     dtstep, u10, v10, alt, dz8w, xland, xlat, xlong, chem,             &
     emit_crustal,do_isorropia,mgai,mgaj,caai,caaj,kai,kaj,             &
     ids,ide, jds,jde, kds,kde,                                         &
     ims,ime, jms,jme, kms,kme,                                         &
     its,ite, jts,jte, kts,kte                                          )
!
! Routine to calculate seasalt emissions for SOA_VBS over the time
! dtstep...
! William.Gustafson@pnl.gov; 10-May-2007
!------------------------------------------------------------------------

   USE module_mosaic_addemiss, only:    seasalt_emitfactors_1bin
   USE module_data_soa2

   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: ids,ide, jds,jde, kds,kde,            &
                                  ims,ime, jms,jme, kms,kme,            &
                                  its,ite, jts,jte, kts,kte

   REAL, INTENT(IN   ) ::    dtstep

! 10-m wind speed components (m/s)
   REAL,  DIMENSION( ims:ime , jms:jme ),                               &
          INTENT(IN   ) ::   u10, v10, xland, xlat, xlong

! trace species mixing ratios (aerosol mass = ug/kg; number = #/kg)
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),              &
         INTENT(INOUT ) ::   chem
! non-transported species
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),              &
         INTENT(INOUT ) ::   mgai,mgaj,caai,caaj,kai,kaj

! alt  = 1.0/(dry air density) in (m3/kg)
! dz8w = layer thickness in (m)
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                     &
         INTENT(IN   ) ::   alt, dz8w

   LOGICAL, INTENT(IN ) :: do_isorropia,emit_crustal
   REAL, PARAMETER :: SSTUNE = 0.75

! local variables
   integer :: i, j, k, l, l_na, l_cl, n
    integer :: p1st

    real :: dum, dumdlo, dumdhi, dumoceanfrac, dumspd10
    real :: factaa, factbb, fraccl, fracna
!liqy   
        real :: fracca, frack,  fracmg, fracso4
!liqy-20140709 

    real :: ssemfact_numb_i, ssemfact_numb_j, ssemfact_numb_c
    real :: ssemfact_mass_i, ssemfact_mass_j, ssemfact_mass_c


!   Nope, we won't because the parameterization is only valid down to
!   0.1 microns.
!   Setup in units of cm.
!    dumdlo = 0.039e-4
!    dumdhi = 0.078e-4
    ssemfact_numb_i = 0.
    ssemfact_mass_i = 0.

!   Compute emissions factors for the accumulation mode...
!   Potentially, we could go down to 0.078 microns to match the bin
!   boundary for MOSAIC, but MOSAIC is capped at 0.1 too. The upper end
!   has been chosen to match the MOSAIC bin boundary closest to two
!   standard deviations from the default bin mean diameter for the coarse
!   mode.
    dumdlo = 0.1e-4
    dumdhi = 1.250e-4
    call seasalt_emitfactors_1bin( 1, dumdlo, dumdhi,   &
         ssemfact_numb_j, dum, ssemfact_mass_j )

!   Compute emissions factors for the coarse mode...
    dumdlo = 1.25e-4
    dumdhi = 10.0e-4
    call seasalt_emitfactors_1bin( 1, dumdlo, dumdhi,   &
         ssemfact_numb_c, dum, ssemfact_mass_c )

!   Convert mass emissions factor from (g/m2/s) to (ug/m2/s)
    ssemfact_mass_i = ssemfact_mass_i*1.0e6
    ssemfact_mass_j = ssemfact_mass_j*1.0e6
    ssemfact_mass_c = ssemfact_mass_c*1.0e6

!   Loop over i,j and apply seasalt emissions
    k = kts
    do j = jts, jte
    do i = its, ite
    !Skip this point if over land. xland=1 for land and 2 for water.
    !Also, there is no way to differentiate fresh from salt water.
    !Currently, this assumes all water is salty.
       if( xland(i,j) < 1.5 .or. &
           ( xlat(i,j) .gt. 41.0 .and. xlat(i,j) .lt. 49.5 .and. &        ! No emissions over GreAT LAKES
             xlong(i,j) .gt. -93.0 .and. xlong(i,j) .lt. -75.0 ) ) cycle
    !wig: As far as I can tell, only real.exe knows the fractional breakdown
    !     of land use. So, in wrf.exe, dumoceanfrac will always be 1.
       dumoceanfrac = 1. !fraction of grid i,j that is salt water
       dumspd10 = dumoceanfrac* &
            ( (u10(i,j)*u10(i,j) + v10(i,j)*v10(i,j))**(0.5*3.41) )
!   factaa is (s*m2/kg-air)
!   factaa*ssemfact_mass(n) is (s*m2/kg-air)*(ug/m2/s) = ug/kg-air
!   factaa*ssemfact_numb(n) is (s*m2/kg-air)*( #/m2/s) =  #/kg-air
       factaa = (dtstep/dz8w(i,k,j))*alt(i,k,j)
       factbb = factaa * dumspd10
    if ( emit_crustal ) then
       fracna = mw_na_aer / (mw_na_aer + mw_cl_aer)
       fraccl = 1.0 - fracna
    else
       fracna = 10.7838/35.171
       fraccl = 19.3529/35.171
       fracca =  0.4121/35.171
       frack  =  0.3991/35.171
       fracmg =  1.2837/35.171
    endif
!   Add the emissions into the chem array...
       chem(i,k,j,p_naai) = chem(i,k,j,p_naai) +   &
                            factbb * ssemfact_mass_i * fracna * SSTUNE
       chem(i,k,j,p_clai) = chem(i,k,j,p_clai) +   &
                            factbb * ssemfact_mass_i * fraccl * SSTUNE
!     if (do_isorropia .and. emit_crustal) then
!       mgai(i,k,j)        = mgai(i,k,j) +         &
!                            factbb * ssemfact_mass_i * fracmg * SSTUNE
!       caai(i,k,j)        = caai(i,k,j) +         &
!                            factbb * ssemfact_mass_i * fracca * SSTUNE
!       kai(i,k,j)         = kai(i,k,j)  +         &
!                            factbb * ssemfact_mass_i * frack  * SSTUNE
!     endif
        chem(i,k,j,p_nu0)  = chem(i,k,j,p_nu0) +   &
                            factbb * ssemfact_numb_i*  SSTUNE
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------   
       chem(i,k,j,p_naaj) = chem(i,k,j,p_naaj) +   &
                            factbb * ssemfact_mass_j * fracna * SSTUNE
       chem(i,k,j,p_claj) = chem(i,k,j,p_claj) +   &
                            factbb * ssemfact_mass_j * fraccl * SSTUNE
!    if (do_isorropia .and. emit_crustal ) then
!       mgaj(i,k,j)        = mgaj(i,k,j) +         &
!                            factbb * ssemfact_mass_j * fracmg * SSTUNE
!       caaj(i,k,j)        = caaj(i,k,j) +         &
!                            factbb * ssemfact_mass_j * fracca * SSTUNE
!       kaj(i,k,j)         = kaj(i,k,j)  +         &
!                            factbb * ssemfact_mass_j * frack  * SSTUNE
!     endif
       chem(i,k,j,p_ac0)  = chem(i,k,j,p_ac0) +   &
                            factbb * ssemfact_numb_j * SSTUNE
!-------------------------------------------------------------------------
       chem(i,k,j,p_seas) = chem(i,k,j,p_seas) +   &
                            factbb * ssemfact_mass_c * SSTUNE * 0.33
       chem(i,k,j,p_corn) = chem(i,k,j,p_corn) +   &
                            factbb * ssemfact_numb_c * SSTUNE * 0.33
!liqy-20140709
    end do !i
    end do !j
END SUBROUTINE soa2_seasalt_emiss


SUBROUTINE soa2_seasalt_emiss_feuntes(                                  & 
     dtstep, u10, v10, alt, dz8w, xland, chem,                          &
     ids,ide, jds,jde, kds,kde,                                         &
     ims,ime, jms,jme, kms,kme,                                         &
     its,ite, jts,jte, kts,kte                                          )


   USE module_mosaic_addemiss, only:    seasalt_emit_feuntes_1bin
   USE module_data_soa2

   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: ids,ide, jds,jde, kds,kde,            &
                                  ims,ime, jms,jme, kms,kme,            &
                                  its,ite, jts,jte, kts,kte

   REAL, INTENT(IN   ) ::    dtstep

! 10-m wind speed components (m/s)
   REAL,  DIMENSION( ims:ime , jms:jme ),                               &
          INTENT(IN   ) ::   u10, v10, xland

! trace species mixing ratios (aerosol mass = ug/kg; number = #/kg)
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),              &
         INTENT(INOUT ) ::   chem

! alt  = 1.0/(dry air density) in (m3/kg)
! dz8w = layer thickness in (m)
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                     &
         INTENT(IN   ) ::   alt, dz8w

! local variables
   integer :: i, j, k, l, l_na, l_cl, n
    integer :: p1st

    real :: dum, dumdlo, dumdhi, dumoceanfrac, dumspd10
    real :: factaa, factbb, fraccl, fracna
!liqy   
        real :: fracca, frack,  fracmg, fracso4
!liqy-20140709 

    real :: ssemfact_numb_i, ssemfact_numb_j, ssemfact_numb_c
    real :: ssemfact_mass_i, ssemfact_mass_j, ssemfact_mass_c


! For fuentes seasalt
     real, save :: alpha_f1(4) =   &
                 (/ 12.328, 38.077, 102.31, 281.65 /)
    real, save :: alpha_f2(4) =   &
                 (/ 2.2958, 8.0935, 25.251, 46.80 /)
    real, save :: alpha_f3(4) =   &
                 (/ 0.00452, 0.00705, 0.00080, 0.000761 /)
    real, save :: beta_f1(4) =   &
                 (/ 0.0311, -0.031633, 0.013154, -0.0017762 /)
    real, save :: beta_f2(4) =   &
                 (/ -13.916, 35.73, -9.7651, 1.1665 /)
    real, save :: beta_f3(4) =   &
                 (/ 4747.8, 12920.0, 7313.4, 6610.0 /)
    real :: nti(4), dp0gi(4)

    real, save :: oc02um(2) = (/ 0.0, 280.0 /)              

        ! if using Feuntes et al (2011) then calculate the seawater OC dependent
        ! factors
                do i=1,4
                        nti(i)   = beta_f1(i) * oc02um(1)**2.0 & 
                                 + beta_f2(i) * oc02um(1) + beta_f3(i)
                        dp0gi(i) = alpha_f1(i) + alpha_f2(i) & 
                                 * exp(-alpha_f3(i)*oc02um(1))
                end do


!   Nope, we won't because the parameterization is only valid down to
!   0.1 microns.
!   Setup in units of cm.
!    dumdlo = 0.039e-4
!    dumdhi = 0.078e-4
    ssemfact_numb_i = 0.
    ssemfact_mass_i = 0.

!   Compute emissions factors for the accumulation mode...
!   Potentially, we could go down to 0.078 microns to match the bin
!   boundary for MOSAIC, but MOSAIC is capped at 0.1 too. The upper end
!   has been chosen to match the MOSAIC bin boundary closest to two
!   standard deviations from the default bin mean diameter for the coarse
!   mode.
    dumdlo = 0.1e-4
    dumdhi = 1.250e-4
    call seasalt_emit_feuntes_1bin( dumdlo, dumdhi,   &
         ssemfact_numb_j, dum, ssemfact_mass_j, nti, dp0gi, oc02um(1) )

!   Compute emissions factors for the coarse mode...
    dumdlo = 1.25e-4
    dumdhi = 10.0e-4
    call seasalt_emit_feuntes_1bin( dumdlo, dumdhi,   &
         ssemfact_numb_c, dum, ssemfact_mass_c, nti, dp0gi, oc02um(1) )

!   Convert mass emissions factor from (g/m2/s) to (ug/m2/s)
    ssemfact_mass_i = ssemfact_mass_i*1.0e6
    ssemfact_mass_j = ssemfact_mass_j*1.0e6
    ssemfact_mass_c = ssemfact_mass_c*1.0e6

!   Loop over i,j and apply seasalt emissions
    k = kts
    do j = jts, jte
    do i = its, ite
    !Skip this point if over land. xland=1 for land and 2 for water.
    !Also, there is no way to differentiate fresh from salt water.
    !Currently, this assumes all water is salty.
       if( xland(i,j) < 1.5 ) cycle
    !wig: As far as I can tell, only real.exe knows the fractional breakdown
    !     of land use. So, in wrf.exe, dumoceanfrac will always be 1.
       dumoceanfrac = 1. !fraction of grid i,j that is salt water
       dumspd10 = dumoceanfrac* &
            ( (u10(i,j)*u10(i,j) + v10(i,j)*v10(i,j))**(0.5*3.41) )
!   factaa is (s*m2/kg-air)
!   factaa*ssemfact_mass(n) is (s*m2/kg-air)*(ug/m2/s) = ug/kg-air
!   factaa*ssemfact_numb(n) is (s*m2/kg-air)*( #/m2/s) =  #/kg-air
       factaa = (dtstep/dz8w(i,k,j))*alt(i,k,j)
       factbb = factaa * dumspd10
       fracna = mw_na_aer / (mw_na_aer + mw_cl_aer)
       fraccl = 1.0 - fracna
!   Add the emissions into the chem array...
       chem(i,k,j,p_naai) = chem(i,k,j,p_naai) +   &
                            factbb * ssemfact_mass_i * fracna
       chem(i,k,j,p_clai) = chem(i,k,j,p_clai) +   &
                            factbb * ssemfact_mass_i * fraccl
        chem(i,k,j,p_nu0)  = chem(i,k,j,p_nu0) +   &
                            factbb * ssemfact_numb_i
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------   
       chem(i,k,j,p_naaj) = chem(i,k,j,p_naaj) +   &
                            factbb * ssemfact_mass_j * fracna
       chem(i,k,j,p_claj) = chem(i,k,j,p_claj) +   &
                            factbb * ssemfact_mass_j * fraccl
       chem(i,k,j,p_ac0)  = chem(i,k,j,p_ac0) +   &
                            factbb * ssemfact_numb_j
!-------------------------------------------------------------------------
       chem(i,k,j,p_seas) = chem(i,k,j,p_seas) +   &
                            factbb * ssemfact_mass_c
       chem(i,k,j,p_corn) = chem(i,k,j,p_corn) +   &
                            factbb * ssemfact_numb_c
!liqy-20140709
    end do !i
    end do !j


END SUBROUTINE soa2_seasalt_emiss_feuntes
!====================================================================================
!add another dust emission scheme following GOCART mechanism  --czhao
!09/17/2009
!====================================================================================
  subroutine soa2_dust_gocartemis (ktau,dt,num_soil_layers,alt,u_phy,    &
         v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,                   &
         ivgtyp,isltyp,xland,dx,g,                                         &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
  USE module_data_gocart_dust
  USE module_configure
  USE module_state_description
  USE module_model_constants, ONLY: mwdry
  USE module_data_mosaic_asect
  IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: ktau, num_soil_layers,           &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   INTEGER,DIMENSION( ims:ime , jms:jme )                  ,               &
          INTENT(IN   ) ::                                                 &
                                                     ivgtyp,               &
                                                     isltyp
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                                   chem
  REAL, DIMENSION( ims:ime, num_soil_layers, jms:jme ) ,      &
      INTENT(INOUT) ::                               smois
   REAL,  DIMENSION( ims:ime , jms:jme, 3 )                   ,               &
          INTENT(IN   ) ::    erod
   REAL,  DIMENSION( ims:ime , jms:jme )                   ,               &
          INTENT(IN   ) ::                                                 &
                                                     u10,                  &
                                                     v10,                  &
                                                     xland
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) ::                                                 &
                                                        alt,               &
                                                     dz8w,p8w,             &
                                              u_phy,v_phy,rho_phy

  REAL, INTENT(IN   ) :: dt,dx,g
!
! local variables
!

  integer :: nmx,i,j,k,ndt,imx,jmx,lmx
  integer ilwi, start_month
  real*8, DIMENSION (3) :: erodin
  real*8, DIMENSION (5) :: bems
  real*8  w10m,gwet,airden,airmas
  real*8  cdustemis,jdustemis,cdustcon,jdustcon
  real*8  cdustdens,jdustdens,mass1part,jdustdiam,cdustdiam,dp_meanvol_tmp
  real*8  dxy
  real*8  conver,converi
  real dttt
  real soilfacj,rhosoilj,rhosoilc
  real totalemis,accfrac,corfrac,rscale1,rscale2

  accfrac=0.07              ! assign 7% to accumulation mode
  corfrac=0.93              ! assign 93% to coarse mode
  rscale1=1.00  ! to account for the dust larger than 10um in radius
  rscale2=1.02  ! to account for the dust larger than 10um in radius
  accfrac=accfrac*rscale1
  corfrac=corfrac*rscale2

  rhosoilj=2.5e3
  rhosoilc=2.6e3
  soilfacj=soilfac*rhosoilj/rhosoilc

  conver=1.e-9
  converi=1.e9
!
! number of dust bins
  nmx=5
  k=kts
  do j=jts,jte
  do i=its,ite
!
! don't do dust over water!!!
     if(xland(i,j).lt.1.5)then

     ilwi=1
     start_month = 3   ! it doesn't matter, ch_dust is not a month dependent now, a constant
     w10m=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
     airmas=-(p8w(i,kts+1,j)-p8w(i,kts,j))*dx*dx/g   ! kg 

! we don't trust the u10,v10 values, if model layers are very thin near surface
     if(dz8w(i,kts,j).lt.12.)w10m=sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))
    !erodin(1)=erod(i,j,1)/dx/dx   ! czhao erod shouldn't be scaled to the area,
    !because it's a fraction
    !erodin(2)=erod(i,j,2)/dx/dx
    !erodin(3)=erod(i,j,3)/dx/dx
     erodin(1)=erod(i,j,1)
     erodin(2)=erod(i,j,2)
     erodin(3)=erod(i,j,3)
!
!  volumetric soil moisture over porosity
     gwet=smois(i,1,j)/porosity(isltyp(i,j))
     ndt=ifix(dt)
     airden=rho_phy(i,kts,j)
     dxy=dx*dx

     call soa_vbs_source_du( nmx, dt,i,j, &
                            erodin, ilwi, dxy, w10m, gwet, airden, airmas, &
                            bems,start_month,g)

!bems: kg/timestep/cell
    !sum up the dust emission from 0.1-10 um in radius 
    ! unit change from kg/timestep/cell to ug/m2/s
    totalemis=(sum(bems(1:5))/dt)*converi/dxy
     ! to account for the particles larger than 10 um radius
     ! based on assumed size distribution
    jdustemis = totalemis*accfrac   ! accumulation mode
    cdustemis = totalemis*corfrac   ! coarse mode 

         cdustcon = sum(bems(1:5))*corfrac/airmas  ! kg/kg-dryair
         cdustcon = cdustcon * converi   ! ug/kg-dryair
         jdustcon = sum(bems(1:5))*accfrac/airmas  ! kg/kg-dryair
         jdustcon = jdustcon * converi   ! ug/kg-dryair

         chem(i,k,j,p_p25j)=chem(i,k,j,p_p25j) + jdustcon
         chem(i,k,j,p_soila)=chem(i,k,j,p_soila) + cdustcon

! czhao doing dust number emission following pm10
! use soilfacj for accumulation mode because GOCART assign a less dense dust in
! accumulation mode
       chem(i,k,j,p_ac0) =  chem(i,k,j,p_ac0) + jdustcon * factnuma*soilfacj
       chem(i,k,j,p_corn) =  chem(i,k,j,p_corn) + cdustcon * factnumc*soilfac

     endif
  enddo
  enddo

end subroutine soa2_dust_gocartemis

  SUBROUTINE soa_vbs_source_du( nmx, dt1,i,j, &
                     erod, ilwi, dxy, w10m, gwet, airden, airmas, &
                     bems,month,g0)

! ****************************************************************************
! *  Evaluate the source of each dust particles size classes  (kg/m3)        
! *  by soil emission.
! *  Input:
! *         EROD      Fraction of erodible grid cell                (-)
! *                   for 1: Sand, 2: Silt, 3: Clay
! *         DUSTDEN   Dust density                                  (kg/m3)
! *         DXY       Surface of each grid cell                     (m2)
! *         AIRVOL    Volume occupy by each grid boxes              (m3)
! *         NDT1      Time step                                     (s)
! *         W10m      Velocity at the anemometer level (10meters)   (m/s)
! *         u_tresh   Threshold velocity for particule uplifting    (m/s)
! *         CH_dust   Constant to fudge the total emission of dust  (s2/m2)
! *      
! *  Output:
! *         DSRC      Source of each dust type           (kg/timestep/cell) 
! *
! *  Working:
! *         SRC       Potential source                   (kg/m/timestep/cell)
! *
! ****************************************************************************

 USE module_data_gocart_dust

  INTEGER, INTENT(IN)    :: nmx
  REAL*8,    INTENT(IN)  :: erod(ndcls)
  INTEGER, INTENT(IN)    :: ilwi,month

  REAL*8,    INTENT(IN)    :: w10m, gwet
  REAL*8,    INTENT(IN)    :: dxy
  REAL*8,    INTENT(IN)    :: airden, airmas
  REAL*8,    INTENT(OUT)   :: bems(nmx)

  REAL*8    :: den(nmx), diam(nmx)
  REAL*8    :: tsrc, u_ts0, cw, u_ts, dsrc, srce
  REAL, intent(in)    :: g0
  REAL    :: rhoa, g,dt1
  INTEGER :: i, j, n, m, k

  ! default is 1 ug s2 m-5 == 1.e-9 kg s2 m-5
  !ch_dust(:,:)=0.8D-9   ! ch_dust is defined here instead of in the
  !chemics_ini.F if with SOA_VBS  -czhao
   ch_dust(:,:)=1.0D-9  ! default 
  !ch_dust(:,:)=0.65D-9   ! ch_dust is tuned to match MODIS and AERONET
  !measurements over Sahara 
  !ch_dust(:,:)=1.0D-9*0.36  ! ch_dust is scaled to soa_vbs total dust emission

  ! executable statemenst

  DO n = 1, nmx
     ! Threshold velocity as a function of the dust density and the diameter
     ! from Bagnold (1941)
     den(n) = den_dust(n)*1.0D-3
     diam(n) = 2.0*reff_dust(n)*1.0D2
     g = g0*1.0E2
     ! Pointer to the 3 classes considered in the source data files
     m = ipoint(n)
     tsrc = 0.0
              rhoa = airden*1.0D-3
              u_ts0 = 0.13*1.0D-2*SQRT(den(n)*g*diam(n)/rhoa)* &
                   SQRT(1.0+0.006/den(n)/g/(diam(n))**2.5)/ &
                   SQRT(1.928*(1331.0*(diam(n))**1.56+0.38)**0.092-1.0)

              ! Case of surface dry enough to erode
             IF (gwet < 0.5) THEN  !  Pete's modified value
!              IF (gwet < 0.2) THEN
                 u_ts = MAX(0.0D+0,u_ts0*(1.2D+0+2.0D-1*LOG10(MAX(1.0D-3,gwet))))
              ELSE
                 ! Case of wet surface, no erosion
                 u_ts = 100.0
              END IF
              srce = frac_s(n)*erod(m)*dxy  ! (m2)
              IF (ilwi == 1 ) THEN
                 dsrc = ch_dust(n,month)*srce*w10m**2 &
                      * (w10m - u_ts)*dt1  ! (kg)
              ELSE
                 dsrc = 0.0
              END IF
              IF (dsrc < 0.0) dsrc = 0.0

              ! Update dust mixing ratio at first model level.
              !tc(n) = tc(n) + dsrc / airmas    !kg/kg-dryair -czhao
              bems(n) = dsrc     ! kg/timestep/cell

  ENDDO

END SUBROUTINE soa_vbs_source_du

! RAR, 05/30/2013: The following subroutines are removed here, because they aren't used
! soa_vbs_dust_emiss
! soa_vbs_dust_gocartemis
! soa_vbs_source_du

END Module module_made_soa2

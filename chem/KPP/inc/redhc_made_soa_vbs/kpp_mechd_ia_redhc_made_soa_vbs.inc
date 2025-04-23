! RAR: modified to handle the new SOA mechanism based on the VBS approach and multi-generational
! VOC oxidation mechanism
! Correspondence between RACM and SAPRC-99 SOA precursors
! 1) OLT -> OLE1
! 2) OLI -> OLE2
! 3) TOL -> ARO1
! 4) XYL -> ARO2
! 5) CSL -> ARO2
! 6) HC5 -> ALK4
! 7) HC8 -> ALK5
! 8) ISO -> ISO
! 9) API -> TERP
!10) LIM -> TERP
!11) SESQ-> SESQ
!
if(p_nu0.gt.1)then


     ! OLE only, RACM_SOA is OLI and OLT
     roleho  = RCONST(48)
     roleo3  = RCONST(49)
     roleno3 = RCONST(50)

     ! TOL
     rtolho =   RCONST(67)

     ! NO XYL

     ! NO CSL

     ! PAR, RACM_SOA is HC5 and HC8
     rsoaalkho  =  RCONST(88)

     ! ISO
     risoho  =  RCONST(54)
     risoo3  =  RCONST(55)
     risono3 =  RCONST(56)

     ! C10H16, RACM is seperate for API and LIM
     rc10h16ho  = RCONST(69)
     rc10h16o3  = RCONST(68)
     rc10h16no3 = RCONST(70)

     ! SESQ
     rsesqho  =  RCONST(72)
     rsesqo3  =  RCONST(71)
     rsesqno3 =  RCONST(73)

     ! IVOC
     rivocho = RCONST(87)

     
     ! NO Isoprene radical, RACM uses ISOP + HO, NO, MO2, ACO3 to deterimine branching ratio, what whould we use?
     ! Based on Lane et al., 2008, just need RO2/NO/HO2 reactions...
     rxo2no  = RCONST(59)
     rxo2xo2 = RCONST(60)
     rxo2ho2 = RCONST(62) 
! ---------------------    
  
     ! production from anthropogenic VOCs
     PRDROG(PSOAALK)= rsoaalkho*var(ind_soaalk)*var(ind_oh_nt)

     PRDROG(POLE)= roleho*var(ind_ole)*var(ind_oh_nt) + roleo3*var(ind_ole)*var(ind_o3) + roleno3*var(ind_ole)*var(ind_no3)

     PRDROG(PARO)= rtolho*var(ind_toluene)*var(ind_oh_nt)

     PRDROG(PIVOC)= rivocho*var(ind_ivoc)*var(ind_oh_nt)


     ! Biogenic
     PRDROG(PISOP)= risoho*var(ind_isop)*var(ind_oh_nt) + risoo3*var(ind_isop)*var(ind_o3) + risono3*var(ind_isop)*var(ind_no3)

     PRDROG(PTERP)= rc10h16ho*var(ind_c10h16)*var(ind_oh_nt) + rc10h16o3*var(ind_c10h16)*var(ind_o3) + rc10h16no3*var(ind_c10h16)*var(ind_no3)

     PRDROG(PSESQ)= rsesqho*var(ind_sesq)*var(ind_oh_nt) + rsesqo3*var(ind_sesq)*var(ind_o3) + rsesqno3*var(ind_sesq)*var(ind_no3)

! RAR: to calculate the branching ratios to determine high NOx versus low NOx

     PRDROG(PBRCH)= rxo2no*var(ind_no)

!     ! VDROG carrying the branching ratios
!     if (PRDROG(PBRCH)>1.E-12) then
!        ro2loss= PRDROG(PBRCH) + rxo2xo2*var(ind_xo2_nt) + rxo2ho2*var(ind_ho2_nt)
!        VDROG3_SOA2( i,k,j,PBRCH )= MIN( 1.D0,(PRDROG(PBRCH)/ro2loss) )
!     else
!        VDROG3_SOA2( i,k,j,PBRCH )= 0.
!     end if
     if (PRDROG(PBRCH) > 1.E-12) then
       if (ind_xo2_nt > 0 .and. ind_ho2_nt > 0 .and. ind_xo2_nt <= size(var) .and. ind_ho2_nt <= size(var)) then
         ro2loss = PRDROG(PBRCH) + rxo2xo2*var(ind_xo2_nt) + rxo2ho2*var(ind_ho2_nt)
         if (ro2loss /= 0.0) then
           VDROG3_SOA2(i,k,j,PBRCH) = MIN(1.D0, (PRDROG(PBRCH)/ro2loss))
         else
           VDROG3_SOA2(i,k,j,PBRCH) = 0.
         end if
       else
         VDROG3_SOA2(i,k,j,PBRCH) = 0.
       end if
     else
       VDROG3_SOA2(i,k,j,PBRCH) = 0.
     end if


     DO n = 1, ldrog_soa2-2
        VDROG3_SOA2( i,k,j,n ) =  oconv* PRDROG( n ) * DTSTEPC
        VDROG3_SOA2( i,k,j,n ) =  MAX( 0., VDROG3_SOA2( i,k,j,n ) )
     ENDDO
! All the produced OCV at bin-0 (CVASOA0) goes to the aerosol phase
! Below is the VDROG3_SOA2 component carrying CVASOA0 production
    VDROG3_SOA2(i,k,j,POCV0)= MAX( 0.0_dp, oconv* var(ind_cvasoa0) )
    var(ind_cvasoa0)= 0.0_dp

endif

! RAR: debugging
!if (i==8 .AND. j==18) then
!   if (k==1) then
!        write(*,*)'rhch5ho',rhc5ho,'rhc8ho',rhc8ho,'rhc8ho',roltho,'roliho',roliho, &
!                  'rtolho',rtolho,'rxylho',rxylho,'rsesqno3',rsesqno3
!        write(*,*)'ind_toluene',ind_toluene,'var(ind_toluene)',var(ind_toluene)
!        write(*,*)'ind_oh_nt',ind_oh_nt,'var(ind_oh_nt)',var(ind_oh_nt)
!        write(*,*)'ind_isop',ind_isop,'risoho',risoho
!        write(*,*)'PRDROG(PBRCH)', PRDROG(PBRCH),'ro2loss=',ro2loss        
!        write(*,*)'VDROG3(8,1,18,:)', VDROG3(i,k,j,:)
!   end if
!end if
!
!if (j==18 .AND. k==1) then
!        write(*,*)'VDROG3(:,18,1,:)', VDROG3(i,k,j,:)
!end if

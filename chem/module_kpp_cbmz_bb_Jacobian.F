! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! The ODE Jacobian of Chemical Model File
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
! File                 : cbmz_bb_Jacobian.f90
! Time                 : Wed Jan 29 03:51:10 2025
! Working directory    : /scratch1/BMC/wrf-chem/minsusee/WRFV4.6/chem/KPP/mechanisms/cbmz_bb
! Equation file        : cbmz_bb.kpp
! Output root filename : cbmz_bb
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE cbmz_bb_Jacobian

  USE cbmz_bb_Parameters
  USE cbmz_bb_JacobianSP

  IMPLICIT NONE

CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Jac_SP_Vec - function for sparse multiplication: sparse Jacobian times vector
!   Arguments :
!      JVS       - sparse Jacobian of variables
!      UV        - User vector for variables
!      JUV       - Jacobian times user vector
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Jac_SP_Vec ( JVS, UV, JUV )

! JVS - sparse Jacobian of variables
  REAL(kind=dp) :: JVS(LU_NONZERO)
! UV - User vector for variables
  REAL(kind=dp) :: UV(NVAR)
! JUV - Jacobian times user vector
  REAL(kind=dp) :: JUV(NVAR)

  JUV(1) = JVS(1)*UV(1)+JVS(2)*UV(7)+JVS(3)*UV(50)
  JUV(2) = JVS(4)*UV(2)
  JUV(3) = JVS(5)*UV(3)
  JUV(4) = JVS(6)*UV(4)+JVS(7)*UV(23)+JVS(8)*UV(34)+JVS(9)*UV(51)
  JUV(5) = JVS(10)*UV(5)+JVS(11)*UV(34)+JVS(12)*UV(35)+JVS(13)*UV(39)+JVS(14)*UV(45)+JVS(15)*UV(51)+JVS(16)*UV(53)&
             &+JVS(17)*UV(55)
  JUV(6) = JVS(18)*UV(6)+JVS(19)*UV(51)
  JUV(7) = JVS(20)*UV(7)+JVS(21)*UV(50)
  JUV(8) = JVS(22)*UV(8)+JVS(23)*UV(50)
  JUV(9) = JVS(24)*UV(9)+JVS(25)*UV(50)+JVS(26)*UV(55)
  JUV(10) = JVS(27)*UV(10)+JVS(28)*UV(52)+JVS(29)*UV(53)
  JUV(11) = JVS(30)*UV(11)+JVS(31)*UV(50)
  JUV(12) = JVS(32)*UV(12)+JVS(33)*UV(52)+JVS(34)*UV(54)
  JUV(13) = JVS(35)*UV(13)+JVS(36)*UV(50)
  JUV(14) = JVS(37)*UV(14)+JVS(38)*UV(34)+JVS(39)*UV(39)+JVS(40)*UV(50)+JVS(41)*UV(51)
  JUV(15) = JVS(42)*UV(15)+JVS(43)*UV(26)+JVS(44)*UV(50)+JVS(45)*UV(52)+JVS(46)*UV(54)
  JUV(16) = JVS(47)*UV(16)+JVS(48)*UV(50)+JVS(49)*UV(52)+JVS(50)*UV(55)
  JUV(17) = JVS(51)*UV(17)+JVS(52)*UV(34)+JVS(53)*UV(39)+JVS(54)*UV(40)+JVS(55)*UV(50)+JVS(56)*UV(51)
  JUV(18) = JVS(57)*UV(11)+JVS(58)*UV(13)+JVS(59)*UV(18)+JVS(60)*UV(50)+JVS(61)*UV(56)
  JUV(19) = JVS(62)*UV(19)+JVS(63)*UV(29)+JVS(64)*UV(34)+JVS(65)*UV(39)+JVS(66)*UV(41)+JVS(67)*UV(47)+JVS(68)*UV(48)&
              &+JVS(69)*UV(49)+JVS(70)*UV(50)+JVS(71)*UV(51)+JVS(72)*UV(54)+JVS(73)*UV(56)
  JUV(20) = JVS(74)*UV(20)+JVS(75)*UV(44)+JVS(76)*UV(50)+JVS(77)*UV(55)
  JUV(21) = JVS(78)*UV(21)+JVS(79)*UV(40)+JVS(80)*UV(50)+JVS(81)*UV(55)
  JUV(22) = JVS(82)*UV(22)+JVS(83)*UV(50)+JVS(84)*UV(52)+JVS(85)*UV(55)+JVS(86)*UV(56)
  JUV(23) = JVS(87)*UV(23)+JVS(88)*UV(50)+JVS(89)*UV(51)
  JUV(24) = JVS(90)*UV(24)+JVS(91)*UV(34)+JVS(92)*UV(39)+JVS(93)*UV(44)+JVS(94)*UV(50)+JVS(95)*UV(51)
  JUV(25) = JVS(96)*UV(6)+JVS(97)*UV(25)+JVS(98)*UV(51)+JVS(99)*UV(52)+JVS(100)*UV(54)+JVS(101)*UV(56)
  JUV(26) = JVS(102)*UV(11)+JVS(103)*UV(13)+JVS(104)*UV(26)+JVS(105)*UV(50)+JVS(106)*UV(54)
  JUV(27) = JVS(107)*UV(12)+JVS(108)*UV(26)+JVS(109)*UV(27)+JVS(110)*UV(36)+JVS(111)*UV(42)+JVS(112)*UV(43)+JVS(113)&
              &*UV(45)+JVS(114)*UV(50)+JVS(115)*UV(52)+JVS(116)*UV(54)+JVS(117)*UV(55)
  JUV(28) = JVS(118)*UV(23)+JVS(119)*UV(28)+JVS(120)*UV(30)+JVS(121)*UV(33)+JVS(122)*UV(34)+JVS(123)*UV(35)+JVS(124)&
              &*UV(36)+JVS(125)*UV(39)+JVS(126)*UV(42)+JVS(127)*UV(43)+JVS(128)*UV(45)+JVS(129)*UV(50)+JVS(130)*UV(51)&
              &+JVS(131)*UV(54)+JVS(132)*UV(56)
  JUV(29) = JVS(133)*UV(13)+JVS(134)*UV(19)+JVS(135)*UV(29)+JVS(136)*UV(31)+JVS(137)*UV(32)+JVS(141)*UV(45)+JVS(145)&
              &*UV(50)+JVS(147)*UV(54)+JVS(148)*UV(55)+JVS(149)*UV(56)
  JUV(30) = JVS(150)*UV(18)+JVS(151)*UV(26)+JVS(152)*UV(30)+JVS(153)*UV(50)+JVS(154)*UV(51)+JVS(156)*UV(56)
  JUV(31) = JVS(157)*UV(31)+JVS(158)*UV(35)+JVS(159)*UV(54)+JVS(160)*UV(55)+JVS(161)*UV(56)
  JUV(32) = JVS(162)*UV(32)+JVS(163)*UV(35)+JVS(164)*UV(50)+JVS(165)*UV(55)+JVS(166)*UV(56)
  JUV(33) = JVS(167)*UV(33)+JVS(168)*UV(45)+JVS(169)*UV(50)+JVS(170)*UV(55)+JVS(171)*UV(56)
  JUV(34) = JVS(172)*UV(34)+JVS(173)*UV(50)+JVS(174)*UV(51)+JVS(175)*UV(54)
  JUV(35) = JVS(176)*UV(35)+JVS(177)*UV(50)+JVS(178)*UV(51)+JVS(179)*UV(54)
  JUV(36) = JVS(180)*UV(20)+JVS(181)*UV(23)+JVS(182)*UV(24)+JVS(183)*UV(30)+JVS(184)*UV(32)+JVS(185)*UV(33)+JVS(186)&
              &*UV(34)+JVS(187)*UV(35)+JVS(188)*UV(36)+JVS(190)*UV(41)+JVS(191)*UV(44)+JVS(192)*UV(45)+JVS(193)*UV(46)&
              &+JVS(194)*UV(50)+JVS(195)*UV(51)+JVS(196)*UV(54)+JVS(198)*UV(56)
  JUV(37) = JVS(199)*UV(11)+JVS(200)*UV(13)+JVS(201)*UV(23)+JVS(202)*UV(26)+JVS(203)*UV(30)+JVS(204)*UV(34)+JVS(205)&
              &*UV(35)+JVS(206)*UV(37)+JVS(207)*UV(39)+JVS(208)*UV(42)+JVS(209)*UV(45)+JVS(210)*UV(47)+JVS(211)*UV(48)&
              &+JVS(212)*UV(49)+JVS(213)*UV(50)+JVS(214)*UV(51)+JVS(215)*UV(54)+JVS(216)*UV(55)+JVS(217)*UV(56)
  JUV(38) = JVS(218)*UV(33)+JVS(219)*UV(38)+JVS(220)*UV(39)+JVS(221)*UV(45)+JVS(222)*UV(46)+JVS(223)*UV(47)+JVS(224)&
              &*UV(48)+JVS(225)*UV(49)+JVS(226)*UV(50)+JVS(227)*UV(51)+JVS(228)*UV(54)+JVS(230)*UV(56)
  JUV(39) = JVS(231)*UV(39)+JVS(232)*UV(50)+JVS(233)*UV(51)+JVS(234)*UV(54)
  JUV(40) = JVS(235)*UV(17)+JVS(236)*UV(21)+JVS(237)*UV(34)+JVS(238)*UV(39)+JVS(239)*UV(40)+JVS(240)*UV(47)+JVS(241)&
              &*UV(48)+JVS(242)*UV(49)+JVS(243)*UV(50)+JVS(244)*UV(51)+JVS(245)*UV(54)+JVS(246)*UV(55)+JVS(247)*UV(56)
  JUV(41) = JVS(248)*UV(34)+JVS(249)*UV(39)+JVS(250)*UV(41)+JVS(251)*UV(49)+JVS(252)*UV(50)+JVS(254)*UV(54)+JVS(255)&
              &*UV(55)+JVS(256)*UV(56)
  JUV(42) = JVS(257)*UV(13)+JVS(258)*UV(30)+JVS(259)*UV(33)+JVS(260)*UV(34)+JVS(261)*UV(39)+JVS(262)*UV(42)+JVS(263)&
              &*UV(45)+JVS(264)*UV(46)+JVS(265)*UV(47)+JVS(266)*UV(50)+JVS(267)*UV(51)+JVS(268)*UV(54)+JVS(270)*UV(56)
  JUV(43) = JVS(271)*UV(8)+JVS(272)*UV(21)+JVS(273)*UV(23)+JVS(274)*UV(30)+JVS(275)*UV(31)+JVS(276)*UV(33)+JVS(277)&
              &*UV(34)+JVS(278)*UV(35)+JVS(279)*UV(39)+JVS(280)*UV(40)+JVS(281)*UV(41)+JVS(282)*UV(43)+JVS(283)*UV(45)&
              &+JVS(284)*UV(47)+JVS(285)*UV(48)+JVS(286)*UV(49)+JVS(287)*UV(50)+JVS(288)*UV(51)+JVS(289)*UV(54)+JVS(291)&
              &*UV(56)
  JUV(44) = JVS(292)*UV(14)+JVS(293)*UV(20)+JVS(294)*UV(34)+JVS(295)*UV(38)+JVS(296)*UV(39)+JVS(297)*UV(43)+JVS(298)&
              &*UV(44)+JVS(299)*UV(45)+JVS(304)*UV(50)+JVS(305)*UV(51)+JVS(306)*UV(53)+JVS(307)*UV(54)+JVS(308)*UV(55)&
              &+JVS(309)*UV(56)
  JUV(45) = JVS(310)*UV(31)+JVS(311)*UV(32)+JVS(312)*UV(35)+JVS(313)*UV(45)+JVS(314)*UV(50)+JVS(315)*UV(51)+JVS(316)&
              &*UV(54)+JVS(318)*UV(56)
  JUV(46) = JVS(319)*UV(38)+JVS(320)*UV(39)+JVS(322)*UV(46)+JVS(326)*UV(50)+JVS(327)*UV(51)+JVS(328)*UV(54)+JVS(329)&
              &*UV(55)+JVS(330)*UV(56)
  JUV(47) = JVS(331)*UV(32)+JVS(332)*UV(33)+JVS(335)*UV(46)+JVS(336)*UV(47)+JVS(337)*UV(48)+JVS(339)*UV(50)+JVS(342)&
              &*UV(55)
  JUV(48) = JVS(344)*UV(29)+JVS(347)*UV(34)+JVS(349)*UV(39)+JVS(352)*UV(47)+JVS(353)*UV(48)+JVS(355)*UV(50)+JVS(356)&
              &*UV(51)+JVS(357)*UV(54)+JVS(358)*UV(55)+JVS(359)*UV(56)
  JUV(49) = JVS(360)*UV(15)+JVS(361)*UV(18)+JVS(363)*UV(31)+JVS(364)*UV(32)+JVS(366)*UV(41)+JVS(367)*UV(45)+JVS(368)&
              &*UV(48)+JVS(369)*UV(49)+JVS(370)*UV(50)+JVS(372)*UV(52)+JVS(373)*UV(54)+JVS(374)*UV(55)+JVS(375)*UV(56)
  JUV(50) = JVS(376)*UV(6)+JVS(377)*UV(7)+JVS(378)*UV(8)+JVS(379)*UV(9)+JVS(380)*UV(11)+JVS(381)*UV(13)+JVS(382)*UV(14)&
              &+JVS(383)*UV(16)+JVS(384)*UV(17)+JVS(385)*UV(20)+JVS(386)*UV(21)+JVS(387)*UV(22)+JVS(388)*UV(23)+JVS(389)&
              &*UV(24)+JVS(390)*UV(26)+JVS(391)*UV(27)+JVS(392)*UV(28)+JVS(393)*UV(29)+JVS(394)*UV(30)+JVS(398)*UV(34)&
              &+JVS(399)*UV(35)+JVS(400)*UV(36)+JVS(401)*UV(38)+JVS(402)*UV(39)+JVS(405)*UV(42)+JVS(406)*UV(43)+JVS(408)&
              &*UV(45)+JVS(410)*UV(47)+JVS(412)*UV(49)+JVS(413)*UV(50)+JVS(414)*UV(51)+JVS(415)*UV(52)+JVS(417)*UV(54)&
              &+JVS(418)*UV(55)+JVS(419)*UV(56)
  JUV(51) = JVS(420)*UV(23)+JVS(421)*UV(25)+JVS(422)*UV(30)+JVS(423)*UV(34)+JVS(424)*UV(35)+JVS(425)*UV(39)+JVS(426)&
              &*UV(45)+JVS(427)*UV(50)+JVS(428)*UV(51)+JVS(429)*UV(52)+JVS(430)*UV(53)+JVS(432)*UV(55)+JVS(433)*UV(56)
  JUV(52) = JVS(434)*UV(10)+JVS(435)*UV(12)+JVS(436)*UV(15)+JVS(437)*UV(16)+JVS(438)*UV(18)+JVS(439)*UV(22)+JVS(440)&
              &*UV(25)+JVS(442)*UV(27)+JVS(443)*UV(31)+JVS(444)*UV(32)+JVS(445)*UV(33)+JVS(448)*UV(37)+JVS(450)*UV(40)&
              &+JVS(451)*UV(41)+JVS(454)*UV(44)+JVS(456)*UV(46)+JVS(458)*UV(48)+JVS(459)*UV(49)+JVS(460)*UV(50)+JVS(461)&
              &*UV(51)+JVS(462)*UV(52)+JVS(463)*UV(53)+JVS(464)*UV(54)+JVS(465)*UV(55)+JVS(466)*UV(56)
  JUV(53) = JVS(467)*UV(10)+JVS(468)*UV(30)+JVS(469)*UV(34)+JVS(470)*UV(35)+JVS(471)*UV(38)+JVS(472)*UV(39)+JVS(473)&
              &*UV(42)+JVS(474)*UV(43)+JVS(475)*UV(45)+JVS(476)*UV(46)+JVS(480)*UV(50)+JVS(481)*UV(51)+JVS(482)*UV(52)&
              &+JVS(483)*UV(53)+JVS(484)*UV(54)+JVS(485)*UV(55)+JVS(486)*UV(56)
  JUV(54) = JVS(487)*UV(12)+JVS(488)*UV(25)+JVS(489)*UV(26)+JVS(490)*UV(27)+JVS(491)*UV(34)+JVS(492)*UV(35)+JVS(493)&
              &*UV(36)+JVS(494)*UV(37)+JVS(495)*UV(39)+JVS(496)*UV(40)+JVS(497)*UV(41)+JVS(498)*UV(42)+JVS(499)*UV(43)&
              &+JVS(500)*UV(44)+JVS(501)*UV(45)+JVS(502)*UV(46)+JVS(504)*UV(48)+JVS(506)*UV(50)+JVS(507)*UV(51)+JVS(508)&
              &*UV(52)+JVS(509)*UV(53)+JVS(510)*UV(54)+JVS(511)*UV(55)+JVS(512)*UV(56)
  JUV(55) = JVS(513)*UV(7)+JVS(514)*UV(8)+JVS(515)*UV(9)+JVS(516)*UV(11)+JVS(517)*UV(13)+JVS(518)*UV(16)+JVS(519)*UV(18)&
              &+JVS(520)*UV(20)+JVS(521)*UV(21)+JVS(522)*UV(23)+JVS(523)*UV(24)+JVS(524)*UV(26)+JVS(525)*UV(28)+JVS(526)&
              &*UV(30)+JVS(527)*UV(31)+JVS(528)*UV(32)+JVS(529)*UV(33)+JVS(530)*UV(34)+JVS(531)*UV(35)+JVS(532)*UV(36)&
              &+JVS(533)*UV(37)+JVS(534)*UV(39)+JVS(535)*UV(40)+JVS(536)*UV(41)+JVS(537)*UV(42)+JVS(538)*UV(43)+JVS(539)&
              &*UV(44)+JVS(540)*UV(45)+JVS(541)*UV(46)+JVS(542)*UV(47)+JVS(543)*UV(48)+JVS(544)*UV(49)+JVS(545)*UV(50)&
              &+JVS(546)*UV(51)+JVS(547)*UV(52)+JVS(548)*UV(53)+JVS(549)*UV(54)+JVS(550)*UV(55)+JVS(551)*UV(56)
  JUV(56) = JVS(552)*UV(18)+JVS(553)*UV(22)+JVS(554)*UV(25)+JVS(555)*UV(31)+JVS(556)*UV(32)+JVS(557)*UV(33)+JVS(559)&
              &*UV(37)+JVS(561)*UV(40)+JVS(562)*UV(41)+JVS(564)*UV(44)+JVS(566)*UV(46)+JVS(568)*UV(48)+JVS(570)*UV(50)&
              &+JVS(571)*UV(51)+JVS(572)*UV(52)+JVS(573)*UV(53)+JVS(574)*UV(54)+JVS(575)*UV(55)+JVS(576)*UV(56)
      
END SUBROUTINE Jac_SP_Vec

! End of Jac_SP_Vec function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! JacTR_SP_Vec - sparse multiplication: sparse Jacobian transposed times vector
!   Arguments :
!      JVS       - sparse Jacobian of variables
!      UV        - User vector for variables
!      JTUV      - Jacobian transposed times user vector
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE JacTR_SP_Vec ( JVS, UV, JTUV )

! JVS - sparse Jacobian of variables
  REAL(kind=dp) :: JVS(LU_NONZERO)
! UV - User vector for variables
  REAL(kind=dp) :: UV(NVAR)
! JTUV - Jacobian transposed times user vector
  REAL(kind=dp) :: JTUV(NVAR)

  JTUV(1) = JVS(1)*UV(1)
  JTUV(2) = JVS(4)*UV(2)
  JTUV(3) = JVS(5)*UV(3)
  JTUV(4) = JVS(6)*UV(4)
  JTUV(5) = JVS(10)*UV(5)
  JTUV(6) = JVS(18)*UV(6)+JVS(96)*UV(25)+JVS(376)*UV(50)
  JTUV(7) = JVS(2)*UV(1)+JVS(20)*UV(7)+JVS(377)*UV(50)+JVS(513)*UV(55)
  JTUV(8) = JVS(22)*UV(8)+JVS(271)*UV(43)+JVS(378)*UV(50)+JVS(514)*UV(55)
  JTUV(9) = JVS(24)*UV(9)+JVS(379)*UV(50)+JVS(515)*UV(55)
  JTUV(10) = JVS(27)*UV(10)+JVS(434)*UV(52)+JVS(467)*UV(53)
  JTUV(11) = JVS(30)*UV(11)+JVS(57)*UV(18)+JVS(102)*UV(26)+JVS(199)*UV(37)+JVS(380)*UV(50)+JVS(516)*UV(55)
  JTUV(12) = JVS(32)*UV(12)+JVS(107)*UV(27)+JVS(435)*UV(52)+JVS(487)*UV(54)
  JTUV(13) = JVS(35)*UV(13)+JVS(58)*UV(18)+JVS(103)*UV(26)+JVS(133)*UV(29)+JVS(200)*UV(37)+JVS(257)*UV(42)+JVS(381)&
               &*UV(50)+JVS(517)*UV(55)
  JTUV(14) = JVS(37)*UV(14)+JVS(292)*UV(44)+JVS(382)*UV(50)
  JTUV(15) = JVS(42)*UV(15)+JVS(360)*UV(49)+JVS(436)*UV(52)
  JTUV(16) = JVS(47)*UV(16)+JVS(383)*UV(50)+JVS(437)*UV(52)+JVS(518)*UV(55)
  JTUV(17) = JVS(51)*UV(17)+JVS(235)*UV(40)+JVS(384)*UV(50)
  JTUV(18) = JVS(59)*UV(18)+JVS(150)*UV(30)+JVS(361)*UV(49)+JVS(438)*UV(52)+JVS(519)*UV(55)+JVS(552)*UV(56)
  JTUV(19) = JVS(62)*UV(19)+JVS(134)*UV(29)
  JTUV(20) = JVS(74)*UV(20)+JVS(180)*UV(36)+JVS(293)*UV(44)+JVS(385)*UV(50)+JVS(520)*UV(55)
  JTUV(21) = JVS(78)*UV(21)+JVS(236)*UV(40)+JVS(272)*UV(43)+JVS(386)*UV(50)+JVS(521)*UV(55)
  JTUV(22) = JVS(82)*UV(22)+JVS(387)*UV(50)+JVS(439)*UV(52)+JVS(553)*UV(56)
  JTUV(23) = JVS(7)*UV(4)+JVS(87)*UV(23)+JVS(118)*UV(28)+JVS(181)*UV(36)+JVS(201)*UV(37)+JVS(273)*UV(43)+JVS(388)*UV(50)&
               &+JVS(420)*UV(51)+JVS(522)*UV(55)
  JTUV(24) = JVS(90)*UV(24)+JVS(182)*UV(36)+JVS(389)*UV(50)+JVS(523)*UV(55)
  JTUV(25) = JVS(97)*UV(25)+JVS(421)*UV(51)+JVS(440)*UV(52)+JVS(488)*UV(54)+JVS(554)*UV(56)
  JTUV(26) = JVS(43)*UV(15)+JVS(104)*UV(26)+JVS(108)*UV(27)+JVS(151)*UV(30)+JVS(202)*UV(37)+JVS(390)*UV(50)+JVS(489)&
               &*UV(54)+JVS(524)*UV(55)
  JTUV(27) = JVS(109)*UV(27)+JVS(391)*UV(50)+JVS(442)*UV(52)+JVS(490)*UV(54)
  JTUV(28) = JVS(119)*UV(28)+JVS(392)*UV(50)+JVS(525)*UV(55)
  JTUV(29) = JVS(63)*UV(19)+JVS(135)*UV(29)+JVS(344)*UV(48)+JVS(393)*UV(50)
  JTUV(30) = JVS(120)*UV(28)+JVS(152)*UV(30)+JVS(183)*UV(36)+JVS(203)*UV(37)+JVS(258)*UV(42)+JVS(274)*UV(43)+JVS(394)&
               &*UV(50)+JVS(422)*UV(51)+JVS(468)*UV(53)+JVS(526)*UV(55)
  JTUV(31) = JVS(136)*UV(29)+JVS(157)*UV(31)+JVS(275)*UV(43)+JVS(310)*UV(45)+JVS(363)*UV(49)+JVS(443)*UV(52)+JVS(527)&
               &*UV(55)+JVS(555)*UV(56)
  JTUV(32) = JVS(137)*UV(29)+JVS(162)*UV(32)+JVS(184)*UV(36)+JVS(311)*UV(45)+JVS(331)*UV(47)+JVS(364)*UV(49)+JVS(444)&
               &*UV(52)+JVS(528)*UV(55)+JVS(556)*UV(56)
  JTUV(33) = JVS(121)*UV(28)+JVS(167)*UV(33)+JVS(185)*UV(36)+JVS(218)*UV(38)+JVS(259)*UV(42)+JVS(276)*UV(43)+JVS(332)&
               &*UV(47)+JVS(445)*UV(52)+JVS(529)*UV(55)+JVS(557)*UV(56)
  JTUV(34) = JVS(8)*UV(4)+JVS(11)*UV(5)+JVS(38)*UV(14)+JVS(52)*UV(17)+JVS(64)*UV(19)+JVS(91)*UV(24)+JVS(122)*UV(28)&
               &+JVS(172)*UV(34)+JVS(186)*UV(36)+JVS(204)*UV(37)+JVS(237)*UV(40)+JVS(248)*UV(41)+JVS(260)*UV(42)+JVS(277)&
               &*UV(43)+JVS(294)*UV(44)+JVS(347)*UV(48)+JVS(398)*UV(50)+JVS(423)*UV(51)+JVS(469)*UV(53)+JVS(491)*UV(54)&
               &+JVS(530)*UV(55)
  JTUV(35) = JVS(12)*UV(5)+JVS(123)*UV(28)+JVS(158)*UV(31)+JVS(163)*UV(32)+JVS(176)*UV(35)+JVS(187)*UV(36)+JVS(205)&
               &*UV(37)+JVS(278)*UV(43)+JVS(312)*UV(45)+JVS(399)*UV(50)+JVS(424)*UV(51)+JVS(470)*UV(53)+JVS(492)*UV(54)&
               &+JVS(531)*UV(55)
  JTUV(36) = JVS(110)*UV(27)+JVS(124)*UV(28)+JVS(188)*UV(36)+JVS(400)*UV(50)+JVS(493)*UV(54)+JVS(532)*UV(55)
  JTUV(37) = JVS(206)*UV(37)+JVS(448)*UV(52)+JVS(494)*UV(54)+JVS(533)*UV(55)+JVS(559)*UV(56)
  JTUV(38) = JVS(219)*UV(38)+JVS(295)*UV(44)+JVS(319)*UV(46)+JVS(401)*UV(50)+JVS(471)*UV(53)
  JTUV(39) = JVS(13)*UV(5)+JVS(39)*UV(14)+JVS(53)*UV(17)+JVS(65)*UV(19)+JVS(92)*UV(24)+JVS(125)*UV(28)+JVS(207)*UV(37)&
               &+JVS(220)*UV(38)+JVS(231)*UV(39)+JVS(238)*UV(40)+JVS(249)*UV(41)+JVS(261)*UV(42)+JVS(279)*UV(43)+JVS(296)&
               &*UV(44)+JVS(320)*UV(46)+JVS(349)*UV(48)+JVS(402)*UV(50)+JVS(425)*UV(51)+JVS(472)*UV(53)+JVS(495)*UV(54)&
               &+JVS(534)*UV(55)
  JTUV(40) = JVS(54)*UV(17)+JVS(79)*UV(21)+JVS(239)*UV(40)+JVS(280)*UV(43)+JVS(450)*UV(52)+JVS(496)*UV(54)+JVS(535)&
               &*UV(55)+JVS(561)*UV(56)
  JTUV(41) = JVS(66)*UV(19)+JVS(190)*UV(36)+JVS(250)*UV(41)+JVS(281)*UV(43)+JVS(366)*UV(49)+JVS(451)*UV(52)+JVS(497)&
               &*UV(54)+JVS(536)*UV(55)+JVS(562)*UV(56)
  JTUV(42) = JVS(111)*UV(27)+JVS(126)*UV(28)+JVS(208)*UV(37)+JVS(262)*UV(42)+JVS(405)*UV(50)+JVS(473)*UV(53)+JVS(498)&
               &*UV(54)+JVS(537)*UV(55)
  JTUV(43) = JVS(112)*UV(27)+JVS(127)*UV(28)+JVS(282)*UV(43)+JVS(297)*UV(44)+JVS(406)*UV(50)+JVS(474)*UV(53)+JVS(499)&
               &*UV(54)+JVS(538)*UV(55)
  JTUV(44) = JVS(75)*UV(20)+JVS(93)*UV(24)+JVS(191)*UV(36)+JVS(298)*UV(44)+JVS(454)*UV(52)+JVS(500)*UV(54)+JVS(539)&
               &*UV(55)+JVS(564)*UV(56)
  JTUV(45) = JVS(14)*UV(5)+JVS(113)*UV(27)+JVS(128)*UV(28)+JVS(141)*UV(29)+JVS(168)*UV(33)+JVS(192)*UV(36)+JVS(209)&
               &*UV(37)+JVS(221)*UV(38)+JVS(263)*UV(42)+JVS(283)*UV(43)+JVS(299)*UV(44)+JVS(313)*UV(45)+JVS(367)*UV(49)&
               &+JVS(408)*UV(50)+JVS(426)*UV(51)+JVS(475)*UV(53)+JVS(501)*UV(54)+JVS(540)*UV(55)
  JTUV(46) = JVS(193)*UV(36)+JVS(222)*UV(38)+JVS(264)*UV(42)+JVS(322)*UV(46)+JVS(335)*UV(47)+JVS(456)*UV(52)+JVS(476)&
               &*UV(53)+JVS(502)*UV(54)+JVS(541)*UV(55)+JVS(566)*UV(56)
  JTUV(47) = JVS(67)*UV(19)+JVS(210)*UV(37)+JVS(223)*UV(38)+JVS(240)*UV(40)+JVS(265)*UV(42)+JVS(284)*UV(43)+JVS(336)&
               &*UV(47)+JVS(352)*UV(48)+JVS(410)*UV(50)+JVS(542)*UV(55)
  JTUV(48) = JVS(68)*UV(19)+JVS(211)*UV(37)+JVS(224)*UV(38)+JVS(241)*UV(40)+JVS(285)*UV(43)+JVS(337)*UV(47)+JVS(353)&
               &*UV(48)+JVS(368)*UV(49)+JVS(458)*UV(52)+JVS(504)*UV(54)+JVS(543)*UV(55)+JVS(568)*UV(56)
  JTUV(49) = JVS(69)*UV(19)+JVS(212)*UV(37)+JVS(225)*UV(38)+JVS(242)*UV(40)+JVS(251)*UV(41)+JVS(286)*UV(43)+JVS(369)&
               &*UV(49)+JVS(412)*UV(50)+JVS(459)*UV(52)+JVS(544)*UV(55)
  JTUV(50) = JVS(3)*UV(1)+JVS(21)*UV(7)+JVS(23)*UV(8)+JVS(25)*UV(9)+JVS(31)*UV(11)+JVS(36)*UV(13)+JVS(40)*UV(14)+JVS(44)&
               &*UV(15)+JVS(48)*UV(16)+JVS(55)*UV(17)+JVS(60)*UV(18)+JVS(70)*UV(19)+JVS(76)*UV(20)+JVS(80)*UV(21)+JVS(83)&
               &*UV(22)+JVS(88)*UV(23)+JVS(94)*UV(24)+JVS(105)*UV(26)+JVS(114)*UV(27)+JVS(129)*UV(28)+JVS(145)*UV(29)&
               &+JVS(153)*UV(30)+JVS(164)*UV(32)+JVS(169)*UV(33)+JVS(173)*UV(34)+JVS(177)*UV(35)+JVS(194)*UV(36)+JVS(213)&
               &*UV(37)+JVS(226)*UV(38)+JVS(232)*UV(39)+JVS(243)*UV(40)+JVS(252)*UV(41)+JVS(266)*UV(42)+JVS(287)*UV(43)&
               &+JVS(304)*UV(44)+JVS(314)*UV(45)+JVS(326)*UV(46)+JVS(339)*UV(47)+JVS(355)*UV(48)+JVS(370)*UV(49)+JVS(413)&
               &*UV(50)+JVS(427)*UV(51)+JVS(460)*UV(52)+JVS(480)*UV(53)+JVS(506)*UV(54)+JVS(545)*UV(55)+JVS(570)*UV(56)
  JTUV(51) = JVS(9)*UV(4)+JVS(15)*UV(5)+JVS(19)*UV(6)+JVS(41)*UV(14)+JVS(56)*UV(17)+JVS(71)*UV(19)+JVS(89)*UV(23)&
               &+JVS(95)*UV(24)+JVS(98)*UV(25)+JVS(130)*UV(28)+JVS(154)*UV(30)+JVS(174)*UV(34)+JVS(178)*UV(35)+JVS(195)&
               &*UV(36)+JVS(214)*UV(37)+JVS(227)*UV(38)+JVS(233)*UV(39)+JVS(244)*UV(40)+JVS(267)*UV(42)+JVS(288)*UV(43)&
               &+JVS(305)*UV(44)+JVS(315)*UV(45)+JVS(327)*UV(46)+JVS(356)*UV(48)+JVS(414)*UV(50)+JVS(428)*UV(51)+JVS(461)&
               &*UV(52)+JVS(481)*UV(53)+JVS(507)*UV(54)+JVS(546)*UV(55)+JVS(571)*UV(56)
  JTUV(52) = JVS(28)*UV(10)+JVS(33)*UV(12)+JVS(45)*UV(15)+JVS(49)*UV(16)+JVS(84)*UV(22)+JVS(99)*UV(25)+JVS(115)*UV(27)&
               &+JVS(372)*UV(49)+JVS(415)*UV(50)+JVS(429)*UV(51)+JVS(462)*UV(52)+JVS(482)*UV(53)+JVS(508)*UV(54)+JVS(547)&
               &*UV(55)+JVS(572)*UV(56)
  JTUV(53) = JVS(16)*UV(5)+JVS(29)*UV(10)+JVS(306)*UV(44)+JVS(430)*UV(51)+JVS(463)*UV(52)+JVS(483)*UV(53)+JVS(509)&
               &*UV(54)+JVS(548)*UV(55)+JVS(573)*UV(56)
  JTUV(54) = JVS(34)*UV(12)+JVS(46)*UV(15)+JVS(72)*UV(19)+JVS(100)*UV(25)+JVS(106)*UV(26)+JVS(116)*UV(27)+JVS(131)&
               &*UV(28)+JVS(147)*UV(29)+JVS(159)*UV(31)+JVS(175)*UV(34)+JVS(179)*UV(35)+JVS(196)*UV(36)+JVS(215)*UV(37)&
               &+JVS(228)*UV(38)+JVS(234)*UV(39)+JVS(245)*UV(40)+JVS(254)*UV(41)+JVS(268)*UV(42)+JVS(289)*UV(43)+JVS(307)&
               &*UV(44)+JVS(316)*UV(45)+JVS(328)*UV(46)+JVS(357)*UV(48)+JVS(373)*UV(49)+JVS(417)*UV(50)+JVS(464)*UV(52)&
               &+JVS(484)*UV(53)+JVS(510)*UV(54)+JVS(549)*UV(55)+JVS(574)*UV(56)
  JTUV(55) = JVS(17)*UV(5)+JVS(26)*UV(9)+JVS(50)*UV(16)+JVS(77)*UV(20)+JVS(81)*UV(21)+JVS(85)*UV(22)+JVS(117)*UV(27)&
               &+JVS(148)*UV(29)+JVS(160)*UV(31)+JVS(165)*UV(32)+JVS(170)*UV(33)+JVS(216)*UV(37)+JVS(246)*UV(40)+JVS(255)&
               &*UV(41)+JVS(308)*UV(44)+JVS(329)*UV(46)+JVS(342)*UV(47)+JVS(358)*UV(48)+JVS(374)*UV(49)+JVS(418)*UV(50)&
               &+JVS(432)*UV(51)+JVS(465)*UV(52)+JVS(485)*UV(53)+JVS(511)*UV(54)+JVS(550)*UV(55)+JVS(575)*UV(56)
  JTUV(56) = JVS(61)*UV(18)+JVS(73)*UV(19)+JVS(86)*UV(22)+JVS(101)*UV(25)+JVS(132)*UV(28)+JVS(149)*UV(29)+JVS(156)&
               &*UV(30)+JVS(161)*UV(31)+JVS(166)*UV(32)+JVS(171)*UV(33)+JVS(198)*UV(36)+JVS(217)*UV(37)+JVS(230)*UV(38)&
               &+JVS(247)*UV(40)+JVS(256)*UV(41)+JVS(270)*UV(42)+JVS(291)*UV(43)+JVS(309)*UV(44)+JVS(318)*UV(45)+JVS(330)&
               &*UV(46)+JVS(359)*UV(48)+JVS(375)*UV(49)+JVS(419)*UV(50)+JVS(433)*UV(51)+JVS(466)*UV(52)+JVS(486)*UV(53)&
               &+JVS(512)*UV(54)+JVS(551)*UV(55)+JVS(576)*UV(56)
      
END SUBROUTINE JacTR_SP_Vec

! End of JacTR_SP_Vec function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE cbmz_bb_Jacobian


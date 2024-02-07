subroutine mg321_secs (mxd, x, v)
integer, intent (in) :: mxd
real , intent (in) :: x(0:mg321_nr-1)
real , intent (out) :: v(0:mg321_ivs(mxd)-1)
!! Note: We stop at degree 8 for now.
!-----------------------------------------------------------------------
integer, parameter :: m=3, n=2, nk=6, npv=43, &
  m2=m*(m-1), n2=n*(n-1), mn=m*n, m2n=m*(m-1)*n, mn2=m*n*(n-1), &
  m2n2=m2n*(n-1), m3n=m*(m-1)*(m-2)*n, mn3=m*n*(n-1)*(n-2)
integer :: i0, i1, j0, j1, k0
real  :: pv(0:npv-1), d(0:nk-1,0:nk-1), d2(0:nk-1,0:nk-1), &
  d3(0:nk-1,0:nk-1), d4(0:nk-1,0:nk-1), d5(0:nk-1,0:nk-1), &
  d6(0:nk-1,0:nk-1), d7(0:nk-1,0:nk-1)
call mg321_setd ()
pv = 0
k0 = m+n
do i0 = 0, m-1
 do i1 = 0, m-1
 if (i1.ne.i0) then
  if (2.le.mxd) then
   call mg321_deg2_i1k0 ()
  endif
  if (3.le.mxd) then
   call mg321_deg3_i1k0 ()
  endif
 endif
 enddo
enddo
do j0 = m, m+n-1
 do i0 = 0, m-1
  if (2.le.mxd) then
   call mg321_deg2_i0j0k0 ()
  endif
  if (3.le.mxd) then
   call mg321_deg3_i0j0 ()
   call mg321_deg3_i0j0k0 ()
  endif
  if (4.le.mxd) then
   call mg321_deg4_i0j0 ()
   call mg321_deg4_i0j0k0 ()
  endif
  if (5.le.mxd) then
   call mg321_deg5_i0j0 ()
   call mg321_deg5_i0j0k0 ()
  endif
  do i1 = 0, m-1
  if (i1.ne.i0) then
   if (2.le.mxd) then
    call mg321_deg2_i1j0 ()
   endif
   if (3.le.mxd) then
    call mg321_deg3_i1j0 ()
    call mg321_deg3_i1j0k0 ()
   endif
   if (4.le.mxd) then
    call mg321_deg4_i1j0 ()
    call mg321_deg4_i1j0k0 ()
   endif
   if (5.le.mxd) then
    call mg321_deg5_i1j0 ()
   endif
  endif
  enddo
 enddo
 do j1 = m, m+n-1
 if (j1.ne.j0) then
  do i0 = 0, m-1
   if (3.le.mxd) then
    call mg321_deg3_i0j1k0 ()
   endif
   if (4.le.mxd) then
    call mg321_deg4_i0j1 ()
    call mg321_deg4_i0j1k0 ()
   endif
   do i1 = 0, m-1
   if (i1.ne.i0) then
    if (3.le.mxd) then
     call mg321_deg3_i1j1 ()
    endif
   endif
   enddo
  enddo
 endif
 enddo
enddo
v(0) = 1
if (2.le.mxd) then
 v(1:4) = pv(0:3)
endif
if (3.le.mxd) then
 v(5:20) = pv(4:19)
endif
if (4.le.mxd) then
 v(21) = pv(0)*pv(0)
 v(22) = pv(0)*pv(1)
 v(23) = pv(1)*pv(1)
 v(24) = pv(0)*pv(2)
 v(25) = pv(1)*pv(2)
 v(26) = pv(2)*pv(2)
 v(27) = pv(0)*pv(3)
 v(28) = pv(1)*pv(3)
 v(29) = pv(2)*pv(3)
 v(30:49) = pv(20:39)
endif
if (5.le.mxd) then
 v(50) = pv(3)*pv(4)
 v(51) = pv(1)*pv(5)
 v(52) = pv(0)*pv(6)
 v(53) = pv(1)*pv(6)
 v(54) = pv(2)*pv(6)
 v(55) = pv(0)*pv(7)
 v(56) = pv(1)*pv(7)
 v(57) = pv(2)*pv(7)
 v(58) = pv(3)*pv(7)
 v(59) = pv(0)*pv(8)
 v(60) = pv(1)*pv(8)
 v(61) = pv(2)*pv(8)
 v(62) = pv(3)*pv(8)
 v(63) = pv(0)*pv(9)
 v(64) = pv(1)*pv(9)
 v(65) = pv(2)*pv(9)
 v(66) = pv(3)*pv(9)
 v(67) = pv(0)*pv(10)
 v(68) = pv(3)*pv(10)
 v(69) = pv(3)*pv(11)
 v(70) = pv(0)*pv(12)
 v(71) = pv(0)*pv(13)
 v(72) = pv(1)*pv(13)
 v(73) = pv(2)*pv(13)
 v(74) = pv(3)*pv(13)
 v(75) = pv(0)*pv(14)
 v(76) = pv(1)*pv(14)
 v(77) = pv(2)*pv(14)
 v(78) = pv(3)*pv(14)
 v(79) = pv(0)*pv(15)
 v(80) = pv(3)*pv(15)
 v(81) = pv(0)*pv(16)
 v(82) = pv(1)*pv(16)
 v(83) = pv(3)*pv(16)
 v(84) = pv(0)*pv(17)
 v(85) = pv(1)*pv(17)
 v(86) = pv(2)*pv(17)
 v(87) = pv(0)*pv(18)
 v(88) = pv(1)*pv(18)
 v(89) = pv(2)*pv(18)
 v(90) = pv(0)*pv(19)
 v(91) = pv(1)*pv(19)
 v(92) = pv(2)*pv(19)
 v(93:95) = pv(40:42)
endif
if (6.le.mxd) then
 v(96) = pv(8)*pv(8)
 v(97) = pv(4)*pv(9)
 v(98) = pv(5)*pv(9)
 v(99) = pv(6)*pv(9)
 v(100) = pv(7)*pv(9)
 v(101) = pv(8)*pv(9)
 v(102) = pv(9)*pv(9)
 v(103) = pv(4)*pv(11)
 v(104) = pv(4)*pv(12)
 v(105) = pv(4)*pv(13)
 v(106) = pv(9)*pv(13)
 v(107) = pv(4)*pv(14)
 v(108) = pv(5)*pv(14)
 v(109) = pv(6)*pv(14)
 v(110) = pv(7)*pv(14)
 v(111) = pv(8)*pv(14)
 v(112) = pv(9)*pv(14)
 v(113) = pv(14)*pv(14)
 v(114) = pv(4)*pv(15)
 v(115) = pv(7)*pv(15)
 v(116) = pv(9)*pv(15)
 v(117) = pv(10)*pv(15)
 v(118) = pv(4)*pv(16)
 v(119) = pv(5)*pv(16)
 v(120) = pv(6)*pv(16)
 v(121) = pv(7)*pv(16)
 v(122) = pv(8)*pv(16)
 v(123) = pv(9)*pv(16)
 v(124) = pv(10)*pv(16)
 v(125) = pv(11)*pv(16)
 v(126) = pv(14)*pv(16)
 v(127) = pv(4)*pv(18)
 v(128) = pv(9)*pv(18)
 v(129) = pv(4)*pv(19)
 v(130) = pv(7)*pv(19)
 v(131) = pv(8)*pv(19)
 v(132) = pv(9)*pv(19)
 v(133) = pv(10)*pv(19)
 v(134) = pv(11)*pv(19)
 v(135) = pv(14)*pv(19)
 v(136) = pv(0)*pv(20)
 v(137) = pv(1)*pv(20)
 v(138) = pv(3)*pv(21)
 v(139) = pv(0)*pv(22)
 v(140) = pv(2)*pv(22)
 v(141) = pv(0)*pv(23)
 v(142) = pv(1)*pv(23)
 v(143) = pv(2)*pv(23)
 v(144) = pv(3)*pv(23)
 v(145) = pv(1)*pv(24)
 v(146) = pv(2)*pv(24)
 v(147) = pv(0)*pv(25)
 v(148) = pv(1)*pv(25)
 v(149) = pv(2)*pv(25)
 v(150) = pv(3)*pv(25)
 v(151) = pv(0)*pv(26)
 v(152) = pv(1)*pv(26)
 v(153) = pv(2)*pv(26)
 v(154) = pv(0)*pv(27)
 v(155) = pv(0)*pv(28)
 v(156) = pv(0)*pv(29)
 v(157) = pv(2)*pv(29)
 v(158) = pv(3)*pv(29)
 v(159) = pv(0)*pv(30)
 v(160) = pv(1)*pv(30)
 v(161) = pv(3)*pv(30)
 v(162) = pv(0)*pv(31)
 v(163) = pv(1)*pv(31)
 v(164) = pv(2)*pv(31)
 v(165) = pv(0)*pv(32)
 v(166) = pv(1)*pv(32)
 v(167) = pv(2)*pv(32)
 v(168) = pv(3)*pv(32)
 v(169) = pv(0)*pv(33)
 v(170) = pv(1)*pv(33)
 v(171) = pv(2)*pv(33)
 v(172) = pv(0)*pv(34)
 v(173) = pv(1)*pv(34)
 v(174) = pv(2)*pv(34)
 v(175) = pv(0)*pv(35)
 v(176) = pv(1)*pv(35)
 v(177) = pv(2)*pv(35)
 v(178) = pv(0)*pv(36)
 v(179) = pv(1)*pv(36)
 v(180) = pv(2)*pv(36)
 v(181) = pv(0)*pv(37)
 v(182) = pv(1)*pv(37)
 v(183) = pv(2)*pv(37)
 v(184) = pv(0)*pv(38)
 v(185) = pv(1)*pv(38)
 v(186) = pv(2)*pv(38)
 v(187) = pv(0)*pv(39)
 v(188) = pv(1)*pv(39)
 v(189) = pv(2)*pv(39)
endif
if (7.le.mxd) then
 v(190) = pv(0)*pv(0)*pv(15)
 v(191) = pv(0)*pv(0)*pv(16)
 v(192) = pv(0)*pv(1)*pv(16)
 v(193) = pv(16)*pv(20)
 v(194) = pv(5)*pv(23)
 v(195) = pv(9)*pv(23)
 v(196) = pv(13)*pv(23)
 v(197) = pv(14)*pv(23)
 v(198) = pv(16)*pv(23)
 v(199) = pv(4)*pv(24)
 v(200) = pv(7)*pv(24)
 v(201) = pv(9)*pv(24)
 v(202) = pv(4)*pv(25)
 v(203) = pv(8)*pv(25)
 v(204) = pv(9)*pv(25)
 v(205) = pv(4)*pv(26)
 v(206) = pv(7)*pv(26)
 v(207) = pv(8)*pv(26)
 v(208) = pv(9)*pv(26)
 v(209) = pv(10)*pv(26)
 v(210) = pv(14)*pv(26)
 v(211) = pv(15)*pv(26)
 v(212) = pv(4)*pv(27)
 v(213) = pv(9)*pv(27)
 v(214) = pv(10)*pv(27)
 v(215) = pv(8)*pv(28)
 v(216) = pv(4)*pv(29)
 v(217) = pv(5)*pv(29)
 v(218) = pv(6)*pv(29)
 v(219) = pv(7)*pv(29)
 v(220) = pv(8)*pv(29)
 v(221) = pv(9)*pv(29)
 v(222) = pv(10)*pv(29)
 v(223) = pv(11)*pv(29)
 v(224) = pv(14)*pv(29)
 v(225) = pv(4)*pv(30)
 v(226) = pv(6)*pv(30)
 v(227) = pv(7)*pv(30)
 v(228) = pv(9)*pv(30)
 v(229) = pv(10)*pv(30)
 v(230) = pv(11)*pv(30)
 v(231) = pv(14)*pv(30)
 v(232) = pv(4)*pv(31)
 v(233) = pv(7)*pv(31)
 v(234) = pv(8)*pv(31)
 v(235) = pv(9)*pv(31)
 v(236) = pv(10)*pv(31)
 v(237) = pv(11)*pv(31)
 v(238) = pv(14)*pv(31)
 v(239) = pv(4)*pv(32)
 v(240) = pv(5)*pv(32)
 v(241) = pv(6)*pv(32)
 v(242) = pv(7)*pv(32)
 v(243) = pv(8)*pv(32)
 v(244) = pv(9)*pv(32)
 v(245) = pv(10)*pv(32)
 v(246) = pv(11)*pv(32)
 v(247) = pv(12)*pv(32)
 v(248) = pv(14)*pv(32)
 v(249) = pv(8)*pv(33)
 v(250) = pv(9)*pv(33)
 v(251) = pv(13)*pv(33)
 v(252) = pv(14)*pv(33)
 v(253) = pv(16)*pv(33)
 v(254) = pv(4)*pv(34)
 v(255) = pv(9)*pv(34)
 v(256) = pv(10)*pv(34)
 v(257) = pv(11)*pv(34)
 v(258) = pv(14)*pv(34)
 v(259) = pv(4)*pv(35)
 v(260) = pv(8)*pv(35)
 v(261) = pv(9)*pv(35)
 v(262) = pv(10)*pv(35)
 v(263) = pv(11)*pv(35)
 v(264) = pv(15)*pv(35)
 v(265) = pv(4)*pv(36)
 v(266) = pv(5)*pv(36)
 v(267) = pv(7)*pv(36)
 v(268) = pv(8)*pv(36)
 v(269) = pv(9)*pv(36)
 v(270) = pv(10)*pv(36)
 v(271) = pv(11)*pv(36)
 v(272) = pv(12)*pv(36)
 v(273) = pv(13)*pv(36)
 v(274) = pv(14)*pv(36)
 v(275) = pv(15)*pv(36)
 v(276) = pv(16)*pv(36)
 v(277) = pv(4)*pv(37)
 v(278) = pv(8)*pv(37)
 v(279) = pv(9)*pv(37)
 v(280) = pv(10)*pv(37)
 v(281) = pv(11)*pv(37)
 v(282) = pv(14)*pv(37)
 v(283) = pv(4)*pv(38)
 v(284) = pv(5)*pv(38)
 v(285) = pv(7)*pv(38)
 v(286) = pv(8)*pv(38)
 v(287) = pv(9)*pv(38)
 v(288) = pv(10)*pv(38)
 v(289) = pv(11)*pv(38)
 v(290) = pv(14)*pv(38)
 v(291) = pv(15)*pv(38)
 v(292) = pv(4)*pv(39)
 v(293) = pv(5)*pv(39)
 v(294) = pv(6)*pv(39)
 v(295) = pv(7)*pv(39)
 v(296) = pv(8)*pv(39)
 v(297) = pv(9)*pv(39)
 v(298) = pv(10)*pv(39)
 v(299) = pv(11)*pv(39)
 v(300) = pv(13)*pv(39)
 v(301) = pv(14)*pv(39)
 v(302) = pv(0)*pv(40)
 v(303) = pv(1)*pv(40)
 v(304) = pv(2)*pv(40)
 v(305) = pv(3)*pv(40)
 v(306) = pv(0)*pv(41)
 v(307) = pv(1)*pv(41)
 v(308) = pv(2)*pv(41)
 v(309) = pv(3)*pv(41)
 v(310) = pv(0)*pv(42)
 v(311) = pv(1)*pv(42)
 v(312) = pv(2)*pv(42)
 v(313) = pv(3)*pv(42)
endif
if (8.le.mxd) then
 v(314) = pv(3)*pv(4)*pv(16)
 v(315) = pv(1)*pv(2)*pv(23)
 v(316) = pv(1)*pv(1)*pv(24)
 v(317) = pv(1)*pv(1)*pv(26)
 v(318) = pv(1)*pv(2)*pv(26)
 v(319) = pv(0)*pv(0)*pv(28)
 v(320) = pv(0)*pv(0)*pv(31)
 v(321) = pv(0)*pv(1)*pv(31)
 v(322) = pv(0)*pv(0)*pv(32)
 v(323) = pv(0)*pv(1)*pv(32)
 v(324) = pv(1)*pv(1)*pv(34)
 v(325) = pv(0)*pv(1)*pv(36)
 v(326) = pv(1)*pv(1)*pv(36)
 v(327) = pv(1)*pv(2)*pv(36)
 v(328) = pv(0)*pv(1)*pv(37)
 v(329) = pv(1)*pv(1)*pv(37)
 v(330) = pv(0)*pv(0)*pv(38)
 v(331) = pv(1)*pv(1)*pv(38)
 v(332) = pv(0)*pv(2)*pv(38)
 v(333) = pv(0)*pv(0)*pv(39)
 v(334) = pv(0)*pv(1)*pv(39)
 v(335) = pv(23)*pv(26)
 v(336) = pv(23)*pv(28)
 v(337) = pv(20)*pv(29)
 v(338) = pv(23)*pv(29)
 v(339) = pv(26)*pv(29)
 v(340) = pv(27)*pv(29)
 v(341) = pv(29)*pv(29)
 v(342) = pv(20)*pv(30)
 v(343) = pv(23)*pv(30)
 v(344) = pv(21)*pv(31)
 v(345) = pv(23)*pv(31)
 v(346) = pv(20)*pv(32)
 v(347) = pv(21)*pv(32)
 v(348) = pv(23)*pv(32)
 v(349) = pv(24)*pv(32)
 v(350) = pv(26)*pv(32)
 v(351) = pv(27)*pv(32)
 v(352) = pv(29)*pv(32)
 v(353) = pv(29)*pv(33)
 v(354) = pv(30)*pv(33)
 v(355) = pv(32)*pv(33)
 v(356) = pv(20)*pv(34)
 v(357) = pv(23)*pv(34)
 v(358) = pv(29)*pv(34)
 v(359) = pv(32)*pv(34)
 v(360) = pv(20)*pv(35)
 v(361) = pv(23)*pv(35)
 v(362) = pv(20)*pv(36)
 v(363) = pv(23)*pv(36)
 v(364) = pv(24)*pv(36)
 v(365) = pv(26)*pv(36)
 v(366) = pv(29)*pv(36)
 v(367) = pv(32)*pv(36)
 v(368) = pv(20)*pv(37)
 v(369) = pv(21)*pv(37)
 v(370) = pv(23)*pv(37)
 v(371) = pv(25)*pv(37)
 v(372) = pv(26)*pv(37)
 v(373) = pv(29)*pv(37)
 v(374) = pv(30)*pv(37)
 v(375) = pv(32)*pv(37)
 v(376) = pv(20)*pv(38)
 v(377) = pv(21)*pv(38)
 v(378) = pv(23)*pv(38)
 v(379) = pv(24)*pv(38)
 v(380) = pv(26)*pv(38)
 v(381) = pv(29)*pv(38)
 v(382) = pv(32)*pv(38)
 v(383) = pv(20)*pv(39)
 v(384) = pv(21)*pv(39)
 v(385) = pv(22)*pv(39)
 v(386) = pv(23)*pv(39)
 v(387) = pv(24)*pv(39)
 v(388) = pv(25)*pv(39)
 v(389) = pv(26)*pv(39)
 v(390) = pv(27)*pv(39)
 v(391) = pv(28)*pv(39)
 v(392) = pv(29)*pv(39)
 v(393) = pv(4)*pv(40)
 v(394) = pv(7)*pv(40)
 v(395) = pv(9)*pv(40)
 v(396) = pv(10)*pv(40)
 v(397) = pv(11)*pv(40)
 v(398) = pv(13)*pv(40)
 v(399) = pv(14)*pv(40)
 v(400) = pv(16)*pv(40)
 v(401) = pv(17)*pv(40)
 v(402) = pv(18)*pv(40)
 v(403) = pv(19)*pv(40)
 v(404) = pv(4)*pv(41)
 v(405) = pv(5)*pv(41)
 v(406) = pv(6)*pv(41)
 v(407) = pv(7)*pv(41)
 v(408) = pv(9)*pv(41)
 v(409) = pv(10)*pv(41)
 v(410) = pv(11)*pv(41)
 v(411) = pv(14)*pv(41)
 v(412) = pv(15)*pv(41)
 v(413) = pv(17)*pv(41)
 v(414) = pv(18)*pv(41)
 v(415) = pv(19)*pv(41)
 v(416) = pv(4)*pv(42)
 v(417) = pv(5)*pv(42)
 v(418) = pv(6)*pv(42)
 v(419) = pv(7)*pv(42)
 v(420) = pv(8)*pv(42)
 v(421) = pv(9)*pv(42)
 v(422) = pv(10)*pv(42)
 v(423) = pv(11)*pv(42)
 v(424) = pv(12)*pv(42)
 v(425) = pv(13)*pv(42)
 v(426) = pv(14)*pv(42)
 v(427) = pv(15)*pv(42)
 v(428) = pv(16)*pv(42)
 v(429) = pv(17)*pv(42)
 v(430) = pv(18)*pv(42)
 v(431) = pv(19)*pv(42)
endif
return
contains
subroutine mg321_setd ()
 integer :: i, j, k, l0, l1
 integer, parameter, dimension(0:3) :: lb=(/ 0, m, m+n, m+n+1 /)
 k = 0
 do l1 = 0, 2
  do l0 = 0, l1
   do j = lb(l1), lb(l1+1)-1
    do i = lb(l0), min(j-1,lb(l0+1)-1)
     d(i,j) = x(k)
     d(j,i) = d(i,j)
     d2(i,j) = x(k)**2
     d2(j,i) = d2(i,j)
     d3(i,j) = x(k)**3
     d3(j,i) = d3(i,j)
     d4(i,j) = x(k)**4
     d4(j,i) = d4(i,j)
     d5(i,j) = x(k)**5
     d5(j,i) = d5(i,j)
     d6(i,j) = x(k)**6
     d6(j,i) = d6(i,j)
     d7(i,j) = x(k)**7
     d7(j,i) = d7(i,j)
     k = k+1
    enddo
   enddo
  enddo
 enddo
 if (k.ne.mg321_nr) then
  stop 'mg321_setd: bad count'
 endif
 do i = 0, nk-1
  d(i,i) = 0
  d2(i,i) = 0
  d3(i,i) = 0
  d4(i,i) = 0
  d5(i,i) = 0
  d6(i,i) = 0
  d7(i,i) = 0
 enddo
end subroutine mg321_setd
subroutine mg321_deg2_i1j0 ()
 pv(0) = pv(0)+d(i0,i1)*d(i0,j0)/m2n
end subroutine mg321_deg2_i1j0
subroutine mg321_deg2_i1k0 ()
 pv(1) = pv(1)+d(i0,i1)*d(i0,k0)/m2
end subroutine mg321_deg2_i1k0
subroutine mg321_deg2_i0j0k0 ()
 pv(2) = pv(2)+d(i0,j0)*d(i0,k0)/mn
 pv(3) = pv(3)+d(i0,j0)*d(j0,k0)/mn
end subroutine mg321_deg2_i0j0k0
subroutine mg321_deg3_i0j0 ()
 pv(6) = pv(6)+d3(i0,j0)/mn
end subroutine mg321_deg3_i0j0
subroutine mg321_deg3_i1j0 ()
 pv(4) = pv(4)+d2(i0,i1)*d(i0,j0)/m2n
 pv(5) = pv(5)+d(i0,i1)*d2(i0,j0)/m2n
 pv(7) = pv(7)+d(i0,i1)*d(i0,j0)*d(i1,j0)/m2n
 pv(8) = pv(8)+d2(i0,j0)*d(i1,j0)/m2n
end subroutine mg321_deg3_i1j0
subroutine mg321_deg3_i1j1 ()
 pv(9) = pv(9)+d(i0,i1)*d(i0,j0)*d(i0,j1)/m2n2
end subroutine mg321_deg3_i1j1
subroutine mg321_deg3_i1k0 ()
 pv(10) = pv(10)+d2(i0,i1)*d(i0,k0)/m2
 pv(15) = pv(15)+d(i0,i1)*d2(i0,k0)/m2
end subroutine mg321_deg3_i1k0
subroutine mg321_deg3_i0j0k0 ()
 pv(12) = pv(12)+d2(i0,j0)*d(i0,k0)/mn
 pv(16) = pv(16)+d(i0,j0)*d2(i0,k0)/mn
 pv(18) = pv(18)+d2(i0,j0)*d(j0,k0)/mn
 pv(19) = pv(19)+d(i0,j0)*d(i0,k0)*d(j0,k0)/mn
end subroutine mg321_deg3_i0j0k0
subroutine mg321_deg3_i1j0k0 ()
 pv(11) = pv(11)+d(i0,i1)*d(i0,j0)*d(i0,k0)/m2n
 pv(13) = pv(13)+d(i0,j0)*d(i1,j0)*d(i0,k0)/m2n
 pv(17) = pv(17)+d(i0,i1)*d(i0,j0)*d(j0,k0)/m2n
end subroutine mg321_deg3_i1j0k0
subroutine mg321_deg3_i0j1k0 ()
 pv(14) = pv(14)+d(i0,j0)*d(i0,j1)*d(i0,k0)/mn2
end subroutine mg321_deg3_i0j1k0
subroutine mg321_deg4_i0j0 ()
 pv(22) = pv(22)+d4(i0,j0)/mn
end subroutine mg321_deg4_i0j0
subroutine mg321_deg4_i1j0 ()
 pv(20) = pv(20)+d2(i0,i1)*d2(i0,j0)/m2n
 pv(21) = pv(21)+d(i0,i1)*d3(i0,j0)/m2n
 pv(23) = pv(23)+d2(i0,i1)*d(i0,j0)*d(i1,j0)/m2n
 pv(24) = pv(24)+d(i0,i1)*d2(i0,j0)*d(i1,j0)/m2n
 pv(25) = pv(25)+d3(i0,j0)*d(i1,j0)/m2n
end subroutine mg321_deg4_i1j0
subroutine mg321_deg4_i0j1 ()
 pv(26) = pv(26)+d3(i0,j0)*d(i0,j1)/mn2
end subroutine mg321_deg4_i0j1
subroutine mg321_deg4_i0j0k0 ()
 pv(28) = pv(28)+d3(i0,j0)*d(i0,k0)/mn
 pv(31) = pv(31)+d2(i0,j0)*d2(i0,k0)/mn
 pv(35) = pv(35)+d3(i0,j0)*d(j0,k0)/mn
 pv(38) = pv(38)+d2(i0,j0)*d(i0,k0)*d(j0,k0)/mn
 pv(39) = pv(39)+d(i0,j0)*d2(i0,k0)*d(j0,k0)/mn
end subroutine mg321_deg4_i0j0k0
subroutine mg321_deg4_i1j0k0 ()
 pv(27) = pv(27)+d(i0,i1)*d2(i0,j0)*d(i0,k0)/m2n
 pv(29) = pv(29)+d(i0,i1)*d(i0,j0)*d(i1,j0)*d(i0,k0)/m2n
 pv(30) = pv(30)+d2(i0,j0)*d(i1,j0)*d(i0,k0)/m2n
 pv(32) = pv(32)+d(i0,j0)*d(i1,j0)*d2(i0,k0)/m2n
 pv(33) = pv(33)+d2(i0,i1)*d(i0,j0)*d(j0,k0)/m2n
 pv(34) = pv(34)+d(i0,i1)*d2(i0,j0)*d(j0,k0)/m2n
 pv(37) = pv(37)+d(i0,i1)*d(i0,j0)*d(i0,k0)*d(j0,k0)/m2n
end subroutine mg321_deg4_i1j0k0
subroutine mg321_deg4_i0j1k0 ()
 pv(36) = pv(36)+d2(i0,j0)*d(i0,j1)*d(j0,k0)/mn2
end subroutine mg321_deg4_i0j1k0
subroutine mg321_deg5_i0j0 ()
 pv(41) = pv(41)+d5(i0,j0)/mn
end subroutine mg321_deg5_i0j0
subroutine mg321_deg5_i1j0 ()
 pv(40) = pv(40)+d(i0,i1)*d4(i0,j0)/m2n
end subroutine mg321_deg5_i1j0
subroutine mg321_deg5_i0j0k0 ()
 pv(42) = pv(42)+d4(i0,j0)*d(i0,k0)/mn
end subroutine mg321_deg5_i0j0k0
end subroutine mg321_secs

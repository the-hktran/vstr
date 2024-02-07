module inv_mg321_d
!..use and access
use inv_share_d
implicit none
private
public :: mg321_gens, mg321_prims, mg321_secs, mg321_base
!..types
!..data
integer, parameter, public :: &
  mg321_nr=15, mg321_ngrp=12, mg321_ngen=3, &
  mg321_dvp(0:9) = (/ 0, 5, 6, 3, 0, 0, 1, 0, 0, 0 /), &
  mg321_ivp(0:9) = (/ 0, 5, 11, 14, 14, 14, 15, 15, 15, 15 /), &
  mg321_dvs(0:9) = (/ 1, 0, 4, 16, 29, 46, 94, 124, 118, 118 /), &
  mg321_ivs(0:9) = (/ 1, 1, 5, 21, 50, 96, 190, 314, 432, 550 /), &
  mg321_dvb(0:9) = (/ 1, 5, 25, 104, 389, 1303, 4008, 11363, 30061, &
    74702 /), &
  mg321_ivb(0:9) = (/ 1, 6, 31, 135, 524, 1827, 5835, 17198, 47259, &
    121961 /)
!..procedures
contains
include "mg321_gens.f90"
include "mg321_prims.f90"
include "mg321_secs.f90"
include "mg321_base.f90"
end module inv_mg321_d

module inv_mg411_d
!..use and access
use inv_share_d
implicit none
private
public :: mg411_gens, mg411_prims, mg411_secs, mg411_base
!..types
!..data
integer, parameter, public :: &
  mg411_nr=15, mg411_ngrp=24, mg411_ngen=2, &
  mg411_dvp(0:9) = (/ 0, 4, 4, 4, 3, 0, 0, 0, 0, 0 /), &
  mg411_ivp(0:9) = (/ 0, 4, 8, 12, 15, 15, 15, 15, 15, 15 /), &
  mg411_dvs(0:9) = (/ 1, 0, 3, 13, 32, 62, 129, 221, 335, 442 /), &
  mg411_ivs(0:9) = (/ 1, 1, 4, 17, 49, 111, 240, 461, 796, 1238 /), &
  mg411_dvb(0:9) = (/ 1, 4, 17, 65, 230, 736, 2197, 6093, 15864, &
    38960 /), &
  mg411_ivb(0:9) = (/ 1, 5, 22, 87, 317, 1053, 3250, 9343, 25207, &
    64167 /)
!..procedures
contains
include "mg411_gens.f90"
include "mg411_prims.f90"
include "mg411_secs.f90"
include "mg411_base.f90"
end module inv_mg411_d

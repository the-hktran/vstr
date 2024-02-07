module inv_share_d
!..use and access
implicit none
private
public :: inv_base, mgx_mk1d, mgx_mk2d
!..types
!..data
!..procedures
contains
include "inv_base.f90"
include "mgx_mk1d.f90"
include "mgx_mk2d.f90"
end module inv_share_d

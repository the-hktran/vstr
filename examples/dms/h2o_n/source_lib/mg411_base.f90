subroutine mg411_base (mxd, x, w)
integer, intent (in) :: mxd
real , intent (in) :: x(0:mg411_nr-1)
real , intent (out) :: w(0:mg411_ivb(mxd)-1)
!-----------------------------------------------------------------------
real  :: u(0:mg411_nr-1), v(0:mg411_ivs(mxd)-1)
call mg411_prims (x, u)
call mg411_secs (mxd, x, v)
call inv_base (mxd, mg411_ivp, mg411_ivs, mg411_ivb, u, v, w)
end subroutine mg411_base

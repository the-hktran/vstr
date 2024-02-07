subroutine mg321_base (mxd, x, w)
integer, intent (in) :: mxd
real , intent (in) :: x(0:mg321_nr-1)
real , intent (out) :: w(0:mg321_ivb(mxd)-1)
!-----------------------------------------------------------------------
real  :: u(0:mg321_nr-1), v(0:mg321_ivs(mxd)-1)
call mg321_prims (x, u)
call mg321_secs (mxd, x, v)
call inv_base (mxd, mg321_ivp, mg321_ivs, mg321_ivb, u, v, w)
end subroutine mg321_base

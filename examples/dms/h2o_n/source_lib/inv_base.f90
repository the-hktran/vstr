subroutine inv_base (mxd, ivp, ivs, ivb, u, v, w)
integer, intent (in) :: mxd, ivp(0:), ivs(0:), ivb(0:)
real , intent (in) :: u(0:ivp(mxd)-1), v(0:ivs(mxd)-1)
real , intent (out) :: w(0:ivb(mxd)-1)
!-----------------------------------------------------------------------
integer :: l(0:mxd,0:ivp(mxd)), ind, i, k, d, inc
l(0,0:ivp(mxd)) = 0
w(0) = v(0)
ind = 1
do d = 1, mxd
 do k = 1, d
  do i = ivp(k-1), ivp(k)-1
   l(d,i) = ind
   inc = l(d+1-k,0)-l(d-k,i)
   w(ind:ind+inc-1) = u(i)*w(l(d-k,i):l(d-k,i)+inc-1)
   ind = ind+inc
  enddo
 enddo
 l(d,ivp(d):ivp(mxd)) = ind
 w(ind:ind+ivs(d)-ivs(d-1)-1) = v(ivs(d-1):ivs(d)-1)
 ind = ind+ivs(d)-ivs(d-1)
enddo
end subroutine inv_base

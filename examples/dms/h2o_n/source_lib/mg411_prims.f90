subroutine mg411_prims (x, u)
real , intent (in) :: x(0:mg411_nr-1)
real , intent (out) :: u(0:mg411_nr-1)
!-----------------------------------------------------------------------
integer, parameter :: m=4, m2=m*(m-1)/2
real  :: x0(0:m2-1), x1(0:m-1), x2(0:m-1), x3(0:0), t0(0:m-1)
x0 = x(0:m2-1)
x1 = x(m2:m2+m-1)
x2 = x(m2+m:m2+2*m-1)
x3(0) = x(m2+2*m)
t0(0) = (x0(0)+x0(1)+x0(3))/3
t0(1) = (x0(0)+x0(2)+x0(4))/3
t0(2) = (x0(1)+x0(2)+x0(5))/3
t0(3) = (x0(3)+x0(4)+x0(5))/3
u(0) = sum(x0)/size(x0)
u(1) = sum(x1)/size(x1)
u(2) = sum(x2)/size(x2)
u(3) = sum(x3)/size(x3)
u(4) = sum(t0**2)/size(t0)
u(5) = sum(x0**2)/size(x0)
u(6) = sum(x1**2)/size(x1)
u(7) = sum(x2**2)/size(x2)
u(8) = sum(t0**3)/size(t0)
u(9) = sum(x0**3)/size(x0)
u(10) = sum(x1**3)/size(x1)
u(11) = sum(x2**3)/size(x2)
u(12) = sum(t0**4)/size(t0)
u(13) = sum(x1**4)/size(x1)
u(14) = sum(x2**4)/size(x2)
return
end subroutine mg411_prims

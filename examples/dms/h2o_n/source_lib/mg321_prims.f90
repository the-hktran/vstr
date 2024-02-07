subroutine mg321_prims (x, u)
real , intent (in) :: x(0:mg321_nr-1)
real , intent (out) :: u(0:mg321_nr-1)
!-----------------------------------------------------------------------
integer, parameter :: m=3, n=2, m2=m*(m-1)/2, n2=n*(n-1)/2
integer :: i, j
real  :: x0(0:m2-1), x1(0:m*n-1), x2(0:n2-1), x3(0:m-1), &
  x4(0:n-1), t0(0:m-1), e(0:m-1), f(0:n-1)
x0 = x(0:m2-1)
x1 = x(m2:m2+m*n-1)
x2 = x(m2+m*n:m2+m*n+n2-1)
x3 = x(m2+m*n+n2:m2+m*n+n2+m-1)
x4 = x(m2+m*n+n2+m:m2+m*n+n2+m+n-1)
t0(0) = (x0(0)+x0(1))/2
t0(1) = (x0(0)+x0(2))/2
t0(2) = (x0(1)+x0(2))/2
do i = 0, m-1
 e(i) = sum(x1(i:i+m*(n-1):m))/n
enddo
do j = 0, n-1
 f(j) = sum(x1(m*j:m*(j+1)-1))/m
enddo
u(0) = sum(x0)/size(x0)
u(1) = sum(x1)/size(x1)
u(2) = sum(x2)/size(x2)
u(3) = sum(x3)/size(x3)
u(4) = sum(x4)/size(x4)
u(5) = sum(t0**2)/size(t0)
u(6) = sum(e**2)/size(e)
u(7) = sum(f**2)/size(f)
u(8) = sum(x1**2)/size(x1)
u(9) = sum(x3**2)/size(x3)
u(10) = sum(x4**2)/size(x4)
u(11) = sum(t0**3)/size(t0)
u(12) = sum(e**3)/size(e)
u(13) = sum(x3**3)/size(x3)
u(14) = sum(x1**6)/size(x1)
return
end subroutine mg321_prims

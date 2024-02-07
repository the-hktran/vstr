subroutine getd0 (nk, r0, d0)
implicit none
integer :: nk
real :: r0(0:nk-1,0:nk-1), d0(0:nk-1,0:nk-1)
integer :: i, j
do i = 0, nk-1
 d0(i,i) = 0
 do j = i+1, nk-1
  d0(i,j) = 3.d0*exp(-r0(i,j)/3)
  d0(j,i) = d0(i,j)
 enddo
enddo
return
end subroutine getd0

subroutine getr0 (nk, xn, r0)
implicit none
integer nk
real xn(0:2,0:nk-1), r0(0:nk-1,0:nk-1)
integer i, j
do i = 0, nk-1
 r0(i,i) = 0
 do j = i+1, nk-1
  r0(i,j) = sqrt((xn(0,j)-xn(0,i))**2+(xn(1,j)-xn(1,i))**2+ &
    (xn(2,j)-xn(2,i))**2)
  r0(j,i) = r0(i,j)
 enddo
enddo
return
end subroutine getr0

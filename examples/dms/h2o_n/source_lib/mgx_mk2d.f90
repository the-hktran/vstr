subroutine mgx_mk2d (nk, nr, la, x, d)
integer, intent (in) :: nk, nr, la(0:)
real , intent (in) :: x(0:nr-1)
real , intent (out) :: d(0:nk-1,0:nk-1)
integer :: i, j, k, l0, l1, lb(0:size(la))
lb(0) = 0
do l0 = 1, size(la)
 lb(l0) = lb(l0-1)+la(l0-1)
enddo
if (lb(size(la)).ne.nk) then
 stop 'mgx_mk2d: bad count nk'
endif
k = 0
do l1 = 0, size(lb)-1
 do l0 = 0, l1
  do j = lb(l1), lb(l1+1)-1
   do i = lb(l0), min(j-1,lb(l0+1)-1)
    d(i,j) = x(k)
    d(j,i) = d(i,j)
    k = k+1
   enddo
  enddo
 enddo
enddo
if (k.ne.nr) then
 stop 'mgx_mk2d: bad count nr'
endif
do i = 0, nk-1
 d(i,i) = 0
enddo
end subroutine mgx_mk2d

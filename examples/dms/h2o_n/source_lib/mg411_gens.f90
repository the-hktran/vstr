subroutine mg411_gens (ind, x, y)
integer, intent (in) :: ind
real , intent (in) :: x(0:mg411_nr-1)
real , intent (out) :: y(0:mg411_nr-1)
!-----------------------------------------------------------------------
integer :: iord(0:mg411_nr-1)
! Numbering (block revlex):
! ( 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14)
! (01, 02, 12, 03, 13, 23, 04, 14, 24, 34, 05, 15, 25, 35, 45)
select case (ind)
case (0, -1)
! permutation (0,1)
 iord = (/ 0, 2, 1, 4, 3, 5, 7, 6, 8, 9, 11, 10, 12, 13, 14 /)
case (1, -2)
! permutation (0,1,2,3)
 iord = (/ 2, 4, 5, 0, 1, 3, 7, 8, 9, 6, 11, 12, 13, 10, 14 /)
case default
 stop 'mg411_gens: invalid index'
end select
if (0.le.ind) then
 y(iord(:)) = x(:)
else
 y(:) = x(iord(:))
endif
return
end subroutine mg411_gens

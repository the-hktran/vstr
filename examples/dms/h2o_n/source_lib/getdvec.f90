subroutine getdvec (ma, mb, xn, vec, chg)
use inv_share_d
use inv_mg321_d
use inv_mg411_d
implicit none
! version for H4O2
integer, parameter :: nk=6, nk2=nk*(nk-1)/2
integer :: ma, mb
real  :: xn(0:2,0:nk-1), vec(0:3,0:ma+mb-1)
integer :: i, j, k, dma, dmb, ind(0:nk-1)
real  :: r0(0:nk-1,0:nk-1), d0(0:nk-1,0:nk-1), &
  d1(0:nk-1,0:nk-1), x0(0:nk2-1), v0(0:ma+mb-1), t0
!add by yimin to get charges
real::chg(0:5,0:ma+mb-1)
chg=0.d0
!add by yimin to get charges
!-----------------------------------------------------------------------
dma = -1 ; dmb = -1
do k = 9, 0, -1
 if (ma.eq.mg321_ivb(k)) then
  dma = k
 endif
 if (mb.eq.mg411_ivb(k)) then
  dmb = k
 endif
enddo
if (dma.eq.-1.or.dmb.eq.-1) then
 stop 'getvec: bad ma, mb'
endif
vec = 0
call getr0 (nk, xn, r0)
call getd0 (nk, r0, d0)
! vector factor xn(0:2,0) (species X)
ind = (/ 1, 2, 3, 4, 5, 0 /)
call sortd0 (nk, d0, ind, d1)
call mgx_mk1d (nk, nk2, (/3,2,1/), d1, x0)
call mg321_base (dma, x0, v0)
do k = 0, ma-1
 vec(0:2,k) = vec(0:2,k)+xn(0:2,0)*v0(k)
 !add by yimin to get charges
 chg(0,k)=v0(k)
 !add by yimin to get charges
 vec(3,k) = vec(3,k)+v0(k)
enddo
! vector factor xn(0:2,1) (species X)
ind = (/ 0, 2, 3, 4, 5, 1 /)
call sortd0 (nk, d0, ind, d1)
call mgx_mk1d (nk, nk2, (/3,2,1/), d1, x0)
call mg321_base (dma, x0, v0)
do k = 0, ma-1
 vec(0:2,k) = vec(0:2,k)+xn(0:2,1)*v0(k)
 !add by yimin to get charges
 chg(1,k)=v0(k)
 !add by yimin to get charges
 vec(3,k) = vec(3,k)+v0(k)
enddo
! vector factor xn(0:2,2) (species X)
ind = (/ 0, 1, 3, 4, 5, 2 /)
call sortd0 (nk, d0, ind, d1)
call mgx_mk1d (nk, nk2, (/3,2,1/), d1, x0)
call mg321_base (dma, x0, v0)
do k = 0, ma-1
 vec(0:2,k) = vec(0:2,k)+xn(0:2,2)*v0(k)
 !add by yimin to get charges
 chg(2,k)=v0(k)
 !add by yimin to get charges
 vec(3,k) = vec(3,k)+v0(k)
enddo
! vector factor xn(0:2,3) (species X)
ind = (/ 0, 1, 2, 4, 5, 3 /)
call sortd0 (nk, d0, ind, d1)
call mgx_mk1d (nk, nk2, (/3,2,1/), d1, x0)
call mg321_base (dma, x0, v0)
do k = 0, ma-1
 vec(0:2,k) = vec(0:2,k)+xn(0:2,3)*v0(k)
 !add by yimin to get charges
 chg(3,k)=v0(k)
 !add by yimin to get charges
 vec(3,k) = vec(3,k)+v0(k)
enddo
! vector factor xn(0:2,4) (species Y)
ind = (/ 0, 1, 2, 3, 5, 4 /)
call sortd0 (nk, d0, ind, d1)
call mgx_mk1d (nk, nk2, (/4,1,1/), d1, x0)
call mg411_base (dmb, x0, v0)
do k = 0, mb-1
 vec(0:2,ma+k) = vec(0:2,ma+k)+xn(0:2,4)*v0(k)
 !add by yimin to get charges
 chg(4,ma+k)=v0(k)
 !add by yimin to get charges
 vec(3,ma+k) = vec(3,ma+k)+v0(k)
enddo
! vector factor xn(0:2,5) (species Y)
ind = (/ 0, 1, 2, 3, 4, 5 /)
call sortd0 (nk, d0, ind, d1)
call mgx_mk1d (nk, nk2, (/4,1,1/), d1, x0)
call mg411_base (dmb, x0, v0)
do k = 0, mb-1
 vec(0:2,ma+k) = vec(0:2,ma+k)+xn(0:2,5)*v0(k)
 !add by yimin to get charges
 chg(5,ma+k)=v0(k)
 !add by yimin to get charges
 vec(3,ma+k) = vec(3,ma+k)+v0(k)
enddo
return
contains
subroutine sortd0 (nk, d0, ind, d1)
integer, intent (in) :: nk
real , intent (in) :: d0(0:nk-1,0:nk-1)
integer, intent (in) :: ind(0:nk-1)
real , intent (out) :: d1(0:nk-1,0:nk-1)
do j = 0, nk-1
 do i = 0, nk-1
  d1(i,j) = d0(ind(i),ind(j))
 enddo
enddo
end subroutine sortd0
end

SUBROUTINE calc_h2o_n_dip(natm,x,dip)
use dms_shell
!implicit none

integer, intent(in) :: natm
double precision, dimension(natm, 3), intent(in) :: x
integer::i,j,dim
double precision, intent(out)::dip(3)
double precision, dimension(3*natm)::xx

dim=3*natm
call dms_init(dim/9)
  do i=1,3
    do j=1,natm
      xx(3*(j-1)+i) = x(j, i) !x(natm*(i-1)+j)
    end do
  end do
!x=x/auang
call dipole(xx,dip)

end

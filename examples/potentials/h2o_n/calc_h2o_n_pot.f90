  SUBROUTINE calc_h2o_n_pot(natm,x,eng,mu)
  use pes_shell
  !implicit none
  integer, intent(in)::natm
  double precision,dimension(natm, 3),intent(in)::x
  real,dimension(3*natm)::tmpx
  double precision,dimension(3,natm)::gd1,gd2
  real,dimension(3*natm,3*natm)::H
  integer::i,j,ierr
  character::symb
  real::time_start,time_end,eng_ccsd,engf,engb
  real::eps
  double precision, intent(out)::eng
  double precision,dimension(3), intent(out)::mu
  double precision, dimension(3*natm)::xx
 
  call pes_init(natm/3)
  eps = 0.001
  do i=1,3
    do j=1,natm
      xx(3*(j-1)+i) = x(j, i) !x(natm*(i-1)+j)
    end do
  end do
  !x = x / auang
    
  call pot_gd(xx,eng,gd1,mu)
  end

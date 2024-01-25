  SUBROUTINE calc_h2o_n_pot(natm,x,eng,mu)
  use pes_shell
  !implicit none
  integer, intent(in)::natm
  double precision,dimension(3*natm),intent(in)::x
  real,dimension(3*natm)::tmpx
  double precision,dimension(3,natm)::gd1,gd2
  real,dimension(3*natm,3*natm)::H
  integer::i,j,ierr
  character::symb
  real::time_start,time_end,eng_ccsd,engf,engb
  real::eps
  double precision, intent(out)::eng
  double precision,dimension(3), intent(out)::mu
 
  call pes_init(natm/3)
  eps = 0.001

     !x = x / auang
    
  call pot_gd(x,eng,gd1,mu)
  end

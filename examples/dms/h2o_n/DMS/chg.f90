!==================================================!
! LTP2011 DMS                              !
!==================================================!
subroutine dip_ltp2011(xx,dp,chg)
  real,dimension(3,3),intent(in)::xx
  real,dimension(3),intent(out)::dp
  real,dimension(3),intent(out)::chg
  ! ::::::::::::::::::::
  !real,dimension(1:size(x,1),3)::xn
  !integer::natm
  real::r1,r2,cos_theta,mux,muy
  real::htheta,sin_htheta,cos_htheta
  real,dimension(9)::x1d

  call mon_rtheta(xx,r1,r2,cos_theta)
  call DIPS(muy,mux,r1,r2,cos_theta)
  
  htheta=0.5*acos(cos_theta)
  sin_htheta=sin(htheta)
  cos_htheta=cos(htheta)

  chg(1)=(muy*cos_htheta + mux*sin_htheta) /(2*r1*sin_htheta*cos_htheta)
  chg(2)=(mux*sin_htheta - muy*cos_htheta) /(2*r2*sin_htheta*cos_htheta)
  chg(3)= -chg(1)-chg(2)

  dp=matmul(xx,chg)

  return
end subroutine dip_ltp2011

!=================================================
! calculate r1, r2 and cos(theta) of a monomer, for the use of LTP2011
!=================================================
subroutine mon_rtheta(xx,r1,r2,cos_theta)
  real,dimension(3,3),intent(in)::xx
  real,intent(out)::r1,r2,cos_theta
  real,dimension(3)::vec1,vec2
  real::r1r2
  
  vec1=xx(:,3)-xx(:,1)
  r1=dsqrt(dot_product(vec1,vec1))
  vec2=xx(:,3)-xx(:,2)
  r2=dsqrt(dot_product(vec2,vec2))

  r1r2=dot_product(vec1,vec2)
  cos_theta=r1r2/(r1*r2)

  return
end subroutine mon_rtheta

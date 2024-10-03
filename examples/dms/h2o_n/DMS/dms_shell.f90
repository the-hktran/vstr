module dms_shell
  implicit none
  ! constants
  real,parameter::auang=0.5291772083
  real,parameter::aucm=219474.6313710
  real,parameter::aukcal=627.51

  ! some global variables
  integer,dimension(:,:),allocatable::idx_2b

contains
  !=================================================!
  ! Initializing HBB water potential                !
  !=================================================!
  subroutine dms_init(nw)
    integer,intent(in)::nw
    !::::::::::::::::::::

    ! 2-body init
    call predip('/insomnia001/home/hkt2112/Code/vstr/examples/dms/h2o_n/coef/h4o2.dms2b.coeff.dat')
    if(allocated(idx_2b)) deallocate(idx_2b)
    allocate(idx_2b(6,nw*(nw-1)/2))
    call map_2b(nw,idx_2b)


    return
  end subroutine dms_init

  !==================================================!
  ! water dipole moment                              !
  !==================================================!
  subroutine dipole(x,dp)
    real,dimension(:),intent(in)::x
    real,dimension(3),intent(out)::dp
    ! ::::::::::::::::::::
    real,dimension(3,size(x)/3)::xn
    integer::natm

    natm=size(x)/3
    xn=reshape(x,(/3,natm/))

    call dip12bhbb(natm,xn,dp)

    return
  end subroutine dipole

!==================================================
! dipole moment using HBB DMS fit
!==================================================
subroutine dip12bhbb(natm,xx,dp)
  integer,intent(in)::natm
  real,dimension(3,natm),intent(in)::xx
  real,dimension(1:3),intent(inout)::dp
  !::::::::::::::::::::
  integer,dimension(1:6)::idx
  real,dimension(1:natm)::chg_m,chg_d
  real,dimension(3,6)::x2
  real,dimension(3,3)::x1,xw
  real,dimension(3)::vect,dp2,chgm,dpm
  real,dimension(6)::chg1,chg1p,chg2
  integer::i,j,nw

  nw=natm/3

  !get monomer charges
  chg_m=0.d0
  do i=1,nw
     xw(:,3)=xx(:,2*nw+i)
     xw(:,1:2)=xx(:,2*i-1:2*i)
     call dip_ltp2011(xw,dpm,chgm)
     chg_m(2*nw+i)=chgm(3)
     chg_m(2*i-1:2*i)=chgm(1:2)
   end do


 !get intrinsic atomic charges of each dimer
  chg_d=0.d0
  do i=1,size(idx_2b,2)
     idx=idx_2b(:,i)
     x2(:,(/(j,j=1,6)/))=xx(:,idx)
     call calcdip(dp2,chg2,transpose(x2))
     chg_d(idx)=chg_d(idx)+chg2
  end do

  chg_d=chg_d+chg_m
  dp=matmul(xx,chg_d)


  return
end subroutine dip12bhbb

!==================================================
! Mapping 2body
!==================================================
subroutine map_2b(nw,idx)
  integer,intent(in)::nw
  integer,dimension(:,:),intent(inout)::idx
  !::::::::::::::::::::
  integer::i,j,vec(6),fo,l

  fo=nw*2
  l=1
  do i=1,nw-1
     vec((/1,2,3/))=(/fo+i,i*2-1,i*2/)
     do j=i+1,nw
        vec((/4,5,6/))=(/fo+j,j*2-1,j*2/)
        idx(:,l)=vec
        l=l+1
     end do
  end do

  return
end subroutine map_2b
end module dms_shell

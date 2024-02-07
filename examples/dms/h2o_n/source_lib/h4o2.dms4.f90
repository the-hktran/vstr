        subroutine predip(dname)

        implicit double precision (a-h,o-z)
        implicit integer (i-n)

        character(len=*),intent(in)::dname
        double precision V, cart_in(6,3), v0
        double precision dc0(0:5,0:5),dw0(0:5,0:5), coef(0:2879)

        integer i,j,k,i1,j1,k1

        common/NCOEdip/ms,mr
        common/h4o2coefdip/dc0,dw0,coef

        ms=1827 ; mr=1053

        open(20,file=dname,status='old')
        read(20,*)
        read(20,*)
        read(20,*)(coef(i1),i1=0,ms+mr-1)
        close(20)

        return
        end

!***************************************************************
        subroutine calcdip(dip,chg,cart_in)

        implicit none

        integer ms,mr
        double precision dc0(0:5,0:5),dw0(0:5,0:5),coef(0:2879)

        common/NCOEdip/ms,mr
        common/h4o2coefdip/dc0,dw0,coef
        double precision V, cart_in(6,3), dip(1:3),chg(1:6),chg_ord(0:5)

        double precision cart0(6,3),cart1(6,3)
        integer i,j,k,l,i1,j1,k1,l1,i2,j2,k2,l2

        double precision rvec(0:3),d0(0:5,0:5),r0(0:5,0:5)
        double precision xnuc(0:2,0:5),vec(0:3,0:2879),chgb(0:5,0:2879)

         do j=1,3
          xnuc(j-1,0:1)=cart_in(2:3,j)
          xnuc(j-1,2:3)=cart_in(5:6,j)
          xnuc(j-1,4)=cart_in(1,j)
          xnuc(j-1,5)=cart_in(4,j)
        end do

        call getdvec (ms, mr, xnuc(0:2,0:5), vec,chgb)

        chg_ord=0.d0
        do i=0,5
          chg_ord(i)=dot_product(chgb(i,:),coef)
        end do
     
        chg(2:3)=chg_ord(0:1)
        chg(5:6)=chg_ord(2:3)
        chg(1)=chg_ord(4)
        chg(4)=chg_ord(5)
        
        return
        end subroutine calcdip

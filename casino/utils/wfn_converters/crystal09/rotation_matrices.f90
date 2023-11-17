 SUBROUTINE gsym11(xyv)
  USE numbers
  USE lmaxxx
  USE parinf_module
  USE paral1_module
  USE shell_info
  USE rotmatrix
  IMPLICIT REAL(float) (a-h,o-z)
  IMPLICIT INTEGER (i-n)
  DIMENSION xyv(3,3,48)
  DIMENSION tv(lmax_dft13,lmax_dft13)
  PARAMETER (toll=1e-6_float)
  LOGICAL lpr9
  REAL(float),ALLOCATABLE :: tto0(:)
  mvf=inf(2)
  mxtp=inf(175)
  lpr9=lprint(9)/=0.and.iameq0
  allocate(tto0(shnaosq(mxtp)*mvf),stat=ierr)
  if(ierr/=0)call errvrs(0,'gsym11','tto0 allocation')
  limrot1=0
  i0=0
  do i=0,mxtp
   n1=shnao(i)
   lvalue=i-1
   do mv=1,mvf
    select case (i)
     case (0)
      tv(1,1)=1._float
     case (1)
      tv(1,1)=1._float
      tv(2:4,1)=0._float
      tv(1,2:4)=0._float
      tv(2:4,2:4)=xyv(1:3,1:3,mv)
     case (2)
      tv(1:3,1:3)=xyv(1:3,1:3,mv)
     case default
      call gsym33(xyv(1,1,mv),tv,lvalue)
    endselect
    do j=1,n1
     tto0(i0+1:i0+n1)=tv(1:n1,j)
     i0=i0+n1
    enddo
    if(lpr9)then
     write(iout,"(/' Rotation matrix for symmop',i3,' shell type= ',a)")mv,&
      &shtxt(i)
     call matout(tv,n1,n1)
    endif
   enddo
  enddo
  limrot1=count(abs(tto0)>=toll)
  if(allocated(tto))then
   if(limrot1>size(tto))then
    deallocate(mmom)
    deallocate(mmo)
    deallocate(tto)
   endif
  endif
  if(.not.allocated(tto))then
   allocate(tto(limrot1),stat=ierr)
   if(ierr/=0)call errvrs(0,'gsym11','tto allocation')
   allocate(mmo(limrot1),stat=ierr)
   if(ierr/=0)call errvrs(0,'gsym11','mmo allocation')
   allocate(mmom(limrot1),stat=ierr)
   if(ierr/=0)call errvrs(0,'gsym11','mmom allocation')
   if(lpr9)call errnic(2,limrot1,'gsym11','limrot1')
  endif
  minz(1)=0
  kk=0
  ii=0
  i0=0
  do i=0,mxtp
   n1=shnao(i)
   do mv=1,mvf
    do k=1,n1
     do j=1,n1
      i0=i0+1
      if(abs(tto0(i0))<toll)cycle
      ii=ii+1
      tto(ii)=tto0(i0)
      mmo(ii)=j
      mmom(ii)=k
     enddo
    enddo
    kk=kk+1
    minz(kk+1)=ii
   enddo
  enddo
  deallocate(tto0)
  return
 END SUBROUTINE gsym11


 SUBROUTINE gsym22
  USE numbers
  USE lmaxxx
  USE parinf_module
  USE paral1_module
  USE shell_info
  USE rotmatrix
  IMPLICIT REAL(float) (a-h,o-z)
  IMPLICIT INTEGER (i-n)
  LOGICAL lpr9
  DIMENSION igath(lmax_dft13*lmax_dft13)
  lpr9=lprint(9)/=0.and.iameq0
  mvf=inf(2)
  mxtp=inf(21)
  limrotx=0
  limroty=0
  ibase=0

  do i=0,mxtp
   jbase=0
   do j=0,i
    if (inf(106+i)+inf(106+j).eq.2) then
     do mv=1,mvf
      ii2=minz(ibase+mv+1)
      ii1=minz(ibase+mv)
      limrotx=limrotx+(ii2-ii1)*(minz(jbase+mv+1)-minz(jbase+mv))
      if(i/=j)cycle
      do m=ii1+1,ii2
       do n=ii1+1,ii2
        if (mmo(m).lt.mmo(n)) limroty=limroty+1
       enddo
      enddo
     enddo
    endif
    jbase=jbase+mvf
   enddo
   ibase=ibase+mvf
  enddo

  limrot2=limrotx+limrotx-limroty

  if(allocated(zo34))then
   if(limrot2>size(zo34))then
    deallocate(ko34)
    deallocate(nom34)
    deallocate(no34)
    deallocate(zo34)
   endif
  endif

  if(.not.allocated(zo34))then
   allocate(zo34(limrot2),stat=ierr)
   if(ierr/=0)call errvrs(0,'gsym22','zo34 allocation')
   allocate(no34(limrot2),stat=ierr)
   if(ierr/=0)call errvrs(0,'gsym22','no34 allocation')
   allocate(nom34(limrot2),stat=ierr)
   if(ierr/=0)call errvrs(0,'gsym22','nom34 allocation')
   allocate(ko34(limrot2),stat=ierr)
   if(ierr/=0)call errvrs(0,'gsym22','ko34 allocation')
   if(lpr9)call errnic(2,limrot2,'gsym22','limrot2')
  endif

  ii=0
  j4=0
  do lat4=0,mxtp
   if(inf(106+lat4)==1)then
    id4=shnao(lat4)
    j3=0
    do lat3=0,lat4
     if(inf(106+lat3)==1)then
      id3=shnao(lat3)
      novf34(0,lat3,lat4)=ii
      jjx=ii
      do mv=1,mvf
       m=0
       do i=1,id3
        n=i
        do j=1,id4
         m=m+1
         igath(m)=n
         n=n+id3
        enddo
       enddo
       do i=minz(j3+mv)+1,minz(j3+mv+1)
        mo=(mmo(i)-1)*id4
        mom=(mmom(i)-1)*id4
        do j=minz(j4+mv)+1,minz(j4+mv+1)
         ii=ii+1
         no34(ii)=mmo(j)+mo
         ko34(ii)=igath(no34(ii))
         nom34(ii)=mmom(j)+mom
         zo34(ii)=tto(i)*tto(j)
        enddo
       enddo
       novf34(mv,lat3,lat4)=ii
      enddo
      if(lat3/=lat4)then
       m=ii-jjx
       novf34(0:mvf,lat4,lat3)=novf34(0:mvf,lat3,lat4)+m
       m=ii
       do i=jjx+1,m
        ii=ii+1
        no34(ii)=ko34(i)
        ko34(ii)=no34(i)
        nom34(ii)=igath(nom34(i))
        zo34(ii)=zo34(i)
       enddo
      endif
     endif
     j3=j3+mvf
    enddo
   endif
   j4=j4+mvf
  enddo
  inf(76)=ii
  j4=0
  do lat4=0,mxtp
  if(inf(106+lat4)==1)then
   id4=shnao(lat4)
   lovf34(0,lat4)=ii
   do mv=1,mvf
    do i=minz(j4+mv)+1,minz(j4+mv+1)
     mo3=mmo(i)
     mo=iky(mo3-1)
     mom=(mmom(i)-1)*id4
     do j=minz(j4+mv)+1,minz(j4+mv+1)
      mo4=mmo(j)
      if(mo3>=mo4)then
       ii=ii+1
       no34(ii)=mo+mo4
       nom34(ii)=mmom(j)+mom
       zo34(ii)=tto(i)*tto(j)
      endif
     enddo
    enddo
    lovf34(mv,lat4)=ii
   enddo
  endif
  j4=j4+mvf
  enddo
  return
 END SUBROUTINE gsym22


 SUBROUTINE gsym33(rot,rotx,lmax)
  USE numbers
  USE lmaxxx
  USE sphfac_module
  USE shell_info
  IMPLICIT REAL(float) (a-h,o-z)
  IMPLICIT INTEGER (i-n)
  DIMENSION rot(3,3)
  DIMENSION rotx(lmax_dft13,lmax_dft13)
  DIMENSION dufact(0:lmax_dft6),dubas(lmax_dft91)
  DIMENSION coez(lmax_dft28),coep(lmax_dft28),coem(lmax_dft28),&
   &proj(lmax_dft28,lmax_dft13)
  proj(1:3,1)=rot(3,1:3)
  proj(1:3,2)=rot(1,1:3)
  proj(1:3,3)=rot(2,1:3)
  xx=rot(1,1)
  yx=rot(2,1)
  zx=rot(3,1)
  xy=rot(1,2)
  yy=rot(2,2)
  zy=rot(3,2)
  xz=rot(1,3)
  yz=rot(2,3)
  zz=rot(3,3)
! Recursions in l
  nbase=6
  do l=2,lmax
   bilbo=l
   l1=l-1
! m=0 and m=l cases
   coez(:nbase)=0._float
   coep(:nbase)=0._float
   coem(:nbase)=0._float
   fa=1._float/bilbo
   mpos=l1+l1
   lbase=0
   k=0
   do i=0,l1
    mbase=lbase
    lbase=lbase+i+1
    do j=0,i
     k=k+1
     fb=proj(k,1)*fa
     fc=proj(k,mpos)
     fd=proj(k,mpos+1)
     kbase=mbase+j
     coez(kbase+1)=fb*zx+coez(kbase+1)
     coep(kbase+1)=fc*xx-fd*yx+coep(kbase+1)
     coem(kbase+1)=fd*xx+fc*yx+coem(kbase+1)
     kbase=lbase+j
     coez(kbase+1)=fb*zy+coez(kbase+1)
     coep(kbase+1)=fc*xy-fd*yy+coep(kbase+1)
     coem(kbase+1)=fd*xy+fc*yy+coem(kbase+1)
     coez(kbase+2)=fb*zz+coez(kbase+2)
     coep(kbase+2)=fc*xz-fd*yz+coep(kbase+2)
     coem(kbase+2)=fd*xz+fc*yz+coem(kbase+2)
    enddo
   enddo
   proj(:nbase,1)=coez(:nbase)
   proj(:nbase,mpos+2)=coep(:nbase)
   proj(:nbase,mpos+3)=coem(:nbase)
! m finite but less than l
   mpos=2
   do m=1,l1
    mneg=mpos+1
    coep(:nbase)=0._float
    coem(:nbase)=0._float
    bilbo=bilbo-1._float
    fa=1._float/bilbo
    k=0
    lbase=0
    do i=0,l1
     mbase=lbase
     lbase=lbase+i+1
     do j=0,i
      k=k+1
      fb=proj(k,mpos)*fa
      fc=proj(k,mneg)*fa
      kbase=mbase+j
      coep(kbase+1)=fb*zx+coep(kbase+1)
      coem(kbase+1)=fc*zx+coem(kbase+1)
      kbase=lbase+j
      coep(kbase+1)=fb*zy+coep(kbase+1)
      coem(kbase+1)=fc*zy+coem(kbase+1)
      coep(kbase+2)=fb*zz+coep(kbase+2)
      coem(kbase+2)=fc*zz+coem(kbase+2)
     enddo
    enddo
    proj(:nbase,mpos)=coep(:nbase)
    proj(:nbase,mneg)=coem(:nbase)
    mpos=mpos+2
   enddo
   nbase=nbase+l+2
  enddo
  lterm=lmax+lmax+1
! Generate overlap factors between Cartesian functions
  coez(1)=1._float
  dufact(0)=1._float
  fb=1._float
  hobbit=2._float
  bilbo=lmax
  fc=bilbo
  m=2
  do i=1,lmax
   bilbo=bilbo+1._float
   hobbit=hobbit/(bilbo*fc)
   fc=fc-1._float
   fa=sqrt(hobbit)
   coez(m)=fa
   coez(m+1)=fa
   m=m+2
   dufact(i)=dufact(i-1)*fb
   fb=fb+2._float
  enddo
  dubas(1:iky(lterm))=0._float
  m=0
  do i=0,lmax
   bilbo=dufact(lmax-i)
   do j=0,i
    dubas(m+j+j+1)=dufact(i-j)*dufact(j)*bilbo
   enddo
   m=i*4+m+3
  enddo
! Generate rotation matrix by projection
  m=lmax*lmax
  do mpos=m+1,m+lterm
   k=ietap(mpos)+1
   l=ietap(mpos+1)
   do mneg=1,lterm
    bilbo=0._float
    do n=k,l
     jx=net1(n)
     mbase=net3(n)
     fc=0._float
     nbase=0
     do i=0,lmax
      do j=0,i
       nbase=nbase+1
       fc=dubas(mbase+j)*proj(nbase,mneg)+fc
      enddo
      mbase=mbase+jx+i+1
     enddo
     bilbo=fsph(n)*fc+bilbo
    enddo
    rotx(mneg,mpos-m)=coez(mneg)*bilbo
   enddo
  enddo
 END SUBROUTINE gsym33


 SUBROUTINE alloc_rot12
  USE rotmatrix
  IMPLICIT INTEGER (i-n)
  allocate (mom12(limrot2),stat=ierr)
  if(ierr.ne.0)call errvrs(0,'alloc_rot12','mom12 allocation')
  allocate (no12(limrot2),stat=ierr)
  if(ierr/=0)call errvrs(0,'alloc_rot12','no12 allocation')
 END SUBROUTINE alloc_rot12


 SUBROUTINE dealloc_rot12
  USE rotmatrix
  deallocate (mom12)
  deallocate (no12)
 END SUBROUTINE dealloc_rot12


 SUBROUTINE alloc_rot34
  USE rotmatrix
  IMPLICIT INTEGER (i-n)
  allocate (mo34(limrot2),stat=ierr)
  if(ierr/=0)call errvrs(0,'alloc_rot34','mo34 allocation')
 END SUBROUTINE alloc_rot34


 SUBROUTINE dealloc_rot34
  USE rotmatrix
  deallocate (mo34)
 END SUBROUTINE dealloc_rot34

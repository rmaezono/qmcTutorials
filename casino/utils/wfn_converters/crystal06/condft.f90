 SUBROUTINE condft
  USE numbers
  USE lmaxxx
  USE sphfac_module
  USE shell_info
  IMPLICIT REAL(FLOAT) (A-H,O-Z)
  IMPLICIT INTEGER (i-n)
  COMMON/LOCO/lpol,npo,&
   &apol(lmax_dft28,lmax_dft13,2),ipol(lmax_dft28,lmax_dft13,2),&
   &jpol(lmax_dft28,lmax_dft13,2),mpol(lmax_dft13,2)

  do loop=1,4
   ietap(loop+1)=loop
   mpol(loop,2)=1
   net1(loop)=0
   net2(loop)=0
   ipol(1,loop,2)=0
   jpol(1,loop,2)=0
   fsph(loop)=1._float
   apol(1,loop,2)=1._float
  enddo
  net1(3)=1
  net1(4)=1
  net2(4)=1
  apol(1,1,1)=1._float
  ipol(1,1,1)=0
  jpol(1,1,1)=0
  ietap(1)=0
  lmqu=0
  lpol=4
  lazy=4
  mpol(1,1)=1
  ipol(1,2,2)=1
  jpol(1,2,2)=1
  jpol(1,3,2)=1
  iold=1
  inew=2
  do l=2,lmax_dft6
   zzv=l
   yyv=1._float-zzv
   ufac=zzv
   xx=zzv-yyv
   sfac=1._float
   lpoll=lpol
   mm=1
 !...  m . lt . l  cases
   do m=1,l
    lmqu=lmqu+1
    bilbo=xx/zzv
    frodo=yyv/zzv
    yyv=yyv-1._float
    do mmm=1,2
     npo=0
     if(m.ne.l)then
      do loop=1,mpol(mm,iold)
       zz=apol(loop,mm,iold)*frodo
       i=ipol(loop,mm,iold)
       j=jpol(loop,mm,iold)
       call convrs(i+2,j+2,zz,fsph,net1,net2)
       call convrs(i,j+2,zz,fsph,net1,net2)
       call convrs(i,j,zz,fsph,net1,net2)
      enddo
     endif
     do loop=1,mpol(mm,inew)
      call convrs(ipol(loop,mm,inew),jpol(loop,mm,inew),&
       &apol(loop,mm,inew)*bilbo,fsph,net1,net2)
     enddo
     lazy=lazy+1
     mpol(mm,iold)=npo
     vfac=sqrt(sfac)
     do loop=1,npo
      lpol=lpol+1
      ipol(loop,mm,iold)=net1(lpol)
      jpol(loop,mm,iold)=net2(lpol)
      apol(loop,mm,iold)=fsph(lpol)
      fsph(lpol)=fsph(lpol)*vfac
     enddo
     ietap(lazy+1)=lpol
     mm=mm+1
     if(m.ne.1)cycle
     sfac=2._float
     exit
    enddo
    ufac=ufac+1._float
    sfac=sfac/(zzv*ufac)
    zzv=zzv-1._float
   enddo
   lmqu=lmqu+1
   mmp1=mm+1
   mneg=mm-1
   mpos=mm-2
   mpo=mpol(mpos,inew)
   mne=mpol(mneg,inew)
 !...   m  =  l   case
   npo=0
   do loop=1,mpo
    call convrs(ipol(loop,mpos,inew)+1,jpol(loop,mpos,inew)+1,&
     &apol(loop,mpos,inew)*xx,fsph,net1,net2)
   enddo
   do loop=1,mne
    call convrs(ipol(loop,mneg,inew),jpol(loop,mneg,inew)+1,&
     &-apol(loop,mneg,inew)*xx,fsph,net1,net2)
   enddo
   lazy=lazy+1
   mpol(mm,iold)=npo
   vfac=sqrt(sfac)
   do loop=1,npo
    lpol=lpol+1
    ipol(loop,mm,iold)=net1(lpol)
    jpol(loop,mm,iold)=net2(lpol)
    apol(loop,mm,iold)=fsph(lpol)
    fsph(lpol)=fsph(lpol)*vfac
   enddo
   ietap(lazy+1)=lpol
 !...  -m  =  l   case
   npo=0
   do loop=1,mpo
    call convrs(ipol(loop,mpos,inew),jpol(loop,mpos,inew)+1,&
     &apol(loop,mpos,inew)*xx,fsph,net1,net2)
   enddo
   do loop=1,mne
    call convrs(ipol(loop,mneg,inew)+1,jpol(loop,mneg,inew)+1,&
     &apol(loop,mneg,inew)*xx,fsph,net1,net2)
   enddo
   mpol(mmp1,iold)=npo
   do loop=1,npo
    lpol=lpol+1
    ipol(loop,mmp1,iold)=net1(lpol)
    jpol(loop,mmp1,iold)=net2(lpol)
    apol(loop,mmp1,iold)=fsph(lpol)
    fsph(lpol)=fsph(lpol)*vfac
   enddo
   lazy=lazy+1
   do loop=lpoll+1,lpol
    net1(loop)=l-net1(loop)
    net2(loop)=l-net2(loop)
   enddo
   ietap(lazy+1)=lpol
   loop=iold
   iold=inew
   inew=loop
  enddo
  do i=1,lpol
   net3(i)=iky(net1(i))+net2(i)+1
  enddo
  do i=0,lmax_dft7
   nu3p(i)=i*i
  enddo
  return
 END SUBROUTINE condft

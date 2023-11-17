 program uf2dwfn
!----------------------------------------------------------------------!
! 2DHF to CASINO interface.                                            !
! Converts numerical Hartree-Fock orbitals provided by 2DHF into       !
! the format read by CASINO.                                           !
! Note that 2DHF is only applicable to diatomic molecules, and         !
! that multideterminant expansions (analogous to configuration         !
! interaction) may be constructed by hand by editing the casino input  !
! file.                                                                !
! Output is formatted for casino.                                      !
! J.R. Trail, TCM group, Cavendish Laboratory, Cambridge, UK 2013      !
! usage ~/bin/uf2dwfn                                                  !
! files required: inp.data inp.orb                                     !
! files provided: dwfn.data                                            !
!----------------------------------------------------------------------!
 implicit none
 real(kind=kind(0.d0)) pi
 integer               iin,length1
 integer               i,j,k,k1,k2,imu,in,inioff
 integer               ishift,iorb,igp
 real(kind=kind(0.d0)) vni(600),vmu(900)
 real(kind=kind(0.d0)) hni,hmu
 character(3) isconfig
 character(15) orbinfo
 character(15) orbinfo1(6)
 character(120*7) names_in_order
 integer indx1,indx2,nindx1,nindx2
 integer sp(2,120),type(2,120),mval(2,120),orb_num(2,120)

!---- vars for read v
 character(80)         header,datetime
 integer               ngrids,nni,nmu
 real(kind=kind(0.d0)) r,rgrid,z1,z2
 real(kind=kind(0.d0)) area(60)
 integer               norb,nel,nexch
 integer                i1b(60),i2b(60),i3b(1830)
 integer i1e(60),i2e(60),i3e(1830)
 integer i1si(60),i2si(60),i3si(1830)
 integer i1ng(60),i2ng(60),i3ng(1830)
 integer i1mu(60),i2mu(60),i3mu(1830)
 real(kind=kind(0.d0)),ALLOCATABLE :: cw_orb(:)
!---- vars for read ^

 iin=20
 pi=3.1415926535897931159979634685442d0

!---- read in data about grids, and then the orbitals
 open(iin,file='inp.orb',form='unformatted')

 read(iin) header
 read(iin) datetime
 read(iin) ngrids,nni,nmu
 read(iin) r,rgrid
 read(iin) z1,z2
 read(iin) norb,nel,nexch
 read(iin) i1b,i2b,i3b,i1e,i2e,i3e,          &
  &        i1si,i2si,i3si,i1ng,i2ng,i3ng,    &
  &        i1mu,i2mu,i3mu
 length1=norb*nni*nmu
 allocate(cw_orb(length1))
 write(*,*) 'Reading orbitals ...'
 do i=norb,1,-1          		! Orbitals in reverse order
  read(iin) cw_orb(i1b(i):i1b(i)+i1si(i)-1)
 enddo
 write(*,*) 'done.'
 read(iin) area
 do i=norb,1,-1
  write(*,*) 'area :',i,area(i)
 enddo

 close(iin)

!write(*,*) ngrids,nni,nmu
!write(*,*) r,rgrid
!write(*,*) z1,z2
!write(*,*) norb,nel,nexch
! nni,nmu - grid parameters
! r       - nucleus-nucleus distance
! norb    - number of orbitals
! nel     - number of electrons

!----- make the ni,mu indices
 hmu=rgrid
 vmu(1)=0.d0
 do i=2,nmu
  vmu(i)=vmu(i-1)+hmu
 enddo

 hni=pi/dble(nni-1)
 do i=1,nni
  vni(i)=dble((i-1))*hni
 enddo

 open(iin,file='inp.data')

 read(iin,*)				! title
 read(iin,'(a20)') header(1:20)		! name header
 write(*,*) 'ok.'
 do 					! read until config specs start
  read(iin,'(a3)') isconfig
  if(isconfig=='con')exit
  if(isconfig=='CON')exit
 enddo

 j=0
 sp=0;type=0;mval=0;orb_num=0;indx1=1;indx2=1
 do i=norb,1,-1 		        ! orbitals in reverse order

  if(j>1)then				! duplicate last read in specs
   j=j-1
  else					! ... or not (usual).
   read(iin,'(i2,1x,a15)') j,orbinfo
   k=0 ; k1=1 ; k2=0 ; orbinfo1(:)=' '
   do k1=1,15
    if(orbinfo(k1:k1)==' ')cycle
    if(k1<=k2)cycle
    k=k+1
    do k2=k1,14
      if(orbinfo(k2+1:k2+1)==' ')exit
    enddo
    orbinfo1(k)=orbinfo(k1:k2)
   enddo
   if(trim(orbinfo1(2))/='u'.and.trim(orbinfo1(2))/='g')then
    do k=6,3,-1
      orbinfo1(k)=orbinfo1(k-1)
    enddo
    orbinfo1(2)=' '
   endif
  endif

  names_in_order(1+7*i:7+7*i)=trim(orbinfo1(1))//' '//trim(orbinfo1(2))
  write(*,*) i,j,orbinfo
  if(trim(orbinfo1(1))=='sigma')then
   if(orbinfo1(3)=='+'.or.orbinfo1(4)=='+')then
    sp     (1,indx1)=1
    type   (1,indx1)=0
    mval   (1,indx1)=0
    orb_num(1,indx1)=i
    indx1=indx1+1
   endif
   if(orbinfo1(3)=='-'.or.orbinfo1(4)=='-')then
    sp     (2,indx2)=1
    type   (2,indx2)=0
    mval   (2,indx2)=0
    orb_num(2,indx2)=i
    indx2=indx2+1
   endif
  endif

  if(trim(orbinfo1(1))=='pi')then
   if(orbinfo1(3)=='+'.or.orbinfo1(4)=='+')then
    sp     (1,indx1)=1
    type   (1,indx1)=1
    mval   (1,indx1)=+1
    orb_num(1,indx1)=i
    indx1=indx1+1
   endif
   if(orbinfo1(3)=='-'.or.orbinfo1(4)=='-')then
    sp     (2,indx2)=1
    type   (2,indx2)=1
    mval   (2,indx2)=+1
    orb_num(2,indx2)=i
    indx2=indx2+1
   endif
   if(orbinfo1(5)=='+'.or.orbinfo1(6)=='+')then
    sp     (1,indx1)=1
    type   (1,indx1)=1
    mval   (1,indx1)=-1
    orb_num(1,indx1)=i
    indx1=indx1+1
   endif
   if(orbinfo1(5)=='-'.or.orbinfo1(6)=='-')then
    sp     (2,indx2)=1
    type   (2,indx2)=1
    mval   (2,indx2)=-1
    orb_num(2,indx2)=i
    indx2=indx2+1
   endif
  endif

  if(trim(orbinfo1(1))=='delta')then
   if(orbinfo1(3)=='+'.or.orbinfo1(4)=='+')then
    sp     (1,indx1)=1
    type   (1,indx1)=2
    mval   (1,indx1)=+2
    orb_num(1,indx1)=i
    indx1=indx1+1
   endif
   if(orbinfo1(3)=='-'.or.orbinfo1(4)=='-')then
    sp     (2,indx2)=1
    type   (2,indx2)=2
    mval   (2,indx2)=+2
    orb_num(2,indx2)=i
    indx2=indx2+1
   endif
   if(orbinfo1(5)=='+'.or.orbinfo1(6)=='+')then
    sp     (1,indx1)=1
    type   (1,indx1)=2
    mval   (1,indx1)=-2
    orb_num(1,indx1)=i
    indx1=indx1+1
   endif
   if(orbinfo1(5)=='-'.or.orbinfo1(6)=='-')then
    sp     (2,indx2)=1
    type   (2,indx2)=2
    mval   (2,indx2)=-2
    orb_num(2,indx2)=i
    indx2=indx2+1
   endif
  endif


  if(trim(orbinfo1(1))=='phi')then
   if(orbinfo1(3)=='+'.or.orbinfo1(4)=='+')then
    sp     (1,indx1)=1
    type   (1,indx1)=3
    mval   (1,indx1)=+3
    orb_num(1,indx1)=i
    indx1=indx1+1
   endif
   if(orbinfo1(3)=='-'.or.orbinfo1(4)=='-')then
    sp     (2,indx2)=1
    type   (2,indx2)=3
    mval   (2,indx2)=+3
    orb_num(2,indx2)=i
    indx2=indx2+1
   endif
   if(orbinfo1(5)=='+'.or.orbinfo1(6)=='+')then
    sp     (1,indx1)=1
    type   (1,indx1)=3
    mval   (1,indx1)=-3
    orb_num(1,indx1)=i
    indx1=indx1+1
   endif
   if(orbinfo1(5)=='-'.or.orbinfo1(6)=='-')then
    sp     (2,indx2)=1
    type   (2,indx2)=3
    mval   (2,indx2)=-3
    orb_num(2,indx2)=i
    indx2=indx2+1
   endif
  endif


 enddo
 nindx1=indx1-1 ; nindx2=indx2-1

 close(iin)

 open(iin,file='dwfn.data')
 write(iin,*)           'Tabulated dimer wave functions in real space'
 write(iin,'(a20)')      header(1:20)
 write(iin,*)           'Nuclei-Nuclei distance'
 write(iin,'(f11.6)')    r
 write(iin,*)           'Atomic numbers'
 write(iin,'(2(1x,i3))') nint(z2),nint(z1)				! swap is correct.
 write(iin,*)           'Total number of orbitals'
 write(iin,'(1(1x,i3))') norb
 write(iin,*)           'Number up/downspin electrons, determinants'
 write(iin,'(3(1x,i3))') nindx1,nindx2,1
 write(iin,*)           'States'
 j=1
 do indx1=nindx1,1,-1     		! Orbitals in reverse order
  write(iin,'(4(1x,i2))') j,orb_num(1,indx1),type(1,indx1),mval(1,indx1)
  j=j+1
 enddo
 j=1
 do indx2=nindx2,1,-1     		! Orbitals in reverse order
  write(iin,'(4(1x,i2))') j,orb_num(2,indx2),type(2,indx2),mval(2,indx2)
  j=j+1
 enddo
 write(iin,*)               'Grid, niXmu'
 write(iin,'(2(1x,i4))')     nni,nmu
 write(iin,'(2(1x,e22.15))') hni,hmu

 do iorb=1,norb
  write(iin,'(a14,1x,i3,1x,a10)') 'Orbital number',iorb,                                &
   &     '('//trim(names_in_order(1+7*iorb:7+7*iorb))//')'
   ishift=i1b(iorb)-1
   do imu=1,nmu
     inioff=(imu-1)*nni
     do in=1,nni
      igp=ishift+inioff+in
      write(iin,'(e25.15)') cw_orb(igp)
     enddo
   enddo
 enddo

 close(iin)
 
 end program

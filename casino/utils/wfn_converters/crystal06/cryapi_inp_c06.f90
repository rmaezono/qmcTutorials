!-----------------------------------------------------------------------------!
! Official CRYSTAL06 routines to read the output of CRYAPI_OUT.               !
!                                                                             !
! Changes                                                                     !
! -------                                                                     !
! 9.2008 - MDT Hideous Torino code tidied and formatted properly.             !
! 9.2008 - MDT Stopped it printing out everything it reads.                   !
!-----------------------------------------------------------------------------!


 MODULE numbers
  IMPLICIT NONE
  INTEGER,PARAMETER :: float=selected_real_kind(13,100)
 END MODULE numbers


 MODULE text_module
  INTEGER,PARAMETER,PUBLIC :: nangecp=5
  CHARACTER(len=24) :: nspin(2)
  CHARACTER(len=24) :: mspin(2)
  CHARACTER(len=4)  :: ndn(6)
  CHARACTER(len=6)  :: xmot(nangecp)
  CHARACTER(len=1)  :: tipo(0:1)
  CHARACTER(len=24) :: hftype(0:2)
  CHARACTER(len=23) :: nomc(8),nomx(5)
  CHARACTER(len=2)  :: symbat(0:93)
  DATA nspin/'    ALPHA+BETA ELECTRONS','    ALPHA-BETA ELECTRONS'/
  DATA mspin/'    ALPHA      ELECTRONS','    BETA       ELECTRONS'/
  DATA ndn/' S',' SP',' P',' D',' F',' G'/
  DATA xmot/ 'W0 TMS','P0 TMS','P1 TMS','P2 TMS','P3 TMS'/
  DATA tipo/'C','R'/
  DATA hftype/' RESTRICTED CLOSED SHELL',&
   &' RESTRICTED OPEN SHELL  ',' UNRESTRICTED OPEN SHELL'/
  DATA nomx/'DIRAC-SLATER LDA','VON BARTH-HEDIN','BECKE',&
   &'PERDEW-BURKE-ERNZERHOF','PERDEW-WANG GGA'/
  DATA nomc/'PERDEW-WANG LSD','PERDEW-ZUNGER','VOSKO-WILK-NUSAIR',&
   &'VON BARTH-HEDIN','PERDEW 86','PERDEW-WANG GGA','LEE-YANG-PARR',&
   &'PERDEW-BURKE-ERNZERHOF'/
  DATA symbat/   'XX','H ','HE','LI','BE','B ','C ','N ',&
   &'O ','F ','NE','NA','MG','AL','SI','P ','S ','CL','AR',&
   &'K ','CA','SC','TI','V ','CR','MN','FE','CO','NI','CU',&
   &'ZN','GA','GE','AS','SE','BR','KR','RB','SR','Y ','ZR',&
   &'NB','MO','TC','RU','RH','PD','AG','CD','IN','SN','SB',&
   &'TE','I ','XE','CS','BA','LA','CE','PR','ND','PM','SM',&
   &'EU','GD','TB','DY','HO','ER','TM','YB','LU','HF','TA',&
   &'W ','RE','OS','IR','PT','AU','HG','TL','PB','BI','PO',&
   &'AT','RN','FR','RA','AC','TH','PA','U ','QC'/
 END MODULE text_module


 MODULE parame_module
  USE numbers
  IMPLICIT REAL(float) (a-h,o-z)
  IMPLICIT INTEGER (i-n)
  PARAMETER (lim001=417,lim002=10000,lim005=2400,lim015=3200,&
 &lim016=1200,lim007=3500,lim006=lim007+lim007-1,&
 &lim018=lim002*lim001,lim042=lim002*2,lim086=lim001*4)
 END MODULE parame_module


 MODULE parinf_module
  USE numbers
  IMPLICIT NONE
  PRIVATE
  INTEGER,PARAMETER,PUBLIC :: limftn=100
  INTEGER,PARAMETER,PUBLIC :: limprn=050
  REAL(float),DIMENSION(:),ALLOCATABLE,PUBLIC :: par
  INTEGER,DIMENSION(1:limprn),PUBLIC :: lprint
  INTEGER,DIMENSION(1:limftn),PUBLIC :: iunit
  INTEGER,DIMENSION(:),ALLOCATABLE,PUBLIC :: itol
  INTEGER,DIMENSION(:),ALLOCATABLE,PUBLIC :: inf
  CHARACTER(len=80),PUBLIC :: itit
  INTEGER,PUBLIC :: iin,iout
 END MODULE parinf_module


 MODULE basato_module
  USE numbers
  USE parame_module
  IMPLICIT REAL(float)(a-h,o-z)
  IMPLICIT INTEGER (i-n)
  COMMON/basato/aznuc(lim016),xa(3,lim016),&
   &che(lim015),exad(lim015+1),xl(3,lim015),&
   &exx(lim042),c1(lim042),c2(lim042),c3(lim042),&
   &cmax(lim042),c2w(lim042),c3w(lim042),&
   &nat(lim016),nshpri(lim016+1),ipseud(lim016),&
   &laa(lim015+1),lan(lim015),lat(lim015),latao(lim015),&
   &ndq(lim015+1),latoat(lim015)!iasymmet(lim016),iprim(lim016)
 END MODULE basato_module


 MODULE infpot_module
  USE numbers
  USE parame_module
  INTEGER,PARAMETER :: nang=5
  REAL(float),DIMENSION(lim042) :: apot,cpot
  INTEGER,DIMENSION(lim042) :: npot
  INTEGER,DIMENSION(lim016*nang) :: nbtyp
  INTEGER,DIMENSION(lim016+1) :: nsom
 END MODULE infpot_module


 MODULE gvect_module
  USE numbers
  USE parame_module
  REAL(float) :: paret(3,3),w1r(3,3)
  REAL(float),DIMENSION(lim007+1) :: gmodus
  REAL(float),DIMENSION(3,lim006) :: xg
  INTEGER,DIMENSION(lim007+1) :: nm,mn
  INTEGER,DIMENSION(lim006) :: nn1
  INTEGER,DIMENSION(3,lim006) :: lg
 END MODULE gvect_module


 MODULE retic_module
  USE numbers
  USE parame_module
  REAL(float),DIMENSION(3,3),PUBLIC :: bret
  REAL(float),DIMENSION(lim001),PUBLIC :: wpj
  REAL(float),DIMENSION(lim086),PUBLIC :: cossma,sinsma
  INTEGER,DIMENSION(lim001),PUBLIC :: latvrs
  INTEGER,DIMENSION(3,lim001),PUBLIC :: jj
  INTEGER,DIMENSION(3,3,48),PUBLIC :: irr
  INTEGER, PUBLIC :: is,isp,nkf,nkif,is1,is2,is3
 END MODULE retic_module


 MODULE epesi_module
  USE numbers
  USE parame_module
  REAL(float),DIMENSION(lim018) :: ene,alfa
 END MODULE epesi_module


 MODULE molbar_module
  USE numbers
  USE parame_module
  REAL(float),DIMENSION(3,lim016),PUBLIC :: xbar
  INTEGER,DIMENSION(lim016),PUBLIC :: n1mol
 END MODULE molbar_module


 MODULE xyvdim_module
  USE numbers
  USE parame_module
  REAL(float),DIMENSION(3,3,48),PUBLIC :: xyv
  REAL(float),DIMENSION(3,48),PUBLIC :: trasv
  INTEGER,DIMENSION(48),PUBLIC :: ninv
  INTEGER,DIMENSION(48,48),PUBLIC :: multab
 END MODULE xyvdim_module


 SUBROUTINE cryread
  USE numbers
  USE parame_module
  USE parinf_module
  USE basato_module
  USE infpot_module
  USE gvect_module
  USE xyvdim_module
  USE text_module
  USE molbar_module
  USE memory_screen
  IMPLICIT REAL(float) (a-h,o-z)
  IMPLICIT INTEGER (i-n)
  LOGICAL iprat(0:399),exist
  CHARACTER(len=7) :: znamz='CRYREAD'
  INTEGER,DIMENSION(:),ALLOCATABLE :: iky
  INTEGER,DIMENSION(:),ALLOCATABLE :: ncf,inzvlb,idime,idimf,idmcou
  INTEGER,DIMENSION(:),ALLOCATABLE :: la3,la4,jncdu,irof,jrof
  INTEGER,DIMENSION(:),ALLOCATABLE :: iccat,iccs3,kvrsp,iccv,icct,icc
  INTEGER,DIMENSION(:),ALLOCATABLE :: ina12,nlana
  INTEGER,DIMENSION(:),ALLOCATABLE :: ila12t,jpoint,iccs1,nla21t
  INTEGER,DIMENSION(:),ALLOCATABLE :: nnnc,nnnc2,la34x,la34v,ilana
  INTEGER,DIMENSION(:),ALLOCATABLE :: nshg,ngshg,nqgshg,nstatg
  INTEGER,DIMENSION(:),ALLOCATABLE :: nngi,nshgi
  REAL(float),DIMENSION(:),ALLOCATABLE :: qtot,sg,fg,pg,a
!------------------------------------!
! Suppress printing by default (MDT) !
  LOGICAL :: PRINT=.false.
!------------------------------------!
  iin=5
  iout=6
  iuni=97

  inquire(file='GRED.DAT',exist=exist)
  if(.not.exist)then
   call errnic(0,iuni,znamz,'File GRED.DAT not found - no wave function data.')
  else
   open(iuni,file='GRED.DAT',form='formatted',status='old')
  endif
  do i=1,limftn
   iunit(i)=i
  enddo
  read(iuni,107)itit

  inquire(file='fort.98',exist=exist)
  if(exist)then
   open(98,file='fort.98',form='formatted',status='old')
   read(98,'(a80)')itit
   close(98)
  endif

  read(iuni,106)luminf,lumtol,lumpar
  allocate(inf(1:luminf),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of inf')
  allocate(itol(1:lumtol),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of itol')
  allocate(par(1:lumpar),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of par')
  call ireadf(iuni,inf,luminf)
  call ireadf(iuni,itol,lumtol)
  call rreadf(iuni,par,lumpar)

  if(PRINT)then
   write(iout,203)itit,hftype(inf(27)+inf(64))
   if(inf(170)/=0)then
    write(iout,1007)
    if((inf(172)/=0).and.(inf(171)/=0))then
     write(iout,9203)trim(nomx(inf(172))),trim(nomc(inf(171)))
     if(abs(par(39)-1._float)>=1e-11_float)write(iout,2011)par(39)
     if(abs(par(38)-1._float)>=1e-11_float)write(iout,2012)par(38)
    else
     if(inf(172)/=0)then
      write(iout,9201)trim(nomx(inf(172)))
      if(abs(par(39)-1._float)>=1e-11_float)write(iout,2011)par(39)
     endif
     if(inf(171)/=0)then
      write(iout,9202)trim(nomc(inf(171)))
      if(abs(par(38)-1._float)>=1e-11_float)write(iout,2012)par(38)
     endif
    endif
    if(inf(173)/=0)write(iout,1445)par(48)
   else
    write(iout,1008)
   endif
  endif ! PRINT

1007 format(' Kohn-Sham Hamiltonian')
1008 format(' Hartree-Fock Hamiltonian')
9201 format(/' The exchange functional    ',a,' is active')
9202 format(/' The correlation functional ',a,' is active')
9203 format(/' (Exchange)[correlation] functional:(',a,')[',a,']'/)
1445 format(/' Hybrid exchange - percetage of Fock exchangeE =',t48,f10.4)
2011 format(' Non-local weighting factor (exchange) =',t48,f10.4)
2012 format(' Non-local weighting factor [correlation] =',t48,f10.4)

! Check of dimensions of static arrays
  laf=inf(20)
  if(inf(20)>lim015)call errnic(0,inf(20),znamz,&
   &'Too many shells - increase lim015 to')
  laf1=laf+1
  laf3=3*laf
  naf=inf(24)
  if(inf(24)>lim016)call errnic(0,inf(24),znamz,&
   &'Too many atoms - increase lim016 to')
  naf1=naf+1
  naf3=3*naf
  ndf=inf(7)
  if(inf(7)>lim002)call errnic(0,inf(7),znamz,&
   &'Too many ao - increase lim002 to')
  nprim=inf(75)
  if(inf(75)>lim042)call errnic(0,inf(75),znamz,&
   &'Too many primitives - increase lim042 to')
  mvlaf=inf(56)
  la34f=inf(73)
  nspsta=inf(64)+1
  allocate(qtot(1:naf),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of qtot')

! COMMON xyv,gvect
! Cartesian components of lattice parameters
  call rreadf(iuni,paret,9)
! Transformation matrix crystallographic <=> primitive cell
  call rreadf(iuni,w1r,9)
! Inverse of symmetry operators
  call ireadf(iuni,ninv,inf(2))
! Multiplication table
  call ireadf(iuni,multab,48*48)
! Rotational part of symmetry operators in cartesian reference
  call rreadf(iuni,xyv,9*inf(2))
! Translational part of symmetry operators
  call rreadf(iuni,trasv,3*inf(2))
! Print of cartesian components of lattice parameters
  if(inf(10)/=0.and.PRINT)write(iout,1130)((paret(i,j),j=1,3),i=1,3)
1130 format(/' Structure and symmetry information '//&
   &' Direct lattice vector components (bohr)'//&
   &t13,'x', t28,'y',t43,'z'/&
   &' b1',3f15.9/' b2',3f15.9/' b3',3f15.9/)
! Print of symmetry operators
  mvf=inf(2)
  if(PRINT)then
   write(iout,200)
   do mv1=1,mvf,2
    mv9=min(mv1+1,mvf)
    write(iout,201)(i,ninv(i),i=mv1,mv9)
    do i=1,3
     write(iout,202)((xyv(i,j,iv),j=1,3),trasv(i,iv),iv=mv1,mv9)
    enddo
   enddo
  endif ! PRINT

! Print of multiplication table
  if(PRINT)then
   write(iout,100)
   call matint_d(multab,mvf,48)
  endif ! PRINT
100 format(/' Point group multiplication table')
200 format(/' ***** Symmops - translators in bohr')
201 format(/2('   No.',i3,' Inverse',i3,20x))
202 format(2(1x,3f8.3,2x,f8.3,5x))
  inf5=inf(5)+1
  inf79=inf(79)
  inf793=inf79*3
! Information on direct lattice - number of direct lattice vectors
  if(inf79>lim006)call errnic(0,inf79,znamz,&
   &'Too many g - increase lim006 to')
! Information on direct lattice - number of stars of direct lattice vectors
  if(inf(5)>lim007)call errnic(0,inf(5),znamz,'Too many g - increase lim007 to')
! Square of modulus of direct lattice stars of vectors
  call rreadf(iuni,gmodus,inf5)
! Cartesian coordinates of direct lattice vectors
  call rreadf(iuni,xg,inf793)
  call ireadf(iuni,nm,inf5)
! Number of vectors in each star
  call ireadf(iuni,mn,inf5)
! Inverse vector
  call ireadf(iuni,nn1,inf79)
! Coordinates of lattice vectors in crystallographic units
  call ireadf(iuni,lg,inf793)
! Print the stars of direct lattice vectors
  nstar=inf(5)
! Number of stars of direct lattice vectors to be printed
  lprint(1)=3
  if(inf(10)/=0.and.PRINT)then
   write(iout,101)inf(79),inf(5),sqrt(gmodus(nstar))
   jmax=min(lprint(1),nstar)
   do 14 lsh=1,jmax
    n=mn(lsh)
    j=mn(lsh+1)
    m=j-n
    r=sqrt(gmodus(lsh))
    write(iout,102)lsh,m,r
    do m=n,j-1
     write(iout,103)m,nn1(m),(lg(i,m),i=1,3),(xg(i,m),i=1,3)
    enddo
14 continue
  endif
101 format(/' No.of g vectors',i5,' stars',i5,' rmax',f11.5,' au')
102 format(/' Star n. ',i3,' N. of vectors ',i3,' r=',1pe15.7/&
     &'   lg i(lg) lx   ly   lz         x           y          z')
103 format(5i5,4x,3f12.6,2i6)
204 format(/2x,a24)
! End printing direct lattice vectors

! COMMON basato
! Nuclear charge
  call rreadf(iuni,aznuc,naf)
! Cartesian coordinates of atoms in the reference cell
  call rreadf(iuni,xa,naf3)
! Formal charge attributed to each shell
  call rreadf(iuni,che,laf)
! Adjoined gaussian of each shell
  call rreadf(iuni,exad,laf)
! Cartesian components of shell coordinates
  call rreadf(iuni,xl,laf3)
! Primitive gaussian exponents
  call rreadf(iuni,exx,nprim)
! Contraction coefficients of gaussians in crystal normalization
  call rreadf(iuni,c1,nprim)
  call rreadf(iuni,c2,nprim)
  call rreadf(iuni,c3,nprim)
  call rreadf(iuni,cmax,nprim)
! Contraction coefficients of gaussians in old normalization
  call rreadf(iuni,c2w,nprim)
  call rreadf(iuni,c3w,nprim)
! Formal atomic number of atoms (z=mod(nat,100))
  call ireadf(iuni,nat,naf)
! First shell of each atom
  call ireadf(iuni,nshpri,naf1)
! First primitive of each shell
  call ireadf(iuni,laa,laf1)
! Number of primitives in each shell
  call ireadf(iuni,lan,laf)
! Shell type (0=s, 1=sp, 2=p, 3=d, 4=f)
  call ireadf(iuni,lat,laf)
! Number of ao per shell
  call ireadf(iuni,latao,laf)
! First ao of each shell
  call ireadf(iuni,ndq,laf1)
! Atoms to which each shell belongs
  call ireadf(iuni,latoat,laf)
  if(PRINT)then
   write(iout,1072)
   do na=1,inf(24)
    iam=nat(na)
    ia=mod(iam,100)
    write(iout,1082)na,iam,symbat(ia),(xa(k,na),k=1,3)
   enddo
   write(iout,1071)
1082 format(i4,2x,i4,1x,a,1p,3e20.12)
1072 format(/' Atom cartesian coordinates (bohr) - primitive cell'/&
     &1x,79('*')/' *      atom',t22,'x',t42,'y',t62,'z'/1x,79('*'))
  endif ! PRINT
! Printing of basis set
  if(PRINT)write(iout,1069)
1069 format(' Variational basis set'/)
  iprat(0:399)=.false.
  if(PRINT)then
   write(iout,1070)
   do na=1,inf(24)
    iam=nat(na)
    ia=mod(iam,100)
    write(iout,1080)na,symbat(ia),(xa(k,na),k=1,3)
    if(iprat(iam))cycle
    iprat(iam)=.true.
    do i=nshpri(na),nshpri(na+1)-1
     lati=lat(i)+1
     k=ndq(i)+1
     if(lati==1)then
      write(iout,1090)k,ndn(lati)
     else
      write(iout,1100)k,ndq(i+1),ndn(lati)
     endif
     do k=laa(i),laa(i+1)-1
      write(iout,1120)exx(k),c1(k),c2(k),c3(k)
     enddo
    enddo
   enddo
   write(iout,1071)
  endif ! PRINT
1080 format(i4,1x,a,3f7.3)
1070 format(1x,79('*')/' Local atomic functions basis set'/1x,79('*')/&
      &'   atom  x(au)  y(au)  z(au)    no. type  exponent',&
      &'  s coef   p coef   d/f/g coef'/1x,79('*'))
1090 format(31x,i4,a)
1100 format(26x,i4,'-',i4,a)
1120 format(40x,1p,4e10.3)
1071 format(1x,79('*')/)

! COMMON infpot - information on ECP
  if(inf(31)/=0)then
   call ireadf(iuni,ipseud,naf)
   read(iuni,108)imax,jmax,itypse,(npot(i),i=1,imax),(nbtyp(i),i=1,jmax),&
    &(nsom(i),i=1,itypse)
   read(iuni,109)(apot(i),i=1,imax),(cpot(i),i=1,imax)
! ECP printing
   if(PRINT)then
    write(iout,1111)
    jmax=0
    itypse=0
    do na=1,naf
     if(ipseud(na)==0)cycle
     jmax=jmax+nangecp
     itypse=itypse+1
     write(iout,8730)mod(nat(na),100),aznuc(na)
     ir=jmax-nangecp
     k=nsom(itypse)
     do j=1,nangecp
      ir=ir+1
      ngau=nbtyp(ir)
      if(ngau==0)cycle
      nb=min(ngau,2)
      write(iout,8750)xmot(j),(apot(i+k),cpot(i+k),npot(i+k),i=1,nb)
      do l=nb+1,ngau,2
       lb=l+1
       if(l==ngau)lb=ngau
       write(iout,8760)(apot(i+k),cpot(i+k),npot(i+k),i=l,lb)
      enddo
      k=k+ngau
     enddo
    enddo
   endif ! PRINT
  endif ! inf(31)
1111 format(/1x,79('*')/' *** Pseudopotential information ***'/&
    &1x,79('*'))
8730  format(/' Atomic number',i4,', Nuclear charge',f7.3/&
     &' Pseudopotential'/'   Type       Exponent      Coeff.    N',&
     &'      Exponent      Coeff.    N')
 8750 format(2x,a6,4(f14.7,f13.7,i4))
 8760 format(8x,4(f14.7,f13.7,i4))

! COMMON molbar - for molsplit only
  if(inf(34)/=0)then
   read(iuni,108)(n1mol(i),i=1,inf(92)+1)
   call rreadf(iuni,xbar,inf(92)*3)
   if(PRINT)then
    write(iout,110)inf(92)
    write(iout,111)(j,(xbar(i,j),i=1,3),j=1,inf(92))
    write(iout,112)(i,n1mol(i),i=1,inf(92))
   endif
111 format(/' Coordinates of the barycentres of the molecules'//&
   &(i6,3e20.12))
110 format(/' Lattice of ',i3,' molecules not interacting')
112 format(/' First atom of the molecule:'/(8(i5,i4,',')))
  endif
! Electronic charge according to Mulliken population analysis
  call rreadf(iuni,qtot,naf)
  if(PRINT)write(iout,1079)
1079 format(/' Mulliken atomic charges'/)
1081 format(2(i4,1x,a,' znuc',f7.3,' qtot',f7.3))
  if(PRINT)write(iout,1081)(na,symbat(mod(nat(na),100)),aznuc(na),qtot(na),&
   &na=1,naf)
106 format(8i10)
107 format(a80)
203 format(/13(' CRY06')//' Information from CRYSTAL06'//1x,80a1//&
   &' Type of calculation : ',a24/)
108 format(/(8i10))
109 format(/(4e20.13))

! COMMON lavlaf
  mvlaf1=mvlaf+1
  laf1=laf+1
  allocate(ncf(1:mvlaf1),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of ncf')
  allocate(inzvlb(1:mvlaf),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of inzvlb')
  allocate(idime(1:mvlaf),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of idime')
  allocate(idimf(1:mvlaf),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of idimf')
  allocate(idmcou(1:mvlaf),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of idmcou')
! Symmetry information - for programmers only
  call ireadf(iuni,ncf,mvlaf1)
  inf793=ncf(mvlaf1)
  call ireadf(iuni,inzvlb,mvlaf)
  call ireadf(iuni,idime,mvlaf)
  call ireadf(iuni,idimf,mvlaf)
  call ireadf(iuni,idmcou,mvlaf)
  allocate(jncdu(1:mvlaf),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of jncdu')
  allocate(irof(1:mvlaf),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of irof')
  allocate(jrof(1:mvlaf),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of jrof')
  allocate(lav(1:laf,1:mvf),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of lav')
  allocate(mgnav(1:naf,1:mvf),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of mgnav')
  allocate(la3(1:inf793),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of la3')
  allocate(la4(1:inf793),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of la4')
  call ireadf(iuni,lav,laf*mvf)
  call ireadf(iuni,mgnav,naf*mvf)
  call ireadf(iuni,la3,inf793)
  call ireadf(iuni,la4,inf793)
  call ireadf(iuni,jncdu,mvlaf)
  call ireadf(iuni,irof,mvlaf)
  call ireadf(iuni,jrof,mvlaf)
  read(iuni,*)isizeiccs3,isizeiccs1,isizeilana
  allocate(iccat(1:naf1),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of iccat')
  allocate(iccs3(1:isizeiccs3),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of iccs3')
  allocate(kvrsp(1:laf*laf),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of kvrsp')
  allocate(iccv(1:laf),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of iccv')
  allocate(icct(1:laf1),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of icct')
  allocate(icc(1:laf),stat=ierror)
  if(ierror/=0)call errnic(0,ierror,znamz,'Allocation of icc')
  call ireadf(iuni,iccat,naf1)
  call ireadf(iuni,iccs3,isizeiccs3)
  call ireadf(iuni,kvrsp,laf*laf)
  call ireadf(iuni,iccv,laf)
  icclaf=iccat(naf1)
  allocate(ina12(1:icclaf-1),stat=ierror)
  if(ierror/=0)call errnic(0,icclaf-1,znamz,'Allocation of ina12')
  allocate(nlana(1:icclaf),stat=ierror)
  if(ierror/=0)call errnic(0,icclaf,znamz,'Allocation of nlana')
  call ireadf(iuni,ina12,icclaf-1)
  call ireadf(iuni,nlana,icclaf)
  call ireadf(iuni,icc,laf)
  call ireadf(iuni,icct,laf1)
  icclaf=icct(laf1)
  allocate(ila12t(1:icclaf),stat=ierror)
  if(ierror/=0)call errnic(0,icclaf,znamz,'Allocation of ila12t')
  allocate(jpoint(1:icclaf),stat=ierror)
  if(ierror/=0)call errnic(0,icclaf,znamz,'Allocation of jpoint')
  allocate(iccs1(1:isizeiccs1),stat=ierror)
  if(ierror/=0)call errnic(0,isizeiccs1,znamz,'Allocation of iccs1')
  allocate(nla21t(1:icclaf),stat=ierror)
  if(ierror/=0)call errnic(0,icclaf,znamz,'Allocation of nla21t')
  allocate(nnnc(1:icclaf),stat=ierror)
  if(ierror/=0)call errnic(0,icclaf,znamz,'Allocation of nnnc')
  allocate(nnnc2(1:icclaf),stat=ierror)
  if(ierror/=0)call errnic(0,icclaf,znamz,'Allocation of nnnc2')
  allocate(la34x(1:icclaf),stat=ierror)
  if(ierror/=0)call errnic(0,icclaf,znamz,'Allocation of la34x')
  allocate(la34v(1:icclaf),stat=ierror)
  if(ierror/=0)call errnic(0,icclaf,znamz,'allocation of la34v')
  allocate(ilana(1:isizeilana),stat=ierror)
  if(ierror/=0)call errnic(0,isizeilana,znamz,'Allocation of ilana')
  call ireadf(iuni,ila12t,icclaf)
  call ireadf(iuni,jpoint,icclaf)
  call ireadf(iuni,iccs1,isizeiccs1)
  call ireadf(iuni,nla21t,icclaf)
  call ireadf(iuni,nnnc,icclaf)
  call ireadf(iuni,nnnc2,icclaf)
  call ireadf(iuni,la34x,icclaf)
  call ireadf(iuni,la34v,icclaf)
  call ireadf(iuni,ilana,isizeilana)
  nngidmf=inf(39)*inf(133)
  allocate(nshg(1:la34f),stat=ierror)
  if(ierror/=0)call errnic(0,la34f,znamz,'Allocation of nshg')
  allocate(ngshg(1:la34f*inf(37)),stat=ierror)
  if(ierror/=0)call errnic(0,la34f*inf(37),znamz,'Allocation of ngshg')
  allocate(nqgshg(1:la34f*inf(145)),stat=ierror)
  if(ierror/=0)call errnic(0,la34f*inf(145),znamz,'Allocation of nqgshg')
  allocate(nstatg(1:mvlaf1),stat=ierror)
  if(ierror/=0)call errnic(0,mvlaf1,znamz,'Allocation of nqgshg')
  allocate(nngi(1:nngidmf),stat=ierror)
  if(ierror/=0)call errnic(0,nngidmf,znamz,'Allocation of nngi')
  allocate(nshgi(1:nngidmf),stat=ierror)
  if(ierror/=0)call errnic(0,nngidmf,znamz,'Allocation of nshgi')
  call ireadf(iuni,nshg,la34f)
  call ireadf(iuni,ngshg,la34f*inf(37))
  call ireadf(iuni,nqgshg,la34f*inf(145))

! COMMON inoetc
  call ireadf(iuni,nstatg,mvlaf1)
  call ireadf(iuni,nngi,nngidmf)
  call ireadf(iuni,nshgi,nngidmf)
  npgt=inf(12)
  nfgt=inf(11)
  allocate(sg(1:npgt),stat=ierror)
  if(ierror/=0)call errnic(0,npgt,znamz,'Allocation of sg')
  allocate(fg(1:nfgt*nspsta),stat=ierror)
  if(ierror/=0)call errnic(0,nfgt,znamz,'Allocation of fg')
  allocate(pg(1:npgt*nspsta),stat=ierror)
  if(ierror/=0)call errnic(0,npgt,znamz,'Allocation of pg')
  allocate(a(1:ndf*ndf),stat=ierror)
  if(ierror/=0)call errnic(0,ndf*ndf,znamz,'Allocation of a')
! IKY computed to print Hamiltonian matrix
  allocate(iky(1:ndf+1),stat=ierror)
  if(ierror/=0)call errnic(0,ndf,znamz,'Allocation of iky')
! fort.3 overlap matrix

  call rreadf(iuni,sg,npgt) !!!!!!!!!!!!!!!!!!!!!! CRYSTAL06 version
!  call rreadf(iuni,sg,nfgt)  !!!!!!!!!!!!!!!!!!!!!! CRYSTAL09 version

! Unrestricted open shell:
! Fock matrix: 2 arrays, alpha and beta
! fort.11 fock/ks matrix
  call rreadf(iuni,fg,nfgt)
  if(nspsta/=1)call rreadf(iuni,fg(nfgt+1),nfgt)
! Unrestricted open shell:
! Density matrix: 1 array, alpha+beta, alfa-beta (total, spin)
! fort.13 density matrix
  call rreadf(iuni,pg,npgt*nspsta)

! Printing of density matrix
  if(PRINT)then

   write(iout,180)
180 format(/' Density matrix direct lattice - g=0'/)
   mgmg=1
! MATPMN - density matrix direct space - mgmg direct lattice vector
! Modify the value of mgmg to print the density matrix relative to the
! first mgmg direct lattice vectors.
   ndf=inf(7)
   laf=inf(20)
   nstap=0
   do isigma=1,nspsta
    write(iout,204)nspin(isigma)
    a(1:ndf*ndf)=0._float
    k=0
    do l1=1,laf
     la1=latao(l1)
     do 123 ll2=icct(l1)+1,icct(l1+1)
      l2=ila12t(ll2)
      ksh=idime(jpoint(ll2))
      if(ksh==0)goto 123
      ippop=la34v(ll2)
      inz=ngshg(ippop+1)+1
      do mg=inz,ngshg(ippop+ksh+1)
       if(nqgshg(mg)/=mgmg)cycle
       la2=latao(l2)
       ii=ndq(l2)+k
       m=(mg-inz)*la1*la2+nnnc(ll2)
       do i=1,la1
        do j=1,la2
         m=m+1
         a(ii+j)=pg(m+nstap)
        enddo
        ii=ii+ndf
       enddo
       goto 123
      enddo
123  continue
     k=la1*ndf+k
    enddo

    do m=1,ndf,10
     k=min(m+9,ndf)
     write(iout,40) (j,j=m,k)
     k=k-m
     inz1=m
     do i=1,ndf
      inz2=inz1+k
      write(iout,50)i,(a(j),j=inz1,inz2)
      inz1=inz1+ndf
     enddo
    enddo
    nstap=npgt
   enddo

  endif ! PRINT

50 format(i4,3x,1p,10e12.4)
40 format(/3x,10(8x,i4)/)
  mgmg=1

! MATFG ini - Hamiltonian matrix in direct space - mgmg direct lattice vector
! Modify the value of mgmg to print the density matrix relative to the
! first mgmg direct lattice vectors.
  if(PRINT)then
   write(iout,181)
181 format(/' Hamiltonian matrix direct lattice - g=0')
   ndf=inf(7)
   laf=inf(20)
   nstap=0
   do isigma=1,nspsta
    if(nspsta/=1)write(iout,204)mspin(isigma)
    ii=0
    do i=1,ndf+1
     iky(i)=ii
     ii=ii+i
    enddo
    a(1:iky(ndf+1))=0._float
    do l1=1,laf
     la1=latao(l1)
     inz1=ndq(l1)+1
     ifn1=ndq(l1+1)
     do l2=l1,laf
      ll2=kvrsp(iccv(l1)+l2)
      if(ll2==0)cycle
      ksh=idimf(jpoint(ll2))
      if(ksh==0)cycle
      ippop=iccs3(ll2)
      mstart=nnnc2(ll2)
      inz=ngshg(ippop+1)+1
      if(l1==l2)then
       mmult=iky(la1+1)
      else
       la2=latao(l2)
       inz2=ndq(l2)
       mmult=la1*la2
      endif
      do mg=inz,ngshg(ippop+ksh+1)
       if(nqgshg(mg)/=mgmg)cycle
       m=(mg-inz)*mmult+mstart
       if(l1==l2)then
        do i=inz1,ifn1
         ii=iky(i)
         do j=inz1,i
          m=m+1
          a(ii+j)=fg(m)
         enddo
        enddo
       else
        do i=1,la2
         ii=iky(i+inz2)
         mm=m+i
         do j=inz1,ifn1
          a(ii+j)=fg(mm+nstap)
          mm=mm+la2
         enddo
        enddo
       endif
       exit
      enddo
     enddo
    enddo
    ii=1
    do m=1,ndf,10
     i=min(m+9,ndf)
     write(iout,40)(j,j=m,i)
     do i=m,ndf
      inz1=iky(i)+ii
      inz2=inz1+min(i-m,9)
      write(iout,50)i,(a(j),j=inz1,inz2)
     enddo
     ii=ii+10
    enddo
    nstap=nfgt
   enddo
! MATFG end
  endif ! PRINT

  if(inf(139)/=0)then
   if(PRINT)write(iout,1068)
1068 format(/' Localization data')
   do isigma=1,nspsta
    if(nspsta/=1)write(iout,204)mspin(isigma)
    call readloc(ndf,naf,iuni,PRINT)
   enddo
  endif
  close(iuni)
  deallocate(iky)
  deallocate(a)
  deallocate(pg)
  deallocate(fg)
  deallocate(sg)
  deallocate(nshgi)
  deallocate(nngi)
  deallocate(nstatg)
  deallocate(nqgshg)
  deallocate(ngshg)
  deallocate(nshg)
  deallocate(ilana)
  deallocate(la34v)
  deallocate(la34x)
  deallocate(nnnc2)
  deallocate(nnnc)
  deallocate(nla21t)
  deallocate(iccs1)
  deallocate(jpoint)
  deallocate(ila12t)
  deallocate(nlana)
  deallocate(ina12)
  deallocate(icc)
  deallocate(icct)
  deallocate(iccv)
  deallocate(kvrsp)
  deallocate(iccs3)
  deallocate(iccat)
  deallocate(la4)
  deallocate(la3)
  deallocate(jrof)
  deallocate(irof)
  deallocate(jncdu)
  deallocate(idmcou)
  deallocate(idimf)
  deallocate(idime)
  deallocate(inzvlb)
  deallocate(ncf)
  deallocate(qtot)
 END SUBROUTINE cryread


 SUBROUTINE readloc(ndf,naf,iu,PRINT)
  USE numbers
  USE parame_module
  USE parinf_module
  USE gvect_module
  IMPLICIT REAL(float)(a-h,o-z)
  IMPLICIT INTEGER (i-n)
  REAL(float),DIMENSION(:,:,:),ALLOCATABLE :: wfun
  REAL(float),DIMENSION(:,:),ALLOCATABLE :: wcent
  REAL(float),DIMENSION(:),ALLOCATABLE :: poploc
  INTEGER,DIMENSION(:),ALLOCATABLE :: lbands
  LOGICAL PRINT
  CHARACTER(len=7) :: znamz='readloc'

  read(iu,106)mmgg,mmmg0,limcll,nocc
106 format(8i10)
! mmgg, mmmg0 used by external program
  allocate(lbands(1:nocc),stat=ierror)
  if(ierror/=0)call errnic(0,nocc,znamz,'Memory allocation for lbands - failed')
  call ireadf(iu,lbands,nocc)
  if(PRINT)write(iout,1005)(lbands(i),i=1,nocc)
1005 format(/' Wannier functions - list of active bands'/(10i7/))

! Read Wannier functions
  allocate(wfun(1:ndf,1:nocc,1:mmmg0),stat=ierror)
  if(ierror/=0)call errnic(0,mmmg0*nocc*ndf,znamz,'Memory allocation for wfun &
   &- failed')
  call rreadf(iu,wfun,mmmg0*nocc*ndf)

! Print wf for the first 3 lattice vectors
  if(PRINT)then
   do mg=1,min(3,mmmg0)
    write(iout,1000)(lg(i,mg),i=1,3)
    call cmacol(wfun(1,1,mg),ndf,1,nocc,1)
   enddo
  endif
1000 format(/'g = (',3i4,')')

! Read wfs centroids
  allocate(wcent(1:3,1:nocc),stat=ierror)
  if(ierror/=0)call errnic(0,nocc*3,znamz,'Memory allocation for wcent - &
   &failed')
  call rreadf(iu,wcent,nocc*3)
  if(PRINT)then
   write(iout,1004)
   do ic=1,nocc
    write(iout,1006)ic,(wcent(n,ic),n=1,3)
   enddo
  endif
1004 format(//' Wannier',t26,'centroid''s coordinates (bohr):'/)
1006 format(i6,3e20.12)

! Read atomic population from unit 80
  lnm=naf*nocc*mmmg0
  allocate(poploc(lnm),stat=ierror)
  if(ierror/=0)call errnic(0,lnm,znamz,'Memory allocation for poploc - failed')
  call rreadf(iu,poploc,lnm)
  deallocate(poploc)
  deallocate(wcent)
  deallocate(wfun)
  deallocate(lbands)

 END SUBROUTINE readloc


 SUBROUTINE cmacol(a,n,niniz,nifin,latv)
  USE numbers
  USE parame_module
  USE parinf_module
  IMPLICIT REAL(float)(a-h,o-z)
  IMPLICIT INTEGER (i-n)
  DIMENSION a(*),fit(10)
  n10=n*10
  if(latv==0)then
   n2=n+n
   mbase=(niniz-1)*n2
   do m=niniz,nifin,5
    k=min(m+4,nifin)
    write(iout,40)(j,j,j=m,k)
    do i=1,n
     nbase=mbase+i+i
     l=0
     do l1=m,k
      fit1=a(nbase-1)
      fit2=a(nbase)
      if(abs(fit1)<(1.e-7))fit1=0.
      if(abs(fit2)<(1.e-7))fit2=0.
      fit(l+1)=fit1
      fit(l+2)=fit2
      nbase=nbase+n2
      l=l+2
     enddo
     write(iout,50)i,(fit(l1),l1=1,l)
    enddo
    mbase=mbase+n10
   enddo
  else
   mbase=(niniz-1)*n
   do m=niniz,nifin,10
    k=min(m+9,nifin)
    write(iout,40)(j,j=m,k)
    do i=1,n
     nbase=mbase+i
     l=0
     do l1=m,k
      l=l+1
      fit1=a(nbase)
      if(abs(fit1)<(1.e-7))fit1=0.
      fit(l)=fit1
      nbase=nbase+n
     enddo
     write(iout,50)i,(fit(l1),l1=1,l)
    enddo
    mbase=mbase+n10
   enddo
  endif
50 format(i5,3x,10g12.4)
40 format(/4x,10(8x,i4)/)
 END SUBROUTINE cmacol


 SUBROUTINE kredin
!------------------------------------------------------------!
! This routine was derived from bmat/smat routine.           !
! Reads eigenvectors in IBZ and generate eigenvectors        !
! in the full brillouin zone, applying symmetry operators    !
! and time reversal symmetry.                                !
!------------------------------------------------------------!
  USE numbers
  USE parinf_module
  USE retic_module
  USE epesi_module
  USE xyvdim_module
  USE text_module
  IMPLICIT REAL(float) (a-h,o-z)
  IMPLICIT INTEGER (i-n)
  LOGICAL exist
  REAL(float),DIMENSION(:),ALLOCATABLE :: ar
  REAL(float),DIMENSION(3):: akxyz
  INTEGER,DIMENSION(:),ALLOCATABLE :: norder
  INTEGER,DIMENSION(3):: kr
!-----------------------------!
! Suppress printing (MDT)     !
  LOGICAL :: PRINT=.false.
!-----------------------------!
  CHARACTER(len=6) :: znamz='kredin'

  mvf=inf(2)
  ndf=inf(7)
  nspsta=inf(64)+1
  ndfsqr2=ndf*ndf*2
  allocate(ar(ndfsqr2),stat=ierror )
  if(ierror/=0)call errnic(0,ndfsqr2,znamz,'AR memory allocation')
  nband=ndf*nspsta
  io30=iunit(30)
  inquire(file='KRED.DAT',exist=exist)

  if(.not.exist)then
   call errnic(0,ndf,znamz,'File KRED.DAT not found - eigenvectors')
  else
   open(unit=io30,file='KRED.DAT',form='formatted',status='old')
  endif

  read(io30,101)is1,is2,is3,nkf,bret
  if(nkf>lim001)call errnic(0,nkf,znamz,'Too many k points - increase &
   &lim001 to')
  is1is2=is1*is2
  allocate(norder(is1*is2*is3),stat=ierror )
  if(ierror/=0)call errnic(0,is1*is2*is3,znamz,'Norder memory allocation')
! nkf          number of k points in ibz
! inf(1)       number of symmetry operators extended to inversion
! inf(2)       number of symmetry operators
! inf(7)       number of basis FUNCTIONs (ao)
! inf(10)      translational symmetry: 0 (0d) 1 (1d) 2 (2d) 3 (3d)
! is1(is2,is3) monkhorst net shrinking factors
! inf(64)      0 (restricted closed shell) / 1 (unrestricted open shell)
  ntot=(inf(64)+1)*inf(7)*nkf
! Coordinates of k points in lattice vectors units
  call ireadf(io30,jj,3*nkf)
! Type of k point: COMPLEX (0), REAL (1)
  call ireadf(io30,latvrs,nkf)
! Symmetry operators in lattice vectors units
  call ireadf(io30,irr,9*48)
! Geometrical weight of k points
  call rreadf(io30,wpj,nkf)
  if(ntot>lim018)call errnic(0,ntot,znamz,&
   &'No. of k points x bands too large - increase lim018 to')
! Eigenvalues - unrestricted: alpha(all k points),beta(all k points)
  call rreadf(io30,ene,ntot)
! Alfa - weight of the eigenvalues - computed from fermi energy calculation
  call rreadf(io30,alfa,ntot)

  if(PRINT)write(iout,100)
  if(inf(10)/=0)then
   if(PRINT)then
    write(iout,1121)is1,is2,is3,nkf,is1,(nk,tipo(latvrs(nk)),(jj(i,nk),i=1,3),&
     &nk=1,nkf)
1121 format(' Shrinking factor (monkhorst net)',t35,3i3,&
      &t46,'Number of k points in the ibz',t76,i4/&
      &1x,79('*')/' *** k points coordinates (oblique coordinates in unit&
      &s of is =',i3,')'/(4(i4,'-',a1,'(',3i3,')')))
    write(iout,1103)(wpj(nk),nk=1,nkf)
1103 format(/' Geometrical weight of k points - monkhorst'//(08f10.5))
1104 format(/' Weight of k points for each band'/)
1105 format(' Band n',i4/(15f7.4))
   endif ! PRINT
   nband=nspsta*inf(7)
   ndinz=nkf*nband
   if(PRINT)then
    write(iout,1104)
    do n=1,nband
     write(iout,1105)n,(alfa(nkn),nkn=n,ndinz,nband)
    enddo
    write(iout,1129)
    do i=1,3
     write(iout,1130)(bret(i,j),j=1,3)
    enddo
   endif ! PRINT
  endif ! inf(10)
  jniz=0
  nband=ndf*nspsta
  if(PRINT)then
   do isigma=1,nspsta
    if(nspsta/=1)write(iout,204)mspin(isigma)
    iniz=jniz
    do k=1,nkf
     latw=latvrs(k)
     write(iout,105)k,jj(1,k),jj(2,k),jj(3,k),(ene(iniz+i),i=1,ndf)
     iniz=iniz+nband
    enddo
    jniz=ndf
   enddo
  endif ! PRINT
105 format(/' Eigenvalues (au) - k=',i4,' (',3i3,')',1p/(10e12.4))
1129 format(/' Reciprocal lattice vectors components. (a.u.) '//&
      &6x,1(2x,'x',12x,'y',12x,'z',12x))
1130 format(3f13.7,2x,3f13.7)
100 format(39(' k')//' Eigenvalues and related information are read',&
     &' from file KRED.DAT (full bz)'//39(' k')/)
101 format(4i4/1p,(3e21.13))

  is10=is1*16
  is20=is2*16
  is30=is3*16
  vrsiss=1._float/is1
  iniz=0
  if(PRINT)then
   if(inf(10)==0)then
    write(iout,205)
   else
    write(iout,104)
    write(iout,200)
    do mv1=1,mvf,4
     mv9=min(mv1+3,mvf)
     write(iout,201)(i,i=mv1,mv9)
      do i=1,3
      write(iout,202)((irr(i,j,iv),j=1,3),iv=mv1,mv9)
     enddo
    enddo
   endif
  endif ! PRINT
200 format(/' symmops - reciprocal lattice'/)
201 format(/4(7x,' no.',i3,6x))
202 format(4(3x,3i4,5x))
204 format(/2x,a24)

  jniz=0
  do isigma=1,nspsta
  norder(1:is1is2*is3)=0
  if(nspsta>1.and.PRINT)write(iout,204)mspin(isigma)
  iniz=jniz

  do 1010 k=1,nkf

   latw=latvrs(k)
   if(PRINT)write(iout,115)k,tipo(latw),jj(1,k),jj(2,k),jj(3,k),&
    &(ene(iniz+i),i=1,ndf)
104 format(//' hamiltonian eigenvectors in the full brillouin zone'/)
205 format(//' hamiltonian eigenvectors'/)
115 format(/' eigenvalues - k=',i4,'-',a1,'(',3i3,') (ibz)',1p/(10e12.4))
106 format(/' eigenvectors - k=',i4,'-',a1,'(',3i3,') - symmop ',i4,&
     &' (',3i3,') (coordinates units is=',i3,')'/&
     &'(',3f12.8,' a.u.)')
107 format(/' time reversal symmetry applied ')
   iniz=iniz+nband
   jr1=jj(1,k)
   jr2=jj(2,k)
   jr3=jj(3,k)
   if(latvrs(k)/=0)goto 889

!........................ COMPLEX k-point ............................

   do 776 mv=1,mvf
    mr1=mod(irr(1,1,mv)*jr1+irr(1,2,mv)*jr2+irr(1,3,mv)*jr3+is10,is1)
    mr2=mod(irr(2,1,mv)*jr1+irr(2,2,mv)*jr2+irr(2,3,mv)*jr3+is20,is2)
    mr3=mod(irr(3,1,mv)*jr1+irr(3,2,mv)*jr2+irr(3,3,mv)*jr3+is30,is3)
    kr(1)=mr1
    kr(2)=mr2
    kr(3)=mr3
    nrec=mr1+1+mr2*is1+mr3*is1is2
    if(norder(nrec)/=0)goto 776
    do i=1,3
     akxyz(i)=(bret(1,i)*kr(1)+bret(2,i)*kr(2)+bret(3,i)*kr(3))*vrsiss
    enddo
    if(PRINT)write(iout,106)k,tipo(latw),jr1,jr2,jr3,mv,kr,is1,akxyz
    call ireadf(io30,kr,3)
! Check: k generated vs k written by crystal06
    norder(nrec)=2
    call rreadf(io30,ar,ndf*ndf*(2-latw))
    if(PRINT)call cmacol(ar,ndf,1,ndf,latw)
776 continue

! No time reversal symmetry - skip
   if(mvf==inf(1))goto 777
! Time reversal symmetry applied
   jr1=mod(is1-jr1,is1)
   jr2=mod(is2-jr2,is2)
   jr3=mod(is3-jr3,is3)
   do 778 mv=1,mvf
    mr1=mod(irr(1,1,mv)*jr1+irr(1,2,mv)*jr2+irr(1,3,mv)*jr3+is10,is1)
    mr2=mod(irr(2,1,mv)*jr1+irr(2,2,mv)*jr2+irr(2,3,mv)*jr3+is20,is2)
    mr3=mod(irr(3,1,mv)*jr1+irr(3,2,mv)*jr2+irr(3,3,mv)*jr3+is30,is3)
    kr(1)=mr1
    kr(2)=mr2
    kr(3)=mr3
    nrec=mr1+1+mr2*is1+mr3*is1is2
    if(norder(nrec)/=0)goto 778
    do i=1,3
     akxyz(i)=(bret(1,i)*kr(1)+bret(2,i)*kr(2)+bret(3,i)*kr(3))*vrsiss
    enddo
    if(PRINT)then
     write(iout,107)
     write(iout,106)k,tipo(latw),jj(1,k),jj(2,k),jj(3,k),mv,kr,is1,akxyz
    endif ! PRINT
    call ireadf(io30,kr,3)
    norder(nrec)=2
    call rreadf(io30,ar,ndf*ndf*(2-latw))
    if(PRINT)call cmacol(ar,ndf,1,ndf,latw)
778 continue
777 goto 1010

!.......................... REAL k-point .............................
889 continue
   do 886 mv=1,mvf
    mr1=mod(irr(1,1,mv)*jr1+irr(1,2,mv)*jr2+irr(1,3,mv)*jr3+is10,is1)
    mr2=mod(irr(2,1,mv)*jr1+irr(2,2,mv)*jr2+irr(2,3,mv)*jr3+is20,is2)
    mr3=mod(irr(3,1,mv)*jr1+irr(3,2,mv)*jr2+irr(3,3,mv)*jr3+is30,is3)
    kr(1)=mr1
    kr(2)=mr2
    kr(3)=mr3
    nrec=mr1+1+mr2*is1+mr3*is1is2
    if(norder(nrec)/=0)goto 886
    do i=1,3
     akxyz(i)=(bret(1,i)*kr(1)+bret(2,i)*kr(2)+bret(3,i)*kr(3))*vrsiss
    enddo
    if(PRINT)write(iout,106)k,tipo(latw),jr1,jr2,jr3,mv,kr,is1,akxyz
    call ireadf(io30,kr,3)
!    print 1943,'kr',kr
1943 format(a,9i3)
    norder(nrec)=1
    call rreadf(io30,ar,ndf*ndf*(2-latw))
    if(PRINT)call cmacol(ar,ndf,1,ndf,latw)
886 continue

1010 continue ! end of loop over k points
   iniz=ndf
   enddo
   close(io30)
  END SUBROUTINE kredin


 SUBROUTINE ireadf(jout,ia,n)
  USE numbers
  IMPLICIT REAL(float) (a-h,o-z)
  IMPLICIT INTEGER (i-n)
  DIMENSION ia(n)
  read(jout,*)ia
 END SUBROUTINE ireadf

 SUBROUTINE rreadf(jin,a,n)
  USE numbers
  IMPLICIT REAL(float) (a-h,o-z)
  IMPLICIT INTEGER (i-n)
  DIMENSION a(n)
  read(jin,*)a
 END SUBROUTINE rreadf


 SUBROUTINE errnic(ierr,inic,namz,mzss)
  USE numbers
  USE parinf_module
  IMPLICIT REAL(float) (a-h,o-z)
  IMPLICIT INTEGER (i-n)
  CHARACTER(*) namz
  CHARACTER(*) mzss
  if(ierr==0)then
   write(iout,1)trim(namz),trim(mzss),inic
1  format(' Error **** ',a,' **** ',a,i10)
   stop
  elseif(ierr==1) then
   write(iout,2)trim(namz),trim(mzss),inic
2  format(' Warning **** ',a,' **** ',a,i10)
  else
   write(iout,3)trim(namz),trim(mzss),inic
3  format(' Information **** ',a,' **** ',a,i10)
  endif
 END SUBROUTINE errnic


 SUBROUTINE matint_d(l,nr,ndim)
  USE numbers
  USE parinf_module
  IMPLICIT REAL(float) (a-h,o-z)
  IMPLICIT INTEGER (i-n)
  DIMENSION l(ndim,ndim)
  do m=1,nr,20
   k=min(m+19,nr)
   write(iout,40)(j,j=m,k)
   do i=1,nr
    write(iout,50)i,(l(i,j),j=m,k)
   enddo
  enddo
40 format(/7x,20i6)
50 format(i4,3x,20(i6))
 END SUBROUTINE matint_d

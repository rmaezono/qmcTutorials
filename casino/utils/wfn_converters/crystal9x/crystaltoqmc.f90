!------------------------------------------------------------------------!
! CRYSTALTOQMC format converter                                          !
!                                                                        !
! Read unformatted fort.12, fort.10 and fort.30 files from CRYSTAL95/98  !
! carrying out reduction of k point net ('plucking') if necessary.       !
! Merge these three files and write formatted gwfn.data containing       !
! appropriate data for QMC calculation with CASINO.                      !
!                                                                        !
! NB: CRYSTAL2003 is capable of writing gwfn.data files directly.        !
!                                                                        !
! MDT 1997                                                               !
! MDT Jun 1999 - modified to write new format of gwfn.data file.         !
! MDT Nov 2000 - fully fortran90-ized.                                   !
! MDT Dec 2000 - bug fix to make real/complex k reordering correct in    !
!                big supercells                                          !
! MDT Feb 2001 - added facility to read eigenvalues from fort.30 file    !
!                and write them to gwfn.data so we can do metals.        !
! MDT Mar 2001 - modified to work with CRYSTAL98 as well as CRYSTAL95    !
! MDT Oct 2003 - stopped it asking about k points and supercells if the  !
!                system is finite.                                       !
!------------------------------------------------------------------------!

MODULE utils
 IMPLICIT NONE
 INTEGER,PARAMETER :: dp=kind(1.d0)


 CONTAINS


 SUBROUTINE rread(jin,a,n)
 IMPLICIT NONE
 INTEGER,INTENT(in) :: jin,n
 REAL(dp),INTENT(out) :: a(n)
 read(jin)a
 END SUBROUTINE rread


 SUBROUTINE iread(jin,ia,n)
 IMPLICIT NONE
 INTEGER,INTENT(in) :: jin,n
 INTEGER,INTENT(out) :: ia(n)
 read(jin)ia
 END SUBROUTINE iread


 SUBROUTINE err1(ierr,namz,mzss)
 IMPLICIT NONE
 INTEGER,INTENT(in) :: ierr
 CHARACTER(*),INTENT(in) :: namz,mzss
 if(ierr==0)then
  write(6,1)namz,mzss
1 format(' ERROR **** ',a,' **** ',a)
  stop
 else
  write(6,2)namz,mzss
2 format(' WARNING **** ',a,' **** ',a)
 endif
 END SUBROUTINE err1


 SUBROUTINE minv3(p,pinv,det)
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: p(3,3)
 REAL(dp),INTENT(out) :: pinv(3,3),det
 REAL(dp) :: f1,f2,f3,deti
 REAL(dp),PARAMETER :: toll=1d-16
 f1=p(2,2)*p(3,3)-p(2,3)*p(3,2)
 f2=p(2,3)*p(3,1)-p(2,1)*p(3,3)
 f3=p(2,1)*p(3,2)-p(2,2)*p(3,1)
 det=p(1,1)*f1+p(1,2)*f2+p(1,3)*f3
 if(abs(det)<toll)then
  call err1(1,'minv3 ','Determinant of direct lattice vectors is very small.')
 else
  deti=1.d0/det
  pinv(1,1)=f1*deti
  pinv(2,1)=f2*deti
  pinv(3,1)=f3*deti
  f1=p(1,1)*deti
  f2=p(1,2)*deti
  f3=p(1,3)*deti
  pinv(1,2)=f3*p(3,2)-f2*p(3,3)
  pinv(2,2)=f1*p(3,3)-f3*p(3,1)
  pinv(3,2)=f2*p(3,1)-f1*p(3,2)
  pinv(1,3)=f2*p(2,3)-f3*p(2,2)
  pinv(2,3)=f3*p(2,1)-f1*p(2,3)
  pinv(3,3)=f1*p(2,2)-f2*p(2,1)
 endif
 return
 END SUBROUTINE minv3


 SUBROUTINE mxmb(a,nca,nra,b,ncb,nrb,r,ncr,nrr,ncol,nlink,nrow)
 IMPLICIT NONE
 INTEGER,INTENT(in) :: nca,nra,ncb,nrb,ncr,nrr,ncol,nlink,nrow
 REAL(dp),INTENT(in) :: a(*),b(*)
 REAL(dp),INTENT(inout) :: r(*)
 INTEGER i,j,k,ia,ib,ir,iaa,ibb,irr
 REAL(dp) fac
 ir=1 ; ib=1
 do j=1,nrow
  ibb=ib ; ia=1
  do k=1,nlink
   fac=b(ibb)
   if(fac/=0.d0)then
    irr=ir ; iaa=ia
    do i=1,ncol
     r(irr)=fac*a(iaa)+r(irr)
     irr=irr+ncr ; iaa=iaa+nca
    enddo
   endif
   ibb=ibb+ncb
   ia=ia+nra
  enddo
  ir=ir+nrr ; ib=ib+nrb
 enddo
 END SUBROUTINE mxmb


 CHARACTER(20) FUNCTION i2s(n)
!-----------------------------------------------------------------------!
! I2S                                                                   !
! ===                                                                   !
! Convert integers to left justified strings that can be printed in the !
! middle of a sentence without introducing large amounts of white space.!
!                                                                       !
! Calling routine is intended to include something like:                !
! USE utilities                                                         !
! INTEGER i                                                             !
! i=12                                                                  !
! write(6,*)'Integer number ',trim(i2s(i)),' with words at the end.'    !
!-----------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: n
  INTEGER i,j
  CHARACTER tmp,sign

  if(n==0)then
   i2s='0' ; return
  endif
  sign=' ' ; if(n<0)sign='-'

  do i=1,len(i2s)
   i2s(i:i)=' '
  enddo

  i=abs(n)
  do j=1,len(i2s)
   if(i==0)exit
   i2s(j:j)=achar(ichar('0')+mod(i,10))
   i=i/10
  enddo

  i=1 ; j=len_trim(i2s)
  do
   if(i>=j)exit
   tmp=i2s(j:j)
   i2s(j:j)=i2s(i:i)
   i2s(i:i)=tmp
   i=i+1
   j=j-1
  enddo

  i2s=trim(sign)//i2s

 END FUNCTION i2s


END MODULE utils


PROGRAM crystaltoqmc
!---------------------------!
! Main program starts here. !
!---------------------------!
 USE utils
 IMPLICIT NONE

! Emulate CRYSTAL9X f77 common blocks:
! parinf common
 INTEGER,PARAMETER :: limpar=50,limtol=50,liminf=200
 REAL(dp) :: par(limpar)
 INTEGER itol(limtol),inf(liminf),iout
 CHARACTER(4) :: itit(20)
 CHARACTER,ALLOCATABLE :: newtit(:)
! xyvdim common
 REAL(dp) xyv(9,48),trasv(3,48)
 INTEGER ninv(48)
! gvect common
 REAL(dp) paret(3,3),w1r(3,3)
 REAL(dp),ALLOCATABLE :: gmodus(:),xg(:,:),xa(:,:)
 INTEGER,ALLOCATABLE :: nm(:),mn(:),nn1(:),lg(:,:)
! basato common
 REAL(dp),ALLOCATABLE :: aznuc(:),che(:),exad(:),xl(:,:),exx(:),     &
  &c1(:),c2(:),c3(:),cmax(:),c2w(:),c3w(:)
 INTEGER,ALLOCATABLE :: nat(:),nshpri(:),ipseud(:),laa(:),lan(:),lat(:),     &
  &latao(:),ndq(:),latoat(:)
! psiijj common
 INTEGER ii,jj,itypse
! infpot common
 REAL(dp),ALLOCATABLE :: apot(:),cpot(:)
 INTEGER,ALLOCATABLE :: npot(:),nbtyp(:),nsom(:)
! spinor common
 INTEGER,ALLOCATABLE :: ispin(:),ichmod(:)
! molbar common
 REAL(dp),ALLOCATABLE :: xbar(:,:)
 INTEGER,ALLOCATABLE :: n1mol(:)
! basold common
 REAL(dp) parold(3,3),trasvo(3,48)
 REAL(dp),ALLOCATABLE :: exaold(:),xold(:,:),xlold(:,:),hmodus(:),xgold(:,:)
! knetout output
 INTEGER,ALLOCATABLE :: latvrs2(:),jjj(:,:)
 REAL(dp),ALLOCATABLE :: wpj(:),ene(:)

! Local stuff
 REAL(dp),PARAMETER :: toll=1.d-6
 INTEGER io10,io12,idum,laf,naf,nprim,laf3,laf1,naf3,naf1,inf5,inf79,inf793, &
  &ndfrf,ndfcf,k,ncell(3),is1,is2,is3,mvf,irr(3,3,48),i,j,iv,jv,l1,&
  &l2,l3,iswp,lll,laz(48),lbz(48),lcz(48),lay(48),lby(48),lcy(48),mv,is1m1,  &
  &is2m1,is10,is20,is21,is30,k1,k2,k3,isigma,nc(3),iq1,iq2,iq3,nkfqmc,nkf2,  &
  &ichoose,nkf,norb,n,jbase,kbase,nb,titlength,titsize,ind,&
  &nelec,na,kr,kc,ik,ck_size,j3(1),ialloc,kmax,nk,ntot,io30
 INTEGER,ALLOCATABLE :: mjj(:,:),latvrs(:),latqmc(:),mqmc(:,:),kpoint(:),    &
  &jk(:,:),ak_pointer(:,:),kir(:),kic(:),kir_plucked(:),kic_plucked(:)
 REAL(dp) brot(3,3),bret(3,3),crot(3,3),drot(3,3),det,ris1,ris2,     &
  &ris3,rtmp(3),rem1,rem2,rem3,root3,harmonic_coeff_times_norm(28,8),  &
  &fraction_dftx
 REAL(dp),ALLOCATABLE :: akt(:),kvec(:),ktemp(:,:),kvec_reordered(:),&
  &c_prim(:),ak(:),ckba(:)
 LOGICAL input,eigenvalues_present,periodic
 LOGICAL,ALLOCATABLE :: issp(:)
 CHARACTER(13),PARAMETER :: reals3="(3(1pe20.13))"
 CHARACTER(13),PARAMETER :: reals4="(4(1pe20.13))"
 CHARACTER(13),PARAMETER :: int8="(8i10)"
 CHARACTER(5) exch_functional(5),corr_functional(8)
 DATA exch_functional/'LDA','VBH','Becke','PBE','PWGGA'/
 DATA corr_functional/'PWLSD','PZ','VWN','VBH','P86','PWGGA','LYP','PBE'/

!---------------------------------------------------------------------------!

 ialloc=0

! Open the relevant files
 io10=10 ; io12=12 ; iout=13 ; io30=30
 open(io10,status='unknown',form='unformatted',err=100)
 open(io12,status='unknown',form='unformatted',err=12)
 open(iout,file='gwfn.data',status='unknown',form='formatted')

! Check if it's a molecule in which case we don't need to worry about k points
 read(io12)(idum,i=1,32)
 read(io12)par,itol,inf
 if(inf(10)==0)then
  periodic=.false.
 else
  periodic=.true.
 endif

! Ask about supercell sizes etc. if we have PBCs.
 if(periodic)then

  do ! until 1 or 2 is input
   write(6,*)'Choose:'
   write(6,*)'(1) Standard k space net reduction (plucking)'
   write(6,*)'(2) Manual selection of k points'
   write(6,*)
   read(5,*)ichoose
   if(ichoose==1.or.ichoose==2)exit
  enddo

  if(ichoose==1)then
   input=.false.
   write(6,*)'Monkhorst net shrinking factors used in CRYSTAL calculation?'
   write(6,*)'(e.g. 8 8 8)'
   read(5,*)nc(1:3)
   write(6,*)'Supercell size to be used in CASINO QMC calculation?'
   write(6,*)'(e.g. 2 2 2)'
   read(5,*)ncell(1:3)
  else ! ichoose==2
   input=.true.
   write(6,*)'Monkhorst net shrinking factors used in CRYSTAL calculation?'
   read(5,*)ncell(1:3)
   nc=ncell
   write(6,*)'How many k points do you want?'
   read(5,*)nkfqmc
   allocate(kpoint(nkfqmc),stat=ialloc) ; if(ialloc/=0)goto 1
   kpoint=0.d0
   write(6,*)'Enter the k point sequence numbers:'
   write(6,*)'(see CRYSTAL output file)'
   do i=1,nkfqmc
    read(5,*)kpoint(i)
   enddo
  endif

 else ! finite system

  input=.false.
  nc(1:3)=1
  ncell(1:3)=1

 endif ! periodic or not

 kmax=nc(1)*nc(2)*nc(3)

! Check QMC net is a subset of the CRYSTAL one.
 is1=ncell(1) ; is2=ncell(2) ; is3=ncell(3)
 is21=is1*is2
 iswp=is21*is3
 is1m1=is1-1
 is2m1=is2-1
 is10=is1*16
 is20=is2*16
 is30=is3*16
 ris1=1.d0/dble(is1) ; ris2=1.d0/dble(is2) ; ris3=1.d0/dble(is3)
 rem1=mod(dble(nc(1)),dble(is1))
 rem2=mod(dble(nc(2)),dble(is2))
 rem3=mod(dble(nc(3)),dble(is3))
 if(rem1>toll.or.rem2>toll.or.rem3>toll)then
  write(6,*)
  write(6,*)'Tough.'
  write(6,*)'The CRYSTAL MP net does not contain the necessary k points'
  write(6,*)'to do a QMC calculation with a supercell of this size.'
  write(6,*)
  stop
 endif
 iq1=nc(1)/is1 ; iq2=nc(2)/is2 ; iq3=nc(3)/is3
 allocate(latqmc(kmax),issp(kmax),mqmc(3,kmax),mjj(3,kmax),latvrs(kmax), &
  &kvec(3*kmax),kvec_reordered(3*kmax),kir(kmax),kic(kmax),ktemp(3,kmax),&
  &jk(kmax,2),ak_pointer(kmax,2),kir_plucked(kmax),kic_plucked(kmax),    &
  &stat=ialloc) ; if(ialloc/=0)goto 2
 latqmc=0
 if(.not.input)nkfqmc=0

! Start to read in the CRYSTAL files:

! read 32 lum0xx dimensions which we don't care about
! read(io12)(idum,i=1,32)

! read parinf/gvect/xyvdim commons
! read(io12)par,itol,inf

 read(io12)itit
 titlength=len(itit)
 titsize=size(itit)
 titlength=titlength*titsize
 allocate(newtit(titlength))
 backspace 12
 read(io12)newtit
 read(io12)paret,w1r,xyv,trasv,ninv

! Prepare symmetry stuff (parts of the following lifted from CRYSTAL)
 mvf=inf(1)
 det=1.d0/par(34)
 do i=1,3
  do j=1,3
   brot(j,i)=paret(i,j)*det
  enddo
 enddo
 call minv3(brot,bret,det)
 do iv=1,mvf
  do i=1,3
   do j=1,3
    crot(j,i)=0.d0
    drot(j,i)=0.d0
   enddo
  enddo
  call mxmb(bret,1,3,xyv(1,iv),1,3,crot,1,3,3,3,3)
  call mxmb(crot,1,3,brot,1,3,drot,1,3,3,3,3)
  jv=ninv(iv)
  do i=1,3
   do j=1,3
    irr(j,i,jv)=nint(drot(i,j))
   enddo
  enddo
 enddo

! Map the QMC supercell k points onto the CRYSTAL Monkhorst-Pack net
! (BZ symmetry stuff not necessary for now but leave it for possible future use)
 mvf=1 ! number of symmops = 1
 issp=.false.
 nkf=0
 lll=0
 do l3=0,is3-1
  laz(1:mvf)=irr(1,3,1:mvf)*l3+is10
  lbz(1:mvf)=irr(2,3,1:mvf)*l3+is20
  lcz(1:mvf)=irr(3,3,1:mvf)*l3+is30
  do l2=0,is2m1
   lay(1:mvf)=irr(1,2,1:mvf)*l2+laz(1:mvf)
   lby(1:mvf)=irr(2,2,1:mvf)*l2+lbz(1:mvf)
   lcy(1:mvf)=irr(3,2,1:mvf)*l2+lcz(1:mvf)
   do l1=0,is1m1
    lll=lll+1
    if(issp(lll))cycle
    nkf=nkf+1
    mqmc(1,nkf)=l1*iq1 ; mqmc(2,nkf)=l2*iq2 ; mqmc(3,nkf)=l3*iq3
    do mv=1,mvf
     k3=mod(irr(3,1,mv)*l1+lcy(mv),is3)
     k2=mod(irr(2,1,mv)*l1+lby(mv),is2)
     k1=mod(irr(1,1,mv)*l1+lay(mv),is1)
     k=k3*is21+k2*is1+k1
     if(issp(k+1))cycle
     issp(k+1)=.true.
     if(k1/=0)k1=is1-k1
     if(k2/=0)k2=is2-k2
     if(k3/=0)k3=is3-k3
     k=k3*is21+k2*is1+k1
     if(issp(k+1))cycle
     issp(k+1)=.true.
    enddo ! mv (symmops)
   enddo ! l1
  enddo ! l2
 enddo ! l3

 if(input)then
  if(any(kpoint(:)>nkf))then
   write(6,*)
   write(6,'(1x,a,i5,a)')'ERROR: Selected k point numbers outside range (1 : ',&
    &nkf,')'
   stop
  endif
 endif

! Carry on reading and writing the CRYSTAL files (now we know what inf(25) is)
 if(.not.input)nkfqmc=nkf
 inf(25)=nkfqmc
! read basato common
 laf=inf(20)   ! number of shells
 naf=inf(24)   ! number of atoms
 nprim=inf(75) ! number of Gaussian primitives in the basis
 laf3=laf*3
 laf1=laf+1
 naf3=naf*3
 naf1=naf+1
 allocate(aznuc(naf),xa(3,naf),che(laf),exad(laf),xl(3,laf),exx(nprim),      &
  &c1(nprim),c2(nprim),c3(nprim),cmax(nprim),c2w(nprim),c3w(nprim),nat(naf), &
  &nshpri(naf1),ipseud(naf),laa(laf1),lan(laf),lat(laf),latao(laf),ndq(laf1),&
  &latoat(laf),c_prim(nprim),stat=ialloc) ; if(ialloc/=0)goto 3
 call rread(io12,aznuc,naf)
 call rread(io12,xa,naf3)
 call rread(io12,che,laf)
 call rread(io12,exad,laf)
 call rread(io12,xl,laf3)
 call rread(io12,exx,nprim)
 call rread(io12,c1,nprim)
 call rread(io12,c2,nprim)
 call rread(io12,c3,nprim)
 call rread(io12,cmax,nprim)
 call rread(io12,c2w,nprim)
 call rread(io12,c3w,nprim)
 call iread(io12,nat,naf)
 call iread(io12,nshpri,naf1)
 call iread(io12,ipseud,naf)
 call iread(io12,laa,laf1)
 call iread(io12,lan,laf)
 call iread(io12,lat,laf)
 call iread(io12,latao,laf)
 call iread(io12,ndq,laf1)
 call iread(io12,latoat,laf)
 deallocate(exad,cmax,c2w,c3w,ipseud,latoat,ndq)

! read psiijj/infpot commons
 if(inf(31)/=0)then
  read(io12)ii,jj,itypse
  allocate(apot(ii),cpot(ii),npot(ii),nbtyp(jj),nsom(itypse),stat=ialloc)
  if(ialloc/=0)goto 4
  call rread(io12,apot,ii)
  call rread(io12,cpot,ii)
  call iread(io12,npot,ii)
  call iread(io12,nbtyp,jj)
  call iread(io12,nsom,itypse)
  deallocate(apot,cpot,npot,nbtyp,nsom)
 endif

! read spinor common
 allocate(ispin(naf),ichmod(naf),stat=ialloc) ; if(ialloc/=0)goto 5
 call iread(io12,ispin,naf)
 call iread(io12,ichmod,naf)
 deallocate(ispin,ichmod)

! read molbar common
 if(inf(34)/=0)then
  allocate(xbar(3,naf),n1mol(naf),stat=ialloc) ; if(ialloc/=0)goto 6
  call rread(io12,xbar,naf3)
  call iread(io12,n1mol,naf)
  deallocate(xbar,n1mol)
 endif

! read basold common
 inf5=inf(5)+1
 if(inf(15)/=0)then
  allocate(exaold(laf),xold(3,naf),xlold(3,laf),hmodus(inf5), &
   &xgold(3,inf(79)*3),stat=ialloc) ; if(ialloc/=0)goto 7
  read(io12)parold,trasvo
  call rread(io12,exaold,laf)
  call rread(io12,xold,naf3)
  call rread(io12,xlold,laf3)
  call rread(io12,hmodus,inf5)
  call rread(io12,xgold,inf(79)*3)
  deallocate(exaold,xold,xlold,hmodus,xgold)
 endif

! read gvect common
 inf79=inf(79)
 inf793=inf79*3
 allocate(gmodus(inf5),xg(3,inf79),nm(inf5),mn(inf5),nn1(inf79),lg(3,inf79), &
  &stat=ialloc) ; if(ialloc/=0)goto 8
 call rread(io12,gmodus,inf5)
 call rread(io12,xg,inf793)
 call iread(io12,nm,inf5)
 call iread(io12,mn,inf5)
 call iread(io12,nn1,inf79)
 call iread(io12,lg,inf793)
 deallocate(gmodus,xg,nm,mn,nn1,lg)

! Play with reciprocal lattice
! ----------------------------

! Generate the Monkhorst-Pack net data so we can figure out which points are
! real and which are complex i.e. whether lie on the BZ edge or not.

 is1=nc(1) ; is2=nc(2) ; is3=nc(3)
 is21=is1*is2
 iswp=is21*is3
 is1m1=is1-1
 is2m1=is2-1
 is10=is1*16
 is20=is2*16
 is30=is3*16
 ris1=1.d0/is1 ; ris2=1.d0/is2 ; ris3=1.d0/is3
 ind=1
 kvec=0.d0
 mvf=1 ! should set to inf(50) if you ever need to account for point symmetry
 issp(1:iswp)=.false.
 nkf=0
 lll=0
 do l3=0,is3-1
  laz(1:mvf)=irr(1,3,1:mvf)*l3
  lbz(1:mvf)=irr(2,3,1:mvf)*l3
  lcz(1:mvf)=irr(3,3,1:mvf)*l3
  do l2=0,is2m1
   lay(1:mvf)=irr(1,2,1:mvf)*l2+laz(1:mvf)
   lby(1:mvf)=irr(2,2,1:mvf)*l2+lbz(1:mvf)
   lcy(1:mvf)=irr(3,2,1:mvf)*l2+lcz(1:mvf)
   do l1=0,is1m1
    lll=lll+1
    if(issp(lll))cycle
    nkf=nkf+1

! Define which k points are part of both the QMC net and the CRYSTAL net
    if(input)then
     do i=1,nkfqmc
      latqmc(kpoint(i))=1
     enddo
    else
     do i=1,nkfqmc
      if(l1==mqmc(1,i).and.l2==mqmc(2,i).and.l3==mqmc(3,i))then
       latqmc(nkf)=1
       exit ! i.e. goto mjj statement
      endif
     enddo
    endif

    mjj(1,nkf)=l1 ; mjj(2,nkf)=l2 ; mjj(3,nkf)=l3
    rtmp(1)=l1*ris1 ; rtmp(2)=l2*ris2 ; rtmp(3)=l3*ris3
    call mxmb(bret,3,1,rtmp,1,1,kvec(ind),1,1,3,3,1)
    ind=ind+3
    do mv=1,mvf
     k3=mod(irr(3,1,mv)*l1+lcy(mv),is3)
     k2=mod(irr(2,1,mv)*l1+lby(mv),is2)
     k1=mod(irr(1,1,mv)*l1+lay(mv),is1)
     k=k3*is21+k2*is1+k1
     if(issp(k+1))cycle
     issp(k+1)=.true.
     if(k1/=0)k1=is1-k1
     if(k2/=0)k2=is2-k2
     if(k3/=0)k3=is3-k3
     k=k3*is21+k2*is1+k1
     if(issp(k+1))cycle
     issp(k+1)=.true.
    enddo ! mv (symmops)
   enddo ! l1
  enddo ! l2
 enddo ! l3

 do k=1,nkf
  latvrs(k)=0
  if(mod(mjj(1,k)*2,is1)/=0)cycle
  if(mod(mjj(2,k)*2,is2)/=0)cycle
  if(mod(mjj(3,k)*2,is3)==0)latvrs(k)=1
 enddo

 deallocate(issp,mqmc,mjj)

 if(periodic)then

  if(nkfqmc==1.and.nkf==1)then
   write(6,'(1x,4a)')trim(i2s(nkfqmc)),' QMC k point generated from ',&
    &trim(i2s(nkf)),' CRYSTAL k point.'
  elseif(nkfqmc==1.and.nkf/=1)then
   write(6,'(1x,4a)')trim(i2s(nkfqmc)),' QMC k point generated from ',&
    &trim(i2s(nkf)),' CRYSTAL k points.'
  else
   write(6,'(1x,4a)')trim(i2s(nkfqmc)),' QMC k points generated from ',&
    &trim(i2s(nkf)),' CRYSTAL k points.'
  endif

  if(nkfqmc==1)then
   ind=0
   do k=1,nkf
    if(latqmc(k)==1)then
     write(6,*)
     write(6,*)'Wrote out single k point information'
     write(6,*)'Coordinates (Angstrom):'
     write(6,'(4(1pe20.13))')(kvec(ind+i),i=1,3)
     write(6,*)'Complex k point flag (1=real, 0=complex):'
     write(6,'(8i10)')latvrs(k)
     write(6,*)
    endif
    ind=ind+3
   enddo
  endif

 endif ! periodic

! Now read the eigenvalues from the fort.30 file (periodic systems only)
! This bit introduced 2.2001 MDT (N.B. now reading in stuff that we
! have just had to laboriously calculate.. Sigh.)

 if(periodic)then
  inquire(file='fort.30',exist=eigenvalues_present)
  if(eigenvalues_present)then
   open(io30,status='unknown',form='unformatted',err=30)
   read(io30)idum,nkf2,idum,bret
   if(nkf2/=nkf)then
    write(6,*)'Calculated NKF  :',nkf
    write(6,*)'NKF from fort.30: ',nkf2
    call err1(0,'CRYSTALTOQMC','Calculated NKF value different&
     & from that given in the fort.30 file.')
   endif
   ntot=(inf(64)+1)*inf(7)*nkf2
   allocate(jjj(3,nkf2),latvrs2(nkf2),wpj(nkf2),ene(ntot),stat=ialloc)
   call iread(io30,jjj,3*nkf2)
   call iread(io30,latvrs2,nkf2)
   call iread(io30,irr,9*48)
   call rread(io30,wpj,nkf2)
   call rread(io30,ene,ntot)
   deallocate(jjj,latvrs2,wpj)
  endif
 endif

! Now write out the gwfn.data file for CASINO to read

!=================Title section=============================

! Hack the length of the title. This is necessary because of the boring
! f77 method of writing out an unformatted variable length title that
! CRYSTAL uses.
 do i=1,titlength
  j3=iachar(newtit(i))
  if(j3(1)<32.or.j3(1)>122)then ! we have a non-ascii character (ugh!)
   titlength=i-1
   exit
  endif
 enddo
 write(iout,fmt='(80a1)')(newtit(i),i=1,titlength)

!=================Basic information=========================

 write(iout,fmt="(/'BASIC INFO'/'----------')")
! replace_me_with_name_of_code is done by the crysgen script
 write(iout,fmt="('Generated by:'/'replace_me_with_name_of_code')")
 write(iout,fmt="('Method:')")
 if(inf(170)/=0)then
  write(iout,fmt="('DFT')")
 else
  if(inf(27)==0.and.inf(64)==0)then
   write(iout,fmt="('RHF')")
  elseif(inf(27)==1.and.inf(64)==1)then
   write(iout,fmt="('UHF')")
  elseif(inf(27)==0.and.inf(64)==1)then
   write(iout,fmt="('ROHF')")
  endif
 endif
 write(iout,fmt="('DFT functional:')")

 if(inf(170)/=0)then ! DFT
  if(inf(172)/=0)then
   fraction_dftx=1.d0-par(48)*0.01d0
   if(fraction_dftx<1.d0)then ! hybrid xc functional
    write(iout,fmt="('Fock exchange * ',f9.5,' + ',a,' exchange * ',f9.5, &
     &' + ',a,' correlation')")1.d0-fraction_dftx,                        &
     &trim(adjustl(exch_functional(inf(172)))),fraction_dftx,             &
     &trim(adjustl(corr_functional(inf(171))))
   else ! non-hybrid xc functional
    write(iout,fmt="(a,' correlation, ',a,' exchange')")                  &
     &trim(adjustl(corr_functional(inf(171)))),                           &
     &trim(adjustl(exch_functional(inf(172))))
   endif
  else ! correlation-only functional
   write(iout,fmt="('Hartree-Fock exchange plus ',a,' correlation.')")    &
    &trim(adjustl(corr_functional(inf(171))))
  endif
 else ! Hartree-Fock
  write(iout,fmt="('none')")
 endif

 write(iout,fmt="('Periodicity:'/i1)")inf(10)
 write(iout,fmt="('Spin unrestricted:')")
 if(inf(64)/=0)then
  write(iout,fmt="('.true.')")
 else
  write(iout,fmt="('.false.')")
 endif
 write(iout,fmt="('Nuclear repulsion energy (au/atom):')")
 write(iout,*)par(18)/dble(naf)

 nelec=0
 do na=1,naf
  do l1=nshpri(na),nshpri(na+1)-1
   nelec=nelec+int(che(l1))
  enddo
 enddo

 deallocate(che)

 write(iout,fmt="('Number of electrons per primitive cell:')")
 write(iout,*)nelec

!=================Geometry input==============================

 write(iout,fmt="(/'GEOMETRY'/'--------')")

 write(iout,fmt="('Number of atoms:')")
 write(iout,*)naf

 write(iout,fmt="('Atomic positions (au):')")
 write(iout,fmt=reals3)((xa(i,j),i=1,3),j=1,naf)

 write(iout,fmt="('Atomic numbers for each atom:')")
 write(iout,fmt=int8)(nat(i),i=1,naf)

 write(iout,fmt="('Valence charges for each atom:')")
 write(iout,fmt=reals4)(aznuc(i),i=1,naf)

 deallocate(aznuc,xa,nat)

 if(inf(10)==0)then ! molecules

  kr=1 ; kc=0 ; kir_plucked(1)=1

 else ! solids,polymers,slabs

  write(iout,fmt="('Primitive lattice vectors (au):')")
  write(iout,fmt=reals3)((paret(i,j),j=1,3),i=1,3) ! i.e. don't write(")paret

!=================K SPACE NET=================================

  write(iout,fmt="(/'K SPACE NET'/'-----------')")

  write(iout,fmt="('Number of k points')")
  write(iout,*)nkfqmc

  ik=1 ; nk=0 ; kr=0 ; kc=0
  do k=1,nkf
   if(latqmc(k)==1)then
    nk=nk+1
    ktemp(1,k)=kvec(ik)
    ktemp(2,k)=kvec(ik+1)
    ktemp(3,k)=kvec(ik+2)
    if(latvrs(k)==1)then
     kr=kr+1
     kir(kr)=k
     kir_plucked(kr)=nk
    else
     kc=kc+1
     kic(kc)=k
     kic_plucked(kc)=nk
    endif
   endif
   ik=ik+3
  enddo
  ik=1
  do i=1,kr
   k=kir(i)
   kvec_reordered(ik)=ktemp(1,k)
   kvec_reordered(ik+1)=ktemp(2,k)
   kvec_reordered(ik+2)=ktemp(3,k)
   ik=ik+3
  enddo
  do i=1,kc
   k=kic(i)
   kvec_reordered(ik)=ktemp(1,k)
   kvec_reordered(ik+1)=ktemp(2,k)
   kvec_reordered(ik+2)=ktemp(3,k)
   ik=ik+3
  enddo

  write(iout,fmt="('Number of ''real'' k points on BZ edge')")
  write(iout,*)kr

  write(iout,fmt="('k point coordinates (au)')")
  write(iout,fmt=reals3)(kvec_reordered(k),k=1,nkfqmc*3)

  deallocate(ktemp)
 endif ! polymers/slabs/solids

!=================Basis set===================================

 write(iout,fmt="(/'BASIS SET'/'---------')")

 write(iout,fmt="('Number of Gaussian centres')")
 write(iout,*)naf

 write(iout,fmt="('Number of shells per primitive cell')")
 write(iout,*)laf

 write(iout,fmt="('Number of basis functions (''AO'') per primitive cell')")
 write(iout,*)inf(7)

 write(iout,fmt="('Number of Gaussian primitives per primitive cell')")
 write(iout,*)nprim

 j=maxval(lat(1:laf))

 write(iout,fmt="('Highest shell angular momentum (s/p/d/f/g... 1/2/3/4/5...)'&
  &)")
 if(j==0)then
  write(iout,*)j+1
 elseif(j==1.or.j==2)then
  write(iout,*)2
 elseif(j==3)then
  write(iout,*)j
 else
  call err1(0,'CRYSTALTOQMC','Incorrect value of LAT.')
 endif

 write(iout,fmt="('Code for shell types (s/sp/p/d/f... 1/2/3/4/5...)')")
 write(iout,fmt=int8)(lat(i)+1,i=1,laf)

 write(iout,fmt="('Number of primitive Gaussians in each shell')")
 write(iout,fmt=int8)(lan(i),i=1,laf)

 write(iout,fmt="('Sequence number of first shell on each centre')")
 write(iout,fmt=int8)(nshpri(i),i=1,naf+1)

 write(iout,fmt="('Exponents of Gaussian primitives')")
 write(iout,fmt=reals4)(exx(i),i=1,nprim)

 write(iout,fmt="('Correctly normalised contraction coefficients')")
 do i=1,laf
  if(lat(i)==0.or.lat(i)==1)then
   c_prim(laa(i):laa(i+1)-1)=c1(laa(i):laa(i+1)-1)
  elseif(lat(i)==2)then
   c_prim(laa(i):laa(i+1)-1)=c2(laa(i):laa(i+1)-1)
  elseif(lat(i)==3)then
   c_prim(laa(i):laa(i+1)-1)=c3(laa(i):laa(i+1)-1)
  else
   write(iout,fmt=*)'Duh.'
  endif
 enddo
 write(iout,fmt=reals4)(c_prim(i),i=1,nprim)

 if(any(lat(1:laf)==1))then
  write(iout,fmt="('2nd contraction coefficients (p coeff.,&
   & for sp shells, 0.0 otherwise)')")
  write(iout,fmt=reals4)(c2(i),i=1,nprim)
 endif

 write(iout,fmt="('Position of each shell (au)')")
 write(iout,fmt=reals3)((xl(i,j),i=1,3),j=1,laf)

 deallocate(xl,exx,c1,c2,c3,laa,lan,c_prim,nshpri)

!=================Multideterminant information============================

 write(iout,fmt="(/'MULTIDETERMINANT INFORMATION'/&
  &'----------------------------')")
 write(iout,fmt="('GS')") ! no CRYSTAL excited state calcs, sadly.

!=================Eigenvector coefficients================================

 write(iout,fmt="(/'ORBITAL COEFFICIENTS'/'-----------------------')")

! Read in eigenvector coefficients from fort.10, discarding those from
! k points we don't want and store in temporary vector AK

 ck_size=(kr*inf(7)*inf(7)+(nkfqmc-kr)*2*inf(7)*inf(7))*(inf(64)+1)

 allocate(ak(ck_size),ckba(ck_size),akt(inf(7)*inf(7)*2),stat=ialloc)
 if(ialloc/=0)goto 10

 ind=0
 ndfrf=inf(7)*inf(7)
 ndfcf=ndfrf+ndfrf
 do isigma=0,inf(64)
  do k=1,nkf
   if(latvrs(k)==1)then
    call rread(io10,akt,ndfrf)
    if(latqmc(k)==1)then
     do i=1,ndfrf
      if(abs(akt(i))>1.d-12)then
       ak(ind+i)=akt(i)
      else
       ak(ind+i)=0.d0
      endif
     enddo
     ind=ind+ndfrf
    endif
   else
    call rread(io10,akt,ndfcf)
    if(latqmc(k)==1)then
     do i=1,ndfcf
      if(abs(akt(i))>1.d-12)then
       ak(ind+i)=akt(i)
      else
       ak(ind+i)=0.d0
      endif
     enddo
     ind=ind+ndfcf
    endif
   endif
  enddo
 enddo

 deallocate(akt)

! Reorder eigenvector coefficients in ak --> ckba
! Also multiply coefficients of d basis functions by the
! harmonic coefficients * a missing normalization piece
! ao        harmonic coefficient      normalization
! 3z2-r2         0.5                    1
! xz             3.0                    1/root3
! yz             3.0                    1/root3
! x2-y2          3.0                    1/(2*root3)
! xy             6.0                    1/(2*root3)

 n=0
 do isigma=1,inf(64)+1
  nk=0
  do k=1,nkf
   if(latqmc(k)==1)then
    nk=nk+1
    ak_pointer(nk,isigma)=n
    if(latvrs(k)==1)then
     n=n+ndfrf
    else
     n=n+ndfcf
    endif
   endif
  enddo
 enddo
 do isigma=1,inf(64)+1
  do k=1,kr
   jk(k,isigma)=ak_pointer(kir_plucked(k),isigma)
  enddo
  do k=1,kc
   jk(kr+k,isigma)=ak_pointer(kic_plucked(k),isigma)
  enddo
 enddo

 deallocate(ak_pointer,kir_plucked,kic_plucked)

! s
 harmonic_coeff_times_norm(1,1)=1.d0
!sp
 harmonic_coeff_times_norm(1:4,2)=1.d0
! p
 harmonic_coeff_times_norm(1:3,3)=1.d0
!d (harmonic)
 root3=sqrt(3.d0)
 harmonic_coeff_times_norm(1,4)=0.5d0
 harmonic_coeff_times_norm(2,4)=root3
 harmonic_coeff_times_norm(3,4)=root3
 harmonic_coeff_times_norm(4,4)=0.5d0*root3
 harmonic_coeff_times_norm(5,4)=root3

 j=1
 naf=inf(7) ! total number of bands (occ/unocc) and number of basis functions
 do isigma=1,inf(64)+1
  do k=1,kr
   i=1
   jbase=jk(k,isigma)
   do nb=1,naf
    do l1=1,laf
     norb=latao(l1)
     kbase=lat(l1)+1
     do n=1,norb
      ckba(j)=ak(jbase+i)*harmonic_coeff_times_norm(n,kbase)
      i=i+1
      j=j+1
     enddo ! ao in shell
    enddo ! shells
   enddo ! bands
  enddo ! real k
  do k=kr+1,nkfqmc
   i=1
   jbase=jk(k,isigma)
   do nb=1,naf
    do l1=1,laf
     norb=latao(l1)
     kbase=lat(l1)+1
     do n=1,norb
      ckba(j)=ak(jbase+i)*harmonic_coeff_times_norm(n,kbase)
      ckba(j+1)=ak(jbase+i+1)*harmonic_coeff_times_norm(n,kbase)
      i=i+2
      j=j+2
     enddo ! ao in shell
    enddo ! shells
   enddo ! bands
  enddo ! real k
 enddo ! spins

 deallocate(ak,jk,latao,lat)

 write(iout,fmt=reals4)(ckba(i),i=1,ck_size)

 deallocate(ckba)

!====================Write out band energies for periodic systems=========

 if(inf(10)/=0.and.eigenvalues_present)then
  write(iout,fmt="(/'EIGENVALUES'/'-----------')")
! First write out the band energies at all the 'real' k points
  ik=0 ; k=0
  n=inf(7)*(inf(64)+1)
  do j=1,nkf
   if(latqmc(j)==1.and.latvrs(j)==1)then
    k=k+1
    jbase=(j-1)*n
    do isigma=1,inf(64)+1
     if(inf(64)==0)then
      write(iout,'(a,i6,3f14.8)')'k ',k,kvec_reordered(ik+1),&
      &kvec_reordered(ik+2),kvec_reordered(ik+3)
     else
      write(iout,'(a,i2,a,i6,3f14.8)')'spin ',isigma,' k ',k,&
      &kvec_reordered(ik+1),kvec_reordered(ik+2),kvec_reordered(ik+3)
     endif
     if(isigma==1)then
      write(iout,fmt=reals4)(ene(i),i=jbase+1,jbase+inf(7))
     else
      write(iout,fmt=reals4)(ene(i),i=jbase+inf(7)+1,jbase+inf(7)+inf(7))
     endif
    enddo
    ik=ik+3
   endif
  enddo ! k
! Now write out the band energies for the 'complex' k points
  do j=1,nkf
   if(latqmc(j)==1.and.latvrs(j)==0)then ! 'complex' k points
    k=k+1
    jbase=(j-1)*n
    do isigma=1,inf(64)+1
     if(inf(64)==0)then
      write(iout,'(a,i6,3f14.8)')'k ',k,kvec_reordered(ik+1),&
      &kvec_reordered(ik+2),kvec_reordered(ik+3)
     else
      write(iout,'(a,i2,a,i6,3f14.8)')'spin ',isigma,' k ',k,&
      &kvec_reordered(ik+1),kvec_reordered(ik+2),kvec_reordered(ik+3)
     endif
     if(isigma==1)then
      write(iout,fmt=reals4)(ene(i),i=jbase+1,jbase+inf(7))
     else
      write(iout,fmt=reals4)(ene(i),i=jbase+inf(7)+1,jbase+inf(7)+inf(7))
     endif
    enddo
    ik=ik+3
   endif
  enddo ! k
 endif

!====================Append any other info. to the input file===========

 write(iout,fmt=&
  & "(/'Input and output files for this calculation (not read by CASINO)')")
 write(iout,fmt=&
  & "('=====================================================',&
  & '===========================')")

 write(iout,*)

 write(6,*)'Written gwfn.data file.'

 close(io10)
 close(io12)
 close(iout)

 stop

1  call err1(0,'CRYSTALTOQMC','Allocation problem (1).')
2  call err1(0,'CRYSTALTOQMC','Allocation problem (2).')
3  call err1(0,'CRYSTALTOQMC','Allocation problem (3).')
4  call err1(0,'CRYSTALTOQMC','Allocation problem (4).')
5  call err1(0,'CRYSTALTOQMC','Allocation problem (5).')
6  call err1(0,'CRYSTALTOQMC','Allocation problem (6).')
7  call err1(0,'CRYSTALTOQMC','Allocation problem (7).')
8  call err1(0,'CRYSTALTOQMC','Allocation problem (8).')
10 call err1(0,'CRYSTALTOQMC','Allocation problem (10).')
12 call err1(0,'CRYSTALTOQMC','Error opening fort.12 file')
30 call err1(0,'CRYSTALTOQMC','Error opening fort.30 file')
100 call err1(0,'CRYSTALTOQMC','Error opening fort.10 file')

END PROGRAM crystaltoqmc

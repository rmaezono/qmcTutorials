!-------------------------------------------------------------------!
! Nonorthogonal wave function localization program                  !
! Dario Alfe` 2.2004                                                !
!                                                                   !
! Given a set of orbitals {psi} this program finds the linear       !
! combination that maximizes the localization of the orbitals on    !
! chosen spheres or parallelepipeds in the cell. The orbitals {psi} !
! do not need to be orthogonal.  The original and localized         !
! orbitals are both represented in a plane-wave basis.              !
!                                                                   !
! INPUT: files 'pwfn.data' and 'centres.dat'                        !
! OUTPUT: file 'pwfn.data.localized'                                !
!                                                                   !
! Format of 'centres.dat':                                          !
!                                                                   !
! Number of centres                                                 !
! 4                                                                 !
! Display coefficients of linear transformation (0=NO; 1=YES)       !
! 1                                                                 !
! Use spherical (1) or parallelepiped (2) localization regions      !
! 1                                                                 !
! x,y & z coords of centres ; radius ; no. orbs on centre (up & dn) !
! 0.0 0.0 0.0     5.0     2  2                                      !
! 0.4 1.3 2.9     3.0     1  1                                      !
! 8.0 3.9 1.4     5.0     2  2                                      !
! 7.2 3.4 5.6     6.0     3  3                                      !
!                                                                   !
! The plane-wave files are in the standard CASINO format.  Note     !
! only one k point is allowed.                                      !
!                                                                   !
! External routines: ZHEGV from lapack (gen. eigenproblem solver).  !
!                                                                   !
! Changes                                                           !
! =======                                                           !
! 12.2005  NDD  Tidied and rearranged.  Changed format of input.    !
! 12.2005  NDD  Fixed bugs for spin-polarized systems.              !
! 12.2005  NDD  Made it possible to omit some bands from transf.    !
!  1.2006  NDD  More tidying.  Parallelepiped localization regions. !
!               Added grid_mult.                                    !
!  8.2008  NDD  Warn if transformation is (nearly) singular.        !
!-------------------------------------------------------------------!

MODULE utils
 IMPLICIT NONE
 INTEGER,PARAMETER :: dp=kind(1.d0)
 REAL(dp),PARAMETER :: pi=3.14159265358979324d0,twopi=2.d0*pi, &
  &rec_twopi=1.d0/(2.d0*pi)


CONTAINS


 SUBROUTINE localize
!-----------------------------------------------------------------------------!
! This subroutine generates localized plane wave orbitals, by forming         !
! a linear combination of the orbitals in pwfn.data.  The plane-wave orbitals !
! are assumed to be complex, and only to exist at Gamma.                      !
!-----------------------------------------------------------------------------!
  USE singleton, ONLY : fftn
  IMPLICIT NONE
  INTEGER i,j,l,nr(3),ibnd,ibnd1,ibnd2,ig,npwx,l1,l2,l3,nks,nat,na,atomno, &
   &ibnd_new,nbnd(2),spin,num_spins,nelec,icount,nbnd_trans,ncentres,lwork, &
   &info,ntot(2),iprint,icut,ierr,ialloc,gmax(3),lower_lim_1,fft_ndata, &
   &upper_lim_1,lower_lim_2,upper_lim_2,lower_lim_3,upper_lim_3,l1p,l2p,l3p
  INTEGER,PARAMETER :: io_plwfn=10,io_pwfn=11,io_centres=12
  INTEGER,ALLOCATABLE :: norbs_at_centre(:,:),g(:,:)
  REAL(dp) r(3),at(3,3),bg(3,3),tau(3),energycpt,ecut,xk(3),r1(3),rdiff(3),&
   &radius2,tempr,gtemp(3)
  REAL(dp),ALLOCATABLE :: g2(:),centre(:,:),radius(:),lambda(:),rwork(:),&
   &et(:),radius_lv(:,:)
  COMPLEX(dp) :: det_tmat
  COMPLEX(dp),ALLOCATABLE :: evc(:,:),evcl(:),evcr(:,:,:,:),work(:), &
   &amat(:,:),bmat(:,:),bmatwork(:,:),fftdata(:),tmat(:,:)
  CHARACTER(80) char80
  LOGICAL lsda,pwreal
! Tolerance for comparing numbers etc.
  REAL(dp),PARAMETER :: eps=1.d-10
! Grid multiplicity, for increasing fineness of real-space integration grid.
  REAL(dp),PARAMETER :: grid_mult=1.d0

  INTERFACE
   SUBROUTINE zhegv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info)
    IMPLICIT NONE
    CHARACTER,INTENT(in) :: jobz,uplo
    INTEGER,INTENT(in) :: itype,lda,ldb,lwork,n
    COMPLEX(kind(1.d0)),INTENT(inout) :: a(lda,*),b(ldb,*)
    COMPLEX(kind(1.d0)),INTENT(out) :: work(*)
    INTEGER,INTENT(out) :: info
    REAL(kind(1.d0)),INTENT(out) :: rwork(*),w(*)
   END SUBROUTINE zhegv
  END INTERFACE

! Read in data from centres.dat
  open(io_centres,file='centres.dat',status='old',iostat=ierr)
  if(ierr/=0)then
   write(6,*)'Problem opening centres.dat.'
   stop
  endif ! ierr/=0
  call skipio(io_centres,1)
  read(io_centres,*,err=10,end=20)ncentres
  if(ncentres<=0)then
   write(6,*)'Number of centres must be positive.'
   stop
  endif ! ncentres<=0
  allocate(centre(3,ncentres),radius(ncentres),norbs_at_centre(ncentres,2), &
   &stat=ialloc)
  if(ialloc/=0)then
   write(6,*)'Allocation error.'
   stop
  endif ! ialloc/=0
  call skipio(io_centres,1)
  read(io_centres,*,err=10,end=20)iprint
  call skipio(io_centres,1)
  read(io_centres,*,err=10,end=20)icut
  if(icut/=1.and.icut/=2)then
   write(6,*)'Must choose either a spherical (1) or parallelepiped (2) &
    &localization region.'
   stop
  endif ! icut
  if(icut==2)then
   allocate(radius_lv(3,ncentres),stat=ialloc)
   if(ialloc/=0)then
    write(6,*)'Allocation error: radius_lv.'
    stop
   endif ! ialloc/=0
  endif ! icut=2
  ntot=0
  call skipio(io_centres,1)
  do i=1,ncentres
   read(io_centres,*,err=10,end=20)centre(:,i),radius(i),norbs_at_centre(i,1:2)
   if(radius(i)<=0.d0)then
    write(6,*)'Found a negative radius in centres.dat!'
    stop
   endif ! radius<=0
   if(any(norbs_at_centre(i,:)<0))then
    write(6,*)'A centre has a negative number of orbitals in centres.dat!'
    stop
   elseif(all(norbs_at_centre(i,:)==0))then
    write(6,*)'Warning: found a centre with no orbitals localized on it.'
   endif ! norbs_at_centre<0
   ntot=ntot+norbs_at_centre(i,:)
  enddo ! i
  close(io_centres)

  open(io_pwfn,file='pwfn.data',status='old',iostat=ierr)
  if(ierr/=0)then
   write(6,*)'Problem opening pwfn.data.'
   stop
  endif ! ierr/=0

  open(io_plwfn,file='pwfn.data.localized',status='replace',iostat=ierr)
  if(ierr/=0)then
   write(6,*)'Problem opening pwfn.data.localized.'
   stop
  endif ! ierr/=0

  read(io_pwfn,'(a)',err=30,end=40)char80
  char80=trim(adjustl(char80))//' (loc)'
  write(io_plwfn,'(a)')trim(char80)
  call skipio(io_pwfn,4)
  write(io_plwfn,*)
  write(io_plwfn,*)'BASIC INFO'
  write(io_plwfn,*)'----------'
  write(io_plwfn,*)'Generated by:'
  read(io_pwfn,'(a)',err=30,end=40)char80
  write(io_plwfn,*)trim(adjustl(char80))
  call skipio(io_pwfn,1)
  write(io_plwfn,*)'Method:'
  read(io_pwfn,'(a)',err=30,end=40)char80
  write(io_plwfn,*)trim(adjustl(char80))
  call skipio(io_pwfn,1)
  write(io_plwfn,*)'DFT Functional:'
  read(io_pwfn,'(a)',err=30,end=40)char80
  write(io_plwfn,*)trim(adjustl(char80))
  call skipio(io_pwfn,1)
  write(io_plwfn,*)'Pseudopotential'
  read(io_pwfn,'(a)',err=30,end=40)char80
  write(io_plwfn,*)trim(adjustl(char80))
  call skipio(io_pwfn,1)
  write(io_plwfn,*)'Plane wave cutoff (au)'
  read(io_pwfn,*,err=30,end=40)ecut
  write(io_plwfn,*)ecut
  call skipio(io_pwfn,1)
  write(io_plwfn,*)'Spin polarized:'
  read(io_pwfn,*,err=30,end=40)lsda
  write(io_plwfn,*)lsda
  call skipio(io_pwfn,1)
  write(io_plwfn,*)'Total energy (au per primitive cell)'
  read(io_pwfn,*,err=30,end=40)energycpt
  write(io_plwfn,*)energycpt
  call skipio(io_pwfn,1)
  write(io_plwfn,*)'Kinetic energy (au per primitive cell)'
  read(io_pwfn,*,err=30,end=40)energycpt
  write(io_plwfn,*)energycpt
  call skipio(io_pwfn,1)
  write(io_plwfn,*)'Local potential energy (au per primitive cell)'
  read(io_pwfn,*,err=30,end=40)energycpt
  write(io_plwfn,*)energycpt
  call skipio(io_pwfn,1)
  write(io_plwfn,*)'Non local potential energy(au per primitive cell)'
  read(io_pwfn,*,err=30,end=40)energycpt
  write(io_plwfn,*)energycpt
  call skipio(io_pwfn,3)
  write(io_plwfn,*)'Electron electron energy (au per primitive cell)'
  write(io_plwfn,*)0.d0
  write(io_plwfn,*)'Ion-ion energy (au per primitive cell)'
  read(io_pwfn,*,err=30,end=40)energycpt
  write(io_plwfn,*)energycpt
  call skipio(io_pwfn,1)
  write(io_plwfn,*)'Number of electrons per primitive cell'
  read(io_pwfn,*,err=30,end=40)nelec
  write(io_plwfn,*)nelec
  call skipio(io_pwfn,4)
  write(io_plwfn,*)
  write(io_plwfn,*)'GEOMETRY'
  write(io_plwfn,*)'-------- '
  write(io_plwfn,*)'Number of atoms per primitive cell '
  read(io_pwfn,*,err=30,end=40)nat
  write(io_plwfn,*)nat
  call skipio(io_pwfn,1)
  write(io_plwfn,*)'Atomic number and position of the atoms(au)'
  write(6,*)'Atomic number and position of the atoms(au):'
  do na=1,nat
   read(io_pwfn,*,err=30,end=40)atomno,tau(:)
   write(io_plwfn,'(i3,3(1x,es24.17))')atomno,tau(:)
   write(6,'(1x,i3,2x,3(1x,f20.12))')atomno,tau(:)
  enddo ! na
  write(6,*)
  call skipio(io_pwfn,1)
  write(io_plwfn,*)'Primitive lattice vectors (au)'
  write(6,*)'Primitive lattice vectors (au):'
  read(io_pwfn,*,err=30,end=40)at(:,1)
  read(io_pwfn,*,err=30,end=40)at(:,2)
  read(io_pwfn,*,err=30,end=40)at(:,3)
  write(io_plwfn,'(1x,3(es24.17,1x))')at(:,1)
  write(io_plwfn,'(1x,3(es24.17,1x))')at(:,2)
  write(io_plwfn,'(1x,3(es24.17,1x))')at(:,3)
  write(6,'(1x,3(f20.13,1x))')at(:,1)
  write(6,'(1x,3(f20.13,1x))')at(:,2)
  write(6,'(1x,3(f20.13,1x))')at(:,3)
  write(6,*)

  if(.not.lsda.and.any(norbs_at_centre(:,2)/=norbs_at_centre(:,1)))then
   write(6,*)'Have requested different localisations for up and down spins.'
   write(6,*)'However, your data is supposed to be non-spin-polarised.'
   stop
  endif ! non-spin polar, but different localization for up & down.

  if(icut==2)then
   write(6,*)'Parallelepiped cutoff'
  else
   write(6,*)'Spherical cutoff'
  endif ! icut
  write(6,*)

! Calculate reciprocal lattice vectors (columns of bg), without factor of 2pi.
  call inve(at,bg)
  bg=transpose(bg)

! Work out "radii" of centres in terms of lattice vectors if a parallelepiped
! shaped cell is to be used.
  if(icut==2)then
   do i=1,ncentres
    radius_lv(1,i)=radius(i)*sqrt(dot_product(bg(:,1),bg(:,1)))
    radius_lv(2,i)=radius(i)*sqrt(dot_product(bg(:,2),bg(:,2)))
    radius_lv(3,i)=radius(i)*sqrt(dot_product(bg(:,3),bg(:,3)))
   enddo ! i
  endif ! icut=2

  call skipio(io_pwfn,4)
  write(io_plwfn,*)
  write(io_plwfn,*)'G VECTORS'
  write(io_plwfn,*)'---------'
  write(io_plwfn,*)'Number of G-vectors'
  read(io_pwfn,*,err=30,end=40)npwx
  write(io_plwfn,*)npwx
  call skipio(io_pwfn,1)
  write(io_plwfn,*)'Gx Gy Gz (au)'

  allocate(g(3,npwx),g2(npwx),stat=ialloc)
  if(ialloc/=0)then
   write(6,*)'Allocation error: G vectors.'
   stop
  endif ! ialloc/=0
  gmax=0
  do ig=1,npwx
   read(io_pwfn,*,err=30,end=40)gtemp
   write(io_plwfn,'(1x,es24.17,1x,es24.17,1x,es24.17)')gtemp
   g2(ig)=dot_product(gtemp,gtemp)
   g(:,ig)=nint(matmul(gtemp,at)*rec_twopi)
   if(abs(g(1,ig))>gmax(1))gmax(1)=abs(g(1,ig))
   if(abs(g(2,ig))>gmax(2))gmax(2)=abs(g(2,ig))
   if(abs(g(3,ig))>gmax(3))gmax(3)=abs(g(3,ig))
  enddo ! ig
  nr=2*ceiling(dble(gmax)*grid_mult)+2
  write(6,*)'Number of plane waves: '//trim(i2s(npwx))
  write(6,*)'FFT grid             : '//trim(i2s(nr(1)))//' by ' &
   &//trim(i2s(nr(2)))//' by '//trim(i2s(nr(3)))
  write(6,*)

  call skipio(io_pwfn,4)
  write(io_plwfn,*)
  write(io_plwfn,*)'WAVE FUNCTION'
  write(io_plwfn,*)'-------------'
  write(io_plwfn,*)'Number of k-points'
  read(io_pwfn,*,err=30,end=40)nks
  if(nks/=1)then
   write(6,*)'There should only be one k point in pwfn.data: Gamma.'
   stop
  endif ! nks/=1
  write(io_plwfn,*)nks
  if(lsda)then
   num_spins=2
  else
   num_spins=1
  endif ! num_spins

  fft_ndata=nr(1)*nr(2)*nr(3)
  allocate(fftdata(fft_ndata),stat=ialloc)
  if(ialloc/=0)then
   write(6,*)'Allocation error: FFT data.'
   stop
  endif ! ialloc

  if(any(g(:,1)/=0))then
   write(6,*)'First G vector in pwfn.data must be 0.'
   stop
  endif ! G(1)/=0

  call skipio(io_pwfn,1)
  write(io_plwfn,*)'k-point # ; # of bands (up spin/down spin) ; &
   &k-point coords (au)'
  read(io_pwfn,*,err=30,end=40)i,nbnd(1),nbnd(2),xk(:)
  if(abs(xk(1))>eps.or.abs(xk(2))>eps.or.abs(xk(3))>eps)then
   write(6,*)'Non-Gamma k points are not allowed at present.'
   stop
  endif ! k/=0
  if(any(ntot(1:num_spins)>nbnd(1:num_spins)))then
   write(6,*)'Have requested more orbitals to be localized than are available!'
   stop
  endif ! ntot>nbnd

  write(io_plwfn,'(1x,i1,1x,1x,2(i5,1x),3(es24.17,1x))')1,nbnd(1:2),xk(:)

  do spin=1,num_spins

   if(lsda)then
    write(6,*)'SPIN '//trim(i2s(spin))
    write(6,*)'======'
    write(6,*)
   endif ! spin_polarized

! Number of bands included in the transformation
   nbnd_trans=min(nbnd(spin),ntot(spin))

   lwork=2*nbnd_trans-1
   allocate(evcr(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1,nbnd_trans),et(nbnd_trans),&
    &evc(npwx,nbnd_trans),evcl(npwx),lambda(nbnd_trans), &
    &amat(nbnd_trans,nbnd_trans),bmat(nbnd_trans,nbnd_trans), &
    &bmatwork(nbnd_trans,nbnd_trans),work(lwork),rwork(3*nbnd_trans-2), &
    &tmat(nbnd_trans,nbnd_trans),stat=ialloc)
   if(ialloc/=0)then
    write(6,*)'Allocation problem: work arrays.'
    stop
   endif ! ialloc/=0

   do ibnd=1,nbnd_trans
    call skipio(io_pwfn,1)
    read(io_pwfn,*,err=30,end=40)i,i,et(ibnd)
    call skipio(io_pwfn,1)
! Read in the wavefunctions.  Treat them as complex, even if they happen
! to be real.
    read(io_pwfn,'(a)',err=30,end=40)char80
    if(index(char80,',')>0)then
     pwreal=.false.
    else
     pwreal=.true.
    endif ! Coefficient contains ",".
    backspace(io_pwfn)
    if(pwreal)then
     do ig=1,npwx
      read(io_pwfn,*,err=30,end=40)tempr
      evc(ig,ibnd)=cmplx(tempr,0.d0,dp)
     enddo ! ig
    else
     do ig=1,npwx
      read(io_pwfn,*,err=30,end=40)evc(ig,ibnd)
     enddo ! ig
    endif ! pwreal
   enddo ! ibnd

   do ibnd=1,nbnd_trans

    fftdata=cmplx(0.d0,0.d0,dp)
    do ig=1,npwx
     icount=1+g(1,ig)+nr(1)*(indexfn(g(1,ig))+g(2,ig)+nr(2)*(indexfn(g(2,ig)) &
      &+g(3,ig)+nr(3)*indexfn(g(3,ig))))
     fftdata(icount)=evc(ig,ibnd)
    enddo ! ig

! Transform to real space
    call fftn(fftdata,nr,inv=.true.)
    fftdata=fftdata*sqrt(dble(fft_ndata))
    icount=1
    do l3=0,nr(3)-1
     do l2=0,nr(2)-1
      do l1=0,nr(1)-1
       evcr(l1,l2,l3,ibnd)=fftdata(icount)
       icount=icount+1
      enddo ! l1
     enddo ! l2
    enddo ! l3

   enddo  ! bands

! Calculate bmat(ij) = <psi_i|psi_j>.
   do ibnd2=1,nbnd_trans
    bmat(ibnd2,ibnd2)=cmplx(0.d0,0.d0,dp)
    do ig=1,npwx
     bmat(ibnd2,ibnd2)=bmat(ibnd2,ibnd2)+conjg(evc(ig,ibnd2))*evc(ig,ibnd2)
    enddo ! ig
    do ibnd1=ibnd2+1,nbnd_trans
     bmat(ibnd1,ibnd2)=cmplx(0.d0,0.d0,dp)
     do ig=1,npwx
      bmat(ibnd1,ibnd2)=bmat(ibnd1,ibnd2)+conjg(evc(ig,ibnd1))*evc(ig,ibnd2)
     enddo ! ig
     bmat(ibnd2,ibnd1)=conjg(bmat(ibnd1,ibnd2))
    enddo ! ibnd1
   enddo ! ibnd2

! Now for each centre...
   ibnd_new=0
   do i=1,ncentres
    radius2=radius(i)*radius(i)
    if(4.d0*radius2*dot_product(bg(:,1),bg(:,1))>1.d0.or.&
     &4.d0*radius2*dot_product(bg(:,2),bg(:,2))>1.d0.or.&
     &4.d0*radius2*dot_product(bg(:,3),bg(:,3))>1.d0)then
     write(6,*)'Localization radius is too big for simulation cell.'
     stop
    endif ! radius extends outside cell.
    amat=cmplx(0.d0,0.d0,dp)

    if(icut==1)then
! Spherical cutoff.

     do l1=0,nr(1)-1
      r1(1)=dble(l1)/dble(nr(1))
      do l2=0,nr(2)-1
       r1(2)=dble(l2)/dble(nr(2))
       do l3=0,nr(3)-1
        r1(3)=dble(l3)/dble(nr(3))
        r=matmul(at,r1)-centre(:,i)
        call min_image(r,at,bg,rdiff)
! Calculate the matrix amat(ij) = <psi_i|psi_j>_gamma, where <|>_gamma
! is the integral over the sphere centred on centre(:,i) with radius radius(i)
        if(dot_product(rdiff,rdiff)<radius2)then
         do ibnd2=1,nbnd_trans
          do ibnd1=1,ibnd2
           amat(ibnd1,ibnd2)=amat(ibnd1,ibnd2)+conjg(evcr(l1,l2,l3,ibnd1)) &
            &*evcr(l1,l2,l3,ibnd2)
          enddo ! ibnd1
         enddo ! ibnd2
        endif
       enddo ! l3
      enddo ! l2
     enddo ! l1

    else
! Parallelepiped cutoff.

     lower_lim_1=ceiling((dot_product(centre(:,i),bg(:,1))-radius_lv(1,i)) &
      &*dble(nr(1)))
     upper_lim_1=floor((dot_product(centre(:,i),bg(:,1))+radius_lv(1,i)) &
      &*dble(nr(1)))
     lower_lim_2=ceiling((dot_product(centre(:,i),bg(:,2))-radius_lv(2,i)) &
      &*dble(nr(2)))
     upper_lim_2=floor((dot_product(centre(:,i),bg(:,2))+radius_lv(2,i)) &
      &*dble(nr(2)))
     lower_lim_3=ceiling((dot_product(centre(:,i),bg(:,3))-radius_lv(3,i)) &
      &*dble(nr(3)))
     upper_lim_3=floor((dot_product(centre(:,i),bg(:,3))+radius_lv(3,i)) &
      &*dble(nr(3)))
     do l1p=lower_lim_1,upper_lim_1
      l1=modulo(l1p,nr(1))
      do l2p=lower_lim_2,upper_lim_2
       l2=modulo(l2p,nr(2))
       do l3p=lower_lim_3,upper_lim_3
        l3=modulo(l3p,nr(3))
! Calculate the matrix amat(ij) = <psi_i|psi_j>_gamma, where <|>_gamma
! is the integral over the sphere centred on centre(:,i) with radius radius(i)
        do ibnd2=1,nbnd_trans
         do ibnd1=1,ibnd2
          amat(ibnd1,ibnd2)=amat(ibnd1,ibnd2)+conjg(evcr(l1,l2,l3,ibnd1)) &
           &*evcr(l1,l2,l3,ibnd2)
         enddo ! ibnd1
        enddo ! ibnd2
       enddo ! l3
      enddo ! l2p
     enddo ! l1p

    endif ! icut=1
    amat=amat/dble(nr(1)*nr(2)*nr(3))

! Solve the generalized eigenvalue problem.
    bmatwork=bmat
    call zhegv(1,'v','u',nbnd_trans,amat,nbnd_trans,bmatwork,nbnd_trans, &
     &lambda,work,lwork,rwork,info)
    if(info>0)then
     write(6,*)'WARNING, zhegv info = '//trim(i2s(info))
     write(6,*)
    endif ! info>0

! Note that new orbitals are the ones with the highest eigenvalues; hence
! we count down from nbnd_trans when assigning orbitals to a particular centre.
    write(6,*)'Centre '//trim(i2s(i))
    write(6,*)'-----------'
    write(6,'(" ",a,3(f10.6,1x))')'Position    : ',centre(:,i)
    write(6,*)'No. orbitals: ',norbs_at_centre(i,spin)
    write(6,*)'Trunc. rad. : ',radius(i)
    write(6,*)
    write(6,*)'Localization fraction (eigenvalue) for each orbital at centre:'
    do j=nbnd_trans,nbnd_trans-norbs_at_centre(i,spin)+1,-1
     write(6,*)'  Orbital '//trim(i2s(nbnd_trans+1-j))//'  ',lambda(j)
    enddo ! j
    write(6,*)
    do j=1,norbs_at_centre(i,spin)
     ibnd_new=ibnd_new+1
     if(iprint/=0)then
      write(6,*)'Coefficients of linear combination for orbital # ' &
       &//trim(i2s(j))//' of centre '//trim(i2s(i))//':'
      do l=1,nbnd_trans
       write(6,*)'  Coeff. of orig. band '//trim(i2s(l))//'  ', &
        &amat(l,nbnd_trans-j+1)
      enddo ! l
     endif ! iprint/=0
     tmat(1:nbnd_trans,ibnd_new)=amat(1:nbnd_trans,nbnd_trans-j+1)
     evcl=cmplx(0.d0,0.d0,dp)
     do l=1,nbnd_trans
      do ig=1,npwx
       evcl(ig)=evcl(ig)+amat(l,nbnd_trans-j+1)*evc(ig,l)
      enddo ! ig
     enddo ! l
     write(io_plwfn,*)'Band, spin, eigenvalue (au), centre, radius, orb, &
      &norm^2'
     write(io_plwfn,'(1x,i5,1x,i1,1x,5(es24.17,1x),i3,1x,es24.17)')ibnd_new, &
      &spin,et(ibnd_new),centre(:,i),radius(i),j,lambda(nbnd_trans-j+1)
     write(io_plwfn,*)'Eigenvector coefficients'
     do ig=1,npwx
      write(io_plwfn,*)evcl(ig)
     enddo ! ig
    enddo ! j
    if(iprint/=0)write(6,*)
   enddo ! i

! Work out if tmat (matrix of transformation) is singular.
   call eval_det(nbnd_trans,tmat,det_tmat)
   if(abs(det_tmat)<1.d-7)then
    write(6,*)'Warning: determinant of orbital-transformation matrix is small.'
    write(6,*)'Determinant: ',det_tmat
    write(6,*)'Perhaps consider using different localization regions.'
    write(6,*)
   endif ! |T| is small.

   deallocate(evc,et,evcl,evcr,amat,bmat,bmatwork,lambda,work,rwork,tmat)

! Copy out the eigenfunctions for bands not included in the transformation.
   do ibnd=nbnd_trans+1,nbnd(spin)
    call skipio(io_pwfn,1)
    write(io_plwfn,*)'Band, spin, eigenvalue (au)'
    read(io_pwfn,*,err=30,end=40)i,j,tempr
    write(io_plwfn,'(1x,i5,1x,i1,1x,es24.17)')i,j,tempr
    call skipio(io_pwfn,1)
    write(io_plwfn,*)'Eigenvector coefficients'
    do ig=1,npwx
     read(io_pwfn,'(a)',err=30,end=40)char80
     write(io_plwfn,*)trim(adjustl(char80))
    enddo ! ig
   enddo ! ibnd

  enddo ! spin

  deallocate(fftdata,norbs_at_centre,centre,radius,g,g2)
  if(allocated(radius_lv))deallocate(radius_lv)

  close(io_pwfn)
  close(io_plwfn)

  write(6,*)'Have generated pwfn.data.localized file.'
  write(6,*)

  return

10 write(6,*)'Error reading centres.dat.'
  stop
20 write(6,*)'Unexpectedly reached end of centres.dat.'
  stop
30 write(6,*)'Error reading pwfn.data.'
  stop
40 write(6,*)'Unexpectedly reached end of pwfn.data.'
  stop
 END SUBROUTINE localize


 SUBROUTINE skipio(io,n)
!----------------------------------------------------------------!
! This subroutine skips n lines from the file opened on unit io. !
!----------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: io,n
  INTEGER i,ierr
  do i=1,n
   read(io,*,iostat=ierr)
   if(ierr<0)then
    write(6,*)'Unexpectedly reached end of file whilst skipping lines.'
    stop
   elseif(ierr>0)then
    write(6,*)'Error reading file whilst skipping lines.'
    stop
   endif ! ierr
  enddo ! i
 END SUBROUTINE skipio


 INTEGER FUNCTION indexfn(i)
!-----------------------------------------------------------!
! This function returns 1 if i is negative and 0 otherwise. !
!-----------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: i
  if(i>=0)then
   indexfn=0
  else
   indexfn=1
  endif ! i>=0
 END FUNCTION indexfn


 SUBROUTINE min_image(a,latt_vect,rec_vect,b)
!--------------------------------------------------------------!
! This subroutine computes b as the minimum image vector of    !
! vector a with respect to the lattice specified by latt_vect. !
! So -b is the vector from a to its closest lattice point.     !
!--------------------------------------------------------------!
  IMPLICIT NONE
  REAL(dp),INTENT(in) :: a(3),latt_vect(3,3),rec_vect(3,3)
  REAL(dp),INTENT(out) :: b(3)
  REAL(dp) Delta(3),mag_b_sq,dist_sq
  INTEGER n0(3),n(3),i,j,k
  INTEGER,PARAMETER :: check_shell=2
  n0(1)=floor(dot_product(a(1:3),rec_vect(1:3,1)))
  n0(2)=floor(dot_product(a(1:3),rec_vect(1:3,2)))
  n0(3)=floor(dot_product(a(1:3),rec_vect(1:3,3)))
  mag_b_sq=a(1)**2+a(2)**2+a(3)**2
  b(1:3)=a(1:3)
  do i=-check_shell,check_shell+1
   n(1)=n0(1)+i
   do j=-check_shell,check_shell+1
    n(2)=n0(2)+j
    do k=-check_shell,check_shell+1
     n(3)=n0(3)+k
     Delta(1:3)=a(1:3)-n(1)*latt_vect(1:3,1)-n(2)*latt_vect(1:3,2) &
      &-n(3)*latt_vect(1:3,3)
     dist_sq=Delta(1)**2+Delta(2)**2+Delta(3)**2
     if(dist_sq<mag_b_sq)then
      mag_b_sq=dist_sq
      b(1:3)=Delta(1:3)
     endif ! dist_sq<mag_b_sq
    enddo ! k
   enddo ! j
  enddo ! i
 END SUBROUTINE min_image


 SUBROUTINE inve(v,inv)
!----------------------!
! Inverts 3x3 matrices !
!----------------------!
  IMPLICIT NONE
  REAL(dp),INTENT(in) :: v(3,3)
  REAL(dp),INTENT(out) :: inv(3,3)
  REAL(dp) d
  d=v(1,1)*(v(2,2)*v(3,3)-v(2,3)*v(3,2))+ &
   &v(2,1)*(v(3,2)*v(1,3)-v(1,2)*v(3,3))+ &
   &v(3,1)*(v(1,2)*v(2,3)-v(1,3)*v(2,2))
  if(d==0.d0)then
   write(6,*)'Trying to invert a singular matrix...'
   stop
  endif ! d=0
  d=1.d0/d
  inv(1,1)=(v(2,2)*v(3,3)-v(2,3)*v(3,2))*d
  inv(1,2)=(v(3,2)*v(1,3)-v(1,2)*v(3,3))*d
  inv(1,3)=(v(1,2)*v(2,3)-v(1,3)*v(2,2))*d
  inv(2,1)=(v(3,1)*v(2,3)-v(2,1)*v(3,3))*d
  inv(2,2)=(v(1,1)*v(3,3)-v(3,1)*v(1,3))*d
  inv(2,3)=(v(2,1)*v(1,3)-v(1,1)*v(2,3))*d
  inv(3,1)=(v(2,1)*v(3,2)-v(2,2)*v(3,1))*d
  inv(3,2)=(v(3,1)*v(1,2)-v(1,1)*v(3,2))*d
  inv(3,3)=(v(1,1)*v(2,2)-v(1,2)*v(2,1))*d
 END SUBROUTINE inve


 CHARACTER(12) FUNCTION i2s(n)
!-----------------------------------------------------------------------!
! Convert integers to left justified strings that can be printed in the !
! middle of a sentence without introducing large amounts of white space.!
!-----------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: n
  INTEGER i,j
  INTEGER,PARAMETER :: ichar0=ichar('0')
  i2s=''
  i=abs(n)
  do j=len(i2s),1,-1
   i2s(j:j)=achar(ichar0+mod(i,10))
   i=i/10 ; if(i==0)exit
  enddo ! j
  if(n<0)then
   i2s='-'//adjustl(i2s)
  else
   i2s=adjustl(i2s)
  endif ! n<0
 END FUNCTION i2s


 SUBROUTINE lu_decom_cmplx(a,piv,n)
!------------------------------------!
! LU decomposition. A overwritten on !
! output. Complex version.           !
!------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n
 COMPLEX(dp),INTENT(inout) :: a(n,n)
 INTEGER,INTENT(inout) :: piv(n)
 INTEGER ierr
 INTERFACE
  SUBROUTINE zgetrf(m,n,a,lda,ipiv,info)
   IMPLICIT NONE
   INTEGER,INTENT(in) :: lda,m,n
   COMPLEX(kind(1.d0)),INTENT(inout) :: a(lda,*)
   INTEGER,INTENT(out) :: info,ipiv(*)
  END SUBROUTINE zgetrf
 END INTERFACE
 call zgetrf(n,n,a,n,piv,ierr)
 if(ierr==0)return
 if(ierr<0)write(6,*)'ZGETRF says parameter #',trim(i2s(-ierr)),' has an &
  &illegal value.'
 if(ierr>0)write(6,*)'ZGETRF says the ',trim(i2s(ierr)),' diagonal element &
  &is exactly zero.'
 stop
 END SUBROUTINE lu_decom_cmplx


 SUBROUTINE eval_det(n,a,det)
!----------------------------------------------------------------------!
! This subroutine evaluates the determinant of a matrix a, using an LU !
! decomposition algorithm.                                             !
!----------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: n
  COMPLEX(dp),INTENT(in) :: a(n,n)
  COMPLEX(dp),INTENT(out) :: det
  COMPLEX(dp) bjj
  COMPLEX(dp),ALLOCATABLE :: b(:,:)
  INTEGER j,ialloc
  INTEGER,ALLOCATABLE :: piv(:)
  LOGICAL flip
  allocate(b(n,n),piv(n),stat=ialloc)
  if(ialloc/=0)then
   write(*,*)'Allocation problem in EVAL_DET.'
   stop
  endif ! ialloc
  b=a
  call lu_decom_cmplx(b,piv,n)
  flip=.false.
  det=cmplx(1.d0,0.d0,dp)
  do j=1,n
   bjj=b(j,j)
   det=det*bjj
   if(j/=piv(j))flip=.not.flip
  enddo ! j
  if(flip)det=-det
  deallocate(b)
 END SUBROUTINE eval_det


END MODULE utils


PROGRAM localizer
!---------------------------!
! Main program starts here. !
!---------------------------!
 USE utils
 IMPLICIT NONE

 write(6,*)
 write(6,*)'LOCALIZER'
 write(6,*)'========='
 write(6,*)

 call localize

 write(6,*)'Program finished.'
 write(6,*)

END PROGRAM localizer



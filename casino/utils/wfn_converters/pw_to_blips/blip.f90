!-------------------------------------------------------------------!
! Plane waves -> blip conversion utility                            !
! Dario Alfe` 3.2001                                                !
!                                                                   !
! This program requires a CASINO-format pwfn.data file containing   !
! plane-wave orbitals.  It produces a CASINO-format bwfn.data file  !
! with the same orbitals expressed in terms of blip functions.  If  !
! the user supplies a centres.dat file then the blip orbitals will  !
! be truncated at the edge of the specified localization regions.   !
! The format of the centres.dat file is:                            !
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
! Skin thickness (au)                <- OPTIONAL LINE               !
! 0.0                                <- OPTIONAL LINE               !
!                                                                   !
! The centres.dat file will only be read if the orbitals exist at   !
! Gamma only.  Note that if the orbitals exist at Gamma only then   !
! real orbitals will be constructed, whereas if other k points are  !
! present then the original complex orbitals will be written out and!
! the construction of real orbitals will take place in CASINO       !
! itself.                                                           !
!                                                                   !
! Changes                                                           !
! -------                                                           !
!  5.2002 MDT - Tidied up. Made conform to f90 standard.            !
! 11.2003 DA  - Added fast fourier transforms                       !
! 10.2004 NDD - Tried to make the program a bit more user-friendly. !
! 10.2004 NDD - Some modifications for spin-polarized systems.      !
!  1.2005 NDD - Extensive changes.  Merged with blipl.              !
! 10.2006 NDD - Added alternative approximation for gamma array.    !
!  7.2007 NDD - Fixed bug in overlap test for large parallelepiped  !
!               trunc. localised orbitals (don't use min_image).    !
!               Some tidying.                                       !
!  2.2008 NDD - Report error bars on overlaps.                      !
!  8.2008 NDD - New random number generator.                        !
!  1.2010 NN  - Single-pass reading of pwfn.data input file         !
! 12.2011 NDD - Translate molecules into centre of blip grid.       !
!-------------------------------------------------------------------!

MODULE blip_utils
!---------------------------------------------------------------!
! Various subroutines for performing PW -> blip transformation. !
!---------------------------------------------------------------!
 USE helpers
 IMPLICIT NONE

! Runtime parameters, may be set via blip.in file
! Grid multiplicity
 REAL(dp) :: xmul=1.d0
! Number of random points for overlap test.
 INTEGER :: n_points_for_test=0
! Evaluate KE and norm2?
 LOGICAL :: norm2_ke_calc=.false.
! Number of times to carry out overlap test (to get error bars).
 INTEGER :: n_overlap_tests=12
! Which approx for "gamma" (see the blip paper / NDD's notes) to be used?
! 1=the approximation suggested in the blip paper; 2=NDD's approx.
 INTEGER :: gamma_approx=1
! Intended periodicity.  If this is less than 3 then the atoms will be
! translated to lie close to the centre of the blip grid.
 INTEGER :: intended_periodicity=3
 NAMELIST /blipin/ &
  & xmul,              &
  & n_points_for_test, &
  & norm2_ke_calc,     &
  & n_overlap_tests,   &
  & gamma_approx
! Blip grid dimensions (integer and real, and squared)
 INTEGER nr(3)
 REAL(dp) rnr(3),rnr2(3)
! Number of G vectors, number of k vectors and number of spin types
 INTEGER NGvec,Nkvec,no_spin_types
! Unit numbers for I/O.
 INTEGER,PARAMETER :: io_bwfn=10,io_pwfn=11,io_centres=12
! G vectors
 INTEGER,ALLOCATABLE :: g(:,:)
! Number of bands for a given k point and spin
 INTEGER,ALLOCATABLE :: nband(:,:)
 REAL(dp),ALLOCATABLE :: kvec(:,:)
 REAL(dp),ALLOCATABLE :: g2(:)
! Lattice vectors and reciprocal lattice vectors (w/o 2pi) in columns.
 REAL(dp) at(3,3),bg(3,3)
! Magnitudes of reciprocal lattice vectors
 REAL(dp) mag_rlv(3)
! Products of reciprocal lattice vectors
 REAL(dp) lvp(6)
! gamma array.
 REAL(dp),ALLOCATABLE :: gamma(:)
! Real / complex blip coefficients
 REAL(dp),ALLOCATABLE :: avc(:,:,:)
 COMPLEX(dp),ALLOCATABLE :: cavc(:,:,:)
! Is the system spin-polarized?
 LOGICAL blipreal,spin_polarized
! Variables to pass to external functions for integration.
 REAL(dp) centrepos_lv_pass(3),r_pass,r_sin_theta_pass, &
  &r_cos_theta_pass,r2_sin_theta_pass,rvec_lv_pass(3),rlv_min_pass(3), &
  &rlv_max_pass(3),kvec_pass(3),kvec2_pass
! Offset to atom positions
 REAL(dp) :: R_offset(3)
! Array of exp(-iG.R_offset).
 COMPLEX(dp),ALLOCATABLE :: expmigdotR_offset(:)


 PUBLIC


 CONTAINS


 SUBROUTINE get_args
!----------------------------------------!
! Ask the user for the input parameters. !
!----------------------------------------!
  IMPLICIT NONE
  INTEGER ierr
  CHARACTER(1) char1
  LOGICAL blipin_exists

! Read in from a blip.in file.  Why not just use "cat input_parameters > blip"?
  inquire(file='blip.in',exist=blipin_exists)
  if(blipin_exists)then
   write(6,*)'Reading options from blip.in'
   open(8,file='blip.in',status="old",iostat=ierr)
   if(ierr/=0)call errstop('GET_ARGS','Error opening blip.in.')
   read(8,blipin,iostat=ierr)
   if(ierr/=0)call errstop('GET_ARGS','Error reading blip.in.')
   close(8)
   write(6,*)
   write(6,*)'Using the following settings:'
   write(6,blipin)
   write(6,*)
   return
  endif ! blipin_exists

  do
   write(6,*)'Please enter the grid multiplicity (suggested values >=1.0).'
   read(5,*,iostat=ierr)xmul
   if(ierr/=0)xmul=-1.d0
   if(xmul<=0.d0)then
    write(6,*)'Please choose a positive number.'
    write(6,*)
   elseif(xmul<1.d0)then
    write(6,*)'Warning: usually best to choose a multiplicity greater than 1.'
    exit
   else
    exit
   endif ! xmul
  enddo
  write(6,*)
  do
   write(6,*)'Number of random points for overlap test (choose 0 to skip test).'
   read(5,*,iostat=ierr)n_points_for_test
   if(ierr/=0)n_points_for_test=-1
   if(n_points_for_test>=0)exit
   write(6,*)'Please try again.  Choose a non-negative number.'
  enddo
  write(6,*)
  do
   write(6,*)'Calculate the KE of the blip orbitals? (y/n)'
   read(5,*,iostat=ierr)char1
   if(ierr/=0)char1='z'
   if(char1=='y'.or.char1=='Y')then
    norm2_ke_calc=.true.
    exit
   elseif(char1=='n'.or.char1=='N')then
    norm2_ke_calc=.false.
    exit
   else
    write(6,*)'Please try again.'
   endif ! ierr
  enddo
  write(6,*)
  do
   call wordwrap('Translate the atoms into the centre of the blip region (for &
    &reduced-periodicity calculations)?  (y/n)')
   read(5,*,iostat=ierr)char1
   if(ierr/=0)char1='z'
   if(char1=='y'.or.char1=='Y')then
    write(6,*)
    call wordwrap('Please enter the intended periodicity (0, 1, 2 or 3).  &
     &Note that this program merely translates the atoms and orbitals; to &
     &use the intended periodicity in CASINO you must specify the &
     &BLIP_PERIODICITY keyword in the input file.')
    do
     read(5,*,iostat=ierr)intended_periodicity
     if(ierr/=0)intended_periodicity=-1
     if(intended_periodicity>=0.and.intended_periodicity<=3)exit
     write(6,*)'Please enter an integer between 0 and 3.'
    enddo
    exit
   elseif(char1=='n'.or.char1=='N')then
    intended_periodicity=3
    exit
   else
    write(6,*)'Please try again.'
   endif ! ierr
  enddo
  write(6,*)

  write(6,*)'Using the following settings:'
  write(6,*)'  Grid multiplicity           : ',xmul
  write(6,*)'  Number of overlap tests     : '//trim(i2s(n_overlap_tests))
  write(6,*)'  Number of points for test   : '//trim(i2s(n_points_for_test))
  write(6,*)'  Norm^2 & kinetic-energy test: '//l2s(norm2_ke_calc)
  write(6,*)'  Choice of blip method       : '//trim(i2s(gamma_approx))
  write(6,*)'  Intended periodicity        : '//trim(i2s(intended_periodicity))
  write(6,*)

 END SUBROUTINE get_args


 SUBROUTINE initialize
!-----------------------------!
! Open input and output file  !
!-----------------------------!
  IMPLICIT NONE
  INTEGER ierr
  open(io_pwfn,file='pwfn.data',status='old',action='read',iostat=ierr)
  if(ierr/=0)call errstop('INITIALIZE', &
   &'Problem opening pwfn.data for reading. Stopping.')
  open(io_bwfn,file='bwfn.data',status='unknown',action='write',iostat=ierr)
  if(ierr/=0)call errstop('INITIALIZE','Problem opening bwfn.data for writing. &
   &Stopping.')
 END SUBROUTINE initialize


 SUBROUTINE finalize
!-----------------------------!
! Close input and output file !
!-----------------------------!
  IMPLICIT NONE
  close(io_bwfn)
  close(io_pwfn)
  if(allocated(kvec))deallocate(kvec)
  if(allocated(nband))deallocate(nband)
 END SUBROUTINE finalize


 SUBROUTINE bwfn_header
!-----------------------------------------------------------------------!
! In this subroutine, information from the top of the pwfn.data file is !
! transferred into the bwfn.data file.  The geometry and G vectors are  !
! also read in.  The gamma array is constructed.                        !
!-----------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER nbasis,na,nelec,gmax(3),ialloc,ig,sum_atomno,i
  REAL(dp) da(3),k,k2,k4,cosk,molcentre(3)
  CHARACTER(80) title,code,method,functional,pseudo_type
  REAL(dp) plane_wave_cutoff,total_energy,kinetic_energy,&
   &local_potential_energy,non_local_potential_energy,electron_electron_energy,&
   &eionion,lvdotR_offset
  INTEGER,ALLOCATABLE :: atomno(:)
  REAL(dp),ALLOCATABLE :: atompos(:,:),g_raw(:,:)
  COMPLEX(dp),ALLOCATABLE :: gfac1(:),gfac2(:),gfac3(:)

  read(io_pwfn,'(a)',err=10,end=20)title                 ;call skipio(io_pwfn,4)
  read(io_pwfn,'(a)',err=10,end=20)code                  ;call skipio(io_pwfn,1)
  read(io_pwfn,'(a)',err=10,end=20)method                ;call skipio(io_pwfn,1)
  read(io_pwfn,'(a)',err=10,end=20)functional            ;call skipio(io_pwfn,1)
  read(io_pwfn,'(a)',err=10,end=20)pseudo_type           ;call skipio(io_pwfn,1)
  read(io_pwfn,*,err=10,end=20)plane_wave_cutoff         ;call skipio(io_pwfn,1)
  read(io_pwfn,*,err=10,end=20)spin_polarized            ;call skipio(io_pwfn,1)
  read(io_pwfn,*,err=10,end=20)total_energy              ;call skipio(io_pwfn,1)
  read(io_pwfn,*,err=10,end=20)kinetic_energy            ;call skipio(io_pwfn,1)
  read(io_pwfn,*,err=10,end=20)local_potential_energy    ;call skipio(io_pwfn,1)
  read(io_pwfn,*,err=10,end=20)non_local_potential_energy;call skipio(io_pwfn,1)
  read(io_pwfn,*,err=10,end=20)electron_electron_energy  ;call skipio(io_pwfn,1)
  read(io_pwfn,*,err=10,end=20)eionion                   ;call skipio(io_pwfn,1)
  read(io_pwfn,*,err=10,end=20)nelec                     ;call skipio(io_pwfn,4)
  read(io_pwfn,*,err=10,end=20)nbasis
  allocate(atomno(nbasis),atompos(3,nbasis),stat=ialloc)
  if(ialloc/=0)call errstop('BWFN_HEADER','Allocation problem: atoms.')
  call skipio(io_pwfn,1)
  do na=1,nbasis
   read(io_pwfn,*,err=10,end=20)atomno(na),atompos(:,na)
  enddo ! na
  call skipio(io_pwfn,1)
  read(io_pwfn,*,err=10,end=20)at(:,1)
  read(io_pwfn,*,err=10,end=20)at(:,2)
  read(io_pwfn,*,err=10,end=20)at(:,3)
  call skipio(io_pwfn,4)
  read(io_pwfn,*,err=10,end=20)NGvec
  call skipio(io_pwfn,1)
  allocate(g_raw(3,NGvec),g(3,NGvec),g2(NGvec),stat=ialloc)
  if(ialloc/=0)call errstop('BWFN_HEADER','Allocation problem: G vectors.')
  do ig=1,NGvec
   read(io_pwfn,*,err=10,end=20)g_raw(:,ig)
  enddo ! ig
  call skipio(io_pwfn,4)
  read(io_pwfn,*,err=10,end=20)Nkvec

! Compute reciprocal lattice vectors (cols of bg), w/o factor of 2pi.
  call inve(at,bg)
  bg=transpose(bg)
  mag_rlv(1)=sqrt(sum(bg(1:3,1)**2))
  mag_rlv(2)=sqrt(sum(bg(1:3,2)**2))
  mag_rlv(3)=sqrt(sum(bg(1:3,3)**2))
  lvp(1)=bg(1,1)**2+bg(2,1)**2+bg(3,1)**2
  lvp(2)=bg(1,2)**2+bg(2,2)**2+bg(3,2)**2
  lvp(3)=bg(1,3)**2+bg(2,3)**2+bg(3,3)**2
  lvp(4)=2.d0*(bg(1,1)*bg(1,2)+bg(2,1)*bg(2,2)+bg(3,1)*bg(3,2))
  lvp(5)=2.d0*(bg(1,2)*bg(1,3)+bg(2,2)*bg(2,3)+bg(3,2)*bg(3,3))
  lvp(6)=2.d0*(bg(1,3)*bg(1,1)+bg(2,3)*bg(2,1)+bg(3,3)*bg(3,1))

  do ig=1,NGvec
   g2(ig)=sum(g_raw(:,ig)**2)
   g(:,ig)=nint(matmul(g_raw(:,ig),at)*rec_twopi)
  enddo ! ig
  gmax(1)=maxval(abs(g(1,:)))
  gmax(2)=maxval(abs(g(2,:)))
  gmax(3)=maxval(abs(g(3,:)))

! Work out blip grid.
  nr(:)=2*ceiling(dble(gmax(:))*xmul)+2
  rnr(:)=dble(nr(:))
  rnr2(:)=rnr(:)*rnr(:)

  write(6,*)'Blip grid: '//trim(i2s(nr(1)))//' by '//trim(i2s(nr(2))) &
   &//' by '//trim(i2s(nr(3)))
  write(6,*)

  allocate(kvec(nkvec,3),nband(nkvec,2),stat=ialloc)
  if(ialloc/=0)call errstop('BWFN_HEADER','Allocation error: kvec,nband.')

  if(spin_polarized)then
   no_spin_types=2
  else
   no_spin_types=1
  endif ! spin_polarized

! Calculate gamma.
  allocate(gamma(NGvec),stat=ialloc)
  if(ialloc/=0)call errstop('BWFN_HEADER','Allocation error: gamma.')
  da(:)=twopi/rnr(:)
  if(gamma_approx==1)then
   do ig=1,NGvec
    if(g(1,ig)/=0)then
     k=da(1)*g(1,ig) ; cosk=cos(k) ; k2=k*k ; k4=k2*k2
     gamma(ig)=k4/(6.d0*((cosk-2.d0)*cosk+1.d0))
    else
     gamma(ig)=two_thirds
    endif ! k/=0
    if(g(2,ig)/=0)then
     k=da(2)*g(2,ig) ; cosk=cos(k) ; k2=k*k ; k4=k2*k2
     gamma(ig)=gamma(ig)*k4/(6.d0*((cosk-2.d0)*cosk+1.d0))
    else
     gamma(ig)=gamma(ig)*two_thirds
    endif ! k/=0
    if(g(3,ig)/=0)then
     k=da(3)*g(3,ig) ; cosk=cos(k) ; k2=k*k ; k4=k2*k2
     gamma(ig)=gamma(ig)*k4/(6.d0*((cosk-2.d0)*cosk+1.d0))
    else
     gamma(ig)=gamma(ig)*two_thirds
    endif ! k/=0
   enddo ! ig
  elseif(gamma_approx==2)then
   do ig=1,NGvec
    gamma(ig)=1.d0/((1.d0+0.5d0*cos(da(1)*g(1,ig))) &
     &*(1.d0+0.5d0*cos(da(2)*g(2,ig)))*(1.d0+0.5d0*cos(da(3)*g(3,ig))))
   enddo ! ig
  else
   call errstop('BWFN_HEADER','Bug: bad gamma_approx.')
  endif ! gamma_approx

! Set up offset of atoms to put them in the middle of the blip grid in finite
! systems.
  R_offset=0.d0
  sum_atomno=sum(atomno)
  if(intended_periodicity<3.and.sum_atomno>0)then
! "Centre" of "molecule": average position of atoms weighted by the atomic
! number.
   molcentre=0.d0
   do na=1,nbasis
    molcentre=molcentre+dble(atomno(na))*atompos(:,na)
   enddo ! na
   molcentre=molcentre/dble(sum_atomno)
! Offset by (1/2)c-(R.c/|c|^2)c, where R is the centre of the molecule and c
! is the third lattice vector.
   R_offset=R_offset+(0.5d0-dot_product(molcentre,at(1:3,3)) &
     &/dot_product(at(1:3,3),at(1:3,3)))*at(1:3,3)
   if(intended_periodicity<2)then
! Offset by (1/2)b-(R.b/|b|^2)b, where R is the centre of the molecule and b
! is the second lattice vector.
    R_offset=R_offset+(0.5d0-dot_product(molcentre,at(1:3,2)) &
      &/dot_product(at(1:3,2),at(1:3,2)))*at(1:3,2)
    if(intended_periodicity==0)then
! Offset by (1/2)a-(R.a/|a|^2)a, where R is the centre of the molecule and a
! is the first lattice vector.
     R_offset=R_offset+(0.5d0-dot_product(molcentre,at(1:3,1)) &
       &/dot_product(at(1:3,1),at(1:3,1)))*at(1:3,1)
    endif ! 0D-periodic
   endif ! 0D or 1D-periodic

   write(6,*)'All atom positions etc. will be displaced by:'
   write(6,*)R_offset
   write(6,*)

! Translate the atoms.
   do na=1,nbasis
    atompos(1:3,na)=atompos(1:3,na)+R_offset
   enddo ! na

   allocate(expmigdotR_offset(NGvec),stat=ialloc)
   if(ialloc/=0)call errstop('BWFN_HEADER','Allocation error: &
    &expmigdotR_offset.')

   allocate(gfac1(-gmax(1):gmax(1)),gfac2(-gmax(2):gmax(2)), &
    &gfac3(-gmax(3):gmax(3)),stat=ialloc)
   if(ialloc/=0)call errstop('BWFN_HEADER','Allocation error: gfac.')

   gfac1(0)=cmplx(1.d0,0.d0,dp)
   if(gmax(1)>0)then
    lvdotR_offset=twopi*dot_product(bg(1:3,1),R_offset)
    gfac1(1)=cmplx(cos(lvdotR_offset),-sin(lvdotR_offset),dp)
    gfac1(-1)=conjg(gfac1(1))
    do i=2,gmax(1)
     gfac1(i)=gfac1(i-1)*gfac1(1)
     gfac1(-i)=conjg(gfac1(i))
    enddo ! i
   endif ! gmax(1)>0

   gfac2(0)=cmplx(1.d0,0.d0,dp)
   if(gmax(2)>0)then
    lvdotR_offset=twopi*dot_product(bg(1:3,2),R_offset)
    gfac2(1)=cmplx(cos(lvdotR_offset),-sin(lvdotR_offset),dp)
    gfac2(-1)=conjg(gfac2(1))
    do i=2,gmax(2)
     gfac2(i)=gfac2(i-1)*gfac2(1)
     gfac2(-i)=conjg(gfac2(i))
    enddo ! i
   endif ! gmax(2)>0

   gfac3(0)=cmplx(1.d0,0.d0,dp)
   if(gmax(3)>0)then
    lvdotR_offset=twopi*dot_product(bg(1:3,3),R_offset)
    gfac3(1)=cmplx(cos(lvdotR_offset),-sin(lvdotR_offset),dp)
    gfac3(-1)=conjg(gfac3(1))
    do i=2,gmax(3)
     gfac3(i)=gfac3(i-1)*gfac3(1)
     gfac3(-i)=conjg(gfac3(i))
    enddo ! i
   endif ! gmax(3)>0

   do ig=1,NGvec
    expmigdotR_offset(ig)=gfac1(g(1,ig))*gfac2(g(2,ig))*gfac3(g(3,ig))
   enddo ! ig

   deallocate(gfac1,gfac2,gfac3)

  endif ! intended_periodicity<3

! Write header of bwfn.data.
  write(io_bwfn,'(a)')trim(title)
  write(io_bwfn,*)
  write(io_bwfn,*)'BASIC INFO'
  write(io_bwfn,*)'----------'
  write(io_bwfn,*)'Generated by:'
  write(io_bwfn,*)trim(code)
  write(io_bwfn,*)'Method:'
  write(io_bwfn,*)trim(method)
  write(io_bwfn,*)'DFT Functional:'
  write(io_bwfn,*)trim(functional)
  write(io_bwfn,*)'Pseudopotential'
  write(io_bwfn,*)trim(pseudo_type)
  write(io_bwfn,*)'Plane wave cutoff (au)'
  write(io_bwfn,*)plane_wave_cutoff
  write(io_bwfn,*)'Spin polarized:'
  write(io_bwfn,*)spin_polarized
  write(io_bwfn,*)'Total energy (au per primitive cell)'
  write(io_bwfn,*)total_energy
  write(io_bwfn,*)'Kinetic energy (au per primitive cell)'
  write(io_bwfn,*)kinetic_energy
  write(io_bwfn,*)'Local potential energy (au per primitive cell)'
  write(io_bwfn,*)local_potential_energy
  write(io_bwfn,*)'Non local potential energy(au per primitive cell)'
  write(io_bwfn,*)non_local_potential_energy
  write(io_bwfn,*)'Electron electron energy (au per primitive cell)'
  write(io_bwfn,*)electron_electron_energy
  write(io_bwfn,*)'Ion-ion energy (au per primitive cell)'
  write(io_bwfn,*)eionion
  write(io_bwfn,*)'Number of electrons per primitive cell'
  write(io_bwfn,*)nelec
  write(io_bwfn,*)
  write(io_bwfn,*)'GEOMETRY'
  write(io_bwfn,*)'-------- '
  write(io_bwfn,*)'Number of atoms per primitive cell '
  write(io_bwfn,*)nbasis
  write(io_bwfn,*)'Atomic number and position of the atoms(au) '
  do na=1,nbasis
   write(io_bwfn,'(1x,a,3(1x,es24.17))')trim(i2s(atomno(na))),atompos(:,na)
  enddo ! na
  write(io_bwfn,*)'Primitive lattice vectors (au)'
  write(io_bwfn,'(3(1x,es24.17))')at(:,1)
  write(io_bwfn,'(3(1x,es24.17))')at(:,2)
  write(io_bwfn,'(3(1x,es24.17))')at(:,3)
  write(io_bwfn,*)
  write(io_bwfn,*)'G VECTORS'
  write(io_bwfn,*)'---------'
  write(io_bwfn,*)'Number of G-vectors'
  write(io_bwfn,*)NGvec
  write(io_bwfn,*)'Gx Gy Gz (au)'

  do ig=1,NGvec
   write(io_bwfn,'(3(1x,es24.17))')g_raw(:,ig)
  enddo ! ig

  write(io_bwfn,*)'Blip grid'
  write(io_bwfn,*)trim(i2s(nr(1)))//" "//trim(i2s(nr(2)))//" " &
   &//trim(i2s(nr(3)))

! Read in number of k points.
  write(io_bwfn,*)
  write(io_bwfn,*)'WAVE FUNCTION'
  write(io_bwfn,*)'-------------'
  write(io_bwfn,*)'Number of k-points'
  write(io_bwfn,*)Nkvec

  return
10 call errstop('BWFN_HEADER','Error reading pwfn.data.')
20 call errstop('BWFN_HEADER','Unexpectedly reached end of pwfn.data.')
 END SUBROUTINE bwfn_header


 SUBROUTINE bwfn_transform
!-----------------------------------------------------------------------!
! Carry out PW->blip transformation.                                    !
! The header line of the very first k-point is read here, because it is !
! needed for the decision between blipp and blipk.                      !
!-----------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER ik

  call skipio(io_pwfn,1)
  write(io_bwfn,*)'k-point # ; # of bands (up spin/down spin) ; &
   &k-point coords (au)'
  read(io_pwfn,*,err=10,end=20)ik,nband(1,1:2),kvec(1,:)
  write(io_bwfn,'(1x,a,3(1x,es24.17))')trim(i2s(ik))//" "//&
   &trim(i2s(nband(1,1)))//" "//trim(i2s(nband(1,2))),kvec(1,:)

  if(nkvec==1.and.all(abs(kvec(1,:))<eps))then
   call blipp
  else
   call blipk
  endif ! Just Gamma present
  if(intended_periodicity<3)deallocate(expmigdotR_offset)
  write(6,*)

  ! Check that k vectors are appropriate for the intended periodicity.
  if(intended_periodicity<3)then
   if(any(kvec(:,3)/=0.d0))then
    call wordwrap('WARNING: the k vectors should be zero in the directions of &
     &reduced periodicity (z).')
    write(6,*)
   endif ! Multiple k in "z" direction.
   if(intended_periodicity<2)then
    if(any(kvec(:,2)/=0.d0))then
     call wordwrap('WARNING: the k vectors should be zero in the directions of &
      &reduced periodicity (y).')
     write(6,*)
    endif ! Multiple k in "y" direction.
    if(intended_periodicity==0)then
     if(any(kvec(:,1)/=0.d0))then
      call wordwrap('WARNING: the k vectors should be zero in the directions &
       &of reduced periodicity (x).')
      write(6,*)
     endif ! Multiple k in "x" direction.
    endif ! 0D-periodic
   endif ! 0D or 1D-periodic
  endif ! Reduced periodicity

  return
10 call errstop('BWFN_TRANSFORM','Error reading pwfn.data.')
20 call errstop('BWFN_TRANSFORM','Unexpectedly reached end of pwfn.data.')
 END SUBROUTINE bwfn_transform


! Subroutine blipp is in a separate file.
 INCLUDE "blipp.f90"


! Subroutine blipk is in a separate file.
 INCLUDE "blipk.f90"


 SUBROUTINE write_overlaps(ibnd,av_overlap,av_overlap2)
!-------------------------------------------------------------------------!
! Write out the overlaps of the value, gradient and Laplacian of the blip !
! orbitals.  Give error bars where possible.                              !
!-------------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: ibnd
  REAL(dp),INTENT(in) :: av_overlap(5),av_overlap2(5)
  INTEGER k
  REAL(dp) err_overlap(5)
  CHARACTER(4) char4
  CHARACTER(12) char12_arr(5)
  if(n_overlap_tests<2)call errstop('WRITE_OVERLAPS','Error: need at least &
   &two overlap tests, to estimate error bars.')
  err_overlap=sqrt(max(av_overlap2-av_overlap**2,0.d0)/dble(n_overlap_tests-1))
  char4=trim(i2s(ibnd)) ; char4=adjustr(char4)
  do k=1,5
   char12_arr(k)=trim(write_mean(av_overlap(k),err_overlap(k)))
! Not room to quote error bar.  Just quote mean.
   if(index(char12_arr(k),')')==0)write(char12_arr(k),'(f12.9)')av_overlap(k)
  enddo ! k
  write(6,'(1x,a,2(1x,a),2x,3(1x,a))')char4,char12_arr(1:5)
 END SUBROUTINE write_overlaps


END MODULE blip_utils


PROGRAM blip
!---------------------------!
! Main program starts here. !
!---------------------------!
 USE blip_utils,ONLY : get_args,initialize,bwfn_header,bwfn_transform,finalize
 IMPLICIT NONE

 write(6,*)
 write(6,*)'PLANE WAVE -> BLIP CONVERSION UTILITY'
 write(6,*)'====================================='
 write(6,*)

! Ask user for xmul etc.
 call get_args
! Open input and output file.
 call initialize

! Transfer information in header of pwfn.data into bwfn.data and read in
! geometry and G vectors.
 call bwfn_header

! Carry out PW->blip transformation.
 call bwfn_transform

! Close input and output file.
 call finalize

 write(6,*)'Program finished.  bwfn.data has been generated.'
 write(6,*)

END PROGRAM blip


! The following functions need to be external, because they are used as
! arguments to adapt_simpson

REAL(kind=kind(1.d0)) FUNCTION norm2_integrand_phi(phi)
!-------------------------------------------------------------------------!
! This function returns psi^2.r^2.sin(theta) for a given point relative   !
! to the centre of psi specified by spherical polar coordinates r, theta  !
! and phi.                                                                !
!-------------------------------------------------------------------------!
 USE blip_utils, ONLY : bg,r_sin_theta_pass,r_cos_theta_pass, &
  &r2_sin_theta_pass,centrepos_lv_pass,blip3d
 USE helpers, ONLY : dp
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: phi
 REAL(dp) rvec(3),r1(3),rpsi,lap,grad(3)
 rvec=(/r_sin_theta_pass*cos(phi),r_sin_theta_pass*sin(phi), &
  &r_cos_theta_pass/)+centrepos_lv_pass
 r1=matmul(rvec,bg)+centrepos_lv_pass
 call blip3d(rpsi,lap,grad,r1,.true.,.false.)
 norm2_integrand_phi=rpsi**2*r2_sin_theta_pass
END FUNCTION norm2_integrand_phi


REAL(kind=kind(1.d0)) FUNCTION norm2_integrand_theta(theta)
!--------------------------------------------------------------------------!
! This function returns the integral of psi^2.r^2.sin(theta) with respect  !
! to phi for a given theta and r, where r is the distance from the centre. !
!--------------------------------------------------------------------------!
 USE blip_utils, ONLY : r_sin_theta_pass,r_cos_theta_pass, &
  &r2_sin_theta_pass,r_pass
 USE helpers, ONLY : adapt_simpson,twopi,dp
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: theta
 INTERFACE
  REAL(kind=kind(1.d0)) FUNCTION norm2_integrand_phi(phi)
   USE helpers,ONLY : dp
   REAL(dp),INTENT(in) :: phi
  END FUNCTION norm2_integrand_phi
 END INTERFACE
 r_sin_theta_pass=r_pass*sin(theta)
 r2_sin_theta_pass=r_pass*r_sin_theta_pass
 r_cos_theta_pass=r_pass*cos(theta)
 call adapt_simpson(norm2_integrand_phi,0.d0,twopi,norm2_integrand_theta)
END FUNCTION norm2_integrand_theta


REAL(kind=kind(1.d0)) FUNCTION norm2_integrand_r(r)
!--------------------------------------------------------------------------!
! This function returns the integral of psi^2.r^2.sin(theta) with respect  !
! to theta and phi for a given r, where r is the distance from the centre. !
!--------------------------------------------------------------------------!
 USE blip_utils, ONLY : r_pass
 USE helpers, ONLY : adapt_simpson,pi,dp
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: r
 INTERFACE
  REAL(kind=kind(1.d0)) FUNCTION norm2_integrand_theta(theta)
   USE helpers,ONLY : dp
   REAL(dp),INTENT(in) :: theta
  END FUNCTION norm2_integrand_theta
 END INTERFACE
 r_pass=r
 call adapt_simpson(norm2_integrand_theta,0.d0,pi,norm2_integrand_r)
END FUNCTION norm2_integrand_r


REAL(kind=kind(1.d0)) FUNCTION norm2_integrand_xlv(xlv)
!--------------------------------------------------------------------------!
! This function returns psi^2 at a point whose coordinates in terms of the !
! lattice vectors are xlv, ylv and zlv.                                    !
!--------------------------------------------------------------------------!
 USE blip_utils, ONLY : blip3d,rvec_lv_pass
 USE helpers,ONLY : dp
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: xlv
 REAL(dp) rpsi,lap,grad(3)
 rvec_lv_pass(1)=xlv
 call blip3d(rpsi,lap,grad,rvec_lv_pass,.true.,.false.)
 norm2_integrand_xlv=rpsi**2
END FUNCTION norm2_integrand_xlv


REAL(kind=kind(1.d0)) FUNCTION norm2_integrand_ylv(ylv)
!---------------------------------------------------------------------------!
! This function returns the integral of psi^2 with respect to xlv for a     !
! given ylv and zlv, where xlv, ylv and zlv are the coordinates of a point  !
! in terms of the lattice vectors.                                          !
!---------------------------------------------------------------------------!
 USE blip_utils, ONLY : rvec_lv_pass,rlv_min_pass,rlv_max_pass
 USE helpers, ONLY : adapt_simpson,dp
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: ylv
 INTERFACE
  REAL(kind=kind(1.d0)) FUNCTION norm2_integrand_xlv(xlv)
   USE helpers,ONLY : dp
   REAL(dp),INTENT(in) :: xlv
  END FUNCTION norm2_integrand_xlv
 END INTERFACE
 rvec_lv_pass(2)=ylv
 call adapt_simpson(norm2_integrand_xlv,rlv_min_pass(1),rlv_max_pass(1), &
  &norm2_integrand_ylv)
END FUNCTION norm2_integrand_ylv


REAL(kind=kind(1.d0)) FUNCTION norm2_integrand_zlv(zlv)
!-------------------------------------------------------------------------!
! This function returns the integral of psi^2 with respect to xlv and ylv !
! for a given zlv, where xlv, ylv and zlv are the coordinates of a point  !
! in terms of the lattice vectors.                                        !
!-------------------------------------------------------------------------!
 USE blip_utils, ONLY : rvec_lv_pass,rlv_min_pass,rlv_max_pass
 USE helpers, ONLY : adapt_simpson,dp
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: zlv
 INTERFACE
  REAL(kind=kind(1.d0)) FUNCTION norm2_integrand_ylv(ylv)
   USE helpers,ONLY : dp
   REAL(dp),INTENT(in) :: ylv
  END FUNCTION norm2_integrand_ylv
 END INTERFACE
 rvec_lv_pass(3)=zlv
 call adapt_simpson(norm2_integrand_ylv,rlv_min_pass(2),rlv_max_pass(2), &
  &norm2_integrand_zlv)
END FUNCTION norm2_integrand_zlv


REAL(kind=kind(1.d0)) FUNCTION c_norm2_integrand_xlv(xlv)
!----------------------------------------------------------------------------!
! This function returns |psi|^2 at a point whose coordinates in terms of the !
! lattice vectors are xlv, ylv and zlv.                                      !
!----------------------------------------------------------------------------!
 USE blip_utils, ONLY : blip3dk,rvec_lv_pass
 USE helpers, ONLY : dp
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: xlv
 COMPLEX(dp) rpsi,lap,grad(3)
 rvec_lv_pass(1)=xlv
 call blip3dk(rpsi,lap,grad,rvec_lv_pass,.true.,.false.)
 c_norm2_integrand_xlv=dble(conjg(rpsi)*rpsi)
END FUNCTION c_norm2_integrand_xlv


REAL(kind=kind(1.d0)) FUNCTION c_norm2_integrand_ylv(ylv)
!---------------------------------------------------------------------------!
! This function returns the integral of |psi|^2 with respect to xlv for a   !
! given ylv and zlv, where xlv, ylv and zlv are the coordinates of a point  !
! in terms of the lattice vectors.                                          !
!---------------------------------------------------------------------------!
 USE blip_utils, ONLY : rvec_lv_pass,rlv_min_pass,rlv_max_pass
 USE helpers, ONLY : adapt_simpson,dp
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: ylv
 INTERFACE
  REAL(kind=kind(1.d0)) FUNCTION c_norm2_integrand_xlv(xlv)
   USE helpers,ONLY : dp
   REAL(dp),INTENT(in) :: xlv
  END FUNCTION c_norm2_integrand_xlv
 END INTERFACE
 rvec_lv_pass(2)=ylv
 call adapt_simpson(c_norm2_integrand_xlv,rlv_min_pass(1),rlv_max_pass(1), &
  &c_norm2_integrand_ylv)
END FUNCTION c_norm2_integrand_ylv


REAL(kind=kind(1.d0)) FUNCTION c_norm2_integrand_zlv(zlv)
!---------------------------------------------------------------------------!
! This function returns the integral of |psi|^2 with respect to xlv and ylv !
! for a given zlv, where xlv, ylv and zlv are the coordinates of a point    !
! in terms of the lattice vectors.                                          !
!---------------------------------------------------------------------------!
 USE blip_utils, ONLY : rvec_lv_pass,rlv_min_pass,rlv_max_pass
 USE helpers, ONLY : adapt_simpson,dp
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: zlv
 INTERFACE
  REAL(kind=kind(1.d0)) FUNCTION c_norm2_integrand_ylv(ylv)
   USE helpers,ONLY : dp
   REAL(dp),INTENT(in) :: ylv
  END FUNCTION c_norm2_integrand_ylv
 END INTERFACE
 rvec_lv_pass(3)=zlv
 call adapt_simpson(c_norm2_integrand_ylv,rlv_min_pass(2),rlv_max_pass(2), &
  &c_norm2_integrand_zlv)
END FUNCTION c_norm2_integrand_zlv


REAL(kind=kind(1.d0)) FUNCTION ke_integrand_phi(phi)
!------------------------------------------------------------------------!
! This function returns psi.lap(psi).r^2.sin(theta) for a given point    !
! relative to the centre of psi specified by spherical polar coordinates !
! r, theta and phi.                                                      !
! (One can instead return -|grad(psi)|^2.r^2.sin(theta) - should give    !
! the same KE when integrated.  Uncomment the indicated lines.)          !
!------------------------------------------------------------------------!
 USE blip_utils, ONLY : bg,r_sin_theta_pass,r_cos_theta_pass, &
  &r2_sin_theta_pass,centrepos_lv_pass,blip3d
 USE helpers,ONLY : dp
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: phi
 REAL(dp) rvec(3),r1(3),rpsi,lap,grad(3)
 rvec=(/r_sin_theta_pass*cos(phi),r_sin_theta_pass*sin(phi), &
  &r_cos_theta_pass/)+centrepos_lv_pass
 r1=matmul(rvec,bg)+centrepos_lv_pass
 call blip3d(rpsi,lap,grad,r1,.true.,.true.)                 ! psi.lap(psi)
 ke_integrand_phi=rpsi*lap*r2_sin_theta_pass                 ! psi.lap(psi)
! call blip3d(rpsi,lap,grad,r1,.false.,.true.)               ! -|grad(psi)|^2
! ke_integrand_phi=-dot_product(grad,grad)*r2_sin_theta_pass ! -|grad(psi)|^2
END FUNCTION ke_integrand_phi


REAL(kind=kind(1.d0)) FUNCTION ke_integrand_theta(theta)
!------------------------------------------------------------------------!
! This function returns the integral of psi.lap(psi).r^2.sin(theta) with !
! respect to phi for a given theta and r, where r is the distance from   !
! the centre.                                                            !
!------------------------------------------------------------------------!
 USE blip_utils, ONLY : r_sin_theta_pass,r_cos_theta_pass, &
  &r2_sin_theta_pass,r_pass
 USE helpers, ONLY : adapt_simpson,twopi,dp
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: theta
 INTERFACE
  REAL(kind=kind(1.d0)) FUNCTION ke_integrand_phi(phi)
   USE helpers,ONLY : dp
   REAL(dp),INTENT(in) :: phi
  END FUNCTION ke_integrand_phi
 END INTERFACE
 r_sin_theta_pass=r_pass*sin(theta)
 r2_sin_theta_pass=r_pass*r_sin_theta_pass
 r_cos_theta_pass=r_pass*cos(theta)
 call adapt_simpson(ke_integrand_phi,0.d0,twopi,ke_integrand_theta)
END FUNCTION ke_integrand_theta


REAL(kind=kind(1.d0)) FUNCTION ke_integrand_r(r)
!------------------------------------------------------------------------!
! This function returns the integral of psi.lap(psi).r^2.sin(theta) with !
! respect to theta and phi for a given r, where r is the distance from   !
! the centre.                                                            !
!------------------------------------------------------------------------!
 USE blip_utils, ONLY : r_pass
 USE helpers, ONLY : adapt_simpson,pi,dp
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: r
 INTERFACE
  REAL(kind=kind(1.d0)) FUNCTION ke_integrand_theta(theta)
   USE helpers,ONLY : dp
   REAL(dp),INTENT(in) :: theta
  END FUNCTION ke_integrand_theta
 END INTERFACE
 r_pass=r
 call adapt_simpson(ke_integrand_theta,0.d0,pi,ke_integrand_r)
END FUNCTION ke_integrand_r


REAL(kind=kind(1.d0)) FUNCTION ke_integrand_xlv(xlv)
!--------------------------------------------------------------------------!
! This function returns psi.lap(psi) at a point whose coordinates in terms !
! of the lattice vectors are xlv, ylv and zlv.                             !
! (One can instead return -|grad(psi)|^2 - should give the same KE when    !
! integrated.  Uncomment the indicated lines.)                             !
!--------------------------------------------------------------------------!
 USE blip_utils, ONLY : blip3d,rvec_lv_pass
 USE helpers, ONLY : dp
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: xlv
 REAL(dp) rpsi,lap,grad(3)
 rvec_lv_pass(1)=xlv
 call blip3d(rpsi,lap,grad,rvec_lv_pass,.true.,.true.)     ! psi.lap(psi)
 ke_integrand_xlv=rpsi*lap                                 ! psi.lap(psi)
! call blip3d(rpsi,lap,grad,rvec_lv_pass,.false.,.true.)   ! -|grad(psi)|^2
! ke_integrand_xlv=-dot_product(grad,grad)                 ! -|grad(psi)|^2
END FUNCTION ke_integrand_xlv


REAL(kind=kind(1.d0)) FUNCTION ke_integrand_ylv(ylv)
!----------------------------------------------------------------------------!
! This function returns the integral of psi.lap(psi) with respect to xlv for !
! a given ylv and zlv, where xlv, ylv and zlv are the coordinates of a point !
! in terms of the lattice vectors.                                           !
!----------------------------------------------------------------------------!
 USE blip_utils, ONLY : rvec_lv_pass,rlv_min_pass,rlv_max_pass
 USE helpers, ONLY : adapt_simpson,dp
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: ylv
 INTERFACE
  REAL(kind=kind(1.d0)) FUNCTION ke_integrand_xlv(xlv)
   USE helpers,ONLY : dp
   REAL(dp),INTENT(in) :: xlv
  END FUNCTION ke_integrand_xlv
 END INTERFACE
 rvec_lv_pass(2)=ylv
 call adapt_simpson(ke_integrand_xlv,rlv_min_pass(1),rlv_max_pass(1),ke_integrand_ylv)
END FUNCTION ke_integrand_ylv


REAL(kind=kind(1.d0)) FUNCTION ke_integrand_zlv(zlv)
!----------------------------------------------------------------------------!
! This function returns the integral of psi.lap(psi) with respect to xlv and !
! ylv for a given zlv, where xlv, ylv and zlv are the coordinates of a point !
! in terms of the lattice vectors.                                           !
!----------------------------------------------------------------------------!
 USE blip_utils, ONLY : rvec_lv_pass,rlv_min_pass,rlv_max_pass
 USE helpers, ONLY : adapt_simpson,dp
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: zlv
 INTERFACE
  REAL(kind=kind(1.d0)) FUNCTION ke_integrand_ylv(ylv)
   USE helpers,ONLY : dp
   REAL(dp),INTENT(in) :: ylv
  END FUNCTION ke_integrand_ylv
 END INTERFACE
 rvec_lv_pass(3)=zlv
 call adapt_simpson(ke_integrand_ylv,rlv_min_pass(2),rlv_max_pass(2),ke_integrand_zlv)
END FUNCTION ke_integrand_zlv


REAL(kind=kind(1.d0)) FUNCTION c_ke_integrand_xlv(xlv)
!-------------------------------------------------------------------------!
! This function returns Re(psi*.lap(psi)) at a point whose coordinates in !
! terms of the lattice vectors are xlv, ylv and zlv.                      !
! (One can instead return -|grad(psi)|^2 - should give the same KE when   !
! integrated.  Uncomment the indicated lines.)                            !
!-------------------------------------------------------------------------!
 USE blip_utils, ONLY : blip3dk,rvec_lv_pass,kvec2_pass,kvec_pass
 USE helpers, ONLY : dp,twoi !,iunity
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: xlv
 COMPLEX(dp) rpsi,lap,grad(3)
 rvec_lv_pass(1)=xlv
 call blip3dk(rpsi,lap,grad,rvec_lv_pass,.true.,.true.)
 c_ke_integrand_xlv=dble(conjg(rpsi)*(-kvec2_pass*rpsi &     ! psi*.lap(psi)
  &+twoi*dot_product(kvec_pass,grad)+lap))                   ! psi*.lap(psi)
! c_ke_integrand_xlv=-dble(kvec2_pass*conjg(rpsi)*rpsi &     ! -|grad(psi)|^2
!  &-iunity*conjg(rpsi)*dot_product(kvec_pass,grad) &        ! -|grad(psi)|^2
!  &+iunity*rpsi*conjg(dot_product(kvec_pass,grad)) &        ! -|grad(psi)|^2
!  &+dot_product(grad,grad))                                 ! -|grad(psi)|^2
!               (NB dot_product(a,b)=a*.b.)
END FUNCTION c_ke_integrand_xlv


REAL(kind=kind(1.d0)) FUNCTION c_ke_integrand_ylv(ylv)
!-----------------------------------------------------------------------------!
! This function returns the integral of Re(psi*.lap(psi)) with respect to xlv !
! for a given ylv and zlv, where xlv, ylv and zlv are the coordinates of a    !
! point in terms of the lattice vectors.                                      !
!-----------------------------------------------------------------------------!
 USE blip_utils, ONLY : rvec_lv_pass,rlv_min_pass,rlv_max_pass
 USE helpers, ONLY : adapt_simpson,dp
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: ylv
 INTERFACE
  REAL(kind=kind(1.d0)) FUNCTION c_ke_integrand_xlv(xlv)
   USE helpers,ONLY : dp
   REAL(dp),INTENT(in) :: xlv
  END FUNCTION c_ke_integrand_xlv
 END INTERFACE
 rvec_lv_pass(2)=ylv
 call adapt_simpson(c_ke_integrand_xlv,rlv_min_pass(1),rlv_max_pass(1), &
  &c_ke_integrand_ylv)
END FUNCTION c_ke_integrand_ylv


REAL(kind=kind(1.d0)) FUNCTION c_ke_integrand_zlv(zlv)
!----------------------------------------------------------------------------!
! This function returns the integral of Re(psi*.lap(psi)) with respect to    !
! xlv and ylv for a given zlv, where xlv, ylv and zlv are the coordinates of !
! a point in terms of the lattice vectors.                                   !
!----------------------------------------------------------------------------!
 USE blip_utils, ONLY : rvec_lv_pass,rlv_min_pass,rlv_max_pass
 USE helpers, ONLY : adapt_simpson,dp
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: zlv
 INTERFACE
  REAL(kind=kind(1.d0)) FUNCTION c_ke_integrand_ylv(ylv)
   USE helpers,ONLY : dp
   REAL(dp),INTENT(in) :: ylv
  END FUNCTION c_ke_integrand_ylv
 END INTERFACE
 rvec_lv_pass(3)=zlv
 call adapt_simpson(c_ke_integrand_ylv,rlv_min_pass(2),rlv_max_pass(2), &
  &c_ke_integrand_zlv)
END FUNCTION c_ke_integrand_zlv

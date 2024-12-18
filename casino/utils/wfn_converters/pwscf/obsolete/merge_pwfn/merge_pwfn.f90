!-----------------------------------------------------------------------------!
! This program merges pwfn.datax [1=1,2,...] files produced by a parallel     !
! PWSCF run into one pwfn.data file.                                          !
!                                                                             !
! Usage:                                                                      !
! ./pwfn_merge                                                                !
!                                                                             !
! D. Alfe                                                                     !
!                                                                             !
! Changes                                                                     !
! -------                                                                     !
! 1.2006 MDT - CASINO formatted.                                              !
! 2.2008 MDT - Removed unused variable declarations.                          !
!-----------------------------------------------------------------------------!

MODULE utils
 IMPLICIT NONE


 CONTAINS


 SUBROUTINE skipio(io,n)
!---------------------------------------!
! Skip n lines from file on channel io. !
!---------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: io,n
 INTEGER i
 do i=1,n
  read(io,*)
 enddo
 END SUBROUTINE skipio


 SUBROUTINE create_index(y,x_index)
!-----------------------------------------------------------------------------!
! This subroutine creates an index array x_index for the n items of data in   !
! the array y.  Adapted from Numerical Recipes.                               !
!-----------------------------------------------------------------------------!
 IMPLICIT NONE
 DOUBLE PRECISION,INTENT(in) :: y(:)
 INTEGER,INTENT(out) :: x_index(:)
 INTEGER,PARAMETER :: ins_sort_thresh=7,stacksize=80
 INTEGER n,i,x_indexj,ir,itemp,j,jstack,k,l,lp1,istack(stacksize)
 DOUBLE PRECISION yj
 n=size(x_index)
 do j=1,n
  x_index(j)=j
 enddo ! j
 if(n<=1)return
 jstack=0
 l=1
 ir=n
 do
  if(ir-l<ins_sort_thresh)then
jloop : do j=l+1,ir
    x_indexj=x_index(j) ; yj=y(x_indexj)
    do i=j-1,l,-1
     if(y(x_index(i))<=yj)then
      x_index(i+1)=x_indexj
      cycle jloop
     endif ! y(x_index(i))<=yj
     x_index(i+1)=x_index(i)
    enddo ! i
    x_index(l)=x_indexj
   enddo jloop ! j
   if(jstack==0)return
   ir=istack(jstack)
   l=istack(jstack-1)
   jstack=jstack-2
  else
   k=(l+ir)/2
   lp1=l+1
   itemp=x_index(k)    ; x_index(k)=x_index(lp1)  ; x_index(lp1)=itemp
   if(y(x_index(l))>y(x_index(ir)))then
    itemp=x_index(l)   ; x_index(l)=x_index(ir)   ; x_index(ir)=itemp
   endif
   if(y(x_index(lp1))>y(x_index(ir)))then
    itemp=x_index(lp1) ; x_index(lp1)=x_index(ir) ; x_index(ir)=itemp
   endif
   if(y(x_index(l))>y(x_index(lp1)))then
    itemp=x_index(l)   ; x_index(l)=x_index(lp1)  ; x_index(lp1)=itemp
   endif
   i=lp1
   j=ir
   x_indexj=x_index(lp1)
   yj=y(x_indexj)
   do
    do
     i=i+1
     if(y(x_index(i))>=yj)exit
    enddo ! i
    do
     j=j-1
     if(y(x_index(j))<=yj)exit
    enddo ! j
    if(j<i)exit
    itemp=x_index(i) ; x_index(i)=x_index(j) ; x_index(j)=itemp
   enddo
   x_index(lp1)=x_index(j)
   x_index(j)=x_indexj
   jstack=jstack+2
   if(jstack>stacksize)then
    write(6,*)'stacksize is too small.'
    stop
   endif ! jstack>stacksize
   if(ir-i+1>=j-l)then
    istack(jstack)=ir
    istack(jstack-1)=i
    ir=j-1
   else
    istack(jstack)=j-1
    istack(jstack-1)=l
    l=i
   endif ! ir-i+1>=j-l
  endif ! ir-l<ins_sort_thresh
 enddo
 END SUBROUTINE create_index


END MODULE utils


PROGRAM pwfn_merge
!---------------------------!
! Main program starts here. !
!---------------------------!
  USE utils
  IMPLICIT NONE
  INTEGER i,ibnd,ig,npwx,nks,ik,nproc,nat,na,atomno, &
   &nbnd(2),ispin,nspin,nelec
  INTEGER :: io=91
  INTEGER,PARAMETER :: dp=kind(1.d0)
  INTEGER, ALLOCATABLE :: npw(:),indx(:)
  REAL(dp) at(3,3),g(3),tau(3),etot,ewld,ecut,xk(3),et
  REAL(dp),ALLOCATABLE :: g1(:,:),g2(:)
  COMPLEX(dp) evc
  COMPLEX(dp),ALLOCATABLE :: evc1(:)
  LOGICAL lsda
  CHARACTER(3) nd_nmbr
  CHARACTER(80) title

  print*,'MERGE_PWFN'
  print*,'=========='
  print*,'How many files?'
  read(*,*) nproc
  allocate(npw(nproc))

  if(nproc<1)then
   write(*,*)'Must have at least one file.'
   stop
  elseif(nproc<10)then
   nd_nmbr='1'
  elseif(nproc<100)then
   nd_nmbr='01'
  elseif(nproc<1000)then
   nd_nmbr='001'
  else
   write(*,*)'Need to update merge_pwfn for more than 1000 processors.'
   stop
  endif ! nproc

  open(11,file='pwfn.data'//nd_nmbr,status='old')
  open(io,file='tmp.data',status='unknown')

  read(11,'(a)')title
  write(io,'(a)')title
  call skipio(11,11)
  write(io,'(a)')
  write(io,'(a)')'BASIC INFO'
  write(io,'(a)')'----------'
  write(io,'(a)')'Generated by:'
  write(io,'(a)')' PWSCF'
  write(io,'(a)')'method:'
  write(io,'(a)')' DFT'
  write(io,'(a)')'DFT Functional:'
  write(io,'(a)')' unknown'
  write(io,'(a)')'Pseudopotential'
  write(io,'(a)')' unknown'
  call skipio(11,1)
  write(io,'(a)')'Plane wave cutoff (au)'
  read(11,*)ecut
  write(io,*)ecut
  call skipio(11,1)
  write(io,'(a)')'Spin polarized:'
  read(11,*)lsda
  write(io,*)lsda
  if(lsda)then
   nspin=2
  else
   nspin=1
  endif
  call skipio(11,1)
  write(io,'(a)')'Total energy (au per primitive cell)'
  read(11,*)etot
  write(io,*)etot
  call skipio(11,1)
  write(io,'(a)')'Kinetic energy (au per primitive cell)'
  read(11,*)etot
  write(io,*)etot
  call skipio(11,1)
  write(io,'(a)')'Local potential energy (au per primitive cell)'
  read(11,*)etot
  write(io,*)etot
  call skipio(11,1)
  write(io,'(a)')'Non local potential energy(au per primitive cell)'
  read(11,*)etot
  write(io,*)etot
  call skipio(11,1)
  write(io,'(a)')'Electron-electron energy (au per primitive cell)'
  read(11,*)etot
  write(io,*)etot
  call skipio(11,1)
  write(io,'(a)')'Ion-ion energy (au per primitive cell)'
  read(11,*)ewld
  write(io,*)ewld
  call skipio(11,1)
  write(io,'(a)')'Number of electrons per primitive cell'
  read(11,*)nelec
  write(io,*)nelec
  call skipio(11,4)
  write(io,'(a)')' '
  write(io,'(a)')'GEOMETRY'
  write(io,'(a)')'-------- '
  write(io,'(a)')'Number of atoms per primitive cell '
  read(11,*)nat
  write(io,*)nat
  call skipio(11,1)
  write(io,'(a)')'Atomic number and position of the atoms(au) '
  do na=1,nat
   read(11,*)atomno,tau(:)
   write(io,'(i6,3f20.12)')atomno,tau(:)
  enddo
  call skipio(11,1)
  write(io,'(a)')'primitive lattice vectors (au) '
  read(11,*)at(:,1)
  read(11,*)at(:,2)
  read(11,*)at(:,3)
  write(io,'(3f20.12)')at(:,1)
  write(io,'(3f20.12)')at(:,2)
  write(io,'(3f20.12)')at(:,3)

  call skipio(11,4)
  write(io,'(a)')' '
  write(io,'(a)')'G VECTORS'
  write(io,'(a)')'---------'
  write(io,'(a)')'Number of g-vectors'


  read(11,*)npw(1)
  npwx=npw(1)
  print*,'File ',nd_nmbr,' contains ',npwx,' plane waves.'
  nd_nmbr='   '
  do i=2,nproc
   if(nproc<10)then
    write(nd_nmbr(1:1),'(i1)')i
   elseif(nproc<100)then
    if(i<10)then
     nd_nmbr='0'
     write(nd_nmbr(2:2),'(i1)')i
    else
     write(nd_nmbr(1:2),'(i2)')i
    endif
   else
    if(i<10)then
     nd_nmbr='00'
     write(nd_nmbr(3:3),'(i1)')i
    elseif(i<100)then
     nd_nmbr='0'
     write(nd_nmbr(2:3),'(i2)')i
    else
     write(nd_nmbr,'(i3)')i
    endif
   endif
   open(i+10,file='pwfn.data'//nd_nmbr,status='old')
   call skipio(i+10,36+nat+8)
   read(i+10,*)npw(i)
   npwx=npwx+npw(i)
   print*,'File ',nd_nmbr,' contains ',npw(i),' plane waves, &
    &total up to now: ',npwx
  enddo

  write(io,*)npwx
  write(io,'(a)')'gx gy gz (au)'

  do i=1,nproc
   call skipio(10+i,1)
   do ig=1,npw(i)
    read(i+10,*)g(:)
    write(io,'(3f20.12)')g(:)
   enddo
  enddo

  write(io,'(a)')' '
  write(io,'(a)')'WAVE FUNCTION'
  write(io,'(a)')'-------------'
  write(io,'(a)')'Number of k-points'
  call skipio(11,4)
  read(11,*)nks
  write(io,*)nks

  do i=2,nproc
   call skipio(i+10,5)
  enddo
  do ik=1,nks
   do i=2,nproc
    call skipio(i+10,2)
   enddo
   call skipio(11,1)
   write(io,'(a)')'k-point # ; # of bands (up spin/down spin) ; &
    &k-point coords (au)'
   read(11,*)i,nbnd(1),nbnd(2),xk(:)
    write(io,'(3i4,3f20.16)') ik, nbnd(1), nbnd(2), xk(:)
   do ispin=1,nspin
    do ibnd=1,nbnd(ispin)
     call skipio(11,1)
     write(io,'(a)')'band, spin, eigenvalue (au)'
     read(11,*)i,i,et
     write(io,'(2i5,f20.12)')ibnd,ispin,et
     call skipio(11,1)
! Read in the wavefunctions
     write(io,'(a)') 'eigenvectors coefficients'
     do ig=1, npw(1)
      read(11,*) evc
      write(io,*) evc
     enddo
     do i=2,nproc
      call skipio(i+10,3)
      do ig=1,npw(i)
       read(i+10,*)evc
       write(io,*)evc
      enddo
     enddo
    enddo
   enddo
  enddo

  close(io)
  close(11)

  print*,'Reordering G vectors...'

  open(11,file='tmp.data',status='old')
  open(io,file='pwfn.data',status='unknown')

  read(11,'(a)') title
  write(io,'(a)') title
  call skipio(11,11)
  write(io,'(a)')
  write(io,'(a)')'BASIC INFO'
  write(io,'(a)')'----------'
  write(io,'(a)')'Generated by:'
  write(io,'(a)')' PWSCF'
  write(io,'(a)')'Method:'
  write(io,'(a)')' DFT'
  write(io,'(a)')'DFT Functional:'
  write(io,'(a)')' unknown'
  write(io,'(a)')'Pseudopotential'
  write(io,'(a)')' unknown'
  call skipio(11,1)
  write(io,'(a)')'Plane wave cutoff (au)'
  read(11,*)ecut
  write(io,*)ecut
  call skipio(11,1)
  write(io,'(a)')'Spin polarized:'
  read(11,*)lsda
  write(io,*)lsda
  if(lsda) then
   nspin=2
  else
   nspin=1
  endif
  call skipio(11,1)
  write(io,'(a)')'Total energy (au per primitive cell)'
  read(11,*) etot
  write(io,*)etot
  call skipio(11,1)
  write(io,'(a)')'Kinetic energy (au per primitive cell)'
  read(11,*)etot
  write(io,*)etot
  call skipio(11,1)
  write(io,'(a)')'Local potential energy (au per primitive cell)'
  read(11,*)etot
  write(io,*)etot
  call skipio(11,1)
  write(io,'(a)')'Non local potential energy(au per primitive cell)'
  read(11,*)etot
  write(io,*)etot
  call skipio(11,1)
  write(io,'(a)')'Electron-electron energy (au per primitive cell)'
  read(11,*)etot
  write(io,*)etot
  call skipio(11,1)
  write(io,'(a)')'Ion-ion energy (au per primitive cell)'
  read(11,*)ewld
  write(io,*)ewld
  call skipio(11,1)
  write(io,'(a)')'number of electrons per primitive cell'
  read(11,*)nelec
  write(io,*)nelec
  call skipio(11,4)
  write(io,'(a)')' '
  write(io,'(a)')'GEOMETRY'
  write(io,'(a)')'-------- '
  write(io,'(a)')'Number of atoms per primitive cell '
  read(11,*)nat
  write(io,*)nat
  call skipio(11,1)
  write(io,'(a)')'Atomic number and position of the atoms(au) '
  do na=1,nat
   read(11,*)atomno,tau(:)
   write(io,'(i6,3f20.12)')atomno,tau(:)
  enddo
  call skipio(11,1)
  write(io,'(a)')'Primitive lattice vectors (au) '
  read(11,*) at(:,1)
  read(11,*) at(:,2)
  read(11,*) at(:,3)
  write(io,'(3f20.12)')at(:,1)
  write(io,'(3f20.12)')at(:,2)
  write(io,'(3f20.12)')at(:,3)

  call skipio(11,4)
  write(io,'(a)')' '
  write(io,'(a)')'G VECTORS'
  write(io,'(a)')'---------'
  write(io,'(a)')'Number of G-vectors'

  read(11,*)npwx
  write(io,*)npwx
  write(io,'(a)')'gx gy gz (au)'
  call skipio(11,1)

  allocate(g1(3,npwx),g2(npwx),indx(npwx),evc1(npwx) )

  do ig=1,npwx
   read(11,*)g1(:,ig)
   g2(ig)=dot_product(g1(:,ig),g1(:,ig))
  enddo
  call create_index(g2,indx)

  do ig=1,npwx
   write(io,'(3f20.12)')g1(:,indx(ig))
  enddo

  write(io,'(a)')' '
  write(io,'(a)')'WAVE FUNCTION'
  write(io,'(a)')'-------------'
  write(io,'(a)')'Number of k-points'
  call skipio(11,4)
  read(11,*)nks
  write(io,*)nks

  do ik=1,nks
   call skipio(11,1)
   write(io,'(a)')'k-point # ; # of bands (up spin/down spin) ; k-point &
    &coords (au)'
   read(11,*)i,nbnd(1),nbnd(2),xk(:)
    write(io,'(3i4,3f20.16)')ik,nbnd(1),nbnd(2),xk(:)
   do ispin=1,nspin
    do ibnd=1,nbnd(ispin)
     call skipio(11,1)
     write(io,'(a)')'band, spin, eigenvalue (au)'
     read(11,*)i,i,et
     write(io,'(2i5,f20.12)')ibnd,ispin,et
     call skipio(11,1)
! Read in the wave functions
     write(io,'(a)') 'EIGENVECTOR COEFFICIENTS'
     do ig=1, npwx
      read(11,*) evc1(ig)
     enddo
     do ig=1, npwx
      write(io,*) evc1(indx(ig))
     enddo
    enddo
   enddo
  enddo

  close(11,status='delete')

END PROGRAM pwfn_merge


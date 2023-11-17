 PROGRAM plot_mpc
!-----------------------------------------------------------------------------!
! PLOT_MPC                                                                    !
! Plotting program for density.data.                                          !
!                                                                             !
! Plots in real space the following quantities stored in Fourier space:       !
!                                                                             !
! 1. the SCF density calculated from the input trial wave function            !
! 2. 1/r (NOT CURRENTLY IMPLEMENTED)                                          !
!                                                                             !
! MDT 5.2006                                                                  !
!-----------------------------------------------------------------------------!
 USE dsp
 USE esdf
 USE run_control
 USE format_utils, ONLY : i2s
 IMPLICIT NONE

 INTEGER ialloc,ierr,npoint1,npoint2,npoint3,iplot,iset,n1,n2,n3
 INTEGER,PARAMETER :: io=10
 REAL(dp) treal,timag
 REAL(dp),PARAMETER :: pi=3.14159265358979324d0,twopi=2.d0*pi,&
  &one_over_twopi=1.d0/twopi
 LOGICAL file_present

! Basic info
 INTEGER nbasis
 REAL(dp) pa1(3),pa2(3),pa3(3),pb1(3),pb2(3),pb3(3),volume
 REAL(dp),ALLOCATABLE :: basis(:,:)

! Plot description

 INTEGER nlines,plot_dim
 REAL(dp) apoint(3),bpoint(3),cpoint(3),dpoint(3),xprod(3),r_ab(3),r_ac(3),&
  &r_ad(3),step_ab(3),step_ac(3),step_ad(3),r(3),r1(3),step_ab_length

! G vector sets
 INTEGER den_grange,den_ngvec
 INTEGER,ALLOCATABLE :: den_pgmap(:,:)
 REAL(dp),ALLOCATABLE :: den_gvec(:,:)
 COMPLEX(dp),ALLOCATABLE :: den_expigdotr(:),den_mwork(:,:)

! Density
 INTEGER ntypes
 COMPLEX(dp),ALLOCATABLE :: den_sc(:,:)

!----------------------------------------------------------------------------

  write(6,*)
  write(6,*)'PLOT_MPC : visualization of quantities from mpc.data file'
  write(6,*)'========================================================='
  write(6,*)

  inquire(file='mpc.data',exist=file_present)
  if(.not.file_present)then
   write(6,*)"No file 'mpc.data' containing the expectation value data."
   write(6,*)
   call errstop('PLOT_MPC','Quitting.')
  endif

  call read_mpc

  inquire(file='input',exist=file_present)
  if(.not.file_present)then
   write(6,*)"No file 'input' containing the controlling plot_expval block."
   write(6,*)
   call errstop('PLOT_MPC','Quitting.')
  endif

  call read_input

  call talk_to_user

  call write_mpc

  stop


CONTAINS


 SUBROUTINE read_mpc
!-----------------------------------------------------------------------------!
! Call routine to determine what data blocks mpc.data contains, then          !
! call the relevant additional routines to read these blocks.                 !
!-----------------------------------------------------------------------------!
 IMPLICIT NONE

 INTEGER i,j,k,ialloc,ierr

 write(6,*)'Found mpc.data file. Reading data.'
 write(6,*)

 open(io,file='mpc.data',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_MPC','Problem opening mpc.data file.')

! First read in DENSITY DATA part of mpc.data

! Read in and check the primitive lattice translation vectors :
 call skip(io,9)
 read(io,*,err=6,end=6)pa1
 read(io,*,err=6,end=6)pa2
 read(io,*,err=6,end=6)pa3
 volume=pa1(1)*pa2(2)*pa3(3)+pa1(2)*pa2(3)*pa3(1)+pa1(3)*pa2(1)*pa3(2)-&
  &pa1(3)*pa2(2)*pa3(1)-pa1(1)*pa2(3)*pa3(2)-pa1(2)*pa2(1)*pa3(3)
 pb1(1)=twopi*(pa2(2)*pa3(3)-pa2(3)*pa3(2))/volume
 pb1(2)=twopi*(pa2(3)*pa3(1)-pa2(1)*pa3(3))/volume
 pb1(3)=twopi*(pa2(1)*pa3(2)-pa3(1)*pa2(2))/volume
 pb2(1)=twopi*(pa3(2)*pa1(3)-pa1(2)*pa3(3))/volume
 pb2(2)=twopi*(pa1(1)*pa3(3)-pa3(1)*pa1(3))/volume
 pb2(3)=twopi*(pa3(1)*pa1(2)-pa1(1)*pa3(2))/volume
 pb3(1)=twopi*(pa1(2)*pa2(3)-pa2(2)*pa1(3))/volume
 pb3(2)=twopi*(pa2(1)*pa1(3)-pa1(1)*pa2(3))/volume
 pb3(3)=twopi*(pa1(1)*pa2(2)-pa1(2)*pa2(1))/volume

! Read in number of atoms in basis
 read(io,*,err=6,end=6)
 read(io,*,err=6,end=6)nbasis

! Read in the positions of the basis ions ( in atomic units ).
 read(io,*,err=6,end=6)
 allocate(basis(3,nbasis),stat=ialloc)
 if(ialloc/=0)call errstop('READ_MPC','Allocation problem with basis.')
 do i=1,nbasis
  read(io,*,err=6,end=6)k,(basis(j,i),j=1,3)
 enddo

! Read in number of particle types (1=electrons,2=electrons/holes)

 call skip(io,3)
 read(io,*,err=6,end=6)ntypes

! Read in number of G vectors

 read(io,*,err=6,end=6)
 read(io,*,err=6,end=6)den_ngvec

! Allocate arrays dependent on number of G vectors
 allocate(den_gvec(3,den_ngvec),den_sc(den_ngvec,ntypes),stat=ialloc)
 if(ialloc/=0)call errstop('READ_MPC','Allocation problem for density arrays.')

! Read in G vectors
 read(io,*,err=6,end=6)
 do i=1,den_ngvec
  read(io,*,err=6,end=6)den_gvec(1,i),den_gvec(2,i),den_gvec(3,i)
 enddo

! Read in charge density for each particle type

 do j=1,ntypes
  call skip(io,2)
  read(io,*,err=6,end=6)i
  if(i/=j)write(6,*)'Error reading particle type of set ',trim(i2s(j)),'.'
  read(io,*,err=6,end=6)
  do i=1,den_ngvec
   read(io,*,err=6,end=6)treal,timag
   den_sc(i,j)=cmplx(treal,timag,kind=dp)
  enddo
 enddo ! particle types

 call skip(io,3)

! Now read in EEPOT DATA part of mpc.data

!---------------------------------------------------------------------------!
! TO BE IMPLEMENTED                                                         !
!---------------------------------------------------------------------------!

 return

6 call errstop('READ_MPC','Problem reading DENSITY DATA section of mpc.data &
   &file.')

 END SUBROUTINE read_mpc


 SUBROUTINE read_input
 IMPLICIT NONE
 INTEGER j
 CHARACTER(16) plot_dim_string

 call esdf_init('input') ! i.e. read the input file
 call esdf_warnout

! Read the plot_expval block
 npoint1=0 ; npoint2=0 ; npoint3=0
 apoint=0.d0 ; bpoint=0.d0 ; cpoint=0.d0 ; dpoint=0.d0

 if(esdf_block('plot_expval',nlines))then

  write(6,*)'Found input file. Reading plot_expval block.'
  write(6,*)

! Item to plot
! Dimensionality
  read(block_data(1),*,err=2,end=3)plot_dim_string
  select case(trim(adjustl(plot_dim_string)))
   case('line','LINE','1','1D','1d','1-D','1-d')
    plot_dim=1
   case('plane','PLANE','2','2D','2d','2-D','2-d')
    plot_dim=2
   case('volume','VOLUME','3','3D','3d','3-D','3-d')
    plot_dim=3
   case default
    write(6,*)'The second line in the plot_expval block must be one of:'
    write(6,*)' line, plane, volume (or similar).'
    call errstop('READ_INPUT','Stopping.')
  end select
! No. of points to plot, starting point and ending point A & B [& C [& D]]
  select case(plot_dim)
   case(1)
    read(block_data(2),*,err=2,end=3)npoint1
    read(block_data(3),*,err=2,end=3)(apoint(j),j=1,3)
    read(block_data(4),*,err=2,end=3)(bpoint(j),j=1,3)
   case(2)
    read(block_data(2),*,err=2,end=3)npoint1,npoint2
    read(block_data(3),*,err=2,end=3)(apoint(j),j=1,3)
    read(block_data(4),*,err=2,end=3)(bpoint(j),j=1,3)
    read(block_data(5),*,err=2,end=3)(cpoint(j),j=1,3)
   case(3)
    read(block_data(2),*,err=2,end=3)npoint1,npoint2,npoint3
    read(block_data(3),*,err=2,end=3)(apoint(j),j=1,3)
    read(block_data(4),*,err=2,end=3)(bpoint(j),j=1,3)
    read(block_data(5),*,err=2,end=3)(cpoint(j),j=1,3)
    read(block_data(6),*,err=2,end=3)(dpoint(j),j=1,3)
  end select

! Check line lengths and collinearity
  r_ab=bpoint-apoint ; treal=sqrt(sum(r_ab(:)**2))
  if(treal==0.d0)call errstop('READ_INPUT','Length AB zero in plot_expval &
   &block.')
  if(plot_dim>1)then
   r_ac=cpoint-apoint ; treal=sqrt(sum(r_ac(:)**2))
   if(treal==0.d0)call errstop('READ_INPUT','Length AC zero in plot_expval &
    &block.')
   xprod(1)=r_ab(2)*r_ac(3)-r_ab(3)*r_ac(2)
   xprod(2)=r_ab(3)*r_ac(1)-r_ab(1)*r_ac(3)
   xprod(3)=r_ab(1)*r_ac(2)-r_ab(2)*r_ac(1)
   if(all(abs(xprod)<1.d-13))call errstop('READ_INPUT',&
    &'Collinear points A, B and C in plot_expval block.')
  endif
  if(plot_dim>2)then
   r_ad=dpoint-apoint ; treal=sqrt(sum(r_ad(:)**2))
   if(treal==0.d0)call errstop('READ_INPUT','Length AD zero in plot_expval &
    &block.')
   xprod(1)=r_ab(2)*r_ad(3)-r_ab(3)*r_ad(2)
   xprod(2)=r_ab(3)*r_ad(1)-r_ab(1)*r_ad(3)
   xprod(3)=r_ab(1)*r_ad(2)-r_ab(2)*r_ad(1)
   if(all(abs(xprod)<1.d-13))call errstop('READ_INPUT',&
    &'Collinear points A, B and D in qmc_plot.')
   xprod(1)=r_ac(2)*r_ad(3)-r_ac(3)*r_ad(2)
   xprod(2)=r_ac(3)*r_ad(1)-r_ac(1)*r_ad(3)
   xprod(3)=r_ac(1)*r_ad(2)-r_ac(2)*r_ad(1)
   if(all(abs(xprod)<1.d-13))call errstop('READ_INPUT',&
    &'Collinear points A, C and D in plot_expval block.')
  endif

! Define steps along each line to get the required number of points
  step_ab=r_ab/dble(npoint1-1)
  if(plot_dim>1)step_ac=r_ac/dble(npoint2-1)
  if(plot_dim>2)step_ad=r_ad/dble(npoint3-1)

  select case(plot_dim)
   case(1)
    step_ab_length=sqrt(step_ab(1)**2+step_ab(2)**2+step_ab(3)**2)
   case default
    step_ab_length=0.d0
  end select

 else

  call errstop('READ_INPUT','The input file does not contain a plot_expval &
   &block')

 endif

 return

2 call errstop('READ_INPUT','Error reading plot_expval block.')
3 call errstop('READ_INPUT','Reached end of line reading plot_expval block.')

 END SUBROUTINE read_input


 SUBROUTINE talk_to_user
!-----------------------------------------------------------------------------!
! Discuss with the user precisely what he wants to do.                        !
!-----------------------------------------------------------------------------!
 IMPLICIT NONE
 LOGICAL :: mike_hasnt_removed_this=.true.

 if(mike_hasnt_removed_this)then

  write(6,*)'Only plot of density currently implemented. Assuming this is what&
   & you want.'
  write(6,*)
  iplot=1

 endif

! Do required setup based on user responses

 call setup_fourier_basis

 END SUBROUTINE talk_to_user


 SUBROUTINE write_mpc
 IMPLICIT NONE
 write(6,*)'Plotting density from mpc.data file.'
 write(6,*)
 select case(plot_dim)
  case(1)
   open(10,file='lineplot.dat',status='unknown',action='write',iostat=ierr)
   write(6,*)'Data written to lineplot.dat.'
  case(2)
   open(10,file='2Dplot.dat',status='unknown',action='write',iostat=ierr)
   write(6,*)'Data written to 2Dplot.dat.'
  case(3)
   open(10,file='3Dplot.dat',status='unknown',action='write',iostat=ierr)
   write(6,*)'Data written to 3Dplot.dat.'
  case default
   call errstop('WRITE_MPC','Invalid plot_dim code - bug.')
 end select
 if(ierr/=0)call errstop('WRITE_MPC','Error opening plot file for writing.')
 select case(iplot)
  case(1)
   call plot_density
  case(2)
   call plot_eepot
  case default
   call errstop('WRITE_MPC','Invalid code for MPC plot selection - bug.')
 end select
 END SUBROUTINE write_mpc


 SUBROUTINE plot_density
!-----------------------------------------------------------------------------!
! Plot the DENSITY block from mpc.data                                        !
!-----------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,k,g
 REAL(dp) x,den

 iset=1 ! Assuming only 1 particle type for the moment.

 select case (plot_dim)
  case(1) ! line

   if(r_ab(2)==0.d0.and.r_ab(3)==0.d0)then
    x=apoint(1)
   elseif(r_ab(1)==0.d0.and.r_ab(3)==0.d0)then
    x=apoint(2)
   elseif(r_ab(1)==0.d0.and.r_ab(2)==0.d0)then
    x=apoint(3)
   else
    x=0.d0
   endif
   r(:)=apoint(:)
   do i=1,npoint1
    den=0.d0
    call compute_fourier_basis(r)
    do g=1,den_ngvec
     den=den+real(den_sc(g,iset)*den_expigdotr(g),dp)
    enddo
    write(10,'(e18.10,1x,e18.10)')x,den
    r(:)=r(:)+step_ab(:)
    x=x+step_ab_length
   enddo

  case(2) ! plane

   r(:)=apoint(:) ; r1(:)=apoint(:)
   do j=1,npoint2
    do i=1,npoint1
     den=0.d0
     call compute_fourier_basis(r)
     do g=1,den_ngvec
      den=den+real(den_sc(g,iset)*den_expigdotr(g),dp)
     enddo
     write(10,'(e18.10,3(1x,e18.10))')r(:),den
     r(:)=r(:)+step_ab(:)
    enddo
    write(10,*)
    r1(:)=r1(:)+step_ac(:) ; r(:)=r1(:)
   enddo

  case(3) ! volume

   do k=0,npoint3-1
    do j=0,npoint2-1
     do i=1,npoint1
      r(:)=apoint(:)+dble(i-1)*step_ab(:)+dble(j)*step_ac(:)+dble(k)*step_ad(:)
      den=0.d0
      call compute_fourier_basis(r)
      do g=1,den_ngvec
       den=den+real(den_sc(g,iset)*den_expigdotr(g),dp)
      enddo
      write(10,'(e18.10,3(1x,e18.10))')r(:),den
     enddo
     write(10,*)
    enddo
   enddo

  case default
   call errstop('PLOT_DENSITY','Invalid dimensionality for plot.')

 end select

 END SUBROUTINE plot_density


 SUBROUTINE plot_eepot
 IMPLICIT NONE
 call errstop('PLOT_EEPOT','Not implemented yet.')
 END SUBROUTINE plot_eepot


 SUBROUTINE setup_fourier_basis
!-----------------------------------------------------------------------------!
! Perform required setup for calculation of Fourier basis.                    !
!-----------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i

 allocate(den_pgmap(3,den_ngvec),stat=ialloc)
 if(ialloc/=0)call errstop('SETUP_FOURIER_BASIS','Allocation problem : &
  &den_pgmap/grange arrays.')
 allocate(den_expigdotr(den_ngvec),stat=ialloc)
 if(ialloc/=0)call errstop('SETUP_FOURIER_BASIS','Allocation error: &
  &den_expigdotr array.')

! Map G vectors onto repeats of primitive lattice
 den_grange=0

 do i=1,den_ngvec

  n1=nint(one_over_twopi*(den_gvec(1,i)*pa1(1)+ &
   &den_gvec(2,i)*pa1(2)+den_gvec(3,i)*pa1(3)))
  if(abs(n1)>den_grange)den_grange=abs(n1)
  den_pgmap(1,i)=n1

  n2=nint(one_over_twopi*(den_gvec(1,i)*pa2(1)+ &
   &den_gvec(2,i)*pa2(2)+den_gvec(3,i)*pa2(3)))
  if(abs(n2)>den_grange)den_grange=abs(n2)
  den_pgmap(2,i)=n2

  n3=nint(one_over_twopi*(den_gvec(1,i)*pa3(1)+ &
   &den_gvec(2,i)*pa3(2)+den_gvec(3,i)*pa3(3)))
  if(abs(n3)>den_grange)den_grange=abs(n3)
  den_pgmap(3,i)=n3

 enddo ! G vectors in set

 if(den_grange<1)then
  call errstop('SETUP_FOURIER_BASIS','Density G vector sets in density.data &
   &must include non-zero G.')
 endif

 allocate(den_mwork(3,-den_grange:den_grange),stat=ialloc)
 if(ialloc/=0)call errstop('SETUP_FOURIER_BASIS','Allocation problem : &
  &den_mwork')

 END SUBROUTINE setup_fourier_basis


 SUBROUTINE compute_fourier_basis(rvec)
!-----------------------------------------------------------------------------!
! Compute the required exp(-iG.r) terms used in the reciprocal space          !
! representations of expectation values.                                      !
!-----------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i
 REAL(dp),INTENT(in) :: rvec(3)
 REAL(dp) p(3),cos1,cos2,cos3,sin1,sin2,sin3

 p(1)=-pb1(1)*rvec(1)-pb1(2)*rvec(2)-pb1(3)*rvec(3)
 p(2)=-pb2(1)*rvec(1)-pb2(2)*rvec(2)-pb2(3)*rvec(3)
 p(3)=-pb3(1)*rvec(1)-pb3(2)*rvec(2)-pb3(3)*rvec(3)

 cos1=cos(p(1)) ; cos2=cos(p(2)) ; cos3=cos(p(3))
 sin1=sin(p(1)) ; sin2=sin(p(2)) ; sin3=sin(p(3))

 den_mwork(1:3,0)=1.d0
! Must have G/=0 in basis or will crash here
 den_mwork(1,1)=cmplx(cos1,sin1,kind=dp)
 den_mwork(2,1)=cmplx(cos2,sin2,kind=dp)
 den_mwork(3,1)=cmplx(cos3,sin3,kind=dp)
 den_mwork(1,-1)=cmplx(cos1,-sin1,kind=dp)
 den_mwork(2,-1)=cmplx(cos2,-sin2,kind=dp)
 den_mwork(3,-1)=cmplx(cos3,-sin3,kind=dp)

 do i=2,den_grange
  den_mwork(1,i)=den_mwork(1,1)*den_mwork(1,i-1)
  den_mwork(2,i)=den_mwork(2,1)*den_mwork(2,i-1)
  den_mwork(3,i)=den_mwork(3,1)*den_mwork(3,i-1)
  den_mwork(1,-i)=conjg(den_mwork(1,i))
  den_mwork(2,-i)=conjg(den_mwork(2,i))
  den_mwork(3,-i)=conjg(den_mwork(3,i))
 enddo

! Insert alternative faster version for real arrays - xxxx
 do i=1,den_ngvec
  n1=den_pgmap(1,i); n2=den_pgmap(2,i); n3=den_pgmap(3,i)
  den_expigdotr(i)=den_mwork(1,n1)*den_mwork(2,n2)*den_mwork(3,n3)
 enddo

 END SUBROUTINE compute_fourier_basis


END PROGRAM plot_mpc

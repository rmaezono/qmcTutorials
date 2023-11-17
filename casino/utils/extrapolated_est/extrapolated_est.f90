PROGRAM extrapolated_est
 !-------------------------------------------------------------!
 ! This program calculates the extrapolated estimator from 2   !
 ! data files, one DMC and one VMC. If both sets have          !
 ! errorbars, these will be combined in the usual way.         !
 !                                                             !
 ! One can produce the required input files by using           !
 ! PLOT_EXPVAL in a directory containing expval.data           !
 !-------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER, PARAMETER :: dp=KIND(1.d0)
 INTEGER :: io=8,ierr,i,ialloc,fileno
 INTEGER, ALLOCATABLE :: nlines(:)
 REAL(dp), ALLOCATABLE :: rbin(:,:),weight(:,:),new_weight(:),&
  &errorbars(:,:),extrap_err(:)
 REAL(dp) dummy(3)
 REAL(dp), PARAMETER :: tol=1.d-8
 CHARACTER(len=100) infile,tmpstring,dmcfilename,vmcfilename
 LOGICAL have_errorbars(2)

 allocate(nlines(2),stat=ialloc)
 if(ialloc/=0)then
  write(*,*)'Problem allocating array nlines.'
  stop
 endif
 have_errorbars(:)=.false.
 do fileno=1,2
  if(fileno==1)then
   write(*,*)'Please enter the name of the file containing the DMC &
    &data:'
   read(*,*) dmcfilename
   infile=dmcfilename
  endif
  if(fileno==2)then
   write(*,*)'Please enter the name of the file containing the VMC &
    &data:'
   read(*,*) vmcfilename
   infile=vmcfilename
  endif
  open(unit=io,file=trim(infile),status='old',iostat=ierr)
  if(ierr/=0)then
   write(*,*)'Problem reading file '//trim(infile)//'.'
   stop
  endif
  i=0
  do
   read(io,*,iostat=ierr)tmpstring
   if(ierr<0)exit
   i=i+1
   if(ierr>0)then
    write(*,*)'Problem reading file '//trim(infile)//'.'
    stop
   endif
  enddo
  nlines(fileno)=i  ! infile(fileno) thus has nlines
  if(nlines(fileno)/=nlines(1))then
   write(*,*)'The two files must have equal numbers of lines of data.'
   stop
  endif
  if(fileno==1)then
   allocate(rbin(2,nlines(fileno)),weight(2,nlines(fileno)),&
    &new_weight(nlines(fileno)),errorbars(2,nlines(fileno)),&
    &extrap_err(nlines(fileno)),stat=ialloc)
   if(ialloc/=0)then
    write(*,*)'Problem allocating arrays rbin, weight, errorbars&
     & and extrap_err.'
    stop
   endif
  endif ! fileno==1
  rewind(io)
  i=0
  do
   i=i+1
   read(io,*,iostat=ierr)dummy(1:3)           ! See if there are 2 or 3
   if(ierr/=0)then                            ! columns (i.e. errorbars?)
    have_errorbars(fileno)=.false.
    exit
   endif
   if(i==nlines(1).and.ierr==0)then
    have_errorbars(fileno)=.true.
    exit
   endif
  enddo ! First loop over lines
  rewind(io)
  do i=1,nlines(fileno)
   if(have_errorbars(fileno))then
    read(io,*,iostat=ierr)rbin(fileno,i),weight(fileno,i),&
     &errorbars(fileno,i)
   else
    read(io,*,iostat=ierr)rbin(fileno,i),weight(fileno,i)
    errorbars(fileno,i)=0.d0
   endif
   if(ierr/=0)then
    write(*,*)'Problem reading file '//trim(infile)//'. (3)'
    stop
   endif
  enddo ! Second loop over lines
  close(io)
 enddo ! Loop over different files
 if(have_errorbars(1).neqv.have_errorbars(2))then
  write(*,*)'Only one of the files has a complete set of errorbars &
   &specified.'
  stop
 endif
 write(*,*)
 if(have_errorbars(1))write(*,*)'Computing errorbars.'
 if(.not.have_errorbars(1))write(*,*)'Not computing any errorbars.'
 write(*,*)
 write(*,*)'Evaluating extrapolated estimator &
  &using files '//trim(dmcfilename)//' and '//trim(vmcfilename)//' for &
  &DMC and VMC data, respectively.'
 do i=1,nlines(1)
  new_weight(i)=0.d0
  if(((rbin(1,i)-rbin(2,i))**2)>tol**2)then
   write(*,*)'x-components of data (i.e., bin positions) &
    &should be the same across both files for extrapolation. The two &
    &sets differ on line number:', i
   write(*,*)'Difference in bin positions:',rbin(1,i)-rbin(2,i)
   write(*,*)'Tolerance                  :',tol
   stop
  endif
  new_weight(i)=2*weight(1,i)-weight(2,i)
  ! Use the standard method for combining errors
  if(have_errorbars(1))extrap_err(i)=sqrt(4.d0*(errorbars(1,i)**2)&
   &+errorbars(2,i)**2)
 enddo ! Loop over lines

 open(unit=io,file='lineplot_ext.dat',status='new',iostat=ierr)
 if(ierr/=0)then
  write(*,*)'Problem opening lineplot_ext.dat for output (maybe the &
   &file already exists?).'
  stop
 endif
 do i=1,nlines(1)
  if(have_errorbars(1))then
   write(io,*,iostat=ierr)rbin(1,i),new_weight(i),extrap_err(i)
  else
   write(io,*,iostat=ierr)rbin(1,i),new_weight(i)
  endif
  if(ierr/=0)then
   write(*,*)'Problem writing to lineplot_ext.dat.'
   stop
  endif
 enddo
 close(io)
 write(*,*)'Output written to lineplot_ext.dat.'
 write(*,*)
 deallocate(nlines,rbin,weight,new_weight,errorbars,extrap_err)
END PROGRAM


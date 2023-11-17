PROGRAM findmax
! Find maximum abs. value in column ncol of file filename
! Empty lines are skipped
 IMPLICIT NONE
 INTEGER,PARAMETER :: dp=selected_real_kind(kind(1.d0))
 INTEGER ncol,ierr
 REAL(dp) maxv,currv
 CHARACTER(80) filename
 CHARACTER(20),ALLOCATABLE :: unused(:)

 read(*,*)filename,ncol
 if(ncol<1)then
  write(6,*)'ERROR' ; stop
 endif
 allocate(unused(ncol-1))

 open(unit=10,file=filename,status='old',iostat=ierr)
 if(ierr/=0)then
  write(6,*)'ERROR' ; stop
 endif

 maxv=0.d0
 do
  if(ncol==1)then
   read(10,*,err=5,end=10)currv
  else
   read(10,*,err=5,end=10)unused(1:ncol-1),currv
  endif
  currv=abs(currv)
  if(currv>maxv)maxv=currv
5 cycle
 enddo

10 close(10)
 write(6,'(es19.12)')maxv

END PROGRAM findmax

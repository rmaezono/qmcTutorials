!-------------------------------------------------------------------------!
! QUAD_FIT by Neil Drummond, Nov 2002.                                    !
!                                                                         !
! Program for carrying out a quadratic fit to a set of data, in order     !
! to find a local extremum.                                               !
!                                                                         !
! Data file should consist of two or three columns of data, which         !
! are the (x,y) or (x,y +/- delta_y) points to which the fit is           !
! made. If two columns are given then the program will perform a least-   !
! squares fit to the (x,y) data; if three columns are given then a        !
! chi^2 fit to the (x,y +/- delta_y) will be performed.                   !
!-------------------------------------------------------------------------!


MODULE qf_utils
!----------------------------------------------------------------!
! This module contains some subprograms used by quad_fit to read !
! from the data file and carry out the chi^2 fit.                !
!----------------------------------------------------------------!
 IMPLICIT NONE
 PRIVATE
 PUBLIC check_file,read_file,perform_fit,display_fitting_info,rescale_errors, &
  &i2s,dp,errstop,find_centre
 INTEGER,PARAMETER :: dp=kind(1.d0)


CONTAINS


 SUBROUTINE check_file(io,fname,have_errors,N)
!-----------------------------------------------------------------------------!
! This subroutine looks at the supplied data file and works out how many data !
! items there are and whether they have error bars.                           !
!-----------------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: io
  CHARACTER(*),INTENT(in) :: fname
  INTEGER,INTENT(out) :: N
  LOGICAL,INTENT(out) :: have_errors
  INTEGER ierr
  REAL(dp) dummy(3)
  CHARACTER(180) char180
  N=0
  do
   read(io,'(a)',iostat=ierr)char180
   if(ierr>0)call errstop('CHECK_FILE','Error reading '//trim(fname)//'.')
   if(ierr<0)exit
   char180=adjustl(char180)
   if(index(char180,'#')/=1)then
    N=N+1
    read(char180,*,iostat=ierr)dummy(1:3)
    if(N==1)then
     have_errors=(ierr==0)
    else
     if((ierr==0).neqv.have_errors)call errstop('CHECK_FILE', &
      &'Some lines of data have error bars.')
    endif ! N=1
   endif ! Not a comment.
  enddo
  rewind(io)
 END SUBROUTINE check_file


 SUBROUTINE read_file(io,N,x,y,rec_delta_y_sq,fname,have_errors)
!-----------------------------------------------------------------------------!
! Read in data from file.                                                     !
! If we're doing a least-squares fit then just make all the error bars in the !
! chi^2 fit constant, equal to 1. Then the chi^2 procedure reduces to the     !
! least-squares procedure.                                                    !
!-----------------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: io,N
  LOGICAL,INTENT(in) :: have_errors
  CHARACTER(*),INTENT(in) :: fname
  REAL(dp),INTENT(out) :: x(N),y(N),rec_delta_y_sq(N)
  INTEGER ierr,i
  do i=1,N
   if(have_errors)then
    read(io,*,iostat=ierr)x(i),y(i),rec_delta_y_sq(i)
    if(ierr/=0)call errstop('READ_FILE','Error reading file '//trim(fname)//'.')
    if(rec_delta_y_sq(i)<=0.d0)call errstop('READ_FILE', &
     &'Error bars are not all positive in file '//trim(fname)//'.')
    rec_delta_y_sq(i)=1.d0/rec_delta_y_sq(i)**2
   else
    read(io,*,iostat=ierr)x(i),y(i)
    if(ierr/=0)call errstop('READ_FILE','Error reading file '//trim(fname)//'.')
   endif ! have_errors
  enddo ! i
  if(.not.have_errors)rec_delta_y_sq=1.d0
 END SUBROUTINE read_file


 SUBROUTINE perform_fit(N,x,y,rec_delta_y_sq,a,Minv,chi2)
!-------------------------------------------------------------------------!
! This subroutine carries out a chi^2 fit of a quadratic to the N         !
! (x,y +/- delta_y) data points. It returns errors in the fitting         !
! coefficients as well as the coefficients themselves and the (minimised) !
! chi^2 function.  See accompanying notes for info on the calculation.    !
!-------------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: N
  REAL(dp),INTENT(in) :: x(N),y(N),rec_delta_y_sq(N)
  REAL(dp),INTENT(out) :: a(3),Minv(3,3),chi2
  INTEGER i
  REAL(dp) M(3,3),c(3),x_ov_dy2,x2_ov_dy2,x3_ov_dy2,x4_ov_dy2,yfit

! Construct the matrix M and vector c.
  M=0.d0 ; c=0.d0
  do i=1,N
   x_ov_dy2=x(i)*rec_delta_y_sq(i) ; x2_ov_dy2=x(i)*x_ov_dy2
   x3_ov_dy2=x(i)*x2_ov_dy2        ; x4_ov_dy2=x(i)*x3_ov_dy2
   M(1,1)=M(1,1)+x4_ov_dy2
   M(1,2)=M(1,2)+x3_ov_dy2
   M(1,3)=M(1,3)+x2_ov_dy2
   M(2,3)=M(2,3)+x_ov_dy2
   M(3,3)=M(3,3)+rec_delta_y_sq(i)
   c(1)=c(1)+y(i)*x2_ov_dy2
   c(2)=c(2)+y(i)*x_ov_dy2
   c(3)=c(3)+y(i)*rec_delta_y_sq(i)
  enddo ! i
  M(2,1)=M(1,2)
  M(2,2)=M(1,3) ; M(3,1)=M(1,3)
  M(3,2)=M(2,3)

! Invert matrix M.
  call inv33(M,Minv)

! Evaluate the polynomial coefficients.
  a=matmul(Minv,c)

! Evaluate chi^2 function for the optimized parameters.
  chi2=0.d0
  do i=1,N
   yfit=(a(1)*x(i)+a(2))*x(i)+a(3)
   chi2=chi2+(y(i)-yfit)**2*rec_delta_y_sq(i)
  enddo ! i

 END SUBROUTINE perform_fit


 SUBROUTINE display_fitting_info(a,Minv,chi2,x_c)
!----------------------------------------------------------!
! Write out the fitting parameters, etc., with error bars. !
!----------------------------------------------------------!
  IMPLICIT NONE
  REAL(dp),INTENT(in) :: a(3),Minv(3,3),chi2,x_c
  INTEGER j
  REAL(dp) delta_a,z,delta_z,y0,delta_y0,ea(5)
  CHARACTER(25) :: as(3),eas(5),x_cs
  write(6,*)'Fitted quadratic is y(x) = a1.(x-x_c)^2 + a2.(x-x_c) + a3, where:'
  write(6,*)
  do j=1,3
   delta_a=sqrt(Minv(j,j))
   write(6,'(1x,a,2(es24.16,a))')' a'//trim(i2s(j))//' = ',a(j),' +/- ', &
    &delta_a,'.'
  enddo ! j
  write(6,*)
  write(6,'(1x,a,es24.16,a)')'The chi^2 function is: ',chi2,'.'
  write(6,*)
  if(a(1)/=0.d0)then
   call get_extremum(a,Minv,z,delta_z,y0,delta_y0)
   if(a(1)<0.d0)then
    write(6,*)'Fitted quadratic has a maximum at: '
   else
    write(6,*)'Fitted quadratic has a minimum at: '
   endif ! a(1)<0
   write(6,'(1x,es24.16,a,es24.16)')z+x_c,' +/- ',delta_z
   write(6,*)
   write(6,*)'The function value at this point is: '
   write(6,'(1x,es24.16,a,es24.16)')y0,' +/- ',delta_y0
  else
   write(6,*)'You appear to have a straight line.'
   write(6,*)'Hence a local extremum can''t be found.'
  endif ! a(1)/=0
  write(6,*)

! Write out curve to paste into XMGrace.
  write(x_cs,'(es24.16)')x_c
  x_cs=adjustl(x_cs)
  do j=1,3
   write(as(j),'(es24.16)')a(j)
   if(j>1.and.a(j)>=0.d0)then
    as(j)='+'//adjustl(as(j))
   else
    as(j)=adjustl(as(j))
   endif ! need to add "+".
  enddo ! j
  write(6,*)'Fitted quadratic (in a suitable format for pasting into XMGrace):'
  write(6,*)
  write(6,*)'y=('//trim(as(1))//'*(x-'//trim(x_cs)//')'//trim(as(2)) &
   &//')*(x-'//trim(x_cs)//')'//trim(as(3))
  write(6,*)

! Write out error bars to paste into XMGrace.
  ea(1)=Minv(1,1) ; ea(2)=2.d0*Minv(1,2) ; ea(3)=Minv(2,2)+2.d0*Minv(1,3)
  ea(4)=2.d0*Minv(2,3) ; ea(5)=Minv(3,3)
  do j=1,5
   write(eas(j),'(es24.16)')ea(j)
   if(j>1.and.ea(j)>=0.d0)then
    eas(j)='+'//adjustl(eas(j))
   else
    eas(j)=adjustl(eas(j))
   endif ! need to add "+".
  enddo ! j
  write(6,*)'Standard error (in a suitable format for pasting into XMGrace):'
  write(6,*)
  write(6,*)'y1=sqrt(((('//trim(eas(1))//'*(x-'//trim(x_cs)//')'//trim(eas(2)) &
   &//')*(x-'//trim(x_cs)//')'//trim(eas(3))//')*(x-'//trim(x_cs)//')' &
   &//trim(eas(4))//')*(x-'//trim(x_cs)//')'//trim(eas(5)),')'
  write(6,*)

 END SUBROUTINE display_fitting_info


 SUBROUTINE get_extremum(a,Minv,z,delta_z,y0,delta_y0)
!-----------------------------------------------------------------------------!
! Evaluate the extremum of the quadratic, with error bars.  z is the location !
! of the extremum and y0 is the value of the quadratic at the extremum.       !
!-----------------------------------------------------------------------------!
  IMPLICIT NONE
  REAL(dp),INTENT(in) :: a(3),Minv(3,3)
  REAL(dp),INTENT(out) :: z,delta_z,y0,delta_y0
  REAL(dp) p(3)
  if(a(1)==0.d0)call errstop('GET_EXTREMUM', &
   &'Error: fitted "quadratic" is in fact linear.')
  z=-a(2)/(2.d0*a(1))
  p=(/0.5d0*a(2)/a(1)**2,-0.5d0/a(1),0.d0/) ! See notes.
  delta_z=sqrt(pT_Minv_p(p,Minv))
  y0=-0.25d0*a(2)**2/a(1)+a(3)
  p=(/0.25d0*a(2)**2/a(1)**2,-0.5d0*a(2)/a(1),1.d0/) ! See notes.
  delta_y0=sqrt(pT_Minv_p(p,Minv))
 END SUBROUTINE get_extremum


 REAL(dp) FUNCTION pT_Minv_p(p,Minv)
!---------------------------------------------!
! Evaluate p^T M^(-1) p.  NB, M is symmetric. !
!---------------------------------------------!
  IMPLICIT NONE
  REAL(dp),INTENT(in) :: p(3),Minv(3,3)
  pT_Minv_p=p(1)**2*Minv(1,1)+p(2)**2*Minv(2,2)+p(3)**2*Minv(3,3) &
   &+2.d0*(p(1)*(p(2)*Minv(1,2)+p(3)*Minv(1,3))+p(2)*p(3)*Minv(2,3))
 END FUNCTION pT_Minv_p


 SUBROUTINE rescale_errors(N,rec_delta_y_sq,chi2,Minv)
!--------------------------------------------------------------------------!
! Rescale the error bars so that the chi^2 function is equal to N.  Update !
! chi^2 function and Minv matrix accordingly.                              !
!--------------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: N
  REAL(dp),INTENT(inout) :: rec_delta_y_sq(N),chi2,Minv(3,3)
  rec_delta_y_sq(1:N)=rec_delta_y_sq(1:N)*(dble(N)/chi2)
  Minv=Minv*(chi2/dble(N))
  chi2=dble(N)
 END SUBROUTINE rescale_errors


 SUBROUTINE find_centre(N,x,rec_delta_y_sq,x_c)
!--------------------------------!
! Determine the average x value. !
!--------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: N
  REAL(dp),INTENT(in) :: rec_delta_y_sq(N)
  REAL(dp),INTENT(inout) :: x(N)
  REAL(dp),INTENT(out) :: x_c
  x_c=dot_product(rec_delta_y_sq,x)/sum(rec_delta_y_sq)
  x=x-x_c
 END SUBROUTINE find_centre


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


 SUBROUTINE inv33(M,Minv)
!----------------------!
! Inverts 3x3 matrices !
!----------------------!
  IMPLICIT NONE
  REAL(dp),INTENT(in) :: M(3,3)
  REAL(dp),INTENT(out) :: Minv(3,3)
  REAL(dp) d
  d=M(1,1)*(M(2,2)*M(3,3)-M(2,3)*M(3,2))+ &
   &M(2,1)*(M(3,2)*M(1,3)-M(1,2)*M(3,3))+ &
   &M(3,1)*(M(1,2)*M(2,3)-M(1,3)*M(2,2))
  if(d==0.d0)call errstop('INV33','Trying to invert a singular matrix...')
  d=1.d0/d
  Minv(1,1)=(M(2,2)*M(3,3)-M(2,3)*M(3,2))*d
  Minv(1,2)=(M(3,2)*M(1,3)-M(1,2)*M(3,3))*d
  Minv(1,3)=(M(1,2)*M(2,3)-M(1,3)*M(2,2))*d
  Minv(2,1)=(M(3,1)*M(2,3)-M(2,1)*M(3,3))*d
  Minv(2,2)=(M(1,1)*M(3,3)-M(3,1)*M(1,3))*d
  Minv(2,3)=(M(2,1)*M(1,3)-M(1,1)*M(2,3))*d
  Minv(3,1)=(M(2,1)*M(3,2)-M(2,2)*M(3,1))*d
  Minv(3,2)=(M(3,1)*M(1,2)-M(1,1)*M(3,2))*d
  Minv(3,3)=(M(1,1)*M(2,2)-M(1,2)*M(2,1))*d
 END SUBROUTINE inv33


 SUBROUTINE errstop(sub,message)
!---------------------------!
! Report an error and stop. !
!---------------------------!
  IMPLICIT NONE
  CHARACTER(*),INTENT(in) :: sub,message
  write(6,*)
  write(6,*)'ERROR in subroutine '//trim(adjustl(sub))//'.'
  write(6,*)
  call wordwrap(trim(adjustl(message)))
  write(6,*)
  stop
 END SUBROUTINE errstop


 SUBROUTINE wordwrap(text,unit_in,linelength_in)
!-------------------------------------------------------------------------!
! This subroutine prints out the contents of the character string 'text', !
! ensuring that line breaks only occur at space characters.  The output   !
! is written to unit unit_in if this parameter is supplied; otherwise the !
! output is written to unit 6.  The maximum length of each line is given  !
! by linelength_in if this is supplied; otherwise the default line length !
! is 79 characters.                                                       !
!-------------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in),OPTIONAL :: unit_in,linelength_in
  CHARACTER(*),INTENT(in) :: text
  CHARACTER(260) :: temp
  INTEGER :: i,unit,lentext,startpoint,stoppoint,lastpos,linelength
  if(present(unit_in))then
   unit=unit_in
  else
   unit=6
  endif ! unit supplied.
  lentext=len(trim(text))
  if(lentext<1)then
   write(unit,*)
   return
  endif ! No text
  if(present(linelength_in))then
   if(linelength_in>=2)then
    linelength=linelength_in
   else
    linelength=2
   endif ! sensible line-length supplied.
  else
   linelength=79
  endif ! linelength present.
  startpoint=1
  i=0
  do
   i=i+1
   stoppoint=startpoint+linelength-1
   if(stoppoint<=lentext)then
    lastpos=index(trim(text(startpoint:stoppoint))," ",.true.)
    if(lastpos>0)stoppoint=startpoint+lastpos-1
   else
    stoppoint=lentext
   endif ! stoppoint <= length of text
   if(i==1)then
! Allow the user to indent the first line, if (s)he wishes.
    temp=text(startpoint:stoppoint) ! or else pathscale f90 fails to compile
    write(unit,*)trim(temp)
   else
    temp=text(startpoint:stoppoint) ! or else pathscale f90 fails to compile
    write(unit,*)trim(adjustl(temp))
   endif ! i=1
   if(stoppoint==lentext)then
    exit
   else
    startpoint=stoppoint+1
   endif ! Finished text?
  enddo ! Loop over lines.
 END SUBROUTINE wordwrap


END MODULE qf_utils



PROGRAM quad_fit
!---------------------------!
! Main program starts here. !
!---------------------------!
 USE qf_utils,ONLY : check_file,read_file,perform_fit,display_fitting_info, &
  &rescale_errors,i2s,dp,errstop,find_centre
 IMPLICIT NONE
 INTEGER N,ierr
 INTEGER,PARAMETER :: io=8
 REAL(dp),ALLOCATABLE :: x(:),y(:),rec_delta_y_sq(:)
 REAL(dp) chi2,a(3),Minv(3,3),x_c
 CHARACTER(72) fname
 LOGICAL have_errors,file_exists

 write(6,*)
 write(6,*)'QUAD_FIT'
 write(6,*)'========'
 write(6,*)

! Get user to enter name of data file.
 do
  write(6,*)'Please enter name of data file:'
  read(5,*,iostat=ierr)fname
  if(ierr/=0)fname=''
  fname=adjustl(fname)
  write(6,*)
  inquire(file=trim(fname),exist=file_exists)
  if(file_exists)exit
  write(6,*)'That file does not seem to exist.  Please try again.'
 enddo

! Open and read the data file.
 open(unit=io,file=trim(fname),iostat=ierr,status='old')
 if(ierr/=0)call errstop('MAIN','Cannot open file '//trim(fname)//'.  Sorry.')
 call check_file(io,fname,have_errors,N)
 write(6,*)'There seem to be '//trim(i2s(N))//' lines in file '//trim(fname) &
  &//'.'
 if(have_errors)then
  write(6,*)'Format appears to be x, y, delta_y.'
 else
  write(6,*)'Format appears to be x, y.'
 endif ! have_errors
 allocate(x(N),y(N),rec_delta_y_sq(N),stat=ierr)
 if(ierr/=0)call errstop('MAIN','Allocation error.')
 call read_file(io,N,x,y,rec_delta_y_sq,fname,have_errors)
 close(io)
 write(6,*)'Data read in.'
 write(6,*)

! Determine the "centre" of the data.  Offset the data accordingly.
 call find_centre(N,x,rec_delta_y_sq,x_c)
 write(6,*)'Offset x_c: ',x_c
 write(6,*)

! Carry out the chi^2 fit to determine the parameters of the
! fitted quadratic and their error bars.
 call perform_fit(N,x,y,rec_delta_y_sq,a,Minv,chi2)

! If we don't have error bars, choose them so that chi^2=N.  Instead of
! performing fit again, just rescale Minv etc appropriately.
 if(.not.have_errors)call rescale_errors(N,rec_delta_y_sq,chi2,Minv)

! Display resulting information about the fitted polynomial.
 call display_fitting_info(a,Minv,chi2,x_c)

 deallocate(x,y,rec_delta_y_sq)

 write(6,*)'Program finished.'
 write(6,*)

END PROGRAM quad_fit




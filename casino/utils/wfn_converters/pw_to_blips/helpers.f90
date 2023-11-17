MODULE helpers
!---------------------------------------------------------------!
! Various helper routines needed by the blip program            !
!---------------------------------------------------------------!
 IMPLICIT NONE
! Double precision kind.
 INTEGER,PARAMETER :: dp=kind(1.d0)

! Some constants.
 REAL(dp),PARAMETER :: pi=3.14159265358979324d0
 REAL(dp),PARAMETER :: twopi=2.d0*pi
 REAL(dp),PARAMETER :: rec_twopi=1.d0/twopi
 REAL(dp),PARAMETER :: eps=1.d-10
 REAL(dp),PARAMETER :: third=1.d0/3.d0
 REAL(dp),PARAMETER :: two_thirds=2.d0/3.d0

 COMPLEX(dp),PARAMETER :: czero=(0.d0,0.d0)
 COMPLEX(dp),PARAMETER :: iunity=(0.d0,1.d0)
 COMPLEX(dp),PARAMETER :: twoi=(0.d0,2.d0)

! Precision for Simpson's rule - can be changed before calls to multiple
! integrals.  Precision for the KE and norm^2 integrands.
 REAL(dp) :: eps_simp=1.d-8
 REAL(dp),PARAMETER :: eps_ke=1.d-3
 REAL(dp),PARAMETER :: eps_norm2=1.d-4

 PUBLIC

CONTAINS


 SUBROUTINE skipio(io,n)
!-----------------------------------------------!
! Skip n lines from the file opened on unit io. !
!-----------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: io,n
  INTEGER i,ierr
  do i=1,n
   read(io,*,iostat=ierr)
   if(ierr>0)then
    call errstop('SKIPIO','Error reading file to skip lines.')
   elseif(ierr<0)then
    call errstop('SKIPIO','Unexpectedly reached the end of the file whilst &
     &skipping lines.')
   endif ! ierr
  enddo ! i
 END SUBROUTINE skipio


 SUBROUTINE inve(v,inv)
!-----------------------!
! Inverts 3x3 matrices. !
!-----------------------!
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(in) :: v(3,3)
  DOUBLE PRECISION,INTENT(out) :: inv(3,3)
  DOUBLE PRECISION d
  d=v(1,1)*(v(2,2)*v(3,3)-v(2,3)*v(3,2))+ &
   &v(2,1)*(v(3,2)*v(1,3)-v(1,2)*v(3,3))+ &
   &v(3,1)*(v(1,2)*v(2,3)-v(1,3)*v(2,2))
  if(d==0.d0)call errstop('INVE','Trying to invert a singular determinant.')
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


 RECURSIVE SUBROUTINE adapt_simpson(f,a,b,s)
!--------------------------------------------------------------------------!
! This subroutine integrates a function f from a to b using Simpson's      !
! rule (adaptive).  The result is stored in s.                             !
!--------------------------------------------------------------------------!
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(in) :: a,b
  DOUBLE PRECISION,INTENT(out) :: s
  INTEGER :: j,it,jj
  DOUBLE PRECISION :: ps,pst,st,delta_x,sumf,tnm,x
  INTEGER,PARAMETER :: jmax=2000
  INTERFACE
   DOUBLE PRECISION FUNCTION f(x)
    DOUBLE PRECISION,INTENT(in) :: x
   END FUNCTION f
  END INTERFACE
  pst=0.5d0*huge(1.d0)
  ps=pst
  do j=1,jmax
   if(j==1)then
    st=0.5d0*(b-a)*(f(a)+f(b))
    it=1
   else
    tnm=dble(it)
    delta_x=(b-a)/tnm
    x=a+0.5d0*delta_x
    sumf=0.d0
    do jj=1,it
     sumf=sumf+f(x)
     x=x+delta_x
    enddo ! jj
    st=0.5d0*(st+(b-a)*sumf/tnm)
    it=it*2
   endif
   s=third*(4.d0*st-pst)
   if(j>5.and.(abs(s-ps)<eps_simp*abs(ps).or.(s==0.d0.and.ps==0.d0)))return
   ps=s
   pst=st
  enddo ! j
  call errstop('ADAPT_SIMPSON','Have not been able to converge integral in &
   &adapt_simpson.')
 END SUBROUTINE adapt_simpson


 DOUBLE PRECISION FUNCTION det_33(A)
!-----------------------------------------------------------------!
! This function returns the determinant of a real, 3x3 matrix, A. !
!-----------------------------------------------------------------!
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(in) :: A(3,3)
  det_33=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
   &+A(2,1)*(A(3,2)*A(1,3)-A(1,2)*A(3,3)) &
   &+A(3,1)*(A(1,2)*A(2,3)-A(1,3)*A(2,2))
 END FUNCTION det_33


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


 CHARACTER(1) FUNCTION l2s(l)
!-------------------------------------------------------------------------!
! Convert logical variable to a string of length 1.                       !
!-------------------------------------------------------------------------!
 IMPLICIT NONE
 LOGICAL,INTENT(in) :: l
 if(l)then
  l2s='T'
 else
  l2s='F'
 endif
 END FUNCTION l2s


 CHARACTER(72) FUNCTION write_mean(av,std_err_in_mean,err_prec_in)
!-----------------------------------------------------------------------------!
! Write out a mean value with the standard error in the mean in the form      !
! av(std_err_in_mean), e.g. 0.123546(7).  err_prec_in specifies the number of !
! digits of precision to which the error should be quoted (by default 1).     !
!-----------------------------------------------------------------------------!
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(in) :: av,std_err_in_mean
  INTEGER,INTENT(in),OPTIONAL :: err_prec_in
  INTEGER lowest_digit_to_quote,err_quote,err_prec,int_part,dec_part,i
  INTEGER,PARAMETER :: err_prec_default=1
  DOUBLE PRECISION av_quote
  CHARACTER(1) sgn
  CHARACTER(72) zero_pad

  if(std_err_in_mean<=0.d0)then
   write_mean='ERROR: NON-POSITIVE ERROR BAR!!!'
   return
  endif ! Error is negative

  if(present(err_prec_in))then
   if(err_prec_in>=1)then
    err_prec=err_prec_in
   else
    write_mean='ERROR: NON-POSITIVE PRECISION!!!'
    return
   endif ! err_prec_in sensible.
  else
   err_prec=err_prec_default
  endif ! Accuracy of error supplied.

! Work out lowest digit of precision that should be retained in the
! mean (i.e. the digit in terms of which the error is specified).
! Calculate the error in terms of this digit and round.
  lowest_digit_to_quote=floor(log(std_err_in_mean)/log(10.d0))+1-err_prec
  err_quote=nint(std_err_in_mean*10.d0**dble(-lowest_digit_to_quote))
  if(err_quote==10**err_prec)then
   lowest_digit_to_quote=lowest_digit_to_quote+1
   err_quote=err_quote/10
  endif ! err_quote rounds up to next figure.

  if(err_quote>=10**err_prec.or.err_quote<10**(err_prec-1))then
   write_mean='ERROR: BUG IN WRITE_MEAN!!!'
   return
  endif ! Check error is in range.

! Truncate the mean to the relevant precision.  Establish its sign,
! then take the absolute value and work out the integer part.
  av_quote=anint(av*10.d0**dble(-lowest_digit_to_quote)) &
   &*10.d0**dble(lowest_digit_to_quote)
  if(av_quote<0.d0)then
   sgn='-'
   av_quote=-av_quote
  else
   sgn=''
  endif ! Sign
  if(aint(av_quote)>dble(huge(1)))then
   write_mean='ERROR: NUMBERS ARE TOO LARGE IN WRITE_MEAN!'
   return
  endif ! Vast number
  int_part=floor(av_quote)

  if(lowest_digit_to_quote<0)then
! If the error is in a decimal place then construct string using
! integer part and decimal part, noting that the latter may need to
! be padded with zeros, e.g. if we want "0001" rather than "1".
   if(anint((av_quote-dble(int_part)) &
    &*10.d0**dble(-lowest_digit_to_quote))>dble(huge(1)))then
    write_mean='ERROR: NUMBERS ARE TOO LARGE IN WRITE_MEAN!'
    return
   endif ! Vast number
   dec_part=nint((av_quote-dble(int_part))*10.d0**dble(-lowest_digit_to_quote))
   zero_pad=''
   if(dec_part<0)then
    write_mean='ERROR: BUG IN WRITE_MEAN! (2)'
    return
   endif ! dec
   do i=1,-lowest_digit_to_quote-no_digits_int(dec_part)
    zero_pad(i:i)='0'
   enddo ! i
   write_mean=sgn//trim(i2s(int_part))//'.'//trim(zero_pad) &
    &//trim(i2s(dec_part))//'('//trim(i2s(err_quote))//')'
  else
! If the error is in a figure above the decimal point then, of
! course, we don't have to worry about a decimal part.
   write_mean=sgn//trim(i2s(int_part))//'(' &
    &//trim(i2s(err_quote*10**lowest_digit_to_quote))//')'
  endif ! lowest_digit_to_quote<0

 END FUNCTION write_mean


 INTEGER FUNCTION no_digits_int(i)
!----------------------------------------------------------------------!
! Calculate the number of digits in integer i.  For i>0 this should be !
! floor(log(i)/log(10))+1, but sometimes rounding errors cause this    !
! expression to give the wrong result.                                 !
!----------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: i
  INTEGER j,k
  j=i ; k=1
  do
   j=j/10
   if(j==0)exit
   k=k+1
  enddo
  no_digits_int=k
 END FUNCTION no_digits_int


 SUBROUTINE min_image(a,lat_vec,rec_vec,b,mag_rec_vec)
!----------------------------------------------------------------------------!
! This subroutine computes b as the minimum-image vector of vector a with    !
! respect to the lattice specified by lat_vec.  So -b is the vector from a   !
! to its closest lattice point.                                              !
! lat_vec holds the lattice vectors in columns, while rec_vec holds the      !
! reciprocal-lattice vectors without the factor of 2pi in columns.           !
! The vertex-checking algorithm is used, which could fail in some cases.     !
!----------------------------------------------------------------------------!
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(in) :: a(3),lat_vec(3,3),rec_vec(3,3),mag_rec_vec(3)
  DOUBLE PRECISION,INTENT(out) :: b(3)
  DOUBLE PRECISION Delta(3,8),mag_b_sq,mag_b,dist_sq(8),dist2,Delta1(3), &
   &Delta2(3),Delta3(3)
  INTEGER n(3),i,j,k,check_shell(3),kk(1)

! Work out which parallelepiped-shaped unit cell contains point a.
  n(1)=floor(dot_product(a,rec_vec(1:3,1)))
  n(2)=floor(dot_product(a,rec_vec(1:3,2)))
  n(3)=floor(dot_product(a,rec_vec(1:3,3)))

! Use the vertex-checking algorithm.
  Delta(1:3,1)=a-n(1)*lat_vec(1:3,1)-n(2)*lat_vec(1:3,2)-n(3)*lat_vec(1:3,3)
  dist_sq(1)=sum(Delta(1:3,1)**2)
  Delta(1:3,2)=Delta(1:3,1)-lat_vec(1:3,1)
  dist_sq(2)=sum(Delta(1:3,2)**2)
  Delta(1:3,3)=Delta(1:3,1)-lat_vec(1:3,2)
  dist_sq(3)=sum(Delta(1:3,3)**2)
  Delta(1:3,4)=Delta(1:3,2)-lat_vec(1:3,2)
  dist_sq(4)=sum(Delta(1:3,4)**2)
  Delta(1:3,5)=Delta(1:3,1)-lat_vec(1:3,3)
  dist_sq(5)=sum(Delta(1:3,5)**2)
  Delta(1:3,6)=Delta(1:3,5)-lat_vec(1:3,1)
  dist_sq(6)=sum(Delta(1:3,6)**2)
  Delta(1:3,7)=Delta(1:3,5)-lat_vec(1:3,2)
  dist_sq(7)=sum(Delta(1:3,7)**2)
  Delta(1:3,8)=Delta(1:3,7)-lat_vec(1:3,1)
  dist_sq(8)=sum(Delta(1:3,8)**2)
  kk=minloc(dist_sq(1:8))
  b(1:3)=Delta(1:3,kk(1))

! Calculate the number of additional shells that must be checked.
! All shells that lie within a sphere of the radius of the minimum-image
! distance found by the vertex-checking algorithm must be tested.
  mag_b_sq=sum(b(1:3)**2) ; mag_b=sqrt(mag_b_sq)
  check_shell(1:3)=floor(mag_b*mag_rec_vec(1:3))

! Check all the additional shells of lattice points (if necessary).
  if(any(check_shell>0))then
   do i=n(1)-check_shell(1),n(1)+check_shell(1)+1
    Delta1=a-i*lat_vec(1:3,1)
    do j=n(2)-check_shell(2),n(2)+check_shell(2)+1
     Delta2=Delta1-j*lat_vec(1:3,2)
     do k=n(3)-check_shell(3),n(3)+check_shell(3)+1
      Delta3=Delta2-k*lat_vec(1:3,3)
      dist2=sum(Delta3(1:3)**2)
      if(dist2<mag_b_sq)then
       mag_b_sq=dist2
       b=Delta3
      endif ! dist_sq<mag_b_sq
     enddo ! k
    enddo ! j
   enddo ! i
  endif ! check_shell>0

 END SUBROUTINE min_image


 SUBROUTINE errstop(sub,message)
!---------------------------!
! Report an error and stop. !
!---------------------------!
  IMPLICIT NONE
  CHARACTER(*),INTENT(in) :: sub,message
  write(*,*)
  write(*,*)'ERROR in subroutine '//trim(adjustl(sub))//'.'
  write(*,*)
  call wordwrap(trim(adjustl(message)))
  write(*,*)
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


END MODULE helpers

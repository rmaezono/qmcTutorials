!------------------------------------------------------------------------!
! This program enables the extrapolation of DMC energies to zero         !
! time step by fitting a polynomial form to the DMC energy as a          !
! function of time step, i.e. E = k1*tau^n1 + k2*tau^n2 + ...            !
! Note that the exponents n1, n2, ... do not need to be integers.        !
! One of them should be zero, however.                                   !
!                                                                        !
! The program will ask the user for the name of the file containing the  !
! energy against time step data.  The file should have either two or     !
! three columns holding, respectively, the time step, the DMC energy and !
! the statistical error bar in the DMC energy.                           !
!                                                                        !
! NDD, 21.01.02                                                          !
!                                                                        !
! Changes                                                                !
! =======                                                                !
! NDD  09.2006  Extensive tidying.  Simplified format of input file.     !
!               Added printout of fitting function.                      !
! NDD  01.2009  Simplified calculation of error bars.                    !
!------------------------------------------------------------------------!


MODULE utils
!----------------------------!
! Miscellaneous subroutines. !
!----------------------------!
 IMPLICIT NONE
 PRIVATE
 PUBLIC calc_parameters,chi_squared,get_exponents,get_file,check_file, &
  &read_file,i2s,construct_graphstr,construct_poly,write_mean,no_points, &
  &no_terms,tau_point,e_point,n_exp,rec_sigma_sq,dp
! Double precision type
 INTEGER,PARAMETER :: dp=kind(1.d0)
! Number of (tau,E) data points
 INTEGER no_points,no_terms
! Time steps, energies, exponents in fitting function and errors in energies.
 REAL(dp),ALLOCATABLE :: tau_point(:),e_point(:),n_exp(:),rec_sigma_sq(:)


CONTAINS


 SUBROUTINE calc_parameters(a,sigma_a)
!-------------------------------------------------------------!
! Evaluate the parameters using the method given in my notes. !
!-------------------------------------------------------------!
  IMPLICIT NONE
  REAL(dp),INTENT(out) :: a(no_terms)
! Following should just be INTENT(out), but NAG compiler chokes.
  REAL(dp),INTENT(inout),OPTIONAL :: sigma_a(no_terms)
  INTEGER info,i,j,ipiv(no_terms),lwork,ialloc
  REAL(dp) M(no_terms,no_terms),Minv(no_terms,no_terms),c(no_terms),tempr(1)
  REAL(dp),ALLOCATABLE :: work(:)
  LOGICAL errorbars_present
  INTERFACE
   SUBROUTINE dsytrf(UPLO,N,A,LDA,IPIV,WORK,LWORK,INFO)
    CHARACTER(1),INTENT(in) :: UPLO
    INTEGER,INTENT(in) :: N,LDA,LWORK
    REAL(kind(1.d0)),INTENT(out) :: A(LDA,*),WORK(*)
    INTEGER,INTENT(out) :: IPIV(*)
    INTEGER,INTENT(out) :: INFO
   END SUBROUTINE dsytrf
   SUBROUTINE dsytri(UPLO,N,A,LDA,IPIV,WORK,INFO)
    CHARACTER(1),INTENT(in) :: UPLO
    INTEGER,INTENT(in) :: N,LDA,IPIV(*)
    REAL(kind(1.d0)),INTENT(inout) :: A(LDA,*)
    REAL(kind(1.d0)),INTENT(out) :: WORK(*)
    INTEGER,INTENT(out) :: INFO
   END SUBROUTINE dsytri
  END INTERFACE

  errorbars_present=present(sigma_a)

! Construct vector c.
  call construct_c(c,errorbars_present)

! Construct matrix M.
  call construct_M(M,errorbars_present)

! Invert matrix M.
  Minv=M
  call dsytrf('L',no_terms,Minv(1,1),no_terms,ipiv(1),tempr(1),-1,info)
  if(info/=0)then
   write(6,*)'Matrix inversion failed (1).'
   stop
  endif ! info/=0
  lwork=nint(tempr(1))
  allocate(work(lwork),stat=ialloc)
  if(ialloc/=0)then
   write(6,*)'Allocation error: WORK (1).'
   stop
  endif ! ialloc/=0
  call dsytrf('L',no_terms,Minv(1,1),no_terms,ipiv(1),work(1),lwork,info)
  if(info/=0)then
   write(6,*)'Matrix inversion failed (2).'
   stop
  endif ! info/=0
  deallocate(work)
  allocate(work(no_terms),stat=ialloc)
  if(ialloc/=0)then
   write(6,*)'Allocation error: WORK (2).'
   stop
  endif ! ialloc/-0
  call dsytri('L',no_terms,Minv(1,1),no_terms,ipiv(1),work(1),info)
  if(info/=0)then
   write(6,*)'Matrix inversion failed (3).'
   stop
  endif ! info/=0
  deallocate(work)

! Complete the upper triangle of Minv.
  do i=1,no_terms-1
   do j=i+1,no_terms
    Minv(i,j)=Minv(j,i)
   enddo ! j
  enddo ! i

! Hence evaluate the coefficients of the terms in the polynomial.
  a=matmul(Minv,c)

! Evaluate the standard errors in the coefficients.
  if(errorbars_present)then
   do j=1,no_terms
    sigma_a(j)=sqrt(Minv(j,j))
   enddo ! j
  endif ! errorbars_present

 END SUBROUTINE calc_parameters


 REAL(dp) FUNCTION chi_squared(a,errorbars_present)
!--------------------------------------------------------------------------!
! Evaluate the chi-squared value of the fit.  If error bars are not given, !
! return the least-squares function.                                       !
!--------------------------------------------------------------------------!
  IMPLICIT NONE
  REAL(dp),INTENT(in) :: a(no_terms)
  LOGICAL,INTENT(in) :: errorbars_present
  INTEGER i,k
  REAL(dp) e_fit
  chi_squared=0.d0
  do i=1,no_points
   e_fit=0.d0
   do k=1,no_terms
    e_fit=e_fit+a(k)*tau_point(i)**n_exp(k)
   enddo ! k
   if(errorbars_present)then
    chi_squared=chi_squared+(e_point(i)-e_fit)**2*rec_sigma_sq(i)
   else
    chi_squared=chi_squared+(e_point(i)-e_fit)**2
   endif ! errorbars_present
  enddo ! i
 END FUNCTION chi_squared


 SUBROUTINE construct_c(c,errorbars_present)
!----------------------------------------!
! Construct the vector c (see my notes). !
!----------------------------------------!
  IMPLICIT NONE
  REAL(dp),INTENT(out) :: c(no_terms)
  LOGICAL,INTENT(in) :: errorbars_present
  INTEGER i,j
  if(errorbars_present)then
   do j=1,no_terms
    c(j)=0.d0
    do i=1,no_points
     c(j)=c(j)+e_point(i)*tau_point(i)**n_exp(j)*rec_sigma_sq(i)
    enddo ! i
   enddo ! j
  else
   do j=1,no_terms
    c(j)=0.d0
    do i=1,no_points
     c(j)=c(j)+e_point(i)*tau_point(i)**n_exp(j)
    enddo ! i
   enddo ! j
  endif ! errorbars_present
 END SUBROUTINE construct_c


 SUBROUTINE construct_M(M,errorbars_present)
!----------------------------------------!
! Construct the matrix M (see my notes). !
! Only need lower triangular part.       !
!----------------------------------------!
  IMPLICIT NONE
  REAL(dp),INTENT(out) :: M(no_terms,no_terms)
  LOGICAL,INTENT(in) :: errorbars_present
  INTEGER i,j,k
  M=0.d0
  if(errorbars_present)then
   do k=1,no_terms
    do j=k,no_terms
     do i=1,no_points
      M(j,k)=M(j,k)+tau_point(i)**(n_exp(j)+n_exp(k))*rec_sigma_sq(i)
     enddo ! i
    enddo ! j
   enddo ! k
  else
   do k=1,no_terms
    do j=k,no_terms
     do i=1,no_points
      M(j,k)=M(j,k)+tau_point(i)**(n_exp(j)+n_exp(k))
     enddo ! i
    enddo ! j
   enddo ! k
  endif ! errorbars_present
 END SUBROUTINE construct_M


 SUBROUTINE get_exponents
!-----------------------------------------------------------------------!
! Ask the user to enter the exponents of the terms in the fitting poly. !
!-----------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER i,j,ierr,ialloc
  REAL(dp) temp
  REAL(dp),PARAMETER :: tol_zero=1.d-8
  LOGICAL test_flag

  do
   write(6,*)'Please enter the number of terms in the polynomial to be fitted.'
   read(5,*,iostat=ierr)no_terms
   if(ierr/=0)then
    write(6,*)'Error - please try again.'
   else
    if(no_terms<1)then
     write(6,*)'Interpolating polynomial must have at least one term!'
    else
     if(no_terms>no_points)then
      write(6,*)'Insufficient data to determine polynomial coefficients.'
     else
      exit
     endif ! no_terms>no_points
    endif ! no_terms<1
   endif ! ierr/=0
  enddo

  allocate(n_exp(no_terms),stat=ialloc)
  if(ialloc/=0)then
   write(6,*)'Allocation error.'
   stop
  endif ! ialloc/=0

  do
   write(6,*)'Please enter the '//trim(i2s(no_terms)) &
    &//' exponents in the fitting polynomial.'
   read(5,*,iostat=ierr)n_exp(1:no_terms)
   if(ierr/=0)then
    write(6,*)'Error - please try again.'
   else
    if(any(n_exp<0))then
     write(6,*)'Found negative exponent; please try again.'
    else
     test_flag=.false.
ol:  do i=1,no_terms
      do j=1,i-1
       if(abs(n_exp(i)-n_exp(j))<=2.d0*tol_zero)then
        test_flag=.true.
        exit ol
       endif ! Identical exponents
      enddo ! j
     enddo ol ! i
     if(test_flag)then
      write(6,*)'Two exponents appear to be identical.  Please try again.'
     else
      if(all(abs(n_exp)>=tol_zero))then
       write(6,*)'No constant term in fitting polynomial.'
      else
       exit
      endif ! No constant term.
     endif ! Identical exponents
    endif ! Negative exponents
   endif ! Error
  enddo

! Sort exponents into ascending order
  do i=1,no_terms-1
   do j=i+1,no_terms
    if(n_exp(j)<n_exp(i))then
     temp=n_exp(i)
     n_exp(i)=n_exp(j)
     n_exp(j)=temp
    endif ! Swap needed
   enddo ! j
  enddo ! i

! n_exp(1) is within tol_zero of zero; make it rigorously zero.
  n_exp(1)=0.d0

 END SUBROUTINE get_exponents


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


 SUBROUTINE get_file(in_file)
!------------------------!
! Find out the filename. !
!------------------------!
  IMPLICIT NONE
  CHARACTER(*),INTENT(out) :: in_file
  INTEGER ierr
  LOGICAL file_exists
  do
   write(6,*)'Please enter name of data file.'
   read(5,'(a)',iostat=ierr)in_file
   if(ierr/=0)then
    write(6,*)'Error - please try again.'
   else
    in_file=adjustl(in_file)
    inquire(file=trim(in_file),exist=file_exists)
    if(.not.file_exists)then
     write(6,*)'File does not appear to exist. Please try again.'
    else
     exit
    endif ! File nonexistent
   endif ! ierr/=0
  enddo
  write(6,*)
 END SUBROUTINE get_file


 SUBROUTINE check_file(io,in_file,no_lines,errorbars_present)
!------------------------------------!
! Count the lines in the input file. !
!------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: io
  CHARACTER(*),INTENT(in) :: in_file
  INTEGER,INTENT(out) :: no_lines
  LOGICAL,INTENT(out) :: errorbars_present
  INTEGER ierr,ipos,phys_lines
  REAL(dp) dummy(3)
  CHARACTER(1024) line
  errorbars_present=.true.
  no_lines=0 ; phys_lines=0
! Look for first line
  do
   read(io,'(a)',iostat=ierr)line
   if(ierr<0)exit
   phys_lines=phys_lines+1
   if(ierr>0)then
    write(6,*)'Error reading '//trim(in_file)//' at line ' &
     &//trim(i2s(phys_lines))//'.'
    stop
   endif
   line=adjustl(line)
   ipos=scan(line,'#!')
   if(ipos==1)cycle
   if(ipos>1)line=line(1:ipos-1)
   if(len_trim(line)==0)cycle
   if(errorbars_present)then
    read(line,*,iostat=ierr)dummy(1:3)
    if(ierr/=0)errorbars_present=.false.
   endif
   if(.not.errorbars_present)then
    read(line,*,iostat=ierr)dummy(1:2)
    if(ierr/=0)then
     write(6,*)'Error reading '//trim(in_file)//' at line '//&
      &trim(i2s(phys_lines))//'.'
     stop
    endif
   endif
   no_lines=no_lines+1
  enddo
  rewind(io)
  if(no_lines==0)then
   write(6,*)'File '//trim(in_file)//' doesn''t contain any data.'
   stop
  endif
 END SUBROUTINE check_file


 SUBROUTINE read_file(io,in_file,errorbars_present)
!-------------------------------------!
! Read in the data in the input file. !
!-------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: io
  CHARACTER(*),INTENT(in) :: in_file
  LOGICAL,INTENT(in) :: errorbars_present
  INTEGER i,ierr,ipos,phys_lines
  CHARACTER(1024) line
  do i=1,no_points
   phys_lines=0
   do
    read(io,'(a)',iostat=ierr)line
    if(ierr<0)then
     write(6,*)'Error reading '//trim(in_file)//' at line '//&
      &trim(i2s(phys_lines))//': unexpected end of file.'
     stop
    endif
    phys_lines=phys_lines+1
    if(ierr>0)then
     write(6,*)'Error reading '//trim(in_file)//' at line '//&
      &trim(i2s(phys_lines))//'.'
     stop
    endif
    line=adjustl(line)
    ipos=scan(line,'#!')
    if(ipos==1)cycle
    if(ipos>1)line=line(1:ipos-1)
    if(len_trim(line)==0)cycle
    exit
   enddo
   if(errorbars_present)then
    read(line,*,iostat=ierr)tau_point(i),e_point(i),rec_sigma_sq(i)
   else
    read(line,*,iostat=ierr)tau_point(i),e_point(i)
   endif
   if(ierr/=0)then
    write(6,*)'Error reading '//trim(in_file)//'.'
    stop
   endif ! ierr/=0
   if(errorbars_present)then
    if(rec_sigma_sq(i)<=0.d0)then
     write(*,*)'Non-positive energy error bar at data line '//trim(i2s(i))//'.'
     stop
    endif
    rec_sigma_sq(i)=1.d0/rec_sigma_sq(i)**2
   endif ! errorbars_present
  enddo ! i
 END SUBROUTINE read_file


 SUBROUTINE construct_graphstr(a,graphstr)
!------------------------------------------------------------------------!
! This subroutine returns the fitted polynomial in a suitable format for !
! pasting into xmgrace.                                                  !
!------------------------------------------------------------------------!
  IMPLICIT NONE
  REAL(dp),INTENT(in) :: a(no_terms)
  CHARACTER(*),INTENT(out) :: graphstr
  INTEGER j
  REAL(dp),PARAMETER :: tol_zero=1.d-8
  CHARACTER(1) plusstr
  CHARACTER(30) coeffstr
  CHARACTER(36) pwstr
  graphstr='y='
  do j=1,no_terms
   if(abs(anint(n_exp(j))-n_exp(j))<tol_zero)then
    if(nint(n_exp(j))==0)then
     pwstr=''
    elseif(nint(n_exp(j))==1)then
     pwstr='*x'
    else
     pwstr='*x^'//trim(i2s(nint(n_exp(j))))
    endif ! n_exp=0
   else
    write(pwstr,*)n_exp(j)
    pwstr='*x^'//trim(adjustl(pwstr))
   endif ! Power is an integer
   if(j>1.and.a(j)>=0.d0)then
    plusstr='+'
   else
    plusstr=''
   endif ! "+" needed
   write(coeffstr,*)a(j)
   coeffstr=adjustl(coeffstr)
   graphstr=trim(graphstr)//trim(plusstr)//trim(coeffstr)//trim(pwstr)
  enddo ! j
! Change "d" to "e" in graphstr.
  do j=1,len_trim(graphstr)
   if(graphstr(j:j)=='d'.or.graphstr(j:j)=='D')graphstr(j:j)='E'
  enddo ! j
 END SUBROUTINE construct_graphstr


 SUBROUTINE construct_poly(polystr)
!------------------------------------------------------------------------!
! This subroutine returns the fitted polynomial in a suitable format for !
! pasting into xmgrace.                                                  !
!------------------------------------------------------------------------!
  IMPLICIT NONE
  CHARACTER(*),INTENT(out) :: polystr
  INTEGER j
  REAL(dp),PARAMETER :: tol_zero=1.d-8
  CHARACTER(36) pwstr
  polystr='E ='
  do j=1,no_terms
   if(abs(anint(n_exp(j))-n_exp(j))<tol_zero)then
    if(nint(n_exp(j))==0)then
     pwstr=''
    elseif(nint(n_exp(j))==1)then
     pwstr='*tau'
    else
     pwstr='*tau^'//trim(i2s(nint(n_exp(j))))
    endif ! n_exp=0
   else
    write(pwstr,*)n_exp(j)
    pwstr='*tau^'//trim(adjustl(pwstr))
   endif ! Power is an integer
   if(j>1)then
    polystr=trim(polystr)//' + k_'//trim(i2s(j))//trim(pwstr)
   else
    polystr=trim(polystr)//' k_'//trim(i2s(j))//trim(pwstr)
   endif ! "+" needed
  enddo ! j
 END SUBROUTINE construct_poly


 CHARACTER(72) FUNCTION write_mean(av,std_err_in_mean,err_prec_in)
!--------------------------------------------------------------------!
! Write out a mean value with the standard error in the mean in the  !
! form av(std_err_in_mean), e.g. 0.123546(7).  err_prec_in specifies !
! the number of digits of precision to which the error should be     !
! quoted (by default 1).                                             !
!--------------------------------------------------------------------!
  IMPLICIT NONE
  REAL(dp),INTENT(in) :: av,std_err_in_mean
  INTEGER,INTENT(in),OPTIONAL :: err_prec_in
  INTEGER lowest_digit_to_quote,err_quote,err_prec,int_part,dec_part,i
  INTEGER,PARAMETER :: err_prec_default=1
  REAL(dp) av_quote
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
!------------------------------------------------------------------!
! Calculate the number of digits in integer i.  For i>0 this       !
! should be floor(log(i)/log(10))+1, but sometimes rounding errors !
! cause this expression to give the wrong result.                  !
!------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: i
  INTEGER j,k
  j=i
  k=1
  do
   j=j/10
   if(j==0)exit
   k=k+1
  enddo
  no_digits_int=k
 END FUNCTION no_digits_int


END MODULE utils


PROGRAM extrapolate_tau
!---------------------------!
! Main program starts here. !
!---------------------------!
 USE utils
 IMPLICIT NONE
 INTEGER j,ierr,ialloc,io
 REAL(dp),ALLOCATABLE :: a(:),sigma_a(:)
 CHARACTER(3) indstr
 CHARACTER(80) in_file
 CHARACTER(320) tempstr
 LOGICAL errorbars_present

 write(6,*)
 write(6,*)'DMC time-step extrapolator'
 write(6,*)'=========================='
 write(6,*)

! Read in (tau,energy) data from file specified by user.
 call get_file(in_file)
 io=8
 write(6,*)'Reading from "'//trim(in_file)//'".'
 open(unit=io,file=trim(in_file),status='old',iostat=ierr)
 if(ierr/=0)then
  write(6,*)'Error opening "'//trim(in_file)//'".'
  stop
 endif ! ierr
 call check_file(io,in_file,no_points,errorbars_present)
 write(6,*)'Number of data points supplied: '//trim(i2s(no_points))
 if(no_points<1)then
  write(6,*)'More DMC energy points needed.'
  stop
 endif
 if(errorbars_present)then
  write(6,*)'Error bars on energy data are supplied.'
  allocate(tau_point(no_points),e_point(no_points),rec_sigma_sq(no_points), &
   &stat=ialloc)
 else
  write(6,*)'No error bars are supplied for energy data.'
  allocate(tau_point(no_points),e_point(no_points),stat=ialloc)
 endif ! errorbars_present
 if(ialloc/=0)then
  write(6,*)'Allocation error.'
  stop
 endif ! ialloc/=0
 call read_file(io,in_file,errorbars_present)
 if(any(tau_point<=0.d0))write(6,*) &
  &'Warning: negative time step found in input file.'
 write(6,*)
 close(io)

! Ask user to enter exponents for fitting polynomial.
 call get_exponents
 if(errorbars_present)then
  allocate(a(no_terms),sigma_a(no_terms),stat=ialloc)
 else
  allocate(a(no_terms),stat=ialloc)
 endif ! errorbars_present
 if(ialloc/=0)then
  write(6,*)'Allocation error (2).'
  stop
 endif ! ialloc/=0
 call construct_poly(tempstr)
 write(6,*)'Fitting polynomial is of form:'
 write(6,*)trim(tempstr)
 write(6,*)

! We solve Ma=c for vector of coefficients a (see my notes).
 if(errorbars_present)then
  write(6,*)'Performing chi-squared fit...'
  call calc_parameters(a,sigma_a)
 else
  write(6,*)'Performing least-squares fit...'
  call calc_parameters(a)
 endif ! errorbars_present
 write(6,*)'Done.'
 write(6,*)

! Evaluate chi-squared.
 if(errorbars_present)then
  write(6,*)'Chi-squared value: ',chi_squared(a,errorbars_present)
 else
  write(6,*)'Least-squares function: ',chi_squared(a,errorbars_present)
 endif ! errorbars_present
 write(6,*)

 if(errorbars_present)then
  do j=1,no_terms
   indstr=adjustl(trim(i2s(j)))
   write(6,'(" ",a," = ",es24.16," +/- ",es24.16)')'k_'//indstr,a(j),sigma_a(j)
  enddo ! j
 else
  do j=1,no_terms
   indstr=adjustl(trim(i2s(j)))
   write(6,'(" ",a," = ",es24.16)')'k_'//indstr,a(j)
  enddo ! j
 endif ! errorbars_present
 write(6,*)

! Evaluate the energy extrapolated to time step of zero as the
! coefficient of tau**0.0 (NB, exponents have been sorted).
 if(errorbars_present)then
  write(6,*)'DMC energy at zero time step: '//trim(write_mean(a(1),sigma_a(1)))
 else
  write(6,*)'DMC energy at zero time step: ',a(1)
 endif ! errorbars_present
 write(6,*)

! Write out fitted polynomial in a form that can be pasted into xmgrace.
 call construct_graphstr(a,tempstr)
 write(6,*)'Fitted polynomial in y(x) form (for pasting into XMGrace):'
 write(6,*)trim(tempstr)
 write(6,*)

 deallocate(tau_point,e_point,n_exp,a)
 if(errorbars_present)deallocate(rec_sigma_sq,sigma_a)

 write(6,*)'Program finished.'
 write(6,*)

END PROGRAM extrapolate_tau

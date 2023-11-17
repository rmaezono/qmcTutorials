!------------------------------------------------------------------------------!
! TAHELPER        EM/NDD  01.2014                                              !
!                                                                              !
! This program fits a function to a set of data                                !
! The set of data are "i,y(i),error in y(i),z(i)" in E_v_twist.dat file.       !
! In E_v_twist.dat file, i is the number of the twist.                         !
! The fitting function is of "y(i)=a0+a1.z(i)" form.                           !
! Here y(i) is the DMC energy for each twist.                                  !
!      a0 is final twist-averaged energy.                                      !
! This program is not intended to be called directly by users.  It is called by!
! the twistanalysis_castep script.                                             !
!------------------------------------------------------------------------------!

MODULE utils
!----------------------------!
! Miscellaneous subroutines. !
!----------------------------!
 IMPLICIT NONE
 INTEGER,PARAMETER :: dp=KIND(1.d0)
 INTEGER nlines
 REAL(dp),ALLOCATABLE :: rec_sigma_y_sq(:),z(:),y(:)
 REAL(dp) twe_fit,err_fit,chi2
 PRIVATE
 PUBLIC rec_sigma_y_sq,z,y,nlines,twe_fit,err_fit,chi2,read_file,fit_parameter


CONTAINS


 SUBROUTINE read_file
!-------------------------------------!
! Reads E_v_twist.dat file. file.     !
!-------------------------------------!
  IMPLICIT NONE
  INTEGER i,j,ierr
  CHARACTER(50) char50
  open(unit=8,file='E_v_twist.dat',status='old',iostat=ierr)
  if(ierr/=0)then
   write(6,*)'Cannot open E_v_twist.dat file.'
   stop
  endif ! ierr
! Read number of data lines in E_v_twist.dat file
  nlines=0
  do 
   read(8,*,iostat=ierr)char50
   if(ierr>0)write(6,*)'Problem reading E_v_twist.dat file.'
   if(ierr<0)exit
   nlines=nlines+1
  enddo
  rewind(8)
! Read all data set value in E_v_twist.dat file 
  allocate(y(nlines),rec_sigma_y_sq(nlines),z(nlines))  
  do i=1,nlines
   read(8,*,iostat=ierr)j,y(i),rec_sigma_y_sq(i),z(i) 
   if(j/=i)then
    write(*,*)'Unexpected data order in E_v_twist.dat file.'
    stop
   endif ! j/=i
   if(rec_sigma_y_sq(i)<=0)then
    write(*,*)'Non-positive error bar in E_v_twist.dat file.'
    stop
   endif
   rec_sigma_y_sq(i)=1.d0/rec_sigma_y_sq(i)**2
  enddo ! i
  close(8)
 END SUBROUTINE read_file


 SUBROUTINE fit_parameter
!-------------------------------------------------------------!
! Evaluates matrix components Ma=C.                           ! 
! Performs least-squares fit. See NDD's notes.                !
!-------------------------------------------------------------!
  IMPLICIT NONE  
  REAL(dp),PARAMETER :: tol=1.d-8
  REAL(dp) det_M
  REAL(dp) M(2,2),M_inv(2,2),c(2)
  REAL(dp) a(2)
  INTEGER i,j,k
! Define 2 by 2 matrix M(j,k).
  do k=1,2
   do j=1,2
    M(j,k)=0.d0
    do i=1,nlines
     M(j,k)=M(j,k)+z(i)**(j+k-2)*rec_sigma_y_sq(i)
    enddo ! i
   enddo ! j
  enddo ! k
! Define determinant M(j,k).
  det_M=M(1,1)*M(2,2)-M(2,1)*M(1,2)
! Find the inverse of M(j,k).
  if(abs(det_M)<tol*(sum(abs(M))*0.25d0))then
   write(*,*)'Matrix M is singular.'
   stop
  endif
  M_inv(1,1)=M(2,2)/det_M
  M_inv(1,2)=-M(1,2)/det_M
  M_inv(2,1)=M_inv(1,2)
  M_inv(2,2)=M(1,1)/det_M
! Define vector C(j)
  do j=1,2
   c(j)=0.d0
   do i=1,nlines
    c(j)=c(j)+y(i)*z(i)**(j-1)*rec_sigma_y_sq(i)
   enddo ! i
  enddo ! j
! Evaluate the coefficients of the terms in the function.
  a=matmul(M_inv,c)
  twe_fit=a(1)  
  err_fit=SQRT(M_inv(1,1))
  chi2=chi_squared(a)
 END SUBROUTINE fit_parameter


 REAL(dp) FUNCTION chi_squared(n)
! Evaluate the chi-squared value of the fit.
  IMPLICIT NONE
  REAL(dp),INTENT(in) :: n(2)
  REAL(dp) y_fit
  INTEGER i,k
  chi_squared=0.d0
  do i=1,nlines
   y_fit=0.d0
   do k=1,2
    y_fit=y_fit+n(k)*z(i)**(k-1)
   enddo ! k
   chi_squared=chi_squared+(y(i)-y_fit)**2*rec_sigma_y_sq(i)
  enddo ! i
 END FUNCTION chi_squared


END MODULE utils


PROGRAM tahelper
!---------------------------!
! Main program starts here. !
!---------------------------!
 USE utils
 IMPLICIT NONE
! Reads in the existing E_v_twist.dat file.
 call read_file
! Reports fitted parameter and the error bar.
 call fit_parameter
 write(6,*)'Final twist-averaged energy   : ',twe_fit
 write(6,*)'Error in twist-averaged energy: ',err_fit
 write(6,*)'Chi-squared value             : ',chi2
 deallocate (y,rec_sigma_y_sq,z)
END PROGRAM tahelper

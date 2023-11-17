 MODULE numerical
!---------------------------------------------------------------------!
! Module containing standard numerical functions.                     !
!                                                                     !
! Currently contains:                                                 !
!  - straight_line_fit                                                !
!  - interp_nev (interpolation with Neville's algorithm)              !
!---------------------------------------------------------------------!
 USE dsp
 IMPLICIT NONE
 PRIVATE
 PUBLIC straight_line_fit,interp_nev

 CONTAINS


  SUBROUTINE straight_line_fit(x,y,n,a,b,siga,sigb,chi2)
  IMPLICIT NONE
  INTEGER,INTENT(in) :: n
  REAL(dp),INTENT(in) :: x(n),y(n)
  REAL(dp),INTENT(out) :: a,b,siga,sigb,chi2
  REAL(dp) sigdat,ss,sx,sxoss,sy,st2
  REAL(dp) :: t(n)

  ss=real(n,dp)
  sx=sum(x)
  sy=sum(y)

  sxoss=sx/ss
  t(:)=x(:)-sxoss

  b=dot_product(t,y)

  st2=dot_product(t,t)
  b=b/st2
  a=(sy-sx*b)/ss
  siga=sqrt((1.d0+sx*sx/(ss*st2))/ss)
  sigb=sqrt(1.d0/st2)
  t(:)=y(:)-a-b*x(:)

  chi2=dot_product(t,t)
  sigdat=sqrt(chi2/real(n-2,dp))
  siga=siga*sigdat
  sigb=sigb*sigdat

  END SUBROUTINE straight_line_fit


  SUBROUTINE interp_nev(points,val,npoints,x,val_interp,err_interp)
!-------------------------------------------------------------!
! INTERP_NEV: an implementation of Neville's algorithm which  !
! constructs the unique interpolating polynomial of degree    !
! n-1 through n points. The technique is based on the Newton  !
! form of the interpolating polynomial and the recursion      !
! relation for the divided differences.                       !
!                                                             !
! Reference: http://en.wikipedia.org/wiki/Neville's_algorithm !
!                                                             !
! Supplied with (1) NPOINTS function values VAL on a set of   !
! POINTS and (2) an interpolation point X, then the routine   !
! returns an interpolated value VAL_INTERP with an estimated  !
! error ERR_INTERP.                                           !
!-------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: npoints
  REAL(dp),INTENT(in) :: points(npoints),val(npoints),x
  REAL(dp),INTENT(out) :: val_interp,err_interp
  INTEGER i,j,n
  INTEGER,PARAMETER :: npoly=10 ! npoints>npoly deliberately not checked.
  REAL(dp) c(npoly),d(npoly),sep1,sep2,sep3,cmd

  sep1=abs(x-points(1))
  n=1
  do i=1,npoints
   sep2=abs(x-points(i))
   if(sep2<sep1)then
    sep1=sep2
    n=i
   endif
   c(i)=val(i) ; d(i)=val(i)
  enddo
! n is now the index of the grid point closest to the interpolation point x
  val_interp=val(n)

! Perform Neville's algorithm
  n=n-1
  do i=1,npoints-1 ! loop over tableau columns
   do j=1,npoints-i
    sep1=points(j)-x ; sep2=points(i+j)-x
    cmd=c(j+1)-d(j)
! c and d are the differences between parents and daughters in the tableau
! as we move from left to right.
    sep3=sep1-sep2
    c(j)=sep1*cmd/sep3 ; d(j)=sep2*cmd/sep3
   enddo
   if(n*2<npoints-i)then
    err_interp=c(n+1)
   else
    err_interp=d(n)
    n=n-1
   endif
   val_interp=val_interp+err_interp ! interpolated value
  enddo

  END SUBROUTINE interp_nev


 END MODULE numerical

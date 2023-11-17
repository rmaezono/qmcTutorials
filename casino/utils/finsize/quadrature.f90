 SUBROUTINE gauss_quadrature(x1,x2,n,func,integral)
!--------------------------------------------------------------------!
! Gauss quadrature of function func between points a and b using     !
! n points with abscissae x and weights w.                           !
!--------------------------------------------------------------------!
  USE dsp
  USE run_control, ONLY : errstop
  INTEGER,INTENT(in) :: n
  REAL(dp),INTENT(in) :: x1,x2
  REAL(dp),INTENT(out) :: integral
  INTEGER j,ialloc
  REAL(dp),ALLOCATABLE :: x(:),w(:)
  REAL(dp),PARAMETER :: pi=3.14159265358979324d0
  INTERFACE
   FUNCTION func(x)
    USE dsp,ONLY : dp
    REAL(dp) :: func
    REAL(dp),INTENT(in) :: x
   END FUNCTION func
  END INTERFACE

  allocate(x(n),w(n),stat=ialloc)
  if(ialloc/=0)call errstop('GAUSS_QUADRATURE','Allocation problem.')

! Get abscissae and weights for Gauss-Legendre (adapt for other schemes)
  call gauss_legendre_xw

  integral=0.d0
  do j=1,n
   integral=integral+w(j)*func(x(j))
  enddo


 CONTAINS


  SUBROUTINE gauss_legendre_xw
!--------------------------------------------------------------------------!
! Calculate abscissae x and weights w for Gauss-Legendre quadrature.       !
!--------------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER its,j,m
  INTEGER,PARAMETER :: maxit=10
  REAL(dp) xl,xm
  REAL(dp),PARAMETER :: eps=1.d-14
  REAL(dp),DIMENSION((size(x)+1)/2) :: p1,p2,p3,pp,z,z1
  LOGICAL,DIMENSION((size(x)+1)/2) :: unfinished

  if(size(x)/=size(w))call errstop('GAULEG','Number of abscissa not equal to &
   &number of weights.')

  m=(n+1)/2
  xm=0.5d0*(x2+x1)
  xl=0.5d0*(x2-x1)
  z=cos(pi*(arith_prog(1,1,m)-0.25d0)/(n+0.5d0))
  unfinished=.true.

  do its=1,maxit
   where(unfinished)
    p1=1.d0
    p2=0.d0
   end where
   do j=1,n
    where(unfinished)
     p3=p2
     p2=p1
     p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
    end where
   enddo
   where(unfinished)
    pp=n*(z*p1-p2)/(z*z-1.d0)
    z1=z
    z=z1-p1/pp
    unfinished=(abs(z-z1)>eps)
   end where
   if(.not.any(unfinished))exit
  enddo

  if(its==maxit+1)call errstop('GAULEG','Too many iterations.')
  x(1:m)=xm-xl*z
  x(n:n-m+1:-1)=xm+xl*z
  w(1:m)=2.d0*xl/((1.d0-z**2)*pp**2)
  w(n:n-m+1:-1)=w(1:m)

  END SUBROUTINE gauss_legendre_xw


  FUNCTION arith_prog(first,increment,n)
!-----------------------------------------------------!
! Array function returning an arithmetic progression  !
!-----------------------------------------------------!
  INTEGER,INTENT(in) :: first,increment,n
  INTEGER k,k2,temp
  INTEGER,PARAMETER :: npar_arth=16,npar2_arth=8
  INTEGER,DIMENSION(n) :: arith_prog
  if(n>0)arith_prog(1)=first
  if(n<=npar_arth)then
   do k=2,n
    arith_prog(k)=arith_prog(k-1)+increment
   enddo
  else
   do k=2,npar2_arth
    arith_prog(k)=arith_prog(k-1)+increment
   enddo
   temp=increment*npar2_arth
   k=npar2_arth
   do
    if(k>=n)exit
    k2=k+k
    arith_prog(k+1:min(k2,n))=temp+arith_prog(1:min(k,n-k))
    temp=temp+temp
    k=k2
   enddo
  endif
  END FUNCTION arith_prog


 END SUBROUTINE gauss_quadrature

MODULE runge_kutta
 USE dsp
 USE run_control, ONLY : errstop_master
 USE store, ONLY : dim
 IMPLICIT NONE
 PRIVATE
 PUBLIC stepper_rk


CONTAINS


 SUBROUTINE stepper_rk(x,vel,t,htry,eps,xscal,hdid,hnext,velocity)
 IMPLICIT NONE
 REAL(dp),INTENT(inout) :: x(dim)
 REAL(dp),INTENT(inout) :: t
 REAL(dp),INTENT(in) :: vel(dim),xscal(dim)
 REAL(dp),INTENT(in) :: htry,eps
 REAL(dp),INTENT(out) :: hdid,hnext

 INTERFACE
  SUBROUTINE velocity(t,x,vel)
   USE dsp
   USE store, ONLY : dim
   IMPLICIT NONE
   REAL(dp),INTENT(in) :: t
   REAL(dp),INTENT(in) :: x(dim)
   REAL(dp),INTENT(out) :: vel(dim)
   END SUBROUTINE velocity
 END INTERFACE
 
 REAL(dp) errmax,h,htemp,tnew
 REAL(dp) xerr(dim),xtemp(dim)

 h=htry ! set stepsize to initial value

 do
  call runge_kutta_step(x,vel,t,h,xtemp,xerr,velocity)
  errmax=maxval(abs(xerr(:)/xscal(:)))/eps ! scale relative to required tol
  if(errmax<=1.d0)exit ! Step succeeded. Go compute size of next step.
! Truncation error too large - reduce stepsize
  htemp=0.9d0*h*(errmax**(-0.25d0))
  h=sign(max(abs(htemp),0.1d0*abs(h)),h) ! No more than a factor of 10
  tnew=t+h
  if(tnew==t)call errstop_master('STEPPER_RK','Stepsize underflow.')
 enddo

! Compute size of next step. No more than a factor of 5 increase.
 if(errmax>1.89d-4)then ! where 1.89d-4 = (5/0.9)**(1/0.2)
  hnext=0.9d0*h*(errmax**(-0.2d0))
 else
  hnext=5.d0*h
 endif

 hdid=h
 t=t+h
 x(:)=xtemp(:)

 END SUBROUTINE stepper_rk


 SUBROUTINE runge_kutta_step(x,vel,t,h,xout,xerr,velocity)
!------------------------------------------------------------------------------!
! RUNGE_KUTTA_STEP                                                             !
! ================                                                             !
!                                                                              !
! Do one Runge-Kutta-Fehlberg step from time t to t+h, starting at position x. !
!                                                                              !
! Output : xout(1:dim) - final position                                        !
!          xerr        - estimated error in final position                     !
!------------------------------------------------------------------------------!
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: x(dim),vel(dim)
 REAL(dp),INTENT(in) :: t,h
 REAL(dp),INTENT(out) :: xout(dim),xerr(dim)

 INTERFACE
  SUBROUTINE velocity(t,x,vel)
   USE dsp
   USE store, ONLY : dim
   IMPLICIT NONE
   REAL(dp),INTENT(in) :: t
   REAL(dp),INTENT(in) :: x(dim)
   REAL(dp),INTENT(out) :: vel(dim)
  END SUBROUTINE velocity
 END INTERFACE

 REAL(dp),DIMENSION(dim) :: ak2,ak3,ak4,ak5,ak6,xtemp
 REAL(dp),PARAMETER :: a2=0.2d0,a3=0.3d0,a4=0.6d0,a5=1.d0,a6=0.875d0,&
  &b21=0.2d0,b31=3.d0/40.d0,b32=9.d0/40.d0,b41=0.3d0,b42=-0.9d0,b43=1.2d0,&
  &b51=-11.d0/54.d0,b52=2.5d0,b53=-70.d0/27.d0,b54=35.d0/27.d0,& 
  &b61=1631.d0/55296.d0,b62=175.d0/512.d0,b63=575.d0/13824.d0,&
  &b64=44275.d0/110592.d0,b65=253.d0/4096.d0,c1=37.d0/378.d0,&
  &c3=250.d0/621.d0,c4=125.d0/594.d0,c6=512.d0/1771.d0,&
  &dc1=c1-2825.d0/27648.d0,dc3=c3-18575.d0/48384.d0,&
  &dc4=c4-13525.d0/55296.d0,dc5=-277.d0/14336.d0,dc6=c6-0.25d0

 xtemp=x+b21*h*vel                                       ! first step

 call velocity(t+a2*h,xtemp,ak2)                         ! second step
 xtemp=x+h*(b31*vel+b32*ak2)                             

 call velocity(t+a3*h,xtemp,ak3)                         ! third step
 xtemp=x+h*(b41*vel+b42*ak2+b43*ak3)                     

 call velocity(t+a4*h,xtemp,ak4)                         ! fourth step
 xtemp=x+h*(b51*vel+b52*ak2+b53*ak3+b54*ak4)

 call velocity(t+a5*h,xtemp,ak5)                         ! fifth step
 xtemp=x+h*(b61*vel+b62*ak2+b63*ak3+b64*ak4+b65*ak5)

 call velocity(t+a6*h,xtemp,ak6)                         ! sixth step
 xout=x+h*(c1*vel+c3*ak3+c4*ak4+c6*ak6) ! accumulate increments with  
                                        ! proper weights

 xerr=h*(dc1*vel+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6)

 END SUBROUTINE runge_kutta_step


END MODULE runge_kutta

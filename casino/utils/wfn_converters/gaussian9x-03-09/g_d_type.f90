FUNCTION g_d_type(rvec)
!------------------------------------------------------------------!
! Evaluate a normalized (?) 3D d-type primitive Gaussian at the    !
! position, rvec.                                                  !
!------------------------------------------------------------------!
 USE paramfile
 IMPLICIT none
 REAL(KIND=dp),INTENT(in) :: rvec(3)
 INTEGER,SAVE :: icall=0 ! First call of function or not
 REAL(KIND=dp) g_d_type
 REAL(KIND=dp) rvec2(3) ! (x^2, y^2, z^2)
 REAL(KIND=dp) rmag2 ! x^2 + y^2 + z^2
 REAL(KIND=dp),PARAMETER :: alpha=0.057348d0, alpha7=alpha**7
! Constant parts of normalization factors for d functions
 REAL(KIND=dp),SAVE :: dnorm2 !,dnorm1,dnorm3

! Initialisation
 if(icall==0)then
  icall=1
!  dnorm1=(32.d0*oneover_pi_cubed*alpha7)**0.25d0
  dnorm2=(128.d0*oneover_pi_cubed*alpha7)**0.25d0
!  dnorm3=((2048.d0*oneover_pi_cubed*alpha7)**0.25d0)/sqrt(3.d0)
 endif

 rvec2=rvec*rvec
 rmag2=sum(rvec2)

! Calculate the contributions to a 'pure' d function...
! g_d_type=dnorm1*(rvec2(3)-rmag2) & ! z^2 - r^2
 g_d_type=dnorm2*(rvec2(1)-rvec2(2)) ! x^2 - y^2
! ! The last three parts have a different normalization...
! +dnorm3*(rvec(1)*rvec(2) &! xy
! +rvec(1)*rvec(3)         &! xz
! +rvec(2)*rvec(3))         ! yz

 g_d_type=g_d_type*exp(-alpha*rmag2)

END FUNCTION g_d_type

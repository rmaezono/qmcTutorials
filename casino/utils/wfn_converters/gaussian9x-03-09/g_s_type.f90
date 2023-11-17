FUNCTION g_s_type(rvec)
!----------------------------------------------------------------!
! Evaluate a normalized (?) 3D s-type primitive Gaussian at the  !
! position, rvec.                                                !
!----------------------------------------------------------------!
 USE paramfile
 IMPLICIT NONE
 REAL(KIND=dp),INTENT(in) :: rvec(3)
 REAL(KIND=dp) g_s_type
 REAL(KIND=dp) rmag2 ! x^2 + y^2 + z^2
 REAL(KIND=dp) anorm ! Normalization factor
 REAL(KIND=dp),PARAMETER :: alpha=1.36529000d+00,alpha3=alpha**3

 rmag2=sum(rvec*rvec)

 anorm=(8.d0*alpha3*oneover_pi_cubed)**0.25d0
 g_s_type=anorm*exp(-alpha*rmag2)

END FUNCTION g_s_type

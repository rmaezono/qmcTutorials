MODULE eval_density

 USE dsp
 IMPLICIT NONE
 PRIVATE
 PUBLIC eval_density_t0_1d,eval_density_t0_2d,eval_density_t0_3d

CONTAINS 

 SUBROUTINE eval_density_t0_1d(dentype,x,rho)
!----------------------------------------------------------------------------!
! Evaluate the initial density rho in 1D at t=0 and position x.              !
! The dentype parameter may be used to select different density types.       !
!----------------------------------------------------------------------------!
 USE run_control,ONLY : errstop
 USE store, ONLY : two_over_pi,four_over_pi,half_cell_x
 INTEGER,INTENT(in) :: dentype
 REAL(dp),INTENT(in) :: x(1)
 REAL(dp),INTENT(out) :: rho

 select case(dentype)
  case(0)
   rho=two_over_pi*sin(x(1))**2
  case(1)
   if(x(1)<half_cell_x)then
    rho=four_over_pi*sin(2.d0*x(1))**2 ! CHECK NORMALIZATION
   else
    rho=0.d0
   endif
  case default
   call errstop('EVAL_DENSITY_T0_1D','Invalid value of DENTYPE in input.')
 end select

 END SUBROUTINE eval_density_t0_1d


 SUBROUTINE eval_density_t0_2d(dentype,x,rho)
!----------------------------------------------------------------------------!
! Evaluate the initial density rho in 2D at t=0 and position x.              !
! The dentype parameter may be used to select different density types.       !
!----------------------------------------------------------------------------!
 USE run_control,ONLY : errstop
 USE choose_wfn,ONLY : eval_psi_squared_t0_2d
 USE store, ONLY : four_over_pi_squared,sixteen_over_pi_squared,&
  &pi_over_two,half_cell_x,half_cell_y
 INTEGER,INTENT(in) :: dentype
 REAL(dp),INTENT(inout) :: x(2)
 REAL(dp),INTENT(out) :: rho
 REAL(dp) psi2

 select case(dentype)
  case(0)
   rho=four_over_pi_squared*sin(x(1))**2*sin(x(2))**2
  case(1)
   if(x(1)<=half_cell_x.and.x(2)<=half_cell_y)then
    rho=sixteen_over_pi_squared*(sin(2.d0*x(1))**2*sin(2.0d0*x(2))**2)
   else
    rho=0.d0
   endif
  case(2)
   if(x(1)>=half_cell_x.and.x(2)<=half_cell_y)then
    rho=sixteen_over_pi_squared*(sin(2.d0*(x(1)-pi_over_two))**2*&
     &sin(2.d0*x(2))**2)
   else
    rho=0.d0
   endif
  case(3)
   if(x(1)<=half_cell_x.and.x(2)>half_cell_y)then
    rho=sixteen_over_pi_squared*(sin(2.d0*x(1))**2*&
     &sin(2.d0*(x(2)-pi_over_two))**2)
   else
    rho=0.d0
   endif
  case(4)
   if(x(1)>=half_cell_x.and.x(2)>half_cell_y)then
    rho=sixteen_over_pi_squared*(sin(2.d0*(x(1)-pi_over_two))**2* &
     &sin(2.d0*(x(2)-pi_over_two))**2)
   else
    rho=0.d0
   endif
  case(5)
   call eval_psi_squared_t0_2d(x,psi2)
   rho=four_over_pi_squared*sin(x(1))**2*sin(x(2))**2
   rho=rho*psi2
  case default
   call errstop('EVAL_DENSITY_T0_2D','Invalid value of DENTYPE in input.')
 end select

 END SUBROUTINE eval_density_t0_2d


 SUBROUTINE eval_density_t0_3d(dentype,x,rho)
!----------------------------------------------------------------------------!
! Evaluate the initial density rho in 3D at t=0 and position x.              !
! The dentype parameter may be used to select different density types.       !
!----------------------------------------------------------------------------!
! USE run_control,ONLY : errstop
 USE store, ONLY : eight_over_pi_cubed
 INTEGER,INTENT(in) :: dentype
 REAL(dp),INTENT(in) :: x(3)
 REAL(dp),INTENT(out) :: rho

! select case(dentype)
!  case(1)

 rho=eight_over_pi_cubed*sin(x(1))**2*sin(x(2))**2*sin(x(3))**2

!  case default
!   call errstop('EVAL_DENSITY_T0_3D','Invalid value of DENTYPE in input.')
! end select

 END SUBROUTINE eval_density_t0_3d

END MODULE eval_density

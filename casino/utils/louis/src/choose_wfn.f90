MODULE choose_wfn
 USE dsp
 USE sin_wfn
 USE scaled_sin_wfn
 IMPLICIT NONE
 PRIVATE
 PUBLIC eval_psi_squared_t0_1d,eval_psi_squared_t0_2d,eval_psi_squared_t0_3d,&
  &eval_psi_squared_1d,eval_psi_squared_2d,eval_psi_squared_3d,&
  &energy_and_variance

CONTAINS

 SUBROUTINE eval_psi_squared_t0_1d(x,rho)
 USE store,ONLY : wfn_type
 REAL(dp),INTENT(inout) :: x(1)
 REAL(dp),INTENT(inout) :: rho

 select case(trim(adjustl(wfn_type)))
  case('sine_wave')
   call sw_eval_psisq_t0_1d(x,rho)
  case('scaled_sine')
   call ss_eval_psisq_t0_1d(x,rho)
 end select

 END SUBROUTINE eval_psi_squared_t0_1d

 SUBROUTINE eval_psi_squared_t0_2d(x,rho)
 USE store,ONLY : wfn_type
 REAL(dp),INTENT(inout) :: x(2)
 REAL(dp),INTENT(inout) :: rho

 select case(trim(adjustl(wfn_type)))
  case('sine_wave')
   call sw_eval_psisq_t0_2d(x,rho)
  case('scaled_sine')
   call ss_eval_psisq_t0_2d(x,rho)
 end select

 END SUBROUTINE eval_psi_squared_t0_2d

 SUBROUTINE eval_psi_squared_t0_3d(x,rho)
 USE store,ONLY : wfn_type
 REAL(dp),INTENT(inout) :: x(3)
 REAL(dp),INTENT(inout) :: rho

 select case(trim(adjustl(wfn_type)))
  case('sine_wave')
   call sw_eval_psisq_t0_3d(x,rho)
  case('scaled_sine')
   call ss_eval_psisq_t0_3d(x,rho)
 end select

 END SUBROUTINE eval_psi_squared_t0_3d

 SUBROUTINE eval_psi_squared_1d(t,x,rho)
 USE store,ONLY : wfn_type
 REAL(dp),INTENT(inout) :: t
 REAL(dp),INTENT(inout) :: x(1)
 REAL(dp),INTENT(inout) :: rho

 select case(trim(adjustl(wfn_type)))
  case('sine_wave')
   call sw_eval_psisq_1d(t,x,rho)
  case('scaled_sine')
   call ss_eval_psisq_1d(t,x,rho)
 end select

 END SUBROUTINE eval_psi_squared_1d

 SUBROUTINE eval_psi_squared_2d(t,x,rho)
 USE store,ONLY : wfn_type
 REAL(dp),INTENT(inout) :: t
 REAL(dp),INTENT(inout) :: x(2)
 REAL(dp),INTENT(inout) :: rho

 select case(trim(adjustl(wfn_type)))
  case('sine_wave')
   call sw_eval_psisq_2d(t,x,rho)
  case('scaled_sine')
   call ss_eval_psisq_2d(t,x,rho)
 end select

 END SUBROUTINE eval_psi_squared_2d

 SUBROUTINE eval_psi_squared_3d(t,x,rho)
 USE store,ONLY : wfn_type
 REAL(dp),INTENT(inout) :: t
 REAL(dp),INTENT(inout) :: x(3)
 REAL(dp),INTENT(inout) :: rho

 select case(trim(adjustl(wfn_type)))
  case('sine_wave')
   call sw_eval_psisq_3d(t,x,rho)
  case('scaled_sine')
   call ss_eval_psisq_3d(t,x,rho)
 end select

 END SUBROUTINE eval_psi_squared_3d

 SUBROUTINE energy_and_variance
 USE store, ONLY : wfn_type
 select case(trim(adjustl(wfn_type)))
  case('sine_wave')
   call sw_energy_and_variance
  case('scaled_sine')
   call ss_energy_and_variance
 end select

 END SUBROUTINE energy_and_variance

END MODULE choose_wfn


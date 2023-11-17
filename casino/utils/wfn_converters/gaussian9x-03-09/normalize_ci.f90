SUBROUTINE normalize_ci
!----------------------------------------------------------------------------!
! Re-normalize the CIS expansion of each excited state to allow for the      !
! fact that we may have a rather incomplete expansion (how incomplete will   !
! depend on the value of N in Iop(9/40=N) in the Gaussian94 route section).  !
!----------------------------------------------------------------------------!
 USE cis_data
 USE g94_wavefunction, ONLY: SPIN, ispin_lim
 IMPLICIT NONE
 INTEGER ispin
 REAL(KIND=dp) wnorm,truenorm,sqrt_half

 if(icis_out/=0)then ! Only check normalization of CIS expansion if we
                     ! are dealing with an excited state rather than
                     ! the ground state

  sqrt_half=sqrt(0.5d0)

! Normalize the CIS coefficients for the selected excited state, icis_out

  wnorm=0.0d0
  do ispin=1,ispin_lim
   wnorm=wnorm+sum(ci_coeff(1:Nconfig(ispin),ispin)* &
    &ci_coeff(1:Nconfig(ispin),ispin))
  enddo
  write(*,fmt="('Old normalization of excited state ',i2,' = ',e12.5)") &
   &icis_out,wnorm

  ci_norm=wnorm ! Store normalization of original expansion for later
                ! analysis of the wave function
  truenorm=1.d0/sqrt(wnorm)
! If this is a spin-restricted calculation for an excited state
! then we are only explicitly dealing with half of the wave function
! and must take this into account when we normalize it.
  if(.not.spin)then
! ('singlet' not set unless it is a spin-restricted calc.)
!if(singlet)truenorm=sqrt_half*truenorm
   ci_norm=2.0d0*ci_norm
   truenorm=sqrt_half*truenorm
  endif

  do ispin=1,ispin_lim
   ci_coeff(1:Nconfig(ispin),ispin)=truenorm*ci_coeff(1:Nconfig(ispin),ispin)
  enddo

  wnorm=0
  do ispin=1,ispin_lim,1
   wnorm=wnorm + SUM(ci_coeff(1:Nconfig(ispin),ispin)* &
    &ci_coeff(1:Nconfig(ispin),ispin))
  enddo
  write(*,fmt="('New normalization of excited state ',i2,' = ',e12.5)") &
   &icis_out,wnorm
 endif

END SUBROUTINE normalize_ci

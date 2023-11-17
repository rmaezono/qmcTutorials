MODULE integ_params
 USE paramfile
 IMPLICIT none

! MO to output for plotting, it is also the one whose normalization is tested.
 INTEGER,SAVE :: mo_ref
! Spin of MO to plot
 INTEGER,SAVE :: mo_spin
! If (CIS) then this is the excited state to output
 INTEGER,SAVE :: Nex_plt
! Number of data points to scan in each dimension
 INTEGER,PARAMETER :: irange=2000
! Step size between points (same for each dimension)
 REAL(KIND=dp),PARAMETER :: dx=0.005d0
! Associated volume element for integration (in normalization_check)
 REAL(KIND=dp),PARAMETER :: dv=dx*dx*dx

END MODULE integ_params

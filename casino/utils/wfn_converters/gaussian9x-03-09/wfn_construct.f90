SUBROUTINE wfn_construct
!-------------------------------------------------------------------!
! Calculate the wave function using knowledge of the coefficients   !
! defining the  molecular orbital and the angular momentum of each  !
! basis function.                                                   !
!-------------------------------------------------------------------!
 USE cis_data
 USE g94_wavefunction
 USE integ_params
 IMPLICIT NONE

 INTEGER ix
 REAL(KIND=dp) x,rtemp
 REAL(KIND=dp) rvec(3) ! Position vector of point at which to calculate wfn
 REAL(KIND=dp) origin(3) ! Origin of coordinates
 REAL(KIND=dp) Radial ! The value of the radial wfn
 REAL(KIND=dp) psi ! Function to calculate wfn at rvec
 LOGICAL plotradial
 CHARACTER(1) inchar

 write(6,fmt="(/'Plot (R)adial or (A)bsolute MO?')")
 read(*,fmt="(a1)")inchar
 call capitalise(inchar,1)

 if(inchar=='R')then
  plotradial=.true.
 elseif(inchar=='A')then
  plotradial=.false.
 else
  write(6,*) 'Defaulting to absolute plot...'
  plotradial=.false.
 endif

 open(unit=24,file='wfn.out',status='unknown')
 write(24,fmt="('# Gaussian wave function; x,Psi:')")

! Set coordinates of starting point for scan
! origin=(/0.05d0,0.1d0,-3.18d0/)
 origin=(/0.d0,0.d0,0.d0/)
! The parameters determining the extent and granularity of the plot
! are set in the integ_params module.
 do ix=0,irange

! Scan along x axis...
! x=real(ix,dp)*dx
! rvec=(/x,0.d0,0.d0/)

! Scan along (1,0,1)/sqrt(2) axis
  x=real(ix)*dx*0.707106781d0
  rvec=(/x,0.d0,x/)
  rvec=rvec+origin

! Calculate the value of the MO at this point
  rtemp = psi(rvec,mo_ref,mo_spin)

! Round to zero if answer is very small to prevent plotting programs
! getting confused
  if(abs(rtemp)<1.d-20)rtemp=0.d0

! This part for simple x axis scan...
! if(plotradial)then
!  Radial=two_rtpi*x*rtemp ! Radial wave function
! else
!  Radial=rtemp             ! Full wave function
! endif
! write(24,fmt="(e13.4,2x,e13.4)")x,Radial

! This part needed if scanning along (1,1,0)/sqrt(2)...
  if(plotradial)then
   Radial=two_rtpi*x*1.414213562373d0*rtemp ! Radial wave function
  else
   Radial=rtemp ! Full wave function
  endif
  write(24,fmt="(e13.4,2x,e13.4)")1.414213562373d0*x,Radial

 enddo

 close(unit=24)

END SUBROUTINE wfn_construct

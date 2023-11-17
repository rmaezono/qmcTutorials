SUBROUTINE normalization_check
!----------------------------------------!
! Check the normalization of a function. !
!----------------------------------------!
 USE cis_data
 USE g94_wavefunction
 USE integ_params
 IMPLICIT NONE

 INTEGER i1,j1,k1
 REAL(KIND=dp) x1,y1,z1,ainteg,tmp
 REAL(KIND=dp) posn(3) ! position vector of current point
 REAL(KIND=dp) psi ! The function to be integrated

 open(unit=21,file='norm.log',status='unknown')

 ainteg=0.d0

! The parameters determining the range and step size for the integration are
! set in the integ_params module.

 do i1=-irange,irange
  if(i1/=0)then
! This bit ensures that the point at which the wave function is calculated
! is at the centre of the volume element.
   x1=sign((real(abs(i1),dp)-0.5d0),real(i1,dp))*dx
   do j1=-irange,irange
    if(j1/=0)then
     y1=sign((real(abs(j1),dp)-0.5d0),real(j1,dp))*dx
     do k1=-irange,irange
      if(k1/=0)then
       z1=sign((real(abs(k1),dp)-0.5d0),real(k1,dp))*dx

       posn=(/x1,y1,z1/)

       tmp=psi(posn,mo_ref,mo_spin)

       if(abs(tmp)>1.d0)then
        write(6,*)'psi = ',tmp
        write(6,*)'mo_ref = ',mo_ref,' mo_spin = ',mo_spin
        write(6,*)posn
       endif

       ainteg=ainteg+tmp*tmp*dv

      endif
     enddo
    endif
   enddo
  endif
 enddo

 if(CIS)then
  write(21,fmt="(/'Integration of state ',i2,' over interval ',&
   &f6.1,' to ',f6.1,' with dx = ',f8.4,' gives '&
   &,e14.6)")Nex_plt,-dx*irange,dx*irange,dx,ainteg
  else
   write(21,fmt="(/'Integration of MO ',i2,', spin ',i1,' over interval '&
    &,f6.1,' to ',f6.1,' with'/'dx = ',f8.4,' gives '&
    &,e14.6)")mo_ref,mo_spin,-dx*irange,dx*irange,dx,ainteg
  endif

 close(unit=21)

END SUBROUTINE normalization_check

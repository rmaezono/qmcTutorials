SUBROUTINE con_coeffs
!------------------------------------------------------------!
! Multiply the part of the normalization factor that is      !
! common to each shell into the contraction coefficients and !
! store them in an accessible form.                          !
!                                                            !
! Changes                                                    !
! -------                                                    !
! 3.2015 KDD - Fixed issue with harmonic g functions.        !
!------------------------------------------------------------!
 USE g94_wavefunction
 IMPLICIT none

 INTEGER ic,iprim,ishell
 REAL(KIND=dp) anorm,tv1,tv2,tv3,tv4,tv5,tv6,tv7

 tv1=(8.d0*oneover_pi_cubed)**0.25d0   ! (2^3 / PI^3)**.25
 tv2=2.d0*tv1                          ! (2^7 / PI^3)**.25
 tv3=2.d0*tv2                          ! (2^11 / PI^3)**.25
 tv4=(512.d0*oneover_pi_cubed)**0.25d0
 tv5=2.0d0*tv3                         ! (2^15 / PI^3)**.25
! tv6=tv2/sqrt(3.d0)                    ! (2^7 / 9*PI^3)**.25
 tv6=(2048.d0*oneover_pi_cubed)**0.25d0
 tv7=2.0d0*tv3/sqrt(3.d0)

 ic=0
 do ishell=1,Nshells
  do iprim=1,Nprim(ishell)
   ic=ic+1

! Identify shell type and hence the part of the normalization factor common to
! all associated basis functions
   shelltype: select case(Lshell(ishell))

    case(0) ! S shell
     anorm=tv1*(shexpnt(ic)**0.75d0)

    case(-1) ! SP shell but this is the con. coeff. for S part
     anorm=tv1*(shexpnt(ic)**0.75d0)

    case(1) ! P shell
     anorm=tv2*(shexpnt(ic)**1.25d0)

    case(-2) ! Harmonic D shell
     anorm=tv2*(shexpnt(ic)**1.75d0)

    case(2) ! Cartesian D shell
     anorm=tv3*(shexpnt(ic)**1.75d0)

    case(-3) ! Harmonic F shell.
     anorm=tv4*(shexpnt(ic)**2.25d0)

    case(3) ! Cartesian F shell
     anorm=tv5*(shexpnt(ic)**2.25d0)

    case(-4) ! Harmonic G shell
     anorm=tv6*(shexpnt(ic)**2.75d0)

    case(4) ! Cartesian G shell
     anorm=tv7*(shexpnt(ic)**2.75d0)

    case default
     write(6,fmt="('Shells higher than g not supported.')")

   end select shelltype

! Multiply this common normalization factor into the contraction
! coefficients and store the latter more conveniently
   c_prim(ic)=anorm*c_prim(ic)
   con_coeff(iprim,ishell)=c_prim(ic)

  enddo
 enddo

 if(num_sp>0)then
! We also have SP shells so do the same for them...
  ic=0
  do ishell=1,Nshells
   do iprim=1,Nprim(ishell)
    ic=ic+1
    anorm=tv2*(shexpnt(ic)**1.25d0)
    c2_prim(ic)=anorm*c2_prim(ic)
    sp_coeff(iprim,ishell)=c2_prim(ic)
   enddo
  enddo

 endif

END SUBROUTINE con_coeffs

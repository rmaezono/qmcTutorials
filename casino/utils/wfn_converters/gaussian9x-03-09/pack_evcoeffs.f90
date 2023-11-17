SUBROUTINE pack_evcoeffs(temp,Ndim)
!----------------------------------------------------------------------------!
! Store the eigenvector coefficients in a directly accessible form with      !
! alpha and beta MO coeff.'s stored separately.  At this stage also          !
! multiply in the remaining normalization factors - required because the     !
! normalization of d and higher shells is not the same for all basis         !
! functions. The common parts of their normalization are included in the     !
! contraction coefficients.                                                  !
!                                                                            !
! The multiplication of the normalization factors and ev coeffs could be     !
! linearized at some stage - store the normalization factors for each        !
! basis function type in arrays and send to BLAS.                            !
!                                                                            !
! Changes                                                                    !
! -------                                                                    !
! 3.2015 KDD - Fixed issue with harmonic g functions.                        !
!----------------------------------------------------------------------------!
 USE g94_wavefunction
 IMPLICIT none

 INTEGER,INTENT(in) :: Ndim
 REAL(KIND=dp),INTENT(in) :: temp(Ndim)
 INTEGER ispin,imo,ifun,ic,ishell
 REAL(KIND=dp) oneover_rt3,oneover_rt5,oneover_rt35,oneover_rt315,sqrt2,tv1,&
  &tv2,tv3,tv4,tv5,tv6,tv7,harm_f_coeffs(7),harm_g_coeffs(9)

! Precalculate some constant factors used below
 oneover_rt3=1.d0/sqrt(3.d0)
 oneover_rt5=1.d0/sqrt(5.d0)
 oneover_rt35=1.d0/sqrt(35.d0)
 oneover_rt315=1.d0/sqrt(315.d0)
 sqrt2=sqrt(2.d0)
 tv1=2.d0**1.5d0
 tv2=2.d0*oneover_rt3
 tv3=sqrt(8.d0/15.d0)
 tv4=2.d0/sqrt(7.d0)
 tv5=1.d0/sqrt(15.d0)
 tv6=4.d0*sqrt(2.d0)*oneover_rt5
 tv7=sqrt(8.d0)*oneover_rt315

 harm_f_coeffs(1)=sqrt(8.d0/15.d0)        ! m =  0
 harm_f_coeffs(2)=2.d0*oneover_rt5/3.d0   ! m =  1
 harm_f_coeffs(3)=harm_f_coeffs(2)        ! m = -1
 harm_f_coeffs(4)=sqrt2/15.d0             ! m =  2
 harm_f_coeffs(5)=(2.d0**1.5d0)/30.d0     ! m = -2
 harm_f_coeffs(6)=oneover_rt3/15.d0       ! m =  3
 harm_f_coeffs(7)=harm_f_coeffs(6)        ! m = -3

! harm_g_coeffs(1)=8.d0*oneover_rt35       ! m =  0
! harm_g_coeffs(2)=0.4d0*sqrt2*tv4         ! m =  1
! harm_g_coeffs(3)=harm_g_coeffs(2)        ! m = -1
! harm_g_coeffs(4)=4.d0/(15.d0*sqrt(7.d0)) ! m =  2
! harm_g_coeffs(5)=harm_g_coeffs(4)        ! m = -2
! harm_g_coeffs(6)=(2.d0**1.5d0)/105.d0    ! m =  3
! harm_g_coeffs(7)=harm_g_coeffs(6)        ! m = -3
! harm_g_coeffs(8)=1.d0/105.d0             ! m =  4
! harm_g_coeffs(9)=1.d0/420.d0             ! m = -4
! KDD fix 3.2015
 harm_g_coeffs(1)=sqrt(16.d0/105.d0)              ! m =  0
 harm_g_coeffs(2)=harm_g_coeffs(1)/sqrt(10.d0)    ! m =  1
 harm_g_coeffs(3)=harm_g_coeffs(2)                ! m = -1
 harm_g_coeffs(4)=harm_g_coeffs(1)/sqrt(180.d0)   ! m =  2
 harm_g_coeffs(5)=harm_g_coeffs(4)                ! m = -2
 harm_g_coeffs(6)=harm_g_coeffs(1)/sqrt(2520.d0)  ! m =  3
 harm_g_coeffs(7)=harm_g_coeffs(6)                ! m = -3
 harm_g_coeffs(8)=harm_g_coeffs(1)/sqrt(20160.d0) ! m =  4
 harm_g_coeffs(9)=harm_g_coeffs(8)                ! m = -4

 ic=0
 do ispin=1,ispin_lim
  do imo=1,Nmo ! Loop over all MOs of spin ispin_lim
   ifun=0 ! Counter for current basis function

   do ishell=1,Nshells
    ifun=ifun+1

    shelltype2: select case(Lshell(ishell))
     case(0) ! S shell
      ic=ic+1
      evcoeff1(ifun,imo,ispin)=temp(ic)
     case(-1) ! SP shell
      evcoeff1(ifun:ifun+3,imo,ispin)=temp(ic+1:ic+4)
      ifun=ifun+3
      ic=ic+4
     case(1) ! P shell
      evcoeff1(ifun:ifun+2,imo,ispin)=temp(ic+1:ic+3)
      ifun=ifun+2
      ic=ic+3
     case(-2) ! Harmonic D shell
      evcoeff1(ifun,imo,ispin)  =oneover_rt3*temp(ic+1)
      evcoeff1(ifun+1,imo,ispin)=2.0d0*temp(ic+2)
      evcoeff1(ifun+2,imo,ispin)=2.0d0*temp(ic+3)
      evcoeff1(ifun+3,imo,ispin)=temp(ic+4)
      evcoeff1(ifun+4,imo,ispin)=2.0d0*temp(ic+5)
      ifun=ifun+4
      ic=ic+5
     case(2) ! Cartesian D shell
      evcoeff1(ifun:ifun+2,imo,ispin)  =oneover_rt3*temp(ic+1:ic+3)
      evcoeff1(ifun+3:ifun+5,imo,ispin)=temp(ic+4:ic+6)
      ifun=ifun+5
      ic=ic+6
     case(-3) ! Harmonic F shell
      evcoeff1(ifun:ifun+6,imo,ispin)=harm_f_coeffs*temp(ic+1:ic+7)
      ifun=ifun+6
      ic=ic+7
     case(3) ! Cartesian F shell
      ic=ic+1
      evcoeff1(ifun:ifun+2,imo,ispin)  =tv5*temp(ic:ic+2)
      ic=ic+3
      evcoeff1(ifun+3:ifun+8,imo,ispin)=oneover_rt3*temp(ic:ic+5)
      ic=ic+6
      evcoeff1(ifun+9,imo,ispin)=temp(ic)
      ifun=ifun+9
     case(-4) ! Harmonic G shell
      evcoeff1(ifun:ifun+8,imo,ispin)=harm_g_coeffs*temp(ic+1:ic+9)
      ifun=ifun+8
      ic=ic+9
     case(4) ! Cartesian G shell
      ic=ic+1
      evcoeff1(ifun,imo,ispin)   = oneover_rt35*temp(ic); ic=ic+1
      evcoeff1(ifun+1,imo,ispin) = oneover_rt5*temp(ic) ; ic=ic+1
      evcoeff1(ifun+2,imo,ispin) = oneover_rt3*temp(ic) ; ic=ic+1
      evcoeff1(ifun+3,imo,ispin) = oneover_rt5*temp(ic) ; ic=ic+1
      evcoeff1(ifun+4,imo,ispin) = oneover_rt35*temp(ic); ic=ic+1
      evcoeff1(ifun+5,imo,ispin) = oneover_rt5*temp(ic) ; ic=ic+1
      evcoeff1(ifun+6,imo,ispin) =             temp(ic) ; ic=ic+1
      evcoeff1(ifun+7,imo,ispin) =             temp(ic) ; ic=ic+1
      evcoeff1(ifun+8,imo,ispin) = oneover_rt5*temp(ic) ; ic=ic+1
      evcoeff1(ifun+9,imo,ispin) = oneover_rt3*temp(ic) ; ic=ic+1
      evcoeff1(ifun+10,imo,ispin) =            temp(ic) ; ic=ic+1
      evcoeff1(ifun+11,imo,ispin) = oneover_rt3*temp(ic); ic=ic+1
      evcoeff1(ifun+12:ifun+13,imo,ispin) = oneover_rt5*temp(ic:ic+1)
      ic=ic+2
      evcoeff1(ifun+14,imo,ispin) = oneover_rt35*temp(ic)
      ifun=ifun+14
     case default
      write(6,fmt="('h and higher functions not implemented in &
       &pack_evcoeffs.')")

    end select shelltype2

   enddo
  enddo
 enddo

END SUBROUTINE pack_evcoeffs

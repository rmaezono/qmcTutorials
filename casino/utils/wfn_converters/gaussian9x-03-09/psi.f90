FUNCTION psi(rvec,imo,ispin)
!------------------------------------------------------------------------------!
! Calculate the value of MO `imo' (spin ispin) at rvec=(x,y,z). All            !
! normalization factors are already included in the eigenvector and            !
! contraction coefficients (done in pack_evcoeffs and con_coeffs respectively).!
!------------------------------------------------------------------------------!
 USE g94_wavefunction
 IMPLICIT none
 INTEGER,INTENT(in) :: imo   ! Which MO to calculate
 INTEGER,INTENT(in) :: ispin ! =1 if alpha, 2 if beta MO/band
 REAL(KIND=dp),INTENT(in) :: rvec(3) ! Point at which to calculate

 REAL(KIND=dp) psi

! Just counters: ibasfun keeps track of which basis fn we are on and ishell
! which (contracted) shell
 INTEGER i,ibasfun,ishell
! Loop over contractions
 INTEGER icntr
! Counter for which primitive we're on
 INTEGER iprim
! Exponent - only need one since all functions in a given shell have the
! same exponent.
 REAL(KIND=dp) alpha
! Initially the full exponent of the Gaussian and then, if not too big, the
! value of the Gaussian itself
! REAL(KIND=dp) gaussn
! Position vector of the point at which the wave function is being evaluated
! wrt to each shell centre.
 REAL(KIND=dp) relvec(3,Nshells)
! Ditto but components squared
 REAL(KIND=dp) relvec2(3,Nshells)
! Magnitude of relative position vector squared
 REAL(KIND=dp) relmag2(Nshells)

! Temporary variable for holding product of normalization factor and a
! primitive Gaussian
 REAL(KIND=dp) gval
! Used for summing up the basis function contributions to a primitive Gaussian.
! e.g. for a p-type there are 3 basis functions.  For a Cartesian i-type
! function there are 28.
 REAL(KIND=dp) cntrn(28)
! Holds the individual basis function values before multiplication by the
! Gaussian and contraction coefficient followed by addition to sum over
! basis functions
 REAL(KIND=dp) primfn(28)
! Flag that identifies whether a Gaussian is big enough to be worth evaluating
! LOGICAL evaluate

 psi=0.d0

 do i=1,Nshells
  relvec(:,i)=rvec(:)-shll_posn(:,i)
  relvec2(:,i)=relvec(:,i)*relvec(:,i)
  relmag2(i)=sum(relvec2(:,i))
 enddo

 iprim=0
 ibasfun=0

 do ishell=1,Nshells
  ibasfun=ibasfun+1

! Apply some crude screening...
! Would be more efficient if exponents were sorted such that they were
! always ordered in decreasing size.
! gaussn=alpha*relmag2(ishell)
! if(gaussn<15.d0)then
!  evaluate=.true.
!  gaussn=exp(-gaussn)
! else
!  evaluate=.false.
! endif

  function_type: select case(Lshell(ishell))

   case(0) ! s type
    cntrn(1)=0.d0
    do icntr=1,Nprim(ishell)
     iprim=iprim+1      ! Counter for current primitive
     alpha=shexpnt(iprim) ! Exponent of this primitive
     cntrn(1)=cntrn(1)+con_coeff(icntr,ishell)*exp(-alpha*relmag2(ishell))
    enddo

! Add the basis function constructed from this contraction, multiplied
! by the appropriate MO expansion coefficient
    psi=psi+evcoeff1(ibasfun,imo,ispin)*cntrn(1)

   case(1) ! p type

    cntrn(1:3)=0.d0
    do icntr=1,Nprim(ishell)
     iprim=iprim+1
     alpha=shexpnt(iprim)
     gval=con_coeff(icntr,ishell)*exp(-alpha*relmag2(ishell))
     cntrn(1:3)=cntrn(1:3)+relvec(:,ishell)*gval
    enddo
    psi=psi+sum(evcoeff1(ibasfun:ibasfun+2,imo,ispin)*cntrn(1:3))
    ibasfun=ibasfun+2

   case(-1) ! sp type

    cntrn(1:4)=0.d0
    do icntr=1,Nprim(ishell)
     iprim=iprim+1
     alpha=shexpnt(iprim)
     gval=exp(-alpha*relmag2(ishell))
     cntrn(1)=cntrn(1)+con_coeff(icntr,ishell)*gval
     gval=gval*sp_coeff(icntr,ishell)
     cntrn(2:4)=cntrn(2:4)+relvec(:,ishell)*gval
    enddo

    psi=psi+sum(evcoeff1(ibasfun:ibasfun+3,imo,ispin)*cntrn(1:4))
    ibasfun=ibasfun+3

   case(-2) ! 'pure' d type - 5 basis functions

    cntrn(1:5)=0.d0
    do icntr=1,Nprim(ishell)
     iprim=iprim+1
     gval=shexpnt(iprim)
! Calculate the contributions to a 'pure' d function...
! Gaussian94 outputs the coefficients for these in the following order:
! 3z^2-r^2, xz, yz, x^2-y^2, xy, which corresponds to m=0,m=1,m=-1,...
! The harmonic coefficients for these functions, unlike those for the
! functions of higher AM, are not present here as they are included in
! the eigenvector coefficients.  This is because it is what the QMC
! code expects.
     primfn(1)=3.d0*relvec2(3,ishell)-relmag2(ishell) ! 3z^2 - r^2
     primfn(2)=relvec(1,ishell)*relvec(3,ishell) ! xz
     primfn(3)=relvec(2,ishell)*relvec(3,ishell) ! yz
     primfn(4)=relvec2(1,ishell)-relvec2(2,ishell) ! x^2 - y^2
     primfn(5)=relvec(1,ishell)*relvec(2,ishell) ! xy
     gval=con_coeff(icntr,ishell)*exp(-gval*relmag2(ishell))
     cntrn(1:5)=cntrn(1:5)+gval*primfn(1:5)
    enddo

    psi=psi+sum(evcoeff1(ibasfun:ibasfun+4,imo,ispin)*cntrn(1:5))
    ibasfun=ibasfun+4

   case(2) ! Cartesian d type - 6 basis functions

    cntrn=0.d0
    do icntr=1,Nprim(ishell)
     iprim=iprim+1
     gval=shexpnt(iprim)
     gval=con_coeff(icntr,ishell)*exp(-gval*relmag2(ishell))
! Calculate the six contributions to a cartesian d function...
! The coefficients for these are output by G94 in the following order:
! x^2, y^2, z^2, xy, xz, yz
     primfn(1:3)=relvec2(1:3,ishell) ! x^2, y^2 and z^2
     primfn(4)=relvec(1,ishell)*relvec(2,ishell) ! xy
     primfn(5)=relvec(1,ishell)*relvec(3,ishell) ! xz
     primfn(6)=relvec(2,ishell)*relvec(3,ishell) ! yz
     cntrn(1:6)=cntrn(1:6)+gval*primfn(1:6)
    enddo

    psi=psi+sum(evcoeff1(ibasfun:ibasfun+5,imo,ispin)*cntrn(1:6))
    ibasfun=ibasfun+5

   case(-3) ! Harmonic f function - 7 basis functions

    ! m=0
    primfn(1)=(5.d0*relvec2(3,ishell)-3.d0*relmag2(ishell))*relvec(3,ishell)
    primfn(1)=0.5d0*primfn(1)
    ! m=1, m=-1
    primfn(2:3)=(5.d0*relvec2(3,ishell)-relmag2(ishell))*relvec(1:2,ishell)
    primfn(2:3)=1.5d0*primfn(2:3)
    ! m=2
    primfn(4)=relvec(3,ishell)*(relvec2(1,ishell)-relvec2(2,ishell))
    primfn(4)=15.d0*primfn(4)
    ! m=-2
    primfn(5)=relvec(1,ishell)*relvec(2,ishell)*relvec(3,ishell)
    primfn(5)=30.d0*primfn(5)
    ! m=3
    primfn(6)=relvec(1,ishell)*(15.d0*relvec2(1,ishell)-4.5d1*relvec2(2,ishell))
    ! m=-3
    primfn(7)=relvec(2,ishell)*(4.5d1*relvec2(1,ishell)-15.d0*relvec2(2,ishell))
    cntrn=0.d0
    do icntr=1,Nprim(ishell)
     iprim=iprim+1
     gval=shexpnt(iprim)
     gval=con_coeff(icntr,ishell)*exp(-gval*relmag2(ishell))
     cntrn(1:7)=cntrn(1:7)+gval*primfn(1:7)
    enddo

    psi=psi+sum(evcoeff1(ibasfun:ibasfun+6,imo,ispin)*cntrn(1:7))
    ibasfun=ibasfun+6

   case(-4) ! Harmonic g function - 9 basis functions

    ! m=0
    primfn(1)=3.d0*relmag2(ishell)*(relmag2(ishell)-10.d0*relvec2(3,ishell))+ &
     &35.d0*relvec2(3,ishell)*relvec2(3,ishell)
    primfn(1)=0.125d0*primfn(1)
    ! m=1 and m=-1
    primfn(2:3)=2.5d0*relvec(1:2,ishell)*relvec(3,ishell)* &
     &(7.d0*relvec2(3,ishell)-3.d0*relmag2(ishell))
    gval=7.d0*relvec2(3,ishell)-relmag2(ishell)
    ! m=2
    primfn(4)=7.5d0*gval*(relvec2(1,ishell)-relvec2(2,ishell))
    ! m=-2
    primfn(5)=15.d0*relvec(1,ishell)*relvec(2,ishell)*gval
    ! m=3
    primfn(6)=relvec(1,ishell)*relvec(3,ishell)* &
     &(105.d0*relvec2(1,ishell)-315.d0*relvec2(2,ishell))
    ! m=-3
    primfn(7)=relvec(2,ishell)*relvec(3,ishell)* &
     &(315.d0*relvec2(1,ishell)-105.d0*relvec2(2,ishell))
    ! m=4
    primfn(8)=relvec2(1,ishell)*(relvec2(1,ishell)- &
     &6.0d0*relvec2(2,ishell))+relvec2(2,ishell)*relvec2(2,ishell)
    primfn(8)=105.d0*primfn(8)
    ! m=-4
    primfn(9)=420.d0*relvec(1,ishell)*relvec(2,ishell)* &
     &(relvec2(1,ishell)-relvec2(2,ishell))

    cntrn=0.d0
    do icntr=1,Nprim(ishell)
     iprim=iprim+1
     gval=con_coeff(icntr,ishell)*exp(-shexpnt(iprim)*relmag2(ishell))
     cntrn(1:9)=cntrn(1:9)+gval*primfn(1:9)
    enddo

    psi=psi+sum(evcoeff1(ibasfun:ibasfun+8,imo,ispin)*cntrn(1:9))
    ibasfun=ibasfun+8

   case default
    write(*,fmt="('Oh dear, L(',i2,') = ',i3)") ibasfun,imo,Lshell(ishell)
    stop

  end select function_type

 enddo

END FUNCTION psi

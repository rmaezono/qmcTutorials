MODULE random_numbers
!------------------------------------------------------------------------!
! RANDOM_NUMBERS                                                         !
! --------------                                                         !
! Module to handle all tasks pertaining to the generation of             !
! pseudo-random number sequences.                                        !
!                                                                        !
! Uses the RANLUX algorithm.                                             !
! (see header notes in routine 'ranlux' below).                          !
!                                                                        !
! MDT 8.2002                                                             !
!------------------------------------------------------------------------!
 USE dsp
 USE parallel
 USE run_control, ONLY : errstop
 USE store,       ONLY : o,twopi
 IMPLICIT NONE
 PRIVATE
 PUBLIC ranx,ranx_twopi,ranx_pm,initialize_random,default_seed,randomseed

 INTEGER,SAVE :: current_ran
 INTEGER,PARAMETER :: num_ranx=63,ranluxlevel=3
 REAL r(num_ranx)

! Variables used in the RANLUX algorithm
 INTEGER,PARAMETER :: maxlev=4,lxdflt=3,itwo24=2**24,icons=2147483563
 INTEGER :: ndskip(0:maxlev)=(/0,24,73,199,365/),igiga=1000000000,i24=24,  &
  &j24=10,in24=0,kount=0,mkount=0,default_seed=31415927
 INTEGER,SAVE :: next(24),luxlev=lxdflt,nskip,jseed,randomseed
 REAL,SAVE :: seeds(24),carry=0.,twom24,twom12
 LOGICAL,SAVE :: notyet=.true.


CONTAINS


 SUBROUTINE initialize_random
!------------------------------------------------------------------------!
! Generate unique seed for each node. Use Park-Miller generator to       !
! generate random default seeds for multiple nodes from the master seed  !
! JSDFLT. If restarting from an earlier calc the generator should        !
! instead reconstruct its position in the sequence from the saved state  !
! in the qmc.ran file.                                                   !
!------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i
 REAL(dp) temp

! Initialize generator on each node with the unique seed.
 call ranlux_init(ranluxlevel,randomseed,0,0,.true.,.false.)
! Flush random number buffer.
 current_ran=0
 do i=1,24
  temp=ranx() ! otherwise first 24 numbers used are independent of ranluxlevel.
 enddo

 END SUBROUTINE initialize_random


 REAL(dp) FUNCTION ranx()
!--------------------------------------------------------------------------!
! Generate uniformly distributed random number in the range 0 --> 1 .      !
! Make generator return random numbers num_ranx at a time for efficiency.  !
!--------------------------------------------------------------------------!
 IMPLICIT NONE
 if(current_ran==0)then
  call ranlux(r,num_ranx)
  current_ran=num_ranx
 endif
 ranx=dble(r(current_ran))
 current_ran=current_ran-1
 END FUNCTION ranx


 REAL(dp) FUNCTION ranx_twopi()
!--------------------------------------------------------------------------!
! Generate uniformly distributed random number in the range 0 --> 2pi.     !
! Make generator return random numbers num_ranx at a time for efficiency.  !
!--------------------------------------------------------------------------!
 IMPLICIT NONE
 if(current_ran==0)then
  call ranlux(r,num_ranx)
  current_ran=num_ranx
 endif
 ranx_twopi=dble(r(current_ran))*twopi
 current_ran=current_ran-1
 END FUNCTION ranx_twopi


 REAL(dp) FUNCTION ranx_pm()
!--------------------------------------------------------------------------!
! Generate uniformly distributed random number in the range -1 --> +1.     !
! Make generator return random numbers num_ranx at a time for efficiency.  !
!--------------------------------------------------------------------------!
 IMPLICIT NONE
 
 if(ranx()>0.5d0)then
  ranx_pm=ranx() 
 else
  ranx_pm=-ranx() 
 endif

 END FUNCTION ranx_pm


 SUBROUTINE ranlux(rvec,lenv)
!--------------------------------------------------------------------------!
! RANLUX                                                                   !
! ------                                                                   !
! Call ranlux(rvec,lenv) returns a vector rvec of lenv 32-bit random       !
! floating point numbers between zero and one (not including the endpoints)!
!--------------------------------------------------------------------------!

!--------------------------------------------------------------------------!
!                               NOTES                                      !
!                               =====                                      !
!  Subtract-and-borrow random number generator proposed by Marsaglia and   !
!  Zaman, implemented by F. James with the name RCARRY in 1991, and later  !
!  improved by Martin Luescher in 1993 to produce "Luxury Pseudorandom     !
!  Numbers".  Fortran 77 coded by F. James, 1993                           !
!                                                                          !
!  References:                                                             !
!  M. Luscher, Computer Physics Communications  79 (1994) 100              !
!  F. James, Computer Physics Communications 79 (1994) 111                 !
!                                                                          !
!  As implemented in CASINO, this generator is based on a Fortran 90       !
!  version written by Alan Miller (alan@mel.dms.csiro.au).                 !
!                                                                          !
!   LUXURY LEVELS.                                                         !
!   ------ ------      The available luxury levels are:                    !
!                                                                          !
!  level 0  (p=24): equivalent to the original RCARRY of Marsaglia         !
!           and Zaman, very long period, but fails many tests.             !
!  level 1  (p=48): considerable improvement in quality over level 0,      !
!           now passes the gap test, but still fails spectral test.        !
!  level 2  (p=97): Passes all known tests, but  theoretically still       !
!           defective.                                                     !
!  level 3  (p=223):  Any theoretically possible correlations have         !
!           very small chance of being observed [DEFAULT VALUE].           !
!  level 4  (p=389): highest possible luxury, all 24 bits chaotic.         !
!                                                                          !
!  Luxury Level     0   1   2   3     4                                    !
!    ndskip        /0, 24, 73, 199, 365/                                   !
! Corresponds to p=24  48  97  223  389                                    !
!     time factor   1   2   3    6   10   on slow workstation              !
!                   1 1.5   2    3    5   on fast mainframe                !
!                   1 1.5 2.5    5  8.5   on PC using LF90                 !
!                                                                          !
!--------------------------------------------------------------------------!

 IMPLICIT NONE
 INTEGER,INTENT(in) :: lenv
 REAL,INTENT(out) :: rvec(lenv)
 INTEGER ivec,isk
 REAL uni

 if(notyet)call errstop('RANLUX','RANLUX called without initialization.')

! The main generator: "Subtract-with-borrow", as proposed by Marsaglia and
! Zaman, Florida State University, March, 1989

 do ivec=1,lenv

  uni=seeds(j24)-seeds(i24)-carry
  if(uni<0.)then
   uni=uni+1.0
   carry=twom24
  else
   carry=0.
  endif
  seeds(i24)=uni
  i24=next(i24)
  j24=next(j24)
  rvec(ivec)=uni

! Small numbers (with less than 12 "significant" bits) are "padded".
  if(uni<twom12)then
   rvec(ivec)=rvec(ivec)+twom24*seeds(j24)

! And zero is forbidden in case someone takes a logarithm.
   if(rvec(ivec)==0.)rvec(ivec)=twom24*twom24
  endif

! Skipping to luxury. As proposed by Martin Luscher.
  in24=in24+1

  if(in24==24)then
   in24=0
   kount=kount+nskip
   do isk=1,nskip
    uni=seeds(j24)-seeds(i24)-carry
    if(uni<0.)then
     uni=uni+1.0
     carry=twom24
    else
     carry=0.
    endif
    seeds(i24)=uni
    i24=next(i24)
    j24=next(j24)
   enddo ! isk
  endif

 enddo ! ivec

 kount=kount+lenv
 if(kount>=igiga)then
  mkount=mkount+1
  kount=kount-igiga
 endif

 END SUBROUTINE ranlux


 SUBROUTINE ranlux_init(lux,ins,k1,k2,use_default,print_ran_init)
!------------------------------------------------------------!
! Subroutine to initialize from one or three integers.       !
! Call ranlux_init(LUX,INT,K1,K2) initializes the generator  !
! from one 32-bit integer INT and sets luxury level LUX      !
! which is integer between zero and MAXLEV. The k1 and k2    !
! parameters should be set to zero unless restarting at a    !
! break point given by output of RLUXAT (see RLUXAT).        !
!------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: lux,ins,k1,k2
 LOGICAL,INTENT(in) :: use_default,print_ran_init
 INTEGER i,iouter,iseeds(24),isk,k,inner,izip,izip2,p
 REAL uni

 luxlev=lux
 nskip=ndskip(luxlev)
 in24=0

 if(am_master.and.print_ran_init)then
  write(o,*)'Initialize random number generator'
  write(o,*)'=================================='
  write(o,*)'Generator                                 :  RANLUX'
  write(o,'(1x,a,i2)')'RANLUX luxury level                       : ',luxlev
  p=nskip+24
  if(p<10)then
   write(o,'(1x,a,i2)')'p value                                   : ',p
  elseif(p>10.and.p<100)then
   write(o,'(1x,a,i3)')'p value                                   : ',p
  elseif(p>100.and.p<1000)then
   write(o,'(1x,a,i4)')'p value                                   : ',p
  elseif(p>1000.and.p<10000)then
   write(o,'(1x,a,i5)')'p value                                   : ',p
  endif
  write(o,*)
 endif

 if(ins<0)call errstop('RANLUX_INIT', &
  &'Illegal initialization : negative input seed.')

 if(ins>0)then
  jseed=ins
 else
  jseed=randomseed
 endif
 if(print_ran_init.and.use_default)then
  if(nnodes==1)then
   write(o,'(1x,a,i10)')'Initialized from seed                     :  ',jseed
   write(o,*)
  elseif(am_master)then
   write(o,'(1x,a,i10)')'Initialized from default seed on node 0   :  ',jseed
   write(o,*)
  endif
 endif

 notyet=.false.
 twom24=1.
 do i=1,24
  twom24=twom24*0.5
  k=jseed/53668
  jseed=40014*(jseed-k*53668)-k*12211
  if(jseed<0)jseed=jseed+icons
  iseeds(i)=mod(jseed,itwo24)
 enddo

 twom12=twom24*4096.
 do i=1,24
  seeds(i)=real(iseeds(i))*twom24
  next(i)=i-1
 enddo

 next(1)=24
 i24=24
 j24=10
 carry=0.
 if(seeds(24)==0.)carry=twom24

! If restarting at a break point, skip K1 + IGIGA*K2
! Note that this is the number of numbers delivered to
! the user PLUS the number skipped (if luxury > 0).
 kount=k1
 mkount=k2
 if(k1+k2/=0)then

  do iouter=1,k2+1
   inner=igiga
   if(iouter==k2+1)inner=k1
   do isk=1,inner
    uni=seeds(j24)-seeds(i24)-carry
    if(uni<0.)then
     uni=uni+1.0
     carry=twom24
    else
     carry=0.
    endif
    seeds(i24)=uni
    i24=next(i24)
    j24=next(j24)
   enddo ! isk
  enddo ! iouter

! Get the right value of IN24 by direct calculation.
  in24=mod(kount,nskip+24)
  if(mkount>0)then
   izip=mod(igiga,nskip+24)
   izip2=mkount*izip+in24
   in24=mod(izip2,nskip+24)
  endif

! Now IN24 had better be between zero and 23 inclusive.
  if(in24>23)then
   write(o,'(a/a,3i11,a,i5)') &
    &'  Error in restarting with RLUXGO:','  The values',ins,&
    &k1,k2,' cannot occur at luxury level',luxlev
    in24=0
  endif
 endif

 END SUBROUTINE ranlux_init


END MODULE random_numbers

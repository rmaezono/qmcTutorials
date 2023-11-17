!----------------------------------------------------------------------------!
! TWISTOFFSET_PWCSF                                                          !
!                                                                            !
! MDT 11.2011 (based on NDD's 'twistoffset_castep').                         !
!                                                                            !
! This program applies a random offset to the list of k vectors listed       !
! in a file linked to standard input, and is used by the twistav_pwscf       !
! script when doing twist-averaged VMC/DMC runs with CASINO and PWSCF.       !
! The offset data are written to a file called 'pwscf_kpoints.out'.          !
!                                                                            !
! CHANGES                                                                    !
! =======                                                                    !
! None yet.                                                                  !
!----------------------------------------------------------------------------!


MODULE dsp
!------------------------!
! Double precision type. !
!------------------------!
 IMPLICIT NONE
 INTEGER,PARAMETER :: dp=KIND(1.d0)
END MODULE dsp


MODULE rng
!----------------------------------------------------------------------------!
! Random number generator, using the method suggested by D.E. Knuth in       !
! Seminumerical Algorithms (vol 2 of The Art of Computer Programming).       !
! The method is based on lagged Fibonacci sequences with subtraction.        !
!----------------------------------------------------------------------------!
 USE dsp,ONLY : dp
 IMPLICIT NONE
 PRIVATE
 PUBLIC ranx
 INTEGER,PARAMETER :: KK=100,LL=37 ! Leave these.
 REAL(dp) :: ranstate(kk)  ! Determines output of gen_ran_array.


CONTAINS


 REAL(dp) FUNCTION ranx()
!------------------------------------------------------------------------------!
! Return a random number uniformly distributed in [0,1).                       !
! Uses M. Luescher's suggestion: generate 1009 random numbers at a time using  !
! Knuth's algorithm, but only use the first 100.                               !
!------------------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,PARAMETER :: seed=310952  ! Random seed, betw. 0 & 2^30-3.
  INTEGER,PARAMETER :: Nran=1009,Nkeep=100 ! See comment on p. 188 of Knuth.
  INTEGER,SAVE :: i=-1
  REAL(dp),SAVE :: ran_array(Nran)
  if(i==-1)then
   call init_rng(seed) ! Initialize the RNG.
   i=Nkeep
  endif ! First call.
  if(i==Nkeep)then
   call gen_ran_array(ran_array,Nran) ! Generate a new array of random nos.
   i=0
  endif ! i=Nkeep
  i=i+1
  ranx=ran_array(i)
 END FUNCTION ranx


 SUBROUTINE gen_ran_array(ran_array,N)
!---------------------------------------------------------------!
! Generate an array of N random numbers: see Knuth's ran_array. !
!---------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: N
  REAL(dp),INTENT(out) :: ran_array(N)
  INTEGER j
  ran_array(1:KK)=ranstate(1:KK)
  do j=KK+1,N
   ran_array(j)=mod(ran_array(j-KK)+ran_array(j-LL),1.d0)
  enddo ! j
  do j=1,LL
   ranstate(j)=mod(ran_array(N+j-KK)+ran_array(N+j-LL),1.d0)
  enddo ! j
  do j=LL+1,KK
   ranstate(j)=mod(ran_array(N+j-KK)+ranstate(j-LL),1.d0)
  enddo ! j
 END SUBROUTINE gen_ran_array


 SUBROUTINE init_rng(seed)
!--------------------------------------------!
! Initialize the RNG: see Knuth's ran_start. !
!--------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: seed
  INTEGER j,s,t,sseed
  INTEGER,PARAMETER :: MM=2**30,TT=70
  REAL(dp) ss,x(KK+KK-1)
  REAL(dp),PARAMETER :: ULP=1.d0/2.d0**52,ULP2=2.d0*ULP
  if(seed<0)then
   sseed=MM-1-mod(-1-seed,MM)
  else
   sseed=mod(seed,MM)
  endif ! seed<0
  ss=ULP2*dble(sseed+2)
  do j=1,KK
   x(j)=ss
   ss=ss+ss
   if(ss>=1.d0)ss=ss-1.d0+ULP2
  enddo ! j
  x(2)=x(2)+ULP
  s=sseed
  t=TT-1
  do
   do j=KK,2,-1
    x(j+j-1)=x(j)
    x(j+j-2)=0.d0
   enddo ! j
   do j=KK+KK-1,KK+1,-1
    x(j-(KK-LL))=mod(x(j-(KK-LL))+x(j),1.d0)
    x(j-KK)=mod(x(j-KK)+x(j),1.d0)
   enddo ! j
   if(mod(s,2)==1)then
    do j=KK,1,-1
     x(j+1)=x(j)
    enddo ! j
    x(1)=x(KK+1)
    x(LL+1)=mod(x(LL+1)+x(KK+1),1.d0)
   endif ! s odd
   if(s/=0)then
    s=s/2
   else
    t=t-1
   endif ! s/=0
   if(t<=0)exit
  enddo
  ranstate(1+KK-LL:KK)=x(1:LL)
  ranstate(1:KK-LL)=x(LL+1:KK)
  do j=1,10
   call gen_ran_array(x,KK+KK-1)
  enddo ! j
 END SUBROUTINE init_rng


END MODULE rng


MODULE utils
!--------------------------!
! Miscellaneous utilities. !
!--------------------------!
 USE dsp
 IMPLICIT NONE


CONTAINS


 SUBROUTINE errstop(msg)
!------------------------------------!
! Stop and display an error message. !
!------------------------------------!
  IMPLICIT NONE
  CHARACTER(*),INTENT(in) :: msg
  write(6,'(a)')'TWISTOFFSET ERROR: '//trim(msg)
  stop
 END SUBROUTINE errstop


END MODULE utils


PROGRAM twistoffset
!---------------------------!
! Main program starts here. !
!---------------------------!
 USE dsp,ONLY : dp
 USE rng,ONLY : ranx
 USE utils,ONLY : errstop
 IMPLICIT NONE
 INTEGER ierr,ik,nk,twist,i
 REAL(dp) koffset(3)
! REAL(dp),PARAMETER :: tol=1.d-8
 REAL(dp),ALLOCATABLE :: kpoint(:,:),weight(:)

! Call random number generator an appropriate number of times, so we don't
! use the same offset at each twist.
! Find k offset.  Each component is uniform in [0,1].
 read(5,*,iostat=ierr)twist
 if(ierr/=0)call errstop('Error reading twist number.')
 koffset=0.d0
 do i=1,twist
  koffset=(/ranx(),ranx(),ranx()/)
 enddo

! Read in number of k points.
 read(5,*,iostat=ierr)nk
 if(ierr/=0)call errstop('Error reading number of k points.')
 if(nk<1)call errstop('Incorrect value of nk read from pwscf_kpoints.in&
  & via standard input.')

 allocate(kpoint(3,nk),weight(nk),stat=ierr)
 if(ierr/=0)call errstop('Error allocating k vectors.')

! Read in k points
 do ik=1,nk
  read(5,*,end=1,err=1)kpoint(1:3,ik),weight(ik)
 enddo

! if(any(abs(weight(1:nk)-1.d0/real(nk,dp))>tol))call errstop('The k-point &
!  &weights should all be the same, and should add up to 1.')
!
! Except they don't in PWSCF (check?) so skip this test.. :-)

! Open output file
 open(unit=8,file='pwscf_kpoints.out',status='replace',iostat=ierr)
 if(ierr/=0)call errstop('Cannot open pwscf_kpoints.out file.')

! Apply offset to k points.  Put them all in [0,1].  Write out new k points.
 write(8,*)twist+1
 write(8,*)nk
 do ik=1,nk
!  write(8,'(4(es23.16,1x))')modulo(kpoint(1:3,ik)+koffset,1.d0),1.d0/dble(nk)
  write(8,'(4(es23.16,1x))')modulo(kpoint(1:3,ik)+koffset,1.d0),weight(ik)
 enddo ! ik

 stop

1 call errstop('Error reading k points from pwscf_kpoints.in via standard&
   & input.')

END PROGRAM twistoffset

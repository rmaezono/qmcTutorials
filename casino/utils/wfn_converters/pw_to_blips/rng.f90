MODULE rng
!----------------------------------------------------------------------------!
! Random number generator, using the method suggested by D.E. Knuth in       !
! Seminumerical Algorithms (vol 2 of The Art of Computer Programming).       !
! The method is based on lagged Fibonacci sequences with subtraction.        !
!----------------------------------------------------------------------------!
 IMPLICIT NONE
 PRIVATE
 PUBLIC ranx,init_rng
 INTEGER,PARAMETER :: KK=100,LL=37 ! Leave these.
 DOUBLE PRECISION :: ranstate(kk)  ! Determines output of gen_ran_array.

 INTEGER,PARAMETER :: default_seed=310952  ! Random seed, betw. 0 & 2^30-3.
 INTEGER,PARAMETER :: Nran=1009,Nkeep=100 ! See comment on p. 188 of Knuth.
 INTEGER :: ran_array_idx=-1
 DOUBLE PRECISION :: ran_array(Nran)


CONTAINS


 DOUBLE PRECISION FUNCTION ranx()
!------------------------------------------------------------------------------!
! Return a random number uniformly distributed in [0,1).                       !
! Uses M. Luescher's suggestion: generate 1009 random numbers at a time using  !
! Knuth's algorithm, but only use the first 100.                               !
!------------------------------------------------------------------------------!
  IMPLICIT NONE
  if(ran_array_idx==-1)then
   call init_rng(default_seed) ! Initialize the RNG.
  endif ! First call.
  if(ran_array_idx==Nkeep)then
   call gen_ran_array(ran_array,Nran) ! Generate a new array of random nos.
   ran_array_idx=0
  endif ! i=Nkeep
  ran_array_idx=ran_array_idx+1
  ranx=ran_array(ran_array_idx)
 END FUNCTION ranx


 SUBROUTINE gen_ran_array(ran_array,N)
!---------------------------------------------------------------!
! Generate an array of N random numbers: see Knuth's ran_array. !
!---------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: N
  DOUBLE PRECISION,INTENT(out) :: ran_array(N)
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
  DOUBLE PRECISION ss,x(KK+KK-1)
  DOUBLE PRECISION,PARAMETER :: ULP=1.d0/2.d0**52,ULP2=2.d0*ULP
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
  ran_array_idx=Nkeep
 END SUBROUTINE init_rng


END MODULE rng

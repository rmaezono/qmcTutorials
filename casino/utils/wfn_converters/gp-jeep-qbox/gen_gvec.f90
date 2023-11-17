MODULE gen_gvec_module
IMPLICIT NONE
CONTAINS


SUBROUTINE gen_gvec(ecut,alat,blat,clat,ng,gvec_int,bCount_only)
!-------------------------------------------------------------------!
! This subroutine generates a list of gvectors equivalent to those  !
! generated in the JEEP code (like the genbasis utility)            !
!-------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(inout) :: ng,gvec_int(3,ng)
 DOUBLE PRECISION,INTENT(in) :: ecut,alat,blat,clat
 LOGICAL,INTENT(in) :: bCount_only
 INTEGER,ALLOCATABLE :: index(:),gvec_temp(:,:)
 DOUBLE PRECISION,ALLOCATABLE :: esort(:)
 INTEGER kxmax,kymax,kzmax,nx,ny,nz,i,j,k
 DOUBLE PRECISION sqrt_ecut,twopi,fac1,fac2,fac3,e,d0,d1,d2,&
  &kn0,kn1,kn2,dnorm,sknorm,sp

 twopi=8.d0*atan(1.d0)

 sqrt_ecut=sqrt(ecut)
 kxmax=int(alat*sqrt_ecut/twopi)
 kymax=int(blat*sqrt_ecut/twopi)
 kzmax=int(clat*sqrt_ecut/twopi)

 nx=2*kxmax+2
 do
  nx=nx+2
  if(factorizable(nx))exit
 enddo
 ny=2*kymax+2
 do
  ny=ny+2
  if(factorizable(ny))exit
 enddo
 nz=2*kzmax+2
 do
  nz=nz+2
  if(factorizable(nz))exit
 enddo

 fac1=twopi/alat ; fac2=twopi/blat ; fac3=twopi/clat

! Define norms used to split degeneracy within a shell
 d0=(2*kymax+1)*(2*kzmax+1)
 d1=(2*kzmax+1)
 d2=1
 dnorm=sqrt(d0*d0+d1*d1+d2*d2)
 d0=d0/dnorm; d1=d1/dnorm; d2=d2/dnorm

! Count how many g-vectors to keep
! sp is used to only keep half the Kx=0 plane.
 ng=1
 do i=0,kxmax
  do j=-kymax,kymax
   do k=-kzmax,kzmax
    sp=i*d0+j*d1+k*d2
    if(sp>0)then
     e=0.5d0*((fac1*dble(i))**2+(fac2*dble(j))**2+(fac3*dble(k))**2)
     if(e<0.5d0*ecut)ng=ng+1
    endif
   enddo
  enddo
 enddo

 if(bCount_only)return

! Allocate the set of g-vectors
 if(allocated(esort))deallocate(esort)
 allocate(esort(ng))

 ng=1
 gvec_int(1,ng)=0 ; gvec_int(2,ng)=0 ; gvec_int(3,ng)=0
 esort(ng)=0.d0
 do i=0,kxmax
  do j=-kymax,kymax
   do k=-kzmax,kzmax
    e=0.5d0*((fac1*dble(i))**2+(fac2*dble(j))**2+(fac3*dble(k))**2)
    sp=i*d0+j*d1+k*d2
    if(sp>0)then
     if(e<0.5d0*ecut)then
      ng=ng+1
      gvec_int(1,ng)=i ; gvec_int(2,ng)=j ; gvec_int(3,ng)=k
      sknorm=sqrt (dble(i)**2+dble(j)**2+dble(k)**2)
      kn0=dble(i)/sknorm
      kn1=dble(j)/sknorm
      kn2=dble(k)/sknorm
      esort(ng)=1000000*e-(kn0*d0+kn1*d1+kn2*d2)
     endif
    endif
   enddo
  enddo
 enddo

! Sort the gvectors according to Francois's criteria

! Insertion sort from Num. Recip. (N^2 also) - 28 seconds in test
!  do j=2,ng
!   etemp=esort(j)
!   gtemp=gvec_int(:,j)
!   do i=j-1,1,-1
!    if(esort(i)<etemp)exit
!    esort(i+1)=esort(i)
!    gvec_int(:,i+1)=gvec_int(:,i)
!   enddo
!   esort(i+1)=etemp
!   gvec_int(:,i+1)=gtemp
!  enddo

! Using a canned qsort routine - 1.5 seconds in test
 allocate(index(ng),gvec_temp(3,ng))
 call sortrx(ng,esort,index)
 do i=1,ng
  gvec_temp(:,i)=gvec_int(:,index(i))
 enddo
 gvec_int=gvec_temp
 deallocate(index,gvec_temp)

END SUBROUTINE gen_gvec


LOGICAL FUNCTION factorizable(i)
 IMPLICIT NONE
 INTEGER,INTENT(in) :: i
 INTEGER j

 j=i
 factorizable=.false.
 do
  if(mod(j,2)==0)then
   j=j/2
  else
   exit
  endif
 enddo
 do
  if(mod(j,3)==0)then
   j=j/3
  else
   exit
  endif
 enddo
 do
  if(mod(j,5)==0)then
   j=j/5
  else
   exit
  endif
 enddo
 if(j==1)factorizable=.true.

END FUNCTION factorizable


SUBROUTINE SORTRX(N,DATA,INDEX)
!----------------------------------------------------------------------!
!                                                                      !
!     SORTRX -- SORT, Real input, indeX output                         !
!                                                                      !
!                                                                      !
!     Input:  N     INTEGER                                            !
!             DATA  REAL*8                                             !
!                                                                      !
!     Output: INDEX INTEGER (DIMENSION N)                              !
!                                                                      !
! This routine performs an in-memory sort of the first N elements of   !
! array DATA, returning into array INDEX the indices of elements of    !
! DATA arranged in ascending order.  Thus,                             !
!                                                                      !
!    DATA(INDEX(1)) will be the smallest number in array DATA;         !
!    DATA(INDEX(N)) will be the largest number in DATA.                !
!                                                                      !
! The original data is not physically rearranged.  The original order  !
! of equal input values is not necessarily preserved.                  !
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
! SORTRX uses a hybrid QuickSort algorithm, based on several           !
! suggestions in Knuth, Volume 3, Section 5.2.2.  In particular, the   !
! "pivot key" [my term] for dividing each subsequence is chosen to be  !
! the median of the first, last, and middle values of the subsequence; !
! and the QuickSort is cut off when a subsequence has 9 or fewer       !
! elements, and a straight insertion sort of the entire array is done  !
! at the end.  The result is comparable to a pure insertion sort for   !
! very short arrays, and very fast for very large arrays (of order 12  !
! micro-sec/element on the 3081K for arrays of 10K elements).  It is   !
! also not subject to the poor performance of the pure QuickSort on    !
! partially ordered data.                                              !
!                                                                      !
! Created:  15 Jul 1986  Len Moss                                      !
!                                                                      !
!----------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: N
 DOUBLE PRECISION,INTENT(in) :: DATA(N)
 INTEGER,INTENT(out) :: INDEX(N)
 INTEGER LSTK(31),RSTK(31),ISTK
 INTEGER L,R,I,J,P,INDEXP,INDEXT
 DOUBLE PRECISION DATAP
! Quit QuickSort-ing when a subsequence contains M or fewer
! elements and finish off at end with straight insertion sort.
! According to Knuth, V.3, the optimum value of M is around 9.
 INTEGER,PARAMETER :: M=9

! Make initial guess for INDEX
 do I=1,N
  index(I)=I
 enddo ! I

! If array is short, skip QuickSort and go directly to the straight
! insertion sort.
 if(N<=M)goto 900

! The "Qn:"s correspond roughly to steps in Algorithm Q, Knuth, V.3,
! PP.116-117, modified to select the median of the first, last, and middle
! elements as the "pivot key" (in Knuth's notation, "K").  Also modified to
! leave data in place and produce an INDEX array.  To simplify comments, let
! DATA[I]=DATA(INDEX(I)).

! Q1: Initialize
 ISTK=0
 L=1
 R=N

200 continue

! Q2: Sort the subsequence DATA[L]..DATA[R].
!
! At this point, DATA[l] <= DATA[m] <= DATA[r] for all l < L,
! r > R, and L <= m <= R.  (First time through, there is no
! DATA for l < L or r > R.)

 I=L
 J=R

! Q2.5: Select pivot key
!
! Let the pivot, P, be the midpoint of this subsequence,
! P=(L+R)/2; then rearrange INDEX(L), INDEX(P), and INDEX(R)
! so the corresponding DATA values are in increasing order.
! The pivot key, DATAP, is then DATA[P].

 P=(L+R)/2
 INDEXP=index(P)
 DATAP=data(INDEXP)

 if(data(index(L))>DATAP)then
  index(P)=index(L)
  index(L)=INDEXP
  INDEXP=index(P)
  DATAP=data(INDEXP)
 endif

 if(DATAP>data(index(R)))then
  if(data(index(L))>data(index(R)))then
   index(P)=index(L)
   index(L)=index(R)
  else
   index(P)=index(R)
  endif
  index(R)=INDEXP
  INDEXP=index(P)
  DATAP=data(INDEXP)
 endif

! Now we swap values between the right and left sides and/or
! move DATAP until all smaller values are on the left and all
! larger values are on the right.  Neither the left or right
! side will be internally ordered yet; however, DATAP will be
! in its final position.

300 continue

! Q3: Search for datum on left >= DATAP
!
! At this point, DATA[L] <= DATAP.  We can therefore start scanning
! up from L, looking for a value >= DATAP (this scan is guaranteed
! to terminate since we initially placed DATAP near the middle of
! the subsequence).

 I=I+1
 if(data(index(I))<DATAP)goto 300

400 continue

! Q4: Search for datum on right <= DATAP
!
! At this point, DATA[R] >= DATAP.  We can therefore start scanning
! down from R, looking for a value <= DATAP (this scan is guaranteed
! to terminate since we initially placed DATAP near the middle of
! the subsequence).

 J=J-1
 if(data(index(J))>DATAP)goto 400

! Q5: Have the two scans collided?

 if(I<J)then

! Q6: No, interchange DATA[I] <--> DATA[J] and continue

  INDEXT=index(I)
  index(I)=index(J)
  index(J)=INDEXT
  goto 300
 else

! Q7: Yes, select next subsequence to sort
!
! At this point, I >= J and DATA[l] <= DATA[I] == DATAP <= DATA[r],
! for all L <= l < I and J < r <= R.  If both subsequences are
! more than M elements long, push the longer one on the stack and
! go back to QuickSort the shorter; if only one is more than M
! elements long, go back and QuickSort it; otherwise, pop a
! subsequence off the stack and QuickSort it.

  if(R-J>=I-L.and.I-L>M)then
   ISTK=ISTK+1
   LSTK(ISTK)=J+1
   RSTK(ISTK)=R
   R=I-1
  elseif(I-L>R-J.and.R-J>M)then
   ISTK=ISTK+1
   LSTK(ISTK)=L
   RSTK(ISTK)=I-1
   L=J+1
  elseif(R-J>M)then
   L=J+1
  elseif(I-L>M)then
   R=I-1
  else
! Q8: Pop the stack, or terminate QuickSort if empty
   if(ISTK<1)goto 900
   L=LSTK(ISTK)
   R=RSTK(ISTK)
   ISTK=ISTK-1
  endif
  goto 200
 endif

900 continue

! Q9: Straight Insertion sort
 do I=2,N
  if(data(index(I-1))>data(index(I)))then
   INDEXP=index(I)
   DATAP=data(INDEXP)
   P=I-1
920 continue
   index(P+1)=index(P)
   P=P-1
   if(P>0)then
    if(data(index(P))>DATAP)goto 920
   endif
   index(P+1)=INDEXP
  endif
 enddo

END SUBROUTINE SORTRX


END MODULE gen_gvec_module

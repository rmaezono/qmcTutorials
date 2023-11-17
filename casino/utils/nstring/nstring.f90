PROGRAM nstring
 !-------------------------------------------------------------------------!
 ! NSTRING                                                                 !
 ! -------                                                                 !
 ! Generate integer number sequence from a to b step c.'                   !
 ! Useful for strings required in correlation.data files if you're using   !
 ! Type 1 labelling of the sets of atoms in chi and f functions.           !
 !                                                                         !
 ! MDT 5.2002                                                              !
 !-------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER a,b,c,i
 INTEGER,PARAMETER :: maxn=10000
 INTEGER d(maxn)
 write(6,*)'NSTRING'
 write(6,*)'-------'
 write(6,*)
 write(6,*)'Generate integer number sequence from a to b step c.'
 write(6,*)'a?'
 do
  read(5,*)a
  if(a>0.and.a<maxn)exit
  write(6,*)
  write(6,*)'a must be +ve and less than ',maxn
  write(6,*)
 enddo
 do
  write(6,*)'b?'
  read(5,*)b
  if(b>0.and.b<maxn)exit
  write(6,*)
  write(6,*)'b must be +ve and less than ',maxn
  write(6,*)
 enddo
 do
  write(6,*)'c?'
  read(5,*)c
  if(c>0.and.b<a)then
   write(6,*)
   write(6,*)'c is +ve but b < a)'
   write(6,*)
  elseif(c<0.and.b>a)then
   write(6,*)
   write(6,*)'c is -ve but b > a)'
   write(6,*)
  elseif(c==0)then
   write(6,*)
   write(6,*)'c should not be zero'
   write(6,*)
  else
   exit
  endif
 enddo
 write(6,*)

 do i=a,b,c
  d(i)=i
 enddo
 write(6,'(3x,32767(a,1x))')(trim(i2s(d(i))),i=a,b,c)


CONTAINS


 CHARACTER(20) FUNCTION i2s(n)
!-----------------------------------------------------------------------!
! I2S                                                                   !
! ===                                                                   !
! Convert integers to left justified strings that can be printed in the !
! middle of a sentence without introducing large amounts of white space.!
!                                                                       !
! Calling routine is intended to include something like:                !
! INTEGER i                                                             !
! i=12                                                                  !
! write(6,*)'Integer number ',trim(i2s(i)),' with words at the end.'    !
!-----------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: n
  INTEGER i,j
  CHARACTER tmp,sign

  if(n==0)then
   i2s='0' ; return
  endif
  sign=' ' ; if(n<0)sign='-'

  do i=1,len(i2s)
   i2s(i:i)=' '
  enddo

  i=abs(n)
  do j=1,len(i2s)
   if(i==0)exit
   i2s(j:j)=achar(ichar('0')+mod(i,10))
   i=i/10
  enddo

  i=1 ; j=len_trim(i2s)
  do
   if(i>=j)exit
   tmp=i2s(j:j)
   i2s(j:j)=i2s(i:i)
   i2s(i:i)=tmp
   i=i+1
   j=j-1
  enddo

  i2s=trim(sign)//i2s

 END FUNCTION i2s


END PROGRAM nstring

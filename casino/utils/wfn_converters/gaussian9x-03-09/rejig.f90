MODULE rejig
!---------------------------------------------------------------!
! Given n values in d, this routine sorts them into descending  !
! order and reorders the contents of reg accordingly.  The      !
! method is straight insertion.                                 !
!---------------------------------------------------------------!

CONTAINS

 SUBROUTINE numsrt(d,n,np,reg)
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n,np
 INTEGER,INTENT(inout) :: d(np)
 INTEGER,INTENT(inout),OPTIONAL :: reg(np)
 INTEGER i,j,k,p,iq

 do i=1,n-1
  k=i
  p=d(i)
  do j=i+1,n
   if(d(j)>=p)then
    k=j
    p=d(j)
   endif
  enddo
  if(k/=i)then
! We have swapped two elements
   d(k)=d(i)
   d(i)=p
   if(PRESENT(reg))then
! Apply same change to elements of reg
    iq=reg(k)
    reg(k)=reg(i)
    reg(i)=iq
   endif
  endif
 enddo

 END SUBROUTINE numsrt

END MODULE rejig

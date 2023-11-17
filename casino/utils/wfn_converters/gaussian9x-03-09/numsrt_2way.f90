SUBROUTINE numsrt_2way(d,n,np,reg,dirn)
!----------------------------------------------------------------!
! Given n values in d, this routine sorts them into descending   !
! or ascending order (depending on whether dirn is 'd' or 'a')   !
! and reorders the contents of reg accordingly. The method       !
! is straight insertion.                                         !
!----------------------------------------------------------------!
 INTEGER,INTENT(in) :: n,np
 INTEGER,INTENT(inout) :: d(np),reg(np)
 CHARACTER(1),INTENT(in) :: dirn
 INTEGER i,j,k,p,iq

 if(dirn=='d'.or.dirn=='D')then
! Sort into descending order...
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
! Apply same change to elements of reg
    iq=reg(k)
    reg(k)=reg(i)
    reg(i)=iq
   endif
  enddo
 else
! Sort into ascending order...
  do i=1,n-1
   k=i
   p=d(i)
   do j=i+1,n
    if(d(j)<=p)then
     k=j
     p=d(j)
    endif
   enddo
   if(k/=i)then
! We have swapped two elements
    d(k)=d(i)
    d(i)=p
! Apply same change to elements of reg
    iq=reg(k)
    reg(k)=reg(i)
    reg(i)=iq
   endif
  enddo
 endif

END SUBROUTINE numsrt_2way

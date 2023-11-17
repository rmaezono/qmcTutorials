 SUBROUTINE mcm(is1,is2,is3,is,isj1,isj2,isj3)
  USE numbers
  IMPLICIT REAL(float) (a-h,o-z)
  IMPLICIT INTEGER (i-n)
  is=is1
  n=is2
1 do l=1,n
   m=is*l
   if(mod(m,n)==0)exit
  enddo
  is=m
  if(n==is3)goto 5
  n=is3
  goto 1
5 isj1=m/is1
  isj2=m/is2
  isj3=m/n
 END SUBROUTINE mcm

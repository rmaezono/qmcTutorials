 SUBROUTINE matout(a,nr,nc)
  USE numbers
  USE parinf_module
  IMPLICIT REAL(float) (a-h,o-z)
  IMPLICIT INTEGER (i-n)
  DIMENSION a(nr,*)
  do m=1,nc,10
   k=min(m+9,nc)
   write(iout,40)(j,j=m,k)
   do i=1,nr
    write(iout,50)i,(a(i,j),j=m,k)
   enddo
  enddo
40 format(/7x,10(5x,i3,4x)/)
50 format(i4,3x,1p,10e12.4)
 END SUBROUTINE matout

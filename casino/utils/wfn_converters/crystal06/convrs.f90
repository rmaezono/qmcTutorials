 SUBROUTINE convrs(I,J,Z,APOLE,IPOLE,JPOLE)
  USE numbers
  USE lmaxxx
  IMPLICIT REAL(FLOAT) (A-H,O-Z)
  IMPLICIT INTEGER (I-N)
  DIMENSION APOLE(*),IPOLE(*),JPOLE(*)
  COMMON/LOCO/LPOL,NPO

  lponpo=lpol+npo
  do ll=lpol+1,lponpo
   if(i==ipole(ll).and.j==jpole(ll))goto 1
  enddo
  npo=npo+1
  apole(lponpo+1)=z
  ipole(lponpo+1)=i
  jpole(lponpo+1)=j
  return
1 bilbo=apole(ll)+z
  if(abs(bilbo)<(1e-9_float))then
   do l=ll+1,lponpo
    apole(l-1)=apole(l)
    ipole(l-1)=ipole(l)
    jpole(l-1)=jpole(l)
   enddo
   npo=npo-1
  else
   apole(ll)=bilbo
  endif
  return
 END SUBROUTINE convrs

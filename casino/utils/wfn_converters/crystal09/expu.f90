 SUBROUTINE expu(j1,j2,j3)
  USE numbers
  USE parame_module
  USE parinf_module
  USE gvect_module
  USE retic_module
  USE expo_module
  IMPLICIT REAL(float) (a-h,o-z)
  IMPLICIT INTEGER (i-n)
  j1vrs=j1*inf(89)
  j2vrs=j2*inf(90)
  j3vrs=j3*inf(91)
  ism=inf(88)
  iq=ism*8192
  ex(1,1)=1._float
  ex(2,1)=0._float
  do mg=2,inf(28),2
   jl=mod(j1vrs*lg(1,mg)+j2vrs*lg(2,mg)+j3vrs*lg(3,mg)+iq,ism)
   vrs=cossma(jl+1)
   wrs=sinsma(jl+1)
   ex(1,mg)=vrs
   ex(2,mg)=wrs
   ex(1,mg+1)=vrs
   ex(2,mg+1)=-wrs
  enddo
 END SUBROUTINE expu

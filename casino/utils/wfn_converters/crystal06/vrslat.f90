 SUBROUTINE vrslat
  USE numbers
  USE parame_module
  USE parinf_module
  USE retic_module
  IMPLICIT REAL(float) (a-h,o-z)
  IMPLICIT INTEGER (i-n)
  latsum=0
  do k=1,nkf
  latvrs(k)=0
  if(mod(jj(1,k)*2,is1)/=0)cycle
  if(mod(jj(2,k)*2,is2)/=0)cycle
  if(mod(jj(3,k)*2,is3)==0)then
   latvrs(k)=1
   latsum=latsum+1
  endif
  enddo
  if(latsum==nkf)then
   inf(81)=1
  else
   inf(81)=0
  endif
  call mcm(is1,is2,is3,ism,isj1,isj2,isj3)
  if(ism>lim086)call errnic(0,ism,'vrslat',&
   &'ISM for Monkhorst net too large - increase LIM086')
  inf(88)=ism
  inf(89)=isj1
  inf(90)=isj2
  inf(91)=isj3
  vrs=par(34)/ism
  cossma(1)=1._float
  sinsma(1)=0._float
  do k=1,ism-1
   wrs=k*vrs
   cossma(k+1)=cos(wrs)
   sinsma(k+1)=sin(wrs)
  enddo
 END SUBROUTINE vrslat

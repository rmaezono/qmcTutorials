 SUBROUTINE estroe(a,ar,mv,nbini,nbfi,mbands)
  USE numbers
  USE parame_module
  USE parinf_module
  USE memory_screen
  USE rotmatrix
  USE basato_module
  USE expo_module
  IMPLICIT REAL(float) (a-h,o-z)
  IMPLICIT INTEGER (i-n)
  DIMENSION a(*),ar(*)
  ndf=inf(7)
  indbas=(nbini-1)*ndf
  ndf=ndf+ndf
  mvf=inf(2)
  ar(1:mbands*ndf)=0._float
  do lavrs=1,inf(24)
   i=mgnav(lavrs,mv)
   eco=ex(1,i)
   esi=ex(2,i)
   do la=nshpri(lavrs),nshpri(lavrs+1)-1
    ica=ndq(lav(la,mv))+indbas
    inf3=lat(la)*mvf+mv
    ico=ndq(la)
    do i=minz(inf3)+1,minz(inf3+1)
     vsi=tto(i)
     vco=eco*vsi
     vsi=esi*vsi
     ica1=(mmo(i)+ica)*2
     ica2=(mmom(i)+ico)*2
     do ind=nbini,nbfi
      a1=a(ica1-1)
      a2=a(ica1)
      ar(ica2-1)=a1*vco+a2*vsi+ar(ica2-1)
      ar(ica2)=a1*vsi-a2*vco+ar(ica2)
      ica1=ica1+ndf
      ica2=ica2+ndf
     enddo
    enddo
   enddo
  enddo
 END SUBROUTINE estroe


 SUBROUTINE estrof(a,ar,mv,nbini,nbfi,mbands)
  USE numbers
  USE parame_module
  USE parinf_module
  USE memory_screen
  USE rotmatrix
  USE basato_module
  USE expo_module
  IMPLICIT REAL(float) (a-h,o-z)
  IMPLICIT INTEGER (i-n)
  DIMENSION a(*),ar(*)
  ndf=inf(7)
  indbas=(nbini-1)*ndf
  ndf=ndf+ndf
  mvf=inf(2)
  ar(1:mbands*ndf)=0._float
  do lavrs=1,inf(24)
   i=mgnav(lavrs,mv)
   eco=ex(1,i)
   esi=ex(2,i)
   do la=nshpri(lavrs),nshpri(lavrs+1)-1
    ica=ndq(lav(la,mv))+indbas
    inf3=lat(la)*mvf+mv
    ico=ndq(la)
    do i=minz(inf3)+1,minz(inf3+1)
     vsi=tto(i)
     vco=eco*vsi
     vsi=esi*vsi
     ica1=(mmo(i)+ica)*2
     ica2=(mmom(i)+ico)*2
     do ind=nbini,nbfi
      a1=a(ica1-1)
      a2=a(ica1)
      ar(ica2-1)=a1*vco-a2*vsi+ar(ica2-1)
      ar(ica2)=a1*vsi+a2*vco+ar(ica2)
      ica1=ica1+ndf
      ica2=ica2+ndf
     enddo
    enddo
   enddo
  enddo
 END SUBROUTINE estrof


 SUBROUTINE estrog(a,ar,mv,nbini,nbfi,mbands)
  USE numbers
  USE parame_module
  USE parinf_module
  USE memory_screen
  USE rotmatrix
  USE basato_module
  USE expo_module
  IMPLICIT REAL(float) (a-h,o-z)
  IMPLICIT INTEGER (i-n)
  DIMENSION a(*),ar(*)
  ndf=inf(7)
  mvf=inf(2)
  indbas=(nbini-1)*ndf
  ar(1:mbands*ndf)=0._float
  do lavrs=1,inf(24)
   eco=ex(1,mgnav(lavrs,mv))
   do la=nshpri(lavrs),nshpri(lavrs+1)-1
    ica=ndq(lav(la,mv))+indbas
    ico=ndq(la)
    inf3=lat(la)*mvf+mv
    do i=minz(inf3)+1,minz(inf3+1)
     vco=tto(i)*eco
     ica1=mmo(i)+ica
     ica2=mmom(i)+ico
     do ind=nbini,nbfi
      ar(ica2)=a(ica1)*vco+ar(ica2)
      ica1=ica1+ndf
      ica2=ica2+ndf
     enddo
    enddo
   enddo
  enddo
 END SUBROUTINE estrog



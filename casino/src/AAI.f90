module slaarnaai
use dsp
use slaarnaat
use slaarnabg
use slaarnabp
use store
use slaarnaag,   only : czero,zi
use slaarnabt,   only : mxmb,mxmbc
use run_control, only : errstop
implicit none
private
public wfdet_s,wfdet_s_bf
real(dp) xyzzyaaaa1(1),xyzzyaaab1(1),xyzzyaaac1(1)
complex(dp) xyzzyaaad1(1),xyzzyaaae1(1),xyzzyaaaf1(1)
contains
subroutine wfdet_s(r,igl,lap_convert,orbital,ion,jspin,band,nk,ion_s,v&
&irtual,spin_type_in,rs,grs,ls)
implicit none
integer,intent(in) :: igl,orbital,ion,jspin,band,nk,ion_s,spin_type_in
real(dp),intent(in) :: r,lap_convert
real(dp),intent(inout) :: rs(2),grs(2),ls(2)
logical,intent(in) :: virtual
if(.not.isperiodic)then
if(igl==0)then
call xyzzyaaaa2
else
call xyzzyaaab2
endif
else
if(igl==0)then
call xyzzyaaac2
else
call xyzzyaaad2
endif
endif
contains
subroutine xyzzyaaaa2
implicit none
integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3,xyzzyaaad3,xyzzyaaae3,xyzzyaa&
&af3
real(dp) xyzzyaaag3,xyzzyaaah3
ao_m=0.d0
xyzzyaaag3=r*r
do xyzzyaaaa3=first_shell(ion),first_shell(ion+1)-1
xyzzyaaab3=shell_am(xyzzyaaaa3)
if(xyzzyaaab3>2)cycle
xyzzyaaah3=min_exponent(xyzzyaaaa3)*xyzzyaaag3
if(xyzzyaaah3>screening_tolerance)cycle
xyzzyaaac3=first_ao(xyzzyaaaa3)
xyzzyaaad3=primitive(xyzzyaaaa3)
xyzzyaaae3=primitive(xyzzyaaaa3+1)-1
do xyzzyaaaf3=xyzzyaaad3,xyzzyaaae3
xyzzyaaah3=exponent(xyzzyaaaf3)*xyzzyaaag3
ao_m(xyzzyaaac3)=ao_m(xyzzyaaac3)+c_prim(xyzzyaaaf3)*exp(-xyzzyaaah3)
enddo
enddo
if(.not.virtual)then
xyzzyaaaa1(1)=0.d0
call mxmb(rck(orbital,1:num_ao,1,spin_type_in),num_ao,1,1,ao_m,1,1,xyz&
&zyaaaa1,1,1,1,num_ao,1)
rs(1)=xyzzyaaaa1(1)*rnorm
else
rs(1)=sum(ao_m(:)*rck_ex(:,orbital,jspin))*rnorm
endif
end subroutine xyzzyaaaa2
subroutine xyzzyaaab2
implicit none
integer xyzzyaaaa4,xyzzyaaab4,xyzzyaaac4,xyzzyaaad4,xyzzyaaae4,xyzzyaa&
&af4
real(dp) xyzzyaaag4,xyzzyaaah4,xyzzyaaai4,xyzzyaaaj4,xyzzyaaak4,xyzzya&
&aal4,xyzzyaaam4,xyzzyaaan4,xyzzyaaao4
ao_m=0.d0
alap_m=0.d0
agra1_m=0.d0
xyzzyaaag4=r*r
do xyzzyaaaa4=first_shell(ion),first_shell(ion+1)-1
xyzzyaaab4=shell_am(xyzzyaaaa4)
if(xyzzyaaab4>2)cycle
xyzzyaaaj4=min_exponent(xyzzyaaaa4)
xyzzyaaai4=xyzzyaaaj4*xyzzyaaag4
if(xyzzyaaai4>screening_tolerance)cycle
xyzzyaaac4=first_ao(xyzzyaaaa4)
xyzzyaaad4=primitive(xyzzyaaaa4)
xyzzyaaae4=primitive(xyzzyaaaa4+1)-1
do xyzzyaaaf4=xyzzyaaad4,xyzzyaaae4
xyzzyaaaj4=exponent(xyzzyaaaf4)
xyzzyaaai4=xyzzyaaaj4*xyzzyaaag4
xyzzyaaak4=exp(-xyzzyaaai4)
xyzzyaaah4=-xyzzyaaaj4-xyzzyaaaj4
xyzzyaaam4=c_prim(xyzzyaaaf4)*xyzzyaaak4
xyzzyaaal4=xyzzyaaam4*xyzzyaaah4
xyzzyaaao4=r*xyzzyaaal4
xyzzyaaai4=lap_convert+xyzzyaaah4*xyzzyaaag4
xyzzyaaan4=xyzzyaaai4*xyzzyaaal4
ao_m(xyzzyaaac4)=ao_m(xyzzyaaac4)+xyzzyaaam4
alap_m(xyzzyaaac4)=alap_m(xyzzyaaac4)+xyzzyaaan4
agra1_m(xyzzyaaac4)=agra1_m(xyzzyaaac4)+xyzzyaaao4
enddo
enddo
if(.not.virtual)then
xyzzyaaaa1(1)=0.d0
xyzzyaaab1(1)=0.d0
xyzzyaaac1(1)=0.d0
call mxmb(rck(orbital,1:num_ao,1,spin_type_in),num_ao,1,1,ao_m,1,1,xyz&
&zyaaaa1,1,1,1,num_ao,1)
call mxmb(rck(orbital,1:num_ao,1,spin_type_in),num_ao,1,1,alap_m,1,1,x&
&yzzyaaab1,1,1,1,num_ao,1)
call mxmb(rck(orbital,1:num_ao,1,spin_type_in),num_ao,1,1,agra1_m,1,1,&
&xyzzyaaac1,1,1,1,num_ao,1)
rs(1)=xyzzyaaaa1(1)*rnorm
ls(1)=xyzzyaaab1(1)*rnorm
grs(1)=xyzzyaaac1(1)*rnorm
else
rs(1)=rnorm*sum(ao_m(:)*rck_ex(:,orbital,jspin))
ls(1)=rnorm*sum(alap_m(:)*rck_ex(:,orbital,jspin))
grs(1)=rnorm*sum(agra1_m(:)*rck_ex(:,orbital,jspin))
endif
end subroutine xyzzyaaab2
subroutine xyzzyaaac2
implicit none
integer xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5,xyzzyaaad5,xyzzyaaae5,xyzzyaa&
&af5,xyzzyaaag5
real(dp) xyzzyaaah5,xyzzyaaai5,xyzzyaaaj5,xyzzyaaak5,xyzzyaaal5,xyzzya&
&aam5,xyzzyaaan5,xyzzyaaao5,xyzzyaaap5,xyzzyaaaq5,xyzzyaaar5,xyzzyaaas&
&5,xyzzyaaat5
complex(dp) xyzzyaaau5
ao_m=0.d0
xyzzyaaah5=r*r
do xyzzyaaab5=first_shell(ion),first_shell(ion+1)-1
xyzzyaaac5=shell_am(xyzzyaaab5)
if(xyzzyaaac5>2)cycle
xyzzyaaak5=min_exponent(xyzzyaaab5)
xyzzyaaai5=xyzzyaaak5*xyzzyaaah5
if(xyzzyaaai5>screening_tolerance)cycle
xyzzyaaad5=first_ao(xyzzyaaab5)
xyzzyaaae5=primitive(xyzzyaaab5)
xyzzyaaaf5=primitive(xyzzyaaab5+1)-1
do xyzzyaaag5=xyzzyaaae5,xyzzyaaaf5
xyzzyaaak5=exponent(xyzzyaaag5)
xyzzyaaai5=xyzzyaaak5*xyzzyaaah5
xyzzyaaal5=exp(-xyzzyaaai5)
select case (xyzzyaaac5)
case(1)
xyzzyaaaj5=c_prim(xyzzyaaag5)*xyzzyaaal5
case(2)
xyzzyaaaj5=c_prim(xyzzyaaag5)*xyzzyaaal5
case default
xyzzyaaaj5=0.d0
end select
ao_m(xyzzyaaad5)=ao_m(xyzzyaaad5)+xyzzyaaaj5
enddo
enddo
if(.not.virtual)then
xyzzyaaam5=kvec(1,nk)
xyzzyaaan5=kvec(2,nk)
xyzzyaaao5=kvec(3,nk)
xyzzyaaap5=gpcell(1,ion_s)
xyzzyaaaq5=gpcell(2,ion_s)
xyzzyaaar5=gpcell(3,ion_s)
xyzzyaaas5=xyzzyaaam5*xyzzyaaap5+xyzzyaaan5*xyzzyaaaq5+xyzzyaaao5*xyzz&
&yaaar5
if(nk<=num_real_k)then
xyzzyaaaa1(1)=0.d0
xyzzyaaat5=cos(xyzzyaaas5)
rbf(1:num_ao,nk)=ao_m(1:num_ao)*xyzzyaaat5
call mxmb(rck(band,1:num_ao,nk,spin_type_in),num_ao,1,1,rbf(1,nk),1,1,&
&xyzzyaaaa1,1,1,1,num_ao,1)
rs(1)=xyzzyaaaa1(1)*rnorm
else
xyzzyaaad1(1)=czero
xyzzyaaau5=exp(xyzzyaaas5*zi)
do xyzzyaaaa5=1,num_ao
cbf(xyzzyaaaa5,nk)=cmplx(ao_m(xyzzyaaaa5),0.d0,kind=dp)*xyzzyaaau5
enddo
call mxmbc(cck(band,1:num_ao,nk,spin_type_in),num_ao,1,1,cbf(1,nk),1,1&
&,xyzzyaaad1,1,1,1,num_ao,1)
rs(1)=2.d0*dble(xyzzyaaad1(1))*rnorm
rs(2)=2.d0*aimag(xyzzyaaad1(1))*rnorm
endif
else
call errstop('TEMPORARY','No cusp corrected virtual orbitals for perio&
&dic Gaussian orbitals')
endif
end subroutine xyzzyaaac2
subroutine xyzzyaaad2
implicit none
integer xyzzyaaaa6,xyzzyaaab6,xyzzyaaac6,xyzzyaaad6,xyzzyaaae6,xyzzyaa&
&af6
real(dp) xyzzyaaag6,xyzzyaaah6,xyzzyaaai6,xyzzyaaaj6,xyzzyaaak6,xyzzya&
&aal6,xyzzyaaam6,xyzzyaaan6,xyzzyaaao6,xyzzyaaap6,xyzzyaaaq6,xyzzyaaar&
&6,xyzzyaaas6,xyzzyaaat6,xyzzyaaau6,xyzzyaaav6,xyzzyaaaw6
complex(dp) xyzzyaaax6
ao_m=0.d0
alap_m=0.d0
agra1_m=0.d0
xyzzyaaag6=r*r
do xyzzyaaaa6=first_shell(ion),first_shell(ion+1)-1
xyzzyaaab6=shell_am(xyzzyaaaa6)
if(xyzzyaaab6>2)cycle
xyzzyaaaj6=min_exponent(xyzzyaaaa6)
xyzzyaaai6=xyzzyaaaj6*xyzzyaaag6
if(xyzzyaaai6>screening_tolerance)cycle
xyzzyaaac6=first_ao(xyzzyaaaa6)
xyzzyaaad6=primitive(xyzzyaaaa6)
xyzzyaaae6=primitive(xyzzyaaaa6+1)-1
do xyzzyaaaf6=xyzzyaaad6,xyzzyaaae6
xyzzyaaaj6=exponent(xyzzyaaaf6)
xyzzyaaai6=xyzzyaaaj6*xyzzyaaag6
xyzzyaaak6=exp(-xyzzyaaai6)
xyzzyaaah6=-xyzzyaaaj6-xyzzyaaaj6
select case (xyzzyaaab6)
case(1)
xyzzyaaam6=c_prim(xyzzyaaaf6)*xyzzyaaak6
xyzzyaaal6=xyzzyaaam6*xyzzyaaah6
xyzzyaaao6=r*xyzzyaaal6
xyzzyaaai6=lap_convert+xyzzyaaah6*xyzzyaaag6
xyzzyaaan6=xyzzyaaai6*xyzzyaaal6
case(2)
xyzzyaaal6=c_prim(xyzzyaaaf6)*xyzzyaaak6
xyzzyaaam6=xyzzyaaal6
xyzzyaaao6=xyzzyaaah6*r*xyzzyaaal6
xyzzyaaai6=lap_convert+xyzzyaaah6*xyzzyaaag6
xyzzyaaan6=xyzzyaaah6*xyzzyaaai6*xyzzyaaal6
case default
xyzzyaaam6=0.d0
xyzzyaaao6=0.d0
xyzzyaaan6=0.d0
end select
ao_m(xyzzyaaac6)=ao_m(xyzzyaaac6)+xyzzyaaam6
alap_m(xyzzyaaac6)=alap_m(xyzzyaaac6)+xyzzyaaan6
agra1_m(xyzzyaaac6)=agra1_m(xyzzyaaac6)+xyzzyaaao6
enddo
enddo
if(.not.virtual)then
xyzzyaaap6=kvec(1,nk)
xyzzyaaaq6=kvec(2,nk)
xyzzyaaar6=kvec(3,nk)
xyzzyaaas6=gpcell(1,ion_s)
xyzzyaaat6=gpcell(2,ion_s)
xyzzyaaau6=gpcell(3,ion_s)
xyzzyaaav6=xyzzyaaap6*xyzzyaaas6+xyzzyaaaq6*xyzzyaaat6+xyzzyaaar6*xyzz&
&yaaau6
if(nk<=num_real_k)then
xyzzyaaaa1(1)=0.d0
xyzzyaaab1(1)=0.d0
xyzzyaaac1(1)=0.d0
xyzzyaaaw6=cos(xyzzyaaav6)
rbf(1:num_ao,nk)=ao_m(1:num_ao)*xyzzyaaaw6
rblap(1:num_ao,nk)=alap_m(1:num_ao)*xyzzyaaaw6
rbgra1(1:num_ao,nk)=agra1_m(1:num_ao)*xyzzyaaaw6
call mxmb(rck(band,1:num_ao,nk,spin_type_in),num_ao,1,1,rbf(1,nk),1,1,&
&xyzzyaaaa1,1,1,1,num_ao,1)
call mxmb(rck(band,1:num_ao,nk,spin_type_in),num_ao,1,1,rblap(1,nk),1,&
&1,xyzzyaaab1,1,1,1,num_ao,1)
call mxmb(rck(band,1:num_ao,nk,spin_type_in),num_ao,1,1,rbgra1(1,nk),1&
&,1,xyzzyaaac1,1,1,1,num_ao,1)
rs(1)=xyzzyaaaa1(1)*rnorm
ls(1)=xyzzyaaab1(1)*rnorm
grs(1)=xyzzyaaac1(1)*rnorm
else
xyzzyaaad1(1)=czero
xyzzyaaae1(1)=czero
xyzzyaaaf1(1)=czero
xyzzyaaax6=exp(xyzzyaaav6*zi)
cbf(1:num_ao,nk)=cmplx(ao_m(1:num_ao),0.d0,kind=dp)*xyzzyaaax6
cblap(1:num_ao,nk)=cmplx(alap_m(1:num_ao),0.d0,kind=dp)*xyzzyaaax6
cbgra1(1:num_ao,nk)=cmplx(agra1_m(1:num_ao),0.d0,kind=dp)*xyzzyaaax6
call mxmbc(cck(band,1:num_ao,nk,spin_type_in),num_ao,1,1,cbf(1,nk),1,1&
&,xyzzyaaad1,1,1,1,num_ao,1)
call mxmbc(cck(band,1:num_ao,nk,spin_type_in),num_ao,1,1,cblap(1,nk),1&
&,1,xyzzyaaae1,1,1,1,num_ao,1)
call mxmbc(cck(band,1:num_ao,nk,spin_type_in),num_ao,1,1,cbgra1(1,nk),&
&1,1,xyzzyaaaf1,1,1,1,num_ao,1)
xyzzyaaai6=2.d0*rnorm
rs(1)=dble(xyzzyaaad1(1))*xyzzyaaai6
rs(2)=aimag(xyzzyaaad1(1))*xyzzyaaai6
ls(1)=dble(xyzzyaaae1(1))*xyzzyaaai6
ls(2)=aimag(xyzzyaaae1(1))*xyzzyaaai6
grs(1)=dble(xyzzyaaaf1(1))*xyzzyaaai6
grs(2)=aimag(xyzzyaaaf1(1))*xyzzyaaai6
endif
else
call errstop('TEMPORARY','No cusp corrected virtual orbitals for perio&
&dic Gaussian orbitals')
endif
end subroutine xyzzyaaad2
end subroutine wfdet_s
subroutine wfdet_s_bf(r,rsquared,xx,yy,zz,xy,xz,yz,lap_convert,orbital&
&,ion,jspin,band,nk,ion_s,virtual,spin_type_in,rs,grs,sds)
implicit none
integer,intent(in) :: orbital,ion,jspin,band,nk,ion_s,spin_type_in
real(dp),intent(in) :: r,rsquared,xx,yy,zz,xy,xz,yz,lap_convert
real(dp),intent(inout) :: rs(2),grs(2),sds(6,2)
logical,intent(in) :: virtual
if(.not.isperiodic)then
call xyzzyaaaa7
else
call errstop('WFDET_S_BF','Routine not coded for periodic case')
endif
contains
subroutine xyzzyaaaa7
implicit none
integer xyzzyaaaa8,xyzzyaaab8,xyzzyaaac8,xyzzyaaad8,xyzzyaaae8,xyzzyaa&
&af8
real(dp) xyzzyaaag8,xyzzyaaah8,xyzzyaaai8,xyzzyaaaj8,xyzzyaaak8,xyzzya&
&aal8,xyzzyaaam8(1),xyzzyaaan8(1),xyzzyaaao8(1),xyzzyaaap8(1),xyzzyaaa&
&q8(1),xyzzyaaar8(1),xyzzyaaas8(6)
ao_m=0.d0
agra1_m=0.d0
asderiv1_m=0.d0
asderiv2_m=0.d0
asderiv3_m=0.d0
asderiv4_m=0.d0
asderiv5_m=0.d0
asderiv6_m=0.d0
do xyzzyaaac8=first_shell(ion),first_shell(ion+1)-1
xyzzyaaaf8=shell_am(xyzzyaaac8)
if(xyzzyaaaf8>2)cycle
xyzzyaaai8=min_exponent(xyzzyaaac8)
xyzzyaaah8=xyzzyaaai8*rsquared
if(xyzzyaaah8>screening_tolerance)cycle
xyzzyaaaa8=first_ao(xyzzyaaac8)
xyzzyaaad8=primitive(xyzzyaaac8)
xyzzyaaae8=primitive(xyzzyaaac8+1)-1
do xyzzyaaab8=xyzzyaaad8,xyzzyaaae8
xyzzyaaai8=exponent(xyzzyaaab8)
xyzzyaaah8=xyzzyaaai8*rsquared
xyzzyaaak8=c_prim(xyzzyaaab8)*exp(-xyzzyaaah8)
xyzzyaaag8=-xyzzyaaai8-xyzzyaaai8
xyzzyaaaj8=xyzzyaaak8*xyzzyaaag8
xyzzyaaal8=r*xyzzyaaaj8
xyzzyaaas8(1)=(1.d0+xyzzyaaag8*xx)*xyzzyaaaj8
xyzzyaaas8(2)=(1.d0+xyzzyaaag8*yy)*xyzzyaaaj8
xyzzyaaas8(3)=(1.d0+xyzzyaaag8*zz)*xyzzyaaaj8
xyzzyaaaj8=xyzzyaaaj8*xyzzyaaag8
xyzzyaaas8(4)=xy*xyzzyaaaj8
xyzzyaaas8(5)=xz*xyzzyaaaj8
xyzzyaaas8(6)=yz*xyzzyaaaj8
ao_m(xyzzyaaaa8)=ao_m(xyzzyaaaa8)+xyzzyaaak8
agra1_m(xyzzyaaaa8)=agra1_m(xyzzyaaaa8)+xyzzyaaal8
asderiv1_m(xyzzyaaaa8)=asderiv1_m(xyzzyaaaa8)+xyzzyaaas8(1)
asderiv2_m(xyzzyaaaa8)=asderiv2_m(xyzzyaaaa8)+xyzzyaaas8(2)
asderiv3_m(xyzzyaaaa8)=asderiv3_m(xyzzyaaaa8)+xyzzyaaas8(3)
asderiv4_m(xyzzyaaaa8)=asderiv4_m(xyzzyaaaa8)+xyzzyaaas8(4)
asderiv5_m(xyzzyaaaa8)=asderiv5_m(xyzzyaaaa8)+xyzzyaaas8(5)
asderiv6_m(xyzzyaaaa8)=asderiv6_m(xyzzyaaaa8)+xyzzyaaas8(6)
enddo
enddo
if(.not.virtual)then
xyzzyaaaa1(1)=0.d0
xyzzyaaab1(1)=0.d0
xyzzyaaac1(1)=0.d0
xyzzyaaam8(1)=0.d0
xyzzyaaan8(1)=0.d0
xyzzyaaao8(1)=0.d0
xyzzyaaap8(1)=0.d0
xyzzyaaaq8(1)=0.d0
xyzzyaaar8(1)=0.d0
call mxmb(rck(orbital,1:num_ao,1,spin_type_in),num_ao,1,1,ao_m,1,1,xyz&
&zyaaaa1,1,1,1,num_ao,1)
call mxmb(rck(orbital,1:num_ao,1,spin_type_in),num_ao,1,1,agra1_m,1,1,&
&xyzzyaaac1,1,1,1,num_ao,1)
call mxmb(rck(orbital,1:num_ao,1,spin_type_in),num_ao,1,1,asderiv1_m,1&
&,1,xyzzyaaam8,1,1,1,num_ao,1)
call mxmb(rck(orbital,1:num_ao,1,spin_type_in),num_ao,1,1,asderiv2_m,1&
&,1,xyzzyaaan8,1,1,1,num_ao,1)
call mxmb(rck(orbital,1:num_ao,1,spin_type_in),num_ao,1,1,asderiv3_m,1&
&,1,xyzzyaaao8,1,1,1,num_ao,1)
call mxmb(rck(orbital,1:num_ao,1,spin_type_in),num_ao,1,1,asderiv4_m,1&
&,1,xyzzyaaap8,1,1,1,num_ao,1)
call mxmb(rck(orbital,1:num_ao,1,spin_type_in),num_ao,1,1,asderiv5_m,1&
&,1,xyzzyaaaq8,1,1,1,num_ao,1)
call mxmb(rck(orbital,1:num_ao,1,spin_type_in),num_ao,1,1,asderiv6_m,1&
&,1,xyzzyaaar8,1,1,1,num_ao,1)
rs(1)=xyzzyaaaa1(1)*rnorm
grs(1)=xyzzyaaac1(1)*rnorm
sds(1,1)=xyzzyaaam8(1)*rnorm
sds(2,1)=xyzzyaaan8(1)*rnorm
sds(3,1)=xyzzyaaao8(1)*rnorm
sds(4,1)=xyzzyaaap8(1)*rnorm
sds(5,1)=xyzzyaaaq8(1)*rnorm
sds(6,1)=xyzzyaaar8(1)*rnorm
else
rs(1)=rnorm*sum(ao_m(:)*rck_ex(:,orbital,jspin))
grs(1)=rnorm*sum(agra1_m(:)*rck_ex(:,orbital,jspin))
sds(1,1)=rnorm*sum(asderiv1_m(:)*rck_ex(:,orbital,jspin))
sds(2,1)=rnorm*sum(asderiv2_m(:)*rck_ex(:,orbital,jspin))
sds(3,1)=rnorm*sum(asderiv3_m(:)*rck_ex(:,orbital,jspin))
sds(4,1)=rnorm*sum(asderiv4_m(:)*rck_ex(:,orbital,jspin))
sds(5,1)=rnorm*sum(asderiv5_m(:)*rck_ex(:,orbital,jspin))
sds(6,1)=rnorm*sum(asderiv6_m(:)*rck_ex(:,orbital,jspin))
endif
end subroutine xyzzyaaaa7
end subroutine wfdet_s_bf
end module slaarnaai

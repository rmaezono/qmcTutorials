module slaarnaba
use dsp
use slaarnaat
use slaarnabg
use slaarnabp
use store
use slaarnaag, only : czero,zi
use slaarnabi,      only : use_gpcc,cusp_wfdet_cmplx
use slaarnabt, only : gmult,g3mult,dcopy,multi_ddot_s,multi_zdotu_s
implicit none
private
public gauss_per_orb_eval
integer xyzzyaaaa1,xyzzyaaab1
contains
subroutine gauss_per_orb_eval(rvec,jspin,norb,orbmask,norbc,corbmask,v&
&al,fsd,orbval,orbgrad,orblap)
implicit none
integer,intent(in) :: jspin,norb,norbc
real(dp),intent(in) :: rvec(3)
real(dp),intent(inout) :: orbval(norb,real1_complex2),orbgrad(3,norb,r&
&eal1_complex2),orblap(norb,real1_complex2)
logical,intent(in) :: val,fsd,orbmask(norb),corbmask(norbc)
integer xyzzyaaaa2
real(dp) xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2,xyzzyaaaf2,xyzzya&
&aag2,xyzzyaaah2,xyzzyaaai2(3)
select case (periodicity)
case(3)
xyzzyaaai2=anint(matmul(rvec,painv))
new_origin=matmul(xyzzyaaai2,pamat)
case(2)
xyzzyaaai2(1:2)=anint(matmul(rvec(1:2),painv(1:2,1:2)))
new_origin(1:2)=matmul(xyzzyaaai2(1:2),pamat(1:2,1:2))
new_origin(3)=0.d0
case(1)
xyzzyaaai2(1)=anint(rvec(1)*painv(1,1))
new_origin(1)=xyzzyaaai2(1)*pamat(1,1)
new_origin(2:3)=0.d0
end select
new_origin_to_rvec=rvec-new_origin
xyzzyaaab2=new_origin(1)
xyzzyaaac2=new_origin(2)
xyzzyaaad2=new_origin(3)
do xyzzyaaaa2=1,num_real_k
xyzzyaaae2=kvec(1,xyzzyaaaa2)
xyzzyaaaf2=kvec(2,xyzzyaaaa2)
xyzzyaaag2=kvec(3,xyzzyaaaa2)
xyzzyaaah2=xyzzyaaae2*xyzzyaaab2+xyzzyaaaf2*xyzzyaaac2+xyzzyaaag2*xyzz&
&yaaad2
coskdotg_shift(xyzzyaaaa2)=cos(xyzzyaaah2)
enddo
do xyzzyaaaa2=num_real_k_plus_1,num_k
xyzzyaaae2=kvec(1,xyzzyaaaa2)
xyzzyaaaf2=kvec(2,xyzzyaaaa2)
xyzzyaaag2=kvec(3,xyzzyaaaa2)
xyzzyaaah2=xyzzyaaae2*xyzzyaaab2+xyzzyaaaf2*xyzzyaaac2+xyzzyaaag2*xyzz&
&yaaad2
expikdotg_shift(xyzzyaaaa2)=exp(xyzzyaaah2*zi)
enddo
xyzzyaaaa1=which_ssingle(jspin,spin_dep_gs)
xyzzyaaab1=which_ssingle(jspin,spin_dep_in)
if(.not.fsd)then
call xyzzyaaac1(rvec,jspin,norb,orbmask,norbc,corbmask,orbval)
elseif(.not.val)then
call xyzzyaaad1(rvec,jspin,norb,orbmask,norbc,corbmask,orbgrad,orblap)
else
call xyzzyaaae1(rvec,jspin,norb,orbmask,norbc,corbmask,orbval,orbgrad,&
&orblap)
endif
end subroutine gauss_per_orb_eval
subroutine xyzzyaaac1(rvec,jspin,norb,orbmask,norbc,corbmask,orbval)
implicit none
integer,intent(in) :: jspin,norb,norbc
real(dp),intent(in) :: rvec(3)
real(dp),intent(inout) :: orbval(norb,real1_complex2)
logical,intent(in) :: orbmask(norb),corbmask(norbc)
integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3,xyzzyaaad3,xyzzyaaae3,xyzzyaa&
&af3,xyzzyaaag3,xyzzyaaah3,xyzzyaaai3,xyzzyaaaj3,xyzzyaaak3,xyzzyaaal3&
&,xyzzyaaam3,xyzzyaaan3,xyzzyaaao3,xyzzyaaap3,xyzzyaaaq3,xyzzyaaar3,xy&
&zzyaaas3(num_ao),xyzzyaaat3,xyzzyaaau3
real(dp) xyzzyaaav3,xyzzyaaaw3,xyzzyaaax3,xyzzyaaay3,xyzzyaaaz3,xyzzya&
&aba3,xyzzyaabb3,xyzzyaabc3,xyzzyaabd3,xyzzyaabe3,xyzzyaabf3(7),xyzzya&
&abg3,xyzzyaabh3,xyzzyaabi3,xyzzyaabj3,xyzzyaabk3,xyzzyaabl3,xyzzyaabm&
&3,xyzzyaabn3,xyzzyaabo3,xyzzyaabp3,xyzzyaabq3,xyzzyaabr3
logical xyzzyaabs3,xyzzyaabt3(num_shells)
xyzzyaaaa3=0
rbf=0.d0
xyzzyaabt3=.false.
if(complex_states)cbf=czero
do xyzzyaaab3=1,num_cell
do xyzzyaaac3=1,num_real_k
phase(xyzzyaaac3)=coskdotg(xyzzyaaac3,xyzzyaaab3)*coskdotg_shift(xyzzy&
&aaac3)
enddo
do xyzzyaaac3=num_real_k_plus_1,num_k
cphase(xyzzyaaac3)=expikdotg(xyzzyaaac3,xyzzyaaab3)*expikdotg_shift(xy&
&zzyaaac3)
enddo
do xyzzyaaad3=1,num_centres_in_cell(xyzzyaaab3)
xyzzyaaav3=xpos_in_cell(xyzzyaaad3,xyzzyaaab3)
xyzzyaaaw3=ypos_in_cell(xyzzyaaad3,xyzzyaaab3)
xyzzyaaax3=zpos_in_cell(xyzzyaaad3,xyzzyaaab3)
xyzzyaaav3=new_origin_to_rvec(1)-xyzzyaaav3
xyzzyaaaw3=new_origin_to_rvec(2)-xyzzyaaaw3
xyzzyaaax3=new_origin_to_rvec(3)-xyzzyaaax3
xyzzyaaay3=xyzzyaaav3*xyzzyaaav3
xyzzyaaaz3=xyzzyaaaw3*xyzzyaaaw3
xyzzyaaba3=xyzzyaaax3*xyzzyaaax3
xyzzyaabb3=xyzzyaaay3+xyzzyaaaz3+xyzzyaaba3
xyzzyaaaf3=num_shells_on_centre(xyzzyaaad3,xyzzyaaab3)
do xyzzyaaag3=1,xyzzyaaaf3
xyzzyaaaa3=xyzzyaaaa3+1
xyzzyaaah3=shell_sequence_number(xyzzyaaaa3)
xyzzyaabd3=min_exponent(xyzzyaaah3)
xyzzyaabc3=xyzzyaabd3*xyzzyaabb3
if(xyzzyaabc3>screening_tolerance)cycle
xyzzyaaai3=shell_am(xyzzyaaah3)
xyzzyaaaj3=numao_in_shell(xyzzyaaah3)
do xyzzyaaaq3=1,xyzzyaaaj3
xyzzyaabf3(xyzzyaaaq3)=0.d0
enddo
xyzzyaaal3=primitive(xyzzyaaah3)
xyzzyaaam3=primitive(xyzzyaaah3+1)-1
xyzzyaabs3=.false.
do xyzzyaaan3=xyzzyaaal3,xyzzyaaam3
xyzzyaabd3=exponent(xyzzyaaan3)
xyzzyaabc3=xyzzyaabd3*xyzzyaabb3
if(xyzzyaabc3>30.d0)cycle
xyzzyaabs3=.true.
xyzzyaabt3(xyzzyaaah3)=.true.
xyzzyaabe3=exp(-xyzzyaabc3)
xyzzyaabc3=c_prim(xyzzyaaan3)*xyzzyaabe3
select case (xyzzyaaai3)
case(1)
xyzzyaabf3(1)=xyzzyaabf3(1)+xyzzyaabc3
case(2)
xyzzyaabf3(1)=xyzzyaabf3(1)+xyzzyaabc3
xyzzyaabc3=c_prim2(xyzzyaaan3)*xyzzyaabe3
xyzzyaabf3(2)=xyzzyaabf3(2)+xyzzyaabc3*xyzzyaaav3
xyzzyaabf3(3)=xyzzyaabf3(3)+xyzzyaabc3*xyzzyaaaw3
xyzzyaabf3(4)=xyzzyaabf3(4)+xyzzyaabc3*xyzzyaaax3
case(3)
xyzzyaabf3(1)=xyzzyaabf3(1)+xyzzyaabc3*xyzzyaaav3
xyzzyaabf3(2)=xyzzyaabf3(2)+xyzzyaabc3*xyzzyaaaw3
xyzzyaabf3(3)=xyzzyaabf3(3)+xyzzyaabc3*xyzzyaaax3
case(4)
xyzzyaabf3(1)=xyzzyaabf3(1)+xyzzyaabc3*(3.d0*xyzzyaaba3-xyzzyaabb3)
xyzzyaabf3(2)=xyzzyaabf3(2)+xyzzyaabc3*xyzzyaaav3*xyzzyaaax3
xyzzyaabf3(3)=xyzzyaabf3(3)+xyzzyaabc3*xyzzyaaaw3*xyzzyaaax3
xyzzyaabf3(4)=xyzzyaabf3(4)+xyzzyaabc3*(xyzzyaaay3-xyzzyaaaz3)
xyzzyaabf3(5)=xyzzyaabf3(5)+xyzzyaabc3*xyzzyaaav3*xyzzyaaaw3
case(5)
xyzzyaabg3=xyzzyaaay3*xyzzyaaav3
xyzzyaabh3=xyzzyaaaz3*xyzzyaaaw3
xyzzyaabi3=xyzzyaaba3*xyzzyaaax3
xyzzyaabj3=xyzzyaaav3*xyzzyaaaw3
xyzzyaabk3=xyzzyaaav3*xyzzyaaax3
xyzzyaabl3=xyzzyaaaw3*xyzzyaaax3
xyzzyaabm3=xyzzyaabj3*xyzzyaaav3
xyzzyaabn3=xyzzyaabj3*xyzzyaaaw3
xyzzyaabq3=xyzzyaabk3*xyzzyaaav3
xyzzyaabr3=xyzzyaabk3*xyzzyaaax3
xyzzyaabo3=xyzzyaabl3*xyzzyaaaw3
xyzzyaabp3=xyzzyaabl3*xyzzyaaax3
xyzzyaabf3(1)=xyzzyaabf3(1)+xyzzyaabc3*(xyzzyaabi3-1.5d0*(xyzzyaabq3+x&
&yzzyaabo3))
xyzzyaabf3(2)=xyzzyaabf3(2)+xyzzyaabc3*(6.d0*xyzzyaabr3-1.5d0*(xyzzyaa&
&bg3+xyzzyaabn3))
xyzzyaabf3(3)=xyzzyaabf3(3)+xyzzyaabc3*(6.d0*xyzzyaabp3-1.5d0*(xyzzyaa&
&bm3+xyzzyaabh3))
xyzzyaabf3(4)=xyzzyaabf3(4)+xyzzyaabc3*15.d0*(xyzzyaabq3-xyzzyaabo3)
xyzzyaabf3(5)=xyzzyaabf3(5)+xyzzyaabc3*30.d0*xyzzyaabj3*xyzzyaaax3
xyzzyaabf3(6)=xyzzyaabf3(6)+xyzzyaabc3*(15.d0*xyzzyaabg3-45.d0*xyzzyaa&
&bn3)
xyzzyaabf3(7)=xyzzyaabf3(7)+xyzzyaabc3*(45.d0*xyzzyaabm3-15.d0*xyzzyaa&
&bh3)
end select
enddo
if(xyzzyaabs3)then
xyzzyaaak3=first_ao(xyzzyaaah3)-1
if(num_real_k_gt_1)then
do xyzzyaaac3=1,num_real_k-1,2
if(xyzzyaaaj3==4)then
do xyzzyaaaq3=1,4
rbf(xyzzyaaak3+xyzzyaaaq3,xyzzyaaac3)=rbf(xyzzyaaak3+xyzzyaaaq3,xyzzya&
&aac3)+xyzzyaabf3(xyzzyaaaq3)*phase(xyzzyaaac3)
rbf(xyzzyaaak3+xyzzyaaaq3,xyzzyaaac3+1)=rbf(xyzzyaaak3+xyzzyaaaq3,xyzz&
&yaaac3+1)+xyzzyaabf3(xyzzyaaaq3)*phase(xyzzyaaac3+1)
enddo
else
do xyzzyaaaq3=1,xyzzyaaaj3
rbf(xyzzyaaak3+xyzzyaaaq3,xyzzyaaac3)=rbf(xyzzyaaak3+xyzzyaaaq3,xyzzya&
&aac3)+xyzzyaabf3(xyzzyaaaq3)*phase(xyzzyaaac3)
rbf(xyzzyaaak3+xyzzyaaaq3,xyzzyaaac3+1)=rbf(xyzzyaaak3+xyzzyaaaq3,xyzz&
&yaaac3+1)+xyzzyaabf3(xyzzyaaaq3)*phase(xyzzyaaac3+1)
enddo
endif
enddo
endif
if(num_real_k_odd)then
do xyzzyaaaq3=1,xyzzyaaaj3
rbf(xyzzyaaak3+xyzzyaaaq3,num_real_k)=rbf(xyzzyaaak3+xyzzyaaaq3,num_re&
&al_k)+xyzzyaabf3(xyzzyaaaq3)*phase(num_real_k)
enddo
endif
do xyzzyaaac3=num_real_k_plus_1,num_k
do xyzzyaaaq3=1,xyzzyaaaj3
cbf(xyzzyaaak3+xyzzyaaaq3,xyzzyaaac3)=cbf(xyzzyaaak3+xyzzyaaaq3,xyzzya&
&aac3)+cmplx(xyzzyaabf3(xyzzyaaaq3),0.d0,kind=dp)*cphase(xyzzyaaac3)
enddo
enddo
endif
enddo
enddo
enddo
xyzzyaaar3=0
do xyzzyaaah3=1,num_shells
if(.not.xyzzyaabt3(xyzzyaaah3))cycle
xyzzyaaak3=first_ao(xyzzyaaah3)-1
do xyzzyaaaq3=1,numao_in_shell(xyzzyaaah3)
xyzzyaaar3=xyzzyaaar3+1
xyzzyaaas3(xyzzyaaar3)=xyzzyaaak3+xyzzyaaaq3
enddo
enddo
if(use_gpcc)call cusp_wfdet_cmplx(rvec,jspin,.true.,.false.,cusp_val_c&
&,cusp_grad_c,cusp_lap_c)
xyzzyaaao3=0
xyzzyaaap3=0
xyzzyaaat3=gauss_offset(jspin)
xyzzyaaau3=gauss_offsetc(jspin)
do xyzzyaaac3=1,num_real_k
xyzzyaaaf3=nband(xyzzyaaac3,xyzzyaaaa1)
orb1(xyzzyaaao3+1:xyzzyaaao3+xyzzyaaaf3)=multi_ddot_s(xyzzyaaaf3,xyzzy&
&aaar3,xyzzyaaas3,corbmask(xyzzyaaap3+xyzzyaaau3),rck(:,:,xyzzyaaac3,x&
&yzzyaaab1),maxb*num_ao,maxb,rbf(1,xyzzyaaac3),rnorm)
xyzzyaaao3=xyzzyaaao3+xyzzyaaaf3
xyzzyaaap3=xyzzyaaap3+xyzzyaaaf3
enddo
if(use_gpcc)then
do xyzzyaaaq3=1,xyzzyaaao3
if(.not.corbmask(xyzzyaaaq3+xyzzyaaau3-1))cycle
orb1(xyzzyaaaq3)=orb1(xyzzyaaaq3)+dble(cusp_val_c(xyzzyaaaq3+xyzzyaaau&
&3-1))
enddo
endif
if(complex_states)then
do xyzzyaaac3=num_real_k_plus_1,num_k
xyzzyaaaf3=nband(xyzzyaaac3,xyzzyaaaa1)
if(complex_wf)then
corb1(1:xyzzyaaaf3)=multi_zdotu_s(xyzzyaaaf3,xyzzyaaar3,xyzzyaaas3,cor&
&bmask(xyzzyaaap3+xyzzyaaau3),cck(:,:,xyzzyaaac3,xyzzyaaab1),maxb*num_&
&ao,maxb,cbf(1,xyzzyaaac3),rnorm)
else
corb1(1:xyzzyaaaf3)=multi_zdotu_s(xyzzyaaaf3,xyzzyaaar3,xyzzyaaas3,cor&
&bmask(xyzzyaaap3+xyzzyaaau3),cck(:,:,xyzzyaaac3,xyzzyaaab1),maxb*num_&
&ao,maxb,cbf(1,xyzzyaaac3),two_rnorm)
endif
do xyzzyaaaq3=1,xyzzyaaaf3
xyzzyaaao3=xyzzyaaao3+1
xyzzyaaap3=xyzzyaaap3+1
if(.not.corbmask(xyzzyaaap3+xyzzyaaau3-1))cycle
if(use_gpcc)then
corb1(xyzzyaaaq3)=corb1(xyzzyaaaq3)+cusp_val_c(xyzzyaaap3+xyzzyaaau3-1&
&)
endif
orb1(xyzzyaaao3)=dble(corb1(xyzzyaaaq3))
if(xyzzyaaao3+1<=nuc_nele(jspin))then
xyzzyaaao3=xyzzyaaao3+1
orb1(xyzzyaaao3)=aimag(corb1(xyzzyaaaq3))
endif
enddo
enddo
endif
call dcopy(xyzzyaaao3,orb1(1),1,orbval(xyzzyaaat3,1),1)
if(excite)then
do xyzzyaaaq3=1,gauss_ex_norb
xyzzyaaae3=gauss_gs_norb+xyzzyaaaq3
if(.not.orbmask(xyzzyaaae3))cycle
xyzzyaaac3=virtual_k(xyzzyaaaq3,jspin)
if(xyzzyaaac3<=num_real_k)then
psi(xyzzyaaaq3)=sum(rbf(:,xyzzyaaac3)*rck_ex(:,full2vrt(xyzzyaaaq3,jsp&
&in),jspin))*rnorm
else
psi(xyzzyaaaq3)=sum(dble(cbf(:,xyzzyaaac3)*cck_ex(:,full2vrt(xyzzyaaaq&
&3,jspin),jspin)))*two_rnorm
endif
if(use_gpcc)then
psi(xyzzyaaaq3)=psi(xyzzyaaaq3)+dble(cusp_val_c(ridx2zidx(xyzzyaaae3))&
&)
endif
orbval(xyzzyaaae3,1)=psi(xyzzyaaaq3)
enddo
endif
end subroutine xyzzyaaac1
subroutine xyzzyaaad1(rvec,jspin,norb,orbmask,norbc,corbmask,orbgrad,o&
&rblap)
implicit none
integer,intent(in) :: jspin,norb,norbc
real(dp),intent(in) :: rvec(3)
real(dp),intent(inout) :: orbgrad(3,norb,real1_complex2),orblap(norb,r&
&eal1_complex2)
logical,intent(in) :: orbmask(norb),corbmask(norbc)
integer xyzzyaaaa4,xyzzyaaab4,xyzzyaaac4,xyzzyaaad4,xyzzyaaae4,xyzzyaa&
&af4,xyzzyaaag4,xyzzyaaah4,xyzzyaaai4,xyzzyaaaj4,xyzzyaaak4,xyzzyaaal4&
&,xyzzyaaam4,xyzzyaaan4,xyzzyaaao4,xyzzyaaap4,xyzzyaaaq4,xyzzyaaar4,xy&
&zzyaaas4,xyzzyaaat4,xyzzyaaau4,xyzzyaaav4,xyzzyaaaw4(num_ao)
real(dp) xyzzyaaax4,xyzzyaaay4,xyzzyaaaz4,xyzzyaaba4,xyzzyaabb4,xyzzya&
&abc4,xyzzyaabd4,xyzzyaabe4,xyzzyaabf4,xyzzyaabg4,xyzzyaabh4,xyzzyaabi&
&4,xyzzyaabj4,xyzzyaabk4,xyzzyaabl4,xyzzyaabm4,xyzzyaabn4,xyzzyaabo4,x&
&yzzyaabp4,xyzzyaabq4,xyzzyaabr4,xyzzyaabs4,xyzzyaabt4,xyzzyaabu4,xyzz&
&yaabv4,xyzzyaabw4,xyzzyaabx4,xyzzyaaby4,xyzzyaabz4,xyzzyaaca4,xyzzyaa&
&cb4,xyzzyaacc4,xyzzyaacd4,xyzzyaace4,xyzzyaacf4,xyzzyaacg4,xyzzyaach4&
&,xyzzyaaci4,xyzzyaacj4,xyzzyaack4,xyzzyaacl4(7),xyzzyaacm4(21),xyzzya&
&acn4(7),xyzzyaaco4(21),xyzzyaacp4,xyzzyaacq4,xyzzyaacr4,xyzzyaacs4,xy&
&zzyaact4,xyzzyaacu4,xyzzyaacv4,xyzzyaacw4,xyzzyaacx4
complex(dp) xyzzyaacy4
logical xyzzyaacz4,xyzzyaada4(num_shells)
xyzzyaaas4=0
rblap=0.d0
rbgra1=0.d0
rbgra2=0.d0
rbgra3=0.d0
xyzzyaada4=.false.
if(complex_states)then
cblap=czero
cbgra1=czero
cbgra2=czero
cbgra3=czero
endif
do xyzzyaaaa4=1,num_cell
do xyzzyaaab4=1,num_real_k
phase(xyzzyaaab4)=coskdotg(xyzzyaaab4,xyzzyaaaa4)*coskdotg_shift(xyzzy&
&aaab4)
enddo
do xyzzyaaab4=num_real_k_plus_1,num_k
cphase(xyzzyaaab4)=expikdotg(xyzzyaaab4,xyzzyaaaa4)*expikdotg_shift(xy&
&zzyaaab4)
enddo
do xyzzyaaac4=1,num_centres_in_cell(xyzzyaaaa4)
xyzzyaaay4=xpos_in_cell(xyzzyaaac4,xyzzyaaaa4)
xyzzyaaaz4=ypos_in_cell(xyzzyaaac4,xyzzyaaaa4)
xyzzyaaba4=zpos_in_cell(xyzzyaaac4,xyzzyaaaa4)
xyzzyaaay4=new_origin_to_rvec(1)-xyzzyaaay4
xyzzyaaaz4=new_origin_to_rvec(2)-xyzzyaaaz4
xyzzyaaba4=new_origin_to_rvec(3)-xyzzyaaba4
xyzzyaabb4=xyzzyaaay4*xyzzyaaay4
xyzzyaabc4=xyzzyaaaz4*xyzzyaaaz4
xyzzyaabd4=xyzzyaaba4*xyzzyaaba4
xyzzyaaax4=xyzzyaabb4+xyzzyaabc4+xyzzyaabd4
xyzzyaaad4=num_shells_on_centre(xyzzyaaac4,xyzzyaaaa4)
do xyzzyaaae4=1,xyzzyaaad4
xyzzyaaas4=xyzzyaaas4+1
xyzzyaaaf4=shell_sequence_number(xyzzyaaas4)
xyzzyaacg4=min_exponent(xyzzyaaaf4)
xyzzyaabp4=xyzzyaacg4*xyzzyaaax4
if(xyzzyaabp4>screening_tolerance)cycle
xyzzyaaag4=shell_am(xyzzyaaaf4)
xyzzyaaah4=numao_in_shell(xyzzyaaaf4)
xyzzyaaai4=xyzzyaaah4*3
do xyzzyaaan4=1,xyzzyaaah4
xyzzyaacn4(xyzzyaaan4)=0.d0
enddo
do xyzzyaaan4=1,xyzzyaaai4
xyzzyaaco4(xyzzyaaan4)=0.d0
enddo
xyzzyaaaq4=primitive(xyzzyaaaf4)
xyzzyaaar4=primitive(xyzzyaaaf4+1)-1
xyzzyaacz4=.false.
do xyzzyaaap4=xyzzyaaaq4,xyzzyaaar4
xyzzyaacg4=exponent(xyzzyaaap4)
xyzzyaabp4=xyzzyaacg4*xyzzyaaax4
if(xyzzyaabp4>30.d0)cycle
xyzzyaacz4=.true.
xyzzyaada4(xyzzyaaaf4)=.true.
xyzzyaach4=exp(-xyzzyaabp4)
xyzzyaabl4=-xyzzyaacg4-xyzzyaacg4
select case (xyzzyaaag4)
case(1)
xyzzyaaci4=c_prim(xyzzyaaap4)*xyzzyaach4*xyzzyaabl4
xyzzyaacm4(1)=xyzzyaaay4*xyzzyaaci4
xyzzyaacm4(2)=xyzzyaaaz4*xyzzyaaci4
xyzzyaacm4(3)=xyzzyaaba4*xyzzyaaci4
xyzzyaabp4=3.d0+xyzzyaabl4*xyzzyaaax4
xyzzyaacl4(1)=xyzzyaabp4*xyzzyaaci4
case(2)
xyzzyaacm4(1)=xyzzyaabl4*xyzzyaaay4
xyzzyaacm4(2)=xyzzyaabl4*xyzzyaaaz4
xyzzyaacm4(3)=xyzzyaabl4*xyzzyaaba4
xyzzyaacm4(4)=1.d0+xyzzyaacm4(1)*xyzzyaaay4
xyzzyaacm4(5)=xyzzyaacm4(1)*xyzzyaaaz4
xyzzyaacm4(6)=xyzzyaacm4(1)*xyzzyaaba4
xyzzyaacm4(8)=1.d0+xyzzyaacm4(2)*xyzzyaaaz4
xyzzyaacm4(9)=xyzzyaacm4(2)*xyzzyaaba4
xyzzyaacm4(7)=xyzzyaacm4(5)
xyzzyaacm4(10)=xyzzyaacm4(6)
xyzzyaacm4(11)=xyzzyaacm4(9)
xyzzyaacm4(12)=1.d0+xyzzyaacm4(3)*xyzzyaaba4
xyzzyaaci4=c_prim(xyzzyaaap4)*xyzzyaach4
xyzzyaacj4=c_prim2(xyzzyaaap4)*xyzzyaach4
xyzzyaacm4(1:3)=xyzzyaacm4(1:3)*xyzzyaaci4
xyzzyaacm4(4:12)=xyzzyaacm4(4:12)*xyzzyaacj4
xyzzyaabp4=3.d0+xyzzyaabl4*xyzzyaaax4
xyzzyaabq4=xyzzyaabl4*xyzzyaabp4
xyzzyaabr4=(xyzzyaabq4+2.d0*xyzzyaabl4)*xyzzyaacj4
xyzzyaacl4(1)=xyzzyaabq4*xyzzyaaci4
xyzzyaacl4(2)=xyzzyaabr4*xyzzyaaay4
xyzzyaacl4(3)=xyzzyaabr4*xyzzyaaaz4
xyzzyaacl4(4)=xyzzyaabr4*xyzzyaaba4
case(3)
xyzzyaabm4=xyzzyaabl4*xyzzyaaay4
xyzzyaabn4=xyzzyaabl4*xyzzyaaaz4
xyzzyaabo4=xyzzyaabl4*xyzzyaaba4
xyzzyaabs4=xyzzyaabl4*(5+xyzzyaabl4*xyzzyaaax4)
xyzzyaacm4(1)=1.d0+xyzzyaabm4*xyzzyaaay4
xyzzyaacm4(2)=xyzzyaabm4*xyzzyaaaz4
xyzzyaacm4(3)=xyzzyaabm4*xyzzyaaba4
xyzzyaacm4(5)=1.d0+xyzzyaabn4*xyzzyaaaz4
xyzzyaacm4(6)=xyzzyaabn4*xyzzyaaba4
xyzzyaacm4(9)=1.d0+xyzzyaabo4*xyzzyaaba4
xyzzyaacm4(4)=xyzzyaacm4(2)
xyzzyaacm4(7)=xyzzyaacm4(3)
xyzzyaacm4(8)=xyzzyaacm4(6)
xyzzyaacl4(1)=xyzzyaabs4*xyzzyaaay4
xyzzyaacl4(2)=xyzzyaabs4*xyzzyaaaz4
xyzzyaacl4(3)=xyzzyaabs4*xyzzyaaba4
xyzzyaacj4=c_prim(xyzzyaaap4)*xyzzyaach4
xyzzyaacm4(1:9)=xyzzyaacm4(1:9)*xyzzyaacj4
xyzzyaacl4(1:3)=xyzzyaacl4(1:3)*xyzzyaacj4
case(4)
xyzzyaabe4=xyzzyaaay4*xyzzyaaaz4
xyzzyaabf4=xyzzyaaaz4*xyzzyaaba4
xyzzyaabg4=xyzzyaaay4*xyzzyaaba4
xyzzyaabh4=xyzzyaabb4-xyzzyaabc4
xyzzyaabi4=3.d0*xyzzyaabd4-xyzzyaaax4
xyzzyaabm4=xyzzyaabl4*xyzzyaaay4
xyzzyaabn4=xyzzyaabl4*xyzzyaaaz4
xyzzyaabo4=xyzzyaabl4*xyzzyaaba4
xyzzyaabp4=1.d0+xyzzyaabm4*xyzzyaaay4
xyzzyaabq4=1.d0+xyzzyaabn4*xyzzyaaaz4
xyzzyaabr4=1.d0+xyzzyaabo4*xyzzyaaba4
xyzzyaabs4=xyzzyaabl4*xyzzyaabi4-2.d0
xyzzyaabt4=xyzzyaabl4*xyzzyaabh4
xyzzyaacm4(1)=xyzzyaabs4*xyzzyaaay4
xyzzyaacm4(2)=xyzzyaabs4*xyzzyaaaz4
xyzzyaacm4(3)=(xyzzyaabs4+6.d0)*xyzzyaaba4
xyzzyaacm4(4)=xyzzyaabp4*xyzzyaaba4
xyzzyaacm4(5)=xyzzyaabm4*xyzzyaabf4
xyzzyaacm4(6)=xyzzyaabr4*xyzzyaaay4
xyzzyaacm4(7)=xyzzyaacm4(5)
xyzzyaacm4(8)=xyzzyaabq4*xyzzyaaba4
xyzzyaacm4(9)=xyzzyaabr4*xyzzyaaaz4
xyzzyaacm4(10)=(2.d0+xyzzyaabl4*xyzzyaabh4)*xyzzyaaay4
xyzzyaacm4(11)=(xyzzyaabt4-2.d0)*xyzzyaaaz4
xyzzyaacm4(12)=xyzzyaabt4*xyzzyaaba4
xyzzyaacm4(13)=xyzzyaabp4*xyzzyaaaz4
xyzzyaacm4(14)=xyzzyaabq4*xyzzyaaay4
xyzzyaacm4(15)=xyzzyaacm4(5)
xyzzyaack4=c_prim(xyzzyaaap4)*xyzzyaach4
xyzzyaacm4(1:15)=xyzzyaacm4(1:15)*xyzzyaack4
xyzzyaabp4=7.d0+xyzzyaabl4*xyzzyaaax4
xyzzyaabp4=xyzzyaabl4*xyzzyaabp4*xyzzyaack4
xyzzyaacl4(1)=xyzzyaabp4*xyzzyaabi4
xyzzyaacl4(2)=xyzzyaabp4*xyzzyaabg4
xyzzyaacl4(3)=xyzzyaabp4*xyzzyaabf4
xyzzyaacl4(4)=xyzzyaabp4*xyzzyaabh4
xyzzyaacl4(5)=xyzzyaabp4*xyzzyaabe4
case(5)
xyzzyaabe4=xyzzyaaay4*xyzzyaaaz4
xyzzyaabg4=xyzzyaaay4*xyzzyaaba4
xyzzyaabf4=xyzzyaaaz4*xyzzyaaba4
xyzzyaacx4=xyzzyaabe4*xyzzyaaba4
xyzzyaabj4=xyzzyaabb4+xyzzyaabc4
xyzzyaabk4=xyzzyaabc4-xyzzyaabb4
xyzzyaacp4=2.d0*xyzzyaacg4
xyzzyaacq4=3.d0*xyzzyaacg4
xyzzyaacr4=xyzzyaacp4*xyzzyaabb4
xyzzyaacs4=xyzzyaacp4*xyzzyaabc4
xyzzyaact4=xyzzyaacp4*xyzzyaabd4
xyzzyaacu4=4.d0*xyzzyaabd4
xyzzyaacf4=2.d0*xyzzyaacg4*xyzzyaaax4-9.d0
xyzzyaacv4=3.d0*xyzzyaabb4
xyzzyaacw4=3.d0*xyzzyaabc4
xyzzyaabr4=xyzzyaacv4+xyzzyaacw4-(2.d0*xyzzyaabd4)
xyzzyaabq4=xyzzyaacg4*xyzzyaabr4-3.d0
xyzzyaacm4(1)=xyzzyaabg4*xyzzyaabq4
xyzzyaacm4(2)=xyzzyaabf4*xyzzyaabq4
xyzzyaacl4(1)=(-xyzzyaacg4)*xyzzyaaba4*xyzzyaabr4*xyzzyaacf4
xyzzyaabr4=6.d0*xyzzyaabd4
xyzzyaabs4=xyzzyaacg4*xyzzyaabr4-3.d0
xyzzyaabt4=2.d0*xyzzyaact4*xyzzyaabd4
xyzzyaacm4(3)=0.5d0*(xyzzyaabr4-xyzzyaabt4+xyzzyaabj4*xyzzyaabs4)
xyzzyaabu4=xyzzyaacr4*xyzzyaabb4-xyzzyaabc4+xyzzyaacu4
xyzzyaabv4=xyzzyaacp4*(xyzzyaabc4-xyzzyaacu4)-3
xyzzyaacm4(4)=1.5d0*(xyzzyaabu4+xyzzyaabb4*xyzzyaabv4)
xyzzyaabv4=xyzzyaacq4*(xyzzyaabj4-xyzzyaacu4)
xyzzyaacm4(5)=xyzzyaabe4*(xyzzyaabv4-3.d0)
xyzzyaacm4(6)=xyzzyaabg4*(xyzzyaabv4+12.d0)
xyzzyaacl4(2)=(-xyzzyaaay4)*xyzzyaabv4*xyzzyaacf4
xyzzyaacm4(7)=xyzzyaacm4(5)
xyzzyaacm4(9)=xyzzyaabf4*(xyzzyaabv4+12.d0)
xyzzyaabw4=xyzzyaacs4*xyzzyaabc4+xyzzyaabb4*(xyzzyaacs4-1.d0)
xyzzyaabx4=xyzzyaabc4*(3.d0+8.d0*xyzzyaacg4*xyzzyaabd4)
xyzzyaacm4(8)=1.5d0*(xyzzyaabw4-xyzzyaabx4+xyzzyaacu4)
xyzzyaacl4(3)=(-xyzzyaaaz4)*xyzzyaabv4*xyzzyaacf4
xyzzyaaby4=xyzzyaacg4*xyzzyaabk4
xyzzyaacm4(10)=30.d0*xyzzyaabg4*(xyzzyaaby4+1.d0)
xyzzyaacm4(11)=30.d0*xyzzyaabf4*(xyzzyaaby4-1.d0)
xyzzyaacm4(12)=15.d0*xyzzyaabk4*(xyzzyaact4-1.d0)
xyzzyaacl4(4)=-30.d0*xyzzyaaba4*xyzzyaaby4*xyzzyaacf4
xyzzyaacm4(13)=30.d0*(1.d0-xyzzyaacr4)*xyzzyaabf4
xyzzyaacm4(14)=30.d0*(1.d0-xyzzyaacs4)*xyzzyaabg4
xyzzyaacm4(15)=30.d0*(1.d0-xyzzyaact4)*xyzzyaabe4
xyzzyaacl4(5)=30.d0*xyzzyaacp4*xyzzyaacx4*xyzzyaacf4
xyzzyaabz4=xyzzyaacv4*(1+xyzzyaacs4)
xyzzyaaca4=xyzzyaacr4*xyzzyaabb4+xyzzyaacw4
xyzzyaacm4(16)=-15.d0*(xyzzyaaca4-xyzzyaabz4)
xyzzyaacb4=(-30.d0)*xyzzyaacg4*(xyzzyaabb4-xyzzyaacw4)
xyzzyaacm4(17)=xyzzyaabe4*(xyzzyaacb4-90.d0)
xyzzyaacm4(18)=xyzzyaabg4*xyzzyaacb4
xyzzyaacl4(6)=(-xyzzyaaay4)*xyzzyaacb4*xyzzyaacf4
xyzzyaacc4=30.d0*xyzzyaacg4*(xyzzyaabc4-xyzzyaacv4)
xyzzyaacd4=xyzzyaabb4*(6.d0*xyzzyaacg4*xyzzyaabc4-3.d0)
xyzzyaace4=xyzzyaabc4*(xyzzyaacs4-3.d0)
xyzzyaacm4(19)=xyzzyaabe4*(90.d0+xyzzyaacc4)
xyzzyaacm4(20)=15.d0*(xyzzyaace4-xyzzyaacd4)
xyzzyaacm4(21)=xyzzyaabf4*xyzzyaacc4
xyzzyaacl4(7)=(-xyzzyaaaz4)*xyzzyaacc4*xyzzyaacf4
xyzzyaack4=c_prim(xyzzyaaap4)*xyzzyaach4
do xyzzyaaan4=1,21
xyzzyaacm4(xyzzyaaan4)=xyzzyaacm4(xyzzyaaan4)*xyzzyaack4
enddo
do xyzzyaaan4=1,7
xyzzyaacl4(xyzzyaaan4)=xyzzyaacl4(xyzzyaaan4)*xyzzyaack4
enddo
end select
do xyzzyaaan4=1,xyzzyaaah4
xyzzyaacn4(xyzzyaaan4)=xyzzyaacn4(xyzzyaaan4)+xyzzyaacl4(xyzzyaaan4)
enddo
do xyzzyaaan4=1,xyzzyaaai4
xyzzyaaco4(xyzzyaaan4)=xyzzyaaco4(xyzzyaaan4)+xyzzyaacm4(xyzzyaaan4)
enddo
enddo
if(xyzzyaacz4)then
xyzzyaaat4=first_ao(xyzzyaaaf4)
call gmult(rblap(xyzzyaaat4,1),num_ao,num_real_k,xyzzyaacn4,phase,xyzz&
&yaaah4)
call g3mult(rbgra1(xyzzyaaat4,1),rbgra2(xyzzyaaat4,1),rbgra3(xyzzyaaat&
&4,1),num_ao,num_real_k,xyzzyaaco4,phase,xyzzyaaah4)
xyzzyaaau4=xyzzyaaat4+xyzzyaaah4-1
do xyzzyaaab4=num_real_k_plus_1,num_k
xyzzyaacy4=cphase(xyzzyaaab4)
cblap (xyzzyaaat4:xyzzyaaau4,xyzzyaaab4)=cblap (xyzzyaaat4:xyzzyaaau4,&
&xyzzyaaab4)+xyzzyaacn4(1:xyzzyaaah4)*xyzzyaacy4
cbgra1(xyzzyaaat4:xyzzyaaau4,xyzzyaaab4)=cbgra1(xyzzyaaat4:xyzzyaaau4,&
&xyzzyaaab4)+xyzzyaaco4(1:xyzzyaaai4:3)*xyzzyaacy4
cbgra2(xyzzyaaat4:xyzzyaaau4,xyzzyaaab4)=cbgra2(xyzzyaaat4:xyzzyaaau4,&
&xyzzyaaab4)+xyzzyaaco4(2:xyzzyaaai4:3)*xyzzyaacy4
cbgra3(xyzzyaaat4:xyzzyaaau4,xyzzyaaab4)=cbgra3(xyzzyaaat4:xyzzyaaau4,&
&xyzzyaaab4)+xyzzyaaco4(3:xyzzyaaai4:3)*xyzzyaacy4
enddo
endif
enddo
enddo
enddo
xyzzyaaav4=0
do xyzzyaaaf4=1,num_shells
if(.not.xyzzyaada4(xyzzyaaaf4))cycle
xyzzyaaat4=first_ao(xyzzyaaaf4)-1
do xyzzyaaan4=1,numao_in_shell(xyzzyaaaf4)
xyzzyaaav4=xyzzyaaav4+1
xyzzyaaaw4(xyzzyaaav4)=xyzzyaaat4+xyzzyaaan4
enddo
enddo
if(use_gpcc)call cusp_wfdet_cmplx(rvec,jspin,.false.,.true.,cusp_val_c&
&,cusp_grad_c,cusp_lap_c)
xyzzyaaal4=0
xyzzyaaam4=0
xyzzyaaaj4=gauss_offset(jspin)
xyzzyaaak4=gauss_offsetc(jspin)
do xyzzyaaab4=1,num_real_k
xyzzyaaad4=nband(xyzzyaaab4,xyzzyaaaa1)
orb2(xyzzyaaal4+1:xyzzyaaal4+xyzzyaaad4)=multi_ddot_s(xyzzyaaad4,xyzzy&
&aaav4,xyzzyaaaw4,corbmask(xyzzyaaam4+xyzzyaaak4),rck(:,:,xyzzyaaab4,x&
&yzzyaaab1),maxb*num_ao,maxb,rblap(1,xyzzyaaab4),rnorm)
orb3(xyzzyaaal4+1:xyzzyaaal4+xyzzyaaad4)=multi_ddot_s(xyzzyaaad4,xyzzy&
&aaav4,xyzzyaaaw4,corbmask(xyzzyaaam4+xyzzyaaak4),rck(:,:,xyzzyaaab4,x&
&yzzyaaab1),maxb*num_ao,maxb,rbgra1(1,xyzzyaaab4),rnorm)
orb4(xyzzyaaal4+1:xyzzyaaal4+xyzzyaaad4)=multi_ddot_s(xyzzyaaad4,xyzzy&
&aaav4,xyzzyaaaw4,corbmask(xyzzyaaam4+xyzzyaaak4),rck(:,:,xyzzyaaab4,x&
&yzzyaaab1),maxb*num_ao,maxb,rbgra2(1,xyzzyaaab4),rnorm)
orb5(xyzzyaaal4+1:xyzzyaaal4+xyzzyaaad4)=multi_ddot_s(xyzzyaaad4,xyzzy&
&aaav4,xyzzyaaaw4,corbmask(xyzzyaaam4+xyzzyaaak4),rck(:,:,xyzzyaaab4,x&
&yzzyaaab1),maxb*num_ao,maxb,rbgra3(1,xyzzyaaab4),rnorm)
xyzzyaaal4=xyzzyaaal4+xyzzyaaad4
xyzzyaaam4=xyzzyaaam4+xyzzyaaad4
enddo
if(use_gpcc)then
do xyzzyaaan4=1,xyzzyaaal4
if(.not.corbmask(xyzzyaaan4+xyzzyaaak4-1))cycle
orb2(xyzzyaaan4)=orb2(xyzzyaaan4)+dble(cusp_lap_c(xyzzyaaan4+xyzzyaaak&
&4-1))
orb3(xyzzyaaan4)=orb3(xyzzyaaan4)+dble(cusp_grad_c(1,xyzzyaaan4+xyzzya&
&aak4-1))
orb4(xyzzyaaan4)=orb4(xyzzyaaan4)+dble(cusp_grad_c(2,xyzzyaaan4+xyzzya&
&aak4-1))
orb5(xyzzyaaan4)=orb5(xyzzyaaan4)+dble(cusp_grad_c(3,xyzzyaaan4+xyzzya&
&aak4-1))
enddo
endif
if(complex_states)then
do xyzzyaaab4=num_real_k_plus_1,num_k
xyzzyaaad4=nband(xyzzyaaab4,xyzzyaaaa1)
corb2(1:xyzzyaaad4)=multi_zdotu_s(xyzzyaaad4,xyzzyaaav4,xyzzyaaaw4,cor&
&bmask(xyzzyaaam4+xyzzyaaak4),cck(:,:,xyzzyaaab4,xyzzyaaab1),maxb*num_&
&ao,maxb,cblap(1,xyzzyaaab4),two_rnorm)
corb3(1:xyzzyaaad4)=multi_zdotu_s(xyzzyaaad4,xyzzyaaav4,xyzzyaaaw4,cor&
&bmask(xyzzyaaam4+xyzzyaaak4),cck(:,:,xyzzyaaab4,xyzzyaaab1),maxb*num_&
&ao,maxb,cbgra1(1,xyzzyaaab4),two_rnorm)
corb4(1:xyzzyaaad4)=multi_zdotu_s(xyzzyaaad4,xyzzyaaav4,xyzzyaaaw4,cor&
&bmask(xyzzyaaam4+xyzzyaaak4),cck(:,:,xyzzyaaab4,xyzzyaaab1),maxb*num_&
&ao,maxb,cbgra2(1,xyzzyaaab4),two_rnorm)
corb5(1:xyzzyaaad4)=multi_zdotu_s(xyzzyaaad4,xyzzyaaav4,xyzzyaaaw4,cor&
&bmask(xyzzyaaam4+xyzzyaaak4),cck(:,:,xyzzyaaab4,xyzzyaaab1),maxb*num_&
&ao,maxb,cbgra3(1,xyzzyaaab4),two_rnorm)
do xyzzyaaan4=1,xyzzyaaad4
xyzzyaaal4=xyzzyaaal4+1
xyzzyaaam4=xyzzyaaam4+1
if(.not.corbmask(xyzzyaaam4+xyzzyaaak4-1))cycle
if(use_gpcc)then
corb2(xyzzyaaan4)=corb2(xyzzyaaan4)+cusp_lap_c(xyzzyaaam4+xyzzyaaak4-1&
&)
corb3(xyzzyaaan4)=corb3(xyzzyaaan4)+cusp_grad_c(1,xyzzyaaam4+xyzzyaaak&
&4-1)
corb4(xyzzyaaan4)=corb4(xyzzyaaan4)+cusp_grad_c(2,xyzzyaaam4+xyzzyaaak&
&4-1)
corb5(xyzzyaaan4)=corb5(xyzzyaaan4)+cusp_grad_c(3,xyzzyaaam4+xyzzyaaak&
&4-1)
endif
orb2(xyzzyaaal4)=dble(corb2(xyzzyaaan4))
orb3(xyzzyaaal4)=dble(corb3(xyzzyaaan4))
orb4(xyzzyaaal4)=dble(corb4(xyzzyaaan4))
orb5(xyzzyaaal4)=dble(corb5(xyzzyaaan4))
if(xyzzyaaal4+1<=nuc_nele(jspin))then
xyzzyaaal4=xyzzyaaal4+1
orb2(xyzzyaaal4)=aimag(corb2(xyzzyaaan4))
orb3(xyzzyaaal4)=aimag(corb3(xyzzyaaan4))
orb4(xyzzyaaal4)=aimag(corb4(xyzzyaaan4))
orb5(xyzzyaaal4)=aimag(corb5(xyzzyaaan4))
endif
enddo
enddo
endif
call dcopy(xyzzyaaal4,orb2(1),1,orblap(xyzzyaaaj4,1),1)
call dcopy(xyzzyaaal4,orb3(1),1,orbgrad(1,xyzzyaaaj4,1),3)
call dcopy(xyzzyaaal4,orb4(1),1,orbgrad(2,xyzzyaaaj4,1),3)
call dcopy(xyzzyaaal4,orb5(1),1,orbgrad(3,xyzzyaaaj4,1),3)
if(excite)then
do xyzzyaaan4=1,gauss_ex_norb
xyzzyaaao4=gauss_gs_norb+xyzzyaaan4
if(.not.orbmask(xyzzyaaao4))cycle
xyzzyaaab4=virtual_k(xyzzyaaan4,jspin)
xyzzyaaad4=full2vrt(xyzzyaaan4,jspin)
if(xyzzyaaab4<=num_real_k)then
flap(xyzzyaaan4)=sum(rblap(:,xyzzyaaab4)*rck_ex(:,xyzzyaaad4,jspin))*r&
&norm
fgra1(xyzzyaaan4)=sum(rbgra1(:,xyzzyaaab4)*rck_ex(:,xyzzyaaad4,jspin))&
&*rnorm
fgra2(xyzzyaaan4)=sum(rbgra2(:,xyzzyaaab4)*rck_ex(:,xyzzyaaad4,jspin))&
&*rnorm
fgra3(xyzzyaaan4)=sum(rbgra3(:,xyzzyaaab4)*rck_ex(:,xyzzyaaad4,jspin))&
&*rnorm
else
flap(xyzzyaaan4)=sum(dble(cblap(:,xyzzyaaab4)*cck_ex(:,xyzzyaaad4,jspi&
&n)))*two_rnorm
fgra1(xyzzyaaan4)=sum(dble(cbgra1(:,xyzzyaaab4)*cck_ex(:,xyzzyaaad4,js&
&pin)))*two_rnorm
fgra2(xyzzyaaan4)=sum(dble(cbgra2(:,xyzzyaaab4)*cck_ex(:,xyzzyaaad4,js&
&pin)))*two_rnorm
fgra3(xyzzyaaan4)=sum(dble(cbgra3(:,xyzzyaaab4)*cck_ex(:,xyzzyaaad4,js&
&pin)))*two_rnorm
endif
if(use_gpcc)then
flap(xyzzyaaan4)=flap(xyzzyaaan4)+dble(cusp_lap_c(ridx2zidx(xyzzyaaao4&
&)))
fgra1(xyzzyaaan4)=fgra1(xyzzyaaan4)+dble(cusp_grad_c(1,ridx2zidx(xyzzy&
&aaao4)))
fgra2(xyzzyaaan4)=fgra2(xyzzyaaan4)+dble(cusp_grad_c(2,ridx2zidx(xyzzy&
&aaao4)))
fgra3(xyzzyaaan4)=fgra3(xyzzyaaan4)+dble(cusp_grad_c(3,ridx2zidx(xyzzy&
&aaao4)))
endif
orblap(xyzzyaaao4,1)=flap(xyzzyaaan4)
orbgrad(1,xyzzyaaao4,1)=fgra1(xyzzyaaan4)
orbgrad(2,xyzzyaaao4,1)=fgra2(xyzzyaaan4)
orbgrad(3,xyzzyaaao4,1)=fgra3(xyzzyaaan4)
enddo
endif
end subroutine xyzzyaaad1
subroutine xyzzyaaae1(rvec,jspin,norb,orbmask,norbc,corbmask,orbval,or&
&bgrad,orblap)
implicit none
integer,intent(in) :: jspin,norb,norbc
logical,intent(in) :: orbmask(norb),corbmask(norbc)
real(dp),intent(in) :: rvec(3)
real(dp),intent(inout) :: orbval(norb,real1_complex2),orbgrad(3,norb,r&
&eal1_complex2),orblap(norb,real1_complex2)
integer xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5,xyzzyaaad5,xyzzyaaae5,xyzzyaa&
&af5,xyzzyaaag5,xyzzyaaah5,xyzzyaaai5,xyzzyaaaj5,xyzzyaaak5,xyzzyaaal5&
&,xyzzyaaam5,xyzzyaaan5,xyzzyaaao5,xyzzyaaap5,xyzzyaaaq5,xyzzyaaar5,xy&
&zzyaaas5,xyzzyaaat5,xyzzyaaau5,xyzzyaaav5,xyzzyaaaw5(num_ao)
logical xyzzyaaax5,xyzzyaaay5(num_shells)
real(dp) xyzzyaaaz5,xyzzyaaba5,xyzzyaabb5,xyzzyaabc5,xyzzyaabd5,xyzzya&
&abe5,xyzzyaabf5,xyzzyaabg5,xyzzyaabh5,xyzzyaabi5,xyzzyaabj5,xyzzyaabk&
&5,xyzzyaabl5,xyzzyaabm5,xyzzyaabn5,xyzzyaabo5,xyzzyaabp5,xyzzyaabq5,x&
&yzzyaabr5,xyzzyaabs5,xyzzyaabt5,xyzzyaabu5,xyzzyaabv5,xyzzyaabw5,xyzz&
&yaabx5,xyzzyaaby5,xyzzyaabz5,xyzzyaaca5,xyzzyaacb5,xyzzyaacc5,xyzzyaa&
&cd5,xyzzyaace5,xyzzyaacf5,xyzzyaacg5,xyzzyaach5,xyzzyaaci5,xyzzyaacj5&
&,xyzzyaack5,xyzzyaacl5,xyzzyaacm5,xyzzyaacn5,xyzzyaaco5,xyzzyaacp5,xy&
&zzyaacq5(7),xyzzyaacr5(7),xyzzyaacs5(21),xyzzyaact5(7),xyzzyaacu5(7),&
&xyzzyaacv5(21),xyzzyaacw5,xyzzyaacx5,xyzzyaacy5,xyzzyaacz5,xyzzyaada5&
&,xyzzyaadb5,xyzzyaadc5,xyzzyaadd5,xyzzyaade5,xyzzyaadf5,xyzzyaadg5,xy&
&zzyaadh5,xyzzyaadi5,xyzzyaadj5,xyzzyaadk5
complex(dp) xyzzyaadl5
xyzzyaaar5=0
rbf=0.d0
rblap=0.d0
rbgra1=0.d0
rbgra2=0.d0
rbgra3=0.d0
xyzzyaaay5=.false.
if(complex_states)then
cbf=czero
cblap=czero
cbgra1=czero
cbgra2=czero
cbgra3=czero
endif
do xyzzyaaaa5=1,num_cell
do xyzzyaaab5=1,num_real_k
phase(xyzzyaaab5)=coskdotg(xyzzyaaab5,xyzzyaaaa5)*coskdotg_shift(xyzzy&
&aaab5)
enddo
do xyzzyaaab5=num_real_k_plus_1,num_k
cphase(xyzzyaaab5)=expikdotg(xyzzyaaab5,xyzzyaaaa5)*expikdotg_shift(xy&
&zzyaaab5)
enddo
do xyzzyaaaq5=1,num_centres_in_cell(xyzzyaaaa5)
xyzzyaaba5=xpos_in_cell(xyzzyaaaq5,xyzzyaaaa5)
xyzzyaabb5=ypos_in_cell(xyzzyaaaq5,xyzzyaaaa5)
xyzzyaabc5=zpos_in_cell(xyzzyaaaq5,xyzzyaaaa5)
xyzzyaaba5=new_origin_to_rvec(1)-xyzzyaaba5
xyzzyaabb5=new_origin_to_rvec(2)-xyzzyaabb5
xyzzyaabc5=new_origin_to_rvec(3)-xyzzyaabc5
xyzzyaabd5=xyzzyaaba5*xyzzyaaba5
xyzzyaabe5=xyzzyaabb5*xyzzyaabb5
xyzzyaabf5=xyzzyaabc5*xyzzyaabc5
xyzzyaaaz5=xyzzyaabd5+xyzzyaabe5+xyzzyaabf5
xyzzyaaau5=num_shells_on_centre(xyzzyaaaq5,xyzzyaaaa5)
do xyzzyaaac5=1,xyzzyaaau5
xyzzyaaar5=xyzzyaaar5+1
xyzzyaaap5=shell_sequence_number(xyzzyaaar5)
xyzzyaacl5=min_exponent(xyzzyaaap5)
xyzzyaabu5=xyzzyaacl5*xyzzyaaaz5
if(xyzzyaabu5>screening_tolerance)cycle
xyzzyaaad5=shell_am(xyzzyaaap5)
xyzzyaaae5=numao_in_shell(xyzzyaaap5)
xyzzyaaaf5=xyzzyaaae5*3
do xyzzyaaai5=1,xyzzyaaae5
xyzzyaact5(xyzzyaaai5)=0.d0
enddo
do xyzzyaaai5=1,xyzzyaaae5
xyzzyaacu5(xyzzyaaai5)=0.d0
enddo
do xyzzyaaai5=1,xyzzyaaaf5
xyzzyaacv5(xyzzyaaai5)=0.d0
enddo
xyzzyaaan5=primitive(xyzzyaaap5)
xyzzyaaao5=primitive(xyzzyaaap5+1)-1
xyzzyaaax5=.false.
do xyzzyaaak5=xyzzyaaan5,xyzzyaaao5
xyzzyaacl5=exponent(xyzzyaaak5)
xyzzyaabu5=xyzzyaacl5*xyzzyaaaz5
if(xyzzyaabu5>30.d0)cycle
xyzzyaaax5=.true.
xyzzyaaay5(xyzzyaaap5)=.true.
xyzzyaacm5=exp(-xyzzyaabu5)
xyzzyaabq5=-xyzzyaacl5-xyzzyaacl5
select case (xyzzyaaad5)
case(1)
xyzzyaacq5(1)=c_prim(xyzzyaaak5)*xyzzyaacm5
xyzzyaacn5=xyzzyaacq5(1)*xyzzyaabq5
xyzzyaacs5(1)=xyzzyaaba5*xyzzyaacn5
xyzzyaacs5(2)=xyzzyaabb5*xyzzyaacn5
xyzzyaacs5(3)=xyzzyaabc5*xyzzyaacn5
xyzzyaabu5=3.d0+xyzzyaabq5*xyzzyaaaz5
xyzzyaacr5(1)=xyzzyaabu5*xyzzyaacn5
case(2)
xyzzyaacn5=c_prim(xyzzyaaak5)*xyzzyaacm5
xyzzyaaco5=c_prim2(xyzzyaaak5)*xyzzyaacm5
xyzzyaacq5(1)=xyzzyaacn5
xyzzyaacq5(2)=xyzzyaaco5*xyzzyaaba5
xyzzyaacq5(3)=xyzzyaaco5*xyzzyaabb5
xyzzyaacq5(4)=xyzzyaaco5*xyzzyaabc5
xyzzyaacs5(1)=xyzzyaabq5*xyzzyaaba5
xyzzyaacs5(2)=xyzzyaabq5*xyzzyaabb5
xyzzyaacs5(3)=xyzzyaabq5*xyzzyaabc5
xyzzyaacs5(4)=1.d0+xyzzyaacs5(1)*xyzzyaaba5
xyzzyaacs5(5)=xyzzyaacs5(1)*xyzzyaabb5
xyzzyaacs5(6)=xyzzyaacs5(1)*xyzzyaabc5
xyzzyaacs5(8)=1.d0+xyzzyaacs5(2)*xyzzyaabb5
xyzzyaacs5(9)=xyzzyaacs5(2)*xyzzyaabc5
xyzzyaacs5(7)=xyzzyaacs5(5)
xyzzyaacs5(10)=xyzzyaacs5(6)
xyzzyaacs5(11)=xyzzyaacs5(9)
xyzzyaacs5(12)=1.d0+xyzzyaacs5(3)*xyzzyaabc5
xyzzyaacs5(1:3)=xyzzyaacs5(1:3)*xyzzyaacn5
xyzzyaacs5(4:12)=xyzzyaacs5(4:12)*xyzzyaaco5
xyzzyaabu5=3.d0+xyzzyaabq5*xyzzyaaaz5
xyzzyaabv5=xyzzyaabq5*xyzzyaabu5
xyzzyaabw5=(xyzzyaabv5+2.d0*xyzzyaabq5)*xyzzyaaco5
xyzzyaacr5(1)=xyzzyaabv5*xyzzyaacn5
xyzzyaacr5(2)=xyzzyaabw5*xyzzyaaba5
xyzzyaacr5(3)=xyzzyaabw5*xyzzyaabb5
xyzzyaacr5(4)=xyzzyaabw5*xyzzyaabc5
case(3)
xyzzyaaco5=c_prim(xyzzyaaak5)*xyzzyaacm5
xyzzyaacq5(1)=xyzzyaaco5*xyzzyaaba5
xyzzyaacq5(2)=xyzzyaaco5*xyzzyaabb5
xyzzyaacq5(3)=xyzzyaaco5*xyzzyaabc5
xyzzyaabr5=xyzzyaabq5*xyzzyaaba5
xyzzyaabs5=xyzzyaabq5*xyzzyaabb5
xyzzyaabt5=xyzzyaabq5*xyzzyaabc5
xyzzyaabx5=xyzzyaabq5*(5+xyzzyaabq5*xyzzyaaaz5)
xyzzyaacs5(1)=1.d0+xyzzyaabr5*xyzzyaaba5
xyzzyaacs5(2)=xyzzyaabr5*xyzzyaabb5
xyzzyaacs5(3)=xyzzyaabr5*xyzzyaabc5
xyzzyaacs5(5)=1.d0+xyzzyaabs5*xyzzyaabb5
xyzzyaacs5(6)=xyzzyaabs5*xyzzyaabc5
xyzzyaacs5(9)=1.d0+xyzzyaabt5*xyzzyaabc5
xyzzyaacs5(4)=xyzzyaacs5(2)
xyzzyaacs5(7)=xyzzyaacs5(3)
xyzzyaacs5(8)=xyzzyaacs5(6)
xyzzyaacr5(1)=xyzzyaabx5*xyzzyaaba5
xyzzyaacr5(2)=xyzzyaabx5*xyzzyaabb5
xyzzyaacr5(3)=xyzzyaabx5*xyzzyaabc5
do xyzzyaaai5=1,9
xyzzyaacs5(xyzzyaaai5)=xyzzyaacs5(xyzzyaaai5)*xyzzyaaco5
enddo
do xyzzyaaai5=1,3
xyzzyaacr5(xyzzyaaai5)=xyzzyaacr5(xyzzyaaai5)*xyzzyaaco5
enddo
case(4)
xyzzyaabg5=xyzzyaaba5*xyzzyaabb5
xyzzyaabh5=xyzzyaabb5*xyzzyaabc5
xyzzyaabi5=xyzzyaaba5*xyzzyaabc5
xyzzyaabk5=xyzzyaabd5-xyzzyaabe5
xyzzyaabj5=3.d0*xyzzyaabf5-xyzzyaaaz5
xyzzyaacq5(1)=xyzzyaabj5
xyzzyaacq5(2)=xyzzyaabi5
xyzzyaacq5(3)=xyzzyaabh5
xyzzyaacq5(4)=xyzzyaabk5
xyzzyaacq5(5)=xyzzyaabg5
xyzzyaabr5=xyzzyaabq5*xyzzyaaba5
xyzzyaabs5=xyzzyaabq5*xyzzyaabb5
xyzzyaabt5=xyzzyaabq5*xyzzyaabc5
xyzzyaabu5=1.d0+xyzzyaabr5*xyzzyaaba5
xyzzyaabv5=1.d0+xyzzyaabs5*xyzzyaabb5
xyzzyaabw5=1.d0+xyzzyaabt5*xyzzyaabc5
xyzzyaabx5=xyzzyaabq5*xyzzyaabj5-2.d0
xyzzyaaby5=xyzzyaabq5*xyzzyaabk5
xyzzyaacs5(1)=xyzzyaabx5*xyzzyaaba5
xyzzyaacs5(2)=xyzzyaabx5*xyzzyaabb5
xyzzyaacs5(3)=(xyzzyaabx5+6.d0)*xyzzyaabc5
xyzzyaacs5(4)=xyzzyaabu5*xyzzyaabc5
xyzzyaacs5(5)=xyzzyaabr5*xyzzyaabh5
xyzzyaacs5(6)=xyzzyaabw5*xyzzyaaba5
xyzzyaacs5(7)=xyzzyaacs5(5)
xyzzyaacs5(8)=xyzzyaabv5*xyzzyaabc5
xyzzyaacs5(9)=xyzzyaabw5*xyzzyaabb5
xyzzyaacs5(10)=(2.d0+xyzzyaabq5*xyzzyaabk5)*xyzzyaaba5
xyzzyaacs5(11)=(xyzzyaaby5-2.d0)*xyzzyaabb5
xyzzyaacs5(12)=xyzzyaaby5*xyzzyaabc5
xyzzyaacs5(13)=xyzzyaabu5*xyzzyaabb5
xyzzyaacs5(14)=xyzzyaabv5*xyzzyaaba5
xyzzyaacs5(15)=xyzzyaacs5(5)
xyzzyaacp5=c_prim(xyzzyaaak5)*xyzzyaacm5
xyzzyaacs5(1:15)=xyzzyaacs5(1:15)*xyzzyaacp5
xyzzyaacq5(1:5)=xyzzyaacq5(1:5)*xyzzyaacp5
xyzzyaabu5=7.d0+xyzzyaabq5*xyzzyaaaz5
xyzzyaabu5=xyzzyaabq5*xyzzyaabu5*xyzzyaacp5
xyzzyaacr5(1)=xyzzyaabu5*xyzzyaabj5
xyzzyaacr5(2)=xyzzyaabu5*xyzzyaabi5
xyzzyaacr5(3)=xyzzyaabu5*xyzzyaabh5
xyzzyaacr5(4)=xyzzyaabu5*xyzzyaabk5
xyzzyaacr5(5)=xyzzyaabu5*xyzzyaabg5
case(5)
xyzzyaabn5=xyzzyaabd5*xyzzyaaba5
xyzzyaabo5=xyzzyaabe5*xyzzyaabb5
xyzzyaabp5=xyzzyaabf5*xyzzyaabc5
xyzzyaabg5=xyzzyaaba5*xyzzyaabb5
xyzzyaabi5=xyzzyaaba5*xyzzyaabc5
xyzzyaabh5=xyzzyaabb5*xyzzyaabc5
xyzzyaade5=xyzzyaabg5*xyzzyaaba5
xyzzyaadf5=xyzzyaabg5*xyzzyaabb5
xyzzyaadg5=xyzzyaabi5*xyzzyaaba5
xyzzyaadh5=xyzzyaabi5*xyzzyaabc5
xyzzyaadj5=xyzzyaabh5*xyzzyaabb5
xyzzyaadi5=xyzzyaabh5*xyzzyaabc5
xyzzyaadk5=xyzzyaabg5*xyzzyaabc5
xyzzyaabl5=xyzzyaabd5+xyzzyaabe5
xyzzyaabm5=xyzzyaabe5-xyzzyaabd5
xyzzyaadc5=3.d0*xyzzyaabd5
xyzzyaadd5=3.d0*xyzzyaabe5
xyzzyaacw5=2.d0*xyzzyaacl5
xyzzyaacx5=3.d0*xyzzyaacl5
xyzzyaacy5=xyzzyaacw5*xyzzyaabd5
xyzzyaacz5=xyzzyaacw5*xyzzyaabe5
xyzzyaada5=xyzzyaacw5*xyzzyaabf5
xyzzyaadb5=4.d0*xyzzyaabf5
xyzzyaacq5(1)=xyzzyaabp5-1.5d0*(xyzzyaadg5+xyzzyaadj5)
xyzzyaacq5(2)=6.d0*xyzzyaadh5-1.5d0*(xyzzyaabn5+xyzzyaadf5)
xyzzyaacq5(3)=6.d0*xyzzyaadi5-1.5d0*(xyzzyaade5+xyzzyaabo5)
xyzzyaacq5(4)=15.d0*(xyzzyaadg5-xyzzyaadj5)
xyzzyaacq5(5)=30.d0*xyzzyaadk5
xyzzyaacq5(6)=15.d0*xyzzyaabn5-45.d0*xyzzyaadf5
xyzzyaacq5(7)=45.d0*xyzzyaade5-15.d0*xyzzyaabo5
xyzzyaacp5=c_prim(xyzzyaaak5)*xyzzyaacm5
xyzzyaacq5(1:7)=xyzzyaacq5(1:7)*xyzzyaacp5
xyzzyaack5=xyzzyaacw5*(xyzzyaacw5*xyzzyaaaz5-9.d0)
xyzzyaacr5(1:7)=xyzzyaacq5(1:7)*xyzzyaack5
xyzzyaabv5=xyzzyaadc5+xyzzyaadd5-(2.d0*xyzzyaabf5)
xyzzyaabv5=xyzzyaacl5*xyzzyaabv5-3.d0
xyzzyaacs5(1)=xyzzyaabi5*xyzzyaabv5
xyzzyaacs5(2)=xyzzyaabh5*xyzzyaabv5
xyzzyaabw5=6.d0*xyzzyaabf5
xyzzyaabx5=xyzzyaacl5*xyzzyaabw5-3.d0
xyzzyaaby5=2.d0*xyzzyaada5*xyzzyaabf5
xyzzyaacs5(3)=0.5d0*(xyzzyaabw5-xyzzyaaby5+xyzzyaabl5*xyzzyaabx5)
xyzzyaabz5=xyzzyaacy5*xyzzyaabd5-xyzzyaabe5+xyzzyaadb5
xyzzyaaca5=xyzzyaacw5*(xyzzyaabe5-xyzzyaadb5)-3
xyzzyaacs5(4)=1.5d0*(xyzzyaabz5+xyzzyaabd5*xyzzyaaca5)
xyzzyaaca5=xyzzyaacx5*(xyzzyaabl5-xyzzyaadb5)
xyzzyaacs5(5)=xyzzyaabg5*(xyzzyaaca5-3.d0)
xyzzyaacs5(6)=xyzzyaabi5*(xyzzyaaca5+12.d0)
xyzzyaacs5(7)=xyzzyaacs5(5)
xyzzyaacs5(9)=xyzzyaabh5*(xyzzyaaca5+12.d0)
xyzzyaacb5=xyzzyaacz5*xyzzyaabe5+xyzzyaabd5*(xyzzyaacz5-1.d0)
xyzzyaacc5=xyzzyaabe5*(3.d0+8.d0*xyzzyaacl5*xyzzyaabf5)
xyzzyaacs5(8)=1.5d0*(xyzzyaacb5-xyzzyaacc5+xyzzyaadb5)
xyzzyaacd5=xyzzyaacl5*xyzzyaabm5
xyzzyaacs5(10)=30.d0*xyzzyaabi5*(xyzzyaacd5+1.d0)
xyzzyaacs5(11)=30.d0*xyzzyaabh5*(xyzzyaacd5-1.d0)
xyzzyaacs5(12)=15.d0*xyzzyaabm5*(xyzzyaada5-1.d0)
xyzzyaacs5(13)=30.d0*(1.d0-xyzzyaacy5)*xyzzyaabh5
xyzzyaacs5(14)=30.d0*(1.d0-xyzzyaacz5)*xyzzyaabi5
xyzzyaacs5(15)=30.d0*(1.d0-xyzzyaada5)*xyzzyaabg5
xyzzyaace5=xyzzyaadc5*(1+xyzzyaacz5)
xyzzyaacf5=xyzzyaacy5*xyzzyaabd5+xyzzyaadd5
xyzzyaacs5(16)=-15.d0*(xyzzyaacf5-xyzzyaace5)
xyzzyaacg5=(-30.d0)*xyzzyaacl5*(xyzzyaabd5-xyzzyaadd5)
xyzzyaacs5(17)=xyzzyaabg5*(xyzzyaacg5-90.d0)
xyzzyaacs5(18)=xyzzyaabi5*xyzzyaacg5
xyzzyaach5=30.d0*xyzzyaacl5*(xyzzyaabe5-xyzzyaadc5)
xyzzyaaci5=xyzzyaabd5*(6.d0*xyzzyaacl5*xyzzyaabe5-3.d0)
xyzzyaacj5=xyzzyaabe5*(xyzzyaacz5-3.d0)
xyzzyaacs5(19)=xyzzyaabg5*(90.d0+xyzzyaach5)
xyzzyaacs5(20)=15.d0*(xyzzyaacj5-xyzzyaaci5)
xyzzyaacs5(21)=xyzzyaabh5*xyzzyaach5
xyzzyaacs5(1:21)=xyzzyaacs5(1:21)*xyzzyaacp5
end select
do xyzzyaaai5=1,xyzzyaaae5
xyzzyaact5(xyzzyaaai5)=xyzzyaact5(xyzzyaaai5)+xyzzyaacq5(xyzzyaaai5)
enddo
do xyzzyaaai5=1,xyzzyaaae5
xyzzyaacu5(xyzzyaaai5)=xyzzyaacu5(xyzzyaaai5)+xyzzyaacr5(xyzzyaaai5)
enddo
do xyzzyaaai5=1,xyzzyaaaf5
xyzzyaacv5(xyzzyaaai5)=xyzzyaacv5(xyzzyaaai5)+xyzzyaacs5(xyzzyaaai5)
enddo
enddo
if(xyzzyaaax5)then
xyzzyaaas5=first_ao(xyzzyaaap5)
call gmult(rbf(xyzzyaaas5,1),num_ao,num_real_k,xyzzyaact5,phase,xyzzya&
&aae5)
call gmult(rblap(xyzzyaaas5,1),num_ao,num_real_k,xyzzyaacu5,phase,xyzz&
&yaaae5)
call g3mult(rbgra1(xyzzyaaas5,1),rbgra2(xyzzyaaas5,1),rbgra3(xyzzyaaas&
&5,1),num_ao,num_real_k,xyzzyaacv5,phase,xyzzyaaae5)
xyzzyaaat5=xyzzyaaas5+xyzzyaaae5-1
do xyzzyaaab5=num_real_k_plus_1,num_k
xyzzyaadl5=cphase(xyzzyaaab5)
cbf(xyzzyaaas5:xyzzyaaat5,xyzzyaaab5)=cbf(xyzzyaaas5:xyzzyaaat5,xyzzya&
&aab5)+xyzzyaact5(1:xyzzyaaae5)*xyzzyaadl5
cblap(xyzzyaaas5:xyzzyaaat5,xyzzyaaab5)=cblap(xyzzyaaas5:xyzzyaaat5,xy&
&zzyaaab5)+xyzzyaacu5(1:xyzzyaaae5)*xyzzyaadl5
cbgra1(xyzzyaaas5:xyzzyaaat5,xyzzyaaab5)=cbgra1(xyzzyaaas5:xyzzyaaat5,&
&xyzzyaaab5)+xyzzyaacv5(1:xyzzyaaaf5:3)*xyzzyaadl5
cbgra2(xyzzyaaas5:xyzzyaaat5,xyzzyaaab5)=cbgra2(xyzzyaaas5:xyzzyaaat5,&
&xyzzyaaab5)+xyzzyaacv5(2:xyzzyaaaf5:3)*xyzzyaadl5
cbgra3(xyzzyaaas5:xyzzyaaat5,xyzzyaaab5)=cbgra3(xyzzyaaas5:xyzzyaaat5,&
&xyzzyaaab5)+xyzzyaacv5(3:xyzzyaaaf5:3)*xyzzyaadl5
enddo
endif
enddo
enddo
enddo
xyzzyaaav5=0
do xyzzyaaap5=1,num_shells
if(.not.xyzzyaaay5(xyzzyaaap5))cycle
xyzzyaaas5=first_ao(xyzzyaaap5)-1
do xyzzyaaai5=1,numao_in_shell(xyzzyaaap5)
xyzzyaaav5=xyzzyaaav5+1
xyzzyaaaw5(xyzzyaaav5)=xyzzyaaas5+xyzzyaaai5
enddo
enddo
if(use_gpcc)call cusp_wfdet_cmplx(rvec,jspin,.true.,.true.,cusp_val_c,&
&cusp_grad_c,cusp_lap_c)
xyzzyaaag5=0
xyzzyaaah5=0
xyzzyaaal5=gauss_offset(jspin)
xyzzyaaam5=gauss_offsetc(jspin)
do xyzzyaaab5=1,num_real_k
xyzzyaaau5=nband(xyzzyaaab5,xyzzyaaaa1)
orb1(xyzzyaaag5+1:xyzzyaaag5+xyzzyaaau5)=multi_ddot_s(xyzzyaaau5,xyzzy&
&aaav5,xyzzyaaaw5,corbmask(xyzzyaaah5+xyzzyaaam5),rck(:,:,xyzzyaaab5,x&
&yzzyaaab1),maxb*num_ao,maxb,rbf(1,xyzzyaaab5),rnorm)
orb2(xyzzyaaag5+1:xyzzyaaag5+xyzzyaaau5)=multi_ddot_s(xyzzyaaau5,xyzzy&
&aaav5,xyzzyaaaw5,corbmask(xyzzyaaah5+xyzzyaaam5),rck(:,:,xyzzyaaab5,x&
&yzzyaaab1),maxb*num_ao,maxb,rblap(1,xyzzyaaab5),rnorm)
orb3(xyzzyaaag5+1:xyzzyaaag5+xyzzyaaau5)=multi_ddot_s(xyzzyaaau5,xyzzy&
&aaav5,xyzzyaaaw5,corbmask(xyzzyaaah5+xyzzyaaam5),rck(:,:,xyzzyaaab5,x&
&yzzyaaab1),maxb*num_ao,maxb,rbgra1(1,xyzzyaaab5),rnorm)
orb4(xyzzyaaag5+1:xyzzyaaag5+xyzzyaaau5)=multi_ddot_s(xyzzyaaau5,xyzzy&
&aaav5,xyzzyaaaw5,corbmask(xyzzyaaah5+xyzzyaaam5),rck(:,:,xyzzyaaab5,x&
&yzzyaaab1),maxb*num_ao,maxb,rbgra2(1,xyzzyaaab5),rnorm)
orb5(xyzzyaaag5+1:xyzzyaaag5+xyzzyaaau5)=multi_ddot_s(xyzzyaaau5,xyzzy&
&aaav5,xyzzyaaaw5,corbmask(xyzzyaaah5+xyzzyaaam5),rck(:,:,xyzzyaaab5,x&
&yzzyaaab1),maxb*num_ao,maxb,rbgra3(1,xyzzyaaab5),rnorm)
xyzzyaaag5=xyzzyaaag5+xyzzyaaau5
xyzzyaaah5=xyzzyaaah5+xyzzyaaau5
enddo
if(use_gpcc)then
do xyzzyaaai5=1,xyzzyaaag5
if(.not.corbmask(xyzzyaaai5+xyzzyaaam5-1))cycle
orb1(xyzzyaaai5)=orb1(xyzzyaaai5)+dble(cusp_val_c(xyzzyaaai5+xyzzyaaam&
&5-1))
orb2(xyzzyaaai5)=orb2(xyzzyaaai5)+dble(cusp_lap_c(xyzzyaaai5+xyzzyaaam&
&5-1))
orb3(xyzzyaaai5)=orb3(xyzzyaaai5)+dble(cusp_grad_c(1,xyzzyaaai5+xyzzya&
&aam5-1))
orb4(xyzzyaaai5)=orb4(xyzzyaaai5)+dble(cusp_grad_c(2,xyzzyaaai5+xyzzya&
&aam5-1))
orb5(xyzzyaaai5)=orb5(xyzzyaaai5)+dble(cusp_grad_c(3,xyzzyaaai5+xyzzya&
&aam5-1))
enddo
endif
if(complex_states)then
do xyzzyaaab5=num_real_k_plus_1,num_k
xyzzyaaau5=nband(xyzzyaaab5,xyzzyaaaa1)
corb1(1:xyzzyaaau5)=multi_zdotu_s(xyzzyaaau5,xyzzyaaav5,xyzzyaaaw5,cor&
&bmask(xyzzyaaah5+xyzzyaaam5),cck(:,:,xyzzyaaab5,xyzzyaaab1),maxb*num_&
&ao,maxb,cbf(1,xyzzyaaab5),two_rnorm)
corb2(1:xyzzyaaau5)=multi_zdotu_s(xyzzyaaau5,xyzzyaaav5,xyzzyaaaw5,cor&
&bmask(xyzzyaaah5+xyzzyaaam5),cck(:,:,xyzzyaaab5,xyzzyaaab1),maxb*num_&
&ao,maxb,cblap(1,xyzzyaaab5),two_rnorm)
corb3(1:xyzzyaaau5)=multi_zdotu_s(xyzzyaaau5,xyzzyaaav5,xyzzyaaaw5,cor&
&bmask(xyzzyaaah5+xyzzyaaam5),cck(:,:,xyzzyaaab5,xyzzyaaab1),maxb*num_&
&ao,maxb,cbgra1(1,xyzzyaaab5),two_rnorm)
corb4(1:xyzzyaaau5)=multi_zdotu_s(xyzzyaaau5,xyzzyaaav5,xyzzyaaaw5,cor&
&bmask(xyzzyaaah5+xyzzyaaam5),cck(:,:,xyzzyaaab5,xyzzyaaab1),maxb*num_&
&ao,maxb,cbgra2(1,xyzzyaaab5),two_rnorm)
corb5(1:xyzzyaaau5)=multi_zdotu_s(xyzzyaaau5,xyzzyaaav5,xyzzyaaaw5,cor&
&bmask(xyzzyaaah5+xyzzyaaam5),cck(:,:,xyzzyaaab5,xyzzyaaab1),maxb*num_&
&ao,maxb,cbgra3(1,xyzzyaaab5),two_rnorm)
do xyzzyaaai5=1,xyzzyaaau5
xyzzyaaag5=xyzzyaaag5+1
xyzzyaaah5=xyzzyaaah5+1
if(.not.corbmask(xyzzyaaah5+xyzzyaaam5-1))cycle
if(use_gpcc)then
corb1(xyzzyaaai5)=corb1(xyzzyaaai5)+cusp_val_c(xyzzyaaah5+xyzzyaaam5-1&
&)
corb2(xyzzyaaai5)=corb2(xyzzyaaai5)+cusp_lap_c(xyzzyaaah5+xyzzyaaam5-1&
&)
corb3(xyzzyaaai5)=corb3(xyzzyaaai5)+cusp_grad_c(1,xyzzyaaah5+xyzzyaaam&
&5-1)
corb4(xyzzyaaai5)=corb4(xyzzyaaai5)+cusp_grad_c(2,xyzzyaaah5+xyzzyaaam&
&5-1)
corb5(xyzzyaaai5)=corb5(xyzzyaaai5)+cusp_grad_c(3,xyzzyaaah5+xyzzyaaam&
&5-1)
endif
orb1(xyzzyaaag5)=dble(corb1(xyzzyaaai5))
orb2(xyzzyaaag5)=dble(corb2(xyzzyaaai5))
orb3(xyzzyaaag5)=dble(corb3(xyzzyaaai5))
orb4(xyzzyaaag5)=dble(corb4(xyzzyaaai5))
orb5(xyzzyaaag5)=dble(corb5(xyzzyaaai5))
if(xyzzyaaag5+1<=nuc_nele(jspin))then
xyzzyaaag5=xyzzyaaag5+1
orb1(xyzzyaaag5)=aimag(corb1(xyzzyaaai5))
orb2(xyzzyaaag5)=aimag(corb2(xyzzyaaai5))
orb3(xyzzyaaag5)=aimag(corb3(xyzzyaaai5))
orb4(xyzzyaaag5)=aimag(corb4(xyzzyaaai5))
orb5(xyzzyaaag5)=aimag(corb5(xyzzyaaai5))
endif
enddo
enddo
endif
call dcopy(xyzzyaaag5,orb1(1),1,orbval(xyzzyaaal5,1),1)
call dcopy(xyzzyaaag5,orb2(1),1,orblap(xyzzyaaal5,1),1)
call dcopy(xyzzyaaag5,orb3(1),1,orbgrad(1,xyzzyaaal5,1),3)
call dcopy(xyzzyaaag5,orb4(1),1,orbgrad(2,xyzzyaaal5,1),3)
call dcopy(xyzzyaaag5,orb5(1),1,orbgrad(3,xyzzyaaal5,1),3)
if(excite)then
do xyzzyaaai5=1,gauss_ex_norb
xyzzyaaaj5=gauss_gs_norb+xyzzyaaai5
if(.not.orbmask(xyzzyaaaj5))cycle
xyzzyaaab5=virtual_k(xyzzyaaai5,jspin)
xyzzyaaau5=full2vrt(xyzzyaaai5,jspin)
if(xyzzyaaab5<=num_real_k)then
psi(xyzzyaaai5)=sum(rbf(:,xyzzyaaab5)*rck_ex(:,xyzzyaaau5,jspin))*rnor&
&m
flap(xyzzyaaai5)=sum(rblap(:,xyzzyaaab5)*rck_ex(:,xyzzyaaau5,jspin))*r&
&norm
fgra1(xyzzyaaai5)=sum(rbgra1(:,xyzzyaaab5)*rck_ex(:,xyzzyaaau5,jspin))&
&*rnorm
fgra2(xyzzyaaai5)=sum(rbgra2(:,xyzzyaaab5)*rck_ex(:,xyzzyaaau5,jspin))&
&*rnorm
fgra3(xyzzyaaai5)=sum(rbgra3(:,xyzzyaaab5)*rck_ex(:,xyzzyaaau5,jspin))&
&*rnorm
else
psi(xyzzyaaai5)=sum(dble(cbf(:,xyzzyaaab5)*cck_ex(:,xyzzyaaau5,jspin))&
&)*two_rnorm
flap(xyzzyaaai5)=sum(dble(cblap(:,xyzzyaaab5)*cck_ex(:,xyzzyaaau5,jspi&
&n)))*two_rnorm
fgra1(xyzzyaaai5)=sum(dble(cbgra1(:,xyzzyaaab5)*cck_ex(:,xyzzyaaau5,js&
&pin)))*two_rnorm
fgra2(xyzzyaaai5)=sum(dble(cbgra2(:,xyzzyaaab5)*cck_ex(:,xyzzyaaau5,js&
&pin)))*two_rnorm
fgra3(xyzzyaaai5)=sum(dble(cbgra3(:,xyzzyaaab5)*cck_ex(:,xyzzyaaau5,js&
&pin)))*two_rnorm
endif
if(use_gpcc)then
psi(xyzzyaaai5)=psi(xyzzyaaai5)+dble(cusp_val_c(ridx2zidx(xyzzyaaaj5))&
&)
flap(xyzzyaaai5)=flap(xyzzyaaai5)+dble(cusp_lap_c(ridx2zidx(xyzzyaaaj5&
&)))
fgra1(xyzzyaaai5)=fgra1(xyzzyaaai5)+dble(cusp_grad_c(1,ridx2zidx(xyzzy&
&aaaj5)))
fgra2(xyzzyaaai5)=fgra2(xyzzyaaai5)+dble(cusp_grad_c(2,ridx2zidx(xyzzy&
&aaaj5)))
fgra3(xyzzyaaai5)=fgra3(xyzzyaaai5)+dble(cusp_grad_c(3,ridx2zidx(xyzzy&
&aaaj5)))
endif
orbval(xyzzyaaaj5,1)=psi(xyzzyaaai5)
orblap(xyzzyaaaj5,1)=flap(xyzzyaaai5)
orbgrad(1,xyzzyaaaj5,1)=fgra1(xyzzyaaai5)
orbgrad(2,xyzzyaaaj5,1)=fgra2(xyzzyaaai5)
orbgrad(3,xyzzyaaaj5,1)=fgra3(xyzzyaaai5)
enddo
endif
end subroutine xyzzyaaae1
end module slaarnaba

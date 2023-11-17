module slaarnabi
use slaarnaag,   only : czero,c_one
use dsp,         only : dp
use format_utils,only : wout,i2s,wordwrap
use slaarnabg,    only : rion,iontype,nbasis,isperiodic,painv,pamat
use slaarnabt,   only : dcopy,zcopy,ddot,exp_protect,bracket,gsbracket&
&,brentmin
use slaarnaca,       only : zion,is_ae
use parallel,    only : am_master
use run_control, only : errstop,errwarn,check_alloc
use store,       only : which_ssingle,no_ssingles,pcharge,nspin
implicit none
private
public cusp_wfdet_real,cusp_wfdet_cmplx,initialize_gpcc,setup_gpcc,use&
&_gpcc,determine_poly_outer,determine_poly_inner,naeions_prim,nradgrid&
&,nsphgrid,radgrid,sphgrid,ae_index,spherical_av_real,spherical_av_cmp&
&lx,gpcc_init_buffer,gpcc_finalize_buffer,gpcc_load_from_buffer,gpcc_s&
&ave_to_buffer,eval_spline_reg,setup_spline_reg
integer xyzzyaaaa1,xyzzyaaab1
integer xyzzyaaac1
logical,allocatable :: xyzzyaaad1(:,:)
real(dp),parameter :: xyzzyaaae1=3.25819d0,xyzzyaaaf1=-15.0126d0,xyzzy&
&aaag1=33.7308d0,xyzzyaaah1=-42.8705d0,xyzzyaaai1=31.2276d0,xyzzyaaaj1&
&=-12.1316d0,xyzzyaaak1=1.94692d0
real(dp),allocatable :: xyzzyaaal1(:,:,:),xyzzyaaam1(:,:),xyzzyaaan1(:&
&,:),xyzzyaaao1(:)
real(dp),allocatable :: xyzzyaaap1(:)
integer,allocatable :: xyzzyaaaq1(:)
real(dp),allocatable :: xyzzyaaar1(:,:)
complex(dp),allocatable :: xyzzyaaas1(:,:)
integer xyzzyaaat1,xyzzyaaau1
real(dp) xyzzyaaav1,xyzzyaaaw1,xyzzyaaax1,xyzzyaaay1,xyzzyaaaz1,xyzzya&
&aba1,xyzzyaabb1
integer,parameter :: nradgrid=100
real(dp),allocatable :: radgrid(:,:)
real(dp),allocatable :: xyzzyaabc1(:,:,:),xyzzyaabd1(:,:,:)
integer nsphgrid
real(dp),allocatable :: sphgrid(:,:),xyzzyaabe1(:)
integer,parameter :: xyzzyaabf1=2
logical use_gpcc
integer naeions_prim
integer,allocatable :: ae_index(:)
real(dp),parameter :: xyzzyaabg1=50.d0
integer,parameter :: xyzzyaabh1=nradgrid+1
logical xyzzyaabi1
logical xyzzyaabj1
real(dp),parameter :: xyzzyaabk1=0.9d0
real(dp),allocatable :: xyzzyaabl1(:,:)
real(dp) xyzzyaabm1(3,3),xyzzyaabn1(3)
real(dp),allocatable :: xyzzyaabo1(:)
real(dp),allocatable :: xyzzyaabp1(:,:),xyzzyaabq1(:,:),xyzzyaabr1(:,:&
&),xyzzyaabs1(:,:),xyzzyaabt1(:,:),xyzzyaabu1(:,:),xyzzyaabv1(:,:)
complex(dp),allocatable :: xyzzyaabw1(:,:)
integer xyzzyaabx1,xyzzyaaby1,xyzzyaabz1
contains
subroutine cusp_wfdet_real(rvec,jspin,val,fsd,cusp_val,cusp_grad,cusp_&
&lap,cusp_sderivs)
use slaarnabq,only : minimum_image
implicit none
integer,intent(in) :: jspin
real(dp),intent(in) :: rvec(3)
real(dp),intent(out) :: cusp_val(xyzzyaaac1),cusp_grad(3,xyzzyaaac1),c&
&usp_lap(xyzzyaaac1)
real(dp),intent(inout),optional :: cusp_sderivs(6,xyzzyaaac1)
logical,intent(in) :: val,fsd
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2
real(dp) xyzzyaaae2(3),xyzzyaaaf2,xyzzyaaag2,xyzzyaaah2,xyzzyaaai2(4),&
&xyzzyaaaj2(4),xyzzyaaak2(3),xyzzyaaal2,xyzzyaaam2(6),xyzzyaaan2,xyzzy&
&aaao2,xyzzyaaap2,xyzzyaaaq2,xyzzyaaar2(3),xyzzyaaas2,xyzzyaaat2,xyzzy&
&aaau2,xyzzyaaav2,xyzzyaaaw2,xyzzyaaax2(4,1)
if(val)then
!$omp workshare
cusp_val(:)=0.d0
!$omp end workshare
endif
if(fsd)then
!$omp workshare
cusp_grad(:,:)=0.d0
cusp_lap(:)=0.d0
!$omp end workshare
if(present(cusp_sderivs))then
!$omp workshare
cusp_sderivs(:,:)=0.d0
!$omp end workshare
endif
endif
if(pcharge(jspin)/=-1.d0)return
xyzzyaaad2=which_ssingle(jspin,xyzzyaaaa1)
do xyzzyaaaa2=1,naeions_prim
xyzzyaaax2(1:3,1)=rvec-rion(1:3,ae_index(xyzzyaaaa2))
xyzzyaaax2(4,1)=0.d0
if(isperiodic)call minimum_image(4,1,xyzzyaaax2)
xyzzyaaag2=dot_product(xyzzyaaax2(1:3,1),xyzzyaaax2(1:3,1))
if(xyzzyaaag2<xyzzyaaao1(xyzzyaaaa2))then
xyzzyaaae2=xyzzyaaax2(1:3,1)
xyzzyaaaf2=sqrt(xyzzyaaag2)
xyzzyaaai2=(/xyzzyaaaf2,xyzzyaaag2,xyzzyaaaf2*xyzzyaaag2,xyzzyaaag2*xy&
&zzyaaag2/)
if(fsd)then
if(xyzzyaaaf2/=0.d0)then
xyzzyaaah2=1.d0/xyzzyaaaf2
else
xyzzyaaah2=0.d0
endif
xyzzyaaaj2=(/xyzzyaaah2,2.d0,3.d0*xyzzyaaaf2,4.d0*xyzzyaaag2/)
xyzzyaaak2=(/2.d0,6.d0*xyzzyaaaf2,12.d0*xyzzyaaag2/)
if(present(cusp_sderivs))then
xyzzyaaal2=xyzzyaaah2*xyzzyaaah2
xyzzyaaam2(1:3)=xyzzyaaae2**2*xyzzyaaal2
xyzzyaaam2(4)=xyzzyaaae2(1)*xyzzyaaae2(2)*xyzzyaaal2
xyzzyaaam2(5)=xyzzyaaae2(1)*xyzzyaaae2(3)*xyzzyaaal2
xyzzyaaam2(6)=xyzzyaaae2(2)*xyzzyaaae2(3)*xyzzyaaal2
endif
endif
!$omp do private(xyzzyaaab2,xyzzyaaat2,xyzzyaaau2,xyzzyaaaw2,xyzzyaaan&
!$omp &2,xyzzyaaao2,xyzzyaaav2,xyzzyaaap2,xyzzyaaaq2,xyzzyaaar2,xyzzya&
!$omp &aas2)
do xyzzyaaab2=1,xyzzyaaac1
if(.not.xyzzyaaad1(xyzzyaaab2,xyzzyaaad2))cycle
if(xyzzyaaaf2>=xyzzyaaan1(xyzzyaaab2,xyzzyaaaa2))cycle
call eval_spline_reg(nradgrid,radgrid(1:nradgrid,xyzzyaaaa2),xyzzyaabc&
&1(1:nradgrid,xyzzyaaab2,xyzzyaaaa2),xyzzyaabd1(1:nradgrid,xyzzyaaab2,&
&xyzzyaaaa2),xyzzyaaaf2,val,fsd,xyzzyaaat2,xyzzyaaau2,xyzzyaaaw2)
xyzzyaaan2=xyzzyaaar1(xyzzyaaab2,xyzzyaaaa2)
xyzzyaaao2=exp_protect(xyzzyaaal1(0,xyzzyaaab2,xyzzyaaaa2)+ddot(4,xyzz&
&yaaal1(1,xyzzyaaab2,xyzzyaaaa2),1,xyzzyaaai2(1),1))
if(val)cusp_val(xyzzyaaab2)=xyzzyaaan2*(xyzzyaaam1(xyzzyaaab2,xyzzyaaa&
&a2)+xyzzyaaao2-xyzzyaaat2)
if(fsd)then
xyzzyaaav2=xyzzyaaah2*xyzzyaaau2
xyzzyaaap2=ddot(4,xyzzyaaal1(1,xyzzyaaab2,xyzzyaaaa2),1,xyzzyaaaj2(1),&
&1)
xyzzyaaaq2=ddot(3,xyzzyaaal1(2,xyzzyaaab2,xyzzyaaaa2),1,xyzzyaaak2(1),&
&1)
xyzzyaaar2=xyzzyaaap2*xyzzyaaae2
xyzzyaaas2=xyzzyaaaq2+2.d0*xyzzyaaap2
cusp_grad(1:3,xyzzyaaab2)=xyzzyaaan2*(xyzzyaaao2*xyzzyaaar2-xyzzyaaav2&
&*xyzzyaaae2)
cusp_lap(xyzzyaaab2)=xyzzyaaan2*(xyzzyaaao2*(xyzzyaaas2+dot_product(xy&
&zzyaaar2,xyzzyaaar2)) -(xyzzyaaaw2+2.d0*xyzzyaaav2))
if(present(cusp_sderivs))then
do xyzzyaaac2=1,3
cusp_sderivs(xyzzyaaac2,xyzzyaaab2)=xyzzyaaan2*(xyzzyaaao2*(xyzzyaaam2&
&(xyzzyaaac2)*(xyzzyaaaq2-xyzzyaaap2) +xyzzyaaap2+xyzzyaaar2(xyzzyaaac&
&2)**2) -(xyzzyaaav2*(1.d0-xyzzyaaam2(xyzzyaaac2))+xyzzyaaaw2*xyzzyaaa&
&m2(xyzzyaaac2)))
enddo
cusp_sderivs(4,xyzzyaaab2)=xyzzyaaan2*(xyzzyaaao2*(xyzzyaaam2(4)*(xyzz&
&yaaaq2-xyzzyaaap2) +xyzzyaaar2(1)*xyzzyaaar2(2)) +(xyzzyaaav2*xyzzyaa&
&am2(4)-xyzzyaaaw2*xyzzyaaam2(4)))
cusp_sderivs(5,xyzzyaaab2)=xyzzyaaan2*(xyzzyaaao2*(xyzzyaaam2(5)*(xyzz&
&yaaaq2-xyzzyaaap2) +xyzzyaaar2(1)*xyzzyaaar2(3)) +(xyzzyaaav2*xyzzyaa&
&am2(5)-xyzzyaaaw2*xyzzyaaam2(5)))
cusp_sderivs(6,xyzzyaaab2)=xyzzyaaan2*(xyzzyaaao2*(xyzzyaaam2(6)*(xyzz&
&yaaaq2-xyzzyaaap2) +xyzzyaaar2(2)*xyzzyaaar2(3)) +(xyzzyaaav2*xyzzyaa&
&am2(6)-xyzzyaaaw2*xyzzyaaam2(6)))
endif
endif
enddo
!$omp end do
exit
endif
enddo
end subroutine cusp_wfdet_real
subroutine cusp_wfdet_cmplx(rvec,jspin,val,fsd,cusp_val,cusp_grad,cusp&
&_lap,cusp_sderivs)
use slaarnabq,only : min_image_brute_force
implicit none
integer,intent(in) :: jspin
real(dp),intent(in) :: rvec(3)
complex(dp),intent(out) :: cusp_val(xyzzyaaac1),cusp_grad(3,xyzzyaaac1&
&),cusp_lap(xyzzyaaac1)
complex(dp),intent(inout),optional :: cusp_sderivs(6,xyzzyaaac1)
logical,intent(in) :: val,fsd
integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3,xyzzyaaad3
real(dp) xyzzyaaae3(3),xyzzyaaaf3,xyzzyaaag3,xyzzyaaah3,xyzzyaaai3(4),&
&xyzzyaaaj3(4),xyzzyaaak3(3),xyzzyaaal3,xyzzyaaam3(6),xyzzyaaan3,xyzzy&
&aaao3,xyzzyaaap3(3),xyzzyaaaq3,xyzzyaaar3,xyzzyaaas3,xyzzyaaat3,xyzzy&
&aaau3,xyzzyaaav3,xyzzyaaaw3,xyzzyaaax3(3)
complex(dp) xyzzyaaay3
if(val)then
!$omp workshare
cusp_val(:)=czero
!$omp end workshare
endif
if(fsd)then
!$omp workshare
cusp_grad(:,:)=czero
cusp_lap(:)=czero
!$omp end workshare
if(present(cusp_sderivs))then
!$omp workshare
cusp_sderivs(:,:)=czero
!$omp end workshare
endif
endif
if(pcharge(jspin)/=-1.d0)return
xyzzyaaad3=which_ssingle(jspin,xyzzyaaaa1)
do xyzzyaaaa3=1,naeions_prim
xyzzyaaax3=rvec-rion(1:3,ae_index(xyzzyaaaa3))
if(isperiodic)then
call min_image_brute_force(3,xyzzyaaax3,xyzzyaabm1,painv,xyzzyaaae3,xy&
&zzyaabn1)
else
xyzzyaaae3=xyzzyaaax3
endif
xyzzyaaag3=dot_product(xyzzyaaae3,xyzzyaaae3)
if(xyzzyaaag3<xyzzyaaao1(xyzzyaaaa3))then
xyzzyaaax3=xyzzyaaax3-xyzzyaaae3
xyzzyaaaf3=sqrt(xyzzyaaag3)
xyzzyaaai3=(/xyzzyaaaf3,xyzzyaaag3,xyzzyaaaf3*xyzzyaaag3,xyzzyaaag3*xy&
&zzyaaag3/)
if(fsd)then
if(xyzzyaaaf3/=0.d0)then
xyzzyaaah3=1.d0/xyzzyaaaf3
else
xyzzyaaah3=0.d0
endif
xyzzyaaaj3=(/xyzzyaaah3,2.d0,3.d0*xyzzyaaaf3,4.d0*xyzzyaaag3/)
xyzzyaaak3=(/2.d0,6.d0*xyzzyaaaf3,12.d0*xyzzyaaag3/)
if(present(cusp_sderivs))then
xyzzyaaal3=xyzzyaaah3*xyzzyaaah3
xyzzyaaam3(1:3)=xyzzyaaae3**2*xyzzyaaal3
xyzzyaaam3(4)=xyzzyaaae3(1)*xyzzyaaae3(2)*xyzzyaaal3
xyzzyaaam3(5)=xyzzyaaae3(1)*xyzzyaaae3(3)*xyzzyaaal3
xyzzyaaam3(6)=xyzzyaaae3(2)*xyzzyaaae3(3)*xyzzyaaal3
endif
endif
!$omp do private(xyzzyaaab3,xyzzyaaaw3,xyzzyaaay3,xyzzyaaar3,xyzzyaaas&
!$omp &3,xyzzyaaau3,xyzzyaaav3,xyzzyaaat3,xyzzyaaan3,xyzzyaaao3,xyzzya&
!$omp &aap3,xyzzyaaaq3)
do xyzzyaaab3=1,xyzzyaaac1
if(.not.xyzzyaaad1(xyzzyaaab3,xyzzyaaad3))cycle
if(xyzzyaaaf3>=xyzzyaaan1(xyzzyaaab3,xyzzyaaaa3))cycle
if(isperiodic)then
xyzzyaaaw3=dot_product(xyzzyaabl1(1:3,xyzzyaaab3),xyzzyaaax3)
xyzzyaaay3=xyzzyaaas1(xyzzyaaab3,xyzzyaaaa3)*cmplx(cos(xyzzyaaaw3),sin&
&(xyzzyaaaw3),dp)
else
xyzzyaaay3=xyzzyaaas1(xyzzyaaab3,xyzzyaaaa3)
endif
call eval_spline_reg(nradgrid,radgrid(1:nradgrid,xyzzyaaaa3),xyzzyaabc&
&1(1:nradgrid,xyzzyaaab3,xyzzyaaaa3),xyzzyaabd1(1:nradgrid,xyzzyaaab3,&
&xyzzyaaaa3),xyzzyaaaf3,val,fsd,xyzzyaaar3,xyzzyaaas3,xyzzyaaau3)
xyzzyaaav3=exp_protect(xyzzyaaal1(0,xyzzyaaab3,xyzzyaaaa3)+ddot(4,xyzz&
&yaaal1(1,xyzzyaaab3,xyzzyaaaa3),1,xyzzyaaai3(1),1))
if(val)cusp_val(xyzzyaaab3)=xyzzyaaay3*(xyzzyaaam1(xyzzyaaab3,xyzzyaaa&
&a3)+xyzzyaaav3-xyzzyaaar3)
if(fsd)then
xyzzyaaat3=xyzzyaaah3*xyzzyaaas3
xyzzyaaan3=ddot(4,xyzzyaaal1(1,xyzzyaaab3,xyzzyaaaa3),1,xyzzyaaaj3(1),&
&1)
xyzzyaaao3=ddot(3,xyzzyaaal1(2,xyzzyaaab3,xyzzyaaaa3),1,xyzzyaaak3(1),&
&1)
xyzzyaaap3=xyzzyaaan3*xyzzyaaae3
xyzzyaaaq3=xyzzyaaao3+2.d0*xyzzyaaan3
cusp_grad(1:3,xyzzyaaab3)=xyzzyaaay3*(xyzzyaaav3*xyzzyaaap3-xyzzyaaat3&
&*xyzzyaaae3)
cusp_lap(xyzzyaaab3)=xyzzyaaay3*(xyzzyaaav3*(xyzzyaaaq3+dot_product(xy&
&zzyaaap3,xyzzyaaap3))-(xyzzyaaau3+2.d0*xyzzyaaat3))
if(present(cusp_sderivs))then
do xyzzyaaac3=1,3
cusp_sderivs(xyzzyaaac3,xyzzyaaab3)=xyzzyaaay3*(xyzzyaaav3*(xyzzyaaam3&
&(xyzzyaaac3)*(xyzzyaaao3-xyzzyaaan3)+xyzzyaaan3+xyzzyaaap3(xyzzyaaac3&
&)**2)-(xyzzyaaat3*(1.d0-xyzzyaaam3(xyzzyaaac3))+xyzzyaaau3*xyzzyaaam3&
&(xyzzyaaac3)))
enddo
cusp_sderivs(4,xyzzyaaab3)=xyzzyaaay3*(xyzzyaaav3*(xyzzyaaam3(4)*(xyzz&
&yaaao3-xyzzyaaan3)+xyzzyaaap3(1)*xyzzyaaap3(2))+(xyzzyaaat3*xyzzyaaam&
&3(4)-xyzzyaaau3*xyzzyaaam3(4)))
cusp_sderivs(5,xyzzyaaab3)=xyzzyaaay3*(xyzzyaaav3*(xyzzyaaam3(5)*(xyzz&
&yaaao3-xyzzyaaan3)+xyzzyaaap3(1)*xyzzyaaap3(3))+(xyzzyaaat3*xyzzyaaam&
&3(5)-xyzzyaaau3*xyzzyaaam3(5)))
cusp_sderivs(6,xyzzyaaab3)=xyzzyaaay3*(xyzzyaaav3*(xyzzyaaam3(6)*(xyzz&
&yaaao3-xyzzyaaan3)+xyzzyaaap3(2)*xyzzyaaap3(3))+(xyzzyaaat3*xyzzyaaam&
&3(6)-xyzzyaaau3*xyzzyaaam3(6)))
endif
endif
enddo
!$omp end do
exit
endif
enddo
end subroutine cusp_wfdet_cmplx
subroutine initialize_gpcc(spin_dep_gpcc_in,norb_in,orbmask_in,gamma_o&
&nly_in)
use slaarnabq,only : min_image_brute_force
implicit none
integer,intent(in) :: spin_dep_gpcc_in,norb_in
logical,intent(in) :: orbmask_in(norb_in,nspin),gamma_only_in
integer xyzzyaaaa4,xyzzyaaab4,xyzzyaaac4,xyzzyaaad4,xyzzyaaae4,xyzzyaa&
&af4,xyzzyaaag4,xyzzyaaah4
integer,parameter :: xyzzyaaai4=6
real(dp) xyzzyaaaj4,xyzzyaaak4,xyzzyaaal4(3),xyzzyaaam4,xyzzyaaan4
real(dp),parameter :: xyzzyaaao4=1.25d0
real(dp),parameter :: xyzzyaaap4=1.d-8
if(am_master.and.xyzzyaabf1>3)then
call wout('GPCC INITIALIZATION')
call wout('===================')
endif
xyzzyaaaa1=spin_dep_gpcc_in
xyzzyaaab1=no_ssingles(xyzzyaaaa1)
xyzzyaaac1=norb_in
allocate(xyzzyaaad1(xyzzyaaac1,xyzzyaaab1),xyzzyaabo1(xyzzyaaab1),stat&
&=xyzzyaaae4)
call check_alloc(xyzzyaaae4,'INITIALIZE_GPCC','orbmask_gpcc')
do xyzzyaaag4=1,nspin
xyzzyaaad1(:,which_ssingle(xyzzyaaag4,xyzzyaaaa1))=orbmask_in(:,xyzzya&
&aag4)
enddo
do xyzzyaaag4=1,nspin
xyzzyaaaf4=which_ssingle(xyzzyaaag4,xyzzyaaaa1)
xyzzyaabo1(xyzzyaaaf4)=pcharge(xyzzyaaag4)
do xyzzyaaah4=xyzzyaaag4+1,nspin
if(which_ssingle(xyzzyaaah4,xyzzyaaaa1)==xyzzyaaaf4.and.abs(pcharge(xy&
&zzyaaah4)-pcharge(xyzzyaaag4))>xyzzyaaap4)call errstop('INITIALIZE_GP&
&CC','Spins '//trim(i2s(xyzzyaaag4))//' and '//trim(i2s(xyzzyaaah4))//&
&' are supposed to have the same orbitals, but they have different cha&
&rges.')
enddo
enddo
xyzzyaabj1=gamma_only_in
xyzzyaabm1=transpose(pamat)
xyzzyaabn1(1)=sqrt(dot_product(painv(1:3,1),painv(1:3,1)))
xyzzyaabn1(2)=sqrt(dot_product(painv(1:3,2),painv(1:3,2)))
xyzzyaabn1(3)=sqrt(dot_product(painv(1:3,3),painv(1:3,3)))
naeions_prim=0
do xyzzyaaaa4=1,nbasis
if(is_ae(xyzzyaaaa4).and.zion(iontype(xyzzyaaaa4))>0.d0)naeions_prim=n&
&aeions_prim+1
enddo
if(am_master.and.naeions_prim==0)call errstop('INITIALIZE_GPCC','Cusp &
&correction activated, but no all-electron ions are present.')
allocate(ae_index(naeions_prim),stat=xyzzyaaae4)
call check_alloc(xyzzyaaae4,'INITIALIZE_GPCC','ae_index')
xyzzyaaab4=0
do xyzzyaaaa4=1,nbasis
if(is_ae(xyzzyaaaa4).and.zion(iontype(xyzzyaaaa4))>0.d0)then
xyzzyaaab4=xyzzyaaab4+1
ae_index(xyzzyaaab4)=xyzzyaaaa4
endif
enddo
allocate(radgrid(nradgrid,naeions_prim),xyzzyaaap1(naeions_prim),xyzzy&
&aaaq1(naeions_prim),stat=xyzzyaaae4)
call check_alloc(xyzzyaaae4,'INITIALIZE_GPCC','radgrid')
do xyzzyaaab4=1,naeions_prim
xyzzyaaap1(xyzzyaaab4)=1.d0/zion(iontype(ae_index(xyzzyaaab4)))
enddo
do xyzzyaaab4=1,naeions_prim-1
xyzzyaaam4=1.d0/zion(iontype(ae_index(xyzzyaaab4)))
do xyzzyaaac4=xyzzyaaab4+1,naeions_prim
xyzzyaaan4=1.d0/zion(iontype(ae_index(xyzzyaaac4)))
if(isperiodic)then
call min_image_brute_force(3,rion(1:3,ae_index(xyzzyaaac4))-rion(1:3,a&
&e_index(xyzzyaaab4)),xyzzyaabm1,painv,xyzzyaaal4,xyzzyaabn1)
else
xyzzyaaal4=rion(1:3,ae_index(xyzzyaaac4))-rion(1:3,ae_index(xyzzyaaab4&
&))
endif
xyzzyaaak4=dot_product(xyzzyaaal4,xyzzyaaal4)
if(xyzzyaaak4<(xyzzyaaam4+xyzzyaaan4)**2)then
xyzzyaaak4=sqrt(xyzzyaaak4)
xyzzyaaap1(xyzzyaaab4)=min(xyzzyaaap1(xyzzyaaab4),xyzzyaaam4*xyzzyaaak&
&4/(xyzzyaaam4+xyzzyaaan4))
xyzzyaaap1(xyzzyaaac4)=min(xyzzyaaap1(xyzzyaaac4),xyzzyaaan4*xyzzyaaak&
&4/(xyzzyaaam4+xyzzyaaan4))
endif
enddo
enddo
do xyzzyaaab4=1,naeions_prim
xyzzyaaaj4=xyzzyaaao4*xyzzyaaap1(xyzzyaaab4)/dble(nradgrid-1)
xyzzyaaaq1(xyzzyaaab4)=floor(xyzzyaaap1(xyzzyaaab4)/xyzzyaaaj4)+1
do xyzzyaaad4=1,nradgrid
radgrid(xyzzyaaad4,xyzzyaaab4)=dble(xyzzyaaad4-1)*xyzzyaaaj4
enddo
if(am_master.and.xyzzyaabf1>3)call wout('rcmax for ion '//trim(i2s(ae_&
&index(xyzzyaaab4)))//' : ',xyzzyaaap1(xyzzyaaab4),rfmt='(es24.16)')
enddo
if(am_master.and.xyzzyaabf1>3)call wout('Rule for spherical integratio&
&n: '//trim(i2s(xyzzyaaai4)))
call xyzzyaaca1(xyzzyaaai4)
if(am_master.and.xyzzyaabf1>3)call wout('Spherical grid set up.')
if(xyzzyaabj1)then
allocate(xyzzyaabc1(nradgrid,xyzzyaaac1,naeions_prim),xyzzyaaar1(xyzzy&
&aaac1,naeions_prim),stat=xyzzyaaae4)
else
allocate(xyzzyaabc1(nradgrid,xyzzyaaac1,naeions_prim),xyzzyaaas1(xyzzy&
&aaac1,naeions_prim),stat=xyzzyaaae4)
endif
call check_alloc(xyzzyaaae4,'SETUP_GPCC','orbbar_radgrid')
xyzzyaabc1=0.d0
if(xyzzyaabj1)then
xyzzyaaar1=0.d0
else
xyzzyaaas1=czero
endif
if(am_master.and.xyzzyaabf1>3)then
call wout('Finished initialization.')
call wout()
endif
end subroutine initialize_gpcc
subroutine xyzzyaaca1(irule)
use slaarnaag,only : one_over_root_two,twentyoneth,fifteenth,thirtieth&
&,twelfth,sixth,eleventh,third
implicit none
integer,intent(in) :: irule
integer xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5,xyzzyaaad5,xyzzyaaae5,xyzzyaa&
&af5
real(dp) xyzzyaaag5,xyzzyaaah5,xyzzyaaai5,xyzzyaaaj5,xyzzyaaak5
xyzzyaaag5=sqrt(third)
xyzzyaaah5=sqrt(eleventh)
xyzzyaaai5=3.d0*xyzzyaaah5
xyzzyaaaj5=sqrt((5.d0-sqrt(5.d0))*0.1d0)
xyzzyaaak5=sqrt((5.d0+sqrt(5.d0))*0.1d0)
select case (irule)
case(1)
nsphgrid=1
allocate(sphgrid(3,nsphgrid),xyzzyaabe1(nsphgrid),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_SPHERICAL_GRID','')
sphgrid(1:3,1)=(/0.d0,0.d0,1.d0/)
xyzzyaabe1(1)=1.d0
case(2)
nsphgrid=4
allocate(sphgrid(3,nsphgrid),xyzzyaabe1(nsphgrid),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_SPHERICAL_GRID','')
sphgrid(1:3,1)=(/xyzzyaaag5,xyzzyaaag5,xyzzyaaag5/)
sphgrid(1:3,2)=(/xyzzyaaag5,-xyzzyaaag5,-xyzzyaaag5/)
sphgrid(1:3,3)=(/-xyzzyaaag5,xyzzyaaag5,-xyzzyaaag5/)
sphgrid(1:3,4)=(/-xyzzyaaag5,-xyzzyaaag5,xyzzyaaag5/)
xyzzyaabe1(1:4)=0.25d0
case(3)
nsphgrid=6
allocate(sphgrid(3,nsphgrid),xyzzyaabe1(nsphgrid),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_SPHERICAL_GRID','')
sphgrid(1:3,1)=(/1.d0,0.d0,0.d0/)
sphgrid(1:3,3)=(/0.d0,1.d0,0.d0/)
sphgrid(1:3,5)=(/0.d0,0.d0,1.d0/)
do xyzzyaaaa5=1,5,2
sphgrid(1:3,xyzzyaaaa5+1)=-sphgrid(1:3,xyzzyaaaa5)
enddo
xyzzyaabe1(1:6)=sixth
case(4)
nsphgrid=12
allocate(sphgrid(3,nsphgrid),xyzzyaabe1(nsphgrid),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_SPHERICAL_GRID','')
sphgrid(1:3,1)=(/0.d0,xyzzyaaaj5,xyzzyaaak5/)
sphgrid(1:3,3)=(/0.d0,-xyzzyaaaj5,xyzzyaaak5/)
sphgrid(1:3,5)=(/-xyzzyaaak5,0.d0,xyzzyaaaj5/)
sphgrid(1:3,7)=(/xyzzyaaak5,0.d0,xyzzyaaaj5/)
sphgrid(1:3,9)=(/-xyzzyaaaj5,xyzzyaaak5,0.d0/)
sphgrid(1:3,11)=(/xyzzyaaaj5,xyzzyaaak5,0.d0/)
do xyzzyaaaa5=1,11,2
sphgrid(1:3,xyzzyaaaa5+1)=-sphgrid(1:3,xyzzyaaaa5)
enddo
xyzzyaabe1(1:12)=twelfth
case(5)
nsphgrid=18
allocate(sphgrid(3,nsphgrid),xyzzyaabe1(nsphgrid),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_SPHERICAL_GRID','')
sphgrid(1:3,1)=(/1.d0,0.d0,0.d0/)
sphgrid(1:3,3)=(/0.d0,1.d0,0.d0/)
sphgrid(1:3,5)=(/0.d0,0.d0,1.d0/)
do xyzzyaaaa5=1,5,2
sphgrid(1:3,xyzzyaaaa5+1)=-sphgrid(1:3,xyzzyaaaa5)
enddo
xyzzyaabe1(1:6)=thirtieth
do xyzzyaaaa5=7,15,4
if(xyzzyaaaa5==7)then
xyzzyaaac5=1
xyzzyaaad5=2
xyzzyaaae5=3
elseif(xyzzyaaaa5==11)then
xyzzyaaac5=2
xyzzyaaad5=3
xyzzyaaae5=1
else
xyzzyaaac5=3
xyzzyaaad5=1
xyzzyaaae5=2
endif
sphgrid(xyzzyaaac5,xyzzyaaaa5)=one_over_root_two
sphgrid(xyzzyaaad5,xyzzyaaaa5)=one_over_root_two
sphgrid(xyzzyaaae5,xyzzyaaaa5)=0.d0
sphgrid(xyzzyaaac5,xyzzyaaaa5+1)=one_over_root_two
sphgrid(xyzzyaaad5,xyzzyaaaa5+1)=-one_over_root_two
sphgrid(xyzzyaaae5,xyzzyaaaa5+1)=0.d0
sphgrid(xyzzyaaac5,xyzzyaaaa5+2)=-one_over_root_two
sphgrid(xyzzyaaad5,xyzzyaaaa5+2)=one_over_root_two
sphgrid(xyzzyaaae5,xyzzyaaaa5+2)=0.d0
sphgrid(xyzzyaaac5,xyzzyaaaa5+3)=-one_over_root_two
sphgrid(xyzzyaaad5,xyzzyaaaa5+3)=-one_over_root_two
sphgrid(xyzzyaaae5,xyzzyaaaa5+3)=0.d0
enddo
xyzzyaabe1(7:18)=fifteenth
case(6)
nsphgrid=26
allocate(sphgrid(3,nsphgrid),xyzzyaabe1(nsphgrid),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_SPHERICAL_GRID','')
sphgrid(1:3,1)=(/1.d0,0.d0,0.d0/)
sphgrid(1:3,3)=(/0.d0,1.d0,0.d0/)
sphgrid(1:3,5)=(/0.d0,0.d0,1.d0/)
do xyzzyaaaa5=1,5,2
sphgrid(1:3,xyzzyaaaa5+1)=-sphgrid(1:3,xyzzyaaaa5)
enddo
xyzzyaabe1(1:6)=twentyoneth
do xyzzyaaaa5=7,15,4
if(xyzzyaaaa5==7)then
xyzzyaaac5=1
xyzzyaaad5=2
xyzzyaaae5=3
elseif(xyzzyaaaa5==11)then
xyzzyaaac5=2
xyzzyaaad5=3
xyzzyaaae5=1
else
xyzzyaaac5=3
xyzzyaaad5=1
xyzzyaaae5=2
endif
sphgrid(xyzzyaaac5,xyzzyaaaa5)=one_over_root_two
sphgrid(xyzzyaaad5,xyzzyaaaa5)=one_over_root_two
sphgrid(xyzzyaaae5,xyzzyaaaa5)=0.d0
sphgrid(xyzzyaaac5,xyzzyaaaa5+1)=one_over_root_two
sphgrid(xyzzyaaad5,xyzzyaaaa5+1)=-one_over_root_two
sphgrid(xyzzyaaae5,xyzzyaaaa5+1)=0.d0
sphgrid(xyzzyaaac5,xyzzyaaaa5+2)=-one_over_root_two
sphgrid(xyzzyaaad5,xyzzyaaaa5+2)=one_over_root_two
sphgrid(xyzzyaaae5,xyzzyaaaa5+2)=0.d0
sphgrid(xyzzyaaac5,xyzzyaaaa5+3)=-one_over_root_two
sphgrid(xyzzyaaad5,xyzzyaaaa5+3)=-one_over_root_two
sphgrid(xyzzyaaae5,xyzzyaaaa5+3)=0.d0
enddo
xyzzyaabe1(7:18)=4.d0/105.d0
do xyzzyaaaa5=0,7
do xyzzyaaab5=1,3
if(btest(xyzzyaaaa5,xyzzyaaab5-1))then
sphgrid(xyzzyaaab5,19+xyzzyaaaa5)=-xyzzyaaag5
else
sphgrid(xyzzyaaab5,19+xyzzyaaaa5)=xyzzyaaag5
endif
enddo
enddo
xyzzyaabe1(19:26)=27.d0/840.d0
case(7)
nsphgrid=50
allocate(sphgrid(3,nsphgrid),xyzzyaabe1(nsphgrid),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_SPHERICAL_GRID','')
sphgrid(1:3,1)=(/1.d0,0.d0,0.d0/)
sphgrid(1:3,3)=(/0.d0,1.d0,0.d0/)
sphgrid(1:3,5)=(/0.d0,0.d0,1.d0/)
do xyzzyaaaa5=1,5,2
sphgrid(1:3,xyzzyaaaa5+1)=-sphgrid(1:3,xyzzyaaaa5)
enddo
xyzzyaabe1(1:6)=4.d0/315.d0
do xyzzyaaaa5=7,15,4
if(xyzzyaaaa5==7)then
xyzzyaaac5=1
xyzzyaaad5=2
xyzzyaaae5=3
elseif(xyzzyaaaa5==11)then
xyzzyaaac5=2
xyzzyaaad5=3
xyzzyaaae5=1
else
xyzzyaaac5=3
xyzzyaaad5=1
xyzzyaaae5=2
endif
sphgrid(xyzzyaaac5,xyzzyaaaa5)=one_over_root_two
sphgrid(xyzzyaaad5,xyzzyaaaa5)=one_over_root_two
sphgrid(xyzzyaaae5,xyzzyaaaa5)=0.d0
sphgrid(xyzzyaaac5,xyzzyaaaa5+1)=one_over_root_two
sphgrid(xyzzyaaad5,xyzzyaaaa5+1)=-one_over_root_two
sphgrid(xyzzyaaae5,xyzzyaaaa5+1)=0.d0
sphgrid(xyzzyaaac5,xyzzyaaaa5+2)=-one_over_root_two
sphgrid(xyzzyaaad5,xyzzyaaaa5+2)=one_over_root_two
sphgrid(xyzzyaaae5,xyzzyaaaa5+2)=0.d0
sphgrid(xyzzyaaac5,xyzzyaaaa5+3)=-one_over_root_two
sphgrid(xyzzyaaad5,xyzzyaaaa5+3)=-one_over_root_two
sphgrid(xyzzyaaae5,xyzzyaaaa5+3)=0.d0
enddo
xyzzyaabe1(7:18)=64.d0/2835.d0
sphgrid(1:3,19:26)=xyzzyaaag5
do xyzzyaaaa5=0,7
do xyzzyaaab5=1,3
if(btest(xyzzyaaaa5,xyzzyaaab5-1))sphgrid(xyzzyaaab5,19+xyzzyaaaa5)=-s&
&phgrid(xyzzyaaab5,19+xyzzyaaaa5)
enddo
enddo
xyzzyaabe1(19:26)=27.d0/1280.d0
sphgrid(1,27:34)=xyzzyaaah5
sphgrid(2,27:34)=xyzzyaaah5
sphgrid(3,27:34)=xyzzyaaai5
sphgrid(1,35:42)=xyzzyaaah5
sphgrid(2,35:42)=xyzzyaaai5
sphgrid(3,35:42)=xyzzyaaah5
sphgrid(1,43:50)=xyzzyaaai5
sphgrid(2,43:50)=xyzzyaaah5
sphgrid(3,43:50)=xyzzyaaah5
do xyzzyaaaa5=0,7
do xyzzyaaab5=1,3
if(btest(xyzzyaaaa5,xyzzyaaab5-1))then
sphgrid(xyzzyaaab5,27+xyzzyaaaa5)=-sphgrid(xyzzyaaab5,27+xyzzyaaaa5)
sphgrid(xyzzyaaab5,35+xyzzyaaaa5)=-sphgrid(xyzzyaaab5,35+xyzzyaaaa5)
sphgrid(xyzzyaaab5,43+xyzzyaaaa5)=-sphgrid(xyzzyaaab5,43+xyzzyaaaa5)
endif
enddo
enddo
xyzzyaabe1(27:50)=14641.d0/725760.d0
case default
call errstop('SETUP_SPHERICAL_GRID','Unsupported grid code: ' //trim(i&
&2s(irule))//'.')
end select
end subroutine xyzzyaaca1
subroutine spherical_av_real(spin,jion,igrid,orb_sphgrid_r)
implicit none
integer,intent(in) :: spin,igrid,jion
real(dp),intent(in) :: orb_sphgrid_r(nsphgrid,xyzzyaaac1)
integer xyzzyaaaa6
if(.not.xyzzyaabj1)call errstop('SPHERICAL_AV_REAL','Erroneous call.')
if(igrid==1)then
do xyzzyaaaa6=1,xyzzyaaac1
if(xyzzyaaad1(xyzzyaaaa6,spin))then
xyzzyaabc1(1,xyzzyaaaa6,jion)=orb_sphgrid_r(1,xyzzyaaaa6)
if(xyzzyaabc1(1,xyzzyaaaa6,jion)>=0.d0)then
xyzzyaaar1(xyzzyaaaa6,jion)=1.d0
else
xyzzyaaar1(xyzzyaaaa6,jion)=-1.d0
xyzzyaabc1(1,xyzzyaaaa6,jion)=-xyzzyaabc1(1,xyzzyaaaa6,jion)
endif
endif
enddo
else
do xyzzyaaaa6=1,xyzzyaaac1
if(xyzzyaaad1(xyzzyaaaa6,spin))then
if(xyzzyaaar1(xyzzyaaaa6,jion)>0.d0)then
xyzzyaabc1(igrid,xyzzyaaaa6,jion)=ddot(nsphgrid,xyzzyaabe1(1),1,orb_sp&
&hgrid_r(1,xyzzyaaaa6),1)
else
xyzzyaabc1(igrid,xyzzyaaaa6,jion)=-ddot(nsphgrid,xyzzyaabe1(1),1,orb_s&
&phgrid_r(1,xyzzyaaaa6),1)
endif
endif
enddo
endif
end subroutine spherical_av_real
subroutine spherical_av_cmplx(spin,jion,igrid,orb_sphgrid_c)
implicit none
integer,intent(in) :: spin,igrid,jion
complex(dp),intent(in) :: orb_sphgrid_c(nsphgrid,xyzzyaaac1)
integer xyzzyaaaa7,xyzzyaaab7
real(dp) xyzzyaaac7,xyzzyaaad7,xyzzyaaae7
if(xyzzyaabj1)call errstop('SPHERICAL_AV_CMPLX','Erroneous call.')
if(igrid==1)then
do xyzzyaaaa7=1,xyzzyaaac1
if(xyzzyaaad1(xyzzyaaaa7,spin))then
xyzzyaabc1(1,xyzzyaaaa7,jion)=abs(orb_sphgrid_c(1,xyzzyaaaa7))
if(xyzzyaabc1(1,xyzzyaaaa7,jion)>0.d0)then
xyzzyaaas1(xyzzyaaaa7,jion)=orb_sphgrid_c(1,xyzzyaaaa7)/xyzzyaabc1(1,x&
&yzzyaaaa7,jion)
else
xyzzyaaas1(xyzzyaaaa7,jion)=c_one
endif
endif
enddo
else
do xyzzyaaaa7=1,xyzzyaaac1
if(xyzzyaaad1(xyzzyaaaa7,spin))then
xyzzyaaad7=dble(xyzzyaaas1(xyzzyaaaa7,jion))
xyzzyaaae7=aimag(xyzzyaaas1(xyzzyaaaa7,jion))
xyzzyaaac7=0.d0
do xyzzyaaab7=1,nsphgrid
xyzzyaaac7=xyzzyaaac7+xyzzyaabe1(xyzzyaaab7)*(xyzzyaaad7*dble(orb_sphg&
&rid_c(xyzzyaaab7,xyzzyaaaa7))+xyzzyaaae7*aimag(orb_sphgrid_c(xyzzyaaa&
&b7,xyzzyaaaa7)))
enddo
xyzzyaabc1(igrid,xyzzyaaaa7,jion)=xyzzyaaac7
endif
enddo
endif
end subroutine spherical_av_cmplx
subroutine setup_gpcc(kvec_in)
implicit none
real(dp),intent(in),optional :: kvec_in(3,xyzzyaaac1)
integer xyzzyaaaa8,xyzzyaaab8,xyzzyaaac8,xyzzyaaad8,xyzzyaaae8,xyzzyaa&
&af8
integer,parameter :: xyzzyaaag8=10
real(dp) xyzzyaaah8,xyzzyaaai8,xyzzyaaaj8,xyzzyaaak8,xyzzyaaal8,xyzzya&
&aam8,xyzzyaaan8
real(dp),parameter :: xyzzyaaao8=2.d0
real(dp),parameter :: xyzzyaaap8=1.d-3
real(dp),parameter :: xyzzyaaaq8=0.2d0
real(dp),parameter :: xyzzyaaar8=0.1d0
real(dp),allocatable :: xyzzyaaas8(:)
logical xyzzyaaat8
logical,save :: xyzzyaaau8=.true.
real(dp) objfn_outer
external objfn_outer
if(xyzzyaabf1>1.and.am_master.and.xyzzyaaau8)then
call wout('GPCC setup')
call wout('==========')
endif
if(xyzzyaaau8)then
allocate(xyzzyaaan1(xyzzyaaac1,naeions_prim),xyzzyaaao1(naeions_prim),&
&xyzzyaaal1(0:4,xyzzyaaac1,naeions_prim),xyzzyaabd1(nradgrid,xyzzyaaac&
&1,naeions_prim),xyzzyaaam1(xyzzyaaac1,naeions_prim),stat=xyzzyaaad8)
call check_alloc(xyzzyaaad8,'SETUP_GPCC','rc_arr')
xyzzyaaan1=0.d0
xyzzyaaao1=0.d0
xyzzyaaal1=0.d0
xyzzyaabd1=0.d0
xyzzyaaam1=0.d0
endif
allocate(xyzzyaaas8(xyzzyaabh1),stat=xyzzyaaad8)
call check_alloc(xyzzyaaad8,'SETUP_GPCC','rnode')
do xyzzyaaab8=1,xyzzyaaab1
if(xyzzyaabo1(xyzzyaaab8)/=-1.d0)cycle
if(am_master.and.xyzzyaabf1>3.and.xyzzyaaau8.and.xyzzyaaab1>1)call wou&
&t('SPIN TYPE '//trim(i2s(xyzzyaaab8)))
do xyzzyaaaa8=1,naeions_prim
if(am_master.and.xyzzyaabf1>3.and.xyzzyaaau8)then
call wout(' ION '//trim(i2s(ae_index(xyzzyaaaa8))))
call wout('  rcmax : ',xyzzyaaap1(xyzzyaaaa8))
endif
do xyzzyaaac8=1,xyzzyaaac1
if(.not.xyzzyaaad1(xyzzyaaac8,xyzzyaaab8))cycle
if(am_master.and.xyzzyaabf1>3.and.xyzzyaaau8)call wout('  ORBITAL ' //&
&trim(i2s(xyzzyaaac8)))
if(all(xyzzyaabc1(1:nradgrid,xyzzyaaac8,xyzzyaaaa8)==0.d0))then
xyzzyaaan1(xyzzyaaac8,xyzzyaaaa8)=0.d0
if(am_master.and.xyzzyaabf1>3.and.xyzzyaaau8)call wout('   Orbital is &
&zero at nucleus.  Ignoring.')
cycle
endif
call setup_spline_reg(nradgrid,radgrid(1:nradgrid,xyzzyaaaa8),xyzzyaab&
&c1(1:nradgrid,xyzzyaaac8,xyzzyaaaa8),xyzzyaabd1(1:nradgrid,xyzzyaaac8&
&,xyzzyaaaa8),0.d0)
if(am_master.and.xyzzyaabf1>3.and.xyzzyaaau8)call wout('   Spline gene&
&rated.')
if(minval(xyzzyaabc1(1:xyzzyaaaq1(xyzzyaaaa8),xyzzyaaac8,xyzzyaaaa8))<&
&0.d0)then
xyzzyaaam1(xyzzyaaac8,xyzzyaaaa8)=minval(xyzzyaabc1(1:xyzzyaaaq1(xyzzy&
&aaaa8),xyzzyaaac8,xyzzyaaaa8))*xyzzyaaao8
else
xyzzyaaam1(xyzzyaaac8,xyzzyaaaa8)=0.d0
endif
if(am_master.and.xyzzyaabf1>3.and.xyzzyaaau8)then
if(xyzzyaabj1)then
if(xyzzyaaar1(xyzzyaaac8,xyzzyaaaa8)>=0.d0)then
call wout('   Orbital is positive at nucleus.')
else
call wout('   Orbital is negative at nucleus.')
endif
else
call wout('   Phase at nucleus: ',xyzzyaaas1(xyzzyaaac8,xyzzyaaaa8))
endif
call wout('   C (inc. sign)   : ',xyzzyaaam1(xyzzyaaac8,xyzzyaaaa8))
endif
call xyzzyaach1(xyzzyaaac8,xyzzyaaaa8,xyzzyaaau8,xyzzyaaaf8,xyzzyaaas8&
&)
call xyzzyaacg1(xyzzyaaac8,xyzzyaaaa8,xyzzyaaaf8,xyzzyaaas8,xyzzyaaai8&
&)
if(am_master.and.xyzzyaabf1>3.and.xyzzyaaau8)call wout('   Initial rc &
&  : ',xyzzyaaai8)
xyzzyaaba1=xyzzyaaam1(xyzzyaaac8,xyzzyaaaa8)
xyzzyaabb1=xyzzyaabc1(1,xyzzyaaac8,xyzzyaaaa8)
xyzzyaaav1=zion(iontype(ae_index(xyzzyaaaa8)))
xyzzyaaat1=xyzzyaaac8
xyzzyaaau1=xyzzyaaaa8
if(xyzzyaaaf8>0)then
xyzzyaabi1=any(xyzzyaaas8(1:xyzzyaaaf8)<xyzzyaaar8*xyzzyaaap1(xyzzyaaa&
&a8))
else
xyzzyaabi1=.false.
endif
if(am_master.and.xyzzyaabf1>3.and.xyzzyaaau8.and.xyzzyaabi1) call wout&
&('   Nodes lie close to r=0.')
if(.not.xyzzyaabi1)then
xyzzyaaaj8=(1.d0-xyzzyaaaq8)*xyzzyaaai8
xyzzyaaak8=min((1.d0+xyzzyaaaq8)*xyzzyaaai8,xyzzyaaap1(xyzzyaaaa8))
if(xyzzyaaaf8>0)xyzzyaaak8=min(xyzzyaaak8,xyzzyaabk1*xyzzyaaas8(1))
call bracket(xyzzyaaaj8,xyzzyaaak8,objfn_outer,xyzzyaaag8,xyzzyaaal8,x&
&yzzyaaam8,xyzzyaaan8,xyzzyaaat8)
if(xyzzyaaat8)then
call brentmin(xyzzyaaal8,xyzzyaaam8,xyzzyaaan8,objfn_outer,xyzzyaaai8,&
&xyzzyaaah8,xyzzyaaap8)
if(am_master.and.xyzzyaabf1>3.and.xyzzyaaau8)call wout('   Bracketing &
&of min of obj fn with respect to rc successful.')
elseif(objfn_outer(xyzzyaaaj8)<=objfn_outer(xyzzyaaak8))then
xyzzyaaai8=xyzzyaaaj8
else
xyzzyaaai8=xyzzyaaak8
endif
if(.not.xyzzyaaat8.and.am_master.and.xyzzyaabf1>3.and.xyzzyaaau8)call &
&wout('   Bracketing of min of obj fn with respect to rc unsuccessful.&
&')
if(am_master.and.xyzzyaabf1>3.and.xyzzyaaau8)call wout('   Optimized r&
&c : ',xyzzyaaai8)
endif
xyzzyaaan1(xyzzyaaac8,xyzzyaaaa8)=xyzzyaaai8
call determine_poly_outer(xyzzyaaai8,xyzzyaaau8,xyzzyaaah8,xyzzyaaal1(&
&0:4,xyzzyaaac8,xyzzyaaaa8))
if(am_master.and.xyzzyaabf1>3.and.xyzzyaaau8)then
call wout('   Polynomial coefficients evaluated.')
do xyzzyaaae8=0,4
call wout('   alpha_'//trim(i2s(xyzzyaaae8))//' : ',xyzzyaaal1(xyzzyaa&
&ae8,xyzzyaaac8,xyzzyaaaa8))
enddo
endif
enddo
xyzzyaaao1(xyzzyaaaa8)=maxval(xyzzyaaan1(1:xyzzyaaac1,xyzzyaaaa8))**2
enddo
enddo
deallocate(xyzzyaaas8)
if(.not.xyzzyaabj1.and.xyzzyaaau8)then
if(.not.present(kvec_in))call errstop('SETUP_GPCC','No k vectors.')
allocate(xyzzyaabl1(3,xyzzyaaac1),stat=xyzzyaaad8)
call check_alloc(xyzzyaaad8,'SETUP_GPCC','kvec_gpcc')
call dcopy(3*xyzzyaaac1,kvec_in(1,1),1,xyzzyaabl1(1,1),1)
endif
if(am_master.and.xyzzyaabf1>1.and.xyzzyaaau8)then
if(xyzzyaabf1>3)call wout()
call wordwrap('Cusp correction successfully constructed for each orbit&
&al at each ion.')
call wout()
endif
xyzzyaaau8=.false.
end subroutine setup_gpcc
subroutine determine_poly_outer(rc,first_setup,objfn_val,alpha)
implicit none
real(dp),intent(in) :: rc
real(dp),intent(out) :: objfn_val
real(dp),intent(out),optional :: alpha(0:4)
logical,intent(in) :: first_setup
real(dp) xyzzyaaaa9,xyzzyaaab9,xyzzyaaac9,xyzzyaaad9,xyzzyaaae9,xyzzya&
&aaf9,xyzzyaaag9
real(dp),parameter :: xyzzyaaah9=1.d-8
real(dp) objfn_inner
external objfn_inner
call eval_spline_reg(nradgrid,radgrid(1:nradgrid,xyzzyaaau1),xyzzyaabc&
&1(1:nradgrid,xyzzyaaat1,xyzzyaaau1),xyzzyaabd1(1:nradgrid,xyzzyaaat1,&
&xyzzyaaau1),rc,.true.,.true.,xyzzyaaab9,xyzzyaaac9,xyzzyaaad9)
if(present(alpha).and.am_master.and.xyzzyaabf1>4.and.first_setup)then
call wout('   Uncorr s-part of orb at r=0               : ',xyzzyaabb1&
&)
call wout('   Uncorr s-part of orb at r=rc              : ',xyzzyaaab9&
&)
call wout('   Deriv of uncorr s-part of orb at r=rc     : ',xyzzyaaac9&
&)
call wout('   2nd deriv of uncorr s-part of orb at r=rc : ',xyzzyaaad9&
&)
endif
xyzzyaaaz1=rc
if(xyzzyaaab9==xyzzyaaba1)call errstop('DETERMINE_POLY_OUTER','Orbital&
& value at rc is equal to C.')
xyzzyaaaw1=log(xyzzyaaab9-xyzzyaaba1)
xyzzyaaax1=xyzzyaaac9/(xyzzyaaab9-xyzzyaaba1)
xyzzyaaay1=xyzzyaaad9/(xyzzyaaab9-xyzzyaaba1)
if(present(alpha).and.am_master.and.xyzzyaabf1>4.and.first_setup)then
call wout('   X1 : ',xyzzyaaaw1)
call wout('   X2 : ',xyzzyaaax1)
call wout('   X3 : ',xyzzyaaay1)
endif
if(.not.xyzzyaabi1)then
xyzzyaaae9=log(xyzzyaabb1-xyzzyaaba1)
xyzzyaaaf9=log(2.d0*xyzzyaabb1-xyzzyaaba1)
if(xyzzyaaae9==xyzzyaaaf9)xyzzyaaaf9=xyzzyaaae9+xyzzyaaah9
call gsbracket(xyzzyaaae9,xyzzyaaaf9,xyzzyaaag9,objfn_inner)
if(present(alpha).and.am_master.and.xyzzyaabf1>4.and.first_setup)then
call wout('   p(0) (A); obj fn : ',(/xyzzyaaae9,objfn_inner(xyzzyaaae9&
&)/),rfmt='(es20.12)')
call wout('   p(0) (B); obj fn : ',(/xyzzyaaaf9,objfn_inner(xyzzyaaaf9&
&)/),rfmt='(es20.12)')
call wout('   p(0) (C); obj fn : ',(/xyzzyaaag9,objfn_inner(xyzzyaaag9&
&)/),rfmt='(es20.12)')
call wout('    (These values should bracket the minimum of the objecti&
&ve fn.)')
endif
call brentmin(xyzzyaaae9,xyzzyaaaf9,xyzzyaaag9,objfn_inner,xyzzyaaaa9,&
&objfn_val,xyzzyaaah9)
else
xyzzyaaaa9=log(xyzzyaabb1-xyzzyaaba1)
endif
if(present(alpha))then
if(am_master.and.xyzzyaabf1>4.and.first_setup)then
if(xyzzyaabi1)then
call wout('   Orbital value at r=0 chosen.')
call wout('   p(0): ',xyzzyaaaa9,rfmt='(es20.12)')
else
call wout('   Brent minimization finished.')
call wout('   p(0); obj fn: ',(/xyzzyaaaa9,objfn_val/),rfmt='(es20.12)&
&')
endif
endif
call determine_poly_inner(xyzzyaaaa9,first_setup,objfn_val,alpha)
endif
end subroutine determine_poly_outer
subroutine determine_poly_inner(p0,first_setup,objfn_val,alpha_out)
implicit none
real(dp),intent(in) :: p0
real(dp),intent(out) :: objfn_val
real(dp),intent(out),optional :: alpha_out(0:4)
logical,intent(in) :: first_setup
real(dp) xyzzyaaaa10,xyzzyaaab10,xyzzyaaac10(0:4)
xyzzyaaaa10=-xyzzyaaav1*(xyzzyaaba1*exp_protect(-p0)+1.d0)
xyzzyaaab10=p0
if(present(alpha_out).and.am_master.and.xyzzyaabf1>4.and.first_setup)t&
&hen
call wout('   X4 : ',xyzzyaaaa10)
call wout('   X5 : ',xyzzyaaab10)
endif
call xyzzyaacb1(xyzzyaaaw1,xyzzyaaax1,xyzzyaaay1,xyzzyaaaa10,xyzzyaaab&
&10,xyzzyaaaz1,xyzzyaaac10)
if(xyzzyaabi1)then
objfn_val=0.d0
else
objfn_val=xyzzyaace1(xyzzyaaav1,xyzzyaaac10,xyzzyaaba1,xyzzyaaaz1)
endif
if(present(alpha_out))alpha_out=xyzzyaaac10
end subroutine determine_poly_inner
subroutine xyzzyaacb1(x1,x2,x3,x4,x5,rc,alpha)
implicit none
real(dp),intent(in) :: x1,x2,x3,x4,x5,rc
real(dp),intent(out) :: alpha(0:4)
real(dp) xyzzyaaaa11,xyzzyaaab11,xyzzyaaac11
if(rc>0.d0)then
xyzzyaaaa11=1.d0/rc
else
xyzzyaaaa11=0.d0
endif
xyzzyaaab11=x3-x2*x2
xyzzyaaac11=x1-x5
alpha(0)=x5
alpha(1)=x4
alpha(2)=0.5d0*xyzzyaaab11+3.d0*xyzzyaaaa11*(-x2-x4+2.d0*xyzzyaaaa11*x&
&yzzyaaac11)
alpha(3)=xyzzyaaaa11*(-xyzzyaaab11+xyzzyaaaa11*(5.d0*x2+3.d0*x4-8.d0*x&
&yzzyaaaa11*xyzzyaaac11))
alpha(4)=xyzzyaaaa11*xyzzyaaaa11*(0.5d0*xyzzyaaab11+xyzzyaaaa11*(-x4-2&
&.d0*x2+3.d0*xyzzyaaaa11*xyzzyaaac11))
end subroutine xyzzyaacb1
real(dp) function xyzzyaacc1(z,beta0,r)
implicit none
real(dp),intent(in) :: z,beta0,r
xyzzyaacc1=z*z*(beta0+r*r*(xyzzyaaae1+r*(xyzzyaaaf1+r*(xyzzyaaag1+r*(x&
&yzzyaaah1+r*(xyzzyaaai1+r*(xyzzyaaaj1+r*xyzzyaaak1)))))))
end function xyzzyaacc1
real(dp) function xyzzyaacd1(z,alpha,c,r)
implicit none
real(dp),intent(in) :: z,alpha(0:4),r,c
integer xyzzyaaaa13
real(dp) xyzzyaaab13,xyzzyaaac13,xyzzyaaad13,xyzzyaaae13,xyzzyaaaf13,x&
&yzzyaaag13,xyzzyaaah13
if(r>0.d0)then
xyzzyaaab13=alpha(0)
xyzzyaaac13=0.d0
xyzzyaaad13=0.d0
xyzzyaaaf13=1.d0
xyzzyaaag13=0.d0
do xyzzyaaaa13=1,4
xyzzyaaad13=xyzzyaaad13+alpha(xyzzyaaaa13)*dble(xyzzyaaaa13)*xyzzyaaag&
&13
xyzzyaaag13=dble(xyzzyaaaa13)*xyzzyaaaf13
xyzzyaaac13=xyzzyaaac13+alpha(xyzzyaaaa13)*xyzzyaaag13
xyzzyaaaf13=xyzzyaaaf13*r
xyzzyaaab13=xyzzyaaab13+alpha(xyzzyaaaa13)*xyzzyaaaf13
enddo
xyzzyaaah13=exp_protect(xyzzyaaab13)
xyzzyaaae13=1.d0/r
if(xyzzyaaah13+c/=0.d0)then
xyzzyaacd1=-0.5d0*xyzzyaaah13/(xyzzyaaah13+c)*(xyzzyaaac13*(2.d0*xyzzy&
&aaae13+xyzzyaaac13)+xyzzyaaad13)-z*xyzzyaaae13
else
xyzzyaacd1=0.d0
endif
else
xyzzyaaah13=exp_protect(alpha(0))
if(xyzzyaaah13+c/=0.d0)then
xyzzyaacd1=-0.5d0*xyzzyaaah13/(c+xyzzyaaah13)*(6.d0*alpha(2)+alpha(1)*&
&alpha(1))
else
xyzzyaacd1=0.d0
endif
endif
end function xyzzyaacd1
real(dp) function xyzzyaace1(z,alpha,c,rc)
implicit none
real(dp),intent(in) :: z,alpha(0:4),rc,c
integer xyzzyaaaa14,xyzzyaaab14
integer,parameter :: xyzzyaaac14=100,xyzzyaaad14=100
real(dp) xyzzyaaae14,xyzzyaaaf14,xyzzyaaag14,xyzzyaaah14,xyzzyaaai14,x&
&yzzyaaaj14(0:8),xyzzyaaak14(0:8),xyzzyaaal14(0:8),xyzzyaaam14,xyzzyaa&
&an14,xyzzyaaao14,xyzzyaaap14,xyzzyaaaq14,xyzzyaaar14,xyzzyaaas14,xyzz&
&yaaat14,xyzzyaaau14,xyzzyaaav14,xyzzyaaaw14,xyzzyaaax14,xyzzyaaay14,x&
&yzzyaaaz14,xyzzyaaba14,xyzzyaabb14,xyzzyaabc14,xyzzyaabd14,xyzzyaabe1&
&4,xyzzyaabf14,xyzzyaabg14,xyzzyaabh14,xyzzyaabi14,xyzzyaabj14(0:8),xy&
&zzyaabk14,xyzzyaabl14
real(dp),parameter :: xyzzyaabm14=1.d-12
logical xyzzyaabn14
xyzzyaabn14=(nint(z)==1)
if(xyzzyaabn14)then
xyzzyaaaf14=xyzzyaacd1(z,alpha,c,rc)
xyzzyaabj14=(/xyzzyaaaf14,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/)
else
xyzzyaaae14=(xyzzyaacd1(z,alpha,c,rc)-xyzzyaacc1(z,0.d0,rc))/z**2
xyzzyaabj14=z**2*(/xyzzyaaae14,0.d0,xyzzyaaae1,xyzzyaaaf1,xyzzyaaag1,x&
&yzzyaaah1,xyzzyaaai1,xyzzyaaaj1,xyzzyaaak1/)
endif
xyzzyaace1=0.d0
xyzzyaaai14=-1.d0
xyzzyaabk14=rc/dble(xyzzyaaac14-1)
do xyzzyaaaa14=0,xyzzyaaac14-1
xyzzyaaah14=dble(xyzzyaaaa14)*xyzzyaabk14
if(.not.xyzzyaabn14)xyzzyaaaf14=xyzzyaacc1(z,xyzzyaaae14,xyzzyaaah14)
xyzzyaaag14=(xyzzyaacd1(z,alpha,c,xyzzyaaah14)-xyzzyaaaf14)**2
if(xyzzyaaag14>xyzzyaace1)then
xyzzyaaai14=xyzzyaaah14
xyzzyaace1=xyzzyaaag14
endif
enddo
if(xyzzyaaai14<0.d0)call errstop('DELTA_E_L_SQ_MAX','Cannot locate max&
&imum by brute force.')
if(xyzzyaaai14>0.d0)then
xyzzyaaah14=xyzzyaaai14
else
xyzzyaaah14=xyzzyaabk14
endif
xyzzyaabi14=xyzzyaabm14*rc
do xyzzyaaab14=1,xyzzyaaad14
xyzzyaaaj14(0)=1.d0
xyzzyaaak14(0)=0.d0
xyzzyaaal14(0)=0.d0
do xyzzyaaaa14=1,8
xyzzyaaaj14(xyzzyaaaa14)=xyzzyaaaj14(xyzzyaaaa14-1)*xyzzyaaah14
xyzzyaaak14(xyzzyaaaa14)=xyzzyaaaj14(xyzzyaaaa14-1)*dble(xyzzyaaaa14)
xyzzyaaal14(xyzzyaaaa14)=xyzzyaaak14(xyzzyaaaa14-1)*dble(xyzzyaaaa14)
enddo
xyzzyaaam14=ddot(7,xyzzyaabj14(2),1,xyzzyaaak14(2),1)
xyzzyaaan14=ddot(7,xyzzyaabj14(2),1,xyzzyaaal14(2),1)
xyzzyaaao14=ddot(5,alpha(0),1,xyzzyaaaj14(0),1)
xyzzyaaap14=ddot(4,alpha(1),1,xyzzyaaak14(1),1)
xyzzyaaaq14=ddot(3,alpha(2),1,xyzzyaaal14(2),1)
xyzzyaaar14=6.d0*alpha(3)+24.d0*alpha(4)*xyzzyaaah14
xyzzyaaas14=24.d0*alpha(4)
xyzzyaabl14=c*exp_protect(-xyzzyaaao14)+1.d0
if(xyzzyaabl14/=0.d0)then
xyzzyaaat14=1.d0/xyzzyaabl14
else
xyzzyaaat14=0.d0
endif
xyzzyaaau14=xyzzyaaat14*(1.d0-xyzzyaaat14)*xyzzyaaap14
xyzzyaaav14=(1.d0-2.d0*xyzzyaaat14)*xyzzyaaau14*xyzzyaaap14+xyzzyaaat1&
&4*(1.d0-xyzzyaaat14)*xyzzyaaaq14
xyzzyaaaw14=1.d0/xyzzyaaah14
xyzzyaaax14=xyzzyaaaw14*xyzzyaaaw14
xyzzyaaay14=xyzzyaaax14*xyzzyaaaw14
xyzzyaaaz14=2.d0*xyzzyaaaw14+xyzzyaaap14
xyzzyaaba14=-2.d0*xyzzyaaax14+xyzzyaaaq14
xyzzyaabb14=4.d0*xyzzyaaay14+xyzzyaaar14
xyzzyaabc14=xyzzyaaap14*xyzzyaaaz14+xyzzyaaaq14
xyzzyaabd14=xyzzyaaaq14*xyzzyaaaz14+xyzzyaaap14*xyzzyaaba14+xyzzyaaar1&
&4
xyzzyaabe14=xyzzyaaar14*xyzzyaaaz14+2.d0*xyzzyaaaq14*xyzzyaaba14+xyzzy&
&aaap14*xyzzyaabb14+xyzzyaaas14
xyzzyaabf14=-0.5d0*(xyzzyaaat14*xyzzyaabd14+xyzzyaaau14*xyzzyaabc14)+z&
&*xyzzyaaax14
xyzzyaabg14=-0.5d0*(xyzzyaaat14*xyzzyaabe14+xyzzyaaav14*xyzzyaabc14)-x&
&yzzyaaau14*xyzzyaabd14-2.d0*z*xyzzyaaay14
if(xyzzyaaan14==xyzzyaabg14)exit
xyzzyaabh14=xyzzyaaah14
xyzzyaaah14=xyzzyaaah14-(xyzzyaaam14-xyzzyaabf14)/(xyzzyaaan14-xyzzyaa&
&bg14)
if(xyzzyaaah14>rc)xyzzyaaah14=rc
if(xyzzyaaah14<=0.d0)xyzzyaaah14=0.5d0*xyzzyaabh14
if(abs(xyzzyaaah14-xyzzyaabh14)<xyzzyaabi14)exit
enddo
if(.not.xyzzyaabn14)xyzzyaaaf14=xyzzyaacc1(z,xyzzyaaae14,xyzzyaaah14)
xyzzyaaag14=(xyzzyaacd1(z,alpha,c,xyzzyaaah14)-xyzzyaaaf14)**2
if(xyzzyaace1-xyzzyaaag14>xyzzyaabm14*xyzzyaace1.and.xyzzyaaai14>0.d0 &
&.and.xyzzyaabf1>4.and.am_master)then
call wout('Brute force: DE2=',xyzzyaace1)
call wout('NR:          DE2=',xyzzyaaag14)
call errwarn('DELTA_E_L_SQ_MAX','Newton-Raphson iteration failed.')
endif
xyzzyaace1=max(xyzzyaaag14,xyzzyaace1)
end function xyzzyaace1
real(dp) function xyzzyaacf1(iorb,jion,r)
implicit none
integer,intent(in) :: iorb,jion
real(dp),intent(in) :: r
real(dp) xyzzyaaaa15,xyzzyaaab15,xyzzyaaac15,xyzzyaaad15
if(r>0.d0)then
xyzzyaaaa15=1.d0/r
call eval_spline_reg(nradgrid,radgrid(1:nradgrid,jion),xyzzyaabc1(1:nr&
&adgrid,iorb,jion),xyzzyaabd1(1:nradgrid,iorb,jion),r,.true.,.true.,xy&
&zzyaaab15,xyzzyaaac15,xyzzyaaad15)
if(xyzzyaaab15/=0.d0)then
xyzzyaacf1=(-0.5d0*xyzzyaaad15-xyzzyaaac15*xyzzyaaaa15)/xyzzyaaab15-zi&
&on(iontype(ae_index(jion)))*xyzzyaaaa15
else
xyzzyaacf1=0.d0
endif
else
xyzzyaacf1=0.d0
endif
end function xyzzyaacf1
subroutine xyzzyaacg1(iorb,jion,no_nodes,rnode,rc)
implicit none
integer,intent(in) :: iorb,jion,no_nodes
real(dp),intent(in) :: rnode(xyzzyaabh1)
real(dp),intent(out) :: rc
integer xyzzyaaaa16
integer,parameter :: xyzzyaaab16=100,xyzzyaaac16=20
real(dp) xyzzyaaad16,xyzzyaaae16,xyzzyaaaf16,xyzzyaaag16,xyzzyaaah16,x&
&yzzyaaai16,xyzzyaaaj16,xyzzyaaak16,xyzzyaaal16
real(dp),parameter :: xyzzyaaam16=1.d-8
logical xyzzyaaan16
xyzzyaaad16=zion(iontype(ae_index(jion)))
xyzzyaaan16=(nint(xyzzyaaad16)==1)
if(xyzzyaaan16)then
xyzzyaaae16=xyzzyaacf1(iorb,jion,xyzzyaaap1(jion))
else
xyzzyaaaf16=(xyzzyaacf1(iorb,jion,xyzzyaaap1(jion))-xyzzyaacc1(xyzzyaa&
&ad16,0.d0,xyzzyaaap1(jion)))/xyzzyaaad16**2
endif
if(no_nodes==0)then
xyzzyaaal16=xyzzyaaap1(jion)
else
xyzzyaaal16=min(xyzzyaaap1(jion),xyzzyaabk1*minval(rnode(1:no_nodes)))
endif
xyzzyaaag16=xyzzyaaad16*xyzzyaaad16/xyzzyaabg1
xyzzyaaah16=xyzzyaaal16/dble(xyzzyaaab16)
iloop: do xyzzyaaaa16=xyzzyaaab16,1,-1
rc=dble(xyzzyaaaa16)*xyzzyaaah16
if(.not.xyzzyaaan16)xyzzyaaae16=xyzzyaacc1(xyzzyaaad16,xyzzyaaaf16,rc)
if(abs(xyzzyaacf1(iorb,jion,rc)-xyzzyaaae16)>=xyzzyaaag16)exit
enddo iloop
xyzzyaaak16=xyzzyaaam16*xyzzyaaap1(jion)
if(rc<xyzzyaaal16-xyzzyaaak16)then
xyzzyaaai16=rc
xyzzyaaaj16=min(xyzzyaaap1(jion),xyzzyaaai16+xyzzyaaah16)
do xyzzyaaaa16=1,xyzzyaaac16
rc=0.5d0*(xyzzyaaai16+xyzzyaaaj16)
if(xyzzyaaaj16-xyzzyaaai16<xyzzyaaak16)exit
if(.not.xyzzyaaan16)xyzzyaaae16=xyzzyaacc1(xyzzyaaad16,xyzzyaaaf16,rc)
if(abs(xyzzyaacf1(iorb,jion,rc)-xyzzyaaae16)>=xyzzyaaag16)then
xyzzyaaai16=rc
else
xyzzyaaaj16=rc
endif
enddo
endif
end subroutine xyzzyaacg1
subroutine xyzzyaach1(iorb,jion,first_setup,no_nodes,rnode)
implicit none
integer,intent(in) :: iorb,jion
integer,intent(out) :: no_nodes
real(dp),intent(out) :: rnode(xyzzyaabh1)
logical,intent(in) :: first_setup
integer xyzzyaaaa17,xyzzyaaab17,xyzzyaaac17
integer,parameter :: xyzzyaaad17=1000,xyzzyaaae17=20
real(dp) xyzzyaaaf17,xyzzyaaag17,xyzzyaaah17,xyzzyaaai17,xyzzyaaaj17,x&
&yzzyaaak17,xyzzyaaal17,xyzzyaaam17,xyzzyaaan17,xyzzyaaao17
real(dp),parameter :: xyzzyaaap17=1.d-8
logical xyzzyaaaq17
xyzzyaaan17=xyzzyaaap17*xyzzyaaap1(jion)
xyzzyaaaf17=xyzzyaaap1(jion)/dble(xyzzyaaad17)
xyzzyaaaa17=0
xyzzyaaah17=xyzzyaabc1(1,iorb,jion)
do xyzzyaaab17=0,xyzzyaaad17
xyzzyaaag17=dble(xyzzyaaab17)*xyzzyaaaf17
call eval_spline_reg(nradgrid,radgrid(1:nradgrid,jion),xyzzyaabc1(1:nr&
&adgrid,iorb,jion),xyzzyaabd1(1:nradgrid,iorb,jion),xyzzyaaag17,.true.&
&,.false.,xyzzyaaai17,xyzzyaaaj17,xyzzyaaak17)
if((xyzzyaaah17>0.d0.and.xyzzyaaai17<=0.d0).or.(xyzzyaaah17<0.d0.and.x&
&yzzyaaai17>=0.d0))then
xyzzyaaaa17=xyzzyaaaa17+1
if(xyzzyaaaa17>xyzzyaabh1)call errstop('DETERMINE_NODES','Need to incr&
&ease the maxnodes parameter.')
if(xyzzyaaai17/=xyzzyaaah17)then
xyzzyaaal17=xyzzyaaag17-xyzzyaaaf17-xyzzyaaah17*xyzzyaaaf17/(xyzzyaaai&
&17-xyzzyaaah17)
else
xyzzyaaal17=xyzzyaaag17-0.5d0*xyzzyaaaf17
endif
xyzzyaaah17=xyzzyaaai17
do xyzzyaaac17=1,xyzzyaaae17
call eval_spline_reg(nradgrid,radgrid(1:nradgrid,jion),xyzzyaabc1(1:nr&
&adgrid,iorb,jion),xyzzyaabd1(1:nradgrid,iorb,jion),xyzzyaaal17,.true.&
&,.true.,xyzzyaaai17,xyzzyaaaj17,xyzzyaaak17)
if(xyzzyaaaj17==0.d0)exit
xyzzyaaam17=-xyzzyaaai17/xyzzyaaaj17
xyzzyaaal17=xyzzyaaal17+xyzzyaaam17
if(abs(xyzzyaaam17)<xyzzyaaan17)exit
enddo
rnode(xyzzyaaaa17)=max(xyzzyaaal17,0.d0)
else
xyzzyaaah17=xyzzyaaai17
endif
enddo
no_nodes=xyzzyaaaa17
do
xyzzyaaaq17=.true.
do xyzzyaaab17=1,no_nodes-1
if(rnode(xyzzyaaab17)>rnode(xyzzyaaab17+1))then
xyzzyaaaq17=.false.
xyzzyaaao17=rnode(xyzzyaaab17)
rnode(xyzzyaaab17)=rnode(xyzzyaaab17+1)
rnode(xyzzyaaab17+1)=xyzzyaaao17
endif
enddo
if(xyzzyaaaq17)exit
enddo
if(am_master.and.xyzzyaabf1>3.and.first_setup.and.xyzzyaaaa17>0)then
call wout('   This orbital has '//trim(i2s(xyzzyaaaa17))//' nodes in [&
&0,rcmax].')
if(xyzzyaabf1>4)then
do xyzzyaaab17=1,no_nodes
call eval_spline_reg(nradgrid,radgrid(1:nradgrid,jion),xyzzyaabc1(1:nr&
&adgrid,iorb,jion),xyzzyaabd1(1:nradgrid,iorb,jion),rnode(xyzzyaaab17)&
&,.true.,.false.,xyzzyaaai17,xyzzyaaaj17,xyzzyaaak17)
call wout('   Node ; orbbar : ',(/rnode(xyzzyaaab17),xyzzyaaai17/),rfm&
&t='(es20.12)')
enddo
endif
endif
end subroutine xyzzyaach1
subroutine setup_spline_reg(n,xarr,yarr,d2yarr,yp1,ypn)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: xarr(n),yarr(n)
real(dp),intent(in),optional :: yp1,ypn
real(dp),intent(out) :: d2yarr(n)
integer xyzzyaaaa18,xyzzyaaab18
real(dp) xyzzyaaac18,xyzzyaaad18,xyzzyaaae18,xyzzyaaaf18
real(dp),parameter :: xyzzyaaag18=1.d-8
real(dp),allocatable :: xyzzyaaah18(:)
xyzzyaaac18=xarr(2)-xarr(1)
xyzzyaaad18=xyzzyaaag18*xyzzyaaac18
do xyzzyaaaa18=3,n
if(abs(xarr(xyzzyaaaa18)-xarr(xyzzyaaaa18-1)-xyzzyaaac18)>xyzzyaaad18)&
&call errstop('SETUP_SPLINE_REG','Bug: x data are not regular.')
enddo
xyzzyaaae18=1.d0/xyzzyaaac18
xyzzyaaaf18=xyzzyaaae18*xyzzyaaae18
allocate(xyzzyaaah18(n),stat=xyzzyaaab18)
call check_alloc(xyzzyaaab18,'SETUP_SPLINE','')
if(present(yp1))then
xyzzyaaah18(1)=3.d0*xyzzyaaae18*((yarr(2)-yarr(1))*xyzzyaaae18-yp1)
d2yarr(1)=-0.5d0
else
xyzzyaaah18(1)=0.d0
d2yarr(1)=0.d0
endif
do xyzzyaaaa18=2,n-1
d2yarr(xyzzyaaaa18)=-1.d0/(d2yarr(xyzzyaaaa18-1)+4.d0)
xyzzyaaah18(xyzzyaaaa18)=(6.d0*(yarr(xyzzyaaaa18)+yarr(xyzzyaaaa18)-ya&
&rr(xyzzyaaaa18+1)-yarr(xyzzyaaaa18-1))*xyzzyaaaf18+xyzzyaaah18(xyzzya&
&aaa18-1))*d2yarr(xyzzyaaaa18)
enddo
if(present(ypn))then
xyzzyaaah18(n)=3.d0*xyzzyaaae18*(ypn-(yarr(n)-yarr(n-1))*xyzzyaaae18)
d2yarr(n)=(2.d0*xyzzyaaah18(n)-xyzzyaaah18(n-1))/(d2yarr(n-1)+2.d0)
else
xyzzyaaah18(n)=0.d0
d2yarr(n)=0.d0
endif
do xyzzyaaaa18=n-1,1,-1
d2yarr(xyzzyaaaa18)=d2yarr(xyzzyaaaa18)*d2yarr(xyzzyaaaa18+1)+xyzzyaaa&
&h18(xyzzyaaaa18)
enddo
deallocate(xyzzyaaah18)
end subroutine setup_spline_reg
subroutine eval_spline_reg(n,xarr,yarr,d2yarr,x,val,fsd,y,dy,d2y)
use slaarnaag,only : sixth
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: xarr(n),x,yarr(n),d2yarr(n)
real(dp),intent(out) :: y,dy,d2y
logical,intent(in) :: val,fsd
integer xyzzyaaaa19,xyzzyaaab19
real(dp) xyzzyaaac19,xyzzyaaad19,xyzzyaaae19,xyzzyaaaf19
xyzzyaaae19=xarr(2)-xarr(1)
xyzzyaaaf19=1.d0/xyzzyaaae19
xyzzyaaab19=min(max(int((x-xarr(1))*xyzzyaaaf19)+1,1),n-1)
xyzzyaaaa19=xyzzyaaab19+1
xyzzyaaac19=(xarr(xyzzyaaaa19)-x)*xyzzyaaaf19
xyzzyaaad19=1.d0-xyzzyaaac19
if(val)y=xyzzyaaac19*yarr(xyzzyaaab19)+xyzzyaaad19*yarr(xyzzyaaaa19)+(&
&xyzzyaaac19*(xyzzyaaac19*xyzzyaaac19-1.d0)*d2yarr(xyzzyaaab19) +xyzzy&
&aaad19*(xyzzyaaad19*xyzzyaaad19-1.d0)*d2yarr(xyzzyaaaa19))*xyzzyaaae1&
&9*xyzzyaaae19*sixth
if(fsd)then
dy=(yarr(xyzzyaaaa19)-yarr(xyzzyaaab19))*xyzzyaaaf19+xyzzyaaae19*((0.5&
&d0*xyzzyaaad19*xyzzyaaad19-sixth)*d2yarr(xyzzyaaaa19) -(0.5d0*xyzzyaa&
&ac19*xyzzyaaac19-sixth)*d2yarr(xyzzyaaab19))
d2y=xyzzyaaad19*d2yarr(xyzzyaaaa19)+xyzzyaaac19*d2yarr(xyzzyaaab19)
endif
end subroutine eval_spline_reg
subroutine gpcc_init_buffer(nbuf)
implicit none
integer,intent(in) :: nbuf
integer xyzzyaaaa20
xyzzyaabx1=xyzzyaaac1*naeions_prim
xyzzyaaby1=nradgrid*xyzzyaabx1
xyzzyaabz1=5*xyzzyaabx1
allocate(xyzzyaabp1(naeions_prim,0:nbuf),xyzzyaabq1(xyzzyaabx1,0:nbuf)&
&,xyzzyaabr1(xyzzyaaby1,0:nbuf),xyzzyaabs1(xyzzyaaby1,0:nbuf),xyzzyaab&
&u1(xyzzyaabz1,0:nbuf),xyzzyaabv1(xyzzyaabx1,0:nbuf),stat=xyzzyaaaa20)
call check_alloc(xyzzyaaaa20,'GPCC_INIT_BUFFER','max_rc_sq_buf')
xyzzyaabp1=0.d0
xyzzyaabq1=0.d0
xyzzyaabr1=0.d0
xyzzyaabs1=0.d0
xyzzyaabu1=0.d0
xyzzyaabv1=0.d0
if(xyzzyaabj1)then
allocate(xyzzyaabt1(xyzzyaabx1,0:nbuf),stat=xyzzyaaaa20)
call check_alloc(xyzzyaaaa20,'GPCC_INIT_BUFFER','sgn_buf')
xyzzyaabt1=0.d0
else
allocate(xyzzyaabw1(xyzzyaabx1,nbuf),stat=xyzzyaaaa20)
call check_alloc(xyzzyaaaa20,'GPCC_INIT_BUFFER','phase_buf')
xyzzyaabw1=czero
endif
end subroutine gpcc_init_buffer
subroutine gpcc_finalize_buffer
implicit none
deallocate(xyzzyaabp1,xyzzyaabq1,xyzzyaabr1,xyzzyaabs1,xyzzyaabu1,xyzz&
&yaabv1)
if(xyzzyaabj1)then
deallocate(xyzzyaabt1)
else
deallocate(xyzzyaabw1)
endif
end subroutine gpcc_finalize_buffer
subroutine gpcc_load_from_buffer(ibuf)
implicit none
integer,intent(in) :: ibuf
call dcopy(naeions_prim,xyzzyaabp1(1,ibuf),1,xyzzyaaao1(1),1)
call dcopy(xyzzyaabx1,xyzzyaabq1(1,ibuf),1,xyzzyaaan1(1,1),1)
call dcopy(xyzzyaaby1,xyzzyaabr1(1,ibuf),1,xyzzyaabc1(1,1,1),1)
call dcopy(xyzzyaaby1,xyzzyaabs1(1,ibuf),1,xyzzyaabd1(1,1,1),1)
if(xyzzyaabj1)then
call dcopy(xyzzyaabx1,xyzzyaabt1(1,ibuf),1,xyzzyaaar1(1,1),1)
else
call zcopy(xyzzyaabx1,xyzzyaabw1(1,ibuf),1,xyzzyaaas1(1,1),1)
endif
call dcopy(xyzzyaabz1,xyzzyaabu1(1,ibuf),1,xyzzyaaal1(0,1,1),1)
call dcopy(xyzzyaabx1,xyzzyaabv1(1,ibuf),1,xyzzyaaam1(1,1),1)
end subroutine gpcc_load_from_buffer
subroutine gpcc_save_to_buffer(ibuf)
implicit none
integer,intent(in) :: ibuf
call dcopy(naeions_prim,xyzzyaaao1(1),1,xyzzyaabp1(1,ibuf),1)
call dcopy(xyzzyaabx1,xyzzyaaan1(1,1),1,xyzzyaabq1(1,ibuf),1)
call dcopy(xyzzyaaby1,xyzzyaabc1(1,1,1),1,xyzzyaabr1(1,ibuf),1)
call dcopy(xyzzyaaby1,xyzzyaabd1(1,1,1),1,xyzzyaabs1(1,ibuf),1)
if(xyzzyaabj1)then
call dcopy(xyzzyaabx1,xyzzyaaar1(1,1),1,xyzzyaabt1(1,ibuf),1)
else
call zcopy(xyzzyaabx1,xyzzyaaas1(1,1),1,xyzzyaabw1(1,ibuf),1)
endif
call dcopy(xyzzyaabz1,xyzzyaaal1(0,1,1),1,xyzzyaabu1(1,ibuf),1)
call dcopy(xyzzyaabx1,xyzzyaaam1(1,1),1,xyzzyaabv1(1,ibuf),1)
end subroutine gpcc_save_to_buffer
end module slaarnabi
function objfn_inner(p0)
use dsp, only : dp
use slaarnabi,only : determine_poly_inner
implicit none
real(dp),intent(in) :: p0
real(dp) objfn_inner
call determine_poly_inner(p0,.false.,objfn_inner)
end function objfn_inner
function objfn_outer(rc)
use dsp, only : dp
use slaarnabi,only : determine_poly_outer
implicit none
real(dp),intent(in) :: rc
real(dp) objfn_outer
call determine_poly_outer(rc,.false.,objfn_outer)
end function objfn_outer

module slaarnaas
use dsp
use slaarnaag,   only : pi,twopi,fourpi,third,half,one_over_root_two,e&
&uler,two_euler,root_three_over_two
use file_utils,  only : open_units
use format_utils,only : wout,r2s,r2s2,i2s,l2s,wordwrap,write_list_int,&
&capitalize
use slaarnaat,   only : gautol
use slaarnabg,    only : nbasis,periodicity,a1,a2,a3,scell_matrix,scel&
&l_mat_inv,npcells,pa1,pa2,pa3,b1,b2,b3,isperiodic,pamat,pbmat,pb1,pb2&
&,pb3,amat,ainv,dimensionality,area,wigner_seitz_radius,pvolume,nitot,&
&rion,volume,bmat,binv
use slaarnabp,        only : wf_nd,wf_np,wf_nm,wf_d,mdet_max_mods
use slaarnabt,   only : dcopy,inverse3,dscal,ddot,quicksort,erfc,lu_de&
&com,lu_solve_once
use parallel,    only : am_master
use run_control, only : errstop,errstop_master,timer,check_alloc
use slaarnaci,         only : print_sdwdata,sdw_setup,sdw_kvec,sdw_orb&
&_theta
use store,       only : complex_wf,heg_orbtype,ndet,nspin,levels_ssing&
&les,nemax,heg_nele,constant_energy,netot,real1_complex2,nele,use_back&
&flow,self_term,inv_pmass,pcharge,which_ssingle,update_by_column,heg_l&
&ayer,popp_spin,pspin,no_ssingles,interaction,mc_twist_freq,nspin_ndet&
&,mc_twist_av,pmass,have_biex3pot,fix_holes,orb_norm,k_offset,site_pos&
&_in_cell,spin_dep_wc,wc_gauss_exp,type_wc,heg_nlayers,heg_ylayer,heg_&
&zlayer,which_spin,heg_slatt
implicit none
private
public init_free_orbs,free_biex1_orb_eval,free_crystal_orb_eval,free_p&
&airing_orb_eval,free_fluid_orb_eval,heg_mpc_long,half_mpc_constant,ex&
&cite_heg,write_heg,corr_heg_required,heg_nbasis,heg_cell,heg_wigner_b&
&asis,heg_crystal_type,r_s,rho,heg_ferro,heg_nosites,heg_crystal_sites&
&,heg_repeat,calc_hf_energies,mc_twist_offset,me_biex3,mh_biex3,mu_bie&
&x3,xx_sep,setup_freeorb_params,finish_freeorb_params,get_freeorb_para&
&ms,put_freeorb_params,harmwire_b,get_free_orbmap,get_free_orbdesc,get&
&_free_ndesc,free_norb,gs_kvec,eval_finite_hf_energies,ppmcta_partial_&
&sum,ppmcta_fit
integer xyzzyaaaa1,heg_nbasis,xyzzyaaab1,heg_repeat,xyzzyaaac1,xyzzyaa&
&ad1,xyzzyaaae1(3),xyzzyaaaf1,xyzzyaaag1(3),free_norb
integer,allocatable :: xyzzyaaah1(:),xyzzyaaai1(:),heg_nosites(:),xyzz&
&yaaaj1(:,:,:),xyzzyaaak1(:),xyzzyaaal1(:,:),xyzzyaaam1(:,:),xyzzyaaan&
&1(:),xyzzyaaao1(:,:,:)
real(dp) xyzzyaaap1,half_mpc_constant,xyzzyaaaq1,xyzzyaaar1,xyzzyaaas1&
&,xyzzyaaat1,xyzzyaaau1,xyzzyaaav1,r_s,rho,heg_cell(3,3)
real(dp) :: me_biex3=-1.d0,mh_biex3=-1.d0,mu_biex3=-1.d0,xx_sep=-1.d0,&
&harmwire_b=-1.d0
real(dp),allocatable :: gs_kvec(:,:),xyzzyaaaw1(:),heg_wigner_basis(:,&
&:),xyzzyaaax1(:,:),xyzzyaaay1(:,:,:),heg_crystal_sites(:,:,:),xyzzyaa&
&az1(:,:),xyzzyaaba1(:),xyzzyaabb1(:,:)
complex(dp),allocatable :: xyzzyaabc1(:),xyzzyaabd1(:),xyzzyaabe1(:)
logical corr_heg_required,xyzzyaabf1,xyzzyaabg1,xyzzyaabh1,xyzzyaabi1,&
&xyzzyaabj1,xyzzyaabk1,xyzzyaabl1,xyzzyaabm1,xyzzyaabn1,xyzzyaabo1
logical,allocatable :: heg_ferro(:)
character(11) heg_crystal_type
character(80) title,tmpr,tmpr2
integer xyzzyaabp1,xyzzyaabq1,xyzzyaabr1,xyzzyaabs1,xyzzyaabt1,xyzzyaa&
&bu1,xyzzyaabv1,xyzzyaabw1,xyzzyaabx1,xyzzyaaby1,xyzzyaabz1,xyzzyaaca1&
&,xyzzyaacb1
integer,allocatable :: xyzzyaacc1(:,:),xyzzyaacd1(:,:),xyzzyaace1(:),x&
&yzzyaacf1(:,:),xyzzyaacg1(:,:),xyzzyaach1(:,:),xyzzyaaci1(:,:),xyzzya&
&acj1(:,:),xyzzyaack1(:,:),xyzzyaacl1(:,:,:),xyzzyaacm1(:,:),xyzzyaacn&
&1(:,:)
real(dp) xyzzyaaco1,xyzzyaacp1,xyzzyaacq1,xyzzyaacr1,xyzzyaacs1(3,3),x&
&yzzyaact1(3)
real(dp),allocatable :: xyzzyaacu1(:,:,:),xyzzyaacv1(:,:),xyzzyaacw1(:&
&,:,:),xyzzyaacx1(:,:),xyzzyaacy1(:),xyzzyaacz1(:),xyzzyaada1(:),xyzzy&
&aadb1(:,:,:),xyzzyaadc1(:,:,:),xyzzyaadd1(:,:),xyzzyaade1(:,:),xyzzya&
&adf1(:,:,:,:),xyzzyaadg1(:,:),xyzzyaadh1(:,:)
logical xyzzyaadi1,xyzzyaadj1,xyzzyaadk1,xyzzyaadl1,xyzzyaadm1
logical :: xyzzyaadn1=.false.
integer xyzzyaado1,xyzzyaadp1
integer,allocatable :: xyzzyaadq1(:),xyzzyaadr1(:),xyzzyaads1(:)
real(dp),allocatable :: xyzzyaadt1(:,:),xyzzyaadu1(:),xyzzyaadv1(:,:)
logical excite_heg,xyzzyaadw1,xyzzyaadx1
logical,allocatable :: xyzzyaady1(:)
integer,allocatable :: xyzzyaadz1(:,:),xyzzyaaea1(:)
real(dp),allocatable :: xyzzyaaeb1(:,:),xyzzyaaec1(:)
integer,parameter :: xyzzyaaed1=9999
integer xyzzyaaee1,xyzzyaaef1
integer xyzzyaaeg1,xyzzyaaeh1,xyzzyaaei1,xyzzyaaej1,xyzzyaaek1,xyzzyaa&
&el1,xyzzyaaem1,xyzzyaaen1,xyzzyaaeo1
real(dp),allocatable :: xyzzyaaep1(:,:,:,:),xyzzyaaeq1(:,:,:),xyzzyaae&
&r1(:,:,:,:),xyzzyaaes1(:),xyzzyaaet1(:),xyzzyaaeu1(:),xyzzyaaev1(:,:,&
&:,:),xyzzyaaew1(:,:,:,:),xyzzyaaex1(:,:,:),xyzzyaaey1(:,:,:),xyzzyaae&
&z1(:,:,:,:,:),xyzzyaafa1(:,:,:),xyzzyaafb1(:,:,:),xyzzyaafc1(:,:),xyz&
&zyaafd1(:,:),xyzzyaafe1(:,:,:)
integer xyzzyaaff1,xyzzyaafg1
integer,allocatable :: xyzzyaafh1(:,:)
contains
subroutine init_free_orbs
implicit none
integer xyzzyaaaa2,xyzzyaaab2,d,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2,xyzzy&
&aaaf2,xyzzyaaag2,xyzzyaaah2,xyzzyaaai2,xyzzyaaaj2,xyzzyaaak2,xyzzyaaa&
&l2,xyzzyaaam2,xyzzyaaan2
real(dp) xyzzyaaao2,xyzzyaaap2,xyzzyaaaq2
logical xyzzyaaar2,xyzzyaaas2,xyzzyaaat2
character(80) term_label,char_80,temp1,temp2,temp3,temp4,temp5
xyzzyaabl1=complex_wf
type_wc=0
xyzzyaaat2=.false.
xyzzyaadi1=.false.
xyzzyaadj1=.false.
xyzzyaadk1=.false.
xyzzyaadl1=.false.
xyzzyaabk1=.true.
xyzzyaadm1=.false.
xyzzyaabf1=any(heg_orbtype==1)
xyzzyaabg1=any(heg_orbtype<0)
xyzzyaabh1=any(heg_orbtype==2)
xyzzyaabj1=any(heg_orbtype==3)
xyzzyaabm1=any(heg_orbtype==4)
xyzzyaabn1=any(heg_orbtype==5)
xyzzyaabo1=any(heg_orbtype==6)
xyzzyaabi1=(dimensionality==0)
excite_heg=any(wf_nd(:,:)>0)
xyzzyaadx1=excite_heg.or.(free_norb>maxval(heg_nele))
if(any(wf_np(:,:)>0).or.any(wf_nm(:,:)>0))call errstop_master('INIT_FR&
&EE_ORBS','Subtractions and additions are not supported.  Set complex_&
&wf=T, supply the approriate number of electrons, and make any necessa&
&ry excitations.')
xyzzyaadw1=(ndet>1)
if(am_master)then
call wout('Free-particle orbitals')
call wout('======================')
write(tmpr,'(a,t10,a,t15,a,t25,a,t40,a)')'MD term','Det','Particle','N&
&o. particles','Type of free orbitals'
call wout(tmpr)
call wout(repeat('-',64))
do xyzzyaaad2=1,ndet
temp1=trim(i2s(xyzzyaaad2))
do xyzzyaaaa2=1,nspin
xyzzyaaas2=.true.
temp2=trim(i2s(xyzzyaaaa2))
temp3=trim(i2s(xyzzyaaaa2))
temp4=trim(i2s(heg_nele(xyzzyaaaa2)))
select case(heg_orbtype(xyzzyaaaa2,xyzzyaaad2))
case(0)
temp5='None'
case(1)
if(xyzzyaabl1)then
temp5='Complex plane waves'
else
temp5='Real plane waves'
endif
if(wf_nd(xyzzyaaad2,xyzzyaaaa2)>0)temp5=trim(temp5)//' - excited(' //t&
&rim(i2s(wf_nd(xyzzyaaad2,xyzzyaaaa2)))//')'
case(2)
temp5='Wigner Crystal'
case(3)
temp5='Spin-density waves'
case(4)
temp5='Biexciton (biex1)'
case(5)
temp5='Biexciton (biex2)'
case(6)
temp5='Biexciton (biex3)'
case(7)
temp5='Exciton molecule'
case(101:)
temp5='Expot orbitals - SET #'//trim(i2s(heg_orbtype(xyzzyaaaa2,xyzzya&
&aad2)-100))
case default
temp5='Pairing'
if(xyzzyaaaa2<-heg_orbtype(xyzzyaaaa2,xyzzyaaad2))then
temp3=trim(i2s(xyzzyaaaa2))//' & '//trim(i2s(-heg_orbtype(xyzzyaaaa2,x&
&yzzyaaad2)))
temp4=trim(i2s(heg_nele(xyzzyaaaa2)))//' + '//trim(i2s(heg_nele(-heg_o&
&rbtype(xyzzyaaaa2,xyzzyaaad2))))
else
xyzzyaaas2=.false.
endif
end select
if(xyzzyaaas2)then
write(tmpr,'(a,t10,a,t15,a,t25,a,t40,a)')trim(temp1),trim(temp2),trim(&
&temp3),trim(temp4),trim(temp5)
call wout(tmpr)
temp1=''
endif
enddo
enddo
call wout(repeat('-',64))
call wout()
endif
if(xyzzyaabf1.and.any(heg_orbtype(:,1)/=1.and.heg_nele>0))call errstop&
&_master('INIT_FREE_ORBS','If a fluid phase is to be simulated then de&
&terminant 1 should contain the ground-state fluid-phase orbitals.  If&
& this is a problem then it can easily be changed.')
if(mc_twist_av.and..not.xyzzyaabf1)call errstop_master('INIT_FREE_ORBS&
&','Monte Carlo twist averaging can only be performed for fluid phases&
&.')
inquire(file='correlation.data',exist=xyzzyaaar2)
if(.not.xyzzyaaar2.and.corr_heg_required)call errstop_master('INIT_FRE&
&E_ORBS','Cannot find correlation.data file, which is required.')
if(xyzzyaaar2.and.corr_heg_required)then
call open_units(xyzzyaaaa1,xyzzyaaac2)
open(unit=xyzzyaaaa1,file='correlation.data',status='old',iostat=xyzzy&
&aaac2)
if(xyzzyaaac2/=0)call errstop_master('INIT_FREE_ORBS','Problem opening&
& correlation.data .')
if(am_master)then
call wout('Reading free orbitals from correlation.data file.')
call wout()
endif
do
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaac2)char_80
if(trim(adjustl(char_80))=='START FREE_ORBS')exit
if(xyzzyaaac2<0)call errstop_master('INIT_FREE_ORBS','Cannot find "STA&
&RT FREE_ORBS" in correlation.data.')
if(xyzzyaaac2>0)call errstop_master('INIT_FREE_ORBS','Problem reading &
&correlation.data.')
enddo
read(xyzzyaaaa1,*,err=100,end=101)
read(xyzzyaaaa1,'(a)',err=100,end=101)title
title=trim(adjustl(title))
if(am_master)call wout(' Title                  : '//trim(adjustl(titl&
&e)))
mainloop: do
read(xyzzyaaaa1,'(a)',err=100,end=101)term_label
term_label=adjustl(term_label)
select case(trim(term_label))
case('START PAIRING')
if(.not.xyzzyaabg1)call errstop_master('INIT_FREE_ORBS','Found pairing&
& term in correlation.data but no pairing in use.')
if(xyzzyaadj1.or.xyzzyaadi1.or.xyzzyaadk1.or.xyzzyaadl1)call errstop_m&
&aster('INIT_FREE_ORBS','Two PAIRING sets in correlation.data ?')
if(am_master)call wout(' Pairing orbitals:')
read(xyzzyaaaa1,*,err=100,end=101)
read(xyzzyaaaa1,*,err=100,end=101)xyzzyaabu1
if(xyzzyaabu1/=0)call errstop_master('INIT_FREE_ORBS','Spin dependence&
& for pairing parameters may only be 0 currently.')
select case(xyzzyaabu1)
case(0)
xyzzyaabv1=1
end select
allocate(xyzzyaada1(xyzzyaabv1),stat=xyzzyaaae2)
call check_alloc(xyzzyaaae2,'INIT_FREE_ORBS','pair_norm')
if(am_master)call wout('  Particle-pair types    : '//trim(i2s(xyzzyaa&
&bv1)))
pairing_loop: do
read(xyzzyaaaa1,'(a)',err=100,end=101)term_label
term_label=adjustl(term_label)
select case(trim(term_label))
case('START GAUSSIAN TERM')
if(xyzzyaadj1)call errstop_master('INIT_FREE_ORBS','More than one Gaus&
&sian term found.')
xyzzyaadj1=.true.
if(am_master)call wout('  Gaussian term:')
read(xyzzyaaaa1,*,err=100,end=101)
read(xyzzyaaaa1,*,err=100,end=101)xyzzyaabq1
read(xyzzyaaaa1,*,err=100,end=101)
read(xyzzyaaaa1,*,err=100,end=101)xyzzyaadm1
allocate(xyzzyaacu1(real1_complex2,xyzzyaabq1,xyzzyaabv1),xyzzyaacc1(x&
&yzzyaabq1,xyzzyaabv1),xyzzyaacv1(xyzzyaabq1,xyzzyaabv1),xyzzyaacd1(xy&
&zzyaabq1,xyzzyaabv1),stat=xyzzyaaae2)
call check_alloc(xyzzyaaae2,'INIT_FREE_ORBS','pair_gauss_*')
xyzzyaacu1=0.d0
xyzzyaacc1=1
xyzzyaacv1=1.d0
xyzzyaacd1=1
if(xyzzyaadm1)then
call xyzzyaaga1
xyzzyaaco1=1.d0
xyzzyaabt1=1
if(am_master)call wout('   Number of parameters  : '//trim(i2s(xyzzyaa&
&bv1)))
read(xyzzyaaaa1,*,err=100,end=101)
read(xyzzyaaaa1,*,iostat=xyzzyaaac2)xyzzyaaco1,xyzzyaabt1
tmpr=r2s2(xyzzyaaco1,'(f18.12)')
if(xyzzyaaac2/=0)then
backspace xyzzyaaaa1
if(am_master)then
if(xyzzyaabt1==1)then
call wout('   Rex     (optimizable) : '//trim(tmpr)//' (default)')
else
call wout('   Rex           (fixed) : '//trim(tmpr)//' (default)')
endif
endif
elseif(am_master)then
if(xyzzyaabt1==1)then
call wout('   Rex     (optimizable) : '//trim(tmpr))
else
call wout('   Rex           (fixed) : '//trim(tmpr))
endif
endif
if(xyzzyaaco1<1.1d-8)call errstop_master('INIT_FREE_ORBS','Rex paramet&
&er too small!')
xyzzyaaaq2=1.d0/xyzzyaaco1
xyzzyaaap2=1.d0/sqrt(sum(xyzzyaacz1(:)**2))
do xyzzyaaab2=1,xyzzyaabv1
xyzzyaacu1(1,:,xyzzyaaab2)=xyzzyaacz1(:)*xyzzyaaap2
if(xyzzyaabl1)xyzzyaacu1(2,:,xyzzyaaab2)=0.d0
xyzzyaacv1(:,xyzzyaaab2)=(xyzzyaaaq2*xyzzyaacy1(:))**2
enddo
else
if(am_master)call wout('  Number of parameters  : '//trim(i2s(xyzzyaab&
&v1*xyzzyaabq1*(1+real1_complex2))))
read(xyzzyaaaa1,*,err=100,end=101)
read_pair_gauss: do xyzzyaaab2=1,xyzzyaabv1
do xyzzyaaaa2=1,xyzzyaabq1
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaac2)char_80
if(xyzzyaaac2==0)then
read(char_80,*,iostat=xyzzyaaac2)xyzzyaacu1(1:real1_complex2,xyzzyaaaa&
&2,xyzzyaaab2),xyzzyaacc1(xyzzyaaaa2,xyzzyaaab2)
if(complex_wf.and.xyzzyaaac2/=0)then
xyzzyaacu1(2,xyzzyaaaa2,xyzzyaaab2)=0.d0
read(char_80,*,iostat=xyzzyaaac2)xyzzyaacu1(1,xyzzyaaaa2,xyzzyaaab2),x&
&yzzyaacc1(xyzzyaaaa2,xyzzyaaab2)
endif
endif
if(xyzzyaaac2/=0)then
backspace xyzzyaaaa1
if(am_master)call wout('   Unspecified params    : Set to default (g=0&
&, w=1)')
goto 501
endif
tmpr=r2s2(xyzzyaacu1(1,xyzzyaaaa2,xyzzyaaab2),'(f18.12)')
if(xyzzyaabl1)tmpr=trim(tmpr)//','//trim(r2s2(xyzzyaacu1(2,xyzzyaaaa2,&
&xyzzyaaab2),'(f18.12)'))
if(am_master)then
if(xyzzyaacc1(xyzzyaaaa2,xyzzyaaab2)==1)then
call wout('   g_'//trim(i2s(xyzzyaaaa2))//','//trim(i2s(xyzzyaaab2))//&
&'     (optimizable) : '//trim(tmpr))
else
call wout('   g_'//trim(i2s(xyzzyaaaa2))//','//trim(i2s(xyzzyaaab2))//&
&'           (fixed) : '//trim(tmpr))
endif
endif
read(xyzzyaaaa1,*,iostat=xyzzyaaac2)xyzzyaacv1(xyzzyaaaa2,xyzzyaaab2),&
&xyzzyaacd1(xyzzyaaaa2,xyzzyaaab2)
if(xyzzyaaac2/=0)then
backspace xyzzyaaaa1
if(am_master)call wout('   Unspecified params    : Set to default (c=0&
&, w=1)')
goto 501
endif
tmpr=r2s2(xyzzyaacv1(xyzzyaaaa2,xyzzyaaab2),'(f18.12)')
if(am_master)then
if(xyzzyaacd1(xyzzyaaaa2,xyzzyaaab2)==1)then
call wout('   w_'//trim(i2s(xyzzyaaaa2))//','//trim(i2s(xyzzyaaab2))//&
&'     (optimizable) : '//trim(tmpr))
else
call wout('   w_'//trim(i2s(xyzzyaaaa2))//','//trim(i2s(xyzzyaaab2))//&
&'           (fixed) : '//trim(tmpr))
endif
endif
enddo
enddo read_pair_gauss
501    continue
endif
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaac2)char_80
if(xyzzyaaac2/=0)call errstop_master('INIT_FREE_ORBS','String "END GAU&
&SSIAN TERM" not found.')
if(trim(adjustl(char_80))/='END GAUSSIAN TERM')call errstop_master('IN&
&IT_FREE_ORBS','String "END GAUSSIAN TERM" not found.')
case('START PLANE-WAVE TERM')
if(xyzzyaadi1)call errstop_master('INIT_FREE_ORBS','More than one plan&
&e-wave term found.')
xyzzyaadi1=.true.
if(am_master)call wout('  Plane-wave term:')
read(xyzzyaaaa1,*,err=100,end=101)
read(xyzzyaaaa1,*,err=100,end=101)xyzzyaabr1
allocate(xyzzyaacw1(real1_complex2,xyzzyaabr1,xyzzyaabv1),xyzzyaacf1(x&
&yzzyaabr1,xyzzyaabv1),stat=xyzzyaaae2)
call check_alloc(xyzzyaaae2,'INIT_FREE_ORBS','pair_pw_c')
xyzzyaacw1=0.d0
xyzzyaacf1=1
if(am_master)call wout('   Number of parameters  : '//trim(i2s(xyzzyaa&
&bv1*xyzzyaabr1*real1_complex2)))
read(xyzzyaaaa1,*,err=100,end=101)
read_pair_pw: do xyzzyaaab2=1,xyzzyaabv1
do xyzzyaaaa2=1,xyzzyaabr1
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaac2)char_80
if(xyzzyaaac2==0)then
read(char_80,*,iostat=xyzzyaaac2)xyzzyaacw1(1:real1_complex2,xyzzyaaaa&
&2,xyzzyaaab2),xyzzyaacf1(xyzzyaaaa2,xyzzyaaab2)
if(complex_wf.and.xyzzyaaac2/=0)then
xyzzyaacw1(2,xyzzyaaaa2,xyzzyaaab2)=0.d0
read(char_80,*,iostat=xyzzyaaac2)xyzzyaacw1(1,xyzzyaaaa2,xyzzyaaab2),x&
&yzzyaacf1(xyzzyaaaa2,xyzzyaaab2)
endif
endif
if(xyzzyaaac2/=0)then
backspace xyzzyaaaa1
if(am_master)call wout('   Unspecified params    : Set to default')
goto 502
endif
xyzzyaabk1=.false.
tmpr=r2s2(xyzzyaacw1(1,xyzzyaaaa2,xyzzyaaab2),'(f18.12)')
if(xyzzyaabl1)tmpr=trim(tmpr)//','//trim(r2s2(xyzzyaacw1(2,xyzzyaaaa2,&
&xyzzyaaab2),'(f18.12)'))
if(am_master)then
if(xyzzyaacf1(xyzzyaaaa2,xyzzyaaab2)==1)then
call wout('   c_'//trim(i2s(xyzzyaaaa2))//','//trim(i2s(xyzzyaaab2))//&
&'   (optimizable) : '//trim(tmpr))
else
call wout('   c_'//trim(i2s(xyzzyaaaa2))//','//trim(i2s(xyzzyaaab2))//&
&'         (fixed) : '//trim(tmpr))
endif
endif
enddo
enddo read_pair_pw
502   continue
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaac2)char_80
if(xyzzyaaac2/=0)call errstop_master('INIT_FREE_ORBS','String "END PLA&
&NE-WAVE TERM" not found.')
if(trim(adjustl(char_80))/='END PLANE-WAVE TERM')call errstop_master('&
&INIT_FREE_ORBS','String "END PLANE-WAVE TERM" not found.')
case('START POLYNOMIAL TERM')
if(xyzzyaadk1)call errstop_master('INIT_FREE_ORBS','More than one poly&
&nomial term found.')
xyzzyaadk1=.true.
if(am_master)call wout('  Polynomial term:')
read(xyzzyaaaa1,*,err=100,end=101)
read(xyzzyaaaa1,*,err=100,end=101)xyzzyaabx1
if(xyzzyaabx1<1)call errstop_master('INIT_FREE_ORBS','The order N_p of&
& polynomial pairing orbitals cannot be less than 1.')
read(xyzzyaaaa1,*,err=100,end=101)
read(xyzzyaaaa1,*,err=100,end=101)xyzzyaaby1
if(xyzzyaaby1<0)call errstop_master('INIT_FREE_ORBS','The truncation o&
&rder C_p of polynomial pairing orbitals cannot be less than 0.')
allocate(xyzzyaadb1(real1_complex2,0:xyzzyaabx1,xyzzyaabv1),xyzzyaach1&
&(0:xyzzyaabx1,xyzzyaabv1),stat=xyzzyaaae2)
call check_alloc(xyzzyaaae2,'INIT_FREE_ORBS','poly')
xyzzyaadb1=0.d0
xyzzyaach1=1
xyzzyaach1(1,:)=-1
if(am_master)then
call wout('   Order of polynomials  : '//trim(i2s(xyzzyaabx1)))
call wout('   Truncation, C_p       : '//trim(i2s(xyzzyaaby1)))
call wout('   Number of parameters  : '//trim(i2s(1+xyzzyaabv1*(xyzzya&
&abx1+1)*real1_complex2)))
endif
read(xyzzyaaaa1,*,err=100,end=101)
read(xyzzyaaaa1,*,err=100,end=101)xyzzyaacp1,xyzzyaabz1
if(xyzzyaacp1<0.d0)call errstop_master('INIT_FREE_ORBS','The cutoff L_&
&p of polynomial pairing orbitals must be positive (set to zero for de&
&fault).')
if(am_master)then
if(xyzzyaacp1==0.d0)then
tmpr="0.99 * L_WS [not yet known, see below]"
else
tmpr=r2s2(xyzzyaacp1,'(f18.12)')
endif
endif
select case(xyzzyaabz1)
case(0)
if(am_master)call wout('   Cutoff (au)   (fixed) : '//trim(tmpr))
case(1)
xyzzyaadn1=.false.
if(am_master)call wout('   Cutoff (au)     (opt) : '//trim(tmpr))
case(2)
xyzzyaabz1=1
xyzzyaadn1=.true.
if(am_master)call wout('   Cutoff (au) (lim opt) : '//trim(tmpr))
case default
call errstop_master('INIT_FREE_ORBS','Unknown L_p_optable value')
end select
read_pair_poly: do xyzzyaaab2=1,xyzzyaabv1
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaac2)char_80
if(xyzzyaaac2==0)then
read(char_80,*,iostat=xyzzyaaac2)xyzzyaadb1(1:real1_complex2,0,xyzzyaa&
&ab2),xyzzyaach1(0,xyzzyaaab2)
if(complex_wf.and.xyzzyaaac2/=0)then
xyzzyaadb1(2,0,xyzzyaaab2)=0.d0
read(char_80,*,iostat=xyzzyaaac2)xyzzyaadb1(1,0,xyzzyaaab2),xyzzyaach1&
&(0,xyzzyaaab2)
endif
endif
if(xyzzyaaac2/=0)then
backspace xyzzyaaaa1
if(am_master)call wout('   Unspecified params    : Set to default (alp&
&ha=0, optable=1)')
goto 503
endif
if(xyzzyaacp1==0.d0)then
xyzzyaadb1(:,1,xyzzyaaab2)=0.d0
else
xyzzyaadb1(:,1,xyzzyaaab2)=xyzzyaadb1(:,0,xyzzyaaab2)*real(xyzzyaaby1,&
&dp)/xyzzyaacp1
endif
if(am_master)then
tmpr=r2s2(xyzzyaadb1(1,0,xyzzyaaab2),'(f18.12)')
tmpr2=''
if(xyzzyaabl1)tmpr2=','//trim(r2s2(xyzzyaadb1(2,0,xyzzyaaab2),'(f18.12&
&)'))
if(xyzzyaach1(0,xyzzyaaab2)==1)then
call wout('   alpha_0,'//trim(i2s(xyzzyaaab2))//'       (opt) : '//tri&
&m(tmpr)//trim(tmpr2))
else
call wout('   alpha_0,'//trim(i2s(xyzzyaaab2))//'     (fixed) : '//tri&
&m(tmpr)//trim(tmpr2))
endif
if(xyzzyaacp1==0.d0.and.xyzzyaaby1/=0)then
tmpr='[to be calculated when L_p known]'
tmpr2=''
else
tmpr=r2s2(xyzzyaadb1(1,1,xyzzyaaab2),'(f18.12)')
tmpr2=''
if(xyzzyaabl1)tmpr2=','//trim(r2s2(xyzzyaadb1(2,1,xyzzyaaab2),'(f18.12&
&)'))
endif
call wout('   alpha_1,'//trim(i2s(xyzzyaaab2))//'    (constr) : '//tri&
&m(tmpr)//trim(tmpr2))
endif
do xyzzyaaaa2=2,xyzzyaabx1
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaac2)char_80
if(xyzzyaaac2==0)then
read(char_80,*,iostat=xyzzyaaac2)xyzzyaadb1(1:real1_complex2,xyzzyaaaa&
&2,xyzzyaaab2),xyzzyaach1(xyzzyaaaa2,xyzzyaaab2)
if(complex_wf.and.xyzzyaaac2/=0)then
xyzzyaadb1(2,xyzzyaaaa2,xyzzyaaab2)=0.d0
read(char_80,*,iostat=xyzzyaaac2)xyzzyaadb1(1,xyzzyaaaa2,xyzzyaaab2),x&
&yzzyaach1(xyzzyaaaa2,xyzzyaaab2)
endif
endif
if(xyzzyaaac2/=0)then
backspace xyzzyaaaa1
if(am_master)call wout('   Unspecified params    : Set to default (alp&
&ha=0, optable=1)')
exit read_pair_poly
endif
if(am_master)then
tmpr=r2s2(xyzzyaadb1(1,xyzzyaaaa2,xyzzyaaab2),'(f18.12)')
tmpr2=''
if(xyzzyaabl1)tmpr2=','//trim(r2s2(xyzzyaadb1(2,xyzzyaaaa2,xyzzyaaab2)&
&,'(f18.12)'))
if(xyzzyaach1(xyzzyaaaa2,xyzzyaaab2)==1)then
call wout('   alpha_'//trim(i2s(xyzzyaaaa2))//','//trim(i2s(xyzzyaaab2&
&))//'       (opt) : '//trim(tmpr)//trim(tmpr2))
else
call wout('   alpha_'//trim(i2s(xyzzyaaaa2))//','//trim(i2s(xyzzyaaab2&
&))//'     (fixed) : '//trim(tmpr)//trim(tmpr2))
endif
endif
enddo
enddo read_pair_poly
503   continue
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaac2)char_80
if(xyzzyaaac2/=0)call errstop_master('INIT_FREE_ORBS','String "END POL&
&YNOMIAL TERM" not found.')
if(trim(adjustl(char_80))/='END POLYNOMIAL TERM')call errstop_master('&
&INIT_FREE_ORBS','String "END POLYNOMIAL TERM" not found.')
case('START SLATER TERM')
if(xyzzyaadl1)call errstop_master('INIT_FREE_ORBS','More than one Slat&
&er term found.')
xyzzyaadl1=.true.
if(am_master)call wout('  Slater term:')
read(xyzzyaaaa1,*,err=100,end=101)
read(xyzzyaaaa1,*,err=100,end=101)xyzzyaaca1
read(xyzzyaaaa1,*,err=100,end=101)
read(xyzzyaaaa1,*,err=100,end=101)xyzzyaacb1
if(xyzzyaaca1+xyzzyaacb1<1)call errstop_master('INIT_FREE_ORBS','Slate&
&r-type pairing orbital is empty!')
allocate(xyzzyaadc1(real1_complex2,xyzzyaaca1,xyzzyaabv1),xyzzyaaci1(x&
&yzzyaaca1,xyzzyaabv1),              xyzzyaadd1(xyzzyaaca1,xyzzyaabv1)&
&,                      xyzzyaade1(xyzzyaaca1,xyzzyaabv1),            &
&          xyzzyaack1(xyzzyaaca1,xyzzyaabv1),              xyzzyaacj1(&
&xyzzyaaca1,xyzzyaabv1),              xyzzyaadf1(real1_complex2,3,xyzz&
&yaacb1,xyzzyaabv1),     xyzzyaacl1(3,xyzzyaacb1,xyzzyaabv1),         &
&   xyzzyaadg1(xyzzyaacb1,xyzzyaabv1),                      xyzzyaacm1&
&(xyzzyaacb1,xyzzyaabv1),              xyzzyaadh1(xyzzyaacb1,xyzzyaabv&
&1),                      xyzzyaacn1(xyzzyaacb1,xyzzyaabv1),          &
&    stat=xyzzyaaae2)
call check_alloc(xyzzyaaae2,'INIT_FREE_ORBS','Slater')
xyzzyaadc1=0.d0
xyzzyaadc1(1,:,:)=1.d0
xyzzyaaci1=1
xyzzyaade1=1.d0
xyzzyaack1=1
xyzzyaadf1=0.d0
xyzzyaacl1=1
xyzzyaadh1=1.d0
xyzzyaacn1=1
xyzzyaacj1=1
xyzzyaacm1=1
do xyzzyaaaa2=1,xyzzyaaca1
xyzzyaadd1(xyzzyaaaa2,:)=real(xyzzyaaaa2,dp)
enddo
do xyzzyaaaa2=1,xyzzyaacb1
xyzzyaadg1(xyzzyaaaa2,:)=real(xyzzyaaaa2,dp)
enddo
if(am_master)then
call wout('   Order of S term       : '//trim(i2s(xyzzyaaca1)))
call wout('   Order of P term       : '//trim(i2s(xyzzyaacb1)))
call wout('   Number of parameters  : '//trim(i2s(3*xyzzyaaca1+(dimens&
&ionality+1)*xyzzyaacb1*real1_complex2)))
endif
read(xyzzyaaaa1,*,err=100,end=101)
read_pair_slater: do xyzzyaaab2=1,xyzzyaabv1
do xyzzyaaaa2=1,xyzzyaaca1
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaac2)char_80
if(xyzzyaaac2==0)then
read(char_80,*,iostat=xyzzyaaac2)xyzzyaadc1(1:real1_complex2,xyzzyaaaa&
&2,xyzzyaaab2),xyzzyaaci1(xyzzyaaaa2,xyzzyaaab2)
if(complex_wf.and.xyzzyaaac2/=0)then
xyzzyaadc1(2,xyzzyaaaa2,xyzzyaaab2)=0.d0
read(char_80,*,iostat=xyzzyaaac2)xyzzyaadc1(1,xyzzyaaaa2,xyzzyaaab2),x&
&yzzyaaci1(xyzzyaaaa2,xyzzyaaab2)
endif
endif
if(xyzzyaaac2/=0)then
backspace xyzzyaaaa1
if(am_master)call wout('   Unspecified params    : Set to default (c=1&
&, a>=1, b=1, d=1, optable=1)')
goto 504
endif
if(am_master)then
tmpr=r2s2(xyzzyaadc1(1,xyzzyaaaa2,xyzzyaaab2),'(f18.12)')
tmpr2=''
if(xyzzyaabl1)tmpr2=','//trim(r2s2(xyzzyaadc1(2,xyzzyaaaa2,xyzzyaaab2)&
&,'(f18.12)'))
if(xyzzyaaci1(xyzzyaaaa2,xyzzyaaab2)==1)then
call wout('   S_c_'//trim(i2s(xyzzyaaaa2))//','//trim(i2s(xyzzyaaab2))&
&//'         (opt) : '//trim(tmpr)//trim(tmpr2))
else
call wout('   S_c_'//trim(i2s(xyzzyaaaa2))//','//trim(i2s(xyzzyaaab2))&
&//'       (fixed) : '//trim(tmpr)//trim(tmpr2))
endif
endif
read(xyzzyaaaa1,*,iostat=xyzzyaaac2)xyzzyaadd1(xyzzyaaaa2,xyzzyaaab2),&
&xyzzyaacj1(xyzzyaaaa2,xyzzyaaab2)
if(xyzzyaaac2/=0)then
backspace xyzzyaaaa1
if(am_master)call wout('   Unspecified params    : Set to default (c=1&
&, a>=1, b=1, d=1, optable=1)')
goto 504
endif
if(am_master)then
tmpr=r2s2(xyzzyaadd1(xyzzyaaaa2,xyzzyaaab2),'(f18.12)')
if(xyzzyaacj1(xyzzyaaaa2,xyzzyaaab2)==1)then
call wout('   S_a_'//trim(i2s(xyzzyaaaa2))//','//trim(i2s(xyzzyaaab2))&
&//'         (opt) : '//trim(tmpr))
else
call wout('   S_a_'//trim(i2s(xyzzyaaaa2))//','//trim(i2s(xyzzyaaab2))&
&//'       (fixed) : '//trim(tmpr))
endif
endif
read(xyzzyaaaa1,*,iostat=xyzzyaaac2)xyzzyaade1(xyzzyaaaa2,xyzzyaaab2),&
&xyzzyaack1(xyzzyaaaa2,xyzzyaaab2)
if(xyzzyaaac2/=0)then
backspace xyzzyaaaa1
if(am_master)call wout('   Unspecified params    : Set to default (c=1&
&, a>=1, b=1, d=1, optable=1)')
goto 504
endif
if(am_master)then
tmpr=r2s2(xyzzyaade1(xyzzyaaaa2,xyzzyaaab2),'(f18.12)')
if(xyzzyaack1(xyzzyaaaa2,xyzzyaaab2)==1)then
call wout('   S_b_'//trim(i2s(xyzzyaaaa2))//','//trim(i2s(xyzzyaaab2))&
&//'         (opt) : '//trim(tmpr))
else
call wout('   S_b_'//trim(i2s(xyzzyaaaa2))//','//trim(i2s(xyzzyaaab2))&
&//'       (fixed) : '//trim(tmpr))
endif
endif
enddo
do xyzzyaaaa2=1,xyzzyaacb1
do d=1,dimensionality
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaac2)char_80
if(xyzzyaaac2==0)then
read(char_80,*,iostat=xyzzyaaac2)xyzzyaadf1(1:real1_complex2,d,xyzzyaa&
&aa2,xyzzyaaab2),xyzzyaacl1(d,xyzzyaaaa2,xyzzyaaab2)
if(complex_wf.and.xyzzyaaac2/=0)then
xyzzyaadf1(2,d,xyzzyaaaa2,xyzzyaaab2)=0.d0
read(char_80,*,iostat=xyzzyaaac2)xyzzyaadf1(1,d,xyzzyaaaa2,xyzzyaaab2)&
&,xyzzyaacl1(d,xyzzyaaaa2,xyzzyaaab2)
endif
endif
if(xyzzyaaac2/=0)then
backspace xyzzyaaaa1
if(am_master)call wout('   Unspecified params    : Set to default (c=1&
&, a>=1, b=1, optable=1)')
goto 504
endif
if(am_master)then
tmpr=r2s2(xyzzyaadf1(1,d,xyzzyaaaa2,xyzzyaaab2),'(f18.12)')
tmpr2=''
if(xyzzyaabl1)tmpr2=','//trim(r2s2(xyzzyaadf1(2,d,xyzzyaaaa2,xyzzyaaab&
&2),'(f18.12)'))
if(xyzzyaacl1(d,xyzzyaaaa2,xyzzyaaab2)==1)then
call wout('   P_c_'//trim(i2s(d))//','//trim(i2s(xyzzyaaaa2))//','//tr&
&im(i2s(xyzzyaaab2))//'       (opt) : '//trim(tmpr)//trim(tmpr2))
else
call wout('   P_c_'//trim(i2s(d))//','//trim(i2s(xyzzyaaaa2))//','//tr&
&im(i2s(xyzzyaaab2))//'     (fixed) : '//trim(tmpr)//trim(tmpr2))
endif
endif
enddo
read(xyzzyaaaa1,*,iostat=xyzzyaaac2)xyzzyaadg1(xyzzyaaaa2,xyzzyaaab2),&
&xyzzyaacm1(xyzzyaaaa2,xyzzyaaab2)
if(xyzzyaaac2/=0)then
backspace xyzzyaaaa1
if(am_master)call wout('   Unspecified params    : Set to default (c=1&
&, a>=1, b=1, optable=1)')
goto 504
endif
if(am_master)then
tmpr=r2s2(xyzzyaadg1(xyzzyaaaa2,xyzzyaaab2),'(f18.12)')
if(xyzzyaacm1(xyzzyaaaa2,xyzzyaaab2)==1)then
call wout('   P_a_'//trim(i2s(xyzzyaaaa2))//','//trim(i2s(xyzzyaaab2))&
&//'         (opt) : '//trim(tmpr))
else
call wout('   P_a_'//trim(i2s(xyzzyaaaa2))//','//trim(i2s(xyzzyaaab2))&
&//'       (fixed) : '//trim(tmpr))
endif
endif
read(xyzzyaaaa1,*,iostat=xyzzyaaac2)xyzzyaadh1(xyzzyaaaa2,xyzzyaaab2),&
&xyzzyaacn1(xyzzyaaaa2,xyzzyaaab2)
if(xyzzyaaac2/=0)then
backspace xyzzyaaaa1
if(am_master)call wout('   Unspecified params    : Set to default (c=1&
&, a>=1, b=1, optable=1)')
goto 504
endif
if(am_master)then
tmpr=r2s2(xyzzyaadh1(xyzzyaaaa2,xyzzyaaab2),'(f18.12)')
if(xyzzyaacn1(xyzzyaaaa2,xyzzyaaab2)==1)then
call wout('   P_b_'//trim(i2s(xyzzyaaaa2))//','//trim(i2s(xyzzyaaab2))&
&//'         (opt) : '//trim(tmpr))
else
call wout('   P_b_'//trim(i2s(xyzzyaaaa2))//','//trim(i2s(xyzzyaaab2))&
&//'       (fixed) : '//trim(tmpr))
endif
endif
enddo
enddo read_pair_slater
504  continue
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaac2)char_80
if(xyzzyaaac2/=0)call errstop_master('INIT_FREE_ORBS','String "END SLA&
&TER TERM" not found.')
if(trim(adjustl(char_80))/='END SLATER TERM')call errstop_master('INIT&
&_FREE_ORBS','String "END SLATER TERM" not found.')
case('END PAIRING')
if(.not.(xyzzyaadj1.or.xyzzyaadi1.or.xyzzyaadk1.or.xyzzyaadl1))call er&
&rstop_master('INIT_FREE_ORBS','No known pairing term found in the PAI&
&RING block!')
do xyzzyaaab2=1,xyzzyaabv1
xyzzyaaap2=0.d0
if(xyzzyaadj1)xyzzyaaap2=xyzzyaaap2+sum(xyzzyaacu1(:,:,xyzzyaaab2)**2)
if(xyzzyaadi1)xyzzyaaap2=xyzzyaaap2+sum(xyzzyaacw1(:,:,xyzzyaaab2)**2)
if(xyzzyaadk1)xyzzyaaap2=xyzzyaaap2+sum(xyzzyaadb1(:,:,xyzzyaaab2)**2)
if(xyzzyaadl1)xyzzyaaap2=xyzzyaaap2+sum(xyzzyaadc1(:,:,xyzzyaaab2)**2)&
&+sum(xyzzyaadf1(:,:,:,xyzzyaaab2)**2)
if(xyzzyaaap2<=0.d0)call errstop_master('INIT_FREE_ORBS','Pairing orbi&
&tal is zero.')
xyzzyaada1(xyzzyaaab2)=1.d0/sqrt(xyzzyaaap2)
enddo
goto 505
case default
call errstop_master('INIT_FREE_ORBS','Unknown pairing term found in PA&
&IRING block.')
end select
enddo pairing_loop
505 continue
case('START WIGNER CRYSTAL')
if(.not.xyzzyaabh1)call errstop_master('INIT_FREE_ORBS','Found Wigner &
&crystal term in correlation.data but no WC in use.')
if(type_wc/=0)call errstop_master('INIT_FREE_ORBS','Two WIGNER CRYSTAL&
& sets in correlation.data?')
if(am_master)call wout(' Wigner Crystal orbitals:')
read(xyzzyaaaa1,*,err=100,end=101)
read(xyzzyaaaa1,*,err=100,end=101)char_80
char_80=adjustl(char_80)
call capitalize(char_80)
if(trim(char_80)=='GAUSSIAN')then
type_wc=1
elseif(trim(char_80)=='PLANE-WAVE')then
type_wc=2
elseif(trim(char_80)=='GAUSSIAN+PLANE-WAVE'.or.trim(char_80)=='HYBRID'&
&)then
type_wc=3
else
read(char_80,*,err=100,end=101)type_wc
endif
if(type_wc<1.or.type_wc>3)call errstop_master('INIT_FREE_ORBS','Type o&
&f WC orbitals must be 1, 2 or 3.')
if(am_master)then
select case(type_wc)
case(1)
call wout('  Type of orbitals      : Gaussian')
case(2)
call wout('  Type of orbitals      : Plane-wave expansion')
case(3)
call wout('  Type of orbitals      : Hybrid')
end select
endif
read(xyzzyaaaa1,*,err=100,end=101)
read(xyzzyaaaa1,*,err=100,end=101)spin_dep_wc
if(spin_dep_wc<0.or.spin_dep_wc>levels_ssingles)call errstop_master('I&
&NIT_FREE_ORBS','Spin dependence of WC parameters must be 0 -- '//trim&
&(i2s(levels_ssingles))//'.')
xyzzyaabw1=no_ssingles(spin_dep_wc)
if(am_master)then
call wout('  Spin types            : '//trim(i2s(xyzzyaabw1)))
endif
if(type_wc==1.or.type_wc==3)then
allocate(wc_gauss_exp(xyzzyaabw1),xyzzyaace1(xyzzyaabw1),stat=xyzzyaaa&
&e2)
call check_alloc(xyzzyaaae2,'INIT_FREE_ORBS','wc_gauss_exp')
wc_gauss_exp=1.d0
xyzzyaace1=1
if(am_master)then
if(type_wc==1)then
call wout('  Number of parameters  : '//trim(i2s(xyzzyaabw1)))
else
call wout('  No. of Gaussian params: '//trim(i2s(xyzzyaabw1)))
endif
endif
read(xyzzyaaaa1,*,err=100,end=101)
read_wc_gauss: do xyzzyaaab2=1,xyzzyaabw1
read(xyzzyaaaa1,*,iostat=xyzzyaaac2)wc_gauss_exp(xyzzyaaab2),xyzzyaace&
&1(xyzzyaaab2)
if(xyzzyaaac2/=0)then
backspace xyzzyaaaa1
if(am_master)call wout('  Unspecified params    : Set to default')
goto 506
endif
if(wc_gauss_exp(xyzzyaaab2)<0.d0)call errstop_master('INIT_FREE_ORBS',&
&'Gaussian exponents for Wigner Crystal must be >= 0 .')
tmpr=r2s2(wc_gauss_exp(xyzzyaaab2),'(f18.12)')
if(am_master)then
if(xyzzyaace1(xyzzyaaab2)==1)then
call wout('  w_'//trim(i2s(xyzzyaaab2))//'     (optimizable) : '//trim&
&(tmpr))
else
call wout('  w_'//trim(i2s(xyzzyaaab2))//'           (fixed) : '//trim&
&(tmpr))
endif
endif
enddo read_wc_gauss
506  continue
endif
if(type_wc==2.or.type_wc==3)then
read(xyzzyaaaa1,*,err=100,end=101)
read(xyzzyaaaa1,*,err=100,end=101)xyzzyaabs1
read(xyzzyaaaa1,*,err=100,end=101)
allocate(xyzzyaacx1(xyzzyaabs1,xyzzyaabw1),xyzzyaacg1(xyzzyaabs1,xyzzy&
&aabw1),stat=xyzzyaaae2)
call check_alloc(xyzzyaaae2,'INIT_FREE_ORBS','wc_pw_c')
xyzzyaacx1=0.d0
xyzzyaacg1=1
if(type_wc==2)then
xyzzyaacx1(1,:)=1.d0
xyzzyaacg1(1,:)=0
endif
if(am_master)then
if(type_wc==2)then
call wout('  Number of parameters  : '//trim(i2s(xyzzyaabw1*xyzzyaabs1&
&)))
else
call wout('  Number of PW params   : '//trim(i2s(xyzzyaabw1*xyzzyaabs1&
&)))
endif
endif
read_wc_pw: do xyzzyaaab2=1,xyzzyaabw1
do xyzzyaaaa2=1,xyzzyaabs1
read(xyzzyaaaa1,*,iostat=xyzzyaaac2)xyzzyaacx1(xyzzyaaaa2,xyzzyaaab2),&
&xyzzyaacg1(xyzzyaaaa2,xyzzyaaab2)
if(xyzzyaaac2/=0)then
backspace xyzzyaaaa1
if(am_master)call wout('  Unspecified params    : Set to default')
xyzzyaacx1(xyzzyaaaa2,xyzzyaaab2)=0.d0
xyzzyaacg1(xyzzyaaaa2,xyzzyaaab2)=1
if(type_wc==2.and.xyzzyaaaa2==1)then
xyzzyaacx1(1,xyzzyaaab2)=1.d0
xyzzyaacg1(1,xyzzyaaab2)=0
endif
goto 507
endif
tmpr=r2s2(xyzzyaacx1(xyzzyaaaa2,xyzzyaaab2),'(f18.12)')
if(am_master)then
if(xyzzyaacg1(xyzzyaaaa2,xyzzyaaab2)==1)then
call wout('  c_'//trim(i2s(xyzzyaaaa2))//','//trim(i2s(xyzzyaaab2))//'&
&     (optimizable) : '//trim(tmpr))
else
call wout('  c_'//trim(i2s(xyzzyaaaa2))//','//trim(i2s(xyzzyaaab2))//'&
&           (fixed) : '//trim(tmpr))
endif
endif
enddo
enddo read_wc_pw
507  continue
if(type_wc==2)xyzzyaacg1(1,1:xyzzyaabw1)=0
endif
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaac2)char_80
if(xyzzyaaac2/=0)call errstop_master('INIT_FREE_ORBS','String "END WIG&
&NER CRYSTAL" not found.')
if(trim(adjustl(char_80))/='END WIGNER CRYSTAL')call errstop_master('I&
&NIT_FREE_ORBS','String "END WIGNER CRYSTAL" not found.')
case('START SPIN DENSITY WAVE')
if(.not.xyzzyaabj1)call errstop_master('INIT_FREE_ORBS','Found SDW sec&
&tion in correlation.data but no SDW orbitals in use.')
if(xyzzyaaat2)call errstop_master('INIT_FREE_ORBS','Two Spin Density W&
&ave sections in correlation.data ?')
xyzzyaaat2=.true.
allocate(sdw_kvec(3,nele(1)),sdw_orb_theta(nele(1)),stat=xyzzyaaae2)
call check_alloc(xyzzyaaae2,'INIT_FREE_ORBS','k-vectors')
read(xyzzyaaaa1,*,err=100,end=101)
read(xyzzyaaaa1,*,err=100,end=101)xyzzyaaaf2
if(am_master)then
call wout('  Reading '//trim(i2s(xyzzyaaaf2))//' k-vectors.')
endif
read(xyzzyaaaa1,*,err=100,end=101)
xyzzyaaah2=0
do xyzzyaaaa2=1,xyzzyaaaf2
xyzzyaaah2=xyzzyaaah2+1
if(xyzzyaaah2>nele(1))call errstop_master('INIT_FREE_ORBS','Number of &
&SDW orbitals defined in correlation.data is incorrect')
read(xyzzyaaaa1,*,err=100,end=101)sdw_kvec(1:3,xyzzyaaah2)
read(xyzzyaaaa1,*,err=100,end=101)xyzzyaaao2,xyzzyaaag2
if(xyzzyaaag2==1)then
sdw_orb_theta(xyzzyaaah2)=xyzzyaaao2
elseif(xyzzyaaag2==2)then
xyzzyaaah2=xyzzyaaah2+1
if(xyzzyaaah2>nele(1))call errstop_master('INIT_FREE_ORBS','Number of &
&SDW orbitals defined in correlation.data is incorrect')
sdw_kvec(1:3,xyzzyaaah2)=sdw_kvec(1:3,xyzzyaaah2-1)
sdw_orb_theta(xyzzyaaah2-1)=xyzzyaaao2
sdw_orb_theta(xyzzyaaah2)=xyzzyaaao2+pi
else
call errstop_master('INIT_FREE_ORBS','Incorrect occupation number')
endif
enddo
if(xyzzyaaah2<nele(1))call errstop_master('INIT_FREE_ORBS','Number of &
&SDW orbitals defined in correlation.data is incorrect')
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaac2)char_80
if(xyzzyaaac2/=0)call errstop_master('INIT_FREE_ORBS','String "END SPI&
&N DENSITY WAVE" not found.')
if(trim(adjustl(char_80))/='END SPIN DENSITY WAVE')call errstop_master&
&('INIT_FREE_ORBS','String "END SPIN DENSITY WAVE" not found.')
case('END FREE_ORBS')
goto 601
case default
call errstop_master('INIT_FREE_ORBS','Unrecognized label '//trim(term_&
&label)//' reading FREE_ORBS term in correlation.data .')
end select
enddo mainloop
601 continue
close(xyzzyaaaa1)
if(am_master)then
call wout()
call wout('Finished reading free-particle block from correlation.data.&
&')
call wout()
endif
if(xyzzyaabg1.and..not.(xyzzyaadj1.or.xyzzyaadi1.or.xyzzyaadk1.or.xyzz&
&yaadl1))call errstop_master('INIT_FREE_ORBS','No (known) pairing para&
&meters found in correlation.data .')
if(xyzzyaabh1.and.type_wc==0)call errstop_master('INIT_FREE_ORBS','WC &
&parameters not found in correlation.data .')
endif
xyzzyaaab1=0
if(xyzzyaadi1)xyzzyaaab1=xyzzyaabr1
if(type_wc==2.or.type_wc==3)xyzzyaaab1=max(xyzzyaabs1,xyzzyaaab1)
if((xyzzyaaab1>0.or.xyzzyaabf1).and.periodicity==0.and.(any(nele>1) .o&
&r.complex_wf))call errstop_master('INIT_FREE_ORBS','Cannot have plane&
& waves in aperiodic systems (unless we defined a cell size, but we do&
& not.)')
if(xyzzyaabg1)then
xyzzyaaff1=0
xyzzyaafg1=0
allocate(xyzzyaaao1(nemax,nspin,ndet),stat=xyzzyaaae2)
call check_alloc(xyzzyaaae2,'INIT_FREE_ORBS','free_orbmap')
xyzzyaaao1=0
do xyzzyaaai2=1,nspin
do xyzzyaaad2=1,ndet
xyzzyaaaj2=-heg_orbtype(xyzzyaaai2,xyzzyaaad2)
xyzzyaaak2=0
if(xyzzyaaaj2>1)xyzzyaaak2=sum(nele(1:xyzzyaaaj2-1))
if(update_by_column(xyzzyaaai2,xyzzyaaad2))then
xyzzyaaal2=heg_nele(xyzzyaaaj2)
else
xyzzyaaal2=nele(xyzzyaaaj2)
endif
do xyzzyaaam2=1,xyzzyaaal2
xyzzyaaao1(xyzzyaaam2,xyzzyaaai2,xyzzyaaad2)=xyzzyaaak2+xyzzyaaam2
enddo
enddo
enddo
free_norb=netot
elseif(xyzzyaabf1)then
xyzzyaaff1=1
xyzzyaafg1=0
allocate(xyzzyaaao1(maxval(heg_nele),nspin,ndet),stat=xyzzyaaae2)
call check_alloc(xyzzyaaae2,'INIT_FREE_ORBS','free_orbmap')
xyzzyaaao1=0
do xyzzyaaai2=1,nspin
do xyzzyaaam2=1,heg_nele(xyzzyaaai2)
xyzzyaaao1(xyzzyaaam2,xyzzyaaai2,1)=xyzzyaaam2
enddo
do xyzzyaaad2=2,ndet
xyzzyaaao1(:,xyzzyaaai2,xyzzyaaad2)=xyzzyaaao1(:,xyzzyaaai2,1)
enddo
enddo
if(free_norb==0)free_norb=maxval(heg_nele)
elseif(xyzzyaabh1)then
xyzzyaaff1=0
xyzzyaafg1=0
allocate(xyzzyaaao1(maxval(heg_nele),nspin,ndet),stat=xyzzyaaae2)
call check_alloc(xyzzyaaae2,'INIT_FREE_ORBS','free_orbmap')
xyzzyaaao1=0
do xyzzyaaai2=1,nspin
do xyzzyaaad2=1,ndet
xyzzyaaan2=heg_slatt(xyzzyaaai2,xyzzyaaad2)
do xyzzyaaam2=1,heg_nele(xyzzyaaai2)
xyzzyaaao1(xyzzyaaam2,xyzzyaaai2,xyzzyaaad2)=(xyzzyaaan2-1)*maxval(heg&
&_nele)+xyzzyaaam2
enddo
enddo
enddo
free_norb=maxval(heg_slatt)*maxval(heg_nele)
elseif(xyzzyaabm1.or.xyzzyaabn1.or.xyzzyaabo1)then
xyzzyaaff1=0
xyzzyaafg1=0
allocate(xyzzyaaao1(1,nspin,1),stat=xyzzyaaae2)
call check_alloc(xyzzyaaae2,'INIT_FREE_ORBS','free_orbmap')
xyzzyaaao1=0
do xyzzyaaai2=1,nspin
xyzzyaaao1(1,xyzzyaaai2,1)=1
enddo
free_norb=1
endif
if(.not.xyzzyaabi1)then
call xyzzyaafi1
else
heg_cell(1,:)=a1
heg_cell(2,:)=a2
heg_cell(3,:)=a3
r_s=0.d0
if(am_master)then
call wout('System is a real system (has nuclei)')
call wout()
endif
endif
call xyzzyaafj1
return
100 call errstop_master('INIT_FREE_ORBS','Error reading correlation.da&
&ta .')
101 call errstop_master('INIT_FREE_ORBS','File correlation.data ended &
&unexpectedly.')
end subroutine init_free_orbs
subroutine xyzzyaafi1
use slaarnaan, only : init_geometry,netot_nitot_products,print_geometr&
&y
use slaarnabq,     only : min_image_brute_force
use slaarnaca,         only : read_ppots
implicit none
integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3
real(dp) xyzzyaaad3,xyzzyaaae3,xyzzyaaaf3,xyzzyaaag3,xyzzyaaah3(3),xyz&
&zyaaai3
real(dp),allocatable :: xyzzyaaaj3(:)
character(80) tmpr
if(am_master)then
call wout('Model system geometry setup')
call wout('===========================')
endif
a1=heg_cell(1,:)
a2=heg_cell(2,:)
a3=heg_cell(3,:)
xyzzyaaag3=real(sum(heg_nele),dp)/real(heg_nlayers,dp)
select case(periodicity)
case(3)
rho=3.d0/(fourpi*r_s**3)
case(2)
rho=1.d0/(pi*r_s**2)
case(1)
rho=0.5d0/r_s
case(0)
rho=0.d0
end select
select case(periodicity)
case(3)
xyzzyaaaf3=abs(a1(1)*a2(2)*a3(3)+a1(2)*a2(3)*a3(1)+a1(3)*a2(1)*a3(2)-a&
&1(3)*a2(2)*a3(1)-a1(1)*a2(3)*a3(2)-a1(2)*a2(1)*a3(3))
if(xyzzyaaaf3==0.d0)call errstop_master('SETUP_HEG_GEOMETRY','Cell vec&
&tors are linearly dependent.')
xyzzyaaag3=xyzzyaaag3/xyzzyaaaf3
xyzzyaaap1=(xyzzyaaag3/rho)**third
a1=a1*xyzzyaaap1
a2=a2*xyzzyaaap1
a3=a3*xyzzyaaap1
case(2)
xyzzyaaae3=abs(a1(1)*a2(2)-a1(2)*a2(1))
if(xyzzyaaae3==0.d0)call errstop_master('SETUP_HEG_GEOMETRY','Cell vec&
&tors are collinear.')
xyzzyaaag3=xyzzyaaag3/xyzzyaaae3
xyzzyaaap1=(xyzzyaaag3/rho)**half
a1=a1*xyzzyaaap1
a2=a2*xyzzyaaap1
case(1)
xyzzyaaad3=abs(a1(1))
if(xyzzyaaad3==0.d0)call errstop_master('SETUP_HEG_GEOMETRY','Cell vec&
&tor appears to be zero.')
xyzzyaaag3=xyzzyaaag3/xyzzyaaad3
xyzzyaaap1=xyzzyaaag3/rho
a1=a1*xyzzyaaap1
end select
do xyzzyaaaa3=1,3
if(scell_matrix(xyzzyaaaa3,xyzzyaaaa3)==0)scell_matrix(xyzzyaaaa3,xyzz&
&yaaaa3)=1
enddo
scell_mat_inv=inverse3(dble(scell_matrix),determinant=xyzzyaaai3)
npcells=abs(nint(xyzzyaaai3))
pa1=scell_mat_inv(1,1)*a1+scell_mat_inv(1,2)*a2+scell_mat_inv(1,3)*a3
pa2=scell_mat_inv(2,1)*a1+scell_mat_inv(2,2)*a2+scell_mat_inv(2,3)*a3
pa3=scell_mat_inv(3,1)*a1+scell_mat_inv(3,2)*a2+scell_mat_inv(3,3)*a3
call init_geometry
if(nbasis>0)then
call read_ppots
call netot_nitot_products
if(am_master)call print_geometry
endif
allocate(xyzzyaaaj3(nspin),stat=xyzzyaaab3)
call check_alloc(xyzzyaaab3,'SETUP_HEG_GEOMETRY','spin_polarization')
do xyzzyaaaa3=1,nspin
xyzzyaaaj3(xyzzyaaaa3)=-20.d0
if(popp_spin(xyzzyaaaa3)>xyzzyaaaa3)then
xyzzyaaaj3(xyzzyaaaa3)=real(heg_nele(xyzzyaaaa3)-heg_nele(popp_spin(xy&
&zzyaaaa3)),dp)/real(heg_nele(xyzzyaaaa3)+heg_nele(popp_spin(xyzzyaaaa&
&3)),dp)
if(pspin(xyzzyaaaa3)<0.d0)xyzzyaaaj3(xyzzyaaaa3)=-xyzzyaaaj3(xyzzyaaaa&
&3)
endif
enddo
if(isperiodic)then
xyzzyaacs1=transpose(bmat)
xyzzyaact1(1)=dot_product(binv(1:3,1),binv(1:3,1))
xyzzyaact1(2)=dot_product(binv(1:3,2),binv(1:3,2))
xyzzyaact1(3)=dot_product(binv(1:3,3),binv(1:3,3))
call min_image_brute_force(periodicity,k_offset,xyzzyaacs1,binv,xyzzya&
&aah3,xyzzyaact1)
k_offset=xyzzyaaah3
endif
if(am_master)then
call wout('Dimensionality                                 : ' //trim(i&
&2s(dimensionality)))
call wout('Periodicity                                    : ' //trim(i&
&2s(periodicity)))
if(isperiodic)then
tmpr=r2s(r_s,'(f20.12)')
call wout('r_s parameter                                  : '//trim(tm&
&pr))
tmpr=r2s(rho,'(f20.12)')
call wout('Density                                        : '//trim(tm&
&pr))
endif
do xyzzyaaaa3=1,nspin
if(xyzzyaaaj3(xyzzyaaaa3)<-1.d0)cycle
call wout('Spin polarization '//trim(i2s(xyzzyaaaa3))//'-'//trim(i2s(p&
&opp_spin(xyzzyaaaa3))) //'                          : ',xyzzyaaaj3(xy&
&zzyaaaa3),rfmt='(f20.12)',adjust=.true.)
enddo
if(isperiodic)then
call wout('Direct lattice vectors                         :')
call wout(' A1 =  ',a1(1:periodicity),rfmt='(e22.12)')
if(periodicity>1)then
call wout(' A2 =  ',a2(1:periodicity),rfmt='(e22.12)')
if(periodicity==3)call wout(' A3 =  ',a3(1:periodicity),rfmt='(e22.12)&
&')
endif
endif
if(heg_nlayers>1)then
do xyzzyaaac3=1,heg_nlayers
if(dimensionality==2)then
call wout('z (Layer '//trim(i2s(xyzzyaaac3)) //')                     &
&               : ',heg_zlayer(xyzzyaaac3),rfmt='(f20.12)',adjust=.tru&
&e.)
else
call wout('y (Wire '//trim(i2s(xyzzyaaac3)) //')                      &
&               : ',heg_ylayer(xyzzyaaac3),rfmt='(f20.12)',adjust=.tru&
&e.)
call wout('z (Wire '//trim(i2s(xyzzyaaac3)) //')                      &
&               : ',heg_zlayer(xyzzyaaac3),rfmt='(f20.12)',adjust=.tru&
&e.)
endif
enddo
endif
if(periodicity>0)then
call wout('Reciprocal lattice vectors                     :')
call wout(' B1 =  ',b1(1:periodicity),rfmt='(e22.12)')
if(periodicity>1)then
call wout(' B2 =  ',b2(1:periodicity),rfmt='(e22.12)')
if(periodicity==3)call wout(' B3 =  ',b3(1:periodicity),rfmt='(e22.12)&
&')
endif
endif
if(periodicity==3)then
tmpr=r2s(volume,'(es20.12)')
call wout('Simulation cell volume                         : '//trim(tm&
&pr))
elseif(periodicity==2)then
tmpr=r2s(area,'(es20.12)')
call wout('Simulation cell area                           : '//trim(tm&
&pr))
elseif(periodicity==1)then
tmpr=r2s(abs(a1(1)),'(es20.12)')
call wout('Simulation cell length                         : '//trim(tm&
&pr))
endif
if(periodicity>0)then
tmpr=r2s(wigner_seitz_radius,'(f20.13)')
if(periodicity==1)then
call wout('Radius of line inscribed in Wigner-Seitz cell  : '//trim(tm&
&pr))
elseif(periodicity==2)then
call wout('Radius of circle inscribed in Wigner-Seitz cell: '//trim(tm&
&pr))
else
call wout('Radius of sphere inscribed in Wigner-Seitz cell: '//trim(tm&
&pr))
endif
endif
if(xyzzyaabl1.and.(xyzzyaabf1.or.xyzzyaabg1))then
if(mc_twist_av)then
if(mc_twist_freq>0)then
call wout('Monte Carlo integration over k-vector offset will be perfor&
&med.')
call wout('Frequency of offset changes                    : ' //trim(i&
&2s(mc_twist_freq)))
elseif(mc_twist_freq<0)then
call errstop('SETUP_HEG_GEOMETRY','Monte Carlo twist frequency should &
&be non-negative.')
endif
call wout('Initial offset to k-vector lattice             : ')
else
call wout('Offset to k-vector lattice                     : ')
endif
call wout(' k_s = ',k_offset(1:periodicity))
endif
if(xyzzyaabo1)then
tmpr=r2s(xx_sep,'(f20.13)')
if(.not.fix_holes)then
call wout('Exciton-exciton separation                     : '//trim(tm&
&pr))
tmpr=r2s(me_biex3,'(f20.13)')
call wout('(Real) electron mass                           : '//trim(tm&
&pr))
tmpr=r2s(mh_biex3,'(f20.13)')
call wout('(Real) hole mass                               : '//trim(tm&
&pr))
else
call wout('Hole-hole separation                           : '//trim(tm&
&pr))
endif
endif
call wout()
endif
deallocate(xyzzyaaaj3)
if(xyzzyaabo1)then
nitot=1
allocate(rion(3,nitot),stat=xyzzyaaab3)
call check_alloc(xyzzyaaab3,'SETUP_HEG_GEOMETRY','rion')
rion=0.d0
call netot_nitot_products
endif
end subroutine xyzzyaafi1
subroutine xyzzyaafj1
implicit none
integer xyzzyaaaa4,xyzzyaaab4,xyzzyaaac4,xyzzyaaad4,xyzzyaaae4,xyzzyaa&
&af4,xyzzyaaag4,xyzzyaaah4,xyzzyaaai4,xyzzyaaaj4,xyzzyaaak4,xyzzyaaal4&
&,xyzzyaaam4,xyzzyaaan4,xyzzyaaao4,xyzzyaaap4,xyzzyaaaq4,xyzzyaaar4
real(dp) xyzzyaaas4,xyzzyaaat4,xyzzyaaau4,xyzzyaaav4,xyzzyaaaw4,xyzzya&
&aax4(3),xyzzyaaay4(heg_nlayers),xyzzyaaaz4
real(dp),parameter :: xyzzyaaba4=1.d-8
logical xyzzyaabb4,xyzzyaabc4
if(.not.any(heg_orbtype/=0.and.heg_orbtype<101))return
if(am_master)then
call wout('Free-particle orbital setup')
call wout('===========================')
endif
if(xyzzyaabm1.or.xyzzyaabn1)then
if(nspin/=4.or.any(nele/=1))call errstop_master('SETUP_HEG_ORBITALS','&
&Can only use biexciton with 1 spin-up electron, 1 spin-down electron,&
&1 spin-up hole and 1 spin-down hole.')
if(periodicity>0)call errstop_master('SETUP_HEG_ORBITALS','Biexciton s&
&ystem should not be periodic.')
if(abs(inv_pmass(2)-inv_pmass(1))>xyzzyaaba4*abs(inv_pmass(1)))call er&
&rstop_master('SETUP_HEG_ORBITALS','Masses of the two electrons should&
& be the same.')
if(abs(inv_pmass(4)-inv_pmass(3))>xyzzyaaba4*abs(inv_pmass(3)))call er&
&rstop_master('SETUP_HEG_ORBITALS','Masses of the two holes should be &
&the same.')
if(abs(pcharge(2)-pcharge(1))>xyzzyaaba4*abs(pcharge(1)))call errstop_&
&master('SETUP_HEG_ORBITALS','Charges of the two electrons should be t&
&he same.')
if(abs(pcharge(4)-pcharge(3))>xyzzyaaba4*abs(pcharge(3)))call errstop_&
&master('SETUP_HEG_ORBITALS','Charges of the two holes should be the s&
&ame.')
endif
if(xyzzyaabm1)then
if(dimensionality/=2)call errstop_master('SETUP_HEG_ORBITALS','Biexcit&
&on system should be 2D.')
if(heg_nlayers>1)then
if(heg_nlayers>2)call errstop_master('SETUP_HEG_ORBITALS','More than t&
&wo layers.')
if(heg_layer(1)==2.or.heg_layer(2)==2.or.heg_layer(3)==1 .or.heg_layer&
&(4)==1)call errstop_master('SETUP_HEG_ORBITALS','Particles 3 and 4 sh&
&ould be in the top layer.')
endif
if(any(heg_orbtype/=4))call errstop_master('SETUP_HEG_ORBITALS','Canno&
&t mix biexciton orbitals with other HEG orbitals.')
endif
if(xyzzyaabo1)then
if(nspin/=2.or.any(nele/=1))call errstop_master('SETUP_HEG_ORBITALS','&
&Can only use effective biexciton with 1 spin-up particle and 1 spin-d&
&own particle.')
if(periodicity>0)call errstop_master('SETUP_HEG_ORBITALS','Biexciton s&
&ystem should not be periodic.')
if(.not.fix_holes)then
if(me_biex3<=0.d0.or.mh_biex3<=0.d0)call errstop_master('SETUP_HEG_ORB&
&ITALS','Please define the me_biex3 and mh_biex3 parameters in the fre&
&e_particles input block.')
mu_biex3=me_biex3*mh_biex3/(me_biex3+mh_biex3)
if(abs(inv_pmass(1)-inv_pmass(2))>xyzzyaaba4*abs(inv_pmass(1)) .or.abs&
&(pmass(1)-mu_biex3)>xyzzyaaba4*abs(pmass(1)))call errstop_master('SET&
&UP_HEG_ORBITALS','Masses of the two effective particles should be the&
& reduced mass of an exciton.')
endif
if(any(pcharge/=0.d0))call errstop_master('SETUP_HEG_ORBITALS','Charge&
&s of the two particles should be zero.')
if(dimensionality/=2)call errstop_master('SETUP_HEG_ORBITALS','Effecti&
&ve biexciton system should be 2D.')
if(any(heg_layer/=1))call errstop_master('SETUP_HEG_ORBITALS','Should &
&just have a single layer in effective biexciton.')
if(any(heg_orbtype/=6))call errstop_master('SETUP_HEG_ORBITALS','Canno&
&t mix biexciton orbitals with other HEG orbitals.')
if(xx_sep<0.d0)call errstop_master('SETUP_HEG_ORBITALS','The xx_sep pa&
&rameter must be defined in the free_particles input block.')
have_biex3pot=.true.
endif
xyzzyaaav1=real(gautol,dp)*log(0.1d0)
xyzzyaaaw4=1.d0
do xyzzyaaag4=2,nemax
xyzzyaaaw4=xyzzyaaaw4*real(xyzzyaaag4,dp)**(-0.5d0/real(nemax,dp))
enddo
xyzzyaaaw4=orb_norm*xyzzyaaaw4
xyzzyaaaq1=xyzzyaaaw4
xyzzyaaat1=xyzzyaaaw4
if(type_wc==1.or.type_wc==3)then
xyzzyaaas1=orb_norm*exp(0.25d0)
else
xyzzyaaas1=xyzzyaaaw4
endif
xyzzyaaar1=xyzzyaaaq1**third
if(type_wc==2)then
do xyzzyaaai4=1,xyzzyaabw1
if(xyzzyaacx1(1,xyzzyaaai4)==0.d0)call errstop_master('SETUP_HEG_ORBIT&
&ALS','Constant term in Wigner-crystal plane-wave expansion should be &
&nonzero.')
call dscal(xyzzyaabs1-1,xyzzyaaas1/xyzzyaacx1(1,xyzzyaaai4),xyzzyaacx1&
&(2,xyzzyaaai4),1)
xyzzyaacx1(1,xyzzyaaai4)=xyzzyaaas1
enddo
endif
if(xyzzyaabf1.or.xyzzyaadi1.or.type_wc==2.or.type_wc==3)then
if(.not.isperiodic.and.(any(heg_nele>1).or.complex_wf))call errstop_ma&
&ster('SETUP_HEG_ORBITALS','Need a periodic system to generate plane w&
&aves.')
if(xyzzyaabf1)call xyzzyaafu1
if(xyzzyaadi1.or.type_wc==2.or.type_wc==3)call xyzzyaafw1
if(.not.xyzzyaabl1)then
do xyzzyaaao4=1,nspin
if(.not.(any(heg_orbtype(xyzzyaaao4,:)==1)))cycle
if(heg_nele(xyzzyaaao4)==0.or.any(xyzzyaaah1==heg_nele(xyzzyaaao4)))cy&
&cle
if(am_master)then
call wout()
call wordwrap('Number of particles of each type (e.g. up-spin electron&
&s) must be a "magic number" (occupying k-vectors by full stars of inc&
&reasing length) in order for the wave function to be real-valued. Set&
& COMPLEX_WF to T in the input file, or set the number of particles to&
& a magic number, which for the current lattice type are:')
xyzzyaaak4=size(xyzzyaaah1)-1
call write_list_int(xyzzyaaak4,xyzzyaaah1(1:xyzzyaaak4),12,5,1)
endif
call errstop_master('SETUP_HEG_ORBITALS','Wrong number of particles fo&
&r real-valued wave function.')
enddo
endif
if(xyzzyaaab1>0)then
xyzzyaaag4=0
do xyzzyaaao4=1,nspin
xyzzyaaag4=max(xyzzyaaag4,maxloc(xyzzyaaah1,1,xyzzyaaah1==heg_nele(xyz&
&zyaaao4)))
enddo
if(all(heg_orbtype==1))then
do xyzzyaaao4=1,nspin
if(.not.((any(heg_orbtype(xyzzyaaao4,:)==2).and.type_wc==2).or.(any(he&
&g_orbtype(xyzzyaaao4,:)<0).and.xyzzyaadi1)))cycle
if(heg_nele(xyzzyaaao4)<=xyzzyaaah1(xyzzyaaab1))cycle
call errstop_master('SETUP_HEG_ORBITALS','Number of coefficients in pl&
&ane-wave expansion must be equal to or larger than the number of occu&
&pied stars: '//trim(i2s(xyzzyaaag4))//'.')
enddo
endif
if(xyzzyaadi1.and.xyzzyaabk1)xyzzyaacw1(1,1:xyzzyaaag4,:)=1.d0
endif
endif
if(periodicity>0)then
if(.not.xyzzyaabh1)then
npcells=1
xyzzyaaav4=1.d0
else
npcells=heg_repeat**periodicity
xyzzyaaav4=real(heg_repeat,dp)
endif
scell_matrix=0
scell_matrix(1,1)=nint(xyzzyaaav4)
scell_matrix(2,2)=nint(xyzzyaaav4)
scell_matrix(3,3)=nint(xyzzyaaav4)
pa1=a1/xyzzyaaav4
pa2=a2/xyzzyaaav4
pa3=a3/xyzzyaaav4
pamat(1,1:3)=pa1(1:3)
pamat(2,1:3)=pa2(1:3)
pamat(3,1:3)=pa3(1:3)
pbmat=twopi*transpose(inverse3(pamat,determinant=pvolume))
pvolume=abs(pvolume)
pb1(1:3)=pbmat(1,1:3)
pb2(1:3)=pbmat(2,1:3)
pb3(1:3)=pbmat(3,1:3)
endif
if(xyzzyaabh1)then
xyzzyaabb4=.true.
xyzzyaabc4=.true.
select case(trim(heg_crystal_type))
case('cubic')
xyzzyaabc4=(periodicity==3)
xyzzyaabb4=xyzzyaabb4.and.a1(2)==0.d0.and.a1(3)==0.d0
xyzzyaabb4=xyzzyaabb4.and.a2(1)==0.d0.and.a2(3)==0.d0
xyzzyaabb4=xyzzyaabb4.and.a3(1)==0.d0.and.a3(2)==0.d0
case('bcc')
xyzzyaabc4=(periodicity==3)
xyzzyaabb4=xyzzyaabb4.and.-a1(1)==a1(2).and.-a1(1)==a1(3)
xyzzyaabb4=xyzzyaabb4.and.-a2(2)==a2(1).and.-a2(2)==a2(3)
xyzzyaabb4=xyzzyaabb4.and.-a3(3)==a3(1).and.-a3(3)==a3(2)
case('fcc')
xyzzyaabc4=(periodicity==3)
xyzzyaabb4=xyzzyaabb4.and.a1(1)==0.d0.and.a1(2)==a1(3)
xyzzyaabb4=xyzzyaabb4.and.a2(2)==0.d0.and.a2(1)==a2(3)
xyzzyaabb4=xyzzyaabb4.and.a3(3)==0.d0.and.a3(1)==a3(2)
case('manual')
case('rectangular')
xyzzyaabc4=(periodicity==2)
xyzzyaabb4=xyzzyaabb4.and.a1(2)==0.d0.and.a2(1)==0.d0
case('hexagonal')
xyzzyaabc4=(periodicity==2)
xyzzyaabb4=xyzzyaabb4.and.a1(2)==0.d0
xyzzyaabb4=xyzzyaabb4.and.abs(a1(1)*0.5d0-a2(1))<xyzzyaaba4
xyzzyaabb4=xyzzyaabb4.and.abs(a1(1)*root_three_over_two-a2(2))<xyzzyaa&
&ba4
case default
call errstop_master('SETUP_HEG_ORBITALS','Lattice not recognized.')
end select
if(.not.xyzzyaabc4)then
call errstop_master('SETUP_HEG_ORBITALS','Bad periodicity for '//trim(&
&heg_crystal_type)//' lattice.')
elseif(.not.xyzzyaabb4)then
call errstop_master('SETUP_HEG_ORBITALS','Lattice given in input not r&
&ecognized as '//trim(heg_crystal_type)//'.')
endif
if(heg_crystal_type/='manual')then
if(periodicity<2)call errstop_master('SETUP_HEG_ORBITALS','Need a 2D- &
&or 3D-periodic system to use pre-defined WC lattices.')
allocate(heg_nosites(heg_nbasis),stat=xyzzyaaac4)
call check_alloc(xyzzyaaac4,'SETUP_HEG_ORBITALS','heg_nosites')
heg_nosites=0
xyzzyaaac1=0
do xyzzyaaah4=1,ndet
do xyzzyaaai4=1,nspin
xyzzyaaaq4=heg_slatt(xyzzyaaai4,xyzzyaaah4)
if(xyzzyaaaq4==0)cycle
if(.not.heg_ferro(xyzzyaaaq4))then
if(trim(heg_crystal_type)=='fcc')call errstop_master('SETUP_HEG_ORBITA&
&LS','Antiferromagnet is frustrated in fcc lattices.')
if(popp_spin(xyzzyaaai4)==0)call errstop_master('SETUP_HEG_ORBITALS','&
&Cannot generate antiferromagnet if opposite-spin particle is not pres&
&ent.')
if(heg_slatt(popp_spin(xyzzyaaai4),xyzzyaaah4)/=xyzzyaaaq4)call errsto&
&p_master('SETUP_HEG_ORBITALS','Antiferromagnet requested for particle&
& '//trim(i2s(xyzzyaaai4))//' but opposite spin particle ('//trim(i2s(&
&popp_spin(xyzzyaaai4)))//') has not been assigned the same sublattice&
&.')
if(mod(npcells,2)/=0)call errstop_master('SETUP_HEG_ORBITALS','Cannot &
&generate antiferromagnet for an odd value of "repeat".')
if(heg_nele(xyzzyaaai4)/=npcells/2.and.am_master)call errstop_master('&
&SETUP_HEG_ORBITALS','Mismatch between number of particles ('//trim(i2&
&s(heg_nele(xyzzyaaai4)))//') and sites ('//trim(i2s(npcells/2))//') f&
&or sublattice '//trim(i2s(xyzzyaaaq4))//'.')
else
if(heg_nele(xyzzyaaai4)/=npcells)call errstop_master('SETUP_HEG_ORBITA&
&LS','Mismatch between number of particles ('//trim(i2s(heg_nele(xyzzy&
&aaai4)))//') and sites ('//trim(i2s(npcells))//') for sublattice '//t&
&rim(i2s(xyzzyaaaq4))//'.')
endif
xyzzyaaaj4=npcells
if(.not.heg_ferro(xyzzyaaaq4))xyzzyaaaj4=npcells/2
heg_nosites(xyzzyaaaq4)=xyzzyaaaj4
xyzzyaaac1=max(xyzzyaaac1,xyzzyaaaj4)
enddo
enddo
xyzzyaaap4=maxval(heg_slatt)
allocate(xyzzyaaay1(3,xyzzyaaac1,xyzzyaaap4),stat=xyzzyaaac4)
call check_alloc(xyzzyaaac4,'SETUP_HEG_ORBITALS','site_pos')
xyzzyaaay1=0.d0
if(periodicity==2)then
xyzzyaaar4=0
else
xyzzyaaar4=nint(xyzzyaaav4)-1
endif
do xyzzyaaah4=1,ndet
do xyzzyaaai4=1,nspin
xyzzyaaaq4=heg_slatt(xyzzyaaai4,xyzzyaaah4)
if(xyzzyaaaq4==0)cycle
if(.not.heg_ferro(xyzzyaaaq4).and.popp_spin(xyzzyaaai4)<xyzzyaaai4)cyc&
&le
xyzzyaaaa4=0
xyzzyaaab4=0
do xyzzyaaad4=0,nint(xyzzyaaav4)-1
do xyzzyaaae4=0,nint(xyzzyaaav4)-1
do xyzzyaaaf4=0,xyzzyaaar4
xyzzyaaax4=real((/xyzzyaaad4,xyzzyaaae4,xyzzyaaaf4/),dp)
xyzzyaaal4=xyzzyaaad4+xyzzyaaae4+xyzzyaaaf4
xyzzyaaas4=sum((xyzzyaaax4(:)+heg_wigner_basis(:,xyzzyaaaq4))*pamat(:,&
&1))
if(periodicity>1)then
xyzzyaaat4=sum((xyzzyaaax4(:)+heg_wigner_basis(:,xyzzyaaaq4))*pamat(:,&
&2))
else
xyzzyaaat4=heg_ylayer(heg_layer(xyzzyaaai4))
endif
if(periodicity==3)then
xyzzyaaau4=sum((xyzzyaaax4(:)+heg_wigner_basis(:,xyzzyaaaq4))*pamat(:,&
&3))
else
if(.not.heg_ferro(xyzzyaaaq4).and.mod(xyzzyaaal4,2)/=0)then
xyzzyaaau4=heg_zlayer(heg_layer(popp_spin(xyzzyaaai4)))
else
xyzzyaaau4=heg_zlayer(heg_layer(xyzzyaaai4))
endif
endif
if(.not.heg_ferro(xyzzyaaaq4))then
if(mod(xyzzyaaal4,2)==0)then
xyzzyaaaa4=xyzzyaaaa4+1
xyzzyaaay1(1:3,xyzzyaaaa4,xyzzyaaaq4)=(/xyzzyaaas4,xyzzyaaat4,xyzzyaaa&
&u4/)
else
xyzzyaaab4=xyzzyaaab4+1
xyzzyaaay1(1:3,xyzzyaaab4,heg_slatt(popp_spin(xyzzyaaai4),xyzzyaaah4))&
&=(/xyzzyaaas4,xyzzyaaat4,xyzzyaaau4/)
endif
else
xyzzyaaaa4=xyzzyaaaa4+1
xyzzyaaay1(1:3,xyzzyaaaa4,xyzzyaaaq4)=(/xyzzyaaas4,xyzzyaaat4,xyzzyaaa&
&u4/)
endif
enddo
enddo
enddo
enddo
enddo
else
xyzzyaaac1=0
do xyzzyaaah4=1,ndet
do xyzzyaaai4=1,nspin
xyzzyaaaq4=heg_slatt(xyzzyaaai4,xyzzyaaah4)
if(xyzzyaaaq4==0)cycle
xyzzyaaaj4=heg_nosites(xyzzyaaaq4)
if(xyzzyaaaj4/=heg_nele(xyzzyaaai4))call errstop_master('SETUP_HEG_ORB&
&ITALS','Mismatch between number of particles ('//trim(i2s(heg_nele(xy&
&zzyaaai4)))//') and sites ('//trim(i2s(xyzzyaaaj4))//') for sublattic&
&e '//trim(i2s(xyzzyaaaq4))//'.')
xyzzyaaac1=max(xyzzyaaac1,xyzzyaaaj4)
enddo
enddo
xyzzyaaap4=maxval(heg_slatt)
allocate(xyzzyaaay1(3,xyzzyaaac1,xyzzyaaap4),stat=xyzzyaaac4)
call check_alloc(xyzzyaaac4,'SETUP_HEG_ORBITALS','site_pos')
xyzzyaaay1=0.d0
do xyzzyaaah4=1,ndet
do xyzzyaaai4=1,nspin
xyzzyaaaq4=heg_slatt(xyzzyaaai4,xyzzyaaah4)
if(xyzzyaaaq4==0)cycle
xyzzyaaaj4=heg_nosites(xyzzyaaaq4)
do xyzzyaaaa4=1,xyzzyaaaj4
do xyzzyaaag4=1,3
if(xyzzyaaag4<=periodicity)then
xyzzyaaay1(xyzzyaaag4,xyzzyaaaa4,xyzzyaaaq4)=sum(heg_crystal_sites(1:3&
&,xyzzyaaaa4,xyzzyaaaq4)*amat(1:3,xyzzyaaag4))
else
xyzzyaaay1(xyzzyaaag4,xyzzyaaaa4,xyzzyaaaq4)=heg_crystal_sites(xyzzyaa&
&ag4,xyzzyaaaa4,xyzzyaaaq4)
endif
enddo
enddo
enddo
enddo
endif
call xyzzyaafy1
if(nitot==0)call xyzzyaafk1
endif
if(xyzzyaabj1)call sdw_setup
if(xyzzyaabg1)then
if(xyzzyaadk1)then
if(xyzzyaacp1==0.d0)then
xyzzyaacp1=0.99d0*wigner_seitz_radius
do xyzzyaaai4=1,xyzzyaabv1
xyzzyaadb1(:,1,xyzzyaaai4)=xyzzyaadb1(:,0,xyzzyaaai4)*real(xyzzyaaby1,&
&dp)/xyzzyaacp1
xyzzyaaaz4=sum(xyzzyaadb1(:,:,xyzzyaaai4)**2)
if(xyzzyaadj1)xyzzyaaaz4=xyzzyaaaz4+sum(xyzzyaacu1(:,:,xyzzyaaai4)**2)
if(xyzzyaadi1)xyzzyaaaz4=xyzzyaaaz4+sum(xyzzyaacw1(:,:,xyzzyaaai4)**2)
if(xyzzyaadl1)xyzzyaaaz4=xyzzyaaaz4+sum(xyzzyaadc1(:,:,xyzzyaaai4)**2)&
&+sum(xyzzyaadf1(:,:,:,xyzzyaaai4)**2)
if(xyzzyaaaz4<=0.d0)call errstop_master('SETUP_HEG_ORBITALS','Pairing &
&orbital is zero.')
xyzzyaada1(xyzzyaaai4)=1.d0/sqrt(xyzzyaaaz4)
enddo
endif
if(xyzzyaadn1.and.xyzzyaacp1>wigner_seitz_radius)call errstop_master('&
&SETUP_HEG_ORBITALS','L_p is limited to < L_WS, but is > L_WS.')
xyzzyaacq1=1.d0/xyzzyaacp1
xyzzyaacr1=xyzzyaacq1*xyzzyaacq1
endif
call xyzzyaafz1
endif
xyzzyaaau1=0.d0
if(isperiodic.and.periodicity<3.and.heg_nlayers>1)then
xyzzyaaay4=0.d0
do xyzzyaaai4=1,nspin
xyzzyaaay4(heg_layer(xyzzyaaai4))=xyzzyaaay4(heg_layer(xyzzyaaai4)) +d&
&ble(nele(xyzzyaaai4))*pcharge(xyzzyaaai4)
enddo
xyzzyaaau1=0.d0
if(dimensionality==2)then
do xyzzyaaam4=1,heg_nlayers-1
do xyzzyaaan4=xyzzyaaam4+1,heg_nlayers
xyzzyaaau1=xyzzyaaau1+abs(heg_zlayer(xyzzyaaan4)-heg_zlayer(xyzzyaaam4&
&)) *xyzzyaaay4(xyzzyaaam4)*xyzzyaaay4(xyzzyaaan4)
enddo
enddo
xyzzyaaau1=twopi*xyzzyaaau1/area
else
do xyzzyaaam4=1,heg_nlayers-1
do xyzzyaaan4=xyzzyaaam4+1,heg_nlayers
xyzzyaaau1=xyzzyaaau1+xyzzyaaay4(xyzzyaaam4)*xyzzyaaay4(xyzzyaaan4) *l&
&og((heg_ylayer(xyzzyaaam4)-heg_ylayer(xyzzyaaan4))**2 +(heg_zlayer(xy&
&zzyaaam4)-heg_zlayer(xyzzyaaan4))**2)
enddo
enddo
xyzzyaaau1=xyzzyaaau1/abs(a1(1))
endif
endif
constant_energy=constant_energy+xyzzyaaau1
if(am_master)then
if(xyzzyaabh1)then
call wout('Type of crystal             : '//trim(heg_crystal_type))
call wout('Number of sublattices       : '//trim(i2s(heg_nbasis)))
if(isperiodic)then
call wout('Primitive cell              :')
call wout('',pa1(1:periodicity),rfmt='(f20.12)',rsep=' ')
if(periodicity>1)then
call wout('',pa2(1:periodicity),rfmt='(f20.12)',rsep=' ')
if(periodicity==3)call wout('',pa3(1:periodicity),rfmt='(f20.12)',rsep&
&=' ')
endif
endif
if(type_wc==1.or.type_wc==3)call wout('Range of WC Gaussians       : '&
&//trim(i2s(maxval(xyzzyaaak1(:))))//' cells (max)')
if(nitot>0)call wout('Number of centres ("ions")  : '//trim(i2s(nitot)&
&))
endif
if(xyzzyaadj1)call wout('Range of pairing Gaussians  : '//trim(i2s(xyz&
&zyaaad1))//' cells')
if(xyzzyaadk1)call wout('Range of pairing polynomial : '//trim(i2s(xyz&
&zyaaad1))//' cells')
if(xyzzyaabj1)call print_sdwdata()
call wout('Free-particle orbitals set up.')
call wout()
endif
end subroutine xyzzyaafj1
subroutine xyzzyaafk1
use slaarnaan, only : netot_nitot_products
use slaarnabg, only : ignore_ionic_interactions
implicit none
integer xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5,xyzzyaaad5,xyzzyaaae5
real(dp),allocatable :: xyzzyaaaf5(:,:)
real(dp),parameter :: xyzzyaaag5=1.d-7
logical xyzzyaaah5
allocate(xyzzyaaaf5(3,xyzzyaaac1*nspin_ndet),stat=xyzzyaaae5)
call check_alloc(xyzzyaaae5,'SET_UP_INHOMOGENEITY','1')
nitot=0
xyzzyaaaa5=maxval(heg_slatt)
do xyzzyaaab5=1,xyzzyaaaa5
do xyzzyaaac5=1,heg_nosites(xyzzyaaab5)
xyzzyaaah5=.true.
do xyzzyaaad5=1,nitot
if(all(abs(xyzzyaaaf5(:,xyzzyaaad5)-xyzzyaaay1(:,xyzzyaaac5,xyzzyaaab5&
&)) <xyzzyaaag5*r_s))then
xyzzyaaah5=.false.
exit
endif
enddo
if(xyzzyaaah5)then
nitot=nitot+1
xyzzyaaaf5(:,nitot)=xyzzyaaay1(:,xyzzyaaac5,xyzzyaaab5)
endif
enddo
enddo
allocate(rion(3,nitot),stat=xyzzyaaae5)
call check_alloc(xyzzyaaae5,'SET_UP_INHOMOGENEITY','2')
rion(1:3,1:nitot)=xyzzyaaaf5(1:3,1:nitot)
deallocate(xyzzyaaaf5)
ignore_ionic_interactions=.true.
call netot_nitot_products
end subroutine xyzzyaafk1
subroutine write_heg(correlation_name)
implicit none
character(20),intent(in) :: correlation_name
integer xyzzyaaaa6,xyzzyaaab6,xyzzyaaac6,xyzzyaaad6,xyzzyaaae6,xyzzyaa&
&af6,xyzzyaaag6
logical xyzzyaaah6
if(.not.am_master)return
inquire(file=trim(correlation_name),exist=xyzzyaaah6)
if(xyzzyaaah6)then
open(unit=xyzzyaaaa1,file=trim(correlation_name),status='old',position&
&='append',iostat=xyzzyaaaa6)
else
open(unit=xyzzyaaaa1,file=trim(correlation_name),status='replace',iost&
&at=xyzzyaaaa6)
endif
write(xyzzyaaaa1,*)'START FREE_ORBS'
write(xyzzyaaaa1,*)'Title'
write(xyzzyaaaa1,*)trim(adjustl(title))
if(xyzzyaabg1)then
write(xyzzyaaaa1,*)'START PAIRING'
write(xyzzyaaaa1,*)'Spin-pair dependence (CURRENTLY HAS NO EFFECT)'
write(xyzzyaaaa1,'(3x,a)')trim(i2s(xyzzyaabu1))
if(xyzzyaadj1)then
write(xyzzyaaaa1,*)'START GAUSSIAN TERM'
write(xyzzyaaaa1,*)'Number of gaussians'
write(xyzzyaaaa1,'(3x,a)')trim(i2s(xyzzyaabq1))
write(xyzzyaaaa1,*)'Fit to exponential (T/F)'
write(xyzzyaaaa1,'(3x,a)')trim(l2s(xyzzyaadm1))
write(xyzzyaaaa1,*)'Parameter ;          Optimizable (0=NO; 1=YES)'
if(xyzzyaadm1)then
do xyzzyaaac6=1,xyzzyaabv1
write(xyzzyaaaa1,*)xyzzyaaco1,xyzzyaabt1,'      ! Rex'
enddo
endif
do xyzzyaaac6=1,xyzzyaabv1
do xyzzyaaab6=1,xyzzyaabq1
write(xyzzyaaaa1,*)xyzzyaacu1(1:real1_complex2,xyzzyaaab6,xyzzyaaac6)*&
&xyzzyaada1(xyzzyaaac6),xyzzyaacc1(xyzzyaaab6,xyzzyaaac6),'      ! g_'&
&,trim(i2s(xyzzyaaab6)),',',trim(i2s(xyzzyaaac6))
write(xyzzyaaaa1,*)xyzzyaacv1(xyzzyaaab6,xyzzyaaac6),xyzzyaacd1(xyzzya&
&aab6,xyzzyaaac6),'      ! exp_',trim(i2s(xyzzyaaab6)),',',trim(i2s(xy&
&zzyaaac6))
enddo
enddo
write(xyzzyaaaa1,*)'END GAUSSIAN TERM'
endif
if(xyzzyaadi1)then
write(xyzzyaaaa1,*)'START PLANE-WAVE TERM'
write(xyzzyaaaa1,*)'Number of plane waves'
write(xyzzyaaaa1,'(3x,a)')trim(i2s(xyzzyaabr1))
write(xyzzyaaaa1,*)'Parameter ;          Optimizable (0=NO; 1=YES)'
do xyzzyaaac6=1,xyzzyaabv1
do xyzzyaaab6=1,xyzzyaabr1
write(xyzzyaaaa1,*)xyzzyaacw1(1:real1_complex2,xyzzyaaab6,xyzzyaaac6)*&
&xyzzyaada1(xyzzyaaac6),xyzzyaacf1(xyzzyaaab6,xyzzyaaac6),'      ! c_'&
&,trim(i2s(xyzzyaaab6)),',',trim(i2s(xyzzyaaac6))
enddo
enddo
write(xyzzyaaaa1,*)'END PLANE-WAVE TERM'
endif
if(xyzzyaadk1)then
write(xyzzyaaaa1,*)'START POLYNOMIAL TERM'
write(xyzzyaaaa1,*)'Order of polynomials'
write(xyzzyaaaa1,'(3x,a)')trim(i2s(xyzzyaabx1))
write(xyzzyaaaa1,*)'Truncation order of polynomial'
write(xyzzyaaaa1,'(3x,a)')trim(i2s(xyzzyaaby1))
xyzzyaaag6=xyzzyaabz1
if(xyzzyaadn1)xyzzyaaag6=2
write(xyzzyaaaa1,*)'Parameter ;          Optimizable (0=NO; 1=YES)'
write(xyzzyaaaa1,*)xyzzyaacp1,xyzzyaaag6,'      ! L_p'
do xyzzyaaac6=1,xyzzyaabv1
do xyzzyaaab6=0,xyzzyaabx1
if(xyzzyaach1(xyzzyaaab6,xyzzyaaac6)<0)cycle
write(xyzzyaaaa1,*)xyzzyaadb1(1:real1_complex2,xyzzyaaab6,xyzzyaaac6)*&
&xyzzyaada1(xyzzyaaac6),xyzzyaach1(xyzzyaaab6,xyzzyaaac6),'      ! alp&
&ha_',trim(i2s(xyzzyaaab6)),',',trim(i2s(xyzzyaaac6))
enddo
enddo
write(xyzzyaaaa1,*)'END POLYNOMIAL TERM'
endif
if(xyzzyaadl1)then
write(xyzzyaaaa1,*)'START SLATER TERM'
write(xyzzyaaaa1,*)'Order of S-type expansion'
write(xyzzyaaaa1,'(3x,a)')trim(i2s(xyzzyaaca1))
write(xyzzyaaaa1,*)'Order of P-type expansion'
write(xyzzyaaaa1,'(3x,a)')trim(i2s(xyzzyaacb1))
write(xyzzyaaaa1,*)'Parameter ;          Optimizable (0=NO; 1=YES)'
do xyzzyaaac6=1,xyzzyaabv1
do xyzzyaaab6=1,xyzzyaaca1
write(xyzzyaaaa1,*)xyzzyaadc1(1:real1_complex2,xyzzyaaab6,xyzzyaaac6)*&
&xyzzyaada1(xyzzyaaac6),xyzzyaaci1(xyzzyaaab6,xyzzyaaac6),'      ! S_c&
&_',trim(i2s(xyzzyaaab6)),',',trim(i2s(xyzzyaaac6))
write(xyzzyaaaa1,*)xyzzyaadd1(xyzzyaaab6,xyzzyaaac6),xyzzyaacj1(xyzzya&
&aab6,xyzzyaaac6),'      ! S_a_',trim(i2s(xyzzyaaab6)),',',trim(i2s(xy&
&zzyaaac6))
write(xyzzyaaaa1,*)xyzzyaade1(xyzzyaaab6,xyzzyaaac6),xyzzyaack1(xyzzya&
&aab6,xyzzyaaac6),'      ! S_b_',trim(i2s(xyzzyaaab6)),',',trim(i2s(xy&
&zzyaaac6))
enddo
do xyzzyaaab6=1,xyzzyaacb1
do xyzzyaaad6=1,dimensionality
write(xyzzyaaaa1,*)xyzzyaadf1(1:real1_complex2,xyzzyaaad6,xyzzyaaab6,x&
&yzzyaaac6)*xyzzyaada1(xyzzyaaac6),xyzzyaacl1(xyzzyaaad6,xyzzyaaab6,xy&
&zzyaaac6),'      ! P_c_',trim(i2s(xyzzyaaad6)),',',   trim(i2s(xyzzya&
&aab6)),',',trim(i2s(xyzzyaaac6))
enddo
write(xyzzyaaaa1,*)xyzzyaadg1(xyzzyaaab6,xyzzyaaac6),xyzzyaacm1(xyzzya&
&aab6,xyzzyaaac6),'      ! P_a_',trim(i2s(xyzzyaaab6)),',',trim(i2s(xy&
&zzyaaac6))
write(xyzzyaaaa1,*)xyzzyaadh1(xyzzyaaab6,xyzzyaaac6),xyzzyaacn1(xyzzya&
&aab6,xyzzyaaac6),'      ! P_b_',trim(i2s(xyzzyaaab6)),',',trim(i2s(xy&
&zzyaaac6))
enddo
enddo
write(xyzzyaaaa1,*)'END SLATER TERM'
endif
write(xyzzyaaaa1,*)'END PAIRING'
endif
if(xyzzyaabh1)then
write(xyzzyaaaa1,*)'START WIGNER CRYSTAL'
write(xyzzyaaaa1,*)'Type of WC orbitals ("gaussian", "plane-wave", or &
&"hybrid")'
if(type_wc==1)then
write(xyzzyaaaa1,*)'  gaussian'
elseif(type_wc==2)then
write(xyzzyaaaa1,*)'  plane-wave'
elseif(type_wc==3)then
write(xyzzyaaaa1,*)'  hybrid'
else
call errstop('WRITE_HEG','Bug.')
endif
write(xyzzyaaaa1,*)'Spin dependence (0 -> eu=ed ; 1 -> eu/=ed)'
write(xyzzyaaaa1,'(3x,a)')trim(i2s(spin_dep_wc))
if(type_wc==1.or.type_wc==3)then
write(xyzzyaaaa1,*)'Gaussian exponent(s) ; Optimizable (0=NO; 1=YES)'
do xyzzyaaac6=1,xyzzyaabw1
write(xyzzyaaaa1,*)wc_gauss_exp(xyzzyaaac6),xyzzyaace1(xyzzyaaac6),'      !&
& exp_' //trim(i2s(xyzzyaaac6))
enddo
endif
if(type_wc==2.or.type_wc==3)then
write(xyzzyaaaa1,*)'Number of plane waves'
write(xyzzyaaaa1,'(3x,a)')trim(i2s(xyzzyaabs1))
write(xyzzyaaaa1,*)'PW coefficients      ; Optimizable (0=NO; 1=YES)'
do xyzzyaaac6=1,xyzzyaabw1
do xyzzyaaab6=1,xyzzyaabs1
write(xyzzyaaaa1,*)xyzzyaacx1(xyzzyaaab6,xyzzyaaac6),xyzzyaacg1(xyzzya&
&aab6,xyzzyaaac6),'      ! c_' //trim(i2s(xyzzyaaab6))//','//trim(i2s(&
&xyzzyaaac6))
enddo
enddo
endif
write(xyzzyaaaa1,*)'END WIGNER CRYSTAL'
endif
if(xyzzyaabj1)then
write(xyzzyaaaa1,*)'START SPIN DENSITY WAVE'
write(xyzzyaaaa1,*)'Number of k-vectors'
xyzzyaaae6=1
do xyzzyaaab6=2,nele(1)
if(any(sdw_kvec(:,xyzzyaaab6)/=sdw_kvec(:,xyzzyaaab6-1)))then
xyzzyaaae6=xyzzyaaae6+1
endif
enddo
write(xyzzyaaaa1,*)xyzzyaaae6
write(xyzzyaaaa1,*)'k-vector (3 components) ; theta, occupancy'
do xyzzyaaab6=1,nele(1)
xyzzyaaaf6=1
if(xyzzyaaab6<nele(1))then
if(.not.any(sdw_kvec(:,xyzzyaaab6)/=sdw_kvec(:,xyzzyaaab6+1)))then
xyzzyaaaf6=2
endif
endif
if(xyzzyaaab6>1)then
if(.not.any(sdw_kvec(:,xyzzyaaab6)/=sdw_kvec(:,xyzzyaaab6-1)))then
cycle
endif
endif
write(xyzzyaaaa1,*)sdw_kvec(:,xyzzyaaab6)
write(xyzzyaaaa1,*)sdw_orb_theta(xyzzyaaab6),xyzzyaaaf6
enddo
write(xyzzyaaaa1,*)'END SPIN DENSITY WAVE'
endif
write(xyzzyaaaa1,*)'END FREE_ORBS'
close(xyzzyaaaa1)
end subroutine write_heg
subroutine setup_freeorb_params(nparam)
implicit none
integer,intent(inout) :: nparam
integer xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7
nparam=0
if(xyzzyaabg1)then
if(xyzzyaadj1)then
if(xyzzyaadm1)then
do xyzzyaaab7=1,xyzzyaabv1
if(xyzzyaabt1==1)nparam=nparam+1
enddo
else
do xyzzyaaab7=1,xyzzyaabv1
do xyzzyaaaa7=1,xyzzyaabq1
if(xyzzyaacc1(xyzzyaaaa7,xyzzyaaab7)==1)nparam=nparam+real1_complex2
if(xyzzyaacd1(xyzzyaaaa7,xyzzyaaab7)==1)nparam=nparam+1
enddo
enddo
endif
endif
if(xyzzyaadi1)then
do xyzzyaaab7=1,xyzzyaabv1
do xyzzyaaaa7=1,xyzzyaabr1
if(xyzzyaacf1(xyzzyaaaa7,xyzzyaaab7)==1)nparam=nparam+real1_complex2
enddo
enddo
endif
if(xyzzyaadk1)then
if(xyzzyaabz1==1)nparam=nparam+1
do xyzzyaaab7=1,xyzzyaabv1
do xyzzyaaaa7=0,xyzzyaabx1
if(xyzzyaach1(xyzzyaaaa7,xyzzyaaab7)==1)nparam=nparam+real1_complex2
enddo
enddo
endif
if(xyzzyaadl1)then
do xyzzyaaab7=1,xyzzyaabv1
do xyzzyaaaa7=1,xyzzyaaca1
if(xyzzyaaci1(xyzzyaaaa7,xyzzyaaab7)==1)nparam=nparam+real1_complex2
if(xyzzyaacj1(xyzzyaaaa7,xyzzyaaab7)==1)nparam=nparam+1
if(xyzzyaack1(xyzzyaaaa7,xyzzyaaab7)==1)nparam=nparam+1
enddo
do xyzzyaaaa7=1,xyzzyaacb1
do xyzzyaaac7=1,dimensionality
if(xyzzyaacl1(xyzzyaaac7,xyzzyaaaa7,xyzzyaaab7)==1)nparam=nparam+real1&
&_complex2
enddo
if(xyzzyaacm1(xyzzyaaaa7,xyzzyaaab7)==1)nparam=nparam+1
if(xyzzyaacn1(xyzzyaaaa7,xyzzyaaab7)==1)nparam=nparam+1
enddo
enddo
endif
endif
if(xyzzyaabh1)then
if(type_wc==1.or.type_wc==3)then
do xyzzyaaab7=1,xyzzyaabw1
if(xyzzyaace1(xyzzyaaab7)==1)nparam=nparam+1
enddo
endif
if(type_wc==2.or.type_wc==3)then
do xyzzyaaab7=1,xyzzyaabw1
do xyzzyaaaa7=1,xyzzyaabs1
if(xyzzyaacg1(xyzzyaaaa7,xyzzyaaab7)==1)nparam=nparam+1
enddo
enddo
endif
endif
xyzzyaabp1=nparam
call xyzzyaafl1
end subroutine setup_freeorb_params
subroutine finish_freeorb_params
implicit none
call xyzzyaafm1
end subroutine finish_freeorb_params
subroutine xyzzyaafl1
implicit none
integer xyzzyaaaa9,xyzzyaaab9
xyzzyaaab9=xyzzyaabp1
if(xyzzyaabg1)then
if(xyzzyaadj1)then
allocate(xyzzyaaep1(real1_complex2,xyzzyaabq1,xyzzyaabv1,0:xyzzyaaab9)&
&,xyzzyaaeq1(xyzzyaabq1,xyzzyaabv1,0:xyzzyaaab9),stat=xyzzyaaaa9)
call check_alloc(xyzzyaaaa9,'SETUP_FREEORB_PBUFFER','pairing_gauss')
xyzzyaaep1=0.d0
xyzzyaaeq1=0.d0
xyzzyaaeg1=real1_complex2*xyzzyaabq1*xyzzyaabv1
xyzzyaaeh1=xyzzyaabq1*xyzzyaabv1
endif
if(xyzzyaadi1)then
allocate(xyzzyaaer1(real1_complex2,xyzzyaabr1,xyzzyaabv1,0:xyzzyaaab9)&
&,stat=xyzzyaaaa9)
call check_alloc(xyzzyaaaa9,'SETUP_FREEORB_PBUFFER','pairing_pw')
xyzzyaaer1=0.d0
xyzzyaaei1=real1_complex2*xyzzyaabr1*xyzzyaabv1
endif
if(xyzzyaadk1)then
allocate(xyzzyaaes1(0:xyzzyaaab9),xyzzyaaet1(0:xyzzyaaab9),           &
&      xyzzyaaeu1(0:xyzzyaaab9),xyzzyaaev1(real1_complex2,0:xyzzyaabx1&
&,xyzzyaabv1,0:xyzzyaaab9),stat=xyzzyaaaa9)
call check_alloc(xyzzyaaaa9,'SETUP_FREEORB_PBUFFER','pairing_poly')
xyzzyaaes1=0.d0
xyzzyaaet1=0.d0
xyzzyaaeu1=0.d0
xyzzyaaev1=0.d0
xyzzyaaej1=real1_complex2*(xyzzyaabx1+1)*xyzzyaabv1
endif
if(xyzzyaadl1)then
if(xyzzyaaca1>0)then
allocate(xyzzyaaew1(real1_complex2,xyzzyaaca1,xyzzyaabv1,0:xyzzyaaab9)&
&,xyzzyaaex1(xyzzyaaca1,xyzzyaabv1,0:xyzzyaaab9),  xyzzyaaey1(xyzzyaac&
&a1,xyzzyaabv1,0:xyzzyaaab9),stat=xyzzyaaaa9)
call check_alloc(xyzzyaaaa9,'SETUP_FREEORB_PBUFFER','pairing_slater_s'&
&)
xyzzyaaew1=0.d0
xyzzyaaex1=0.d0
xyzzyaaey1=0.d0
xyzzyaaek1=real1_complex2*xyzzyaaca1*xyzzyaabv1
xyzzyaael1=xyzzyaaca1*xyzzyaabv1
endif
if(xyzzyaacb1>0)then
allocate(xyzzyaaez1(real1_complex2,3,xyzzyaacb1,xyzzyaabv1,0:xyzzyaaab&
&9),xyzzyaafa1(xyzzyaacb1,xyzzyaabv1,0:xyzzyaaab9),    xyzzyaafb1(xyzz&
&yaacb1,xyzzyaabv1,0:xyzzyaaab9),stat=xyzzyaaaa9)
call check_alloc(xyzzyaaaa9,'SETUP_FREEORB_PBUFFER','pairing_slater_p'&
&)
xyzzyaaez1=0.d0
xyzzyaafa1=0.d0
xyzzyaafb1=0.d0
xyzzyaaem1=real1_complex2*3*xyzzyaacb1*xyzzyaabv1
xyzzyaaen1=xyzzyaacb1*xyzzyaabv1
endif
endif
allocate(xyzzyaafc1(xyzzyaabv1,0:xyzzyaaab9),stat=xyzzyaaaa9)
call check_alloc(xyzzyaaaa9,'SETUP_FREEORB_PBUFFER','pair_norm')
xyzzyaafc1=0.d0
endif
if(xyzzyaabh1)then
if(type_wc==1.or.type_wc==3)then
allocate(xyzzyaafd1(xyzzyaabw1,0:xyzzyaaab9),stat=xyzzyaaaa9)
call check_alloc(xyzzyaaaa9,'SETUP_FREEORB_PBUFFER','wc_gauss_exp')
xyzzyaafd1=0.d0
endif
if(type_wc==2.or.type_wc==3)then
allocate(xyzzyaafe1(xyzzyaabs1,xyzzyaabw1,0:xyzzyaaab9),stat=xyzzyaaaa&
&9)
call check_alloc(xyzzyaaaa9,'SETUP_FREEORB_PBUFFER','wc_pw_c')
xyzzyaafe1=0.d0
xyzzyaaeo1=xyzzyaabs1*xyzzyaabw1
endif
endif
end subroutine xyzzyaafl1
subroutine xyzzyaafm1
implicit none
if(xyzzyaabg1)then
if(xyzzyaadj1)deallocate(xyzzyaaep1,xyzzyaaeq1)
if(xyzzyaadi1)deallocate(xyzzyaaer1)
if(xyzzyaadk1)deallocate(xyzzyaaes1,xyzzyaaet1,xyzzyaaeu1,xyzzyaaev1)
if(xyzzyaadl1)then
if(xyzzyaaca1>0)deallocate(xyzzyaaew1,xyzzyaaex1,xyzzyaaey1)
if(xyzzyaacb1>0)deallocate(xyzzyaaez1,xyzzyaafa1,xyzzyaafb1)
endif
deallocate(xyzzyaafc1)
endif
if(xyzzyaabh1)then
if(type_wc==1.or.type_wc==3)deallocate(xyzzyaafd1)
if(type_wc==2.or.type_wc==3)deallocate(xyzzyaafe1)
endif
end subroutine xyzzyaafm1
subroutine get_freeorb_params(params,has_lolim,lolim,has_hilim,hilim,i&
&s_shallow,is_redundant,is_linear,is_loglinear,has_aderiv,affect_map,l&
&abel)
implicit none
real(dp),intent(inout) :: params(xyzzyaabp1),lolim(xyzzyaabp1),hilim(x&
&yzzyaabp1)
logical,intent(inout) :: has_lolim(xyzzyaabp1),has_hilim(xyzzyaabp1),i&
&s_shallow(xyzzyaabp1),is_redundant(xyzzyaabp1),is_linear(xyzzyaabp1),&
&is_loglinear(xyzzyaabp1),has_aderiv(xyzzyaabp1),affect_map(xyzzyaabp1&
&,xyzzyaabp1)
character(2),intent(inout) :: label(xyzzyaabp1)
integer xyzzyaaaa11,xyzzyaaab11,xyzzyaaac11,xyzzyaaad11,xyzzyaaae11
real(dp) xyzzyaaaf11,xyzzyaaag11
has_lolim=.false.
lolim=0.d0
has_hilim=.false.
hilim=0.d0
is_shallow=.false.
is_redundant=.false.
is_linear=.false.
is_loglinear=.false.
has_aderiv=.false.
affect_map=.false.
do xyzzyaaaa11=1,xyzzyaabp1
affect_map(xyzzyaaaa11,xyzzyaaaa11)=.true.
enddo
label='OP'
xyzzyaaaf11=1.1d-8
if(isperiodic)xyzzyaaag11=0.999999d0*wigner_seitz_radius
xyzzyaaaa11=0
if(xyzzyaabg1)then
do xyzzyaaac11=1,xyzzyaabv1
if(xyzzyaadj1)xyzzyaacu1(:,:,xyzzyaaac11)=xyzzyaacu1(:,:,xyzzyaaac11)*&
&xyzzyaada1(xyzzyaaac11)
if(xyzzyaadi1)xyzzyaacw1(:,:,xyzzyaaac11)=xyzzyaacw1(:,:,xyzzyaaac11)*&
&xyzzyaada1(xyzzyaaac11)
if(xyzzyaadk1)xyzzyaadb1(:,:,xyzzyaaac11)=xyzzyaadb1(:,:,xyzzyaaac11)*&
&xyzzyaada1(xyzzyaaac11)
if(xyzzyaadl1)xyzzyaadc1(:,:,xyzzyaaac11)=xyzzyaadc1(:,:,xyzzyaaac11)*&
&xyzzyaada1(xyzzyaaac11)
xyzzyaada1(xyzzyaaac11)=1.d0
enddo
if(xyzzyaadj1)then
if(xyzzyaadm1)then
do xyzzyaaac11=1,xyzzyaabv1
if(xyzzyaabt1==1)then
xyzzyaaaa11=xyzzyaaaa11+1
params(xyzzyaaaa11)=xyzzyaaco1
has_lolim(xyzzyaaaa11)=.true.
lolim(xyzzyaaaa11)=xyzzyaaaf11
endif
enddo
else
do xyzzyaaac11=1,xyzzyaabv1
do xyzzyaaab11=1,xyzzyaabq1
if(xyzzyaacc1(xyzzyaaab11,xyzzyaaac11)==1)then
do xyzzyaaad11=1,real1_complex2
xyzzyaaaa11=xyzzyaaaa11+1
params(xyzzyaaaa11)=xyzzyaacu1(xyzzyaaad11,xyzzyaaab11,xyzzyaaac11)
enddo
endif
if(xyzzyaacd1(xyzzyaaab11,xyzzyaaac11)==1)then
xyzzyaaaa11=xyzzyaaaa11+1
params(xyzzyaaaa11)=xyzzyaacv1(xyzzyaaab11,xyzzyaaac11)
has_lolim(xyzzyaaaa11)=.true.
lolim(xyzzyaaaa11)=xyzzyaaaf11
endif
enddo
enddo
endif
endif
if(xyzzyaadi1)then
do xyzzyaaac11=1,xyzzyaabv1
do xyzzyaaab11=1,xyzzyaabr1
if(xyzzyaacf1(xyzzyaaab11,xyzzyaaac11)==1)then
do xyzzyaaad11=1,real1_complex2
xyzzyaaaa11=xyzzyaaaa11+1
params(xyzzyaaaa11)=xyzzyaacw1(xyzzyaaad11,xyzzyaaab11,xyzzyaaac11)
enddo
endif
enddo
enddo
endif
if(xyzzyaadk1)then
if(xyzzyaabz1==1)then
xyzzyaaaa11=xyzzyaaaa11+1
params(xyzzyaaaa11)=xyzzyaacp1
has_lolim(xyzzyaaaa11)=.true.
lolim(xyzzyaaaa11)=xyzzyaaaf11
if(xyzzyaadn1)then
has_hilim(xyzzyaaaa11)=.true.
hilim(xyzzyaaaa11)=xyzzyaaag11
endif
is_shallow(xyzzyaaaa11)=.true.
is_redundant(xyzzyaaaa11)=all(xyzzyaadb1(:,:,:)==0.d0)
endif
do xyzzyaaac11=1,xyzzyaabv1
do xyzzyaaab11=0,xyzzyaabx1
if(xyzzyaach1(xyzzyaaab11,xyzzyaaac11)==1)then
do xyzzyaaad11=1,real1_complex2
xyzzyaaaa11=xyzzyaaaa11+1
params(xyzzyaaaa11)=xyzzyaadb1(xyzzyaaad11,xyzzyaaab11,xyzzyaaac11)
enddo
endif
enddo
enddo
endif
if(xyzzyaadl1)then
do xyzzyaaac11=1,xyzzyaabv1
do xyzzyaaab11=1,xyzzyaaca1
if(xyzzyaaci1(xyzzyaaab11,xyzzyaaac11)==1)then
do xyzzyaaad11=1,real1_complex2
xyzzyaaaa11=xyzzyaaaa11+1
params(xyzzyaaaa11)=xyzzyaadc1(xyzzyaaad11,xyzzyaaab11,xyzzyaaac11)
enddo
endif
if(xyzzyaacj1(xyzzyaaab11,xyzzyaaac11)==1)then
xyzzyaaaa11=xyzzyaaaa11+1
params(xyzzyaaaa11)=xyzzyaadd1(xyzzyaaab11,xyzzyaaac11)
has_lolim(xyzzyaaaa11)=.true.
lolim(xyzzyaaaa11)=xyzzyaaaf11
endif
if(xyzzyaack1(xyzzyaaab11,xyzzyaaac11)==1)then
xyzzyaaaa11=xyzzyaaaa11+1
params(xyzzyaaaa11)=xyzzyaade1(xyzzyaaab11,xyzzyaaac11)
has_lolim(xyzzyaaaa11)=.true.
lolim(xyzzyaaaa11)=0.d0
endif
enddo
do xyzzyaaab11=1,xyzzyaacb1
do xyzzyaaae11=1,dimensionality
xyzzyaadf1(:,xyzzyaaae11,xyzzyaaab11,xyzzyaaac11)=xyzzyaadf1(:,xyzzyaa&
&ae11,xyzzyaaab11,xyzzyaaac11)*xyzzyaada1(xyzzyaaac11)
if(xyzzyaacl1(xyzzyaaae11,xyzzyaaab11,xyzzyaaac11)==1)then
do xyzzyaaad11=1,real1_complex2
xyzzyaaaa11=xyzzyaaaa11+1
params(xyzzyaaaa11)=xyzzyaadf1(xyzzyaaad11,xyzzyaaae11,xyzzyaaab11,xyz&
&zyaaac11)
enddo
endif
enddo
if(xyzzyaacm1(xyzzyaaab11,xyzzyaaac11)==1)then
xyzzyaaaa11=xyzzyaaaa11+1
params(xyzzyaaaa11)=xyzzyaadg1(xyzzyaaab11,xyzzyaaac11)
has_lolim(xyzzyaaaa11)=.true.
lolim(xyzzyaaaa11)=xyzzyaaaf11
endif
if(xyzzyaacn1(xyzzyaaab11,xyzzyaaac11)==1)then
xyzzyaaaa11=xyzzyaaaa11+1
params(xyzzyaaaa11)=xyzzyaadh1(xyzzyaaab11,xyzzyaaac11)
has_lolim(xyzzyaaaa11)=.true.
lolim(xyzzyaaaa11)=0.d0
endif
enddo
enddo
endif
endif
if(xyzzyaabh1)then
if(type_wc==1.or.type_wc==3)then
do xyzzyaaac11=1,xyzzyaabw1
if(xyzzyaace1(xyzzyaaac11)==1)then
xyzzyaaaa11=xyzzyaaaa11+1
params(xyzzyaaaa11)=wc_gauss_exp(xyzzyaaac11)
has_lolim(xyzzyaaaa11)=.true.
lolim(xyzzyaaaa11)=xyzzyaaaf11
endif
enddo
endif
if(type_wc==2.or.type_wc==3)then
do xyzzyaaac11=1,xyzzyaabw1
do xyzzyaaab11=1,xyzzyaabs1
if(xyzzyaacg1(xyzzyaaab11,xyzzyaaac11)==1)then
xyzzyaaaa11=xyzzyaaaa11+1
params(xyzzyaaaa11)=xyzzyaacx1(xyzzyaaab11,xyzzyaaac11)
endif
enddo
enddo
endif
endif
end subroutine get_freeorb_params
subroutine put_freeorb_params(params,ignore,iparam_buffer,prestore,bad&
&_params)
implicit none
integer,intent(in) :: iparam_buffer
real(dp),intent(inout) :: params(xyzzyaabp1)
logical,intent(in) :: ignore(xyzzyaabp1),prestore
logical,intent(out) :: bad_params
integer xyzzyaaaa12,xyzzyaaab12,xyzzyaaac12,xyzzyaaad12,xyzzyaaae12
real(dp) xyzzyaaaf12,xyzzyaaag12,xyzzyaaah12,xyzzyaaai12
logical xyzzyaaaj12,xyzzyaaak12
bad_params=.false.
if(prestore)then
call xyzzyaafo1(iparam_buffer)
return
endif
xyzzyaaaj12=.false.
xyzzyaaak12=.false.
xyzzyaaae12=0
if(xyzzyaabg1)then
if(xyzzyaadj1)then
xyzzyaaag12=minval(xyzzyaacv1(:,:))
if(xyzzyaadm1)then
do xyzzyaaac12=1,xyzzyaabv1
if(xyzzyaabt1==1)then
xyzzyaaae12=xyzzyaaae12+1
if(.not.ignore(xyzzyaaae12))then
xyzzyaaco1=params(xyzzyaaae12)
xyzzyaaai12=1.d0/xyzzyaaco1
xyzzyaacv1(:,1)=(xyzzyaaai12*xyzzyaacy1(:))**2
endif
endif
enddo
else
do xyzzyaaac12=1,xyzzyaabv1
do xyzzyaaaa12=1,xyzzyaabq1
if(xyzzyaacc1(xyzzyaaaa12,xyzzyaaac12)==1)then
do xyzzyaaab12=1,real1_complex2
xyzzyaaae12=xyzzyaaae12+1
if(.not.ignore(xyzzyaaae12))xyzzyaacu1(xyzzyaaab12,xyzzyaaaa12,xyzzyaa&
&ac12)=params(xyzzyaaae12)
enddo
endif
if(xyzzyaacd1(xyzzyaaaa12,xyzzyaaac12)==1)then
xyzzyaaae12=xyzzyaaae12+1
if(.not.ignore(xyzzyaaae12))xyzzyaacv1(xyzzyaaaa12,xyzzyaaac12)=params&
&(xyzzyaaae12)
endif
enddo
enddo
endif
xyzzyaaah12=minval(xyzzyaacv1(:,:))
if(xyzzyaaag12/=xyzzyaaah12)xyzzyaaaj12=.true.
endif
if(xyzzyaadi1)then
do xyzzyaaac12=1,xyzzyaabv1
do xyzzyaaaa12=1,xyzzyaabr1
if(xyzzyaacf1(xyzzyaaaa12,xyzzyaaac12)==1)then
do xyzzyaaab12=1,real1_complex2
xyzzyaaae12=xyzzyaaae12+1
if(.not.ignore(xyzzyaaae12))xyzzyaacw1(xyzzyaaab12,xyzzyaaaa12,xyzzyaa&
&ac12)=params(xyzzyaaae12)
enddo
endif
enddo
enddo
endif
if(xyzzyaadk1)then
if(xyzzyaabz1==1)then
xyzzyaaae12=xyzzyaaae12+1
if(.not.ignore(xyzzyaaae12))then
xyzzyaacp1=params(xyzzyaaae12)
xyzzyaacq1=1.d0/xyzzyaacp1
xyzzyaacr1=xyzzyaacq1*xyzzyaacq1
xyzzyaaaj12=.true.
endif
endif
do xyzzyaaac12=1,xyzzyaabv1
do xyzzyaaaa12=0,xyzzyaabx1
if(xyzzyaach1(xyzzyaaaa12,xyzzyaaac12)==1)then
do xyzzyaaab12=1,real1_complex2
xyzzyaaae12=xyzzyaaae12+1
if(.not.ignore(xyzzyaaae12))xyzzyaadb1(xyzzyaaab12,xyzzyaaaa12,xyzzyaa&
&ac12)=params(xyzzyaaae12)
enddo
endif
enddo
xyzzyaadb1(:,1,xyzzyaaac12)=xyzzyaadb1(:,0,xyzzyaaac12)*real(xyzzyaaby&
&1,dp)*xyzzyaacq1
enddo
endif
if(xyzzyaadl1)then
do xyzzyaaac12=1,xyzzyaabv1
do xyzzyaaaa12=1,xyzzyaaca1
if(xyzzyaaci1(xyzzyaaaa12,xyzzyaaac12)==1)then
do xyzzyaaab12=1,real1_complex2
xyzzyaaae12=xyzzyaaae12+1
if(.not.ignore(xyzzyaaae12))xyzzyaadc1(xyzzyaaab12,xyzzyaaaa12,xyzzyaa&
&ac12)=params(xyzzyaaae12)
enddo
endif
if(xyzzyaacj1(xyzzyaaaa12,xyzzyaaac12)==1)then
xyzzyaaae12=xyzzyaaae12+1
if(.not.ignore(xyzzyaaae12))xyzzyaadd1(xyzzyaaaa12,xyzzyaaac12)=params&
&(xyzzyaaae12)
endif
if(xyzzyaack1(xyzzyaaaa12,xyzzyaaac12)==1)then
xyzzyaaae12=xyzzyaaae12+1
if(.not.ignore(xyzzyaaae12))xyzzyaade1(xyzzyaaaa12,xyzzyaaac12)=params&
&(xyzzyaaae12)
endif
enddo
do xyzzyaaaa12=1,xyzzyaacb1
do xyzzyaaad12=1,dimensionality
if(xyzzyaacl1(xyzzyaaad12,xyzzyaaaa12,xyzzyaaac12)==1)then
do xyzzyaaab12=1,real1_complex2
xyzzyaaae12=xyzzyaaae12+1
if(.not.ignore(xyzzyaaae12))xyzzyaadf1(xyzzyaaab12,xyzzyaaad12,xyzzyaa&
&aa12,xyzzyaaac12)=params(xyzzyaaae12)
enddo
endif
enddo
if(xyzzyaacm1(xyzzyaaaa12,xyzzyaaac12)==1)then
xyzzyaaae12=xyzzyaaae12+1
if(.not.ignore(xyzzyaaae12))xyzzyaadg1(xyzzyaaaa12,xyzzyaaac12)=params&
&(xyzzyaaae12)
endif
if(xyzzyaacn1(xyzzyaaaa12,xyzzyaaac12)==1)then
xyzzyaaae12=xyzzyaaae12+1
if(.not.ignore(xyzzyaaae12))xyzzyaadh1(xyzzyaaaa12,xyzzyaaac12)=params&
&(xyzzyaaae12)
endif
enddo
enddo
endif
do xyzzyaaac12=1,xyzzyaabv1
xyzzyaaaf12=0.d0
if(xyzzyaadj1)xyzzyaaaf12=xyzzyaaaf12+sum(xyzzyaacu1(:,:,xyzzyaaac12)*&
&*2)
if(xyzzyaadi1)xyzzyaaaf12=xyzzyaaaf12+sum(xyzzyaacw1(:,:,xyzzyaaac12)*&
&*2)
if(xyzzyaadk1)xyzzyaaaf12=xyzzyaaaf12+sum(xyzzyaadb1(:,:,xyzzyaaac12)*&
&*2)
if(xyzzyaadl1)xyzzyaaaf12=xyzzyaaaf12+sum(xyzzyaadc1(:,:,xyzzyaaac12)*&
&*2)+sum(xyzzyaadf1(:,:,:,xyzzyaaac12)**2)
if(xyzzyaaaf12<=0.d0)call errstop_master('PUT_HEG_PARAMS','Pairing orb&
&ital is zero.')
xyzzyaada1(xyzzyaaac12)=1.d0/sqrt(xyzzyaaaf12)
enddo
if(xyzzyaaaj12.and.iparam_buffer==0)call xyzzyaafz1
endif
if(xyzzyaabh1)then
if(type_wc==1.or.type_wc==3)then
do xyzzyaaac12=1,xyzzyaabw1
if(xyzzyaace1(xyzzyaaac12)==1)then
xyzzyaaae12=xyzzyaaae12+1
if(.not.ignore(xyzzyaaae12))then
wc_gauss_exp(xyzzyaaac12)=params(xyzzyaaae12)
xyzzyaaaj12=.true.
endif
endif
enddo
endif
if(type_wc==2.or.type_wc==3)then
do xyzzyaaac12=1,xyzzyaabw1
do xyzzyaaaa12=1,xyzzyaabs1
if(xyzzyaacg1(xyzzyaaaa12,xyzzyaaac12)==1)then
xyzzyaaae12=xyzzyaaae12+1
if(.not.ignore(xyzzyaaae12))xyzzyaacx1(xyzzyaaaa12,xyzzyaaac12)=params&
&(xyzzyaaae12)
endif
enddo
enddo
endif
if(xyzzyaaak12.and.iparam_buffer==0)call xyzzyaafy1
endif
call xyzzyaafn1(iparam_buffer)
end subroutine put_freeorb_params
subroutine xyzzyaafn1(indx)
implicit none
integer,intent(in) :: indx
if(xyzzyaabg1)then
if(xyzzyaadj1)then
call dcopy(xyzzyaaeg1,xyzzyaacu1(1,1,1),1,xyzzyaaep1(1,1,1,indx),1)
call dcopy(xyzzyaaeh1,xyzzyaacv1(1,1),1,xyzzyaaeq1(1,1,indx),1)
endif
if(xyzzyaadi1)then
call dcopy(xyzzyaaei1,xyzzyaacw1(1,1,1),1,xyzzyaaer1(1,1,1,indx),1)
endif
if(xyzzyaadk1)then
xyzzyaaes1(indx)=xyzzyaacp1
xyzzyaaet1(indx)=xyzzyaacq1
xyzzyaaeu1(indx)=xyzzyaacr1
call dcopy(xyzzyaaej1,xyzzyaadb1(1,0,1),1,xyzzyaaev1(1,0,1,indx),1)
endif
if(xyzzyaadl1)then
if(xyzzyaaca1>0)then
call dcopy(xyzzyaaek1,xyzzyaadc1(1,1,1),1,xyzzyaaew1(1,1,1,indx),1)
call dcopy(xyzzyaael1,xyzzyaadd1(1,1),1,xyzzyaaex1(1,1,indx),1)
call dcopy(xyzzyaael1,xyzzyaade1(1,1),1,xyzzyaaey1(1,1,indx),1)
endif
if(xyzzyaacb1>0)then
call dcopy(xyzzyaaem1,xyzzyaadf1(1,1,1,1),1,xyzzyaaez1(1,1,1,1,indx),1&
&)
call dcopy(xyzzyaaen1,xyzzyaadg1(1,1),1,xyzzyaafa1(1,1,indx),1)
call dcopy(xyzzyaaen1,xyzzyaadh1(1,1),1,xyzzyaafb1(1,1,indx),1)
endif
endif
call dcopy(xyzzyaabv1,xyzzyaada1(1),1,xyzzyaafc1(1,indx),1)
endif
if(xyzzyaabh1)then
if(type_wc==1.or.type_wc==3)then
call dcopy(xyzzyaabw1,wc_gauss_exp(1),1,xyzzyaafd1(1,indx),1)
endif
if(type_wc==2.or.type_wc==3)then
call dcopy(xyzzyaaeo1,xyzzyaacx1(1,1),1,xyzzyaafe1(1,1,indx),1)
endif
endif
end subroutine xyzzyaafn1
subroutine xyzzyaafo1(indx)
implicit none
integer,intent(in) :: indx
if(xyzzyaabg1)then
if(xyzzyaadj1)then
call dcopy(xyzzyaaeg1,xyzzyaaep1(1,1,1,indx),1,xyzzyaacu1(1,1,1),1)
call dcopy(xyzzyaaeh1,xyzzyaaeq1(1,1,indx),1,xyzzyaacv1(1,1),1)
endif
if(xyzzyaadi1)then
call dcopy(xyzzyaaei1,xyzzyaaer1(1,1,1,indx),1,xyzzyaacw1(1,1,1),1)
endif
if(xyzzyaadk1)then
xyzzyaacp1=xyzzyaaes1(indx)
xyzzyaacq1=xyzzyaaet1(indx)
xyzzyaacr1=xyzzyaaeu1(indx)
call dcopy(xyzzyaaej1,xyzzyaaev1(1,0,1,indx),1,xyzzyaadb1(1,0,1),1)
endif
if(xyzzyaadl1)then
if(xyzzyaaca1>0)then
call dcopy(xyzzyaaek1,xyzzyaaew1(1,1,1,indx),1,xyzzyaadc1(1,1,1),1)
call dcopy(xyzzyaael1,xyzzyaaex1(1,1,indx),1,xyzzyaadd1(1,1),1)
call dcopy(xyzzyaael1,xyzzyaaey1(1,1,indx),1,xyzzyaade1(1,1),1)
endif
if(xyzzyaacb1>0)then
call dcopy(xyzzyaaem1,xyzzyaaez1(1,1,1,1,indx),1,xyzzyaadf1(1,1,1,1),1&
&)
call dcopy(xyzzyaaen1,xyzzyaafa1(1,1,indx),1,xyzzyaadg1(1,1),1)
call dcopy(xyzzyaaen1,xyzzyaafb1(1,1,indx),1,xyzzyaadh1(1,1),1)
endif
endif
call dcopy(xyzzyaabv1,xyzzyaafc1(1,indx),1,xyzzyaada1(1),1)
endif
if(xyzzyaabh1)then
if(type_wc==1.or.type_wc==3)then
call dcopy(xyzzyaabw1,xyzzyaafd1(1,indx),1,wc_gauss_exp(1),1)
endif
if(type_wc==2.or.type_wc==3)then
call dcopy(xyzzyaaeo1,xyzzyaafe1(1,1,indx),1,xyzzyaacx1(1,1),1)
endif
endif
end subroutine xyzzyaafo1
subroutine free_fluid_orb_eval(rvec,jspin,lnorb,norb,orbmask,val,fsd,o&
&rbval,orbgrad,orblap,orbsderivs)
implicit none
integer,intent(in) :: jspin,lnorb,norb
real(dp),intent(in) :: rvec(3)
real(dp),intent(inout) :: orbval(lnorb,*),orbgrad(3,lnorb,*),orblap(ln&
&orb,*)
real(dp),intent(inout),optional :: orbsderivs(6,lnorb,*)
logical,intent(in) :: val,fsd,orbmask(*)
integer xyzzyaaaa15,xyzzyaaab15,xyzzyaaac15
real(dp) xyzzyaaad15,xyzzyaaae15,xyzzyaaaf15,xyzzyaaag15,xyzzyaaah15
complex(dp) xyzzyaaai15,xyzzyaaaj15,xyzzyaaak15,xyzzyaaal15,xyzzyaaam1&
&5
if(heg_nele(jspin)==0)return
if(periodicity==0)then
orbval(:,1)=1.d0
orbgrad(:,:,1)=0.d0
orblap(:,1)=0.d0
if(present(orbsderivs))orbsderivs(:,:,1)=0.d0
return
endif
xyzzyaaah15=dot_product(b1,rvec)
xyzzyaaaj15=cmplx(cos(xyzzyaaah15),sin(xyzzyaaah15),dp)
xyzzyaaah15=dot_product(b2,rvec)
xyzzyaaak15=cmplx(cos(xyzzyaaah15),sin(xyzzyaaah15),dp)
xyzzyaaah15=dot_product(b3,rvec)
xyzzyaaal15=cmplx(cos(xyzzyaaah15),sin(xyzzyaaah15),dp)
xyzzyaabc1(0)=cmplx(xyzzyaaar1,0.d0,dp)
xyzzyaabd1(0)=xyzzyaabc1(0)
xyzzyaabe1(0)=xyzzyaabc1(0)
if(xyzzyaabl1)then
do xyzzyaaab15=1,xyzzyaaae1(1)
xyzzyaabc1(xyzzyaaab15)=xyzzyaabc1(xyzzyaaab15-1)*xyzzyaaaj15
xyzzyaabc1(-xyzzyaaab15)=conjg(xyzzyaabc1(xyzzyaaab15))
enddo
else
do xyzzyaaab15=1,xyzzyaaae1(1)
xyzzyaabc1(xyzzyaaab15)=xyzzyaabc1(xyzzyaaab15-1)*xyzzyaaaj15
enddo
endif
do xyzzyaaab15=1,xyzzyaaae1(2)
xyzzyaabd1(xyzzyaaab15)=xyzzyaabd1(xyzzyaaab15-1)*xyzzyaaak15
xyzzyaabd1(-xyzzyaaab15)=conjg(xyzzyaabd1(xyzzyaaab15))
enddo
do xyzzyaaab15=1,xyzzyaaae1(3)
xyzzyaabe1(xyzzyaaab15)=xyzzyaabe1(xyzzyaaab15-1)*xyzzyaaal15
xyzzyaabe1(-xyzzyaaab15)=conjg(xyzzyaabe1(xyzzyaaab15))
enddo
if(xyzzyaabl1)then
xyzzyaaah15=dot_product(k_offset,rvec)
xyzzyaaam15=cmplx(cos(xyzzyaaah15),sin(xyzzyaaah15),dp)
do xyzzyaaab15=-xyzzyaaae1(3),xyzzyaaae1(3)
xyzzyaabe1(xyzzyaaab15)=xyzzyaabe1(xyzzyaaab15)*xyzzyaaam15
enddo
xyzzyaaac15=0
if(.not.fsd)then
do xyzzyaaab15=1,xyzzyaaai1(jspin)
xyzzyaaac15=xyzzyaaac15+1
if(.not.orbmask(xyzzyaaac15))cycle
xyzzyaaai15=xyzzyaabc1(xyzzyaaam1(1,xyzzyaaab15))*xyzzyaabd1(xyzzyaaam&
&1(2,xyzzyaaab15))*xyzzyaabe1(xyzzyaaam1(3,xyzzyaaab15))
orbval(xyzzyaaac15,1)=dble(xyzzyaaai15)
orbval(xyzzyaaac15,2)=aimag(xyzzyaaai15)
enddo
else
if(.not.val)then
if(present(orbsderivs))then
do xyzzyaaab15=1,xyzzyaaai1(jspin)
xyzzyaaac15=xyzzyaaac15+1
if(.not.orbmask(xyzzyaaac15))cycle
xyzzyaaai15=xyzzyaabc1(xyzzyaaam1(1,xyzzyaaab15))*xyzzyaabd1(xyzzyaaam&
&1(2,xyzzyaaab15))*xyzzyaabe1(xyzzyaaam1(3,xyzzyaaab15))
xyzzyaaaf15=dble(xyzzyaaai15)
xyzzyaaae15=-aimag(xyzzyaaai15)
xyzzyaaag15=-xyzzyaaaf15
orbgrad(1:3,xyzzyaaac15,1)=gs_kvec(1:3,xyzzyaaab15)*xyzzyaaae15
orbgrad(1:3,xyzzyaaac15,2)=gs_kvec(1:3,xyzzyaaab15)*xyzzyaaaf15
orblap(xyzzyaaac15,1)=xyzzyaaaw1(xyzzyaaab15)*xyzzyaaag15
orblap(xyzzyaaac15,2)=xyzzyaaaw1(xyzzyaaab15)*xyzzyaaae15
orbsderivs(1:6,xyzzyaaac15,1)=xyzzyaabb1(1:6,xyzzyaaab15)*xyzzyaaag15
orbsderivs(1:6,xyzzyaaac15,2)=xyzzyaabb1(1:6,xyzzyaaab15)*xyzzyaaae15
enddo
else
do xyzzyaaab15=1,xyzzyaaai1(jspin)
xyzzyaaac15=xyzzyaaac15+1
if(.not.orbmask(xyzzyaaac15))cycle
xyzzyaaai15=xyzzyaabc1(xyzzyaaam1(1,xyzzyaaab15))*xyzzyaabd1(xyzzyaaam&
&1(2,xyzzyaaab15))*xyzzyaabe1(xyzzyaaam1(3,xyzzyaaab15))
xyzzyaaaf15=dble(xyzzyaaai15)
xyzzyaaae15=-aimag(xyzzyaaai15)
orbgrad(1:3,xyzzyaaac15,1)=gs_kvec(1:3,xyzzyaaab15)*xyzzyaaae15
orbgrad(1:3,xyzzyaaac15,2)=gs_kvec(1:3,xyzzyaaab15)*xyzzyaaaf15
orblap(xyzzyaaac15,1)=-xyzzyaaaw1(xyzzyaaab15)*xyzzyaaaf15
orblap(xyzzyaaac15,2)=xyzzyaaaw1(xyzzyaaab15)*xyzzyaaae15
enddo
endif
else
if(present(orbsderivs))then
do xyzzyaaab15=1,xyzzyaaai1(jspin)
xyzzyaaac15=xyzzyaaac15+1
if(.not.orbmask(xyzzyaaac15))cycle
xyzzyaaai15=xyzzyaabc1(xyzzyaaam1(1,xyzzyaaab15))*xyzzyaabd1(xyzzyaaam&
&1(2,xyzzyaaab15))*xyzzyaabe1(xyzzyaaam1(3,xyzzyaaab15))
xyzzyaaaf15=dble(xyzzyaaai15)
xyzzyaaad15=aimag(xyzzyaaai15)
xyzzyaaag15=-xyzzyaaaf15
xyzzyaaae15=-xyzzyaaad15
orbval(xyzzyaaac15,1)=xyzzyaaaf15
orbval(xyzzyaaac15,2)=xyzzyaaad15
orbgrad(1:3,xyzzyaaac15,1)=gs_kvec(1:3,xyzzyaaab15)*xyzzyaaae15
orbgrad(1:3,xyzzyaaac15,2)=gs_kvec(1:3,xyzzyaaab15)*xyzzyaaaf15
orblap(xyzzyaaac15,1)=xyzzyaaaw1(xyzzyaaab15)*xyzzyaaag15
orblap(xyzzyaaac15,2)=xyzzyaaaw1(xyzzyaaab15)*xyzzyaaae15
orbsderivs(1:6,xyzzyaaac15,1)=xyzzyaabb1(1:6,xyzzyaaab15)*xyzzyaaag15
orbsderivs(1:6,xyzzyaaac15,2)=xyzzyaabb1(1:6,xyzzyaaab15)*xyzzyaaae15
enddo
else
do xyzzyaaab15=1,xyzzyaaai1(jspin)
xyzzyaaac15=xyzzyaaac15+1
if(.not.orbmask(xyzzyaaac15))cycle
xyzzyaaai15=xyzzyaabc1(xyzzyaaam1(1,xyzzyaaab15))*xyzzyaabd1(xyzzyaaam&
&1(2,xyzzyaaab15))*xyzzyaabe1(xyzzyaaam1(3,xyzzyaaab15))
xyzzyaaaf15=dble(xyzzyaaai15)
xyzzyaaad15=aimag(xyzzyaaai15)
xyzzyaaae15=-aimag(xyzzyaaai15)
orbval(xyzzyaaac15,1)=xyzzyaaaf15
orbval(xyzzyaaac15,2)=xyzzyaaad15
orbgrad(1:3,xyzzyaaac15,1)=gs_kvec(1:3,xyzzyaaab15)*xyzzyaaae15
orbgrad(1:3,xyzzyaaac15,2)=gs_kvec(1:3,xyzzyaaab15)*xyzzyaaaf15
orblap(xyzzyaaac15,1)=-xyzzyaaaw1(xyzzyaaab15)*xyzzyaaaf15
orblap(xyzzyaaac15,2)=xyzzyaaaw1(xyzzyaaab15)*xyzzyaaae15
enddo
endif
endif
endif
else
xyzzyaaac15=1
if(.not.fsd)then
if(orbmask(xyzzyaaac15))orbval(xyzzyaaac15,1)=xyzzyaaaq1*one_over_root&
&_two
do xyzzyaaab15=2,xyzzyaaai1(jspin)
if(.not.orbmask(xyzzyaaac15+1).and..not.orbmask(xyzzyaaac15+2))then
xyzzyaaac15=xyzzyaaac15+2
cycle
endif
xyzzyaaai15=xyzzyaabc1(xyzzyaaam1(1,xyzzyaaab15))*xyzzyaabd1(xyzzyaaam&
&1(2,xyzzyaaab15))*xyzzyaabe1(xyzzyaaam1(3,xyzzyaaab15))
xyzzyaaac15=xyzzyaaac15+1
if(orbmask(xyzzyaaac15))orbval(xyzzyaaac15,1)=aimag(xyzzyaaai15)
xyzzyaaac15=xyzzyaaac15+1
if(orbmask(xyzzyaaac15))orbval(xyzzyaaac15,1)=dble(xyzzyaaai15)
enddo
else
if(orbmask(xyzzyaaac15))then
orbgrad(1:3,xyzzyaaac15,1)=0.d0
orblap(xyzzyaaac15,1)=0.d0
endif
if(.not.val)then
if(present(orbsderivs))then
if(orbmask(xyzzyaaac15))orbsderivs(1:6,xyzzyaaac15,1)=0.d0
do xyzzyaaab15=2,xyzzyaaai1(jspin)
if(.not.orbmask(xyzzyaaac15+1).and..not.orbmask(xyzzyaaac15+2))then
xyzzyaaac15=xyzzyaaac15+2
cycle
endif
xyzzyaaai15=xyzzyaabc1(xyzzyaaam1(1,xyzzyaaab15))*xyzzyaabd1(xyzzyaaam&
&1(2,xyzzyaaab15))*xyzzyaabe1(xyzzyaaam1(3,xyzzyaaab15))
xyzzyaaaf15=dble(xyzzyaaai15)
xyzzyaaae15=-aimag(xyzzyaaai15)
xyzzyaaag15=-xyzzyaaaf15
xyzzyaaac15=xyzzyaaac15+1
if(orbmask(xyzzyaaac15))then
orbgrad(1:3,xyzzyaaac15,1)=gs_kvec(1:3,xyzzyaaab15)*xyzzyaaaf15
orblap(xyzzyaaac15,1)=xyzzyaaaw1(xyzzyaaab15)*xyzzyaaae15
orbsderivs(1:6,xyzzyaaac15,1)=xyzzyaabb1(1:6,xyzzyaaab15)*xyzzyaaae15
endif
xyzzyaaac15=xyzzyaaac15+1
if(orbmask(xyzzyaaac15))then
orbgrad(1:3,xyzzyaaac15,1)=gs_kvec(1:3,xyzzyaaab15)*xyzzyaaae15
orblap(xyzzyaaac15,1)=xyzzyaaaw1(xyzzyaaab15)*xyzzyaaag15
orbsderivs(1:6,xyzzyaaac15,1)=xyzzyaabb1(1:6,xyzzyaaab15)*xyzzyaaag15
endif
enddo
else
do xyzzyaaab15=2,xyzzyaaai1(jspin)
if(.not.orbmask(xyzzyaaac15+1).and..not.orbmask(xyzzyaaac15+2))then
xyzzyaaac15=xyzzyaaac15+2
cycle
endif
xyzzyaaai15=xyzzyaabc1(xyzzyaaam1(1,xyzzyaaab15))*xyzzyaabd1(xyzzyaaam&
&1(2,xyzzyaaab15))*xyzzyaabe1(xyzzyaaam1(3,xyzzyaaab15))
xyzzyaaaf15=dble(xyzzyaaai15)
xyzzyaaae15=-aimag(xyzzyaaai15)
xyzzyaaac15=xyzzyaaac15+1
if(orbmask(xyzzyaaac15))then
orbgrad(1:3,xyzzyaaac15,1)=gs_kvec(1:3,xyzzyaaab15)*xyzzyaaaf15
orblap(xyzzyaaac15,1)=xyzzyaaaw1(xyzzyaaab15)*xyzzyaaae15
endif
xyzzyaaac15=xyzzyaaac15+1
if(orbmask(xyzzyaaac15))then
orbgrad(1:3,xyzzyaaac15,1)=gs_kvec(1:3,xyzzyaaab15)*xyzzyaaae15
orblap(xyzzyaaac15,1)=-xyzzyaaaw1(xyzzyaaab15)*xyzzyaaaf15
endif
enddo
endif
else
orbval(xyzzyaaac15,1)=xyzzyaaaq1*one_over_root_two
if(present(orbsderivs))then
orbsderivs(:,xyzzyaaac15,1)=0.d0
do xyzzyaaab15=2,xyzzyaaai1(jspin)
if(.not.orbmask(xyzzyaaac15+1).and..not.orbmask(xyzzyaaac15+2))then
xyzzyaaac15=xyzzyaaac15+2
cycle
endif
xyzzyaaai15=xyzzyaabc1(xyzzyaaam1(1,xyzzyaaab15))*xyzzyaabd1(xyzzyaaam&
&1(2,xyzzyaaab15))*xyzzyaabe1(xyzzyaaam1(3,xyzzyaaab15))
xyzzyaaaf15=dble(xyzzyaaai15)
xyzzyaaad15=aimag(xyzzyaaai15)
xyzzyaaag15=-xyzzyaaaf15
xyzzyaaae15=-xyzzyaaad15
xyzzyaaac15=xyzzyaaac15+1
if(orbmask(xyzzyaaac15))then
orbval(xyzzyaaac15,1)=xyzzyaaad15
orbgrad(1:3,xyzzyaaac15,1)=gs_kvec(1:3,xyzzyaaab15)*xyzzyaaaf15
orblap(xyzzyaaac15,1)=xyzzyaaaw1(xyzzyaaab15)*xyzzyaaae15
orbsderivs(1:6,xyzzyaaac15,1)=xyzzyaabb1(1:6,xyzzyaaab15)*xyzzyaaae15
endif
xyzzyaaac15=xyzzyaaac15+1
if(orbmask(xyzzyaaac15))then
orbval(xyzzyaaac15,1)=xyzzyaaaf15
orbgrad(1:3,xyzzyaaac15,1)=gs_kvec(1:3,xyzzyaaab15)*xyzzyaaae15
orblap(xyzzyaaac15,1)=xyzzyaaaw1(xyzzyaaab15)*xyzzyaaag15
orbsderivs(1:6,xyzzyaaac15,1)=xyzzyaabb1(1:6,xyzzyaaab15)*xyzzyaaag15
endif
enddo
else
do xyzzyaaab15=2,xyzzyaaai1(jspin)
if(.not.orbmask(xyzzyaaac15+1).and..not.orbmask(xyzzyaaac15+2))then
xyzzyaaac15=xyzzyaaac15+2
cycle
endif
xyzzyaaai15=xyzzyaabc1(xyzzyaaam1(1,xyzzyaaab15))*xyzzyaabd1(xyzzyaaam&
&1(2,xyzzyaaab15))*xyzzyaabe1(xyzzyaaam1(3,xyzzyaaab15))
xyzzyaaaf15=dble(xyzzyaaai15)
xyzzyaaad15=aimag(xyzzyaaai15)
xyzzyaaae15=-xyzzyaaad15
xyzzyaaac15=xyzzyaaac15+1
if(orbmask(xyzzyaaac15))then
orbval(xyzzyaaac15,1)=xyzzyaaad15
orbgrad(1:3,xyzzyaaac15,1)=gs_kvec(1:3,xyzzyaaab15)*xyzzyaaaf15
orblap(xyzzyaaac15,1)=xyzzyaaaw1(xyzzyaaab15)*xyzzyaaae15
endif
xyzzyaaac15=xyzzyaaac15+1
if(orbmask(xyzzyaaac15))then
orbval(xyzzyaaac15,1)=xyzzyaaaf15
orbgrad(1:3,xyzzyaaac15,1)=gs_kvec(1:3,xyzzyaaab15)*xyzzyaaae15
orblap(xyzzyaaac15,1)=-xyzzyaaaw1(xyzzyaaab15)*xyzzyaaaf15
endif
enddo
endif
endif
endif
endif
if(xyzzyaadx1)then
do xyzzyaaaa15=1,xyzzyaadp1
xyzzyaaac15=xyzzyaado1+xyzzyaaaa15
if(.not.orbmask(xyzzyaaac15))cycle
xyzzyaaah15=dot_product(xyzzyaadt1(1:3,xyzzyaaaa15),rvec)
xyzzyaaaf15=xyzzyaaaq1*cos(xyzzyaaah15)
xyzzyaaad15=xyzzyaaaq1*sin(xyzzyaaah15)
xyzzyaaag15=-xyzzyaaaf15
xyzzyaaae15=-xyzzyaaad15
if(xyzzyaabl1)then
if(val)then
orbval(xyzzyaaac15,1)=xyzzyaaaf15
orbval(xyzzyaaac15,2)=xyzzyaaad15
endif
if(fsd)then
orbgrad(1:3,xyzzyaaac15,1)=xyzzyaadt1(1:3,xyzzyaaaa15)*xyzzyaaae15
orbgrad(1:3,xyzzyaaac15,2)=xyzzyaadt1(1:3,xyzzyaaaa15)*xyzzyaaaf15
orblap(xyzzyaaac15,1)=xyzzyaadu1(xyzzyaaaa15)*xyzzyaaag15
orblap(xyzzyaaac15,2)=xyzzyaadu1(xyzzyaaaa15)*xyzzyaaae15
if(present(orbsderivs))then
orbsderivs(1:6,xyzzyaaac15,1)=xyzzyaadv1(1:6,xyzzyaaaa15)*xyzzyaaag15
orbsderivs(1:6,xyzzyaaac15,2)=xyzzyaadv1(1:6,xyzzyaaaa15)*xyzzyaaae15
endif
endif
elseif(xyzzyaady1(xyzzyaaaa15))then
if(val)orbval(xyzzyaaac15,1)=xyzzyaaaf15
if(fsd)then
orbgrad(1:3,xyzzyaaac15,1)=xyzzyaadt1(1:3,xyzzyaaaa15)*xyzzyaaae15
orblap(xyzzyaaac15,1)=xyzzyaadu1(xyzzyaaaa15)*xyzzyaaag15
if(present(orbsderivs))orbsderivs(1:6,xyzzyaaac15,1)=xyzzyaadv1(1:6,xy&
&zzyaaaa15)*xyzzyaaag15
endif
else
if(val)orbval(xyzzyaaac15,1)=xyzzyaaad15
if(fsd)then
orbgrad(1:3,xyzzyaaac15,1)=xyzzyaadt1(1:3,xyzzyaaaa15)*xyzzyaaaf15
orblap(xyzzyaaac15,1)=xyzzyaadu1(xyzzyaaaa15)*xyzzyaaae15
if(present(orbsderivs))orbsderivs(1:6,xyzzyaaac15,1)=xyzzyaadv1(1:6,xy&
&zzyaaaa15)*xyzzyaaae15
endif
endif
enddo
endif
if(complex_wf.and..not.xyzzyaabl1)then
if(val)orbval(1:norb,2)=0.d0
if(fsd)then
orbgrad(1:3,1:norb,2)=0.d0
orblap(1:norb,2)=0.d0
if(present(orbsderivs))orbsderivs(1:6,1:norb,2)=0.d0
endif
endif
end subroutine free_fluid_orb_eval
subroutine free_crystal_orb_eval(rvec,jspin,lnorb,norb,orbmask,val,fsd&
&,orbval,orbgrad,orblap,orbsderivs)
implicit none
integer,intent(in) :: jspin,lnorb,norb
real(dp),intent(in) :: rvec(3)
real(dp),intent(inout) :: orbval(lnorb,*),orbgrad(3,lnorb,*),orblap(ln&
&orb,*)
real(dp),intent(inout),optional :: orbsderivs(6,lnorb,*)
logical,intent(in) :: val,fsd,orbmask(*)
integer xyzzyaaaa16,xyzzyaaab16,xyzzyaaac16,xyzzyaaad16,xyzzyaaae16,xy&
&zzyaaaf16,xyzzyaaag16,xyzzyaaah16,xyzzyaaai16,xyzzyaaaj16,xyzzyaaak16
real(dp) xyzzyaaal16,xyzzyaaam16,xyzzyaaan16,xyzzyaaao16,xyzzyaaap16(3&
&),xyzzyaaaq16(3),xyzzyaaar16,xyzzyaaas16,xyzzyaaat16,xyzzyaaau16(6),x&
&yzzyaaav16,xyzzyaaaw16,xyzzyaaax16(3)
xyzzyaaai16=maxval(heg_slatt)
xyzzyaaak16=maxval(heg_nele)
xyzzyaaag16=which_ssingle(jspin,spin_dep_wc)
if(isperiodic)then
if(dimensionality==3)then
xyzzyaaaq16=matmul(rvec,ainv)
xyzzyaaaq16=modulo(xyzzyaaaq16,1.d0)
xyzzyaaax16=matmul(xyzzyaaaq16,amat)
elseif(dimensionality==2)then
xyzzyaaaq16(1:2)=matmul(rvec(1:2),ainv(1:2,1:2))
xyzzyaaaq16(1:2)=modulo(xyzzyaaaq16(1:2),1.d0)
xyzzyaaax16(1:2)=matmul(xyzzyaaaq16(1:2),amat(1:2,1:2))
xyzzyaaax16(3)=0.d0
else
xyzzyaaax16(1)=modulo(rvec(1),amat(1,1))
xyzzyaaax16(2:3)=0.d0
endif
else
xyzzyaaax16=rvec
endif
if(val)orbval(1:norb,1:real1_complex2)=0.d0
if(fsd)then
orbgrad(1:3,1:norb,1:real1_complex2)=0.d0
orblap(1:norb,1:real1_complex2)=0.d0
if(present(orbsderivs))orbsderivs(1:6,1:norb,1:real1_complex2)=0.d0
endif
if(type_wc==1.or.type_wc==3)then
xyzzyaaat16=wc_gauss_exp(xyzzyaaag16)
do xyzzyaaaj16=1,xyzzyaaai16
if(.not.any(orbmask((xyzzyaaaj16-1)*xyzzyaaak16+1:xyzzyaaaj16*xyzzyaaa&
&k16)))cycle
do xyzzyaaab16=1,xyzzyaaak1(xyzzyaaaj16)
do xyzzyaaaa16=1,xyzzyaaal1(xyzzyaaab16,xyzzyaaaj16)
xyzzyaaam16=xyzzyaaax16(1)-site_pos_in_cell(1,xyzzyaaaa16,xyzzyaaab16,&
&xyzzyaaaj16)
if(dimensionality==3)then
xyzzyaaan16=xyzzyaaax16(2)-site_pos_in_cell(2,xyzzyaaaa16,xyzzyaaab16,&
&xyzzyaaaj16)
xyzzyaaao16=xyzzyaaax16(3)-site_pos_in_cell(3,xyzzyaaaa16,xyzzyaaab16,&
&xyzzyaaaj16)
xyzzyaaal16=xyzzyaaam16*xyzzyaaam16+xyzzyaaan16*xyzzyaaan16+xyzzyaaao1&
&6*xyzzyaaao16
elseif(dimensionality==2)then
xyzzyaaan16=xyzzyaaax16(2)-site_pos_in_cell(2,xyzzyaaaa16,xyzzyaaab16,&
&xyzzyaaaj16)
xyzzyaaal16=xyzzyaaam16*xyzzyaaam16+xyzzyaaan16*xyzzyaaan16
else
xyzzyaaal16=xyzzyaaam16*xyzzyaaam16
endif
xyzzyaaar16=-xyzzyaaat16*xyzzyaaal16
if(xyzzyaaar16>xyzzyaaav1)then
xyzzyaaar16=xyzzyaaas1*exp(xyzzyaaar16)
xyzzyaaah16=xyzzyaaaj1(xyzzyaaaa16,xyzzyaaab16,xyzzyaaaj16)
if(val)orbval(xyzzyaaah16,1)=orbval(xyzzyaaah16,1)+xyzzyaaar16
if(fsd)then
xyzzyaaaw16=xyzzyaaat16+xyzzyaaat16
xyzzyaaas16=xyzzyaaaw16*xyzzyaaar16
orbgrad(1,xyzzyaaah16,1)=orbgrad(1,xyzzyaaah16,1)-xyzzyaaam16*xyzzyaaa&
&s16
if(dimensionality==3)then
orbgrad(3,xyzzyaaah16,1)=orbgrad(3,xyzzyaaah16,1)-xyzzyaaao16*xyzzyaaa&
&s16
orbgrad(2,xyzzyaaah16,1)=orbgrad(2,xyzzyaaah16,1)-xyzzyaaan16*xyzzyaaa&
&s16
orblap(xyzzyaaah16,1)=orblap(xyzzyaaah16,1)+xyzzyaaas16*(-3.d0+xyzzyaa&
&aw16*xyzzyaaal16)
elseif(dimensionality==2)then
orbgrad(2,xyzzyaaah16,1)=orbgrad(2,xyzzyaaah16,1)-xyzzyaaan16*xyzzyaaa&
&s16
orblap(xyzzyaaah16,1)=orblap(xyzzyaaah16,1)+xyzzyaaas16*(-2.d0+xyzzyaa&
&aw16*xyzzyaaal16)
else
orblap(xyzzyaaah16,1)=orblap(xyzzyaaah16,1)+xyzzyaaas16*(-1.d0+xyzzyaa&
&aw16*xyzzyaaal16)
endif
if(present(orbsderivs))then
orbsderivs(1,xyzzyaaah16,1)=orbsderivs(1,xyzzyaaah16,1)+xyzzyaaas16*(x&
&yzzyaaaw16*xyzzyaaam16*xyzzyaaam16-1.d0)
if(dimensionality>1)then
orbsderivs(2,xyzzyaaah16,1)=orbsderivs(2,xyzzyaaah16,1)+xyzzyaaas16*(x&
&yzzyaaaw16*xyzzyaaan16*xyzzyaaan16-1.d0)
orbsderivs(4,xyzzyaaah16,1)=orbsderivs(4,xyzzyaaah16,1)+xyzzyaaas16*xy&
&zzyaaaw16*xyzzyaaam16*xyzzyaaan16
if(dimensionality==3)then
orbsderivs(3,xyzzyaaah16,1)=orbsderivs(3,xyzzyaaah16,1)+xyzzyaaas16*(x&
&yzzyaaaw16*xyzzyaaao16*xyzzyaaao16-1.d0)
orbsderivs(5,xyzzyaaah16,1)=orbsderivs(5,xyzzyaaah16,1)+xyzzyaaas16*xy&
&zzyaaaw16*xyzzyaaam16*xyzzyaaao16
orbsderivs(6,xyzzyaaah16,1)=orbsderivs(6,xyzzyaaah16,1)+xyzzyaaas16*xy&
&zzyaaaw16*xyzzyaaan16*xyzzyaaao16
endif
endif
endif
endif
endif
enddo
enddo
enddo
endif
if(type_wc==2.or.type_wc==3)then
do xyzzyaaaj16=1,xyzzyaaai16
if(.not.any(orbmask((xyzzyaaaj16-1)*xyzzyaaak16+1:xyzzyaaaj16*xyzzyaaa&
&k16)))cycle
if(val.and..not.fsd)then
do xyzzyaaah16=1,heg_nele(jspin)
orbval(xyzzyaaah16,1)=orbval(xyzzyaaah16,1)+xyzzyaacx1(1,xyzzyaaag16)
xyzzyaaaq16(1:periodicity)=xyzzyaaax16(1:periodicity) -site_pos_in_cel&
&l(1:periodicity,xyzzyaaah16,1,xyzzyaaaj16)
xyzzyaaad16=1
do xyzzyaaae16=2,xyzzyaabs1
xyzzyaaac16=xyzzyaaad16+1
xyzzyaaad16=xyzzyaaan1(xyzzyaaae16)
xyzzyaaas16=0.d0
do xyzzyaaaf16=xyzzyaaac16,xyzzyaaad16
xyzzyaaas16=xyzzyaaas16+cos(sum(xyzzyaaaz1(1:periodicity,xyzzyaaaf16)*&
&xyzzyaaaq16(1:periodicity)))
enddo
orbval(xyzzyaaah16,1)=orbval(xyzzyaaah16,1)+xyzzyaacx1(xyzzyaaae16,xyz&
&zyaaag16)*xyzzyaaas16
enddo
enddo
elseif(fsd.and..not.val)then
do xyzzyaaah16=1,heg_nele(jspin)
xyzzyaaaq16(1:periodicity)=xyzzyaaax16(1:periodicity) -site_pos_in_cel&
&l(1:periodicity,xyzzyaaah16,1,xyzzyaaaj16)
xyzzyaaad16=1
do xyzzyaaae16=2,xyzzyaabs1
xyzzyaaac16=xyzzyaaad16+1
xyzzyaaad16=xyzzyaaan1(xyzzyaaae16)
xyzzyaaau16(1:6)=0.d0
xyzzyaaap16(1:3)=0.d0
xyzzyaaas16=0.d0
do xyzzyaaaf16=xyzzyaaac16,xyzzyaaad16
xyzzyaaar16=sum(xyzzyaaaz1(1:periodicity,xyzzyaaaf16)*xyzzyaaaq16(1:pe&
&riodicity))
xyzzyaaap16(1:periodicity)=xyzzyaaap16(1:periodicity) -xyzzyaaaz1(1:pe&
&riodicity,xyzzyaaaf16)*sin(xyzzyaaar16)
xyzzyaaav16=cos(xyzzyaaar16)
xyzzyaaas16=xyzzyaaas16+xyzzyaaav16
if(present(orbsderivs))then
xyzzyaaau16(1)=xyzzyaaau16(1)-xyzzyaaaz1(1,xyzzyaaaf16)*xyzzyaaaz1(1,x&
&yzzyaaaf16)*xyzzyaaav16
if(dimensionality>1)then
xyzzyaaau16(2)=xyzzyaaau16(2)-xyzzyaaaz1(2,xyzzyaaaf16)*xyzzyaaaz1(2,x&
&yzzyaaaf16)*xyzzyaaav16
xyzzyaaau16(4)=xyzzyaaau16(4)-xyzzyaaaz1(1,xyzzyaaaf16)*xyzzyaaaz1(2,x&
&yzzyaaaf16)*xyzzyaaav16
if(dimensionality==3)then
xyzzyaaau16(3)=xyzzyaaau16(3)-xyzzyaaaz1(3,xyzzyaaaf16)*xyzzyaaaz1(3,x&
&yzzyaaaf16)*xyzzyaaav16
xyzzyaaau16(5)=xyzzyaaau16(5)-xyzzyaaaz1(1,xyzzyaaaf16)*xyzzyaaaz1(3,x&
&yzzyaaaf16)*xyzzyaaav16
xyzzyaaau16(6)=xyzzyaaau16(6)-xyzzyaaaz1(2,xyzzyaaaf16)*xyzzyaaaz1(3,x&
&yzzyaaaf16)*xyzzyaaav16
endif
endif
endif
enddo
xyzzyaaar16=xyzzyaacx1(xyzzyaaae16,xyzzyaaag16)
orbgrad(1:3,xyzzyaaah16,1)=orbgrad(1:3,xyzzyaaah16,1)+xyzzyaaar16*xyzz&
&yaaap16(1:3)
orblap(xyzzyaaah16,1)=orblap(xyzzyaaah16,1)-xyzzyaaar16*xyzzyaaba1(xyz&
&zyaaac16)*xyzzyaaas16
if(present(orbsderivs))orbsderivs(1:6,xyzzyaaah16,1)=orbsderivs(1:6,xy&
&zzyaaah16,1)+xyzzyaaar16*xyzzyaaau16(1:6)
enddo
enddo
else
do xyzzyaaah16=1,heg_nele(jspin)
orbval(xyzzyaaah16,1)=orbval(xyzzyaaah16,1)+xyzzyaacx1(1,xyzzyaaag16)
xyzzyaaaq16(1:periodicity)=xyzzyaaax16(1:periodicity) -site_pos_in_cel&
&l(1:periodicity,xyzzyaaah16,1,xyzzyaaaj16)
xyzzyaaad16=1
do xyzzyaaae16=2,xyzzyaabs1
xyzzyaaau16(1:6)=0.d0
xyzzyaaap16(1:3)=0.d0
xyzzyaaas16=0.d0
xyzzyaaac16=xyzzyaaad16+1
xyzzyaaad16=xyzzyaaan1(xyzzyaaae16)
do xyzzyaaaf16=xyzzyaaac16,xyzzyaaad16
xyzzyaaar16=sum(xyzzyaaaz1(1:periodicity,xyzzyaaaf16)*xyzzyaaaq16(1:pe&
&riodicity))
xyzzyaaap16(1:periodicity)=xyzzyaaap16(1:periodicity) -xyzzyaaaz1(1:pe&
&riodicity,xyzzyaaaf16)*sin(xyzzyaaar16)
xyzzyaaav16=cos(xyzzyaaar16)
xyzzyaaas16=xyzzyaaas16+xyzzyaaav16
if(present(orbsderivs))then
xyzzyaaau16(1)=xyzzyaaau16(1)-xyzzyaaaz1(1,xyzzyaaaf16)*xyzzyaaaz1(1,x&
&yzzyaaaf16)*xyzzyaaav16
if(dimensionality>1)then
xyzzyaaau16(2)=xyzzyaaau16(2)-xyzzyaaaz1(2,xyzzyaaaf16)*xyzzyaaaz1(2,x&
&yzzyaaaf16)*xyzzyaaav16
xyzzyaaau16(4)=xyzzyaaau16(4)-xyzzyaaaz1(1,xyzzyaaaf16)*xyzzyaaaz1(2,x&
&yzzyaaaf16)*xyzzyaaav16
if(dimensionality==3)then
xyzzyaaau16(3)=xyzzyaaau16(3)-xyzzyaaaz1(3,xyzzyaaaf16)*xyzzyaaaz1(3,x&
&yzzyaaaf16)*xyzzyaaav16
xyzzyaaau16(5)=xyzzyaaau16(5)-xyzzyaaaz1(1,xyzzyaaaf16)*xyzzyaaaz1(3,x&
&yzzyaaaf16)*xyzzyaaav16
xyzzyaaau16(6)=xyzzyaaau16(6)-xyzzyaaaz1(2,xyzzyaaaf16)*xyzzyaaaz1(3,x&
&yzzyaaaf16)*xyzzyaaav16
endif
endif
endif
enddo
xyzzyaaar16=xyzzyaacx1(xyzzyaaae16,xyzzyaaag16)
orbval(xyzzyaaah16,1)=orbval(xyzzyaaah16,1)+xyzzyaaar16*xyzzyaaas16
orbgrad(1:3,xyzzyaaah16,1)=orbgrad(1:3,xyzzyaaah16,1)+xyzzyaaar16*xyzz&
&yaaap16(1:3)
orblap(xyzzyaaah16,1)=orblap(xyzzyaaah16,1)-xyzzyaaar16*xyzzyaaba1(xyz&
&zyaaac16)*xyzzyaaas16
if(present(orbsderivs))orbsderivs(1:6,xyzzyaaah16,1)=orbsderivs(1:6,xy&
&zzyaaah16,1)+xyzzyaaar16*xyzzyaaau16(1:6)
enddo
enddo
endif
enddo
endif
end subroutine free_crystal_orb_eval
subroutine free_pairing_orb_eval(eevecs,jspin,lnorb,norb,orbmask,val,f&
&sd,orbval,orbgrad,orblap,orbsderivs)
implicit none
integer,intent(in) :: jspin,lnorb,norb
real(dp),intent(in) :: eevecs(4,netot)
real(dp),intent(inout) :: orbval(lnorb,*),orbgrad(3,lnorb,*),orblap(ln&
&orb,*)
real(dp),intent(inout),optional :: orbsderivs(6,lnorb,*)
logical,intent(in) :: val,fsd,orbmask(*)
integer xyzzyaaaa17,xyzzyaaab17,xyzzyaaac17,xyzzyaaad17,xyzzyaaae17,xy&
&zzyaaaf17,xyzzyaaag17,xyzzyaaah17,xyzzyaaai17,xyzzyaaaj17
real(dp) xyzzyaaak17,xyzzyaaal17,xyzzyaaam17(3),xyzzyaaan17(3),xyzzyaa&
&ao17,xyzzyaaap17,xyzzyaaaq17,xyzzyaaar17,xyzzyaaas17,xyzzyaaat17,xyzz&
&yaaau17,xyzzyaaav17,xyzzyaaaw17,xyzzyaaax17,xyzzyaaay17,xyzzyaaaz17,x&
&yzzyaaba17,xyzzyaabb17,xyzzyaabc17,xyzzyaabd17,xyzzyaabe17,xyzzyaabf1&
&7,xyzzyaabg17(3),xyzzyaabh17,xyzzyaabi17,xyzzyaabj17,xyzzyaabk17,xyzz&
&yaabl17,xyzzyaabm17,xyzzyaabn17,xyzzyaabo17,xyzzyaabp17,xyzzyaabq17,x&
&yzzyaabr17(3),xyzzyaabs17(3),xyzzyaabt17(6),xyzzyaabu17(6),xyzzyaabv1&
&7,xyzzyaabw17,xyzzyaabx17,xyzzyaaby17,xyzzyaabz17,xyzzyaaca17,xyzzyaa&
&cb17(3),xyzzyaacc17(3),xyzzyaacd17(3),xyzzyaace17(3),xyzzyaacf17(6),x&
&yzzyaacg17(6),xyzzyaach17(6),xyzzyaaci17(6),xyzzyaacj17(6),xyzzyaack1&
&7(6),xyzzyaacl17,xyzzyaacm17,xyzzyaacn17,xyzzyaaco17,xyzzyaacp17,xyzz&
&yaacq17(6)
real(dp),allocatable,save :: xyzzyaacr17(:,:)
logical xyzzyaacs17,xyzzyaact17
logical,save :: xyzzyaacu17=.true.
if(xyzzyaacu17)then
if(xyzzyaabl1.and.any(k_offset/=0.d0))then
allocate(xyzzyaacr17(norb,real1_complex2),stat=xyzzyaaaj17)
call check_alloc(xyzzyaaaj17,'FREE_PAIRING_ORB_EVAL','orbval_save')
xyzzyaacr17=0.d0
endif
xyzzyaacu17=.false.
endif
xyzzyaacs17=val.or.(xyzzyaabl1.and.any(k_offset/=0.d0))
xyzzyaact17=fsd
select case(xyzzyaabu1)
case(0)
xyzzyaaaf17=1
end select
xyzzyaaaw17=xyzzyaaat1*xyzzyaada1(xyzzyaaaf17)
if(xyzzyaacs17)then
if(.not.val)xyzzyaacr17(1:norb,1:real1_complex2)=orbval(1:norb,1:real1&
&_complex2)
orbval(1:norb,1:real1_complex2)=0.d0
endif
if(xyzzyaact17)then
orbgrad(1:3,1:norb,1:real1_complex2)=0.d0
orblap(1:norb,1:real1_complex2)=0.d0
if(present(orbsderivs))orbsderivs(1:6,1:norb,1:real1_complex2)=0.d0
endif
if(xyzzyaadj1)then
xyzzyaaap17=-minval(xyzzyaacv1(:,xyzzyaaaf17))
if(xyzzyaacs17)then
do xyzzyaaah17=1,netot
if(.not.orbmask(xyzzyaaah17))cycle
if(jspin>which_spin(xyzzyaaah17))then
xyzzyaaav17=-1.d0
else
xyzzyaaav17=1.d0
endif
xyzzyaaak17=eevecs(4,xyzzyaaah17)*eevecs(4,xyzzyaaah17)
xyzzyaaao17=xyzzyaaak17*xyzzyaaap17
if(xyzzyaaao17>xyzzyaaav1)then
do xyzzyaaaa17=1,xyzzyaabq1
xyzzyaabn17=xyzzyaaaw17*exp(-xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzy&
&aaak17)
do xyzzyaaai17=1,real1_complex2
orbval(xyzzyaaah17,xyzzyaaai17)=orbval(xyzzyaaah17,xyzzyaaai17)+xyzzya&
&acu1(xyzzyaaai17,xyzzyaaaa17,xyzzyaaaf17)*xyzzyaabn17
enddo
enddo
endif
enddo
endif
if(xyzzyaact17)then
if(dimensionality==3)then
do xyzzyaaah17=1,netot
if(.not.orbmask(xyzzyaaah17))cycle
if(jspin>which_spin(xyzzyaaah17))then
xyzzyaaav17=-1.d0
else
xyzzyaaav17=1.d0
endif
xyzzyaaak17=eevecs(4,xyzzyaaah17)*eevecs(4,xyzzyaaah17)
xyzzyaaao17=xyzzyaaak17*xyzzyaaap17
if(xyzzyaaao17>xyzzyaaav1)then
do xyzzyaaaa17=1,xyzzyaabq1
xyzzyaabn17=xyzzyaaaw17*2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*exp(-&
&xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaak17)
do xyzzyaaai17=1,real1_complex2
xyzzyaaal17=xyzzyaabn17*xyzzyaacu1(xyzzyaaai17,xyzzyaaaa17,xyzzyaaaf17&
&)
orbgrad(1:3,xyzzyaaah17,xyzzyaaai17)=orbgrad(1:3,xyzzyaaah17,xyzzyaaai&
&17)-xyzzyaaav17*eevecs(1:3,xyzzyaaah17)*xyzzyaaal17
orblap(xyzzyaaah17,xyzzyaaai17)=orblap(xyzzyaaah17,xyzzyaaai17)+(2.d0*&
&xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaak17-3.d0)*xyzzyaaal17
if(present(orbsderivs))then
xyzzyaaaq17=xyzzyaaav17*eevecs(1,xyzzyaaah17)
xyzzyaaar17=xyzzyaaav17*eevecs(2,xyzzyaaah17)
xyzzyaaas17=xyzzyaaav17*eevecs(3,xyzzyaaah17)
orbsderivs(1,xyzzyaaah17,xyzzyaaai17)=orbsderivs(1,xyzzyaaah17,xyzzyaa&
&ai17)+xyzzyaaal17*(2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaaq&
&17*xyzzyaaaq17-1.d0)
orbsderivs(2,xyzzyaaah17,xyzzyaaai17)=orbsderivs(2,xyzzyaaah17,xyzzyaa&
&ai17)+xyzzyaaal17*(2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaar&
&17*xyzzyaaar17-1.d0)
orbsderivs(3,xyzzyaaah17,xyzzyaaai17)=orbsderivs(3,xyzzyaaah17,xyzzyaa&
&ai17)+xyzzyaaal17*(2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaas&
&17*xyzzyaaas17-1.d0)
orbsderivs(4,xyzzyaaah17,xyzzyaaai17)=orbsderivs(4,xyzzyaaah17,xyzzyaa&
&ai17)+xyzzyaaal17*2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaaq1&
&7*xyzzyaaar17
orbsderivs(5,xyzzyaaah17,xyzzyaaai17)=orbsderivs(5,xyzzyaaah17,xyzzyaa&
&ai17)+xyzzyaaal17*2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaaq1&
&7*xyzzyaaas17
orbsderivs(6,xyzzyaaah17,xyzzyaaai17)=orbsderivs(6,xyzzyaaah17,xyzzyaa&
&ai17)+xyzzyaaal17*2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaar1&
&7*xyzzyaaas17
endif
enddo
enddo
endif
enddo
elseif(dimensionality==2)then
do xyzzyaaah17=1,netot
if(.not.orbmask(xyzzyaaah17))cycle
if(jspin>which_spin(xyzzyaaah17))then
xyzzyaaav17=-1.d0
else
xyzzyaaav17=1.d0
endif
xyzzyaaak17=eevecs(4,xyzzyaaah17)*eevecs(4,xyzzyaaah17)
xyzzyaaao17=xyzzyaaak17*xyzzyaaap17
if(xyzzyaaao17>xyzzyaaav1)then
do xyzzyaaaa17=1,xyzzyaabq1
xyzzyaabn17=xyzzyaaaw17*2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*exp(-&
&xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaak17)
do xyzzyaaai17=1,real1_complex2
xyzzyaaal17=xyzzyaabn17*xyzzyaacu1(xyzzyaaai17,xyzzyaaaa17,xyzzyaaaf17&
&)
orbgrad(1:2,xyzzyaaah17,xyzzyaaai17)=orbgrad(1:2,xyzzyaaah17,xyzzyaaai&
&17)-xyzzyaaav17*eevecs(1:2,xyzzyaaah17)*xyzzyaaal17
orblap(xyzzyaaah17,xyzzyaaai17)=orblap(xyzzyaaah17,xyzzyaaai17)+(2.d0*&
&xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaak17-2.d0)*xyzzyaaal17
if(present(orbsderivs))then
xyzzyaaaq17=xyzzyaaav17*eevecs(1,xyzzyaaah17)
xyzzyaaar17=xyzzyaaav17*eevecs(2,xyzzyaaah17)
orbsderivs(1,xyzzyaaah17,xyzzyaaai17)=orbsderivs(1,xyzzyaaah17,xyzzyaa&
&ai17)+xyzzyaaal17*(2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaaq&
&17*xyzzyaaaq17-1.d0)
orbsderivs(2,xyzzyaaah17,xyzzyaaai17)=orbsderivs(2,xyzzyaaah17,xyzzyaa&
&ai17)+xyzzyaaal17*(2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaar&
&17*xyzzyaaar17-1.d0)
orbsderivs(4,xyzzyaaah17,xyzzyaaai17)=orbsderivs(4,xyzzyaaah17,xyzzyaa&
&ai17)+xyzzyaaal17*2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaaq1&
&7*xyzzyaaar17
endif
enddo
enddo
endif
enddo
else
do xyzzyaaah17=1,netot
if(.not.orbmask(xyzzyaaah17))cycle
if(jspin>which_spin(xyzzyaaah17))then
xyzzyaaav17=-1.d0
else
xyzzyaaav17=1.d0
endif
xyzzyaaak17=eevecs(4,xyzzyaaah17)*eevecs(4,xyzzyaaah17)
xyzzyaaao17=xyzzyaaak17*xyzzyaaap17
if(xyzzyaaao17>xyzzyaaav1)then
do xyzzyaaaa17=1,xyzzyaabq1
xyzzyaaal17=xyzzyaaaw17*2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*exp(-&
&xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaak17)
do xyzzyaaai17=1,real1_complex2
xyzzyaaal17=xyzzyaabn17*xyzzyaacu1(xyzzyaaai17,xyzzyaaaa17,xyzzyaaaf17&
&)
orbgrad(1,xyzzyaaah17,xyzzyaaai17)=orbgrad(1,xyzzyaaah17,xyzzyaaai17)-&
&xyzzyaaav17*eevecs(1,xyzzyaaah17)*xyzzyaaal17
orblap(xyzzyaaah17,xyzzyaaai17)=orblap(xyzzyaaah17,xyzzyaaai17)+(2.d0*&
&xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaak17-1.d0)*xyzzyaaal17
if(present(orbsderivs))orbsderivs(1,xyzzyaaah17,xyzzyaaai17)=orbsderiv&
&s(1,xyzzyaaah17,xyzzyaaai17)+xyzzyaaal17*(2.d0*xyzzyaacv1(xyzzyaaaa17&
&,xyzzyaaaf17)*eevecs(1,xyzzyaaah17)*eevecs(1,xyzzyaaah17)-1.d0)
enddo
enddo
endif
enddo
endif
endif
do xyzzyaaag17=2,xyzzyaaad1
xyzzyaaam17(1:3)=xyzzyaaax1(1:3,xyzzyaaag17)
if(xyzzyaacs17)then
do xyzzyaaah17=1,netot
if(.not.orbmask(xyzzyaaah17))cycle
if(jspin>which_spin(xyzzyaaah17))then
xyzzyaaav17=-1.d0
else
xyzzyaaav17=1.d0
endif
xyzzyaaan17(1:dimensionality)=xyzzyaaav17*eevecs(1:dimensionality,xyzz&
&yaaah17)+xyzzyaaam17(1:dimensionality)
xyzzyaaak17=sum(xyzzyaaan17(1:dimensionality)**2)
xyzzyaaao17=xyzzyaaak17*xyzzyaaap17
if(xyzzyaaao17>xyzzyaaav1)then
do xyzzyaaaa17=1,xyzzyaabq1
xyzzyaabn17=xyzzyaaaw17*exp(-xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzy&
&aaak17)
do xyzzyaaai17=1,real1_complex2
orbval(xyzzyaaah17,xyzzyaaai17)=orbval(xyzzyaaah17,xyzzyaaai17)+xyzzya&
&acu1(xyzzyaaai17,xyzzyaaaa17,xyzzyaaaf17)*xyzzyaabn17
enddo
enddo
endif
enddo
endif
if(xyzzyaact17)then
if(dimensionality==3)then
do xyzzyaaah17=1,netot
if(.not.orbmask(xyzzyaaah17))cycle
if(jspin>which_spin(xyzzyaaah17))then
xyzzyaaav17=-1.d0
else
xyzzyaaav17=1.d0
endif
xyzzyaaan17(1:3)=xyzzyaaav17*eevecs(1:3,xyzzyaaah17)+xyzzyaaam17(1:3)
xyzzyaaak17=sum(xyzzyaaan17**2)
xyzzyaaao17=xyzzyaaak17*xyzzyaaap17
if(xyzzyaaao17>xyzzyaaav1)then
do xyzzyaaaa17=1,xyzzyaabq1
xyzzyaabn17=xyzzyaaaw17*2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*exp(-&
&xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaak17)
do xyzzyaaai17=1,real1_complex2
xyzzyaaal17=xyzzyaabn17*xyzzyaacu1(xyzzyaaai17,xyzzyaaaa17,xyzzyaaaf17&
&)
orbgrad(1:3,xyzzyaaah17,xyzzyaaai17)=orbgrad(1:3,xyzzyaaah17,xyzzyaaai&
&17)-xyzzyaaan17(1:3)*xyzzyaaal17
orblap(xyzzyaaah17,xyzzyaaai17)=orblap(xyzzyaaah17,xyzzyaaai17)+(2.d0*&
&xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaak17-3.d0)*xyzzyaaal17
if(present(orbsderivs))then
xyzzyaaaq17=xyzzyaaan17(1)
xyzzyaaar17=xyzzyaaan17(2)
xyzzyaaas17=xyzzyaaan17(3)
orbsderivs(1,xyzzyaaah17,xyzzyaaai17)=orbsderivs(1,xyzzyaaah17,xyzzyaa&
&ai17)+xyzzyaaal17*(2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaaq&
&17*xyzzyaaaq17-1.d0)
orbsderivs(2,xyzzyaaah17,xyzzyaaai17)=orbsderivs(2,xyzzyaaah17,xyzzyaa&
&ai17)+xyzzyaaal17*(2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaar&
&17*xyzzyaaar17-1.d0)
orbsderivs(3,xyzzyaaah17,xyzzyaaai17)=orbsderivs(3,xyzzyaaah17,xyzzyaa&
&ai17)+xyzzyaaal17*(2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaas&
&17*xyzzyaaas17-1.d0)
orbsderivs(4,xyzzyaaah17,xyzzyaaai17)=orbsderivs(4,xyzzyaaah17,xyzzyaa&
&ai17)+xyzzyaaal17*2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaaq1&
&7*xyzzyaaar17
orbsderivs(5,xyzzyaaah17,xyzzyaaai17)=orbsderivs(5,xyzzyaaah17,xyzzyaa&
&ai17)+xyzzyaaal17*2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaaq1&
&7*xyzzyaaas17
orbsderivs(6,xyzzyaaah17,xyzzyaaai17)=orbsderivs(6,xyzzyaaah17,xyzzyaa&
&ai17)+xyzzyaaal17*2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaar1&
&7*xyzzyaaas17
endif
enddo
enddo
endif
enddo
elseif(dimensionality==2)then
do xyzzyaaah17=1,netot
if(.not.orbmask(xyzzyaaah17))cycle
if(jspin>which_spin(xyzzyaaah17))then
xyzzyaaav17=-1.d0
else
xyzzyaaav17=1.d0
endif
xyzzyaaan17(1:2)=xyzzyaaav17*eevecs(1:2,xyzzyaaah17)+xyzzyaaam17(1:2)
xyzzyaaak17=xyzzyaaan17(1)*xyzzyaaan17(1)+xyzzyaaan17(2)*xyzzyaaan17(2&
&)
xyzzyaaao17=xyzzyaaak17*xyzzyaaap17
if(xyzzyaaao17>xyzzyaaav1)then
do xyzzyaaaa17=1,xyzzyaabq1
xyzzyaabn17=xyzzyaaaw17*2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*exp(-&
&xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaak17)
do xyzzyaaai17=1,real1_complex2
xyzzyaaal17=xyzzyaabn17*xyzzyaacu1(xyzzyaaai17,xyzzyaaaa17,xyzzyaaaf17&
&)
orbgrad(1:2,xyzzyaaah17,xyzzyaaai17)=orbgrad(1:2,xyzzyaaah17,xyzzyaaai&
&17)-xyzzyaaan17(1:2)*xyzzyaaal17
orblap(xyzzyaaah17,xyzzyaaai17)=orblap(xyzzyaaah17,xyzzyaaai17)+(2.d0*&
&xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaak17-2.d0)*xyzzyaaal17
if(present(orbsderivs))then
xyzzyaaaq17=xyzzyaaan17(1)
xyzzyaaar17=xyzzyaaan17(2)
orbsderivs(1,xyzzyaaah17,xyzzyaaai17)=orbsderivs(1,xyzzyaaah17,xyzzyaa&
&ai17)+xyzzyaaal17*(2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaaq&
&17*xyzzyaaaq17-1.d0)
orbsderivs(2,xyzzyaaah17,xyzzyaaai17)=orbsderivs(2,xyzzyaaah17,xyzzyaa&
&ai17)+xyzzyaaal17*(2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaar&
&17*xyzzyaaar17-1.d0)
orbsderivs(4,xyzzyaaah17,xyzzyaaai17)=orbsderivs(4,xyzzyaaah17,xyzzyaa&
&ai17)+xyzzyaaal17*2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaaq1&
&7*xyzzyaaar17
endif
enddo
enddo
endif
enddo
else
do xyzzyaaah17=1,netot
if(.not.orbmask(xyzzyaaah17))cycle
if(jspin>which_spin(xyzzyaaah17))then
xyzzyaaav17=-1.d0
else
xyzzyaaav17=1.d0
endif
xyzzyaaan17(1)=xyzzyaaav17*eevecs(1,xyzzyaaah17)+xyzzyaaam17(1)
xyzzyaaak17=xyzzyaaan17(1)*xyzzyaaan17(1)
xyzzyaaao17=xyzzyaaak17*xyzzyaaap17
if(xyzzyaaao17>xyzzyaaav1)then
do xyzzyaaaa17=1,xyzzyaabq1
xyzzyaabn17=xyzzyaaaw17*2.d0*xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*exp(-&
&xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaak17)
do xyzzyaaai17=1,real1_complex2
xyzzyaaal17=xyzzyaabn17*xyzzyaacu1(xyzzyaaai17,xyzzyaaaa17,xyzzyaaaf17&
&)
orbgrad(1,xyzzyaaah17,xyzzyaaai17)=orbgrad(1,xyzzyaaah17,xyzzyaaai17)-&
&xyzzyaaan17(1)*xyzzyaaal17
orblap(xyzzyaaah17,xyzzyaaai17)=orblap(xyzzyaaah17,xyzzyaaai17)+(2.d0*&
&xyzzyaacv1(xyzzyaaaa17,xyzzyaaaf17)*xyzzyaaak17-1.d0)*xyzzyaaal17
if(present(orbsderivs))orbsderivs(1,xyzzyaaah17,xyzzyaaai17)=orbsderiv&
&s(1,xyzzyaaah17,xyzzyaaai17)+xyzzyaaal17*(2.d0*xyzzyaacv1(xyzzyaaaa17&
&,xyzzyaaaf17)*xyzzyaaak17-1.d0)
enddo
enddo
endif
enddo
endif
endif
enddo
endif
if(xyzzyaadi1)then
if(.not.xyzzyaact17)then
do xyzzyaaah17=1,netot
if(.not.orbmask(xyzzyaaah17))cycle
if(jspin>which_spin(xyzzyaaah17))then
xyzzyaaav17=-1.d0
else
xyzzyaaav17=1.d0
endif
xyzzyaaan17=eevecs(1:3,xyzzyaaah17)*xyzzyaaav17
orbval(xyzzyaaah17,1)=orbval(xyzzyaaah17,1)+xyzzyaaaw17*xyzzyaacw1(1,1&
&,xyzzyaaaf17)
if(xyzzyaabl1)orbval(xyzzyaaah17,2)=orbval(xyzzyaaah17,2)+xyzzyaaaw17*&
&xyzzyaacw1(2,1,xyzzyaaaf17)
xyzzyaaae17=1
do xyzzyaaab17=2,xyzzyaabr1
xyzzyaabq17=0.d0
xyzzyaaad17=xyzzyaaae17+1
xyzzyaaae17=xyzzyaaan1(xyzzyaaab17)
do xyzzyaaac17=xyzzyaaad17,xyzzyaaae17
xyzzyaaal17=ddot(dimensionality,xyzzyaaaz1(1,xyzzyaaac17),1,xyzzyaaan1&
&7(1),1)
xyzzyaabq17=xyzzyaabq17+cos(xyzzyaaal17)
enddo
xyzzyaabv17=xyzzyaaaw17*xyzzyaabq17
orbval(xyzzyaaah17,1)=orbval(xyzzyaaah17,1)+xyzzyaacw1(1,xyzzyaaab17,x&
&yzzyaaaf17)*xyzzyaabv17
if(xyzzyaabl1)orbval(xyzzyaaah17,2)=orbval(xyzzyaaah17,2)+xyzzyaacw1(2&
&,xyzzyaaab17,xyzzyaaaf17)*xyzzyaabv17
enddo
enddo
else
do xyzzyaaah17=1,netot
if(.not.orbmask(xyzzyaaah17))cycle
if(jspin>which_spin(xyzzyaaah17))then
xyzzyaaav17=-1.d0
else
xyzzyaaav17=1.d0
endif
xyzzyaaan17=eevecs(1:3,xyzzyaaah17)*xyzzyaaav17
if(xyzzyaacs17)then
orbval(xyzzyaaah17,1)=orbval(xyzzyaaah17,1)+xyzzyaaaw17*xyzzyaacw1(1,1&
&,xyzzyaaaf17)
if(xyzzyaabl1)orbval(xyzzyaaah17,2)=orbval(xyzzyaaah17,2)+xyzzyaaaw17*&
&xyzzyaacw1(2,1,xyzzyaaaf17)
endif
xyzzyaaae17=1
do xyzzyaaab17=2,xyzzyaabr1
xyzzyaabr17(:)=0.d0
xyzzyaabq17=0.d0
xyzzyaabt17=0.d0
xyzzyaaad17=xyzzyaaae17+1
xyzzyaaae17=xyzzyaaan1(xyzzyaaab17)
do xyzzyaaac17=xyzzyaaad17,xyzzyaaae17
xyzzyaaal17=ddot(dimensionality,xyzzyaaaz1(1,xyzzyaaac17),1,xyzzyaaan1&
&7(1),1)
xyzzyaaat17=cos(xyzzyaaal17)
xyzzyaaau17=sin(xyzzyaaal17)
xyzzyaabq17=xyzzyaabq17+xyzzyaaat17
xyzzyaabr17(1:dimensionality)=xyzzyaabr17(1:dimensionality)-xyzzyaaaz1&
&(1:dimensionality,xyzzyaaac17)*xyzzyaaau17
if(present(orbsderivs))then
xyzzyaabt17(1)=xyzzyaabt17(1)-xyzzyaaaz1(1,xyzzyaaac17)*xyzzyaaaz1(1,x&
&yzzyaaac17)*xyzzyaaat17
if(dimensionality>1)then
xyzzyaabt17(2)=xyzzyaabt17(2)-xyzzyaaaz1(2,xyzzyaaac17)*xyzzyaaaz1(2,x&
&yzzyaaac17)*xyzzyaaat17
xyzzyaabt17(4)=xyzzyaabt17(4)-xyzzyaaaz1(1,xyzzyaaac17)*xyzzyaaaz1(2,x&
&yzzyaaac17)*xyzzyaaat17
if(dimensionality==3)then
xyzzyaabt17(3)=xyzzyaabt17(3)-xyzzyaaaz1(3,xyzzyaaac17)*xyzzyaaaz1(3,x&
&yzzyaaac17)*xyzzyaaat17
xyzzyaabt17(5)=xyzzyaabt17(5)-xyzzyaaaz1(1,xyzzyaaac17)*xyzzyaaaz1(3,x&
&yzzyaaac17)*xyzzyaaat17
xyzzyaabt17(6)=xyzzyaabt17(6)-xyzzyaaaz1(2,xyzzyaaac17)*xyzzyaaaz1(3,x&
&yzzyaaac17)*xyzzyaaat17
endif
endif
endif
enddo
xyzzyaabo17=xyzzyaaaw17*xyzzyaacw1(1,xyzzyaaab17,xyzzyaaaf17)
if(xyzzyaabl1)xyzzyaabp17=xyzzyaaaw17*xyzzyaacw1(2,xyzzyaaab17,xyzzyaa&
&af17)
if(xyzzyaacs17)then
orbval(xyzzyaaah17,1)=orbval(xyzzyaaah17,1)+xyzzyaabo17*xyzzyaabq17
if(xyzzyaabl1)orbval(xyzzyaaah17,2)=orbval(xyzzyaaah17,2)+xyzzyaabp17*&
&xyzzyaabq17
endif
orbgrad(1:dimensionality,xyzzyaaah17,1)=orbgrad(1:dimensionality,xyzzy&
&aaah17,1)+xyzzyaabo17*xyzzyaabr17(1:dimensionality)
if(xyzzyaabl1)orbgrad(1:dimensionality,xyzzyaaah17,2)=orbgrad(1:dimens&
&ionality,xyzzyaaah17,2)+xyzzyaabp17*xyzzyaabr17(1:dimensionality)
if(.not.present(orbsderivs))then
xyzzyaabv17=xyzzyaaba1(xyzzyaaad17)*xyzzyaabq17
orblap(xyzzyaaah17,1)=orblap(xyzzyaaah17,1)-xyzzyaabo17*xyzzyaabv17
if(xyzzyaabl1)orblap(xyzzyaaah17,2)=orblap(xyzzyaaah17,2)-xyzzyaabp17*&
&xyzzyaabv17
else
orbsderivs(1:6,xyzzyaaah17,1)=orbsderivs(1:6,xyzzyaaah17,1)+xyzzyaabo1&
&7*xyzzyaabt17
if(xyzzyaabl1)orbsderivs(1:6,xyzzyaaah17,2)=orbsderivs(1:6,xyzzyaaah17&
&,2) +xyzzyaabp17*xyzzyaabt17
endif
enddo
enddo
endif
endif
if(xyzzyaadk1)then
if(.not.xyzzyaact17)then
do xyzzyaaah17=1,netot
if(.not.orbmask(xyzzyaaah17))cycle
if(jspin>which_spin(xyzzyaaah17))then
xyzzyaaav17=-1.d0
else
xyzzyaaav17=1.d0
endif
xyzzyaaax17=eevecs(4,xyzzyaaah17)
if(xyzzyaaax17>=xyzzyaacp1)cycle
call xyzzyaafp1(xyzzyaaax17,xyzzyaaaf17,xyzzyaaay17,xyzzyaabb17)
orbval(xyzzyaaah17,1)=orbval(xyzzyaaah17,1)+xyzzyaaaw17*xyzzyaaay17
if(xyzzyaabl1)orbval(xyzzyaaah17,2)=orbval(xyzzyaaah17,2)+xyzzyaaaw17*&
&xyzzyaabb17
enddo
elseif(.not.xyzzyaacs17)then
do xyzzyaaah17=1,netot
if(.not.orbmask(xyzzyaaah17))cycle
if(jspin>which_spin(xyzzyaaah17))then
xyzzyaaav17=-1.d0
else
xyzzyaaav17=1.d0
endif
xyzzyaaax17=eevecs(4,xyzzyaaah17)
xyzzyaaan17(1:dimensionality)=xyzzyaaav17*eevecs(1:dimensionality,xyzz&
&yaaah17)
if(xyzzyaaax17>=xyzzyaacp1)cycle
call xyzzyaafq1(xyzzyaaax17,xyzzyaaaf17,xyzzyaaay17,xyzzyaabb17,xyzzya&
&aaz17,xyzzyaabc17,xyzzyaaba17,xyzzyaabd17)
call xyzzyaaft1(xyzzyaaax17,xyzzyaaan17,xyzzyaabg17,xyzzyaabe17,xyzzya&
&abf17)
xyzzyaaaz17=xyzzyaaaz17*xyzzyaaaw17
xyzzyaaba17=xyzzyaaba17*xyzzyaaaw17
orbgrad(1:dimensionality,xyzzyaaah17,1)=orbgrad(1:dimensionality,xyzzy&
&aaah17,1)+xyzzyaabg17(1:dimensionality)*xyzzyaaaz17
if(.not.present(orbsderivs))then
orblap(xyzzyaaah17,1)=orblap(xyzzyaaah17,1)+xyzzyaabe17*xyzzyaaba17+xy&
&zzyaabf17*xyzzyaaaz17
else
xyzzyaabh17=xyzzyaaaz17/xyzzyaaax17
xyzzyaabi17=xyzzyaaba17-xyzzyaabh17
orbsderivs(1,xyzzyaaah17,1)=orbsderivs(1,xyzzyaaah17,1)+xyzzyaabg17(1)&
&*xyzzyaabg17(1)*xyzzyaabi17+xyzzyaabh17
if(dimensionality>1)then
orbsderivs(2,xyzzyaaah17,1)=orbsderivs(2,xyzzyaaah17,1)+xyzzyaabg17(2)&
&*xyzzyaabg17(2)*xyzzyaabi17+xyzzyaabh17
orbsderivs(4,xyzzyaaah17,1)=orbsderivs(4,xyzzyaaah17,1)+xyzzyaabg17(1)&
&*xyzzyaabg17(2)*xyzzyaabi17
if(dimensionality>2)then
orbsderivs(3,xyzzyaaah17,1)=orbsderivs(3,xyzzyaaah17,1)+xyzzyaabg17(3)&
&*xyzzyaabg17(3)*xyzzyaabi17+xyzzyaabh17
orbsderivs(5,xyzzyaaah17,1)=orbsderivs(5,xyzzyaaah17,1)+xyzzyaabg17(1)&
&*xyzzyaabg17(3)*xyzzyaabi17
orbsderivs(6,xyzzyaaah17,1)=orbsderivs(6,xyzzyaaah17,1)+xyzzyaabg17(2)&
&*xyzzyaabg17(3)*xyzzyaabi17
endif
endif
endif
if(xyzzyaabl1)then
xyzzyaabc17=xyzzyaabc17*xyzzyaaaw17
xyzzyaabd17=xyzzyaabd17*xyzzyaaaw17
orbgrad(1:dimensionality,xyzzyaaah17,2)=orbgrad(1:dimensionality,xyzzy&
&aaah17,2)+xyzzyaabg17(1:dimensionality)*xyzzyaabc17
if(.not.present(orbsderivs))then
orblap(xyzzyaaah17,2)=orblap(xyzzyaaah17,2)+xyzzyaabe17*xyzzyaabd17+xy&
&zzyaabf17*xyzzyaabc17
else
xyzzyaabh17=xyzzyaabc17/xyzzyaaax17
xyzzyaabi17=xyzzyaabd17-xyzzyaabh17
orbsderivs(1,xyzzyaaah17,2)=orbsderivs(1,xyzzyaaah17,2)+xyzzyaabg17(1)&
&*xyzzyaabg17(1)*xyzzyaabi17+xyzzyaabh17
if(dimensionality>1)then
orbsderivs(2,xyzzyaaah17,2)=orbsderivs(2,xyzzyaaah17,2)+xyzzyaabg17(2)&
&*xyzzyaabg17(2)*xyzzyaabi17+xyzzyaabh17
orbsderivs(4,xyzzyaaah17,2)=orbsderivs(4,xyzzyaaah17,2)+xyzzyaabg17(1)&
&*xyzzyaabg17(2)*xyzzyaabi17
if(dimensionality>2)then
orbsderivs(3,xyzzyaaah17,2)=orbsderivs(3,xyzzyaaah17,2)+xyzzyaabg17(3)&
&*xyzzyaabg17(3)*xyzzyaabi17+xyzzyaabh17
orbsderivs(5,xyzzyaaah17,2)=orbsderivs(5,xyzzyaaah17,2)+xyzzyaabg17(1)&
&*xyzzyaabg17(3)*xyzzyaabi17
orbsderivs(6,xyzzyaaah17,2)=orbsderivs(6,xyzzyaaah17,2)+xyzzyaabg17(2)&
&*xyzzyaabg17(3)*xyzzyaabi17
endif
endif
endif
endif
enddo
else
do xyzzyaaah17=1,netot
if(.not.orbmask(xyzzyaaah17))cycle
if(jspin>which_spin(xyzzyaaah17))then
xyzzyaaav17=-1.d0
else
xyzzyaaav17=1.d0
endif
xyzzyaaax17=eevecs(4,xyzzyaaah17)
xyzzyaaan17(1:dimensionality)=xyzzyaaav17*eevecs(1:dimensionality,xyzz&
&yaaah17)
if(xyzzyaaax17>=xyzzyaacp1)cycle
call xyzzyaafq1(xyzzyaaax17,xyzzyaaaf17,xyzzyaaay17,xyzzyaabb17,xyzzya&
&aaz17,xyzzyaabc17,xyzzyaaba17,xyzzyaabd17)
call xyzzyaaft1(xyzzyaaax17,xyzzyaaan17,xyzzyaabg17,xyzzyaabe17,xyzzya&
&abf17)
xyzzyaaaz17=xyzzyaaaz17*xyzzyaaaw17
xyzzyaaba17=xyzzyaaba17*xyzzyaaaw17
orbval(xyzzyaaah17,1)=orbval(xyzzyaaah17,1)+xyzzyaaay17*xyzzyaaaw17
orbgrad(1:dimensionality,xyzzyaaah17,1)=orbgrad(1:dimensionality,xyzzy&
&aaah17,1)+xyzzyaabg17(1:dimensionality)*xyzzyaaaz17
if(.not.present(orbsderivs))then
orblap(xyzzyaaah17,1)=orblap(xyzzyaaah17,1)+xyzzyaabe17*xyzzyaaba17+xy&
&zzyaabf17*xyzzyaaaz17
else
xyzzyaabh17=xyzzyaaaz17/xyzzyaaax17
xyzzyaabi17=xyzzyaaba17-xyzzyaabh17
orbsderivs(1,xyzzyaaah17,1)=orbsderivs(1,xyzzyaaah17,1)+xyzzyaabg17(1)&
&*xyzzyaabg17(1)*xyzzyaabi17+xyzzyaabh17
if(dimensionality>1)then
orbsderivs(2,xyzzyaaah17,1)=orbsderivs(2,xyzzyaaah17,1)+xyzzyaabg17(2)&
&*xyzzyaabg17(2)*xyzzyaabi17+xyzzyaabh17
orbsderivs(4,xyzzyaaah17,1)=orbsderivs(4,xyzzyaaah17,1)+xyzzyaabg17(1)&
&*xyzzyaabg17(2)*xyzzyaabi17
if(dimensionality>2)then
orbsderivs(3,xyzzyaaah17,1)=orbsderivs(3,xyzzyaaah17,1)+xyzzyaabg17(3)&
&*xyzzyaabg17(3)*xyzzyaabi17+xyzzyaabh17
orbsderivs(5,xyzzyaaah17,1)=orbsderivs(5,xyzzyaaah17,1)+xyzzyaabg17(1)&
&*xyzzyaabg17(3)*xyzzyaabi17
orbsderivs(6,xyzzyaaah17,1)=orbsderivs(6,xyzzyaaah17,1)+xyzzyaabg17(2)&
&*xyzzyaabg17(3)*xyzzyaabi17
endif
endif
endif
if(xyzzyaabl1)then
xyzzyaabc17=xyzzyaabc17*xyzzyaaaw17
xyzzyaabd17=xyzzyaabd17*xyzzyaaaw17
orbval(xyzzyaaah17,2)=orbval(xyzzyaaah17,2)+xyzzyaabb17*xyzzyaaaw17
orbgrad(1:dimensionality,xyzzyaaah17,2)=orbgrad(1:dimensionality,xyzzy&
&aaah17,2)+xyzzyaabg17(1:dimensionality)*xyzzyaabc17
if(.not.present(orbsderivs))then
orblap(xyzzyaaah17,2)=orblap(xyzzyaaah17,2)+xyzzyaabe17*xyzzyaabd17+xy&
&zzyaabf17*xyzzyaabc17
else
xyzzyaabh17=xyzzyaabc17/xyzzyaaax17
xyzzyaabi17=xyzzyaabd17-xyzzyaabh17
orbsderivs(1,xyzzyaaah17,2)=orbsderivs(1,xyzzyaaah17,2)+xyzzyaabg17(1)&
&*xyzzyaabg17(1)*xyzzyaabi17+xyzzyaabh17
if(dimensionality>1)then
orbsderivs(2,xyzzyaaah17,2)=orbsderivs(2,xyzzyaaah17,2)+xyzzyaabg17(2)&
&*xyzzyaabg17(2)*xyzzyaabi17+xyzzyaabh17
orbsderivs(4,xyzzyaaah17,2)=orbsderivs(4,xyzzyaaah17,2)+xyzzyaabg17(1)&
&*xyzzyaabg17(2)*xyzzyaabi17
if(dimensionality>2)then
orbsderivs(3,xyzzyaaah17,2)=orbsderivs(3,xyzzyaaah17,2)+xyzzyaabg17(3)&
&*xyzzyaabg17(3)*xyzzyaabi17+xyzzyaabh17
orbsderivs(5,xyzzyaaah17,2)=orbsderivs(5,xyzzyaaah17,2)+xyzzyaabg17(1)&
&*xyzzyaabg17(3)*xyzzyaabi17
orbsderivs(6,xyzzyaaah17,2)=orbsderivs(6,xyzzyaaah17,2)+xyzzyaabg17(2)&
&*xyzzyaabg17(3)*xyzzyaabi17
endif
endif
endif
endif
enddo
endif
do xyzzyaaag17=2,xyzzyaaad1
xyzzyaaam17(1:dimensionality)=xyzzyaaax1(1:dimensionality,xyzzyaaag17)
if(.not.xyzzyaact17)then
do xyzzyaaah17=1,netot
if(.not.orbmask(xyzzyaaah17))cycle
if(jspin>which_spin(xyzzyaaah17))then
xyzzyaaav17=-1.d0
else
xyzzyaaav17=1.d0
endif
xyzzyaaan17(1:dimensionality)=xyzzyaaav17*eevecs(1:dimensionality,xyzz&
&yaaah17)+xyzzyaaam17(1:dimensionality)
xyzzyaaax17=sqrt(sum(xyzzyaaan17(1:dimensionality)**2))
if(xyzzyaaax17>=xyzzyaacp1)cycle
call xyzzyaafp1(xyzzyaaax17,xyzzyaaaf17,xyzzyaaay17,xyzzyaabb17)
orbval(xyzzyaaah17,1)=orbval(xyzzyaaah17,1)+xyzzyaaaw17*xyzzyaaay17
if(xyzzyaabl1)orbval(xyzzyaaah17,2)=orbval(xyzzyaaah17,2)+xyzzyaaaw17*&
&xyzzyaabb17
enddo
elseif(.not.xyzzyaacs17)then
do xyzzyaaah17=1,netot
if(.not.orbmask(xyzzyaaah17))cycle
if(jspin>which_spin(xyzzyaaah17))then
xyzzyaaav17=-1.d0
else
xyzzyaaav17=1.d0
endif
xyzzyaaan17(1:dimensionality)=xyzzyaaav17*eevecs(1:dimensionality,xyzz&
&yaaah17)+xyzzyaaam17(1:dimensionality)
xyzzyaaax17=sqrt(sum(xyzzyaaan17(1:dimensionality)**2))
if(xyzzyaaax17>=xyzzyaacp1)cycle
call xyzzyaafq1(xyzzyaaax17,xyzzyaaaf17,xyzzyaaay17,xyzzyaabb17,xyzzya&
&aaz17,xyzzyaabc17,xyzzyaaba17,xyzzyaabd17)
call xyzzyaaft1(xyzzyaaax17,xyzzyaaan17,xyzzyaabg17,xyzzyaabe17,xyzzya&
&abf17)
xyzzyaaaz17=xyzzyaaaz17*xyzzyaaaw17
xyzzyaaba17=xyzzyaaba17*xyzzyaaaw17
orbgrad(1:dimensionality,xyzzyaaah17,1)=orbgrad(1:dimensionality,xyzzy&
&aaah17,1)+xyzzyaabg17(1:dimensionality)*xyzzyaaaz17
if(.not.present(orbsderivs))then
orblap(xyzzyaaah17,1)=orblap(xyzzyaaah17,1)+xyzzyaabe17*xyzzyaaba17+xy&
&zzyaabf17*xyzzyaaaz17
else
xyzzyaabh17=xyzzyaaaz17/xyzzyaaax17
xyzzyaabi17=xyzzyaaba17-xyzzyaabh17
orbsderivs(1,xyzzyaaah17,1)=orbsderivs(1,xyzzyaaah17,1)+xyzzyaabg17(1)&
&*xyzzyaabg17(1)*xyzzyaabi17+xyzzyaabh17
if(dimensionality>1)then
orbsderivs(2,xyzzyaaah17,1)=orbsderivs(2,xyzzyaaah17,1)+xyzzyaabg17(2)&
&*xyzzyaabg17(2)*xyzzyaabi17+xyzzyaabh17
orbsderivs(4,xyzzyaaah17,1)=orbsderivs(4,xyzzyaaah17,1)+xyzzyaabg17(1)&
&*xyzzyaabg17(2)*xyzzyaabi17
if(dimensionality>2)then
orbsderivs(3,xyzzyaaah17,1)=orbsderivs(3,xyzzyaaah17,1)+xyzzyaabg17(3)&
&*xyzzyaabg17(3)*xyzzyaabi17+xyzzyaabh17
orbsderivs(5,xyzzyaaah17,1)=orbsderivs(5,xyzzyaaah17,1)+xyzzyaabg17(1)&
&*xyzzyaabg17(3)*xyzzyaabi17
orbsderivs(6,xyzzyaaah17,1)=orbsderivs(6,xyzzyaaah17,1)+xyzzyaabg17(2)&
&*xyzzyaabg17(3)*xyzzyaabi17
endif
endif
endif
if(xyzzyaabl1)then
xyzzyaabc17=xyzzyaabc17*xyzzyaaaw17
xyzzyaabd17=xyzzyaabd17*xyzzyaaaw17
orbgrad(1:dimensionality,xyzzyaaah17,2)=orbgrad(1:dimensionality,xyzzy&
&aaah17,2)+xyzzyaabg17(1:dimensionality)*xyzzyaabc17
if(.not.present(orbsderivs))then
orblap(xyzzyaaah17,2)=orblap(xyzzyaaah17,2)+xyzzyaabe17*xyzzyaabd17+xy&
&zzyaabf17*xyzzyaabc17
else
xyzzyaabh17=xyzzyaabc17/xyzzyaaax17
xyzzyaabi17=xyzzyaabd17-xyzzyaabh17
orbsderivs(1,xyzzyaaah17,2)=orbsderivs(1,xyzzyaaah17,2)+xyzzyaabg17(1)&
&*xyzzyaabg17(1)*xyzzyaabi17+xyzzyaabh17
if(dimensionality>1)then
orbsderivs(2,xyzzyaaah17,2)=orbsderivs(2,xyzzyaaah17,2)+xyzzyaabg17(2)&
&*xyzzyaabg17(2)*xyzzyaabi17+xyzzyaabh17
orbsderivs(4,xyzzyaaah17,2)=orbsderivs(4,xyzzyaaah17,2)+xyzzyaabg17(1)&
&*xyzzyaabg17(2)*xyzzyaabi17
if(dimensionality>2)then
orbsderivs(3,xyzzyaaah17,2)=orbsderivs(3,xyzzyaaah17,2)+xyzzyaabg17(3)&
&*xyzzyaabg17(3)*xyzzyaabi17+xyzzyaabh17
orbsderivs(5,xyzzyaaah17,2)=orbsderivs(5,xyzzyaaah17,2)+xyzzyaabg17(1)&
&*xyzzyaabg17(3)*xyzzyaabi17
orbsderivs(6,xyzzyaaah17,2)=orbsderivs(6,xyzzyaaah17,2)+xyzzyaabg17(2)&
&*xyzzyaabg17(3)*xyzzyaabi17
endif
endif
endif
endif
enddo
else
do xyzzyaaah17=1,netot
if(.not.orbmask(xyzzyaaah17))cycle
if(jspin>which_spin(xyzzyaaah17))then
xyzzyaaav17=-1.d0
else
xyzzyaaav17=1.d0
endif
xyzzyaaan17(1:dimensionality)=xyzzyaaav17*eevecs(1:dimensionality,xyzz&
&yaaah17)+xyzzyaaam17(1:dimensionality)
xyzzyaaax17=sqrt(sum(xyzzyaaan17(1:dimensionality)**2))
if(xyzzyaaax17>=xyzzyaacp1)cycle
call xyzzyaafq1(xyzzyaaax17,xyzzyaaaf17,xyzzyaaay17,xyzzyaabb17,xyzzya&
&aaz17,xyzzyaabc17,xyzzyaaba17,xyzzyaabd17)
call xyzzyaaft1(xyzzyaaax17,xyzzyaaan17,xyzzyaabg17,xyzzyaabe17,xyzzya&
&abf17)
xyzzyaaaz17=xyzzyaaaz17*xyzzyaaaw17
xyzzyaaba17=xyzzyaaba17*xyzzyaaaw17
orbval(xyzzyaaah17,1)=orbval(xyzzyaaah17,1)+xyzzyaaay17*xyzzyaaaw17
orbgrad(1:dimensionality,xyzzyaaah17,1)=orbgrad(1:dimensionality,xyzzy&
&aaah17,1)+xyzzyaabg17(1:dimensionality)*xyzzyaaaz17
if(.not.present(orbsderivs))then
orblap(xyzzyaaah17,1)=orblap(xyzzyaaah17,1)+xyzzyaabe17*xyzzyaaba17+xy&
&zzyaabf17*xyzzyaaaz17
else
xyzzyaabh17=xyzzyaaaz17/xyzzyaaax17
xyzzyaabi17=xyzzyaaba17-xyzzyaabh17
orbsderivs(1,xyzzyaaah17,1)=orbsderivs(1,xyzzyaaah17,1)+xyzzyaabg17(1)&
&*xyzzyaabg17(1)*xyzzyaabi17+xyzzyaabh17
if(dimensionality>1)then
orbsderivs(2,xyzzyaaah17,1)=orbsderivs(2,xyzzyaaah17,1)+xyzzyaabg17(2)&
&*xyzzyaabg17(2)*xyzzyaabi17+xyzzyaabh17
orbsderivs(4,xyzzyaaah17,1)=orbsderivs(4,xyzzyaaah17,1)+xyzzyaabg17(1)&
&*xyzzyaabg17(2)*xyzzyaabi17
if(dimensionality>2)then
orbsderivs(3,xyzzyaaah17,1)=orbsderivs(3,xyzzyaaah17,1)+xyzzyaabg17(3)&
&*xyzzyaabg17(3)*xyzzyaabi17+xyzzyaabh17
orbsderivs(5,xyzzyaaah17,1)=orbsderivs(5,xyzzyaaah17,1)+xyzzyaabg17(1)&
&*xyzzyaabg17(3)*xyzzyaabi17
orbsderivs(6,xyzzyaaah17,1)=orbsderivs(6,xyzzyaaah17,1)+xyzzyaabg17(2)&
&*xyzzyaabg17(3)*xyzzyaabi17
endif
endif
endif
if(xyzzyaabl1)then
xyzzyaabc17=xyzzyaabc17*xyzzyaaaw17
xyzzyaabd17=xyzzyaabd17*xyzzyaaaw17
orbval(xyzzyaaah17,2)=orbval(xyzzyaaah17,2)+xyzzyaabb17*xyzzyaaaw17
orbgrad(1:dimensionality,xyzzyaaah17,2)=orbgrad(1:dimensionality,xyzzy&
&aaah17,2)+xyzzyaabg17(1:dimensionality)*xyzzyaabc17
if(.not.present(orbsderivs))then
orblap(xyzzyaaah17,2)=orblap(xyzzyaaah17,2)+xyzzyaabe17*xyzzyaabd17+xy&
&zzyaabf17*xyzzyaabc17
else
xyzzyaabh17=xyzzyaabc17/xyzzyaaax17
xyzzyaabi17=xyzzyaabd17-xyzzyaabh17
orbsderivs(1,xyzzyaaah17,2)=orbsderivs(1,xyzzyaaah17,2)+xyzzyaabg17(1)&
&*xyzzyaabg17(1)*xyzzyaabi17+xyzzyaabh17
if(dimensionality>1)then
orbsderivs(2,xyzzyaaah17,2)=orbsderivs(2,xyzzyaaah17,2)+xyzzyaabg17(2)&
&*xyzzyaabg17(2)*xyzzyaabi17+xyzzyaabh17
orbsderivs(4,xyzzyaaah17,2)=orbsderivs(4,xyzzyaaah17,2)+xyzzyaabg17(1)&
&*xyzzyaabg17(2)*xyzzyaabi17
if(dimensionality>2)then
orbsderivs(3,xyzzyaaah17,2)=orbsderivs(3,xyzzyaaah17,2)+xyzzyaabg17(3)&
&*xyzzyaabg17(3)*xyzzyaabi17+xyzzyaabh17
orbsderivs(5,xyzzyaaah17,2)=orbsderivs(5,xyzzyaaah17,2)+xyzzyaabg17(1)&
&*xyzzyaabg17(3)*xyzzyaabi17
orbsderivs(6,xyzzyaaah17,2)=orbsderivs(6,xyzzyaaah17,2)+xyzzyaabg17(2)&
&*xyzzyaabg17(3)*xyzzyaabi17
endif
endif
endif
endif
enddo
endif
enddo
endif
if(xyzzyaadl1)then
if(.not.xyzzyaact17)then
do xyzzyaaah17=1,netot
if(.not.orbmask(xyzzyaaah17))cycle
if(jspin>which_spin(xyzzyaaah17))then
xyzzyaaav17=-1.d0
else
xyzzyaaav17=1.d0
endif
xyzzyaaax17=eevecs(4,xyzzyaaah17)
xyzzyaaan17(1:dimensionality)=xyzzyaaav17*eevecs(1:dimensionality,xyzz&
&yaaah17)
call xyzzyaafr1(xyzzyaaax17,xyzzyaaan17,xyzzyaaaf17,xyzzyaabj17,xyzzya&
&abk17)
orbval(xyzzyaaah17,1)=orbval(xyzzyaaah17,1)+xyzzyaaaw17*xyzzyaabj17
if(xyzzyaabl1)orbval(xyzzyaaah17,2)=orbval(xyzzyaaah17,2)+xyzzyaaaw17*&
&xyzzyaabk17
enddo
else
do xyzzyaaah17=1,netot
if(.not.orbmask(xyzzyaaah17))cycle
if(jspin>which_spin(xyzzyaaah17))then
xyzzyaaav17=-1.d0
else
xyzzyaaav17=1.d0
endif
xyzzyaaax17=eevecs(4,xyzzyaaah17)
xyzzyaaan17(1:dimensionality)=xyzzyaaav17*eevecs(1:dimensionality,xyzz&
&yaaah17)
if(present(orbsderivs))then
call xyzzyaafs1(xyzzyaaax17,xyzzyaaan17,xyzzyaaaf17,xyzzyaabj17,xyzzya&
&abk17,xyzzyaabr17,xyzzyaabs17,xyzzyaabl17,xyzzyaabm17,xyzzyaabt17,xyz&
&zyaabu17)
orbsderivs(1:6,xyzzyaaah17,1)=orbsderivs(1:6,xyzzyaaah17,1)+xyzzyaabt1&
&7(:)*xyzzyaaaw17
if(xyzzyaabl1)orbsderivs(1:6,xyzzyaaah17,2)=orbsderivs(1:6,xyzzyaaah17&
&,2)+xyzzyaabu17(:)*xyzzyaaaw17
else
call xyzzyaafs1(xyzzyaaax17,xyzzyaaan17,xyzzyaaaf17,xyzzyaabj17,xyzzya&
&abk17,xyzzyaabr17,xyzzyaabs17,xyzzyaabl17,xyzzyaabm17)
orblap(xyzzyaaah17,1)=orblap(xyzzyaaah17,1)+xyzzyaabl17*xyzzyaaaw17
if(xyzzyaabl1)orblap(xyzzyaaah17,2)=orblap(xyzzyaaah17,2)+xyzzyaabm17*&
&xyzzyaaaw17
endif
if(xyzzyaacs17)then
orbval(xyzzyaaah17,1)=orbval(xyzzyaaah17,1)+xyzzyaabj17*xyzzyaaaw17
if(xyzzyaabl1)orbval(xyzzyaaah17,2)=orbval(xyzzyaaah17,2)+xyzzyaabk17*&
&xyzzyaaaw17
endif
orbgrad(1:dimensionality,xyzzyaaah17,1)=orbgrad(1:dimensionality,xyzzy&
&aaah17,1)+xyzzyaabr17(1:dimensionality)*xyzzyaaaw17
if(xyzzyaabl1)orbgrad(1:dimensionality,xyzzyaaah17,2)=orbgrad(1:dimens&
&ionality,xyzzyaaah17,2)+xyzzyaabs17(1:dimensionality)*xyzzyaaaw17
enddo
endif
endif
if(xyzzyaabl1.and.any(k_offset/=0.d0))then
if(present(orbsderivs))then
xyzzyaacq17(1:6)=(/k_offset(1)*k_offset(1),k_offset(2)*k_offset(2),k_o&
&ffset(3)*k_offset(3),k_offset(1)*k_offset(2),k_offset(1)*k_offset(3),&
&k_offset(2)*k_offset(3)/)
else
xyzzyaacp17=ddot(3,k_offset(1),1,k_offset(1),1)
endif
do xyzzyaaah17=1,netot
if(.not.orbmask(xyzzyaaah17))cycle
if(jspin>which_spin(xyzzyaaah17))then
xyzzyaaav17=-1.d0
else
xyzzyaaav17=1.d0
endif
xyzzyaaan17(1:3)=xyzzyaaav17*eevecs(1:3,xyzzyaaah17)
xyzzyaabw17=ddot(3,xyzzyaaan17(1),1,k_offset(1),1)
xyzzyaabx17=cos(xyzzyaabw17)
xyzzyaaby17=sin(xyzzyaabw17)
xyzzyaabz17=orbval(xyzzyaaah17,1)
xyzzyaaca17=orbval(xyzzyaaah17,2)
orbval(xyzzyaaah17,1)=xyzzyaabx17*xyzzyaabz17-xyzzyaaby17*xyzzyaaca17
orbval(xyzzyaaah17,2)=xyzzyaabx17*xyzzyaaca17+xyzzyaaby17*xyzzyaabz17
if(xyzzyaact17)then
xyzzyaacb17(1:3)=-xyzzyaaby17*k_offset(1:3)
xyzzyaacc17(1:3)=xyzzyaabx17*k_offset(1:3)
xyzzyaacd17(1:3)=orbgrad(1:3,xyzzyaaah17,1)
xyzzyaace17(1:3)=orbgrad(1:3,xyzzyaaah17,2)
orbgrad(1:3,xyzzyaaah17,1)=xyzzyaacb17(1:3)*xyzzyaabz17-xyzzyaacc17(1:&
&3)*xyzzyaaca17+xyzzyaabx17*xyzzyaacd17(1:3)-xyzzyaaby17*xyzzyaace17(1&
&:3)
orbgrad(1:3,xyzzyaaah17,2)=xyzzyaacb17(1:3)*xyzzyaaca17+xyzzyaacc17(1:&
&3)*xyzzyaabz17+xyzzyaabx17*xyzzyaace17(1:3)+xyzzyaaby17*xyzzyaacd17(1&
&:3)
if(present(orbsderivs))then
xyzzyaacf17(1:6)=-xyzzyaabx17*xyzzyaacq17(1:6)
xyzzyaacg17(1:6)=-xyzzyaaby17*xyzzyaacq17(1:6)
xyzzyaach17(1:6)=orbsderivs(1:6,xyzzyaaah17,1)
xyzzyaaci17(1:6)=orbsderivs(1:6,xyzzyaaah17,2)
xyzzyaacj17(1)=2*(xyzzyaacb17(1)*xyzzyaacd17(1)-xyzzyaacc17(1)*xyzzyaa&
&ce17(1))
xyzzyaacj17(2)=2*(xyzzyaacb17(2)*xyzzyaacd17(2)-xyzzyaacc17(2)*xyzzyaa&
&ce17(2))
xyzzyaacj17(3)=2*(xyzzyaacb17(3)*xyzzyaacd17(3)-xyzzyaacc17(3)*xyzzyaa&
&ce17(3))
xyzzyaacj17(4)=xyzzyaacb17(1)*xyzzyaacd17(2)+xyzzyaacb17(2)*xyzzyaacd1&
&7(1)-xyzzyaacc17(1)*xyzzyaace17(2)-xyzzyaacc17(2)*xyzzyaace17(1)
xyzzyaacj17(5)=xyzzyaacb17(1)*xyzzyaacd17(3)+xyzzyaacb17(3)*xyzzyaacd1&
&7(1)-xyzzyaacc17(1)*xyzzyaace17(3)-xyzzyaacc17(3)*xyzzyaace17(1)
xyzzyaacj17(6)=xyzzyaacb17(2)*xyzzyaacd17(3)+xyzzyaacb17(3)*xyzzyaacd1&
&7(2)-xyzzyaacc17(2)*xyzzyaace17(3)-xyzzyaacc17(3)*xyzzyaace17(2)
xyzzyaack17(1)=2*(xyzzyaacb17(1)*xyzzyaace17(1)+xyzzyaacc17(1)*xyzzyaa&
&cd17(1))
xyzzyaack17(2)=2*(xyzzyaacb17(2)*xyzzyaace17(2)+xyzzyaacc17(2)*xyzzyaa&
&cd17(2))
xyzzyaack17(3)=2*(xyzzyaacb17(3)*xyzzyaace17(3)+xyzzyaacc17(3)*xyzzyaa&
&cd17(3))
xyzzyaack17(4)=xyzzyaacb17(1)*xyzzyaace17(2)+xyzzyaacb17(2)*xyzzyaace1&
&7(1)+xyzzyaacc17(1)*xyzzyaacd17(2)+xyzzyaacc17(2)*xyzzyaacd17(1)
xyzzyaack17(5)=xyzzyaacb17(1)*xyzzyaace17(3)+xyzzyaacb17(3)*xyzzyaace1&
&7(1)+xyzzyaacc17(1)*xyzzyaacd17(3)+xyzzyaacc17(3)*xyzzyaacd17(1)
xyzzyaack17(6)=xyzzyaacb17(2)*xyzzyaace17(3)+xyzzyaacb17(3)*xyzzyaace1&
&7(2)+xyzzyaacc17(2)*xyzzyaacd17(3)+xyzzyaacc17(3)*xyzzyaacd17(2)
orbsderivs(1:6,xyzzyaaah17,1)=xyzzyaacf17(1:6)*xyzzyaabz17-xyzzyaacg17&
&(1:6)*xyzzyaaca17+xyzzyaacj17(1:6)+xyzzyaabx17*xyzzyaach17(1:6)-xyzzy&
&aaby17*xyzzyaaci17(1:6)
orbsderivs(1:6,xyzzyaaah17,2)=xyzzyaacf17(1:6)*xyzzyaaca17+xyzzyaacg17&
&(1:6)*xyzzyaabz17+xyzzyaack17(1:6)+xyzzyaabx17*xyzzyaaci17(1:6)+xyzzy&
&aaby17*xyzzyaach17(1:6)
else
xyzzyaacl17=-xyzzyaabx17*xyzzyaacp17
xyzzyaacm17=-xyzzyaaby17*xyzzyaacp17
xyzzyaacn17=orblap(xyzzyaaah17,1)
xyzzyaaco17=orblap(xyzzyaaah17,2)
orblap(xyzzyaaah17,1)=xyzzyaacl17*xyzzyaabz17-xyzzyaacm17*xyzzyaaca17+&
&2.d0*(ddot(3,xyzzyaacb17(1),1,xyzzyaacd17(1),1)-ddot(3,xyzzyaacc17(1)&
&,1,xyzzyaace17(1),1))+xyzzyaabx17*xyzzyaacn17-xyzzyaaby17*xyzzyaaco17
orblap(xyzzyaaah17,2)=xyzzyaacl17*xyzzyaaca17+xyzzyaacm17*xyzzyaabz17+&
&2.d0*(ddot(3,xyzzyaacb17(1),1,xyzzyaace17(1),1)+ddot(3,xyzzyaacc17(1)&
&,1,xyzzyaacd17(1),1))+xyzzyaabx17*xyzzyaaco17+xyzzyaaby17*xyzzyaacn17
endif
endif
enddo
endif
if(xyzzyaacs17.and..not.val)orbval(1:norb,1:real1_complex2)=xyzzyaacr1&
&7(1:norb,1:real1_complex2)
end subroutine free_pairing_orb_eval
subroutine free_biex1_orb_eval(rvec,jspin,lnorb,norb,orbmask,val,fsd,o&
&rbval,orbgrad,orblap,orbsderivs)
implicit none
integer,intent(in) :: jspin,lnorb,norb
real(dp),intent(in) :: rvec(3)
real(dp),intent(inout) :: orbval(lnorb,*),orbgrad(3,lnorb,*),orblap(ln&
&orb,*)
real(dp),intent(inout),optional :: orbsderivs(6,lnorb,*)
logical,intent(in) :: val,fsd,orbmask(*)
if(val)orbval(1:norb,1:real1_complex2)=1.d0
if(fsd)then
orbgrad(1:3,1:norb,1:real1_complex2)=0.d0
orblap(1:norb,1:real1_complex2)=0.d0
endif
end subroutine free_biex1_orb_eval
subroutine xyzzyaafp1(r,s,polyr,polyi)
implicit none
integer,intent(in) :: s
real(dp),intent(in) :: r
real(dp),intent(out):: polyr,polyi
integer xyzzyaaaa19
real(dp) xyzzyaaab19,xyzzyaaac19,xyzzyaaad19
if(.not.xyzzyaabl1)then
xyzzyaaab19=r
polyr=xyzzyaadb1(1,0,s)
polyi=0.d0
do xyzzyaaaa19=1,xyzzyaabx1
polyr=polyr+xyzzyaadb1(1,xyzzyaaaa19,s)*xyzzyaaab19
xyzzyaaab19=xyzzyaaab19*r
enddo
else
xyzzyaaab19=r
polyr=xyzzyaadb1(1,0,s)
polyi=xyzzyaadb1(2,0,s)
do xyzzyaaaa19=1,xyzzyaabx1
polyr=polyr+xyzzyaadb1(1,xyzzyaaaa19,s)*xyzzyaaab19
polyi=polyi+xyzzyaadb1(2,xyzzyaaaa19,s)*xyzzyaaab19
xyzzyaaab19=xyzzyaaab19*r
enddo
endif
if(xyzzyaaby1==0)then
xyzzyaaad19=1.d0
else
xyzzyaaac19=1.d0-r*xyzzyaacq1
select case(xyzzyaaby1)
case(1)
xyzzyaaad19=xyzzyaaac19
case(2)
xyzzyaaad19=xyzzyaaac19*xyzzyaaac19
case default
xyzzyaaad19=xyzzyaaac19**xyzzyaaby1
end select
endif
polyr=polyr*xyzzyaaad19
polyi=polyi*xyzzyaaad19
end subroutine xyzzyaafp1
subroutine xyzzyaafq1(r,s,polyr,polyi,dpolyr,dpolyi,d2polyr,d2polyi)
implicit none
integer,intent(in) :: s
real(dp),intent(in) :: r
real(dp),intent(out):: polyr,polyi,dpolyr,dpolyi,d2polyr,d2polyi
integer xyzzyaaaa20
real(dp) xyzzyaaab20,xyzzyaaac20,xyzzyaaad20,xyzzyaaae20,xyzzyaaaf20,x&
&yzzyaaag20,xyzzyaaah20,xyzzyaaai20,xyzzyaaaj20
if(.not.xyzzyaabl1)then
polyr=xyzzyaadb1(1,0,s)+xyzzyaadb1(1,1,s)*r
dpolyr=xyzzyaadb1(1,1,s)
d2polyr=0.d0
polyi=0.d0
dpolyi=0.d0
d2polyi=0.d0
xyzzyaaad20=0.d0
xyzzyaaac20=1.d0
xyzzyaaab20=r
do xyzzyaaaa20=2,xyzzyaabx1
xyzzyaaad20=xyzzyaaac20*xyzzyaaaa20
xyzzyaaac20=xyzzyaaab20*xyzzyaaaa20
xyzzyaaab20=xyzzyaaab20*r
polyr=polyr+xyzzyaadb1(1,xyzzyaaaa20,s)*xyzzyaaab20
dpolyr=dpolyr+xyzzyaadb1(1,xyzzyaaaa20,s)*xyzzyaaac20
d2polyr=d2polyr+xyzzyaadb1(1,xyzzyaaaa20,s)*xyzzyaaad20
enddo
else
polyr=xyzzyaadb1(1,0,s)+xyzzyaadb1(1,1,s)*r
polyi=xyzzyaadb1(2,0,s)+xyzzyaadb1(2,1,s)*r
dpolyr=xyzzyaadb1(1,1,s)
dpolyi=xyzzyaadb1(2,1,s)
d2polyr=0.d0
d2polyi=0.d0
xyzzyaaad20=0.d0
xyzzyaaac20=1.d0
xyzzyaaab20=r
do xyzzyaaaa20=2,xyzzyaabx1
xyzzyaaad20=xyzzyaaac20*xyzzyaaaa20
xyzzyaaac20=xyzzyaaab20*xyzzyaaaa20
xyzzyaaab20=xyzzyaaab20*r
polyr=polyr+xyzzyaadb1(1,xyzzyaaaa20,s)*xyzzyaaab20
polyi=polyi+xyzzyaadb1(2,xyzzyaaaa20,s)*xyzzyaaab20
dpolyr=dpolyr+xyzzyaadb1(1,xyzzyaaaa20,s)*xyzzyaaac20
dpolyi=dpolyi+xyzzyaadb1(2,xyzzyaaaa20,s)*xyzzyaaac20
d2polyr=d2polyr+xyzzyaadb1(1,xyzzyaaaa20,s)*xyzzyaaad20
d2polyi=d2polyi+xyzzyaadb1(2,xyzzyaaaa20,s)*xyzzyaaad20
enddo
endif
if(xyzzyaaby1==0)then
xyzzyaaah20=1.d0
xyzzyaaai20=0.d0
xyzzyaaaj20=0.d0
else
xyzzyaaae20=1.d0-r*xyzzyaacq1
select case(xyzzyaaby1)
case(1)
xyzzyaaah20=xyzzyaaae20
xyzzyaaai20=-xyzzyaacq1
xyzzyaaaj20=0.d0
case(2)
xyzzyaaah20=xyzzyaaae20*xyzzyaaae20
xyzzyaaai20=-2*xyzzyaacq1*xyzzyaaae20
xyzzyaaaj20=2*xyzzyaacr1
case(3)
xyzzyaaaf20=xyzzyaaae20*xyzzyaaae20
xyzzyaaah20=xyzzyaaaf20*xyzzyaaae20
xyzzyaaai20=-3*xyzzyaacq1*xyzzyaaaf20
xyzzyaaaj20=6*xyzzyaacr1*xyzzyaaae20
case default
xyzzyaaaf20=xyzzyaaae20*xyzzyaaae20
xyzzyaaag20=xyzzyaaae20**(xyzzyaaby1-2)
xyzzyaaah20=xyzzyaaaf20*xyzzyaaag20
xyzzyaaai20=-real(xyzzyaaby1,dp)*xyzzyaacq1*xyzzyaaag20*xyzzyaaae20
xyzzyaaaj20=real(xyzzyaaby1*(xyzzyaaby1-1),dp)*xyzzyaacr1*xyzzyaaag20
end select
endif
d2polyr=polyr*xyzzyaaaj20+2*dpolyr*xyzzyaaai20+d2polyr*xyzzyaaah20
d2polyi=polyi*xyzzyaaaj20+2*dpolyi*xyzzyaaai20+d2polyi*xyzzyaaah20
dpolyr=polyr*xyzzyaaai20+dpolyr*xyzzyaaah20
dpolyi=polyi*xyzzyaaai20+dpolyi*xyzzyaaah20
polyr=polyr*xyzzyaaah20
polyi=polyi*xyzzyaaah20
end subroutine xyzzyaafq1
subroutine xyzzyaafr1(r,rvec,s,fsr,fsi)
implicit none
integer,intent(in) :: s
real(dp),intent(in) :: r,rvec(3)
real(dp),intent(out) :: fsr,fsi
integer xyzzyaaaa21
real(dp) xyzzyaaab21,xyzzyaaac21,xyzzyaaad21,xyzzyaaae21,xyzzyaaaf21,x&
&yzzyaaag21,xyzzyaaah21
fsr=0.d0
fsi=0.d0
if(.not.xyzzyaabl1)then
xyzzyaaab21=r*r
do xyzzyaaaa21=1,xyzzyaaca1
xyzzyaaaf21=xyzzyaadd1(xyzzyaaaa21,s)*xyzzyaaab21
xyzzyaaag21=xyzzyaade1(xyzzyaaaa21,s)*r
xyzzyaaah21=1.d0/(1.d0+xyzzyaaag21)
xyzzyaaac21=exp(-xyzzyaaaf21*xyzzyaaah21)
fsr=fsr+xyzzyaadc1(1,xyzzyaaaa21,s)*xyzzyaaac21
enddo
do xyzzyaaaa21=1,xyzzyaacb1
xyzzyaaad21=ddot(3,rvec(1),1,xyzzyaadf1(1,1,xyzzyaaaa21,s),real1_compl&
&ex2)
xyzzyaaaf21=xyzzyaadg1(xyzzyaaaa21,s)*xyzzyaaab21
xyzzyaaag21=xyzzyaadh1(xyzzyaaaa21,s)*r
xyzzyaaah21=1.d0/(1.d0+xyzzyaaag21)
xyzzyaaac21=exp(-xyzzyaaaf21*xyzzyaaah21)
fsr=fsr+xyzzyaaad21*xyzzyaaac21
enddo
else
xyzzyaaab21=r*r
do xyzzyaaaa21=1,xyzzyaaca1
xyzzyaaaf21=xyzzyaadd1(xyzzyaaaa21,s)*xyzzyaaab21
xyzzyaaag21=xyzzyaade1(xyzzyaaaa21,s)*r
xyzzyaaah21=1.d0/(1.d0+xyzzyaaag21)
xyzzyaaac21=exp(-xyzzyaaaf21*xyzzyaaah21)
fsr=fsr+xyzzyaadc1(1,xyzzyaaaa21,s)*xyzzyaaac21
fsi=fsi+xyzzyaadc1(2,xyzzyaaaa21,s)*xyzzyaaac21
enddo
do xyzzyaaaa21=1,xyzzyaacb1
xyzzyaaad21=ddot(3,rvec(1),1,xyzzyaadf1(1,1,xyzzyaaaa21,s),real1_compl&
&ex2)
xyzzyaaae21=ddot(3,rvec(1),1,xyzzyaadf1(2,1,xyzzyaaaa21,s),real1_compl&
&ex2)
xyzzyaaaf21=xyzzyaadg1(xyzzyaaaa21,s)*xyzzyaaab21
xyzzyaaag21=xyzzyaadh1(xyzzyaaaa21,s)*r
xyzzyaaah21=1.d0/(1.d0+xyzzyaaag21)
xyzzyaaac21=exp(-xyzzyaaaf21*xyzzyaaah21)
fsr=fsr+xyzzyaaad21*xyzzyaaac21
fsi=fsi+xyzzyaaae21*xyzzyaaac21
enddo
endif
end subroutine xyzzyaafr1
subroutine xyzzyaafs1(r,rvec,s,fsr,fsi,grad_fsr,grad_fsi,lap_fsr,lap_f&
&si,sderivs_fsr,sderivs_fsi)
implicit none
integer,intent(in) :: s
real(dp),intent(in) :: r,rvec(3)
real(dp),intent(out) :: fsr,fsi,grad_fsr(3),grad_fsi(3),lap_fsr,lap_fs&
&i
real(dp),intent(out),optional :: sderivs_fsr(6),sderivs_fsi(6)
integer xyzzyaaaa22
real(dp) xyzzyaaab22,xyzzyaaac22,xyzzyaaad22(3),xyzzyaaae22,xyzzyaaaf2&
&2,xyzzyaaag22,xyzzyaaah22,xyzzyaaai22,xyzzyaaaj22,xyzzyaaak22,xyzzyaa&
&al22,xyzzyaaam22,xyzzyaaan22,xyzzyaaao22(3),xyzzyaaap22(3),xyzzyaaaq2&
&2,xyzzyaaar22,xyzzyaaas22,xyzzyaaat22,xyzzyaaau22
fsr=0.d0
fsi=0.d0
grad_fsr=0.d0
grad_fsi=0.d0
lap_fsr=0.d0
lap_fsi=0.d0
if(present(sderivs_fsr))then
sderivs_fsr=0.d0
sderivs_fsi=0.d0
endif
if(.not.xyzzyaabl1)then
xyzzyaaab22=r*r
xyzzyaaac22=0.d0
if(r/=0.d0)xyzzyaaac22=1.d0/r
xyzzyaaad22=rvec*xyzzyaaac22
do xyzzyaaaa22=1,xyzzyaaca1
xyzzyaaan22=xyzzyaadd1(xyzzyaaaa22,s)
xyzzyaaaj22=xyzzyaade1(xyzzyaaaa22,s)*r
xyzzyaaak22=1.d0/(1.d0+xyzzyaaaj22)
xyzzyaaal22=xyzzyaaak22*xyzzyaaak22
xyzzyaaaq22=-xyzzyaaan22*xyzzyaaab22*xyzzyaaak22
xyzzyaaas22=-xyzzyaaan22*(2.d0+xyzzyaaaj22)*xyzzyaaal22
xyzzyaaar22=xyzzyaaas22*r
xyzzyaaat22=-2.d0*xyzzyaaan22*xyzzyaaal22*xyzzyaaak22
xyzzyaaau22=xyzzyaaat22+xyzzyaaar22*xyzzyaaar22
xyzzyaaae22=exp(xyzzyaaaq22)
xyzzyaaaf22=xyzzyaadc1(1,xyzzyaaaa22,s)*xyzzyaaae22
fsr=fsr+xyzzyaaaf22
grad_fsr=grad_fsr+xyzzyaaaf22*xyzzyaaar22*xyzzyaaad22
if(present(sderivs_fsr))then
sderivs_fsr(1)=sderivs_fsr(1)+xyzzyaaaf22*(xyzzyaaau22*xyzzyaaad22(1)*&
&xyzzyaaad22(1)+xyzzyaaas22*(1.d0-xyzzyaaad22(1)*xyzzyaaad22(1)))
sderivs_fsr(2)=sderivs_fsr(2)+xyzzyaaaf22*(xyzzyaaau22*xyzzyaaad22(2)*&
&xyzzyaaad22(2)+xyzzyaaas22*(1.d0-xyzzyaaad22(2)*xyzzyaaad22(2)))
sderivs_fsr(3)=sderivs_fsr(3)+xyzzyaaaf22*(xyzzyaaau22*xyzzyaaad22(3)*&
&xyzzyaaad22(3)+xyzzyaaas22*(1.d0-xyzzyaaad22(3)*xyzzyaaad22(3)))
sderivs_fsr(4)=sderivs_fsr(4)+xyzzyaaaf22*(xyzzyaaau22-xyzzyaaas22)*xy&
&zzyaaad22(1)*xyzzyaaad22(2)
sderivs_fsr(5)=sderivs_fsr(5)+xyzzyaaaf22*(xyzzyaaau22-xyzzyaaas22)*xy&
&zzyaaad22(1)*xyzzyaaad22(3)
sderivs_fsr(6)=sderivs_fsr(6)+xyzzyaaaf22*(xyzzyaaau22-xyzzyaaas22)*xy&
&zzyaaad22(2)*xyzzyaaad22(3)
else
lap_fsr=lap_fsr+xyzzyaaaf22*(xyzzyaaau22+xyzzyaaas22*real(dimensionali&
&ty-1,dp))
endif
enddo
do xyzzyaaaa22=1,xyzzyaacb1
xyzzyaaan22=xyzzyaadg1(xyzzyaaaa22,s)
xyzzyaaaj22=xyzzyaadh1(xyzzyaaaa22,s)*r
xyzzyaaao22=xyzzyaadf1(1,:,xyzzyaaaa22,s)
xyzzyaaah22=ddot(3,rvec(1),1,xyzzyaaao22(1),1)
xyzzyaaak22=1.d0/(1.d0+xyzzyaaaj22)
xyzzyaaal22=xyzzyaaak22*xyzzyaaak22
xyzzyaaaq22=-xyzzyaaan22*xyzzyaaab22*xyzzyaaak22
xyzzyaaas22=-xyzzyaaan22*(2.d0+xyzzyaaaj22)*xyzzyaaal22
xyzzyaaar22=xyzzyaaas22*r
xyzzyaaat22=-2.d0*xyzzyaaan22*xyzzyaaal22*xyzzyaaak22
xyzzyaaau22=xyzzyaaat22+xyzzyaaar22*xyzzyaaar22
xyzzyaaae22=exp(xyzzyaaaq22)
xyzzyaaaf22=xyzzyaaah22*xyzzyaaae22
fsr=fsr+xyzzyaaaf22
grad_fsr=grad_fsr+(xyzzyaaao22+xyzzyaaah22*xyzzyaaar22*xyzzyaaad22)*xy&
&zzyaaae22
if(present(sderivs_fsr))then
sderivs_fsr(1)=sderivs_fsr(1)+xyzzyaaae22*(xyzzyaaas22*(xyzzyaaah22*(1&
&.d0-xyzzyaaad22(1)*xyzzyaaad22(1))+2.d0*xyzzyaaao22(1)*rvec(1))+xyzzy&
&aaah22*xyzzyaaad22(1)*xyzzyaaad22(1)*xyzzyaaau22)
sderivs_fsr(2)=sderivs_fsr(2)+xyzzyaaae22*(xyzzyaaas22*(xyzzyaaah22*(1&
&.d0-xyzzyaaad22(2)*xyzzyaaad22(2))+2.d0*xyzzyaaao22(2)*rvec(2))+xyzzy&
&aaah22*xyzzyaaad22(2)*xyzzyaaad22(2)*xyzzyaaau22)
sderivs_fsr(3)=sderivs_fsr(3)+xyzzyaaae22*(xyzzyaaas22*(xyzzyaaah22*(1&
&.d0-xyzzyaaad22(3)*xyzzyaaad22(3))+2.d0*xyzzyaaao22(3)*rvec(3))+xyzzy&
&aaah22*xyzzyaaad22(3)*xyzzyaaad22(3)*xyzzyaaau22)
sderivs_fsr(4)=sderivs_fsr(4)+xyzzyaaae22*(xyzzyaaas22*(-xyzzyaaah22*x&
&yzzyaaad22(1)*xyzzyaaad22(2)+xyzzyaaao22(1)*rvec(2)+xyzzyaaao22(2)*rv&
&ec(1))+xyzzyaaah22*xyzzyaaad22(1)*xyzzyaaad22(2)*xyzzyaaau22)
sderivs_fsr(5)=sderivs_fsr(5)+xyzzyaaae22*(xyzzyaaas22*(-xyzzyaaah22*x&
&yzzyaaad22(1)*xyzzyaaad22(3)+xyzzyaaao22(1)*rvec(3)+xyzzyaaao22(3)*rv&
&ec(1))+xyzzyaaah22*xyzzyaaad22(1)*xyzzyaaad22(3)*xyzzyaaau22)
sderivs_fsr(6)=sderivs_fsr(6)+xyzzyaaae22*(xyzzyaaas22*(-xyzzyaaah22*x&
&yzzyaaad22(2)*xyzzyaaad22(3)+xyzzyaaao22(2)*rvec(3)+xyzzyaaao22(3)*rv&
&ec(2))+xyzzyaaah22*xyzzyaaad22(2)*xyzzyaaad22(3)*xyzzyaaau22)
else
lap_fsr=lap_fsr+xyzzyaaae22*xyzzyaaah22*(xyzzyaaas22*real(1+dimensiona&
&lity,dp)+xyzzyaaau22)
endif
enddo
else
xyzzyaaab22=r*r
xyzzyaaac22=0.d0
if(r/=0.d0)xyzzyaaac22=1.d0/r
xyzzyaaad22=rvec*xyzzyaaac22
do xyzzyaaaa22=1,xyzzyaaca1
xyzzyaaan22=xyzzyaadd1(xyzzyaaaa22,s)
xyzzyaaaj22=xyzzyaade1(xyzzyaaaa22,s)*r
xyzzyaaak22=1.d0/(1.d0+xyzzyaaaj22)
xyzzyaaal22=xyzzyaaak22*xyzzyaaak22
xyzzyaaaq22=-xyzzyaaan22*xyzzyaaab22*xyzzyaaak22
xyzzyaaas22=-xyzzyaaan22*(2.d0+xyzzyaaaj22)*xyzzyaaal22
xyzzyaaar22=xyzzyaaas22*r
xyzzyaaat22=-2.d0*xyzzyaaan22*xyzzyaaal22*xyzzyaaak22
xyzzyaaau22=xyzzyaaat22+xyzzyaaar22*xyzzyaaar22
xyzzyaaae22=exp(xyzzyaaaq22)
xyzzyaaaf22=xyzzyaadc1(1,xyzzyaaaa22,s)*xyzzyaaae22
xyzzyaaag22=xyzzyaadc1(2,xyzzyaaaa22,s)*xyzzyaaae22
fsr=fsr+xyzzyaaaf22
fsi=fsi+xyzzyaaag22
grad_fsr=grad_fsr+xyzzyaaaf22*xyzzyaaar22*xyzzyaaad22
grad_fsi=grad_fsi+xyzzyaaag22*xyzzyaaar22*xyzzyaaad22
if(present(sderivs_fsr))then
xyzzyaaam22=xyzzyaaau22*xyzzyaaad22(1)*xyzzyaaad22(1)+xyzzyaaas22*(1.d&
&0-xyzzyaaad22(1)*xyzzyaaad22(1))
sderivs_fsr(1)=sderivs_fsr(1)+xyzzyaaaf22*xyzzyaaam22
sderivs_fsi(1)=sderivs_fsi(1)+xyzzyaaag22*xyzzyaaam22
xyzzyaaam22=xyzzyaaau22*xyzzyaaad22(2)*xyzzyaaad22(2)+xyzzyaaas22*(1.d&
&0-xyzzyaaad22(2)*xyzzyaaad22(2))
sderivs_fsr(2)=sderivs_fsr(2)+xyzzyaaaf22*xyzzyaaam22
sderivs_fsi(2)=sderivs_fsi(2)+xyzzyaaag22*xyzzyaaam22
xyzzyaaam22=xyzzyaaau22*xyzzyaaad22(3)*xyzzyaaad22(3)+xyzzyaaas22*(1.d&
&0-xyzzyaaad22(3)*xyzzyaaad22(3))
sderivs_fsr(3)=sderivs_fsr(3)+xyzzyaaaf22*xyzzyaaam22
sderivs_fsi(3)=sderivs_fsi(3)+xyzzyaaag22*xyzzyaaam22
xyzzyaaam22=(xyzzyaaau22-xyzzyaaas22)*xyzzyaaad22(1)*xyzzyaaad22(2)
sderivs_fsr(4)=sderivs_fsr(4)+xyzzyaaaf22*xyzzyaaam22
sderivs_fsi(4)=sderivs_fsi(4)+xyzzyaaag22*xyzzyaaam22
xyzzyaaam22=(xyzzyaaau22-xyzzyaaas22)*xyzzyaaad22(1)*xyzzyaaad22(3)
sderivs_fsr(5)=sderivs_fsr(5)+xyzzyaaaf22*xyzzyaaam22
sderivs_fsi(5)=sderivs_fsi(5)+xyzzyaaag22*xyzzyaaam22
xyzzyaaam22=(xyzzyaaau22-xyzzyaaas22)*xyzzyaaad22(2)*xyzzyaaad22(3)
sderivs_fsr(6)=sderivs_fsr(6)+xyzzyaaaf22*xyzzyaaam22
sderivs_fsi(6)=sderivs_fsi(6)+xyzzyaaag22*xyzzyaaam22
else
xyzzyaaam22=xyzzyaaau22+xyzzyaaas22*real(dimensionality-1,dp)
lap_fsr=lap_fsr+xyzzyaaaf22*xyzzyaaam22
lap_fsi=lap_fsi+xyzzyaaag22*xyzzyaaam22
endif
enddo
do xyzzyaaaa22=1,xyzzyaacb1
xyzzyaaan22=xyzzyaadg1(xyzzyaaaa22,s)
xyzzyaaaj22=xyzzyaadh1(xyzzyaaaa22,s)*r
xyzzyaaao22=xyzzyaadf1(1,:,xyzzyaaaa22,s)
xyzzyaaah22=ddot(3,rvec(1),1,xyzzyaaao22(1),1)
xyzzyaaap22=xyzzyaadf1(2,:,xyzzyaaaa22,s)
xyzzyaaai22=ddot(3,rvec(1),1,xyzzyaaap22(1),1)
xyzzyaaak22=1.d0/(1.d0+xyzzyaaaj22)
xyzzyaaal22=xyzzyaaak22*xyzzyaaak22
xyzzyaaaq22=-xyzzyaaan22*xyzzyaaab22*xyzzyaaak22
xyzzyaaas22=-xyzzyaaan22*(2.d0+xyzzyaaaj22)*xyzzyaaal22
xyzzyaaar22=xyzzyaaas22*r
xyzzyaaat22=-2.d0*xyzzyaaan22*xyzzyaaal22*xyzzyaaak22
xyzzyaaau22=xyzzyaaat22+xyzzyaaar22*xyzzyaaar22
xyzzyaaae22=exp(xyzzyaaaq22)
xyzzyaaaf22=xyzzyaaah22*xyzzyaaae22
xyzzyaaag22=xyzzyaaai22*xyzzyaaae22
fsr=fsr+xyzzyaaaf22
fsi=fsi+xyzzyaaag22
grad_fsr=grad_fsr+(xyzzyaaao22+xyzzyaaah22*xyzzyaaar22*xyzzyaaad22)*xy&
&zzyaaae22
grad_fsi=grad_fsi+(xyzzyaaap22+xyzzyaaai22*xyzzyaaar22*xyzzyaaad22)*xy&
&zzyaaae22
if(present(sderivs_fsr))then
sderivs_fsr(1)=sderivs_fsr(1)+xyzzyaaae22*(xyzzyaaas22*(xyzzyaaah22*(1&
&.d0-xyzzyaaad22(1)*xyzzyaaad22(1))+2.d0*xyzzyaaao22(1)*rvec(1))+xyzzy&
&aaah22*xyzzyaaad22(1)*xyzzyaaad22(1)*xyzzyaaau22)
sderivs_fsr(2)=sderivs_fsr(2)+xyzzyaaae22*(xyzzyaaas22*(xyzzyaaah22*(1&
&.d0-xyzzyaaad22(2)*xyzzyaaad22(2))+2.d0*xyzzyaaao22(2)*rvec(2))+xyzzy&
&aaah22*xyzzyaaad22(2)*xyzzyaaad22(2)*xyzzyaaau22)
sderivs_fsr(3)=sderivs_fsr(3)+xyzzyaaae22*(xyzzyaaas22*(xyzzyaaah22*(1&
&.d0-xyzzyaaad22(3)*xyzzyaaad22(3))+2.d0*xyzzyaaao22(3)*rvec(3))+xyzzy&
&aaah22*xyzzyaaad22(3)*xyzzyaaad22(3)*xyzzyaaau22)
sderivs_fsr(4)=sderivs_fsr(4)+xyzzyaaae22*(xyzzyaaas22*(-xyzzyaaah22*x&
&yzzyaaad22(1)*xyzzyaaad22(2)+xyzzyaaao22(1)*rvec(2)+xyzzyaaao22(2)*rv&
&ec(1))+xyzzyaaah22*xyzzyaaad22(1)*xyzzyaaad22(2)*xyzzyaaau22)
sderivs_fsr(5)=sderivs_fsr(5)+xyzzyaaae22*(xyzzyaaas22*(-xyzzyaaah22*x&
&yzzyaaad22(1)*xyzzyaaad22(3)+xyzzyaaao22(1)*rvec(3)+xyzzyaaao22(3)*rv&
&ec(1))+xyzzyaaah22*xyzzyaaad22(1)*xyzzyaaad22(3)*xyzzyaaau22)
sderivs_fsr(6)=sderivs_fsr(6)+xyzzyaaae22*(xyzzyaaas22*(-xyzzyaaah22*x&
&yzzyaaad22(2)*xyzzyaaad22(3)+xyzzyaaao22(2)*rvec(3)+xyzzyaaao22(3)*rv&
&ec(2))+xyzzyaaah22*xyzzyaaad22(2)*xyzzyaaad22(3)*xyzzyaaau22)
sderivs_fsi(1)=sderivs_fsi(1)+xyzzyaaae22*(xyzzyaaas22*(xyzzyaaai22*(1&
&.d0-xyzzyaaad22(1)*xyzzyaaad22(1))+2.d0*xyzzyaaap22(1)*rvec(1))+xyzzy&
&aaai22*xyzzyaaad22(1)*xyzzyaaad22(1)*xyzzyaaau22)
sderivs_fsi(2)=sderivs_fsi(2)+xyzzyaaae22*(xyzzyaaas22*(xyzzyaaai22*(1&
&.d0-xyzzyaaad22(2)*xyzzyaaad22(2))+2.d0*xyzzyaaap22(2)*rvec(2))+xyzzy&
&aaai22*xyzzyaaad22(2)*xyzzyaaad22(2)*xyzzyaaau22)
sderivs_fsi(3)=sderivs_fsi(3)+xyzzyaaae22*(xyzzyaaas22*(xyzzyaaai22*(1&
&.d0-xyzzyaaad22(3)*xyzzyaaad22(3))+2.d0*xyzzyaaap22(3)*rvec(3))+xyzzy&
&aaai22*xyzzyaaad22(3)*xyzzyaaad22(3)*xyzzyaaau22)
sderivs_fsi(4)=sderivs_fsi(4)+xyzzyaaae22*(xyzzyaaas22*(-xyzzyaaai22*x&
&yzzyaaad22(1)*xyzzyaaad22(2)+xyzzyaaap22(1)*rvec(2)+xyzzyaaap22(2)*rv&
&ec(1))+xyzzyaaai22*xyzzyaaad22(1)*xyzzyaaad22(2)*xyzzyaaau22)
sderivs_fsi(5)=sderivs_fsi(5)+xyzzyaaae22*(xyzzyaaas22*(-xyzzyaaai22*x&
&yzzyaaad22(1)*xyzzyaaad22(3)+xyzzyaaap22(1)*rvec(3)+xyzzyaaap22(3)*rv&
&ec(1))+xyzzyaaai22*xyzzyaaad22(1)*xyzzyaaad22(3)*xyzzyaaau22)
sderivs_fsi(6)=sderivs_fsi(6)+xyzzyaaae22*(xyzzyaaas22*(-xyzzyaaai22*x&
&yzzyaaad22(2)*xyzzyaaad22(3)+xyzzyaaap22(2)*rvec(3)+xyzzyaaap22(3)*rv&
&ec(2))+xyzzyaaai22*xyzzyaaad22(2)*xyzzyaaad22(3)*xyzzyaaau22)
else
xyzzyaaam22=xyzzyaaas22*real(1+dimensionality,dp)+xyzzyaaau22
lap_fsr=lap_fsr+xyzzyaaaf22*xyzzyaaam22
lap_fsi=lap_fsi+xyzzyaaag22*xyzzyaaam22
endif
enddo
endif
end subroutine xyzzyaafs1
subroutine xyzzyaaft1(r,vecr,grad,grad2,lap)
implicit none
real(dp),intent(in) :: r,vecr(:)
real(dp),intent(out) :: grad(3),grad2,lap
real(dp) xyzzyaaaa23
grad2=1.d0
if(r==0.d0)then
grad(1:3)=0.d0
lap=0.d0
return
endif
select case(dimensionality)
case(3)
xyzzyaaaa23=1.d0/r
lap=2*xyzzyaaaa23
grad(1:3)=vecr(1:3)*xyzzyaaaa23
case(2)
xyzzyaaaa23=1.d0/r
lap=xyzzyaaaa23
grad(1:2)=vecr(1:2)*xyzzyaaaa23
grad(3)=0.d0
case(1)
grad(1)=sign(1.d0,vecr(1))
grad(2:3)=0.d0
lap=0.d0
end select
end subroutine xyzzyaaft1
subroutine xyzzyaafu1
use slaarnaan, only : lattice_generator
implicit none
integer xyzzyaaaa24,xyzzyaaab24,xyzzyaaac24,xyzzyaaad24,xyzzyaaae24,xy&
&zzyaaaf24,xyzzyaaag24,xyzzyaaah24,xyzzyaaai24,xyzzyaaaj24(3),xyzzyaaa&
&k24,xyzzyaaal24,xyzzyaaam24,xyzzyaaan24
integer,allocatable :: xyzzyaaao24(:),xyzzyaaap24(:),xyzzyaaaq24(:),xy&
&zzyaaar24(:)
real(dp) xyzzyaaas24(3,3),xyzzyaaat24(3),xyzzyaaau24(3),xyzzyaaav24(3)&
&,xyzzyaaaw24(3),xyzzyaaax24,xyzzyaaay24,xyzzyaaaz24,xyzzyaaba24,xyzzy&
&aabb24,xyzzyaabc24(3),xyzzyaabd24(3),xyzzyaabe24,xyzzyaabf24,xyzzyaab&
&g24,xyzzyaabh24,xyzzyaabi24,xyzzyaabj24,xyzzyaabk24,xyzzyaabl24,xyzzy&
&aabm24,xyzzyaabn24,xyzzyaabo24,xyzzyaabp24,xyzzyaabq24
real(dp),parameter :: xyzzyaabr24=1.d-12
real(dp),allocatable :: xyzzyaabs24(:)
logical xyzzyaabt24
logical,parameter :: xyzzyaabu24=.false.
character(59) char59
if((any(k_offset/=0.d0).or.mc_twist_av).and..not.xyzzyaabl1)call errst&
&op_master('HEG_FLUID_SETUP','Must use a complex wave function if the &
&k-vector offset is nonzero or is to be averaged over.')
allocate(xyzzyaaai1(nspin),stat=xyzzyaaab24)
call check_alloc(xyzzyaaab24,'HEG_FLUID_SETUP','nk_gs')
if(xyzzyaabl1)then
xyzzyaaai1=heg_nele
else
xyzzyaaai1=(heg_nele+1)/2
endif
xyzzyaaaf1=maxval(xyzzyaaai1)
xyzzyaado1=maxval(heg_nele)
xyzzyaadp1=free_norb-xyzzyaado1
if(excite_heg)then
if(any(wf_d(1,:,:,:)>1).or.any(wf_d(3,:,:,:)>1).or.any(wf_d(1,:,:,:)<0&
&).or.any(wf_d(3,:,:,:)<0))call errstop_master('HEG_FLUID_SETUP','Band&
& indices in excited state calculations must equal 1 (where in e.g. DE&
&T 1 1 PR X Y X'' Y'' the band indices are X and X'').  You should spe&
&cify the virtual orbitals using the k index in the MDET block in corr&
&elation.data.')
do xyzzyaaad24=1,nspin
if(any(wf_d(2,:,:,xyzzyaaad24)>heg_nele(xyzzyaaad24)))call errstop_mas&
&ter('HEG_FLUID_SETUP','Orbital to excite from is not occupied in the &
&ground state.')
do xyzzyaaae24=1,ndet
do xyzzyaaaf24=1,mdet_max_mods-1
if((wf_d(4,xyzzyaaaf24,xyzzyaaae24,xyzzyaaad24)<=heg_nele(xyzzyaaad24)&
&.and.wf_d(4,xyzzyaaaf24,xyzzyaaae24,xyzzyaaad24)>0)) call errstop_mas&
&ter('HEG_FLUID_SETUP','Virtual orbital already occupied in the ground&
& state.')
if(wf_d(2,xyzzyaaaf24,xyzzyaaae24,xyzzyaaad24)>0)then
if(any(wf_d(2,xyzzyaaaf24+1:mdet_max_mods,xyzzyaaae24,xyzzyaaad24)==wf&
&_d(2,xyzzyaaaf24,xyzzyaaae24,xyzzyaaad24))) call errstop_master('HEG_&
&FLUID_SETUP','You are attempting to promote two electrons from the sa&
&me orbital.')
endif
if(wf_d(4,xyzzyaaaf24,xyzzyaaae24,xyzzyaaad24)>0)then
if(any(wf_d(4,xyzzyaaaf24+1:mdet_max_mods,xyzzyaaae24,xyzzyaaad24)==wf&
&_d(4,xyzzyaaaf24,xyzzyaaae24,xyzzyaaad24))) call errstop_master('HEG_&
&FLUID_SETUP','You are attempting to promote two electrons into the sa&
&me orbital.')
endif
enddo
enddo
enddo
allocate(xyzzyaadq1(mdet_max_mods*ndet*nspin),stat=xyzzyaaab24)
call check_alloc(xyzzyaaab24,'HEG_FLUID_SETUP','')
allocate(xyzzyaadr1(mdet_max_mods*ndet*nspin),xyzzyaads1(mdet_max_mods&
&*ndet*nspin),stat=xyzzyaaab24)
call check_alloc(xyzzyaaab24,'HEG_FLUID_SETUP','')
xyzzyaaag24=0
do xyzzyaaad24=1,nspin
do xyzzyaaae24=1,ndet
if(heg_orbtype(xyzzyaaad24,xyzzyaaae24)/=1)cycle
do xyzzyaaaf24=1,wf_nd(xyzzyaaae24,xyzzyaaad24)
xyzzyaaag24=xyzzyaaag24+1
xyzzyaadq1(xyzzyaaag24)=wf_d(2,xyzzyaaaf24,xyzzyaaae24,xyzzyaaad24)
xyzzyaabt24=.false.
do xyzzyaaah24=1,xyzzyaadp1
if(xyzzyaads1(xyzzyaaah24)==wf_d(4,xyzzyaaaf24,xyzzyaaae24,xyzzyaaad24&
&))then
xyzzyaabt24=.true.
exit
endif
enddo
if(xyzzyaabt24)then
xyzzyaadr1(xyzzyaaag24)=xyzzyaaah24
else
xyzzyaadp1=xyzzyaadp1+1
xyzzyaadr1(xyzzyaaag24)=xyzzyaadp1
xyzzyaads1(xyzzyaadp1)=wf_d(4,xyzzyaaaf24,xyzzyaaae24,xyzzyaaad24)
endif
enddo
enddo
enddo
xyzzyaaaa24=xyzzyaadp1
allocate(xyzzyaadt1(3,xyzzyaaaa24),xyzzyaadu1(xyzzyaaaa24),stat=xyzzya&
&aab24)
call check_alloc(xyzzyaaab24,'HEG_FLUID_SETUP','')
if(use_backflow)then
allocate(xyzzyaadv1(6,xyzzyaaaa24),stat=xyzzyaaab24)
call check_alloc(xyzzyaaab24,'HEG_FLUID_SETUP','')
endif
if(.not.xyzzyaabl1)then
allocate(xyzzyaady1(xyzzyaaaa24),stat=xyzzyaaab24)
call check_alloc(xyzzyaaab24,'HEG_FLUID_SETUP','')
endif
elseif(xyzzyaadx1)then
allocate(xyzzyaads1(xyzzyaadp1),xyzzyaadt1(3,xyzzyaadp1),xyzzyaadu1(xy&
&zzyaadp1),xyzzyaady1(xyzzyaadp1),stat=xyzzyaaab24)
call check_alloc(xyzzyaaab24,'HEG_FLUID_SETUP','virtual_k*')
xyzzyaads1=(/(xyzzyaaaa24,xyzzyaaaa24=xyzzyaado1+1,xyzzyaado1+xyzzyaad&
&p1)/)
if(use_backflow)then
allocate(xyzzyaadv1(6,xyzzyaadp1),stat=xyzzyaaab24)
call check_alloc(xyzzyaaab24,'HEG_FLUID_SETUP','virtual_kprod (in have&
&_virtual)')
endif
endif
allocate(gs_kvec(3,xyzzyaaaf1),xyzzyaaaw1(xyzzyaaaf1),xyzzyaaam1(3,xyz&
&zyaaaf1),stat=xyzzyaaab24)
call check_alloc(xyzzyaaab24,'HEG_FLUID_SETUP','gs_kvec_int,...')
if(use_backflow)then
allocate(xyzzyaabb1(6,xyzzyaaaf1),stat=xyzzyaaab24)
call check_alloc(xyzzyaaab24,'HEG_FLUID_SETUP','gs_kprod')
endif
allocate(xyzzyaaeb1(3,xyzzyaaed1),xyzzyaaec1((xyzzyaaed1+3)/2),xyzzyaa&
&dz1(3,xyzzyaaed1),xyzzyaaao24((xyzzyaaed1+3)/2),stat=xyzzyaaab24)
call check_alloc(xyzzyaaab24,'HEG_FLUID_SETUP','')
if(isperiodic)then
call lattice_generator(xyzzyaaed1,periodicity,b1,b2,b3,xyzzyaaeb1,xyzz&
&yaaec1,xyzzyaaao24,xyzzyaaas24,xyzzyaadz1,xyzzyaaai24)
else
xyzzyaaai24=1
xyzzyaaeb1=0.d0
xyzzyaadz1=0
xyzzyaaao24(1)=1
xyzzyaaao24(2)=2
endif
deallocate(xyzzyaaec1)
if(allocated(xyzzyaaah1))deallocate(xyzzyaaah1,xyzzyaaan1)
allocate(xyzzyaaah1(xyzzyaaai24),xyzzyaaan1(xyzzyaaai24),stat=xyzzyaaa&
&b24)
call check_alloc(xyzzyaaab24,'HEG_FLUID_SETUP','last_k')
xyzzyaaah1(1:xyzzyaaai24)=xyzzyaaao24(2:xyzzyaaai24+1)-1
xyzzyaaan1(1:xyzzyaaai24)=(xyzzyaaah1(1:xyzzyaaai24)+1)/2
deallocate(xyzzyaaao24)
xyzzyaaae1=0
if(xyzzyaabl1)then
allocate(xyzzyaaea1(xyzzyaaed1),xyzzyaaec1(xyzzyaaed1),stat=xyzzyaaab2&
&4)
call check_alloc(xyzzyaaab24,'HEG_FLUID_SETUP','k index')
if(.not.mc_twist_av)then
do xyzzyaaac24=1,xyzzyaaed1
xyzzyaaec1(xyzzyaaac24)=sum((xyzzyaaeb1(1:3,xyzzyaaac24)+k_offset)**2)
enddo
call quicksort(xyzzyaaed1,xyzzyaaec1(1),xyzzyaaea1(1))
allocate(xyzzyaabs24(xyzzyaaed1),xyzzyaaap24(xyzzyaaed1),xyzzyaaaq24(x&
&yzzyaaed1),stat=xyzzyaaab24)
call check_alloc(xyzzyaaab24,'HEG_FLUID_SETUP','kangle')
xyzzyaaaa24=1
do while(xyzzyaaaa24<xyzzyaaed1)
xyzzyaabk24=xyzzyaaec1(xyzzyaaea1(xyzzyaaaa24))
xyzzyaabl24=xyzzyaabr24*xyzzyaabk24
xyzzyaabq24=xyzzyaabk24+xyzzyaabk24
xyzzyaaak24=1
do xyzzyaaal24=xyzzyaaaa24+1,xyzzyaaed1
if(xyzzyaaec1(xyzzyaaea1(xyzzyaaal24))-xyzzyaabk24>xyzzyaabl24)exit
xyzzyaaak24=xyzzyaaak24+1
enddo
if(xyzzyaaak24>1)then
if(periodicity==1)then
do xyzzyaaal24=1,xyzzyaaak24
xyzzyaabs24(xyzzyaaal24)=-xyzzyaaeb1(1,xyzzyaaea1(xyzzyaaaa24+xyzzyaaa&
&l24-1))
enddo
elseif(periodicity==2)then
do xyzzyaaal24=1,xyzzyaaak24
xyzzyaaam24=xyzzyaaea1(xyzzyaaaa24+xyzzyaaal24-1)
xyzzyaabg24=xyzzyaaeb1(1,xyzzyaaam24)+k_offset(1)
xyzzyaabh24=xyzzyaaeb1(2,xyzzyaaam24)+k_offset(2)
if(xyzzyaabh24>=0.d0)then
if(xyzzyaabg24>=0.d0)then
xyzzyaabs24(xyzzyaaal24)=xyzzyaabh24*xyzzyaabh24-xyzzyaabk24
else
xyzzyaabs24(xyzzyaaal24)=xyzzyaabg24*xyzzyaabg24
endif
else
if(xyzzyaabg24<0.d0)then
xyzzyaabs24(xyzzyaaal24)=xyzzyaabk24+xyzzyaabh24*xyzzyaabh24
else
xyzzyaabs24(xyzzyaaal24)=xyzzyaabq24+xyzzyaabg24*xyzzyaabg24
endif
endif
enddo
else
do xyzzyaaal24=1,xyzzyaaak24
xyzzyaabs24(xyzzyaaal24)=xyzzyaaeb1(3,xyzzyaaea1(xyzzyaaaa24+xyzzyaaal&
&24-1))
enddo
call quicksort(xyzzyaaak24,xyzzyaabs24(1),xyzzyaaap24(1))
xyzzyaabp24=2.d0*abs(xyzzyaabs24(xyzzyaaap24(xyzzyaaak24)))
xyzzyaabo24=xyzzyaabr24*abs(xyzzyaabs24(xyzzyaaap24(xyzzyaaak24)))
xyzzyaabi24=xyzzyaabs24(xyzzyaaap24(1))
do xyzzyaaal24=1,xyzzyaaak24-1
xyzzyaabj24=xyzzyaabs24(xyzzyaaap24(xyzzyaaal24+1))
xyzzyaabn24=xyzzyaabj24-xyzzyaabi24
if(xyzzyaabn24<xyzzyaabp24.and.xyzzyaabn24>xyzzyaabo24)xyzzyaabp24=xyz&
&zyaabn24
xyzzyaabi24=xyzzyaabj24
enddo
if(xyzzyaabp24>0.d0)then
xyzzyaabm24=5.d0*xyzzyaabk24/xyzzyaabp24
else
xyzzyaabm24=0.d0
endif
do xyzzyaaal24=1,xyzzyaaak24
xyzzyaaam24=xyzzyaaea1(xyzzyaaaa24+xyzzyaaal24-1)
xyzzyaabg24=xyzzyaaeb1(1,xyzzyaaam24)+k_offset(1)
xyzzyaabh24=xyzzyaaeb1(2,xyzzyaaam24)+k_offset(2)
xyzzyaabi24=xyzzyaaeb1(3,xyzzyaaam24)+k_offset(3)
if(xyzzyaabh24>=0.d0)then
if(xyzzyaabg24>=0.d0)then
xyzzyaabs24(xyzzyaaal24)=xyzzyaabh24*xyzzyaabh24-xyzzyaabk24
else
xyzzyaabs24(xyzzyaaal24)=xyzzyaabg24*xyzzyaabg24
endif
else
if(xyzzyaabg24<0.d0)then
xyzzyaabs24(xyzzyaaal24)=xyzzyaabk24+xyzzyaabh24*xyzzyaabh24
else
xyzzyaabs24(xyzzyaaal24)=xyzzyaabq24+xyzzyaabg24*xyzzyaabg24
endif
endif
xyzzyaabs24(xyzzyaaal24)=xyzzyaabs24(xyzzyaaal24)-xyzzyaabm24*xyzzyaab&
&i24
enddo
endif
call quicksort(xyzzyaaak24,xyzzyaabs24(1),xyzzyaaap24(1))
do xyzzyaaal24=1,xyzzyaaak24
xyzzyaaaq24(xyzzyaaaa24+xyzzyaaal24-1)=xyzzyaaea1(xyzzyaaaa24+xyzzyaaa&
&p24(xyzzyaaal24)-1)
enddo
xyzzyaaea1(xyzzyaaaa24:xyzzyaaaa24+xyzzyaaak24-1)=xyzzyaaaq24(xyzzyaaa&
&a24:xyzzyaaaa24+xyzzyaaak24-1)
endif
xyzzyaaaa24=xyzzyaaaa24+xyzzyaaak24
enddo
deallocate(xyzzyaabs24,xyzzyaaap24,xyzzyaaaq24)
do xyzzyaaaa24=1,xyzzyaaaf1
xyzzyaaac24=xyzzyaaea1(xyzzyaaaa24)
xyzzyaaam1(1:3,xyzzyaaaa24)=xyzzyaadz1(1:3,xyzzyaaac24)
gs_kvec(1:3,xyzzyaaaa24)=xyzzyaaeb1(1:3,xyzzyaaac24)+k_offset
xyzzyaaaw1(xyzzyaaaa24)=xyzzyaaec1(xyzzyaaac24)
if(use_backflow)call xyzzyaafv1(gs_kvec(1:3,xyzzyaaaa24),xyzzyaabb1(1:&
&6,xyzzyaaaa24))
if(abs(xyzzyaaam1(1,xyzzyaaaa24))>xyzzyaaae1(1))xyzzyaaae1(1)=abs(xyzz&
&yaaam1(1,xyzzyaaaa24))
if(abs(xyzzyaaam1(2,xyzzyaaaa24))>xyzzyaaae1(2))xyzzyaaae1(2)=abs(xyzz&
&yaaam1(2,xyzzyaaaa24))
if(abs(xyzzyaaam1(3,xyzzyaaaa24))>xyzzyaaae1(3))xyzzyaaae1(3)=abs(xyzz&
&yaaam1(3,xyzzyaaaa24))
enddo
if(excite_heg.or.xyzzyaadx1)then
do xyzzyaaaa24=1,xyzzyaadp1
xyzzyaaac24=xyzzyaaea1(xyzzyaads1(xyzzyaaaa24))
xyzzyaadt1(1:3,xyzzyaaaa24)=xyzzyaaeb1(1:3,xyzzyaaac24)+k_offset
xyzzyaadu1(xyzzyaaaa24)=xyzzyaaec1(xyzzyaaac24)
if(use_backflow)call xyzzyaafv1(xyzzyaadt1(1:3,xyzzyaaaa24),xyzzyaadv1&
&(1:6,xyzzyaaaa24))
enddo
endif
if(excite_heg.and.xyzzyaabu24.and.am_master)then
call wout('Detailed excitation information')
do xyzzyaaad24=1,nspin
if(nspin>1.and.any(wf_nd(:,xyzzyaaad24)>0))call wout('SPIN '//trim(i2s&
&(xyzzyaaad24)))
do xyzzyaaae24=1,ndet
if(ndet>1.and.wf_nd(xyzzyaaae24,xyzzyaaad24)>0)call wout(' DETERMINANT&
& ' //trim(i2s(xyzzyaaae24)))
do xyzzyaaaf24=1,wf_nd(xyzzyaaae24,xyzzyaaad24)
if(wf_nd(xyzzyaaae24,xyzzyaaad24)>1)call wout('  EXCITATION '//trim(i2&
&s(xyzzyaaaf24)))
xyzzyaabd24=xyzzyaaeb1(1:3,xyzzyaaea1(wf_d(2,xyzzyaaaf24,xyzzyaaae24,x&
&yzzyaaad24)))+k_offset
xyzzyaabf24=sqrt(xyzzyaaec1(xyzzyaaea1(wf_d(2,xyzzyaaaf24,xyzzyaaae24,&
&xyzzyaaad24))))
write(char59,'(es19.12,1x,es19.12,1x,es19.12)')xyzzyaabd24
call wout('   PROMO. fm k  = '//char59)
call wout('            |k| = ',xyzzyaabf24)
xyzzyaabc24=xyzzyaaeb1(1:3,xyzzyaaea1(wf_d(4,xyzzyaaaf24,xyzzyaaae24,x&
&yzzyaaad24)))+k_offset
xyzzyaabe24=sqrt(xyzzyaaec1(xyzzyaaea1(wf_d(4,xyzzyaaaf24,xyzzyaaae24,&
&xyzzyaaad24))))
write(char59,'(es19.12,1x,es19.12,1x,es19.12)')xyzzyaabc24
call wout('          to k  = '//char59)
call wout('            |k| = ',xyzzyaabe24)
if(xyzzyaabe24>0.d0.and.xyzzyaabf24>0.d0)call wout('   Angle between e&
&lectron and hole: ',acos(dot_product(xyzzyaabc24,xyzzyaabd24)/(xyzzya&
&abe24*xyzzyaabf24)))
enddo
enddo
enddo
call wout()
endif
deallocate(xyzzyaaeb1,xyzzyaadz1,xyzzyaaec1)
allocate(xyzzyaabc1(-xyzzyaaae1(1):xyzzyaaae1(1)),xyzzyaabd1(-xyzzyaaa&
&e1(2):xyzzyaaae1(2)),xyzzyaabe1(-xyzzyaaae1(3):xyzzyaaae1(3)),stat=xy&
&zzyaaab24)
call check_alloc(xyzzyaaab24,'HEG_FLUID_SETUP','EXPIKDOTR_B arrays')
else
do xyzzyaaac24=1,xyzzyaaed1
xyzzyaaea1(xyzzyaaac24)=xyzzyaaac24
enddo
if(periodicity==3)then
xyzzyaaat24=b1+b2+b3
xyzzyaaau24=b1+b2-b3
xyzzyaaav24=b1-b2+b3
xyzzyaaaw24=b1-b2-b3
xyzzyaaay24=0.5d0*sqrt(max(dot_product(xyzzyaaat24,xyzzyaaat24),dot_pr&
&oduct(xyzzyaaau24,xyzzyaaau24),dot_product(xyzzyaaav24,xyzzyaaav24),d&
&ot_product(xyzzyaaaw24,xyzzyaaaw24)))
elseif(periodicity==2)then
xyzzyaaat24=b1+b2
xyzzyaaau24=b1-b2
xyzzyaaay24=0.5d0*sqrt(max(dot_product(xyzzyaaat24,xyzzyaaat24),dot_pr&
&oduct(xyzzyaaau24,xyzzyaaau24)))
else
xyzzyaaay24=0.5d0*abs(b1(1))
endif
xyzzyaaba24=sqrt(dot_product(xyzzyaaeb1(1:3,xyzzyaaaf1),xyzzyaaeb1(1:3&
&,xyzzyaaaf1)))
xyzzyaaax24=(xyzzyaaba24+xyzzyaaay24)**2
xyzzyaabb24=0.d0
xyzzyaaaj24=0
xyzzyaaag1=0
do xyzzyaaaa24=1,xyzzyaaed1
xyzzyaaaz24=dot_product(xyzzyaaeb1(1:3,xyzzyaaaa24),xyzzyaaeb1(1:3,xyz&
&zyaaaa24))
if(xyzzyaaaz24>=xyzzyaaax24)exit
if(abs(xyzzyaadz1(1,xyzzyaaaa24))>xyzzyaaaj24(1))xyzzyaaaj24(1)=abs(xy&
&zzyaadz1(1,xyzzyaaaa24))
if(abs(xyzzyaadz1(2,xyzzyaaaa24))>xyzzyaaaj24(2))xyzzyaaaj24(2)=abs(xy&
&zzyaadz1(2,xyzzyaaaa24))
if(abs(xyzzyaadz1(3,xyzzyaaaa24))>xyzzyaaaj24(3))xyzzyaaaj24(3)=abs(xy&
&zzyaadz1(3,xyzzyaaaa24))
if(xyzzyaaaz24<=xyzzyaabb24)then
xyzzyaaee1=xyzzyaaaa24+1
xyzzyaaag1=xyzzyaaaj24
endif
enddo
xyzzyaaef1=xyzzyaaaa24-1
allocate(xyzzyaabc1(-xyzzyaaaj24(1):xyzzyaaaj24(1)),xyzzyaabd1(-xyzzya&
&aaj24(2):xyzzyaaaj24(2)),xyzzyaabe1(-xyzzyaaaj24(3):xyzzyaaaj24(3)),s&
&tat=xyzzyaaab24)
call check_alloc(xyzzyaaab24,'HEG_FLUID_SETUP','EXPIKDOTR_B arrays')
do xyzzyaaaa24=1,xyzzyaaee1-1
xyzzyaaam1(1:3,xyzzyaaaa24)=xyzzyaadz1(1:3,xyzzyaaaa24)
gs_kvec(1:3,xyzzyaaaa24)=xyzzyaaeb1(1:3,xyzzyaaaa24)+k_offset
xyzzyaaaw1(xyzzyaaaa24)=sum(gs_kvec(1:3,xyzzyaaaa24)**2)
if(use_backflow)call xyzzyaafv1(gs_kvec(1:3,xyzzyaaaa24),xyzzyaabb1(1:&
&6,xyzzyaaaa24))
enddo
do xyzzyaaac24=xyzzyaaee1,xyzzyaaef1
xyzzyaaec1(xyzzyaaac24)=sum((xyzzyaaeb1(1:3,xyzzyaaac24)+k_offset)**2)
enddo
call quicksort(xyzzyaaef1-xyzzyaaee1+1,xyzzyaaec1(xyzzyaaee1),xyzzyaae&
&a1(xyzzyaaee1))
xyzzyaaea1(xyzzyaaee1:xyzzyaaef1)=xyzzyaaea1(xyzzyaaee1:xyzzyaaef1)+(x&
&yzzyaaee1-1)
xyzzyaaae1=xyzzyaaag1
do xyzzyaaaa24=xyzzyaaee1,xyzzyaaaf1
xyzzyaaac24=xyzzyaaea1(xyzzyaaaa24)
xyzzyaaam1(1:3,xyzzyaaaa24)=xyzzyaadz1(1:3,xyzzyaaac24)
xyzzyaaaw1(xyzzyaaaa24)=xyzzyaaec1(xyzzyaaac24)
gs_kvec(1:3,xyzzyaaaa24)=xyzzyaaeb1(1:3,xyzzyaaac24)+k_offset
if(use_backflow)call xyzzyaafv1(gs_kvec(1:3,xyzzyaaaa24),xyzzyaabb1(1:&
&6,xyzzyaaaa24))
if(abs(xyzzyaaam1(1,xyzzyaaaa24))>xyzzyaaae1(1))xyzzyaaae1(1)=abs(xyzz&
&yaaam1(1,xyzzyaaaa24))
if(abs(xyzzyaaam1(2,xyzzyaaaa24))>xyzzyaaae1(2))xyzzyaaae1(2)=abs(xyzz&
&yaaam1(2,xyzzyaaaa24))
if(abs(xyzzyaaam1(3,xyzzyaaaa24))>xyzzyaaae1(3))xyzzyaaae1(3)=abs(xyzz&
&yaaam1(3,xyzzyaaaa24))
enddo
if(excite_heg)then
do xyzzyaaaa24=1,xyzzyaadp1
xyzzyaaac24=xyzzyaaea1(xyzzyaads1(xyzzyaaaa24))
xyzzyaadu1(xyzzyaaaa24)=xyzzyaaec1(xyzzyaaac24)
xyzzyaadt1(1:3,xyzzyaaaa24)=xyzzyaaeb1(1:3,xyzzyaaac24)+k_offset
if(use_backflow)call xyzzyaafv1(xyzzyaadt1(1:3,xyzzyaaaa24),xyzzyaadv1&
&(1:6,xyzzyaaaa24))
enddo
endif
endif
else
xyzzyaaam1(1:3,1)=0
gs_kvec(1:3,1)=0.d0
xyzzyaaaw1(1)=0.d0
if(use_backflow)xyzzyaabb1(1:3,1)=0.d0
xyzzyaaac24=0
do xyzzyaaaa24=2,xyzzyaaaf1
xyzzyaaac24=xyzzyaaac24+2
if(xyzzyaadz1(1,xyzzyaaac24)<0)then
xyzzyaaam1(1:3,xyzzyaaaa24)=-xyzzyaadz1(1:3,xyzzyaaac24)
gs_kvec(1:3,xyzzyaaaa24)=-xyzzyaaeb1(1:3,xyzzyaaac24)
else
xyzzyaaam1(1:3,xyzzyaaaa24)=xyzzyaadz1(1:3,xyzzyaaac24)
gs_kvec(1:3,xyzzyaaaa24)=xyzzyaaeb1(1:3,xyzzyaaac24)
endif
xyzzyaaaw1(xyzzyaaaa24)=sum(gs_kvec(1:3,xyzzyaaaa24)**2)
if(use_backflow)call xyzzyaafv1(gs_kvec(1:3,xyzzyaaaa24),xyzzyaabb1(1:&
&6,xyzzyaaaa24))
if(xyzzyaaam1(1,xyzzyaaaa24)>xyzzyaaae1(1))xyzzyaaae1(1)=xyzzyaaam1(1,&
&xyzzyaaaa24)
if(abs(xyzzyaaam1(2,xyzzyaaaa24))>xyzzyaaae1(2))xyzzyaaae1(2)=abs(xyzz&
&yaaam1(2,xyzzyaaaa24))
if(abs(xyzzyaaam1(3,xyzzyaaaa24))>xyzzyaaae1(3))xyzzyaaae1(3)=abs(xyzz&
&yaaam1(3,xyzzyaaaa24))
enddo
if(xyzzyaadx1)then
do xyzzyaaaa24=1,xyzzyaadp1
xyzzyaaac24=xyzzyaads1(xyzzyaaaa24)
if(xyzzyaadz1(1,xyzzyaaac24)<0)then
xyzzyaadt1(1:3,xyzzyaaaa24)=-xyzzyaaeb1(1:3,xyzzyaaac24)
else
xyzzyaadt1(1:3,xyzzyaaaa24)=xyzzyaaeb1(1:3,xyzzyaaac24)
endif
xyzzyaadu1(xyzzyaaaa24)=sum(xyzzyaadt1(1:3,xyzzyaaaa24)**2)
if(use_backflow)call xyzzyaafv1(xyzzyaadt1(1:3,xyzzyaaaa24),xyzzyaadv1&
&(1:6,xyzzyaaaa24))
xyzzyaady1(xyzzyaaaa24)=(mod(xyzzyaaac24,2)==1)
enddo
endif
deallocate(xyzzyaaeb1,xyzzyaadz1)
allocate(xyzzyaabc1(0:xyzzyaaae1(1)),xyzzyaabd1(-xyzzyaaae1(2):xyzzyaa&
&ae1(2)),xyzzyaabe1(-xyzzyaaae1(3):xyzzyaaae1(3)),stat=xyzzyaaab24)
call check_alloc(xyzzyaaab24,'HEG_FLUID_SETUP','EXPIKDOTR_B arrays')
endif
allocate(xyzzyaaar24(xyzzyaado1+xyzzyaadp1),stat=xyzzyaaab24)
call check_alloc(xyzzyaaab24,'HEG_FLUID_SETUP','orb_kvec_idx')
xyzzyaaar24(:)=0
if(xyzzyaabl1)then
do xyzzyaaac24=1,xyzzyaado1
xyzzyaaar24(xyzzyaaac24)=xyzzyaaac24
enddo
else
xyzzyaaar24(1)=1
xyzzyaaan24=1
do xyzzyaaac24=2,xyzzyaaaf1
xyzzyaaan24=xyzzyaaan24+1
xyzzyaaar24(xyzzyaaan24)=xyzzyaaac24
if(xyzzyaaan24==xyzzyaado1)exit
xyzzyaaan24=xyzzyaaan24+1
xyzzyaaar24(xyzzyaaan24)=xyzzyaaac24
enddo
endif
if(xyzzyaadx1)then
do xyzzyaaan24=1,xyzzyaadp1
if(xyzzyaabl1)then
xyzzyaaar24(xyzzyaaan24)=xyzzyaaea1(xyzzyaads1(xyzzyaaan24))
else
xyzzyaaar24(xyzzyaaan24)=xyzzyaads1(xyzzyaaan24)
endif
enddo
if(excite_heg)then
xyzzyaaag24=0
do xyzzyaaad24=1,nspin
do xyzzyaaae24=1,ndet
do xyzzyaaaf24=1,wf_nd(xyzzyaaae24,xyzzyaaad24)
xyzzyaaag24=xyzzyaaag24+1
xyzzyaaan24=xyzzyaadr1(xyzzyaaag24)
if(.not.any(xyzzyaaao1==xyzzyaado1+xyzzyaaan24))free_norb=free_norb+1
xyzzyaaao1(xyzzyaadq1(xyzzyaaag24),xyzzyaaad24,xyzzyaaae24)=xyzzyaado1&
&+xyzzyaaan24
enddo
enddo
enddo
endif
endif
if(xyzzyaado1+xyzzyaadp1/=free_norb)call errstop('HEG_FLUID_SETUP','Di&
&sagreement in total number of orbitals.')
allocate(xyzzyaafh1(xyzzyaaff1,free_norb),stat=xyzzyaaab24)
call check_alloc(xyzzyaaab24,'HEG_FLUID_SETUP','free_orbdesc_int')
xyzzyaafh1(1,1:free_norb)=xyzzyaaar24(1:free_norb)
if(xyzzyaabl1.and..not.mc_twist_av)deallocate(xyzzyaaea1)
end subroutine xyzzyaafu1
subroutine mc_twist_offset(k_offset_in)
use slaarnabq, only : min_image_brute_force
use slaarnacc, only : ranx
implicit none
real(dp),intent(in),optional :: k_offset_in(3)
integer xyzzyaaaa25,xyzzyaaab25
real(dp) xyzzyaaac25(3)
if(.not.present(k_offset_in))then
if(dimensionality==3)then
xyzzyaaac25=ranx()*b1+ranx()*b2+ranx()*b3
elseif(dimensionality==2)then
xyzzyaaac25=ranx()*b1+ranx()*b2
else
xyzzyaaac25=ranx()*b1
endif
call min_image_brute_force(periodicity,xyzzyaaac25,xyzzyaacs1,binv,k_o&
&ffset,xyzzyaact1)
else
k_offset=k_offset_in
endif
do xyzzyaaab25=1,xyzzyaaee1-1
gs_kvec(1:3,xyzzyaaab25)=xyzzyaaeb1(1:3,xyzzyaaab25)+k_offset
xyzzyaaaw1(xyzzyaaab25)=sum(gs_kvec(1:3,xyzzyaaab25)**2)
if(use_backflow)call xyzzyaafv1(gs_kvec(1:3,xyzzyaaab25),xyzzyaabb1(1:&
&6,xyzzyaaab25))
enddo
do xyzzyaaaa25=xyzzyaaee1,xyzzyaaef1
xyzzyaaec1(xyzzyaaaa25)=sum((xyzzyaaeb1(1:3,xyzzyaaaa25)+k_offset)**2)
enddo
call quicksort(xyzzyaaef1-xyzzyaaee1+1,xyzzyaaec1(xyzzyaaee1),xyzzyaae&
&a1(xyzzyaaee1))
xyzzyaaea1(xyzzyaaee1:xyzzyaaef1) =xyzzyaaea1(xyzzyaaee1:xyzzyaaef1)+(&
&xyzzyaaee1-1)
xyzzyaaae1=xyzzyaaag1
do xyzzyaaab25=xyzzyaaee1,xyzzyaaaf1
xyzzyaaaa25=xyzzyaaea1(xyzzyaaab25)
xyzzyaaam1(1:3,xyzzyaaab25)=xyzzyaadz1(1:3,xyzzyaaaa25)
xyzzyaaaw1(xyzzyaaab25)=xyzzyaaec1(xyzzyaaaa25)
gs_kvec(1:3,xyzzyaaab25)=xyzzyaaeb1(1:3,xyzzyaaaa25)+k_offset
if(use_backflow)call xyzzyaafv1(gs_kvec(1:3,xyzzyaaab25),xyzzyaabb1(1:&
&6,xyzzyaaab25))
if(abs(xyzzyaaam1(1,xyzzyaaab25))>xyzzyaaae1(1))xyzzyaaae1(1)=abs(xyzz&
&yaaam1(1,xyzzyaaab25))
if(abs(xyzzyaaam1(2,xyzzyaaab25))>xyzzyaaae1(2))xyzzyaaae1(2)=abs(xyzz&
&yaaam1(2,xyzzyaaab25))
if(abs(xyzzyaaam1(3,xyzzyaaab25))>xyzzyaaae1(3))xyzzyaaae1(3)=abs(xyzz&
&yaaam1(3,xyzzyaaab25))
enddo
if(excite_heg)then
do xyzzyaaab25=1,xyzzyaadp1
xyzzyaaaa25=xyzzyaaea1(xyzzyaads1(xyzzyaaab25))
xyzzyaadu1(xyzzyaaab25)=xyzzyaaec1(xyzzyaaaa25)
xyzzyaadt1(1:3,xyzzyaaab25)=xyzzyaaeb1(1:3,xyzzyaaaa25)+k_offset
if(use_backflow)call xyzzyaafv1(xyzzyaadt1(1:3,xyzzyaaab25),xyzzyaadv1&
&(1:6,xyzzyaaab25))
enddo
endif
end subroutine mc_twist_offset
subroutine ppmcta_partial_sum(n,e,k,x,partial_sum,w)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: e(n),k(n),x(n)
real(dp),intent(in),optional :: w(n)
real(dp),intent(out) :: partial_sum(9)
if(present(w))then
partial_sum(1)=sum(w)
partial_sum(2)=sum(w*e)
partial_sum(3)=sum(w*k)
partial_sum(4)=sum(w*k**2)
partial_sum(5)=sum(w*e*k)
partial_sum(6)=sum(w*x)
partial_sum(7)=sum(w*x**2)
partial_sum(8)=sum(w*e*x)
partial_sum(9)=sum(w*k*x)
else
partial_sum(1)=dble(n)
partial_sum(2)=sum(e)
partial_sum(3)=sum(k)
partial_sum(4)=sum(k**2)
partial_sum(5)=sum(e*k)
partial_sum(6)=sum(x)
partial_sum(7)=sum(x**2)
partial_sum(8)=sum(e*x)
partial_sum(9)=sum(k*x)
endif
end subroutine ppmcta_partial_sum
subroutine ppmcta_fit(partial_sum,a,b)
implicit none
real(dp),intent(in) :: partial_sum(9)
real(dp),intent(out) :: a,b
real(dp) xyzzyaaaa27,xyzzyaaab27,xyzzyaaac27,xyzzyaaad27,xyzzyaaae27,x&
&yzzyaaaf27,xyzzyaaag27,xyzzyaaah27,xyzzyaaai27,xyzzyaaaj27,xyzzyaaak2&
&7,xyzzyaaal27,xyzzyaaam27,xyzzyaaan27,xyzzyaaao27
logical xyzzyaaap27
xyzzyaaap27=partial_sum(7)/=0.d0
xyzzyaaaa27=1.d0/partial_sum(1)
xyzzyaaab27=partial_sum(2)*xyzzyaaaa27
xyzzyaaac27=partial_sum(3)*xyzzyaaaa27
xyzzyaaad27=partial_sum(4)*xyzzyaaaa27
xyzzyaaae27=partial_sum(5)*xyzzyaaaa27
if(xyzzyaaap27)then
xyzzyaaaf27=partial_sum(6)*xyzzyaaaa27
xyzzyaaag27=partial_sum(7)*xyzzyaaaa27
xyzzyaaah27=partial_sum(8)*xyzzyaaaa27
xyzzyaaai27=partial_sum(9)*xyzzyaaaa27
endif
xyzzyaaaj27=xyzzyaaad27-xyzzyaaac27**2
xyzzyaaal27=xyzzyaaae27-xyzzyaaab27*xyzzyaaac27
if(xyzzyaaap27)then
xyzzyaaak27=xyzzyaaag27-xyzzyaaaf27**2
xyzzyaaam27=xyzzyaaah27-xyzzyaaab27*xyzzyaaaf27
xyzzyaaan27=xyzzyaaai27-xyzzyaaac27*xyzzyaaaf27
xyzzyaaao27=xyzzyaaan27**2-xyzzyaaaj27*xyzzyaaak27
if(xyzzyaaao27==0.d0)then
a=-1.d0
b=0.d0
else
a=(xyzzyaaak27*xyzzyaaal27-xyzzyaaam27*xyzzyaaan27)/xyzzyaaao27
b=(xyzzyaaaj27*xyzzyaaam27-xyzzyaaal27*xyzzyaaan27)/xyzzyaaao27
endif
else
if(xyzzyaaaj27==0.d0)then
a=-1.d0
else
a=-xyzzyaaal27/xyzzyaaaj27
endif
b=0.d0
endif
end subroutine ppmcta_fit
subroutine xyzzyaafv1(k,kprod)
implicit none
real(dp),intent(in) :: k(3)
real(dp),intent(out) :: kprod(6)
kprod(1)=k(1)*k(1)
kprod(2)=k(2)*k(2)
kprod(3)=k(3)*k(3)
kprod(4)=k(1)*k(2)
kprod(5)=k(1)*k(3)
kprod(6)=k(2)*k(3)
end subroutine xyzzyaafv1
subroutine xyzzyaafw1
use slaarnaan, only : lattice_generator
implicit none
integer,allocatable :: xyzzyaaaa29(:)
integer xyzzyaaab29,xyzzyaaac29,xyzzyaaad29,xyzzyaaae29,xyzzyaaaf29,xy&
&zzyaaag29,xyzzyaaah29
real(dp) xyzzyaaai29(3,3)
allocate(xyzzyaaeb1(3,xyzzyaaed1),xyzzyaaec1((xyzzyaaed1+3)/2),xyzzyaa&
&dz1(3,xyzzyaaed1),xyzzyaaaa29((xyzzyaaed1+3)/2),stat=xyzzyaaab29)
call check_alloc(xyzzyaaab29,'SETUP_HEG_K_LATTICE','')
call lattice_generator(xyzzyaaed1,periodicity,b1,b2,b3,xyzzyaaeb1,xyzz&
&yaaec1,xyzzyaaaa29,xyzzyaaai29,xyzzyaadz1,xyzzyaaah29)
deallocate(xyzzyaadz1)
if(allocated(xyzzyaaah1))deallocate(xyzzyaaah1,xyzzyaaan1)
allocate(xyzzyaaah1(xyzzyaaah29),xyzzyaaan1(xyzzyaaah29),stat=xyzzyaaa&
&b29)
call check_alloc(xyzzyaaab29,'SETUP_HEG_K_LATTICE','last_k')
xyzzyaaah1(1:xyzzyaaah29)=xyzzyaaaa29(2:xyzzyaaah29+1)-1
xyzzyaaan1(1:xyzzyaaah29)=(xyzzyaaah1(1:xyzzyaaah29)+1)/2
if(xyzzyaaab1>xyzzyaaah29)call errstop_master('SETUP_HEG_K_LATTICE','N&
&umber of generated stars too small for required plane-wave expansion.&
&')
xyzzyaaag29=max(maxval(heg_nele),xyzzyaaah1(max(xyzzyaaab1,1)))
xyzzyaaaf29=(xyzzyaaag29+1)/2
allocate(xyzzyaaaz1(3,xyzzyaaaf29),xyzzyaaba1(xyzzyaaaf29+1),stat=xyzz&
&yaaab29)
call check_alloc(xyzzyaaab29,'SETUP_HEG_K_LATTICE','')
xyzzyaaaz1(1:3,1)=0.d0
xyzzyaaac29=1
mu_loop: do xyzzyaaae29=2,xyzzyaaah29
do xyzzyaaad29=xyzzyaaaa29(xyzzyaaae29),xyzzyaaah1(xyzzyaaae29),2
xyzzyaaac29=xyzzyaaac29+1
if(xyzzyaaac29>xyzzyaaaf29)exit mu_loop
xyzzyaaaz1(1:3,xyzzyaaac29)=xyzzyaaeb1(1:3,xyzzyaaad29)
xyzzyaaba1(xyzzyaaac29)=xyzzyaaec1(xyzzyaaae29)
enddo
enddo mu_loop
deallocate(xyzzyaaeb1,xyzzyaaec1,xyzzyaaaa29)
end subroutine xyzzyaafw1
subroutine heg_mpc_long
implicit none
integer xyzzyaaaa30(3)
integer,parameter :: xyzzyaaab30(3)=(/128,128,128/)
real(dp) xyzzyaaac30,xyzzyaaad30,xyzzyaaae30,xyzzyaaaf30,xyzzyaaag30,x&
&yzzyaaah30,xyzzyaaai30,xyzzyaaaj30,xyzzyaaak30,xyzzyaaal30,xyzzyaaam3&
&0,xyzzyaaan30,xyzzyaaao30
call timer('MPC_LONG',.true.)
xyzzyaaaa30=xyzzyaaab30
xyzzyaaag30=1.d0/real(xyzzyaaaa30(1)**2,dp)
call xyzzyaafx1(xyzzyaaaa30,xyzzyaaad30)
xyzzyaaaa30=xyzzyaaaa30*2
xyzzyaaah30=1.d0/real(xyzzyaaaa30(1)**2,dp)
call xyzzyaafx1(xyzzyaaaa30,xyzzyaaae30)
xyzzyaaaa30=xyzzyaaaa30*2
xyzzyaaai30=1.d0/real(xyzzyaaaa30(1)**2,dp)
call xyzzyaafx1(xyzzyaaaa30,xyzzyaaaf30)
xyzzyaaaj30=xyzzyaaah30*xyzzyaaai30/((xyzzyaaag30-xyzzyaaah30)*(xyzzya&
&aag30-xyzzyaaai30))
xyzzyaaak30=xyzzyaaag30*xyzzyaaai30/((xyzzyaaah30-xyzzyaaag30)*(xyzzya&
&aah30-xyzzyaaai30))
xyzzyaaal30=xyzzyaaag30*xyzzyaaah30/((xyzzyaaai30-xyzzyaaag30)*(xyzzya&
&aai30-xyzzyaaah30))
xyzzyaaac30=xyzzyaaaj30*xyzzyaaad30+xyzzyaaak30*xyzzyaaae30+xyzzyaaal3&
&0*xyzzyaaaf30
if(am_master)then
xyzzyaaam30=xyzzyaaai30/(xyzzyaaai30-xyzzyaaah30)
xyzzyaaan30=xyzzyaaah30/(xyzzyaaah30-xyzzyaaai30)
xyzzyaaao30=xyzzyaaam30*xyzzyaaae30+xyzzyaaan30*xyzzyaaaf30
call wout()
call wout('Analysis of extrapolation of Fourier coefficient of f(r):')
call wout()
call wout('    1/[FFT grid(1)]^2             f_{G=0}')
call wout('-------------------------------------------------')
call wout('',(/xyzzyaaag30,xyzzyaaad30/),rfmt='(es24.16)')
call wout('',(/xyzzyaaah30,xyzzyaaae30/),rfmt='(es24.16)')
call wout('',(/xyzzyaaai30,xyzzyaaaf30/),rfmt='(es24.16)')
call wout('     0 (lin. extrap.)   ',xyzzyaaao30,rfmt='(es24.16)')
call wout('     0 (quad. extrap.)  ',xyzzyaaac30,rfmt='(es24.16)')
call wout()
call wout('Please verify that the extrapolated coefficient looks corre&
&ct.')
endif
half_mpc_constant=-xyzzyaaac30*rho*0.5d0
call timer('MPC_LONG',.false.)
end subroutine heg_mpc_long
subroutine xyzzyaafx1(ngrid,fg_0)
use slaarnabq, only : minimum_image
implicit none
integer,intent(in) :: ngrid(3)
real(dp),intent(out) :: fg_0
integer xyzzyaaaa31,xyzzyaaab31,xyzzyaaac31
real(dp) xyzzyaaad31,xyzzyaaae31,xyzzyaaaf31,xyzzyaaag31,xyzzyaaah31,x&
&yzzyaaai31,xyzzyaaaj31,xyzzyaaak31,xyzzyaaal31,xyzzyaaam31(3),xyzzyaa&
&an31
xyzzyaaae31=-0.5d0/wigner_seitz_radius**3
xyzzyaaaf31=1.5d0/wigner_seitz_radius
xyzzyaaan31=wigner_seitz_radius**2
if(periodicity==3)then
xyzzyaaaj31=1.d0/real(ngrid(1),dp)
xyzzyaaak31=1.d0/real(ngrid(2),dp)
xyzzyaaal31=1.d0/real(ngrid(3),dp)
fg_0=0.d0
xyzzyaaai31=0.d0
do xyzzyaaac31=1,ngrid(3)
xyzzyaaai31=xyzzyaaai31+xyzzyaaal31
xyzzyaaah31=0.d0
do xyzzyaaab31=1,ngrid(2)
xyzzyaaah31=xyzzyaaah31+xyzzyaaak31
xyzzyaaag31=0.d0
do xyzzyaaaa31=1,ngrid(1)
xyzzyaaag31=xyzzyaaag31+xyzzyaaaj31
xyzzyaaam31=xyzzyaaag31*a1+xyzzyaaah31*a2+xyzzyaaai31*a3
call minimum_image(3,1,xyzzyaaam31)
xyzzyaaad31=sum(xyzzyaaam31**2)
if(xyzzyaaad31<xyzzyaaan31)then
fg_0=fg_0+xyzzyaaae31*xyzzyaaad31+xyzzyaaaf31
else
fg_0=fg_0+1.d0/sqrt(xyzzyaaad31)
endif
enddo
enddo
enddo
fg_0=fg_0*volume/real(ngrid(1)*ngrid(2)*ngrid(3),dp)+twopi/5.d0*xyzzya&
&aan31
else
xyzzyaaaj31=1.d0/real(ngrid(1),dp)
xyzzyaaak31=1.d0/real(ngrid(2),dp)
fg_0=0.d0
xyzzyaaah31=0.d0
do xyzzyaaab31=1,ngrid(2)
xyzzyaaah31=xyzzyaaah31+xyzzyaaak31
xyzzyaaag31=0.d0
do xyzzyaaaa31=1,ngrid(1)
xyzzyaaag31=xyzzyaaag31+xyzzyaaaj31
xyzzyaaam31(1:2)=xyzzyaaag31*a1(1:2)+xyzzyaaah31*a2(1:2)
xyzzyaaam31(3)=0.d0
call minimum_image(3,1,xyzzyaaam31)
xyzzyaaad31=sum(xyzzyaaam31(1:2)**2)
if(xyzzyaaad31<xyzzyaaan31)then
fg_0=fg_0+xyzzyaaae31*xyzzyaaad31+xyzzyaaaf31
else
fg_0=fg_0+1.d0/sqrt(xyzzyaaad31)
endif
enddo
enddo
fg_0=fg_0*area/real(ngrid(1)*ngrid(2),dp)+3.d0*pi/4.d0*wigner_seitz_ra&
&dius
endif
end subroutine xyzzyaafx1
subroutine xyzzyaafy1
use slaarnabg, only : num_g,s_lattice
implicit none
integer xyzzyaaaa32,xyzzyaaab32,xyzzyaaac32,xyzzyaaad32,xyzzyaaae32,xy&
&zzyaaaf32,xyzzyaaag32(xyzzyaaac1),xyzzyaaah32,xyzzyaaai32,xyzzyaaaj32&
&,xyzzyaaak32,xyzzyaaal32,xyzzyaaam32,xyzzyaaan32
real(dp) xyzzyaaao32,xyzzyaaap32,xyzzyaaaq32,xyzzyaaar32,xyzzyaaas32,x&
&yzzyaaat32,xyzzyaaau32,xyzzyaaav32,xyzzyaaaw32,xyzzyaaax32,xyzzyaaay3&
&2,xyzzyaaaz32,xyzzyaaba32,xyzzyaabb32,xyzzyaabc32(xyzzyaaac1),xyzzyaa&
&bd32(xyzzyaaac1),xyzzyaabe32(xyzzyaaac1),xyzzyaabf32,xyzzyaabg32,y,xy&
&zzyaabh32,xyzzyaabi32,xyzzyaabj32,xyzzyaabk32,xyzzyaabl32,xyzzyaabm32&
&,xyzzyaabn32,xyzzyaabo32,xyzzyaabp32,xyzzyaabq32,xyzzyaabr32
real(dp),allocatable :: xyzzyaabs32(:,:),xyzzyaabt32(:,:),xyzzyaabu32(&
&:,:),xyzzyaabv32(:),xyzzyaabw32(:,:,:),xyzzyaabx32(:,:,:),xyzzyaaby32&
&(:,:)
logical,save :: xyzzyaabz32=.true.
if(xyzzyaabz32)then
xyzzyaaak32=maxval(heg_slatt)
allocate(xyzzyaaak1(xyzzyaaak32),stat=xyzzyaaaf32)
call check_alloc(xyzzyaaaf32,'PRE_SCREENING_WC','num_ps_cells')
xyzzyaaak1=0
endif
if(isperiodic.and.periodicity/=dimensionality)call errstop_master('PRE&
&_SCREENING_WC','At present Wigner-crystal calculations can only be ca&
&rried out if periodicity=dimensionality or for non-periodic systems.'&
&)
if(isperiodic.and.(type_wc==1.or.type_wc==3))then
if(xyzzyaabz32)then
allocate(xyzzyaaal1(num_g,xyzzyaaak32),site_pos_in_cell(3,xyzzyaaac1,n&
&um_g,xyzzyaaak32),xyzzyaaaj1(xyzzyaaac1,num_g,xyzzyaaak32),stat=xyzzy&
&aaaf32)
call check_alloc(xyzzyaaaf32,'PRE_SCREENING_WC','num_sites_in_cell,...&
&')
endif
xyzzyaaal1=0
site_pos_in_cell=0.d0
xyzzyaaaj1=0
if(periodicity==1)then
xyzzyaabr32=0.5d0*a1(1)
elseif(periodicity==2)then
allocate(xyzzyaabs32(4,2),xyzzyaabu32(4,2),xyzzyaabv32(4),stat=xyzzyaa&
&af32)
call check_alloc(xyzzyaaaf32,'PRE_SCREENING_WC','point,...')
xyzzyaabs32(1,:)=0.d0
xyzzyaabs32(2,:)=a1(1:2)
xyzzyaabs32(3,:)=a1(1:2)+a2(1:2)
xyzzyaabs32(4,:)=a2(1:2)
do xyzzyaaaa32=1,4
xyzzyaaah32=xyzzyaaaa32+1
if(xyzzyaaaa32==4)xyzzyaaah32=1
xyzzyaabg32=xyzzyaabs32(xyzzyaaah32,1)-xyzzyaabs32(xyzzyaaaa32,1)
y=xyzzyaabs32(xyzzyaaah32,2)-xyzzyaabs32(xyzzyaaaa32,2)
xyzzyaabf32=sqrt(xyzzyaabg32*xyzzyaabg32+y*y)
xyzzyaabv32(xyzzyaaaa32)=xyzzyaabf32
xyzzyaabu32(xyzzyaaaa32,1)=xyzzyaabg32/xyzzyaabf32
xyzzyaabu32(xyzzyaaaa32,2)=y/xyzzyaabf32
enddo
elseif(periodicity==3)then
allocate(xyzzyaabt32(6,3),xyzzyaabw32(6,4,3),xyzzyaabx32(6,4,3),xyzzya&
&aby32(6,4),stat=xyzzyaaaf32)
call check_alloc(xyzzyaaaf32,'PRE_SCREENING_WC','normal,...')
xyzzyaabw32(1,1,:)=0.d0
xyzzyaabw32(1,2,:)=a1
xyzzyaabw32(1,3,:)=a1+a3
xyzzyaabw32(1,4,:)=a3
xyzzyaabw32(2,1,:)=0.d0
xyzzyaabw32(2,2,:)=a1
xyzzyaabw32(2,3,:)=a1+a2
xyzzyaabw32(2,4,:)=a2
xyzzyaabw32(3,1,:)=0.d0
xyzzyaabw32(3,2,:)=a2
xyzzyaabw32(3,3,:)=a2+a3
xyzzyaabw32(3,4,:)=a3
xyzzyaabw32(4,1,:)=a1
xyzzyaabw32(4,2,:)=a1+a2
xyzzyaabw32(4,3,:)=a1+a2+a3
xyzzyaabw32(4,4,:)=a1+a3
xyzzyaabw32(5,1,:)=a2
xyzzyaabw32(5,2,:)=a1+a2
xyzzyaabw32(5,3,:)=a1+a2+a3
xyzzyaabw32(5,4,:)=a2+a3
xyzzyaabw32(6,1,:)=a3
xyzzyaabw32(6,2,:)=a1+a3
xyzzyaabw32(6,3,:)=a1+a2+a3
xyzzyaabw32(6,4,:)=a2+a3
do xyzzyaaaa32=1,6
do xyzzyaaah32=1,4
xyzzyaaai32=xyzzyaaah32+1
if(xyzzyaaah32==4)xyzzyaaai32=1
xyzzyaabg32=xyzzyaabw32(xyzzyaaaa32,xyzzyaaai32,1)-xyzzyaabw32(xyzzyaa&
&aa32,xyzzyaaah32,1)
y=xyzzyaabw32(xyzzyaaaa32,xyzzyaaai32,2)-xyzzyaabw32(xyzzyaaaa32,xyzzy&
&aaah32,2)
xyzzyaabh32=xyzzyaabw32(xyzzyaaaa32,xyzzyaaai32,3)-xyzzyaabw32(xyzzyaa&
&aa32,xyzzyaaah32,3)
xyzzyaabf32=sqrt(xyzzyaabg32*xyzzyaabg32+y*y+xyzzyaabh32*xyzzyaabh32)
xyzzyaaby32(xyzzyaaaa32,xyzzyaaah32)=xyzzyaabf32
xyzzyaabx32(xyzzyaaaa32,xyzzyaaah32,1)=xyzzyaabg32/xyzzyaabf32
xyzzyaabx32(xyzzyaaaa32,xyzzyaaah32,2)=y/xyzzyaabf32
xyzzyaabx32(xyzzyaaaa32,xyzzyaaah32,3)=xyzzyaabh32/xyzzyaabf32
enddo
enddo
xyzzyaabt32(1,1)=a1(2)*a3(3)-a1(3)*a3(2)
xyzzyaabt32(1,2)=-a1(1)*a3(3)+a1(3)*a3(1)
xyzzyaabt32(1,3)=a1(1)*a3(2)-a1(2)*a3(1)
xyzzyaabf32=sqrt(sum(xyzzyaabt32(1,:)**2))
xyzzyaabt32(1,:)=xyzzyaabt32(1,:)/xyzzyaabf32
xyzzyaabt32(5,:)=xyzzyaabt32(1,:)
xyzzyaabt32(2,1)=a1(2)*a2(3)-a1(3)*a2(2)
xyzzyaabt32(2,2)=-a1(1)*a2(3)+a1(3)*a2(1)
xyzzyaabt32(2,3)=a1(1)*a2(2)-a1(2)*a2(1)
xyzzyaabf32=sqrt(sum(xyzzyaabt32(2,:)**2))
xyzzyaabt32(2,:)=xyzzyaabt32(2,:)/xyzzyaabf32
xyzzyaabt32(6,:)=xyzzyaabt32(2,:)
xyzzyaabt32(3,1)=a2(2)*a3(3)-a2(3)*a3(2)
xyzzyaabt32(3,2)=-a2(1)*a3(3)+a2(3)*a3(1)
xyzzyaabt32(3,3)=a2(1)*a3(2)-a2(2)*a3(1)
xyzzyaabf32=sqrt(sum(xyzzyaabt32(3,:)**2))
xyzzyaabt32(3,:)=xyzzyaabt32(3,:)/xyzzyaabf32
xyzzyaabt32(4,:)=xyzzyaabt32(3,:)
endif
do xyzzyaaal32=1,xyzzyaaak32
xyzzyaaak1(xyzzyaaal32)=0
do xyzzyaaab32=1,num_g
xyzzyaaar32=s_lattice(1,xyzzyaaab32)
xyzzyaaas32=s_lattice(2,xyzzyaaab32)
xyzzyaaat32=s_lattice(3,xyzzyaaab32)
xyzzyaaac32=0
do xyzzyaaae32=1,heg_nosites(xyzzyaaal32)
xyzzyaaao32=xyzzyaaay1(1,xyzzyaaae32,xyzzyaaal32)
xyzzyaaap32=xyzzyaaay1(2,xyzzyaaae32,xyzzyaaal32)
xyzzyaaaq32=xyzzyaaay1(3,xyzzyaaae32,xyzzyaaal32)
xyzzyaaao32=xyzzyaaao32+xyzzyaaar32
xyzzyaaap32=xyzzyaaap32+xyzzyaaas32
xyzzyaaaq32=xyzzyaaaq32+xyzzyaaat32
xyzzyaaau32=1.d100
if(periodicity==1)then
if(xyzzyaaao32>xyzzyaabr32)then
xyzzyaaau32=xyzzyaaao32-xyzzyaabr32
elseif(xyzzyaaao32<-xyzzyaabr32)then
xyzzyaaau32=-xyzzyaabr32-xyzzyaaao32
else
xyzzyaaau32=0.d0
endif
elseif(periodicity==2)then
do xyzzyaaaa32=1,4
xyzzyaaah32=xyzzyaaaa32+1
if(xyzzyaaaa32==4)xyzzyaaah32=1
xyzzyaabg32=xyzzyaaao32-xyzzyaabs32(xyzzyaaaa32,1)
y=xyzzyaaap32-xyzzyaabs32(xyzzyaaaa32,2)
xyzzyaabf32=abs(xyzzyaabg32*xyzzyaabu32(xyzzyaaaa32,2)-y*xyzzyaabu32(x&
&yzzyaaaa32,1))
xyzzyaabi32=abs(xyzzyaabg32*xyzzyaabu32(xyzzyaaaa32,1)+y*xyzzyaabu32(x&
&yzzyaaaa32,2))
xyzzyaabg32=xyzzyaaao32-xyzzyaabs32(xyzzyaaah32,1)
y=xyzzyaaap32-xyzzyaabs32(xyzzyaaah32,2)
xyzzyaabj32=abs(xyzzyaabg32*xyzzyaabu32(xyzzyaaaa32,1)+y*xyzzyaabu32(x&
&yzzyaaaa32,2))
if(max(xyzzyaabi32,xyzzyaabj32)>xyzzyaabv32(xyzzyaaaa32))then
xyzzyaabm32=min(xyzzyaabi32,xyzzyaabj32)
xyzzyaabk32=sqrt(xyzzyaabf32*xyzzyaabf32+xyzzyaabm32*xyzzyaabm32)
else
xyzzyaabk32=xyzzyaabf32
endif
if(xyzzyaabk32<xyzzyaaau32)xyzzyaaau32=xyzzyaabk32
enddo
else
do xyzzyaaaa32=1,6
xyzzyaabg32=xyzzyaaao32-xyzzyaabw32(xyzzyaaaa32,1,1)
y=xyzzyaaap32-xyzzyaabw32(xyzzyaaaa32,1,2)
xyzzyaabh32=xyzzyaaaq32-xyzzyaabw32(xyzzyaaaa32,1,3)
xyzzyaabf32=abs(xyzzyaabg32*xyzzyaabt32(xyzzyaaaa32,1)+y*xyzzyaabt32(x&
&yzzyaaaa32,2)+xyzzyaabh32*xyzzyaabt32(xyzzyaaaa32,3))
xyzzyaabg32=xyzzyaaao32+xyzzyaabf32*xyzzyaabt32(xyzzyaaaa32,1)
xyzzyaaaw32=xyzzyaabg32-xyzzyaabw32(xyzzyaaaa32,1,1)
y=xyzzyaaap32+xyzzyaabf32*xyzzyaabt32(xyzzyaaaa32,2)
xyzzyaaax32=y-xyzzyaabw32(xyzzyaaaa32,1,2)
xyzzyaabh32=xyzzyaaaq32+xyzzyaabf32*xyzzyaabt32(xyzzyaaaa32,3)
xyzzyaaay32=xyzzyaabh32-xyzzyaabw32(xyzzyaaaa32,1,3)
xyzzyaabk32=abs(xyzzyaaaw32*xyzzyaabt32(xyzzyaaaa32,1)+xyzzyaaax32*xyz&
&zyaabt32(xyzzyaaaa32,2)+xyzzyaaay32*xyzzyaabt32(xyzzyaaaa32,3))
if(xyzzyaabk32>xyzzyaabf32)then
xyzzyaabg32=xyzzyaaao32-xyzzyaabf32*xyzzyaabt32(xyzzyaaaa32,1)
xyzzyaaaw32=xyzzyaabg32-xyzzyaabw32(xyzzyaaaa32,1,1)
y=xyzzyaaap32-xyzzyaabf32*xyzzyaabt32(xyzzyaaaa32,2)
xyzzyaaax32=y-xyzzyaabw32(xyzzyaaaa32,1,2)
xyzzyaabh32=xyzzyaaaq32-xyzzyaabf32*xyzzyaabt32(xyzzyaaaa32,3)
xyzzyaaay32=xyzzyaabh32-xyzzyaabw32(xyzzyaaaa32,1,3)
endif
xyzzyaaaz32=xyzzyaaaw32*ainv(1,1)+xyzzyaaax32*ainv(2,1)+xyzzyaaay32*ai&
&nv(3,1)
xyzzyaaba32=xyzzyaaaw32*ainv(1,2)+xyzzyaaax32*ainv(2,2)+xyzzyaaay32*ai&
&nv(3,2)
xyzzyaabb32=xyzzyaaaw32*ainv(1,3)+xyzzyaaax32*ainv(2,3)+xyzzyaaay32*ai&
&nv(3,3)
xyzzyaabq32=1.d-8
if(xyzzyaaaz32>=-xyzzyaabq32.and.xyzzyaaaz32<=(1.d0+xyzzyaabq32).and.x&
&yzzyaaba32>=-xyzzyaabq32.and.xyzzyaaba32<=(1.d0+xyzzyaabq32).and.xyzz&
&yaabb32>=-xyzzyaabq32.and.xyzzyaabb32<=(1.d0+xyzzyaabq32))then
xyzzyaabk32=xyzzyaabf32
if(xyzzyaabk32<xyzzyaaau32)xyzzyaaau32=xyzzyaabk32
else
do xyzzyaaah32=1,4
xyzzyaaai32=xyzzyaaah32+1
if(xyzzyaaah32==4)xyzzyaaai32=1
xyzzyaaaw32=xyzzyaabg32-xyzzyaabw32(xyzzyaaaa32,xyzzyaaah32,1)
xyzzyaaax32=y-xyzzyaabw32(xyzzyaaaa32,xyzzyaaah32,2)
xyzzyaaay32=xyzzyaabh32-xyzzyaabw32(xyzzyaaaa32,xyzzyaaah32,3)
xyzzyaabn32=xyzzyaaax32*xyzzyaabx32(xyzzyaaaa32,xyzzyaaah32,3)-xyzzyaa&
&ay32*xyzzyaabx32(xyzzyaaaa32,xyzzyaaah32,2)
xyzzyaabo32=-xyzzyaaaw32*xyzzyaabx32(xyzzyaaaa32,xyzzyaaah32,3)+xyzzya&
&aay32*xyzzyaabx32(xyzzyaaaa32,xyzzyaaah32,1)
xyzzyaabp32=xyzzyaaaw32*xyzzyaabx32(xyzzyaaaa32,xyzzyaaah32,2)-xyzzyaa&
&ax32*xyzzyaabx32(xyzzyaaaa32,xyzzyaaah32,1)
xyzzyaabl32=sqrt(xyzzyaabn32*xyzzyaabn32+xyzzyaabo32*xyzzyaabo32+xyzzy&
&aabp32*xyzzyaabp32)
xyzzyaabi32=abs(xyzzyaaaw32*xyzzyaabx32(xyzzyaaaa32,xyzzyaaah32,1)+xyz&
&zyaaax32*xyzzyaabx32(xyzzyaaaa32,xyzzyaaah32,2)+xyzzyaaay32*xyzzyaabx&
&32(xyzzyaaaa32,xyzzyaaah32,3))
xyzzyaaaw32=xyzzyaabg32-xyzzyaabw32(xyzzyaaaa32,xyzzyaaai32,1)
xyzzyaaax32=y-xyzzyaabw32(xyzzyaaaa32,xyzzyaaai32,2)
xyzzyaaay32=xyzzyaabh32-xyzzyaabw32(xyzzyaaaa32,xyzzyaaai32,3)
xyzzyaabj32=abs(xyzzyaaaw32*xyzzyaabx32(xyzzyaaaa32,xyzzyaaah32,1)+xyz&
&zyaaax32*xyzzyaabx32(xyzzyaaaa32,xyzzyaaah32,2)+xyzzyaaay32*xyzzyaabx&
&32(xyzzyaaaa32,xyzzyaaah32,3))
if(max(xyzzyaabi32,xyzzyaabj32)>xyzzyaaby32(xyzzyaaaa32,xyzzyaaah32))t&
&hen
xyzzyaabm32=min(xyzzyaabi32,xyzzyaabj32)
xyzzyaabk32=sqrt(xyzzyaabf32*xyzzyaabf32+xyzzyaabl32*xyzzyaabl32+xyzzy&
&aabm32*xyzzyaabm32)
else
xyzzyaabk32=sqrt(xyzzyaabf32*xyzzyaabf32+xyzzyaabl32*xyzzyaabl32)
endif
if(xyzzyaabk32<xyzzyaaau32)xyzzyaaau32=xyzzyaabk32
enddo
endif
enddo
endif
do xyzzyaaad32=1,nspin
xyzzyaaam32=0
do xyzzyaaan32=1,ndet
if(heg_slatt(xyzzyaaad32,xyzzyaaan32)==xyzzyaaal32)then
xyzzyaaam32=xyzzyaaal32
exit
endif
enddo
if(xyzzyaaam32==0)cycle
xyzzyaaaj32=which_ssingle(xyzzyaaad32,spin_dep_wc)
xyzzyaaav32=-wc_gauss_exp(xyzzyaaaj32)*xyzzyaaau32*xyzzyaaau32
if(xyzzyaaav32>xyzzyaaav1.or.xyzzyaaab32==1)then
xyzzyaaac32=xyzzyaaac32+1
xyzzyaabc32(xyzzyaaac32)=xyzzyaaao32
xyzzyaabd32(xyzzyaaac32)=xyzzyaaap32
xyzzyaabe32(xyzzyaaac32)=xyzzyaaaq32
xyzzyaaag32(xyzzyaaac32)=xyzzyaaae32
exit
endif
enddo
enddo
if(xyzzyaaac32>0)then
xyzzyaaak1(xyzzyaaal32)=xyzzyaaak1(xyzzyaaal32)+1
xyzzyaaal1(xyzzyaaak1(xyzzyaaal32),xyzzyaaal32)=xyzzyaaac32
do xyzzyaaae32=1,xyzzyaaac32
site_pos_in_cell(1:3,xyzzyaaae32,xyzzyaaak1(xyzzyaaal32),xyzzyaaal32)=&
&(/xyzzyaabc32(xyzzyaaae32),xyzzyaabd32(xyzzyaaae32),xyzzyaabe32(xyzzy&
&aaae32)/)
xyzzyaaaj1(xyzzyaaae32,xyzzyaaak1(xyzzyaaal32),xyzzyaaal32)=xyzzyaaag3&
&2(xyzzyaaae32)+(xyzzyaaal32-1)*maxval(heg_nele)
enddo
endif
enddo
if(xyzzyaabz32.and.any(heg_nele(:)/=0).and.xyzzyaaak1(xyzzyaaal32)==0)&
&call errstop_master('PRE_SCREENING_WC','Orbital screening problem: no&
& cells contain significant centres')
enddo
elseif(xyzzyaabz32)then
allocate(xyzzyaaal1(1,xyzzyaaak32),site_pos_in_cell(3,xyzzyaaac1,1,xyz&
&zyaaak32),xyzzyaaaj1(xyzzyaaac1,1,xyzzyaaak32),stat=xyzzyaaaf32)
call check_alloc(xyzzyaaaf32,'PRE_SCREENING_WC','num_sites_in_cell (2)&
&')
xyzzyaaak1(:)=1
xyzzyaaal1=0
do xyzzyaaal32=1,xyzzyaaak32
xyzzyaaal1(1,xyzzyaaal32)=heg_nosites(xyzzyaaal32)
do xyzzyaaae32=1,heg_nosites(xyzzyaaal32)
site_pos_in_cell(1:3,xyzzyaaae32,1,xyzzyaaal32)=xyzzyaaay1(1:3,xyzzyaa&
&ae32,xyzzyaaal32)
xyzzyaaaj1(xyzzyaaae32,1,xyzzyaaal32)=xyzzyaaae32+(xyzzyaaal32-1)*maxv&
&al(heg_nele)
enddo
enddo
endif
if(allocated(xyzzyaabs32))deallocate(xyzzyaabs32)
if(allocated(xyzzyaabt32))deallocate(xyzzyaabt32)
if(allocated(xyzzyaabu32))deallocate(xyzzyaabu32)
if(allocated(xyzzyaabv32))deallocate(xyzzyaabv32)
if(allocated(xyzzyaabw32))deallocate(xyzzyaabw32)
if(allocated(xyzzyaabx32))deallocate(xyzzyaabx32)
if(allocated(xyzzyaaby32))deallocate(xyzzyaaby32)
xyzzyaabz32=.false.
end subroutine xyzzyaafy1
subroutine xyzzyaafz1
use slaarnabg, only : num_g,s_lattice
implicit none
integer xyzzyaaaa33,xyzzyaaab33,xyzzyaaac33,xyzzyaaad33,xyzzyaaae33,xy&
&zzyaaaf33
real(dp) xyzzyaaag33,xyzzyaaah33,xyzzyaaai33,xyzzyaaaj33,xyzzyaaak33,y&
&,xyzzyaaal33,xyzzyaaam33(3),xyzzyaaan33,xyzzyaaao33,xyzzyaaap33,xyzzy&
&aaaq33,xyzzyaaar33,xyzzyaaas33,xyzzyaaat33,xyzzyaaau33,xyzzyaaav33,xy&
&zzyaaaw33,xyzzyaaax33,xyzzyaaay33,xyzzyaaaz33,xyzzyaaba33,xyzzyaabb33&
&,xyzzyaabc33,xyzzyaabd33
real(dp),allocatable :: xyzzyaabe33(:,:),xyzzyaabf33(:,:),xyzzyaabg33(&
&:,:),xyzzyaabh33(:),xyzzyaabi33(:,:,:),xyzzyaabj33(:,:,:),xyzzyaabk33&
&(:,:),xyzzyaabl33(:,:),xyzzyaabm33(:,:,:)
logical xyzzyaabn33
logical,save :: xyzzyaabo33=.true.
select case(xyzzyaabu1)
case(0)
xyzzyaaaf33=1
end select
if(periodicity==1)call errstop_master('PRE_SCREENING_PAIRING','Pairing&
& orbitals cannot (yet) be used for 1D HEGs.')
if(isperiodic.and.periodicity/=dimensionality)call errstop_master('PRE&
&_SCREENING_PAIRING','This subroutine assumes that periodicity=dimensi&
&onality.')
if(isperiodic.and.(xyzzyaadj1.or.xyzzyaadk1))then
if(xyzzyaabo33)allocate(xyzzyaaax1(3,num_g))
xyzzyaaad1=0
if(periodicity==2)then
allocate(xyzzyaabe33(4,2),xyzzyaabg33(4,2),xyzzyaabh33(4),xyzzyaabl33(&
&4,2),stat=xyzzyaaab33)
call check_alloc(xyzzyaaab33,'PRE_SCREENING_PAIRING','1')
xyzzyaabe33(1,:)=0.d0
xyzzyaabe33(2,:)=a1(1:2)
xyzzyaabe33(3,:)=a1(1:2)+a2(1:2)
xyzzyaabe33(4,:)=a2(1:2)
do xyzzyaaac33=1,4
xyzzyaaad33=xyzzyaaac33+1
if(xyzzyaaac33==4)xyzzyaaad33=1
xyzzyaaak33=xyzzyaabe33(xyzzyaaad33,1)-xyzzyaabe33(xyzzyaaac33,1)
y=xyzzyaabe33(xyzzyaaad33,2)-xyzzyaabe33(xyzzyaaac33,2)
xyzzyaaan33=sqrt(xyzzyaaak33*xyzzyaaak33+y*y)
xyzzyaabh33(xyzzyaaac33)=xyzzyaaan33
xyzzyaabg33(xyzzyaaac33,1)=xyzzyaaak33/xyzzyaaan33
xyzzyaabg33(xyzzyaaac33,2)=y/xyzzyaaan33
enddo
xyzzyaaam33(1)=0.5d0*(a1(1)+a2(1))
xyzzyaaam33(2)=0.5d0*(a1(2)+a2(2))
xyzzyaaam33(3)=0.d0
else
allocate(xyzzyaabf33(6,3),xyzzyaabi33(6,4,3),xyzzyaabj33(6,4,3),xyzzya&
&abk33(6,4),xyzzyaabm33(6,4,3),stat=xyzzyaaab33)
call check_alloc(xyzzyaaab33,'PRE_SCREENING_PAIRING','2')
xyzzyaabi33(1,1,:)=0.d0
xyzzyaabi33(1,2,:)=a1
xyzzyaabi33(1,3,:)=a1+a3
xyzzyaabi33(1,4,:)=a3
xyzzyaabi33(2,1,:)=0.d0
xyzzyaabi33(2,2,:)=a1
xyzzyaabi33(2,3,:)=a1+a2
xyzzyaabi33(2,4,:)=a2
xyzzyaabi33(3,1,:)=0.d0
xyzzyaabi33(3,2,:)=a2
xyzzyaabi33(3,3,:)=a2+a3
xyzzyaabi33(3,4,:)=a3
xyzzyaabi33(4,1,:)=a1
xyzzyaabi33(4,2,:)=a1+a2
xyzzyaabi33(4,3,:)=a1+a2+a3
xyzzyaabi33(4,4,:)=a1+a3
xyzzyaabi33(5,1,:)=a2
xyzzyaabi33(5,2,:)=a1+a2
xyzzyaabi33(5,3,:)=a1+a2+a3
xyzzyaabi33(5,4,:)=a2+a3
xyzzyaabi33(6,1,:)=a3
xyzzyaabi33(6,2,:)=a1+a3
xyzzyaabi33(6,3,:)=a1+a2+a3
xyzzyaabi33(6,4,:)=a2+a3
do xyzzyaaac33=1,6
do xyzzyaaad33=1,4
xyzzyaaae33=xyzzyaaad33+1
if(xyzzyaaad33==4)xyzzyaaae33=1
xyzzyaaak33=xyzzyaabi33(xyzzyaaac33,xyzzyaaae33,1)-xyzzyaabi33(xyzzyaa&
&ac33,xyzzyaaad33,1)
y=xyzzyaabi33(xyzzyaaac33,xyzzyaaae33,2)-xyzzyaabi33(xyzzyaaac33,xyzzy&
&aaad33,2)
xyzzyaaal33=xyzzyaabi33(xyzzyaaac33,xyzzyaaae33,3)-xyzzyaabi33(xyzzyaa&
&ac33,xyzzyaaad33,3)
xyzzyaaan33=sqrt(xyzzyaaak33*xyzzyaaak33+y*y+xyzzyaaal33*xyzzyaaal33)
xyzzyaabk33(xyzzyaaac33,xyzzyaaad33)=xyzzyaaan33
xyzzyaabj33(xyzzyaaac33,xyzzyaaad33,1)=xyzzyaaak33/xyzzyaaan33
xyzzyaabj33(xyzzyaaac33,xyzzyaaad33,2)=y/xyzzyaaan33
xyzzyaabj33(xyzzyaaac33,xyzzyaaad33,3)=xyzzyaaal33/xyzzyaaan33
enddo
enddo
xyzzyaabf33(1,1)=a1(2)*a3(3)-a1(3)*a3(2)
xyzzyaabf33(1,2)=-a1(1)*a3(3)+a1(3)*a3(1)
xyzzyaabf33(1,3)=a1(1)*a3(2)-a1(2)*a3(1)
xyzzyaaan33=sqrt(sum(xyzzyaabf33(1,:)**2))
xyzzyaabf33(1,:)=xyzzyaabf33(1,:)/xyzzyaaan33
xyzzyaabf33(5,:)=xyzzyaabf33(1,:)
xyzzyaabf33(2,1)=a1(2)*a2(3)-a1(3)*a2(2)
xyzzyaabf33(2,2)=-a1(1)*a2(3)+a1(3)*a2(1)
xyzzyaabf33(2,3)=a1(1)*a2(2)-a1(2)*a2(1)
xyzzyaaan33=sqrt(sum(xyzzyaabf33(2,:)**2))
xyzzyaabf33(2,:)=xyzzyaabf33(2,:)/xyzzyaaan33
xyzzyaabf33(6,:)=xyzzyaabf33(2,:)
xyzzyaabf33(3,1)=a2(2)*a3(3)-a2(3)*a3(2)
xyzzyaabf33(3,2)=-a2(1)*a3(3)+a2(3)*a3(1)
xyzzyaabf33(3,3)=a2(1)*a3(2)-a2(2)*a3(1)
xyzzyaaan33=sqrt(sum(xyzzyaabf33(3,:)**2))
xyzzyaabf33(3,:)=xyzzyaabf33(3,:)/xyzzyaaan33
xyzzyaabf33(4,:)=xyzzyaabf33(3,:)
xyzzyaaam33(:)=0.5d0*(a1(:)+a2(:)+a3(:))
endif
do xyzzyaaaa33=1,num_g
xyzzyaaag33=s_lattice(1,xyzzyaaaa33)
xyzzyaaah33=s_lattice(2,xyzzyaaaa33)
xyzzyaaai33=s_lattice(3,xyzzyaaaa33)
if(periodicity==2)then
xyzzyaabl33(:,1)=xyzzyaabe33(:,1)+xyzzyaaag33
xyzzyaabl33(:,2)=xyzzyaabe33(:,2)+xyzzyaaah33
else
xyzzyaabm33(:,:,1)=xyzzyaabi33(:,:,1)+xyzzyaaag33
xyzzyaabm33(:,:,2)=xyzzyaabi33(:,:,2)+xyzzyaaah33
xyzzyaabm33(:,:,3)=xyzzyaabi33(:,:,3)+xyzzyaaai33
endif
xyzzyaaao33=1.d100
if(periodicity==2)then
do xyzzyaaac33=1,4
xyzzyaaad33=xyzzyaaac33+1
if(xyzzyaaac33==4)xyzzyaaad33=1
xyzzyaaak33=xyzzyaaam33(1)-xyzzyaabl33(xyzzyaaac33,1)
y=xyzzyaaam33(2)-xyzzyaabl33(xyzzyaaac33,2)
xyzzyaaan33=abs(xyzzyaaak33*xyzzyaabg33(xyzzyaaac33,2)-y*xyzzyaabg33(x&
&yzzyaaac33,1))
xyzzyaaba33=abs(xyzzyaaak33*xyzzyaabg33(xyzzyaaac33,1)+y*xyzzyaabg33(x&
&yzzyaaac33,2))
xyzzyaaak33=xyzzyaaam33(1)-xyzzyaabl33(xyzzyaaad33,1)
y=xyzzyaaam33(2)-xyzzyaabl33(xyzzyaaad33,2)
xyzzyaabb33=abs(xyzzyaaak33*xyzzyaabg33(xyzzyaaac33,1)+y*xyzzyaabg33(x&
&yzzyaaac33,2))
if(max(xyzzyaaba33,xyzzyaabb33)>xyzzyaabh33(xyzzyaaac33))then
xyzzyaabc33=min(xyzzyaaba33,xyzzyaabb33)
xyzzyaaav33=sqrt(xyzzyaaan33*xyzzyaaan33+xyzzyaabc33*xyzzyaabc33)
else
xyzzyaaav33=xyzzyaaan33
endif
if(xyzzyaaav33<xyzzyaaao33)xyzzyaaao33=xyzzyaaav33
enddo
else
do xyzzyaaac33=1,6
xyzzyaaak33=xyzzyaaam33(1)-xyzzyaabm33(xyzzyaaac33,1,1)
y=xyzzyaaam33(2)-xyzzyaabm33(xyzzyaaac33,1,2)
xyzzyaaal33=xyzzyaaam33(3)-xyzzyaabm33(xyzzyaaac33,1,3)
xyzzyaaan33=abs(xyzzyaaak33*xyzzyaabf33(xyzzyaaac33,1)+y*xyzzyaabf33(x&
&yzzyaaac33,2)+xyzzyaaal33*xyzzyaabf33(xyzzyaaac33,3))
xyzzyaaak33=xyzzyaaam33(1)+xyzzyaaan33*xyzzyaabf33(xyzzyaaac33,1)
xyzzyaaap33=xyzzyaaak33-xyzzyaabm33(xyzzyaaac33,1,1)
y=xyzzyaaam33(2)+xyzzyaaan33*xyzzyaabf33(xyzzyaaac33,2)
xyzzyaaaq33=y-xyzzyaabm33(xyzzyaaac33,1,2)
xyzzyaaal33=xyzzyaaam33(3)+xyzzyaaan33*xyzzyaabf33(xyzzyaaac33,3)
xyzzyaaar33=xyzzyaaal33-xyzzyaabm33(xyzzyaaac33,1,3)
xyzzyaaav33=abs(xyzzyaaap33*xyzzyaabf33(xyzzyaaac33,1)+xyzzyaaaq33*xyz&
&zyaabf33(xyzzyaaac33,2)+xyzzyaaar33*xyzzyaabf33(xyzzyaaac33,3))
if(xyzzyaaav33>xyzzyaaan33)then
xyzzyaaak33=xyzzyaaam33(1)-xyzzyaaan33*xyzzyaabf33(xyzzyaaac33,1)
xyzzyaaap33=xyzzyaaak33-xyzzyaabm33(xyzzyaaac33,1,1)
y=xyzzyaaam33(2)-xyzzyaaan33*xyzzyaabf33(xyzzyaaac33,2)
xyzzyaaaq33=y-xyzzyaabm33(xyzzyaaac33,1,2)
xyzzyaaal33=xyzzyaaam33(3)-xyzzyaaan33*xyzzyaabf33(xyzzyaaac33,3)
xyzzyaaar33=xyzzyaaal33-xyzzyaabm33(xyzzyaaac33,1,3)
endif
xyzzyaaax33=xyzzyaaap33*ainv(1,1)+xyzzyaaaq33*ainv(2,1)+xyzzyaaar33*ai&
&nv(3,1)
xyzzyaaay33=xyzzyaaap33*ainv(1,2)+xyzzyaaaq33*ainv(2,2)+xyzzyaaar33*ai&
&nv(3,2)
xyzzyaaaz33=xyzzyaaap33*ainv(1,3)+xyzzyaaaq33*ainv(2,3)+xyzzyaaar33*ai&
&nv(3,3)
xyzzyaabd33=1.d-8
if(xyzzyaaax33>=-xyzzyaabd33.and.xyzzyaaax33<=(1.d0+xyzzyaabd33).and.x&
&yzzyaaay33>=-xyzzyaabd33.and.xyzzyaaay33<=(1.d0+xyzzyaabd33).and.xyzz&
&yaaaz33>=-xyzzyaabd33.and.xyzzyaaaz33<=(1.d0+xyzzyaabd33))then
xyzzyaaav33=xyzzyaaan33
if(xyzzyaaav33<xyzzyaaao33)xyzzyaaao33=xyzzyaaav33
else
do xyzzyaaad33=1,4
xyzzyaaae33=xyzzyaaad33+1
if(xyzzyaaad33==4)xyzzyaaae33=1
xyzzyaaap33=xyzzyaaak33-xyzzyaabm33(xyzzyaaac33,xyzzyaaad33,1)
xyzzyaaaq33=y-xyzzyaabm33(xyzzyaaac33,xyzzyaaad33,2)
xyzzyaaar33=xyzzyaaal33-xyzzyaabm33(xyzzyaaac33,xyzzyaaad33,3)
xyzzyaaas33=xyzzyaaaq33*xyzzyaabj33(xyzzyaaac33,xyzzyaaad33,3)-xyzzyaa&
&ar33*xyzzyaabj33(xyzzyaaac33,xyzzyaaad33,2)
xyzzyaaat33=-xyzzyaaap33*xyzzyaabj33(xyzzyaaac33,xyzzyaaad33,3)+xyzzya&
&aar33*xyzzyaabj33(xyzzyaaac33,xyzzyaaad33,1)
xyzzyaaau33=xyzzyaaap33*xyzzyaabj33(xyzzyaaac33,xyzzyaaad33,2)-xyzzyaa&
&aq33*xyzzyaabj33(xyzzyaaac33,xyzzyaaad33,1)
xyzzyaaaw33=sqrt(xyzzyaaas33*xyzzyaaas33+xyzzyaaat33*xyzzyaaat33+xyzzy&
&aaau33*xyzzyaaau33)
xyzzyaaba33=abs(xyzzyaaap33*xyzzyaabj33(xyzzyaaac33,xyzzyaaad33,1)+xyz&
&zyaaaq33*xyzzyaabj33(xyzzyaaac33,xyzzyaaad33,2)+xyzzyaaar33*xyzzyaabj&
&33(xyzzyaaac33,xyzzyaaad33,3))
xyzzyaaap33=xyzzyaaak33-xyzzyaabm33(xyzzyaaac33,xyzzyaaae33,1)
xyzzyaaaq33=y-xyzzyaabm33(xyzzyaaac33,xyzzyaaae33,2)
xyzzyaaar33=xyzzyaaal33-xyzzyaabm33(xyzzyaaac33,xyzzyaaae33,3)
xyzzyaabb33=abs(xyzzyaaap33*xyzzyaabj33(xyzzyaaac33,xyzzyaaad33,1)+xyz&
&zyaaaq33*xyzzyaabj33(xyzzyaaac33,xyzzyaaad33,2)+xyzzyaaar33*xyzzyaabj&
&33(xyzzyaaac33,xyzzyaaad33,3))
if(max(xyzzyaaba33,xyzzyaabb33)>xyzzyaabk33(xyzzyaaac33,xyzzyaaad33))t&
&hen
xyzzyaabc33=min(xyzzyaaba33,xyzzyaabb33)
xyzzyaaav33=sqrt(xyzzyaaan33*xyzzyaaan33+xyzzyaaaw33*xyzzyaaaw33+xyzzy&
&aabc33*xyzzyaabc33)
else
xyzzyaaav33=sqrt(xyzzyaaan33*xyzzyaaan33+xyzzyaaaw33*xyzzyaaaw33)
endif
if(xyzzyaaav33<xyzzyaaao33)xyzzyaaao33=xyzzyaaav33
enddo
endif
enddo
endif
xyzzyaabn33=xyzzyaaaa33==1
if(xyzzyaadj1)then
xyzzyaaaj33=-xyzzyaaao33*xyzzyaaao33*minval(xyzzyaacv1(:,xyzzyaaaf33))
xyzzyaabn33=xyzzyaabn33.or.xyzzyaaaj33>xyzzyaaav1
endif
if(xyzzyaadk1)xyzzyaabn33=xyzzyaabn33.or.xyzzyaaao33<xyzzyaacp1
if(xyzzyaabn33)then
xyzzyaaad1=xyzzyaaad1+1
if(periodicity==2)xyzzyaaai33=0.d0
xyzzyaaax1(1:3,xyzzyaaad1)=(/xyzzyaaag33,xyzzyaaah33,xyzzyaaai33/)
endif
enddo
if(xyzzyaaad1==0.and.xyzzyaabo33)call errstop_master('PRE_SCREENING_PA&
&IRING','Orbital screening problem: no cells are significant.')
elseif(xyzzyaabo33)then
allocate(xyzzyaaax1(3,1),stat=xyzzyaaab33)
call check_alloc(xyzzyaaab33,'PRE_SCREENING_PAIRING','3')
xyzzyaaad1=1
xyzzyaaax1(1:3,1)=0.d0
endif
if(allocated(xyzzyaabe33))deallocate(xyzzyaabe33)
if(allocated(xyzzyaabf33))deallocate(xyzzyaabf33)
if(allocated(xyzzyaabg33))deallocate(xyzzyaabg33)
if(allocated(xyzzyaabh33))deallocate(xyzzyaabh33)
if(allocated(xyzzyaabi33))deallocate(xyzzyaabi33)
if(allocated(xyzzyaabj33))deallocate(xyzzyaabj33)
if(allocated(xyzzyaabk33))deallocate(xyzzyaabk33)
if(allocated(xyzzyaabl33))deallocate(xyzzyaabl33)
if(allocated(xyzzyaabm33))deallocate(xyzzyaabm33)
xyzzyaabo33=.false.
end subroutine xyzzyaafz1
subroutine xyzzyaaga1
implicit none
integer xyzzyaaaa34,xyzzyaaab34,xyzzyaaac34,xyzzyaaad34
integer,allocatable :: xyzzyaaae34(:)
real(dp) xyzzyaaaf34,xyzzyaaag34,xyzzyaaah34,xyzzyaaai34,xyzzyaaaj34
real(dp),allocatable :: xyzzyaaak34(:,:),xyzzyaaal34(:)
logical xyzzyaaam34
xyzzyaaaa34=xyzzyaabq1
allocate(xyzzyaaak34(xyzzyaaaa34,xyzzyaaaa34),xyzzyaaae34(xyzzyaaaa34)&
&,xyzzyaacy1(xyzzyaaaa34),xyzzyaaal34(xyzzyaaaa34),xyzzyaacz1(xyzzyaaa&
&a34),stat=xyzzyaaad34)
call check_alloc(xyzzyaaad34,'GAUSS_FIT_COEFF','')
do xyzzyaaab34=1,xyzzyaaaa34
xyzzyaaaf34=1.d0/(log(real(xyzzyaaaa34+1,dp))-log(real(xyzzyaaab34,dp)&
&))
xyzzyaacy1(xyzzyaaab34)=xyzzyaaaf34
xyzzyaaal34(xyzzyaaab34)=exp(0.25d0/xyzzyaaaf34**2)*erfc(0.5d0/xyzzyaa&
&af34)
enddo
do xyzzyaaab34=1,xyzzyaaaa34
xyzzyaaaf34=xyzzyaacy1(xyzzyaaab34)
xyzzyaaag34=xyzzyaaaf34**2
xyzzyaaak34(xyzzyaaab34,xyzzyaaab34)=one_over_root_two
do xyzzyaaac34=xyzzyaaab34+1,xyzzyaaaa34
xyzzyaaah34=xyzzyaacy1(xyzzyaaac34)
xyzzyaaai34=xyzzyaaah34**2
xyzzyaaaj34=1.d0/sqrt(xyzzyaaag34+xyzzyaaai34)
xyzzyaaak34(xyzzyaaab34,xyzzyaaac34)=xyzzyaaaf34*xyzzyaaaj34
xyzzyaaak34(xyzzyaaac34,xyzzyaaab34)=xyzzyaaah34*xyzzyaaaj34
enddo
enddo
call lu_decom(xyzzyaaak34,xyzzyaaae34,xyzzyaaaa34,xyzzyaaaa34,xyzzyaaa&
&m34)
if(xyzzyaaam34)call errstop_master('GAUSS_FIT_COEFF','Singular matrix &
&found in LU decomposition. Fit to exponential failed. Too many Gaussi&
&ans?')
call dcopy(xyzzyaaaa34,xyzzyaaal34,1,xyzzyaacz1,1)
call lu_solve_once(xyzzyaaak34,xyzzyaaae34,xyzzyaacz1,xyzzyaaaa34,xyzz&
&yaaaa34)
deallocate(xyzzyaaak34,xyzzyaaae34,xyzzyaaal34)
do xyzzyaaab34=2,xyzzyaaaa34
xyzzyaacz1(xyzzyaaab34)=xyzzyaacz1(xyzzyaaab34)/xyzzyaacz1(1)
enddo
xyzzyaacz1(1)=1.d0
end subroutine xyzzyaaga1
subroutine eval_finite_hf_energies(hf_ke,hf_ex)
implicit none
real(dp),intent(out) :: hf_ke
real(dp),intent(out),optional :: hf_ex
integer xyzzyaaaa35,xyzzyaaab35,xyzzyaaac35
real(dp) :: xyzzyaaad35
if(.not.xyzzyaabf1.or.xyzzyaabg1.or.xyzzyaabh1.or.xyzzyaabi1.or.xyzzya&
&adw1.or.excite_heg)call errstop('EVAL_HF_KE','MC twist averaging is o&
&nly for fluid ground states.')
hf_ke=0.d0
do xyzzyaaaa35=1,nspin
if(nele(xyzzyaaaa35)>0)hf_ke=hf_ke+sum(xyzzyaaaw1(1:xyzzyaaai1(xyzzyaa&
&aa35)))*inv_pmass(xyzzyaaaa35)
enddo
if(xyzzyaabl1)then
hf_ke=0.5d0*hf_ke/real(netot,dp)
else
hf_ke=hf_ke/real(netot,dp)
endif
if(present(hf_ex))then
hf_ex=0.d0
do xyzzyaaaa35=1,nspin
if(nele(xyzzyaaaa35)>0)then
xyzzyaaad35=0.d0
if(xyzzyaabl1)then
do xyzzyaaab35=1,xyzzyaaai1(xyzzyaaaa35)-1
do xyzzyaaac35=xyzzyaaab35+1,xyzzyaaai1(xyzzyaaaa35)
xyzzyaaad35=xyzzyaaad35+xyzzyaagb1(sum((gs_kvec(1:3,xyzzyaaab35)-gs_kv&
&ec(1:3,xyzzyaaac35))**2))
enddo
enddo
else
do xyzzyaaab35=1,xyzzyaaai1(xyzzyaaaa35)
do xyzzyaaac35=xyzzyaaab35,xyzzyaaai1(xyzzyaaaa35)
if(xyzzyaaab35/=xyzzyaaac35)then
xyzzyaaad35=xyzzyaaad35+2*xyzzyaagb1(sum((gs_kvec(1:3,xyzzyaaab35)-gs_&
&kvec(1:3,xyzzyaaac35))**2))
endif
if(xyzzyaaab35/=1)then
if(xyzzyaaab35==xyzzyaaac35)then
xyzzyaaad35=xyzzyaaad35+xyzzyaagb1(sum((gs_kvec(1:3,xyzzyaaab35)+gs_kv&
&ec(:,xyzzyaaac35))**2))
else
xyzzyaaad35=xyzzyaaad35+2*xyzzyaagb1(sum((gs_kvec(1:3,xyzzyaaab35)+gs_&
&kvec(1:3,xyzzyaaac35))**2))
endif
endif
enddo
enddo
endif
hf_ex=hf_ex+pcharge(xyzzyaaaa35)**2*xyzzyaaad35
endif
enddo
if(periodicity==3)then
hf_ex=-hf_ex/(dble(netot)*volume)
elseif(periodicity==2)then
hf_ex=-hf_ex/(dble(netot)*area)
else
hf_ex=-hf_ex/(dble(netot)*abs(a1(1)))
endif
if(xyzzyaabl1)then
do xyzzyaaaa35=1,nspin
if(nele(xyzzyaaaa35)>0)hf_ex=hf_ex+dble(xyzzyaaai1(xyzzyaaaa35))*self_&
&term*pcharge(xyzzyaaaa35)**2/dble(netot)
enddo
else
do xyzzyaaaa35=1,nspin
if(nele(xyzzyaaaa35)>0)hf_ex=hf_ex+dble(2*xyzzyaaai1(xyzzyaaaa35)-1)*s&
&elf_term*pcharge(xyzzyaaaa35)**2/dble(netot)
enddo
endif
endif
end subroutine eval_finite_hf_energies
subroutine calc_hf_energies
implicit none
integer xyzzyaaaa36,xyzzyaaab36
real(dp) xyzzyaaac36,xyzzyaaad36,xyzzyaaae36,xyzzyaaaf36,xyzzyaaag36,x&
&yzzyaaah36,xyzzyaaai36,xyzzyaaaj36
logical xyzzyaaak36
if(.not.xyzzyaabf1.or.xyzzyaabg1.or.xyzzyaabh1.or.xyzzyaabi1.or.xyzzya&
&adw1.or.excite_heg.or..not.am_master)return
call timer('HF_ENERGIES',.true.)
call wout('Squared magnitude of Fermi wave vector (au)')
call wout('===========================================')
do xyzzyaaaa36=1,nspin
if(nele(xyzzyaaaa36)>0)then
tmpr=r2s(maxval(xyzzyaaaw1(1:xyzzyaaai1(xyzzyaaaa36))),'(f20.12)')
call wout('|k_F|^2 for spin '//trim(i2s(xyzzyaaaa36))//' : '//trim(tmp&
&r))
endif
enddo
call wout()
xyzzyaaak36=.false.
select case(trim(interaction))
case('ewald','mpc','ewald_mpc','mpc_ewald')
xyzzyaaak36=.true.
end select
if(xyzzyaaak36)then
call eval_finite_hf_energies(xyzzyaaac36,xyzzyaaad36)
else
call eval_finite_hf_energies(xyzzyaaac36)
endif
call wout('Hartree-Fock energies of finite system in au per particle')
call wout('=========================================================')
if(xyzzyaaak36)call wout('(Using Ewald interaction.)')
tmpr=r2s2(xyzzyaaac36,'(f20.12)')
call wout('HF kinetic energy         : '//trim(tmpr))
if(xyzzyaaak36)then
tmpr=r2s2(xyzzyaaad36,'(f20.12)')
call wout('HF exchange energy        : '//trim(tmpr))
tmpr=r2s2(xyzzyaaac36+xyzzyaaad36,'(f20.12)')
call wout('HF energy                 : '//trim(tmpr))
endif
call wout()
if(mc_twist_av.and.xyzzyaaak36)then
xyzzyaaab36=max(100,nint(2.d7/dble(netot**2)))
call xyzzyaagd1(xyzzyaaab36,xyzzyaaac36,xyzzyaaad36,xyzzyaaai36,xyzzya&
&aaj36)
call wout('Twist-averaged Hartree-Fock energies of finite system in au&
& per particle')
call wout('===========================================================&
&=============')
call wout('(Using Ewald interaction.)')
call wout('Number of random twists   : '//trim(i2s(xyzzyaaab36)))
write(tmpr,'(a,f20.12,a,f20.12)')'HF kinetic energy         : ',xyzzya&
&aac36,' +/- ',xyzzyaaai36
call wout(tmpr)
write(tmpr,'(a,f20.12,a,f20.12)')'HF exchange energy        : ',xyzzya&
&aad36,' +/- ',xyzzyaaaj36
call wout(tmpr)
write(tmpr,'(a,f20.12,a,f20.12)')'HF energy                 : ',xyzzya&
&aac36+xyzzyaaad36,' +/- ',sqrt(xyzzyaaai36**2+xyzzyaaaj36**2)
call wout(tmpr)
call wout()
endif
if(isperiodic)then
call wout('Hartree-Fock energies of infinite system in au per particle&
&')
call wout('===========================================================&
&')
if(periodicity==3)then
xyzzyaaaf36=(6.d0*pi*pi/volume)**third
xyzzyaaac36=0.3d0
xyzzyaaad36=-0.75d0/pi
elseif(periodicity==2)then
xyzzyaaaf36=sqrt(pi/area)
xyzzyaaac36=1.d0
xyzzyaaad36=-8.d0*third/pi
else
xyzzyaaaf36=1.d0/abs(a1(1))
xyzzyaaac36=pi**2/6.d0
xyzzyaaad36=1.d0
endif
xyzzyaaae36=1.d0/dble(netot)
xyzzyaaag36=sum(nele(:)**(dble(periodicity+2)/dble(periodicity))*inv_p&
&mass(:))
xyzzyaaac36=xyzzyaaac36*xyzzyaaae36*xyzzyaaaf36**2*xyzzyaaag36
if(xyzzyaaak36)then
if(periodicity>1)then
xyzzyaaah36=sum(nele(:)**(dble(periodicity+1)/dble(periodicity))*pchar&
&ge(:)**2)
else
xyzzyaaah36=0.d0
do xyzzyaaaa36=1,nspin
if(nele(xyzzyaaaa36)>0)xyzzyaaah36=xyzzyaaah36+(dble(nele(xyzzyaaaa36)&
&)*pcharge(xyzzyaaaa36))**2 *(log(pi*dble(nele(xyzzyaaaa36))/abs(a1(1)&
&))+euler-1.5d0)
enddo
endif
xyzzyaaad36=xyzzyaaad36*xyzzyaaae36*xyzzyaaaf36*xyzzyaaah36
endif
tmpr=r2s2(xyzzyaaac36,'(f20.12)')
call wout('HF kinetic energy         : '//trim(tmpr))
if(xyzzyaaak36)then
tmpr=r2s2(xyzzyaaad36,'(f20.12)')
call wout('HF exchange energy        : '//trim(tmpr))
tmpr=r2s2(xyzzyaaac36+xyzzyaaad36,'(f20.12)')
call wout('HF energy                 : '//trim(tmpr))
endif
call wout()
endif
if(xyzzyaaau1/=0.d0)then
tmpr=r2s2(xyzzyaaau1/dble(netot),'(f20.12)')
call wout('Capacitor energy (au/part): '//trim(tmpr))
call wout()
endif
call timer('HF_ENERGIES',.false.)
end subroutine calc_hf_energies
real(dp) function xyzzyaagb1(k2)
use slaarnabt,only : exp_int
implicit none
real(dp),intent(in) :: k2
if(periodicity==3)then
xyzzyaagb1=fourpi/k2
elseif(periodicity==2)then
xyzzyaagb1=twopi/sqrt(k2)
else
if(harmwire_b<=0.d0)then
xyzzyaagb1=-log(0.25d0*k2)-two_euler
else
xyzzyaagb1=exp_int(harmwire_b**2*k2)*exp(harmwire_b**2*k2)
endif
endif
end function xyzzyaagb1
subroutine xyzzyaagc1(hfta_k_offset,hf_ke,hf_ex)
implicit none
real(dp),intent(in) :: hfta_k_offset(3)
real(dp),intent(out) :: hf_ke,hf_ex
integer xyzzyaaaa38,xyzzyaaab38,xyzzyaaac38,xyzzyaaad38
real(dp) xyzzyaaae38
do xyzzyaaac38=1,xyzzyaaef1
xyzzyaaec1(xyzzyaaac38)=sum((xyzzyaaeb1(1:3,xyzzyaaac38)+hfta_k_offset&
&(1:3))**2)
enddo
call quicksort(xyzzyaaef1-xyzzyaaee1+1,xyzzyaaec1(xyzzyaaee1),xyzzyaae&
&a1(xyzzyaaee1))
xyzzyaaea1(xyzzyaaee1:xyzzyaaef1)=xyzzyaaea1(xyzzyaaee1:xyzzyaaef1)+(x&
&yzzyaaee1-1)
hf_ke=0.d0
hf_ex=0.d0
do xyzzyaaaa38=1,nspin
xyzzyaaae38=0.d0
do xyzzyaaac38=1,xyzzyaaai1(xyzzyaaaa38)
xyzzyaaae38=xyzzyaaae38+xyzzyaaec1(xyzzyaaea1(xyzzyaaac38))
enddo
hf_ke=hf_ke+xyzzyaaae38*inv_pmass(xyzzyaaaa38)
xyzzyaaae38=0.d0
if(periodicity==3)then
do xyzzyaaab38=1,xyzzyaaai1(xyzzyaaaa38)-1
xyzzyaaad38=xyzzyaaea1(xyzzyaaab38)
do xyzzyaaac38=xyzzyaaab38+1,xyzzyaaai1(xyzzyaaaa38)
xyzzyaaae38=xyzzyaaae38+fourpi/sum((xyzzyaaeb1(1:3,xyzzyaaad38)-xyzzya&
&aeb1(1:3,xyzzyaaea1(xyzzyaaac38)))**2)
enddo
enddo
elseif(periodicity==2)then
do xyzzyaaab38=1,xyzzyaaai1(xyzzyaaaa38)-1
xyzzyaaad38=xyzzyaaea1(xyzzyaaab38)
do xyzzyaaac38=xyzzyaaab38+1,xyzzyaaai1(xyzzyaaaa38)
xyzzyaaae38=xyzzyaaae38+twopi/sqrt(sum((xyzzyaaeb1(1:2,xyzzyaaad38)-xy&
&zzyaaeb1(1:2,xyzzyaaea1(xyzzyaaac38)))**2))
enddo
enddo
else
do xyzzyaaab38=1,xyzzyaaai1(xyzzyaaaa38)-1
xyzzyaaad38=xyzzyaaea1(xyzzyaaab38)
do xyzzyaaac38=xyzzyaaab38+1,xyzzyaaai1(xyzzyaaaa38)
xyzzyaaae38=xyzzyaaae38-log(0.25d0*(xyzzyaaeb1(1,xyzzyaaad38)-xyzzyaae&
&b1(1,xyzzyaaea1(xyzzyaaac38)))**2)
enddo
enddo
xyzzyaaae38=xyzzyaaae38-euler*dble((xyzzyaaai1(xyzzyaaaa38)*(xyzzyaaai&
&1(xyzzyaaaa38)-1)))
endif
hf_ex=hf_ex+pcharge(xyzzyaaaa38)**2*xyzzyaaae38
enddo
end subroutine xyzzyaagc1
subroutine xyzzyaagd1(no_twists_hfta,hfta_ke,hfta_ex,hfta_ke_err,hfta_&
&ex_err)
use slaarnabq, only : min_image_brute_force
use slaarnacc, only : ranx
implicit none
integer,intent(in) :: no_twists_hfta
real(dp),intent(out) :: hfta_ke,hfta_ex,hfta_ke_err,hfta_ex_err
integer xyzzyaaaa39,xyzzyaaab39
real(dp) xyzzyaaac39(3),xyzzyaaad39,xyzzyaaae39,xyzzyaaaf39(3),xyzzyaa&
&ag39,xyzzyaaah39
hfta_ke=0.d0
hfta_ex=0.d0
xyzzyaaag39=0.d0
xyzzyaaah39=0.d0
do xyzzyaaaa39=1,no_twists_hfta
if(periodicity==3)then
xyzzyaaac39=ranx()*b1+ranx()*b2+ranx()*b3
elseif(periodicity==2)then
xyzzyaaac39=ranx()*b1+ranx()*b2
else
xyzzyaaac39=ranx()*b1
endif
call min_image_brute_force(periodicity,xyzzyaaac39,xyzzyaacs1,binv,xyz&
&zyaaaf39,xyzzyaact1)
call xyzzyaagc1(xyzzyaaaf39,xyzzyaaad39,xyzzyaaae39)
hfta_ke=hfta_ke+xyzzyaaad39
hfta_ex=hfta_ex+xyzzyaaae39
xyzzyaaag39=xyzzyaaag39+xyzzyaaad39**2
xyzzyaaah39=xyzzyaaah39+xyzzyaaae39**2
enddo
hfta_ke=hfta_ke/dble(no_twists_hfta)
hfta_ex=hfta_ex/dble(no_twists_hfta)
xyzzyaaag39=xyzzyaaag39/dble(no_twists_hfta)
xyzzyaaah39=xyzzyaaah39/dble(no_twists_hfta)
hfta_ke_err=sqrt(max((xyzzyaaag39-hfta_ke**2)/dble(no_twists_hfta-1),0&
&.d0))
hfta_ex_err=sqrt(max((xyzzyaaah39-hfta_ex**2)/dble(no_twists_hfta-1),0&
&.d0))
hfta_ke=0.5d0*hfta_ke/dble(netot)
hfta_ke_err=0.5d0*hfta_ke_err/dble(netot)
if(periodicity==3)then
hfta_ex=-hfta_ex/(dble(netot)*volume)
hfta_ex_err=hfta_ex_err/(dble(netot)*volume)
elseif(periodicity==2)then
hfta_ex=-hfta_ex/(dble(netot)*area)
hfta_ex_err=hfta_ex_err/(dble(netot)*area)
else
hfta_ex=-hfta_ex/(dble(netot)*abs(a1(1)))
hfta_ex_err=hfta_ex_err/(dble(netot)*abs(a1(1)))
endif
do xyzzyaaab39=1,nspin
if(nele(xyzzyaaab39)>0)hfta_ex=hfta_ex+dble(xyzzyaaai1(xyzzyaaab39))*s&
&elf_term*pcharge(xyzzyaaab39)**2 /dble(netot)
enddo
end subroutine xyzzyaagd1
subroutine get_free_orbmap(row_offset,norb,orbmap)
implicit none
integer,intent(inout) :: row_offset(nspin),norb,orbmap(nemax,nspin,nde&
&t)
integer xyzzyaaaa40,xyzzyaaab40,xyzzyaaac40
if(xyzzyaabg1)then
do xyzzyaaaa40=1,nspin
do xyzzyaaab40=1,ndet
xyzzyaaac40=max(heg_nele(xyzzyaaaa40),heg_nele(-heg_orbtype(xyzzyaaaa4&
&0,xyzzyaaab40)))
orbmap(row_offset(xyzzyaaaa40)+1:row_offset(xyzzyaaaa40)+xyzzyaaac40,x&
&yzzyaaaa40,xyzzyaaab40)=norb+xyzzyaaao1(1:xyzzyaaac40,xyzzyaaaa40,xyz&
&zyaaab40)
row_offset(xyzzyaaaa40)=row_offset(xyzzyaaaa40)+xyzzyaaac40
enddo
enddo
norb=norb+free_norb
elseif(xyzzyaabh1.or.xyzzyaabf1)then
do xyzzyaaaa40=1,nspin
orbmap(row_offset(xyzzyaaaa40)+1:row_offset(xyzzyaaaa40)+heg_nele(xyzz&
&yaaaa40),xyzzyaaaa40,:)=norb+xyzzyaaao1(1:heg_nele(xyzzyaaaa40),xyzzy&
&aaaa40,:)
row_offset(xyzzyaaaa40)=row_offset(xyzzyaaaa40)+heg_nele(xyzzyaaaa40)
enddo
norb=norb+free_norb
elseif(xyzzyaabm1.or.xyzzyaabn1.or.xyzzyaabo1)then
do xyzzyaaaa40=1,nspin
orbmap(row_offset(xyzzyaaaa40)+1,xyzzyaaaa40,:)=norb+xyzzyaaao1(1,xyzz&
&yaaaa40,:)
row_offset(xyzzyaaaa40)=row_offset(xyzzyaaaa40)+1
enddo
norb=norb+free_norb
endif
end subroutine get_free_orbmap
subroutine get_free_ndesc(ndesc_int,ndesc_dp)
implicit none
integer,intent(inout) :: ndesc_int,ndesc_dp
ndesc_int=xyzzyaaff1
ndesc_dp=xyzzyaafg1
end subroutine get_free_ndesc
subroutine get_free_orbdesc(norb,ndesc_int,ndesc_dp,orbdesc_int,orbdes&
&c_dp)
implicit none
integer,intent(in) :: norb,ndesc_int,ndesc_dp
integer,intent(inout) :: orbdesc_int(ndesc_int,norb)
real(dp),intent(inout) :: orbdesc_dp(ndesc_dp,norb)
if(xyzzyaabf1)then
orbdesc_int(1:xyzzyaaff1,1:free_norb)=xyzzyaafh1(1:xyzzyaaff1,1:free_n&
&orb)
endif
end subroutine get_free_orbdesc
end module slaarnaas

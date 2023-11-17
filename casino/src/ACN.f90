module slaarnacm_linjas
use slaarnaaf
use dsp
use slaarnaam
use parallel
use slaarnach
use slaarnacs
use casl,          only : write_casl
use slaarnaag,     only : twopi,third,ninth,rec_54
use slaarnaah,   only : write_correlation_header
use format_utils,  only : wout,i2s,byte2human,labelled_list
use slaarnabl,       only : get_linear_basis
use slaarnabt,     only : choose,ddot,daxpy,dscal,dadd,dcopy
use slaarnacc,only : ranx,get_random_state,put_random_state
use run_control,   only : errstop_master,errstop,errwarn,timer,loop_ti&
&me_estimate,check_alloc
use store,         only : netot,nspin,inv_pmass,three_netot,opt_cycle,&
&chkpoint_level,isopt,isopt_vmc,isvmc_opt,opt_cycles,rng_restart_safe,&
&complex_wf,which_spin
implicit none
private
public varmin_linjas_main,lsfval,lsfgrad,vm_linjas_method,vm_linjas_it&
&s
integer xyzzyaaaa1,xyzzyaaab1
integer xyzzyaaac1
logical,allocatable :: xyzzyaaad1(:)
character(2),allocatable :: xyzzyaaae1(:)
real(dp) xyzzyaaaf1
real(dp),allocatable :: xyzzyaaag1(:),xyzzyaaah1(:),xyzzyaaai1(:),xyzz&
&yaaaj1(:)
integer xyzzyaaak1,xyzzyaaal1,xyzzyaaam1,xyzzyaaan1
real(dp),allocatable :: xyzzyaaao1(:),xyzzyaaap1(:,:),xyzzyaaaq1(:,:),&
&xyzzyaaar1(:,:)
integer xyzzyaaas1,xyzzyaaat1,xyzzyaaau1,xyzzyaaav1
real(dp) xyzzyaaaw1,xyzzyaaax1,xyzzyaaay1,xyzzyaaaz1
real(dp),allocatable :: xyzzyaaba1(:),xyzzyaabb1(:),xyzzyaabc1(:),xyzz&
&yaabd1(:),xyzzyaabe1(:),xyzzyaabf1(:)
integer xyzzyaabg1,xyzzyaabh1,xyzzyaabi1,xyzzyaabj1
character(10) :: vm_linjas_method
integer vm_linjas_its
real(dp),allocatable :: xyzzyaabk1(:),xyzzyaabl1(:),xyzzyaabm1(:,:),xy&
&zzyaabn1(:),xyzzyaabo1(:,:),xyzzyaabp1(:),xyzzyaabq1(:)
real(dp),allocatable :: xyzzyaabr1(:)
contains
subroutine varmin_linjas_main
implicit none
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2
integer(i64) xyzzyaaae2
real(dp),allocatable :: xyzzyaaaf2(:)
logical xyzzyaaag2(0),xyzzyaaah2(1)
logical :: xyzzyaaai2=.true.
character(20) req_item(5),opt_item(0),req_extra(0),opt_extra(1),dum_co&
&nfig_item(0),extra_item(1)
character(512) errmsg
real(dp),allocatable :: xyzzyaaaj2(:),xyzzyaaak2(:)
logical,allocatable :: xyzzyaaal2(:),xyzzyaaam2(:),xyzzyaaan2(:),xyzzy&
&aaao2(:),xyzzyaaap2(:)
integer xyzzyaaaq2,xyzzyaaar2
logical xyzzyaaas2,xyzzyaaat2
character(10) main_method_code
character(80) method_name
call timer('VARMIN',.true.)
call timer('SETUP',.true.)
if(am_master)then
call wout('Performing variance minimization using Varmin_linjas')
call wout('====================================================')
if(complex_wf)call errstop('VARMIN_LINJAS_MAIN','Varmin_linjas can onl&
&y be used for real wave functions at present.  Please either (i) choo&
&se a different OPT_METHOD, (ii) optimize your wave function at a k_s &
&such that the wave function is real or (iii) ask Neil to implement Va&
&rmin_linjas for complex wave functions.')
endif
call setup_wfn_params(xyzzyaaaa1)
allocate(xyzzyaaaj2(xyzzyaaaa1),xyzzyaaak2(xyzzyaaaa1),xyzzyaaad1(xyzz&
&yaaaa1),xyzzyaaal2(xyzzyaaaa1),xyzzyaaam2(xyzzyaaaa1),xyzzyaaan2(xyzz&
&yaaaa1),xyzzyaaao2(xyzzyaaaa1),xyzzyaaap2(xyzzyaaaa1),xyzzyaaae1(xyzz&
&yaaaa1),stat=xyzzyaaab2)
call check_alloc(xyzzyaaab2,'SETUP_VARMIN','opt_params,...','Reduce nu&
&mber of parameters.')
xyzzyaaaj2=0.d0
xyzzyaaad1=.false.
call get_params(xyzzyaaaj2,xyzzyaaal2,xyzzyaaam2,xyzzyaaan2,xyzzyaaao2&
&,xyzzyaaap2,xyzzyaaae1,rescale=.false.)
xyzzyaaad1=.not.xyzzyaaao2
xyzzyaaab1=count(xyzzyaaao2)
if(xyzzyaaab1==0)call errstop_master('VARMIN_LINJAS_MAIN','No linear J&
&astrow parameters to optimize.')
deallocate(xyzzyaaal2,xyzzyaaam2,xyzzyaaan2,xyzzyaaao2,xyzzyaaap2)
if(am_master)call wout('Number of linear parameters in Jastrow : ' //t&
&rim(i2s(xyzzyaaab1)))
allocate(xyzzyaaaf2(xyzzyaaab1),stat=xyzzyaaab2)
call check_alloc(xyzzyaaab2,'VARMIN_LINJAS_MAIN','')
xyzzyaaaf2=pack(xyzzyaaaj2,.not.xyzzyaaad1)
deallocate(xyzzyaaaj2)
call xyzzyaabs1(xyzzyaaae2)
if(am_master)call wout('G arrays allocated with size           : '//tr&
&im(byte2human(8.d0*dble(xyzzyaaae2))))
xyzzyaaac2=-1
req_item=(/'RELE           ','FIDET          ','LAPDET         ',  'LO&
&CAL_POTENTIAL','NLTOT          '/)
opt_extra='RANDOM'
call load_configs(xyzzyaaac2,'VMC',req_item,opt_item,xyzzyaaag2,req_ex&
&tra,opt_extra,xyzzyaaah2)
if(xyzzyaaah2(1))call put_random_state(random_state_config)
if(am_master)then
call wout('Number of configs loaded               : '//trim(i2s(xyzzya&
&aac2*nnodes)))
endif
call xyzzyaacf1(xyzzyaaac2)
call dismantle_configs
if(am_master)then
if(nnodes>1)then
call wout('Averaging G arrays over nodes.')
else
call wout('Averaging G arrays.')
endif
endif
call xyzzyaach1
xyzzyaaar2=10000
main_method_code=trim(vm_linjas_method)
xyzzyaaas2=.false.
xyzzyaaat2=.true.
select case(trim(vm_linjas_method))
case('MC')
method_name='Monte Carlo line minimization'
xyzzyaaat2=.false.
case('LM')
method_name='simple line minimization'
xyzzyaaar2=40000
xyzzyaaat2=.false.
case('GN')
method_name='Gauss-Newton'
case('SD')
method_name='steepest descents'
case('CG')
method_name='conjugate gradients'
case('BFGS')
method_name='BFGS'
case('SD_MC')
method_name='steepest descents and Monte Carlo'
main_method_code='SD'
xyzzyaaas2=.true.
case('CG_MC')
method_name='conjugate gradients and Monte Carlo'
main_method_code='CG'
xyzzyaaas2=.true.
case('BFGS_MC')
method_name='BFGS and Monte Carlo'
main_method_code='BFGS'
xyzzyaaas2=.true.
case('GN_MC')
method_name='Gauss-Newton and Monte Carlo'
main_method_code='GN'
xyzzyaaas2=.true.
case default
call errstop_master('VARMIN_LINJAS_MAIN','Unknown minimization method &
&"'//trim(vm_linjas_method)//'".')
end select
xyzzyaaaq2=1
if(xyzzyaaas2)xyzzyaaaq2=10
if(vm_linjas_its<0)vm_linjas_its=xyzzyaaar2
if(am_master)then
if(xyzzyaaat2)then
call wout('Constructing Gamma and Lambda arrays.')
else
call wout('Constructing Gamma arrays.')
endif
endif
call xyzzyaabu1(xyzzyaaaa2)
if(am_master)call wout('Gamma arrays allocated with size: '//trim(byte&
&2human(8.d0*dble(xyzzyaaaa2))))
if(xyzzyaaat2)then
call xyzzyaabw1(xyzzyaaaa2)
if(am_master)call wout('Lambda arrays allocated with size: '//trim(byt&
&e2human(8.d0*dble(xyzzyaaaa2))))
endif
call xyzzyaabt1
call timer('SETUP',.false.)
if(am_master)then
call wout('Minimization strategy: '//trim(method_name)//'.')
call wout()
call wout('Parameter set at start of process:')
call put_params(xyzzyaaaf2,xyzzyaaad1,.false.,params_print=xyzzyaaak2)
call labelled_list(xyzzyaaak2(1:xyzzyaaab1),pack(xyzzyaaae1,.not.xyzzy&
&aaad1))
call wout('LSF at start of process: ',lsfval(xyzzyaaaf2))
call wout()
endif
call put_params(xyzzyaaaf2,xyzzyaaad1,.false.,params_print=xyzzyaaak2)
do xyzzyaaad2=1,xyzzyaaaq2
select case(trim(main_method_code))
case('MC')
call xyzzyaaci1(vm_linjas_its,xyzzyaaaf2)
case('LM')
call xyzzyaacn1(vm_linjas_its,xyzzyaaaf2)
case('GN')
call xyzzyaacj1(vm_linjas_its,xyzzyaaaf2)
case('SD')
call xyzzyaack1(vm_linjas_its,xyzzyaaaf2)
case('CG')
call xyzzyaacl1(vm_linjas_its,xyzzyaaaf2)
case('BFGS')
call xyzzyaacm1(vm_linjas_its,xyzzyaaaf2)
end select
if(xyzzyaaas2.and.xyzzyaaad2<xyzzyaaaq2)call xyzzyaaci1(vm_linjas_its,&
&xyzzyaaaf2)
call put_params(xyzzyaaaf2,xyzzyaaad1,.false.,params_print=xyzzyaaak2)
call write_correlation_header(.false.,.true.)
call update_wfn_casl
if(am_master)then
call write_casl(':parameters.casl','parameters.'//trim(i2s(max(1,opt_c&
&ycle)))//'.casl',errmsg)
if(len_trim(errmsg)>0)call errstop_master('VARMIN_LINJAS_MAIN',trim(er&
&rmsg))
endif
enddo
if(am_master)then
call wout('Parameter set at end of process:')
call labelled_list(xyzzyaaak2(1:xyzzyaaab1),pack(xyzzyaaae1,.not.xyzzy&
&aaad1))
call wout('LSF at end of process: ',lsfval(xyzzyaaaf2))
call wout()
call wout('Optimization complete.')
call wout()
endif
call finish_wfn_params
deallocate(xyzzyaaad1,xyzzyaaae1)
deallocate(xyzzyaaaf2)
if(allocated(xyzzyaaag1))call xyzzyaabv1
if(allocated(xyzzyaaap1))call xyzzyaabx1
extra_item(1)='RANDOM'
call init_config_accumulation('OPT',0,dum_config_item,extra_item)
call get_random_state(random_state_config)
call end_config_accumulation(.false.)
if(.not.rng_restart_safe)xyzzyaaai2=.false.
select case(chkpoint_level)
case(-1)
continue
case(0)
if(isopt.or.((isvmc_opt.or.isopt_vmc).and.opt_cycle==opt_cycles))call &
&write_configs(xyzzyaaai2)
case(1,2)
call write_configs(xyzzyaaai2)
end select
call timer('VARMIN',.false.)
end subroutine varmin_linjas_main
subroutine xyzzyaabs1(storage64)
implicit none
integer(i64),intent(out) :: storage64
integer xyzzyaaaa3
integer(i64) xyzzyaaab3,xyzzyaaac3,xyzzyaaad3,xyzzyaaae3,xyzzyaaaf3
xyzzyaaaf3=int(xyzzyaaab1,i64)
xyzzyaaab3=xyzzyaaaf3
if(xyzzyaaab3>int(huge(1),i64))call errstop_master('ALLOCATE_G_DATA','&
&Size of G1 array exceeds largest representable integer. Cannot alloca&
&te arrays - too many parameters.')
xyzzyaabg1=int(xyzzyaaab3)
xyzzyaaac3=(xyzzyaaab1*(xyzzyaaab1+1_i64))/2_i64
if(xyzzyaaac3>int(huge(1),i64))call errstop_master('ALLOCATE_G_DATA','&
&Size of G2 array exceeds largest representable integer. Cannot alloca&
&te arrays - too many parameters.')
xyzzyaabh1=int(xyzzyaaac3)
xyzzyaaad3=(xyzzyaaab1**2_i64*(xyzzyaaab1+1))/2_i64
if(xyzzyaaad3>int(huge(1),i64))call errstop_master('ALLOCATE_G_DATA','&
&Size of G3 array exceeds largest representable integer. Cannot alloca&
&te arrays - too many parameters.')
xyzzyaabi1=int(xyzzyaaad3)
xyzzyaaae3=(((xyzzyaaab1*(xyzzyaaab1+1_i64))/2_i64) *((xyzzyaaab1*(xyz&
&zyaaab1+1_i64))/2_i64+1_i64))/2_i64
if(xyzzyaaae3>int(huge(1),i64))call errstop_master('ALLOCATE_G_DATA','&
&Size of G4 array exceeds largest representable integer. Cannot alloca&
&te arrays - too many parameters.')
xyzzyaabj1=int(xyzzyaaae3)
allocate(xyzzyaaba1(xyzzyaabg1),xyzzyaabb1(xyzzyaabh1),xyzzyaabc1(xyzz&
&yaabg1),xyzzyaabd1(xyzzyaabh1),xyzzyaabe1(xyzzyaabi1),xyzzyaabf1(xyzz&
&yaabj1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'ALLOCATE_G_DATA','G arrays')
storage64=3+2*(xyzzyaaab3+xyzzyaaac3)+xyzzyaaad3+xyzzyaaae3
end subroutine xyzzyaabs1
subroutine xyzzyaabt1
implicit none
deallocate(xyzzyaaba1,xyzzyaabb1,xyzzyaabc1,xyzzyaabd1,xyzzyaabe1,xyzz&
&yaabf1)
end subroutine xyzzyaabt1
subroutine xyzzyaabu1(storage)
implicit none
integer,intent(out) :: storage
integer xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5,xyzzyaaad5,xyzzyaaae5,xyzzyaa&
&af5,xyzzyaaag5,xyzzyaaah5,xyzzyaaai5,xyzzyaaaj5
real(dp) xyzzyaaak5
if(.not.allocated(xyzzyaaag1))then
xyzzyaaak1=xyzzyaaab1
xyzzyaaal1=choose(xyzzyaaab1+1,2)
xyzzyaaam1=choose(xyzzyaaab1+2,3)
xyzzyaaan1=choose(xyzzyaaab1+3,4)
allocate(xyzzyaaag1(xyzzyaaak1),xyzzyaaah1(xyzzyaaal1),xyzzyaaai1(xyzz&
&yaaam1),xyzzyaaaj1(xyzzyaaan1),stat=xyzzyaaaj5)
call check_alloc(xyzzyaaaj5,'CONSTRUCT_GAMMA','Gamma <1>')
else
if(xyzzyaaak1/=xyzzyaaab1)then
call xyzzyaabv1
xyzzyaaak1=xyzzyaaab1
xyzzyaaal1=choose(xyzzyaaab1+1,2)
xyzzyaaam1=choose(xyzzyaaab1+2,3)
xyzzyaaan1=choose(xyzzyaaab1+3,4)
allocate(xyzzyaaag1(xyzzyaaak1),xyzzyaaah1(xyzzyaaal1),xyzzyaaai1(xyzz&
&yaaam1),xyzzyaaaj1(xyzzyaaan1),stat=xyzzyaaaj5)
call check_alloc(xyzzyaaaj5,'CONSTRUCT_GAMMA','Gamma <2>')
endif
endif
storage=1+xyzzyaaak1+xyzzyaaal1+xyzzyaaam1+xyzzyaaan1
xyzzyaaak5=xyzzyaaaz1-0.5d0*xyzzyaaax1
xyzzyaaaf1=xyzzyaaay1-xyzzyaaak5**2
call dcopy(xyzzyaaab1,xyzzyaabc1(1),1,xyzzyaaag1(1),1)
call daxpy(xyzzyaaab1,xyzzyaaak5,xyzzyaaba1(1),1,xyzzyaaag1(1),1)
xyzzyaaah1=0.d0
xyzzyaaai1=0.d0
xyzzyaaaj1=0.d0
do xyzzyaaaa5=1,xyzzyaaab1
do xyzzyaaab5=1,xyzzyaaab1
xyzzyaaae5=xyzzyaaby1(xyzzyaaaa5,xyzzyaaab5)
xyzzyaaaf5=xyzzyaaae5
xyzzyaaah1(xyzzyaaae5)=xyzzyaaah1(xyzzyaaae5)+xyzzyaabd1(xyzzyaaaf5)-0&
&.25d0*xyzzyaaba1(xyzzyaaaa5)*xyzzyaaba1(xyzzyaaab5)+xyzzyaabb1(xyzzya&
&aaf5)*xyzzyaaak5
do xyzzyaaac5=1,xyzzyaaab1
xyzzyaaae5=xyzzyaabz1(xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5)
xyzzyaaah5=xyzzyaacb1(xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5)
xyzzyaaai1(xyzzyaaae5)=xyzzyaaai1(xyzzyaaae5)+xyzzyaabe1(xyzzyaaah5)-0&
&.5d0*xyzzyaabb1(xyzzyaaaf5)*xyzzyaaba1(xyzzyaaac5)
do xyzzyaaad5=1,xyzzyaaab1
xyzzyaaae5=xyzzyaaca1(xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5,xyzzyaaad5)
xyzzyaaai5=xyzzyaacc1(xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5,xyzzyaaad5)
xyzzyaaag5=xyzzyaaby1(xyzzyaaac5,xyzzyaaad5)
xyzzyaaaj1(xyzzyaaae5)=xyzzyaaaj1(xyzzyaaae5)+xyzzyaabf1(xyzzyaaai5)-0&
&.25d0*xyzzyaabb1(xyzzyaaaf5)*xyzzyaabb1(xyzzyaaag5)
enddo
enddo
enddo
enddo
end subroutine xyzzyaabu1
subroutine xyzzyaabv1
implicit none
deallocate(xyzzyaaag1,xyzzyaaah1,xyzzyaaai1,xyzzyaaaj1)
xyzzyaaak1=0
xyzzyaaal1=0
xyzzyaaam1=0
xyzzyaaan1=0
end subroutine xyzzyaabv1
subroutine xyzzyaabw1(storage)
implicit none
integer,intent(out) :: storage
integer xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7,xyzzyaaad7,xyzzyaaae7,xyzzyaa&
&af7,xyzzyaaag7,xyzzyaaah7,xyzzyaaai7,xyzzyaaaj7,xyzzyaaak7,xyzzyaaal7&
&,xyzzyaaam7,xyzzyaaan7,xyzzyaaao7
real(dp) xyzzyaaap7
if(.not.allocated(xyzzyaaap1))then
xyzzyaaat1=xyzzyaaab1
xyzzyaaau1=choose(xyzzyaaab1+1,2)
xyzzyaaav1=choose(xyzzyaaab1+2,3)
allocate(xyzzyaaao1(xyzzyaaab1),xyzzyaaap1(xyzzyaaab1,xyzzyaaat1),xyzz&
&yaaaq1(xyzzyaaab1,xyzzyaaau1),xyzzyaaar1(xyzzyaaab1,xyzzyaaav1),stat=&
&xyzzyaaao7)
call check_alloc(xyzzyaaao7,'CONSTRUCT_LAMBDA','Lambda arrays <1>')
xyzzyaaas1=xyzzyaaab1
xyzzyaaat1=xyzzyaaab1*xyzzyaaat1
xyzzyaaau1=xyzzyaaab1*xyzzyaaau1
xyzzyaaav1=xyzzyaaab1*xyzzyaaav1
else
if(xyzzyaaas1/=xyzzyaaab1)then
call xyzzyaabx1
xyzzyaaat1=xyzzyaaab1
xyzzyaaau1=choose(xyzzyaaab1+1,2)
xyzzyaaav1=choose(xyzzyaaab1+2,3)
allocate(xyzzyaaao1(xyzzyaaab1),xyzzyaaap1(xyzzyaaab1,xyzzyaaat1),xyzz&
&yaaaq1(xyzzyaaab1,xyzzyaaau1),xyzzyaaar1(xyzzyaaab1,xyzzyaaav1),stat=&
&xyzzyaaao7)
call check_alloc(xyzzyaaao7,'CONSTRUCT_LAMBDA','Lambda arrays <2>')
xyzzyaaas1=xyzzyaaab1
xyzzyaaat1=xyzzyaaab1*xyzzyaaat1
xyzzyaaau1=xyzzyaaab1*xyzzyaaau1
xyzzyaaav1=xyzzyaaab1*xyzzyaaav1
endif
endif
storage=xyzzyaaas1+xyzzyaaat1+xyzzyaaau1+xyzzyaaav1
xyzzyaaap7=xyzzyaaaz1-0.5d0*xyzzyaaax1
xyzzyaaaq1=0.d0
xyzzyaaar1=0.d0
do xyzzyaaaa7=1,xyzzyaaab1
do xyzzyaaab7=1,xyzzyaaab1
xyzzyaaah7=xyzzyaaby1(xyzzyaaaa7,xyzzyaaab7)
do xyzzyaaac7=1,xyzzyaaab1
xyzzyaaae7=xyzzyaabz1(xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7)
xyzzyaaam7=xyzzyaaby1(xyzzyaaab7,xyzzyaaac7)
do xyzzyaaad7=1,xyzzyaaab1
xyzzyaaak7=xyzzyaacc1(xyzzyaaad7,xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7)
xyzzyaaal7=xyzzyaacc1(xyzzyaaaa7,xyzzyaaab7,xyzzyaaad7,xyzzyaaac7)
xyzzyaaag7=xyzzyaaby1(xyzzyaaad7,xyzzyaaaa7)
xyzzyaaan7=xyzzyaaby1(xyzzyaaad7,xyzzyaaac7)
xyzzyaaar1(xyzzyaaad7,xyzzyaaae7)=xyzzyaaar1(xyzzyaaad7,xyzzyaaae7)+2*&
&(xyzzyaabf1(xyzzyaaak7) +xyzzyaabf1(xyzzyaaal7))-0.5d0*(xyzzyaabb1(xy&
&zzyaaag7)*xyzzyaabb1(xyzzyaaam7) +xyzzyaabb1(xyzzyaaah7)*xyzzyaabb1(x&
&yzzyaaan7))
enddo
enddo
xyzzyaaae7=xyzzyaaby1(xyzzyaaaa7,xyzzyaaab7)
do xyzzyaaad7=1,xyzzyaaab1
xyzzyaaai7=xyzzyaacb1(xyzzyaaad7,xyzzyaaaa7,xyzzyaaab7)
xyzzyaaaj7=xyzzyaacb1(xyzzyaaaa7,xyzzyaaab7,xyzzyaaad7)
xyzzyaaag7=xyzzyaaby1(xyzzyaaad7,xyzzyaaaa7)
xyzzyaaaq1(xyzzyaaad7,xyzzyaaae7)=xyzzyaaaq1(xyzzyaaad7,xyzzyaaae7)+2*&
&xyzzyaabe1(xyzzyaaai7) +xyzzyaabe1(xyzzyaaaj7)-xyzzyaabb1(xyzzyaaag7)&
&*xyzzyaaba1(xyzzyaaab7) -0.5d0*xyzzyaabb1(xyzzyaaah7)*xyzzyaaba1(xyzz&
&yaaad7)
enddo
enddo
do xyzzyaaad7=1,xyzzyaaab1
xyzzyaaaf7=xyzzyaaby1(xyzzyaaaa7,xyzzyaaad7)
xyzzyaaap1(xyzzyaaad7,xyzzyaaaa7)=2*(xyzzyaabd1(xyzzyaaaf7)+xyzzyaabb1&
&(xyzzyaaaf7)*xyzzyaaap7) -0.5d0*xyzzyaaba1(xyzzyaaaa7)*xyzzyaaba1(xyz&
&zyaaad7)
enddo
enddo
call dcopy(xyzzyaaab1,xyzzyaabc1(1),1,xyzzyaaao1(1),1)
call daxpy(xyzzyaaab1,xyzzyaaap7,xyzzyaaba1(1),1,xyzzyaaao1(1),1)
end subroutine xyzzyaabw1
subroutine xyzzyaabx1
implicit none
deallocate(xyzzyaaao1,xyzzyaaap1,xyzzyaaaq1,xyzzyaaar1)
xyzzyaaas1=0
xyzzyaaat1=0
xyzzyaaau1=0
xyzzyaaav1=0
end subroutine xyzzyaabx1
integer function xyzzyaaby1(i,j)
implicit none
integer,intent(in) :: i,j
integer xyzzyaaaa9,xyzzyaaab9
if(i<j)then
xyzzyaaaa9=i
xyzzyaaab9=j
else
xyzzyaaaa9=j
xyzzyaaab9=i
endif
xyzzyaaby1=((xyzzyaaaa9-1)*(2*xyzzyaaab1+2-xyzzyaaaa9))/2+xyzzyaaab9-x&
&yzzyaaaa9+1
end function xyzzyaaby1
integer function xyzzyaabz1(i,j,k)
implicit none
integer,intent(in) :: i,j,k
integer xyzzyaaaa10,xyzzyaaab10,xyzzyaaac10
if(k>j)then
if(j>=i)then
xyzzyaaaa10=i
xyzzyaaab10=j
xyzzyaaac10=k
else
if(i>k)then
xyzzyaaaa10=j
xyzzyaaab10=k
xyzzyaaac10=i
else
xyzzyaaaa10=j
xyzzyaaab10=i
xyzzyaaac10=k
endif
endif
else
if(i>=j)then
xyzzyaaaa10=k
xyzzyaaab10=j
xyzzyaaac10=i
else
if(i>k)then
xyzzyaaaa10=k
xyzzyaaab10=i
xyzzyaaac10=j
else
xyzzyaaaa10=i
xyzzyaaab10=k
xyzzyaaac10=j
endif
endif
endif
xyzzyaabz1=(xyzzyaaaa10*(3*xyzzyaaab1*(xyzzyaaab1+2)+2 +xyzzyaaaa10*(x&
&yzzyaaaa10-3*xyzzyaaab1-3))+3*(xyzzyaaab1*(2*xyzzyaaab10-xyzzyaaab1-3&
&) +xyzzyaaab10*(1-xyzzyaaab10)+2*xyzzyaaac10))/6
end function xyzzyaabz1
integer function xyzzyaaca1(i,j,k,l)
implicit none
integer,intent(in) :: i,j,k,l
integer xyzzyaaaa11,xyzzyaaab11,xyzzyaaac11,xyzzyaaad11
if(i<=j)then
if(j<=k)then
if(k<=l)then
xyzzyaaaa11=i
xyzzyaaab11=j
xyzzyaaac11=k
xyzzyaaad11=l
else
if(j<=l)then
xyzzyaaaa11=i
xyzzyaaab11=j
xyzzyaaac11=l
xyzzyaaad11=k
else
if(i<l)then
xyzzyaaaa11=i
xyzzyaaab11=l
xyzzyaaac11=j
xyzzyaaad11=k
else
xyzzyaaaa11=l
xyzzyaaab11=i
xyzzyaaac11=j
xyzzyaaad11=k
endif
endif
endif
else
if(i<=k)then
if(j<=l)then
xyzzyaaaa11=i
xyzzyaaab11=k
xyzzyaaac11=j
xyzzyaaad11=l
else
if(k<=l)then
xyzzyaaaa11=i
xyzzyaaab11=k
xyzzyaaac11=l
xyzzyaaad11=j
else
if(i<l)then
xyzzyaaaa11=i
xyzzyaaab11=l
xyzzyaaac11=k
xyzzyaaad11=j
else
xyzzyaaaa11=l
xyzzyaaab11=i
xyzzyaaac11=k
xyzzyaaad11=j
endif
endif
endif
else
if(j<=l)then
xyzzyaaaa11=k
xyzzyaaab11=i
xyzzyaaac11=j
xyzzyaaad11=l
else
if(i<=l)then
xyzzyaaaa11=k
xyzzyaaab11=i
xyzzyaaac11=l
xyzzyaaad11=j
else
if(k<l)then
xyzzyaaaa11=k
xyzzyaaab11=l
xyzzyaaac11=i
xyzzyaaad11=j
else
xyzzyaaaa11=l
xyzzyaaab11=k
xyzzyaaac11=i
xyzzyaaad11=j
endif
endif
endif
endif
endif
else
if(i<=k)then
if(k<=l)then
xyzzyaaaa11=j
xyzzyaaab11=i
xyzzyaaac11=k
xyzzyaaad11=l
else
if(i<=l)then
xyzzyaaaa11=j
xyzzyaaab11=i
xyzzyaaac11=l
xyzzyaaad11=k
else
if(j<l)then
xyzzyaaaa11=j
xyzzyaaab11=l
xyzzyaaac11=i
xyzzyaaad11=k
else
xyzzyaaaa11=l
xyzzyaaab11=j
xyzzyaaac11=i
xyzzyaaad11=k
endif
endif
endif
else
if(j<=k)then
if(i<=l)then
xyzzyaaaa11=j
xyzzyaaab11=k
xyzzyaaac11=i
xyzzyaaad11=l
else
if(k<=l)then
xyzzyaaaa11=j
xyzzyaaab11=k
xyzzyaaac11=l
xyzzyaaad11=i
else
if(j<l)then
xyzzyaaaa11=j
xyzzyaaab11=l
xyzzyaaac11=k
xyzzyaaad11=i
else
xyzzyaaaa11=l
xyzzyaaab11=j
xyzzyaaac11=k
xyzzyaaad11=i
endif
endif
endif
else
if(i<=l)then
xyzzyaaaa11=k
xyzzyaaab11=j
xyzzyaaac11=i
xyzzyaaad11=l
else
if(j<=l)then
xyzzyaaaa11=k
xyzzyaaab11=j
xyzzyaaac11=l
xyzzyaaad11=i
else
if(k<l)then
xyzzyaaaa11=k
xyzzyaaab11=l
xyzzyaaac11=j
xyzzyaaad11=i
else
xyzzyaaaa11=l
xyzzyaaab11=k
xyzzyaaac11=j
xyzzyaaad11=i
endif
endif
endif
endif
endif
endif
xyzzyaaca1=(xyzzyaaab1*(xyzzyaaab1*((4*xyzzyaaaa11-4)*xyzzyaaab1 +6*xy&
&zzyaaaa11*(3-xyzzyaaaa11)-24+12*xyzzyaaab11)+12*xyzzyaaab11*(2-xyzzya&
&aab11)+24*xyzzyaaac11-44+xyzzyaaaa11*(22+xyzzyaaaa11*(4*xyzzyaaaa11-1&
&8))) +24*xyzzyaaad11+12*xyzzyaaac11*(1-xyzzyaaac11)+4*xyzzyaaab11*(xy&
&zzyaaab11*(xyzzyaaab11-3)+2)+xyzzyaaaa11*(xyzzyaaaa11*(xyzzyaaaa11*(6&
&-xyzzyaaaa11)-11)+6))/24
end function xyzzyaaca1
integer function xyzzyaacb1(i,j,kp)
implicit none
integer,intent(in) :: i,j,kp
integer xyzzyaaaa12,xyzzyaaab12
if(j<i)then
xyzzyaaaa12=j
xyzzyaaab12=i
else
xyzzyaaaa12=i
xyzzyaaab12=j
endif
xyzzyaacb1=(xyzzyaaab1*(2*(xyzzyaaaa12-1)*xyzzyaaab1-2+2*xyzzyaaab12+x&
&yzzyaaaa12*(1-xyzzyaaaa12)))/2+kp
end function xyzzyaacb1
integer function xyzzyaacc1(i,j,k,l)
implicit none
integer,intent(in) :: i,j,k,l
integer xyzzyaaaa13,xyzzyaaab13,xyzzyaaac13,xyzzyaaad13,xyzzyaaae13,xy&
&zzyaaaf13,xyzzyaaag13
if(j<i)then
xyzzyaaaa13=j
xyzzyaaab13=i
else
xyzzyaaaa13=i
xyzzyaaab13=j
endif
xyzzyaaae13=((xyzzyaaaa13-1)*(2*xyzzyaaab1+2-xyzzyaaaa13))/2+xyzzyaaab&
&13-xyzzyaaaa13+1
if(l<k)then
xyzzyaaac13=l
xyzzyaaad13=k
else
xyzzyaaac13=k
xyzzyaaad13=l
endif
xyzzyaaaf13=((xyzzyaaac13-1)*(2*xyzzyaaab1+2-xyzzyaaac13))/2+xyzzyaaad&
&13-xyzzyaaac13+1
if(xyzzyaaaf13<xyzzyaaae13)then
xyzzyaaag13=xyzzyaaaa13
xyzzyaaaa13=xyzzyaaac13
xyzzyaaac13=xyzzyaaag13
xyzzyaaag13=xyzzyaaab13
xyzzyaaab13=xyzzyaaad13
xyzzyaaad13=xyzzyaaag13
endif
xyzzyaacc1=(xyzzyaaab1*(xyzzyaaab1*(((xyzzyaaaa13-1)*xyzzyaaab1+xyzzya&
&aab13)*4 +xyzzyaaaa13*(14-6*xyzzyaaaa13)-12)+xyzzyaaaa13*(xyzzyaaaa13&
&*(4*xyzzyaaaa13-10)+10-8*xyzzyaaab13)-16+12*xyzzyaaab13+8*xyzzyaaac13&
&) +xyzzyaaaa13*(xyzzyaaaa13*(xyzzyaaaa13*(2-xyzzyaaaa13)+4*xyzzyaaab1&
&3-3)+2-4*xyzzyaaab13)+4*(2*xyzzyaaad13+xyzzyaaac13*(1-xyzzyaaac13)+xy&
&zzyaaab13*(1-xyzzyaaab13)))/8
end function xyzzyaacc1
real(dp) function xyzzyaacd1(a1,a2,a3,a4)
implicit none
real(dp),intent(in) :: a1,a2,a3,a4
real(dp) xyzzyaaaa14,xyzzyaaab14,xyzzyaaac14,xyzzyaaad14,xyzzyaaae14,x&
&yzzyaaaf14,xyzzyaaag14,xyzzyaaah14,xyzzyaaai14,xyzzyaaaj14,xyzzyaaak1&
&4,xyzzyaaal14,xyzzyaaam14,xyzzyaaan14,xyzzyaaao14,xyzzyaaap14
if(a4<=0.d0)call errstop('min_quartic','Quartic has no lower bound!')
xyzzyaaaa14=1.d0/a4
xyzzyaaab14=0.75d0*a3*xyzzyaaaa14
xyzzyaaac14=0.5d0*a2*xyzzyaaaa14
xyzzyaaad14=0.25d0*a1*xyzzyaaaa14
xyzzyaaae14=ninth*(xyzzyaaab14**2-3*xyzzyaaac14)
xyzzyaaaf14=rec_54*(2*xyzzyaaab14**3-9*xyzzyaaab14*xyzzyaaac14+27*xyzz&
&yaaad14)
if(xyzzyaaaf14**2<xyzzyaaae14**3)then
xyzzyaaag14=sqrt(xyzzyaaae14)
xyzzyaaap14=acos(xyzzyaaaf14/xyzzyaaag14**3)
xyzzyaaah14=-2*xyzzyaaag14*cos(third*xyzzyaaap14)-third*xyzzyaaab14
xyzzyaaai14=-2*xyzzyaaag14*cos(third*(xyzzyaaap14+twopi))-third*xyzzya&
&aab14
xyzzyaaaj14=-2*xyzzyaaag14*cos(third*(xyzzyaaap14-twopi))-third*xyzzya&
&aab14
xyzzyaaak14=(((a4*xyzzyaaah14+a3)*xyzzyaaah14+a2)*xyzzyaaah14+a1)*xyzz&
&yaaah14
xyzzyaaal14=(((a4*xyzzyaaai14+a3)*xyzzyaaai14+a2)*xyzzyaaai14+a1)*xyzz&
&yaaai14
xyzzyaaam14=(((a4*xyzzyaaaj14+a3)*xyzzyaaaj14+a2)*xyzzyaaaj14+a1)*xyzz&
&yaaaj14
if(xyzzyaaal14<xyzzyaaam14)then
if(xyzzyaaak14<xyzzyaaal14)then
xyzzyaacd1=xyzzyaaah14
else
xyzzyaacd1=xyzzyaaai14
endif
else
if(xyzzyaaak14<xyzzyaaam14)then
xyzzyaacd1=xyzzyaaah14
else
xyzzyaacd1=xyzzyaaaj14
endif
endif
else
xyzzyaaan14=-sign(1.d0,xyzzyaaaf14)*(abs(xyzzyaaaf14)+sqrt(xyzzyaaaf14&
&**2-xyzzyaaae14**3))**third
if(xyzzyaaan14/=0.d0)then
xyzzyaaao14=xyzzyaaae14/xyzzyaaan14
else
xyzzyaaao14=0.d0
endif
xyzzyaacd1=(xyzzyaaan14+xyzzyaaao14)-third*xyzzyaaab14
endif
end function xyzzyaacd1
subroutine xyzzyaace1(a_vect,b_vect,omega0,omega1,omega2,omega3,omega4&
&)
implicit none
real(dp),intent(in) :: a_vect(xyzzyaaab1),b_vect(xyzzyaaab1)
real(dp),intent(out) :: omega0,omega1,omega2,omega3,omega4
integer xyzzyaaaa15,xyzzyaaab15,xyzzyaaac15,xyzzyaaad15,xyzzyaaae15,xy&
&zzyaaaf15,xyzzyaaag15
real(dp) xyzzyaaah15,xyzzyaaai15,xyzzyaaaj15,xyzzyaaak15,xyzzyaaal15,x&
&yzzyaaam15,xyzzyaaan15,xyzzyaaao15,xyzzyaaap15,xyzzyaaaq15,xyzzyaaar1&
&5,xyzzyaaas15,xyzzyaaat15,xyzzyaaau15,xyzzyaaav15,xyzzyaaaw15,xyzzyaa&
&ax15,xyzzyaaay15,xyzzyaaaz15,xyzzyaaba15,xyzzyaabb15,xyzzyaabc15,xyzz&
&yaabd15,xyzzyaabe15,xyzzyaabf15,xyzzyaabg15,xyzzyaabh15,xyzzyaabi15
omega4=0.d0
omega0=xyzzyaaaf1
xyzzyaaan15=0.d0
xyzzyaaao15=0.d0
xyzzyaaap15=0.d0
xyzzyaaaq15=0.d0
xyzzyaaaz15=0.d0
xyzzyaaba15=0.d0
xyzzyaabb15=0.d0
xyzzyaabc15=0.d0
xyzzyaabd15=0.d0
xyzzyaabe15=0.d0
xyzzyaabf15=0.d0
xyzzyaabg15=0.d0
xyzzyaabh15=0.d0
xyzzyaabi15=0.d0
xyzzyaaae15=0
xyzzyaaaf15=0
xyzzyaaag15=0
do xyzzyaaaa15=1,xyzzyaaab1
xyzzyaaah15=0.d0
xyzzyaaak15=xyzzyaaag1(xyzzyaaaa15)
xyzzyaaar15=0.d0
xyzzyaaas15=0.d0
xyzzyaaat15=0.d0
xyzzyaaau15=0.d0
xyzzyaaav15=0.d0
xyzzyaaaw15=0.d0
do xyzzyaaab15=xyzzyaaaa15,xyzzyaaab1
xyzzyaaae15=xyzzyaaae15+1
xyzzyaaai15=0.d0
xyzzyaaal15=xyzzyaaah1(xyzzyaaae15)
xyzzyaaax15=0.d0
xyzzyaaay15=0.d0
do xyzzyaaac15=xyzzyaaab15,xyzzyaaab1
xyzzyaaaf15=xyzzyaaaf15+1
xyzzyaaaj15=0.d0
xyzzyaaam15=xyzzyaaai1(xyzzyaaaf15)
do xyzzyaaad15=xyzzyaaac15,xyzzyaaab1
xyzzyaaag15=xyzzyaaag15+1
xyzzyaaaj15=xyzzyaaaj15+xyzzyaaaj1(xyzzyaaag15)*b_vect(xyzzyaaad15)
xyzzyaaam15=xyzzyaaam15+xyzzyaaaj1(xyzzyaaag15)*a_vect(xyzzyaaad15)
enddo
xyzzyaaai15=xyzzyaaai15+xyzzyaaaj15*b_vect(xyzzyaaac15)
xyzzyaaal15=xyzzyaaal15+xyzzyaaam15*a_vect(xyzzyaaac15)
xyzzyaaax15=xyzzyaaax15+xyzzyaaaj15*a_vect(xyzzyaaac15)
xyzzyaaay15=xyzzyaaay15+xyzzyaaam15*b_vect(xyzzyaaac15)
enddo
xyzzyaaah15=xyzzyaaah15+xyzzyaaai15*b_vect(xyzzyaaab15)
xyzzyaaak15=xyzzyaaak15+xyzzyaaal15*a_vect(xyzzyaaab15)
xyzzyaaar15=xyzzyaaar15+xyzzyaaai15*a_vect(xyzzyaaab15)
xyzzyaaas15=xyzzyaaas15+xyzzyaaax15*b_vect(xyzzyaaab15)
xyzzyaaat15=xyzzyaaat15+xyzzyaaay15*b_vect(xyzzyaaab15)
xyzzyaaau15=xyzzyaaau15+xyzzyaaax15*a_vect(xyzzyaaab15)
xyzzyaaav15=xyzzyaaav15+xyzzyaaay15*a_vect(xyzzyaaab15)
xyzzyaaaw15=xyzzyaaaw15+xyzzyaaal15*b_vect(xyzzyaaab15)
enddo
omega4=omega4+xyzzyaaah15*b_vect(xyzzyaaaa15)
omega0=omega0+xyzzyaaak15*a_vect(xyzzyaaaa15)
xyzzyaaan15=xyzzyaaan15+xyzzyaaah15*a_vect(xyzzyaaaa15)
xyzzyaaao15=xyzzyaaao15+xyzzyaaar15*b_vect(xyzzyaaaa15)
xyzzyaaap15=xyzzyaaap15+xyzzyaaas15*b_vect(xyzzyaaaa15)
xyzzyaaaq15=xyzzyaaaq15+xyzzyaaat15*b_vect(xyzzyaaaa15)
xyzzyaaaz15=xyzzyaaaz15+xyzzyaaar15*a_vect(xyzzyaaaa15)
xyzzyaaba15=xyzzyaaba15+xyzzyaaas15*a_vect(xyzzyaaaa15)
xyzzyaabb15=xyzzyaabb15+xyzzyaaat15*a_vect(xyzzyaaaa15)
xyzzyaabc15=xyzzyaabc15+xyzzyaaau15*b_vect(xyzzyaaaa15)
xyzzyaabd15=xyzzyaabd15+xyzzyaaav15*b_vect(xyzzyaaaa15)
xyzzyaabe15=xyzzyaabe15+xyzzyaaaw15*b_vect(xyzzyaaaa15)
xyzzyaabf15=xyzzyaabf15+xyzzyaaak15*b_vect(xyzzyaaaa15)
xyzzyaabg15=xyzzyaabg15+xyzzyaaaw15*a_vect(xyzzyaaaa15)
xyzzyaabh15=xyzzyaabh15+xyzzyaaav15*a_vect(xyzzyaaaa15)
xyzzyaabi15=xyzzyaabi15+xyzzyaaau15*a_vect(xyzzyaaaa15)
enddo
omega3=xyzzyaaan15+xyzzyaaao15+xyzzyaaap15+xyzzyaaaq15
omega2=xyzzyaaaz15+xyzzyaaba15+xyzzyaabb15+xyzzyaabc15+xyzzyaabd15+xyz&
&zyaabe15
omega1=xyzzyaabf15+xyzzyaabg15+xyzzyaabh15+xyzzyaabi15
end subroutine xyzzyaace1
real(dp) function lsfval(alpha)
implicit none
real(dp),intent(in) :: alpha(xyzzyaaab1)
integer xyzzyaaaa16,xyzzyaaab16,xyzzyaaac16,xyzzyaaad16,xyzzyaaae16,xy&
&zzyaaaf16,xyzzyaaag16
real(dp) xyzzyaaah16,xyzzyaaai16,xyzzyaaaj16
lsfval=xyzzyaaaf1
xyzzyaaae16=0
xyzzyaaaf16=0
xyzzyaaag16=0
do xyzzyaaaa16=1,xyzzyaaab1
xyzzyaaah16=xyzzyaaag1(xyzzyaaaa16)
do xyzzyaaab16=xyzzyaaaa16,xyzzyaaab1
xyzzyaaae16=xyzzyaaae16+1
xyzzyaaai16=xyzzyaaah1(xyzzyaaae16)
do xyzzyaaac16=xyzzyaaab16,xyzzyaaab1
xyzzyaaaf16=xyzzyaaaf16+1
xyzzyaaaj16=xyzzyaaai1(xyzzyaaaf16)
do xyzzyaaad16=xyzzyaaac16,xyzzyaaab1
xyzzyaaag16=xyzzyaaag16+1
xyzzyaaaj16=xyzzyaaaj16+xyzzyaaaj1(xyzzyaaag16)*alpha(xyzzyaaad16)
enddo
xyzzyaaai16=xyzzyaaai16+xyzzyaaaj16*alpha(xyzzyaaac16)
enddo
xyzzyaaah16=xyzzyaaah16+xyzzyaaai16*alpha(xyzzyaaab16)
enddo
lsfval=lsfval+xyzzyaaah16*alpha(xyzzyaaaa16)
enddo
end function lsfval
subroutine lsfgrad(alpha,dlsf_dalpha)
implicit none
real(dp),intent(in) :: alpha(xyzzyaaab1)
real(dp),intent(out) :: dlsf_dalpha(xyzzyaaab1)
integer xyzzyaaaa17,xyzzyaaab17,xyzzyaaac17,xyzzyaaad17,xyzzyaaae17
real(dp) :: xyzzyaaaf17(xyzzyaaab1),xyzzyaaag17(xyzzyaaab1)
call dcopy(xyzzyaaab1,xyzzyaaao1(1),1,dlsf_dalpha(1),1)
xyzzyaaad17=0
xyzzyaaae17=0
do xyzzyaaaa17=1,xyzzyaaab1
call dcopy(xyzzyaaab1,xyzzyaaap1(1,xyzzyaaaa17),1,xyzzyaaaf17(1),1)
do xyzzyaaab17=xyzzyaaaa17,xyzzyaaab1
xyzzyaaad17=xyzzyaaad17+1
call dcopy(xyzzyaaab1,xyzzyaaaq1(1,xyzzyaaad17),1,xyzzyaaag17(1),1)
do xyzzyaaac17=xyzzyaaab17,xyzzyaaab1
xyzzyaaae17=xyzzyaaae17+1
call daxpy(xyzzyaaab1,alpha(xyzzyaaac17),xyzzyaaar1(1,xyzzyaaae17),1,x&
&yzzyaaag17(1),1)
enddo
call daxpy(xyzzyaaab1,alpha(xyzzyaaab17),xyzzyaaag17(1),1,xyzzyaaaf17(&
&1),1)
enddo
call daxpy(xyzzyaaab1,alpha(xyzzyaaaa17),xyzzyaaaf17(1),1,dlsf_dalpha(&
&1),1)
enddo
end subroutine lsfgrad
subroutine xyzzyaacf1(opt_nconfig)
implicit none
integer,intent(in) :: opt_nconfig
integer xyzzyaaaa18,xyzzyaaab18
xyzzyaaaw1=0.d0
xyzzyaaax1=0.d0
xyzzyaaba1=0.d0
xyzzyaabb1=0.d0
xyzzyaaay1=0.d0
xyzzyaabc1=0.d0
xyzzyaabd1=0.d0
xyzzyaabe1=0.d0
xyzzyaabf1=0.d0
xyzzyaaaz1=0.d0
xyzzyaaac1=0
call scratch_protect(xyzzyaaac1)
call energy_scratch_request(xyzzyaaac1)
call ederiv_scratch_request(xyzzyaaac1,xyzzyaaad1)
call setup_scratch
call which_scratch(xyzzyaaac1)
call setup_wfn_utils
call setup_energy_utils
call setup_storage(opt_nconfig,xyzzyaaad1)
call setup_storage_energy(opt_nconfig)
allocate(xyzzyaabm1(three_netot,xyzzyaaab1),xyzzyaabk1(xyzzyaabg1),xyz&
&zyaabl1(xyzzyaabh1),xyzzyaabn1(three_netot),xyzzyaabo1(3,netot),xyzzy&
&aabp1(three_netot),xyzzyaabq1(three_netot),stat=xyzzyaaaa18)
call check_alloc(xyzzyaaaa18,'CONSTRUCT_COEFFS_FROM_CONFIGS','grad_f,.&
&..')
xyzzyaabo1=0.d0
xyzzyaabp1=0.d0
xyzzyaabq1=0.d0
if(any(inv_pmass/=1.d0).and..not.allocated(xyzzyaabr1))then
allocate(xyzzyaabr1(nspin),stat=xyzzyaaaa18)
call check_alloc(xyzzyaaaa18,'CONSTRUCT_COEFFS_FROM_CONFIGS','rec_root&
&,...')
xyzzyaabr1(1:nspin)=sqrt(inv_pmass(1:nspin))
endif
do xyzzyaaab18=1,opt_nconfig
call loop_time_estimate('Constructing G arrays.',xyzzyaaab18,opt_nconf&
&ig)
call load_from_storage(xyzzyaaac1,xyzzyaaab18)
call load_from_storage_energy(xyzzyaaac1,xyzzyaaab18)
call xyzzyaacg1(is=xyzzyaaac1)
enddo
call loop_time_estimate('Done.')
if(opt_nconfig/=nint(xyzzyaaaw1))call errstop('CONSTRUCT_COEFFS_FROM_C&
&ONFIGS','Problem: weight differs from number of configs.')
deallocate(xyzzyaabm1,xyzzyaabk1,xyzzyaabl1,xyzzyaabn1,xyzzyaabo1,xyzz&
&yaabp1,xyzzyaabq1)
call finish_storage_energy
call finish_storage
call finish_energy_utils
call finish_wfn_utils
call finish_scratch
end subroutine xyzzyaacf1
subroutine xyzzyaacg1(is)
implicit none
integer,intent(in) :: is
integer xyzzyaaaa19,xyzzyaaab19,xyzzyaaac19,xyzzyaaad19,xyzzyaaae19,xy&
&zzyaaaf19,xyzzyaaag19,xyzzyaaah19
real(dp) v,xyzzyaaai19,xyzzyaaaj19(xyzzyaaab1),xyzzyaaak19,xyzzyaaal19&
&,xyzzyaaam19,xyzzyaaan19,ecomps(n_ecomp)
complex(dp) xyzzyaaao19(3),xyzzyaaap19
call eval_local_energy(is,ecomps=ecomps)
v=ecomps(i_pei)
xyzzyaabo1=0.d0
xyzzyaaal19=0.d0
do xyzzyaaaf19=1,netot
call wfn_loggrad(xyzzyaaaf19,is,1,xyzzyaaao19,prefetch_sd=.true.,prefe&
&tch_aderiv=.true.)
call wfn_loglap(xyzzyaaaf19,is,1,xyzzyaaap19)
xyzzyaabo1(1:3,xyzzyaaaf19)=dble(xyzzyaaao19(1:3))*sqrt(inv_pmass(whic&
&h_spin(xyzzyaaaf19)))
xyzzyaaal19=xyzzyaaal19+dble(xyzzyaaap19)*inv_pmass(which_spin(xyzzyaa&
&af19))
enddo
xyzzyaaan19=ddot(three_netot,xyzzyaabo1(1,1),1,xyzzyaabo1(1,1),1)
xyzzyaaal19=xyzzyaaal19+xyzzyaaan19
call get_linear_basis(is,xyzzyaaab1,xyzzyaaaa1,xyzzyaaad1,xyzzyaabm1,x&
&yzzyaaaj19,xyzzyaabn1,xyzzyaaak19)
call dcopy(three_netot,xyzzyaabn1(1),1,xyzzyaabp1(1),1)
call dadd(three_netot,xyzzyaabo1(1,1),1,xyzzyaabp1(1),1)
call dcopy(three_netot,xyzzyaabp1(1),1,xyzzyaabq1(1),1)
call dadd(three_netot,xyzzyaabo1(1,1),1,xyzzyaabq1(1),1)
xyzzyaaai19=xyzzyaaak19+xyzzyaaal19+ddot(three_netot,xyzzyaabn1(1),1,x&
&yzzyaabq1(1),1)
xyzzyaaac19=0
do xyzzyaaaa19=1,xyzzyaaab1
xyzzyaabk1(xyzzyaaaa19)=2*ddot(three_netot,xyzzyaabm1(1,xyzzyaaaa19),1&
&,xyzzyaabp1(1),1)+xyzzyaaaj19(xyzzyaaaa19)
do xyzzyaaab19=xyzzyaaaa19,xyzzyaaab1
xyzzyaaac19=xyzzyaaac19+1
xyzzyaabl1(xyzzyaaac19)=ddot(three_netot,xyzzyaabm1(1,xyzzyaaaa19),1,x&
&yzzyaabm1(1,xyzzyaaab19),1)
enddo
enddo
xyzzyaaam19=0.5d0*xyzzyaaai19-v
xyzzyaaay1=xyzzyaaay1+xyzzyaaam19**2
call daxpy(xyzzyaaab1,xyzzyaaam19,xyzzyaabk1(1),1,xyzzyaabc1(1),1)
xyzzyaaac19=1
xyzzyaaad19=1
xyzzyaaae19=1
xyzzyaaag19=xyzzyaaab1+1
do xyzzyaaaa19=1,xyzzyaaab1
xyzzyaaag19=xyzzyaaag19-1
call daxpy(xyzzyaaag19,0.25d0*xyzzyaabk1(xyzzyaaaa19),xyzzyaabk1(xyzzy&
&aaaa19),1,xyzzyaabd1(xyzzyaaac19),1)
call daxpy(xyzzyaaag19,xyzzyaaam19,xyzzyaabl1(xyzzyaaac19),1,xyzzyaabd&
&1(xyzzyaaac19),1)
xyzzyaaah19=((xyzzyaaag19+1)*xyzzyaaag19)/2+1
do xyzzyaaab19=xyzzyaaaa19,xyzzyaaab1
call daxpy(xyzzyaaab1,0.5d0*xyzzyaabl1(xyzzyaaac19),xyzzyaabk1(1),1,xy&
&zzyaabe1(xyzzyaaad19),1)
xyzzyaaad19=xyzzyaaad19+xyzzyaaab1
xyzzyaaah19=xyzzyaaah19-1
call daxpy(xyzzyaaah19,0.25d0*xyzzyaabl1(xyzzyaaac19),xyzzyaabl1(xyzzy&
&aaac19),1,xyzzyaabf1(xyzzyaaae19),1)
xyzzyaaae19=xyzzyaaae19+xyzzyaaah19
xyzzyaaac19=xyzzyaaac19+1
enddo
enddo
xyzzyaaaw1=xyzzyaaaw1+1.d0
xyzzyaaax1=xyzzyaaax1+xyzzyaaai19
call dadd(xyzzyaabg1,xyzzyaabk1(1),1,xyzzyaaba1(1),1)
call dadd(xyzzyaabh1,xyzzyaabl1(1),1,xyzzyaabb1(1),1)
xyzzyaaaz1=xyzzyaaaz1+v
end subroutine xyzzyaacg1
subroutine xyzzyaach1
implicit none
integer xyzzyaaaa20,xyzzyaaab20
real(dp) xyzzyaaac20,xyzzyaaad20,xyzzyaaae20,xyzzyaaaf20,xyzzyaaag20
real(dp),allocatable :: xyzzyaaah20(:)
if(nnodes>1)then
if(.not.am_master)allocate(xyzzyaaah20(1))
call mpi_reduce(xyzzyaaaw1,xyzzyaaac20,1,mpi_double_precision,mpi_sum,&
&0,mpi_comm_world,xyzzyaaaa20)
call checkmpi(xyzzyaaaa20,'Sum weight in average_G_data')
call mpi_reduce(xyzzyaaax1,xyzzyaaad20,1,mpi_double_precision,mpi_sum,&
&0,mpi_comm_world,xyzzyaaaa20)
call checkmpi(xyzzyaaaa20,'Sum g0 in average_G_data')
call mpi_reduce(xyzzyaaay1,xyzzyaaae20,1,mpi_double_precision,mpi_sum,&
&0,mpi_comm_world,xyzzyaaaa20)
call checkmpi(xyzzyaaaa20,'Sum BG0 in average_G_data')
call mpi_reduce(xyzzyaaaz1,xyzzyaaaf20,1,mpi_double_precision,mpi_sum,&
&0,mpi_comm_world,xyzzyaaaa20)
call checkmpi(xyzzyaaaa20,'Sum V in average_G_data')
if(am_master)then
xyzzyaaaw1=xyzzyaaac20
xyzzyaaax1=xyzzyaaad20
xyzzyaaay1=xyzzyaaae20
xyzzyaaaz1=xyzzyaaaf20
allocate(xyzzyaaah20(xyzzyaabg1),stat=xyzzyaaab20)
call check_alloc(xyzzyaaab20,'AVERAGE_G_DATA','1')
else
xyzzyaaaw1=0.d0
xyzzyaaax1=0.d0
xyzzyaaay1=0.d0
xyzzyaaaz1=0.d0
endif
call mpi_reduce(xyzzyaaba1,xyzzyaaah20,xyzzyaabg1,mpi_double_precision&
&,mpi_sum,0,mpi_comm_world,xyzzyaaaa20)
call checkmpi(xyzzyaaaa20,'Sum g1 in average_G_data')
if(am_master)then
call dcopy(xyzzyaabg1,xyzzyaaah20(1),1,xyzzyaaba1(1),1)
else
xyzzyaaba1=0.d0
endif
call mpi_reduce(xyzzyaabc1,xyzzyaaah20,xyzzyaabg1,mpi_double_precision&
&,mpi_sum,0,mpi_comm_world,xyzzyaaaa20)
call checkmpi(xyzzyaaaa20,'Sum BG1 in average_G_data')
if(am_master)then
call dcopy(xyzzyaabg1,xyzzyaaah20(1),1,xyzzyaabc1(1),1)
deallocate(xyzzyaaah20)
allocate(xyzzyaaah20(xyzzyaabh1),stat=xyzzyaaab20)
call check_alloc(xyzzyaaab20,'AVERAGE_G_DATA','2')
else
xyzzyaabc1=0.d0
endif
call mpi_reduce(xyzzyaabb1,xyzzyaaah20,xyzzyaabh1,mpi_double_precision&
&,mpi_sum,0,mpi_comm_world,xyzzyaaaa20)
call checkmpi(xyzzyaaaa20,'Sum g2 in average_G_data')
if(am_master)then
call dcopy(xyzzyaabh1,xyzzyaaah20(1),1,xyzzyaabb1(1),1)
else
xyzzyaabb1=0.d0
endif
call mpi_reduce(xyzzyaabd1,xyzzyaaah20,xyzzyaabh1,mpi_double_precision&
&,mpi_sum,0,mpi_comm_world,xyzzyaaaa20)
call checkmpi(xyzzyaaaa20,'Sum BG2 in average_G_data')
if(am_master)then
call dcopy(xyzzyaabh1,xyzzyaaah20(1),1,xyzzyaabd1(1),1)
deallocate(xyzzyaaah20)
allocate(xyzzyaaah20(xyzzyaabi1),stat=xyzzyaaab20)
call check_alloc(xyzzyaaab20,'AVERAGE_G_DATA','3')
else
xyzzyaabd1=0.d0
endif
call mpi_reduce(xyzzyaabe1,xyzzyaaah20,xyzzyaabi1,mpi_double_precision&
&,mpi_sum,0,mpi_comm_world,xyzzyaaaa20)
call checkmpi(xyzzyaaaa20,'Sum BG3 in average_G_data')
if(am_master)then
call dcopy(xyzzyaabi1,xyzzyaaah20(1),1,xyzzyaabe1(1),1)
deallocate(xyzzyaaah20)
allocate(xyzzyaaah20(xyzzyaabj1),stat=xyzzyaaab20)
call check_alloc(xyzzyaaab20,'AVERAGE_G_DATA','4')
else
xyzzyaabe1=0.d0
endif
call mpi_reduce(xyzzyaabf1,xyzzyaaah20,xyzzyaabj1,mpi_double_precision&
&,mpi_sum,0,mpi_comm_world,xyzzyaaaa20)
call checkmpi(xyzzyaaaa20,'Sum BG4 in average_G_data')
if(am_master)then
call dcopy(xyzzyaabj1,xyzzyaaah20(1),1,xyzzyaabf1(1),1)
deallocate(xyzzyaaah20)
else
xyzzyaabf1=0.d0
endif
if(.not.am_master)deallocate(xyzzyaaah20)
call mpi_bcast(xyzzyaaaw1,1,mpi_double_precision,0,mpi_comm_world,xyzz&
&yaaaa20)
call mpi_bcast(xyzzyaaax1,1,mpi_double_precision,0,mpi_comm_world,xyzz&
&yaaaa20)
call mpi_bcast(xyzzyaaba1,xyzzyaabg1,mpi_double_precision,0,mpi_comm_w&
&orld,xyzzyaaaa20)
call mpi_bcast(xyzzyaabb1,xyzzyaabh1,mpi_double_precision,0,mpi_comm_w&
&orld,xyzzyaaaa20)
call mpi_bcast(xyzzyaaay1,1,mpi_double_precision,0,mpi_comm_world,xyzz&
&yaaaa20)
call mpi_bcast(xyzzyaabc1,xyzzyaabg1,mpi_double_precision,0,mpi_comm_w&
&orld,xyzzyaaaa20)
call mpi_bcast(xyzzyaabd1,xyzzyaabh1,mpi_double_precision,0,mpi_comm_w&
&orld,xyzzyaaaa20)
call mpi_bcast(xyzzyaabe1,xyzzyaabi1,mpi_double_precision,0,mpi_comm_w&
&orld,xyzzyaaaa20)
call mpi_bcast(xyzzyaabf1,xyzzyaabj1,mpi_double_precision,0,mpi_comm_w&
&orld,xyzzyaaaa20)
call mpi_bcast(xyzzyaaaz1,1,mpi_double_precision,0,mpi_comm_world,xyzz&
&yaaaa20)
endif
if(xyzzyaaaw1>=0.5d0)then
xyzzyaaag20=1.d0/xyzzyaaaw1
xyzzyaaax1=xyzzyaaax1*xyzzyaaag20
call dscal(xyzzyaabg1,xyzzyaaag20,xyzzyaaba1,1)
call dscal(xyzzyaabh1,xyzzyaaag20,xyzzyaabb1,1)
xyzzyaaay1=xyzzyaaay1*xyzzyaaag20
call dscal(xyzzyaabg1,xyzzyaaag20,xyzzyaabc1,1)
call dscal(xyzzyaabh1,xyzzyaaag20,xyzzyaabd1,1)
call dscal(xyzzyaabi1,xyzzyaaag20,xyzzyaabe1,1)
call dscal(xyzzyaabj1,xyzzyaaag20,xyzzyaabf1,1)
xyzzyaaaz1=xyzzyaaaz1*xyzzyaaag20
endif
xyzzyaaaw1=1.d0
end subroutine xyzzyaach1
subroutine xyzzyaaci1(no_lines,alpha)
implicit none
integer,intent(in) :: no_lines
real(dp),intent(inout) :: alpha(xyzzyaaab1)
integer xyzzyaaaa21,xyzzyaaab21,xyzzyaaac21,xyzzyaaad21,xyzzyaaae21
real(dp) xyzzyaaaf21,xyzzyaaag21,xyzzyaaah21,xyzzyaaai21,xyzzyaaaj21,x&
&yzzyaaak21,xyzzyaaal21,xyzzyaaam21
real(dp),allocatable :: xyzzyaaan21(:),xyzzyaaao21(:)
allocate(xyzzyaaan21(xyzzyaaab1),stat=xyzzyaaab21)
call check_alloc(xyzzyaaab21,'MINIMIZE_LSF_MC','dir array')
do xyzzyaaaa21=1,no_lines
do xyzzyaaac21=1,xyzzyaaab1
xyzzyaaan21(xyzzyaaac21)=ranx()-0.5d0
enddo
call xyzzyaace1(alpha,xyzzyaaan21,xyzzyaaag21,xyzzyaaah21,xyzzyaaai21,&
&xyzzyaaaj21,xyzzyaaak21)
if(xyzzyaaak21>0.d0)then
xyzzyaaaf21=xyzzyaacd1(xyzzyaaah21,xyzzyaaai21,xyzzyaaaj21,xyzzyaaak21&
&)
else
xyzzyaaaf21=0.d0
endif
xyzzyaaam21=(((xyzzyaaak21*xyzzyaaaf21+xyzzyaaaj21)*xyzzyaaaf21+xyzzya&
&aai21)*xyzzyaaaf21+xyzzyaaah21)*xyzzyaaaf21+xyzzyaaag21
if(xyzzyaaam21<0.d0)cycle
call daxpy(xyzzyaaab1,xyzzyaaaf21,xyzzyaaan21(1),1,alpha(1),1)
enddo
deallocate(xyzzyaaan21)
xyzzyaaal21=lsfval(alpha)
if(nnodes>1)then
if(am_master)then
allocate(xyzzyaaao21(0:nnodes-1),stat=xyzzyaaab21)
else
allocate(xyzzyaaao21(1),stat=xyzzyaaab21)
endif
call check_alloc(xyzzyaaab21,'MINIMIZE_LSF_MC','')
call mpi_gather(xyzzyaaal21,1,mpi_double_precision,xyzzyaaao21,1,mpi_d&
&ouble_precision,0,mpi_comm_world,xyzzyaaae21)
call checkmpi(xyzzyaaae21,'Gather LSF values in minimize_LSF_MC.')
if(am_master)then
xyzzyaaad21=0
do xyzzyaaac21=1,nnodes-1
if(xyzzyaaao21(xyzzyaaac21)<xyzzyaaal21)then
xyzzyaaal21=xyzzyaaao21(xyzzyaaac21)
xyzzyaaad21=xyzzyaaac21
endif
enddo
endif
deallocate(xyzzyaaao21)
call mpi_bcast(xyzzyaaad21,1,mpi_integer,0,mpi_comm_world,xyzzyaaae21)
call checkmpi(xyzzyaaae21,'Broadcasting best_node in minimize_LSF_MC.'&
&)
call mpi_bcast(alpha,xyzzyaaab1,mpi_double_precision,xyzzyaaad21,mpi_c&
&omm_world,xyzzyaaae21)
call checkmpi(xyzzyaaae21,'Broadcasting alpha in minimize_LSF_MC.')
endif
end subroutine xyzzyaaci1
subroutine xyzzyaacj1(maxiter,alpha)
use toms573, only : nl2sol,dfault
implicit none
integer,intent(in) :: maxiter
real(dp),intent(inout) :: alpha(xyzzyaaab1)
integer xyzzyaaaa22,xyzzyaaab22
integer,allocatable :: xyzzyaaac22(:)
real(dp) xyzzyaaad22,xyzzyaaae22
real(dp),allocatable :: xyzzyaaaf22(:),xyzzyaaag22(:)
external calcr_varmin_linjas,calcj_varmin_linjas
if(am_master)then
allocate(xyzzyaaag22(xyzzyaaab1),stat=xyzzyaaaa22)
call check_alloc(xyzzyaaaa22,'MINIMIZE_LSF_GAUSSNEWTON','alpha_init')
xyzzyaaad22=lsfval(alpha)
call dcopy(xyzzyaaab1,alpha(1),1,xyzzyaaag22(1),1)
allocate(xyzzyaaac22(xyzzyaaab1+60),xyzzyaaaf22(93+xyzzyaaab1*(xyzzyaa&
&ab1+3)+xyzzyaaab1*(3*xyzzyaaab1+33)/2),stat=xyzzyaaaa22)
call check_alloc(xyzzyaaaa22,'MINIMIZE_LSF_GAUSSNEWTON','NL2SOL arrays&
&')
call dfault(xyzzyaaac22,xyzzyaaaf22)
xyzzyaaac22(14:15)=0
xyzzyaaac22(17)=500000
xyzzyaaac22(18)=maxiter
xyzzyaaac22(19:24)=0
call nl2sol(xyzzyaaab1,xyzzyaaab1,alpha,calcr_varmin_linjas,calcj_varm&
&in_linjas,xyzzyaaac22,xyzzyaaaf22)
if(xyzzyaaac22(1)==3)then
call wout('X-convergence.  Gauss-Newton process complete.')
elseif(xyzzyaaac22(1)==4)then
call wout('Relative-function convergence.  Gauss-Newton process comple&
&te.')
elseif(xyzzyaaac22(1)==5)then
call wout('Both X and relative-function convergence.  Gauss-Newton pro&
&cess complete.')
elseif(xyzzyaaac22(1)==6)then
call wout('Absolute function convergence.  Gauss-Newton process comple&
&te.')
elseif(xyzzyaaac22(1)==7)then
call wout('Singular convergence.  Gauss-Newton process complete.')
elseif(xyzzyaaac22(1)==8)then
call wout('False convergence.  Gauss-Newton process complete.')
elseif(xyzzyaaac22(1)==9.or.xyzzyaaac22(1)==10)then
call wout('Max. no. of iterations reached.  Gauss-Newton process compl&
&ete.')
else
call wout('NL2SOL return code: '//trim(i2s(xyzzyaaac22(1))))
call errwarn('MINIMIZE_LSF_GAUSSNEWTON','Gauss-Newton process was unsu&
&ccessful.')
endif
deallocate(xyzzyaaac22,xyzzyaaaf22)
xyzzyaaae22=lsfval(alpha)
if(xyzzyaaae22>xyzzyaaad22)call dcopy(xyzzyaaab1,xyzzyaaag22(1),1,alph&
&a(1),1)
deallocate(xyzzyaaag22)
endif
if(nnodes>1)then
call mpi_bcast(alpha,xyzzyaaab1,mpi_double_precision,0,mpi_comm_world,&
&xyzzyaaab22)
call checkmpi(xyzzyaaab22,'Broadcasting alpha in minimize_LSF_GaussNew&
&ton.')
endif
end subroutine xyzzyaacj1
subroutine xyzzyaack1(maxiter,alpha)
implicit none
integer,intent(in) :: maxiter
real(dp),intent(inout) :: alpha(xyzzyaaab1)
integer xyzzyaaaa23,xyzzyaaab23,xyzzyaaac23
real(dp) xyzzyaaad23,xyzzyaaae23,xyzzyaaaf23,xyzzyaaag23,xyzzyaaah23,x&
&yzzyaaai23,xyzzyaaaj23
real(dp),allocatable :: xyzzyaaak23(:)
if(am_master)then
allocate(xyzzyaaak23(xyzzyaaab1),stat=xyzzyaaaa23)
call check_alloc(xyzzyaaaa23,'MINIMIZE_LSF_SD','gradient')
do xyzzyaaab23=1,maxiter
call lsfgrad(alpha,xyzzyaaak23)
call xyzzyaace1(alpha,xyzzyaaak23,xyzzyaaae23,xyzzyaaaf23,xyzzyaaag23,&
&xyzzyaaah23,xyzzyaaai23)
if(xyzzyaaai23>0.d0)then
xyzzyaaad23=xyzzyaacd1(xyzzyaaaf23,xyzzyaaag23,xyzzyaaah23,xyzzyaaai23&
&)
else
exit
endif
if(xyzzyaaad23==0.d0)exit
xyzzyaaaj23=(((xyzzyaaai23*xyzzyaaad23+xyzzyaaah23)*xyzzyaaad23+xyzzya&
&aag23)*xyzzyaaad23+xyzzyaaaf23)*xyzzyaaad23+xyzzyaaae23
if(xyzzyaaaj23<0.d0)exit
call daxpy(xyzzyaaab1,xyzzyaaad23,xyzzyaaak23(1),1,alpha(1),1)
enddo
deallocate(xyzzyaaak23)
endif
if(nnodes>1)then
call mpi_bcast(alpha,xyzzyaaab1,mpi_double_precision,0,mpi_comm_world,&
&xyzzyaaac23)
call checkmpi(xyzzyaaac23,'Broadcasting alpha in minimize_LSF_SD.')
endif
end subroutine xyzzyaack1
subroutine xyzzyaacl1(maxiter,alpha)
implicit none
integer,intent(in) :: maxiter
real(dp),intent(inout) :: alpha(xyzzyaaab1)
integer xyzzyaaaa24,xyzzyaaab24,xyzzyaaac24,xyzzyaaad24
real(dp) xyzzyaaae24,xyzzyaaaf24,xyzzyaaag24,xyzzyaaah24,xyzzyaaai24,x&
&yzzyaaaj24,xyzzyaaak24,xyzzyaaal24,xyzzyaaam24,xyzzyaaan24,xyzzyaaao2&
&4
real(dp),allocatable :: xyzzyaaap24(:),xyzzyaaaq24(:),xyzzyaaar24(:)
logical xyzzyaaas24
if(am_master)then
allocate(xyzzyaaap24(xyzzyaaab1),xyzzyaaaq24(xyzzyaaab1),xyzzyaaar24(x&
&yzzyaaab1),stat=xyzzyaaaa24)
call check_alloc(xyzzyaaaa24,'MINIMIZE_LSF_CG','dir,...')
xyzzyaaao24=lsfval(alpha)
call lsfgrad(alpha,xyzzyaaap24)
do xyzzyaaac24=1,xyzzyaaab1
xyzzyaaaq24(xyzzyaaac24)=-xyzzyaaap24(xyzzyaaac24)
xyzzyaaar24(xyzzyaaac24)=xyzzyaaaq24(xyzzyaaac24)
xyzzyaaap24(xyzzyaaac24)=xyzzyaaar24(xyzzyaaac24)
enddo
xyzzyaaas24=.false.
do xyzzyaaab24=1,maxiter
call xyzzyaace1(alpha,xyzzyaaap24,xyzzyaaaf24,xyzzyaaag24,xyzzyaaah24,&
&xyzzyaaai24,xyzzyaaaj24)
if(xyzzyaaaj24>0.d0)then
xyzzyaaae24=xyzzyaacd1(xyzzyaaag24,xyzzyaaah24,xyzzyaaai24,xyzzyaaaj24&
&)
else
xyzzyaaae24=0.d0
endif
xyzzyaaan24=(((xyzzyaaaj24*xyzzyaaae24+xyzzyaaai24)*xyzzyaaae24+xyzzya&
&aah24)*xyzzyaaae24+xyzzyaaag24)*xyzzyaaae24+xyzzyaaaf24
if(xyzzyaaan24<0.d0)then
xyzzyaaas24=.true.
exit
endif
if(xyzzyaaan24>=xyzzyaaao24)then
if(xyzzyaaas24)exit
xyzzyaaas24=.true.
xyzzyaaao24=lsfval(alpha)
call lsfgrad(alpha,xyzzyaaap24)
do xyzzyaaac24=1,xyzzyaaab1
xyzzyaaaq24(xyzzyaaac24)=-xyzzyaaap24(xyzzyaaac24)
xyzzyaaar24(xyzzyaaac24)=xyzzyaaaq24(xyzzyaaac24)
xyzzyaaap24(xyzzyaaac24)=xyzzyaaar24(xyzzyaaac24)
enddo
cycle
endif
xyzzyaaao24=xyzzyaaan24
call daxpy(xyzzyaaab1,xyzzyaaae24,xyzzyaaap24(1),1,alpha(1),1)
call lsfgrad(alpha,xyzzyaaap24)
xyzzyaaal24=ddot(xyzzyaaab1,xyzzyaaaq24(1),1,xyzzyaaaq24(1),1)
xyzzyaaam24=0.d0
do xyzzyaaac24=1,xyzzyaaab1
xyzzyaaam24=xyzzyaaam24+(xyzzyaaap24(xyzzyaaac24)+xyzzyaaaq24(xyzzyaaa&
&c24))*xyzzyaaap24(xyzzyaaac24)
enddo
if(xyzzyaaal24==0.d0)exit
xyzzyaaak24=xyzzyaaam24/xyzzyaaal24
do xyzzyaaac24=1,xyzzyaaab1
xyzzyaaaq24(xyzzyaaac24)=-xyzzyaaap24(xyzzyaaac24)
xyzzyaaar24(xyzzyaaac24)=xyzzyaaaq24(xyzzyaaac24)+xyzzyaaak24*xyzzyaaa&
&r24(xyzzyaaac24)
xyzzyaaap24(xyzzyaaac24)=xyzzyaaar24(xyzzyaaac24)
enddo
enddo
deallocate(xyzzyaaap24,xyzzyaaaq24,xyzzyaaar24)
endif
if(nnodes>1)then
call mpi_bcast(alpha,xyzzyaaab1,mpi_double_precision,0,mpi_comm_world,&
&xyzzyaaad24)
call checkmpi(xyzzyaaad24,'Broadcasting alpha in minimize_LSF_CG.')
endif
end subroutine xyzzyaacl1
subroutine xyzzyaacm1(maxiter,alpha)
implicit none
integer,intent(in) :: maxiter
real(dp),intent(inout) :: alpha(xyzzyaaab1)
integer xyzzyaaaa25,xyzzyaaab25,xyzzyaaac25,xyzzyaaad25,xyzzyaaae25
real(dp) xyzzyaaaf25,xyzzyaaag25,xyzzyaaah25,xyzzyaaai25,xyzzyaaaj25,x&
&yzzyaaak25,xyzzyaaal25,xyzzyaaam25,xyzzyaaan25,xyzzyaaao25,xyzzyaaap2&
&5
real(dp),allocatable :: xyzzyaaaq25(:),xyzzyaaar25(:),xyzzyaaas25(:),x&
&yzzyaaat25(:),xyzzyaaau25(:,:)
logical xyzzyaaav25
if(am_master)then
allocate(xyzzyaaar25(xyzzyaaab1),xyzzyaaas25(xyzzyaaab1),xyzzyaaaq25(x&
&yzzyaaab1),xyzzyaaat25(xyzzyaaab1),xyzzyaaau25(xyzzyaaab1,xyzzyaaab1)&
&,stat=xyzzyaaaa25)
call check_alloc(xyzzyaaaa25,'MINIMIZE_LSF_CG','')
xyzzyaaal25=lsfval(alpha)
call lsfgrad(alpha,xyzzyaaar25)
do xyzzyaaab25=1,xyzzyaaab1
xyzzyaaau25(1:xyzzyaaab1,xyzzyaaab25)=0.d0
xyzzyaaau25(xyzzyaaab25,xyzzyaaab25)=1.d0
xyzzyaaaq25(xyzzyaaab25)=-xyzzyaaar25(xyzzyaaab25)
enddo
xyzzyaaav25=.false.
do xyzzyaaae25=1,maxiter
call xyzzyaace1(alpha,xyzzyaaaq25,xyzzyaaag25,xyzzyaaah25,xyzzyaaai25,&
&xyzzyaaaj25,xyzzyaaak25)
if(xyzzyaaak25>0.d0)then
xyzzyaaaf25=xyzzyaacd1(xyzzyaaah25,xyzzyaaai25,xyzzyaaaj25,xyzzyaaak25&
&)
else
xyzzyaaaf25=0.d0
endif
xyzzyaaap25=xyzzyaaal25
xyzzyaaal25=(((xyzzyaaak25*xyzzyaaaf25+xyzzyaaaj25)*xyzzyaaaf25+xyzzya&
&aai25)*xyzzyaaaf25+xyzzyaaah25)*xyzzyaaaf25+xyzzyaaag25
if(xyzzyaaal25<0.d0)then
xyzzyaaav25=.true.
exit
endif
if(xyzzyaaap25<=xyzzyaaal25)then
if(xyzzyaaav25)exit
xyzzyaaav25=.true.
do xyzzyaaab25=1,xyzzyaaab1
xyzzyaaau25(1:xyzzyaaab1,xyzzyaaab25)=0.d0
xyzzyaaau25(xyzzyaaab25,xyzzyaaab25)=1.d0
xyzzyaaaq25(xyzzyaaab25)=-xyzzyaaar25(xyzzyaaab25)
enddo
cycle
endif
call dscal(xyzzyaaab1,xyzzyaaaf25,xyzzyaaaq25(1),1)
call dadd(xyzzyaaab1,xyzzyaaaq25(1),1,alpha(1),1)
call dcopy(xyzzyaaab1,xyzzyaaar25(1),1,xyzzyaaas25(1),1)
call lsfgrad(alpha,xyzzyaaar25)
do xyzzyaaab25=1,xyzzyaaab1
xyzzyaaas25(xyzzyaaab25)=xyzzyaaar25(xyzzyaaab25)-xyzzyaaas25(xyzzyaaa&
&b25)
enddo
do xyzzyaaab25=1,xyzzyaaab1
xyzzyaaat25(xyzzyaaab25)=ddot(xyzzyaaab1,xyzzyaaau25(1,xyzzyaaab25),1,&
&xyzzyaaas25(1),1)
enddo
xyzzyaaam25=ddot(xyzzyaaab1,xyzzyaaas25(1),1,xyzzyaaaq25(1),1)
xyzzyaaan25=ddot(xyzzyaaab1,xyzzyaaas25(1),1,xyzzyaaat25(1),1)
if(xyzzyaaam25>0.d0.and.xyzzyaaan25/=0.d0)then
xyzzyaaam25=1.d0/xyzzyaaam25
xyzzyaaao25=1.d0/xyzzyaaan25
do xyzzyaaab25=1,xyzzyaaab1
xyzzyaaas25(xyzzyaaab25)=xyzzyaaam25*xyzzyaaaq25(xyzzyaaab25)-xyzzyaaa&
&o25*xyzzyaaat25(xyzzyaaab25)
enddo
do xyzzyaaab25=1,xyzzyaaab1
xyzzyaaau25(xyzzyaaab25,xyzzyaaab25)=xyzzyaaau25(xyzzyaaab25,xyzzyaaab&
&25)+xyzzyaaam25*xyzzyaaaq25(xyzzyaaab25)*xyzzyaaaq25(xyzzyaaab25)-xyz&
&zyaaao25*xyzzyaaat25(xyzzyaaab25)*xyzzyaaat25(xyzzyaaab25)+xyzzyaaan2&
&5*xyzzyaaas25(xyzzyaaab25)*xyzzyaaas25(xyzzyaaab25)
do xyzzyaaac25=xyzzyaaab25+1,xyzzyaaab1
xyzzyaaau25(xyzzyaaab25,xyzzyaaac25)=xyzzyaaau25(xyzzyaaab25,xyzzyaaac&
&25)+xyzzyaaam25*xyzzyaaaq25(xyzzyaaab25)*xyzzyaaaq25(xyzzyaaac25)-xyz&
&zyaaao25*xyzzyaaat25(xyzzyaaab25)*xyzzyaaat25(xyzzyaaac25)+xyzzyaaan2&
&5*xyzzyaaas25(xyzzyaaab25)*xyzzyaaas25(xyzzyaaac25)
xyzzyaaau25(xyzzyaaac25,xyzzyaaab25)=xyzzyaaau25(xyzzyaaab25,xyzzyaaac&
&25)
enddo
enddo
endif
do xyzzyaaab25=1,xyzzyaaab1
xyzzyaaaq25(xyzzyaaab25)=-ddot(xyzzyaaab1,xyzzyaaau25(1,xyzzyaaab25),1&
&,xyzzyaaar25(1),1)
enddo
enddo
deallocate(xyzzyaaar25,xyzzyaaas25,xyzzyaaaq25,xyzzyaaat25,xyzzyaaau25&
&)
endif
if(nnodes>1)then
call mpi_bcast(alpha,xyzzyaaab1,mpi_double_precision,0,mpi_comm_world,&
&xyzzyaaad25)
call checkmpi(xyzzyaaad25,'Broadcasting alpha in minimize_LSF_CG.')
endif
end subroutine xyzzyaacm1
subroutine xyzzyaacn1(maxiter,alpha)
implicit none
integer,intent(in) :: maxiter
real(dp),intent(inout) :: alpha(xyzzyaaab1)
integer xyzzyaaaa26,xyzzyaaab26,xyzzyaaac26,xyzzyaaad26
real(dp) xyzzyaaae26,xyzzyaaaf26,xyzzyaaag26,xyzzyaaah26,xyzzyaaai26,x&
&yzzyaaaj26,xyzzyaaak26,xyzzyaaal26
real(dp),allocatable :: xyzzyaaam26(:)
if(am_master)then
allocate(xyzzyaaam26(xyzzyaaab1),stat=xyzzyaaaa26)
call check_alloc(xyzzyaaaa26,'MINIMIZE_LSF_MC','dir array')
xyzzyaaam26(1:xyzzyaaab1)=0.d0
do xyzzyaaab26=1,maxiter
xyzzyaaaf26=0.d0
do xyzzyaaac26=1,xyzzyaaab1
xyzzyaaam26(xyzzyaaac26)=1.d0
call xyzzyaace1(alpha,xyzzyaaam26,xyzzyaaag26,xyzzyaaah26,xyzzyaaai26,&
&xyzzyaaaj26,xyzzyaaak26)
if(xyzzyaaak26>0.d0)then
xyzzyaaae26=xyzzyaacd1(xyzzyaaah26,xyzzyaaai26,xyzzyaaaj26,xyzzyaaak26&
&)
else
xyzzyaaae26=0.d0
endif
xyzzyaaal26=(((xyzzyaaak26*xyzzyaaae26+xyzzyaaaj26)*xyzzyaaae26+xyzzya&
&aai26)*xyzzyaaae26+xyzzyaaah26)*xyzzyaaae26+xyzzyaaag26
if(xyzzyaaal26<0.d0)exit
if(xyzzyaaae26>xyzzyaaaf26)xyzzyaaaf26=xyzzyaaae26
alpha(xyzzyaaac26)=alpha(xyzzyaaac26)+xyzzyaaae26
xyzzyaaam26(xyzzyaaac26)=0.d0
enddo
if(xyzzyaaaf26==0.d0)exit
enddo
deallocate(xyzzyaaam26)
endif
if(nnodes>1)then
call mpi_bcast(alpha,xyzzyaaab1,mpi_double_precision,0,mpi_comm_world,&
&xyzzyaaad26)
call checkmpi(xyzzyaaad26,'Broadcasting alpha in minimize_LSF_SD.')
endif
end subroutine xyzzyaacn1
end module slaarnacm_linjas
subroutine calcr_varmin_linjas(n,p,x,nf,r,uiparm,urparm,ufparm)
use dsp
use slaarnacm_linjas,only : lsfval
implicit none
integer,intent(in) :: n,p
integer,intent(inout) :: nf
integer,intent(inout),optional :: uiparm(:)
real(dp),intent(in) :: x(:)
real(dp),intent(out) :: r(:)
real(dp),intent(inout),optional :: urparm(:),ufparm
r(1)=lsfval(x)
r(2:n)=0.d0
end subroutine calcr_varmin_linjas
subroutine calcj_varmin_linjas(n,p,x,nf,j,uiparm,urparm,ufparm)
use dsp
use slaarnacm_linjas,only : lsfgrad
implicit none
integer,intent(in) :: n,p
integer,intent(inout) :: nf
integer,intent(inout),optional :: uiparm(:)
real(dp),intent(in) :: x(:)
real(dp),intent(inout),optional :: urparm(:),ufparm
real(dp),intent(out) :: j(:)
integer xyzzyaaaa28,xyzzyaaab28,xyzzyaaac28
real(dp) xyzzyaaad28(p)
call lsfgrad(x,xyzzyaaad28)
xyzzyaaac28=0
do xyzzyaaab28=1,p
xyzzyaaac28=xyzzyaaac28+1
j(xyzzyaaac28)=xyzzyaaad28(xyzzyaaab28)
do xyzzyaaaa28=2,n
xyzzyaaac28=xyzzyaaac28+1
j(xyzzyaaac28)=0.d0
enddo
enddo
end subroutine calcj_varmin_linjas

module slaarnaaf
use dsp
use parallel
use store
use file_utils,     only : open_units
use format_utils,   only : wout,i2s,int2log,log2int
use slaarnabt,      only : dcopy,zcopy
use slaarnacc, only : ranx_max,ranx_gauss_max,ranx_buffer,ranx_gauss_b&
&uffer,ranx_indx,ranx_gauss_indx
use run_control,    only : errstop,errstop_master,errwarn,timer,check_&
&alloc,get_total_time,total_time_previous
implicit none
private
public make_chkpoint_groups,load_configs,write_configs,init_config_acc&
&umulation,add_config,delete_config_vmc,end_config_accumulation,disman&
&tle_configs,shift_config_files,backup_config_file,copy_config_file,lo&
&gp_config
integer,parameter,private :: xyzzyaaaa1=1,xyzzyaaab1=2,xyzzyaaac1=3,xy&
&zzyaaad1=4,xyzzyaaae1=5,xyzzyaaaf1=6,xyzzyaaag1=7,xyzzyaaah1=8,xyzzya&
&aai1=9,xyzzyaaaj1=10,xyzzyaaak1=11,xyzzyaaal1=12,xyzzyaaam1=13,xyzzya&
&aan1=14,xyzzyaaao1=15,xyzzyaaap1=15
character(20),private :: xyzzyaaaq1(xyzzyaaap1)
data xyzzyaaaq1 /'RELE','SELE','ETOT','LOGDET','FIDET','FI_PROD_DET','&
&LAPDET','PROD_LAPDET','LOCAL_POTENTIAL','NLTOT','STOT','WDMC','VALJAS&
&','LOGP','TWIST'/
logical,private :: xyzzyaaar1=.false.,xyzzyaaas1=.false.,xyzzyaaat1(1:&
&xyzzyaaap1)=.false.
integer,private :: xyzzyaaau1=0,xyzzyaaav1=0,xyzzyaaaw1=0,xyzzyaaax1=0&
&,xyzzyaaay1=0
integer,allocatable,private :: xyzzyaaaz1(:)
character(20),allocatable,private :: xyzzyaaba1(:),xyzzyaabb1(:)
character(20),private :: xyzzyaabc1='NONE'
character(160),public :: con_in='config.in',con_out='config.out',con_b&
&ackup='config.backup'
character(146),public :: con_loc
integer,allocatable,public :: sele_config(:,:)
real(dp),allocatable,public :: rele_config(:,:,:),lapdet_config(:),eto&
&t_config(:),stot_config(:),local_potential_config(:),nltot_config(:),&
&fidet_config(:,:,:,:),fi_prod_det_config(:,:,:,:,:,:),prod_lapdet_con&
&fig(:,:,:,:,:),wdmc_config(:),valjas_config(:),twist_config(:,:)
real(dp),allocatable :: logp_config(:)
complex(dp),allocatable,public :: logdet_config(:,:,:)
integer,public :: vmc_steps_config,no_difftypes_config
integer,allocatable,public :: sele_vmc_config(:)
real(dp),public :: final_vmce_config,final_vmcde_config,final_vmcdeu_c&
&onfig,final_vmcvar_config,final_vmcmove_config,final_vmctau_config,fi&
&nal_vmce2_config,final_vmcde2_config,final_vmcde2u_config
real(dp),allocatable,public :: dtvmc_array_config(:),rele_vmc_config(:&
&,:)
integer xyzzyaabd1,xyzzyaabe1
integer,public :: dmcequil_steps_config,dmcstats_steps_config,tpdmc_co&
&nfig
real(dp),public :: ebest_config,ebest_init_config,eref_config,dteff_eb&
&est_init_config,dteff_best_config,numerator_wt2_config,denominator_wt&
&_config,denominator_wt2_config,log_pi_wt_config,log_pi_wt2_config,num&
&er_expect_config(13)
real(dp),allocatable,public :: log_pi_wt_array_config(:),log_pi_wt_arr&
&ay2_config(:)
logical,public :: lwdmc_config,growth_estimator_config
integer xyzzyaabf1
integer,public :: reblock_nbs_config
integer,allocatable,public :: reblock_block_length_config(:)
integer,allocatable,public :: reblock_s_u_f_config(:)
real(dp),allocatable,public :: reblock_sum_w2_closed_config(:)
real(dp),allocatable,public :: reblock_sum_w_open_config(:)
real(dp),allocatable,public :: reblock_sum_o2_closed_config(:,:)
real(dp),allocatable,public :: reblock_sum_o_open_config(:,:)
integer,public :: reblock_nstep_config,reblock_nobs_config
real(dp),public :: reblock_sum_w_config
real(dp),allocatable,public :: reblock_sum_o_config(:)
logical,public :: popstats_config
real(dp),public :: reblock_sum_w2_config
real(dp),allocatable,public :: reblock_sum_ow0_config(:)
real(dp),allocatable,public :: reblock_sum_ow2_config(:)
real(dp),allocatable,public :: reblock_sum_o2w0_config(:)
real(dp),allocatable,public :: reblock_sum_o2w2_config(:)
integer,public :: random_state_config(25)
integer,public :: nitot_config
real(dp),allocatable,public :: rion_config(:,:)
integer,public :: checkpoint_ncpu,mpicomm_group,ckmymaster
integer,allocatable :: xyzzyaabg1(:)
logical,public :: am_master_group
contains
subroutine make_chkpoint_groups
implicit none
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2,xyzzyaa&
&af2,xyzzyaaag2
integer,allocatable :: xyzzyaaah2(:),xyzzyaaai2(:),xyzzyaaaj2(:),xyzzy&
&aaak2(:),xyzzyaaal2(:)
logical xyzzyaaam2
if(checkpoint_ncpu<=0.or.checkpoint_ncpu>nnodes)checkpoint_ncpu=nnodes
xyzzyaaae2=nnodes/checkpoint_ncpu
allocate(xyzzyaaah2(xyzzyaaae2),xyzzyaaai2(xyzzyaaae2),xyzzyaaaj2(xyzz&
&yaaae2),xyzzyaaak2(xyzzyaaae2+1),xyzzyaaal2(xyzzyaaae2+1),stat=xyzzya&
&aag2)
call check_alloc(xyzzyaaag2,'MAKE_CHKPOINT_GROUPS','chkpoint_nnodes et&
&c.')
xyzzyaaaa2=(nnodes-xyzzyaaae2*checkpoint_ncpu)/xyzzyaaae2
xyzzyaaah2(:)=checkpoint_ncpu+xyzzyaaaa2
xyzzyaaab2=mod(nnodes-xyzzyaaae2*checkpoint_ncpu,xyzzyaaae2)
do xyzzyaaac2=1,xyzzyaaab2
xyzzyaaah2(xyzzyaaac2)=xyzzyaaah2(xyzzyaaac2)+1
enddo
xyzzyaaai2(1)=0
xyzzyaaaj2(1)=xyzzyaaah2(1)-1
xyzzyaaad2=1
do xyzzyaaaa2=2,xyzzyaaae2
xyzzyaaai2(xyzzyaaaa2)=xyzzyaaaj2(xyzzyaaaa2-1)+1
xyzzyaaaj2(xyzzyaaaa2)=xyzzyaaai2(xyzzyaaaa2)+xyzzyaaah2(xyzzyaaaa2)-1
if(my_node>=xyzzyaaai2(xyzzyaaaa2))xyzzyaaad2=xyzzyaaaa2
enddo
if(any(my_node==xyzzyaaai2(:)))then
xyzzyaaam2=.true.
ckmymaster=my_node
else
xyzzyaaam2=.false.
ckmymaster=xyzzyaaai2(xyzzyaaad2)
endif
call mpi_comm_group(mpi_comm_world,xyzzyaaaf2,ierror)
call checkmpi(ierror,'mpi_comm_group in setup_redist_groups')
do xyzzyaaaa2=1,xyzzyaaae2
allocate(xyzzyaabg1(xyzzyaaah2(xyzzyaaaa2)),stat=xyzzyaaag2)
call check_alloc(xyzzyaaag2,'MAKE_CHKPOINT_GROUPS','chkpoint_ranks')
xyzzyaaac2=0
do xyzzyaaab2=xyzzyaaai2(xyzzyaaaa2),xyzzyaaaj2(xyzzyaaaa2)
xyzzyaaac2=xyzzyaaac2+1
xyzzyaabg1(xyzzyaaac2)=xyzzyaaab2
enddo
call mpi_group_incl(xyzzyaaaf2,xyzzyaaah2(xyzzyaaaa2),xyzzyaabg1,xyzzy&
&aaak2(xyzzyaaaa2),ierror)
call checkmpi(ierror,'mpi_group_incl in make_checkpoint_groups')
call mpi_comm_create(mpi_comm_world,xyzzyaaak2(xyzzyaaaa2),xyzzyaaal2(&
&xyzzyaaaa2),ierror)
call checkmpi(ierror,'mpi_comm_create in make_checkpoint_groups')
deallocate(xyzzyaabg1)
enddo
checkpoint_ncpu=xyzzyaaah2(xyzzyaaad2)
mpicomm_group=xyzzyaaal2(xyzzyaaad2)
am_master_group=xyzzyaaam2
if(am_master_group)then
allocate(xyzzyaabg1(checkpoint_ncpu-1),stat=xyzzyaaag2)
xyzzyaaac2=0
do xyzzyaaab2=xyzzyaaai2(xyzzyaaad2)+1,xyzzyaaaj2(xyzzyaaad2)
xyzzyaaac2=xyzzyaaac2+1
xyzzyaabg1(xyzzyaaac2)=xyzzyaaab2
enddo
endif
deallocate(xyzzyaaal2,xyzzyaaak2,xyzzyaaah2,xyzzyaaai2,xyzzyaaaj2)
end subroutine make_chkpoint_groups
subroutine load_configs(nconfig_inout,gen_by,required_item,optional_it&
&em,present_citem,required_extra,optional_extra,present_eitem,from_bac&
&kup)
implicit none
integer,intent(inout) :: nconfig_inout
logical,optional,intent(in) :: from_backup
logical,dimension(:),intent(out) :: present_citem,present_eitem
character(*),intent(in) :: gen_by
character(20),dimension(:),intent(in) :: required_item,optional_item,r&
&equired_extra,optional_extra
character(20) gen_by_expected
character(160) con_file
call timer('LOAD_CONFIGS',.true.)
gen_by_expected=gen_by
con_file=con_in
if(present(from_backup))then
if(from_backup)con_file=con_backup
endif
if(xyzzyaaas1)call errstop_master('LOAD_CONFIGS','Have not called END_&
&CONFIG_ACCUMULATION.')
if(.not.xyzzyaaar1)call xyzzyaabh1(con_file,nconfig_inout)
call xyzzyaabi1(nconfig_inout,gen_by_expected,required_item,optional_i&
&tem,present_citem,required_extra,optional_extra,present_eitem)
call timer('LOAD_CONFIGS',.false.)
end subroutine load_configs
subroutine xyzzyaabh1(con_file,nconfig_inout)
use slaarnabg, only : atom_basis_type
implicit none
integer,intent(in) :: nconfig_inout
character(160),intent(in) :: con_file
integer xyzzyaaaa4,xyzzyaaab4,xyzzyaaac4,xyzzyaaad4,xyzzyaaae4,xyzzyaa&
&af4,xyzzyaaag4,xyzzyaaah4,xyzzyaaai4,xyzzyaaaj4,xyzzyaaak4,xyzzyaaal4
logical xyzzyaaam4
character(20) label
integer xyzzyaaan4,xyzzyaaao4,xyzzyaaap4,xyzzyaaaq4,xyzzyaaar4,xyzzyaa&
&as4
integer,allocatable :: xyzzyaaat4(:)
character(20) interaction_in,atom_basis_type_in
integer xyzzyaaau4(10),xyzzyaaav4,xyzzyaaaw4,xyzzyaaax4,xyzzyaaay4(25)&
&,xyzzyaaaz4,xyzzyaaba4
real xyzzyaabb4(ranx_max),xyzzyaabc4(ranx_gauss_max)
real(dp) xyzzyaabd4(23)
logical xyzzyaabe4
integer,allocatable :: xyzzyaabf4(:)
real(dp),allocatable :: xyzzyaabg4(:,:)
integer,allocatable :: xyzzyaabh4(:,:)
real(dp),allocatable :: xyzzyaabi4(:,:,:),xyzzyaabj4(:),xyzzyaabk4(:),&
&xyzzyaabl4(:),xyzzyaabm4(:),xyzzyaabn4(:),xyzzyaabo4(:,:,:,:),xyzzyaa&
&bp4(:,:,:,:,:,:),xyzzyaabq4(:,:,:,:,:),xyzzyaabr4(:),xyzzyaabs4(:),xy&
&zzyaabt4(:),xyzzyaabu4(:,:)
complex(dp),allocatable :: xyzzyaabv4(:,:,:)
xyzzyaabc1='NONE'
xyzzyaaaw1=0
xyzzyaaay1=0
xyzzyaaap4=0
xyzzyaaan4=-1
xyzzyaaao4=-1
ebest_config=0.d0
ebest_init_config=huge(1.d0)
eref_config=0.d0
dteff_ebest_init_config=0.d0
dteff_best_config=0.d0
vmc_steps_config=0
dmcequil_steps_config=0
dmcstats_steps_config=0
xyzzyaaar4=0
xyzzyaaas4=0
numer_expect_config(1:13)=0.d0
denominator_wt_config=0.d0
log_pi_wt_config=0.d0
log_pi_wt2_config=0.d0
numerator_wt2_config=0.d0
denominator_wt2_config=0.d0
tpdmc_config=0
lwdmc_config=.false.
growth_estimator_config=.false.
if(am_master_group)then
inquire(file=trim(con_file),exist=xyzzyaaam4)
call mpi_bcast(xyzzyaaam4,1,mpi_logical,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast con_file in read_configs')
if(.not.xyzzyaaam4)return
call open_units(xyzzyaaab4,xyzzyaaac4)
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Could not find free i/o &
&unit number on node '//trim(i2s(my_node))//'.')
open(unit=xyzzyaaab4,file=trim(con_file),form='unformatted',status='ol&
&d',iostat=xyzzyaaac4)
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem opening config f&
&ile on node '//trim(i2s(my_node))//'.')
read(xyzzyaaab4,iostat=xyzzyaaac4)label
if(xyzzyaaac4<0)call errstop('READ_CONFIGS','Config file appears to be&
& empty?')
if(xyzzyaaac4>0)call errstop('READ_CONFIGS','Problem reading first lab&
&el in config file. Perhaps this config file is in an older format? Tr&
&y using the update_config utility.')
if(trim(label)/='INFO')call errstop('READ_CONFIGS','First label in con&
&fig file is wrong. Perhaps this config file is in old format? Try usi&
&ng the update_config utility.')
rewind(xyzzyaaab4)
do
read(xyzzyaaab4,iostat=xyzzyaaac4)label
call mpi_bcast(xyzzyaaac4,1,mpi_integer,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast ierr in read_configs [master]')
if(xyzzyaaac4<0)then
if(nconfig_inout<=0)then
close(xyzzyaaab4)
open_unit(xyzzyaaab4)=.false.
xyzzyaaau1=0
xyzzyaaav1=0
xyzzyaaar1=.true.
if(xyzzyaaay1>0)then
call mpi_bcast(xyzzyaabb1(1)(1:1),len(xyzzyaabb1)*xyzzyaaay1,mpi_chara&
&cter,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast extra_item in read_configs [master]')
endif
return
endif
call errstop('READ_CONFIGS','End-of-file reached while reading section&
& name.')
endif
if(xyzzyaaac4>0)call errstop('READ_CONFIGS','Problem reading section n&
&ame.')
call mpi_bcast(label,len(label),mpi_character,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast label in read_configs [master]')
select case(trim(label))
case('INFO')
do
read(xyzzyaaab4,iostat=xyzzyaaac4)label
if(xyzzyaaac4<0)call errstop('READ_CONFIGS','End-of-file reached while&
& reading section names.')
if(xyzzyaaac4>0)call errstop('READ_CONFIGS','Problem reading section n&
&ame.')
if(trim(label)=='END INFO')goto 10
select case(trim(label))
case('NEXTRA')
read(xyzzyaaab4,iostat=xyzzyaaac4)xyzzyaaay1
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading NEXTRA.'&
&)
allocate(xyzzyaabb1(xyzzyaaay1),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','extra_item')
case('GEN_BY')
read(xyzzyaaab4,iostat=xyzzyaaac4)xyzzyaabc1
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading GEN_BY.'&
&)
case('NSPIN')
read(xyzzyaaab4,iostat=xyzzyaaac4)xyzzyaaan4
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading NSPIN.')
if(xyzzyaaan4/=nspin)call errstop('READ_CONFIGS','Wrong NSPIN in confi&
&g file.')
case('NELE')
if(xyzzyaaan4<0)call errstop('READ_CONFIGS','NELE found before NSPIN i&
&n config file. This is a bug.')
if(allocated(xyzzyaaat4))call errstop('READ_CONFIGS','NELE found twice&
& in config file.')
allocate(xyzzyaaat4(nspin),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','NELE_IN')
read(xyzzyaaab4,iostat=xyzzyaaac4)xyzzyaaat4
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading NELE.')
if(any(xyzzyaaat4/=nele))call errstop('READ_CONFIGS','Wrong NELE in co&
&nfig file.')
case('NDET')
read(xyzzyaaab4,iostat=xyzzyaaac4)xyzzyaaao4
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading NDET.')
if(xyzzyaaao4/=ndet)call errstop('READ_CONFIGS','Wrong NDET in config &
&file.')
case('COMPLEX_WF')
read(xyzzyaaab4,iostat=xyzzyaaac4)xyzzyaaar4
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading COMPLEX_&
&WF.')
if(int2log(xyzzyaaar4).neqv.complex_wf)call errstop('READ_CONFIGS','Wr&
&ong COMPLEX_WF in config file.')
case('NONCOLL_SPIN')
read(xyzzyaaab4,iostat=xyzzyaaac4)xyzzyaaas4
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading NONCOLL_&
&SPIN.')
if(int2log(xyzzyaaas4).neqv.noncoll_spin)call errstop('READ_CONFIGS','&
&Wrong NONCOLL_SPIN in config file.')
case('INTERACTION')
read(xyzzyaaab4,iostat=xyzzyaaac4)interaction_in
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading INTERACT&
&ION.')
if(trim(interaction_in)/=trim(interaction))call errwarn('READ_CONFIGS'&
&,'INTERACTION in config file differs from current value.')
case('ATOM_BASIS_TYPE')
read(xyzzyaaab4,iostat=xyzzyaaac4)atom_basis_type_in
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading ATOM_BAS&
&IS_TYPE.')
if(trim(atom_basis_type_in)/=trim(atom_basis_type))call errstop('READ_&
&CONFIGS','Wrong ATOM_BASIS_TYPE in config file.')
case('TOTAL_TIME')
read(xyzzyaaab4,iostat=xyzzyaaac4)total_time_previous
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading TOTAL_TI&
&ME.')
case default
call errstop('READ_CONFIGS',"Unknown info '"//trim(label)//"' found in&
& config file.")
end select
enddo
10 continue
call mpi_bcast(xyzzyaaay1,1,mpi_integer,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast nextra in read_configs [master]')
call mpi_bcast(xyzzyaabc1,20,mpi_character,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast generated_by in read_configs [master]')
case('DEFINE_CONFIGS')
read(xyzzyaaab4,iostat=xyzzyaaac4)xyzzyaaav1
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading number o&
&f configs.')
if(xyzzyaaav1<1)call errstop('READ_CONFIGS','Config file defines confi&
&gs but contains none?')
read(xyzzyaaab4,iostat=xyzzyaaac4)xyzzyaaaw1
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading number o&
&f items per config.')
allocate(xyzzyaaba1(xyzzyaaaw1),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','config_item')
xyzzyaaat1=.false.
do xyzzyaaaf4=1,xyzzyaaaw1
read(xyzzyaaab4,iostat=xyzzyaaac4)label
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading config c&
&ontent label.')
xyzzyaaba1(xyzzyaaaf4)=trim(label)
xyzzyaaai4=0
do xyzzyaaaj4=1,xyzzyaaap1
if(trim(xyzzyaaba1(xyzzyaaaf4))==trim(xyzzyaaaq1(xyzzyaaaj4)))then
xyzzyaaai4=xyzzyaaaj4
goto 20
endif
enddo
20 continue
if(xyzzyaaai4==0)call errstop('READ_CONFIGS',"Unrecognized item '"//tr&
&im(label)//"'.")
if(any(xyzzyaaat1(xyzzyaaai4:)))call errstop('READ_CONFIGS','Config it&
&ems in wrong order or repeated.')
xyzzyaaat1(xyzzyaaai4)=.true.
enddo
read(xyzzyaaab4,iostat=xyzzyaaac4)label
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS',"Problem reading 'END DEF&
&INE_CONFIGS'.")
if(trim(label)/='END DEFINE_CONFIGS')call errstop('READ_CONFIGS',"Expe&
&cted to find 'END DEFINE_CONFIGS' but didn't.")
call mpi_bcast(xyzzyaaav1,1,mpi_integer,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast nconfig_total in read_configs [master]')
call mpi_bcast(xyzzyaaat1,xyzzyaaap1,mpi_logical,0,mpicomm_group,ierro&
&r)
call checkmpi(ierror,'bcast lconfig_item in read_configs [master]')
case('VMC_SAVED_STATE')
xyzzyaaap4=xyzzyaaap4+1
if(xyzzyaaap4>xyzzyaaay1)call errstop('READ_CONFIGS','Extra block VMC_&
&SAVED_STATE found in config file but not expected.')
xyzzyaabb1(xyzzyaaap4)=trim(label)
do
read(xyzzyaaab4,iostat=xyzzyaaac4)label
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS',"Problem reading label in&
& VMC_SAVED_STATE section in config file.")
call mpi_bcast(label,len(label),mpi_character,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast label in read_configs [master]')
select case(trim(label))
case('NO_DIFFTYPES')
read(xyzzyaaab4,iostat=xyzzyaaac4)no_difftypes_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading NO_DIFFT&
&YPES.')
if(no_difftypes_config/=no_difftypes)call errstop('READ_CONFIGS','Wron&
&g NO_DIFFTYPES in config file.')
case('DTVMC')
allocate(dtvmc_array_config(no_difftypes),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','DTVMC_ARRAY_CONFIG')
read(xyzzyaaab4,iostat=xyzzyaaac4)dtvmc_array_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading DTVMC.')
call mpi_bcast(dtvmc_array_config,no_difftypes,mpi_double_precision,0,&
&mpicomm_group,ierror)
call checkmpi(ierror,'bcast dtvmc_array_config in read_configs [master&
&]')
case('VMC_STEPS')
read(xyzzyaaab4,iostat=xyzzyaaac4)vmc_steps_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading VMC_STEP&
&S.')
call mpi_bcast(vmc_steps_config,1,mpi_integer,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast vmc_steps_config in read_configs [master]'&
&)
case('NNODES')
read(xyzzyaaab4,iostat=xyzzyaaac4)xyzzyaaaq4
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading NNODES.'&
&)
call mpi_bcast(xyzzyaaaq4,1,mpi_integer,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast nnodes_in in read_configs [master]')
case('RELE')
allocate(rele_vmc_config(3,netot),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','RELE_VMC_CONFIG')
rele_vmc_config=0.d0
if(checkpoint_ncpu>1)then
allocate(xyzzyaabg4(3,netot),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','RELE_VMC_BUFF')
xyzzyaabg4=0.d0
endif
xyzzyaaav4=0
do xyzzyaaah4=0,xyzzyaaaq4-1
if(xyzzyaaah4==my_node)then
read(xyzzyaaab4,iostat=xyzzyaaac4)rele_vmc_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading RELE on &
&node '//trim(i2s(my_node))//'.')
elseif(xyzzyaaah4>my_node.and.xyzzyaaah4<(my_node+checkpoint_ncpu))the&
&n
read(xyzzyaaab4,iostat=xyzzyaaac4)xyzzyaabg4
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading RELE for&
& other nodes on node '//trim(i2s(my_node))//'.')
xyzzyaaav4=xyzzyaaav4+1
call qmpi_ssend(xyzzyaabg4,xyzzyaabg1(xyzzyaaav4),1,'READ_CONFIGS','RE&
&LE_VMC_BUFF')
call checkmpi(ierror,'ssend rele_vmc_buff in read_configs')
else
read(xyzzyaaab4,iostat=xyzzyaaac4)
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem skipping over ot&
&her nodes'' RELE on node '//trim(i2s(my_node))//'.')
endif
enddo
if(checkpoint_ncpu>1)deallocate(xyzzyaabg4)
case('SELE')
allocate(sele_vmc_config(netot),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','SELE_VMC_CONFIG')
sele_vmc_config(:)=0
if(checkpoint_ncpu>1)then
allocate(xyzzyaabf4(netot),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','SELE_VMC_BUFF')
xyzzyaabf4=0
endif
xyzzyaaav4=0
do xyzzyaaah4=0,xyzzyaaaq4-1
if(xyzzyaaah4==my_node)then
read(xyzzyaaab4,iostat=xyzzyaaac4)sele_vmc_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading SELE on &
&node '//trim(i2s(my_node))//'.')
elseif(xyzzyaaah4>my_node.and.xyzzyaaah4<(my_node+checkpoint_ncpu))the&
&n
read(xyzzyaaab4,iostat=xyzzyaaac4)xyzzyaabf4
xyzzyaaav4=xyzzyaaav4+1
call qmpi_ssend(xyzzyaabf4,xyzzyaabg1(xyzzyaaav4),1,'READ_CONFIGS','SE&
&LE_VMC_BUFF')
call checkmpi(ierror,'ssend sele_vmc_buff in read_configs')
else
read(xyzzyaaab4,iostat=xyzzyaaac4)
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem skipping over ot&
&her nodes'' SELE on node '//trim(i2s(my_node))//'.')
endif
enddo
if(checkpoint_ncpu>1)deallocate(xyzzyaabf4)
case('END VMC_SAVED_STATE')
goto 30
case default
call errstop('READ_CONFIGS','Unrecognized information '//trim(label)//&
&' in VMC_SAVED_STATE section of config file.')
end select
enddo
30 continue
case('FINAL_VMC_RESULT')
xyzzyaaap4=xyzzyaaap4+1
if(xyzzyaaap4>xyzzyaaay1)call errstop('READ_CONFIGS','Extra block FINA&
&L_VMC_RESULT found in config file but not expected.')
xyzzyaabb1(xyzzyaaap4)=trim(label)
final_vmce_config=0.d0
final_vmcde_config=0.d0
final_vmcdeu_config=0.d0
final_vmcvar_config=0.d0
final_vmctau_config=0.d0
final_vmcmove_config=0.d0
final_vmce2_config=0.d0
final_vmcde2_config=0.d0
final_vmcde2u_config=0.d0
do
read(xyzzyaaab4,iostat=xyzzyaaac4)label
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS',"Problem reading label in&
& FINAL_VMC_RESULT section in config file.")
select case(trim(label))
case('ENERGY')
read(xyzzyaaab4,iostat=xyzzyaaac4)final_vmce_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading ENERGY')
case('ERRORBAR')
read(xyzzyaaab4,iostat=xyzzyaaac4)final_vmcde_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading ERRORBAR&
&')
case('ERRORBARU')
read(xyzzyaaab4,iostat=xyzzyaaac4)final_vmcdeu_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading ERRORBAR&
&U')
case('VARIANCE')
read(xyzzyaaab4,iostat=xyzzyaaac4)final_vmcvar_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading VARIANCE&
&')
case('CORRTIME')
read(xyzzyaaab4,iostat=xyzzyaaac4)final_vmctau_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading CORRTIME&
&')
case('TOTMOVE')
read(xyzzyaaab4,iostat=xyzzyaaac4)final_vmcmove_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading TOTMOVE'&
&)
case('ENERGY2')
read(xyzzyaaab4,iostat=xyzzyaaac4)final_vmce2_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading ENERGY2'&
&)
case('ERRORBAR2')
read(xyzzyaaab4,iostat=xyzzyaaac4)final_vmcde2_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading ERRORBAR&
&2')
case('ERRORBAR2U')
read(xyzzyaaab4,iostat=xyzzyaaac4)final_vmcde2u_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading ERRORBAR&
&2U')
case('END FINAL_VMC_RESULT')
goto 40
case default
call errstop('READ_CONFIGS','Unrecognized information '//trim(label)//&
&' in FINAL_VMC_RESULT section of config file.')
end select
enddo
40 continue
xyzzyaabd4(1)=final_vmce_config
xyzzyaabd4(2)=final_vmcde_config
xyzzyaabd4(3)=final_vmcdeu_config
xyzzyaabd4(4)=final_vmcvar_config
xyzzyaabd4(5)=final_vmctau_config
xyzzyaabd4(6)=final_vmcmove_config
xyzzyaabd4(7)=final_vmce2_config
xyzzyaabd4(8)=final_vmcde2_config
xyzzyaabd4(9)=final_vmcde2u_config
call mpi_bcast(xyzzyaabd4,9,mpi_double_precision,0,mpicomm_group,ierro&
&r)
call checkmpi(ierror,'bcast vmc_saved_state''s dbuff in read_configs [&
&master].')
case('DMC_SAVED_STATE')
xyzzyaaap4=xyzzyaaap4+1
if(xyzzyaaap4>xyzzyaaay1)call errstop('READ_CONFIGS','Extra block DMC_&
&SAVED_STATE found in config file but not expected.')
xyzzyaabb1(xyzzyaaap4)=trim(label)
xyzzyaaau4(4)=0
xyzzyaaau4(5)=0
do
read(xyzzyaaab4,iostat=xyzzyaaac4)label
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS',"Problem reading label in&
& DMC_SAVED_STATE section in config file.")
select case(trim(label))
case('EBEST')
read(xyzzyaaab4,iostat=xyzzyaaac4)ebest_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading EBEST.')
case('EBEST_INIT')
read(xyzzyaaab4,iostat=xyzzyaaac4)ebest_init_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading EBEST_IN&
&IT.')
case('EREF')
read(xyzzyaaab4,iostat=xyzzyaaac4)eref_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading EREF.')
case('DTEFF_EBEST_INIT')
read(xyzzyaaab4,iostat=xyzzyaaac4)dteff_ebest_init_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading DTEFF_EB&
&EST_INIT.')
case('DTEFF_BEST')
read(xyzzyaaab4,iostat=xyzzyaaac4)dteff_best_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading DTEFF_BE&
&ST.')
case('DMCEQUIL_STEPS')
read(xyzzyaaab4,iostat=xyzzyaaac4)dmcequil_steps_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading DMCEQUIL&
&_STEPS.')
case('DMCSTATS_STEPS')
read(xyzzyaaab4,iostat=xyzzyaaac4)dmcstats_steps_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading DMCSTATS&
&_STEPS.')
case('TPDMC')
read(xyzzyaaab4,iostat=xyzzyaaac4)tpdmc_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading TPDMC.')
case('LWDMC')
read(xyzzyaaab4,iostat=xyzzyaaac4)xyzzyaabd1
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading LWDMC.')
lwdmc_config=int2log(xyzzyaabd1)
case('GROWTH_ESTIMATOR')
read(xyzzyaaab4,iostat=xyzzyaaac4)xyzzyaabe1
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading GROWTH_E&
&STIMATOR.')
growth_estimator_config=int2log(xyzzyaabe1)
case('NUMER_EXPECT')
read(xyzzyaaab4,iostat=xyzzyaaac4)numer_expect_config(1:13)
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading NUMER_EX&
&PECT.')
case('DENOMINATOR_WT')
read(xyzzyaaab4,iostat=xyzzyaaac4)denominator_wt_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading DENOMINA&
&TOR_WT.')
case('LOG_PI_WT')
read(xyzzyaaab4,iostat=xyzzyaaac4)log_pi_wt_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading LOG_PI_W&
&T.')
case('LOG_PI_WT2')
read(xyzzyaaab4,iostat=xyzzyaaac4)log_pi_wt2_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading LOG_PI_W&
&T2.')
case('LOG_PI_WT_ARRAY')
xyzzyaaau4(4)=1
if(tpdmc_config<0)call errstop('READ_CONFIGS','LOG_PI_WT_ARRAY found b&
&efore TPDMC in config file. This is a bug.')
if(tpdmc_config==0)call errstop('READ_CONFIGS','LOG_PI_WT_ARRAY found &
&in config file even though the stored TPDMC is zero. This is a bug.')
if(allocated(log_pi_wt_array_config))deallocate(log_pi_wt_array_config&
&)
allocate(log_pi_wt_array_config(0:tpdmc_config-1),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','LOG_PI_WT_ARRAY')
read(xyzzyaaab4,iostat=xyzzyaaac4)log_pi_wt_array_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading LOG_PI_W&
&T_ARRAY.')
case('LOG_PI_WT_ARRAY2')
xyzzyaaau4(5)=1
if(tpdmc_config<0)call errstop('READ_CONFIGS','LOG_PI_WT_ARRAY2 found &
&before TPDMC in config file. This is a bug.')
if(allocated(log_pi_wt_array2_config))deallocate(log_pi_wt_array2_conf&
&ig)
allocate(log_pi_wt_array2_config(0:tpdmc_config),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','LOG_PI_WT_ARRAY2')
read(xyzzyaaab4,iostat=xyzzyaaac4)log_pi_wt_array2_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading LOG_PI_W&
&T_ARRAY2.')
case('NUMERATOR_WT2')
read(xyzzyaaab4,iostat=xyzzyaaac4)numerator_wt2_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading NUMERATO&
&R_WT2.')
case('DENOMINATOR_WT2')
read(xyzzyaaab4,iostat=xyzzyaaac4)denominator_wt2_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading DENOMINA&
&TOR_WT2.')
case('END DMC_SAVED_STATE')
goto 50
case default
call errstop('READ_CONFIGS','Unrecognized information '//trim(label)//&
&' in DMC_SAVED_STATE section of config file.')
end select
enddo
50 continue
if(ebest_init_config>=huge(1.d0))ebest_init_config=ebest_config
xyzzyaabd4(1)=ebest_config
xyzzyaabd4(2)=eref_config
xyzzyaabd4(3)=dteff_ebest_init_config
xyzzyaabd4(4)=dteff_best_config
xyzzyaabd4(5:17)=numer_expect_config(1:13)
xyzzyaabd4(18)=denominator_wt_config
xyzzyaabd4(19)=log_pi_wt_config
xyzzyaabd4(20)=log_pi_wt2_config
xyzzyaabd4(21)=numerator_wt2_config
xyzzyaabd4(22)=denominator_wt2_config
xyzzyaabd4(23)=ebest_init_config
call mpi_bcast(xyzzyaabd4,23,mpi_double_precision,0,mpicomm_group,ierr&
&or)
call checkmpi(ierror,'bcast dmc_saved_state''s dbuff in read_configs [&
&master].')
xyzzyaaau4(1)=dmcequil_steps_config
xyzzyaaau4(2)=dmcstats_steps_config
xyzzyaaau4(3)=tpdmc_config
xyzzyaaau4(6)=0
if(lwdmc_config)xyzzyaaau4(6)=1
xyzzyaaau4(7)=0
if (growth_estimator_config)xyzzyaaau4(7)=1
call mpi_bcast(xyzzyaaau4,7,mpi_integer,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast dmc_saved_state''s ibuff in read_configs [&
&master]')
if(xyzzyaaau4(4)/=0)then
call mpi_bcast(log_pi_wt_array_config(0),tpdmc_config,mpi_double_preci&
&sion,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast log_pi_wt_array_config in read_configs [ma&
&ster]')
endif
if(xyzzyaaau4(5)/=0)then
call mpi_bcast(log_pi_wt_array2_config(0),tpdmc_config+1,mpi_double_pr&
&ecision,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast log_pi_wt_array2_config in read_configs [m&
&aster]')
endif
case('REBLOCK_DATA')
xyzzyaaap4=xyzzyaaap4+1
if(xyzzyaaap4>xyzzyaaay1)call errstop('READ_CONFIGS','Extra block REBL&
&OCK_DATA found in config file but not expected.')
xyzzyaabb1(xyzzyaaap4)=trim(label)
read(xyzzyaaab4,iostat=xyzzyaaac4)reblock_nbs_config,reblock_nstep_con&
&fig,reblock_nobs_config,reblock_sum_w_config,xyzzyaabf1
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading # of reb&
&lock blocks.')
popstats_config=int2log(xyzzyaabf1)
if(allocated(reblock_block_length_config))deallocate(reblock_block_len&
&gth_config,reblock_s_u_f_config,reblock_sum_w2_closed_config,reblock_&
&sum_w_open_config,reblock_sum_o2_closed_config,reblock_sum_o_open_con&
&fig,reblock_sum_o_config)
allocate(reblock_block_length_config(reblock_nbs_config),reblock_s_u_f&
&_config(reblock_nbs_config),reblock_sum_w2_closed_config(reblock_nbs_&
&config),reblock_sum_w_open_config(reblock_nbs_config),reblock_sum_o2_&
&closed_config(reblock_nobs_config,reblock_nbs_config),reblock_sum_o_o&
&pen_config(reblock_nobs_config,reblock_nbs_config),reblock_sum_o_conf&
&ig(reblock_nobs_config),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','reblock_block_length_confi&
&g')
read(xyzzyaaab4,iostat=xyzzyaaac4)reblock_block_length_config,reblock_&
&s_u_f_config,reblock_sum_w2_closed_config,reblock_sum_w_open_config,r&
&eblock_sum_o2_closed_config,reblock_sum_o_open_config,reblock_sum_o_c&
&onfig
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading reblock &
&blocks.')
if(popstats_config)then
if(allocated(reblock_sum_ow0_config))deallocate(reblock_sum_ow0_config&
&,reblock_sum_ow2_config,reblock_sum_o2w0_config,reblock_sum_o2w2_conf&
&ig)
allocate(reblock_sum_ow0_config(reblock_nobs_config),reblock_sum_ow2_c&
&onfig(reblock_nobs_config),reblock_sum_o2w0_config(reblock_nobs_confi&
&g),reblock_sum_o2w2_config(reblock_nobs_config),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','reblock_sum_ow0_config')
read(xyzzyaaab4,iostat=xyzzyaaac4)reblock_sum_w2_config,reblock_sum_ow&
&0_config,reblock_sum_ow2_config,reblock_sum_o2w0_config,reblock_sum_o&
&2w2_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading popstats&
&.')
endif
read(xyzzyaaab4,iostat=xyzzyaaac4)label
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS',"Problem reading 'END REB&
&LOCK_DATA' label in config file.")
if(trim(label)/='END REBLOCK_DATA')call errstop('READ_CONFIGS',"Wrong &
&label read in when expecting 'END REBLOCK_DATA'.")
case('RANDOM')
xyzzyaaap4=xyzzyaaap4+1
if(xyzzyaaap4>xyzzyaaay1)call errstop('READ_CONFIGS','Extra block RAND&
&OM found in config file but not expected.')
xyzzyaabb1(xyzzyaaap4)=trim(label)
read(xyzzyaaab4,iostat=xyzzyaaac4)label
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading label in&
& RANDOM section in config file.')
if(trim(label)/='NNODES')call errstop('READ_CONFIGS','Wrong label for &
&NNODES in RANDOM section.')
read(xyzzyaaab4,iostat=xyzzyaaac4)xyzzyaaaq4
call mpi_bcast(xyzzyaaaq4,1,mpi_integer,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast nnodes_in in read_configs [master]')
random_state_config=0
if((checkpoint_ncpu>1).and.(my_node<xyzzyaaaq4-1))xyzzyaaay4=0
xyzzyaaav4=0
do xyzzyaaah4=0,xyzzyaaaq4-1
if(xyzzyaaah4==my_node)then
read(xyzzyaaab4,iostat=xyzzyaaac4)random_state_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading RANDOM_S&
&TATE on node '//trim(i2s(my_node))//'.')
elseif(xyzzyaaah4>my_node.and.xyzzyaaah4<(my_node+checkpoint_ncpu))the&
&n
xyzzyaaav4=xyzzyaaav4+1
read(xyzzyaaab4,iostat=xyzzyaaac4)xyzzyaaay4
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading RANDOM_S&
&TATE on node '//trim(i2s(my_node))//'.')
call qmpi_ssend(xyzzyaaay4,xyzzyaabg1(xyzzyaaav4),1,'READ_CONFIGS','RA&
&NDOM_STATE_BUFF')
call checkmpi(ierror,'ssend random_state_buff in read_configs on node &
&'//trim(i2s(my_node))//'.')
else
read(xyzzyaaab4,iostat=xyzzyaaac4)
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS',"Problem skipping over ot&
&her nodes' RANDOM_STATE on node "//trim(i2s(my_node))//'.')
endif
enddo
read(xyzzyaaab4,iostat=xyzzyaaac4)label
if(trim(label)=='BUFFER')then
xyzzyaabe4=.true.
else
xyzzyaabe4=.false.
endif
call mpi_bcast(xyzzyaabe4,1,mpi_logical,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast found_rng_buffer in read_configs [master]'&
&)
if(xyzzyaabe4)then
read(xyzzyaaab4)xyzzyaaak4,xyzzyaaal4
if(xyzzyaaak4/=ranx_max.and.xyzzyaaal4/=ranx_gauss_max)call errstop('R&
&EAD_CONFIGS','Config file was generated with random-number buffers of&
& a different size to those currently defined internally in CASINO.')
xyzzyaaav4=0
do xyzzyaaah4=0,xyzzyaaaq4-1
if(xyzzyaaah4==my_node)then
read(xyzzyaaab4,iostat=xyzzyaaac4)ranx_indx,ranx_gauss_indx
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading RANX_IND&
&X on node '//trim(i2s(my_node))//'.')
elseif(xyzzyaaah4>my_node.and.xyzzyaaah4<(my_node+checkpoint_ncpu))the&
&n
xyzzyaaav4=xyzzyaaav4+1
read(xyzzyaaab4,iostat=xyzzyaaac4)xyzzyaaaz4,xyzzyaaba4
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading RANX_IND&
&X_BUFF on node '//trim(i2s(my_node))//'.')
call qmpi_ssend(xyzzyaaaz4,xyzzyaabg1(xyzzyaaav4),1,'READ_CONFIGS','RA&
&NX_INDX_BUFF')
call checkmpi(ierror,'ssend ranx_indx_buff in read_configs on node '//&
&trim(i2s(my_node))//'.')
call qmpi_ssend(xyzzyaaba4,xyzzyaabg1(xyzzyaaav4),1,'READ_CONFIGS','RA&
&NX_GAUSS_INDX_BUFF')
call checkmpi(ierror,'ssend ranx_gauss_indx_buff in read_configs on no&
&de '//trim(i2s(my_node))//'.')
else
read(xyzzyaaab4,iostat=xyzzyaaac4)
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS',"Problem skipping over ot&
&her nodes' RANX_INDX on node "//trim(i2s(my_node))//'.')
endif
enddo
xyzzyaaav4=0
do xyzzyaaah4=0,xyzzyaaaq4-1
if(xyzzyaaah4==my_node)then
read(xyzzyaaab4,iostat=xyzzyaaac4)ranx_buffer
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading RANX_BUF&
&FER on node '//trim(i2s(my_node))//'.')
elseif(xyzzyaaah4>my_node.and.xyzzyaaah4<(my_node+checkpoint_ncpu))the&
&n
xyzzyaaav4=xyzzyaaav4+1
read(xyzzyaaab4,iostat=xyzzyaaac4)xyzzyaabb4
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading RANX_BUF&
&FER_BUFF on node '//trim(i2s(my_node))//'.')
call qmpi_ssend(xyzzyaabb4,xyzzyaabg1(xyzzyaaav4),1,'READ_CONFIGS','RA&
&NX_BUFFER_BUFF')
call checkmpi(ierror,'ssend ranx_buffer_buff in read_configs on node '&
&//trim(i2s(my_node))//'.')
else
read(xyzzyaaab4,iostat=xyzzyaaac4)
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS',"Problem skipping over ot&
&her nodes' RANX_BUFFER on node "//trim(i2s(my_node))//'.')
endif
enddo
xyzzyaaav4=0
do xyzzyaaah4=0,xyzzyaaaq4-1
if(xyzzyaaah4==my_node)then
read(xyzzyaaab4,iostat=xyzzyaaac4)ranx_gauss_buffer
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading RANX_GAU&
&SS_BUFFER on node '//trim(i2s(my_node))//'.')
elseif(xyzzyaaah4>my_node.and.xyzzyaaah4<(my_node+checkpoint_ncpu))the&
&n
xyzzyaaav4=xyzzyaaav4+1
read(xyzzyaaab4,iostat=xyzzyaaac4)xyzzyaabc4
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading RANX_GAU&
&SS_BUFFER_BUFF on node '//trim(i2s(my_node))//'.')
call qmpi_ssend(xyzzyaabc4,xyzzyaabg1(xyzzyaaav4),1,'READ_CONFIGS','RA&
&NX_GAUSS_BUFFER_BUFF')
call checkmpi(ierror,'ssend ranx_gauss_buffer_buff in read_configs on &
&node '//trim(i2s(my_node))//'.')
else
read(xyzzyaaab4,iostat=xyzzyaaac4)
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS',"Problem skipping over ot&
&her nodes' RANX_GAUSS_BUFFER on node "//trim(i2s(my_node))//'.')
endif
enddo
read(xyzzyaaab4,iostat=xyzzyaaac4)label
endif
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS',"Problem reading 'END RAN&
&DOM' label in config file.")
if(trim(label)/='END RANDOM')call errstop('READ_CONFIGS',"Wrong label &
&read in when expecting 'END RANDOM'.")
case('GEOMETRY')
xyzzyaaap4=xyzzyaaap4+1
if(xyzzyaaap4>xyzzyaaay1)call errstop('READ_CONFIGS','Extra block GEOM&
&ETRY found in config file but not expected.')
xyzzyaabb1(xyzzyaaap4)=trim(label)
read(xyzzyaaab4,iostat=xyzzyaaac4)nitot_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading NITOT.')
call mpi_bcast(nitot_config,1,mpi_integer,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast nitot_config in read_configs [master]')
if(allocated(rion_config))deallocate(rion_config)
if(nitot_config>0)then
allocate(rion_config(3,nitot_config),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','RION_CONFIG')
read(xyzzyaaab4,iostat=xyzzyaaac4)rion_config
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS','Problem reading RION.')
call mpi_bcast(rion_config,3*nitot_config,mpi_double_precision,0,mpico&
&mm_group,ierror)
call checkmpi(ierror,'bcast rion_config in read_configs [master]')
endif
read(xyzzyaaab4,iostat=xyzzyaaac4)label
if(xyzzyaaac4/=0)call errstop('READ_CONFIGS',"Problem reading 'END GEO&
&METRY' label in config file.")
if(trim(label)/='END GEOMETRY')call errstop('READ_CONFIGS',"Wrong labe&
&l read in when expecting 'END GEOMETRY'.")
case('CONFIGS')
goto 3
case default
call errstop('READ_CONFIGS','Unrecognized section '//trim(label)//' in&
& config file.')
end select
enddo
3 continue
call mpi_bcast(xyzzyaabb1(1)(1:1),len(xyzzyaabb1)*xyzzyaaay1,mpi_chara&
&cter,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast extra_item in read_configs [master]')
else
call mpi_bcast(xyzzyaaam4,1,mpi_logical,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast con_exists in read_configs [slave]')
if(.not.xyzzyaaam4)return
do
call mpi_bcast(xyzzyaaac4,1,mpi_integer,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast ierr in read_configs [slave]')
if(xyzzyaaac4<0)then
if(nconfig_inout<=0)then
xyzzyaaau1=0
xyzzyaaav1=0
xyzzyaaar1=.true.
if(xyzzyaaay1>0)then
call mpi_bcast(xyzzyaabb1(1)(1:1),len(xyzzyaabb1)*xyzzyaaay1,mpi_chara&
&cter,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast extra_item in read_configs [slave]')
endif
return
endif
endif
call mpi_bcast(label,len(label),mpi_character,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast label in read_configs [slave]')
select case(trim(label))
case("INFO")
call mpi_bcast(xyzzyaaay1,1,mpi_integer,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast nextra in read_configs [slave]')
call mpi_bcast(xyzzyaabc1,20,mpi_character,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast generated_by in read_configs [slave]')
allocate(xyzzyaabb1(xyzzyaaay1),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','extra_item')
case('DEFINE_CONFIGS')
call mpi_bcast(xyzzyaaav1,1,mpi_integer,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast nconfig_total in read_configs [slave]')
call mpi_bcast(xyzzyaaat1,xyzzyaaap1,mpi_logical,0,mpicomm_group,ierro&
&r)
call checkmpi(ierror,'bcast lconfig_item in read_configs [master]')
xyzzyaaaw1=count(xyzzyaaat1)
allocate(xyzzyaaba1(xyzzyaaaw1),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','config_item')
xyzzyaaba1=''
xyzzyaaaf4=0
do xyzzyaaaj4=1,xyzzyaaap1
if(xyzzyaaat1(xyzzyaaaj4))then
xyzzyaaaf4=xyzzyaaaf4+1
xyzzyaaba1(xyzzyaaaf4)=trim(xyzzyaaaq1(xyzzyaaaj4))
endif
enddo
case('VMC_SAVED_STATE')
do
call mpi_bcast(label,len(label),mpi_character,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast label in read_configs [slave]')
select case(trim(label))
case('DTVMC')
allocate(dtvmc_array_config(no_difftypes),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','DTVMC_ARRAY_CONFIG')
call mpi_bcast(dtvmc_array_config,no_difftypes,mpi_double_precision,0,&
&mpicomm_group,ierror)
call checkmpi(ierror,'bcast dtvmc_array_config in read_configs [slave]&
&')
case('VMC_STEPS')
call mpi_bcast(vmc_steps_config,1,mpi_integer,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast vmc_steps_config in read_configs [slave]')
case('NNODES')
call mpi_bcast(xyzzyaaaq4,1,mpi_integer,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast nnodes_in in read_configs [slave]')
case('RELE')
allocate(rele_vmc_config(3,netot),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','RELE_VMC_CONFIG')
rele_vmc_config=0.d0
if(my_node<xyzzyaaaq4)then
call qmpi_recv(rele_vmc_config,ckmymaster,1,'READ_CONFIGS','RELE_VMC_C&
&ONFIG')
endif
case('SELE')
allocate(sele_vmc_config(netot),stat=xyzzyaaad4)
sele_vmc_config(:)=0
call check_alloc(xyzzyaaad4,'READ_CONFIGS','SELE_VMC_CONFIG')
if(my_node<xyzzyaaaq4)then
call qmpi_recv(sele_vmc_config,ckmymaster,1,'READ_CONFIGS','SELE_VMC_C&
&ONFIG')
endif
case('END VMC_SAVED_STATE')
goto 60
end select
enddo
60  continue
case('FINAL_VMC_RESULT')
call mpi_bcast(xyzzyaabd4,9,mpi_double_precision,0,mpicomm_group,ierro&
&r)
call checkmpi(ierror,'bcast final_vmc_result''s dbuff in read_configs &
&[slave]')
final_vmce_config=xyzzyaabd4(1)
final_vmcde_config=xyzzyaabd4(2)
final_vmcdeu_config=xyzzyaabd4(3)
final_vmcvar_config=xyzzyaabd4(4)
final_vmctau_config=xyzzyaabd4(5)
final_vmcmove_config=xyzzyaabd4(6)
final_vmce2_config=xyzzyaabd4(7)
final_vmcde2_config=xyzzyaabd4(8)
final_vmcde2u_config=xyzzyaabd4(9)
case('DMC_SAVED_STATE')
call mpi_bcast(xyzzyaabd4,23,mpi_double_precision,0,mpicomm_group,ierr&
&or)
call checkmpi(ierror,'bcast dmc_saved_state''s dbuff in read_configs [&
&slave]')
ebest_config=xyzzyaabd4(1)
eref_config=xyzzyaabd4(2)
dteff_ebest_init_config=xyzzyaabd4(3)
dteff_best_config=xyzzyaabd4(4)
numer_expect_config(1:13)=xyzzyaabd4(5:17)
denominator_wt_config=xyzzyaabd4(18)
log_pi_wt_config=xyzzyaabd4(19)
log_pi_wt2_config=xyzzyaabd4(20)
numerator_wt2_config=xyzzyaabd4(21)
denominator_wt2_config=xyzzyaabd4(22)
ebest_init_config=xyzzyaabd4(23)
call mpi_bcast(xyzzyaaau4,7,mpi_integer,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast dmc_saved_state''s ibuff in read_configs [&
&slave]')
dmcequil_steps_config=xyzzyaaau4(1)
dmcstats_steps_config=xyzzyaaau4(2)
tpdmc_config=xyzzyaaau4(3)
lwdmc_config=.false.
if(xyzzyaaau4(6)/=0)lwdmc_config=.true.
growth_estimator_config=.false.
if(xyzzyaaau4(7)/=0)growth_estimator_config=.true.
if(xyzzyaaau4(4)/=0)then
if(allocated(log_pi_wt_array_config))deallocate(log_pi_wt_array_config&
&)
allocate(log_pi_wt_array_config(0:tpdmc_config-1),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','LOG_PI_WT_ARRAY')
call mpi_bcast(log_pi_wt_array_config(0),tpdmc_config,mpi_double_preci&
&sion,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast log_pi_wt_array_config in read_configs [sl&
&ave]')
endif
if(xyzzyaaau4(5)/=0)then
if(allocated(log_pi_wt_array2_config))deallocate(log_pi_wt_array2_conf&
&ig)
allocate(log_pi_wt_array2_config(0:tpdmc_config),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','LOG_PI_WT_ARRAY2')
call mpi_bcast(log_pi_wt_array2_config(0),tpdmc_config+1,mpi_double_pr&
&ecision,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast log_pi_wt_array2_config in read_configs [s&
&lave]')
endif
case('REBLOCK_DATA')
continue
case('RANDOM')
call mpi_bcast(xyzzyaaaq4,1,mpi_integer,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast nnodes_in in read_configs [slave]')
random_state_config=0
if(my_node<xyzzyaaaq4)then
call qmpi_recv(random_state_config,ckmymaster,1,'READ_CONFIGS','RANDOM&
&_STATE_CONFIG')
endif
call mpi_bcast(xyzzyaabe4,1,mpi_logical,0,mpicomm_group,ierror)
if(my_node<xyzzyaaaq4)then
if(xyzzyaabe4)then
call qmpi_recv(ranx_indx,ckmymaster,1,'READ_CONFIGS','RANX_INDX')
call qmpi_recv(ranx_gauss_indx,ckmymaster,1,'READ_CONFIGS','RANX_GAUSS&
&_INDX')
call qmpi_recv(ranx_buffer,ckmymaster,1,'READ_CONFIGS','RANX_BUFFER')
call qmpi_recv(ranx_gauss_buffer,ckmymaster,1,'READ_CONFIGS','RANX_GAU&
&SS_BUFFER')
endif
endif
case('GEOMETRY')
call mpi_bcast(nitot_config,1,mpi_integer,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast nitot_config in read_configs [slave]')
if(allocated(rion_config))deallocate(rion_config)
if(nitot_config>0)then
allocate(rion_config(3,nitot_config),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','RION_CONFIG')
call mpi_bcast(rion_config,3*nitot_config,mpi_double_precision,0,mpico&
&mm_group,ierror)
call checkmpi(ierror,'bcast rion_config in read_configs [slave]')
endif
case('CONFIGS')
goto 70
end select
enddo
70 continue
call mpi_bcast(xyzzyaabb1(1)(1:1),len(xyzzyaabb1)*xyzzyaaay1,mpi_chara&
&cter,0,mpicomm_group,ierror)
call checkmpi(ierror,'bcast extra_item in read_configs [slave]')
endif
if(am_master_group)then
if(xyzzyaaap4/=xyzzyaaay1)call errstop('READ_CONFIGS','Have read '//tr&
&im(i2s(xyzzyaaap4))//' extra items, but the config file claims it has&
& '//trim(i2s(xyzzyaaay1))//'.')
endif
if(xyzzyaaaw1==0)call errstop('READ_CONFIGS','Config contents undefine&
&d when about to read configs.')
allocate(xyzzyaaaz1(0:nnodes-1),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','nconfig_node')
xyzzyaaaz1(:)=xyzzyaaav1/nnodes
if(nconfig_inout<-1)then
xyzzyaaaa4=mod(xyzzyaaav1,nnodes)-1
xyzzyaaaz1(0:xyzzyaaaa4)=xyzzyaaaz1(0:xyzzyaaaa4)+1
elseif(nconfig_inout>0)then
xyzzyaaaz1(:)=min(xyzzyaaaz1(:),nconfig_inout)
endif
xyzzyaaau1=xyzzyaaaz1(my_node)
xyzzyaaav1=sum(xyzzyaaaz1)
if(xyzzyaaav1>0)then
if(am_master_group)then
xyzzyaaae4=0
if(my_node>0)xyzzyaaae4=sum(xyzzyaaaz1(0:my_node-1))*xyzzyaaaw1
do xyzzyaaaa4=1,xyzzyaaae4
read(xyzzyaaab4,iostat=xyzzyaaac4)
if(xyzzyaaac4<0)call errstop('READ_CONFIGS','End-of-file found when sk&
&ipping over other nodes data on node '//trim(i2s(my_node))//'.')
if(xyzzyaaac4>0)call errstop('READ_CONFIGS','Problem skipping over oth&
&er nodes data on node '//trim(i2s(my_node))//'.')
enddo
endif
if(xyzzyaaat1(xyzzyaaaa1))then
allocate(rele_config(3,netot,xyzzyaaau1),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','rele')
rele_config=0.d0
endif
if(xyzzyaaat1(xyzzyaaab1))then
allocate(sele_config(netot,xyzzyaaau1),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','sele')
sele_config=0
endif
if(xyzzyaaat1(xyzzyaaac1))then
allocate(etot_config(xyzzyaaau1),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','etot')
etot_config=0.d0
endif
if(xyzzyaaat1(xyzzyaaad1))then
allocate(logdet_config(nspin,ndet,xyzzyaaau1),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','logdet')
logdet_config=cmplx(0.d0,0.d0,dp)
endif
if(xyzzyaaat1(xyzzyaaae1))then
allocate(fidet_config(3,netot,real1_complex2,xyzzyaaau1),stat=xyzzyaaa&
&d4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','fidet')
fidet_config=0.d0
endif
if(xyzzyaaat1(xyzzyaaaf1))then
allocate(fi_prod_det_config(3,ndet,nemax,real1_complex2,nspin,xyzzyaaa&
&u1),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','fi_prod_det')
fi_prod_det_config=0.d0
endif
if(xyzzyaaat1(xyzzyaaag1))then
allocate(lapdet_config(xyzzyaaau1),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','lapdet')
lapdet_config=0.d0
endif
if(xyzzyaaat1(xyzzyaaah1))then
allocate(prod_lapdet_config(ndet,nemax,real1_complex2,nspin,xyzzyaaau1&
&),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','prod_lapdet')
prod_lapdet_config=0.d0
endif
if(xyzzyaaat1(xyzzyaaai1))then
allocate(local_potential_config(xyzzyaaau1),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','local_potential')
local_potential_config=0.d0
endif
if(xyzzyaaat1(xyzzyaaaj1))then
allocate(nltot_config(xyzzyaaau1),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','nltot')
nltot_config=0.d0
endif
if(xyzzyaaat1(xyzzyaaak1))then
allocate(stot_config(xyzzyaaau1),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','stot')
stot_config=0.d0
endif
if(xyzzyaaat1(xyzzyaaal1))then
allocate(wdmc_config(xyzzyaaau1),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','wdmc')
wdmc_config=0.d0
endif
if(xyzzyaaat1(xyzzyaaam1))then
allocate(valjas_config(xyzzyaaau1),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','valjas')
valjas_config=0.d0
endif
if(xyzzyaaat1(xyzzyaaan1))then
allocate(logp_config(xyzzyaaau1),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','logp')
logp_config=0.d0
endif
if(xyzzyaaat1(xyzzyaaao1))then
allocate(twist_config(5,xyzzyaaau1),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','twist')
twist_config=0.d0
endif
if(am_master_group)then
do xyzzyaaag4=1,xyzzyaaau1
if(xyzzyaaat1(xyzzyaaaa1))read(xyzzyaaab4,err=1,end=2)rele_config(:,:,&
&xyzzyaaag4)
if(xyzzyaaat1(xyzzyaaab1))read(xyzzyaaab4,err=1,end=2)sele_config(:,xy&
&zzyaaag4)
if(xyzzyaaat1(xyzzyaaac1))read(xyzzyaaab4,err=1,end=2)etot_config(xyzz&
&yaaag4)
if(xyzzyaaat1(xyzzyaaad1))read(xyzzyaaab4,err=1,end=2)logdet_config(:,&
&:,xyzzyaaag4)
if(xyzzyaaat1(xyzzyaaae1))read(xyzzyaaab4,err=1,end=2)fidet_config(:,:&
&,:,xyzzyaaag4)
if(xyzzyaaat1(xyzzyaaaf1))read(xyzzyaaab4,err=1,end=2)fi_prod_det_conf&
&ig(:,:,:,:,:,xyzzyaaag4)
if(xyzzyaaat1(xyzzyaaag1))read(xyzzyaaab4,err=1,end=2)lapdet_config(xy&
&zzyaaag4)
if(xyzzyaaat1(xyzzyaaah1))read(xyzzyaaab4,err=1,end=2)prod_lapdet_conf&
&ig(:,:,:,:,xyzzyaaag4)
if(xyzzyaaat1(xyzzyaaai1))read(xyzzyaaab4,err=1,end=2)local_potential_&
&config(xyzzyaaag4)
if(xyzzyaaat1(xyzzyaaaj1))read(xyzzyaaab4,err=1,end=2)nltot_config(xyz&
&zyaaag4)
if(xyzzyaaat1(xyzzyaaak1))read(xyzzyaaab4,err=1,end=2)stot_config(xyzz&
&yaaag4)
if(xyzzyaaat1(xyzzyaaal1))read(xyzzyaaab4,err=1,end=2)wdmc_config(xyzz&
&yaaag4)
if(xyzzyaaat1(xyzzyaaam1))read(xyzzyaaab4,err=1,end=2)valjas_config(xy&
&zzyaaag4)
if(xyzzyaaat1(xyzzyaaan1))read(xyzzyaaab4,err=1,end=2)logp_config(xyzz&
&yaaag4)
if(xyzzyaaat1(xyzzyaaao1))read(xyzzyaaab4,err=1,end=2)twist_config(:,x&
&yzzyaaag4)
enddo
if(checkpoint_ncpu>1)then
xyzzyaaax4=maxval(xyzzyaaaz1(my_node+1:my_node+checkpoint_ncpu-1))
if(xyzzyaaat1(xyzzyaaaa1))then
allocate(xyzzyaabi4(3,netot,xyzzyaaax4),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','rele_buff')
xyzzyaabi4=0.d0
endif
if(xyzzyaaat1(xyzzyaaab1))then
allocate(xyzzyaabh4(netot,xyzzyaaax4),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','sele_buff')
xyzzyaabh4=0
endif
if(xyzzyaaat1(xyzzyaaac1))then
allocate(xyzzyaabk4(xyzzyaaax4),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','etot_buff')
xyzzyaabk4=0.d0
endif
if(xyzzyaaat1(xyzzyaaad1))then
allocate(xyzzyaabv4(nspin,ndet,xyzzyaaax4),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','logdet_buff')
xyzzyaabv4=cmplx(0.d0,0.d0,dp)
endif
if(xyzzyaaat1(xyzzyaaae1))then
allocate(xyzzyaabo4(3,netot,real1_complex2,xyzzyaaax4),stat=xyzzyaaad4&
&)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','fidet_buff')
xyzzyaabo4=0.d0
endif
if(xyzzyaaat1(xyzzyaaaf1))then
allocate(xyzzyaabp4(3,ndet,nemax,real1_complex2,nspin,xyzzyaaax4),stat&
&=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','fi_prod_det_buff')
xyzzyaabp4=0.d0
endif
if(xyzzyaaat1(xyzzyaaag1))then
allocate(xyzzyaabj4(xyzzyaaax4),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','lapdet_buff')
xyzzyaabj4=0.d0
endif
if(xyzzyaaat1(xyzzyaaah1))then
allocate(xyzzyaabq4(ndet,nemax,real1_complex2,nspin,xyzzyaaax4),stat=x&
&yzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','prod_lapdet_buff')
xyzzyaabq4=0.d0
endif
if(xyzzyaaat1(xyzzyaaai1))then
allocate(xyzzyaabm4(xyzzyaaax4),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','local_potential_buff')
xyzzyaabm4=0.d0
endif
if(xyzzyaaat1(xyzzyaaaj1))then
allocate(xyzzyaabn4(xyzzyaaax4),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','nltot_buff')
xyzzyaabn4=0.d0
endif
if(xyzzyaaat1(xyzzyaaak1))then
allocate(xyzzyaabl4(xyzzyaaax4),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','stot_buff')
xyzzyaabl4=0.d0
endif
if(xyzzyaaat1(xyzzyaaal1))then
allocate(xyzzyaabr4(xyzzyaaax4),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','wdmc_buff')
xyzzyaabr4=0.d0
endif
if(xyzzyaaat1(xyzzyaaam1))then
allocate(xyzzyaabs4(xyzzyaaax4),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','valjas_buff')
xyzzyaabs4=0.d0
endif
if(xyzzyaaat1(xyzzyaaan1))then
allocate(xyzzyaabt4(xyzzyaaax4),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','logp_buff')
xyzzyaabt4=0.d0
endif
if(xyzzyaaat1(xyzzyaaao1))then
allocate(xyzzyaabu4(5,xyzzyaaax4),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'READ_CONFIGS','twist_buff')
xyzzyaabu4=0.d0
endif
endif
do xyzzyaaaa4=1,checkpoint_ncpu-1
xyzzyaaav4=xyzzyaabg1(xyzzyaaaa4)
xyzzyaaaw4=xyzzyaaaz1(xyzzyaaav4)
do xyzzyaaag4=1,xyzzyaaaw4
if(xyzzyaaat1(xyzzyaaaa1))read(xyzzyaaab4,err=1,end=2)xyzzyaabi4(:,:,x&
&yzzyaaag4)
if(xyzzyaaat1(xyzzyaaab1))read(xyzzyaaab4,err=1,end=2)xyzzyaabh4(:,xyz&
&zyaaag4)
if(xyzzyaaat1(xyzzyaaac1))read(xyzzyaaab4,err=1,end=2)xyzzyaabk4(xyzzy&
&aaag4)
if(xyzzyaaat1(xyzzyaaad1))read(xyzzyaaab4,err=1,end=2)xyzzyaabv4(:,:,x&
&yzzyaaag4)
if(xyzzyaaat1(xyzzyaaae1))read(xyzzyaaab4,err=1,end=2)xyzzyaabo4(:,:,:&
&,xyzzyaaag4)
if(xyzzyaaat1(xyzzyaaaf1))read(xyzzyaaab4,err=1,end=2)xyzzyaabp4(:,:,:&
&,:,:,xyzzyaaag4)
if(xyzzyaaat1(xyzzyaaag1))read(xyzzyaaab4,err=1,end=2)xyzzyaabj4(xyzzy&
&aaag4)
if(xyzzyaaat1(xyzzyaaah1))read(xyzzyaaab4,err=1,end=2)xyzzyaabq4(:,:,:&
&,:,xyzzyaaag4)
if(xyzzyaaat1(xyzzyaaai1))read(xyzzyaaab4,err=1,end=2)xyzzyaabm4(xyzzy&
&aaag4)
if(xyzzyaaat1(xyzzyaaaj1))read(xyzzyaaab4,err=1,end=2)xyzzyaabn4(xyzzy&
&aaag4)
if(xyzzyaaat1(xyzzyaaak1))read(xyzzyaaab4,err=1,end=2)xyzzyaabl4(xyzzy&
&aaag4)
if(xyzzyaaat1(xyzzyaaal1))read(xyzzyaaab4,err=1,end=2)xyzzyaabr4(xyzzy&
&aaag4)
if(xyzzyaaat1(xyzzyaaam1))read(xyzzyaaab4,err=1,end=2)xyzzyaabs4(xyzzy&
&aaag4)
if(xyzzyaaat1(xyzzyaaan1))read(xyzzyaaab4,err=1,end=2)xyzzyaabt4(xyzzy&
&aaag4)
if(xyzzyaaat1(xyzzyaaao1))read(xyzzyaaab4,err=1,end=2)xyzzyaabu4(:,xyz&
&zyaaag4)
enddo
if(xyzzyaaat1(xyzzyaaaa1))then
call qmpi_ssend(xyzzyaabi4(:,:,1:xyzzyaaaw4),xyzzyaaav4,xyzzyaaaa1,'RE&
&AD_CONFIGS','RELE_BUFF')
endif
if(xyzzyaaat1(xyzzyaaab1))then
call qmpi_ssend(xyzzyaabh4(:,1:xyzzyaaaw4),xyzzyaaav4,xyzzyaaab1,'READ&
&_CONFIGS','SELE_BUFF')
endif
if(xyzzyaaat1(xyzzyaaac1))then
call qmpi_ssend(xyzzyaabk4(1:xyzzyaaaw4),xyzzyaaav4,xyzzyaaac1,'READ_C&
&ONFIGS','ETOT_BUFF')
endif
if(xyzzyaaat1(xyzzyaaad1))then
call qmpi_ssend(xyzzyaabv4(:,:,1:xyzzyaaaw4),xyzzyaaav4,xyzzyaaad1,'RE&
&AD_CONFIGS','LOGDET_BUFF')
endif
if(xyzzyaaat1(xyzzyaaae1))then
call qmpi_ssend(xyzzyaabo4(:,:,:,1:xyzzyaaaw4),xyzzyaaav4,xyzzyaaae1,'&
&READ_CONFIGS','FIDET_BUFF')
endif
if(xyzzyaaat1(xyzzyaaaf1))then
call qmpi_ssend(xyzzyaabp4(:,:,:,:,:,1:xyzzyaaaw4),xyzzyaaav4,xyzzyaaa&
&f1,'READ_CONFIGS','FI_PROD_DET_BUFF')
endif
if(xyzzyaaat1(xyzzyaaag1))then
call qmpi_ssend(xyzzyaabj4(1:xyzzyaaaw4),xyzzyaaav4,xyzzyaaag1,'READ_C&
&ONFIGS','LAPDET_BUFF')
endif
if(xyzzyaaat1(xyzzyaaah1))then
call qmpi_ssend(xyzzyaabq4(:,:,:,:,1:xyzzyaaaw4),xyzzyaaav4,xyzzyaaah1&
&,'READ_CONFIGS','PROD_LAPDET_BUFF')
endif
if(xyzzyaaat1(xyzzyaaai1))then
call qmpi_ssend(xyzzyaabm4(1:xyzzyaaaw4),xyzzyaaav4,xyzzyaaai1,'READ_C&
&ONFIGS','LOCAL_POTENTIAL_BUFF')
endif
if(xyzzyaaat1(xyzzyaaaj1))then
call qmpi_ssend(xyzzyaabn4(1:xyzzyaaaw4),xyzzyaaav4,xyzzyaaaj1,'READ_C&
&ONFIGS','NLTOT_BUFF')
endif
if(xyzzyaaat1(xyzzyaaak1))then
call qmpi_ssend(xyzzyaabl4(1:xyzzyaaaw4),xyzzyaaav4,xyzzyaaak1,'READ_C&
&ONFIGS','STOT_BUFF')
endif
if(xyzzyaaat1(xyzzyaaal1))then
call qmpi_ssend(xyzzyaabr4(1:xyzzyaaaw4),xyzzyaaav4,xyzzyaaal1,'READ_C&
&ONFIGS','WDMC_BUFF')
endif
if(xyzzyaaat1(xyzzyaaam1))then
call qmpi_ssend(xyzzyaabs4(1:xyzzyaaaw4),xyzzyaaav4,xyzzyaaam1,'READ_C&
&ONFIGS','VALJAS_BUFF')
endif
if(xyzzyaaat1(xyzzyaaan1))then
call qmpi_ssend(xyzzyaabt4(1:xyzzyaaaw4),xyzzyaaav4,xyzzyaaan1,'READ_C&
&ONFIGS','LOGP_BUFF')
call checkmpi(ierror,'ssend logp_buff in read_configs')
endif
if(xyzzyaaat1(xyzzyaaao1))then
call qmpi_ssend(xyzzyaabu4(:,1:xyzzyaaaw4),xyzzyaaav4,xyzzyaaao1,'READ&
&_CONFIGS','TWIST_BUFF')
call checkmpi(ierror,'ssend twist_buff in read_configs')
endif
enddo
if(allocated(xyzzyaabi4))deallocate(xyzzyaabi4)
if(allocated(xyzzyaabh4))deallocate(xyzzyaabh4)
if(allocated(xyzzyaabk4))deallocate(xyzzyaabk4)
if(allocated(xyzzyaabv4))deallocate(xyzzyaabv4)
if(allocated(xyzzyaabo4))deallocate(xyzzyaabo4)
if(allocated(xyzzyaabp4))deallocate(xyzzyaabp4)
if(allocated(xyzzyaabj4))deallocate(xyzzyaabj4)
if(allocated(xyzzyaabq4))deallocate(xyzzyaabq4)
if(allocated(xyzzyaabm4))deallocate(xyzzyaabm4)
if(allocated(xyzzyaabn4))deallocate(xyzzyaabn4)
if(allocated(xyzzyaabl4))deallocate(xyzzyaabl4)
if(allocated(xyzzyaabr4))deallocate(xyzzyaabr4)
if(allocated(xyzzyaabs4))deallocate(xyzzyaabs4)
if(allocated(xyzzyaabt4))deallocate(xyzzyaabt4)
if(allocated(xyzzyaabu4))deallocate(xyzzyaabu4)
else
if(xyzzyaaat1(xyzzyaaaa1))then
call qmpi_recv(rele_config,ckmymaster,xyzzyaaaa1,'READ_CONFIGS','RELE_&
&CONFIG')
endif
if(xyzzyaaat1(xyzzyaaab1))then
call qmpi_recv(sele_config,ckmymaster,xyzzyaaab1,'READ_CONFIGS','SELE_&
&CONFIG')
endif
if(xyzzyaaat1(xyzzyaaac1))then
call qmpi_recv(etot_config,ckmymaster,xyzzyaaac1,'READ_CONFIGS','ETOT_&
&CONFIG')
endif
if(xyzzyaaat1(xyzzyaaad1))then
call qmpi_recv(logdet_config,ckmymaster,xyzzyaaad1,'READ_CONFIGS','LOG&
&DET_CONFIG')
endif
if(xyzzyaaat1(xyzzyaaae1))then
call qmpi_recv(fidet_config,ckmymaster,xyzzyaaae1,'READ_CONFIGS','FIDE&
&T_CONFIG')
endif
if(xyzzyaaat1(xyzzyaaaf1))then
call qmpi_recv(fi_prod_det_config,ckmymaster,xyzzyaaaf1,'READ_CONFIGS'&
&,'FI_PROD_DET_CONFIG')
endif
if(xyzzyaaat1(xyzzyaaag1))then
call qmpi_recv(lapdet_config,ckmymaster,xyzzyaaag1,'READ_CONFIGS','LAP&
&DET_CONFIG')
endif
if(xyzzyaaat1(xyzzyaaah1))then
call qmpi_recv(prod_lapdet_config,ckmymaster,xyzzyaaah1,'READ_CONFIGS'&
&,'PROD_LAPDET_CONFIG')
endif
if(xyzzyaaat1(xyzzyaaai1))then
call qmpi_recv(local_potential_config,ckmymaster,xyzzyaaai1,'READ_CONF&
&IGS','LOCAL_POTENTIAL_CONFIG')
endif
if(xyzzyaaat1(xyzzyaaaj1))then
call qmpi_recv(nltot_config,ckmymaster,xyzzyaaaj1,'READ_CONFIGS','NLTO&
&T_CONFIG')
endif
if(xyzzyaaat1(xyzzyaaak1))then
call qmpi_recv(stot_config,ckmymaster,xyzzyaaak1,'READ_CONFIGS','STOT_&
&CONFIG')
endif
if(xyzzyaaat1(xyzzyaaal1))then
call qmpi_recv(wdmc_config,ckmymaster,xyzzyaaal1,'READ_CONFIGS','WDMC_&
&CONFIG')
endif
if(xyzzyaaat1(xyzzyaaam1))then
call qmpi_recv(valjas_config,ckmymaster,xyzzyaaam1,'READ_CONFIGS','VAL&
&JAS_CONFIG')
endif
if(xyzzyaaat1(xyzzyaaan1))then
call qmpi_recv(logp_config,ckmymaster,xyzzyaaan1,'READ_CONFIGS','LOGP_&
&CONFIG')
endif
if(xyzzyaaat1(xyzzyaaao1))then
call qmpi_recv(twist_config,ckmymaster,xyzzyaaao1,'READ_CONFIGS','TWIS&
&T_CONFIG')
endif
endif
endif
if(am_master_group)then
close(xyzzyaaab4)
open_unit(xyzzyaaab4)=.false.
endif
xyzzyaaar1=.true.
if(allocated(xyzzyaaat4))deallocate(xyzzyaaat4)
return
1 call errstop('READ_CONFIGS','Problem reading config '//trim(i2s(xyzz&
&yaaag4))//' on node '//trim(i2s(my_node))//'.')
2 call errstop('READ_CONFIGS','End-of-file found when reading config '&
&//trim(i2s(xyzzyaaag4))//' on node '//trim(i2s(my_node))//'.')
end subroutine xyzzyaabh1
subroutine xyzzyaabi1(nconfig_inout,gen_by_expected,required_item,opti&
&onal_item,present_citem,required_extra,optional_extra,present_eitem)
implicit none
integer,intent(inout) :: nconfig_inout
logical,dimension(:),intent(out) :: present_citem,present_eitem
character(20),intent(in) :: gen_by_expected
character(20),dimension(:),intent(in) :: required_item,optional_item,r&
&equired_extra,optional_extra
integer xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5,xyzzyaaad5,xyzzyaaae5,xyzzyaa&
&af5,xyzzyaaag5,xyzzyaaah5,xyzzyaaai5,xyzzyaaaj5
logical,allocatable :: xyzzyaaak5(:),xyzzyaaal5(:),xyzzyaaam5(:),xyzzy&
&aaan5(:)
character(20) label
xyzzyaaaa5=size(required_item)
xyzzyaaab5=size(optional_item)
xyzzyaaac5=size(required_extra)
xyzzyaaad5=size(optional_extra)
xyzzyaaaf5=0
if(.not.xyzzyaaar1.and.(xyzzyaaaa5/=0.or.xyzzyaaac5/=0))call errstop_m&
&aster('CHECK_CONFIGS','Previous QMC data requested, but no data in me&
&mory and config file not found.')
if(xyzzyaaau1<nconfig_inout)call errstop('CHECK_CONFIGS','Have loaded &
&fewer configs than required on node '//trim(i2s(my_node))//'.')
if(nconfig_inout==-1.and.xyzzyaaav1/=xyzzyaaau1*nnodes)call errstop('C&
&HECK_CONFIGS','Requested an equal number of configs per node but some&
&how have loaded the wrong number on node '//trim(i2s(my_node))//'.')
nconfig_inout=xyzzyaaau1
xyzzyaaaj5=0
xyzzyaaaj5=index(gen_by_expected,trim(xyzzyaabc1))
if(xyzzyaaaj5==0.and.trim(gen_by_expected)/='ANY'.and.trim(xyzzyaabc1)&
&/='NONE')call errstop_master('CHECK_CONFIGS','Loaded configs were gen&
&erated by '//trim(xyzzyaabc1)//' instead of '//trim(gen_by_expected)/&
&/'. Quitting.')
if(xyzzyaaaa5+xyzzyaaab5>0.and.xyzzyaaav1>0.and.xyzzyaaar1)then
allocate(xyzzyaaak5(xyzzyaaaw1),xyzzyaaal5(xyzzyaaaw1),xyzzyaaam5(xyzz&
&yaaaa5),xyzzyaaan5(xyzzyaaab5),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'CHECK_CONFIGS','L_ITEM')
xyzzyaaak5=.false.
xyzzyaaal5=.false.
xyzzyaaam5=.false.
xyzzyaaan5=.false.
do xyzzyaaae5=1,xyzzyaaaw1
label=xyzzyaaba1(xyzzyaaae5)
do xyzzyaaag5=1,xyzzyaaaa5
if(trim(label)==trim(required_item(xyzzyaaag5)))then
if(xyzzyaaam5(xyzzyaaag5))call errstop_master('CHECK_CONFIGS','Config &
&item '//trim(label)//' requested more than once.')
xyzzyaaak5(xyzzyaaae5)=.true.
xyzzyaaam5(xyzzyaaag5)=.true.
endif
enddo
do xyzzyaaah5=1,xyzzyaaab5
if(trim(label)==trim(optional_item(xyzzyaaah5)))then
if(xyzzyaaan5(xyzzyaaah5))call errstop_master('CHECK_CONFIGS','Optiona&
&l config item '//trim(label)//' requested more than once.')
if(xyzzyaaak5(xyzzyaaae5))call errstop_master('CHECK_CONFIGS','Config &
&item '//trim(label)//' requested more than once (at least once as req&
&uired and once as optional).')
xyzzyaaal5(xyzzyaaae5)=.true.
xyzzyaaan5(xyzzyaaah5)=.true.
endif
enddo
enddo
if(any(.not.xyzzyaaam5))then
if(am_master)then
call wout('The following config data is needed but is not present in t&
&he config file:')
do xyzzyaaag5=1,xyzzyaaaa5
if(.not.xyzzyaaam5(xyzzyaaag5))call wout(trim(required_item(xyzzyaaag5&
&)))
enddo
call wout()
endif
call errstop_master('CHECK_CONFIGS','Stopping.')
endif
do xyzzyaaae5=1,xyzzyaaaw1
label=xyzzyaaba1(xyzzyaaae5)
if(xyzzyaaak5(xyzzyaaae5).or.xyzzyaaal5(xyzzyaaae5))cycle
select case(trim(label))
case('RELE')
deallocate(rele_config)
case('SELE')
deallocate(sele_config)
case('ETOT')
deallocate(etot_config)
case('LOGDET')
deallocate(logdet_config)
case('FIDET')
deallocate(fidet_config)
case('FI_PROD_DET')
deallocate(fi_prod_det_config)
case('LAPDET')
deallocate(lapdet_config)
case('PROD_LAPDET')
deallocate(prod_lapdet_config)
case('LOCAL_POTENTIAL')
deallocate(local_potential_config)
case('NLTOT')
deallocate(nltot_config)
case('STOT')
deallocate(stot_config)
case('WDMC')
deallocate(wdmc_config)
case('VALJAS')
deallocate(valjas_config)
case('LOGP')
deallocate(logp_config)
case('TWIST')
deallocate(twist_config)
end select
enddo
do xyzzyaaah5=1,xyzzyaaab5
label=optional_item(xyzzyaaah5)
select case(trim(label))
case('RELE')
if(.not.allocated(rele_config))allocate(rele_config(3,netot,xyzzyaaau1&
&),stat=xyzzyaaaf5)
case('SELE')
if(.not.allocated(sele_config))allocate(sele_config(netot,xyzzyaaau1),&
&stat=xyzzyaaaf5)
case('ETOT')
if(.not.allocated(etot_config))allocate(etot_config(xyzzyaaau1),stat=x&
&yzzyaaaf5)
case('LOGDET')
if(.not.allocated(logdet_config))allocate(logdet_config(nspin,ndet,xyz&
&zyaaau1),stat=xyzzyaaaf5)
case('FIDET')
if(.not.allocated(fidet_config))allocate(fidet_config(3,netot,real1_co&
&mplex2,xyzzyaaau1),stat=xyzzyaaaf5)
case('FI_PROD_DET')
if(.not.allocated(fi_prod_det_config))allocate(fi_prod_det_config(3,nd&
&et,nemax,real1_complex2,nspin,xyzzyaaau1),stat=xyzzyaaaf5)
case('LAPDET')
if(.not.allocated(lapdet_config))allocate(lapdet_config(xyzzyaaau1),st&
&at=xyzzyaaaf5)
case('PROD_LAPDET')
if(.not.allocated(prod_lapdet_config))allocate(prod_lapdet_config(ndet&
&,nemax,real1_complex2,nspin,xyzzyaaau1),stat=xyzzyaaaf5)
case('LOCAL_POTENTIAL')
if(.not.allocated(local_potential_config))allocate(local_potential_con&
&fig(xyzzyaaau1),stat=xyzzyaaaf5)
case('NLTOT')
if(.not.allocated(nltot_config))allocate(nltot_config(xyzzyaaau1),stat&
&=xyzzyaaaf5)
case('STOT')
if(.not.allocated(stot_config))allocate(stot_config(xyzzyaaau1),stat=x&
&yzzyaaaf5)
case('WDMC')
if(.not.allocated(wdmc_config))allocate(wdmc_config(xyzzyaaau1),stat=x&
&yzzyaaaf5)
case('VALJAS')
if(.not.allocated(valjas_config))then
allocate(valjas_config(xyzzyaaau1),stat=xyzzyaaaf5)
if(xyzzyaaaf5==0)valjas_config=0.d0
endif
case('LOGP')
if(.not.allocated(logp_config))allocate(logp_config(xyzzyaaau1),stat=x&
&yzzyaaaf5)
case('TWIST')
if(.not.allocated(twist_config))allocate(twist_config(5,xyzzyaaau1),st&
&at=xyzzyaaaf5)
end select
call check_alloc(xyzzyaaaf5,'CHECK_CONFIGS',trim(xyzzyaaba1(xyzzyaaah5&
&)))
enddo
if(xyzzyaaab5>0)present_citem=xyzzyaaan5
deallocate(xyzzyaaak5,xyzzyaaal5,xyzzyaaam5,xyzzyaaan5)
elseif(xyzzyaaav1>0)then
if(allocated(rele_config))deallocate(rele_config)
if(allocated(sele_config))deallocate(sele_config)
if(allocated(logdet_config))deallocate(logdet_config)
if(allocated(lapdet_config))deallocate(lapdet_config)
if(allocated(etot_config))deallocate(etot_config)
if(allocated(stot_config))deallocate(stot_config)
if(allocated(local_potential_config))deallocate(local_potential_config&
&)
if(allocated(nltot_config))deallocate(nltot_config)
if(allocated(fidet_config))deallocate(fidet_config)
if(allocated(fi_prod_det_config))deallocate(fi_prod_det_config)
if(allocated(prod_lapdet_config))deallocate(prod_lapdet_config)
if(allocated(wdmc_config))deallocate(wdmc_config)
if(allocated(valjas_config))deallocate(valjas_config)
if(allocated(logp_config))deallocate(logp_config)
if(allocated(twist_config))deallocate(twist_config)
if(allocated(xyzzyaaba1))deallocate(xyzzyaaba1)
if(allocated(xyzzyaaaz1))deallocate(xyzzyaaaz1)
xyzzyaaau1=0
xyzzyaaav1=0
if(xyzzyaaab5>0)present_eitem=.false.
endif
if(xyzzyaaac5+xyzzyaaad5>0.and.xyzzyaaar1)then
allocate(xyzzyaaak5(xyzzyaaay1),xyzzyaaal5(xyzzyaaay1),xyzzyaaam5(xyzz&
&yaaac5),xyzzyaaan5(xyzzyaaad5),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'CHECK_CONFIGS','L_ITEM <2>')
xyzzyaaak5=.false.
xyzzyaaal5=.false.
xyzzyaaam5=.false.
xyzzyaaan5=.false.
do xyzzyaaai5=1,xyzzyaaay1
label=xyzzyaabb1(xyzzyaaai5)
do xyzzyaaag5=1,xyzzyaaac5
if(trim(label)==trim(required_extra(xyzzyaaag5)))then
if(xyzzyaaam5(xyzzyaaag5))call errstop_master('CHECK_CONFIGS','Extra i&
&tem '//trim(label)//' requested more than once.')
xyzzyaaak5(xyzzyaaai5)=.true.
xyzzyaaam5(xyzzyaaag5)=.true.
endif
enddo
do xyzzyaaah5=1,xyzzyaaad5
if(trim(label)==trim(optional_extra(xyzzyaaah5)))then
if(xyzzyaaan5(xyzzyaaah5))call errstop_master('CHECK_CONFIGS','Optiona&
&l extra item '//trim(label)//' requested more than once.')
if(xyzzyaaak5(xyzzyaaai5))call errstop_master('CHECK_CONFIGS','Extra i&
&tem '//trim(label)//' requested more than once (at least once as requ&
&ired and once as optional).')
xyzzyaaal5(xyzzyaaai5)=.true.
xyzzyaaan5(xyzzyaaah5)=.true.
endif
enddo
enddo
if(any(.not.xyzzyaaam5))then
if(am_master)then
call wout('The following extra data is needed but is not present in th&
&e config file:')
do xyzzyaaag5=1,xyzzyaaac5
if(.not.xyzzyaaam5(xyzzyaaag5))call wout(trim(required_extra(xyzzyaaag&
&5)))
enddo
call wout()
endif
call errstop_master('CHECK_CONFIGS','Stopping.')
endif
do xyzzyaaai5=1,xyzzyaaay1
if(xyzzyaaak5(xyzzyaaai5).or.xyzzyaaal5(xyzzyaaai5))cycle
label=xyzzyaabb1(xyzzyaaai5)
select case(trim(label))
case('VMC_SAVED_STATE')
if(allocated(dtvmc_array_config))deallocate(dtvmc_array_config)
if(allocated(rele_vmc_config))deallocate(rele_vmc_config)
if(allocated(sele_vmc_config))deallocate(sele_vmc_config)
case('DMC_SAVED_STATE')
if(allocated(log_pi_wt_array_config))deallocate(log_pi_wt_array_config&
&)
if(allocated(log_pi_wt_array2_config))deallocate(log_pi_wt_array2_conf&
&ig)
end select
enddo
do xyzzyaaah5=1,xyzzyaaad5
label=optional_extra(xyzzyaaah5)
select case(trim(label))
case('VMC_SAVED_STATE')
if(.not.allocated(dtvmc_array_config))allocate(dtvmc_array_config(no_d&
&ifftypes),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'CHECK_CONFIGS','DTVMC_ARRAY_CONFIG')
if(.not.allocated(rele_vmc_config))allocate(rele_vmc_config(3,netot),s&
&tat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'CHECK_CONFIGS','RELE_VMC_CONFIG')
if(.not.allocated(sele_vmc_config))allocate(sele_vmc_config(netot),sta&
&t=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'CHECK_CONFIGS','SELE_VMC_CONFIG')
case('DMC_SAVED_STATE')
if(tpdmc_config>0)then
if(.not.allocated(log_pi_wt_array_config))allocate(log_pi_wt_array_con&
&fig(0:tpdmc_config-1),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'CHECK_CONFIGS','LOG_PI_WT_ARRAY_CONFIG')
endif
if(.not.allocated(log_pi_wt_array2_config))allocate(log_pi_wt_array2_c&
&onfig(0:tpdmc_config),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'CHECK_CONFIGS','LOG_PI_WT_ARRAY2_CONFIG')
end select
enddo
if(xyzzyaaad5>0)present_eitem=xyzzyaaan5
deallocate(xyzzyaaak5,xyzzyaaal5,xyzzyaaam5,xyzzyaaan5)
else
do xyzzyaaai5=1,xyzzyaaay1
label=xyzzyaabb1(xyzzyaaai5)
select case(trim(label))
case('VMC_SAVED_STATE')
if(allocated(dtvmc_array_config))deallocate(dtvmc_array_config)
if(allocated(rele_vmc_config))deallocate(rele_vmc_config)
if(allocated(sele_vmc_config))deallocate(sele_vmc_config)
case('DMC_SAVED_STATE')
if(allocated(log_pi_wt_array_config))deallocate(log_pi_wt_array_config&
&)
if(allocated(log_pi_wt_array2_config))deallocate(log_pi_wt_array2_conf&
&ig)
end select
enddo
if(xyzzyaaad5>0)present_eitem=.false.
endif
end subroutine xyzzyaabi1
subroutine write_configs(rngbuf)
use slaarnabg, only : atom_basis_type,nitot,rion
use shalloc, only : shallocate,deshallocate,shallocate_barrier,need_sh&
&m,nnpsmp,nsmps,smp_masters,my_smpproc,am_smpmaster
implicit none
logical,intent(in) :: rngbuf
integer xyzzyaaaa6,xyzzyaaab6,xyzzyaaac6,xyzzyaaad6,xyzzyaaae6,xyzzyaa&
&af6,xyzzyaaag6,xyzzyaaah6,xyzzyaaai6,xyzzyaaaj6,xyzzyaaak6,xyzzyaaal6&
&,xyzzyaaam6(0:nnodes-1)
integer,allocatable :: xyzzyaaan6(:,:)
integer,pointer :: xyzzyaaao6(:,:)=>null(),xyzzyaaap6(:,:)=>null(),xyz&
&zyaaaq6(:)=>null(),xyzzyaaar6(:)=>null()
real(sp) xyzzyaaas6
real(sp),pointer :: xyzzyaaat6(:)=>null()
real(dp),allocatable :: xyzzyaaau6(:,:,:),xyzzyaaav6(:),xyzzyaaaw6(:),&
&xyzzyaaax6(:),xyzzyaaay6(:),xyzzyaaaz6(:),xyzzyaaba6(:,:,:,:),xyzzyaa&
&bb6(:,:,:,:,:,:),xyzzyaabc6(:,:,:,:,:),xyzzyaabd6(:),xyzzyaabe6(:),xy&
&zzyaabf6(:),xyzzyaabg6(:,:)
real(dp),pointer :: xyzzyaabh6(:,:,:)=>null()
complex(dp),allocatable :: xyzzyaabi6(:,:,:)
character(20) label
call timer('WRITE_CONFIGS',.true.)
xyzzyaaam6(:)=0
if(xyzzyaaav1>0)then
xyzzyaaam6(:)=xyzzyaaaz1(:)
if(xyzzyaaas1)xyzzyaaam6(:)=xyzzyaaax1
endif
xyzzyaaaj6=xyzzyaaam6(my_node)
xyzzyaaak6=sum(xyzzyaaam6)
if(am_master)then
call open_units(xyzzyaaac6,xyzzyaaad6)
if(xyzzyaaad6/=0)call errstop('WRITE_CONFIGS','Could not find free i/o&
& unit number.')
open(unit=xyzzyaaac6,file=trim(con_out),form='unformatted',status='rep&
&lace',iostat=xyzzyaaad6)
if(xyzzyaaad6/=0)call errstop('WRITE_CONFIGS','Problem opening config &
&file.')
call get_total_time(xyzzyaaas6)
call xyzzyaabj1(xyzzyaaac6,'INFO')
call xyzzyaabj1(xyzzyaaac6,'NEXTRA')
write(xyzzyaaac6)xyzzyaaay1
call xyzzyaabj1(xyzzyaaac6,'GEN_BY')
write(xyzzyaaac6)xyzzyaabc1
call xyzzyaabj1(xyzzyaaac6,'NSPIN')
write(xyzzyaaac6)nspin
call xyzzyaabj1(xyzzyaaac6,'NELE')
write(xyzzyaaac6)nele
call xyzzyaabj1(xyzzyaaac6,'NDET')
write(xyzzyaaac6)ndet
call xyzzyaabj1(xyzzyaaac6,'COMPLEX_WF')
write(xyzzyaaac6)log2int(complex_wf)
call xyzzyaabj1(xyzzyaaac6,'NONCOLL_SPIN')
write(xyzzyaaac6)log2int(noncoll_spin)
call xyzzyaabj1(xyzzyaaac6,'INTERACTION')
write(xyzzyaaac6)interaction
call xyzzyaabj1(xyzzyaaac6,'ATOM_BASIS_TYPE')
write(xyzzyaaac6)atom_basis_type
call xyzzyaabj1(xyzzyaaac6,'TOTAL_TIME')
write(xyzzyaaac6)xyzzyaaas6
call xyzzyaabj1(xyzzyaaac6,'END INFO')
if(xyzzyaaak6>0)then
call xyzzyaabj1(xyzzyaaac6,'DEFINE_CONFIGS')
write(xyzzyaaac6)xyzzyaaak6
write(xyzzyaaac6)xyzzyaaaw1
do xyzzyaaae6=1,xyzzyaaaw1
call xyzzyaabj1(xyzzyaaac6,xyzzyaaba1(xyzzyaaae6))
enddo
call xyzzyaabj1(xyzzyaaac6,'END DEFINE_CONFIGS')
endif
do xyzzyaaaf6=1,xyzzyaaay1
label=xyzzyaabb1(xyzzyaaaf6)
call xyzzyaabj1(xyzzyaaac6,(trim(label)))
select case(trim(label))
case('VMC_SAVED_STATE')
call xyzzyaabj1(xyzzyaaac6,'NO_DIFFTYPES')
write(xyzzyaaac6)no_difftypes
call xyzzyaabj1(xyzzyaaac6,'DTVMC')
write(xyzzyaaac6)dtvmc_array_config
call xyzzyaabj1(xyzzyaaac6,'VMC_STEPS')
write(xyzzyaaac6)vmc_steps_config
call xyzzyaabj1(xyzzyaaac6,'NNODES')
write(xyzzyaaac6)nnodes
call xyzzyaabj1(xyzzyaaac6,'RELE')
if(need_shm)then
call shallocate(xyzzyaabh6,3,netot,nnpsmp,stat=xyzzyaaai6)
call check_alloc(xyzzyaaai6,'write_configs','rele_vmc_config_all')
write(xyzzyaaac6)rele_vmc_config(:,:)
call shallocate_barrier
do xyzzyaaag6=2,nnpsmp
write(xyzzyaaac6)xyzzyaabh6(:,:,xyzzyaaag6)
enddo
do xyzzyaaaa6=2,nsmps
call qmpi_recv(xyzzyaabh6,smp_masters(xyzzyaaaa6),load_msg,'WRITE_CONF&
&IGS','rele_vmc_config_all')
call checkmpi(ierror,'receiving rele_vmc_config_all in write_configs')
do xyzzyaaag6=1,nnpsmp
write(xyzzyaaac6)xyzzyaabh6(:,:,xyzzyaaag6)
enddo
enddo
call shallocate_barrier
call deshallocate(xyzzyaabh6)
else
allocate(xyzzyaabh6(3,netot,nnodes),stat=xyzzyaaai6)
call check_alloc(xyzzyaaai6,'write_configs','rele_vmc_config_all')
call mpi_gather(rele_vmc_config,three_netot,mpi_double_precision,xyzzy&
&aabh6,three_netot,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gather rele_vmc in write_configs')
do xyzzyaaag6=1,nnodes
write(xyzzyaaac6)xyzzyaabh6(:,:,xyzzyaaag6)
enddo
deallocate(xyzzyaabh6)
endif
if(noncoll_spin)then
call xyzzyaabj1(xyzzyaaac6,'SELE')
if(need_shm)then
call shallocate(xyzzyaaao6,netot,nnpsmp,stat=xyzzyaaai6)
call check_alloc(xyzzyaaai6,'write_configs','sele_vmc_config_all')
write(xyzzyaaac6)sele_vmc_config(:)
call shallocate_barrier
do xyzzyaaag6=2,nnpsmp
write(xyzzyaaac6)xyzzyaaao6(:,xyzzyaaag6)
enddo
do xyzzyaaaa6=2,nsmps
call qmpi_recv(xyzzyaaao6,smp_masters(xyzzyaaaa6),load_msg,'WRITE_CONF&
&IGS','sele_vmc_config_all')
call checkmpi(ierror,'receiving sele_vmc_config_all in write_configs')
do xyzzyaaag6=1,nnpsmp
write(xyzzyaaac6)xyzzyaaao6(:,xyzzyaaag6)
enddo
enddo
call shallocate_barrier
call deshallocate(xyzzyaaao6)
else
allocate(xyzzyaaao6(netot,nnodes),stat=xyzzyaaai6)
call check_alloc(xyzzyaaai6,'write_configs','sele_vmc_config_all')
call mpi_gather(sele_vmc_config,netot,mpi_integer,xyzzyaaao6,netot,mpi&
&_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gather sele_vmc in write_configs')
do xyzzyaaag6=1,nnodes
write(xyzzyaaac6)xyzzyaaao6(:,xyzzyaaag6)
enddo
deallocate(xyzzyaaao6)
endif
endif
case('FINAL_VMC_RESULT')
call xyzzyaabj1(xyzzyaaac6,'ENERGY')
write(xyzzyaaac6)final_vmce_config
call xyzzyaabj1(xyzzyaaac6,'ERRORBAR')
write(xyzzyaaac6)final_vmcde_config
call xyzzyaabj1(xyzzyaaac6,'ERRORBARU')
write(xyzzyaaac6)final_vmcdeu_config
call xyzzyaabj1(xyzzyaaac6,'VARIANCE')
write(xyzzyaaac6)final_vmcvar_config
call xyzzyaabj1(xyzzyaaac6,'TOTMOVE')
write(xyzzyaaac6)final_vmcmove_config
call xyzzyaabj1(xyzzyaaac6,'CORRTIME')
write(xyzzyaaac6)final_vmctau_config
if(trim(interaction)=='ewald_mpc'.or.trim(interaction)=='mpc_ewald')th&
&en
call xyzzyaabj1(xyzzyaaac6,'ENERGY2')
write(xyzzyaaac6)final_vmce2_config
call xyzzyaabj1(xyzzyaaac6,'ERRORBAR2')
write(xyzzyaaac6)final_vmcde2_config
call xyzzyaabj1(xyzzyaaac6,'ERRORBAR2U')
write(xyzzyaaac6)final_vmcde2u_config
endif
case('DMC_SAVED_STATE')
call xyzzyaabj1(xyzzyaaac6,'EBEST')
write(xyzzyaaac6)ebest_config
call xyzzyaabj1(xyzzyaaac6,'EBEST_INIT')
write(xyzzyaaac6)ebest_init_config
call xyzzyaabj1(xyzzyaaac6,'EREF')
write(xyzzyaaac6)eref_config
call xyzzyaabj1(xyzzyaaac6,'DTEFF_EBEST_INIT')
write(xyzzyaaac6)dteff_ebest_init_config
call xyzzyaabj1(xyzzyaaac6,'DTEFF_BEST')
write(xyzzyaaac6)dteff_best_config
call xyzzyaabj1(xyzzyaaac6,'DMCEQUIL_STEPS')
write(xyzzyaaac6)dmcequil_steps_config
call xyzzyaabj1(xyzzyaaac6,'DMCSTATS_STEPS')
write(xyzzyaaac6)dmcstats_steps_config
call xyzzyaabj1(xyzzyaaac6,'TPDMC')
write(xyzzyaaac6)tpdmc_config
call xyzzyaabj1(xyzzyaaac6,'LWDMC')
write(xyzzyaaac6)log2int(lwdmc_config)
call xyzzyaabj1(xyzzyaaac6,'GROWTH_ESTIMATOR')
write(xyzzyaaac6)log2int(growth_estimator_config)
call xyzzyaabj1(xyzzyaaac6,'NUMER_EXPECT')
write(xyzzyaaac6)numer_expect_config
call xyzzyaabj1(xyzzyaaac6,'DENOMINATOR_WT')
write(xyzzyaaac6)denominator_wt_config
call xyzzyaabj1(xyzzyaaac6,'LOG_PI_WT')
write(xyzzyaaac6)log_pi_wt_config
if(tpdmc_config>0)then
call xyzzyaabj1(xyzzyaaac6,'LOG_PI_WT_ARRAY')
write(xyzzyaaac6)log_pi_wt_array_config
endif
if(growth_estimator_config)then
call xyzzyaabj1(xyzzyaaac6,'NUMERATOR_WT2')
write(xyzzyaaac6)numerator_wt2_config
call xyzzyaabj1(xyzzyaaac6,'DENOMINATOR_WT2')
write(xyzzyaaac6)denominator_wt2_config
call xyzzyaabj1(xyzzyaaac6,'LOG_PI_WT2')
write(xyzzyaaac6)log_pi_wt2_config
call xyzzyaabj1(xyzzyaaac6,'LOG_PI_WT_ARRAY2')
write(xyzzyaaac6)log_pi_wt_array2_config
endif
case('REBLOCK_DATA')
write(xyzzyaaac6)reblock_nbs_config,reblock_nstep_config,reblock_nobs_&
&config,reblock_sum_w_config,log2int(popstats_config)
write(xyzzyaaac6)reblock_block_length_config,reblock_s_u_f_config,rebl&
&ock_sum_w2_closed_config,reblock_sum_w_open_config,reblock_sum_o2_clo&
&sed_config,reblock_sum_o_open_config,reblock_sum_o_config
if(popstats_config)then
write(xyzzyaaac6)reblock_sum_w2_config,reblock_sum_ow0_config,reblock_&
&sum_ow2_config,reblock_sum_o2w0_config,reblock_sum_o2w2_config
endif
case('RANDOM')
call xyzzyaabj1(xyzzyaaac6,'NNODES')
write(xyzzyaaac6)nnodes
if(need_shm)then
call shallocate(xyzzyaaap6,25,nnpsmp,stat=xyzzyaaai6)
call check_alloc(xyzzyaaai6,'write_configs','random_state_config_all')
write(xyzzyaaac6)random_state_config(:)
call shallocate_barrier
do xyzzyaaag6=2,nnpsmp
write(xyzzyaaac6)xyzzyaaap6(:,xyzzyaaag6)
enddo
do xyzzyaaaa6=2,nsmps
call qmpi_recv(xyzzyaaap6,smp_masters(xyzzyaaaa6),load_msg,'WRITE_CONF&
&IGS','random_state_config_all')
call checkmpi(ierror,'receiving random_state_config_all in write_confi&
&gs')
do xyzzyaaag6=1,nnpsmp
write(xyzzyaaac6)xyzzyaaap6(:,xyzzyaaag6)
enddo
enddo
call shallocate_barrier
call deshallocate(xyzzyaaap6)
else
allocate(xyzzyaaap6(25,nnodes),stat=xyzzyaaai6)
call check_alloc(xyzzyaaai6,'write_configs','random_state_config_all')
call mpi_gather(random_state_config,25,mpi_integer,xyzzyaaap6,25,mpi_i&
&nteger,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gather random_state in write_configs')
do xyzzyaaag6=1,nnodes
write(xyzzyaaac6)xyzzyaaap6(:,xyzzyaaag6)
enddo
deallocate(xyzzyaaap6)
endif
if(rngbuf)then
call xyzzyaabj1(xyzzyaaac6,'BUFFER')
write(xyzzyaaac6)ranx_max,ranx_gauss_max
if(need_shm)then
call shallocate(xyzzyaaaq6,nnpsmp,stat=xyzzyaaai6)
call check_alloc(xyzzyaaai6,'write_configs','ranx_indx_all')
call shallocate(xyzzyaaar6,nnpsmp,stat=xyzzyaaai6)
call check_alloc(xyzzyaaai6,'write_configs','ranx_gauss_indx_all')
write(xyzzyaaac6)ranx_indx,ranx_gauss_indx
call shallocate_barrier
do xyzzyaaag6=2,nnpsmp
write(xyzzyaaac6)xyzzyaaaq6(xyzzyaaag6),xyzzyaaar6(xyzzyaaag6)
enddo
do xyzzyaaaa6=2,nsmps
call qmpi_recv(xyzzyaaaq6,smp_masters(xyzzyaaaa6),load_msg,'WRITE_CONF&
&IGS','ranx_indx_all')
call checkmpi(ierror,'receiving ranx_indx_all in write_configs')
call qmpi_recv(xyzzyaaar6,smp_masters(xyzzyaaaa6),load_msg,'WRITE_CONF&
&IGS','ranx_gauss_indx_all')
call checkmpi(ierror,'receiving ranx_gauss_indx_all in write_configs')
do xyzzyaaag6=1,nnpsmp
write(xyzzyaaac6)xyzzyaaaq6(xyzzyaaag6),xyzzyaaar6(xyzzyaaag6)
enddo
enddo
call shallocate_barrier
call deshallocate(xyzzyaaaq6)
call deshallocate(xyzzyaaar6)
else
allocate(xyzzyaaaq6(nnodes),stat=xyzzyaaai6)
call check_alloc(xyzzyaaai6,'write_configs','ranx_indx_all')
allocate(xyzzyaaar6(nnodes),stat=xyzzyaaai6)
call check_alloc(xyzzyaaai6,'write_configs','ranx_gauss_indx_all')
call mpi_gather(ranx_indx,1,mpi_integer,xyzzyaaaq6,1,mpi_integer,0,mpi&
&_comm_world,ierror)
call checkmpi(ierror,'gather ranx_indx in write_configs')
call mpi_gather(ranx_gauss_indx,1,mpi_integer,xyzzyaaar6,1,mpi_integer&
&,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gather ranx_gauss_indx in write_configs')
do xyzzyaaag6=1,nnodes
write(xyzzyaaac6)xyzzyaaaq6(xyzzyaaag6),xyzzyaaar6(xyzzyaaag6)
enddo
deallocate(xyzzyaaaq6,xyzzyaaar6)
endif
if(need_shm)then
call shallocate(xyzzyaaat6,ranx_max*nnpsmp,stat=xyzzyaaai6)
call check_alloc(xyzzyaaai6,'write_configs','rng_buffer_all')
write(xyzzyaaac6)ranx_buffer(:)
call shallocate_barrier
xyzzyaaab6=ranx_max+1
do xyzzyaaag6=2,nnpsmp
write(xyzzyaaac6)xyzzyaaat6(xyzzyaaab6:xyzzyaaab6+ranx_max-1)
xyzzyaaab6=xyzzyaaab6+ranx_max
enddo
do xyzzyaaaa6=2,nsmps
call qmpi_recv(xyzzyaaat6,smp_masters(xyzzyaaaa6),load_msg,'WRITE_CONF&
&IGS','rng_buffer_all')
call checkmpi(ierror,'receiving rng_buffer_all in write_configs')
xyzzyaaab6=1
do xyzzyaaag6=1,nnpsmp
write(xyzzyaaac6)xyzzyaaat6(xyzzyaaab6:xyzzyaaab6+ranx_max-1)
xyzzyaaab6=xyzzyaaab6+ranx_max
enddo
enddo
call shallocate_barrier
call deshallocate(xyzzyaaat6)
else
allocate(xyzzyaaat6(ranx_max*nnodes),stat=xyzzyaaai6)
call check_alloc(xyzzyaaai6,'write_configs','rng_buffer_all')
call mpi_gather(ranx_buffer,ranx_max,mpi_real,xyzzyaaat6,ranx_max,mpi_&
&real,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gather ranx_buffer in write_configs')
xyzzyaaaa6=1
do xyzzyaaag6=1,nnodes
write(xyzzyaaac6)xyzzyaaat6(xyzzyaaaa6:xyzzyaaaa6+ranx_max-1)
xyzzyaaaa6=xyzzyaaaa6+ranx_max
enddo
deallocate(xyzzyaaat6)
endif
if(need_shm)then
call shallocate(xyzzyaaat6,ranx_gauss_max*nnpsmp,stat=xyzzyaaai6)
call check_alloc(xyzzyaaai6,'write_configs','rng_buffer_all <2>')
write(xyzzyaaac6)ranx_gauss_buffer(:)
call shallocate_barrier
xyzzyaaab6=ranx_gauss_max+1
do xyzzyaaag6=2,nnpsmp
write(xyzzyaaac6)xyzzyaaat6(xyzzyaaab6:xyzzyaaab6+ranx_gauss_max-1)
xyzzyaaab6=xyzzyaaab6+ranx_gauss_max
enddo
do xyzzyaaaa6=2,nsmps
call qmpi_recv(xyzzyaaat6,smp_masters(xyzzyaaaa6),load_msg,'WRITE_CONF&
&IGS','rng_buffer_all')
call checkmpi(ierror,'receiving rng_buffer_all in write_configs <2>')
xyzzyaaab6=1
do xyzzyaaag6=1,nnpsmp
write(xyzzyaaac6)xyzzyaaat6(xyzzyaaab6:xyzzyaaab6+ranx_gauss_max-1)
xyzzyaaab6=xyzzyaaab6+ranx_gauss_max
enddo
enddo
call shallocate_barrier
call deshallocate(xyzzyaaat6)
else
allocate(xyzzyaaat6(ranx_gauss_max*nnodes),stat=xyzzyaaai6)
call check_alloc(xyzzyaaai6,'write_configs','rng_buffer_all <2>')
call mpi_gather(ranx_gauss_buffer,ranx_gauss_max,mpi_real,xyzzyaaat6,r&
&anx_gauss_max,mpi_real,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gather ranx_gauss_buffer in write_configs')
xyzzyaaaa6=1
do xyzzyaaag6=1,nnodes
write(xyzzyaaac6)xyzzyaaat6(xyzzyaaaa6:xyzzyaaaa6+ranx_gauss_max-1)
xyzzyaaaa6=xyzzyaaaa6+ranx_gauss_max
enddo
deallocate(xyzzyaaat6)
endif
endif
case('GEOMETRY')
write(xyzzyaaac6)nitot
if(nitot>0)write(xyzzyaaac6)rion
end select
call xyzzyaabj1(xyzzyaaac6,'END '//trim(label))
enddo
if(xyzzyaaak6>0)then
call xyzzyaabj1(xyzzyaaac6,'CONFIGS')
do xyzzyaaah6=1,xyzzyaaaj6
if(xyzzyaaat1(xyzzyaaaa1))write(xyzzyaaac6)rele_config(:,:,xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaab1))write(xyzzyaaac6)sele_config(:,xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaac1))write(xyzzyaaac6)etot_config(xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaad1))write(xyzzyaaac6)logdet_config(:,:,xyzzyaaah&
&6)
if(xyzzyaaat1(xyzzyaaae1))write(xyzzyaaac6)fidet_config(:,:,:,xyzzyaaa&
&h6)
if(xyzzyaaat1(xyzzyaaaf1))write(xyzzyaaac6)fi_prod_det_config(:,:,:,:,&
&:,xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaag1))write(xyzzyaaac6)lapdet_config(xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaah1))write(xyzzyaaac6)prod_lapdet_config(:,:,:,:,&
&xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaai1))write(xyzzyaaac6)local_potential_config(xyzz&
&yaaah6)
if(xyzzyaaat1(xyzzyaaaj1))write(xyzzyaaac6)nltot_config(xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaak1))write(xyzzyaaac6)stot_config(xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaal1))write(xyzzyaaac6)wdmc_config(xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaam1))write(xyzzyaaac6)valjas_config(xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaan1))write(xyzzyaaac6)logp_config(xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaao1))write(xyzzyaaac6)twist_config(:,xyzzyaaah6)
enddo
call qmc_barrier
do xyzzyaaag6=1,nnodes-1
xyzzyaaal6=xyzzyaaam6(xyzzyaaag6)
if(xyzzyaaal6>0)then
do xyzzyaaae6=1,xyzzyaaaw1
select case(trim(xyzzyaaba1(xyzzyaaae6)))
case('RELE')
allocate(xyzzyaaau6(3,netot,xyzzyaaal6),stat=xyzzyaaai6)
case('SELE')
allocate(xyzzyaaan6(netot,xyzzyaaal6),stat=xyzzyaaai6)
case('ETOT')
allocate(xyzzyaaaw6(xyzzyaaal6),stat=xyzzyaaai6)
case('LOGDET')
allocate(xyzzyaabi6(nspin,ndet,xyzzyaaal6),stat=xyzzyaaai6)
case('FIDET')
allocate(xyzzyaaba6(3,netot,real1_complex2,xyzzyaaal6),stat=xyzzyaaai6&
&)
case('FI_PROD_DET')
allocate(xyzzyaabb6(3,ndet,nemax,real1_complex2,nspin,xyzzyaaal6),stat&
&=xyzzyaaai6)
case('LAPDET')
allocate(xyzzyaaav6(xyzzyaaal6),stat=xyzzyaaai6)
case('PROD_LAPDET')
allocate(xyzzyaabc6(ndet,nemax,real1_complex2,nspin,xyzzyaaal6),stat=x&
&yzzyaaai6)
case('LOCAL_POTENTIAL')
allocate(xyzzyaaay6(xyzzyaaal6),stat=xyzzyaaai6)
case('NLTOT')
allocate(xyzzyaaaz6(xyzzyaaal6),stat=xyzzyaaai6)
case('STOT')
allocate(xyzzyaaax6(xyzzyaaal6),stat=xyzzyaaai6)
case('WDMC')
allocate(xyzzyaabd6(xyzzyaaal6),stat=xyzzyaaai6)
case('VALJAS')
allocate(xyzzyaabe6(xyzzyaaal6),stat=xyzzyaaai6)
case('LOGP')
allocate(xyzzyaabf6(xyzzyaaal6),stat=xyzzyaaai6)
case('TWIST')
allocate(xyzzyaabg6(5,xyzzyaaal6),stat=xyzzyaaai6)
end select
call check_alloc(xyzzyaaai6,'WRITE_CONFIGS',trim(xyzzyaaba1(xyzzyaaae6&
&))//'_NODE')
enddo
do xyzzyaaae6=1,xyzzyaaaw1
select case(trim(xyzzyaaba1(xyzzyaaae6)))
case('RELE')
call qmpi_recv(xyzzyaaau6,xyzzyaaag6,load_msg,'WRITE_CONFIGS','RELE_NO&
&DE')
case('SELE')
call qmpi_recv(xyzzyaaan6,xyzzyaaag6,load_msg,'WRITE_CONFIGS','SELE_NO&
&DE')
case('ETOT')
call qmpi_recv(xyzzyaaaw6,xyzzyaaag6,load_msg,'WRITE_CONFIGS','ETOT_NO&
&DE')
case('LOGDET')
call qmpi_recv(xyzzyaabi6,xyzzyaaag6,load_msg,'WRITE_CONFIGS','LOGDET_&
&NODE')
case('FIDET')
call qmpi_recv(xyzzyaaba6,xyzzyaaag6,load_msg,'WRITE_CONFIGS','FIDET_N&
&ODE')
case('FI_PROD_DET')
call qmpi_recv(xyzzyaabb6,xyzzyaaag6,load_msg,'WRITE_CONFIGS','FI_PROD&
&_DET_NODE')
case('LAPDET')
call qmpi_recv(xyzzyaaav6,xyzzyaaag6,load_msg,'WRITE_CONFIGS','LAPDET_&
&NODE')
case('PROD_LAPDET')
call qmpi_recv(xyzzyaabc6,xyzzyaaag6,load_msg,'WRITE_CONFIGS','PROD_LA&
&PDET_NODE')
case('LOCAL_POTENTIAL')
call qmpi_recv(xyzzyaaay6,xyzzyaaag6,load_msg,'WRITE_CONFIGS','LOCAL_P&
&OTENTIAL_NODE')
case('NLTOT')
call qmpi_recv(xyzzyaaaz6,xyzzyaaag6,load_msg,'WRITE_CONFIGS','NLTOT_N&
&ODE')
case('STOT')
call qmpi_recv(xyzzyaaax6,xyzzyaaag6,load_msg,'WRITE_CONFIGS','STOT_NO&
&DE')
case('WDMC')
call qmpi_recv(xyzzyaabd6,xyzzyaaag6,load_msg,'WRITE_CONFIGS','WDMC_NO&
&DE')
case('VALJAS')
call qmpi_recv(xyzzyaabe6,xyzzyaaag6,load_msg,'WRITE_CONFIGS','VALJAS_&
&NODE')
case('LOGP')
call qmpi_recv(xyzzyaabf6,xyzzyaaag6,load_msg,'WRITE_CONFIGS','LOGP_NO&
&DE')
case('TWIST')
call qmpi_recv(xyzzyaabg6,xyzzyaaag6,load_msg,'WRITE_CONFIGS','TWIST_N&
&ODE')
end select
call checkmpi(ierror,'recv '//trim(xyzzyaaba1(xyzzyaaae6))//' in write&
&_configs')
enddo
do xyzzyaaah6=1,xyzzyaaal6
if(xyzzyaaat1(xyzzyaaaa1))write(xyzzyaaac6)xyzzyaaau6(:,:,xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaab1))write(xyzzyaaac6)xyzzyaaan6(:,xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaac1))write(xyzzyaaac6)xyzzyaaaw6(xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaad1))write(xyzzyaaac6)xyzzyaabi6(:,:,xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaae1))write(xyzzyaaac6)xyzzyaaba6(:,:,:,xyzzyaaah6&
&)
if(xyzzyaaat1(xyzzyaaaf1))write(xyzzyaaac6)xyzzyaabb6(:,:,:,:,:,xyzzya&
&aah6)
if(xyzzyaaat1(xyzzyaaag1))write(xyzzyaaac6)xyzzyaaav6(xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaah1))write(xyzzyaaac6)xyzzyaabc6(:,:,:,:,xyzzyaaa&
&h6)
if(xyzzyaaat1(xyzzyaaai1))write(xyzzyaaac6)xyzzyaaay6(xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaaj1))write(xyzzyaaac6)xyzzyaaaz6(xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaak1))write(xyzzyaaac6)xyzzyaaax6(xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaal1))write(xyzzyaaac6)xyzzyaabd6(xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaam1))write(xyzzyaaac6)xyzzyaabe6(xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaan1))write(xyzzyaaac6)xyzzyaabf6(xyzzyaaah6)
if(xyzzyaaat1(xyzzyaaao1))write(xyzzyaaac6)xyzzyaabg6(:,xyzzyaaah6)
enddo
if(allocated(xyzzyaaau6))deallocate(xyzzyaaau6)
if(allocated(xyzzyaaan6))deallocate(xyzzyaaan6)
if(allocated(xyzzyaabi6))deallocate(xyzzyaabi6)
if(allocated(xyzzyaaav6))deallocate(xyzzyaaav6)
if(allocated(xyzzyaaaw6))deallocate(xyzzyaaaw6)
if(allocated(xyzzyaaax6))deallocate(xyzzyaaax6)
if(allocated(xyzzyaaay6))deallocate(xyzzyaaay6)
if(allocated(xyzzyaaaz6))deallocate(xyzzyaaaz6)
if(allocated(xyzzyaaba6))deallocate(xyzzyaaba6)
if(allocated(xyzzyaabb6))deallocate(xyzzyaabb6)
if(allocated(xyzzyaabc6))deallocate(xyzzyaabc6)
if(allocated(xyzzyaabd6))deallocate(xyzzyaabd6)
if(allocated(xyzzyaabe6))deallocate(xyzzyaabe6)
if(allocated(xyzzyaabf6))deallocate(xyzzyaabf6)
if(allocated(xyzzyaabg6))deallocate(xyzzyaabg6)
endif
call qmc_barrier
enddo
endif
close(xyzzyaaac6)
open_unit(xyzzyaaac6)=.false.
else
do xyzzyaaaf6=1,xyzzyaaay1
label=xyzzyaabb1(xyzzyaaaf6)
select case(trim(label))
case('VMC_SAVED_STATE')
if(need_shm)then
call shallocate(xyzzyaabh6,3,netot,nnpsmp,stat=xyzzyaaai6)
xyzzyaabh6(:,:,my_smpproc)=rele_vmc_config(:,:)
call shallocate_barrier
if(am_smpmaster)call qmpi_ssend(xyzzyaabh6,0,load_msg,'WRITE_CONFIGS',&
&'rele_vmc_config_all')
call shallocate_barrier
call deshallocate(xyzzyaabh6)
else
allocate(xyzzyaabh6(1,1,1))
call mpi_gather(rele_vmc_config,three_netot,mpi_double_precision,xyzzy&
&aabh6,three_netot,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gather rele_vmc in write_configs')
deallocate(xyzzyaabh6)
endif
if(noncoll_spin)then
if(need_shm)then
call shallocate(xyzzyaaao6,netot,nnpsmp,stat=xyzzyaaai6)
xyzzyaaao6(:,my_smpproc)=sele_vmc_config(:)
call shallocate_barrier
if(am_smpmaster)call qmpi_ssend(xyzzyaaao6,0,load_msg,'WRITE_CONFIGS',&
&'sele_vmc_config_all')
call shallocate_barrier
call deshallocate(xyzzyaaao6)
else
allocate(xyzzyaaao6(1,1))
call mpi_gather(sele_vmc_config,netot,mpi_integer,xyzzyaaao6,netot,mpi&
&_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gather sele_vmc in write_configs')
deallocate(xyzzyaaao6)
endif
endif
case('REBLOCK_DATA')
continue
case('RANDOM')
if(need_shm)then
call shallocate(xyzzyaaap6,25,nnpsmp,stat=xyzzyaaai6)
xyzzyaaap6(:,my_smpproc)=random_state_config(:)
call shallocate_barrier
if(am_smpmaster)call qmpi_ssend(xyzzyaaap6,0,load_msg,'WRITE_CONFIGS',&
&'random_state_config_all')
call shallocate_barrier
call deshallocate(xyzzyaaap6)
else
allocate(xyzzyaaap6(1,1))
call mpi_gather(random_state_config,25,mpi_integer,xyzzyaaap6,25,mpi_i&
&nteger,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gather random_state_config_all in write_configs'&
&)
deallocate(xyzzyaaap6)
endif
if(rngbuf)then
if(need_shm)then
call shallocate(xyzzyaaaq6,nnpsmp,stat=xyzzyaaai6)
xyzzyaaaq6(my_smpproc)=ranx_indx
call shallocate(xyzzyaaar6,nnpsmp,stat=xyzzyaaai6)
xyzzyaaar6(my_smpproc)=ranx_gauss_indx
call shallocate_barrier
if(am_smpmaster)then
call qmpi_ssend(xyzzyaaaq6,0,load_msg,'WRITE_CONFIGS','ranx_indx_all')
call qmpi_ssend(xyzzyaaar6,0,load_msg,'WRITE_CONFIGS','ranx_gauss_indx&
&_all')
endif
call shallocate_barrier
call deshallocate(xyzzyaaaq6)
call deshallocate(xyzzyaaar6)
else
allocate(xyzzyaaaq6(1))
allocate(xyzzyaaar6(1))
call mpi_gather(ranx_indx,1,mpi_integer,xyzzyaaaq6,1,mpi_integer,0,mpi&
&_comm_world,ierror)
call checkmpi(ierror,'gather ranx_indx in write_configs')
call mpi_gather(ranx_gauss_indx,1,mpi_integer,xyzzyaaar6,1,mpi_integer&
&,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gather ranx_gauss_indx in write_configs')
deallocate(xyzzyaaaq6,xyzzyaaar6)
endif
if(need_shm)then
call shallocate(xyzzyaaat6,ranx_max*nnpsmp,stat=xyzzyaaai6)
xyzzyaaaa6=(my_smpproc-1)*ranx_max+1
xyzzyaaat6(xyzzyaaaa6:xyzzyaaaa6+ranx_max-1)=ranx_buffer(:)
call shallocate_barrier
if(am_smpmaster)call qmpi_ssend(xyzzyaaat6,0,load_msg,'WRITE_CONFIGS',&
&'rng_buffer_all')
call shallocate_barrier
call deshallocate(xyzzyaaat6)
else
allocate(xyzzyaaat6(1))
call mpi_gather(ranx_buffer,ranx_max,mpi_real,xyzzyaaat6,ranx_max,mpi_&
&real,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gather ranx_buffer in write_configs')
deallocate(xyzzyaaat6)
endif
if(need_shm)then
call shallocate(xyzzyaaat6,ranx_gauss_max*nnpsmp,stat=xyzzyaaai6)
xyzzyaaaa6=(my_smpproc-1)*ranx_gauss_max+1
xyzzyaaat6(xyzzyaaaa6:xyzzyaaaa6+ranx_gauss_max-1)=ranx_gauss_buffer(:&
&)
call shallocate_barrier
if(am_smpmaster)call qmpi_ssend(xyzzyaaat6,0,load_msg,'WRITE_CONFIGS',&
&'rng_buffer_all <2>')
call shallocate_barrier
call deshallocate(xyzzyaaat6)
else
allocate(xyzzyaaat6(1))
call mpi_gather(ranx_gauss_buffer,ranx_gauss_max,mpi_real,xyzzyaaat6,r&
&anx_gauss_max,mpi_real,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gather ranx_gauss_buffer in write_configs')
deallocate(xyzzyaaat6)
endif
endif
case('GEOMETRY')
continue
end select
enddo
if(xyzzyaaak6>0)then
do xyzzyaaag6=0,nnodes-1
if(my_node==xyzzyaaag6.and.xyzzyaaaj6>0)then
do xyzzyaaae6=1,xyzzyaaaw1
select case(trim(xyzzyaaba1(xyzzyaaae6)))
case('RELE')
call qmpi_ssend(rele_config(:,:,1:xyzzyaaaj6),0,load_msg,'WRITE_CONFIG&
&S','RELE_CONFIG')
case('SELE')
call qmpi_ssend(sele_config(:,1:xyzzyaaaj6),0,load_msg,'WRITE_CONFIGS'&
&,'SELE_CONFIG')
case('ETOT')
call qmpi_ssend(etot_config(1:xyzzyaaaj6),0,load_msg,'WRITE_CONFIGS','&
&ETOT_CONFIG')
case('LOGDET')
call qmpi_ssend(logdet_config(:,:,1:xyzzyaaaj6),0,load_msg,'WRITE_CONF&
&IGS','LOGDET_CONFIG')
case('FIDET')
call qmpi_ssend(fidet_config(:,:,:,1:xyzzyaaaj6),0,load_msg,'WRITE_CON&
&FIGS','FIDET_CONFIG')
case('FI_PROD_DET')
call qmpi_ssend(fi_prod_det_config(:,:,:,:,:,1:xyzzyaaaj6),0,load_msg,&
&'WRITE_CONFIG','FI_PROD_DET_CONFIG')
case('LAPDET')
call qmpi_ssend(lapdet_config(1:xyzzyaaaj6),0,load_msg,'WRITE_CONFIGS'&
&,'LAPDET_CONFIG')
case('PROD_LAPDET')
call qmpi_ssend(prod_lapdet_config(:,:,:,:,1:xyzzyaaaj6),0,load_msg,'W&
&RITE_CONFIG','PROD_LAPDET_CONFIG')
case('LOCAL_POTENTIAL')
call qmpi_ssend(local_potential_config(1:xyzzyaaaj6),0,load_msg,'WRITE&
&_CONFIGS','LOCAL_POTENTIAL_CONFIG')
case('NLTOT')
call qmpi_ssend(nltot_config(1:xyzzyaaaj6),0,load_msg,'WRITE_CONFIGS',&
&'NLTOT_CONFIG')
case('STOT')
call qmpi_ssend(stot_config(1:xyzzyaaaj6),0,load_msg,'WRITE_CONFIGS','&
&STOT_CONFIG')
case('WDMC')
call qmpi_ssend(wdmc_config(1:xyzzyaaaj6),0,load_msg,'WRITE_CONFIGS','&
&WDMC_CONFIG')
case('VALJAS')
call qmpi_ssend(valjas_config(1:xyzzyaaaj6),0,load_msg,'WRITE_CONFIGS'&
&,'VALJAS_CONFIG')
case('LOGP')
call qmpi_ssend(logp_config(1:xyzzyaaaj6),0,load_msg,'WRITE_CONFIGS','&
&LOGP_CONFIG')
case('TWIST')
call qmpi_ssend(twist_config(:,1:xyzzyaaaj6),0,load_msg,'WRITE_CONFIGS&
&','TWIST_CONFIG')
end select
call checkmpi(ierror,'ssend '//trim(xyzzyaaba1(xyzzyaaae6))//' in writ&
&e_configs')
enddo
endif
call qmc_barrier
enddo
endif
endif
call timer('WRITE_CONFIGS',.false.)
end subroutine write_configs
subroutine xyzzyaabj1(io,label)
implicit none
integer,intent(in) :: io
character(*),intent(in) :: label
integer xyzzyaaaa7
character(20) label20
label20=label
write(io,iostat=xyzzyaaaa7)label20
if(xyzzyaaaa7/=0)call errstop('WRITE_LABEL','Problem writing label: '/&
&/trim(label)//'.')
end subroutine
subroutine init_config_accumulation(gen_by,ncfg,cfg_items,xtra_items)
implicit none
integer,intent(in) :: ncfg
character(*),intent(in) :: gen_by
character(20),dimension(:),intent(inout) :: cfg_items,xtra_items
integer xyzzyaaaa8,xyzzyaaab8,xyzzyaaac8,xyzzyaaad8,xyzzyaaae8,xyzzyaa&
&af8,xyzzyaaag8
integer,allocatable :: xyzzyaaah8(:,:)
real(dp),allocatable :: xyzzyaaai8(:,:,:),xyzzyaaaj8(:),xyzzyaaak8(:),&
&xyzzyaaal8(:),xyzzyaaam8(:),xyzzyaaan8(:),xyzzyaaao8(:,:,:,:),xyzzyaa&
&ap8(:,:,:,:,:,:),xyzzyaaaq8(:,:,:,:,:),xyzzyaaar8(:),xyzzyaaas8(:),xy&
&zzyaaat8(:),xyzzyaaau8(:,:)
complex(dp),allocatable :: xyzzyaaav8(:,:,:)
logical xyzzyaaaw8,xyzzyaaax8(xyzzyaaap1)
if(xyzzyaaas1)call errstop_master('INIT_CONFIG_ACCUMULATION','Already &
&initialized. Bug.')
select case(trim(gen_by))
case('VMC','DMC','OPT')
continue
case default
call errstop_master('INIT_CONFIG_ACCUMULATION','Unknown generation met&
&hod '//trim(gen_by)//'.')
end select
if(.not.xyzzyaaar1.or.xyzzyaaav1==0)then
xyzzyaabc1=gen_by
xyzzyaaaw1=size(cfg_items)
allocate(xyzzyaaba1(xyzzyaaaw1),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','config_item')
xyzzyaaba1=cfg_items
xyzzyaaat1=.false.
do xyzzyaaab8=1,xyzzyaaaw1
xyzzyaaaf8=0
do xyzzyaaag8=1,xyzzyaaap1
if(trim(xyzzyaaba1(xyzzyaaab8))==trim(xyzzyaaaq1(xyzzyaaag8)))then
xyzzyaaaf8=xyzzyaaag8
exit
endif
enddo
if(xyzzyaaaf8==0)call errstop_master('INIT_CONFIG_ACCUMULATION',"Unrec&
&ognized item '"//trim(xyzzyaaba1(xyzzyaaab8))//"'.")
if(xyzzyaaat1(xyzzyaaaf8))call errstop_master('INIT_CONFIG_ACCUMULATIO&
&N','Config items repeated.')
xyzzyaaat1(xyzzyaaaf8)=.true.
enddo
xyzzyaaab8=0
do xyzzyaaag8=1,xyzzyaaap1
if(xyzzyaaat1(xyzzyaaag8))then
xyzzyaaab8=xyzzyaaab8+1
xyzzyaaba1(xyzzyaaab8)=xyzzyaaaq1(xyzzyaaag8)
endif
enddo
if(xyzzyaaab8/=xyzzyaaaw1)call errstop_master('INIT_CONFIG_ACCUMULATIO&
&N','Can''t count.')
xyzzyaaau1=ncfg
xyzzyaaax1=0
allocate(xyzzyaaaz1(0:nnodes-1),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','nconfig_node')
call qmpi_allgather(xyzzyaaau1,xyzzyaaaz1(0:nnodes-1),mpi_comm_world,'&
&INIT_CONFIG_ACCUMULATION','nconfig')
xyzzyaaav1=sum(xyzzyaaaz1)
do xyzzyaaab8=1,xyzzyaaaw1
select case(trim(xyzzyaaba1(xyzzyaaab8)))
case('RELE')
allocate(rele_config(3,netot,xyzzyaaau1),stat=xyzzyaaaa8)
if(xyzzyaaaa8==0)rele_config=0.d0
case('SELE')
allocate(sele_config(netot,xyzzyaaau1),stat=xyzzyaaaa8)
if(xyzzyaaaa8==0)sele_config=0
case('ETOT')
allocate(etot_config(xyzzyaaau1),stat=xyzzyaaaa8)
if(xyzzyaaaa8==0)etot_config=0.d0
case('LOGDET')
allocate(logdet_config(nspin,ndet,xyzzyaaau1),stat=xyzzyaaaa8)
if(xyzzyaaaa8==0)logdet_config=cmplx(0.d0,0.d0,dp)
case('FIDET')
allocate(fidet_config(3,netot,real1_complex2,xyzzyaaau1),stat=xyzzyaaa&
&a8)
if(xyzzyaaaa8==0)fidet_config=0.d0
case('FI_PROD_DET')
allocate(fi_prod_det_config(3,ndet,nemax,real1_complex2,nspin,xyzzyaaa&
&u1),stat=xyzzyaaaa8)
if(xyzzyaaaa8==0)fi_prod_det_config=0.d0
case('LAPDET')
allocate(lapdet_config(xyzzyaaau1),stat=xyzzyaaaa8)
if(xyzzyaaaa8==0)lapdet_config=0.d0
case('PROD_LAPDET')
allocate(prod_lapdet_config(ndet,nemax,real1_complex2,nspin,xyzzyaaau1&
&),stat=xyzzyaaaa8)
if(xyzzyaaaa8==0)prod_lapdet_config=0.d0
case('LOCAL_POTENTIAL')
allocate(local_potential_config(xyzzyaaau1),stat=xyzzyaaaa8)
if(xyzzyaaaa8==0)local_potential_config=0.d0
case('NLTOT')
allocate(nltot_config(xyzzyaaau1),stat=xyzzyaaaa8)
if(xyzzyaaaa8==0)nltot_config=0.d0
case('STOT')
allocate(stot_config(xyzzyaaau1),stat=xyzzyaaaa8)
if(xyzzyaaaa8==0)stot_config=0.d0
case('WDMC')
allocate(wdmc_config(xyzzyaaau1),stat=xyzzyaaaa8)
if(xyzzyaaaa8==0)wdmc_config=0.d0
case('VALJAS')
allocate(valjas_config(xyzzyaaau1),stat=xyzzyaaaa8)
if(xyzzyaaaa8==0)valjas_config=0.d0
case('LOGP')
allocate(logp_config(xyzzyaaau1),stat=xyzzyaaaa8)
if(xyzzyaaaa8==0)logp_config=0.d0
case('TWIST')
allocate(twist_config(5,xyzzyaaau1),stat=xyzzyaaaa8)
if(xyzzyaaaa8==0)twist_config=0.d0
case default
call errstop_master('INIT_CONFIG_ACCUMULATION','Unknown config item ty&
&pe '//trim(xyzzyaaba1(xyzzyaaab8))//'.')
end select
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION',trim(xyzzyaaba1&
&(xyzzyaaab8)))
enddo
else
if(trim(xyzzyaabc1)/=trim(gen_by))call errstop_master('INIT_CONFIG_ACC&
&UMULATION','Trying to continue '//trim(xyzzyaabc1)//' config accumula&
&tion using '//trim(gen_by)//'.')
if(xyzzyaaaw1/=size(cfg_items))call errstop_master('INIT_CONFIG_ACCUMU&
&LATION','Trying to continue config accumulation, but different config&
& items are present.')
xyzzyaaax8=.false.
do xyzzyaaab8=1,xyzzyaaaw1
xyzzyaaaf8=0
do xyzzyaaag8=1,xyzzyaaap1
if(trim(cfg_items(xyzzyaaab8))==trim(xyzzyaaaq1(xyzzyaaag8)))then
xyzzyaaaf8=xyzzyaaag8
exit
endif
enddo
if(xyzzyaaaf8==0)call errstop_master('INIT_CONFIG_ACCUMULATION',"Unrec&
&ognized item '"//trim(cfg_items(xyzzyaaab8))//"'.")
if(xyzzyaaax8(xyzzyaaaf8))call errstop_master('INIT_CONFIG_ACCUMULATIO&
&N','Config items repeated.')
xyzzyaaax8(xyzzyaaaf8)=.true.
enddo
xyzzyaaab8=0
do xyzzyaaag8=1,xyzzyaaap1
if(xyzzyaaax8(xyzzyaaag8))then
xyzzyaaab8=xyzzyaaab8+1
cfg_items(xyzzyaaab8)=xyzzyaaaq1(xyzzyaaag8)
endif
enddo
if(xyzzyaaab8/=xyzzyaaaw1)call errstop_master('INIT_CONFIG_ACCUMULATIO&
&N','Can''t count.')
if(any(cfg_items/=xyzzyaaba1))call errstop_master('INIT_CONFIG_ACCUMUL&
&ATION','Trying to continue config accumulation, but different config &
&items are present.')
xyzzyaaad8=xyzzyaaau1
xyzzyaaau1=xyzzyaaau1+ncfg
xyzzyaaax1=xyzzyaaad8
call mpi_gather(xyzzyaaau1,1,mpi_integer,xyzzyaaaz1,1,mpi_integer,0,mp&
&i_comm_world,ierror)
call checkmpi(ierror,'gather nconfig in init_config_accumulation')
call mpi_bcast(xyzzyaaaz1(0:nnodes-1),nnodes,mpi_integer,0,mpi_comm_wo&
&rld,ierror)
call checkmpi(ierror,'bcast nconfig_node in init_config_accumulation')
xyzzyaaav1=sum(xyzzyaaaz1)
if(any(xyzzyaaaz1(:)/=xyzzyaaaz1(0)))call errstop_master('INIT_CONFIG_&
&ACCUMULATION','Trying to continue config accumulation with a differen&
&t number of configs per node, which is unsupported.')
do xyzzyaaab8=1,xyzzyaaaw1
select case(trim(xyzzyaaba1(xyzzyaaab8)))
case('RELE')
allocate(xyzzyaaai8(3,netot,xyzzyaaad8),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','RELE_PREV')
call dcopy(three_netot*xyzzyaaad8,rele_config(1,1,1),1,xyzzyaaai8(1,1,&
&1),1)
deallocate(rele_config)
allocate(rele_config(3,netot,xyzzyaaau1),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','RELE_CONFIG')
call dcopy(three_netot*xyzzyaaad8,xyzzyaaai8(1,1,1),1,rele_config(1,1,&
&1),1)
deallocate(xyzzyaaai8)
case('SELE')
allocate(xyzzyaaah8(netot,xyzzyaaad8),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','SELE_PREV')
xyzzyaaah8=sele_config
allocate(sele_config(netot,xyzzyaaau1),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','SELE_CONFIG')
sele_config(:,1:xyzzyaaad8)=xyzzyaaah8(:,1:xyzzyaaad8)
deallocate(xyzzyaaah8)
case('ETOT')
allocate(xyzzyaaak8(xyzzyaaad8),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','ETOT_PREV')
call dcopy(xyzzyaaad8,etot_config(1),1,xyzzyaaak8(1),1)
deallocate(etot_config)
allocate(etot_config(xyzzyaaau1),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','ETOT_CONFIG')
call dcopy(xyzzyaaad8,xyzzyaaak8(1),1,etot_config(1),1)
deallocate(xyzzyaaak8)
case('LOGDET')
allocate(xyzzyaaav8(nspin,ndet,xyzzyaaad8),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','LOGDET_PREV')
call zcopy(size_det*xyzzyaaad8,logdet_config(1,1,1),1,xyzzyaaav8(1,1,1&
&),1)
deallocate(logdet_config)
allocate(logdet_config(nspin,ndet,xyzzyaaau1),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','LOGDET_CONFIG'&
&)
call zcopy(size_det*xyzzyaaad8,xyzzyaaav8(1,1,1),1,logdet_config(1,1,1&
&),1)
deallocate(xyzzyaaav8)
case('FIDET')
allocate(xyzzyaaao8(3,netot,real1_complex2,xyzzyaaad8),stat=xyzzyaaaa8&
&)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','FIDET_PREV')
call dcopy(size_fidet*xyzzyaaad8,fidet_config(1,1,1,1),1,xyzzyaaao8(1,&
&1,1,1),1)
deallocate(fidet_config)
allocate(fidet_config(3,netot,real1_complex2,xyzzyaaau1),stat=xyzzyaaa&
&a8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','FIDET_CONFIG')
call dcopy(size_fidet*xyzzyaaad8,xyzzyaaao8(1,1,1,1),1,fidet_config(1,&
&1,1,1),1)
deallocate(xyzzyaaao8)
case('FI_PROD_DET')
allocate(xyzzyaaap8(3,ndet,nemax,real1_complex2,nspin,xyzzyaaad8),stat&
&=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','FI_PROD_DET_PR&
&EV')
call dcopy(size_fi_prod_det*xyzzyaaad8,fi_prod_det_config(1,1,1,1,1,1)&
&,1,xyzzyaaap8(1,1,1,1,1,1),1)
deallocate(fi_prod_det_config)
allocate(fi_prod_det_config(3,ndet,nemax,real1_complex2,nspin,xyzzyaaa&
&u1),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','FI_PROD_DET_CO&
&NFIG')
call dcopy(size_fi_prod_det*xyzzyaaad8,xyzzyaaap8(1,1,1,1,1,1),1,fi_pr&
&od_det_config(1,1,1,1,1,1),1)
deallocate(xyzzyaaap8)
case('LAPDET')
allocate(xyzzyaaaj8(xyzzyaaad8),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','LAPDET_PREV')
call dcopy(xyzzyaaad8,lapdet_config(1),1,xyzzyaaaj8(1),1)
deallocate(lapdet_config)
allocate(lapdet_config(xyzzyaaau1),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','LAPDET_CONFIG'&
&)
call dcopy(xyzzyaaad8,xyzzyaaaj8(1),1,lapdet_config(1),1)
deallocate(xyzzyaaaj8)
case('PROD_LAPDET')
allocate(xyzzyaaaq8(ndet,nemax,real1_complex2,nspin,xyzzyaaad8),stat=x&
&yzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','PROD_LAPDET_PR&
&EV')
call dcopy(size_prod_lapdet*xyzzyaaad8,prod_lapdet_config(1,1,1,1,1),1&
&,xyzzyaaaq8(1,1,1,1,1),1)
deallocate(prod_lapdet_config)
allocate(prod_lapdet_config(ndet,nemax,real1_complex2,nspin,xyzzyaaau1&
&),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','PROD_LAPDET_CO&
&NFIG')
call dcopy(size_prod_lapdet*xyzzyaaad8,xyzzyaaaq8(1,1,1,1,1),1,prod_la&
&pdet_config(1,1,1,1,1),1)
deallocate(xyzzyaaaq8)
case('LOCAL_POTENTIAL')
allocate(xyzzyaaam8(xyzzyaaad8),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','LOCAL_POTENTIA&
&L_PREV')
call dcopy(xyzzyaaad8,local_potential_config(1),1,xyzzyaaam8(1),1)
deallocate(local_potential_config)
allocate(local_potential_config(xyzzyaaau1),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','LOCAL_POTENTIA&
&L_CONFIG')
call dcopy(xyzzyaaad8,xyzzyaaam8(1),1,local_potential_config(1),1)
deallocate(xyzzyaaam8)
case('NLTOT')
allocate(xyzzyaaan8(xyzzyaaad8),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','NLTOT_PREV')
call dcopy(xyzzyaaad8,nltot_config(1),1,xyzzyaaan8(1),1)
deallocate(nltot_config)
allocate(nltot_config(xyzzyaaau1),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','NLTOT_CONFIG')
call dcopy(xyzzyaaad8,xyzzyaaan8(1),1,nltot_config(1),1)
deallocate(xyzzyaaan8)
case('STOT')
allocate(xyzzyaaal8(xyzzyaaad8),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','STOT_PREV')
call dcopy(xyzzyaaad8,stot_config(1),1,xyzzyaaal8(1),1)
deallocate(stot_config)
allocate(stot_config(xyzzyaaau1),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','STOT_CONFIG')
call dcopy(xyzzyaaad8,xyzzyaaal8(1),1,stot_config(1),1)
deallocate(xyzzyaaal8)
case('WDMC')
allocate(xyzzyaaar8(xyzzyaaad8),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','WDMC_PREV')
call dcopy(xyzzyaaad8,wdmc_config(1),1,xyzzyaaar8(1),1)
deallocate(wdmc_config)
allocate(wdmc_config(xyzzyaaau1),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','WDMC_CONFIG')
call dcopy(xyzzyaaad8,xyzzyaaar8(1),1,wdmc_config(1),1)
deallocate(xyzzyaaar8)
case('VALJAS')
allocate(xyzzyaaas8(xyzzyaaad8),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','VALJAS_PREV')
call dcopy(xyzzyaaad8,valjas_config(1),1,xyzzyaaas8(1),1)
deallocate(valjas_config)
allocate(valjas_config(xyzzyaaau1),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','VALJAS_CONFIG'&
&)
valjas_config=0.d0
call dcopy(xyzzyaaad8,xyzzyaaas8(1),1,valjas_config(1),1)
deallocate(xyzzyaaas8)
case('LOGP')
allocate(xyzzyaaat8(xyzzyaaad8),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','LOGP_PREV')
call dcopy(xyzzyaaad8,logp_config(1),1,xyzzyaaat8(1),1)
deallocate(logp_config)
allocate(logp_config(xyzzyaaau1),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','LOGP_CONFIG')
logp_config=0.d0
call dcopy(xyzzyaaad8,xyzzyaaat8(1),1,logp_config(1),1)
deallocate(xyzzyaaat8)
case('TWIST')
allocate(xyzzyaaau8(5,xyzzyaaad8),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','TWIST_PREV')
call dcopy(5*xyzzyaaad8,twist_config(1,1),1,xyzzyaaau8(1,1),1)
deallocate(twist_config)
allocate(twist_config(5,xyzzyaaau1),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','TWIST_CONFIG')
twist_config=0.d0
call dcopy(5*xyzzyaaad8,xyzzyaaau8(1,1),1,twist_config(1,1),1)
deallocate(xyzzyaaau8)
case default
call errstop_master('INIT_CONFIG_ACCUMULATION','Unknown config item ty&
&pe '//trim(xyzzyaaba1(xyzzyaaab8))//'.')
end select
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION',trim(xyzzyaaba1&
&(xyzzyaaab8)))
enddo
endif
if(.not.xyzzyaaar1)then
xyzzyaaay1=size(xtra_items)
allocate(xyzzyaabb1(xyzzyaaay1),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','extra_item')
xyzzyaabb1=xtra_items
do xyzzyaaac8=1,xyzzyaaay1
select case(trim(xyzzyaabb1(xyzzyaaac8)))
case('VMC_SAVED_STATE')
allocate(rele_vmc_config(3,netot),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','RELE_VMC_CONFI&
&G')
allocate(sele_vmc_config(netot),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','SELE_VMC_CONFI&
&G')
allocate(dtvmc_array_config(no_difftypes),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','DTVMC_ARRAY_CO&
&NFIG')
case('FINAL_VMC_RESULT')
continue
case('DMC_SAVED_STATE')
if(tpdmc_config>0)then
allocate(log_pi_wt_array_config(0:tpdmc_config-1),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','LOG_PI_WT_ARRA&
&Y_CONFIG')
endif
allocate(log_pi_wt_array2_config(0:tpdmc_config),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','LOG_PI_WT_ARRA&
&Y2_CONFIG')
case('REBLOCK_DATA')
continue
case('RANDOM')
continue
case('GEOMETRY')
continue
case default
call errstop_master('INIT_CONFIG_ACCUMULATION','Unknown extra item typ&
&e '//trim(xyzzyaabb1(xyzzyaaac8))//'.')
end select
enddo
else
if(xyzzyaaay1>0)then
do xyzzyaaac8=1,xyzzyaaay1
xyzzyaaaw8=.false.
do xyzzyaaae8=1,size(xtra_items)
if(trim(xtra_items(xyzzyaaae8))==trim(xyzzyaabb1(xyzzyaaac8)))then
xyzzyaaaw8=.true.
exit
endif
enddo
if(xyzzyaaaw8)cycle
select case(trim(xyzzyaabb1(xyzzyaaac8)))
case('VMC_SAVED_STATE')
if(allocated(rele_vmc_config))deallocate(rele_vmc_config)
if(allocated(sele_vmc_config))deallocate(sele_vmc_config)
if(allocated(dtvmc_array_config))deallocate(dtvmc_array_config)
case('DMC_SAVED_STATE')
if(allocated(log_pi_wt_array_config))deallocate(log_pi_wt_array_config&
&)
if(allocated(log_pi_wt_array2_config))deallocate(log_pi_wt_array2_conf&
&ig)
end select
enddo
endif
do xyzzyaaae8=1,size(xtra_items)
xyzzyaaaw8=.false.
do xyzzyaaac8=1,xyzzyaaay1
if(trim(xtra_items(xyzzyaaae8))==trim(xyzzyaabb1(xyzzyaaac8)))then
xyzzyaaaw8=.true.
exit
endif
enddo
if(xyzzyaaaw8)cycle
select case(trim(xyzzyaabb1(xyzzyaaac8)))
case('VMC_SAVED_STATE')
allocate(rele_vmc_config(3,netot),sele_vmc_config(netot),dtvmc_array_c&
&onfig(no_difftypes),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','VMC_SAVED_STAT&
&E arrays')
case('DMC_SAVED_STATE')
if(tpdmc_config>0)then
allocate(log_pi_wt_array_config(0:tpdmc_config-1),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','LOG_PI_WT_ARRA&
&Y_CONFIG')
endif
allocate(log_pi_wt_array2_config(0:tpdmc_config),stat=xyzzyaaaa8)
call check_alloc(xyzzyaaaa8,'INIT_CONFIG_ACCUMULATION','LOG_PI_WT_ARRA&
&Y2_CONFIG')
end select
enddo
if(xyzzyaaay1>0)deallocate(xyzzyaabb1)
xyzzyaaay1=size(xtra_items)
allocate(xyzzyaabb1(xyzzyaaay1))
xyzzyaabb1=xtra_items
endif
xyzzyaaar1=.false.
xyzzyaaas1=.true.
end subroutine init_config_accumulation
subroutine add_config(modify,rele,sele,etot,logdet,fidet,fi_prod_det,l&
&apdet,prod_lapdet,local_potential,nltot,stot,wdmc,valjas,logp,twist)
implicit none
integer,intent(inout),optional :: sele(netot)
real(dp),intent(inout),optional :: rele(3,netot),etot,fidet(3,netot,re&
&al1_complex2),fi_prod_det(3,ndet,nemax,real1_complex2,nspin),lapdet,p&
&rod_lapdet(ndet,nemax,real1_complex2,nspin),local_potential,nltot,sto&
&t,wdmc,valjas,logp,twist(5)
complex(dp),intent(inout),optional :: logdet(nspin,ndet)
logical,intent(in),optional :: modify
integer xyzzyaaaa9
call timer('ADD_CONFIG',.true.)
if(.not.xyzzyaaas1)call errstop_master('ADD_CONFIG','Have not called I&
&NIT_CONFIG_ACCUMULATION.')
if(present(modify))then
if(.not.modify)xyzzyaaax1=xyzzyaaax1+1
else
xyzzyaaax1=xyzzyaaax1+1
endif
xyzzyaaaa9=xyzzyaaax1
if(xyzzyaaaa9>xyzzyaaau1)call errstop('ADD_CONFIG','One config too man&
&y.')
if(xyzzyaaat1(xyzzyaaaa1).and.present(rele))call dcopy(three_netot,rel&
&e(1,1),1,rele_config(1,1,xyzzyaaaa9),1)
if(xyzzyaaat1(xyzzyaaab1).and.present(sele))sele_config(:,xyzzyaaaa9)=&
&sele
if(xyzzyaaat1(xyzzyaaac1).and.present(etot))etot_config(xyzzyaaaa9)=et&
&ot
if(xyzzyaaat1(xyzzyaaad1).and.present(logdet))call zcopy(size_det,logd&
&et(1,1),1,logdet_config(1,1,xyzzyaaaa9),1)
if(xyzzyaaat1(xyzzyaaae1).and.present(fidet))call dcopy(size_fidet,fid&
&et(1,1,1),1,fidet_config(1,1,1,xyzzyaaaa9),1)
if(xyzzyaaat1(xyzzyaaaf1).and.present(fi_prod_det))call dcopy(size_fi_&
&prod_det,fi_prod_det(1,1,1,1,1),1,fi_prod_det_config(1,1,1,1,1,xyzzya&
&aaa9),1)
if(xyzzyaaat1(xyzzyaaag1).and.present(lapdet))lapdet_config(xyzzyaaaa9&
&)=lapdet
if(xyzzyaaat1(xyzzyaaah1).and.present(prod_lapdet))call dcopy(size_pro&
&d_lapdet,prod_lapdet(1,1,1,1),1,prod_lapdet_config(1,1,1,1,xyzzyaaaa9&
&),1)
if(xyzzyaaat1(xyzzyaaai1).and.present(local_potential))local_potential&
&_config(xyzzyaaaa9)=local_potential
if(xyzzyaaat1(xyzzyaaaj1).and.present(nltot))nltot_config(xyzzyaaaa9)=&
&nltot
if(xyzzyaaat1(xyzzyaaak1).and.present(stot))stot_config(xyzzyaaaa9)=st&
&ot
if(xyzzyaaat1(xyzzyaaal1).and.present(wdmc))wdmc_config(xyzzyaaaa9)=wd&
&mc
if(xyzzyaaat1(xyzzyaaam1).and.present(valjas))valjas_config(xyzzyaaaa9&
&)=valjas
if(xyzzyaaat1(xyzzyaaan1).and.present(logp))logp_config(xyzzyaaaa9)=lo&
&gp
if(xyzzyaaat1(xyzzyaaao1).and.present(twist))twist_config(:,xyzzyaaaa9&
&)=twist(:)
call timer('ADD_CONFIG',.false.)
end subroutine add_config
subroutine delete_config_vmc(n)
implicit none
integer,intent(in) :: n
if(n<0)return
xyzzyaaax1=xyzzyaaax1-n
end subroutine delete_config_vmc
subroutine end_config_accumulation(exceeded)
implicit none
logical,intent(in) :: exceeded
if(xyzzyaaax1/=xyzzyaaau1.and..not.exceeded)call errstop('END_CONFIG_A&
&CCUMULATION','Have not reached the total number of configs yet ('//tr&
&im(i2s(xyzzyaaax1))//'/'//trim(i2s(xyzzyaaau1))//').')
xyzzyaaar1=.true.
xyzzyaaas1=.false.
end subroutine end_config_accumulation
subroutine dismantle_configs
implicit none
if(xyzzyaaas1)call errstop('DISMANTLE_CONFIGS','Inappropriate call: st&
&ill doing accumulation. Bug.')
if(.not.xyzzyaaar1)call errstop('DISMANTLE_CONFIGS','Inappropriate cal&
&l: no data in memory. Bug.')
if(allocated(rele_config))deallocate(rele_config)
if(allocated(sele_config))deallocate(sele_config)
if(allocated(logdet_config))deallocate(logdet_config)
if(allocated(lapdet_config))deallocate(lapdet_config)
if(allocated(etot_config))deallocate(etot_config)
if(allocated(stot_config))deallocate(stot_config)
if(allocated(local_potential_config))deallocate(local_potential_config&
&)
if(allocated(nltot_config))deallocate(nltot_config)
if(allocated(fidet_config))deallocate(fidet_config)
if(allocated(fi_prod_det_config))deallocate(fi_prod_det_config)
if(allocated(prod_lapdet_config))deallocate(prod_lapdet_config)
if(allocated(wdmc_config))deallocate(wdmc_config)
if(allocated(valjas_config))deallocate(valjas_config)
if(allocated(logp_config))deallocate(logp_config)
if(allocated(twist_config))deallocate(twist_config)
if(allocated(xyzzyaaba1))deallocate(xyzzyaaba1)
if(allocated(xyzzyaaaz1))deallocate(xyzzyaaaz1)
xyzzyaaau1=0
xyzzyaaav1=0
xyzzyaaaw1=0
xyzzyaaat1=.false.
if(allocated(rele_vmc_config))deallocate(rele_vmc_config)
if(allocated(sele_vmc_config))deallocate(sele_vmc_config)
if(allocated(dtvmc_array_config))deallocate(dtvmc_array_config)
if(allocated(log_pi_wt_array_config))deallocate(log_pi_wt_array_config&
&)
if(allocated(log_pi_wt_array2_config))deallocate(log_pi_wt_array2_conf&
&ig)
if(allocated(xyzzyaabb1))deallocate(xyzzyaabb1)
random_state_config=0
xyzzyaaay1=0
xyzzyaaar1=.false.
end subroutine dismantle_configs
subroutine shift_config_files
implicit none
call copy_config_file(con_out,con_in,move=.true.)
end subroutine shift_config_files
subroutine backup_config_file
implicit none
call copy_config_file(con_out,con_backup,move=.false.)
end subroutine backup_config_file
subroutine copy_config_file(con_from,con_to,move)
implicit none
logical,intent(in),optional :: move
character(160),intent(in) :: con_from,con_to
integer xyzzyaaaa15,xyzzyaaab15,xyzzyaaac15,xyzzyaaad15,xyzzyaaae15,xy&
&zzyaaaf15,xyzzyaaag15,xyzzyaaah15,xyzzyaaai15
logical xyzzyaaaj15,xyzzyaaak15,xyzzyaaal15(xyzzyaaap1)
character(20) label
integer xyzzyaaam15,xyzzyaaan15,xyzzyaaao15,xyzzyaaap15,xyzzyaaaq15,xy&
&zzyaaar15,xyzzyaaas15,xyzzyaaat15,xyzzyaaau15,xyzzyaaav15,xyzzyaaaw15&
&,xyzzyaaax15,xyzzyaaay15,xyzzyaaaz15,xyzzyaaba15,xyzzyaabb15,xyzzyaab&
&c15,xyzzyaabd15,xyzzyaabe15(25),xyzzyaabf15,xyzzyaabg15
integer,allocatable :: xyzzyaabh15(:),xyzzyaabi15(:),xyzzyaabj15(:)
real(sp) xyzzyaabk15
real(dp) xyzzyaabl15,xyzzyaabm15,xyzzyaabn15,xyzzyaabo15,xyzzyaabp15,x&
&yzzyaabq15,xyzzyaabr15,xyzzyaabs15,xyzzyaabt15(5),xyzzyaabu15,xyzzyaa&
&bv15,xyzzyaabw15,xyzzyaabx15,xyzzyaaby15,xyzzyaabz15,xyzzyaaca15,xyzz&
&yaacb15,xyzzyaacc15,xyzzyaacd15,xyzzyaace15(13),xyzzyaacf15,xyzzyaacg&
&15,xyzzyaach15,xyzzyaaci15,xyzzyaacj15,xyzzyaack15,xyzzyaacl15,xyzzya&
&acm15,xyzzyaacn15
real(dp),allocatable :: xyzzyaaco15(:),xyzzyaacp15(:,:),xyzzyaacq15(:,&
&:),xyzzyaacr15(:,:,:),xyzzyaacs15(:,:,:,:,:),xyzzyaact15(:,:,:,:),xyz&
&zyaacu15(:),xyzzyaacv15(:)
complex(dp),allocatable :: xyzzyaacw15(:,:)
character(20) gen_by_in,interaction_in,atom_basis_type_in
character(20),allocatable :: xyzzyaacx15(:)
integer xyzzyaacy15,xyzzyaacz15,xyzzyaada15,xyzzyaadb15,xyzzyaadc15
integer,allocatable :: xyzzyaadd15(:),xyzzyaade15(:)
real(dp) xyzzyaadf15,xyzzyaadg15
real(dp),allocatable :: xyzzyaadh15(:),xyzzyaadi15(:),xyzzyaadj15(:,:)&
&,xyzzyaadk15(:,:),xyzzyaadl15(:),xyzzyaadm15(:),xyzzyaadn15(:),xyzzya&
&ado15(:),xyzzyaadp15(:),xyzzyaadq15(:,:)
integer xyzzyaadr15,xyzzyaads15
real,allocatable :: xyzzyaadt15(:),xyzzyaadu15(:)
if(am_master)then
call timer('COPY_CONFIG_FILE',.true.)
xyzzyaaaj15=.false.
if(present(move))xyzzyaaaj15=move
gen_by_in='NONE'
xyzzyaaaq15=0
xyzzyaaaw15=-1
xyzzyaaam15=-1
xyzzyaaan15=-1
xyzzyaabd15=0
xyzzyaaay15=0
xyzzyaaav15=1
xyzzyaaaz15=0
xyzzyaaat15=-1
inquire(file=trim(con_from),exist=xyzzyaaak15)
if(.not.xyzzyaaak15)then
inquire(file=trim(con_to),exist=xyzzyaaak15)
if(.not.xyzzyaaak15)then
call qmc_barrier
call timer('COPY_CONFIG_FILE',.false.)
return
endif
call open_units(xyzzyaaab15,xyzzyaaac15)
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Could not find free&
& i/o unit number for deleting existing file.')
open(unit=xyzzyaaab15,file=trim(con_to),form='unformatted',status='old&
&',iostat=xyzzyaaac15)
if(xyzzyaaac15==0)close(xyzzyaaab15,status='delete')
open_unit(xyzzyaaab15)=.false.
call qmc_barrier
call timer('COPY_CONFIG_FILE',.false.)
return
endif
call open_units(xyzzyaaaa15,xyzzyaaac15)
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Could not find free&
& i/o unit number (FROM).')
open(unit=xyzzyaaaa15,file=trim(con_from),form='unformatted',status='o&
&ld',iostat=xyzzyaaac15)
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem opening '//&
&trim(con_from)//' for copy.')
call open_units(xyzzyaaab15,xyzzyaaac15)
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Could not find free&
& i/o unit number (TO).')
open(unit=xyzzyaaab15,file=trim(con_to),form='unformatted',status='rep&
&lace',iostat=xyzzyaaac15)
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem opening '//&
&trim(con_to)//'.')
read(xyzzyaaaa15,iostat=xyzzyaaac15)label
if(xyzzyaaac15<0)call errstop('COPY_CONFIG_FILE','Config file appears &
&to be empty?')
if(xyzzyaaac15>0)call errstop('COPY_CONFIG_FILE','Problem reading firs&
&t label in config file. Perhaps this config file is in old format? Tr&
&y using the update_config utility.')
if(trim(label)/='INFO')call errstop('COPY_CONFIG_FILE','First label in&
& config file is wrong. Perhaps this config file is in old format? Try&
& using the update_config utility.')
rewind(xyzzyaaaa15)
do
read(xyzzyaaaa15,iostat=xyzzyaaac15)label
if(xyzzyaaac15<0)then
if(xyzzyaaaw15<0)then
xyzzyaaaw15=0
exit
endif
call errstop('COPY_CONFIG_FILE','End-of-file reached while reading sec&
&tion name.')
endif
if(xyzzyaaac15>0)call errstop('COPY_CONFIG_FILE','Problem reading sect&
&ion name.')
write(xyzzyaaab15)label
select case(trim(label))
case('INFO')
do
read(xyzzyaaaa15,iostat=xyzzyaaac15)label
if(xyzzyaaac15<0)call errstop('COPY_CONFIG_FILE','End-of-file reached &
&while reading info label.')
if(xyzzyaaac15>0)call errstop('COPY_CONFIG_FILE','Problem reading info&
& label.')
write(xyzzyaaab15)label
if(trim(label)=='END INFO')goto 80
select case(trim(label))
case('NEXTRA')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaaas15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading NEX&
&TRA.')
write(xyzzyaaab15)xyzzyaaas15
case('GEN_BY')
read(xyzzyaaaa15,iostat=xyzzyaaac15)gen_by_in
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading GEN&
&_BY.')
write(xyzzyaaab15)gen_by_in
case('NSPIN')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaaam15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading NSP&
&IN.')
write(xyzzyaaab15)xyzzyaaam15
case('NELE')
if(xyzzyaaam15<0)call errstop('COPY_CONFIG_FILE','NELE found before NS&
&PIN in config file. This is a bug.')
allocate(xyzzyaabh15(xyzzyaaam15),stat=xyzzyaaad15)
call check_alloc(xyzzyaaad15,'COPY_CONFIG_FILE','NELE_IN')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabh15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading NEL&
&E.')
write(xyzzyaaab15)xyzzyaabh15
xyzzyaaao15=maxval(xyzzyaabh15)
xyzzyaaap15=sum(xyzzyaabh15)
deallocate(xyzzyaabh15)
case('NDET')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaaan15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading NDE&
&T.')
write(xyzzyaaab15)xyzzyaaan15
case('COMPLEX_WF')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaaay15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading COM&
&PLEX_WF.')
if(int2log(xyzzyaaay15))xyzzyaaav15=2
write(xyzzyaaab15)xyzzyaaay15
case('NONCOLL_SPIN')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaaaz15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading NON&
&COLL_SPIN.')
write(xyzzyaaab15)xyzzyaaaz15
case('INTERACTION')
read(xyzzyaaaa15,iostat=xyzzyaaac15)interaction_in
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading INT&
&ERACTION.')
write(xyzzyaaab15)interaction_in
case('ATOM_BASIS_TYPE')
read(xyzzyaaaa15,iostat=xyzzyaaac15)atom_basis_type_in
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading ATO&
&M_BASIS_TYPE.')
write(xyzzyaaab15)atom_basis_type_in
case('TOTAL_TIME')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabk15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading TOT&
&AL_TIME_PREVIOUS.')
write(xyzzyaaab15)xyzzyaabk15
case default
call errstop('COPY_CONFIG_FILE','Found unknown info '//trim(label)//'.&
&')
end select
enddo
80 continue
case('DEFINE_CONFIGS')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaaaw15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading num&
&ber of configs.')
write(xyzzyaaab15)xyzzyaaaw15
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaaaq15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading num&
&ber of items per config.')
write(xyzzyaaab15)xyzzyaaaq15
allocate(xyzzyaacx15(xyzzyaaaq15),stat=xyzzyaaad15)
call check_alloc(xyzzyaaad15,'COPY_CONFIG_FILE','config_item')
xyzzyaaal15=.false.
do xyzzyaaar15=1,xyzzyaaaq15
read(xyzzyaaaa15,iostat=xyzzyaaac15)label
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading con&
&fig content label.')
write(xyzzyaaab15)label
xyzzyaacx15(xyzzyaaar15)=trim(label)
xyzzyaaaf15=0
do xyzzyaaag15=1,xyzzyaaap1
if(trim(xyzzyaacx15(xyzzyaaar15))==trim(xyzzyaaaq1(xyzzyaaag15)))then
xyzzyaaaf15=xyzzyaaag15
goto 90
endif
enddo
90   continue
if(xyzzyaaaf15==0)call errstop_master('COPY_CONFIG_FILE',"Unrecognized&
& item '"//trim(label)//"'.")
if(any(xyzzyaaal15(xyzzyaaaf15:)))call errstop_master('COPY_CONFIG_FIL&
&E','Config items in wrong order or repeated.')
xyzzyaaal15(xyzzyaaaf15)=.true.
enddo
read(xyzzyaaaa15,iostat=xyzzyaaac15)label
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE',"Problem reading 'EN&
&D DEFINE_CONFIGS'.")
if(trim(label)/='END DEFINE_CONFIGS')call errstop('COPY_CONFIG_FILE',"&
&Expected to find 'END DEFINE_CONFIGS' but didn't.")
write(xyzzyaaab15)label
case('VMC_SAVED_STATE')
do
read(xyzzyaaaa15,iostat=xyzzyaaac15)label
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE',"Problem reading lab&
&el in VMC_SAVED_STATE section in config file.")
write(xyzzyaaab15)label
select case(trim(label))
case('NO_DIFFTYPES')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaaat15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading NO_&
&DIFFTYPES.')
write(xyzzyaaab15)xyzzyaaat15
case('DTVMC')
if(xyzzyaaat15<1)call errstop('COPY_CONFIG_FILE','DTVMC found before N&
&O_DIFFTYPES in config file.')
allocate(xyzzyaaco15(xyzzyaaat15),stat=xyzzyaaad15)
call check_alloc(xyzzyaaad15,'COPY_CONFIG_FILE','DTVMC_ARRAY_IN')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaaco15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading DTV&
&MC.')
write(xyzzyaaab15)xyzzyaaco15
deallocate(xyzzyaaco15)
case('VMC_STEPS')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaaau15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading VMC&
&_STEPS.')
write(xyzzyaaab15)xyzzyaaau15
case('NNODES')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabd15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading NNO&
&DES.')
write(xyzzyaaab15)xyzzyaabd15
case('RELE')
if(xyzzyaaap15<1)call errstop('COPY_CONFIG_FILE','NETOT undefined when&
& reading RELE VMC_SAVED_STATE block. This is a bug.')
if(xyzzyaabd15<1)call errstop('COPY_CONFIG_FILE','NNODES undefined whe&
&n reading RELE VMC_SAVED_STATE block. This is a bug.')
allocate(xyzzyaacp15(3,xyzzyaaap15),stat=xyzzyaaad15)
call check_alloc(xyzzyaaad15,'COPY_CONFIG_FILE','RELE_VMC_IN')
do xyzzyaaae15=1,xyzzyaabd15
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaacp15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading REL&
&E_VMC.')
write(xyzzyaaab15)xyzzyaacp15
enddo
deallocate(xyzzyaacp15)
case('SELE')
if(xyzzyaaap15<1)call errstop('COPY_CONFIG_FILE','NETOT undefined when&
& reading SELE in VMC_SAVED_STATE block. This is a bug.')
if(xyzzyaabd15<1)call errstop('COPY_CONFIG_FILE','NNODES undefined whe&
&n reading RELE VMC_SAVED_STATE block. This is a bug.')
allocate(xyzzyaabi15(xyzzyaaap15),stat=xyzzyaaad15)
call check_alloc(xyzzyaaad15,'COPY_CONFIG_FILE','SELE_VMC_CONFIG')
do xyzzyaaae15=1,xyzzyaabd15
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabi15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading SEL&
&E.')
write(xyzzyaaab15)xyzzyaabi15
enddo
deallocate(xyzzyaabi15)
case('END VMC_SAVED_STATE')
exit
case default
call errstop('COPY_CONFIG_FILE','Unrecognized information '//trim(labe&
&l)//' in VMC_SAVED_STATE section of config file.')
end select
enddo
case('FINAL_VMC_RESULT')
do
read(xyzzyaaaa15,iostat=xyzzyaaac15)label
if(xyzzyaaac15/=0)call errstop_master('COPY_CONFIG_FILE',"Problem read&
&ing label in FINAL_VMC_RESULT section in config file.")
write(xyzzyaaab15)label
select case(trim(label))
case('ENERGY')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaacf15
if(xyzzyaaac15/=0)call errstop_master('COPY_CONFIG_FILE','Problem read&
&ing ENERGY')
write(xyzzyaaab15)xyzzyaacf15
case('ERRORBAR')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaacg15
if(xyzzyaaac15/=0)call errstop_master('COPY_CONFIG_FILE','Problem read&
&ing ERRORBAR')
write(xyzzyaaab15)xyzzyaacg15
case('ERRORBARU')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaach15
if(xyzzyaaac15/=0)call errstop_master('COPY_CONFIG_FILE','Problem read&
&ing ERRORBARU')
write(xyzzyaaab15)xyzzyaach15
case('VARIANCE')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaaci15
if(xyzzyaaac15/=0)call errstop_master('COPY_CONFIG_FILE','Problem read&
&ing VARIANCE')
write(xyzzyaaab15)xyzzyaaci15
case('CORRTIME')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaack15
if(xyzzyaaac15/=0)call errstop_master('COPY_CONFIG_FILE','Problem read&
&ing CORRTIME')
write(xyzzyaaab15)xyzzyaack15
case('TOTMOVE')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaacj15
if(xyzzyaaac15/=0)call errstop_master('COPY_CONFIG_FILE','Problem read&
&ing TOTMOVE')
write(xyzzyaaab15)xyzzyaacj15
case('ENERGY2')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaacl15
if(xyzzyaaac15/=0)call errstop_master('COPY_CONFIG_FILE','Problem read&
&ing ENERGY2')
write(xyzzyaaab15)xyzzyaacl15
case('ERRORBAR2')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaacm15
if(xyzzyaaac15/=0)call errstop_master('COPY_CONFIG_FILE','Problem read&
&ing ERRORBAR2')
write(xyzzyaaab15)xyzzyaacm15
case('ERRORBAR2U')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaacn15
if(xyzzyaaac15/=0)call errstop_master('COPY_CONFIG_FILE','Problem read&
&ing ERRORBAR2U')
write(xyzzyaaab15)xyzzyaacn15
case('END FINAL_VMC_RESULT')
exit
case default
call errstop_master('COPY_CONFIG_FILE','Unrecognized information '//tr&
&im(label)//' in FINAL_VMC_RESULT section of config file.')
end select
enddo
case('DMC_SAVED_STATE')
do
read(xyzzyaaaa15,iostat=xyzzyaaac15)label
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE',"Problem reading lab&
&el in DMC_SAVED_STATE section in config file.")
write(xyzzyaaab15)label
select case(trim(label))
case('EBEST')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabu15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading EBE&
&ST.')
write(xyzzyaaab15)xyzzyaabu15
case('EBEST_INIT')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabv15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading EBE&
&ST_INIT.')
write(xyzzyaaab15)xyzzyaabv15
case('EREF')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabw15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading ERE&
&F.')
write(xyzzyaaab15)xyzzyaabw15
case('DTEFF_EBEST_INIT')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabx15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading DTE&
&FF_EBEST_INIT.')
write(xyzzyaaab15)xyzzyaabx15
case('DTEFF_BEST')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaaby15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading DTE&
&FF_BEST.')
write(xyzzyaaab15)xyzzyaaby15
case('DMCEQUIL_STEPS')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaaba15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading DMC&
&EQUIL_STEPS.')
write(xyzzyaaab15)xyzzyaaba15
case('DMCSTATS_STEPS')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabb15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading DMC&
&STATS_STEPS.')
write(xyzzyaaab15)xyzzyaabb15
case('TPDMC')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabc15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading TPD&
&MC.')
write(xyzzyaaab15)xyzzyaabc15
case('LWDMC')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabf15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading LWD&
&MC.')
write(xyzzyaaab15)xyzzyaabf15
case('GROWTH_ESTIMATOR')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabg15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading GRO&
&WTH_ESTIMATOR.')
write(xyzzyaaab15)xyzzyaabg15
case('NUMER_EXPECT')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaace15(1:13)
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading NUM&
&ER_EXPECT.')
write(xyzzyaaab15)xyzzyaace15(1:13)
case('DENOMINATOR_WT')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaaca15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading DEN&
&OMINATOR_WT.')
write(xyzzyaaab15)xyzzyaaca15
case('LOG_PI_WT')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaacc15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading LOG&
&_PI_WT.')
write(xyzzyaaab15)xyzzyaacc15
case('LOG_PI_WT2')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaacd15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading LOG&
&_PI_WT2.')
write(xyzzyaaab15)xyzzyaacd15
case('LOG_PI_WT_ARRAY')
if(xyzzyaabc15<0)call errstop('COPY_CONFIG_FILE','LOG_PI_WT_ARRAY foun&
&d before TPDMC in config file. This is a bug.')
if(xyzzyaabc15==0)call errstop('COPY_CONFIG_FILE','LOG_PI_WT_ARRAY fou&
&nd in config file even though the stored TPDMC is zero. This is a bug&
&.')
allocate(xyzzyaacu15(0:xyzzyaabc15-1),stat=xyzzyaaad15)
call check_alloc(xyzzyaaad15,'COPY_CONFIG_FILE','LOG_PI_WT_ARRAY')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaacu15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading LOG&
&_PI_WT_ARRAY.')
write(xyzzyaaab15)xyzzyaacu15
deallocate(xyzzyaacu15)
case('LOG_PI_WT_ARRAY2')
if(xyzzyaabc15<0)call errstop('COPY_CONFIG_FILE','LOG_PI_WT_ARRAY2 fou&
&nd before TPDMC in config file. This is a bug.')
allocate(xyzzyaacv15(0:xyzzyaabc15),stat=xyzzyaaad15)
call check_alloc(xyzzyaaad15,'COPY_CONFIG_FILE','LOG_PI_WT_ARRAY2')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaacv15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading LOG&
&_PI_WT_ARRAY2.')
write(xyzzyaaab15)xyzzyaacv15
deallocate(xyzzyaacv15)
case('NUMERATOR_WT2')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabz15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading NUM&
&ERATOR_WT2.')
write(xyzzyaaab15)xyzzyaabz15
case('DENOMINATOR_WT2')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaacb15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading DEN&
&OMINATOR_WT2.')
write(xyzzyaaab15)xyzzyaacb15
case('END DMC_SAVED_STATE')
exit
case default
call errstop('COPY_CONFIG_FILE','Unrecognized information '//trim(labe&
&l)//' in DMC_SAVED_STATE section of config file.')
end select
enddo
case('REBLOCK_DATA')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaacy15,xyzzyaacz15,xyzzyaada1&
&5,xyzzyaadf15,xyzzyaadb15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading # o&
&f reblock blocks.')
write(xyzzyaaab15)xyzzyaacy15,xyzzyaacz15,xyzzyaada15,xyzzyaadf15,xyzz&
&yaadb15
allocate(xyzzyaadd15(xyzzyaacy15),xyzzyaade15(xyzzyaacy15),xyzzyaadh15&
&(xyzzyaacy15),xyzzyaadi15(xyzzyaacy15),xyzzyaadj15(xyzzyaada15,xyzzya&
&acy15),xyzzyaadk15(xyzzyaada15,xyzzyaacy15),xyzzyaadl15(xyzzyaada15),&
&stat=xyzzyaaad15)
call check_alloc(xyzzyaaad15,'COPY_CONFIG_FILE','reblock_block_length_&
&config')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaadd15,xyzzyaade15,xyzzyaadh1&
&5,xyzzyaadi15,xyzzyaadj15,xyzzyaadk15,xyzzyaadl15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading reb&
&lock blocks.')
write(xyzzyaaab15)xyzzyaadd15,xyzzyaade15,xyzzyaadh15,xyzzyaadi15,xyzz&
&yaadj15,xyzzyaadk15,xyzzyaadl15
deallocate(xyzzyaadd15,xyzzyaade15,xyzzyaadh15,xyzzyaadi15,xyzzyaadj15&
&,xyzzyaadk15,xyzzyaadl15)
if(int2log(xyzzyaadb15))then
allocate( xyzzyaadm15(xyzzyaada15),xyzzyaadn15(xyzzyaada15),xyzzyaado1&
&5(xyzzyaada15),xyzzyaadp15(xyzzyaada15),stat=xyzzyaaad15)
call check_alloc(xyzzyaaad15,'COPY_CONFIG_FILE','reblock_sum_ow0_in')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaadg15,xyzzyaadm15,xyzzyaadn1&
&5,xyzzyaado15,xyzzyaadp15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading pop&
&stats.')
write(xyzzyaaab15)xyzzyaadg15,xyzzyaadm15,xyzzyaadn15,xyzzyaado15,xyzz&
&yaadp15
deallocate( xyzzyaadm15,xyzzyaadn15,xyzzyaado15,xyzzyaadp15)
endif
read(xyzzyaaaa15,iostat=xyzzyaaac15)label
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading lab&
&el END REBLOCK_DATA in config file.')
if(trim(label)/='END REBLOCK_DATA')call errstop('COPY_CONFIG_FILE','Wr&
&ong label instead of END REBLOCK_DATA in config file.')
write(xyzzyaaab15)label
case('RANDOM')
read(xyzzyaaaa15,iostat=xyzzyaaac15)label
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading lab&
&el NNODES in RANDOM section of config file.')
if(trim(label)/='NNODES')call errstop('COPY_CONFIG_FILE','Wrong label &
&instead of NNODES in RANDOM section of config file.')
write(xyzzyaaab15)label
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabd15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading NNO&
&DES in RANDOM section of config file.')
write(xyzzyaaab15)xyzzyaabd15
do xyzzyaaae15=1,xyzzyaabd15
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabe15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading RAN&
&DOM_STATE for node '//trim(i2s(xyzzyaaae15))//' in RANDOM section of &
&config file.')
write(xyzzyaaab15)xyzzyaabe15
enddo
read(xyzzyaaaa15,iostat=xyzzyaaac15)label
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading lab&
&el END RANDOM/BUFFER in config file.')
if(trim(label)=='BUFFER')then
write(xyzzyaaab15)label
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaaah15,xyzzyaaai15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading RAN&
&X_MAX and RANX_GAUSS_MAX in RANDOM section of config file.')
if(xyzzyaaah15/=ranx_max.and.xyzzyaaai15/=ranx_gauss_max)call errstop(&
&'COPY_CONFIG_FILE','Config file was generated with random-number buff&
&ers of a different size to those currently defined internally in CASI&
&NO.')
write(xyzzyaaab15)xyzzyaaah15,xyzzyaaai15
do xyzzyaaah15=1,xyzzyaabd15
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaadr15,xyzzyaads15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading RAN&
&X_INDX in RANDOM section of config file.')
write(xyzzyaaab15)xyzzyaadr15,xyzzyaads15
enddo
allocate(xyzzyaadt15(ranx_max))
do xyzzyaaah15=1,xyzzyaabd15
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaadt15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading RAN&
&X_BUFFER in RANDOM section of config file.')
write(xyzzyaaab15)xyzzyaadt15
enddo
deallocate(xyzzyaadt15)
allocate(xyzzyaadu15(ranx_gauss_max))
do xyzzyaaah15=1,xyzzyaabd15
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaadu15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading RAN&
&X_GAUSS_BUFFER in RANDOM section of config file.')
write(xyzzyaaab15)xyzzyaadu15
enddo
deallocate(xyzzyaadu15)
read(xyzzyaaaa15,iostat=xyzzyaaac15)label
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading lab&
&el END RANDOM in config file <2>.')
endif
if(trim(label)/='END RANDOM')call errstop('COPY_CONFIG_FILE','Wrong la&
&bel instead of END RANDOM in config file.')
write(xyzzyaaab15)label
case('GEOMETRY')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaadc15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading NIT&
&OT.')
write(xyzzyaaab15)xyzzyaadc15
if(xyzzyaadc15>0)then
allocate(xyzzyaadq15(3,xyzzyaadc15),stat=xyzzyaaad15)
call check_alloc(xyzzyaaad15,'COPY_CONFIG_FILE','RION_IN')
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaadq15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading RIO&
&N.')
write(xyzzyaaab15)xyzzyaadq15
deallocate(xyzzyaadq15)
endif
read(xyzzyaaaa15,iostat=xyzzyaaac15)label
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE',"Problem reading 'EN&
&D GEOMETRY' label in config file.")
if(trim(label)/='END GEOMETRY')call errstop('COPY_CONFIG_FILE',"Wrong &
&label read in when expecting 'END GEOMETRY'.")
write(xyzzyaaab15)label
case('CONFIGS')
goto 10
case default
call errstop('COPY_CONFIG_FILE','Unrecognized section '//trim(label)//&
&' in config file.')
end select
enddo
10 continue
if(xyzzyaaaw15<0)call errstop('COPY_CONFIG_FILE','Number of configs un&
&defined when about to read configs.')
if(xyzzyaaaw15>0)then
if(xyzzyaaaq15==0)call errstop('COPY_CONFIG_FILE','Config contents und&
&efined when about to read configs.')
do xyzzyaaar15=1,xyzzyaaaq15
select case(trim(xyzzyaacx15(xyzzyaaar15)))
case('RELE')
allocate(xyzzyaacq15(3,xyzzyaaap15),stat=xyzzyaaad15)
case('SELE')
allocate(xyzzyaabj15(xyzzyaaap15),stat=xyzzyaaad15)
case('LOGDET')
allocate(xyzzyaacw15(xyzzyaaam15,xyzzyaaan15),stat=xyzzyaaad15)
case('FIDET')
allocate(xyzzyaacr15(3,xyzzyaaap15,xyzzyaaav15),stat=xyzzyaaad15)
case('FI_PROD_DET')
allocate(xyzzyaacs15(3,xyzzyaaan15,xyzzyaaao15,xyzzyaaav15,xyzzyaaam15&
&),stat=xyzzyaaad15)
case('PROD_LAPDET')
allocate(xyzzyaact15(xyzzyaaan15,xyzzyaaao15,xyzzyaaav15,xyzzyaaam15),&
&stat=xyzzyaaad15)
case default
xyzzyaaad15=0
end select
call check_alloc(xyzzyaaad15,'COPY_CONFIG_FILE',trim(xyzzyaacx15(xyzzy&
&aaar15)))
enddo
do xyzzyaaax15=1,xyzzyaaaw15
if(xyzzyaaal15(xyzzyaaaa1))then
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaacq15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading REL&
&E for config '//trim(i2s(xyzzyaaax15))//'.')
write(xyzzyaaab15)xyzzyaacq15
endif
if(xyzzyaaal15(xyzzyaaab1))then
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabj15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading SEL&
&E for config '//trim(i2s(xyzzyaaax15))//'.')
write(xyzzyaaab15)xyzzyaabj15
endif
if(xyzzyaaal15(xyzzyaaac1))then
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabm15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading ETO&
&T for config '//trim(i2s(xyzzyaaax15))//'.')
write(xyzzyaaab15)xyzzyaabm15
endif
if(xyzzyaaal15(xyzzyaaad1))then
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaacw15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading LOG&
&DET for config '//trim(i2s(xyzzyaaax15))//'.')
write(xyzzyaaab15)xyzzyaacw15
endif
if(xyzzyaaal15(xyzzyaaae1))then
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaacr15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading FID&
&ET for config '//trim(i2s(xyzzyaaax15))//'.')
write(xyzzyaaab15)xyzzyaacr15
endif
if(xyzzyaaal15(xyzzyaaaf1))then
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaacs15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading FI_&
&PROD_DET for config '//trim(i2s(xyzzyaaax15))//'.')
write(xyzzyaaab15)xyzzyaacs15
endif
if(xyzzyaaal15(xyzzyaaag1))then
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabl15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading LAP&
&DET for config '//trim(i2s(xyzzyaaax15))//'.')
write(xyzzyaaab15)xyzzyaabl15
endif
if(xyzzyaaal15(xyzzyaaah1))then
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaact15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading PRO&
&D_LAPDET for config '//trim(i2s(xyzzyaaax15))//'.')
write(xyzzyaaab15)xyzzyaact15
endif
if(xyzzyaaal15(xyzzyaaai1))then
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabn15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading LOC&
&AL_POTENTIAL for config '//trim(i2s(xyzzyaaax15))//'.')
write(xyzzyaaab15)xyzzyaabn15
endif
if(xyzzyaaal15(xyzzyaaaj1))then
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabo15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading NLT&
&OT for config '//trim(i2s(xyzzyaaax15))//'.')
write(xyzzyaaab15)xyzzyaabo15
endif
if(xyzzyaaal15(xyzzyaaak1))then
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabp15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading STO&
&T for config '//trim(i2s(xyzzyaaax15))//'.')
write(xyzzyaaab15)xyzzyaabp15
endif
if(xyzzyaaal15(xyzzyaaal1))then
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabq15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading WDM&
&C for config '//trim(i2s(xyzzyaaax15))//'.')
write(xyzzyaaab15)xyzzyaabq15
endif
if(xyzzyaaal15(xyzzyaaam1))then
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabr15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading VAL&
&JAS for config '//trim(i2s(xyzzyaaax15))//'.')
write(xyzzyaaab15)xyzzyaabr15
endif
if(xyzzyaaal15(xyzzyaaan1))then
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabs15
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading LOG&
&P for config '//trim(i2s(xyzzyaaax15))//'.')
write(xyzzyaaab15)xyzzyaabs15
endif
if(xyzzyaaal15(xyzzyaaao1))then
read(xyzzyaaaa15,iostat=xyzzyaaac15)xyzzyaabt15(:)
if(xyzzyaaac15/=0)call errstop('COPY_CONFIG_FILE','Problem reading TWI&
&ST for config '//trim(i2s(xyzzyaaax15))//'.')
write(xyzzyaaab15)xyzzyaabt15(:)
endif
enddo
do xyzzyaaar15=1,xyzzyaaaq15
select case(trim(xyzzyaacx15(xyzzyaaar15)))
case('RELE')
deallocate(xyzzyaacq15)
case('SELE')
deallocate(xyzzyaabj15)
case('LOGDET')
deallocate(xyzzyaacw15)
case('FIDET')
deallocate(xyzzyaacr15)
case('FI_PROD_DET')
deallocate(xyzzyaacs15)
case('PROD_LAPDET')
deallocate(xyzzyaact15)
end select
enddo
endif
if(xyzzyaaaj15)then
close(xyzzyaaaa15,status='delete')
else
close(xyzzyaaaa15)
endif
close(xyzzyaaab15)
open_unit(xyzzyaaaa15)=.false.
open_unit(xyzzyaaab15)=.false.
call timer('COPY_CONFIG_FILE',.false.)
endif
if(allocated(xyzzyaacx15))deallocate(xyzzyaacx15)
call qmc_barrier
end subroutine copy_config_file
end module slaarnaaf

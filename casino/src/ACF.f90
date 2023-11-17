module slaarnacf
use dsp
use slaarnaam
use slaarnach
use store
use slaarnacs
use slaarnaab,only : forces_scratch_request,setup_forces_accum,eval_lo&
&cal_forces,nfcomps
use slaarnaan,only : points
use file_utils,   only : open_units
use format_utils, only : wout,i2s
use slaarnabg,     only : dimensionality
use slaarnabt,    only : dnrm2,ddot
use parallel,     only : am_master
use slaarnaca,        only : have_ppots
use run_control,  only : errstop,check_alloc
use slaarnabs,    only : v_non_local_tmove
implicit none
private
public rmc_main
integer xyzzyaaaa1
integer xyzzyaaab1
type rmc_config
type(config_geom),pointer :: pt_geom=>null()
type(config_wfn),pointer :: pt_wfn=>null()
real(dp),pointer :: pt_drift(:,:)=>null(),pt_obs(:)=>null()
real(dp) eloc
real(dp) mag2_drift
integer(i64) id
end type
type(rmc_config),allocatable :: xyzzyaaac1(:)
real(dp),allocatable,target :: xyzzyaaad1(:,:,:),xyzzyaaae1(:,:)
integer xyzzyaaaf1,xyzzyaaag1,xyzzyaaah1,xyzzyaaai1,xyzzyaaaj1,xyzzyaa&
&ak1,xyzzyaaal1,xyzzyaaam1,xyzzyaaan1,xyzzyaaao1,xyzzyaaap1,xyzzyaaaq1&
&,xyzzyaaar1,xyzzyaaas1
integer(i64) xyzzyaaat1
logical xyzzyaaau1
real(dp),parameter :: xyzzyaaav1=0.5d0
logical,parameter :: xyzzyaaaw1=.true.
contains
subroutine rmc_main(iaccum,bounce,dtrmc,nmove_rmc,nblock_rmc,rmc_rep_l&
&ength,rmc_move_length,corper_rmc,rmc_meas_pos,aborted)
use slaarnacc, only : ranx
implicit none
integer,intent(in) :: nmove_rmc,nblock_rmc,rmc_rep_length,rmc_move_len&
&gth,corper_rmc
real(dp),intent(in) :: dtrmc
logical,intent(in) :: iaccum,bounce,rmc_meas_pos
logical,intent(out) :: aborted
integer i,xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2(netot),xyzzyaaae&
&2,xyzzyaaaf2
integer(i64) xyzzyaaag2
integer(i64),parameter :: xyzzyaaah2=0,xyzzyaaai2=-1,xyzzyaaaj2=-2
real(dp) xyzzyaaak2,xyzzyaaal2,xyzzyaaam2(3,netot),xyzzyaaan2(3,netot)&
&,xyzzyaaao2,xyzzyaaap2,xyzzyaaaq2,xyzzyaaar2,xyzzyaaas2
complex(dp) xyzzyaaat2,xyzzyaaau2
logical xyzzyaaav2(netot),xyzzyaaaw2
logical,parameter :: xyzzyaaax2=.true.
character(256) tmpr
if(am_master)then
call wout()
call wout('BEGIN RMC CALCULATION')
call wout('=====================')
call wout()
endif
aborted=.false.
if(xyzzyaaax2)then
call open_units(xyzzyaaaa1,xyzzyaaae2)
if(xyzzyaaae2/=0)call errstop('RMC_MAIN','Ran out of I/O units <1>.')
open(unit=xyzzyaaaa1,file='rmc.rep',form='unformatted',status='unknown&
&',iostat=xyzzyaaae2)
if(xyzzyaaae2/=0)call errstop('RMC_MAIN','Problem opening rmc.rep.')
call open_units(xyzzyaaab1,xyzzyaaae2)
if(xyzzyaaae2/=0)call errstop('RMC_MAIN','Ran out of I/O units <2>.')
open(unit=xyzzyaaab1,file='rmc.config',form='unformatted',status='unkn&
&own',iostat=xyzzyaaae2)
if(xyzzyaaae2/=0)call errstop('RMC_MAIN','Problem opening rmc.config.'&
&)
endif
xyzzyaaak1=0
xyzzyaaam1=0
call scratch_protect(xyzzyaaak1)
call scratch_protect(xyzzyaaam1)
call scratch_request(ratiocfg_from=xyzzyaaak1,ratiocfg_to=xyzzyaaam1)
call energy_scratch_request(xyzzyaaak1)
call energy_scratch_request(xyzzyaaam1)
if(forces)then
call forces_scratch_request(xyzzyaaak1)
call forces_scratch_request(xyzzyaaam1)
endif
call scratch_request(drift=xyzzyaaak1)
call scratch_request(drift=xyzzyaaam1)
call energy_scratch_request(xyzzyaaak1)
call setup_scratch
call which_scratch(xyzzyaaak1)
call which_scratch(xyzzyaaam1)
call setup_wfn_utils
call setup_energy_utils
if(forces)call setup_forces_accum
xyzzyaaav2(:)=.false.
xyzzyaaan2=0.d0
call points(xyzzyaaam2,xyzzyaaad2,xyzzyaaan2,xyzzyaaav2,printout=.true&
&.)
call xyzzyaabd1(rmc_rep_length,rmc_move_length,xyzzyaaam2,dtrmc,rmc_me&
&as_pos)
if(xyzzyaaax2)then
do i=1,rmc_rep_length
call xyzzyaaax1(xyzzyaaac1(xyzzyaabf1(i)))
enddo
endif
xyzzyaaaf2=0
xyzzyaaac2=0
xyzzyaaas2=0.d0
do xyzzyaaaa2=1,nblock_rmc
do xyzzyaaab2=1,nmove_rmc
if(.not.bounce)then
if(ranx()<0.5d0)then
if(xyzzyaaax2)write(xyzzyaaaa1)xyzzyaaah2,xyzzyaaaj2
xyzzyaaaf1=xyzzyaabf1(xyzzyaaah1+1)
xyzzyaaag1=-xyzzyaaag1
endif
endif
xyzzyaaap2=1.d0
xyzzyaaaq2=1.d0
xyzzyaaak2=0.d0
do i=1,rmc_move_length
call xyzzyaabb1(xyzzyaaac1(xyzzyaabh1(i-1)),xyzzyaaac1(xyzzyaabh1(i)),&
&dtrmc,xyzzyaaar2)
xyzzyaaap2=xyzzyaaap2*xyzzyaaar2
xyzzyaaaq2=xyzzyaaaq2*xyzzyaaay1(xyzzyaaac1(xyzzyaabg1(i+1)),xyzzyaaac&
&1(xyzzyaabg1(i)),dtrmc)
xyzzyaaak2=xyzzyaaak2+xyzzyaaaz1(xyzzyaaac1(xyzzyaabh1(i-1)),xyzzyaaac&
&1(xyzzyaabh1(i)),dtrmc)-xyzzyaaaz1(xyzzyaaac1(xyzzyaabg1(i)),xyzzyaaa&
&c1(xyzzyaabg1(i+1)),dtrmc)
if(xyzzyaaax2)call xyzzyaaax1(xyzzyaaac1(xyzzyaabh1(i)))
enddo
xyzzyaaat2=xyzzyaaba1(xyzzyaaac1(xyzzyaabh1(0)),xyzzyaaac1(xyzzyaabh1(&
&rmc_move_length)))
xyzzyaaau2=xyzzyaaba1(xyzzyaaac1(xyzzyaabg1(1)),xyzzyaaac1(xyzzyaabg1(&
&rmc_move_length+1)))
xyzzyaaal2=abs(xyzzyaaat2*xyzzyaaau2)*xyzzyaaaq2/xyzzyaaap2*exp(xyzzya&
&aak2)
xyzzyaaaw2=(real(xyzzyaaat2,dp)<0.d0)
if(.not.xyzzyaaaw2.and.xyzzyaaal2>ranx())then
if(xyzzyaaax2)then
do i=1,rmc_move_length
write(xyzzyaaaa1)xyzzyaaac1(xyzzyaabg1(i))%id,xyzzyaaac1(xyzzyaabh1(i)&
&)%id
enddo
endif
xyzzyaaac2=xyzzyaaac2+1
xyzzyaaaf1=xyzzyaabg1(1)
else
if(xyzzyaaaw2)xyzzyaaaf2=xyzzyaaaf2+1
if(bounce)then
if(xyzzyaaax2)write(xyzzyaaaa1)xyzzyaaah2,xyzzyaaaj2
xyzzyaaaf1=xyzzyaabf1(xyzzyaaah1+1)
xyzzyaaag1=-xyzzyaaag1
endif
if(xyzzyaaax2)write(xyzzyaaaa1)xyzzyaaah2,xyzzyaaai2
endif
xyzzyaaas2=xyzzyaaas2+xyzzyaaac1(xyzzyaabf1(1))%eloc
enddo
enddo
call xyzzyaabe1
call finish_energy_utils
call finish_wfn_utils
call finish_scratch
if(am_master)then
xyzzyaaag2=xyzzyaabi1()-1_i64
write(tmpr,*)xyzzyaaag2
call wout('Last ID used            : '//trim(adjustl(tmpr)))
call wout('RMC energy              : ',xyzzyaaas2/dble(nblock_rmc*nmov&
&e_rmc))
xyzzyaaao2=dble(xyzzyaaac2)/dble(nblock_rmc*nmove_rmc)
call wout('Acceptance              : ',xyzzyaaao2)
if(xyzzyaaao2<1.d0)then
call wout('Average bounce length   : ',1.d0/(1.d0-xyzzyaaao2))
else
call wout('Average bounce length   : Infinity')
endif
call wout('Number of node crossings: '//trim(i2s(xyzzyaaaf2)))
endif
if(xyzzyaaax2)then
close(xyzzyaaaa1)
open_unit(xyzzyaaaa1)=.false.
close(xyzzyaaab1)
open_unit(xyzzyaaab1)=.false.
endif
end subroutine rmc_main
subroutine xyzzyaaax1(cf)
type(rmc_config),intent(in) :: cf
integer xyzzyaaaa3(netot)
real(dp) xyzzyaaab3(3,netot)
call load_from_pt(xyzzyaaak1,cf%pt_geom,cf%pt_wfn,xyzzyaaab3,xyzzyaaaa&
&3)
write(xyzzyaaab1)cf%id,cf%pt_obs
end subroutine xyzzyaaax1
real(dp) function xyzzyaaay1(c_from,c_to,dtrmc)
implicit none
real(dp),intent(in) :: dtrmc
type(rmc_config),intent(in) :: c_from,c_to
integer xyzzyaaaa4(3,netot),xyzzyaaab4(3,netot)
real(dp) xyzzyaaac4(3,netot),xyzzyaaad4(3,netot),xyzzyaaae4
call load_from_pt(xyzzyaaak1,c_from%pt_geom,c_from%pt_wfn,xyzzyaaac4,x&
&yzzyaaaa4)
call load_from_pt(xyzzyaaam1,c_to%pt_geom,c_to%pt_wfn,xyzzyaaad4,xyzzy&
&aaab4)
xyzzyaaae4=dnrm2(3*netot,xyzzyaaad4-xyzzyaaac4-2.d0*xyzzyaaav1*dtrmc*c&
&_from%pt_drift,1)
xyzzyaaay1=exp(-xyzzyaaae4**2/(4.d0*xyzzyaaav1*dtrmc))
end function xyzzyaaay1
real(dp) function xyzzyaaaz1(c_from,c_to,dtrmc)
implicit none
real(dp),intent(in) :: dtrmc
type(rmc_config),intent(in) :: c_from,c_to
integer xyzzyaaaa5(3,netot),xyzzyaaab5(3,netot)
real(dp) xyzzyaaac5(3,netot),xyzzyaaad5(3,netot),xyzzyaaae5,xyzzyaaaf5
call load_from_pt(xyzzyaaak1,c_from%pt_geom,c_from%pt_wfn,xyzzyaaac5,x&
&yzzyaaaa5)
call load_from_pt(xyzzyaaam1,c_to%pt_geom,c_to%pt_wfn,xyzzyaaad5,xyzzy&
&aaab5)
xyzzyaaae5=dnrm2(3*netot,xyzzyaaad5-xyzzyaaac5,1)
xyzzyaaaf5=ddot(3*netot,xyzzyaaad5-xyzzyaaac5,1,c_to%pt_drift-c_from%p&
&t_drift,1)
xyzzyaaaz1=-0.5d0*dtrmc*(c_from%eloc+c_to%eloc) -0.5d0*dtrmc*xyzzyaaav&
&1*(c_from%mag2_drift+c_to%mag2_drift) -xyzzyaaae5**2/(4.d0*xyzzyaaav1&
&*dtrmc)-0.5d0*xyzzyaaaf5
end function xyzzyaaaz1
complex(dp) function xyzzyaaba1(c_from,c_to)
implicit none
type(rmc_config),intent(in) :: c_from,c_to
integer xyzzyaaaa6(netot)
real(dp) xyzzyaaab6(3,netot)
logical isnan,isinf
call load_from_pt(xyzzyaaak1,c_from%pt_geom,c_from%pt_wfn,xyzzyaaab6,x&
&yzzyaaaa6)
call load_from_pt(xyzzyaaam1,c_to%pt_geom,c_to%pt_wfn,xyzzyaaab6,xyzzy&
&aaaa6)
call wfn_ratio(xyzzyaaak1,xyzzyaaam1,0,ratio=xyzzyaaba1,isnan=isnan,is&
&inf=isinf)
end function xyzzyaaba1
subroutine xyzzyaabb1(c_from,c_to,dtrmc,prob)
use slaarnabg,      only : dimensionality
use slaarnacc,only : ranx_gaussian
implicit none
real(dp),intent(in) :: dtrmc
real(dp),intent(out) :: prob
type(rmc_config),intent(in):: c_from
type(rmc_config),intent(inout):: c_to
integer xyzzyaaaa7(netot),i,xyzzyaaab7
real(dp) xyzzyaaac7(3,netot),xyzzyaaad7(3,netot),xyzzyaaae7(3,netot)
call load_from_pt(xyzzyaaak1,c_from%pt_geom,c_from%pt_wfn,xyzzyaaac7,x&
&yzzyaaaa7)
do i=1,netot
do xyzzyaaab7=1,dimensionality
xyzzyaaad7(xyzzyaaab7,i)=ranx_gaussian(sqrt(2.d0*dtrmc*xyzzyaaav1))
enddo
enddo
xyzzyaaae7=xyzzyaaac7+2.d0*xyzzyaaav1*dtrmc*c_from%pt_drift+xyzzyaaad7
prob=exp(-dnrm2(3*netot,xyzzyaaad7,1)**2/(4.d0*xyzzyaaav1*dtrmc))
call define_config(xyzzyaaak1,xyzzyaaae7,xyzzyaaaa7)
call save_to_pt(xyzzyaaak1,c_to%pt_geom,c_to%pt_wfn)
call xyzzyaabc1(c_to,.true.)
c_to%id=xyzzyaabi1()
end subroutine xyzzyaabb1
subroutine xyzzyaabc1(config,do_meas)
implicit none
logical,intent(in) :: do_meas
type(rmc_config),intent(inout) :: config
integer xyzzyaaaa8(netot),xyzzyaaab8,xyzzyaaac8
real(dp) xyzzyaaad8(3,netot),ecomps(n_ecomp),xyzzyaaae8,xyzzyaaaf8
complex(dp) xyzzyaaag8(3)
logical isnan,isinf
call load_from_pt(xyzzyaaak1,config%pt_geom,config%pt_wfn,xyzzyaaad8,x&
&yzzyaaaa8)
do xyzzyaaab8=1,netot
call wfn_loggrad(xyzzyaaab8,xyzzyaaak1,0,xyzzyaaag8,isnan=isnan,isinf=&
&isinf)
config%pt_drift(:,xyzzyaaab8)=real(xyzzyaaag8,dp)
enddo
config%mag2_drift=dnrm2(3*netot,config%pt_drift,1)**2
call eval_local_energy(xyzzyaaak1,etot=config%eloc,ecomps=ecomps,isnan&
&=isnan,isinf=isinf)
if(do_meas)then
if(xyzzyaaaw1)config%pt_obs(xyzzyaaao1)=config%eloc
if(xyzzyaaau1)then
do xyzzyaaab8=1,netot
do xyzzyaaac8=1,dimensionality
config%pt_obs(xyzzyaaaq1+dimensionality*(xyzzyaaab8-1)+xyzzyaaac8-1)=x&
&yzzyaaad8(xyzzyaaac8,xyzzyaaab8)
enddo
enddo
endif
if(forces)then
xyzzyaaae8=0.d0
xyzzyaaaf8=0.d0
if(have_ppots)then
if(.not.use_tmove)then
xyzzyaaae8=ecomps(i_potinl)
else
call v_non_local_tmove(xyzzyaaae8,xyzzyaaaf8)
endif
endif
call eval_local_forces(xyzzyaaak1,config%eloc,ecomps(i_kei),xyzzyaaae8&
&,xyzzyaaaf8,config%pt_obs(xyzzyaaas1:(xyzzyaaas1+xyzzyaaar1-1)))
endif
endif
end subroutine xyzzyaabc1
subroutine xyzzyaabd1(body_length,prop_move_length,initial_rele,dtrmc,&
&meas_pos_set)
use slaarnabt, only : dcopy
implicit none
integer,intent(in) :: body_length,prop_move_length
real(dp),intent(in) :: dtrmc,initial_rele(3,netot)
logical,intent(in)  :: meas_pos_set
integer xyzzyaaaa9(netot),i,xyzzyaaab9
real(dp) xyzzyaaac9
xyzzyaaau1=meas_pos_set
xyzzyaaaa9=0
xyzzyaaag1=1
xyzzyaaah1=body_length
xyzzyaaai1=prop_move_length
xyzzyaaaj1=xyzzyaaah1+xyzzyaaai1
xyzzyaaaf1=xyzzyaaaj1-1
xyzzyaaat1=0
xyzzyaaal1=0
if(xyzzyaaaw1)then
xyzzyaaan1=1
xyzzyaaao1=xyzzyaaal1+1
xyzzyaaal1=xyzzyaaal1+xyzzyaaan1
endif
if(xyzzyaaau1)then
xyzzyaaap1=dimensionality*netot
xyzzyaaaq1=xyzzyaaal1+1
xyzzyaaal1=xyzzyaaal1+xyzzyaaap1
endif
if(forces)then
xyzzyaaar1=nfcomps
xyzzyaaas1=xyzzyaaal1+1
xyzzyaaal1=xyzzyaaal1+xyzzyaaar1
endif
allocate(xyzzyaaac1(xyzzyaaaj1),xyzzyaaad1(3,netot,xyzzyaaaj1),xyzzyaa&
&ae1(xyzzyaaal1,xyzzyaaaj1),stat=xyzzyaaab9)
call check_alloc(xyzzyaaab9,'INIT_REPTILE','')
do i=1,xyzzyaaaj1
call gen_config(xyzzyaaac1(i)%pt_geom,xyzzyaaac1(i)%pt_wfn)
xyzzyaaac1(i)%pt_drift=>xyzzyaaad1(:,:,i)
xyzzyaaac1(i)%pt_obs=>xyzzyaaae1(:,i)
enddo
call define_config(xyzzyaaak1,initial_rele,xyzzyaaaa9)
call dcopy(three_netot,initial_rele(1,1),1,rele_scr(1,1,xyzzyaaak1),1)
call save_to_pt(xyzzyaaak1,xyzzyaaac1(xyzzyaabf1(1))%pt_geom,xyzzyaaac&
&1(xyzzyaabf1(1))%pt_wfn)
xyzzyaaac1(xyzzyaabf1(1))%id=xyzzyaabi1()
call xyzzyaabc1(xyzzyaaac1(xyzzyaabf1(1)),.true.)
do i=2,xyzzyaaah1
call xyzzyaabb1(xyzzyaaac1(xyzzyaabf1(i-1)),xyzzyaaac1(xyzzyaabf1(i)),&
&dtrmc,xyzzyaaac9)
enddo
end subroutine xyzzyaabd1
subroutine xyzzyaabe1
implicit none
deallocate(xyzzyaaac1,xyzzyaaad1,xyzzyaaae1)
end subroutine xyzzyaabe1
integer function xyzzyaabf1(i)
implicit none
integer,intent(in) :: i
xyzzyaabf1=1+modulo(xyzzyaaaf1+xyzzyaaag1*i-1,xyzzyaaaj1)
end function xyzzyaabf1
integer function xyzzyaabg1(i)
implicit none
integer,intent(in) :: i
xyzzyaabg1=xyzzyaabf1(xyzzyaaah1-i+1)
end function xyzzyaabg1
integer function xyzzyaabh1(i)
implicit none
integer,intent(in) :: i
xyzzyaabh1=xyzzyaabf1(1-i)
end function xyzzyaabh1
integer(i64) function xyzzyaabi1()
implicit none
xyzzyaaat1=xyzzyaaat1+1
xyzzyaabi1=xyzzyaaat1
end function xyzzyaabi1
end module slaarnacf

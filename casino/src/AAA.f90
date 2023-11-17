module slaarnaaa
use slaarnaag
use dsp
use parallel
use slaarnacg
use slaarnach
use store
use slaarnacq
use slaarnacs
use slaarnaam,only : eval_local_energy
use format_utils,only : wout,r2s
use run_control, only : errstop,check_alloc
implicit none
private
public setup_alt_utils,finish_alt_utils,alt_accept,alt_reject,wfn_altr&
&atio,alt_halflogpdf,weight_config_vmc,fieller,complextosimple,simplet&
&ocomplex
public f_scr,logp_scr,f_valid,logp_valid
real(dp),allocatable :: f_scr(:),logp_scr(:)
logical,allocatable :: f_valid(:),logp_valid(:)
logical xyzzyaaaa1,xyzzyaaab1,xyzzyaaac1
contains
subroutine complextosimple(is0,is1)
implicit none
integer,intent(in),optional :: is0,is1
if(simplepdf==1)then
if(.not.isitcomplex)call errstop('COMPLEXTOSIMPLE','Called in simple.'&
&)
if(present(is0))call empty_scratch_wfn(is0)
if(present(is1))call empty_scratch_wfn(is1)
use_backflow=.false.
use_jastrow=.false.
use_orbmods=.false.
if(altsamp==2)then
detstart=1
detstop=ndet_smp
else
detstart=1
detstop=1
endif
isitcomplex=.false.
wfdet_orbmask=wfdet_orbmask_scr2
endif
end subroutine complextosimple
subroutine simpletocomplex(is0,is1)
implicit none
integer,intent(in),optional :: is0,is1
if(simplepdf==1)then
if(isitcomplex)call errstop('SIMPLETOCOMPLEX','Called in complex.')
detstart=1
detstop=ndet
use_backflow=xyzzyaaab1
use_jastrow=xyzzyaaaa1
use_orbmods=xyzzyaaac1
if(present(is0))call empty_scratch_wfn(is0)
if(present(is1))call empty_scratch_wfn(is1)
isitcomplex=.true.
wfdet_orbmask=wfdet_orbmask_scr1
endif
end subroutine simpletocomplex
subroutine setup_alt_utils
implicit none
integer xyzzyaaaa4
allocate(f_scr(nscratch),f_valid(nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_ALT_UTILS','f_scr')
f_scr=0.0
f_valid=.false.
allocate(logp_scr(nscratch),logp_valid(nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_ALT_UTILS','logp_scr')
logp_scr=0.d0
logp_valid=.false.
allocate(wfdet_orbmask_scr1(wfdet_norb,nspin),wfdet_orbmask_scr2(wfdet&
&_norb,nspin),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_ALT_UTILS','wfdet_orbmask_scr*')
wfdet_orbmask_scr1=wfdet_orbmask
wfdet_orbmask_scr2=wfdet_orbmask
if(use_altsamp)then
if(simplepdf==1)then
wfdet_orbmask_scr2=.false.
detstart=1
detstop=ndet_smp
call get_wfdet_orbmask(wfdet_norb,wfdet_orbmap,wfdet_orbmask_scr2)
detstart=1
detstop=ndet
elseif(simplepdf==2)then
call setup_sampling_utils
endif
endif
xyzzyaaab1=use_backflow
xyzzyaaaa1=use_jastrow
xyzzyaaac1=use_orbmods
end subroutine setup_alt_utils
subroutine finish_alt_utils
implicit none
deallocate(f_valid,logp_valid)
deallocate(f_scr,logp_scr)
deallocate(wfdet_orbmask_scr1,wfdet_orbmask_scr2)
if(use_altsamp.and.simplepdf==2)call finish_sampling_utils
end subroutine finish_alt_utils
subroutine alt_accept(is,js)
implicit none
integer,intent(in) :: is,js
f_valid(is)=f_valid(js)
logp_valid(is)=logp_valid(js)
f_scr(is)=f_scr(js)
logp_scr(is)=logp_scr(js)
f_valid(js)=.false.
logp_valid(js)=.false.
if(use_altsamp.and.simplepdf==2)call sampling_accept_move(is,js)
end subroutine alt_accept
subroutine alt_reject(js)
implicit none
integer,intent(in) :: js
f_valid(js)=.false.
logp_valid(js)=.false.
if(use_altsamp.and.simplepdf==2)call clear_scratch_sampling(js)
end subroutine alt_reject
subroutine wfn_altratio(is,js,ilevel,relprob,truncprob,isnan,isinf)
implicit none
integer,intent(in) :: is,js,ilevel
real(dp),intent(out),optional :: relprob,truncprob
logical,intent(out),optional :: isnan,isinf
integer xyzzyaaaa8
real(dp) xyzzyaaab8,xyzzyaaac8,xyzzyaaad8
if(altsamp/=1.or..not.use_jastrow.or.simplepdf==2)then
xyzzyaaaa8=ilevel
if(ilevel==2.and.(.not.use_jastrow.or.simplepdf==2))then
if(present(relprob))relprob=1.d0
if(present(truncprob))truncprob=1.d0
return
endif
else
if(ilevel==2)xyzzyaaaa8=1
if(ilevel==1)xyzzyaaaa8=2
endif
if(simplepdf==2)then
call sampling_ratio(is,js,relprob=xyzzyaaab8,isnan=isnan,isinf=isinf)
else
call wfn_ratio(is,js,xyzzyaaaa8,relprob=xyzzyaaab8,isnan=isnan,isinf=i&
&sinf)
endif
if(xyzzyaaaa8==1)then
call xyzzyaaad1(is,xyzzyaaac8,isnan,isinf)
call xyzzyaaad1(js,xyzzyaaad8,isnan,isinf)
xyzzyaaab8=xyzzyaaab8*xyzzyaaad8/xyzzyaaac8
endif
if(present(relprob))relprob=xyzzyaaab8
if(present(truncprob))truncprob=min(xyzzyaaab8,1.d0)
end subroutine wfn_altratio
subroutine xyzzyaaad1(is,fis,isnan,isinf)
implicit none
integer,intent(in) :: is
real(dp),intent(out) :: fis
logical,intent(out),optional :: isnan,isinf
real(dp) xyzzyaaaa9,xyzzyaaab9
logical nan,inf
if(present(isnan))isnan=.false.
if(present(isinf))isinf=.false.
if(f_valid(is))then
fis=f_scr(is)
return
endif
if(altsamp==1)then
call eval_local_energy(is,etot=xyzzyaaaa9,fix_nl_grid=.true.,isnan=nan&
&,isinf=inf,scr_nl_override=scr_nl_2)
xyzzyaaab9=abs(xyzzyaaaa9-vmc_optimum_e0)
fis=sqrt(xyzzyaaab9**2+2.d0*vmc_optimum_ew**2)
if(present(isnan))isnan=nan
if(present(isinf))isinf=inf
else
fis=1.d0
endif
f_scr(is)=fis
f_valid(is)=.true.
end subroutine xyzzyaaad1
subroutine alt_halflogpdf(is,logpval)
implicit none
integer,intent(in) :: is
real(dp),intent(out) :: logpval
real(dp) xyzzyaaaa10
complex(dp) xyzzyaaab10
logical iszero,isnan,isinf
if(simplepdf==1.and.isitcomplex)call errstop('ALT_HALFLOGPDF','Called &
&in complex.')
if(logp_valid(is))then
logpval=logp_scr(is)
return
endif
call wfn_logval(is,xyzzyaaab10,iszero=iszero)
if(simplepdf==2)then
call sampling_logval(is,xyzzyaaab10,iszero=iszero)
else
call wfn_logval(is,xyzzyaaab10,iszero=iszero)
endif
call xyzzyaaad1(is,xyzzyaaaa10,isnan=isnan,isinf=isinf)
logpval=dble(xyzzyaaab10)+0.5d0*log(xyzzyaaaa10)
logp_scr(is)=logpval
logp_valid(is)=.true.
end subroutine alt_halflogpdf
subroutine weight_config_vmc(is,weight)
implicit none
integer,intent(in) :: is
real(dp),intent(out) :: weight
complex(dp) xyzzyaaaa11,xyzzyaaab11
logical iszero
if(simplepdf==1.and..not.isitcomplex)call errstop('WEIGHT_CONFIG_VMC',&
&'Called in simple.')
if(.not.logp_valid(is))call errstop('WEIGHT_CONFIG_VMC','PDF appears t&
&o not be available in scratch.')
call wfn_logval(is,xyzzyaaab11,iszero=iszero)
if(iszero)then
weight=0.d0
else
xyzzyaaaa11=exp(xyzzyaaab11-cmplx(logp_scr(is),0.d0,dp))
weight=dble(xyzzyaaaa11)**2+aimag(xyzzyaaaa11)**2
endif
end subroutine weight_config_vmc
subroutine fieller(w,x,n,rn,writeout,dist)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: rn,x(:),w(:)
logical,intent(in) :: writeout,dist
real(dp) xyzzyaaaa12,xyzzyaaab12,xyzzyaaac12,xyzzyaaad12,xyzzyaaae12,x&
&yzzyaaaf12,xyzzyaaag12,xyzzyaaah12,xyzzyaaai12,xyzzyaaaj12,xyzzyaaak1&
&2,xyzzyaaal12,xyzzyaaam12,xyzzyaaan12,xyzzyaaao12,xyzzyaaap12,xyzzyaa&
&aq12,xyzzyaaar12,xyzzyaaas12,xyzzyaaat12,xyzzyaaau12,xyzzyaaav12,xyzz&
&yaaaw12
character(80) tmpr,tmpr2,tmpr3
xyzzyaaaf12=1.d0
xyzzyaaad12=1.d0/rn
xyzzyaaae12=1.d0/(rn*(rn-1.d0))
xyzzyaaag12=maxval(w(1:n))
xyzzyaaah12=minval(w(1:n))
if(dist)then
call mpi_reduce(xyzzyaaag12,xyzzyaaav12,1,mpi_double_precision,mpi_max&
&,0,mpi_comm_world,ierror)
call mpi_reduce(xyzzyaaah12,xyzzyaaaw12,1,mpi_double_precision,mpi_min&
&,0,mpi_comm_world,ierror)
else
xyzzyaaav12=xyzzyaaag12
xyzzyaaaw12=xyzzyaaah12
endif
if(am_master.and.writeout)then
tmpr=r2s(xyzzyaaav12,'(es12.4)')
tmpr2=r2s(xyzzyaaaw12,'(es12.4)')
tmpr3=r2s(xyzzyaaav12/xyzzyaaaw12,'(es12.4)')
call wout(' Max,Min,Max/Min weights : '//trim(tmpr)//', '//trim(tmpr2)&
&//', '//trim(tmpr3))
call wout()
endif
xyzzyaaah12=sum(x(1:n))
if(dist)then
call mpi_reduce(xyzzyaaah12,xyzzyaaaj12,1,mpi_double_precision,mpi_sum&
&,0,mpi_comm_world,ierror)
call mpi_reduce(xyzzyaaag12,xyzzyaaak12,1,mpi_double_precision,mpi_sum&
&,0,mpi_comm_world,ierror)
else
xyzzyaaaj12=xyzzyaaah12
xyzzyaaak12=xyzzyaaag12
endif
xyzzyaaak12=1.d0
xyzzyaaaj12=xyzzyaaaj12*xyzzyaaad12
if(dist)then
call mpi_bcast(xyzzyaaak12,1,mpi_double_precision,0,mpi_comm_world,ier&
&ror)
call mpi_bcast(xyzzyaaaj12,1,mpi_double_precision,0,mpi_comm_world,ier&
&ror)
endif
xyzzyaaaa12=xyzzyaaaj12
xyzzyaaag12=sum((x(1:n)-xyzzyaaaj12)**2)
xyzzyaaah12=sum((x(1:n)-xyzzyaaaj12)*(w(1:n)-xyzzyaaak12))
xyzzyaaai12=sum((w(1:n)-xyzzyaaak12)**2)
if(dist)then
call mpi_reduce(xyzzyaaag12,xyzzyaaal12,1,mpi_double_precision,mpi_sum&
&,0,mpi_comm_world,ierror)
call mpi_reduce(xyzzyaaah12,xyzzyaaam12,1,mpi_double_precision,mpi_sum&
&,0,mpi_comm_world,ierror)
call mpi_reduce(xyzzyaaai12,xyzzyaaan12,1,mpi_double_precision,mpi_sum&
&,0,mpi_comm_world,ierror)
else
xyzzyaaal12=xyzzyaaag12
xyzzyaaam12=xyzzyaaah12
xyzzyaaan12=xyzzyaaai12
endif
xyzzyaaan12=xyzzyaaan12*xyzzyaaae12
xyzzyaaam12=xyzzyaaam12*xyzzyaaae12
xyzzyaaal12=xyzzyaaal12*xyzzyaaae12
if(am_master)then
xyzzyaaao12=(xyzzyaaak12*xyzzyaaak12-xyzzyaaaf12**2*xyzzyaaan12)*0.5d0
xyzzyaaap12=-(xyzzyaaak12*xyzzyaaaj12-xyzzyaaaf12**2*xyzzyaaam12)
xyzzyaaaq12=(xyzzyaaaj12*xyzzyaaaj12-xyzzyaaaf12**2*xyzzyaaal12)*0.5d0
xyzzyaaag12=0.5d0/xyzzyaaao12
xyzzyaaas12=(-xyzzyaaap12-sqrt(xyzzyaaap12**2-4.0d0*xyzzyaaao12*xyzzya&
&aaq12))*xyzzyaaag12
xyzzyaaar12=(-xyzzyaaap12+sqrt(xyzzyaaap12**2-4.0d0*xyzzyaaao12*xyzzya&
&aaq12))*xyzzyaaag12
xyzzyaaac12=xyzzyaaas12
xyzzyaaab12=xyzzyaaar12
if(writeout)then
if(dist)then
tmpr=r2s(xyzzyaaaa12,'(f15.10)')
tmpr2=r2s(0.5d0*(xyzzyaaar12-xyzzyaaas12),'(f15.10)')
call wout(' Est[E_tot] (|l2-l1|/2)  : '//trim(tmpr)//' ( '//trim(tmpr2&
&)//' )')
tmpr=r2s(xyzzyaaas12,'(f15.10)')
tmpr2=r2s(xyzzyaaar12,'(f15.10)')
call wout(' (E)  Cnf. interval      : '//trim(tmpr)//'    '//trim(tmpr&
&2))
tmpr=r2s(xyzzyaaal12,'(es15.8)')
tmpr2=r2s(xyzzyaaam12,'(es15.8)')
tmpr3=r2s(xyzzyaaan12,'(es15.8)')
call wout(' (E)  c22,c12,c11        : '//trim(tmpr)//'    '//trim(tmpr&
&2)//'    '//trim(tmpr3))
else
tmpr=r2s(xyzzyaaaa12,'(f15.10)')
tmpr2=r2s(0.5d0*(xyzzyaaar12-xyzzyaaas12),'(f15.10)')
call wout('|Est[E_tot] (|l2-l1|/2)  : '//trim(tmpr)//' ( '//trim(tmpr2&
&)//' )')
tmpr=r2s(xyzzyaaas12,'(f15.10)')
tmpr2=r2s(xyzzyaaar12,'(f15.10)')
call wout('|(E)  Cnf. interval      : '//trim(tmpr)//'    '//trim(tmpr&
&2))
tmpr=r2s(xyzzyaaal12,'(es15.8)')
tmpr2=r2s(xyzzyaaam12,'(es15.8)')
tmpr3=r2s(xyzzyaaan12,'(es15.8)')
call wout('|(E)  c22,c12,c11        : '//trim(tmpr)//'    '//trim(tmpr&
&2)//'    '//trim(tmpr3))
endif
if(xyzzyaaap12**2<4.d0*xyzzyaaao12*xyzzyaaaq12)call wout('WARNING!    &
&             : Complex l1,l2 for Est[E_tot]')
call wout()
endif
endif
if(dist)then
call mpi_bcast(xyzzyaaac12,1,mpi_double_precision,0,mpi_comm_world,ier&
&ror)
call mpi_bcast(xyzzyaaab12,1,mpi_double_precision,0,mpi_comm_world,ier&
&ror)
endif
xyzzyaaah12=sum(w(1:n)*(x(1:n)/w(1:n)-xyzzyaaaa12)**2)
if(dist)then
call mpi_reduce(xyzzyaaah12,xyzzyaaaj12,1,mpi_double_precision,mpi_sum&
&,0,mpi_comm_world,ierror)
else
xyzzyaaaj12=xyzzyaaah12
endif
xyzzyaaag12=1.d0/(rn-1.d0)
xyzzyaaaj12=xyzzyaaaj12*xyzzyaaag12
xyzzyaaag12=1.d0/xyzzyaaak12
xyzzyaaat12=xyzzyaaaj12*xyzzyaaag12
if(dist)then
call mpi_bcast(xyzzyaaaj12,1,mpi_double_precision,0,mpi_comm_world,ier&
&ror)
call mpi_bcast(xyzzyaaat12,1,mpi_double_precision,0,mpi_comm_world,ier&
&ror)
endif
xyzzyaaah12=sum(w(1:n)*abs(x(1:n)/w(1:n)-xyzzyaaaa12))
if(dist)then
call mpi_reduce(xyzzyaaah12,xyzzyaaaj12,1,mpi_double_precision,mpi_sum&
&,0,mpi_comm_world,ierror)
else
xyzzyaaaj12=xyzzyaaah12
endif
xyzzyaaag12=1.d0/rn
xyzzyaaaj12=xyzzyaaaj12*xyzzyaaag12
xyzzyaaag12=1.d0/xyzzyaaak12
xyzzyaaau12=xyzzyaaaj12*xyzzyaaag12
if(dist)then
call mpi_bcast(xyzzyaaaj12,1,mpi_double_precision,0,mpi_comm_world,ier&
&ror)
call mpi_bcast(xyzzyaaau12,1,mpi_double_precision,0,mpi_comm_world,ier&
&ror)
endif
if(am_master)then
if(writeout)then
tmpr=r2s(xyzzyaaaf12/sqrt(rn)*sqrt(xyzzyaaat12),'(f15.10)')
tmpr2=r2s(xyzzyaaaf12/sqrt(rn)*xyzzyaaau12,'(f15.10)')
if(dist)then
call wout('    Std. E_tot    error  : '//trim(tmpr))
call wout('    Optimum E_tot error  : '//trim(tmpr2))
else
call wout('|   Std. E_tot    error  : '//trim(tmpr))
call wout('|   Optimum E_tot error  : '//trim(tmpr2))
endif
endif
endif
end subroutine fieller
end module slaarnaaa

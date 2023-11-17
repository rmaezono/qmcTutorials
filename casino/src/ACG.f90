module slaarnacg
use dsp
use slaarnach
use run_control, only : check_alloc,timer
use store,       only : netot,nele,nemax,nspin,complex_wf,which_ie,whi&
&ch_spin,which_ii,real1_complex2,pairing_wf
use slaarnacq, only : wfdet,wfdet_orbmask,wfdet_norb
implicit none
private
public setup_sampling_utils,finish_sampling_utils,sampling_ratio,sampl&
&ing_logval,sampling_accept_move,reset_config_sampling,clear_scratch_s&
&ampling
real(dp),allocatable :: xyzzyaaaa1(:,:)
logical,allocatable :: xyzzyaaab1(:,:)
real(dp),allocatable :: xyzzyaaac1(:,:,:,:),xyzzyaaad1(:,:)
logical,allocatable :: xyzzyaaae1(:,:,:,:),xyzzyaaaf1(:,:)
logical,allocatable :: xyzzyaaag1(:),xyzzyaaah1(:)
logical,parameter :: xyzzyaaai1=.true.
real(dp),parameter :: xyzzyaaaj1=.5d0
integer xyzzyaaak1
integer,allocatable :: xyzzyaaal1(:)
real(dp),allocatable,target :: xyzzyaaam1(:,:),xyzzyaaan1(:,:,:),xyzzy&
&aaao1(:,:)
contains
subroutine setup_sampling_utils
implicit none
integer xyzzyaaaa2,xyzzyaaab2
allocate(xyzzyaaaa1(netot,nscratch),xyzzyaaab1(netot,nscratch),stat=xy&
&zzyaaaa2)
call check_alloc(xyzzyaaaa2,'SETUP_SAMPLING_UTILS','orbmax_scr')
xyzzyaaaa1=0.d0
if(xyzzyaaai1)then
allocate(xyzzyaaac1(nemax,nemax,nspin,nscratch),xyzzyaaae1(nemax,nemax&
&,nspin,nscratch),xyzzyaaag1(nscratch),stat=xyzzyaaaa2)
call check_alloc(xyzzyaaaa2,'SETUP_SAMPLING_UTILS','coalcut_scr')
xyzzyaaac1=0.d0
xyzzyaaae1=.false.
allocate(xyzzyaaad1(nemax,nscratch),xyzzyaaaf1(nemax,nscratch),xyzzyaa&
&ah1(nscratch),stat=xyzzyaaaa2)
call check_alloc(xyzzyaaaa2,'SETUP_SAMPLING_UTILS','coalcut1_chscr')
xyzzyaaad1=0.d0
xyzzyaaaf1=.false.
endif
do xyzzyaaab2=1,nscratch
call clear_scratch_sampling(xyzzyaaab2)
enddo
call xyzzyaaap1
end subroutine setup_sampling_utils
subroutine finish_sampling_utils
implicit none
deallocate(xyzzyaaaa1,xyzzyaaab1)
if(xyzzyaaai1)then
deallocate(xyzzyaaac1,xyzzyaaae1,xyzzyaaag1)
deallocate(xyzzyaaad1,xyzzyaaaf1,xyzzyaaah1)
endif
call xyzzyaaaq1
end subroutine finish_sampling_utils
subroutine xyzzyaaap1
implicit none
integer xyzzyaaaa4
allocate(xyzzyaaam1(wfdet_norb,real1_complex2),xyzzyaaan1(3,wfdet_norb&
&,real1_complex2),xyzzyaaao1(wfdet_norb,real1_complex2),xyzzyaaal1(wfd&
&et_norb),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'ALLOCATE_TEMPS_SAMPLING','orbval_t')
xyzzyaaam1=0.d0
xyzzyaaan1=0.d0
xyzzyaaao1=0.d0
xyzzyaaal1=0
end subroutine xyzzyaaap1
subroutine xyzzyaaaq1
implicit none
deallocate(xyzzyaaam1,xyzzyaaan1,xyzzyaaao1,xyzzyaaal1)
end subroutine xyzzyaaaq1
subroutine sampling_ratio(is,js,ratio,relprob,truncprob,isinf,isnan)
implicit none
integer,intent(in) :: is,js
real(dp),intent(out),optional :: relprob,truncprob
complex(dp),intent(out),optional :: ratio
logical,intent(out),optional :: isnan,isinf
integer xyzzyaaaa6,xyzzyaaab6,xyzzyaaac6,xyzzyaaad6
real(dp) xyzzyaaae6,xyzzyaaaf6,xyzzyaaag6
logical xyzzyaaah6,xyzzyaaai6,xyzzyaaaj6,xyzzyaaak6
call timer('SAMPLING_RATIO',.true.)
xyzzyaaai6=.false.
xyzzyaaah6=.false.
if(buffer_move1_from(js)==is)then
xyzzyaaaa6=buffer_move1_from_ii(js)
call xyzzyaaar1(xyzzyaaaa6,is)
xyzzyaaae6=xyzzyaaaa1(xyzzyaaaa6,is)
call xyzzyaaar1(xyzzyaaaa6,js)
xyzzyaaaf6=xyzzyaaaa1(xyzzyaaaa6,js)
if(xyzzyaaae6==0.d0.and.xyzzyaaaf6==0.d0)then
xyzzyaaah6=.true.
elseif(xyzzyaaae6==0.d0)then
xyzzyaaai6=.true.
elseif(xyzzyaaaf6==0.d0)then
xyzzyaaag6=0.d0
else
xyzzyaaag6=xyzzyaaaf6/xyzzyaaae6
endif
if(xyzzyaaai1.and.xyzzyaaag6/=0.d0.and..not.(xyzzyaaah6.or.xyzzyaaai6)&
&)then
call xyzzyaaas1(is)
call xyzzyaaat1(xyzzyaaaa6,js)
xyzzyaaab6=which_ie(xyzzyaaaa6)
xyzzyaaad6=which_spin(xyzzyaaaa6)
do xyzzyaaac6=1,nele(xyzzyaaad6)
if(xyzzyaaac6==xyzzyaaab6)cycle
xyzzyaaae6=xyzzyaaac1(xyzzyaaac6,xyzzyaaab6,xyzzyaaad6,is)
xyzzyaaaj6=xyzzyaaae1(xyzzyaaac6,xyzzyaaab6,xyzzyaaad6,is)
xyzzyaaaf6=xyzzyaaad1(xyzzyaaac6,js)
xyzzyaaak6=xyzzyaaaf1(xyzzyaaac6,js)
if(xyzzyaaae6==0.d0.and.xyzzyaaaf6==0.d0)then
xyzzyaaah6=.true.
exit
elseif(xyzzyaaae6==0.d0)then
xyzzyaaai6=.true.
exit
elseif(xyzzyaaaf6==0.d0)then
xyzzyaaag6=0.d0
exit
elseif(xyzzyaaaj6.and.xyzzyaaak6)then
xyzzyaaag6=xyzzyaaag6*(xyzzyaaaf6/xyzzyaaae6)
elseif(xyzzyaaaj6)then
xyzzyaaag6=xyzzyaaag6/xyzzyaaae6
elseif(xyzzyaaak6)then
xyzzyaaag6=xyzzyaaag6*xyzzyaaaf6
endif
enddo
endif
else
xyzzyaaag6=1.d0
do xyzzyaaaa6=1,netot
call xyzzyaaar1(xyzzyaaaa6,is)
xyzzyaaae6=xyzzyaaaa1(xyzzyaaaa6,is)
call xyzzyaaar1(xyzzyaaaa6,js)
xyzzyaaaf6=xyzzyaaaa1(xyzzyaaaa6,js)
if(xyzzyaaae6==0.d0.and.xyzzyaaaf6==0.d0)then
xyzzyaaah6=.true.
exit
elseif(xyzzyaaae6==0.d0)then
xyzzyaaai6=.true.
exit
elseif(xyzzyaaaf6==0.d0)then
xyzzyaaag6=0.d0
exit
else
xyzzyaaag6=xyzzyaaag6*(xyzzyaaaf6/xyzzyaaae6)
endif
enddo
if(xyzzyaaai1.and.xyzzyaaag6/=0.d0.and..not.(xyzzyaaah6.or.xyzzyaaai6)&
&)then
call xyzzyaaas1(is)
call xyzzyaaas1(js)
ispin_loop: do xyzzyaaad6=1,nspin
do xyzzyaaab6=1,nele(xyzzyaaad6)
do xyzzyaaac6=xyzzyaaab6+1,nele(xyzzyaaad6)
xyzzyaaae6=xyzzyaaac1(xyzzyaaac6,xyzzyaaab6,xyzzyaaad6,is)
xyzzyaaaj6=xyzzyaaae1(xyzzyaaac6,xyzzyaaab6,xyzzyaaad6,is)
xyzzyaaaf6=xyzzyaaac1(xyzzyaaac6,xyzzyaaab6,xyzzyaaad6,js)
xyzzyaaak6=xyzzyaaae1(xyzzyaaac6,xyzzyaaab6,xyzzyaaad6,js)
if(xyzzyaaae6==0.d0.and.xyzzyaaaf6==0.d0)then
xyzzyaaah6=.true.
exit ispin_loop
elseif(xyzzyaaae6==0.d0)then
xyzzyaaai6=.true.
exit ispin_loop
elseif(xyzzyaaaf6==0.d0)then
xyzzyaaag6=0.d0
exit ispin_loop
elseif(xyzzyaaaj6.and.xyzzyaaak6)then
xyzzyaaag6=xyzzyaaag6*(xyzzyaaaf6/xyzzyaaae6)
elseif(xyzzyaaaj6)then
xyzzyaaag6=xyzzyaaag6/xyzzyaaae6
elseif(xyzzyaaak6)then
xyzzyaaag6=xyzzyaaag6*xyzzyaaaf6
endif
enddo
enddo
enddo ispin_loop
endif
endif
if(present(ratio))ratio=cmplx(xyzzyaaag6,0.d0,dp)
if(present(relprob))relprob=xyzzyaaag6**2
if(present(truncprob))truncprob=min(xyzzyaaag6**2,1.d0)
if(present(isnan))isnan=xyzzyaaah6
if(present(isinf))isinf=xyzzyaaai6
call timer('SAMPLING_RATIO',.false.)
end subroutine sampling_ratio
subroutine sampling_logval(is,logval,iszero)
implicit none
integer,intent(in) :: is
complex(dp),intent(out) :: logval
logical,intent(out),optional :: iszero
integer xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7,xyzzyaaad7
real(dp) xyzzyaaae7,xyzzyaaaf7
logical xyzzyaaag7
call timer('SAMPLING_LOGVAL',.true.)
xyzzyaaae7=0.d0
xyzzyaaag7=.false.
do xyzzyaaaa7=1,netot
call xyzzyaaar1(xyzzyaaaa7,is)
xyzzyaaaf7=xyzzyaaaa1(xyzzyaaaa7,is)
if(xyzzyaaaf7==0.d0)then
xyzzyaaag7=.true.
exit
else
xyzzyaaae7=xyzzyaaae7+log(xyzzyaaaf7)
endif
enddo
if(xyzzyaaai1.and..not.xyzzyaaag7)then
call xyzzyaaas1(is)
ispin_loop: do xyzzyaaac7=1,nspin
do xyzzyaaab7=1,nele(xyzzyaaac7)
do xyzzyaaad7=xyzzyaaab7+1,nele(xyzzyaaac7)
if(xyzzyaaae1(xyzzyaaad7,xyzzyaaab7,xyzzyaaac7,is))then
xyzzyaaaf7=xyzzyaaac1(xyzzyaaad7,xyzzyaaab7,xyzzyaaac7,is)
if(xyzzyaaaf7==0.d0)then
xyzzyaaag7=.true.
exit ispin_loop
else
xyzzyaaae7=xyzzyaaae7+log(xyzzyaaaf7)
endif
endif
enddo
enddo
enddo ispin_loop
endif
if(xyzzyaaag7)xyzzyaaae7=-1000.d0
logval=cmplx(xyzzyaaae7,0.d0,dp)
if(present(iszero))iszero=xyzzyaaag7
call timer('SAMPLING_LOGVAL',.false.)
end subroutine sampling_logval
subroutine sampling_accept_move(is,js)
implicit none
integer,intent(in) :: is,js
integer xyzzyaaaa8,xyzzyaaab8,xyzzyaaac8
call timer('SAMPLING_ACCEPT_MOVE',.true.)
if(buffer_move1_from(js)==is)then
xyzzyaaaa8=buffer_move1_from_ii(js)
if(xyzzyaaab1(xyzzyaaaa8,js))xyzzyaaaa1(xyzzyaaaa8,is)=xyzzyaaaa1(xyzz&
&yaaaa8,js)
xyzzyaaab1(xyzzyaaaa8,is)=xyzzyaaab1(xyzzyaaaa8,js)
if(xyzzyaaai1)then
if(xyzzyaaah1(js).and.xyzzyaaag1(is))then
xyzzyaaab8=which_ie(xyzzyaaaa8)
xyzzyaaac8=which_spin(xyzzyaaaa8)
xyzzyaaac1(1:nele(xyzzyaaac8),xyzzyaaab8,xyzzyaaac8,is)=xyzzyaaad1(1:n&
&ele(xyzzyaaac8),js)
xyzzyaaae1(1:nele(xyzzyaaac8),xyzzyaaab8,xyzzyaaac8,is)=xyzzyaaaf1(1:n&
&ele(xyzzyaaac8),js)
xyzzyaaac1(xyzzyaaab8,1:nele(xyzzyaaac8),xyzzyaaac8,is)=xyzzyaaad1(1:n&
&ele(xyzzyaaac8),js)
xyzzyaaae1(xyzzyaaab8,1:nele(xyzzyaaac8),xyzzyaaac8,is)=xyzzyaaaf1(1:n&
&ele(xyzzyaaac8),js)
xyzzyaaag1(is)=.true.
else
xyzzyaaag1(is)=.false.
endif
endif
else
do xyzzyaaaa8=1,netot
if(xyzzyaaab1(xyzzyaaaa8,js))xyzzyaaaa1(xyzzyaaaa8,is)=xyzzyaaaa1(xyzz&
&yaaaa8,js)
xyzzyaaab1(xyzzyaaaa8,is)=xyzzyaaab1(xyzzyaaaa8,js)
enddo
if(xyzzyaaai1)then
if(xyzzyaaag1(js))then
xyzzyaaac1(:,:,:,is)=xyzzyaaac1(:,:,:,js)
xyzzyaaae1(:,:,:,is)=xyzzyaaae1(:,:,:,js)
xyzzyaaag1(is)=.true.
else
xyzzyaaag1(is)=.false.
endif
xyzzyaaah1(is)=.false.
endif
endif
xyzzyaaab1(:,js)=.false.
if(xyzzyaaai1)then
xyzzyaaag1(js)=.false.
xyzzyaaah1(js)=.false.
endif
call timer('SAMPLING_ACCEPT_MOVE',.false.)
end subroutine sampling_accept_move
subroutine reset_config_sampling(js)
implicit none
integer,intent(in) :: js
call clear_scratch_sampling(js)
end subroutine reset_config_sampling
subroutine clear_scratch_sampling(is)
implicit none
integer,intent(in) :: is
xyzzyaaab1(:,is)=.false.
if(xyzzyaaai1)then
xyzzyaaag1(is)=.false.
xyzzyaaah1(is)=.false.
endif
end subroutine clear_scratch_sampling
subroutine xyzzyaaar1(ii,is)
implicit none
integer,intent(in) :: ii,is
integer xyzzyaaaa11,xyzzyaaab11
real(dp) xyzzyaaac11,xyzzyaaad11
if(xyzzyaaab1(ii,is))return
xyzzyaaab11=which_spin(ii)
if(buffer_move1_from(is)/=0)then
if(pairing_wf)call get_eevecs1_ch(ii,is)
call wfdet(rele1_chscr(1,is),sele1_chscr(is),xyzzyaaab11,wfdet_norb,wf&
&det_orbmask(1,xyzzyaaab11),.true.,.false.,xyzzyaaam1,xyzzyaaan1,xyzzy&
&aaao1,eevecs1_chscr(1,1,is),orb_m=xyzzyaaak1,orb_rmap=xyzzyaaal1)
else
call get_rsele(is)
if(pairing_wf)call get_eevecs(is)
call wfdet(rele_scr(1,ii,is),sele_scr(ii,is),xyzzyaaab11,wfdet_norb,wf&
&det_orbmask(1,xyzzyaaab11),.true.,.false.,xyzzyaaam1,xyzzyaaan1,xyzzy&
&aaao1,eevecs_scr(1,1,ii,is),orb_m=xyzzyaaak1,orb_rmap=xyzzyaaal1)
endif
xyzzyaaad11=0.d0
do xyzzyaaaa11=1,wfdet_norb
if(.not.wfdet_orbmask(xyzzyaaaa11,xyzzyaaab11))cycle
if(complex_wf)then
xyzzyaaac11=sqrt(xyzzyaaam1(xyzzyaaaa11,1)*xyzzyaaam1(xyzzyaaaa11,1)+x&
&yzzyaaam1(xyzzyaaaa11,2)*xyzzyaaam1(xyzzyaaaa11,2))
else
xyzzyaaac11=abs(xyzzyaaam1(xyzzyaaaa11,1))
endif
if(xyzzyaaac11>xyzzyaaad11)xyzzyaaad11=xyzzyaaac11
enddo
xyzzyaaaa1(ii,is)=xyzzyaaad11
xyzzyaaab1(ii,is)=.true.
end subroutine xyzzyaaar1
recursive subroutine xyzzyaaas1(xyzzyaaaa12)
implicit none
integer,intent(in) :: xyzzyaaaa12
integer xyzzyaaab12,xyzzyaaac12,xyzzyaaad12,xyzzyaaae12,xyzzyaaaf12,xy&
&zzyaaag12
real(dp) xyzzyaaah12
logical xyzzyaaai12,xyzzyaaaj12
if(xyzzyaaag1(xyzzyaaaa12))return
xyzzyaaai12=.false.
xyzzyaaab12=buffer_move1_from(xyzzyaaaa12)
if(xyzzyaaab12/=0)xyzzyaaai12=.true.
if(xyzzyaaai12)then
call xyzzyaaas1(xyzzyaaab12)
xyzzyaaac12=buffer_move1_from_ii(xyzzyaaaa12)
xyzzyaaae12=which_ie(xyzzyaaac12)
xyzzyaaag12=which_spin(xyzzyaaac12)
call xyzzyaaat1(xyzzyaaac12,xyzzyaaaa12)
xyzzyaaac1(:,:,:,xyzzyaaaa12)=xyzzyaaac1(:,:,:,xyzzyaaab12)
xyzzyaaae1(:,:,:,xyzzyaaaa12)=xyzzyaaae1(:,:,:,xyzzyaaab12)
xyzzyaaac1(1:nele(xyzzyaaag12),xyzzyaaae12,xyzzyaaag12,xyzzyaaaa12)=xy&
&zzyaaad1(1:nele(xyzzyaaag12),xyzzyaaaa12)
xyzzyaaae1(1:nele(xyzzyaaag12),xyzzyaaae12,xyzzyaaag12,xyzzyaaaa12)=xy&
&zzyaaaf1(1:nele(xyzzyaaag12),xyzzyaaaa12)
xyzzyaaac1(xyzzyaaae12,1:nele(xyzzyaaag12),xyzzyaaag12,xyzzyaaaa12)=xy&
&zzyaaad1(1:nele(xyzzyaaag12),xyzzyaaaa12)
xyzzyaaae1(xyzzyaaae12,1:nele(xyzzyaaag12),xyzzyaaag12,xyzzyaaaa12)=xy&
&zzyaaaf1(1:nele(xyzzyaaag12),xyzzyaaaa12)
else
call get_eevecs(xyzzyaaaa12)
do xyzzyaaac12=1,netot
xyzzyaaae12=which_ie(xyzzyaaac12)
xyzzyaaag12=which_spin(xyzzyaaac12)
xyzzyaaac1(xyzzyaaae12,xyzzyaaae12,xyzzyaaag12,xyzzyaaaa12)=1.d0
xyzzyaaae1(xyzzyaaae12,xyzzyaaae12,xyzzyaaag12,xyzzyaaaa12)=.false.
xyzzyaaad12=xyzzyaaac12
do xyzzyaaaf12=xyzzyaaae12+1,nele(xyzzyaaag12)
xyzzyaaad12=xyzzyaaad12+1
call xyzzyaaau1(eevecs_scr(4,xyzzyaaac12,xyzzyaaad12,xyzzyaaaa12),xyzz&
&yaaah12,xyzzyaaaj12)
xyzzyaaac1(xyzzyaaaf12,xyzzyaaae12,xyzzyaaag12,xyzzyaaaa12)=xyzzyaaah1&
&2
xyzzyaaae1(xyzzyaaaf12,xyzzyaaae12,xyzzyaaag12,xyzzyaaaa12)=xyzzyaaaj1&
&2
xyzzyaaac1(xyzzyaaae12,xyzzyaaaf12,xyzzyaaag12,xyzzyaaaa12)=xyzzyaaah1&
&2
xyzzyaaae1(xyzzyaaae12,xyzzyaaaf12,xyzzyaaag12,xyzzyaaaa12)=xyzzyaaaj1&
&2
enddo
enddo
endif
xyzzyaaag1(xyzzyaaaa12)=.true.
end subroutine xyzzyaaas1
subroutine xyzzyaaat1(ii,is)
implicit none
integer,intent(in) :: ii,is
integer xyzzyaaaa13,xyzzyaaab13,xyzzyaaac13,xyzzyaaad13
if(xyzzyaaah1(is))return
if(xyzzyaaag1(is))then
xyzzyaaab13=which_ie(ii)
xyzzyaaad13=which_spin(ii)
xyzzyaaad1(1:nele(xyzzyaaad13),is)=xyzzyaaac1(1:nele(xyzzyaaad13),xyzz&
&yaaab13,xyzzyaaad13,is)
xyzzyaaaf1(1:nele(xyzzyaaad13),is)=xyzzyaaae1(1:nele(xyzzyaaad13),xyzz&
&yaaab13,xyzzyaaad13,is)
else
xyzzyaaab13=which_ie(ii)
xyzzyaaad13=which_spin(ii)
call get_eevecs1_ch(ii,is)
xyzzyaaaa13=which_ii(1,xyzzyaaad13)-1
do xyzzyaaac13=1,nele(xyzzyaaad13)
xyzzyaaaa13=xyzzyaaaa13+1
if(xyzzyaaab13==xyzzyaaac13)then
xyzzyaaad1(xyzzyaaab13,is)=1.d0
xyzzyaaaf1(xyzzyaaab13,is)=.false.
else
call xyzzyaaau1(eevecs1_chscr(4,xyzzyaaaa13,is),xyzzyaaad1(xyzzyaaac13&
&,is),xyzzyaaaf1(xyzzyaaac13,is))
endif
enddo
endif
xyzzyaaah1(is)=.true.
end subroutine xyzzyaaat1
subroutine xyzzyaaau1(r,f,in_radius)
implicit none
real(dp),intent(in) :: r
real(dp),intent(out) :: f
logical,intent(out) :: in_radius
if(r>xyzzyaaaj1)then
f=1.d0
in_radius=.false.
else
f=r/xyzzyaaaj1
in_radius=.true.
endif
end subroutine xyzzyaaau1
end module slaarnacg

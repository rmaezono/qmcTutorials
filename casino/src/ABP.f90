module slaarnabp
use slaarnaag
use dsp
use parallel
use store
use file_utils,  only : open_units
use format_utils,only : wout,i2s
use slaarnabg,    only : atom_basis_type
use slaarnabt,   only : dcopy
use run_control, only : errstop,errstop_master,check_alloc
implicit none
private
public print_determinant_info,read_mdet,read_mdet_wfn,write_mdet,broad&
&cast_mdet,setup_mdet_params,finish_mdet_params,get_mdet_params,put_md&
&et_params
public mdet_max_mods,detcoef,wf_np,wf_nm,wf_nd,wf_d,wf_pm,mods,modifie&
&d_mdet,no_orb_phases,orb_phase_band,orb_phase_kpoint,orb_phase_spin,o&
&rb_phase_det,orb_phase,detcoef_label,detcoef_optable,orig_ndet,orig_d&
&etcoef,mdet_title
integer, parameter :: mdet_max_mods=10000
integer orig_ndet
real(dp),allocatable :: detcoef(:),orig_detcoef(:)
integer,allocatable :: wf_np(:,:),wf_nm(:,:),wf_nd(:,:)
integer,allocatable :: wf_pm(:,:,:,:),wf_d(:,:,:,:)
integer xyzzyaaaa1
integer xyzzyaaab1,xyzzyaaac1,xyzzyaaad1
logical xyzzyaaae1,modified_mdet,mods
logical xyzzyaaaf1
integer,parameter :: xyzzyaaag1=100
integer :: no_orb_phases=0
integer,allocatable :: orb_phase_band(:),orb_phase_kpoint(:),orb_phase&
&_spin(:),orb_phase_det(:)
real(dp),allocatable :: orb_phase(:)
integer xyzzyaaah1
integer,allocatable :: xyzzyaaai1(:),xyzzyaaaj1(:),xyzzyaaak1(:),detco&
&ef_label(:),detcoef_optable(:)
real(dp),allocatable :: xyzzyaaal1(:,:)
character(80) :: mdet_title='No title'
character(2) mdet_type
character(80),allocatable :: xyzzyaaam1(:)
real(dp),allocatable :: xyzzyaaan1(:,:)
contains
subroutine print_determinant_info
implicit none
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2
character(22) char22
character(80) tmpr
if(.not.am_master)return
if(atom_basis_type=='non_int_he')return
if(orig_ndet>1)then
call wout('Determinants :')
do xyzzyaaab2=1,nspin
call wout(' '//trim(i2s(xyzzyaaab2))//': '//trim(i2s(nele(xyzzyaaab2))&
&)//' x '//trim(pname(xyzzyaaab2)))
enddo
call wout(trim(i2s(orig_ndet))//' terms in multideterminant expansion.&
&')
else
call wout('Single determinant :')
endif
call wout()
if(dble(ndet)*dble(nspin+1)>1000.d0)then
call wout('(Output too long, omitted.  See correlation.data / xwfn.dat&
&a instead.)')
call wout()
return
endif
do xyzzyaaaa2=1,orig_ndet
do xyzzyaaab2=1,nspin
if(allocated(heg_orbtype))then
if(heg_orbtype(xyzzyaaab2,xyzzyaaaa2)<0.and.-heg_orbtype(xyzzyaaab2,xy&
&zzyaaaa2)<xyzzyaaab2)cycle
endif
char22='Det '//trim(i2s(xyzzyaaaa2))//' spin '//trim(i2s(xyzzyaaab2))
if(wf_np(xyzzyaaaa2,xyzzyaaab2)==0.and.wf_nm(xyzzyaaaa2,xyzzyaaab2)==0&
&.and.wf_nd(xyzzyaaaa2,xyzzyaaab2)==0)then
tmpr=char22//' : ground state'
call wout(tmpr)
else
if(wf_np(xyzzyaaaa2,xyzzyaaab2)/=0)then
do xyzzyaaac2=1,wf_np(xyzzyaaaa2,xyzzyaaab2)
tmpr=char22//' : e added to band '//trim(i2s(wf_pm(1,xyzzyaaac2,xyzzya&
&aaa2,xyzzyaaab2)))//' kp '//trim(i2s(wf_pm(2,xyzzyaaac2,xyzzyaaaa2,xy&
&zzyaaab2)))
call wout(tmpr)
enddo
endif
if(wf_nm(xyzzyaaaa2,xyzzyaaab2)/=0)then
do xyzzyaaac2=1,wf_nm(xyzzyaaaa2,xyzzyaaab2)
tmpr=char22//' : e removed from band '//trim(i2s(wf_pm(3,xyzzyaaac2,xy&
&zzyaaaa2,xyzzyaaab2)))//' kp '//trim(i2s(wf_pm(4,xyzzyaaac2,xyzzyaaaa&
&2,xyzzyaaab2)))
call wout(tmpr)
enddo
endif
if(wf_nd(xyzzyaaaa2,xyzzyaaab2)/=0)then
do xyzzyaaac2=1,wf_nd(xyzzyaaaa2,xyzzyaaab2)
tmpr=char22//' : e promoted from band '//trim(i2s(wf_d(1,xyzzyaaac2,xy&
&zzyaaaa2,xyzzyaaab2)))//' kp '//trim(i2s(wf_d(2,xyzzyaaac2,xyzzyaaaa2&
&,xyzzyaaab2)))//' to band '//trim(i2s(wf_d(3,xyzzyaaac2,xyzzyaaaa2,xy&
&zzyaaab2)))//' kp '//trim(i2s(wf_d(4,xyzzyaaac2,xyzzyaaaa2,xyzzyaaab2&
&)))
call wout(tmpr)
enddo
endif
endif
enddo
char22='Det '//trim(i2s(xyzzyaaaa2))//' coefficient'
write(tmpr,'(f16.8)')detcoef(xyzzyaaaa2)
tmpr=char22//' : '//trim(adjustl(tmpr))
call wout(tmpr)
enddo
call wout()
end subroutine print_determinant_info
subroutine read_mdet
implicit none
integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3,xyzzyaaad3,xyzzyaaae3
logical xyzzyaaaf3,xyzzyaaag3,xyzzyaaah3
character(80) rline,tmpr
modified_mdet=.false.
inquire(file="correlation.data",exist=xyzzyaaaf1)
if(xyzzyaaaf1)then
call open_units(xyzzyaaaa1,xyzzyaaab3)
if(xyzzyaaab3/=0)call errstop('READ_MDET','Unable to find free i/o uni&
&t.')
open(unit=xyzzyaaaa1,file='correlation.data',status='old',iostat=xyzzy&
&aaab3)
if(xyzzyaaab3/=0)call errstop('READ_MDET','Problem opening correlation&
&.data')
do
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaab3)rline
if(xyzzyaaab3>0.and.am_master)call errstop('READ_MDET','Problem readin&
&g correlation.data. Please check this file.')
if(xyzzyaaab3<0)then
xyzzyaaaf1=.false.
exit
endif
if(trim(adjustl(rline))=='START MDET')then
xyzzyaaaf1=.true.
exit
endif
enddo
if(.not.xyzzyaaaf1)close(xyzzyaaaa1)
endif
if(xyzzyaaaf1)then
if(am_master)then
call wout('Reading multideterminant/excitation data from correlation.d&
&ata.')
call wout('===========================================================&
&====')
endif
read(xyzzyaaaa1,*,err=1,end=1)
read(xyzzyaaaa1,'(a)',err=1,end=1)mdet_title
if(am_master)call wout('Title: '//trim(adjustl(mdet_title)))
read(xyzzyaaaa1,*,err=1,end=1)
read(xyzzyaaaa1,*,err=1,end=1)mdet_type
if(mdet_type=='GS')then
ndet=1
mods=.false.
xyzzyaaae1=.false.
if(am_master)call wout('Ground state as defined in xwfn.data will be u&
&sed.')
elseif(mdet_type=='SD')then
ndet=1
mods=.true.
xyzzyaaae1=.false.
modified_mdet=.true.
if(am_master)call wout('A single determinant with excitations will be &
&used.')
elseif(mdet_type=='MD')then
xyzzyaaae1=.true.
mods=.not.pairing_wf
modified_mdet=.true.
if(am_master)call wout('Multiple determinants will be used.')
read(xyzzyaaaa1,*,err=1,end=1)ndet
if(ndet<1.and.am_master)call errstop('READ_MDET','Number of determinan&
&ts should be positive in correlation.data.')
if(am_master)call wout('Number of determinants:  '//trim(i2s(ndet)))
elseif(am_master)then
call errstop('READ_MDET','First line of WAVEFUNCTION set in correlatio&
&n.data must be "GS", "SD" or "MD".')
endif
if(allocated(detcoef))deallocate(detcoef,wf_np,wf_nm,wf_nd)
allocate(detcoef(ndet),wf_np(ndet,nspin),wf_nm(ndet,nspin),wf_nd(ndet,&
&nspin),stat=xyzzyaaac3)
call check_alloc(xyzzyaaac3,'READ_MDET','1')
detcoef=1.d0
wf_np=0
wf_nm=0
wf_nd=0
if(mods)then
if(allocated(xyzzyaaam1))deallocate(xyzzyaaam1)
allocate(xyzzyaaam1(mdet_max_mods),stat=xyzzyaaac3)
call check_alloc(xyzzyaaac3,'READ_MDET','2.5')
endif
if(xyzzyaaae1)then
if(allocated(detcoef_label))deallocate(detcoef_label,detcoef_optable)
allocate(detcoef_label(ndet),detcoef_optable(ndet),stat=xyzzyaaac3)
call check_alloc(xyzzyaaac3,'READ_MDET','3')
xyzzyaaad3=1
xyzzyaaac1=0
xyzzyaaae3=0
xyzzyaaah3=ndet>1000
if(am_master)then
if(.not.xyzzyaaah3)then
call wout('Determinant ;   Coefficient ;     Label ; Optimizable')
else
call wout('(Output too long, omitted.  See correlation.data instead.)'&
&)
endif
endif
do xyzzyaaaa3=1,ndet
read(xyzzyaaaa1,*,err=1,end=1)detcoef(xyzzyaaaa3),detcoef_label(xyzzya&
&aaa3),detcoef_optable(xyzzyaaaa3)
if(am_master.and..not.xyzzyaaah3)then
if(detcoef_optable(xyzzyaaaa3)==0)then
write(tmpr,'(i8,2x,es19.10,3x,i5,5x,a)')xyzzyaaaa3,detcoef(xyzzyaaaa3)&
&,detcoef_label(xyzzyaaaa3),"Fixed"
elseif(detcoef_optable(xyzzyaaaa3)==1)then
write(tmpr,'(i8,2x,es19.10,3x,i5,5x,a)')xyzzyaaaa3,detcoef(xyzzyaaaa3)&
&,detcoef_label(xyzzyaaaa3),"Free"
else
call errstop('READ_MDET','Optimizable flag must be 0 (fixed) or 1 (fre&
&e).')
endif
call wout(tmpr)
endif
if(detcoef_label(xyzzyaaaa3)==xyzzyaaad3+1)then
if(xyzzyaaae3==0.and.am_master)call errstop('READ_MDET','No determinan&
&ts with label '//trim(i2s(xyzzyaaad3))//'.')
xyzzyaaad3=xyzzyaaad3+1
if(xyzzyaaae3>xyzzyaaac1)xyzzyaaac1=xyzzyaaae3
xyzzyaaae3=1
else
xyzzyaaae3=xyzzyaaae3+1
endif
if(detcoef_label(xyzzyaaaa3)/=xyzzyaaad3)call errstop('READ_MDET','Lab&
&els of determinant coefficients must be monotonically increasing natu&
&ral numbers.')
enddo
if(xyzzyaaae3==0.and.am_master)call errstop('READ_MDET','No determinan&
&ts with label '//trim(i2s(xyzzyaaad3))//'.')
if(xyzzyaaae3>xyzzyaaac1)xyzzyaaac1=xyzzyaaae3
xyzzyaaab1=xyzzyaaad3
if(am_master)call wout()
if(allocated(xyzzyaaai1))deallocate(xyzzyaaai1,xyzzyaaaj1,xyzzyaaak1,x&
&yzzyaaal1)
allocate(xyzzyaaai1(xyzzyaaab1),xyzzyaaaj1(xyzzyaaab1),xyzzyaaak1(xyzz&
&yaaab1),xyzzyaaal1(xyzzyaaac1,xyzzyaaab1),stat=xyzzyaaac3)
call check_alloc(xyzzyaaac3,'READ_MDET','4')
xyzzyaaal1(:,:)=0.d0
xyzzyaaad3=0
do xyzzyaaaa3=1,ndet
if(xyzzyaaad3/=detcoef_label(xyzzyaaaa3))then
xyzzyaaad3=detcoef_label(xyzzyaaaa3)
xyzzyaaak1(xyzzyaaad3)=1
xyzzyaaaj1(xyzzyaaad3)=xyzzyaaaa3
xyzzyaaal1(1,xyzzyaaad3)=1.d0
xyzzyaaai1(xyzzyaaad3)=detcoef_optable(xyzzyaaaa3)
else
xyzzyaaak1(xyzzyaaad3)=xyzzyaaak1(xyzzyaaad3)+1
if(detcoef(xyzzyaaaj1(xyzzyaaad3))/=0.d0)then
xyzzyaaal1(xyzzyaaak1(xyzzyaaad3),xyzzyaaad3)=detcoef(xyzzyaaaa3)/detc&
&oef(xyzzyaaaj1(xyzzyaaad3))
else
detcoef(xyzzyaaaa3)=0.d0
xyzzyaaal1(xyzzyaaak1(xyzzyaaad3),xyzzyaaad3)=0.d0
endif
if(xyzzyaaai1(xyzzyaaad3)/=detcoef_optable(xyzzyaaaa3).and.am_master)c&
&all errstop('READ_MDET','All determinant coefficients belonging to th&
&e same set must have the same optimizable flag.')
endif
enddo
if(count(xyzzyaaai1==1)>=xyzzyaaab1.and.am_master)call errstop('READ_M&
&DET','At least one set of determinant coefficients must be fixed.')
if(all(detcoef==0.d0).and.am_master)call errstop('READ_MDET','At least&
& one determinant coefficient must be nonzero.')
endif
if(mods)then
if(allocated(wf_pm))deallocate(wf_pm,wf_d)
allocate(wf_pm(4,mdet_max_mods,ndet,nspin),wf_d(4,mdet_max_mods,ndet,n&
&spin),stat=xyzzyaaac3)
call check_alloc(xyzzyaaac3,'READ_MDET','5')
wf_pm=0
wf_d=0
call xyzzyaaaq1(xyzzyaaaa1,.true.)
else
read(xyzzyaaaa1,'(a)',err=1,end=1)rline
if(am_master.and.trim(adjustl(rline))/='END MDET')call errstop('READ_M&
&DET','Was expecting to find "END MDET" in correlation.data.')
endif
close(xyzzyaaaa1)
else
if(am_master)then
call wout('No multideterminant/excitation data supplied.')
if(trim(atom_basis_type)=='none')then
call wout('Using ground state defined in input file.')
elseif(allocated(heg_orbtype))then
call wout('Using ground state defined in xwfn.data and input files.')
else
call wout('Using ground state defined in xwfn.data file.')
endif
endif
call xyzzyaaao1
endif
if(am_master)call wout()
detstart=1
detstop=ndet
if(use_altsamp.and.simplepdf==1)then
if(allocated(detcoef_label))then
xyzzyaaaf3=(detcoef_label(ndet)==1.or.ndet==1)
else
xyzzyaaaf3=ndet==1
endif
xyzzyaaag3=(.not.use_backflow.and..not.use_orbmods.and.ndet==1)
if(altsamp==2.and.xyzzyaaaf3.and..not.xyzzyaaag3)call errstop_master('&
&READ_MDET','Need more than 1 CSF for alt. sampling with SIMPLEPDF.')
if(ndet_smp==0)then
if(.not.xyzzyaaaf1)then
ndet_smp=1
else
if(altsamp==2)then
do xyzzyaaaa3=1,ndet-1
if(detcoef_label(xyzzyaaaa3)==2.and.detcoef_label(xyzzyaaaa3+1)==3)nde&
&t_smp=xyzzyaaaa3
enddo
if(detcoef_label(ndet)==2)ndet_smp=ndet
else
do xyzzyaaaa3=1,ndet-1
if(detcoef_label(xyzzyaaaa3)==1.and.detcoef_label(xyzzyaaaa3+1)==2)nde&
&t_smp=xyzzyaaaa3
enddo
if(detcoef_label(ndet)==1)ndet_smp=ndet
endif
endif
endif
if(ndet_smp==0)call errstop_master('READ_MDET','Problem with determini&
&ng wave function for SIMPLEPDF.')
endif
call xyzzyaaap1
return
1 if(am_master)call errstop('READ_MDET','Problem reading MDET block in&
& correlation.data.  Please check the file.')
end subroutine read_mdet
subroutine xyzzyaaao1
implicit none
integer xyzzyaaaa4
ndet=1
if(allocated(detcoef))deallocate(detcoef,wf_np,wf_nm,wf_nd)
allocate(detcoef(ndet),wf_np(ndet,nspin),wf_nm(ndet,nspin),wf_nd(ndet,&
&nspin),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'DEFAULT_MDET','detcoef')
detcoef=1.d0
wf_np=0
wf_nm=0
wf_nd=0
call xyzzyaaap1
end subroutine xyzzyaaao1
subroutine read_mdet_wfn(iunit)
implicit none
integer,intent(in) :: iunit
integer xyzzyaaaa5,xyzzyaaab5
character(80) rline
read(iunit,fmt='(a)',err=666,end=666)rline
if(trim(adjustl(rline))=='GS')then
ndet=1
mods=.false.
xyzzyaaae1=.false.
elseif(trim(adjustl(rline))=='SD')then
ndet=1
mods=.true.
xyzzyaaae1=.false.
modified_mdet=.false.
elseif(trim(adjustl(rline))=='MD')then
read(iunit,*,err=666,end=666)ndet
xyzzyaaae1=.true.
mods=.true.
modified_mdet=.false.
if(pairing_wf)mods=.false.
else
call errstop('READ_MDET_WFN','Unknown state-type in wave function file&
&: ' //trim(adjustl(rline))//'.')
endif
if(allocated(detcoef))deallocate(detcoef,wf_np,wf_nm,wf_nd)
allocate(detcoef(ndet),wf_np(ndet,nspin),wf_nm(ndet,nspin),wf_nd(ndet,&
&nspin),stat=xyzzyaaab5)
call check_alloc(xyzzyaaab5,'READ_MDET_WFN','1')
detcoef=1.d0
wf_np=0
wf_nm=0
wf_nd=0
detstart=1
detstop=ndet
if(pairing_wf.and.xyzzyaaae1)then
do xyzzyaaaa5=1,ndet
read(iunit,*,err=666,end=666)detcoef(xyzzyaaaa5)
enddo
endif
if(mods)then
if(allocated(wf_pm))deallocate(wf_pm,wf_d)
allocate(wf_pm(4,mdet_max_mods,ndet,nspin),wf_d(4,mdet_max_mods,ndet,n&
&spin),stat=xyzzyaaab5)
call check_alloc(xyzzyaaab5,'READ_MDET_WFN','3')
wf_pm=0
wf_d=0
if(xyzzyaaae1)then
do xyzzyaaaa5=1,ndet
read(iunit,*,err=666,end=666)detcoef(xyzzyaaaa5)
enddo
endif
call xyzzyaaaq1(iunit,.false.)
endif
call xyzzyaaap1
return
666 if(am_master)call errstop('READ_MDET_WFN','Problem reading multide&
&terminant data in xwfn.data.')
end subroutine read_mdet_wfn
subroutine broadcast_mdet
implicit none
integer xyzzyaaaa6
call mpi_bcast(ndet,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting ndet in broadcast_mdet')
call mpi_bcast(mods,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting mods in broadcast_mdet.')
call mpi_bcast(xyzzyaaae1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting scal in broadcast_mdet.')
if(.not.am_master)then
if(allocated(detcoef))deallocate(detcoef,wf_np,wf_nm,wf_nd)
allocate(detcoef(ndet),wf_np(ndet,nspin),wf_nm(ndet,nspin),wf_nd(ndet,&
&nspin),stat=xyzzyaaaa6)
call check_alloc(xyzzyaaaa6,'BROADCAST_MDET','1')
endif
call mpi_bcast(detcoef,ndet,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'Broadcasting detcoef in broadcast_mdet.')
call mpi_bcast(wf_nd,ndet*2,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting wf_nd in broadcast_mdet.')
call mpi_bcast(wf_np,ndet*2,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting wf_np in broadcast_mdet.')
call mpi_bcast(wf_nm,ndet*2,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting wf_nm in broadcast_mdet.')
call mpi_bcast(detstart,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting detstart in broadcast_mdet')
call mpi_bcast(detstop,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting detstop in broadcast_mdet')
if(mods)then
if(.not.am_master)then
if(allocated(wf_pm))deallocate(wf_pm,wf_d)
allocate(wf_pm(4,mdet_max_mods,ndet,nspin),wf_d(4,mdet_max_mods,ndet,n&
&spin),stat=xyzzyaaaa6)
call check_alloc(xyzzyaaaa6,'BROADCAST_MDET','3')
endif
call mpi_bcast(wf_pm,4*mdet_max_mods*ndet*nspin,mpi_integer,0,mpi_comm&
&_world,ierror)
call checkmpi(ierror,'Broadcasting wf_pm in broadcast_mdet.')
call mpi_bcast(wf_d,4*mdet_max_mods*ndet*nspin,mpi_integer,0,mpi_comm_&
&world,ierror)
call checkmpi(ierror,'Broadcasting wf_d in broadcast_mdet.')
endif
call xyzzyaaap1
end subroutine broadcast_mdet
subroutine xyzzyaaap1
implicit none
integer xyzzyaaaa7
if(allocated(orig_detcoef))deallocate(orig_detcoef)
orig_ndet=ndet
allocate(orig_detcoef(orig_ndet),stat=xyzzyaaaa7)
call check_alloc(xyzzyaaaa7,'READ_MDET','orig_detcoef')
orig_detcoef=detcoef
end subroutine xyzzyaaap1
subroutine xyzzyaaaq1(io_wfn,store_tokens)
implicit none
integer,intent(in) :: io_wfn
logical,intent(in) :: store_tokens
integer xyzzyaaaa8,xyzzyaaab8,xyzzyaaac8,xyzzyaaad8,xyzzyaaae8,xyzzyaa&
&af8,xyzzyaaag8,xyzzyaaah8,xyzzyaaai8
character(80) rline,rtok
if(store_tokens.and.am_master)call wout('Excitation specifications (se&
&e manual):')
xyzzyaaah8=0
do
xyzzyaaah8=xyzzyaaah8+1
read(io_wfn,fmt='(a)',err=666,end=666)rline
if(store_tokens)then
if(xyzzyaaah8>mdet_max_mods.and.am_master)call errstop('READ_TOKENS','&
&Too many tokens supplied. Please increase the mdet_max_mods parameter&
& in mdet.f90.')
xyzzyaaam1(xyzzyaaah8)=rline
endif
call xyzzyaaar1(rtok,rline)
if(rtok=='DET')then
call xyzzyaaar1(rtok,rline)
read(rtok,*)xyzzyaaai8
call xyzzyaaar1(rtok,rline)
read(rtok,*)xyzzyaaaa8
if((xyzzyaaaa8<1.or.xyzzyaaaa8>nspin.or.xyzzyaaai8<1.or.xyzzyaaai8>nde&
&t).and.am_master)then
call wout( 'Determinant no/spin out of range.')
call wout( 'Det '//trim(i2s(xyzzyaaai8))//'; Spin '//trim(i2s(xyzzyaaa&
&a8)))
call errstop('READ_TOKENS','Quitting.')
endif
call xyzzyaaar1(rtok,rline)
if(rtok=='PR')then
call xyzzyaaar1(rtok,rline)
read(rtok,*)xyzzyaaab8
call xyzzyaaar1(rtok,rline)
read(rtok,*)xyzzyaaac8
call xyzzyaaar1(rtok,rline)
read(rtok,*)xyzzyaaad8
call xyzzyaaar1(rtok,rline)
read(rtok,*)xyzzyaaae8
xyzzyaaaf8=wf_nd(xyzzyaaai8,xyzzyaaaa8)+1
wf_nd(xyzzyaaai8,xyzzyaaaa8)=xyzzyaaaf8
if(xyzzyaaaf8>mdet_max_mods.and.am_master)then
call wout('Too many promotions/additions/subtractions.')
call wout('Increase parameter MDET_MAX_MODS in mdet.f90.')
call errstop('READ_TOKENS','Quitting.')
endif
wf_d(1,xyzzyaaaf8,xyzzyaaai8,xyzzyaaaa8)=xyzzyaaab8
wf_d(2,xyzzyaaaf8,xyzzyaaai8,xyzzyaaaa8)=xyzzyaaac8
wf_d(3,xyzzyaaaf8,xyzzyaaai8,xyzzyaaaa8)=xyzzyaaad8
wf_d(4,xyzzyaaaf8,xyzzyaaai8,xyzzyaaaa8)=xyzzyaaae8
elseif(rtok=='PL')then
call xyzzyaaar1(rtok,rline)
read(rtok,*)xyzzyaaad8
call xyzzyaaar1(rtok,rline)
read(rtok,*)xyzzyaaae8
xyzzyaaaf8=wf_np(xyzzyaaai8,xyzzyaaaa8)+1
wf_np(xyzzyaaai8,xyzzyaaaa8)=xyzzyaaaf8
if(xyzzyaaaf8>mdet_max_mods.and.am_master)then
call wout('Too many promotions/additions/subtractions.')
call wout('Increase parameter MDET_MAX_MODS in mdet.f90.')
call errstop('READ_TOKENS','Quitting.')
endif
wf_pm(1,xyzzyaaaf8,xyzzyaaai8,xyzzyaaaa8)=xyzzyaaad8
wf_pm(2,xyzzyaaaf8,xyzzyaaai8,xyzzyaaaa8)=xyzzyaaae8
elseif(rtok=='MI')then
call xyzzyaaar1(rtok,rline)
read(rtok,*)xyzzyaaab8
call xyzzyaaar1(rtok,rline)
read(rtok,*)xyzzyaaac8
xyzzyaaaf8=wf_nm(xyzzyaaai8,xyzzyaaaa8)+1
wf_nm(xyzzyaaai8,xyzzyaaaa8)=xyzzyaaaf8
if(xyzzyaaaf8>mdet_max_mods.and.am_master)then
call wout('Too many promotions/additions/subtractions.')
call wout('Increase parameter MDET_MAX_MODS in mdet.f90.')
call errstop('READ_TOKENS','Quitting.')
endif
wf_pm(3,xyzzyaaaf8,xyzzyaaai8,xyzzyaaaa8)=xyzzyaaab8
wf_pm(4,xyzzyaaaf8,xyzzyaaai8,xyzzyaaaa8)=xyzzyaaac8
else
call errstop_master('READ_TOKENS','Unknown specifier '//trim(adjustl(r&
&tok)) //' in state spec. See READ_MDET for options.')
endif
elseif(rtok=='BLIP_PHASE'.or.rtok=="ORB_PHASE")then
if(.not.allocated(orb_phase))then
allocate(orb_phase_band(xyzzyaaag1),orb_phase_kpoint(xyzzyaaag1),orb_p&
&hase_spin(xyzzyaaag1),orb_phase_det(xyzzyaaag1),orb_phase(xyzzyaaag1)&
&,stat=xyzzyaaag8)
call check_alloc(xyzzyaaag8,'READ_TOKENS','11')
endif
no_orb_phases=no_orb_phases+1
if(no_orb_phases>xyzzyaaag1)call errstop_master('READ_TOKENS','Too man&
&y phases are specified.  Increase the max_orb_phases parameter in mde&
&t.f90.')
call xyzzyaaar1(rtok,rline)
read(rtok,*)orb_phase_det(no_orb_phases)
call xyzzyaaar1(rtok,rline)
read(rtok,*)orb_phase_spin(no_orb_phases)
call xyzzyaaar1(rtok,rline)
read(rtok,*)orb_phase_band(no_orb_phases)
call xyzzyaaar1(rtok,rline)
read(rtok,*)orb_phase_kpoint(no_orb_phases)
call xyzzyaaar1(rtok,rline)
read(rtok,*)orb_phase(no_orb_phases)
else
exit
endif
if(store_tokens.and.am_master)call wout(trim(adjustl(xyzzyaaam1(xyzzya&
&aah8))))
enddo
xyzzyaaad1=xyzzyaaah8-1
if(am_master.and.store_tokens.and.xyzzyaaad1<1)call wout('No excitatio&
&ns specifed.')
return
666 if(am_master)call errstop('READ_TOKENS','Problem reading tokens in&
& multideterminant data.')
end subroutine xyzzyaaaq1
subroutine xyzzyaaar1(t,s)
implicit none
character(len=*),intent(out) :: t
character(len=*),intent(inout) :: s
integer xyzzyaaaa9
s=adjustl(s)
xyzzyaaaa9=scan(s," ,")
t=s(1:xyzzyaaaa9)
s=s(xyzzyaaaa9:len(s))
end subroutine xyzzyaaar1
subroutine write_mdet(correlation_name)
implicit none
character(20),intent(in) :: correlation_name
integer xyzzyaaaa10,xyzzyaaab10
logical xyzzyaaac10
if(am_master)then
inquire(file=trim(correlation_name),exist=xyzzyaaac10)
if(xyzzyaaac10)then
open(unit=xyzzyaaaa1,file=trim(correlation_name),status='old',position&
&='append',iostat=xyzzyaaaa10)
else
open(unit=xyzzyaaaa1,file=trim(correlation_name),status='replace',iost&
&at=xyzzyaaaa10)
endif
if(xyzzyaaaa10/=0)call errstop('WRITE_MDET','Problem opening ' //trim(&
&correlation_name)//'.')
write(xyzzyaaaa1,*)'START MDET'
write(xyzzyaaaa1,*)'Title'
write(xyzzyaaaa1,*)trim(adjustl(mdet_title))
write(xyzzyaaaa1,*)'Multideterminant/excitation specification (see man&
&ual)'
write(xyzzyaaaa1,*)mdet_type
if(xyzzyaaae1)then
write(xyzzyaaaa1,'(3x,a)')trim(i2s(orig_ndet))//'                                  !&
& Number of dets'
do xyzzyaaab10=1,orig_ndet
write(xyzzyaaaa1,*)orig_detcoef(xyzzyaaab10),"  "//trim(i2s(detcoef_la&
&bel(xyzzyaaab10)))//"   "//trim(i2s(detcoef_optable(xyzzyaaab10)))//"      !&
& c_"//trim(i2s(xyzzyaaab10))//" ; label ; opt-flag"
enddo
endif
if(mods)then
do xyzzyaaab10=1,xyzzyaaad1
write(xyzzyaaaa1,*)trim(xyzzyaaam1(xyzzyaaab10))
enddo
endif
write(xyzzyaaaa1,*)'END MDET'
write(xyzzyaaaa1,*)
close(xyzzyaaaa1)
endif
end subroutine write_mdet
subroutine setup_mdet_params(nparam)
implicit none
integer,intent(inout) :: nparam
integer xyzzyaaaa11
nparam=0
xyzzyaaah1=0
if(.not.opt_det_coeff)return
if(modified_mdet.and.xyzzyaaae1)then
do xyzzyaaaa11=1,xyzzyaaab1
if(xyzzyaaai1(xyzzyaaaa11)==1)nparam=nparam+1
enddo
endif
xyzzyaaah1=nparam
call xyzzyaaas1
end subroutine setup_mdet_params
subroutine finish_mdet_params
implicit none
if(.not.opt_det_coeff)return
call xyzzyaaat1
end subroutine finish_mdet_params
subroutine xyzzyaaas1
implicit none
integer xyzzyaaaa13
allocate(xyzzyaaan1(orig_ndet,0:xyzzyaaah1),stat=xyzzyaaaa13)
call check_alloc(xyzzyaaaa13,'SETUP_MDET_PBUFFER','')
xyzzyaaan1=0.d0
end subroutine xyzzyaaas1
subroutine xyzzyaaat1
implicit none
deallocate(xyzzyaaan1)
end subroutine xyzzyaaat1
subroutine get_mdet_params(params,has_lolim,lolim,has_hilim,hilim,is_s&
&hallow,is_redundant,is_linear,is_loglinear,has_aderiv,affect_map,labe&
&l)
implicit none
real(dp),intent(inout) :: params(xyzzyaaah1),lolim(xyzzyaaah1),hilim(x&
&yzzyaaah1)
logical,intent(inout) :: has_lolim(xyzzyaaah1),has_hilim(xyzzyaaah1),i&
&s_shallow(xyzzyaaah1),is_redundant(xyzzyaaah1),is_linear(xyzzyaaah1),&
&is_loglinear(xyzzyaaah1),has_aderiv(xyzzyaaah1),affect_map(xyzzyaaah1&
&,xyzzyaaah1)
character(2),intent(inout) :: label(xyzzyaaah1)
integer xyzzyaaaa15,xyzzyaaab15
has_lolim=.false.
lolim=0.d0
has_hilim=.false.
hilim=0.d0
is_shallow=.false.
is_redundant=.false.
is_linear=.true.
is_loglinear=.false.
has_aderiv=.false.
affect_map=.false.
do xyzzyaaaa15=1,xyzzyaaah1
affect_map(xyzzyaaaa15,xyzzyaaaa15)=.true.
enddo
label='DC'
xyzzyaaaa15=0
do xyzzyaaab15=1,xyzzyaaab1
if(xyzzyaaai1(xyzzyaaab15)==1)then
xyzzyaaaa15=xyzzyaaaa15+1
params(xyzzyaaaa15)=orig_detcoef(xyzzyaaaj1(xyzzyaaab15))
endif
enddo
end subroutine get_mdet_params
subroutine put_mdet_params(params,ignore,iparam_buffer,prestore,bad_pa&
&rams)
implicit none
integer,intent(in) :: iparam_buffer
real(dp),intent(inout) :: params(xyzzyaaah1)
logical,intent(in) :: ignore(xyzzyaaah1),prestore
logical,intent(out) :: bad_params
integer xyzzyaaaa16,xyzzyaaab16,xyzzyaaac16,xyzzyaaad16
bad_params=.false.
if(prestore)then
call xyzzyaaav1(iparam_buffer)
return
endif
xyzzyaaaa16=0
do xyzzyaaab16=1,xyzzyaaab1
if(xyzzyaaai1(xyzzyaaab16)==1)then
xyzzyaaaa16=xyzzyaaaa16+1
if(.not.ignore(xyzzyaaaa16))then
orig_detcoef(xyzzyaaaj1(xyzzyaaab16))=params(xyzzyaaaa16)
xyzzyaaac16=1
do xyzzyaaad16=xyzzyaaaj1(xyzzyaaab16)+1,xyzzyaaaj1(xyzzyaaab16)+xyzzy&
&aaak1(xyzzyaaab16)-1
xyzzyaaac16=xyzzyaaac16+1
orig_detcoef(xyzzyaaad16)=xyzzyaaal1(xyzzyaaac16,xyzzyaaab16)*params(x&
&yzzyaaaa16)
enddo
endif
endif
enddo
call xyzzyaaau1(iparam_buffer)
end subroutine put_mdet_params
subroutine xyzzyaaau1(indx)
implicit none
integer,intent(in) :: indx
call dcopy(orig_ndet,orig_detcoef(1),1,xyzzyaaan1(1,indx),1)
end subroutine xyzzyaaau1
subroutine xyzzyaaav1(indx)
implicit none
integer,intent(in) :: indx
call dcopy(orig_ndet,xyzzyaaan1(1,indx),1,orig_detcoef(1),1)
end subroutine xyzzyaaav1
end module slaarnabp

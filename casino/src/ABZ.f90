module slaarnabyter
use dsp
use slaarnaam
use parallel
use slaarnach
use slaarnacs
use slaarnaai,      only : wfdet_s
use slaarnaan, only : ee_distances,ei_distances,points,map_to_simcell
use slaarnaap,         only : eval_expot
use file_utils,    only : open_units
use format_utils,  only : wout,i2s,wordwrap,capitalize
use slaarnaat,     only : spin_dep_gs,s_plot
use slaarnabg,      only : nitot,isperiodic,atom_basis_type,model_syst&
&em,npcells,rion,dimensionality
use run_control,   only : errstop,errstop_master,errwarn,timer,check_a&
&lloc
use store,         only : netot,nemax,real1_complex2,nspin,nele,ndet,p&
&airing_wf,nuc_nele,three_netot,heg_orbtype,use_jastrow,open_unit,use_&
&backflow,vmc_cfg_by_cfg,no_difftypes,which_spin,which_ssingle
use slaarnacp,           only : equilibration
use slaarnacq,   only : wfdet,wfdet_norb,wfdet_orbmap,wfdet_orbmask,co&
&py_orb_to_det
implicit none
private
public init_plotter,setup_plotter,qmc_plotter
integer,parameter :: nthings_to_plot=10,xyzzyaaaa1=1,xyzzyaaab1=2,xyzz&
&yaaac1=3,xyzzyaaad1=4,xyzzyaaae1=5,xyzzyaaaf1=6,xyzzyaaag1=7,xyzzyaaa&
&h1=8,xyzzyaaai1=9,xyzzyaaaj1=10
character(20),parameter :: xyzzyaaak1(nthings_to_plot)=(/'orb      ','&
&orb_gradx','orb_grady','orb_gradz','orb_lap  ',  'eipot    ','energy &
&  ','wfn      ','nodes    ','expot    '/)
character(64),parameter :: plot_name(nthings_to_plot)=(/'selected orbi&
&tals              ','x-gradient of selected orbitals',  'y-gradient o&
&f selected orbitals','z-gradient of selected orbitals',  'Laplacian o&
&f selected orbitals ','electron-ion local potential   ',  'local ener&
&gy                   ','wave function                  ',  'wave func&
&tion and nodes        ','external potential             '/)
character(64),parameter :: xyzzyaaal1(3)=(/'lineplot.dat','2Dplot.dat &
& ','3Dplot.dat  '/)
character(64),parameter :: xyzzyaaam1(3)=(/'1Dnodes.dat','2Dnodes.dat'&
&,'3Dnodes.dat'/)
integer :: xyzzyaaan1=0
character(180),allocatable :: xyzzyaaao1(:)
integer what_to_plot,plot_dim,grid_npoint_ab,grid_npoint_ac,grid_npoin&
&t_ad,plot_nfunctions,plot_ncomponents,plot_nfix_elec
integer,allocatable :: xyzzyaaap1(:),xyzzyaaaq1(:),xyzzyaaar1(:),xyzzy&
&aaas1(:)
real(dp) xyzzyaaat1(3),xyzzyaaau1(3),xyzzyaaav1(3),d(3)
real(dp),allocatable :: xyzzyaaaw1(:,:)
contains
subroutine init_plotter(nlines,input_block)
implicit none
integer,intent(in) :: nlines
character(180),intent(in) :: input_block(nlines)
integer xyzzyaaaa2,xyzzyaaab2
if(xyzzyaaan1>0.or.allocated(xyzzyaaao1))call errstop_master('INIT_PLO&
&TTER','Reading QMC_PLOT input block twice.')
xyzzyaaan1=nlines
if(nlines<1)return
allocate(xyzzyaaao1(xyzzyaaan1),stat=xyzzyaaaa2)
call check_alloc(xyzzyaaaa2,'INIT_PLOTTER','qmc_plot_block')
do xyzzyaaab2=1,nlines
xyzzyaaao1(xyzzyaaab2)=trim(adjustl(input_block(xyzzyaaab2)))
enddo
end subroutine init_plotter
subroutine setup_plotter
implicit none
integer xyzzyaaaa3,xyzzyaaab3,indx,ialloc,ierr
real(dp) xyzzyaaac3,xyzzyaaad3(3),xyzzyaaae3(3),xyzzyaaaf3(3),xyzzyaaa&
&g3(3)
real(dp),parameter :: xyzzyaaah3=1.d-13
character(20) what_to_plot_string,plot_dim_string
if(xyzzyaaan1<1.or..not.allocated(xyzzyaaao1))call errstop_master('SET&
&UP_PLOTTER','QMC_PLOT input block not found in input file, or not rec&
&eived in INIT_PLOTTER.')
what_to_plot=0
grid_npoint_ab=0
grid_npoint_ac=0
grid_npoint_ad=0
xyzzyaaat1=0.d0
xyzzyaaau1=0.d0
xyzzyaaav1=0.d0
d=0.d0
plot_nfunctions=1
plot_ncomponents=1
indx=0
call xyzzyaaai3
read(xyzzyaaao1(indx),*,iostat=ierr)what_to_plot_string
call xyzzyaaaj3
what_to_plot_string=trim(adjustl(what_to_plot_string))
call capitalize(what_to_plot_string,decapitalize_in=.true.)
do xyzzyaaaa3=1,nthings_to_plot
if(trim(what_to_plot_string)==trim(xyzzyaaak1(xyzzyaaaa3)))then
what_to_plot=xyzzyaaaa3
exit
endif
enddo
if(what_to_plot==0)then
if(am_master)then
call wout('The first line in the QMC_PLOT input block must be one of:'&
&)
do xyzzyaaaa3=1,nthings_to_plot
call wout(' - '//trim(xyzzyaaak1(xyzzyaaaa3)))
enddo
endif
call errstop_master('SETUP_PLOTTER','Quitting.')
endif
call xyzzyaaai3
read(xyzzyaaao1(indx),*,iostat=ierr)plot_dim_string
call xyzzyaaaj3
call capitalize(plot_dim_string,decapitalize_in=.true.)
select case(trim(adjustl(plot_dim_string)))
case('line','1','1d','1-d')
plot_dim=1
case('plane','2','2d','2-d')
plot_dim=2
case('volume','3','3d','3-d')
plot_dim=3
case default
if(am_master)then
call wout('The second line in the qmc_plot block must be one of:')
call wout(' - line')
call wout(' - plane')
call wout(' - volume')
endif
call errstop_master('SETUP_PLOTTER','Stopping.')
end select
select case(plot_dim)
case(1)
call xyzzyaaai3
read(xyzzyaaao1(indx),*,iostat=ierr)grid_npoint_ab
call xyzzyaaaj3
call xyzzyaaai3
read(xyzzyaaao1(indx),*,iostat=ierr)xyzzyaaat1(1:dimensionality)
call xyzzyaaaj3
call xyzzyaaai3
read(xyzzyaaao1(indx),*,iostat=ierr)xyzzyaaau1(1:dimensionality)
call xyzzyaaaj3
case(2)
call xyzzyaaai3
read(xyzzyaaao1(indx),*,iostat=ierr)grid_npoint_ab,grid_npoint_ac
call xyzzyaaaj3
call xyzzyaaai3
read(xyzzyaaao1(indx),*,iostat=ierr)xyzzyaaat1(1:dimensionality)
call xyzzyaaaj3
call xyzzyaaai3
read(xyzzyaaao1(indx),*,iostat=ierr)xyzzyaaau1(1:dimensionality)
call xyzzyaaaj3
call xyzzyaaai3
read(xyzzyaaao1(indx),*,iostat=ierr)xyzzyaaav1(1:dimensionality)
call xyzzyaaaj3
case(3)
call xyzzyaaai3
read(xyzzyaaao1(indx),*,iostat=ierr)grid_npoint_ab,grid_npoint_ac,grid&
&_npoint_ad
call xyzzyaaaj3
call xyzzyaaai3
read(xyzzyaaao1(indx),*,iostat=ierr)xyzzyaaat1(1:dimensionality)
call xyzzyaaaj3
call xyzzyaaai3
read(xyzzyaaao1(indx),*,iostat=ierr)xyzzyaaau1(1:dimensionality)
call xyzzyaaaj3
call xyzzyaaai3
read(xyzzyaaao1(indx),*,iostat=ierr)xyzzyaaav1(1:dimensionality)
call xyzzyaaaj3
call xyzzyaaai3
read(xyzzyaaao1(indx),*,iostat=ierr)d(1:dimensionality)
call xyzzyaaaj3
end select
xyzzyaaad3=xyzzyaaat1-xyzzyaaau1
xyzzyaaac3=sqrt(sum(xyzzyaaad3(:)**2))
if(xyzzyaaac3==0.d0)call errstop_master('SETUP_PLOTTER','Zero line len&
&gth AB in QMC_PLOT input block.')
if(plot_dim>1)then
xyzzyaaae3=xyzzyaaat1-xyzzyaaav1
xyzzyaaac3=sqrt(sum(xyzzyaaae3(:)**2))
if(xyzzyaaac3==0.d0)call errstop_master('SETUP_PLOTTER','Zero line len&
&gth AC in QMC_PLOT input block.')
xyzzyaaag3(1)=xyzzyaaad3(2)*xyzzyaaae3(3)-xyzzyaaad3(3)*xyzzyaaae3(2)
xyzzyaaag3(2)=xyzzyaaad3(3)*xyzzyaaae3(1)-xyzzyaaad3(1)*xyzzyaaae3(3)
xyzzyaaag3(3)=xyzzyaaad3(1)*xyzzyaaae3(2)-xyzzyaaad3(2)*xyzzyaaae3(1)
if(all(abs(xyzzyaaag3)<xyzzyaaah3))call errstop_master('SETUP_PLOTTER'&
&,'Points A, B and C in QMC_PLOT input block are collinear.')
endif
if(plot_dim>2)then
xyzzyaaaf3=xyzzyaaat1-d
xyzzyaaac3=sqrt(sum(xyzzyaaaf3(:)**2))
if(xyzzyaaac3==0.d0)call errstop_master('SETUP_PLOTTER','Zero line len&
&gth AD for qmc_plot.')
xyzzyaaag3(1)=xyzzyaaad3(2)*xyzzyaaaf3(3)-xyzzyaaad3(3)*xyzzyaaaf3(2)
xyzzyaaag3(2)=xyzzyaaad3(3)*xyzzyaaaf3(1)-xyzzyaaad3(1)*xyzzyaaaf3(3)
xyzzyaaag3(3)=xyzzyaaad3(1)*xyzzyaaaf3(2)-xyzzyaaad3(2)*xyzzyaaaf3(1)
if(all(abs(xyzzyaaag3)<xyzzyaaah3))call errstop_master('SETUP_PLOTTER'&
&,'Points A, B and D in QMC_PLOT input block are collinear.')
xyzzyaaag3(1)=xyzzyaaae3(2)*xyzzyaaaf3(3)-xyzzyaaae3(3)*xyzzyaaaf3(2)
xyzzyaaag3(2)=xyzzyaaae3(3)*xyzzyaaaf3(1)-xyzzyaaae3(1)*xyzzyaaaf3(3)
xyzzyaaag3(3)=xyzzyaaae3(1)*xyzzyaaaf3(2)-xyzzyaaae3(2)*xyzzyaaaf3(1)
if(all(abs(xyzzyaaag3)<xyzzyaaah3))call errstop_master('SETUP_PLOTTER'&
&,'Points A, C and D in QMC_PLOT input block are collinear.')
endif
select case(what_to_plot)
case(xyzzyaaaa1,xyzzyaaab1,xyzzyaaac1,xyzzyaaad1,xyzzyaaae1)
call xyzzyaaai3
read(xyzzyaaao1(indx),*,iostat=ierr)plot_nfunctions
call xyzzyaaaj3
if(plot_nfunctions<1)call errstop_master('SETUP_PLOTTER','Number of or&
&bitals to plot is zero in the QMC_PLOT input block. Quitting.')
allocate(xyzzyaaap1(plot_nfunctions),xyzzyaaaq1(plot_nfunctions),stat=&
&ialloc)
call check_alloc(ialloc,'SETUP_PLOTTER','plot_orb_indx')
xyzzyaaap1=0
xyzzyaaaq1=0
call xyzzyaaai3
read(xyzzyaaao1(indx),*,iostat=ierr)xyzzyaaap1(1:plot_nfunctions)
call xyzzyaaaj3
call xyzzyaaai3
read(xyzzyaaao1(indx),*,iostat=ierr)xyzzyaaaq1(1:plot_nfunctions)
call xyzzyaaaj3
case(xyzzyaaag1,xyzzyaaah1,xyzzyaaai1)
call xyzzyaaai3
read(xyzzyaaao1(indx),*,iostat=ierr)plot_nfix_elec
call xyzzyaaaj3
allocate(xyzzyaaar1(plot_nfix_elec),xyzzyaaas1(plot_nfix_elec),xyzzyaa&
&aw1(3,plot_nfix_elec),stat=ialloc)
call check_alloc(ialloc,'SETUP_PLOTTER','plot_fix_rele')
do xyzzyaaab3=1,plot_nfix_elec
call xyzzyaaai3
read(xyzzyaaao1(indx),*,iostat=ierr)xyzzyaaar1(xyzzyaaab3),xyzzyaaas1(&
&xyzzyaaab3),xyzzyaaaw1(1:dimensionality,xyzzyaaab3)
call xyzzyaaaj3
enddo
end select
if(xyzzyaaan1/=indx)call errstop_master('SETUP_PLOTTER','Unexpected ex&
&tra lines at the end of the QMC_PLOT input block.')
deallocate(xyzzyaaao1)
xyzzyaaan1=0
contains
subroutine xyzzyaaai3()
implicit none
indx=indx+1
if(xyzzyaaan1<indx)call errstop_master('SETUP_PLOTTER','Expecting to f&
&ind at least '//trim(i2s(indx))//' lines in the QMC_PLOT input block.&
&')
end subroutine xyzzyaaai3
subroutine xyzzyaaaj3()
implicit none
if(ierr<0)then
call errstop_master('SETUP_PLOTTER','Problem reading line '//trim(i2s(&
&indx))//' in the QMC_PLOTTER input block.')
elseif(ierr>0)then
call errstop_master('SETUP_PLOTTER','Reached end of line '//trim(i2s(i&
&ndx))//' in the QMC_PLOTTER input block without finding all the data &
&it should contain.')
endif
end subroutine xyzzyaaaj3
end subroutine setup_plotter
subroutine qmc_plotter(nequil,dtvmc_array,dtvmc_shift_array,opt_dtvmc_&
&array,initial_rele,initial_rele_set)
implicit none
integer,intent(in) :: nequil,opt_dtvmc_array(no_difftypes)
real(dp),intent(in) :: initial_rele(3,netot)
real(dp),intent(inout) :: dtvmc_array(no_difftypes),dtvmc_shift_array(&
&no_difftypes)
logical,intent(in) :: initial_rele_set(netot)
integer io,xyzzyaaaa6,xyzzyaaab6,xyzzyaaac6,xyzzyaaad6,xyzzyaaae6,xyzz&
&yaaaf6,xyzzyaaag6,xyzzyaaah6,tot_npoint,xyzzyaaai6,xyzzyaaaj6,xyzzyaa&
&ak6,xyzzyaaal6,xyzzyaaam6(2,2,2),ifn,xyzzyaaan6,xyzzyaaao6,iplotmin,i&
&plotmax,his_node,xyzzyaaap6,ialloc,xyzzyaaaq6
integer,allocatable :: map_point(:,:,:),xyzzyaaar6(:)
real(dp) xyzzyaaas6,xyzzyaaat6,xyzzyaaau6(3),xyzzyaaav6(3),xyzzyaaaw6(&
&3),xyzzyaaax6,xyzzyaaay6,xyzzyaaaz6(3),xyzzyaaba6,etot,xyzzyaabb6,xyz&
&zyaabc6(2),xyzzyaabd6(2),xyzzyaabe6(2),xyzzyaabf6(no_difftypes)
real(dp) :: xyzzyaabg6=1.d0
real(dp),allocatable :: grid_point(:,:),xyzzyaabh6(:,:),xyzzyaabi6(:,:&
&),xyzzyaabj6(:,:,:),xyzzyaabk6(:,:,:),xyzzyaabl6(:,:,:,:),xyzzyaabm6(&
&:,:,:),xyzzyaabn6(:,:),plot_data(:,:,:),xyzzyaabo6(:,:),xyzzyaabp6(:,&
&:,:),xyzzyaabq6(:,:)
complex(dp) xyzzyaabr6
logical xyzzyaabs6,xyzzyaabt6,xyzzyaabu6,xyzzyaabv6,xyzzyaabw6,isnan,i&
&sinf,xyzzyaabx6,xyzzyaaby6
logical,allocatable :: valid_plot_point(:)
character(64) filename
character(80) tmpr
call timer('QMC_PLOTTER',.true.)
xyzzyaaau6(1:3)=(xyzzyaaau1(1:3)-xyzzyaaat1(1:3))/dble(grid_npoint_ab-&
&1)
if(plot_dim>1)then
xyzzyaaav6(1:3)=(xyzzyaaav1(1:3)-xyzzyaaat1(1:3))/dble(grid_npoint_ac-&
&1)
else
xyzzyaaav6=0.d0
grid_npoint_ac=1
endif
if(plot_dim>2)then
xyzzyaaaw6(1:3)=(d(1:3)-xyzzyaaat1(1:3))/dble(grid_npoint_ad-1)
else
xyzzyaaaw6=0.d0
grid_npoint_ad=1
endif
tot_npoint=grid_npoint_ab*grid_npoint_ac*grid_npoint_ad
allocate(plot_data(plot_ncomponents,tot_npoint,plot_nfunctions),valid_&
&plot_point(tot_npoint),grid_point(3,tot_npoint),map_point(grid_npoint&
&_ab,grid_npoint_ac,grid_npoint_ad),stat=ialloc)
call check_alloc(ialloc,'QMC_PLOTTER','plot_data')
valid_plot_point=.true.
xyzzyaaai6=0
do xyzzyaaal6=1,grid_npoint_ad
do xyzzyaaak6=1,grid_npoint_ac
do xyzzyaaaj6=1,grid_npoint_ab
xyzzyaaai6=xyzzyaaai6+1
map_point(xyzzyaaaj6,xyzzyaaak6,xyzzyaaal6)=xyzzyaaai6
grid_point(1:3,xyzzyaaai6)=xyzzyaaat1(1:3)+real(xyzzyaaaj6-1,dp)*xyzzy&
&aaau6(1:3)+real(xyzzyaaak6-1,dp)*xyzzyaaav6(1:3)+real(xyzzyaaal6-1,dp&
&)*xyzzyaaaw6(1:3)
enddo
enddo
enddo
if(tot_npoint<nnodes)call errstop_master('QMC_PLOTTER','Number of poin&
&ts to plot is less than the number of processors.')
iplotmin=nint(dble(tot_npoint)/dble(nnodes)*dble(my_node))+1
if(my_node==nnodes-1)then
iplotmax=tot_npoint
else
iplotmax=nint(dble(tot_npoint)/dble(nnodes)*dble(my_node+1))
endif
if(am_master)then
filename=xyzzyaaal1(plot_dim)
call wout(trim(i2s(plot_dim))//'-DIMENSIONAL QMC_PLOT REQUESTED')
call wout('================================')
call wout()
call wout("Plotting "//trim(plot_name(what_to_plot))//" to '"//trim(fi&
&lename)//"'.")
call wout(trim(i2s(tot_npoint))//' points to scan.')
call wout()
endif
allocate(xyzzyaabi6(3,netot),xyzzyaaar6(netot),stat=ialloc)
call check_alloc(ialloc,'QMC_PLOTTER','rele')
xyzzyaabi6=0.d0
xyzzyaaar6=0
select case(what_to_plot)
case(xyzzyaaaa1,xyzzyaaab1,xyzzyaaac1,xyzzyaaad1,xyzzyaaae1)
if(any(xyzzyaaap1(:)>nemax).or.any(xyzzyaaap1(:)<1))then
if(am_master)call wordwrap('Some of the orbital sequence numbers for t&
&he line orbital plot are out of range, e.g., > maximum no of electron&
&s of up or down spin. Change the input.')
call errstop_master('QMC_PLOTTER','Quitting')
endif
if(any(xyzzyaaaq1<1).or.any(xyzzyaaaq1>nspin))call errstop_master('QMC&
&_PLOTTER','Spins for orbital plots should be between 1 and '//trim(i2&
&s(nspin))//'.')
allocate(xyzzyaabk6(nemax,real1_complex2,ndet),xyzzyaabl6(3,nemax,real&
&1_complex2,ndet),xyzzyaabm6(nemax,real1_complex2,ndet),xyzzyaabh6(4,n&
&etot),stat=ialloc)
call check_alloc(ialloc,'QMC_PLOTTER','rpsi')
xyzzyaabk6=0.d0
xyzzyaabl6=0.d0
xyzzyaabm6=0.d0
xyzzyaabh6=0.d0
allocate(xyzzyaabo6(wfdet_norb,real1_complex2),xyzzyaabp6(3,wfdet_norb&
&,real1_complex2),xyzzyaabq6(wfdet_norb,real1_complex2),stat=ialloc)
call check_alloc(ialloc,'QMC_PLOTTER','orbval')
xyzzyaabo6=0.d0
xyzzyaabp6=0.d0
xyzzyaabq6=0.d0
if(trim(atom_basis_type)/='gaussian')s_plot=.false.
inquire(file='orbitals.in',exist=xyzzyaabu6)
if(s_plot)then
if(xyzzyaabu6)then
call open_units(io,xyzzyaaap6)
if(xyzzyaaap6/=0)call errstop_master('QMC_PLOTTER','Unable to find fre&
&e I/O unit.')
open(io,file='orbitals.in',status='unknown',iostat=xyzzyaaap6)
if(xyzzyaaap6/=0)call errstop('QMC_PLOTTER','Unable to open orbitals.i&
&n.')
read(io,*,iostat=xyzzyaaap6)xyzzyaaaf6,xyzzyaaag6,xyzzyaaah6
if(((xyzzyaaah6/=1.and.xyzzyaaah6/=2).or.xyzzyaaaf6<1.or.xyzzyaaag6<1.&
&or.xyzzyaaaf6>nele(xyzzyaaah6).or.xyzzyaaag6>nitot.or.xyzzyaaap6/=0).&
&and.am_master)call errwarn('QMC_PLOTTER','Problem with orbitals.in fi&
&le')
if(xyzzyaaaf6>nuc_nele(xyzzyaaah6))then
xyzzyaabw6=.true.
xyzzyaaaf6=xyzzyaaaf6-nuc_nele(xyzzyaaah6)
else
xyzzyaabw6=.false.
endif
close(io)
open_unit(io)=.false.
else
call errstop_master('QMC_PLOTTER','Keyword s_plot set to T, but no orb&
&itals.in file exists.')
endif
endif
if(pairing_wf)then
call points(xyzzyaabi6,xyzzyaaar6,initial_rele,initial_rele_set)
call mpi_bcast(xyzzyaabi6,three_netot,mpi_double_precision,0,mpi_comm_&
&world,ierror)
call checkmpi(ierror,'Broadcasting rele in qmc_plotter.')
call mpi_bcast(xyzzyaaar6,netot,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting sele in qmc_plotter.')
do ifn=1,plot_nfunctions
do xyzzyaaad6=1,nspin
xyzzyaaac6=-heg_orbtype(xyzzyaaad6,1)
if(nuc_nele(xyzzyaaad6)<xyzzyaaap1(ifn).and.xyzzyaaac6>0)then
xyzzyaaaa6=xyzzyaaap1(ifn)-nuc_nele(xyzzyaaad6)
if(xyzzyaaac6>1)xyzzyaaaa6=xyzzyaaaa6+sum(nele(1:xyzzyaaac6-1))
xyzzyaabi6(1:3,xyzzyaaaa6)=xyzzyaaat1(1:3)
endif
enddo
enddo
endif
case(xyzzyaaaf1)
allocate(xyzzyaabj6(4,nitot,netot),xyzzyaabn6(3,nitot),stat=ialloc)
call check_alloc(ialloc,'QMC_PLOTTER','2')
case(xyzzyaaag1,xyzzyaaah1,xyzzyaaai1)
if(nele(1)<1)call errstop_master('QMC_PLOTTER','Need particles of type&
& 1 in order to plot the energy or wave function.')
if(any(xyzzyaaar1>nspin).or.any(xyzzyaaar1<1))call errstop_master('QMC&
&_PLOTTER','Spin of fixed particle is out of range.')
xyzzyaabs6=use_jastrow
if(isperiodic.or..not.model_system)use_jastrow=.false.
xyzzyaabt6=use_backflow
if(isperiodic.or..not.model_system)use_backflow=.false.
call init_wfn
xyzzyaaan6=0
xyzzyaaao6=0
vmc_cfg_by_cfg=.false.
call scratch_request(ratio1_from=xyzzyaaan6,ratio1_to=xyzzyaaao6)
call setup_scratch
call which_scratch(xyzzyaaan6)
call which_scratch(xyzzyaaao6)
call setup_wfn_utils
call points(xyzzyaabi6,xyzzyaaar6,initial_rele,initial_rele_set)
xyzzyaaaq6=1
xyzzyaabf6=sqrt(dtvmc_array)
call equilibration(xyzzyaaan6,xyzzyaaao6,xyzzyaabi6,xyzzyaaar6,dtvmc_a&
&rray,xyzzyaabf6,dtvmc_shift_array,opt_dtvmc_array,nequil,xyzzyaaaq6)
call mpi_bcast(xyzzyaabi6,three_netot,mpi_double_precision,0,mpi_comm_&
&world,ierror)
call checkmpi(ierror,'Broadcasting rele in qmc_plotter.')
call mpi_bcast(xyzzyaaar6,netot,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting sele in qmc_plotter.')
call finish_wfn_utils
call finish_scratch
use_jastrow=xyzzyaabs6
use_backflow=xyzzyaabt6
call init_wfn
do xyzzyaaae6=1,plot_nfix_elec
xyzzyaaaa6=xyzzyaaas1(xyzzyaaae6)
if(xyzzyaaar1(xyzzyaaae6)>1)xyzzyaaaa6=xyzzyaaaa6+sum(nele(1:xyzzyaaar&
&1(xyzzyaaae6)-1))
if(xyzzyaaaa6==1)call errstop_master('QMC_PLOTTER','You have requested&
& to fix particle number 1, but it needs to be moved to do the plot. R&
&emove particle number 1 from the list of fixed particles in the QMC_P&
&LOT block')
if(xyzzyaaaa6>netot)call errstop_master('QMC_PLOTTER','List of fixed p&
&articles in qmc_plot block includes non-existing particles?')
xyzzyaabi6(1:dimensionality,xyzzyaaaa6)=xyzzyaaaw1(1:dimensionality,xy&
&zzyaaae6)
enddo
if(isperiodic)call map_to_simcell(3,netot,xyzzyaabi6)
if(am_master)then
call wout('Particle positions for plot:')
write(tmpr,'(a2,1x,a3,3(1x,a20))')'s','i','x','y','z'
call wout(tmpr)
call wout(repeat('-',70))
do xyzzyaaac6=1,nspin
do xyzzyaaab6=1,nele(xyzzyaaac6)
xyzzyaaaa6=xyzzyaaab6
if(xyzzyaaac6>1)xyzzyaaaa6=xyzzyaaaa6+sum(nele(1:xyzzyaaac6-1))
if(xyzzyaaaa6==1)cycle
write(tmpr,'(i2,1x,i3,3(1x,f20.12))')xyzzyaaac6,xyzzyaaab6,xyzzyaabi6(&
&1:3,xyzzyaaaa6)
call wout(tmpr)
enddo
enddo
call wout()
endif
xyzzyaaan6=0
xyzzyaaao6=0
call scratch_request(ratio1_from=xyzzyaaan6,ratio1_to=xyzzyaaao6)
if(what_to_plot==xyzzyaaag1)call energy_scratch_request(xyzzyaaan6)
call setup_scratch
call which_scratch(xyzzyaaan6)
call which_scratch(xyzzyaaao6)
call setup_wfn_utils
if(what_to_plot==xyzzyaaag1)then
call setup_energy_utils
if(isperiodic.and.model_system)then
xyzzyaabb6=1.d0/dble(netot)
else
xyzzyaabb6=1.d0/dble(npcells)
endif
endif
call define_config(xyzzyaaan6,xyzzyaabi6,xyzzyaaar6)
end select
xyzzyaaai6=0
outer: do xyzzyaaal6=1,grid_npoint_ad
do xyzzyaaak6=1,grid_npoint_ac
do xyzzyaaaj6=1,grid_npoint_ab
xyzzyaaai6=xyzzyaaai6+1
if(xyzzyaaai6<iplotmin)cycle
if(xyzzyaaai6>iplotmax)exit outer
xyzzyaaaz6(1:3)=grid_point(1:3,xyzzyaaai6)
xyzzyaabi6(1:3,1)=xyzzyaaaz6(1:3)
select case(what_to_plot)
case(xyzzyaaaa1,xyzzyaaab1,xyzzyaaac1,xyzzyaaad1,xyzzyaaae1)
if(all(xyzzyaaaz6(:)==0.d0).and.what_to_plot/=xyzzyaaaa1)then
valid_plot_point(xyzzyaaai6)=.false.
cycle
endif
if(pairing_wf)call ee_distances(netot,xyzzyaabi6(1:3,1),xyzzyaabi6,xyz&
&zyaabh6)
if(s_plot)xyzzyaaat6=sqrt((xyzzyaaaz6(1)-rion(1,1))**2+(xyzzyaaaz6(2)-&
&rion(2,1))**2+(xyzzyaaaz6(3)-rion(3,1))**2)
do xyzzyaaad6=1,nspin
if(all(xyzzyaaaq1/=xyzzyaaad6))cycle
if(what_to_plot==xyzzyaaaa1)then
if(.not.s_plot)then
call wfdet(xyzzyaaaz6,0,xyzzyaaad6,wfdet_norb,wfdet_orbmask(1,xyzzyaaa&
&d6),.true.,.false.,xyzzyaabo6,xyzzyaabp6,xyzzyaabq6,xyzzyaabh6)
call copy_orb_to_det(1,xyzzyaaad6,wfdet_orbmap,xyzzyaabo6,xyzzyaabk6)
else
call wfdet_s(xyzzyaaat6,0,1.d0,xyzzyaaaf6,xyzzyaaag6,xyzzyaaah6,1,1,xy&
&zzyaaag6,xyzzyaabw6,which_ssingle(xyzzyaaad6,spin_dep_gs),xyzzyaabc6,&
&xyzzyaabd6,xyzzyaabe6)
endif
else
if(.not.s_plot)then
call wfdet(xyzzyaaaz6,0,xyzzyaaad6,wfdet_norb,wfdet_orbmask(1,xyzzyaaa&
&d6),.false.,.true.,xyzzyaabo6,xyzzyaabp6,xyzzyaabq6,xyzzyaabh6)
call copy_orb_to_det(3,xyzzyaaad6,wfdet_orbmap,xyzzyaabp6,xyzzyaabl6)
call copy_orb_to_det(1,xyzzyaaad6,wfdet_orbmap,xyzzyaabq6,xyzzyaabm6)
else
call wfdet_s(xyzzyaaat6,1,1.d0,xyzzyaaaf6,xyzzyaaag6,xyzzyaaah6,1,1,xy&
&zzyaaag6,xyzzyaabw6,which_ssingle(xyzzyaaad6,spin_dep_gs),xyzzyaabc6,&
&xyzzyaabd6,xyzzyaabe6)
endif
endif
do ifn=1,plot_nfunctions
if(xyzzyaaaq1(ifn)/=xyzzyaaad6)cycle
if(xyzzyaaap1(ifn)>nele(xyzzyaaad6))call errstop_master('QMC_PLOTTER',&
&'Orbital to be plotted out of range for selected spin.')
select case (what_to_plot)
case(xyzzyaaaa1)
if(.not.s_plot)then
plot_data(1,xyzzyaaai6,ifn)=xyzzyaabk6(xyzzyaaap1(ifn),1,1)
else
plot_data(1,xyzzyaaai6,ifn)=xyzzyaabc6(1)
endif
case(xyzzyaaab1)
if(.not.s_plot)then
plot_data(1,xyzzyaaai6,ifn)=xyzzyaabl6(1,xyzzyaaap1(ifn),1,1)
else
plot_data(1,xyzzyaaai6,ifn)=xyzzyaabd6(1)
endif
case(xyzzyaaac1)
if(.not.s_plot)then
plot_data(1,xyzzyaaai6,ifn)=xyzzyaabl6(2,xyzzyaaap1(ifn),1,1)
else
plot_data(1,xyzzyaaai6,ifn)=xyzzyaabd6(1)
endif
case(xyzzyaaad1)
if(.not.s_plot)then
plot_data(1,xyzzyaaai6,ifn)=xyzzyaabl6(3,xyzzyaaap1(ifn),1,1)
else
plot_data(1,xyzzyaaai6,ifn)=xyzzyaabd6(1)
endif
case(xyzzyaaae1)
if(.not.s_plot)then
plot_data(1,xyzzyaaai6,ifn)=xyzzyaabm6(xyzzyaaap1(ifn),1,1)
else
plot_data(1,xyzzyaaai6,ifn)=xyzzyaabe6(1)
endif
end select
enddo
enddo
case(xyzzyaaaf1)
call ei_distances(xyzzyaaaz6,nitot,rion,xyzzyaabj6(1,1,1))
call ionic_potential(1,xyzzyaabj6,xyzzyaaax6,xyzzyaaay6,xyzzyaaba6,xyz&
&zyaabn6,xyzzyaabx6,xyzzyaaby6,isinf)
if(xyzzyaabx6.or.xyzzyaaby6.or.isinf)then
plot_data(1,xyzzyaaai6,1)=0.d0
valid_plot_point(xyzzyaaai6)=.false.
else
plot_data(1,xyzzyaaai6,1)=xyzzyaaax6
endif
case(xyzzyaaag1)
call define_config_oneelec(1,xyzzyaaan6,xyzzyaaao6,xyzzyaaaz6,xyzzyaaa&
&r6(1))
call accept_move(xyzzyaaan6,xyzzyaaao6)
xyzzyaabi6(:,1)=xyzzyaaaz6(:)
call eval_local_energy(xyzzyaaan6,etot=etot,fix_nl_grid=.true.,isnan=i&
&snan,isinf=isinf)
if(isnan.or.isinf)then
plot_data(1,xyzzyaaai6,1)=0.d0
valid_plot_point(xyzzyaaai6)=.false.
else
plot_data(1,xyzzyaaai6,1)=etot*xyzzyaabb6
endif
case(xyzzyaaah1,xyzzyaaai1)
if(all(xyzzyaaaz6(:)==0.d0))then
valid_plot_point(xyzzyaaai6)=.false.
cycle
endif
call define_config_oneelec(1,xyzzyaaan6,xyzzyaaao6,xyzzyaaaz6,xyzzyaaa&
&r6(1))
call wfn_ratio(xyzzyaaan6,xyzzyaaao6,0,ratio=xyzzyaabr6,isnan=isnan,is&
&inf=isinf)
if(isnan.or.isinf)then
plot_data(1,xyzzyaaai6,1)=0.d0
valid_plot_point(xyzzyaaai6)=.false.
else
plot_data(1,xyzzyaaai6,1)=dble(xyzzyaabr6)
endif
case(xyzzyaaaj1)
call eval_expot(xyzzyaaaz6,which_spin(1),xyzzyaaas6,isinf)
if(isinf)then
plot_data(1,xyzzyaaai6,1)=0.d0
valid_plot_point(xyzzyaaai6)=.false.
else
plot_data(1,xyzzyaaai6,1)=xyzzyaaas6
endif
case default
call errstop_master('QMC_PLOTTER','Don''t know what I''m plotting. Bug&
&.')
end select
enddo
enddo
enddo outer
if(am_master)then
do his_node=1,nnodes-1
iplotmin=nint(dble(tot_npoint)/dble(nnodes)*dble(his_node))+1
if(his_node==nnodes-1)then
iplotmax=tot_npoint
else
iplotmax=nint(dble(tot_npoint)/dble(nnodes)*dble(his_node+1))
endif
call mpi_recv(valid_plot_point(iplotmin:iplotmax),iplotmax-iplotmin+1,&
&mpi_logical,his_node,1,mpi_comm_world,status,ierror)
call checkmpi(ierror,'recv valid_plot_point in qmc_plotter.')
call mpi_recv(plot_data(1:plot_ncomponents,iplotmin:iplotmax,1:plot_nf&
&unctions),plot_ncomponents*(iplotmax-iplotmin+1)*plot_nfunctions,mpi_&
&double_precision,his_node,2,mpi_comm_world,status,ierror)
call checkmpi(ierror,'recv plot_data in qmc_plotter.')
enddo
else
call mpi_ssend(valid_plot_point(iplotmin:iplotmax),iplotmax-iplotmin+1&
&,mpi_logical,0,1,mpi_comm_world,ierror)
call checkmpi(ierror,'ssend valid_plot_point in qmc_plotter.')
call mpi_ssend(plot_data(1:plot_ncomponents,iplotmin:iplotmax,1:plot_n&
&functions),plot_ncomponents*(iplotmax-iplotmin+1)*plot_nfunctions,mpi&
&_double_precision,0,2,mpi_comm_world,ierror)
endif
if(am_master)then
if(all(.not.valid_plot_point))then
call errwarn('QMC_PLOTTER','Apparently none of the plot points are val&
&id.  No output produced.  Please check that both the path you plot an&
&d the positions of any fixed particles are sensible.')
else
call open_units(io,xyzzyaaap6)
if(xyzzyaaap6/=0)call errstop_master('QMC_PLOTTER','Unable to find fre&
&e I/O unit for plot file.')
open(io,file=filename,status='unknown',form='formatted',iostat=xyzzyaa&
&ap6)
if(xyzzyaaap6/=0)call errstop('QMC_PLOTTER','Unable to open '//trim(fi&
&lename)//'.')
if(plot_dim==1)then
xyzzyaaaz6(:)=xyzzyaaat1(:)-xyzzyaaau1(:)
xyzzyaabg6=sqrt(sum((xyzzyaaat1(:)-xyzzyaaau1(:))**2))/real(grid_npoin&
&t_ab-1,dp)
if(xyzzyaaaz6(2)==0.d0.and.xyzzyaaaz6(3)==0.d0)then
xyzzyaaaz6(1)=xyzzyaaat1(1)
elseif(xyzzyaaaz6(1)==0.d0.and.xyzzyaaaz6(3)==0.d0)then
xyzzyaaaz6(1)=xyzzyaaat1(2)
elseif(xyzzyaaaz6(1)==0.d0.and.xyzzyaaaz6(2)==0.d0)then
xyzzyaaaz6(1)=xyzzyaaat1(3)
else
xyzzyaaaz6(1)=0.d0
endif
endif
do ifn=1,plot_nfunctions
if(ifn/=1)write(io,*)'&'
xyzzyaaai6=0
do xyzzyaaal6=1,grid_npoint_ad
if(xyzzyaaal6>1)write(io,*)
do xyzzyaaak6=1,grid_npoint_ac
if(xyzzyaaak6>1)write(io,*)
do xyzzyaaaj6=1,grid_npoint_ab
xyzzyaaai6=xyzzyaaai6+1
if(.not.valid_plot_point(xyzzyaaai6).and.what_to_plot/=xyzzyaaaa1)cycl&
&e
if(plot_dim==1)then
write(io,*)xyzzyaaaz6(1)+real(xyzzyaaaj6-1,dp)*xyzzyaabg6,plot_data(1:&
&plot_ncomponents,xyzzyaaai6,ifn)
else
write(io,'(e18.10,3(1x,e18.10))')grid_point(1:3,xyzzyaaai6),plot_data(&
&1:plot_ncomponents,xyzzyaaai6,ifn)
endif
enddo
enddo
enddo
enddo
close(io)
open_unit(io)=.false.
if(plot_dim==1)then
call wout('Plot finished. Use e.g. xmgrace or gnuplot to view.')
elseif(plot_dim==2)then
call wout('Plot finished. Use e.g. plot_2D utility to view.')
else
call wout('Plot finished.')
endif
endif
endif
if(what_to_plot==xyzzyaaai1.and.am_master)then
call open_units(io,xyzzyaaap6)
if(xyzzyaaap6/=0)call errstop_master('QMC_PLOTTER','Unable to find fre&
&e I/O unit for node file.')
filename=trim(xyzzyaaam1(plot_dim))
open(io,file=filename,status='unknown',form='formatted',iostat=xyzzyaa&
&ap6)
if(xyzzyaaap6/=0)call errstop('QMC_PLOTTER','Unable to open '//trim(fi&
&lename)//'.')
xyzzyaaam6=0
select case(plot_dim)
case(1)
do xyzzyaaaj6=1,grid_npoint_ab-1
xyzzyaabv6=.false.
xyzzyaaam6(1,1,1)=nint(sign(1.d0,plot_data(1,map_point(xyzzyaaaj6,1,1)&
&,1)))
xyzzyaaam6(2,1,1)=nint(sign(1.d0,plot_data(1,map_point(xyzzyaaaj6+1,1,&
&1),1)))
xyzzyaabv6=(xyzzyaaam6(2,1,1)/=xyzzyaaam6(1,1,1))
if(xyzzyaabv6)then
xyzzyaaaz6=0.5d0*(grid_point(1:3,map_point(xyzzyaaaj6,1,1))+grid_point&
&(1:3,map_point(xyzzyaaaj6+1,1,1)))
write(io,*)xyzzyaaaz6
endif
enddo
call wout('Nodal points plotted to '//trim(filename)//'.')
case(2)
do xyzzyaaak6=1,grid_npoint_ac-1
do xyzzyaaaj6=1,grid_npoint_ab-1
xyzzyaabv6=.false.
xyzzyaaam6(1,1,1)=nint(sign(1.d0,plot_data(1,map_point(xyzzyaaaj6,xyzz&
&yaaak6,1),1)))
xyzzyaaam6(2,1,1)=nint(sign(1.d0,plot_data(1,map_point(xyzzyaaaj6+1,xy&
&zzyaaak6,1),1)))
xyzzyaaam6(1,2,1)=nint(sign(1.d0,plot_data(1,map_point(xyzzyaaaj6,xyzz&
&yaaak6+1,1),1)))
xyzzyaaam6(2,2,1)=nint(sign(1.d0,plot_data(1,map_point(xyzzyaaaj6+1,xy&
&zzyaaak6+1,1),1)))
xyzzyaabv6=any(xyzzyaaam6(:,:,1)/=xyzzyaaam6(1,1,1))
if(xyzzyaabv6)then
xyzzyaaaz6=0.5d0*(grid_point(1:3,map_point(xyzzyaaaj6,xyzzyaaak6,1))+g&
&rid_point(1:3,map_point(xyzzyaaaj6+1,xyzzyaaak6+1,1)))
write(io,*)xyzzyaaaz6
endif
enddo
enddo
call wout('Nodal curve plotted to '//trim(filename)//'.')
case(3)
do xyzzyaaal6=1,grid_npoint_ad-1
do xyzzyaaak6=1,grid_npoint_ac-1
do xyzzyaaaj6=1,grid_npoint_ab-1
xyzzyaabv6=.false.
xyzzyaaam6(1,1,1)=nint(sign(1.d0,plot_data(1,map_point(xyzzyaaaj6,xyzz&
&yaaak6,xyzzyaaal6),1)))
xyzzyaaam6(2,1,1)=nint(sign(1.d0,plot_data(1,map_point(xyzzyaaaj6+1,xy&
&zzyaaak6,xyzzyaaal6),1)))
xyzzyaaam6(1,2,1)=nint(sign(1.d0,plot_data(1,map_point(xyzzyaaaj6,xyzz&
&yaaak6+1,xyzzyaaal6),1)))
xyzzyaaam6(1,1,2)=nint(sign(1.d0,plot_data(1,map_point(xyzzyaaaj6,xyzz&
&yaaak6,xyzzyaaal6+1),1)))
xyzzyaaam6(2,2,1)=nint(sign(1.d0,plot_data(1,map_point(xyzzyaaaj6+1,xy&
&zzyaaak6+1,xyzzyaaal6),1)))
xyzzyaaam6(2,1,2)=nint(sign(1.d0,plot_data(1,map_point(xyzzyaaaj6+1,xy&
&zzyaaak6,xyzzyaaal6+1),1)))
xyzzyaaam6(1,2,2)=nint(sign(1.d0,plot_data(1,map_point(xyzzyaaaj6,xyzz&
&yaaak6+1,xyzzyaaal6+1),1)))
xyzzyaaam6(2,2,2)=nint(sign(1.d0,plot_data(1,map_point(xyzzyaaaj6+1,xy&
&zzyaaak6+1,xyzzyaaal6+1),1)))
xyzzyaabv6=any(xyzzyaaam6(:,:,:)/=xyzzyaaam6(1,1,1))
if(xyzzyaabv6)then
xyzzyaaaz6=0.5d0*(grid_point(1:3,map_point(xyzzyaaaj6,xyzzyaaak6,xyzzy&
&aaal6))+grid_point(1:3,map_point(xyzzyaaaj6+1,xyzzyaaak6+1,xyzzyaaal6&
&+1)))
write(io,*)xyzzyaaaz6
endif
enddo
enddo
enddo
call wout('Nodal surface plotted to '//trim(filename)//'.')
end select
close(io)
open_unit(io)=.false.
endif
deallocate(plot_data,grid_point,map_point,valid_plot_point)
deallocate(xyzzyaabi6,xyzzyaaar6)
select case(what_to_plot)
case(xyzzyaaaa1,xyzzyaaab1,xyzzyaaac1,xyzzyaaad1,xyzzyaaae1)
deallocate(xyzzyaabk6,xyzzyaabl6,xyzzyaabm6,xyzzyaabh6)
case(xyzzyaaaf1)
deallocate(xyzzyaabj6,xyzzyaabn6)
case(xyzzyaaag1,xyzzyaaah1,xyzzyaaai1)
if(what_to_plot==xyzzyaaag1)call finish_energy_utils
call finish_wfn_utils
call finish_scratch
end select
call timer('QMC_PLOTTER',.false.)
end subroutine qmc_plotter
end module slaarnabyter

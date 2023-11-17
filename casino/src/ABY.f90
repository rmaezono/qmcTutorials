module slaarnaby
use dsp
use slaarnaam
use parallel
use slaarnach
use slaarnacs
use slaarnaan, only : points,map_to_simcell
use file_utils,    only : open_units
use format_utils,  only : wout,i2s,l2s,capitalize,get_field,pad_int,wo&
&rdwrap
use slaarnabg,      only : dimensionality,isperiodic
use run_control,   only : errstop_master,timer,check_alloc
use store,         only : netot,nspin,nele,three_netot,open_unit,vmc_c&
&fg_by_cfg,no_difftypes,which_ii,which_ie,which_spin
use slaarnacp,           only : equilibration
implicit none
private
public init_plot,setup_plot,plot_main
integer :: xyzzyaaaa1=0
character(180),allocatable :: xyzzyaaab1(:)
integer what_to_plot,what_to_plot_sec,plot_dim,grid_npoint_ab,grid_npo&
&int_ac,grid_npoint_ad,plot_nfunctions,plot_ncomponents,plot_ii,plot_r&
&ank
real(dp) xyzzyaaac1(3),xyzzyaaad1(3),xyzzyaaae1(3),xyzzyaaaf1(3)
real(dp),allocatable :: xyzzyaaag1(:,:)
logical xyzzyaaah1,xyzzyaaai1,xyzzyaaaj1,xyzzyaaak1
logical,allocatable :: xyzzyaaal1(:)
character(64) now_plotting_name
character(64),allocatable :: function_name(:)
contains
subroutine init_plot(nlines,input_block)
implicit none
integer,intent(in) :: nlines
character(180),intent(in) :: input_block(nlines)
integer xyzzyaaaa2,xyzzyaaab2
if(xyzzyaaaa1>0.or.allocated(xyzzyaaab1))call errstop_master('INIT_PLO&
&T','Reading PLOT input block twice.')
xyzzyaaaa1=nlines
if(nlines<1)return
allocate(xyzzyaaab1(xyzzyaaaa1),stat=xyzzyaaaa2)
call check_alloc(xyzzyaaaa2,'INIT_PLOT','plot_block')
do xyzzyaaab2=1,nlines
xyzzyaaab1(xyzzyaaab2)=trim(adjustl(input_block(xyzzyaaab2)))
enddo
end subroutine init_plot
subroutine setup_plot
implicit none
integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3,ialloc,xyzzyaaad3,xyzzyaaae3(&
&2),xyzzyaaaf3,xyzzyaaag3,xyzzyaaah3,xyzzyaaai3,xyzzyaaaj3
integer,allocatable :: xyzzyaaak3(:)
real(dp) xyzzyaaal3,xyzzyaaam3(3),xyzzyaaan3(3),xyzzyaaao3(3),xyzzyaaa&
&p3(3)
real(dp),parameter :: xyzzyaaaq3=1.d-13
character(180) line,field
character(64),allocatable :: xyzzyaaar3(:),xyzzyaaas3(:)
what_to_plot=0
what_to_plot_sec=0
plot_ii=0
plot_rank=0
xyzzyaaah1=.false.
xyzzyaaai1=.false.
plot_ncomponents=0
plot_nfunctions=0
grid_npoint_ab=0
grid_npoint_ac=0
grid_npoint_ad=0
xyzzyaaac1=0.d0
xyzzyaaad1=0.d0
xyzzyaaae1=0.d0
xyzzyaaaf1=0.d0
xyzzyaaae3=0
call enumerate_plot_wfn(xyzzyaaae3(1))
call enumerate_plot_energy(xyzzyaaae3(2))
xyzzyaaag3=sum(xyzzyaaae3(1:2))
allocate(xyzzyaaak3(xyzzyaaag3),stat=ialloc)
call check_alloc(ialloc,'SETUP_PLOT','which_plot_sec')
do xyzzyaaah3=1,2
xyzzyaaak3(sum(xyzzyaaae3(1:xyzzyaaah3-1))+1:sum(xyzzyaaae3(1:xyzzyaaa&
&h3)))=xyzzyaaah3
enddo
allocate(xyzzyaaar3(xyzzyaaag3),xyzzyaaas3(xyzzyaaag3),stat=ialloc)
call check_alloc(ialloc,'SETUP_PLOT','plot_keyword,...')
call enumerate_plot_wfn(xyzzyaaae3(1),xyzzyaaar3(1:xyzzyaaae3(1)),xyzz&
&yaaas3(1:xyzzyaaae3(1)))
call enumerate_plot_energy(xyzzyaaae3(2),xyzzyaaar3(xyzzyaaae3(1)+1:su&
&m(xyzzyaaae3(1:2))),xyzzyaaas3(xyzzyaaae3(1)+1:sum(xyzzyaaae3(1:2))))
xyzzyaaac3=0
call xyzzyaaat3
read(xyzzyaaab1(xyzzyaaac3),'(a)',iostat=xyzzyaaad3)line
call xyzzyaaau3
line=trim(adjustl(line))
call capitalize(line,decapitalize_in=.true.)
do xyzzyaaaa3=1,xyzzyaaag3
if(trim(line)==trim(xyzzyaaar3(xyzzyaaaa3)))then
what_to_plot=xyzzyaaaa3
exit
endif
enddo
if(what_to_plot==0)then
if(am_master)then
call wout('For this system the following objects can be plotted (keywo&
&rd - description):')
do xyzzyaaaa3=1,xyzzyaaag3
call wout(' * '//trim(xyzzyaaar3(xyzzyaaaa3))//' - '//trim(xyzzyaaas3(&
&xyzzyaaaa3)))
enddo
if(trim(line)/='help')then
call wout()
call wout('The provided keyword "'//trim(line)//'" is not valid.')
endif
endif
call errstop_master('SETUP_PLOT','Quitting.')
endif
now_plotting_name=trim(xyzzyaaas3(xyzzyaaaa3))
call xyzzyaaat3
read(xyzzyaaab1(xyzzyaaac3),'(a)',iostat=xyzzyaaad3)line
call xyzzyaaau3
if(.not.((xyzzyaaav3(1,line,'electron').and.xyzzyaaav3(3,line,'spin'))&
&.or.(xyzzyaaav3(1,line,'particle').and.xyzzyaaav3(3,line,'type'))))ca&
&ll errstop_master('SETUP_PLOT','In line '//trim(i2s(xyzzyaaac3))//' o&
&f PLOT input block: expected "electron <ie> spin <ispin>" or "particl&
&e <ie> type <ispin>".')
field=get_field(2,line)
read(field,*,iostat=xyzzyaaad3)xyzzyaaai3
call xyzzyaaau3
field=get_field(4,line)
read(field,*,iostat=xyzzyaaad3)xyzzyaaaj3
call xyzzyaaau3
if(xyzzyaaaj3<1.or.xyzzyaaaj3>nspin)call errstop_master('SETUP_PLOT','&
&In line '//trim(i2s(xyzzyaaac3))//' of PLOT input block: <ispin> out &
&of range.')
if(xyzzyaaai3<1.or.xyzzyaaai3>nele(xyzzyaaaj3))call errstop_master('SE&
&TUP_PLOT','In line '//trim(i2s(xyzzyaaac3))//' of PLOT input block: <&
&ie> out of range.')
plot_ii=which_ii(xyzzyaaai3,xyzzyaaaj3)
what_to_plot_sec=xyzzyaaak3(what_to_plot)
what_to_plot=what_to_plot-sum(xyzzyaaae3(1:what_to_plot_sec-1))
select case(what_to_plot_sec)
case(1)
call query_plot_wfn(what_to_plot,plot_ii,plot_rank,xyzzyaaah1,xyzzyaaa&
&i1,xyzzyaaaj1,xyzzyaaak1,plot_nfunctions)
case(2)
call query_plot_energy(what_to_plot,plot_rank,xyzzyaaah1,xyzzyaaai1,xy&
&zzyaaaj1,xyzzyaaak1,plot_nfunctions)
end select
if(plot_nfunctions<1)call errstop_master('SETUP_PLOT','No functions to&
& plot according to routine. Possible bug.')
allocate(function_name(plot_nfunctions),stat=ialloc)
call check_alloc(ialloc,'SETUP_PLOT','function_name,...')
select case(what_to_plot_sec)
case(1)
call query_plot_wfn(what_to_plot,plot_ii,plot_rank,xyzzyaaah1,xyzzyaaa&
&i1,xyzzyaaaj1,xyzzyaaak1,plot_nfunctions,function_name)
case(2)
call query_plot_energy(what_to_plot,plot_rank,xyzzyaaah1,xyzzyaaai1,xy&
&zzyaaaj1,xyzzyaaak1,plot_nfunctions,function_name)
end select
plot_ncomponents=dimensionality**plot_rank
if(xyzzyaaah1)plot_ncomponents=plot_ncomponents*2
if(xyzzyaaai1)plot_ncomponents=plot_ncomponents*2
if(plot_ncomponents<1)call errstop_master('SETUP_PLOT','No components &
&to plot according to routine. Bug.')
call xyzzyaaat3
read(xyzzyaaab1(xyzzyaaac3),'(a)',iostat=xyzzyaaad3)line
call xyzzyaaau3
call capitalize(line,decapitalize_in=.true.)
select case(trim(adjustl(line)))
case('line','1','1d')
plot_dim=1
case('plane','2','2d')
plot_dim=2
case('volume','3','3d')
plot_dim=3
case default
if(am_master)then
call wout('Line '//trim(i2s(xyzzyaaac3))//' in the PLOT input block mu&
&st be one of:')
call wout(' - "line", "1D", "1"')
call wout(' - "plane", "2D", "2"')
call wout(' - "volume", "3D", "3"')
endif
call errstop_master('SETUP_PLOT','Stopping.')
end select
if(plot_dim>dimensionality)call errstop_master('SETUP_PLOT','Plot dime&
&nsionality greater than system dimensionality.')
grid_npoint_ab=0
grid_npoint_ac=0
grid_npoint_ad=0
call xyzzyaaat3
read(xyzzyaaab1(xyzzyaaac3),'(a)',iostat=xyzzyaaad3)line
call xyzzyaaau3
if(.not.xyzzyaaav3(1,line,'grid'))call errstop_master('SETUP_PLOT','In&
& line '//trim(i2s(xyzzyaaac3))//' of PLOT input block: expected "grid&
& <n_AB> [<n_AC> [<n_AD>]]".')
field=get_field(2,line)
read(field,*,iostat=xyzzyaaad3)grid_npoint_ab
call xyzzyaaau3
if(grid_npoint_ab<1)call errstop_master('SETUP_PLOT','In line '//trim(&
&i2s(xyzzyaaac3))//' of PLOT input block: <n_AB> out of range.')
if(plot_dim>1)then
field=get_field(3,line)
read(field,*,iostat=xyzzyaaad3)grid_npoint_ac
call xyzzyaaau3
if(grid_npoint_ac<1)call errstop_master('SETUP_PLOT','In line '//trim(&
&i2s(xyzzyaaac3))//' of PLOT input block: <n_AC> out of range.')
endif
if(plot_dim>2)then
field=get_field(4,line)
read(field,*,iostat=xyzzyaaad3)grid_npoint_ad
call xyzzyaaau3
if(grid_npoint_ad<1)call errstop_master('SETUP_PLOT','In line '//trim(&
&i2s(xyzzyaaac3))//' of PLOT input block: <n_AD> out of range.')
endif
xyzzyaaac1=0.d0
xyzzyaaad1=0.d0
xyzzyaaae1=0.d0
xyzzyaaaf1=0.d0
call xyzzyaaat3
read(xyzzyaaab1(xyzzyaaac3),'(a)',iostat=xyzzyaaad3)line
call xyzzyaaau3
if(.not.xyzzyaaav3(1,line,'a'))call errstop_master('SETUP_PLOT','In li&
&ne '//trim(i2s(xyzzyaaac3))//' of PLOT input block: expected "A <x> [&
&<y> [<z>]]".')
do xyzzyaaaf3=1,dimensionality
field=get_field(1+xyzzyaaaf3,line)
read(field,*,iostat=xyzzyaaad3)xyzzyaaac1(xyzzyaaaf3)
call xyzzyaaau3
enddo
call xyzzyaaat3
read(xyzzyaaab1(xyzzyaaac3),'(a)',iostat=xyzzyaaad3)line
call xyzzyaaau3
if(.not.xyzzyaaav3(1,line,'b'))call errstop_master('SETUP_PLOT','In li&
&ne '//trim(i2s(xyzzyaaac3))//' of PLOT input block: expected "B <x> [&
&<y> [<z>]]".')
do xyzzyaaaf3=1,dimensionality
field=get_field(1+xyzzyaaaf3,line)
read(field,*,iostat=xyzzyaaad3)xyzzyaaad1(xyzzyaaaf3)
call xyzzyaaau3
enddo
if(plot_dim>1)then
call xyzzyaaat3
read(xyzzyaaab1(xyzzyaaac3),'(a)',iostat=xyzzyaaad3)line
call xyzzyaaau3
if(.not.xyzzyaaav3(1,line,'c'))call errstop_master('SETUP_PLOT','In li&
&ne '//trim(i2s(xyzzyaaac3))//' of PLOT input block: expected "C <x> [&
&<y> [<z>]]".')
do xyzzyaaaf3=1,dimensionality
field=get_field(1+xyzzyaaaf3,line)
read(field,*,iostat=xyzzyaaad3)xyzzyaaae1(xyzzyaaaf3)
call xyzzyaaau3
enddo
endif
if(plot_dim>2)then
call xyzzyaaat3
read(xyzzyaaab1(xyzzyaaac3),'(a)',iostat=xyzzyaaad3)line
call xyzzyaaau3
if(.not.xyzzyaaav3(1,line,'d'))call errstop_master('SETUP_PLOT','In li&
&ne '//trim(i2s(xyzzyaaac3))//' of PLOT input block: expected "D <x> [&
&<y> [<z>]]".')
do xyzzyaaaf3=1,dimensionality
field=get_field(1+xyzzyaaaf3,line)
read(field,*,iostat=xyzzyaaad3)xyzzyaaaf1(xyzzyaaaf3)
call xyzzyaaau3
enddo
endif
xyzzyaaam3=xyzzyaaac1-xyzzyaaad1
xyzzyaaal3=sqrt(sum(xyzzyaaam3(:)**2))
if(xyzzyaaal3==0.d0)call errstop_master('SETUP_PLOT','Zero line length&
& AB in PLOT input block.')
if(plot_dim>1)then
xyzzyaaan3=xyzzyaaac1-xyzzyaaae1
xyzzyaaal3=sqrt(sum(xyzzyaaan3(:)**2))
if(xyzzyaaal3==0.d0)call errstop_master('SETUP_PLOT','Zero line length&
& AC in PLOT input block.')
xyzzyaaap3(1)=xyzzyaaam3(2)*xyzzyaaan3(3)-xyzzyaaam3(3)*xyzzyaaan3(2)
xyzzyaaap3(2)=xyzzyaaam3(3)*xyzzyaaan3(1)-xyzzyaaam3(1)*xyzzyaaan3(3)
xyzzyaaap3(3)=xyzzyaaam3(1)*xyzzyaaan3(2)-xyzzyaaam3(2)*xyzzyaaan3(1)
if(all(abs(xyzzyaaap3)<xyzzyaaaq3))call errstop_master('SETUP_PLOT','P&
&oints A, B and C in PLOT input block are collinear.')
endif
if(plot_dim>2)then
xyzzyaaao3=xyzzyaaac1-xyzzyaaaf1
xyzzyaaal3=sqrt(sum(xyzzyaaao3(:)**2))
if(xyzzyaaal3==0.d0)call errstop_master('SETUP_PLOT','Zero line length&
& AD for plot.')
xyzzyaaap3(1)=xyzzyaaam3(2)*xyzzyaaao3(3)-xyzzyaaam3(3)*xyzzyaaao3(2)
xyzzyaaap3(2)=xyzzyaaam3(3)*xyzzyaaao3(1)-xyzzyaaam3(1)*xyzzyaaao3(3)
xyzzyaaap3(3)=xyzzyaaam3(1)*xyzzyaaao3(2)-xyzzyaaam3(2)*xyzzyaaao3(1)
if(all(abs(xyzzyaaap3)<xyzzyaaaq3))call errstop_master('SETUP_PLOT','P&
&oints A, B and D in PLOT input block are collinear.')
xyzzyaaap3(1)=xyzzyaaan3(2)*xyzzyaaao3(3)-xyzzyaaan3(3)*xyzzyaaao3(2)
xyzzyaaap3(2)=xyzzyaaan3(3)*xyzzyaaao3(1)-xyzzyaaan3(1)*xyzzyaaao3(3)
xyzzyaaap3(3)=xyzzyaaan3(1)*xyzzyaaao3(2)-xyzzyaaan3(2)*xyzzyaaao3(1)
if(all(abs(xyzzyaaap3)<xyzzyaaaq3))call errstop_master('SETUP_PLOT','P&
&oints A, C and D in PLOT input block are collinear.')
endif
allocate(xyzzyaaag1(3,netot),xyzzyaaal1(netot),stat=ialloc)
call check_alloc(ialloc,'SETUP_PLOT','plot_fix_rele')
xyzzyaaag1=0.d0
xyzzyaaal1=.false.
do while(xyzzyaaaa1>xyzzyaaac3)
call xyzzyaaat3
read(xyzzyaaab1(xyzzyaaac3),'(a)',iostat=xyzzyaaad3)line
call xyzzyaaau3
field=get_field(1,line)
select case(trim(field))
case('fix')
if(.not.(((xyzzyaaav3(2,line,'electron').and.xyzzyaaav3(4,line,'spin')&
&).or.(xyzzyaaav3(2,line,'particle').and.xyzzyaaav3(4,line,'type'))).a&
&nd.xyzzyaaav3(6,line,'@')))call errstop_master('SETUP_PLOT','In line &
&'//trim(i2s(xyzzyaaac3))//' of PLOT input block: expected "fix electr&
&on <ie> spin <ispin> @ <x> [<y> [<z>]]" or "fix particle <ie> type <i&
&spin> @ <x> [<y> [<z>]]".')
field=get_field(3,line)
read(field,*,iostat=xyzzyaaad3)xyzzyaaai3
call xyzzyaaau3
field=get_field(5,line)
read(field,*,iostat=xyzzyaaad3)xyzzyaaaj3
call xyzzyaaau3
if(xyzzyaaaj3<1.or.xyzzyaaaj3>nspin)call errstop_master('SETUP_PLOT','&
&In line '//trim(i2s(xyzzyaaac3))//' of PLOT input block: <ispin> out &
&of range.')
if(xyzzyaaai3<1.or.xyzzyaaai3>nele(xyzzyaaaj3))call errstop_master('SE&
&TUP_PLOT','In line '//trim(i2s(xyzzyaaac3))//' of PLOT input block: <&
&ie> out of range.')
xyzzyaaab3=which_ii(xyzzyaaai3,xyzzyaaaj3)
if(xyzzyaaal1(xyzzyaaab3))call errstop_master('SETUP_PLOT','In line '/&
&/trim(i2s(xyzzyaaac3))//' of PLOT input block: this particle is alrea&
&dy fixed.')
xyzzyaaal1(xyzzyaaab3)=.true.
do xyzzyaaaf3=1,dimensionality
field=get_field(6+xyzzyaaaf3,line)
read(field,*,iostat=xyzzyaaad3)xyzzyaaag1(xyzzyaaaf3,xyzzyaaab3)
call xyzzyaaau3
enddo
case default
call errstop_master('SETUP_PLOT','Extra line '//trim(i2s(xyzzyaaac3))/&
&/' in PLOT input block not understood.')
end select
enddo
deallocate(xyzzyaaas3,xyzzyaaar3)
deallocate(xyzzyaaab1)
xyzzyaaaa1=0
contains
subroutine xyzzyaaat3()
implicit none
xyzzyaaac3=xyzzyaaac3+1
if(xyzzyaaaa1<xyzzyaaac3)call errstop_master('SETUP_PLOT','Expecting t&
&o find at least '//trim(i2s(xyzzyaaac3))//' lines in the PLOT input b&
&lock.')
end subroutine xyzzyaaat3
subroutine xyzzyaaau3()
implicit none
if(xyzzyaaad3>0)then
call errstop_master('SETUP_PLOT','Problem reading line '//trim(i2s(xyz&
&zyaaac3))//' in the PLOT input block.')
elseif(xyzzyaaad3<0)then
call errstop_master('SETUP_PLOT','Reached end of line '//trim(i2s(xyzz&
&yaaac3))//' in the PLOT input block without finding all the data it s&
&hould contain.')
endif
end subroutine xyzzyaaau3
logical function xyzzyaaav3(ifield,line,string)
implicit none
integer,intent(in) :: ifield
character(*),intent(in) :: line,string
character(180) tstring
tstring=line
call capitalize(tstring,decapitalize_in=.true.)
xyzzyaaav3=trim(get_field(ifield,tstring))==trim(string)
end function xyzzyaaav3
end subroutine setup_plot
subroutine plot_main(nequil,dtvmc_array,dtvmc_shift_array,opt_dtvmc_ar&
&ray,initial_rele,initial_rele_set)
implicit none
integer,intent(in) :: nequil,opt_dtvmc_array(no_difftypes)
real(dp),intent(in) :: initial_rele(3,netot)
real(dp),intent(inout) :: dtvmc_array(no_difftypes),dtvmc_shift_array(&
&no_difftypes)
logical,intent(in) :: initial_rele_set(netot)
integer xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7,xyzzyaaad7,tot_npoint,xyzzyaa&
&ae7,xyzzyaaaf7,xyzzyaaag7,xyzzyaaah7,ifn,xyzzyaaai7,xyzzyaaaj7,iplotm&
&in,iplotmax,inode,xyzzyaaak7,ialloc,xyzzyaaal7
integer,allocatable :: map_point(:,:,:),xyzzyaaam7(:)
real(dp) xyzzyaaan7(3),xyzzyaaao7(3),xyzzyaaap7(3),xyzzyaaaq7(3),xyzzy&
&aaar7(no_difftypes)
real(dp),allocatable :: grid_point(:,:),xyzzyaaas7(:,:),plot_data(:,:,&
&:)
character(64) filename
character(80) tmpr
character(512) line
call timer('PLOT_MAIN',.true.)
xyzzyaaan7(1:3)=(xyzzyaaad1(1:3)-xyzzyaaac1(1:3))/dble(grid_npoint_ab-&
&1)
if(plot_dim>1)then
xyzzyaaao7(1:3)=(xyzzyaaae1(1:3)-xyzzyaaac1(1:3))/dble(grid_npoint_ac-&
&1)
else
xyzzyaaao7=0.d0
grid_npoint_ac=1
endif
if(plot_dim>2)then
xyzzyaaap7(1:3)=(xyzzyaaaf1(1:3)-xyzzyaaac1(1:3))/dble(grid_npoint_ad-&
&1)
else
xyzzyaaap7=0.d0
grid_npoint_ad=1
endif
tot_npoint=grid_npoint_ab*grid_npoint_ac*grid_npoint_ad
allocate(plot_data(plot_ncomponents,plot_nfunctions,tot_npoint),grid_p&
&oint(3,tot_npoint),map_point(grid_npoint_ab,grid_npoint_ac,grid_npoin&
&t_ad),stat=ialloc)
call check_alloc(ialloc,'PLOT_MAIN','plot_data')
xyzzyaaae7=0
do xyzzyaaah7=1,grid_npoint_ad
do xyzzyaaag7=1,grid_npoint_ac
do xyzzyaaaf7=1,grid_npoint_ab
xyzzyaaae7=xyzzyaaae7+1
map_point(xyzzyaaaf7,xyzzyaaag7,xyzzyaaah7)=xyzzyaaae7
grid_point(1:3,xyzzyaaae7)=xyzzyaaac1(1:3)+real(xyzzyaaaf7-1,dp)*xyzzy&
&aaan7(1:3)+real(xyzzyaaag7-1,dp)*xyzzyaaao7(1:3)+real(xyzzyaaah7-1,dp&
&)*xyzzyaaap7(1:3)
enddo
enddo
enddo
if(tot_npoint<nnodes)call errstop_master('PLOT_MAIN','Number of points&
& to plot is less than the number of processors.')
iplotmin=nint(dble(tot_npoint)/dble(nnodes)*dble(my_node))+1
if(my_node==nnodes-1)then
iplotmax=tot_npoint
else
iplotmax=nint(dble(tot_npoint)/dble(nnodes)*dble(my_node+1))
endif
if(am_master)then
call wout(trim(i2s(plot_dim))//'-DIMENSIONAL PLOT REQUESTED')
call wout('============================')
call wout()
call wout('Preparing system for plot.')
call wout()
endif
allocate(xyzzyaaas7(3,netot),xyzzyaaam7(netot),stat=ialloc)
call check_alloc(ialloc,'PLOT_MAIN','rele')
xyzzyaaas7=0.d0
xyzzyaaam7=0
if(count(xyzzyaaal1)<netot)then
xyzzyaaai7=0
xyzzyaaaj7=0
vmc_cfg_by_cfg=.false.
call scratch_request(ratio1_from=xyzzyaaai7,ratio1_to=xyzzyaaaj7)
call setup_scratch
call which_scratch(xyzzyaaai7)
call which_scratch(xyzzyaaaj7)
call setup_wfn_utils
call points(xyzzyaaas7,xyzzyaaam7,initial_rele,initial_rele_set)
xyzzyaaal7=1
xyzzyaaar7=sqrt(dtvmc_array)
call equilibration(xyzzyaaai7,xyzzyaaaj7,xyzzyaaas7,xyzzyaaam7,dtvmc_a&
&rray,xyzzyaaar7,dtvmc_shift_array,opt_dtvmc_array,nequil,xyzzyaaal7)
call mpi_bcast(xyzzyaaas7,three_netot,mpi_double_precision,0,mpi_comm_&
&world,ierror)
call checkmpi(ierror,'Broadcasting rele in plot_main.')
call mpi_bcast(xyzzyaaam7,netot,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting sele in plot_main.')
call finish_wfn_utils
call finish_scratch
else
if(am_master)then
call wout('All particles fixed, no equilibration required.')
call wout()
endif
endif
do xyzzyaaab7=1,netot
if(xyzzyaaal1(xyzzyaaab7))xyzzyaaas7(1:3,xyzzyaaab7)=xyzzyaaag1(1:3,xy&
&zzyaaab7)
enddo
if(isperiodic)call map_to_simcell(3,netot,xyzzyaaas7)
if(am_master)then
call wout('Particle positions for plot (paste into PLOT block to fix):&
&')
xyzzyaaab7=0
do xyzzyaaad7=1,nspin
do xyzzyaaac7=1,nele(xyzzyaaad7)
xyzzyaaab7=xyzzyaaab7+1
write(tmpr,'(1x,a,'//trim(i2s(dimensionality))//'(1x,es15.8))')'fix el&
&ectron '//trim(i2s(xyzzyaaac7))//' spin '//trim(i2s(xyzzyaaad7))//' @&
&',xyzzyaaas7(1:dimensionality,xyzzyaaab7)
call wout(tmpr)
enddo
enddo
call wout()
call wordwrap('Evaluating '//trim(now_plotting_name)//' at '//trim(i2s&
&(tot_npoint))//' positions of particle '//trim(i2s(which_ie(plot_ii))&
&)//' of type '//trim(i2s(which_spin(plot_ii)))//' in a '//trim(i2s(pl&
&ot_dim))//'-dimensional region.')
call wout()
endif
xyzzyaaai7=0
xyzzyaaaj7=0
call scratch_request(ratio1_from=xyzzyaaai7,ratio1_to=xyzzyaaaj7)
call energy_scratch_request(xyzzyaaai7)
call energy_scratch_request(xyzzyaaaj7)
call scratch_request(kinetic_detail=xyzzyaaai7)
call scratch_request(kinetic_detail=xyzzyaaaj7)
call scratch_request(wfn_detail=xyzzyaaai7)
call scratch_request(wfn_detail=xyzzyaaaj7)
call setup_scratch
call which_scratch(xyzzyaaai7)
call which_scratch(xyzzyaaaj7)
call setup_wfn_utils
call setup_energy_utils
call define_config(xyzzyaaai7,xyzzyaaas7,xyzzyaaam7)
do xyzzyaaae7=iplotmin,iplotmax
xyzzyaaaq7(1:3)=grid_point(1:3,xyzzyaaae7)
xyzzyaaas7(1:3,plot_ii)=xyzzyaaaq7(1:3)
select case(what_to_plot_sec)
case(1)
call define_config_oneelec(plot_ii,xyzzyaaai7,xyzzyaaaj7,xyzzyaaaq7,xy&
&zzyaaam7(1))
call get_plot_wfn(what_to_plot,plot_ii,xyzzyaaai7,xyzzyaaaj7,plot_data&
&(1,1,xyzzyaaae7))
case(2)
call define_config_oneelec(plot_ii,xyzzyaaai7,xyzzyaaaj7,xyzzyaaaq7,xy&
&zzyaaam7(1))
call accept_move(xyzzyaaai7,xyzzyaaaj7)
call get_plot_energy(what_to_plot,xyzzyaaai7,plot_data(1,1,xyzzyaaae7)&
&)
case default
call errstop_master('PLOT_MAIN','Don''t know what I''m plotting. Bug.'&
&)
end select
enddo
if(am_master)then
do inode=1,nnodes-1
iplotmin=nint(dble(tot_npoint)/dble(nnodes)*dble(inode))+1
if(inode==nnodes-1)then
iplotmax=tot_npoint
else
iplotmax=nint(dble(tot_npoint)/dble(nnodes)*dble(inode+1))
endif
call mpi_recv(plot_data(1:plot_ncomponents,1:plot_nfunctions,iplotmin:&
&iplotmax),plot_ncomponents*plot_nfunctions*(iplotmax-iplotmin+1),mpi_&
&double_precision,inode,2,mpi_comm_world,status,ierror)
call checkmpi(ierror,'recv plot_data in plot_main.')
enddo
else
call mpi_ssend(plot_data(1:plot_ncomponents,1:plot_nfunctions,iplotmin&
&:iplotmax),plot_ncomponents*plot_nfunctions*(iplotmax-iplotmin+1),mpi&
&_double_precision,0,2,mpi_comm_world,ierror)
endif
if(am_master)then
call wout('Writing functions:')
do ifn=1,plot_nfunctions
if(plot_nfunctions==1)then
filename='casinoplot.dat'
else
filename='casinoplot_'//trim(pad_int(ifn,plot_nfunctions,'0'))//'.dat'
endif
call wordwrap('* '//trim(filename)//': '//trim(function_name(ifn)))
call open_units(xyzzyaaaa7,xyzzyaaak7)
if(xyzzyaaak7/=0)call errstop_master('PLOT_MAIN','Unable to find free &
&I/O unit for plot file.')
open(xyzzyaaaa7,file=filename,status='unknown',form='formatted')
write(xyzzyaaaa7,'(a)')'# CASINO PLOT'
call xyzzyaaap1('plot-title',function_name(ifn))
call xyzzyaaap1('system-dimensionality',i2s(dimensionality))
call xyzzyaaap1('function-rank',i2s(plot_rank))
if(plot_rank>0)then
call xyzzyaaap1('function-rotate-as-tensor',l2s(xyzzyaaaj1))
endif
if(plot_rank==1)then
call xyzzyaaap1('function-translation-invariant',l2s(.not.xyzzyaaak1))
endif
call xyzzyaaap1('function-complex',l2s(xyzzyaaah1))
call xyzzyaaap1('function-stderr',l2s(xyzzyaaai1))
call xyzzyaaap1('region-dimensionality',i2s(plot_dim))
write(line,*)xyzzyaaac1(1:dimensionality)
call xyzzyaaap1('region-a',line)
write(line,*)xyzzyaaad1(1:dimensionality)
call xyzzyaaap1('region-b',line)
if(plot_dim>1)then
write(line,*)xyzzyaaae1(1:dimensionality)
call xyzzyaaap1('region-c',adjustl(line))
if(plot_dim>2)then
write(line,*)xyzzyaaaf1(1:dimensionality)
call xyzzyaaap1('region-d',adjustl(line))
endif
endif
call xyzzyaaap1('region-ngrid-ab',i2s(grid_npoint_ab))
if(plot_dim>1)then
call xyzzyaaap1('region-ngrid-ac',i2s(grid_npoint_ac))
if(plot_dim>2)then
call xyzzyaaap1('region-ngrid-ad',i2s(grid_npoint_ad))
endif
endif
xyzzyaaae7=0
do xyzzyaaah7=1,grid_npoint_ad
do xyzzyaaag7=1,grid_npoint_ac
do xyzzyaaaf7=1,grid_npoint_ab
xyzzyaaae7=xyzzyaaae7+1
write(line,*)grid_point(1:dimensionality,xyzzyaaae7),plot_data(1:plot_&
&ncomponents,ifn,xyzzyaaae7)
write(xyzzyaaaa7,'(1x,a)')trim(line)
enddo
enddo
enddo
close(xyzzyaaaa7)
open_unit(xyzzyaaaa7)=.false.
enddo
call wout()
call wout('Finished.')
call wout()
endif
deallocate(plot_data,grid_point,map_point)
deallocate(xyzzyaaas7,xyzzyaaam7)
call finish_energy_utils
call finish_wfn_utils
call finish_scratch
call timer('PLOT_MAIN',.false.)
contains
subroutine xyzzyaaap1(keyword,value)
implicit none
character(*),intent(in) :: keyword,value
write(xyzzyaaaa7,'(a)')'# '//trim(adjustl(keyword))//': '//trim(adjust&
&l(value))
end subroutine xyzzyaaap1
end subroutine plot_main
end module slaarnaby

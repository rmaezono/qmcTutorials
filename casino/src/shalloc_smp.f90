module shalloc
use comms
use dsp
use format_utils,only   : i2s,i2s64,d2s,wout
use parallel,only       : nnodes,am_master,checkmpi
use slaarnacc,only : ranx_max
use run_control,only    : check_alloc,errstop,errstop_master
use store,only          : chkpoint_level,netot
implicit none
private
public shallocate,deshallocate,shallocate_barrier,need_shm,init_shm,am&
&_smpmaster,smp_nodes,nsmps,nnpsmp,shallocate_blip,deshallocate_blip,s&
&hallocate_gauss,shm_size,shm_size_nproc,shalloc_clean,my_smpgroup,my_&
&smpmaster,my_smpproc,smp_masters
logical,parameter :: need_shm=.true.
integer xyzzyaaaa1,xyzzyaaab1,xyzzyaaac1
integer xyzzyaaad1
logical am_smpmaster
logical,parameter :: xyzzyaaae1=.false.
integer nsmps
integer nnpsmp
integer,allocatable :: smp_nodes(:)
integer,allocatable :: smp_masters(:)
integer my_smpgroup,my_smpmaster,my_smpproc
integer(kind=i64) xyzzyaaaf1,xyzzyaaag1
integer shm_size_nproc
character(256) tmpr
interface shallocate
module procedure xyzzyaaah1,xyzzyaaal1,xyzzyaabf1,xyzzyaaap1,xyzzyaaat&
&1,xyzzyaaax1,xyzzyaabj1,xyzzyaabb1,xyzzyaabr1,xyzzyaabn1
end interface
interface deshallocate
module procedure xyzzyaaaj1,xyzzyaaan1,xyzzyaabh1,xyzzyaaar1,xyzzyaaav&
&1,xyzzyaaaz1,xyzzyaabl1,xyzzyaabd1,xyzzyaabt1,xyzzyaabp1
end interface
interface shallocate_blip
module procedure xyzzyaabv1,xyzzyaabx1,xyzzyaabz1,xyzzyaacb1
end interface
interface deshallocate_blip
module procedure xyzzyaaar1,xyzzyaaav1,xyzzyaaaz1,xyzzyaabl1,xyzzyaabd&
&1,xyzzyaabt1,xyzzyaabp1
end interface
interface shallocate_gauss
module procedure xyzzyaacd1
end interface
external alloc_shm
external dealloc_shm
external get_smp_list
external get_nnpsmp
external set_shm_debug
contains
subroutine init_shm(nnodes,my_node)
implicit none
integer,intent(in) :: nnodes,my_node
integer xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7,xyzzyaaad7,xyzzyaaae7
integer,allocatable :: xyzzyaaaf7(:)
xyzzyaaad1=my_node
if(xyzzyaaae1)call set_shm_debug()
call get_nnpsmp(xyzzyaaab7)
if(xyzzyaaab7>0)then
xyzzyaaaa7=nnodes/xyzzyaaab7
if(xyzzyaaaa7*xyzzyaaab7/=nnodes)call errstop_master('INIT_SHM','Total&
& number of MPI processes must be exactly divisible by the number of p&
&rocesses over which to share memory, but it is not. Change the latter&
& using the --shm/--shmem flags to runqmc or by setting an environment&
& variable CASINO_NUMABLK.)')
allocate(smp_masters(xyzzyaaaa7),smp_nodes(xyzzyaaab7),stat=xyzzyaaae7&
&)
else
xyzzyaaaa7=nnodes
xyzzyaaab7=nnodes
allocate(smp_masters(xyzzyaaaa7),smp_nodes(xyzzyaaab7),stat=xyzzyaaae7&
&)
endif
call check_alloc(xyzzyaaae7,'INIT_SHM','smp_masters, smp_nodes')
call get_smp_list(nsmps,smp_masters,nnpsmp,smp_nodes)
am_smpmaster=my_node==smp_nodes(1)
if(nsmps*nnpsmp/=nnodes)call errstop_master('INIT_SHM','Total number o&
&f MPI processes must be exactly divisible by the number of processes &
&over which to share memory, but it is not. In this case -- in the abs&
&ence of user- or environment-supplied information (you''re not using &
&runqmc, right?) -- the grouping of processor cores into nodes was ded&
&uced from the processor names obtained using mpi_get_processor_name. &
&You may change this by using the --shm/--shmem flags to runqmc, or by&
& setting an environment variable CASINO_NUMABLK.)')
if(nsmps<nnodes)then
allocate(xyzzyaaaf7(nsmps),stat=xyzzyaaae7)
call check_alloc(xyzzyaaae7,'INIT_SHM','shm_temp')
xyzzyaaaf7(1:nsmps)=smp_masters(1:nsmps)
deallocate(smp_masters)
allocate(smp_masters(nsmps),stat=xyzzyaaae7)
call check_alloc(xyzzyaaae7,'INIT_SHM','smp_masters <2>')
smp_masters=xyzzyaaaf7
deallocate(xyzzyaaaf7)
endif
if(nnpsmp<nnodes)then
allocate(xyzzyaaaf7(nnpsmp),stat=xyzzyaaae7)
call check_alloc(xyzzyaaae7,'INIT_SHM','shm_temp <2>')
xyzzyaaaf7(1:nnpsmp)=smp_nodes(1:nnpsmp)
deallocate(smp_nodes)
allocate(smp_nodes(nnpsmp),stat=xyzzyaaae7)
call check_alloc(xyzzyaaae7,'INIT_SHM','smp_nodes <2>')
smp_nodes=xyzzyaaaf7
deallocate(xyzzyaaaf7)
endif
my_smpgroup=(my_node/nnpsmp)+1
my_smpmaster=(my_smpgroup-1)*nnpsmp
if(my_smpmaster==0)my_smpproc=my_node+1
if(my_smpmaster/=0)my_smpproc=mod(my_node,my_smpmaster)+1
call mpi_comm_group(mpi_comm_world,xyzzyaaaa1,xyzzyaaae7)
call checkmpi(xyzzyaaae7,'mpi_comm_group in init_shm')
call mpi_group_incl(xyzzyaaaa1,nsmps,smp_masters,xyzzyaaac7,xyzzyaaae7&
&)
call checkmpi(xyzzyaaae7,'mpi_group_incl in init_shm')
call mpi_comm_create(mpi_comm_world,xyzzyaaac7,xyzzyaaab1,xyzzyaaae7)
call checkmpi(xyzzyaaae7,'mpi_comm_create in init_shm')
call mpi_group_incl(xyzzyaaaa1,nnpsmp,smp_nodes,xyzzyaaad7,xyzzyaaae7)
call checkmpi(xyzzyaaae7,'mpi_group_incl <2> in init_shm')
call mpi_comm_create(mpi_comm_world,xyzzyaaad7,xyzzyaaac1,xyzzyaaae7)
call checkmpi(xyzzyaaae7,'mpi_comm_create <2> in init_shm')
xyzzyaaaf1=0_i64
xyzzyaaag1=0_i64
end subroutine init_shm
subroutine shallocate_barrier
implicit none
integer xyzzyaaaa8
call mpi_barrier(xyzzyaaac1,xyzzyaaaa8)
call checkmpi(xyzzyaaaa8,'mpi_barrier in shallocate_barrier')
end subroutine shallocate_barrier
subroutine xyzzyaaah1(a,i,stat)
implicit none
integer,pointer :: a(:)
integer,intent(in) :: i
integer,intent(out),optional :: stat
integer(mpi_address_kind) :: xyzzyaaaa9
integer xyzzyaaab9(i)
pointer (p_aa,xyzzyaaab9)
integer xyzzyaaac9
integer(kind=i64) :: xyzzyaaad9
xyzzyaaac9=0
if(present(stat))xyzzyaaac9=1
xyzzyaaad9=int(i,i64)
if(am_master)then
xyzzyaaaf1=xyzzyaaaf1+xyzzyaaad9*4_i64
xyzzyaaag1=max(xyzzyaaag1,xyzzyaaaf1)
endif
call alloc_shm(p_aa,xyzzyaaad9,mpi_integer,xyzzyaaac9)
call xyzzyaaae9(xyzzyaaab9)
if(present(stat))stat=xyzzyaaac9
if(xyzzyaaae1)then
xyzzyaaaa9=loc(a(1))
write(tmpr,"('pe',i2,': i_1: p_aa=',z16,', l_a=',z16)")xyzzyaaad1,p_aa&
&,xyzzyaaaa9
call wout(tmpr)
endif
contains
subroutine xyzzyaaae9(xyzzyaaab9)
implicit none
integer,target :: xyzzyaaab9(i)
a=>xyzzyaaab9
end subroutine xyzzyaaae9
end subroutine xyzzyaaah1
subroutine xyzzyaaaj1(a)
implicit none
integer,pointer :: a(:)
integer(mpi_address_kind) xyzzyaaaa11
call xyzzyaaab11(a)
call dealloc_shm(xyzzyaaaa11,int(size(a),i64),mpi_integer)
if(am_master)xyzzyaaaf1=xyzzyaaaf1-int(size(a),i64)*4_i64
nullify(a)
contains
subroutine xyzzyaaab11(a)
implicit none
integer,intent(in) :: a(:)
xyzzyaaaa11=loc(a)
end subroutine xyzzyaaab11
end subroutine xyzzyaaaj1
subroutine xyzzyaaal1(a,i,j,stat)
implicit none
integer,pointer :: a(:,:)
integer,intent(in) :: i,j
integer,intent(out),optional :: stat
integer(mpi_address_kind) :: xyzzyaaaa13
integer xyzzyaaab13(i,j)
pointer (p_aa,xyzzyaaab13)
integer xyzzyaaac13
integer(kind=i64) :: xyzzyaaad13
xyzzyaaac13=0
if(present(stat))xyzzyaaac13=1
xyzzyaaad13=int(i,i64)*int(j,i64)
if(am_master)then
xyzzyaaaf1=xyzzyaaaf1+xyzzyaaad13*4_i64
xyzzyaaag1=max(xyzzyaaag1,xyzzyaaaf1)
endif
call alloc_shm(p_aa,xyzzyaaad13,mpi_integer,xyzzyaaac13)
call xyzzyaaae13(xyzzyaaab13)
if(present(stat))stat=xyzzyaaac13
if(xyzzyaaae1)then
xyzzyaaaa13=loc(a(1,1))
write(tmpr,"('pe',i2,': i_2: p_aa=',z16,', l_a=',z16)")xyzzyaaad1,p_aa&
&,xyzzyaaaa13
call wout(tmpr)
endif
contains
subroutine xyzzyaaae13(xyzzyaaab13)
implicit none
integer,target :: xyzzyaaab13(i,j)
a=>xyzzyaaab13
end subroutine xyzzyaaae13
end subroutine xyzzyaaal1
subroutine xyzzyaaan1(a)
implicit none
integer,pointer :: a(:,:)
integer(mpi_address_kind) xyzzyaaaa15
call xyzzyaaab15(a)
call dealloc_shm(xyzzyaaaa15,int(size(a),i64),mpi_integer)
if(am_master)xyzzyaaaf1=xyzzyaaaf1-int(size(a),i64)*4_i64
nullify(a)
contains
subroutine xyzzyaaab15(a)
implicit none
integer,intent(in) :: a(:,:)
xyzzyaaaa15=loc(a)
end subroutine xyzzyaaab15
end subroutine xyzzyaaan1
subroutine xyzzyaaap1(a,i,stat)
implicit none
real(dp),pointer :: a(:)
integer,intent(in) :: i
integer,intent(out),optional :: stat
integer(mpi_address_kind) :: xyzzyaaaa17
real(dp) xyzzyaaab17(i)
pointer (p_aa,xyzzyaaab17)
integer xyzzyaaac17
xyzzyaaac17=0
if(present(stat))xyzzyaaac17=1
if(am_master)then
xyzzyaaaf1=xyzzyaaaf1+int(i,i64)*8_i64
xyzzyaaag1=max(xyzzyaaag1,xyzzyaaaf1)
endif
call alloc_shm(p_aa,int(i,i64),mpi_double_precision,xyzzyaaac17)
call xyzzyaaad17(xyzzyaaab17)
if(present(stat))stat=xyzzyaaac17
if(xyzzyaaae1)then
xyzzyaaaa17=loc(a)
write(tmpr,"('pe',i2,': rdp_1: p_aa=',z16,', l_a=',z16)")xyzzyaaad1,p_&
&aa,xyzzyaaaa17
call wout(tmpr)
endif
contains
subroutine xyzzyaaad17(xyzzyaaab17)
implicit none
real(dp),target :: xyzzyaaab17(i)
a=>xyzzyaaab17
end subroutine xyzzyaaad17
end subroutine xyzzyaaap1
subroutine xyzzyaaar1(a)
implicit none
real(dp),pointer :: a(:)
integer(mpi_address_kind) xyzzyaaaa19
call xyzzyaaab19(a)
if(am_master)xyzzyaaaf1=xyzzyaaaf1-int(size(a),i64)*8_i64
call dealloc_shm(xyzzyaaaa19,int(size(a),i64),mpi_double_precision)
nullify(a)
contains
subroutine xyzzyaaab19(a)
implicit none
real(dp),intent(in) :: a(:)
xyzzyaaaa19=loc(a)
end subroutine xyzzyaaab19
end subroutine xyzzyaaar1
subroutine xyzzyaaat1(a,i,j,stat)
implicit none
real(dp),pointer :: a(:,:)
integer,intent(in) :: i,j
integer,intent(out),optional :: stat
integer(mpi_address_kind) :: xyzzyaaaa21
real(dp) xyzzyaaab21(i,j)
pointer (p_aa,xyzzyaaab21)
integer xyzzyaaac21
integer(kind=i64) :: xyzzyaaad21
xyzzyaaac21=0
if(present(stat))xyzzyaaac21=1
xyzzyaaad21=int(i,i64)*int(j,i64)
if(am_master)then
xyzzyaaaf1=xyzzyaaaf1+xyzzyaaad21*8_i64
xyzzyaaag1=max(xyzzyaaag1,xyzzyaaaf1)
endif
call alloc_shm(p_aa,xyzzyaaad21,mpi_double_precision,xyzzyaaac21)
call xyzzyaaae21(xyzzyaaab21)
if(present(stat))stat=xyzzyaaac21
if(xyzzyaaae1)then
xyzzyaaaa21=loc(a(1,1))
write(tmpr,"('pe',i2,': rdp_2: p_aa=',z16,', l_a=',z16)")xyzzyaaad1,p_&
&aa,xyzzyaaaa21
call wout(tmpr)
endif
contains
subroutine xyzzyaaae21(xyzzyaaab21)
implicit none
real(dp),target :: xyzzyaaab21(i,j)
a=>xyzzyaaab21
end subroutine xyzzyaaae21
end subroutine xyzzyaaat1
subroutine xyzzyaaav1(a)
implicit none
real(dp),pointer :: a(:,:)
integer(mpi_address_kind) xyzzyaaaa23
call xyzzyaaab23(a)
call dealloc_shm(xyzzyaaaa23,int(size(a),i64),mpi_double_precision)
if(am_master)xyzzyaaaf1=xyzzyaaaf1-int(size(a),i64)*8_i64
nullify(a)
contains
subroutine xyzzyaaab23(a)
implicit none
real(dp),intent(in) :: a(:,:)
xyzzyaaaa23=loc(a)
end subroutine xyzzyaaab23
end subroutine xyzzyaaav1
subroutine xyzzyaaax1(a,i,j,k,stat)
implicit none
real(dp),pointer :: a(:,:,:)
integer,intent(in) :: i,j,k
integer,intent(out),optional :: stat
integer(mpi_address_kind) :: xyzzyaaaa25
real(dp) xyzzyaaab25(i,j,k)
pointer (p_aa,xyzzyaaab25)
integer xyzzyaaac25
integer(kind=i64) :: xyzzyaaad25
xyzzyaaac25=0
if(present(stat))xyzzyaaac25=1
xyzzyaaad25=int(i,i64)*int(j,i64)*int(k,i64)
if(am_master)then
xyzzyaaaf1=xyzzyaaaf1+xyzzyaaad25*8_i64
xyzzyaaag1=max(xyzzyaaag1,xyzzyaaaf1)
endif
call alloc_shm(p_aa,xyzzyaaad25,mpi_double_precision,xyzzyaaac25)
call xyzzyaaae25(xyzzyaaab25)
if(present(stat))stat=xyzzyaaac25
if(xyzzyaaae1)then
xyzzyaaaa25=loc(a(1,1,1))
write(tmpr,"('pe',i2,': rdp_3: p_aa=',z16,', l_a=',z16)")xyzzyaaad1,p_&
&aa,xyzzyaaaa25
call wout(tmpr)
endif
contains
subroutine xyzzyaaae25(xyzzyaaab25)
implicit none
real(dp),target :: xyzzyaaab25(i,j,k)
a=>xyzzyaaab25
end subroutine xyzzyaaae25
end subroutine xyzzyaaax1
subroutine xyzzyaaaz1(a)
implicit none
real(dp),pointer :: a(:,:,:)
integer(mpi_address_kind) xyzzyaaaa27
call xyzzyaaab27(a)
call dealloc_shm(xyzzyaaaa27,int(size(a),i64),mpi_double_precision)
if(am_master)xyzzyaaaf1=xyzzyaaaf1-int(size(a),i64)*8_i64
nullify(a)
contains
subroutine xyzzyaaab27(a)
implicit none
real(dp),intent(in) :: a(:,:,:)
xyzzyaaaa27=loc(a)
end subroutine xyzzyaaab27
end subroutine xyzzyaaaz1
subroutine xyzzyaabb1(a,i,j,k,l,stat)
implicit none
real(dp),pointer :: a(:,:,:,:)
integer,intent(in) :: i,j,k,l
integer,intent(out),optional :: stat
integer(mpi_address_kind) :: xyzzyaaaa29
real(dp) xyzzyaaab29(i,j,k,l)
pointer (p_aa,xyzzyaaab29)
integer xyzzyaaac29
integer(kind=i64) :: xyzzyaaad29
xyzzyaaac29=0
if(present(stat))xyzzyaaac29=1
xyzzyaaad29=int(i,i64)*int(j,i64)*int(k,i64)*int(l,i64)
if(am_master)then
xyzzyaaaf1=xyzzyaaaf1+xyzzyaaad29*8_i64
xyzzyaaag1=max(xyzzyaaag1,xyzzyaaaf1)
endif
call alloc_shm(p_aa,xyzzyaaad29,mpi_double_precision,xyzzyaaac29)
call xyzzyaaae29(xyzzyaaab29)
if(present(stat))stat=xyzzyaaac29
if(xyzzyaaae1)then
xyzzyaaaa29=loc(a(1,1,1,1))
write(tmpr,"('pe',i2,': rdp_4: p_aa=',z16,', l_a=',z16)")xyzzyaaad1,p_&
&aa,xyzzyaaaa29
call wout(tmpr)
endif
contains
subroutine xyzzyaaae29(xyzzyaaab29)
implicit none
real(dp),target :: xyzzyaaab29(i,j,k,l)
a=>xyzzyaaab29
end subroutine xyzzyaaae29
end subroutine xyzzyaabb1
subroutine xyzzyaabd1(a)
implicit none
real(dp),pointer :: a(:,:,:,:)
integer(mpi_address_kind) xyzzyaaaa31
call xyzzyaaab31(a)
call dealloc_shm(xyzzyaaaa31,int(size(a),i64),mpi_double_precision)
if(am_master)xyzzyaaaf1=xyzzyaaaf1-int(size(a),i64)*8_i64
nullify(a)
contains
subroutine xyzzyaaab31(a)
implicit none
real(dp),intent(in) :: a(:,:,:,:)
xyzzyaaaa31=loc(a)
end subroutine xyzzyaaab31
end subroutine xyzzyaabd1
subroutine xyzzyaabf1(a,i,stat)
implicit none
real(sp),pointer :: a(:)
integer,intent(in) :: i
integer,intent(out),optional :: stat
integer(mpi_address_kind) :: xyzzyaaaa33
real(sp) xyzzyaaab33(i)
pointer (p_aa,xyzzyaaab33)
integer xyzzyaaac33
xyzzyaaac33=0
if(present(stat))xyzzyaaac33=1
if(am_master)then
xyzzyaaaf1=xyzzyaaaf1+int(i,i64)*4_i64
xyzzyaaag1=max(xyzzyaaag1,xyzzyaaaf1)
endif
call alloc_shm(p_aa,int(i,i64),mpi_real,xyzzyaaac33)
call xyzzyaaad33(xyzzyaaab33)
if(present(stat))stat=xyzzyaaac33
if(xyzzyaaae1)then
xyzzyaaaa33=loc(a)
write(tmpr,"('pe',i2,': rsp_1: p_aa=',z16,', l_a=',z16)")xyzzyaaad1,p_&
&aa,xyzzyaaaa33
call wout(tmpr)
endif
contains
subroutine xyzzyaaad33(xyzzyaaab33)
implicit none
real(sp),target :: xyzzyaaab33(i)
a=>xyzzyaaab33
end subroutine xyzzyaaad33
end subroutine xyzzyaabf1
subroutine xyzzyaabh1(a)
implicit none
real(sp),pointer :: a(:)
integer(mpi_address_kind) xyzzyaaaa35
call xyzzyaaab35(a)
if(am_master)xyzzyaaaf1=xyzzyaaaf1-int(size(a),i64)*4_i64
call dealloc_shm(xyzzyaaaa35,int(size(a),i64),mpi_real)
nullify(a)
contains
subroutine xyzzyaaab35(a)
implicit none
real(sp),intent(in) :: a(:)
xyzzyaaaa35=loc(a)
end subroutine xyzzyaaab35
end subroutine xyzzyaabh1
subroutine xyzzyaabj1(a,i,j,k,l,stat)
implicit none
real(sp),pointer :: a(:,:,:,:)
integer,intent(in) :: i,j,k,l
integer,intent(out),optional :: stat
integer(mpi_address_kind) :: xyzzyaaaa37
real(sp) xyzzyaaab37(i,j,k,l)
pointer (p_aa,xyzzyaaab37)
integer xyzzyaaac37
integer(kind=i64) :: xyzzyaaad37
xyzzyaaac37=0
if(present(stat))xyzzyaaac37=1
xyzzyaaad37=int(i,i64)*int(j,i64)*int(k,i64)*int(l,i64)
if(am_master)then
xyzzyaaaf1=xyzzyaaaf1+xyzzyaaad37*4_i64
xyzzyaaag1=max(xyzzyaaag1,xyzzyaaaf1)
endif
call alloc_shm(p_aa,xyzzyaaad37,mpi_real,xyzzyaaac37)
call xyzzyaaae37(xyzzyaaab37)
if(present(stat))stat=xyzzyaaac37
if(xyzzyaaae1)then
xyzzyaaaa37=loc(a(1,1,1,1))
write(tmpr,"('pe',i2,': rsp_4: p_aa=',z16,', l_a=',z16)")xyzzyaaad1,p_&
&aa,xyzzyaaaa37
call wout(tmpr)
endif
contains
subroutine xyzzyaaae37(xyzzyaaab37)
implicit none
real(sp),target :: xyzzyaaab37(i,j,k,l)
a=>xyzzyaaab37
end subroutine xyzzyaaae37
end subroutine xyzzyaabj1
subroutine xyzzyaabl1(a)
implicit none
real(sp),pointer :: a(:,:,:,:)
integer(mpi_address_kind) xyzzyaaaa39
call xyzzyaaab39(a)
call dealloc_shm(xyzzyaaaa39,int(size(a),i64),mpi_real)
if(am_master)xyzzyaaaf1=xyzzyaaaf1-int(size(a),i64)*4_i64
nullify(a)
contains
subroutine xyzzyaaab39(a)
real(sp),intent(in) :: a(:,:,:,:)
xyzzyaaaa39=loc(a)
end subroutine xyzzyaaab39
end subroutine xyzzyaabl1
subroutine xyzzyaabn1(a,i,j,k,l,stat)
implicit none
complex(dp),pointer :: a(:,:,:,:)
integer,intent(in) :: i,j,k,l
integer,intent(out),optional :: stat
integer(mpi_address_kind)::xyzzyaaaa41
complex(dp) xyzzyaaab41(i,j,k,l)
pointer (p_aa,xyzzyaaab41)
integer xyzzyaaac41
integer(kind=i64) :: xyzzyaaad41
xyzzyaaac41=0
if(present(stat))xyzzyaaac41=1
xyzzyaaad41=int(i,i64)*int(j,i64)*int(k,i64)*int(l,i64)
if(am_master)then
xyzzyaaaf1=xyzzyaaaf1+xyzzyaaad41*16_i64
xyzzyaaag1=max(xyzzyaaag1,xyzzyaaaf1)
endif
call alloc_shm(p_aa,xyzzyaaad41,mpi_double_complex,xyzzyaaac41)
call xyzzyaaae41(xyzzyaaab41)
if(present(stat))stat=xyzzyaaac41
if(xyzzyaaae1)then
xyzzyaaaa41=loc(a(1,1,1,1))
write(tmpr,"('pe',i2,': cdp_4: p_aa=',z16,', l_a=',z16)")xyzzyaaad1,p_&
&aa,xyzzyaaaa41
call wout(tmpr)
endif
contains
subroutine xyzzyaaae41(xyzzyaaab41)
implicit none
complex(dp),target :: xyzzyaaab41(i,j,k,l)
a=>xyzzyaaab41
end subroutine xyzzyaaae41
end subroutine xyzzyaabn1
subroutine xyzzyaabp1(a)
implicit none
complex(dp),pointer :: a(:,:,:,:)
integer(mpi_address_kind) xyzzyaaaa43
call xyzzyaaab43(a)
call dealloc_shm(xyzzyaaaa43,int(size(a),i64),mpi_double_complex)
if(am_master)xyzzyaaaf1=xyzzyaaaf1-int(size(a),i64)*16_i64
nullify(a)
contains
subroutine xyzzyaaab43(a)
implicit none
complex(dp),intent(in) :: a(:,:,:,:)
xyzzyaaaa43=loc(a)
end subroutine xyzzyaaab43
end subroutine xyzzyaabp1
subroutine xyzzyaabr1(a,i,j,k,l,stat)
implicit none
complex(sp),pointer :: a(:,:,:,:)
integer,intent(in) :: i,j,k,l
integer,intent(out),optional :: stat
integer(mpi_address_kind) :: xyzzyaaaa45
complex(sp) xyzzyaaab45(i,j,k,l)
pointer (p_aa,xyzzyaaab45)
integer xyzzyaaac45
integer(kind=i64) :: xyzzyaaad45
xyzzyaaac45=0
if(present(stat))xyzzyaaac45=1
xyzzyaaad45=int(i,i64)*int(j,i64)*int(k,i64)*int(l,i64)
if(am_master)then
xyzzyaaaf1=xyzzyaaaf1+xyzzyaaad45*8_i64
xyzzyaaag1=max(xyzzyaaag1,xyzzyaaaf1)
endif
call alloc_shm(p_aa,xyzzyaaad45,mpi_complex,xyzzyaaac45)
call xyzzyaaae45(xyzzyaaab45)
if(present(stat))stat=xyzzyaaac45
if(xyzzyaaae1)then
xyzzyaaaa45=loc(a(1,1,1,1))
write(tmpr,"('pe',i2,': csp_4: p_aa=',z16,', l_a=',z16)")xyzzyaaad1,p_&
&aa,xyzzyaaaa45
call wout(tmpr)
endif
contains
subroutine xyzzyaaae45(xyzzyaaab45)
implicit none
complex(sp),target :: xyzzyaaab45(i,j,k,l)
a=>xyzzyaaab45
end subroutine xyzzyaaae45
end subroutine xyzzyaabr1
subroutine xyzzyaabt1(a)
implicit none
complex(sp),pointer :: a(:,:,:,:)
integer(mpi_address_kind) xyzzyaaaa47
call xyzzyaaab47(a)
call dealloc_shm(xyzzyaaaa47,int(size(a),i64),mpi_complex)
if(am_master)xyzzyaaaf1=xyzzyaaaf1-int(size(a),i64)*8_i64
nullify(a)
contains
subroutine xyzzyaaab47(a)
implicit none
complex(sp),intent(in) :: a(:,:,:,:)
xyzzyaaaa47=loc(a)
end subroutine xyzzyaaab47
end subroutine xyzzyaabt1
subroutine xyzzyaabv1(a,i,j,k,l,stat)
implicit none
real(sp),pointer :: a(:,:,:,:)
integer,intent(in) :: i,j,k,l
integer,intent(out),optional :: stat
integer(mpi_address_kind) :: xyzzyaaaa49
real(sp) xyzzyaaab49(i,j,k,l)
pointer (p_aa,xyzzyaaab49)
integer xyzzyaaac49
integer(kind=i64) :: xyzzyaaad49
xyzzyaaac49=0
if(present(stat))xyzzyaaac49=1
xyzzyaaad49=int(i,i64)*int(j,i64)*int(k,i64)*int(l,i64)
if(am_master)then
xyzzyaaaf1=xyzzyaaaf1+xyzzyaaad49*4_i64
xyzzyaaag1=max(xyzzyaaag1,xyzzyaaaf1)
endif
call alloc_shm(p_aa,xyzzyaaad49,mpi_real,xyzzyaaac49)
call xyzzyaaae49(xyzzyaaab49)
if(present(stat))stat=xyzzyaaac49
if(xyzzyaaae1)then
xyzzyaaaa49=loc(a(1,0,0,0))
write(tmpr,"('pe',i2,': rsp_4: p_aa=',z16,', l_a=',z16)")xyzzyaaad1,p_&
&aa,xyzzyaaaa49
call wout(tmpr)
endif
contains
subroutine xyzzyaaae49(xyzzyaaab49)
implicit none
real(sp),target :: xyzzyaaab49(i,0:j-1,0:k-1,0:l-1)
a=>xyzzyaaab49
end subroutine xyzzyaaae49
end subroutine xyzzyaabv1
subroutine xyzzyaabx1(a,i,j,k,l,stat)
implicit none
real(dp),pointer :: a(:,:,:,:)
integer,intent(in) :: i,j,k,l
integer,intent(out),optional :: stat
integer(mpi_address_kind) :: xyzzyaaaa51
real(dp) xyzzyaaab51(i,j,k,l)
pointer (p_aa,xyzzyaaab51)
integer xyzzyaaac51
integer(kind=i64) :: xyzzyaaad51
xyzzyaaac51=0
if(present(stat))xyzzyaaac51=1
xyzzyaaad51=int(i,i64)*int(j,i64)*int(k,i64)*int(l,i64)
if(am_master)then
xyzzyaaaf1=xyzzyaaaf1+xyzzyaaad51*8_i64
xyzzyaaag1=max(xyzzyaaag1,xyzzyaaaf1)
endif
call alloc_shm(p_aa,xyzzyaaad51,mpi_double_precision,xyzzyaaac51)
call xyzzyaaae51(xyzzyaaab51)
if(present(stat))stat=xyzzyaaac51
if(xyzzyaaae1)then
xyzzyaaaa51=loc(a(1,0,0,0))
write(tmpr,"('pe',i2,': rdp_4: p_aa=',z16,', l_a=',z16)")xyzzyaaad1,p_&
&aa,xyzzyaaaa51
call wout(tmpr)
endif
contains
subroutine xyzzyaaae51(xyzzyaaab51)
implicit none
real(dp),target :: xyzzyaaab51(i,0:j-1,0:k-1,0:l-1)
a=>xyzzyaaab51
end subroutine xyzzyaaae51
end subroutine xyzzyaabx1
subroutine xyzzyaabz1(a,i,j,k,l,stat)
implicit none
complex(sp),pointer :: a(:,:,:,:)
integer,intent(in) :: i,j,k,l
integer,intent(out),optional :: stat
integer(mpi_address_kind) :: xyzzyaaaa53
complex(sp) xyzzyaaab53(i,j,k,l)
pointer (p_aa,xyzzyaaab53)
integer xyzzyaaac53
integer(kind=i64) :: xyzzyaaad53
xyzzyaaac53=0
if(present(stat))xyzzyaaac53=1
xyzzyaaad53=int(i,i64)*int(j,i64)*int(k,i64)*int(l,i64)
if(am_master)then
xyzzyaaaf1=xyzzyaaaf1+xyzzyaaad53*8_i64
xyzzyaaag1=max(xyzzyaaag1,xyzzyaaaf1)
endif
call alloc_shm(p_aa,xyzzyaaad53,mpi_complex,xyzzyaaac53)
call xyzzyaaae53(xyzzyaaab53)
if(present(stat))stat=xyzzyaaac53
if(xyzzyaaae1)then
xyzzyaaaa53=loc(a(1,0,0,0))
write(tmpr,"('pe',i2,': csp_4: p_aa=',z16,', l_a=',z16)")xyzzyaaad1,p_&
&aa,xyzzyaaaa53
call wout(tmpr)
endif
contains
subroutine xyzzyaaae53(xyzzyaaab53)
implicit none
complex(sp),target :: xyzzyaaab53(i,0:j-1,0:k-1,0:l-1)
a=>xyzzyaaab53
end subroutine xyzzyaaae53
end subroutine xyzzyaabz1
subroutine xyzzyaacb1(a,i,j,k,l,stat)
implicit none
complex(dp),pointer :: a(:,:,:,:)
integer,intent(in) :: i,j,k,l
integer,intent(out),optional :: stat
integer(mpi_address_kind)::xyzzyaaaa55
complex(dp) xyzzyaaab55(i,j,k,l)
pointer (p_aa,xyzzyaaab55)
integer xyzzyaaac55
integer(kind=i64) :: xyzzyaaad55
xyzzyaaac55=0
if(present(stat))xyzzyaaac55=1
xyzzyaaad55=int(i,i64)*int(j,i64)*int(k,i64)*int(l,i64)
if(am_master)then
xyzzyaaaf1=xyzzyaaaf1+xyzzyaaad55*16_i64
xyzzyaaag1=max(xyzzyaaag1,xyzzyaaaf1)
endif
call alloc_shm(p_aa,xyzzyaaad55,mpi_double_complex,xyzzyaaac55)
call xyzzyaaae55(xyzzyaaab55)
if(present(stat))stat=xyzzyaaac55
if(xyzzyaaae1)then
xyzzyaaaa55=loc(a(1,0,0,0))
write(tmpr,"('pe',i2,': cdp_4: p_aa=',z16,', l_a=',z16)")xyzzyaaad1,p_&
&aa,xyzzyaaaa55
call wout(tmpr)
endif
contains
subroutine xyzzyaaae55(xyzzyaaab55)
implicit none
complex(dp),target :: xyzzyaaab55(i,0:j-1,0:k-1,0:l-1)
a=>xyzzyaaab55
end subroutine xyzzyaaae55
end subroutine xyzzyaacb1
subroutine xyzzyaacd1(a,i,j,k1,k2,l,stat)
implicit none
complex(dp),pointer :: a(:,:,:,:)
integer,intent(in) :: i,j,k1,k2,l
integer,intent(out),optional :: stat
integer(mpi_address_kind)::xyzzyaaaa57
complex(dp) xyzzyaaab57(i,j,k1:k2,l)
pointer (p_aa,xyzzyaaab57)
integer xyzzyaaac57
integer(kind=i64) :: xyzzyaaad57
xyzzyaaac57=0
if(present(stat))xyzzyaaac57=1
xyzzyaaad57=int(i,i64)*int(j,i64)*int((k2-k1+1),i64)*int(l,i64)
if(am_master)then
xyzzyaaaf1=xyzzyaaaf1+xyzzyaaad57*16_i64
xyzzyaaag1=max(xyzzyaaag1,xyzzyaaaf1)
endif
call alloc_shm(p_aa,xyzzyaaad57,mpi_double_complex,xyzzyaaac57)
call xyzzyaaae57(xyzzyaaab57)
if(present(stat))stat=xyzzyaaac57
if(xyzzyaaae1)then
xyzzyaaaa57=loc(a(1,1,k1,1))
write(tmpr,"('pe',i2,': cdp_4: p_aa=',z16,', l_a=',z16)")xyzzyaaad1,p_&
&aa,xyzzyaaaa57
call wout(tmpr)
endif
contains
subroutine xyzzyaaae57(xyzzyaaab57)
implicit none
complex(dp),target :: xyzzyaaab57(i,j,k1:k2,l)
a=>xyzzyaaab57
end subroutine xyzzyaaae57
end subroutine xyzzyaacd1
subroutine shm_size
implicit none
integer xyzzyaaaa59
integer(kind=i64) xyzzyaaab59,xyzzyaaac59
if(shm_size_nproc==0)then
xyzzyaaaa59=nnpsmp
else
xyzzyaaaa59=shm_size_nproc
endif
if(chkpoint_level/=-1)then
xyzzyaaab59=int(max(4*ranx_max*xyzzyaaaa59,8*3*netot*xyzzyaaaa59),i64)
else
xyzzyaaab59=0_i64
endif
xyzzyaaac59=max(xyzzyaaag1,xyzzyaaaf1+xyzzyaaab59)
call wout('Shared memory')
call wout('=============')
call wout()
if(chkpoint_level==-1)then
call wout('Config writes turned off by CHECKPOINT keyword; Shm size re&
&duced.')
call wout()
else
if(shm_size_nproc>0)then
call wout('Current number of processes/node overridden by SHM_SIZE_NPR&
&OC in input.')
call wout('Assuming '//trim(i2s(xyzzyaaaa59))//' MPI processes/node.')
call wout()
endif
endif
if(xyzzyaaac59/=0_i64)then
call wout('Maximum shared memory required : '//trim(i2s64(xyzzyaaac59)&
&)//' bytes = '//trim(d2s(real(xyzzyaaac59,dp)*1.d-6))//' MB')
call wout('(BG_SHAREDMEMSIZE via --user.shmemsize on Blue Gene machine&
&s)')
else
call wout('For this system with this input no shared memory will be us&
&ed.')
endif
call wout()
end subroutine shm_size
subroutine shalloc_clean
implicit none
call clean_shm()
end subroutine shalloc_clean
end module shalloc

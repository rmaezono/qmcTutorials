module shalloc
use dsp
implicit none
private
public shallocate,shallocate_blip,shallocate_gauss,deshallocate,deshal&
&locate_blip,shallocate_barrier,need_shm,init_shm,am_smpmaster,smp_nod&
&es,nsmps,nnpsmp,shm_size,shm_size_nproc,shalloc_clean,my_smpgroup,my_&
&smpmaster,my_smpproc,smp_masters
logical,parameter :: need_shm=.false.
logical :: am_smpmaster=.true.
integer nsmps
integer :: nnpsmp=1
integer,allocatable :: smp_nodes(:)
integer,allocatable :: smp_masters(:)
integer shm_size_nproc
integer :: my_smpgroup,my_smpmaster,my_smpproc=0
interface shallocate
module procedure xyzzyaaaa1,xyzzyaaac1,xyzzyaaae1,xyzzyaaai1,xyzzyaaak&
&1,xyzzyaaam1,xyzzyaaag1,xyzzyaaao1,xyzzyaaaq1,xyzzyaaas1
end interface
interface deshallocate
module procedure xyzzyaaab1,xyzzyaaad1,xyzzyaaaf1,xyzzyaaaj1,xyzzyaaal&
&1,xyzzyaaan1,xyzzyaaah1,xyzzyaaap1,xyzzyaaar1,xyzzyaaat1
end interface
interface shallocate_blip
module procedure xyzzyaaau1,xyzzyaaav1,xyzzyaaaw1,xyzzyaaax1
end interface
interface deshallocate_blip
module procedure xyzzyaaaj1,xyzzyaaal1,xyzzyaaan1,xyzzyaaah1,xyzzyaaap&
&1,xyzzyaaar1,xyzzyaaat1
end interface
interface shallocate_gauss
module procedure xyzzyaaay1
end interface
contains
subroutine init_shm(nnodes,my_node)
implicit none
integer,intent(in) :: nnodes,my_node
nsmps=nnodes
my_smpgroup=my_node
my_smpmaster=my_node
end subroutine init_shm
subroutine shallocate_barrier
implicit none
end subroutine shallocate_barrier
subroutine xyzzyaaaa1(a,i,stat)
implicit none
integer,pointer :: a(:)
integer,intent(in) :: i
integer,intent(out),optional :: stat
integer xyzzyaaaa9
allocate(a(i),stat=xyzzyaaaa9)
if(present(stat))stat=xyzzyaaaa9
end subroutine xyzzyaaaa1
subroutine xyzzyaaab1(a)
implicit none
integer,pointer :: a(:)
deallocate(a)
end subroutine xyzzyaaab1
subroutine xyzzyaaac1(a,i,j,stat)
implicit none
integer,pointer :: a(:,:)
integer,intent(in) :: i,j
integer,intent(out),optional :: stat
integer xyzzyaaaa11
allocate(a(i,j),stat=xyzzyaaaa11)
if(present(stat))stat=xyzzyaaaa11
end subroutine xyzzyaaac1
subroutine xyzzyaaad1(a)
implicit none
integer,pointer :: a(:,:)
deallocate(a)
end subroutine xyzzyaaad1
subroutine xyzzyaaae1(a,i,stat)
implicit none
real(sp),pointer :: a(:)
integer,intent(in) :: i
integer,intent(out),optional :: stat
integer xyzzyaaaa13
allocate(a(i),stat=xyzzyaaaa13)
if(present(stat))stat=xyzzyaaaa13
end subroutine xyzzyaaae1
subroutine xyzzyaaaf1(a)
implicit none
real(sp),pointer :: a(:)
deallocate(a)
end subroutine xyzzyaaaf1
subroutine xyzzyaaag1(a,i,j,k,l,stat)
implicit none
real(sp),pointer :: a(:,:,:,:)
integer,intent(in) :: i,j,k,l
integer,intent(out),optional :: stat
integer xyzzyaaaa15
allocate(a(i,j,k,l),stat=xyzzyaaaa15)
if(present(stat))stat=xyzzyaaaa15
end subroutine xyzzyaaag1
subroutine xyzzyaaah1(a)
implicit none
real(sp),pointer :: a(:,:,:,:)
deallocate(a)
end subroutine xyzzyaaah1
subroutine xyzzyaaai1(a,i,stat)
implicit none
real(dp),pointer :: a(:)
integer,intent(in) :: i
integer,intent(out),optional :: stat
integer xyzzyaaaa17
allocate(a(i),stat=xyzzyaaaa17)
if(present(stat))stat=xyzzyaaaa17
end subroutine xyzzyaaai1
subroutine xyzzyaaaj1(a)
implicit none
real(dp),pointer :: a(:)
deallocate(a)
end subroutine xyzzyaaaj1
subroutine xyzzyaaak1(a,i,j,stat)
implicit none
real(dp),pointer :: a(:,:)
integer,intent(in) :: i,j
integer,intent(out),optional :: stat
integer xyzzyaaaa19
allocate(a(i,j),stat=xyzzyaaaa19)
if(present(stat))stat=xyzzyaaaa19
end subroutine xyzzyaaak1
subroutine xyzzyaaal1(a)
implicit none
real(dp),pointer :: a(:,:)
deallocate(a)
end subroutine xyzzyaaal1
subroutine xyzzyaaam1(a,i,j,k,stat)
implicit none
real(dp),pointer :: a(:,:,:)
integer,intent(in) :: i,j,k
integer,intent(out),optional :: stat
integer xyzzyaaaa21
allocate(a(i,j,k),stat=xyzzyaaaa21)
if(present(stat))stat=xyzzyaaaa21
end subroutine xyzzyaaam1
subroutine xyzzyaaan1(a)
implicit none
real(dp),pointer :: a(:,:,:)
deallocate(a)
end subroutine xyzzyaaan1
subroutine xyzzyaaao1(a,i,j,k,l,stat)
implicit none
real(dp),pointer :: a(:,:,:,:)
integer,intent(in) :: i,j,k,l
integer,intent(out),optional :: stat
integer xyzzyaaaa23
allocate(a(i,j,k,l),stat=xyzzyaaaa23)
if(present(stat))stat=xyzzyaaaa23
end subroutine xyzzyaaao1
subroutine xyzzyaaap1(a)
implicit none
real(dp),pointer :: a(:,:,:,:)
deallocate(a)
end subroutine xyzzyaaap1
subroutine xyzzyaaaq1(a,i,j,k,l,stat)
implicit none
integer,intent(in) :: i,j,k,l
integer,intent(out),optional :: stat
complex(sp),pointer :: a(:,:,:,:)
integer xyzzyaaaa25
allocate(a(i,j,k,l),stat=xyzzyaaaa25)
if(present(stat))stat=xyzzyaaaa25
end subroutine xyzzyaaaq1
subroutine xyzzyaaar1(a)
implicit none
complex(sp),pointer :: a(:,:,:,:)
deallocate(a)
end subroutine xyzzyaaar1
subroutine xyzzyaaas1(a,i,j,k,l,stat)
implicit none
integer,intent(in) :: i,j,k,l
integer,intent(out),optional :: stat
complex(dp),pointer :: a(:,:,:,:)
integer xyzzyaaaa27
allocate(a(i,j,k,l),stat=xyzzyaaaa27)
if(present(stat))stat=xyzzyaaaa27
end subroutine xyzzyaaas1
subroutine xyzzyaaat1(a)
implicit none
complex(dp),pointer :: a(:,:,:,:)
deallocate(a)
end subroutine xyzzyaaat1
subroutine xyzzyaaau1(a,i,j,k,l,stat)
implicit none
real(sp),pointer :: a(:,:,:,:)
integer,intent(in) :: i,j,k,l
integer,intent(out),optional :: stat
integer xyzzyaaaa29
allocate(a(i,0:j-1,0:k-1,0:l-1),stat=xyzzyaaaa29)
if(present(stat))stat=xyzzyaaaa29
end subroutine xyzzyaaau1
subroutine xyzzyaaav1(a,i,j,k,l,stat)
implicit none
real(dp),pointer :: a(:,:,:,:)
integer,intent(in) :: i,j,k,l
integer,intent(out),optional :: stat
integer xyzzyaaaa30
allocate(a(i,0:j-1,0:k-1,0:l-1),stat=xyzzyaaaa30)
if(present(stat))stat=xyzzyaaaa30
end subroutine xyzzyaaav1
subroutine xyzzyaaaw1(a,i,j,k,l,stat)
implicit none
integer,intent(in) :: i,j,k,l
integer,intent(out),optional :: stat
complex(sp),pointer :: a(:,:,:,:)
integer xyzzyaaaa31
allocate(a(i,0:j-1,0:k-1,0:l-1),stat=xyzzyaaaa31)
if(present(stat))stat=xyzzyaaaa31
end subroutine xyzzyaaaw1
subroutine xyzzyaaax1(a,i,j,k,l,stat)
implicit none
integer,intent(in) :: i,j,k,l
integer,intent(out),optional :: stat
complex(dp),pointer :: a(:,:,:,:)
integer xyzzyaaaa32
allocate(a(i,0:j-1,0:k-1,0:l-1),stat=xyzzyaaaa32)
if(present(stat))stat=xyzzyaaaa32
end subroutine xyzzyaaax1
subroutine xyzzyaaay1(a,i,j,k1,k2,l,stat)
implicit none
integer,intent(in) :: i,j,k1,k2,l
integer,intent(out),optional :: stat
complex(dp),pointer :: a(:,:,:,:)
integer xyzzyaaaa33
allocate(a(i,j,k1:k2,l),stat=xyzzyaaaa33)
if(present(stat))stat=xyzzyaaaa33
end subroutine xyzzyaaay1
subroutine shm_size
implicit none
return
end subroutine shm_size
subroutine shalloc_clean
implicit none
return
end subroutine shalloc_clean
end module shalloc

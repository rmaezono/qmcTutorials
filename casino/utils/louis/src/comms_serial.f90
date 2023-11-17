module comms
use dsp
implicit none
private
integer,public :: mpi_bottom
integer,parameter,public :: mpi_success=0,mpi_err_buffer=1,mpi_err_cou&
&nt=2,mpi_err_type=3,mpi_err_tag=4,mpi_err_comm=5,mpi_err_rank=6,mpi_e&
&rr_root=7,mpi_err_group=8,mpi_err_op=9,mpi_err_topology=10,mpi_err_di&
&ms=11,mpi_err_arg=12,mpi_err_unknown=13,mpi_err_truncate=14,mpi_err_o&
&ther=15,mpi_err_intern=16,mpi_err_in_status=17,mpi_err_pending=18,mpi&
&_err_request=19,mpi_err_lastcode=4114,mpi_undefined=(-32766),mpi_grap&
&h=1,mpi_cart=2,mpi_proc_null=(-1),mpi_bsend_overhead=512,mpi_source=2&
&,mpi_tag=3,mpi_error=4,mpi_status_size=4,mpi_max_processor_name=256,m&
&pi_max_error_string=512,mpi_max_name_string=63,mpi_comm_null=0,mpi_da&
&tatype_null=0,mpi_errhandler_null=0,mpi_group_null=0,mpi_keyval_inval&
&id=0,mpi_request_null=0,mpi_ident=0,mpi_congruent=1,mpi_similar=2,mpi&
&_unequal=3,mpi_errors_are_fatal=119,mpi_errors_return=120,mpi_complex&
&=23,mpi_double_complex=24,mpi_logical=25,mpi_real=26,mpi_double_preci&
&sion=27,mpi_integer=28,mpi_2integer=29,mpi_2complex=30,mpi_2double_co&
&mplex=31,mpi_2real=32,mpi_2double_precision=33,mpi_character=1,mpi_by&
&te=3,mpi_ub=16,mpi_lb=15,mpi_packed=14,mpi_integer1=0,mpi_integer2=0,&
&mpi_integer4=0,mpi_real4=0,mpi_real8=0,mpi_max=100,mpi_min=101,mpi_su&
&m=102,mpi_prod=103,mpi_land=104,mpi_band=105,mpi_lor=106,mpi_bor=107,&
&mpi_lxor=108,mpi_bxor=109,mpi_minloc=110,mpi_maxloc=111,mpi_op_null=0&
&,mpi_group_empty=90,mpi_comm_world=91,mpi_comm_self=92,mpi_tag_ub=80,&
&mpi_host=82,mpi_io=84,mpi_wtime_is_global=86,mpi_any_source=(-2),mpi_&
&any_tag=(-1),mpi_version=1,mpi_subversion=1,mpi_offset_kind=4,mpi_mod&
&e_create=1,mpi_mode_wronly=1,mpi_mode_rdonly=1,mpi_info_null=0,mpi_se&
&ek_end=0,mpi_integer8=0,mpi_address_kind=8,mpi_status_ignore=1,mpi_st&
&atuses_ignore=1
public :: mpi_bcast
interface mpi_bcast
module procedure xyzzyaaaa1,xyzzyaaab1,xyzzyaaac1,xyzzyaaad1,xyzzyaaae&
&1,xyzzyaaaf1,xyzzyaaag1,xyzzyaaah1,xyzzyaaai1,xyzzyaaaj1,xyzzyaaak1,x&
&yzzyaaal1,xyzzyaaam1,xyzzyaaan1,xyzzyaaao1,xyzzyaaap1,xyzzyaaaq1,xyzz&
&yaaar1,xyzzyaaas1,xyzzyaaat1,xyzzyaaau1,xyzzyaaav1,xyzzyaaaw1,xyzzyaa&
&ax1,xyzzyaaay1,xyzzyaaaz1,xyzzyaaba1,xyzzyaabc1,xyzzyaabd1,xyzzyaabe1&
&,xyzzyaabf1,xyzzyaabg1,xyzzyaabb1
end interface
public :: mpi_gather
interface mpi_gather
module procedure xyzzyaabh1,xyzzyaabi1,xyzzyaabj1,xyzzyaabk1,xyzzyaabl&
&1,xyzzyaabm1,xyzzyaabn1,xyzzyaabo1,xyzzyaabp1,xyzzyaabq1,xyzzyaabr1,x&
&yzzyaabs1,xyzzyaabt1,xyzzyaabu1,xyzzyaabv1,xyzzyaabw1,xyzzyaabx1,xyzz&
&yaabz1,xyzzyaaca1,xyzzyaacb1,xyzzyaacc1,xyzzyaacd1,xyzzyaaby1
end interface
public :: mpi_gather_in_place
interface mpi_gather_in_place
module procedure xyzzyaacf1,xyzzyaace1
end interface
public :: mpi_recv
interface mpi_recv
module procedure xyzzyaacg1,xyzzyaach1,xyzzyaaci1,xyzzyaacj1,xyzzyaack&
&1,xyzzyaacl1,xyzzyaacm1,xyzzyaacn1,xyzzyaaco1,xyzzyaacp1,xyzzyaacq1,x&
&yzzyaacr1,xyzzyaacs1,xyzzyaact1,                                     &
&                      xyzzyaacu1,xyzzyaacv1,xyzzyaacw1,xyzzyaacx1,xyz&
&zyaacy1,    xyzzyaacz1,xyzzyaada1,xyzzyaadb1
end interface
public :: mpi_irecv
interface mpi_irecv
module procedure xyzzyaadc1,xyzzyaadd1,xyzzyaade1,xyzzyaadf1,xyzzyaadg&
&1,xyzzyaadh1,xyzzyaadi1,xyzzyaadj1,xyzzyaadk1,xyzzyaadl1,xyzzyaadm1,x&
&yzzyaadn1,xyzzyaado1,xyzzyaadp1,xyzzyaadq1,xyzzyaadr1,xyzzyaads1,xyzz&
&yaadt1,xyzzyaadu1,xyzzyaadv1
end interface
public :: mpi_file_read
interface mpi_file_read
module procedure xyzzyaadw1,xyzzyaadx1,xyzzyaady1,xyzzyaadz1,xyzzyaaea&
&1,xyzzyaaeb1,xyzzyaaec1,xyzzyaaed1,xyzzyaaee1,xyzzyaaef1,xyzzyaaeg1,x&
&yzzyaaeh1,xyzzyaaei1,xyzzyaaej1,xyzzyaaek1,xyzzyaael1,xyzzyaaem1,xyzz&
&yaaen1,xyzzyaaeo1,xyzzyaaep1,xyzzyaaeq1,xyzzyaaer1,xyzzyaaes1,xyzzyaa&
&et1
end interface
public :: mpi_file_read_at
interface mpi_file_read_at
module procedure xyzzyaaeu1,xyzzyaaev1,xyzzyaaew1,xyzzyaaex1,xyzzyaaey&
&1,xyzzyaaez1,xyzzyaafa1,xyzzyaafb1,xyzzyaafc1,xyzzyaafd1,xyzzyaafe1,x&
&yzzyaaff1,xyzzyaafg1,xyzzyaafh1,xyzzyaafi1,xyzzyaafj1,xyzzyaafk1,xyzz&
&yaafl1,xyzzyaafm1,xyzzyaafn1,xyzzyaafo1,xyzzyaafp1
end interface
public :: mpi_file_read_all
interface mpi_file_read_all
module procedure xyzzyaafq1,xyzzyaafr1,xyzzyaafs1,xyzzyaaft1
end interface
public :: mpi_file_write
interface mpi_file_write
module procedure xyzzyaafu1,xyzzyaafv1,xyzzyaafw1,xyzzyaafx1,xyzzyaafy&
&1,xyzzyaafz1,xyzzyaaga1,xyzzyaagb1,xyzzyaagc1,xyzzyaagd1,xyzzyaage1,x&
&yzzyaagf1,xyzzyaagg1,xyzzyaagh1,xyzzyaagi1,xyzzyaagj1,xyzzyaagk1,xyzz&
&yaagl1,xyzzyaagm1,xyzzyaagn1,xyzzyaago1,xyzzyaagp1,xyzzyaagq1
end interface
public :: mpi_file_write_at
interface mpi_file_write_at
module procedure xyzzyaagr1,xyzzyaags1,xyzzyaagt1,xyzzyaagu1,xyzzyaagv&
&1,xyzzyaagw1,xyzzyaagx1,xyzzyaagy1,xyzzyaagz1,xyzzyaaha1,xyzzyaahb1,x&
&yzzyaahc1,xyzzyaahd1,xyzzyaahe1,xyzzyaahf1,xyzzyaahg1,xyzzyaahh1,xyzz&
&yaahi1,xyzzyaahj1,xyzzyaahk1,xyzzyaahl1,xyzzyaahm1
end interface
public :: mpi_file_write_all
interface mpi_file_write_all
module procedure xyzzyaahn1,xyzzyaaho1,xyzzyaahp1,xyzzyaahq1
end interface
public :: mpi_reduce
interface mpi_reduce
module procedure xyzzyaahr1,xyzzyaahs1,xyzzyaaht1,xyzzyaahu1,xyzzyaahv&
&1,xyzzyaahw1,xyzzyaahx1,xyzzyaahy1,xyzzyaahz1,xyzzyaaia1,xyzzyaaic1,x&
&yzzyaaid1,xyzzyaaie1,xyzzyaaif1,xyzzyaaig1,xyzzyaaib1
end interface
public :: mpi_allreduce
interface mpi_allreduce
module procedure xyzzyaaih1,xyzzyaaii1,xyzzyaaij1,xyzzyaaik1,xyzzyaail&
&1,xyzzyaaim1,xyzzyaain1
end interface
public :: mpi_scatter
interface mpi_scatter
module procedure xyzzyaaio1,xyzzyaaip1,xyzzyaaiq1,xyzzyaair1,xyzzyaais&
&1,xyzzyaait1,xyzzyaaiu1,xyzzyaaiv1,xyzzyaaiw1,xyzzyaaix1,xyzzyaaiy1,x&
&yzzyaaiz1
end interface
public :: mpi_ssend
interface mpi_ssend
module procedure xyzzyaaja1,xyzzyaajb1,xyzzyaajc1,xyzzyaajd1,xyzzyaaje&
&1,xyzzyaajf1,xyzzyaajg1,xyzzyaajh1,xyzzyaaji1,xyzzyaajj1,xyzzyaajk1,x&
&yzzyaajl1,xyzzyaajm1,xyzzyaajn1,xyzzyaajo1,xyzzyaajp1,xyzzyaajq1,xyzz&
&yaajr1,xyzzyaajs1
end interface
public :: mpi_isend
interface mpi_isend
module procedure xyzzyaajt1,xyzzyaaju1,xyzzyaajv1,xyzzyaajw1,xyzzyaajx&
&1,xyzzyaajy1,xyzzyaajz1,xyzzyaaka1,xyzzyaakb1,xyzzyaakc1,xyzzyaakd1,x&
&yzzyaake1,xyzzyaakf1,xyzzyaakg1,xyzzyaakh1,xyzzyaaki1,xyzzyaakj1,xyzz&
&yaakk1
end interface
public :: mpi_issend
interface mpi_issend
module procedure xyzzyaakl1,xyzzyaakm1,xyzzyaakn1,xyzzyaako1,xyzzyaakp&
&1,xyzzyaakq1,xyzzyaakr1,xyzzyaaks1,xyzzyaakt1,xyzzyaaku1,xyzzyaakv1,x&
&yzzyaakw1,xyzzyaakx1,xyzzyaaky1,xyzzyaakz1,xyzzyaala1,xyzzyaalb1,xyzz&
&yaalc1
end interface
public :: mpi_alltoallv
interface mpi_alltoallv
module procedure xyzzyaalg1,xyzzyaalh1,xyzzyaali1,xyzzyaalj1
end interface
public :: mpi_address
interface mpi_address
module procedure xyzzyaalk1,xyzzyaall1
end interface
public :: mpi_allgather
interface mpi_allgather
module procedure xyzzyaalm1
end interface
public :: mpi_alltoall
interface mpi_alltoall
module procedure xyzzyaald1,xyzzyaale1,xyzzyaalf1
end interface
public :: mpi_abort,mpi_barrier,mpi_iprobe,mpi_comm_rank,mpi_comm_size&
&,mpi_comm_split,mpi_error_string,mpi_init,mpi_finalize,mpi_file_close&
&,mpi_file_get_byte_offset,mpi_file_get_position,mpi_file_open,mpi_fil&
&e_seek,mpi_file_set_view,mpi_type_commit,mpi_type_contiguous,mpi_type&
&_free,mpi_type_size,mpi_type_struct,mpi_type_vector,mpi_wtime,mpi_wai&
&t,mpi_waitall,mpi_comm_group,mpi_group_incl,mpi_comm_create,mpi_waits&
&ome
contains
subroutine mpi_abort(ipm_comm,ipm_in,ipm_info)
integer ipm_comm,ipm_in,ipm_info
ipm_info=mpi_success
end subroutine mpi_abort
subroutine mpi_barrier(ipm_comm,ipm_info)
integer ipm_comm,ipm_info
ipm_info=mpi_success
end subroutine mpi_barrier
subroutine mpi_iprobe(source,tag,comm,flag,status,ierror)
integer source,tag,comm,status(mpi_status_size),ierror
logical flag
ierror=mpi_success
flag=.false.
end subroutine mpi_iprobe
subroutine xyzzyaaaa1(buffer,count,datatype,root,comm,ierror)
complex(dp) buffer
integer count,datatype,root,comm,ierror
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaaaa1
subroutine xyzzyaaab1(buffer,count,datatype,root,comm,ierror)
complex(dp) buffer(:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaaab1
subroutine xyzzyaaac1(buffer,count,datatype,root,comm,ierror)
complex(dp) buffer(:,:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaaac1
subroutine xyzzyaaad1(buffer,count,datatype,root,comm,ierror)
complex(dp) buffer(:,:,:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaaad1
subroutine xyzzyaaae1(buffer,count,datatype,root,comm,ierror)
complex(dp) buffer(:,:,:,:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaaae1
subroutine xyzzyaaaf1(buffer,count,datatype,root,comm,ierror)
complex(dp) buffer(:,:,:,:,:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaaaf1
subroutine xyzzyaaag1(buffer,count,datatype,root,comm,ierror)
complex(dp) buffer(:,:,:,:,:,:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaaag1
subroutine xyzzyaaah1(buffer,count,datatype,root,comm,ierror)
complex(sp) buffer(:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_complex)ierror=mpi_err_type
end subroutine xyzzyaaah1
subroutine xyzzyaaai1(buffer,count,datatype,root,comm,ierror)
complex(sp) buffer(:,:,:,:,:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_complex)ierror=mpi_err_type
end subroutine xyzzyaaai1
subroutine xyzzyaaaj1(buffer,count,datatype,root,comm,ierror)
real(dp) buffer
integer count,datatype,root,comm,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaaaj1
subroutine xyzzyaaak1(buffer,count,datatype,root,comm,ierror)
real(dp) buffer(:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaaak1
subroutine xyzzyaaal1(buffer,count,datatype,root,comm,ierror)
real(dp) buffer(:,:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaaal1
subroutine xyzzyaaam1(buffer,count,datatype,root,comm,ierror)
real(dp) buffer(:,:,:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaaam1
subroutine xyzzyaaan1(buffer,count,datatype,root,comm,ierror)
real(dp) buffer(:,:,:,:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaaan1
subroutine xyzzyaaao1(buffer,count,datatype,root,comm,ierror)
real(dp) buffer(:,:,:,:,:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaaao1
subroutine xyzzyaaap1(buffer,count,datatype,root,comm,ierror)
integer buffer
integer count,datatype,root,comm,ierror
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaaap1
subroutine xyzzyaaaq1(buffer,count,datatype,root,comm,ierror)
integer buffer(:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaaaq1
subroutine xyzzyaaar1(buffer,count,datatype,root,comm,ierror)
integer buffer(:,:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaaar1
subroutine xyzzyaaas1(buffer,count,datatype,root,comm,ierror)
integer buffer(:,:,:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaaas1
subroutine xyzzyaaat1(buffer,count,datatype,root,comm,ierror)
integer buffer(:,:,:,:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaaat1
subroutine xyzzyaaau1(buffer,count,datatype,root,comm,ierror)
logical buffer
integer count,datatype,root,comm,ierror
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaaau1
subroutine xyzzyaaav1(buffer,count,datatype,root,comm,ierror)
logical buffer(:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaaav1
subroutine xyzzyaaaw1(buffer,count,datatype,root,comm,ierror)
logical buffer(:,:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaaaw1
subroutine xyzzyaaax1(buffer,count,datatype,root,comm,ierror)
logical buffer(:,:,:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaaax1
subroutine xyzzyaaay1(buffer,count,datatype,root,comm,ierror)
character buffer
integer count,datatype,root,comm,ierror
if(datatype/=mpi_character)ierror=mpi_err_type
end subroutine xyzzyaaay1
subroutine xyzzyaaaz1(buffer,count,datatype,root,comm,ierror)
character buffer(:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_character)ierror=mpi_err_type
end subroutine xyzzyaaaz1
subroutine xyzzyaaba1(buffer,count,datatype,root,comm,ierror)
character buffer(:,:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_character)ierror=mpi_err_type
end subroutine xyzzyaaba1
subroutine xyzzyaabb1(buffer,count,datatype,root,comm,ierror)
real(sp) buffer
integer count,datatype,root,comm,ierror
if(datatype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaabb1
subroutine xyzzyaabc1(buffer,count,datatype,root,comm,ierror)
real(sp) buffer(:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaabc1
subroutine xyzzyaabd1(buffer,count,datatype,root,comm,ierror)
real(sp) buffer(:,:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaabd1
subroutine xyzzyaabe1(buffer,count,datatype,root,comm,ierror)
real(sp) buffer(:,:,:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaabe1
subroutine xyzzyaabf1(buffer,count,datatype,root,comm,ierror)
real(sp) buffer(:,:,:,:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaabf1
subroutine xyzzyaabg1(buffer,count,datatype,root,comm,ierror)
real(sp) buffer(:,:,:,:,:)
integer count,datatype,root,comm,ierror
if(datatype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaabg1
subroutine mpi_comm_rank(ipm_comm,ipm_rank,ipm_info)
integer ipm_comm,ipm_rank,ipm_info
ipm_info=mpi_success
ipm_rank=0
end subroutine mpi_comm_rank
subroutine mpi_comm_size(ipm_comm,ipm_node,ipm_info)
integer ipm_comm,ipm_info,ipm_node
ipm_info=mpi_success
ipm_node=1
end subroutine mpi_comm_size
subroutine mpi_error_string(ipm_error,ipm_string,ipm_len,ipm_info)
integer ipm_error,ipm_len,ipm_info
character(mpi_max_error_string) ipm_string
end subroutine mpi_error_string
subroutine mpi_finalize(ipm_info)
integer ipm_info
ipm_info=mpi_success
end subroutine mpi_finalize
subroutine xyzzyaabh1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
complex(dp) sendbuf,recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf(1)=sendbuf
if(sendtype/=mpi_double_complex)ierror=mpi_err_type
if(recvtype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaabh1
subroutine xyzzyaabi1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
complex(dp) sendbuf(:),recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf=sendbuf
if(sendtype/=mpi_double_complex)ierror=mpi_err_type
if(recvtype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaabi1
subroutine xyzzyaabj1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
complex(dp) sendbuf(:),recvbuf(:,:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf(:,1)=sendbuf(:)
if(sendtype/=mpi_double_complex)ierror=mpi_err_type
if(recvtype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaabj1
subroutine xyzzyaabk1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
complex(dp) sendbuf(:,:),recvbuf(:,:,:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf(:,:,1)=sendbuf(:,:)
if(sendtype/=mpi_double_complex)ierror=mpi_err_type
if(recvtype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaabk1
subroutine xyzzyaabl1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
complex(dp) sendbuf(:,:,:),recvbuf(:,:,:,:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf(:,:,:,1)=sendbuf(:,:,:)
if(sendtype/=mpi_double_complex)ierror=mpi_err_type
if(recvtype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaabl1
subroutine xyzzyaabm1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
complex(dp) sendbuf(:,:,:,:),recvbuf(:,:,:,:,:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf(:,:,:,:,1)=sendbuf(:,:,:,:)
if(sendtype/=mpi_double_complex)ierror=mpi_err_type
if(recvtype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaabm1
subroutine xyzzyaabn1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
real sendbuf(:),recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf=sendbuf
if(sendtype/=mpi_real)ierror=mpi_err_type
if(recvtype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaabn1
subroutine xyzzyaabo1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
real(dp) sendbuf,recvbuf
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf=sendbuf
if(sendtype/=mpi_double_precision)ierror=mpi_err_type
if(recvtype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaabo1
subroutine xyzzyaabp1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
real(dp) sendbuf,recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf(1)=sendbuf
if(sendtype/=mpi_double_precision)ierror=mpi_err_type
if(recvtype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaabp1
subroutine xyzzyaabq1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
real(dp) sendbuf(:),recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf=sendbuf
if(sendtype/=mpi_double_precision)ierror=mpi_err_type
if(recvtype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaabq1
subroutine xyzzyaabr1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
real(dp) sendbuf(:),recvbuf(:,:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf(:,1)=sendbuf(:)
if(sendtype/=mpi_double_precision)ierror=mpi_err_type
if(recvtype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaabr1
subroutine xyzzyaabs1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
real(dp) sendbuf(:,:),recvbuf(:,:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf(:,:)=sendbuf(:,:)
if(sendtype/=mpi_double_precision)ierror=mpi_err_type
if(recvtype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaabs1
subroutine xyzzyaabt1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
real(dp) sendbuf(:,:),recvbuf(:,:,:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf(:,:,1)=sendbuf(:,:)
if(sendtype/=mpi_double_precision)ierror=mpi_err_type
if(recvtype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaabt1
subroutine xyzzyaabu1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
real(dp) sendbuf(:,:,:),recvbuf(:,:,:,:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf(:,:,:,1)=sendbuf(:,:,:)
if(sendtype/=mpi_double_precision)ierror=mpi_err_type
if(recvtype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaabu1
subroutine xyzzyaabv1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
real(dp) sendbuf(:,:,:,:),recvbuf(:,:,:,:,:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf(:,:,:,:,1)=sendbuf(:,:,:,:)
if(sendtype/=mpi_double_precision)ierror=mpi_err_type
if(recvtype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaabv1
subroutine xyzzyaabw1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
real(dp) sendbuf(:,:,:,:,:),recvbuf(:,:,:,:,:,:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf(:,:,:,:,:,1)=sendbuf(:,:,:,:,:)
if(sendtype/=mpi_double_precision)ierror=mpi_err_type
if(recvtype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaabw1
subroutine xyzzyaabx1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
integer sendbuf,recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf(1)=sendbuf
if(sendtype/=mpi_integer)ierror=mpi_err_type
if(recvtype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaabx1
subroutine xyzzyaaby1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
integer(i64) sendbuf,recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf(1)=sendbuf
if(sendtype/=mpi_integer)ierror=mpi_err_type
if(recvtype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaaby1
subroutine xyzzyaabz1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
integer sendbuf(:),recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf=sendbuf
if(sendtype/=mpi_integer)ierror=mpi_err_type
if(recvtype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaabz1
subroutine xyzzyaaca1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
integer sendbuf(:),recvbuf(:,:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf(:,1)=sendbuf(:)
if(sendtype/=mpi_integer)ierror=mpi_err_type
if(recvtype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaaca1
subroutine xyzzyaacb1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
logical sendbuf,recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf(1)=sendbuf
if(sendtype/=mpi_logical)ierror=mpi_err_type
if(recvtype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaacb1
subroutine xyzzyaacc1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
logical sendbuf(:),recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf=sendbuf
if(sendtype/=mpi_logical)ierror=mpi_err_type
if(recvtype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaacc1
subroutine xyzzyaacd1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
logical sendbuf(:),recvbuf(:,:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf(:,1)=sendbuf(:)
if(sendtype/=mpi_logical)ierror=mpi_err_type
if(recvtype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaacd1
subroutine xyzzyaace1(sendcount,sendtype,recvbuf,recvcount,recvtype,ro&
&ot,comm,ierror)
real(dp) recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
if(sendtype/=mpi_double_precision)ierror=mpi_err_type
if(recvtype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaace1
subroutine xyzzyaacf1(sendcount,sendtype,recvbuf,recvcount,recvtype,ro&
&ot,comm,ierror)
integer recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
if(sendtype/=mpi_integer)ierror=mpi_err_type
if(recvtype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaacf1
subroutine mpi_init(ipm_info)
integer ipm_info
ipm_info=mpi_success
end subroutine mpi_init
subroutine xyzzyaacg1(buf,count,datatype,source,tag,comm,status,ierror&
&)
complex(dp) buf
integer count,datatype,source,tag,comm,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaacg1
subroutine xyzzyaach1(buf,count,datatype,source,tag,comm,status,ierror&
&)
complex(dp) buf(:)
integer count,datatype,source,tag,comm,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaach1
subroutine xyzzyaaci1(buf,count,datatype,source,tag,comm,status,ierror&
&)
complex(dp) buf(:,:)
integer count,datatype,source,tag,comm,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaaci1
subroutine xyzzyaacj1(buf,count,datatype,source,tag,comm,status,ierror&
&)
complex(dp) buf(:,:,:)
integer count,datatype,source,tag,comm,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaacj1
subroutine xyzzyaack1(buf,count,datatype,source,tag,comm,status,ierror&
&)
complex(dp) buf(:,:,:,:)
integer count,datatype,source,tag,comm,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaack1
subroutine xyzzyaacl1(buf,count,datatype,source,tag,comm,status,ierror&
&)
complex(dp) buf(:,:,:,:,:)
integer count,datatype,source,tag,comm,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaacl1
subroutine xyzzyaacm1(buf,count,datatype,source,tag,comm,status,ierror&
&)
real buf(:)
integer count,datatype,source,tag,comm,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaacm1
subroutine xyzzyaacn1(buf,count,datatype,source,tag,comm,status,ierror&
&)
real(dp) buf
integer count,datatype,source,tag,comm,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaacn1
subroutine xyzzyaaco1(buf,count,datatype,source,tag,comm,status,ierror&
&)
real(dp) buf(:)
integer count,datatype,source,tag,comm,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaaco1
subroutine xyzzyaacp1(buf,count,datatype,source,tag,comm,status,ierror&
&)
real(dp) buf(:,:)
integer count,datatype,source,tag,comm,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaacp1
subroutine xyzzyaacq1(buf,count,datatype,source,tag,comm,status,ierror&
&)
real(dp) buf(:,:,:)
integer count,datatype,source,tag,comm,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaacq1
subroutine xyzzyaacr1(buf,count,datatype,source,tag,comm,status,ierror&
&)
real(dp) buf(:,:,:,:)
integer count,datatype,source,tag,comm,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaacr1
subroutine xyzzyaacs1(buf,count,datatype,source,tag,comm,status,ierror&
&)
real(dp) buf(:,:,:,:,:)
integer count,datatype,source,tag,comm,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaacs1
subroutine xyzzyaact1(buf,count,datatype,source,tag,comm,status,ierror&
&)
real(dp) buf(:,:,:,:,:,:)
integer count,datatype,source,tag,comm,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaact1
subroutine xyzzyaacu1(buf,count,datatype,source,tag,comm,status,ierror&
&)
integer buf
integer count,datatype,source,tag,comm,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaacu1
subroutine xyzzyaacv1(buf,count,datatype,source,tag,comm,status,ierror&
&)
integer buf(:)
integer count,datatype,source,tag,comm,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaacv1
subroutine xyzzyaacw1(buf,count,datatype,source,tag,comm,status,ierror&
&)
integer buf(:)
integer count,datatype,source,tag,comm,status(:,:),ierror
ierror=mpi_success
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaacw1
subroutine xyzzyaacx1(buf,count,datatype,source,tag,comm,status,ierror&
&)
integer buf(:,:)
integer count,datatype,source,tag,comm,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaacx1
subroutine xyzzyaacy1(buf,count,datatype,source,tag,comm,status,ierror&
&)
integer buf(:,:,:)
integer count,datatype,source,tag,comm,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaacy1
subroutine xyzzyaacz1(buf,count,datatype,source,tag,comm,status,ierror&
&)
logical buf
integer count,datatype,source,tag,comm,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaacz1
subroutine xyzzyaada1(buf,count,datatype,source,tag,comm,status,ierror&
&)
logical buf(:)
integer count,datatype,source,tag,comm,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaada1
subroutine xyzzyaadb1(buf,count,datatype,source,tag,comm,status,ierror&
&)
logical buf(:,:)
integer count,datatype,source,tag,comm,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaadb1
subroutine xyzzyaadc1(buf,count,datatype,source,tag,comm,request,ierro&
&r)
complex(dp) buf
integer count,datatype,source,tag,comm,request,ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaadc1
subroutine xyzzyaadd1(buf,count,datatype,source,tag,comm,request,ierro&
&r)
complex(dp) buf(:)
integer count,datatype,source,tag,comm,request,ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaadd1
subroutine xyzzyaade1(buf,count,datatype,source,tag,comm,request,ierro&
&r)
complex(dp) buf(:,:)
integer count,datatype,source,tag,comm,request,ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaade1
subroutine xyzzyaadf1(buf,count,datatype,source,tag,comm,request,ierro&
&r)
complex(dp) buf(:,:,:)
integer count,datatype,source,tag,comm,request,ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaadf1
subroutine xyzzyaadg1(buf,count,datatype,source,tag,comm,request,ierro&
&r)
complex(dp) buf(:,:,:,:)
integer count,datatype,source,tag,comm,request,ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaadg1
subroutine xyzzyaadh1(buf,count,datatype,source,tag,comm,request,ierro&
&r)
complex(dp) buf(:,:,:,:,:)
integer count,datatype,source,tag,comm,request,ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaadh1
subroutine xyzzyaadi1(buf,count,datatype,source,tag,comm,request,ierro&
&r)
real(dp) buf
integer count,datatype,source,tag,comm,request,ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaadi1
subroutine xyzzyaadj1(buf,count,datatype,source,tag,comm,request,ierro&
&r)
real(dp) buf(:)
integer count,datatype,source,tag,comm,request,ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaadj1
subroutine xyzzyaadk1(buf,count,datatype,source,tag,comm,request,ierro&
&r)
real(dp) buf(:,:)
integer count,datatype,source,tag,comm,request,ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaadk1
subroutine xyzzyaadl1(buf,count,datatype,source,tag,comm,request,ierro&
&r)
real(dp) buf(:,:,:)
integer count,datatype,source,tag,comm,request,ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaadl1
subroutine xyzzyaadm1(buf,count,datatype,source,tag,comm,request,ierro&
&r)
real(dp) buf(:,:,:,:)
integer count,datatype,source,tag,comm,request,ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaadm1
subroutine xyzzyaadn1(buf,count,datatype,source,tag,comm,request,ierro&
&r)
real(dp) buf(:,:,:,:,:)
integer count,datatype,source,tag,comm,request,ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaadn1
subroutine xyzzyaado1(buf,count,datatype,source,tag,comm,request,ierro&
&r)
real(dp) buf(:,:,:,:,:,:)
integer count,datatype,source,tag,comm,request,ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaado1
subroutine xyzzyaadp1(buf,count,datatype,source,tag,comm,request,ierro&
&r)
integer buf
integer count,datatype,source,tag,comm,request,ierror
ierror=mpi_success
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaadp1
subroutine xyzzyaadq1(buf,count,datatype,source,tag,comm,request,ierro&
&r)
integer buf(:)
integer count,datatype,source,tag,comm,request,ierror
ierror=mpi_success
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaadq1
subroutine xyzzyaadr1(buf,count,datatype,source,tag,comm,request,ierro&
&r)
integer buf(:,:)
integer count,datatype,source,tag,comm,request,ierror
ierror=mpi_success
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaadr1
subroutine xyzzyaads1(buf,count,datatype,source,tag,comm,request,ierro&
&r)
integer buf(:,:,:)
integer count,datatype,source,tag,comm,request,ierror
ierror=mpi_success
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaads1
subroutine xyzzyaadt1(buf,count,datatype,source,tag,comm,request,ierro&
&r)
logical buf
integer count,datatype,source,tag,comm,request,ierror
ierror=mpi_success
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaadt1
subroutine xyzzyaadu1(buf,count,datatype,source,tag,comm,request,ierro&
&r)
logical buf(:)
integer count,datatype,source,tag,comm,request,ierror
ierror=mpi_success
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaadu1
subroutine xyzzyaadv1(buf,count,datatype,source,tag,comm,request,ierro&
&r)
logical buf(:,:)
integer count,datatype,source,tag,comm,request,ierror
ierror=mpi_success
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaadv1
subroutine mpi_file_open(comm,filename,amode,info,fh,ierror)
integer comm,amode,info,fh,ierror
character(*) filename
end subroutine mpi_file_open
subroutine mpi_file_close(fh,ierror)
integer fh,ierror
ierror=mpi_success
end subroutine mpi_file_close
subroutine mpi_file_get_byte_offset(fh,offset,disp,ierror)
integer fh,ierror
integer(mpi_offset_kind) offset,disp
ierror=mpi_success
end subroutine mpi_file_get_byte_offset
subroutine mpi_file_get_position(fh,offset,ierror)
integer fh,ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
end subroutine mpi_file_get_position
subroutine xyzzyaadw1(fh,buf,count,datatype,status,ierror)
complex(dp) buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaadw1
subroutine xyzzyaadx1(fh,buf,count,datatype,status,ierror)
complex(dp) buf(:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaadx1
subroutine xyzzyaady1(fh,buf,count,datatype,status,ierror)
complex(dp) buf(:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaady1
subroutine xyzzyaadz1(fh,buf,count,datatype,status,ierror)
complex(dp) buf(:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaadz1
subroutine xyzzyaaea1(fh,buf,count,datatype,status,ierror)
complex(dp) buf(:,:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaaea1
subroutine xyzzyaaeb1(fh,buf,count,datatype,status,ierror)
complex(dp) buf(:,:,:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaaeb1
subroutine xyzzyaaec1(fh,buf,count,datatype,status,ierror)
complex(sp) buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_complex)ierror=mpi_err_type
end subroutine xyzzyaaec1
subroutine xyzzyaaed1(fh,buf,count,datatype,status,ierror)
real(dp) buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaaed1
subroutine xyzzyaaee1(fh,buf,count,datatype,status,ierror)
real(dp) buf(:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaaee1
subroutine xyzzyaaef1(fh,buf,count,datatype,status,ierror)
real(dp) buf(:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaaef1
subroutine xyzzyaaeg1(fh,buf,count,datatype,status,ierror)
real(dp) buf(:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaaeg1
subroutine xyzzyaaeh1(fh,buf,count,datatype,status,ierror)
real(dp) buf(:,:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaaeh1
subroutine xyzzyaaei1(fh,buf,count,datatype,status,ierror)
real(dp) buf(:,:,:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaaei1
subroutine xyzzyaaej1(fh,buf,count,datatype,status,ierror)
real(dp) buf(:,:,:,:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaaej1
subroutine xyzzyaaek1(fh,buf,count,datatype,status,ierror)
real buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaaek1
subroutine xyzzyaael1(fh,buf,count,datatype,status,ierror)
real buf(:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaael1
subroutine xyzzyaaem1(fh,buf,count,datatype,status,ierror)
real buf(:,:,:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaaem1
subroutine xyzzyaaen1(fh,buf,count,datatype,status,ierror)
integer buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaaen1
subroutine xyzzyaaeo1(fh,buf,count,datatype,status,ierror)
integer buf(:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaaeo1
subroutine xyzzyaaep1(fh,buf,count,datatype,status,ierror)
integer buf(:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaaep1
subroutine xyzzyaaeq1(fh,buf,count,datatype,status,ierror)
logical buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaaeq1
subroutine xyzzyaaer1(fh,buf,count,datatype,status,ierror)
logical buf(:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaaer1
subroutine xyzzyaaes1(fh,buf,count,datatype,status,ierror)
logical buf(:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaaes1
subroutine xyzzyaaet1(fh,buf,count,datatype,status,ierror)
character(*) buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_character)ierror=mpi_err_type
end subroutine xyzzyaaet1
subroutine xyzzyaaeu1(fh,offset,buf,count,datatype,status,ierror)
complex(dp) buf
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaaeu1
subroutine xyzzyaaev1(fh,offset,buf,count,datatype,status,ierror)
complex(dp) buf(:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaaev1
subroutine xyzzyaaew1(fh,offset,buf,count,datatype,status,ierror)
complex(dp) buf(:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaaew1
subroutine xyzzyaaex1(fh,offset,buf,count,datatype,status,ierror)
complex(dp) buf(:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaaex1
subroutine xyzzyaaey1(fh,offset,buf,count,datatype,status,ierror)
complex(dp) buf(:,:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaaey1
subroutine xyzzyaaez1(fh,offset,buf,count,datatype,status,ierror)
complex(dp) buf(:,:,:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaaez1
subroutine xyzzyaafa1(fh,offset,buf,count,datatype,status,ierror)
real(dp) buf
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaafa1
subroutine xyzzyaafb1(fh,offset,buf,count,datatype,status,ierror)
real(dp) buf(:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaafb1
subroutine xyzzyaafc1(fh,offset,buf,count,datatype,status,ierror)
real(dp) buf(:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaafc1
subroutine xyzzyaafd1(fh,offset,buf,count,datatype,status,ierror)
real(dp) buf(:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaafd1
subroutine xyzzyaafe1(fh,offset,buf,count,datatype,status,ierror)
real(dp) buf(:,:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaafe1
subroutine xyzzyaaff1(fh,offset,buf,count,datatype,status,ierror)
real(dp) buf(:,:,:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaaff1
subroutine xyzzyaafg1(fh,offset,buf,count,datatype,status,ierror)
real(dp) buf(:,:,:,:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaafg1
subroutine xyzzyaafh1(fh,offset,buf,count,datatype,status,ierror)
real buf
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaafh1
subroutine xyzzyaafi1(fh,offset,buf,count,datatype,status,ierror)
real buf(:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaafi1
subroutine xyzzyaafj1(fh,offset,buf,count,datatype,status,ierror)
integer buf
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaafj1
subroutine xyzzyaafk1(fh,offset,buf,count,datatype,status,ierror)
integer buf(:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaafk1
subroutine xyzzyaafl1(fh,offset,buf,count,datatype,status,ierror)
integer buf(:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaafl1
subroutine xyzzyaafm1(fh,offset,buf,count,datatype,status,ierror)
logical buf
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaafm1
subroutine xyzzyaafn1(fh,offset,buf,count,datatype,status,ierror)
logical buf(:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaafn1
subroutine xyzzyaafo1(fh,offset,buf,count,datatype,status,ierror)
logical buf(:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaafo1
subroutine xyzzyaafp1(fh,offset,buf,count,datatype,status,ierror)
character(*) buf
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_character)ierror=mpi_err_type
end subroutine xyzzyaafp1
subroutine xyzzyaafq1(fh,buf,count,datatype,status,ierror)
complex(dp) buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaafq1
subroutine xyzzyaafr1(fh,buf,count,datatype,status,ierror)
complex(sp) buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_complex)ierror=mpi_err_type
end subroutine xyzzyaafr1
subroutine xyzzyaafs1(fh,buf,count,datatype,status,ierror)
real(dp) buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaafs1
subroutine xyzzyaaft1(fh,buf,count,datatype,status,ierror)
real(sp) buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaaft1
subroutine mpi_file_seek(fh,offset,whence,ierror)
integer fh,whence,ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
end subroutine mpi_file_seek
subroutine mpi_file_set_view(fh,disp,etype,filetype,datarep,info,ierro&
&r)
integer fh,etype,filetype,info,ierror
integer(mpi_offset_kind) disp
character(*) datarep
ierror=mpi_success
end subroutine mpi_file_set_view
subroutine xyzzyaafu1(fh,buf,count,datatype,status,ierror)
complex(dp) buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaafu1
subroutine xyzzyaafv1(fh,buf,count,datatype,status,ierror)
complex(dp) buf(:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaafv1
subroutine xyzzyaafw1(fh,buf,count,datatype,status,ierror)
complex(dp) buf(:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaafw1
subroutine xyzzyaafx1(fh,buf,count,datatype,status,ierror)
complex(dp) buf(:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaafx1
subroutine xyzzyaafy1(fh,buf,count,datatype,status,ierror)
complex(dp) buf(:,:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaafy1
subroutine xyzzyaafz1(fh,buf,count,datatype,status,ierror)
complex(dp) buf(:,:,:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaafz1
subroutine xyzzyaaga1(fh,buf,count,datatype,status,ierror)
complex(sp) buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_complex)ierror=mpi_err_type
end subroutine xyzzyaaga1
subroutine xyzzyaagb1(fh,buf,count,datatype,status,ierror)
real(dp) buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaagb1
subroutine xyzzyaagc1(fh,buf,count,datatype,status,ierror)
real(dp) buf(:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaagc1
subroutine xyzzyaagd1(fh,buf,count,datatype,status,ierror)
real(dp) buf(:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaagd1
subroutine xyzzyaage1(fh,buf,count,datatype,status,ierror)
real(dp) buf(:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaage1
subroutine xyzzyaagf1(fh,buf,count,datatype,status,ierror)
real(dp) buf(:,:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaagf1
subroutine xyzzyaagg1(fh,buf,count,datatype,status,ierror)
real(dp) buf(:,:,:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaagg1
subroutine xyzzyaagh1(fh,buf,count,datatype,status,ierror)
real(dp) buf(:,:,:,:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaagh1
subroutine xyzzyaagi1(fh,buf,count,datatype,status,ierror)
real buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaagi1
subroutine xyzzyaagj1(fh,buf,count,datatype,status,ierror)
real buf(:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaagj1
subroutine xyzzyaagk1(fh,buf,count,datatype,status,ierror)
integer buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaagk1
subroutine xyzzyaagl1(fh,buf,count,datatype,status,ierror)
integer buf(:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaagl1
subroutine xyzzyaagm1(fh,buf,count,datatype,status,ierror)
integer buf(:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaagm1
subroutine xyzzyaagn1(fh,buf,count,datatype,status,ierror)
logical buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaagn1
subroutine xyzzyaago1(fh,buf,count,datatype,status,ierror)
logical buf(:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaago1
subroutine xyzzyaagp1(fh,buf,count,datatype,status,ierror)
logical buf(:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaagp1
subroutine xyzzyaagq1(fh,buf,count,datatype,status,ierror)
character(*) buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_character)ierror=mpi_err_type
end subroutine xyzzyaagq1
subroutine xyzzyaagr1(fh,offset,buf,count,datatype,status,ierror)
complex(dp) buf
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaagr1
subroutine xyzzyaags1(fh,offset,buf,count,datatype,status,ierror)
complex(dp) buf(:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaags1
subroutine xyzzyaagt1(fh,offset,buf,count,datatype,status,ierror)
complex(dp) buf(:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaagt1
subroutine xyzzyaagu1(fh,offset,buf,count,datatype,status,ierror)
complex(dp) buf(:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaagu1
subroutine xyzzyaagv1(fh,offset,buf,count,datatype,status,ierror)
complex(dp) buf(:,:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaagv1
subroutine xyzzyaagw1(fh,offset,buf,count,datatype,status,ierror)
complex(dp) buf(:,:,:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaagw1
subroutine xyzzyaagx1(fh,offset,buf,count,datatype,status,ierror)
real(dp) buf
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaagx1
subroutine xyzzyaagy1(fh,offset,buf,count,datatype,status,ierror)
real(dp) buf(:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaagy1
subroutine xyzzyaagz1(fh,offset,buf,count,datatype,status,ierror)
real(dp) buf(:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaagz1
subroutine xyzzyaaha1(fh,offset,buf,count,datatype,status,ierror)
real(dp) buf(:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaaha1
subroutine xyzzyaahb1(fh,offset,buf,count,datatype,status,ierror)
real(dp) buf(:,:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaahb1
subroutine xyzzyaahc1(fh,offset,buf,count,datatype,status,ierror)
real(dp) buf(:,:,:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaahc1
subroutine xyzzyaahd1(fh,offset,buf,count,datatype,status,ierror)
real(dp) buf(:,:,:,:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaahd1
subroutine xyzzyaahe1(fh,offset,buf,count,datatype,status,ierror)
real buf
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaahe1
subroutine xyzzyaahf1(fh,offset,buf,count,datatype,status,ierror)
real buf(:,:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaahf1
subroutine xyzzyaahg1(fh,offset,buf,count,datatype,status,ierror)
integer buf
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaahg1
subroutine xyzzyaahh1(fh,offset,buf,count,datatype,status,ierror)
integer buf(:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaahh1
subroutine xyzzyaahi1(fh,offset,buf,count,datatype,status,ierror)
integer buf(:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaahi1
subroutine xyzzyaahj1(fh,offset,buf,count,datatype,status,ierror)
logical buf
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaahj1
subroutine xyzzyaahk1(fh,offset,buf,count,datatype,status,ierror)
logical buf(:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaahk1
subroutine xyzzyaahl1(fh,offset,buf,count,datatype,status,ierror)
logical buf(:,:)
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaahl1
subroutine xyzzyaahm1(fh,offset,buf,count,datatype,status,ierror)
character(*) buf
integer fh,count,datatype,status(mpi_status_size),ierror
integer(mpi_offset_kind) offset
ierror=mpi_success
if(datatype/=mpi_character)ierror=mpi_err_type
end subroutine xyzzyaahm1
subroutine xyzzyaahn1(fh,buf,count,datatype,status,ierror)
complex(dp) buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaahn1
subroutine xyzzyaaho1(fh,buf,count,datatype,status,ierror)
complex(sp) buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_complex)ierror=mpi_err_type
end subroutine xyzzyaaho1
subroutine xyzzyaahp1(fh,buf,count,datatype,status,ierror)
real(dp) buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaahp1
subroutine xyzzyaahq1(fh,buf,count,datatype,status,ierror)
real(sp) buf
integer fh,count,datatype,status(mpi_status_size),ierror
ierror=mpi_success
if(datatype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaahq1
subroutine xyzzyaahr1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_root,ipm_comm,ipm_info)
complex(dp) ipm_sendbuf,ipm_recvbuf
integer ipm_count,ipm_datatype,ipm_op,ipm_root,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_double_complex)ipm_info=mpi_err_type
end subroutine xyzzyaahr1
subroutine xyzzyaahs1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_root,ipm_comm,ipm_info)
complex(dp) ipm_sendbuf(:),ipm_recvbuf(:)
integer ipm_count,ipm_datatype,ipm_op,ipm_root,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_double_complex)ipm_info=mpi_err_type
end subroutine xyzzyaahs1
subroutine xyzzyaaht1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_root,ipm_comm,ipm_info)
complex(dp) ipm_sendbuf(:,:),ipm_recvbuf(:,:)
integer ipm_count,ipm_datatype,ipm_op,ipm_root,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_double_complex)ipm_info=mpi_err_type
end subroutine xyzzyaaht1
subroutine xyzzyaahu1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_root,ipm_comm,ipm_info)
complex(dp) ipm_sendbuf(:,:,:),ipm_recvbuf(:,:,:)
integer ipm_count,ipm_datatype,ipm_op,ipm_root,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_double_complex)ipm_info=mpi_err_type
end subroutine xyzzyaahu1
subroutine xyzzyaahv1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_root,ipm_comm,ipm_info)
real(dp) ipm_sendbuf,ipm_recvbuf
integer ipm_count,ipm_datatype,ipm_op,ipm_root,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_double_precision)ipm_info=mpi_err_type
end subroutine xyzzyaahv1
subroutine xyzzyaahw1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_root,ipm_comm,ipm_info)
real(dp) ipm_sendbuf(:),ipm_recvbuf(:)
integer ipm_count,ipm_datatype,ipm_op,ipm_root,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_double_precision)ipm_info=mpi_err_type
end subroutine xyzzyaahw1
subroutine xyzzyaahx1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_root,ipm_comm,ipm_info)
real(dp) ipm_sendbuf(:,:),ipm_recvbuf(:,:)
integer ipm_count,ipm_datatype,ipm_op,ipm_root,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_double_precision)ipm_info=mpi_err_type
end subroutine xyzzyaahx1
subroutine xyzzyaahy1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_root,ipm_comm,ipm_info)
real(dp) ipm_sendbuf(:,:,:),ipm_recvbuf(:,:,:)
integer ipm_count,ipm_datatype,ipm_op,ipm_root,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_double_precision)ipm_info=mpi_err_type
end subroutine xyzzyaahy1
subroutine xyzzyaahz1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_root,ipm_comm,ipm_info)
real(dp) ipm_sendbuf(:,:,:,:),ipm_recvbuf(:,:,:,:)
integer ipm_count,ipm_datatype,ipm_op,ipm_root,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_double_precision)ipm_info=mpi_err_type
end subroutine xyzzyaahz1
subroutine xyzzyaaia1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_root,ipm_comm,ipm_info)
integer ipm_sendbuf,ipm_recvbuf
integer ipm_count,ipm_datatype,ipm_op,ipm_root,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_integer)ipm_info=mpi_err_type
end subroutine xyzzyaaia1
subroutine xyzzyaaib1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_root,ipm_comm,ipm_info)
integer(i64) ipm_sendbuf,ipm_recvbuf
integer ipm_count,ipm_datatype,ipm_op,ipm_root,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_integer8)ipm_info=mpi_err_type
end subroutine xyzzyaaib1
subroutine xyzzyaaic1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_root,ipm_comm,ipm_info)
integer ipm_sendbuf(:),ipm_recvbuf(:)
integer ipm_count,ipm_datatype,ipm_op,ipm_root,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_integer)ipm_info=mpi_err_type
end subroutine xyzzyaaic1
subroutine xyzzyaaid1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_root,ipm_comm,ipm_info)
integer ipm_sendbuf(:,:),ipm_recvbuf(:,:)
integer ipm_count,ipm_datatype,ipm_op,ipm_root,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_integer)ipm_info=mpi_err_type
end subroutine xyzzyaaid1
subroutine xyzzyaaie1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_root,ipm_comm,ipm_info)
integer ipm_sendbuf(:,:,:),ipm_recvbuf(:,:,:)
integer ipm_count,ipm_datatype,ipm_op,ipm_root,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_integer)ipm_info=mpi_err_type
end subroutine xyzzyaaie1
subroutine xyzzyaaif1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_root,ipm_comm,ipm_info)
logical ipm_sendbuf,ipm_recvbuf
integer ipm_count,ipm_datatype,ipm_op,ipm_root,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_logical)ipm_info=mpi_err_type
end subroutine xyzzyaaif1
subroutine xyzzyaaig1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_root,ipm_comm,ipm_info)
logical ipm_sendbuf(:,:,:),ipm_recvbuf(:,:,:)
integer ipm_count,ipm_datatype,ipm_op,ipm_root,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_logical)ipm_info=mpi_err_type
end subroutine xyzzyaaig1
subroutine xyzzyaaih1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_comm,ipm_info)
real(dp) ipm_sendbuf,ipm_recvbuf
integer ipm_count,ipm_datatype,ipm_op,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_double_precision)ipm_info=mpi_err_type
end subroutine xyzzyaaih1
subroutine xyzzyaaii1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_comm,ipm_info)
real(dp) ipm_sendbuf(:),ipm_recvbuf(:)
integer ipm_count,ipm_datatype,ipm_op,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_double_precision)ipm_info=mpi_err_type
end subroutine xyzzyaaii1
subroutine xyzzyaaij1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_comm,ipm_info)
real(dp) ipm_sendbuf(:,:),ipm_recvbuf(:,:)
integer ipm_count,ipm_datatype,ipm_op,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_double_precision)ipm_info=mpi_err_type
end subroutine xyzzyaaij1
subroutine xyzzyaaik1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_comm,ipm_info)
real(dp) ipm_sendbuf(:,:,:),ipm_recvbuf(:,:,:)
integer ipm_count,ipm_datatype,ipm_op,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_double_precision)ipm_info=mpi_err_type
end subroutine xyzzyaaik1
subroutine xyzzyaail1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_comm,ipm_info)
real(dp) ipm_sendbuf(:,:,:,:),ipm_recvbuf(:,:,:,:)
integer ipm_count,ipm_datatype,ipm_op,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_double_precision)ipm_info=mpi_err_type
end subroutine xyzzyaail1
subroutine xyzzyaaim1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_comm,ipm_info)
logical ipm_sendbuf(:),ipm_recvbuf(:)
integer ipm_count,ipm_datatype,ipm_op,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_logical)ipm_info=mpi_err_type
end subroutine xyzzyaaim1
subroutine xyzzyaain1(ipm_sendbuf,ipm_recvbuf,ipm_count,ipm_datatype,i&
&pm_op,ipm_comm,ipm_info)
integer ipm_sendbuf,ipm_recvbuf
integer ipm_count,ipm_datatype,ipm_op,ipm_comm,ipm_info
ipm_recvbuf=ipm_sendbuf
ipm_info=mpi_success
if(ipm_datatype/=mpi_integer)ipm_info=mpi_err_type
end subroutine xyzzyaain1
subroutine xyzzyaaio1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
complex(dp) sendbuf,recvbuf
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf=sendbuf
if(sendtype/=mpi_double_complex)ierror=mpi_err_type
if(recvtype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaaio1
subroutine xyzzyaaip1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
complex(dp) sendbuf(:),recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf=sendbuf
if(sendtype/=mpi_double_complex)ierror=mpi_err_type
if(recvtype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaaip1
subroutine xyzzyaaiq1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
complex(dp) sendbuf(:,:),recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf(:)=sendbuf(:,1)
if(sendtype/=mpi_double_complex)ierror=mpi_err_type
if(recvtype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaaiq1
subroutine xyzzyaair1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
real(dp) sendbuf,recvbuf
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf=sendbuf
if(sendtype/=mpi_double_precision)ierror=mpi_err_type
if(recvtype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaair1
subroutine xyzzyaais1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
real(dp) sendbuf(:),recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf=sendbuf
if(sendtype/=mpi_double_precision)ierror=mpi_err_type
if(recvtype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaais1
subroutine xyzzyaait1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
real(dp) sendbuf(:,:),recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf(:)=sendbuf(:,1)
if(sendtype/=mpi_double_precision)ierror=mpi_err_type
if(recvtype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaait1
subroutine xyzzyaaiu1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
integer sendbuf,recvbuf
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf=sendbuf
if(sendtype/=mpi_integer)ierror=mpi_err_type
if(recvtype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaaiu1
subroutine xyzzyaaiv1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
integer sendbuf(:),recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf=sendbuf
if(sendtype/=mpi_integer)ierror=mpi_err_type
if(recvtype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaaiv1
subroutine xyzzyaaiw1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
integer sendbuf(:,:),recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf(:)=sendbuf(:,1)
if(sendtype/=mpi_integer)ierror=mpi_err_type
if(recvtype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaaiw1
subroutine xyzzyaaix1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
logical sendbuf,recvbuf
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf=sendbuf
if(sendtype/=mpi_logical)ierror=mpi_err_type
if(recvtype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaaix1
subroutine xyzzyaaiy1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
logical sendbuf(:),recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf=sendbuf
if(sendtype/=mpi_logical)ierror=mpi_err_type
if(recvtype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaaiy1
subroutine xyzzyaaiz1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,root,comm,ierror)
logical sendbuf(:,:),recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
recvbuf(:)=sendbuf(:,1)
if(sendtype/=mpi_logical)ierror=mpi_err_type
if(recvtype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaaiz1
subroutine xyzzyaaja1(buf,count,datatype,dest,tag,comm,ierror)
complex(dp) buf
integer count,datatype,dest,tag,comm,ierror
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaaja1
subroutine xyzzyaajb1(buf,count,datatype,dest,tag,comm,ierror)
complex(dp) buf(:)
integer count,datatype,dest,tag,comm,ierror
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaajb1
subroutine xyzzyaajc1(buf,count,datatype,dest,tag,comm,ierror)
complex(dp) buf(:,:)
integer count,datatype,dest,tag,comm,ierror
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaajc1
subroutine xyzzyaajd1(buf,count,datatype,dest,tag,comm,ierror)
complex(dp) buf(:,:,:)
integer count,datatype,dest,tag,comm,ierror
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaajd1
subroutine xyzzyaaje1(buf,count,datatype,dest,tag,comm,ierror)
real buf(:)
integer count,datatype,dest,tag,comm,ierror
if(datatype/=mpi_real)ierror=mpi_err_type
end subroutine xyzzyaaje1
subroutine xyzzyaajf1(buf,count,datatype,dest,tag,comm,ierror)
real(dp) buf
integer count,datatype,dest,tag,comm,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaajf1
subroutine xyzzyaajg1(buf,count,datatype,dest,tag,comm,ierror)
real(dp) buf(:)
integer count,datatype,dest,tag,comm,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaajg1
subroutine xyzzyaajh1(buf,count,datatype,dest,tag,comm,ierror)
real(dp) buf(:,:)
integer count,datatype,dest,tag,comm,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaajh1
subroutine xyzzyaaji1(buf,count,datatype,dest,tag,comm,ierror)
real(dp) buf(:,:,:)
integer count,datatype,dest,tag,comm,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaaji1
subroutine xyzzyaajj1(buf,count,datatype,dest,tag,comm,ierror)
real(dp) buf(:,:,:,:)
integer count,datatype,dest,tag,comm,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaajj1
subroutine xyzzyaajk1(buf,count,datatype,dest,tag,comm,ierror)
real(dp) buf(:,:,:,:,:)
integer count,datatype,dest,tag,comm,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaajk1
subroutine xyzzyaajl1(buf,count,datatype,dest,tag,comm,ierror)
real(dp) buf(:,:,:,:,:,:)
integer count,datatype,dest,tag,comm,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaajl1
subroutine xyzzyaajm1(buf,count,datatype,dest,tag,comm,ierror)
integer buf
integer count,datatype,dest,tag,comm,ierror
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaajm1
subroutine xyzzyaajn1(buf,count,datatype,dest,tag,comm,ierror)
integer buf(:)
integer count,datatype,dest,tag,comm,ierror
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaajn1
subroutine xyzzyaajo1(buf,count,datatype,dest,tag,comm,ierror)
integer buf(:,:)
integer count,datatype,dest,tag,comm,ierror
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaajo1
subroutine xyzzyaajp1(buf,count,datatype,dest,tag,comm,ierror)
integer buf(:,:,:)
integer count,datatype,dest,tag,comm,ierror
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaajp1
subroutine xyzzyaajq1(buf,count,datatype,dest,tag,comm,ierror)
logical buf
integer count,datatype,dest,tag,comm,ierror
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaajq1
subroutine xyzzyaajr1(buf,count,datatype,dest,tag,comm,ierror)
logical buf(:)
integer count,datatype,dest,tag,comm,ierror
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaajr1
subroutine xyzzyaajs1(buf,count,datatype,dest,tag,comm,ierror)
logical buf(:,:)
integer count,datatype,dest,tag,comm,ierror
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaajs1
subroutine xyzzyaajt1(buf,count,datatype,dest,tag,comm,request,ierror)
complex(dp) buf
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaajt1
subroutine xyzzyaaju1(buf,count,datatype,dest,tag,comm,request,ierror)
complex(dp) buf(:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaaju1
subroutine xyzzyaajv1(buf,count,datatype,dest,tag,comm,request,ierror)
complex(dp) buf(:,:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaajv1
subroutine xyzzyaajw1(buf,count,datatype,dest,tag,comm,request,ierror)
complex(dp) buf(:,:,:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaajw1
subroutine xyzzyaajx1(buf,count,datatype,dest,tag,comm,request,ierror)
real(dp) buf
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaajx1
subroutine xyzzyaajy1(buf,count,datatype,dest,tag,comm,request,ierror)
real(dp) buf(:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaajy1
subroutine xyzzyaajz1(buf,count,datatype,dest,tag,comm,request,ierror)
real(dp) buf(:,:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaajz1
subroutine xyzzyaaka1(buf,count,datatype,dest,tag,comm,request,ierror)
real(dp) buf(:,:,:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaaka1
subroutine xyzzyaakb1(buf,count,datatype,dest,tag,comm,request,ierror)
real(dp) buf(:,:,:,:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaakb1
subroutine xyzzyaakc1(buf,count,datatype,dest,tag,comm,request,ierror)
real(dp) buf(:,:,:,:,:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaakc1
subroutine xyzzyaakd1(buf,count,datatype,dest,tag,comm,request,ierror)
real(dp) buf(:,:,:,:,:,:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaakd1
subroutine xyzzyaake1(buf,count,datatype,dest,tag,comm,request,ierror)
integer buf
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaake1
subroutine xyzzyaakf1(buf,count,datatype,dest,tag,comm,request,ierror)
integer buf(:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaakf1
subroutine xyzzyaakg1(buf,count,datatype,dest,tag,comm,request,ierror)
integer buf(:,:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaakg1
subroutine xyzzyaakh1(buf,count,datatype,dest,tag,comm,request,ierror)
integer buf(:,:,:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaakh1
subroutine xyzzyaaki1(buf,count,datatype,dest,tag,comm,request,ierror)
logical buf
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaaki1
subroutine xyzzyaakj1(buf,count,datatype,dest,tag,comm,request,ierror)
logical buf(:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaakj1
subroutine xyzzyaakk1(buf,count,datatype,dest,tag,comm,request,ierror)
logical buf(:,:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaakk1
subroutine xyzzyaakl1(buf,count,datatype,dest,tag,comm,request,ierror)
complex(dp) buf
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaakl1
subroutine xyzzyaakm1(buf,count,datatype,dest,tag,comm,request,ierror)
complex(dp) buf(:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaakm1
subroutine xyzzyaakn1(buf,count,datatype,dest,tag,comm,request,ierror)
complex(dp) buf(:,:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaakn1
subroutine xyzzyaako1(buf,count,datatype,dest,tag,comm,request,ierror)
complex(dp) buf(:,:,:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_complex)ierror=mpi_err_type
end subroutine xyzzyaako1
subroutine xyzzyaakp1(buf,count,datatype,dest,tag,comm,request,ierror)
real(dp) buf
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaakp1
subroutine xyzzyaakq1(buf,count,datatype,dest,tag,comm,request,ierror)
real(dp) buf(:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaakq1
subroutine xyzzyaakr1(buf,count,datatype,dest,tag,comm,request,ierror)
real(dp) buf(:,:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaakr1
subroutine xyzzyaaks1(buf,count,datatype,dest,tag,comm,request,ierror)
real(dp) buf(:,:,:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaaks1
subroutine xyzzyaakt1(buf,count,datatype,dest,tag,comm,request,ierror)
real(dp) buf(:,:,:,:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaakt1
subroutine xyzzyaaku1(buf,count,datatype,dest,tag,comm,request,ierror)
real(dp) buf(:,:,:,:,:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaaku1
subroutine xyzzyaakv1(buf,count,datatype,dest,tag,comm,request,ierror)
real(dp) buf(:,:,:,:,:,:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_double_precision)ierror=mpi_err_type
end subroutine xyzzyaakv1
subroutine xyzzyaakw1(buf,count,datatype,dest,tag,comm,request,ierror)
integer buf
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaakw1
subroutine xyzzyaakx1(buf,count,datatype,dest,tag,comm,request,ierror)
integer buf(:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaakx1
subroutine xyzzyaaky1(buf,count,datatype,dest,tag,comm,request,ierror)
integer buf(:,:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaaky1
subroutine xyzzyaakz1(buf,count,datatype,dest,tag,comm,request,ierror)
integer buf(:,:,:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaakz1
subroutine xyzzyaala1(buf,count,datatype,dest,tag,comm,request,ierror)
logical buf
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaala1
subroutine xyzzyaalb1(buf,count,datatype,dest,tag,comm,request,ierror)
logical buf(:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaalb1
subroutine xyzzyaalc1(buf,count,datatype,dest,tag,comm,request,ierror)
logical buf(:,:)
integer count,datatype,dest,tag,comm,request,ierror
if(datatype/=mpi_logical)ierror=mpi_err_type
end subroutine xyzzyaalc1
subroutine xyzzyaald1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,comm,ierror)
integer sendbuf,recvbuf
integer sendcount,sendtype,recvcount,recvtype,comm,ierror
recvbuf=sendbuf
if(sendtype/=mpi_integer)ierror=mpi_err_type
if(recvtype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaald1
subroutine xyzzyaale1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,comm,ierror)
integer sendbuf(:),recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,comm,ierror
recvbuf=sendbuf
if(sendtype/=mpi_integer)ierror=mpi_err_type
if(recvtype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaale1
subroutine xyzzyaalf1(sendbuf,sendcount,sendtype,recvbuf,recvcount,rec&
&vtype,comm,ierror)
integer sendbuf(:,:),recvbuf(:,:)
integer sendcount,sendtype,recvcount,recvtype,comm,ierror
recvbuf=sendbuf
if(sendtype/=mpi_integer)ierror=mpi_err_type
if(recvtype/=mpi_integer)ierror=mpi_err_type
end subroutine xyzzyaalf1
subroutine mpi_comm_split(comm1,id,node,comm,ierror)
integer comm1,id,node,comm,ierror
end subroutine mpi_comm_split
subroutine mpi_type_commit(mpi_input_type,ierror)
integer mpi_input_type,ierror
end subroutine mpi_type_commit
subroutine mpi_type_contiguous(count,oldtype,newtype,ierror)
integer count,oldtype,newtype,ierror
ierror=mpi_success
end subroutine mpi_type_contiguous
subroutine mpi_type_free(datatype,ierror)
integer datatype,ierror
ierror=mpi_success
end subroutine mpi_type_free
subroutine mpi_type_struct(count,array_of_blocklengths,array_of_displa&
&cements,array_of_types,newtype,ierror)
integer count,array_of_blocklengths(*),array_of_displacements(*),array&
&_of_types(*),newtype,ierror
end subroutine mpi_type_struct
subroutine mpi_type_size(datatype,size,ierror)
integer datatype,size,ierror
ierror=mpi_success
end subroutine mpi_type_size
subroutine mpi_type_vector(count,bl,stride,oldtype,newtype,ierror)
integer count,bl,stride,oldtype,newtype,ierror
ierror=mpi_success
end subroutine mpi_type_vector
subroutine xyzzyaalg1(sendbuf,sendcounts,sdispls,sendtype,recvbuf,recv&
&counts,rdispls,recvtype,comm,ierror)
integer sendbuf,recvbuf
integer sendcounts(:),sdispls(:),sendtype,recvcounts(:),rdispls(:),rec&
&vtype,comm,ierror
recvbuf=sendbuf
end subroutine xyzzyaalg1
subroutine xyzzyaalh1(sendbuf,sendcounts,sdispls,sendtype,recvbuf,recv&
&counts,rdispls,recvtype,comm,ierror)
real(dp) sendbuf,recvbuf
integer sendcounts(:),sdispls(:),sendtype,recvcounts(:),rdispls(:),rec&
&vtype,comm,ierror
recvbuf=sendbuf
end subroutine xyzzyaalh1
subroutine xyzzyaali1(sendbuf,sendcounts,sdispls,sendtype,recvbuf,recv&
&counts,rdispls,recvtype,comm,ierror)
integer sendbuf(:,:),recvbuf(:,:)
integer sendcounts(*),sdispls(*),sendtype,recvcounts(*),rdispls(*),rec&
&vtype,comm,ierror
recvbuf=sendbuf
end subroutine xyzzyaali1
subroutine xyzzyaalj1(sendbuf,sendcounts,sdispls,sendtype,recvbuf,recv&
&counts,rdispls,recvtype,comm,ierror)
real(dp) sendbuf(:,:),recvbuf(:,:)
integer sendcounts(*),sdispls(*),sendtype,recvcounts(*),rdispls(*),rec&
&vtype,comm,ierror
recvbuf=sendbuf
end subroutine xyzzyaalj1
subroutine xyzzyaalk1(location,address,ierror)
integer location,address,ierror
end subroutine xyzzyaalk1
subroutine xyzzyaall1(location,address,ierror)
logical location
integer address,ierror
end subroutine xyzzyaall1
subroutine xyzzyaalm1(sendbuf,sendcount,sendtype,recvbuf,recvcounts,re&
&cvtype,comm,ierror)
integer sendbuf,recvbuf(:)
integer sendcount,sendtype,recvcounts,recvtype,comm,ierror
recvbuf=sendbuf
end subroutine xyzzyaalm1
subroutine mpi_wait(request,ignore,ierror)
integer request,ignore,ierror
ierror=mpi_success
end subroutine mpi_wait
subroutine mpi_waitall(n,request,ignore,ierror)
integer n,request(n),ignore,ierror
ierror=mpi_success
end subroutine mpi_waitall
subroutine mpi_waitsome(n,request,outcount,which_nodes,ignore,ierror)
integer n,request(n),ignore,ierror,outcount,which_nodes(n)
outcount=1
which_nodes(1)=1
ierror=mpi_success
end subroutine mpi_waitsome
subroutine mpi_comm_group(comm_world,group_world,ierror)
integer comm_world,group_world,ierror
group_world=comm_world
ierror=mpi_success
end subroutine mpi_comm_group
subroutine mpi_group_incl(group_world,n,ranks,group_out,ierror)
integer group_world,n,ranks(n),group_out,ierror
group_out=group_world
ierror=mpi_success
end subroutine mpi_group_incl
subroutine mpi_comm_create(comm_world,n,comm,ierror)
integer comm_world,n,comm,ierror
comm=comm_world
ierror=mpi_success
end subroutine mpi_comm_create
real(dp) function mpi_wtime()
mpi_wtime=-1.d0
end function mpi_wtime
end module comms

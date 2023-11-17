module parallel
use comms
use dsp,          only: dp,i64
use format_utils, only: wout,i2s
implicit none
public
private dp,i64,wout,i2s
integer my_node,nnodes,ierror,status(mpi_status_size)
integer,parameter :: aveid=103,move_msg=107,instruct_msg=108,ncon_msg=&
&109,load_msg=112
logical am_master,am_slave,use_timer
integer,parameter :: max_no_messages=5000
integer nbreq(max_no_messages)
integer,allocatable :: rg_comm(:,:)
interface mpi_bcast_safe
module procedure mpi_bcast_safe_c1,mpi_bcast_safe_c2,mpi_bcast_safe_c3&
&,mpi_bcast_safe_c4,mpi_bcast_safe_c5,mpi_bcast_safe_c6,mpi_bcast_safe&
&_sc1,mpi_bcast_safe_sc4,mpi_bcast_safe_d1,mpi_bcast_safe_d2,mpi_bcast&
&_safe_d3,mpi_bcast_safe_d4,mpi_bcast_safe_d5,mpi_bcast_safe_i1,mpi_bc&
&ast_safe_i2,mpi_bcast_safe_i3,mpi_bcast_safe_i4,mpi_bcast_safe_l1,mpi&
&_bcast_safe_l2,mpi_bcast_safe_l3,mpi_bcast_safe_s1,mpi_bcast_safe_s2,&
&mpi_bcast_safe_r1,mpi_bcast_safe_r2,mpi_bcast_safe_r3,mpi_bcast_safe_&
&r4,mpi_bcast_safe_r5
end interface
interface qmpi_ssend
module procedure qmpi_ssend_i,qmpi_ssend_i1,qmpi_ssend_i2,qmpi_ssend_r&
&1,qmpi_ssend_d,qmpi_ssend_d1,qmpi_ssend_d2,qmpi_ssend_d3,qmpi_ssend_d&
&4,qmpi_ssend_d5,qmpi_ssend_d6,qmpi_ssend_c,qmpi_ssend_c1,qmpi_ssend_c&
&2,qmpi_ssend_c3
end interface
interface qmpi_isend
module procedure qmpi_isend_i,qmpi_isend_i1,qmpi_isend_i2,qmpi_isend_d&
&,qmpi_isend_d1,qmpi_isend_d2
end interface
interface qmpi_issend
module procedure qmpi_issend_i,qmpi_issend_i1,qmpi_issend_i2,qmpi_isse&
&nd_d,qmpi_issend_d1,qmpi_issend_d2
end interface
interface qmpi_recv
module procedure qmpi_recv_i,qmpi_recv_i1,qmpi_recv_i2,qmpi_recv_r1,qm&
&pi_recv_d,qmpi_recv_d1,qmpi_recv_d2,qmpi_recv_d3,qmpi_recv_d4,qmpi_re&
&cv_d5,qmpi_recv_d6,qmpi_recv_c,qmpi_recv_c1,qmpi_recv_c2,qmpi_recv_c3
end interface
interface qmpi_irecv
module procedure qmpi_irecv_i,qmpi_irecv_i1,qmpi_irecv_i2,qmpi_irecv_d&
&,qmpi_irecv_d1,qmpi_irecv_d2
end interface
interface qmpi_bcast
module procedure qmpi_bcast_i,qmpi_bcast_i1,qmpi_bcast_l,qmpi_bcast_d,&
&qmpi_bcast_d1
end interface
interface qmpi_gather
module procedure qmpi_gather_i,qmpi_gather_i1,qmpi_gather_d,qmpi_gathe&
&r_d1
end interface
interface qmpi_allgather
module procedure qmpi_allgather_i
end interface
interface qmpi_gather_in_place
module procedure qmpi_gather_in_place_i1,qmpi_gather_in_place_d1
end interface
interface qmpi_reduce
module procedure qmpi_reduce_i,qmpi_reduce_i1,qmpi_reduce_d,qmpi_reduc&
&e_d1,qmpi_reduce_d2
end interface
integer,private,parameter :: nsafe_bytes=33554432
private mpi_bcast_safe_c,mpi_bcast_safe_sc,mpi_bcast_safe_d,mpi_bcast_&
&safe_i,mpi_bcast_safe_l,mpi_bcast_safe_s,mpi_bcast_safe_r,mpi_bcast_s&
&afe_c1,mpi_bcast_safe_c2,mpi_bcast_safe_c3,mpi_bcast_safe_c4,mpi_bcas&
&t_safe_c5,mpi_bcast_safe_c6,mpi_bcast_safe_sc1,mpi_bcast_safe_sc4,mpi&
&_bcast_safe_d1,mpi_bcast_safe_d2,mpi_bcast_safe_d3,mpi_bcast_safe_d4,&
&mpi_bcast_safe_d5,mpi_bcast_safe_i1,mpi_bcast_safe_i2,mpi_bcast_safe_&
&i3,mpi_bcast_safe_i4,mpi_bcast_safe_l1,mpi_bcast_safe_l2,mpi_bcast_sa&
&fe_l3,mpi_bcast_safe_s1,mpi_bcast_safe_s2,mpi_bcast_safe_r1,mpi_bcast&
&_safe_r2,mpi_bcast_safe_r3,mpi_bcast_safe_r4,mpi_bcast_safe_r5
contains
subroutine init_parallel
use store, only : wout_inhibit_node
implicit none
call mpi_init(ierror)
call mpi_comm_size(mpi_comm_world,nnodes,ierror)
call mpi_comm_rank(mpi_comm_world,my_node,ierror)
if(my_node==0)then
am_master=.true.
am_slave=.false.
use_timer=.true.
else
am_master=.false.
am_slave=.true.
use_timer=.false.
endif
wout_inhibit_node=am_slave
end subroutine init_parallel
subroutine end_parallel
implicit none
call mpi_finalize(ierror)
end subroutine end_parallel
subroutine checkmpi(ie,errmesg)
implicit none
integer,intent(in) :: ie
character(*),intent(in) ::  errmesg
integer xyzzyaaaa15,xyzzyaaab15,xyzzyaaac15
character(mpi_max_error_string) string
if(ie/=mpi_success)then
call wout()
call wout('MPI returns with error code : '//trim(i2s(ie)))
call wout(' => CASINO error message is : '//trim(errmesg))
call mpi_error_string(ie,string,xyzzyaaab15,xyzzyaaaa15)
call wout(' => MPI error message is    : '//trim(string))
call wout()
call mpi_abort(mpi_comm_world,-1,xyzzyaaac15)
stop
endif
end subroutine checkmpi
subroutine qmpi_ssend_i(src,node,tag,rt,id)
implicit none
integer,intent(in) :: src
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_ssend(src,1,mpi_integer,node,tag,mpi_comm_world,ierror)
call checkmpi(ierror,'s-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_ssend_i
subroutine qmpi_ssend_i1(src,node,tag,rt,id)
implicit none
integer,intent(in) :: src(:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_ssend(src,size(src,1),mpi_integer,node,tag,mpi_comm_world,ier&
&ror)
call checkmpi(ierror,'s-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_ssend_i1
subroutine qmpi_ssend_i2(src,node,tag,rt,id)
implicit none
integer,intent(in) :: src(:,:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_ssend(src,size(src,1)*size(src,2),mpi_integer,node,tag,mpi_co&
&mm_world,ierror)
call checkmpi(ierror,'s-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_ssend_i2
subroutine qmpi_isend_i(src,node,tag,request,rt,id)
implicit none
integer,intent(in) :: src
integer,intent(in) :: node,tag
integer,intent(out) :: request
character(*),intent(in) :: rt,id
call mpi_isend(src,1,mpi_integer,node,tag,mpi_comm_world,request,ierro&
&r)
call checkmpi(ierror,'i-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_isend_i
subroutine qmpi_isend_i1(src,node,tag,request,rt,id)
implicit none
integer,intent(in) :: src(:)
integer,intent(in) :: node,tag
integer,intent(out) :: request
character(*),intent(in) :: rt,id
call mpi_isend(src,size(src,1),mpi_integer,node,tag,mpi_comm_world,req&
&uest,ierror)
call checkmpi(ierror,'i-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_isend_i1
subroutine qmpi_isend_i2(src,node,tag,request,rt,id)
implicit none
integer,intent(in) :: src(:,:)
integer,intent(in) :: node,tag
integer,intent(out) :: request
character(*),intent(in) :: rt,id
call mpi_isend(src,size(src,1)*size(src,2),mpi_integer,node,tag,mpi_co&
&mm_world,request,ierror)
call checkmpi(ierror,'i-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_isend_i2
subroutine qmpi_issend_i(src,node,tag,request,rt,id)
implicit none
integer,intent(in) :: src
integer,intent(in) :: node,tag
integer,intent(out) :: request
character(*),intent(in) :: rt,id
call mpi_issend(src,1,mpi_integer,node,tag,mpi_comm_world,request,ierr&
&or)
call checkmpi(ierror,'is-sending '//id//' ('//trim(i2s(tag))//') in '/&
&/rt)
end subroutine qmpi_issend_i
subroutine qmpi_issend_i1(src,node,tag,request,rt,id)
implicit none
integer,intent(in) :: src(:)
integer,intent(in) :: node,tag
integer,intent(out) :: request
character(*),intent(in) :: rt,id
call mpi_issend(src,size(src,1),mpi_integer,node,tag,mpi_comm_world,re&
&quest,ierror)
call checkmpi(ierror,'is-sending '//id//' ('//trim(i2s(tag))//') in '/&
&/rt)
end subroutine qmpi_issend_i1
subroutine qmpi_issend_i2(src,node,tag,request,rt,id)
implicit none
integer,intent(in) :: src(:,:)
integer,intent(in) :: node,tag
integer,intent(out) :: request
character(*),intent(in) :: rt,id
call mpi_issend(src,size(src,1)*size(src,2),mpi_integer,node,tag,mpi_c&
&omm_world,request,ierror)
call checkmpi(ierror,'is-sending '//id//' ('//trim(i2s(tag))//') in '/&
&/rt)
end subroutine qmpi_issend_i2
subroutine qmpi_ssend_r1(src,node,tag,rt,id)
implicit none
real,intent(in) :: src(:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_ssend(src,size(src,1),mpi_real,node,tag,mpi_comm_world,ierror&
&)
call checkmpi(ierror,'s-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_ssend_r1
subroutine qmpi_ssend_d(src,node,tag,rt,id)
implicit none
real(dp),intent(in) :: src
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_ssend(src,1,mpi_double_precision,node,tag,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'s-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_ssend_d
subroutine qmpi_ssend_d1(src,node,tag,rt,id)
implicit none
real(dp),intent(in) :: src(:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_ssend(src,size(src,1),mpi_double_precision,node,tag,mpi_comm_&
&world,ierror)
call checkmpi(ierror,'s-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_ssend_d1
subroutine qmpi_ssend_d2(src,node,tag,rt,id)
implicit none
real(dp),intent(in) :: src(:,:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_ssend(src,size(src,1)*size(src,2),mpi_double_precision,node,t&
&ag,mpi_comm_world,ierror)
call checkmpi(ierror,'s-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_ssend_d2
subroutine qmpi_ssend_d3(src,node,tag,rt,id)
implicit none
real(dp),intent(in) :: src(:,:,:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_ssend(src,size(src,1)*size(src,2)*size(src,3),mpi_double_prec&
&ision,node,tag,mpi_comm_world,ierror)
call checkmpi(ierror,'s-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_ssend_d3
subroutine qmpi_ssend_d4(src,node,tag,rt,id)
implicit none
real(dp),intent(in) :: src(:,:,:,:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_ssend(src,size(src,1)*size(src,2)*size(src,3)*size(src,4),mpi&
&_double_precision,node,tag,mpi_comm_world,ierror)
call checkmpi(ierror,'s-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_ssend_d4
subroutine qmpi_ssend_d5(src,node,tag,rt,id)
implicit none
real(dp),intent(in) :: src(:,:,:,:,:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_ssend(src,size(src,1)*size(src,2)*size(src,3)*size(src,4)*siz&
&e(src,5),mpi_double_precision,node,tag,mpi_comm_world,ierror)
call checkmpi(ierror,'s-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_ssend_d5
subroutine qmpi_ssend_d6(src,node,tag,rt,id)
implicit none
real(dp),intent(in) :: src(:,:,:,:,:,:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_ssend(src,size(src,1)*size(src,2)*size(src,3)*size(src,4)*siz&
&e(src,5)*size(src,6),mpi_double_precision,node,tag,mpi_comm_world,ier&
&ror)
call checkmpi(ierror,'s-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_ssend_d6
subroutine qmpi_isend_d(src,node,tag,request,rt,id)
implicit none
real(dp),intent(in) :: src
integer,intent(in) :: node,tag
integer,intent(out) :: request
character(*),intent(in) :: rt,id
call mpi_isend(src,1,mpi_double_precision,node,tag,mpi_comm_world,requ&
&est,ierror)
call checkmpi(ierror,'i-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_isend_d
subroutine qmpi_isend_d1(src,node,tag,request,rt,id)
implicit none
real(dp),intent(in) :: src(:)
integer,intent(in) :: node,tag
integer,intent(out) :: request
character(*),intent(in) :: rt,id
call mpi_isend(src,size(src,1),mpi_double_precision,node,tag,mpi_comm_&
&world,request,ierror)
call checkmpi(ierror,'i-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_isend_d1
subroutine qmpi_isend_d2(src,node,tag,request,rt,id)
implicit none
real(dp),intent(in) :: src(:,:)
integer,intent(in) :: node,tag
integer,intent(out) :: request
character(*),intent(in) :: rt,id
call mpi_isend(src,size(src,1)*size(src,2),mpi_double_precision,node,t&
&ag,mpi_comm_world,request,ierror)
call checkmpi(ierror,'i-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_isend_d2
subroutine qmpi_issend_d(src,node,tag,request,rt,id)
implicit none
real(dp),intent(in) :: src
integer,intent(in) :: node,tag
integer,intent(out) :: request
character(*),intent(in) :: rt,id
call mpi_issend(src,1,mpi_double_precision,node,tag,mpi_comm_world,req&
&uest,ierror)
call checkmpi(ierror,'i-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_issend_d
subroutine qmpi_issend_d1(src,node,tag,request,rt,id)
implicit none
real(dp),intent(in) :: src(:)
integer,intent(in) :: node,tag
integer,intent(out) :: request
character(*),intent(in) :: rt,id
call mpi_issend(src,size(src,1),mpi_double_precision,node,tag,mpi_comm&
&_world,request,ierror)
call checkmpi(ierror,'i-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_issend_d1
subroutine qmpi_issend_d2(src,node,tag,request,rt,id)
implicit none
real(dp),intent(in) :: src(:,:)
integer,intent(in) :: node,tag
integer,intent(out) :: request
character(*),intent(in) :: rt,id
call mpi_issend(src,size(src,1)*size(src,2),mpi_double_precision,node,&
&tag,mpi_comm_world,request,ierror)
call checkmpi(ierror,'i-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_issend_d2
subroutine qmpi_ssend_c(src,node,tag,rt,id)
implicit none
complex(dp),intent(in) :: src
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_ssend(src,1,mpi_double_complex,node,tag,mpi_comm_world,ierror&
&)
call checkmpi(ierror,'s-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_ssend_c
subroutine qmpi_ssend_c1(src,node,tag,rt,id)
implicit none
complex(dp),intent(in) :: src(:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_ssend(src,size(src,1),mpi_double_complex,node,tag,mpi_comm_wo&
&rld,ierror)
call checkmpi(ierror,'s-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_ssend_c1
subroutine qmpi_ssend_c2(src,node,tag,rt,id)
implicit none
complex(dp),intent(in) :: src(:,:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_ssend(src,size(src,1)*size(src,2),mpi_double_complex,node,tag&
&,mpi_comm_world,ierror)
call checkmpi(ierror,'s-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_ssend_c2
subroutine qmpi_ssend_c3(src,node,tag,rt,id)
implicit none
complex(dp),intent(in) :: src(:,:,:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_ssend(src,size(src,1)*size(src,2)*size(src,3),mpi_double_comp&
&lex,node,tag,mpi_comm_world,ierror)
call checkmpi(ierror,'s-sending '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_ssend_c3
subroutine qmpi_recv_i(dst,node,tag,rt,id)
implicit none
integer,intent(out) :: dst
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_recv(dst,1,mpi_integer,node,tag,mpi_comm_world,status,ierror)
call checkmpi(ierror,'receiving '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_recv_i
subroutine qmpi_recv_i1(dst,node,tag,rt,id)
implicit none
integer,intent(out) :: dst(:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_recv(dst,size(dst,1),mpi_integer,node,tag,mpi_comm_world,stat&
&us,ierror)
call checkmpi(ierror,'receiving '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_recv_i1
subroutine qmpi_recv_i2(dst,node,tag,rt,id)
implicit none
integer,intent(out) :: dst(:,:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_recv(dst,size(dst,1)*size(dst,2),mpi_integer,node,tag,mpi_com&
&m_world,status,ierror)
call checkmpi(ierror,'receiving '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_recv_i2
subroutine qmpi_recv_r1(dst,node,tag,rt,id)
implicit none
real,intent(out) :: dst(:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_recv(dst,size(dst,1),mpi_real,node,tag,mpi_comm_world,status,&
&ierror)
call checkmpi(ierror,'receiving '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_recv_r1
subroutine qmpi_recv_d(dst,node,tag,rt,id)
implicit none
real(dp),intent(out) :: dst
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_recv(dst,1,mpi_double_precision,node,tag,mpi_comm_world,statu&
&s,ierror)
call checkmpi(ierror,'receiving '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_recv_d
subroutine qmpi_recv_d1(dst,node,tag,rt,id)
implicit none
real(dp),intent(out) :: dst(:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_recv(dst,size(dst,1),mpi_double_precision,node,tag,mpi_comm_w&
&orld,status,ierror)
call checkmpi(ierror,'receiving '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_recv_d1
subroutine qmpi_recv_d2(dst,node,tag,rt,id)
implicit none
real(dp),intent(out) :: dst(:,:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_recv(dst,size(dst,1)*size(dst,2),mpi_double_precision,node,ta&
&g,mpi_comm_world,status,ierror)
call checkmpi(ierror,'receiving '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_recv_d2
subroutine qmpi_recv_d3(dst,node,tag,rt,id)
implicit none
real(dp),intent(out) :: dst(:,:,:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_recv(dst,size(dst,1)*size(dst,2)*size(dst,3),mpi_double_preci&
&sion,node,tag,mpi_comm_world,status,ierror)
call checkmpi(ierror,'receiving '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_recv_d3
subroutine qmpi_recv_d4(dst,node,tag,rt,id)
implicit none
real(dp),intent(out) :: dst(:,:,:,:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_recv(dst,size(dst,1)*size(dst,2)*size(dst,3)*size(dst,4),mpi_&
&double_precision,node,tag,mpi_comm_world,status,ierror)
call checkmpi(ierror,'receiving '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_recv_d4
subroutine qmpi_recv_d5(dst,node,tag,rt,id)
implicit none
real(dp),intent(out) :: dst(:,:,:,:,:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_recv(dst,size(dst,1)*size(dst,2)*size(dst,3)*size(dst,4)*size&
&(dst,5),mpi_double_precision,node,tag,mpi_comm_world,status,ierror)
call checkmpi(ierror,'receiving '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_recv_d5
subroutine qmpi_recv_d6(dst,node,tag,rt,id)
implicit none
real(dp),intent(out) :: dst(:,:,:,:,:,:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_recv(dst,size(dst,1)*size(dst,2)*size(dst,3)*size(dst,4)*size&
&(dst,5)*size(dst,6),mpi_double_precision,node,tag,mpi_comm_world,stat&
&us,ierror)
call checkmpi(ierror,'receiving '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_recv_d6
subroutine qmpi_recv_c(dst,node,tag,rt,id)
implicit none
complex(dp),intent(out) :: dst
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_recv(dst,1,mpi_double_complex,node,tag,mpi_comm_world,status,&
&ierror)
call checkmpi(ierror,'receiving '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_recv_c
subroutine qmpi_recv_c1(dst,node,tag,rt,id)
implicit none
complex(dp),intent(out) :: dst(:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_recv(dst,size(dst,1),mpi_double_complex,node,tag,mpi_comm_wor&
&ld,status,ierror)
call checkmpi(ierror,'receiving '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_recv_c1
subroutine qmpi_recv_c2(dst,node,tag,rt,id)
implicit none
complex(dp),intent(out) :: dst(:,:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_recv(dst,size(dst,1)*size(dst,2),mpi_double_complex,node,tag,&
&mpi_comm_world,status,ierror)
call checkmpi(ierror,'receiving '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_recv_c2
subroutine qmpi_recv_c3(dst,node,tag,rt,id)
implicit none
complex(dp),intent(out) :: dst(:,:,:)
integer,intent(in) :: node,tag
character(*),intent(in) :: rt,id
call mpi_recv(dst,size(dst,1)*size(dst,2)*size(dst,3),mpi_double_compl&
&ex,node,tag,mpi_comm_world,status,ierror)
call checkmpi(ierror,'receiving '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_recv_c3
subroutine qmpi_irecv_i(dst,node,tag,request,rt,id)
implicit none
integer,intent(out) :: dst
integer,intent(in) :: node,tag
integer,intent(out) :: request
character(*),intent(in) :: rt,id
call mpi_irecv(dst,1,mpi_integer,node,tag,mpi_comm_world,request,ierro&
&r)
call checkmpi(ierror,'receiving '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_irecv_i
subroutine qmpi_irecv_i1(dst,node,tag,request,rt,id)
implicit none
integer,intent(out) :: dst(:)
integer,intent(in) :: node,tag
integer,intent(out) :: request
character(*),intent(in) :: rt,id
call mpi_irecv(dst,size(dst,1),mpi_integer,node,tag,mpi_comm_world,req&
&uest,ierror)
call checkmpi(ierror,'receiving '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_irecv_i1
subroutine qmpi_irecv_i2(dst,node,tag,request,rt,id)
implicit none
integer,intent(out) :: dst(:,:)
integer,intent(in) :: node,tag
integer,intent(out) :: request
character(*),intent(in) :: rt,id
call mpi_irecv(dst,size(dst,1)*size(dst,2),mpi_integer,node,tag,mpi_co&
&mm_world,request,ierror)
call checkmpi(ierror,'receiving '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_irecv_i2
subroutine qmpi_irecv_d(dst,node,tag,request,rt,id)
implicit none
real(dp),intent(out) :: dst
integer,intent(in) :: node,tag
integer,intent(out) :: request
character(*),intent(in) :: rt,id
call mpi_irecv(dst,1,mpi_double_precision,node,tag,mpi_comm_world,requ&
&est,ierror)
call checkmpi(ierror,'receiving '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_irecv_d
subroutine qmpi_irecv_d1(dst,node,tag,request,rt,id)
implicit none
real(dp),intent(out) :: dst(:)
integer,intent(in) :: node,tag
integer,intent(out) :: request
character(*),intent(in) :: rt,id
call mpi_irecv(dst,size(dst,1),mpi_double_precision,node,tag,mpi_comm_&
&world,request,ierror)
call checkmpi(ierror,'receiving '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_irecv_d1
subroutine qmpi_irecv_d2(dst,node,tag,request,rt,id)
implicit none
real(dp),intent(out) :: dst(:,:)
integer,intent(in) :: node,tag
integer,intent(out) :: request
character(*),intent(in) :: rt,id
call mpi_irecv(dst,size(dst,1)*size(dst,2),mpi_double_precision,node,t&
&ag,mpi_comm_world,request,ierror)
call checkmpi(ierror,'receiving '//id//' ('//trim(i2s(tag))//') in '//&
&rt)
end subroutine qmpi_irecv_d2
subroutine qmpi_bcast_i(val,rt,id)
implicit none
integer,intent(inout) :: val
character(*),intent(in) :: rt,id
call mpi_bcast(val,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting '//id//' in '//rt)
end subroutine qmpi_bcast_i
subroutine qmpi_bcast_i1(val,rt,id)
implicit none
integer,intent(inout) :: val(:)
character(*),intent(in) :: rt,id
call mpi_bcast(val,size(val,1),mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting '//id//' in '//rt)
end subroutine qmpi_bcast_i1
subroutine qmpi_bcast_l(val,rt,id)
implicit none
logical,intent(inout) :: val
character(*),intent(in) :: rt,id
call mpi_bcast(val,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting '//id//' in '//rt)
end subroutine qmpi_bcast_l
subroutine qmpi_bcast_d(val,rt,id)
implicit none
real(dp),intent(inout) :: val
character(*),intent(in) :: rt,id
call mpi_bcast(val,1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting '//id//' in '//rt)
end subroutine qmpi_bcast_d
subroutine qmpi_bcast_d1(val,rt,id)
implicit none
real(dp),intent(inout) :: val(:)
character(*),intent(in) :: rt,id
call mpi_bcast(val,size(val,1),mpi_double_precision,0,mpi_comm_world,i&
&error)
call checkmpi(ierror,'Broadcasting '//id//' in '//rt)
end subroutine qmpi_bcast_d1
subroutine qmpi_gather_i(src,dst,rt,id)
implicit none
integer,intent(in) :: src
integer,intent(out) :: dst(nnodes)
character(*),intent(in) :: rt,id
call mpi_gather(src,1,mpi_integer,dst,1,mpi_integer,0,mpi_comm_world,i&
&error)
call checkmpi(ierror,'Gathering '//id//' in '//rt)
end subroutine qmpi_gather_i
subroutine qmpi_gather_i1(src,dst,rt,id)
implicit none
integer,intent(in) :: src(:)
integer,intent(out) :: dst(size(src,1),nnodes)
character(*),intent(in) :: rt,id
call mpi_gather(src,size(src,1),mpi_integer,dst,size(src,1),mpi_intege&
&r,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Gathering '//id//' in '//rt)
end subroutine qmpi_gather_i1
subroutine qmpi_gather_d(src,dst,rt,id)
implicit none
real(dp),intent(in) :: src
real(dp),intent(out) :: dst(nnodes)
character(*),intent(in) :: rt,id
call mpi_gather(src,1,mpi_integer,dst,1,mpi_integer,0,mpi_comm_world,i&
&error)
call checkmpi(ierror,'Gathering '//id//' in '//rt)
end subroutine qmpi_gather_d
subroutine qmpi_gather_d1(src,dst,rt,id)
implicit none
real(dp),intent(in) :: src(:)
real(dp),intent(out) :: dst(size(src,1),nnodes)
character(*),intent(in) :: rt,id
call mpi_gather(src,size(src,1),mpi_integer,dst,size(src,1),mpi_intege&
&r,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Gathering '//id//' in '//rt)
end subroutine qmpi_gather_d1
subroutine qmpi_allgather_i(src,dst,comm,rt,id)
implicit none
integer,intent(in) :: src,comm
integer,intent(out) :: dst(:)
character(*),intent(in) :: rt,id
call mpi_allgather(src,1,mpi_integer,dst,1,mpi_integer,comm,ierror)
call checkmpi(ierror,'Gathering '//id//' in '//rt)
end subroutine qmpi_allgather_i
subroutine qmpi_gather_in_place_i1(buf,rt,id)
implicit none
integer,intent(inout) :: buf(:)
character(*),intent(in) :: rt,id
integer xyzzyaaaa74(1)
if(am_master)then
call mpi_gather_in_place(1,mpi_integer,buf,1,mpi_integer,0,mpi_comm_wo&
&rld,ierror)
else
call mpi_gather(buf,1,mpi_integer,xyzzyaaaa74,1,mpi_integer,0,mpi_comm&
&_world,ierror)
endif
call checkmpi(ierror,'In-place gathering '//id//' in '//rt)
end subroutine qmpi_gather_in_place_i1
subroutine qmpi_gather_in_place_d1(buf,rt,id)
implicit none
real(dp),intent(inout) :: buf(:)
character(*),intent(in) :: rt,id
real(dp) xyzzyaaaa75(1)
if(am_master)then
call mpi_gather_in_place(1,mpi_double_precision,buf,1,mpi_double_preci&
&sion,0,mpi_comm_world,ierror)
else
call mpi_gather(buf,1,mpi_double_precision,xyzzyaaaa75,1,mpi_double_pr&
&ecision,0,mpi_comm_world,ierror)
endif
call checkmpi(ierror,'In-place gathering '//id//' in '//rt)
end subroutine qmpi_gather_in_place_d1
subroutine qmpi_reduce_i(src,dst,typ,rt,id)
implicit none
integer,intent(in) :: src
integer,intent(out) :: dst
integer,intent(in) :: typ
character(*),intent(in) :: rt,id
call mpi_reduce(src,dst,1,mpi_integer,typ,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Reducing '//id//' in '//rt)
end subroutine qmpi_reduce_i
subroutine qmpi_reduce_i1(src,dst,typ,rt,id)
implicit none
integer,intent(in) :: src(:)
integer,intent(out) :: dst(size(src,1))
integer,intent(in) :: typ
character(*),intent(in) :: rt,id
call mpi_reduce(src,dst,size(src,1),mpi_integer,typ,0,mpi_comm_world,i&
&error)
call checkmpi(ierror,'Reducing '//id//' in '//rt)
end subroutine qmpi_reduce_i1
subroutine qmpi_reduce_d(src,dst,typ,rt,id)
implicit none
real(dp),intent(in) :: src
real(dp),intent(out) :: dst
integer,intent(in) :: typ
character(*),intent(in) :: rt,id
call mpi_reduce(src,dst,1,mpi_double_precision,typ,0,mpi_comm_world,ie&
&rror)
call checkmpi(ierror,'Reducing '//id//' in '//rt)
end subroutine qmpi_reduce_d
subroutine qmpi_reduce_d1(src,dst,typ,rt,id)
implicit none
real(dp),intent(in) :: src(:)
real(dp),intent(out) :: dst(size(src,1))
integer,intent(in) :: typ
character(*),intent(in) :: rt,id
call mpi_reduce(src,dst,size(src,1),mpi_double_precision,typ,0,mpi_com&
&m_world,ierror)
call checkmpi(ierror,'Reducing '//id//' in '//rt)
end subroutine qmpi_reduce_d1
subroutine qmpi_reduce_d2(src,dst,typ,rt,id)
implicit none
real(dp),intent(in) :: src(:,:)
real(dp),intent(out) :: dst(size(src,1),size(src,2))
integer,intent(in) :: typ
character(*),intent(in) :: rt,id
call mpi_reduce(src,dst,size(src,1)*size(src,2),mpi_double_precision,t&
&yp,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Reducing '//id//' in '//rt)
end subroutine qmpi_reduce_d2
subroutine qmc_barrier
implicit none
integer xyzzyaaaa81
logical xyzzyaaab81
character(256) tmpr
if(nnodes<=128)then
call mpi_iprobe(mpi_any_source,mpi_any_tag,mpi_comm_world,xyzzyaaab81,&
&status,ierror)
call checkmpi(ierror,'QMC barrier IPROBE')
if(xyzzyaaab81)then
call wout()
call wout("MPI queue not empty at QMC barrier.")
tmpr=''
do xyzzyaaaa81=1,mpi_status_size
tmpr=trim(tmpr)//' '//trim(i2s(status(xyzzyaaaa81)))
enddo
call wout("Status field of message"//trim(tmpr))
call wout()
call mpi_abort(mpi_comm_world,-1,ierror)
stop
endif
endif
call mpi_barrier(mpi_comm_world,ierror)
call checkmpi(ierror,'QMC barrier')
end subroutine qmc_barrier
subroutine mpi_bcast_safe_c1(buffer,count,datatype,root,comm,ierror)
complex(dp) buffer(:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_c(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_c1
subroutine mpi_bcast_safe_c2(buffer,count,datatype,root,comm,ierror)
complex(dp) buffer(:,:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_c(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_c2
subroutine mpi_bcast_safe_c3(buffer,count,datatype,root,comm,ierror)
complex(dp) buffer(:,:,:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_c(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_c3
subroutine mpi_bcast_safe_c4(buffer,count,datatype,root,comm,ierror)
complex(dp) buffer(:,:,:,:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_c(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_c4
subroutine mpi_bcast_safe_c5(buffer,count,datatype,root,comm,ierror)
complex(dp) buffer(:,:,:,:,:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_c(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_c5
subroutine mpi_bcast_safe_c6(buffer,count,datatype,root,comm,ierror)
complex(dp) buffer(:,:,:,:,:,:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_c(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_c6
subroutine mpi_bcast_safe_sc1(buffer,count,datatype,root,comm,ierror)
complex buffer(:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_sc(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_sc1
subroutine mpi_bcast_safe_sc4(buffer,count,datatype,root,comm,ierror)
complex buffer(:,:,:,:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_sc(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_sc4
subroutine mpi_bcast_safe_d1(buffer,count,datatype,root,comm,ierror)
real(dp) buffer(:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_d(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_d1
subroutine mpi_bcast_safe_d2(buffer,count,datatype,root,comm,ierror)
real(dp) buffer(:,:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_d(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_d2
subroutine mpi_bcast_safe_d3(buffer,count,datatype,root,comm,ierror)
real(dp) buffer(:,:,:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_d(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_d3
subroutine mpi_bcast_safe_d4(buffer,count,datatype,root,comm,ierror)
real(dp) buffer(:,:,:,:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_d(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_d4
subroutine mpi_bcast_safe_d5(buffer,count,datatype,root,comm,ierror)
real(dp) buffer(:,:,:,:,:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_d(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_d5
subroutine mpi_bcast_safe_r1(buffer,count,datatype,root,comm,ierror)
real buffer(:)
integer(i64) count
integer datatype,root,comm,ierror
call mpi_bcast_safe_r(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_r1
subroutine mpi_bcast_safe_r2(buffer,count,datatype,root,comm,ierror)
real buffer(:,:)
integer(i64) count
integer datatype,root,comm,ierror
call mpi_bcast_safe_r(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_r2
subroutine mpi_bcast_safe_r3(buffer,count,datatype,root,comm,ierror)
real buffer(:,:,:)
integer(i64) count
integer datatype,root,comm,ierror
call mpi_bcast_safe_r(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_r3
subroutine mpi_bcast_safe_r4(buffer,count,datatype,root,comm,ierror)
real buffer(:,:,:,:)
integer(i64) count
integer datatype,root,comm,ierror
call mpi_bcast_safe_r(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_r4
subroutine mpi_bcast_safe_r5(buffer,count,datatype,root,comm,ierror)
real buffer(:,:,:,:,:)
integer(i64) count
integer datatype,root,comm,ierror
call mpi_bcast_safe_r(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_r5
subroutine mpi_bcast_safe_i1(buffer,count,datatype,root,comm,ierror)
integer buffer(:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_i(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_i1
subroutine mpi_bcast_safe_i2(buffer,count,datatype,root,comm,ierror)
integer buffer(:,:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_i(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_i2
subroutine mpi_bcast_safe_i3(buffer,count,datatype,root,comm,ierror)
integer buffer(:,:,:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_i(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_i3
subroutine mpi_bcast_safe_i4(buffer,count,datatype,root,comm,ierror)
integer buffer(:,:,:,:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_i(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_i4
subroutine mpi_bcast_safe_l1(buffer,count,datatype,root,comm,ierror)
logical buffer(:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_l(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_l1
subroutine mpi_bcast_safe_l2(buffer,count,datatype,root,comm,ierror)
logical buffer(:,:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_l(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_l2
subroutine mpi_bcast_safe_l3(buffer,count,datatype,root,comm,ierror)
logical buffer(:,:,:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_l(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_l3
subroutine mpi_bcast_safe_s1(buffer,count,datatype,root,comm,ierror)
character buffer(:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_s(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_s1
subroutine mpi_bcast_safe_s2(buffer,count,datatype,root,comm,ierror)
character buffer(:,:)
integer count,datatype,root,comm,ierror
call mpi_bcast_safe_s(buffer,count,datatype,root,comm,ierror)
end subroutine mpi_bcast_safe_s2
subroutine mpi_bcast_safe_c(buffer,count,datatype,root,comm,ierror)
implicit none
complex(dp) buffer(*)
integer count,datatype,root,comm,ierror
integer xyzzyaaaa109,xyzzyaaab109,xyzzyaaac109
xyzzyaaaa109=0
xyzzyaaab109=nsafe_bytes/16
do while(xyzzyaaaa109<count)
xyzzyaaac109=min(count-xyzzyaaaa109,xyzzyaaab109)
call mpi_bcast(buffer(xyzzyaaaa109+1:xyzzyaaaa109+xyzzyaaac109),xyzzya&
&aac109,datatype,root,comm,ierror)
xyzzyaaaa109=xyzzyaaaa109+xyzzyaaac109
enddo
end subroutine mpi_bcast_safe_c
subroutine mpi_bcast_safe_sc(buffer,count,datatype,root,comm,ierror)
implicit none
complex buffer(*)
integer count,datatype,root,comm,ierror
integer xyzzyaaaa110,xyzzyaaab110,xyzzyaaac110
xyzzyaaaa110=0
xyzzyaaab110=nsafe_bytes/8
do while(xyzzyaaaa110<count)
xyzzyaaac110=min(count-xyzzyaaaa110,xyzzyaaab110)
call mpi_bcast(buffer(xyzzyaaaa110+1:xyzzyaaaa110+xyzzyaaac110),xyzzya&
&aac110,datatype,root,comm,ierror)
xyzzyaaaa110=xyzzyaaaa110+xyzzyaaac110
enddo
end subroutine mpi_bcast_safe_sc
subroutine mpi_bcast_safe_d(buffer,count,datatype,root,comm,ierror)
implicit none
real(dp) buffer(*)
integer count,datatype,root,comm,ierror
integer xyzzyaaaa111,xyzzyaaab111,xyzzyaaac111
xyzzyaaaa111=0
xyzzyaaab111=nsafe_bytes/8
do while(xyzzyaaaa111<count)
xyzzyaaac111=min(count-xyzzyaaaa111,xyzzyaaab111)
call mpi_bcast(buffer(xyzzyaaaa111+1:xyzzyaaaa111+xyzzyaaac111),xyzzya&
&aac111,datatype,root,comm,ierror)
xyzzyaaaa111=xyzzyaaaa111+xyzzyaaac111
enddo
end subroutine mpi_bcast_safe_d
subroutine mpi_bcast_safe_i(buffer,count,datatype,root,comm,ierror)
implicit none
integer buffer(*)
integer count,datatype,root,comm,ierror
integer xyzzyaaaa112,xyzzyaaab112,xyzzyaaac112
xyzzyaaaa112=0
xyzzyaaab112=nsafe_bytes/4
do while(xyzzyaaaa112<count)
xyzzyaaac112=min(count-xyzzyaaaa112,xyzzyaaab112)
call mpi_bcast(buffer(xyzzyaaaa112+1:xyzzyaaaa112+xyzzyaaac112),xyzzya&
&aac112,datatype,root,comm,ierror)
xyzzyaaaa112=xyzzyaaaa112+xyzzyaaac112
enddo
end subroutine mpi_bcast_safe_i
subroutine mpi_bcast_safe_l(buffer,count,datatype,root,comm,ierror)
implicit none
logical buffer(*)
integer count,datatype,root,comm,ierror
integer xyzzyaaaa113,xyzzyaaab113,xyzzyaaac113
xyzzyaaaa113=0
xyzzyaaab113=nsafe_bytes/4
do while(xyzzyaaaa113<count)
xyzzyaaac113=min(count-xyzzyaaaa113,xyzzyaaab113)
call mpi_bcast(buffer(xyzzyaaaa113+1:xyzzyaaaa113+xyzzyaaac113),xyzzya&
&aac113,datatype,root,comm,ierror)
xyzzyaaaa113=xyzzyaaaa113+xyzzyaaac113
enddo
end subroutine mpi_bcast_safe_l
subroutine mpi_bcast_safe_s(buffer,count,datatype,root,comm,ierror)
implicit none
character buffer(*)
integer count,datatype,root,comm,ierror
integer xyzzyaaaa114,xyzzyaaab114,xyzzyaaac114
xyzzyaaaa114=0
xyzzyaaab114=nsafe_bytes
do while(xyzzyaaaa114<count)
xyzzyaaac114=min(count-xyzzyaaaa114,xyzzyaaab114)
call mpi_bcast(buffer(xyzzyaaaa114+1:xyzzyaaaa114+xyzzyaaac114),xyzzya&
&aac114,datatype,root,comm,ierror)
xyzzyaaaa114=xyzzyaaaa114+xyzzyaaac114
enddo
end subroutine mpi_bcast_safe_s
subroutine mpi_bcast_safe_r(buffer,count,datatype,root,comm,ierror)
implicit none
real buffer(*)
integer datatype,root,comm,ierror,xyzzyaaaa115
integer(i64) count,xyzzyaaab115,xyzzyaaac115
xyzzyaaab115=0
xyzzyaaac115=nsafe_bytes/4
do while(xyzzyaaab115<count)
xyzzyaaaa115=int(min(count-xyzzyaaab115,xyzzyaaac115))
call mpi_bcast(buffer(xyzzyaaab115+1),xyzzyaaaa115,datatype,root,comm,&
&ierror)
xyzzyaaab115=xyzzyaaab115+xyzzyaaaa115
enddo
end subroutine mpi_bcast_safe_r
end module parallel

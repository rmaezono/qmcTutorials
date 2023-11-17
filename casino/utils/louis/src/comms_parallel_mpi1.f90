module comms
implicit none
include 'mpif.h'
interface mpi_gather_in_place
module procedure mpi_gather_in_place_i1,mpi_gather_in_place_d1
end interface mpi_gather_in_place
contains
subroutine mpi_gather_in_place_i1(sendcount,sendtype,recvbuf,recvcount&
&,recvtype,root,comm,ierror)
implicit none
integer recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
integer,allocatable :: xyzzyaaaa3(:)
integer xyzzyaaab3(1),xyzzyaaac3,xyzzyaaad3
call mpi_comm_rank(comm,xyzzyaaad3,ierror)
if(ierror/=0)return
if(xyzzyaaad3==root)then
allocate(xyzzyaaaa3(sendcount),stat=xyzzyaaac3)
if(xyzzyaaac3/=0)then
ierror=mpi_err_unknown
return
endif
xyzzyaaaa3(1:sendcount)=recvbuf(1:sendcount)
call mpi_gather(xyzzyaaaa3,sendcount,sendtype,recvbuf,recvcount,recvty&
&pe,root,comm,ierror)
deallocate(xyzzyaaaa3)
else
call mpi_gather(recvbuf,sendcount,sendtype,xyzzyaaab3(1),recvcount,rec&
&vtype,recvtype,root,comm,ierror)
endif
end subroutine mpi_gather_in_place_i1
subroutine mpi_gather_in_place_d1(sendcount,sendtype,recvbuf,recvcount&
&,recvtype,root,comm,ierror)
implicit none
double precision recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
double precision,allocatable :: xyzzyaaaa4(:)
integer xyzzyaaab4,xyzzyaaac4
double precision xyzzyaaad4(1)
call mpi_comm_rank(comm,xyzzyaaac4,ierror)
if(ierror/=0)return
if(xyzzyaaac4==root)then
allocate(xyzzyaaaa4(sendcount),stat=xyzzyaaab4)
if(xyzzyaaab4/=0)then
ierror=mpi_err_unknown
return
endif
xyzzyaaaa4(1:sendcount)=recvbuf(1:sendcount)
call mpi_gather(xyzzyaaaa4,sendcount,sendtype,recvbuf,recvcount,recvty&
&pe,root,comm,ierror)
deallocate(xyzzyaaaa4)
else
call mpi_gather(recvbuf,sendcount,sendtype,xyzzyaaad4(1),recvcount,rec&
&vtype,recvtype,root,comm,ierror)
endif
end subroutine mpi_gather_in_place_d1
end module comms

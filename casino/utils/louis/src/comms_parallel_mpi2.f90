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
call mpi_gather(mpi_in_place,sendcount,sendtype,recvbuf,recvcount,recv&
&type,root,comm,ierror)
end subroutine mpi_gather_in_place_i1
subroutine mpi_gather_in_place_d1(sendcount,sendtype,recvbuf,recvcount&
&,recvtype,root,comm,ierror)
implicit none
double precision recvbuf(:)
integer sendcount,sendtype,recvcount,recvtype,root,comm,ierror
call mpi_gather(mpi_in_place,sendcount,sendtype,recvbuf,recvcount,recv&
&type,root,comm,ierror)
end subroutine mpi_gather_in_place_d1
end module comms

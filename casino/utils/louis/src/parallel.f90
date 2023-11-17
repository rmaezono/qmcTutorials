MODULE parallel
!-----------------------------------------------------------------------------!
! Module for parallel things.                                                 !
!                                                                             !
! MDT 12/2009                                                                 !
!                                                                             !
! Changes                                                                     !
! -------                                                                     !
! None yet.                                                                   !
!-----------------------------------------------------------------------------!
 USE comms
 USE dsp
 IMPLICIT NONE
 INTEGER my_node,nnodes,ierror,status(mpi_status_size)
 INTEGER,PARAMETER :: aveid=103,move_msg=107,instruct_msg=108,ncon_msg=109,&
  &load_msg=112
 LOGICAL am_master,am_slave,use_timer

! Safe broadcast routine: it does several small broadcasts rather than one
! so big that something overflows inside the MPI library.
! NB, scalars not allowed.
 INTERFACE mpi_bcast_safe
  MODULE PROCEDURE &
   &mpi_bcast_safe_c1,mpi_bcast_safe_c2,mpi_bcast_safe_c3,&
   &mpi_bcast_safe_c4,mpi_bcast_safe_c5,mpi_bcast_safe_c6,&
   &mpi_bcast_safe_sc1,mpi_bcast_safe_sc5,&
   &mpi_bcast_safe_d1,mpi_bcast_safe_d2,mpi_bcast_safe_d3,&
   &mpi_bcast_safe_d4,mpi_bcast_safe_d5,&
   &mpi_bcast_safe_i1,mpi_bcast_safe_i2,mpi_bcast_safe_i3,&
   &mpi_bcast_safe_i4,&
   &mpi_bcast_safe_l1,mpi_bcast_safe_l2,mpi_bcast_safe_l3,&
   &mpi_bcast_safe_s1,mpi_bcast_safe_s2,mpi_bcast_safe_r1,&
   &mpi_bcast_safe_r2,mpi_bcast_safe_r3,mpi_bcast_safe_r4,mpi_bcast_safe_r5
 END INTERFACE
 INTEGER,PRIVATE,PARAMETER :: nsafe_bytes=33554432 ! 32 MiB (32 * 1024 * 1024 B)
 PRIVATE mpi_bcast_safe_c,mpi_bcast_safe_sc,mpi_bcast_safe_d,&
  &mpi_bcast_safe_i,mpi_bcast_safe_l,mpi_bcast_safe_s,mpi_bcast_safe_r,&
  &mpi_bcast_safe_c1,mpi_bcast_safe_c2,mpi_bcast_safe_c3,mpi_bcast_safe_c4,&
  &mpi_bcast_safe_c5,mpi_bcast_safe_c6,mpi_bcast_safe_sc1,mpi_bcast_safe_sc5,&
  &mpi_bcast_safe_d1,mpi_bcast_safe_d2,mpi_bcast_safe_d3,mpi_bcast_safe_d4,&
  &mpi_bcast_safe_d5,mpi_bcast_safe_i1,mpi_bcast_safe_i2,mpi_bcast_safe_i3,&
  &mpi_bcast_safe_i4,mpi_bcast_safe_l1,mpi_bcast_safe_l2,mpi_bcast_safe_l3,&
  &mpi_bcast_safe_s1,mpi_bcast_safe_s2,mpi_bcast_safe_r1,mpi_bcast_safe_r2,&
  &mpi_bcast_safe_r3,mpi_bcast_safe_r4,mpi_bcast_safe_r5


CONTAINS


 SUBROUTINE init_parallel
 IMPLICIT NONE
 call mpi_init(ierror)
 call mpi_comm_size(mpi_comm_world,nnodes,ierror)
 call mpi_comm_rank(mpi_comm_world,my_node,ierror)
 if(my_node==0)then
  am_master=.true.
  am_slave=.false.
  use_timer=.true. ! can turn this off in input
 else
  am_master=.false.
  am_slave=.true.
  use_timer=.false.
 endif
 END SUBROUTINE init_parallel


 SUBROUTINE end_parallel
 IMPLICIT NONE
 call mpi_finalize(ierror)
 END SUBROUTINE end_parallel


 SUBROUTINE checkmpi(ie,errmesg)
!-------------------------------------------------!
! Checks the value of ierror after each MPI call. !
!-------------------------------------------------!
 USE store,ONLY : o
 IMPLICIT NONE
 INTEGER,INTENT(in) :: ie
 CHARACTER(*),INTENT(in) ::  errmesg
 INTEGER ierr,resultlen,ierror
 CHARACTER(mpi_max_error_string) string
 if(ie/=mpi_success)then
  write(o,*)
  write(o,*)'MPI returns with error code : ',ie
  write(o,*)' => LOUIS error message is : ',errmesg
  call mpi_error_string(ie,string,resultlen,ierr)
  write(o,*)' => MPI error message is    : ',string
  write(o,*)
  call mpi_abort(mpi_comm_world,-1,ierror)
  call mpi_finalize(ierror)
 endif
 END SUBROUTINE checkmpi


 SUBROUTINE barrier
 IMPLICIT NONE
 call mpi_barrier(mpi_comm_world,ierror)
 END SUBROUTINE barrier


! Overloading routines for 'safe' broadcasts
 SUBROUTINE mpi_bcast_safe_c1(buffer,count,datatype,root,comm,ierror)
 COMPLEX(dp) buffer(:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_c(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_c1
 SUBROUTINE mpi_bcast_safe_c2(buffer,count,datatype,root,comm,ierror)
 COMPLEX(dp) buffer(:,:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_c(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_c2
 SUBROUTINE mpi_bcast_safe_c3(buffer,count,datatype,root,comm,ierror)
 COMPLEX(dp) buffer(:,:,:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_c(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_c3
 SUBROUTINE mpi_bcast_safe_c4(buffer,count,datatype,root,comm,ierror)
 COMPLEX(dp) buffer(:,:,:,:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_c(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_c4
 SUBROUTINE mpi_bcast_safe_c5(buffer,count,datatype,root,comm,ierror)
 COMPLEX(dp) buffer(:,:,:,:,:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_c(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_c5
 SUBROUTINE mpi_bcast_safe_c6(buffer,count,datatype,root,comm,ierror)
 COMPLEX(dp) buffer(:,:,:,:,:,:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_c(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_c6
 SUBROUTINE mpi_bcast_safe_sc1(buffer,count,datatype,root,comm,ierror)
 COMPLEX buffer(:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_sc(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_sc1
 SUBROUTINE mpi_bcast_safe_sc5(buffer,count,datatype,root,comm,ierror)
 COMPLEX buffer(:,:,:,:,:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_sc(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_sc5
 SUBROUTINE mpi_bcast_safe_d1(buffer,count,datatype,root,comm,ierror)
 REAL(dp) buffer(:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_d(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_d1
 SUBROUTINE mpi_bcast_safe_d2(buffer,count,datatype,root,comm,ierror)
 REAL(dp) buffer(:,:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_d(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_d2
 SUBROUTINE mpi_bcast_safe_d3(buffer,count,datatype,root,comm,ierror)
 REAL(dp) buffer(:,:,:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_d(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_d3
 SUBROUTINE mpi_bcast_safe_d4(buffer,count,datatype,root,comm,ierror)
 REAL(dp) buffer(:,:,:,:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_d(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_d4
 SUBROUTINE mpi_bcast_safe_d5(buffer,count,datatype,root,comm,ierror)
 REAL(dp) buffer(:,:,:,:,:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_d(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_d5
 SUBROUTINE mpi_bcast_safe_r1(buffer,count,datatype,root,comm,ierror)
 REAL buffer(:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_r(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_r1
 SUBROUTINE mpi_bcast_safe_r2(buffer,count,datatype,root,comm,ierror)
 REAL buffer(:,:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_r(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_r2
 SUBROUTINE mpi_bcast_safe_r3(buffer,count,datatype,root,comm,ierror)
 REAL buffer(:,:,:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_r(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_r3
 SUBROUTINE mpi_bcast_safe_r4(buffer,count,datatype,root,comm,ierror)
 REAL buffer(:,:,:,:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_r(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_r4
 SUBROUTINE mpi_bcast_safe_r5(buffer,count,datatype,root,comm,ierror)
 REAL buffer(:,:,:,:,:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_r(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_r5
 SUBROUTINE mpi_bcast_safe_i1(buffer,count,datatype,root,comm,ierror)
 INTEGER buffer(:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_i(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_i1
 SUBROUTINE mpi_bcast_safe_i2(buffer,count,datatype,root,comm,ierror)
 INTEGER buffer(:,:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_i(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_i2
 SUBROUTINE mpi_bcast_safe_i3(buffer,count,datatype,root,comm,ierror)
 INTEGER buffer(:,:,:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_i(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_i3
 SUBROUTINE mpi_bcast_safe_i4(buffer,count,datatype,root,comm,ierror)
 INTEGER buffer(:,:,:,:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_i(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_i4
 SUBROUTINE mpi_bcast_safe_l1(buffer,count,datatype,root,comm,ierror)
 LOGICAL buffer(:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_l(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_l1
 SUBROUTINE mpi_bcast_safe_l2(buffer,count,datatype,root,comm,ierror)
 LOGICAL buffer(:,:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_l(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_l2
 SUBROUTINE mpi_bcast_safe_l3(buffer,count,datatype,root,comm,ierror)
 LOGICAL buffer(:,:,:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_l(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_l3
 SUBROUTINE mpi_bcast_safe_s1(buffer,count,datatype,root,comm,ierror)
 CHARACTER buffer(:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_s(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_s1
 SUBROUTINE mpi_bcast_safe_s2(buffer,count,datatype,root,comm,ierror)
 CHARACTER buffer(:,:)
 INTEGER count,datatype,root,comm,ierror
 call mpi_bcast_safe_s(buffer,count,datatype,root,comm,ierror)
 END SUBROUTINE mpi_bcast_safe_s2


! Actual routines for 'safe' broadcasts


 SUBROUTINE mpi_bcast_safe_c(buffer,count,datatype,root,comm,ierror)
 IMPLICIT NONE
 COMPLEX(dp) buffer(*)
 INTEGER count,datatype,root,comm,ierror
 INTEGER ioffset,nsafe_elements,nsize
 ioffset=0 ; nsafe_elements=nsafe_bytes/16 ! bytes per dp complex
 do while(ioffset<count)
  nsize=min(count-ioffset,nsafe_elements)
  call mpi_bcast(buffer(ioffset+1:ioffset+nsize),nsize,datatype,root,comm,&
   &ierror)
  ioffset=ioffset+nsize
 enddo
 END SUBROUTINE mpi_bcast_safe_c

 SUBROUTINE mpi_bcast_safe_sc(buffer,count,datatype,root,comm,ierror)
 IMPLICIT NONE
 COMPLEX buffer(*)
 INTEGER count,datatype,root,comm,ierror
 INTEGER ioffset,nsafe_elements,nsize
 ioffset=0 ; nsafe_elements=nsafe_bytes/8 ! bytes per sp complex
 do while(ioffset<count)
  nsize=min(count-ioffset,nsafe_elements)
  call mpi_bcast(buffer(ioffset+1:ioffset+nsize),nsize,datatype,root,comm,&
   &ierror)
  ioffset=ioffset+nsize
 enddo
 END SUBROUTINE mpi_bcast_safe_sc

 SUBROUTINE mpi_bcast_safe_d(buffer,count,datatype,root,comm,ierror)
 IMPLICIT NONE
 REAL(dp) buffer(*)
 INTEGER count,datatype,root,comm,ierror
 INTEGER ioffset,nsafe_elements,nsize
 ioffset=0 ; nsafe_elements=nsafe_bytes/8 ! bytes per dp real
 do while(ioffset<count)
  nsize=min(count-ioffset,nsafe_elements)
  call mpi_bcast(buffer(ioffset+1:ioffset+nsize),nsize,datatype,root,comm,&
   &ierror)
  ioffset=ioffset+nsize
 enddo
 END SUBROUTINE mpi_bcast_safe_d

 SUBROUTINE mpi_bcast_safe_i(buffer,count,datatype,root,comm,ierror)
 IMPLICIT NONE
 INTEGER buffer(*)
 INTEGER count,datatype,root,comm,ierror
 INTEGER ioffset,nsafe_elements,nsize
 ioffset=0 ; nsafe_elements=nsafe_bytes/4 ! bytes per integer
 do while(ioffset<count)
  nsize=min(count-ioffset,nsafe_elements)
  call mpi_bcast(buffer(ioffset+1:ioffset+nsize),nsize,datatype,root,comm,&
   &ierror)
  ioffset=ioffset+nsize
 enddo
 END SUBROUTINE mpi_bcast_safe_i

 SUBROUTINE mpi_bcast_safe_l(buffer,count,datatype,root,comm,ierror)
 IMPLICIT NONE
 LOGICAL buffer(*)
 INTEGER count,datatype,root,comm,ierror
 INTEGER ioffset,nsafe_elements,nsize
 ioffset=0 ; nsafe_elements=nsafe_bytes/4 ! bytes per boolean
 do while(ioffset<count)
  nsize=min(count-ioffset,nsafe_elements)
  call mpi_bcast(buffer(ioffset+1:ioffset+nsize),nsize,datatype,root,comm,&
   &ierror)
  ioffset=ioffset+nsize
 enddo
 END SUBROUTINE mpi_bcast_safe_l

 SUBROUTINE mpi_bcast_safe_s(buffer,count,datatype,root,comm,ierror)
 IMPLICIT NONE
 CHARACTER buffer(*)
 INTEGER count,datatype,root,comm,ierror
 INTEGER ioffset,nsafe_elements,nsize
 ioffset=0 ; nsafe_elements=nsafe_bytes ! 1 byte per character (assume ascii)
 do while(ioffset<count)
  nsize=min(count-ioffset,nsafe_elements)
  call mpi_bcast(buffer(ioffset+1:ioffset+nsize),nsize,datatype,root,comm,&
   &ierror)
  ioffset=ioffset+nsize
 enddo
 END SUBROUTINE mpi_bcast_safe_s

 SUBROUTINE mpi_bcast_safe_r(buffer,count,datatype,root,comm,ierror)
 IMPLICIT NONE
 REAL buffer(*)
 INTEGER count,datatype,root,comm,ierror
 INTEGER ioffset,nsafe_elements,nsize
 ioffset=0 ; nsafe_elements=nsafe_bytes/4 ! bytes per sp real
 do while(ioffset<count)
  nsize=min(count-ioffset,nsafe_elements)
  call mpi_bcast(buffer(ioffset+1:ioffset+nsize),nsize,datatype,root,comm,&
   &ierror)
  ioffset=ioffset+nsize
 enddo
 END SUBROUTINE mpi_bcast_safe_r


END MODULE parallel

module openmp_base
!$ use omp_lib
!$ use parallel, only : am_master
implicit none
private
public openmp_setup,use_openmp
!$ public get_omp_index_range
!$ public nthreads,my_thread
!$ public ne_min_ee,nitot_min,newald_min
!$ integer,parameter:: ne_min_ee=100
!$ integer,parameter:: nitot_min=200
!$ integer,parameter:: newald_min=100
!$ integer nthreads,my_thread
!$omp threadprivate(my_thread)
logical :: use_openmp=.false.
contains
subroutine openmp_setup()
implicit none
!$ use_openmp=.true.
!$ call omp_set_dynamic(.false.)
!$omp parallel default(none) shared(am_master,nthreads)
!$ my_thread=omp_get_thread_num()
!$ nthreads=omp_get_num_threads()
!$omp end parallel
end subroutine openmp_setup
!$ subroutine get_omp_index_range(istart,iend)
!$ implicit none
!$ integer,intent(inout) :: istart,iend
!$ integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3,xyzzyaaad3
!$ xyzzyaaaa3=iend-istart+1
!$ xyzzyaaac3=xyzzyaaaa3/nthreads
!$ xyzzyaaad3=mod(xyzzyaaaa3,nthreads)
!$ istart=istart+xyzzyaaac3*my_thread
!$ iend=istart+xyzzyaaac3-1
!$ if(xyzzyaaad3/=0)then
!$  xyzzyaaab3=nthreads-xyzzyaaad3
!$  if(my_thread>=xyzzyaaab3)then
!$   istart=istart+my_thread-xyzzyaaab3
!$   iend=iend+my_thread-xyzzyaaab3+1
!$  endif
!$ endif
!$ end subroutine get_omp_index_range
end module openmp_base

program main
use openmp_base
use parallel
use format_utils,only : wout,i2s,global_time_heading,print_centred_lin&
&e,write_list_int
use run_control, only : timer_start,timer_end
use store,       only : o,output_file
use shalloc,     only : init_shm,need_shm,smp_nodes,nnpsmp,shalloc_cle&
&an
implicit none
interface
subroutine monte_carlo
end subroutine monte_carlo
end interface
integer xyzzyaaaa1
character(256) version
include 'VERSION'
call init_parallel
call timer_start
call init_shm(nnodes,my_node)
call openmp_setup
o=7
if(am_master)then
output_file='out'
else
output_file='.out_node'//trim(i2s(my_node))
endif
if(am_master)then
call wout(repeat('-',78))
call wout()
call print_centred_line(' #####                                     ',&
&80)
call print_centred_line('##   ##    ##     ####   ##  #   ##   #### ',&
&80)
call print_centred_line('##        ####   ##      ##  ##  ##  ##  ##',&
&80)
call print_centred_line('##       ##  ##   ####   ##  ### ##  ##  ##',&
&80)
call print_centred_line('##       ######      ##  ##  ## ###  ##  ##',&
&80)
call print_centred_line('##   ##  ##  ##  ##  ##  ##  ##  ##  ##  ##',&
&80)
call print_centred_line(' #####   ##  ##   ####   ##  ##   #   #### ',&
&80)
call wout()
call print_centred_line('The Cambridge Quantum Monte Carlo Code',80)
call print_centred_line('CASINO '//trim(version),80)
call wout()
call print_centred_line('Main Authors : R.J.Needs, M.D.Towler, N.D.Dru&
&mmond and P.Lopez Rios',80)
call wout()
call print_centred_line('Theory of Condensed Matter Group, Cavendish L&
&aboratory,',80)
call print_centred_line('University of Cambridge, Cambridge, CB3 OHE, &
&UK.',80)
call wout()
call print_centred_line('CASINO web page: vallico.net/casinoqmc',80)
call wout()
call print_centred_line('Current contact: mdt26 at cam.ac.uk',80)
call wout()
call wout(repeat('-',78))
call global_time_heading(.true.)
call wout()
if(nnodes==1)then
call wout('Sequential run: not using MPI.')
else
call wout('Running in parallel using '//trim(i2s(nnodes))//' MPI proce&
&sses.')
endif
!$ if(nthreads==1)then
!$  call wout('Running a single thread per process.')
!$ else
!$  call wout('Running '//trim(i2s(nthreads))//' threads per process.'&
!$ &)
!$ endif
call wout()
if(need_shm)then
call wout('Using shared memory within nodes. Ranks of processes on nod&
&e 1:')
xyzzyaaaa1=floor(log10(dble(max(maxval(smp_nodes(1:nnpsmp)),1))))+1
call write_list_int(nnpsmp,smp_nodes(1:nnpsmp),78/(xyzzyaaaa1+1),xyzzy&
&aaaa1,1)
call wout()
deallocate(smp_nodes)
endif
endif
call monte_carlo
if(am_master)then
call timer_end
call global_time_heading(.false.)
endif
call shalloc_clean
call end_parallel
end program main

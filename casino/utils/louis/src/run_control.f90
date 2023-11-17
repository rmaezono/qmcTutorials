MODULE run_control
!---------------------------------------------------------------------!
! Flow control module, comprising:                                    !
! - Error-handling routines:                                          !
!   errstop,errstop2,errstop_master,errstop_quiet,errwarn,all_stop    !
! - Timing routines:                                                  !
!   timer,timer_start,timer_end,tcputime,walltime,exceeds_time_limit  !
! - Other routines:                                                   !
!   reallocate                                                        !
!                                                                     !
! The timing routines store the routine-calling structure using a     !
! linked-tree hierarchical structure. Items in the tree can be added, !
! closed and re-opened by calling TIMER below.                        !
!---------------------------------------------------------------------!
 USE dsp
 USE format_utils, ONLY : i2s,r2ss,wordwrap,time2human
 USE parallel,     ONLY : am_master
 USE store,        ONLY : o,errstop_skip,max_cpu_time,max_real_time,ialloc
 IMPLICIT NONE
 INTERFACE reallocate
  MODULE PROCEDURE reallocate_rv,reallocate_rm,reallocate_iv,reallocate_im,&
   &reallocate_hv
 END INTERFACE
 PRIVATE
! Public routines:
 PUBLIC timer,timer_start,timer_end,tcputime,walltime,exceeds_time_limit,&
  &errstop,errstop_master,errstop2,errstop_quiet,errwarn,all_stop,&
  &loop_time_estimate,reallocate

 REAL(sp),PARAMETER :: inv_wall_rate=0.001
! Hard limit on number of nested timers, preventing the code from eating
! RAM away forever when timers are not deactivated.
 INTEGER,PARAMETER :: hard_max_level=20
! Traceback facility for deactivated timers
 CHARACTER(30) :: tb_name(1:hard_max_level)='(unknown)'
! Flag to force showing the whole structure of 'collapsed' routines.
 LOGICAL,PARAMETER :: ignore_collapse=.false.
! Display the routine structure *as they are called* (for debugging)
 LOGICAL,PARAMETER :: debug_structure=.false.
 REAL(sp) start_time(3),start_wtime

 TYPE timed_routine
! This derived type contains the timing information of a particular
! routine in the hierarchy, as well as pointers to its parent, its
! first child and its contiguous sibling.
  REAL(sp) time(3),start(3),rtime(3)
  CHARACTER(30) name
  LOGICAL error,collapse
  INTEGER(i64) hits
  TYPE(timed_routine),POINTER :: first_sublevel,next,parent
 END TYPE timed_routine
 TYPE(timed_routine),POINTER :: level_zero,pt_current
 INTEGER :: ilevel_current=0


CONTAINS


 SUBROUTINE timer(label,activate,collapse)
!-------------------------------------------------------------!
! Keep track of the time spent by routine LABEL, starting the !
! timer if ACTIVATE or stopping it if not. The optional flag  !
! COLLAPSE will prevent routines nested inside LABEL from     !
! being shown in the final timer report.                      !
! PLR 01.2007                                                 !
!-------------------------------------------------------------!
 USE parallel, ONLY : use_timer,am_master
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: label
 LOGICAL,INTENT(in) :: activate
 LOGICAL,INTENT(in),OPTIONAL :: collapse
 LOGICAL err
 REAL(sp) current_time(3),time_diff(3)
 TYPE(timed_routine),POINTER :: pt_search

 if(debug_structure.and.am_master)then
  if(activate)write(o,'(a)')'+'//repeat(' ',ilevel_current)//trim(label)
 endif

 if(.not.use_timer)then
  if(activate)then
   ilevel_current=min(ilevel_current+1,hard_max_level)
   tb_name(ilevel_current)=trim(label)
  else
   ilevel_current=max(ilevel_current-1,0)
  endif
  return
 endif

 if(activate)then
  err=.false.
  if(ilevel_current>=hard_max_level)return
  ilevel_current=ilevel_current+1
  tb_name(ilevel_current)=trim(label)
! Determine whether this is a new routine or an already-known one
  pt_search=>pt_current%first_sublevel
  if(associated(pt_search,pt_current))then
! Current level has no children: create first child
   call create_timed_routine(pt_current,pt_search,label,.true.,err)
  else
! Search for label in tree after timed routine
   do
    if(trim(pt_search%name)==label)exit ! found
    if(associated(pt_search%next,pt_current))then
! End of children from current level: add a new child
     call create_timed_routine(pt_current,pt_search,label,.false.,err)
     exit
    endif
    pt_search=>pt_search%next ! next child
   enddo
  endif
  if(err)return ! ignore allocation errors
! Point pt_current at the new routine, and update its contents.
  pt_current=>pt_search
  pt_current%hits=pt_current%hits+1_i64
  if(present(collapse))pt_current%collapse=pt_current%collapse.or.&
   &collapse
  call cputime(current_time)
  pt_current%start(:)=current_time(:)
 else ! deactivate
  call cputime(current_time)
  if(trim(pt_current%name)==label)then
! Deactivate current timer
   time_diff(:)=current_time(:)-pt_current%start(:)
   pt_current%time(:)=pt_current%time(:)+time_diff(:)
   pt_current=>pt_current%parent
   ilevel_current=ilevel_current-1
  else
! This is an error. Try finding timer label in parent levels.
   pt_search=>pt_current
   do
    pt_search=>pt_search%parent
    if(associated(pt_search,level_zero))return ! timer not found - no action
    if(trim(pt_search%name)==label)exit ! timer found
   enddo
! Timer found, deactivate all timers in the middle and flag the error
   do
    time_diff(:)=current_time(:)-pt_current%start(:)
    pt_current%time(:)=pt_current%time(:)+time_diff(:)
    if(associated(pt_current,pt_search))exit
    pt_current%error=.true.
    pt_current=>pt_current%parent
    ilevel_current=ilevel_current-1
   enddo
  endif ! correct deactivation or not
 endif ! activate or deactivate timer

 END SUBROUTINE timer


 SUBROUTINE create_timed_routine(my_parent,new_level,label,first,err)
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: label
 LOGICAL,INTENT(in) :: first
 LOGICAL,INTENT(out) :: err
 TYPE(timed_routine),POINTER :: my_parent,new_level,pt_creation
 INTEGER ialloc
 err=.false. ; allocate(pt_creation,stat=ialloc)
 if(ialloc/=0)then
  err=.true. ; return
 endif
! Links to the new item
 if(first)then
  my_parent%first_sublevel=>pt_creation
 else
  new_level%next=>pt_creation
 endif
 new_level=>pt_creation
 new_level%name=label
 new_level%time=0.0
 new_level%rtime=0.0
 new_level%start=0.0
 new_level%hits=0
 new_level%error=.false.
 new_level%collapse=.false.
! Links from new item
 new_level%first_sublevel=>new_level
 new_level%next=>my_parent
 new_level%parent=>my_parent
 END SUBROUTINE create_timed_routine


 SUBROUTINE timer_start
!------------------------------------------!
! Start the timers for the whole program.  !
!------------------------------------------!
 USE parallel, ONLY : am_master
 IMPLICIT NONE
 if(.not.am_master)return
 allocate(level_zero)
 pt_current=>level_zero
 level_zero%name='LEVEL_ZERO'
 level_zero%time=0.0
 level_zero%rtime=0.0
 level_zero%start=0.0
 level_zero%hits=1
 level_zero%error=.false.
 level_zero%collapse=.false.
 level_zero%first_sublevel=>level_zero
 level_zero%next=>level_zero
 level_zero%parent=>level_zero
 call cputime(start_time)
 start_wtime=walltime()
 END SUBROUTINE timer_start


 SUBROUTINE timer_end
!------------------------------------------!
! Stop the timers for the whole program.   !
!------------------------------------------!
 USE parallel, ONLY : use_timer
 IMPLICIT NONE
 INTEGER level
 INTEGER(i64) hits
 REAL(sp) total_time,total_wtime,cpu_ratio_branch,cpu_ratio_tot,tv1(3),&
  &total_time_min,total_time_hr,total_wtime_min,total_wtime_hr,system_ratio
 LOGICAL errors
 CHARACTER(40) label
 CHARACTER(80) line
 CHARACTER(24) tmpr,tmpr2,tmpr3
 TYPE(timed_routine),POINTER :: pt_current,pt_temp

 line=' '
 total_time=tcputime()-start_time(3)
 total_time_min=total_time/60.0_sp
 total_time_hr=total_time_min/60.0_sp
 total_wtime=walltime()-start_wtime
 total_wtime_min=total_wtime/60.0_sp
 total_wtime_hr=total_wtime_min/60.0_sp
 write(o,*)
 tmpr=r2ss(total_time,'(f13.4)') 
 tmpr2=r2ss(total_time_min,'(f13.4)')
 tmpr3=r2ss(total_time_hr,'(f13.4)')
 if(total_time<60.0_sp)then
  write(o,'(1x,a)')'Total LOUIS CPU time  : : : '//trim(tmpr)//' sec'
 elseif(total_time<3600.0_sp)then
  write(o,'(1x,a)')'Total LOUIS CPU time  : : : '//trim(tmpr)//' sec ('&
   &//trim(tmpr2)//' min)'
 else
  write(o,'(1x,a)')'Total LOUIS CPU time  : : : '//trim(tmpr)//' sec ('&
   &//trim(tmpr2)//' min, '//trim(tmpr3)//' hrs)'
 endif
 tmpr=r2ss(total_wtime,'(f13.4)') 
 tmpr2=r2ss(total_wtime_min,'(f13.4)')
 tmpr3=r2ss(total_wtime_hr,'(f13.4)')
 if(total_wtime<60.0_sp)then
  write(o,'(1x,a)')'Total LOUIS real time : : : '//trim(tmpr)//' sec'
 elseif(total_wtime<3600.0_sp)then
  write(o,'(1x,a)')'Total LOUIS real time : : : '//trim(tmpr)//' sec ('&
   &//trim(tmpr2)//' min)'
 else
  write(o,'(1x,a)')'Total LOUIS real time : : : '//trim(tmpr)//' sec ('&
   &//trim(tmpr2)//' min, '//trim(tmpr3)//' hrs)'
 endif
 write(o,*)

 if(use_timer)then

  level=0 ; pt_current=>level_zero ; errors=.false.
  pt_current%time(3)=total_time
  write(o,*)'Timing information:'
  write(o,*)'----------------------------------CPU---Sys/CPU--CPU/brn--&
   &CPU/tot--Hits--'

  mainprint: do
! Find next label to print
   if(associated(pt_current%first_sublevel,pt_current).or.&
    &(pt_current%collapse.and..not.ignore_collapse))then
! No deeper level in this branch
    if(associated(pt_current%next,pt_current%parent))then
! No next same-level from parent. Print unaccounted-for time for parent, then
! go up the 'parent' hierarchy until a 'next' pointer can be followed.
     do
      if(associated(pt_current%parent,pt_current))exit mainprint
      pt_current=>pt_current%parent ; level=level-1
      if(pt_current%rtime(3)>0.)then
       label=repeat('-',level)//'[rest]'
       system_ratio=0. ; cpu_ratio_tot=0. ; cpu_ratio_branch=0.
       if(pt_current%rtime(3)>0.)system_ratio=min(1.,&
        &max(0.,pt_current%rtime(2)/pt_current%rtime(3)))
       if(total_time>0.)cpu_ratio_tot=min(1.,max(0.,pt_current%rtime(3)/&
        &total_time))
       if(pt_current%time(3)>0.)cpu_ratio_branch=min(1.,max(0.,&
        &pt_current%rtime(3)/pt_current%time(3)))
       call format_timer_output(label,system_ratio,max(0.0,&
        &pt_current%rtime(3)),cpu_ratio_branch,cpu_ratio_tot,-1_i64,line)
       write(o,'(a)')trim(line)
      endif
      if(associated(pt_current%next,pt_current%parent))cycle
      pt_current=>pt_current%next
      exit
     enddo
    else
     pt_current=>pt_current%next
    endif
   else
! Gather total times for children
    tv1(:)=0. ; pt_temp=>pt_current%first_sublevel
    do
     tv1=tv1+pt_temp%time
     if(associated(pt_temp%next,pt_temp%parent))exit
     pt_temp=>pt_temp%next
    enddo
    pt_current%rtime=pt_current%time-tv1
! Move on to next level
    pt_current=>pt_current%first_sublevel ; level=level+1
   endif
   system_ratio=0. ; cpu_ratio_tot=0. ; cpu_ratio_branch=0.
   if(pt_current%time(3)>0.)system_ratio=min(1.,max(0.,pt_current%time(2)/&
    &pt_current%time(3)))
   if(total_time>0.)cpu_ratio_tot=min(1.,max(0.,pt_current%time(3)/total_time))
   pt_temp=>pt_current%parent
   if(pt_temp%time(3)>0.)cpu_ratio_branch=min(1.,max(0.,pt_current%time(3)/&
    &pt_temp%time(3)))
   if(pt_current%collapse.and.level>1)then
    label=repeat('-',level-2)//'+'//trim(pt_current%name)
   else
    label=repeat('-',level-1)//trim(pt_current%name)
   endif
   if(pt_current%error)then
    label=trim(label)//'*' ; errors=.true.
   endif
   hits=pt_current%hits
   if(total_time>0.)cpu_ratio_tot=min(1.,max(0.,pt_current%time(3)/total_time))
   call format_timer_output(label,system_ratio,max(0.0,pt_current%time(3)),&
    &cpu_ratio_branch,cpu_ratio_tot,hits,line)
   write(o,'(a)')trim(line)
  enddo mainprint
  write(o,"(1x,73('-'))")
  write(o,*)'      CPU: CPU time (seconds)     Sys/CPU: system-to-CPU &
   &time ratio'
  write(o,*)'   CPU/brn: CPU time (% of branch)      CPU/tot: CPU time &
   &(% of total)'
  if(errors)write(o,*)'             *: routines without correct timer &
   &deactivation.'
  write(o,"(1x,73('='))")

 else ! .not.use_timer
!  write(o,*)
!  write(o,*)'Subroutine timers deactivated (use TIMING_INFO input keyword)'
!  write(o,*)
  write(o,"(1x,73('='))")
 endif

 END SUBROUTINE timer_end


 SUBROUTINE format_timer_output(label,s_ratio,c_time,c_ratio_branch,&
  &c_ratio_tot,hits,line)
 IMPLICIT NONE
 INTEGER(i64),INTENT(in) :: hits
 REAL(sp),INTENT(in) :: s_ratio,c_time,c_ratio_branch,c_ratio_tot
 CHARACTER(40),INTENT(in) :: label
 CHARACTER(80),INTENT(out) :: line
 INTEGER temp_power,hits10
 INTEGER(i64) temp_hits
 REAL(sp) rhits
 CHARACTER(40) hits_string,suffix
 if(hits>=0_i64)then
! Divide hits by 10^(multiple of 3) to display in human-readable form.
  temp_hits=hits ; temp_power=0
  do
   if(temp_hits>=1000_i64)then
    temp_hits=temp_hits/1000_i64 ; temp_power=temp_power+3
   else
    exit
   endif
  enddo
! Display two digits.
  if(temp_hits<10.and.temp_power>0)then
   rhits=real(hits)*10.**(-temp_power)
   if(nint(rhits)>9)then
    hits_string='10'
   else
    hits10=int(hits/10_i64**(temp_power-1))
    hits_string=trim(i2s(hits10/10))//'.'//trim(i2s(mod(hits10,10)))
   endif
  else
   hits_string=trim(i2s(int(temp_hits)))
  endif
  select case(temp_power)
  case(0) ; suffix=''  ; case(3)  ; suffix='k' ; case(6) ; suffix='M'
  case(9) ; suffix='G' ; case(12) ; suffix='T'
  case default
   suffix='E'//trim(i2s(temp_power))
  end select
 else
  hits_string='' ; suffix=''
 endif
 write(line,"(1x,a,t30,f10.2,t42,f6.2,'%',t51,f6.2,'%',t60,f6.2,'%',t69,a3,a)")&
  &trim(label),c_time,100.*s_ratio,100.*c_ratio_branch,&
  &100.*c_ratio_tot,trim(hits_string),trim(suffix)
 END SUBROUTINE format_timer_output


 LOGICAL FUNCTION exceeds_time_limit(first_block)
!--------------------------------------------------------------------------!
! Check if the calculation has exceeded the imposed time limit, or if it's !
! going to do so in the next block. Broadcast result to slave nodes.       !
!--------------------------------------------------------------------------!
 USE parallel
 USE store, ONLY : max_cpu_time,max_real_time
 IMPLICIT NONE
 LOGICAL,INTENT(in) :: first_block
 REAL(sp) test_cpu(3),test_real
 INTEGER,SAVE :: nblock
 REAL(sp),SAVE :: ave_block_cpu_time,first_cpu_time,ave_block_real_time,&
  &first_real_time
 LOGICAL,SAVE :: first_call=.true.,have_max_cpu,have_max_real,have_max

 exceeds_time_limit=.false.
 if(first_call)then
  have_max_cpu=max_cpu_time>0. ; have_max_real=max_real_time>0.
  have_max=have_max_cpu.or.have_max_real
  first_call=.false.
 endif

 if(am_master)then
  if(have_max)then
   if(first_block)then
    nblock=0
   else
    nblock=nblock+1
   endif
  endif
  if(have_max_cpu)then
   call cputime(test_cpu)
   if(first_block)then
    first_cpu_time=test_cpu(3) ; ave_block_cpu_time=0.
   else
    ave_block_cpu_time=(test_cpu(3)-first_cpu_time)/real(nblock,sp)
   endif
   if(test_cpu(3)-start_time(3)+ave_block_cpu_time*1.2>max_cpu_time)&
    &exceeds_time_limit=.true.
  endif
  if(have_max_real)then
   test_real=walltime()
   if(first_block)then
    first_real_time=test_real ; ave_block_real_time=0.
   else
    ave_block_real_time=(test_real-first_real_time)/real(nblock,sp)
   endif
   if(test_real-start_wtime+ave_block_real_time*1.2>max_real_time)&
    &exceeds_time_limit=.true.
  endif
 endif

 if(have_max)then
  call mpi_bcast(exceeds_time_limit,1,mpi_logical,0,mpi_comm_world,ierror)
  call checkmpi(ierror,'Broadcasting exceeds_time_limit in exceeds_time_limit')
 endif

 END FUNCTION exceeds_time_limit


 SUBROUTINE cputime(time)
!------------------------------------------------------------!
! CPUTIME                                                    !
! -------                                                    !
!                                                            !
! MDT 4.2002                                                 !
!                                                            !
! Returns the amount of cpu time used since some previous    !
! time, in seconds. The previous time is arbitrary, but may  !
! be the time that the job or the program started, for       !
! example. The difference between two separate calls of this !
! routine is then (approximately) the cpu time used between  !
! the calls. The f90 intrinsic date_and_time is used to      !
! calculate wall time - being careful that the result does   !
! overflow the largest available integer (which the f90      !
! system_clock intrinsic will do very easily), while the     !
! common fortran extension etime is used to calculate CPU    !
! time (and its breakdown into user and system bits). On     !
! systems not supporting the etime extension (e.g. NAG) a    !
! separate etime.c routine can be built into LOUIS.          !
!------------------------------------------------------------!
 IMPLICIT NONE
 REAL(sp),INTENT(out) :: time(3)
 REAL(sp) user_and_system_time(2)
 INTERFACE
  REAL(kind(1.0)) FUNCTION etime(user_and_system_time)
   REAL(kind(1.0)),INTENT(inout) :: user_and_system_time(2)
  END FUNCTION etime
 END INTERFACE
 time(3)=etime(user_and_system_time) ! total
 time(1)=user_and_system_time(1)     ! user
 time(2)=user_and_system_time(2)     ! system
 END SUBROUTINE cputime


 REAL(sp) FUNCTION tcputime()
 IMPLICIT NONE
 REAL(sp) user_and_system_time(2)
 REAL(sp),EXTERNAL :: etime
 tcputime=etime(user_and_system_time) ! total only
 END FUNCTION tcputime


 REAL(sp) FUNCTION walltime()
 IMPLICIT NONE
 INTEGER(i64) time
 INTEGER(i64),SAVE :: time_offset
 LOGICAL,SAVE :: first_call=.true.
 call clock_long(time)
 if(first_call)then
  time_offset=time
  first_call=.false.
 endif
 walltime=real(time-time_offset,sp)*inv_wall_rate
 END FUNCTION walltime


 SUBROUTINE clock_long(t)
!------------------------------------------------------!
! Returns the number of milliseconds that have elapsed !
! since 1st January 2000 00:00. The output t has type  !
! long integer (i64) to make sure it can hold a        !
! sufficiently large number.                           !
!------------------------------------------------------!
 IMPLICIT NONE
 INTEGER(i64),INTENT(out) :: t
 INTEGER i,time_values(8)
! 'time_values' are (/ year, month, day, timezone, hour, min, sec, millisec /)
 t=0_i64 ; call date_and_time(values=time_values)
! Year (in days)
 do i=2000,time_values(1)-1
  t=t+days_upto_month(i,12)
 enddo
 t=t+days_upto_month(time_values(1),time_values(2)-1) ! Month (in days)
 t=t+time_values(3)-1_i64 ! Days
 t=t*24_i64 ; t=t+time_values(5) ! Hours
 t=t*60_i64 ; t=t+time_values(6) ! Minutes
 t=t*60_i64 ; t=t+time_values(7) ! Seconds
 t=t*1000_i64 ; t=t+time_values(8) ! Milliseconds
 END SUBROUTINE clock_long


 INTEGER(i64) FUNCTION days_upto_month(year,month)
!----------------------------------------------------------------------------!
! Returns the total number of days that are in the months from January up to !
! "month", in the year "year".                                               !
!----------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: year,month
 INTEGER(i64),PARAMETER :: days_cumm(12)=(/31,59,90,120,151,181,212,243,273,&
  &304,334,365/)
 days_upto_month=0
 if(month<=0.or.month>=13)return
 days_upto_month=days_cumm(month)
! If year is a leap-year, February has 29 days
 if(month>1.and.((mod(year,4)==0.and.mod(year,100)/=0).or.mod(year,400)==0))&
  &days_upto_month=days_upto_month+1_i64
 END FUNCTION days_upto_month


 SUBROUTINE loop_time_estimate(description,iter,num_iters)
!-------------------------------------------------------------------------!
! Print estimated time for completing a given loop. The routine should be !
! called at every iteration of the loop, providing the current iteration  !
! number ITER and total number of iterations NUM_ITERS. The estimate will !
! be printed once the loop has gone beyond the second iteration AND the   !
! code has spent over REPORT_ONSET CPU seconds in the loop.               !
! Additionally, if there is a time limit defined, an error is raised if   !
! the loop will take longer than the remaining time (CPU or wall).        !
! Optionally, a final summary can be printed if the routine is called     !
! after the loop without specifying ITER or NUM_ITERS.                    !
!                                                                         !
! NN 12.2008                                                              !
!-------------------------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: description
 INTEGER,INTENT(in),OPTIONAL :: iter,num_iters
 REAL(sp),SAVE :: ctime1=0.,ctime2=0.,wtime2=0.
 LOGICAL,SAVE :: estimate_reported=.false.
 REAL(sp) ctime_diff,ctime_est,wtime_diff,wtime_est,wt
 REAL(sp),PARAMETER :: report_onset=60.
 if(.not.am_master)return
 if(.not.present(iter))then ! final call
  ctime_diff=tcputime()-ctime1
  write(o,*)trim(description),' [total CPU time: ',&
   &trim(time2human(ctime_diff)),']'
  write(o,*)
 elseif(iter==1)then ! 1st iteration
  ctime1=tcputime() ; estimate_reported=.false.
  write(o,*)trim(description)
 elseif(iter==2)then ! 2nd iteration
  ctime2=tcputime() ; wtime2=walltime()
 elseif(.not.estimate_reported)then ! enough iterations to report
  ctime_diff=tcputime()-ctime2
  if(ctime_diff>=report_onset)then ! enough elapsed time to report
   ctime_est=real(num_iters-iter+1)*ctime_diff/real(iter-2)
   write(o,*)' [CPU time: ',trim(time2human(ctime_diff)),' elapsed, ',&
    &trim(time2human(ctime_est)),' remaining]'
   estimate_reported=.true.
! Check time limits
   if(max_cpu_time>0.)then
    if(max_cpu_time+start_time(3)<tcputime()+ctime_est)call errstop_master&
     &('LOOP_TIME_ESTIMATE','Estimated CPU time exceeds available time.')
   endif
   if(max_real_time>0.)then
    wt=walltime()
    wtime_diff=wt-wtime2
    wtime_est=real(num_iters-iter+1)*wtime_diff/real(iter-2)
    if(max_real_time+start_wtime<wt+wtime_est)call errstop_master&
     &('LOOP_TIME_ESTIMATE','Estimated walltime exceeds available time.')
   endif
  endif ! ctime_diff>=report_onset
 endif ! presence/values of ITER
 END SUBROUTINE loop_time_estimate


 SUBROUTINE traceback
!-----------------------------------------------------------------------!
! Print recorded position in the code at time of encountering an error. !
!-----------------------------------------------------------------------!
 USE parallel,ONLY : use_timer
 IMPLICIT NONE
 INTEGER i
 TYPE(timed_routine),POINTER :: pt_temp
 write(o,*)'LOUIS internal traceback:'
 if(.not.use_timer)then
  if(ilevel_current>0)then
   write(o,*)' Problem detected at ',trim(tb_name(ilevel_current))
  else
   write(o,*)' Problem detected at an unknown routine'
  endif
  do i=ilevel_current-1,1,-1
   write(o,*)' Called from ',trim(tb_name(i))
  enddo
 else
  write(o,*)' Problem detected at ',trim(pt_current%name)
  pt_temp=>pt_current
  do
   if(associated(pt_temp%parent,level_zero))exit
   pt_temp=>pt_temp%parent
   write(o,*)' Called from ',trim(pt_temp%name)
  enddo
 endif
 write(o,*)' Called from MAIN'
 END SUBROUTINE traceback


 SUBROUTINE errstop(routine,error)
!-------------------------------------------------------!
! Write out routine name and error message then stop.   !
!-------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: routine,error
 if(.not.errstop_skip)then
  write(o,*)
  write(o,*)'ERROR : '//trim(routine)
  call wordwrap(trim(error))
  write(o,*)
  call traceback
  write(o,*)
  write(o,'(1x,72(''=''))')
  write(o,*)
  call all_abort
 endif
 END SUBROUTINE errstop


 SUBROUTINE errstop_master(routine,error)
!--------------------------------------------------------!
! Write out routine name and error message (master node) !
! and then stop (all nodes).                             !
!--------------------------------------------------------!
 USE parallel, ONLY : am_master,barrier
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: routine,error
 if(errstop_skip)return
 if(am_master)then
  write(o,*)
  write(o,*)'ERROR : '//trim(routine)
  call wordwrap(trim(error))
  write(o,*)
  call traceback
  write(o,*)
  write(o,'(1x,72(''=''))')
  write(o,*)
 else ! ensure no node stops before master does
  write(o,*)'Hit global error in '//trim(routine)//':'
  call wordwrap(trim(error))
  write(o,*)
  write(o,*)'Waiting for master node to exit.'
  call barrier
 endif
 call all_abort
 END SUBROUTINE errstop_master


 SUBROUTINE errstop2(routine,error,int1)
!--------------------------------------------------------------------------!
! Write out routine name, error message and associated integer then stop.  !
!--------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: int1
 CHARACTER(*),INTENT(in) :: routine,error
 if(.not.errstop_skip)then
  write(o,'(/1x,''ERROR : '',a,/1x,a,1x,a/)')routine,error,trim(i2s(int1))
  write(o,*)
  call traceback
  write(o,*)
  write(o,'(1x,72(''=''))')
  write(o,*)
  call all_abort
 endif
 END SUBROUTINE errstop2


 SUBROUTINE errstop_quiet(routine,error)
!-------------------------------------------------------------------------!
! Write out routine name and error message then stop without whingeing    !
! about typewriters. Often won't work on parallel machines!               !
!-------------------------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: routine,error
 if(.not.errstop_skip)then
  write(o,1)routine,error
1 format(/1x,'ERROR : ',a,/1x,a/)
  call traceback
  write(o,*)
  write(o,'(1x,78(''-''))')
  call all_stop
 endif
 END SUBROUTINE errstop_quiet


 SUBROUTINE errwarn(routine,warning)
!--------------------------!
! Print a warning message. !
!--------------------------!
 USE parallel
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: routine,warning
 if(am_master)then
  call wordwrap('Warning: ['//trim(routine)//'] '//trim(warning))
  write(o,*)
 endif
 END SUBROUTINE errwarn


 SUBROUTINE all_abort
 USE parallel
 IMPLICIT NONE
 close(o)
 call mpi_abort(mpi_comm_world,-1,ierror)
 stop
 END SUBROUTINE all_abort


 SUBROUTINE all_stop
 USE parallel
 IMPLICIT NONE
 close(o)
 call mpi_finalize(ierror)
 stop
 END SUBROUTINE all_stop


 FUNCTION reallocate_rv(p,n)
 REAL(dp),DIMENSION(:),POINTER :: p,reallocate_rv
 INTEGER,INTENT(in) :: n
 INTEGER nold
 allocate(reallocate_rv(n),stat=ialloc)
 if(ialloc/=0)call errstop('REALLOCATE_RV','Allocation problem.')
 if(.not.associated(p))return
 nold=size(p)
 reallocate_rv(1:min(nold,n))=p(1:min(nold,n))
 deallocate(p)
 END FUNCTION reallocate_rv

 FUNCTION reallocate_iv(p,n)
 INTEGER,DIMENSION(:),POINTER :: p,reallocate_iv
 INTEGER,INTENT(in) :: n
 INTEGER nold
 allocate(reallocate_iv(n),stat=ialloc)
 if(ialloc/=0)call errstop('REALLOCATE_IV','Allocation problem.')
 if(.not.associated(p))return
 nold=size(p)
 reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
 deallocate(p)
 END FUNCTION reallocate_iv

 FUNCTION reallocate_hv(p,n)
 CHARACTER(1),DIMENSION(:),POINTER :: p,reallocate_hv
 INTEGER,INTENT(in) :: n
 INTEGER nold
 allocate(reallocate_hv(n),stat=ialloc)
 if(ialloc/=0)call errstop('REALLOCATE_HV','Allocation problem.')
 if(.not.associated(p))return
 nold=size(p)
 reallocate_hv(1:min(nold,n))=p(1:min(nold,n))
 deallocate(p)
 END FUNCTION reallocate_hv

 FUNCTION reallocate_rm(p,n,m)
 REAL(dp),DIMENSION(:,:),POINTER :: p,reallocate_rm
 INTEGER,INTENT(in) :: n,m
 INTEGER nold,mold
 allocate(reallocate_rm(n,m),stat=ialloc)
 if(ialloc/=0)call errstop('REALLOCATE_RM','Allocation problem.')
 if(.not.associated(p))return
 nold=size(p,1)
 mold=size(p,2)
 reallocate_rm(1:min(nold,n),1:min(mold,m))=&
 p(1:min(nold,n),1:min(mold,m))
 deallocate(p)
 END FUNCTION reallocate_rm

 FUNCTION reallocate_im(p,n,m)
 INTEGER,DIMENSION(:,:),POINTER :: p,reallocate_im
 INTEGER,INTENT(in) :: n,m
 INTEGER nold,mold
 allocate(reallocate_im(n,m),stat=ialloc)
 if(ialloc/=0)call errstop('REALLOCATE_IM','Allocation problem.')
 if(.not.associated(p))return
 nold=size(p,1)
 mold=size(p,2)
 reallocate_im(1:min(nold,n),1:min(mold,m))=&
 p(1:min(nold,n),1:min(mold,m))
 deallocate(p)
 END FUNCTION reallocate_im


END MODULE run_control

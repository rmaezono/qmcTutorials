!--------------------------------------------------------------------!
! QUICKBLOCK                                                         !
!                                                                    !
! Neil Drummond, 5/2004                                              !
!                                                                    !
! This is a simple reblocking utility.  It is intended for           !
! analysing very large dmc.hist files.                               !
! If run in a directory containing (d/v)mc.hist, then it will reblock!
! the total energy data in the (d/v)mc.hist file.  Otherwise, it will!
! ask the user for the filename.  If asked to reblock a file other   !
! than (d/v)mc.hist then it will ask the user which column contains  !
! data, whether weighted averages are to be used, and which column   !
! contains the weights.                                              !
!--------------------------------------------------------------------!

MODULE utils
!-------------------------------!
! Miscellaneous utilities, etc. !
!-------------------------------!
 IMPLICIT NONE
 INTEGER,PARAMETER :: dp=kind(1.d0)


CONTAINS


 CHARACTER(12) FUNCTION i2s(n)
!-----------------------------------------------------------------------!
! I2S                                                                   !
! ===                                                                   !
! Convert integers to left justified strings that can be printed in the !
! middle of a sentence without introducing large amounts of white space.!
!                                                                       !
! Calling routine is intended to include something like:                !
! USE utilities                                                         !
! INTEGER i                                                             !
! i=12                                                                  !
! write(*,*)'Integer number ',trim(i2s(i)),' with words at the end.'    !
!-----------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: n
  INTEGER i,j
  INTEGER,PARAMETER :: ichar0=ichar('0')
  i2s=''
  i=abs(n)
  do j=len(i2s),1,-1
   i2s(j:j)=achar(ichar0+mod(i,10))
   i=i/10 ; if(i==0)exit
  enddo ! j
  if(n<0)then
   i2s='-'//adjustl(i2s)
  else
   i2s=adjustl(i2s)
  endif ! n<0
 END FUNCTION i2s


 SUBROUTINE errstop(sub,message)
!---------------------------!
! Report an error and stop. !
!---------------------------!
  IMPLICIT NONE
  CHARACTER(*),INTENT(in) :: sub,message
  write(*,*)
  write(*,*)'ERROR in subroutine '//trim(adjustl(sub))//'.'
  write(*,*)
  call wordwrap(trim(adjustl(message)))
  write(*,*)
  stop
 END SUBROUTINE errstop


 SUBROUTINE wordwrap(text,unit_in,linelength_in)
!-------------------------------------------------------------------------!
! This subroutine prints out the contents of the character string 'text', !
! ensuring that line breaks only occur at space characters.  The output   !
! is written to unit unit_in if this parameter is supplied; otherwise the !
! output is written to unit 6.  The maximum length of each line is given  !
! by linelength_in if this is supplied; otherwise the default line length !
! is 79 characters.                                                       !
!-------------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in),OPTIONAL :: unit_in,linelength_in
  CHARACTER(*),INTENT(in) :: text
  CHARACTER(260) :: temp
  INTEGER :: i,unit,lentext,startpoint,stoppoint,lastpos,linelength
  if(present(unit_in))then
   unit=unit_in
  else
   unit=6
  endif ! unit supplied.
  lentext=len(trim(text))
  if(lentext<1)then
   write(unit,*)
   return
  endif ! No text
  if(present(linelength_in))then
   if(linelength_in>=2)then
    linelength=linelength_in
   else
    linelength=2
   endif ! sensible line-length supplied.
  else
   linelength=79
  endif ! linelength present.
  startpoint=1
  i=0
  do
   i=i+1
   stoppoint=startpoint+linelength-1
   if(stoppoint<=lentext)then
    lastpos=index(trim(text(startpoint:stoppoint))," ",.true.)
    if(lastpos>0)stoppoint=startpoint+lastpos-1
   else
    stoppoint=lentext
   endif ! stoppoint <= length of text
   if(i==1)then
! Allow the user to indent the first line, if (s)he wishes.
    temp=text(startpoint:stoppoint) ! or else pathscale f90 fails to compile
    write(unit,*)trim(temp)
   else
    temp=text(startpoint:stoppoint) ! or else pathscale f90 fails to compile
    write(unit,*)trim(adjustl(temp))
   endif ! i=1
   if(stoppoint==lentext)then
    exit
   else
    startpoint=stoppoint+1
   endif ! Finished text?
  enddo ! Loop over lines.
 END SUBROUTINE wordwrap


END MODULE utils


MODULE analyse_data
!-------------------------------------------------------------!
! Miscellaneous subroutines for reading & analysing the data. !
!-------------------------------------------------------------!
 USE utils
 IMPLICIT NONE
 INTEGER no_pts ! Number of data points
! Arrays with the data points and (optionally) weights.
 REAL(dp),ALLOCATABLE :: data_array(:),weight_array(:)


CONTAINS


 SUBROUTINE get_file(in_file)
!-------------------------------------------------!
! Find out the filename by one method or another. !
!-------------------------------------------------!
  IMPLICIT NONE
  CHARACTER(80),INTENT(out) :: in_file
  INTEGER ierr
  LOGICAL file_present,dmc_exists,vmc_exists
  inquire(file='dmc.hist',exist=dmc_exists)
  inquire(file='vmc.hist',exist=vmc_exists)
! Ask the user for the file to be analysed.
  do
   if(dmc_exists)then
    write(*,*)'Please enter name of data file (press enter to choose dmc.hist).'
   elseif(vmc_exists)then
    write(*,*)'Please enter name of data file (press enter to choose vmc.hist).'
   else
    write(*,*)'Please enter name of data file.'
   endif ! dmc_exists
   read(*,'(a)',iostat=ierr)in_file
   if(ierr/=0)in_file=''
   if(len_trim(in_file)==0)then
    if(dmc_exists)then
     in_file="dmc.hist"
     exit
    elseif(vmc_exists)then
     in_file="vmc.hist"
     exit
    endif ! dmc_exists
   endif ! Nothing entered.
   in_file=adjustl(in_file)
   inquire(file=trim(in_file),exist=file_present)
   if(file_present)then
    write(*,*)
    exit
   else
    write(*,*)'Cannot find this file.  Please give another filename.'
   endif ! file_present
  enddo ! Loop asking for filename
 END SUBROUTINE get_file


 SUBROUTINE read_data(weighted)
!------------------------------------------------------------!
! Prepare the arrays containing the data points and weights. !
!------------------------------------------------------------!
  IMPLICIT NONE
  LOGICAL,INTENT(out) :: weighted
  INTEGER ierr,N,Nequil,wt_or_not,data_column,weight_column,i,skip_lines,&
   &ialloc,yorn
  REAL(dp),ALLOCATABLE :: row(:)
  LOGICAL default
  CHARACTER(1) char1
  CHARACTER(80) in_file

! Find out the name of the data file.
  call get_file(in_file)

! Open the data file.
  open(unit=8,file=trim(in_file),status='old',iostat=ierr)
  if(ierr/=0)call errstop('READ_DATA','Sorry, cannot open file '&
   &//trim(in_file)//'.')

! Count the number of lines in the file.
  N=0
  Nequil=0
  do
   read(8,*,iostat=ierr)char1
   if(ierr>0)call errstop('CHECK_FILE','Error counting data in '&
    &//trim(in_file)//'.')
   if(ierr<0)exit
   if(char1/='#')then
    N=N+1
   elseif(N>0)then
    Nequil=N
   endif
  enddo
  rewind(8)
  if(N<2)call errstop('READ_DATA','Data file should have at least two lines.')

! Establish whether "default" behaviour is to be assumed.
  if(trim(in_file)=='dmc.hist'.or.trim(in_file)=='dmc.hist')then
   do
    write(*,*)'Use the default behaviour for '//trim(in_file)//'?  0=NO, 1=YES.'
    read(*,*,iostat=ierr)yorn
    if(ierr/=0)yorn=-1
    write(*,*)
    if(yorn==1)then
     default=.true. ; exit
    elseif(yorn==0)then
     default=.false. ; exit
    endif ! char1
   enddo
   if(default)then
    if(trim(in_file)=='dmc.hist')then
! Default behaviour for analysing a dmc.hist file.
     write(*,*)'Analysing total energy in dmc.hist.'
     write(*,*)'Will use weights.'
     weighted=.true.
     data_column=4
     weight_column=2
    else
! Default behaviour for analysing a vmc.hist file.
     write(*,*)'Analysing total energy in vmc.hist.'
     weighted=.false.
     data_column=2
     weight_column=-1
    endif ! VMC / DMC
   endif ! default
  else
   default=.false.
  endif ! dmc.hist or vmc.hist.

  if(.not.default)then
! Find out if the user wants to compute weighted averages and find
! out which column of the data file contains the data and which
! contains the weights.
   do
    write(*,*)'Analysing a file called '//trim(in_file)//'.'
    write(*,*)
    write(*,*)'Do you want to use weighted averages?  0=NO, 1=YES.'
    read(*,*,iostat=ierr)wt_or_not
    if(ierr/=0)wt_or_not=-1
    if(wt_or_not==0)then
     weighted=.false.
     exit
    elseif(wt_or_not==1)then
     weighted=.true.
     exit
    endif
   enddo ! Loop asking for weighted averages.
   do
    write(*,*)'Which column of '//trim(in_file)//' contains the data?'
    read(*,*,iostat=ierr)data_column
    if(ierr/=0)data_column=-1
    if(data_column<0)then
     write(*,*)'This should be greater than 0.'
    else
     exit
    endif ! data column<0
   enddo ! Loop asking for data column.
   if(weighted)then
    do
     write(*,*)'Which column of '//trim(in_file)// &
      &' contains the weights?'
     read(*,*,iostat=ierr)weight_column
     if(ierr/=0)weight_column=-1
     if(weight_column<0.or.data_column==weight_column)then
      write(*,*)'Try again.'
     else
      exit
     endif ! weight_column<0 etc
    enddo ! Loop asking for weight column
   endif ! weighted
  endif ! Do we need to ask the user about these things.
  write(*,*)

! Allocate the array for reading in a row.
  if(weighted)then
   allocate(row(max(data_column,weight_column)),stat=ialloc)
  else
   allocate(row(data_column),stat=ialloc)
  endif ! weighted
  if(ialloc/=0)call errstop('READ_DATA','Allocation problem: row.')

! Find out how many lines are to be skipped.
  if(Nequil==0)then
   write(*,*)trim(in_file)//' contains '//trim(i2s(N))//' lines.'
  else
   write(*,*)trim(in_file)//' contains '//trim(i2s(N))//' lines, including '//&
    &trim(i2s(Nequil))//' equilibration lines.'
  endif
  do
   write(*,*)'How many initial lines to discard (-1 for detected &
    &equilibration length)?'
   read(*,*,iostat=ierr)skip_lines
   if(ierr/=0)skip_lines=-2
   if(skip_lines==-1)then
    if(Nequil>N-2)then
     write(*,*)'No data left to analyze if equilibration lines are discarded.'
    else
     skip_lines=Nequil
     exit
    endif
   elseif(skip_lines<0.or.skip_lines>N-2)then
    write(*,*)'Number of lines to skip must be between 0 and '// &
     &trim(i2s(N-2))//'.'
   else
    exit
   endif ! Problem with skip_lines
  enddo ! Loop asking for skip_lines

! Number of points to be evaluated
  no_pts=N-skip_lines

  write(*,*)'Averages will be taken over '//trim(i2s(no_pts))//' data points.'
  write(*,*)

! Allocate the data and weight arrays.
  allocate(data_array(no_pts),stat=ialloc)
  if(ialloc/=0)call errstop('READ_DATA','Allocation problem: data_array.')
  if(weighted)then
   allocate(weight_array(no_pts),stat=ialloc)
   if(ialloc/=0)call errstop('READ_DATA','Allocation problem: weight_array.')
  endif ! weighted

! Read in the data and weights.
  i=0
  do
   if(i>=skip_lines)exit
   read(8,*,iostat=ierr)char1
   if(ierr/=0)call errstop('READ_DATA','Error reading '//trim(in_file)//'.')
   if(char1/='#')i=i+1
  enddo ! i
  i=0
  do
   if(i>=no_pts)exit
   read(8,*,iostat=ierr)row(:)
   if(ierr==0)then
    i=i+1
    data_array(i)=row(data_column)
    if(weighted)then
     weight_array(i)=row(weight_column)
     if(weight_array(i)<=0.d0)call errstop('READ_DATA',&
      &'Found a non-positive weight at line '//trim(i2s(skip_lines+i))//' of '&
      &//trim(in_file)//'.')
    endif ! weighted
   elseif(ierr<0)then
    call errstop('READ_DATA','Unexpectedly reached the end of '//trim(in_file)&
     &//'.')
   endif ! line not a comment.
  enddo ! i

  deallocate(row)

  close(8)

 END SUBROUTINE read_data


 SUBROUTINE compute_stats_weighted
!--------------------------------------------------------------!
! Compute the weighted average of the data, and calculate the  !
! error bar as a function of reblocking transformation number. !
!--------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER i,k,no_rtns,rtn,no_blocks,block_length,no_pts_in_last_block,j,ialloc
  REAL(dp) average,var,std_err,tot_weight, &
   &tot_weight_sq,block_average,red_tot_weight,block_weight, &
   &delta_std_err,eff_no_blocks,bl,ncorr
  REAL(dp),ALLOCATABLE :: std_err_vector(:)

! Compute average of data.  Multiply each data point by its weight.
  average=0.d0
  tot_weight=0.d0
  do i=1,no_pts
   data_array(i)=weight_array(i)*data_array(i)
   average=average+data_array(i)
   tot_weight=tot_weight+weight_array(i)
  enddo ! i
  if(tot_weight<=0.d0)call errstop('COMPUTE_STATS_WEIGHTED',&
   &'Total weight should be positive.')
  average=average/tot_weight

  write(*,*)'Average of data = ',average
  write(*,*)
  write(*,*)'Reblock transformation   ;   Std. error   ;  Error in std. error'

! Number of reblocking transformations
  no_rtns=floor(log(dble(no_pts))/log(2.d0))
  allocate(std_err_vector(no_rtns),stat=ialloc)
  if(ialloc/=0)call errstop('COMPUTE_STATS_WEIGHTED',&
   &'Allocation error: std_err_vector.')

! Initial block length
  block_length=1

  do rtn=0,no_rtns-1

! Number of blocks
   no_blocks=no_pts/block_length

! Evaluate the sum of the squares of the deviations from the average.
! Last, incomplete block has fewer data points and hence a smaller weight.
   var=0.d0
   tot_weight_sq=0.d0
   k=0
   do i=1,no_blocks
    block_average=0.d0
    block_weight=0.d0
    do j=1,block_length
     k=k+1
     block_average=block_average+data_array(k)
     block_weight=block_weight+weight_array(k)
    enddo ! j
    block_average=block_average/block_weight
    var=var+(block_average-average)**2*block_weight
    tot_weight_sq=tot_weight_sq+block_weight**2
   enddo ! i
   block_average=0.d0
   block_weight=0.d0
   no_pts_in_last_block=0
   do
    k=k+1
    if(k>no_pts)exit
    no_pts_in_last_block=no_pts_in_last_block+1
    block_average=block_average+data_array(k)
    block_weight=block_weight+weight_array(k)
   enddo ! k
   if(no_pts_in_last_block>0)then
    block_average=block_average/dble(block_weight)
    var=var+(block_average-average)**2*block_weight
    tot_weight_sq=tot_weight_sq+block_weight**2
   endif ! last block nonzero

! Evaluate variance, standard error in mean and error in standard error.
   red_tot_weight=tot_weight-tot_weight_sq/tot_weight
   var=var/red_tot_weight

   eff_no_blocks=dble(no_blocks)+dble(no_pts_in_last_block)/dble(block_length)

   std_err=sqrt(var/eff_no_blocks)
   delta_std_err=std_err/sqrt(2.d0*(eff_no_blocks-1.d0))
   std_err_vector(rtn+1)=std_err

   write(*,*)rtn,std_err,delta_std_err

! Double the block length for the next reblock.
   block_length=2*block_length

  enddo ! rtn

  write(*,*)

! Analyse reblock plot to suggest reblocked error bar.
  if(std_err_vector(1)>0.d0)then
   bl=1.d0
   do rtn=0,no_rtns-1
    ncorr=(std_err_vector(rtn+1)/std_err_vector(1))**2
    if(bl*bl*bl>=dble(no_pts+no_pts)*ncorr*ncorr)then
     write(*,*)'Correlation-corrected std. error = ',std_err_vector(rtn+1)
     exit
    endif
    bl=bl+bl
   enddo ! rtn
   if(rtn>=no_rtns)write(*,*)'Correlation-corrected std. error =  N/A'
  else
   write(*,*)'Correlation-corrected std. error =  N/A'
  endif
  write(*,*)

  deallocate(data_array,weight_array,std_err_vector)

 END SUBROUTINE compute_stats_weighted


 SUBROUTINE compute_stats_unweighted
!---------------------------------------------------------------!
! Compute the unweighted average of the data, and calculate the !
! error bar as a function of reblocking transformation number.  !
!---------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER i,k,no_rtns,rtn,no_blocks,block_length,no_pts_in_last_block,j,ialloc
  REAL(dp) average,last_block_weight,var,std_err,tot_weight, &
   &tot_weight_sq,block_average,red_tot_weight,rec_block_length, &
   &delta_std_err,bl,ncorr
  REAL(dp),ALLOCATABLE :: std_err_vector(:)

! Compute average of data.
  average=0.d0
  do i=1,no_pts
   average=average+data_array(i)
  enddo ! i
  average=average/dble(no_pts)

  write(*,*)'Average of data = ',average
  write(*,*)
  write(*,*)'Reblock transformation   ;   Std. error   ;  Error in std. error'

! Number of reblocking transformations
  no_rtns=floor(log(dble(no_pts))/log(2.d0))
  allocate(std_err_vector(no_rtns),stat=ialloc)
  if(ialloc/=0)call errstop('COMPUTE_STATS_UNWEIGHTED',&
   &'Allocation error: std_err_vector.')

! Initial block length
  block_length=1

  do rtn=0,no_rtns-1

! Number of blocks
   rec_block_length=1.d0/dble(block_length)
   no_blocks=no_pts/block_length

! Evaluate the sum of the squares of the deviations from the average.
! Weight the last, incomplete block by its size as a fraction of the others.
   var=0.d0
   k=0
   do i=1,no_blocks
    block_average=0.d0
    do j=1,block_length
     k=k+1
     block_average=block_average+data_array(k)
    enddo ! j
    block_average=block_average*rec_block_length
    var=var+(block_average-average)**2
   enddo ! i
   block_average=0.d0
   no_pts_in_last_block=0
   do
    k=k+1
    if(k>no_pts)exit
    no_pts_in_last_block=no_pts_in_last_block+1
    block_average=block_average+data_array(k)
   enddo ! k
   last_block_weight=dble(no_pts_in_last_block)*rec_block_length
   if(no_pts_in_last_block>0)then
    block_average=block_average/dble(no_pts_in_last_block)
    var=var+(block_average-average)**2*last_block_weight
   endif ! last block nonzero

! Evaluate variance, standard error in mean and error in standard error.
   tot_weight=dble(no_blocks)+last_block_weight
   tot_weight_sq=dble(no_blocks)+last_block_weight**2
   red_tot_weight=tot_weight-tot_weight_sq/tot_weight
   var=var/red_tot_weight
   std_err=sqrt(var/tot_weight)
   delta_std_err=std_err/sqrt(2.d0*(tot_weight-1.d0))
   std_err_vector(rtn+1)=std_err

   write(*,*)rtn,std_err,delta_std_err

! Double the block length for the next reblock.
   block_length=2*block_length

  enddo ! rtn

  write(*,*)

! Analyse reblock plot to suggest reblocked error bar.
  if(std_err_vector(1)>0.d0)then
   bl=1.d0
   do rtn=0,no_rtns-1
    ncorr=(std_err_vector(rtn+1)/std_err_vector(1))**2
    if(bl*bl*bl>=dble(no_pts+no_pts)*ncorr*ncorr)then
     write(*,*)'Correlation-corrected std. error = ',std_err_vector(rtn+1)
     exit
    endif
    bl=bl+bl
   enddo ! rtn
   if(rtn>=no_rtns)write(*,*)'Correlation-corrected std. error =  N/A'
  else
   write(*,*)'Correlation-corrected std. error =  N/A'
  endif
  write(*,*)

  deallocate(data_array,std_err_vector)

 END SUBROUTINE compute_stats_unweighted


END MODULE analyse_data


PROGRAM quickblock
!---------------------------!
! Main program starts here. !
!---------------------------!
 USE analyse_data
 IMPLICIT NONE
 LOGICAL weighted

 write(*,*)
 write(*,*)'Quickblock'
 write(*,*)

 call read_data(weighted)

 if(weighted)then
  call compute_stats_weighted
 else
  call compute_stats_unweighted
 endif ! weighted

 write(*,*)'Program finished.'
 write(*,*)

END PROGRAM quickblock

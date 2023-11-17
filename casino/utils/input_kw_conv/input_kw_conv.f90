!--------------------------------------------------------------------------!
! CONV_TO_STD   NDD   05/08/2012                                           !
!                                                                          !
! This program allows the user to convert CASINO input files with the old  !
! keyword set to the new keyword set (and vice versa).                     !
!                                                                          !
! MDT NOTE: the 'old keywords' have now (Jan 2014) been flagged as         !
! 'severely deprecated'. Is is intended to remove support for them in      !
! early 2015, at which point this script is likely to be very useful.      !
!                                                                          !
!--------------------------------------------------------------------------!

MODULE utils
!--------------------------!
! Miscellaneous utilities. !
!--------------------------!
 IMPLICIT NONE
 INTEGER,PARAMETER :: dp=KIND(1.d0)


CONTAINS


 SUBROUTINE capitalize(string,decapitalize_in)
!-----------------------------------------------------------------------------!
! This subroutine converts the lower-case characters in string to upper-case  !
! characters (or vice versa if the optional decapitalize arg. is set to T).   !
!-----------------------------------------------------------------------------!
  IMPLICIT NONE
  LOGICAL,INTENT(in),OPTIONAL :: decapitalize_in
  CHARACTER(*),INTENT(inout) :: string
  INTEGER i,ichar_string
  INTEGER,PARAMETER :: ichara=ichar('a'),icharz=ichar('z'),&
   &icharBA=ichar('A'),icharBZ=ichar('Z'),deltaichar=icharBA-ichara
  LOGICAL decapitalize
  if(present(decapitalize_in))then
   decapitalize=decapitalize_in
  else
   decapitalize=.false.
  endif ! decapitalize_in supplied
  if(decapitalize)then
   do i=1,len_trim(string) ! Upper case -> lower case.
    ichar_string=ichar(string(i:i))
    if(ichar_string>=icharBA.and.ichar_string<=icharBZ) &
     &string(i:i)=achar(ichar_string-deltaichar)
   enddo ! i
  else ! Lower case -> upper case.
   do i=1,len_trim(string)
    ichar_string=ichar(string(i:i))
    if(ichar_string>=ichara.and.ichar_string<=icharz) &
     &string(i:i)=achar(ichar_string+deltaichar)
   enddo ! i
  endif ! decapitalize
 END SUBROUTINE capitalize


 SUBROUTINE errstop(sub,message)
!---------------------------!
! Report an error and stop. !
!---------------------------!
  IMPLICIT NONE
  CHARACTER(*),INTENT(in) :: sub,message
  write(6,*)
  write(6,*)'ERROR in subroutine '//trim(adjustl(sub))//'.'
  write(6,*)
  call wordwrap(trim(adjustl(message)))
  write(6,*)
  stop
 END SUBROUTINE errstop


 CHARACTER(12) FUNCTION i2s(n)
!------------------------------------------------------------------------!
! Convert integers to left justified strings that can be printed in the  !
! middle of a sentence without introducing large amounts of white space. !
!------------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: n
  INTEGER i,j
  INTEGER,PARAMETER :: ichar0=ICHAR('0')
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


 SUBROUTINE wordwrap(text,unit_in,linelength_in)
!---------------------------------------------------------------------------!
! This subroutine prints out the contents of the character string 'text',   !
! ensuring that line breaks only occur at space characters.  The output     !
! is written to unit unit_in if this parameter is supplied; otherwise the   !
! output is written to unit o.  The maximum length of each line is given    !
! by linelength_in if this is supplied; otherwise the default line length   !
! is 79 characters.
!---------------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in),OPTIONAL :: unit_in,linelength_in
  CHARACTER(*),INTENT(in) :: text
  INTEGER i,unit,lentext,startpoint,stoppoint,lastpos,linelength
  CHARACTER(260) temp
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
   i=1+1
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


MODULE change_input
!---------------------------------------------------------------!
! Subroutines for converting the keyword set in the input file. !
!---------------------------------------------------------------!
 USE utils
 IMPLICIT NONE
 PRIVATE
 PUBLIC read_input,change_keywords,write_input,finalise
 INTEGER nline_max,nline,keywordset,nkeywords(2)
 CHARACTER(200),ALLOCATABLE :: inlines(:)
 CHARACTER(50),ALLOCATABLE :: keywordlist(:,:)


CONTAINS


 SUBROUTINE setup_keywordlist
!---------------------------------------------------------!
! Set up the keyword lists (for checking the input file). !
!---------------------------------------------------------!
  IMPLICIT NONE
  INTEGER ialloc

  nkeywords(1)=22 ! Number of old keywords.
  nkeywords(2)=22 ! Number of new keywords.

  allocate(keywordlist(maxval(nkeywords),2),stat=ialloc)
  if(ialloc/=0)call errstop('SETUP_KEYWORDLIST','Allocation error.')

! Old VMC keywords
  keywordlist(1,1)='nequil'
  keywordlist(2,1)='nmove'
  keywordlist(3,1)='nblock'
  keywordlist(4,1)='corper'
  keywordlist(5,1)='nvmcave'
  keywordlist(6,1)='nwrcon'
  keywordlist(7,1)='vmc_twist_av'
  keywordlist(8,1)='nequil_ta'

! Old DMC keywords
  keywordlist(9,1)='nmove_dmc_equil'
  keywordlist(10,1)='nmove_dmc_stats'
  keywordlist(11,1)='nblock_dmc_equil'
  keywordlist(12,1)='nblock_dmc_stats'
  keywordlist(13,1)='nmove_dmcmd_equil'
  keywordlist(14,1)='nmove_dmcmd_stats'
  keywordlist(15,1)='nconfig'
  keywordlist(16,1)='nconfig_prelim'
  keywordlist(17,1)='trip_popn'
  keywordlist(18,1)='corper_dmc'
  keywordlist(19,1)='ndmcave'
  keywordlist(20,1)='num_dmc_twists'
  keywordlist(21,1)='nmove_dmct_equil'
  keywordlist(22,1)='nblock_dmct_equil'

! New VMC keywords
  keywordlist(1,2)='vmc_equil_nstep'
  keywordlist(2,2)='vmc_nstep'
  keywordlist(3,2)='vmc_nblock'
  keywordlist(4,2)='vmc_nconfig_write'
  keywordlist(5,2)='vmc_decorr_period'
  keywordlist(6,2)='vmc_ave_period'
  keywordlist(7,2)='vmc_reequil_nstep'
  keywordlist(8,2)='vmc_ntwist'
  
! New DMC keywords
  keywordlist(9,2)='dmc_equil_nstep'
  keywordlist(10,2)='dmc_equil_nblock'
  keywordlist(11,2)='dmc_stats_nstep'
  keywordlist(12,2)='dmc_stats_nblock'
  keywordlist(13,2)='dmcmd_stats_nstep'
  keywordlist(14,2)='dmcmd_stats_nblock'
  keywordlist(15,2)='dmc_target_weight'
  keywordlist(16,2)='dmc_ave_period'
  keywordlist(17,2)='dmc_decorr_period'
  keywordlist(18,2)='dmc_trip_weight'
  keywordlist(19,2)='dmc_nconf_prelim'
  keywordlist(20,2)='dmc_ntwist'
  keywordlist(21,2)='dmc_reequil_nstep'
  keywordlist(22,2)='dmc_reequil_nblock'

 END SUBROUTINE setup_keywordlist


 SUBROUTINE read_input
!-------------------------!
! Read in the input file. !
!-------------------------!
  IMPLICIT NONE
  INTEGER ierr,i,j,ialloc
  CHARACTER(1) separator
  CHARACTER(50) keyword
  CHARACTER(120) char120

  call setup_keywordlist

  open(unit=8,file='input',status='old',iostat=ierr)
  if(ierr/=0)call errstop('READ_INPUT','Cannot open input.')
  nline=0
  keywordset=0
  do
   read(8,'(a)',iostat=ierr)char120
   if(ierr>0)call errstop('READ_INPUT','Error reading input (1).')
   if(ierr<0)exit
   nline=nline+1
   read(char120,*,iostat=ierr)keyword,separator
   if(ierr==0)then
    keyword=adjustl(keyword)
    if(separator==':')then
     call capitalize(keyword,.true.)
ksets : do j=1,2
      do i=1,nkeywords(j)
       if(trim(keyword)==trim(keywordlist(i,j)))then
        if(keywordset/=j.and.keywordset/=0)call errstop('READ_INPUT', &
         &'Input file seems to contain both old and new keywords.')
        keywordset=j
        exit ksets
       endif ! keyword belong to set j
      enddo ! i
     enddo ksets ! j
    endif ! Separator=":"
   endif ! ierr/=0
  enddo ! lines

  deallocate(keywordlist)

  if(keywordset==1)then
   call wordwrap('The input file uses the old keyword set.  &
    &Will convert to the new keyword set.')
  elseif(keywordset==2)then
   call wordwrap('The input file uses the new keyword set.  &
    &Will convert to the old keyword set.')
  else
   call errstop('READ_INPUT','Unable to determine whether the &
    &input file contains the old or new keyword set.')
  endif ! keywordset
  write(6,*)

  rewind(8)

  nline_max=nline+maxval(nkeywords)
  allocate(inlines(nline_max),stat=ialloc)
  if(ialloc/=0)call errstop('READ_INPUT','Allocation error: inlines.')

  do i=1,nline
   read(8,'(a)',iostat=ierr)inlines(i)
   if(ierr/=0)call errstop('READ_INPUT','Error reading input file (2).')
   inlines(i)=adjustl(inlines(i))
  enddo ! i

  close(8)

 END SUBROUTINE read_input


 SUBROUTINE change_keywords
!--------------------------------------!
! Make the changes to the keyword set. !
!--------------------------------------!
  IMPLICIT NONE
  INTEGER ierr,nnodes,tempi,iline,nvmcave,nblock,nmove,nwrcon,vmc_ntwist,&
   &iline_vmc_ntwist,iline_vmc_twist_av,vmc_nstep,vmc_nconfig_write,&
   &nblock_dmc_equil,nblock_dmct_equil,nblock_dmc_stats,corper_dmc,&
   &ndmcave,nmove_dmc_equil,nmove_dmct_equil,nmove_dmc_stats,nconfig_prelim,&
   &dmc_ave_period,dmc_equil_nstep,dmc_reequil_nstep,dmc_stats_nstep,&
   &dmc_nconf_prelim
  REAL(dp) nconfig,trip_popn,dmc_target_weight,dmc_trip_weight
  LOGICAL vmc_twist_av

  do
   write(6,*)'Please enter the number of cores on which you intend to run.'
   read(5,*,iostat=ierr)nnodes
   if(ierr/=0)nnodes=-1
   if(nnodes<1)then
    write(6,*)'Please try again.  Enter a positive integer.'
   else
    exit
   endif ! nnodes<1
  enddo
  write(6,*)

  if(keywordset==1)then
! OLD -> NEW

! nvmcave -> vmc_ave_period
   call get_keyword_value_int('nvmcave',nvmcave,iline)
   if(iline>0)then
    call replace_keyword_int(iline,'vmc_ave_period',nvmcave,&
     &'#*! vmc.hist reduction factor (Integer)')
   else
    nvmcave=1
   endif ! iline>0

! corper->vmc_decorr_period
   call get_keyword_value_int('corper',tempi,iline)
   if(iline>0)call replace_keyword_int(iline,'vmc_decorr_period',tempi,&
    &'#*! Decorrelation loop length (Integer)')

! nequil->vmc_equil_nstep
   call get_keyword_value_int('nequil',tempi,iline)
   if(iline>0)call replace_keyword_int(iline,'vmc_equil_nstep',tempi,&
    &'#*! Number of equilibration steps (Integer)')

! nequil_ta->vmc_reequil_nstep
   call get_keyword_value_int('nequil_ta',tempi,iline)
   if(iline>0)call replace_keyword_int(iline,'vmc_reequil_nstep',tempi,&
    &'#*! No. of post-twist-change equil (Integer)')

! Insert vmc_twist_av if nec.
   call get_keyword_value_logical('vmc_twist_av',vmc_twist_av,&
    &iline_vmc_twist_av)
   if(iline_vmc_twist_av==0)vmc_twist_av=.false.

   if(vmc_twist_av)then

! nmove,nblock -> vmc_nstep,vmc_nblock
    call get_keyword_value_int('nblock',nblock,iline)
    if(iline>0)then
     call replace_keyword_int(iline,'vmc_nblock',1,&
      &'#*! Number of checkpoints (Integer)')
    else
     nblock=1
    endif ! nblock
    call replace_keyword_int(iline_vmc_twist_av,'vmc_ntwist',nblock,&
     &'#*! No. of VMC twist angles (Integer)')
    call get_keyword_value_int('nmove',nmove,iline)
    if(iline>0)call replace_keyword_int(iline,'vmc_nstep',&
     &nmove*nnodes*nvmcave,'#*! Number of steps (Integer)')

   else

! nmove,nblock -> vmc_nstep,vmc_nblock
    call get_keyword_value_int('nblock',nblock,iline)
    if(iline>0)then
     call replace_keyword_int(iline,'vmc_nblock',nblock,&
      &'#*! Number of checkpoints (Integer)')
    else
     nblock=1
    endif ! nblock
    call get_keyword_value_int('nmove',nmove,iline)
    if(iline>0)call replace_keyword_int(iline,'vmc_nstep', &
     &nmove*nnodes*nblock*nvmcave,'#*! Number of steps (Integer)')

   endif ! vmc_twist_av

! nwrcon->vmc_nconfig_write
   call get_keyword_value_int('nwrcon',nwrcon,iline)
   if(iline>0)call replace_keyword_int(iline,'vmc_nconfig_write',&
    &nnodes*nwrcon,'#*! Number of configs to write (Integer)')

! nmove_dmcmd_equil->dmcmd_equil_nstep
   call get_keyword_value_int('nmove_dmcmd_equil',tempi,iline)
   if(iline>0)call replace_keyword_int(iline,'dmcmd_equil_nstep',tempi,&
    &'#*! Number of equil steps in DMC-MD (Integer)')

! nmove_dmcmd_stats->dmcmd_stats_nstep
   call get_keyword_value_int('nmove_dmcmd_stats',tempi,iline)
   if(iline>0)call replace_keyword_int(iline,'dmcmd_stats_nstep',tempi,&
    &'#*! Number of stats steps in DMC-MD (Integer)')

! num_dmc_twists->dmc_ntwist
   call get_keyword_value_int('num_dmc_twists',tempi,iline)
   if(iline>0)call replace_keyword_int(iline,'dmc_ntwist',tempi,&
    &'#*! No. of DMC twist angles (Integer)')

! nblock_dmc_equil->dmc_equil_nblock
   call get_keyword_value_int('nblock_dmc_equil',nblock_dmc_equil,iline)
   if(iline>0)then
    call replace_keyword_int(iline,'dmc_equil_nblock',&
     &nblock_dmc_equil,'#*! Number of checkpoints (Integer)')
   else
    nblock_dmc_equil=1
   endif ! iline>0

! nblock_dmct_equil->dmc_reequil_nblock
   call get_keyword_value_int('nblock_dmct_equil',nblock_dmct_equil,iline)
   if(iline>0)then
    call replace_keyword_int(iline,'dmc_reequil_nblock',&
     &nblock_dmct_equil,'#*! No. of post twist equil blocks (Integer)')
   else
    nblock_dmct_equil=1
   endif ! iline>0

! nblock_dmc_stats->dmc_stats_nblock
   call get_keyword_value_int('nblock_dmc_stats',nblock_dmc_stats,iline)
   if(iline>0)then
    call replace_keyword_int(iline,'dmc_stats_nblock',&
     &nblock_dmc_stats,'#*! Number of checkpoints (Integer)')
   else
    nblock_dmc_stats=1
   endif ! iline>0

! corper_dmc->dmc_decorr_period
   call get_keyword_value_int('corper_dmc',corper_dmc,iline)
   if(iline>0)then
    call replace_keyword_int(iline,'dmc_decorr_period',&
     &corper_dmc,'#*! DMC correlation period (Integer)')
   else
    corper_dmc=1
   endif ! iline>0

! ndmcave->dmc_ave_period
   call get_keyword_value_int('ndmcave',ndmcave,iline)
   if(iline>0)then
    call replace_keyword_int(iline,'dmc_ave_period',&
     &max(ndmcave/corper_dmc,1),'#*! dmc.hist reduction factor (Integer)')
   else
    ndmcave=1
   endif ! iline>0

! nmove_dmc_equil->dmc_equil_nstep
   call get_keyword_value_int('nmove_dmc_equil',nmove_dmc_equil,iline)
   if(iline>0)then
    dmc_equil_nstep=nblock_dmc_equil*nmove_dmc_equil*ndmcave/corper_dmc
    if(mod(nblock_dmc_equil*nmove_dmc_equil*ndmcave,corper_dmc)>0)&
     &dmc_equil_nstep=dmc_equil_nstep+1
    call replace_keyword_int(iline,'dmc_equil_nstep',dmc_equil_nstep,&
     &'#*! Number of equil steps (Integer)')
   endif ! iline>0

! nmove_dmct_equil->dmc_reequil_nstep
   call get_keyword_value_int('nmove_dmct_equil',nmove_dmct_equil,iline)
   if(iline>0)then
    dmc_reequil_nstep=nblock_dmct_equil*nmove_dmct_equil*ndmcave/corper_dmc
    if(mod(nblock_dmct_equil*nmove_dmct_equil*ndmcave,corper_dmc)>0)&
     &dmc_reequil_nstep=dmc_reequil_nstep+1
    call replace_keyword_int(iline,'dmc_reequil_nstep',&
     &dmc_reequil_nstep,'#*! No. of post twist equil moves (Integer)')
   endif ! iline>0

! nmove_dmc_stats->dmc_stats_nstep
   call get_keyword_value_int('nmove_dmc_stats',nmove_dmc_stats,iline)
   if(iline>0)then
    dmc_stats_nstep=nblock_dmc_stats*nmove_dmc_stats*ndmcave/corper_dmc
    if(mod(nblock_dmc_stats*nmove_dmc_stats*ndmcave,corper_dmc)>0)&
     &dmc_stats_nstep=dmc_stats_nstep+1
    call replace_keyword_int(iline,'dmc_stats_nstep',dmc_stats_nstep,&
     &'#*! Number of stats accum steps (Integer)')
   endif ! iline>0

! nconfig->dmc_target_weight
   call get_keyword_value_real('nconfig',nconfig,iline)
   if(iline>0)call replace_keyword_real(iline,'dmc_target_weight',&
    &nconfig*dble(nnodes),'#*! Total target weight in DMC (Real)')

! trip_popn->dmc_trip_weight
   call get_keyword_value_real('trip_popn',trip_popn,iline)
   if(iline>0)call replace_keyword_real(iline,'dmc_trip_weight',&
    &trip_popn*dble(nnodes),'#*! Catastrophe threshold (Real)')

! nconfig_prelim->dmc_nconf_prelim
   call get_keyword_value_int('nconfig_prelim',nconfig_prelim,iline)
   if(iline>0)call replace_keyword_int(iline,'dmc_nconf_prelim',&
    &nconfig_prelim*nnodes,'#*! No of preliminary DMC configs (Integer)')

  elseif(keywordset==2)then
! NEW -> OLD

! vmc_ave_period -> nvmcave
   call get_keyword_value_int('vmc_ave_period',nvmcave,iline)
   if(iline>0)then
    call replace_keyword_int(iline,'nvmcave',nvmcave,&
     &'#*! # succ VMC pts to average (Integer)')
   else
    nvmcave=1
   endif ! iline

! vmc_decorr_period->corper
   call get_keyword_value_int('vmc_decorr_period',tempi,iline)
   if(iline>0)call replace_keyword_int(iline,'corper',tempi,&
    &'#*! Decorrelation loop length (Integer)')

! vmc_equil_nstep->nequil
   call get_keyword_value_int('vmc_equil_nstep',tempi,iline)
   if(iline>0)call replace_keyword_int(iline,'nequil',tempi,&
    &'#*! Number of equilibration steps (Integer)')

! vmc_reequil_nstep->nequil_ta
   call get_keyword_value_int('vmc_reequil_nstep',tempi,iline)
   if(iline>0)call replace_keyword_int(iline,'nequil_ta',tempi,&
    &'#*! No. of post-twist-change equil (Integer)')

! vmc_ntwist->vmc_twist_av
   call get_keyword_value_int('vmc_ntwist',vmc_ntwist,iline_vmc_ntwist)
   if(iline_vmc_ntwist>0)then
    vmc_twist_av=(vmc_ntwist>0)
   else
    vmc_twist_av=.false.
   endif ! vmc_ntwist

   if(vmc_twist_av)then
      
    call replace_keyword_logical(iline_vmc_ntwist,'vmc_twist_av',.true., &
     &'#*! Perform VMC twist averaging (Boolean)')

! nmove,nblock -> vmc_nstep,vmc_nblock
    call get_keyword_value_int('vmc_nblock',nblock,iline)
    if(iline==0)then
     call insert_keyword_int(iline_vmc_ntwist,'nblock',vmc_ntwist,&
      &'#*! Number of blocks (Integer)')
    else
     if(nblock/=1)call errstop('CHANGE_INPUT',&
      &'Should have vmc_nblock=1 in input.')
     call replace_keyword_int(iline,'nblock',vmc_ntwist,&
      &'#*! No. of VMC twist angles (Integer)')
    endif ! iline=0

    call get_keyword_value_int('vmc_nstep',vmc_nstep,iline)
    if(iline>0)then
     nmove=vmc_nstep/(nnodes*nvmcave)
     if(mod(vmc_nstep,nnodes*nvmcave)>0)nmove=nmove+1
     call replace_keyword_int(iline,'nmove',nmove,&
      &'#*! Number of moves per block (Integer)')
    endif ! iline>0

   else

! nmove,nblock -> vmc_nstep,vmc_nblock
    call get_keyword_value_int('vmc_nblock',nblock,iline)
    if(iline>0)then
     call replace_keyword_int(iline,'nblock',nblock,&
      &'#*! Number of blocks (Integer)')
    else
     nblock=1
    endif ! nblock
    call get_keyword_value_int('vmc_nstep',vmc_nstep,iline)
    if(iline>0)then
     nmove=vmc_nstep/(nnodes*nblock*nvmcave)
     if(mod(vmc_nstep,nnodes*nblock*nvmcave)>0)nmove=nmove+1
     call replace_keyword_int(iline,'nmove',nmove,&
      &'#*! Number of moves per block (Integer)')
    endif ! iline>0

   endif ! vmc_twist_av

! vmc_nconfig_write->nwrcon
   call get_keyword_value_int('vmc_nconfig_write',vmc_nconfig_write,iline)
   if(iline>0)then
    nwrcon=vmc_nconfig_write/nnodes
    if(mod(vmc_nconfig_write,nnodes)>0)nwrcon=nwrcon+1
    call replace_keyword_int(iline,'nwrcon',nwrcon,&
     &'#*! # configs to write (Integer)')
   endif ! iline>0

! dmcmd_equil_nstep->nmove_dmcmd_equil
   call get_keyword_value_int('dmcmd_equil_nstep',tempi,iline)
   if(iline>0)call replace_keyword_int(iline,'nmove_dmcmd_equil',tempi,&
    &'#*! Number of equil moves in DMC-MD (Integer)')

! dmcmd_stats_nstep->nmove_dmcmd_stats
   call get_keyword_value_int('dmcmd_stats_nstep',tempi,iline)
   if(iline>0)call replace_keyword_int(iline,'nmove_dmcmd_stats',tempi,&
    &'#*! Number of stats moves in DMC-MD (Integer)')

! dmc_ntwist->num_dmc_twists
   call get_keyword_value_int('dmc_ntwist',tempi,iline)
   if(iline>0)call replace_keyword_int(iline,'num_dmc_twists',tempi,&
    &'#*! No. of DMC twist angles (Integer)')

! dmc_equil_nblock->nblock_dmc_equil
   call get_keyword_value_int('dmc_equil_nblock',nblock_dmc_equil,iline)
   if(iline>0)then
    call replace_keyword_int(iline,'nblock_dmc_equil',&
     &nblock_dmc_equil,'#*! Number of blocks (Integer)')
   else
    nblock_dmc_equil=1
   endif ! iline>0

! dmc_reequil_nblock->nblock_dmct_equil
   call get_keyword_value_int('dmc_reequil_nblock',nblock_dmct_equil,iline)
   if(iline>0)then
    call replace_keyword_int(iline,'nblock_dmct_equil',&
     &nblock_dmct_equil,'#*! No. of post twist equil blocks (Integer)')
   else
    nblock_dmct_equil=1
   endif ! iline>0

! dmc_stats_nblock->nblock_dmc_stats
   call get_keyword_value_int('dmc_stats_nblock',nblock_dmc_stats,iline)
   if(iline>0)then
    call replace_keyword_int(iline,'nblock_dmc_stats',nblock_dmc_stats,&
     &'#*! Number of blocks (Integer)')
   else
    nblock_dmc_stats=1
   endif ! iline>0

! dmc_decorr_period->corper_dmc
   call get_keyword_value_int('dmc_decorr_period',corper_dmc,iline)
   if(iline>0)then
    call replace_keyword_int(iline,'corper_dmc',corper_dmc,&
     &'#*! DMC correlation period (Integer)')
   else
    corper_dmc=1
   endif ! iline>0

! dmc_ave_period->ndmcave
   call get_keyword_value_int('dmc_ave_period',dmc_ave_period,iline)
   if(iline>0)then
    ndmcave=dmc_ave_period*corper_dmc
    call replace_keyword_int(iline,'ndmcave',ndmcave,&
     &'#*! dmc.hist reduction factor (Integer)')
   else
    ndmcave=corper_dmc
   endif ! iline>0

! dmc_equil_nstep->nmove_dmc_equil
   call get_keyword_value_int('dmc_equil_nstep',dmc_equil_nstep,iline)
   if(iline>0)call replace_keyword_int(iline,'nmove_dmc_equil',&
    &1+(dmc_equil_nstep-1)/(nblock_dmc_equil*max(ndmcave/corper_dmc,1)),&
    &'#*! Number of moves (Integer)')

! dmc_reequil_nstep->nmove_dmct_equil
   call get_keyword_value_int('dmc_reequil_nstep',dmc_reequil_nstep,iline)
   if(iline>0)call replace_keyword_int(iline,'nmove_dmct_equil', &
    &1+(dmc_reequil_nstep-1)/(nblock_dmct_equil*max(ndmcave/corper_dmc,1)),&
    &'#*! No. of post twist equil moves (Integer)')

! dmc_stats_nstep->nmove_dmc_stats
   call get_keyword_value_int('dmc_stats_nstep',dmc_stats_nstep,iline)
   if(iline>0)call replace_keyword_int(iline,'nmove_dmc_stats', &
    &1+(dmc_stats_nstep-1)/(nblock_dmc_stats*max(ndmcave/corper_dmc,1)),&
    &'#*! Number of moves (Integer)')

! dmc_target_weight->nconfig
   call get_keyword_value_real('dmc_target_weight',dmc_target_weight,iline)
   if(iline>0)call replace_keyword_real(iline,'nconfig',&
    &dmc_target_weight/dble(nnodes),'#*! # of configs (DMC) (Real)')

! dmc_trip_weight->trip_popn
   call get_keyword_value_real('dmc_trip_weight',dmc_trip_weight,iline)
   if(iline>0)call replace_keyword_real(iline,'trip_popn',&
    &dmc_trip_weight/dble(nnodes),'#*! Catastrophe threshold (Real)')

! dmc_nconf_prelim->nconfig_prelim
   call get_keyword_value_int('dmc_nconf_prelim',dmc_nconf_prelim,iline)
   if(iline>0)call replace_keyword_int(iline,'nconfig_prelim',&
    &max(dmc_nconf_prelim/nnodes,1), &
    &'#*! No of preliminary DMC configs (Integer)')

  else
   call errstop('CHANGE_KEYWORDS','Keyword-set confusion.')
  endif ! keywordset

 END SUBROUTINE change_keywords


 SUBROUTINE get_keyword_value_int(keywordseek,kvalue,iline)
!---------------------------------------------------!
! Determine the numerical value of a given keyword. !
!---------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(out) :: kvalue,iline
  CHARACTER(*),INTENT(in) :: keywordseek
  INTEGER i,ierr
  CHARACTER(1) separator
  CHARACTER(50) keywordin
  iline=0
  do i=1,nline
   read(inlines(i),*,iostat=ierr)keywordin,separator,kvalue
   if(ierr==0)then
    if(separator==':'.and.trim(keywordin)==keywordseek)then
     iline=i
     return
    endif ! keyword found
   endif ! ierr/=0
  enddo ! i
 END SUBROUTINE get_keyword_value_int


 SUBROUTINE get_keyword_value_logical(keywordseek,kvalue,iline)
!-------------------------------------------------!
! Determine the value of a given logical keyword. !
!-------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(out) :: iline
  LOGICAL,INTENT(out) :: kvalue
  CHARACTER(*),INTENT(in) :: keywordseek
  INTEGER i,ierr
  CHARACTER(1) separator
  CHARACTER(50) keywordin
  iline=0
  do i=1,nline
   read(inlines(i),*,iostat=ierr)keywordin,separator,kvalue
   if(ierr==0)then
    if(separator==':'.and.trim(keywordin)==keywordseek)then
     iline=i
     return
    endif ! keyword found
   endif ! ierr/=0
  enddo ! i
 END SUBROUTINE get_keyword_value_logical


 SUBROUTINE get_keyword_value_real(keywordseek,kvalue,iline)
!----------------------------------------------!
! Determine the value of a given real keyword. !
!----------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(out) :: iline
  REAL(dp),INTENT(out) :: kvalue
  CHARACTER(*),INTENT(in) :: keywordseek
  INTEGER i,ierr
  CHARACTER(1) separator
  CHARACTER(50) keywordin
  iline=0
  do i=1,nline
   read(inlines(i),*,iostat=ierr)keywordin,separator,kvalue
   if(ierr==0)then
    if(separator==':'.and.trim(keywordin)==keywordseek)then
     iline=i
     return
    endif ! keyword found
   endif ! ierr/=0
  enddo ! i
 END SUBROUTINE get_keyword_value_real


 SUBROUTINE replace_keyword_int(iline,keyword,kvalue,kcomment)
!-------------------------------------------------!
! Replace iline with the given keyword and value. !
!-------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: iline,kvalue
  CHARACTER(*),INTENT(in) :: keyword,kcomment
  CHARACTER(15) kvaluestring
  CHARACTER(17) keywordstring
  keywordstring=trim(adjustl(keyword))
  kvaluestring=trim(i2s(kvalue))
  write(inlines(iline),'(a)')keywordstring//' : '//kvaluestring//trim(kcomment)
 END SUBROUTINE replace_keyword_int


 SUBROUTINE replace_keyword_logical(iline,keyword,kvalue,kcomment)
!-------------------------------------------------!
! Replace iline with the given keyword and value. !
!-------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: iline
  LOGICAL,INTENT(in) :: kvalue
  CHARACTER(*),INTENT(in) :: keyword,kcomment
  CHARACTER(15) kvaluestring
  CHARACTER(17) keywordstring
  keywordstring=trim(adjustl(keyword))
  if(kvalue)then
   kvaluestring='T'
  else
   kvaluestring='F'
  endif ! kvalue
  write(inlines(iline),'(a)')keywordstring//' : '//adjustl(kvaluestring)&
   &//trim(kcomment)
 END SUBROUTINE replace_keyword_logical


 SUBROUTINE replace_keyword_real(iline,keyword,kvalue,kcomment)
!-------------------------------------------------!
! Replace iline with the given keyword and value. !
!-------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: iline
  REAL(dp),INTENT(in) :: kvalue
  CHARACTER(*),INTENT(in) :: keyword,kcomment
  CHARACTER(15) kvaluestring
  CHARACTER(17) keywordstring
  keywordstring=trim(adjustl(keyword))
  write(kvaluestring,'(f15.7)')kvalue
  write(inlines(iline),'(a)')keywordstring//' : '//adjustl(kvaluestring)&
   &//trim(kcomment)
 END SUBROUTINE replace_keyword_real


 SUBROUTINE insert_keyword_int(iline,keyword,kvalue,kcomment)
!-----------------------------------------------------------------!
! Insert a new line after iline with the given keyword and value. !
!-----------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: iline,kvalue
  CHARACTER(*),INTENT(in) :: keyword,kcomment
  INTEGER jline
  CHARACTER(15) kvaluestring
  CHARACTER(17) keywordstring
  if(iline<0.or.iline>nline)call errstop('INSERT_KEYWORD_LOGICAL', &
   &'Bug: line '//trim(i2s(iline))//' is out of range.')
  if(nline>=nline_max)call errstop('INSERT_KEYWORD_INT', &
   &'NLINE_MAX parameter exceeded.')
  nline=nline+1
  do jline=nline,iline+2
   inlines(jline)=inlines(jline-1)
  enddo ! jline
  keywordstring=trim(adjustl(keyword))
  kvaluestring=trim(i2s(kvalue))
  write(inlines(iline+1),'(a)')keywordstring//' : '//kvaluestring//&
   &trim(kcomment)
 END SUBROUTINE insert_keyword_int


 SUBROUTINE write_input
!-------------------------------------!
! Write out the converted input file. !
!-------------------------------------!
  IMPLICIT NONE
  INTEGER ierr,i
  open(unit=9,file='input_converted',status='replace',iostat=ierr)
  if(ierr/=0)call errstop('WRITE_INPUT','Unable to open input_converted.')
  do i=1,nline
   write(9,'(a)')trim(inlines(i))
  enddo ! i
  close(9)
  write(6,*)'Have written converted input file to input_converted.'
  write(6,*)
 END SUBROUTINE write_input


 SUBROUTINE finalise
!--------------------!
! Deallocate arrays. !
!--------------------!
  IMPLICIT NONE
  if(allocated(inlines))deallocate(inlines)
 END SUBROUTINE finalise


END MODULE change_input


PROGRAM input_kw_conv
!--------------------------!
! Main program stats here. !
!--------------------------!
 USE change_input,ONLY : read_input,change_keywords,write_input,finalise
 IMPLICIT NONE

 write(6,*)
 write(6,*)'NEW/OLD KEYWORD CONVERTER'
 write(6,*)'=============================='
 write(6,*)

! Read in the existing input file and determine the keyword set used.
 call read_input

! Change the new keyword set to the old set or vice versa.
 call change_keywords

! Write out the converted input file.
 call write_input

! Deallocate arrays.
 call finalise

 write(6,*)'Program finished.'
 write(6,*)

END PROGRAM input_kw_conv

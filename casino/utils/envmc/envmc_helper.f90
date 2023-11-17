PROGRAM envmc_helper
!------------------------------------------------------------!
! Read a CASINO output file, compute and print relevant data !
! for the ENVMC utility to print out.                        !
!------------------------------------------------------------!
IMPLICIT NONE
INTEGER io,ierr,nmove,nblock,nrun,opt_cycles,totvmc,svmc,nfig,iET,iET_alt,itime&
 &,itau,iVar,iEV,iEVee,iEVei,iEVnl,iEK,iEKf,iEKt,iER,iERmp,iERmv,iERend,iEReed,&
 &iERret,iERtot,iEfinal,iEfinal_ct,iEfinal_otfr,iEfinal2,iEfinal2_ct,&
 &iEfinal2_otfr
INTEGER,PARAMETER :: maxrun=10000,sp=selected_real_kind(kind(1.0)),&
 &dp=selected_real_kind(kind(1.d0))
INTEGER nblock_run(maxrun)
REAL(dp) total_time,Efinal,dEfinal,dEfinal_ct,dEfinal_otfr,Efinal2,dEfinal2,&
 &dEfinal2_ct,dEfinal2_otfr
REAL(dp),ALLOCATABLE :: ET(:),dET(:),ET_alt(:),dET_alt(:),time(:),tau(:),&
 &dtau(:),EV(:),dEV(:),EVee(:),dEVee(:),EVei(:),dEVei(:),EVnl(:),dEVnl(:),&
 &EK(:),dEK(:),EKf(:),dEKf(:),EKt(:),dEKt(:),Var(:),dVar(:),ER(:),dER(:),&
 &ERmp(:),dERmp(:),ERmv(:),dERmv(:),ERend(:),dERend(:),EReed(:),dEReed(:),&
 &ERret(:),dERret(:),ERtot(:),dERtot(:),nmove_block(:)
LOGICAL isopt,isdmc,vmc1,vmcf,iteracsrch,varerror,use_blocktime
CHARACTER(80) runtype,iterac_sword

call open_file
call read_header
call read_data_and_print_out


CONTAINS


 SUBROUTINE open_file()
!-----------------------!
! Open the output file. !
!-----------------------!
 IMPLICIT NONE
 CHARACTER(512) out_file
 io=10
 write(0,*)'File name:'
 read(5,*,iostat=ierr)out_file
 if(ierr/=0)then ; write(6,*)'not_found' ; stop ; endif
 open(unit=io,file=out_file,status='old',form='formatted',iostat=ierr)
 if(ierr/=0)then ; write(6,*)'not_found' ; stop ; endif
 write(0,*)'Number of figures in errorbars (0 disables conversion):'
 read(5,*,iostat=ierr)nfig
 if(ierr/=0)then ; write(6,*)'nf_error' ; stop ; endif
 if(nfig<0)then  ; write(6,*)'nf_error' ; stop ; endif
 END SUBROUTINE open_file


 SUBROUTINE read_header()
!-----------------------------------------!
! Parse beginning of file, and gather the !
! values of various important keywords.   !
!-----------------------------------------!
 IMPLICIT NONE
 INTEGER vmc_nstep,vmc_nblock,vmc_ntwist
 REAL(sp) block_time
 LOGICAL old_input,new_input
 CHARACTER(80) f1,f2,fn
 CHARACTER(512) rline
 CHARACTER(2048) warning

! Set appropriate defaults.
 runtype='notfound' ; iterac_sword='' ; iteracsrch=.false.
 isopt=.false. ; vmc1=.false. ; vmcf=.false. ; isdmc=.false.
 nmove=0 ; nblock=0 ; opt_cycles=0 ; use_blocktime=.false.
 varerror=.true.

! Go to general input parameters section.
 do
  read(io,'(a)',iostat=ierr)rline
  if(ierr/=0)then
   write(6,*)'not_valid 1'
   stop
  endif
  if(trim(adjustl(rline))=="General input parameters")exit
! Check for warnings.
  f1=field(1,rline)
  if(trim(adjustl(f1))=='Warning'.or.trim(adjustl(f1))=='Warning:')then
   warning=trim(adjustl(rline))
   do ! read until first blank line
    read(io,'(a)',iostat=ierr)rline
    if(ierr/=0)then
     write(6,*)trim(warning) ; write(6,*)'not_valid 2'
     stop
    endif
    if(trim(adjustl(rline))=='')exit
    warning=trim(warning)//' '//trim(adjustl(rline))
   enddo
   write(6,'(a)')trim(warning)
  endif
 enddo

! Read general parameters.
 do
  read(io,'(a)',iostat=ierr)rline
  if(ierr/=0)then
   write(6,*)'not_valid 3'
   stop
  endif
  f1=field(1,rline)
  if(len(trim(adjustl(f1)))==0)exit ! exit at first blank line
  select case(trim(adjustl(f1)))
   case('RUNTYPE')
    runtype=field(-1,rline)
    select case(trim(adjustl(runtype)))
     case('vmc') ; continue
     case('vmc_opt') ; isopt=.true. ; vmc1=.true.
     case('opt_vmc') ; isopt=.true. ; vmc1=.false.
     case('vmc_dmc') ; isdmc=.true.
     case default ; write(6,*)'not_valid 4' ; stop
    end select
   case('ITERAC')
    fn=field(-1,rline)
    select case(trim(adjustl(fn)))
     case('3') ; iterac_sword='Ewald' ; iteracsrch=.true.
     case('4') ; iterac_sword='MPC' ; iteracsrch=.true.
    end select
   case('INTERACTION')
    fn=field(-1,rline)
    select case(trim(adjustl(fn)))
     case('ewald_mpc') ; iterac_sword='Ewald' ; iteracsrch=.true.
     case('mpc_ewald') ; iterac_sword='MPC' ; iteracsrch=.true.
    end select
  end select
 enddo
 if(trim(adjustl(runtype))=='notfound')then
  write(6,*)'not_valid 5' ; stop
 endif

! Go to specific input parameters section.
 do
  read(io,'(a)',iostat=ierr)rline
  if(ierr/=0)then
   write(6,*)'not_valid 6' ; stop
  endif
  select case(trim(adjustl(rline)))
   case('VMC input parameters')
    if(isopt.or.isdmc)then
     write(6,*)'not_valid 7'
     stop
    endif
    exit
   case('VMC/Variance minimization input parameters',&
    &'VMC/Variance minimization (varmin_linjas) input parameters',&
    &'VMC/Energy minimization input parameters',&
    &'VMC/MAD minimization input parameters',&
    &'VMC/optimization input parameters')
    if(.not.isopt)then
     write(6,*)'not_valid 8'
     stop
    endif
    exit
   case('VMC/DMC input parameters')
    if(.not.isdmc)then
     write(6,*)'not_valid 9'
     stop
    endif
    exit
  end select
 enddo

! Read specific parameters.
 old_input=.false. ; new_input=.false. ; vmc_ntwist=0
 do
  read(io,'(a)',iostat=ierr)rline
  if(ierr/=0)then
   write(6,*)'not_valid 10'
   stop
  endif
  f1=field(1,rline)
  if(len(trim(adjustl(f1)))==0)exit
  select case(trim(adjustl(f1)))
   case('NMOVE')
    old_input=.true.
    fn=field(-1,rline) ; if(trim(fn)=='node)')fn=field(-4,rline)
    read(fn,*)nmove
   case('NBLOCK')
    old_input=.true.
    fn=field(-1,rline)
    read(fn,*)nblock
   case('POSTFIT_VMC')
    fn=field(-1,rline)
    read(fn,*)vmcf
   case('OPT_CYCLES')
    fn=field(-1,rline)
    read(fn,*)opt_cycles
   case('VMC_NSTEP')
    new_input=.true.
    fn=field(-1,rline)
    read(fn,*)vmc_nstep
   case('VMC_NBLOCK')
    new_input=.true.
    fn=field(-1,rline)
    read(fn,*)vmc_nblock
   case('VMC_NTWIST')
    new_input=.true.
    fn=field(-1,rline)
    read(fn,*)vmc_ntwist
   case('BLOCK_TIME')
    fn=field(-1,rline)
    read(fn,*)block_time
    if(block_time>0.)use_blocktime=.true.
  end select
 enddo

 if(.not.old_input.and..not.new_input)then
  write(6,*)'not_valid 11'
  stop
 endif

 if(use_blocktime)then ! No nblock from input params. Count blocks instead.

  nrun=0
  nblock=0
  do
   read(io,'(a)',iostat=ierr)rline
   if(ierr/=0)exit
   f1=field(1,rline)
   f2=field(2,rline)
   if(trim(adjustl(f1))=='Block'.and.trim(adjustl(f2))=='average')&
    &nblock=nblock+1
   if(trim(adjustl(f1))=='FINAL'.and.trim(adjustl(f2))=='RESULT:')then
    nrun=nrun+1
    if(nrun>maxrun)then
     write(6,*)'not_valid 12'
     stop
    endif
    nblock_run(nrun)=nblock
    nblock=0
   endif
  enddo
  rewind(io)

 else ! not use_blocktime

  if(.not.old_input)then
   nblock=vmc_nblock ; nmove=vmc_nstep/vmc_nblock
   if(vmc_ntwist>0)then
    vmc_nblock=1
    nblock=vmc_ntwist
   endif
  endif
  nrun=1 ; nblock_run(1)=nblock

 endif

 if(any(nblock_run(1:nrun)==0))then
  write(6,*)'not_valid 13'
  stop
 endif

! Setup number of runs.
 totvmc=1 ; svmc=1
 if(isopt)then
  totvmc=opt_cycles ; if(vmcf)totvmc=totvmc+1
  if(.not.vmc1)svmc=2
  if(use_blocktime)then
   if(nrun/=totvmc)then
    write(6,*)'not_valid 14'
    stop
   endif
  else
   nblock_run(2:totvmc)=nblock
  endif
 endif

! Print relevant data.
 write(6,*)'TOTVMC ',trim(i2s(totvmc))
 write(6,*)'STARTAT ',trim(i2s(svmc-1))

 END SUBROUTINE read_header


 SUBROUTINE allocate_arrays()
!--------------------------!
! Allocate storage arrays. !
!--------------------------!
 IMPLICIT NONE
 INTEGER ialloc

 ialloc=0

 allocate(ET(nblock),dET(nblock),ET_alt(nblock),dET_alt(nblock),&
  &nmove_block(nblock),time(nblock),tau(nblock),dtau(nblock),Var(nblock),&
  &dVar(nblock),stat=ialloc)
 if(ialloc/=0)then
  write(6,*)'alloc_error 1'
  stop
 endif

 allocate(EV(nblock),dEV(nblock),EVee(nblock),dEVee(nblock),&
  &EVei(nblock),dEVei(nblock),EVnl(nblock),dEVnl(nblock),stat=ialloc)
 if(ialloc/=0)then
  write(6,*)'alloc_error 2'
  stop
 endif

 allocate(EK(nblock),dEK(nblock),EKf(nblock),dEKf(nblock),EKt(nblock),&
  &dEKt(nblock),stat=ialloc)
 if(ialloc/=0)then
  write(6,*)'alloc_error 3'
  stop
 endif

 allocate(ER(nblock),dER(nblock),ERmp(nblock),dERmp(nblock),ERmv(nblock),&
  &dERmv(nblock),ERend(nblock),dERend(nblock),EReed(nblock),dEReed(nblock)&
  &,ERret(nblock),dERret(nblock),ERtot(nblock),dERtot(nblock),stat=ialloc)
 if(ialloc/=0)then
  write(6,*)'alloc_error 4'
  stop
 endif

 END SUBROUTINE allocate_arrays


 SUBROUTINE deallocate_arrays()
!----------------------------!
! Deallocate storage arrays. !
!----------------------------!
 IMPLICIT NONE

 deallocate(ET,dET,ET_alt,dET_alt,nmove_block,time,tau,dtau,Var,dVar,EV,dEV,&
  &EVee,dEVee,EVei,dEVei,EVnl,dEVnl,EK,dEK,EKf,dEKf,EKt,dEKt,ER,dER,ERmp,&
  &dERmp,ERmv,dERmv,ERend,dERend,EReed,dEReed,ERret,dERret,ERtot,dERtot)

 END SUBROUTINE deallocate_arrays


 SUBROUTINE read_data_and_print_out()
!-----------------------------------------!
! Read the energy data in the output file !
! and print it out.                       !
!-----------------------------------------!
 IMPLICIT NONE
 INTEGER svmc,ivmc,totvmc,iblock,irun
 REAL(dp) t1
 LOGICAL ignore_next_variance,unitsdone
 CHARACTER(80) start,finish,f1,f2,f3,f4,fn,units
 CHARACTER(512) rline,rline2
 CHARACTER(2048) warning

 total_time=0.d0 ; ierr=0 ; unitsdone=.false.

! Setup number of runs.
 totvmc=1 ; svmc=1
 if(isopt)then
  totvmc=opt_cycles ; if(vmcf)totvmc=totvmc+1
  if(.not.vmc1)svmc=2
 endif

! Loop over VMC cycles.
 irun=0
 main_loop: do ivmc=svmc,totvmc
  irun=irun+1
  if(irun>1)call deallocate_arrays
  nblock=nblock_run(irun)
  call allocate_arrays
  write(6,*)'VMC ',trim(i2s(ivmc))
  iET=0 ; ET=0.d0 ; dET=0.d0
  iET_alt=0 ; ET_alt=0.d0 ; dET_alt=0.d0
  itime=0 ; time=0.d0
  itau=0 ; tau=0.d0 ; dtau=0.d0
  iVar=0 ; Var=0.d0 ; dVar=0.d0
  iEV=0 ; EV=0.d0 ; dEV=0.d0
  iEVee=0 ; EVee=0.d0 ; dEVee=0.d0
  iEVei=0 ; EVei=0.d0 ; dEVei=0.d0
  iEVnl=0 ; EVnl=0.d0 ; dEVnl=0.d0
  iEK=0 ; EK=0.d0 ; dEK=0.d0
  iEKf=0 ; EKf=0.d0 ; dEKf=0.d0
  iEKt=0 ; EKt=0.d0 ; dEKt=0.d0

  iER=0 ; ER=0.d0 ; dER=0.d0
  iERmp=0 ; ERmp=0.d0 ; dERmp=0.d0
  iERmv=0 ; ERmv=0.d0 ; dERmv=0.d0
  iERend=0 ; ERend=0.d0 ; dERend=0.d0
  iEReed=0 ; EReed=0.d0 ; dEReed=0.d0
  iERret=0 ; ERret=0.d0 ; dERret=0.d0
  iERtot=0 ; ERtot=0.d0 ; dERtot=0.d0

  iEfinal=0 ; Efinal=0.d0 ; dEfinal=0.d0
  iEfinal_ct=0 ; dEfinal_ct=0.d0
  iEfinal_otfr=0 ; dEfinal_otfr=0.d0

  iEfinal2=0 ; Efinal2=0.d0 ; dEfinal=0.d0
  iEfinal2_ct=0 ; dEfinal_ct=0.d0
  iEfinal2_otfr=0 ; dEfinal_otfr=0.d0

  if(.not.use_blocktime)nmove_block(:)=real(nmove,dp)

  iblock=0 ; ignore_next_variance=.false.
! Read up to beginning of VMC cycle.
  finish='*     *     *     *     *     *     *     *     *     *     *     *'
  if(isopt)then
   start='PERFORMING VMC CONFIGURATION-GENERATION CALCULATION No. '//&
    &trim(i2s(ivmc))
   if(ivmc==totvmc)then
    finish='abcdefghijklmnopqrstuvwxyz'
    if(vmcf)start='PERFORMING POST-FIT VMC CALCULATION.'
   endif
  elseif(isdmc)then
   start='PERFORMING A VMC CONFIGURATION-GENERATION CALCULATION.'
  else
   start='PERFORMING A SINGLE VMC CALCULATION.'
   finish='abcdefghijklmnopqrstuvwxyz'
  endif
  do
   read(io,'(a)',iostat=ierr)rline
   if(ierr/=0)then ; write(6,*)'EOF' ; return ; endif
   if(trim(adjustl(rline))==trim(adjustl(start)))exit

! Check for warnings.
   f1=field(1,rline)
   if(trim(adjustl(f1))=='Warning'.or.trim(adjustl(f1))=='Warning:')then
    warning=trim(adjustl(rline))
    do ! read until first blank line
     read(io,'(a)',iostat=ierr)rline
     if(ierr/=0)then
      write(6,'(a)')trim(warning) ; write(6,*)'EOF'
      return
     endif
     if(trim(adjustl(rline))=='')exit
     warning=trim(warning)//' '//trim(adjustl(rline))
    enddo
    write(6,'(a)')trim(warning)
   endif
  enddo

! Start reading data.
  do
   read(io,'(a)',iostat=ierr)rline
   if(ierr/=0)then
    call analyse_data ; write(6,*)'EOF' ; return
   endif
   if(trim(adjustl(rline))==trim(adjustl(finish)))then
    call analyse_data ; cycle main_loop
   endif
   f1=field(1,rline)

   select case(trim(adjustl(f1)))

    case('Warning')
     warning=trim(adjustl(rline))
     do ! read until first blank line
      read(io,'(a)',iostat=ierr)rline
      if(ierr/=0)then
       call analyse_data
       write(6,'(a)')trim(warning)
       write(6,*)'EOF'
       return
      endif
      if(trim(adjustl(rline))=='')exit
      warning=trim(warning)//' '//trim(adjustl(rline))
     enddo
     write(6,'(a)')trim(warning)

    case('Warning:')
     warning=trim(adjustl(rline))
     do ! read until first blank line
      read(io,'(a)',iostat=ierr)rline
      if(ierr/=0)then
       call analyse_data
       write(6,'(a)')trim(warning)
       write(6,*)'EOF'
       return
      endif
      if(trim(adjustl(rline))=='')exit
      warning=trim(warning)//' '//trim(adjustl(rline))
     enddo
     write(6,'(a)')trim(warning)

    case('In')
     f2=field(2,rline) ; if(trim(adjustl(f2))/='block')cycle
     f4=field(4,rline)
     if(trim(adjustl(f4))/=trim(i2s(iblock+1)))then
      write(6,*)'iblock_error' ; return
     endif
     iblock=iblock+1

    case('Correlation')
     f2=field(2,rline) ; if(trim(adjustl(f2))/='time')cycle
     if(.not.increase_index(itau,iblock))cycle
     call parse_xdx(field(-3,rline),field(-1,rline),itau,tau,dtau)

    case('No.')
     if(iblock==0)cycle
     fn=field(-1,rline) ; read(fn,*,iostat=ierr)nmove
     if(ierr/=0)cycle
     nmove_block(iblock)=real(nmove,dp)

    case('Number')
     if(iblock==0)cycle
     fn=field(-1,rline) ; read(fn,*,iostat=ierr)nmove
     if(ierr/=0)cycle
     nmove_block(iblock)=real(nmove,dp)

    case('Block')
     if(unitsdone)cycle
     unitsdone=.true.
     select case(trim(adjustl(rline)))
      case('Block average energies (au per simulation cell)')
       units='au/sim.cell'
      case('Block average energies (au per primitive cell)')
       units='au/prim.cell'
      case('Block average energies (au per particle)')
       units='au/particle'
      case('Block average energies (au)')
       units='au'
      case default
       if(trim(units)=='')unitsdone=.false.
     end select
     if(unitsdone)write(6,*)'UNITS ',trim(adjustl(units))

    case('Total')
     f2=field(2,rline)
     select case(trim(adjustl(f2)))
      case('energy')
       f3=field(3,rline)
       select case(trim(adjustl(f3)))
        case('(relativistic)')
         if(.not.increase_index(iER,iblock))cycle
         read(io,'(a)',iostat=ierr)rline2
         if(ierr/=0)then
          call analyse_data ; write(6,*)'EOF' ; return
         endif
         call parse_xdx(field(-1,rline),field(-1,rline2),iER,ER,dER)
        case default
         if(iteracsrch)then
          f4=field(4,rline)
          if(trim(adjustl(f4))/=trim(adjustl(iterac_sword)))then
           if(.not.increase_index(iET_alt,iblock))cycle
           read(io,'(a)',iostat=ierr)rline2
           if(ierr/=0)then
            call analyse_data ; write(6,*)'EOF' ; return
           endif
           call parse_xdx(field(-1,rline),field(-1,rline2),iET_alt,ET_alt,&
            &dET_alt)
           cycle
          endif
         endif
         if(.not.increase_index(iET,iblock))cycle
         read(io,'(a)',iostat=ierr)rline2
         if(ierr/=0)then
          call analyse_data ; write(6,*)'EOF' ; return
         endif
         call parse_xdx(field(-1,rline),field(-1,rline2),iET,ET,dET)
       end select
      case('relativistic')
       if(.not.increase_index(iERtot,iblock))cycle
       read(io,'(a)',iostat=ierr)rline2
       if(ierr/=0)then
        call analyse_data ; write(6,*)'EOF' ; return
       endif
       call parse_xdx(field(-1,rline),field(-1,rline2),iERtot,ERtot,dERtot)
      case('CASINO')
       f3=field(3,rline) ; if(trim(adjustl(f3))/='CPU')cycle
       fn=field(-1,rline) ; read(fn,*,iostat=ierr)t1
       if(ierr/=0)cycle
       total_time=t1
     end select

    case('Potential')
     if(.not.increase_index(iEV,iblock))cycle
     read(io,'(a)',iostat=ierr)rline2
     if(ierr/=0)then
      call analyse_data ; write(6,*)'EOF' ; return
     endif
     call parse_xdx(field(-1,rline),field(-1,rline2),iEV,EV,dEV)

    case('Kinetic')
     f3=field(3,rline)
     select case(trim(adjustl(f3)))
      case('KEI')
       if(.not.increase_index(iEK,iblock))cycle
       read(io,'(a)',iostat=ierr)rline2
       if(ierr/=0)then
        call analyse_data ; write(6,*)'EOF' ; return
       endif
       call parse_xdx(field(-1,rline),field(-1,rline2),iEK,EK,dEK)
      case('TI')
       if(.not.increase_index(iEKt,iblock))cycle
       read(io,'(a)',iostat=ierr)rline2
       if(ierr/=0)then
        call analyse_data ; write(6,*)'EOF' ; return
       endif
       call parse_xdx(field(-1,rline),field(-1,rline2),iEKt,EKt,dEKt)
      case('FISQ')
       if(.not.increase_index(iEKf,iblock))cycle
       read(io,'(a)',iostat=ierr)rline2
       if(ierr/=0)then
        call analyse_data ; write(6,*)'EOF' ; return
       endif
       call parse_xdx(field(-1,rline),field(-1,rline2),iEKf,EKf,dEKf)
     end select

    case('Mass')
     f2=field(2,rline)
     select case(trim(adjustl(f2)))
      case('polarization')
       if(.not.increase_index(iERmp,iblock))cycle
       read(io,'(a)',iostat=ierr)rline2
       if(ierr/=0)then
        call analyse_data ; write(6,*)'EOF' ; return
       endif
       call parse_xdx(field(-1,rline),field(-1,rline2),iERmp,ERmp,dERmp)
      case('velocity')
       if(.not.increase_index(iERmv,iblock))cycle
       read(io,'(a)',iostat=ierr)rline2
       if(ierr/=0)then
        call analyse_data ; write(6,*)'EOF' ; return
       endif
       call parse_xdx(field(-1,rline),field(-1,rline2),iERmv,ERmv,dERmv)
     end select

    case('Electron-nucleus')
     if(.not.increase_index(iERend,iblock))cycle
     read(io,'(a)',iostat=ierr)rline2
     if(ierr/=0)then
      call analyse_data ; write(6,*)'EOF' ; return
     endif
     call parse_xdx(field(-1,rline),field(-1,rline2),iERend,ERend,dERend)

    case('Electron-electron')
     if(.not.increase_index(iEReed,iblock))cycle
     read(io,'(a)',iostat=ierr)rline2
     if(ierr/=0)then
      call analyse_data ; write(6,*)'EOF' ; return
     endif
     call parse_xdx(field(-1,rline),field(-1,rline2),iEReed,EReed,dEReed)

    case('Retardation')
     if(.not.increase_index(iERret,iblock))cycle
     read(io,'(a)',iostat=ierr)rline2
     if(ierr/=0)then
      call analyse_data ; write(6,*)'EOF' ; return
     endif
     call parse_xdx(field(-1,rline),field(-1,rline2),iERret,ERret,dERret)

    case('Ewald','e-e','e-i','e-n')
     f2=field(2,rline)
     select case(trim(adjustl(f1))//' '//trim(adjustl(f2)))
      case('Ewald e-e','e-e interaction')
       if(trim(adjustl(iterac_sword))=='MPC')cycle
       if(.not.increase_index(iEVee,iblock))cycle
       read(io,'(a)',iostat=ierr)rline2
       if(ierr/=0)then
        call analyse_data ; write(6,*)'EOF' ; return
       endif
       call parse_xdx(field(-1,rline),field(-1,rline2),iEVee,EVee,dEVee)
      case('Ewald e-i','Ewald e-n','e-i int.','e-i interaction',&
       &'e-n interaction')
       f3=field(3,rline) ; f4=field(4,rline)
       if(trim(adjustl(f3))=='(non-local)'.or.&
        &trim(adjustl(f4))=='(non-local)')then
        if(.not.increase_index(iEVnl,iblock))cycle
        read(io,'(a)',iostat=ierr)rline2
        if(ierr/=0)then
         call analyse_data ; write(6,*)'EOF' ; return
        endif
        call parse_xdx(field(-1,rline),field(-1,rline2),iEVnl,EVnl,dEVnl)
       else
        if(.not.increase_index(iEVei,iblock))cycle
        read(io,'(a)',iostat=ierr)rline2
        if(ierr/=0)then
         call analyse_data ; write(6,*)'EOF' ; return
        endif
        call parse_xdx(field(-1,rline),field(-1,rline2),iEVei,EVei,dEVei)
       endif
     end select

    case('MPC')
     if(trim(adjustl(iterac_sword))=='Ewald')cycle
     if(.not.increase_index(iEVee,iblock))cycle
     read(io,'(a)',iostat=ierr)rline2
     if(ierr/=0)then
      call analyse_data ; write(6,*)'EOF' ; return
     endif
     call parse_xdx(field(-1,rline),field(-1,rline2),iEVee,EVee,dEVee)

    case('FINAL') ! New format out file (3/2014)
     do
      read(io,'(a)',iostat=ierr)rline
      if(ierr/=0)then
       write(6,'(a)')trim(warning) ; write(6,*)'EOF' ; return
      endif
      if(trim(adjustl(rline))=='')cycle
      fn=field(-2,rline)
      if(trim(adjustl(fn))=='No')then
       iEfinal=1
       call parse_r(field(1,rline),iEfinal,Efinal)
       call parse_r(field(3,rline),iEfinal,dEfinal)
      endif
      if(trim(adjustl(fn))=='time')then
       f1=field(1,rline)
       if(trim(adjustl(f1))/='Insufficient')then
        iEfinal_ct=1
        call parse_r(field(3,rline),iEfinal_ct,dEfinal_ct)
       endif
      endif
      if(trim(adjustl(fn))=='reblocking')then
       f1=field(1,rline)
       if(trim(adjustl(f1))/='Not'.and.trim(adjustl(f1))/='Insufficient')then
        iEfinal_otfr=1
        call parse_r(field(3,rline),iEfinal_otfr,dEfinal_otfr)
       endif
       exit
      endif
     enddo
     if(iteracsrch)then ! get second energy
      do
       read(io,'(a)',iostat=ierr)rline
       if(ierr/=0)then
        write(6,'(a)')trim(warning) ; write(6,*)'EOF' ; return
       endif
       if(trim(adjustl(rline))=='')cycle
       fn=field(-2,rline)
       if(trim(adjustl(fn))=='No')then
        iEfinal2=1
        call parse_r(field(1,rline),iEfinal2,Efinal2)
        call parse_r(field(3,rline),iEfinal2,dEfinal2)
       endif
       if(trim(adjustl(fn))=='time')then
        f1=field(1,rline)
        if(trim(adjustl(f1))/='Insufficient')then
         iEfinal2_ct=1
         call parse_r(field(3,rline),iEfinal2_ct,dEfinal2_ct)
        endif
       endif
       if(trim(adjustl(fn))=='reblocking')then
        f1=field(1,rline)
        if(trim(adjustl(f1))/='Not'.and.trim(adjustl(f1))/='Insufficient')then
         iEfinal2_otfr=1
         call parse_r(field(3,rline),iEfinal2_otfr,dEfinal2_otfr)
        endif
        exit
       endif
      enddo
     endif

    case('CALC_VARIANCE')
     f2=field(2,rline)
     if(trim(adjustl(f2))=='output:')ignore_next_variance=.true.

    case('Variance','Var.')
     if(.not.increase_index(iVar,iblock))cycle
     if(ignore_next_variance)then
      ignore_next_variance=.false. ; cycle
     endif
     if(varerror)then
      read(io,'(a)',iostat=ierr)rline2
      if(ierr/=0)then
       call analyse_data ; write(6,*)'EOF' ; return
      endif
      call parse_xdx(field(-1,rline),field(-1,rline2),iVar,Var,dVar)
     else
      call parse_x(field(-1,rline),iVar,Var)
     endif

    case('Time')
     if(.not.increase_index(itime,iblock))cycle
     call parse_x(field(-1,rline),itime,time)

   end select

  enddo

 enddo main_loop

 if(isdmc)then ! in this case, total time not read in above loop
  do
   read(io,'(a)',iostat=ierr)rline
   if(ierr/=0)then ; write(6,*)'EOF' ; return ; endif
   f1=field(1,rline) ; if(trim(adjustl(f1))/='Total')cycle
   f2=field(2,rline) ; if(trim(adjustl(f2))/='CASINO')cycle
   f3=field(3,rline) ; if(trim(adjustl(f3))/='CPU')cycle
   fn=field(-1,rline) ; read(fn,*,iostat=ierr)t1 ; if(ierr/=0)cycle
   total_time=t1
   if(total_time/=0.d0)call write_out_short('TOTAL_TIME',total_time)
   write(6,*)'EOF' ; return
  enddo
 endif

 END SUBROUTINE read_data_and_print_out


 SUBROUTINE analyse_data()
!------------------------------------------------!
! Analyse statistics of data gathered so far and !
! print results for current VMC block.           !
!------------------------------------------------!
 IMPLICIT NONE
 REAL(dp) cfactor,t1,t2,ke,dke,v,dv

! Process correlation time.
 cfactor=1.d0
 select case(itau)
  case(-1) ; write(6,*)'TAU NaN' ; case(0) ; write(6,*)'TAU N/A'
  case default
   call reaverage_w(itau,tau,dtau,nmove_block,cfactor,t1,t2)
   if(t1<0.d0)then
    write(6,*)'TAU NaN'
   else
    cfactor=sqrt(t1)
    call write_out('TAU',t1,t2)
   endif
 end select

! Process total energy.
 select case(iET)
  case(-1) ; write(6,*)'E_TOTAL NaN' ; case(0) ; write(6,*)'E_TOTAL N/A'
  case default
   call reaverage_w(iET,ET,dET,nmove_block,cfactor,t1,t2)
   call write_out('E_TOTAL',t1,t2)
 end select

! Process variance.
 select case(iVar)
  case(-1) ; write(6,*)'VAR NaN' ; case(0) ; write(6,*)'VAR N/A'
  case default
   if(varerror)then
    call reaverage_w(iVar,Var,dVar,nmove_block,cfactor,t1,t2)
   else
    call average(iVar,Var,t1,t2)
   endif
   t1=max(t1,0.d0)
   call write_out('VAR',t1,t2)
 end select

! Process alternative total energy (MPC/Ewald).
 select case(iET_alt)
  case(-1) ; write(6,*)'E_ALT NaN' ; case(0) ; write(6,*)'E_ALT N/A'
  case default
   call reaverage_w(iET_alt,ET_alt,dET_alt,nmove_block,cfactor,t1,t2)
   call write_out('E_ALT',t1,t2)
 end select

! Process kinetic energy.
 select case(iEK)
  case(-1) ; write(6,*)'E_KEI NaN' ; case(0) ; write(6,*)'E_KEI N/A'
  case default
   call reaverage_w(iEK,EK,dEK,nmove_block,cfactor,t1,t2)
   call write_out('E_KEI',t1,t2)
   ke=t1 ; dke=t2
 end select
 select case(iEKt)
  case(-1) ; write(6,*)'E_TI NaN' ; case(0) ; write(6,*)'E_TI N/A'
  case default
   call reaverage_w(iEKt,EKt,dEKt,nmove_block,cfactor,t1,t2)
   call write_out('E_TI',t1,t2)
 end select
 select case(iEKf)
  case(-1) ; write(6,*)'E_FISQ NaN' ; case(0) ; write(6,*)'E_FISQ N/A'
  case default
   call reaverage_w(iEKf,EKf,dEKf,nmove_block,cfactor,t1,t2)
   call write_out('E_FISQ',t1,t2)
 end select

! Process potential energy.
 select case(iEV)
  case(-1) ; write(6,*)'E_V NaN' ; case(0) ; write(6,*)'E_V N/A'
  case default
   call reaverage_w(iEV,EV,dEV,nmove_block,cfactor,t1,t2)
   call write_out('E_V',t1,t2)
   v=t1 ; dv=t2
 end select
 select case(iEVee)
  case(-1) ; write(6,*)'E_Vee NaN' ; case(0) ; write(6,*)'E_Vee N/A'
  case default
   call reaverage_w(iEVee,EVee,dEVee,nmove_block,cfactor,t1,t2)
   call write_out('E_Vee',t1,t2)
 end select
 select case(iEVei)
  case(-1) ; write(6,*)'E_Vei NaN' ; case(0) ; write(6,*)'E_Vei N/A'
  case default
   call reaverage_w(iEVei,EVei,dEVei,nmove_block,cfactor,t1,t2)
   call write_out('E_Vei',t1,t2)
 end select
 select case(iEVnl)
  case(-1) ; write(6,*)'E_Vnl NaN' ; case(0) ; write(6,*)'E_Vnl N/A'
  case default
   call reaverage_w(iEVnl,EVnl,dEVnl,nmove_block,cfactor,t1,t2)
   call write_out('E_Vnl',t1,t2)
 end select

! Process virial ratio.
 if(iEK==-1.or.iEV==-1)then
  write(6,*)'VR NaN'
 elseif(iEK==0.or.iEV==0)then
  write(6,*)'VR N/A'
 else
  t1=v/ke
  t2=abs(t1*(dv/v+dke/ke))
  call write_out('VR',t1,t2)
 endif

! Process relativistic terms.
 select case(iER)
  case(-1) ; write(6,*)'E_R NaN' ; case(0) ; write(6,*)'E_R N/A'
  case default
   call reaverage_w(iER,ER,dER,nmove_block,cfactor,t1,t2)
   call write_out('E_R',t1,t2)
 end select
 select case(iERmp)
  case(-1) ; write(6,*)'E_Rmp NaN' ; case(0) ; write(6,*)'E_Rmp N/A'
  case default
   call reaverage_w(iERmp,ERmp,dERmp,nmove_block,cfactor,t1,t2)
   call write_out('E_Rmp',t1,t2)
 end select
 select case(iERmv)
  case(-1) ; write(6,*)'E_Rmv NaN' ; case(0) ; write(6,*)'E_Rmv N/A'
  case default
   call reaverage_w(iERmv,ERmv,dERmv,nmove_block,cfactor,t1,t2)
   call write_out('E_Rmv',t1,t2)
 end select
 select case(iERend)
  case(-1) ; write(6,*)'E_Rend NaN' ; case(0) ; write(6,*)'E_Rend N/A'
  case default
   call reaverage_w(iERend,ERend,dERend,nmove_block,cfactor,t1,t2)
   call write_out('E_Rend',t1,t2)
 end select
 select case(iEReed)
  case(-1) ; write(6,*)'E_Reed NaN' ; case(0) ; write(6,*)'E_Reed N/A'
  case default
   call reaverage_w(iEReed,EReed,dEReed,nmove_block,cfactor,t1,t2)
   call write_out('E_Reed',t1,t2)
 end select
 select case(iERret)
  case(-1) ; write(6,*)'E_Rret NaN' ; case(0) ; write(6,*)'E_Rret N/A'
  case default
   call reaverage_w(iERret,ERret,dERret,nmove_block,cfactor,t1,t2)
   call write_out('E_Rret',t1,t2)
 end select
 select case(iERtot)
  case(-1) ; write(6,*)'E_Rtot NaN' ; case(0) ; write(6,*)'E_Rtot N/A'
  case default
   call reaverage_w(iERtot,ERtot,dERtot,nmove_block,cfactor,t1,t2)
   call write_out('E_Rtot',t1,t2)
 end select

! Process block time.
 select case(itime)
  case(-1) ; write(6,*)'BLOCK_TIME NaN' ; case(0) ; write(6,*)'BLOCK_TIME N/A'
  case default
   call average(itime,time,t1,t2)
   call write_out('BLOCK_TIME',t1,t2)
 end select

! Process total time.
 if(total_time/=0.d0)call write_out_short('TOTAL_TIME',total_time)

! Process final result
 select case(iEfinal)
  case(-1) ; write(6,*)'EFINAL NaN'
  case(0) ; write(6,*)'EFINAL N/A'
  case default ; call write_out('EFINAL',Efinal,dEfinal)
 end select
 select case(iEfinal_ct)
  case(-1) ; write(6,*)'EFINAL_CT NaN'
  case(0) ; write(6,*)'EFINAL_CT N/A'
  case default ; call write_out('EFINAL_CT',Efinal,dEfinal_ct)
 end select
 select case(iEfinal_otfr)
  case(-1) ; write(6,*)'EFINAL_OTFR NaN'
  case(0) ; write(6,*)'EFINAL_OTFR N/A'
  case default ; call write_out('EFINAL_OTFR',Efinal,dEfinal_otfr)
 end select
 select case(iEfinal2)
  case(-1) ; write(6,*)'EFINAL2 NaN'
  case(0) ; write(6,*)'EFINAL2 N/A'
  case default
   call write_out('EFINAL2',Efinal2,dEfinal2)
 end select
 select case(iEfinal2_ct)
  case(-1) ; write(6,*)'EFINAL2_CT NaN'
  case(0) ; write(6,*)'EFINAL2_CT N/A'
  case default ; call write_out('EFINAL2_CT',Efinal2,dEfinal2_ct)
 end select
 select case(iEfinal2_otfr)
  case(-1) ; write(6,*)'EFINAL2_OTFR NaN'
  case(0) ; write(6,*)'EFINAL2_OTFR N/A'
  case default ; call write_out('EFINAL2_OTFR',Efinal,dEfinal2_otfr)
 end select

 END SUBROUTINE analyse_data


 CHARACTER(80) FUNCTION field(n,line)
!--------------------------------------------------------!
! Return the N-th field of string LINE, where the fields !
! are separated by one or more spaces.                   !
! If N is negative, return the |N|-th field of LINE from !
! the end.                                               !
!--------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n
 CHARACTER(512),INTENT(in) :: line
 CHARACTER(512) tline
 INTEGER i,k,absn
 LOGICAL back
 if(n==0)then
  field='' ; return
 endif
 absn=abs(n) ; tline=trim(adjustl(line)) ; back=(n<0)
 do i=1,absn-1
  k=scan(trim(adjustl(tline)),' ',back)
  if(k==0)then
   field='' ; return
  endif
  if(back)then
   tline=trim(adjustl(tline(1:k-1)))
  else
   tline=trim(adjustl(tline(k+1:)))
  endif
 enddo
 k=scan(trim(adjustl(tline)),' ',back)
 if(k==0)then
  field=trim(adjustl(tline))
 elseif(back)then
  field=trim(adjustl(tline(k+1:)))
 else
  field=trim(adjustl(tline(1:k-1)))
 endif
 END FUNCTION field


 CHARACTER(20) FUNCTION i2s(n)
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n
 INTEGER i,j
 CHARACTER tmp,sign
 if(n==0)then
  i2s='0' ; return
 endif
 sign=' ' ; if(n<0)sign='-'
 do i=1,len(i2s)
  i2s(i:i)=' '
 enddo
 i=abs(n)
 do j=1,len(i2s)
  if(i==0)exit
  i2s(j:j)=achar(ichar('0')+mod(i,10)) ; i=i/10
 enddo
 i=1 ; j=len_trim(i2s)
 do
  if(i>=j)exit
  tmp=i2s(j:j)
  i2s(j:j)=i2s(i:i)
  i2s(i:i)=tmp
  i=i+1 ; j=j-1
 enddo
 i2s=trim(sign)//i2s
 END FUNCTION i2s


! SUBROUTINE reaverage(n,vector,dvector,prefactor,mean,std_dev)
!------------------------------------------------------------!
! Given a set of data vector(1:n) and associated error bars  !
! dvector(1:n), compute the average and the global error bar,!
! knowing that each datum was generated by a set of NMOVE    !
! data (with NMOVE constant).                                !
!                                                            !
! MEAN is just the average of block energies. STD_DEV**2 is  !
! given by (<O^2>-<O>^2+(nmove-1)<dO^2>)/(n*nmove-1),        !
! PREFACTOR should be 1.d0 for decorrelated data, or it can  !
! be the correlation time for serially correlated data.      !
!                                                            !
! NO LONGER USED (MDT 5.2014). Replaced with reaverage_w.    !
!                                                            !
!------------------------------------------------------------!
!
! IMPLICIT NONE
! INTEGER,INTENT(in) :: n
! REAL(dp),INTENT(in) :: vector(:),dvector(:),prefactor
! REAL(dp),INTENT(out) :: mean,std_dev
! REAL(dp) sumv,sumv2,sumdv2,invn,invnt1,t1
! mean=0.d0 ; std_dev=0.d0
! if(n<1)then
!  return
! elseif(n==1)then
!  mean=vector(1) ; std_dev=dvector(1)*prefactor
! else
!  invn=1.d0/dble(n) ; invnt1=1.d0/dble(n*nmove-1) ; sumv=sum(vector(1:n))
!  sumv2=sum(vector(1:n)**2) ; sumdv2=sum(dvector(1:n)**2)
!  mean=sumv*invn ; t1=sumv2+dble(nmove-1)*sumdv2 ; t1=t1*invn-mean*mean
!  t1=t1*invnt1 ; if(t1<0.d0)t1=0.d0
!  std_dev=prefactor*sqrt(t1)
! endif
! END SUBROUTINE reaverage


 SUBROUTINE reaverage_w(n,vector,dvector,weight,prefactor,mean,std_dev)
!---------------------------------------------------------------!
! Given a set of data vector(1:n) and associated error bars     !
! dvector(1:n), compute the average and the global error bar,   !
! knowing that each datum was generated by a set of WEIGHT      !
! samples (with WEIGHT potentially different for each datum).   !
!                                                               !
! MEAN is just the average of block energies. STD_DEV**2 is     !
! given by :                                                    !
!                                                               !
! <O^2>_w - <O>_w)^2 +  sum_i weight_i*(weight_i-1) std_dev_i^2 !
!                       --------------------------------------- !
!                                  totmove                      !
! ------------------------------------------------------------- !
!                       totmove - 1                             !
!                                                               !
! where                                                         !
!                                                               !
! totmove = sum_i weight_i                                      !
!                                                               !
! <O^2>_2 = sum_i weight_i data_i^2                             !
!           -----------------------                             !
!                   totmove                                     !
!                                                               !
! <O>_2   = sum_i weight_i data_i                               !
!           ---------------------                               !
!                   totmove                                     !
!                                                               !
! PREFACTOR should be 1.d0 for decorrelated data, or it can     !
! be the correlation time for serially correlated data.         !
!                                                               !
! MDT 5.2014                                                    !
!---------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n
 REAL(dp),INTENT(in) :: vector(:),dvector(:),weight(:),prefactor
 REAL(dp),INTENT(out) :: mean,std_dev
 REAL(dp) sumv,sumv2,sumdv2,invn,invnt1,t1,totmove

 mean=0.d0 ; std_dev=0.d0
 if(n<1)then
  return
 elseif(n==1)then
  mean=vector(1) ; std_dev=dvector(1)*prefactor
 else
  totmove=sum(weight(1:n))
  invn=1.d0/totmove
  invnt1=1.d0/(totmove-1.d0)
  sumv=dot_product(weight(1:n),vector(1:n))
  sumv2=dot_product(weight(1:n),vector(1:n)**2)
  sumdv2=dot_product(weight(1:n)*(weight(1:n)-1.d0),dvector(1:n)**2)
  mean=sumv*invn
  t1=sumv2*invn-mean*mean
  t1=t1+sumdv2*invn
  t1=t1*invnt1
  if(t1<0.d0)t1=0.d0
  std_dev=prefactor*sqrt(t1)
 endif
 END SUBROUTINE reaverage_w


 SUBROUTINE average(n,vector,mean,std_dev)
!------------------------------------------------------------!
! Given a set of data vector(1:n), compute the average and   !
! the error bar.  STD_DEV**2 is given by (<O^2>-<O>^2)/(n-1).!
!------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n
 REAL(dp),INTENT(in) :: vector(:)
 REAL(dp),INTENT(out) :: mean,std_dev
 REAL(dp) sumv,sumv2,invn,invn1,t1
 mean=0.d0 ; std_dev=0.d0
 if(n<1)then
  return
 elseif(n==1)then
  mean=vector(1) ; std_dev=0.d0
 else
  invn=1.d0/dble(n) ; invn1=1.d0/dble(n-1) ; sumv=sum(vector(1:n))
  sumv2=sum(vector(1:n)**2) ; mean=sumv*invn ; t1=sumv2*invn-mean*mean
  t1=t1*invn1 ; if(t1<0.d0)t1=0.d0
  std_dev=sqrt(t1)
 endif
 END SUBROUTINE average


 LOGICAL FUNCTION increase_index(indx,iblk)
 IMPLICIT NONE
 INTEGER,INTENT(in) :: indx,iblk
 if(indx==-1)then
  increase_index=.false.
 else
  increase_index=indx+1==iblk
 endif
 END FUNCTION increase_index


 SUBROUTINE write_out(label,x,dx)
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: x,dx
 CHARACTER(*),INTENT(in) :: label
 CHARACTER(80+2*nfig) char_80
 char_80=epmd2ed(x,dx)
 write(6,'(a)')trim(label)//' '//trim(char_80)
 END SUBROUTINE write_out


 SUBROUTINE write_out_short(label,x)
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: x
 CHARACTER(*),INTENT(in) :: label
 CHARACTER(80) char_80
 write(char_80,'(f21.4)')x
 write(6,'(a)')trim(label)//' '//trim(char_80)
 END SUBROUTINE write_out_short


 CHARACTER(80+2*nfig) FUNCTION epmd2ed(x,dx)
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: x,dx
 CHARACTER(80+nfig) c1,c2,c3
 INTEGER oom
 REAL(dp) resc_factor,resc_x,resc_dx
 if(dx==0.d0)then
  if(x==0.d0)then
   epmd2ed="0"
  else
   write(c1,'(f21.12)')x
   epmd2ed=trim(adjustl(c1))
  endif
  return
 endif
 if(nfig==0)then
  write(c1,'(f21.12)')x
  write(c2,'(f21.12)')dx
  epmd2ed=trim(adjustl(c1))//' +/- '//trim(adjustl(c2))
  return
 endif
 oom=nfig-1+ceiling(-log10(dx))
 resc_factor=10.d0**oom
 resc_dx=anint(resc_factor*dx)
 if(abs(resc_dx/10.d0**nfig-1.d0)<1.d-5)then
  oom=oom-1
  resc_factor=resc_factor*0.1d0
  resc_dx=resc_dx*0.1d0
 endif
 if(oom<0)resc_dx=resc_dx*10.d0**(-oom)
 resc_x=anint(resc_factor*x)/resc_factor
 if(oom>0)then
  c2='(f'//trim(i2s(oom+9))//'.'//trim(i2s(oom))//')'
  write(c1,c2)resc_x
 else
  c1=trim(i2s(int(resc_x)))
 endif
 write(c3,'(f'//trim(i2s(nfig+9))//'.0)')resc_dx
 c3=adjustl(c3)
 c3=c3(1:len_trim(c3)-1)
! c3=i2s(resc_dx)
 epmd2ed=trim(adjustl(c1))//'('//trim(c3)//')'
 END FUNCTION epmd2ed


 SUBROUTINE parse_r(fx,i,x)
!--------------------------------------------------------!
! Parse text field fx into x, with i returning unchanged !
! if successful, or as -1 if not.                        !
!--------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: fx
 INTEGER,INTENT(inout) :: i
 REAL(dp),INTENT(inout) :: x
 INTEGER ierr
 i=1
 if(trim(fx)=='NaN')then
  i=-1
  return
 endif
 read(fx,*,iostat=ierr)x
 if(ierr/=0)then
  i=-1
  return
 endif
 END SUBROUTINE parse_r


 SUBROUTINE parse_x(fx,i,x)
!----------------------------------------------------------!
! Parse text field fx into x(i+1), with i returning as i+1 !
! if successful, or as -1 if not.                          !
!----------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: fx
 INTEGER,INTENT(inout) :: i
 REAL(dp),INTENT(inout) :: x(*)
 INTEGER ierr
 i=i+1
 if(trim(fx)=='NaN')then
  i=-1
  return
 endif
 read(fx,*,iostat=ierr)x(i)
 if(ierr/=0)then
  i=-1
  return
 endif
 END SUBROUTINE parse_x


 SUBROUTINE parse_xdx(fx,fdx,i,x,dx)
!------------------------------------------------------------!
! Parse text fields fx and fdx into x(i+1) and dx(i+1), with !
! i returning as i+1 if successful, or as -1 if not.         !
!------------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: fx,fdx
 INTEGER,INTENT(inout) :: i
 REAL(dp),INTENT(inout) :: x(*),dx(*)
 INTEGER ierr
 i=i+1
 if(trim(fx)=='NaN')then
  i=-1
  return
 endif
 read(fx,*,iostat=ierr)x(i)
 if(ierr/=0)then
  i=-1
  return
 endif
 if(len_trim(fdx)==0)then
  dx(i)=0.d0
 else
  if(trim(fdx)=='NaN')then
   i=-1
   return
  endif
  read(fdx,*,iostat=ierr)dx(i)
  if(ierr/=0)then
   i=-1
   return
  endif
 endif
 END SUBROUTINE parse_xdx


END PROGRAM envmc_helper

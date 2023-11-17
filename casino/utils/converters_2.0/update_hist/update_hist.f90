!----------------------------------------------------------------------------!
! UPDATE_HIST                                                                !
! ===========                                                                !
!                                                                            !
! Neil Drummond, 8.2005                                                      !
!                                                                            !
! This utility (i) converts the vmc.hist file produced by CASINO 1.xxx into  !
! the qmc.hist format used by CASINO 2.xxx, and (ii) converts CASINO 1.xxx   !
! dmc.hist and dmc.hist2 files into CASINO 2.xxx qmc.hist format.            !
!----------------------------------------------------------------------------!

MODULE utils
 IMPLICIT NONE
! Tolerance for energy check etc.
 DOUBLE PRECISION,PARAMETER :: tol=1.d-6


 CONTAINS


 SUBROUTINE convert_vmc
!------------------------------------------------------------!
! This subroutine converts vmc.hist (old) -> vmc.hist (new). !
!------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER i,ierr,no_cols_vmc_hist,vmc_block_length,dum_i,ned,neu,btype, &
   &netot,nbasis,npcells,nhu,nhd,iterac,periodicity,no_data_cols
  DOUBLE PRECISION eionion,eionion_temp,eionion_exact,etot_temp,dum_real1, &
   &dum_real2
  DOUBLE PRECISION,ALLOCATABLE :: row_data(:)
  LOGICAL mpc,vcpp,relativistic,electron_gas,unsure,awfn,bwfn,gwfn,pwfn,ewald
  LOGICAL :: eionion_const=.true.

  open(unit=8,file='vmc.hist',status='old',iostat=ierr)
  if(ierr/=0)then
   write(6,*)' Problem opening vmc.hist.'
   stop
  endif

! Read information in header of old vmc.hist file.
  read(8,*,iostat=ierr)
  if(ierr/=0)then
   write(6,*)' Problem reading vmc.hist.'
   stop
  endif
  read(8,*,iostat=ierr)dum_i,nbasis,neu,ned,periodicity,no_cols_vmc_hist, &
   &eionion_exact
  if(ierr/=0)then
   write(6,*)' Problem reading vmc.hist.'
   stop
  endif
  if(no_cols_vmc_hist==15)then
   write(6,*)' Relativistic terms are not present in vmc.hist.'
   relativistic=.false.
  elseif(no_cols_vmc_hist==21)then
   write(6,*)' Relativistic terms are present in vmc.hist.'
   relativistic=.true.
  else
   write(6,*)' Sorry, format of vmc.hist file is wrong.'
   write(6,*)' No. of columns in vmc.hist: '//trim(i2s(no_cols_vmc_hist))//'.'
   stop
  endif
  write(6,*)
  if(nbasis<0.or.neu<0.or.ned<0.or.periodicity<0.or.periodicity>3)then
   write(6,*)'Information in header of old vmc.hist doesn''t make sense.'
   stop
  endif

! Electron(-hole) gas?
  read(8,*,iostat=ierr)dum_real1,dum_real2
  if(ierr/=0)then
   electron_gas=.false.
   nhu=0 ; nhd=0
   write(6,*)' Your system is not an electron(-hole) gas.'
  else
   if((abs(dum_real1)<tol.or.abs(dum_real1-anint(dum_real1)) &
    &<tol*abs(dum_real1)).and.(abs(dum_real2)<tol &
    &.or.abs(dum_real2-anint(dum_real2))<tol*abs(dum_real2)))then
    electron_gas=.true.
    nhu=nint(dum_real1)
    nhd=nint(dum_real2)
    write(6,*)' Your system is an electron(-hole) gas.'
   else
    electron_gas=.false.
    nhu=0 ; nhd=0
    write(6,*)' Your system is not an electron(-hole) gas.'
   endif ! ierr
  endif ! ierr/=0
  write(6,*)
  rewind(8)

  allocate(row_data(1:no_cols_vmc_hist),stat=ierr)
  if(ierr/=0)then
   write(6,*)' Problem allocating row_data.'
   stop
  endif

  read(8,*,iostat=ierr)
  if(ierr/=0)then
   write(6,*)' Problem reading vmc.hist.'
   stop
  endif
  read(8,*,iostat=ierr)vmc_block_length
  if(ierr/=0)then
   write(6,*)' Problem reading block length in vmc.hist.'
   stop
  endif
  if(vmc_block_length<1)then
   write(6,*)' VMC block length is '//trim(i2s(vmc_block_length)) &
    &//', which is less than 1.'
   stop
  endif
  if(electron_gas)then
   read(8,*,iostat=ierr)
   if(ierr/=0)then
    write(6,*)' Problem reading vmc.hist.'
    stop
   endif
  endif ! electron gas
  read(8,fmt="(4(1PE20.13))",iostat=ierr)row_data(1:no_cols_vmc_hist)
  if(ierr/=0)then
   write(6,*)' Problem reading vmc.hist.'
   stop
  endif

  if(row_data(12)/=0.d0.or.row_data(13)/=0.d0.or.row_data(14)/=0.d0)then
   write(6,*)' Core-polarization pseudopotential data are present in vmc.hist.'
   vcpp=.true.
  else
   write(6,*)' Core-polarization pseudopotential data are not present in &
    &vmc.hist.'
   vcpp=.false.
  endif ! VCPP present.
  write(6,*)

  if(row_data(10)/=0.d0.or.row_data(11)/=0.d0)then
   write(6,*)' MPC data are present in vmc.hist.'
   mpc=.true.
  else
   write(6,*)' MPC data are not present in vmc.hist.'
   mpc=.false.
  endif ! MPC present.
  write(6,*)
  if(row_data(6)/=0.d0)then
   write(6,*)' Ewald data are present in vmc.hist.'
   ewald=.true.
  else
   write(6,*)' Ewald data are not present in vmc.hist.'
   ewald=.false.
  endif ! Ewald present.
  write(6,*)
  rewind(8)

! Work out total number of particles.
  netot=neu+ned+nhu+nhd
  if(netot/=1)then
   write(6,*)' There are '//trim(i2s(netot))//' particles in the simulation.'
  else
   write(6,*)' There is 1 particle in the simulation.'
  endif ! singular/plural.
  write(6,*)

! Work out btype value.
  if(electron_gas)then
   btype=0
  else
! Look for xwfn.data files.  If none exist or more than one exist then
! ask the user to supply btype.
   inquire(file='pwfn.data',exist=pwfn)
   inquire(file='bwfn.data',exist=bwfn)
   inquire(file='gwfn.data',exist=gwfn)
   inquire(file='awfn.data',exist=awfn)
   unsure=.not.(pwfn.or.bwfn.or.gwfn.or.awfn).or.(pwfn.and.bwfn) &
    &.or.(pwfn.and.gwfn).or.(pwfn.and.awfn).or.(bwfn.and.gwfn) &
    &.or.(bwfn.and.awfn).or.(gwfn.and.awfn)
   if(unsure)then
    do
     write(6,*)' Please enter the basis type (btype input parameter)'
     read(5,*,iostat=ierr)btype
     if(ierr/=0)btype=-1
     if(btype<1.or.btype>4)then
      write(6,*)' Should have 1<=btype<=4'
      write(6,*)
     else
      exit
     endif
    enddo
    write(6,*)
   else
    if(pwfn)btype=1
    if(bwfn)btype=4
    if(gwfn)btype=2
    if(awfn)btype=3
   endif ! unsure
  endif ! electron_gas
  write(6,*)' The btype input parameter is '//trim(i2s(btype))//'.'
  write(6,*)

! Work out iterac value.
  if(periodicity>0)then
   if(ewald.and..not.mpc)then
    iterac=1
   elseif(mpc.and..not.ewald)then
    iterac=2
   else
    do
     write(6,*)' Please enter the value of the iterac input parameter.'
     read(5,*,iostat=ierr)iterac
     if(ierr/=0)iterac=-1
     if(iterac<1.or.iterac>4)then
      write(6,*)' iterac must be between 1 and 4.'
      write(6,*)
     else
      exit
     endif
    enddo ! Get iterac value.
    write(6,*)
   endif
  else
   iterac=1
  endif ! periodicity
  write(6,*)' The iterac input parameter is '//trim(i2s(iterac))//'.'
  write(6,*)

  if(periodicity>0)then
   do
    write(6,*)' Please enter the number of primitive cells in your simulation.'
    read(5,*,iostat=ierr)npcells
    if(ierr/=0)npcells=-1
    if(npcells<1)then
     write(6,*)' Please try again.'
     write(6,*)
    else
     exit
    endif
   enddo
   write(6,*)
  else
   npcells=1
  endif
  if(npcells/=1)then
   write(6,*)' There are '//trim(i2s(npcells))//' primitive cells.'
  else
   write(6,*)' There is 1 primitive cell.'
  endif ! singular / plural.
  write(6,*)

  open(unit=10,file='vmc.hist.new',status='replace',iostat=ierr)
  if(ierr/=0)then
   write(6,*)' Problem opening vmc.hist.new.  Stopping.'
   stop
  endif ! ierr

! Produce header for top of vmc.hist.new.
  write(10,'(a)')'# Title'
  write(10,'(a)')'# vmc.hist file converted from old format.'
  write(10,'(a)')'# File version number'
  write(10,'(a)')'# 1'
  write(10,'(a)')'# QMC method (VMC, DMC, PIMC, AFMC or RMC)'
  write(10,'(a)')'# VMC'
  write(10,'(a)')'# Electron-electron interaction type (iterac keyword)'
  write(10,'(a)')'# '//trim(i2s(iterac))
  write(10,'(a)')'# Constant (ion-ion) energy'
  write(10,'(a,es26.16)')'# ',eionion_exact
  write(10,'(a)')'# Number of electrons (and other particles) in simulation'
  write(10,'(a)')'# '//trim(i2s(netot))
  write(10,'(a)')'# Number of atoms per primitive cell'
  write(10,'(a)')'# '//trim(i2s(nbasis))
  write(10,'(a)')'# Number of primitive cells'
  write(10,'(a)')'# '//trim(i2s(npcells))
  write(10,'(a)')'# Basis type (btype keyword)'
  write(10,'(a)')'# '//trim(i2s(btype))
  write(10,'(a)')'# Periodic system? (0=NO; 1=YES)'
  if(periodicity==0)then
   write(10,'(a)')'# 0'
  else
   write(10,'(a)')'# 1'
  endif
  write(10,'(a)')'# Number of data columns (excluding iteration number)'
  no_data_cols=6
  if(mpc)no_data_cols=no_data_cols+2
  if(vcpp)no_data_cols=no_data_cols+3
  if(relativistic)no_data_cols=no_data_cols+5
  write(10,'(a)')'# '//trim(i2s(no_data_cols))
  write(10,'(a)')'# Data items (in order)'
  write(10,'(a)')'# ETOT'
  write(10,'(a)')'# KEI'
  write(10,'(a)')'# TI'
  write(10,'(a)')'# EWALD'
  write(10,'(a)')'# LOCAL'
  write(10,'(a)')'# NONLOCAL'
  if(mpc)then
   write(10,'(a)')'# SHORT'
   write(10,'(a)')'# LONG'
  endif ! mpc
  if(vcpp)then
   write(10,'(a)')'# CPPEI'
   write(10,'(a)')'# CPPE'
   write(10,'(a)')'# CPPEE'
  endif ! vcpp
  if(relativistic)then
   write(10,'(a)')'# MASSPOL'
   write(10,'(a)')'# MASSVEL'
   write(10,'(a)')'# DARWINEN'
   write(10,'(a)')'# DARWINEE'
   write(10,'(a)')'# RETARD'
  endif ! relativistic.
  write(10,'(a)')'# Raw QMC data.'

! Read data from vmc.hist and write to vmc.hist.new, line by line.
! See the CASINO 1.xxx manual for the definition of each data item.
  i=0
  do
   i=i+1
   if(mod(i-1,vmc_block_length)==0)then
    read(8,*,iostat=ierr)
    if(ierr>0)then
     write(6,*)' Problem reading block header in vmc.hist (1).'
     stop
    elseif(ierr<0)then
     write(6,*)' End of vmc.hist occurs at line number '//trim(i2s(i-1))//'.'
     exit
    endif
    read(8,*,iostat=ierr)vmc_block_length
    if(ierr/=0)then
     write(6,*)' Problem reading block header in vmc.hist (2).'
     stop
    endif
    if(electron_gas)then
     read(8,*,iostat=ierr)
     if(ierr/=0)then
      write(6,*)' Problem reading vmc.hist.'
      stop
     endif
    endif ! electron gas
   endif ! Block header
   read(8,fmt="(4(1PE20.13))",iostat=ierr)row_data(1:no_cols_vmc_hist)
   if(ierr>0)then
     write(6,*)' Problem reading block header in vmc.hist (1).'
     stop
    elseif(ierr<0)then
    write(6,*)' End of vmc.hist occurs at line number '//trim(i2s(i-1))//'.'
    if(mod(i-1,vmc_block_length)/=0)write(6,*)' WARNING: this is not the &
     &end of a block.'
    exit
   endif

! Check that the energy components add up as expected.
   if(eionion_const)then
    etot_temp=row_data(3)+row_data(7)+row_data(8)+row_data(12) &
     &+row_data(13)+row_data(14)
    if(iterac==1.or.iterac==3)then
! Ewald used in total energy.
     etot_temp=etot_temp+row_data(6)
    else
! MPC used in total energy
     etot_temp=etot_temp+row_data(9)
    endif ! iterac
    if(relativistic)etot_temp=etot_temp+row_data(16)+row_data(17) &
     &+row_data(18)+row_data(19)+row_data(20)
    eionion_temp=row_data(1)-etot_temp
    if(i==1)eionion=eionion_temp
    if(abs(eionion_temp-eionion)>tol)then
     write(6,*)' WARNING: a component of energy is unaccounted for.'
     write(6,*)' First evaluation of ion-ion energy: ',eionion
     write(6,*)' Later evaluation of ion-ion energy: ',eionion_temp
     eionion_const=.false.
    endif ! eionion not constant
   endif ! eionion_const

   if(.not.mpc.and..not.vcpp.and..not.relativistic)then
    write(10,'(a,1x,6(es26.16,1x))')trim(i2s(i)),row_data(1),row_data(3), &
     &row_data(4),row_data(6),row_data(7),row_data(8)
   elseif(.not.mpc.and..not.vcpp.and.relativistic)then
    write(10,'(a,1x,11(es26.16,1x))')trim(i2s(i)),row_data(1),row_data(3), &
     &row_data(4),row_data(6),row_data(7),row_data(8),row_data(16), &
     &row_data(17),row_data(18),row_data(19),row_data(20)
   elseif(.not.mpc.and.vcpp.and..not.relativistic)then
    write(10,'(a,1x,9(es26.16,1x))')trim(i2s(i)),row_data(1),row_data(3), &
     &row_data(4),row_data(6),row_data(7),row_data(8),row_data(12), &
     &row_data(13),row_data(14)
   elseif(.not.mpc.and.vcpp.and.relativistic)then
    write(10,'(a,1x,14(es26.16,1x))')trim(i2s(i)),row_data(1),row_data(3), &
     &row_data(4),row_data(6),row_data(7),row_data(8),row_data(12), &
     &row_data(13),row_data(14),row_data(16),row_data(17),row_data(18), &
     &row_data(19),row_data(20)
   elseif(mpc.and..not.vcpp.and..not.relativistic)then
    write(10,'(a,1x,8(es26.16,1x))')trim(i2s(i)),row_data(1),row_data(3), &
     &row_data(4),row_data(6),row_data(7),row_data(8),row_data(10),row_data(11)
   elseif(mpc.and..not.vcpp.and.relativistic)then
    write(10,'(a,1x,13(es26.16,1x))')trim(i2s(i)),row_data(1),row_data(3), &
     &row_data(4),row_data(6),row_data(7),row_data(8),row_data(10), &
     &row_data(11),row_data(16),row_data(17),row_data(18),row_data(19), &
     &row_data(20)
   elseif(mpc.and.vcpp.and..not.relativistic)then
    write(10,'(a,1x,11(es26.16,1x))')trim(i2s(i)),row_data(1),row_data(3), &
     &row_data(4),row_data(6),row_data(7),row_data(8),row_data(10), &
     &row_data(11),row_data(12),row_data(13),row_data(14)
   else
    write(10,'(a,1x,16(es26.16,1x))')trim(i2s(i)),row_data(1),row_data(3), &
    &row_data(4),row_data(6),row_data(7),row_data(8),row_data(10), &
    &row_data(11),row_data(12),row_data(13),row_data(14),row_data(16), &
     &row_data(17),row_data(18),row_data(19),row_data(20)
   endif ! MPC/VCPP/Relativistic.

  enddo ! i

  close(8)
  close(10)
  deallocate(row_data)

  if(abs(eionion-eionion_exact)>tol.and.eionion_const)then
   write(6,*)' WARNING: ion-ion energy is different from the one quoted &
    &in the block header.'
   write(6,*)' This might be a problem.'
  endif ! eionion const, but wrong const.
  write(6,*)

 END SUBROUTINE convert_vmc


 SUBROUTINE convert_dmc
!---------------------------------------------------------------!
! This subroutine converts dmc.hist(2) (old) -> dmc.hist (new). !
!---------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER i,ialloc,ierr,no_cols_dmc_hist,no_cols_dmc_hist2,netot,btype, &
   &npcells,nbasis,iterac,periodicity,nequil,no_data_cols
  DOUBLE PRECISION eionion,eionion_temp,etot_temp
  DOUBLE PRECISION,ALLOCATABLE :: row_data1(:),row_data2(:),tempdata(:)
  CHARACTER(1) y_or_n
  CHARACTER(512) bigstring
  LOGICAL mpc,vcpp,relativistic,dmc_hist,dmc_hist2,ewald,awfn,bwfn,gwfn,pwfn, &
   &unsure
  LOGICAL :: eionion_const=.true.,warn_given=.false.

! Are both dmc.hist and dmc.hist2 present?
  inquire(file='dmc.hist',exist=dmc_hist)
  if(.not.dmc_hist)then
   write(6,*)' Please try to find the dmc.hist file.'
   stop
  endif
  inquire(file='dmc.hist2',exist=dmc_hist2)

  open(unit=8,file='dmc.hist',status='old',iostat=ierr)
  if(ierr/=0)then
   write(6,*)' Problem opening dmc.hist.'
   stop
  endif

! Find out how many items per line there are in dmc.hist. (As a check.)
  bigstring=''
  read(8,fmt='(a)',iostat=ierr)bigstring
  if(ierr/=0)then
   write(6,*)' Problem reading 1st line of dmc.hist file.'
   stop
  endif
  no_cols_dmc_hist=0
  do
   allocate(tempdata(no_cols_dmc_hist+1),stat=ialloc)
   if(ialloc/=0)then
    write(6,*)' Allocation error.'
    stop
   endif
   read(bigstring,*,iostat=ierr)tempdata(1:no_cols_dmc_hist+1)
   deallocate(tempdata)
   if(ierr/=0)exit
   no_cols_dmc_hist=no_cols_dmc_hist+1
  enddo ! no_cols_dmc_hist
  if(no_cols_dmc_hist==0)then
   write(6,*)' No data on first line of dmc.hist2? I/O error.'
   stop
  endif
  rewind(8)
  if(no_cols_dmc_hist/=13)then
   write(6,*)' Sorry, format of dmc.hist file is wrong.'
   write(6,*)' No. of columns in dmc.hist: '//trim(i2s(no_cols_dmc_hist))//'.'
   stop
  endif

  allocate(row_data1(1:no_cols_dmc_hist),stat=ierr)
  if(ierr/=0)then
   write(6,*)' Problem allocating row_data1.'
   stop
  endif

! Get number of atoms in primitive cell and periodicity.
  read(8,*)row_data1(1:no_cols_dmc_hist)
  nbasis=nint(row_data1(10))
  periodicity=nint(row_data1(11))
  if(nbasis<0.or.periodicity<0.or.periodicity>3)then
   write(6,*)' Problem with nbasis or periodicity in dmc.hist.'
   stop
  endif
  rewind(8)

  if(dmc_hist2)then

   open(unit=9,file='dmc.hist2',status='old',iostat=ierr)
   if(ierr/=0)then
    write(6,*)' Problem opening dmc.hist2.'
    stop
   endif

! Find out how many items per line are in dmc.hist2.
   bigstring=''
   read(9,fmt='(a)',iostat=ierr)bigstring
   if(ierr/=0)then
    write(6,*)' Problem reading 1st line of dmc.hist file.'
    stop
   endif
   no_cols_dmc_hist2=0
   do
    allocate(tempdata(no_cols_dmc_hist2+1),stat=ialloc)
    if(ialloc/=0)then
     write(6,*)' Allocation error.'
     stop
    endif
    read(bigstring,*,iostat=ierr)tempdata(1:no_cols_dmc_hist2+1)
    deallocate(tempdata)
    if(ierr/=0)exit
    no_cols_dmc_hist2=no_cols_dmc_hist2+1
   enddo ! no_cols_dmc_hist2
   if(no_cols_dmc_hist2==0)then
    write(6,*)' No data on first line of dmc.hist2? I/O error.'
    stop
   endif
   rewind(9)

   if(no_cols_dmc_hist2==17)then
    write(6,*)' Relativistic data are not present in dmc.hist2.'
    relativistic=.false.
   elseif(no_cols_dmc_hist2==23)then
    write(6,*)' Relativistic data are present in dmc.hist2.'
    relativistic=.true.
   else
    write(6,*)' Sorry, format of dmc.hist2 file is wrong.'
    write(6,*)' No. of columns in dmc.hist2: '//trim(i2s(no_cols_dmc_hist2)) &
     &//'.'
    stop
   endif
   write(6,*)

   allocate(row_data2(1:no_cols_dmc_hist2),stat=ierr)
   if(ierr/=0)then
    write(6,*)' Problem allocating row_data2.'
    stop
   endif

   read(9,*,iostat=ierr)row_data2(1:no_cols_dmc_hist2)
   if(row_data2(9)/=0.d0.or.row_data2(10)/=0.d0.or.row_data2(11)/=0.d0)then
    write(6,*)' Core-polarization pseudopotential data are present in &
     &dmc.hist2.'
    vcpp=.true.
   else
    write(6,*)' Core-polarization pseudopotential data are not present in &
     &dmc.hist2.'
    vcpp=.false.
   endif ! VCPP present.
   write(6,*)

   if(row_data2(5)/=0.d0.or.row_data2(6)/=0.d0)then
    write(6,*)' MPC data are present in dmc.hist2.'
    mpc=.true.
   else
    write(6,*)' MPC data are not present in dmc.hist2.'
    mpc=.false.
   endif ! MPC present.
   write(6,*)

   if(row_data2(4)/=0.d0)then
    write(6,*)' Ewald data are present in dmc.hist2.'
    ewald=.true.
   else
    write(6,*)' Ewald data are not present in dmc.hist2.'
    ewald=.false.
   endif ! MPC present.
   write(6,*)

   rewind(9)

  else

   do
    write(6,*)' Are you accumulating relativistic data (y/n)?'
    read(5,*,iostat=ierr)y_or_n
    if(ierr/=0)y_or_n='a'
    if(y_or_n=='y')y_or_n='Y'
    if(y_or_n=='n')y_or_n='N'
    if(y_or_n=='N')then
     relativistic=.false.
     exit
    elseif(y_or_n=='Y')then
     relativistic=.true.
     exit
    else
     write(6,*)' Please try again.  Enter "y" or "n".'
     write(6,*)
    endif
   enddo
   write(6,*)

   do
    write(6,*)' Are you accumulating core-polarization pseudopotential &
     &data (y/n)?'
    read(5,*,iostat=ierr)y_or_n
    if(ierr/=0)y_or_n='a'
    if(y_or_n=='y')y_or_n='Y'
    if(y_or_n=='n')y_or_n='N'
    if(y_or_n=='N')then
     vcpp=.false.
     exit
    elseif(y_or_n=='Y')then
     vcpp=.true.
     exit
    else
     write(6,*)' Please try again.  Enter "y" or "n".'
     write(6,*)
    endif
   enddo
   write(6,*)

   do
    write(6,*)' Are you accumulating MPC data (y/n)?'
    read(5,*,iostat=ierr)y_or_n
    if(ierr/=0)y_or_n='a'
    if(y_or_n=='y')y_or_n='Y'
    if(y_or_n=='n')y_or_n='N'
    if(y_or_n=='N')then
     mpc=.false.
     exit
    elseif(y_or_n=='Y')then
     mpc=.true.
     exit
    else
     write(6,*)' Please try again.  Enter "y" or "n".'
     write(6,*)
    endif
   enddo
   write(6,*)

   do
    write(6,*)' Are you accumulating Ewald data (y/n)?'
    read(5,*,iostat=ierr)y_or_n
    if(ierr/=0)y_or_n='a'
    if(y_or_n=='y')y_or_n='Y'
    if(y_or_n=='n')y_or_n='N'
    if(y_or_n=='N')then
     ewald=.false.
     exit
    elseif(y_or_n=='Y')then
     ewald=.true.
     exit
    else
     write(6,*)' Please try again.  Enter "y" or "n".'
     write(6,*)
    endif
   enddo
   write(6,*)

   do
    write(6,*)' Please enter the constant ion-ion energy.'
    read(5,*,iostat=ierr)eionion
    if(ierr/=0)then
     write(6,*)' Please try again.'
     write(6,*)
    else
     exit
    endif
   enddo
   write(6,*)

! Set the data in dmc.hist2 to zero if dmc.hist2 is not present.
   if(relativistic)then
    no_cols_dmc_hist2=22
   else
    no_cols_dmc_hist2=17
   endif ! relativistic
   allocate(row_data2(1:no_cols_dmc_hist2),stat=ierr)
   if(ierr/=0)then
    write(6,*)' Problem allocating row_data2.'
    stop
   endif
   row_data2(:)=0.d0

  endif ! dmc_hist2

  if(nbasis==0)then
   btype=0
  else
! Look for xwfn.data files.  If none exist or more than one exist then
! ask the user to supply btype.
   inquire(file='pwfn.data',exist=pwfn)
   inquire(file='bwfn.data',exist=bwfn)
   inquire(file='gwfn.data',exist=gwfn)
   inquire(file='awfn.data',exist=awfn)
   unsure=.not.(pwfn.or.bwfn.or.gwfn.or.awfn).or.(pwfn.and.bwfn) &
    &.or.(pwfn.and.gwfn).or.(pwfn.and.awfn).or.(bwfn.and.gwfn) &
    &.or.(bwfn.and.awfn).or.(gwfn.and.awfn)
   if(unsure)then
    do
     write(6,*)' Please enter the basis type (btype input parameter)'
     read(5,*,iostat=ierr)btype
     if(ierr/=0)btype=-1
     if(btype<1.or.btype>4)then
      write(6,*)' Should have 1<=btype<=4'
      write(6,*)
     else
      exit
     endif
    enddo
    write(6,*)
   else
    if(pwfn)btype=1
    if(bwfn)btype=4
    if(gwfn)btype=2
    if(awfn)btype=3
   endif ! unsure
  endif ! electron_gas
  write(6,*)' The btype input parameter is '//trim(i2s(btype))//'.'
  write(6,*)

! Work out iterac value.
  if(periodicity>0)then
   if(ewald.and..not.mpc)then
    iterac=1
   elseif(mpc.and..not.ewald)then
    iterac=2
   else
    do
     write(6,*)' Please enter the value of the iterac input parameter.'
     read(5,*,iostat=ierr)iterac
     if(ierr/=0)iterac=-1
     if(iterac<1.or.iterac>4)then
      write(6,*)' iterac must be between 1 and 4.'
      write(6,*)
     else
      exit
     endif
    enddo ! Get iterac value.
    write(6,*)
   endif
  else
   iterac=1
  endif ! periodicity
  write(6,*)' The iterac input parameter is '//trim(i2s(iterac))//'.'
  write(6,*)

  if(periodicity>0)then
   do
    write(6,*)' Please enter the number of primitive cells in your simulation.'
    read(5,*,iostat=ierr)npcells
    if(ierr/=0)npcells=-1
    if(npcells<1)then
     write(6,*)' Please try again.'
     write(6,*)
    else
     exit
    endif
   enddo
   write(6,*)
  else
   npcells=1
  endif
  if(npcells/=1)then
   write(6,*)' There are '//trim(i2s(npcells))//' primitive cells.'
  else
   write(6,*)' There is 1 primitive cell.'
  endif ! singular/plural
  write(6,*)

  do
   write(6,*)' Please enter the total number of particles in your simulation.'
   read(5,*,iostat=ierr)netot
   if(ierr/=0)netot=-1
   if(netot<1)then
    write(6,*)' Please try again.'
    write(6,*)
   else
    exit
   endif
  enddo
  write(6,*)
  if(netot/=1)then
   write(6,*)' There are '//trim(i2s(netot))//' particles in the simulation.'
  else
   write(6,*)' There is 1 particle in the simulation.'
  endif ! singular/plural
  write(6,*)

  do
   write(6,*)' Please enter the number of equilibration moves.'
   read(5,*,iostat=ierr)nequil
   if(ierr/=0)nequil=-1
   if(nequil<0)then
    write(6,*)' Please try again.'
    write(6,*)
   else
    exit
   endif
  enddo
  write(6,*)
  if(nequil/=1)then
   write(6,*)' The first '//trim(i2s(nequil)) &
    &//' moves are equilibration moves.'
  else
   write(6,*)' The first move is an equilibration move.'
  endif ! singular/plural
  write(6,*)

  open(unit=10,file='dmc.hist.new',status='replace',iostat=ierr)
  if(ierr/=0)then
   write(6,*)' Problem opening dmc.hist.new.  Stopping.'
   stop
  endif ! ierr

! Read data from dmc.hist(2) files and write to dmc.hist.new, line by line.
! See the CASINO 1.xxx manual for the definition of each data item.
  i=0
  do
   i=i+1
   read(8,*,iostat=ierr)row_data1(1:no_cols_dmc_hist)
   if(ierr<0)then
    write(6,*)' End of dmc.hist occurs at line number '//trim(i2s(i-1))//'.'
    exit
   elseif(ierr>0)then
    write(6,*)' Error in dmc.hist at line number '//trim(i2s(i-1))//'.'
    exit
   endif

   if(dmc_hist2)then
    read(9,*,iostat=ierr)row_data2(1:no_cols_dmc_hist2)
    if(ierr<0)then
     write(6,*)' End of dmc.hist2 occurs at line number '//trim(i2s(i-1))//'.'
     exit
    elseif(ierr>0)then
     write(6,*)' Error in dmc.hist2 at line number '//trim(i2s(i-1))//'.'
     exit
    endif
   endif ! dmc_hist2

   if(.not.warn_given.and.nint(row_data1(1))/=i)then
    write(6,*)' WARNING: the iteration numbers in dmc.hist behave oddly &
     &at line '//trim(i2s(i))//'.'
    write(6,*)' You might want to check this.'
    warn_given=.true.
   endif ! Line number unexpected.

   if(dmc_hist2.and.eionion_const)then
    etot_temp=row_data2(2)+row_data2(8)+row_data2(7)+row_data2(9) &
     &+row_data2(10)+row_data2(11)
    if(iterac==1.or.iterac==3)then
! Ewald used in total energy.
     etot_temp=etot_temp+row_data2(4)
    else
! MPC used in total energy.
     etot_temp=etot_temp+row_data2(5)+row_data2(6)
    endif
    if(relativistic)etot_temp=etot_temp+row_data2(18)+row_data2(19) &
     &+row_data2(20)+row_data2(21)+row_data2(22)
    eionion_temp=row_data1(3)-etot_temp
    if(i==1)then
     eionion=eionion_temp
     if(nbasis==0.or.(nbasis==1.and.periodicity==0))eionion=0.d0
    endif ! i=1
    if(abs(eionion_temp-eionion)>tol)then
     write(6,*)' WARNING: a component of energy is unaccounted for.'
     write(6,*)' First calculation of eionion: ',eionion
     write(6,*)' Later calculation of eionion: ',eionion_temp
     write(6,*)' NOTE THAT THE ION-ION ENERGY IN dmc.hist.new IS INCORRECT.'
     write(6,*)' Please update the header manually.'
     eionion_const=.false.
    endif ! eionion not constant
   endif ! Can check eionion.

! Produce header for top of dmc.hist.new.
   if(i==1)then
    write(10,'(a)')'# Title'
    write(10,'(a)')'# dmc.hist file converted from old format.'
    write(10,'(a)')'# File version number'
    write(10,'(a)')'# 1'
    write(10,'(a)')'# QMC method (VMC, DMC, PIMC, AFMC or RMC)'
    write(10,'(a)')'# DMC'
    write(10,'(a)')'# Electron-electron interaction type (iterac keyword)'
    write(10,'(a)')'# '//trim(i2s(iterac))
    write(10,'(a)')'# Constant (ion-ion) energy'
    write(10,'(a,es26.16)')'# ',eionion
    write(10,'(a)')'# Number of electrons (and other particles) in simulation'
    write(10,'(a)')'# '//trim(i2s(netot))
    write(10,'(a)')'# Number of atoms per primitive cell'
    write(10,'(a)')'# '//trim(i2s(nbasis))
    write(10,'(a)')'# Number of primitive cells'
    write(10,'(a)')'# '//trim(i2s(npcells))
    write(10,'(a)')'# Basis type (btype keyword)'
    write(10,'(a)')'# '//trim(i2s(btype))
    write(10,'(a)')'# Periodic (0=NO; 1=YES)'
    if(periodicity==0)then
     write(10,'(a)')'# 0'
    else
     write(10,'(a)')'# 1'
    endif
    write(10,'(a)')'# Number of data columns (excluding iteration number)'
    no_data_cols=12
    if(mpc)no_data_cols=no_data_cols+2
    if(vcpp)no_data_cols=no_data_cols+3
    if(relativistic)no_data_cols=no_data_cols+5
    write(10,'(a)')'# '//trim(i2s(no_data_cols))
    write(10,'(a)')'# Data items (in order)'
    write(10,'(a)')'# WEIGHT'
    write(10,'(a)')'# NCONF'
    write(10,'(a)')'# ETOT'
    write(10,'(a)')'# EREF'
    write(10,'(a)')'# EBEST'
    write(10,'(a)')'# ACC'
    write(10,'(a)')'# TEFF'
    write(10,'(a)')'# KEI'
    write(10,'(a)')'# TI'
    write(10,'(a)')'# EWALD'
    write(10,'(a)')'# LOCAL'
    write(10,'(a)')'# NONLOCAL'
    if(mpc)then
     write(10,'(a)')'# SHORT'
     write(10,'(a)')'# LONG'
    endif ! mpc
    if(vcpp)then
     write(10,'(a)')'# CPPEI'
     write(10,'(a)')'# CPPE'
     write(10,'(a)')'# CPPEE'
    endif ! vcpp
    if(relativistic)then
     write(10,'(a)')'# MASSPOL'
     write(10,'(a)')'# MASSVEL'
     write(10,'(a)')'# DARWINEN'
     write(10,'(a)')'# DARWINEE'
     write(10,'(a)')'# RETARD'
    endif ! relativistic
    write(10,'(a)')'# Raw QMC data'
    if(nequil==0)write(10,'(a)')'#### START STATS'
   endif ! i=1

   if(.not.mpc.and..not.vcpp.and..not.relativistic)then
    write(10,'(a,1x,es26.16,1x,a,1x,10(es26.16,1x))') &
     &trim(i2s(nint(row_data1(1)))),row_data1(13), &
     &trim(i2s(nint(row_data1(2)))),row_data1(3),row_data1(4), &
     &row_data1(5),row_data1(6),row_data2(17),row_data2(2), &
     &row_data2(1),row_data2(4),row_data2(8),row_data2(7)
   elseif(.not.mpc.and..not.vcpp.and.relativistic)then
    write(10,'(a,1x,es26.16,1x,a,1x,15(es26.16,1x))') &
     &trim(i2s(nint(row_data1(1)))),row_data1(13), &
     &trim(i2s(nint(row_data1(2)))),row_data1(3),row_data1(4), &
     &row_data1(5),row_data1(6),row_data2(17),row_data2(2), &
     &row_data2(1),row_data2(4),row_data2(8),row_data2(7), &
     &row_data2(18),row_data2(19),row_data2(20),row_data2(21),row_data2(22)
   elseif(.not.mpc.and.vcpp.and..not.relativistic)then
    write(10,'(a,1x,es26.16,1x,a,1x,13(es26.16,1x))') &
     &trim(i2s(nint(row_data1(1)))),row_data1(13), &
     &trim(i2s(nint(row_data1(2)))),row_data1(3),row_data1(4), &
     &row_data1(5),row_data1(6),row_data2(17),row_data2(2), &
     &row_data2(1),row_data2(4),row_data2(8),row_data2(7),row_data2(9), &
     &row_data2(10),row_data2(11)
   elseif(.not.mpc.and.vcpp.and.relativistic)then
    write(10,'(a,1x,es26.16,1x,a,1x,18(es26.16,1x))') &
     &trim(i2s(nint(row_data1(1)))),row_data1(13), &
     &trim(i2s(nint(row_data1(2)))),row_data1(3),row_data1(4), &
     &row_data1(5),row_data1(6),row_data2(17),row_data2(2), &
     &row_data2(1),row_data2(4),row_data2(8),row_data2(7),row_data2(9), &
     &row_data2(10),row_data2(11),row_data2(18),row_data2(19), &
     &row_data2(20),row_data2(21),row_data2(22)
   elseif(mpc.and..not.vcpp.and..not.relativistic)then
    write(10,'(a,1x,es26.16,1x,a,1x,12(es26.16,1x))') &
     &trim(i2s(nint(row_data1(1)))),row_data1(13), &
     &trim(i2s(nint(row_data1(2)))),row_data1(3),row_data1(4), &
     &row_data1(5),row_data1(6),row_data2(17),row_data2(2), &
     &row_data2(1),row_data2(4),row_data2(8),row_data2(7), &
     &row_data2(5),row_data2(6)
   elseif(mpc.and..not.vcpp.and.relativistic)then
    write(10,'(a,1x,es26.16,1x,a,1x,17(es26.16,1x))') &
     &trim(i2s(nint(row_data1(1)))),row_data1(13), &
     &trim(i2s(nint(row_data1(2)))),row_data1(3),row_data1(4), &
     &row_data1(5),row_data1(6),row_data2(17),row_data2(2), &
     &row_data2(1),row_data2(4),row_data2(8),row_data2(7), &
     &row_data2(5),row_data2(6),row_data2(18),row_data2(19), &
     &row_data2(20),row_data2(21),row_data2(22)
   elseif(mpc.and.vcpp.and..not.relativistic)then
    write(10,'(a,1x,es26.16,1x,a,1x,15(es26.16,1x))') &
     &trim(i2s(nint(row_data1(1)))),row_data1(13), &
     &trim(i2s(nint(row_data1(2)))),row_data1(3),row_data1(4), &
     &row_data1(5),row_data1(6),row_data2(17),row_data2(2), &
     &row_data2(1),row_data2(4),row_data2(8),row_data2(7), &
     &row_data2(5),row_data2(6),row_data2(9),row_data2(10),row_data2(11)
   else
    write(10,'(a,1x,es26.16,1x,a,1x,20(es26.16,1x))') &
     &trim(i2s(nint(row_data1(1)))),row_data1(13), &
     &trim(i2s(nint(row_data1(2)))),row_data1(3),row_data1(4), &
     &row_data1(5),row_data1(6),row_data2(17),row_data2(2), &
     &row_data2(1),row_data2(4),row_data2(8),row_data2(7),row_data2(5), &
     &row_data2(6),row_data2(9),row_data2(10),row_data2(11),row_data2(18), &
     &row_data2(19),row_data2(20),row_data2(21),row_data2(22)
   endif ! MPC/VCPP/Relativistic.

   if(i==nequil)write(10,'(a)')'#### START STATS'

  enddo ! i
  write(6,*)

  close(8)
  if(dmc_hist2)close(9)
  close(10)
  deallocate(row_data1,row_data2)

  if(dmc_hist2.and.eionion_const)then
   write(6,*)' The ion-ion interaction energy is: ',eionion
   write(6,*)' Please check that this is correct.  If it is not then please &
    &update the value'
   write(6,*)' in the header of dmc.hist.new.'
   write(6,*)
  endif

 END SUBROUTINE convert_dmc


 CHARACTER(20) FUNCTION i2s(n)
!-----------------------------------------------------------------------!
! I2S                                                                   !
! ===                                                                   !
! Convert integers to left-justified strings that can be printed in the !
! middle of a sentence without introducing large amounts of white space.!
!                                                                       !
! Calling routine is intended to include something like:                !
! USE utilities                                                         !
! INTEGER i                                                             !
! i=12                                                                  !
! write(6,*)'Integer number ',trim(i2s(i)),' with words at the end.'    !
!-----------------------------------------------------------------------!
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
   i2s(j:j)=achar(ichar('0')+mod(i,10))
   i=i/10
  enddo
  i=1 ; j=len_trim(i2s)
  do
   if(i>=j)exit
   tmp=i2s(j:j)
   i2s(j:j)=i2s(i:i)
   i2s(i:i)=tmp
   i=i+1
   j=j-1
  enddo
  i2s=trim(sign)//i2s
 END FUNCTION i2s


END MODULE utils


PROGRAM update_hist
!----------------------------!
! Main program starts here.  !
!----------------------------!
 USE utils
 IMPLICIT NONE
 LOGICAL vmc_hist,dmc_hist,dmc_hist2

 write(6,*)
 write(6,*)'UPDATE_HIST'
 write(6,*)'==========='
 write(6,*)

! Work out what file we are to read: vmc.hist, dmc.hist or dmc.hist2.
 inquire(file='vmc.hist',exist=vmc_hist)
 inquire(file='dmc.hist',exist=dmc_hist)
 inquire(file='dmc.hist2',exist=dmc_hist2)

 if(vmc_hist)then
  write(6,*)'Converting a vmc.hist file...'
  write(6,*)
  call convert_vmc
  write(6,*)'Done.  A new-format vmc.hist file has been produced.'
  write(6,*)'It is called vmc.hist.new.'
  write(6,*)
 endif ! vmc.hist file

 if(dmc_hist)then
  if(.not.dmc_hist2)write(6,*)'Converting a dmc.hist file...'
  if(dmc_hist2)write(6,*)'Converting dmc.hist and dmc.hist2 files...'
  write(6,*)
  call convert_dmc
  write(6,*)'Done.  A new-format dmc.hist file has been produced.'
  write(6,*)'It is called dmc.hist.new.'
  write(6,*)
 endif ! dmc.hist(2) file

 if(dmc_hist2.and..not.dmc_hist)then
  write(6,*)'An old-format dmc.hist2 file is present, but this cannot be &
   &converted without'
  write(6,*)'the accompanying data in dmc.hist.'
  write(6,*)
 endif ! dmc.hist2 but no dmc.hist.

 if(.not.vmc_hist.and..not.dmc_hist.and..not.dmc_hist2)then
  write(6,*)'There do not seem to be any .hist files in this directory.'
  write(6,*)'No .hist.new files have been generated.'
  write(6,*)
 endif ! no .hist files.

 write(6,*)'Program finished.'
 write(6,*)

END PROGRAM update_hist

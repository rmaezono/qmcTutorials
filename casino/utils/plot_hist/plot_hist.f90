!------------------------------------------------------------------------!
! PLOT_HIST                                                              !
! ========                                                               !
!                                                                        !
! Neil Drummond, 8.2005                                                  !
!                                                                        !
! This utility reads in a qmc.hist, vmc.hist or dmc.hist file and        !
! allows the user to plot his or her choice of data against move number. !
! The resulting (q/v/d)mc.hist.plot files can be read by                 !
! XMGrace.  Several items of data may be plotted.                        !
!------------------------------------------------------------------------!

MODULE utils
!---------------------------------!
! Miscellaneous subroutines, etc. !
!---------------------------------!
 IMPLICIT NONE
 PRIVATE
 PUBLIC plot_qmc
! Tags for the columns of the data file, specifying where each data item
! is held.  If a tag is negative, the data item isn't present.
 INTEGER tag_step,tag_energy,tag_K,tag_T,tag_fisq,tag_Ewald,tag_local,&
  &tag_nonlocal,tag_short,tag_long,tag_cppei,tag_cppe,tag_cppee,tag_masspol,&
  &tag_massvel,tag_darwinen,tag_darwinee,tag_retard,tag_dipole1,tag_dipole2, &
  &tag_dipole3,tag_dipolesq,tag_contact_den,tag_weight,tag_nconf,tag_eref, &
  &tag_ebest,tag_acc,tag_teff,tag_hf_ke,tag_hf_ex
! Filename of qmc.hist and qmc.hist.plot.
 CHARACTER(10) filename
 CHARACTER(15) filename_plot
! Information from the header of qmc.hist.
 INTEGER no_cols_qmc,iterac,isper_flag,version


CONTAINS


 SUBROUTINE plot_qmc
!---------------------------------------------------------------------------!
! This subroutine reads in data from (q/v/d)mc.hist, offers the user a menu !
! of columns to plot, and plots the selected data to (q/v/d)mc.hist.plot.   !
!---------------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER :: ierr,ialloc,choice,i,which_file
  DOUBLE PRECISION,ALLOCATABLE :: row_data(:)
  LOGICAL qmc,vmc,dmc
  LOGICAL :: first_time=.true.
  LOGICAL :: warned=.false.
  CHARACTER(1) temp
  CHARACTER(80) datastring
  CHARACTER(640) char640

! Establish whether qmc.hist, vmc.hist and dmc.hist should be plotted.
  inquire(file='qmc.hist',exist=qmc)
  inquire(file='dmc.hist',exist=dmc)
  inquire(file='vmc.hist',exist=vmc)
  if(qmc)then
   dmc=.false.
   vmc=.false.
  endif
  if(vmc.and.dmc)then
   do
    write(6,*)'Both vmc.hist and dmc.hist files exist in this directory.'
    write(6,*)'Would you like to plot data from (1) vmc.hist or (2) dmc.hist?'
    read(5,*,iostat=ierr)which_file
    if(ierr/=0)which_file=-1
    if(which_file==1)then
     dmc=.false.
     exit
    elseif(which_file==2)then
     vmc=.false.
     exit
    else
     write(6,*)'Please try again.'
     write(6,*)
    endif
   enddo
   write(6,*)
  endif ! Both VMC and DMC present.
  if(qmc)then
   write(6,*)'PLOT DATA FROM QMC.HIST'
   write(6,*)'======================'
   filename='qmc.hist'
  elseif(vmc)then
   write(6,*)'PLOT DATA FROM VMC.HIST'
   write(6,*)'======================'
   filename='vmc.hist'
  elseif(dmc)then
   write(6,*)'PLOT DATA FROM DMC.HIST'
   write(6,*)'======================'
   filename='dmc.hist'
  else
   write(6,*)'No qmc.hist, vmc.hist or dmc.hist files exist in this directory.'
   stop
  endif ! dmc etc
  write(6,*)

! Open the input and output files.
  open(unit=8,file=trim(filename),status='old',iostat=ierr)
  if(ierr/=0)then
   write(6,*)'Problem opening '//trim(filename)//'.'
   stop
  endif ! ierr/=0
  filename_plot=trim(adjustl(filename))//'.plot'
  open(unit=10,file=trim(filename_plot),status='replace',iostat=ierr)
  if(ierr/=0)then
   write(6,*)'Problem opening '//trim(filename_plot)//'.'
   stop
  endif ! ierr

! Read header from qmc.hist.
! Skip title.
  read(8,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  read(8,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
! Skip version number.
  read(8,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  read(8,*,iostat=ierr)temp,version
  call check_ierr(ierr)
  call check_hash(temp)
  if(version/=1)then
   write(6,*)'File version should be 1 in header.  Stopping.'
   stop
  endif ! version
! Skip QMC method.
  read(8,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  read(8,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
! Get interaction-type (iterac).
  read(8,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  read(8,*,iostat=ierr)temp,datastring
  call check_ierr(ierr)
  call check_hash(temp)
  datastring=adjustl(datastring)
  if(trim(datastring)=='none')then
   iterac=0
  elseif(trim(datastring)=='coulomb'.or.trim(datastring)=='ewald' &
   &.or.trim(datastring)=='manual')then
   iterac=1
  elseif(trim(datastring)=='mpc')then
   iterac=2
  elseif(trim(datastring)=='ewald_mpc')then
   iterac=3
  elseif(trim(datastring)=='mpc_ewald')then
   iterac=4
  else
! Check for old format.
   read(datastring,*,iostat=ierr)iterac
   call check_ierr(ierr)
   if(iterac<0.or.iterac>4)then
    write(6,*)'ITERAC value should be between 0 and 4.  Stopping.'
    stop
   endif ! iterac value out of range.
  endif ! interaction-type
! Skip constant (ion-ion) energy.
  read(8,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  read(8,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
! Skip total number of electrons.
  read(8,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  read(8,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
! Skip number of atoms per primitive cell.
  read(8,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  read(8,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
! Skip number of primitive cells.
  read(8,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  read(8,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
! Skip basis-type keyword.
  read(8,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  read(8,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
! Get periodicity.
  read(8,*,iostat=ierr)
  call check_ierr(ierr)
  call check_hash(temp)
  read(8,*,iostat=ierr)temp,isper_flag
  call check_ierr(ierr)
  call check_hash(temp)
  if(isper_flag/=0.and.isper_flag/=1)then
   write(6,*)'Periodicity flag must be 0 or 1.'
   stop
  endif ! periodicity.
  if(isper_flag==0.and.iterac/=0.and.iterac/=1)then
   write(6,*)'Interaction-type should be 0 or 1 for finite systems.  &
    &Contradiction in header.'
   stop
  endif
! Get number of data columns.  Increase it by 1, since the line-numbers will
! also be read.
  read(8,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  read(8,*,iostat=ierr)temp,no_cols_qmc
  call check_ierr(ierr)
  call check_hash(temp)
  if(no_cols_qmc<1)then
   write(6,*)'No data to plot.  Stopping.'
   stop
  endif
  no_cols_qmc=no_cols_qmc+1
! Get items in qmc.hist
  tag_step=1         ! Move number
  tag_energy=-1      ! Total energy
  tag_K=-1           ! KEI kinetic-energy estimator
  tag_T=-1           ! TI kinetic-energy estimator
  tag_fisq=-1        ! FISQ kinetic-energy estimator
  tag_Ewald=-1       ! 1/r or Ewald e-e interaction
  tag_local=-1       ! Local electron-ion energy
  tag_nonlocal=-1    ! Nonlocal electron-ion energy
  tag_short=-1       ! Short-range part of MPC
  tag_long=-1        ! Long-range part of MPC
  tag_cppei=-1       ! Electron-ion CPP term
  tag_cppe=-1        ! Electron CPP term
  tag_cppee=-1       ! Electron-electron CPP term
  tag_masspol=-1     ! Mass-polarization term
  tag_massvel=-1     ! Mass-velocity term
  tag_darwinen=-1    ! Darwin e-n term
  tag_darwinee=-1    ! Darwin e-e term
  tag_retard=-1      ! Retardation term
  tag_dipole1=-1     ! Dipole moment
  tag_dipole2=-1     !    "     "
  tag_dipole3=-1     !    "     "
  tag_dipolesq=-1    ! Square of dipole moment
  tag_contact_den=-1 ! Electron-positron contact density
  tag_weight=-1      ! Total weight of configs
  tag_nconf=-1       ! Number of configs
  tag_eref=-1        ! Reference energy
  tag_ebest=-1       ! Best estimate of energy
  tag_acc=-1         ! Acceptance ratio
  tag_teff=-1        ! Effective time step
  tag_hf_ke=-1       ! HEG HF KE
  tag_hf_ex=-1       ! HEG HF exchange energy
  read(8,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  do i=2,no_cols_qmc
   read(8,*,iostat=ierr)temp,datastring
   call check_ierr(ierr)
   call check_hash(temp)
   datastring=adjustl(datastring)
   if(trim(datastring)=='ETOT')then
    call check_tag_free(tag_energy)
    tag_energy=i
   elseif(trim(datastring)=='KEI')then
    call check_tag_free(tag_K)
    tag_K=i
   elseif(trim(datastring)=='TI')then
    call check_tag_free(tag_T)
    tag_T=i
   elseif(trim(datastring)=='FISQ')then
    call check_tag_free(tag_fisq)
    tag_fisq=i
   elseif(trim(datastring)=='EWALD')then
    call check_tag_free(tag_Ewald)
    tag_Ewald=i
   elseif(trim(datastring)=='LOCAL')then
    call check_tag_free(tag_local)
    tag_local=i
   elseif(trim(datastring)=='NONLOCAL')then
    call check_tag_free(tag_nonlocal)
    tag_nonlocal=i
   elseif(trim(datastring)=='SHORT')then
    call check_tag_free(tag_short)
    tag_short=i
   elseif(trim(datastring)=='LONG')then
    call check_tag_free(tag_long)
    tag_long=i
   elseif(trim(datastring)=='CPPEI')then
    call check_tag_free(tag_cppei)
    tag_cppei=i
   elseif(trim(datastring)=='CPPE')then
    call check_tag_free(tag_cppe)
    tag_cppe=i
   elseif(trim(datastring)=='CPPEE')then
    call check_tag_free(tag_cppee)
    tag_cppee=i
   elseif(trim(datastring)=='MASSPOL')then
    call check_tag_free(tag_masspol)
    tag_masspol=i
   elseif(trim(datastring)=='MASSVEL')then
    call check_tag_free(tag_massvel)
    tag_massvel=i
   elseif(trim(datastring)=='DARWINEN')then
    call check_tag_free(tag_darwinen)
    tag_darwinen=i
   elseif(trim(datastring)=='DARWINEE')then
    call check_tag_free(tag_darwinee)
    tag_darwinee=i
   elseif(trim(datastring)=='RETARD')then
    call check_tag_free(tag_retard)
    tag_retard=i
   elseif(trim(datastring)=='DIPOLE1')then
    call check_tag_free(tag_dipole1)
    tag_dipole1=i
   elseif(trim(datastring)=='DIPOLE2')then
    call check_tag_free(tag_dipole2)
    tag_dipole2=i
   elseif(trim(datastring)=='DIPOLE3')then
    call check_tag_free(tag_dipole3)
    tag_dipole3=i
   elseif(trim(datastring)=='DIPOLESQ')then
    call check_tag_free(tag_dipolesq)
    tag_dipolesq=i
   elseif(trim(datastring)=='CONTACT_DEN')then
    call check_tag_free(tag_contact_den)
    tag_contact_den=i
   elseif(trim(datastring)=='WEIGHT')then
    call check_tag_free(tag_weight)
    tag_weight=i
   elseif(trim(datastring)=='NCONF')then
    call check_tag_free(tag_nconf)
    tag_nconf=i
   elseif(trim(datastring)=='EREF')then
    call check_tag_free(tag_eref)
    tag_eref=i
   elseif(trim(datastring)=='EBEST')then
    call check_tag_free(tag_ebest)
    tag_ebest=i
   elseif(trim(datastring)=='ACC')then
    call check_tag_free(tag_acc)
    tag_acc=i
   elseif(trim(datastring)=='TEFF')then
    call check_tag_free(tag_teff)
    tag_teff=i
   elseif(trim(datastring)=='HF_KE')then
    call check_tag_free(tag_hf_ke)
    tag_hf_ke=i
   elseif(trim(datastring)=='HF_EX')then
    call check_tag_free(tag_hf_ex)
    tag_hf_ex=i
   else
    write(6,*)'Column label not recognised.'
    write(6,*)'Label is: '//trim(datastring)
    stop
   endif ! Label
  enddo ! i

  allocate(row_data(no_cols_qmc),stat=ialloc)
  if(ialloc/=0)then
   write(6,*)'Allocation error.'
   stop
  endif ! ialloc/=0

! Menu loop
  do

   write(6,*)'Please select one of the following items to plot:'
   do i=1,no_cols_qmc
    if(tag_step==i)write(6,'(" ",i2,a)')i, &
     &': Step counter;'
    if(tag_weight==i)write(6,'(" ",i2,a)')i, &
     &': Total weight of all configurations;'
    if(tag_nconf==i)write(6,'(" ",i2,a)')i, &
     &': Number of configurations;'
    if(iterac==2.or.iterac==4)then
     if(tag_energy==i)write(6,'(" ",i2,a)')i, &
      &': Total energy (using MPC);'
    else
     if(isper_flag==0.or.iterac==0)then
      if(tag_energy==i)write(6,'(" ",i2,a)')i, &
       &': Total energy;'
     else
      if(tag_energy==i)write(6,'(" ",i2,a)')i, &
       &': Total energy (using Ewald);'
     endif ! isperiodic
    endif ! iterac
    if(tag_eref==i)write(6,'(" ",i2,a)')i, &
     &': Reference energy;'
    if(tag_ebest==i)write(6,'(" ",i2,a)')i, &
     &': Best estimate of ground-state energy;'
    if(tag_acc==i)write(6,'(" ",i2,a)')i, &
     &': Move-acceptance ratio;'
    if(tag_teff==i)write(6,'(" ",i2,a)')i, &
     &': Effective time step;'
    if(tag_K==i)write(6,'(" ",i2,a)')i, &
     &': Kinetic energy (K);'
    if(tag_T==i)write(6,'(" ",i2,a)')i, &
     &': Kinetic energy (T);'
    if(tag_fisq==i)write(6,'(" ",i2,a)')i, &
     &': Kinetic energy (F);'
    if(isper_flag==0)then
     if(tag_Ewald==i)write(6,'(" ",i2,a)')i, &
      &': Interaction energy (Coulomb);'
    else
     if(tag_Ewald==i)write(6,'(" ",i2,a)')i, &
      &': Interaction energy (Ewald);'
    endif ! isperiodic
    if(tag_local==i)write(6,'(" ",i2,a)')i, &
     &': Local electron-ion energy (inc. external energy);'
    if(tag_nonlocal==i)write(6,'(" ",i2,a)')i, &
     &': Nonlocal electron-ion energy;'
    if(tag_short==i)write(6,'(" ",i2,a)')i, &
     &': Short-range part of MPC;'
    if(tag_long==i)write(6,'(" ",i2,a)')i, &
     &': Long-range part of MPC;'
    if(tag_cppei==i)write(6,'(" ",i2,a)')i, &
     &': Electron-ion CPP term;'
    if(tag_cppe==i)write(6,'(" ",i2,a)')i, &
     &': Electron CPP term;'
    if(tag_cppee==i)write(6,'(" ",i2,a)')i, &
     &': Electron-electron CPP term;'
    if(tag_masspol==i)write(6,'(" ",i2,a)')i, &
     &': Mass-polarization term;'
    if(tag_massvel==i)write(6,'(" ",i2,a)')i, &
     &': Mass-velocity term;'
    if(tag_darwinen==i)write(6,'(" ",i2,a)')i, &
     &': Electron-nucleus Darwin term;'
    if(tag_darwinee==i)write(6,'(" ",i2,a)')i, &
     &': Electron-electron Darwin term;'
    if(tag_retard==i)write(6,'(" ",i2,a)')i, &
     &': Retardation term;'
    if(tag_dipole1==i)write(6,'(" ",i2,a)')i, &
     &': Dipole moment (x-component);'
    if(tag_dipole2==i)write(6,'(" ",i2,a)')i, &
     &': Dipole moment (y-component);'
    if(tag_dipole3==i)write(6,'(" ",i2,a)')i, &
     &': Dipole moment (z-component);'
    if(tag_dipolesq==i)write(6,'(" ",i2,a)')i, &
     &': Square of dipole moment;'
    if(tag_contact_den==i)write(6,'(" ",i2,a)')i, &
     &': Electron-positron contact density;'
    if(tag_hf_ke==i)write(6,'(" ",i2,a)')i, &
     &': HEG Hartree-Fock kinetic energy;'
    if(tag_hf_ex==i)write(6,'(" ",i2,a)')i, &
     &': HEG Hartree-Fock exchange energy;'
   enddo ! i
   write(6,*)'or enter 0 to exit.'
   read(5,*,iostat=ierr)choice
   if(ierr/=0)choice=-1

   if(choice==0)then
! Exit.
    write(6,*)
    exit
   elseif(choice>=1.and.choice<=no_cols_qmc)then
! Make a plot.
    if(.not.first_time)write(10,*)'&'
    i=0
    do
     read(8,'(a)',iostat=ierr)char640
     if(ierr<0)then
      exit
     elseif(ierr>0)then
      write(6,*)'Problem reading dmc.hist.'
      stop
     endif
     if(index(char640,'#')==0)then
      i=i+1
      read(char640,*,iostat=ierr)row_data(1:no_cols_qmc)
      if(ierr/=0)then
       write(6,*)'Problem reading '//trim(filename)//'.'
       stop
      endif
      if(.not.warned.and.nint(row_data(1))/=i)then
       write(6,*)'Warning: line label wrong at line '//trim(i2s(i))//'.'
       warned=.true.
      endif
      write(10,*)i,row_data(choice)
     endif ! Line not a comment.
    enddo ! i
    rewind(8)
    first_time=.false.
    write(6,*)
    write(6,*)'The data have been plotted.'
   else
! User is clearly an idiot & cannot follow simple instructions.
    write(6,*)'Please try again.'
   endif ! choice
   write(6,*)

  enddo ! choice loop

  close(8)
  deallocate(row_data)

  write(6,*)'Finished plotting data.'
  if(.not.first_time)then
   write(6,*)'Data has been written to '//trim(filename_plot)//'.'
   close(10)
  else
   close(10,status='delete')
  endif ! Data written
  write(6,*)

 END SUBROUTINE plot_qmc


 SUBROUTINE check_hash(char)
!---------------------------------------------------------------------------!
! This sub is used to check that the 1st char in each header line is a "#". !
!---------------------------------------------------------------------------!
  IMPLICIT NONE
  CHARACTER(1),INTENT(in) :: char
  if(char/='#')then
   write(6,*)'Header line does not have a "#" in front.  Stopping.'
   stop
  endif
 END SUBROUTINE check_hash


 SUBROUTINE check_ierr(ierr)
!------------------------------------------------------!
! Complain if there has been a problem reading a file. !
!------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: ierr
  if(ierr/=0)then
   write(6,*)'Problem reading '//trim(filename)//'.'
   stop
  endif
 END SUBROUTINE check_ierr


 SUBROUTINE check_tag_free(tag)
!----------------------------------------------!
! Complain if a tag has already been assigned. !
!----------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: tag
  if(tag/=-1)then
   write(6,*)'Tag assigned twice.  Two column labels must be the same.'
   stop
  endif
 END SUBROUTINE check_tag_free


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
! write(6,*)'Integer number ',trim(i2s(i)),' with words at the end.'    !
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


END MODULE utils


PROGRAM plot_hist
!---------------------------!
! Main program starts here. !
!---------------------------!
 USE utils
 IMPLICIT NONE

 write(6,*)
 write(6,*)'O-----------O'
 write(6,*)'| PLOT_HIST |'
 write(6,*)'O-----------O'
 write(6,*)

 call plot_qmc

 write(6,*)'Program finished.'
 write(6,*)

END PROGRAM plot_hist


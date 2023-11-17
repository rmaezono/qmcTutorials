SUBROUTINE read_G9xout
!-------------------------------------------------------------------------!
! Identify whether the Gaussian run was simply a G.S. calculation or      !
! whether it produced excited states, either via CIS, TDDFT or CASSCF.    !
! If a CIS or TDDFT calc. was done then the user is given the choice of   !
! outputting the g.s. or one of the excited states.  The CASSCF state to  !
! output is taken from the 'state of interest' specified in the Gaussian  !
! route section.  If an excited state is required, the configurations and !
! coefficients are read from the output file.                             !
!-------------------------------------------------------------------------!
 USE awk_like
 USE cis_data
 USE g94_wavefunction
 IMPLICIT NONE

 INTEGER i,ifield,ifail
 INTEGER iconfigA ! Keeps count of configurations for each alpha excited state
 INTEGER iconfigB ! Ditto for beta excited states
 INTEGER iend ! Find where to cut-off the string containing n-n repulsion energy
 INTEGER ispin ! Loop over spins
 INTEGER iocca ! Which MO currently exciting from in CIS
 INTEGER ioccb
 INTEGER nfields
 INTEGER Nunpaired
 INTEGER iseof ! negative if EOF reached
 INTEGER ie,icon ! Loop counters
 INTEGER icas_label ! Label of CAS config. read from file
 INTEGER :: icheck=0 ! Used to test success of memory allocation/file handling
 INTEGER,PARAMETER :: io_in=33 ! Unit no. to read from
 INTEGER,PARAMETER :: io_out1=34 ! Units to output details of
 INTEGER,PARAMETER :: io_out2=35 ! excitations to
 INTEGER,ALLOCATABLE :: ilabels(:) ! Used to keep track of sorted CAS configs
 REAL(KIND=dp) weight ! Weight of excitation from a given MO in CIS
 REAL(KIND=dp) coeff_copy(Ncas_config_max) ! Used in reordering det coeffs
 LOGICAL done ! For control of reading HF eigenvalues
 LOGICAL search ! Flag to control scanning of output file for excited states.
 CHARACTER(22) new_string
 CHARACTER(35) input_file
 CHARACTER(35) filename ! For constructing names of out files
 CHARACTER(130) instring ! Storage for each line of out file, prior to parsing
 CHARACTER(130) instring_cp ! A copy of line with n-n repulsion energy
 CHARACTER(130) elements(NF_max) ! Same as fields in 'awk_like'

! RAndC: introduce a couple of variables
 CHARACTER(20) :: prevectext
 REAL(KIND=dp) :: tempc2dub=0.d0

 input_file=trim(g94_file)//'.out'

 open(unit=io_in,file=input_file,status='old',iostat=icheck)

 if(icheck>0)then
  write(*,fmt="(/'Cannot open file ',a35,'.')") input_file
  stop
 endif

! Identify what sort of Gaussian calculation was done. Take special
! action to read-in excited states if it was a CIS/TDDFT or CASSCF job.
 search=.true.
 do while(search)

  read(io_in,fmt="(a75)")instring
  icheck=index(string=trim(instring),substring="#")

  if(icheck/=0)then

! A # may also appear in the output file if it is combined with the log
! file so we have to check for that too...
   icheck=index(string=trim(instring),substring="LSBATCH:")
   if(icheck==0)then

    search=.false.

    call capitalise(instring,len(instring))
    icheck=index(string=trim(instring),substring="CIS")
    if(icheck==0)then

     icheck=index(string=trim(instring),substring=" TD")
     if(icheck==0)then

! This is not a CIS or TDDFT calculation output file so
! we must check to see whether it is a CAS output
      CIS=.false.
      icheck=index(string=trim(instring),substring="CAS")
      if(icheck==0)then ! This is not a CAS calculation output file
       CAS=.false.
      else ! Identify the 'state of interest' in the CAS calculation
       call scan_string(instring,fields,NF,NF_max," ")
       do ifield=1,NF
        icheck=index(string=fields(ifield),substring="NROOT")
        if(icheck/=0)exit
       enddo

       if(icheck==0)then ! Nroot not specified so takes default value of 1
        Nroot=1
       else ! Take 'n' where field(ifield) contains NROOT=n
        read(fields(ifield)(icheck+6:icheck+6),fmt="(i1)")Nroot
       endif
      endif

     else ! We are dealing with a TDDFT output file
      CAS=.false.
     endif

    else ! We are dealing with a CIS output file
     CAS=.false.
    endif
   endif
  endif
 enddo

! Get the spin multiplicity of the CASSCF wave function
 if(CAS)then
  call getline(io_in,"Multiplicity",fields,NF,NF_max,ifail,.true.)
  if(ifail<0)then
   write(*,fmt="('Error searching for ''Multiplicity''.')")
   stop
  endif
  read(fields(6),*) Multiplicity
  write(*,fmt="(/'Spin multiplicity = ',i1)")Multiplicity
 endif

! Get the nuclear-nuclear repulsion energy from the output file
 write(*,fmt="(/'Reading n-n repulsion energy from ',a35)")input_file
 search=.true.
 iseof=0
 instring_cp=" "
 if(Natom>1)then
! Search for last occurence of n-n repulsion energy
  do while(iseof==0)

   read(io_in,fmt="(a75)", iostat=iseof)instring
   if(iseof==0)then
    if(instring(8:31)=="nuclear repulsion energy")then
     instring_cp=trim(instring)
    endif 
   endif
  enddo
  if(len_trim(instring_cp)>0)then
   iend=index(string=instring_cp,substring="Hartrees")

   new_string=instring_cp(32:(iend-1))
! adjustl appears to fail if it is asked to act on a substring
! hence the need for new_string
   new_string=adjustl(new_string)
   read(new_string,*) Eionion
! Divide by the no. of atoms in the system
   Eionion=Eionion/real(Natom,dp)
  else
   Eionion=0.d0
  endif
 else
  Eionion=0.d0
 endif

! Set up the parameters describing the ground state in case there are no
! excited states
 Nconfig=0
 Nexcite=0
! The G.S. consists of a product of one alpha determinant and 1 beta det.
! for a spin polarised calc (provided that we have at least 1 alpha and
! at least 1 beta electron, respectively).
 if(Nspin(1)>0)Nconfig(1)=1
 if(SPIN.and.(Nspin(2)>0))Nconfig(2)=1

 ci_coeff(1,1)=1.d0 ! No expansion coeffs for G.S.

 if(CIS)then
! Next section gets the excitations and their associated prefactors in
! the CIS expansion.
  write(*,fmt="(/'Scanning CIS/TDDFT excitations from ',a35)")input_file
  search=.true.

  do while(search)

   read(io_in,fmt="(a75)",iostat=icheck)instring

   if(icheck==0)then

    if(instring(1:14)==" Excited State")then
     Nexcite=Nexcite+1
     write(6,*)trim(instring)
    endif

   elseif(icheck<0)then
    search=.false. ! We have reached the end of the file
    if(Nexcite==0)then
     write(*,fmt="('ERROR - CIS specified in route section but failed to&
      & find any excited states in out file.'//'Outputting the ground state.')")
     CIS=.false.
    endif

   else
    write(*,fmt="('Error reading file ',a35)")input_file
    stop
   endif
  enddo

  if(Nexcite>0)then
   write(*,fmt="(/'Construct CASINO file for an excited state (y/n)?')")
   read(*,fmt="(a1)")instring

   if((instring=='y').or.(instring=='Y'))then

    if(Nexcite>1)then
     write(*,fmt="(/'There are ',i2,' excited states, which one do you &
      &want'/'to output?')")Nexcite
     read(*,fmt="(i2)")icis_out
    else
     write(*,fmt="(/'Only one excited state so outputting that one.')")
     icis_out=1
    endif
   endif
  else ! Output the ground state
   icis_out=0
  endif

  if(icis_out>0)then
! Now that we know which excited state we're going to need we can go back
! and get it...
   rewind(io_in)

   write(*,fmt="(/'Reading CIS excitation No. ',i2)")icis_out
   search=.true.
   Nexcite=0
   iocca=0 ! Current alpha occ. MO being excited from
   ioccb=0 !   "     beta   "    "  "      "      "

   do while(search)

    read(io_in,fmt="(a75)")instring
    if(instring(1:14)==" Excited State")then

     Nexcite=Nexcite+1
     if(Nexcite==icis_out)then

! If this is a spin-restricted calculation then we must identify
! the spin state of the excitation so that we can reconstruct it
! correctly for the QMC calculation
      if(.not.SPIN)then
       singlet=.true.
       if(instring(23:29)=="Triplet")singlet=.false.
      endif

      iconfiga=0 ! Zero the counters for the number of alpha/beta
      iconfigb=0 ! configurations contributing to this excited state

      read(io_in,fmt="(a75)")instring

      do while(instring(10:11)=="->")

! Check to see whether this is an alpha or beta orbital and
! store this information.  If it is not a spin polarised
! calc. then treat the orbitals as alpha orbitals
       if((.not.SPIN).or.(instring(8:8)=='A'))then
        iconfiga=iconfiga+1 ! Increment alpha config counter for this state
        if(iconfiga>Ncon_max)then
         write(*,fmt="('Only ',i4,' configurations per excited state allowed - &
          &increase Ncon_max in cis_data.f90.')")Ncon_max
         stop
        endif

        if(.not.SPIN)then ! Formatting slightly different in this case...
         read(instring(6:8),fmt="(i3)")Orbitals(1,iconfiga,1)
        else
         read(instring(5:7),fmt="(i3)")Orbitals(1,iconfiga,1)
        endif
! The MO excited into
        read(instring(12:14),fmt="(i3)")Orbitals(2,iconfiga,1)

! Finally, the CIS expansion coeff. for this configuration
        read(instring(23:),fmt=*)ci_coeff(iconfiga,1)

! Write all excitations from alpha orbital I to fromI.1.dat as MO excited into,
! percentage weight..
        if(Orbitals(1,iconfiga,1)/=iocca)then
! Update which occupied MO we are currently exciting from
         iocca=Orbitals(1,iconfiga,1)
! Open new file to store details of excitations from this MO
         if(iocca/=0)close(io_out1)
         write(filename,*)iocca
         filename="from"//trim(adjustl(filename))//".1.dat"
         open(unit=io_out1,file=trim(filename),status='unknown')
        endif

! We may not have 100% of the wave function (it depends on what tolerance is
! used to decide whether or not to output a configuration) so don't convert
! to a percentage yet...
!       weight=100.d0*ci_coeff(iconfiga,1)*ci_coeff(iconfiga,1)
        weight=ci_coeff(iconfiga,1)*ci_coeff(iconfiga,1)
! Gaussian only outputs spin-up excitation if it is a spin-restricted calc. so
! need to multiply by two to account for the spin-down excitation too.
        if(.not.SPIN)weight=2.d0*weight
        write(io_out1,*)Orbitals(2,iconfiga,1),weight

       else ! This is a beta orbital...
        iconfigb=iconfigb+1 ! Increment beta config'n counter for this state
        if(iconfigb>Ncon_max)then
         write (*,fmt="('Only ',i4,' configurations per excited&
          & state allowed - increase Ncon_max in cis_data.f90.')")Ncon_max
         stop
        endif
! The MO excited from
        read(instring(5:7),fmt="(i3)")Orbitals(1,iconfigb,2)
! The MO excited into
        read(instring(12:14),fmt="(i3)")Orbitals(2,iconfigb,2)
! Associated expansion coefficient
        read(instring(23:32),fmt="(f10.7)")ci_coeff(iconfigb,2)

! Write all excitations from beta orbital I to fromI.2.dat as MO excited into,
! percentage weight..
        if(Orbitals(1,iconfigb,2)/=ioccb)then
! Update which occupied MO we are currently exciting from
         ioccb=Orbitals(1,iconfigb,2)
! Open new file to store details of excitations from this MO
         if(ioccb/=0)close(io_out2)
         write(filename,*)ioccb
         filename="from"//trim(adjustl(filename))//".2.dat"
         open(unit=io_out2,file=trim(filename),status='unknown')
        endif

        weight=100.d0*ci_coeff(iconfigb,2)*ci_coeff(iconfigb,2)
! Gaussian only outputs spin-up excitation if it is a spin-
! restricted calc. so need to multiply by two to account for
! the spin-down excitation too.
        if(.not.SPIN)weight=2.d0*weight
        write(io_out2,*)Orbitals(2,iconfigb,2),weight

       endif

       read(io_in,fmt="(a75)")instring ! Get the next line of the file
      enddo

      search=.false. ! We've found the state we were after so that's it
     endif
    endif
   enddo ! End of 'search' loop looking for excited states

   Nconfig(1)=iconfiga ! Store the no. of configs in this state
   if(SPIN)then
    Nconfig(2)=iconfigb
    write(*,fmt="('...it consists of ',i4,' alpha and ',i4,' beta &
     &configurations.')")iconfiga,iconfigb
   else
    write(*,fmt="('...it consists of ',i4,' configurations.')")iconfiga
   endif

  endif

 elseif(CAS)then
  rewind(io_in)
! Read-in the CASSCF expansion of the state CONVERGED TO in the Gaussian
! calculation - this is specified by the NROOT keyword in the root section
  write(*,fmt="(/'Reading CASSCF state number ',i1,' from ',a35/)") &
   &Nroot,input_file

! Get number of active electrons and orbitals. The order in which these are
! given depends on what IOp's (if any) are used...
  call getline(io_in,'active',fields,NF,NF_max,ifail,.true.)
  do i=1,NF
   icheck=index(string=trim(fields(i)),substring="orbitals")
   if(icheck/=0)then ! This is 'standard' output with orbitals given first
    read(fields(NF),*)nact_orbs
    call getline(io_in,' ',fields,NF,NF_max,ifail,.true.)
    read(fields(NF),*)nact_elecs
    exit
   endif
  enddo

  if(icheck==0)then ! This is IOp-forced output with electrons first
   read(fields(NF),*)nact_elecs
   call getline(io_in,' ',fields,NF,NF_max,ifail,.true.)
   read(fields(NF),*)nact_orbs
  endif

! Determine the number of alpha and beta electrons in the active
! space.  In addition, determine the number of frozen, occupied MO so
! that we can convert from Gaussian's numbering of the active
! orbitals back to the 'ROHF' numbering we will need to specify the
! excitations for CASINO.
  Nunpaired=abs(Nspin(1)-Nspin(2))
  if(Nunpaired<=nact_elecs)then
   Nspin_cas(1)=(nact_elecs-Nunpaired)/2+Nunpaired
   Nspin_cas(2)=nact_elecs-Nspin_cas(1)
   Nfrozen=(Nelec-nact_elecs)/2
  else
   Nspin_cas(1)=nact_elecs
   Nspin_cas(2)=0
   Nfrozen=Nunpaired-Nact_elecs ! No. of unpaired electrons in frozen MO
! No. of occupied, frozen MOs = No. of singly occupied + no. of
! doubly occupied frozen MOs...
   Nfrozen=Nfrozen+(Nelec-nact_elecs-Nfrozen)/2
  endif

  write(*,fmt="('No. of active spin-up = ',i2,/'No. of active spin-down &
   &= ',i2)")Nspin_cas
  write(*,fmt="('No. of occupied, frozen MOs = ',i2)")Nfrozen

! Check whether or not the CASSCF calculation used a Slater Determinant
! basis - read the line following that with the number of electrons on.
  call getline(io_in,'SLATER',fields,NF,NF_max,ifail,.false.)
  if(ifail<0)then
   slater=.false.
   rewind(io_in)
  else
   slater=.true.
  endif

! Currently G98 only outputs a maximum of 50 configurations and coefficients
! in each CAS eigenvector and so that's the value of Ncas_config_max which is
! set in cis_data.f90

! How many basis functions there are in the CASSCF calc.
! If it is less than Ncas_config_max then it will be the number of
! elements of the eigenvector that we will have to deal with.
  call getline(io_in,'DELETED',fields,NF,NF_max,ifail,.true.)
  if(ifail>0)then
   read(fields(6),*)Ncas_basis
  else
   rewind(io_in)
! No easy way to find how many basis functions we have, therefore
! ASSUME that this is a big calculation and therefore Gaussian only
! outputs Ncas_config_max (=50 for G98 rev. A9).
   Ncas_basis=Ncas_config_max
  endif

  allocate(icasdet(Nact_elecs,Ncas_basis),ilabels(Ncas_config_max),&
   &icas_config(Ncas_config_max),cas_coeffs(Ncas_config_max),&
   &icas_ref(Nact_elecs),stat=icheck)
  if(icheck/=0)stop 'Could not allocate memory in read_G9xout.'

! First we read the states and associated coefficients for the CASSCF
! expansion of the state of interest (specified by Nroot).

! RAndC b: later versions of G03 (checked upto E.01) put the search string
!          "final printing" in small case, just to be awkward.
  if(trim(adjustl(code_used))=="Gaussian 03, Revision C or later" &
     & .OR. trim(adjustl(code_used))=="Gaussian 09")then
   prevectext="final printing."
  else
   prevectext="FINAL PRINTING"
  endif ! RAndC e

! RAndC b: comment oot this chunk
!  call getline(io_in,'FINAL PRINTING',fields,NF,NF_max,ifail,.true.)
!  if(ifail<0)then
!   write(*,fmt="('Error searching for ''FINAL PRINTING''.')")
!   stop
!  endif
!  call getline(io_in,' ',fields,NF,NF_max,ifail,.true.)
!  if(ifail<0)then
!   write(*,fmt="('Unexpected end of file after ''FINAL PRINTING''.')")
!   stop
!  endif
! RAndC e

! RAndC b: replace the chunk immediately above
  call getline(io_in,prevectext,fields,NF,NF_max,ifail,.true.)
  if(ifail<0)then
   write(*,fmt=*) " Error searching for ", prevectext,"."
   stop
  endif
!   call getline(io_in,' ',fields,NF,NF_max,ifail,.true.)
   call getline(io_in,'EIGENVALUES ',fields,NF,NF_max,ifail,.true.)
  if(ifail<0)then
   write(*,fmt="('Unexpected end of file after ''FINAL PRINTING''.')")
   write(*,fmt="('> Was P specified in the Gaussian input route section, &
    &i.e., #P?')")
   stop
  endif
! RAndC e

! RAndC b: cope with the case-change.
! *NB* Not using this part of the IF statement in G03
  if(trim(fields(1))=="FINAL".or.trim(fields(1))=="final")then
!  if (trim(fields(1))=="FINAL")then
! RAndC e

! This is a calculation where the coefficients of ALL of the configurations
! are printed out in order

! Look for the line containing the start of the Nroot'th vector..
   write(instring,*)Nroot
   call getline(io_in,' '//trim(adjustl(instring))//' ',fields,NF, &
    &NF_max,ifail,.true.)
   if(ifail<0)then
    write(*,fmt="('Error reading CAS vector.')")
    stop
   endif

   icon=0 ! Zero counter of no. of config.'s  in CAS expansion

! First two numbers on the line are the root no. and eigenvalue
! respectively so skip these
   do i=3,NF
    icon=icon+1
! Configurations are output in ascending order so indexing
! is easy but still has to be done...
    icas_config(icon)=icon
! Read associated coefficient
    read(fields(i),*)cas_coeffs(icon)
   enddo

! Continue reading the eigenvector until we have all coefficients

   do while(icon<Ncas_basis)
! Read next line - this will contain only coefficients
    call getline(io_in,' ',fields,NF,NF_max,ifail,.true.)
    if(ifail<0)call fatal('Hit EOF reading CAS coefficients.')

    do i=1,NF
     icon=icon+1
     icas_config(icon)=icon
! Read associated coefficient
     read(fields(i),*)cas_coeffs(icon)
    enddo
   enddo

   Ncas_basis=icon
! RAndC b: see modification above around line 457
  elseif(trim(fields(1))=="EIGENVALUES")then
!  elseif(trim(fields(NF-1)) == "EIGENVALUE")then
! RAndC e

! 'Standard' G98 output where only the first 50 (say) coefficients
! are output and therefore have to be labelled

! Move forward to the root of interest
!   do i=1,Nroot-1
    call getline(io_in,'EIGENVALUE',fields,NF,NF_max,ifail,.true.)
    if(ifail<0)call fatal('Error searching for ''EIGENVALUE''.')
!   enddo

   icon=0 ! Zero counter of number of configurations  in CAS expansion

   search=.true.
   do while(search)
    read(io_in,fmt="(a130)")instring
    call scan_string(trim(adjustl(instring)),fields,NF,NF_max," ")
    if(fields(1)(1:1)=="(")then
     call scan_string(trim(adjustl(instring)),fields,NF,NF_max,"(")

     nfields=NF
     do i=1,nfields
      icon=icon+1
      call scan_string(trim(adjustl(fields(i))),elements,NF,NF_max,")")
! For big CAS calculations Gaussian's output is broken and get
! '*****' instead of configuration number. Have to skip these cases..
! RAndC b: read the value of CSF coefficient for later use
! convert from read string to double.
      read(elements(2),*)tempc2dub
! RAndC e
      if(elements(1)(1:1)=='*')then
       icon=icon-1
       cycle
! RAndC b:  Ignore coefficients close to zero.
      elseif(DABS(tempc2dub)<=1.d-7)then
       write(*,*)
       write(*,*)' > WARNING: CSFs with coefficients < 1.0D-07 have been &
        &dropped.'
       icon=icon-1
       cycle
! RAndC e
      endif
      read(elements(1),*)icas_config(icon)
      read(elements(2),*)cas_coeffs(icon)
     enddo
    else
     search=.false.
    endif
   enddo

! Correct the number of basis functions such that it is now just the
! number that are used explicitly in the wave function expansion
   Ncas_basis=icon
! RAndC b: error message says it all
   if(Ncas_basis==50.and.abs(cas_coeffs(50))>1.1d-7)then
    write(*,*)
    write(*,*)' > Gaussian doesn''t currently print out more than 50 CSF &
     &coefficients.'
    write(*,*)' > Perhaps this feature is hidden in an obscure IOp?'
    write(*,*)' > Anyway, the sum of the squares of the coefficients will &
     &not equal 1.0.'
   endif
! RAndC e

  else
   write(6,fmt="('Unrecognised format for CAS coefficent output.')")
   stop
  endif

! Now we read only those configurations actually involved in the CAS
! eigenvector expansion
! We must rewind the G98/4 output file and read the configurations that
! correspond to the labels read from the CAS expansion vector
  rewind(io_in)

! ilabels will keep track of the reordering that occurs when the list
! of config.'s in the expansion is sorted into ascending order
  do i=1,Ncas_config_max
   ilabels(i)=i
  enddo
!!$ do i=1,Ncas_basis
!!$  write(66,*)icas_config(i),ilabels(i)
!!$ enddo
  call numsrt_2way(icas_config,Ncas_basis,Ncas_config_max,ilabels,'a')
!!$ write(66,*) 'After sort..'
!!$ do i=1,Ncas_basis
!!$  write(66,*)icas_config(i),ilabels(i)
!!$ enddo
! Now reorder coefficients too so as to avoid confusion later
  coeff_copy(1:Ncas_basis)=cas_coeffs(ilabels(1:Ncas_basis))
  cas_coeffs(1:Ncas_basis)=coeff_copy(1:Ncas_basis)

  call getline(io_in,'PRIMARY',fields,NF,NF_max,ifail,.true.)
  if(ifail>0)then
! Output is 'standard' with reference config. labelled separately
! First we read the ground state/reference configuration...
   if(trim(adjustl(code_used))=="Gaussian 98")then
    do ie=1,Nact_elecs
! RAndC b: use of the "2" isn't necessary in release G98 A.11.
!    read(fields(2+ie),fmt="(i1)")icas_ref(ie)
     read(fields(3+ie),fmt="(i1)")icas_ref(ie)
! RAndC e
    enddo
   else
    do ie=1,Nact_elecs
     read(fields(3+ie),fmt="(i1)")icas_ref(ie)
    enddo
   endif

! Now we look for the relevant configurations...
   search=.true.
   if(icas_config(1)==1)then
! CAS expansion includes reference config. which needs to be skipped
    icasdet(:,1)=icas_ref
    icon=2
   else
    icon=1
   endif

   do while(search)
! Read the label for this configuration
    call getline(io_in,' ',fields,NF,NF_max,ifail,.true.)
    if(ifail<0)call fatal('Hit end of file in reading config. def.''s')
    read(fields(1),*)icas_label

! Check to see whether it is on our list
    if(icas_label==icas_config(icon))then
! Alpha MOs
     call getline(io_in,' ',fields,NF,NF_max,ifail,.true.)
     if(ifail<0)call fatal('Hit end of file in reading config. def.''s')
     do ie=1,Nspin_cas(1)
      read(fields(ie),*)icasdet(ie,icon)
     enddo
! Beta MOs
     call getline(io_in,' ',fields,NF,NF_max,ifail,.true.)
     if(ifail < 0)call fatal('Hit end of file in reading config. def.''s')
     do ie=1,Nspin_cas(2)
      read(fields(ie),*)icasdet(ie+Nspin_cas(1),icon)
     enddo

     icon=icon+1
! Have we now found all of the necessary configurations?
     if(icon>Ncas_basis)search=.false.
    else ! It isn't so skip it..
     read(io_in,*)
     read(io_in,*)
    endif
   enddo

  else
! Output is forced with IOps and reference config is not labelled
! separately - it is simply no. 1 in the list
   REWIND(io_in)
   call getline(io_in,'Slater Determinants',fields,NF,NF_max,ifail,.true.)
   if(ifail<0)call fatal('Could not identify start of configs list.')
   read(io_in,*)
   read(io_in,*)
   call getline(io_in,' ',fields,NF,NF_max,ifail,.true.)
   if(ifail<0)call fatal('Hit EOF whilst looking for reference config.')
   do ie=1,Nspin_cas(1)
    read(fields(ie),*)icas_ref(ie)
   enddo
   call getline(io_in,' ',fields,NF,NF_max,ifail,.true.)
   if(ifail<0)call fatal('Hit EOF whilst looking for reference config.')
   do i=1,Nspin_cas(2)
    ie=i+Nspin_cas(1)
    read(fields(i),*)icas_ref(ie)
   enddo

! Does CAS expansion include reference config.?
   if(icas_config(1)==1)icasdet(:,1)=icas_ref

! Now we look for the (rest of the) relevant configurations...
   search=.true.
   icon=2
   do while(search)
! Have a blank line between configurations in this format...
    read(io_in,*)
! Read the label for this configuration
    call getline(io_in,' ',fields,NF,NF_max,ifail,.true.)
    if(ifail<0)call fatal('Hit end of file in reading config. definitions.')
    read(fields(1),*)icas_label

! Check to see whether it is on our list
    if(icas_label==icas_config(icon))then
! Alpha MOs
     call getline(io_in,' ',fields,NF,NF_max,ifail,.true.)
     if(ifail<0)call fatal('Hit end of file in reading config. definitions.')
     do ie=1,Nspin_cas(1)
      read(fields(ie),*)icasdet(ie,icon)
     enddo
! Beta MOs
     call getline(io_in,' ',fields,NF,NF_max,ifail,.true.)
     if(ifail < 0)call fatal('Hit end of file in reading config. definitions.')
     do ie=1,Nspin_cas(2),1
      read(fields(ie),*) icasdet(ie+Nspin_cas(1),icon)
     enddo

     icon=icon+1
! Have we now found all of the necessary configurations?
     if(icon>Ncas_basis)search=.false.
    else ! It isn't so skip it...
     read(io_in,*)
     read(io_in,*)
    endif
   enddo

  endif

! At this point we can convert from Gaussian's numbering for the
! active orbitals back to the original numbering scheme.  To do this
! we use the number of frozen, occupied orbitals (Nfrozen) determined
! earlier.

  icas_ref(:)=icas_ref(:)+Nfrozen
  icasdet(:,1:Ncas_basis)=icasdet(:,1:Ncas_basis)+Nfrozen

 endif ! End of CAS input section

! This is an entirely separate part for analysing the percentage breakdown
! of the CIS excited state
 if(CIS)then
  write(*,fmt="(/'Reading HF spectrum from output file...'/)")
  rewind(io_in)
! How many distinct occupied orbitals do we have?
  if(SPIN)then ! spin-polarised calc.
   Nocc_mo=Nspin(1)+Nspin(2)
  else ! spin-restricted calc. - cannot do an "ROCIS" calculation so
       ! Nspin(1) must = Nspin(2) but we'd better check...
   if(Nspin(1)/=Nspin(2))then
    write(*,fmt="('WARNING: is a spin-restricted CIS calculation but &
     &Nspin(1) is not equal to Nspin(2) - skipping analysis of CIS state')")
    analyze_cis=.false.
   endif
   Nocc_mo=Nspin(1)
  endif

  allocate(hf_spectrum(Nmo,2),stat=icheck)
  if(icheck<0)then
   write(6,fmt="('Failed to allocate memory for eigenvalue &
    &spectrum...skipping analysis')")
   analyze_cis=.false.
  else

! Read the eigenvalues of the occupied orbitals, both alpha and beta
! (if a spin-polarised calculation)
   close(io_out1)
   open(unit=io_out1,file='hf_spectrum.dat',status='unknown')
   do ispin=1,ispin_lim
    ie=0
    if(ispin==1)then
     call getline(io_in,'Alpha  occ. eigenvalues',fields,NF,NF_max,ifail,.true.)
    else
     call getline(io_in,'Beta  occ. eigenvalues',fields,NF,NF_max,ifail,.true.)
    endif
    if(ifail<0)exit

    done=.false.
    do while(.not.done)
     do i=5,NF
      ie=ie+1
      read(fields(i),*,iostat=icheck)hf_spectrum(ie,ispin)
      if(icheck/=0)then
       ie=ie-1
       done=.true.
       exit
      endif
     enddo
     if(ispin==1)then
      call getline(io_in,'Alpha ',fields,NF,NF_max,ifail,.true.)
     else
      call getline(io_in,'Beta ',fields,NF,NF_max,ifail,.true.)
     endif
     if(ifail<0)done=.true.
    enddo
! The number of HF eigenvalues successfully read from output file - will
! be needed in sum_degen_excite.f90
    Nmo_read(ispin)=ie
    write(io_out1,fmt="('Spin ',i1,' HF eigenvalues:')")ispin
    write(io_out1,*)(i,hf_spectrum(i,ispin),i=1,ie)
   enddo
   if(ie==0)then
    write(6,fmt="('Failed to read eigenvalue spectrum...skipping analysis')")
    analyze_cis=.false.
   endif

  endif

 endif

 close(io_out1)
 close(io_out2)
 close(unit=io_in)

END SUBROUTINE read_G9xout

SUBROUTINE user_control(test)
!----------------------------------------------------------------------------!
! Ask the user which Gaussian job they want us to convert, read it, and then !
! act according to its contents and the user's responses.                    !
!----------------------------------------------------------------------------!
 USE cis_data
 USE g94_wavefunction
 IMPLICIT NONE
 LOGICAL,INTENT(inout) :: test
 CHARACTER(1) instring

 write(*,fmt="(/'Gaussian94/98/03 -> QMC wave function conversion utility.'/&
  &/'This program assumes that the Gaussian job file, output file'/&
  &'and formatted checkpoint file are all in the current directory'/&
  &'and share the same root name.  The suffixes of the output and'/&
  &'checkpoint files are assumed to be ''.out'' and ''.fchk'' (or ''.Fchk''),&
  & respectively.'/'The resulting ''gwfn.data'' file is in &
  &''[root name].qmc''.')")
 write(*,fmt="(/'Please enter the root of the Gaussian94/98/03 input file &
  &name:')")
 read(*,fmt="(A30)") g94_file
! g94_file="triplet1_slater"

! Obtain from the output file the version of Gaussian (94,98,03) that was used
 call get_gauss_version

! Read-in the wave function parameters from the formatted checkpoint file
 call read_fchk

! Read the output file produced by Gaussian94 to get the nuclear-nuclear
! repulsion energy.  If we are interested in the CIS wave function then
! also read in the configurations defining the excited states together
! with the associated CIS expansion coefficients
 call read_G9xout

! Since we may be dealing with a relatively incomplete CIS expansion
! (depending on the tolerance used in Gaussian to identify the configs
! to output) we need to ensure that the CIS wave function is properly
! normalized.  This is not necessary for QMC but it keeps things tidy.
 if(CIS)call normalize_ci

 if(CAS)call cas_wfn

 if(CIS)then
  if(maxval(Nconfig(:))>1)then
   write(*,fmt="(/'Do you want to re-sum the determinants in the &
    &expansion of this state (y/n)?'/'(You will also &
    &need to do this if you wish to perform an &
    &analysis'/'of the composition of the excited state.)')")
   read(*,fmt="(a1)")instring

   if((instring=='y').or.(instring=='Y'))then
    resummed=.true.
    write(*,fmt="(/'Resumming excited state ',i2,'...')")icis_out
    call re_sum
   endif

! re_sum needs to have been called because it calculates the number of
! excitations out of each occupied MO
   if(analyze_cis.and.resummed)call analyze_cis_state
  endif

 elseif(CAS)then
  write(*,fmt="(/'Outputting the CASSCF wave function...')")
 else
! We only have the ground state so that's the one we send to the QMC
! input file
  write(*,fmt="(/'Outputting the ground-state wave function...')")
  icis_out=0
 endif

! The initialisation of test in g9xtoqmc determines whether or not
! we ask the user if they want to test the MOs
 if(test)then
  write(*,fmt="(/'Do you wish to test the normalization of/plot the &
   &molecular orbitals (y/n)?')")
  read(*,fmt="(a1)") instring
  if((instring=='y').or.(instring=='Y'))then
   test=.true.
  else
   test=.false.
  endif
 endif

END SUBROUTINE user_control

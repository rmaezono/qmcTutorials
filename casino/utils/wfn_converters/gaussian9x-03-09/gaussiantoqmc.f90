PROGRAM gauss_wfn
 IMPLICIT none
!--------------------------------------------------------------------!
! Utility to read the results of a Gaussian94/98/03 calculation and  !
! output the resulting wavefunction in a form compatible with        !
! CASINO. The code REQUIRES the existence of a formatted             !
! checkpoint file (produced by putting FormCheck=(MO,Basis) in       !
! the route section of the Gaussian job file).  It expects this      !
! file to have a '.Fchk' suffix.  The output file of the Gaussian    !
! job is also REQUIRED.  It is assumed that this has a '.out'        !
! suffix.  If the  original Gaussian job file is present then it     !
! will be appended to the end of the CASINO gwfn.data file.          !
!                                                                    !
! ARP 2001                                                           !
!                                                                    !
! Changes                                                            !
! -------                                                            !
! 5.2003  MDT - Rewrite and tidy. New Makefile.                      !
! 10.2003 MDT - Added Gaussian03 support.                            !
!--------------------------------------------------------------------!

 LOGICAL :: test=.false.! Flag to plot out & possibly
                         ! test the normalization of MOs.  User
                         ! is only given the choice if this is
                         ! set to true here.

! Ask the user for the name of the job file and read-in the
! wave function from the corresponding .Fchk file.  If the user
! wants a CIS wave function, read in the coefficients from the
! .out file.
  call user_control(test)

! Additional bit for testing
  if(test)call wfn_test

! Create an input file for CASINO using this G94/98/03 wave function
  call qmc_write

END PROGRAM gauss_wfn

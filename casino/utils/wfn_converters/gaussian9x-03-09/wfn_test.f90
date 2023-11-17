SUBROUTINE wfn_test
!---------------------------------------------------------------------------!
! Ask the user which function/wave function to plot out and whether or      !
! not to test its normalization.                                            !
!---------------------------------------------------------------------------!
 USE g94_wavefunction, ONLY: Nmo,Nspin,SPIN
 USE integ_params
 IMPLICIT NONE
 INTEGER mo_plot
 CHARACTER(1) instring

! Ask the user which MO to plot
 if(SPIN)then
  write(*,fmt="(/'We have ',i3,' alpha and ',i3,' beta MOs of which ',i2/&
   &'and ',i2,' are occupied, respectively.')")Nmo,Nmo,Nspin
 else
  write(*,fmt="(/'Spin restricted calculation. We have ',i3,' alpha MOs&
   & of which ',i2/'are occupied.')")Nmo,Nspin(1)
 endif

 write(*,fmt="(/'The alpha MOs come before the beta MOs.  Which MO'/&
  &'do you want to plot/test the normalization of?')")
 READ (*,*) mo_plot

! Since MOs are referenced according to their spin
! we have to shift their index, hence in general mo_ref /= mo_plot
 if(mo_plot<=Nmo)then
! Alpha MO
  mo_spin=1
  mo_ref=mo_plot
 else
  mo_spin=2
  mo_ref=mo_plot-Nmo
 endif

 write(*,fmt="(/'Test the normalization of this wave function (y/n)?')")
 read(*,fmt="(a1)") instring

! Perform crude integration to check normalization of this function
 if((instring=='y').or.(instring=='Y'))call normalization_check

! Reconstruct the single MO for plotting.
 call wfn_construct

END SUBROUTINE wfn_test

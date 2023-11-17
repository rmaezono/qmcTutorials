SUBROUTINE cas_wfn
!--------------------------------------------------------------------------!
! Manipulate the CAS wave function in order to make it easier to  output.  !
!--------------------------------------------------------------------------!
 USE cis_data
 USE g94_wavefunction
 IMPLICIT none

 INTEGER icount,isign,idet,imo,iref_up(Nspin_cas(1)),iref_dn(Nspin_cas(2)), &
  &iconfig_up(Nspin_cas(1)),iconfig_dn(Nspin_cas(2))
 REAL(KIND=dp) sumsq
 CHARACTER(1) instring
 CHARACTER(2) tmpstring
 CHARACTER(60) fmtstring(2)
 LOGICAL free_mo(Ncoeffs) ! free_mo(imo)=.true. in MO imo is not involved in
                          ! any CAS configurations

! Construct the format string for printing out each of the configurations in
! the CAS expansion, and that for the reference configuration
 write(fmtstring(1),fmt="(i2)")Nspin_cas(1)
 if(Nspin_cas(2)>0)then
  write(tmpstring,fmt="(I2)")Nspin_cas(2)
  fmtstring(2)="('Config. ',I2,': [',"//trim(adjustl(fmtstring(1)))//&
   &"i2,'] [',"//trim(adjustl(tmpstring))//"i2,'], Isign = ',i2)"
  fmtstring(1)="('Reference state: [',"//trim(adjustl(fmtstring(1)))//&
   &"i2,'] [',"//trim(adjustl(tmpstring))//"i2,']')"
 else
  fmtstring(2)="('Config. ',I2,': [',"//trim(adjustl(fmtstring(1)))//&
   &"i2,'] [ ], Isign = ',I2)"
  fmtstring(1)="('Reference state: [',"//trim(adjustl(fmtstring(1)))//&
   &"i2,'] [ ]')"
 endif

! Bring each of the configurations into maximum coincidence with
! the reference configuration
 iref_up=icas_ref(1:Nspin_cas(1):1)
 iref_dn=icas_ref(Nspin_cas(1)+1:Nact_elecs:1)

! For outputting details of the CAS configurations
!  write(77,fmt=fmtstring(1))icas_ref(:)

 sumsq=0.d0
! Initialise all elements of free_mo to .true.
 free_mo=.true.
 do idet=1,Ncas_basis,1
  isign=1
  iconfig_up=icasdet(1:Nspin_cas(1),idet)
  iconfig_dn=icasdet(Nspin_cas(1)+1:Nact_elecs,idet)

! Spin-up determinant
  call max_coincidence(iref_up,iconfig_up,Nspin_cas(1),isign)

! Spin-down determinant
  call max_coincidence(iref_dn,iconfig_dn,Nspin_cas(2),isign)

! Apply isign to determinant coefficient
  cas_coeffs(idet)=real(isign,dp)*cas_coeffs(idet)

! Restore configuration
  icasdet(1:Nspin_cas(1),idet)=iconfig_up
  icasdet(Nspin_cas(1)+1:Nact_elecs,idet)=iconfig_dn

! For outputting details of the CAS configurations
!  write(77,fmt=trim(fmtstring(2)))idet,icasdet(:,idet),isign

! Check the normalization of the expansion coefficients
  sumsq=sumsq+cas_coeffs(idet)*cas_coeffs(idet)

! Search for MOs that are not involved in any of the CAS
! configurations.  These will be needed later in order to store
! resummed MOs.
  free_mo(icasdet(1,idet):icasdet(Nact_elecs,idet))=.false.
 enddo

! Now create an array that lists each of the MOs that are free for us
! to store resummed MOs in
 icount=0
 do imo=1,Nmo
  if(free_mo(imo))then
   icount=icount+1
   list_of_free_mos(icount)=imo
  endif
 enddo
! Finally, store the number of MOs that are available for storage
 Nfree_mo=icount-1

 write(*,fmt="(/'Sum of squares of CAS coefficients = ',F9.5)")sumsq

 write(*,fmt="(/'Do you want to re-sum the determinants in the &
  &expansion of this state (y/n)?')")
 read(*,fmt="(A1)") instring
! instring="y"

 if((instring=='y').or.(instring=='Y'))then
  resummed=.true.
  call resum_cas
 else
  resummed=.false.
 endif

END SUBROUTINE cas_wfn

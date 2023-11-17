SUBROUTINE analyze_cis_state
!---------------------------------------------------------------------!
! Take a CIS excited state and examine the MOs that are excited from. !
! For those that are degenerate, add the percentage weights for       !
! excitations into each final MO.                                     !
!---------------------------------------------------------------------!
 USE cis_data
 USE g94_wavefunction
 IMPLICIT NONE
 INTEGER ishift,imo,ispin
 INTEGER,PARAMETER :: max_degen=4
 REAL(KIND=dp) diff
! Tolerance for deciding whether two MOs are degenerate or not
 REAL(KIND=dp),PARAMETER :: degen_tol=1.d-6

 imo=1
 do ispin=1,ispin_lim

  do while(imo<=nspin(ispin))

   if(npromote(imo,ispin)>0)then
    ishift=1
    do
     diff=abs(hf_spectrum(imo+ishift,ispin)-hf_spectrum(imo,ispin))
     if(diff<degen_tol)then
      ishift=ishift+1
      if(ishift>max_degen)then
       write(*,fmt="('A degeneracy greater than ',i1,' appears &
        &to be present!')")max_degen
        stop
      elseif((imo+ishift)>nspin(ispin))then
! have reached the end of all mos of spin ispin so the
! current set imo -> imo+ishift-1 must be degenerate
       call sum_degen_excite(degen_tol,max_degen,imo,imo+ishift-1,ispin)
       imo=imo+ishift
       exit
      endif
      cycle
     else
! the current state is not degenerate with the previous ones
! so have reached end of current set of degenerate orbitals
      call sum_degen_excite(degen_tol,max_degen,imo,imo+ishift-1,ispin)
      imo=imo+ishift
      exit
     endif
    enddo
   else
    imo=imo+1
   endif
  enddo
 enddo

END SUBROUTINE analyze_cis_state

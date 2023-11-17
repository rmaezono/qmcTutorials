SUBROUTINE cas_write(iounit)
!-----------------------------------------------------------------!
! Output the CASSCF wave function according to the rules for the  !
! multi-determinant specification section of the gwfn.data file.  !
!-----------------------------------------------------------------!
 USE cis_data
 IMPLICIT none

 INTEGER,INTENT(in) :: iounit ! Unit number of output file
 INTEGER idet,imo,icount
! INTEGER iSz ! z component of spin projection
 INTEGER num_det ! Number of determinants left after resumming

 write(iounit,fmt="('MD')")

 if(slater)then ! we have standard Slater determinants in the expansion - easier

  if(.not.resummed)then
! Number of determinants in wfn expansion
   write(iounit,*)Ncas_basis

! Expansion coefficients
   do idet=1,Ncas_basis
    write(iounit,*)cas_coeffs(idet)
   enddo

! Define the excitations
   do idet=1,Ncas_basis
! Spin-up
    do imo=1,Nspin_cas(1)
     if(icasdet(imo,idet)/=icas_ref(imo))then
! MO is not the same as that in reference state; must be an excitation
      write(iounit,fmt="('DET ',i3,1x,i1,1x,'PR',2(1x,i3,1x,i1))") &
       &idet,1,icas_ref(imo),1,icasdet(imo,idet),1
     endif
    enddo
! Spin-down
    do imo=(Nspin_cas(1)+1),Nact_elecs
     if(icasdet(imo,idet)/=icas_ref(imo))then
      write(iounit,fmt="('DET ',i3,1x,i1,1x,'PR',2(1x,i3,1x,i1))") &
       &idet,2,icas_ref(imo),1,icasdet(imo,idet),1
     endif
    enddo
   enddo

  else
! Expansion has been (partially) resummed

   num_det=Ncas_basis-num_det_deleted
! Number of determinants in wfn expansion is thus reduced
   write(iounit,*)num_det

! Expansion coefficients
   do idet=1,Ncas_basis
    if(.not.deleted_det(idet))then
     if(resummed_det(idet))then
! Expansion coefficients are absorbed into the resummed MO
      write(iounit,*)1.d0
     else
! This determinant remained unaffected by resumming so output its coefficient
      write(iounit,*)cas_coeffs(idet)
     endif
    endif
   enddo

! Define the excitations
   icount = 0
   do idet=1,Ncas_basis

    if(.not.deleted_det(idet))then
     icount=icount+1
! Spin-up
     do imo=1,Nspin_cas(1)
      if(icasdet(imo,idet)/=icas_ref(imo))then
! MO is not the same as that in reference state; must be an excitation
       write(iounit,fmt="('DET ',i3,1x,i1,1x,'PR',2(1x,i3,1x,i1))") &
        &icount,1,icas_ref(imo),1,icasdet(imo,idet),1
      endif
     enddo
! Spin-down
     do imo=(Nspin_cas(1)+1),Nact_elecs
      if(icasdet(imo,idet)/=icas_ref(imo))then
       write(iounit,fmt="('DET ',i3,1x,i1,1x,'PR',2(1x,i3,1x,i1))") &
        &icount,2,icas_ref(imo),1,icasdet(imo,idet),1
      endif
     enddo
    endif

   enddo
  endif

 else

  write(*,fmt="(/'I''m afraid that I have yet to work out how to express&
   & Gaussian''s spin'/'configurations in terms of Slater determinants and&
   & therefore I''d suggest'/'you re-run this CAS calculation with the&
   & ''SlaterDet'' or IOp(4/46=3) option.')")
   stop

!!$  ! CAS expansion consists of spin configurations therefore we have to
!!$  ! take into account the multiplicity of the wave function
!!$  if(multiplicity==3)then
!!$  ! Calculate z spin projection in reference state
!!$   iSz=abs(Nspin_cas(1)-Nspin_cas(2))+1
!!$
!!$   if(iSz/=0)then
!!$  ! Triplet wave functions with Sz /= 0 consist of a single
!!$  ! determinant so we can treat the spin configuration as though it
!!$  ! were a Slater determinant.
!!$
!!$  ! Number of determinants in wfn expansion
!!$    write(iounit,*)Ncas_basis
!!$
!!$  ! Expansion coefficients
!!$    do idet=1,Ncas_basis
!!$     write(iounit,*) cas_coeffs(idet)
!!$    enddo
!!$
!!$    do idet=1,Ncas_basis
!!$  ! Spin-up
!!$     do imo=1,Nspin_cas(1)
!!$      if(icasdet(imo,idet)/=icas_ref(imo))then
!!$  ! MO is not the same as that in reference state; must be an excitation
!!$       write(iounit,fmt="('DET ',i3,1x,i1,1x,'PR',2(1x,i3,1x,i1))") &
!!$        &idet,1,icas_ref(imo),1,icasdet(imo,idet),1
!!$      endif
!!$     enddo
!!$  ! Spin-down
!!$     do imo=(Nspin_cas(1)+1),Nact_elecs
!!$      if(icasdet(imo,idet)/=icas_ref(imo))then
!!$       write(iounit,fmt="('DET ',i3,1x,i1,1x,'PR',2(1x,i3,1x,i1))") &
!!$        &idet,2,icas_ref(imo),1,icasdet(imo,idet),1
!!$      endif
!!$     enddo
!!$    enddo
!!$
!!$   else
!!$  ! Triplet wave functions for a single excitation consist of
!!$  ! two determinants; spin-down excitation minus equivalent
!!$  ! spin-up excitation.  However, for double, triple etc.
!!$  ! excitations life is hard and I haven't worked it out yet...
!!$    write(*,fmt="('Ah - we have a triplet state of non-zero &
!!$     &multiplicity expressed in terms of spin configurations.'/&
!!$     &'I have yet to work out how to convert this to Slater &
!!$     &determinants so my advice would be to rerun the CAS with &
!!$     &the ''SlaterDet'' option.')")
!!$    stop
!!$
!!$   endif
!!$
!!$  elseif(multiplicity==1)then
!!$  ! Singlet state - have two Slater determinants for each spin
!!$  ! configuration in the expansion for single excitations but
!!$  ! I haven't worked out the details for triples and higher.
!!$  ! Hence we chicken out at this point...
!!$
!!$   write(*,FMT="(/'Ah - a singlet state in terms of spin &
!!$    &configurations eh? My advice to'/'you would be to rerun &
!!$    &that CAS calculation with the ''SlaterDet'' option.')")
!!$   stop
!!$
!!$  else
!!$   write(*,fmt="('Unsupported multiplicity, Multiplicity = ',i1)") &
!!$    &Multiplicity
!!$   stop
!!$  endif
 endif

END SUBROUTINE cas_write

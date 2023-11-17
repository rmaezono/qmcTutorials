SUBROUTINE sum_degen_excite(degen_tol,max_degen,istart,ifinish,ispin)
!--------------------------------------------------------------------------!
! Read weights for excitations out of each MO from istart to ifinish, sum  !
! them, and output to a new file.                                          !
!--------------------------------------------------------------------------!
 USE cis_data
 USE g94_wavefunction
 IMPLICIT NONE
 INTEGER,INTENT(in) :: max_degen ! Range of occ. MOs that are degenerate
 INTEGER,INTENT(in) :: ispin ! Their spin (1=up, 2=down)
 INTEGER,INTENT(in) :: istart,ifinish
 REAL(KIND=dp),INTENT(in) :: degen_tol
 INTEGER icheck,ishift,jmo,ifinal
 INTEGER :: io_out=25
 REAL(KIND=dp) temp_weight,diff,wght_to_percent
 REAL(KIND=dp),PARAMETER :: tol_zero=1.d-12
! weight(i) holds percentage weight of excitation into MO i
 REAL(KIND=dp),ALLOCATABLE :: weight(:)
 CHARACTER(30) :: tmpstring,filename
 CHARACTER(50) :: fmt_string

 allocate(weight(Nmo),stat=icheck)
 if(icheck<0)then
  write(*,fmt="('Failed to allocate memory in sum_degen_excite')")
  return
 endif

 weight=0.d0
! Conversion factor to give percentage contribution of given state to
! CIS expansion (allows for the fact that the latter might not be complete)
 wght_to_percent=100.d0/ci_norm
 do jmo=istart,ifinish
! For excitations from MO jmo, read weights associated with
! each MO excited into and add to weights for the other degenerate MOs
  write(tmpstring,*)jmo
  filename='from'//trim(adjustl(tmpstring))//"."
  write(tmpstring,*)ispin
  filename=trim(filename)//trim(adjustl(tmpstring))//".dat"
  open(unit=io_out,file=trim(adjustl(filename)),status='old')
  do
   read(io_out,*,iostat=icheck)ifinal,temp_weight
   if(icheck==0)then
    weight(ifinal)=weight(ifinal)+wght_to_percent*temp_weight
   else
    exit
   endif
  enddo
  close(io_out)
 enddo

! Output the summed weights for excitation out of istart->ifinish
! However, some of the MOs excited into will be degenerate too so
! need to check for this and add the weights of any that are
 write(tmpstring,*) istart
 filename='from'//trim(adjustl(tmpstring))//"_"
 write(tmpstring,*) ifinish
 filename=trim(filename)//trim(adjustl(tmpstring))//"."
 write(tmpstring,*) ispin
 filename=trim(filename)//trim(adjustl(tmpstring))//".dat"
 write(*,fmt="('Writing summed weights to ',(a))") filename
 open(unit=io_out,file=trim(filename),status='unknown')
 write(io_out,fmt="('# Sum of excitations out of MOs ',i3,' to ',i3)") &
  &istart,ifinish
 write(io_out,fmt="('#',15x,'Final MO')")
 write(io_out,fmt="('# Degeneracy  Number (LUMO=1)  Energy rel. HOMO (eV)&
  &       Weight')")
 fmt_string="(4x,i1,10x,i4,14x,f9.5,13x,f11.7)"

 jmo=Nocc_mo+1
 do while(jmo<Nmo_read(ispin))

  if(weight(jmo)>tol_zero)then
   ishift=1
   do
    diff=abs(hf_spectrum(jmo+ishift,ispin)-hf_spectrum(jmo,ispin))
    if(diff<degen_tol)then
     ishift=ishift+1
     if(ishift>max_degen)then
      write(*,fmt="('A degeneracy greater than ',i1,' appears &
       &to be present (sum_degen_excite)...')")max_degen
      write(*,fmt="('...skipping analysis of CIS wave function.')")
      write(*,fmt="('jmo = ',i3)")jmo
      return
     elseif((jmo+ishift)>Nmo)then
! Have reached the end of all MOs so the current set jmo -> jmo+ishift-1
! must be degenerate
      write(io_out,fmt=trim(adjustl(fmt_string))) ishift,jmo-Nocc_mo, &
       &eV*(hf_spectrum(jmo,ispin)-hf_spectrum(Nocc_mo,ispin)),&
       &sum(weight(jmo:jmo+ishift-1))

      jmo=jmo+ishift
      exit
     endif

     cycle
    else
! The current state is not degenerate with the previous ones
! so have reached end of current set of degenerate orbitals
     write(io_out,fmt=trim(adjustl(fmt_string))) ishift, jmo-Nocc_mo, &
      &eV*(hf_spectrum(jmo,ispin)-hf_spectrum(Nocc_mo,ispin)),&
      &sum(weight(jmo:jmo+ishift-1))

     jmo=jmo+ishift
     exit
    endif
   enddo
  else
   jmo=jmo+1
  endif
 enddo

 close(io_out)

END SUBROUTINE sum_degen_excite

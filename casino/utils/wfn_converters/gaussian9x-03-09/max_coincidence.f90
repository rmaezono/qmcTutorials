SUBROUTINE max_coincidence(iconfig_ref,iconfig,nlen,isign)
!------------------------------------------------------------------!
! Subroutine to bring determinant (list of orbitals) iconfig into  !
! maximum coincidence with the determinant stored in iconfig_ref.  !
! Keeps track of sign changes associated with swapping rows of a   !
! determinant.                                                     !
!------------------------------------------------------------------!
 IMPLICIT none
 INTEGER,INTENT(in) :: nlen
 INTEGER,INTENT(in) :: iconfig_ref(nlen)
 INTEGER,INTENT(inout) :: iconfig(nlen)
 INTEGER,INTENT(inout) :: isign
 INTEGER i,j,iref_orb,itemp

 do i=1,nlen
  iref_orb=iconfig_ref(i)
  do j=1,nlen,1

   if((iconfig(j)==iref_orb).and.(i/=j))then
! Swap the orbitals at positions i and j in iconfig()
    itemp=iconfig(j)
    iconfig(j)=iconfig(i)
    iconfig(i)=itemp
! Account for the sign change
    isign=-1*isign
    exit ! As these are spin orbitals they can only be singly
         ! occupied so now we've found it we can stop looking.
   endif

  enddo
 enddo

END SUBROUTINE max_coincidence

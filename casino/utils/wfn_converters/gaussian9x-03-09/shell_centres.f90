SUBROUTINE shell_centres(icentre,ifirst_shll)
!---------------------------------------------------------------------!
! Identify the positions of the unique shell centres and store the    !
! first shell index corresponding to each. (This assumes that         !
! Gaussian94 outputs them in some sort of order rather than jumping   !
! back and forth between centres.)                                    !
!---------------------------------------------------------------------!
 USE g94_wavefunction
 IMPLICIT NONE
 INTEGER,INTENT(out) :: icentre ! Returns the number of unique shell centres
 INTEGER,INTENT(out) :: ifirst_shll(Nshells+1) ! Index of the first shell at
                                               ! each unique centre
 INTEGER ishell,iuniq

 REAL(KIND=dp) centre(3,Nshells) ! Holds the unique shell centres
 REAL(KIND=dp) diff(3) ! Diff between current shell centre and
                       !previously stored centre
 LOGICAL newcentre ! Flag for the discovery of a new shell centre

 icentre=1
 centre(:,1)=shll_posn(:,1)
 ifirst_shll(1)=1

 do ishell=2,Nshells

  newcentre=.true.

  do iuniq=1,icentre
   diff=abs(centre(:,iuniq) - shll_posn(:,ishell))
   if(sum(diff)<1.d-6)newcentre=.false. ! This is not a new centre
  enddo

  if(newcentre)then
   icentre=icentre+1
   centre(:,icentre)=shll_posn(:,ishell)
   ifirst_shll(icentre)=ishell
  endif

 enddo

! False data to allow loops of form:
! do n=1,num_centres
!  do shell=first_shell(n),first_shell(n+1)-1
 ifirst_shll(icentre+1)=Nshells+1

!DBG
! write(*,fmt="('There are ',i2,' unique shell centres at:')") icentre
! write(*,fmt="(3(1pe20.13))") centre(:,1:icentre)
!DBG
END SUBROUTINE shell_centres

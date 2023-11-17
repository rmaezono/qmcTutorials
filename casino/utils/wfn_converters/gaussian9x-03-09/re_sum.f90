SUBROUTINE re_sum
!---------------------------------------------------------------------!
! Re-sum the determinants in the expansion of the CIS/TD-DFT          !
! excited state to be output, icis_out, so that there are only as     !
! many determinants as there are distinct occupied orbitals excited   !
! from.                                                               !
!---------------------------------------------------------------------!
 USE cis_data
 USE g94_wavefunction
 IMPLICIT NONE

 INTEGER iconfig ! Loops over configurations within excited state icis_out
 INTEGER ivirt ! Stores current virtual MO being excited into during resumming
 INTEGER Nvirtual ! Index for storage of extra virtual MOs created by resumming
 INTEGER iocc ! The occupied MO being excited from in the current configuration
 INTEGER ipromo_count ! Counts the no. of promotions from each occupied orbital
 INTEGER iex ! Loops over excitations from a given MO
 INTEGER memtest ! Check on success of memory allocation
 INTEGER ispin ! Loops over spin - 1=alpha, 2=beta
 REAL(KIND=dp) :: ci1 ! Temporary variable to hold CIS expansion coefficients

 iconfig=max(Nspin(1),Nspin(2)) ! Acts as a temporary variable here
 allocate(Npromote(iconfig,2),stat=memtest )
 if(memtest/=0)stop 'Could not allocate memory in re_sum.'

 do ispin=1,ispin_lim

  Npromote(:,ispin)=0 ! Initialise the storage of how many excitations there
                      ! are from each occupied MO of spin `ispin'

  if(Nconfig(ispin)>0)then
   Ndet(ispin)=1
   iocc=Orbitals(1,1,ispin) ! The first configuration is an
   ipromo_count=1           ! excitation from this state - count it
  else
   Ndet(ispin)=0
   ipromo_count=0 ! There are no configs of this spin contributing to icis_out
   Npromote(:,ispin)=0
  endif

  do iconfig=2,Nconfig(ispin)

! Gaussian outputs the excitations from one occupied MO consecutively
! and it also loops over the occupied orbitals in order too.  This makes
! it easy to find the next (different) occ. orbital to be excited from:
   if(Orbitals(1,iconfig,ispin)/=iocc)then
! This configuration involves an excitation from a different MO
! so store the no. of excitations from the last MO
    Npromote(iocc,ispin)=ipromo_count
! Increment the no. of occupied MOs that are excited from to form
! the icis_out excited state
    Ndet(ispin)=Ndet(ispin)+1
! Update which occupied MO we are currently exciting from
    iocc=Orbitals(1,iconfig,ispin)
! Initialise the count of the no. of excitations from this MO
    ipromo_count=1
   else
! Another excitation from the iocc'th MO
    ipromo_count=ipromo_count+1
   endif

  enddo
! That's all of the excitations so store the count of the no. from the
! last distinct occupied MO
  if(Nconfig(ispin)>0)Npromote(iocc,ispin)=ipromo_count

  if(ispin==1)then
   write(*,fmt="(/'There are ',i3,' resummed alpha determinants:')")Ndet(1)
   write(*,fmt="((i4,' excitations from alpha orbital ',i2))") &
    &(Npromote(iocc,1),iocc,iocc=1,Nspin(1))
  else
   if(Nconfig(ispin)>0)then
    write(*,fmt="(/'There are ',i3,' resummed beta determinants:')")Ndet(2)
    write(*,fmt="((i4,' excitations from beta orbital ',i2))") &
     &(Npromote(iocc,2),iocc,iocc=1,Nspin(2))
   else
    write(*,fmt="(/'There are no resummed beta determinants.')")
   endif
  endif
  write(*,*)

! Now do the resumming itself

  iconfig=0 ! Configuration index - used to reference CIS expansion coeff.
            ! and the virtual orbital being excited into
! Nvirtual=Ncoeffs ! Set index to the current no. of orbitals so that
                   ! those created by resumming are stored after them
  Nvirtual=Nmo

! Loop over occupied orbitals
  do iocc=1,Nspin(ispin)

   if(Npromote(iocc,ispin)>0)then
    Nvirtual=Nvirtual+1 ! Create a new resummed determinant

! Loop over excitations from the current occ. orbital
    do iex=1,Npromote(iocc,ispin)
     iconfig=iconfig+1
     ivirt=Orbitals(2,iconfig,ispin)
     ci1=ci_coeff(iconfig,ispin)
! Loop over basis functions - adding to evcoeff1 is OK since
! it is initialised to zero in read_fchk
     evcoeff1(:,Nvirtual,ispin)=evcoeff1(:,Nvirtual,ispin) &
      &+ci1*evcoeff1(:,ivirt,ispin)
    enddo

   endif
  enddo

! Overwrite the first Ndet virtual orbitals with these resummed orbitals.
! This is permitted since CASINO will not attempt any further
! promotions or subtractions with the determinants given it for a
! multideterminant calculation.
  Nvirtual=Nmo
  do ivirt=(Nspin(ispin)+1),(Nspin(ispin)+Ndet(ispin))
   Nvirtual=Nvirtual+1
   evcoeff1(:,ivirt,ispin)=evcoeff1(:,Nvirtual,ispin)
  enddo

 enddo

END SUBROUTINE re_sum

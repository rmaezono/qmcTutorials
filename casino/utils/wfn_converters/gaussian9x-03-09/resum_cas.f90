SUBROUTINE resum_cas
!-------------------------------!
! Resum the CAS wave function.  !
! Last modified, 28/9/2000.     !
!-------------------------------!
 USE cis_data
 USE g94_wavefunction
 USE rejig
 IMPLICIT NONE

! G98 always does ROHF-based CASSCF calculations and therefore we
! only have to deal with spin-up eigenvectors - these are all that
! are stored.  ispin_rohf ensures that it is these eigenvectors
! that we resum, independently of the spin of the determinant.
 INTEGER,PARAMETER :: ispin_rohf=1
 INTEGER ispin,isign
! Loop counters
 INTEGER idet,jdet,kdet,imo,jmo,ig,imem
! Control loop limits according to spin (for shorter source code)
 INTEGER lim,lim1,lim2,mo_shift
 INTEGER num_excite
! Counter for storing resummed MOs
 INTEGER ivirt,ifree_mo
 INTEGER num_del ! Number of deleted dets.
 INTEGER num_remain ! Number of remaining dets.
! Ngroups(idet,ispin) - number of different resummable groups associated
! with determinant idet with spin part ispin factored out
 INTEGER Ngroups(Ncas_basis,2)
! num_members(ig,idet,ispin) - number of members in group ig, of
! excitations from idet, when spin ispin part factors out
 INTEGER num_members(Ncas_basis,Ncas_basis,2)
! members(i,ig,idet,ispin) - identity of member i, of group ig for
! excitations from idet with ispin part common to both dets
 INTEGER members(Ncas_basis,Ncas_basis,Ncas_basis,2)
! Keep a count of determinants resummed and removed
 INTEGER icount_resum
! List of dets that have been resummed and removed
 INTEGER resummed_list(Ncas_basis)
! Temp. space for bringing resummable dets into max. coincidence
 INTEGER iconfig(Nact_elecs)
 INTEGER ncommon ! No. of MOs in common between two dets
 INTEGER mo_common(Nact_elecs,Ncas_basis)
 INTEGER ibig, ibig_grp
! Store which det (& spin) we will resum next
 INTEGER idet_resum, ispin_resum
! Position of MO to be resummed in list making up icasdet
 INTEGER imo_resum
! Spin of the determinants that factor out of resummation
 INTEGER ispin_cmmn
! No. of dets to resum in current group: used as check to see that it's not zero
 INTEGER num_to_sum
! Used to test success of memory allocation
 INTEGER icheck
 REAL(KIND=dp) ci1 ! Temporary copy of Gaussian coefficient in an MO
 LOGICAL disagree, excite, present
! True if this member has been deleted due to its already having
! been resummed. Initialised to .false. below.
 LOGICAL deleted_member(Ncas_basis,Ncas_basis,Ncas_basis,2)
! For error-checking output
 CHARACTER(2) tmpstring
 CHARACTER(60) fmtstring(2)

! Construct format strings for error-checking output
 write(fmtstring(1),fmt="(i2)") Nspin_cas(1)
 if(Nspin_cas(2)>0)then
  write(tmpstring,fmt="(i2)") Nspin_cas(2)
  fmtstring(2)="('Config. ',i2,': [',"//trim(adjustl(fmtstring(1)))//&
   &"i2,'] [',"//trim(adjustl(tmpstring))//"i2,'], Isign = ',i2)"

  fmtstring(1)="('Reference state: [',"//trim(adjustl(fmtstring(1)))//&
   &"i2,'] [',"//trim(adjustl(tmpstring))//"i2,']')"
 else
  fmtstring(2)="('Config. ',i2,': [',"//trim(adjustl(fmtstring(1)))//&
   &"i2,'] [ ], Isign = ',i2)"

  fmtstring(1)="('Reference state: [',"//trim(adjustl(fmtstring(1)))//&
   &"i2,'] [ ]')"
 endif

 allocate(resummed_det(Ncas_basis), deleted_det(Ncas_basis), STAT=icheck)
 if(icheck/=0)call fatal('Could not allocate memory in resum_cas.')

! Initialise excitation counters etc. to zero
 Ngroups=0 ; Num_members=0
 icount_resum=0
 resummed_list=0 ; num_det_deleted=0
 deleted_member=.false. ; resummed_det=.false.
 deleted_det=.false.

! Resummed MOs will be stored in evcoeff(:,ifree_mo,ispin_rohf) where
! ifree_mo identifies a spare MO as listed in list_of_free_mos and CAS
! is ROHF so ispin_rohf=1.
 Nmo_resummed=0 ! Counter for resummed MOs created - used to move
                ! through list_of_free_mos

! First, we identify the possible excitations from each determinant..
 do idet=1,Ncas_basis
  do jdet=idet+1,Ncas_basis

   do ispin=1,2
    if(ispin==1)then
     mo_shift=0
    else
     mo_shift=Nspin_cas(1)
    endif

    do imo=(mo_shift+1),(mo_shift+Nspin_cas(ispin))
     disagree=.true.
     do jmo=(mo_shift+1),(mo_shift+Nspin_cas(ispin))
      if(icasdet(imo,idet)==icasdet(jmo,jdet))disagree=.false.
     enddo
     if(disagree)exit ! spin ispin parts must differ because imo
    enddo             ! not found in jdet

! If the ispin parts of idet and jdet are the same then
! look to see whether the opposite spin parts are resummable
    if(.not.disagree)then
     num_excite=0
! Alters the loop limits according to spin - limits run
! over the MOs with spin OPPOSITE to ispin
     lim=2-ispin
     lim1=1+lim*Nspin_cas(1)
     lim2=Nspin_cas(1)+lim*Nspin_cas(2)
     do imo=lim1,lim2
      excite=.true.
      do jmo=lim1,lim2
       if(icasdet(imo,idet)==icasdet(jmo,jdet))excite=.false.
      enddo
      if(excite)then
       num_excite=num_excite+1
       if(num_excite>1)exit ! Can't resum if dets differ by >1 MO
      endif
     enddo

     if(num_excite==1)then
! Identify the MOs common to both determinants (as these
! allow us to distinguish different single excitations)

! Temporarily increment the counter for groups of excitations
! from this determinant in order to give us some space whilst
! we check whether it is in a new group or not
      Ngroups(idet,ispin)=Ngroups(idet,ispin)+1
      ncommon=0
      do imo=lim1,lim2
       do jmo=lim1,lim2
        if(icasdet(imo,idet)==icasdet(jmo,jdet))then
         ncommon=ncommon+1
         mo_common(ncommon,Ngroups(idet,ispin))=icasdet(imo,idet)
        endif
       enddo
      enddo
! Sort this list of MOs that are common to the 2 dets into
! descending order for ease of comparison with others.  The
! array that can be reordered as mo_common is sorted is an
! OPTIONAL argument and is omitted here
      call numsrt(mo_common(1:ncommon,Ngroups(idet,ispin)),ncommon,Ncas_basis)

! Check whether this excitation can be resummed together with
! any of the others out of idet
! Check to see whether this excitation is similar to any of
! the others that we've found so far for idet,ispin...
      present=.false.
      do ig=1,Ngroups(idet,ispin)-1
       if(all(mo_common(1:ncommon,Ngroups(idet,ispin))==&
        &mo_common(1:ncommon,ig)))then
        present=.true.
        exit
       endif
      enddo
      if(present)then
! This excitation is in fact a member of the ig'th group so let's include it
       num_members(ig,idet,ispin)=num_members(ig,idet,ispin)+1
       members(num_members(ig,idet,ispin),ig,idet,ispin)=jdet
! ...and correct the number of groups
       Ngroups(idet,ispin)=Ngroups(idet,ispin)-1
      else
! This is a new type of excitation so initialise its membership record
       num_members(Ngroups(idet,ispin),idet,ispin)=1
       members(1,Ngroups(idet,ispin),idet,ispin)=jdet
      endif

     endif
    endif
   enddo ! Loop over spin
  enddo ! Inner loop over determinants - jdet

! Output of resummable groups of determinants for testing purposes
! if(any(Ngroups(idet,1:2)>0))then
!  write(66,fmt="(/'Det. ',i2,',   up-spin: ')")idet
!  do ig=1,Ngroups(idet,2)
!   write(66,fmt="('Group ',i2,', ',i2,' &
!    &members: ',5(i2,1x))") ig,num_members(ig,idet,2), &
!    &members(1:num_members(ig,idet,2),ig,idet,2)
!  enddo
!
!  write(66,fmt="(9x,'down-spin: ')")
!  do ig=1,Ngroups(idet,1)
!   write(66,fmt="('Group ',i2,', ',i2,' &
!    &members: ',5(i2,1x))")ig,num_members(ig,idet,1), &
!    &members(1:num_members(ig,idet,1),ig,idet,1)
!  enddo
!
! endif
 enddo

! Keep looping until all possible groups have been resummed, starting
! with the largest...
 do

! Find the biggest group...
  ibig=0
  do ispin=1,2
   do idet=1,Ncas_basis
! Check to see whether idet has already been resummed
    if(any(idet==resummed_list(1:icount_resum)))cycle
    if(Ngroups(idet,ispin)>0)then
! Check each group for whether it contains dets that have already been resummed
     do ig=1,Ngroups(idet,ispin)
      num_del=0
      do imem=1,Num_members(ig,idet,ispin)
       if(deleted_member(imem,ig,idet,ispin))then
        num_del=num_del+1
        cycle
       elseif(any(members(imem,ig,idet,ispin)==resummed_list(1:icount_resum)))&
        &then
! Member imem of group ig has already been resummed and is therefore deleted
        deleted_member(imem,ig,idet,ispin)=.true.
        num_del=num_del+1
       endif
      enddo
! See how many determinants that leaves us with
      num_remain=Num_members(ig,idet,ispin)-num_del
      if(num_remain>ibig)then
! This is the best group we have found so far - store how many dets.
! it has, which det. (and spin) it is for and which group
       ibig=num_remain
       idet_resum=idet
       ispin_cmmn=ispin
       ibig_grp=ig
      endif
     enddo

    endif
   enddo
  enddo

! If ibig=0 then the largest group contains a single determinant which
! means that it can't be resummed.  Hence we can do no more...
  if(ibig==0)exit

! ispin_cmmn is the spin of the set of determinants COMMON to all those
! about to be resummed.  Therefore, the determinants to be resummed are
! of opposite spin...
  if(ispin_cmmn==1)then
   ispin_resum=2
  else
   ispin_resum=1
  endif

! Error-checking output...
! write(78,fmt="('About to resum the following determinants...')")
! write(78,fmt=trim(fmtstring(2)))idet_resum,icasdet(:,idet_resum),1

! Bring the group into maximum coincidence with
! icasdet(:,idet_resum) where ':' loops over electrons of spin
! ispin_cmmn
  num_to_sum=0
  do ig=1,Num_members(ibig_grp,idet_resum,ispin_cmmn)

! Check to see whether this member has already been resummed
   if(deleted_member(ig,ibig_grp,idet_resum,ispin_cmmn))cycle

   kdet=members(ig,ibig_grp,idet_resum,ispin_cmmn)
   num_to_sum=num_to_sum+1
! Make temporary copy of determinant of spin ispin_resum in order
! to bring it into maximum coincidence
   lim=ispin_resum-1
   lim1=1+lim*Nspin_cas(1)
   lim2=Nspin_cas(1)+lim*Nspin_cas(2)

   iconfig(1:Nspin_cas(ispin_resum))=icasdet(lim1:lim2,kdet)
   isign=1
   call max_coincidence(icasdet(lim1:lim2,idet_resum),iconfig,&
    &Nspin_cas(ispin_resum),isign)
! Apply isign to determinant coefficient
   cas_coeffs(kdet)=real(isign,dp)*cas_coeffs(kdet)
! Restore the determinant
   icasdet(lim1:lim2,kdet)=iconfig(1:Nspin_cas(ispin_resum))
! Error-checking output...
!  write(78,fmt=trim(fmtstring(2)))kdet,icasdet(:,kdet),isign
  enddo

! Paranoid error checking
  if(num_to_sum==0)then
   call fatal('Attempting to resum no determinants in resum_cas - BUG!')
  endif

! Almost ready, just need to identify which of the MOs is to be
! resummed - compare first two determinants and find the position
! of the MOs that are not common to both of them - need some care
! as the first n of the group might already have been resummed...
  do ig=1,Num_members(ibig_grp,idet_resum,ispin_cmmn)
   if(.not.deleted_member(ig,ibig_grp,idet_resum,ispin_cmmn))exit
  enddo
! Now we have a determinant that hasn't been resummed we can do the comparison..
  do imo_resum=lim1,lim2
   if(icasdet(imo_resum,idet_resum)/=icasdet(imo_resum,kdet))exit
  enddo

!==================== Now the resumming itself =========================
! Keep a count of resummed MOs and identify the next free MO on the
! list in order to store the eigenvector we're about to create
  Nmo_resummed=Nmo_resummed+1
  if(Nmo_resummed>Nfree_mo)then
   call fatal("No. of resummed MOs > greater than no. of MOs free for&
    & storage - cannot resum wave function.")
  endif
  ifree_mo=list_of_free_mos(Nmo_resummed)
! Initialise that eigenvector ready to make new one...
  evcoeff1(:,ifree_mo,ispin_rohf)=0.d0

! First for the determinant that labels the rest of the group
! Keep a record of determinants resummed
  icount_resum=icount_resum+1
  if(icount_resum>Ncas_basis)then
   call fatal("Too many determinants resummed in resum_cas?")
  endif

  resummed_list(icount_resum)=idet_resum
! Resummed but not deleted unlike the determinants below
  resummed_det(idet_resum)=.true.
! Coefficient of the excited determinant
  ci1=cas_coeffs(idet_resum)
! The MO 'excited into' in this determinant
  ivirt=icasdet(imo_resum,idet_resum)
! Loop over basis functions - adding to evcoeff1 is OK since
! it is initialised to zero in read_fchk
  evcoeff1(:,ifree_mo,ispin_rohf)=evcoeff1(:,ifree_mo,ispin_rohf) &
   &+ci1*evcoeff1(:,ivirt,ispin_rohf)
! Now the group itself...
  do ig=1,Num_members(ibig_grp,idet_resum,ispin_cmmn)

! Check to see whether this determinant has already been resummed
   if(deleted_member(ig,ibig_grp,idet_resum,ispin_cmmn))cycle

   kdet=members(ig,ibig_grp,idet_resum,ispin_cmmn)
! Keep a record of determinants resummed
   icount_resum=icount_resum+1
   if(icount_resum>Ncas_basis)then
    write(*,*)'icount_resum=',icount_resum
    call fatal("Too many determinants resummed in resum_cas?")
   endif
   resummed_list(icount_resum)=kdet
   num_det_deleted=num_det_deleted+1
! Resummed and deleted...
   resummed_det(kdet)=.true.
   deleted_det(kdet)=.true.

   ci1=cas_coeffs(kdet)
! imo_resum holds the position of the MO being resummed
   ivirt=icasdet(imo_resum,kdet)

! Loop over basis functions
   evcoeff1(:,ifree_mo,ispin_rohf)=evcoeff1(:,ifree_mo,ispin_rohf) &
    &+ci1*evcoeff1(:,ivirt,ispin_rohf)
  enddo
! Replace the MO in the original determinant with the new, resummed one
  icasdet(imo_resum,idet_resum)=ifree_mo

 enddo ! End of loop to search for determinants to resum

END SUBROUTINE resum_cas

MODULE casino_interface
!-------------------------------------------------------------------------!
! Interface to the CASINO Quantum Monte Carlo code                        !
! See www.tcm.phy.cam.ac.uk/~mdt26/casino.html                            !
! 12.2002 Mike Towler, University of Cambridge                            !
!                                                                         !
! Changes                                                                 !
! -------                                                                 !
! 11.2003 MDT - Tidied and updated.                                       !
!-------------------------------------------------------------------------!
 USE basato_module
 USE expo_module
 USE gvect_module
 USE numbers
 USE parinf_module
 USE retic_module
 USE xyvdin_module

 IMPLICIT NONE
 PRIVATE
 PUBLIC qmc_main

 INTEGER ialloc,mvf,ndf,nk(3),number_of_supercells
 INTEGER,PARAMETER :: r2s_length=80
 CHARACTER(r2s_length) tmpr
 LOGICAL periodic


CONTAINS


 SUBROUTINE qmc_main
!-------------------------------------------------------------------------!
! Read the relevant part of the properties input file then loop over      !
! the different supercell sizes required. For each cell size, regenerate  !
! the k point net then write out the appropriate gwfn.data file which     !
! contains the trial wave function in a format understood by the CASINO   !
! program.                                                                !
!-------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,ierr
 INTEGER,ALLOCATABLE :: nc(:,:)
 REAL(FLOAT) r1,r2,r3
 REAL(FLOAT),PARAMETER :: tol=1.e-6_float

 ialloc=0

! Read the properties input file

 periodic=inf(10)/=0

 write(6,*)
 write(6,'(1x,a)')'CRYSTAL-CASINO INTERFACE FOR QUANTUM MONTE CARLO &
  &CALCULATIONS'
 write(6,'(1x,a)')'-------------------------------------------------&
  &------------'
 write(6,*)

 if(periodic)then

  read(iin,*,iostat=ierr)number_of_supercells
  if(ierr/=0)call errvrs(0,'CASINO_INTERFACE','Error reading properties input.')

  allocate(nc(3,number_of_supercells),stat=ialloc)
  if(ialloc/=0)call errvrs(0,'CASINO_INTERFACE','Allocation problem <1>.')

  do i=1,number_of_supercells
   read(iin,*,iostat=ierr)nc(1,i),nc(2,i),nc(3,i)
  enddo
  if(ierr/=0)call errvrs(0,'CASINO_INTERFACE','Error reading properties input.')

  nk(1)=inf(53) ; nk(2)=inf(54) ; nk(3)=inf(55) ! no. of MP k points

  do i=1,number_of_supercells
 ! Check QMC supercell size is same as or subset of original CRYSTAL MP k grid.
   r1=mod(dble(nk(1)),dble(nc(1,i)))
   r2=mod(dble(nk(2)),dble(nc(2,i)))
   r3=mod(dble(nk(3)),dble(nc(3,i)))
   if(r1>tol.or.r2>tol.or.r3>tol)call errvrs(0,'CASINO_INTERFACE',&
    &'Mismatch between CRYSTAL k grid and QMC supercell.')
   call qmc(nc(1,i)) ! Generate gwfn.data for this cell size
  enddo

 else ! atoms/molecules - assume no cell info in properties input.

  number_of_supercells=1
  allocate(nc(3,1),stat=ialloc)
  if(ialloc/=0)call errvrs(0,'CASINO_INTERFACE','Allocation problem <2>.')
  nc(1:3,1)=1 ; nk(1)=1 ; nk(2)=1 ; nk(3)=1
  call qmc(nc(1,1))

 endif
 deallocate(nc)

 END SUBROUTINE qmc_main


 SUBROUTINE qmc(nc)
!-------------------------------------------------------------------------!
! Write a gwfn.data file for the given supercell size nc                  !
!-------------------------------------------------------------------------!
 IMPLICIT NONE

 INTEGER,INTENT(in) :: nc(3)
 INTEGER idum,laf,naf,nprim,ndfrf,ndfcf,k,mvf,i,j,l1,                        &
  &l2,l3,lll,laz(48),lbz(48),lcz(48),lay(48),lby(48),lcy(48),mv,             &
  &is10,is20,is21,is30,k1,k2,k3,spin,iq1,iq2,iq3,nkfqmc,nkf_irred,           &
  &norb,n,jbase,kbase,nb,ind,ink,iop,jop,nelec,na,kr,kc,ik,ck_size,j3(1),    &
  &kmax,ntot,io,ierr,nkfqmc_irred,lric(48),symmop,mr1,mr2,mr3,nr1,nr2,nr3,   &
  &nkf_irred2,kcmax
 INTEGER,ALLOCATABLE :: latqmc(:),mqmc(:,:),latvrs_local(:),jk(:,:),         &
  &ak_pointer(:,:),kir(:),kic(:),kir_plucked(:),kic_plucked(:),latvrs_qmc(:),&
  &latvrs_qmc_irred(:),latqmc_irred(:),latvrs_full(:),kmap(:),found_k(:),    &
  &irr_local(:,:,:),jjj_local(:,:)
 REAL(FLOAT) brot(3,3),det,ris1,ris2,ris3,rtmp(3),root3,                     &
  &harmonic_coeff_times_norm(28,8),fraction_dftx
 REAL(FLOAT),ALLOCATABLE :: akt1(:),akt2(:),kvec(:),ktemp(:,:),              &
  &kvec_reordered(:),c_prim(:),ak(:),ckba(:)
! knetout output
 REAL(FLOAT),ALLOCATABLE :: wpj_local(:),eigenvalues(:)
 LOGICAL eigenvalues_present,lvet(48),op
 LOGICAL,ALLOCATABLE :: issp(:)
 CHARACTER(5) exch_functional(5),corr_functional(8)
 CHARACTER(80) gwfn_filename
 CHARACTER(13),PARAMETER :: reals3="(3(1pe20.13))"
 CHARACTER(13),PARAMETER :: reals4="(4(1pe20.13))"
 CHARACTER(13),PARAMETER :: int8="(8i10)"
 DATA exch_functional/'LDA','VBH','Becke','PBE','PWGGA'/
 DATA corr_functional/'PWLSD','PZ','VWN','VBH','P86','PWGGA','LYP','PBE'/

 if(number_of_supercells==1)then
  gwfn_filename='gwfn.data'
 else
  gwfn_filename='gwfn.data_'//trim(i2s(nc(1)))//trim(i2s(nc(2)))//&
   &trim(i2s(nc(3)))
 endif

! Basic definitions

 ndf=inf(7)    ! number of AOs / number of bands
 laf=inf(20)   ! number of shells
 naf=inf(24)   ! number of atoms
 nprim=inf(75) ! number of Gaussian primitives
 ndfrf=ndf*ndf
 ndfcf=ndfrf+ndfrf
 allocate(c_prim(nprim),stat=ialloc)
 if(ialloc/=0)call errvrs(0,'CASINO_INTERFACE','Allocation problem <3>.')
 brot=transpose(paret)/par(34)
 call minv3(brot,bret,det)

! Allocate arrays which depend on QMC supercell size

 kmax=nk(1)*nk(2)*nk(3)  ; kcmax=nc(1)*nc(2)*nc(3)
 allocate(latqmc(kmax),issp(kmax),mqmc(3,kcmax),kmap(kmax),                  &
  &kvec(3*kcmax),kvec_reordered(3*kcmax),kir(kcmax),kic(kcmax),ktemp(3,kmax),&
  &jk(kcmax,2),ak_pointer(kcmax,2),kir_plucked(kcmax),kic_plucked(kcmax),    &
  &latvrs_qmc(kcmax),latvrs_qmc_irred(kcmax),latqmc_irred(kmax),             &
  &latvrs_full(kmax),found_k(kmax),stat=ialloc)
 if(ialloc/=0)call errvrs(0,'CASINO_INTERFACE','Allocation problem <4>.')

 latqmc_irred(:)=0 ; latqmc(:)=0 ; latvrs_full(:)=0

! Basic properties of the irreducible QMC grid

 is1=nc(1)               ; is2=nc(2)               ; is3=nc(3)
 is10=is1*16             ; is20=is2*16             ; is30=is3*16
 iq1=nk(1)/is1           ; iq2=nk(2)/is2           ; iq3=nk(3)/is3
 is21=is1*is2
 mvf=inf(50)
 issp(:)=.false.
 lll=0
 nkfqmc_irred=0
 do l3=0,is3-1
  laz(1:mvf)=irr(1,3,1:mvf)*l3+is10
  lbz(1:mvf)=irr(2,3,1:mvf)*l3+is20
  lcz(1:mvf)=irr(3,3,1:mvf)*l3+is30
  do l2=0,is2-1
   lay(1:mvf)=irr(1,2,1:mvf)*l2+laz(1:mvf)
   lby(1:mvf)=irr(2,2,1:mvf)*l2+lbz(1:mvf)
   lcy(1:mvf)=irr(3,2,1:mvf)*l2+lcz(1:mvf)
   do l1=0,is1-1
    lll=lll+1
    if(issp(lll))cycle
    nkfqmc_irred=nkfqmc_irred+1                                ! OUTPUT
    do mv=1,mvf ! loop over symmops
     k3=mod(irr(3,1,mv)*l1+lcy(mv),is3)
     k2=mod(irr(2,1,mv)*l1+lby(mv),is2)
     k1=mod(irr(1,1,mv)*l1+lay(mv),is1)
     mqmc(1,nkfqmc_irred)=k1*iq1                               ! OUTPUT
     mqmc(2,nkfqmc_irred)=k2*iq2                               ! OUTPUT
     mqmc(3,nkfqmc_irred)=k3*iq3                               ! OUTPUT
     k=k3*is21+k2*is1+k1
     if(issp(k+1))cycle
     issp(k+1)=.true.
     if(k1/=0)k1=is1-k1
     if(k2/=0)k2=is2-k2
     if(k3/=0)k3=is3-k3
     k=k3*is21+k2*is1+k1
     if(issp(k+1))cycle
     issp(k+1)=.true.
    enddo ! mvf (symmops)
   enddo ! l1
  enddo ! l2
 enddo ! l3

 do k=1,nkfqmc_irred
  latvrs_qmc_irred(k)=0
  if(mod(mqmc(1,k)*2,nk(1))/=0)cycle
  if(mod(mqmc(2,k)*2,nk(2))/=0)cycle
  if(mod(mqmc(3,k)*2,nk(3))==0)latvrs_qmc_irred(k)=1
 enddo

! Basic properties of the irreducible MP grid

 is1=nk(1) ; is2=nk(2) ; is3=nk(3)
 is21=is1*is2
 is10=is1*16
 is20=is2*16
 is30=is3*16
 kvec(:)=0._float
 issp(:)=.false.
 nkf_irred=0
 lll=0
 do l3=0,is3-1
  laz(1:mvf)=irr(1,3,1:mvf)*l3+is10
  lbz(1:mvf)=irr(2,3,1:mvf)*l3+is20
  lcz(1:mvf)=irr(3,3,1:mvf)*l3+is30
  do l2=0,is2-1
   lay(1:mvf)=irr(1,2,1:mvf)*l2+laz(1:mvf)
   lby(1:mvf)=irr(2,2,1:mvf)*l2+lbz(1:mvf)
   lcy(1:mvf)=irr(3,2,1:mvf)*l2+lcz(1:mvf)
   do l1=0,is1-1
    lll=lll+1
    if(issp(lll))cycle
    nkf_irred=nkf_irred+1                                         ! OUTPUT
    do mv=1,mvf
     k3=mod(irr(3,1,mv)*l1+lcy(mv),is3)
     k2=mod(irr(2,1,mv)*l1+lby(mv),is2)
     k1=mod(irr(1,1,mv)*l1+lay(mv),is1)
     k=k3*is21+k2*is1+k1
     if(issp(k+1))cycle
     do i=1,nkfqmc_irred
      if(k1==mqmc(1,i).and.k2==mqmc(2,i).and.k3==mqmc(3,i))then
       latqmc_irred(nkf_irred)=1                                  ! OUTPUT
       exit
      endif
     enddo
     issp(k+1)=.true.
     if(k1/=0)k1=is1-k1
     if(k2/=0)k2=is2-k2
     if(k3/=0)k3=is3-k3
     k=k3*is21+k2*is1+k1
     if(issp(k+1))cycle
     do i=1,nkfqmc_irred
      if(k1==mqmc(1,i).and.k2==mqmc(2,i).and.k3==mqmc(3,i))then
       latqmc_irred(nkf_irred)=1                                  ! OUTPUT
       exit
      endif
     enddo
     issp(k+1)=.true.
    enddo ! mv (symmops)
   enddo ! l1
  enddo ! l2
 enddo ! l3

! Calculate number of real and complex k points in the full QMC grid, plus the
! coordinates of the corresponding k points.

 nkfqmc=0 ; kr=0 ; kc=0 ; ind=1
 ris1=1._float/(nc(1)*iq1)
 ris2=1._float/(nc(2)*iq2)
 ris3=1._float/(nc(3)*iq3)
a:do k=1,nkfqmc_irred
  iop=0
  l1=mqmc(1,k) ; l2=mqmc(2,k) ; l3=mqmc(3,k) ! integer coords of QMC irred k
! Establish invariant subgroup
  if(latvrs_qmc_irred(k)==1)then ! real k point
   lvet(1:mvf)=.true.
   do mv=1,mvf ! symmops
    laz(mv)=mod(irr(1,1,mv)*l1+irr(1,2,mv)*l2+irr(1,3,mv)*l3+is10,is1)
    lbz(mv)=mod(irr(2,1,mv)*l1+irr(2,2,mv)*l2+irr(2,3,mv)*l3+is20,is2)
    lcz(mv)=mod(irr(3,1,mv)*l1+irr(3,2,mv)*l2+irr(3,3,mv)*l3+is30,is3)
    if(abs(laz(mv)-l1)+abs(lbz(mv)-l2)+abs(lcz(mv)-l3)/=0)cycle
    iop=iop+1
    lric(iop)=mv
   enddo
  else ! complex k point
   nr1=mod(is1-l1,is1) ; nr2=mod(is2-l2,is2) ; nr3=mod(is3-l3,is3)
   lvet(1:mvf)=.true.
   do mv=1,mvf ! symmops
    mr1=mod(irr(1,1,mv)*l1+irr(1,2,mv)*l2+irr(1,3,mv)*l3+is10,is1)
    mr2=mod(irr(2,1,mv)*l1+irr(2,2,mv)*l2+irr(2,3,mv)*l3+is20,is2)
    mr3=mod(irr(3,1,mv)*l1+irr(3,2,mv)*l2+irr(3,3,mv)*l3+is30,is3)
    if(abs(mr1-l1)+abs(mr2-l2)+abs(mr3-l3)==0)then
     iop=iop+1
     lric(iop)=mv
    else
     if(abs(mr1-nr1)+abs(mr2-nr2)+abs(mr3-nr3)==0)then
      iop=iop+1
      lric(iop)=-mv
     endif
    endif
    laz(mv)=mr1 ; lbz(mv)=mr2 ; lcz(mv)=mr3
   enddo
  endif ! real or complex k point
1 continue
  do symmop=1,mvf
   if(lvet(symmop))then
    nkfqmc=nkfqmc+1                                               ! OUTPUT
    rtmp(1)=real(laz(symmop),kind=kind(0._float))*ris1
    rtmp(2)=real(lbz(symmop),kind=kind(0._float))*ris2
    rtmp(3)=real(lcz(symmop),kind=kind(0._float))*ris3
    call mxmb(bret,3,1,rtmp,1,1,kvec(ind),1,1,3,3,1)              ! OUTPUT
    ind=ind+3
    if(latvrs_qmc_irred(k)==1)then
     kr=kr+1                                                      ! OUTPUT
    else
     kc=kc+1                                                      ! OUTPUT
    endif
    do jop=1,iop
     if(latvrs_qmc_irred(k)==1)then
      mv=multab(symmop,lric(jop))
      lvet(mv)=.false.
     else
      i=lric(jop)
      j=abs(i)
      mv=multab(symmop,j)
      lvet(mv)=.false.
     endif
    enddo
    goto 1
   endif ! lvet(symmop)
  enddo ! symmops
 enddo a

! Read in IBZ orbital coefficients from fort.8, discarding those from
! k points we don't want, unfold the eigenvector coeffs into the full
! Brillouin zone, and store them in temporary vector AK.

 ck_size=(kr*ndfrf+(nkfqmc-kr)*ndfcf)*(inf(64)+1)

 allocate(ak(ck_size),ckba(ck_size),stat=ialloc)
 if(ialloc/=0)call errvrs(0,'CASINO_INTERFACE','Allocation problem <5>.')
 if(kc==0)then
  allocate(akt1(ndfrf),stat=ialloc)
 else
  allocate(akt1(ndfcf),stat=ialloc)
 endif
 if(ialloc/=0)call errvrs(0,'CASINO_INTERFACE','Allocation problem <6>.')
 if(periodic)then
  allocate(akt2(ndfcf),stat=ialloc)
  if(ialloc/=0)call errvrs(0,'CASINO_INTERFACE','Allocation problem <7>.')
 endif

 ind=0
 mvf=inf(2)

 do spin=0,inf(64)

  ik=0 ; nkf=0
b:do k=1,nkf_irred

   l1=jj(1,k) ; l2=jj(2,k) ; l3=jj(3,k) ! jj are integer coords of MP IBZ k

! Read eigenvectors off disk
   if(latvrs(k)==1)then ! real k point
    call rread(8,akt1,ndfrf)
   else                 ! complex k point
    call rread(8,akt1,ndfcf)
   endif

! Transform the orbital coefficients from IBZ k to new k
   if(periodic)then

    found_k(:)=0

    if(latvrs(k)==1)then ! REAL

     do mv=1,mvf

      mr1=mod(irr(1,1,mv)*l1+irr(1,2,mv)*l2+irr(1,3,mv)*l3+is10,is1)
      mr2=mod(irr(2,1,mv)*l1+irr(2,2,mv)*l2+irr(2,3,mv)*l3+is20,is2)
      mr3=mod(irr(3,1,mv)*l1+irr(3,2,mv)*l2+irr(3,3,mv)*l3+is30,is3)
      i=mr1+1+mr2*is1+mr3*is21
      if(found_k(i)/=0)cycle
      if(latqmc_irred(k)==1)call expt(mr1,mr2,mr3) ! calc e(ik.g) for lots of g
      found_k(i)=1

      nkf=nkf+1                                                 ! OUTPUT
      kmap(nkf)=k                                               ! OUTPUT

      latvrs_full(nkf)=1                                        ! OUTPUT

      if(latqmc_irred(k)==1)then ! CASINO wants this k point
       latqmc(nkf)=1                                            ! OUTPUT
       call estrog(akt1,akt2,ninv(mv),1,ndf,ndf)
       do i=1,ndfrf
        if(abs(akt2(i))>1.e-12_float)then
         ak(ind+i)=akt2(i)                                      ! OUTPUT
        else
         ak(ind+i)=0._float
        endif
       enddo
       ind=ind+ndfrf
      endif

     enddo

    else ! COMPLEX

     do mv=1,mvf

      mr1=mod(irr(1,1,mv)*l1+irr(1,2,mv)*l2+irr(1,3,mv)*l3+is10,is1)
      mr2=mod(irr(2,1,mv)*l1+irr(2,2,mv)*l2+irr(2,3,mv)*l3+is20,is2)
      mr3=mod(irr(3,1,mv)*l1+irr(3,2,mv)*l2+irr(3,3,mv)*l3+is30,is3)
      i=mr1+1+mr2*is1+mr3*is21
      if(found_k(i)/=0)cycle
      mr1=mod(is1-mr1,is1) ; mr2=mod(is2-mr2,is2) ; mr3=mod(is3-mr3,is3)
      i=mr1+1+mr2*is1+mr3*is21
      if(found_k(i)/=0)cycle
      if(latqmc_irred(k)==1)call expu(mr1,mr2,mr3)
      found_k(i)=2

      nkf=nkf+1                                                 ! OUTPUT
      kmap(nkf)=k                                               ! OUTPUT

      if(latqmc_irred(k)==1)then ! CASINO wants this k point
       latqmc(nkf)=1                                            ! OUTPUT
       call estrof(akt1,akt2,ninv(mv),1,ndf,ndf)
       do i=1,ndfcf
        if(abs(akt2(i))>1.e-12_float)then
         ak(ind+i)=akt2(i)                                      ! OUTPUT
        else
         ak(ind+i)=0._float
        endif
       enddo
       ind=ind+ndfcf
      endif ! CASINO wants this k point

     enddo

     if(mvf==inf(1))cycle b
     l1=mod(is1-l1,is1) ; l2=mod(is2-l2,is2) ; l3=mod(is3-l3,is3)

     do mv=1,mvf

      mr1=mod(irr(1,1,mv)*l1+irr(1,2,mv)*l2+irr(1,3,mv)*l3+is10,is1)
      mr2=mod(irr(2,1,mv)*l1+irr(2,2,mv)*l2+irr(2,3,mv)*l3+is20,is2)
      mr3=mod(irr(3,1,mv)*l1+irr(3,2,mv)*l2+irr(3,3,mv)*l3+is30,is3)
      i=mr1+1+mr2*is1+mr3*is21
      if(found_k(i)/=0)cycle
      if(latqmc_irred(k)==1)call expu(mr1,mr2,mr3)
      found_k(i)=2

      nkf=nkf+1                                                 ! OUTPUT
      kmap(nkf)=k                                               ! OUTPUT

      if(latqmc_irred(k)==1)then ! CASINO wants this k point
       latqmc(nkf)=1                                            ! OUTPUT
       call estroe(akt1,akt2,ninv(mv),1,ndf,ndf)
       do i=1,ndfcf
        if(abs(akt2(i))>1.e-12_float)then
         ak(ind+i)=akt2(i)                                     ! OUTPUT
        else
         ak(ind+i)=0._float
        endif
       enddo
       ind=ind+ndfcf
      endif ! CASINO wants this k point

     enddo

    endif ! real or complex k

   else ! atom/molecule

    nkf=1
    kmap(1)=1
    latqmc(1)=1
    latvrs_full(1)=1
    n=ndfrf*spin
    do i=1,ndfrf
     if(abs(akt1(i))>1.e-12_float)then
      ak(n+i)=akt1(i)
     else
      ak(n+i)=0._float
     endif
    enddo

   endif ! periodic or not

  enddo b ! irreducible k

 enddo ! spin

 deallocate(akt1)
 if(periodic)deallocate(akt2)
 rewind(8)

! Mess about with k points

 if(periodic)then

  ik=1 ; ink=0 ; kr=0 ; kc=0
  do k=1,nkf
   if(latqmc(k)==1)then
    ink=ink+1
    ktemp(1,ink)=kvec(ik)
    ktemp(2,ink)=kvec(ik+1)
    ktemp(3,ink)=kvec(ik+2)
    if(latvrs_full(k)==1)then
     kr=kr+1
     kir(kr)=k ! not used
     kir_plucked(kr)=ink
    else
     kc=kc+1
     kic(kc)=k ! not used
     kic_plucked(kc)=ink
    endif
    ik=ik+3
   endif
  enddo
  ik=1
  do i=1,kr
   k=kir_plucked(i)
   kvec_reordered(ik)=ktemp(1,k)
   kvec_reordered(ik+1)=ktemp(2,k)
   kvec_reordered(ik+2)=ktemp(3,k)
   ik=ik+3
  enddo
  do i=1,kc
   k=kic_plucked(i)
   kvec_reordered(ik)=ktemp(1,k)
   kvec_reordered(ik+1)=ktemp(2,k)
   kvec_reordered(ik+2)=ktemp(3,k)
   ik=ik+3
  enddo

 else

  kr=1 ; kc=0 ; kir_plucked(1)=1
  kir(1)=1 ! not used

 endif

! Reorder orbital coefficients in ak --> ckba
! Also multiply coefficients of d basis functions by the
! harmonic coefficients * a missing normalization piece
! ao        harmonic coefficient      normalization
! 3z2-r2         0.5                    1
! xz             3.0                    1/root3
! yz             3.0                    1/root3
! x2-y2          3.0                    1/(2*root3)
! xy             6.0                    1/(2*root3)

 n=0
 do spin=1,inf(64)+1
  ink=0
  do k=1,nkf
   if(latqmc(k)==1)then
    ink=ink+1
    ak_pointer(ink,spin)=n
    if(latvrs_full(k)==1)then
     n=n+ndfrf
    else
     n=n+ndfcf
    endif
   endif
  enddo
 enddo
 do spin=1,inf(64)+1
  do k=1,kr
   jk(k,spin)=ak_pointer(kir_plucked(k),spin)
  enddo
  do k=1,kc
   jk(kr+k,spin)=ak_pointer(kic_plucked(k),spin)
  enddo
 enddo

 deallocate(ak_pointer)

! s
 harmonic_coeff_times_norm(1,1)=1._float
!sp
 harmonic_coeff_times_norm(1:4,2)=1._float
! p
 harmonic_coeff_times_norm(1:3,3)=1._float
!d (harmonic)
 root3=sqrt(3._float)
 harmonic_coeff_times_norm(1,4)=0.5_float
 harmonic_coeff_times_norm(2,4)=root3
 harmonic_coeff_times_norm(3,4)=root3
 harmonic_coeff_times_norm(4,4)=0.5_float*root3
 harmonic_coeff_times_norm(5,4)=root3

 j=1
 do spin=1,inf(64)+1
  do k=1,kr ! real k points in full QMC grid
   i=1
   jbase=jk(k,spin)
   do nb=1,ndf
    do l1=1,laf
     norb=latao(l1)
     kbase=lat(l1)+1
     do n=1,norb
      ckba(j)=ak(jbase+i)*harmonic_coeff_times_norm(n,kbase)
      i=i+1
      j=j+1
     enddo ! ao in shell
    enddo ! shells
   enddo ! bands
  enddo ! real k
  do k=kr+1,nkfqmc ! complex k points in full QMC grid
   i=1
   jbase=jk(k,spin)
   do nb=1,ndf
    do l1=1,laf
     norb=latao(l1)
     kbase=lat(l1)+1
     do n=1,norb
      ckba(j)=ak(jbase+i)*harmonic_coeff_times_norm(n,kbase)
      ckba(j+1)=ak(jbase+i+1)*harmonic_coeff_times_norm(n,kbase)
      i=i+2
      j=j+2
     enddo ! ao in shell
    enddo ! shells
   enddo ! bands
  enddo ! real k
 enddo ! spins

 deallocate(ak,jk)

! Now read the eigenvalues from the fort.30 file (periodic systems only)

 if(periodic)then
  inquire(file='fort.30',exist=eigenvalues_present)
  if(eigenvalues_present)then
   read(30,iostat=ierr)idum,nkf_irred2,idum,bret
   if(ierr/=0)call errvrs(0,'CASINO_INTERFACE','Error reading fort.30 file.')
   if(nkf_irred2/=nkf_irred)then
    write(6,*)'Calculated NKF  :',nkf_irred
    write(6,*)'NKF from fort.30: ',nkf_irred2
    call errvrs(0,'CASINO_INTERFACE','Calculated NKF value different&
     & from that given in the fort.30 file.')
   endif
   ntot=(inf(64)+1)*ndf*nkf_irred
   allocate(jjj_local(3,nkf_irred),latvrs_local(nkf_irred), &
    &wpj_local(nkf_irred),eigenvalues(ntot),irr_local(3,3,48),stat=ialloc)
   if(ialloc/=0)call errvrs(0,'CASINO_INTERFACE','Allocation error <6>')
   call iread(30,jjj_local,3*nkf_irred)
   call iread(30,latvrs_local,nkf_irred)
   call iread(30,irr_local,9*48)
   call rread(30,wpj_local,nkf_irred)
   call rread(30,eigenvalues,ntot)
   deallocate(jjj_local,latvrs_local,wpj_local,irr_local)
  endif
  rewind(30)
 endif

!===========================================================

! Now write out the gwfn.data file for CASINO to read

 do io=60,100
  inquire(unit=io,opened=op)
  if(.not.op)exit
 enddo
 if(io>99)call errvrs(0,'CASINO_INTERFACE','Problem finding unused unit number')

 open(io,file=trim(adjustl(gwfn_filename)),status='unknown',form='formatted',&
  &iostat=ierr)
 if(ierr/=0)call errvrs(0,'CASINO_INTERFACE','Error opening gwfn.data file.')

!=================Basic information=========================

 write(io,fmt='(18a4)')(itit(i),i=1,18)
 write(io,fmt="(/'BASIC INFO'/'----------')")
 write(io,fmt="('Generated by:'/' CRYSTAL2003')")

 write(io,fmt="('Method:')")
 if(inf(170)/=0)then
  write(io,fmt="(' DFT')")
 else
  if(inf(27)==0.and.inf(64)==0)then
   write(io,fmt="(' RHF')")
  elseif(inf(27)==1.and.inf(64)==1)then
   write(io,fmt="(' UHF')")
  elseif(inf(27)==0.and.inf(64)==1)then
   write(io,fmt="(' ROHF')")
  endif
 endif

 write(io,fmt="('DFT functional:')")
 if(inf(170)/=0)then ! DFT
  if(inf(172)/=0)then
   fraction_dftx=1._float-par(48)*0.01_float
   if(fraction_dftx < 1._float)then ! hybrid xc functional
    write(io,fmt="('1x,Fock exchange * ',f7.5,' + ',a,' exchange * ',f7.5,&
     &' + ',a,' correlation')")1._float-fraction_dftx,                    &
     &trim(adjustl(exch_functional(inf(172)))),fraction_dftx,             &
     &trim(adjustl(corr_functional(inf(171))))
   else ! non-hybrid xc functional
    write(io,fmt="(1x,a,' correlation, ',a,' exchange')")                 &
     &trim(adjustl(corr_functional(inf(171)))),                           &
     &trim(adjustl(exch_functional(inf(172))))
   endif
  else ! correlation-only functional
   write(io,fmt="(' Hartree-Fock exchange plus ',a,' correlation.')")     &
    &trim(adjustl(corr_functional(inf(171))))
  endif
 else ! Hartree-Fock
  write(io,fmt="(' None')")
 endif

 write(io,fmt="('Periodicity:'/1x,a)")trim(i2s(inf(10)))
 write(io,fmt="('Spin unrestricted:')")
 if(inf(64)/=0)then
  write(io,fmt="(' .true.')")
 else
  write(io,fmt="(' .false.')")
 endif

 write(io,fmt="('Nuclear repulsion energy (au/atom):')")
 tmpr=r2s(par(18)/real(naf,kind=kind(0._float)),'(f24.15)')
 write(io,'(1x,a)')trim(tmpr)

 nelec=0
 do na=1,naf
  do l1=nshpri(na),nshpri(na+1)-1
   nelec=nelec+int(che(l1))
  enddo
 enddo

 if(periodic)then
  write(io,fmt="('Number of electrons per primitive cell')")
 else
  write(io,fmt="('Number of electrons')")
 endif
 write(io,'(1x,a)')trim(i2s(nelec))

!=================Geometry input==============================

 write(io,fmt="(/'GEOMETRY'/'--------')")

 write(io,fmt="('Number of atoms')")
 write(io,'(1x,a)')trim(i2s(naf))

 write(io,fmt="('Atomic positions (au)')")
 write(io,fmt=reals3)((xa(i,j),i=1,3),j=1,naf)

 write(io,fmt="('Atomic numbers for each atom')")
 write(io,fmt=int8)(nat(i),i=1,naf)

 write(io,fmt="('Valence charges for each atom')")
 write(io,fmt=reals4)(aznuc(i),i=1,naf)

 if(periodic)then

  write(io,fmt="('Primitive lattice vectors (au)')")
  write(io,fmt=reals3)((paret(i,j),j=1,3),i=1,3) ! i.e. don't write(")paret

!=================K SPACE NET=================================

  write(io,fmt="(/'K SPACE NET'/'-----------')")

  write(io,fmt="('Number of k points')")
  write(io,'(1x,a)')trim(i2s(nkfqmc))

  write(io,fmt="('Number of ''real'' k points on BZ edge')")
  write(io,'(1x,a)')trim(i2s(kr))

  write(io,fmt="('k point coordinates (au)')")
  write(io,fmt=reals3)(kvec_reordered(k),k=1,nkfqmc*3)

 endif ! polymers/slabs/solids

!=================Basis set===================================

 write(io,fmt="(/'BASIS SET'/'---------')")

 write(io,fmt="('Number of Gaussian centres')")
 write(io,'(1x,a)')trim(i2s(naf))

 if(periodic)then
  write(io,fmt="('Number of shells per primitive cell')")
 else
  write(io,fmt="('Number of shells')")
 endif
 write(io,'(1x,a)')trim(i2s(laf))

 if(periodic)then
  write(io,fmt="('Number of basis functions (''AO'') per primitive cell')")
 else
  write(io,fmt="('Number of basis functions (''AO'')')")
 endif
 write(io,'(1x,a)')trim(i2s(ndf))

 if(periodic)then
  write(io,fmt="('Number of Gaussian primitives per primitive cell')")
 else
  write(io,fmt="('Number of Gaussian primitives')")
 endif
 write(io,'(1x,a)')trim(i2s(nprim))

 j=maxval(lat(1:laf))

 write(io,fmt="('Highest shell angular momentum (s/p/d/f/g... 1/2/3/4/5...)')")
 if(j==0)then
  write(io,'(a)')' 1'
 elseif(j==1.or.j==2)then
  write(io,'(a)')' 2'
 elseif(j>=3)then
  write(io,'(1x,a)')trim(i2s(j))
 endif

 write(io,fmt="('Code for shell types (s/sp/p/d/f... 1/2/3/4/5...)')")
 write(io,fmt=int8)(lat(i)+1,i=1,laf)

 write(io,fmt="('Number of primitive Gaussians in each shell')")
 write(io,fmt=int8)(lan(i),i=1,laf)

 write(io,fmt="('Sequence number of first shell on each centre')")
 write(io,fmt=int8)(nshpri(i),i=1,naf+1)

 write(io,fmt="('Exponents of Gaussian primitives')")
 write(io,fmt=reals4)(exx(i),i=1,nprim)

 write(io,fmt="('Correctly normalised contraction coefficients')")
 do i=1,laf
  if(lat(i)==0.or.lat(i)==1)then
   c_prim(laa(i):laa(i+1)-1)=c1(laa(i):laa(i+1)-1)
  elseif(lat(i)==2)then
   c_prim(laa(i):laa(i+1)-1)=c2(laa(i):laa(i+1)-1)
  elseif(lat(i)==3)then
   c_prim(laa(i):laa(i+1)-1)=c3(laa(i):laa(i+1)-1)
  else
   call errvrs(0,'CASINO_INTERFACE','Don''t know where to find contraction &
    &coefficients for f and g functions.')
  endif
 enddo
 write(io,fmt=reals4)(c_prim(i),i=1,nprim)
 deallocate(c_prim)

 if(any(lat(1:laf)==1))then
  write(io,fmt="('2nd contraction coefficients (p coeff.,&
   & for sp shells, 0.0 otherwise)')")
  write(io,fmt=reals4)(c2(i),i=1,nprim)
 endif

 write(io,fmt="('Position of each shell (au)')")
 write(io,fmt=reals3)((xl(i,j),i=1,3),j=1,laf)

!=================Multideterminant information============================

 write(io,fmt="(/'MULTIDETERMINANT INFORMATION'/&
  &'----------------------------')")
 write(io,fmt="('GS')") ! only single determinant possible with CRYSTAL

!=================Orbital coefficients====================================

 write(io,fmt="(/'ORBITAL COEFFICIENTS'/'---------------------------')")

 write(io,fmt=reals4)(ckba(i),i=1,ck_size)

 deallocate(ckba)

!====================Write out band energies for periodic systems=========

 if(periodic.and.eigenvalues_present)then
  write(io,fmt="(/'EIGENVALUES'/'-----------')")
! First write out the band energies at all the 'real' k points
  ik=0 ; k=0 ; kr=0 ; kc=0
  n=ndf*(inf(64)+1)
  do j=1,nkf
   if(latqmc(j)==1.and.latvrs_full(j)==1)then
    k=k+1 ; kr=kr+1
    jbase=(kmap(kir_plucked(kr))-1)*n
    do spin=1,inf(64)+1
     if(inf(64)==0)then
      write(io,'(a,i6,3f14.8)')'k ',k,kvec_reordered(ik+1),&
      &kvec_reordered(ik+2),kvec_reordered(ik+3)
     else
      write(io,'(a,i2,a,i6,3f14.8)')'spin ',spin,' k ',k,&
      &kvec_reordered(ik+1),kvec_reordered(ik+2),kvec_reordered(ik+3)
     endif
     if(spin==1)then
      write(io,fmt=reals4)(eigenvalues(i),i=jbase+1,jbase+ndf)
     else
      write(io,fmt=reals4)(eigenvalues(i),i=jbase+ndf+1,jbase+ndf+ndf)
     endif
    enddo
    ik=ik+3
   endif
  enddo ! k
! Now write out the band energies for the 'complex' k points
  do j=1,nkf
   if(latqmc(j)==1.and.latvrs_full(j)==0)then ! 'complex' k points
    k=k+1 ; kc=kc+1
    jbase=(kmap(kic_plucked(kc))-1)*n
    do spin=1,inf(64)+1
     if(inf(64)==0)then
      write(io,'(a,i6,3f14.8)')'k ',k,kvec_reordered(ik+1),&
      &kvec_reordered(ik+2),kvec_reordered(ik+3)
     else
      write(io,'(a,i2,a,i6,3f14.8)')'spin ',spin,' k ',k,&
      &kvec_reordered(ik+1),kvec_reordered(ik+2),kvec_reordered(ik+3)
     endif
     if(spin==1)then
      write(io,fmt=reals4)(eigenvalues(i),i=jbase+1,jbase+ndf)
     else
      write(io,fmt=reals4)(eigenvalues(i),i=jbase+ndf+1,jbase+ndf+ndf)
     endif
    enddo
    ik=ik+3
   endif
  enddo ! k
 endif

 deallocate(kir_plucked,kic_plucked)
 if(periodic)deallocate(eigenvalues)

!====================Append any other info. to the input file===========

 write(io,fmt=&
  & "(/'Input and output files for this calculation (not read by CASINO)')")
 write(io,fmt=&
  & "('=====================================================',&
  & '===========================')")

 write(io,*)
 close(io)

 if(periodic)then
  write(6,*)'CASINO gwfn.data FILE GENERATED FOR ',trim(i2s(nc(1))),'x',&
   &trim(i2s(nc(2))),'x',trim(i2s(nc(3))),' PERIODIC SYSTEM.'
 else
  write(6,*)'CASINO gwfn.data FILE GENERATED FOR FINITE SYSTEM.'
 endif
 write(6,*)

 deallocate(latqmc,issp,mqmc,kmap,kvec,kvec_reordered,kir,kic,ktemp, &
  &latvrs_qmc,latvrs_qmc_irred,latqmc_irred,latvrs_full,found_k,stat=ialloc)


 END SUBROUTINE qmc


 CHARACTER(20) FUNCTION i2s(n)
!-----------------------------------------------------------------------!
! I2S                                                                   !
! ===                                                                   !
! Convert integers to left justified strings that can be printed in the !
! middle of a sentence without introducing large amounts of white space.!
!                                                                       !
! Calling routine is intended to include something like:                !
! INTEGER i                                                             !
! i=12                                                                  !
! write(6,*)'Integer number ',trim(i2s(i)),' with words at the end.'    !
!-----------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,n
 CHARACTER tmp,asign

 if(n==0)then
  i2s='0' ; return
 endif
 asign=' ' ; if(n<0)asign='-'

 do i=1,len(i2s)
  i2s(i:i)=' '
 enddo

 i=abs(n)
 do j=1,len(i2s)
  if(i==0)exit
  i2s(j:j)=achar(ichar('0')+mod(i,10))
  i=i/10
 enddo

 i=1 ; j=len_trim(i2s)
 do
  if(i>=j)exit
  tmp=i2s(j:j)
  i2s(j:j)=i2s(i:i)
  i2s(i:i)=tmp
  i=i+1
  j=j-1
 enddo

 i2s=trim(asign)//i2s

 END FUNCTION i2s


 CHARACTER(r2s_length) FUNCTION r2s(r,real_format)
!-------------------------------------------------------------------------!
! Converts real variable with arbitrary format to string that can be      !
! trimmed and printed in the middle of a sentence without introducing     !
! large amounts of white space, as you would if you did                   !
! write(6,'(f12.6)')12.0 or similar. Note you need to pass through the    !
! format string e.g. f12.6 .                                              !
!                                                                         !
! Calling routine is intended to include something like:                  !
! REAL(float) r                                                           !
! r=12.d0                                                                 !
! tmpr=r2s(r,'(f12.6)')                                                   !
! write(6,*)'Real number ',trim(tmpr),' with words at the end.'           !
!                                                                         !
! Note : DON'T USE R2S IN A WRITE STATEMENT SINCE THIS IS ILLEGAL         !
! IN FORTRAN90 (ALTHOUGH NOT IN FORTRAN200X). IF ANYONE HAS TIME, FEEL    !
! FREE TO WRITE A VERSION OF THIS WHICH ISN'T ILLEGAL - SIMILAR TO        !
! I2S ABOVE - SO THAT PEOPLE WHO HAVEN'T READ THIS NOTE DON'T FEEL        !
! TEMPTED TO CALL R2S IN A WRITE STATEMENT.                               !
!-------------------------------------------------------------------------!

 IMPLICIT NONE
 REAL(float),INTENT(in) :: r
 CHARACTER(*),INTENT(in) :: real_format

 if(len(r2s)>0)then
  write(r2s,real_format)r
  r2s=adjustl(r2s)
 endif

 END FUNCTION r2s


END MODULE casino_interface

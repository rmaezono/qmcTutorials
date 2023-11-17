SUBROUTINE read_fchk
!------------------------------------------------------------!
! Subroutine to read a formatted checkpoint file produced by !
! Gaussian94/98/03/09 and extract all information necessary  !
! to reconstruct the wave function.                          !
!------------------------------------------------------------!
 USE awk_like
 USE cis_data
 USE g94_wavefunction
 IMPLICIT NONE

 INTEGER icheck,i,ia,ic,nline
 INTEGER memtest ! Flag for checking success of memory allocation
 INTEGER Ncoeff2 ! Ncoeffs^2
 INTEGER Natom_x_3 ! Number of atoms * 3
 REAL(KIND=dp),ALLOCATABLE :: Temp(:) ! Temporary storage for MO coeff read
 REAL(KIND=dp) dummy ! Dummy variable for reading values from file
 LOGICAL search ! Flag to control scanning of the Fchk file
 LOGICAL lfound
 CHARACTER(35) fchk_file
 CHARACTER(65) instring ! Character string for taking lines of the Fchk file

 fchk_file=trim(g94_file)//'.fchk'
 open(unit=30,file=fchk_file,status='old',iostat=icheck)
 if(icheck>0)then 
  fchk_file=trim(g94_file)//'.Fchk'
  open(unit=30,file=fchk_file,status='old',iostat=icheck)
 endif
 if(icheck>0)then
  write(*,fmt="(/'Cannot open the formatted checkpoint file (.fchk or .Fchk)' )")
  stop
 else
  write(*,fmt="(/'Reading file ',a35)")fchk_file
 endif

 read(30,fmt="(a71)")title_txt
 read(30,fmt="(a71)")job_txt

! Identify whether this is spin-restricted or unrestricted calculation.
 call scan_string(job_txt,fields,NF,NF_max," ")

 if((fields(2)(1:1)=='R').or.(trim(fields(2))=="CASSCF"))then
  SPIN=.false.
  ispin_lim=1 ! Reset limit for loops over spin so as to exclude spin down
 endif

! Get the number of atoms
 read(30,fmt="(58x,i3)")Natom
 write(*,fmt="(/i3,' atoms.')")Natom

! Look for the specification of the number of electrons
 search=.true.
 do while(search)
  read(30,fmt="(a65)")instring
  call scan_string(instring,fields,NF,NF_max," ")
  if(fields(3)=="electrons")search=.false.
 enddo
 read(fields(NF),*)Nelec

! The number of alpha and beta electrons - use to determine which MOs
! are occupied in the ground state.
 call getval(30,'alpha',3,6,dummy) ; Nspin(1)=int(dummy)
 call getval(30,'beta',3,6,dummy) ; Nspin(2)=int(dummy)
 write(*,fmt="(i3,' alpha and ',i3,' beta electrons.')")Nspin

! The number of basis functions
 call getval(30,'basis',3,6,dummy) ; Ncoeffs=int(dummy)
 Ncoeff2=Ncoeffs*Ncoeffs
 write(*,fmt="(i4,' basis functions.')")Ncoeffs
  call checkval(30, 'independant',3, lfound)
  if (lfound) then
    call getval(30,'independant',3,6,dummy) ; Nmo=int(dummy)
    write(*,fmt="(i4,' independent functions.')")Nmo
    Ncoeff2=Nmo*Ncoeffs
  else
    call checkval(30, 'independent',3, lfound)
    if (lfound) then 
      call getval(30,'independent',3,6,dummy) ; Nmo=int(dummy)
      write(*,fmt="(i4,' independent functions.')")Nmo
      Ncoeff2=Nmo*Ncoeffs
    else 
    write(*,*) "none"
    Nmo=Ncoeffs
    Ncoeff2=Ncoeffs*Ncoeffs
    endif
  endif

! Nmo_total is the number of MOs we allow room for in evcoeff1.  This
! must be large enough to provide temporary storage for resummed orbitals.
! For CIS, Nmo_total need only be equal to Nmo + the number of electrons
! but for CASSCF it must be larger.  The factor of 5 below is empirical.
 Nmo_total=Nmo+5*Nelec

!===============Memory allocation (1)=================================

! The number of shells (may or may not be contracted)
  call getval(30,'contracted',3,6,dummy) ; Nshells=int(dummy)
  write(*,fmt="(i4,' contracted shells.')")Nshells
! Highest degree of angular momentum in the basis
  call getval(30,'angular',2,5,dummy) ; Max_AM=int(dummy)
  write(*,fmt="('Max. AM of any shell is ',i2,'.')")Max_AM
! Largest degree of contraction
  call getval(30,'degree',2,6,dummy) ; Max_contract=int(dummy)
  write(*,fmt="('Largest  degree of contraction is ',i2,'.')")Max_contract
! The number of primitive Gaussian type functions
  call getval(30,'primitive',3,6,dummy) ; Nprimgtf=int(dummy)
  write(*,fmt="(i4,' primitive shells.')")Nprimgtf
! Allocations (done later for Gaussian 09)
  allocate(Lshell(Nshells),shexpnt(Nprimgtf),shll_posn(3,Nshells), &
   &ish_map(Nshells),evcoeff1(Ncoeffs,Nmo_total,2), &
   &con_coeff(Max_contract,Nshells),Nprim(Nshells), &
   &sp_coeff(Max_contract,Nshells), &
   &c_prim(Nprimgtf),c2_prim(Nprimgtf),list_of_free_mos(Nmo),stat=memtest)
  if(memtest/=0)stop 'Could not allocate memory in read_fchk (1a).'
  Lshell=0 ; shexpnt=0.d0 ; shll_posn=0.d0 ; ish_map=0
  evcoeff1=0.d0 ; con_coeff=0.d0
  Nprim=0 ; sp_coeff=0.d0 ; c_prim=0.d0 ; c2_prim=0.d0

 allocate(Temp(2*Ncoeff2),atom_posn(3,Natom),Natomic_no(Natom), &
  &Atomic_chrg(Natom),stat=memtest)

 if(memtest/=0)stop 'Could not allocate memory in read_fchk (2).'
 Temp=0.d0 ; atom_posn=0.d0 ; Natomic_no=0 ; Atomic_chrg=0.d0

! ==================Atomic numbers====================================

 call getval(30,'Atomic',1,5,dummy)
 if(Natom/=int(dummy))then
  write(6,fmt="('ERROR - number of atoms does not agree with number of&
   & atomic numbers specified in Fchk file.')")
  stop
 endif
 write(*,fmt="('Reading atomic numbers...')")

 nline=max((Natom/6),1) ! If only have 1 data point still have 1 line to read
 ic=0
 do i=1,nline
  read(30,*)Natomic_no(ic+1:min(ic+6,Natom))
  ic=ic+6
 enddo
! However, there may be a last line with less than 6 numbers..
 if((mod(Natom,6)/=0).and.(Natom>6))then
  read(30,*)Natomic_no(ic+1:Natom)
 endif

! ==================Nuclear charges================================

 write(*,fmt="('Reading valence charges...')")

 call getval2(30,'Nuclear',1,'charges',2,5,dummy)
 if(Natom/=int(dummy))then
  write(6,fmt="('ERROR - number of atoms does not agree with number of&
   & nuclear charges specified in Fchk file.')")
  stop
 endif

 nline=max((Natom/5),1) ! If only have 1 data point still have 1 line to read
 ic=0
 do i=1,nline
  read(30,*)Atomic_chrg(ic+1:min(ic+5,Natom))
  ic=ic+5
 enddo
! However, there may be a last line with less than 5 numbers...
 if((mod(Natom,5)/=0).and.(Natom>5))then
  read(30,*)Atomic_chrg(ic+1:Natom)
 endif

! ==================Cartesian coordinates of each atomic centre=====

 write(*,fmt="('Reading coordinates of each centre...')")

 call getval(30,'cartesian',2,6,dummy)
 if(3*Natom/=int(dummy))then
  write(6,fmt="('ERROR - 3x number of atoms does not agree with number of&
   & cartesian coordinates specified in Fchk file.')")
  stop
 endif
! Store x,y and z sequentially in temp and then transfer to directly
! referenceable form in atom_posn(x_i,iatom)

 Natom_x_3=3*Natom
 nline=max(Natom_x_3/5,1)

 ic=0
 do i=1,nline
  read(30,*)Temp(ic+1:min(Natom_x_3,ic+5))
  ic=ic+5
 enddo

! There may be a last line with less than 5 numbers..
 if((mod(Natom_x_3,5)/=0).and.(Natom_x_3>5))then
  read(30,*)temp(ic+1:Natom_x_3)
 endif

! Store in directly referenceable form..
 ic=0
 do ia=1,Natom
  atom_posn(1:3,ia)=temp(ic+1:ic+3)
  ic=ic+3
 enddo

! =================AM of each shell/ Shell type========================

 write(*,fmt="('Reading AM of each shell...')")

 call getval2(30,'Shell',1,'types',2,5,dummy)
 if(Nshells/=int(dummy))then
  write(6,fmt="('ERROR - no. of contracted shells does not agree with number of&
   & shell types specified in Fchk file.')")
  stop
 endif

! There are six fields on each line and Nshell numbers in total
! so the number of complete lines is:
 nline=max(Nshells/6,1)

 ic=0
 do i=1,nline
  read(30,*)Lshell(ic+1:min(ic+6,Nshells))
  ic=ic+6
 enddo
! However, there may be a last line with less than 6 numbers...
 if((mod(Nshells,6)/=0).and.(Nshells>6))then
  read(30,*)Lshell(ic+1:Nshells)
 endif

! L=-1 indicates an SP shell (whereas L=1 is a P shell) and L=+2
! indicates the use of cartesian type D shells rather than the 'pure'
! type with only five basis functions.
 num_sp=0
 do i=1,Nshells
  if(Lshell(i)==-1)then
   num_sp=num_sp+1
  elseif(Lshell(i)>3)then
   write(*,fmt="('h and higher angular momentum shells are not supported.')")
   stop
  endif
 enddo
 write(*,fmt="('We have ',i2,' SP shells.')")num_sp

! ==============The number of primitives/shell===============================

 write(*,fmt="('Reading the no. of primitives/shell...')")

 call getval(30,'primitives',3,8,dummy)
 if(Nshells/=int(dummy))then
  write(6,fmt="('ERROR - number of contracted shells does not agree with number&
   & of shells for which number of primitives is specified in Fchk file.')")
  stop
 endif

! nline is still Nshells/6 so no need to calculate it again...
 ic=0
 do i=1,nline
  read(30,*)Nprim(ic+1:min(ic+6,Nshells))
  ic=ic+6
 enddo
! However, there may be a last line with less than 6 numbers...
 if((mod(Nshells,6)/=0).and.(Nshells>6))then
  read(30,*)Nprim(ic+1:Nshells)
 endif

! ======================Shell to atom map======================================

 write(*,fmt="('Reading the shell to atom map...')")

 call getval(30,'map',4,7,dummy)
 if(Nshells/=int(dummy))then
  write(6,fmt="('ERROR - number of contracted shells does not agree with number&
   & of shells for which atom mapping is specified in Fchk file.')")
  stop
 endif

! Same no. of data points and hence lines as in the last section
 ic=0
 do i=1,nline
  read(30,*)ish_map(ic+1:min(ic+6,Nshells))
  ic=ic+6
 enddo
! However, there may be a last line with less than 6 numbers...
 if((mod(Nshells,6)/=0).and.(Nshells>6))then
  read(30,*)ish_map(ic+1:Nshells)
 endif

!=======================Primitive exponents=============================

 write(*,fmt="('Reading the primitive exponents...')")

 call getval(30,'exponents',2,5,dummy)
 if(Nprimgtf/=int(dummy))then
  write(6,fmt="('ERROR - number of primitive Gaussians does not agree with&
   & number of exponents specified in Fchk file.')")
  stop
 endif

! Primitive exponents - one for each shell
! There are now only five fields per line so recalculate the number
! of lines of data to read in...
 nline=max(Nprimgtf/5,1)

 ic=0
 do i=1,nline
  read(30,*)shexpnt(ic+1:min(ic+5,Nprimgtf))
  ic=ic+5
 enddo
! However, there may be a last line with less than 5 numbers...
 if((mod(Nprimgtf,5)/=0).and.(Nprimgtf>5))then
  read(30,*)shexpnt(ic+1:Nprimgtf)
 endif

!========================Contraction coefficients =======================

 write(*,fmt="('Reading the contraction coefficients...')")

 call getval(30,'Contraction',1,5,dummy)
 if(Nprimgtf/=int(dummy))then
  write(6,fmt="('ERROR - number of primitive Gaussians does not agree with&
   & number of contraction coefficients specified in Fchk file.')")
  stop
 endif

! Contraction coefficients - have as many of these as we have primitive
! Gaussian-type functions.

 ic=0
 do i=1,nline
  read(30,*)c_prim(ic+1:min(ic+5,Nprimgtf))
  ic=ic+5
 enddo
! However, there may be a last line with less than 5 numbers...
 if((mod(Nprimgtf,5)/=0).and.(Nprimgtf>5))then
  read(30,*)c_prim(ic+1:Nprimgtf)
 endif

! If there are SP shells then there is another block of information
! giving the contraction coefficients associated with the P functions
! of the SP shell.  (Those for the S are included with the normal
! contraction coefficients info.)
 if(num_sp>0)then
  read(30,*)

  ic=0
  do i=1,nline
   read(30,*)c2_prim(ic+1:min(ic+5,Nprimgtf))
   ic=ic+5
  enddo
! However, there may be a last line with less than 5 numbers...
  if((mod(Nprimgtf,5)/=0).and.(Nprimgtf>5))then
   read(30,*)c2_prim(ic+1:Nprimgtf)
  endif
 endif

! Include common normalization factors in the contraction coefficients
! and store in a referenceable form...
 call con_coeffs

!=====================Coordinates of each shell============================

 write(*,fmt="('Reading the Coordinates of each shell...')")

 call getval(30,'Coordinates',1,7,dummy)
 if((3*Nshells)/=int(dummy))then
  write(6,fmt="('ERROR - number of contracted shells does not agree with&
   & number of shell centres specified in Fchk file.')")
  stop
 endif

! Coordinates of each shell - (x,y,z) for `Nshells' shells
 nline=max((3*Nshells)/5,1)
 ic=0
 do i=1,nline
  read(30,*)Temp(ic+1:min(ic+5,3*Nshells))
  ic=ic+5
 enddo
! However, there may be a last line with less than 5 numbers...
 if((mod((3*Nshells),5)/=0).and.((3*Nshells)>5))then
  read(30,*)Temp(ic+1:(3*Nshells))
 endif

! Store in directly referenceable form...
 ic=0
 do i=1,Nshells
  shll_posn(1:3,i)=Temp(ic+1:ic+3)
  ic=ic+3
 enddo

!=========================Alpha MO expansion coefficients=================

 write(*,fmt="('Reading the alpha MO exp''n coeffs...')")

 call getval(30,'MO',2,6,dummy)
 if(Ncoeff2/=int(dummy))then
  write(6,fmt="('ERROR: (number of basis functions * no. of MO) does not agree&
   & with number of eigenvector coefficients  specified in Fchk file.')")
  stop
 endif

! Now the Alpha eigenvector coefficients.  We have Ncoeffs basis functions and
! Nmo MOs.  Therefore the total number of coefficients is Ncoeffs*Nmo=Ncoeff2,
! with five to a line.
 nline=max(Ncoeff2/5,1)
 ic=0
 do i=1,nline
  read(30,*)temp(ic+1:min(ic+5,Ncoeff2))
  ic=ic+5
 enddo
! However, there may be a last line with less than 5 numbers...
 if((mod(Ncoeff2,5)/=0).and.(Ncoeff2>5))then
  read(30,*)temp(ic+1:Ncoeff2)
 endif

!=========================Beta MO expansion coefficients=============

 if(SPIN)then
  write(*,fmt="('Reading the beta MO expansion coefficients...')")

  read(30,*)
! Exactly the same but for the Beta MO coefficients...
  ic=Ncoeff2
  do i=1,nline
   read(30,*)temp(ic+1:min(ic+5,(2*Ncoeff2)))
   ic=ic+5
  enddo
! However, there may be a last line with less than 5 numbers...
  if((mod(Ncoeff2,5)/=0).and.(Ncoeff2>5))then
   read(30,*)temp(ic+1:(2*Ncoeff2))
  endif
 endif

! Now store in a directly accessible form with alpha and beta MO coefficients
! stored separately.  Also multiply in the remaining
! normalization factors - required because the normalization of d and higher
! shells is not the same for all basis functions.  The common parts of their
! normalization are included in the contraction coefficients.
 call pack_evcoeffs(temp,2*Ncoeff2)

 close(unit=30)

END SUBROUTINE read_fchk

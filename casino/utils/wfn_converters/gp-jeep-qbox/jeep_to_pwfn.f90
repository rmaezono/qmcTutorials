!----------------------------------------------------------------------------!
! JEEP_TO_PWFN                                                               !
! ------------                                                               !
! Neil Drummond, 10.2003                                                     !
!                                                                            !
! This code reads in a jeep.wf JEEP wavefunction file and atoms.sys file and !
! generates a CASINO pwfn.data file.                                         !
!----------------------------------------------------------------------------!


MODULE constants
 IMPLICIT NONE
 DOUBLE PRECISION,PARAMETER :: pi=3.14159265358979324d0
 DOUBLE PRECISION,PARAMETER :: twopi=3.14159265358979324d0*2.d0
END MODULE constants


MODULE pseudo
!---------------------------------------------------------!
! Define the pseudo atoms used in the JEEP calculations.  !
!---------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,PARAMETER :: no_elements=92 ! For present purposes...
 CHARACTER(12) element_names(no_elements),element_names_caps(no_elements)
 DATA element_names/'hydrogen','helium','lithium','beryllium','boron','carbon',&
  &'nitrogen','oxygen','fluorine','neon','sodium','magnesium','aluminium',     &
  &'silicon','phosphorus','sulphur','chlorine','argon','potassium','calcium',  &
  &'scandium','titanium','vanadium','chromium','manganese','iron','cobalt',    &
  &'nickel','copper','zinc','gallium','germanium','astatine','selenium',       &
  &'bromine','krypton','rubidium','strontium','yttrium','zirconium','niobium', &
  &'molybdenum','technetium','ruthenium','rhodium','palladium','silver',       &
  &'cadmium','indium','tin','antimony','tellurium','iodine','xenon','caesium', &
  &'barium','lanthanum','cerium','praesodymium','neodymium','promethium',      &
  &'samarium','europium','gadolinium','terbium','dysprosium','holmium',        &
  &'erbium','thulium','ytterbium','lutetium','hafnium','tantalum','tungsten',  &
  &'rhenium','osmium','iridium','platinum','gold','mercury','thallium','lead', &
  &'bismuth','polonium','astatine','radon','francium','radium','actinium',     &
  &'thorium','protactinium','uranium'/
 DATA element_names_caps /'Hydrogen','Helium','Lithium','Beryllium','Boron',   &
  &'Carbon','Nitrogen','Oxygen','Fluorine','Neon','Sodium','Magnesium',        &
  &'Aluminium','Silicon','Phosphorus','Sulphur','Chlorine','Argon','Potassium',&
  &'Calcium','Scandium','Titanium','Vanadium','Chromium','Manganese','Iron',   &
  &'Cobalt','Nickel','Copper','Zinc','Gallium','Germanium','Astatine',         &
  &'Selenium','Bromine','Krypton','Rubidium','Strontium','Yttrium','Zirconium',&
  &'Niobium','Molybdenum','Technetium','Ruthenium','Rhodium','Palladium',      &
  &'Silver','Cadmium','Indium','Tin','Antimony','Tellurium','Iodine','Xenon',  &
  &'Caesium','Barium','Lanthanum','Cerium','Praesodymium','Neodymium',         &
  &'Promethium','Samarium','Europium','Gadolinium','Terbium','Dysprosium',     &
  &'Holmium','Erbium','Thulium','Ytterbium','Lutetium','Hafnium','Tantalum',   &
  &'Tungsten','Rhenium','Osmium','Iridium','Platinum','Gold','Mercury',        &
  &'Thallium','Lead','Bismuth','Polonium','Astatine','Radon','Francium',       &
  &'Radium','Actinium','Thorium','Protactinium','Uranium'/

END MODULE pseudo


PROGRAM jeep_to_pwfn
!-----------------!
! Main program.   !
!-----------------!
 USE constants
 USE gen_gvec_module
 USE pseudo
 IMPLICIT NONE
 INTEGER i,j,ng_temp,ierr,natoms,nstates,ng_half,istate
 INTEGER,ALLOCATABLE :: atom_Z(:),g_int(:,:)
 DOUBLE PRECISION b1(3),b2(3),b3(3),rtemp(3),ecut,ecut_Ry,ref_a,ref_b, &
  &ref_c,a,b,c
 DOUBLE PRECISION,ALLOCATABLE :: atom_r(:,:)
 COMPLEX(KIND=KIND(0.d0)),ALLOCATABLE :: pw_coeff(:)
 CHARACTER(4) ctemp,ctemp7,ctemp2
 CHARACTER(12) ptemp
 LOGICAL :: jeep_wf_exists
 DOUBLE PRECISION,PARAMETER :: ecut_max=500

 write(6,*)
 write(6,*)'JEEP_TO_PWFN'
 write(6,*)

 inquire(file='jeep.wf',exist=jeep_wf_exists)
 if(.not.jeep_wf_exists)then
  write(6,*)'You need to supply a jeep.wf file with the jeep wfn.'
  stop
 endif

! Find out how many states there should be and how many to write out.
 write(6,*)'Reading jeep.wf'
 call readwf1(nstates,ng_temp)

 write(6,'(" There are ",I4," states, with ",I8," G-vectors")') &
  &nstates,ng_temp
 write(6,*)

! Get the atomic positions and types from the JEEP atoms.sys file
 write(6,*)'Reading atoms.sys.'
 write(6,*)
 open(12,file='atoms.sys',form='formatted',status='old',iostat=ierr)
 if(ierr/=0)then
  write(6,*)'Problem opening atoms.sys.'
  stop
 endif
 read(12,*)ctemp7,a,b,c
 write(6,'(" Dimensions of simulation box: ",f8.3," by ",f8.3," by ",f8.3)')a,&
  &b,c

! Count atoms in atoms.sys.
 natoms=0
 do
  read(12,*,iostat=i)ctemp,ctemp2,ptemp,rtemp
  if(i/=0)exit
  natoms=natoms+1
 enddo

 write(6,'(" Number of atoms in atoms.sys: ",i4)')natoms
 write(6,*)

 rewind(12)
 read(12,*)ctemp7,a,b,c

 allocate(atom_r(3,natoms),atom_Z(natoms),stat=ierr)
 if(ierr/=0)then
  write(6,*)'Allocation problem [3.5].'
  stop
 endif

! Read atomic positions and subtract off CoM. Work out what
! atomic number the atoms have.
 write(6,*)'Atomic positions:'
 do i=1,natoms
  read(12,*)ctemp,ctemp2,ptemp,atom_r(1:3,i)
  atom_Z(i)=-1
  do j=1,no_elements
   if(index(ptemp,trim(adjustl(element_names(j))))/=0)then
    if(atom_Z(i)==-1)then
     atom_Z(i)=j
    else
     if(index(trim(adjustl(element_names(j))), &
      &trim(adjustl(element_names(atom_Z(i)))))/=0)then
      atom_Z(i)=j
     elseif(index(trim(adjustl(element_names(atom_Z(i)))), &
      &trim(adjustl(element_names(j))))==0)then
      write(6,*)'Can''t decide what element atom ',ctemp2,' is.'
      write(6,*)'It contains more than one element name!'
      stop
     endif
    endif
   endif
   if(index(ptemp,trim(adjustl(element_names_caps(j))))/=0)then
    if(atom_Z(i)==-1)then
     atom_Z(i)=j
    else
     if(index(trim(adjustl(element_names_caps(j))), &
      &trim(adjustl(element_names(atom_Z(i)))))/=0)then
      atom_Z(i)=j
     elseif(index(trim(adjustl(element_names_caps(atom_Z(i)))), &
      &trim(adjustl(element_names(j))))==0)then
      write(6,*)'Can''t decide what element atom ',ctemp2,' is.'
      write(6,*)'It contains more than one element name!'
      stop
     endif
    endif
   endif
  enddo
  if(atom_Z(i)==-1)then
   write(6,*)'Haven''t been able to find the atomic number for atom ',ctemp2
   write(6,*)'Check the name of your pseudopotential. It should contain'
   write(6,*)'the name of the desired element in lower case.'
   stop
  endif
  write(6,'(" atom ",A4," ",A12,3(1X,F15.8))')ctemp2,ptemp, &
   &atom_r(1:3,i)
 enddo
 write(6,*)

 close(12)

 ref_a=a ; ref_b=b ; ref_c=c

! Guess the cutoff and generate the g-vector list
 ecut_Ry=0.d0
 allocate(g_int(3,ng_temp),stat=ierr)
 if(ierr/=0)then
  write(6,*)'Allocation problem [4].'
  stop
 endif

 do i=1,int(ecut_max)
  call gen_gvec(dble(i),ref_a,ref_b,ref_c,ng_half,g_int,.true.)
  if(ng_half==ng_temp)then
   ecut_Ry=dble(i)
   write(6,'(" It looks like the Ecut is ",F6.1," Rydbergs.")')ecut_Ry
   exit
  endif
 enddo
 if(ecut_Ry==0.d0)then
  write(6,*)'Ecut value cannot be determined.'
  stop
 endif

 call gen_gvec(ecut_Ry,ref_a,ref_b,ref_c,ng_half,g_int,.false.)
 write(6,'(" There are ",I8," g-vectors in the original box")')ng_half
 write(6,*)

! Convert ecut to au
 ecut=0.5d0*ecut_Ry

! Reciprocal lattice vectors.
 b1(1:3)=0.d0 ; b2(1:3)=0.d0 ; b3(1:3)=0.d0
 b1(1)=twopi/a
 b2(2)=twopi/b
 b3(3)=twopi/c

! Allocate the storage of the G-vector coefficients for one spline
 allocate(pw_coeff(ng_half),stat=ierr)
 if(ierr/=0)then
  write(6,*)'Allocation problem [6].'
  stop
 endif

 open(unit=9,file='pwfn.data',status='replace',iostat=ierr)
 if(ierr/=0)then
  write(6,*)'Problem opening pwfn.data.'
  stop
 endif

 write(9,*)'Please replace this line with a title.'
 write(9,*)
 write(9,*)'BASIC INFO'
 write(9,*)'----------'
 write(9,*)'Generated by:'
 write(9,*)'JEEP'
 write(9,*)'Method:'
 write(9,*)'DFT'
 write(9,*)'DFT Functional'
 write(9,*)'unknown'
 write(9,*)'Pseudopotential'
 write(9,*)'unknown'
 write(9,*)'Plane wave cutoff (au)'
 write(9,*)ecut
 write(9,*)'Spin polarized:'
 write(9,*)'F'
 write(9,*)'Total energy (au per primitive cell)'
 write(9,*)0.d0
 write(9,*)'Kinetic energy (au per primitive cell)'
 write(9,*)0.d0
 write(9,*)'Local potential energy (au per primitive cell)'
 write(9,*)0.d0
 write(9,*)'Non-local potential energy (au per primitive cell)'
 write(9,*)0.d0
 write(9,*)'Electron-electron energy (au per primitive cell)'
 write(9,*)0.d0
 write(9,*)'Ion-ion energy (au per primitive cell)'
 write(9,*)0.d0
 write(9,*)'Number of electrons per primitive cell'
 write(9,*)0
 write(9,*)
 write(9,*)'GEOMETRY'
 write(9,*)'--------'
 write(9,*)'Number of atoms per primitive cell'
 write(9,*)natoms
 write(9,*)'Atomic numbers and positions of atoms (au)'
 do i=1,natoms
  write(9,'(" ",i2,"  ",f15.8," ",f15.8," ",f15.8)')atom_Z(i),atom_r(1:3,i)
 enddo
 write(9,*)'Primitive lattice vectors (au)'
 write(9,'(" ",f15.8," ",f15.8," ",f15.8)')a,0.d0,0.d0
 write(9,'(" ",f15.8," ",f15.8," ",f15.8)')0.d0,b,0.d0
 write(9,'(" ",f15.8," ",f15.8," ",f15.8)')0.d0,0.d0,c
 write(9,*)
 write(9,*)'G VECTORS'
 write(9,*)'---------'
 write(9,*)'Number of G-vectors'
 write(9,*)2*ng_half-1
 write(9,*)'Gx Gy Gz (au)'
 write(9,'(" ",f15.8," ",f15.8," ",f15.8)')0.d0,0.d0,0.d0
 do i=2,ng_half
  write(9,'(" ",f15.8," ",f15.8," ",f15.8)')g_int(1,i)*b1(1:3)+ &
   &g_int(2,i)*b2(1:3)+g_int(3,i)*b3(1:3)
  write(9,'(" ",f15.8," ",f15.8," ",f15.8)')-g_int(1,i)*b1(1:3) &
   &-g_int(2,i)*b2(1:3)-g_int(3,i)*b3(1:3)
 enddo
 write(9,*)
 write(9,*)'WAVE FUNCTION'
 write(9,*)'-------------'
 write(9,*)'Number of k-points'
 write(9,*)1
 write(9,*)'k-point # ; # of bands (up spin/down spin) ; k-point coords (au)'
 write(9,'(" ",i1," ",i4," ",i1," ",f6.1," ",f6.1," ",f6.1)')1,nstates,0, &
  &0.d0,0.d0,0.d0
 do istate=1,nstates
  write(9,*)'Band,spin,eigenvalue (au)'
  write(9,'(" ",i4," ",i1," ",f15.8)')istate,1,dble(istate)
  write(9,*)'Eigenvector coefficients'
  call readwf2(pw_coeff,istate)
  write(9,*)pw_coeff(1)
  do i=2,ng_half
   write(9,*)pw_coeff(i)
   write(9,*)conjg(pw_coeff(i))
  enddo
 enddo

 write(6,*)'The pwfn.data file is generated.'
 write(6,*)'Please insert energies and titles etc at the top of this file.'
 write(6,*)'Note that the "eigenvalues" aren''t the actual eigenvalues!'
 write(6,*)
 write(6,*)'Program finished.'
 write(6,*)

END PROGRAM jeep_to_pwfn

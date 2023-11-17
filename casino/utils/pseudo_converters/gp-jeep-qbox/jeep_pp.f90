!------------------------------------------------------------!
! JEEP_PP                                                    !
! -------                                                    !
!                                                            !
! AJW, 2003                                                  !
!                                                            !
! Convert a pseudopotential file used by JEEP into the       !
! format used by CASINO.                                     !
!                                                            !
! Changes:                                                   !
! --------                                                   !
! 10.2003 NDD - Generalized and updated.                     !
! 10.2004 NDD - Write out [element symbol]_pp.data.          !
!------------------------------------------------------------!


MODULE pseudo
!---------------------------------------------------------!
! Define the pseudo atoms used in the JEEP calculations.  !
!---------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,PARAMETER :: no_elements=92
 CHARACTER(12) element_names(no_elements),element_names_caps(no_elements)
 CHARACTER(2) periodic_table_nocap(no_elements)
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
 DATA periodic_table_nocap/ &
  & 'h ','he','li','be','b ','c ','n ',                    &
  & 'o ','f ','ne','na','mg','al','si','p ','s ','cl','ar',&
  & 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',&
  & 'zn','ga','ge','as','se','br','kr','rb','sr','y ','zr',&
  & 'nb','mo','tc','ru','rh','pd','ag','cd','in','sn','sb',&
  & 'te','i ','xe','cs','ba','la','ce','pr','nd','pm','sm',&
  & 'eu','gd','tb','dy','ho','er','tm','yb','lu','hf','ta',&
  & 'w ','re','os','ir','pt','au','hg','tl','pb','bi','po',&
  & 'at','rn','fr','ra','ac','th','pa','u '/
END MODULE pseudo


PROGRAM jeep_pp
 USE pseudo
 IMPLICIT NONE

 INTEGER i,j,num_comment_lines,num_grid,valence_charge,l_local,ialloc,atomic_no
 DOUBLE PRECISION,ALLOCATABLE :: r(:),v_s(:),v_p(:),v_d(:)
 DOUBLE PRECISION temp_r,temp_phi,rcut,mass
 CHARACTER(1) char
 CHARACTER(72) in_file,out_file,comment

 write(6,*)
 write(6,*)'JEEP_TO_CASINO_PSEUDO'
 write(6,*)
 write(6,*)'Please enter name of input (JEEP pseudopotential) file.'
 read(5,*)in_file
 write(6,*)
 open(unit=8,file=trim(in_file),status='old',iostat=i)
 if(i/=0)then
  write(6,*)'Can''t open file ',trim(adjustl(in_file)),'.'
  stop
 endif

 atomic_no=-1
 do j=1,no_elements
  if(index(in_file,trim(adjustl(element_names(j))))/=0)then
   if(atomic_no==-1)then
    atomic_no=j
   else
    if(index(trim(adjustl(element_names(j))), &
     &trim(adjustl(element_names(atomic_no))))/=0)then
     atomic_no=j
    elseif(index(trim(adjustl(element_names(atomic_no))), &
     &trim(adjustl(element_names(j))))==0)then
     write(6,*)'Can''t decide which element this pseudopotential is for.'
     write(6,*)'It contains more than one element name!'
     stop
    endif
   endif
  endif
  if(index(in_file,trim(adjustl(element_names_caps(j))))/=0)then
   if(atomic_no==-1)then
    atomic_no=j
   else
    if(index(trim(adjustl(element_names_caps(j))), &
     &trim(adjustl(element_names(atomic_no))))/=0)then
     atomic_no=j
    elseif(index(trim(adjustl(element_names_caps(atomic_no))), &
     &trim(adjustl(element_names(j))))==0)then
     write(6,*)'Can''t decide which element this pseudopotential is for.'
     write(6,*)'It contains more than one element name!'
     stop
    endif
   endif
  endif
 enddo
 if(atomic_no==-1)then
  write(6,*)'Haven''t been able to find the atomic number.'
  write(6,*)'Check the name of your pseudopotential. It should contain'
  write(6,*)'the name of the desired element in lower case.'
  stop
 endif

 out_file=trim(periodic_table_nocap(atomic_no))//'_pp.data'
 write(6,*)'Writing pseudopotential to ',trim(out_file),'.'
 write(6,*)
 open(unit=9,file=trim(out_file),status='new',iostat=i)
 if(i/=0)then
  write(6,*)'Can''t open file ',trim(adjustl(out_file)),'.'
  write(6,*)'NB, won''t overwrite existing files.'
  stop
 endif

! Read comments at top of jeep pseudopotential.
 num_comment_lines=0
 do
  read(8,*)char
  if(char=="#")then
   num_comment_lines=num_comment_lines+1
  else
   exit
  endif
 enddo
 rewind(8)
 write(6,*)'Number of comment lines : ',num_comment_lines
 do i=1,num_comment_lines
  read(8,*)
 enddo

 read(8,*)num_grid,valence_charge,rcut,mass,l_local
 l_local=l_local-1
 write(6,*)'l_local=',l_local

 allocate(r(num_grid),v_s(num_grid),v_p(num_grid),v_d(num_grid),stat=ialloc)
 if(ialloc/=0)then
  write(6,*)'Allocation problem'
  stop
 endif

! Read s part.
 do i=1,num_grid
  read(8,*)r(i),v_s(i)
 enddo
 read(8,*,end=101)
 do i=1,num_grid
  read(8,*)temp_r,temp_phi
 enddo

! Read p part.
 read(8,*,end=101)
 do i=1,num_grid
  read(8,*)r(i),v_p(i)
 enddo
 read(8,*)
 do i=1,num_grid
  read(8,*)temp_r,temp_phi
 enddo

! read d part
 read(8,*,end=102)
 do i=1,num_grid
  read(8,*)r(i),v_d(i)
 enddo

101 continue
102 continue

 write(6,*)'Please give header/description for CASINO pseudopot file.'
 write(6,*)'No more than 72 characters, please.'
 read(5,'(a)')comment
 write(6,*)

 write(9,*)trim(adjustl(comment))
 write(9,*)'Atomic number and pseudo-charge'
 write(9,'(" ",i3,"  ",f6.1)')atomic_no,dble(valence_charge)
 write(9,*)'Energy units (rydberg/hartree/ev):'
 write(9,*)'rydberg'
 write(9,*)'Angular momentum of local component (0=s,1=p,2=d..) for DFT &
  &and QMC'
 write(9,'(i2,1x,i2)')l_local,l_local
 write(9,*)'NLRULE override (1) VMC/DMC (2) config gen &
  &(0 ==> input/default value)'
 write(9,*)'0 0'
 write(9,*)'Number of grid points'
 write(9,'(i12)')num_grid
 write(9,*)'R(i) in atomic units'

! Print real-space grid
 do i=1,num_grid
  write(9,*)r(i)
 enddo

 if(l_local==0)then

! Write the s part for s
  write(9,*)'r*potential (L=0) in Ry'
  do i=1,num_grid
   write(9,*)2.d0*r(i)*v_s(i)
!            ^^^^
!       Convert to Ry
  enddo
! Write the s part for p
  write(9,*)'r*potential (L=1) in Ry'
  do i=1,num_grid
   write(9,*)2.d0*r(i)*v_s(i)
!            ^^^^
!       Convert to Ry
  enddo
! Write the s part for d
  write(9,*)'r*potential (L=2) in Ry'
  do i=1,num_grid
   write(9,*)2.d0*r(i)*v_s(i)
!            ^^^^
!       Convert to Ry
  enddo

 elseif(l_local==1)then

! Write the s part for s
  write(9,*)'r*potential (L=0) in Ry'
  do i=1,num_grid
   write(9,*)2.d0*r(i)*v_s(i)
!            ^^^^
!       Convert to Ry
  enddo
! Write the p part for p
  write(9,*)'r*potential (L=1) in Ry'
  do i=1,num_grid
   write(9,*)2.d0*r(i)*v_p(i)
!            ^^^^
!       Convert to Ry
  enddo
! Write the p part for d
  write(9,*)'r*potential (L=2) in Ry'
  do i=1,num_grid
   write(9,*)2.d0*r(i)*v_p(i)
!            ^^^^
!       Convert to Ry
  enddo

 elseif(l_local==2)then

! Write the s part for s
  write(9,*)'r*potential (L=0) in Ry'
  do i=1,num_grid
   write(9,*)2.d0*r(i)*v_s(i)
!            ^^^^
!       Convert to Ry
  enddo
! Write the p part for p
  write(9,*)'r*potential (L=1) in Ry'
  do i=1,num_grid
   write(9,*)2.d0*r(i)*v_p(i)
!            ^^^^
!       Convert to Ry
  enddo
! Write the d part for d
  write(9,*)'r*potential (L=2) in Ry'
  do i=1,num_grid
   write(9,*)2.d0*r(i)*v_d(i)
!            ^^^^
!       Convert to Ry
  enddo

 else

  write(6,*)'Problem with l_local = ',l_local
  stop

 endif ! value of l_local

 deallocate(r,v_s,v_p,v_d)

 write(6,*)'Program finished.'
 write(6,*)

END PROGRAM jeep_pp

PROGRAM ion_dist
!------------------------------------------------------------!
! Program for automating the generation of the edist_by_ion  !
! block in the input file for CASINO. Works for              !
! antiferromagnetic Wigner crystals.                         !
!------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,n,no_electrons_prim,no_prim_cells,total_no_k_points,no_real_k_points
 INTEGER,ALLOCATABLE :: spin_dat(:)
 CHARACTER(72) dummy
 CHARACTER(25) dummy_25
! CHARACTER(5+13*no_electrons_prim) fmt_spin_in
 CHARACTER(10000) fmt_spin_in
 LOGICAL have_gwfn

 write(6,*)
 write(6,*)'ION_DIST'
 write(6,*)

 inquire(file='gwfn.data',exist=have_gwfn)
 if(.not.have_gwfn)then
  write(6,*)'The gwfn.data file is missing. Stopping.'
  stop
 endif

 open(unit=8,file='gwfn.data',status='old',iostat=i)
 if(i/=0)then
  write(6,*)'Sorry, can''t open gwfn.data.'
  stop
 endif

! Scan to find line "K SPACE NET". NB Can't use Monkhorst-Pack mesh
! given in crystal output because data may have been "plucked".
 do
  read(8,'(a)',err=90,end=90)dummy
  if(index(dummy,'K SPACE NET')/=0)exit
 enddo
 read(8,*)dummy
 read(8,*)dummy
 read(8,*)total_no_k_points
 read(8,*)dummy
 read(8,*)no_real_k_points
 no_prim_cells=2*(total_no_k_points-no_real_k_points)+no_real_k_points

 write(6,10)'Supercell consists of ',no_prim_cells,' primitive cells.'

! Scan to find line "NUMBER OF A", just before no. electrons in prim cell.
 do
  read(8,'(a)',err=100,end=100)dummy
  if(index(dummy,'NUMBER OF AO')/=0)exit
 enddo

 read(8,'(a,i5)')dummy_25,no_electrons_prim
 write(6,10)'Primitive cell contains ',no_electrons_prim,' electrons.'
 allocate(spin_dat(no_electrons_prim))

! Scan to find line "ATOMIC SPINS SET TO", just before line with spin data.
 do
  read(8,'(a)',err=110,end=110)dummy
  if(index(dummy,'ATOMIC SPINS SET TO')/=0)exit
 enddo

 call read_in_spins
 write(6,*)'Read in spins of electrons in primitive lattice.'

 close(8)

! Write out data in form required for edist_by_ion input block.
 open(unit=9,file='ion_dist.out')
 do n=1,no_prim_cells
  do i=1,no_electrons_prim
   if(spin_dat(i)==1)then
    write(9,*)(n-1)*no_electrons_prim+i,1,0
   else
    write(9,*)(n-1)*no_electrons_prim+i,0,1
   endif
  enddo
 enddo
 close(9)

 write(6,*)'edist_by_ion data placed in "ion_dist.out".'

 deallocate(spin_dat)

 write(6,*)
 write(6,*)'Program finished.'
 write(6,*)

 stop

10 format(' ',a,i4,a)
90 write(6,*)'Couldn''t find line "K SPACE NET" in gwfn.data.'
 stop
100 write(6,*)'Couldn''t find line "NUMBER OF AO" in gwfn.data.'
 stop
110 write(6,*)'Couldn''t find line "ATOMIC SPINS SET TO" in gwfn.data.'
 stop


CONTAINS


 SUBROUTINE read_in_spins
! Read in (ATOM, AT. N., SPIN) data.
 IMPLICIT NONE
 INTEGER i
 INTEGER,ALLOCATABLE :: spin_dat_in(:)

 fmt_spin_in='('
 do i=1,no_electrons_prim
  fmt_spin_in=trim(fmt_spin_in)//'" ",i3,i4,i2,'
 enddo
 fmt_spin_in=trim(fmt_spin_in)//'" ")'

 allocate(spin_dat_in(3*no_electrons_prim))
 read(8,trim(fmt_spin_in))spin_dat_in(1:3*no_electrons_prim)
 do i=1,no_electrons_prim
  spin_dat(i)=spin_dat_in(3*i)
 enddo
 deallocate(spin_dat_in)

 END SUBROUTINE read_in_spins


END PROGRAM ion_dist


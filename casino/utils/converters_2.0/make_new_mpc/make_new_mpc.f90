PROGRAM make_new_mpc
!--------------------------------------------------------------------------!
! Construct new format mpc.data file by merging old density.data and       !
! eepot.data files with slight format changes.                             !
!                                                                          !
! MDT 5.2005                                                               !
!--------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,ierr,nat,den_ngvec,eepot_ngvec,io_density,io_mpc,io_eepot,n1,n2,n3
 INTEGER,ALLOCATABLE :: atno(:)
 INTEGER,PARAMETER :: dp=kind(1.d0),file_version=1,r2s_length=80
 REAL(dp) pa1(3),pa2(3),pa3(3),a1(3),a2(3),a3(3),treal,timag
 REAL(dp),ALLOCATABLE :: den_gvec(:,:),eepot_gvec(:,:),basis(:,:),eepot_sc(:)
 COMPLEX(KIND=KIND(0d0)),ALLOCATABLE :: den_sc(:)
 LOGICAL eepot_present,density_present
 CHARACTER(80) :: ctext
 CHARACTER(r2s_length) tmpr,tmpr2,tmpr3,tmpr4

! Check existence of relevant files and open them

 inquire(file='density.data',exist=density_present)
 inquire(file='eepot.data',exist=eepot_present)
 if(.not.density_present.or..not.eepot_present)then
  write(6,*)'MAKE_NEW_MPC will create a new format mpc.data from the old format'
  write(6,*)'eepot.data and density.data files, but only if they exist.'
  write(6,*)
  stop
 endif

 io_density=10 ; io_eepot=11 ; io_mpc=12
 open(io_density,file='density.data',status='old',action='read',err=10)
 open(io_eepot,file='eepot.data',status='old',action='read',err=20)
 open(io_mpc,file='mpc.data',status='unknown',action='write',err=30)

! Read the density.data file

! Primitive lattice translation vectors
 read(io_density,'(a)',err=40,end=40)ctext
 call skip(io_density,2)
 read(io_density,*,err=40)pa1
 read(io_density,*,err=40)pa2
 read(io_density,*,err=40)pa3
 call skip(io_density,5)

! Number of atoms in basis
 read(io_density,*,err=40)nat
 call skip(io_density,1)

! Positions of basis atoms (a.u.)
 allocate(basis(3,nat),atno(nat),stat=ierr)
 do i=1,nat
  read(io_density,*,err=40)atno(i),(basis(j,i),j=1,3)
 enddo
 call skip(io_density,1)

! Number of primitive cell G-vectors.
 read(io_density,*,err=40)den_ngvec

! G-vectors

 allocate(den_gvec(3,den_ngvec),den_sc(den_ngvec),stat=ierr)
 if(ierr/=0)then
  write(6,*)'Allocation problem.' ; stop
 endif

 call skip(io_density,1)
 do i=1,den_ngvec
  read(io_density,*,err=40)den_gvec(1,i),den_gvec(2,i),den_gvec(3,i)
 enddo

! Self-consistent density

 call skip(io_density,1)
 do i=1,den_ngvec
  read(io_density,*,err=40,end=40)treal,timag
  den_sc(i)=cmplx(treal,timag,kind=dp)
 enddo

 close(io_density)

! Write the density into the mpc.data file

 write(io_mpc,'(a)')'START MPC DATA'
 write(io_mpc,'(a)')'Title'
 write(io_mpc,'(a)')'  No title given'
 write(io_mpc,'(a)')'File version'
 write(io_mpc,'(1x,a)')trim(i2s(file_version))
 write(io_mpc,*)

 write(io_mpc,'(a)')'START DENSITY DATA'
 write(io_mpc,'(a)')'Self consistent charge density in reciprocal space'
 write(io_mpc,'(a)')'Real space primitive cell translation vectors (au)'

 tmpr=r2s(pa1(1),'(f20.15)')
 tmpr2=r2s(pa1(2),'(f20.15)')
 tmpr3=r2s(pa1(3),'(f20.15)')
 write(io_mpc,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)
 tmpr=r2s(pa2(1),'(f20.15)')
 tmpr2=r2s(pa2(2),'(f20.15)')
 tmpr3=r2s(pa2(3),'(f20.15)')
 write(io_mpc,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)
 tmpr=r2s(pa3(1),'(f20.15)')
 tmpr2=r2s(pa3(2),'(f20.15)')
 tmpr3=r2s(pa3(3),'(f20.15)')
 write(io_mpc,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)

 write(io_mpc,'(a)')'Number of atoms in the primitive cell'
 write(io_mpc,'(1x,a)')trim(i2s(nat))
 write(io_mpc,'(a)')'Positions of atoms (au)'
 do i=1,nat
  tmpr=r2s(basis(1,i),'(f20.15)')
  tmpr2=r2s(basis(2,i),'(f20.15)')
  tmpr3=r2s(basis(3,i),'(f20.15)')
  write(io_mpc,'(1x,4(1x,a))')trim(i2s(atno(i))),trim(tmpr),trim(tmpr2), &
   &trim(tmpr3)
 enddo

 write(io_mpc,'(a)')'Energy cutoff used for G-vectors (au)'
 write(io_mpc,'(a)')'  Unknown'
 write(io_mpc,'(a)')'Number of particle types (1=electrons,2=electrons/holes)'
 write(io_mpc,*)' 1'
 write(io_mpc,'(a)')'Number of G-vectors'
 write(io_mpc,'(1x,a)')trim(i2s(den_ngvec))
 write(io_mpc,'(a)')'G-vectors (au)'

 do i=1,den_ngvec
  tmpr=r2s(den_gvec(1,i),'(f20.14)')
  tmpr2=r2s(den_gvec(2,i),'(f20.14)')
  tmpr3=r2s(den_gvec(3,i),'(f20.14)')
  write(io_mpc,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)
 enddo

 write(io_mpc,'(a)')'START SET 1'
 write(io_mpc,'(a)')'Particle type'
 write(io_mpc,'(a)')'  1'
 write(io_mpc,'(a)')'Complex charge density (real part, imaginary part)'
 do i=1,den_ngvec
  tmpr=r2s(dble(den_sc(i)),'(f20.14)')
  tmpr2=r2s(aimag(den_sc(i)),'(f20.14)')
  write(io_mpc,'(2(1x,a))')trim(tmpr),trim(tmpr2)
 enddo
 write(io_mpc,'(a)')'END SET 1'
 write(io_mpc,'(a)')'END DENSITY DATA'
 write(io_mpc,*)

 deallocate(basis,den_gvec,den_sc)

! Read the eepot.data file

! Simulation cell translation vectors
 read(io_eepot,fmt='(a)',err=50,end=50)ctext
 call skip(io_eepot,1)
 read(io_eepot,*,err=50)a1
 read(io_eepot,*,err=50)a2
 read(io_eepot,*,err=50)a3
 call skip(io_eepot,1)

! Multiples of primitive translation vectors
 read(io_eepot,*,err=50)n1,n2,n3
 call skip(io_eepot,4)

! Number of simulation cell G-vectors.
 read(io_eepot,*,err=50)eepot_ngvec
 call skip(io_eepot,1)

! G vectors and coefficients

 allocate(eepot_gvec(3,eepot_ngvec),eepot_sc(eepot_ngvec),stat=ierr)
 if(ierr/=0)then
  write(6,*)'Allocation problem <2>.' ; stop
 endif

 do i=1,eepot_ngvec
  read(io_eepot,*,err=50)eepot_gvec(1,i),eepot_gvec(2,i),eepot_gvec(3,i), &
   &eepot_sc(i)
 enddo

 close(io_eepot)

! Write eepot.data into the mpc.data file

 write(io_mpc,'(a)')'START EEPOT DATA'
 write(io_mpc,'(a)')'Real space simulation cell translation vectors (au)'
 tmpr=r2s(a1(1),'(f23.15)')
 tmpr2=r2s(a1(2),'(f23.15)')
 tmpr3=r2s(a1(3),'(f23.15)')
 write(io_mpc,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)
 tmpr=r2s(a2(1),'(f23.15)')
 tmpr2=r2s(a2(2),'(f23.15)')
 tmpr3=r2s(a2(3),'(f23.15)')
 write(io_mpc,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)
 tmpr=r2s(a3(1),'(f23.15)')
 tmpr2=r2s(a3(2),'(f23.15)')
 tmpr3=r2s(a3(3),'(f23.15)')
 write(io_mpc,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)
 write(io_mpc,'(a)')'Multiples of primitive translation vectors'
 write(io_mpc,'(1x,3(1x,a))')trim(i2s(n1)),trim(i2s(n2)),trim(i2s(n3))
 write(io_mpc,'(a)')'Energy cutoff used for G-vectors (au)'
 write(io_mpc,'(a)')'  Unknown'
 write(io_mpc,'(a)')'Number of G-vectors'
 write(io_mpc,'(1x,a)')trim(i2s(eepot_ngvec))
 write(io_mpc,'(a)')'G-vectors (au) and Fourier coefficients'
 do i=1,eepot_ngvec
  tmpr=r2s(eepot_gvec(1,i),'(f17.12)')
  tmpr2=r2s(eepot_gvec(2,i),'(f17.12)')
  tmpr3=r2s(eepot_gvec(3,i),'(f17.12)')
  tmpr4=r2s(eepot_sc(i),'(f23.12)')
  write(io_mpc,'(4(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3),trim(tmpr4)
 enddo
 write(io_mpc,'(a)')'END EEPOT DATA'
 write(io_mpc,'(a)')

 write(io_mpc,'(a)')'END MPC DATA'

 close(io_mpc)

 stop

! Error handling

10 write(6,*)'Error opening density.data file' ; stop
20 write(6,*)'Error opening eepot.data file' ; stop
30 write(6,*)'Error opening new mpc.data file' ; stop
40 write(6,*)'Error reading density.data file' ; stop
50 write(6,*)'Error reading eepot.data file' ; stop


CONTAINS



 CHARACTER(20) FUNCTION i2s(n)
!-----------------------------------------------------------------------!
! I2S                                                                   !
! ===                                                                   !
! Convert integers to left justified strings that can be printed in the !
! middle of a sentence without introducing large amounts of white space.!
!                                                                       !
! Calling routine is intended to include something like:                !
! USE utilities                                                         !
! INTEGER i                                                             !
! i=12                                                                  !
! write(6,*)'Integer number ',trim(i2s(i)),' with words at the end.'    !
!-----------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: n
  INTEGER i,j
  CHARACTER tmp,sign

  if(n==0)then
   i2s='0' ; return
  endif
  sign=' ' ; if(n<0)sign='-'

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

  i2s=trim(sign)//i2s

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
! USE utilities                                                           !
! REAL(dp) r                                                              !
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
 REAL(dp),INTENT(in) :: r
 CHARACTER(*),INTENT(in) :: real_format

 write(r2s,real_format)r
 if(r<0)then
  r2s=adjustl(r2s)
 else
  r2s=' '//adjustl(r2s)
 endif

 END FUNCTION r2s


 SUBROUTINE skip(iunit,nskip)
!-------------------------------------!
! Skip records in a free format file. !
!-------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: iunit,nskip
 INTEGER i
 do i=1,nskip
  read(iunit,*)
 enddo
 END SUBROUTINE skip


END PROGRAM make_new_mpc


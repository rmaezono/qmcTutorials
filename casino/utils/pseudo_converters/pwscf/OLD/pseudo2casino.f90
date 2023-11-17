! ----------------------------------------------
! pseudopotential converter PWSCF to CASINO
! Dario Alfe`
!-----------------------------------------------
program pseudo2casino

implicit none

character*20 :: file_pseudo, file_casino, cdum
logical :: existing, l, l1
real*8, allocatable :: s(:)
real*8 :: dx, r, x0, z, xval
integer :: ir, nx, i, n, nloc

5    print '('' file pseudopot. > '',$)'
     read '(a)', file_pseudo
     inquire (file=file_pseudo,exist=existing)
     if (.not. existing) then
       print '('' File '',a,'' not found'')', file_pseudo
       go to 5
     end if
     print '('' file casino > '',$)'
     read '(a)', file_casino

open (1, file=file_pseudo, status='old')
open (2, file=file_casino, status='unknown')

read(1,*)cdum
write(2,*)cdum
write(2,'('' Atomic number and pseudo-charge'')')
read(1,*)cdum,xval,n,n,n,l,nloc,l
read(1,*)z,x0,dx,nx,n
allocate(s(nx))
write(2,'(1x,i2,f4.0)') int(z), xval
write(2,'('' Energy units (rydberg/hartree/ev):'')')
write(2,'('' rydberg'')')
write(2,'('' Angular momentum of local component (0=s,1=p,2=d..)'')')
write(2,'(1x,i1)')nloc
write(2,'('' NLRULE override (1) VMC/DMC (2) config gen (0 ==> input/default value)'')')
write(2,'('' 0 0'')')
write(2,'('' Number of grid points'')')
write(2,'(i12)') nx+1
write(2,'('' R(i) in atomic units'')')
r = 0.0d0
write(2,'(g23.16)') r
do ir = 1, nx
   r = x0 + (ir-1)*dx
   r = exp(r) / z
   write(2,'(g23.16)') r
enddo
write(2,'('' r*s pseudopotential in Ry'')')
read(1,*)
read(1,*) (s(i),i=1,nx)
write(2,'(g23.16)') 0.0d0
do ir = 1, nx
   r = x0 + (ir-1)*dx
   r = exp(r) / z
   write(2,'(g23.16)') r*s(ir)
enddo
write(2,'('' r*p pseudopotential in Ry'')')
read(1,*)
read(1,*) (s(i),i=1,nx)
write(2,'(g23.16)') 0.0d0
do ir = 1, nx
   r = x0 + (ir-1)*dx
   r = exp(r) / z
   write(2,'(g23.16)') r*s(ir)
enddo
write(2,'('' r*d pseudopotential in Ry'')')
read(1,*)
read(1,*) (s(i),i=1,nx)
write(2,'(g23.16)') 0.0d0
do ir = 1, nx
   r = x0 + (ir-1)*dx
   r = exp(r) / z
   write(2,'(g23.16)') r*s(ir)
enddo

end

module slaarnabm
use dsp,         only : dp
use format_utils,only : wout
use run_control, only : errstop
implicit none
integer xyzzyaaaa1
real(dp) xyzzyaaab1(3),xyzzyaaac1(3,3),xyzzyaaad1
logical use_magnetic_field
private
public use_magnetic_field,read_magnetic_field,eval_magnetic_energy
contains
subroutine read_magnetic_field
use file_utils,only : open_units
use parallel,only : am_master
use store,only : complex_wf,open_unit
implicit none
integer xyzzyaaaa2,xyzzyaaab2
real(dp) xyzzyaaac2(3)
real(dp),parameter :: xyzzyaaad2=1.d-8
character(80) char80,tmpr
if(am_master)then
call wout('External magnetic fields')
call wout('========================')
call wout('Reading magnetic field from expot.data.')
endif
if(.not.complex_wf.and.am_master)call errstop('READ_MAGNETIC_FIELD','C&
&OMPLEX_WF must be T if magnetic fields are present.')
xyzzyaaaa1=0
call open_units(xyzzyaaaa2,xyzzyaaab2)
if(xyzzyaaab2/=0)call errstop('READ_MAGNETIC_FIELD','Cannot find avail&
&able unit for reading expot.data.')
open(unit=xyzzyaaaa2,file='expot.data',status='old',action='read',iost&
&at=xyzzyaaab2)
if(xyzzyaaab2/=0)call errstop('READ_MAGNETIC_FIELD','Unable to open ex&
&pot.data.')
do
read(xyzzyaaaa2,'(a)',iostat=xyzzyaaab2)char80
if(xyzzyaaab2>0)then
call errstop('READ_MAGNETIC_FIELD','Error reading expot.data.')
elseif(xyzzyaaab2<0)then
exit
endif
char80=adjustl(char80)
if(trim(char80)=="START UNIFORM MAGNETIC FIELD")then
if(xyzzyaaaa1/=0.and.am_master)call errstop('READ_MAGNETIC_FIELD','Mor&
&e than one type of magnetic field is specified in expot.data.')
xyzzyaaaa1=1
if(am_master)call wout('Uniform external magnetic field.')
read(xyzzyaaaa2,*,err=666,end=666)
read(xyzzyaaaa2,*,err=666,end=666)xyzzyaaab1(1:3)
if(am_master)then
call wout(' Magnetic vector potential is A = A0 + A1.r, where:')
write(tmpr,'(a,es20.12,",",es20.12,",",es20.12,")")')'  A0 = (',xyzzya&
&aab1(1:3)
call wout(tmpr)
endif
read(xyzzyaaaa2,*,err=666,end=666)
read(xyzzyaaaa2,*,err=666,end=666)xyzzyaaac1(1,1:3)
read(xyzzyaaaa2,*,err=666,end=666)xyzzyaaac1(2,1:3)
read(xyzzyaaaa2,*,err=666,end=666)xyzzyaaac1(3,1:3)
if(am_master)then
write(tmpr,'(a,es20.12," ",es20.12," ",es20.12,")")')'  A1 = (',xyzzya&
&aac1(1,1:3)
call wout(tmpr)
write(tmpr,'(a,es20.12," ",es20.12," ",es20.12,")")')'       (',xyzzya&
&aac1(2,1:3)
call wout(tmpr)
write(tmpr,'(a,es20.12," ",es20.12," ",es20.12,")")')'       (',xyzzya&
&aac1(3,1:3)
call wout(tmpr)
xyzzyaaac2(1)=xyzzyaaac1(3,2)-xyzzyaaac1(2,3)
xyzzyaaac2(2)=xyzzyaaac1(1,3)-xyzzyaaac1(3,1)
xyzzyaaac2(3)=xyzzyaaac1(2,1)-xyzzyaaac1(1,2)
call wout(' Uniform magnetic field B = curl(A) is:')
write(tmpr,'(a,es20.12,",",es20.12,",",es20.12,")")')'  B = (',xyzzyaa&
&ac2(1:3)
call wout(tmpr)
xyzzyaaad1=xyzzyaaac1(1,1)+xyzzyaaac1(2,2)+xyzzyaaac1(3,3)
if(abs(xyzzyaaad1)<xyzzyaaad2*sum(abs(xyzzyaaac1))/9.d0)then
call wout(' Coulomb gauge has been used.')
else
call wout(' Unknown gauge.')
endif
call wout()
endif
read(xyzzyaaaa2,'(a)',iostat=xyzzyaaab2)char80
if(xyzzyaaab2/=0)char80=''
if(trim(adjustl(char80))/='END UNIFORM MAGNETIC FIELD'.and.am_master)c&
&all errstop('READ_MAGNETIC_FIELD','Was expecting to find "END UNIFORM&
& MAGNETIC FIELD".')
endif
enddo
close(xyzzyaaaa2)
open_unit(xyzzyaaaa2)=.false.
if(xyzzyaaaa1==0.and.am_master)call errstop('READ_MAGNETIC_FIELD','Hav&
&e not found a magnetic field in expot.data.  At present only the UNIF&
&ORM MAGNETIC FIELD block is allowed.')
return
666 call errstop('READ_MAGNETIC_FIELD','Error reading expot.data.')
end subroutine read_magnetic_field
subroutine eval_magnetic_energy(rvec,ispin,drift_r,drift_i,emagnetic_r&
&,emagnetic_i)
use store,only : pcharge,inv_pmass,electron_system
implicit none
integer,intent(in) :: ispin
real(dp),intent(in) :: rvec(3),drift_r(3),drift_i(3)
real(dp),intent(out) :: emagnetic_r,emagnetic_i
real(dp) xyzzyaaaa3(3),xyzzyaaab3
call xyzzyaaae1(rvec,xyzzyaaaa3,xyzzyaaab3)
if(electron_system)then
emagnetic_r=0.5d0*dot_product(xyzzyaaaa3,xyzzyaaaa3)+dot_product(xyzzy&
&aaaa3,drift_i)
emagnetic_i=-0.5d0*xyzzyaaab3-dot_product(xyzzyaaaa3,drift_r)
else
emagnetic_r=(0.5d0*pcharge(ispin)*dot_product(xyzzyaaaa3,xyzzyaaaa3)-d&
&ot_product(xyzzyaaaa3,drift_i))*pcharge(ispin)*inv_pmass(ispin)
emagnetic_i=(0.5d0*xyzzyaaab3+dot_product(xyzzyaaaa3,drift_r))*pcharge&
&(ispin)*inv_pmass(ispin)
endif
end subroutine eval_magnetic_energy
subroutine xyzzyaaae1(rvec,a,div_a)
implicit none
real(dp),intent(in) :: rvec(3)
real(dp),intent(out) :: a(3),div_a
if(xyzzyaaaa1==1)then
a=xyzzyaaab1+matmul(xyzzyaaac1,rvec)
div_a=xyzzyaaad1
else
call errstop('EVAL_MAGNETIC_VECTOR_POT','Only uniform magnetic fields &
&are currently implemented.')
endif
end subroutine xyzzyaaae1
end module slaarnabm

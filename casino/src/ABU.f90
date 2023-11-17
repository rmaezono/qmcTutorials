module slaarnabu
use dsp,         only : dp
use file_utils,  only : open_units
use format_utils,only : wout,i2s,display_param
use slaarnabt,   only : dcopy
use parallel,    only : am_master
use slaarnaca,       only : have_ae,zion
use run_control, only : errstop,errwarn,check_alloc
use store,       only : levels_ssingles,no_ssingles,which_ssingle
implicit none
private
public read_orbmods,write_orbmods,eval_orbmods,setup_orbmod_params,fin&
&ish_orbmod_params,get_orbmod_params,put_orbmod_params,setup_orbmods
integer xyzzyaaaa1
integer xyzzyaaab1,xyzzyaaac1
integer xyzzyaaad1
integer xyzzyaaae1
integer,allocatable :: xyzzyaaaf1(:)
integer,allocatable :: xyzzyaaag1(:),xyzzyaaah1(:)
real(dp),allocatable :: xyzzyaaai1(:,:),xyzzyaaaj1(:,:),xyzzyaaak1(:,:&
&,:)
integer,allocatable :: xyzzyaaal1(:,:),xyzzyaaam1(:,:),xyzzyaaan1(:,:,&
&:)
integer,allocatable :: xyzzyaaao1(:)
character(80) :: xyzzyaaap1
real(dp),allocatable :: xyzzyaaaq1(:,:,:),xyzzyaaar1(:,:,:),xyzzyaaas1&
&(:,:,:,:)
integer xyzzyaaat1,xyzzyaaau1
contains
subroutine read_orbmods
implicit none
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2,xyzzyaa&
&af2(3),i,j,xyzzyaaag2,total_no_params
integer,allocatable :: xyzzyaaah2(:,:,:)
real(dp),allocatable :: xyzzyaaai2(:,:,:)
character(80) :: xyzzyaaaj2
call open_units(xyzzyaaaa1,xyzzyaaaa2)
if(xyzzyaaaa2/=0)call errstop('READ_ORBMODS','Cannot find free IO unit&
&.')
open(unit=xyzzyaaaa1,file='correlation.data',status='old',iostat=xyzzy&
&aaaa2)
if(xyzzyaaaa2/=0)call errstop('READ_ORBMODS','Problem opening correlat&
&ion.data')
if(am_master)then
call wout()
call wout('Atomic orbital modification functions')
call wout('=====================================')
call wout('Reading orbital modifications from correlation.data file.')
call wout()
endif
do
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaaa2)xyzzyaaaj2
if(xyzzyaaaa2>0.and.am_master)call errstop('READ_ORBMODS','Problem rea&
&ding correlation.data. Check this file.')
if(xyzzyaaaa2<0.and.am_master)call errstop('READ_ORBMODS','Could not f&
&ind "START ORBMODS" in correlation.data.')
if(trim(adjustl(xyzzyaaaj2))=='START ORBMODS')exit
enddo
read(xyzzyaaaa1,*,err=666,end=666)
read(xyzzyaaaa1,'(a)',err=666,end=666)xyzzyaaap1
if(am_master)call wout('Title: '//trim(adjustl(xyzzyaaap1)))
read(xyzzyaaaa1,*,err=666,end=666)
read(xyzzyaaaa1,*,err=666,end=666)xyzzyaaab1
if(am_master)then
call wout('Spin dependence: '//trim(i2s(xyzzyaaab1)))
if(xyzzyaaab1<0.or.xyzzyaaab1>levels_ssingles)call errstop('READ_ORBMO&
&DS','Spin-dep should be 0 -- '//trim(i2s(levels_ssingles))//'.')
endif
xyzzyaaac1=no_ssingles(xyzzyaaab1)
read(xyzzyaaaa1,*,err=666,end=666)
read(xyzzyaaaa1,*,err=666,end=666)xyzzyaaad1
if(am_master)then
call wout('Number of modification functions: '//trim(i2s(xyzzyaaad1)))
if(xyzzyaaad1<0)call errstop('READ_ORBMODS','No. of mod fns is <0.')
endif
allocate(xyzzyaaag1(xyzzyaaad1),xyzzyaaah1(xyzzyaaad1),xyzzyaaai1(xyzz&
&yaaad1,xyzzyaaac1),xyzzyaaal1(xyzzyaaad1,xyzzyaaac1),xyzzyaaaj1(xyzzy&
&aaad1,xyzzyaaac1),xyzzyaaam1(xyzzyaaad1,xyzzyaaac1),xyzzyaaaf1(xyzzya&
&aad1),stat=xyzzyaaab2)
call check_alloc(xyzzyaaab2,'READ_ORBMODS','1')
xyzzyaaai1=1.d0
xyzzyaaal1=1
xyzzyaaaj1=0.d0
xyzzyaaam1=1
total_no_params=0
do xyzzyaaac2=1,xyzzyaaad1
do
read(xyzzyaaaa1,'(a)',err=666,end=666)xyzzyaaaj2
if(trim(adjustl(xyzzyaaaj2))=='START MODFN '//trim(i2s(xyzzyaaac2)))ex&
&it
if(trim(adjustl(xyzzyaaaj2))/=''.and.am_master)call errstop('READ_ORBM&
&ODS','Was expecting to find "START MODFN '//trim(i2s(xyzzyaaac2))//'"&
&.')
enddo
if(am_master)call wout('Orbital modification no. '//trim(i2s(xyzzyaaac&
&2))//':')
read(xyzzyaaaa1,*,err=666,end=666)
read(xyzzyaaaa1,*,err=666,end=666)xyzzyaaag1(xyzzyaaac2),xyzzyaaah1(xy&
&zzyaaac2)
if(am_master)then
call wout(' n quantum number                  :  '//trim(i2s(xyzzyaaag&
&1(xyzzyaaac2))))
call wout(' l quantum number                  :  '//trim(i2s(xyzzyaaah&
&1(xyzzyaaac2))))
if(xyzzyaaag1(xyzzyaaac2)<1)call errstop('READ_ORBMODS','Should have n&
&>0.')
if(xyzzyaaah1(xyzzyaaac2)<0.or.xyzzyaaah1(xyzzyaaac2)>=xyzzyaaag1(xyzz&
&yaaac2))call errstop('READ_ORBMODS','Should have 0<=l<n.')
endif
read(xyzzyaaaa1,*,err=666,end=666)
read(xyzzyaaaa1,*,err=666,end=666)xyzzyaaaf1(xyzzyaaac2)
if(am_master)then
call wout(' Expansion order                   :  '//trim(i2s(xyzzyaaaf&
&1(xyzzyaaac2))))
if(xyzzyaaaf1(xyzzyaaac2)<1)call errstop('READ_ORBMODS','Expansion ord&
&er should be at least 1.')
endif
if(xyzzyaaac2==1)then
allocate(xyzzyaaak1(0:xyzzyaaaf1(xyzzyaaac2),xyzzyaaad1,xyzzyaaac1),xy&
&zzyaaan1(0:xyzzyaaaf1(xyzzyaaac2),xyzzyaaad1,xyzzyaaac1),stat=xyzzyaa&
&ab2)
call check_alloc(xyzzyaaab2,'READ_ORBMODS','2')
xyzzyaaak1=0.d0
xyzzyaaan1=1
elseif(xyzzyaaaf1(xyzzyaaac2)>maxval(xyzzyaaaf1(1:xyzzyaaac2-1)))then
xyzzyaaaf2=shape(xyzzyaaak1)
allocate(xyzzyaaai2(0:xyzzyaaaf2(1)-1,xyzzyaaaf2(2),xyzzyaaaf2(3)),xyz&
&zyaaah2(0:xyzzyaaaf2(1)-1,xyzzyaaaf2(2),xyzzyaaaf2(3)),stat=xyzzyaaab&
&2)
call check_alloc(xyzzyaaab2,'READ_ORBMODS','3')
xyzzyaaai2=xyzzyaaak1
xyzzyaaah2=xyzzyaaan1
deallocate(xyzzyaaak1,xyzzyaaan1)
allocate(xyzzyaaak1(0:xyzzyaaaf1(xyzzyaaac2),xyzzyaaad1,xyzzyaaac1),xy&
&zzyaaan1(0:xyzzyaaaf1(xyzzyaaac2),xyzzyaaad1,xyzzyaaac1),stat=xyzzyaa&
&ab2)
call check_alloc(xyzzyaaab2,'READ_ORBMODS','4')
xyzzyaaak1=0.d0
xyzzyaaan1=1
xyzzyaaak1(0:xyzzyaaaf2(1)-1,xyzzyaaaf2(2),xyzzyaaaf2(3))=xyzzyaaai2(0&
&:xyzzyaaaf2(1)-1,xyzzyaaaf2(2),xyzzyaaaf2(3))
xyzzyaaan1(0:xyzzyaaaf2(1)-1,xyzzyaaaf2(2),xyzzyaaaf2(3))=xyzzyaaah2(0&
&:xyzzyaaaf2(1)-1,xyzzyaaaf2(2),xyzzyaaaf2(3))
deallocate(xyzzyaaai2,xyzzyaaah2)
endif
read(xyzzyaaaa1,*,err=666,end=666)
param_loop: do xyzzyaaad2=1,xyzzyaaac1
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)xyzzyaaai1(xyzzyaaac2,xyzzyaaad2),&
&xyzzyaaal1(xyzzyaaac2,xyzzyaaad2)
if(xyzzyaaaa2/=0)then
backspace xyzzyaaaa1
if(am_master)call wout(' Not all coefficients supplied: rest given def&
&ault values.')
xyzzyaaai1(xyzzyaaac2,xyzzyaaad2)=1.d0
xyzzyaaal1(xyzzyaaac2,xyzzyaaad2)=1
exit param_loop
endif
if(am_master)then
call display_param(xyzzyaaai1(xyzzyaaac2,xyzzyaaad2),xyzzyaaal1(xyzzya&
&aac2,xyzzyaaad2),'A_'//trim(i2s(xyzzyaaac2))//','//trim(i2s(xyzzyaaad&
&2)),indent_in=2)
if(xyzzyaaal1(xyzzyaaac2,xyzzyaaad2)/=0.and.xyzzyaaal1(xyzzyaaac2,xyzz&
&yaaad2)/=1)call errstop('READ_ORBMODS','Optimizable flag should be 0 &
&(fixed) or 1 (optimizable).')
if(xyzzyaaai1(xyzzyaaac2,xyzzyaaad2)<=0.d0)call errstop('READ_ORBMODS'&
&,'The A parameter should be strictly positive.')
endif
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)xyzzyaaaj1(xyzzyaaac2,xyzzyaaad2),&
&xyzzyaaam1(xyzzyaaac2,xyzzyaaad2)
if(xyzzyaaaa2/=0)then
backspace xyzzyaaaa1
if(am_master)call wout(' Not all coefficients supplied: rest given def&
&ault values.')
xyzzyaaaj1(xyzzyaaac2,xyzzyaaad2)=0.d0
xyzzyaaam1(xyzzyaaac2,xyzzyaaad2)=1
exit param_loop
endif
if(am_master)then
call display_param(xyzzyaaaj1(xyzzyaaac2,xyzzyaaad2),xyzzyaaam1(xyzzya&
&aac2,xyzzyaaad2),'B_'//trim(i2s(xyzzyaaac2))//','//trim(i2s(xyzzyaaad&
&2)),indent_in=2)
if(xyzzyaaam1(xyzzyaaac2,xyzzyaaad2)/=0.and.xyzzyaaam1(xyzzyaaac2,xyzz&
&yaaad2)/=1)call errstop('READ_ORBMODS','Optimizable flag should be 0 &
&(fixed) or 1 (optimizable).')
if(xyzzyaaaj1(xyzzyaaac2,xyzzyaaad2)<0.d0)call errstop('READ_ORBMODS',&
&'The B parameter should be non-negative.')
endif
do xyzzyaaae2=0,xyzzyaaaf1(xyzzyaaac2)
if(xyzzyaaae2/=1.or.xyzzyaaah1(xyzzyaaac2)/=0)then
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)xyzzyaaak1(xyzzyaaae2,xyzzyaaac2,x&
&yzzyaaad2),xyzzyaaan1(xyzzyaaae2,xyzzyaaac2,xyzzyaaad2)
if(xyzzyaaaa2/=0)then
backspace xyzzyaaaa1
if(am_master)call wout(' Not all coefficients supplied: rest given def&
&ault values.')
xyzzyaaak1(xyzzyaaae2,xyzzyaaac2,xyzzyaaad2)=0.d0
xyzzyaaan1(xyzzyaaae2,xyzzyaaac2,xyzzyaaad2)=1
exit param_loop
endif
if(am_master)then
call display_param(xyzzyaaak1(xyzzyaaae2,xyzzyaaac2,xyzzyaaad2),xyzzya&
&aan1(xyzzyaaae2,xyzzyaaac2,xyzzyaaad2),'c_'//trim(i2s(xyzzyaaae2))//'&
&,'//trim(i2s(xyzzyaaac2))//','//trim(i2s(xyzzyaaad2)),indent_in=2)
if(xyzzyaaan1(xyzzyaaae2,xyzzyaaac2,xyzzyaaad2)/=0.and.xyzzyaaan1(xyzz&
&yaaae2,xyzzyaaac2,xyzzyaaad2)/=1)call errstop('READ_ORBMODS','Optimiz&
&able flag should be 0 (fixed) or 1 (optimizable).')
endif
endif
enddo
enddo param_loop
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaaa2)xyzzyaaaj2
if(xyzzyaaaa2/=0.and.am_master)call errstop('READ_ORBMODS','Was expect&
&ing to find "END MODFN '//trim(i2s(xyzzyaaac2))//'".')
if(trim(adjustl(xyzzyaaaj2))/='END MODFN '//trim(i2s(xyzzyaaac2)))call&
& errstop('READ_ORBMODS','Was expecting to find "END MODFN ' //trim(i2&
&s(xyzzyaaac2))//'".')
if(xyzzyaaah1(xyzzyaaac2)==0)then
xyzzyaaag2=xyzzyaaac1*(2+xyzzyaaaf1(xyzzyaaac2))
else
xyzzyaaag2=xyzzyaaac1*(3+xyzzyaaaf1(xyzzyaaac2))
endif
total_no_params=total_no_params+xyzzyaaag2
if(am_master)call wout(' There are '//trim(i2s(xyzzyaaag2))//' paramet&
&ers in this mod. fn., including A and B.')
enddo
do
read(xyzzyaaaa1,'(a)',err=666,end=666)xyzzyaaaj2
if(trim(adjustl(xyzzyaaaj2))=='END ORBMODS')exit
if(trim(adjustl(xyzzyaaaj2))/='')call errstop('READ_ORBMODS','Was expe&
&cting to find "END ORBMODS".')
enddo
do
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaaa2)xyzzyaaaj2
if(xyzzyaaaa2<0)exit
if(xyzzyaaaa2>0)call errstop('READ_ORBMODS','Problem reading correlati&
&on.data. Check this file.')
if(trim(adjustl(xyzzyaaaj2))=='START ORBMODS')call errstop('READ_ORBMO&
&DS','There seems to be more than one set of atomic orbital modificati&
&on data in correlation.data.')
enddo
close(xyzzyaaaa1)
if(am_master)then
call wout('Finished reading atomic orbital modifications from correlat&
&ion.data.')
call wout('In total, there are '//trim(i2s(total_no_params))//' parame&
&ters in the modification functions.')
call wout()
endif
if(am_master)then
do i=1,xyzzyaaad1-1
do j=i+1,xyzzyaaad1
if(xyzzyaaag1(i)==xyzzyaaag1(j).and.xyzzyaaah1(i)==xyzzyaaah1(j))call &
&errstop('READ_ORBMODS','Modification functions '//trim(i2s(i))//' and&
& '//trim(i2s(j))//' have the same quantum numbers.')
enddo
enddo
endif
return
666 call errstop('READ_ORBMODS','Problem reading orbital modifications&
& in correlation.data. Check this file.')
end subroutine read_orbmods
subroutine setup_orbmods(aw_norb_rad,aw_orbdat_rad)
implicit none
integer,intent(in) :: aw_norb_rad,aw_orbdat_rad(3,aw_norb_rad)
call xyzzyaaav1
call xyzzyaaaw1(aw_norb_rad,aw_orbdat_rad)
end subroutine setup_orbmods
subroutine write_orbmods(correlation_name)
implicit none
character(20),intent(in) :: correlation_name
integer xyzzyaaaa4,xyzzyaaab4,xyzzyaaac4,xyzzyaaad4
logical xyzzyaaae4
character(64) blurb
if(am_master)then
inquire(file=trim(correlation_name),exist=xyzzyaaae4)
if(xyzzyaaae4)then
open(unit=xyzzyaaaa1,file=trim(correlation_name),status='old',position&
&='append',iostat=xyzzyaaad4)
else
open(unit=xyzzyaaaa1,file=trim(correlation_name),status='replace',iost&
&at=xyzzyaaad4)
endif
if(xyzzyaaad4/=0)call errstop('WRITE_ORBMODS','Problem opening '//trim&
&(correlation_name)//'.')
write(xyzzyaaaa1,*)'START ORBMODS'
write(xyzzyaaaa1,*)'Title'
write(xyzzyaaaa1,*)trim(xyzzyaaap1)
write(xyzzyaaaa1,*)'Spin dependence (0->u=d; 1->u/=d)'
write(xyzzyaaaa1,'(3x,a)')trim(i2s(xyzzyaaab1))
write(xyzzyaaaa1,*)'Number of modification functions'
write(xyzzyaaaa1,'(3x,a)')trim(i2s(xyzzyaaad1))
do xyzzyaaaa4=1,xyzzyaaad1
write(xyzzyaaaa1,*)'START MODFN ',trim(i2s(xyzzyaaaa4))
write(xyzzyaaaa1,*)'Quantum numbers n and l'
write(xyzzyaaaa1,'(3x,a,2x,a)')trim(i2s(xyzzyaaag1(xyzzyaaaa4))),trim(&
&i2s(xyzzyaaah1(xyzzyaaaa4)))
write(xyzzyaaaa1,*)'Expansion order'
write(xyzzyaaaa1,'(3x,a)')trim(i2s(xyzzyaaaf1(xyzzyaaaa4)))
write(xyzzyaaaa1,*)'Parameters in modification function     ;  Optimiz&
&able (0=NO; 1=YES)'
do xyzzyaaab4=1,xyzzyaaac1
blurb='      ! A_'//trim(i2s(xyzzyaaab4))
write(xyzzyaaaa1,*)xyzzyaaai1(xyzzyaaaa4,xyzzyaaab4),xyzzyaaal1(xyzzya&
&aaa4,xyzzyaaab4),trim(blurb)
blurb='      ! B_'//trim(i2s(xyzzyaaab4))
write(xyzzyaaaa1,*)xyzzyaaaj1(xyzzyaaaa4,xyzzyaaab4),xyzzyaaam1(xyzzya&
&aaa4,xyzzyaaab4),trim(blurb)
do xyzzyaaac4=0,xyzzyaaaf1(xyzzyaaaa4)
if(xyzzyaaac4/=1.or.xyzzyaaah1(xyzzyaaaa4)/=0)then
blurb='      ! c_'//trim(i2s(xyzzyaaac4))//','//trim(i2s(xyzzyaaab4))
write(xyzzyaaaa1,*)xyzzyaaak1(xyzzyaaac4,xyzzyaaaa4,xyzzyaaab4),xyzzya&
&aan1(xyzzyaaac4,xyzzyaaaa4,xyzzyaaab4),trim(blurb)
endif
enddo
enddo
write(xyzzyaaaa1,*)'END MODFN '//trim(i2s(xyzzyaaaa4))
enddo
write(xyzzyaaaa1,*)'END ORBMODS'
write(xyzzyaaaa1,*)
close(xyzzyaaaa1)
endif
end subroutine write_orbmods
subroutine xyzzyaaav1
implicit none
integer xyzzyaaaa5,xyzzyaaab5
if(have_ae)then
do xyzzyaaaa5=1,xyzzyaaad1
if(xyzzyaaah1(xyzzyaaaa5)==0)then
do xyzzyaaab5=1,xyzzyaaac1
xyzzyaaak1(1,xyzzyaaaa5,xyzzyaaab5)=-zion(1)*xyzzyaaak1(0,xyzzyaaaa5,x&
&yzzyaaab5)
xyzzyaaan1(1,xyzzyaaaa5,xyzzyaaab5)=-1
enddo
endif
enddo
else
do xyzzyaaaa5=1,xyzzyaaad1
if(xyzzyaaah1(xyzzyaaaa5)==0)then
do xyzzyaaab5=1,xyzzyaaac1
xyzzyaaak1(1,xyzzyaaaa5,xyzzyaaab5)=0.d0
xyzzyaaan1(1,xyzzyaaaa5,xyzzyaaab5)=-1
enddo
endif
enddo
endif
end subroutine xyzzyaaav1
subroutine xyzzyaaaw1(aw_norb_rad,aw_orbdat_rad)
implicit none
integer,intent(in) :: aw_norb_rad,aw_orbdat_rad(3,aw_norb_rad)
integer xyzzyaaaa6,xyzzyaaab6,xyzzyaaac6
logical xyzzyaaad6
if(.not.allocated(xyzzyaaao1))then
allocate(xyzzyaaao1(aw_norb_rad),stat=xyzzyaaac6)
call check_alloc(xyzzyaaac6,'CONSTRUCT_INDICES','modfn_index')
endif
do xyzzyaaaa6=1,aw_norb_rad
xyzzyaaad6=.false.
do xyzzyaaab6=1,xyzzyaaad1
if(xyzzyaaag1(xyzzyaaab6)==aw_orbdat_rad(2,xyzzyaaaa6).and.xyzzyaaah1(&
&xyzzyaaab6)==aw_orbdat_rad(3,xyzzyaaaa6))then
if(xyzzyaaad6.and.am_master)call errstop('CONSTRUCT_INDICES','Mod-fn a&
&mbiguity.')
xyzzyaaao1(xyzzyaaaa6)=xyzzyaaab6
xyzzyaaad6=.true.
endif
enddo
if(.not.xyzzyaaad6)then
xyzzyaaao1(xyzzyaaaa6)=0
call errwarn('CONSTRUCT_INDICES','no mod. fn. found for orb. '//trim(i&
&2s(xyzzyaaaa6))//'. This orbital has quantum numbers n='//trim(i2s(aw&
&_orbdat_rad(2,xyzzyaaaa6)))//', l='//trim(i2s(aw_orbdat_rad(3,xyzzyaa&
&aa6)))//'.')
endif
enddo
if(am_master)then
do xyzzyaaab6=1,xyzzyaaad1
xyzzyaaad6=.false.
orb_loop : do xyzzyaaaa6=1,aw_norb_rad
if(xyzzyaaao1(xyzzyaaaa6)==xyzzyaaab6)then
xyzzyaaad6=.true.
exit orb_loop
endif
enddo orb_loop
if(.not.xyzzyaaad6)call errwarn('CONSTRUCT_INDICES','mod. fn. '//trim(&
&i2s(xyzzyaaab6))//' appears in correlation.data, but is not used.')
enddo
endif
end subroutine xyzzyaaaw1
subroutine eval_orbmods(r,jspin,iorb_rad,w_r_to_l,d_w_r_to_l,d2_w_r_to&
&_l)
implicit none
integer,intent(in) :: jspin,iorb_rad
real(dp),intent(in) :: r
real(dp),intent(out) :: w_r_to_l,d_w_r_to_l,d2_w_r_to_l
integer xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7,xyzzyaaad7,xyzzyaaae7
real(dp) xyzzyaaaf7,xyzzyaaag7,xyzzyaaah7,xyzzyaaai7,xyzzyaaaj7,xyzzya&
&aak7,xyzzyaaal7,xyzzyaaam7,xyzzyaaan7,xyzzyaaao7,xyzzyaaap7,xyzzyaaaq&
&7,xyzzyaaar7,xyzzyaaas7,xyzzyaaat7
xyzzyaaaa7=xyzzyaaao1(iorb_rad)
if(xyzzyaaaa7>0)then
xyzzyaaab7=which_ssingle(jspin,xyzzyaaab1)
xyzzyaaac7=xyzzyaaaf1(xyzzyaaaa7)
xyzzyaaaf7=exp(-xyzzyaaai1(xyzzyaaaa7,xyzzyaaab7)*r*r/(1.d0+xyzzyaaaj1&
&(xyzzyaaaa7,xyzzyaaab7)*r))
xyzzyaaag7=-xyzzyaaaf7*xyzzyaaai1(xyzzyaaaa7,xyzzyaaab7)*r*(2.d0+xyzzy&
&aaaj1(xyzzyaaaa7,xyzzyaaab7)*r)/((1.d0+xyzzyaaaj1(xyzzyaaaa7,xyzzyaaa&
&b7)*r)*(1.d0+xyzzyaaaj1(xyzzyaaaa7,xyzzyaaab7)*r))
xyzzyaaah7=xyzzyaaaf7*xyzzyaaai1(xyzzyaaaa7,xyzzyaaab7)*(-2.d0+r*(-2.d&
&0*xyzzyaaaj1(xyzzyaaaa7,xyzzyaaab7)+r*xyzzyaaai1(xyzzyaaaa7,xyzzyaaab&
&7)*(4.d0+r*xyzzyaaaj1(xyzzyaaaa7,xyzzyaaab7)*(4.d0+xyzzyaaaj1(xyzzyaa&
&aa7,xyzzyaaab7)*r))))/(1.d0+xyzzyaaaj1(xyzzyaaaa7,xyzzyaaab7)*r)**4
xyzzyaaai7=xyzzyaaak1(0,xyzzyaaaa7,xyzzyaaab7)+xyzzyaaak1(1,xyzzyaaaa7&
&,xyzzyaaab7)*r
xyzzyaaaj7=xyzzyaaak1(1,xyzzyaaaa7,xyzzyaaab7)
xyzzyaaak7=0.d0
if(xyzzyaaac7>=2)then
xyzzyaaan7=1.d0
xyzzyaaam7=r
xyzzyaaal7=r*r
do xyzzyaaad7=2,xyzzyaaac7-1
xyzzyaaai7=xyzzyaaai7+xyzzyaaak1(xyzzyaaad7,xyzzyaaaa7,xyzzyaaab7)*xyz&
&zyaaal7
xyzzyaaam7=dble(xyzzyaaad7)*xyzzyaaam7
xyzzyaaaj7=xyzzyaaaj7+xyzzyaaak1(xyzzyaaad7,xyzzyaaaa7,xyzzyaaab7)*xyz&
&zyaaam7
xyzzyaaak7=xyzzyaaak7+xyzzyaaak1(xyzzyaaad7,xyzzyaaaa7,xyzzyaaab7)*dbl&
&e(xyzzyaaad7)*xyzzyaaan7
xyzzyaaan7=xyzzyaaam7
xyzzyaaam7=xyzzyaaal7
xyzzyaaal7=xyzzyaaal7*r
enddo
xyzzyaaai7=xyzzyaaai7+xyzzyaaak1(xyzzyaaad7,xyzzyaaaa7,xyzzyaaab7)*xyz&
&zyaaal7
xyzzyaaaj7=xyzzyaaaj7+xyzzyaaak1(xyzzyaaac7,xyzzyaaaa7,xyzzyaaab7)*dbl&
&e(xyzzyaaac7)*xyzzyaaam7
xyzzyaaak7=xyzzyaaak7+xyzzyaaak1(xyzzyaaac7,xyzzyaaaa7,xyzzyaaab7)*dbl&
&e(xyzzyaaac7)*xyzzyaaan7
endif
xyzzyaaao7=xyzzyaaai7*xyzzyaaaf7
xyzzyaaap7=xyzzyaaaj7*xyzzyaaaf7+xyzzyaaag7*xyzzyaaai7
xyzzyaaaq7=xyzzyaaak7*xyzzyaaaf7+2.d0*xyzzyaaaj7*xyzzyaaag7+xyzzyaaai7&
&*xyzzyaaah7
xyzzyaaae7=xyzzyaaah1(xyzzyaaaa7)
if(xyzzyaaae7==0)then
w_r_to_l=xyzzyaaao7
d_w_r_to_l=xyzzyaaap7
d2_w_r_to_l=xyzzyaaaq7
elseif(xyzzyaaae7==1)then
w_r_to_l=r*xyzzyaaao7
d_w_r_to_l=r*xyzzyaaap7+xyzzyaaao7
d2_w_r_to_l=r*xyzzyaaaq7+2.d0*xyzzyaaap7
elseif(xyzzyaaae7==2)then
w_r_to_l=r*r*xyzzyaaao7
d_w_r_to_l=r*(2.d0*xyzzyaaao7+r*xyzzyaaap7)
d2_w_r_to_l=2.d0*xyzzyaaao7+r*(4.d0*xyzzyaaap7+r*xyzzyaaaq7)
else
xyzzyaaar7=r**(xyzzyaaae7-2)
xyzzyaaas7=xyzzyaaar7*r
xyzzyaaat7=xyzzyaaas7*r
w_r_to_l=xyzzyaaat7*xyzzyaaao7
d_w_r_to_l=xyzzyaaat7*xyzzyaaap7+dble(xyzzyaaae7)*xyzzyaaas7*xyzzyaaao&
&7
d2_w_r_to_l=xyzzyaaat7*xyzzyaaaq7+2.d0*dble(xyzzyaaae7)*xyzzyaaas7*xyz&
&zyaaap7 +dble(xyzzyaaae7)*(dble(xyzzyaaae7)-1.d0)*xyzzyaaar7*xyzzyaaa&
&o7
endif
else
w_r_to_l=0.d0
d_w_r_to_l=0.d0
d2_w_r_to_l=0.d0
endif
end subroutine eval_orbmods
subroutine setup_orbmod_params(nparam)
implicit none
integer,intent(inout) :: nparam
integer xyzzyaaaa8,xyzzyaaab8,xyzzyaaac8
nparam=0
do xyzzyaaaa8=1,xyzzyaaac1
do xyzzyaaab8=1,xyzzyaaad1
if(xyzzyaaal1(xyzzyaaab8,xyzzyaaaa8)==1)nparam=nparam+1
if(xyzzyaaam1(xyzzyaaab8,xyzzyaaaa8)==1)nparam=nparam+1
do xyzzyaaac8=0,xyzzyaaaf1(xyzzyaaab8)
if(xyzzyaaan1(xyzzyaaac8,xyzzyaaab8,xyzzyaaaa8)==1)nparam=nparam+1
enddo
enddo
enddo
xyzzyaaae1=nparam
call xyzzyaaax1
end subroutine setup_orbmod_params
subroutine finish_orbmod_params
implicit none
call xyzzyaaay1
end subroutine finish_orbmod_params
subroutine xyzzyaaax1
implicit none
integer xyzzyaaaa10
allocate(xyzzyaaaq1(xyzzyaaad1,xyzzyaaac1,0:xyzzyaaae1),xyzzyaaar1(xyz&
&zyaaad1,xyzzyaaac1,0:xyzzyaaae1),xyzzyaaas1(0:maxval(xyzzyaaaf1),xyzz&
&yaaad1,xyzzyaaac1,0:xyzzyaaae1),stat=xyzzyaaaa10)
call check_alloc(xyzzyaaaa10,'SETUP_ORBMOD_PBUFFER','')
xyzzyaaaq1=0.d0
xyzzyaaar1=0.d0
xyzzyaaas1=0.d0
xyzzyaaat1=xyzzyaaad1*xyzzyaaac1
xyzzyaaau1=(maxval(xyzzyaaaf1)+1)*xyzzyaaad1*xyzzyaaac1
end subroutine xyzzyaaax1
subroutine xyzzyaaay1
implicit none
deallocate(xyzzyaaaq1,xyzzyaaar1,xyzzyaaas1)
end subroutine xyzzyaaay1
subroutine get_orbmod_params(params,has_lolim,lolim,has_hilim,hilim,is&
&_shallow,is_redundant,is_linear,is_loglinear,has_aderiv,affect_map,la&
&bel)
implicit none
real(dp),intent(inout) :: params(xyzzyaaae1),lolim(xyzzyaaae1),hilim(x&
&yzzyaaae1)
logical,intent(inout) :: has_lolim(xyzzyaaae1),has_hilim(xyzzyaaae1),i&
&s_shallow(xyzzyaaae1),is_redundant(xyzzyaaae1),is_linear(xyzzyaaae1),&
&is_loglinear(xyzzyaaae1),has_aderiv(xyzzyaaae1),affect_map(xyzzyaaae1&
&,xyzzyaaae1)
character(2),intent(inout) :: label(xyzzyaaae1)
integer xyzzyaaaa12,xyzzyaaab12,xyzzyaaac12,xyzzyaaad12
has_lolim=.false.
lolim=0.d0
has_hilim=.false.
hilim=0.d0
is_shallow=.false.
is_redundant=.false.
is_linear=.false.
is_loglinear=.false.
has_aderiv=.false.
affect_map=.false.
do xyzzyaaaa12=1,xyzzyaaae1
affect_map(xyzzyaaaa12,xyzzyaaaa12)=.true.
enddo
label='OP'
xyzzyaaaa12=0
do xyzzyaaab12=1,xyzzyaaac1
do xyzzyaaac12=1,xyzzyaaad1
if(xyzzyaaal1(xyzzyaaac12,xyzzyaaab12)==1)then
xyzzyaaaa12=xyzzyaaaa12+1
params(xyzzyaaaa12)=xyzzyaaai1(xyzzyaaac12,xyzzyaaab12)
has_lolim(xyzzyaaaa12)=.true.
lolim(xyzzyaaaa12)=0.d0
endif
if(xyzzyaaam1(xyzzyaaac12,xyzzyaaab12)==1)then
xyzzyaaaa12=xyzzyaaaa12+1
params(xyzzyaaaa12)=xyzzyaaaj1(xyzzyaaac12,xyzzyaaab12)
has_lolim(xyzzyaaaa12)=.true.
lolim(xyzzyaaaa12)=0.d0
endif
do xyzzyaaad12=0,xyzzyaaaf1(xyzzyaaac12)
if(xyzzyaaan1(xyzzyaaad12,xyzzyaaac12,xyzzyaaab12)==1)then
xyzzyaaaa12=xyzzyaaaa12+1
params(xyzzyaaaa12)=xyzzyaaak1(xyzzyaaad12,xyzzyaaac12,xyzzyaaab12)
endif
enddo
enddo
enddo
end subroutine get_orbmod_params
subroutine put_orbmod_params(params,ignore,iparam_buffer,prestore,bad_&
&params)
implicit none
integer,intent(in) :: iparam_buffer
real(dp),intent(inout) :: params(xyzzyaaae1)
logical,intent(in) :: ignore(xyzzyaaae1),prestore
logical,intent(out) :: bad_params
integer xyzzyaaaa13,xyzzyaaab13,xyzzyaaac13,xyzzyaaad13
bad_params=.false.
if(prestore)then
call xyzzyaaba1(iparam_buffer)
return
endif
xyzzyaaaa13=0
do xyzzyaaab13=1,xyzzyaaac1
do xyzzyaaad13=1,xyzzyaaad1
if(xyzzyaaal1(xyzzyaaad13,xyzzyaaab13)==1)then
xyzzyaaaa13=xyzzyaaaa13+1
if(.not.ignore(xyzzyaaaa13))xyzzyaaai1(xyzzyaaad13,xyzzyaaab13)=params&
&(xyzzyaaaa13)
endif
if(xyzzyaaam1(xyzzyaaad13,xyzzyaaab13)==1)then
xyzzyaaaa13=xyzzyaaaa13+1
if(.not.ignore(xyzzyaaaa13))xyzzyaaaj1(xyzzyaaad13,xyzzyaaab13)=params&
&(xyzzyaaaa13)
endif
do xyzzyaaac13=0,xyzzyaaaf1(xyzzyaaad13)
if(xyzzyaaan1(xyzzyaaac13,xyzzyaaad13,xyzzyaaab13)==1)then
xyzzyaaaa13=xyzzyaaaa13+1
if(.not.ignore(xyzzyaaaa13))xyzzyaaak1(xyzzyaaac13,xyzzyaaad13,xyzzyaa&
&ab13)=params(xyzzyaaaa13)
endif
enddo
enddo
enddo
call xyzzyaaav1
call xyzzyaaaz1(iparam_buffer)
end subroutine put_orbmod_params
subroutine xyzzyaaaz1(indx)
implicit none
integer,intent(in) :: indx
call dcopy(xyzzyaaat1,xyzzyaaai1(1,1),1,xyzzyaaaq1(1,1,indx),1)
call dcopy(xyzzyaaat1,xyzzyaaaj1(1,1),1,xyzzyaaar1(1,1,indx),1)
call dcopy(xyzzyaaau1,xyzzyaaak1(0,1,1),1,xyzzyaaas1(0,1,1,indx),1)
end subroutine xyzzyaaaz1
subroutine xyzzyaaba1(indx)
implicit none
integer,intent(in) :: indx
call dcopy(xyzzyaaat1,xyzzyaaaq1(1,1,indx),1,xyzzyaaai1(1,1),1)
call dcopy(xyzzyaaat1,xyzzyaaar1(1,1,indx),1,xyzzyaaaj1(1,1),1)
call dcopy(xyzzyaaau1,xyzzyaaas1(0,1,1,indx),1,xyzzyaaak1(0,1,1),1)
end subroutine xyzzyaaba1
end module slaarnabu

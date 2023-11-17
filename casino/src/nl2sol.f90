module machine_constants
use run_control, only : errstop
implicit none
integer,parameter :: dp=selected_real_kind(precision(1.d0),range(1.d0)&
&)
contains
function imdcon(k) result(ival)
implicit none
integer,intent(in) :: k
integer ival
integer :: xyzzyaaaa2(3)=(/6,0,5/)
ival=xyzzyaaaa2(k)
end function imdcon
function rmdcon(k) result(fn_val)
implicit none
integer,intent(in) :: k
real(dp) fn_val
select case (k)
case(1) 
fn_val=tiny(1._dp)
case(2) 
fn_val=sqrt(tiny(1._dp))
case(3) 
fn_val=epsilon(1._dp)
case(4) 
fn_val=sqrt(epsilon(1._dp))
case(5) 
fn_val=sqrt(huge(1._dp))
case(6) 
fn_val=huge(1._dp)
case default 
fn_val=0.d0
call errstop('RMDCON','Illegal argument to function RMDCON in NL2SOL m&
&odule.')
end select
end function rmdcon
end module machine_constants
module toms573
use machine_constants
implicit none
logical :: stop_nl2sol=.false.
contains
subroutine nl2sol(n,p,x,calcr,calcj,iv,v,uiparm,urparm,ufparm)
integer,intent(in) :: n, p
integer,intent(inout) :: iv(:)
integer,intent(inout),optional :: uiparm(:)
real(dp),intent(inout) :: x(:),v(:)
real(dp),intent(inout),optional :: urparm(:),ufparm
interface
subroutine calcr(n,p,x,nf,r,uiparm,urparm,ufparm)
use machine_constants
implicit none
integer,intent(in) :: n,p
integer,intent(inout) :: nf
integer,intent(inout),optional :: uiparm(:)
real(dp),intent(in) :: x(:)
real(dp),intent(inout),optional :: urparm(:),ufparm
real(dp),intent(out) :: r(:)
end subroutine calcr
subroutine calcj(n,p,x,nf,j,uiparm,urparm,ufparm)
use machine_constants
implicit none
integer,intent(in) :: n, p
integer,intent(inout) :: nf
integer,intent(inout),optional :: uiparm(:)
real(dp),intent(in) :: x(:)
real(dp),intent(inout),optional :: urparm(:),ufparm
real(dp),intent(out) :: j(:)
end subroutine calcj
end interface
logical xyzzyaaaa5
integer xyzzyaaab5,xyzzyaaac5,nf,xyzzyaaad5
integer,parameter :: xyzzyaaae5=6,xyzzyaaaf5=7,xyzzyaaag5=2
integer,parameter :: xyzzyaaah5=27,j=33,r=50
xyzzyaaab5=94+2*n+p*(3*p+31)/2
iv(xyzzyaaah5)=xyzzyaaab5
xyzzyaaad5=xyzzyaaab5+p
iv(r)=xyzzyaaad5
xyzzyaaac5=xyzzyaaad5+n
iv(j)=xyzzyaaac5
xyzzyaaaa5=.true.
if(iv(1)/=0.and.iv(1)/=12)goto 40
xyzzyaaaa5=.false.
iv(xyzzyaaae5)=1
iv(xyzzyaaaf5)=1
10 nf=iv(xyzzyaaae5)
call calcr(n,p,x,nf,v(xyzzyaaad5:xyzzyaaad5+n-1),uiparm,urparm,ufparm)
if(xyzzyaaaa5)goto 20
if(nf>0)goto 30
iv(1)=13
goto 60
20 if(nf<=0)iv(xyzzyaaag5)=1
goto 40
30 call calcj(n,p,x,iv(xyzzyaaaf5),v(xyzzyaaac5:xyzzyaaac5+n*p-1),uipa&
&rm,urparm,ufparm)
if(iv(xyzzyaaaf5)==0)goto 50
xyzzyaaaa5=.true.
40 call nl2itr(v(xyzzyaaab5:xyzzyaaab5+p-1),iv,v(xyzzyaaac5:xyzzyaaac5&
&+n*p-1),n,n,p,v(xyzzyaaad5:xyzzyaaad5+n-1),v(1:xyzzyaaab5-1),x)
if(iv(1)-2<0)then
goto 10
elseif(iv(1)-2==0)then
goto 30
else
goto 70
endif
50 iv(1)=15
60 call itsmry(v(xyzzyaaab5:),iv,p,v,x)
70 return
end subroutine nl2sol
subroutine nl2sno(n,p,x,calcr,iv,v,uiparm,urparm,ufparm)
integer,intent(in) :: n,p
integer,intent(inout) :: iv(:)
integer,intent(inout),optional :: uiparm(:)
real(dp),intent(inout) :: x(:),v(:)
real(dp),intent(inout),optional :: urparm(:),ufparm
real(dp) x_h(p),xyzzyaaaa8(p)
interface
subroutine calcr(n,p,x,nf,r,nl2sol_iteration,been_rst,stop_opt,x_h,bad&
&_point)
use machine_constants
implicit none
integer,intent(in) :: n,p,nl2sol_iteration,been_rst
integer,intent(inout) :: nf
real(dp),intent(in) :: x(:)
real(dp),intent(in),optional :: x_h(:)
real(dp),intent(out) :: r(:)
logical,intent(inout) :: stop_opt
logical,intent(inout),optional :: bad_point
end subroutine calcr
end interface
logical xyzzyaaab8,bad_point
integer xyzzyaaac8,xyzzyaaad8,xyzzyaaae8,xyzzyaaaf8,xyzzyaaag8,xyzzyaa&
&ah8,nf,xyzzyaaai8,xyzzyaaaj8
real(dp) h,xyzzyaaak8
real(dp),parameter :: xyzzyaaal8=1.d0,xyzzyaaam8=0.d0
integer,parameter :: xyzzyaaan8=14,xyzzyaaao8=15,xyzzyaaap8=27,xyzzyaa&
&aq8=16,xyzzyaaar8=33,xyzzyaaas8=6,xyzzyaaat8=7,r=50,xyzzyaaau8=2
integer,parameter :: xyzzyaaav8=36,xyzzyaaaw8=38
iv(31)=0
xyzzyaaad8=94+2*n+p*(3*p+31)/2
iv(xyzzyaaap8)=xyzzyaaad8
xyzzyaaaj8=xyzzyaaad8+p
iv(r)=xyzzyaaaj8
xyzzyaaaf8=xyzzyaaaj8+n
iv(xyzzyaaar8)=xyzzyaaaf8
xyzzyaaai8=xyzzyaaaf8-1
if(iv(1)==0)call dfault(iv,v)
iv(xyzzyaaao8)=-abs(iv(xyzzyaaao8))
if(iv(xyzzyaaan8)/=0.and.iv(xyzzyaaao8)==0)iv(xyzzyaaao8)=-1
xyzzyaaab8=.true.
if(iv(1)/=12)goto 80
xyzzyaaab8=.false.
iv(xyzzyaaas8)=1
iv(xyzzyaaat8)=1
if(iv(xyzzyaaaq8)>0)call vscopy(p,v(xyzzyaaad8:xyzzyaaad8+p-1),xyzzyaa&
&al8)
if(v(xyzzyaaaw8)>xyzzyaaam8)call vscopy(p,v(xyzzyaaad8:xyzzyaaad8+p-1)&
&,v(xyzzyaaaw8))
10 nf=iv(xyzzyaaas8)
bad_point=.false.
call calcr(n,p,x,nf,v(xyzzyaaaj8:xyzzyaaaj8+n-1),iv(31),iv(9),stop_nl2&
&sol,bad_point=bad_point)
if(bad_point)iv(2)=1
if(xyzzyaaab8)goto 20
if(nf>0)goto 30
iv(1)=13
goto 90
20 if(nf<=0)iv(xyzzyaaau8)=1
goto 80
30 xyzzyaaag8=xyzzyaaaf8
xyzzyaaac8=xyzzyaaad8
do xyzzyaaah8=1,p
xyzzyaaak8=x(xyzzyaaah8)
h=v(xyzzyaaav8)*max(abs(xyzzyaaak8),xyzzyaaal8/v(xyzzyaaac8))
xyzzyaaac8=xyzzyaaac8+1
xyzzyaaaa8(xyzzyaaah8)=h
x_h(xyzzyaaah8)=xyzzyaaak8+h
enddo
nf=iv(xyzzyaaat8)
call calcr(n,p,x,nf,v(xyzzyaaaf8:xyzzyaaaf8+n*p-1),iv(31),iv(9),stop_n&
&l2sol,x_h=x_h)
xyzzyaaag8=xyzzyaaaf8
do xyzzyaaah8=1,p
h=xyzzyaaaa8(xyzzyaaah8)
do xyzzyaaae8=xyzzyaaaj8,xyzzyaaai8
v(xyzzyaaag8)=(v(xyzzyaaag8)-v(xyzzyaaae8))/h
xyzzyaaag8=xyzzyaaag8+1
enddo
enddo
xyzzyaaab8=.true.
80 call nl2itr(v(xyzzyaaad8:xyzzyaaad8+p-1),iv,v(xyzzyaaaf8:xyzzyaaaf8&
&+n*p-1),n,n,p,v(xyzzyaaaj8:xyzzyaaaj8+n-1),v(1:xyzzyaaad8-1),x)
if(iv(1)-2<0)then
goto 10
elseif(iv(1)-2==0)then
goto 30
else
goto 100
endif
90 call itsmry(v(xyzzyaaad8:),iv,p,v,x)
100 return
end subroutine nl2sno
subroutine nl2itr(d,iv,jac,n,nn,p,r,v,x)
integer,intent(in) :: n, nn, p
integer,intent(inout) :: iv(:)
real(dp),intent(inout) :: d(:),jac(:),r(:),v(:),x(:)
integer xyzzyaaaa10,xyzzyaaab10,xyzzyaaac10,xyzzyaaad10,xyzzyaaae10,xy&
&zzyaaaf10,xyzzyaaag10,xyzzyaaah10,xyzzyaaai10,xyzzyaaaj10,xyzzyaaak10&
&,xyzzyaaal10,xyzzyaaam10,xyzzyaaan10,xyzzyaaao10,xyzzyaaap10,xyzzyaaa&
&q10,xyzzyaaar10,xyzzyaaas10,xyzzyaaat10,xyzzyaaau10,xyzzyaaav10,xyzzy&
&aaaw10,xyzzyaaax10,xyzzyaaay10,xyzzyaaaz10,xyzzyaaba10,xyzzyaabb10,xy&
&zzyaabc10,xyzzyaabd10,xyzzyaabe10,xyzzyaabf10,xyzzyaabg10
real(dp) xyzzyaabh10,xyzzyaabi10,xyzzyaabj10,xyzzyaabk10,xyzzyaabl10
real(dp),allocatable :: xyzzyaabm10(:,:)
integer,parameter :: xyzzyaabn10=34, xyzzyaabo10=26, xyzzyaabp10=14, x&
&yzzyaabq10=15, xyzzyaabr10=43, xyzzyaabs10=16, xyzzyaabt10=28, xyzzya&
&abu10=44, xyzzyaabv10=32, xyzzyaabw10=25, xyzzyaabx10=61,  xyzzyaaby1&
&0=60, xyzzyaabz10=3, xyzzyaaca10=35, xyzzyaacb10=36, xyzzyaacc10=37, &
&xyzzyaacd10=58, xyzzyaace10=38, xyzzyaacf10=5, xyzzyaacg10=17, xyzzya&
&ach10=18, xyzzyaaci10=6,    xyzzyaacj10=7, xyzzyaack10=40, xyzzyaacl1&
&0=41, xyzzyaacm10=30, xyzzyaacn10=31,   xyzzyaaco10=49, xyzzyaacp10=8&
&, xyzzyaacq10=51, xyzzyaacr10=9, xyzzyaacs10=52, xyzzyaact10=53,   xy&
&zzyaacu10=55, xyzzyaacv10=11, xyzzyaacw10=56, xyzzyaacx10=57, xyzzyaa&
&cy10=12,  xyzzyaacz10=2, xyzzyaada10=59, xyzzyaadb10=13, xyzzyaadc10=&
&60
integer,parameter :: xyzzyaadd10=43, xyzzyaade10=1, xyzzyaadf10=38, xy&
&zzyaadg10=2, xyzzyaadh10=37,  xyzzyaadi10=10, xyzzyaadj10=11, xyzzyaa&
&dk10=45, xyzzyaadl10=13, xyzzyaadm10=4, xyzzyaadn10=23,  xyzzyaado10=&
&39, xyzzyaadp10=87, xyzzyaadq10=35, xyzzyaadr10=9, xyzzyaads10=21,  x&
&yzzyaadt10=7, xyzzyaadu10=16, xyzzyaadv10=8, xyzzyaadw10=9, xyzzyaadx&
&10=42,    xyzzyaady10=47, xyzzyaadz10=5, xyzzyaaea10=29, xyzzyaaeb10=&
&30, xyzzyaaec10=78,  xyzzyaaed10=48
real(dp),parameter :: xyzzyaaee10=0.5d0,xyzzyaaef10=-1.d0,xyzzyaaeg10=&
&1.d0,xyzzyaaeh10=0.d0
allocate(xyzzyaabm10(nn,p))
xyzzyaaan10=0
do xyzzyaaal10=1,p
do xyzzyaaaf10=1,nn
xyzzyaaan10=xyzzyaaan10+1
xyzzyaabm10(xyzzyaaaf10,xyzzyaaal10)=jac(xyzzyaaan10)
enddo
enddo
xyzzyaaaf10=iv(1)
if(xyzzyaaaf10==1)goto 20
if(xyzzyaaaf10==2)goto 50
call parchk(iv,n,nn,p,v)
xyzzyaaaf10=iv(1)-2
select case(xyzzyaaaf10)
case (1:6)
goto 300
case (7,9)
goto 150
case (8)
goto 100
case (10)
goto 10
case default
goto 590
end select
10 iv(xyzzyaacn10)=0
iv(xyzzyaaci10)=1
iv(xyzzyaacm10)=1
iv(xyzzyaacj10)=1
iv(xyzzyaace10)=-1
iv(xyzzyaacv10)=2
iv(xyzzyaacz10)=0
iv(xyzzyaabn10)=0
iv(xyzzyaabo10)=0
iv(xyzzyaack10)=0
iv(xyzzyaacl10)=0
iv(xyzzyaacb10)=-1
iv(xyzzyaacp10)=0
iv(xyzzyaact10)=xyzzyaadp10+2*p
xyzzyaaas10=p*(p+1)/2
iv(xyzzyaadc10)=iv(xyzzyaact10)+xyzzyaaas10
iv(xyzzyaacu10)=iv(xyzzyaadc10)+p
iv(xyzzyaacw10)=iv(xyzzyaacu10)+p
iv(xyzzyaabr10)=iv(xyzzyaacw10)+p
iv(xyzzyaabt10)=iv(xyzzyaabr10)+p
iv(xyzzyaacc10)=iv(xyzzyaabt10)+p
iv(xyzzyaacq10)=iv(xyzzyaacc10)+p
iv(xyzzyaacs10)=iv(xyzzyaacq10)+p
iv(xyzzyaaco10)=iv(xyzzyaacs10)+n
iv(xyzzyaabu10)=iv(xyzzyaaco10)+n
iv(xyzzyaada10)=iv(xyzzyaabu10)+xyzzyaaas10
iv(xyzzyaacd10)=iv(xyzzyaada10)+4*p+7
if(v(xyzzyaadf10)>=xyzzyaaeh10)call vscopy(p,d(1:p),v(xyzzyaadf10))
if(v(xyzzyaado10)>xyzzyaaeh10)call vscopy(p,v(xyzzyaadp10:xyzzyaadp10+&
&p-1),v(xyzzyaado10))
xyzzyaaaf10=xyzzyaadp10+p
if(v(xyzzyaadh10)>xyzzyaaeh10)call vscopy(p,v(xyzzyaaaf10:xyzzyaaaf10+&
&p-1),v(xyzzyaadh10))
v(xyzzyaadw10)=xyzzyaaeh10
v(xyzzyaadz10)=xyzzyaaeh10
v(xyzzyaadv10)=v(xyzzyaadq10)/(xyzzyaaeg10+v(xyzzyaads10))
iv(xyzzyaacf10)=1
if(iv(xyzzyaabw10)==2)iv(xyzzyaacf10)=2
xyzzyaabc10=iv(xyzzyaact10)
if(iv(xyzzyaabw10)==0)call vscopy(xyzzyaaas10,v(xyzzyaabc10:xyzzyaabc1&
&0+xyzzyaaas10-1),xyzzyaaeh10)
20 xyzzyaabk10=v2norm(n,r)
if(xyzzyaabk10>v(xyzzyaadx10))iv(xyzzyaacz10)=1
if(iv(xyzzyaacz10)/=0)goto 30
v(xyzzyaadi10)=xyzzyaaee10*xyzzyaabk10**2
30 if(iv(xyzzyaace10)<0)then
goto 40
elseif(iv(xyzzyaace10)==0)then
goto 300
else
goto 540
endif
40 if(iv(xyzzyaacz10)==0)goto 60
iv(1)=13
goto 580
50 if(iv(xyzzyaacj10)/=0)goto 60
iv(1)=15
goto 580
60 iv(xyzzyaacb10)=-1
xyzzyaaab10=iv(xyzzyaabt10)
do xyzzyaaaf10=1,p
v(xyzzyaaab10)=dotprd(n,r,xyzzyaabm10(:,xyzzyaaaf10))
xyzzyaaab10=xyzzyaaab10+1
enddo
if(iv(xyzzyaace10)>0) goto 520
if(iv(xyzzyaabs10)>0)call dupdat(d(1:p),iv,xyzzyaabm10,n,p,v)
xyzzyaaax10=iv(xyzzyaacs10)
call vcopy (n,v(xyzzyaaax10:xyzzyaaax10+n-1),r)
xyzzyaaat10=iv(xyzzyaaco10)
call vcopy(n,v(xyzzyaaat10:xyzzyaaat10+n-1),r)
xyzzyaaab10=iv(xyzzyaabt10)
xyzzyaaaa10=iv(xyzzyaabr10)
xyzzyaaal10=xyzzyaaaa10
do xyzzyaaaf10=1,p
v(xyzzyaaal10)=v(xyzzyaaab10)/d(xyzzyaaaf10)
xyzzyaaal10=xyzzyaaal10+1
xyzzyaaab10=xyzzyaaab10+1
enddo
v(xyzzyaade10)=v2norm(p,v(xyzzyaaaa10:))
if(iv(xyzzyaabn10)/=0)goto 510
if(iv(xyzzyaace10)==0)goto 460
iv(xyzzyaace10)=0
90 call itsmry(d,iv,p,v,x)
100 xyzzyaaal10=iv(xyzzyaacn10)
if(xyzzyaaal10<iv(xyzzyaach10))goto 110
iv(1)=10
goto 580
110 iv(xyzzyaacn10)=xyzzyaaal10+1
if(xyzzyaaal10==0) goto 130
xyzzyaaba10=iv(xyzzyaacu10)
do xyzzyaaaf10=1,p
v(xyzzyaaba10)=d(xyzzyaaaf10)*v(xyzzyaaba10)
xyzzyaaba10=xyzzyaaba10+1
enddo
xyzzyaaba10=iv(xyzzyaacu10)
v(xyzzyaadv10)=v(xyzzyaadu10)*v2norm(p,v(xyzzyaaba10:))
130 xyzzyaabg10=iv(xyzzyaadc10)
v(xyzzyaadl10)=v(xyzzyaadi10)
iv(xyzzyaaca10)=-1
iv(xyzzyaabz10)=4
iv(xyzzyaabu10)=-abs(iv(xyzzyaabu10))
iv(xyzzyaacx10)=iv(xyzzyaacf10)
call vcopy(p,v(xyzzyaabg10:xyzzyaabg10+p-1),x)
140 if(.not.stopx())goto 160
iv(1)=11
goto 170
150 if(v(xyzzyaadi10)>=v(xyzzyaadl10))goto 160
v(xyzzyaadu10)=xyzzyaaeg10
xyzzyaaal10=iv(xyzzyaacn10)
goto 110
160 if(iv(xyzzyaaci10)<iv(xyzzyaacg10)+iv(xyzzyaack10))goto 180
iv(1)=9
170 if(v(xyzzyaadi10)>=v(xyzzyaadl10))goto 580
iv(xyzzyaabn10)=iv(1)
goto 450
180 xyzzyaaba10=iv(xyzzyaacu10)
xyzzyaabf10=iv(xyzzyaada10)
if(iv(xyzzyaacf10)==2)goto 220
xyzzyaaat10=iv(xyzzyaaco10)
if(iv(xyzzyaacb10)>=0)goto 190
xyzzyaaaw10=iv(xyzzyaacq10)
if(-1==iv(xyzzyaacb10))call qrfact(n,p,xyzzyaabm10,v(xyzzyaaaw10:xyzzy&
&aaaw10+p-1),iv(xyzzyaabx10:xyzzyaabx10+p-1),iv(xyzzyaabv10),0,v(xyzzy&
&aabf10:xyzzyaabf10+p-1))
call qapply(n,p,xyzzyaabm10,v(xyzzyaaat10:),iv(xyzzyaabv10))
190 xyzzyaaae10=iv(xyzzyaabu10)
if(xyzzyaaae10>0)goto 210
xyzzyaaae10=-xyzzyaaae10
iv(xyzzyaabu10)=xyzzyaaae10
xyzzyaaal10=xyzzyaaae10
xyzzyaaaw10=iv(xyzzyaacq10)
v(xyzzyaaal10)=v(xyzzyaaaw10)
if(p==1)goto 210
do xyzzyaaaf10=2,p
call vcopy(xyzzyaaaf10-1,v(xyzzyaaal10+1:xyzzyaaal10+1+xyzzyaaaf10-1-1&
&),xyzzyaabm10(:,xyzzyaaaf10))
xyzzyaaal10=xyzzyaaal10+xyzzyaaaf10
xyzzyaaaw10=xyzzyaaaw10+1
v(xyzzyaaal10)=v(xyzzyaaaw10)
enddo
210 xyzzyaaab10=iv(xyzzyaabt10)
call lmstep(d,v(xyzzyaaab10:xyzzyaaab10+p-1),iv(xyzzyaabv10),iv(xyzzya&
&abx10:xyzzyaabx10+p-1),iv(xyzzyaacb10),p,v(xyzzyaaat10:xyzzyaaat10+n-&
&1),v(xyzzyaaae10:xyzzyaaae10+(p*(p+1))/2-1),v(xyzzyaaba10:xyzzyaaba10&
&+p-1),v(1:86),v(xyzzyaabf10:))
goto 290
220 if(iv(xyzzyaabu10)>0)goto 280
xyzzyaaae10=-iv(xyzzyaabu10)
iv(xyzzyaabu10)=xyzzyaaae10
xyzzyaabc10=iv(xyzzyaact10)
if(-1/=iv(xyzzyaacb10))goto 250
do xyzzyaaaf10=1,p
xyzzyaabk10=xyzzyaaeg10/d(xyzzyaaaf10)
do xyzzyaaal10=1,xyzzyaaaf10
v(xyzzyaaae10)=xyzzyaabk10*(dotprd(n,xyzzyaabm10(:,xyzzyaaaf10),xyzzya&
&abm10(:,xyzzyaaal10))+v(xyzzyaabc10))/d(xyzzyaaal10)
xyzzyaaae10=xyzzyaaae10+1
xyzzyaabc10=xyzzyaabc10+1
enddo
enddo
goto 280
250 xyzzyaaay10=xyzzyaabc10-xyzzyaaae10
xyzzyaaad10=xyzzyaaae10-1
xyzzyaaaj10=iv(xyzzyaabx10)
xyzzyaabl10=xyzzyaaeg10/d(xyzzyaaaj10)
xyzzyaaav10=iv(xyzzyaacq10)-1
xyzzyaabi10=v(xyzzyaaav10+1)
do xyzzyaaaf10=1,p
xyzzyaaan10=xyzzyaaby10+xyzzyaaaf10
xyzzyaaah10=iv(xyzzyaaan10)
xyzzyaaae10=xyzzyaaad10+xyzzyaaah10*(xyzzyaaah10-1)/2
xyzzyaaan10=xyzzyaaae10+xyzzyaaah10
xyzzyaaar10=xyzzyaaan10+xyzzyaaay10
xyzzyaabk10=xyzzyaaeg10/d(xyzzyaaah10)
xyzzyaaau10=xyzzyaaav10+xyzzyaaaf10
xyzzyaabh10=v(xyzzyaaau10)**2
if(xyzzyaaaf10>1)xyzzyaabh10=xyzzyaabh10+dotprd(xyzzyaaaf10-1,xyzzyaab&
&m10(:,xyzzyaaaf10),xyzzyaabm10(:,xyzzyaaaf10))
v(xyzzyaaan10)=(xyzzyaabh10+v(xyzzyaaar10))*xyzzyaabk10**2
if(xyzzyaaaf10==1)cycle
xyzzyaaan10=xyzzyaaae10+xyzzyaaaj10
if(xyzzyaaah10<xyzzyaaaj10)xyzzyaaan10=xyzzyaaan10+((xyzzyaaaj10-xyzzy&
&aaah10)*(xyzzyaaaj10+xyzzyaaah10-3))/2
xyzzyaaar10=xyzzyaaan10+xyzzyaaay10
v(xyzzyaaan10)=xyzzyaabk10*(xyzzyaabi10*xyzzyaabm10(1,xyzzyaaaf10)+v(x&
&yzzyaaar10))*xyzzyaabl10
if(xyzzyaaaf10==2)cycle
xyzzyaaag10=xyzzyaaaf10-1
do xyzzyaaal10=2,xyzzyaaag10
xyzzyaaak10=xyzzyaaby10+xyzzyaaal10
xyzzyaaai10=iv(xyzzyaaak10)
xyzzyaaan10=xyzzyaaae10+xyzzyaaai10
if(xyzzyaaah10<xyzzyaaai10)xyzzyaaan10=xyzzyaaan10+((xyzzyaaai10-xyzzy&
&aaah10)*(xyzzyaaai10+xyzzyaaah10-3))/2
xyzzyaaar10=xyzzyaaan10+xyzzyaaay10
xyzzyaaam10=xyzzyaaal10-1
xyzzyaaau10=xyzzyaaav10+xyzzyaaal10
v(xyzzyaaan10)=xyzzyaabk10*(dotprd(xyzzyaaam10,xyzzyaabm10(:,xyzzyaaaf&
&10),xyzzyaabm10(:,xyzzyaaal10))+v(xyzzyaaau10)*xyzzyaabm10(xyzzyaaal1&
&0,xyzzyaaaf10)+v(xyzzyaaar10))/d(xyzzyaaai10)
enddo
enddo
280 xyzzyaaae10=iv(xyzzyaabu10)
xyzzyaaaa10=iv(xyzzyaabr10)
xyzzyaaap10=iv(xyzzyaacd10)
call gqtstp(d,v(xyzzyaaaa10:xyzzyaaaa10+p-1),v(xyzzyaaae10:xyzzyaaae10&
&+(p*(p+1))/2-1),iv(xyzzyaaca10),v(xyzzyaaap10:xyzzyaaap10+(p*(p+1))/2&
&-1),p,v(xyzzyaaba10:xyzzyaaba10+p-1),v(1:86),v(xyzzyaabf10:xyzzyaabf1&
&0+4*p+7-1))
290 if(iv(xyzzyaabz10)==6)goto 300
xyzzyaabg10=iv(xyzzyaadc10)
xyzzyaaba10=iv(xyzzyaacu10)
call vaxpy(p,x(1:p),xyzzyaaeg10,v(xyzzyaaba10:xyzzyaaba10+p-1),v(xyzzy&
&aabg10:xyzzyaabg10+p-1))
iv(xyzzyaaci10)=iv(xyzzyaaci10)+1
iv(1)=1
iv(xyzzyaacz10)=0
goto 590
300 xyzzyaaba10=iv(xyzzyaacu10)
xyzzyaaaq10=iv(xyzzyaacw10)
xyzzyaabg10=iv(xyzzyaadc10)
call assess(d,iv,p,v(xyzzyaaba10:xyzzyaaba10+p-1),v(xyzzyaaaq10:xyzzya&
&aaq10+p-1),v(1:86),x,v(xyzzyaabg10:xyzzyaabg10+p-1))
if(iv(xyzzyaacy10)==0)goto 310
iv(xyzzyaabu10)=-abs(iv(xyzzyaabu10))
iv(xyzzyaacx10)=iv(xyzzyaacx10)+2
call vcopy(xyzzyaadr10,v(1:xyzzyaadr10),v(xyzzyaaec10:))
310 if(iv(xyzzyaacr10)==0)goto 320
xyzzyaaax10=iv(xyzzyaacs10)
call vcopy(n,r(1:n),v(xyzzyaaax10:))
320 xyzzyaaan10=iv(xyzzyaabz10)-4
xyzzyaabb10=iv(xyzzyaacf10)
if(xyzzyaaan10>0)then
select case(xyzzyaaan10)
case (1)
goto 340
case (2)
goto 360
case (3:8)
goto 370
case (9)
goto 500
case (10)
goto 460
end select
endif
xyzzyaabh10=v(xyzzyaadt10)-v(xyzzyaadj10)
xyzzyaaaz10=iv(xyzzyaacc10)
xyzzyaabc10=iv(xyzzyaact10)
call slvmul(p,v(xyzzyaaaz10:xyzzyaaaz10+p*(p+1)/2-1),v(xyzzyaabc10:xyz&
&zyaabc10+(p*(p+1))/2-1),v(xyzzyaaba10:xyzzyaaba10+p-1))
xyzzyaabj10=xyzzyaaee10*dotprd(p,v(xyzzyaaba10:),v(xyzzyaaaz10:))
if(iv(xyzzyaacf10)==1)xyzzyaabj10=-xyzzyaabj10
if(abs(xyzzyaabh10+xyzzyaabj10)*v(xyzzyaadk10)>=abs(xyzzyaabh10))goto &
&330
iv(xyzzyaacf10)=3-iv(xyzzyaacf10)
if(iv(xyzzyaacf10)==1)iv(xyzzyaaca10)=-1
if(iv(xyzzyaacf10)==2.and.iv(xyzzyaacb10)>0)iv(xyzzyaacb10)=0
if(-2<xyzzyaaan10) goto 380
iv(xyzzyaabu10)=-abs(iv(xyzzyaabu10))
iv(xyzzyaacx10)=iv(xyzzyaacx10)+2
call vcopy(xyzzyaadr10,v(xyzzyaaec10:xyzzyaaec10+xyzzyaadr10-1),v)
goto 350
330 if(-3<xyzzyaaan10)goto 380
v(xyzzyaadv10)=v(xyzzyaadu10)*v(xyzzyaadg10)
goto 140
340 v(xyzzyaadv10)=v(xyzzyaadu10)*v(xyzzyaadg10)
350 if(v(xyzzyaadi10)>=v(xyzzyaadl10))goto 140
xyzzyaaax10=iv(xyzzyaacs10)
call vcopy(n,v(xyzzyaaax10:xyzzyaaax10+n-1),r)
goto 140
360 v(xyzzyaadv10)=v(xyzzyaadq10)
goto 180
370 iv(xyzzyaabn10)=xyzzyaaan10
if(v(xyzzyaadi10)>=v(xyzzyaadl10))goto 510
if(iv(xyzzyaadb10)==14)goto 510
iv(xyzzyaadb10)=14
380 iv(xyzzyaabo10)=0
xyzzyaaao10=iv(xyzzyaacc10)
if(iv(xyzzyaacb10)>=0)goto 400
do xyzzyaaaf10=1,p
v(xyzzyaaao10)=dotprd(n,xyzzyaabm10(:,xyzzyaaaf10),r)
xyzzyaaao10=xyzzyaaao10+1
enddo
goto 410
400 xyzzyaaat10=iv(xyzzyaaco10)
call vcopy(n,v(xyzzyaaat10:xyzzyaaat10+n-1),r)
call qapply(n,p,xyzzyaabm10,v(xyzzyaaat10:xyzzyaaat10+n-1),iv(xyzzyaab&
&v10))
xyzzyaaaw10=iv(xyzzyaacq10)
xyzzyaabd10=iv(xyzzyaacw10)
call rptmul(3,iv(xyzzyaabx10:xyzzyaabx10+p-1),xyzzyaabm10,p,v(xyzzyaaa&
&w10:xyzzyaaaw10+p-1),v(xyzzyaaat10:xyzzyaaat10+p-1),v(xyzzyaaao10:xyz&
&zyaaao10+p-1),v(xyzzyaabd10:xyzzyaabd10+p-1))
410 if(iv(xyzzyaabz10)/=3)goto 450
xyzzyaaba10=iv(xyzzyaacu10)
xyzzyaabd10=iv(xyzzyaacw10)
xyzzyaabe10=iv(xyzzyaadc10)
if(xyzzyaabb10==2)goto 420
xyzzyaaaw10=iv(xyzzyaacq10)
call rptmul(2,iv(xyzzyaabx10:xyzzyaabx10+p-1),xyzzyaabm10,p,v(xyzzyaaa&
&w10:xyzzyaaaw10+p-1),v(xyzzyaaba10:xyzzyaaba10+p-1),v(xyzzyaabd10:xyz&
&zyaabd10+p-1),v(xyzzyaabe10:xyzzyaabe10+p-1))
goto 450
420 xyzzyaaae10=iv(xyzzyaabu10)
xyzzyaaal10=xyzzyaabe10
do xyzzyaaaf10=1,p
v(xyzzyaaal10)=d(xyzzyaaaf10)*v(xyzzyaaba10)
xyzzyaaal10=xyzzyaaal10+1
xyzzyaaba10=xyzzyaaba10+1
enddo
call slvmul(p,v(xyzzyaabd10:xyzzyaabd10+p*(p+1)/2-1),v(xyzzyaaae10:),v&
&(xyzzyaabe10:))
do xyzzyaaaf10=1,p
v(xyzzyaabd10)=d(xyzzyaaaf10)*v(xyzzyaabd10)
xyzzyaabd10=xyzzyaabd10+1
enddo
450 iv(xyzzyaacm10)=iv(xyzzyaacm10)+1
xyzzyaaab10=iv(xyzzyaabt10)
xyzzyaaac10=iv(xyzzyaada10)
call vcopy(p,v(xyzzyaaac10:xyzzyaaac10+p-1),v(xyzzyaaab10:))
iv(1)=2
goto 590
460 xyzzyaaac10=iv(xyzzyaada10)
xyzzyaaab10=iv(xyzzyaabt10)
call vaxpy(p,v(xyzzyaaac10:xyzzyaaac10+p-1),xyzzyaaef10,v(xyzzyaaac10:&
&),v(xyzzyaaab10:))
xyzzyaaba10=iv(xyzzyaacu10)
xyzzyaabd10=iv(xyzzyaacw10)
xyzzyaabe10=iv(xyzzyaadc10)
if(iv(xyzzyaabz10)/=3)goto 490
xyzzyaaal10=xyzzyaabd10
xyzzyaaan10=xyzzyaaac10
do xyzzyaaaf10=1,p
v(xyzzyaaal10)=(v(xyzzyaaal10)-v(xyzzyaaan10))/d(xyzzyaaaf10)
xyzzyaaal10=xyzzyaaal10+1
xyzzyaaan10=xyzzyaaan10+1
enddo
if(v2norm(p,v(xyzzyaabd10:))<=v(xyzzyaade10)*v(xyzzyaaea10))goto 480
if(dotprd(p,v(xyzzyaaab10:),v(xyzzyaaba10:))>=v(xyzzyaadm10)*v(xyzzyaa&
&eb10))goto 490
480 v(xyzzyaadu10)=v(xyzzyaadn10)
490 xyzzyaaao10=iv(xyzzyaacc10)
call vaxpy (p,v(xyzzyaaao10:xyzzyaaao10+p-1),xyzzyaaef10,v(xyzzyaaao10&
&:),v(xyzzyaaab10:))
xyzzyaabc10=iv(xyzzyaact10)
call slvmul(p,v(xyzzyaabd10:xyzzyaabd10+(p*(p+1)/2)-1),v(xyzzyaabc10:)&
&,v(xyzzyaaba10:))
xyzzyaabl10=abs(dotprd(p,v(xyzzyaaba10:),v(xyzzyaabd10:)))
xyzzyaabk10=abs(dotprd(p,v(xyzzyaaba10:),v(xyzzyaaao10:)))
v(xyzzyaady10)=xyzzyaaeg10
if(xyzzyaabk10<xyzzyaabl10)v(xyzzyaady10)=xyzzyaabk10/xyzzyaabl10
call slupdt(v(xyzzyaabc10:xyzzyaabc10+(p*(p+1))/2-1),v(xyzzyaadd10),p,&
&v(xyzzyaady10),v(xyzzyaaba10:xyzzyaaba10+p-1),v(xyzzyaabd10:xyzzyaabd&
&10+p-1),v(xyzzyaabe10:xyzzyaabe10+p-1),v(xyzzyaaac10:xyzzyaaac10+p-1)&
&,v(xyzzyaaed10),v(xyzzyaaao10:xyzzyaaao10+p-1))
iv(1)=2
goto 90
500 iv(1)=14
goto 580
510 if(iv(xyzzyaabq10)==0.and.iv(xyzzyaabp10)==0)goto 570
if(iv(xyzzyaabo10)/=0)goto 570
if(iv(xyzzyaabn10)>=7)goto 570
iv(xyzzyaace10)=0
520 call covclc(xyzzyaaaf10,d,iv,xyzzyaabm10,n,p,r,v,x)
select case(xyzzyaaaf10)
case(1:2)
iv(xyzzyaack10)=iv(xyzzyaack10)+1
iv(xyzzyaaci10)=iv(xyzzyaaci10)+1
iv(xyzzyaacr10)=xyzzyaaaf10
iv(1)=1
goto 590
case(3)
goto 550
case(4)
iv(xyzzyaace10)=0
if(iv(xyzzyaacn10)==0)iv(xyzzyaace10)=-1
goto 570
end select
540 if(iv(xyzzyaacr10)==1.or.iv(xyzzyaacz10)/=0)goto 520
iv(xyzzyaacj10)=iv(xyzzyaaci10)
550 iv(xyzzyaacl10)=iv(xyzzyaacl10)+1
iv(xyzzyaacm10)=iv(xyzzyaacm10)+1
iv(1)=2
goto 590
570 iv(1)=iv(xyzzyaabn10)
iv(xyzzyaabn10)=0
580 call itsmry(d,iv,p,v,x)
590 xyzzyaaan10=0
do xyzzyaaal10=1,p
do xyzzyaaaf10=1,nn
xyzzyaaan10=xyzzyaaan10+1
jac(xyzzyaaan10)=xyzzyaabm10(xyzzyaaaf10,xyzzyaaal10)
enddo
enddo
deallocate(xyzzyaabm10)
return
end subroutine nl2itr
subroutine assess(d,iv,p,step,stlstg,v,x,x0)
integer,intent(in) :: p
integer,intent(inout) :: iv(:)
real(dp),intent(in) :: d(:),x0(:)
real(dp),intent(inout) :: step(:),stlstg(:),v(:),x(:)
logical xyzzyaaaa11
integer xyzzyaaab11,xyzzyaaac11
real(dp) xyzzyaaad11,xyzzyaaae11,xyzzyaaaf11,xyzzyaaag11,xyzzyaaah11
real(dp),parameter :: xyzzyaaai11=0.5d0,xyzzyaaaj11=1.d0,xyzzyaaak11=2&
&.d0,xyzzyaaal11=0.d0
integer,parameter :: xyzzyaaam11=3,xyzzyaaan11=4,xyzzyaaao11=5,xyzzyaa&
&ap11=6,xyzzyaaaq11=7,xyzzyaaar11=8, xyzzyaaas11=9,xyzzyaaat11=10,xyzz&
&yaaau11=11,xyzzyaaav11=12,xyzzyaaaw11=2,xyzzyaaax11=13
integer,parameter :: xyzzyaaay11=31,xyzzyaaaz11=22,xyzzyaaba11=2,xyzzy&
&aabb11=3,xyzzyaabc11=18,           xyzzyaabd11=10,xyzzyaabe11=11,xyzz&
&yaabf11=12,xyzzyaabg11=13,xyzzyaabh11=14,xyzzyaabi11=4,xyzzyaabj11=23&
&,xyzzyaabk11=35,xyzzyaabl11=6, xyzzyaabm11=15,xyzzyaabn11=7,xyzzyaabo&
&11=16,xyzzyaabp11=24,xyzzyaabq11=25,xyzzyaabr11=17,xyzzyaabs11=32,xyz&
&zyaabt11=5,xyzzyaabu11=26,xyzzyaabv11=27,xyzzyaabw11=28,xyzzyaabx11=3&
&3,xyzzyaaby11=34
xyzzyaaac11=iv(xyzzyaaap11)
iv(xyzzyaaav11)=0
iv(xyzzyaaas11)=0
xyzzyaaag11=xyzzyaaaj11
xyzzyaaaa11=.true.
xyzzyaaab11=iv(xyzzyaaam11)
select case (xyzzyaaab11)
case (1)
goto 20
case (2)
goto 30
case (3:4)
goto 10
case (5)
goto 40
case (6)
goto 280
case (7:11)
goto 220
case (12)
goto 170
case default
iv(xyzzyaaam11)=13
goto 300
end select
10 iv(xyzzyaaat11)=1
iv(xyzzyaaar11)=0
v(xyzzyaabf11)=v(xyzzyaabg11)
if(iv(xyzzyaaaw11)==0)goto 90
iv(xyzzyaaat11)=-1
iv(xyzzyaaax11)=xyzzyaaab11
goto 60
20 if(iv(xyzzyaaao11)/=iv(xyzzyaaan11))goto 30
iv(xyzzyaaat11)=iv(xyzzyaaau11)
iv(xyzzyaaar11)=-1
goto 90
30 iv(xyzzyaaat11)=iv(xyzzyaaat11)+1
40 if(iv(xyzzyaaat11)>0)goto 50
if(iv(xyzzyaaaw11)/=0)goto 60
iv(xyzzyaaat11)=-iv(xyzzyaaat11)
xyzzyaaab11=iv(xyzzyaaax11)
select case(xyzzyaaab11)
case(1)
goto 20
case(2)
goto 30
case(3:4)
goto 90
case(5)
goto 70
end select
50 if(iv(xyzzyaaaw11)==0)goto 70
if(iv(xyzzyaaar11)>0)goto 80
iv(xyzzyaaat11)=-iv(xyzzyaaat11)
iv(xyzzyaaax11)=iv(xyzzyaaam11)
60 v(xyzzyaabo11)=v(xyzzyaaaz11)
iv(xyzzyaaar11)=iv(xyzzyaaar11)-1
iv(xyzzyaaam11)=5
goto 300
70 if(v(xyzzyaabd11)<v(xyzzyaabf11))goto 90
if(iv(xyzzyaaao11)==iv(xyzzyaaan11))goto 80
iv(xyzzyaaao11)=iv(xyzzyaaan11)
iv(xyzzyaaav11)=1
80 if(v(xyzzyaabf11)>=v(xyzzyaabg11))goto 90
iv(xyzzyaaas11)=1
v(xyzzyaabd11)=v(xyzzyaabf11)
v(xyzzyaabn11)=v(xyzzyaabm11)
v(xyzzyaabi11)=v(xyzzyaabh11)
if(iv(xyzzyaaav11)==0)xyzzyaaag11=v(xyzzyaaba11)/v(xyzzyaabc11)
v(xyzzyaaba11)=v(xyzzyaabc11)
xyzzyaaac11=iv(xyzzyaaaq11)
xyzzyaaaa11=.false.
90 xyzzyaaaf11=reldst(p,d,x,x0)
if(xyzzyaaaa11)goto 110
do xyzzyaaab11=1,p
step(xyzzyaaab11)=stlstg(xyzzyaaab11)
x(xyzzyaaab11)=x0(xyzzyaaab11)+stlstg(xyzzyaaab11)
enddo
110 v(xyzzyaabe11)=v(xyzzyaabg11)-v(xyzzyaabd11)
if(v(xyzzyaabe11)>v(xyzzyaabv11)*v(xyzzyaabn11))goto 140
v(xyzzyaabr11)=xyzzyaaaf11
if(v(xyzzyaabd11)<v(xyzzyaabg11))goto 120
iv(xyzzyaaan11)=iv(xyzzyaaao11)
v(xyzzyaabf11)=v(xyzzyaabd11)
v(xyzzyaabd11)=v(xyzzyaabg11)
call vcopy(p,x(1:p),x0)
iv(xyzzyaaas11)=1
goto 130
120 iv(xyzzyaaaq11)=xyzzyaaac11
130 iv(xyzzyaaam11)=1
if(iv(xyzzyaaat11)<iv(xyzzyaaau11))goto 160
iv(xyzzyaaam11)=5
iv(xyzzyaaar11)=iv(xyzzyaaar11)-1
goto 160
140 iv(xyzzyaaaq11)=xyzzyaaac11
xyzzyaaag11=xyzzyaaaj11
if(xyzzyaaaa11)v(xyzzyaabr11)=xyzzyaaaf11
v(xyzzyaabc11)=v(xyzzyaaba11)
if(v(xyzzyaabe11)>v(xyzzyaabn11)*v(xyzzyaabu11))goto 190
if(iv(xyzzyaaat11)>=iv(xyzzyaaau11))goto 150
iv(xyzzyaaam11)=2
goto 160
150 iv(xyzzyaaam11)=4
160 iv(xyzzyaaax11)=iv(xyzzyaaam11)
xyzzyaaad11=v(xyzzyaabi11)+v(xyzzyaabe11)
v(xyzzyaabo11)=xyzzyaaai11*xyzzyaaag11
if(xyzzyaaad11<v(xyzzyaabi11))v(xyzzyaabo11)=xyzzyaaag11*max(v(xyzzyaa&
&bp11),xyzzyaaai11*v(xyzzyaabi11)/xyzzyaaad11)
170 if(v(xyzzyaabr11)<=v(xyzzyaaby11))goto 180
iv(xyzzyaaam11)=iv(xyzzyaaax11)
if(v(xyzzyaabd11)<v(xyzzyaabg11))goto 200
goto 230
180 iv(xyzzyaaam11)=12
goto 240
190 if(v(xyzzyaabe11)<(-v(xyzzyaabw11)*v(xyzzyaabi11)))goto 210
if(iv(xyzzyaaar11)<0)goto 210
if(iv(xyzzyaaas11)==1)goto 210
v(xyzzyaabo11)=v(xyzzyaabq11)
xyzzyaaae11=v(xyzzyaabi11)
if(v(xyzzyaabe11)<(xyzzyaaai11/v(xyzzyaabo11)-xyzzyaaaj11)*xyzzyaaae11&
&) v(xyzzyaabo11)=max(v(xyzzyaabj11),xyzzyaaai11*xyzzyaaae11/ (xyzzyaa&
&ae11+v(xyzzyaabe11)))
iv(xyzzyaaam11)=4
if(v(xyzzyaabt11)==xyzzyaaal11)goto 230
iv(xyzzyaaam11)=5
iv(xyzzyaaar11)=iv(xyzzyaaar11)+1
200 v(xyzzyaabf11)=v(xyzzyaabd11)
iv(xyzzyaaan11)=iv(xyzzyaaao11)
call vcopy(p,stlstg(1:p),step)
v(xyzzyaabc11)=v(xyzzyaaba11)
iv(xyzzyaaaq11)=xyzzyaaac11
v(xyzzyaabm11)=v(xyzzyaabn11)
v(xyzzyaabh11)=v(xyzzyaabi11)
goto 230
210 v(xyzzyaabo11)=xyzzyaaaj11
iv(xyzzyaaam11)=3
goto 230
220 iv(xyzzyaaam11)=iv(xyzzyaaax11)
if(v(xyzzyaabc11)>=xyzzyaaal11)goto 240
iv(xyzzyaaam11)=12
goto 240
230 iv(xyzzyaaax11)=iv(xyzzyaaam11)
240 if(abs(v(xyzzyaabd11))<v(xyzzyaaay11))iv(xyzzyaaam11)=10
if(xyzzyaaai11*v(xyzzyaabe11)>v(xyzzyaabn11))goto 300
xyzzyaaad11=v(xyzzyaabs11)*abs(v(xyzzyaabg11))
if(v(xyzzyaaba11)>v(xyzzyaabk11).and.v(xyzzyaabn11)<=xyzzyaaad11)iv(xy&
&zzyaaam11)=11
if(v(xyzzyaabb11)<xyzzyaaal11)goto 250
xyzzyaaab11=0
if((v(xyzzyaabl11)>xyzzyaaal11.and.v(xyzzyaabl11)<=xyzzyaaad11).or. (v&
&(xyzzyaabl11)==xyzzyaaal11.and.v(xyzzyaabn11)==xyzzyaaal11))xyzzyaaab&
&11=2
if(v(xyzzyaabt11)==xyzzyaaal11.and.v(xyzzyaabr11)<=v(xyzzyaabx11).and.&
&xyzzyaaaa11)xyzzyaaab11=xyzzyaaab11+1
if(xyzzyaaab11>0)iv(xyzzyaaam11)=xyzzyaaab11+6
250 if(abs(iv(xyzzyaaam11)-3)>2.and.iv(xyzzyaaam11)/=12)goto 300
if(v(xyzzyaaba11)>v(xyzzyaabk11))goto 260
if(v(xyzzyaabn11)>=xyzzyaaad11)goto 300
if(v(xyzzyaabb11)<=xyzzyaaal11)goto 270
if(xyzzyaaai11*v(xyzzyaabb11)<=v(xyzzyaabk11)) goto 300
goto 270
260 if(xyzzyaaai11*v(xyzzyaaba11)<=v(xyzzyaabk11))goto 300
xyzzyaaah11=v(xyzzyaabk11)/v(xyzzyaaba11)
if(xyzzyaaah11*(xyzzyaaak11-xyzzyaaah11)*v(xyzzyaabn11)>=xyzzyaaad11)g&
&oto 300
270 if(v(xyzzyaabl11)<xyzzyaaal11)goto 290
v(xyzzyaabh11)=v(xyzzyaabi11)
v(xyzzyaabc11)=v(xyzzyaaba11)
if(iv(xyzzyaaam11)==12)v(xyzzyaabc11)=-v(xyzzyaabc11)
v(xyzzyaabm11)=v(xyzzyaabn11)
iv(xyzzyaaam11)=6
call vcopy(p,stlstg(1:p),step)
goto 300
280 v(xyzzyaabi11)=v(xyzzyaabh11)
v(xyzzyaaba11)=abs(v(xyzzyaabc11))
call vcopy(p,step(1:p),stlstg)
iv(xyzzyaaam11)=iv(xyzzyaaax11)
if(v(xyzzyaabc11)<=xyzzyaaal11)iv(xyzzyaaam11)=12
v(xyzzyaabl11)=-v(xyzzyaabn11)
v(xyzzyaabn11)=v(xyzzyaabm11)
290 if(-v(xyzzyaabl11)<=v(xyzzyaabs11)*abs(v(xyzzyaabg11)))iv(xyzzyaaa&
&m11)=11
300 return
end subroutine assess
subroutine covclc(covirc,d,iv,j,n,p,r,v,x)
integer,intent(in) :: n, p
integer,intent(inout) :: covirc,iv(:)
real(dp),intent(in) :: d(:)
real(dp),intent(inout) :: v(:),x(:),j(:,:)
real(dp),intent(inout) :: r(:)
integer xyzzyaaaa12,xyzzyaaab12,xyzzyaaac12,xyzzyaaad12,xyzzyaaae12,xy&
&zzyaaaf12,xyzzyaaag12,xyzzyaaah12,xyzzyaaai12,xyzzyaaaj12,xyzzyaaak12&
&,xyzzyaaal12,xyzzyaaam12,xyzzyaaan12,xyzzyaaao12,xyzzyaaap12,xyzzyaaa&
&q12,xyzzyaaar12,xyzzyaaas12,xyzzyaaat12,xyzzyaaau12,xyzzyaaav12,xyzzy&
&aaaw12,xyzzyaaax12,xyzzyaaay12,xyzzyaaaz12,xyzzyaaba12,xyzzyaabb12
integer,save :: xyzzyaabc12
real(dp) xyzzyaabd12,xyzzyaabe12,xyzzyaabf12
logical xyzzyaabg12
real(dp),parameter :: xyzzyaabh12=0.5d0,xyzzyaabi12=-0.5d0,xyzzyaabj12&
&=1.d0,xyzzyaabk12=2.d0,xyzzyaabl12=0.d0
integer,parameter :: xyzzyaabm12=26, xyzzyaabn12=15, xyzzyaabo12=50, x&
&yzzyaabp12=44, xyzzyaabq12=40,  xyzzyaabr12=10, xyzzyaabs12=46, xyzzy&
&aabt12=28, xyzzyaabu12=44, xyzzyaabv12=32, xyzzyaabw12=61, xyzzyaabx1&
&2=60, xyzzyaaby12=35, xyzzyaabz12=36, xyzzyaaca12=58, xyzzyaacb12=38,&
& xyzzyaacc12=7, xyzzyaacd12=49, xyzzyaace12=51, xyzzyaacf12=52, xyzzy&
&aacg12=54, xyzzyaach12=12, xyzzyaaci12=2, xyzzyaacj12=59,  xyzzyaack1&
&2=49
covirc=4
xyzzyaaan12=iv(xyzzyaabn12)
xyzzyaaaq12=iv(xyzzyaacb12)
if(xyzzyaaaq12>0)goto 10
iv(xyzzyaaby12)=-1
if(iv(xyzzyaabz12)>0)iv(xyzzyaabz12)=0
if(abs(xyzzyaaan12)>=3)goto 310
v(xyzzyaabs12)=v(xyzzyaabr12)
xyzzyaaam12=iv(xyzzyaacf12)
call vcopy (n,v(xyzzyaaam12:xyzzyaaam12+n-1),r)
10 if(xyzzyaaaq12>p)goto 220
if(xyzzyaaan12<0)goto 110
xyzzyaaab12=iv(xyzzyaacj12)+p
xyzzyaaac12=iv(xyzzyaabt12)
if(xyzzyaaaq12>0)goto 20
call vcopy(p,v(xyzzyaaab12:xyzzyaaab12+p-1),v(xyzzyaaac12:))
iv(xyzzyaach12)=iv(xyzzyaacc12)
goto 90
20 xyzzyaabd12=v(xyzzyaabo12)
x(xyzzyaaaq12)=v(xyzzyaack12)
if(iv(xyzzyaaci12)==0)goto 40
if(xyzzyaabd12*x(xyzzyaaaq12)>xyzzyaabl12)goto 30
iv(xyzzyaabm12)=-2
goto 210
30 xyzzyaabd12=xyzzyaabi12*xyzzyaabd12
goto 100
40 xyzzyaabc12=iv(xyzzyaaca12)
xyzzyaaaa12=xyzzyaaac12+p-1
do xyzzyaaah12=xyzzyaaac12,xyzzyaaaa12
v(xyzzyaaah12)=(v(xyzzyaaah12)-v(xyzzyaaab12))/xyzzyaabd12
xyzzyaaab12=xyzzyaaab12+1
enddo
xyzzyaaam12=xyzzyaabc12+xyzzyaaaq12*(xyzzyaaaq12-1)/2
xyzzyaaap12=xyzzyaaam12+xyzzyaaaq12-2
if(xyzzyaaaq12==1)goto 70
do xyzzyaaah12=xyzzyaaam12,xyzzyaaap12
v(xyzzyaaah12)=xyzzyaabh12*(v(xyzzyaaah12)+v(xyzzyaaac12))
xyzzyaaac12=xyzzyaaac12+1
enddo
70 xyzzyaaap12=xyzzyaaap12+1
do xyzzyaaah12=xyzzyaaaq12,p
v(xyzzyaaap12)=v(xyzzyaaac12)
xyzzyaaap12=xyzzyaaap12+xyzzyaaah12
xyzzyaaac12=xyzzyaaac12+1
enddo
90 xyzzyaaaq12=xyzzyaaaq12+1
iv(xyzzyaacb12)=xyzzyaaaq12
if(xyzzyaaaq12>p)goto 210
xyzzyaabd12=v(xyzzyaabp12)*max(xyzzyaabj12/d(xyzzyaaaq12),abs(x(xyzzya&
&aaq12)))
if(x(xyzzyaaaq12)<xyzzyaabl12)xyzzyaabd12=-xyzzyaabd12
v(xyzzyaack12)=x(xyzzyaaaq12)
100 x(xyzzyaaaq12)=x(xyzzyaaaq12)+xyzzyaabd12
v(xyzzyaabo12)=xyzzyaabd12
covirc=2
goto 390
110 xyzzyaaay12=iv(xyzzyaacj12)+p-1
xyzzyaaar12=xyzzyaaaq12-1
xyzzyaaas12=xyzzyaaaq12*xyzzyaaar12/2
if(xyzzyaaaq12>0)goto 120
iv(xyzzyaacg12)=0
goto 200
120 xyzzyaaah12=iv(xyzzyaacg12)
if(xyzzyaaah12>0)goto 180
if(iv(xyzzyaaci12)==0)goto 140
xyzzyaaax12=xyzzyaaay12+xyzzyaaaq12
xyzzyaabd12=v(xyzzyaaax12)
if(xyzzyaabd12*x(xyzzyaack12)>xyzzyaabl12)goto 130
iv(xyzzyaabm12)=-2
goto 390
130 xyzzyaabd12=xyzzyaabi12*xyzzyaabd12
x(xyzzyaaaq12)=x(xyzzyaack12)+xyzzyaabd12
v(xyzzyaaax12)=xyzzyaabd12
covirc=1
goto 390
140 xyzzyaaat12=p*(p-1)/2
xyzzyaabc12=iv(xyzzyaaca12)
xyzzyaaag12=xyzzyaabc12+xyzzyaaat12+xyzzyaaar12
v(xyzzyaaag12)=v(xyzzyaabr12)
xyzzyaaae12=xyzzyaabc12+xyzzyaaas12
if(xyzzyaaar12==0)goto 160
xyzzyaaaf12=xyzzyaabc12+xyzzyaaat12
do xyzzyaaah12=1,xyzzyaaar12
v(xyzzyaaae12)=v(xyzzyaabs12)-(v(xyzzyaabr12)+v(xyzzyaaaf12))
xyzzyaaae12=xyzzyaaae12+1
xyzzyaaaf12=xyzzyaaaf12+1
enddo
160 v(xyzzyaaae12)=v(xyzzyaabr12)-xyzzyaabk12*v(xyzzyaabs12)
xyzzyaaah12=1
170 iv(xyzzyaacg12)=xyzzyaaah12
xyzzyaaaw12=xyzzyaaay12+xyzzyaaah12
v(xyzzyaabo12)=x(xyzzyaaah12)
x(xyzzyaaah12)=x(xyzzyaaah12)+v(xyzzyaaaw12)
if(xyzzyaaah12==xyzzyaaaq12)x(xyzzyaaah12)=v(xyzzyaack12)-v(xyzzyaaaw1&
&2)
covirc=1
goto 390
180 x(xyzzyaaah12)=v(xyzzyaabo12)
if(iv(xyzzyaaci12)==0)goto 190
iv(xyzzyaabm12)=-2
goto 390
190 xyzzyaaaw12=xyzzyaaay12+xyzzyaaah12
xyzzyaaae12=xyzzyaabc12+xyzzyaaas12+xyzzyaaah12-1
xyzzyaaax12=xyzzyaaay12+xyzzyaaaq12
v(xyzzyaaae12)=(v(xyzzyaaae12)+v(xyzzyaabr12))/(v(xyzzyaaaw12)*v(xyzzy&
&aaax12))
xyzzyaaah12=xyzzyaaah12+1
if(xyzzyaaah12<=xyzzyaaaq12)goto 170
iv(xyzzyaacg12)=0
x(xyzzyaaaq12)=v(xyzzyaack12)
200 xyzzyaaaq12=xyzzyaaaq12+1
iv(xyzzyaacb12)=xyzzyaaaq12
if(xyzzyaaaq12>p)goto 210
xyzzyaabd12=v(xyzzyaabq12)*max(xyzzyaabj12/d(xyzzyaaaq12),abs(x(xyzzya&
&aaq12)))
if(x(xyzzyaaaq12)<xyzzyaabl12)xyzzyaabd12=-xyzzyaabd12
v(xyzzyaack12)=x(xyzzyaaaq12)
x(xyzzyaaaq12)=x(xyzzyaaaq12)+xyzzyaabd12
xyzzyaaax12=xyzzyaaay12+xyzzyaaaq12
v(xyzzyaaax12)=xyzzyaabd12
covirc=1
goto 390
210 xyzzyaaam12=iv(xyzzyaacf12)
call vcopy(n,r(1:n),v(xyzzyaaam12:))
v(xyzzyaabr12)=v(xyzzyaabs12)
if(xyzzyaaan12<0)goto 220
iv(xyzzyaacc12)=iv(xyzzyaach12)
xyzzyaaau12=iv(xyzzyaacd12)
call vcopy(n,v(xyzzyaaau12:xyzzyaaau12+n-1),r)
if(iv(xyzzyaabm12)<0)goto 390
covirc=3
goto 390
220 xyzzyaabc12=iv(xyzzyaaca12)
xyzzyaaad12=xyzzyaabc12
if(abs(xyzzyaaan12)==2)goto 230
xyzzyaaad12=abs(iv(xyzzyaabu12))
iv(xyzzyaabu12)=-xyzzyaaad12
230 call lsqrt(1,p,v(xyzzyaaad12:xyzzyaaad12+(p*(p+1)/2)-1),v(xyzzyaab&
&c12:),xyzzyaaal12)
iv(xyzzyaabm12)=-1
if(xyzzyaaal12/=0)goto 390
xyzzyaabb12=iv(xyzzyaacj12)+p
if(abs(xyzzyaaan12)>1)goto 340
call vscopy(p*(p+1)/2,v(xyzzyaabc12:xyzzyaabc12+p*(p+1)/2-1),xyzzyaabl&
&12)
xyzzyaabg12=iv(xyzzyaabz12)==(-1)
xyzzyaaaq12=p
if(xyzzyaabg12)xyzzyaaaq12=n
xyzzyaaba12=xyzzyaabb12-1
xyzzyaaav12=iv(xyzzyaace12)
do xyzzyaaah12=1,xyzzyaaaq12
if(xyzzyaabg12) goto 250
call vscopy(p,v(xyzzyaabb12:xyzzyaabb12+p-1),xyzzyaabl12)
xyzzyaaai12=xyzzyaabx12+xyzzyaaah12
xyzzyaaap12=xyzzyaaba12+iv(xyzzyaaai12)
v(xyzzyaaap12)=v(xyzzyaaav12)
xyzzyaaav12=xyzzyaaav12+1
if(xyzzyaaah12==p)goto 270
xyzzyaaak12=xyzzyaaah12+1
do xyzzyaaam12=xyzzyaaak12,p
xyzzyaaaj12=xyzzyaabx12+xyzzyaaam12
xyzzyaaap12=xyzzyaaba12+iv(xyzzyaaaj12)
v(xyzzyaaap12)=j(xyzzyaaah12,xyzzyaaam12)
enddo
goto 270
250 xyzzyaaap12=xyzzyaaba12
do xyzzyaaam12=1,p
xyzzyaaap12=xyzzyaaap12+1
v(xyzzyaaap12)=j(xyzzyaaah12,xyzzyaaam12)
enddo
270 call livmul(p,v(xyzzyaabb12:xyzzyaabb12+p-1),v(xyzzyaaad12:),v(xyz&
&zyaabb12:))
call litvmu(p,v(xyzzyaabb12:xyzzyaabb12+p-1),v(xyzzyaaad12:),v(xyzzyaa&
&bb12:))
xyzzyaaao12=xyzzyaabc12
do xyzzyaaam12=1,p
xyzzyaaap12=xyzzyaaba12+xyzzyaaam12
xyzzyaabf12=v(xyzzyaaap12)
do xyzzyaaap12=1,xyzzyaaam12
xyzzyaaaz12=xyzzyaaba12+xyzzyaaap12
v(xyzzyaaao12)=v(xyzzyaaao12)+xyzzyaabf12*v(xyzzyaaaz12)
xyzzyaaao12=xyzzyaaao12+1
enddo
enddo
enddo
goto 370
310 xyzzyaaav12=iv(xyzzyaace12)
if(iv(xyzzyaabz12)/=(-1))goto 320
xyzzyaaau12=iv(xyzzyaacd12)
call vcopy(n,v(xyzzyaaau12:xyzzyaaau12+n-1),r)
xyzzyaabb12=iv(xyzzyaacj12)+p
call qrfact(n,p,j,v(xyzzyaaav12:xyzzyaaav12+p-1),iv(xyzzyaabw12:xyzzya&
&abw12+p-1),iv(xyzzyaabv12),0,v(xyzzyaabb12:xyzzyaabb12+p-1))
iv(xyzzyaabz12)=-2
320 iv(xyzzyaabm12)=-1
if(iv(xyzzyaabv12)/=0)goto 390
xyzzyaabc12=iv(xyzzyaaca12)
xyzzyaaad12=abs(iv(xyzzyaabu12))
iv(xyzzyaabu12)=-xyzzyaaad12
xyzzyaaap12=xyzzyaaad12
do xyzzyaaah12=1,p
if(xyzzyaaah12>1)call vcopy(xyzzyaaah12-1,v(xyzzyaaap12:xyzzyaaap12+xy&
&zzyaaah12-1-1),j(:,xyzzyaaah12))
xyzzyaaap12=xyzzyaaap12+xyzzyaaah12-1
v(xyzzyaaap12)=v(xyzzyaaav12)
xyzzyaaap12=xyzzyaaap12+1
xyzzyaaav12=xyzzyaaav12+1
enddo
340 call linvrt(p,v(xyzzyaaad12:xyzzyaaad12+p*(p+1)/2-1),v(xyzzyaaad12&
&:))
call ltsqar(p,v(xyzzyaaad12:xyzzyaaad12+p*(p+1)/2-1),v(xyzzyaaad12:))
if(xyzzyaaad12==xyzzyaabc12)goto 370
do xyzzyaaah12=1,p
xyzzyaaaq12=xyzzyaabx12+xyzzyaaah12
xyzzyaaai12=iv(xyzzyaaaq12)
xyzzyaaao12=xyzzyaabc12-1+xyzzyaaai12*(xyzzyaaai12-1)/2
do xyzzyaaam12=1,xyzzyaaah12
xyzzyaaaq12=xyzzyaabx12+xyzzyaaam12
xyzzyaaaj12=iv(xyzzyaaaq12)
xyzzyaaap12=xyzzyaaao12+xyzzyaaaj12
if(xyzzyaaaj12>xyzzyaaai12)xyzzyaaap12=xyzzyaaap12+(xyzzyaaaj12-xyzzya&
&aai12)*(xyzzyaaaj12+xyzzyaaai12-3)/2
v(xyzzyaaap12)=v(xyzzyaaad12)
xyzzyaaad12=xyzzyaaad12+1
enddo
enddo
370 iv(xyzzyaabm12)=xyzzyaabc12
xyzzyaabe12=v(xyzzyaabr12)/(xyzzyaabh12*real(max(1,n-p)))
xyzzyaaam12=xyzzyaabc12-1+p*(p+1)/2
v(xyzzyaabc12:xyzzyaaam12)=xyzzyaabe12*v(xyzzyaabc12:xyzzyaaam12)
390 return
end subroutine covclc
subroutine dfault(iv,v)
integer,intent(inout) :: iv(:)
real(dp),intent(inout) :: v(:)
real(dp) :: xyzzyaaaa13,xyzzyaaab13,xyzzyaaac13
real(dp),parameter :: xyzzyaaad13=1.d0,xyzzyaaae13=3.d0
integer,parameter :: xyzzyaaaf13=14,xyzzyaaag13=15,xyzzyaaah13=16,xyzz&
&yaaai13=25,xyzzyaaaj13=17,xyzzyaaak13=18,xyzzyaaal13=19,xyzzyaaam13=2&
&0,xyzzyaaan13=21,xyzzyaaao13=22,xyzzyaaap13=23,xyzzyaaaq13=24
integer,parameter :: xyzzyaaar13=31,xyzzyaaas13=43,xyzzyaaat13=22,xyzz&
&yaaau13=44,xyzzyaaav13=41,xyzzyaaaw13=38,xyzzyaaax13=40,xyzzyaaay13=3&
&6,xyzzyaaaz13=37,xyzzyaaba13=19,xyzzyaabb13=45,xyzzyaabc13=23,xyzzyaa&
&bd13=39,xyzzyaabe13=35,xyzzyaabf13=20,xyzzyaabg13=21,xyzzyaabh13=24,x&
&yzzyaabi13=25,xyzzyaabj13=32,xyzzyaabk13=42,xyzzyaabl13=26,xyzzyaabm1&
&3=27,xyzzyaabn13=28,xyzzyaabo13=29,xyzzyaabp13=30,xyzzyaabq13=33,xyzz&
&yaabr13=34
iv(1)=12
iv(xyzzyaaaf13)=1
iv(xyzzyaaag13)=1
iv(xyzzyaaah13)=1
iv(xyzzyaaai13)=0
iv(xyzzyaaaj13)=200
iv(xyzzyaaak13)=150
iv(xyzzyaaal13)=1
iv(xyzzyaaam13)=1
iv(xyzzyaaan13)=imdcon(1)
iv(xyzzyaaao13)=1
iv(xyzzyaaap13)=1
iv(xyzzyaaaq13)=1
xyzzyaaaa13=rmdcon(3)
v(xyzzyaaar13)=1.d-20
if(xyzzyaaaa13>1.d-10)v(xyzzyaaar13)=xyzzyaaaa13**2
v(xyzzyaaas13)=max(1.d-6,1.d2*xyzzyaaaa13)
v(xyzzyaaat13)=0.5d0
xyzzyaaac13=rmdcon(4)
v(xyzzyaaau13)=xyzzyaaac13
v(xyzzyaaav13)=0.6d0
v(xyzzyaaaw13)=0.d0
xyzzyaaab13=xyzzyaaaa13**(xyzzyaaad13/xyzzyaaae13)
v(xyzzyaaax13)=xyzzyaaab13
v(xyzzyaaay13)=xyzzyaaac13
v(xyzzyaaaz13)=1.d0
v(xyzzyaaba13)=0.1d0
v(xyzzyaabb13)=1.5d0
v(xyzzyaabc13)=2.d0
v(xyzzyaabd13)=1.d-6
v(xyzzyaabe13)=100.d0
v(xyzzyaabf13)=-0.1d0
v(xyzzyaabg13)=0.1d0
v(xyzzyaabh13)=0.1d0
v(xyzzyaabi13)=4.d0
v(xyzzyaabj13)=max(1.d-10,xyzzyaaab13**2)
v(xyzzyaabk13)=rmdcon(5)
v(xyzzyaabl13)=0.1d0
v(xyzzyaabm13)=1.d-4
v(xyzzyaabn13)=0.75d0
v(xyzzyaabo13)=0.5d0
v(xyzzyaabp13)=0.75d0
v(xyzzyaabq13)=xyzzyaaac13
v(xyzzyaabr13)=1.d2*xyzzyaaaa13
return
end subroutine dfault
function dotprd(p,x,y) result(fn_val)
integer,intent(in) :: p
real(dp),intent(in) :: x(:),y(:)
real(dp) :: fn_val
integer xyzzyaaaa14
real(dp) xyzzyaaab14
real(dp),parameter :: xyzzyaaac14=1.d0,xyzzyaaad14=0.d0
real(dp),save :: xyzzyaaae14=0.d0
fn_val=xyzzyaaad14
if(p<=0)goto 30
if(xyzzyaaae14==xyzzyaaad14)xyzzyaaae14=rmdcon(2)
do xyzzyaaaa14=1,p
xyzzyaaab14=max(abs(x(xyzzyaaaa14)),abs(y(xyzzyaaaa14)))
if(xyzzyaaab14>xyzzyaaac14)goto 10
if(xyzzyaaab14<xyzzyaaae14)cycle
xyzzyaaab14=(x(xyzzyaaaa14)/xyzzyaaae14)*y(xyzzyaaaa14)
if(abs(xyzzyaaab14)<xyzzyaaae14)cycle
10 fn_val=fn_val+x(xyzzyaaaa14)*y(xyzzyaaaa14)
enddo
30 return
end function dotprd
subroutine dupdat(d,iv,j,n,p,v)
integer,intent(in)   :: iv(:),n,p
real(dp),intent(in)  :: j(:,:),v(:)
real(dp),intent(inout) :: d(:)
integer xyzzyaaaa15,xyzzyaaab15,xyzzyaaac15,xyzzyaaad15
real(dp) xyzzyaaae15,xyzzyaaaf15,xyzzyaaag15
real(dp) :: xyzzyaaah15=0.d0
integer,parameter :: xyzzyaaai15=41,xyzzyaaaj15=16,xyzzyaaak15=86,xyzz&
&yaaal15=31,xyzzyaaam15=53
xyzzyaaab15=iv(xyzzyaaaj15)
if(xyzzyaaab15==1)goto 10
if(iv(xyzzyaaal15)>0)goto 30
10 xyzzyaaag15=v(xyzzyaaai15)
xyzzyaaaa15=xyzzyaaak15+p
xyzzyaaad15=iv(xyzzyaaam15)-1
do xyzzyaaab15=1,p
xyzzyaaad15=xyzzyaaad15+xyzzyaaab15
xyzzyaaae15=v(xyzzyaaad15)
xyzzyaaaf15=v2norm(n,j(:,xyzzyaaab15))
if(xyzzyaaae15>xyzzyaaah15)xyzzyaaaf15=sqrt(xyzzyaaaf15*xyzzyaaaf15+xy&
&zzyaaae15)
xyzzyaaac15=xyzzyaaak15+xyzzyaaab15
xyzzyaaaa15=xyzzyaaaa15+1
if(xyzzyaaaf15<v(xyzzyaaac15))xyzzyaaaf15=max(v(xyzzyaaaa15),v(xyzzyaa&
&ac15))
d(xyzzyaaab15)=max(xyzzyaaag15*d(xyzzyaaab15),xyzzyaaaf15)
enddo
30 return
end subroutine dupdat
subroutine gqtstp(d,dig,dihdi,ka,l,p,step,v,w)
integer,intent(in) :: p
integer,intent(inout) :: ka
real(dp),intent(in) :: d(:),dig(:)
real(dp),intent(inout) :: l(:),v(:),step(:),w(:),dihdi(:)
integer xyzzyaaaa16,xyzzyaaab16,xyzzyaaac16,xyzzyaaad16,xyzzyaaae16,xy&
&zzyaaaf16,xyzzyaaag16,xyzzyaaah16,xyzzyaaai16,xyzzyaaaj16,xyzzyaaak16&
&,xyzzyaaal16,xyzzyaaam16,xyzzyaaan16,xyzzyaaao16, xyzzyaaap16,xyzzyaa&
&aq16,xyzzyaaar16,xyzzyaaas16,xyzzyaaat16,xyzzyaaau16
real(dp) xyzzyaaav16,xyzzyaaaw16,xyzzyaaax16,xyzzyaaay16,xyzzyaaaz16,x&
&yzzyaaba16,xyzzyaabb16,xyzzyaabc16,xyzzyaabd16,xyzzyaabe16,xyzzyaabf1&
&6,xyzzyaabg16,xyzzyaabh16,xyzzyaabi16,xyzzyaabj16,xyzzyaabk16,xyzzyaa&
&bl16,xyzzyaabm16,xyzzyaabn16,xyzzyaabo16,xyzzyaabp16,xyzzyaabq16
logical xyzzyaabr16
integer,parameter :: xyzzyaabs16=1,xyzzyaabt16=2,xyzzyaabu16=3,xyzzyaa&
&bv16=19,xyzzyaabw16=4,   xyzzyaabx16=6,xyzzyaaby16=20,xyzzyaabz16=21,&
&xyzzyaaca16=7,xyzzyaacb16=8,xyzzyaacc16=9,xyzzyaacd16=5
real(dp),parameter :: xyzzyaace16=50.d0,xyzzyaacf16=4.d0,xyzzyaacg16=0&
&.5d0,xyzzyaach16=2.d0, xyzzyaaci16=-1.d0,xyzzyaacj16=1.d0,xyzzyaack16&
&=1.d-3,xyzzyaacl16=6.d0,xyzzyaacm16=3.d0,xyzzyaacn16=2.d0,xyzzyaaco16&
&=0.d0
real(dp),save :: xyzzyaacp16=0.d0
xyzzyaaaz16=0.d0
xyzzyaabd16=0.d0
xyzzyaaav16=0.d0
xyzzyaabp16=0.d0
xyzzyaaaa16=p+1
xyzzyaaae16=xyzzyaaaa16+1
xyzzyaaaf16=xyzzyaaae16+1
xyzzyaaao16=xyzzyaaaf16+1
xyzzyaaap16=xyzzyaaao16+1
xyzzyaaas16=xyzzyaaap16+1
xyzzyaaad16=xyzzyaaas16+1
xyzzyaaac16=xyzzyaaad16
xyzzyaaab16=xyzzyaaac16+1
xyzzyaaar16=xyzzyaaac16+p
xyzzyaaaq16=xyzzyaaar16+1
xyzzyaabh16=v(xyzzyaacb16)
xyzzyaabe16=v(xyzzyaabz16)*xyzzyaabh16
xyzzyaabf16=v(xyzzyaaby16)*xyzzyaabh16
xyzzyaabg16=xyzzyaacn16*v(xyzzyaabv16)/(xyzzyaacm16*(xyzzyaacf16*(v(xy&
&zzyaaby16)+xyzzyaacj16)*(xyzzyaach16+xyzzyaacj16)+  xyzzyaach16+xyzzy&
&aacn16)*xyzzyaabh16**2)
xyzzyaabc16=xyzzyaaco16
xyzzyaaba16=v(xyzzyaabv16)/xyzzyaacl16
xyzzyaaaj16=0
xyzzyaabr16=.false.
xyzzyaaam16=ka+50
if(ka>=0)goto 290
xyzzyaaal16=0
xyzzyaabp16=xyzzyaaci16
ka=0
xyzzyaaam16=50
xyzzyaaak16=0
do xyzzyaaag16=1,p
xyzzyaaak16=xyzzyaaak16+xyzzyaaag16
xyzzyaaan16=xyzzyaaac16+xyzzyaaag16
w(xyzzyaaan16)=dihdi(xyzzyaaak16)
enddo
xyzzyaabo16=xyzzyaaco16
xyzzyaaak16=p*(p+1)/2
do xyzzyaaag16=1,xyzzyaaak16
xyzzyaabm16=abs(dihdi(xyzzyaaag16))
if(xyzzyaabo16<xyzzyaabm16)xyzzyaabo16=xyzzyaabm16
enddo
w(xyzzyaaaa16)=xyzzyaabo16
30 call lsqrt (1, p, l(1:(p*(p+1)/2)), dihdi, xyzzyaaaj16)
if(xyzzyaaaj16==0)goto 50
xyzzyaaak16=xyzzyaaaj16*(xyzzyaaaj16+1)/2
xyzzyaabm16=l(xyzzyaaak16)
l(xyzzyaaak16)=xyzzyaacj16
w(1:xyzzyaaaj16) = xyzzyaaco16
w(xyzzyaaaj16)=xyzzyaacj16
call litvmu(xyzzyaaaj16,w(1:xyzzyaaaj16),l,w)
xyzzyaabo16=v2norm(xyzzyaaaj16,w)
xyzzyaabb16=-xyzzyaabm16/xyzzyaabo16/xyzzyaabo16
v(xyzzyaabu16)=-xyzzyaabb16
if(xyzzyaabr16)goto 200
v(xyzzyaabx16)=xyzzyaaco16
goto 60
50 xyzzyaabb16=xyzzyaaco16
call livmul(p,w(xyzzyaaaq16:xyzzyaaaq16+p-1),l,dig)
v(xyzzyaabx16)=xyzzyaacg16*dotprd(p,w(xyzzyaaaq16:),w(xyzzyaaaq16:))
call litvmu(p,w(xyzzyaaaq16:xyzzyaaaq16+p-1),l,w(xyzzyaaaq16:))
xyzzyaaaz16=v2norm(p,w(xyzzyaaaq16:))
v(xyzzyaabu16)=xyzzyaaaz16
xyzzyaabd16=xyzzyaaaz16-xyzzyaabh16
if(xyzzyaabd16<=xyzzyaabe16)goto 260
if(xyzzyaabr16)goto 200
60 v(xyzzyaabs16)=v2norm(p,dig)
if(v(xyzzyaabs16)==xyzzyaaco16)goto 430
xyzzyaaal16=0
do xyzzyaaag16=1,p
xyzzyaabq16=xyzzyaaco16
if(xyzzyaaag16==1)goto 80
xyzzyaaah16=xyzzyaaag16-1
do xyzzyaaak16=1,xyzzyaaah16
xyzzyaaal16=xyzzyaaal16+1
xyzzyaabm16=abs(dihdi(xyzzyaaal16))
xyzzyaabq16=xyzzyaabq16+xyzzyaabm16
w(xyzzyaaak16)=w(xyzzyaaak16)+xyzzyaabm16
enddo
80 w(xyzzyaaag16)=xyzzyaabq16
xyzzyaaal16=xyzzyaaal16+1
enddo
xyzzyaaal16=1
xyzzyaabo16=w(xyzzyaaab16)-w(1)
if(p<=1)goto 110
do xyzzyaaag16=2,p
xyzzyaaak16=xyzzyaaac16+xyzzyaaag16
xyzzyaabm16=w(xyzzyaaak16)-w(xyzzyaaag16)
if(xyzzyaabm16>=xyzzyaabo16)cycle
xyzzyaabo16=xyzzyaabm16
xyzzyaaal16=xyzzyaaag16
enddo
110 xyzzyaabk16=w(xyzzyaaal16)
xyzzyaaak16=xyzzyaaac16+xyzzyaaal16
xyzzyaaax16=w(xyzzyaaak16)
xyzzyaaan16=xyzzyaaal16*(xyzzyaaal16-1)/2+1
xyzzyaaai16=1
xyzzyaabm16=xyzzyaaco16
do xyzzyaaag16=1,p
if(xyzzyaaag16==xyzzyaaal16)goto 120
xyzzyaaaw16=abs(dihdi(xyzzyaaan16))
xyzzyaabj16=w(xyzzyaaag16)
xyzzyaaak16=xyzzyaaac16+xyzzyaaag16
xyzzyaabo16=xyzzyaacg16*(xyzzyaaax16-w(xyzzyaaak16)+xyzzyaabj16-xyzzya&
&aaw16)
xyzzyaabo16=xyzzyaabo16+sqrt(xyzzyaabo16*xyzzyaabo16+xyzzyaabk16*xyzzy&
&aaaw16)
if(xyzzyaabm16<xyzzyaabo16)xyzzyaabm16=xyzzyaabo16
if(xyzzyaaag16<xyzzyaaal16)goto 130
120 xyzzyaaai16=xyzzyaaag16
130 xyzzyaaan16=xyzzyaaan16+xyzzyaaai16
enddo
w(xyzzyaaaf16)=xyzzyaaax16-xyzzyaabm16
xyzzyaabp16=v(xyzzyaabs16)/xyzzyaabh16-w(xyzzyaaaf16)
xyzzyaaal16=1
xyzzyaabo16=w(xyzzyaaab16)+w(1)
if(p<=1)goto 160
do xyzzyaaag16=2,p
xyzzyaaak16=xyzzyaaac16+xyzzyaaag16
xyzzyaabm16=w(xyzzyaaak16)+w(xyzzyaaag16)
if(xyzzyaabm16<=xyzzyaabo16)cycle
xyzzyaabo16=xyzzyaabm16
xyzzyaaal16=xyzzyaaag16
enddo
160 xyzzyaabk16=w(xyzzyaaal16)
xyzzyaaak16=xyzzyaaac16+xyzzyaaal16
xyzzyaaax16=w(xyzzyaaak16)
xyzzyaaan16=xyzzyaaal16*(xyzzyaaal16-1)/2+1
xyzzyaaai16=1
xyzzyaabm16=xyzzyaaco16
do xyzzyaaag16=1,p
if(xyzzyaaag16==xyzzyaaal16)goto 170
xyzzyaaaw16=abs(dihdi(xyzzyaaan16))
xyzzyaabj16=w(xyzzyaaag16)
xyzzyaaak16=xyzzyaaac16+xyzzyaaag16
xyzzyaabo16=xyzzyaacg16*(w(xyzzyaaak16)+xyzzyaabj16-xyzzyaaaw16-xyzzya&
&aax16)
xyzzyaabo16=xyzzyaabo16+sqrt(xyzzyaabo16*xyzzyaabo16+xyzzyaabk16*xyzzy&
&aaaw16)
if(xyzzyaabm16<xyzzyaabo16)xyzzyaabm16=xyzzyaabo16
if(xyzzyaaag16<xyzzyaaal16)goto 180
170 xyzzyaaai16=xyzzyaaag16
180 xyzzyaaan16=xyzzyaaan16+xyzzyaaai16
enddo
w(xyzzyaaae16)=xyzzyaaax16+xyzzyaabm16
xyzzyaabb16=max(xyzzyaabb16,v(xyzzyaabs16)/xyzzyaabh16-w(xyzzyaaae16))
xyzzyaaav16=abs(v(xyzzyaacd16))*v(xyzzyaacc16)/xyzzyaabh16
if(xyzzyaaaj16/=0)goto 200
call livmul(p,w(1:p),l,w(xyzzyaaaq16:))
xyzzyaabm16=v2norm(p,w)
w(xyzzyaaap16)=xyzzyaaaz16/xyzzyaabm16/xyzzyaabm16
xyzzyaabb16=max(xyzzyaabb16,xyzzyaabd16*w(xyzzyaaap16))
200 ka=ka+1
if(-v(xyzzyaabu16)>=xyzzyaaav16.or.xyzzyaaav16<xyzzyaabb16.or.xyzzyaaa&
&v16>=xyzzyaabp16)xyzzyaaav16=xyzzyaabp16* max(xyzzyaack16,sqrt(xyzzya&
&abb16/xyzzyaabp16))
xyzzyaaal16=0
do xyzzyaaag16=1,p
xyzzyaaal16=xyzzyaaal16+xyzzyaaag16
xyzzyaaak16=xyzzyaaac16+xyzzyaaag16
dihdi(xyzzyaaal16)=w(xyzzyaaak16)+xyzzyaaav16
enddo
call lsqrt(1,p,l(1:p*(p+1)/2),dihdi,xyzzyaaaj16)
if(xyzzyaaaj16==0)goto 230
xyzzyaaak16=(xyzzyaaaj16*(xyzzyaaaj16+1))/2
xyzzyaabm16=l(xyzzyaaak16)
l(xyzzyaaak16)=xyzzyaacj16
do xyzzyaaag16=1,xyzzyaaaj16
w(xyzzyaaag16)=xyzzyaaco16
enddo
w(xyzzyaaaj16)=xyzzyaacj16
call litvmu(xyzzyaaaj16,w(1:xyzzyaaaj16),l,w)
xyzzyaabo16=v2norm(xyzzyaaaj16,w)
xyzzyaabb16=xyzzyaaav16-xyzzyaabm16/xyzzyaabo16/xyzzyaabo16
v(xyzzyaabu16)=-xyzzyaabb16
goto 200
230 call livmul(p,w(xyzzyaaaq16:xyzzyaaaq16+p-1),l,dig)
call litvmu(p,w(xyzzyaaaq16:xyzzyaaaq16+p-1),l,w(xyzzyaaaq16:))
xyzzyaaaz16=v2norm(p,w(xyzzyaaaq16:))
xyzzyaabd16=xyzzyaaaz16-xyzzyaabh16
if(xyzzyaabd16<=xyzzyaabe16.and.xyzzyaabd16>=xyzzyaabf16)goto 270
if(xyzzyaabd16==xyzzyaabc16)goto 270
xyzzyaabc16=xyzzyaabd16
if(xyzzyaabd16>xyzzyaaco16)goto 240
if(v(xyzzyaabu16)>xyzzyaaco16)goto 240
xyzzyaaay16=xyzzyaaav16+v(xyzzyaabu16)
xyzzyaabn16=xyzzyaaav16*xyzzyaaaz16*xyzzyaaaz16+dotprd(p,dig,w(xyzzyaa&
&aq16:))
if(xyzzyaaay16<xyzzyaabg16*xyzzyaabn16)goto 250
240 if(ka>=xyzzyaaam16)goto 270
call livmul(p,w(1:p),l,w(xyzzyaaaq16:))
xyzzyaabo16=v2norm(p,w)
if(xyzzyaabd16<xyzzyaaco16)xyzzyaabp16=min(xyzzyaabp16,xyzzyaaav16)
xyzzyaaav16=xyzzyaaav16+(xyzzyaabd16/xyzzyaabo16)*(xyzzyaaaz16/xyzzyaa&
&bo16)*(xyzzyaaaz16/xyzzyaabh16)
xyzzyaabb16=max(xyzzyaabb16,xyzzyaaav16)
goto 200
250 if(xyzzyaacp16==xyzzyaaco16)xyzzyaacp16=xyzzyaace16*rmdcon(3)
if(xyzzyaaay16>xyzzyaacp16*w(xyzzyaaaa16))goto 330
goto 270
260 xyzzyaaav16=xyzzyaaco16
270 step(1:p)=-w(xyzzyaaar16+1:xyzzyaaar16+p)/d(1:p)
v(xyzzyaabw16)=-dotprd(p,dig,w(xyzzyaaaq16:))
v(xyzzyaaca16)=xyzzyaacg16*(abs(xyzzyaaav16)*xyzzyaaaz16*xyzzyaaaz16-v&
&(xyzzyaabw16))
goto 410
290 if(v(xyzzyaabu16)<=xyzzyaaco16.or.v(xyzzyaabu16)-xyzzyaabh16>xyzzy&
&aabe16)goto 310
xyzzyaabr16=.true.
ka=ka+1
xyzzyaaal16=0
do xyzzyaaag16=1,p
xyzzyaaal16=xyzzyaaal16+xyzzyaaag16
xyzzyaaak16=xyzzyaaac16+xyzzyaaag16
dihdi(xyzzyaaal16)=w(xyzzyaaak16)
enddo
xyzzyaabp16=xyzzyaaci16
goto 30
310 if(ka==0)goto 50
xyzzyaaaz16=w(xyzzyaaad16)
xyzzyaaav16=abs(v(xyzzyaacd16))
xyzzyaabd16=xyzzyaaaz16-xyzzyaabh16
xyzzyaabm16=v(xyzzyaabs16)/xyzzyaabh16
if(xyzzyaabh16>v(xyzzyaacc16))goto 320
xyzzyaabp16=xyzzyaabm16-w(xyzzyaaaf16)
xyzzyaabb16=xyzzyaaco16
if(xyzzyaaav16>xyzzyaaco16)xyzzyaabb16=w(xyzzyaaao16)
xyzzyaabb16=max(xyzzyaabb16,xyzzyaabm16-w(xyzzyaaae16))
if(v(xyzzyaabu16)>xyzzyaaco16)xyzzyaabb16=max(xyzzyaabb16,(v(xyzzyaabu&
&16)-xyzzyaabh16)*w(xyzzyaaap16))
goto 240
320 xyzzyaabp16=xyzzyaabm16-w(xyzzyaaaf16)
if(xyzzyaaav16>xyzzyaaco16)xyzzyaabp16=min(xyzzyaabp16,w(xyzzyaaas16))
xyzzyaabb16=max(xyzzyaaco16,-v(xyzzyaabu16),xyzzyaabm16-w(xyzzyaaae16)&
&)
if(v(xyzzyaabu16)>xyzzyaaco16)xyzzyaabb16=max(xyzzyaabb16,(v(xyzzyaabu&
&16)-xyzzyaabh16)*w(xyzzyaaap16))
goto 240
330 xyzzyaaav16=-xyzzyaaav16
xyzzyaaau16=xyzzyaaar16+p
xyzzyaaat16=xyzzyaaau16+1
xyzzyaaay16=xyzzyaach16*xyzzyaaay16
call lsvmin(p,l,w(xyzzyaaat16:xyzzyaaat16+p-1),w(1:p),xyzzyaabm16)
xyzzyaaal16=0
340 w(1:p)=xyzzyaabm16*w(1:p)
call litvmu(p,w(1:p),l,w)
xyzzyaabo16=xyzzyaacj16/v2norm(p,w)
xyzzyaabm16=xyzzyaabo16*xyzzyaabm16
if(xyzzyaabm16<=xyzzyaaay16)goto 370
if(xyzzyaaal16>30)goto 270
xyzzyaaal16=xyzzyaaal16+1
w(xyzzyaaau16+1:xyzzyaaau16+p)=xyzzyaabo16*w(1:p)
call livmul(p,w(1:p),l,w(xyzzyaaat16:))
xyzzyaabm16=xyzzyaacj16/v2norm(p,w)
goto 340
370 w(1:p)=xyzzyaabo16*w(1:p)
xyzzyaabl16=dotprd(p,w(xyzzyaaaq16:),w)
xyzzyaabo16=(xyzzyaabh16+xyzzyaaaz16)*(xyzzyaabh16-xyzzyaaaz16)
xyzzyaabi16=sqrt(xyzzyaabl16*xyzzyaabl16+xyzzyaabo16)
if(xyzzyaabl16<xyzzyaaco16)xyzzyaabi16=-xyzzyaabi16
xyzzyaabj16=xyzzyaabo16/(xyzzyaabl16+xyzzyaabi16)
v(xyzzyaaca16)=xyzzyaacg16*xyzzyaabn16
xyzzyaabo16=xyzzyaaco16
xyzzyaabm16=xyzzyaabj16*(xyzzyaaav16*xyzzyaabl16-xyzzyaacg16*xyzzyaabj&
&16*(xyzzyaaav16+xyzzyaabm16*dotprd(p,w(xyzzyaaat16:),w)))
if(xyzzyaabm16<xyzzyaaba16*xyzzyaabn16)goto 390
v(xyzzyaaca16)=v(xyzzyaaca16)+xyzzyaabm16
xyzzyaaaz16=xyzzyaabh16
xyzzyaabo16=-xyzzyaabj16
390 do xyzzyaaag16=1,p
xyzzyaaak16=xyzzyaaar16+xyzzyaaag16
w(xyzzyaaak16)=xyzzyaabo16*w(xyzzyaaag16)-w(xyzzyaaak16)
step(xyzzyaaag16)=w(xyzzyaaak16)/d(xyzzyaaag16)
enddo
v(xyzzyaabw16)=dotprd(p,dig,w(xyzzyaaaq16:))
410 v(xyzzyaabt16)=xyzzyaaaz16
v(xyzzyaacd16)=xyzzyaaav16
w(xyzzyaaao16)=xyzzyaabb16
w(xyzzyaaas16)=xyzzyaabp16
v(xyzzyaacc16)=xyzzyaabh16
w(xyzzyaaad16)=xyzzyaaaz16
xyzzyaaak16=0
do xyzzyaaag16=1,p
xyzzyaaak16=xyzzyaaak16+xyzzyaaag16
xyzzyaaal16=xyzzyaaac16+xyzzyaaag16
dihdi(xyzzyaaak16)=w(xyzzyaaal16)
enddo
goto 450
430 v(xyzzyaacd16)=xyzzyaaco16
v(xyzzyaaca16)=xyzzyaaco16
v(xyzzyaabt16)=xyzzyaaco16
v(xyzzyaabw16)=xyzzyaaco16
step(1:p)=xyzzyaaco16
450 return
end subroutine gqtstp
subroutine itsmry (d,iv,p,v,x)
integer,intent(in) :: p
integer,intent(inout) :: iv(:)
real(dp),intent(in) :: d(:),v(:),x(:)
integer xyzzyaaaa17,xyzzyaaab17,xyzzyaaac17,xyzzyaaad17,xyzzyaaae17,xy&
&zzyaaaf17,xyzzyaaag17,nf,xyzzyaaah17,xyzzyaaai17,pu
real(dp) :: xyzzyaaaj17,xyzzyaaak17,xyzzyaaal17,reldf
integer, parameter :: xyzzyaaam17=26,xyzzyaaan17=14,xyzzyaaao17=28,xyz&
&zyaaap17=15,xyzzyaaaq17=39,    xyzzyaaar17=6,xyzzyaaas17=40,xyzzyaaat&
&17=41,xyzzyaaau17=30,xyzzyaaav17=31,  xyzzyaaaw17=19,xyzzyaaax17=48,x&
&yzzyaaay17=21,xyzzyaaaz17=22,        xyzzyaaba17=23,xyzzyaabb17=57,xy&
&zzyaabc17=24
integer, parameter :: xyzzyaabd17=2,f=10,xyzzyaabe17=13,xyzzyaabf17=11&
&,xyzzyaabg17=6,xyzzyaabh17=7,   xyzzyaabi17=17,xyzzyaabj17=47,xyzzyaa&
&bk17=5
real(dp),parameter :: xyzzyaabl17=0.d0
pu=iv(xyzzyaaay17)
if(pu<0)goto 610
xyzzyaaae17=iv(1)
xyzzyaaai17=iv(xyzzyaaaw17)
if(xyzzyaaae17<2.or.xyzzyaaae17>15)goto 320
if(xyzzyaaai17==0)goto 70
if(xyzzyaaae17>=12)goto 70
if(xyzzyaaae17>=10.and.iv(xyzzyaaax17)==0)goto 70
if(xyzzyaaae17>2)goto 10
iv(xyzzyaaax17)=iv(xyzzyaaax17)+1
if(iv(xyzzyaaax17)<abs(xyzzyaaai17))goto 610
10 nf=iv(xyzzyaaar17)-abs(iv(xyzzyaaas17))
iv(xyzzyaaax17)=0
reldf=xyzzyaabl17
xyzzyaaal17=xyzzyaabl17
xyzzyaaak17=v(xyzzyaabe17)
if(xyzzyaaak17<=xyzzyaabl17)goto 20
reldf=v(xyzzyaabf17)/xyzzyaaak17
xyzzyaaal17=v(xyzzyaabh17)/xyzzyaaak17
20 if(xyzzyaaai17>0)goto 40
iv(xyzzyaaaq17)=0
write(pu,60)iv(xyzzyaaav17),nf,v(f),reldf,xyzzyaaal17,v(xyzzyaabi17)
goto 70
40 continue
iv(xyzzyaaaq17)=0
xyzzyaaag17=iv(xyzzyaabb17)
xyzzyaaaj17=xyzzyaabl17
if(xyzzyaaak17>xyzzyaabl17)xyzzyaaaj17=v(xyzzyaabg17)/xyzzyaaak17
write(pu,*)'NL2SOL iteration summary :'
write(pu,*)'--------------------------'
write(pu,'(a,1x,i16)')' Iteration                  :',iv(xyzzyaaav17)
write(pu,'(a,1x,i16)')' No. function evals.        :',nf
write(pu,'(a,1x,e16.8)')' Function value             :',v(f)
write(pu,'(a,1x,e16.8)')' Relative function diff.    :',reldf
write(pu,'(a,1x,e16.8)')' Predicted rel. func. diff. :',xyzzyaaal17
write(pu,'(a,1x,e16.8)')' Rel. parameter change      :',v(xyzzyaabi17)
select case(xyzzyaaag17)
case(1)
write(pu,'(a)')' Quadratic model            :   Gauss-Newton'
case(2)
write(pu,'(a)')' Quadratic model            :   Augmented'
case(3)
write(pu,'(a)')' Quadratic models           :   Gauss-Newton + Augment&
&ed'
case(4)
write(pu,'(a)')' Quadratic models           :   Augmented + Gauss-Newt&
&on'
case(5)
write(pu,'(a)')' Quadratic models           :   Gauss-Newton + Augment&
&ed + Gauss-Newton'
case(6)
write(pu,'(a)')' Quadratic models           :   Augmented + Gauss-Newt&
&on + Augmented'
end select
write(pu,'(a,1x,e16.8)')' Marquardt parameter        :',v(xyzzyaabk17)
write(pu,'(a,1x,e16.8)')' Sizing factor              :',v(xyzzyaabj17)
write(pu,'(a,1x,e16.8)')' Norm of d times step       :',v(xyzzyaabd17)
write(pu,'(a,1x,e16.8)')' Convergence parameter      :',xyzzyaaaj17
write(pu,*)
60 format(' ', i5, i6, 4f11.3, a3, a4, 4f11.3)
70 select case (xyzzyaaae17)
case (1:2)
goto 610
case (3)
write(pu,90)
goto 370
90 format(' NL2SOL: parameter convergence.'/)
case (4)
write(pu,110)
goto 370
110 format(' NL2SOL: relative function convergence.'/)
case (5)
write(pu,130)
goto 370
130 format(' NL2SOL: parameter and relative function convergence.'/)
case (6)
write(pu,150)
goto 370
150 format(' NL2SOL: absolute function convergence.'/)
case (7)
write(pu,170)
goto 370
170 format(' NL2SOL: singular convergence.'/)
case (8)
write(pu,190)
goto 370
190 format(' NL2SOL: false convergence.'/)
case (9)
write(pu,210)
goto 370
210 format(' NL2SOL: function evaluation limit.'/)
case (10)
write(pu,230)
goto 370
230 format(' NL2SOL: iteration limit.'/)
case (11)
write(pu,250)
goto 370
250 format(' NL2SOL: halting due to stopx.'/)
case (12)
goto 340
case (13)
write(pu,270)
goto 340
270 format(' NL2SOL: initial sum of squares overflows.'/)
case (14)
write(pu,290)
goto 610
290 format(' NL2SOL: bad parameters to assess.'/)
case (15)
write(pu,310)
310 format(' NL2SOL: Jacobian could not be computed.'/)
if(iv(xyzzyaaav17)>0)goto 420
goto 340
end select
320 write(pu,330)xyzzyaaae17
330 format(' NL2SOL: iv(1) is ',i5/)
goto 610
340 if(iv(xyzzyaabc17)/=0)then
write(pu,*)'NL2SOL initial parameters :'
write(pu,*)'---------------------------'
write(pu,*)'    i              x(i)           d(i)'
write(pu,*)'----- ----------------- --------------'
write(pu,350)(xyzzyaaac17,x(xyzzyaaac17),d(xyzzyaaac17),xyzzyaaac17=1,&
&p)
write(pu,*)'--------------------------------------'
endif
350 format(1x,i5,1x,f17.6,1x,f14.3)
if(xyzzyaaae17>=13)goto 610
iv(xyzzyaaaq17)=0
iv(xyzzyaaax17)=0
if(xyzzyaaai17==0)goto 610
write(pu,360)v(f)
360 format(' Initial function value: ',e16.8/)
goto 610
370 iv(xyzzyaaaq17)=1
if(iv(xyzzyaaba17)==0)goto 420
xyzzyaaak17=v(xyzzyaabe17)
xyzzyaaal17=xyzzyaabl17
xyzzyaaaj17=xyzzyaabl17
if(xyzzyaaak17<=xyzzyaabl17)goto 380
xyzzyaaal17=v(xyzzyaabh17)/xyzzyaaak17
xyzzyaaaj17=v(xyzzyaabg17)/xyzzyaaak17
380 nf=iv(xyzzyaaar17)-iv(xyzzyaaas17)
xyzzyaaah17=iv(xyzzyaaau17)-iv(xyzzyaaat17)
write(pu,*)'NL2SOL final summary :'
write(pu,*)'----------------------'
write(pu,'(a,1x,e16.8)')' Final function value       :',v(f)
write(pu,'(a,1x,e16.8)')' Rel. parameter change      :',v(xyzzyaabi17)
write(pu,'(a,1x,i16)')' Number of function evals.  :',nf
write(pu,'(a,1x,i16)')' Number of gradient evals.  :',xyzzyaaah17
write(pu,'(a,1x,e16.8)')' Predicted rel. func. diff. :',xyzzyaaal17
write(pu,'(a,1x,e16.8)')' Convergence parameter      :',xyzzyaaaj17
if(iv(xyzzyaaas17)>0)write(pu,400)iv(xyzzyaaas17)
400 format (' Func. evals for covariance : ',i16)
if(iv(xyzzyaaat17)>0)write(pu,410)iv(xyzzyaaat17)
410 format (' Grad. evals for covariance : ',i16)
write(pu,*)
420 if(iv(xyzzyaaaz17)==0)goto 460
iv(xyzzyaaaq17)=1
xyzzyaaab17=iv(xyzzyaaao17)
write(pu,*)'Final parameters :'
write(pu,*)'------------------'
write(pu,*)'    i              x(i)           d(i)           g(i)'
write(pu,*)'----- ----------------- -------------- --------------'
do xyzzyaaac17=1,p
write(pu,450)xyzzyaaac17,x(xyzzyaaac17),d(xyzzyaaac17),v(xyzzyaaab17)
xyzzyaaab17=xyzzyaaab17+1
enddo
write(pu,*)'-----------------------------------------------------'
write(pu,*)
450 format(1x,i5,1x,f17.6,2(1x,f14.3))
460 if(iv(xyzzyaaan17)==0)goto 610
xyzzyaaaa17=iv(xyzzyaaam17)
iv(xyzzyaaaq17)=1
if(xyzzyaaaa17<0.0)then
goto 470
elseif(xyzzyaaaa17==0.0)then
goto 500
else
goto 520
endif
470 if(-1==xyzzyaaaa17)write(pu,480)
480 format(' NL2SOL: indefinite covariance matrix.'/)
if(-2==xyzzyaaaa17)write(pu,490)
490 format(' NL2SOL: oversize steps in computing covariance.'/)
goto 610
500 write(pu,510)
510 format(' NL2SOL: covariance matrix not computed'/)
goto 610
520 xyzzyaaac17=abs(iv(xyzzyaaap17))
write(pu,*)'NL2SOL print-out of covariance matrix :'
write(pu,*)'---------------------------------------'
if(xyzzyaaac17<=1)write(pu,530)
530 format(' Cov. matrix = scale * H**-1 * (J**t * J) * H**-1'/)
if(xyzzyaaac17==2)write(pu,540)
540 format (' Cov. matrix = scale * H**-1'/)
if(xyzzyaaac17>=3)write(pu,550)
550 format(' Cov. matrix = scale * (J**t * J)**-1'/)
xyzzyaaad17=xyzzyaaaa17-1
if(xyzzyaaai17<=0)goto 580
do xyzzyaaac17=1,p
xyzzyaaaf17=xyzzyaaad17+1
xyzzyaaad17=xyzzyaaad17+xyzzyaaac17
write(pu,'(a,i3,a)')'Row ',xyzzyaaac17,' is:'
write(pu,570)v(xyzzyaaaf17:xyzzyaaad17)
enddo
write(pu,*)
570 format(4(1x,e16.8))
goto 610
580 do xyzzyaaac17=1,p
xyzzyaaaf17=xyzzyaaad17+1
xyzzyaaad17=xyzzyaaad17+xyzzyaaac17
write(pu,'(a,i3,a)')'Row ',xyzzyaaac17,' is:'
write(pu,600)v(xyzzyaaaf17:xyzzyaaad17)
enddo
write(pu,*)
600 format(4(1x,e16.8))
610 return
end subroutine itsmry
subroutine linvrt(n,lin,l)
integer,intent(in) :: n
real(dp),intent(in) :: l(:)
real(dp),intent(inout) :: lin(:)
integer xyzzyaaaa18,xyzzyaaab18,xyzzyaaac18,xyzzyaaad18,xyzzyaaae18,xy&
&zzyaaaf18,xyzzyaaag18,xyzzyaaah18,xyzzyaaai18
real(dp) xyzzyaaaj18
real(dp),parameter :: xyzzyaaak18=1.d0,xyzzyaaal18=0.d0
xyzzyaaai18=n+1
xyzzyaaae18=n*(xyzzyaaai18)/2
do xyzzyaaab18=1,n
xyzzyaaaa18=xyzzyaaai18-xyzzyaaab18
lin(xyzzyaaae18)=xyzzyaaak18/l(xyzzyaaae18)
if(xyzzyaaaa18 <= 1) goto 40
xyzzyaaaf18=xyzzyaaae18
xyzzyaaac18=xyzzyaaaa18-1
do xyzzyaaad18=1,xyzzyaaac18
xyzzyaaaj18=xyzzyaaal18
xyzzyaaae18=xyzzyaaaf18
xyzzyaaah18=xyzzyaaaf18-xyzzyaaad18
do xyzzyaaag18=1,xyzzyaaad18
xyzzyaaaj18=xyzzyaaaj18-l(xyzzyaaah18)*lin(xyzzyaaae18)
xyzzyaaae18=xyzzyaaae18-1
xyzzyaaah18=xyzzyaaah18+xyzzyaaag18-xyzzyaaaa18
enddo
lin(xyzzyaaae18)=xyzzyaaaj18/l(xyzzyaaah18)
enddo
xyzzyaaae18=xyzzyaaae18-1
enddo
40 return
end subroutine linvrt
subroutine litvmu (n, x, l, y)
integer,intent(in) :: n
real(dp),intent(in) :: l(:),y(:)
real(dp),intent(inout) :: x(:)
integer xyzzyaaaa19,xyzzyaaab19,xyzzyaaac19,xyzzyaaad19,xyzzyaaae19,xy&
&zzyaaaf19,xyzzyaaag19
real(dp) xyzzyaaah19
real(dp),parameter :: xyzzyaaai19=0.d0
x(1:n)=y(1:n)
xyzzyaaag19=n+1
xyzzyaaae19=n*(n+1)/2
do xyzzyaaab19=1,n
xyzzyaaaa19=xyzzyaaag19-xyzzyaaab19
xyzzyaaah19=x(xyzzyaaaa19)/l(xyzzyaaae19)
x(xyzzyaaaa19)=xyzzyaaah19
if(xyzzyaaaa19<=1)goto 40
xyzzyaaae19=xyzzyaaae19-xyzzyaaaa19
if(xyzzyaaah19==xyzzyaaai19)cycle
xyzzyaaad19=xyzzyaaaa19-1
do xyzzyaaaf19=1,xyzzyaaad19
xyzzyaaac19=xyzzyaaae19+xyzzyaaaf19
x(xyzzyaaaf19)=x(xyzzyaaaf19)-xyzzyaaah19*l(xyzzyaaac19)
enddo
enddo
40 return
end subroutine litvmu
subroutine livmul(n,x,l,y)
integer,intent(in) :: n
real(dp),intent(in) :: l(:),y(:)
real(dp),intent(inout) :: x(:)
integer xyzzyaaaa20,xyzzyaaab20,xyzzyaaac20
real(dp) xyzzyaaad20
real(dp),parameter :: xyzzyaaae20=0.d0
do xyzzyaaac20=1,n
if(y(xyzzyaaac20)/=xyzzyaaae20)goto 20
x(xyzzyaaac20)=xyzzyaaae20
enddo
goto 40
20 xyzzyaaab20=xyzzyaaac20*(xyzzyaaac20+1)/2
x(xyzzyaaac20)=y(xyzzyaaac20)/l(xyzzyaaab20)
if(xyzzyaaac20>=n)goto 40
xyzzyaaac20=xyzzyaaac20+1
do xyzzyaaaa20=xyzzyaaac20,n
xyzzyaaad20=dotprd(xyzzyaaaa20-1,l(xyzzyaaab20+1:),x)
xyzzyaaab20=xyzzyaaab20+xyzzyaaaa20
x(xyzzyaaaa20)=(y(xyzzyaaaa20)-xyzzyaaad20)/l(xyzzyaaab20)
enddo
40 return
end subroutine livmul
subroutine lmstep(d,g,ierr,ipivot,ka,p,qtr,r,step,v,w)
integer,intent(in) :: p
integer,intent(inout) :: ierr,ipivot(:),ka
real(dp),intent(in) :: d(:),g(:),qtr(:),r(:)
real(dp),intent(inout) :: v(:),w(:)
real(dp),intent(inout) :: step(:)
integer xyzzyaaaa21,xyzzyaaab21,xyzzyaaac21,xyzzyaaad21,xyzzyaaae21,xy&
&zzyaaaf21,xyzzyaaag21,xyzzyaaah21,xyzzyaaai21,xyzzyaaaj21,xyzzyaaak21&
&,xyzzyaaal21,xyzzyaaam21,xyzzyaaan21,xyzzyaaao21,xyzzyaaap21
integer,parameter :: xyzzyaaaq21=1,xyzzyaaar21=2,xyzzyaaas21=3,xyzzyaa&
&at21=19,xyzzyaaau21=4,xyzzyaaav21=6,xyzzyaaaw21=20,xyzzyaaax21=21,xyz&
&zyaaay21=7,xyzzyaaaz21=8,xyzzyaaba21=9,xyzzyaabb21=5
real(dp) xyzzyaabc21,xyzzyaabd21,xyzzyaabe21,xyzzyaabf21,xyzzyaabg21,x&
&yzzyaabh21,xyzzyaabi21,xyzzyaabj21,xyzzyaabk21,xyzzyaabl21,xyzzyaabm2&
&1,xyzzyaabn21,xyzzyaabo21,xyzzyaabp21, xyzzyaabq21,xyzzyaabr21,xyzzya&
&abs21,xyzzyaabt21,xyzzyaabu21,xyzzyaabv21,xyzzyaabw21,xyzzyaabx21,xyz&
&zyaaby21
real(dp),parameter :: xyzzyaabz21=256.d0,xyzzyaaca21=8.d0,xyzzyaacb21=&
&0.5d0, xyzzyaacc21=-1.d0,xyzzyaacd21=1.d0,xyzzyaace21=1.d-3,xyzzyaacf&
&21=3.d0,xyzzyaacg21=2.5d0,xyzzyaach21=0.d0
xyzzyaabq21=0.d0
xyzzyaabe21=0.d0
xyzzyaaai21=p+1
xyzzyaaaj21=xyzzyaaai21+1
xyzzyaaap21=xyzzyaaaj21+1
xyzzyaaaa21=xyzzyaaap21+1
xyzzyaaao21=xyzzyaaaa21
xyzzyaaan21=xyzzyaaao21+1
xyzzyaaak21=p*(p+1)/2
xyzzyaaam21=xyzzyaaak21+xyzzyaaao21
xyzzyaaal21=xyzzyaaam21+1
xyzzyaabr21=v(xyzzyaaaz21)
if(xyzzyaabr21>xyzzyaach21)xyzzyaabq21=v(xyzzyaaat21)/((xyzzyaaca21*(v&
&(xyzzyaaaw21)+xyzzyaacd21)+xyzzyaacf21)* xyzzyaabr21**2)
xyzzyaabo21=v(xyzzyaaax21)*xyzzyaabr21
xyzzyaabp21=v(xyzzyaaaw21)*xyzzyaabr21
xyzzyaabi21=xyzzyaacd21/xyzzyaabz21
xyzzyaabg21=xyzzyaabz21*xyzzyaabz21
xyzzyaabm21=xyzzyaach21
xyzzyaabl21=xyzzyaach21
xyzzyaabx21=xyzzyaach21
xyzzyaaag21=ka+12
if(ka<0)then
goto 10
elseif(ka==0)then
goto 20
else
goto 310
endif
10 ka=0
xyzzyaaag21=12
xyzzyaaaf21=p
if(ierr/=0)xyzzyaaaf21=abs(ierr)-1
v(xyzzyaaav21)=xyzzyaacb21*dotprd(xyzzyaaaf21,qtr,qtr)
20 v(xyzzyaaas21)=xyzzyaacc21
if(ierr/=0)goto 50
call litvmu(p,w(1:p),r,qtr)
do xyzzyaaab21=1,p
xyzzyaaae21=ipivot(xyzzyaaab21)
step(xyzzyaaab21)=d(xyzzyaaae21)*w(xyzzyaaab21)
enddo
xyzzyaabh21=v2norm(p,step)
v(xyzzyaaas21)=xyzzyaabh21
xyzzyaabn21=xyzzyaabh21-xyzzyaabr21
if(xyzzyaabn21<=xyzzyaabo21)goto 350
if(ka>0)goto 70
do xyzzyaaab21=1,p
xyzzyaaae21=ipivot(xyzzyaaab21)
step(xyzzyaaab21)=d(xyzzyaaae21)*(step(xyzzyaaab21)/xyzzyaabh21)
enddo
call livmul(p,step(1:p),r,step)
xyzzyaabv21=xyzzyaacd21/v2norm(p,step)
w(xyzzyaaaj21)=(xyzzyaabv21/xyzzyaabh21)*xyzzyaabv21
xyzzyaabl21=xyzzyaabn21*w(xyzzyaaaj21)
50 w(1:p)=g(1:p)/d(1:p)
v(xyzzyaaaq21)=v2norm(p,w)
xyzzyaabx21=v(xyzzyaaaq21)/xyzzyaabr21
if(xyzzyaabx21<=xyzzyaach21)goto 330
xyzzyaabe21=abs(v(xyzzyaabb21))*v(xyzzyaaba21)/xyzzyaabr21
70 ka=ka+1
call vcopy(xyzzyaaak21,w(xyzzyaaan21:xyzzyaaan21+xyzzyaaak21-1),r)
call vcopy(p,w(xyzzyaaal21:xyzzyaaal21+p-1),qtr)
if(xyzzyaabe21<=xyzzyaach21.or.xyzzyaabe21<xyzzyaabl21.or.xyzzyaabe21>&
&=xyzzyaabx21)xyzzyaabe21=xyzzyaabx21*max(xyzzyaace21,sqrt(xyzzyaabl21&
&/xyzzyaabx21))
xyzzyaabu21=sqrt(xyzzyaabe21)
w(1:p)=xyzzyaacd21
do xyzzyaaab21=1,p
xyzzyaaah21=xyzzyaaab21*(xyzzyaaab21+1)/2+xyzzyaaao21
xyzzyaaby21=w(xyzzyaaah21)
xyzzyaabk21=xyzzyaacd21
xyzzyaabj21=w(xyzzyaaab21)
xyzzyaaae21=ipivot(xyzzyaaab21)
xyzzyaabd21=xyzzyaabu21*d(xyzzyaaae21)
if(xyzzyaabd21 >= abs(xyzzyaaby21)) goto 110
90 xyzzyaabc21=xyzzyaabd21/xyzzyaaby21
xyzzyaabf21=xyzzyaabk21*xyzzyaabc21/xyzzyaabj21
xyzzyaabv21=xyzzyaabc21*xyzzyaabf21+xyzzyaacd21
if(xyzzyaabv21 > xyzzyaacg21) goto 110
w(xyzzyaaab21)=xyzzyaabj21/xyzzyaabv21
xyzzyaabk21=xyzzyaabk21/xyzzyaabv21
w(xyzzyaaah21)=xyzzyaabv21*xyzzyaaby21
xyzzyaabc21=-xyzzyaabc21
do xyzzyaaae21=xyzzyaaab21,p
xyzzyaaah21=xyzzyaaah21+xyzzyaaae21
step(xyzzyaaae21)=xyzzyaabc21*w(xyzzyaaah21)
enddo
goto 130
110 xyzzyaabf21=xyzzyaaby21/xyzzyaabd21
xyzzyaabc21=xyzzyaabj21*xyzzyaabf21/xyzzyaabk21
xyzzyaabv21=xyzzyaabc21*xyzzyaabf21+xyzzyaacd21
if(xyzzyaabv21>xyzzyaacg21)goto 90
w(xyzzyaaab21)=xyzzyaabk21/xyzzyaabv21
xyzzyaabk21=xyzzyaabj21/xyzzyaabv21
w(xyzzyaaah21)=xyzzyaabv21*xyzzyaabd21
do xyzzyaaae21=xyzzyaaab21,p
xyzzyaaah21=xyzzyaaah21+xyzzyaaae21
xyzzyaaby21=w(xyzzyaaah21)
step(xyzzyaaae21)=-xyzzyaaby21
w(xyzzyaaah21)=xyzzyaabc21*xyzzyaaby21
enddo
130 if(xyzzyaaab21==p)goto 240
xyzzyaaac21=xyzzyaaab21+1
do xyzzyaaad21=xyzzyaaac21,p
xyzzyaaah21=xyzzyaaad21*(xyzzyaaad21+1)/2+xyzzyaaao21
xyzzyaaby21=w(xyzzyaaah21)
xyzzyaabs21=step(xyzzyaaad21-1)
xyzzyaabj21=w(xyzzyaaad21)
if(xyzzyaabj21>=xyzzyaabi21)goto 150
xyzzyaabj21=xyzzyaabj21*xyzzyaabg21
xyzzyaaby21=xyzzyaaby21/xyzzyaabz21
xyzzyaaaf21=xyzzyaaah21
do xyzzyaaae21=xyzzyaaad21,p
xyzzyaaaf21=xyzzyaaaf21+xyzzyaaae21
w(xyzzyaaaf21)=w(xyzzyaaaf21)/xyzzyaabz21
enddo
150 if(abs(xyzzyaabs21)>abs(xyzzyaaby21))goto 180
if(xyzzyaabs21==xyzzyaach21)cycle
160 xyzzyaabc21=xyzzyaabs21/xyzzyaaby21
xyzzyaabf21=xyzzyaabk21*xyzzyaabc21/xyzzyaabj21
xyzzyaabv21=xyzzyaabc21*xyzzyaabf21+xyzzyaacd21
if(xyzzyaabv21>xyzzyaacg21)goto 180
w(xyzzyaaah21)=xyzzyaabv21*xyzzyaaby21
w(xyzzyaaad21)=xyzzyaabj21/xyzzyaabv21
xyzzyaabk21=xyzzyaabk21/xyzzyaabv21
do xyzzyaaae21=xyzzyaaad21,p
xyzzyaaah21=xyzzyaaah21+xyzzyaaae21
xyzzyaaby21=w(xyzzyaaah21)
xyzzyaabt21=step(xyzzyaaae21)
w(xyzzyaaah21)=xyzzyaaby21+xyzzyaabf21*xyzzyaabt21
step(xyzzyaaae21)=xyzzyaabt21-xyzzyaabc21*xyzzyaaby21
enddo
goto 200
180 xyzzyaabf21=xyzzyaaby21/xyzzyaabs21
xyzzyaabc21=xyzzyaabj21*xyzzyaabf21/xyzzyaabk21
xyzzyaabv21=xyzzyaabc21*xyzzyaabf21+xyzzyaacd21
if(xyzzyaabv21>xyzzyaacg21)goto 160
w(xyzzyaaad21)=xyzzyaabk21/xyzzyaabv21
xyzzyaabk21=xyzzyaabj21/xyzzyaabv21
w(xyzzyaaah21)=xyzzyaabv21*xyzzyaabs21
do xyzzyaaae21=xyzzyaaad21,p
xyzzyaaah21=xyzzyaaah21+xyzzyaaae21
xyzzyaaby21=w(xyzzyaaah21)
xyzzyaabt21=step(xyzzyaaae21)
w(xyzzyaaah21)=xyzzyaabc21*xyzzyaaby21+xyzzyaabt21
step(xyzzyaaae21)=xyzzyaabf21*xyzzyaabt21-xyzzyaaby21
enddo
200 if(xyzzyaabk21>=xyzzyaabi21)cycle
xyzzyaabk21=xyzzyaabk21*xyzzyaabg21
step(xyzzyaaad21:p)=step(xyzzyaaad21:p)/xyzzyaabz21
enddo
enddo
240 call litvmu(p,w(xyzzyaaal21:xyzzyaaal21+p-1),w(xyzzyaaan21:),w(xyz&
&zyaaal21:))
do xyzzyaaab21=1,p
xyzzyaaae21=ipivot(xyzzyaaab21)
xyzzyaaaf21=xyzzyaaam21+xyzzyaaab21
xyzzyaabv21=w(xyzzyaaaf21)
step(xyzzyaaae21)=-xyzzyaabv21
w(xyzzyaaaf21)=xyzzyaabv21*d(xyzzyaaae21)
enddo
xyzzyaabh21=v2norm(p,w(xyzzyaaal21:))
xyzzyaabn21=xyzzyaabh21-xyzzyaabr21
if(xyzzyaabn21<=xyzzyaabo21.and.xyzzyaabn21>=xyzzyaabp21)goto 370
if(xyzzyaabm21==xyzzyaabn21)goto 370
xyzzyaabm21=xyzzyaabn21
if(xyzzyaabn21>xyzzyaach21)goto 270
if(ka>=xyzzyaaag21)goto 370
xyzzyaabw21=xyzzyaabe21*xyzzyaabh21*xyzzyaabh21-dotprd(p,step,g)
if(xyzzyaabe21>=xyzzyaabw21*xyzzyaabq21)goto 270
v(xyzzyaabb21)=-xyzzyaabe21
goto 380
260 if(xyzzyaabn21<xyzzyaach21)xyzzyaabx21=min(xyzzyaabx21,xyzzyaabe21&
&)
goto 280
270 if(xyzzyaabn21<xyzzyaach21)xyzzyaabx21=xyzzyaabe21
280 do xyzzyaaab21=1,p
xyzzyaaae21=ipivot(xyzzyaaab21)
xyzzyaaaf21=xyzzyaaam21+xyzzyaaab21
step(xyzzyaaab21)=d(xyzzyaaae21)*(w(xyzzyaaaf21)/xyzzyaabh21)
enddo
call livmul(p,step(1:p),w(xyzzyaaan21:),step)
step(1:p)=step(1:p)/sqrt(w(1:p))
xyzzyaabv21=xyzzyaacd21/v2norm(p,step)
xyzzyaabe21=xyzzyaabe21+xyzzyaabv21*xyzzyaabn21*xyzzyaabv21/xyzzyaabr2&
&1
xyzzyaabl21=max(xyzzyaabl21,xyzzyaabe21)
goto 70
310 xyzzyaabl21=w(xyzzyaaai21)
xyzzyaabx21=w(xyzzyaaap21)
if(v(xyzzyaaas21)>xyzzyaach21.and.v(xyzzyaaas21)-xyzzyaabr21<=xyzzyaab&
&o21)goto 20
xyzzyaabe21=abs(v(xyzzyaabb21))
xyzzyaabh21=w(xyzzyaaaa21)
xyzzyaabn21=xyzzyaabh21-xyzzyaabr21
xyzzyaabv21=v(xyzzyaaaq21)/xyzzyaabr21
if(xyzzyaabr21>v(xyzzyaaba21))goto 320
xyzzyaabx21=xyzzyaabv21
if(xyzzyaabe21<=xyzzyaach21)xyzzyaabl21=xyzzyaach21
if(v(xyzzyaaas21)>xyzzyaach21)xyzzyaabl21=max(xyzzyaabl21,(v(xyzzyaaas&
&21)-xyzzyaabr21)*w(xyzzyaaaj21))
goto 260
320 if(xyzzyaabe21<=xyzzyaach21.or.xyzzyaabx21>xyzzyaabv21)xyzzyaabx21&
&=xyzzyaabv21
xyzzyaabl21=xyzzyaach21
if(v(xyzzyaaas21)>xyzzyaach21)xyzzyaabl21=max(xyzzyaabl21,(v(xyzzyaaas&
&21)-xyzzyaabr21)*w(xyzzyaaaj21))
goto 260
330 v(xyzzyaabb21)=xyzzyaach21
xyzzyaabh21=xyzzyaach21
xyzzyaabl21=xyzzyaach21
xyzzyaabx21=xyzzyaach21
v(xyzzyaaau21)=xyzzyaach21
v(xyzzyaaay21)=xyzzyaach21
step(1:p)=xyzzyaach21
goto 390
350 xyzzyaabe21=xyzzyaach21
do xyzzyaaab21=1,p
xyzzyaaae21=ipivot(xyzzyaaab21)
step(xyzzyaaae21)=-w(xyzzyaaab21)
enddo
370 v(xyzzyaabb21)=xyzzyaabe21
380 v(xyzzyaaau21)=dotprd(p,step,g)
v(xyzzyaaay21)=xyzzyaacb21*(xyzzyaabe21*xyzzyaabh21*xyzzyaabh21-v(xyzz&
&yaaau21))
390 v(xyzzyaaar21)=xyzzyaabh21
w(xyzzyaaaa21)=xyzzyaabh21
w(xyzzyaaai21)=xyzzyaabl21
w(xyzzyaaap21)=xyzzyaabx21
v(xyzzyaaba21)=xyzzyaabr21
end subroutine lmstep
subroutine lsqrt (n1, n, l, a, irc)
integer,intent(in) :: n1, n
integer,intent(out) :: irc
real(dp),intent(in) :: a(:)
real(dp),intent(inout) :: l(:)
integer xyzzyaaaa22,xyzzyaaab22,xyzzyaaac22,xyzzyaaad22,xyzzyaaae22,xy&
&zzyaaaf22,xyzzyaaag22,xyzzyaaah22,xyzzyaaai22,xyzzyaaaj22
real(dp) :: xyzzyaaak22,xyzzyaaal22,xyzzyaaam22=0.d0
xyzzyaaae22=n1*(n1-1)/2
do xyzzyaaaa22=n1,n
xyzzyaaal22=xyzzyaaam22
if(xyzzyaaaa22==1)goto 40
xyzzyaaai22=0
xyzzyaaad22=xyzzyaaaa22-1
do xyzzyaaaf22=1,xyzzyaaad22
xyzzyaaak22=xyzzyaaam22
if(xyzzyaaaf22==1)goto 20
xyzzyaaah22=xyzzyaaaf22-1
do xyzzyaaaj22=1,xyzzyaaah22
xyzzyaaac22=xyzzyaaae22+xyzzyaaaj22
xyzzyaaag22=xyzzyaaai22+xyzzyaaaj22
xyzzyaaak22=xyzzyaaak22+l(xyzzyaaac22)*l(xyzzyaaag22)
enddo
20 xyzzyaaab22=xyzzyaaae22+xyzzyaaaf22
xyzzyaaai22=xyzzyaaai22+xyzzyaaaf22
xyzzyaaak22=(a(xyzzyaaab22)-xyzzyaaak22)/l(xyzzyaaai22)
l(xyzzyaaab22)=xyzzyaaak22
xyzzyaaal22=xyzzyaaal22+xyzzyaaak22*xyzzyaaak22
enddo
40 xyzzyaaae22=xyzzyaaae22+xyzzyaaaa22
xyzzyaaak22=a(xyzzyaaae22)-xyzzyaaal22
if(xyzzyaaak22<=xyzzyaaam22)goto 60
l(xyzzyaaae22)=sqrt(xyzzyaaak22)
enddo
irc=0
goto 70
60 l(xyzzyaaae22)=xyzzyaaak22
irc=xyzzyaaaa22
70 return
end subroutine lsqrt
subroutine lsvmin (p,l,x,y,fn_val)
integer,intent(in) :: p
real(dp),intent(in) :: l(:)
real(dp),intent(inout) :: x(:), y(:)
real(dp),intent(out) :: fn_val
integer xyzzyaaaa23,xyzzyaaab23,xyzzyaaac23,xyzzyaaad23,xyzzyaaae23,xy&
&zzyaaaf23,xyzzyaaag23,xyzzyaaah23,xyzzyaaai23
real(dp) xyzzyaaaj23,xyzzyaaak23,xyzzyaaal23,xyzzyaaam23,xyzzyaaan23,x&
&yzzyaaao23,xyzzyaaap23
real(dp),parameter :: xyzzyaaaq23=0.5d0,xyzzyaaar23=1.d0,xyzzyaaas23=9&
&973.d0,xyzzyaaat23=0.d0
integer,save :: xyzzyaaau23=2
xyzzyaaab23=0
do xyzzyaaaa23=1,p
x(xyzzyaaaa23)=xyzzyaaat23
xyzzyaaab23=xyzzyaaab23+xyzzyaaaa23
if(l(xyzzyaaab23)==xyzzyaaat23)goto 100
enddo
if(mod(xyzzyaaau23,9973)==0)xyzzyaaau23=2
xyzzyaaai23=p+1
do xyzzyaaaf23=1,p
xyzzyaaac23=xyzzyaaai23-xyzzyaaaf23
xyzzyaaau23=mod(3432*xyzzyaaau23,9973)
xyzzyaaaj23=xyzzyaaaq23*(xyzzyaaar23+real(xyzzyaaau23,dp)/xyzzyaaas23)
xyzzyaaap23=(xyzzyaaaj23-x(xyzzyaaac23))
xyzzyaaao23=(-xyzzyaaaj23-x(xyzzyaaac23))
xyzzyaaam23=abs(xyzzyaaap23)
xyzzyaaal23=abs(xyzzyaaao23)
xyzzyaaag23=xyzzyaaac23-1
xyzzyaaah23=xyzzyaaac23*xyzzyaaag23/2
xyzzyaaae23=xyzzyaaah23+xyzzyaaac23
xyzzyaaap23=xyzzyaaap23/l(xyzzyaaae23)
xyzzyaaao23=xyzzyaaao23/l(xyzzyaaae23)
if(xyzzyaaag23 == 0) goto 30
do xyzzyaaaa23=1,xyzzyaaag23
xyzzyaaad23=xyzzyaaah23+xyzzyaaaa23
xyzzyaaam23=xyzzyaaam23+abs(x(xyzzyaaaa23)+l(xyzzyaaad23)*xyzzyaaap23)
xyzzyaaal23=xyzzyaaal23+abs(x(xyzzyaaaa23)+l(xyzzyaaad23)*xyzzyaaao23)
enddo
30 if(xyzzyaaal23>xyzzyaaam23)xyzzyaaap23=xyzzyaaao23
x(xyzzyaaac23)=xyzzyaaap23
if(xyzzyaaag23==0)cycle
do xyzzyaaaa23=1,xyzzyaaag23
xyzzyaaad23=xyzzyaaah23+xyzzyaaaa23
x(xyzzyaaaa23)=x(xyzzyaaaa23)+l(xyzzyaaad23)*xyzzyaaap23
enddo
enddo
xyzzyaaan23=xyzzyaaar23/v2norm(p,x)
x(1:p)=xyzzyaaan23*x(1:p)
do xyzzyaaac23=1,p
xyzzyaaak23=xyzzyaaat23
xyzzyaaag23=xyzzyaaac23-1
xyzzyaaah23=xyzzyaaac23*xyzzyaaag23/2
if(xyzzyaaag23==0)goto 80
do xyzzyaaaa23=1,xyzzyaaag23
xyzzyaaad23=xyzzyaaah23+xyzzyaaaa23
xyzzyaaak23=xyzzyaaak23+l(xyzzyaaad23)*y(xyzzyaaaa23)
enddo
80 xyzzyaaae23=xyzzyaaah23+xyzzyaaac23
y(xyzzyaaac23)=(x(xyzzyaaac23)-xyzzyaaak23)/l(xyzzyaaae23)
enddo
fn_val=xyzzyaaar23/v2norm(p,y)
goto 110
100 fn_val=xyzzyaaat23
110 return
end subroutine lsvmin
subroutine ltsqar(n,a,l)
integer,intent(in) :: n
real(dp),intent(in) :: l(:)
real(dp),intent(inout) :: a(:)
integer xyzzyaaaa24,xyzzyaaab24,xyzzyaaac24,xyzzyaaad24,xyzzyaaae24,xy&
&zzyaaaf24,xyzzyaaag24
real(dp) xyzzyaaah24,xyzzyaaai24
xyzzyaaab24=0
do xyzzyaaaa24=1,n
xyzzyaaad24=xyzzyaaab24+1
xyzzyaaab24=xyzzyaaab24+xyzzyaaaa24
xyzzyaaag24=1
if(xyzzyaaaa24==1)goto 30
xyzzyaaac24=xyzzyaaab24-1
do xyzzyaaae24=xyzzyaaad24,xyzzyaaac24
xyzzyaaai24=l(xyzzyaaae24)
do xyzzyaaaf24=xyzzyaaad24,xyzzyaaae24
a(xyzzyaaag24)=a(xyzzyaaag24)+xyzzyaaai24*l(xyzzyaaaf24)
xyzzyaaag24=xyzzyaaag24+1
enddo
enddo
30 xyzzyaaah24=l(xyzzyaaab24)
do xyzzyaaae24=xyzzyaaad24,xyzzyaaab24
a(xyzzyaaae24)=xyzzyaaah24*l(xyzzyaaae24)
enddo
enddo
end subroutine ltsqar
subroutine parchk(iv,n,nn,p,v)
integer,intent(in) :: n,nn,p
integer,intent(inout) :: iv(:)
real(dp),intent(inout) :: v(:)
integer xyzzyaaaa25,xyzzyaaab25,xyzzyaaac25,xyzzyaaad25,xyzzyaaae25,xy&
&zzyaaaf25,xyzzyaaag25
real(dp) :: xyzzyaaah25, xyzzyaaai25
integer,parameter :: xyzzyaaaj25=27
real(dp),parameter :: xyzzyaaak25=0.d0
integer,parameter :: xyzzyaaal25=16,xyzzyaaam25=29,xyzzyaaan25=37,xyzz&
&yaaao25=19,xyzzyaaap25=25, xyzzyaaaq25=39,xyzzyaaar25=86,xyzzyaaas25=&
&87,xyzzyaaat25=45,xyzzyaaau25=46,xyzzyaaav25=47,xyzzyaaaw25=20,xyzzya&
&aax25=51,xyzzyaaay25=21
real(dp),save :: xyzzyaaaz25=0.d0,xyzzyaaba25=1.d0
character(len=4),parameter :: xyzzyaabb25(2,27) = reshape( (/  'epsl',&
&'on..', 'phmn','fc..',  'phmx','fc..', 'decf','ac..',  'incf','ac..',&
& 'rdfc','mn..',  'rdfc','mx..', 'tune','r1..',  'tune','r2..', 'tune'&
&,'r3..',  'tune','r4..', 'tune','r5..',  'afct','ol..', 'rfct','ol..'&
&,  'xcto','l...', 'xfto','l...',  'lmax','0...', 'dltf','dj..',  'd0i&
&n','it..', 'dini','t...',  'jtin','it..', 'dltf','dc..',  'dfac','...&
&.', 'rlim','it..',  'cosm','in..', 'delt','a0..',  'fuzz','....' /), &
&(/ 2, 27 /) )
real(dp) :: xyzzyaabc25(27) = (/ 1.0d-3, -0.99d0, 1.0d-3, 1.0d-2, 1.2d&
&0, 1.d-2, 1.2d0, 0.d0, 0.d0, 1.d-3, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.&
&d0, 0.d0, 0.d0,    0.d0, -10.d0, 0.d0, 0.d0, 0.d0, 1.d10, 0.d0, 0.d0,&
& 1.01d0 /)
real(dp) :: xyzzyaabd25(27) = (/ 0.9d0, -1.d-3, 1.d1, 0.8d0,    1.d2, &
&0.8d0, 1.d2, 0.5d0, 0.5d0, 1.d0, 1.d0, 1.d0, 1.d0, 0.1d0, 1.d0, 1.d0,&
& 1.d0, 1.d0,   1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0,    1.d0, 1.d&
&2 /)
character (len=4), parameter :: xyzzyaabe25(3) = (/ '---c', 'hang', 'e&
&d v' /), xyzzyaabf25(3) = (/ 'nond', 'efau', 'lt v' /)
if(iv(1)==0)call dfault(iv,v)
xyzzyaaag25=iv(xyzzyaaay25)
xyzzyaaab25=iv(1)
if(xyzzyaaab25/=12)goto 30
if(nn>=n.and.n>=p.and.p>=1)goto 20
iv(1)=16
if(xyzzyaaag25>=0)write(xyzzyaaag25,10)nn,n,p
10 format(' NL2SOL: Bad nn, n, or p. nn =', i5, ', n =', i5, ', p =', &
&i5/)
goto 300
20 xyzzyaaad25=iv(21)
call dfault(iv(21:),v(33:))
iv(21)=xyzzyaaad25
iv(xyzzyaaam25)=iv(xyzzyaaal25+20)
iv(xyzzyaaat25)=n
iv(xyzzyaaau25)=nn
iv(xyzzyaaav25)=p
goto 80
30 if(n==iv(xyzzyaaat25).and.nn==iv(xyzzyaaau25).and.p==iv(xyzzyaaav25&
&))goto 50
iv(1)=17
if(xyzzyaaag25>=0)write(xyzzyaaag25,40)iv(xyzzyaaau25),iv(xyzzyaaat25)&
&,iv(xyzzyaaav25),nn,n,p
40 format(' NL2SOL: (nn,n,p) changed from (', i5, ',', i5, ',', i3,')'&
&/' to (' , i5, ',', i5, ',', i3, ').'/)
goto 300
50 if(xyzzyaaab25<=11.and.xyzzyaaab25>=1)goto 70
iv(1)=50
if(xyzzyaaag25>=0)write(xyzzyaaag25,60)xyzzyaaab25
60 format(' NL2SOL: iv(1) =', i5, ' should be between 0 and 12.')
goto 300
70 continue
80 if(xyzzyaaaz25>xyzzyaaba25)goto 90
xyzzyaaba25=rmdcon(1)
xyzzyaaah25=rmdcon(3)
xyzzyaaaz25=rmdcon(6)
xyzzyaabc25(12)=xyzzyaaah25
xyzzyaabd25(12)=xyzzyaaaz25
xyzzyaabc25(13)=xyzzyaaba25
xyzzyaabd25(13)=xyzzyaaaz25
xyzzyaabc25(14)=xyzzyaaah25
xyzzyaabc25(17)=xyzzyaaba25
xyzzyaabd25(17)=xyzzyaaaz25
xyzzyaabc25(18)=xyzzyaaah25
xyzzyaabd25(19)=xyzzyaaaz25
xyzzyaabd25(20)=xyzzyaaaz25
xyzzyaabd25(21)=xyzzyaaaz25
xyzzyaabc25(22)=xyzzyaaah25
xyzzyaabd25(24)=rmdcon(5)
xyzzyaabc25(25)=xyzzyaaah25
xyzzyaabc25(26)=xyzzyaaah25
90 xyzzyaaaf25=0
if(iv(xyzzyaaap25)>=0.and.iv(xyzzyaaap25)<=2)goto 110
xyzzyaaaf25=18
if(xyzzyaaag25>=0)write(xyzzyaaag25,100)iv(xyzzyaaap25)
100 format (' NL2SOL: iv(25) =', i4, ' should be between 0 and 2.')
110 xyzzyaaad25=xyzzyaaao25
do xyzzyaaaa25=1,xyzzyaaaj25
xyzzyaaai25=v(xyzzyaaad25)
if(xyzzyaaai25>=xyzzyaabc25(xyzzyaaaa25).and.xyzzyaaai25<=xyzzyaabd25(&
&xyzzyaaaa25))goto 130
xyzzyaaaf25=xyzzyaaad25
if(xyzzyaaag25>=0)write(xyzzyaaag25,120)xyzzyaaad25,xyzzyaaai25,xyzzya&
&abc25(xyzzyaaaa25),xyzzyaabd25(xyzzyaaaa25)
120 format (' NL2SOL: v(',i2,') =',e12.4, ' should be between',e12.4,'&
& and',e12.4)
130 xyzzyaaad25=xyzzyaaad25+1
enddo
if(xyzzyaaab25==12.and.v(xyzzyaaaq25)>xyzzyaaak25)goto 170
xyzzyaaac25=xyzzyaaar25+p
do xyzzyaaaa25=xyzzyaaas25,xyzzyaaac25
if(v(xyzzyaaaa25)>xyzzyaaak25)cycle
xyzzyaaad25=xyzzyaaaa25-xyzzyaaar25
if(xyzzyaaag25>=0)write(xyzzyaaag25,150)xyzzyaaad25,xyzzyaaaa25,v(xyzz&
&yaaaa25)
150 format(' NL2SOL: jtol(',i3, ') = v(',i3, ') =',f11.3,' should be p&
&ositive.')
xyzzyaaaf25=xyzzyaaaa25
enddo
170 if(xyzzyaaaf25==0)goto 180
iv(1)=xyzzyaaaf25
goto 300
180 if(xyzzyaaag25<0.or.iv(xyzzyaaaw25)==0)goto 300
if(xyzzyaaab25/=12.or.iv(xyzzyaaap25)==0)goto 200
xyzzyaaaf25=1
write(xyzzyaaag25,190)iv(xyzzyaaap25)
190 format(' iv(25) =',i3)
200 if(iv(xyzzyaaal25)==iv(xyzzyaaam25))goto 220
xyzzyaaaf25=1
write(xyzzyaaag25,210)iv(xyzzyaaal25)
210 format(' iv(16) =', i3)
220 xyzzyaaad25=xyzzyaaao25
xyzzyaaae25=xyzzyaaax25
do xyzzyaaaa25=1,xyzzyaaaj25
if(v(xyzzyaaad25)==v(xyzzyaaae25))goto 250
xyzzyaaaf25=1
write(xyzzyaaag25,240)xyzzyaaad25,v(xyzzyaaad25)
240 format (' v(',i2,') =',e12.4)
250 xyzzyaaad25=xyzzyaaad25+1
xyzzyaaae25=xyzzyaaae25+1
enddo
iv(xyzzyaaam25)=iv(xyzzyaaal25)
call vcopy (xyzzyaaaj25,v(xyzzyaaax25:xyzzyaaax25+xyzzyaaaj25-1),v(xyz&
&zyaaao25:))
if(xyzzyaaab25/=12)goto 300
if(v(xyzzyaaaq25)>xyzzyaaak25)goto 280
xyzzyaaac25=xyzzyaaar25+p
write(xyzzyaaag25,*)'Initial jtol array:'
write(xyzzyaaag25,270)v(xyzzyaaas25:xyzzyaaac25)
270 format(6(1x,e12.4))
280 if(v(xyzzyaaan25)>xyzzyaaak25)goto 300
xyzzyaaad25=xyzzyaaas25+p
xyzzyaaae25=xyzzyaaad25+p-1
write(xyzzyaaag25,*)'Initial d0 array:'
write(xyzzyaaag25,290)v(xyzzyaaad25:xyzzyaaae25)
290 format (6(1x,e12.4))
300 return
end subroutine parchk
subroutine qapply(n,p,j,r,ierr)
integer,intent(in) :: n, p, ierr
real(dp),intent(in) :: j(:,:)
real(dp),intent(inout) :: r(:)
integer xyzzyaaaa26,xyzzyaaab26,xyzzyaaac26
real(dp) xyzzyaaad26
xyzzyaaaa26=p
if(ierr/=0)xyzzyaaaa26=abs(ierr)-1
do xyzzyaaab26=1,xyzzyaaaa26
xyzzyaaac26=n-xyzzyaaab26+1
xyzzyaaad26=-dotprd(xyzzyaaac26,j(xyzzyaaab26:,xyzzyaaab26),r(xyzzyaaa&
&b26:))
r(xyzzyaaab26:n)=r(xyzzyaaab26:n)+xyzzyaaad26*j(xyzzyaaab26:n,xyzzyaaa&
&b26)
enddo
end subroutine qapply
subroutine qrfact(m,n,qr,alpha,ipivot,ierr,nopivk,sum)
integer,intent(in) :: m, n, nopivk
integer,intent(inout) :: ierr, ipivot(:)
real(dp),intent(inout) :: qr(:,:)
real(dp),intent(inout) :: alpha(:),sum(:)
integer xyzzyaaaa27,xyzzyaaab27,xyzzyaaac27,xyzzyaaad27,xyzzyaaae27,xy&
&zzyaaaf27,xyzzyaaag27
real(dp) xyzzyaaah27,xyzzyaaai27,xyzzyaaaj27,xyzzyaaak27,xyzzyaaal27,x&
&yzzyaaam27,xyzzyaaan27,xyzzyaaao27
real(dp),save :: xyzzyaaap27=0._dp,xyzzyaaaq27=0._dp
real(dp),parameter :: xyzzyaaar27=1.d0,xyzzyaaas27=0.01d0,xyzzyaaat27=&
&0.99d0,xyzzyaaau27=0.d0
if(xyzzyaaaq27>xyzzyaaau27)goto 10
xyzzyaaaq27=rmdcon(1)
xyzzyaaap27=rmdcon(4)
10 ierr=0
xyzzyaaan27=xyzzyaaas27*xyzzyaaap27
do xyzzyaaab27=1,n
sum(xyzzyaaab27)=v2norm(m,qr(:,xyzzyaaab27))
ipivot(xyzzyaaab27)=xyzzyaaab27
enddo
xyzzyaaaf27=min(m,n)
do xyzzyaaad27=1,xyzzyaaaf27
xyzzyaaag27=m-xyzzyaaad27+1
xyzzyaaal27=xyzzyaaau27
xyzzyaaac27=0
if(xyzzyaaad27<=nopivk)goto 50
do xyzzyaaab27=xyzzyaaad27,n
if(xyzzyaaal27>=sum(xyzzyaaab27))cycle
xyzzyaaal27=sum(xyzzyaaab27)
xyzzyaaac27=xyzzyaaab27
enddo
if(xyzzyaaac27 == 0) goto 120
if(xyzzyaaac27 == xyzzyaaad27) goto 50
xyzzyaaaa27=ipivot(xyzzyaaad27)
ipivot(xyzzyaaad27)=ipivot(xyzzyaaac27)
ipivot(xyzzyaaac27)=xyzzyaaaa27
sum(xyzzyaaac27)=sum(xyzzyaaad27)
sum(xyzzyaaad27)=xyzzyaaal27
do xyzzyaaaa27=1,m
xyzzyaaal27=qr(xyzzyaaaa27,xyzzyaaad27)
qr(xyzzyaaaa27,xyzzyaaad27)=qr(xyzzyaaaa27,xyzzyaaac27)
qr(xyzzyaaaa27,xyzzyaaac27)=xyzzyaaal27
enddo
50 xyzzyaaak27=xyzzyaaau27
do xyzzyaaaa27=xyzzyaaad27,m
if(abs(qr(xyzzyaaaa27,xyzzyaaad27))>xyzzyaaak27)xyzzyaaak27=abs(qr(xyz&
&zyaaaa27,xyzzyaaad27))
enddo
if(xyzzyaaak27<xyzzyaaaq27)goto 110
xyzzyaaah27=v2norm(xyzzyaaag27,qr(xyzzyaaad27:,xyzzyaaad27))/xyzzyaaak&
&27
xyzzyaaal27=xyzzyaaah27**2
xyzzyaaaj27=qr(xyzzyaaad27,xyzzyaaad27)
if(xyzzyaaaj27>=xyzzyaaau27)xyzzyaaah27=-xyzzyaaah27
alpha(xyzzyaaad27)=xyzzyaaah27*xyzzyaaak27
xyzzyaaai27=xyzzyaaak27*sqrt(xyzzyaaal27-(xyzzyaaaj27*xyzzyaaah27/xyzz&
&yaaak27))
qr(xyzzyaaad27,xyzzyaaad27)=xyzzyaaaj27-alpha(xyzzyaaad27)
qr(xyzzyaaad27:m,xyzzyaaad27)=qr(xyzzyaaad27:m,xyzzyaaad27)/xyzzyaaai2&
&7
xyzzyaaae27=xyzzyaaad27+1
if(xyzzyaaae27>n)cycle
do xyzzyaaab27=xyzzyaaae27,n
xyzzyaaam27=-dotprd(xyzzyaaag27,qr(xyzzyaaad27:,xyzzyaaad27),qr(xyzzya&
&aad27:,xyzzyaaab27))
call vaxpy (xyzzyaaag27, qr(xyzzyaaad27:xyzzyaaad27+xyzzyaaag27-1,xyzz&
&yaaab27), xyzzyaaam27, qr(xyzzyaaad27:,xyzzyaaad27), qr(xyzzyaaad27:,&
&xyzzyaaab27))
if(xyzzyaaae27>m)cycle
xyzzyaaao27=sum(xyzzyaaab27)
if(xyzzyaaao27<xyzzyaaaq27)cycle
xyzzyaaam27=abs(qr(xyzzyaaad27,xyzzyaaab27)/xyzzyaaao27)
if(xyzzyaaam27<xyzzyaaan27)cycle
if(xyzzyaaam27>=xyzzyaaat27)goto 80
sum(xyzzyaaab27)=xyzzyaaao27*sqrt(xyzzyaaar27-xyzzyaaam27**2)
cycle
80 sum(xyzzyaaab27)=v2norm(m-xyzzyaaad27,qr(xyzzyaaae27:,xyzzyaaab27))
enddo
enddo
goto 150
110 ierr=-xyzzyaaad27
goto 130
120 ierr=xyzzyaaad27
130 do xyzzyaaaa27=xyzzyaaad27,n
alpha(xyzzyaaaa27)=xyzzyaaau27
if(xyzzyaaaa27>xyzzyaaad27)call vscopy(xyzzyaaaa27-xyzzyaaad27,qr(xyzz&
&yaaad27:xyzzyaaad27+xyzzyaaaa27-xyzzyaaad27-1,xyzzyaaaa27),xyzzyaaau2&
&7)
enddo
150 return
end subroutine qrfact
function reldst(p,d,x,x0) result(fn_val)
integer,intent(in) :: p
real(dp),intent(in) :: d(:),x(:),x0(:)
real(dp) :: fn_val
integer xyzzyaaaa28
real(dp) :: xyzzyaaab28,xyzzyaaac28,xyzzyaaad28,xyzzyaaae28=0.d0
xyzzyaaab28=xyzzyaaae28
xyzzyaaad28=xyzzyaaae28
do xyzzyaaaa28=1,p
xyzzyaaac28=abs(d(xyzzyaaaa28)*(x(xyzzyaaaa28)-x0(xyzzyaaaa28)))
if(xyzzyaaab28<xyzzyaaac28)xyzzyaaab28=xyzzyaaac28
xyzzyaaac28=d(xyzzyaaaa28)*(abs(x(xyzzyaaaa28))+abs(x0(xyzzyaaaa28)))
if(xyzzyaaad28<xyzzyaaac28)xyzzyaaad28=xyzzyaaac28
enddo
fn_val=xyzzyaaae28
if(xyzzyaaad28>xyzzyaaae28)fn_val=xyzzyaaab28/xyzzyaaad28
end function reldst
subroutine rptmul(func,ipivot,j,p,rd,x,y,z)
integer,intent(in) :: func,p,ipivot(:)
real(dp),intent(in) :: j(:,:),rd(:)
real(dp),intent(inout) :: x(:),y(:),z(:)
integer xyzzyaaaa29,xyzzyaaab29,xyzzyaaac29,xyzzyaaad29
real(dp) xyzzyaaae29
if(func>2)goto 50
do xyzzyaaaa29=1,p
xyzzyaaac29=ipivot(xyzzyaaaa29)
z(xyzzyaaaa29)=x(xyzzyaaac29)
enddo
y(1)=z(1)*rd(1)
if(p<=1)goto 40
do xyzzyaaac29=2,p
xyzzyaaad29=xyzzyaaac29-1
xyzzyaaae29=z(xyzzyaaac29)
y(1:xyzzyaaad29)=y(1:xyzzyaaad29)+j(1:xyzzyaaad29,xyzzyaaac29)*xyzzyaa&
&ae29
y(xyzzyaaac29)=xyzzyaaae29*rd(xyzzyaaac29)
enddo
40 if(func<=1)goto 110
goto 70
50 y(1:p)=x(1:p)
70 z(1)=y(1)*rd(1)
if(p==1)goto 90
do xyzzyaaaa29=2,p
xyzzyaaab29=xyzzyaaaa29-1
z(xyzzyaaaa29)=y(xyzzyaaaa29)*rd(xyzzyaaaa29)+dotprd(xyzzyaaab29,j(1:,&
&xyzzyaaaa29),y)
enddo
90 do xyzzyaaaa29=1,p
xyzzyaaac29=ipivot(xyzzyaaaa29)
y(xyzzyaaac29)=z(xyzzyaaaa29)
enddo
110 return
end subroutine rptmul
subroutine slupdt(a,cosmin,p,size,step,u,w,wchmtd,wscale,y)
integer,intent(in) :: p
real(dp),intent(in) :: cosmin,size,step(:),wchmtd(:),y(:)
real(dp),intent(inout) :: a(:),u(:),w(:)
real(dp),intent(out) :: wscale
integer xyzzyaaaa30,xyzzyaaab30,xyzzyaaac30
real(dp) xyzzyaaad30,xyzzyaaae30,xyzzyaaaf30,xyzzyaaag30,xyzzyaaah30
real(dp),parameter :: xyzzyaaai30=0.5d0,xyzzyaaaj30=1.d0,xyzzyaaak30=0&
&.d0
xyzzyaaae30=dotprd(p,step,wchmtd)
xyzzyaaad30=cosmin*v2norm(p,step)*v2norm(p,wchmtd)
wscale=xyzzyaaaj30
if(xyzzyaaad30/=xyzzyaaak30)wscale=min(xyzzyaaaj30,abs(xyzzyaaae30/xyz&
&zyaaad30))
xyzzyaaaf30=xyzzyaaak30
if(xyzzyaaae30/=xyzzyaaak30)xyzzyaaaf30=wscale/xyzzyaaae30
w(1:p)=xyzzyaaaf30*wchmtd(1:p)
call slvmul(p,u(1:p),a,step)
xyzzyaaaf30=xyzzyaaai30*(size*dotprd(p,step,u)-dotprd(p,step,y))
u(1:p)=xyzzyaaaf30*w(1:p)+y(1:p)-size*u(1:p)
xyzzyaaac30=1
do xyzzyaaaa30=1,p
xyzzyaaag30=u(xyzzyaaaa30)
xyzzyaaah30=w(xyzzyaaaa30)
do xyzzyaaab30=1,xyzzyaaaa30
a(xyzzyaaac30)=size*a(xyzzyaaac30)+xyzzyaaag30*w(xyzzyaaab30)+xyzzyaaa&
&h30*u(xyzzyaaab30)
xyzzyaaac30=xyzzyaaac30+1
enddo
enddo
end subroutine slupdt
subroutine slvmul(p,y,s,x)
integer,intent(in) :: p
real(dp),intent(in) :: s(:), x(:)
real(dp),intent(inout) :: y(:)
integer xyzzyaaaa31,xyzzyaaab31,xyzzyaaac31,xyzzyaaad31
real(dp) xyzzyaaae31
xyzzyaaac31=1
do xyzzyaaaa31=1,p
y(xyzzyaaaa31)=dotprd(xyzzyaaaa31,s(xyzzyaaac31:),x)
xyzzyaaac31=xyzzyaaac31+xyzzyaaaa31
enddo
if(p<=1)goto 40
xyzzyaaac31=1
do xyzzyaaaa31=2,p
xyzzyaaae31=x(xyzzyaaaa31)
xyzzyaaab31=xyzzyaaaa31-1
xyzzyaaac31=xyzzyaaac31+1
do xyzzyaaad31=1,xyzzyaaab31
y(xyzzyaaad31)=y(xyzzyaaad31)+s(xyzzyaaac31)*xyzzyaaae31
xyzzyaaac31=xyzzyaaac31+1
enddo
enddo
40 return
end subroutine slvmul
function stopx() result(fn_val)
logical :: fn_val
fn_val=stop_nl2sol
stop_nl2sol=.false.
end function stopx
subroutine vaxpy(p,w,a,x,y)
integer,intent(in) :: p
real(dp),intent(in) :: a,x(:),y(:)
real(dp),intent(inout) :: w(:)
w(1:p)=a*x(1:p)+y(1:p)
end subroutine vaxpy
subroutine vcopy (p, y, x)
integer,intent(in) :: p
real(dp),intent(in) :: x(:)
real(dp),intent(inout) :: y(:)
y(1:p)=x(1:p)
end subroutine vcopy
subroutine vscopy (p, y, s)
integer,intent(in) :: p
real(dp),intent(in) :: s
real(dp),intent(inout) :: y(:)
y(1:p)=s
end subroutine vscopy
function v2norm (p, x) result(fn_val)
integer,intent(in) :: p
real(dp),intent(in) :: x(:)
real(dp) :: fn_val
integer xyzzyaaaa36,xyzzyaaab36
real(dp) xyzzyaaac36,xyzzyaaad36,xyzzyaaae36,xyzzyaaaf36
real(dp),parameter :: xyzzyaaag36=1.d0,xyzzyaaah36=0.d0
real(dp),save :: xyzzyaaai36=0.d0
if(p>0)goto 10
fn_val=xyzzyaaah36
goto 70
10 do xyzzyaaaa36=1,p
if(x(xyzzyaaaa36)/=xyzzyaaah36)goto 30
enddo
fn_val=xyzzyaaah36
goto 70
30 xyzzyaaad36=abs(x(xyzzyaaaa36))
if(xyzzyaaaa36<p)goto 40
fn_val=xyzzyaaad36
goto 70
40 xyzzyaaae36=xyzzyaaag36
if(xyzzyaaai36==xyzzyaaah36)xyzzyaaai36=rmdcon(2)
xyzzyaaab36=xyzzyaaaa36+1
do xyzzyaaaa36=xyzzyaaab36,p
xyzzyaaaf36=abs(x(xyzzyaaaa36))
if(xyzzyaaaf36>xyzzyaaad36)goto 50
xyzzyaaac36=xyzzyaaaf36/xyzzyaaad36
if(xyzzyaaac36>xyzzyaaai36)xyzzyaaae36=xyzzyaaae36+xyzzyaaac36*xyzzyaa&
&ac36
cycle
50 xyzzyaaac36=xyzzyaaad36/xyzzyaaaf36
if(xyzzyaaac36<=xyzzyaaai36)xyzzyaaac36=xyzzyaaah36
xyzzyaaae36=xyzzyaaag36+xyzzyaaae36*xyzzyaaac36*xyzzyaaac36
xyzzyaaad36=xyzzyaaaf36
enddo
fn_val=xyzzyaaad36*sqrt(xyzzyaaae36)
70 return
end function v2norm
end module toms573

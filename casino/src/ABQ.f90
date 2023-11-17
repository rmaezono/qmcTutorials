module slaarnabq
use dsp
implicit none
private
public minimum_image,min_image_brute_force
contains
subroutine minimum_image(ldvec,nvec,vec)
use slaarnabg
implicit none
integer,intent(in) :: ldvec,nvec
real(dp),intent(inout) :: vec(ldvec,*)
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2
real(dp) xyzzyaaad2,xyzzyaaae2(8),xyzzyaaaf2(3,8),xyzzyaaag2(3),xyzzya&
&aah2(3),xyzzyaaai2(9),xyzzyaaaj2,xyzzyaaak2,xyzzyaaal2(3,3),xyzzyaaam&
&2
real(dp),parameter :: xyzzyaaan2=1.d-13
logical,save :: xyzzyaaao2=.true.
select case (periodicity)
case(3)
if(xyzzyaaao2)then
xyzzyaaal2(1:3,1)=a1
xyzzyaaal2(1:3,2)=a2
xyzzyaaal2(1:3,3)=a3
call xyzzyaaab1(xyzzyaaal2)
ma1=xyzzyaaal2(1:3,1)
ma2=xyzzyaaal2(1:3,2)
ma3=xyzzyaaal2(1:3,3)
orthogonal=.false.
axis_aligned=.false.
fcc=.false.
bcc=.false.
xyzzyaaag2(1)=abs(ma1(1)*ma2(1)+ma1(2)*ma2(2)+ma1(3)*ma2(3))
xyzzyaaag2(2)=abs(ma1(1)*ma3(1)+ma1(2)*ma3(2)+ma1(3)*ma3(3))
xyzzyaaag2(3)=abs(ma2(1)*ma3(1)+ma2(2)*ma3(2)+ma2(3)*ma3(3))
if((xyzzyaaag2(1)<xyzzyaaan2).and.(xyzzyaaag2(2)<xyzzyaaan2).and.(xyzz&
&yaaag2(3)<xyzzyaaan2))then
orthogonal=.true.
if((abs(ma1(2))<xyzzyaaan2).and.(abs(ma1(3))<xyzzyaaan2).and.(abs(ma2(&
&1))<xyzzyaaan2).and.(abs(ma2(3))<xyzzyaaan2).and.(abs(ma3(1))<xyzzyaa&
&an2).and.(abs(ma3(2))<xyzzyaaan2))then
axis_aligned=.true.
endif
else
if((abs(ma1(1))<xyzzyaaan2).and.(abs(ma2(2))<xyzzyaaan2).and.(abs(ma3(&
&3))<xyzzyaaan2))then
if((abs(ma1(2)-ma1(3))<xyzzyaaan2*abs(ma1(2))).and. (abs(ma2(1)-ma2(3)&
&)<xyzzyaaan2*abs(ma2(1))).and. (abs(ma3(1)-ma3(2))<xyzzyaaan2*abs(ma3&
&(1))).and. (abs(ma1(2)-ma2(1))<xyzzyaaan2*abs(ma1(2))).and. (abs(ma1(&
&2)-ma3(1))<xyzzyaaan2*abs(ma1(2))))fcc=.true.
elseif((abs(ma1(1))<xyzzyaaan2).and.(abs(ma2(3))<xyzzyaaan2).and.(abs(&
&ma3(2))<xyzzyaaan2))then
if((abs(ma1(2)-ma1(3))<xyzzyaaan2*abs(ma1(2))).and. (abs(ma2(1)-ma2(2)&
&)<xyzzyaaan2*abs(ma2(1))).and. (abs(ma3(1)-ma3(3))<xyzzyaaan2*abs(ma3&
&(1))).and. (abs(ma1(2)-ma2(1))<xyzzyaaan2*abs(ma1(2))).and. (abs(ma1(&
&2)-ma3(1))<xyzzyaaan2*abs(ma1(2))))fcc=.true.
elseif((abs(ma1(2))<xyzzyaaan2).and.(abs(ma2(1))<xyzzyaaan2).and.(abs(&
&ma3(3))<xyzzyaaan2))then
if((abs(ma1(1)-ma1(3))<xyzzyaaan2*abs(ma1(1))).and. (abs(ma2(2)-ma2(3)&
&)<xyzzyaaan2*abs(ma2(2))).and. (abs(ma3(1)-ma3(2))<xyzzyaaan2*abs(ma3&
&(1))).and. (abs(ma1(1)-ma2(2))<xyzzyaaan2*abs(ma1(1))).and. (abs(ma1(&
&1)-ma3(1))<xyzzyaaan2*abs(ma1(1))))fcc=.true.
elseif((abs(ma1(2))<xyzzyaaan2).and.(abs(ma2(3))<xyzzyaaan2).and.(abs(&
&ma3(1))<xyzzyaaan2))then
if((abs(ma1(1)-ma1(3))<xyzzyaaan2*abs(ma1(1))).and. (abs(ma2(1)-ma2(2)&
&)<xyzzyaaan2*abs(ma2(1))).and. (abs(ma3(2)-ma3(3))<xyzzyaaan2*abs(ma3&
&(2))).and. (abs(ma1(1)-ma2(1))<xyzzyaaan2*abs(ma1(1))).and. (abs(ma1(&
&1)-ma3(2))<xyzzyaaan2*abs(ma1(1))))fcc=.true.
elseif((abs(ma1(3))<xyzzyaaan2).and.(abs(ma2(1))<xyzzyaaan2).and.(abs(&
&ma3(2))<xyzzyaaan2))then
if((abs(ma1(1)-ma1(2))<xyzzyaaan2*abs(ma1(1))).and. (abs(ma2(2)-ma2(3)&
&)<xyzzyaaan2*abs(ma2(2))).and. (abs(ma3(1)-ma3(3))<xyzzyaaan2*abs(ma3&
&(1))).and. (abs(ma1(1)-ma2(2))<xyzzyaaan2*abs(ma1(1))).and. (abs(ma1(&
&1)-ma3(1))<xyzzyaaan2*abs(ma1(1))))fcc=.true.
elseif((abs(ma1(3))<xyzzyaaan2).and.(abs(ma2(2))<xyzzyaaan2).and.(abs(&
&ma3(1))<xyzzyaaan2))then
if((abs(ma1(1)-ma1(2))<xyzzyaaan2*abs(ma1(1))).and. (abs(ma2(1)-ma2(3)&
&)<xyzzyaaan2*abs(ma2(1))).and. (abs(ma3(2)-ma3(3))<xyzzyaaan2*abs(ma3&
&(2))).and. (abs(ma1(1)-ma2(1))<xyzzyaaan2*abs(ma1(1))).and. (abs(ma1(&
&1)-ma3(2))<xyzzyaaan2*abs(ma1(1))))fcc=.true.
endif
if(fcc)then
if(dabs(ma1(1))>xyzzyaaan2)then
scubelen=2.d0*abs(ma1(1))
else
scubelen=2.d0*abs(ma1(2))
endif
goto 1
endif
xyzzyaaai2(1:3)=abs(ma1(1:3))
xyzzyaaai2(4:6)=abs(ma2(1:3))
xyzzyaaai2(7:9)=abs(ma3(1:3))
do xyzzyaaab2=2,9
if(abs(xyzzyaaai2(1)-xyzzyaaai2(xyzzyaaab2))>xyzzyaaan2*abs(xyzzyaaai2&
&(1)))goto 1
enddo
xyzzyaaah2(1)=ma1(1)*ma1(2)*ma1(3)
xyzzyaaah2(2)=ma2(1)*ma2(2)*ma2(3)
xyzzyaaah2(3)=ma3(1)*ma3(2)*ma3(3)
if((abs(xyzzyaaah2(1)-xyzzyaaah2(2))<xyzzyaaan2*abs(xyzzyaaah2(1))).an&
&d. (abs(xyzzyaaah2(1)-xyzzyaaah2(3))<xyzzyaaan2*abs(xyzzyaaah2(1))).a&
&nd. (abs(xyzzyaaah2(2)-xyzzyaaah2(3))<xyzzyaaan2*abs(xyzzyaaah2(2))))&
&then
bcc=.true.
scubelen=2.d0*abs(ma1(1))
inv_scubelen=1.d0/scubelen
endif
endif
1  xyzzyaaam2=ma1(1)*ma2(2)*ma3(3)+ma1(2)*ma2(3)*ma3(1)+ma1(3)*ma2(1)*&
&ma3(2)-ma1(3)*ma2(2)*ma3(1)-ma1(1)*ma2(3)*ma3(2)-ma1(2)*ma2(1)*ma3(3)
tb(1,1)=(ma2(2)*ma3(3)-ma2(3)*ma3(2))/xyzzyaaam2
tb(1,2)=(ma2(3)*ma3(1)-ma2(1)*ma3(3))/xyzzyaaam2
tb(1,3)=(ma2(1)*ma3(2)-ma3(1)*ma2(2))/xyzzyaaam2
tb(2,1)=(ma3(2)*ma1(3)-ma1(2)*ma3(3))/xyzzyaaam2
tb(2,2)=(ma1(1)*ma3(3)-ma3(1)*ma1(3))/xyzzyaaam2
tb(2,3)=(ma3(1)*ma1(2)-ma1(1)*ma3(2))/xyzzyaaam2
tb(3,1)=(ma1(2)*ma2(3)-ma2(2)*ma1(3))/xyzzyaaam2
tb(3,2)=(ma2(1)*ma1(3)-ma1(1)*ma2(3))/xyzzyaaam2
tb(3,3)=(ma1(1)*ma2(2)-ma1(2)*ma2(1))/xyzzyaaam2
xyzzyaaao2=.false.
endif
if(orthogonal)then
if(axis_aligned)then
do xyzzyaaaa2=1,nvec
xyzzyaaaj2=anint(tb(1,1)*vec(1,xyzzyaaaa2))
vec(1,xyzzyaaaa2)=vec(1,xyzzyaaaa2)-xyzzyaaaj2*ma1(1)
xyzzyaaaj2=anint(tb(2,2)*vec(2,xyzzyaaaa2))
vec(2,xyzzyaaaa2)=vec(2,xyzzyaaaa2)-xyzzyaaaj2*ma2(2)
xyzzyaaaj2=anint(tb(3,3)*vec(3,xyzzyaaaa2))
vec(3,xyzzyaaaa2)=vec(3,xyzzyaaaa2)-xyzzyaaaj2*ma3(3)
enddo
else
do xyzzyaaaa2=1,nvec
xyzzyaaag2(1:3)=tb(1:3,1)*vec(1,xyzzyaaaa2)+tb(1:3,2)*vec(2,xyzzyaaaa2&
&)+tb(1:3,3)*vec(3,xyzzyaaaa2)
vec(1:3,xyzzyaaaa2)=vec(1:3,xyzzyaaaa2)-anint(xyzzyaaag2(1))*ma1(1:3)-&
&anint(xyzzyaaag2(2))*ma2(1:3)-anint(xyzzyaaag2(3))*ma3(1:3)
enddo
endif
else
if(fcc)then
call xyzzyaaaa1(nvec,vec,scubelen)
elseif(bcc)then
do xyzzyaaaa2=1,nvec
xyzzyaaag2(1:3)=tb(1:3,1)*vec(1,xyzzyaaaa2)+tb(1:3,2)*vec(2,xyzzyaaaa2&
&)+tb(1:3,3)*vec(3,xyzzyaaaa2)
do xyzzyaaab2=1,3
if(xyzzyaaag2(xyzzyaaab2)>0.d0)then
xyzzyaaag2(xyzzyaaab2)=aint(xyzzyaaag2(xyzzyaaab2))
else
xyzzyaaag2(xyzzyaaab2)=aint(xyzzyaaag2(xyzzyaaab2)-1)
endif
enddo
vec(1:3,xyzzyaaaa2)=vec(1:3,xyzzyaaaa2)-xyzzyaaag2(1)*ma1(1:3)-xyzzyaa&
&ag2(2)*ma2(1:3)-xyzzyaaag2(3)*ma3(1:3)
vec(1:3,xyzzyaaaa2)=vec(1:3,xyzzyaaaa2)-scubelen*aint(2.d0*vec(1:3,xyz&
&zyaaaa2)*inv_scubelen)
if(((abs(vec(1,xyzzyaaaa2))+abs(vec(2,xyzzyaaaa2))+abs(vec(3,xyzzyaaaa&
&2)))*inv_scubelen)>=0.75d0)then
vec(1:3,xyzzyaaaa2)=vec(1:3,xyzzyaaaa2)-sign(0.5d0*scubelen,vec(1:3,xy&
&zzyaaaa2))
endif
enddo
else
do xyzzyaaaa2=1,nvec
xyzzyaaag2(1:3)=tb(1:3,1)*vec(1,xyzzyaaaa2)+tb(1:3,2)*vec(2,xyzzyaaaa2&
&)+tb(1:3,3)*vec(3,xyzzyaaaa2)
xyzzyaaag2(1)=dble(floor(xyzzyaaag2(1)))
xyzzyaaag2(2)=dble(floor(xyzzyaaag2(2)))
xyzzyaaag2(3)=dble(floor(xyzzyaaag2(3)))
vec(1:3,xyzzyaaaa2)=vec(1:3,xyzzyaaaa2)-xyzzyaaag2(1)*ma1(1:3)-xyzzyaa&
&ag2(2)*ma2(1:3)-xyzzyaaag2(3)*ma3(1:3)
xyzzyaaaf2(1:3,1)=vec(1:3,xyzzyaaaa2)
xyzzyaaaf2(1:3,2)=vec(1:3,xyzzyaaaa2)-ma1(1:3)
xyzzyaaaf2(1:3,3)=vec(1:3,xyzzyaaaa2)-ma2(1:3)
xyzzyaaaf2(1:3,4)=vec(1:3,xyzzyaaaa2)-ma3(1:3)
xyzzyaaaf2(1:3,5)=xyzzyaaaf2(1:3,2)-ma2(1:3)
xyzzyaaaf2(1:3,6)=xyzzyaaaf2(1:3,2)-ma3(1:3)
xyzzyaaaf2(1:3,7)=xyzzyaaaf2(1:3,3)-ma3(1:3)
xyzzyaaaf2(1:3,8)=xyzzyaaaf2(1:3,5)-ma3(1:3)
do xyzzyaaab2=1,8
xyzzyaaae2(xyzzyaaab2)=sum(xyzzyaaaf2(1:3,xyzzyaaab2)**2)
enddo
xyzzyaaad2=xyzzyaaae2(1)
xyzzyaaac2=1
do xyzzyaaab2=2,8
if(xyzzyaaae2(xyzzyaaab2)<xyzzyaaad2)then
xyzzyaaad2=xyzzyaaae2(xyzzyaaab2)
xyzzyaaac2=xyzzyaaab2
endif
enddo
vec(1:3,xyzzyaaaa2)=xyzzyaaaf2(1:3,xyzzyaaac2)
enddo
endif
endif
case(2)
if(xyzzyaaao2)then
xyzzyaaal2(1:3,1)=a1
xyzzyaaal2(1:3,2)=a2
xyzzyaaal2(1:3,3)=0.d0
call xyzzyaaab1(xyzzyaaal2)
ma1=xyzzyaaal2(1:3,1)
ma2=xyzzyaaal2(1:3,2)
ma3=a3
orthogonal=abs(ma1(1)*ma2(1)+ma1(2)*ma2(2))<xyzzyaaan2
xyzzyaaam2=ma1(1)*ma2(2)*ma3(3)+ma1(2)*ma2(3)*ma3(1)+ma1(3)*ma2(1)*ma3&
&(2)-ma1(3)*ma2(2)*ma3(1)-ma1(1)*ma2(3)*ma3(2)-ma1(2)*ma2(1)*ma3(3)
tb(1,1)=(ma2(2)*ma3(3)-ma2(3)*ma3(2))/xyzzyaaam2
tb(1,2)=(ma2(3)*ma3(1)-ma2(1)*ma3(3))/xyzzyaaam2
tb(1,3)=(ma2(1)*ma3(2)-ma3(1)*ma2(2))/xyzzyaaam2
tb(2,1)=(ma3(2)*ma1(3)-ma1(2)*ma3(3))/xyzzyaaam2
tb(2,2)=(ma1(1)*ma3(3)-ma3(1)*ma1(3))/xyzzyaaam2
tb(2,3)=(ma3(1)*ma1(2)-ma1(1)*ma3(2))/xyzzyaaam2
tb(3,1)=(ma1(2)*ma2(3)-ma2(2)*ma1(3))/xyzzyaaam2
tb(3,2)=(ma2(1)*ma1(3)-ma1(1)*ma2(3))/xyzzyaaam2
tb(3,3)=(ma1(1)*ma2(2)-ma1(2)*ma2(1))/xyzzyaaam2
xyzzyaaao2=.false.
endif
if(orthogonal)then
do xyzzyaaaa2=1,nvec
xyzzyaaag2(1:2)=tb(1:2,1)*vec(1,xyzzyaaaa2)+tb(1:2,2)*vec(2,xyzzyaaaa2&
&)
vec(1:2,xyzzyaaaa2)=vec(1:2,xyzzyaaaa2)-anint(xyzzyaaag2(1))*ma1(1:2)-&
&anint(xyzzyaaag2(2))*ma2(1:2)
enddo
else
do xyzzyaaaa2=1,nvec
xyzzyaaag2(1:2)=tb(1:2,1)*vec(1,xyzzyaaaa2)+tb(1:2,2)*vec(2,xyzzyaaaa2&
&)
xyzzyaaag2(1)=dble(floor(xyzzyaaag2(1)))
xyzzyaaag2(2)=dble(floor(xyzzyaaag2(2)))
vec(1:2,xyzzyaaaa2)=vec(1:2,xyzzyaaaa2)-xyzzyaaag2(1)*ma1(1:2)-xyzzyaa&
&ag2(2)*ma2(1:2)
xyzzyaaaf2(1:2,1)=vec(1:2,xyzzyaaaa2)
xyzzyaaaf2(1:2,2)=vec(1:2,xyzzyaaaa2)-ma1(1:2)
xyzzyaaaf2(1:2,3)=vec(1:2,xyzzyaaaa2)-ma2(1:2)
xyzzyaaaf2(1:2,4)=xyzzyaaaf2(1:2,2)-ma2(1:2)
do xyzzyaaab2=1,4
xyzzyaaae2(xyzzyaaab2)=sum(xyzzyaaaf2(1:2,xyzzyaaab2)**2)
enddo
xyzzyaaad2=xyzzyaaae2(1)
xyzzyaaac2=1
do xyzzyaaab2=2,4
if(xyzzyaaae2(xyzzyaaab2)<xyzzyaaad2)then
xyzzyaaad2=xyzzyaaae2(xyzzyaaab2)
xyzzyaaac2=xyzzyaaab2
endif
enddo
vec(1:2,xyzzyaaaa2)=xyzzyaaaf2(1:2,xyzzyaaac2)
enddo
endif
case(1)
xyzzyaaak2=1.d0/a1(1)
do xyzzyaaaa2=1,nvec
vec(1,xyzzyaaaa2)=vec(1,xyzzyaaaa2)-anint(vec(1,xyzzyaaaa2)*xyzzyaaak2&
&,kind=dp)*a1(1)
enddo
end select
end subroutine minimum_image
subroutine xyzzyaaaa1(nvec,vec,scubelen)
implicit none
integer,intent(in) :: nvec
real(dp),intent(in) :: scubelen
real(dp),intent(inout) :: vec(4,nvec)
integer xyzzyaaaa3
real(dp) xyzzyaaab3,xyzzyaaac3,xyzzyaaad3,xyzzyaaae3,xyzzyaaaf3
real(dp),save :: xyzzyaaag3,xyzzyaaah3,xyzzyaaai3,xyzzyaaaj3,xyzzyaaak&
&3,xyzzyaaal3
logical,save :: xyzzyaaam3=.true.
xyzzyaaab3=0.5d0
if(xyzzyaaam3)then
xyzzyaaag3=sqrt(2.d0)
xyzzyaaah3=xyzzyaaab3*xyzzyaaag3
xyzzyaaai3=1.d0/scubelen
xyzzyaaaj3=xyzzyaaag3*xyzzyaaai3
xyzzyaaak3=xyzzyaaab3*scubelen
xyzzyaaal3=scubelen/xyzzyaaag3
xyzzyaaam3=.false.
endif
do xyzzyaaaa3=1,nvec
xyzzyaaad3=(vec(1,xyzzyaaaa3)-vec(2,xyzzyaaaa3))*xyzzyaaai3
xyzzyaaae3=(vec(1,xyzzyaaaa3)+vec(2,xyzzyaaaa3))*xyzzyaaai3
xyzzyaaaf3=vec(3,xyzzyaaaa3)*xyzzyaaaj3
xyzzyaaad3=xyzzyaaad3-nint(xyzzyaaad3)
xyzzyaaae3=xyzzyaaae3-nint(xyzzyaaae3)
xyzzyaaaf3=xyzzyaaaf3-xyzzyaaag3*nint(xyzzyaaah3*xyzzyaaaf3)
xyzzyaaac3=0.5d0*int(abs(xyzzyaaad3)+abs(xyzzyaaae3)+xyzzyaaag3*abs(xy&
&zzyaaaf3))
xyzzyaaad3=xyzzyaaad3-sign(xyzzyaaac3,xyzzyaaad3)
xyzzyaaae3=xyzzyaaae3-sign(xyzzyaaac3,xyzzyaaae3)
xyzzyaaaf3=xyzzyaaaf3-sign(xyzzyaaac3,xyzzyaaaf3)*xyzzyaaag3
vec(1,xyzzyaaaa3)=(xyzzyaaad3+xyzzyaaae3)*xyzzyaaak3
vec(2,xyzzyaaaa3)=(xyzzyaaae3-xyzzyaaad3)*xyzzyaaak3
vec(3,xyzzyaaaa3)=xyzzyaaaf3*xyzzyaaal3
enddo
end subroutine xyzzyaaaa1
subroutine xyzzyaaab1(vecs)
implicit none
real(dp),intent(inout) :: vecs(3,3)
integer xyzzyaaaa4
real(dp) xyzzyaaab4(3,3)
logical xyzzyaaac4
iter: do
xyzzyaaab4=vecs
do xyzzyaaaa4=1,3
vecs(1:3,xyzzyaaaa4)=0.d0
xyzzyaaac4=xyzzyaaac1(vecs)
vecs(1:3,xyzzyaaaa4)=xyzzyaaab4(1:3,xyzzyaaaa4)
if(xyzzyaaac4)cycle iter
enddo
if(xyzzyaaac1(vecs))cycle
exit
enddo iter
end subroutine xyzzyaaab1
logical function xyzzyaaac1(vecs)
implicit none
real(dp),intent(inout) :: vecs(3,3)
integer xyzzyaaaa5
integer :: xyzzyaaab5=1
real(dp) xyzzyaaac5(3,4),xyzzyaaad5,xyzzyaaae5
real(dp),parameter :: xyzzyaaaf5=1.d-13
xyzzyaaad5=0.d0
do xyzzyaaaa5=1,3
xyzzyaaae5=vecs(1,xyzzyaaaa5)**2+vecs(2,xyzzyaaaa5)**2+vecs(3,xyzzyaaa&
&a5)**2
if((xyzzyaaae5-xyzzyaaad5)>xyzzyaaaf5*xyzzyaaad5)then
xyzzyaaad5=xyzzyaaae5
xyzzyaaab5=xyzzyaaaa5
endif
enddo
xyzzyaaac5(1:3,1)=vecs(1:3,1)+vecs(1:3,2)-vecs(1:3,3)
xyzzyaaac5(1:3,2)=vecs(1:3,1)-vecs(1:3,2)+vecs(1:3,3)
xyzzyaaac5(1:3,3)=-vecs(1:3,1)+vecs(1:3,2)+vecs(1:3,3)
xyzzyaaac5(1:3,4)=vecs(1:3,1)+vecs(1:3,2)+vecs(1:3,3)
xyzzyaaac1=.false.
do xyzzyaaaa5=1,4
xyzzyaaae5=xyzzyaaac5(1,xyzzyaaaa5)**2+xyzzyaaac5(2,xyzzyaaaa5)**2+xyz&
&zyaaac5(3,xyzzyaaaa5)**2
if((xyzzyaaae5-xyzzyaaad5)<-xyzzyaaaf5*xyzzyaaad5)then
vecs(1:3,xyzzyaaab5)=xyzzyaaac5(1:3,xyzzyaaaa5)
xyzzyaaac1=.true.
exit
endif
enddo
end function xyzzyaaac1
subroutine min_image_brute_force(periodicity,a,lat_vec,rec_vec,b,mag_r&
&ec_vec)
implicit none
integer,intent(in) :: periodicity
real(dp),intent(in) :: a(3),lat_vec(3,3),rec_vec(3,3)
real(dp),intent(in),optional :: mag_rec_vec(3)
real(dp),intent(out) :: b(3)
integer xyzzyaaaa6(3),xyzzyaaab6,xyzzyaaac6,xyzzyaaad6,xyzzyaaae6(3),x&
&yzzyaaaf6
real(dp) xyzzyaaag6(3,8),xyzzyaaah6,xyzzyaaai6(8),xyzzyaaaj6,xyzzyaaak&
&6,xyzzyaaal6,xyzzyaaam6(3),xyzzyaaan6(3),xyzzyaaao6(3)
if(periodicity==3)then
xyzzyaaaa6(1)=floor(dot_product(a,rec_vec(1:3,1)))
xyzzyaaaa6(2)=floor(dot_product(a,rec_vec(1:3,2)))
xyzzyaaaa6(3)=floor(dot_product(a,rec_vec(1:3,3)))
xyzzyaaag6(1:3,1)=a-xyzzyaaaa6(1)*lat_vec(1:3,1)-xyzzyaaaa6(2)*lat_vec&
&(1:3,2)-xyzzyaaaa6(3)*lat_vec(1:3,3)
xyzzyaaai6(1)=sum(xyzzyaaag6(1:3,1)**2)
xyzzyaaag6(1:3,2)=xyzzyaaag6(1:3,1)-lat_vec(1:3,1)
xyzzyaaai6(2)=sum(xyzzyaaag6(1:3,2)**2)
xyzzyaaag6(1:3,3)=xyzzyaaag6(1:3,1)-lat_vec(1:3,2)
xyzzyaaai6(3)=sum(xyzzyaaag6(1:3,3)**2)
xyzzyaaag6(1:3,4)=xyzzyaaag6(1:3,2)-lat_vec(1:3,2)
xyzzyaaai6(4)=sum(xyzzyaaag6(1:3,4)**2)
xyzzyaaag6(1:3,5)=xyzzyaaag6(1:3,1)-lat_vec(1:3,3)
xyzzyaaai6(5)=sum(xyzzyaaag6(1:3,5)**2)
xyzzyaaag6(1:3,6)=xyzzyaaag6(1:3,5)-lat_vec(1:3,1)
xyzzyaaai6(6)=sum(xyzzyaaag6(1:3,6)**2)
xyzzyaaag6(1:3,7)=xyzzyaaag6(1:3,5)-lat_vec(1:3,2)
xyzzyaaai6(7)=sum(xyzzyaaag6(1:3,7)**2)
xyzzyaaag6(1:3,8)=xyzzyaaag6(1:3,7)-lat_vec(1:3,1)
xyzzyaaai6(8)=sum(xyzzyaaag6(1:3,8)**2)
xyzzyaaaf6=minloc(xyzzyaaai6(1:8),1)
b(1:3)=xyzzyaaag6(1:3,xyzzyaaaf6)
if(present(mag_rec_vec))then
xyzzyaaah6=sum(b(1:3)**2)
xyzzyaaae6(1:3)=floor(sqrt(xyzzyaaah6)*mag_rec_vec(1:3))
if(any(xyzzyaaae6>0))then
do xyzzyaaab6=xyzzyaaaa6(1)-xyzzyaaae6(1),xyzzyaaaa6(1)+xyzzyaaae6(1)+&
&1
xyzzyaaam6=a-xyzzyaaab6*lat_vec(1:3,1)
do xyzzyaaac6=xyzzyaaaa6(2)-xyzzyaaae6(2),xyzzyaaaa6(2)+xyzzyaaae6(2)+&
&1
xyzzyaaan6=xyzzyaaam6-xyzzyaaac6*lat_vec(1:3,2)
do xyzzyaaad6=xyzzyaaaa6(3)-xyzzyaaae6(3),xyzzyaaaa6(3)+xyzzyaaae6(3)+&
&1
xyzzyaaao6=xyzzyaaan6-xyzzyaaad6*lat_vec(1:3,3)
xyzzyaaaj6=sum(xyzzyaaao6(1:3)**2)
if(xyzzyaaaj6<xyzzyaaah6)then
xyzzyaaah6=xyzzyaaaj6
b=xyzzyaaao6
endif
enddo
enddo
enddo
endif
endif
elseif(periodicity==2)then
xyzzyaaaa6(1)=floor(a(1)*rec_vec(1,1)+a(2)*rec_vec(2,1))
xyzzyaaaa6(2)=floor(a(1)*rec_vec(1,2)+a(2)*rec_vec(2,2))
xyzzyaaag6(1:2,1)=a(1:2)-xyzzyaaaa6(1)*lat_vec(1:2,1)-xyzzyaaaa6(2)*la&
&t_vec(1:2,2)
xyzzyaaai6(1)=xyzzyaaag6(1,1)*xyzzyaaag6(1,1)+xyzzyaaag6(2,1)*xyzzyaaa&
&g6(2,1)
xyzzyaaag6(1:2,2)=xyzzyaaag6(1:2,1)-lat_vec(1:2,1)
xyzzyaaai6(2)=xyzzyaaag6(1,2)*xyzzyaaag6(1,2)+xyzzyaaag6(2,2)*xyzzyaaa&
&g6(2,2)
xyzzyaaag6(1:2,3)=xyzzyaaag6(1:2,1)-lat_vec(1:2,2)
xyzzyaaai6(3)=xyzzyaaag6(1,3)*xyzzyaaag6(1,3)+xyzzyaaag6(2,3)*xyzzyaaa&
&g6(2,3)
xyzzyaaag6(1:2,4)=xyzzyaaag6(1:2,3)-lat_vec(1:2,1)
xyzzyaaai6(4)=xyzzyaaag6(1,4)*xyzzyaaag6(1,4)+xyzzyaaag6(2,4)*xyzzyaaa&
&g6(2,4)
xyzzyaaaf6=minloc(xyzzyaaai6(1:4),1)
b(1:2)=xyzzyaaag6(1:2,xyzzyaaaf6)
b(3)=a(3)
else
xyzzyaaab6=floor(a(1)*rec_vec(1,1))
xyzzyaaak6=a(1)-xyzzyaaab6*lat_vec(1,1)
xyzzyaaal6=xyzzyaaak6-lat_vec(1,1)
if(abs(xyzzyaaak6)<abs(xyzzyaaal6))then
b(1)=xyzzyaaak6
else
b(1)=xyzzyaaal6
endif
b(2:3)=a(2:3)
endif
end subroutine min_image_brute_force
end module slaarnabq

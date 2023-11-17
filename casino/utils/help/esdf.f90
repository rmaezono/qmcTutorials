module esdf
use dsp
use esdf_key,   only : load_keywords,kw,numkw
use format_utils,only : wout,wordwrap
use run_control,only : errstop,qmc_stop,errstop_master
implicit none
integer,public,parameter :: llength=180
integer,private,parameter :: nphys=57,ndump=2000
integer,private :: nrecords,nwarns,ndmp
logical keywords_not_in_list
character(llength),private,dimension(:),allocatable :: llist,warns,dum&
&p
character(llength),private,dimension(:,:),allocatable :: tlist
character(llength),public,dimension(:),allocatable :: block_data
type phys_unit
character(10) d,n
real(dp) u
end type phys_unit
type(phys_unit),private,dimension(nphys) :: phy
data phy(1)%d /'m'/
data phy(1)%n /'kg'/
data phy(1)%u /1.d0/
data phy(2)%d /'m'/
data phy(2)%n /'g'/
data phy(2)%u /1.d-3/
data phy(3)%d /'m'/
data phy(3)%n /'amu'/
data phy(3)%u /1.66054d-27/
data phy(4)%d /'l'/
data phy(4)%n /'m'/
data phy(4)%u /1.d0/
data phy(5)%d /'l'/
data phy(5)%n /'nm'/
data phy(5)%u /1.d-9/
data phy(6)%d /'l'/
data phy(6)%n /'ang'/
data phy(6)%u /1.d-10/
data phy(7)%d /'l'/
data phy(7)%n /'bohr'/
data phy(7)%u /0.52917715d-10/
data phy(8)%d /'t'/
data phy(8)%n /'s'/
data phy(8)%u /1.d0/
data phy(9)%d /'t'/
data phy(9)%n /'ns'/
data phy(9)%u /1.d-9/
data phy(10)%d /'t'/
data phy(10)%n /'ps'/
data phy(10)%u /1.d-12/
data phy(11)%d /'t'/
data phy(11)%n /'fs'/
data phy(11)%u /1.d-15/
data phy(12)%d /'e'/
data phy(12)%n /'j'/
data phy(12)%u /1.d0/
data phy(13)%d /'e'/
data phy(13)%n /'erg'/
data phy(13)%u /1.d-7/
data phy(14)%d /'e'/
data phy(14)%n /'ev'/
data phy(14)%u /1.60219d-19/
data phy(15)%d /'e'/
data phy(15)%n /'mev'/
data phy(15)%u /1.60219d-22/
data phy(16)%d /'e'/
data phy(16)%n /'ry'/
data phy(16)%u /2.17991d-18/
data phy(17)%d /'e'/
data phy(17)%n /'mry'/
data phy(17)%u /2.17991d-21/
data phy(18)%d /'e'/
data phy(18)%n /'hartree'/
data phy(18)%u /4.35982d-18/
data phy(19)%d /'e'/
data phy(19)%n /'kcal/mol'/
data phy(19)%u /6.94780d-21/
data phy(20)%d /'e'/
data phy(20)%n /'mhartree'/
data phy(20)%u /4.35982d-21/
data phy(21)%d /'e'/
data phy(21)%n /'kj/mol'/
data phy(21)%u /1.6606d-21/
data phy(22)%d /'e'/
data phy(22)%n /'hz'/
data phy(22)%u /6.6262d-34/
data phy(23)%d /'e'/
data phy(23)%n /'thz'/
data phy(23)%u /6.6262d-22/
data phy(24)%d /'e'/
data phy(24)%n /'cm-1'/
data phy(24)%u /1.986d-23/
data phy(25)%d /'e'/
data phy(25)%n /'cm^-1'/
data phy(25)%u /1.986d-23/
data phy(26)%d /'e'/
data phy(26)%n /'cm**-1'/
data phy(26)%u /1.986d-23/
data phy(27)%d /'f'/
data phy(27)%n /'N'/
data phy(27)%u /1.d0/
data phy(28)%d /'f'/
data phy(28)%n /'ev/ang'/
data phy(28)%u /1.60219d-9/
data phy(29)%d /'f'/
data phy(29)%n /'ry/bohr'/
data phy(29)%u /4.11943d-8/
data phy(30)%d /'l'/
data phy(30)%n /'cm'/
data phy(30)%u /1.d-2/
data phy(31)%d /'p'/
data phy(31)%n /'pa'/
data phy(31)%u /1.d0/
data phy(32)%d /'p'/
data phy(32)%n /'mpa'/
data phy(32)%u /1.d6/
data phy(33)%d /'p'/
data phy(33)%n /'gpa'/
data phy(33)%u /1.d9/
data phy(34)%d /'p'/
data phy(34)%n /'atm'/
data phy(34)%u /1.01325d5/
data phy(35)%d /'p'/
data phy(35)%n /'bar'/
data phy(35)%u /1.d5/
data phy(36)%d /'p'/
data phy(36)%n /'mbar'/
data phy(36)%u /1.d11/
data phy(37)%d /'p'/
data phy(37)%n /'ry/bohr**3'/
data phy(37)%u /1.47108d13/
data phy(38)%d /'p'/
data phy(38)%n /'ev/ang**3'/
data phy(38)%u /1.60219d11/
data phy(39)%d /'c'/
data phy(39)%n /'c'/
data phy(39)%u /1.d0/
data phy(40)%d /'c'/
data phy(40)%n /'e'/
data phy(40)%u /1.602177d-19/
data phy(41)%d /'d'/
data phy(41)%n /'C*m'/
data phy(41)%u /1.d0/
data phy(42)%d /'d'/
data phy(42)%n /'D'/
data phy(42)%u /3.33564d-30/
data phy(43)%d /'d'/
data phy(43)%n /'debye'/
data phy(43)%u /3.33564d-30/
data phy(44)%d /'d'/
data phy(44)%n /'e*bohr'/
data phy(44)%u /8.47835d-30/
data phy(45)%d /'d'/
data phy(45)%n /'e*ang'/
data phy(45)%u /1.602177d-29/
data phy(46)%d /'mom'/
data phy(46)%n /'kg*m**2'/
data phy(46)%u /1.d0/
data phy(47)%d /'mom'/
data phy(47)%n /'ry*fs**2'/
data phy(47)%u /2.1799d-48/
data phy(48)%d /'ef'/
data phy(48)%n /'v/m'/
data phy(48)%u /1.d0/
data phy(49)%d /'ef'/
data phy(49)%n /'v/nm'/
data phy(49)%u /1.d9/
data phy(50)%d /'ef'/
data phy(50)%n /'v/ang'/
data phy(50)%u /1.d10/
data phy(51)%d /'ef'/
data phy(51)%n /'v/bohr'/
data phy(51)%u /1.8897268d10/
data phy(52)%d/'ef'/
data phy(52)%n /'ry/bohr/e'/
data phy(52)%u /2.5711273d11/
data phy(53)%d/'ef'/
data phy(53)%n/'har/bohr/e'/
data phy(53)%u /5.1422546d11/
data phy(54)%d /'e'/
data phy(54)%n /'k'/
data phy(54)%u /1.38066d-23/
data phy(55)%d /'t'/
data phy(55)%n /'hr'/
data phy(55)%u /3600.d0/
data phy(56)%d /'t'/
data phy(56)%n /'min'/
data phy(56)%u /60.d0/
data phy(57)%d /'t'/
data phy(57)%n /'day'/
data phy(57)%u /86400.d0/
contains
subroutine esdf_init(filename)
implicit none
character(*),intent(in) :: filename
integer unit,xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2,xy&
&zzyaaaf2,xyzzyaaag2
integer,parameter :: xyzzyaaah2=3,xyzzyaaai2=3
logical xyzzyaaaj2
character(1) comment(xyzzyaaah2),divide(xyzzyaaai2)
character(llength) cjunk,ctemp
character(llength) sunstring
data comment /'#',';','!'/
data divide /' ','=',':'/
keywords_not_in_list=.false.
call load_keywords
do xyzzyaaab2=1,numkw
ctemp=esdf_reduce(kw(xyzzyaaab2)%label)
kw(xyzzyaaab2)%label=ctemp(1:30)
enddo
unit=0
call esdf_file(unit,filename,xyzzyaaaa2)
cjunk='Unable to open main input file "'//trim(filename)//'"'
if(xyzzyaaaa2==1)then
xyzzyaaag2=0
call errstop('ESDF',trim(cjunk))
else
xyzzyaaag2=huge(1)
endif
nrecords=0
do xyzzyaaab2=1,xyzzyaaag2
read(unit,'(a)',end=100)cjunk
do xyzzyaaac2=1,xyzzyaaah2
xyzzyaaad2=index(cjunk,comment(xyzzyaaac2))
if(xyzzyaaad2>0)cjunk(xyzzyaaad2:)=' '
enddo
do
xyzzyaaad2=index(cjunk,char(9))
if(xyzzyaaad2<1)exit
cjunk(xyzzyaaad2:xyzzyaaad2)=' '
enddo
if(len_trim(cjunk)>0)then
nrecords=nrecords+1
endif
enddo
100 rewind(unit)
allocate(llist(nrecords),block_data(nrecords),tlist(llength,nrecords),&
&warns(nrecords),dump(ndump))
nwarns=0
warns=' '
ndmp=0
dump=' '
nrecords=0
do xyzzyaaab2=1,xyzzyaaag2
read(unit,'(a)',end=101)cjunk
do xyzzyaaac2=1,xyzzyaaah2
xyzzyaaad2=index(cjunk,comment(xyzzyaaac2))
if(xyzzyaaad2>0)cjunk(xyzzyaaad2:)=' '
enddo
do
xyzzyaaad2=index(cjunk,char(9))
if(xyzzyaaad2<1)exit
cjunk(xyzzyaaad2:xyzzyaaad2)=' '
enddo
if(len_trim(cjunk)>0)then
nrecords=nrecords+1
llist(nrecords)=adjustl(cjunk)
endif
enddo
101 close(unit)
tlist=' '
do xyzzyaaab2=1,nrecords
ctemp=llist(xyzzyaaab2)
xyzzyaaae2=0
do
if(len_trim(ctemp)<=0)exit
xyzzyaaad2=minval(index(ctemp,divide),mask=index(ctemp,divide)>0)
if(xyzzyaaad2>1)then
xyzzyaaae2=xyzzyaaae2+1
tlist(xyzzyaaae2,xyzzyaaab2)=adjustl(ctemp(:xyzzyaaad2-1))
endif
ctemp=adjustl(ctemp(xyzzyaaad2+1:))
enddo
enddo
xyzzyaaaj2=.false.
do xyzzyaaab2=1,nrecords
if(esdf_reduce(tlist(1,xyzzyaaab2))=='%block')then
xyzzyaaaj2=.true.
sunstring=esdf_reduce(tlist(2,xyzzyaaab2))
if((count(sunstring==kw%label)==0))then
ctemp='Label "'//trim(esdf_reduce(tlist(2,xyzzyaaab2)))//'" not in key&
&word list'
if(count(ctemp==warns)==0)then
keywords_not_in_list=.true.
call esdf_warn(ctemp)
endif
endif
xyzzyaaaf2=0
do xyzzyaaac2=1,nrecords
if(esdf_reduce(tlist(1,xyzzyaaac2))=='%block')then
if(esdf_reduce(tlist(2,xyzzyaaab2))==esdf_reduce(tlist(2,xyzzyaaac2)))&
&xyzzyaaaf2=xyzzyaaaf2+1
else
if(esdf_reduce(tlist(2,xyzzyaaab2))==esdf_reduce(tlist(1,xyzzyaaac2)))&
&xyzzyaaaf2=xyzzyaaaf2+1
endif
enddo
ctemp='Label "'//trim(esdf_reduce(tlist(2,xyzzyaaab2)))//'" is multipl&
&y defined in the input file. '
if((xyzzyaaaf2>2).and.(count(ctemp==warns)==0))call esdf_warn(ctemp)
endif
sunstring=esdf_reduce(tlist(1,xyzzyaaab2))
if((count(sunstring==kw%label)==0).and.(.not.xyzzyaaaj2).and. sunstrin&
&g/='%endblock')then
ctemp='Label "'//trim(esdf_reduce(tlist(1,xyzzyaaab2)))//'" not in key&
&word list'
if(count(ctemp==warns)==0)then
call esdf_warn(ctemp)
keywords_not_in_list=.true.
endif
endif
if(.not.xyzzyaaaj2)then
xyzzyaaaf2=0
do xyzzyaaac2=1,nrecords
if(esdf_reduce(tlist(1,xyzzyaaac2))=='%block')then
if(esdf_reduce(tlist(1,xyzzyaaab2))==esdf_reduce(tlist(2,xyzzyaaac2)))&
&xyzzyaaaf2=xyzzyaaaf2+1
else
if(esdf_reduce(tlist(1,xyzzyaaab2))==esdf_reduce(tlist(1,xyzzyaaac2)))&
&xyzzyaaaf2=xyzzyaaaf2+1
endif
enddo
ctemp='Label "'//trim(esdf_reduce(tlist(1,xyzzyaaab2)))//'" is multipl&
&y defined in the input file. '
if((xyzzyaaaf2>1).and.(count(ctemp==warns)==0))call esdf_warn(ctemp)
endif
if(esdf_reduce(tlist(1,xyzzyaaab2))=='%endblock')xyzzyaaaj2=.false.
enddo
end subroutine esdf_init
function esdf_string(label,default)
character(*),intent(in) :: label,default
integer xyzzyaaaa3,xyzzyaaab3
character(llength) ctemp,esdf_string
call esdf_lblchk(label,'T')
esdf_string=default
do xyzzyaaaa3=1,nrecords
if(esdf_reduce(tlist(1,xyzzyaaaa3))==esdf_reduce(label))then
esdf_string=trim(adjustl(tlist(2,xyzzyaaaa3)))
do xyzzyaaab3=3,llength
if(len_trim(tlist(xyzzyaaab3,xyzzyaaaa3))==0)exit
esdf_string=trim(esdf_string)//' '//trim(adjustl(tlist(xyzzyaaab3,xyzz&
&yaaaa3)))
enddo
exit
endif
enddo
ndmp=ndmp+1
write(dump(ndmp),*,err=101)trim(esdf_reduce(label)),':',trim(esdf_stri&
&ng)
if(count(dump(ndmp)==dump(1:ndmp-1))>0)ndmp=ndmp-1
return
101 ctemp='Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_str&
&ing'
call esdf_warn(ctemp)
end function esdf_string
function esdf_integer(label,default)
integer,intent(in) :: default
character(*),intent(in) :: label
integer xyzzyaaaa4,esdf_integer
character(llength) ctemp
call esdf_lblchk(label,'I')
esdf_integer=default
do xyzzyaaaa4=1,nrecords
if(esdf_reduce(tlist(1,xyzzyaaaa4))==esdf_reduce(label))then
read(tlist(2,xyzzyaaaa4),*,err=100)esdf_integer
exit
endif
enddo
ndmp=ndmp+1
write(dump(ndmp),*,err=101)trim(esdf_reduce(label)),':',esdf_integer
if(count(dump(ndmp)==dump(1:ndmp-1))>0)ndmp=ndmp-1
return
100 ctemp='Unable to parse "'//trim(esdf_reduce(label))//'" in esdf_in&
&teger'
call esdf_die(ctemp)
101 ctemp='Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_int&
&eger'
call esdf_warn(ctemp)
end function esdf_integer
function esdf_single(label,default)
real(dp),intent(in) :: default
character(*),intent(in) :: label
integer xyzzyaaaa5
real(dp) esdf_single
character(llength) ctemp
call esdf_lblchk(label,'S')
esdf_single=default
do xyzzyaaaa5=1,nrecords
if(esdf_reduce(tlist(1,xyzzyaaaa5))==esdf_reduce(label))then
read(tlist(2,xyzzyaaaa5),*,err=100)esdf_single
exit
endif
enddo
ndmp=ndmp+1
write(dump(ndmp),*,err=101)trim(esdf_reduce(label)),':',esdf_single
if(count(dump(ndmp)==dump(1:ndmp-1))>0)ndmp=ndmp-1
return
100 ctemp='Unable to parse "'//trim(esdf_reduce(label))//'" in esdf_si&
&ngle'
call esdf_die(ctemp)
101 ctemp='Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_sin&
&gle'
call esdf_warn(ctemp)
end function esdf_single
function esdf_double(label,default)
real(dp),intent(in) :: default
character(*),intent(in) :: label
integer xyzzyaaaa6
real(dp) esdf_double
character(llength) ctemp
call esdf_lblchk(label,'D')
esdf_double=default
do xyzzyaaaa6=1,nrecords
if(esdf_reduce(tlist(1,xyzzyaaaa6))==esdf_reduce(label))then
read(tlist(2,xyzzyaaaa6),*,err=100)esdf_double
exit
endif
enddo
ndmp=ndmp+1
write(dump(ndmp),*,err=101)trim(esdf_reduce(label)),':',esdf_double
if(count(dump(ndmp)==dump(1:ndmp-1))>0) ndmp=ndmp-1
return
100 esdf_double=default
ctemp='Unable to parse "'//trim(esdf_reduce(label))//'" in esdf_double&
&'
call esdf_die(ctemp)
101 ctemp='Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_dou&
&ble'
call esdf_warn(ctemp)
end function esdf_double
function esdf_physical(label,default,dunit)
real(dp),intent(in) :: default
character(*),intent(in) :: label,dunit
integer xyzzyaaaa7
real(dp) esdf_physical
character(llength) ctemp,iunit
call esdf_lblchk(label,'P')
esdf_physical=default
do xyzzyaaaa7=1,nrecords
if(esdf_reduce(tlist(1,xyzzyaaaa7))==esdf_reduce(label))then
read(tlist(2,xyzzyaaaa7),*,err=100,end=100)esdf_physical
read(tlist(3,xyzzyaaaa7),*,err=100,end=100)iunit
esdf_physical=esdf_convfac(iunit,dunit)*esdf_physical
exit
endif
enddo
ndmp=ndmp+1
write(dump(ndmp),*,err=101)trim(esdf_reduce(label)),':',esdf_physical,&
&' ',trim(dunit)
if(count(dump(ndmp)==dump(1:ndmp-1))>0)ndmp=ndmp-1
return
100 esdf_physical=default
ctemp='Unable to parse "'//trim(esdf_reduce(label))//'" in esdf_physic&
&al'
call esdf_die(ctemp)
101 ctemp='Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_phy&
&sical'
call esdf_warn(ctemp)
end function esdf_physical
function esdf_defined(label,type)
character(*),intent(in) :: label
character(1),intent(in) :: type
integer xyzzyaaaa8
logical esdf_defined
character(llength) ctemp
call esdf_lblchk(label,type)
esdf_defined=.false.
do xyzzyaaaa8=1,nrecords
if(esdf_reduce(tlist(1,xyzzyaaaa8))==esdf_reduce(label))then
esdf_defined=.true.
exit
endif
enddo
if(esdf_defined)then
ndmp=ndmp+1
write(dump(ndmp),*,err=101)trim(esdf_reduce(label)),':'
if(count(dump(ndmp)==dump(1:ndmp-1))>0)ndmp=ndmp-1
endif
return
101 ctemp='Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_def&
&ined'
call esdf_warn(ctemp)
end function esdf_defined
function esdf_boolean(label,default)
character(*),intent(in) :: label
logical,intent(in) :: default
integer xyzzyaaaa9
logical esdf_boolean
character(llength) ctemp,positive(3),negative(3)
character(llength) sunstring
data positive /'yes','true','t'/
data negative /'no','false','f'/
call esdf_lblchk(label,'L')
esdf_boolean=default
do xyzzyaaaa9=1,nrecords
if(esdf_reduce(tlist(1,xyzzyaaaa9))==esdf_reduce(label))then
if(len_trim(tlist(2,xyzzyaaaa9))==0)then
esdf_boolean=.true.
exit
endif
sunstring=esdf_reduce(tlist(2,xyzzyaaaa9))
if(any(index(positive,sunstring)>0))then
esdf_boolean=.true.
exit
endif
if(any(index(negative,sunstring)>0))then
esdf_boolean=.false.
exit
endif
call esdf_die('Unable to parse boolean value')
endif
enddo
ndmp=ndmp+1
write(dump(ndmp),*,err=101)trim(esdf_reduce(label)),': ',esdf_boolean
if(count(dump(ndmp)==dump(1:ndmp-1))>0) ndmp=ndmp-1
return
101 ctemp='Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_boo&
&lean'
call esdf_warn(ctemp)
end function esdf_boolean
function esdf_block(label,nlines)
implicit none
integer,intent(out) :: nlines
character(*),intent(in) :: label
integer xyzzyaaaa10,xyzzyaaab10
logical esdf_block
character(llength) ctemp
call esdf_lblchk(label,'B')
ctemp='Block "'//trim(esdf_reduce(label))//'" not closed correctly. '
esdf_block=.false.
nlines=0
do xyzzyaaaa10=1,nrecords
if((esdf_reduce(tlist(1,xyzzyaaaa10))==esdf_reduce('%block')).and.(esd&
&f_reduce(tlist(2,xyzzyaaaa10))==esdf_reduce(label)))then
esdf_block=.true.
do
if(esdf_reduce(tlist(1,xyzzyaaaa10+nlines+1))==esdf_reduce('%endblock'&
&))exit
nlines=nlines+1
if(nlines+xyzzyaaaa10>nrecords)call esdf_die(ctemp)
block_data(nlines)=llist(xyzzyaaaa10+nlines)
enddo
if(esdf_reduce(tlist(2,xyzzyaaaa10+nlines+1))/=esdf_reduce(label))call&
& esdf_die(ctemp)
exit
endif
enddo
if(.not.esdf_block)return
ndmp=ndmp+1
write(dump(ndmp),*,err=101)'%block ',trim(esdf_reduce(label)),': '
if(count(dump(ndmp)==dump(1:ndmp-1))>0)then
ndmp=ndmp-1
return
endif
do xyzzyaaab10=1,nlines
ndmp=ndmp+1
if(ndmp>ndump)call errstop('ESDF_BLOCK','Too many lines in block : inc&
&rease the parameter NDUMP in the esdf.f90 module.')
dump(ndmp)=block_data(xyzzyaaab10)
enddo
ndmp=ndmp+1
write(dump(ndmp),*,err=101)'%endblock ',trim(esdf_reduce(label)),': '
return
101 ctemp='Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_blo&
&ck.'
call esdf_warn(ctemp)
end function esdf_block
function esdf_reduce(string_untrimmed)
character(*),intent(in) :: string_untrimmed
integer xyzzyaaaa11,xyzzyaaab11,xyzzyaaac11,xyzzyaaad11,xyzzyaaae11,xy&
&zzyaaaf11
integer,parameter :: xyzzyaaag11=3
character(1) punct(xyzzyaaag11)
character(llength) esdf_reduce,ctemp,string
data punct /'.','_','-'/
xyzzyaaaa11=ichar('A')
xyzzyaaab11=ichar('Z')
xyzzyaaac11=ichar('a')-xyzzyaaaa11
string=adjustl(string_untrimmed)
xyzzyaaaf11=len_trim(string)
esdf_reduce(1:xyzzyaaaf11)=string(1:xyzzyaaaf11)
esdf_reduce(xyzzyaaaf11+1:)=' '
do xyzzyaaae11=1,llength
xyzzyaaad11=ichar(esdf_reduce(xyzzyaaae11:xyzzyaaae11))
if((xyzzyaaad11>=xyzzyaaaa11).and.(xyzzyaaad11<=xyzzyaaab11))esdf_redu&
&ce(xyzzyaaae11:xyzzyaaae11)=char(xyzzyaaac11+xyzzyaaad11)
enddo
do xyzzyaaae11=1,xyzzyaaag11
do
xyzzyaaad11=index(esdf_reduce,punct(xyzzyaaae11))
if(xyzzyaaad11>0)then
ctemp=esdf_reduce
esdf_reduce(xyzzyaaad11:xyzzyaaaf11)=ctemp(xyzzyaaad11+1:xyzzyaaaf11)/&
&/' '
else
exit
endif
enddo
enddo
esdf_reduce=trim(adjustl(esdf_reduce))
end function esdf_reduce
function esdf_convfac(from,to)
character(*),intent(in) :: from,to
integer xyzzyaaaa12,xyzzyaaab12,xyzzyaaac12
real(dp) esdf_convfac
character(llength) ctemp
xyzzyaaab12=0
xyzzyaaac12=0
do xyzzyaaaa12=1,nphys
if(esdf_reduce(from)==phy(xyzzyaaaa12)%n)xyzzyaaab12=xyzzyaaaa12
if(esdf_reduce(to)==phy(xyzzyaaaa12)%n)xyzzyaaac12=xyzzyaaaa12
enddo
if(xyzzyaaab12==0)then
ctemp='Units not recognized in input file : '//trim(esdf_reduce(from))
call esdf_die(ctemp)
endif
if(xyzzyaaac12==0)then
ctemp='Units not recognized in program : '//trim(esdf_reduce(to))
call esdf_die(ctemp)
endif
if(phy(xyzzyaaab12)%d/=phy(xyzzyaaac12)%d)then
ctemp='Dimensions do not match : '//trim(esdf_reduce(from)) //' vs '//&
&trim(esdf_reduce(to))
call esdf_die(ctemp)
endif
esdf_convfac=phy(xyzzyaaab12)%u/phy(xyzzyaaac12)%u
end function esdf_convfac
function esdf_unit(ierr)
integer,intent(out) :: ierr
integer esdf_unit
logical xyzzyaaaa13
ierr=0
do esdf_unit=10,99
inquire(unit=esdf_unit,opened=xyzzyaaaa13,err=100)
if(.not.xyzzyaaaa13)return
enddo
call esdf_warn('Unable to find a free i/o unit using esdf_unit.')
ierr=1
return
100 call esdf_die('Error opening files by esdf_unit.')
end function esdf_unit
subroutine esdf_file(unit,filename,ierr)
integer,intent(out) :: unit,ierr
character(*),intent(in) :: filename
logical xyzzyaaaa14
unit=esdf_unit(ierr)
if(ierr>0)return
inquire(file=trim(filename),exist=xyzzyaaaa14,err=100)
if(.not.xyzzyaaaa14)goto 100
open(unit=unit,file=trim(filename),form='formatted',status='old',actio&
&n='read',err=100)
return
100 ierr=1
end subroutine esdf_file
subroutine esdf_newfile(unit,filename,ierr)
integer,intent(out) :: unit,ierr
character(*),intent(in) :: filename
unit=esdf_unit(ierr)
if(ierr>0)return
open(unit=unit,file=trim(filename),form='formatted',status='replace',e&
&rr=100)
return
100 ierr=1
end subroutine esdf_newfile
subroutine esdf_lblchk(string,typ)
character(*),intent(in) :: string
character(1),intent(in) :: typ
integer xyzzyaaaa16
logical xyzzyaaab16(numkw)
character(llength) ctemp
xyzzyaaab16=(esdf_reduce(string)==kw(:)%label)
xyzzyaaaa16=count(xyzzyaaab16)
if(xyzzyaaaa16==0)then
ctemp='Label "'//trim(esdf_reduce(string))//'" not recognized in keywo&
&rd list.'
call esdf_die(ctemp)
endif
if(xyzzyaaaa16>1)then
ctemp='Label "'//trim(esdf_reduce(string))//'" is multiply defined.'
call esdf_die(ctemp)
endif
where(xyzzyaaab16)
xyzzyaaab16(:)=(typ/=kw(:)%typ(1:1))
endwhere
if(any(xyzzyaaab16))then
ctemp='Label "'//trim(esdf_reduce(string))//'" has been used with wron&
&g type'
call esdf_die(ctemp)
endif
end subroutine esdf_lblchk
subroutine esdf_help(helpword_in,searchword_in)
implicit none
character(*),intent(in) :: helpword_in,searchword_in
integer xyzzyaaaa17,xyzzyaaab17,xyzzyaaac17,xyzzyaaad17
character(1) cl
character(20) ctyp,clev
character(60) title
character(78) ctemp
character(llength) helpword,searchword
character(llength),allocatable :: xyzzyaaae17(:)
helpword=esdf_reduce(helpword_in)
searchword=esdf_reduce(searchword_in)
if(trim(helpword)=='search')then
if(len_trim(searchword)<1)call esdf_die('"searchword" is empty.')
do xyzzyaaaa17=1,numkw
if((index(kw(xyzzyaaaa17)%label,trim(searchword))>0).or.(index(kw(xyzz&
&yaaaa17)%dscrpt,trim(searchword))>0))then
xyzzyaaab17=index(kw(xyzzyaaaa17)%dscrpt,'!*')-1
if(xyzzyaaab17==-1)call esdf_die('Keyword description incorrectly form&
&atted.')
title=kw(xyzzyaaaa17)%dscrpt(1:xyzzyaaab17)
xyzzyaaad17=len_trim(title)
if(xyzzyaaad17>80) call esdf_die('Keyword title too long.')
call wout(kw(xyzzyaaaa17)%label//trim(title))
endif
enddo
call wout()
call qmc_stop
endif
if(trim(helpword)=='all')then
do xyzzyaaaa17=1,numkw
if(len_trim(kw(xyzzyaaaa17)%label)>0)then
xyzzyaaab17=index(kw(xyzzyaaaa17)%dscrpt,'!*')-1
if(xyzzyaaab17==-1)call esdf_die('Keyword description incorrectly form&
&atted.')
title=kw(xyzzyaaaa17)%dscrpt(1:xyzzyaaab17)
xyzzyaaad17=len_trim(title)
if(xyzzyaaad17>80)call esdf_die('Keyword title too long.')
call wout(kw(xyzzyaaaa17)%label//trim(title))
endif
enddo
call wout()
call qmc_stop
endif
if(trim(helpword)=='list')then
do xyzzyaaaa17=1,numkw
if(len_trim(kw(xyzzyaaaa17)%label)>0)then
xyzzyaaab17=index(kw(xyzzyaaaa17)%dscrpt,'!*')-1
if(xyzzyaaab17==-1)call esdf_die('Keyword description incorrectly form&
&atted.')
title=kw(xyzzyaaaa17)%dscrpt(1:xyzzyaaab17)
xyzzyaaad17=len_trim(title)
if(xyzzyaaad17>80)call esdf_die('Keyword title too long.')
call wout(kw(xyzzyaaaa17)%label//" "//kw(xyzzyaaaa17)%typ(1:1)//" "//k&
&w(xyzzyaaaa17)%typ(3:3)//" "//trim(title))
endif
enddo
call wout()
call qmc_stop
endif
if(any((/'basic ','inter ','expert','dummy '/)==trim(helpword)))then
select case(trim(helpword))
case('basic')
cl='B'
case('inter')
cl='I'
case('expert')
cl='E'
case('dummy')
cl='D'
end select
do xyzzyaaaa17=1,numkw
if(kw(xyzzyaaaa17)%typ(3:3)==cl)then
xyzzyaaab17=index(kw(xyzzyaaaa17)%dscrpt,'!*')-1
if(xyzzyaaab17==-1)call esdf_die('Keyword description incorrectly form&
&atted.')
title=kw(xyzzyaaaa17)%dscrpt(1:xyzzyaaab17)
xyzzyaaad17=len_trim(title)
if(xyzzyaaad17>80) call esdf_die('Keyword title too long.')
call wout(kw(xyzzyaaaa17)%label//trim(title))
endif
enddo
call wout()
call qmc_stop
endif
allocate(xyzzyaaae17(numkw))
do xyzzyaaaa17=1,numkw
xyzzyaaae17(xyzzyaaaa17)=esdf_reduce(kw(xyzzyaaaa17)%label)
enddo
if(.not.any(xyzzyaaae17==trim(helpword)))call esdf_die('Keyword not re&
&cognized.')
if(count(xyzzyaaae17==trim(helpword))>1)call esdf_die('Keyword entry d&
&uplicated.')
deallocate(xyzzyaaae17)
do xyzzyaaaa17=1,numkw
if(esdf_reduce(kw(xyzzyaaaa17)%label)==trim(helpword))then
call wout('Keyword : '//trim(kw(xyzzyaaaa17)%label))
xyzzyaaab17=index(kw(xyzzyaaaa17)%dscrpt,'!*')+1
if(xyzzyaaab17==1)then
call wout('Title   : (unknown)')
else
title=trim(adjustl(kw(xyzzyaaaa17)%dscrpt(1:xyzzyaaab17-2)))
xyzzyaaac17=index(title,'*!')
if(xyzzyaaac17>0)title=trim(adjustl(title(xyzzyaaac17+2:)))
call wout('Title   : '//trim(title))
endif
select case(kw(xyzzyaaaa17)%typ(1:1))
case('I')
ctyp='Integer'
case('S')
ctyp='Single Precision'
case('D')
ctyp='Double Precision'
case('P')
ctyp='Physical'
case('T')
ctyp='String'
case('E')
ctyp='Defined'
case('B')
ctyp='Block'
case('L')
ctyp='Boolean'
end select
select case(kw(xyzzyaaaa17)%typ(3:3))
case('B')
clev='Basic'
case('I')
clev='Intermediate'
case('E')
clev='Expert'
case('D')
clev='Dummy'
end select
call wout('Type    : '//trim(ctyp))
call wout('Level   : '//trim(clev))
call wout()
call wout('DESCRIPTION')
call wout('-----------')
xyzzyaaab17=xyzzyaaab17+1
xyzzyaaad17=len_trim(kw(xyzzyaaaa17)%dscrpt)
do while(xyzzyaaab17<xyzzyaaad17)
ctemp=kw(xyzzyaaaa17)%dscrpt(xyzzyaaab17:min(xyzzyaaab17+77,xyzzyaaad1&
&7))
xyzzyaaac17=index(ctemp,'$')-1
if(xyzzyaaac17==-1)xyzzyaaac17=index(ctemp,' ',back=.true.)
write(ctemp,'(a)')adjustl(ctemp(:xyzzyaaac17))
call wout(trim(ctemp))
xyzzyaaab17=xyzzyaaab17+xyzzyaaac17
if(kw(xyzzyaaaa17)%dscrpt(xyzzyaaab17:xyzzyaaab17)=='$')xyzzyaaab17=xy&
&zzyaaab17+1
enddo
endif
enddo
call wout()
call qmc_stop
end subroutine esdf_help
subroutine esdf_die(string)
character(*),intent(in) :: string
call errstop_master('ESDF input',trim(string))
end subroutine esdf_die
subroutine esdf_warn(string)
character(*),intent(in) :: string
nwarns=nwarns+1
warns(nwarns)=string
end subroutine esdf_warn
subroutine esdf_warnout
integer xyzzyaaaa20
do xyzzyaaaa20=1,nwarns
call wout('INPUT ERROR: '//trim(warns(xyzzyaaaa20)))
call wout()
enddo
if(keywords_not_in_list)then
call wordwrap('The presence of unknown keywords in CASINO input files &
&is forbidden. To make a legitimate input, copy any example input file&
& from the examples directory and change the very few system-dependent&
& keywords in the SYSTEM section at the top.')
call errstop('ESDF','Quitting')
endif
end subroutine esdf_warnout
subroutine esdf_close
deallocate(llist,tlist,block_data,warns,dump,kw)
end subroutine esdf_close
subroutine esdf_dump(filename)
character(*),intent(in) :: filename
integer xyzzyaaaa22,xyzzyaaab22,xyzzyaaac22,xyzzyaaad22,xyzzyaaae22
character(llength) cjunk
character(llength),parameter :: xyzzyaaaf22=repeat(' ',llength)
call esdf_newfile(xyzzyaaaa22,filename,xyzzyaaab22)
if(xyzzyaaab22==1)call esdf_die('Unable to open example input file "in&
&put_example"')
xyzzyaaae22=19
do xyzzyaaac22=1,ndmp
xyzzyaaad22=index(dump(xyzzyaaac22),':')
if(xyzzyaaad22>0)then
cjunk=dump(xyzzyaaac22)(1:xyzzyaaad22-1)
dump(xyzzyaaac22)=trim(cjunk)//xyzzyaaaf22(1:xyzzyaaae22-xyzzyaaad22+1&
&)//': ' //trim(adjustl(dump(xyzzyaaac22)(xyzzyaaad22+1:)))//'#'
endif
enddo
xyzzyaaae22=maxval(index(dump(1:ndmp),'#',back=.true.))
xyzzyaaae22=36
do xyzzyaaac22=1,ndmp
xyzzyaaad22=index(dump(xyzzyaaac22),'#',back=.true.)
if(xyzzyaaad22>0)then
dump(xyzzyaaac22)=dump(xyzzyaaac22)(1:xyzzyaaad22-1)//xyzzyaaaf22(1:xy&
&zzyaaae22-xyzzyaaad22+1)//'#'
endif
enddo
do xyzzyaaac22=1,ndmp
xyzzyaaad22=index(dump(xyzzyaaac22),':')
if(xyzzyaaad22>0)then
cjunk=dump(xyzzyaaac22)(1:xyzzyaaad22-1)
do xyzzyaaad22=1,numkw
if(index(cjunk,trim(kw(xyzzyaaad22)%label))>0)exit
enddo
select case(kw(xyzzyaaad22)%typ(1:1))
case('I')
cjunk='Integer'
case('S')
cjunk='Single Precision'
case('D')
cjunk='Double Precision'
case('P')
cjunk='Physical'
case('T')
cjunk='String'
case('E')
cjunk='Defined'
case('B')
cjunk='Block'
case('L')
cjunk='Boolean'
end select
xyzzyaaae22=index(kw(xyzzyaaad22)%dscrpt,'!*')
dump(xyzzyaaac22)=trim(dump(xyzzyaaac22))//trim(kw(xyzzyaaad22)%dscrpt&
&(1:xyzzyaaae22-1))//' ('//trim(adjustl(cjunk))//')'
endif
enddo
do xyzzyaaac22=1,ndmp
write(xyzzyaaaa22,'(a)')adjustl(dump(xyzzyaaac22))
enddo
end subroutine esdf_dump
subroutine help_system(hword,sword)
implicit none
character(*),optional,intent(in) :: hword,sword
character(llength) helpword,searchword
searchword=''
call load_keywords
if(.not.(present(hword)))return
helpword=hword 
searchword=''
if(present(sword))searchword=sword
call wout('CASINO HELP SYSTEM')
call wout('==================')
call wout()
call esdf_help(trim(helpword),trim(searchword))
call wout()
end subroutine help_system
end module esdf

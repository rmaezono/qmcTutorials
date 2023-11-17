module casl
implicit none
private
public read_casl,write_casl,push_casl_context,pop_casl_context,query_c&
&asl_item,get_casl_item,check_unread_casl,delete_casl_item,set_casl_it&
&em,set_casl_block,first_unread_child,unique_casl_string
public casl_keysize,casl_valsize,casl_fullkeysize
type metachar
private
character,pointer :: chars(:)=>null()
end type metachar
type casl_item
private
type(casl_item),pointer :: next=>null(),prev=>null(),parent=>null(),fi&
&rst_child=>null(),last_child=>null()
logical :: is_block=.false.,is_implicit=.false.,is_inline=.false.,been&
&_read=.false.
integer :: ilevel=0,indent=-1,nimplicit=0
type(metachar) :: name,value,unique_name,full_unique_name
type(casl_item),pointer :: avl_left=>null(),avl_right=>null(),avl_pare&
&nt=>null()
integer :: avl_depth=0,avl_imbalance=0
end type casl_item
type(casl_item),pointer :: xyzzyaaaa1=>null()
type(casl_item),pointer :: xyzzyaaab1=>null()
type casl_list
type(casl_item),pointer :: item=>null()
type(casl_list),pointer :: next=>null(),prev=>null()
end type casl_list
type(casl_list),pointer :: xyzzyaaac1=>null()
integer,parameter :: casl_keysize=32
integer,parameter :: casl_fullkeysize=256
integer,parameter :: casl_valsize=2048
integer,parameter :: xyzzyaaad1=2
integer,parameter :: xyzzyaaae1=5
interface get_casl_item
module procedure xyzzyaaas1,xyzzyaaau1,xyzzyaaaw1,xyzzyaaay1,xyzzyaaba&
&1,xyzzyaabc1,xyzzyaabd1
end interface
interface set_casl_item
module procedure xyzzyaabj1,xyzzyaabk1,xyzzyaabl1,xyzzyaabm1,xyzzyaabn&
&1,xyzzyaabo1,xyzzyaabp1
end interface
contains
subroutine xyzzyaaaf1
implicit none
if(.not.associated(xyzzyaaaa1))then
allocate(xyzzyaaaa1)
xyzzyaaaa1%is_block=.true.
endif
if(.not.associated(xyzzyaaac1))then
allocate(xyzzyaaac1)
xyzzyaaac1%item=>xyzzyaaaa1
endif
if(.not.associated(xyzzyaaab1))xyzzyaaab1=>xyzzyaaaa1
end subroutine xyzzyaaaf1
subroutine read_casl(filename,errmsg)
implicit none
character(*),intent(in) :: filename
character(512),intent(inout) :: errmsg
integer xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5,indent,iline,xyzzyaaad5,xyzzy&
&aaae5
logical xyzzyaaaf5,xyzzyaaag5,xyzzyaaah5
character check_eof
type(casl_item),pointer :: xyzzyaaai5,xyzzyaaaj5,item
type(metachar) string,label
type virtual_line
integer :: indent=-1,iline=-1,cont_indent=-1
type(metachar) line
type(virtual_line),pointer :: next=>null()
end type virtual_line
type(virtual_line),pointer :: xyzzyaaak5,xyzzyaaal5,xyzzyaaam5
type(casl_list),pointer :: xyzzyaaan5
errmsg=''
call xyzzyaaaf1
inquire(file=trim(filename),exist=xyzzyaaaf5)
if(.not.xyzzyaaaf5)return
xyzzyaaad5=9
xyzzyaaag5=.true.
do while(xyzzyaaag5)
xyzzyaaad5=xyzzyaaad5+1
inquire(xyzzyaaad5,opened=xyzzyaaag5)
if(xyzzyaaad5>=99)then
errmsg='Could not find free i/o unit.'
return
endif
enddo
open(unit=xyzzyaaad5,file=trim(filename),status='old',iostat=xyzzyaaaa&
&5)
if(xyzzyaaaa5/=0)then
errmsg='Problem opening '//trim(filename)//'.'
return
endif
nullify(xyzzyaaak5,xyzzyaaam5,xyzzyaaal5)
iline=0
do
iline=iline+1
read(xyzzyaaad5,*,iostat=xyzzyaaaa5)check_eof
if(xyzzyaaaa5/=0)exit
backspace(xyzzyaaad5)
call xyzzyaacd1(xyzzyaaad5,string)
do xyzzyaaab5=1,xyzzyaaca1(string)
if(string%chars(xyzzyaaab5)==char(9))string%chars(xyzzyaaab5)=' '
enddo
xyzzyaaab5=xyzzyaacl1(string,'#')
if(xyzzyaaab5>0)call xyzzyaacj1(string,1,xyzzyaaab5-1)
if(xyzzyaacb1(string)==0)cycle
do indent=0,xyzzyaaca1(string)-1
if(string%chars(indent+1)/=' ')exit
enddo
call xyzzyaack1(string)
if(associated(xyzzyaaak5))then
if(indent>xyzzyaaam5%indent)then
xyzzyaaab5=xyzzyaaap5(xyzzyaaci1(xyzzyaaam5%line,xyzzyaaca1(xyzzyaaam5&
&%line)),0,look_for=':')
if(xyzzyaaab5<-1)then
errmsg='Parse pass 1 problem: syntax error at line '//trim(i2s(iline))&
&//' of file '//trim(filename)//'.'
call xyzzyaabw1(string)
call xyzzyaaas5
close(xyzzyaaad5)
return
elseif(xyzzyaaab5/=xyzzyaacb1(xyzzyaaam5%line))then
if(xyzzyaaam5%cont_indent==-1.or.xyzzyaaam5%cont_indent==indent)then
xyzzyaaam5%cont_indent=indent
call xyzzyaace1(xyzzyaaam5%line,' ')
call xyzzyaacf1(xyzzyaaam5%line,string)
call xyzzyaabw1(string)
cycle
else
errmsg='Parse pass 1 problem: continuation line at line '//trim(i2s(il&
&ine))//' of file '//trim(filename)//' has different indentation from &
&previous continuation line(s) of the same line.'
call xyzzyaabw1(string)
call xyzzyaaas5
close(xyzzyaaad5)
return
endif
endif
endif
endif
allocate(xyzzyaaal5)
call xyzzyaabx1(string,xyzzyaaal5%line)
call xyzzyaabw1(string)
xyzzyaaal5%iline=iline
xyzzyaaal5%indent=indent
if(.not.associated(xyzzyaaak5))then
xyzzyaaak5=>xyzzyaaal5
else
xyzzyaaam5%next=>xyzzyaaal5
endif
xyzzyaaal5%next=>xyzzyaaak5
xyzzyaaam5=>xyzzyaaal5
nullify(xyzzyaaal5)
enddo
nullify(xyzzyaaam5)
close(xyzzyaaad5)
if(.not.associated(xyzzyaaak5))return
allocate(xyzzyaaai5)
xyzzyaaai5%is_block=.true.
call xyzzyaaby1(trim(filename),xyzzyaaai5%name)
call xyzzyaaby1(trim(filename),xyzzyaaai5%unique_name)
call xyzzyaaby1(trim(filename),xyzzyaaai5%full_unique_name)
xyzzyaaai5%parent=>xyzzyaaaa1
if(.not.associated(xyzzyaaaa1%first_child))then
xyzzyaaaa1%first_child=>xyzzyaaai5
xyzzyaaaa1%last_child=>xyzzyaaai5
else
xyzzyaaaa1%last_child%next=>xyzzyaaai5
xyzzyaaai5%prev=>xyzzyaaaa1%last_child
xyzzyaaaa1%last_child=>xyzzyaaai5
endif
call xyzzyaabq1(xyzzyaaai5,xyzzyaaah5)
if(xyzzyaaah5)then
errmsg='File '//trim(filename)//' read twice.'
call xyzzyaaas5
return
endif
nullify(xyzzyaaan5)
call push_casl_context(':'//trim(filename))
xyzzyaaae5=0
xyzzyaaal5=>xyzzyaaak5
xyzzyaaaj5=>xyzzyaaai5
do
allocate(item)
item%indent=xyzzyaaal5%indent
if(item%indent>xyzzyaaaj5%indent)then
xyzzyaaaj5%first_child=>item
xyzzyaaaj5%last_child=>item
item%parent=>xyzzyaaaj5
item%ilevel=xyzzyaaaj5%ilevel+1
else
do while(xyzzyaaaj5%indent>item%indent)
if(.not.associated(xyzzyaaaj5%parent))exit
xyzzyaaaj5=>xyzzyaaaj5%parent
enddo
if(xyzzyaaaj5%indent/=item%indent)then
errmsg='Bad indentation at line '//trim(i2s(xyzzyaaal5%iline))//' of f&
&ile '//trim(filename)//': no matching sibling.'
deallocate(item)
call xyzzyaaas5
return
endif
xyzzyaaaj5%next=>item
item%prev=>xyzzyaaaj5
item%parent=>xyzzyaaaj5%parent
item%parent%last_child=>item
item%ilevel=xyzzyaaaj5%ilevel
endif
xyzzyaaaj5=>item
xyzzyaaac5=0
do
xyzzyaaab5=xyzzyaaap5(xyzzyaaci1(xyzzyaaal5%line,xyzzyaaca1(xyzzyaaal5&
&%line)),xyzzyaaac5,look_for=':')
select case(xyzzyaaab5)
case(-3)
errmsg='Internal error.'
call xyzzyaaas5
return
case(-2)
errmsg='Mismatched syntactical delimiters ] ) } found when looking for&
& indicator character in line '//trim(i2s(xyzzyaaal5%iline))//' of fil&
&e '//trim(filename)//'.'
call xyzzyaaas5
return
case(-1)
errmsg='Mismatched syntactical delimiters [ ( { " found when looking f&
&or indicator character in line '//trim(i2s(xyzzyaaal5%iline))//' of f&
&ile '//trim(filename)//'.'
call xyzzyaaas5
return
case(0)
exit
end select
if(xyzzyaaab5==xyzzyaacb1(xyzzyaaal5%line))exit
if(xyzzyaaal5%line%chars(xyzzyaaab5+1)==' ')exit
xyzzyaaac5=xyzzyaaab5
enddo
if(xyzzyaaab5>0)then
call xyzzyaabx1(xyzzyaaal5%line,xyzzyaaaj5%name)
call xyzzyaacj1(xyzzyaaaj5%name,1,xyzzyaaab5-1)
call xyzzyaack1(xyzzyaaaj5%name)
if(xyzzyaacl1(xyzzyaaaj5%name,':')>0)then
errmsg='Missing space after ":" or misnamed item at line '//trim(i2s(x&
&yzzyaaal5%iline))//' of file '//trim(filename)//': an item name canno&
&t contain ":".'
call xyzzyaaas5
return
endif
if(xyzzyaaaj5%name%chars(1)=='%')then
if(xyzzyaaca1(xyzzyaaaj5%name)>1)then
if(xyzzyaaaj5%name%chars(2)=='!')then
call xyzzyaacj1(xyzzyaaaj5%name,3)
call xyzzyaack1(xyzzyaaaj5%name)
call xyzzyaacg1(xyzzyaaaj5%name,'%delete'//trim(i2s(xyzzyaaae5))//'-')
if(.not.associated(xyzzyaaan5))then
allocate(xyzzyaaan5)
else
allocate(xyzzyaaan5%next)
xyzzyaaan5%next%prev=>xyzzyaaan5
xyzzyaaan5=>xyzzyaaan5%next
endif
xyzzyaaan5%item=>xyzzyaaaj5
xyzzyaaae5=xyzzyaaae5+1
else
call xyzzyaaar1(xyzzyaaaj5,label)
errmsg='Illegal item name "'//xyzzyaaci1(label,xyzzyaaca1(label))//'":&
& character "%" at the beginning of the name is reserved for directive&
&s, and this name does not match any known directive.'
call xyzzyaabw1(label)
call xyzzyaaas5
return
endif
else
call xyzzyaaar1(xyzzyaaaj5,label)
errmsg='Illegal item name "'//xyzzyaaci1(label,xyzzyaaca1(label))//'":&
& character "%" at the beginning of the name is reserved for directive&
&s, and this name does not match any known directive.'
call xyzzyaabw1(label)
call xyzzyaaas5
return
endif
endif
call xyzzyaaby1(trim(unique_casl_string(xyzzyaaci1(xyzzyaaaj5%name,xyz&
&zyaaca1(xyzzyaaaj5%name)))),xyzzyaaaj5%unique_name)
call xyzzyaabx1(xyzzyaaaj5%unique_name,xyzzyaaaj5%full_unique_name)
call xyzzyaacg1(xyzzyaaaj5%full_unique_name,':')
call xyzzyaach1(xyzzyaaaj5%full_unique_name,xyzzyaaaj5%parent%full_uni&
&que_name)
if(xyzzyaaab5==xyzzyaacb1(xyzzyaaal5%line))then
if(associated(xyzzyaaal5%next,xyzzyaaak5))then
call xyzzyaabw1(xyzzyaaaj5%value)
else
if(xyzzyaaal5%next%indent>xyzzyaaal5%indent)then
xyzzyaaaj5%is_block=.true.
else
call xyzzyaabw1(xyzzyaaaj5%value)
endif
endif
else
call xyzzyaabx1(xyzzyaaal5%line,xyzzyaaaj5%value)
call xyzzyaacj1(xyzzyaaaj5%value,xyzzyaaab5+1)
call xyzzyaack1(xyzzyaaaj5%value)
endif
else
xyzzyaaaj5%is_implicit=.true.
xyzzyaaaj5%parent%nimplicit=xyzzyaaaj5%parent%nimplicit+1
call xyzzyaaby1('%u'//trim(i2s(xyzzyaaaj5%parent%nimplicit)),xyzzyaaaj&
&5%name)
call xyzzyaaby1('%u'//trim(i2s(xyzzyaaaj5%parent%nimplicit)),xyzzyaaaj&
&5%unique_name)
call xyzzyaabx1(xyzzyaaaj5%unique_name,xyzzyaaaj5%full_unique_name)
call xyzzyaacg1(xyzzyaaaj5%full_unique_name,':')
call xyzzyaach1(xyzzyaaaj5%full_unique_name,xyzzyaaaj5%parent%full_uni&
&que_name)
call xyzzyaabx1(xyzzyaaal5%line,xyzzyaaaj5%value)
call xyzzyaack1(xyzzyaaal5%line)
endif
call xyzzyaabq1(xyzzyaaaj5,xyzzyaaah5)
if(xyzzyaaah5)then
call xyzzyaaar1(xyzzyaaaj5,label)
errmsg='Item '//xyzzyaaci1(label,xyzzyaaca1(label))//' found twice (se&
&cond occurrence at line '//trim(i2s(xyzzyaaal5%iline))//' of file '//&
&trim(filename)//').  This may be due to bad indentation in file '//tr&
&im(filename)//'.'
call xyzzyaaas5
return
endif
if(.not.(xyzzyaaaj5%is_block.or.xyzzyaaaj5%is_implicit))then
if(xyzzyaaca1(xyzzyaaaj5%value)>0)then
if(xyzzyaaaj5%value%chars(1)=='[')then
xyzzyaaaj5%is_block=.true.
xyzzyaaaj5%is_inline=.true.
call xyzzyaaao5(xyzzyaaaj5,errmsg)
if(len_trim(errmsg)>0)then
call xyzzyaaas5
return
endif
endif
endif
endif
if(associated(xyzzyaaal5%next,xyzzyaaak5))exit
xyzzyaaal5=>xyzzyaaal5%next
enddo
call xyzzyaaar5
call xyzzyaaaq5
call pop_casl_context()
contains
recursive subroutine xyzzyaaao5(item,xyzzyaaaa6)
implicit none
type(casl_item),pointer :: item
character(512),intent(inout) :: xyzzyaaaa6
integer xyzzyaaab6,xyzzyaaac6
logical xyzzyaaad6
type(casl_item),pointer :: xyzzyaaae6
type(metachar) label
xyzzyaaaa6=''
xyzzyaaab6=xyzzyaaap5(xyzzyaaci1(item%value,xyzzyaaca1(item%value)),1)
select case(xyzzyaaab6)
case(-3)
xyzzyaaaa6='Internal error.'
return
case(-2)
call xyzzyaaar1(item,label)
xyzzyaaaa6='Mismatched syntactical delimiters ] ) } found when looking&
& for end of inline block "'//xyzzyaaci1(label,xyzzyaaca1(label))//'".&
&'
call xyzzyaabw1(label)
return
case(-1)
call xyzzyaaar1(item,label)
xyzzyaaaa6='Mismatched syntactical delimiters [ ( { " found when looki&
&ng for end of inline block "'//xyzzyaaci1(label,xyzzyaaca1(label))//'&
&".'
call xyzzyaabw1(label)
return
case(0)
call xyzzyaaar1(item,label)
xyzzyaaaa6='Cannot find closing square bracket for inline block "'//xy&
&zzyaaci1(label,xyzzyaaca1(label))//'".'
call xyzzyaabw1(label)
return
end select
if(xyzzyaaab6<xyzzyaaca1(item%value))then
call xyzzyaaar1(item,label)
xyzzyaaaa6='Trailing characters after inline block "'//xyzzyaaci1(labe&
&l,xyzzyaaca1(label))//'".'
call xyzzyaabw1(label)
return
endif
call xyzzyaacj1(item%value,2,xyzzyaaab6-1)
if(xyzzyaacb1(item%value)==0)return
nullify(xyzzyaaae6)
do
if(.not.associated(xyzzyaaae6))then
allocate(xyzzyaaae6)
item%first_child=>xyzzyaaae6
item%last_child=>xyzzyaaae6
else
allocate(xyzzyaaae6%next)
xyzzyaaae6%next%prev=>xyzzyaaae6
xyzzyaaae6=>xyzzyaaae6%next
item%last_child=>xyzzyaaae6
endif
xyzzyaaae6%parent=>item
xyzzyaaae6%ilevel=item%ilevel+1
xyzzyaaae6%indent=item%indent+1
xyzzyaaac6=0
do
xyzzyaaab6=xyzzyaaap5(xyzzyaaci1(item%value,xyzzyaaca1(item%value)),xy&
&zzyaaac6,look_for=',')
select case(xyzzyaaab6)
case(-3)
xyzzyaaaa6='Internal error.'
return
case(-2)
call xyzzyaaar1(item,label)
xyzzyaaaa6='Mismatched syntactical delimiters ] ) } found when looking&
& for separator in inline block "'//xyzzyaaci1(label,xyzzyaaca1(label)&
&)//'".'
call xyzzyaabw1(label)
return
case(-1)
call xyzzyaaar1(item,label)
xyzzyaaaa6='Mismatched syntactical delimiters [ ( { " found when looki&
&ng for separator in inline block "'//xyzzyaaci1(label,xyzzyaaca1(labe&
&l))//'".'
call xyzzyaabw1(label)
return
case(0)
exit
end select
if(xyzzyaaab6==0)exit
if(xyzzyaaab6==xyzzyaacb1(item%value))exit
if(item%value%chars(xyzzyaaab6+1)==' ')exit
xyzzyaaac6=xyzzyaaab6
enddo
call xyzzyaabx1(item%value,xyzzyaaae6%value)
if(xyzzyaaab6<1)then
call xyzzyaabw1(item%value)
else
call xyzzyaacj1(xyzzyaaae6%value,1,xyzzyaaab6-1)
call xyzzyaacj1(item%value,xyzzyaaab6+1)
endif
xyzzyaaab6=xyzzyaacl1(xyzzyaaae6%value,':')
if(xyzzyaaab6>0)then
call xyzzyaabx1(xyzzyaaae6%value,xyzzyaaae6%name)
call xyzzyaacj1(xyzzyaaae6%name,1,xyzzyaaab6-1)
call xyzzyaack1(xyzzyaaae6%name)
if(any(xyzzyaaae6%name%chars(:)==':'))then
xyzzyaaaa6='Missing space or misnamed item at line '//trim(i2s(xyzzyaa&
&al5%iline))//' of file '//trim(filename)//': an item name cannot cont&
&ain ":".'
return
endif
if(xyzzyaaae6%name%chars(1)=='%')then
if(xyzzyaaca1(xyzzyaaae6%name)>1)then
if(xyzzyaaae6%name%chars(2)=='!')then
call xyzzyaacj1(xyzzyaaae6%name,3)
call xyzzyaack1(xyzzyaaae6%name)
call xyzzyaacg1(xyzzyaaae6%name,'%delete'//trim(i2s(xyzzyaaae5))//'-')
if(.not.associated(xyzzyaaan5))then
allocate(xyzzyaaan5)
else
allocate(xyzzyaaan5%next)
xyzzyaaan5%next%prev=>xyzzyaaan5
xyzzyaaan5=>xyzzyaaan5%next
endif
xyzzyaaan5%item=>xyzzyaaae6
xyzzyaaae5=xyzzyaaae5+1
else
call xyzzyaaar1(xyzzyaaae6,label)
xyzzyaaaa6='Illegal item name "'//xyzzyaaci1(label,xyzzyaaca1(label))/&
&/'": character "%" at the beginning of the name is reserved for direc&
&tives, and this name does not match any known directive.'
call xyzzyaabw1(label)
return
endif
else
call xyzzyaaar1(xyzzyaaae6,label)
xyzzyaaaa6='Illegal item name "'//xyzzyaaci1(label,xyzzyaaca1(label))/&
&/'": character "%" at the beginning of the name is reserved for direc&
&tives, and this name does not match any known directive.'
call xyzzyaabw1(label)
return
endif
endif
call xyzzyaaby1(trim(unique_casl_string(xyzzyaaci1(xyzzyaaae6%name,xyz&
&zyaaca1(xyzzyaaae6%name)))),xyzzyaaae6%unique_name)
call xyzzyaabx1(xyzzyaaae6%unique_name,xyzzyaaae6%full_unique_name)
call xyzzyaacg1(xyzzyaaae6%full_unique_name,':')
call xyzzyaach1(xyzzyaaae6%full_unique_name,xyzzyaaae6%parent%full_uni&
&que_name)
call xyzzyaacj1(xyzzyaaae6%value,xyzzyaaab6+1)
call xyzzyaack1(xyzzyaaae6%value)
else
xyzzyaaae6%is_implicit=.true.
item%nimplicit=item%nimplicit+1
call xyzzyaaby1('%u'//trim(i2s(item%nimplicit)),xyzzyaaae6%name)
call xyzzyaaby1('%u'//trim(i2s(item%nimplicit)),xyzzyaaae6%unique_name&
&)
call xyzzyaabx1(xyzzyaaae6%unique_name,xyzzyaaae6%full_unique_name)
call xyzzyaacg1(xyzzyaaae6%full_unique_name,':')
call xyzzyaach1(xyzzyaaae6%full_unique_name,xyzzyaaae6%parent%full_uni&
&que_name)
call xyzzyaack1(xyzzyaaae6%value)
if(xyzzyaacb1(xyzzyaaae6%value)==0)then
call xyzzyaaar1(item,label)
xyzzyaaaa6='Syntax error in contents of "'//xyzzyaaci1(label,xyzzyaaca&
&1(label))//'": unnamed scalar with zero-sized value encountered (i.e.&
&, empty item #'//trim(i2s(item%nimplicit))//' in comma-separated list&
&).'
call xyzzyaabw1(label)
return
endif
endif
call xyzzyaabq1(xyzzyaaae6,xyzzyaaad6)
if(xyzzyaaad6)then
call xyzzyaaar1(xyzzyaaae6,label)
if(.not.xyzzyaaae6%is_implicit)then
xyzzyaaaa6='Item '//xyzzyaaci1(label,xyzzyaaca1(label))//' found twice&
&.'
else
xyzzyaaaa6='Problem with generation of internal names for implicitly n&
&amed items: label "'//xyzzyaaci1(label,xyzzyaaca1(label))//'" generat&
&ed twice. This is a bug.'
endif
call xyzzyaabw1(label)
return
endif
if(.not.xyzzyaaae6%is_implicit)then
if(xyzzyaaca1(xyzzyaaae6%value)>0)then
if(xyzzyaaae6%value%chars(1)=='[')then
xyzzyaaae6%is_block=.true.
xyzzyaaae6%is_inline=.true.
call xyzzyaaao5(xyzzyaaae6,xyzzyaaaa6)
if(len_trim(xyzzyaaaa6)>0)return
endif
endif
endif
if(xyzzyaaca1(item%value)==0)exit
enddo
call xyzzyaabw1(item%value)
end subroutine xyzzyaaao5
recursive function xyzzyaaap5(line,i0,look_for,ignore_other) result(sy&
&n_parse)
implicit none
integer,intent(in) :: i0
logical,intent(in),optional :: ignore_other
character(*),intent(in) :: line
character,intent(in),optional :: look_for
integer syn_parse,i,n
logical no_syntax
character left_delim,right_delim,c1
syn_parse=-1
no_syntax=.false.
if(present(ignore_other))no_syntax=ignore_other
if(.not.present(look_for))then
if(i0<1)then
syn_parse=-3
return
endif
left_delim=line(i0:i0)
select case(left_delim)
case('[')
right_delim=']'
case('(')
right_delim=')'
case('{')
right_delim='}'
case('"')
right_delim='"'
case default
syn_parse=-3
return
end select
else
right_delim=look_for
endif
i=i0
n=len(line)
do
i=i+1
if(i>n)exit
c1=line(i:i)
if(c1==right_delim)then
syn_parse=i
return
elseif(.not.no_syntax)then
select case(c1)
case('[','(','{')
i=xyzzyaaap5(line,i)
case('"')
i=xyzzyaaap5(line,i,ignore_other=.true.)
case(']',')','}')
syn_parse=-2
return
case default
cycle
end select
if(i<0)then
syn_parse=i
return
elseif(i==0)then
syn_parse=-1
return
endif
endif
enddo
syn_parse=0
end function xyzzyaaap5
subroutine xyzzyaaaq5
implicit none
type(virtual_line),pointer :: xyzzyaaaa7,xyzzyaaab7
if(associated(xyzzyaaak5))then
xyzzyaaaa7=>xyzzyaaak5%next
do while(.not.associated(xyzzyaaaa7,xyzzyaaak5))
xyzzyaaab7=>xyzzyaaaa7%next
call xyzzyaabw1(xyzzyaaaa7%line)
deallocate(xyzzyaaaa7)
xyzzyaaaa7=>xyzzyaaab7
enddo
call xyzzyaabw1(xyzzyaaak5%line)
deallocate(xyzzyaaak5)
nullify(xyzzyaaak5)
endif
end subroutine xyzzyaaaq5
subroutine xyzzyaaar5(clean_only)
implicit none
logical,optional :: clean_only
logical xyzzyaaaa8
if(associated(xyzzyaaan5))then
xyzzyaaaa8=.false.
if(present(clean_only))xyzzyaaaa8=clean_only
do
if(.not.xyzzyaaaa8)call xyzzyaabg1(xyzzyaaan5%item)
if(associated(xyzzyaaan5%prev))then
xyzzyaaan5=>xyzzyaaan5%prev
deallocate(xyzzyaaan5%next)
else
deallocate(xyzzyaaan5)
nullify(xyzzyaaan5)
exit
endif
enddo
endif
end subroutine xyzzyaaar5
subroutine xyzzyaaas5
implicit none
call xyzzyaaaq5
call xyzzyaaar5(clean_only=.true.)
if(associated(xyzzyaaai5))call xyzzyaabg1(xyzzyaaai5)
end subroutine xyzzyaaas5
end subroutine read_casl
subroutine write_casl(label,file,errmsg)
implicit none
character(*),intent(in) :: label,file
character(512),intent(inout) :: errmsg
integer xyzzyaaaa10,xyzzyaaab10,xyzzyaaac10
logical xyzzyaaad10
type(casl_item),pointer :: xyzzyaaae10
type(metachar) line
errmsg=''
call xyzzyaaaf1
call xyzzyaabv1(label,xyzzyaaae10)
if(.not.associated(xyzzyaaae10))return
if(.not.associated(xyzzyaaae10%first_child))return
if(.not.associated(xyzzyaaae10%parent,xyzzyaaaa1))return
xyzzyaaaa10=9
xyzzyaaad10=.true.
do while(xyzzyaaad10)
xyzzyaaaa10=xyzzyaaaa10+1
inquire(xyzzyaaaa10,opened=xyzzyaaad10)
if(xyzzyaaaa10>=99)then
errmsg='Could not find free i/o unit.'
return
endif
enddo
open(xyzzyaaaa10,file=trim(file),status='replace',iostat=xyzzyaaab10)
if(xyzzyaaab10/=0)then
errmsg='Problem opening '//trim(file)//'.'
return
endif
xyzzyaaac10=xyzzyaaae10%ilevel+1
xyzzyaaae10=>xyzzyaaae10%first_child
do
call xyzzyaaaf10(xyzzyaaae10,xyzzyaaaa10,xyzzyaaac10,line,.false.,0,.f&
&alse.)
if(.not.associated(xyzzyaaae10%next))exit
write(xyzzyaaaa10,'(a)')''
xyzzyaaae10=>xyzzyaaae10%next
enddo
close(xyzzyaaaa10)
contains
recursive subroutine xyzzyaaaf10(xyzzyaaaa11,xyzzyaaab11,xyzzyaaac11,x&
&yzzyaaag11,xyzzyaaae11,xyzzyaaad11,xyzzyaaaf11)
implicit none
type(casl_item),pointer :: xyzzyaaaa11
integer,intent(in) :: xyzzyaaab11,xyzzyaaac11,xyzzyaaad11
logical,intent(in) :: xyzzyaaae11,xyzzyaaaf11
type(metachar),intent(inout) :: xyzzyaaag11
type(casl_item),pointer :: xyzzyaaah11
if(.not.xyzzyaaae11)then
call xyzzyaabw1(xyzzyaaag11)
if(xyzzyaaaa11%is_block)then
if(associated(xyzzyaaaa11%first_child))then
if(.not.xyzzyaaaa11%is_inline)then
write(xyzzyaaab11,'(a)')repeat(' ',xyzzyaaad1*(xyzzyaaaa11%ilevel-xyzz&
&yaaac11))//xyzzyaaci1(xyzzyaaaa11%name,xyzzyaaca1(xyzzyaaaa11%name))/&
&/':'
xyzzyaaah11=>xyzzyaaaa11%first_child
do
call xyzzyaaaf10(xyzzyaaah11,xyzzyaaab11,xyzzyaaac11,xyzzyaaag11,.fals&
&e.,0,.false.)
if(.not.associated(xyzzyaaah11%next))exit
xyzzyaaah11=>xyzzyaaah11%next
enddo
else
call xyzzyaaby1(repeat(' ',xyzzyaaad1*(xyzzyaaaa11%ilevel-xyzzyaaac11)&
&)//xyzzyaaci1(xyzzyaaaa11%name,xyzzyaaca1(xyzzyaaaa11%name))//': [ ',&
&xyzzyaaag11)
xyzzyaaah11=>xyzzyaaaa11%first_child
do
call xyzzyaaaf10(xyzzyaaah11,xyzzyaaab11,xyzzyaaac11,xyzzyaaag11,.true&
&.,xyzzyaaaa11%ilevel,associated(xyzzyaaah11%next))
if(.not.associated(xyzzyaaah11%next))exit
xyzzyaaah11=>xyzzyaaah11%next
enddo
call xyzzyaace1(xyzzyaaag11,' ]')
write(xyzzyaaab11,'(a)')xyzzyaaci1(xyzzyaaag11,xyzzyaaca1(xyzzyaaag11)&
&)
call xyzzyaabw1(xyzzyaaag11)
endif
else
write(xyzzyaaab11,'(a)')repeat(' ',xyzzyaaad1*(xyzzyaaaa11%ilevel-xyzz&
&yaaac11))//xyzzyaaci1(xyzzyaaaa11%name,xyzzyaaca1(xyzzyaaaa11%name))/&
&/': [ ]'
endif
else
if(.not.xyzzyaaaa11%is_implicit)then
write(xyzzyaaab11,'(a)')repeat(' ',xyzzyaaad1*(xyzzyaaaa11%ilevel-xyzz&
&yaaac11))//xyzzyaaci1(xyzzyaaaa11%name,xyzzyaaca1(xyzzyaaaa11%name))/&
&/': '//xyzzyaaci1(xyzzyaaaa11%value,xyzzyaaca1(xyzzyaaaa11%value))
else
write(xyzzyaaab11,'(a)')repeat(' ',xyzzyaaad1*(xyzzyaaaa11%ilevel-xyzz&
&yaaac11))//xyzzyaaci1(xyzzyaaaa11%value,xyzzyaaca1(xyzzyaaaa11%value)&
&)
endif
endif
else
if(xyzzyaaaa11%is_block)then
call xyzzyaaag10(xyzzyaaag11,xyzzyaaci1(xyzzyaaaa11%name,xyzzyaaca1(xy&
&zzyaaaa11%name))//': [ ',xyzzyaaad11-xyzzyaaac11)
if(associated(xyzzyaaaa11%first_child))then
xyzzyaaah11=>xyzzyaaaa11%first_child
do
call xyzzyaaaf10(xyzzyaaah11,xyzzyaaab11,xyzzyaaac11,xyzzyaaag11,.true&
&.,xyzzyaaad11,associated(xyzzyaaah11%next))
if(.not.associated(xyzzyaaah11%next))exit
xyzzyaaah11=>xyzzyaaah11%next
enddo
endif
if(xyzzyaaaf11)then
call xyzzyaaag10(xyzzyaaag11,' ], ',xyzzyaaad11-xyzzyaaac11)
else
call xyzzyaaag10(xyzzyaaag11,' ]',xyzzyaaad11-xyzzyaaac11)
endif
else
if(.not.xyzzyaaaa11%is_implicit)then
if(xyzzyaaaf11)then
call xyzzyaaag10(xyzzyaaag11,xyzzyaaci1(xyzzyaaaa11%name,xyzzyaaca1(xy&
&zzyaaaa11%name))//': '//xyzzyaaci1(xyzzyaaaa11%value,xyzzyaaca1(xyzzy&
&aaaa11%value))//', ',xyzzyaaad11-xyzzyaaac11)
else
call xyzzyaaag10(xyzzyaaag11,xyzzyaaci1(xyzzyaaaa11%name,xyzzyaaca1(xy&
&zzyaaaa11%name))//': '//xyzzyaaci1(xyzzyaaaa11%value,xyzzyaaca1(xyzzy&
&aaaa11%value)),xyzzyaaad11-xyzzyaaac11)
endif
else
if(xyzzyaaaf11)then
call xyzzyaaag10(xyzzyaaag11,xyzzyaaci1(xyzzyaaaa11%value,xyzzyaaca1(x&
&yzzyaaaa11%value))//', ',xyzzyaaad11-xyzzyaaac11)
else
call xyzzyaaag10(xyzzyaaag11,xyzzyaaci1(xyzzyaaaa11%value,xyzzyaaca1(x&
&yzzyaaaa11%value)),xyzzyaaad11-xyzzyaaac11)
endif
endif
endif
endif
end subroutine xyzzyaaaf10
subroutine xyzzyaaag10(line,string,indent)
implicit none
integer,intent(in) :: indent
character(*),intent(in) :: string
type(metachar),intent(inout) :: line
integer xyzzyaaaa12,xyzzyaaab12,xyzzyaaac12,xyzzyaaad12,xyzzyaaae12
integer,parameter :: xyzzyaaaf12=79
character(len(string)) adjustl_string
xyzzyaaac12=xyzzyaaca1(line)
xyzzyaaad12=len(string)
xyzzyaaae12=xyzzyaaac12+xyzzyaaad12
if(xyzzyaaae12>xyzzyaaaf12)then
write(xyzzyaaaa10,'(a)')xyzzyaaci1(line,xyzzyaaac12)
adjustl_string=adjustl(string)
xyzzyaaaa12=len_trim(string)
xyzzyaaab12=len_trim(adjustl_string)
if(xyzzyaaaa12>xyzzyaaab12)then
call xyzzyaaby1(repeat(' ',xyzzyaaad1*indent+xyzzyaaae1)//adjustl_stri&
&ng(1:xyzzyaaad12-xyzzyaaaa12+xyzzyaaab12),line)
else
call xyzzyaaby1(repeat(' ',xyzzyaaad1*indent+xyzzyaaae1)//string,line)
endif
else
call xyzzyaace1(line,string)
endif
end subroutine xyzzyaaag10
end subroutine write_casl
function unique_casl_string(string) result(unique_string)
implicit none
character(*),intent(in) :: string
integer xyzzyaaaa13,xyzzyaaab13,xyzzyaaac13
character(len(string)) unique_string
do xyzzyaaaa13=1,len_trim(string)
xyzzyaaac13=ichar(string(xyzzyaaaa13:xyzzyaaaa13))
if(xyzzyaaac13>=ichar('A').and.xyzzyaaac13<=ichar('Z'))xyzzyaaac13=xyz&
&zyaaac13+ichar('a')-ichar('A')
unique_string(xyzzyaaaa13:xyzzyaaaa13)=char(xyzzyaaac13)
enddo
xyzzyaaab13=scan(trim(unique_string),' ')
do while(xyzzyaaab13>0)
unique_string(xyzzyaaab13:)=trim(unique_string(xyzzyaaab13+1:))
xyzzyaaab13=scan(trim(unique_string),' ')
enddo
end function unique_casl_string
subroutine push_casl_context(label,respect_spelling)
implicit none
logical,intent(in),optional :: respect_spelling
character(*),intent(in) :: label
logical xyzzyaaaa14,xyzzyaaab14
character(len(label)) crumb
type(casl_item),pointer :: item
xyzzyaaaa14=.false.
if(present(respect_spelling))xyzzyaaaa14=respect_spelling
call xyzzyaaaf1
call xyzzyaabv1(label,item)
xyzzyaaab14=associated(item)
if(xyzzyaaab14)xyzzyaaab14=item%is_block
if(xyzzyaaab14)then
allocate(xyzzyaaac1%next)
xyzzyaaac1%next%prev=>xyzzyaaac1
xyzzyaaac1%next%item=>item
xyzzyaaac1=>xyzzyaaac1%next
if(.not.xyzzyaaaa14.and.len_trim(label)>0)then
call xyzzyaaaq1(trim(label),crumb)
call xyzzyaaby1(trim(crumb),item%name)
endif
endif
end subroutine push_casl_context
subroutine pop_casl_context()
implicit none
call xyzzyaaaf1
if(associated(xyzzyaaac1%prev))then
xyzzyaaac1=>xyzzyaaac1%prev
deallocate(xyzzyaaac1%next)
nullify(xyzzyaaac1%next)
endif
end subroutine pop_casl_context
subroutine xyzzyaaan1(label,crumb)
implicit none
type(metachar),intent(inout) :: label,crumb
integer xyzzyaaaa16
xyzzyaaaa16=xyzzyaacl1(label,':')
call xyzzyaabx1(label,crumb)
if(xyzzyaaaa16>0)then
call xyzzyaacj1(crumb,1,xyzzyaaaa16-1)
call xyzzyaacj1(label,xyzzyaaaa16+1)
else
call xyzzyaabw1(label)
endif
end subroutine xyzzyaaan1
integer function xyzzyaaao1(label)
implicit none
character(*),intent(in) :: label
type(metachar) part_label
integer xyzzyaaaa17
xyzzyaaao1=0
if(len_trim(label)==0)return
xyzzyaaao1=1
call xyzzyaaby1(trim(label),part_label)
do
xyzzyaaaa17=xyzzyaacl1(part_label,':')
if(xyzzyaaaa17<1)exit
xyzzyaaao1=xyzzyaaao1+1
call xyzzyaacj1(part_label,xyzzyaaaa17+1)
enddo
call xyzzyaabw1(part_label)
end function xyzzyaaao1
subroutine xyzzyaaap1(label)
implicit none
type(metachar),intent(inout) :: label
integer xyzzyaaaa18
xyzzyaaaa18=xyzzyaacl1(label,':',back=.true.)
if(xyzzyaaaa18>0)then
call xyzzyaacj1(label,1,xyzzyaaaa18-1)
else
call xyzzyaabw1(label)
endif
end subroutine xyzzyaaap1
subroutine xyzzyaaaq1(label,crumb)
implicit none
character(*),intent(in) :: label
character(len(label)),intent(out) :: crumb
integer xyzzyaaaa19
xyzzyaaaa19=scan(label,':',back=.true.)
if(xyzzyaaaa19>0)then
crumb=label(xyzzyaaaa19+1:)
else
crumb=label
endif
end subroutine xyzzyaaaq1
subroutine xyzzyaaar1(item,label,full,unique)
implicit none
type(casl_item),pointer :: item
type(metachar),intent(inout) :: label
logical,intent(in),optional :: full,unique
logical xyzzyaaaa20,xyzzyaaab20
type(casl_item),pointer :: xyzzyaaac20
xyzzyaaaa20=.true.
if(present(full))xyzzyaaaa20=.not.full
xyzzyaaab20=.false.
if(present(unique))xyzzyaaab20=unique
call xyzzyaabw1(label)
call xyzzyaaaf1
if(.not.associated(item))return
xyzzyaaac20=>item
do
if((xyzzyaaaa20.and.associated(xyzzyaaac1%item,xyzzyaaac20)).or.associ&
&ated(xyzzyaaaa1,xyzzyaaac20))exit
if(xyzzyaaca1(label)==0)then
if(.not.xyzzyaaab20)then
call xyzzyaabx1(xyzzyaaac20%name,label)
else
call xyzzyaabx1(xyzzyaaac20%unique_name,label)
endif
else
if(.not.xyzzyaaab20)then
call xyzzyaacg1(label,trim(xyzzyaaci1(xyzzyaaac20%name,xyzzyaaca1(xyzz&
&yaaac20%name))//':'))
else
call xyzzyaacg1(label,trim(xyzzyaaci1(xyzzyaaac20%unique_name,xyzzyaac&
&a1(xyzzyaaac20%unique_name))//':'))
endif
endif
xyzzyaaac20=>xyzzyaaac20%parent
enddo
end subroutine xyzzyaaar1
subroutine xyzzyaaas1(label,value,ierr,respect_spelling)
implicit none
integer,intent(out) :: ierr
real(kind(1.d0)),intent(out) :: value
logical,intent(in),optional :: respect_spelling
character(*),intent(in) :: label
logical xyzzyaaaa21
character(len(label)) crumb
character(80) tempr
type(casl_item),pointer :: xyzzyaaab21
xyzzyaaaa21=.false.
if(present(respect_spelling))xyzzyaaaa21=respect_spelling
value=0.d0
ierr=-1
call xyzzyaabv1(label,xyzzyaaab21)
if(.not.associated(xyzzyaaab21))return
if(.not.xyzzyaaaa21.and.len_trim(label)>0)then
call xyzzyaaaq1(trim(label),crumb)
call xyzzyaaby1(trim(crumb),xyzzyaaab21%name)
endif
call xyzzyaaat1(xyzzyaaci1(xyzzyaaab21%value,xyzzyaaca1(xyzzyaaab21%va&
&lue)),value,ierr)
call xyzzyaabe1(xyzzyaaab21)
write(tempr,*)value
call xyzzyaaby1(trim(adjustl(tempr)),xyzzyaaab21%value)
end subroutine xyzzyaaas1
subroutine xyzzyaaat1(string,value,ierr)
implicit none
integer,intent(out) :: ierr
real(kind(1.d0)),intent(out) :: value
character(*),intent(in) :: string
real(kind(1.d0)) xyzzyaaaa22
value=0.d0
ierr=0
read(string,*,iostat=ierr)xyzzyaaaa22
if(ierr==0)value=xyzzyaaaa22
end subroutine xyzzyaaat1
subroutine xyzzyaaau1(label,value,ierr,respect_spelling)
implicit none
integer,intent(out) :: ierr
real(kind(1.0)),intent(out) :: value
logical,intent(in),optional :: respect_spelling
character(*),intent(in) :: label
logical xyzzyaaaa23
character(len(label)) crumb
character(80) tempr
type(casl_item),pointer :: xyzzyaaab23
xyzzyaaaa23=.false.
if(present(respect_spelling))xyzzyaaaa23=respect_spelling
value=0.
ierr=-1
call xyzzyaabv1(label,xyzzyaaab23)
if(.not.associated(xyzzyaaab23))return
if(.not.xyzzyaaaa23.and.len_trim(label)>0)then
call xyzzyaaaq1(trim(label),crumb)
call xyzzyaaby1(trim(crumb),xyzzyaaab23%name)
endif
call xyzzyaaav1(xyzzyaaci1(xyzzyaaab23%value,xyzzyaaca1(xyzzyaaab23%va&
&lue)),value,ierr)
call xyzzyaabe1(xyzzyaaab23)
write(tempr,*)value
call xyzzyaaby1(trim(adjustl(tempr)),xyzzyaaab23%value)
end subroutine xyzzyaaau1
subroutine xyzzyaaav1(string,value,ierr)
implicit none
integer,intent(out) :: ierr
real(kind(1.0)),intent(out) :: value
character(*),intent(in) :: string
real(kind(1.0)) xyzzyaaaa24
value=0.0
ierr=0
read(string,*,iostat=ierr)xyzzyaaaa24
if(ierr==0)value=xyzzyaaaa24
end subroutine xyzzyaaav1
subroutine xyzzyaaaw1(label,value,ierr,respect_spelling)
implicit none
integer,intent(out) :: ierr
complex(kind(1.d0)),intent(out) :: value
logical,intent(in),optional :: respect_spelling
character(*),intent(in) :: label
logical xyzzyaaaa25
character(len(label)) crumb
character(80) tempr
type(casl_item),pointer :: xyzzyaaab25
xyzzyaaaa25=.false.
if(present(respect_spelling))xyzzyaaaa25=respect_spelling
value=cmplx(0.d0,0.d0,kind(1.d0))
ierr=-1
call xyzzyaabv1(label,xyzzyaaab25)
if(.not.associated(xyzzyaaab25))return
if(.not.xyzzyaaaa25.and.len_trim(label)>0)then
call xyzzyaaaq1(trim(label),crumb)
call xyzzyaaby1(trim(crumb),xyzzyaaab25%name)
endif
call xyzzyaaax1(xyzzyaaci1(xyzzyaaab25%value,xyzzyaaca1(xyzzyaaab25%va&
&lue)),value,ierr)
call xyzzyaabe1(xyzzyaaab25)
write(tempr,*)value
call xyzzyaaby1(trim(adjustl(tempr)),xyzzyaaab25%value)
end subroutine xyzzyaaaw1
subroutine xyzzyaaax1(string,value,ierr)
implicit none
integer,intent(out) :: ierr
complex(kind(1.d0)),intent(out) :: value
character(*),intent(in) :: string
complex(kind(1.d0)) xyzzyaaaa26
value=cmplx(0.d0,0.d0,kind(1.d0))
ierr=0
read(string,*,iostat=ierr)xyzzyaaaa26
if(ierr==0)value=xyzzyaaaa26
end subroutine xyzzyaaax1
subroutine xyzzyaaay1(label,value,ierr,respect_spelling)
implicit none
integer,intent(out) :: ierr
complex(kind(1.0)),intent(out) :: value
logical,intent(in),optional :: respect_spelling
character(*),intent(in) :: label
logical xyzzyaaaa27
character(len(label)) crumb
character(80) tempr
type(casl_item),pointer :: xyzzyaaab27
xyzzyaaaa27=.false.
if(present(respect_spelling))xyzzyaaaa27=respect_spelling
value=cmplx(0.,0.,kind(1.0))
ierr=-1
call xyzzyaabv1(label,xyzzyaaab27)
if(.not.associated(xyzzyaaab27))return
if(.not.xyzzyaaaa27.and.len_trim(label)>0)then
call xyzzyaaaq1(trim(label),crumb)
call xyzzyaaby1(trim(crumb),xyzzyaaab27%name)
endif
call xyzzyaaaz1(xyzzyaaci1(xyzzyaaab27%value,xyzzyaaca1(xyzzyaaab27%va&
&lue)),value,ierr)
call xyzzyaabe1(xyzzyaaab27)
write(tempr,*)value
call xyzzyaaby1(trim(adjustl(tempr)),xyzzyaaab27%value)
end subroutine xyzzyaaay1
subroutine xyzzyaaaz1(string,value,ierr)
implicit none
integer,intent(out) :: ierr
complex(kind(1.0)),intent(out) :: value
character(*),intent(in) :: string
complex(kind(1.0)) xyzzyaaaa28
value=cmplx(0.0,0.0,kind(1.0))
ierr=0
read(string,*,iostat=ierr)xyzzyaaaa28
if(ierr==0)value=xyzzyaaaa28
end subroutine xyzzyaaaz1
subroutine xyzzyaaba1(label,value,ierr,respect_spelling)
implicit none
integer,intent(out) :: ierr,value
logical,intent(in),optional :: respect_spelling
character(*),intent(in) :: label
logical xyzzyaaaa29
character(len(label)) crumb
type(casl_item),pointer :: xyzzyaaab29
xyzzyaaaa29=.false.
if(present(respect_spelling))xyzzyaaaa29=respect_spelling
value=0
ierr=-1
call xyzzyaabv1(label,xyzzyaaab29)
if(.not.associated(xyzzyaaab29))return
if(.not.xyzzyaaaa29.and.len_trim(label)>0)then
call xyzzyaaaq1(trim(label),crumb)
call xyzzyaaby1(trim(crumb),xyzzyaaab29%name)
endif
call xyzzyaabb1(xyzzyaaci1(xyzzyaaab29%value,xyzzyaaca1(xyzzyaaab29%va&
&lue)),value,ierr)
if(ierr/=0)return
call xyzzyaabe1(xyzzyaaab29)
call xyzzyaaby1(trim(i2s(value)),xyzzyaaab29%value)
end subroutine xyzzyaaba1
subroutine xyzzyaabb1(string,value,ierr)
implicit none
integer,intent(out) :: ierr,value
character(*),intent(in) :: string
integer xyzzyaaaa30
value=0
ierr=0
read(string,*,iostat=ierr)xyzzyaaaa30
if(ierr==0)value=xyzzyaaaa30
end subroutine xyzzyaabb1
subroutine xyzzyaabc1(label,value,ierr,respect_spelling)
implicit none
integer,intent(out) :: ierr
logical,intent(in),optional :: respect_spelling
character(*),intent(in) :: label
character(*),intent(out) :: value
logical xyzzyaaaa31
character(len(label)) crumb
type(casl_item),pointer :: xyzzyaaab31
xyzzyaaaa31=.false.
if(present(respect_spelling))xyzzyaaaa31=respect_spelling
value=""
ierr=-1
call xyzzyaabv1(label,xyzzyaaab31)
if(.not.associated(xyzzyaaab31))return
if(.not.xyzzyaaaa31.and.len_trim(label)>0)then
call xyzzyaaaq1(trim(label),crumb)
call xyzzyaaby1(trim(crumb),xyzzyaaab31%name)
endif
ierr=0
call xyzzyaabz1(xyzzyaaab31%value,value)
call xyzzyaabe1(xyzzyaaab31)
end subroutine xyzzyaabc1
subroutine xyzzyaabd1(label,value,ierr,respect_spelling)
implicit none
integer,intent(out) :: ierr
logical,intent(in),optional :: respect_spelling
logical,intent(out) :: value
character(*),intent(in) :: label
logical xyzzyaaaa32
character(len(label)) crumb
type(casl_item),pointer :: xyzzyaaab32
xyzzyaaaa32=.false.
if(present(respect_spelling))xyzzyaaaa32=respect_spelling
value=.false.
ierr=-1
call xyzzyaabv1(label,xyzzyaaab32)
if(.not.associated(xyzzyaaab32))return
if(.not.xyzzyaaaa32.and.len_trim(label)>0)then
call xyzzyaaaq1(trim(label),crumb)
call xyzzyaaby1(trim(crumb),xyzzyaaab32%name)
endif
select case(unique_casl_string(trim(xyzzyaaci1(xyzzyaaab32%value,xyzzy&
&aaca1(xyzzyaaab32%value)))))
case('.true.','true','t','yes','y','1')
value=.true.
case('.false.','false','f','no','n','0')
value=.false.
case default
return
end select
ierr=0
call xyzzyaabe1(xyzzyaaab32)
if(value)then
call xyzzyaaby1('T',xyzzyaaab32%value)
else
call xyzzyaaby1('F',xyzzyaaab32%value)
endif
end subroutine xyzzyaabd1
subroutine xyzzyaabe1(item)
implicit none
type(casl_item),pointer :: item
type(casl_item),pointer :: xyzzyaaaa33
xyzzyaaaa33=>item
do
if(xyzzyaaaa33%been_read)exit
xyzzyaaaa33%been_read=.true.
if(.not.associated(xyzzyaaaa33%parent))exit
xyzzyaaaa33=>xyzzyaaaa33%parent
enddo
end subroutine xyzzyaabe1
subroutine xyzzyaabf1(item)
implicit none
type(casl_item),pointer :: item
type(casl_item),pointer :: xyzzyaaaa34
item%is_inline=.true.
xyzzyaaaa34=>item
do
if(.not.xyzzyaaaa34%is_inline)exit
xyzzyaaaa34%is_inline=.false.
if(.not.associated(xyzzyaaaa34%parent))exit
xyzzyaaaa34=>xyzzyaaaa34%parent
enddo
end subroutine xyzzyaabf1
subroutine query_casl_item(label,exists,is_block,nchildren,respect_spe&
&lling)
implicit none
integer,intent(out),optional :: nchildren
logical,intent(in),optional :: respect_spelling
logical,intent(out),optional :: exists,is_block
character(*),intent(in) :: label
logical xyzzyaaaa35
character(len(label)) crumb
type(casl_item),pointer :: xyzzyaaab35
call xyzzyaabv1(label,xyzzyaaab35)
if(.not.associated(xyzzyaaab35))then
if(present(exists))exists=.false.
if(present(is_block))is_block=.false.
if(present(nchildren))nchildren=0
else
if(present(exists))exists=.true.
if(present(is_block))is_block=xyzzyaaab35%is_block
xyzzyaaaa35=.false.
if(present(respect_spelling))xyzzyaaaa35=respect_spelling
if(.not.xyzzyaaaa35.and.len_trim(label)>0)then
call xyzzyaaaq1(label,crumb)
call xyzzyaaby1(trim(crumb),xyzzyaaab35%name)
endif
if(present(nchildren))then
nchildren=0
if(xyzzyaaab35%is_block.and.associated(xyzzyaaab35%first_child))then
xyzzyaaab35=>xyzzyaaab35%first_child
do
nchildren=nchildren+1
if(.not.associated(xyzzyaaab35%next))exit
xyzzyaaab35=>xyzzyaaab35%next
enddo
endif
endif
endif
end subroutine query_casl_item
subroutine first_unread_child(label,child_unique_name,ierr,flag_as_rea&
&d)
implicit none
integer,intent(out) :: ierr
logical,intent(in),optional :: flag_as_read
character(*),intent(in) :: label
character(casl_keysize),intent(out) :: child_unique_name
type(casl_item),pointer :: xyzzyaaaa36
child_unique_name=''
ierr=-1
call xyzzyaabv1(label,xyzzyaaaa36)
if(.not.associated(xyzzyaaaa36))return
ierr=-2
if(.not.(xyzzyaaaa36%is_block.and.associated(xyzzyaaaa36%first_child))&
&)return
ierr=-3
xyzzyaaaa36=>xyzzyaaaa36%first_child
do
if(.not.xyzzyaaaa36%been_read)then
ierr=0
call xyzzyaabz1(xyzzyaaaa36%unique_name,child_unique_name)
if(present(flag_as_read))then
if(flag_as_read)xyzzyaaaa36%been_read=.true.
endif
return
endif
if(.not.associated(xyzzyaaaa36%next))exit
xyzzyaaaa36=>xyzzyaaaa36%next
enddo
end subroutine first_unread_child
subroutine check_unread_casl(label,errmsg)
implicit none
character(*),intent(in) :: label
character(*),intent(inout) :: errmsg
integer xyzzyaaaa37
type(casl_item),pointer :: xyzzyaaab37,xyzzyaaac37
xyzzyaaaa37=0
errmsg=''
call xyzzyaabv1(label,xyzzyaaab37)
if(.not.associated(xyzzyaaab37))return
xyzzyaaac37=>xyzzyaaab37
do
if(xyzzyaaac37%is_block.and.associated(xyzzyaaac37%first_child))then
xyzzyaaac37=>xyzzyaaac37%first_child
cycle
endif
if(.not.xyzzyaaac37%been_read)xyzzyaaaa37=xyzzyaaaa37+1
if(associated(xyzzyaaac37,xyzzyaaab37))exit
do while(.not.associated(xyzzyaaac37%next))
xyzzyaaac37=>xyzzyaaac37%parent
if(associated(xyzzyaaac37,xyzzyaaab37))exit
enddo
if(associated(xyzzyaaac37,xyzzyaaab37))exit
xyzzyaaac37=>xyzzyaaac37%next
enddo
if(xyzzyaaaa37>0)errmsg='CASL item "'//label//'" contains '//trim(i2s(&
&xyzzyaaaa37))//' unrecognized entries.'
end subroutine check_unread_casl
subroutine delete_casl_item(label)
character(*),intent(in) :: label
type(casl_item),pointer :: xyzzyaaaa38
call xyzzyaabv1(label,xyzzyaaaa38)
if(.not.associated(xyzzyaaaa38))return
call xyzzyaabg1(xyzzyaaaa38)
end subroutine delete_casl_item
subroutine xyzzyaabg1(item)
implicit none
type(casl_item),pointer :: item
type(casl_item),pointer :: parent,xyzzyaaaa39,xyzzyaaab39
integer xyzzyaaac39,xyzzyaaad39,xyzzyaaae39
logical xyzzyaaaf39
character(32) impl_name
if(item%is_block)then
if(associated(item%first_child))then
xyzzyaaaa39=>item%first_child
do
if(xyzzyaaaa39%is_block.and.associated(xyzzyaaaa39%first_child))then
xyzzyaaaa39=>xyzzyaaaa39%first_child
nullify(xyzzyaaaa39%parent%first_child)
cycle
endif
call xyzzyaabu1(xyzzyaaaa39)
call xyzzyaabw1(xyzzyaaaa39%value)
call xyzzyaabw1(xyzzyaaaa39%name)
call xyzzyaabw1(xyzzyaaaa39%unique_name)
call xyzzyaabw1(xyzzyaaaa39%full_unique_name)
if(associated(xyzzyaaaa39%next))then
xyzzyaaab39=>xyzzyaaaa39%next
else
xyzzyaaab39=>xyzzyaaaa39%parent
endif
deallocate(xyzzyaaaa39)
xyzzyaaaa39=>xyzzyaaab39
if(associated(xyzzyaaaa39,item))exit
enddo
endif
else
call xyzzyaabw1(item%value)
endif
if(associated(item%parent))then
parent=>item%parent
if(associated(parent%first_child,item))then
if(.not.associated(item%next))then
nullify(parent%first_child,parent%last_child)
else
parent%first_child=>item%next
nullify(item%next%prev)
endif
elseif(associated(parent%last_child,item))then
parent%last_child=>item%prev
nullify(item%prev%next)
else
item%prev%next=>item%next
item%next%prev=>item%prev
endif
else
nullify(parent)
endif
call xyzzyaabu1(item)
if(associated(parent))then
call xyzzyaabz1(item%name,impl_name)
if(impl_name(1:2)=='%u')then
parent%nimplicit=parent%nimplicit-1
if(associated(parent%first_child))then
call xyzzyaabb1(trim(impl_name(3:)),xyzzyaaac39,xyzzyaaae39)
if(xyzzyaaac39<=parent%nimplicit)then
xyzzyaaaa39=>parent%first_child
do
call xyzzyaabz1(xyzzyaaaa39%name,impl_name)
if(impl_name(1:2)=='%u')then
call xyzzyaabb1(trim(impl_name(3:)),xyzzyaaad39,xyzzyaaae39)
if(xyzzyaaad39>xyzzyaaac39)then
call xyzzyaabu1(xyzzyaaaa39)
call xyzzyaaby1('%u'//trim(i2s(xyzzyaaad39-1)),xyzzyaaaa39%name)
call xyzzyaaby1('%u'//trim(i2s(xyzzyaaad39-1)),xyzzyaaaa39%unique_name&
&)
call xyzzyaabx1(xyzzyaaaa39%unique_name,xyzzyaaaa39%full_unique_name)
call xyzzyaacg1(xyzzyaaaa39%full_unique_name,':')
call xyzzyaach1(xyzzyaaaa39%full_unique_name,xyzzyaaaa39%parent%full_u&
&nique_name)
call xyzzyaabq1(xyzzyaaaa39,xyzzyaaaf39)
endif
if(xyzzyaaad39==parent%nimplicit)exit
endif
if(.not.associated(xyzzyaaaa39%next))exit
xyzzyaaaa39=>xyzzyaaaa39%next
enddo
endif
endif
endif
endif
call xyzzyaabw1(item%name)
call xyzzyaabw1(item%unique_name)
call xyzzyaabw1(item%full_unique_name)
call xyzzyaabw1(item%value)
deallocate(item)
nullify(item)
end subroutine xyzzyaabg1
subroutine xyzzyaabh1(label,item,errmsg,value,only_if_unset,as_block,r&
&espect_spelling)
implicit none
character(*),intent(in) :: label
type(casl_item),pointer :: item
character(512),intent(inout) :: errmsg
character(*),intent(in),optional :: value
logical,intent(in),optional :: only_if_unset,as_block,respect_spelling
integer xyzzyaaaa40,xyzzyaaab40,xyzzyaaac40,xyzzyaaad40,xyzzyaaae40
logical xyzzyaaaf40,xyzzyaaag40,xyzzyaaah40,xyzzyaaai40
character(len(label)) crumb
type(casl_item),pointer :: parent
type(metachar) mlabel,mlabel2,ulabel,ulabel2
call xyzzyaaaf1
errmsg=''
xyzzyaaag40=.false.
if(present(as_block))xyzzyaaag40=as_block
xyzzyaaah40=.false.
if(present(only_if_unset))xyzzyaaah40=only_if_unset
xyzzyaaai40=.false.
if(present(respect_spelling))xyzzyaaai40=respect_spelling
call xyzzyaabi1(label,mlabel)
call xyzzyaabv1(':'//xyzzyaaci1(mlabel,xyzzyaaca1(mlabel)),item)
if(associated(item))then
if(xyzzyaaag40.neqv.item%is_block)then
errmsg='Block status mismatch for '//label//'.  Check if you have defi&
&ned this item as a scalar/block in the relevant .casl file when it sh&
&ould have been a block/scalar (e.g., you may have tried to define an &
&empty block but forgot the square brackets "[ ]" after the colon, or &
&you may have accidentally indented the line after an empty scalar).'
call xyzzyaabw1(mlabel)
return
endif
if(.not.xyzzyaaai40.and.len_trim(label)>0)then
call xyzzyaaaq1(trim(label),crumb)
call xyzzyaaby1(trim(crumb),item%name)
endif
if(.not.(xyzzyaaah40.or.xyzzyaaag40))call xyzzyaaby1(value,item%value)
call xyzzyaabe1(item)
call xyzzyaabw1(mlabel)
return
endif
call xyzzyaabx1(mlabel,mlabel2)
do while(xyzzyaaca1(mlabel2)>0)
call xyzzyaaap1(mlabel2)
call xyzzyaabv1(':'//xyzzyaaci1(mlabel2,xyzzyaaca1(mlabel2)),parent)
if(associated(parent))exit
enddo
if(xyzzyaaca1(mlabel2)>0)call xyzzyaacj1(mlabel,xyzzyaaca1(mlabel2)+2)
xyzzyaaad40=xyzzyaaao1(xyzzyaaci1(mlabel,xyzzyaaca1(mlabel)))
xyzzyaaae40=xyzzyaaao1(trim(label))
call xyzzyaaby1(trim(label),ulabel)
do xyzzyaaac40=xyzzyaaad40,xyzzyaaae40-1
call xyzzyaaan1(ulabel,ulabel2)
enddo
call xyzzyaaan1(mlabel,mlabel2)
call xyzzyaaan1(ulabel,ulabel2)
if(xyzzyaaca1(mlabel)>0)then
do while(xyzzyaaca1(mlabel)>0)
allocate(item)
call xyzzyaabx1(ulabel2,item%name)
if(associated(parent,xyzzyaaaa1))then
call xyzzyaabx1(ulabel2,item%unique_name)
call xyzzyaabx1(ulabel2,item%full_unique_name)
else
call xyzzyaabx1(mlabel2,item%unique_name)
call xyzzyaabx1(item%unique_name,item%full_unique_name)
call xyzzyaacg1(item%full_unique_name,':')
call xyzzyaach1(item%full_unique_name,parent%full_unique_name)
endif
item%is_block=.true.
item%ilevel=parent%ilevel+1
item%is_inline=parent%is_inline
item%parent=>parent
nullify(item%next)
if(associated(parent%first_child))then
parent%last_child%next=>item
item%prev=>parent%last_child
parent%last_child=>item
else
parent%first_child=>item
parent%last_child=>item
nullify(item%prev)
endif
call xyzzyaabq1(item,xyzzyaaaf40)
if(xyzzyaaaf40)then
errmsg='Failed to insert new item "'//xyzzyaaci1(item%full_unique_name&
&,xyzzyaaca1(item%full_unique_name))//'" in AVL tree since it appears &
&to be already present (ancestor).'
call xyzzyaabw1(mlabel)
call xyzzyaabw1(mlabel2)
call xyzzyaabw1(ulabel)
call xyzzyaabw1(ulabel2)
return
endif
parent=>item
nullify(item)
call xyzzyaaan1(mlabel,mlabel2)
call xyzzyaaan1(ulabel,ulabel2)
enddo
endif
call xyzzyaabw1(mlabel)
call xyzzyaabw1(ulabel)
allocate(item)
call xyzzyaabx1(ulabel2,item%name)
if(associated(parent,xyzzyaaaa1))then
call xyzzyaabx1(ulabel2,item%unique_name)
else
call xyzzyaabx1(mlabel2,item%unique_name)
endif
item%is_block=xyzzyaaag40
item%ilevel=parent%ilevel+1
item%been_read=.true.
item%is_inline=parent%is_inline
item%parent=>parent
if(.not.xyzzyaaag40)then
call xyzzyaabz1(mlabel2,crumb)
if(crumb(1:1)=='%')then
if(crumb(2:2)=='u')then
item%is_implicit=.true.
parent%nimplicit=parent%nimplicit+1
if(trim(crumb)=='%u')then
call xyzzyaace1(mlabel2,trim(i2s(parent%nimplicit)))
call xyzzyaabx1(mlabel2,item%name)
call xyzzyaabx1(mlabel2,item%unique_name)
else
call xyzzyaabb1(crumb(3:),xyzzyaaaa40,xyzzyaaab40)
if(xyzzyaaab40/=0)then
errmsg='Problem reading implicit name index in label.'
call xyzzyaabw1(mlabel2)
call xyzzyaabw1(ulabel2)
return
endif
if(parent%nimplicit/=xyzzyaaaa40)then
errmsg='Problem counting implicitly-named items.'
call xyzzyaabw1(mlabel2)
call xyzzyaabw1(ulabel2)
return
endif
endif
else
errmsg='Cannot create non-implicit item whose name starts with "%".'
call xyzzyaabw1(mlabel2)
call xyzzyaabw1(ulabel2)
return
endif
endif
if(present(value))then
call xyzzyaaby1(trim(value),item%value)
else
call xyzzyaabw1(item%value)
endif
else
if(mlabel2%chars(1)=='%')then
errmsg='Cannot create block item whose name starts with "%".'
call xyzzyaabw1(mlabel2)
call xyzzyaabw1(ulabel2)
return
endif
endif
call xyzzyaabw1(mlabel2)
call xyzzyaabw1(ulabel2)
call xyzzyaabx1(item%unique_name,item%full_unique_name)
if(.not.associated(parent,xyzzyaaaa1))then
call xyzzyaacg1(item%full_unique_name,':')
call xyzzyaach1(item%full_unique_name,parent%full_unique_name)
endif
if(associated(parent%first_child))then
parent%last_child%next=>item
item%prev=>parent%last_child
parent%last_child=>item
else
parent%first_child=>item
parent%last_child=>item
nullify(item%prev)
endif
call xyzzyaabq1(item,xyzzyaaaf40)
if(xyzzyaaaf40)then
errmsg='Failed to insert new item "'//xyzzyaaci1(item%full_unique_name&
&,xyzzyaaca1(item%full_unique_name))//'" in AVL tree since it appears &
&to be already present.'
return
endif
end subroutine xyzzyaabh1
subroutine xyzzyaabi1(label,path)
implicit none
character(*),intent(in) :: label
type(metachar),intent(inout) :: path
type(metachar) tstring,tstring2
call xyzzyaaby1(label,path)
if(xyzzyaaca1(path)>0)then
if(path%chars(1)/=':')then
call xyzzyaacg1(path,':')
call xyzzyaach1(path,xyzzyaaac1%item%full_unique_name)
endif
endif
if(xyzzyaaca1(path)>0)then
do while(path%chars(1)==':')
call xyzzyaacj1(path,2)
if(xyzzyaaca1(path)==0)exit
enddo
endif
if(xyzzyaaca1(path)>0)then
call xyzzyaaan1(path,tstring)
if(xyzzyaaca1(tstring)>0.and.xyzzyaaca1(path)>0)then
call xyzzyaaby1(trim(unique_casl_string(xyzzyaaci1(path,xyzzyaaca1(pat&
&h)))),tstring2)
call xyzzyaacg1(tstring2,':')
call xyzzyaach1(tstring2,tstring)
call xyzzyaabx1(tstring2,path)
call xyzzyaabw1(tstring2)
elseif(xyzzyaaca1(tstring)>0)then
call xyzzyaabx1(tstring,path)
endif
call xyzzyaabw1(tstring)
endif
end subroutine xyzzyaabi1
subroutine xyzzyaabj1(label,value,errmsg,only_if_unset,respect_spellin&
&g)
implicit none
character(*),intent(in) :: label
real(kind(1.d0)),intent(in) :: value
character(512),intent(inout) :: errmsg
logical,intent(in),optional :: only_if_unset,respect_spelling
logical if_unset,xyzzyaaaa42
character(80) tmpr
type(casl_item),pointer :: xyzzyaaab42
errmsg=''
if_unset=.false.
if(present(only_if_unset))if_unset=only_if_unset
xyzzyaaaa42=.false.
if(present(respect_spelling))xyzzyaaaa42=respect_spelling
write(tmpr,*)value
call xyzzyaabh1(label,xyzzyaaab42,errmsg,value=trim(adjustl(tmpr)),onl&
&y_if_unset=if_unset,respect_spelling=xyzzyaaaa42)
end subroutine xyzzyaabj1
subroutine xyzzyaabk1(label,value,errmsg,only_if_unset,respect_spellin&
&g)
implicit none
character(*),intent(in) :: label
real(kind(1.0)),intent(in) :: value
character(512),intent(inout) :: errmsg
logical,intent(in),optional :: only_if_unset,respect_spelling
logical if_unset,xyzzyaaaa43
character(80) tmpr
type(casl_item),pointer :: xyzzyaaab43
errmsg=''
if_unset=.false.
if(present(only_if_unset))if_unset=only_if_unset
xyzzyaaaa43=.false.
if(present(respect_spelling))xyzzyaaaa43=respect_spelling
write(tmpr,*)value
call xyzzyaabh1(label,xyzzyaaab43,errmsg,value=trim(adjustl(tmpr)),onl&
&y_if_unset=if_unset,respect_spelling=xyzzyaaaa43)
end subroutine xyzzyaabk1
subroutine xyzzyaabl1(label,value,errmsg,only_if_unset,respect_spellin&
&g)
implicit none
character(*),intent(in) :: label
complex(kind(1.d0)),intent(in) :: value
character(512),intent(inout) :: errmsg
logical,intent(in),optional :: only_if_unset,respect_spelling
logical if_unset,xyzzyaaaa44
character(80) tmpr
type(casl_item),pointer :: xyzzyaaab44
errmsg=''
if_unset=.false.
if(present(only_if_unset))if_unset=only_if_unset
xyzzyaaaa44=.false.
if(present(respect_spelling))xyzzyaaaa44=respect_spelling
write(tmpr,*)value
call xyzzyaabh1(label,xyzzyaaab44,errmsg,value=trim(adjustl(tmpr)),onl&
&y_if_unset=if_unset,respect_spelling=xyzzyaaaa44)
end subroutine xyzzyaabl1
subroutine xyzzyaabm1(label,value,errmsg,only_if_unset,respect_spellin&
&g)
implicit none
character(*),intent(in) :: label
complex(kind(1.0)),intent(in) :: value
character(512),intent(inout) :: errmsg
logical,intent(in),optional :: only_if_unset,respect_spelling
logical if_unset,xyzzyaaaa45
character(80) tmpr
type(casl_item),pointer :: xyzzyaaab45
errmsg=''
if_unset=.false.
if(present(only_if_unset))if_unset=only_if_unset
xyzzyaaaa45=.false.
if(present(respect_spelling))xyzzyaaaa45=respect_spelling
write(tmpr,*)value
call xyzzyaabh1(label,xyzzyaaab45,errmsg,value=trim(adjustl(tmpr)),onl&
&y_if_unset=if_unset,respect_spelling=xyzzyaaaa45)
end subroutine xyzzyaabm1
subroutine xyzzyaabn1(label,value,errmsg,only_if_unset,respect_spellin&
&g)
implicit none
character(*),intent(in) :: label
integer,intent(in) :: value
character(512),intent(inout) :: errmsg
logical,intent(in),optional :: only_if_unset,respect_spelling
logical if_unset,xyzzyaaaa46
character(80) tmpr
type(casl_item),pointer :: xyzzyaaab46
errmsg=''
if_unset=.false.
if(present(only_if_unset))if_unset=only_if_unset
xyzzyaaaa46=.false.
if(present(respect_spelling))xyzzyaaaa46=respect_spelling
write(tmpr,*)value
call xyzzyaabh1(label,xyzzyaaab46,errmsg,value=trim(adjustl(tmpr)),onl&
&y_if_unset=if_unset,respect_spelling=xyzzyaaaa46)
end subroutine xyzzyaabn1
subroutine xyzzyaabo1(label,value,errmsg,only_if_unset,respect_spellin&
&g)
implicit none
character(*),intent(in) :: label,value
character(512),intent(inout) :: errmsg
logical,intent(in),optional :: only_if_unset,respect_spelling
logical if_unset,xyzzyaaaa47
type(casl_item),pointer :: xyzzyaaab47
errmsg=''
if_unset=.false.
if(present(only_if_unset))if_unset=only_if_unset
xyzzyaaaa47=.false.
if(present(respect_spelling))xyzzyaaaa47=respect_spelling
call xyzzyaabh1(label,xyzzyaaab47,errmsg,value=value,only_if_unset=if_&
&unset,respect_spelling=xyzzyaaaa47)
end subroutine xyzzyaabo1
subroutine xyzzyaabp1(label,value,errmsg,only_if_unset,respect_spellin&
&g)
implicit none
character(*),intent(in) :: label
logical,intent(in) :: value
character(512),intent(inout) :: errmsg
logical,intent(in),optional :: only_if_unset,respect_spelling
logical if_unset,xyzzyaaaa48
character(1) char1
type(casl_item),pointer :: xyzzyaaab48
errmsg=''
if_unset=.false.
if(present(only_if_unset))if_unset=only_if_unset
xyzzyaaaa48=.false.
if(present(respect_spelling))xyzzyaaaa48=respect_spelling
if(value)then
char1='T'
else
char1='F'
endif
call xyzzyaabh1(label,xyzzyaaab48,errmsg,value=char1,only_if_unset=if_&
&unset,respect_spelling=xyzzyaaaa48)
end subroutine xyzzyaabp1
subroutine set_casl_block(label,errmsg,prefer_inline,force_noninline,r&
&espect_spelling,push)
implicit none
character(*),intent(in) :: label
character(512),intent(inout) :: errmsg
logical,intent(in),optional :: prefer_inline,force_noninline,respect_s&
&pelling,push
logical xyzzyaaaa49
type(casl_item),pointer :: item
errmsg=''
xyzzyaaaa49=.false.
if(present(respect_spelling))xyzzyaaaa49=respect_spelling
call xyzzyaabh1(label,item,errmsg,as_block=.true.,respect_spelling=xyz&
&zyaaaa49)
if(present(prefer_inline))then
item%is_inline=prefer_inline
if(associated(item%parent))then
if(item%parent%is_inline)item%is_inline=.true.
endif
endif
if(present(force_noninline))then
if(force_noninline)call xyzzyaabf1(item)
endif
if(present(push))then
if(push)then
allocate(xyzzyaaac1%next)
xyzzyaaac1%next%prev=>xyzzyaaac1
xyzzyaaac1%next%item=>item
xyzzyaaac1=>xyzzyaaac1%next
endif
endif
end subroutine set_casl_block
subroutine xyzzyaabq1(item,already_present)
implicit none
type(casl_item),pointer :: item
logical,intent(out) :: already_present
type(casl_item),pointer :: xyzzyaaaa50
already_present=.false.
xyzzyaaaa50=>xyzzyaaab1
do
select case(compare_metachar(item%full_unique_name,xyzzyaaaa50%full_un&
&ique_name))
case('=')
already_present=.true.
return
case('>')
if(associated(xyzzyaaaa50%avl_right))then
xyzzyaaaa50=>xyzzyaaaa50%avl_right
else
xyzzyaaaa50%avl_right=>item
item%avl_parent=>xyzzyaaaa50
item%avl_depth=1
call xyzzyaabr1(item)
return
endif
case('<')
if(associated(xyzzyaaaa50%avl_left))then
xyzzyaaaa50=>xyzzyaaaa50%avl_left
else
xyzzyaaaa50%avl_left=>item
item%avl_parent=>xyzzyaaaa50
item%avl_depth=1
call xyzzyaabr1(item)
return
endif
end select
enddo
end subroutine xyzzyaabq1
subroutine xyzzyaabr1(item)
implicit none
type(casl_item),pointer :: item
type(casl_item),pointer :: xyzzyaaaa51
xyzzyaaaa51=>item
do
call xyzzyaabt1(xyzzyaaaa51)
call xyzzyaabs1(xyzzyaaaa51)
if(.not.associated(xyzzyaaaa51%avl_parent))then
xyzzyaaab1=>xyzzyaaaa51
exit
endif
xyzzyaaaa51=>xyzzyaaaa51%avl_parent
enddo
end subroutine xyzzyaabr1
subroutine xyzzyaabs1(item)
implicit none
type(casl_item),pointer :: item
type(casl_item),pointer :: xyzzyaaaa52,xyzzyaaab52
integer xyzzyaaac52,xyzzyaaad52
xyzzyaaac52=item%avl_imbalance
if(xyzzyaaac52<-1)then
xyzzyaaaa52=>item%avl_left
xyzzyaaad52=xyzzyaaaa52%avl_imbalance
if(xyzzyaaad52<0)then
if(associated(xyzzyaaaa52%avl_right))then
item%avl_left=>xyzzyaaaa52%avl_right
item%avl_left%avl_parent=>item
else
nullify(item%avl_left)
endif
xyzzyaaaa52%avl_right=>item
if(associated(item%avl_parent))then
xyzzyaaaa52%avl_parent=>item%avl_parent
if(associated(item%avl_parent%avl_left,item))then
item%avl_parent%avl_left=>xyzzyaaaa52
else
item%avl_parent%avl_right=>xyzzyaaaa52
endif
else
nullify(xyzzyaaaa52%avl_parent)
endif
item%avl_parent=>xyzzyaaaa52
call xyzzyaabt1(item)
call xyzzyaabt1(xyzzyaaaa52)
item=>xyzzyaaaa52
elseif(xyzzyaaad52>0)then
xyzzyaaab52=>xyzzyaaaa52%avl_right
if(associated(xyzzyaaab52%avl_left))then
xyzzyaaaa52%avl_right=>xyzzyaaab52%avl_left
xyzzyaaaa52%avl_right%avl_parent=>xyzzyaaaa52
else
nullify(xyzzyaaaa52%avl_right)
endif
if(associated(xyzzyaaab52%avl_right))then
item%avl_left=>xyzzyaaab52%avl_right
item%avl_left%avl_parent=>item
else
nullify(item%avl_left)
endif
xyzzyaaab52%avl_left=>xyzzyaaaa52
xyzzyaaab52%avl_right=>item
if(associated(item%avl_parent))then
xyzzyaaab52%avl_parent=>item%avl_parent
else
nullify(xyzzyaaab52%avl_parent)
endif
if(associated(item%avl_parent))then
xyzzyaaab52%avl_parent=>item%avl_parent
if(associated(item%avl_parent%avl_left,item))then
xyzzyaaab52%avl_parent%avl_left=>xyzzyaaab52
else
xyzzyaaab52%avl_parent%avl_right=>xyzzyaaab52
endif
else
nullify(xyzzyaaaa52%avl_parent)
endif
xyzzyaaab52%avl_left=>xyzzyaaaa52
xyzzyaaaa52%avl_parent=>xyzzyaaab52
xyzzyaaab52%avl_right=>item
item%avl_parent=>xyzzyaaab52
call xyzzyaabt1(xyzzyaaaa52)
call xyzzyaabt1(item)
call xyzzyaabt1(xyzzyaaab52)
item=>xyzzyaaab52
endif
elseif(xyzzyaaac52>1)then
xyzzyaaaa52=>item%avl_right
xyzzyaaad52=xyzzyaaaa52%avl_imbalance
if(xyzzyaaad52<0)then
xyzzyaaab52=>xyzzyaaaa52%avl_left
if(associated(xyzzyaaab52%avl_right))then
xyzzyaaaa52%avl_left=>xyzzyaaab52%avl_right
xyzzyaaaa52%avl_left%avl_parent=>xyzzyaaaa52
else
nullify(xyzzyaaaa52%avl_left)
endif
if(associated(xyzzyaaab52%avl_left))then
item%avl_right=>xyzzyaaab52%avl_left
item%avl_right%avl_parent=>item
else
nullify(item%avl_right)
endif
xyzzyaaab52%avl_right=>xyzzyaaaa52
xyzzyaaab52%avl_left=>item
if(associated(item%avl_parent))then
xyzzyaaab52%avl_parent=>item%avl_parent
else
nullify(xyzzyaaab52%avl_parent)
endif
if(associated(item%avl_parent))then
xyzzyaaab52%avl_parent=>item%avl_parent
if(associated(item%avl_parent%avl_right,item))then
xyzzyaaab52%avl_parent%avl_right=>xyzzyaaab52
else
xyzzyaaab52%avl_parent%avl_left=>xyzzyaaab52
endif
else
nullify(xyzzyaaaa52%avl_parent)
endif
xyzzyaaab52%avl_right=>xyzzyaaaa52
xyzzyaaaa52%avl_parent=>xyzzyaaab52
xyzzyaaab52%avl_left=>item
item%avl_parent=>xyzzyaaab52
call xyzzyaabt1(xyzzyaaaa52)
call xyzzyaabt1(item)
call xyzzyaabt1(xyzzyaaab52)
item=>xyzzyaaab52
elseif(xyzzyaaad52>0)then
if(associated(xyzzyaaaa52%avl_left))then
item%avl_right=>xyzzyaaaa52%avl_left
item%avl_right%avl_parent=>item
else
nullify(item%avl_right)
endif
xyzzyaaaa52%avl_left=>item
if(associated(item%avl_parent))then
xyzzyaaaa52%avl_parent=>item%avl_parent
if(associated(item%avl_parent%avl_right,item))then
item%avl_parent%avl_right=>xyzzyaaaa52
else
item%avl_parent%avl_left=>xyzzyaaaa52
endif
else
nullify(xyzzyaaaa52%avl_parent)
endif
item%avl_parent=>xyzzyaaaa52
call xyzzyaabt1(item)
call xyzzyaabt1(xyzzyaaaa52)
item=>xyzzyaaaa52
endif
endif
end subroutine xyzzyaabs1
subroutine xyzzyaabt1(item)
implicit none
type(casl_item),pointer :: item
integer xyzzyaaaa53,xyzzyaaab53
if(associated(item%avl_left))then
xyzzyaaaa53=item%avl_left%avl_depth
else
xyzzyaaaa53=0
endif
if(associated(item%avl_right))then
xyzzyaaab53=item%avl_right%avl_depth
else
xyzzyaaab53=0
endif
item%avl_depth=max(xyzzyaaaa53,xyzzyaaab53)+1
item%avl_imbalance=xyzzyaaab53-xyzzyaaaa53
end subroutine xyzzyaabt1
subroutine xyzzyaabu1(item)
implicit none
type(casl_item),pointer :: item
type(casl_item),pointer :: xyzzyaaaa54,xyzzyaaab54
nullify(xyzzyaaaa54,xyzzyaaab54)
if(.not.associated(item%avl_left).and..not.associated(item%avl_right))&
&then
if(associated(item%avl_parent))then
if(associated(item%avl_parent%avl_left,item))then
nullify(item%avl_parent%avl_left)
else
nullify(item%avl_parent%avl_right)
endif
xyzzyaaab54=>item%avl_parent
else
nullify(xyzzyaaab1)
endif
elseif(.not.associated(item%avl_left))then
if(associated(item%avl_parent))then
if(associated(item%avl_parent%avl_left,item))then
item%avl_parent%avl_left=>item%avl_right
else
item%avl_parent%avl_right=>item%avl_right
endif
item%avl_right%avl_parent=>item%avl_parent
xyzzyaaab54=>item%avl_right
else
xyzzyaaab1=>item%avl_right
nullify(item%avl_right%avl_parent)
endif
elseif(.not.associated(item%avl_right))then
if(associated(item%avl_parent))then
if(associated(item%avl_parent%avl_left,item))then
item%avl_parent%avl_left=>item%avl_left
else
item%avl_parent%avl_right=>item%avl_left
endif
item%avl_left%avl_parent=>item%avl_parent
xyzzyaaab54=>item%avl_left
else
xyzzyaaab1=>item%avl_left
nullify(item%avl_left%avl_parent)
endif
else
if(item%avl_imbalance>0)then
xyzzyaaaa54=>item%avl_right
do while(associated(xyzzyaaaa54%avl_left))
xyzzyaaaa54=>xyzzyaaaa54%avl_left
enddo
xyzzyaaab54=>xyzzyaaaa54%avl_parent
if(associated(xyzzyaaab54,item))then
nullify(xyzzyaaab54)
else
if(associated(xyzzyaaaa54%avl_right))then
xyzzyaaab54%avl_left=>xyzzyaaaa54%avl_right
xyzzyaaab54%avl_left%avl_parent=>xyzzyaaab54
else
nullify(xyzzyaaab54%avl_left)
endif
endif
if(associated(item%avl_left))then
xyzzyaaaa54%avl_left=>item%avl_left
xyzzyaaaa54%avl_left%avl_parent=>xyzzyaaaa54
else
nullify(xyzzyaaaa54%avl_left)
endif
if(associated(xyzzyaaab54))then
xyzzyaaaa54%avl_right=>item%avl_right
xyzzyaaaa54%avl_right%avl_parent=>xyzzyaaaa54
endif
if(associated(item%avl_parent))then
xyzzyaaaa54%avl_parent=>item%avl_parent
if(associated(item%avl_parent%avl_left,item))then
xyzzyaaaa54%avl_parent%avl_left=>xyzzyaaaa54
else
xyzzyaaaa54%avl_parent%avl_right=>xyzzyaaaa54
endif
else
nullify(xyzzyaaaa54%avl_parent)
xyzzyaaab1=>xyzzyaaaa54
endif
if(.not.associated(xyzzyaaab54))xyzzyaaab54=>xyzzyaaaa54
else
xyzzyaaaa54=>item%avl_left
do while(associated(xyzzyaaaa54%avl_right))
xyzzyaaaa54=>xyzzyaaaa54%avl_right
enddo
xyzzyaaab54=>xyzzyaaaa54%avl_parent
if(associated(xyzzyaaab54,item))then
nullify(xyzzyaaab54)
else
if(associated(xyzzyaaaa54%avl_left))then
xyzzyaaab54%avl_right=>xyzzyaaaa54%avl_left
xyzzyaaab54%avl_right%avl_parent=>xyzzyaaab54
else
nullify(xyzzyaaab54%avl_right)
endif
endif
if(associated(item%avl_right))then
xyzzyaaaa54%avl_right=>item%avl_right
xyzzyaaaa54%avl_right%avl_parent=>xyzzyaaaa54
else
nullify(xyzzyaaaa54%avl_right)
endif
if(associated(xyzzyaaab54))then
xyzzyaaaa54%avl_left=>item%avl_left
xyzzyaaaa54%avl_left%avl_parent=>xyzzyaaaa54
endif
if(associated(item%avl_parent))then
xyzzyaaaa54%avl_parent=>item%avl_parent
if(associated(item%avl_parent%avl_right,item))then
xyzzyaaaa54%avl_parent%avl_right=>xyzzyaaaa54
else
xyzzyaaaa54%avl_parent%avl_left=>xyzzyaaaa54
endif
else
nullify(xyzzyaaaa54%avl_parent)
xyzzyaaab1=>xyzzyaaaa54
endif
if(.not.associated(xyzzyaaab54))xyzzyaaab54=>xyzzyaaaa54
endif
endif
if(associated(xyzzyaaab54))call xyzzyaabr1(xyzzyaaab54)
nullify(item%avl_parent,item%avl_left,item%avl_right)
item%avl_depth=0
item%avl_imbalance=0
end subroutine xyzzyaabu1
subroutine xyzzyaabv1(label,item)
implicit none
character(*),intent(in) :: label
type(casl_item),pointer :: item
type(metachar) mlabel
call xyzzyaabi1(label,mlabel)
if(xyzzyaaca1(mlabel)==0)then
item=>xyzzyaaaa1
return
endif
item=>xyzzyaaab1
do
select case(compare_metachar(mlabel,item%full_unique_name))
case('=')
call xyzzyaabw1(mlabel)
return
case('>')
if(.not.associated(item%avl_right))then
nullify(item)
call xyzzyaabw1(mlabel)
return
endif
item=>item%avl_right
case('<')
if(.not.associated(item%avl_left))then
nullify(item)
call xyzzyaabw1(mlabel)
return
endif
item=>item%avl_left
end select
enddo
end subroutine xyzzyaabv1
subroutine xyzzyaabw1(string)
implicit none
type(metachar),intent(inout) :: string
if(associated(string%chars))then
deallocate(string%chars)
nullify(string%chars)
endif
end subroutine xyzzyaabw1
subroutine xyzzyaabx1(string1,string2)
implicit none
type(metachar),intent(in) :: string1
type(metachar),intent(inout) :: string2
integer xyzzyaaaa57,xyzzyaaab57
call xyzzyaabw1(string2)
xyzzyaaaa57=xyzzyaaca1(string1)
if(xyzzyaaaa57==0)return
allocate(string2%chars(xyzzyaaaa57))
do xyzzyaaab57=1,xyzzyaaaa57
string2%chars(xyzzyaaab57)=string1%chars(xyzzyaaab57)
enddo
end subroutine xyzzyaabx1
subroutine xyzzyaaby1(chars,string)
implicit none
character(*),intent(in) :: chars
type(metachar),intent(inout) :: string
integer xyzzyaaaa58,xyzzyaaab58
call xyzzyaabw1(string)
xyzzyaaaa58=len(chars)
if(xyzzyaaaa58==0)return
allocate(string%chars(xyzzyaaaa58))
do xyzzyaaab58=1,xyzzyaaaa58
string%chars(xyzzyaaab58)=chars(xyzzyaaab58:xyzzyaaab58)
enddo
end subroutine xyzzyaaby1
subroutine xyzzyaabz1(string,chars)
implicit none
type(metachar),intent(in) :: string
character(*),intent(out) :: chars
integer xyzzyaaaa59,xyzzyaaab59,xyzzyaaac59
xyzzyaaab59=len(chars)
if(xyzzyaaab59<1)return
chars=''
if(.not.associated(string%chars))return
xyzzyaaac59=size(string%chars)
do xyzzyaaaa59=1,min(xyzzyaaab59,xyzzyaaac59)
chars(xyzzyaaaa59:xyzzyaaaa59)=string%chars(xyzzyaaaa59)
enddo
end subroutine xyzzyaabz1
integer function xyzzyaaca1(string)
implicit none
type(metachar),intent(in) :: string
xyzzyaaca1=0
if(.not.associated(string%chars))return
xyzzyaaca1=size(string%chars)
end function xyzzyaaca1
integer function xyzzyaacb1(string)
implicit none
type(metachar),intent(in) :: string
integer xyzzyaaaa61,xyzzyaaab61,xyzzyaaac61
xyzzyaacb1=0
xyzzyaaaa61=xyzzyaaca1(string)
if(xyzzyaaaa61==0)return
xyzzyaaab61=1
do xyzzyaaac61=xyzzyaaaa61,1,-1
if(string%chars(xyzzyaaac61)/=' ')exit
enddo
xyzzyaacb1=max(xyzzyaaac61-xyzzyaaab61+1,0)
end function xyzzyaacb1
logical function xyzzyaacc1(string,chars)
implicit none
type(metachar),intent(in) :: string
character(*),intent(in) :: chars
integer xyzzyaaaa62,xyzzyaaab62,xyzzyaaac62
xyzzyaacc1=.false.
xyzzyaaac62=len(chars)
xyzzyaaab62=xyzzyaaca1(string)
if(xyzzyaaac62/=xyzzyaaab62)return
do xyzzyaaaa62=1,xyzzyaaac62
if(string%chars(xyzzyaaaa62)/=chars(xyzzyaaaa62:xyzzyaaaa62))return
enddo
xyzzyaacc1=.true.
end function xyzzyaacc1
character function compare_metachar(string1,string2)
implicit none
type(metachar),intent(in) :: string1,string2
integer xyzzyaaaa63,xyzzyaaab63,xyzzyaaac63
character c1,c2
xyzzyaaaa63=xyzzyaaca1(string1)
xyzzyaaab63=xyzzyaaca1(string2)
if(xyzzyaaaa63==0.or.xyzzyaaab63==0)then
compare_metachar='='
if(xyzzyaaaa63/=0)then
compare_metachar='>'
elseif(xyzzyaaab63/=0)then
compare_metachar='<'
endif
return
endif
do xyzzyaaac63=1,min(xyzzyaaaa63,xyzzyaaab63)
c1=string1%chars(xyzzyaaac63)
c2=string2%chars(xyzzyaaac63)
if(c1/=c2)then
compare_metachar='<'
if(c1>c2)compare_metachar='>'
return
endif
enddo
compare_metachar='='
if(xyzzyaaaa63>xyzzyaaab63)then
compare_metachar='>'
elseif(xyzzyaaaa63<xyzzyaaab63)then
compare_metachar='<'
endif
end function compare_metachar
subroutine xyzzyaacd1(io,string)
implicit none
integer,intent(in) :: io
type(metachar),intent(inout) :: string
integer,parameter :: xyzzyaaaa64=80
character(xyzzyaaaa64) buffer
integer xyzzyaaab64,xyzzyaaac64
call xyzzyaabw1(string)
do
read(io,fmt='(a)',advance='no',eor=9999,size=xyzzyaaab64,iostat=xyzzya&
&aac64)buffer(1:xyzzyaaaa64)
if(xyzzyaaac64/=0)return
call xyzzyaace1(string,buffer)
enddo
9999 continue
call xyzzyaace1(string,buffer(1:xyzzyaaab64))
end subroutine xyzzyaacd1
subroutine xyzzyaace1(string,chars)
implicit none
character(*),intent(in) :: chars
type(metachar),intent(inout) :: string
integer xyzzyaaaa65,xyzzyaaab65,xyzzyaaac65
character,allocatable :: xyzzyaaad65(:)
xyzzyaaaa65=len(chars)
if(xyzzyaaaa65==0)return
xyzzyaaab65=0
if(associated(string%chars))xyzzyaaab65=size(string%chars)
if(xyzzyaaab65>0)then
allocate(xyzzyaaad65(xyzzyaaab65))
xyzzyaaad65(1:xyzzyaaab65)=string%chars(1:xyzzyaaab65)
endif
call xyzzyaabw1(string)
allocate(string%chars(xyzzyaaab65+xyzzyaaaa65))
if(xyzzyaaab65>0)then
string%chars(1:xyzzyaaab65)=xyzzyaaad65(1:xyzzyaaab65)
deallocate(xyzzyaaad65)
endif
do xyzzyaaac65=1,xyzzyaaaa65
string%chars(xyzzyaaab65+xyzzyaaac65)=chars(xyzzyaaac65:xyzzyaaac65)
enddo
end subroutine xyzzyaace1
subroutine xyzzyaacf1(string,string1)
implicit none
type(metachar),intent(in) :: string1
type(metachar),intent(inout) :: string
integer xyzzyaaaa66,xyzzyaaab66
character,allocatable :: xyzzyaaac66(:)
xyzzyaaaa66=xyzzyaaca1(string1)
if(xyzzyaaaa66==0)return
xyzzyaaab66=0
if(associated(string%chars))xyzzyaaab66=size(string%chars)
if(xyzzyaaab66>0)then
allocate(xyzzyaaac66(xyzzyaaab66))
xyzzyaaac66(1:xyzzyaaab66)=string%chars(1:xyzzyaaab66)
endif
call xyzzyaabw1(string)
allocate(string%chars(xyzzyaaab66+xyzzyaaaa66))
if(xyzzyaaab66>0)then
string%chars(1:xyzzyaaab66)=xyzzyaaac66(1:xyzzyaaab66)
deallocate(xyzzyaaac66)
endif
string%chars(xyzzyaaab66+1:xyzzyaaab66+xyzzyaaaa66)=string1%chars(1:xy&
&zzyaaaa66)
end subroutine xyzzyaacf1
subroutine xyzzyaacg1(string,chars)
implicit none
character(*),intent(in) :: chars
type(metachar),intent(inout) :: string
integer xyzzyaaaa67,xyzzyaaab67,xyzzyaaac67
character,allocatable :: xyzzyaaad67(:)
xyzzyaaaa67=len(chars)
if(xyzzyaaaa67==0)return
xyzzyaaab67=0
if(associated(string%chars))xyzzyaaab67=size(string%chars)
if(xyzzyaaab67>0)then
allocate(xyzzyaaad67(xyzzyaaab67))
xyzzyaaad67(1:xyzzyaaab67)=string%chars(1:xyzzyaaab67)
endif
call xyzzyaabw1(string)
allocate(string%chars(xyzzyaaab67+xyzzyaaaa67))
do xyzzyaaac67=1,xyzzyaaaa67
string%chars(xyzzyaaac67)=chars(xyzzyaaac67:xyzzyaaac67)
enddo
if(xyzzyaaab67>0)then
string%chars(xyzzyaaaa67+1:xyzzyaaab67+xyzzyaaaa67)=xyzzyaaad67(1:xyzz&
&yaaab67)
deallocate(xyzzyaaad67)
endif
end subroutine xyzzyaacg1
subroutine xyzzyaach1(string,string1)
implicit none
type(metachar),intent(in) :: string1
type(metachar),intent(inout) :: string
integer xyzzyaaaa68,xyzzyaaab68
character,allocatable :: xyzzyaaac68(:)
xyzzyaaaa68=xyzzyaaca1(string1)
if(xyzzyaaaa68==0)return
xyzzyaaab68=0
if(associated(string%chars))xyzzyaaab68=size(string%chars)
if(xyzzyaaab68>0)then
allocate(xyzzyaaac68(xyzzyaaab68))
xyzzyaaac68(1:xyzzyaaab68)=string%chars(1:xyzzyaaab68)
endif
call xyzzyaabw1(string)
allocate(string%chars(xyzzyaaab68+xyzzyaaaa68))
string%chars(1:xyzzyaaaa68)=string1%chars(1:xyzzyaaaa68)
if(xyzzyaaab68>0)then
string%chars(xyzzyaaaa68+1:xyzzyaaaa68+xyzzyaaab68)=xyzzyaaac68(1:xyzz&
&yaaab68)
deallocate(xyzzyaaac68)
endif
end subroutine xyzzyaach1
function xyzzyaaci1(string,len_chars) result(chars)
implicit none
type(metachar),intent(in) :: string
integer,intent(in) :: len_chars
character(len_chars) chars
integer xyzzyaaaa69,xyzzyaaab69
xyzzyaaaa69=0
chars=''
if(.not.associated(string%chars))return
xyzzyaaaa69=size(string%chars)
do xyzzyaaab69=1,min(len_chars,xyzzyaaaa69)
chars(xyzzyaaab69:xyzzyaaab69)=string%chars(xyzzyaaab69)
enddo
end function xyzzyaaci1
subroutine xyzzyaacj1(string,i0,j1)
implicit none
integer,intent(in) :: i0
integer,intent(in),optional :: j1
type(metachar),intent(inout) :: string
integer xyzzyaaaa70,xyzzyaaab70,xyzzyaaac70
type(metachar) temp_string
xyzzyaaaa70=xyzzyaaca1(string)
if(xyzzyaaaa70==0)return
xyzzyaaac70=xyzzyaaaa70
if(present(j1))xyzzyaaac70=min(j1,xyzzyaaaa70)
if(i0>xyzzyaaac70)then
call xyzzyaabw1(string)
return
endif
xyzzyaaab70=xyzzyaaac70-i0+1
allocate(temp_string%chars(xyzzyaaab70))
temp_string%chars(1:xyzzyaaab70)=string%chars(i0:xyzzyaaac70)
call xyzzyaabw1(string)
string%chars=>temp_string%chars
nullify(temp_string%chars)
end subroutine xyzzyaacj1
subroutine xyzzyaack1(string)
implicit none
type(metachar),intent(inout) :: string
integer xyzzyaaaa71,xyzzyaaab71,xyzzyaaac71
xyzzyaaaa71=xyzzyaaca1(string)
if(xyzzyaaaa71==0)return
do xyzzyaaab71=1,xyzzyaaaa71
if(string%chars(xyzzyaaab71)/=' ')exit
enddo
do xyzzyaaac71=xyzzyaaaa71,1,-1
if(string%chars(xyzzyaaac71)/=' ')exit
enddo
call xyzzyaacj1(string,xyzzyaaab71,xyzzyaaac71)
end subroutine xyzzyaack1
integer function xyzzyaacl1(string,c,back)
implicit none
character,intent(in) :: c
type(metachar),intent(in) :: string
logical,intent(in),optional :: back
integer xyzzyaaaa72
logical xyzzyaaab72
xyzzyaacl1=0
xyzzyaaaa72=xyzzyaaca1(string)
if(xyzzyaaaa72==0)return
xyzzyaaab72=.false.
if(present(back))xyzzyaaab72=back
if(.not.xyzzyaaab72)then
do xyzzyaacl1=1,xyzzyaaaa72
if(string%chars(xyzzyaacl1)==c)return
enddo
else
do xyzzyaacl1=xyzzyaaaa72,1,-1
if(string%chars(xyzzyaacl1)==c)return
enddo
endif
xyzzyaacl1=0
end function xyzzyaacl1
character(20) function i2s(n)
implicit none
integer,intent(in) :: n
integer xyzzyaaaa73,xyzzyaaab73,xyzzyaaac73
i2s=''
xyzzyaaaa73=abs(n)
do while(xyzzyaaaa73>0)
xyzzyaaab73=xyzzyaaaa73/10
xyzzyaaac73=xyzzyaaaa73-xyzzyaaab73*10
if(xyzzyaaac73<0.or.xyzzyaaac73>9)then
i2s='[overflow]'
return
endif
i2s=achar(ichar('0')+xyzzyaaac73)//trim(i2s)
xyzzyaaaa73=xyzzyaaab73
enddo
if(len_trim(i2s)==0)then
i2s='0'
elseif(n<0)then
i2s='-'//trim(i2s)
endif
end function i2s
end module casl

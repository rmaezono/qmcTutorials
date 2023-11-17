PROGRAM combine_cdata
!---------------------------------------------------------------!
! COMBINE_CDATA                                                 !
! =============                                                 !
!                                                               !
! This program reads two or more correlation.data and produces  !
! a correlation.data file on stdout where the parameters are    !
! the average of the parameters in the original files.          !
!                                                               !
! When optimizing a wave function, the first few (often just    !
! the first) iterations move the parameters towards the region  !
! of the minimum, and afterwards the parameters oscillate       !
! around the minimum without any further convergence due to     !
! the numerics of the process. It turns out that the values     !
! each parameter takes in this stage of the optimization are    !
! normally distributed around the optimum value (J.R. Trail, to !
! be published), so it makes sense to simply average them to    !
! attempt to produce the best wave function.                    !
!                                                               !
! This procedure may not be safely applicable when cutoffs or   !
! other non-linear parameters are allowed to vary. Use wisely!  !
!                                                               !
! PLR 08.2009                                                   !
!---------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,PARAMETER :: dp=kind(1.d0)
 CHARACTER(256) line,filename
 CHARACTER(256),POINTER :: ref_line(:)=>null()
 CHARACTER(80),POINTER :: line_fmt(:)=>null()
 CHARACTER(80) f1,f2,in_section,param_fmt,pfmt
 INTEGER o,e,io,nfile,ierr,nparam,nline,delay,np,iline,iparam,ialloc
 LOGICAL exists,has_param
 INTEGER,POINTER :: line_nparam(:)=>null()
 REAL(dp),POINTER :: param(:,:)=>null(),tparam(:)=>null()

 o=6 ; e=0 ; io=10

! Get filename of reference file
 write(e,*)'Enter name of reference correlation.data file'
 read(*,'(a)')filename
 if(len_trim(filename)==0)then
  write(e,*)'Empty string. Quitting.'
  stop
 endif

! Check file exists and open it
 inquire(file=trim(adjustl(filename)),exist=exists)
 if(.not.exists)then
  write(e,*)'File "'//trim(filename)//'" does not exist. Quitting.'
  stop
 endif
 open(io,file=trim(adjustl(filename)),status='old',iostat=ierr)
 if(ierr/=0)then
  write(e,*)'Problem opening file "'//trim(filename)//'". Try again.'
  stop
 endif

! Read file
 nline=0 ; nparam=0 ; in_section='' ; param_fmt='' ; delay=-1
 do
  read(io,'(a)',iostat=ierr)line
  if(ierr/=0)exit
  has_param=.false.
! Check for section change
  f1=field(1,line)
  f2=field(2,line)
  select case(trim(in_section))
  case('')
   if(trim(f1)=='START')then
    select case(trim(f2))
    case('JASTROW','BACKFLOW','MDET','FREE_ORBS','ORBMODS','MOLORBMODS',&
     &'EXMOL')
     in_section=f2
    end select
   endif
  case default
   if(trim(field(1,line))=='END')then
! End of section or end of parameter block
    if(trim(field(2,line))==trim(in_section))in_section=''
    param_fmt='' ; delay=-1
   else
! Check for start of parameter block
    select case(trim(in_section))
    case('JASTROW')
     select case(trim(f1)//' '//trim(f2))
     case('Cutoff (a.u.)','Cutoffs (a.u.)')
      param_fmt='PIx' ; delay=1
     case('Parameter values','Parameter value')
      param_fmt='PIx' ; delay=1
     end select
    case('BACKFLOW')
     select case(trim(f1)//' '//trim(f2))
     case('Cutoff (a.u.)','Cut-off radii')
      param_fmt='PIx' ; delay=1
     case('Parameter ;','Parameter values')
      param_fmt='PIx' ; delay=1
     case('Nucleus ;') ! AE cutoffs
      param_fmt='IIPx' ; delay=1
     end select
    case('FREE_ORBS')
     select case(trim(f1)//' '//trim(f2))
     case('PW coefficients')
      param_fmt='PIx' ; delay=1
     case('Parameter ;')
      param_fmt='PpIx' ; delay=1
     end select
    case('MDET')
     select case(trim(f1)//' '//trim(f2))
     case('Multideterminant/excitation specification')
      param_fmt='PIIx' ; delay=3
     end select
    case('MOLORBMODS','ORBMODS')
     select case(trim(f1)//' '//trim(f2))
     case('Parameter values','Parameters in')
      param_fmt='PIx' ; delay=1
     end select
    case('EXMOL')
     select case(trim(f1)//' '//trim(f2))
     case('Coefficient ;','Parameter ;')
      param_fmt='PIx' ; delay=1
     end select
    end select
   endif
   if(delay>0)then
    delay=delay-1
   elseif(delay==0)then
    has_param=.true.
   endif
  end select
! This line potentially contains parameters. Read them.
  np=0
  if(has_param)then
   nullify(tparam)
   call read_fmt(trim(line),trim(param_fmt),tparam,ierr,actual_fmt=pfmt)
   if(ierr/=0)then
    param_fmt='' ; delay=-1
   elseif(associated(tparam))then
    np=size(tparam,1)
    nparam=nparam+np
    call resize_pointer_dble2((/nparam,1/),param)
    param(nparam-np+1:nparam,1)=tparam(1:np)
    deallocate(tparam)
    nullify(tparam)
   endif
  endif
  if(np==0)pfmt=''
! Store line
  nline=nline+1
  call resize_pointer_char1(256,(/nline/),ref_line)
  call resize_pointer_char1(80,(/nline/),line_fmt)
  call resize_pointer_int1((/nline/),line_nparam)
  ref_line(nline)=trim(line)
  line_nparam(nline)=np
  line_fmt(nline)=pfmt
!(DEBUG)
!if(np==0)then
! write(e,*)trim(line)
!else
! write(e,*)trim(pfmt)
!endif
 enddo

 close(io)

! Check we have read something
 if(nline==0)then
  write(e,*)'File "'//trim(adjustl(filename))//'" is empty.'
  stop
 endif
 if(nparam==0)then
  write(e,*)'File "'//trim(adjustl(filename))//'" contains no parameters.'
  stop
 endif

! Loop over other files
 nfile=1
 do

! Get filename of next file
  write(e,*)'Enter name of correlation.data file #'//trim(i2s(nfile+1))//&
   &', or empty to finish.'
  read(*,'(a)')filename
  if(len_trim(filename)==0)exit

! Check file exists and open it
  inquire(file=trim(adjustl(filename)),exist=exists)
  if(.not.exists)then
   write(e,*)'File "'//trim(filename)//'" does not exist. Skipping.'
   cycle
  endif
  open(io,file=trim(adjustl(filename)),status='old',iostat=ierr)
  if(ierr/=0)then
   write(e,*)'Problem opening file "'//trim(filename)//'". Try again.'
   cycle
  endif

! Loop over lines
  nfile=nfile+1 ; iline=0 ; iparam=0
  call resize_pointer_dble2((/nparam,nfile/),param)
  do
! Read line
   read(io,'(a)',iostat=ierr)line
! Check file is not too short
   if(ierr/=0)then
    if(iline/=nline)then
     write(e,*)'File "'//trim(filename)//'" has fewer lines than reference &
      &file. Ignoring this file.'
     nfile=nfile-1
     call resize_pointer_dble2((/nparam,nfile/),param)
    endif
    exit
   endif
! Check file is not too long
   iline=iline+1
   if(iline>nline)then
    write(e,*)'File "'//trim(filename)//'" has more lines than reference &
     &file. Ignoring this file.'
    nfile=nfile-1
    call resize_pointer_dble2((/nparam,nfile/),param)
    exit
   endif
! Get parameter in this line, if any
   np=line_nparam(iline)
   if(np>0)then
    iparam=iparam+np
    nullify(tparam)
    call read_fmt(trim(line),trim(line_fmt(iline)),tparam,ierr,ref_line=&
     &trim(ref_line(iline)))
    if(ierr/=0)then
     write(e,*)'Error #'//trim(i2s(ierr))//' reading parameters at line #'//&
      &trim(i2s(iline))//' in file "'//trim(filename)//'". Ignoring this file.'
     nfile=nfile-1
     call resize_pointer_dble2((/nparam,nfile/),param)
     exit
    endif
    param(iparam-np+1:iparam,nfile)=tparam(1:np)
    deallocate(tparam)
    nullify(tparam)
   else
! Check full line matches reference
    if(trim(adjustl(line))/=trim(adjustl(ref_line(iline))))then
     write(e,*)'Line #'//trim(i2s(iline))//' in file "'//trim(filename)//&
      &'" does not match reference. Ignoring this file.'
     nfile=nfile-1
     call resize_pointer_dble2((/nparam,nfile/),param)
     exit
    endif
   endif
  enddo ! loop over lines

! Close file
  close(io)

 enddo ! loop over files

! Average parameters
 allocate(tparam(nparam),stat=ialloc)
 if(ialloc/=0)then
  write(e,*)'Allocation error (TPARAM).'
  stop
 endif
 do iparam=1,nparam
  tparam(iparam)=sum(param(iparam,1:nfile))/dble(nfile)
 enddo ! iparam
 deallocate(param)
 nullify(param)

! Produce correlation.data file on stdout
 iparam=0
 do iline=1,nline
  np=line_nparam(iline)
  if(np==0)then
   write(o,'(a)')trim(ref_line(iline))
  else
   iparam=iparam+np
   call write_fmt(trim(ref_line(iline)),trim(line_fmt(iline)),&
    &tparam(iparam-np+1:iparam))
  endif
 enddo ! iline

! Finish
 deallocate(ref_line,line_fmt,line_nparam)
 nullify(ref_line,line_fmt,line_nparam)


CONTAINS


 SUBROUTINE read_fmt(line,line_fmt,param,ierr,actual_fmt,ref_line)
!-------------------------------------------------------------------!
! Read LINE, which is in the format given by LINE_FMT, and extract  !
! the values of however many parameters there are into PARAM(:).    !
! IERR/=0 signals a mismatch between LINE_FMT and LINE. Optionally, !
! return the format encountered in ACTUAL_FMT (differing from       !
! LINE_FMT only in optional fields), and compare non-parameter      !
! fields against a reference line REF_LINE.                         !
!-------------------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: line,line_fmt
 CHARACTER(*),INTENT(in),OPTIONAL :: ref_line
 CHARACTER(80),INTENT(out),OPTIONAL :: actual_fmt
 INTEGER,INTENT(out) :: ierr
 REAL(dp),POINTER :: param(:)
 CHARACTER fm,cfm
 CHARACTER(len_trim(line_fmt)) ofmt
 CHARACTER(len_trim(line)) f1,rf1,char_t1
 INTEGER nfield,ifield,jfield,int_t1,nparam,nshift
 LOGICAL eraseable(len_trim(line_fmt)),erased(len_trim(line_fmt)),&
  &is_param(len_trim(line_fmt))
 REAL(dp) dp_t1

! Initialize
 ierr=0 ; ofmt=''
 nparam=0 ; nfield=len_trim(line_fmt)
 erased=.false. ; is_param=.false.
 do jfield=1,nfield
  fm=line_fmt(jfield:jfield)
  cfm=switch_case(fm)
  if(cfm=='X'.and.jfield<nfield)then ! 'X' can only appear as last item
   ierr=128 ; return
  endif
  eraseable(jfield)=fm/=cfm
 enddo

! Loop over fields in line.
 ifield=1 ; nshift=0
 do
  f1=field(ifield-nshift,line)
! Check if we've got to the end of format string.
  if(ifield>nfield)then
   if(len_trim(f1)==0)exit
   ierr=129
  endif
  if(ierr==0)then
! Get type of this field.
   fm=line_fmt(ifield:ifield)
   cfm=switch_case(fm)
! Check if this is the 'X' field of the format string.
   if(cfm=='X')then
    if(fm=='X'.and.len_trim(f1)==0)ierr=130
    if(ierr==0)then
     if(len_trim(f1)>0)ofmt(nfield-nshift:nfield-nshift)='X'
     if(present(ref_line))then
      jfield=ifield
      do
       f1=field(jfield-nshift,line)
       rf1=field(jfield-nshift,ref_line)
       if(trim(f1)/=trim(rf1))ierr=131
       if(len_trim(f1)==0.or.len_trim(rf1)==0)exit
       jfield=jfield+1
      enddo
     endif
     if(ierr==0)exit
    endif ! ierr==0
   endif ! cfm=='X'
   if(ierr==0)then
! Try to read whatever we need to read.
    f1=field(ifield-nshift,line)
    select case(cfm)
    case('P') ; read(f1,*,iostat=ierr)dp_t1
    case('I') ; read(f1,*,iostat=ierr)int_t1
    case('A') ; read(f1,*,iostat=ierr)char_t1
    end select
    if(ierr==0)then
! Proceed as required for item type.
     select case(cfm)
     case('P')
      nparam=nparam+1
      call resize_pointer_dble1((/nparam/),param)
      param(nparam)=dp_t1
      is_param(ifield)=.true.
     case('I','A')
      if(present(ref_line))then
       rf1=field(ifield-nshift,ref_line)
       if(trim(f1)/=trim(rf1))ierr=132
      endif
     end select ! cfm
    endif ! ierr==0
   endif ! ierr==0
  endif ! ierr==0
! Deal with errors.
  if(ierr/=0)then
   if(any(eraseable(1:ifield)))then
! Set IFIELD to last optional field.
    do jfield=ifield,1,-1
     if(eraseable(jfield))then
      ifield=jfield
      exit
     endif
    enddo ! jfield
! Erase field IFIELD and recalculate the field shift.
    erased(ifield)=.true.
    nshift=count(erased(1:ifield))
    ofmt(ifield-nshift:)=''
! Reset eraseable flags for fields after IFIELD.
    eraseable(ifield)=.false.
    do jfield=ifield+1,nfield
     fm=line_fmt(jfield:jfield)
     cfm=switch_case(fm)
     eraseable(jfield)=fm/=cfm
     erased(jfield)=.false.
    enddo ! jfield
! Delete params found on or after IFIELD.
    nparam=count(is_param(1:ifield-1))
    call resize_pointer_dble1((/nparam/),param)
    ierr=0
   else
! Unrecoverable error (since there are no eraseable fields to resort to).
    return
   endif
  else
! No errors. Append to format and increase field counter.
   ofmt(ifield-nshift:ifield-nshift)=cfm
   ifield=ifield+1
  endif
 enddo ! fields in line

 if(present(actual_fmt))actual_fmt=ofmt

 END SUBROUTINE read_fmt


 SUBROUTINE write_fmt(line,line_fmt,param)
!---------------------------------------------------------!
! Print LINE, which follows the format given in LINE_FMT, !
! replacing the parameters with the elements of PARAM.    !
!---------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: line,line_fmt
 REAL(dp),INTENT(in) :: param(:)
 CHARACTER(256) line_out
 CHARACTER(80) tmpr
 INTEGER ipos0,ipos_start,ipos_end,ifield,iparam
 line_out='' ; ifield=1 ; ipos0=0 ; iparam=0
 do
  if(ipos0>=len_trim(line))exit
  ipos_start=ipos0+verify(line(ipos0+1:),' ')
  if(ipos_start<=ipos0)exit
  ipos_end=ipos_start-1+scan(line(ipos_start:),' ')-1
  if(ipos_end<ipos_start)ipos_end=len_trim(line)
  if(line_fmt(ifield:ifield)=='P')then
   iparam=iparam+1
   write(tmpr,*)param(iparam)
   line_out=trim(line_out)//line(ipos0+1:ipos_start-1)//trim(adjustl(tmpr))
  else
   line_out=trim(line_out)//line(ipos0+1:ipos_end)
  endif
  ipos0=ipos_end
  ifield=ifield+1
 enddo
 write(o,'(a)')trim(line_out)
 END SUBROUTINE write_fmt


 SUBROUTINE resize_pointer_dble1(dims,pt)
!------------------------------------------------------!
! Allocate or resize a first-rank dble pointer PT to   !
! size DIMS, keeping any exisiting data untouched. New !
! elements initialized to zero.                        !
!------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: dims(1)
 REAL(dp),POINTER :: pt(:),pt_new(:)
 INTEGER old_dims(1),ialloc
 if(.not.associated(pt))then
  if(all(dims==0))return
  allocate(pt(dims(1)),stat=ialloc)
  if(ialloc/=0)stop
  pt=0.d0 ; return
 endif
 old_dims=shape(pt)
 if(all(old_dims==dims))return
 if(all(dims==0))then
  deallocate(pt) ; nullify(pt) ; return
 endif
 allocate(pt_new(dims(1)),stat=ialloc)
 if(any(old_dims<dims))pt_new=0
 pt_new(1:min(old_dims(1),dims(1)))=pt(1:min(old_dims(1),dims(1)))
 deallocate(pt)
 pt=>pt_new
 nullify(pt_new)
 END SUBROUTINE resize_pointer_dble1


 SUBROUTINE resize_pointer_dble2(dims,pt)
!------------------------------------------------------!
! Allocate or resize a second-rank dble pointer PT to  !
! size DIMS, keeping any exisiting data untouched. New !
! elements initialized to zero.                        !
!------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: dims(2)
 REAL(dp),POINTER :: pt(:,:),pt_new(:,:)
 INTEGER old_dims(2),ialloc
 if(.not.associated(pt))then
  if(all(dims==0))return
  allocate(pt(dims(1),dims(2)),stat=ialloc)
  if(ialloc/=0)stop
  pt=0.d0 ; return
 endif
 old_dims=shape(pt)
 if(all(old_dims==dims))return
 if(all(dims==0))then
  deallocate(pt) ; nullify(pt) ; return
 endif
 allocate(pt_new(dims(1),dims(2)),stat=ialloc)
 if(any(old_dims<dims))pt_new=0
 pt_new(1:min(old_dims(1),dims(1)),1:min(old_dims(2),dims(2)))=&
  &pt(1:min(old_dims(1),dims(1)),1:min(old_dims(2),dims(2)))
 deallocate(pt)
 pt=>pt_new
 nullify(pt_new)
 END SUBROUTINE resize_pointer_dble2


 SUBROUTINE resize_pointer_int1(dims,pt)
!-------------------------------------------------------!
! Allocate or resize a first-rank integer pointer PT to !
! size DIMS, keeping any exisiting data untouched. New  !
! elements initialized to zero.                         !
!-------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: dims(1)
 INTEGER,POINTER :: pt(:),pt_new(:)
 INTEGER old_dims(1),ialloc
 if(.not.associated(pt))then
  if(all(dims==0))return
  allocate(pt(dims(1)),stat=ialloc)
  if(ialloc/=0)stop
  pt=0 ; return
 endif
 old_dims=shape(pt)
 if(all(old_dims==dims))return
 if(all(dims==0))then
  deallocate(pt) ; nullify(pt) ; return
 endif
 allocate(pt_new(dims(1)),stat=ialloc)
 if(any(old_dims<dims))pt_new=0
 pt_new(1:min(old_dims(1),dims(1)))=pt(1:min(old_dims(1),dims(1)))
 deallocate(pt)
 pt=>pt_new
 nullify(pt_new)
 END SUBROUTINE resize_pointer_int1


 SUBROUTINE resize_pointer_char1(sz,dims,pt)
!------------------------------------------------------!
! Allocate or resize a first-rank character pointer PT !
! to size DIMS, keeping any exisiting data untouched.  !
! New elements initialized to empty string.            !
!------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: sz,dims(1)
 CHARACTER(sz),POINTER :: pt(:),pt_new(:)
 INTEGER old_dims(1),ialloc
 if(.not.associated(pt))then
  if(all(dims==0))return
  allocate(pt(dims(1)),stat=ialloc)
  if(ialloc/=0)stop
  pt='' ; return
 endif
 old_dims=shape(pt)
 if(all(old_dims==dims))return
 if(all(dims==0))then
  deallocate(pt) ; nullify(pt) ; return
 endif
 allocate(pt_new(dims(1)),stat=ialloc)
 if(any(old_dims<dims))pt_new=''
 pt_new(1:min(old_dims(1),dims(1)))=pt(1:min(old_dims(1),dims(1)))
 deallocate(pt)
 pt=>pt_new
 nullify(pt_new)
 END SUBROUTINE resize_pointer_char1


 CHARACTER(20) FUNCTION i2s(n)
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
  i2s(j:j)=achar(ichar('0')+mod(i,10)) ; i=i/10
 enddo
 i=1 ; j=len_trim(i2s)
 do
  if(i>=j)exit
  tmp=i2s(j:j)
  i2s(j:j)=i2s(i:i)
  i2s(i:i)=tmp
  i=i+1 ; j=j-1
 enddo
 i2s=trim(sign)//i2s
 END FUNCTION i2s


 CHARACTER(80) FUNCTION field(n,line)
!--------------------------------------------------------!
! Return the N-th field of string LINE, where the fields !
! are separated by one or more spaces.                   !
! If N is negative, return the |N|-th field of LINE from !
! the end.                                               !
!--------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n
 CHARACTER(*),INTENT(in) :: line
 CHARACTER(len(line)) tline
 INTEGER i,k,absn
 LOGICAL back
 if(n==0)then
  field='' ; return
 endif
 absn=abs(n) ; tline=trim(adjustl(line)) ; back=(n<0)
 do i=1,absn-1
  k=scan(trim(adjustl(tline)),' ',back)
  if(k==0)then
   field='' ; return
  endif
  if(back)then
   tline=trim(adjustl(tline(1:k-1)))
  else
   tline=trim(adjustl(tline(k+1:)))
  endif
 enddo
 k=scan(trim(adjustl(tline)),' ',back)
 if(k==0)then
  field=trim(adjustl(tline))
 elseif(back)then
  field=trim(adjustl(tline(k+1:)))
 else
  field=trim(adjustl(tline(1:k-1)))
 endif
 END FUNCTION field


 FUNCTION switch_case(string,to_lower) RESULT(string_out)
!----------------------------------------------------------!
! Switch case (lower->upper by default, unless TO_LOWER is !
! specified and T) of variable length string STRING.       !
!----------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: string
 LOGICAL,INTENT(in),OPTIONAL :: to_lower
 CHARACTER(len(string)) :: string_out
 INTEGER i,ich,imode
 INTEGER,PARAMETER :: range_start(2)=(/ichar('a'),ichar('A')/),&
  &range_end(2)=(/ichar('z'),ichar('Z')/),                     &
  &range_shift(2)=(/ichar('A')-ichar('a'),ichar('a')-ichar('A')/)
 string_out=string
 imode=1
 if(present(to_lower))then
  if(to_lower)imode=2
 endif
 do i=1,len(string_out)
  ich=ichar(string_out(i:i))
  if(ich>=range_start(imode).and.ich<=range_end(imode))then
   ich=ich+range_shift(imode)
   string_out(i:i)=char(ich)
  endif
 enddo
 END FUNCTION switch_case


END PROGRAM combine_cdata

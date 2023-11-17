! Four F90 units designed to give some of the functionality of AWK
! in finding a pattern and returning the line on which it occurred,
! broken into fields.

 MODULE awk_like
 IMPLICIT NONE
 INTEGER,PARAMETER :: max_line_length=200
 INTEGER,PARAMETER :: NF_max=50 ! Max. no. fields that may be pulled from string
 INTEGER            :: NF ! No. of fields pulled from string
 CHARACTER(LEN=max_line_length) :: fields(NF_max)
 END MODULE awk_like


 SUBROUTINE capitalise(instring,ilen)
!------------------------------------------!
! Convert the string instring to uppercase !
!------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: ilen
 CHARACTER(LEN=ilen),INTENT(inout) :: instring
 INTEGER i, n
 do i=1,ilen,1
  n=iachar(instring(i:i))
  if((n>=97).and.(n<=122))instring(i:i)=achar(n-32)
 enddo
 END SUBROUTINE capitalise


 SUBROUTINE getline(io_unit,keystring,fields,NF,NF_max,ifail,case_sense)
!---------------------------------------------------------------------!
! Scan forward through the file identified by io_unit to find the     !
! string stored in keystring.  Return the fields (items separated by  !
! spaces) from that line.  If keystring contains a single blank space !
! then just take the next line from the file and split it into fields.!
!---------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: io_unit  ! Unit no. of file to read
  INTEGER,INTENT(in) :: NF_max   ! Max no. of fields in a string
  INTEGER,INTENT(out) :: NF      ! Actual no. of fields in a string
  CHARACTER(*),INTENT(inout) :: keystring ! The string to search for
  CHARACTER(*),INTENT(inout) :: fields(NF_max) ! Holds fields pulled from string
  INTEGER,INTENT(out) :: ifail ! set -ve if fail to find string +ve otherwise
  LOGICAL,INTENT(in)  :: case_sense ! Whether search is case-sensitive
  CHARACTER(130) instring
  INTEGER icheck   ! IO checking
  LOGICAL search

  if(keystring==' ')then
! Just get the next line from the file and parse it
   read(io_unit,fmt="(a130)",iostat=icheck)instring
   if(icheck<0)then
    write(6,*) 'Warning - next line is EOF.'
    ifail=-1
   else
    call scan_string(instring,fields,NF,NF_max," ")
    ifail=+1
   endif

  else

   if(case_sense)then
    ! Case-sensitive search...

    ! Search through each line of the file until the search string
    ! is found
    search=.true.
    do while(search)
     read(io_unit,fmt="(a130)",iostat=icheck) instring
     if(icheck<0)then
      write(6,fmt="('Warning - could not find string: ',A)")keystring
      ifail=-1
      exit
     endif

    ! Look to see whether any of the fields contain the string we are
    ! searching for
     if(index(string=instring,substring=keystring)/=0)then
      search=.false.
      ! It does - split this line into fields and return them
      call scan_string(instring,fields,NF,NF_max," ")
      ifail=1
      exit
     endif
    enddo

   else
    ! Search is not case sensitive
    call capitalise(keystring,len(keystring))

    search=.true.
    do while(search)
     read(io_unit,fmt="(A130)",iostat=icheck) instring
     if(icheck<0)then
      write(6,fmt="('Warning - could not find string: ',A)") keystring
      ifail=-1
      exit
     endif

     ! Capitalise each line from file and then try and match with
     ! capitalised search string...
     call capitalise(instring,LEN(instring))

     ! Look to see whether any of the fields contain the string we are
     ! searching for
     if(index(string=instring,substring=keystring)/=0)then
      search=.false.
      ! It does - split this line into fields and return them
      call scan_string(instring,fields,NF,NF_max," ")
      ifail=1
      exit
     endif
    enddo

   endif
  endif

 END SUBROUTINE getline

 SUBROUTINE checkval(io_unit, keystring, istrin, lfound)
!------------------------------------------------------------------------!
! Scan forward through the file identified by io_unit to find if there exist !
! string stored in keystring.  Return .true. in lfound if the string was found, .false. otherwise. !
!------------------------------------------------------------------------!
  USE paramfile
  IMPLICIT none
  INTEGER,INTENT(in) :: io_unit   ! Unit no. of file to read
  CHARACTER(*),INTENT(in) :: keystring ! The string to search for
  INTEGER,INTENT(in) :: istrin ! Which field should contain that string
  INTEGER,PARAMETER :: NF_max=20 ! Max no. of fields in a string
  CHARACTER(130) instring
  CHARACTER(130) fields(NF_max) ! Holds the fields pulled from instring
  INTEGER NF
  INTEGER icheck ! IO checking
  LOGICAL search
  LOGICAL lfound
  rewind(io_unit) !add
  search=.true.
  lfound=.false.
  do while(search)
    read(io_unit,fmt="(A65)",iostat=icheck) instring
    if(icheck<0)then
     search=.false.
    endif
    call scan_string(instring,fields,NF,NF_max," ")
    if(fields(istrin)==keystring)then
     search=.false.
     lfound=.true.
    endif
  enddo
 END SUBROUTINE checkval
 
 SUBROUTINE getval(io_unit,keystring,istrin,ival,value)
!------------------------------------------------------------------------!
! Scan forward through the file identified by io_unit to find the        !
! string stored in keystring.  Return the iposn'th element on this line. !
!------------------------------------------------------------------------!
  USE paramfile
  IMPLICIT none
  INTEGER,INTENT(in) :: io_unit   ! Unit no. of file to read
  CHARACTER(*),INTENT(in) :: keystring ! The string to search for
  INTEGER,INTENT(in) :: istrin ! Which field should contain that string
  INTEGER,INTENT(in) :: ival ! Pos'n of quantity of interest on the line
                             ! on which that string occurs
  REAL(dp),INTENT(out) :: value ! The value of that quantity
  INTEGER,PARAMETER :: NF_max=20 ! Max no. of fields in a string
  CHARACTER(130) instring
  CHARACTER(130) fields(NF_max) ! Holds the fields pulled from instring
  INTEGER NF
  INTEGER icheck ! IO checking
  LOGICAL search
  rewind(io_unit) !add
  search=.true.
  do while(search)
    read(io_unit,fmt="(A65)",iostat=icheck) instring
    if(icheck<0)then
     write(6,*) 'Error - could not find string ',keystring
     stop
    endif
    call scan_string(instring,fields,NF,NF_max," ")
    if(fields(istrin)==keystring)search=.false.
  enddo
  read(fields(ival),*)value
 END SUBROUTINE getval

 SUBROUTINE getval2(io_unit,keystring1,istrin1,keystring2,istrin2,ival,value)
!------------------------------------------------------------------------!
! Scan forward through the file identified by io_unit to find the        !
! string stored in keystring.  Return the iposn'th element on this line. !
!------------------------------------------------------------------------!
  USE paramfile
  IMPLICIT none
  INTEGER,INTENT(in) :: io_unit   ! Unit no. of file to read
  CHARACTER(*),INTENT(in) :: keystring1, keystring2 ! The strings to search for
  INTEGER,INTENT(in) :: istrin1, istrin2 ! Which fields should contain that string
  INTEGER,INTENT(in) :: ival ! Pos'n of quantity of interest on the line
                             ! on which that string occurs
  REAL(dp),INTENT(out) :: value ! The value of that quantity
  INTEGER,PARAMETER :: NF_max=20 ! Max no. of fields in a string
  CHARACTER(130) instring
  CHARACTER(130) fields(NF_max) ! Holds the fields pulled from instring
  INTEGER NF
  INTEGER icheck ! IO checking
  LOGICAL search
  rewind(io_unit) !add
  search=.true.
  do while(search)
    read(io_unit,fmt="(A65)",iostat=icheck) instring
    if(icheck<0)then
     write(6,*) 'Error - could not find string ',keystring1, "or", keystring2
     stop
    endif
    call scan_string(instring,fields,NF,NF_max," ")
    if((fields(istrin1)==keystring1).and.(fields(istrin2)==keystring2))search=.false.
  enddo
  read(fields(ival),*)value
 END SUBROUTINE getval2


 SUBROUTINE scan_string(instring,fields,NF,NF_max,separator)
  IMPLICIT none
  INTEGER,INTENT(out) :: NF
  INTEGER,INTENT(in) :: NF_max
  CHARACTER(*),INTENT(in) :: instring
  CHARACTER(*),INTENT(inout) :: fields(NF_max)
  CHARACTER(1),INTENT(in) :: separator ! Character that separates
                                       ! fields
! Parse the string 'instring' and store the fields (objects separated
! by spaces) in the fields array.

  INTEGER ilength ! Length of string being scanned (excluding
                  ! trailing spaces)
  INTEGER icount  ! Counter for no. of fields found
  INTEGER iposn   ! Current character of string being scanned
  INTEGER ifpt    ! Current character of field being assigned
  LOGICAL newfield! Flag to say whether we are currently within a
                  ! field or between fields (i.e. the last character
                  ! scanned was a space)
  ! String into which to copy instring so as to avoid altering instring
  ! in the following manipulations
  CHARACTER(len=(len(instring))) tempstring

  tempstring=trim(adjustl(instring))
  ilength=len_trim(tempstring)

  iposn=0
  icount=0
  newfield=.true.
  ! Initialise fields as the elements are constructed character by character
  fields=" "

  do while(iposn<ilength)
   iposn=iposn+1
   if(tempstring(iposn:iposn)/=separator)then

    if(newfield)then
     newfield=.false.
     icount=icount + 1
     if(icount>NF_max)then
      write(6,FMT="('Only ',I2,' fields per line allowed.  Increase &
       &NF_max in routine calling scan_string.f90.')") NF_max
      stop
     endif
     fields(icount)(1:1)=tempstring(iposn:iposn)
     ifpt=1
    else
     ifpt=ifpt + 1
     fields(icount)(ifpt:ifpt)=tempstring(iposn:iposn)
    endif

   else
    if(.not.newfield)newfield=.true.
   endif

  enddo

! Return the no. of fields found
  NF=icount

 END SUBROUTINE scan_string

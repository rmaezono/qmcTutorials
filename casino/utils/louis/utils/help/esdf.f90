MODULE esdf
!---------------------------------------------------------------------------!
! ESDF MODULE                                                               !
! Read  the file 'input'.                                                   !
!---------------------------------------------------------------------------!

!---------------------------------------------------------------------------!
!                                                                           !
!  ESDF Summary                                                             !
!  ------------                                                             !
!                                                                           !
!  This module is designed to simplify and enhance the input of data.       !
!  It works from a highly flexible input file.                              !
!                                                                           !
!  DATA is input in a "label <value>" fashion, and input is independent     !
!  of the ordering of the input file. An important feature is the           !
!  requirement that most inputs require default settings to be supplied     !
!  within the main program calling  ESDF. This means that rarely used       !
!  variables will not clutter everyday input files, and, even more          !
!  usefully, "intelligence" may be built into  the main code as the         !
!  defaults may be dependent on other set variables. Block data may also    !
!  be read in. Another important feature is the ability to define           !
!  "physical" values. This means that the input files need not depend on    !
!  the internal physical units used by the main program.                    !
!                                                                           !
!  Usage                                                                    !
!  -----                                                                    !
!                                                                           !
!  First, "USE esdf" wherever you wish to make use of its features. In      !
!  the main program call the initialization routine: call                   !
!  esdf_init('input.esdf'). "input.esdf" is the name of the input file -    !
!  it could be anything. This routine opens the input file, and reads       !
!  into a dynamically allocated storage array. The comments and blank       !
!  lines are stripped out. You are now ready to use the                     !
!  esdf_functions. For example, if you want to read in the number of        !
!  atoms in your calculation, you would use: natom =                        !
!  esdf_integer('NumberOfAtoms',1), where 'NumberOfAtoms' is the label to   !
!  search for in the input file, and '1' is the default. call esdf_close to !
!  deallocate the data arrays. You may then open another input file using   !
!  esdf_init. It is not currently possible to open more that on input       !
!  file at one time.                                                        !
!                                                                           !
!  Syntax                                                                   !
!  ------                                                                   !
!                                                                           !
!  The input file can contain comments. These are defined as anything to    !
!  the right of, and including, '#', ';', or '!'. It is straightforward     !
!  to modify the routine to accept further characters. Blank lines are      !
!  ignored -- use comments and blank lines to make you input file           !
!  readable.                                                                !
!                                                                           !
!  The "labels" are case insensitive (e.g. unitCell is equivalent to        !
!  UnItceLL) and punctuation insensitive (unit.cell is equivalent to        !
!  unit_cell is equivalent to unitcell). Punctuation characters are '.',    !
!  '_', and '-' at the moment. Again - use this feature to improve          !
!  readability.                                                             !
!                                                                           !
!  The following are equivalent ways of defining a physical quantity:       !
!                                                                           !
!  "AgeOfUniverse = 24.d0 s" or "AgeOfUniverse : 24.d0 S" or                !
!  "AgeOfUniverse 24.d0 S"                                                  !
!                                                                           !
!  It would be read in by the main program in the following way:            !
!                                                                           !
!  aou = esdf_physical('ageofuniverse',77.d0,ns)                            !
!                                                                           !
!  "aou" is the double precision variable, 77.d0 is the default number of   !
!  "ns" or nanoseconds. 24s will be converted automatically to its          !
!  equivalent number of nanoseconds.                                        !
!                                                                           !
!  Block data should be placed in the input file as follows:                !
!                                                                           !
!  %block cellvectors                                                       !
!  1.0 1.0 0.0                                                              !
!  0.0 1.0 1.0                                                              !
!  1.0 0.0 1.0                                                              !
!  %endblock cellvectors                                                    !
!                                                                           !
!  And it may be read:                                                      !
!                                                                           !
!    if(esdf_block('CellVectors',nlines))                                   !
!      if(nlines/=3) then (... break out here if the incorrect # of lines)  !
!      do i=1,nlines                                                        !
!       read(block_data(i),*) x,y,z                                         !
!      enddo                                                                !
!    endif                                                                  !
!                                                                           !
!  List of functions                                                        !
!  -----------------                                                        !
!                                                                           !
!  Self explanatory:                                                        !
!                                                                           !
!  esdf_string(label,default)                                               !
!  esdf_integer(label,default)                                              !
!  esdf_single(label,default)                                               !
!  esdf_double(label,default)                                               !
!  esdf_physical(label,default,unit)                                        !
!                                                                           !
!  A little more explanation:                                               !
!                                                                           !
!  esdf_defined(label) is true if "label" found, false otherwise            !
!                                                                           !
!  esdf_boolean(label,default) true if "label yes/true/t (case/punct insens)!
!                              false if"label no/false/f (case/punct insens)!
!                                                                           !
!  The Help feature                                                         !
!  ----------------                                                         !
!                                                                           !
!  The routine "esdf_help(helpword,searchword)" can be used to access the   !
!  information contained within the "esdf_key_mod" module.                  !
!                                                                           !
!  If "helpword" is "search" (case insensitive), then all variables whose   !
!  description contains "searchword" will be output.                        !
!                                                                           !
!  If "helpword" is "basic", "inter", "expert" or "dummy" the variables of  !
!  that type will be displayed.                                             !
!                                                                           !
!  If "helpword" is one of the valid labels, then a description of this     !
!  label will be output.                                                    !
!                                                                           !
!  Finishing off                                                            !
!  -------------                                                            !
!                                                                           !
!  Three routines, "esdf_dump", "esdf_warnout" and "esdf_close", can be     !
!  used to finish the use of ESDF. "esdf_dump" outputs a file               !
!  which could be used as an input file for further runs. "esdf_warnout"    !
!  outputs ESDF warnings to screen, and "esdf_close" deallocates the        !
!  allocated ESDF arrays.                                                   !
!                                                                           !
!---------------------------------------------------------------------------!
 USE dsp
 USE esdf_key
 USE run_control,ONLY : errstop,errstop_quiet,all_stop
 USE store,      ONLY : o

 IMPLICIT NONE

! Set the length of the lines.

 INTEGER,PUBLIC,PARAMETER :: llength=180
 INTEGER,PRIVATE,PARAMETER :: nphys=57,ndump=2000
 INTEGER,PRIVATE :: nrecords,nwarns,ndmp
 LOGICAL keywords_not_in_list
 CHARACTER(llength),PRIVATE,DIMENSION(:),ALLOCATABLE :: llist,warns,dump
 CHARACTER(llength),PRIVATE,DIMENSION(:,:),ALLOCATABLE :: tlist

! The public block data array

 CHARACTER(llength),PUBLIC,DIMENSION(:),ALLOCATABLE :: block_data

! Set the physical units database. Note this should ultimately be merged
! with the one in constants.f90, but it'll do for now to keep esdf separate.

 TYPE phys_unit
  CHARACTER(10) d,n ! d - dimension n - name
  REAL(dp) u ! u - unit
 END TYPE phys_unit

 TYPE(phys_unit),PRIVATE,DIMENSION(nphys) :: phy

!
!     We allow case variations in the units. This could be dangerous
!     (meV --> MeV !!) in real life, but not in this restricted
!     field.
!
! m - mass l - length t - time e - energy f - force p - pressure c- charge
! d - dipole mom - mom inert ef - efield
!

 DATA phy(1)%d /'m'/;DATA phy(1)%n /'kg'/;DATA phy(1)%u /1.d0/
 DATA phy(2)%d /'m'/;DATA phy(2)%n /'g'/;DATA phy(2)%u /1.d-3/
 DATA phy(3)%d /'m'/;DATA phy(3)%n /'amu'/;DATA phy(3)%u /1.66054d-27/
 DATA phy(4)%d /'l'/;DATA phy(4)%n /'m'/;DATA phy(4)%u /1.d0/
 DATA phy(5)%d /'l'/;DATA phy(5)%n /'nm'/;DATA phy(5)%u /1.d-9/
 DATA phy(6)%d /'l'/;DATA phy(6)%n /'ang'/;DATA phy(6)%u /1.d-10/
 DATA phy(7)%d /'l'/;DATA phy(7)%n /'bohr'/;DATA phy(7)%u /0.52917715d-10/
 DATA phy(8)%d /'t'/;DATA phy(8)%n /'s'/;DATA phy(8)%u /1.d0/
 DATA phy(9)%d /'t'/;DATA phy(9)%n /'ns'/;DATA phy(9)%u /1.d-9/
 DATA phy(10)%d /'t'/;DATA phy(10)%n /'ps'/;DATA phy(10)%u /1.d-12/
 DATA phy(11)%d /'t'/;DATA phy(11)%n /'fs'/;DATA phy(11)%u /1.d-15/
 DATA phy(12)%d /'e'/;DATA phy(12)%n /'j'/;DATA phy(12)%u /1.d0/
 DATA phy(13)%d /'e'/;DATA phy(13)%n /'erg'/;DATA phy(13)%u /1.d-7/
 DATA phy(14)%d /'e'/;DATA phy(14)%n /'ev'/;DATA phy(14)%u /1.60219d-19/
 DATA phy(15)%d /'e'/;DATA phy(15)%n /'mev'/;DATA phy(15)%u /1.60219d-22/
 DATA phy(16)%d /'e'/;DATA phy(16)%n /'ry'/;DATA phy(16)%u /2.17991d-18/
 DATA phy(17)%d /'e'/;DATA phy(17)%n /'mry'/;DATA phy(17)%u /2.17991d-21/
 DATA phy(18)%d /'e'/;DATA phy(18)%n /'hartree'/;DATA phy(18)%u /4.35982d-18/
 DATA phy(19)%d /'e'/;DATA phy(19)%n /'kcal/mol'/;DATA phy(19)%u /6.94780d-21/
 DATA phy(20)%d /'e'/;DATA phy(20)%n /'mhartree'/;DATA phy(20)%u /4.35982d-21/
 DATA phy(21)%d /'e'/;DATA phy(21)%n /'kj/mol'/;DATA phy(21)%u /1.6606d-21/
 DATA phy(22)%d /'e'/;DATA phy(22)%n /'hz'/;DATA phy(22)%u /6.6262d-34/
 DATA phy(23)%d /'e'/;DATA phy(23)%n /'thz'/;DATA phy(23)%u /6.6262d-22/
 DATA phy(24)%d /'e'/;DATA phy(24)%n /'cm-1'/;DATA phy(24)%u /1.986d-23/
 DATA phy(25)%d /'e'/;DATA phy(25)%n /'cm^-1'/;DATA phy(25)%u /1.986d-23/
 DATA phy(26)%d /'e'/;DATA phy(26)%n /'cm**-1'/;DATA phy(26)%u /1.986d-23/
 DATA phy(27)%d /'f'/;DATA phy(27)%n /'N'/;DATA phy(27)%u /1.d0/
 DATA phy(28)%d /'f'/;DATA phy(28)%n /'ev/ang'/;DATA phy(28)%u /1.60219d-9/
 DATA phy(29)%d /'f'/;DATA phy(29)%n /'ry/bohr'/;DATA phy(29)%u /4.11943d-8/
 DATA phy(30)%d /'l'/;DATA phy(30)%n /'cm'/;DATA phy(30)%u /1.d-2/
 DATA phy(31)%d /'p'/;DATA phy(31)%n /'pa'/;DATA phy(31)%u /1.d0/
 DATA phy(32)%d /'p'/;DATA phy(32)%n /'mpa'/;DATA phy(32)%u /1.d6/
 DATA phy(33)%d /'p'/;DATA phy(33)%n /'gpa'/;DATA phy(33)%u /1.d9/
 DATA phy(34)%d /'p'/;DATA phy(34)%n /'atm'/;DATA phy(34)%u /1.01325d5/
 DATA phy(35)%d /'p'/;DATA phy(35)%n /'bar'/;DATA phy(35)%u /1.d5/
 DATA phy(36)%d /'p'/;DATA phy(36)%n /'mbar'/;DATA phy(36)%u /1.d11/
 DATA phy(37)%d /'p'/;DATA phy(37)%n /'ry/bohr**3'/;DATA phy(37)%u /1.47108d13/
 DATA phy(38)%d /'p'/;DATA phy(38)%n /'ev/ang**3'/;DATA phy(38)%u /1.60219d11/
 DATA phy(39)%d /'c'/;DATA phy(39)%n /'c'/;DATA phy(39)%u /1.d0/
 DATA phy(40)%d /'c'/;DATA phy(40)%n /'e'/;DATA phy(40)%u /1.602177d-19/
 DATA phy(41)%d /'d'/;DATA phy(41)%n /'C*m'/;DATA phy(41)%u /1.d0/
 DATA phy(42)%d /'d'/;DATA phy(42)%n /'D'/;DATA phy(42)%u /3.33564d-30/
 DATA phy(43)%d /'d'/;DATA phy(43)%n /'debye'/;DATA phy(43)%u /3.33564d-30/
 DATA phy(44)%d /'d'/;DATA phy(44)%n /'e*bohr'/;DATA phy(44)%u /8.47835d-30/
 DATA phy(45)%d /'d'/;DATA phy(45)%n /'e*ang'/;DATA phy(45)%u /1.602177d-29/
 DATA phy(46)%d /'mom'/;DATA phy(46)%n /'kg*m**2'/;DATA phy(46)%u /1.d0/
 DATA phy(47)%d /'mom'/;DATA phy(47)%n /'ry*fs**2'/;DATA phy(47)%u /2.1799d-48/
 DATA phy(48)%d /'ef'/;DATA phy(48)%n /'v/m'/;DATA phy(48)%u /1.d0/
 DATA phy(49)%d /'ef'/;DATA phy(49)%n /'v/nm'/;DATA phy(49)%u /1.d9/
 DATA phy(50)%d /'ef'/;DATA phy(50)%n /'v/ang'/;DATA phy(50)%u /1.d10/
 DATA phy(51)%d /'ef'/;DATA phy(51)%n /'v/bohr'/;DATA phy(51)%u /1.8897268d10/
 DATA phy(52)%d/'ef'/;DATA phy(52)%n /'ry/bohr/e'/;DATA phy(52)%u /2.5711273d11/
 DATA phy(53)%d/'ef'/;DATA phy(53)%n/'har/bohr/e'/;DATA phy(53)%u /5.1422546d11/
 DATA phy(54)%d /'e'/;DATA phy(54)%n /'k'/;DATA phy(54)%u /1.38066d-23/
 DATA phy(55)%d /'t'/;DATA phy(55)%n /'hr'/;DATA phy(55)%u /3600.d0/
 DATA phy(56)%d /'t'/;DATA phy(56)%n /'min'/;DATA phy(56)%u /60.d0/
 DATA phy(57)%d /'t'/;DATA phy(57)%n /'day'/;DATA phy(57)%u /86400.d0/


CONTAINS


 SUBROUTINE esdf_init(filename)
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: filename
 INTEGER unit,ierr,i,j,ic,nt,ndef,nread
 INTEGER,PARAMETER :: ncomm=3,ndiv=3
 LOGICAL inblock
 CHARACTER(1) comment(ncomm),divide(ndiv)
 CHARACTER(llength) cjunk,ctemp
! Added to allow compilation with Sun mpf95 (MDT 1.2003).
 CHARACTER(llength) sunstring

! Define comment characters.
 DATA comment /'#',';','!'/
 DATA divide /' ','=',':'/

 keywords_not_in_list=.false.

! Read the keyword list.
 call load_keywords

! "Reduce" the keyword list for comparison.
 do i=1,numkw
  ctemp=esdf_reduce(kw(i)%label)
  kw(i)%label=ctemp(1:30)
 enddo

! Open the esdf file.
 unit=0 ! or cray compiler complains
 call esdf_file(unit,filename,ierr)
 cjunk='Unable to open main input file "'//trim(filename)//'"'

 if(ierr==1)then
  nread=0
  call errstop('ESDF',trim(cjunk))
 else
  nread=huge(1)
 endif

! Count the number of records (excluding blank lines and commented lines).
 nrecords=0
 do i=1,nread
  read(unit,'(a)',end=100)cjunk
  do j=1,ncomm
   ic=index(cjunk,comment(j))
   if(ic>0)cjunk(ic:)=' '
  enddo
  if(len_trim(cjunk)>0)then
   nrecords=nrecords+1
  endif
 enddo
100 rewind(unit)

! Allocate the array to hold the records and tokens.
 allocate(llist(nrecords),block_data(nrecords),tlist(llength,nrecords),&
  &warns(nrecords),dump(ndump))

! Set the number of warnings to zero.
 nwarns=0 ; warns=' ' ; ndmp=0 ; dump=' '

! Read in the records.
 nrecords=0
 do i=1,nread
  read(unit,'(a)',end=101)cjunk
  do j=1,ncomm
   ic=index(cjunk,comment(j))
   if(ic>0)cjunk(ic:)=' '
  enddo
  if(len_trim(cjunk)>0)then
   nrecords=nrecords+1
   llist(nrecords)=adjustl(cjunk)
  endif
 enddo
101 close(unit)

! Now read in the tokens from llist.
 tlist=' '
 do i=1,nrecords
  ctemp=llist(i)
  nt=0
  do
   if(len_trim(ctemp)<=0)exit
   ic=minval(index(ctemp,divide),mask=index(ctemp,divide)>0)
   if(ic>1)then
    nt=nt+1
    tlist(nt,i)=adjustl(ctemp(:ic-1))
   endif
   ctemp=adjustl(ctemp(ic+1:))
  enddo
 enddo

! Check if any of the "labels" in the input file are unrecognized.
 inblock=.false.
 do i=1,nrecords

! Check if we are in a block.
  if(esdf_reduce(tlist(1,i))=='%block')then
   inblock=.true.
! Check if block label is recognized.
   sunstring=esdf_reduce(tlist(2,i))
   if((count(sunstring==kw%label)==0))then
    ctemp='Label "'//trim(esdf_reduce(tlist(2,i)))//'" not in keyword list'
    if(count(ctemp==warns)==0)then
     keywords_not_in_list=.true.
     call esdf_warn(ctemp)
    endif
   endif
! Check if "label" is multiply defined in the input file.
   ndef=0
   do j=1,nrecords
    if(esdf_reduce(tlist(2,i))==esdf_reduce(tlist(2,j)))ndef=ndef+1
   enddo
   ctemp='Label "'//trim(esdf_reduce(tlist(2,i)))//&
    &'" is multiply defined in the input file. '
   if((ndef>2).and.(count(ctemp==warns)==0))call esdf_warn(ctemp)
  endif

! Check it is in the list of keywords.
  sunstring=esdf_reduce(tlist(1,i))
  if((count(sunstring==kw%label)==0).and.(.not.inblock).and. &
   &sunstring/='%endblock')then
   ctemp='Label "'//trim(esdf_reduce(tlist(1,i)))//'" not in keyword list'
   if(count(ctemp==warns)==0)then
    call esdf_warn(ctemp)
    keywords_not_in_list=.true.
   endif
  endif

  if(.not.inblock)then
! Check if "label" is multiply defined in the input file.
   ndef=0
   do j=1,nrecords
    if(esdf_reduce(tlist(1,i))==esdf_reduce(tlist(1,j)))ndef=ndef+1
   enddo
   ctemp='Label "'//trim(esdf_reduce(tlist(1,i)))//&
    &'" is multiply defined in the input file. '
   if((ndef>1).and.(count(ctemp==warns)==0))call esdf_warn(ctemp)
  endif

! Check if we have left a block.
  if(esdf_reduce(tlist(1,i))=='%endblock')inblock=.false.

 enddo ! i=1,nrecords

 END SUBROUTINE esdf_init


 FUNCTION esdf_string(label,default)
!-----------------------------------------------!
! Return the string attached to the "label".    !
!-----------------------------------------------!
 CHARACTER(*),INTENT(in) :: label,default
 INTEGER i
 CHARACTER(llength) ctemp,esdf_string

! Check "label" is defined.
 call esdf_lblchk(label,'T')

! Set to default.
 esdf_string=default

 do i=1,nrecords
! Search in the first token for "label".
! The first instance is returned.
  if(esdf_reduce(tlist(1,i))==esdf_reduce(label))then
   esdf_string=llist(i)(index(llist(i),trim(tlist(2,i))):)
   exit
  endif
 enddo

! Dump the string used.
 ndmp=ndmp+1
 write(dump(ndmp),*,err=101)trim(esdf_reduce(label)),':',trim(esdf_string)
 if(count(dump(ndmp)==dump(1:ndmp-1))>0)ndmp=ndmp-1

 return

101 ctemp='Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_string'
 call esdf_warn(ctemp)

 END FUNCTION esdf_string


 FUNCTION esdf_integer(label,default)
!-----------------------------------------------!
! Return the integer attached to the "label".   !
!-----------------------------------------------!
 INTEGER,INTENT(in) :: default
 CHARACTER(*),INTENT(in) :: label
 INTEGER i,esdf_integer
 CHARACTER(llength) ctemp

! Check "label" is defined.
 call esdf_lblchk(label,'I')

! Set to default.
 esdf_integer=default

 do i=1,nrecords
! Search in the first token for "label".
! The first instance is returned.
  if(esdf_reduce(tlist(1,i))==esdf_reduce(label))then
   read(tlist(2,i),*,err=100)esdf_integer
   exit
  endif
 enddo

! Dump the value used.
 ndmp=ndmp+1
 write(dump(ndmp),*,err=101)trim(esdf_reduce(label)),':',esdf_integer
 if(count(dump(ndmp)==dump(1:ndmp-1))>0)ndmp=ndmp-1

 return

100 ctemp='Unable to parse "'//trim(esdf_reduce(label))//'" in esdf_integer'
 call esdf_die(ctemp)
101 ctemp='Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_integer'
 call esdf_warn(ctemp)

 END FUNCTION esdf_integer


 FUNCTION esdf_single(label,default)
!----------------------------------------------------------------!
! Return the single precisioned value attached to the "label".   !
!----------------------------------------------------------------!
 REAL(dp),INTENT(in) :: default
 CHARACTER(*),INTENT(in) :: label
 INTEGER i
 REAL(dp) esdf_single
 CHARACTER(llength) ctemp

! Check "label" is defined.
 call esdf_lblchk(label,'S')

! Set to default.
 esdf_single=default

 do i=1,nrecords
! Search in the first token for "label".
! The first instance is returned.
  if(esdf_reduce(tlist(1,i))==esdf_reduce(label))then
   read(tlist(2,i),*,err=100)esdf_single
   exit
  endif
 enddo

! Dump the value used.
 ndmp=ndmp+1
 write(dump(ndmp),*,err=101)trim(esdf_reduce(label)),':',esdf_single
 if(count(dump(ndmp)==dump(1:ndmp-1))>0)ndmp=ndmp-1

 return

100 ctemp='Unable to parse "'//trim(esdf_reduce(label))//'" in esdf_single'
 call esdf_die(ctemp)
101 ctemp='Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_single'
 call esdf_warn(ctemp)

 END FUNCTION esdf_single


 FUNCTION esdf_double(label,default)
!---------------------------------------------------------------!
! Return the double precisioned value attached to the "label".  !
!---------------------------------------------------------------!
 REAL(dp),INTENT(in) :: default
 CHARACTER(*),INTENT(in) :: label
 INTEGER i
 REAL(dp) esdf_double
 CHARACTER(llength) ctemp

! Check "label" is defined.
 call esdf_lblchk(label,'D')

! Set to default.
 esdf_double=default

 do i=1,nrecords
! Search in the first token for "label".
! The first instance is returned.
  if(esdf_reduce(tlist(1,i))==esdf_reduce(label))then
   read(tlist(2,i),*,err=100)esdf_double
   exit
  endif
 enddo

! Dump the value used.
 ndmp=ndmp+1
 write(dump(ndmp),*,err=101)trim(esdf_reduce(label)),':',esdf_double
 if(count(dump(ndmp)==dump(1:ndmp-1))>0) ndmp=ndmp-1

 return

100 esdf_double=default
 ctemp='Unable to parse "'//trim(esdf_reduce(label))//'" in esdf_double'
 call esdf_die(ctemp)
101 ctemp='Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_double'
 call esdf_warn(ctemp)

 END FUNCTION esdf_double


 FUNCTION esdf_physical(label,default,dunit)
!-------------------------------------------------------------------------!
! Return the double precisioned physical value attached to the "label".   !
! Units converted to "dunit".                                             !
!-------------------------------------------------------------------------!
 REAL(dp),INTENT(in) :: default
 CHARACTER(*),INTENT(in) :: label,dunit
 INTEGER i
 REAL(dp) esdf_physical
 CHARACTER(llength) ctemp,iunit

! Check "label" is defined.
 call esdf_lblchk(label,'P')

! Set to default.
 esdf_physical=default

 do i=1,nrecords
! Search in the first token for "label".
! The first instance is returned.
  if(esdf_reduce(tlist(1,i))==esdf_reduce(label))then
   read(tlist(2,i),*,err=100,end=100)esdf_physical
   read(tlist(3,i),*,err=100,end=100)iunit
   esdf_physical=esdf_convfac(iunit,dunit)*esdf_physical
   exit
  endif
 enddo

! Dump the value used.
 ndmp=ndmp+1
 write(dump(ndmp),*,err=101)trim(esdf_reduce(label)),':',esdf_physical,' ',&
 &trim(dunit)
 if(count(dump(ndmp)==dump(1:ndmp-1))>0)ndmp=ndmp-1

 return

100 esdf_physical=default
 ctemp='Unable to parse "'//trim(esdf_reduce(label))//'" in esdf_physical'
 call esdf_die(ctemp)
101 ctemp='Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_physical'
 call esdf_warn(ctemp)

 END FUNCTION esdf_physical


 FUNCTION esdf_defined(label,type)
!-------------------------------------------!
! Is the "label" defined in the input file? !
!-------------------------------------------!
 CHARACTER(*),INTENT(in) :: label
 CHARACTER(1),INTENT(in) :: type
 INTEGER i
 LOGICAL esdf_defined
 CHARACTER(llength) ctemp

! Check "label" is defined.
 call esdf_lblchk(label,type)

! Set to default.
 esdf_defined=.false.

 do i=1,nrecords
! Search in the first token for "label".
! The first instance is returned.
  if(esdf_reduce(tlist(1,i))==esdf_reduce(label))then
   esdf_defined=.true.
   exit
  endif
 enddo

! Dump the value used.
 if(esdf_defined)then
  ndmp=ndmp+1
  write(dump(ndmp),*,err=101)trim(esdf_reduce(label)),':'
  if(count(dump(ndmp)==dump(1:ndmp-1))>0)ndmp=ndmp-1
 endif

 return

101 ctemp='Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_defined'
 call esdf_warn(ctemp)

 END FUNCTION esdf_defined


 FUNCTION esdf_boolean(label,default)
!----------------------------------------------!
! Is the "label" defined in the input file?    !
!----------------------------------------------!
 CHARACTER(*),INTENT(in) :: label
 LOGICAL,INTENT(in) :: default
 INTEGER i
 LOGICAL esdf_boolean
 CHARACTER(llength) ctemp,positive(3),negative(3)
! Added to allow compilation with Sun mpf95 (MDT 1.2003).
 CHARACTER(llength) sunstring
 DATA positive /'yes','true','t'/
 DATA negative /'no','false','f'/

! Check "label" is defined.
 call esdf_lblchk(label,'L')

! Set to default.
 esdf_boolean=default

 do i=1,nrecords
! Search in the first token for "label".
! The first instance is returned.
  if(esdf_reduce(tlist(1,i))==esdf_reduce(label))then
   if(len_trim(tlist(2,i))==0)then
    esdf_boolean=.true.
    exit
   endif
   sunstring=esdf_reduce(tlist(2,i))
!   if(any(index(positive,esdf_reduce(tlist(2,i)))>0))then
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

! Dump the value used.
 ndmp=ndmp+1
 write(dump(ndmp),*,err=101)trim(esdf_reduce(label)),': ',esdf_boolean
 if(count(dump(ndmp)==dump(1:ndmp-1))>0) ndmp=ndmp-1

 return

101 ctemp='Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_boolean'
 call esdf_warn(ctemp)

 END FUNCTION esdf_boolean


 FUNCTION esdf_block(label,nlines)
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: label
 INTEGER,INTENT(out) :: nlines
 INTEGER i,j
 LOGICAL esdf_block
 CHARACTER(llength) ctemp

! Check "label" is defined.
 call esdf_lblchk(label,'B')
 ctemp='Block "'//trim(esdf_reduce(label))//'" not closed correctly '
 esdf_block=.false.
 nlines=0

 do i=1,nrecords
  if((esdf_reduce(tlist(1,i))==esdf_reduce('%block'))&
   &.and.(esdf_reduce(tlist(2,i))==esdf_reduce(label)))then
   esdf_block=.true.
   do
    if(esdf_reduce(tlist(1,i+nlines+1))==esdf_reduce('%endblock'))exit
    nlines=nlines+1
    if(nlines+i>nrecords)call esdf_die(ctemp)
    block_data(nlines)=llist(i+nlines)
   enddo
   if(esdf_reduce(tlist(2,i+nlines+1))/=esdf_reduce(label))call esdf_die(ctemp)
   exit
  endif
 enddo
 if(.not.esdf_block)return

! Dump the block.
 ndmp=ndmp+1
 write(dump(ndmp),*,err=101)'%block ',trim(esdf_reduce(label)),': '
 if(count(dump(ndmp)==dump(1:ndmp-1))>0)then
  ndmp=ndmp-1
  return
 endif
 do j=1,nlines
  ndmp=ndmp+1
  if(ndmp>ndump)call errstop('ESDF_BLOCK','Too many lines in block : increase&
   & the parameter NDUMP in the esdf.f90 module')
  dump(ndmp)=block_data(j)
 enddo
 ndmp=ndmp+1
 write(dump(ndmp),*,err=101)'%endblock ',trim(esdf_reduce(label)),': '

 return

101 ctemp='Unable to dump "'//Trim(esdf_reduce(label))//'" in esdf_block'
 call esdf_warn(ctemp)

 END FUNCTION esdf_block


 FUNCTION esdf_reduce(string_untrimmed)
!---------------------------------------------------------!
! Reduce the string to lower case and remove punctuation. !
!---------------------------------------------------------!
 CHARACTER(*),INTENT(in) :: string_untrimmed
 INTEGER iA,iZ,ishift,ic,i,ln
 INTEGER,PARAMETER :: npunct=3
 CHARACTER(1) punct(npunct)
 CHARACTER(llength) ctemp,string,esdf_reduce
! Define the punctuation to be removed.
 DATA punct /'.','_','-'/

! Initialize system dependant bounds in collating sequence.
 iA=ichar('A') ; iZ=Ichar('Z')
 ishift=ichar('a')-iA

! Initialize output.
 string=adjustl(string_untrimmed)
 ln=len_trim(string)
 esdf_reduce(1:ln)=string(1:ln) ; esdf_reduce(ln+1:)=' '

! Drop all upper case characters to lower case.
 do i=1,llength
  ic=ichar(esdf_reduce(i:i))
  if((ic>=iA).and.(ic<=iZ))esdf_reduce(i:i)=char(ishift+ic)
 enddo

! Now remove punctuation.
 do i=1,npunct
  do
   ic=index(esdf_reduce,punct(i))
   if(ic>0)then
    ctemp=esdf_reduce
    esdf_reduce(ic:ln)=ctemp(ic+1:ln)//' '
   else
    exit
   endif
  enddo
 enddo
 esdf_reduce=trim(adjustl(esdf_reduce))

 END FUNCTION esdf_reduce


 FUNCTION esdf_convfac(from,to)
!-----------------------------------------------------!
! Find the conversion factor between physical units.  !
!-----------------------------------------------------!
 CHARACTER(*),INTENT(in) :: from,to
 INTEGER i,ifrom,ito
 REAL(dp) esdf_convfac
 CHARACTER(llength) ctemp

! Find the index numbers of the from and to units.
 ifrom=0 ; ito=0
 do i=1,nphys
  if(esdf_reduce(from)==phy(i)%n)ifrom=i
  if(esdf_reduce(to)==phy(i)%n)ito=i
 enddo

! Check that the units were recognized.
 if(ifrom==0)then
  ctemp='Units not recognized in input file : '//trim(esdf_reduce(from))
  call esdf_die(ctemp)
 endif
 if(ito==0)then
  ctemp='Units not recognized in program : '//trim(esdf_reduce(to))
  call esdf_die(ctemp)
 endif

! Check that from and to are of the same dimensions.
 if(phy(ifrom)%d/=phy(ito)%d)then
  ctemp='Dimensions do not match : '//trim(esdf_reduce(from))&
   & //' vs '//trim(esdf_reduce(to))
  call esdf_die(ctemp)
 endif

! Set the conversion factor.
 esdf_convfac=phy(ifrom)%u/phy(ito)%u

 END FUNCTION esdf_convfac


 FUNCTION esdf_unit(ierr)
!---------------------------!
! Find an unused i/o unit.  !
!---------------------------!
 INTEGER,INTENT(out) :: ierr
 INTEGER esdf_unit
 LOGICAL op

 ierr=0
 do esdf_unit=10,99
  inquire(unit=esdf_unit,opened=op,err=100)
  if(.not.op)return
 enddo
 call esdf_warn('Unable to find a free i/o unit using esdf_unit.')
 ierr=1
 return

100 call esdf_die('Error opening files by esdf_unit.')

 END FUNCTION esdf_unit


 SUBROUTINE esdf_file(unit,filename,ierr)
!-------------------!
! Open an old file. !
!-------------------!
 INTEGER,INTENT(out) :: unit,ierr
 CHARACTER(*),INTENT(in) :: filename
 LOGICAL ex

 unit=esdf_unit(ierr)
 if(ierr>0)return
 inquire(file=trim(filename),exist=ex,err=100)
 if(.not.ex)goto 100
 open(unit=unit,file=trim(filename),form='formatted',status='old',err=100)
 return

100 ierr=1

 END SUBROUTINE esdf_file


 SUBROUTINE esdf_newfile(unit,filename,ierr)
!--------------------!
! Open a new file.   !
!--------------------!
 CHARACTER(*),INTENT(in) :: filename
 INTEGER,INTENT(out) :: unit,ierr

 unit=esdf_unit(ierr)
 if(ierr>0)return
 open(unit=unit,file=trim(filename),form='formatted',status='replace',err=100)
 return

100 ierr=1

 END SUBROUTINE esdf_newfile


 SUBROUTINE esdf_lblchk(string,typ)
!----------------------------------------------------!
! Check that the label is known, and used correctly. !
!----------------------------------------------------!
 CHARACTER(*),INTENT(in) :: string
 CHARACTER(1),INTENT(in) :: typ
 INTEGER i
 LOGICAL match(numkw)
 CHARACTER(llength) ctemp

! Check if label is recognized.
 match=(esdf_reduce(string)==kw(:)%label)
 i=count(match)
 if(i==0)then
  ctemp='Label "'//trim(esdf_reduce(string))//'" not recognized in keyword list'
  call esdf_die(ctemp)
 endif
 if(i>1)then
  ctemp='Label "'//trim(esdf_reduce(string))//'" is multiply defined'
  call esdf_die(ctemp)
 endif
! compare type of matching keyword.
 where(match)
  match(:)=(typ/=kw(:)%typ(1:1))
 endwhere
 if(any(match))then
  ctemp='Label "'//trim(esdf_reduce(string))//'" has been used with wrong type'
  call esdf_die(ctemp)
 endif

 END SUBROUTINE esdf_lblchk


 SUBROUTINE esdf_help(helpword_in,searchword_in)
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: helpword_in,searchword_in
 INTEGER i,indx,indx2,ln
 CHARACTER(1) cl
 CHARACTER(20) ctyp,clev
 CHARACTER(60) title
 CHARACTER(78) ctemp
 CHARACTER(llength) helpword,searchword
 CHARACTER(llength),ALLOCATABLE :: kcheck(:)

 helpword=esdf_reduce(helpword_in)
 searchword=esdf_reduce(searchword_in)

 if(trim(helpword)=='search')then
  if(len_trim(searchword)<1)call esdf_die('"searchword" is empty.')
! Search for useful keywords.
  do i=1,numkw
   if((index(kw(i)%label,trim(searchword))>0).or.&
    &(index(kw(i)%dscrpt,trim(searchword))>0))then
    indx=index(kw(i)%dscrpt,'!*')-1
    if(indx==-1)call esdf_die('Keyword description incorrectly formatted.')
    title=kw(i)%dscrpt(1:indx)
    ln=len_trim(title)
    if(ln>80) call esdf_die('Keyword title too long.')
    write(o,*)kw(i)%label,Trim(title)
   endif
  enddo
  write(o,*)
  call all_stop
 endif

! All keywords, short description
 if(trim(helpword)=='all')then
  do i=1,numkw
   if(len_trim(kw(i)%label)>0)then
    indx=index(kw(i)%dscrpt,'!*')-1
    if(indx==-1)call esdf_die('Keyword description incorrectly formatted.')
    title=kw(i)%dscrpt(1:indx)
    ln=len_trim(title)
    if(ln>80)call esdf_die('Keyword title too long.')
    write(o,*)kw(i)%label,trim(title)
   endif
  enddo
  write(o,*)
  call all_stop
 endif

! All specific levels of keywords
 if(any((/'basic ','inter ','expert','dummy '/)==trim(helpword)))then
  select case(trim(helpword))
   case('basic')  ; cl='B'
   case('inter')  ; cl='I'
   case('expert') ; cl='E'
   case('dummy')  ; cl='D'
  end select
  do i=1,numkw
   if(kw(i)%typ(3:3)==cl)then
    indx=index(kw(i)%dscrpt,'!*')-1
    if(indx==-1)call esdf_die('Keyword description incorrectly formatted.')
    title=kw(i)%dscrpt(1:indx)
    ln=len_trim(title)
    if(ln>80) call esdf_die('Keyword title too long.')
    write(o,*)kw(i)%label,Trim(title)
   endif
  enddo
  write(o,*)
  call all_stop
 endif

! More information about a specific keyword
! MDT (Original didn't allow punctuation in keywords e.g. vm_debug)
 allocate(kcheck(numkw))
 do i=1,numkw
  kcheck(i)=esdf_reduce(kw(i)%label)
 enddo
 if(.not.any(kcheck==trim(helpword)))call esdf_die('Keyword not recognized.')
 if(count(kcheck==trim(helpword))>1)call esdf_die('Keyword entry duplicated.')
 deallocate(kcheck)
 do i=1,numkw
  if(esdf_reduce(kw(i)%label)==trim(helpword))then
   write(o,*)'Keyword : ',trim(kw(i)%label)
   indx=index(kw(i)%dscrpt,'!*')+1
   if(indx==1)then
    write(o,*)'Title   : (unknown)'
   else
    title=trim(adjustl(kw(i)%dscrpt(1:indx-2)))
    indx2=index(title,'*!')
    if(indx2>0)title=trim(adjustl(title(indx2+2:)))
    write(o,*)'Title   : ',trim(title)
   endif
   select case(kw(i)%typ(1:1))
    case('I') ; ctyp='Integer'
    case('S') ; ctyp='Single Precision'
    case('D') ; ctyp='Double Precision'
    case('P') ; ctyp='Physical'
    case('T') ; ctyp='String'
    case('E') ; ctyp='Defined'
    case('B') ; ctyp='Block'
    case('L') ; ctyp='Boolean'
   end select
   select case(kw(i)%typ(3:3))
    case('B') ; clev='Basic'
    case('I') ; clev='Intermediate'
    case('E') ; clev='Expert'
    case('D') ; clev='Dummy'
   end select
   write(o,*)'Type    : ',trim(ctyp)
   write(o,*)'Level   : ',trim(clev)
   write(o,*)
   write(o,*)'DESCRIPTION'
   write(o,*)'-----------'
   indx=indx+1 ; ln=len_trim(kw(i)%dscrpt)
   do while(indx<ln)
    ctemp=kw(i)%dscrpt(indx:min(indx+77,ln))
    indx2=index(ctemp,'$')-1
    if(indx2==-1)indx2=index(ctemp,' ',back=.true.)
    write(ctemp,'(a)')adjustl(ctemp(:indx2))
    write(o,*)trim(ctemp)
    indx=indx+indx2
    if(kw(i)%dscrpt(indx:indx)=='$')indx=indx+1
   enddo
  endif
 enddo

 write(o,*)
 call all_stop

 END SUBROUTINE esdf_help


 SUBROUTINE esdf_die(string)
!-------------------------------------------------!
! Stop execution due to an error caused by esdf.  !
!-------------------------------------------------!
 CHARACTER(*),INTENT(in) :: string
 call errstop_quiet('ESDF input',trim(string))
 END SUBROUTINE esdf_die


 SUBROUTINE esdf_warn(string)
!----------------------------------------!
! Warning due to an error cause by esdf. !
!----------------------------------------!
 CHARACTER(*),INTENT(in) :: string
 nwarns=nwarns+1
 warns(nwarns)=string
 END SUBROUTINE esdf_warn


 SUBROUTINE esdf_warnout
!---------------------------------------!
! Dump the warnings to screen and quit. !
!---------------------------------------!
 INTEGER i
 do i=1,nwarns
  write(o,*)'INPUT ERROR: '//trim(warns(i))
  write(o,*)
 enddo
 if(keywords_not_in_list)then
  call errstop('ESDF','Unknown keywords detected in LOUIS input file.')
 endif
 END SUBROUTINE esdf_warnout


 SUBROUTINE esdf_close
!-------------------------------------------------------------------!
! Deallocate the data arrays --- call this before re-initializing.  !
!-------------------------------------------------------------------!
 deallocate(llist,tlist,block_data,warns,dump,kw)
 END SUBROUTINE esdf_close


 SUBROUTINE esdf_dump(filename)
!------------------------------------------------------!
! Dump an input file which contains all set variables  !
! including defaults.                                  !
!------------------------------------------------------!
 CHARACTER(*),INTENT(in) :: filename
 INTEGER unit_dump,ierr,i,j,indx
 CHARACTER(llength) cjunk
 CHARACTER(llength),PARAMETER :: spaces=repeat(' ',llength)

! Open the input_example file.
 call esdf_newfile(unit_dump,filename,ierr)
 if(ierr==1)call esdf_die('Unable to open example input file "input_example"')

! indx=maxval(index(dump(1:ndmp),':'))
 indx=19 ! input default

 do i=1,ndmp
  j=index(dump(i),':')
  if(j>0)then
   cjunk=dump(i)(1:j-1)
! MDT Note : this causes IBM compiler to crash.
!   dump(i)=trim(cjunk)//repeat(' ',indx-j+1)//': '&
!   & //trim(adjustl(dump(i)(j+1:)))//'#'
   dump(i)=trim(cjunk)//spaces(1:indx-j+1)//': '&
   & //trim(adjustl(dump(i)(j+1:)))//'#'
  endif
 enddo
 indx=maxval(index(dump(1:ndmp),'#',back=.true.))
 indx=36 ! input default

 do i=1,ndmp
  j=index(dump(i),'#',back=.true.)
  if(j>0)then
! MDT Note : this causes IBM compiler to crash.
!  dump(i)=dump(i)(1:j-1)//repeat(' ',indx-j+1)//'#'
   dump(i)=dump(i)(1:j-1)//spaces(1:indx-j+1)//'#'
  endif
 enddo

 do i=1,ndmp
  j=index(dump(i),':')
  if(j>0)then
   cjunk=dump(i)(1:j-1)
   do j=1,numkw
    if(index(cjunk,trim(kw(j)%label))>0)exit
   enddo
   select case(kw(j)%typ(1:1))
    case('I') ; cjunk='Integer'
    case('S') ; cjunk='Single Precision'
    case('D') ; cjunk='Double Precision'
    case('P') ; cjunk='Physical'
    case('T') ; cjunk='String'
    case('E') ; cjunk='Defined'
    case('B') ; cjunk='Block'
    case('L') ; cjunk='Boolean'
   end select
   indx=index(kw(j)%dscrpt,'!*')
   dump(i)=trim(dump(i))//trim(kw(j)%dscrpt(1:indx-1))//' &
   &('//trim(adjustl(cjunk))//')'
  endif
 enddo

 do i=1,ndmp
  write(unit_dump,'(a)')adjustl(dump(i))
 enddo

 END SUBROUTINE esdf_dump


 SUBROUTINE help_system(hword,sword)
!-------------------------------------------------------------------------!
! Read the temporary file produced by a help script and read the          !
! keywords. Access the esdf help system, output the help, then die.       !
! MDT 8.2000.  Brought here from help.f90 PLR 10.2006.                    !
!-------------------------------------------------------------------------!
  IMPLICIT NONE
  CHARACTER(*),OPTIONAL,INTENT(in) :: hword,sword
  CHARACTER(llength) helpword,searchword

  searchword=''
  call load_keywords
  if(.not.(present(hword)))return
  helpword=hword ; searchword=''
  if(present(sword))searchword=sword
  write(o,*)'LOUIS HELP SYSTEM'
  write(o,*)'================='
  write(o,*)
  call esdf_help(trim(helpword),trim(searchword))
  write(o,*)

 END SUBROUTINE help_system


END MODULE esdf

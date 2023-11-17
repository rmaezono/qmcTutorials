MODULE global
  INTEGER, PARAMETER :: DP=KIND(1.0D0)
  REAL (KIND=DP),PARAMETER :: ha_to_eV=27.2113962d0
END MODULE global

PROGRAM egaussian
  USE global
  IMPLICIT none
! Scan through a G94, G98 or G03 output file and pull-out the SCF energy
! and the KE, PE etc. (if present)

  INTEGER, PARAMETER :: NF_max=20 ! Max no. of fields in a string
  INTEGER, PARAMETER :: io_in=21
  CHARACTER(LEN=30)  :: output_file
  CHARACTER(LEN=80)  :: linestring
  CHARACTER(LEN=80)  :: instring 
  CHARACTER(LEN=80)  :: fields(NF_max) ! Holds the fields pulled from string
  INTEGER            :: NF ! Holds the number of fields found in a string
  INTEGER            :: icheck=0  ! IO checking
  INTEGER            :: length    ! Length of a field string
  REAL(KIND=DP)      :: scf_energy,ke,pe,ee
  LOGICAL            :: verbose
  INTERFACE
   SUBROUTINE pattern_search(io_unit,keystring,instring,ierr)
    INTEGER, INTENT(in) :: io_unit
    CHARACTER(*),INTENT(in) :: keystring
    CHARACTER(*),INTENT(inout) :: instring
    INTEGER,INTENT(out) :: ierr
   END SUBROUTINE pattern_search
   SUBROUTINE pattern_search_last(io_unit,keystring,instring,tempstring,ierr)
    INTEGER, INTENT(in) :: io_unit
    CHARACTER(*),INTENT(in) :: keystring
    CHARACTER(*),INTENT(inout) :: instring
    CHARACTER(*),INTENT(inout) :: tempstring
    INTEGER,INTENT(out) :: ierr
   END SUBROUTINE pattern_search_last
   SUBROUTINE scan_string(instring,fields,NF,NF_max)
    INTEGER,INTENT(out) :: NF
    INTEGER,INTENT(in) :: NF_max
    CHARACTER(*),INTENT(inout) :: instring
    CHARACTER(*),INTENT(inout) :: fields(NF_max)
   END SUBROUTINE scan_string
  END INTERFACE

  WRITE(6,FMT="('Enter full name of Gaussian output file:')")
  READ(5,FMT="(A30)") output_file

  OPEN(unit=io_in,file=TRIM(ADJUSTL(output_file)),status='old',IOSTAT=icheck)

  IF(icheck /= 0)THEN
    WRITE(6,FMT="('Cannot find/open ',A30)") TRIM(ADJUSTL(output_file))
    STOP
  END IF

! Check to see whether Gaussian was told to give verbose output
  CALL pattern_search(io_in,'#P ',linestring,icheck)
  IF(icheck ==0)THEN
    verbose=.true.
  ELSE
    REWIND(io_in)
    CALL pattern_search(io_in,'#p ',linestring,icheck)
    IF(icheck ==0)THEN
     verbose=.true.
    ELSE   
     verbose=.false.
     REWIND(io_in)
    END IF
  END IF

! Find the line on which 'SCF D' occurs
  CALL pattern_search_last(io_in,' SCF D',instring, linestring, icheck)
  IF(LEN_TRIM(instring) /= 0)THEN
    CALL scan_string(instring,fields,NF,NF_max)
    READ(fields(5),*) scf_energy
!    scf_energy=Ha_to_eV*scf_energy ! Convert to eV
    WRITE(6,FMT="(/'SCF tot. energy = ',F15.5,' au')") scf_energy
  ELSE
    WRITE(6,FMT="('End of file while searching for ''SCF D''.')")
    STOP
  END IF
! Now find the line on which the energy break-down is given
  IF(verbose)THEN
    CALL pattern_search_last(io_in,' KE',instring, linestring, icheck)
    IF(LEN_TRIM(instring) /= 0)THEN
      CALL scan_string(instring,fields,NF,NF_max)
      READ(fields(2),*) ke
      length=LEN_TRIM(fields(3))
      READ(fields(3)(4:length),*) pe
      READ(fields(NF),*) ee
! Convert to eV
!      ke=Ha_to_eV*ke; pe=Ha_to_eV*pe; ee=Ha_to_eV*ee
      WRITE(6,FMT="('G.S. KE         = ',F15.5,' au')") ke
      WRITE(6,FMT="('G.S. e-I PE     = ',F15.5,' au')") pe
      WRITE(6,FMT="('G.S. e-e PE     = ',F15.5,' au')") ee
    ELSE
      WRITE(6,FMT="('Cannot find ''KE'' yet verbosity flag is on.')")
      STOP
    END IF
  END IF

END PROGRAM egaussian

SUBROUTINE pattern_search(io_unit,keystring,instring,ierr)
  IMPLICIT none
  INTEGER, INTENT(in) :: io_unit   ! Unit no. of file to read
  CHARACTER(*), INTENT(in)    :: keystring ! The string to search for
  CHARACTER(*), INTENT(inout) :: instring  ! String containing the line
                                           ! on which keystring occurs
  INTEGER, INTENT(out) :: ierr ! Set >0 if end of file occurs

! Scan forward through the file identified by io_unit to find the
! string stored in keystring.  Return the iposn'th element on this
! line.
  INTEGER       :: iposn
  INTEGER       :: ilen      ! Length of instring
  INTEGER       :: icheck    ! IO checking
  LOGICAL       :: search
  CHARACTER(LEN=6) :: fmtstring ! Holds the format string created so
                                ! that we always read-in the correct
                                ! line length to fill instring
  ilen=LEN(instring)
  WRITE(fmtstring,FMT="(I4)") ilen
  fmtstring="(A"//TRIM(ADJUSTL(fmtstring))//")"

  ierr=0
  search=.true.
  DO WHILE(search)
    READ(io_unit,FMT=fmtstring,IOSTAT=icheck) instring
    IF(icheck .lt. 0)THEN
! We have reached the end of the file and not found the desired string
! so flag this and return to calling program unit
      ierr=1
      EXIT
    END IF
    iposn=INDEX(STRING=TRIM(ADJUSTL(instring)),SUBSTRING=keystring)
    IF(iposn /= 0) search=.false.
  END DO

END SUBROUTINE pattern_search

SUBROUTINE pattern_search_last(io_unit,keystring,instring,tempstring,ierr)
  IMPLICIT none
  INTEGER, INTENT(in) :: io_unit   ! Unit no. of file to read
  CHARACTER(*), INTENT(in)    :: keystring ! The string to search for
  CHARACTER(*), INTENT(inout) :: instring  ! String containing the line
                                           ! on which keystring occurs
  CHARACTER(*), INTENT(inout) :: tempstring  ! Temporary string to store lines 
  INTEGER, INTENT(out) :: ierr ! Set >0 if end of file occurs

! Scan forward through the file identified by io_unit to find last occurence
! of the line beginning with keystring (omitting initial spaces))

  INTEGER       :: klen      ! length of keystring
  INTEGER       :: tlen      ! Length of tempstring
  INTEGER       :: icheck    ! IO checking
  LOGICAL       :: search
  CHARACTER(LEN(keystring))   keystring_trim
  CHARACTER(LEN(tempstring))   tempstring_trim
  CHARACTER(LEN=6) :: fmtstring ! Holds the format string created so
                                ! that we always read-in the correct
                                ! line length to fill instring

  REWIND(io_unit)

  instring=" " !to have LEN_TRIM=0 if nothing found

  tlen=LEN(tempstring)
  WRITE(fmtstring,FMT="(I4)") tlen
  fmtstring="(A"//TRIM(ADJUSTL(fmtstring))//")"

  keystring_trim=TRIM(ADJUSTL(keystring))
  klen=LEN_TRIM(keystring_trim)
  ierr=0
  search=.true.
  DO WHILE(search)
    READ(io_unit,FMT=fmtstring,IOSTAT=icheck) tempstring
    tempstring_trim=ADJUSTL(tempstring)
    IF(icheck .lt. 0)THEN
! We have reached the end of the file and not found the desired string
! so flag this and return to calling program unit
      ierr=1
      EXIT
    END IF
!    write(*,*) tempstring_trim(1:klen), klen
    IF(tempstring_trim(1:klen)==keystring_trim)THEN
     instring=tempstring_trim
    END IF
  END DO

END SUBROUTINE pattern_search_last


SUBROUTINE scan_string(instring,fields,NF,NF_max)
  IMPLICIT none
  INTEGER, INTENT(out)      :: NF
  INTEGER, INTENT(in)       :: NF_max
  CHARACTER(*), INTENT(inout) :: instring
  CHARACTER(*), INTENT(inout) :: fields(NF_max)

! Parse the string 'instring' and store the fields (objects separated
! by spaces) in the fields array.

  INTEGER :: ilength ! Length of string being scanned (excluding
                     ! trailing spaces)
  INTEGER :: icount  ! Counter for no. of fields found
  INTEGER :: iposn   ! Current character of string being scanned
  INTEGER :: ifpt    ! Current character of field being assigned
  LOGICAL :: newfield! Flag to say whether we are currently within a
                     ! field or between fields (i.e. the last character
                     ! scanned was a space)

  instring=TRIM(ADJUSTL(instring))
  ilength=LEN_TRIM(instring)

  iposn=0
  icount=0
  newfield=.true.
  ! Initialise fields as the elements are constructed character by character
  fields=" "

  DO WHILE(iposn.lt.ilength)
    iposn=iposn + 1
    IF( instring(iposn:iposn) /= " " )THEN

      IF(newfield)THEN
        newfield=.false.
        icount=icount + 1
        IF(icount .gt. NF_max)THEN
          WRITE(6,FMT="('Only ',I2,' fields per line allowed.  Increase &
                         &NF_max in routine calling scan_string.f90.')") NF_max
          STOP
        END IF
        fields(icount)(1:1)=instring(iposn:iposn)
        ifpt=1
      ELSE
        ifpt=ifpt + 1
        fields(icount)(ifpt:ifpt)=instring(iposn:iposn)
      END IF

    ELSE
      IF(.not. newfield) newfield=.true.
    END IF

  END DO

! Return the no. of fields found
  NF=icount

END SUBROUTINE scan_string



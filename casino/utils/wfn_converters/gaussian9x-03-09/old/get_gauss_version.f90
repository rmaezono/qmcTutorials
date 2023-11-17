SUBROUTINE get_gauss_version
!-------------------------------------------------------------------------!
! Simply gets which version of Gaussian was used to run this particular   !
! job.  Could be done in read_G9Xout but is simpler to do it alone here   !
! for ordering purposes.                                                  !
!-------------------------------------------------------------------------!
  USE awk_like
  USE g94_wavefunction
  IMPLICIT none
  INTEGER, PARAMETER :: io_in = 33 ! Unit no. to read from
  CHARACTER*35       :: input_file
  CHARACTER*130      :: instring ! Storage for each line of output
                                 ! file, prior to parsing
  INTEGER            :: icheck, ifail


  input_file=TRIM(g94_file)//'.out'

  OPEN(unit=io_in,file=input_file,status='old', IOSTAT=icheck)

  IF(icheck .gt. 0)THEN
    WRITE (*,FMT="(/'Cannot open file ',A35,' - where is it?')") input_file
    STOP
  END IF

! Identify whether this was a G94 or G98 calculation using the copyright
! statement that Gaussian always prints
  CALL getline(io_in,'Copyright',fields,NF,NF_max,ifail,.true.)
  IF(ifail < 0)CALL fatal('Failed to read ''Copyright'' from output file')

  instring=fields(3) ! This holds all of the copyright dates
  CALL scan_string(TRIM(instring),fields,NF,NF_max,",")
  ! Gaussian 98 has one more date in its list
  IF(NF .eq. 5)THEN
    code_used="Gaussian 94"
  ELSE IF(NF .eq. 6)THEN
    code_used="Gaussian 98"
  ELSE IF(NF .eq. 7)THEN
    code_used="Gaussian 03"
  ELSE
   write(6,*)"Problem determining which version of Gaussian produced this file."
   stop
  END IF

  CLOSE(io_in)
END SUBROUTINE get_gauss_version

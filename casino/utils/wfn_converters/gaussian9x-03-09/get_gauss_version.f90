SUBROUTINE get_gauss_version
!-------------------------------------------------------------------------!
! Gets which version of Gaussian was used to run this particular job.     !
! Could be done in read_G9Xout but is simpler to do it alone here         !
! for ordering purposes.                                                  !
!-------------------------------------------------------------------------!
 USE awk_like
 USE g94_wavefunction
 IMPLICIT NONE
 INTEGER icheck,ifail
 INTEGER,PARAMETER :: io_in = 33 ! Unit number to read from
 CHARACTER(35) :: input_file
 CHARACTER(130) :: instring ! Storage for each line of out file prior to parsing

 input_file=trim(g94_file)//'.out'

 open(unit=io_in,file=input_file,status='old',iostat=icheck)

 if(icheck>0)then
  write(*,fmt="(/'Cannot open file ',a35,'.')")input_file
  stop
 endif

! Identify whether this was a G94/G98/G03 calculation using the copyright
! statement that Gaussian always prints
 call getline(io_in,'Copyright',fields,NF,NF_max,ifail,.true.)
 if(ifail<0)call fatal('Failed to read ''Copyright'' from output file')

 instring=fields(3) ! This holds all of the copyright dates
 call scan_string(trim(instring),fields,NF,NF_max,",")
! Gaussian 98 has one more date in its list
 if(NF==5)then
  code_used="Gaussian 94"
 elseif(NF==6)then
  code_used="Gaussian 98"
 elseif(NF==7)then
  code_used="Gaussian 03"
! Later revisions of Gaussian 03 have 2004 as the last date in the Copyright
! statement, which brings the total number of fields in this string to 8.
! RAndC b: consider release E.01 (2007).
 elseif(fields(8)=='2004')then
  code_used="Gaussian 03, Revision C or later"
 elseif(fields(8)=='2009')then
  code_used="Gaussian 09"
! RAndC e
 else
  write(6,*)"Problem determining which version of Gaussian produced this file."
  stop
 endif

 close(io_in)

END SUBROUTINE get_gauss_version

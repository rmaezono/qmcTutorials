PROGRAM champer
!-----------------------------------------------------------------!
! CHAMP to CASINO pseudopotential converter                       !
! =========================================                       !
!                                                                 !
! Read CHAMP pseudopotential files and convert them to xx_pp.data !
! files suitable for use with CASINO.                             !
!                                                                 !
! Mike Towler, Oct 2013                                           !
!                                                                 !
! Changes                                                         !
! -------                                                         !
! None yet.                                                       !
!-----------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,n,npotd,zatom,grid_type,ngrid_points,ierr,ialloc
 INTEGER,PARAMETER :: io=10,dp=kind(1.d0)
 REAL(dp) zion,r_asymp,r0,expnt,exp_expnt
 REAL(dp),ALLOCATABLE :: ri(:),psp(:,:)
 LOGICAL file_exists
 CHARACTER(14) outfile
 CHARACTER(80) comment,infile

! Periodic table up to uranium by atomic number, for file names
 CHARACTER(2),PARAMETER :: periodic_table_nocap(0:92)=&
  &(/'wg','h ','he','li','be','b ','c ','n ',&
  & 'o ','f ','ne','na','mg','al','si','p ','s ','cl','ar',&
  & 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',&
  & 'zn','ga','ge','as','se','br','kr','rb','sr','y ','zr',&
  & 'nb','mo','tc','ru','rh','pd','ag','cd','in','sn','sb',&
  & 'te','i ','xe','cs','ba','la','ce','pr','nd','pm','sm',&
  & 'eu','gd','tb','dy','ho','er','tm','yb','lu','hf','ta',&
  & 'w ','re','os','ir','pt','au','hg','tl','pb','bi','po',&
  & 'at','rn','fr','ra','ac','th','pa','u '/)

 write(6,*)
 write(6,*)'CHAMP -> CASINO Pseudopotential Converter'
 write(6,*)
 do
  write(6,*)'Name of CHAMP pseudopotential file to convert:'
  read(5,*,iostat=ierr)infile
  if(ierr==0)exit
  write(6,*)'Invalid name.'     
 enddo
 write(6,*)
 infile=adjustl(infile)
 write(6,*)'Reading file '//trim(infile)//'.'

 open(unit=io,file=infile,form='formatted',status='old',iostat=ierr)
 if(ierr/=0)then
  write(6,*)'Can''t open file '//trim(infile)//'.'
  stop
 endif

! Get atomic number.
 write(6,*)
 write(6,*)'The CHAMP pseudopotential format does not contain the atomic &
  &number,'
 write(6,*)'but this is required by CASINO.'
 write(6,*)
 do
  write(6,*)'Input the atomic number of the atom represented by this pseudo:'
  read(5,*,iostat=ierr)zatom
  if(ierr/=0)zatom=-1
  if(zatom<1.or.zatom>92)then
   write(6,*)'Invalid atomic number.'     
   write(6,*)
  else
   exit
  endif
 enddo

! Read the header of the pseudopotential file.

! Skip comments.
 n=0
 do
  read(io,'(a)',err=10,end=10)comment
  if(comment(1:1)=="#")then
   n=n+1
  else
   exit
  endif
 enddo
 rewind(io)
 do i=1,n
  read(io,*)
 enddo

 read(io,*,err=10,end=10)npotd,zion,r_asymp
 read(io,*,err=10,end=10)grid_type,ngrid_points,r0,expnt,exp_expnt

 allocate(ri(ngrid_points),psp(npotd,ngrid_points),stat=ialloc)
 if(ialloc/=0)then
  write(6,*)'Allocation error.'
  stop
 endif

 do j=1,ngrid_points
  read(io,*)ri(j),(psp(i,j),i=1,npotd)
 enddo 

 outfile=trim(periodic_table_nocap(zatom))//'_pp.data'
 inquire(file=trim(outfile),exist=file_exists)
 if(file_exists)then
  write(6,*)'Output file already exists. Please remove or rename '//&
   &trim(outfile)
  stop
 endif
 open(io,file=outfile,form='formatted',iostat=ierr)
 if(ierr/=0)then
  write(6,*)'Can''t open '//trim(outfile)//' for writing.'
  stop
 endif

! Write the CASINO pseudopotential file.
 write(io,*)'Default title (pseudo converted from CHAMP format)'
 write(io,*)'Atomic number and pseudo-charge'
 write(io,'(1x,a,f16.8)')trim(i2s(zatom)),zion
 write(io,*)'Energy units (rydberg/hartree/ev):'
 write(io,*)'hartree'
 write(io,*)'Angular momentum of local component (0=s,1=p,2=d..) for DFT &
  &and QMC'
 write(io,'(1x,a)')'0 '//trim(i2s(npotd-1)) ! arbitrary
 write(io,'(1x,a)')'NLRULE override (1) VMC/DMC (2) config gen (0 ==> &
  &input/default value)'
 write(io,*)'0 0'
 write(io,*)'Number of grid points'
 write(io,'(1x,a)')trim(i2s(ngrid_points))
 write(io,*)'R(i) in atomic units'
 do i=1,ngrid_points
  write(io,*)ri(i)
 enddo
 do i=1,npotd
  write(io,*)'r*potential (L='//trim(i2s(i-1))//') in hartree'
  do j=1,ngrid_points
   write(io,*)psp(i,j)*ri(j)
  enddo
 enddo

 write(6,*)
 write(6,*)'File '//trim(outfile)//' created.'
 write(6,*)
 write(6,*)'Note: CHAMP specifies the local component of the pseudopotential in'
 write(6,*)'the input file - which this converter does not read; it is here'
 write(6,*)'arbitrarily set to the highest angular momentum. You may wish to'
 write(6,*)'change this - choosing the shallowest channel as the local part is'
 write(6,*)'usually a good idea..'
 write(6,*)

 close(io)

 stop

10 write(*,*)'Error reading file '//trim(infile)//' (header).'
 stop


CONTAINS


 CHARACTER(12) FUNCTION i2s(n)
!------------------------------------------------------------------------!
! Convert integers to left justified strings that can be printed in the  !
! middle of a sentence without introducing large amounts of white space. !
!------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n
 INTEGER i,j
 INTEGER,PARAMETER :: ichar0=ICHAR('0')
 i2s=''
 i=abs(n)
 do j=len(i2s),1,-1
  i2s(j:j)=achar(ichar0+mod(i,10))
  i=i/10 ; if(i==0)exit
 enddo ! j
 if(n<0)then
  i2s='-'//adjustl(i2s)
 else
  i2s=adjustl(i2s)
 endif ! n<0
 END FUNCTION i2s


END PROGRAM champer

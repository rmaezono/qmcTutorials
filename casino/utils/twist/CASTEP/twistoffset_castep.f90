!----------------------------------------------------------------------------!
! TWISTOFFSET_CASTEP    NDD  14/02/08                                        !
!                                                                            !
! This program applies a random offset to the grid of k vectors held in a    !
! CASTEP .cell file, to be used in conjunction with the twistav_castep       !
! script. The offset data are put in a .cell.new file.                       !
!                                                                            !
! CHANGES                                                                    !
! =======                                                                    !
! 08.2008 NDD  New random number generator.                                  !
! 06.2009 NDD  Even newer & better random number generator.                  !
! 11.2011 NDD  Avoid gfortran compiler bug.                                  !
! 07.2013 EM/NDD Added support for reduced-periodicity systems.              !
!----------------------------------------------------------------------------!


MODULE dsp
!------------------------!
! Double precision type. !
!------------------------!
 IMPLICIT NONE
 INTEGER,PARAMETER :: dp=KIND(1.d0)
END MODULE dsp


MODULE rng
!----------------------------------------------------------------------------!
! Random number generator, using the method suggested by D.E. Knuth in       !
! Seminumerical Algorithms (vol 2 of The Art of Computer Programming).       !
! The method is based on lagged Fibonacci sequences with subtraction.        !
!----------------------------------------------------------------------------!
 USE dsp,ONLY : dp
 IMPLICIT NONE
 PRIVATE
 PUBLIC ranx
 INTEGER,PARAMETER :: KK=100,LL=37 ! Leave these.
 REAL(dp) :: ranstate(kk)  ! Determines output of gen_ran_array.


CONTAINS


 REAL(dp) FUNCTION ranx()
!------------------------------------------------------------------------------!
! Return a random number uniformly distributed in [0,1).                       !
! Uses M. Luescher's suggestion: generate 1009 random numbers at a time using  !
! Knuth's algorithm, but only use the first 100.                               !
!------------------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,PARAMETER :: seed=310952  ! Random seed, betw. 0 & 2^30-3.
  INTEGER,PARAMETER :: Nran=1009,Nkeep=100 ! See comment on p. 188 of Knuth.
  INTEGER,SAVE :: i=-1
  REAL(dp),SAVE :: ran_array(Nran)
  if(i==-1)then
   call init_rng(seed) ! Initialize the RNG.
   i=Nkeep
  endif ! First call.
  if(i==Nkeep)then
   call gen_ran_array(ran_array,Nran) ! Generate a new array of random nos.
   i=0
  endif ! i=Nkeep
  i=i+1
  ranx=ran_array(i)
 END FUNCTION ranx


 SUBROUTINE gen_ran_array(ran_array,N)
!---------------------------------------------------------------!
! Generate an array of N random numbers: see Knuth's ran_array. !
!---------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: N
  REAL(dp),INTENT(out) :: ran_array(N)
  INTEGER j
  ran_array(1:KK)=ranstate(1:KK)
  do j=KK+1,N
   ran_array(j)=mod(ran_array(j-KK)+ran_array(j-LL),1.d0)
  enddo ! j
  do j=1,LL
   ranstate(j)=mod(ran_array(N+j-KK)+ran_array(N+j-LL),1.d0)
  enddo ! j
  do j=LL+1,KK
   ranstate(j)=mod(ran_array(N+j-KK)+ranstate(j-LL),1.d0)
  enddo ! j
 END SUBROUTINE gen_ran_array


 SUBROUTINE init_rng(seed)
!--------------------------------------------!
! Initialize the RNG: see Knuth's ran_start. !
!--------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: seed
  INTEGER j,s,t,sseed
  INTEGER,PARAMETER :: MM=2**30,TT=70
  REAL(dp) ss,x(KK+KK-1)
  REAL(dp),PARAMETER :: ULP=1.d0/2.d0**52,ULP2=2.d0*ULP
  if(seed<0)then
   sseed=MM-1-mod(-1-seed,MM)
  else
   sseed=mod(seed,MM)
  endif ! seed<0
  ss=ULP2*dble(sseed+2)
  do j=1,KK
   x(j)=ss
   ss=ss+ss
   if(ss>=1.d0)ss=ss-1.d0+ULP2
  enddo ! j
  x(2)=x(2)+ULP
  s=sseed
  t=TT-1
  do
   do j=KK,2,-1
    x(j+j-1)=x(j)
    x(j+j-2)=0.d0
   enddo ! j
   do j=KK+KK-1,KK+1,-1
    x(j-(KK-LL))=mod(x(j-(KK-LL))+x(j),1.d0)
    x(j-KK)=mod(x(j-KK)+x(j),1.d0)
   enddo ! j
   if(mod(s,2)==1)then
    do j=KK,1,-1
     x(j+1)=x(j)
    enddo ! j
    x(1)=x(KK+1)
    x(LL+1)=mod(x(LL+1)+x(KK+1),1.d0)
   endif ! s odd
   if(s/=0)then
    s=s/2
   else
    t=t-1
   endif ! s/=0
   if(t<=0)exit
  enddo
  ranstate(1+KK-LL:KK)=x(1:LL)
  ranstate(1:KK-LL)=x(LL+1:KK)
  do j=1,10
   call gen_ran_array(x,KK+KK-1)
  enddo ! j
 END SUBROUTINE init_rng


END MODULE rng


MODULE utils
!--------------------------!
! Miscellaneous utilities. !
!--------------------------!
 USE dsp
 IMPLICIT NONE


CONTAINS


 SUBROUTINE errstop(msg)
!------------------------------------!
! Stop and display an error message. !
!------------------------------------!
  IMPLICIT NONE
  CHARACTER(*),INTENT(in) :: msg
  write(6,'(a)')'ERROR: '//trim(msg)
  stop
 END SUBROUTINE errstop


 SUBROUTINE capitalize(string,decapitalize_in)
!-----------------------------------------------------------------------------!
! This subroutine converts the lower-case characters in string to upper-case  !
! characters (or vice versa if the optional decapitalize arg. is set to T).   !
!-----------------------------------------------------------------------------!
  IMPLICIT NONE
  CHARACTER(*),INTENT(inout) :: string
  LOGICAL,INTENT(in),OPTIONAL :: decapitalize_in
  INTEGER :: i,ichar_string
  INTEGER,PARAMETER :: ichara=ICHAR('a'),icharz=ICHAR('z'), &
   &icharBA=ICHAR('A'),icharBZ=ICHAR('Z'),deltaichar=icharBA-ichara
  LOGICAL decapitalize
  if(present(decapitalize_in))then
   decapitalize=decapitalize_in
  else
   decapitalize=.false.
  endif ! decapitalize_in supplied
  if(decapitalize)then
   do i=1,len_trim(string) ! Upper case -> lower case.
    ichar_string=ichar(string(i:i))
    if(ichar_string>=icharBA.and.ichar_string<=icharBZ) &
     &string(i:i)=achar(ichar_string-deltaichar)
   enddo ! i
  else ! Lower case -> upper case.
   do i=1,len_trim(string)
    ichar_string=ichar(string(i:i))
    if(ichar_string>=ichara.and.ichar_string<=icharz) &
     &string(i:i)=achar(ichar_string+deltaichar)
   enddo ! i
  endif ! decapitalize
 END SUBROUTINE capitalize


END MODULE utils


PROGRAM twistoffset
!---------------------------!
! Main program starts here. !
!---------------------------!
 USE dsp,ONLY : dp
 USE rng,ONLY : ranx
 USE utils,ONLY : errstop,capitalize
 IMPLICIT NONE
 INTEGER ierr,ik,Nk,twist,i,periodicity
 INTEGER,PARAMETER :: max_Nk=6*6*6
 REAL(dp) kpoint(3,max_Nk),weight(max_Nk),koffset(3)
 CHARACTER(60) seedname
 CHARACTER(150) char150
 LOGICAL found_kpoints
 REAL(dp),PARAMETER :: tol=1.d-8

! Ask user for seedname of CASTEP file.
 write(6,*)'Seedname?'
 read(5,*,iostat=ierr)seedname
 if(ierr/=0)call errstop('Problem reading seedname.  Bug?')

! Ask user for periodicity.
 write(6,*)'Periodicity?'
 read(5,*,iostat=ierr)periodicity
 if(ierr/=0)periodicity=-1
 if(periodicity<1.or.periodicity>3)call errstop('Error reading periodicity.')

! Call random number generator an appropriate number of times, so we don't
! use the same offset at each twist.
 write(6,*)'Twist number?'
 read(5,*,iostat=ierr)twist
 if(ierr/=0)twist=-1
 if(twist<0)call errstop('Error reading twist.')

 seedname=adjustl(seedname)
 open(unit=8,file=trim(seedname)//'.cell',status='old',iostat=ierr)
 if(ierr/=0)call errstop('Cannot open '//trim(seedname)//'.cell file.')
 open(unit=9,file=trim(seedname)//'.cell.new',status='replace',iostat=ierr)
 if(ierr/=0)call errstop('Cannot open '//trim(seedname)//'.cell.replace file.')

! Find k offset.  Each component is uniform in [0,1].
 koffset=0.d0
 do i=1,twist
  if(periodicity==1)then
   koffset=(/ranx(),0.d0,0.d0/)
  elseif(periodicity==2)then
   koffset=(/ranx(),ranx(),0.d0/)
  elseif(periodicity==3)then
   koffset=(/ranx(),ranx(),ranx()/)
  endif ! periodicity
 enddo ! i

 found_kpoints=.false.
 do
  read(8,'(a)',iostat=ierr)char150
  if(ierr>0)then
   call errstop('Problem reading '//trim(seedname)//'.cell.')
  elseif(ierr<0)then
   exit
  endif ! ierr>0
  write(9,'(a)')trim(char150) ! Copy out each line.
  call capitalize(char150)
  char150=adjustl(char150)

  if(index(char150,'%BLOCK')==1.and.index(char150,'KPOINTS_LIST')>0)then
! Found the kpoint block.
   if(found_kpoints)call errstop('Two sets of k-point data?')
   found_kpoints=.true.

! Read in k points.
   Nk=0
   do
    read(8,'(a)',iostat=ierr)char150
    if(ierr/=0)call errstop('Error reading '//trim(seedname) &
      &//'.cell.replace file.')
    call capitalize(char150)
    char150=adjustl(char150)
    if(index(char150,'%ENDBLOCK')==1.and.index(char150,'KPOINTS_LIST')>0)exit
    Nk=Nk+1
    if(Nk>max_Nk)call errstop('Please increase the max_NK parameter.')
    read(char150,*,iostat=ierr)kpoint(1:3,Nk),weight(Nk)
    if(ierr/=0)call errstop('Error reading '//trim(seedname) &
      &//'.cell.replace file.')
    if(periodicity==1.or.periodicity==2)then
     if(kpoint(3,Nk)/=0.d0)call errstop('k vectors inappropriate for &
      &periodicity.')
     if(periodicity==1.and.kpoint(2,Nk)/=0.d0)call &
      &errstop('k vectors inappropriate for periodicity.')
    endif ! periodicity
   enddo

   backspace(8)

   if(Nk<1)call errstop('Need more k vectors!')
   if(any(abs(weight(1:Nk)-1.d0/dble(Nk))>tol))call errstop('The k-point &
    &weights should all be the same, and should add up to 1.')
  
! Apply offset to k points.  Put them all in [0,1].  Write out new k points.
   do ik=1,Nk
    write(9,'(4(es23.16,1x))')modulo(kpoint(1:3,ik)+koffset,1.d0),1.d0/dble(Nk)
   enddo ! ik

  endif ! kpoints_list

 enddo

 if(.not.found_kpoints)then
  write(9,'("KPOINT_MP_OFFSET ",3(es23.16,1x))')koffset
 endif

END PROGRAM twistoffset

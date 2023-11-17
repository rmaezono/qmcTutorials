PROGRAM ce_calculator
!-------------------------------------------------------------!
! Utility to obtain % of the difference between two reference !
! values achieved by a set of values.                         !
!-------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,PARAMETER :: dp=kind(1.d0)
 REAL(dp) REF0,dREF0,d2REF0,REF100,dREF100,d2REF100,X,dX,d2X,PC,dPC,&
  &a,b,c,a2,b2,c2
 CHARACTER(80) line
 INTEGER ifmt

! Print title.
 write(6,*)'============'
 write(6,*)'% calculator'
 write(6,*)'============'
 write(6,*)

 do

  do

! Read REF0.
   do
    write(6,*)'Enter 0% reference value ["q" to quit]:'
    read(*,'(a)')line
    call process_line(line,ifmt,REF0,dREF0)
    select case(ifmt)
    case(-1)
     write(6,*)'Quitting.'
     write(6,*)
     stop
    case(0,-2)
     write(6,*)'Could not read value. Try again.'
     write(6,*)
     cycle
    end select
    exit
   enddo
   write(6,*)

! Read REF100.
   do
    write(6,*)'Enter 100% reference value ["q" to quit]:'
    read(*,'(a)')line
    call process_line(line,ifmt,REF100,dREF100)
    select case(ifmt)
    case(-1)
     write(6,*)'Quitting.'
     write(6,*)
     stop
    case(0,-2)
     write(6,*)'Could not read value. Try again.'
     write(6,*)
     cycle
    end select
    exit
   enddo
   write(6,*)

! Check references make sense.
   if(REF100==REF0)then
    write(6,*)'Both references are exactly equal. Try again.'
    write(6,*)
    cycle
   endif
   exit

  enddo

! Get squares.
  d2REF0=dREF0*dREF0
  d2REF100=dREF100*dREF100
  a=REF100-REF0 ; a2=a*a

! Read values until user quits.
  do
   write(6,*)'Enter value ["q" to quit, "x" to restart]:'
   read(*,'(a)')line
   call process_line(line,ifmt,X,dX)
   select case(ifmt)
   case(-2)
    write(6,*)'Restarting.'
    write(6,*)
    exit
   case(-1)
    write(6,*)'Quitting.'
    write(6,*)
    stop
   case(0)
    write(6,*)'Could not read value. Try again.'
    write(6,*)
    cycle
   end select
! Compute percentage and print.
   d2X=dX*dX
   b=X-REF0 ; b2=b*b
   c=X-REF100 ; c2=c*c
   PC=b/a
   dPC=sqrt(a2*d2X+b2*d2REF100+c2*d2REF0)/a2
   call print_number(PC*100.d0,dPC*100.d0,ifmt)
  enddo

 enddo


CONTAINS


 SUBROUTINE process_line(line,ifmt,x,dx)
!-----------------------------------------------------!
! Given a line containing a value and an errorbar in  !
! one of several formats, return the value, errorbar  !
! and format X, DX and IFMT, respectively.            !
! Reading errors flagged with IFMT=0; IFMT is set to  !
! -1 if LINE=='q' or 'Q', and IFMT is set to -2 if    !
! LINE=='x' or 'X'.                                   !
!-----------------------------------------------------!
 CHARACTER(80),INTENT(inout) :: line
 INTEGER,INTENT(out) :: ifmt
 REAL(dp),INTENT(out) :: x,dx
 CHARACTER(80) c1
 REAL(dp) temp_x,temp_dx
 INTEGER ich1,ich2,ich3,n,ierr

 ifmt=0
 x=0.d0 ; dx=0.d0
 line=trim(adjustl(line))

! 'Q'?
 if(line=='q'.or.line=='Q')then
  ifmt=-1
  return
 endif

! 'X'?
 if(line=='x'.or.line=='X')then
  ifmt=-2
  return
 endif

! I) Remove LaTeX formatting.

! I.1) Remove '$' (math mode).
 do
  ich1=scan(line,'$')
  if(ich1<1)exit
  if(ich1==1)then
   line=line(2:)
  elseif(ich1==len(line))then
   line=line(:ich1-1)
  else
   line=line(:ich1-1)//line(ich1+1:)
  endif
 enddo

! I.2) Change '&' for '.' (table alignment).
 do
  ich1=scan(line,'&')
  if(ich1<1)exit
  if(ich1==1)then
   line=line(2:)//'.'
  elseif(ich1==len(line))then
   line='.'//line(:ich1-1)
  else
   line=line(:ich1-1)//'.'//line(ich1+1:)
  endif
 enddo

! II) Read number.

! II.1) Try to match 'E dE'.
 read(line,*,iostat=ierr)temp_x,temp_dx
 if(ierr==0)then
  x=temp_x ; dx=temp_dx
  ifmt=1
  return
 endif

! II.2) Try to match 'E +- dE', 'E +/- dE', 'E \pm dE', etc.
 read(line,*,iostat=ierr)temp_x,c1,temp_dx
 if(ierr==0)then
  x=temp_x ; dx=temp_dx
  ifmt=2
  return
 endif

! II.3) Try to match 'E' (no errorbar).
 read(line,*,iostat=ierr)temp_x
 if(ierr==0)then
  x=temp_x ; dx=0.d0
  ifmt=3
  return
 endif

! II.4) Try to match 'E(dE)' (reverse EPMD2ED formatting).
 ich2=scan(line,'(')
 if(ich2>0)then
  c1=line(1:ich2-1)
  read(c1,*,iostat=ierr)temp_x
  if(ierr==0)then
   ich1=scan(c1,'.')
   n=0 ; if(ich1>0)n=ich2-ich1-1
   ich3=scan(line,')')
   if(ich3>ich2+1)then
    c1=line(ich2+1:ich3-1)
    read(c1,*,iostat=ierr)temp_dx
    if(ierr==0)then
     x=temp_x ; dx=temp_dx*10.d0**(-n)
     ifmt=4
     return
    endif
   endif
  endif
 endif

 END SUBROUTINE process_line


 SUBROUTINE print_number(x,dx,ifmt)
!---------------------------------------------------!
! Print the number X +- DX in the specified format. !
!---------------------------------------------------!
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: x,dx
 INTEGER,INTENT(in) :: ifmt
 CHARACTER(80) c1,c2
 write(6,*)
 write(c1,*)x
 if(dx==0.d0)then
  write(6,*)'%: ',trim(adjustl(c1))
 else
  write(c2,*)dx
  select case(ifmt)
  case(1)   ; write(6,*)'%: ',trim(adjustl(c1)),' ',trim(adjustl(c2))
  case(2)   ; write(6,*)'%: ',trim(adjustl(c1)),' +/- ',trim(adjustl(c2))
  case(3,4) ; write(6,*)'%: ',trim(epmd2ed(x,dx,1))
  end select
 endif
 write(6,*)
 END SUBROUTINE print_number


 CHARACTER(80) FUNCTION epmd2ed(x,dx,nfig)
!--------------------------------------------------------------!
! Print e.g. 1.2345(6) from e.g. X=1.2345111 and dX=0.0005923. !
! NFIG determines the number of figures in the errorbar.       !
!--------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: nfig
 REAL(dp),INTENT(in) :: x,dx
 CHARACTER(80) c1,c2
 INTEGER oom,resc_dx
 REAL(dp) resc_factor,resc_x
 if(dx==0.d0)then
  if(x==0.d0)then
   epmd2ed="0"
  else
   write(c1,'(f21.12)')x ; epmd2ed=trim(adjustl(c1))
  endif
  return
 endif
 if(nfig==0)then
  write(c1,'(f21.12)')x ; write(c2,'(f21.12)')dx
  epmd2ed=trim(adjustl(c1))//' +/- '//trim(adjustl(c2))
  return
 endif
 oom=nfig-1+ceiling(-log10(dx))
 resc_factor=10.d0**dble(oom)
 resc_dx=nint(resc_factor*dx)
 if(resc_dx==10**nfig)then
  oom=oom-1
  resc_factor=resc_factor*0.1d0
  resc_dx=resc_dx/10
 endif
 if(oom<0)resc_dx=resc_dx*10**(-oom)
 resc_x=dble(nint(resc_factor*x))/resc_factor
 if(oom>0)then
  c2='(f21.'//trim(i2s(oom))//')'
  write(c1,c2)resc_x
 else
  c1=trim(i2s(int(resc_x)))
 endif
 epmd2ed=trim(adjustl(c1))//'('//trim(i2s(resc_dx))//')'
 END FUNCTION epmd2ed


 CHARACTER(20) FUNCTION i2s(n)
!----------------------------!
! Convert integer to string. !
!----------------------------!
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


END PROGRAM ce_calculator

PROGRAM pspconv
!------------------------------------------------------------------------!
! ABINIT to CASINO pseudopotential converter                             !
!                                                                        !
! Read ABINIT pseudopotential files and convert them to xx_pp.data       !
! files suitable for QMC calculations with CASINO.                       !
!                                                                        !
! Nick Hine Nov 2004                                                     !
!                                                                        !
! Changes                                                                !
! -------                                                                !
! NDMH 11.2004 - Converts Abinit type 1 TM (pspcod=1) Pseudopotentials   !
! NDMH 11.2004 - Also converts FHI (pspcod=6) Pseudopotentials           !
! NDD  07.2012 - Fixed problem with wrapping. Allowed lmax<2. Tidied.    !
!------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER ierr,pspdat,pspcod,pspxc,lmax,lloc,mmax,ii,il
 INTEGER,PARAMETER :: io=19,dp=kind(1.d0)
 REAL(dp) zatom,zion,r2well
 REAL(dp),ALLOCATABLE :: psp(:,:),ri(:)
 CHARACTER(80) title,units,infile
 CHARACTER(14) outfile

! Periodic table up to uranium by atomic number, for file names
 CHARACTER(2),PARAMETER :: periodic_table_nocap(0:92)=     &
  &(/'wg','h ','he','li','be','b ','c ','n ',              &
  & 'o ','f ','ne','na','mg','al','si','p ','s ','cl','ar',&
  & 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',&
  & 'zn','ga','ge','as','se','br','kr','rb','sr','y ','zr',&
  & 'nb','mo','tc','ru','rh','pd','ag','cd','in','sn','sb',&
  & 'te','i ','xe','cs','ba','la','ce','pr','nd','pm','sm',&
  & 'eu','gd','tb','dy','ho','er','tm','yb','lu','hf','ta',&
  & 'w ','re','os','ir','pt','au','hg','tl','pb','bi','po',&
  & 'at','rn','fr','ra','ac','th','pa','u '/)

 write(6,*)
 write(6,*)'PSPCONV: ABINIT -> CASINO Pseudopotential Converter'
 write(6,*)
 do
  write(6,*)'Enter the name of the ABINIT psp file to convert:'
  read(5,*,iostat=ierr)infile
  if(ierr==0)exit
  write(6,*)'Please try again.'
 enddo
 write(6,*)
 infile=adjustl(infile)
 write(6,*)'Reading file '//trim(infile)//'.'

 open(unit=io,file=infile,form='formatted',status='old',iostat=ierr)
 if(ierr/=0)then
  write(6,*)'Can''t open file '//trim(infile)//'.'
  stop
 endif

! Read the header of the psp file
 read(io,'(a)',err=666,end=666)title
 read(io,*,err=666,end=666)zatom,zion,pspdat
 read(io,*,err=666,end=666)pspcod,pspxc,lmax,lloc,mmax,r2well

! Call appropriate subroutine to read the type of pseudopotential found
 select case (pspcod)
  case(1)
   write(6,*)'Troullier-Martins (pspcod=1) pseudopotential.'
   call readpspcod1
  case(2)
   write(6,*)'Goedecker-Teter-Hutter (pspcod=2) pseudopotential.'
   call readpspcod2 ! not supported
  case(3)
   write(6,*)'Hartwigsen-Goedecker-Hutter (pspcod=3) pseudopotential.'
   call readpspcod3 ! not supported
  case(4)
   write(6,*)'Old format pseudopotential (pspcod=4) - not supported.'
   stop
  case(5)
   write(6,*)'Old format pseudopotential (pspcod=5) - not supported.'
   stop
  case(6)
   write(6,*)'Fritz-Haber Institut (pspcod=6) pseudopotential from fhi98pp code.'
   call readpspcod6
  case(8)
   write(6,*)'Flexible format non-standard (pspcod=8) pseudopotential - not &
    &supported.'
   stop
  case default
   write(6,*)'Unknown pseudopotential format found: pspcod=',pspcod
   write(6,*)'Only pspcods 1 and 6 are currently supported.'
   stop
 end select

 close(io)

! Find name of output file
 if(nint(zatom)>92.or.nint(zatom)<0)then
  write(6,*)'Invalid zatom (>92 or <0). zatom = '//trim(i2s(nint(zatom)))//'.'
  stop
 endif

 outfile=trim(periodic_table_nocap(nint(zatom)))//'_pp.data'
 open(io,file=outfile,form='formatted',iostat=ierr)
 if(ierr/=0)then
  write(6,*)'Can''t open '//trim(outfile)//' for writing.'
  stop
 endif

! Write the CASINO pseudopotential file
 write(io,'(a)')trim(title)
 write(io,*)'Atomic number and pseudo-charge'
 write(io,'(1x,i10,f16.8)')nint(zatom),zion
 write(io,*)'Energy units (rydberg/hartree/ev):'
 write(io,*)trim(units)
 write(io,*)'Angular momentum of local component (0=s,1=p,2=d..) for DFT and &
  &QMC'
 write(io,'(1x,i10,1x,i10)')lloc,lloc
 write(io,'(1x,a)')'NLRULE override (1) VMC/DMC (2) config gen (0 ==> &
  &input/default value)'
 write(io,*)'0 0'
 write(io,*)'Number of grid points'
 write(io,'(1x,i10)')mmax
 write(io,*)'R(i) in atomic units'
 do ii=1,mmax
  write(io,*)ri(ii)
 enddo ! ii
 do il=0,lmax
  write(io,*)'r*potential (L='//trim(i2s(il))//') in '//trim(units)
  do ii=1,mmax
   write(io,*)psp(ii,il)*ri(ii)
  enddo ! ii
 enddo ! il

 write(6,*)'File '//trim(outfile)//' created.'
 write(6,*)

 deallocate(psp,ri)

 close(io)

 stop

666 write(*,*)'Error reading file '//trim(infile)//' (header).'
 stop


CONTAINS


 SUBROUTINE readpspcod1
!-----------------------------------------------------------------------!
! Read a pspcod=1 pseudopotential                                       !
! pspcod=1 : Troullier-Martins pseudopotentials, generated by DC Allan  !
! and A Khein. See ~ABINIT/Infos/Psp_infos/psp1.info                    !
!-----------------------------------------------------------------------!
! Locals
 INTEGER ialloc,ii,il,l,nproj,lmax2
 REAL(dp) e990,e999,rcpsp,rms,ekb1,ekb2,epsatm,rchrg,fchrg,qchrg,x

 units='hartree'

 do il=0,lmax
  read(io,*,err=667,end=667)l,e990,e999,nproj,rcpsp
  read(io,*,err=667,end=667)rms,ekb1,ekb2,epsatm
 enddo ! il
 read(io,*,err=667,end=667)rchrg,fchrg,qchrg

! Added by BM to fix core dump by psp files with lmax=1
 lmax2=max(lmax,2)

! Allocate arrays
 allocate(ri(1:mmax),psp(1:mmax,0:lmax2),stat=ialloc)
 if(ialloc/=0)then
  write(6,*)'Error: Unable to allocate arrays for ri and/or psp'
  stop
 endif
 psp=0.d0

! Generate the radial grid according to formula in psp1.info
 do ii=1,mmax
  x=dble(ii-1)/dble(mmax-1)
  ri(ii)=100.d0*(x+1.d-2)**5-1.d-8
 enddo ! ii

! Read the psp data (3 values per line)
 do il=0,lmax
  read(io,*,err=667,end=667)l
  if(l/=il)write(6,*)'Unexpected value '//trim(i2s(l))//' of l found &
   &(expecting '//trim(i2s(il))//').'
  ii=1
  do
   if(ii>mmax)exit ! Check if this works for mmax not a multiple of 3
   read(io,*,err=667,end=667)psp(ii,il),psp(ii+1,il),psp(ii+2,il)
   ii=ii+3
  enddo
 enddo ! il

 return

667 write(6,*)'Error reading '//trim(infile)//' (readpspcod1).'
 stop

 END SUBROUTINE readpspcod1


 SUBROUTINE readpspcod2
!-----------------------------------------------------------------------!
! Read a pspcod=2 pseudopotential                                       !
! pspcod=2 : Goedecker-Teter-Hutter (GTH) pseudopotentials.             !
! See Phys. Rev. B 54, 1703 (1996) if needed                            !
!-----------------------------------------------------------------------!
! This would be rather hard to implement as they're not very suitable for
! transferring onto a grid.
 write(6,*)'Not implemented yet'
 stop
 END SUBROUTINE readpspcod2


 SUBROUTINE readpspcod3
!-----------------------------------------------------------------------!
! Read a pspcod=3 pseudopotential                                       !
! pspcod=3 : Hartwigsen-Goedecker-Hutter pseudopotentials.              !
! See Phys. Rev. B 58, 3641 (1998) if needed, and the file              !
! ~ABINIT/Infos/Psp_infos/psp3.info                                     !
!-----------------------------------------------------------------------!
! So would this.
 write(6,*)'Not implemented yet'
 stop
 END SUBROUTINE readpspcod3


! pspcod=4 or 5 : unsupported old format pseudopotentials
! See ~ABINIT/Infos/Psp_infos/psp45.info

 SUBROUTINE readpspcod6
!-----------------------------------------------------------------------!
! Read a pspcod=6 pseudopotential (Originally by R. Gaudoin)            !
! pspcod=3 : pseudopotentials from the fhi98pp code, see                !
! ~ABINIT/Infos/Psp_infos/psp6.info for details                         !
!-----------------------------------------------------------------------!
! Locals
 INTEGER ialloc,n,di,lmax2
 REAL(dp) f,df1,df2,df3

 units='rydberg'
 do ii=1,15
  read(io,*,err=668,end=668)
 enddo ! ii
 read(io,*,err=668,end=668)n,f
 if(mmax/=n)then
  write(6,*)'Error: mmax is not the same as n - inconsistent psp file.'
  stop
 endif ! mmax/=n

! Added by BM to fix core dump by psp files with lmax=1
 lmax2=max(lmax,2)

! Allocate arrays
 n=n+1
 mmax=n
 allocate(ri(n),psp(n,0:lmax2),stat=ialloc)
 if(ialloc/=0)then
  write(6,*)'Error: Unable to allocate arrays for ri and/or psp'
  stop
 endif ! ialloc/=0
 psp=0.d0

! Read in the radial grid
 ri(1)=0.d0
 do ii=2,n
  read(io,*,err=668,end=668)di,ri(ii),df1,df2
  psp(ii,0)=df2+df2
 enddo ! ii

! Read in the psp data
 do il=1,lmax
  psp(1,il)=0.d0
  do ii=2,n
   read(io,*,err=668,end=668)di,df2,df1,df3
   psp(ii,il)=df3+df3
  enddo ! ii
 enddo ! il

 return

668 write(*,*)'Error reading file '//trim(infile)//' (readpspcod6).'
 stop

 END SUBROUTINE readpspcod6


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


END PROGRAM pspconv

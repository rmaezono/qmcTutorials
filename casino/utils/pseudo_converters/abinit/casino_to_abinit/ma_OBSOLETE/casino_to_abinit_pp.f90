PROGRAM casino2abinit
!------------------------------------------------------------------------!
! CASINO to ABINIT pseudopotential converter                             !
!                                                                        !
! Read CASINO pseudopotential files and convert them to .pspnc           !
! files suitable for DFT calculation with ABINIT. The .pspnc files are   !
! type 1 Troullier-Martins Pseudopotentials.                             !
!                                                                        !
! AM 3.2005                                                              !
!------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER ierr,ialloc,pspcod,pspxc,lmax,lloc,mmax,mmaxa,ii,il,ip
 INTEGER,PARAMETER :: io=19,dp=kind(1.d0),npoly=8
 REAL(dp) zatom,zion,r2well,rel_error,rcpsp(0:2,0:92)
 REAL(dp),ALLOCATABLE :: psp(:,:),ri(:),pspa(:,:),ria(:)
 CHARACTER(80) title,units
 CHARACTER(80) :: infile
 CHARACTER(80) :: outfile
 CHARACTER(8) :: pspdat

! Periodic table up to uranium by atomic number, for file names
 CHARACTER(2) periodic_table_nocap(0:92)
 DATA periodic_table_nocap/ &
  & 'wg','h ','he','li','be','b ','c ','n ',               &
  & 'o ','f ','ne','na','mg','al','si','p ','s ','cl','ar',&
  & 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',&
  & 'zn','ga','ge','as','se','br','kr','rb','sr','y ','zr',&
  & 'nb','mo','tc','ru','rh','pd','ag','cd','in','sn','sb',&
  & 'te','i ','xe','cs','ba','la','ce','pr','nd','pm','sm',&
  & 'eu','gd','tb','dy','ho','er','tm','yb','lu','hf','ta',&
  & 'w ','re','os','ir','pt','au','hg','tl','pb','bi','po',&
  & 'at','rn','fr','ra','ac','th','pa','u '/

! Core radius of CASINO pseudopotentials in a.u. from summary files
  rcpsp=reshape( source=(/0.d0,0.d0,0.d0,0.5d0,0.5d0,0.5d0,&
   &0.6d0,0.6d0,0.6d0,2.19d0,2.37d0,2.37d0,1.88d0,1.96d0,1.96d0,&
   &1.41d0,1.41d0,1.41d0,1.10d0,1.10d0,1.10d0,0.94d0,0.88d0,0.84d0,&
   &0.80d0,0.75d0,0.99d0,0.70d0,0.64d0,0.89d0,0.63d0,0.57d0,0.63d0,&
   &2.70d0,2.85d0,2.85d0,2.38d0,2.38d0,2.38d0,1.94d0,2.28d0,2.28d0,&
   &1.67d0,2.01d0,2.06d0,1.48d0,1.71d0,1.71d0,1.33d0,1.50d0,1.50d0,&
   &1.19d0,1.34d0,1.34d0,1.09d0,1.20d0,1.31d0,3.97d0,3.97d0,3.97d0,&
   &3.26d0,3.26d0,3.26d0,3.02d0,2.90d0,1.01d0,2.85d0,2.73d0,0.90d0,&
   &2.72d0,2.62d0,0.82d0,2.83d0,2.47d0,0.78d0,2.47d0,2.90d0,0.70d0,&
   &2.38d0,2.28d0,0.65d0,2.29d0,2.23d0,0.61d0,2.24d0,2.15d0,0.57d0,&
   &2.34d0,2.07d0,0.55d0,2.09d0,2.48d0,0.51d0,1.84d0,2.37d0,2.82d0,&
   &1.67d0,2.06d0,2.59d0,1.51d0,1.85d0,2.38d0,1.43d0,1.68d0,1.87d0,&
   &1.32d0,1.52d0,1.94d0,1.25d0,1.42d0,1.65d0,4.25d0,4.98d0,4.98d0,&
   &3.57d0,3.97d0,2.07d0,3.25d0,3.39d0,1.67d0,3.09d9,3.17d0,1.50d0,&
   &3.10d0,3.01d0,1.41d0,2.98d0,2.86d0,1.29d0,2.69d0,3.20d0,1.19d0,&
   &2.77d0,2.63d0,1.13d0,2.71d0,2.53d0,1.07d0,2.14d0,2.45d0,1.03d0,&
   &2.56d0,2.39d0,0.97d0,2.31d0,2.75d0,0.91d0,2.06d0,2.70d0,3.25d0,&
   &1.92d0,2.34d0,2.98d0,1.81d0,2.15d0,2.77d0,1.73d0,2.00d0,2.47d0,&
   &1.65d0,1.86d0,2.56d0,1.98d0,1.98d0,2.38d0,4.74d0,4.93d0,4.93d0,&
   &4.13d0,4.53d0,2.36d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,&
   &0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,&
   &0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,&
   &0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,&
   &0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,&
   &3.01d0,3.17d0,1.72d0,2.85d0,3.00d0,1.56d0,2.70d0,2.89d0,1.46d0,&
   &2.59d0,2.63d0,1.37d0,2.52d0,3.00d0,1.30d0,2.43d0,2.59d0,1.25d0,&
   &2.33d0,2.53d0,1.20d0,2.36d0,2.43d0,1.15d0,2.30d0,2.36d0,1.11d0,&
   &2.16d0,2.60d0,1.06d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,&
   &0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,&
   &0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,&
   &0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,&
   &0.00d0,0.00d0,0.00d0/),shape=(/3,93/))

 write(6,*)'WARNING: projectors are required.'
 write(6,*)'================================='
 write(6,*)
 write(6,*)'Currently this utility only reformats the x_pp.data'
 write(6,*)'file into a format readable by ABINIT, but projectors'
 write(6,*)'are required if the pseudopotential is to be used.'

 write(6,*)
 write(6,*)'CASINO to ABINIT Pseudopotential Converter'
 write(6,*)
 write(6,*)'Enter the name of the CASINO pseudopotential file to convert:'
 read(5,*)infile
 write(6,*)
 write(6,*)'Reading file ', trim(infile)

 open(unit=io,file=infile,form='formatted',status='old',iostat=ierr)
 if(ierr/=0)then
  write(6,*)'Can''t open file ', trim(infile)
  stop
 endif

 lmax=2 ! CASINO expects l=0,1,2

! Read the header of CASINO element_pp.data file
 read(io,fmt='(80a)')title
 read(io,*)
 read(io,*)zatom,zion
 read(io,*)
 read(io,*)units
 read(io,*)
 read(io,*)lloc
 read(io,*)
 read(io,*)
 read(io,*)
 read(io,*)mmax
 read(io,*)

 ! Allocate arrays
 allocate(ri(1:mmax),psp(1:mmax,0:lmax),stat=ialloc)
 if(ialloc/=0)then
  write(6,*)'Error: Unable to allocate arrays for ri and/or psp'
  stop
 endif
 ! Read in values of radial grid ri and pseudopotential*ri
 do ii=1,mmax
  read(io,*)ri(ii)
 enddo
 do il=0,lmax
  read(io,*)
  do ii=1,mmax
   read(io,*)psp(ii,il)
   if(ii>1) psp(ii,il)=psp(ii,il)/ri(ii)
  enddo
  ! Linear extrapolation for value of pseudopotential at r=0
  psp(1,il)=psp(2,il)+(psp(3,il)-psp(2,il))/(1.d0-(ri(3)/ri(2)))
 enddo

! Call appropriate subroutine to print the type of pseudopotential

 write(6,*)'Please enter the type of ABINIT pseudopotential:'
 write(6,*)'1 = Troullier-Martins'
 write(6,*)'2 = Goedecker-Teter-Hutter (not supported)'
 write(6,*)'3 = Hartwigsen-Goedecker-Hutter (not supported)'
 write(6,*)'4 = Old format pseudopotential (not supported)'
 write(6,*)'5 = Old format pseudopotential (not supported)'
 write(6,*)'6 = Fritz-Haber Institut (not supported)'
 write(6,*)
 read(5,*)pspcod
 write(6,*)
 pspcod=1
 select case (pspcod)
  case(1)
   write(6,*)'Writing out Troullier-Martins (pspcod=1) pseudopotential'
   call writepspcod1
  case(2)
   write(6,*)'Goedecker-Teter-Hutter (pspcod=2) pseudopotential'
   call writepspcod2 ! not supported
  case(3)
   write(6,*)'Hartwigsen-Goedecker-Hutter (pspcod=3) pseudopotential'
   call writepspcod3 ! not supported
  case(4)
   write(6,*)'Old format pseudopotential (pspcod=4) - Not supported'
   stop
  case(5)
   write(6,*)'Old format pseudopotential (pspcod=5) - Not supported'
   stop
  case(6)
   write(6,*)'Fritz-Haber Institut (pspcod=6) pseudopotential from fhi98pp code'
   call writepspcod6
 end select


CONTAINS


 SUBROUTINE writepspcod1
!-----------------------------------------------------------------------!
! Write a pspcod=1 pseudopotential                                      !
! pspcod=1 : Troullier-Martins pseudopotentials, generated by DC Allan  !
! and A Khein. See ~ABINIT/Infos/Psp_infos/psp1.info                    !
!-----------------------------------------------------------------------!
! Locals
 INTEGER ialloc,ii,il,nproj
 REAL(dp) e990,e999,rms,ekb1,ekb2,epsatm,rchrg,fchrg,qchrg,x

! Find name of output file
 if(zatom>92.or.zatom<0)then
  write(6,*)'Invalid zatom (>92 or <0). zatom = ',zatom
  stop
 endif
 outfile=trim(periodic_table_nocap(nint(zatom)))//'.pspnc'
 open(io,file=outfile,form='formatted',iostat=ierr)
 if(ierr/=0)then
  write(6,*)'Can''t open ',outfile,'for writing'
  stop
 endif

! Convert units to hartree
 if(units=='rydberg')then
  psp=psp*0.5d0
 elseif(units=='ev')then
  psp=psp*0.036749309d0
 endif

! Assign today's date to pspdat
 call date_and_time(pspdat)

! Allocate arrays
 mmaxa=2001 !required by the ABINIT type 1 pseudopotential format
 allocate(ria(1:mmaxa),pspa(1:mmaxa,0:lmax),stat=ialloc)
 if(ialloc/=0)then
  write(6,*)'Error: Unable to allocate arrays for ria and/or pspa'
  stop
 endif

! Generate the radial grid according to formula in psp1.info
 do ii=1,mmaxa
  x=dble(ii-1)/dble(mmaxa-1)
  ria(ii)=100.d0*(x+1.d-2)**5-1.d-8
 enddo

! Assign values to variables required by ABINIT
 pspxc=0        ! No exchange and correlation
 r2well=0.d0    ! No prefactor of a harmonic well (to bind electrons)
 e990= 14.70d0  ! 400eV
 e999= 14.70d0  ! 400eV
 nproj=0        ! No projection functions
 rms=0.d0       ! Measure of pseudopotential quality reflecting the value of
                ! the penalty function in designing the potential
 ekb1=0.d0      ! Kleinmann-Bylander energies for each projection function
 ekb2=0.d0
 epsatm=0.d0    ! Integral Int[0 to Inf] (4*Pi*r*(r*V(r)+Zion))
 rchrg=1.d0     ! the core charge radius for additional core charge used to
                ! match the xc contribution to the hardness matrix
 fchrg=0.d0     ! the prefactor of the core charge expression
 qchrg=0.d0     ! the total (integrated) core charge

! Interpolation to map the CASINO pseudopotential on to the 2001 grid points
! required by ABINIT

 do il=0,lmax
  do ii=1,mmaxa
  call lookup(ri(1),mmax,ria(ii),ip)
  ip=min(max(ip-(npoly-1)/2,1),mmax+1-npoly)
  call interp_nev(ri(ip),psp(ip,il),npoly,ria(ii),pspa(ii,il),rel_error)
  enddo
 enddo

! Write the header of the psp file
 write(io,*)trim(title)
 write(io,'(2f12.8,a12,a25)')zatom,zion,pspdat,'zatom, zion, pspdat'
 write(io,'(a3,a5,a5,a5,a6,f12.8,a40)')trim(i2s(pspcod)),trim(i2s(pspxc)),&
  &trim(i2s(lmax)),trim(i2s(lloc)),trim(i2s(mmaxa)),r2well,&
  &'pspcod,pspxc,lmax,lloc,mmax,r2well'
 do il=0,lmax
  write(io,'(a3,2f8.3,2x,a5,f15.8,a35)')trim(i2s(il)),e990,e999,&
   &trim(i2s(nproj)),rcpsp(il,nint(zatom)),'l,e99.0,e99.9,nproj,rcpsp'
  write(io,'(4f12.8,a23)')rms,ekb1,ekb2,epsatm,'rms,ekb1,ekb2,epsatm'
 enddo
 write(io,'(3f16.12,a20)')rchrg,fchrg,qchrg,'rchrg,fchrg,qchrg'

! Write pseudopotential, 3 values per line
 do il=0,lmax
  write(io,*)' ',trim(i2s(il)),' = l for CASINO pseudopotential'
  ii=1
  do
   if(ii>mmaxa)exit
   write(io,'(3e25.16)')pspa(ii,il),pspa(ii+1,il),pspa(ii+2,il)
   ii=ii+3
  enddo
 enddo

! Write projection functions, 3 values per line
 do il=0,lmax
  write(io,*)' ',trim(i2s(il)),' = l for projection function'
  ii=1
  do
   if(ii>mmaxa)exit
   write(io,'(3e25.16)')0.d0,0.d0,0.d0
   ii=ii+3
  enddo
 enddo

 write(6,*)'done.'

 END SUBROUTINE writepspcod1


 SUBROUTINE writepspcod2
!-----------------------------------------------------------------------!
! Write a pspcod=2 pseudopotential                                      !
! pspcod=2 : Goedecker-Teter-Hutter (GTH) pseudopotentials.             !
! See Phys. Rev. B 54, 1703 (1996) if needed                            !
!-----------------------------------------------------------------------!
! This would be rather hard to implement as they're not very suitable for
! transferring from a grid.
 write(6,*)'Not implemented yet.'
 stop
 END SUBROUTINE writepspcod2


 SUBROUTINE writepspcod3
!-----------------------------------------------------------------------!
! Read a pspcod=3 pseudopotential                                       !
! pspcod=3 : Hartwigsen-Goedecker-Hutter pseudopotentials.              !
! See Phys. Rev. B 58, 3641 (1998) if needed, and the file              !
! ~ABINIT/Infos/Psp_infos/psp3.info                                     !
!-----------------------------------------------------------------------!
! So would this.
 write(6,*)'Not implemented yet.'
 stop
 END SUBROUTINE writepspcod3


! pspcod=4 or 5 : unsupported old format pseudopotentials
! See ~ABINIT/Infos/Psp_infos/psp45.info

 SUBROUTINE writepspcod6
!-----------------------------------------------------------------------!
! Read a pspcod=6 pseudopotential (Originally by R. Gaudoin)            !
! pspcod=3 : pseudopotentials from the fhi98pp code, see                !
! ~ABINIT/Infos/Psp_infos/psp6.info for details                         !
!-----------------------------------------------------------------------!
! Locals
 write(6,*)'Not implemented yet'
 stop
 END SUBROUTINE writepspcod6


 SUBROUTINE lookup(xarr,n,x,jl)
!---------------------------------------------------------------------------!
! This routine works as follows. Given an array xarr(1:n) and given a value !
! x, it returns a value jl such that x is between xarr(jl) and xarr(jl+1).  !
! xarr(1:n) must be monotonic increasing or decreasing and jl=0 or jl=n is  !
! returned to indicate that x is out of range.                              !
!---------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n
 REAL(dp),INTENT(in) :: xarr(*),x
 INTEGER,INTENT(out) :: jl
 INTEGER jm,ju
! Initialize lower and upper limits.
 jl=0 ; ju=n+1
! If not yet done, compute a midpoint and replace either the lower
! or upper limit as appropriate.
 if(xarr(n)>xarr(1))then ! Increasing
  do while(ju-jl>1)
   jm=(ju+jl)/2
   if(x>xarr(jm))then
    jl=jm
   else
    ju=jm
   endif
  enddo ! ju-jl>1
 else ! Decreasing
  do while(ju-jl>1)
   jm=(ju+jl)/2
   if(x<=xarr(jm))then
    jl=jm
   else
    ju=jm
   endif !
  enddo ! ju-jl>1
 endif ! Increasing / decreasing
 END SUBROUTINE lookup


 SUBROUTINE interp_nev(points,val,npoints,x,val_interp,err_interp)
!-------------------------------------------------------------!
! INTERP_NEV: an implementation of Neville's algorithm which  !
! constructs the unique interpolating polynomial of degree    !
! n-1 through n points. The technique is based on the Newton  !
! form of the interpolating polynomial and the recursion      !
! relation for the divided differences.                       !
!                                                             !
! Reference: http://en.wikipedia.org/wiki/Neville's_algorithm !
!                                                             !
! Supplied with (1) NPOINTS function values VAL on a set of   !
! POINTS and (2) an interpolation point X, then the routine   !
! returns an interpolated value VAL_INTERP with an estimated  !
! error ERR_INTERP.                                           !
!-------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: npoints
 REAL(dp),INTENT(in) :: points(npoints),val(npoints),x
 REAL(dp),INTENT(out) :: val_interp,err_interp
 INTEGER i,j,n
 INTEGER,PARAMETER :: npoly=10 ! npoints>npoly deliberately not checked.
 REAL(dp) c(npoly),d(npoly),sep1,sep2,sep3,cmd

 sep1=abs(x-points(1))
 n=1
 do i=1,npoints
  sep2=abs(x-points(i))
  if(sep2<sep1)then
   sep1=sep2
   n=i
  endif
  c(i)=val(i) ; d(i)=val(i)
 enddo
! n is now the index of the grid point closest to the interpolation point x
 val_interp=val(n)

! Perform Neville's algorithm
 n=n-1
 do i=1,npoints-1 ! loop over tableau columns
  do j=1,npoints-i
   sep1=points(j)-x ; sep2=points(i+j)-x
   cmd=c(j+1)-d(j)
! c and d are the differences between parents and daughters in the tableau
! as we move from left to right.
   sep3=sep1-sep2
   c(j)=sep1*cmd/sep3 ; d(j)=sep2*cmd/sep3
  enddo
  if(n*2<npoints-i)then
   err_interp=c(n+1)
  else
   err_interp=d(n)
   n=n-1
  endif
  val_interp=val_interp+err_interp ! interpolated value
 enddo

 END SUBROUTINE interp_nev



 CHARACTER(20) FUNCTION i2s(n)
!-----------------------------------------------------------------------!
! I2S                                                                   !
! ===                                                                   !
! Convert integers to left justified strings that can be printed in the !
! middle of a sentence without introducing large amounts of white space.!
!                                                                       !
! Calling routine is intended to include something like:                !
! USE utilities                                                         !
! INTEGER i                                                             !
! i=12                                                                  !
! write(6,*)'Integer number ',trim(i2s(i)),' with words at the end.'    !
!-----------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER i,j,n
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
   i2s(j:j)=achar(ichar('0')+mod(i,10))
   i=i/10
  enddo

  i=1 ; j=len_trim(i2s)
  do
   if(i>=j)exit
   tmp=i2s(j:j)
   i2s(j:j)=i2s(i:i)
   i2s(i:i)=tmp
   i=i+1
   j=j-1
  enddo

  i2s=trim(sign)//i2s

 END FUNCTION i2s


END PROGRAM casino2abinit

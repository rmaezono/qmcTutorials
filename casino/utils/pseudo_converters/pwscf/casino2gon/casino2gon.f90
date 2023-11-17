!------------------------------------------------------------------------------!
! CASINO2GON                                                                   !
! ----------                                                                   !
! JRT, 2009                                                                    !
!                                                                              !
! Code to convert psedupotentials + projectors from the CASINO website         !
! into a Kleinman-Bylander form in the format used by PWSCF (and others also). !
! Example input file would be (with no !):                                     !
!# atom,Zval,Zgrid,xmin,dx,nrgrid                                              !
!Hf     4.00000000    72.00000000    -7.00000000     0.01250000 1200           !
!0                                                                             !
!   0    2.00000000                                                            !
!   1    0.00000000                                                            !
!   2    2.00000000                                                            !
!"./hf/1/pp.data"                                                              !
!"./hf/1/awfn.data_d2s2_3F"                                                    !
!"./hf/1/awfn.data_d3p1_5G"                                                    !
!                                                                              !
!1         : Does nothing. Just a reminder of what the parameters in the next  !
!            are.                                                              !
!2         : name for pseudopotential, pseudo-Z, and 4 grid parameters for     !
!            pwscf grid.                                                       !
!3         : chosen l of local channel                                         !
!4,5,6     : occupation numbers of each orbital in GS                          !
!7         : path to ppotential input file (can bt pp_gaussian or pp.data)     !
!8,9,10,.. : files from which projectors are taken. projectors are taken from  !
!            files firstcomefirstserved, ie in the above list the d-projector  !
!            comes from awfn.data_d2s2_3F not awfn.data_d3p1_5G.               !
!------------------------------------------------------------------------------

MODULE subs

INTEGER,PARAMETER :: dp=kind(1.d0)

CONTAINS


 SUBROUTINE locate(xx,n,x,j)
!-------------------------------------------------------------!
! Given an array xx(1:n) and given a value x, this routine    !
! returns a value j such that x is between xx(j) and          !
! xx(j+1). xx(1:n) must be monotonic increasing or decreasing !
! and j=0 or j=n is returned to indicate that x is out of     !
! range.                                                      !
!-------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n
 INTEGER,INTENT(out) :: j
 REAL(dp),INTENT(in) :: xx(*),x
 INTEGER jl,jm,ju
! Initialize lower and upper limits.
 jl=0 ; ju=n+1
! If not yet done, compute a midpoint and replace either the lower
! or upper limit as appropriate.
 do while(ju-jl>1)
  jm=(ju+jl)/2
  if((xx(n)>xx(1)).eqv.(x>xx(jm)))then
   jl=jm
  else
   ju=jm
  endif
 enddo
 j=jl
 END SUBROUTINE locate


 SUBROUTINE interp_nev_with_derivs(points,val,npoints,x,val_interp,fd_interp,&
  &sd_interp)
!-------------------------------------------------------------!
! INTERP_NEV_WITH_DERIVS: an implementation of Neville's      !
! algorithm which constructs the unique interpolating         !
! polynomial of degree n-1 through n points. The technique is !
! based on the Newton form of the interpolating polynomial and!
! the recursion relation for the divided differences.         !
!                                                             !
! Reference: http://en.wikipedia.org/wiki/Neville's_algorithm !
!                                                             !
! Supplied with (1) NPOINTS function values VAL on a set of   !
! POINTS and (2) an interpolation point X, then the routine   !
! returns an interpolated value VAL_INTERP together with the  !
! first and second derivatives FD_INTERP and SD_INTERP of the !
! interpolating polynomial. The error estimates ERR_VAL,      !
! ERR_FD, and ERR_SD are not returned.                        !
!-------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: npoints
 REAL(dp),INTENT(in) :: points(npoints),val(npoints),x
 REAL(dp),INTENT(out) :: val_interp,fd_interp,sd_interp
 INTEGER i,j,n
 INTEGER,PARAMETER :: npoly=10
 REAL(dp) sep1,sep2,sep3,cmd,cmd_fd,cmd_sd,err_val,err_fd,err_sd
 REAL(dp) c(npoly),c_fd(npoly),c_sd(npoly),d(npoly),d_fd(npoly),d_sd(npoly)

 sep1=abs(x-points(1))
 n=1
 do i=1,npoints
  sep2=abs(x-points(i))
  if(sep2<sep1)then
   sep1=sep2
   n=i
  endif
  c(i)=val(i) ; c_fd(i)=0.d0 ; c_sd(i)=0.d0
  d(i)=val(i) ; d_fd(i)=0.d0 ; d_sd(i)=0.d0
 enddo
! n is now the index of the grid point closest to the interpolation point x
 val_interp=val(n) ; fd_interp=0.d0 ; sd_interp=0.d0

! Perform Neville's algorithm
 n=n-1
 do i=1,npoints-1 ! loop over tableau columns
  do j=1,npoints-i
   sep1=points(j)-x ; sep2=points(i+j)-x
   cmd=c(j+1)-d(j) ; cmd_fd=c_fd(j+1)-d_fd(j) ; cmd_sd=c_sd(j+1)-d_sd(j)
   sep3=sep1-sep2
! c and d are the differences between parents and daughters in the tableau
! as we move from left to right.
   c(j)=sep1*cmd/sep3
   c_fd(j)=(-cmd+sep1*cmd_fd)/sep3
   c_sd(j)=(-2.d0*cmd_fd+sep1*cmd_sd)/sep3
   d(j)=sep2*cmd/sep3
   d_fd(j)=(-cmd+sep2*cmd_fd)/sep3
   d_sd(j)=(-2.d0*cmd_fd+sep2*cmd_sd)/sep3
  enddo
  if(n*2<npoints-i)then
   err_val=c(n+1) ; err_fd=c_fd(n+1) ; err_sd=c_sd(n+1)
  else
   err_val=d(n) ; err_fd=d_fd(n) ; err_sd=d_sd(n)
   n=n-1
  endif
! Interpolated values
  val_interp=val_interp+err_val
  fd_interp=fd_interp+err_fd
  sd_interp=sd_interp+err_sd
 enddo

 END SUBROUTINE interp_nev_with_derivs


END MODULE subs


PROGRAM casino2gon
 USE subs
 IMPLICIT NONE
 REAL(dp) a,b,zeff,zgrid,xmin,dx
 INTEGER i,j,k,ii,jj,n1,n2,o,l,norb,nelem,nrgrid,nchan,ifail,nrgrid_wf,&
  &norb_thisfile,npoly,lmax,lloc,inp_type,l_index,nrgrid_cas
 INTEGER,ALLOCATABLE :: pwr(:,:),flags(:)
 REAL(dp) :: t1,t2,t3,fa(50),vsd,vpd,occ(3)
 REAL(dp),ALLOCATABLE :: r(:),rab(:),vs(:),vp(:),vd(:),r_cas(:),&
  &vs_cas(:),vp_cas(:),vd_cas(:),r_wf(:),orb_wf(:),orb(:)
 REAL(dp),ALLOCATABLE :: bcoef(:,:),acoef(:,:)
 CHARACTER(2) :: atom
 CHARACTER(6) :: string
 CHARACTER(30) :: filename,filename2

 npoly=6                ! Order of interpolation polynomial
 o=100                  ! Unit number for .gon style output file
 nchan=3                ! Number of channels present in input ppot
 lmax=2                 ! Biggest l value of channels
 occ(:)=0.d0

 open(o,file='ppot.gon')

! Open an input file to identify what i'm doing, and with what.
! atom,zef=(atom name,valence Z), zgrid,xmin,dx,nrgrid=(grid parameters),
! filename=(pp input file)
 open(20,file='./awfn.inp')
 read(20,*)
 read(20,*)atom,zeff,zgrid,xmin,dx,nrgrid
 read(20,*)lloc
 read(20,*)i,occ(1)
 read(20,*)i,occ(2)
 read(20,*)i,occ(3)
 read(20,*)filename
 close(20)

!------------------------------------------------------------------------------
! atom     Name used within *.gon output file
! zeff     Effective Z of pseudopotom
! zgrid    grid parameter for *.gon file
! xmin     grid parameter for *.gon file
! dx       grid parameter for *.gon file
! nrgrid   grid parameter for *.gon file
! lloc     local channel used by pwscf
! i,occ    channel number and its occupation in the ground states
! filename filename of casino format input file
! following files are the wavefunction files used for projectors, selected
! firstcomefirstserved
!------------------------------------------------------------------------------

 allocate(r(nrgrid),rab(nrgrid),vs(nrgrid),vp(nrgrid),vd(nrgrid),orb(nrgrid))
 allocate(flags(nchan))

! Need norb (number of orbitals provided) here, so get it by looking at the
! file list..
 flags(:)=0 ; ifail=0 ; norb=0
 do l_index=0,lmax ! Search for 1st of l=0, write it , then l=1 etc.
  open(20,file='./awfn.inp')
  do i=1,7 ; read(20,*) ; enddo ! Skip first seven lines, to get awfn filenames.
  do
   read(20,*,iostat=ifail)filename2
   if(ifail<0)exit

! Get no. of orbitals in current file.
   open(10,file=filename2)
   do ; read(10,*) string ; if(string=='Total')exit ; enddo
   read(10,*) norb_thisfile

! Read them in, and count.
   do ; read(10,*) string ; if(string=='Radial')exit ; enddo
   read(10,*) nrgrid_wf
   do i=1,nrgrid_wf ; read(10,*) ; enddo ! skip over grid data

   do k=1,norb_thisfile ! loop over orbitals in file
    read(10,*) ; read(10,*) i,j,l         ! Orbital data
    do i=1,nrgrid_wf ; read(10,*) ; enddo ! Skip over orbitals themselves
    if(flags(l+1)==0.and.l_index==l)then  ! Count this orbital if first one
     flags(l+1)=1 ; norb=norb+1
    endif
   enddo
   close(10)

  enddo
  close(20)
 enddo
! norb is now the number of orbitals available, firstcomefirstserved.
! flags(l+1)=1 if orbitals is available, =0 if it is not.

 vs(:)=0.d0 ; vp(:)=0.d0 ; vd(:)=0.d0
 b=dx ; a=exp(xmin)/zgrid
 do i=1,nrgrid
  r(i)=a*(exp(b*dble(i-1))) ; rab(i)=a*b*exp(b*dble(i-1)) ! PWSCF grid
 enddo
 write(*,*) 'Largest grid point, in a.u.    :',r(nrgrid)

! What type is the input file? Flag is inp_type.
 ii=len_trim(filename);inp_type=0
 do i=1,ii-10 ; if(filename(i:i+10)=='pp_gaussian')inp_type=1 ; enddo
 do i=1,ii-6 ; if(filename(i:i+6)=='pp.data')inp_type=2 ; enddo
 if(inp_type==0)then
  write(*,*)'Input ppot type not recognized from filename.' ; stop
 endif
 if(inp_type==1)write(*,*)'Input ppot is Gaussian parameterization.'
 if(inp_type==2)write(*,*)'Input ppot is CASINO tabulation.'

 if(inp_type==1)then ! If Gaussian input, do the following.

   open(10,file=filename)
   do i=1,7 ; read(10,*) ; enddo
   read(10,*)nelem
   allocate(pwr(nchan,nelem),bcoef(nchan,nelem),acoef(nchan,nelem))

   do ii=1,8
    read(10,'(1x,i1,6x,f16.8,6x,f16.8)')pwr(3,ii),bcoef(3,ii),acoef(3,ii)
   enddo
   read(10,*) ; read(10,*)
   do ii=1,8
    read(10,'(1x,i1,6x,f16.8,6x,f16.8)')pwr(1,ii),bcoef(1,ii),acoef(1,ii)
   enddo
   read(10,*) ; read(10,*)
   do ii=1,8
    read(10,'(1x,i1,6x,f16.8,6x,f16.8)')pwr(2,ii),bcoef(2,ii),acoef(2,ii)
   enddo
   close(10)

   do i=1,nrgrid
    vsd=0.d0 ; vpd=0.d0 ; vd(i)=-zeff*r(i)
    do ii=1,8
     vd(i)=vd(i)+acoef(3,ii)*r(i)**pwr(3,ii)*dexp(-bcoef(3,ii)*r(i)*r(i))
     vsd=vsd+acoef(1,ii)*r(i)**pwr(1,ii)*dexp(-bcoef(1,ii)*r(i)*r(i))
     vpd=vpd+acoef(2,ii)*r(i)**pwr(2,ii)*dexp(-bcoef(2,ii)*r(i)*r(i))
    enddo
    vs(i)=(vsd+vd(i))/r(i)**2
    vp(i)=(vpd+vd(i))/r(i)**2
    vd(i)=vd(i)/r(i)**2
   enddo

   vs=vs*2.d0 ; vp=vp*2.d0 ; vd=vd*2.d0 ! Hartrees --> rydbergs

 elseif(inp_type==2)then ! If CASINO input, do the following.

! Note the following dumps r=0 value, since it is meaningless.
   open(10,file=filename)
   do i=1,10 ; read(10,*) ; enddo
   read(10,*)nrgrid_cas ; read(10,*)
   nrgrid_cas=nrgrid_cas-1
   allocate(r_cas(nrgrid_cas),vs_cas(nrgrid_cas),vp_cas(nrgrid_cas),&
    &vd_cas(nrgrid_cas))
   read(10,*) ; do i=1,nrgrid_cas ; read(10,*)  r_cas(i) ; enddo ; read(10,*)
   read(10,*) ; do i=1,nrgrid_cas ; read(10,*) vs_cas(i) ; enddo ; read(10,*)
   read(10,*) ; do i=1,nrgrid_cas ; read(10,*) vp_cas(i) ; enddo ; read(10,*)
   read(10,*) ; do i=1,nrgrid_cas ; read(10,*) vd_cas(i) ; enddo
   close(10)

! Now interpolate to new grid
   do i=1,nrgrid
    if(r(i)<=r_cas(nrgrid_cas).and.r(i)>=r_cas(1))then
     call locate(r_cas(1),nrgrid_cas,r(i),jj)
     jj=min(max(jj-(npoly-1)/2,1),nrgrid_cas+1-npoly)
     call interp_nev_with_derivs(r_cas(jj),vs_cas(jj),npoly,r(i),t1,t2,t3)
     vs(i)=t1
     call interp_nev_with_derivs(r_cas(jj),vp_cas(jj),npoly,r(i),t1,t2,t3)
     vp(i)=t1
     call interp_nev_with_derivs(r_cas(jj),vd_cas(jj),npoly,r(i),t1,t2,t3)
     vd(i)=t1
     vs(i)=vs(i)/r(i)
     vp(i)=vp(i)/r(i)
     vd(i)=vd(i)/r(i)
    else
     write(*,*)'Help! - interpolation failed.'
     stop
    endif
   enddo

 endif ! Done decisions on input type.

! Now vs,vp,vd contains the potentials channels, in rydbergs, on pwscf grid.
! Next output formated for pwscf input in old format  - see
! "espresso-4.0/atomic_doc/pseudo-test/S.gon" and/or
! "http://www.pwscf.org/old/NC_format.txt" .

! First write out ppot channels.

 write(o,'(a3)')'TN '
 write(o,'("''",a2,"''",f8.4,1x,i4,1x,i4,1x,i4,3x,a1,1x,i4,3x,a1,1x,f12.9)')&
  &atom,zeff,lmax,0,0,"F",lloc,"F",0.d0
 write(o,'(1x,f14.10,1x,f14.10,1x,f14.10,1x,i4,i4)') zgrid,xmin,dx,nrgrid,norb

! s,p, then d channel. Actual value, in hartrees, no prefactor of r or anything.
 write(o,*)
 n1=nrgrid/4+1 ; n2=4
 do i=0,n1-1
  if((nrgrid-4*i)<4) n2=nrgrid-4*i
  if(n2/=0)write(o,'(4(1x,e18.11))')(vs(4*i+j+1),j=0,n2-1)
 enddo

 write(o,*)
 n1=nrgrid/4+1 ; n2=4
 do i=0,n1-1
  if((nrgrid-4*i)<4) n2=nrgrid-4*i
  if(n2/=0)write(o,'(4(1x,e18.11))')(vp(4*i+j+1),j=0,n2-1)
 enddo

 write(o,*)
 n1=nrgrid/4+1 ; n2=4
 do i=0,n1-1
  if((nrgrid-4*i)<4)n2=nrgrid-4*i
  if(n2/=0)write(o,'(4(1x,e18.11))')(vd(4*i+j+1),j=0,n2-1)
 enddo

! Next, read in orbitals and interpolate them to the PWSCF grid.

! Loop through orbital files listed in awfn.inp.
! If >1 orbitals exist for a given l, use the first one only.
 flags(:)=0 ; ifail=0

 do l_index=0,lmax ! Search for 1st of l=0, write it, then l=1 etc.

  open(20,file='./awfn.inp')
  do i=1,7 ; read(20,*) ; enddo  ! Skip first seven lines, to get awfn filenames.
  do
   read(20,*,iostat=ifail)filename
   if(ifail<0)exit
   write(*,*) filename

! Get number of orbitals in current file.
   open(10,file=filename)
   do ; read(10,*) string ; if(string=='Total')exit ; enddo
   read(10,*)norb_thisfile

! Read them in, and interpolate to PWSCF grid.
   do ; read(10,*) string ; if(string=='Radial')exit ; enddo
   read(10,*)nrgrid_wf
   allocate(r_wf(nrgrid_wf),orb_wf(nrgrid_wf))
   do i=1,nrgrid_wf
    read(10,*)r_wf(i)
   enddo

   do k=1,norb_thisfile

    read(10,*) ; read(10,*) i,j,l
    do i=1,nrgrid_wf ; read(10,*) orb_wf(i) ; enddo
    do i=1,nrgrid
     if(r(i)<=r_wf(nrgrid_wf))then
      call locate(r_wf(1),nrgrid_wf,r(i),jj)
      jj=min(max(jj-(npoly-1)/2,1),nrgrid_wf+1-npoly)

      if(jj/=1)then
       call interp_nev_with_derivs(r_wf(jj),orb_wf(jj),npoly,r(i),t1,t2,t3)
       orb(i)=t1
      else
       fa(2:npoly-l+1)=orb_wf(2:npoly-l+1)/r_wf(2:npoly-l+1)**l
       call interp_nev_with_derivs(r_wf(2),fa(2),npoly-l,r(i),t1,t2,t3)
       orb(i)=t1*r(i)**l
      endif
     else
      orb(i)=0.d0
     endif
     if(orb(i)<0.d0)orb(i)=0.d0
    enddo
! Current r*(radial part of wf) now in orb, so now we write it out.
! Only use the first l valued orbital in ordered list for projector (one
! allowed per channel).
    if(flags(l+1)==0.and.l_index==l)then
     write(o,*)
     write(o,'(9x,i4,2x,f18.15)')l,occ(l+1)
     n1=nrgrid/4+1 ; n2=4
     do i=0,n1-1
      if((nrgrid-4*i)<4)n2=nrgrid-4*i
      if(n2/=0)write(o,'(4(1x,e18.11))')(orb(4*i+j+1),j=0,n2-1)
     enddo
     flags(l+1)=1
    endif

   enddo
   deallocate(r_wf,orb_wf)
   close(10)

  enddo
  close(20)
 enddo ! Now look for next l_index value...

 write(*,*) 'Finished.'

 close(o)
 stop

END PROGRAM casino2gon

 SUBROUTINE blipk
!------------------------------------------------------------------!
! Read complex PW orbitals at multiple k points from pwfn.data and !
! generate complex blip orbitals.                                  !
!------------------------------------------------------------------!
  USE rng,ONLY : ranx,init_rng
  USE singleton,ONLY : fftn
  IMPLICIT NONE
  INTEGER i,j,k,ibnd,ig,lx,ly,lz,ik,spin,icount,ialloc,maxband,maxband_new
  REAL(dp) r(3),xb(5),xd(5),xbd(5),et,dot_prod,tempr, &
   &norm2_temp,overlap(5),av_overlap(5),av_overlap2(5),modc2,gtemp(3)
  REAL(dp),POINTER :: ke_pworb(:,:,:),ke_bliporb(:,:,:)
  REAL(dp),POINTER :: ke_pworb_new(:,:,:),ke_bliporb_new(:,:,:)
  COMPLEX(dp) :: d,d2,dg(3),blipval,bliplap,blipgrad(3),eigr
  COMPLEX(dp),ALLOCATABLE :: evc(:),fftdata(:)
  CHARACTER(80) char80
  INTERFACE
   REAL(kind=kind(1.d0)) FUNCTION c_norm2_integrand_zlv(zlv)
    USE helpers,ONLY : dp
    REAL(dp),INTENT(in) :: zlv
   END FUNCTION c_norm2_integrand_zlv
   REAL(kind=kind(1.d0)) FUNCTION c_ke_integrand_zlv(zlv)
    USE helpers,ONLY : dp
    REAL(dp),INTENT(in) :: zlv
   END FUNCTION c_ke_integrand_zlv
  END INTERFACE

  if(Nkvec/=1)then
   write(6,*)'Generating complex orbitals at multiple k points.'
  else
   write(6,*)'Generating complex orbitals at a single k point.'
  endif ! Nkvec/=1
  write(6,*)

  allocate(fftdata(nr(1)*nr(2)*nr(3)),evc(NGvec), &
   &cavc(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1),stat=ialloc)
  if(ialloc/=0)call errstop('BLIPK','Allocation problem: blip & PW &
   &coefficients, etc.')

  maxband=0
  nullify(ke_pworb,ke_bliporb)

  do ik=1,Nkvec
   if(ik>1)then
    ! The very first block header was read and writting by bwfn_transform
    ! because it was needed for the decision between blipp and blipk.
    ! All further headers need to be read and written here
    call skipio(io_pwfn,1)
    write(io_bwfn,*)'k-point # ; # of bands (up spin/down spin) ; k-point &
     &coords (au)'
    read(io_pwfn,*,err=10,end=20)i,nband(ik,1),nband(ik,2),kvec(ik,:)
    write(io_bwfn,'(1x,a,3(1x,es24.17))')trim(i2s(ik))//" " &
     &//trim(i2s(nband(ik,1)))//" "//trim(i2s(nband(ik,2))),kvec(ik,:)
   endif

   if(norm2_ke_calc)then
    maxband_new=maxval(nband(:ik,:))
    if(maxband_new>maxband)then
     allocate(ke_pworb_new(maxband_new,Nkvec,no_spin_types), &
      &ke_bliporb_new(maxband_new,Nkvec,no_spin_types),stat=ialloc)
     if(ialloc/=0)call errstop('BLIPK','Allocation error: KE arrays.')
     if(maxband>0)then
      ke_pworb_new(:maxband,:,:)=ke_pworb(:,:,:)
      ke_bliporb_new(:maxband,:,:)=ke_bliporb(:,:,:)
      deallocate(ke_pworb,ke_bliporb)
     endif ! maxband>0
     maxband=maxband_new
     ke_pworb=>ke_pworb_new
     ke_bliporb=>ke_bliporb_new
    endif ! maxband_new>maxband
   endif ! norm2_ke_calc

   do spin=1,no_spin_types
    write(6,'('' k-point '',i5,"  :  (",f10.6,",",f10.6,",",f10.6,")")')&
     &ik,kvec(ik,:)
    if(n_points_for_test>0)then
     if(spin_polarized)write(6,*)'Spin : '//trim(i2s(spin))
     write(6,*)'Band ; alpha ;      alpha (Lap.) ; alpha (gradient)'
    endif ! Overlap test

    do ibnd=1,nband(ik,spin)
     call skipio(io_pwfn,1)
     write(io_bwfn,*)'Band, spin, eigenvalue (au), localized'
     read(io_pwfn,*,err=10,end=20)i,i,et
     write(io_bwfn,'(1x,a,es24.17,a)')trim(i2s(ibnd))//" "//trim(i2s(spin)) &
      &//" ",et," F"
     call skipio(io_pwfn,1)

! Read in the PW coefficients
     do ig=1,NGvec
      read(io_pwfn,'(a)',err=10,end=20)char80
      if(index(char80,",")/=0)then
       read(char80,*)evc(ig)
      else
       read(char80,*)tempr
       evc(ig)=cmplx(tempr,0.d0,dp)
      endif ! Complex coefficients
! Multiply coefficients by exp(-iG.R_offset), to translate the orbitals.
      if(intended_periodicity<3)evc(ig)=evc(ig)*expmigdotR_offset(ig)
     enddo ! ig

     if(norm2_ke_calc)then
! Evaluate the KE and norm2 of the PW orbitals.
      norm2_temp=0.d0
      ke_pworb(ibnd,ik,spin)=0.d0
      do ig=1,NGvec
       gtemp=twopi*matmul(bg,dble(g(:,ig)))
       modc2=dble(evc(ig))**2+aimag(evc(ig))**2
       norm2_temp=norm2_temp+modc2
       ke_pworb(ibnd,ik,spin)=ke_pworb(ibnd,ik,spin)+&
        &sum((gtemp+kvec(ik,:))**2)*modc2
      enddo ! ig
      ke_pworb(ibnd,ik,spin)=0.5d0*ke_pworb(ibnd,ik,spin)/norm2_temp
     endif ! norm2_ke_calc

     fftdata=czero
     do ig=1,NGvec
      icount=1+g(1,ig)+nr(1)*(indexfn(g(1,ig))+g(2,ig) &
       &+nr(2)*(indexfn(g(2,ig))+g(3,ig)+nr(3)*indexfn(g(3,ig))))
      fftdata(icount)=evc(ig)*gamma(ig)
     enddo ! ig

! Carry out FFT to real space
     call fftn(fftdata,nr,inv=.true.)
     fftdata=fftdata*sqrt(dble(nr(1)*nr(2)*nr(3)))

     icount=1
     do lz=0,nr(3)-1
      do ly=0,nr(2)-1
       do lx=0,nr(1)-1
        cavc(lx,ly,lz)=fftdata(icount)
        icount=icount+1
       enddo ! lx
      enddo ! ly
     enddo ! lz

! Write out blip coefficients for extended orbital.
     write(io_bwfn,*)'Complex blip coefficients for extended orbital'
     do lx=0,nr(1)-1
      do ly=0,nr(2)-1
       do lz=0,nr(3)-1
        write(io_bwfn,*)cavc(lx,ly,lz)
       enddo ! lz
      enddo ! ly
     enddo ! lx

     if(norm2_ke_calc)then
! Evaluate KE by numerical integration over the entire sim cell.
! Integrate over the coordinates in terms of the lattice vectors.
! Parallelepiped is a cuboid in this coordinate system.
! Must scale integrals by the cell volume afterwards.
      rlv_min_pass=0.d0
      rlv_max_pass=1.d0
      kvec_pass=kvec(ik,:)
      kvec2_pass=sum(kvec(ik,:)**2)
      eps_simp=eps_norm2
      call adapt_simpson(c_norm2_integrand_zlv,&
       &rlv_min_pass(3),rlv_max_pass(3),norm2_temp)
      eps_simp=eps_ke
      call adapt_simpson(c_ke_integrand_zlv,&
       &rlv_min_pass(3),rlv_max_pass(3),ke_bliporb(ibnd,ik,spin))
      ke_bliporb(ibnd,ik,spin)=-0.5d0*ke_bliporb(ibnd,ik,spin)/norm2_temp
     endif ! norm2_ke_calc

! Carry out the overlap test described in the CASINO manual.
! Repeat the whole test n_overlap_tests times, to compute error bars.
     if(n_points_for_test>0)then
      call init_rng(12345678)

      av_overlap=0.d0 ; av_overlap2=0.d0
      do j=1,n_overlap_tests
       xb=0.d0 ; xd=0.d0 ; xbd=0.d0
       do i=1,n_points_for_test
        r(1)=ranx() ; r(2)=ranx() ; r(3)=ranx()
        call blip3dk(blipval,bliplap,blipgrad,r,.true.,.true.)
        d2=czero ; dg=czero
        d=czero
        do ig=1,NGvec
         dot_prod=twopi*(g(1,ig)*r(1)+g(2,ig)*r(2)+g(3,ig)*r(3))
         eigr=evc(ig)*cmplx(cos(dot_prod),sin(dot_prod),dp)
         d=d+eigr
         d2=d2-eigr*g2(ig)
         dg=dg+(twopi*eigr*iunity)*g(1:3,ig)
        enddo ! ig
        dg=matmul(bg,dg)
! Overlap of orbital values.
        xb(1)=xb(1)+dble(blipval)**2
        xd(1)=xd(1)+dble(d)**2
        xbd(1)=xbd(1)+dble(blipval)*dble(d)
! Overlap of orbital Laplacians.
        xb(2)=xb(2)+dble(bliplap)**2
        xd(2)=xd(2)+dble(d2)**2
        xbd(2)=xbd(2)+dble(bliplap)*dble(d2)
! Overlap of orbital gradients.
        xb(3:5)=xb(3:5)+dble(blipgrad(1:3))**2
        xd(3:5)=xd(3:5)+dble(dg(1:3))**2
        xbd(3:5)=xbd(3:5)+dble(blipgrad(1:3))*dble(dg(1:3))
       enddo ! i
       do k=1,5
        if(xb(k)/=0.d0.and.xd(k)/=0.d0)then
         overlap(k)=xbd(k)**2/(xb(k)*xd(k))
        else
         overlap(k)=0.d0
        endif ! xb & xd nonzero
       enddo ! k
       av_overlap(1:5)=av_overlap(1:5)+overlap(1:5)
       av_overlap2(1:5)=av_overlap2(1:5)+overlap(1:5)**2
      enddo ! j
      av_overlap=av_overlap/dble(n_overlap_tests)
      av_overlap2=av_overlap2/dble(n_overlap_tests)
      call write_overlaps(ibnd,av_overlap,av_overlap2)

     endif ! Overlap test

    enddo ! ibnd

    if(n_points_for_test>0)write(6,*)

   enddo ! spin
  enddo ! ik

! Write out the kinetic energies of the original PW orbitals and the blip
! orbitals.
  if(norm2_ke_calc)then
   write(6,*)'Kinetic energies of orbitals'
   if(spin_polarized)then
    write(6,*)'Spin ;  k point ; band ;   KE of PW orbs (au) ;   &
     &KE of blip orbs (au)'
    do spin=1,no_spin_types
     do ik=1,nkvec
      do ibnd=1,nband(ik,spin)
       write(6,'(3x,i1,2x,i5,2x,i6,3x,f12.8,7x,f12.8)')spin,ik,ibnd, &
        &ke_pworb(ibnd,ik,spin),ke_bliporb(ibnd,ik,spin)
      enddo ! ibnd
     enddo ! ik
    enddo ! spin
   else
    write(6,*)'    k point ;  band ;   KE of PW orbs (au) ;   &
     &KE of blip orbs (au)'
    do ik=1,nkvec
     do ibnd=1,nband(ik,1)
      write(6,'(3x,i5,5x,i6,4x,f12.8,11x,f12.8)')ik,ibnd,ke_pworb(ibnd,ik,1), &
       &ke_bliporb(ibnd,ik,1)
     enddo ! ibnd
    enddo ! ik
   endif ! spin_polarized
   write(6,*)
  endif ! norm2_ke_calc

  deallocate(gamma,fftdata,g,g2,evc,nband)
  deallocate(cavc)
  if(associated(ke_pworb))deallocate(ke_pworb)
  if(associated(ke_bliporb))deallocate(ke_bliporb)

  return

10 call errstop('BLIPK','Error reading pwfn.data.')
20 call errstop('BLIPK','Unexpectedly reached end of pwfn.data.')

 END SUBROUTINE blipk


 SUBROUTINE blip3dk(rpsi,lap,grad,r,lw,lgl)
!----------------------------------------------------------------------------!
! This subroutine evaluates the value of a function, its gradient and its    !
! Laplacian at a vector point r, using the overlapping of blip functions.    !
! The blip grid is defined on a cubic cell, so r should always be given in   !
! units of the crystal lattice vectors.                                      !
!----------------------------------------------------------------------------!
  IMPLICIT NONE
  REAL(dp),INTENT(in) :: r(3)
  LOGICAL,INTENT(in) :: lw,lgl
  COMPLEX(dp),INTENT(out) :: rpsi,lap,grad(3)
  INTEGER ix,iy,iz,ixp,ixpp,ixm,iyp,iypp,iym,izp,izpp,izm
  REAL(dp) t(3),txm,tx,txp,txpp,tym,ty, &
   &typ,typp,tzm,tz,tzp,tzpp,dtxm,dtx,dtxp,dtxpp,dtym,dty,dtyp,dtypp,dtzm,&
   &dtz,dtzp,dtzpp,d2txm,d2tx,d2txp,d2txpp,d2tym,d2ty,d2typ,d2typp,d2tzm,&
   &d2tz,d2tzp,d2tzpp,x,y,z
  COMPLEX(dp) :: d1(16),t1,t2,t3,t4,dt1,dt2,dt3,dt4,d2t1,d2t2,d2t3,d2t4, &
   &templap(6)

  rpsi=czero ; lap=czero ; grad=czero

  ix=modulo(floor(r(1)*rnr(1)),nr(1))
  iy=modulo(floor(r(2)*rnr(2)),nr(2))
  iz=modulo(floor(r(3)*rnr(3)),nr(3))

! The blips are defined as the product of one-dimensional cubic splines, these
! are different from zero only on four adjacent grid points. It follows that
! for any vector r the are only 64 overlapping three-dimensional blips, which
! are the product of all possible combinations of the three one-dimensional
! splines.

! These are the extra 3 coefficients for each dimension needed
  ixp=modulo(ix+1,nr(1))
  ixpp=modulo(ix+2,nr(1))
  ixm=modulo(ix-1,nr(1))
  iyp=modulo(iy+1,nr(2))
  iypp=modulo(iy+2,nr(2))
  iym=modulo(iy-1,nr(2))
  izp=modulo(iz+1,nr(3))
  izpp=modulo(iz+2,nr(3))
  izm=modulo(iz-1,nr(3))

! Now calculate the 12 monodimensional blip functions
  t=modulo(r*rnr,rnr)

  if(.not.lgl)then
! Just need values.
   x=t(1)-dble(ix-1)
   txm=2.d0+x*(-3.d0+x*(1.5d0-0.25d0*x))
   x=t(1)-dble(ix)
   tx=1.d0+x*x*(-1.5d0+0.75d0*x)
   x=t(1)-dble(ix+1)
   txp=1.d0+x*x*(-1.5d0-0.75d0*x)
   x=t(1)-dble(ix+2)
   txpp=2.d0+x*(3.d0+x*(1.5d0+0.25d0*x))
   y=t(2)-dble(iy-1)
   tym=2.d0+y*(-3.d0+y*(1.5d0-0.25d0*y))
   y=t(2)-dble(iy)
   ty=1.d0+y*y*(-1.5d0+0.75d0*y)
   y=t(2)-dble(iy+1)
   typ=1.d0+y*y*(-1.5d0-0.75d0*y)
   y=t(2)-dble(iy+2)
   typp=2.d0+y*(3.d0+y*(1.5d0+0.25d0*y))
   z=t(3)-dble(iz-1)
   tzm=2.d0+z*(-3.d0+z*(1.5d0-0.25d0*z))
   z=t(3)-dble(iz)
   tz=1.d0+z*z*(-1.5d0+0.75d0*z)
   z=t(3)-dble(iz+1)
   tzp=1.d0+z*z*(-1.5d0-0.75d0*z)
   z=t(3)-dble(iz+2)
   tzpp=2.d0+z*(3.d0+z*(1.5d0+0.25d0*z))
  else
! Need first and second derivatives.
   x=t(1)-dble(ix-1)
   txm=2.d0+x*(-3.d0+x*(1.5d0-0.25d0*x))
   dtxm=(-3.d0+x*(3.d0-0.75d0*x))*rnr(1)
   d2txm=(3.d0-1.5d0*x)*rnr2(1)
   x=t(1)-dble(ix)
   tx=1.d0+x*x*(-1.5d0+0.75d0*x)
   dtx=x*(-3.d0+2.25d0*x)*rnr(1)
   d2tx=(-3.d0+4.5d0*x)*rnr2(1)
   x=t(1)-dble(ix+1)
   txp=1.d0+x*x*(-1.5d0-0.75d0*x)
   dtxp=x*(-3.d0-2.25d0*x)*rnr(1)
   d2txp=(-3.d0-4.5d0*x)*rnr2(1)
   x=t(1)-dble(ix+2)
   txpp=2.d0+x*(3.d0+x*(1.5d0+0.25d0*x))
   dtxpp=(3.d0+x*(3.d0+0.75d0*x))*rnr(1)
   d2txpp=(3.d0+1.5d0*x)*rnr2(1)
   y=t(2)-dble(iy-1)
   tym=2.d0+y*(-3.d0+y*(1.5d0-0.25d0*y))
   dtym=(-3.d0+y*(3.d0-0.75d0*y))*rnr(2)
   d2tym=(3.d0-1.5d0*y)*rnr2(2)
   y=t(2)-dble(iy)
   ty=1.d0+y*y*(-1.5d0+0.75d0*y)
   dty=y*(-3.d0+2.25d0*y)*rnr(2)
   d2ty=(-3.d0+4.5d0*y)*rnr2(2)
   y=t(2)-dble(iy+1)
   typ=1.d0+y*y*(-1.5d0-0.75d0*y)
   dtyp=y*(-3.d0-2.25d0*y)*rnr(2)
   d2typ=(-3.d0-4.5d0*y)*rnr2(2)
   y=t(2)-dble(iy+2)
   typp=2.d0+y*(3.d0+y*(1.5d0+0.25d0*y))
   dtypp=(3.d0+y*(3.d0+0.75d0*y))*rnr(2)
   d2typp=(3.d0+1.5d0*y)*rnr2(2)
   z=t(3)-dble(iz-1)
   tzm=2.d0+z*(-3.d0+z*(1.5d0-0.25d0*z))
   dtzm=(-3.d0+z*(3.d0-0.75d0*z))*rnr(3)
   d2tzm=(3.d0-1.5d0*z)*rnr2(3)
   z=t(3)-dble(iz)
   tz=1.d0+z*z*(-1.5d0+0.75d0*z)
   dtz=z*(-3.d0+2.25d0*z)*rnr(3)
   d2tz=(-3.d0+4.5d0*z)*rnr2(3)
   z=t(3)-dble(iz+1)
   tzp=1.d0+z*z*(-1.5d0-0.75d0*z)
   dtzp=z*(-3.d0-2.25d0*z)*rnr(3)
   d2tzp=(-3.d0-4.5d0*z)*rnr2(3)
   z=t(3)-dble(iz+2)
   tzpp=2.d0+z*(3.d0+z*(1.5d0+0.25d0*z))
   dtzpp=(3.d0+z*(3.d0+0.75d0*z))*rnr(3)
   d2tzpp=(3.d0+1.5d0*z)*rnr2(3)
  endif ! not lgl

  d1(1)=cavc(ix,iy,iz)*tz+cavc(ix,iy,izp)*tzp+cavc(ix,iy,izpp)*tzpp+ &
   &cavc(ix,iy,izm)*tzm
  d1(2)=cavc(ix,iyp,iz)*tz+cavc(ix,iyp,izp)*tzp+cavc(ix,iyp,izpp)*tzpp+ &
   &cavc(ix,iyp,izm)*tzm
  d1(3)=cavc(ix,iypp,iz)*tz+cavc(ix,iypp,izp)*tzp+cavc(ix,iypp,izpp)*tzpp+ &
   &cavc(ix,iypp,izm)*tzm
  d1(4)=cavc(ix,iym,iz)*tz+cavc(ix,iym,izp)*tzp+cavc(ix,iym,izpp)*tzpp+ &
   &cavc(ix,iym,izm)*tzm

  d1(5)=cavc(ixp,iy,iz)*tz+cavc(ixp,iy,izp)*tzp+cavc(ixp,iy,izpp)*tzpp+ &
   &cavc(ixp,iy,izm)*tzm
  d1(6)=cavc(ixp,iyp,iz)*tz+cavc(ixp,iyp,izp)*tzp+cavc(ixp,iyp,izpp)*tzpp+ &
   &cavc(ixp,iyp,izm)*tzm
  d1(7)=cavc(ixp,iypp,iz)*tz+cavc(ixp,iypp,izp)*tzp+cavc(ixp,iypp,izpp)*tzpp+ &
   &cavc(ixp,iypp,izm)*tzm
  d1(8)=cavc(ixp,iym,iz)*tz+cavc(ixp,iym,izp)*tzp+cavc(ixp,iym,izpp)*tzpp+ &
   &cavc(ixp,iym,izm)*tzm

  d1(9)=cavc(ixpp,iy,iz)*tz+cavc(ixpp,iy,izp)*tzp+cavc(ixpp,iy,izpp)*tzpp+ &
   &cavc(ixpp,iy,izm)*tzm
  d1(10)=cavc(ixpp,iyp,iz)*tz+cavc(ixpp,iyp,izp)*tzp &
   &+cavc(ixpp,iyp,izpp)*tzpp+cavc(ixpp,iyp,izm)*tzm
  d1(11)=cavc(ixpp,iypp,iz)*tz+cavc(ixpp,iypp,izp)*tzp &
   &+cavc(ixpp,iypp,izpp)*tzpp+cavc(ixpp,iypp,izm)*tzm
  d1(12)=cavc(ixpp,iym,iz)*tz+cavc(ixpp,iym,izp)*tzp &
   &+cavc(ixpp,iym,izpp)*tzpp+cavc(ixpp,iym,izm)*tzm

  d1(13)=cavc(ixm,iy,iz)*tz+cavc(ixm,iy,izp)*tzp+cavc(ixm,iy,izpp)*tzpp+ &
   &cavc(ixm,iy,izm)*tzm
  d1(14)=cavc(ixm,iyp,iz)*tz+cavc(ixm,iyp,izp)*tzp+cavc(ixm,iyp,izpp)*tzpp+ &
   &cavc(ixm,iyp,izm)*tzm
  d1(15)=cavc(ixm,iypp,iz)*tz+cavc(ixm,iypp,izp)*tzp &
   &+cavc(ixm,iypp,izpp)*tzpp+cavc(ixm,iypp,izm)*tzm
  d1(16)=cavc(ixm,iym,iz)*tz+cavc(ixm,iym,izp)*tzp+cavc(ixm,iym,izpp)*tzpp+ &
   &cavc(ixm,iym,izm)*tzm

! Buffer things for subsequent use.
  t1=d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym
  t2=d1(5)*ty+d1(6)*typ+d1(7)*typp+d1(8)*tym
  t3=d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym
  t4=d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym

! The function
  if(lw)then
   rpsi=t1*tx+t2*txp+t3*txpp+t4*txm
   if(.not.lgl)return
  endif ! lw

! Buffer things for subsequent use.
  dt1=d1(1)*dty+d1(2)*dtyp+d1(3)*dtypp+d1(4)*dtym
  dt2=d1(5)*dty+d1(6)*dtyp+d1(7)*dtypp+d1(8)*dtym
  dt3=d1(9)*dty+d1(10)*dtyp+d1(11)*dtypp+d1(12)*dtym
  dt4=d1(13)*dty+d1(14)*dtyp+d1(15)*dtypp+d1(16)*dtym
  d2t1=d1(1)*d2ty+d1(2)*d2typ+d1(3)*d2typp+d1(4)*d2tym
  d2t2=d1(5)*d2ty+d1(6)*d2typ+d1(7)*d2typp+d1(8)*d2tym
  d2t3=d1(9)*d2ty+d1(10)*d2typ+d1(11)*d2typp+d1(12)*d2tym
  d2t4=d1(13)*d2ty+d1(14)*d2typ+d1(15)*d2typp+d1(16)*d2tym

! The gradient, first term involving Theta'(x)
  grad(1)=t1*dtx+t2*dtxp+t3*dtxpp+t4*dtxm
! Second term of the gradient involving Theta'(y)
  grad(2)=dt1*tx+dt2*txp+dt3*txpp+dt4*txm

! The Laplacian: first term involving Theta''(x)
  templap(1)=t1*d2tx+t2*d2txp+t3*d2txpp+t4*d2txm
! Second term of the Laplacian involving Theta''(y)
  templap(2)=d2t1*tx+d2t2*txp+d2t3*txpp+d2t4*txm
! The Laplacian: term involving Theta'(x)Theta'(y)
  templap(4)=dt1*dtx+dt2*dtxp+dt3*dtxpp+dt4*dtxm

! More derivatives.
  d1(1)=cavc(ix,iy,iz)*d2tz+cavc(ix,iy,izp)*d2tzp+cavc(ix,iy,izpp)*d2tzpp+ &
   &cavc(ix,iy,izm)*d2tzm
  d1(2)=cavc(ix,iyp,iz)*d2tz+cavc(ix,iyp,izp)*d2tzp+cavc(ix,iyp,izpp)*d2tzpp+ &
   &cavc(ix,iyp,izm)*d2tzm
  d1(3)=cavc(ix,iypp,iz)*d2tz+cavc(ix,iypp,izp)*d2tzp &
   &+cavc(ix,iypp,izpp)*d2tzpp+cavc(ix,iypp,izm)*d2tzm
  d1(4)=cavc(ix,iym,iz)*d2tz+cavc(ix,iym,izp)*d2tzp+cavc(ix,iym,izpp)*d2tzpp+ &
   &cavc(ix,iym,izm)*d2tzm

  d1(5)=cavc(ixp,iy,iz)*d2tz+cavc(ixp,iy,izp)*d2tzp+cavc(ixp,iy,izpp)*d2tzpp+ &
   &cavc(ixp,iy,izm)*d2tzm
  d1(6)=cavc(ixp,iyp,iz)*d2tz+cavc(ixp,iyp,izp)*d2tzp &
   &+cavc(ixp,iyp,izpp)*d2tzpp+cavc(ixp,iyp,izm)*d2tzm
  d1(7)=cavc(ixp,iypp,iz)*d2tz+cavc(ixp,iypp,izp)*d2tzp &
   &+cavc(ixp,iypp,izpp)*d2tzpp+cavc(ixp,iypp,izm)*d2tzm
  d1(8)=cavc(ixp,iym,iz)*d2tz+cavc(ixp,iym,izp)*d2tzp &
   &+cavc(ixp,iym,izpp)*d2tzpp+cavc(ixp,iym,izm)*d2tzm

  d1(9)=cavc(ixpp,iy,iz)*d2tz+cavc(ixpp,iy,izp)*d2tzp &
   &+cavc(ixpp,iy,izpp)*d2tzpp+cavc(ixpp,iy,izm)*d2tzm
  d1(10)=cavc(ixpp,iyp,iz)*d2tz+cavc(ixpp,iyp,izp)*d2tzp+cavc(ixpp,iyp,izpp)* &
   &d2tzpp+cavc(ixpp,iyp,izm)*d2tzm
  d1(11)=cavc(ixpp,iypp,iz)*d2tz+cavc(ixpp,iypp,izp)*d2tzp &
   &+cavc(ixpp,iypp,izpp)*d2tzpp+cavc(ixpp,iypp,izm)*d2tzm
  d1(12)=cavc(ixpp,iym,iz)*d2tz+cavc(ixpp,iym,izp)*d2tzp+cavc(ixpp,iym,izpp)* &
   &d2tzpp+cavc(ixpp,iym,izm)*d2tzm

  d1(13)=cavc(ixm,iy,iz)*d2tz+cavc(ixm,iy,izp)*d2tzp &
   &+cavc(ixm,iy,izpp)*d2tzpp+cavc(ixm,iy,izm)*d2tzm
  d1(14)=cavc(ixm,iyp,iz)*d2tz+cavc(ixm,iyp,izp)*d2tzp &
   &+cavc(ixm,iyp,izpp)*d2tzpp+cavc(ixm,iyp,izm)*d2tzm
  d1(15)=cavc(ixm,iypp,iz)*d2tz+cavc(ixm,iypp,izp)*d2tzp+cavc(ixm,iypp,izpp)* &
   &d2tzpp+cavc(ixm,iypp,izm)*d2tzm
  d1(16)=cavc(ixm,iym,iz)*d2tz+cavc(ixm,iym,izp)*d2tzp &
   &+cavc(ixm,iym,izpp)*d2tzpp+cavc(ixm,iym,izm)*d2tzm

! Buffer things for subsequent use.
  t1=d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym
  t2=d1(5)*ty+d1(6)*typ+d1(7)*typp+d1(8)*tym
  t3=d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym
  t4=d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym

! Laplacian: term involving Theta''(z)
  templap(3)=t1*tx+t2*txp+t3*txpp+t4*txm

! And the third term of the gradient involving Theta'(z)
  d1(1)=cavc(ix,iy,iz)*dtz+cavc(ix,iy,izp)*dtzp+cavc(ix,iy,izpp)*dtzpp+ &
   &cavc(ix,iy,izm)*dtzm
  d1(2)=cavc(ix,iyp,iz)*dtz+cavc(ix,iyp,izp)*dtzp+cavc(ix,iyp,izpp)*dtzpp+ &
   &cavc(ix,iyp,izm)*dtzm
  d1(3)=cavc(ix,iypp,iz)*dtz+cavc(ix,iypp,izp)*dtzp+cavc(ix,iypp,izpp)*dtzpp+ &
   &cavc(ix,iypp,izm)*dtzm
  d1(4)=cavc(ix,iym,iz)*dtz+cavc(ix,iym,izp)*dtzp+cavc(ix,iym,izpp)*dtzpp+ &
   &cavc(ix,iym,izm)*dtzm

  d1(5)=cavc(ixp,iy,iz)*dtz+cavc(ixp,iy,izp)*dtzp+cavc(ixp,iy,izpp)*dtzpp+ &
   &cavc(ixp,iy,izm)*dtzm
  d1(6)=cavc(ixp,iyp,iz)*dtz+cavc(ixp,iyp,izp)*dtzp+cavc(ixp,iyp,izpp)*dtzpp+ &
   &cavc(ixp,iyp,izm)*dtzm
  d1(7)=cavc(ixp,iypp,iz)*dtz+cavc(ixp,iypp,izp)*dtzp &
   &+cavc(ixp,iypp,izpp)*dtzpp+cavc(ixp,iypp,izm)*dtzm
  d1(8)=cavc(ixp,iym,iz)*dtz+cavc(ixp,iym,izp)*dtzp+cavc(ixp,iym,izpp)*dtzpp+ &
   &cavc(ixp,iym,izm)*dtzm

  d1(9)=cavc(ixpp,iy,iz)*dtz+cavc(ixpp,iy,izp)*dtzp+cavc(ixpp,iy,izpp)*dtzpp+ &
   &cavc(ixpp,iy,izm)*dtzm
  d1(10)=cavc(ixpp,iyp,iz)*dtz+cavc(ixpp,iyp,izp)*dtzp &
   &+cavc(ixpp,iyp,izpp)*dtzpp+cavc(ixpp,iyp,izm)*dtzm
  d1(11)=cavc(ixpp,iypp,iz)*dtz+cavc(ixpp,iypp,izp)*dtzp &
   &+cavc(ixpp,iypp,izpp)*dtzpp+cavc(ixpp,iypp,izm)*dtzm
  d1(12)=cavc(ixpp,iym,iz)*dtz+cavc(ixpp,iym,izp)*dtzp &
   &+cavc(ixpp,iym,izpp)*dtzpp+cavc(ixpp,iym,izm)*dtzm

  d1(13)=cavc(ixm,iy,iz)*dtz+cavc(ixm,iy,izp)*dtzp+cavc(ixm,iy,izpp) &
   &*dtzpp+cavc(ixm,iy,izm)*dtzm
  d1(14)=cavc(ixm,iyp,iz)*dtz+cavc(ixm,iyp,izp)*dtzp+cavc(ixm,iyp,izpp)*dtzpp+&
   &cavc(ixm,iyp,izm)*dtzm
  d1(15)=cavc(ixm,iypp,iz)*dtz+cavc(ixm,iypp,izp)*dtzp &
   &+cavc(ixm,iypp,izpp)*dtzpp+cavc(ixm,iypp,izm)*dtzm
  d1(16)=cavc(ixm,iym,iz)*dtz+cavc(ixm,iym,izp)*dtzp &
   &+cavc(ixm,iym,izpp)*dtzpp+cavc(ixm,iym,izm)*dtzm

! Buffer things for subsequent use
  t1=d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym
  t2=d1(5)*ty+d1(6)*typ+d1(7)*typp+d1(8)*tym
  t3=d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym
  t4=d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym
  dt1=d1(1)*dty+d1(2)*dtyp+d1(3)*dtypp+d1(4)*dtym
  dt2=d1(5)*dty+d1(6)*dtyp+d1(7)*dtypp+d1(8)*dtym
  dt3=d1(9)*dty+d1(10)*dtyp+d1(11)*dtypp+d1(12)*dtym
  dt4=d1(13)*dty+d1(14)*dtyp+d1(15)*dtypp+d1(16)*dtym

! Gradient: term involving Theta'(z)
  grad(3)=t1*tx+t2*txp+t3*txpp+t4*txm

! The Laplacian: term involving Theta'(x)Theta'(z)
  templap(5)=t1*dtx+t2*dtxp+t3*dtxpp+t4*dtxm

! The Laplacian: term involving Theta'(y)Theta'(z)
  templap(6)=dt1*tx+dt2*txp+dt3*txpp+dt4*txm

! Transformation of gradient to the Cartesian grid
  grad(1:3)=matmul(bg,grad(1:3))

! The Laplacian: summing all contributions with appropriate transformation
  lap=templap(1)*lvp(1)+templap(2)*lvp(2)+templap(3)*lvp(3) &
   &+templap(4)*lvp(4)+templap(5)*lvp(6)+templap(6)*lvp(5)

 END SUBROUTINE blip3dk

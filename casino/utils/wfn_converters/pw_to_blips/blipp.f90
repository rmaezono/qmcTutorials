 SUBROUTINE blipp
!------------------------------------------------------------------!
! Read complex PW orbitals at gamma point only from pwfn.data and  !
! generate real blip orbitals.                                     !
!------------------------------------------------------------------!
  USE rng,ONLY : ranx,init_rng
  USE singleton,ONLY : fftn
  IMPLICIT NONE
  INTEGER i,j,k,ibnd,ig,lx,ly,lz,ik,spin,icount,ierr,ialloc,&
   &maxband,maxband_new
  REAL(dp) r(3),xb(5),xd(5),xbd(5),et,dot_prod,tempr, &
   &norm2_temp,overlap(5),av_overlap(5),av_overlap2(5)
  REAL(dp),POINTER :: ke_pworb(:,:,:),ke_bliporb(:,:,:)
  REAL(dp),POINTER :: ke_pworb_new(:,:,:),ke_bliporb_new(:,:,:)
  REAL(dp) :: d,d2,dg(3),blipval,bliplap,blipgrad(3)
  COMPLEX(dp) :: eigr
  COMPLEX(dp),ALLOCATABLE :: evc(:),fftdata(:)
  REAL(dp) norm2_real,norm2_imag,rdiff(3),r1(3),skin_thickness, &
   &cell_volume,norm2_pworb,rdiff_lv(3),ke_real,ke_imag, &
   &half_side_len_lv(3)
  INTEGER ig1,ncentres,max_no_loc_orbs,no_loc_orbs(2),icut,ix,iy,iz,nx,ny,nz, &
   &nxp,nyp,nzp,lx1,ly1,lz1,ibnd2
  INTEGER,ALLOCATABLE :: minus_g(:),norbs_on_centre(:,:)
  REAL(dp),ALLOCATABLE :: centrepos_lv(:,:,:),loc_radius(:), &
   &loc_centre(:,:),trunc_radius(:,:),norm2_frac(:,:)
  LOGICAL lreal,loc_orbs_present,binside
  LOGICAL,ALLOCATABLE :: use_this_g(:)
  CHARACTER(80) char80
  INTERFACE
   REAL(kind=kind(1.d0)) FUNCTION norm2_integrand_r(r)
    USE helpers,ONLY : dp
    REAL(dp),INTENT(in) :: r
   END FUNCTION norm2_integrand_r
   REAL(kind=kind(1.d0)) FUNCTION norm2_integrand_zlv(zlv)
    USE helpers,ONLY : dp
    REAL(dp),INTENT(in) :: zlv
   END FUNCTION norm2_integrand_zlv
   REAL(kind=kind(1.d0)) FUNCTION ke_integrand_r(r)
    USE helpers,ONLY : dp
    REAL(dp),INTENT(in) :: r
   END FUNCTION ke_integrand_r
   REAL(kind=kind(1.d0)) FUNCTION ke_integrand_zlv(zlv)
    USE helpers,ONLY : dp
    REAL(dp),INTENT(in) :: zlv
   END FUNCTION ke_integrand_zlv
  END INTERFACE

  if(Nkvec/=1)then
   call errstop('BLIPP','Bug - there should only be one k point: Gamma.')
  else
   write(6,*)'Generating real orbitals at the Gamma point.'
  endif ! Nkvec/=1
  write(6,*)

  cell_volume=abs(det_33(at))

  inquire(file='centres.dat',exist=loc_orbs_present)

  if(loc_orbs_present)then

! Read in data from the centres.dat file.
   open(io_centres,file='centres.dat',status='old',iostat=ierr)
   if(ierr/=0)call errstop('BLIPP','Cannot open centres.dat.  Stopping.')

   call skipio(io_centres,1)
   read(io_centres,*,err=30,end=40)ncentres
   write(6,*)'Number of centres: '//trim(i2s(ncentres))
   if(ncentres<0)call errstop('BLIPP','Need more centres!')

   if(ncentres>0)then
    allocate(loc_centre(3,ncentres),loc_radius(ncentres), &
     &norbs_on_centre(ncentres,2),stat=ialloc)
    if(ialloc/=0)call errstop('BLIPP','Allocation problem: localization &
     &centres etc.')
   endif ! ncentres

   call skipio(io_centres,3)
   read(io_centres,*,err=30,end=40)icut
   if(icut==2)then
    write(6,*)'A parallelepiped cutoff will be used for localized orbitals.'
   elseif(icut==1)then
    write(6,*)'A spherical cutoff will be used for localized orbitals.'
   else
    call errstop('BLIPP','The cutoff flag should be 1 (spherical cutoff) or 2 &
     &(parallelepiped cutoff).')
   endif ! icut

   call skipio(io_centres,1)

   do i=1,ncentres
    read(io_centres,*,err=30,end=40)loc_centre(1:3,i),loc_radius(i), &
     &norbs_on_centre(i,1:2)
    loc_centre(1:3,i)=loc_centre(1:3,i)+R_offset
    if(loc_radius(i)<=0.d0)call errstop('BLIPP','Cutoff radii in centres.dat &
     &must be positive.')
    if(any(norbs_on_centre(i,:)<0))call errstop('BLIPP','You seem to have a &
     &negative number of orbitals on a centre?')
   enddo ! i

! Total numbers of localized orbitals for up and down spins.
   if(ncentres>0)then
    no_loc_orbs(1)=sum(norbs_on_centre(:,1))
    no_loc_orbs(2)=sum(norbs_on_centre(:,2))
   else
    no_loc_orbs=0
   endif ! ncentres>0
   max_no_loc_orbs=maxval(no_loc_orbs)

! Read optional thickness of shell region; zero by default
   read(io_centres,*,iostat=ierr)
   if(ierr/=0)then
    skin_thickness=0.d0
   else
    read(io_centres,*,iostat=ierr)skin_thickness
    if(ierr/=0)skin_thickness=0.d0
   endif ! ierr/=0
   if(skin_thickness<0.d0)skin_thickness=0.d0
   write(6,*)'Minimum skin thickness        : ',skin_thickness
   write(6,*)

   close(io_centres)

! Check that cutoff radii are not too large
   if((maxval(loc_radius)+skin_thickness)*maxval(mag_rlv)>0.5d0)call &
    &errstop('BLIPP','Localization region does not fit into simulation cell.')

! Centre & radius of each localized orbital
   if(max_no_loc_orbs>0)then
    allocate(centrepos_lv(3,max_no_loc_orbs,no_spin_types), &
     &trunc_radius(max_no_loc_orbs,no_spin_types),stat=ialloc)
    if(ialloc/=0)call errstop('BLIPP','Allocation error: centre positions and &
     &trunc. radii.')
   endif ! max_no_loc_orbs>0
   do spin=1,no_spin_types
    ibnd2=0
    do i=1,ncentres
     do ibnd=1,norbs_on_centre(i,spin)
      ibnd2=ibnd2+1
      centrepos_lv(1:3,ibnd2,spin)=matmul(loc_centre(1:3,i),bg)
      trunc_radius(ibnd2,spin)=loc_radius(i)
     enddo ! band
    enddo ! i
   enddo ! spin
   if(ncentres>0)deallocate(loc_centre,loc_radius,norbs_on_centre)

   if(norm2_ke_calc)then
    allocate(norm2_frac(max_no_loc_orbs,no_spin_types),stat=ialloc)
    if(ialloc/=0)call errstop('BLIPP','Allocation error: norm2_frac.')
   endif ! norm2_ke_calc

  else

   ncentres=0
   no_loc_orbs=0
   max_no_loc_orbs=0

  endif ! loc_orbs_present

  allocate(fftdata(nr(1)*nr(2)*nr(3)),evc(NGvec), &
   &avc(0:nr(1)-1,0:nr(2)-1,0:nr(3)-1),minus_g(NGvec),use_this_g(NGvec),&
   &stat=ialloc)
  if(ialloc/=0)call errstop('BLIPP','Allocation problem: blip & PW &
   &coefficients, etc.')

! Set up to use only +G.
  if(any(g(1:3,1)/=0))call errstop('BLIPP','First G vector in pwfn.data has to &
   &be the zero vector.')
  use_this_g(:)=.true.
  minus_g(:)=0
  minus_g(1)=1
  do ig=2,NGvec-1
   if(use_this_g(ig))then
! Find -G
    do ig1=ig+1,NGvec
     if(all(g(1:3,ig)+g(1:3,ig1)==0))then
      use_this_g(ig1)=.false.
      minus_g(ig)=ig1
      minus_g(ig1)=ig
      exit
     endif ! G(ig)=-G(ig1)
    enddo ! ig1
   endif ! use_this _G
  enddo ! ig

!  write(6,*)'gvec pairing'
!  do ig=1,NGvec
!   write(6,*)ig,g(1:3,ig),minus_g(ig),g(1:3,minus_g(ig))
!  enddo

! There should just be the one k point.
  if(any(abs(kvec(1,:))>eps))call errstop('BLIPP','There should only be one k &
   &vector: Gamma.')

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
    maxband_new = maxval(nband(:ik,:))
    if(maxband_new > maxband)then
     allocate(ke_pworb_new(maxband_new,Nkvec,no_spin_types), &
      &ke_bliporb_new(maxband_new,Nkvec,no_spin_types),stat=ialloc)
     if(ialloc/=0)call errstop('BLIPP','Allocation error: KE arrays.')
     if(maxband>0)then
      ke_pworb_new(:maxband,:,:) = ke_pworb(:,:,:)
      ke_bliporb_new(:maxband,:,:) = ke_bliporb(:,:,:)
      deallocate(ke_pworb,ke_bliporb)
     endif
     maxband = maxband_new
     ke_pworb => ke_pworb_new
     ke_bliporb => ke_bliporb_new
    endif
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
     if(ibnd>no_loc_orbs(spin))then
      write(io_bwfn,'(1x,a,es24.17,a)')trim(i2s(ibnd))//" "//trim(i2s(spin)) &
       &//" ",et," F"
     else
      write(io_bwfn,'(1x,a,es24.17,a)')trim(i2s(ibnd))//" "//trim(i2s(spin)) &
       &//" ",et," T"
     endif ! ibnd>no_loc_orbs
     call skipio(io_pwfn,1)

! Read in the PW coefficients.
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

! Decide whether to take the real or the imaginary part.
     norm2_real=2.d0*dble(evc(1))**2
     norm2_imag=2.d0*aimag(evc(1))**2
     do ig=2,NGvec
      if(use_this_g(ig))then
       if(minus_g(ig)>0)then
        norm2_real=norm2_real+abs(evc(ig)+conjg(evc(minus_g(ig))))**2
        norm2_imag=norm2_imag+abs(evc(ig)-conjg(evc(minus_g(ig))))**2
       else
        norm2_real=norm2_real+abs(evc(ig))**2
        norm2_imag=norm2_imag+abs(evc(ig))**2
       endif ! minus_g(ig)>0
      endif ! use_this_g
     enddo ! ig
     lreal=(norm2_real>=norm2_imag)

     if(norm2_ke_calc)then
! Evaluate the KE and norm2 of the PW orbitals.
      ke_real=0.d0
      ke_imag=0.d0
      do ig=2,NGvec
       if(use_this_g(ig))then
        if(minus_g(ig)>0)then
         ke_real=ke_real+g2(ig)*abs(evc(ig)+conjg(evc(minus_g(ig))))**2
         ke_imag=ke_imag+g2(ig)*abs(evc(ig)-conjg(evc(minus_g(ig))))**2
        else
         ke_real=ke_real+g2(ig)*abs(evc(ig))**2
         ke_imag=ke_imag+g2(ig)*abs(evc(ig))**2
        endif ! minus_g(ig)>0
       endif ! use_this_g
      enddo ! ig
      if(lreal)then
       ke_pworb(ibnd,ik,spin)=0.5d0*ke_real/norm2_real
       norm2_pworb=0.5d0*norm2_real*cell_volume
      else
       ke_pworb(ibnd,ik,spin)=0.5d0*ke_imag/norm2_imag
       norm2_pworb=0.5d0*norm2_imag*cell_volume
      endif ! lreal
     endif ! norm2_ke_calc

     fftdata=czero
     if(lreal)then
 ! avc = \sum_g gamma(g) * Re[ evc(g) * exp( i.g.r ) ]
 ! If the wave function is real then evc(g) = conjg(evc(-g)) and the sum
 ! can be done only on +g because gamma(g) = gamma(-g).
 ! If the WF is not real but we want to take the real part then we take
 ! evc(ig) + conjg( evc(-g) ).
      fftdata(1)=cmplx(dble(evc(1))/3.375d0,0.d0,dp)
      do ig=2,NGvec
       if(use_this_g(ig))then
        icount=1+g(1,ig)+nr(1)*(indexfn(g(1,ig))+g(2,ig) &
         &+nr(2)*(indexfn(g(2,ig))+g(3,ig)+nr(3)*indexfn(g(3,ig))))
        if(minus_g(ig)>0)then
         fftdata(icount)=(evc(ig)+conjg(evc(minus_g(ig))))*gamma(ig)
        else
         fftdata(icount)=2.d0*evc(ig)*gamma(ig)
        endif ! minus_g(ig)>0
       endif ! use this G
      enddo ! ig
     else
 ! If the WF is not imaginary but we want to take the imaginary part then we
 ! take evc(ig)-conjg(evc(-g))
      fftdata(1)=cmplx(0.d0,aimag(evc(1))/3.375d0,dp)
      do ig=2,NGvec
       if(use_this_g(ig))then
        icount=1+g(1,ig)+nr(1)*(indexfn(g(1,ig))+g(2,ig) &
         &+nr(2)*(indexfn(g(2,ig))+g(3,ig)+nr(3)*indexfn(g(3,ig))))
        if(minus_g(ig)>0)then
         fftdata(icount)=(evc(ig)-conjg(evc(minus_g(ig))))*gamma(ig)
        else
         fftdata(icount)=2.d0*evc(ig)*gamma(ig)
        endif ! minus_g>0
       endif ! use this G
      enddo ! ig
     endif ! lreal

! Carry out FFT to real space
     call fftn(fftdata,nr,inv=.true.)
     fftdata=fftdata*sqrt(dble(nr(1)*nr(2)*nr(3)))

     icount=1
     if(lreal)then
      do lz=0,nr(3)-1
       do ly=0,nr(2)-1
        do lx=0,nr(1)-1
         avc(lx,ly,lz)=dble(fftdata(icount))
         icount=icount+1
        enddo ! lx
       enddo ! ly
      enddo ! lz
     else
      do lz=0,nr(3)-1
       do ly=0,nr(2)-1
        do lx=0,nr(1)-1
         avc(lx,ly,lz)=aimag(fftdata(icount))
         icount=icount+1
        enddo ! lx
       enddo ! ly
      enddo ! lz
     endif ! lreal

     if(ibnd>no_loc_orbs(spin))then

! Write out blip coefficients for extended orbital.
      write(io_bwfn,*)'Real blip coefficients for extended orbital'
      do lx=0,nr(1)-1
       do ly=0,nr(2)-1
        do lz=0,nr(3)-1
         write(io_bwfn,*)avc(lx,ly,lz)
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
       eps_simp=eps_norm2
       call adapt_simpson(norm2_integrand_zlv,&
        &rlv_min_pass(3),rlv_max_pass(3),norm2_temp)
       eps_simp=eps_ke
       call adapt_simpson(ke_integrand_zlv,&
        &rlv_min_pass(3),rlv_max_pass(3),ke_bliporb(ibnd,ik,spin))
       ke_bliporb(ibnd,ik,spin)=-0.5d0*ke_bliporb(ibnd,ik,spin)/norm2_temp
      endif ! norm2_ke_calc

     else

 ! Find grid point associated with position of centre.
      ix=floor(centrepos_lv(1,ibnd,spin)*rnr(1))
      iy=floor(centrepos_lv(2,ibnd,spin)*rnr(2))
      iz=floor(centrepos_lv(3,ibnd,spin)*rnr(3))

 ! Work out "start" and "stop" of localized grid, nx and nxp, etc.  Make sure
 ! that localized grid cannot exceed total grid.
      nxp=ceiling(mag_rlv(1)*(trunc_radius(ibnd,spin)+skin_thickness)*rnr(1))
      if(2*nxp+4>nr(1))nxp=(nr(1)-4)/2
      nyp=ceiling(mag_rlv(2)*(trunc_radius(ibnd,spin)+skin_thickness)*rnr(2))
      if(2*nyp+4>nr(2))nyp=(nr(2)-4)/2
      nzp=ceiling(mag_rlv(3)*(trunc_radius(ibnd,spin)+skin_thickness)*rnr(3))
      if(2*nzp+4>nr(3))nzp=(nr(3)-4)/2
      nx=ix-nxp-1
      ny=iy-nyp-1
      nz=iz-nzp-1
      nxp=ix+nxp+2
      nyp=iy+nyp+2
      nzp=iz+nzp+2

 ! Lengths of side of parallelepiped in terms of the lattice vectors.
      if(icut==2)half_side_len_lv(1:3)=mag_rlv(1:3)*trunc_radius(ibnd,spin)

 ! Calculate the fraction of norm^2 inside the sphere/cube.
 ! Note that we calculate the integral over the region inside the inner trunc.
 ! radius, i.e. we do not include the "skin" region.
      if(norm2_ke_calc)then
       centrepos_lv_pass=centrepos_lv(1:3,ibnd,spin)
       if(icut==1)then
        eps_simp=eps_norm2
        call adapt_simpson(norm2_integrand_r,0.d0,trunc_radius(ibnd,spin), &
         &norm2_temp)
        eps_simp=eps_ke
        call adapt_simpson(ke_integrand_r,0.d0,trunc_radius(ibnd,spin), &
         &ke_bliporb(ibnd,ik,spin))
       else
 ! Integrate over the coordinates in terms of the lattice vectors.
 ! Parallelepiped is a cuboid in this coordinate system.
 ! Must scale integrals by the cell volume afterwards.
        rlv_min_pass=centrepos_lv_pass-half_side_len_lv
        rlv_max_pass=centrepos_lv_pass+half_side_len_lv
        eps_simp=eps_norm2
        call adapt_simpson(norm2_integrand_zlv,rlv_min_pass(3),rlv_max_pass(3),&
         &norm2_temp)
        norm2_temp=norm2_temp*cell_volume
        eps_simp=eps_ke
        call adapt_simpson(ke_integrand_zlv,rlv_min_pass(3),rlv_max_pass(3), &
         &ke_bliporb(ibnd,ik,spin))
        ke_bliporb(ibnd,ik,spin)=ke_bliporb(ibnd,ik,spin)*cell_volume
       endif ! icut=1
       ke_bliporb(ibnd,ik,spin)=-0.5d0*ke_bliporb(ibnd,ik,spin)/norm2_temp
 ! Express norm as a fraction of the norm of the PW orbitals.
       norm2_frac(ibnd,spin)=norm2_temp/norm2_pworb
      endif ! norm2_ke_calc

      write(io_bwfn,*)'Centre (au), inner cutoff radius (au), &
       &outer cutoff radius (au) norm2, icut'
      if(norm2_ke_calc)then
       write(io_bwfn,'(3(1x,es24.17),2x,es24.17,1x,es24.17,2x,es24.17,2x,i1)') &
        &matmul(at,centrepos_lv(:,ibnd,spin)),trunc_radius(ibnd,spin), &
        &trunc_radius(ibnd,spin)+skin_thickness,norm2_frac(ibnd,spin),icut
      else
       write(io_bwfn,'(3(1x,es24.17),2x,es24.17,1x,es24.17,a)') &
        &matmul(at,centrepos_lv(:,ibnd,spin)),trunc_radius(ibnd,spin), &
        &trunc_radius(ibnd,spin)+skin_thickness," -1.d0 "//trim(i2s(icut))
      endif ! norm2_ke_calc
      write(io_bwfn,*)'Start of grid; grid dimensions'
      write(io_bwfn,*)trim(i2s(nx))//" "//trim(i2s(ny))//" "//trim(i2s(nz)) &
       &//" "//trim(i2s(nxp-nx+1))//" "//trim(i2s(nyp-ny+1))//" " &
       &//trim(i2s(nzp-nz+1))
      write(io_bwfn,*)'Real blip coefficients for localized orbital'
      do lx1=nx,nxp
       lx=modulo(lx1,nr(1))
       do ly1=ny,nyp
        ly=modulo(ly1,nr(2))
        do lz1=nz,nzp
         lz=modulo(lz1,nr(3))
         write(io_bwfn,*)avc(lx,ly,lz)
        enddo ! lz1
       enddo ! ly1
      enddo ! lx1

     endif ! Localized orbital or not.

! Carry out the overlap test described in the CASINO manual.
! Repeat the whole test n_overlap_tests times, to compute error bars.
     if(n_points_for_test>0)then
      call init_rng(12345678)

      av_overlap=0.d0 ; av_overlap2=0.d0
      do j=1,n_overlap_tests
       xb=0.d0 ; xd=0.d0 ; xbd=0.d0
       do i=1,n_points_for_test
        r(1)=ranx() ; r(2)=ranx() ; r(3)=ranx()
        if(ibnd<=no_loc_orbs(spin))then
         if(icut==1)then
 ! Spherical cutoff.
          r1=matmul(at,r-centrepos_lv(1:3,ibnd,spin))
          call min_image(r1,at,bg,rdiff,mag_rlv)
          binside=(dot_product(rdiff,rdiff)<=trunc_radius(ibnd,spin)**2)
         else
 ! Parallelepiped cutoff.
          rdiff_lv=modulo(r-centrepos_lv(1:3,ibnd,spin)+0.5d0,1.d0)-0.5d0
          binside=(all(abs(rdiff_lv)<=half_side_len_lv))
         endif ! icut
         if(binside)then
          call blip3d(blipval,bliplap,blipgrad,r,.true.,.true.)
         else
          blipval=0.d0 ; bliplap=0.d0 ; blipgrad=0.d0
         endif ! binside
        else
         call blip3d(blipval,bliplap,blipgrad,r,.true.,.true.)
        endif ! icount<no_loc_orbs
        d2=0.d0 ; dg=0.d0
        if(lreal)then
         d=dble(evc(1))
         do ig=2,NGvec
          if(use_this_g(ig))then
           dot_prod=twopi*(g(1,ig)*r(1)+g(2,ig)*r(2)+g(3,ig)*r(3))
           if(minus_g(ig)>0)then
            eigr=(evc(ig)+conjg(evc(minus_g(ig))))* &
             &cmplx(cos(dot_prod),sin(dot_prod),dp)
           else
            eigr=2.d0*evc(ig)*cmplx(cos(dot_prod),sin(dot_prod),dp)
           endif ! minus_G>0
           d=d+dble(eigr)
           d2=d2-dble(eigr)*g2(ig)
           dg=dg-(twopi*aimag(eigr))*g(1:3,ig)
          endif ! use_this_G
         enddo ! ig
        else
         d=aimag(evc(1))
         do ig=2,NGvec
          if(use_this_g(ig))then
           dot_prod=twopi*(g(1,ig)*r(1)+g(2,ig)*r(2)+g(3,ig)*r(3))
           if(minus_g(ig)>0)then
            eigr=(evc(ig)-conjg(evc(minus_g(ig)))) &
             &*cmplx(cos(dot_prod),sin(dot_prod),dp)
           else
            eigr=2.d0*evc(ig)*cmplx(cos(dot_prod),sin(dot_prod),dp)
           endif ! minus_g>0
           d=d+aimag(eigr)
           d2=d2-aimag(eigr)*g2(ig)
           dg=dg+(twopi*dble(eigr))*g(1:3,ig)
          endif ! use_this_G
         enddo ! ig
        endif ! lreal
        dg=matmul(bg,dg)
! Overlap of orbital values.
        xb(1)=xb(1)+blipval**2
        xd(1)=xd(1)+d**2
        xbd(1)=xbd(1)+blipval*d
! Overlap of orbital Laplacians.
        xb(2)=xb(2)+bliplap**2
        xd(2)=xd(2)+d2**2
        xbd(2)=xbd(2)+bliplap*d2
! Overlap of orbital gradients.
        xb(3:5)=xb(3:5)+blipgrad(1:3)**2
        xd(3:5)=xd(3:5)+dg(1:3)**2
        xbd(3:5)=xbd(3:5)+blipgrad(1:3)*dg(1:3)
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

! If localized orbitals are present, write out the fraction of the norm^2
! contained in the localized orbital
  if(max_no_loc_orbs>0.and.norm2_ke_calc)then
   write(6,*)'Fraction of norm^2 in localized orbitals:'
   if(spin_polarized)then
    write(6,*)'Spin ;  Band ;   Frac. norm^2'
    do spin=1,no_spin_types
     do ibnd=1,no_loc_orbs(spin)
      write(6,'(3x,i1,2x,i6,4x,f12.8)')spin,ibnd,norm2_frac(ibnd,spin)
     enddo ! ibnd
    enddo ! spin
   else
    write(6,*)'    Band ;   Frac. norm^2'
    do ibnd=1,no_loc_orbs(1)
     write(6,'(3x,i6,3x,f12.8)')ibnd,norm2_frac(ibnd,1)
    enddo ! ibnd
   endif ! spin_polarized
   write(6,*)
  endif ! Localized orbitals present

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
  deallocate(avc,use_this_g,minus_g)
  if(allocated(centrepos_lv))deallocate(centrepos_lv)
  if(allocated(trunc_radius))deallocate(trunc_radius)
  if(associated(ke_pworb))deallocate(ke_pworb)
  if(associated(ke_bliporb))deallocate(ke_bliporb)
  if(allocated(norm2_frac))deallocate(norm2_frac)

  return

10 call errstop('BLIPP','Error reading pwfn.data.')
20 call errstop('BLIPP','Unexpectedly reached end of pwfn.data.')
30 call errstop('BLIPP','Error reading centres.dat.')
40 call errstop('BLIPP','Unexpectedly reached end of centres.dat.')

 END SUBROUTINE blipp


 SUBROUTINE blip3d(rpsi,lap,grad,r,lw,lgl)
!----------------------------------------------------------------------------!
! This subroutine evaluates the value of a function, its gradient and its    !
! Laplacian at a vector point r, using the overlapping of blip functions.    !
! The blip grid is defined on a cubic cell, so r should always be given in   !
! units of the crystal lattice vectors.                                      !
!----------------------------------------------------------------------------!
  IMPLICIT NONE
  REAL(dp),INTENT(in) :: r(3)
  LOGICAL,INTENT(in) :: lw,lgl
  REAL(dp),INTENT(out) :: rpsi,lap,grad(3)
  INTEGER ix,iy,iz,ixp,ixpp,ixm,iyp,iypp,iym,izp,izpp,izm
  REAL(dp) t(3),txm,tx,txp,txpp,tym,ty, &
   &typ,typp,tzm,tz,tzp,tzpp,dtxm,dtx,dtxp,dtxpp,dtym,dty,dtyp,dtypp,dtzm,&
   &dtz,dtzp,dtzpp,d2txm,d2tx,d2txp,d2txpp,d2tym,d2ty,d2typ,d2typp,d2tzm,&
   &d2tz,d2tzp,d2tzpp,x,y,z
  REAL(dp) d1(16),t1,t2,t3,t4,dt1,dt2,dt3,dt4,d2t1,d2t2,d2t3,d2t4,templap(6)

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

  d1(1)=avc(ix,iy,iz)*tz+avc(ix,iy,izp)*tzp+avc(ix,iy,izpp)*tzpp &
   &+avc(ix,iy,izm)*tzm
  d1(2)=avc(ix,iyp,iz)*tz+avc(ix,iyp,izp)*tzp+avc(ix,iyp,izpp)*tzpp &
   &+avc(ix,iyp,izm)*tzm
  d1(3)=avc(ix,iypp,iz)*tz+avc(ix,iypp,izp)*tzp+avc(ix,iypp,izpp)*tzpp+ &
   &avc(ix,iypp, izm)*tzm
  d1(4)=avc(ix,iym,iz)*tz+avc(ix,iym,izp)*tzp+avc(ix,iym,izpp)*tzpp+ &
   &avc(ix,iym,izm)*tzm

  d1(5)=avc(ixp,iy,iz)*tz+avc(ixp,iy,izp)*tzp+avc(ixp,iy,izpp)*tzpp+ &
   &avc(ixp,iy,izm)*tzm
  d1(6)=avc(ixp,iyp,iz)*tz+avc(ixp,iyp,izp)*tzp+avc(ixp,iyp,izpp)*tzpp+ &
   &avc(ixp, iyp, izm)*tzm
  d1(7)=avc(ixp,iypp,iz)*tz+avc(ixp,iypp,izp)*tzp+avc(ixp,iypp,izpp)*tzpp+ &
   &avc(ixp,iypp,izm)*tzm
  d1(8)=avc(ixp,iym,iz)*tz+avc(ixp,iym,izp)*tzp+avc(ixp,iym,izpp)*tzpp+ &
   &avc(ixp,iym,izm)*tzm

  d1(9)=avc(ixpp,iy,iz)*tz+avc(ixpp,iy,izp)*tzp+avc(ixpp,iy,izpp)*tzpp+ &
   &avc(ixpp,iy,izm)*tzm
  d1(10)=avc(ixpp,iyp,iz)*tz+avc(ixpp,iyp,izp)*tzp+avc(ixpp,iyp,izpp)*tzpp+ &
   &avc(ixpp,iyp,izm)*tzm
  d1(11)=avc(ixpp,iypp,iz)*tz+avc(ixpp,iypp,izp)*tzp &
   &+avc(ixpp,iypp,izpp)*tzpp+avc(ixpp,iypp,izm)*tzm
  d1(12)=avc(ixpp,iym,iz)*tz+avc(ixpp,iym,izp)*tzp+avc(ixpp,iym,izpp)*tzpp+ &
   &avc(ixpp,iym,izm)*tzm

  d1(13)=avc(ixm,iy,iz)*tz+avc(ixm,iy,izp)*tzp+avc(ixm,iy,izpp)*tzpp+ &
   &avc(ixm,iy,izm)*tzm
  d1(14)=avc(ixm,iyp,iz)*tz+avc(ixm,iyp,izp)*tzp+avc(ixm,iyp,izpp)*tzpp+ &
   &avc(ixm,iyp,izm)*tzm
  d1(15)=avc(ixm,iypp,iz)*tz+avc(ixm,iypp,izp)*tzp+avc(ixm,iypp,izpp)*tzpp+ &
   &avc(ixm,iypp,izm)*tzm
  d1(16)=avc(ixm,iym,iz)*tz+avc(ixm,iym,izp)*tzp+avc(ixm,iym,izpp)*tzpp+ &
   &avc(ixm,iym,izm)*tzm

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

! And now the third term of the Laplacian involving Theta''(z)
  d1(1)=avc(ix,iy,iz)*d2tz+avc(ix,iy,izp)*d2tzp+avc(ix,iy,izpp)*d2tzpp+ &
   &avc(ix,iy,izm)*d2tzm
  d1(2)=avc(ix,iyp,iz)*d2tz+avc(ix,iyp,izp)*d2tzp+avc(ix,iyp,izpp)*d2tzpp+ &
   &avc(ix,iyp,izm)*d2tzm
  d1(3)=avc(ix,iypp,iz)*d2tz+avc(ix,iypp,izp)*d2tzp+avc(ix,iypp,izpp)*d2tzpp+ &
   &avc(ix,iypp,izm)*d2tzm
  d1(4)=avc(ix,iym,iz)*d2tz+avc(ix,iym,izp)*d2tzp+avc(ix,iym,izpp)*d2tzpp+ &
   &avc(ix,iym,izm)*d2tzm

  d1(5)=avc(ixp,iy,iz)*d2tz+avc(ixp,iy,izp)*d2tzp+avc(ixp,iy,izpp)*d2tzpp+ &
   &avc(ixp,iy,izm)*d2tzm
  d1(6)=avc(ixp,iyp,iz)*d2tz+avc(ixp,iyp,izp)*d2tzp+avc(ixp,iyp,izpp)*d2tzpp+ &
   &avc(ixp,iyp,izm)*d2tzm
  d1(7)=avc(ixp,iypp,iz)*d2tz+avc(ixp,iypp,izp)*d2tzp+avc(ixp,iypp,izpp)* &
   &d2tzpp+avc(ixp,iypp,izm)*d2tzm
  d1(8)=avc(ixp,iym,iz)*d2tz+avc(ixp,iym,izp)*d2tzp+avc(ixp,iym,izpp)*d2tzpp+ &
   &avc(ixp,iym,izm)*d2tzm

  d1(9)=avc(ixpp,iy,iz)*d2tz+avc(ixpp,iy,izp)*d2tzp+avc(ixpp,iy,izpp)*d2tzpp+ &
   &avc(ixpp,iy,izm)*d2tzm
  d1(10)=avc(ixpp,iyp,iz)*d2tz+avc(ixpp,iyp,izp)*d2tzp+avc(ixpp,iyp,izpp)* &
   &d2tzpp+avc(ixpp,iyp,izm)*d2tzm
  d1(11)=avc(ixpp,iypp,iz)*d2tz+avc(ixpp,iypp,izp)*d2tzp+avc(ixpp,iypp,izpp)* &
   &d2tzpp+avc(ixpp,iypp,izm)*d2tzm
  d1(12)=avc(ixpp,iym,iz)*d2tz+avc(ixpp,iym,izp)*d2tzp+avc(ixpp,iym,izpp)* &
   &d2tzpp+avc(ixpp,iym,izm)*d2tzm

  d1(13)=avc(ixm,iy,iz)*d2tz+avc(ixm,iy,izp)*d2tzp+avc(ixm,iy,izpp)*d2tzpp+ &
   &avc(ixm,iy,izm)*d2tzm
  d1(14)=avc(ixm,iyp,iz)*d2tz+avc(ixm,iyp,izp)*d2tzp &
   &+avc(ixm,iyp,izpp)*d2tzpp+avc(ixm,iyp,izm)*d2tzm
  d1(15)=avc(ixm,iypp,iz)*d2tz+avc(ixm,iypp,izp)*d2tzp+avc(ixm,iypp,izpp)* &
   &d2tzpp+avc(ixm,iypp,izm)*d2tzm
  d1(16)=avc(ixm,iym,iz)*d2tz+avc(ixm,iym,izp)*d2tzp &
   &+avc(ixm,iym,izpp)*d2tzpp+avc(ixm,iym,izm)*d2tzm

! Buffer things for subsequent use.
  t1=d1(1)*ty+d1(2)*typ+d1(3)*typp+d1(4)*tym
  t2=d1(5)*ty+d1(6)*typ+d1(7)*typp+d1(8)*tym
  t3=d1(9)*ty+d1(10)*typ+d1(11)*typp+d1(12)*tym
  t4=d1(13)*ty+d1(14)*typ+d1(15)*typp+d1(16)*tym

! Laplacian: term involving Theta''(z)
  templap(3)=t1*tx+t2*txp+t3*txpp+t4*txm

! And the third term of the gradient involving Theta'(z)
  d1(1)=avc(ix,iy,iz)*dtz+avc(ix,iy,izp)*dtzp+avc(ix,iy,izpp)*dtzpp+ &
   &avc(ix,iy,izm)*dtzm
  d1(2)=avc(ix,iyp,iz)*dtz+avc(ix,iyp,izp)*dtzp+avc(ix,iyp,izpp)*dtzpp+ &
   &avc(ix,iyp,izm)*dtzm
  d1(3)=avc(ix,iypp,iz)*dtz+avc(ix,iypp,izp)*dtzp+avc(ix,iypp,izpp)*dtzpp+ &
   &avc(ix,iypp,izm)*dtzm
  d1(4)=avc(ix,iym,iz)*dtz+avc(ix,iym,izp)*dtzp+avc(ix,iym,izpp)*dtzpp+ &
   &avc(ix, iym, izm)*dtzm

  d1(5)=avc(ixp,iy,iz)*dtz+avc(ixp,iy,izp)*dtzp+avc(ixp,iy,izpp)*dtzpp+ &
   &avc(ixp,iy,izm)*dtzm
  d1(6)=avc(ixp,iyp,iz)*dtz+avc(ixp,iyp,izp)*dtzp+avc(ixp,iyp,izpp)*dtzpp+ &
   &avc(ixp,iyp,izm)*dtzm
  d1(7)=avc(ixp,iypp,iz)*dtz+avc(ixp,iypp,izp)*dtzp+avc(ixp,iypp,izpp)*dtzpp+ &
   &avc(ixp,iypp,izm)*dtzm
  d1(8)=avc(ixp,iym,iz)*dtz+avc(ixp,iym,izp)*dtzp+avc(ixp,iym,izpp)*dtzpp+ &
   &avc(ixp, iym, izm)*dtzm

  d1(9)=avc(ixpp,iy,iz)*dtz+avc(ixpp,iy,izp)*dtzp+avc(ixpp,iy,izpp)*dtzpp+ &
   &avc(ixpp,iy,izm)*dtzm
  d1(10)=avc(ixpp,iyp,iz)*dtz+avc(ixpp,iyp,izp)*dtzp &
   &+avc(ixpp,iyp,izpp)*dtzpp+avc(ixpp,iyp,izm)*dtzm
  d1(11)=avc(ixpp,iypp,iz)*dtz+avc(ixpp,iypp,izp)*dtzp+avc(ixpp,iypp,izpp)* &
   &dtzpp+avc(ixpp,iypp,izm)*dtzm
  d1(12)=avc(ixpp,iym,iz)*dtz+avc(ixpp,iym,izp)*dtzp &
   &+avc(ixpp,iym,izpp)*dtzpp+avc(ixpp,iym,izm)*dtzm

  d1(13)=avc(ixm,iy,iz)*dtz+avc(ixm,iy,izp)*dtzp+avc(ixm,iy,izpp)*dtzpp+ &
   &avc(ixm,iy,izm)*dtzm
  d1(14)=avc(ixm,iyp,iz)*dtz+avc(ixm,iyp,izp)*dtzp+avc(ixm,iyp,izpp)*dtzpp+ &
   &avc(ixm,iyp,izm)*dtzm
  d1(15)=avc(ixm,iypp,iz)*dtz+avc(ixm,iypp,izp)*dtzp &
   &+avc(ixm,iypp,izpp)*dtzpp+avc(ixm,iypp,izm)*dtzm
  d1(16)=avc(ixm,iym,iz)*dtz+avc(ixm,iym,izp)*dtzp+avc(ixm,iym,izpp)*dtzpp+ &
   &avc(ixm,iym,izm)*dtzm

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

 END SUBROUTINE blip3d

PROGRAM supercell
!------------------------------------------------------------------------------!
! SUPERCELL                                                                    !
! =========                                                                    !
! Supercell generator.  Given the primitive lattice vectors of a structure and !
! a target number of primitive cells, this program computes the supercell      !
! lattice vectors that maximize the distance between a point and its closest   !
! periodic image.                                                              !
!                                                                              !
! Note about MATMUL and indexing: we use Proper Fortran Indexing, while        !
! the MATMUL intrinsic function behaves as if indices were not at all          !
! entirely the wrong way around in Fortran or anything.  So, e.g.,             !
! Asim = Smat x Aprim translates into Asim = matmul(Aprim, Smat).              !
!------------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,PARAMETER :: dp=kind(1.d0)
 INTEGER dimensionality,npcells,Smat(3,3),Smat_diag(3,3),i,ierr
 REAL(dp) Aprim(3,3),Aprim_inv(3,3),volume_prim,k_offset(3)
 LOGICAL have_offset,do_diag,do_non_diag,use_inversion
 CHARACTER(80) char80
! Constants.
 REAL(dp),PARAMETER :: two_pi=2.d0*3.14159265358979324d0,inv_twopi=1.d0/two_pi

! Print header.
 write(6,*)
 write(6,*)'SUPERCELL'
 write(6,*)'========='
 write(6,*)'Constructs simulation supercells that maximize the radius of &
  &the sphere'
 write(6,*)'that can be inscribed in the Wigner-Seitz cell of the supercell.'
 write(6,*)'Produces:'
 write(6,*)'* A "supercell matrix" to be given in the "scell_matrix" input &
  &block in CASINO'
 write(6,*)'* k-point grids to be passed to the DFT code for generating the &
  &orbitals.'
 write(6,*)

! Read Aprim.
 Aprim=0.d0
 write(6,*)'Enter the primitive cell vectors (d numbers per line, d lines) &
  &[a.u.]:'
 write(6,*)'A1:'
 read(5,'(a)',iostat=ierr)char80
 if(ierr/=0)then
  write(6,*)'Could not read line from stdin.'
  stop
 endif
 do dimensionality=3,1,-1
  read(char80,*,iostat=ierr)Aprim(1:dimensionality,1)
  if(ierr==0)exit
 enddo
 if(dimensionality<1)then
  write(6,*)'Could not parse A1.'
  stop
 endif
 do i=2,dimensionality
  write(char80,*)i
  write(6,*)'A'//trim(adjustl(char80))//':'
  read(5,'(a)',iostat=ierr)char80
  if(ierr/=0)then
   write(6,*)'Could not read line from stdin.'
   stop
  endif
  read(char80,*,iostat=ierr)Aprim(1:dimensionality,i)
  if(ierr/=0)then
   write(char80,*)dimensionality
   write(6,*)'Could not parse '//trim(adjustl(char80))//'-dimensional vector &
    &from line.'
   stop
  endif
 enddo ! i
 write(6,*)

! Stop if dimensionality is 1.
 if(dimensionality==1)then
  write(6,*)'One-dimensional systems do not need the treatment offered by &
   &this utility.'
  stop
 endif

  ! Invert primitive cell matrix and check it's not singular.
 Aprim_inv=inverse(dimensionality,Aprim,volume_prim)
 if(are_equal(volume_prim,0.d0))then
  write(6,*)'The lattice vectors provided are colinear.'
  stop
 endif

! Read number of primitive cells in supercell.
 write(6,*)'Enter number of primitive cells in the supercell:'
 read(5,*,iostat=ierr)npcells
 if(ierr/=0)then
  write(6,*)'Could not read line from stdin.'
  stop
 endif
 if(npcells<1)then
  write(6,*)'Supercell must contain a positive number of primitive cells.'
  stop
 endif
 write(6,*)

! Read operation mode.
 write(6,*)'Enter one of the following letters (empty for default, "G"):'
 write(6,*)'  G : construct all general non-diagonal supercell matrices'
 write(6,*)'  D : restrict search to diagonal supercell matrices (quick)'
 write(6,*)'  B : do both (useful for diagonal/non-diagonal comparisons)'
 read(5,'(a)',iostat=ierr)char80
 if(ierr/=0)then
  write(6,*)'Could not read line from stdin.'
  stop
 endif
 char80=adjustl(char80)
 select case(char80(1:1))
 case('g','G',' ')
  do_diag=.false.
  do_non_diag=.true.
 case('d','D')
  do_diag=.true.
  do_non_diag=.false.
 case('b','B')
  do_diag=.true.
  do_non_diag=.true.
 case default
  write(6,*)'Invalid option.'
  stop
 end select
 write(6,*)

! Read offset.
 have_offset=.true.
 k_offset=0.d0
 write(6,*)'Enter k-point offset as fractions of the reciprocal lattice &
  &vectors of the'
 write(6,*)'primitive lattice (leave blank to generate all possible k-point &
  &grids'
 write(6,*)'corresponding to real wave functions):'
 read(5,'(a)',iostat=ierr)char80
 if(ierr/=0)then
  write(6,*)'Could not read line from stdin.'
  stop
 endif
 read(char80,*,iostat=ierr)k_offset(1:dimensionality)
 have_offset=ierr==0
 write(6,*)

! Use inversion symmetry or not?
 write(6,*)'Use inversion symmetry to reduce the number of k vectors? (Y/N)'
 read(5,'(a)',iostat=ierr)char80
 if(ierr/=0)then
  write(6,*)'Could not read line from stdin.'
  stop
 endif
 char80=adjustl(char80)
 select case(char80(1:1))
 case('y','Y')
  use_inversion=.true.
 case('n','N')
  use_inversion=.false.
 case default
  write(6,*)'Invalid option.'
  stop
 end select
 write(6,*)

! Report input.
 write(6,*)'Supercell parameters'
 write(6,*)'===================='
 if(volume_prim>0.d0)then
  write(6,*)'Primitive cell vectors (right-handed):'
 else
  write(6,*)'Primitive cell vectors (left-handed):'
 endif
 do i=1,dimensionality
  write(6,*)Aprim(1:dimensionality,i)
 enddo ! i
 write(6,*)
 write(6,'(1x,a,i4)')'Number of primitive cells in supercell: ',npcells
 write(6,*)

! Process diagonal version.
 if(do_diag)then
  write(6,*)'Best geometry with DIAGONAL supercell matrix'
  write(6,*)'============================================'
  call find_Smat(dimensionality,Aprim,npcells,Smat_diag,diagonal=.true.)
  call process_Smat(dimensionality,Smat_diag,Aprim,npcells,have_offset,k_offset)
 endif

! Process non-diagonal version.
 if(do_non_diag)then
  write(6,*)'Best geometry with GENERAL supercell matrix'
  write(6,*)'==========================================='
  call find_Smat(dimensionality,Aprim,npcells,Smat)
  call process_Smat(dimensionality,Smat,Aprim,npcells,have_offset,k_offset)
 endif


CONTAINS


 SUBROUTINE find_Smat(dimensionality,Aprim,npcells,Smat,diagonal)
!-----------------------------------------------------------------------------!
! Find the supercell matrix that produces the most spherical simulation cell. !
!-----------------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: dimensionality,npcells
  REAL(dp),INTENT(in) :: Aprim(3,3)
  INTEGER,INTENT(inout) :: Smat(3,3)
  LOGICAL,INTENT(in),OPTIONAL :: diagonal
  INTEGER sa,sc,sf,sb,sd,se
  REAL(dp) wsr,largest_wsr,Asim(3,3),Asim_inv(3,3),volume_sim,Bsim(3,3)
  LOGICAL diagonal_only

! Initialize
  largest_wsr=0.d0
  Smat=0
  diagonal_only=.false.
  if(present(diagonal))diagonal_only=diagonal

  select case(dimensionality)
  case(3)
! Find factorizations of npcells into three integers.
   do sa=1,npcells
    do sc=1,npcells/sa
     sf=npcells/(sa*sc)
     if(sa*sc*sf/=npcells)cycle
! Loop over valid values of the non-diagonal elements.
     do sb=0,sc-1
      do sd=0,sf-1
       do se=0,sf-1
        if(diagonal_only.and.any((/sb,sd,se/)/=0))cycle
! Get supercell lattice vectors.
        Asim(1:3,1)=sa*Aprim(1:3,1)+sb*Aprim(1:3,2)+sd*Aprim(1:3,3)
        Asim(1:3,2)=                sc*Aprim(1:3,2)+se*Aprim(1:3,3)
        Asim(1:3,3)=                                sf*Aprim(1:3,3)
! Get WSR and store Smat if largest so far.
        Asim_inv=inverse(dimensionality,Asim,volume_sim)
        Bsim=transpose(Asim_inv)*two_pi
        wsr=get_wsr(dimensionality,Asim,Bsim)
        if(wsr>largest_wsr)then
         Smat(1:3,1)=(/sa,sb,sd/)
         Smat(1:3,2)=(/ 0,sc,se/)
         Smat(1:3,3)=(/ 0, 0,sf/)
         largest_wsr=wsr
        endif ! wsr>largest_wsr
       enddo ! se
      enddo ! sd
     enddo ! sb
    enddo ! sc (determines sf)
   enddo ! sa

  case(2)
! Find factorizations of npcells into two integers.
   do sa=1,npcells
    sc=npcells/sa
    if(sa*sc/=npcells)cycle
! Loop over valid values of the non-diagonal elements.
    do sb=0,sc-1
     if(diagonal_only.and.sb/=0)cycle
! Get supercell lattice vectors.
     Asim(1:2,1)=sa*Aprim(1:2,1)+sb*Aprim(1:2,2)
     Asim(1:2,2)=                sc*Aprim(1:2,2)
! Get WSR and store Smat if largest so far.
     Asim_inv=inverse(dimensionality,Asim,volume_sim)
     Bsim=transpose(Asim_inv)*two_pi
     wsr=get_wsr(dimensionality,Asim,Bsim)
     if(wsr>largest_wsr)then
      Smat(1:2,1)=(/sa,sb/)
      Smat(1:2,2)=(/ 0,sc/)
      largest_wsr=wsr
     endif ! wsr>largest_wsr
    enddo ! sb
   enddo ! sa

  end select ! dimensionality

 END SUBROUTINE find_Smat


 SUBROUTINE process_Smat(dimensionality,Smat,Aprim,npcells,have_offset,k_offset)
!-----------------------------------------------------------------!
! Given a supercell matrix, put it in Minkowski-reduced form, and !
! evaluate and print the k points with offset K_SHIFT.            !
!-----------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: dimensionality,npcells
  REAL(dp),INTENT(in) :: Aprim(3,3),k_offset(3)
  LOGICAL,INTENT(in) :: have_offset
  INTEGER,INTENT(inout) :: Smat(3,3)
  INTEGER ialloc,ik,iks,nks,i
  INTEGER,ALLOCATABLE :: k_weight(:)
  REAL(dp) Asim(3,3),Bsim(3,3),Asim_inv(3,3),wsr,volume_sim,Aref(3,3),&
   &Bref(3,3),Aref_inv(3,3),wsr_ref,volume_ref,ks(3),ks_real(3,8),Sinv(3,3),dum
  REAL(dp),ALLOCATABLE :: k_gamma(:,:),k_point(:,:)

! Minkowski-reduce S.
  call minkowski_reduce_prod(dimensionality,Smat,Aprim)

! Get supercell vectors, reciprocals and WSRs.
  Asim=matmul(Aprim,Smat)
  Asim_inv=inverse(dimensionality,Asim,volume_sim)
  Bsim=transpose(Asim_inv)*two_pi
  wsr=get_wsr(dimensionality,Asim,Bsim)

! Construct reference lattice with same volume.
  Aref=0.d0
  select case(dimensionality)
  case(3) ! FCC lattice maximizes WSR
   Aref(1:3,1)=(/0.d0,1.d0,1.d0/)
   Aref(1:3,2)=(/1.d0,0.d0,1.d0/)
   Aref(1:3,3)=(/1.d0,1.d0,0.d0/)
  case(2) ! hexagonal lattice maximizes WSR
   Aref(1:2,1)=(/1.d0,0.d0/)
   Aref(1:2,2)=(/-0.5d0,0.5d0*sqrt(3.d0)/)
  end select
! Normalize to correct volume.
  Aref_inv=inverse(dimensionality,Aref,volume_ref)
  Aref=Aref*abs(volume_sim/volume_ref)**(1.d0/dble(dimensionality))
! Obtain WSR for comparison.
  Aref_inv=inverse(dimensionality,Aref,volume_ref)
  Bref=transpose(Aref_inv)*two_pi
  wsr_ref=get_wsr(dimensionality,Aref,Bref)

! Report results.
  write(6,*)'Supercell matrix:'
  write(char80,*)dimensionality
  do i=1,dimensionality
   write(6,'(1x,'//trim(adjustl(char80))//'(2x,i4))')Smat(1:dimensionality,i)
  enddo ! i
  write(6,*)
  write(6,*)'Wigner-Seitz radius:'
  write(char80,'(f5.1)')100.d0*wsr/wsr_ref
  select case(dimensionality)
  case(3)
   write(6,'(1x,f16.8,a)')wsr,' ('//trim(adjustl(char80))//'% of FCC)'
  case(2)
   write(6,'(1x,f16.8,a)')wsr,' ('//trim(adjustl(char80))//'% of hexagonal)'
  end select
  write(6,*)
  if(volume_sim>0.d0)then
   write(6,*)'Simulation cell vectors (right-handed):'
  else
   write(6,*)'Simulation cell vectors (left-handed):'
  endif
  do i=1,dimensionality
   write(6,*)Asim(1:dimensionality,i)
  enddo ! i
  write(6,*)

! Allocate k point arrays.
  allocate(k_gamma(3,npcells),k_point(3,npcells),k_weight(npcells),stat=ialloc)
  if(ialloc/=0)then
   write(6,*)'Allocation problem (k_point).'
   stop
  endif

! Generate k points at gamma and shift to desired offset.
  call get_k_points(dimensionality,Smat,npcells,k_gamma)

! Construct list of valid offsets.
  if(have_offset)then
   nks=1
   ks_real(1:dimensionality,1)=k_offset(1:dimensionality)
  else
   Sinv=inverse(dimensionality,dble(Smat),dum)
   Sinv=transpose(Sinv)
   select case(dimensionality)
   case(3)
    nks=8
    ks_real(1:3,1)=0.d0
    ks_real(1:3,2)=modulo(0.5d0*Sinv(1:3,1)+0.5d0,1.d0)-0.5d0
    ks_real(1:3,3)=modulo(0.5d0*Sinv(1:3,2)+0.5d0,1.d0)-0.5d0
    ks_real(1:3,4)=modulo(0.5d0*Sinv(1:3,3)+0.5d0,1.d0)-0.5d0
    ks_real(1:3,5)=modulo(0.5d0*(Sinv(1:3,1)+Sinv(1:3,2))+0.5d0,1.d0)-0.5d0
    ks_real(1:3,6)=modulo(0.5d0*(Sinv(1:3,1)+Sinv(1:3,3))+0.5d0,1.d0)-0.5d0
    ks_real(1:3,7)=modulo(0.5d0*(Sinv(1:3,2)+Sinv(1:3,3))+0.5d0,1.d0)-0.5d0
    ks_real(1:3,8)=modulo(0.5d0*(Sinv(1:3,1)+Sinv(1:3,2)+Sinv(1:3,3))+0.5d0,&
     &1.d0)-0.5d0
   case(2)
    nks=4
    ks_real(1:2,1)=0.d0
    ks_real(1:2,2)=modulo(0.5d0*Sinv(1:2,1)+0.5d0,1.d0)-0.5d0
    ks_real(1:2,3)=modulo(0.5d0*Sinv(1:2,2)+0.5d0,1.d0)-0.5d0
    ks_real(1:2,4)=modulo(0.5d0*(Sinv(1:2,1)+Sinv(1:2,2))+0.5d0,1.d0)-0.5d0
   end select
  endif

! Loop over valid offsets and report k points.
  do iks=1,nks
   ks=ks_real(:,iks)
! Shift to desired offset.
   k_point=k_gamma
   call shift_k_points(ks,npcells,k_point,k_weight)
! Report.
   write(char80,*)dimensionality
   write(6,'(1x,a,'//trim(adjustl(char80))//'(1x,f9.6))')'For k-vector offset:',&
    &ks(1:dimensionality)
   write(6,*)' k points and weights:'
   write(char80,*)dimensionality
   do ik=1,npcells
    if(k_weight(ik)>0)write(6,'(2x,'//trim(adjustl(char80))&
     &//'(1x,es24.16),2x,es24.16)')k_point(1:dimensionality,ik),&
     &dble(k_weight(ik))/dble(npcells)
   enddo ! ik
   write(6,*)
  enddo ! iks

 END SUBROUTINE process_Smat


 SUBROUTINE get_k_points(dimensionality,Smat,nk,k_point)
!-------------------------------------------------------!
! Generate list of NK k points in grid centred at Gamma !
! and congruent with the supercell matrix SMAT.         !
!-------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: dimensionality,Smat(3,3),nk
  REAL(dp),INTENT(out) :: k_point(3,nk)
  INTEGER i,j,k,ik
  REAL(dp) Smat_inv(3,3),dum,kfrac(3)

! Initialize.
  k_point=0.d0
  kfrac=0.d0

! Ivert Smat.
  Smat_inv=inverse(dimensionality,dble(Smat),dum)

  select case(dimensionality)
  case(3)
! Loop over reciprocal lattice sites.
   ik=0
   do i=-(nk+1)/2,(nk+1)/2
    do j=-(nk+1)/2,(nk+1)/2
     do k=-(nk+1)/2,(nk+1)/2
      kfrac=i*Smat_inv(1,:)+j*Smat_inv(2,:)+k*Smat_inv(3,:)
! Reject if point is outside first BZ or within tol. of "top" edge;
! accept if point is inside first BZ or within tol. of "bottom" edge.
      if(any(kfrac>0.5d0.or.are_equal(kfrac,0.5d0)).or.&
       &any(kfrac<-0.5d0.and..not.are_equal(kfrac,-0.5d0)))cycle
! Add this point - take care not to run past end of array.
      ik=ik+1
      if(ik>nk)then
       write(6,*)'Found more k points than should have.'
       stop
      endif
      k_point(:,ik)=kfrac
     enddo ! k
    enddo ! j
   enddo ! i
  case(2)
! Loop over reciprocal lattice sites.
   ik=0
   do i=-(nk+1)/2,(nk+1)/2
    do j=-(nk+1)/2,(nk+1)/2
     kfrac(1:2)=i*Smat_inv(1,1:2)+j*Smat_inv(2,1:2)
! Reject if point is outside first BZ or within tol. of "top" edge;
! accept if point is inside first BZ or within tol. of "bottom" edge.
     if(any(kfrac(1:2)>0.5d0.or.are_equal(kfrac(1:2),0.5d0)).or.&
      &any(kfrac(1:2)<-0.5d0.and..not.are_equal(kfrac(1:2),-0.5d0)))cycle
! Add this point - take care not to run past end of array.
     ik=ik+1
     if(ik>nk)then
      write(6,*)'Found more k points than should have.'
      stop
     endif
     k_point(1:2,ik)=kfrac(1:2)
    enddo ! j
   enddo ! i
  end select ! dimensionality

! Check we have all the k points we should have.
  if(ik<nk)then
   write(6,*)'Found fewer k points than we should have.'
   stop
  endif

 END SUBROUTINE get_k_points


 SUBROUTINE shift_k_points(k_shift,nk,k_point,k_weight)
!--------------------------------------------------------------!
! Shift the NK k-points by K_SHIFT, and compute their weights. !
!--------------------------------------------------------------!
  IMPLICIT NONE
  REAL(dp),INTENT(in) :: k_shift(3)
  INTEGER,INTENT(in) :: nk
  REAL(dp),INTENT(inout) :: k_point(3,nk)
  INTEGER,INTENT(inout) :: k_weight(nk)
  INTEGER ik,jk
  REAL(dp) kvec(3)
  k_weight(1:nk)=1
  do ik=1,nk
   k_point(:,ik)=reduce_to_bz(k_point(:,ik)+k_shift)
   if(use_inversion)then
! See if this point forms a pair with a previous point.
    do jk=1,ik-1
     kvec=reduce_to_bz(k_point(:,ik)+k_point(:,jk))
     if(all(are_equal(kvec,0.d0)))exit
    enddo ! jk
    if(jk<ik)then
! It is, so double the weight of the first and zero that of the second.
     k_weight(jk)=k_weight(jk)+k_weight(ik)
     k_weight(ik)=0
    endif ! jk<ik
   endif ! use_inversion
  enddo ! ik
 END SUBROUTINE shift_k_points


 FUNCTION reduce_to_bz(kvec) RESULT(reduced_kvec)
!-------------------------------------------------------!
! Reduce a k vector expressed in fractional coordinates !
! to the first Brillouin zone.                          !
!-------------------------------------------------------!
  IMPLICIT NONE
  REAL(dp),INTENT(in) :: kvec(3)
  REAL(dp) reduced_kvec(3)
! Move to [-0.5, 0.5].
  reduced_kvec=modulo(kvec+0.5d0,1.d0)-0.5d0
! Move to [-0.5, 0.5).
  where(are_equal(reduced_kvec,0.5d0))reduced_kvec=-0.5d0
 END FUNCTION reduce_to_bz


 REAL(dp) FUNCTION get_wsr(dimensionality,A,B)
!--------------------------------!
! Get WSR for the given A and B. !
!--------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: dimensionality
  REAL(dp),INTENT(in) :: A(3,3),B(3,3)
  INTEGER i1,i2,i3,n1,n2,n3
  REAL(dp) lws2,r2

  select case(dimensionality)

  case(3)
   lws2=min(sum(A(1:3,1)**2),sum(A(1:3,2)**2),sum(A(1:3,3)**2))
   n1=int(sqrt(lws2*sum(B(1:3,1)**2))*inv_twopi+0.01d0)
   n2=int(sqrt(lws2*sum(B(1:3,2)**2))*inv_twopi+0.01d0)
   n3=int(sqrt(lws2*sum(B(1:3,3)**2))*inv_twopi+0.01d0)
   do i1=0,n1
    do i2=-n2,n2
     do i3=-n3,n3
      if(i1==0.and.i2==0.and.i3==0)cycle
      r2=sum((/i1,i2,i3/)*A(1,1:3))**2+sum((/i1,i2,i3/)*A(2,1:3))**2+&
       &sum((/i1,i2,i3/)*A(3,1:3))**2
      if(r2<lws2)lws2=r2
     enddo ! i3
    enddo ! i2
   enddo ! i1
   get_wsr=0.5d0*sqrt(lws2)

  case(2)
   lws2=min(sum(A(1:2,1)**2),sum(A(1:2,2)**2))
   n1=int(sqrt(lws2*sum(B(1:2,1)**2))*inv_twopi+0.01d0)
   n2=int(sqrt(lws2*sum(B(1:2,2)**2))*inv_twopi+0.01d0)
   do i1=0,n1
    do i2=-n2,n2
     if(i1==0.and.i2==0)cycle
     r2=sum((/i1,i2/)*A(1,1:2))**2+sum((/i1,i2/)*A(2,1:2))**2
     if(r2<lws2)lws2=r2
    enddo ! i2
   enddo ! i1
   get_wsr=0.5d0*sqrt(lws2)

  case(1)
   get_wsr=0.5d0*A(1,1)

  end select

 END FUNCTION get_wsr


 SUBROUTINE minkowski_reduce_prod(dimensionality,Smat,Aprim)
!---------------------------------------------------------------!
! Minkowski-reduce the product Smat*Aprim by operating on Smat. !
!---------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: dimensionality
  REAL(dp),INTENT(in) :: Aprim(3,3)
  INTEGER,INTENT(inout) :: Smat(3,3)
  INTEGER Umat(3,10),Mmat(3,10),i,j,i_replace,j_replace,j_x2_min
  REAL(dp) Asim(3,3),Xmat(3,10),x2(10),a2(3),a2_max,x2_min

  select case(dimensionality)
    
  case(3)

! Initialize.
   Umat(1:3,1)=(/1,1,0/)
   Umat(1:3,2)=(/1,-1,0/)
   Umat(1:3,3)=(/1,0,1/)
   Umat(1:3,4)=(/1,0,-1/)
   Umat(1:3,5)=(/0,1,1/)
   Umat(1:3,6)=(/0,1,-1/)
   Umat(1:3,7)=(/1,1,1/)
   Umat(1:3,8)=(/1,1,-1/)
   Umat(1:3,9)=(/1,-1,1/)
   Umat(1:3,10)=(/1,-1,-1/)
   Asim=matmul(Aprim,Smat)
   do i=1,3
    a2(i)=sum(Asim(1:3,i)**2)
   enddo

! Loop over iterations.
   do

! Construct M and X.
    Mmat=matmul(Smat,Umat)
    Xmat=matmul(Aprim,Mmat)
    do i=1,10
     x2(i)=sum(Xmat(1:3,i)**2)
    enddo

! Select which vector to replace.
    i_replace=0
    j_replace=0
    a2_max=0.d0
    do i=1,3
     x2_min=a2(i)
     j_x2_min=0
     do j=1,10
      if(Umat(i,j)/=0)then
! Vector j in Xmat can replace vector i in Asim, since it contains vector i
! (otherwise we would end up with linearly dependent vectors).
       if(x2(j)<x2_min)then
        j_x2_min=j
        x2_min=x2(j)
       endif
      endif
     enddo ! j
     if(j_x2_min/=0)then
! Vector i in Asim can be replaced with j_x2_min in Xmat.  See if
! this is the best overall replacement so far.
      if(a2(i)>a2_max)then
       i_replace=i
       j_replace=j_x2_min
       a2_max=a2(i)
      endif
     endif
    enddo ! i

    if(i_replace==0)exit

! Replace longest vector in Asim with shortest vector in Xmat.
    if(Umat(i_replace,j_replace)<0)then
! Keep handedness of vectors.
     Smat(1:3,i_replace)=-Mmat(1:3,j_replace)
     Asim(1:3,i_replace)=-Xmat(1:3,j_replace)
    else
     Smat(1:3,i_replace)=Mmat(1:3,j_replace)
     Asim(1:3,i_replace)=Xmat(1:3,j_replace)
    endif
    a2(i_replace)=x2(j_replace)

   enddo ! Iteration

  case(2)

! Initialize.
   Asim(1:2,1:2)=matmul(Aprim(1:2,1:2),Smat(1:2,1:2))
   a2(1)=sum(Asim(1:2,1)**2)
   a2(2)=sum(Asim(1:2,2)**2)

! Loop over iterations.
   do
    Mmat(1:2,1)=Smat(1:2,1)+Smat(1:2,2)
    Mmat(1:2,2)=Smat(1:2,1)-Smat(1:2,2)
    Xmat(1:2,1:2)=matmul(Aprim(1:2,1:2),Mmat(1:2,1:2))
    x2(1)=sum(Xmat(1:2,1)**2)
    x2(2)=sum(Xmat(1:2,2)**2)

    if(a2(1)<=a2(2))then
     i_replace=2
    else
     i_replace=1
    endif
    if(x2(1)<=x2(2))then
     j_replace=1
    else
     j_replace=2
    endif
    if(x2(j_replace)>=a2(i_replace))exit

    if(i_replace==2.and.j_replace==2)then
! Keep handedness of vectors.
     Smat(1:2,i_replace)=-Mmat(1:2,j_replace)
     Asim(1:2,i_replace)=-Xmat(1:2,j_replace)
    else
     Smat(1:2,i_replace)=Mmat(1:2,j_replace)
     Asim(1:2,i_replace)=Xmat(1:2,j_replace)
    endif
    a2(i_replace)=x2(j_replace)

   enddo ! Iteration

  case(1)
! Nothing to do in 1D.
   continue
  end select

 END SUBROUTINE minkowski_reduce_prod


 FUNCTION inverse(d,a,det)
!------------------------------------------------------------------------------!
! Compute inverse of 1x1/2x2/3x3 matrix contained in A(1:D,1:D), assuming it's !
! not singular.                                                                !
!------------------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: d
  REAL(dp),INTENT(in) :: a(3,3)
  REAL(dp),INTENT(out) :: det
  REAL(dp) f1,f2,f3,inv_det,inverse(3,3)

! Initialize.
  inverse=0.d0

  select case(d)

  case(3)
! Compute determinant.
   f1=a(2,2)*a(3,3)-a(2,3)*a(3,2)
   f2=a(2,3)*a(3,1)-a(2,1)*a(3,3)
   f3=a(2,1)*a(3,2)-a(2,2)*a(3,1)
   det=a(1,1)*f1+a(1,2)*f2+a(1,3)*f3
   if(.not.are_equal(det,0.d0))then
! Compute inverse.
    inv_det=1.d0/det
    inverse(1:3,1)=(/f1,f2,f3/)*inv_det
    f1=a(1,1)*inv_det
    f2=a(1,2)*inv_det
    f3=a(1,3)*inv_det
    inverse(1:3,2)=(/f3*a(3,2)-f2*a(3,3),f1*a(3,3)-f3*a(3,1),&
     &f2*a(3,1)-f1*a(3,2)/)
    inverse(1:3,3)=(/f2*a(2,3)-f3*a(2,2),f3*a(2,1)-f1*a(2,3),&
     &f1*a(2,2)-f2*a(2,1)/)
   endif

  case(2)
! Compute determinant.
   det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
   if(.not.are_equal(det,0.d0))then
! Compute inverse.
    inv_det=1.d0/det
    inverse(1:2,1)=(/a(2,2),-a(2,1)/)*inv_det
    inverse(1:2,2)=(/-a(1,2),a(1,1)/)*inv_det
   endif

  case(1)
! Compute determinant.
   det=a(1,1)
   if(.not.are_equal(det,0.d0))then
! Compute inverse.
    inverse(1,1)=1.d0/det
   endif

  end select

 END FUNCTION inverse


 LOGICAL ELEMENTAL FUNCTION are_equal(x,y)
!------------------------------------------------------------------------------!
! Check if two floating-point numbers are equal within a reasonable tolerance. !
!------------------------------------------------------------------------------!
  IMPLICIT NONE
  REAL(dp),INTENT(in) :: x,y
  REAL(dp) big,small
! Parameters.
  REAL(dp),PARAMETER :: tol_zero=1.d-12,tol_rel=1.d-12
  if(abs(x)<tol_zero.and.abs(y)<tol_zero)then
   are_equal=.true.
  elseif(x>0.d0.eqv.y>0.d0)then
   big=max(abs(x),abs(y))
   small=min(abs(x),abs(y))
   are_equal=1.d0-small/big<tol_rel
  else
   are_equal=.false.
  endif
 END FUNCTION are_equal


END PROGRAM supercell

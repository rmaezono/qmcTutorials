MODULE heg_tools
  !---------------------------------------!
  ! HEG-related routines and other tools. !
  !---------------------------------------!
  IMPLICIT NONE

  ! Public entities explicitly specified.
  PRIVATE
  ! HEG tools.
  PUBLIC compute_self_term,hf_kinetic_energy,hf_exchange_energy
  ! Lattice generator.
  PUBLIC lattice_max_r,lattice_size,gen_lattice_unshifted,&
     &locate_reoccup_range,shift_lattice,overestimate_wsr,lattice_name
  ! Minimum image.
  PUBLIC min_image
  ! Numerical routines.
  PUBLIC isort_dble_preinit,swap1,inverse,are_equal,reblock_weighted
  ! String tools.
  PUBLIC match_line,extract_word


CONTAINS


  SUBROUTINE compute_self_term(dimensionality,volume,A,Ainv,B,Binv,self_term,&
     &ierr)
    !--------------------------------------------------------!
    ! Evaluate the Ewald self term for a system of the given !
    ! dimensionality, volume, and lattice vectors A and B.   !
    !--------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: dimensionality
    DOUBLE PRECISION,INTENT(in) :: volume,A(3,3),Ainv(3,3),B(3,3),Binv(3,3)
    INTEGER,INTENT(out) :: ierr
    DOUBLE PRECISION,INTENT(out) :: self_term
    DOUBLE PRECISION,ALLOCATABLE :: r_lattice(:,:),r2_lattice(:),&
       &k_lattice(:,:),k2_lattice(:)
    INTEGER,ALLOCATABLE :: rindx_lattice(:),kindx_lattice(:)
    DOUBLE PRECISION a11,rsum,ksum,val_gamma,sqrt_gamma,r2_cutoff,k2_cutoff,r,&
       &r2,k,k2,half_inv_sqrt_gamma,quarter_inv_gamma
    INTEGER i,nr1,nr2,nr3,nr,nk1,nk2,nk3,nk,ialloc
    ! Constants.
    DOUBLE PRECISION,PARAMETER :: pi=3.14159265358979324d0
    DOUBLE PRECISION,PARAMETER :: two_pi=2.d0*pi
    DOUBLE PRECISION,PARAMETER :: four_pi=4.d0*pi
    DOUBLE PRECISION,PARAMETER :: four_pi2=4.d0*pi*pi
    DOUBLE PRECISION,PARAMETER :: one_third=1.d0/3.d0

    ierr=0

    select case(dimensionality)

    case(1)

      ! CASINO default number of real-space points - sum done analytically.
      nr=7
      a11=abs(A(1,1))
      self_term=-log(dble(nr)*a11)/a11-(1.d0/(6.d0*dble(nr**2)*a11))*&
         &(1.d0-0.7d0/dble(nr**2))+11.d0/(6.d0*a11)

    case(2)

      ! CASINO value for gamma and cut-offs.
      sqrt_gamma=2.4d0*volume**(-0.5d0)
      r2_cutoff=24.d0/(sqrt_gamma*sqrt_gamma)
      k2_cutoff=four_pi2*10.d0/volume

      half_inv_sqrt_gamma=1.d0/(2.d0*sqrt_gamma)

      ! Construct r and k lattices.
      call lattice_size(dimensionality,Ainv,sqrt(r2_cutoff),nr1,nr2,nr3,nr)
      call lattice_size(dimensionality,Binv,sqrt(k2_cutoff),nk1,nk2,nk3,nk)
      allocate(r_lattice(3,nr),r2_lattice(nr),rindx_lattice(nr),&
         &k_lattice(3,nk),k2_lattice(nk),kindx_lattice(nk),stat=ialloc)
      if(ialloc/=0)then
        ierr=2
        return
      endif
      call gen_lattice_unshifted(dimensionality,A,nr1,nr2,nr3,r_lattice,&
         &r2_lattice,rindx_lattice)
      call gen_lattice_unshifted(dimensionality,B,nk1,nk2,nk3,k_lattice,&
         &k2_lattice,kindx_lattice)

      ! Real space sum.
      rsum=0.d0
      do i=2,nr
        r2=r2_lattice(rindx_lattice(i))
        if(r2>r2_cutoff)exit
        r=sqrt(r2)
        rsum=rsum+erfc(sqrt_gamma*r)/r
      enddo ! i
      if(i>nr)then
        ierr=1
        return
      endif
      ! Add constant term.
      rsum=rsum-2.d0*sqrt_gamma/sqrt(pi)

      ! Reciprocal space sum.
      ksum=0.d0
      do i=2,nk
        k2=k2_lattice(kindx_lattice(i))
        if(k2>k2_cutoff)exit
        k=sqrt(k2)
        ksum=ksum+erfc(k*half_inv_sqrt_gamma)/k
      enddo ! i
      if(i>nk)then
        ierr=1
        return
      endif
      ! Normalize and add constant term.
      ksum=ksum*two_pi/volume-2.d0*sqrt(pi)/(volume*sqrt_gamma)

      ! Total, halved to account for double-counting.
      self_term=.5d0*(rsum+ksum)

      ! Deallocations.
      deallocate(r_lattice,r2_lattice,rindx_lattice,k_lattice,k2_lattice,&
         &kindx_lattice)

    case(3)

      ! CASINO value for gamma and cut-offs.
      sqrt_gamma=2.8d0*volume**(-one_third)
      r2_cutoff=24.d0/(sqrt_gamma*sqrt_gamma)
      k2_cutoff=four_pi2*12.8d0*volume**(-2.d0/3.d0)

      val_gamma=sqrt_gamma*sqrt_gamma
      quarter_inv_gamma=.25d0/val_gamma

      ! Construct r and k lattices.
      call lattice_size(dimensionality,Ainv,sqrt(r2_cutoff),nr1,nr2,nr3,nr)
      call lattice_size(dimensionality,Binv,sqrt(k2_cutoff),nk1,nk2,nk3,nk)
      allocate(r_lattice(3,nr),r2_lattice(nr),rindx_lattice(nr),&
         &k_lattice(3,nk),k2_lattice(nk),kindx_lattice(nk),stat=ialloc)
      if(ialloc/=0)then
        ierr=2
        return
      endif
      call gen_lattice_unshifted(dimensionality,A,nr1,nr2,nr3,r_lattice,&
         &r2_lattice,rindx_lattice)
      call gen_lattice_unshifted(dimensionality,B,nk1,nk2,nk3,k_lattice,&
         &k2_lattice,kindx_lattice)

      ! Real space sum.
      rsum=0.d0
      do i=2,nr
        r2=r2_lattice(rindx_lattice(i))
        if(r2>r2_cutoff)exit
        r=sqrt(r2)
        rsum=rsum+erfc(sqrt_gamma*r)/r
      enddo ! i
      if(i>nr)then
        ierr=1
        return
      endif
      ! Add constant term.
      rsum=rsum-2.d0*sqrt_gamma/sqrt(pi)

      ! Reciprocal space sum.
      ksum=0.d0
      do i=2,nk
        k2=k2_lattice(kindx_lattice(i))
        if(k2>k2_cutoff)exit
        ksum=ksum+exp(-k2*quarter_inv_gamma)/k2
      enddo ! i
      if(i>nk)then
        ierr=1
        return
      endif
      ! Normalize and add constant term.
      ksum=ksum*four_pi/volume-pi/(val_gamma*volume)

      ! Total, halved to account for double-counting.
      self_term=.5d0*(rsum+ksum)

      ! Deallocations.
      deallocate(r_lattice,r2_lattice,rindx_lattice,k_lattice,k2_lattice,&
         &kindx_lattice)

    end select ! dimensionality

  END SUBROUTINE compute_self_term


  DOUBLE PRECISION FUNCTION hf_kinetic_energy(nspin,nele,inv_pmass,inv_netot,&
     &k2_lattice,kindx_lattice)
    !--------------------------------------------------------!
    ! Compute the HF kinetic energy of the configured system !
    ! with the provided set of k-lattice vectors.            !
    !--------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nspin,nele(nspin),kindx_lattice(*)
    DOUBLE PRECISION,INTENT(in) :: k2_lattice(*),inv_pmass(nspin),inv_netot
    INTEGER ispin,ie
    DOUBLE PRECISION t1,t2

    t1=0.d0
    do ispin=1,nspin
      t2=0.d0
      do ie=1,nele(ispin)
        t2=t2+k2_lattice(kindx_lattice(ie))
      enddo ! ie
      t1=t1+t2*inv_pmass(ispin)
    enddo ! ispin
    hf_kinetic_energy=0.5d0*t1*inv_netot

  END FUNCTION hf_kinetic_energy


  DOUBLE PRECISION FUNCTION hf_exchange_energy(dimensionality,nspin,nele,&
     &pcharge2,inv_netot_volume,self_term,k_lattice,kindx_lattice)
    !---------------------------------------------------------!
    ! Compute the HF exchange energy of the configured system !
    ! with the provided set of k-lattice vectors.             !
    !---------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: dimensionality,nspin,nele(nspin),kindx_lattice(*)
    DOUBLE PRECISION,INTENT(in) :: pcharge2(nspin),inv_netot_volume,self_term,&
       &k_lattice(3,*)
    INTEGER ispin,ie,je
    DOUBLE PRECISION t1,t2,kdiff(3),kdiff2,ki(3)

    t1=0.d0
    do ispin=1,nspin
      t2=0.d0
      do ie=1,nele(ispin)-1
        ki=k_lattice(:,kindx_lattice(ie))
        do je=ie+1,nele(ispin)
          kdiff=ki-k_lattice(:,kindx_lattice(je))
          kdiff2=dot_product(kdiff,kdiff)
          t2=t2+coulomb_transform(dimensionality,kdiff2)
        enddo ! je
      enddo ! ie
      t1=t1+t2*pcharge2(ispin)
    enddo ! ispin
    hf_exchange_energy=-t1*inv_netot_volume+self_term

  END FUNCTION hf_exchange_energy


  DOUBLE PRECISION FUNCTION coulomb_transform(dimensionality,k2)
    !------------------------------------------------------!
    ! Evaluate the Coulomb transform of 1/r ar k=sqrt(k2). !
    !------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: dimensionality
    DOUBLE PRECISION,INTENT(in) :: k2
    ! Constants.
    DOUBLE PRECISION,PARAMETER :: pi=3.14159265358979324d0
    DOUBLE PRECISION,PARAMETER :: two_pi=2.d0*pi
    DOUBLE PRECISION,PARAMETER :: four_pi=4.d0*pi
    DOUBLE PRECISION,PARAMETER :: euler_constant=0.577215664901532860606d0

    if(abs(k2)>0.d0)then
      select case(dimensionality)
      case(3)
        coulomb_transform=four_pi/k2
      case(2)
        coulomb_transform=two_pi/sqrt(k2)
      case(1)
        coulomb_transform=-log(0.25d0*k2)-2.d0*euler_constant
      end select
    else
      coulomb_transform=0.d0
    endif

  END FUNCTION coulomb_transform


  DOUBLE PRECISION FUNCTION lattice_max_r(dimensionality,A,n)
    !---------------------------------------------------------------!
    ! Get radius of sphere/circle/segment that holds the N shortest !
    ! lattice points of a cubic/square/linear simulation cell of    !
    ! lattice constant equal to the length of the largest A(:,i)    !
    ! vector.  This serves as an upper bound to the radius of the   !
    ! sphere that contains the N shortest vectors in the cell       !
    ! described by A.                                               !
    !---------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: dimensionality,n
    DOUBLE PRECISION,INTENT(in) :: A(3,3)
    INTEGER i,ir
    DOUBLE PRECISION d,mod_a(3)
    ! Constants.
    DOUBLE PRECISION,PARAMETER :: pi=3.14159265358979324d0
    DOUBLE PRECISION,PARAMETER :: three_over_four_pi=3.d0/(4.d0*pi)
    DOUBLE PRECISION,PARAMETER :: one_third=1.d0/3.d0

    ! Get length of largest b vector.
    do i=1,dimensionality
      mod_a(i)=sqrt(sum(A(1:dimensionality,i)**2))
    enddo ! i
    d=maxval(mod_a(1:dimensionality))

    ! Compute integer prefactor.
    select case(dimensionality)
    case(3)
      ir=ceiling((three_over_four_pi*dble(n))**one_third+0.25d0)
    case(2)
      ! FIXME - adjust
      ir=ceiling((dble(n)/pi)**.5d0+1.d0)
    case(1)
      ir=n/2
    end select

    ! Result.
    lattice_max_r=d*dble(ir)

  END FUNCTION lattice_max_r


  SUBROUTINE lattice_size(dimensionality,Ainv,r,n1,n2,n3,n)
    !-------------------------------------------------------------------!
    ! Construct parallelepiped of reciprocal-space cells that contain a !
    ! sphere of radius R, and return its dimensions N1,N2,N3 and total  !
    ! number of lattice points N.  The latter should be used as an      !
    ! allocation size for the lattice vector array.                     !
    !-------------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: dimensionality
    DOUBLE PRECISION,INTENT(in) :: r,Ainv(3,3)
    INTEGER,INTENT(out) :: n1,n2,n3,n
    DOUBLE PRECISION perp_a(3)
    INTEGER i

    ! Compute dimensions of parallelepiped.
    do i=1,dimensionality
      perp_a(i)=sqrt(sum(Ainv(i,1:dimensionality)**2))
    enddo ! i
    n1=ceiling(r*perp_a(1))
    n2=0
    if(dimensionality>1)n2=ceiling(r*perp_a(2))
    n3=0
    if(dimensionality>2)n3=ceiling(r*perp_a(3))

    ! Total number of vectors.
    n=(2*n1+1)*(2*n2+1)*(2*n3+1)

  END SUBROUTINE lattice_size


  SUBROUTINE gen_lattice_unshifted(dimensionality,A,n1,n2,n3,r_lattice,&
     &r2_lattice,rindx_lattice)
    !------------------------------------------------------!
    ! Generate lattice of N1xN2xN3 vectors in the geometry !
    ! specified by A.                                      !
    !------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: dimensionality,n1,n2,n3
    DOUBLE PRECISION,INTENT(in) :: A(3,3)
    DOUBLE PRECISION,INTENT(inout) :: r_lattice(3,*),r2_lattice(*)
    INTEGER,INTENT(inout) :: rindx_lattice(*)
    INTEGER nvec,i1,i2,i3,i
    DOUBLE PRECISION rvec(3),r2,i1a1(3),i2a2(3),i1a1_i2a2(3)

    ! Initialize counter.
    nvec=0

    ! Add Gamma point.
    nvec=nvec+1
    r_lattice(:,nvec)=0.d0
    r2_lattice(nvec)=0.d0

    if(dimensionality>2)then
      ! Add vectors of the form (0,0,|i3|).(b1,b2,b3)
      rvec=A(:,3)
      do i3=1,n3
        r2=sum(rvec*rvec)
        r_lattice(:,nvec+1)=rvec
        r2_lattice(nvec+1)=r2
        r_lattice(:,nvec+2)=-rvec
        r2_lattice(nvec+2)=r2
        nvec=nvec+2
        rvec=rvec+A(:,3)
      enddo ! i3
    endif

    if(dimensionality>1)then
      ! Add vectors of the form (0,|i2|,i3).(0,b2,b3)
      i2a2=A(:,2)
      do i2=1,n2
        rvec=i2a2-n3*A(:,3)
        do i3=-n3,n3
          r2=sum(rvec*rvec)
          r_lattice(:,nvec+1)=rvec
          r2_lattice(nvec+1)=r2
          r_lattice(:,nvec+2)=-rvec
          r2_lattice(nvec+2)=r2
          nvec=nvec+2
          rvec=rvec+A(:,3)
        enddo ! i3
        i2a2=i2a2+A(:,2)
      enddo ! i2
    endif

    ! Add vectors of the form (|i1|,i2,i3).(0,b2,b3)
    i1a1=A(:,1)
    do i1=1,n1
      i1a1_i2a2=i1a1-n2*A(:,2)
      do i2=-n2,n2
        rvec=i1a1_i2a2-n3*A(:,3)
        do i3=-n3,n3
          r2=sum(rvec*rvec)
          r_lattice(:,nvec+1)=rvec
          r2_lattice(nvec+1)=r2
          r_lattice(:,nvec+2)=-rvec
          r2_lattice(nvec+2)=r2
          nvec=nvec+2
          rvec=rvec+A(:,3)
        enddo ! i3
        i1a1_i2a2=i1a1_i2a2+A(:,2)
      enddo ! i2
      i1a1=i1a1+A(:,1)
    enddo ! i1

    ! Initialize indices for sorting.
    do i=1,nvec
      rindx_lattice(i)=i
    enddo ! i

    ! Sort by length in ascending order using a stable algorithm.
    call isort_dble_preinit(nvec,1,nvec,r2_lattice,rindx_lattice)

  END SUBROUTINE gen_lattice_unshifted


  SUBROUTINE locate_reoccup_range(n,n_occup,gamma_k2_lattice,&
     &gamma_kindx_lattice,max_mod_k_offset,kindx_min,kindx_max)
    !-----------------------------------------------------------------!
    ! Return the indices of the first and last k vectors in the range !
    ! where the occupancy needs to be recalculated when the k-offset  !
    ! changes.                                                        !
    !-----------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n,n_occup,gamma_kindx_lattice(n)
    DOUBLE PRECISION,INTENT(in) :: gamma_k2_lattice(n),max_mod_k_offset
    INTEGER,INTENT(out) :: kindx_min,kindx_max
    DOUBLE PRECISION kF,k2min,k2max
    INTEGER i,ikF
    ikF=n_occup
    kF=sqrt(gamma_k2_lattice(gamma_kindx_lattice(ikF)))
    k2min=max(kF-max_mod_k_offset,0.d0)**2
    k2max=(kF+max_mod_k_offset)**2
    kindx_min=ikF
    do i=ikF-1,1,-1
      if(gamma_k2_lattice(gamma_kindx_lattice(i))<k2min)exit
      kindx_min=i
    enddo ! i
    kindx_max=ikF
    do i=ikF+1,n
      if(gamma_k2_lattice(gamma_kindx_lattice(i))>k2max)exit
      kindx_max=i
    enddo ! i
  END SUBROUTINE locate_reoccup_range


  SUBROUTINE shift_lattice(n,n1,n2,gamma_k_lattice,gamma_kindx_lattice,&
     &k_offset,k_lattice,k2_lattice,kindx_lattice)
    !----------------------------------!
    ! Shift gamma lattice by k_offset. !
    !----------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n,n1,n2,gamma_kindx_lattice(n)
    DOUBLE PRECISION,INTENT(in) :: gamma_k_lattice(3,n),k_offset(3)
    DOUBLE PRECISION,INTENT(inout) :: k_lattice(3,n),k2_lattice(n)
    INTEGER,INTENT(inout) :: kindx_lattice(n)
    INTEGER i,j

    ! Apply the k_offset.
    do i=1,n2
      j=gamma_kindx_lattice(i)
      k_lattice(1:3,j)=gamma_k_lattice(1:3,j)+k_offset(1:3)
      k2_lattice(j)=dot_product(k_lattice(:,j),k_lattice(:,j))
    enddo ! i

    ! Sort, starting from the original ordering, using an algorithm that is
    ! fast for nearly-sorted sequences.
    kindx_lattice=gamma_kindx_lattice
    call isort_dble_preinit(n,n1,n2,k2_lattice,kindx_lattice)

  END SUBROUTINE shift_lattice


  DOUBLE PRECISION FUNCTION overestimate_wsr(dimensionality,A)
    !------------------------------------------------------------------!
    ! Return an overestimate of the radius of the sphere that contains !
    ! the Wigner-Seitz cell of the lattice defined by A.               !
    !------------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: dimensionality
    DOUBLE PRECISION,INTENT(in) :: A(3,3)
    DOUBLE PRECISION diag1_2,diag2_2,diag3_2,diag4_2
    select case(dimensionality)
    case(3)
      diag1_2=sum((A(:,1)+A(:,2)+A(:,3))**2)
      diag2_2=sum((A(:,1)-A(:,2)+A(:,3))**2)
      diag3_2=sum((A(:,1)+A(:,2)-A(:,3))**2)
      diag4_2=sum((A(:,1)-A(:,2)-A(:,3))**2)
      overestimate_wsr=0.5d0*sqrt(max(diag1_2,diag2_2,diag3_2,diag4_2))
    case(2)
      diag1_2=sum((A(:,1)+A(:,2))**2)
      diag2_2=sum((A(:,1)-A(:,2))**2)
      overestimate_wsr=0.5d0*sqrt(max(diag1_2,diag2_2))
    case(1)
      overestimate_wsr=0.5d0*abs(A(1,1))
    end select
  END FUNCTION overestimate_wsr


  CHARACTER(20) FUNCTION lattice_name(dimensionality,A)
    !----------------------------------------------------------!
    ! Return lattice name in a string.  Only recognizes simple !
    ! lattices, and returns "unknown" for others.              !
    !----------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: dimensionality
    DOUBLE PRECISION,INTENT(in) :: A(3,3)
    INTEGER i
    DOUBLE PRECISION mod_a(3),cos_alpha(3)
    ! Constants.
    DOUBLE PRECISION,PARAMETER :: one_third=1.d0/3.d0

    ! Initialize.
    lattice_name=''

    ! Compute length of lattice vectors.
    do i=1,dimensionality
      mod_a(i)=sqrt(sum(A(1:dimensionality,i)**2))
    enddo ! i

    select case(dimensionality)

    case(3)
      ! 3D lattices recognized: sc, fcc, bcc.
      if(are_equal(mod_a(1),mod_a(2)).and.are_equal(mod_a(1),mod_a(3)))then
        cos_alpha(1)=dot_product(A(1:3,2),A(1:3,3))/(mod_a(2)*mod_a(3))
        cos_alpha(2)=dot_product(A(1:3,3),A(1:3,1))/(mod_a(3)*mod_a(1))
        cos_alpha(3)=dot_product(A(1:3,1),A(1:3,2))/(mod_a(1)*mod_a(2))
        if(are_equal(cos_alpha(1),0.d0).and.are_equal(cos_alpha(2),0.d0)&
           &.and.are_equal(cos_alpha(3),0.d0))then
          lattice_name='simple cubic'
        elseif(are_equal(cos_alpha(1),0.5d0).and.&
           &are_equal(cos_alpha(2),0.5d0).and.&
           &are_equal(cos_alpha(3),0.5d0))then
          lattice_name='fcc'
        elseif(are_equal(cos_alpha(1),-one_third).and.&
           &are_equal(cos_alpha(2),-one_third).and.&
           &are_equal(cos_alpha(3),-one_third))then
          lattice_name='bcc'
        endif
      endif
      if(len_trim(lattice_name)==0)lattice_name='unknown 3D'

    case(2)
      ! 2D lattices recognized: square, rectangular, hexagonal.
      cos_alpha(3)=dot_product(A(1:2,1),A(1:2,2))/(mod_a(1)*mod_a(2))
      if(are_equal(cos_alpha(3),0.d0))then
        if(are_equal(mod_a(1),mod_a(2)))then
          lattice_name='square'
        else
          lattice_name='rectangular'
        endif
      elseif(are_equal(mod_a(1),mod_a(2)))then
        if(are_equal(abs(cos_alpha(3)),0.5d0))lattice_name='hexagonal'
      endif
      if(len_trim(lattice_name)==0)lattice_name='unknown 2D'

    case(1)
      lattice_name='linear'

    end select

  END FUNCTION lattice_name


  SUBROUTINE min_image(dimensionality,Amat,Ainv,vec,ierr)
    !----------------------------------------------------------!
    ! Return the shortest periodic image of vec in the lattice !
    ! defined by AMAT.                                         !
    !----------------------------------------------------------!
    INTEGER,INTENT(in) :: dimensionality
    DOUBLE PRECISION,INTENT(in) :: Amat(3,3),Ainv(3,3)
    DOUBLE PRECISION,INTENT(inout) :: vec(3)
    INTEGER,INTENT(out) :: ierr
    INTEGER i,j,k,l,im,subcell_indx(dimensionality),m_max(3),min_im,&
       &mlist(3,8),size_ext,ialloc
    DOUBLE PRECISION dist2,outvec(3),mod2_Ainv_i,min_dist2,mod_outvec,&
       &rel_vec(3,8)
    INTEGER,SAVE :: saved_size_ext=0
    INTEGER,ALLOCATABLE,SAVE :: mlist_ext(:,:)
    DOUBLE PRECISION,ALLOCATABLE,SAVE :: rel_vec_ext(:,:)

    ierr=0
    outvec=0.d0

    ! Get lattice offset.
    do i=1,dimensionality
      subcell_indx(i)=floor(sum(vec(1:dimensionality)*&
         &Ainv(1:dimensionality,i)))
    enddo ! i

    ! Allocate work arrays and construct list of offsets to get to all
    ! vertices of the cell.
    select case(dimensionality)
    case(3)
      l=0
      do k=0,1
        do j=0,1
          do i=0,1
            l=l+1
            mlist(:,l)=(/i,j,k/)
          enddo ! i
        enddo ! j
      enddo ! k
    case(2)
      l=0
      do j=0,1
        do i=0,1
          l=l+1
          mlist(1:2,l)=(/i,j/)
        enddo ! i
      enddo ! j
    case(1)
      l=0
      do i=0,1
        l=l+1
        mlist(1,l)=i
      enddo ! i
    end select

    ! Compute relative vectors from all vertices of the cell.
    min_im=0
    min_dist2=0.d0
    do im=1,l
      do i=1,dimensionality
        rel_vec(i,im)=vec(i)-sum(dble(subcell_indx(1:dimensionality)+&
           &mlist(1:dimensionality,im))*Amat(i,1:dimensionality))
      enddo ! i
      dist2=sum(rel_vec(1:dimensionality,im)**2)
      if(min_im==0.or.min_dist2>dist2)then
        min_im=im
        min_dist2=dist2
      endif
    enddo ! im

    ! Locate minimum distance.
    outvec(1:dimensionality)=rel_vec(1:dimensionality,min_im)

    if(dimensionality==3)then

      ! The single-cell vertex algorithm above is not guaranteed to find the
      ! minimum image in 3 dimensions, so one needs to extend the search to a
      ! number of contiguous cells.

      ! Compute number of contiguous cells to extend the search to in the
      ! direction of each of the lattice vectors.
      mod_outvec=sqrt(min_dist2)
      do i=1,3
        mod2_Ainv_i=sum(Ainv(1:3,i)**2)
        m_max(i)=floor(mod_outvec*mod2_Ainv_i)
      enddo ! i

      ! Proceed only if any additional cells are to be included.
      if(.not.all(m_max==0))then

        ! Allocate work arrays.
        size_ext=product(2*m_max(:)+2)
        if(saved_size_ext<size_ext)then
          if(saved_size_ext>0)deallocate(mlist_ext,rel_vec_ext)
          allocate(mlist_ext(3,size_ext),rel_vec_ext(3,size_ext),stat=ialloc)
          if(ialloc/=0)then
            ierr=2
            return
          endif
          saved_size_ext=size_ext
        endif

        ! Construct list of offsets to get to all relevant vertices.
        l=0
        do k=-m_max(3),m_max(3)+1
          do j=-m_max(2),m_max(2)+1
            do i=-m_max(1),m_max(1)+1
              l=l+1
              mlist_ext(1:3,l)=(/i,j,k/)
            enddo ! i
          enddo ! j
        enddo ! k

        ! Compute relative vectors from all vertices of the cell.
        min_im=0
        do im=1,l
          do i=1,dimensionality
            rel_vec_ext(i,im)=vec(i)-sum(dble(subcell_indx(1:dimensionality)+&
               &mlist_ext(1:dimensionality,im))*Amat(i,1:dimensionality))
          enddo ! i
          dist2=sum(rel_vec_ext(1:dimensionality,im)**2)
          if(min_dist2>dist2)then
            min_im=im
            min_dist2=dist2
          endif
        enddo ! im

        ! Replace result with minimum image vector.
        if(min_im>0)outvec(1:dimensionality)=rel_vec_ext(:,min_im)

      endif ! there are additional cells to check

    endif ! 3D

    ! Return output vector.
    vec=outvec

  END SUBROUTINE min_image


  SUBROUTINE isort_dble_preinit(n,n1,n2,x,indx)
    !-----------------------------------------------------!
    ! Perform insertion sort on a real vector X(1:N) with !
    ! pre-initialized INDX.  This sorting algorithm       !
    ! typically costs ~ N^2, but is stable (preserves the !
    ! order of entries with same X) and becomes order N   !
    ! when X(INDX(:)) is nearly sorted.                   !
    !-----------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n,n1,n2
    DOUBLE PRECISION,INTENT(in) :: x(n)
    INTEGER,INTENT(inout) :: indx(n)
    INTEGER i,j
    DOUBLE PRECISION xi
    ! Loop over elements from n1+1:n2.
    do i=n1+1,n2
      xi=x(indx(i))
      ! Move i-th element upwards until it is greater than previous,
      ! at which point the first I-th elements will be sorted.
      do j=i-1,n1,-1
        if(xi>=x(indx(j)))exit
        call swap1(indx(j),indx(j+1))
      enddo !j
    enddo ! i
  END SUBROUTINE isort_dble_preinit


  SUBROUTINE swap1(x,y)
    !--------------------!
    ! Swap two integers. !
    !--------------------!
    IMPLICIT NONE
    INTEGER,INTENT(inout) :: x,y
    INTEGER z
    z=x
    x=y
    y=z
  END SUBROUTINE swap1


  FUNCTION inverse(d,a,det)
    !----------------------------------------------------------------!
    ! Compute inverse of 1x1/2x2/3x3 matrix contained in A(1:D,1:D), !
    ! assuming it's not singular.                                    !
    !----------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: d
    DOUBLE PRECISION,INTENT(in) :: a(3,3)
    DOUBLE PRECISION,INTENT(out) :: det
    DOUBLE PRECISION f1,f2,f3,inv_det
    DOUBLE PRECISION :: inverse(3,3)
    inverse=0.d0
    select case(d)
    case(3)
      f1=a(2,2)*a(3,3)-a(2,3)*a(3,2)
      f2=a(2,3)*a(3,1)-a(2,1)*a(3,3)
      f3=a(2,1)*a(3,2)-a(2,2)*a(3,1)
      det=a(1,1)*f1+a(1,2)*f2+a(1,3)*f3
      if(.not.are_equal(det,0.d0))then
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
      det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      if(.not.are_equal(det,0.d0))then
        inv_det=1.d0/det
        inverse(1:2,1)=(/a(2,2),-a(2,1)/)*inv_det
        inverse(1:2,2)=(/-a(1,2),a(1,1)/)*inv_det
      endif
    case(1)
      det=a(1,1)
      inverse(1,1)=1.d0/det
    end select
  END FUNCTION inverse


  LOGICAL FUNCTION are_equal(x,y)
    !------------------------------------------------------!
    ! Check if two floating-point numbers are equal within !
    ! a reasonable tolerance.                              !
    !------------------------------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(in) :: x,y
    DOUBLE PRECISION big,small
    ! Parameters.
    DOUBLE PRECISION,PARAMETER :: tol_zero=1.d-50
    DOUBLE PRECISION,PARAMETER :: tol_rel=1.d-12
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


  SUBROUTINE reblock_weighted(m,e,w,mean,err_in_mean)
    !------------------------------------------------------------!
    ! Given a series of M data E(1:M) and weights W(1:M), return !
    ! the mean and reblocked errorbar.                           !
    !------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: m
    DOUBLE PRECISION,INTENT(in) :: e(m),w(m)
    DOUBLE PRECISION,INTENT(out) :: mean,err_in_mean
    INTEGER nrtn,rtn,block_m,nblock,k,i
    DOUBLE PRECISION x(m),sum_x,sum_w,ave_x,var_e,sum_block_w2,block_sum_w,&
       &block_ave_x,eff_nblock,err_e,ncorr,cc_err_e,&
       &err_e_vector(floor(log(dble(m))/log(2.d0)))

    ! Basic quantities.
    x=e*w
    sum_x=sum(x)
    sum_w=sum(w)
    ave_x=sum_x/sum_w
    
    ! Loop over reblocking transformation numbers (RTNs).
    nrtn=floor(log(dble(m))/log(2.d0))
    block_m=1
    do rtn=1,nrtn

      ! Number of blocks at this RTN.
      nblock=int(m/block_m)
      
      ! Evaluate the sum of the squares of the deviations from the average.
      ! Last, incomplete block has fewer data points and hence a smaller
      ! weight.
      var_e=0.d0
      sum_block_w2=0.d0
      k=0
      do i=1,nblock
        block_sum_w=sum(w(k+1:k+block_m))
        block_ave_x=sum(x(k+1:k+block_m))/block_sum_w
        var_e=var_e+block_sum_w*(block_ave_x-ave_x)**2
        sum_block_w2=sum_block_w2+block_sum_w**2
        k=k+block_m
      enddo ! i
      eff_nblock=dble(nblock)
      if(m>k)then
        block_sum_w=sum(w(k+1:m))
        block_ave_x=sum(x(k+1:m))/block_sum_w
        var_e=var_e+block_sum_w*(block_ave_x-ave_x)**2
        sum_block_w2=sum_block_w2+block_sum_w**2
        eff_nblock=eff_nblock+dble(m-k)/dble(block_m)
      endif

      ! Evaluate variance, standard error in mean and error in standard error.
      var_e=var_e/(sum_w-sum_block_w2/sum_w)
      err_e=sqrt(var_e/eff_nblock)
      err_e_vector(rtn)=err_e

      ! Double block length for next reblock.
      block_m=block_m*2

    enddo ! rtn

    ! Analyze reblock plot to obtain correlation-corrected errorbar.
    cc_err_e=maxval(err_e_vector)
    if(err_e_vector(1)>0.d0)then
      block_m=1
      do rtn=1,nrtn
        ncorr=(err_e_vector(rtn)/err_e_vector(1))**2
        if(dble(block_m**3)>=2.d0*dble(m)*ncorr**2)then
          cc_err_e=err_e_vector(rtn)
          exit
        endif
        block_m=block_m*2
      enddo ! rtn
    endif

    ! Return result.
    mean=ave_x
    err_in_mean=cc_err_e

  END SUBROUTINE reblock_weighted


  LOGICAL FUNCTION match_line(line,words)
    !----------------------------------------------------------!
    ! Returns .true. if the initial words in LINE are the same !
    ! as the words in WORDS.                                   !
    !----------------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: line,words
    CHARACTER(len(line)) tline,word1
    CHARACTER(len(words)) twords,word2
    match_line=.true.
    tline=adjustl(line)
    twords=adjustl(words)
    do
      ! Exit if either string is empty.
      if(len_trim(tline)==0)then
        match_line=len_trim(twords)==0
        exit
      elseif(len_trim(twords)==0)then
        match_line=.true.
        exit
      endif
      ! Get next word in tline and twords.
      word1=extract_word(tline,1)
      word2=extract_word(twords,1)
      ! Compare the words.
      if(trim(word1)/=trim(word2))then
        match_line=.false.
        exit
      endif
      ! Remove first word.
      tline=adjustl(tline(len_trim(word1)+1:))
      twords=adjustl(twords(len_trim(word2)+1:))
    enddo
  END FUNCTION match_line


  FUNCTION extract_word(line,iword)
    !------------------------------------!
    ! Returns the IWORD-th word in LINE. !
    !------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: line
    INTEGER,INTENT(in) :: iword
    CHARACTER(len(line)) :: extract_word
    CHARACTER(len(line)) tline
    INTEGER i,ipos
    extract_word=''
    if(iword>0)then
      tline=adjustl(line)
      do i=1,iword-1
        if(len_trim(tline)==0)return
        ipos=scan(tline,' ')
        if(ipos<1)return
        tline=adjustl(tline(ipos+1:))
      enddo ! i
      if(len_trim(tline)==0)return
      ipos=scan(tline,' ')
      if(ipos>0)then
        extract_word=tline(1:ipos-1)
      else
        extract_word=tline
      endif
    elseif(iword<0)then
      tline=adjustr(line)
      do i=1,-iword-1
        if(len_trim(tline)==0)return
        ipos=scan(tline,' ',back=.true.)
        if(ipos<1)return
        ! NB, following must go in two steps, else adjustr operates on string
        ! of length L-IPOS [!].
        tline=tline(:ipos-1)
        tline=adjustr(tline)
      enddo ! i
      if(len_trim(tline)==0)return
      ipos=scan(tline,' ',back=.true.)
      if(ipos>0)then
        extract_word=tline(ipos+1:)
      else
        extract_word=tline
      endif
    endif
  END FUNCTION extract_word


  DOUBLE PRECISION FUNCTION erfc(x)
!-------------------------------------------------------------!
! This function returns an approximation to the complementary !
! error function.  In particular, the fit proposed by P. van  !
! Halen, Elec. Lett. 25, 561 (1989) is used.  The relative    !
! error is less that 2.8E-9 everywhere and the absolute error !
! is less than 1.6E-9 everywhere.                             !
!-------------------------------------------------------------!
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(in) :: x
  DOUBLE PRECISION,PARAMETER :: A1=-1.1283791670955125739d0,&
   &A2=-6.366197121956d-1,&
   &A3=-1.027728162129d-1,A4=1.912427299414d-2,A5=2.401479235527d-4,&
   &A6=-1.786242904258d-3,A7=7.336113173091d-4,A8=-1.655799102866d-4,&
   &A9=2.116490536557d-5,A10=-1.196623630319d-6
  DOUBLE PRECISION t
  t=abs(x)
  erfc=&
   &exp(t*(A1+t*(A2+t*(A3+t*(A4+t*(A5+t*(A6+t*(A7+t*(A8+t*(A9+t*A10))))))))))
  if(x<0.d0)erfc=2.d0-erfc
  END FUNCTION erfc


END MODULE heg_tools

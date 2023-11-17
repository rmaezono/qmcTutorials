MODULE bulirsch_stoer
 USE dsp
 USE run_control, only : errstop_master
 USE store, only : dim
 IMPLICIT NONE
 PRIVATE
 PUBLIC stepper_bs


CONTAINS


 SUBROUTINE stepper_bs(x,vel,t,htry,eps,xscal,hdid,hnext,velocity)
  IMPLICIT NONE
  REAL(dp),INTENT(inout) :: x(dim)
  REAL(dp),INTENT(inout) :: t
  REAL(dp),INTENT(in) :: vel(dim),xscal(dim)
  REAL(dp),INTENT(in) :: htry,eps
  REAL(dp),INTENT(out) :: hdid,hnext

  INTERFACE
   SUBROUTINE velocity(t,x,vel)
    USE dsp
    USE store, ONLY : dim
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: t
    REAL(dp),INTENT(in) :: x(dim)
    REAL(dp),INTENT(out) :: vel(dim)
   END SUBROUTINE velocity
  END INTERFACE

  INTEGER,PARAMETER :: imax=9,kmaxx=imax-1
  REAL(dp),PARAMETER :: safe1=0.25d0,safe2=0.7d0,redmax=1.d-5,&
   &redmin=0.7d0,tiny=1.d-30,scalmx=0.1d0
  INTEGER k,km
  INTEGER,DIMENSION(imax) :: nseq=(/2,4,6,8,10,12,14,16,18/)
  INTEGER,SAVE :: kopt,kmax
  REAL(dp),SAVE :: alf(kmaxx,kmaxx)
  REAL(dp) err(kmaxx)
  REAL(dp),SAVE :: a(imax)
  REAL(dp),SAVE :: epsold=-1.d0,tnew
  REAL(dp) eps1,errmax,fact,h,red,scale,wrkmin,t_est
  REAL(dp) xerr(dim),xsav(dim),xseq(dim)
  LOGICAL reduct
  LOGICAL,SAVE :: first=.true.

  if(eps/=epsold)then
   hnext=-1.d29
   tnew=-1.d29
   eps1=safe1*eps
   a(:)=cumsum(nseq,1)
   where(upper_triangle(kmaxx,kmaxx))alf=eps1**&
    &(outerdiff(a(2:),a(2:))/outerprod(arth(3.d0,2.d0,kmaxx),(a(2:)-a(1)+1.d0)))
   epsold=eps
   do kopt=2,kmaxx-1
    if(a(kopt+1)>a(kopt)*alf(kopt-1,kopt))exit
   enddo
   kmax=kopt
  endif
  h=htry
  xsav(:)=x(:)
  if(h/=hnext.or.t/=tnew)then
   first=.true.
   kopt=kmax
  endif
  reduct=.false.

  main_loop: do
   do k=1,kmax
    tnew=t+h
    if(tnew==t)call errstop_master('STEPPER_BS','Step size underflow.')
    call mmid(xsav,vel,t,h,nseq(k),xseq,velocity)
    t_est=(h/real(nseq(k),dp))**2
    call pzextr(k,t_est,xseq,x,xerr)
    if(k/=1)then
     errmax=maxval(abs(xerr(:)/xscal(:)))
     errmax=max(tiny,errmax)/eps
     km=k-1
     err(km)=(errmax/safe1)**(1.d0/(2*km+1))
    endif
    if(k/=1.and.(k>=kopt-1.or.first))then
     if(errmax<1.d0)exit main_loop
     if(k==kmax.or.k==kopt+1)then
      red=safe2/err(km)
      exit
     elseif(k==kopt)then
      if(alf(kopt-1,kopt)<err(km))then
       red=1.d0/err(km)
       exit
      endif
     elseif(kopt==kmax)then
      if(alf(km,kmax-1)<err(km))then
       red=alf(km,kmax-1)*safe2/err(km)
       exit
      endif
     elseif(alf(km,kopt)<err(km))then
      red=alf(km,kopt-1)/err(km)
      exit
     endif
    endif
   enddo
   red=max(min(red,redmin),redmax)
   h=h*red
   reduct=.true.
  enddo main_loop

  t=tnew
  hdid=h
  first=.false.
  kopt=1+iminloc(a(2:km+1)*max(err(1:km),scalmx))
  scale=max(err(kopt-1),scalmx)
  wrkmin=scale*a(kopt)
  hnext=h/scale
  if(kopt>=k.and.kopt/=kmax.and..not.reduct)then
   fact=max(scale/alf(kopt-1,kopt),scalmx)
   if(a(kopt+1)*fact<=wrkmin)then
    hnext=h/fact
    kopt=kopt+1
   endif
  endif

  END SUBROUTINE stepper_bs


  SUBROUTINE pzextr(iest,t_est,xest,xz,dx)
  IMPLICIT NONE
  INTEGER,INTENT(in) :: iest
  REAL(dp),INTENT(in) :: t_est
  REAL(dp),INTENT(in) :: xest(dim)
  REAL(dp),INTENT(out) :: xz(dim),dx(dim)
  INTEGER,PARAMETER :: iest_max=16
  INTEGER j
  REAL(dp) delta,f1,f2
  REAL(dp) d(dim),tmp(dim),q(dim)
  REAL(dp),DIMENSION(iest_max),SAVE :: t
  REAL(dp),ALLOCATABLE,SAVE :: qcol(:,:)

  if(iest>iest_max)call errstop_master('PZEXTR','Probable misuse, too much &
   &extrapolation')
  if(.not.allocated(qcol))allocate(qcol(dim,iest_max))
  t(iest)=t_est
  dx(:)=xest(:)
  xz(:)=xest(:)
  if(iest==1)then
   qcol(:,1)=xest(:)
  else
   d(:)=xest(:)
   do j=1,iest-1
    delta=1.d0/(t(iest-j)-t_est)
    f1=t_est*delta
    f2=t(iest-j)*delta
    q(:)=qcol(:,j)
    qcol(:,j)=dx(:)
    tmp(:)=d(:)-q(:)
    dx(:)=f1*tmp(:)
    d(:)=f2*tmp(:)
    xz(:)=xz(:)+dx(:)
   enddo
   qcol(:,iest)=dx(:)
  endif

  END SUBROUTINE pzextr


  SUBROUTINE mmid(x,vel,ts,htot,nstep,xout,velocity)
  IMPLICIT NONE
  INTEGER,INTENT(in) :: nstep
  REAL(dp),INTENT(in) :: ts,htot
  REAL(dp),INTENT(in) :: x(dim),vel(dim)
  REAL(dp),INTENT(out) :: xout(dim)

  INTERFACE
   SUBROUTINE velocity(t,x,vel) ! was x,y,vel
    USE dsp
    USE store, ONLY : dim
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: t
    REAL(dp),INTENT(in) :: x(dim)
    REAL(dp),INTENT(out) :: vel(dim)
   END SUBROUTINE velocity
  END INTERFACE

  INTEGER n
  REAL(dp) h,h2,t
  REAL(dp) xm(dim),xn(dim)

  h=htot/real(nstep,dp)
  xm=x
  xn=x+h*vel
  t=ts+h
  call velocity(t,xn,xout)
  h2=2.d0*h
  do n=2,nstep
   call swap(xm,xn)
   xn=xn+h2*xout
   t=t+h
   call velocity(t,xn,xout)
  enddo
  xout=0.5d0*(xm+xn+h*xout)

  END SUBROUTINE mmid


  FUNCTION arth(first,increment,n)
  INTEGER,INTENT(in) :: n
  REAL(dp),INTENT(in) :: first,increment
  INTEGER k,k2
  INTEGER,PARAMETER :: npar_arth=16,npar2_arth=8
  REAL(dp) arth(n),temp
  if(n>0)arth(1)=first
  if(n<=npar_arth)then
   do k=2,n
    arth(k)=arth(k-1)+increment
   enddo
  else
   do k=2,npar2_arth
    arth(k)=arth(k-1)+increment
   enddo
   temp=increment*npar2_arth
   k=npar2_arth
   do
    if(k>=n)exit
    k2=k+k
    arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
    temp=temp+temp
    k=k2
   enddo
  endif
  END FUNCTION arth


  FUNCTION arth_i(first,increment,n)
  INTEGER,INTENT(in) :: first,increment,n
  INTEGER arth_i(n)
  INTEGER k,k2,temp
  INTEGER,PARAMETER :: npar_arth=16,npar2_arth=8
  if(n>0)arth_i(1)=first
  if(n<=npar_arth)then
   do k=2,n
    arth_i(k)=arth_i(k-1)+increment
   enddo
  else
   do k=2,npar2_arth
    arth_i(k)=arth_i(k-1)+increment
   enddo
   temp=increment*npar2_arth
   k=npar2_arth
   do
    if(k>=n)exit
    k2=k+k
    arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
    temp=temp+temp
    k=k2
   enddo
  endif
  END FUNCTION arth_i


  RECURSIVE FUNCTION cumsum(arr,seed) result(ans)
  INTEGER,INTENT(in) :: arr(:)
  INTEGER,OPTIONAL,INTENT(in) :: seed
  INTEGER,DIMENSION(size(arr)) :: ans
  INTEGER n,j,sd
  INTEGER,PARAMETER :: npar_cumsum=16
  n=size(arr)
  if(n==0)return
  sd=0
  if(present(seed))sd=seed
  ans(1)=arr(1)+sd
  if(n<npar_cumsum)then
   do j=2,n
    ans(j)=ans(j-1)+arr(j)
   enddo
  else
   ans(2:n:2)=cumsum(arr(2:n:2)+arr(1:n-1:2),sd)
   ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
  endif
  END FUNCTION cumsum


  FUNCTION iminloc(arr)
  REAL(dp),INTENT(in) :: arr(:)
  INTEGER imin(1)
  INTEGER iminloc
  imin=minloc(arr(:))
  iminloc=imin(1)
  END FUNCTION iminloc


  FUNCTION outerdiff(a,b)
  REAL(dp),INTENT(in) :: a(:),b(:)
  REAL(dp),DIMENSION(size(a),size(b)) :: outerdiff
  outerdiff=spread(a,dim=2,ncopies=size(b))-spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff


  FUNCTION outerdiff_i(a,b)
  INTEGER,INTENT(in) :: a(:),b(:)
  INTEGER,DIMENSION(size(a),size(b)) :: outerdiff_i
  outerdiff_i=spread(a,dim=2,ncopies=size(b))-spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_i


  FUNCTION outerprod(a,b)
  REAL(dp),INTENT(in) :: a(:),b(:)
  REAL(dp),DIMENSION(size(a),size(b)) :: outerprod
  outerprod=spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod


  FUNCTION upper_triangle(j,k)
  INTEGER,INTENT(in) :: j,k
  LOGICAL upper_triangle(j,k)
  INTEGER n
  n=0
  upper_triangle=(outerdiff_i(arth_i(1,1,j),arth_i(1,1,k))<n)
  END FUNCTION upper_triangle


  SUBROUTINE swap(a,b)
  REAL(dp),INTENT(inout) :: a(dim),b(dim)
  REAL(dp) dum(dim)
  dum=a
  a=b
  b=dum
  END SUBROUTINE swap


END MODULE bulirsch_stoer

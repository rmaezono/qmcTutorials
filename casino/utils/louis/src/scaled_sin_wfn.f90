MODULE scaled_sin_wfn
 USE dsp
 IMPLICIT NONE
 PRIVATE
 PUBLIC ss_velocity_deBB_1d,ss_velocity_deBB_2d,ss_velocity_deBB_3d,&
  &ss_velocity_curl1_2d,ss_velocity_curl1_3d,ss_eval_psisq_t0_1d,&
  &ss_eval_psisq_t0_2d,ss_eval_psisq_t0_3d,&
  &ss_eval_psisq_1d,ss_eval_psisq_2d,ss_eval_psisq_3d,&
  &ss_energy_and_variance


CONTAINS


 SUBROUTINE ss_velocity_deBB_1d(t,x,vel)
!--------------------------------------------------------------!
! Calculate the de Broglie velocity ('VELOCITY = "deBB"') for  !
! one-dimensional case.                                        !
!                                                              !
! Output: vel(1)                                               !
!--------------------------------------------------------------!
 USE store, ONLY : dim,nmodes_dim,theta_1d,weights_1d
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: t
 REAL(dp),INTENT(in) :: x(dim)
 REAL(dp),INTENT(out) :: vel(dim)
 INTEGER i
 REAL(dp) ii,t1,t2,amp
 REAL(dp) cos1(nmodes_dim),sin1(nmodes_dim)
 COMPLEX(dp) dfx,psi,phase,exparg

 dfx=(0.d0,0.d0) ; psi=(0.d0,0.d0)

 cos1(1)=cos(x(1)) ; sin1(1)=sin(x(1))
 do i=2,nmodes_dim
  t1=real(i,dp)*x(1) ; cos1(i)=cos(t1) ; sin1(i)=sin(t1)
 enddo

 t2=0.5d0*t
 do i=1,nmodes_dim
  ii=real(i*i,dp)
  exparg=cmplx(0.d0,theta_1d(i)-t2*ii,dp)
  phase=exp(exparg)
  amp=weights_1d(i)
  dfx=dfx+(real(i,dp)*amp*cos1(i))*phase
  psi=psi+sin1(i)*amp*phase
 enddo

! vel(1)=aimag(dfx/psi)
! Optimized version of previous line:
 t1=1.d0/(real(psi,dp)**2+aimag(psi)**2)
 vel(1)=(aimag(dfx)*real(psi,dp)-real(dfx,dp)*aimag(psi))*t1

 END SUBROUTINE ss_velocity_deBB_1d


 SUBROUTINE ss_velocity_deBB_2d(t,x,vel)
!----------------------------------------------------------------!
! Calculate the de Broglie velocity ('VELOCITY = "deBB"').       !
! two-dimensional case.                                          !
!                                                                !
! Rate limiting step : ~77% of time in this routine for 2D case. !
!                                                                !
! Output: vel(1:2)                                               !
!----------------------------------------------------------------!
 USE store, ONLY : dim,nmodes_dim,theta_2d,weights_2d
! USE run_control, ONLY : timer
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: t
 REAL(dp),INTENT(in) :: x(dim)
 REAL(dp),INTENT(out) :: vel(dim)
 INTEGER i,j
 REAL(dp) ij2,t1,t2,t3,amp
 REAL(dp),DIMENSION(nmodes_dim) :: cos1,cos2,sin1,sin2
 COMPLEX(dp) dfx,dfy,psi,phase,exparg
 INTRINSIC exp

!call timer('VELOCITY',.true.)

 dfx=(0.d0,0.d0) ; dfy=(0.d0,0.d0) ; psi=(0.d0,0.d0)

 cos1(1)=cos(x(1)) ; sin1(1)=sin(x(1)) ; cos2(1)=cos(x(2)) ; sin2(1)=sin(x(2))
 do i=2,nmodes_dim
  t1=real(i,dp)*x(1) ; cos1(i)=cos(t1) ; sin1(i)=sin(t1)
  t2=real(i,dp)*x(2) ; cos2(i)=cos(t2) ; sin2(i)=sin(t2)
 enddo

 t3=0.5d0*t
 do j=1,nmodes_dim
  do i=1,nmodes_dim
   ij2=real(i*i+j*j,dp)
   exparg=cmplx(0.d0,theta_2d(i,j)-t3*ij2,dp)
   phase=exp(exparg)
   amp=weights_2d(i,j)
   dfy=dfy+(real(j,dp)*sin1(i)*cos2(j))*amp*phase
   dfx=dfx+(real(i,dp)*cos1(i)*sin2(j))*amp*phase
   psi=psi+sin1(i)*sin2(j)*amp*phase
  enddo
 enddo

! vel(1)=aimag(dfx/psi) ; vel(2)=aimag(dfy/psi)
! Optimized version of previous line:
 t1=1.d0/(real(psi,dp)**2+aimag(psi)**2)
 vel(1)=(aimag(dfx)*real(psi,dp)-real(dfx,dp)*aimag(psi))*t1
 vel(2)=(aimag(dfy)*real(psi,dp)-real(dfy,dp)*aimag(psi))*t1

!call timer('VELOCITY',.false.)

 END SUBROUTINE ss_velocity_deBB_2d


 SUBROUTINE ss_velocity_deBB_3d(t,x,vel)
!--------------------------------------------------------------!
! Calculate the de Broglie velocity ('VELOCITY = "deBB"') for  !
! three-dimensional case.                                      !
!                                                              !
! Output: vel(1:3)                                             !
!--------------------------------------------------------------!
 USE store, ONLY : dim,nmodes_dim,theta_3d,weights_3d
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: t
 REAL(dp),INTENT(in) :: x(dim)
 REAL(dp),INTENT(out) :: vel(dim)
 INTEGER i,j,k
 REAL(dp) ijk2,t1,t2,t3,t4,amp
 REAL(dp),DIMENSION(nmodes_dim) :: cos1,sin1,cos2,sin2,cos3,sin3
 COMPLEX(dp) dfx,dfy,dfz,psi,phase,exparg

 dfx=(0.d0,0.d0) ; dfy=(0.d0,0.d0) ; dfz=(0.d0,0.d0) ; psi=(0.d0,0.d0)

 cos1(1)=cos(x(1)) ; sin1(1)=sin(x(1)) 
 cos2(1)=cos(x(2)) ; sin2(1)=sin(x(2))
 cos3(1)=cos(x(3)) ; sin3(1)=sin(x(3))

 do i=2,nmodes_dim
  t1=real(i,dp)*x(1) ; cos1(i)=cos(t1) ; sin1(i)=sin(t1)
  t2=real(i,dp)*x(2) ; cos2(i)=cos(t2) ; sin2(i)=sin(t2)
  t3=real(i,dp)*x(3) ; cos3(i)=cos(t3) ; sin3(i)=sin(t3)
 enddo

 t4=0.5d0*t
 do k=1,nmodes_dim
  do j=1,nmodes_dim
   do i=1,nmodes_dim
    ijk2=real(i*i+j*j+k*k,dp)
    exparg=cmplx(0.d0,theta_3d(i,j,k)-t4*ijk2,dp)
    phase=exp(exparg)
    amp=weights_3d(i,j,k)
    dfx=dfx+(real(i,dp)*cos1(i)*sin2(j)*sin3(k))*amp*phase
    dfy=dfy+(real(j,dp)*sin1(i)*cos2(j)*sin3(k))*amp*phase
    dfz=dfz+(real(k,dp)*sin1(i)*sin2(j)*cos3(k))*amp*phase
    psi=psi+sin1(i)*sin2(j)*sin3(k)*amp*phase
   enddo
  enddo
 enddo

! vel(1)=aimag(dfx/psi)
! vel(2)=aimag(dfy/psi)
! Replace above by faster version
 t1=1.d0/(real(psi,dp)**2+aimag(psi)**2)
 vel(1)=(aimag(dfx)*real(psi,dp)-real(dfx,dp)*aimag(psi))*t1
 vel(2)=(aimag(dfy)*real(psi,dp)-real(dfy,dp)*aimag(psi))*t1
 vel(3)=(aimag(dfz)*real(psi,dp)-real(dfz,dp)*aimag(psi))*t1

 END SUBROUTINE ss_velocity_deBB_3d


 SUBROUTINE ss_velocity_curl1_2d(t,x,vel)
!----------------------------------------------------------------!
! Calculate the de Broglie velocity + curl psi^2 in 2D.          !
! ('VELOCITY = "curl1"').                                        !
!                                                                !
! Weight of curl term 'curlweight' may be defined in input.      !
!                                                                !
! Output: vel(1:2)                                               !
!----------------------------------------------------------------!
 USE store, ONLY : dim,nmodes_dim,theta_2d,weights_2d,curlweight
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: t
 REAL(dp),INTENT(in) :: x(dim)
 REAL(dp),INTENT(out) :: vel(dim)
 INTEGER i,j
 REAL(dp) ij2,t1,t2,t3,amp
 REAL(dp),DIMENSION(nmodes_dim) :: cos1,cos2,sin1,sin2
 COMPLEX(dp) dfx,dfy,psi,phase,exparg
 INTRINSIC exp

 dfx=(0.d0,0.d0) ; dfy=(0.d0,0.d0) ; psi=(0.d0,0.d0)

 cos1(1)=cos(x(1)) ; sin1(1)=sin(x(1)) ; cos2(1)=cos(x(2)) ; sin2(1)=sin(x(2))
 do i=2,nmodes_dim
  t1=real(i,dp)*x(1) ; cos1(i)=cos(t1) ; sin1(i)=sin(t1)
  t2=real(i,dp)*x(2) ; cos2(i)=cos(t2) ; sin2(i)=sin(t2)
 enddo

 t3=0.5d0*t
 do j=1,nmodes_dim
  do i=1,nmodes_dim
   ij2=real(i*i+j*j,dp)
   exparg=cmplx(0.d0,theta_2d(i,j)-t3*ij2,dp)
   phase=exp(exparg)
   amp=weights_2d(i,j)
   dfy=dfy+real(j,dp)*sin1(i)*cos2(j)*amp*phase
   dfx=dfx+real(i,dp)*cos1(i)*sin2(j)*amp*phase
   psi=psi+sin1(i)*sin2(j)*amp*phase
  enddo
 enddo

! vel(1)=aimag(dfx/psi) ; vel(2)=aimag(dfy/psi)
! Optimized version of previous line:
 t1=1.d0/(real(psi,dp)**2+aimag(psi)**2)
 vel(1)=(aimag(dfx)*real(psi,dp)-real(dfx,dp)*aimag(psi))*t1
 vel(2)=(aimag(dfy)*real(psi,dp)-real(dfy,dp)*aimag(psi))*t1

! Extra curl term
! vel(1)=vel(1)+curlweight*real(dfy/psi,dp)
! vel(2)=vel(2)-curlweight*real(dfx/psi,dp)
! Optimized version of previous two lines:
 t1=1.d0/(real(psi,dp)**2+aimag(psi)**2)
 vel(1)=vel(1)+curlweight*t1*(real(psi,dp)*real(dfy,dp)+aimag(psi)*aimag(dfy)) 
 vel(2)=vel(2)-curlweight*t1*(real(psi,dp)*real(dfx,dp)+aimag(psi)*aimag(dfx)) 

 END SUBROUTINE ss_velocity_curl1_2d


 SUBROUTINE ss_velocity_curl1_3d(t,x,vel)
!--------------------------------------------------------------!
! Calculate the de Broglie velocity + curl psi^2 in 3D.        !
! ('VELOCITY = "curl1"').                                      !
!                                                              !
! Weight of curl term 'curlweight' may be defined in input.    !
!                                                              !
! Output: vel(1:3)                                             !
!--------------------------------------------------------------!
 USE store, ONLY : dim,nmodes_dim,theta_3d,weights_3d,curlweight
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: t
 REAL(dp),INTENT(in) :: x(dim)
 REAL(dp),INTENT(out) :: vel(dim)
 INTEGER i,j,k
 REAL(dp) ijk2,t1,t2,t3,t4,amp
 REAL(dp),DIMENSION(nmodes_dim) :: cos1,sin1,cos2,sin2,cos3,sin3
 COMPLEX(dp) dfx,dfy,dfz,psi,phase,exparg

 dfx=(0.d0,0.d0) ; dfy=(0.d0,0.d0) ; dfz=(0.d0,0.d0) ; psi=(0.d0,0.d0)

 cos1(1)=cos(x(1)) ; sin1(1)=sin(x(1))
 cos2(1)=cos(x(2)) ; sin2(1)=sin(x(2))
 cos3(1)=cos(x(3)) ; sin3(1)=sin(x(3))

 do i=2,nmodes_dim
  t1=real(i,dp)*x(1) ; cos1(i)=cos(t1) ; sin1(i)=sin(t1)
  t2=real(i,dp)*x(2) ; cos2(i)=cos(t2) ; sin2(i)=sin(t2)
  t3=real(i,dp)*x(3) ; cos3(i)=cos(t3) ; sin3(i)=sin(t3)
 enddo

 t4=0.5d0*t
 do k=1,nmodes_dim
  do j=1,nmodes_dim
   do i=1,nmodes_dim
    ijk2=real(i*i+j*j+k*k,dp)
    exparg=cmplx(0.d0,theta_3d(i,j,k)-t4*ijk2,dp)
    phase=exp(exparg)
    amp=weights_3d(i,j,k)
    dfx=dfx+real(i,dp)*cos1(i)*sin2(j)*sin3(k)*amp*phase
    dfy=dfy+real(j,dp)*sin1(i)*cos2(j)*sin3(k)*amp*phase
    dfz=dfz+real(k,dp)*sin1(i)*sin2(j)*cos3(k)*amp*phase
    psi=psi+sin1(i)*sin2(j)*sin3(k)*amp*phase
   enddo
  enddo
 enddo

! vel(1)=aimag(dfx/psi)
! vel(2)=aimag(dfy/psi)
! Replace previous to lines by optimized version
 t1=1.d0/(real(psi,dp)**2+aimag(psi)**2)
 vel(1)=(aimag(dfx)*real(psi,dp)-real(dfx,dp)*aimag(psi))*t1
 vel(2)=(aimag(dfy)*real(psi,dp)-real(dfy,dp)*aimag(psi))*t1
 vel(3)=(aimag(dfz)*real(psi,dp)-real(dfz,dp)*aimag(psi))*t1

! Extra curl term
! vel(1)=vel(1)+curlweight*real((dfy-dfz)/psi,dp)
! vel(2)=vel(2)-curlweight*real((dfz-dfx)/psi,dp)
! vel(3)=vel(3)-curlweight*real((dfx-dfy)/psi,dp)
! Replace previous three lines by optimized version
 t1=1.d0/(real(psi,dp)**2+aimag(psi)**2)
 vel(1)=vel(1)+curlweight*t1*(real(psi,dp)*(real(dfy,dp)-real(dfz,dp))+&
  &aimag(psi)*(aimag(dfy)-aimag(dfz)))
 vel(2)=vel(2)+curlweight*t1*(real(psi,dp)*(real(dfz,dp)-real(dfx,dp))+&
  &aimag(psi)*(aimag(dfz)-aimag(dfx)))
 vel(3)=vel(3)+curlweight*t1*(real(psi,dp)*(real(dfx,dp)-real(dfy,dp))+&
  &aimag(psi)*(aimag(dfx)-aimag(dfy)))

 END SUBROUTINE ss_velocity_curl1_3d

 SUBROUTINE ss_eval_psisq_t0_1d(x,rho)
!----------------------------------------------------------------------------!
! Evaluate psi^2 in 1d at time zero and position x.                          !
!----------------------------------------------------------------------------!
 USE store,ONLY : nmodes_dim,theta_1d,weights_1d
 REAL(dp),INTENT(in) :: x(1)
 REAL(dp),INTENT(out) :: rho
 INTEGER i
 REAL(dp) amp
 COMPLEX(dp) psi,zz

 psi=(0.d0,0.d0)
 do i=1,nmodes_dim
  zz=cmplx(0.d0,theta_1d(i),dp)
  amp=weights_1d(i)
  psi=psi+amp*sin(real(i,dp)*x(1))*exp(zz)
 enddo
 rho=real(psi*conjg(psi),dp)

 END SUBROUTINE ss_eval_psisq_t0_1d


 SUBROUTINE ss_eval_psisq_t0_2d(x,rho)
!----------------------------------------------------------------------------!
! Evaluate psi^2 in 2d at time zero and position x.                          !
!----------------------------------------------------------------------------!
 USE store,ONLY : nmodes_dim,theta_2d,weights_2d
 REAL(dp),INTENT(in) :: x(2)
 REAL(dp),INTENT(out) :: rho
 INTEGER i,j
 REAL(dp) amp
 COMPLEX(dp) psi,zz

 psi=(0.d0,0.d0)
 do j=1,nmodes_dim
  do i=1,nmodes_dim
   zz=cmplx(0.d0,theta_2d(i,j),dp)
   amp=weights_2d(i,j)
   psi=psi+amp*(sin(real(i,dp)*x(1))*sin(real(j,dp)*x(2)))*exp(zz)
  enddo
 enddo
 rho=real(psi*conjg(psi),dp)

 END SUBROUTINE ss_eval_psisq_t0_2d


 SUBROUTINE ss_eval_psisq_t0_3d(x,rho)
!----------------------------------------------------------------------------!
! Evaluate psi^2 in 3d at time zero and position x.                          !
!----------------------------------------------------------------------------!
 USE store,ONLY : nmodes_dim,theta_3d,weights_3d
 REAL(dp),INTENT(in) :: x(3)
 REAL(dp),INTENT(out) :: rho
 INTEGER i,j,k
 REAL(dp) amp
 COMPLEX(dp) psi,zz

 psi=(0.d0,0.d0)
 do k=1,nmodes_dim
  do j=1,nmodes_dim
   do i=1,nmodes_dim
    zz=cmplx(0.d0,theta_3d(i,j,k),dp)
    amp=weights_3d(i,j,k)
    psi=psi+(sin(real(i,dp)*x(1))*sin(real(j,dp)*x(2))*sin(real(k,dp)*x(3)))*&
     &amp*exp(zz)
   enddo
  enddo
 enddo
 rho=real(psi*conjg(psi),dp)

 END SUBROUTINE ss_eval_psisq_t0_3d


 SUBROUTINE ss_eval_psisq_1d(t,x,rho)
!----------------------------------------------------------------------------!
! Evaluate psi^2 in 1d at time t and position x.                             !
!----------------------------------------------------------------------------!
 USE store,ONLY : nmodes_dim,theta_1d,weights_1d
 REAL(dp),INTENT(in) :: t
 REAL(dp),INTENT(in) :: x(1)
 REAL(dp),INTENT(out) :: rho
 INTEGER i
 REAL(dp) t1,amp
 COMPLEX(dp) psi,zz

 psi=(0.d0,0.d0)
 t1=0.5d0*t
 do i=1,nmodes_dim
  zz=cmplx(0.d0,theta_1d(i)-t1*(real(i**2,dp)),dp)
  amp=weights_1d(i)
  psi=psi+sin(real(i,dp)*x(1))*amp*exp(zz)
 enddo
 rho=real(psi*conjg(psi),dp)

 END SUBROUTINE ss_eval_psisq_1d


 SUBROUTINE ss_eval_psisq_2d(t,x,rho)
!----------------------------------------------------------------------------!
! Evaluate psi^2 in 2d at time t and position x.                             !
!----------------------------------------------------------------------------!
 USE store, ONLY : nmodes_dim,theta_2d,weights_2d
 REAL(dp),INTENT(in) :: t
 REAL(dp),INTENT(in) :: x(2)
 REAL(dp),INTENT(out) :: rho
 INTEGER i,j
 REAL(dp) t1,amp
 COMPLEX(dp) psi,zz

 psi=(0.d0,0.d0)
 t1=0.5d0*t
 do j=1,nmodes_dim
  do i=1,nmodes_dim
   zz=cmplx(0.d0,theta_2d(i,j)-t1*(real(i**2+j**2,dp)),dp)
   amp=weights_2d(i,j)
   psi=psi+(sin(real(i,dp)*x(1))*sin(real(j,dp)*x(2)))*amp*exp(zz)
  enddo
 enddo
 rho=real(psi*conjg(psi),dp)

 END SUBROUTINE ss_eval_psisq_2d


 SUBROUTINE ss_eval_psisq_3d(t,x,rho)
!----------------------------------------------------------------------------!
! Evaluate psi^2 in 3d at time t and position x.                             !
!----------------------------------------------------------------------------!
 USE store, ONLY : nmodes_dim,theta_3d,weights_3d
 REAL(dp),INTENT(inout) :: t
 REAL(dp),INTENT(in) :: x(3)
 REAL(dp),INTENT(out) :: rho
 INTEGER i,j,k
 REAL(dp) t1,amp
 COMPLEX(dp) psi,zz

 psi=(0.d0,0.d0)
 t1=0.5d0*t
 do k=1,nmodes_dim
  do j=1,nmodes_dim
   do i=1,nmodes_dim
    zz=cmplx(0.d0,theta_3d(i,j,k)-t1*(real(i**2+j**2+k**2,dp)),dp)
    amp=weights_3d(i,k,k)
    psi=psi+(sin(real(i,dp)*x(1))*sin(real(j,dp)*x(2))*sin(real(k,dp)*x(3)))*&
     &amp*exp(zz)
   enddo
  enddo
 enddo
 rho=real(psi*conjg(psi),dp)

 END SUBROUTINE ss_eval_psisq_3d


 SUBROUTINE ss_energy_and_variance
 USE run_control, ONLY : errstop
 USE store, ONLY : dim,nmodes_dim,mean_energy,variance,cvariance,pi,cell_x
 IMPLICIT NONE
 INTEGER i,j,k
 REAL(dp) t1,esum,esqsum,mean_sq_energy,scale

 esum=0.d0 ; esqsum=0.d0

 scale=pi**2/cell_x**2

 select case(dim)
  case(1)
   do i=1,nmodes_dim
    t1=real(i**2,dp)
    esum=esum+t1
    esqsum=esqsum+t1**2
   enddo
   mean_energy=0.5d0*esum*scale/real(nmodes_dim,dp)
   mean_sq_energy=0.25d0*esqsum*scale**2/real(nmodes_dim,dp)
  case(2)
   do j=1,nmodes_dim
    do i=1,nmodes_dim
     t1=real(i**2+j**2,dp)
     esum=esum+t1
     esqsum=esqsum+t1**2
    enddo
   enddo
   mean_energy=0.5d0*esum*scale/real(nmodes_dim**2,dp)
   mean_sq_energy=0.25d0*esqsum*scale**2/real(nmodes_dim**2,dp)
  case(3)
   do k=1,nmodes_dim
    do j=1,nmodes_dim
     do i=1,nmodes_dim
      t1=real(i**2+j**2+k**2,dp)
      esum=esum+t1
      esqsum=esqsum+t1**2
     enddo
    enddo
   enddo
   mean_energy=0.5d0*esum*scale/real(nmodes_dim**3,dp)
   mean_sq_energy=0.25d0*esqsum*scale**2/real(nmodes_dim**3,dp)
  case default
   call errstop('ENERGY_AND_VARIANCE','Invalid dimension.')
 end select

 variance=mean_sq_energy-mean_energy**2

 esum=0.d0
 esqsum=0.d0
 
 select case(dim)
  case(1)
   do i=1,nmodes_dim
    t1=0.5d0*real(i**2,dp)*scale
    esqsum=esqsum+(t1-mean_energy)**2
    esum=esum+(t1-mean_energy)
   enddo
   cvariance=(esqsum-esum**2/nmodes_dim)/nmodes_dim
  case(2)
   do j=1,nmodes_dim
    do i=1,nmodes_dim
     t1=0.5d0*real(i**2+j**2,dp)*scale
     esqsum=esqsum+(t1-mean_energy)**2
     esum=esum+(t1-mean_energy)
    enddo
   enddo
   cvariance=(esqsum-esum**2/nmodes_dim**2)/nmodes_dim**2
  case(3)
   do k=1,nmodes_dim
    do j=1,nmodes_dim
     do i=1,nmodes_dim
      t1=0.5d0*real(i**2+j**2+k**2,dp)*scale
      esqsum=esqsum+(t1-mean_energy)**2
      esum=esum+(t1-mean_energy)
     enddo
    enddo
   enddo
   cvariance=(esqsum-esum**2/nmodes_dim**3)/nmodes_dim**3
  case default
   call errstop('ENERGY_AND_VARIANCE','Invalid dimension.')
 end select
 
 END SUBROUTINE ss_energy_and_variance


END MODULE scaled_sin_wfn

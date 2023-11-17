module slaarnabc_noopt
use dsp
use format_utils, only : wout,i2s
use slaarnabt,    only : dscal,daxpy
use parallel,     only : am_master
use run_control,  only : errwarn
implicit none
private
public gaussian_elimination,redo_gaussian_elimination,backwards_substi&
&tution,param_dep_list,construct_forced_zero_mask
real(dp),parameter :: xyzzyaaaa1=1.d-10
real(dp),parameter :: xyzzyaaab1=1.d-100
logical,parameter :: xyzzyaaac1=.false.
contains
subroutine gaussian_elimination(nparam,neqn,matrix,rhs,neqn_actual,sol&
&ves_for,eqn_pivot,ierr)
implicit none
integer,intent(in) :: nparam,neqn
integer,intent(out) :: solves_for(neqn),eqn_pivot(neqn),neqn_actual,ie&
&rr
real(dp),intent(inout) :: matrix(nparam,neqn),rhs(neqn)
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2
real(dp) xyzzyaaae2,xyzzyaaaf2
logical xyzzyaaag2(neqn)
do xyzzyaaaa2=1,neqn
xyzzyaaag2(xyzzyaaaa2)=.false.
xyzzyaaae2=maxval(abs(matrix(1:nparam,xyzzyaaaa2)))
if(xyzzyaaae2<xyzzyaaab1)matrix(1:nparam,xyzzyaaaa2)=0.d0
enddo
xyzzyaaaa2=0
do xyzzyaaad2=1,nparam
xyzzyaaab2=0
xyzzyaaaf2=0.d0
do xyzzyaaac2=1,neqn
if(xyzzyaaag2(xyzzyaaac2))cycle
xyzzyaaae2=matrix(xyzzyaaad2,xyzzyaaac2)
if(abs(xyzzyaaae2)>abs(xyzzyaaaf2))then
xyzzyaaab2=xyzzyaaac2
xyzzyaaaf2=xyzzyaaae2
endif
enddo
if(xyzzyaaab2==0)cycle
if(abs(xyzzyaaaf2)<xyzzyaaab1)then
where(.not.xyzzyaaag2(1:neqn))matrix(xyzzyaaad2,1:neqn)=0.d0
cycle
endif
xyzzyaaaa2=xyzzyaaaa2+1
eqn_pivot(xyzzyaaaa2)=xyzzyaaab2
solves_for(xyzzyaaaa2)=xyzzyaaad2
xyzzyaaag2(xyzzyaaab2)=.true.
xyzzyaaae2=1.d0/xyzzyaaaf2
rhs(xyzzyaaab2)=rhs(xyzzyaaab2)*xyzzyaaae2
matrix(xyzzyaaad2,xyzzyaaab2)=1.d0
if(xyzzyaaad2<nparam)then
call dscal(nparam-xyzzyaaad2,xyzzyaaae2,matrix(xyzzyaaad2+1,xyzzyaaab2&
&),1)
where(abs(matrix(xyzzyaaad2+1:nparam,xyzzyaaab2))<xyzzyaaaa1)matrix(xy&
&zzyaaad2+1:nparam,xyzzyaaab2)=0.d0
endif
do xyzzyaaac2=1,neqn
if(xyzzyaaag2(xyzzyaaac2))cycle
xyzzyaaae2=-matrix(xyzzyaaad2,xyzzyaaac2)
matrix(xyzzyaaad2,xyzzyaaac2)=0.d0
if(abs(xyzzyaaae2)<xyzzyaaab1)cycle
rhs(xyzzyaaac2)=rhs(xyzzyaaac2)+xyzzyaaae2*rhs(xyzzyaaab2)
if(xyzzyaaad2<nparam)then
call daxpy(nparam-xyzzyaaad2,xyzzyaaae2,matrix(xyzzyaaad2+1,xyzzyaaab2&
&),1,matrix(xyzzyaaad2+1,xyzzyaaac2),1)
xyzzyaaae2=maxval(abs(matrix(xyzzyaaad2+1:nparam,xyzzyaaac2)))
if(xyzzyaaae2>xyzzyaaaa1)then
xyzzyaaae2=xyzzyaaaa1*xyzzyaaae2
where(abs(matrix(xyzzyaaad2+1:nparam,xyzzyaaac2))<xyzzyaaae2)matrix(xy&
&zzyaaad2+1:nparam,xyzzyaaac2)=0.d0
else
matrix(xyzzyaaad2+1:nparam,xyzzyaaac2)=0.d0
endif
endif
enddo
if(xyzzyaaaa2==neqn)exit
enddo
neqn_actual=xyzzyaaaa2
ierr=0
do xyzzyaaaa2=1,neqn
if(.not.xyzzyaaag2(xyzzyaaaa2).and.abs(rhs(xyzzyaaaa2))>xyzzyaaaa1)ier&
&r=ierr+1
enddo
end subroutine gaussian_elimination
subroutine redo_gaussian_elimination(nparam,neqn,matrix,rhs,neqn_actua&
&l,solves_for,eqn_pivot,bad_params)
implicit none
integer,intent(in) :: nparam,neqn,neqn_actual,solves_for(neqn_actual),&
&eqn_pivot(neqn_actual)
real(dp),intent(inout) :: matrix(nparam,neqn),rhs(neqn)
logical,intent(out),optional :: bad_params
integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3,xyzzyaaad3,xyzzyaaae3,xyzzyaa&
&af3(neqn_actual),xyzzyaaag3(neqn_actual),xyzzyaaah3,xyzzyaaai3
real(dp) xyzzyaaaj3
xyzzyaaai3=0
do xyzzyaaaa3=1,neqn_actual
xyzzyaaae3=solves_for(xyzzyaaaa3)
xyzzyaaab3=eqn_pivot(xyzzyaaaa3)
xyzzyaaaj3=matrix(xyzzyaaae3,xyzzyaaab3)
if(abs(xyzzyaaaj3)>xyzzyaaab1)then
xyzzyaaaj3=1.d0/xyzzyaaaj3
else
xyzzyaaaj3=0.d0
xyzzyaaai3=xyzzyaaai3+1
xyzzyaaaf3(xyzzyaaai3)=xyzzyaaaa3
xyzzyaaag3(xyzzyaaai3)=xyzzyaaae3
endif
rhs(xyzzyaaab3)=rhs(xyzzyaaab3)*xyzzyaaaj3
matrix(xyzzyaaae3,xyzzyaaab3)=1.d0
if(xyzzyaaae3<nparam)then
call dscal(nparam-xyzzyaaae3,xyzzyaaaj3,matrix(xyzzyaaae3+1,xyzzyaaab3&
&),1)
where(abs(matrix(xyzzyaaae3+1:nparam,xyzzyaaab3))<xyzzyaaaa1)matrix(xy&
&zzyaaae3+1:nparam,xyzzyaaab3)=0.d0
endif
do xyzzyaaac3=xyzzyaaaa3+1,neqn_actual
xyzzyaaad3=eqn_pivot(xyzzyaaac3)
xyzzyaaaj3=-matrix(xyzzyaaae3,xyzzyaaad3)
matrix(xyzzyaaae3,xyzzyaaad3)=0.d0
if(abs(xyzzyaaaj3)<xyzzyaaab1)cycle
rhs(xyzzyaaad3)=rhs(xyzzyaaad3)+xyzzyaaaj3*rhs(xyzzyaaab3)
if(xyzzyaaae3<nparam)call daxpy(nparam-xyzzyaaae3,xyzzyaaaj3,matrix(xy&
&zzyaaae3+1,xyzzyaaab3),1,matrix(xyzzyaaae3+1,xyzzyaaad3),1)
enddo
enddo
if(present(bad_params))then
bad_params=xyzzyaaai3>0
elseif(xyzzyaaai3>0)then
call errwarn('REDO_GAUSSIAN_ELIMINATION','non-linear parameters produc&
&e '//trim(i2s(xyzzyaaai3))//' zero(es) in the coefficient matrix for &
&the constraint equations.  We recommend that you stop this calculatio&
&n, change the values of any basis parameters involved to non-integer &
&values, all slightly different from each other, and re-run to avoid p&
&ossible problems.  This run will continue in the hope that the proble&
&m will go away during optimization.  The affected equations and param&
&eters are listed below.')
if(am_master)then
do xyzzyaaah3=1,xyzzyaaai3
call wout('Zero #'//trim(i2s(xyzzyaaah3))//' encountered at equation #&
&'//trim(i2s(xyzzyaaaf3(xyzzyaaah3)))//', parameter #'//trim(i2s(xyzzy&
&aaag3(xyzzyaaah3)))//'.')
enddo
call wout()
endif
endif
end subroutine redo_gaussian_elimination
subroutine backwards_substitution(nparam,neqn,neqn_actual,matrix,rhs,s&
&olves_for,eqn_pivot,param_iszero,param)
implicit none
integer,intent(in) :: nparam,neqn,neqn_actual,solves_for(neqn_actual),&
&eqn_pivot(neqn_actual)
real(dp),intent(in) :: matrix(nparam,neqn),rhs(neqn)
real(dp),intent(inout) :: param(nparam)
logical,intent(in) :: param_iszero(nparam)
integer xyzzyaaaa4,xyzzyaaab4,xyzzyaaac4,xyzzyaaad4
real(dp) xyzzyaaae4,xyzzyaaaf4,xyzzyaaag4,xyzzyaaah4
do xyzzyaaaa4=neqn_actual,1,-1
xyzzyaaac4=solves_for(xyzzyaaaa4)
if(param_iszero(xyzzyaaac4))cycle
xyzzyaaab4=eqn_pivot(xyzzyaaaa4)
xyzzyaaae4=rhs(xyzzyaaab4)
xyzzyaaaf4=0.d0
do xyzzyaaad4=xyzzyaaac4+1,nparam
if(param_iszero(xyzzyaaad4))cycle
xyzzyaaag4=param(xyzzyaaad4)
xyzzyaaaf4=max(abs(xyzzyaaag4),xyzzyaaaf4)
xyzzyaaae4=xyzzyaaae4-matrix(xyzzyaaad4,xyzzyaaab4)*xyzzyaaag4
enddo
xyzzyaaah4=abs(xyzzyaaae4)
if(xyzzyaaah4>0.d0.and.xyzzyaaaf4>0.d0)then
if(xyzzyaaah4<xyzzyaaaf4*xyzzyaaaa1)xyzzyaaae4=0.d0
endif
param(xyzzyaaac4)=xyzzyaaae4
enddo
end subroutine backwards_substitution
subroutine param_dep_list(iparam,nparam,neqn,neqn_actual,matrix,solves&
&_for,eqn_pivot,param_iszero,np,plist,pcoeff)
implicit none
integer,intent(in) :: iparam,nparam,neqn,neqn_actual,solves_for(neqn_a&
&ctual),eqn_pivot(neqn_actual)
integer,intent(inout) :: np,plist(nparam)
real(dp),intent(in) :: matrix(nparam,neqn)
real(dp),intent(inout) :: pcoeff(nparam)
logical,intent(in) :: param_iszero(nparam)
integer xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5,xyzzyaaad5
real(dp) xyzzyaaae5,xyzzyaaaf5,xyzzyaaag5,xyzzyaaah5
np=1
plist(1)=iparam
pcoeff(1)=1.d0
do xyzzyaaaa5=neqn_actual,1,-1
xyzzyaaac5=solves_for(xyzzyaaaa5)
if(param_iszero(xyzzyaaac5))cycle
xyzzyaaab5=eqn_pivot(xyzzyaaaa5)
xyzzyaaaf5=matrix(xyzzyaaac5,xyzzyaaab5)
if(xyzzyaaaf5==0.d0)cycle
xyzzyaaah5=0.d0
do xyzzyaaad5=1,np
if(param_iszero(plist(xyzzyaaad5)))cycle
xyzzyaaae5=matrix(plist(xyzzyaaad5),xyzzyaaab5)
if(xyzzyaaae5==0.d0)cycle
xyzzyaaag5=-xyzzyaaae5/xyzzyaaaf5
xyzzyaaah5=xyzzyaaah5+xyzzyaaag5*pcoeff(xyzzyaaad5)
enddo
if(abs(xyzzyaaah5)>=xyzzyaaaa1)then
np=np+1
plist(np)=xyzzyaaac5
pcoeff(np)=xyzzyaaah5
endif
enddo
end subroutine param_dep_list
subroutine construct_forced_zero_mask(nparam,neqn,neqn_actual,matrix,r&
&hs,solves_for,eqn_pivot,param_iszero)
implicit none
integer,intent(in) :: nparam,neqn,neqn_actual,solves_for(neqn_actual),&
&eqn_pivot(neqn_actual)
real(dp),intent(in) :: matrix(nparam,neqn),rhs(neqn)
logical,intent(inout) :: param_iszero(nparam)
integer xyzzyaaaa6,xyzzyaaab6,xyzzyaaac6
param_iszero=.false.
if(xyzzyaaac1)return
do xyzzyaaaa6=neqn_actual,1,-1
xyzzyaaac6=solves_for(xyzzyaaaa6)
xyzzyaaab6=eqn_pivot(xyzzyaaaa6)
param_iszero(xyzzyaaac6)=rhs(xyzzyaaab6)==0.d0.and.maxval(abs(matrix(x&
&yzzyaaac6+1:nparam,xyzzyaaab6)),.not.param_iszero(xyzzyaaac6+1:nparam&
&))<xyzzyaaaa1
enddo
end subroutine construct_forced_zero_mask
end module slaarnabc_noopt

module slaarnace
use slaarnaag,only : minus_1_over_2c_squared,one_over_4c_squared,amu,n&
&uclear_mass,pi
use dsp,only : dp
use slaarnabg,only : atno,nitot,iontype,isperiodic,dimensionality,peri&
&odicity,volume,nbasis,npcells
use slaarnabk,only : ewald_3d
use slaarnabt,only : ddot,dsum
use slaarnaca,only : zion
use run_control,only : errstop
use store,only : netot,real1_complex2,complex_wf,three_netot,heg_zlaye&
&r,heg_ylayer
implicit none
real(dp) xyzzyaaaa1,darwin_constant,xyzzyaaab1
integer xyzzyaaac1
logical relativistic,eval_maspol,eval_masvel,eval_darwinen,eval_darwin&
&ee,eval_retard,xyzzyaaad1,xyzzyaaae1
private
public relativistic,relativistic_terms,relativity_setup,eval_maspol,ev&
&al_masvel,eval_darwinen,eval_darwinee,eval_retard,darwin_constant
contains
subroutine relativity_setup(isotope_mass,dmc_calc)
use format_utils,only : wout,r2s,wordwrap
use parallel,only : am_master
use run_control,only : errwarn,errwarn_silent
use store,only : nspin,nele,pcharge,pmass,self_term
use slaarnacs,only : psi_s
implicit none
real(dp),intent(in) :: isotope_mass
logical,intent(in) :: dmc_calc
integer xyzzyaaaa2,xyzzyaaab2
real(dp) xyzzyaaac2
character(30) char30
if(am_master)then
call wout('Relativistic corrections setup')
call wout('==============================')
if(complex_wf.and.psi_s/="geminal")call errstop('RELATIVITY_SETUP','Ca&
&nnot use relativistic corrections with complex wave functions at pres&
&ent.')
if(any(pcharge/=-1).or.any(pmass/=1))call errstop('RELATIVITY_SETUP','&
&Cannot use relativistic corrections for non-electron-systems at prese&
&nt.')
if(nspin<2)call errwarn_silent('RELATIVITY_SETUP','Full relativistic c&
&orrections are only available for closed-shell systems. Hope you know&
& what you''re doing.')
if(nspin==2.and.nele(1)/=nele(2))call errwarn_silent('RELATIVITY_SETUP&
&','Full relativistic corrections are only available for closed-shell &
&systems. Hope you know what you''re doing.')
if(nspin>2)call errwarn_silent('RELATIVITY_SETUP','The relativistic co&
&rrection implemented in CASINO are only correct for closed-shell syst&
&ems. Hope you know what you''re doing.')
endif
if(nitot>=1)then
eval_maspol=.true.
if(isotope_mass==0.d0)then
xyzzyaaac2=0.d0
do xyzzyaaaa2=1,nbasis
xyzzyaaac2=xyzzyaaac2+nuclear_mass(atno(xyzzyaaaa2))
enddo
xyzzyaaac2=xyzzyaaac2*dble(npcells)*amu
else
xyzzyaaac2=isotope_mass
endif
xyzzyaaaa1=1.d0/xyzzyaaac2
else
eval_maspol=.false.
endif
eval_darwinee=(dimensionality==3.and.(periodicity==0.or.periodicity==3&
&))
eval_darwinen=(dimensionality==3.and.(periodicity==0.or.periodicity==3&
&)).and.(nitot>0)
xyzzyaaad1=(eval_darwinen.or.eval_darwinee)
if(isperiodic)xyzzyaaab1=dble(netot)*self_term
darwin_constant=0.d0
if(dimensionality==3)then
if(periodicity==3)then
do xyzzyaaaa2=1,nitot
darwin_constant=darwin_constant+zion(iontype(xyzzyaaaa2))
enddo
do xyzzyaaab2=1,nspin
darwin_constant=darwin_constant+pcharge(xyzzyaaab2)*dble(nele(xyzzyaaa&
&b2))
enddo
darwin_constant=darwin_constant/volume*pi*minus_1_over_2c_squared
elseif(isperiodic.and.am_master)then
call errwarn('RELATIVITY_SETUP','Darwin terms not implemented in 1D- a&
&nd 2D-periodic systems.')
endif
elseif(isperiodic.and.am_master)then
call errwarn('RELATIVITY_SETUP','Darwin constant not implemented in 1D&
& and 2D.')
endif
eval_masvel=.true.
eval_retard=.true.
xyzzyaaae1=(all(heg_zlayer==0.d0).and.all(heg_ylayer==0.d0))
if(am_master)then
if(eval_maspol)then
call wout('Mass-polarization term will be evaluated.')
if(nitot>1)call wout(' [This does NOT include all the effects of a fin&
&ite nuclear mass to O(1/M).]')
char30=r2s(xyzzyaaac2,'(es24.16)')
call wout(' Total nuclear mass: '//trim(char30))
else
call wout('Mass-polarization term will not be evaluated.')
endif
if(eval_masvel)then
call wout('Mass-velocity term will be evaluated.')
else
call wout('Mass-velocity term will not be evaluated.')
endif
if(eval_darwinen)then
call wout('Electron-nucleus Darwin term will be evaluated.')
else
call wout('Electron-nucleus Darwin term will not be evaluated.')
endif
if(eval_darwinee)then
call wout('Electron-electron Darwin term will be evaluated.')
else
call wout('Electron-electron Darwin term will not be evaluated.')
endif
if(isperiodic)call wordwrap('It is probably more accurate to calculate&
& the Darwin terms by evaluating the contact pair-correlation function&
& and the charge density at the ions.')
if(darwin_constant/=0.d0)then
char30=r2s(darwin_constant,'(es24.16)')
call wout('Darwin constant (au per particle): '//trim(char30))
endif
if(eval_retard)then
call wout('Retardation term will be evaluated.')
else
call wout('Retardation term will not be evaluated.')
endif
if(dmc_calc)call errwarn('RELATIVITY_SETUP','The estimators for the re&
&lativistic corrections are only approximate in DMC.  Extrapolated est&
&imation will not reduce the order of the error.  This could be fixed &
&if anyone is willing to code up fourth derivatives of the wave functi&
&on.')
call wout()
endif
xyzzyaaac1=max(nitot,1)
end subroutine relativity_setup
subroutine relativistic_terms(dvel,eivecs,relkei,eevecsall,maspol,masv&
&el,darwin_en,darwin_ee,retard,reltot)
implicit none
real(dp),intent(in) :: dvel(3,netot,real1_complex2),eivecs(4,xyzzyaaac&
&1,netot),relkei(netot),eevecsall(4,netot,netot)
real(dp),intent(out) :: maspol,masvel,darwin_en,darwin_ee,retard,relto&
&t
reltot=0.d0
if(eval_maspol)then
call xyzzyaaaf1(dvel,maspol)
reltot=reltot+maspol
endif
if(eval_masvel)then
call xyzzyaaag1(relkei,masvel)
reltot=reltot+masvel
endif
if(xyzzyaaad1)then
call xyzzyaaah1(dvel,eivecs,relkei,eevecsall,darwin_en,darwin_ee)
if(eval_darwinen)reltot=reltot+darwin_en
if(eval_darwinee)reltot=reltot+darwin_ee
endif
if(eval_retard)then
call xyzzyaaai1(dvel,eevecsall,retard)
reltot=reltot+retard
endif
end subroutine relativistic_terms
subroutine xyzzyaaaf1(dvel,maspol)
implicit none
real(dp),intent(in) :: dvel(3,netot,real1_complex2)
real(dp),intent(out) :: maspol
integer xyzzyaaaa4,xyzzyaaab4
maspol=0.d0
do xyzzyaaaa4=1,netot-1
do xyzzyaaab4=xyzzyaaaa4+1,netot
maspol=maspol+dot_product(dvel(1:3,xyzzyaaaa4,1),dvel(1:3,xyzzyaaab4,1&
&))
enddo
enddo
maspol=xyzzyaaaa1*maspol
end subroutine xyzzyaaaf1
subroutine xyzzyaaag1(relkei,masvel)
implicit none
real(dp),intent(in) :: relkei(netot)
real(dp),intent(out) :: masvel
masvel=minus_1_over_2c_squared*ddot(netot,relkei(1),1,relkei(1),1)
end subroutine xyzzyaaag1
subroutine xyzzyaaah1(dvel,eivecs,relkei,eevecsall,darwin_en,darwin_ee&
&)
implicit none
real(dp),intent(in) :: dvel(3,netot,real1_complex2),eivecs(4,xyzzyaaac&
&1,netot),relkei(netot),eevecsall(4,netot,netot)
real(dp),intent(out) :: darwin_en,darwin_ee
integer xyzzyaaaa6,xyzzyaaab6,xyzzyaaac6
real(dp) xyzzyaaad6,xyzzyaaae6,xyzzyaaaf6,xyzzyaaag6,xyzzyaaah6
xyzzyaaad6=ddot(three_netot,dvel(1,1,1),1,dvel(1,1,1),1)-2.d0*dsum(net&
&ot,relkei,1)
if(eval_darwinen)then
xyzzyaaae6=0.d0
if(isperiodic)then
do xyzzyaaaa6=1,netot
do xyzzyaaac6=1,nitot
call ewald_3d(1,eivecs(1,xyzzyaaac6,xyzzyaaaa6),xyzzyaaah6,.false.)
xyzzyaaae6=xyzzyaaae6-zion(iontype(xyzzyaaac6))*xyzzyaaah6
enddo
enddo
else
do xyzzyaaaa6=1,netot
do xyzzyaaac6=1,nitot
xyzzyaaag6=eivecs(4,xyzzyaaac6,xyzzyaaaa6)
if(xyzzyaaag6>0.d0)xyzzyaaae6=xyzzyaaae6-zion(iontype(xyzzyaaac6))/xyz&
&zyaaag6
enddo
enddo
endif
darwin_en=one_over_4c_squared*xyzzyaaad6*xyzzyaaae6
endif
if(eval_darwinee)then
xyzzyaaaf6=0.d0
if(isperiodic)then
do xyzzyaaaa6=1,netot-1
call ewald_3d(netot-xyzzyaaaa6,eevecsall(1,xyzzyaaaa6+1,xyzzyaaaa6),xy&
&zzyaaah6,.false.)
xyzzyaaaf6=xyzzyaaaf6+xyzzyaaah6
enddo
xyzzyaaaf6=xyzzyaaaf6+xyzzyaaab1
else
do xyzzyaaaa6=1,netot-1
do xyzzyaaab6=xyzzyaaaa6+1,netot
xyzzyaaag6=eevecsall(4,xyzzyaaab6,xyzzyaaaa6)
if(xyzzyaaag6>0.d0)xyzzyaaaf6=xyzzyaaaf6+1.d0/xyzzyaaag6
enddo
enddo
endif
darwin_ee=-one_over_4c_squared*xyzzyaaad6*xyzzyaaaf6
endif
end subroutine xyzzyaaah1
subroutine xyzzyaaai1(dvel,eevecsall,retard)
implicit none
real(dp),intent(in) :: dvel(3,netot,real1_complex2),eevecsall(4,netot,&
&netot)
real(dp),intent(out) :: retard
integer xyzzyaaaa7,xyzzyaaab7
real(dp) xyzzyaaac7(3),xyzzyaaad7,xyzzyaaae7,xyzzyaaaf7,xyzzyaaag7
retard=0.d0
if(xyzzyaaae1)then
do xyzzyaaaa7=1,netot-1
do xyzzyaaab7=xyzzyaaaa7+1,netot
xyzzyaaac7=eevecsall(1:3,xyzzyaaab7,xyzzyaaaa7)
xyzzyaaag7=1.d0/eevecsall(4,xyzzyaaab7,xyzzyaaaa7)
xyzzyaaad7=dot_product(xyzzyaaac7,dvel(1:3,xyzzyaaaa7,1))
xyzzyaaae7=dot_product(xyzzyaaac7,dvel(1:3,xyzzyaaab7,1))*xyzzyaaad7
xyzzyaaaf7=dot_product(dvel(1:3,xyzzyaaaa7,1),dvel(1:3,xyzzyaaab7,1))
retard=retard+(xyzzyaaae7*xyzzyaaag7**2+xyzzyaaaf7)*xyzzyaaag7
enddo
enddo
else
do xyzzyaaaa7=1,netot-1
do xyzzyaaab7=xyzzyaaaa7+1,netot
xyzzyaaac7=eevecsall(1:3,xyzzyaaab7,xyzzyaaaa7)
xyzzyaaag7=sqrt(dot_product(xyzzyaaac7,xyzzyaaac7))
if(xyzzyaaag7/=0.d0)xyzzyaaag7=1.d0/xyzzyaaag7
xyzzyaaad7=dot_product(xyzzyaaac7(1:dimensionality),dvel(1:dimensional&
&ity,xyzzyaaaa7,1))
xyzzyaaae7=dot_product(xyzzyaaac7(1:dimensionality),dvel(1:dimensional&
&ity,xyzzyaaab7,1))*xyzzyaaad7
xyzzyaaaf7=dot_product(dvel(1:dimensionality,xyzzyaaaa7,1),dvel(1:dime&
&nsionality,xyzzyaaab7,1))
retard=retard+(xyzzyaaae7*xyzzyaaag7**2+xyzzyaaaf7)*xyzzyaaag7
enddo
enddo
endif
retard=minus_1_over_2c_squared*retard
end subroutine xyzzyaaai1
end module slaarnace

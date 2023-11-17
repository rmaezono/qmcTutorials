module slaarnabj
use dsp
use slaarnaab,only : tagh_forces,nfterms
use file_utils,   only : open_units
use format_utils, only : wout,i2s
use slaarnabg,     only : naxis_forces,nitot_forces
use parallel,     only : am_master
use run_control,  only : errstop,errwarn
use store,        only : use_future,forces,popstats,dmc_md
implicit none
private
public check_hist_header,write_buffered_data,backtrack_hist,insert_his&
&t_marker
integer,public :: tagh_step,tagh_energy,tagh_k,tagh_t,tagh_fisq,tagh_e&
&wald,tagh_local,tagh_nonlocal,tagh_short,tagh_long,tagh_cppei,tagh_cp&
&pe,tagh_cppee,tagh_masspol,tagh_massvel,tagh_darwinen,tagh_darwinee,t&
&agh_retard,tagh_dipole1,tagh_dipole2,tagh_dipole3,tagh_dipole_sq,tagh&
&_weight,tagh_nconf,tagh_eref,tagh_ebest,tagh_acc,tagh_teff,tagh_esqr,&
&tagh_etotalt,tagh_contact_den,tagh_future1,tagh_future2,tagh_future3,&
&tagh_future4,tagh_future5,tagh_future6,tagh_future7,tagh_future8,tagh&
&_future9,tagh_future10,tagh_future0,tagh_hf_ke,tagh_hf_ex
integer,public :: no_cols_qmc
integer xyzzyaaaa1
contains
subroutine check_hist_header(interaction,netot,nbasis,npcells,atom_bas&
&is_type,isperiodic,constant_energy,qmc_method,have_ppots,have_veep,re&
&lativistic,iaccumulate,newrun,eval_dipole_moment,eval_contact_den,wri&
&teout_hist,mc_twist_av)
use slaarnaaq, only : pair_corr,pair_corr_sph,particle_is_fixed
implicit none
integer,intent(in) :: netot,nbasis,npcells
real(dp),intent(in) :: constant_energy
logical,intent(in) :: isperiodic,have_ppots,have_veep,relativistic,iac&
&cumulate,newrun,eval_dipole_moment,eval_contact_den,writeout_hist,mc_&
&twist_av
character(20),intent(in) :: interaction,atom_basis_type
character(3),intent(in) :: qmc_method
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2,xyzzyaa&
&af2,xyzzyaaag2,xyzzyaaah2
real(dp) xyzzyaaai2
real(dp),parameter :: xyzzyaaaj2=1.d-8
logical,parameter :: xyzzyaaak2=.false.
logical xyzzyaaal2,xyzzyaaam2
logical,save :: xyzzyaaan2=.true.
character(1) temp
character(3) qmc_method_h
character(10) filename
character(20) interaction_h,atom_basis_type_h
character(72) datastring
character(80) title
character(500) checkline
if(.not.am_master)return
if(.not.writeout_hist)then
call xyzzyaaac1(interaction,netot,nbasis,npcells,atom_basis_type,isper&
&iodic,constant_energy,qmc_method,have_ppots,have_veep,relativistic,ia&
&ccumulate,eval_dipole_moment,writeout_hist,eval_contact_den,mc_twist_&
&av)
return
endif
if(xyzzyaaan2)then
call open_units(xyzzyaaaa1,xyzzyaaaa2)
if(xyzzyaaaa2/=0)call errstop('CHECK_HIST_HEADER','Cannot find free io&
& unit.')
xyzzyaaan2=.false.
endif
call xyzzyaaad1(qmc_method,filename)
inquire(file=trim(filename),exist=xyzzyaaal2)
if(.not.xyzzyaaal2)then
if(.not.newrun)call errwarn('CHECK_HIST_HEADER','no '//trim(filename)/&
&/' file exists. Will create one.')
call xyzzyaaac1(interaction,netot,nbasis,npcells,atom_basis_type,isper&
&iodic,constant_energy,qmc_method,have_ppots,have_veep,relativistic,ia&
&ccumulate,eval_dipole_moment,writeout_hist,eval_contact_den,mc_twist_&
&av)
else
if(newrun)then
if(.not.dmc_md)then
if(pair_corr.or.pair_corr_sph)then
if(.not.particle_is_fixed)call errstop('CHECK_HIST_HEADER','Input keyw&
&ord newrun = T but a .hist file already exists. Please delete or rena&
&me it.')
else
call errstop('CHECK_HIST_HEADER','Input keyword newrun = T but a .hist&
& file already exists. Please delete or rename it.')
endif
else
call errstop('CHECK_HIST_HEADER','Input keyword newrun = T but a dmc.h&
&ist file already exists. Please delete or rename it. Note that for ea&
&ch restarted nuclear configuration in a DMC-MD trajectory, then a new&
& dmc.hist file should be started - all the relevant file renaming is &
&done by the runqmcmd script, which we really recommend using.')
endif
endif
open(unit=xyzzyaaaa1,file=trim(filename),status='old',iostat=xyzzyaaaa&
&2)
if(xyzzyaaaa2/=0)call errstop('CHECK_HIST_HEADER','Problem opening ' /&
&/trim(filename)//'.')
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaaa2)checkline
call xyzzyaaaf1(xyzzyaaaa2,filename)
checkline=adjustl(checkline)
if(index(checkline,'Block')>0)call errstop('CHECK_HIST_HEADER','You ap&
&pear to be continuing a CASINO version 1 vmc.hist file.  Please use t&
&he UPDATE_HIST utility to update it to the new format.')
if(index(checkline,'#')==0)call errstop('CHECK_HIST_HEADER','Your data&
& file does not seem to start with a header.  This may be because you &
&are using an old-format file.  If this is the case then please use UP&
&DATE_HIST to update your file.')
rewind(xyzzyaaaa1)
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaaa2)datastring
call xyzzyaaaf1(xyzzyaaaa2,filename)
xyzzyaaac2=index(datastring,'#')
if(xyzzyaaac2>0)then
title=datastring(xyzzyaaac2+1:len_trim(datastring))
else
call errstop('CHECK_HIST_HEADER','Header line does not have a "#" in f&
&ront.  Stopping.')
endif
title=adjustl(title)
if(xyzzyaaak2)call wout('Title: '//trim(title))
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp,xyzzyaaae2
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
if(xyzzyaaak2)call wout('File version number is '//trim(i2s(xyzzyaaae2&
&))//'.')
if(xyzzyaaae2/=1)call errstop('CHECK_HIST_HEADER','Version number of '&
& //trim(filename)//' must be 1.')
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp,qmc_method_h
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
qmc_method_h=adjustl(qmc_method_h)
if(xyzzyaaak2)call wout('The data were generated using '//trim(qmc_met&
&hod_h) //'.')
if(trim(qmc_method_h)/='VMC'.and.trim(qmc_method_h)/='DMC')call errsto&
&p('CHECK_HIST_HEADER','Method in '//trim(filename) //' should be eith&
&er VMC or DMC.')
if(trim(filename)=='vmc.hist'.and.trim(qmc_method_h)/='VMC')call errwa&
&rn('CHECK_HIST_HEADER','you appear to have non-VMC data in a file cal&
&led vmc.hist.')
if(trim(filename)=='dmc.hist'.and.trim(qmc_method_h)/='DMC')call errwa&
&rn('CHECK_HIST_HEADER','you appear to have non-DMC data in a file cal&
&led dmc.hist.')
if(trim(qmc_method_h)/=trim(qmc_method))call errstop('CHECK_HIST_HEADE&
&R','QMC method used to generate data in hist file differs from curren&
&t QMC method.')
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp,interaction_h
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
select case(trim(adjustl(interaction_h)))
case('none','coulomb','ewald','ewaldpp','mpc','ewald_mpc','ewaldpp_mpc&
&','mpc_ewald','mpc_ewaldpp','manual')
continue
case('1')
interaction_h='coulomb'
if(isperiodic)interaction_h='ewald'
case('2')
interaction_h='mpc'
case('3')
interaction_h='ewald_mpc'
case('4')
interaction_h='mpc_ewald'
case default
call errstop('CHECK_HIST_HEADER','Value of INTERACTION/ITERAC in .hist&
& file is not valid.')
end select
if(xyzzyaaak2)call wout('The value of the interaction parameter is ' /&
&/trim(adjustl(interaction_h))//'.')
if(trim(interaction_h)/=trim(interaction))call errstop('CHECK_HIST_HEA&
&DER','INTERACTION is different from the value used in previous QMC si&
&mulation.')
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp,xyzzyaaai2
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
if(xyzzyaaak2)call wout('Have got a constant energy contribution.')
if(abs(xyzzyaaai2-constant_energy)>xyzzyaaaj2)call errstop('CHECK_HIST&
&_HEADER','Constant energy has changed since the previous QMC simulati&
&on.')
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp,xyzzyaaaf2
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
if(xyzzyaaak2)then
if(xyzzyaaaf2/=1)then
call wout('There are '//trim(i2s(xyzzyaaaf2)) //' particles in the sim&
&ulation.')
else
call wout('There is 1 particle in the simulation.')
endif
endif
if(xyzzyaaaf2/=netot)call errstop('CHECK_HIST_HEADER','Number of parti&
&cles has changed since the previous QMC simulation.')
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp,xyzzyaaag2
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
if(xyzzyaaak2)then
if(xyzzyaaag2/=1)then
call wout('The primitive cell contains '//trim(i2s(xyzzyaaag2))//' ato&
&ms.')
else
call wout('The primitive cell contains 1 atom.')
endif
endif
if(xyzzyaaag2/=nbasis)call errstop('CHECK_HIST_HEADER','Number of atom&
&s has changed since the previous QMC simulation.')
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp,xyzzyaaah2
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
if(xyzzyaaak2)then
if(xyzzyaaah2/=1)then
call wout('There are '//trim(i2s(xyzzyaaah2))//' primitive cells.')
else
call wout('There is 1 primitive cell.')
endif
endif
if(xyzzyaaah2/=npcells)call errstop('CHECK_HIST_HEADER','Number of pri&
&m. cells has changed since the previous QMC simulation.')
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp,atom_basis_type_h
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
select case(trim(atom_basis_type_h))
case('none','plane-wave','gaussian','slater-type','numerical','blip','&
&dimer','non_int_he','h2')
continue
case('0')
atom_basis_type_h='none'
case('1')
atom_basis_type_h='plane-wave'
case('2')
atom_basis_type_h='gaussian'
case('3')
atom_basis_type_h='numerical'
case('4')
atom_basis_type_h='blip'
case('5')
atom_basis_type_h='non_int_he'
case default
call errstop('CHECK_HIST_HEADER','Value of ATOM_BASIS_TYPE/BTYPE in .h&
&ist file is not valid.')
end select
if(xyzzyaaak2)call wout('The value of the atom_basis_type parameter is&
& ' //trim(atom_basis_type_h)//'.')
if(trim(atom_basis_type_h)/=trim(atom_basis_type))call errstop('CHECK_&
&HIST_HEADER','ATOM_BASIS_TYPE is different from the value used in pre&
&vious QMC simulation.')
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp,xyzzyaaad2
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
if(xyzzyaaad2==1)then
if(xyzzyaaak2)call wout('The system is periodic.')
elseif(xyzzyaaad2==0)then
if(xyzzyaaak2)call wout('The system is not periodic.')
else
call errstop('CHECK_HIST_HEADER','Periodicity flag must be 0 or 1.')
endif
if((xyzzyaaad2/=1.and.isperiodic).or.(xyzzyaaad2/=0.and..not.isperiodi&
&c))call errstop('CHECK_HIST_HEADER','Periodicity has changed since th&
&e previous QMC simulation.')
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp,no_cols_qmc
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
if(xyzzyaaak2)then
if(no_cols_qmc/=1)then
call wout('There are '//trim(i2s(no_cols_qmc))//' columns of data in '&
& //trim(filename)//'.')
else
call wout('There is 1 column of data in '//trim(filename)//'.')
endif
endif
if(no_cols_qmc<1)call errstop('CHECK_HIST_HEADER','No data in .hist fi&
&le.')
no_cols_qmc=no_cols_qmc+1
tagh_step=1
tagh_energy=-1
tagh_etotalt=-1
tagh_k=-1
tagh_t=-1
tagh_fisq=-1
tagh_ewald=-1
tagh_local=-1
tagh_nonlocal=-1
tagh_short=-1
tagh_long=-1
tagh_cppei=-1
tagh_cppe=-1
tagh_cppee=-1
tagh_masspol=-1
tagh_future0=-1
tagh_future1=-1
tagh_future2=-1
tagh_future3=-1
tagh_future4=-1
tagh_future5=-1
tagh_future6=-1
tagh_future7=-1
tagh_future8=-1
tagh_future9=-1
tagh_future10=-1
tagh_massvel=-1
tagh_darwinen=-1
tagh_darwinee=-1
tagh_retard=-1
tagh_dipole1=-1
tagh_dipole2=-1
tagh_dipole3=-1
tagh_dipole_sq=-1
tagh_weight=-1
tagh_nconf=-1
tagh_eref=-1
tagh_ebest=-1
tagh_acc=-1
tagh_teff=-1
tagh_esqr=-1
if(forces)tagh_forces(1:nfterms,1:naxis_forces,1:nitot_forces)=-1
tagh_contact_den=-1
tagh_hf_ke=-1
tagh_hf_ex=-1
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
do xyzzyaaab2=2,no_cols_qmc
read(xyzzyaaaa1,*,iostat=xyzzyaaaa2)temp,datastring
call xyzzyaaaf1(xyzzyaaaa2,filename)
call xyzzyaaae1(temp)
datastring=adjustl(datastring)
if(trim(datastring)=='ETOT')then
call xyzzyaaag1(tagh_energy)
tagh_energy=xyzzyaaab2
elseif(trim(datastring)=='ETOTALT')then
call xyzzyaaag1(tagh_etotalt)
tagh_etotalt=xyzzyaaab2
elseif(trim(datastring)=='KEI')then
call xyzzyaaag1(tagh_k)
tagh_k=xyzzyaaab2
elseif(trim(datastring)=='TI')then
call xyzzyaaag1(tagh_t)
tagh_t=xyzzyaaab2
elseif(trim(datastring)=='FISQ')then
call xyzzyaaag1(tagh_fisq)
tagh_fisq=xyzzyaaab2
elseif(trim(datastring)=='EWALD')then
call xyzzyaaag1(tagh_ewald)
tagh_ewald=xyzzyaaab2
elseif(trim(datastring)=='LOCAL')then
call xyzzyaaag1(tagh_local)
tagh_local=xyzzyaaab2
elseif(trim(datastring)=='NONLOCAL')then
call xyzzyaaag1(tagh_nonlocal)
tagh_nonlocal=xyzzyaaab2
elseif(trim(datastring)=='SHORT')then
call xyzzyaaag1(tagh_short)
tagh_short=xyzzyaaab2
elseif(trim(datastring)=='LONG')then
call xyzzyaaag1(tagh_long)
tagh_long=xyzzyaaab2
elseif(trim(datastring)=='CPPEI')then
call xyzzyaaag1(tagh_cppei)
tagh_cppei=xyzzyaaab2
elseif(trim(datastring)=='CPPE')then
call xyzzyaaag1(tagh_cppe)
tagh_cppe=xyzzyaaab2
elseif(trim(datastring)=='CPPEE')then
call xyzzyaaag1(tagh_cppee)
tagh_cppee=xyzzyaaab2
elseif(trim(datastring)=='MASSPOL')then
call xyzzyaaag1(tagh_masspol)
tagh_masspol=xyzzyaaab2
elseif(trim(datastring)=='MASSVEL')then
call xyzzyaaag1(tagh_massvel)
tagh_massvel=xyzzyaaab2
elseif(trim(datastring)=='DARWINEN')then
call xyzzyaaag1(tagh_darwinen)
tagh_darwinen=xyzzyaaab2
elseif(trim(datastring)=='DARWINEE')then
call xyzzyaaag1(tagh_darwinee)
tagh_darwinee=xyzzyaaab2
elseif(trim(datastring)=='RETARD')then
call xyzzyaaag1(tagh_retard)
tagh_retard=xyzzyaaab2
elseif(trim(datastring)=='DIPOLE1')then
call xyzzyaaag1(tagh_dipole1)
tagh_dipole1=xyzzyaaab2
elseif(trim(datastring)=='DIPOLE2')then
call xyzzyaaag1(tagh_dipole2)
tagh_dipole2=xyzzyaaab2
elseif(trim(datastring)=='DIPOLE3')then
call xyzzyaaag1(tagh_dipole3)
tagh_dipole3=xyzzyaaab2
elseif(trim(datastring)=='DIPOLESQ')then
call xyzzyaaag1(tagh_dipole_sq)
tagh_dipole_sq=xyzzyaaab2
elseif(trim(datastring)=='WEIGHT')then
call xyzzyaaag1(tagh_weight)
tagh_weight=xyzzyaaab2
elseif(trim(datastring)=='NCONF')then
call xyzzyaaag1(tagh_nconf)
tagh_nconf=xyzzyaaab2
elseif(trim(datastring)=='EREF')then
call xyzzyaaag1(tagh_eref)
tagh_eref=xyzzyaaab2
elseif(trim(datastring)=='EBEST')then
call xyzzyaaag1(tagh_ebest)
tagh_ebest=xyzzyaaab2
elseif(trim(datastring)=='ACC')then
call xyzzyaaag1(tagh_acc)
tagh_acc=xyzzyaaab2
elseif(trim(datastring)=='TEFF')then
call xyzzyaaag1(tagh_teff)
tagh_teff=xyzzyaaab2
elseif(trim(datastring(1:3))=='FOR')then
call xyzzyaaab1(datastring,xyzzyaaab2)
elseif(trim(datastring)=='FUTURE0')then
call xyzzyaaag1(tagh_future0)
tagh_future0=xyzzyaaab2
elseif(trim(datastring)=='FUTURE1')then
call xyzzyaaag1(tagh_future1)
tagh_future1=xyzzyaaab2
elseif(trim(datastring)=='FUTURE2')then
call xyzzyaaag1(tagh_future2)
tagh_future2=xyzzyaaab2
elseif(trim(datastring)=='FUTURE3')then
call xyzzyaaag1(tagh_future3)
tagh_future3=xyzzyaaab2
elseif(trim(datastring)=='FUTURE4')then
call xyzzyaaag1(tagh_future4)
tagh_future4=xyzzyaaab2
elseif(trim(datastring)=='FUTURE5')then
call xyzzyaaag1(tagh_future5)
tagh_future5=xyzzyaaab2
elseif(trim(datastring)=='FUTURE6')then
call xyzzyaaag1(tagh_future6)
tagh_future6=xyzzyaaab2
elseif(trim(datastring)=='FUTURE7')then
call xyzzyaaag1(tagh_future7)
tagh_future7=xyzzyaaab2
elseif(trim(datastring)=='FUTURE8')then
call xyzzyaaag1(tagh_future8)
tagh_future8=xyzzyaaab2
elseif(trim(datastring)=='FUTURE9')then
call xyzzyaaag1(tagh_future9)
tagh_future9=xyzzyaaab2
elseif(trim(datastring)=='FUTURE10')then
call xyzzyaaag1(tagh_future10)
tagh_future10=xyzzyaaab2
elseif(trim(datastring)=='ESQR')then
call xyzzyaaag1(tagh_esqr)
tagh_esqr=xyzzyaaab2
elseif(trim(datastring)=='CONTACT_DEN')then
call xyzzyaaag1(tagh_contact_den)
tagh_contact_den=xyzzyaaab2
elseif(trim(datastring)=='HF_KE')then
call xyzzyaaag1(tagh_hf_ke)
tagh_hf_ke=xyzzyaaab2
elseif(trim(datastring)=='HF_EX')then
call xyzzyaaag1(tagh_hf_ex)
tagh_hf_ex=xyzzyaaab2
else
call errstop('CHECK_HIST_HEADER','Column label not recognised.  Label &
&is: '//trim(datastring)//'.')
endif
enddo
if(xyzzyaaak2)then
call wout('Have read in column labels.')
call wout()
endif
if(tagh_energy<=0)call errwarn('CHECK_HIST_HEADER','total energy data &
&are not present!')
if(tagh_k<=0)call errwarn('CHECK_HIST_HEADER','kinetic energy (K) data&
& are not present!')
if((tagh_short>0.or.tagh_long>0).and..not.(trim(interaction_h)=='mpc'.&
&or.trim(interaction_h)=='mpc_ewald'.or.trim(interaction_h)=='ewald_mp&
&c'))call errwarn('CHECK_HIST_HEADER','MPC data are inexplicably prese&
&nt.')
if(tagh_ewald>0.and.(trim(interaction_h)=='mpc'.or.trim(interaction_h)&
&=='none'))call errwarn('CHECK_HIST_HEADER','Ewald/Coulomb data are in&
&explicably present.')
if(tagh_short>0.and.tagh_long<=0)call errwarn('CHECK_HIST_HEADER','onl&
&y have short-ranged part of MPC interaction.')
if(tagh_short<=0.and.tagh_long>0)call errwarn('CHECK_HIST_HEADER','onl&
&y have long-ranged part of MPC interaction.')
select case(trim(interaction_h))
case('coulomb','ewald','ewald_mpc','mpc_ewald','manual')
if(tagh_ewald<=0)call errwarn('CHECK_HIST_HEADER','non-MPC interaction&
& data are not present!')
end select
select case(trim(interaction_h))
case('mpc','ewald_mpc','mpc_ewald')
if(tagh_short<=0)call errwarn('CHECK_HIST_HEADER','short-ranged MPC da&
&ta are not present!')
if(tagh_long<=0)call errwarn('CHECK_HIST_HEADER','long-ranged MPC data&
& are not present!')
end select
if(mc_twist_av.and.(tagh_hf_ke==-1.or.tagh_hf_ex==-1))call errwarn('CH&
&ECK_HIST_HEADER','Missing HF energies.')
if(trim(qmc_method)=='DMC')then
xyzzyaaam2=.false.
do
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaaa2)datastring
if(xyzzyaaaa2<0)exit
if(xyzzyaaaa2>0)call errstop('CHECK_HIST_HEADER','Problem reading hist&
& file.')
if(trim(adjustl(datastring))=='#### START STATS')then
xyzzyaaam2=.true.
exit
endif
enddo
close(xyzzyaaaa1)
if(iaccumulate)then
if(.not.xyzzyaaam2)call insert_hist_marker('DMC')
else
if(xyzzyaaam2)then
if(pair_corr.or.pair_corr_sph)then
if(.not.particle_is_fixed)call errstop('CHECK_HIST_HEADER','A pre-exis&
&ting dmc.hist file contains statistics accumulation data (after the S&
&TART STATS line) but you have requested a period of DMC equilibration&
&. In such a case it is forbidden to append DMC equilibration data to &
&dmc.hist. Change RUNTYPE in input to dmc_stats to avoid this message.&
&')
else
call errstop('CHECK_HIST_HEADER','A pre-existing dmc.hist file contain&
&s statistics accumulation data (after the START STATS line) but you h&
&ave requested a period of DMC equilibration. In such a case it is for&
&bidden to append DMC equilibration data to dmc.hist. Change RUNTYPE i&
&n input to dmc_stats to avoid this message.')
endif
endif
endif
else
close(xyzzyaaaa1)
endif
endif
end subroutine check_hist_header
subroutine xyzzyaaab1(datastring,i)
implicit none
integer,intent(in) :: i
character(72),intent(in) :: datastring
integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3,xyzzyaaad3
character(1) axis(3),atem(22),iaxis_char,item_char
data axis /'X','Y','Z'/
data atem /'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O'&
&,'P','Q','R','S','T','U','V'/
do xyzzyaaac3=1,nitot_forces
read(datastring(6:),*)xyzzyaaad3
if(xyzzyaaad3==xyzzyaaac3)then
do xyzzyaaaa3=1,naxis_forces
iaxis_char=axis(xyzzyaaaa3)
if(datastring(5:5)==iaxis_char)then
do xyzzyaaab3=1,nfterms
item_char=atem(xyzzyaaab3)
if(datastring(4:4)==item_char)then
call xyzzyaaag1(tagh_forces(xyzzyaaab3,xyzzyaaaa3,xyzzyaaac3))
tagh_forces(xyzzyaaab3,xyzzyaaaa3,xyzzyaaac3)=i
endif
enddo
endif
enddo
endif
enddo
end subroutine xyzzyaaab1
subroutine xyzzyaaac1(interaction,netot,nbasis,npcells,atom_basis_type&
&,isperiodic,constant_energy,qmc_method,have_ppots,have_veep,relativis&
&tic,iaccumulate,eval_dipole_moment,writeout_hist,eval_contact_den,mc_&
&twist_av)
use slaarnace,only : eval_maspol,eval_masvel,eval_darwinen,eval_darwin&
&ee,eval_retard
implicit none
integer,intent(in) :: netot,nbasis,npcells
real(dp),intent(in) :: constant_energy
logical,intent(in) :: isperiodic,have_ppots,have_veep,relativistic,iac&
&cumulate,eval_dipole_moment,writeout_hist,eval_contact_den,mc_twist_a&
&v
character(3),intent(in) :: qmc_method
character(20),intent(in) :: interaction,atom_basis_type
integer xyzzyaaaa4,xyzzyaaab4,xyzzyaaac4,xyzzyaaad4,xyzzyaaae4
character(1)  axis(3),atem(22)
character(10) filename
data  axis/'X','Y','Z'/
data  atem/'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O'&
&,'P','Q','R','S','T','U','V'/
if(.not.am_master)return
xyzzyaaab4=1
tagh_step=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
if(trim(adjustl(qmc_method))=='DMC')then
tagh_weight=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_nconf=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
else
tagh_weight=-1
tagh_nconf=-1
endif
tagh_energy=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
select case(trim(interaction))
case('ewald_mpc','mpc_ewald')
tagh_etotalt=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
case default
tagh_etotalt=-1
end select
if(trim(adjustl(qmc_method))=='DMC')then
tagh_eref=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_ebest=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_acc=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_teff=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
else
tagh_eref=-1
tagh_ebest=-1
tagh_acc=-1
tagh_teff=-1
endif
tagh_esqr=-1
if(popstats)then
tagh_esqr=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
endif
tagh_k=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_t=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_fisq=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
select case(trim(interaction))
case('coulomb','ewald','ewald_mpc','mpc_ewald','manual')
tagh_ewald=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
case default
tagh_ewald=-1
end select
tagh_local=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
if(have_ppots)then
tagh_nonlocal=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
else
tagh_nonlocal=-1
endif
select case(trim(interaction))
case('mpc','ewald_mpc','mpc_ewald')
tagh_short=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_long=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
case default
tagh_short=-1
tagh_long=-1
end select
if(have_veep)then
tagh_cppei=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_cppe=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_cppee=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
else
tagh_cppei=-1
tagh_cppe=-1
tagh_cppee=-1
endif
if(relativistic)then
if(eval_maspol)then
tagh_masspol=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
else
tagh_masspol=-1
endif
if(eval_masvel)then
tagh_massvel=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
else
tagh_massvel=-1
endif
if(eval_darwinen)then
tagh_darwinen=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
else
tagh_darwinen=-1
endif
if(eval_darwinee)then
tagh_darwinee=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
else
tagh_darwinee=-1
endif
if(eval_retard)then
tagh_retard=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
else
tagh_retard=-1
endif
else
tagh_masspol=-1
tagh_massvel=-1
tagh_darwinen=-1
tagh_darwinee=-1
tagh_retard=-1
endif
if(eval_dipole_moment)then
tagh_dipole1=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_dipole2=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_dipole3=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_dipole_sq=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
else
tagh_dipole1=-1
tagh_dipole2=-1
tagh_dipole3=-1
tagh_dipole_sq=-1
endif
if(forces)then
do xyzzyaaac4=1,nitot_forces
do xyzzyaaad4=1,naxis_forces
do xyzzyaaae4=1,nfterms
tagh_forces(xyzzyaaae4,xyzzyaaad4,xyzzyaaac4)=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
enddo
enddo
enddo
endif
if(eval_contact_den)then
tagh_contact_den=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
else
tagh_contact_den=-1
endif
if((trim(adjustl(qmc_method))=='DMC').and.(use_future)) then
tagh_future0=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_future1=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_future2=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_future3=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_future4=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_future5=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_future6=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_future7=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_future8=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_future9=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_future10=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
else
tagh_future0=-1
tagh_future1=-1
tagh_future2=-1
tagh_future3=-1
tagh_future4=-1
tagh_future5=-1
tagh_future6=-1
tagh_future7=-1
tagh_future8=-1
tagh_future9=-1
tagh_future10=-1
endif
if(mc_twist_av)then
tagh_hf_ke=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
tagh_hf_ex=xyzzyaaab4
xyzzyaaab4=xyzzyaaab4+1
endif
no_cols_qmc=xyzzyaaab4-1
if(.not.writeout_hist)return
call xyzzyaaad1(qmc_method,filename)
open(unit=xyzzyaaaa1,file=trim(filename),status='replace',iostat=xyzzy&
&aaaa4)
if(xyzzyaaaa4/=0)call errstop('WRITE_HIST_HEADER','Cannot open ' //tri&
&m(filename)//'.')
write(xyzzyaaaa1,'(a)')'# Title'
write(xyzzyaaaa1,'(a)')'# CASINO raw data.'
write(xyzzyaaaa1,'(a)')'# File version number'
write(xyzzyaaaa1,'(a)')'# 1'
write(xyzzyaaaa1,'(a)')'# QMC method (VMC, DMC, PIMC, AFMC or RMC)'
write(xyzzyaaaa1,'(a)')'# '//trim(qmc_method)
write(xyzzyaaaa1,'(a)')'# Electron-electron interaction type (interact&
&ion keyword)'
write(xyzzyaaaa1,'(a)')'# '//trim(interaction)
write(xyzzyaaaa1,'(a)')'# Constant energy contributions'
write(xyzzyaaaa1,'(a,es26.16)')'# ',constant_energy
write(xyzzyaaaa1,'(a)')'# Number of electrons (and other particles) in&
& simulation'
write(xyzzyaaaa1,'(a)')'# '//trim(i2s(netot))
write(xyzzyaaaa1,'(a)')'# Number of atoms per primitive cell'
write(xyzzyaaaa1,'(a)')'# '//trim(i2s(nbasis))
write(xyzzyaaaa1,'(a)')'# Number of primitive cells'
write(xyzzyaaaa1,'(a)')'# '//trim(i2s(npcells))
write(xyzzyaaaa1,'(a)')'# Basis type (atom_basis_type keyword)'
write(xyzzyaaaa1,'(a)')'# '//trim(atom_basis_type)
write(xyzzyaaaa1,'(a)')'# Periodic (0=NO; 1=YES)'
if(isperiodic)then
write(xyzzyaaaa1,'(a)')'# 1'
else
write(xyzzyaaaa1,'(a)')'# 0'
endif
write(xyzzyaaaa1,'(a)')'# Number of data columns (excluding iteration &
&number)'
write(xyzzyaaaa1,'(a)')'# '//trim(i2s(no_cols_qmc-1))
write(xyzzyaaaa1,'(a)')'# Data items (in order)'
do xyzzyaaab4=2,no_cols_qmc
if(tagh_weight==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# WEIGHT'
elseif(tagh_nconf==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# NCONF'
elseif(tagh_energy==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# ETOT'
elseif(tagh_etotalt==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# ETOTALT'
elseif(tagh_eref==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# EREF'
elseif(tagh_ebest==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# EBEST'
elseif(tagh_acc==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# ACC'
elseif(tagh_teff==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# TEFF'
elseif(tagh_esqr==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# ESQR'
elseif(tagh_k==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# KEI'
elseif(tagh_t==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# TI'
elseif(tagh_fisq==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# FISQ'
elseif(tagh_ewald==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# EWALD'
elseif(tagh_local==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# LOCAL'
elseif(tagh_nonlocal==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# NONLOCAL'
elseif(tagh_short==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# SHORT'
elseif(tagh_long==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# LONG'
elseif(tagh_cppei==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# CPPEI'
elseif(tagh_cppe==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# CPPE'
elseif(tagh_cppee==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# CPPEE'
elseif(tagh_masspol==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# MASSPOL'
elseif(tagh_massvel==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# MASSVEL'
elseif(tagh_darwinen==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# DARWINEN'
elseif(tagh_darwinee==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# DARWINEE'
elseif(tagh_retard==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# RETARD'
elseif(tagh_dipole1==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# DIPOLE1'
elseif(tagh_dipole2==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# DIPOLE2'
elseif(tagh_dipole3==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# DIPOLE3'
elseif(tagh_dipole_sq==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# DIPOLESQ'
elseif(tagh_contact_den==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# CONTACT_DEN'
elseif(tagh_future0==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# FUTURE0'
elseif(tagh_future1==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# FUTURE1'
elseif(tagh_future2==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# FUTURE2'
elseif(tagh_future3==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# FUTURE3'
elseif(tagh_future4==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# FUTURE4'
elseif(tagh_future5==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# FUTURE5'
elseif(tagh_future6==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# FUTURE6'
elseif(tagh_future7==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# FUTURE7'
elseif(tagh_future8==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# FUTURE8'
elseif(tagh_future9==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# FUTURE9'
elseif(tagh_future10==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# FUTURE10'
elseif(forces)then
do xyzzyaaac4=1,nitot_forces
do xyzzyaaad4=1,naxis_forces
do xyzzyaaae4=1,nfterms
if(tagh_forces(xyzzyaaae4,xyzzyaaad4,xyzzyaaac4)==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# FOR'//atem(xyzzyaaae4)//axis(xyzzyaaad4)//tr&
&im(i2s(xyzzyaaac4))
endif
enddo
enddo
enddo
elseif(tagh_hf_ke==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# HF_KE'
elseif(tagh_hf_ex==xyzzyaaab4)then
write(xyzzyaaaa1,'(a)')'# HF_EX'
else
call errstop('WRITE_HIST_HEADER','Bug.')
endif
enddo
write(xyzzyaaaa1,'(a)')'# Raw QMC data'
if(iaccumulate.and.trim(qmc_method)=='DMC')write(xyzzyaaaa1,'(a)') '##&
&## START STATS'
close(xyzzyaaaa1)
end subroutine xyzzyaaac1
subroutine insert_hist_marker(qmc_method)
implicit none
character(3),intent(in) :: qmc_method
integer xyzzyaaaa5
character(10) filename
if(.not.am_master)return
call xyzzyaaad1(qmc_method,filename)
open(unit=xyzzyaaaa1,file=trim(filename),status='unknown',position='ap&
&pend',iostat=xyzzyaaaa5)
if(xyzzyaaaa5/=0)call errstop('INSERT_HIST_MARKER','Cannot open '//tri&
&m(filename)//'.')
write(xyzzyaaaa1,'(a)')'#### START STATS'
close(xyzzyaaaa1)
end subroutine insert_hist_marker
subroutine xyzzyaaad1(qmc_method,filename)
implicit none
character(3),intent(in) :: qmc_method
character(10),intent(out) :: filename
if(trim(qmc_method)=='VMC')then
filename='vmc.hist'
elseif(trim(qmc_method)=='DMC')then
filename='dmc.hist'
else
call errstop('GET_FILENAME','QMC method should be VMC or DMC.')
endif
end subroutine xyzzyaaad1
subroutine write_buffered_data(qmc_method,nmove,qmc_hist_buff)
implicit none
integer,intent(in) :: nmove
real(dp),intent(in) :: qmc_hist_buff(no_cols_qmc,nmove)
character(3),intent(in) :: qmc_method
integer xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7,xyzzyaaad7
character(10) filename
character(40) anumber
character(640) writeme
if(.not.am_master)return
call xyzzyaaad1(qmc_method,filename)
open(xyzzyaaaa1,file=trim(filename),status='unknown',position='append'&
&,iostat=xyzzyaaaa7)
if(xyzzyaaaa7/=0)call errstop('WRITE_BUFFERED_DATA','Problem opening '&
&//trim(filename)//'.')
do xyzzyaaab7=1,nmove
writeme=trim(i2s(nint(qmc_hist_buff(tagh_step,xyzzyaaab7))))
do xyzzyaaac7=2,no_cols_qmc
if(xyzzyaaac7==tagh_nconf)then
anumber=trim(i2s(nint(qmc_hist_buff(xyzzyaaac7,xyzzyaaab7))))
else
write(anumber,*)qmc_hist_buff(xyzzyaaac7,xyzzyaaab7)
endif
writeme=trim(writeme)//' '//trim(anumber)
xyzzyaaad7=modulo(xyzzyaaac7,25)
if(xyzzyaaad7==0)then
write(xyzzyaaaa1,'(a)')trim(writeme)
writeme=''
endif
enddo
if(.not.(xyzzyaaad7==0))write(xyzzyaaaa1,'(a)')trim(writeme)
enddo
close(xyzzyaaaa1)
end subroutine write_buffered_data
subroutine backtrack_hist(qmc_method,nmove)
implicit none
integer,intent(in) :: nmove
character(3),intent(in) :: qmc_method
integer xyzzyaaaa8,xyzzyaaab8
character(10) filename
call xyzzyaaad1(qmc_method,filename)
open(xyzzyaaaa1,file=trim(filename),status='unknown',position='append'&
&,iostat=xyzzyaaab8)
if(xyzzyaaab8/=0)call errstop('BACKTRAC_HIST','Cannot open '//trim(fil&
&ename)//'.')
do xyzzyaaaa8=1,nmove
backspace(xyzzyaaaa1)
enddo
endfile(xyzzyaaaa1)
close(xyzzyaaaa1)
end subroutine backtrack_hist
subroutine xyzzyaaae1(char)
implicit none
character(1),intent(in) :: char
if(char/='#')call errstop('CHECK_HASH','Header line does not have a "#&
&" in front.  Stopping.')
end subroutine xyzzyaaae1
subroutine xyzzyaaaf1(ierr,filename)
implicit none
integer,intent(in) :: ierr
character(*),intent(in) :: filename
if(ierr/=0)call errstop('CHECK_IERR','Problem reading '//trim(filename&
&)//'.')
end subroutine xyzzyaaaf1
subroutine xyzzyaaag1(tag)
implicit none
integer,intent(in) :: tag
if(tag/=-1)call errstop('CHECK_TAG_FREE','Tag assigned twice.  Two col&
&umn labels must be the same.')
end subroutine xyzzyaaag1
end module slaarnabj

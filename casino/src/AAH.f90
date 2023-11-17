module slaarnaah
use parallel
use slaarnaad,    only : write_backflow
use slaarnaao,       only : write_exmol
use file_utils,  only : open_units
use format_utils,only : wout,i2s
use slaarnaas,   only : write_heg,corr_heg_required
use slaarnaau,   only : write_gwfmolorb
use slaarnabg,    only : atom_basis_type
use slaarnabl,     only : write_jastrow
use slaarnabp,        only : write_mdet,modified_mdet
use slaarnabn,    only : write_mahan
use slaarnabu,     only : write_orbmods
use run_control, only : errstop,check_alloc
use store,       only : opt_cycle,use_jastrow,use_backflow,use_orbmods
use slaarnacl,    only : write_stowf_params
use slaarnacs,   only : psi_s
implicit none
private
public read_correlation_header,write_correlation_header,jastrow_in_cor&
&r
integer xyzzyaaaa1,xyzzyaaab1,xyzzyaaac1,xyzzyaaad1
integer,parameter :: xyzzyaaae1=10
logical jastrow_in_corr,xyzzyaaaf1,xyzzyaaag1,xyzzyaaah1,xyzzyaaai1,xy&
&zzyaaaj1,xyzzyaaak1,xyzzyaaal1,xyzzyaaam1,xyzzyaaan1
character(80),allocatable :: xyzzyaaao1(:)
contains
subroutine read_correlation_header(verbose)
implicit none
logical,intent(in) :: verbose
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2
logical xyzzyaaad2
logical,save :: xyzzyaaae2=.true.
character(80) char_80
if(am_master)then
call wout('Reading correlation.data')
call wout('========================')
if(xyzzyaaae2)then
call open_units(xyzzyaaaa1,xyzzyaaaa2)
if(xyzzyaaaa2/=0)call errstop('READ_CORRELATION_HEADER','Cannot find f&
&ree io unit.')
call open_units(xyzzyaaab1,xyzzyaaaa2)
if(xyzzyaaaa2/=0)call errstop('READ_CORRELATION_HEADER','Cannot find f&
&ree io unit <2>.')
allocate(xyzzyaaao1(xyzzyaaae1),stat=xyzzyaaab2)
call check_alloc(xyzzyaaab2,'READ_CORRELATION_HEADER','header_text')
xyzzyaaae2=.false.
endif
xyzzyaaac1=0
xyzzyaaaf1=.false.
xyzzyaaah1=.false.
jastrow_in_corr=.false.
xyzzyaaag1=.false.
xyzzyaaai1=.false.
xyzzyaaad1=1
xyzzyaaak1=.false.
xyzzyaaal1=.false.
xyzzyaaam1=.false.
xyzzyaaan1=.false.
inquire(file='correlation.data',exist=xyzzyaaad2)
call mpi_bcast(xyzzyaaad2,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting file_exists in read_correlation_hea&
&der')
if(.not.xyzzyaaad2)then
call wout('No correlation.data file is present.')
call wout()
return
endif
open(unit=xyzzyaaaa1,file='correlation.data',status='old',iostat=xyzzy&
&aaaa2)
if(xyzzyaaaa2/=0)call errstop('READ_CORRELATION_HEADER','Problem openi&
&ng correlation.data.')
do
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaaa2)char_80
if(xyzzyaaaa2>0)then
call errstop('READ_CORRELATION_HEADER','Problem reading correlation.da&
&ta.')
elseif(xyzzyaaaa2<0)then
if(verbose.and.xyzzyaaac1==0)call wout('No header is present in correl&
&ation.data.')
exit
endif
if(trim(adjustl(char_80))=='START HEADER')then
if(xyzzyaaac1>0)call errstop('READ_CORRELATION_HEADER','There would ap&
&pear to be two headers in correlation.data.')
call wout('Header:')
do xyzzyaaac2=1,xyzzyaaae1
read(xyzzyaaaa1,'(a)',err=10,end=10)xyzzyaaao1(xyzzyaaac2)
xyzzyaaao1(xyzzyaaac2)=adjustl(xyzzyaaao1(xyzzyaaac2))
if(trim(adjustl(xyzzyaaao1(xyzzyaaac2)))=='END HEADER')exit
call wout("  "//trim(xyzzyaaao1(xyzzyaaac2)))
xyzzyaaac1=xyzzyaaac1+1
if(xyzzyaaac1>=xyzzyaaae1)call errstop('READ_CORRELATION_HEADER','Head&
&er is more than ' //trim(i2s(xyzzyaaae1-1))//' long. Try to be more c&
&oncise.')
enddo
xyzzyaaaf1=.true.
elseif(trim(adjustl(char_80))=='START VERSION')then
if(xyzzyaaah1)call errstop('READ_CORRELATION_HEADER','More than one ve&
&rsion number in correlation.data.')
xyzzyaaah1=.true.
read(xyzzyaaaa1,*,err=10,end=10)xyzzyaaad1
call wout('Version of correlation.data : '//trim(i2s(xyzzyaaad1)))
read(xyzzyaaaa1,'(a)',err=10,end=10)char_80
if(trim(adjustl(char_80))/='END VERSION')call errstop('READ_CORRELATIO&
&N_HEADER','Had expected to find "END VERSION".')
elseif(trim(adjustl(char_80))=='START JASTROW')then
if(jastrow_in_corr)call errstop('READ_CORRELATION_HEADER','More than o&
&ne Jastrow factor in correlation.data.')
jastrow_in_corr=.true.
call wout('A Jastrow factor is present in correlation.data.')
elseif(trim(adjustl(char_80))=='START BACKFLOW')then
if(xyzzyaaag1)call errstop('READ_CORRELATION_HEADER','More than one ba&
&ckflow function in correlation.data.')
xyzzyaaag1=.true.
call wout('A backflow function is present in correlation.data.')
elseif(trim(adjustl(char_80))=='START MDET')then
if(xyzzyaaai1)call errstop('READ_CORRELATION_HEADER','More than one MD&
&ET set in correlation.data.')
xyzzyaaai1=.true.
call wout('Multideterminant/excitation data is present in correlation.&
&data.')
elseif(trim(adjustl(char_80))=='START FREE_ORBS')then
if(xyzzyaaaj1)call errstop('READ_CORRELATION_HEADER','More than one FR&
&EE_ORBS set in correlation.data.')
xyzzyaaaj1=.true.
call wout('A FREE_ORBS orbital data block is present in correlation.da&
&ta.')
elseif(trim(adjustl(char_80))=='START ORBMODS')then
if(xyzzyaaak1)call errstop('READ_CORRELATION_HEADER','More than one or&
&bital modifications set in correlation.data.')
xyzzyaaak1=.true.
call wout('An orbital modifications data block is present in correlati&
&on.data.')
elseif(trim(adjustl(char_80))=='START MOLORBMODS')then
if(xyzzyaaal1)call errstop('READ_CORRELATION_HEADER','More than one MO&
&LORBMODS set in correlation.data.')
xyzzyaaal1=.true.
call wout('A MOLORBMODS data block is present in correlation.data.')
elseif(trim(adjustl(char_80))=='START EXMOL')then
if(xyzzyaaam1)call errstop('READ_CORRELATION_HEADER','More than one EX&
&MOL set in correlation.data.')
xyzzyaaam1=.true.
call wout('An EXMOL set is present in correlation.data.')
elseif(trim(adjustl(char_80))=='START MAHAN')then
if(xyzzyaaan1)call errstop('READ_CORRELATION_HEADER', 'More than one M&
&AHAN set in correlation.data.')
xyzzyaaan1=.true.
call wout('A MAHAN set is present in correlation.data.')
endif
enddo
close(xyzzyaaaa1)
if(.not.xyzzyaaaf1)then
xyzzyaaaf1=.true.
xyzzyaaac1=1
xyzzyaaao1(1)='No title given.'
endif
if(.not.xyzzyaaah1)call wout('No version number supplied.')
call wout()
endif
if(am_slave)then
call mpi_bcast(xyzzyaaad2,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting file_exists in read_correlation_hea&
&der')
if(.not.xyzzyaaad2)return
endif
call mpi_bcast(jastrow_in_corr,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting jastrow_in_corr in read_correlation&
&_header')
call mpi_bcast(xyzzyaaaf1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting header_in_corr in read_correlation_&
&header')
call mpi_bcast(xyzzyaaag1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting bf_in_corr in read_correlation_head&
&er')
call mpi_bcast(xyzzyaaah1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting version_in_corr in read_correlation&
&_header')
call mpi_bcast(xyzzyaaai1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting mdet_in_corr in read_correlation_he&
&ader')
call mpi_bcast(xyzzyaaaj1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting heg_in_corr in read_correlation_hea&
&der')
call mpi_bcast(xyzzyaaak1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting orbmods_in_corr in read_correlation&
&_header')
call mpi_bcast(xyzzyaaal1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting molorb_in_corr in read_correlation_&
&header')
call mpi_bcast(xyzzyaaam1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting exmol_in_corr in read_correlation_h&
&eader')
call mpi_bcast(xyzzyaaan1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting mahan_in_corr in read_correlation_h&
&eader')
return
10 call errstop('READ_CORRELATION_HEADER','Problem reading correlation&
&.data.')
end subroutine read_correlation_header
subroutine write_correlation_header(verbose,overwrite)
implicit none
logical,intent(in) :: verbose,overwrite
integer xyzzyaaaa3,xyzzyaaab3
logical xyzzyaaac3
character(20) correlation_name
if(am_master)then
if(.not.(jastrow_in_corr.or.xyzzyaaag1.or.xyzzyaaai1.or.xyzzyaaaj1.or.&
&xyzzyaaak1.or.xyzzyaaal1.or.xyzzyaaam1.or.xyzzyaaan1))return
if(opt_cycle>=0)then
correlation_name='correlation.out.'//trim(i2s(opt_cycle))
else
correlation_name='correlation.out'
endif
if(verbose)then
call wout('Writing parameters to '//trim(correlation_name)//'.')
call wout()
endif
inquire(file=trim(correlation_name),exist=xyzzyaaac3)
if(xyzzyaaac3.and..not.overwrite)then
open(unit=xyzzyaaaa1,file=trim(correlation_name),status='old',position&
&='append',iostat=xyzzyaaab3)
else
open(unit=xyzzyaaaa1,file=trim(correlation_name),status='replace',iost&
&at=xyzzyaaab3)
endif
if(xyzzyaaab3/=0)call errstop('WRITE_CORRELATION_HEADER','Cannot open &
&' //trim(correlation_name)//'.')
write(xyzzyaaaa1,*)'START HEADER'
do xyzzyaaaa3=1,xyzzyaaac1
write(xyzzyaaaa1,*)trim(xyzzyaaao1(xyzzyaaaa3))
enddo
write(xyzzyaaaa1,*)'END HEADER'
write(xyzzyaaaa1,*)
write(xyzzyaaaa1,*)'START VERSION'
write(xyzzyaaaa1,'(3x,a)')trim(i2s(xyzzyaaad1))
write(xyzzyaaaa1,*)'END VERSION'
write(xyzzyaaaa1,*)
close(xyzzyaaaa1)
if(jastrow_in_corr)then
if(use_jastrow)then
call write_jastrow(correlation_name)
else
call xyzzyaaap1('correlation.data',correlation_name,'START JASTROW','E&
&ND JASTROW',.false.)
endif
endif
if(xyzzyaaag1)then
if(use_backflow)then
call write_backflow(correlation_name)
else
call xyzzyaaap1('correlation.data',correlation_name,'START BACKFLOW','&
&END BACKFLOW',.false.)
endif
endif
if(xyzzyaaai1)then
if(modified_mdet)then
call write_mdet(correlation_name)
else
call xyzzyaaap1('correlation.data',correlation_name,'START MDET','END &
&MDET',.false.)
endif
endif
if(xyzzyaaaj1)then
if(corr_heg_required)then
call write_heg(correlation_name)
else
call xyzzyaaap1('correlation.data',correlation_name,'START FREE_ORBS',&
&'END FREE_ORBS',.false.)
endif
endif
if(xyzzyaaak1)then
if(use_orbmods)then
call write_orbmods(correlation_name)
else
call xyzzyaaap1('correlation.data',correlation_name,'START ORBMODS','E&
&ND ORBMODS',.false.)
endif
endif
if(xyzzyaaal1)then
if(use_orbmods)then
if(trim(atom_basis_type)=="gaussian")call write_gwfmolorb(correlation_&
&name)
if(trim(atom_basis_type)=="slater-type")call write_stowf_params(correl&
&ation_name)
else
call xyzzyaaap1('correlation.data',correlation_name,'START MOLORMODS',&
&'END MOLORBMODS',.false.)
endif
endif
if(xyzzyaaam1)then
if(trim(psi_s)=='exmol')then
call write_exmol(correlation_name)
else
call xyzzyaaap1('correlation.data',correlation_name,'START EXMOL','END&
& EXMOL',.false.)
endif
endif
if(xyzzyaaan1)then
if(trim(psi_s)=='mahan_ex')then
call write_mahan(correlation_name)
else
call xyzzyaaap1('correlation.data',correlation_name,'START MAHAN', 'EN&
&D MAHAN',.false.)
endif
endif
endif
end subroutine write_correlation_header
subroutine xyzzyaaap1(infile,outfile,starttext,stoptext,overwrite)
implicit none
logical,intent(in) :: overwrite
character(*),intent(in) :: infile,outfile,starttext,stoptext
integer xyzzyaaaa4
logical xyzzyaaab4,xyzzyaaac4
character(200) char_200
if(.not.am_master)return
open(unit=xyzzyaaaa1,file=trim(adjustl(infile)),status='old',iostat=xy&
&zzyaaaa4)
if(xyzzyaaaa4/=0)call errstop('COPY_SET','Problem opening ' //trim(adj&
&ustl(infile))//'.')
inquire(file=trim(adjustl(outfile)),exist=xyzzyaaab4)
if(xyzzyaaab4.and..not.overwrite)then
open(unit=xyzzyaaab1,file=trim(adjustl(outfile)),status='old',position&
&='append',iostat=xyzzyaaaa4)
else
open(unit=xyzzyaaab1,file=trim(adjustl(outfile)),status='replace',iost&
&at=xyzzyaaaa4)
endif
if(xyzzyaaaa4/=0)call errstop('WRITE_CORRELATION_HEADER','Cannot open &
&' //trim(adjustl(outfile))//'.')
xyzzyaaac4=.false.
do
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaaa4)char_200
if(xyzzyaaaa4<0)exit
if(xyzzyaaaa4>0)call errstop('COPY_SET','Problem reading ' //trim(adju&
&stl(infile))//'.')
if(trim(adjustl(char_200))==trim(adjustl(starttext)))xyzzyaaac4=.true.
if(xyzzyaaac4)write(xyzzyaaab1,'(a)')trim(char_200)
if(trim(adjustl(char_200))==trim(adjustl(stoptext)))exit
enddo
if(xyzzyaaac4)write(xyzzyaaab1,*)
close(xyzzyaaaa1)
close(xyzzyaaab1)
end subroutine xyzzyaaap1
end module slaarnaah

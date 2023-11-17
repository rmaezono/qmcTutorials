!{\src2tex{textfont=tt}}
!!****f* ABINIT/driver
!! NAME
!! driver
!!
!! FUNCTION
!! Driver for ground state, response function, and susceptibility calculations.
!! The present routine drives the following operations.
!! An outer loop allows computation related to different data sets.
!! For each data set, either a GS calculation, a RF calculation,
!! or a SUS calculation is made. In both cases, the input variables
!! are transferred in the proper variables, selected big arrays are
!! allocated, then the gstate, respfn or suscep subroutines are called.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2004 ABINIT group (XG,MKV,MM,MT,FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!! codvsn= code version
!! cpui=initial CPU time
!! dtsets(0:ndtset_alloc)=<type datasets_type>contains all input variables
!! filnam(5)=character strings giving file names
!! filstat=character strings giving name of status file
!! mpi_enreg=informations about MPI parallelization
!! ndtset=number of datasets
!! ndtset_alloc=number of datasets, corrected for allocation of at
!!               least one data set.
!! npsp=number of pseudopotentials
!! pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!! walli=initial wall clock time
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  Input/Output
!! dtfil=<type datafiles_type>infos about file names, file unit numbers
!!  (part of which were initialized previously)
!! results_out(0:ndtset_alloc)=<type results_out_type>contains the results
!!   needed for outvars, including evolving variables
!!   Default values are set up in the calling routine
!!
!! NOTES
!! The array filnam is used for the name of input and output files,
!! and roots for generic input, output or temporary files.
!! Pseudopotential file names are set in pspini and pspatm,
!! using another name. The name filstat will be needed beyond gstate to check
!! the appearance of the "exit" flag, to make a hasty exit, as well as
!! in order to output the status of the computation.
!!
!! TODO
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      chkdilatmx,chkexi,getdim_nloc,gstate,leave_new,mkfilename,mkrdim
!!      nonlinear,pstate,respfn,screening,sigma,status,suscep,timab,wrtout
!!      xredxcart
!!
!! SOURCE

 subroutine driver(codvsn,cpui,dtfil,dtsets,filnam,filstat,&
& mpi_enreg,ndtset,ndtset_alloc,npsp,pspheads,results_out,walli)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: ndtset,ndtset_alloc,npsp
 real(dp),intent(in) :: cpui,walli
 character(len=6),intent(in) :: codvsn
 character(len=fnlen),intent(in) :: filstat
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(inout) :: dtfil
 character(len=fnlen),intent(in) :: filnam(5)
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
 type(pspheader_type),intent(in) :: pspheads(npsp)
 type(results_out_type),intent(inout) :: results_out(0:ndtset_alloc)

!Local variables-------------------------------
 integer,parameter :: level=2
 integer,save :: dimekb_old=-1,lmnmax_old=-1,lnmax_old=-1,mqgrid_old=0
 integer,save :: ntypat_old=-1,paw_size_old=-1,usepaw_old=-1
 integer :: angl_size,angl_size_new,basis_size_new,ceksph,getcell,getocc,getvel
 integer :: getxcart,getxred,idtset,iexit,iget,ii,ij_size_new,ilang,ipsp
 integer :: ipspalch,ireadwf,itypalch,itypat,jdtset,jdtset_status,l_max
 integer :: l_max_new,l_size_max,l_size_max_new,l_size_new,lmn2_size_new
 integer :: lmn_size_new,lmnmax,lmnmaxso,lnmax,lnmaxso,mband,mdtset
 integer :: mesh_size_new,mgfft,mk1mem,mkmem,mkqmem,mpsang,mpssoang,mpw
 integer :: mtypalch,mu,n1xccc,natom,natsph,nfft,nkpt,nkptgw,npspalch,nspden
 integer :: nspinor,nsppol,nsym,ntypalch,ntypat,ntyppure,openexit,optdriver
 integer :: paw_size,prtvol,usepaw,will_read
 real(dp) :: etotal
 logical :: test_new
 character(len=2) :: appen
 character(len=4) :: stringfile
 character(len=500) :: message
 character(len=9) :: stringvar
 character(len=fnlen) :: fildens1in,fildensin,filkss,filscr,fnamewff1
 character(len=fnlen) :: fnamewffddk,fnamewffk,fnamewffq
 type(dataset_type) :: dtset
 type(pawang_type) :: pawang
 type(pseudopotential_type) :: psps
 type(results_gs_type) :: results_gs
 integer :: mkmems(3)
 integer,allocatable :: jdtset_(:),npwtot(:)
 real(dp) :: acell(3),rprim(3,3),rprimd(3,3),rprimdget(3,3),strten(6),tsec(2)
 real(dp),allocatable :: occ(:),vel(:,:),xcart(:,:),xred(:,:),xredget(:,:)
 character(len=fnlen) :: filnam_ds(5)
 type(pawrad_type),allocatable :: pawrad(:)
 type(pawtab_type),allocatable :: pawtab(:)

!******************************************************************

!DEBUG ! Do not comment this line, for the time being XG020913
!write(6,*)' driver : enter '
!stop
!ENDDEBUG

 call timab(100,1,tsec)
 call status(0,filstat,iexit,level,'enter         ')

!Structured debugging if prtvol==-level
 prtvol=dtsets(1)%prtvol
 if(prtvol==-level)then
  write(message,'(80a,a,a)')  ('=',ii=1,80),ch10,&
&   ' driver : enter , debug mode '
  call wrtout(06,message,'COLL')
 end if

 mdtset=99

 if(ndtset>mdtset)then
  write(message, '(a,a,a,a,i2,a,i5)' )ch10,&
&  ' driver : BUG ',ch10,&
&  '  The maximal allowed ndtset is ',mdtset,&
&          ' while the input value is ',ndtset,'.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!DEBUG ! Do not comment this line for the time being XG020913
!write(6,*)' driver : before mtypalch '
!ENDDEBUG

!mtypalch=maxval(dtsets(1:ndtset_alloc)%ntypalch) ! Likely troubles with HP compiler
 mtypalch=dtsets(1)%ntypalch
 do ii=1,ndtset_alloc
  mtypalch=max(dtsets(ii)%ntypalch,mtypalch)
 end do

!DEBUG ! Do not comment this line for the time being XG020913
!write(6,*)' driver : before allocate, npsp= ',npsp
!ENDDEBUG

!Allocation of some arrays independent of the dataset
 allocate(psps%filpsp(npsp))
 allocate(psps%pspcod(npsp))
 allocate(psps%pspdat(npsp))
 allocate(psps%pspso(npsp))
 allocate(psps%pspxc(npsp))
 allocate(psps%title(npsp))
 allocate(psps%zionpsp(npsp))
 allocate(psps%znuclpsp(npsp))

 psps%filpsp(1:npsp)=pspheads(1:npsp)%filpsp
 psps%pspcod(1:npsp)=pspheads(1:npsp)%pspcod
 psps%pspdat(1:npsp)=pspheads(1:npsp)%pspdat
 psps%pspso(1:npsp)=pspheads(1:npsp)%pspso
 psps%pspxc(1:npsp)=pspheads(1:npsp)%pspxc
 psps%title(1:npsp)=pspheads(1:npsp)%title
 psps%zionpsp(1:npsp)=pspheads(1:npsp)%zionpsp
 psps%znuclpsp(1:npsp)=pspheads(1:npsp)%znuclpsp

 allocate(jdtset_(0:ndtset))
 if(ndtset/=0)then
  jdtset_(:)=dtsets(0:ndtset)%jdtset
 else
  jdtset_(0)=0
 end if

!*********************************************************************
!Big loop on datasets

!Do loop on idtset (allocate statements are present)
 do idtset=1,ndtset_alloc

  jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=1

  if(ndtset>=2)then
   jdtset_status=jdtset
  else
   jdtset_status=0
  end if

  call status(jdtset_status,filstat,iexit,level,'loop jdtset   ')

  write(message,'(a,80a,a,a,i2,a,66a,a)') ch10,&
&  ('=',mu=1,80),ch10,&
&  '== DATASET ',jdtset,' ',('=',mu=1,66),ch10
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,message,'PERS')     ! PERS is choosen to make debugging easier

  dtset%jdtset=jdtset

! Copy input values
  acell(:)  =dtsets(idtset)%acell_orig(:)
  dtset%acell_orig(:)=dtsets(idtset)%acell_orig(:)

! Determine here wether the calculation is PAW
  usepaw  =0
  if (pspheads(1)%pspcod==7) usepaw=1  ! If paw, all pspcod necessarily are 7 (see iofn2)

! Note that angdeg_orig is NOT YET initialized

  dtset%accesswff=dtsets(idtset)%accesswff
  dtset%cpus     =dtsets(idtset)%cpus
  dtset%diecut   =dtsets(idtset)%diecut
  dtset%dielng   =dtsets(idtset)%dielng
  dtset%diemac   =dtsets(idtset)%diemac
  dtset%diemix   =dtsets(idtset)%diemix
  dtset%diegap   =dtsets(idtset)%diegap
  dtset%dielam   =dtsets(idtset)%dielam

  dtset%td_maxene=dtsets(idtset)%td_maxene

  dtset%gwcalctyp=dtsets(idtset)%gwcalctyp
  dtset%idyson   =dtsets(idtset)%idyson
  dtset%ndyson   =dtsets(idtset)%ndyson
  dtset%intexact =dtsets(idtset)%intexact
  dtset%nbandsus =dtsets(idtset)%nbandsus
  dtset%ldgapp   =dtsets(idtset)%ldgapp

  dtset%boxcutmin=dtsets(idtset)%boxcutmin
  dtset%bxctmindg=dtsets(idtset)%bxctmindg
  dtset%charge   =dtsets(idtset)%charge
  dtset%dedlnn   =dtsets(idtset)%dedlnn
  dtset%dosdeltae=dtsets(idtset)%dosdeltae
  dtset%dtion    =dtsets(idtset)%dtion
  dtset%ecut     =dtsets(idtset)%ecut
  dtset%sciss    =dtsets(idtset)%sciss
  dtset%tsmear   =dtsets(idtset)%tsmear
  dtset%vis      =dtsets(idtset)%vis
  dtset%ecutsm   =dtsets(idtset)%ecutsm
  dtset%effmass  =dtsets(idtset)%effmass
  dtset%dilatmx  =dtsets(idtset)%dilatmx
  dtset%stmbias  =dtsets(idtset)%stmbias
  dtset%strfact  =dtsets(idtset)%strfact
  dtset%freqsusin=dtsets(idtset)%freqsusin
  dtset%freqsuslo=dtsets(idtset)%freqsuslo
  dtset%mdftemp  =dtsets(idtset)%mdftemp
  dtset%mditemp  =dtsets(idtset)%mditemp
  dtset%noseinert=dtsets(idtset)%noseinert
  dtset%optforces=dtsets(idtset)%optforces
  dtset%optnlxccc=dtsets(idtset)%optnlxccc
  dtset%tphysel  =dtsets(idtset)%tphysel
  dtset%kptrlen  =dtsets(idtset)%kptrlen
  dtset%strtarget(1:6)=dtsets(idtset)%strtarget(1:6)
  dtset%strprecon=dtsets(idtset)%strprecon
  dtset%friction =dtsets(idtset)%friction
  dtset%mdwall   =dtsets(idtset)%mdwall
  dtset%fixmom   =dtsets(idtset)%fixmom
  dtset%eshift   =dtsets(idtset)%eshift
  dtset%boxcenter(1:3)=dtsets(idtset)%boxcenter(1:3)
  dtset%ecuteps  =dtsets(idtset)%ecuteps
  dtset%ecutsigx =dtsets(idtset)%ecutsigx
  dtset%ecutwfn  =dtsets(idtset)%ecutwfn
  dtset%ppmfrq   =dtsets(idtset)%ppmfrq
  dtset%zcut     =dtsets(idtset)%zcut
  dtset%soenergy =dtsets(idtset)%soenergy
  dtset%nomegasrd=dtsets(idtset)%nomegasrd
  dtset%omegasrdmax=dtsets(idtset)%omegasrdmax
  dtset%tfkinfunc=dtsets(idtset)%tfkinfunc
  dtset%tfnewton =dtsets(idtset)%tfnewton

  dtset%getocc   =dtsets(idtset)%getocc
  dtset%getwfk   =dtsets(idtset)%getwfk
  dtset%getxcart =dtsets(idtset)%getxcart
  dtset%getxred  =dtsets(idtset)%getxred
  dtset%getden   =dtsets(idtset)%getden
  dtset%getcell  =dtsets(idtset)%getcell
  dtset%getwfq   =dtsets(idtset)%getwfq
  dtset%get1wf   =dtsets(idtset)%get1wf
  dtset%getddk   =dtsets(idtset)%getddk
  dtset%getvel   =dtsets(idtset)%getvel
  dtset%getkss   =dtsets(idtset)%getkss
  dtset%getscr   =dtsets(idtset)%getscr
  dtset%get1den  =dtsets(idtset)%get1den

! Note that kptopt is initialized in both invars1.f and
! invars2.f. So, one does not need dtsets(idtset)%kptopt and kptopt_
! kptopt    =kptopt_(idtset)

  dtset%kptnrm    =dtsets(idtset)%kptnrm
  dtset%dsifkpt(:)=dtsets(idtset)%dsifkpt(:)
  dtset%kptrlatt(:,:)=dtsets(idtset)%kptrlatt(:,:)
  dtset%genafm(:) =dtsets(idtset)%genafm(:)
  mband           =dtsets(idtset)%mband
  dtset%mband     =dtsets(idtset)%mband
  dtset%mffmem    =dtsets(idtset)%mffmem
  mgfft           =dtsets(idtset)%mgfft
  dtset%mgfft     =dtsets(idtset)%mgfft

  dtset%mkmem  =dtsets(idtset)%mkmem
  dtset%mkqmem =dtsets(idtset)%mkqmem
  dtset%mk1mem =dtsets(idtset)%mk1mem
  mkmems(1)=dtset%mkmem
  mkmems(2)=dtset%mkqmem
  mkmems(3)=dtset%mk1mem

  mpw             =dtsets(idtset)%mpw
  dtset%mpw       =dtsets(idtset)%mpw

  psps%mqgrid     =dtsets(idtset)%mqgrid
  if (usepaw==1) psps%mqgrid=dtsets(idtset)%pawmqgrid
  psps%optnlxccc  =dtsets(idtset)%optnlxccc

  natom           =dtsets(idtset)%natom
  dtset%natom     =dtsets(idtset)%natom

  natsph          =dtsets(idtset)%natsph
  dtset%natsph    =dtsets(idtset)%natsph

  dtset%nconeq    =dtsets(idtset)%nconeq

 ! CASINO
  dtset%ndtset    =ndtset
 ! CASINO

  dtset%nelect    =dtsets(idtset)%nelect

  nfft            =dtsets(idtset)%nfft
  dtset%nfft      =dtsets(idtset)%nfft

  dtset%ngfft(:)  =dtsets(idtset)%ngfft(:)
  nkptgw           =dtsets(idtset)%nkptgw
  dtset%nkptgw     =dtsets(idtset)%nkptgw

  nkpt            =dtsets(idtset)%nkpt
  dtset%nkpt      =dtsets(idtset)%nkpt

  dtset%nloalg(:) =dtsets(idtset)%nloalg(:)

  dtset%nline     =dtsets(idtset)%nline
  dtset%nnsclo    =dtsets(idtset)%nnsclo
  dtset%nstep     =dtsets(idtset)%nstep
  dtset%ntime     =dtsets(idtset)%ntime
  dtset%nfreqsus  =dtsets(idtset)%nfreqsus

  dtset%npsp      =npsp

  npspalch        =dtsets(idtset)%npspalch
  dtset%npspalch  =dtsets(idtset)%npspalch

  dtset%nshiftk   =dtsets(idtset)%nshiftk

  nspden          =dtsets(idtset)%nspden
  dtset%nspden    =dtsets(idtset)%nspden

  nspinor         =dtsets(idtset)%nspinor
  dtset%nspinor   =dtsets(idtset)%nspinor

  nsppol          =dtsets(idtset)%nsppol
  dtset%nsppol    =dtsets(idtset)%nsppol

  nsym            =dtsets(idtset)%nsym
  dtset%nsym      =dtsets(idtset)%nsym

  ntypalch        =dtsets(idtset)%ntypalch
  dtset%ntypalch  =dtsets(idtset)%ntypalch

  ntypat           =dtsets(idtset)%ntypat
  dtset%ntypat     =dtsets(idtset)%ntypat

  ntyppure        =dtsets(idtset)%ntyppure
  dtset%ntyppure  =dtsets(idtset)%ntyppure

  dtset%occopt    =dtsets(idtset)%occopt

  dtset%ceksph   =dtsets(idtset)%ceksph
  dtset%enunit   =dtsets(idtset)%enunit
  dtset%exchn2n3 =dtsets(idtset)%exchn2n3
  dtset%frzfermi =dtsets(idtset)%frzfermi
  dtset%ionmov   =dtsets(idtset)%ionmov
  dtset%intxc    =dtsets(idtset)%intxc
  dtset%iprcch   =dtsets(idtset)%iprcch
  dtset%iprcel   =dtsets(idtset)%iprcel
  dtset%iprcfc   =dtsets(idtset)%iprcfc
  dtset%irdwfk   =dtsets(idtset)%irdwfk
  dtset%iscf     =dtsets(idtset)%iscf
  dtset%isecur   =dtsets(idtset)%isecur
  dtset%ixc      =dtsets(idtset)%ixc
  dtset%nqpt     =dtsets(idtset)%nqpt
  dtset%restartxf=dtsets(idtset)%restartxf
  dtset%optcell  =dtsets(idtset)%optcell
  dtset%irdwfq   =dtsets(idtset)%irdwfq
  dtset%ird1wf   =dtsets(idtset)%ird1wf
  dtset%irdddk   =dtsets(idtset)%irdddk
  dtset%irdscr   =dtsets(idtset)%irdscr
  dtset%irdkss   =dtsets(idtset)%irdkss
  dtset%kptopt   =dtsets(idtset)%kptopt
  dtset%chkexit  =dtsets(idtset)%chkexit
  dtset%ikhxc    =dtsets(idtset)%ikhxc
  dtset%nbdbuf   =dtsets(idtset)%nbdbuf
  dtset%localrdwf=dtsets(idtset)%localrdwf
  dtset%efield(1:3)=dtsets(idtset)%efield(1:3)
  dtset%nberry   =dtsets(idtset)%nberry
  dtset%bdberry(1:4)=dtsets(idtset)%bdberry(1:4)
  dtset%delayperm=dtsets(idtset)%delayperm
  dtset%signperm =dtsets(idtset)%signperm
  dtset%nbandkss  =dtsets(idtset)%nbandkss
  dtset%npwkss  =dtsets(idtset)%npwkss
  dtset%berryopt =dtsets(idtset)%berryopt
  dtset%wfoptalg =dtsets(idtset)%wfoptalg
  dtset%nbdblock =dtsets(idtset)%nbdblock
  dtset%ngfftdg(:)=dtsets(idtset)%ngfftdg(:)
  dtset%kssform  =dtsets(idtset)%kssform
  dtset%pawecutdg=dtsets(idtset)%pawecutdg
  dtset%pawlcutd =dtsets(idtset)%pawlcutd
  dtset%pawmixtyp=dtsets(idtset)%pawmixtyp
  dtset%pawntheta=dtsets(idtset)%pawntheta
  dtset%pawnphi  =dtsets(idtset)%pawnphi
  dtset%pawnzlm  =dtsets(idtset)%pawnzlm
  dtset%pawsphmix=dtsets(idtset)%pawsphmix
  dtset%pawvlbox =dtsets(idtset)%pawvlbox
  dtset%pawxcdev =dtsets(idtset)%pawxcdev
  dtset%useylm   =dtsets(idtset)%useylm
  dtset%td_mexcit=dtsets(idtset)%td_mexcit
  dtset%npweps   =dtsets(idtset)%npweps
  dtset%npwsigx   =dtsets(idtset)%npwsigx
  dtset%npwwfn   =dtsets(idtset)%npwwfn
  dtset%nsheps   =dtsets(idtset)%nsheps
  dtset%nshsigx   =dtsets(idtset)%nshsigx
  dtset%nshwfn   =dtsets(idtset)%nshwfn
  dtset%spgroup  =dtsets(idtset)%spgroup
  dtset%ptgroupma=dtsets(idtset)%ptgroupma
  dtset%nctime   =dtsets(idtset)%nctime

  optdriver       =dtsets(idtset)%optdriver
  dtset%optdriver =dtsets(idtset)%optdriver

  dtset%ortalg    =dtsets(idtset)%ortalg

  dtset%prepanl   =dtsets(idtset)%prepanl
  dtset%parareel  =dtsets(idtset)%parareel
  dtset%ecutgros  =dtsets(idtset)%ecutgros
  dtset%npara     =dtsets(idtset)%npara
  dtset%kpara     =dtsets(idtset)%kpara
  dtset%npack     =dtsets(idtset)%npack

  dtset%prteig    =dtsets(idtset)%prteig
  dtset%prtvol    =dtsets(idtset)%prtvol
  dtset%prtden    =dtsets(idtset)%prtden
  dtset%prtpot    =dtsets(idtset)%prtpot
  dtset%prtdos    =dtsets(idtset)%prtdos
  dtset%prtfsurf  =dtsets(idtset)%prtfsurf
  dtset%prtgeo    =dtsets(idtset)%prtgeo
  dtset%prtcml    =dtsets(idtset)%prtcml
  dtset%prtstm    =dtsets(idtset)%prtstm
  dtset%prt1dm    =dtsets(idtset)%prt1dm
  dtset%prtvha    =dtsets(idtset)%prtvha
  dtset%prtvhxc   =dtsets(idtset)%prtvhxc
  dtset%prtvxc    =dtsets(idtset)%prtvxc
  dtset%prtwf     =dtsets(idtset)%prtwf
  dtset%prtbbb    =dtsets(idtset)%prtbbb

  dtset%qprtrb(:) =dtsets(idtset)%qprtrb(:)

  dtset%qpt(:)    =dtsets(idtset)%qpt(:)
  dtset%qptnrm    =dtsets(idtset)%qptnrm
  dtset%qptn(:)   =dtsets(idtset)%qptn(:)

  dtset%ratsph    =dtsets(idtset)%ratsph
  dtset%rfasr     =dtsets(idtset)%rfasr
  dtset%rfatpol(:)=dtsets(idtset)%rfatpol(:)
  dtset%rfdir(:)  =dtsets(idtset)%rfdir(:)
  dtset%rfelfd    =dtsets(idtset)%rfelfd
  dtset%rfmeth    =dtsets(idtset)%rfmeth
  dtset%rfphon    =dtsets(idtset)%rfphon
  dtset%rfstrs    =dtsets(idtset)%rfstrs
  dtset%rfthrd    =dtsets(idtset)%rfthrd
  dtset%rfuser    =dtsets(idtset)%rfuser

  dtset%rf1atpol(:)=dtsets(idtset)%rf1atpol(:)
  dtset%rf1dir(:)  =dtsets(idtset)%rf1dir(:)
  dtset%rf1elfd    =dtsets(idtset)%rf1elfd
  dtset%rf1phon    =dtsets(idtset)%rf1phon

  dtset%rf2atpol(:)=dtsets(idtset)%rf2atpol(:)
  dtset%rf2dir(:)  =dtsets(idtset)%rf2dir(:)
  dtset%rf2elfd    =dtsets(idtset)%rf2elfd
  dtset%rf2phon    =dtsets(idtset)%rf2phon

  dtset%rf3atpol(:)=dtsets(idtset)%rf3atpol(:)
  dtset%rf3dir(:)  =dtsets(idtset)%rf3dir(:)
  dtset%rf3elfd    =dtsets(idtset)%rf3elfd
  dtset%rf3phon    =dtsets(idtset)%rf3phon

  rprim(:,:)      =dtsets(idtset)%rprim_orig(:,:)
  dtset%rprim_orig(:,:)=dtsets(idtset)%rprim_orig(:,:)
  dtset%rprimd_orig(:,:)=dtsets(idtset)%rprimd_orig(:,:)

  dtset%tolmxf    =dtsets(idtset)%tolmxf
  dtset%tolwfr    =dtsets(idtset)%tolwfr
  dtset%toldff    =dtsets(idtset)%toldff
  dtset%toldfe    =dtsets(idtset)%toldfe
  dtset%toldet    =dtsets(idtset)%toldet
  dtset%tolvrs    =dtsets(idtset)%tolvrs

  dtset%useria    =dtsets(idtset)%useria
  dtset%userib    =dtsets(idtset)%userib
  dtset%useric    =dtsets(idtset)%useric
  dtset%userid    =dtsets(idtset)%userid
  dtset%userie    =dtsets(idtset)%userie

  dtset%userra    =dtsets(idtset)%userra
  dtset%userrb    =dtsets(idtset)%userrb
  dtset%userrc    =dtsets(idtset)%userrc
  dtset%userrd    =dtsets(idtset)%userrd
  dtset%userre    =dtsets(idtset)%userre

  dtset%vprtrb(:) =dtsets(idtset)%vprtrb(:)

  allocate(dtset%algalch(ntypalch))
  dtset%algalch(:)=dtsets(idtset)%algalch(1:ntypalch)
  allocate(dtset%amu(ntypat))
  dtset%amu   (:)  =dtsets(idtset)%amu(1:ntypat)
  allocate(dtset%bdgw(2,nkptgw))
  dtset%bdgw  (:,:)=dtsets(idtset)%bdgw(:,1:nkptgw)
  allocate(dtset%densty(ntypat,4))
  dtset%densty(:,:)=dtsets(idtset)%densty(1:ntypat,:)
  allocate(dtset%iatfix(3,natom))
  dtset%iatfix(:,:)=dtsets(idtset)%iatfix(:,1:natom)
  allocate(dtset%iatsph(natsph))
  dtset%iatsph(:)=dtsets(idtset)%iatsph(1:natsph)
  allocate(dtset%istwfk(nkpt))
  dtset%istwfk(:)  =dtsets(idtset)%istwfk(1:nkpt)
  allocate(dtset%kberry(3,dtset%nberry))
  dtset%kberry(:,:)=dtsets(idtset)%kberry(1:3,1:dtset%nberry)
  allocate(dtset%kptgw(3,nkptgw))
  dtset%kptgw (:,:)=dtsets(idtset)%kptgw(:,1:nkptgw)
  allocate(dtset%kpt(3,nkpt))
  dtset%kpt   (:,:)=dtsets(idtset)%kpt(:,1:nkpt)
  allocate(dtset%kptns(3,nkpt))
  dtset%kptns (:,:)=dtsets(idtset)%kptns(:,1:nkpt)
  allocate(dtset%mixalch(npspalch,ntypalch))
  dtset%mixalch(:,:)=dtsets(idtset)%mixalch(1:npspalch,1:ntypalch)
  allocate(dtset%nband(nkpt*nsppol))
  dtset%nband (:)  =dtsets(idtset)%nband(1:nkpt*nsppol)
  allocate(dtset%occ_orig(mband*nkpt*nsppol))
  dtset%occ_orig(:)=dtsets(idtset)%occ_orig(1:mband*nkpt*nsppol)
  allocate(dtset%so_typat(ntypat))
  dtset%so_typat(:) =dtsets(idtset)%so_typat(1:ntypat)
  allocate(dtset%shiftk(3,dtset%nshiftk))
  dtset%shiftk(:,:)=dtsets(idtset)%shiftk(:,1:dtset%nshiftk)
  allocate(dtset%spinat(3,natom))
  dtset%spinat(:,:)=dtsets(idtset)%spinat(:,1:natom)
  allocate(dtset%symafm(nsym))
  dtset%symafm(:)  =dtsets(idtset)%symafm(1:nsym)
  allocate(dtset%symrel(3,3,nsym))
  dtset%symrel(:,:,:)=dtsets(idtset)%symrel(:,:,1:nsym)
  allocate(dtset%tnons(3,nsym))
  dtset%tnons(:,:) =dtsets(idtset)%tnons(:,1:nsym)
  allocate(dtset%typat(natom))
  dtset%typat  (:)  =dtsets(idtset)%typat(1:natom)
  allocate(dtset%vel_orig(3,natom))
  dtset%vel_orig(:,:)=dtsets(idtset)%vel_orig(:,1:natom)

  allocate(dtset%wtatcon(3,natom,dtset%nconeq))
  dtset%wtatcon(:,:,:)=dtsets(idtset)%wtatcon(:,1:natom,1:dtset%nconeq)

  allocate(dtset%wtk(nkpt))
  dtset%wtk   (:)  =dtsets(idtset)%wtk(1:nkpt)
  allocate(dtset%xred_orig(3,natom))
  dtset%xred_orig(:,:)=dtsets(idtset)%xred_orig(:,1:natom)

  allocate(dtset%ziontypat(ntypat))
  dtset%ziontypat(:)  =dtsets(idtset)%ziontypat(1:ntypat)

  allocate(dtset%znucl(npsp))
  dtset%znucl(:)  =dtsets(idtset)%znucl(1:npsp)

! Allocate arrays
  allocate(occ(mband*nkpt*nsppol))
  allocate(vel(3,natom) )
  allocate(xred(3,natom))

  occ   (:)  =dtsets(idtset)%occ_orig(1:mband*nkpt*nsppol)
  vel   (:,:)=dtsets(idtset)%vel_orig(:,1:natom)
  xred  (:,:)=dtsets(idtset)%xred_orig(:,1:natom)

!****************************************************************************
! Treat the file names (get variables)


  filnam_ds(1:5)=filnam(1:5)

! If multi dataset mode, special treatment of filenames 3 and 4 (density and
! wavefunctions input and output, as well as other output files)
  if(ndtset>0)then
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   filnam_ds(3)=trim(filnam(3))//'_DS'//appen
   filnam_ds(4)=trim(filnam(4))//'_DS'//appen
!DEBUG
!  write(6,*)' filnam_ds(3)',filnam_ds(3)
!  write(6,*)' filnam_ds(4)',filnam_ds(4)
!ENDDEBUG
  end if

! According to getwfk and irdwfk, build _WFK file name, referred as fnamewffk
  stringfile='_WFK' ; stringvar='wfk'
  call mkfilename(filnam,fnamewffk,dtset%getwfk,idtset,dtset%irdwfk,jdtset_,&
&  ndtset,stringfile,stringvar,will_read)

  if(optdriver/=1)ireadwf=will_read
  if(ndtset/=0 .and. optdriver==1 .and. will_read==0)then
   write(message, '(a,a,a,a,a,a,a,a,i3,a,a,a,i3,a,i3,a,a,a)' )ch10,&
&   ' driver : ERROR -',ch10,&
&   '  At least one of the input variables irdwfk and getwfk ',ch10,&
&   '  must refer to a valid _WFK file, in the response function',ch10,&
&   '  case, while for idtset=',idtset,',',ch10,&
&   '  they are irdwfk=',dtset%irdwfk,', and getwfk=',dtset%getwfk,'.',ch10,&
&   '  Action : correct irdwfk or getwfk in your input file.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

! Treatment of the other get wavefunction variable, if response function case
! or nonlinear case
  if ((optdriver==1).or.(optdriver==5)) then

!  According to getwfq and irdwfq, build _WFQ file name, referred as fnamewffq
   stringfile='_WFQ' ; stringvar='wfq'
   call mkfilename(filnam,fnamewffq,dtset%getwfq,idtset,dtset%irdwfq,jdtset_,&
&   ndtset,stringfile,stringvar,will_read)
!  If fnamewffq is not initialized thanks to getwfq or irdwfq, use fnamewffk
   if(will_read==0)fnamewffq=fnamewffk

!  According to get1wf and ird1wf, build _1WF file name, referred as fnamewff1
   stringfile='_1WF' ; stringvar='1wf'
   call mkfilename(filnam,fnamewff1,dtset%get1wf,idtset,dtset%ird1wf,jdtset_,&
&   ndtset,stringfile,stringvar,will_read)
   ireadwf=will_read

!  According to getddk and irdddk, build _1WF file name, referred as fnamewffddk
   stringfile='_1WF' ; stringvar='ddk'
   call mkfilename(filnam,fnamewffddk,dtset%getddk,idtset,dtset%irdddk,jdtset_,&
&   ndtset,stringfile,stringvar,will_read)

  end if ! optdriver is 1 or 5

! According to getden, build _DEN file name, referred as fildensin
! A default is available if getden is 0
  stringfile='_DEN' ; stringvar='den'
  call mkfilename(filnam,fildensin,dtset%getden,idtset,0,jdtset_,&
&  ndtset,stringfile,stringvar,will_read)
  if(will_read==0)fildensin=trim(filnam_ds(3))//'_DEN'

! According to get1den, build _DEN file name, referred as fildens1in
! A default is available if get1den is 0
  stringfile='_DEN' ; stringvar='1den'
  call mkfilename(filnam,fildens1in,dtset%get1den,idtset,0,jdtset_,&
&  ndtset,stringfile,stringvar,will_read)
  if(will_read==0)fildens1in=trim(filnam_ds(3))//'_DEN'

  if (optdriver==4) then

!  According to getscr, build _SCR file name, referred as filscr
!  A default is available if getscr is 0
   stringfile='_SCR' ; stringvar='scr'
   call mkfilename(filnam,filscr,dtset%getscr,idtset,dtset%irdscr,jdtset_,&
&   ndtset,stringfile,stringvar,will_read)
   if(will_read==0)filscr=trim(filnam_ds(3))//'_SCR'

  end if

  if ((optdriver==3).or.(optdriver==4)) then

!  According to getkss, build _KSS file name, referred as filkss
!  A default is available if getkss is 0
   stringfile='_KSS' ; stringvar='kss'
   call mkfilename(filnam,filkss,dtset%getkss,idtset,dtset%irdkss,jdtset_,&
&   ndtset,stringfile,stringvar,will_read)
   if(will_read==0)filkss=trim(filnam_ds(3))//'_KSS'

  end if

  dtfil%ireadwf=ireadwf
  dtfil%filnam_ds(1:5)=filnam_ds(1:5)
  dtfil%filscr        =filscr
  dtfil%fildensin     =fildensin
  dtfil%fildens1in    =fildens1in
  dtfil%filkss        =filkss
  dtfil%filstat       =filstat
  dtfil%fnamewffk     =fnamewffk
  dtfil%fnamewffq     =fnamewffq
  dtfil%fnamewffddk   =fnamewffddk
  dtfil%fnamewff1     =fnamewff1

!****************************************************************************
! Treat other get variables

! If multi dataset mode, and already the second dataset,
! treatment of other get variables.
  if( ndtset>1 .and. idtset>1 )then

   getocc=dtset%getocc
   getvel=dtset%getvel
   getxcart=dtset%getxcart
   getxred=dtset%getxred
   getcell=dtset%getcell
!  Should be in chkinp ; should also check that the number of atoms coincide.
   if(getxcart/=0 .and. getxred/=0) then
    write(message, '(a,a,a,a,a,a,i3,a,a,a,i3,a,i3,a,a,a)' )ch10,&
&    ' driver : ERROR -',ch10,&
&    '  The input variables getxcart and getxred cannot be',ch10,&
&    '  simultaneously non-zero, while for idtset=',idtset,',',ch10,&
&    '  they are getxcart=',getxcart,', and getxred=',getxred,'.',ch10,&
&    '  Action : correct getxcart or getxred in your input file.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

!DEBUG
!  write(6,*)' ndtset,idtset,getxred=',ndtset,idtset,getxred
!ENDDEBUG

   if(getocc>0 .or. (getocc<0 .and. idtset+getocc>0) )then
!   In case getocc is a negative number (so must add to idtset)
    if(getocc<0 .and. idtset+getocc>0) iget=idtset+getocc
    if(getocc>0)then
     do iget=1,idtset
      if( dtsets(iget)%jdtset==getocc )exit
     end do
     if(iget==idtset)then
!     The index of the dataset, from which the data ought to be taken,
!     does not correspond to a previous dataset.
      write(message, '(a,a,a,a,i3,a,a,i3,a,a,a,a)' )ch10,&
&      ' driver : ERROR -',ch10,&
&      '  The component number',idtset,' of the input variable getocc,',&
&        ' equal to',getocc,',',ch10,&
&      '  does not correspond to an existing index.',&
&      '  Action : correct getocc or jdtset in your input file.'
      call wrtout(6,message,'COLL')
      call leave_new('COLL')
     end if
    end if
    occ(:)=results_out(iget)%occ(1:mband*nkpt*nsppol)
    write(message, '(a,i3,a,a)' )&
&    ' driver : getocc/=0, take occ from output of dataset with index',&
&         dtsets(iget)%jdtset,'.',ch10
    call wrtout(ab_out,message,'COLL')
    call wrtout(6,message,'COLL')
   end if

!  Getcell has to be done BEFORE getxcart
!  since acell and rprim will be used
   if(getcell>0 .or. (getcell<0 .and. idtset+getcell>0) )then
!   In case getocc is a negative number (so must add to idtset)
    if(getcell<0 .and. idtset+getcell>0) iget=idtset+getcell
    if(getcell>0)then
     do iget=1,idtset
      if( dtsets(iget)%jdtset==getcell )exit
     end do
     if(iget==idtset)then
!     The index of the dataset, from which the data ought to be taken,
!     does not correspond to a previous dataset.
      write(message, '(a,a,a,a,i3,a,a,i3,a,a,a,a)' )ch10,&
&      ' driver : ERROR -',ch10,&
&      '  The component number',idtset,' of the input variable getcell,',&
&        ' equal to',getcell,',',ch10,&
&      '  does not correspond to an existing index.',&
&      '  Action : correct getcell or jdtset in your input file.'
      call wrtout(6,message,'COLL')
      call leave_new('COLL')
     end if
    end if
    acell(:)=results_out(iget)%acell(:)
    rprim(:,:)=results_out(iget)%rprim(:,:)
    write(message, '(a,i3,a,a)' )&
&    ' driver : getcell/=0, take acell and rprim from output of dataset with index',&
&         dtsets(iget)%jdtset,'.',ch10
    call wrtout(ab_out,message,'COLL')
    call wrtout(6,message,'COLL')

!   Check that the new acell and rprim are consistent with the input dilatmx
    call mkrdim(acell,rprim,rprimd)
    call chkdilatmx(dtset%dilatmx,rprimd,dtset%rprimd_orig)

   end if

   if(getxred>0 .or. (getxred<0 .and. idtset+getxred>0) )then
!   In case getxred is a negative number (so must add to idtset)
    if(getxred<0 .and. idtset+getxred>0) iget=idtset+getxred
    if(getxred>0)then
     do iget=1,idtset
      if( dtsets(iget)%jdtset==getxred )exit
     end do
     if(iget==idtset)then
!     The index of the dataset, from which the data ought to be taken,
!     does not correspond to a previous dataset.
      write(message, '(a,a,a,a,i3,a,a,i3,a,a,a,a)' )ch10,&
&      ' driver : ERROR -',ch10,&
&      '  The component number',idtset,' of the input variable getxred,',&
&        ' equal to',getxred,',',ch10,&
&      '  does not correspond to an existing index.',&
&      '  Action : correct getxred or jdtset in your input file.'
      call wrtout(6,message,'COLL')
      call leave_new('COLL')
     end if
    end if
    xred(:,:)=results_out(iget)%xred(:,1:natom)
    write(message, '(a,i3,a,a)' )&
&    ' driver : getxred/=0, take xred from output of dataset with index',&
&         dtsets(iget)%jdtset,'.',ch10
    call wrtout(ab_out,message,'COLL')
    call wrtout(6,message,'COLL')
   end if

   if(getxcart>0 .or. (getxcart<0 .and. idtset+getxcart>0) )then
!   In case getxcart is a negative number (so must add to idtset)
    if(getxcart<0 .and. idtset+getxcart>0) iget=idtset+getxcart
    if(getxcart>0)then
     do iget=1,idtset
      if( dtsets(iget)%jdtset==getxcart )exit
     end do
     if(iget==idtset)then
!     The index of the dataset, from which the data ought to be taken,
!     does not correspond to a previous dataset.
      write(message, '(a,a,a,a,i3,a,a,i3,a,a,a,a)' )ch10,&
&      ' driver : ERROR -',ch10,&
&      '  The component number',idtset,' of the input variable getxcart,',&
&        ' equal to',getxcart,',',ch10,&
&      '  does not correspond to an existing index.',&
&      '  Action : correct getxcart or jdtset in your input file.'
      call wrtout(6,message,'COLL')
      call leave_new('COLL')
     end if
    end if
!   Compute xcart of the previous dataset
    allocate(xcart(3,natom),xredget(3,natom))
    rprimdget(:,:)=results_out(iget)%rprimd(:,:)
    xredget (:,:)=results_out(iget)%xred(:,1:natom)
    call xredxcart(natom,1,rprimdget,xcart,xredget)
!   xcart from previous dataset is computed. Now, produce xred for the new dataset.
!   Now new acell and rprim ...
    call mkrdim(acell,rprim,rprimd)
    call xredxcart(natom,-1,rprimd,xcart,xred)
    deallocate(xcart,xredget)
    write(message, '(a,i3,a,a)' )&
&    ' driver : getxcart/=0, take xcart from output of dataset with index',&
&         dtsets(iget)%jdtset,'.',ch10
    call wrtout(ab_out,message,'COLL')
    call wrtout(6,message,'COLL')
   end if

   if(getvel>0 .or. (getvel<0 .and. idtset+getvel>0) )then
!   In case getvel is a negative number (so must add to idtset)
    if(getvel<0 .and. idtset+getvel>0) iget=idtset+getvel
    if(getvel>0)then
     do iget=1,idtset
      if( dtsets(iget)%jdtset==getvel )exit
     end do
     if(iget==idtset)then
!     The index of the dataset, from which the data ought to be taken,
!     does not correspond to a previous dataset.
      write(message, '(a,a,a,a,i3,a,a,i3,a,a,a,a)' )ch10,&
&      ' driver : ERROR -',ch10,&
&      '  The component number',idtset,' of the input variable getvel,',&
&        ' equal to',getvel,',',ch10,&
&      '  does not correspond to an existing index.',&
&      '  Action : correct getvel or jdtset in your input file.'
      call wrtout(6,message,'COLL')
      call leave_new('COLL')
     end if
    end if
    vel(:,:)=results_out(iget)%vel(:,1:natom)
    write(message, '(a,i3,a,a)' )&
&    ' driver : getvel/=0, take vel from output of dataset with index',&
&         dtsets(iget)%jdtset,'.',ch10
    call wrtout(ab_out,message,'COLL')
    call wrtout(6,message,'COLL')
   end if

  end if

!****************************************************************************
! Treat the pseudopotentials : initialize the psps variable

!Note that mpsang is the max of 1+lmax, with minimal value 1 (even for local psps, at present)
! mpsang=max(maxval(pspheads(1:npsp)%lmax)+1,1) ! might not work with HP compiler
! n1xccc=maxval(pspheads(1:npsp)%xccc)
  mpsang=1
  n1xccc=pspheads(1)%xccc
  do ii=1,npsp
   mpsang=max(pspheads(ii)%lmax+1,mpsang)
   n1xccc=max(pspheads(ii)%xccc,n1xccc)
  end do

! Determine the maximum number of projectors, for the set of pseudo atom
  call getdim_nloc(lmnmax,lmnmaxso,lnmax,lnmaxso,dtset%mixalch,npsp,dtset%npspalch,&
&  ntypat,dtset%ntypalch,pspheads)

  psps%mpsang  =mpsang
  psps%mtypalch=mtypalch
  psps%npsp    =npsp
  psps%npspalch=dtset%npspalch
  psps%ntypat   =ntypat
  psps%ntypalch=dtset%ntypalch
  psps%ntyppure=dtset%ntyppure

  psps%n1xccc  =n1xccc

  psps%usepaw  =usepaw
  psps%useylm  =dtset%useylm

! psps%pspso(:)=pspheads(1:npsp)%pspso   !! This should be the correct coding
  psps%pspso(:)=1
  psps%pspso(1:min(npsp,ntypat))=dtset%so_typat(1:min(npsp,ntypat))    ! XG020717 Strange, will not work with alchemy and spin-orbit

  allocate(psps%algalch(psps%ntypalch))
  allocate(psps%mixalch(psps%npspalch,psps%ntypalch))
  psps%algalch(:)=dtset%algalch(:)
  psps%mixalch(:,:)=dtset%mixalch(:,:)

! Set mpspso
! Warning : mpspso might be different for each dataset.
  psps%mpspso=1
!OCL SCALAR
  do itypat=1,ntypat
   if(dtset%so_typat(itypat)/=1)psps%mpspso=2
  end do
! XG020729 : THIS SHOULD NOT BE NEEDED : only so_typat should determine mpspso.
! Still, the dimensioning is not large enough at present
!OCL SCALAR
  do ipsp=1,npsp
   if(pspheads(ipsp)%pspso/=1)psps%mpspso=2
  end do

! Set mpssoang, lmnmax, lnmax
  if(psps%mpspso==1)then
   psps%mpssoang=psps%mpsang
   psps%lmnmax  =lmnmax
   psps%lnmax   =lnmax
  else
   psps%mpssoang=2*psps%mpsang-1
   psps%lmnmax=lmnmaxso
   psps%lnmax=lnmaxso
  end if
  if (psps%useylm==0) then
   psps%lmnmax=psps%lnmax
  end if

! Set dimekb
  if (psps%usepaw==0) then
   psps%dimekb=psps%lnmax
  else
   psps%dimekb=psps%lmnmax*(psps%lmnmax+1)/2
  end if

! The following arrays are often not deallocated before the end of the dtset loop
! and might keep their content from one dataset to the other,
! if the conditions are fulfilled
  if(dimekb_old/=psps%dimekb .or. ntypat_old/=ntypat .or. &
&    usepaw_old/=psps%usepaw) then
   if(idtset/=1)deallocate(psps%ekb)
   allocate(psps%ekb(psps%dimekb,ntypat*(1-psps%usepaw)))
   dimekb_old=psps%dimekb
   usepaw_old=psps%usepaw
  end if
  if(lmnmax_old/=psps%lmnmax .or. ntypat_old/=ntypat)then
   if(idtset/=1)deallocate(psps%indlmn)
   allocate(psps%indlmn(6,psps%lmnmax,ntypat))
   lmnmax_old=psps%lmnmax
  end if
  if(mqgrid_old/=psps%mqgrid .or. lnmax_old/=psps%lnmax .or. ntypat_old/=ntypat)then
   if(idtset/=1)deallocate(psps%ffspl,psps%qgrid,psps%vlspl)
   allocate(psps%ffspl(psps%mqgrid,2,psps%lnmax,ntypat),psps%qgrid(psps%mqgrid))
   allocate(psps%vlspl(psps%mqgrid,2,ntypat))
   mqgrid_old=psps%mqgrid
   lnmax_old=psps%lnmax
  end if
  if(ntypat_old/=ntypat)then
   if(idtset/=1)deallocate(psps%xcccrc,psps%xccc1d,psps%ziontypat,psps%znucltypat)
   allocate(psps%xcccrc(ntypat),psps%xccc1d(n1xccc,6,ntypat))
   allocate(psps%znucltypat(ntypat))
   allocate(psps%ziontypat(ntypat))
   paw_size=max(ntypat,npsp)
   ntypat_old=ntypat
  end if
  psps%ziontypat(:)=dtset%ziontypat(:)

!****************************************************************************
! PAW allocations.
  if (psps%usepaw==1) then

!  The correct dimension of pawrad/tab is ntypat.
!  In case of paw, no alchemical psp is allowed, so npsp=ntypat
!  However, in case of alchemical psps, pawrad/tab(ipsp) is invoked in
!  pspini. So, in order to avoid any problem, declare pawrad/tab
!  at paw_size=max(ntypat,npsp).
   if (paw_size/=paw_size_old) then
    if(idtset/=1) then
     do itypat=1,ntypat
      deallocate(pawrad(itypat)%rad)
      deallocate(pawrad(itypat)%radfact)
      deallocate(pawtab(itypat)%gnorm)
      deallocate(pawtab(itypat)%shapefunc)
      deallocate(pawtab(itypat)%tphi)
      deallocate(pawtab(itypat)%phi)
      deallocate(pawtab(itypat)%tphitphj)
      deallocate(pawtab(itypat)%phiphj)
      deallocate(pawtab(itypat)%coredens)
      deallocate(pawtab(itypat)%tcoredens)
      deallocate(pawtab(itypat)%qijl)
      deallocate(pawtab(itypat)%eijkl)
      deallocate(pawtab(itypat)%dij0)
      deallocate(pawtab(itypat)%dltij)
      deallocate(pawtab(itypat)%rhoij0)
      deallocate(pawtab(itypat)%sij)
      deallocate(pawrad,pawtab)
     end do
    end if
    allocate(pawrad(paw_size),pawtab(paw_size))
   end if

   l_max_new=mpsang
   l_size_max_new=2*l_max_new-1
   angl_size_new=dtset%pawntheta*dtset%pawnphi

   do itypat=1,ntypat

    basis_size_new=pspheads(itypat)%pawheader%basis_size
    lmn_size_new  =pspheads(itypat)%pawheader%lmn_size
    l_size_new    =pspheads(itypat)%pawheader%l_size
    mesh_size_new =pspheads(itypat)%pawheader%mesh_size
    lmn2_size_new =lmn_size_new*(lmn_size_new+1)/2
    l_size_max_new=max(l_size_new,l_size_max_new)
    ij_size_new   =basis_size_new*(basis_size_new+1)/2

!   Reallocate arrays depending on mesh_size and basis_size
    test_new=(paw_size/=paw_size_old)
    if(.not.test_new)test_new=((mesh_size_new/=pawtab(itypat)%mesh_size)&
&                          .or.(basis_size_new/=pawtab(itypat)%basis_size))
    if (test_new) then
     if((idtset/=1).and.(paw_size==paw_size_old))&
&      deallocate(pawtab(itypat)%tphi,pawtab(itypat)%phi)
     allocate(pawtab(itypat)%tphi(mesh_size_new,basis_size_new))
     allocate(pawtab(itypat)%phi (mesh_size_new,basis_size_new))
     pawtab(itypat)%basis_size=basis_size_new
    end if

!   Reallocate arrays depending on mesh_size and ij_size
    test_new=(paw_size/=paw_size_old)
    if(.not.test_new)test_new=((mesh_size_new/=pawtab(itypat)%mesh_size)&
&                          .or.(ij_size_new/=pawtab(itypat)%ij_size))
    if (test_new) then
     if((idtset/=1).and.(paw_size==paw_size_old))&
&      deallocate(pawtab(itypat)%tphitphj,pawtab(itypat)%phiphj)
     allocate(pawtab(itypat)%tphitphj(mesh_size_new,ij_size_new))
     allocate(pawtab(itypat)%phiphj  (mesh_size_new,ij_size_new))
     pawtab(itypat)%ij_size   =ij_size_new
    end if

!   Reallocate arrays depending on mesh_size and l_size
    test_new=(paw_size/=paw_size_old)
    if(.not.test_new)test_new=((mesh_size_new/=pawtab(itypat)%mesh_size)&
&                         .or.(l_size_new/=pawtab(itypat)%l_size))
    if (test_new) then
     if((idtset/=1).and.(paw_size==paw_size_old))&
&      deallocate(pawtab(itypat)%shapefunc)
     allocate(pawtab(itypat)%shapefunc(mesh_size_new,l_size_new))
    end if

!   Reallocate arrays depending on l_size and lmn2_size
    test_new=(paw_size/=paw_size_old)
    if(.not.test_new)test_new=((l_size_new/=pawtab(itypat)%l_size)&
&                          .or.(lmn2_size_new/=pawtab(itypat)%lmn2_size))
    if (test_new) then
     if((idtset/=1).and.(paw_size==paw_size_old))&
&      deallocate(pawtab(itypat)%qijl)
     allocate(pawtab(itypat)%qijl(l_size_new*l_size_new,lmn2_size_new))
    end if

!   Reallocate arrays depending on l_size
    test_new=(paw_size/=paw_size_old)
    if(.not.test_new)test_new=(l_size_new/=pawtab(itypat)%l_size)
    if (test_new) then
     if((idtset/=1).and.(paw_size==paw_size_old))&
&      deallocate(pawtab(itypat)%gnorm)
     allocate(pawtab(itypat)%gnorm(l_size_new))
     pawtab(itypat)%l_size=l_size_new
    end if

!   Reallocate arrays depending on lmn2_size
    test_new=(paw_size/=paw_size_old)
    if(.not.test_new)test_new=(lmn2_size_new/=pawtab(itypat)%lmn2_size)
    if (test_new) then
     if((idtset/=1).and.(paw_size==paw_size_old))&
&      deallocate(pawtab(itypat)%eijkl,pawtab(itypat)%dij0,&
&                 pawtab(itypat)%dltij,pawtab(itypat)%rhoij0,&
&                 pawtab(itypat)%sij)
     allocate(pawtab(itypat)%eijkl(lmn2_size_new,lmn2_size_new))
     allocate(pawtab(itypat)%dij0(lmn2_size_new))
     allocate(pawtab(itypat)%dltij(lmn2_size_new))
     allocate(pawtab(itypat)%rhoij0(lmn2_size_new))
     allocate(pawtab(itypat)%sij(lmn2_size_new))
     pawtab(itypat)%lmn_size =lmn_size_new
     pawtab(itypat)%lmn2_size=lmn2_size_new
    end if

!   Reallocate arrays depending on mesh_size
    test_new=(paw_size/=paw_size_old)
    if(.not.test_new)test_new=(mesh_size_new/=pawtab(itypat)%mesh_size)
    if (test_new) then
     if((idtset/=1).and.(paw_size==paw_size_old))&
&      deallocate(pawtab(itypat)%coredens,pawtab(itypat)%tcoredens,&
&                 pawrad(itypat)%rad,pawrad(itypat)%radfact)
     allocate(pawtab(itypat)%coredens (mesh_size_new))
     allocate(pawtab(itypat)%tcoredens(mesh_size_new))
     allocate(pawrad(itypat)%rad      (mesh_size_new))
     allocate(pawrad(itypat)%radfact  (mesh_size_new))
     pawtab(itypat)%mesh_size=mesh_size_new
    end if

    paw_size_old=paw_size
   end do ! itypat

   if (idtset==1) then
    angl_size=-1;l_size_max=-1;l_max=-1
   end if

!  Reallocate arrays depending on angl_size and l_size_max
   if ((angl_size_new/=angl_size).or.&
&      (l_size_max_new/=l_size_max)) then
    if(idtset/=1) deallocate(pawang%ylmr)
    allocate(pawang%ylmr(l_size_max_new**2,angl_size_new))
   end if

!  Reallocate arrays depending on l_max and l_size_max
   if ((l_max_new/=l_max).or.&
&      (l_size_max_new/=l_size_max)) then
    if(idtset/=1) deallocate(pawang%zarot)
    allocate(pawang%zarot(l_size_max_new,l_size_max_new,l_max_new,nsym))
    l_size_max=l_size_max_new
   end if

!  Reallocate arrays depending on l_max
   if (l_max_new/=l_max) then
    if(idtset/=1) deallocate(pawang%gntselect)
    allocate(pawang%gntselect((2*l_max_new-1)**2,l_max_new**2,l_max_new**2))
    l_max=l_max_new
   end if

!  Reallocate arrays depending on angl_size
   if (angl_size_new/=angl_size) then
    if(idtset/=1) deallocate(pawang%anginit,pawang%angwgth)
    allocate(pawang%anginit(3,angl_size_new))
    allocate(pawang%angwgth(angl_size_new))
    angl_size=angl_size_new
   end if

  end if

!****************************************************************************
! At this stage, all the data needed for the treatment of one dataset
! have been transferred from multi-dataset arrays.

  mkmem=dtset%mkmem
  iexit=0

! Smaller integer arrays :
  allocate(npwtot(nkpt))

  allocate(results_gs%fcart(3,natom),results_gs%fred(3,natom))
! Also allocate some other components of results_gs,
! even if they are not used in the present routine.
  allocate(results_gs%gresid(3,natom),results_gs%grewtn(3,natom))
  allocate(results_gs%grxc(3,natom),results_gs%synlgr(3,natom))

! Initialize some of these to zero (needed when hasty exit)
  etotal=zero
  strten(:)=zero
  results_gs%eei       =zero
  results_gs%eeig      =zero
  results_gs%eew       =zero
  results_gs%ehart     =zero
  results_gs%eii       =zero
  results_gs%ek        =zero
  results_gs%enl       =zero
  results_gs%entropy   =zero
  results_gs%enxc      =zero
  results_gs%etotal    =zero
  results_gs%fcart(:,:)=zero
  results_gs%fred(:,:) =zero
  results_gs%strten(:) =zero

  if(optdriver==0)then

   call status(jdtset_status,filstat,iexit,level,'call gstate   ')

   if(dtset%parareel==0)then

    mpi_enreg%parareel=0
    call gstate(acell,codvsn,cpui,dtfil,dtset,iexit,&
&    mband,mgfft,mkmem,mpi_enreg,mpw,natom,nfft,nkpt,&
&    npwtot,nspden,nspinor,nsppol,nsym,&
&    occ,pawang,pawrad,pawtab,psps,results_gs,rprim,vel,walli,xred)

   else

!#if defined PARAREEL
    mpi_enreg%parareel=1
    mpi_enreg%paral_level=1
    call pstate(acell,codvsn,cpui,dtfil,dtset,iexit,&
&    mband,mgfft,mkmem,mpi_enreg,mpw,natom,nfft,nkpt,&
&    npwtot,nspden,nspinor,nsppol,nsym,&
&    occ,pawrad,pawtab,psps,results_gs,rprim,vel,walli,xred)
!#endif

   end if

   etotal=results_gs%etotal
   strten(1:6)=results_gs%strten(1:6)

   call status(jdtset_status,filstat,iexit,level,'after gstate  ')

  elseif(optdriver==1)then

   call status(jdtset_status,filstat,iexit,level,'call respfn   ')

   mpi_enreg%paral_level=2
   call respfn(codvsn,cpui,dtfil,dtset,etotal,iexit,&
&   mband,mgfft,mkmem,mkmems,mpi_enreg,mpw,natom,nfft,nkpt,npwtot,nspden,&
&   nspinor,nsppol,nsym,occ,pawrad,pawtab,psps,walli,xred)

   call status(jdtset_status,filstat,iexit,level,'after respfn  ')

  elseif(optdriver==2)then

   call status(jdtset_status,filstat,iexit,level,'call suscep   ')

   mpi_enreg%paral_level=2
   call suscep(dtfil,dtset,iexit,&
&   mband,mgfft,mkmem,mpi_enreg,mpw,natom,nfft,nkpt,&
&   nspden,nspinor,nsppol,nsym,occ,xred)

   call status(jdtset_status,filstat,iexit,level,'after suscep  ')

  elseif(optdriver==3)then

   call status(jdtset_status,filstat,iexit,level,'call screening')

   call screening(acell,dtfil,dtset,rprim)

   call status(jdtset_status,filstat,iexit,level,'after screenin')

  elseif(optdriver==4)then

   call status(jdtset_status,filstat,iexit,level,'call sigma   ')

   mpi_enreg%paral_level=0
   call sigma(acell,dtfil,dtset,mpi_enreg,rprim)

   call status(jdtset_status,filstat,iexit,level,'after sigma  ')

  elseif(optdriver==5)then

   call status(jdtset_status,filstat,iexit,level,'call nonlinear   ')

   mpi_enreg%paral_level=2
   call nonlinear(codvsn,dtfil,dtset,etotal,iexit,&
&   mband,mgfft,mkmem,mpi_enreg,mpw,natom,nfft,nkpt,npwtot,nspden,&
&   nspinor,nsppol,nsym,occ,pawrad,pawtab,psps,xred)

   call status(jdtset_status,filstat,iexit,level,'after nonlinear  ')

  else

!  Error the choice is either 0 -> gstate, 1 -> respfn, 2 -> suscep,
!                             3 -> screening, 4 -> sigma,  5 -> nonlinear
   write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
&   ' driver : ERROR -',ch10,&
&   '  The variable optdriver must be between 0 and 6, but was ',&
&   optdriver,ch10,&
&   '  This is not allowed.  ',ch10,&
&   '  Action : modify optdriver in the input file.'
   call wrtout(06,  message,'COLL')
   call leave_new('COLL')

  end if

!****************************************************************************

!DEBUG
!write(6,*)' driver : after respfn'
!ENDDEBUG

! Transfer of multi dataset outputs from temporaries :
! acell, xred, occ rprim, and vel might be modified from their
! input values
! etotal, fcart, fred, and strten have been computed
! npwtot was already computed before, but is stored only now

  results_out(idtset)%acell(:)          =acell(:)
  results_out(idtset)%etotal            =etotal
  results_out(idtset)%rprim(:,:)        =rprim(:,:)
  call mkrdim(acell,rprim,rprimd)
  results_out(idtset)%rprimd(:,:)       =rprimd(:,:)
  results_out(idtset)%strten(:)         =strten(:)
  results_out(idtset)%fcart(1:3,1:natom)=results_gs%fcart(:,:)
  results_out(idtset)%fred(1:3,1:natom) =results_gs%fred(:,:)
  results_out(idtset)%npwtot(1:nkpt)       =npwtot(1:nkpt)
  results_out(idtset)%occ(1:mband*nkpt*nsppol)=occ(:)
  results_out(idtset)%vel(:,1:natom)    =vel(:,:)
  results_out(idtset)%xred(:,1:natom)   =xred(:,:)

  deallocate(dtset%amu,dtset%algalch,dtset%bdgw)
  deallocate(dtset%densty,dtset%iatfix,dtset%iatsph,dtset%istwfk)
  deallocate(dtset%kberry,dtset%kptgw,dtset%kpt,dtset%kptns)
  deallocate(dtset%mixalch,dtset%nband)
  deallocate(dtset%occ_orig,dtset%so_typat,dtset%shiftk,dtset%spinat)
  deallocate(dtset%symafm,dtset%symrel,dtset%tnons,dtset%typat)
  deallocate(dtset%vel_orig,dtset%wtatcon,dtset%wtk)
  deallocate(dtset%xred_orig,dtset%ziontypat,dtset%znucl)

  deallocate(psps%algalch,psps%mixalch)

  deallocate(occ)
  deallocate(vel,xred)
  deallocate(npwtot)
  deallocate(results_gs%fcart,results_gs%fred)
  deallocate(results_gs%gresid,results_gs%grewtn)
  deallocate(results_gs%grxc,results_gs%synlgr)

  if(iexit/=0)exit

! Check whether exiting was required by the user.
! If found then beat a hasty exit from time steps
  openexit=1 ; if(dtset%chkexit==0) openexit=0
! update for the end of parareel case : chkexi will crash otherwise
  if (mpi_enreg%parareel == 1) then
        mpi_enreg%parareel=0
  end if
  call chkexi(zero,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg,openexit)
  if (iexit/=0)exit

! End do loop on idtset (allocate statements are present -
! an exit statement is present)
 end do

!*********************************************************************

 deallocate(psps%qgrid)
 deallocate(psps%xcccrc)
 deallocate(psps%xccc1d)
 deallocate(psps%vlspl)
 deallocate(psps%ekb)
 deallocate(psps%indlmn)
 deallocate(psps%filpsp)
 deallocate(psps%pspcod)
 deallocate(psps%pspdat)
 deallocate(psps%pspso)
 deallocate(psps%pspxc)
 deallocate(psps%title)
 deallocate(psps%znuclpsp)
 deallocate(psps%znucltypat)
 deallocate(psps%zionpsp)
 deallocate(psps%ziontypat)
 deallocate(psps%ffspl)

 if (psps%usepaw==1) then
  do itypat=1,ntypat
   deallocate(pawrad(itypat)%rad)
   deallocate(pawrad(itypat)%radfact)
   deallocate(pawtab(itypat)%gnorm)
   deallocate(pawtab(itypat)%shapefunc)
   deallocate(pawtab(itypat)%tphi)
   deallocate(pawtab(itypat)%phi)
   deallocate(pawtab(itypat)%tphitphj)
   deallocate(pawtab(itypat)%phiphj)
   deallocate(pawtab(itypat)%coredens)
   deallocate(pawtab(itypat)%tcoredens)
   deallocate(pawtab(itypat)%sij)
   deallocate(pawtab(itypat)%eijkl)
   deallocate(pawtab(itypat)%dij0)
   deallocate(pawtab(itypat)%dltij)
   deallocate(pawtab(itypat)%qijl)
  end do
  deallocate(pawang%anginit,pawang%angwgth)
  deallocate(pawang%ylmr,pawang%zarot)
  deallocate(pawang%gntselect,pawang%realgnt)
  deallocate(pawrad,pawtab)
 end if

 deallocate(jdtset_)

!DEBUG
! write(6,*)' driver : before exit '
! write(6, '(a,9f5.2)' ) ' driver : before exit  , rprim ',rprim (:,:)
!stop
!ENDDEBUG

 call status(0,filstat,iexit,level,'exit          ')
 call timab(100,2,tsec)

 end subroutine driver
!!***

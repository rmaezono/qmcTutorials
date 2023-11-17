!{\src2tex{textfont=tt}}
!!****f* ABINIT/gstate
!! NAME
!! gstate
!!
!! FUNCTION
!! Primary routine for conducting DFT calculations by CG minimization.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2004 ABINIT group (DCA, XG, GMR, JYR, MKV, MT, FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors.
!!
!! INPUTS
!!  codvsn=code version
!!  cpui=initial CPU time
!!  natom =number of atoms in unit cell
!!  mband =maximum number of bands
!!  nfft  =(effective) number of FFT grid points (for this processor)
!!  mgfft =maximum single fft dimension
!!  mkmem =maximum number of k points which can fit in core memory
!!  nkpt  =number of k points
!!  mpw   =maximum number of planewaves in basis sphere (large number)
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=number of channels for spin-polarization (1 or 2)
!!  nsym  =number of symmetry elements in space group
!!  walli=initial wall clock time
!!
!! OUTPUT
!!  npwtot(nkpt) = total number of plane waves at each k point
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation
!!
!! SIDE EFFECTS
!!  acell(3)=unit cell length scales (bohr)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  iexit= exit flag
!!  mpi_enreg=MPI-parallelisation information (some already initialized,
!!   some others to be initialized here)
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   Before entering the first time in gstate, a significant part of
!!   psps has been initialized :
!!   the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,mpsso,mgrid,
!!     ntypat,n1xccc,usepaw,useylm, and the arrays dimensioned to npsp
!!   All the remaining components of psps are to be initialized in the call
!!   to pspini .
!!   The next time the code enters gstate, psps might be identical to the
!!   one of the previous dtset, in which case, no reinitialisation is scheduled
!!   in pspini.f .
!!  rprim(3,3)=dimensionless real space primitive translations
!!  vel(3,natom)=value of velocity
!!  xred(3,natom) = reduced atomic coordinates
!!
!! TODO
!!
!! Not yet possible to use restartxf in parallel when localrdwf==0
!!
!! PARENTS
!!      driver,pstate
!!
!! CHILDREN
!!      blok8,brdmin,bstruct_clean,bstruct_init,chkexi,clnmpi_fft,clnmpi_gs
!!      clnup1,clnup2,delocint,fconv,fixsym,fourdp,getgsc,getph,hdr_clean
!!      hdr_init,hdr_update,initberry,initmpi_fft,initmpi_gs,initro,initylmg
!!      int2char4,inwffil,ioarr,ioddb8,kpgio,leave_new,mkrho,moldyn,move,newocc
!!      outwf,outxfhist,pawinit,prtene,psddb8,pspini,scfcv,setsym,setsymrhoij
!!      setup1,setup2,status,timab,wffclose,wffdelete,wffopen,wffreadskiprec
!!      wrtout,xcomm_world,xme_whoiam,xproc_max
!!
!! SOURCE

 subroutine gstate(acell,codvsn,cpui,dtfil,dtset,iexit,&
& mband,mgfft,mkmem,mpi_enreg,mpw,natom,nfft,nkpt,&
& npwtot,nspden,nspinor,nsppol,nsym,&
& occ,pawang,pawrad,pawtab,psps,results_gs,rprim,vel,walli,xred)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: mband,mgfft,mkmem,mpw,natom,nfft,nkpt,nspden,nspinor
 integer,intent(in) :: nsppol,nsym
 integer,intent(inout) :: iexit
 real(dp),intent(in) :: cpui,walli
 character(len=6),intent(in) :: codvsn
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(inout) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(pawang_type),intent(inout) :: pawang
 type(pseudopotential_type),intent(inout) :: psps
 type(results_gs_type),intent(out) :: results_gs
 integer,intent(out) :: npwtot(nkpt)
 real(dp),intent(inout) :: acell(3),occ(mband*nkpt*nsppol),rprim(3,3)
 real(dp),intent(inout) :: vel(3,natom),xred(3,natom)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!Define file format for different type of files. Presently,
!only one file format is supported for each type of files, but this might
!change soon ...
!2   for wavefunction file, new format (version 2.0 and after)    (fform)
!52  for density rho(r)       (fformr)
!102 for potential V(r) file. (fformv)
 integer,parameter :: fform=2,fformr=52,fformv=102,formeig=0,level=3,response=0
 integer :: accessfil,ask_accurate,bantot,blktyp,choice,fullinit,gscase,iapp
 integer :: iatom,idir,ierr,ifft,ii,index,initialized,ionmov,ios,ir,iscf,isppol
 integer :: itime,itimexit,itypat,ixfh,master,me,mpert,mpsang,msize,mu,mxfh
 integer :: mygroup,nblok,normchoice,nproc,ntime,ntypat,nxfh,openexit,option
 integer :: prtvol,psp_gencond,pwind_alloc,rdwr,restartxf,spaceworld,tim_mkrho
 integer :: tmkmem,toto,vrsddb
 real(dp) :: cpus,ecore,ecut_eff,epulay,etot,fermie,gsqcut_eff,gsqcut_eff_
 real(dp) :: gsqcutbox_eff,residm,rhosum,tolwfr,ucvol
 logical :: ex,od
 character(len=3) :: ipara
 character(len=4) :: tag
 character(len=500) :: message
 character(len=fnlen) :: ddbnm,dscrpt
 type(bandstructure_type) :: bstruct
 type(dens_sym_operator_type) :: densymop_gs
 type(efield_type) :: dtefield
 type(hdr_type) :: hdr
 type(pawene_type) :: pawene
 type(wffile_type) :: wff1,wffnew,wffnow
 integer,allocatable :: atindx(:),atindx1(:),blkflg(:),indsym(:,:,:)
 integer,allocatable :: irrzon(:,:,:),kg(:,:),nattyp(:),npwarr(:),symrec(:,:,:)
 integer,pointer :: pwind(:,:,:)
 real(dp) :: blknrm(3),blkqpt(9),corstr(6),gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp) :: rprimd(3,3),tsec(2)
 real(dp),allocatable :: amass(:),blkval(:,:),cg(:,:),doccde(:),dyfrx2(:,:,:)
 real(dp),allocatable :: eigen(:),gsc(:,:),ph1d(:,:),phnons(:,:,:),resid(:)
 real(dp),allocatable :: rhog(:,:),rhor(:,:),start(:,:),work(:),xfhist(:,:,:,:)
 real(dp),allocatable :: xred_old(:,:),ylm(:,:)
 real(dp),pointer :: pwnsfac(:,:)
 character(len=fnlen) :: tmpfil(8)
!no_abirules
 interface
   subroutine initberry(dtefield,dtfil,dtset,gmet,kg,mband,mkmem,mpi_enreg,&
&              mpw,nkpt,npwarr,nsppol,nsym,occ,pwind,pwind_alloc,pwnsfac,rprimd,symrec)
     use defs_basis
     use defs_datatypes
     use defs_xfuncmpi
     integer,intent(in) :: mband,mkmem,mpw,nkpt,nsppol,nsym
     integer,intent(out):: pwind_alloc
     integer,intent(in) :: kg(3,mpw*mkmem),npwarr(nkpt),symrec(3,3,nsym)
     integer,pointer :: pwind(:,:,:)
     real(dp),pointer :: pwnsfac(:,:)
     real(dp),intent(in) :: gmet(3,3),occ(mband*nkpt*nsppol),rprimd(3,3)
     type(MPI_type),intent(inout) :: mpi_enreg
     type(datafiles_type),intent(in) :: dtfil
     type(dataset_type),intent(inout) :: dtset
     type(efield_type),intent(out) :: dtefield
   end subroutine initberry
 end interface

! ***********************************************************************

!DEBUG
!write(6,*)' gstate : enter'
!stop
!ENDDEBUG

 call timab(32,1,tsec)
 call timab(33,1,tsec)

 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Set up mpi informations from the dataset
 if (mpi_enreg%parareel == 0) then
  mpi_enreg%paral_level=2
  call initmpi_gs(dtset,mpi_enreg)
 else
  mpi_enreg%paral_level=1
 end if

 call initmpi_fft(dtset,mpi_enreg)


!Init spaceworld
 call xcomm_world(spaceworld)
 master =0
!Define me
 call xme_whoiam(me)
!Define nproc
 call xproc_max(nproc,ierr)

!Structured debugging if prtvol==-level
 prtvol=dtset%prtvol
 if(prtvol==-level)then
  write(message,'(80a,a,a)')  ('=',ii=1,80),ch10,&
&   ' gstate : enter , debug mode '
  call wrtout(06,message,'COLL')
 end if

 ntime=dtset%ntime

!Option input variables
 ionmov   =dtset%ionmov
 iscf     =dtset%iscf
 restartxf=dtset%restartxf

!Create names for the temporary files based on dtfil%filnam_ds(5)
!by appending adequate string.
!'_WF1' -> dtfil%unwft1
!'_WF2' -> dtfil%unwft2
!'_KG' ->  dtfil%unkg
!'_DUM' -> tmp_unit (real dummy name)
!'_YLM' -> dtfil%unylm
!'_GSC' -> dtfil%ungsc
!'_PAW' -> dtfil%unpaw
 tmpfil(1)=trim(dtfil%filnam_ds(5))//'_WF1'
 tmpfil(2)=trim(dtfil%filnam_ds(5))//'_WF2'
 tmpfil(3)=trim(dtfil%filnam_ds(5))//'_KG'
 tmpfil(4)=trim(dtfil%filnam_ds(5))//'_DUM'
 tmpfil(6)=trim(dtfil%filnam_ds(5))//'_YLM'
 tmpfil(7)=trim(dtfil%filnam_ds(5))//'_GSC'
 tmpfil(8)=trim(dtfil%filnam_ds(5))//'_PAW'

 if(mpi_enreg%paral_compil_kpt==1)then
! This is the parallel case : the index of the processor must be appended
  call int2char4(mpi_enreg%me,tag)
  do ii=1,8
   tmpfil(ii)=trim(tmpfil(ii))//'_P-'//tag
  end do
 end if

 call status(0,dtfil%filstat,iexit,level,'call setup1   ')

 initialized=0
 ecore=zero ; epulay=zero

 results_gs%grewtn(:,:)=zero
 results_gs%eei        =zero
 results_gs%eeig       =zero
 results_gs%eew        =zero
 results_gs%ehart      =zero
 results_gs%eii        =zero
 results_gs%ek         =zero
 results_gs%enefield   =zero
 results_gs%enl        =zero
 results_gs%enxc       =zero
 results_gs%pel(1:3)   =zero

!Some initializations for PAW
 if (psps%usepaw==1) then
  pawene%etot=zero  ;pawene%etotdc=zero
  pawene%e1=zero    ;pawene%e1dc=zero
  pawene%etild=zero ;pawene%etilddc=zero
  pawene%etild1=zero;pawene%etild1dc=zero
 end if

!Set up for iterations
 ecut_eff= dtset%ecut * (dtset%dilatmx)**2
 allocate(amass(natom))
 call setup1(acell,amass,dtset%amu,bantot,&
& ecut_eff,gmet,gprimd,gsqcut_eff,gsqcutbox_eff,dtset%intxc,ionmov,&
& nfft,natom,dtset%nband,dtset%ngfft,nkpt,dtset%nqpt,nsppol,nsym,psps%ntypat,&
& dtset%qptn,response,rmet,rprim,rprimd,dtset%typat,ucvol)
 call status(0,dtfil%filstat,iexit,level,'call kpgio    ')

!Set up the basis sphere of planewaves
 allocate(kg(3,mpw*mkmem),npwarr(nkpt))
 call kpgio(ecut_eff,dtset%exchn2n3,gmet,dtset%istwfk,kg,tmpfil(3),dtset%kptns,mkmem,&
& dtset%nband,nkpt,'PERS',mpi_enreg,mpw,npwarr,npwtot,nsppol,dtfil%unkg)

!Set up the Ylm for each k point
 mpsang=psps%mpsang
 allocate(ylm(mpw*mkmem,mpsang*mpsang*psps%useylm))
 if (psps%useylm==1) then
  call status(0,dtfil%filstat,iexit,level,'call initylmg ')
  call initylmg(gprimd,kg,dtset%kptns,mkmem,mpi_enreg,mpsang,mpw,nkpt,&
&               npwarr,dtfil%unkg,dtfil%unylm,ylm,tmpfil(6))
 end if

 call timab(33,2,tsec)

!Open and read pseudopotential files
 call status(0,dtfil%filstat,iexit,level,'call pspini   ')
 gsqcut_eff_=gsqcut_eff;if (dtset%pawvlbox==1) gsqcut_eff_=gsqcutbox_eff
 call pspini(ecore,psp_gencond,gsqcut_eff_,iscf,dtset%ixc,level,&
&            natom,pawrad,pawtab,prtvol,psps,dtset%typat)

 call timab(33,1,tsec)

!Initialize band structure datatype
 allocate(doccde(bantot),eigen(bantot))
 doccde(:)=zero ; eigen(:)=zero
 call bstruct_init(bantot,bstruct,doccde,eigen,dtset%istwfk,dtset%kptns,&
& dtset%nband,nkpt,npwarr,nsppol,occ,dtset%wtk)
 deallocate(doccde,eigen)

!Initialize header
 gscase=0
 call hdr_init(bstruct,codvsn,dtset,hdr,pawtab,gscase,psps)

!Update header, with evolving variables, when available
!Here, rprimd, xred and occ are available
 etot=hdr%etot ; fermie=hdr%fermie ; residm=hdr%residm
 call hdr_update(bantot,etot,fermie,hdr,natom,&
&                residm,rprimd,occ,xred)

!Clean band structure datatype (should use it more in the future !)
 call bstruct_clean(bstruct)

 call status(0,dtfil%filstat,iexit,level,'call inwffil  ')

 allocate(cg(2,mpw*nspinor*mband*mkmem*nsppol))
 allocate(eigen(mband*nkpt*nsppol))
 allocate(resid(mband*nkpt*nsppol))
 eigen(:)=0.0d0 ; resid(:)=0.0d0
! mpi_enreg%paralbd=0 ; ask_accurate=0
 ask_accurate=0

!XG 020711 : dtfil should not be reinitialized here !!!
 if (mpi_enreg%parareel == 1) then
  if (mpi_enreg%ipara > 0 ) then
   if (mpi_enreg%jpara == 0) then
    dtfil%ireadwf = 0
   else
    dtfil%ireadwf = 0
    if (mpi_enreg%ipara < 11) write(ipara,'(i1)')mpi_enreg%ipara-1
    if (mpi_enreg%ipara >= 11) write(ipara,'(i2)')mpi_enreg%ipara-1
    if (mpi_enreg%ipara >= 101) write(ipara,'(i3)')mpi_enreg%ipara-1
    dtfil%fnamewffk=trim(dtfil%filnam_ds(4))//'_WFK_'//ipara
   end if
  else
   dtfil%ireadwf = 0
  end if
 end if

!Initialize wavefunctions.
!Warning : ideally, results_gs%fermie and results_gs%residm
!should not be initialized here. One might make them separate variables.

 wff1%unwff=dtfil%unwff1
 call inwffil(ask_accurate,cg,dtset,dtset%ecut,ecut_eff,eigen,dtset%exchn2n3,&
& formeig,gmet,hdr,dtfil%ireadwf,dtset%istwfk,kg,dtset%kptns,&
& dtset%localrdwf,mband,mkmem,mpi_enreg,&
& mpw,dtset%nband,dtset%ngfft,nkpt,npwarr,nspden,nspinor,nsppol,nsym,occ,&
& prtvol,rprimd,dtset%symafm,dtset%symrel,dtset%tnons,&
& dtfil%unkg,wff1,wffnow,dtfil%unwff1,dtfil%unwft1,&
& dtfil%fnamewffk,tmpfil(1))

!Initialize xf history (should be put in inwffil)
 nxfh=0
 if(restartxf>=1 .and. dtfil%ireadwf==1)then

! Should exchange the data about history in parallel localrdwf==0
  if(mpi_enreg%paral_compil_kpt==1 .and. dtset%localrdwf==0)then
   write(message, '(a,a,a,a,a,a)' )ch10,&
&   ' gstate : BUG -',ch10,&
&   '  It is not yet possible to use non-zero restartxf,',ch10,&
&   '  in parallel, when localrdwf=0. Sorry for this ...'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

  allocate(xfhist(3,natom+4,2,0))
  call outxfhist(nxfh,natom,mxfh,xfhist,2,wff1,ios)
  deallocate(xfhist)

  if(ios>0)then
   write(message, '(a,a,a,a,a,a)' )ch10,&
&   ' gstate : BUG -',ch10,&
&   '  An error occurred reading the input wavefunction file,',ch10,&
&   '  with restartxf=1.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  else if(ios==0)then
   write(message, '(a,a,i4,a)' )ch10,&
&   ' gstate : reading',nxfh,' (x,f) history pairs from input wf file.'
   call wrtout(6,message,'COLL')
   call wrtout(ab_out,message,'COLL')
  end if
!WARNING : should check that restartxf is not negative
!WARNING : should check that restartxf /= only when dtfil%ireadwf is activated
 end if

!Allocate the xf history array : takes into account the existing
!pairs, minus those that will be discarded, then those that will
!be computed, governed by ntime, and some additional pairs
!(needed when it will be possible to use xfhist for move.f)
 mxfh=(nxfh-restartxf+1)+ntime+5
 if(mpi_enreg%parareel==1)mxfh=mxfh+500  ! XG020711 : why this value ?
 allocate(xfhist(3,natom+4,2,mxfh))
!WARNING : should check that the number of atoms in the wf file and natom
!are the same

!Initialize the xf history array
 if(nxfh>=restartxf .and. nxfh>0)then
! Eventually skip some of the previous history
  if(restartxf>=2)then
   do ixfh=1,restartxf-1
    call WffReadSkipRec(wff1,ios)
   end do
  end if

! Read and store the relevant history
  call outxfhist(nxfh-restartxf+1,natom,mxfh,xfhist,3,wff1,ios)
 end if

!Close wff1, if it was ever opened (in inwffil)
 if (dtfil%ireadwf==1) then
  call WffClose(wff1,ierr)
 end if

!Initialize second wavefunction file if needed
 if(mkmem==0 .and. dtset%nstep/=0) then
  write(message, '(a,i4,a,a)' )&
&  ' gstate about to open unit',dtfil%unwft2,' for file=',trim(tmpfil(2))
  call wrtout(06,message,'PERS')
  call WffOpen(dtset%accesswff,spaceworld,tmpfil(2),ierr,wffnew,master,me,dtfil%unwft2)
 end if

 call status(0,dtfil%filstat,iexit,level,'call setup2   ')

!Further setup
 allocate(start(3,natom))
 call setup2(dtset%dedlnn,dtset%ecut,epulay,iscf,dtset%istwfk,natom,nkpt,npwtot,&
& start,ucvol,dtset%wtk,xred)

!Allocation for forces and atomic positions
 allocate(xred_old(3,natom))

!Do symmetry stuff only for nsym>1
 allocate(irrzon(nfft**(1-1/nsym),2,nspden/nsppol))
 allocate(phnons(2,nfft**(1-1/nsym),nspden/nsppol))
 irrzon(:,:,:)=0
 allocate(indsym(4,nsym,natom),symrec(3,3,nsym))

 if (nsym>1) then

  call status(0,dtfil%filstat,iexit,level,'call setsym   ')
  call setsym(densymop_gs,indsym,irrzon,iscf,natom,&
&  nfft,dtset%ngfft,nspden,nsppol,nsym,&
&  phnons,dtset%symafm,symrec,dtset%symrel,dtset%tnons,dtset%typat,xred)

! Make sure dtset%iatfix does not break symmetry
  call status(0,dtfil%filstat,iexit,level,'call fixsym   ')
  call fixsym(dtset%iatfix,indsym,natom,nsym)

 else

! The symrec array is used by initberry even in case nsym = 1
  symrec(:,:,1) = 0
  symrec(1,1,1) = 1 ; symrec(2,2,1) = 1 ; symrec(3,3,1) = 1

 end if

!Electric field: initialization
 if ((dtset%berryopt < 0).or.(dtset%berryopt == 4)) then
   nullify(pwind,pwnsfac)
   call initberry(dtefield,dtfil,dtset,gmet,kg,mband,mkmem,mpi_enreg,&
&              mpw,nkpt,npwarr,nsppol,nsym,occ,pwind,pwind_alloc,pwnsfac,rprimd,symrec)
 else
   pwind_alloc = 1
   allocate(pwind(pwind_alloc,2,3),pwnsfac(2,pwind_alloc))
 end if

!Timing for initialisation period
 call timab(33,2,tsec)
 call timab(34,1,tsec)

!Compute new occupation numbers, in case wavefunctions and eigenenergies
!were read from disk, occupation scheme is metallic (this excludes iscf=-1),
!and occupation numbers are required by iscf
 if( dtfil%ireadwf==1 .and. &
&    (dtset%occopt>=3.and.dtset%occopt<=7) .and. &
&    (iscf>0 .or. iscf==-3) ) then

  call status(0,dtfil%filstat,iexit,level,'call newocc   ')
  allocate(doccde(mband*nkpt*nsppol))
! Warning : ideally, results_gs%entropy should not be set up here XG 20011007
! Warning : ideally, results_gs%fermie should not be set up here XG 20011007
! Do not take into account the possible STM bias
  call newocc(doccde,eigen,results_gs%entropy,&
&  results_gs%fermie,&
&  dtset%fixmom,mband,dtset%nband,&
&  dtset%nelect,nkpt,nspinor,nsppol,occ,&
&  dtset%occopt,prtvol,zero,dtset%tphysel,dtset%tsmear,dtset%wtk)
  deallocate(doccde)

 else
! Warning : ideally, results_gs%entropy should not be set up here XG 20011007
  results_gs%entropy=0.0d0
 end if

!Generate an index table of atoms, in order for them to be used
!type after type.
 ntypat=psps%ntypat
 allocate(atindx(natom),atindx1(natom),nattyp(ntypat))
 index=1
 do itypat=1,ntypat
  nattyp(itypat)=0
  do iatom=1,natom
   if(dtset%typat(iatom)==itypat)then
    atindx(iatom)=index
    atindx1(index)=iatom
    index=index+1
    nattyp(itypat)=nattyp(itypat)+1
   end if
  end do
 end do

!PAW: 1- Initialize values for several arrays unchanged during iterations
!     2- initialize <g|S|c> (S=overlap matrix)
!     3- Eventually open temporary storage file
 if(psps%usepaw==1) then
! 1-
  if (psp_gencond==1) then
   call pawinit(psps%indlmn,dtset%ixc,psps%lmnmax,mpsang,dtset%pawnphi,nsym,&
&               dtset%pawntheta,psps%ntypat,pawang,pawrad,pawtab,dtset%pawxcdev)
   call setsymrhoij(gprimd,pawang%l_max-1,nsym,rprimd,&
&                   dtset%symrel,pawang%zarot)
  end if
! 2-
  call status(0,dtfil%filstat,iexit,level,'call getgsc   ')
  allocate(gsc(2,mpw*nspinor*mband*mkmem*nsppol))
  call getgsc(atindx,atindx1,cg,dtfil,dtset,gmet,&
&   gprimd,gsc,tmpfil(7),kg,mband,mgfft,mkmem,mpi_enreg,mpsang,&
&   mpw,natom,nattyp,nkpt,npwarr,nspinor,nsppol,ntypat,pawtab,&
&   psps,ucvol,dtfil%unylm,wffnew,wffnow,xred,ylm)
! 3-
  if(mkmem==0) then
   open(dtfil%unpaw,file=tmpfil(8),form='unformatted',status='unknown')
   rewind(unit=dtfil%unpaw)
  end if
 end if

!Get starting charge density : rhor as well as rhog
 allocate(rhog(2,nfft),rhor(nfft,nspden))
 if (iscf>0) then

  if(dtfil%ireadwf/=0)then

!  Obtain the charge density from wfs that were read previously
   call status(0,dtfil%filstat,iexit,level,'call mkrho    ')
!   tim_mkrho=1 ; mpi_enreg%paralbd=0
   tim_mkrho=1
   call mkrho(cg,densymop_gs,irrzon,dtset%istwfk,kg,mband,mgfft,mkmem,&
&   mpi_enreg,mpw,dtset%nband,nfft,dtset%ngfft,nkpt,npwarr,nspden,nspinor,&
&   nsppol,nsym,occ,phnons,rhog,rhor,dtset%symafm,tim_mkrho,ucvol,&
&   dtfil%unkg,wffnow,dtset%wtk)

  else if(dtfil%ireadwf==0)then

!  Crude, but realistic initialisation of the density
!  There is not point to compute it from random wavefunctions
!  Compute structure factor phases for current atomic pos:
   call status(0,dtfil%filstat,iexit,level,'call getph    ')
   allocate(ph1d(2,3*(2*mgfft+1)*natom))
   call getph(atindx,natom,dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3),&
&   ph1d,xred)
   call status(0,dtfil%filstat,iexit,level,'call initro   ')
   call initro(atindx,dtset%densty,gmet,gsqcut_eff,mgfft,natom,nattyp,&
&   nfft,dtset%ngfft,nspden,ntypat,ph1d,rhog,rhor,&
&   dtset%spinat,ucvol,dtset%ziontypat,dtset%znucl)
   deallocate(ph1d)

  end if

 else if (iscf==-1.or.iscf==-2.or.iscf==-3) then

  call status(0,dtfil%filstat,iexit,level,'call ioarr    ')
! Read rho(r) from a disk file
  rdwr=1
!  set to 1 for netcdf in/output
  accessfil = 0
! Note : results_gs%etotal is read here,
! and might serve in the tddft routine, but it is contrary to the
! intended use of results_gs ...
! Warning : should check the use of results_gs%fermie
! Warning : should check the use of results_gs%residm
! One might make them separate variables.

!DEBUG
!  write(6,*)' gstate : before ioarr, reading the density '
!ENDDEBUG

  call ioarr (accessfil,rhor,results_gs%etotal,fformr,dtfil%fildensin,hdr,&
&  nfft,nspden,rdwr)
! Compute up+down rho(G) by fft
  call status(0,dtfil%filstat,iexit,level,'call fourdp   ')
  allocate(work(nfft))
  work(:)=rhor(:,1)
  call fourdp(1,rhog,work,-1,nfft,dtset%ngfft,0)
  deallocate(work)

 else

! Disallowed value for iscf
  write(message, '(a,a,a,a,i12,a)' )  ch10,&
&   ' gstate : BUG -',ch10,&
&   '  iscf has disallowed value=',iscf,'.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')

 end if

!Debugging : print the different parts of rhor
! MPIWF Warning : this should not be parallelized over space, leave this debugging feature as such.
 if(prtvol==-level)then
  write(message,'(a)') '   ir     rhor(ir)     '
  call wrtout(06,message,'COLL')
  do ir=1,nfft
   if(ir<=11 .or. mod(ir,301)==0 )then
    write(message,'(i5,a,es13.6)')ir,' ',rhor(ir,1)
    call wrtout(06,message,'COLL')
    if(nsppol==2)then
     write(message,'(a,es13.6)')'      ',rhor(ir,2)
     call wrtout(06,message,'COLL')
    end if
   end if
  end do
 end if

 call status(0,dtfil%filstat,iexit,level,'end gstate(1) ')

 if(prtvol==-level)then
  write(message,'(a1,a,a1,a,i1,a)') ch10,&
&   ' gstate : before scfcv, move or brdmin ',&
&   ch10,'  prtvol=-',level,', debugging mode => stop '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 call timab(34,2,tsec)
!Check whether exiting was required by the user.
!If found then do not start minimization steps
!At this first call to chkexi, initialize cpus, if it
!is non-zero (which would mean that no action has to be taken)
!Should do this in driver ...
 cpus=dtset%cpus
 if(abs(cpus)>1.0d-5)cpus=cpus+cpui
 openexit=1 ; if(dtset%chkexit==0) openexit=0
 call chkexi(cpus,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg,openexit)
!If immediate exit, and wavefunctions were not read, must zero eigenvalues
 if (iexit/=0) then
  eigen(:)=zero
 end if
 if (iexit==0) then

! #######################################################################

! If atoms are not being moved, use scfcv directly; else
! call move or brdmin which in turn calls scfcv.

  call timab(35,1,tsec)

  write(message,'(a,80a)')ch10,('=',mu=1,80)
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,message,'COLL')
  if (ionmov==0) then

   call status(0,dtfil%filstat,iexit,level,'call scfcv    ')

!Should merge this call with the call for ionmov==4 and 5
   iapp=0
!   mpi_enreg%paralbd=0
   call scfcv(atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,&
&   ecore,eigen,gsc,hdr,iapp,indsym,initialized,&
&   irrzon,kg,mband,mgfft,mkmem,mpi_enreg,mpsang,mpw,natom,&
&   nattyp,nfft,nkpt,npwarr,nspden,nspinor,nsppol,nsym,occ,&
&   pawang,pawene,pawrad,pawtab,&
&   phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprimd,symrec,&
&   wffnew,wffnow,xred,xred_old,ylm)

  else if (ionmov==1) then
!  Conduct molecular dynamics, with or without viscous damping

   call status(0,dtfil%filstat,iexit,level,'call move     ')
!   mpi_enreg%paralbd=0
    call move(amass,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,&
&   ecore,eigen,gsc,hdr,indsym,initialized,irrzon,&
&   kg,mband,mgfft,mkmem,mpi_enreg,mpsang,mpw,&
&   natom,nattyp,nfft,nkpt,npwarr,nspden,nspinor,nsppol,nsym,occ,&
&   pawang,pawene,pawrad,pawtab,&
&   phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprimd,symrec,&
&   wffnew,wffnow,vel,xred,xred_old,ylm)

  else if (ionmov==2 .or. ionmov==3) then

!  Apply Broyden method for structural optimization, as
!  implemented by Jean-Christophe Charlier (May 1992)

   call status(0,dtfil%filstat,iexit,level,'call brdmin   ')
!   mpi_enreg%paralbd=0

   call brdmin(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,&
&   ecore,eigen,gsc,hdr,indsym,initialized,irrzon,&
&   kg,mband,mgfft,mkmem,mpi_enreg,mpsang,mpw,mxfh,&
&   natom,nattyp,nfft,nkpt,npwarr,nspden,nspinor,nsppol,nsym,nxfh,occ,&
&   pawang,pawene,pawrad,pawtab,&
&   phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprim,symrec,&
&   wffnew,wffnow,vel,xfhist,xred,xred_old,ylm)

  else if (ionmov==4 .or. ionmov==5) then

   do itime=1,ntime

    call status(itime,dtfil%filstat,iexit,level,'call scfcv(mv)')

    if(ionmov==4)then
     if(mod(itime,2)==1)then
      write(message, '(a,a,i3,a)' ) ch10,' STEP NUMBER ',itime,&
&      ' : OPTIMIZE ELECTRONS ------------------------------------'
     else
      write(message, '(a,a,i3,a)' ) ch10,' STEP NUMBER ',itime,&
&      ' : OPTIMIZE ELECTRONS AND IONS ---------------------------'
     end if
    else
     write(message, '(a,a,i3,a)' ) ch10,' STEP NUMBER ',itime,&
&     ' : SIMPLE RELAXATION -------------------------------------'
    end if
    call wrtout(ab_out,message,'COLL')
    call wrtout(06,  message,'COLL')

!   In this case, iapp is simply itime
    iapp=itime
!    mpi_enreg%paralbd=0
    call scfcv (atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,ecore,&
&    eigen,gsc,hdr,iapp,indsym,initialized,&
&    irrzon,kg,mband,mgfft,mkmem,mpi_enreg,mpsang,mpw,natom,&
&    nattyp,nfft,nkpt,npwarr,nspden,nspinor,nsppol,nsym,occ,&
&    pawang,pawene,pawrad,pawtab,&
&    phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprimd,symrec,&
&    wffnew,wffnow,xred,xred_old,ylm)

    if(mod(itime,2)==1)then
!    When the SCF cycle dealt with electrons only,
!    check whether forces are below tolerance; if so, exit
!    from the itime loop
     itimexit=0 ; if(itime==ntime)itimexit=1
     call fconv(results_gs%fcart,dtset%iatfix,itimexit,itime,natom,&
&     ntime,0,1.0d0,dtset%strtarget,results_gs%strten,dtset%tolmxf)
    end if
    if (itimexit/=0) exit

!   Check whether exiting was required by the user.
!   If found then beat a hasty exit from time steps
    openexit=1 ; if(dtset%chkexit==0) openexit=0
    call chkexi(cpus,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg,openexit)
    if (iexit/=0) then
     iexit=0   ! In order not to exit of dataset loop automatically
     exit
    end if

   end do

  else if (ionmov>=6 .and. ionmov<=9) then

!  Molecular dynamics, using Verlet algorithm (ionmov=6)
!  or fake molecular dynamics for minimisation (ionmov=7)
!  or true molecular dynamics with Nose thermostat (ionmov=8)
!  or Langevin dynamics (ionmov=9)

   call status(0,dtfil%filstat,iexit,level,'call moldyn   ')

   call moldyn(acell,amass,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,&
&   dtset,ecore,eigen,gsc,hdr,indsym,initialized,&
&   irrzon,kg,mband,mgfft,mkmem,mpi_enreg,mpsang,mpw,mxfh,&
&   natom,nattyp,nfft,nkpt,npwarr,nspden,nspinor,nsppol,nsym,nxfh,occ,&
&   pawang,pawene,pawrad,pawtab,&
&   phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprim,symrec,&
&   wffnew,wffnow,vel,xfhist,xred,xred_old,ylm)

  else if (ionmov == 10) then

   call delocint(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,&
&   dtset,ecore,eigen,gsc,hdr,indsym,initialized,irrzon,&
&   kg,mband,mgfft,mkmem,mpi_enreg,mpsang,mpw,mxfh,&
&   natom,nattyp,nfft,nkpt,npwarr,nspden,nspinor,nsppol,nsym,nxfh,occ,&
&   pawang,pawene,pawrad,pawtab,&
&   phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprim,symrec,&
&   wffnew,wffnow,vel,xfhist,xred,xred_old,ylm)

  else
!  Not an allowed option
   write(message, '(a,a,a,a,i12,a,a)' ) ch10,&
&   ' gstate : BUG -',ch10,&
&   '  Disallowed value for ionmov=',ionmov,ch10,&
&   '  Allowed values are 0 to 5.'
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if

  call timab(35,2,tsec)

! #####################################################################

!End of the check of hasty exit
 end if

 call timab(36,1,tsec)

 write(message, '(80a,a,a,a,a)' ) ('=',mu=1,80),ch10,ch10,&
&  ' ----iterations are completed or convergence reached----',&
&  ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')

!Close the unneeded temporary data files, if any.
!Other files are closed in clnup1.
 if (mkmem==0) then
  close (unit=dtfil%unkg,status='delete')
  if (psps%useylm==1) close (unit=dtfil%unylm,status='delete')
  if (psps%usepaw==1) close (unit=dtfil%ungsc,status='delete')
  if (psps%usepaw==1) close (unit=dtfil%unpaw,status='delete')
  call WffDelete(wffnew,ierr)
 end if

!Update the header, before using it
 call hdr_update(bantot,results_gs%etotal,results_gs%fermie,hdr,natom,&
&                results_gs%residm,rprimd,occ,xred)

 call status(0,dtfil%filstat,iexit,level,'call outwf    ')

 call outwf(cg,dtfil,dtset,eigen,dtfil%filnam_ds(4),hdr,kg,dtset%kptns,&
& mband,mkmem,mpi_enreg,mpw,mxfh,natom,dtset%nband,nfft,dtset%ngfft,&
& nkpt,npwarr,dtset%nqpt,nspinor,nsppol,dtset%nstep,nxfh,&
& occ,resid,response,wffnow,xfhist)

! CASINO
 call outqmc(cg,dtset,eigen,gprimd,hdr,kg,mband,mkmem,mpw,nkpt,npwarr,&
             &nspinor,nsppol,occ,psps,results_gs)
! CASINO

 call status(0,dtfil%filstat,iexit,level,'call clnup1   ')

 call clnup1(acell,dtset%dosdeltae,eigen,dtset%enunit,&
& results_gs%fermie,dtfil%filnam_ds(4),&
& results_gs%fred,dtset%iatfix,iscf,dtset%kptns,mband,mkmem,mpi_enreg,mpw,&
& natom,dtset%nband,nfft,dtset%ngfft,nkpt,nspden,nspinor,nsppol,dtset%nstep,&
& occ,dtset%occopt,dtset%prtdos,dtset%prteig,dtset%prtstm,prtvol,&
& resid,rhor,rprimd,dtset%tphysel,dtset%tsmear,results_gs%vxcavg,dtset%wtk,xred)

 if (iscf>0 .and. dtset%prtstm==0) then
   call status(0,dtfil%filstat,iexit,level,'call prtene   ')
   call prtene(dtset%berryopt,results_gs%eei,results_gs%eeig,results_gs%eew,&
&       results_gs%ehart,results_gs%eii,results_gs%ek,&
&       results_gs%enefield,results_gs%enl,&
&       results_gs%entropy,results_gs%enxc,pawene%etilddc,&
&       pawene%e1dc-pawene%etild1dc,pawene%etot,pawene%etotdc,&
&       epulay,ab_out,dtset%occopt,dtset%tphysel,dtset%tsmear,psps%usepaw)
 end if

! Open the formatted derivative database file, and write the
! preliminary information
! In the // case, only one processor writes the energy and
! the gradients to the DDB

 if ((psps%usepaw == 0).and.(mpi_enreg%me==0).and.((iscf > 0).or.&
&      (dtset%berryopt == -1).or.(dtset%berryopt) == -3)) then

  call status(0,dtfil%filstat,iexit,level,'call ioddb8   ')
  vrsddb=010929
  dscrpt=' Note : temporary (transfer) database '
  choice=2
  ddbnm=trim(dtfil%filnam_ds(4))//'_DDB'
! tolwfr must be initialized here, but it is a dummy value
  tolwfr=1.0d0
  call ioddb8 (choice,dscrpt,ddbnm,natom,mband,&
&   nkpt,nsym,psps%ntypat,dtfil%unddb,vrsddb,&
&   acell,dtset%amu,dtset%dilatmx,dtset%ecut,dtset%ecutsm,&
&   dtset%intxc,dtset%iscf,dtset%ixc,dtset%kpt,dtset%kptnrm,&
&   natom,dtset%nband,dtset%ngfft,nkpt,nspden,nspinor,&
&   nsppol,nsym,psps%ntypat,occ,dtset%occopt,&
&   rprim,dtset%sciss,dtset%spinat,dtset%symafm,dtset%symrel,&
&   dtset%tnons,tolwfr,dtset%tphysel,dtset%tsmear,&
&   dtset%typat,dtset%wtk,xred,psps%ziontypat,dtset%znucl)

  if (iscf > 0) then
    nblok = 2          ! 1st blok = energy, 2nd blok = gradients
  else
    nblok = 1
  end if
  fullinit = 0
  call psddb8 (choice,psps%dimekb,psps%ekb,fullinit,psps%indlmn,&
&   psps%lmnmax,psps%lnmax,nblok,&
&   psps%ntypat,dtfil%unddb,psps%pspso,psps%usepaw,psps%useylm,vrsddb)

  mpert = natom + 6   ; msize = 3*mpert
  allocate(blkflg(msize),blkval(2,msize))

  blkflg(:) = 0       ; blkval(:,:) = zero
  blkqpt(:) = zero    ; blknrm(:) = one

! Write total energy to the DDB
  if (iscf > 0) then
    blktyp = 0
    blkval(1,1) = results_gs%etotal
    blkflg(1) = 1
    call blok8(blkflg,blknrm,blkqpt,blktyp,blkval,choice,mpert,&
&      msize,natom,dtfil%unddb)
  end if

! Write gradients to the DDB
  blktyp = 4
  blkflg(:) = 0       ; blkval(:,:) = zero
  index = 0
  if (iscf > 0) then
    do iatom = 1, natom
      do idir = 1, 3
        index = index + 1
        blkflg(index) = 1
        blkval(1,index) = results_gs%fred(idir,iatom)
      end do
    end do
  end if

  index = 3*natom + 3
  if ((abs(dtset%berryopt) == 1).or.(abs(dtset%berryopt) == 3)) then
    do idir = 1, 3
      index = index + 1
      if (dtset%rfdir(idir) == 1) then
        blkflg(index) = 1
        blkval(1,index) = results_gs%pel(idir)
      end if
    end do
  end if

  index = 3*natom + 6
  if (iscf > 0) then
    blkflg(index+1:index+6) = 1
    blkval(1,index+1:index+6) = results_gs%strten(1:6)
  end if

  call blok8(blkflg,blknrm,blkqpt,blktyp,blkval,choice,mpert,&
&    msize,natom,dtfil%unddb)

  deallocate(blkflg,blkval)

! Close DDB
  close(dtfil%unddb)

 end if

 if (dtset%nstep>0 .and. dtset%prtstm==0) then
  call status(0,dtfil%filstat,iexit,level,'call clnup2   ')
  call clnup2(psps%n1xccc,results_gs%fred,results_gs%gresid,&
&  results_gs%grewtn,&
&  results_gs%grxc,iscf,natom,prtvol,start,&
&  results_gs%strten,results_gs%synlgr,xred)
 end if

!Deallocate arrays
 deallocate(amass,atindx,atindx1,cg,eigen,indsym)
 deallocate(irrzon,kg,npwarr,nattyp,phnons,resid)
 deallocate(rhog,rhor,start,symrec,xfhist,xred_old,ylm)
 if(psps%usepaw==1) deallocate(gsc)
 if ((dtset%berryopt < 0).or.(dtset%berryopt == 4)) then
   deallocate(pwind,pwnsfac)
   deallocate(dtefield%ikpt_dk,dtefield%idxkstr)
   deallocate(dtefield%sflag,dtefield%cgindex)
   deallocate(dtefield%fkptns,dtefield%indkk_f2ibz,dtefield%i2fbz)
   if (mpi_enreg%paral_compil_kpt == 1) then
    deallocate(mpi_enreg%kptdstrb)
    if (dtset%berryopt == 4) then
     deallocate(mpi_enreg%kptdstrbi,dtefield%cgqindex,dtefield%nneigh)
    end if
   end if
 end if
 if (dtset%berryopt == 4) deallocate(dtefield%smat)

!Clean the header
 call hdr_clean(hdr)

!Clean the MPI informations
 if (mpi_enreg%parareel == 0) then
   call clnmpi_gs(mpi_enreg)
 end if

#         if defined MPI_FFT
         call clnmpi_fft(mpi_enreg)

#         endif

 write(message, '(a,a)' ) ch10,' gstate : exiting '
 call wrtout(06,message,'COLL')

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(36,2,tsec)
 call timab(32,2,tsec)

 end subroutine
!!***

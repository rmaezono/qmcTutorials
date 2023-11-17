!----------------------------------------------------------------------!
! REBLOCK                                                              !
! -------                                                              !
! VMC/DMC statistical analysis code.                                   !
! Averaging and reblocking of vmc.hist and dmc.hist/dmc.hist2 files.   !
! MDT 4.2000                                                           !
!                                                                      !
! Changes                                                              !
! -------                                                              !
! 6.2000  MDT Changed to allow use of MPC or Ewald interactions.       !
! 3.2001  MDT Modified to follow conversion of CASINO output from      !
!         eV/atom to au/primitive cell(or molecule). Added read of     !
!         eionion and electron-hole masses in VMC.                     !
! 3.2001  MDT Added read of dmc.hist2 file and analysis of DMC energy  !
!         components as well as total energy. Changed format of        !
!         dmc.hist by including number of ions in primitive cell and   !
!         periodicity to increase printout flexibility in reblock.     !
! 3.2001  MDT When analysing DMC results, reblock tries to find an     !
!         xwfn.data file (x=g/b/p) and extracts the n-n repulsion      !
!         energy if it can find one. Also tries to find input to       !
!         look for the number of holes.                                !
! 4.2001  MDT With advent of new DMC algorithm, require new way of     !
!         computing number of equilibration steps in a dmc.hist file,  !
!         since eref no longer equal to etrial during accumulation.    !
!         Therefore implemented detection of DMC algorithm type (from  !
!         input file) and a new method of computing nequil if it       !
!         finds dmc_algorithm==3                                       !
! 6.2001  MDT Modified behaviour when dmc.hist2 contains both          !
!         MPC and Ewald data.                                          !
! 2.2002  MDT Modified format of output to more closely match CASINO   !
! 2.2002  MDT Got rid of the different behaviour with different DMC    !
!         algorithms, since now we have only one.                      !
! 3.2002  MDT Made sure reblock understands about core polarization    !
!         energies in periodic systems.                                !
! 3.2002  NDD Allow it to calculate equilibration period for new DMC   !
!         algorithm (4). Gives info on growth and mixed estimators of  !
!         energy.                                                      !
! 7.2002  NDD Utilise full range of data when computing std errors.    !
! 7.2002  NDD Got rid of NAG subroutines and rearranged.               !
! 1.2003  NDD Use total weight at each generation (if available).      !
! 4.2003  NDD Use effective timestep info (if available).              !
! 12.2003 MDT Added analysis of relativistic corrections.              !
! 12.2003 MDT Tidied. Added i2s/r2s for string handling.               !
! 12.2004 MDT Allowed automatic reblocking without the program asking  !
!             questions if reblock.in file is present                  !
! 12.2004 MDT Added kcal/mol as an optional unit for molecules.        !
! 12.2004 MDT Prevented mean energy from being printed out for every   !
!             block length (it doesn't change!).                       !
!  1.2005 AM  Modified analysis of relativistic corrections.           !
!----------------------------------------------------------------------!


MODULE qmcdata
!------------------------!
! Some global variables. !
!------------------------!
 INTEGER nitems,nitems2,nlines,nlines2,nbasis,periodicity,nele(4),       &
  &dmc_algorithm,iterac,nelec,reblock_in_units,reblock_in_blocklength,   &
  &reblock_in_iterac,reblock_in_nstart,reblock_in_nave
 DOUBLE PRECISION,ALLOCATABLE :: items(:,:),items2(:,:)
 DOUBLE PRECISION,PARAMETER :: htoev=27.2113962d0,htokcal=627.507541278d0
 DOUBLE PRECISION eionion,mass_e,mass_h
 LOGICAL vmc,dmc,dmc2,reblock_in,have_total_weights
END MODULE qmcdata


MODULE vmc_records
!---------------------------------------------------!
! What's what in the columns of the vmc.hist file.  !
!---------------------------------------------------!
 INTEGER hist_etot
 INTEGER,PARAMETER :: hist_ewald=1,lpei=2,lkei=3,lti=4,lfisq=5,lpote=6,  &
  &lpoti=7,nlpoti=8,lnewee=9,lshort=10,llong=11,lvcppei=12,lvcppe=13,    &
  &lvcppee=14,hist_mpc=15 ! most of this isn't used any more.
END MODULE vmc_records


MODULE dmc_records
!-------------------------------------------------------------------!
! What's what in the columns of the dmc.hist and dmc.hist2 files.   !
!-------------------------------------------------------------------!
 INTEGER,PARAMETER :: hist_pop=2,hist_etot=3,hist_eref=4,hist_ebest=5,   &
  &hist_acceptr=6,hist_efficiency=9,hist_nbasis=10,hist_periodicity=11,  &
  &hist_egrowth=12,hist_totmult=13
 INTEGER,PARAMETER :: hist2_ti=1,hist2_kei=2,hist2_fi2=3,hist2_pote=4,   &
  &hist2_smpc=5,hist2_lmpc=6,hist2_nloc=7,hist2_loc=8,hist2_vcppei=9,    &
  &hist2_vcppe=10,hist2_vcppee=11,hist2_avwt=12,hist2_minwt=13,          &
  &hist2_maxwt=14,hist2_avage=15,hist2_maxage=16,hist2_dteff=17
END MODULE dmc_records


MODULE stats_calcs
!-------------------------------------------------------------------!
! A collection of subroutines for extracting means, variances, etc  !
! from supplied data.                                               !
!-------------------------------------------------------------------!
 IMPLICIT NONE


CONTAINS


 SUBROUTINE compute_av_var_unweighted(n,data_arr,want_var,av,var)
!----------------------------------------------------------------!
! Compute the mean and (if desired) variance of a set of data.   !
!----------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: n
  DOUBLE PRECISION,INTENT(in) :: data_arr(1:n)
  DOUBLE PRECISION,INTENT(out) :: av,var
  LOGICAL,INTENT(in) :: want_var
  INTEGER i

  if(n>=1)then

! Compute average of data if we have at least 1 point.
   av=sum(data_arr(1:n))/dble(n)

   if(want_var)then
    if(n>=2)then

! If desired, compute variance if we have at least 2 points.
     var=0.d0
     do i=1,n
      var=var+(data_arr(i)-av)**2
     enddo
     var=var/dble(n-1)

    else
     write(6,*)'Can''t compute variance of less than 2 points.'
     stop
    endif ! n>=2
   endif ! want_var

  else
   write(6,*)'Can''t compute average of no data.'
   stop
  endif ! n>=1

 END SUBROUTINE compute_av_var_unweighted


 SUBROUTINE compute_av_var_weighted(n,data_arr,wt_arr,want_var,av,var)
!-----------------------------------------------------------------------!
! Compute the mean and (if desired) variance of a set of weighted data. !
!-----------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: n
  DOUBLE PRECISION,INTENT(in) :: data_arr(1:n),wt_arr(1:n)
  DOUBLE PRECISION,INTENT(out) :: av,var
  LOGICAL,INTENT(in) :: want_var
  INTEGER i
  DOUBLE PRECISION tot_weight,red_tot_weight,sum_wt_sq

  if(n>=1)then

! If we have more than one point, evaluate the total weight.

   if(any(wt_arr<0.d0))then
    write(6,*)'Statistical weights are negative. Stopping.'
    stop
   endif

   tot_weight=sum(wt_arr(1:n))

   if(tot_weight==0.d0)then
    write(6,*)'The total weight is zero. Stopping.'
    stop
   endif

! Evaluate the average of the data.
   av=0.d0
   do i=1,n
    av=av+data_arr(i)*wt_arr(i)
   enddo
   av=av/tot_weight

   if(want_var)then
    if(n>=2)then

! If desired and we have more that 2 data points, compute the variance.

     sum_wt_sq=0.d0
     do i=1,n
      sum_wt_sq=sum_wt_sq+wt_arr(i)**2
     enddo
     red_tot_weight=tot_weight-sum_wt_sq/tot_weight

     var=0.d0
     do i=1,n
      var=var+(data_arr(i)-av)**2*wt_arr(i)
     enddo
     var=var/red_tot_weight

    else
     write(6,*)'Can''t compute variance of fewer than 2 points.'
     stop
    endif ! n>=2
   endif ! want_var

  else
   write(6,*)'Can''t compute average of no data.'
   stop
  endif ! n>=1

 END SUBROUTINE compute_av_var_weighted


 SUBROUTINE compute_av_var_one_wt(n,data_arr,wt_last_item,want_var,av,var)
!-----------------------------------------------------------------------!
! Compute the mean and (if desired) variance of a set of weighted data  !
! in which all the weights except that of the last datum are unity.     !
!-----------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: n
  DOUBLE PRECISION,INTENT(in) :: data_arr(1:n),wt_last_item
  DOUBLE PRECISION,INTENT(out) :: av,var
  LOGICAL,INTENT(in) :: want_var
  INTEGER i
  DOUBLE PRECISION tot_weight,red_tot_weight,sum_wt_sq

  if(n>=1)then

! If we have more than one point, evaluate the total weight.

   if(wt_last_item<0.d0)then
    write(6,*)'Statistical weight is negative. Stopping.'
    stop
   endif

   tot_weight=dble(n-1)+wt_last_item

   if(tot_weight==0.d0)then
    write(6,*)'The total weight is zero. Stopping.'
    stop
   endif

! Evaluate the average of the data.
   av=0.d0
   do i=1,n
    if(i<n)then
     av=av+data_arr(i)
    else
     av=av+data_arr(i)*wt_last_item
    endif
   enddo
   av=av/tot_weight

   if(want_var)then
    if(n>=2)then

! If desired and we have more that 2 data points, compute the variance.

     sum_wt_sq=dble(n-1)+wt_last_item**2
     red_tot_weight=tot_weight-sum_wt_sq/tot_weight

     var=0.d0
     do i=1,n
      if(i<n)then
       var=var+(data_arr(i)-av)**2
      else
       var=var+(data_arr(i)-av)**2*wt_last_item
      endif
     enddo
     var=var/red_tot_weight

    else
     write(6,*)'Can''t compute variance of fewer than 2 points.'
     stop
    endif ! n>=2
   endif ! want_var

  else
   write(6,*)'Can''t compute average of no data.'
   stop
  endif ! n>=1

 END SUBROUTINE compute_av_var_one_wt


 SUBROUTINE compute_stats_unweighted(n,data_arr,av,var,skew,kurt,max_val, &
  &min_val)
!-------------------------------------------------------------------!
! Compute mean, variance, skewness, kurtosis and max and min of a   !
! set of data.                                                      !
!-------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: n
  DOUBLE PRECISION,INTENT(in) :: data_arr(1:n)
  DOUBLE PRECISION,INTENT(out) :: av,var,skew,kurt,max_val,min_val
  INTEGER i
  DOUBLE PRECISION sum_delta_x2,sum_delta_x3,sum_delta_x4

  if(n<2)then
   write(6,*)'Can''t compute variance with fewer than two points.'
   stop
  endif

! Compute average of data
  av=sum(data_arr(1:n))/dble(n)

! Compute max and min.
  max_val=maxval(data_arr(1:n))
  min_val=minval(data_arr(1:n))

! Compute variance, skewness and kurtosis.
  sum_delta_x2=0.d0  ;  sum_delta_x3=0.d0  ;  sum_delta_x4=0.d0
  do i=1,n
   sum_delta_x2=sum_delta_x2+(data_arr(i)-av)**2
   sum_delta_x3=sum_delta_x3+(data_arr(i)-av)**3
   sum_delta_x4=sum_delta_x4+(data_arr(i)-av)**4
  enddo

  var=sum_delta_x2/dble(n-1)
  if(var>0.d0)then
   skew=sum_delta_x3/(dble(n-1)*var**1.5d0)
   kurt=sum_delta_x4/(dble(n-1)*var**2)-3.d0
  else
   skew=0.d0
   kurt=0.d0
  endif

 END SUBROUTINE compute_stats_unweighted

END MODULE stats_calcs


MODULE io_module
!-------------------------------------------------------!
! Subroutines for reading files and displaying results. !
!-------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,PARAMETER :: r2s_length=80
 CHARACTER(r2s_length) tmpr,tmpr2,tmpr3,tmpr4,tmpr5


CONTAINS


 SUBROUTINE read_vmc_data(iu,iout)
!--------------------------------!
! Read data from vmc.hist file.  !
!--------------------------------!
  USE qmcdata
  USE vmc_records
  IMPLICIT NONE
  INTEGER,INTENT(in) :: iu,iout
  INTEGER i,j,k,nmove,nblock,mstart,nmove_old,ialloc
  CHARACTER(13), PARAMETER :: reals4="(4(1PE20.13))"

  ialloc=0
  read(iu,*,end=3)
  rewind(iu)
! Figure out number of items per record (=nave)
  read(iu,*,err=2,end=1)
  read(iu,*,err=2,end=2)nmove,nbasis,nele(1),nele(2),periodicity,nitems,eionion
  rewind(iu)
  allocate(items(nitems,1),stat=ialloc) ! temporarily
  if(ialloc/=0)then
   write(iout,*)'Error allocating items array. Sorry.'
   stop
  endif
! Figure out total number of records
  nblock=0
  do ! till we hit the end of the file which should be at end of block
   read(iu,*,err=2,end=1)
   nblock=nblock+1
   read(iu,*,err=2,end=2)nmove,nbasis,nele(1),nele(2),periodicity,nitems,eionion
   if(nbasis==0)then
    read(iu,*,err=2,end=2)nele(3),nele(4),mass_e,mass_h
   endif
   do j=1,nmove
    read(iu,reals4,err=2,end=2)(items(i,1),i=1,nitems)
   enddo
  enddo

1 write(iout,*)
  write(iout,'(1x,2a)')'Number of data blocks read : ',trim(i2s(nblock))
  write(iout,'(1x,2a)')'Number of moves per block  : ',trim(i2s(nmove))
  write(iout,'(1x,2a)')'Number of records per move : ',trim(i2s(nitems))

  nlines=nblock*nmove
  if(nlines==0)then
   write(6,*)
   write(6,*)'So I''ll stop right there.'
   stop
  endif
  deallocate(items)
  allocate(items(nitems,nlines),stat=ialloc)
  if(ialloc/=0)then
   write(iout,*)'Error allocating items array. Sorry.'
   stop
  endif
  items=0.d0

! Now we've allocated sufficient space, read the data properly
  rewind(iu)
  mstart=0
  do k=1,nblock
   nmove_old=nmove
   read(iu,*,err=2,end=1)
   read(iu,*,err=2,end=1)nmove,nbasis,nele(1),nele(2),periodicity,nitems,eionion
   if(nmove/=nmove_old)then
    write(6,*)'Number of moves per block changed inside vmc.hist file. &
     &Stopping.'
    stop
   endif
   if(nbasis==0)then
    read(iu,*,err=2,end=2)nele(3),nele(4),mass_e,mass_h
   endif
   do j=1,nmove
    read(iu,reals4,err=2,end=1)(items(i,mstart+j),i=1,nitems)
   enddo
   mstart=mstart+nmove
  enddo

  return
2 write(iout,*)'Error reading vmc.hist file. Quitting.'
  stop
3 write(iout,*)
  write(iout,*)'The file is empty. Quitting.'
  stop

 END SUBROUTINE read_vmc_data


 SUBROUTINE read_dmc_data(iu,repu)
!-------------------------------!
! Read data from dmc.hist file. !
!-------------------------------!
  USE dmc_records
  USE qmcdata
  IMPLICIT NONE
  INTEGER,INTENT(in) :: iu,repu
  INTEGER i,j,ierr
  INTEGER,PARAMETER :: maxitem = 128 ! maximum no of items on one line of file
  DOUBLE PRECISION tempdata(maxitem)
  CHARACTER(512) cline

! Deduce number of items per line
  ierr=0
  rewind(iu)
  cline(:)=' '
  read(iu,fmt='(a)',iostat=ierr)cline
  if(ierr/=0)then
   write(repu,*)'Problem reading 1st line of dmc.hist file.'
   stop
  endif
  do nitems=0,maxitem-1,1
   read(cline,fmt=*,err=1,end=1)tempdata(1:nitems+1)
  enddo
1 continue
  if(nitems==0)then
   write(repu,*)'No data on first line of dmc.hist? I/O error.'
   stop
  endif
  write(repu,'(1x,2a)')trim(i2s(nitems)),' items per line of data.'
  allocate(items(nitems,1))
  nlines=1

! Work out how many lines of data we have
  do
   read(iu,fmt=*,end=3,err=3)(items(i,1),i=1,nitems)
   nlines=nlines+1
  enddo
3 continue

! Re-allocate arrays to store this data
  deallocate(items)
  allocate(items(nitems,nlines))
  items=0.d0

  rewind(iu)
  do j=1,nlines
   read(iu,fmt=*,iostat=ierr)(items(i,j),i=1,nitems)
   if(ierr/=0)then
    write(repu,*)'Problem reading dmc.hist at line no.',j
    stop
   endif
  enddo
  write(repu,'(1x,2a)')trim(i2s(nlines)),' lines of data.'

 END SUBROUTINE read_dmc_data


 SUBROUTINE read_dmc2_data(iu,repu)
!---------------------------------!
! Read data from dmc.hist2 file.  !
!---------------------------------!
  USE dmc_records
  USE qmcdata
  IMPLICIT NONE
  INTEGER,INTENT(in) :: iu,repu
  INTEGER i,j,ierr
  INTEGER,PARAMETER :: maxitem=128 ! maximum no of items on one line of file
  DOUBLE PRECISION tempdata(maxitem)
  CHARACTER(512) cline

! Deduce number of items per line
  ierr=0
  rewind(iu)
  cline(:)=' '
  read(iu,fmt='(a)',iostat=ierr)cline
  if(ierr/=0)then
   write(repu,*)'Problem reading 1st line of dmc.hist2 file.'
   stop
  endif
  do nitems2=0,maxitem-1,1
   read(cline,fmt=*,err=1,end=1)tempdata(1:nitems2+1)
  enddo
1 continue
  if(nitems2==0)then
   write(repu,*)'No data on first line of dmc.hist2? I/O error.'
   stop
  endif
  write(repu,'(1x,2a)')trim(i2s(nitems2)),' items per line of data.'
  allocate(items2(nitems2,1))
  nlines2=1

! Work out how many lines of data we have
  do
   read(iu,fmt=*,end=3,err=3)(items2(i,1),i=1,nitems2)
   nlines2=nlines2+1
  enddo
3 continue

! Re-allocate arrays to store this data
  deallocate(items2)
  allocate(items2(nitems2,nlines2))

  rewind(iu)
  do j=1,nlines2
   read(iu,fmt=*,iostat=ierr)(items2(i,j),i=1,nitems2)
   if(ierr/=0)then
    write(repu,*)'Problem reading dmc.hist2 at line no. ',j
    stop
   endif
  enddo
  write(repu,'(1x,2a)')trim(i2s(nlines2)),' lines of data.'

 END SUBROUTINE read_dmc2_data


 SUBROUTINE write_out_vmc(n,no_blocks,blocklength,average,stderr,isperiodic, &
  &e_units)
!---------------------------------------------!
! Print out final reblocked results for VMC.  !
!---------------------------------------------!
  USE qmcdata, ONLY : eionion,nbasis,periodicity
  USE vmc_records, ONLY : hist_etot,hist_mpc,hist_ewald
  IMPLICIT NONE
  INTEGER,INTENT(in) :: n,no_blocks,blocklength
  DOUBLE PRECISION,INTENT(in) :: average(n),stderr(n)
  CHARACTER(9),INTENT(in) :: e_units
  LOGICAL,INTENT(in) :: isperiodic

  write(6,*)
  write(6,*)'FINAL VMC AVERAGES AND ERROR BARS'
  write(6,3)
  write(6,*)'Mean values averaged over all data'
  write(6,1)trim(i2s(no_blocks)),trim(i2s(blocklength))
1 format(1x,'Std. errors from averages over ',a,' blocks of length ',a,'.',&
   &/1x,76('-'))
  write(6,2)
2 format(48x,'Mean',16x,'Std. Error')
3 format(1x,76('-'))
  if (isperiodic) then
   if(hist_etot==hist_ewald) then
    write(6,4)'Total energy using Ewald  (',trim(e_units),&
     &') : ',average(1),stderr(1)
   elseif(hist_etot==hist_mpc)then
    write(6,4)'Total energy using MPC    (',trim(e_units),&
     &') : ',average(15),stderr(15)
   else
    write(6,*)'Type of interaction? Bug.'
    stop
   endif
  else
   write(6,4)'Total energy              (',trim(e_units), &
    &') : ',average(1),stderr(1)
   if(n>15)then
    write(6,4)'Relativistic total energy (',trim(e_units), &
     &') : ',average(1)+average(21),sqrt(stderr(1)**2+stderr(21)**2)
   endif
  endif
  write(6,4)'PEI                       (',trim(e_units),  &
   &') : ',average(2),stderr(2)
  write(6,4)'KEI                       (',trim(e_units),  &
   &') : ',average(3),stderr(3)
  write(6,4)'TI                        (',trim(e_units),  &
   &') : ',average(4),stderr(4)
  write(6,4)'FISQ                      (',trim(e_units),  &
   &') : ',average(5),stderr(5)
  if(nbasis>0)then
   write(6,4)'e-i interaction (local)   (',trim(e_units), &
    &') : ',average(7),stderr(7)
   write(6,4)'                (nonlocal)(',trim(e_units), &
    &') : ',average(8),stderr(8)
  endif
  if(isperiodic) then
   write(6,4)'Ewald e-e interaction     (',trim(e_units), &
    &') : ',average(6),stderr(6)
   write(6,4)'1/r within W.S. cell      (',trim(e_units), &
    &') : ',average(10),stderr(10)
   write(6,4)'Hartree outside W.S. cell (',trim(e_units), &
    &') : ',average(11),stderr(11)
   write(6,4)'MPC interaction           (',trim(e_units), &
    &') : ',average(9),stderr(9)
  else
   write(6,4)'e-e interaction           (',trim(e_units), &
    &') : ',average(6),stderr(6)
  endif
  if(any(average(12:14)/=0.d0))then
   if(periodicity==0)then
    write(6,4)'Core polarization (e-i)   (',trim(e_units),&
     &') : ',average(12),stderr(12)
    write(6,4)'Core polarization (e)     (',trim(e_units),&
     &') : ',average(13),stderr(13)
    write(6,4)'Core polarization (e-e)   (',trim(e_units),&
     &') : ',average(14),stderr(14)
   else
    write(6,4)'Core polarization energy  (',trim(e_units),&
     &') : ',average(12),stderr(12)
   endif
  endif
  if(nbasis>0)then
   write(6,5)'N-N repulsion             (',trim(e_units), &
    &') : ',eionion
  endif
  if(n>15)then ! relativistic components present
   write(6,*)
   write(6,*)'Relativistic corrections'
   write(6,*)
   write(6,4)'Mass polarization         (',trim(e_units),&
    &') : ',average(16),stderr(16)
   write(6,4)'Mass velocity term        (',trim(e_units),&
    &') : ',average(17),stderr(17)
   write(6,4)'Electron-nucleus Darwin   (',trim(e_units),&
    &') : ',average(18),stderr(18)
   write(6,4)'Electron-electron Darwin  (',trim(e_units),&
    &') : ',average(19),stderr(19)
   write(6,4)'Retardation term          (',trim(e_units),&
    &') : ',average(20),stderr(20)
   write(6,4)'Total rel. correction     (',trim(e_units),&
    &') : ',average(21),stderr(21)
  endif

4 format(t2,a,a,a,t38,f21.12,f21.12)
5 format(t2,a,a,a,t38,f21.12)

 END SUBROUTINE write_out_vmc


 SUBROUTINE write_out_dmc(n,no_blocks,blocklength,e,v,e2,v2,average,stderr, &
  &isperiodic,iterac,got_eionion,e_units)
!---------------------------------------------!
! Write out final reblocked results for DMC.  !
!---------------------------------------------!
  USE dmc_records
  USE qmcdata, ONLY : eionion,nbasis,dmc2,periodicity
  IMPLICIT NONE
  INTEGER,INTENT(in) :: n,no_blocks,blocklength,iterac
  DOUBLE PRECISION,INTENT(in) :: average(n),stderr(n),e,v,e2,v2
  CHARACTER(9),INTENT(in) :: e_units
  LOGICAL,INTENT(in) :: isperiodic,got_eionion
  DOUBLE PRECISION temp

  write(6,*)
  write(6,*)'FINAL DMC AVERAGES AND ERROR BARS'
  write(6,3)
  write(6,*)'Mean values averaged over all statistics accumulation data.'
  write(6,1)trim(i2s(no_blocks)),trim(i2s(blocklength))
1 format(1x,'Std. errors from averages over ',a,' blocks of length ',a,'.',&
   &/1x,76('-'))
  write(6,2)
2 format(48x,'Mean',15x,'Std. Error')
3 format(1x,76('-'))
  if(.not.dmc2)then
   write(6,4)'Total energy              (',trim(e_units),') : ',e,v
   write(6,*)
   write(6,*)'(No dmc.hist2 file ==> no energy components)'
   write(6,*)
  else
   if(iterac/=3.and.iterac/=4)then
    write(6,4)'Total energy              (',trim(e_units),') : ',e,v
    if(n>17)then ! relativistic components present
     write(6,4)'Relativistic Total energy (',trim(e_units),') : ', &
      &e+average(23),sqrt(v**2+stderr(23)**2)
    endif
   else
    if(iterac==3)then
     write(6,4)'Total energy (Ewald)      (',trim(e_units),') : ',e,v
     write(6,4)'Total energy (MPC)        (',trim(e_units),') : ',e2,v2
    else ! iterac==4
     write(6,4)'Total energy (MPC)        (',trim(e_units),') : ',e,v
     write(6,4)'Total energy (Ewald)      (',trim(e_units),') : ',e2,v2
    endif
   endif
  endif
  if(dmc2)then
   write(6,4)'KEI                       (',trim(e_units),') : ', &
    &average(hist2_kei),stderr(hist2_kei)
   write(6,4)'TI                        (',trim(e_units),') : ', &
    &average(hist2_ti),stderr(hist2_ti)
   write(6,4)'FISQ                      (',trim(e_units),') : ', &
    &average(hist2_fi2),stderr(hist2_fi2)
   if(nbasis>0)then
    write(6,4)'e-i interaction (local)   (',trim(e_units),') : ', &
     &average(hist2_loc),stderr(hist2_loc)
    write(6,4)'                (nonlocal)(',trim(e_units),') : ', &
     &average(hist2_nloc),stderr(hist2_nloc)
   endif
   if(isperiodic) then
    if(average(hist2_pote)/=0.d0)then
     write(6,4)'Ewald e-e interaction     (',trim(e_units),') : ', &
      &average(hist2_pote),stderr(hist2_pote)
    endif
    if(average(hist2_smpc)/=0.or.average(hist2_lmpc)/=0)then
     write(6,4)'1/r within W.S. cell      (',trim(e_units),') : ', &
      &average(hist2_smpc),stderr(hist2_smpc)
     write(6,4)'Hartree outside W.S. cell (',trim(e_units),') : ', &
      &average(hist2_lmpc),stderr(hist2_lmpc)
     temp=sqrt(stderr(hist2_smpc)**2+stderr(hist2_lmpc)**2)
     write(6,4)'MPC e-e interaction       (',trim(e_units),') : ', &
      &average(hist2_smpc)+average(hist2_lmpc),temp
    endif
   else
    write(6,4)'e-e interaction           (',trim(e_units),') : ', &
     &average(hist2_pote),stderr(hist2_pote)
   endif
   if(any(average(hist2_vcppei:hist2_vcppee)/=0.d0))then
    if(periodicity==0)then
     write(6,4)'e-i polarization          (',trim(e_units),') : ', &
      &average(hist2_vcppei),stderr(hist2_vcppei)
     write(6,4)'e   polarization          (',trim(e_units),') : ', &
      &average(hist2_vcppe),stderr(hist2_vcppe)
     write(6,4)'e-e polarization          (',trim(e_units),') : ', &
      &average(hist2_vcppee),stderr(hist2_vcppee)
    else
     write(6,4)'Core polarization        (',trim(e_units),') : ', &
      &average(hist2_vcppei),stderr(hist2_vcppei)
    endif
   endif
   if(nbasis>1)then
    if(got_eionion)then
     write(6,5)'N-N repulsion             (',trim(e_units),') : ',eionion
    else
     write(6,*)'N-N repulsion: couldn''t find it'
    endif
   endif
   if(n>17)then ! relativistic components present
    write(6,*)
    write(6,*)'Relativistic corrections'
    write(6,*)
    write(6,4)'Mass polarization         (',trim(e_units),&
     &') : ',average(18),stderr(18)
    write(6,4)'Mass velocity term        (',trim(e_units),&
     &') : ',average(19),stderr(19)
    write(6,4)'Electron-nucleus Darwin   (',trim(e_units),&
     &') : ',average(20),stderr(20)
    write(6,4)'Electron-electron Darwin  (',trim(e_units),&
     &') : ',average(21),stderr(21)
    write(6,4)'Retardation term          (',trim(e_units),&
     &') : ',average(22),stderr(22)
    write(6,4)'Total rel. correction     (',trim(e_units),&
     &') : ',average(23),stderr(23)
   endif

  endif

4 format(t2,a,a,a,t38,f21.12,f21.12)
5 format(t2,a,a,a,t38,f21.12)

 END SUBROUTINE write_out_dmc


 CHARACTER(20) FUNCTION i2s(n)
!-----------------------------------------------------------------------!
! I2S                                                                   !
! ===                                                                   !
! Convert integers to left justified strings that can be printed in the !
! middle of a sentence without introducing large amounts of white space.!
!                                                                       !
! Calling routine is intended to include something like:                !
! USE utilities                                                         !
! INTEGER i                                                             !
! i=12                                                                  !
! write(6,*)'Integer number ',trim(i2s(i)),' with words at the end.'    !
!-----------------------------------------------------------------------!
  INTEGER i,j,n
  CHARACTER tmp,sign

  if(n==0)then
   i2s='0' ; return
  endif
  sign=' ' ; if(n<0)sign='-'

  do i=1,len(i2s)
   i2s(i:i)=' '
  enddo

  i=abs(n)
  do j=1,len(i2s)
   if(i==0)exit
   i2s(j:j)=achar(ichar('0')+mod(i,10))
   i=i/10
  enddo

  i=1 ; j=len_trim(i2s)
  do
   if(i>=j)exit
   tmp=i2s(j:j)
   i2s(j:j)=i2s(i:i)
   i2s(i:i)=tmp
   i=i+1
   j=j-1
  enddo

  i2s=trim(sign)//i2s

 END FUNCTION i2s


 CHARACTER(r2s_length) FUNCTION r2s(r,real_format)
!-------------------------------------------------------------------------!
! Converts real variable with arbitrary format to string that can be      !
! trimmed and printed in the middle of a sentence without introducing     !
! large amounts of white space, as you would if you did                   !
! write(6,'(f12.6)')12.0 or similar. Note you need to pass through the    !
! format string e.g. f12.6 .                                              !
!                                                                         !
! Calling routine is intended to include something like:                  !
! DOUBLE PRECISION r                                                      !
! r=12.d0                                                                 !
! tmpr=r2s(r,'(f12.6)')                                                   !
! write(6,*)'Real number ',trim(tmpr),' with words at the end.'           !
!                                                                         !
! Note : DON'T USE R2S IN A WRITE STATEMENT SINCE THIS IS ILLEGAL         !
! IN FORTRAN90 (ALTHOUGH NOT IN FORTRAN200X). IF ANYONE HAS TIME, FEEL    !
! FREE TO WRITE A VERSION OF THIS WHICH ISN'T ILLEGAL - SIMILAR TO        !
! I2S ABOVE - SO THAT PEOPLE WHO HAVEN'T READ THIS NOTE DON'T FEEL        !
! TEMPTED TO CALL R2S IN A WRITE STATEMENT.                               !
!-------------------------------------------------------------------------!

 IMPLICIT NONE
 DOUBLE PRECISION,INTENT(in) :: r
 CHARACTER(*),INTENT(in) :: real_format

 if(len(r2s)>0)then
  write(r2s,real_format)r
  r2s=adjustl(r2s)
 endif

 END FUNCTION r2s


END MODULE io_module


MODULE try_to_find_things
!----------------------------------------------------------------!
! Subroutines for determining the values of parameters using     !
! any CASINO input files etc that may be floating about.         !
!----------------------------------------------------------------!
 USE io_module, ONLY : i2s
 IMPLICIT NONE
 PRIVATE
 PUBLIC try_to_find_number_of_holes,try_to_find_nelec, &
  &try_to_find_dmc_algorithm,try_to_find_iterac,try_to_find_masses, &
  &try_to_find_eionion
 INTEGER,PARAMETER :: max_line_length=200
 INTEGER,PARAMETER :: nf_max=50 ! max # fields that may be pulled from a string
 INTEGER nf ! no. of fields pulled from string
 CHARACTER(len=max_line_length) :: fields(nf_max)


CONTAINS


 SUBROUTINE getline(io_unit,keystring,fields,nf,nf_max,found)
!----------------------------------------------------------------------!
! Scan forward through the file identified by io_unit to find the      !
! string stored in keystring.  Return the fields (items separated by   !
! spaces) from that line.  If keystring contains a single blank space  !
! then just take the next line from the file and split it into fields. !
!----------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: io_unit  ! unit no. of file to read
  INTEGER,INTENT(in) :: nf_max   ! max no. of fields in a string
  INTEGER,INTENT(out) :: nf      ! actual no. of fields in a string
  CHARACTER(*),INTENT(in) :: keystring ! the string to search for
  CHARACTER(*),INTENT(inout) :: fields(nf_max) ! fields pulled from string
  CHARACTER(130) instring
  INTEGER icheck
  LOGICAL search,found

  found=.true.
  if(keystring==' ')then
! Just get the next line from the file and parse it
   read(io_unit,fmt="(a130)",iostat=icheck) instring
   if(icheck < 0)then
    write(6,*)'GETLINE - next line is eof.'
    stop
   else
    call scan_string(instring,fields,nf,nf_max," ")
   endif
  else
! Search through each line of the file until the search string is found
   rewind(io_unit)
   search=.true.
   do while(search)
    read(io_unit,fmt="(a130)",iostat=icheck) instring
    if(icheck<0)then
     found=.false.
     return
    endif

! Look to see whether any of the fields contain the string we want
    if(index(string=instring,substring=keystring) /= 0)then
     search=.false.
! it does - split this line into fields and return them
     call scan_string(instring,fields,nf,nf_max," ")
     exit
    endif
   enddo
  endif

 END SUBROUTINE getline


 SUBROUTINE scan_string(instring,fields,nf,nf_max,separator)
!----------------------------------------------------------------------!
! Parse the string 'instring' and store the fields (objects separated  !
! by spaces) in the fields array.                                      !
!----------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(out) :: nf
  INTEGER,INTENT(in)  :: nf_max
  CHARACTER(*),INTENT(in) :: instring
  CHARACTER(*),INTENT(inout) :: fields(nf_max)
  CHARACTER(1),INTENT(in) :: separator ! character that separates fields
  INTEGER :: ilength ! length of scanned string (excluding trailing spaces)
  INTEGER :: icount  ! counter for no. of fields found
  INTEGER :: iposn   ! current character of string being scanned
  INTEGER :: ifpt    ! current character of field being assigned
  LOGICAL :: newfield! flag to say whether we are currently within a
                     ! field or between fields (i.e. the last character
                     ! scanned was a space)
! String into which to copy instring so as to avoid altering instring
! in the following manipulations
  CHARACTER(len=(len(instring))) :: tempstring
  tempstring=trim(adjustl(instring))
  ilength=len_trim(tempstring)
  iposn=0
  icount=0
  newfield=.true.
! initialise fields as the elements are constructed character by character
  fields=" "
  do while(iposn<ilength)
   iposn=iposn+1
   if(tempstring(iposn:iposn) /= separator)then
    if(newfield)then
     newfield=.false.
     icount=icount+1
     if(icount > nf_max)then
      write(6,fmt="('Only ',i2,' fields per line allowed. Increase &
       &nf_max in routine calling scan_string.f90.')") nf_max
      stop
     endif
     fields(icount)(1:1)=tempstring(iposn:iposn)
     ifpt=1
    else
     ifpt=ifpt+1
     fields(icount)(ifpt:ifpt)=tempstring(iposn:iposn)
    endif
   else
    if(.not.newfield) newfield=.true.
   endif
  enddo
! Return the no. of fields found
  nf=icount

 END SUBROUTINE scan_string


 SUBROUTINE try_to_find_number_of_holes
  USE qmcdata
  IMPLICIT NONE
  INTEGER ierr
  LOGICAL inp,found1,found2

  inquire(file='input',exist=inp)
  if(inp)then
   write(6,*)
   write(6,*)'Found input file.'
   open(9,file='input',status='old',iostat=ierr)
   if(ierr/=0)then
    write(6,*)'Couldn''t open it though.'
   else
    call getline(9,'nhu ',fields,nf,nf_max,found1)
    if(found1)then
     read(fields(3),*)nele(3)
    else
     write(6,*)'Couldn''t find number of up-spin holes. Assuming zero.'
     nele(3)=0
    endif
    call getline(9,'nhd ',fields,nf,nf_max,found2)
    if(found2)then
     read(fields(3),*)nele(4)
    else
     write(6,*)'Couldn''t find number of down-spin holes. Assuming zero.'
     nele(4)=0
    endif
    close(9)
    if(found1.or.found2)then
     write(6,*)'Extracted number of holes.'
     write(6,*)'Number of up-spin holes   : ',nele(3)
     write(6,*)'Number of down-spin holes : ',nele(4)
    endif
   endif
  endif

 END SUBROUTINE try_to_find_number_of_holes


 SUBROUTINE try_to_find_nelec
  USE qmcdata
  IMPLICIT NONE
  INTEGER ierr
  LOGICAL inp,found1,found2

  inquire(file='input',exist=inp)
  if(inp)then
   write(6,*)
   write(6,*)'Found input file.'
   open(9,file='input',status='old',iostat=ierr)
   call getline(9,'neu ',fields,nf,nf_max,found1)
   if(ierr/=0)then
    write(6,*)'Couldn''t open it though..'
   else
    if(found1)then
     read(fields(3),*)nele(1)
    else
     nele(1)=0 ! not really
    endif
    call getline(9,'ned ',fields,nf,nf_max,found2)
    if(found2)then
     read(fields(3),*)nele(2)
    else
     nele(2)=0 ! not really
    endif
    close(9)
    if(found1.or.found2)then
     write(6,*)'Extracted number of electrons.'
     write(6,'(1x,2a)')'Number of up-spin elecs   : ',trim(i2s(nele(1)))
     write(6,'(1x,2a)')'Number of down-spin elecs : ',trim(i2s(nele(2)))
    endif
    nelec=nele(1)+nele(2)
   endif
  endif

 END SUBROUTINE try_to_find_nelec


 SUBROUTINE try_to_find_dmc_algorithm
  USE qmcdata
  IMPLICIT NONE
  INTEGER ierr
  LOGICAL inp,found1,found2,ask_user

  ask_user=.false.

  inquire(file='input',exist=inp)
  if(inp)then
   write(6,*)
   write(6,*)'Found input file.'
   open(9,file='input',status='old',iostat=ierr)
   if(ierr==0)then
! Able to open input file.
    call getline(9,'dmc_algorithm ',fields,nf,nf_max,found1)
    if(found1)then
     read(fields(3),*)dmc_algorithm
    else
     call getline(9,'dmcalgorithm ',fields,nf,nf_max,found2)
     if(found2)then
      read(fields(3),*)dmc_algorithm
     else
      write(*,*)'No information about DMC algorithm in input.'
      write(6,*)'Assuming DMC algorithm is 4 (Latest version).'
      dmc_algorithm=4
     endif
    endif
    if(found1.or.found2)then
     write(6,*)'Extracted dmc_algorithm from input:'
     if(dmc_algorithm==1)then
      write(6,*)'Using value 1 = OLD DMC ALGORITHM to work out nequil.'
     elseif(dmc_algorithm==2)then
      write(6,*)'Using value 2 = NEW DMC ALGORITHM to work out nequil.'
     elseif(dmc_algorithm==3)then
      write(6,*)'Using value 3 = EVEN NEWER DMC ALGORITHM to work out nequil.'
     elseif(dmc_algorithm==4)then
      write(6,*)'Using value 4 = VERY LATEST DMC ALGORITHM to work out nequil.'
     else
      write(6,*)'DMC_ALGORITHM keyword in input file is not equal to 1, 2, 3 or&
       & 4.'
      stop
     endif
    endif
   else
    write(6,*)'Couldn''t open input file.'
    ask_user=.true.
   endif
   close(9)
  else
   write(6,*)'Couldn''t find input file.'
   ask_user=.true.
  endif

  if(ask_user)then
   dmc_algorithm=0
   do while(dmc_algorithm<1.or.dmc_algorithm>4)
    write(6,*)'Please enter DMC algorithm (1--4) (4 being latest algorithm).'
    if(reblock_in)then
     dmc_algorithm=4
     write(6,*)'Non-interactive mode : assuming dmc_algorithm=4'
    else
     read(5,*)dmc_algorithm
    endif
   enddo
  endif

 END SUBROUTINE try_to_find_dmc_algorithm


 SUBROUTINE try_to_find_iterac
  USE qmcdata
  IMPLICIT NONE
  INTEGER ierr
  LOGICAL inp,found1,got_iterac

  got_iterac=.false.
  inquire(file='input',exist=inp)
  if(inp)then
   open(9,file='input',status='old',iostat=ierr)
   if(ierr/=0)then
    write(6,*)'Couldn''t open input file, though it exists..'
   else
    call getline(9,'iterac ',fields,nf,nf_max,found1)
    if(found1)read(fields(3),*)iterac
    close(9)
    if(found1)then
     write(6,'(1x,2a)')'Extracted ITERAC parameter from input: ', &
      &trim(i2s(iterac))
     got_iterac=.true.
     if(iterac==3)then
      write(6,*)'==> Ewald interaction used in DMC propagation'
     elseif(iterac==4)then
      write(6,*)'==> MPC interaction used in DMC propagation'
     else
      write(6,*)'That''s not right. It should be 3 or 4 (?)'
     endif
    endif
   endif
  endif

  if(.not.got_iterac)then
   write(6,*)'Can''t read input file, so please tell me what value of the &
    &ITERAC parameter'
   write(6,*)'you used in this calculation:'
   write(6,*)'3 ==> Ewald interaction used in DMC propagation, or'
   write(6,*)'4 ==> MPC interaction used in DMC propagation?'
   write(6,*)
   do
    read(5,*)iterac
    if(iterac==3.or.iterac==4)then
     exit
    else
     write(6,*)'ITERAC must be 3 or 4 if Ewald and MPC present.'
     write(6,*)
    endif
   enddo
  endif

 END SUBROUTINE try_to_find_iterac


 SUBROUTINE try_to_find_masses
  USE qmcdata
  IMPLICIT NONE
  INTEGER ierr,i
  LOGICAL heg,got_masses

  got_masses=.false.
  inquire(file='heg.data',exist=heg)
  if(heg)then
   write(6,*)
   write(6,*)'Found heg.data file.'
   open(9,file='heg.data',status='old',iostat=ierr)
   if(ierr/=0)then
    write(6,*)'Couldn''t open it though.'
   else
    do i=1,23
     read(9,*)
    enddo
    read(9,*,iostat=ierr)mass_e,mass_h
    close(9)
    if(ierr/=0)then
     write(6,*)'Couldn''t find electron/hole masses though.'
    else
     write(6,*)'Extracted electron/hole masses.'
     got_masses=.true.
    endif
   endif
  endif
  if(.not.got_masses)then
   write(6,*)'Tell me the electron mass: '
   write(6,*)
   read(5,*)mass_e
   write(6,*)'And the hole mass: '
   write(6,*)
   read(5,*)mass_h
  endif

 END SUBROUTINE try_to_find_masses


 SUBROUTINE try_to_find_eionion(got_eionion)
  USE qmcdata
  IMPLICIT NONE
  LOGICAL,INTENT(out) :: got_eionion
  INTEGER ierr,i
  LOGICAL gwfn,pwfn,bwfn

  got_eionion=.false.

  inquire(file='gwfn.data',exist=gwfn)
  inquire(file='pwfn.data',exist=pwfn)
  inquire(file='bwfn.data',exist=bwfn)
  if(gwfn)then
   write(6,*)
   write(6,*)'Found gwfn.data file.'
   open(9,file='gwfn.data',status='old',iostat=ierr)
   if(ierr/=0)then
    write(6,*)'Couldn''t open it though..'
    return
   endif
   do i=1,15
    read(9,*)
   enddo
   read(9,*,err=1,end=1)eionion
   eionion=eionion*nbasis ! au/atom --> au per molecule / au per prim cell
   close(9)
   write(6,*)'Extracted N-N repulsion energy.'
   got_eionion=.true.
   return
1  write(6,*)'Couldn''t find N-N repulsion energy though..'
   return
  endif
  if(pwfn)then
   write(6,*)
   write(6,*)'Found pwfn.data file.'
   open(9,file='pwfn.data',status='old',iostat=ierr)
   if(ierr/=0)then
    write(6,*)'Couldn''t open it though..'
    return
   endif
   do i=1,27
    read(9,*)
   enddo
   read(9,*,err=2,end=2)eionion
   close(9)
   write(6,*)'Extracted N-N repulsion energy.'
   got_eionion=.true.
   return
2  write(6,*)'Couldn''t find N-N repulsion energy though..'
   return
  endif
  if(bwfn)then
   write(6,*)
   write(6,*)'Found bwfn.data file.'
   open(9,file='bwfn.data',status='old',iostat=ierr)
   if(ierr/=0)then
    write(6,*)'Couldn''t open it though..'
    return
   endif
   do i=1,27
    read(9,*)
   enddo
   read(9,*,err=3,end=3)eionion
   close(9)
   write(6,*)'Extracted N-N repulsion energy.'
   got_eionion=.true.
   return
3  write(6,*)'Couldn''t find N-N repulsion energy though..'
   return
  endif
  if(.not.gwfn.and..not.pwfn.and..not.bwfn)then
   write(6,*)
   write(6,*)'No xwfn.data file (x=g/p/b) to extract N-N repulsion from.'
  endif

 END SUBROUTINE try_to_find_eionion

END MODULE try_to_find_things


MODULE stats_analysis
!-----------------------------------------------------------------!
! Subroutines for performing analysis of VMC and DMC .hist files. !
!-----------------------------------------------------------------!
 IMPLICIT NONE
 PRIVATE
 PUBLIC analyse_vmc,analyse_dmc

CONTAINS


 SUBROUTINE analyse_vmc(e_units)
!------------------------------------------------!
! Perform reblocking analysis on the VMC data.   !
!------------------------------------------------!
  USE qmcdata
  USE stats_calcs
  USE vmc_records
  USE io_module,ONLY : i2s,write_out_vmc
  IMPLICIT NONE
  CHARACTER(9),INTENT(in) :: e_units
  INTEGER i,j,jmax,nave,blocklength,no_blocks,ib,&
   &interaction,blockstart,blockstop,ialloc
  INTEGER,PARAMETER :: io_out=23
  DOUBLE PRECISION avetot,vartot,kurt,skew,minval,maxval, &
   &stderr,stderr2,var,last_block_wt,eff_no_blocks,corr_tau,corr_tau_err
  DOUBLE PRECISION,ALLOCATABLE :: block_energy(:),ave(:),err(:)
  LOGICAL isperiodic

! Can only use weights for analysing DMC data.
  have_total_weights=.false.

! Decide whether total energy is for Ewald or MPC (or 1/r) interaction.
  if(items(hist_ewald,1)==0.d0.or.items(hist_mpc,1)==0.d0)then
   if(items(hist_ewald,1)==0.d0.and.items(hist_mpc,1)/=0.d0)then
    hist_etot=hist_mpc
    interaction=2
   endif
   if(items(hist_mpc,1)==0.d0.and.items(hist_ewald,1)/=0.d0)then
    hist_etot=hist_ewald
    interaction=1
   endif
   if(items(hist_mpc,1)==0.d0.and.items(hist_ewald,1)==0.d0)then
    write(6,*)'Total energies with both MPC and Ewald = 0.d0. Problem.'
    stop
   endif
  else
   write(6,*)'Both Ewald and MPC energies present.'
   do
    write(6,*)'Select interaction for use in statistics analysis:'
    write(6,*)'(1) Ewald'
    write(6,*)'(2) MPC'
    write(6,*)
    if(reblock_in)then
     interaction=reblock_in_iterac
     write(6,*)'Using default from reblock.in : ',trim(i2s(interaction))
    else
     read(5,*)interaction
    endif
    if(interaction==1.or.interaction==2)exit
   enddo
   if(interaction==1)hist_etot=hist_ewald
   if(interaction==2)hist_etot=hist_mpc
  endif

! Basic statistical analysis of data from 1-->nlines.
! nave=number of lines to average over (entire data range
! in case of VMC).
  nave=nlines

! Compute lots of information about the total energy data.
  call compute_stats_unweighted(nave,items(hist_etot,1:nave),avetot, &
   &vartot,skew,kurt,maxval,minval)
  call correlation_time(items(hist_etot,1:nave),corr_tau,corr_tau_err)

  write(6,*)
  write(6,*)'Analysis of total energy'
  write(6,*)'------------------------'
  write(6,*)'Minimum energy (',trim(e_units),') = ',minval
  write(6,*)'   Mean energy (',trim(e_units),') = ',avetot
  write(6,*)'Maximum energy (',trim(e_units),') = ',maxval
  write(6,*)'Std. Dev.   s2 (',trim(e_units),') = ',sqrt(vartot)
  write(6,*)'Skewness    s3 (',trim(e_units),') = ',skew
  write(6,*)'Kurtosis    s4 (',trim(e_units),') = ',kurt
  write(6,*)
  write(6,*)'Correlation time            (steps) = ',corr_tau
  write(6,*)'Error in correlation time       +/- = ',corr_tau_err

! Analyse data by thirds - catch very weird unequilibrated runs etc..
  call an_by_thirds(nave,hist_etot,1)

! Reblock data to estimate standard error of mean

! First work out how many block transformations are possible
  do j=1,nave
   if(2**j>=nave)then
    jmax=j-1
    exit
   endif
  enddo

! Array to hold block averages
  allocate(block_energy(nave),stat=ialloc)
  if(ialloc/=0)then
   write(6,*)'Sorry, can''t allocate the block energy array.'
   stop
  endif

  write(6,*)
  write(6,*)'Reblocking of total energy'
  write(6,*)'--------------------------'
  write(6,*)'Performing analysis on set of records from nstart to nstop'
  write(6,*)'which will accomodate exactly blocks with lengths corresponding'
  write(6,*)'to powers of 2 ranging from 0 to maxpower, where'
  write(6,'(1x,7a)')'nstart=','1',',  nstop=',trim(i2s(nave)),' and maxpower=',&
   &trim(i2s(jmax-1)),'.'
  write(6,*)
  write(6,*)'#Bt  #Blocks  Length  Var block means       Std err of mean&
   &    Err in seom'

! File for writing reblocking info for plotting
  open(io_out,file='reblock.plot',status='unknown')
  write(io_out,fmt="('#BT  Std err of mean   Err in seom')")

! Loop over reblock transformation number.
  do j=0,jmax-1

! Block length for this reblock transformation.
   blocklength=2**j

! Number of such blocks (including partially completed block).
   if(mod(nave,blocklength)==0)then
    no_blocks=nave/blocklength
   else
    no_blocks=nave/blocklength+1
   endif

! Number of blocks, including last block as a fraction.
   eff_no_blocks=dble(nave)/dble(blocklength)

! Compute block averages for all but the last block. Weights are 1.
   do ib=1,no_blocks-1
    blockstart=(ib-1)*blocklength+1
    blockstop=ib*blocklength
    call compute_av_var_unweighted(blocklength, &
     &items(hist_etot,blockstart:blockstop),.false.,block_energy(ib),var)
   enddo

! Compute block average for last (possibly incomplete) block.
! Weight is given by ratio of its length to that of the other blocks.
   blockstart=(no_blocks-1)*blocklength+1
   blockstop=nave
   call compute_av_var_unweighted(blockstop-blockstart+1, &
    &items(hist_etot,blockstart:blockstop),.false.,block_energy(no_blocks),var)
   last_block_wt=dble(blockstop-blockstart+1)/dble(blocklength)

! Compute average and variance of block averages.
   call compute_av_var_one_wt(no_blocks,block_energy(1:no_blocks), &
    &last_block_wt,.true.,avetot,vartot)
   stderr=sqrt(vartot)/sqrt(eff_no_blocks)
   stderr2=stderr/sqrt(2.d0*(eff_no_blocks-1.d0))

! And write them out.
   write(6,7)j,no_blocks,blocklength,vartot,stderr,stderr2
7  format(i3,1x,i7,1x,i7,1x,f18.10,1x,f18.10,1x,f18.10)
   write(io_out,*)j,stderr,stderr2

  enddo

  write(6,*)
  write(6,*)'Reblocking info written to reblock.plot file.'
  write(6,*)'(Standard error against reblock transformation number.)'
  write(6,*)'View with plot_reblock utility.'

  close(io_out)

! Now choose block length for computing averages of all quantities.
  do
   write(6,*)
   write(6,*)'Choose block length for computing final statistics: (0 to exit)'
   if(reblock_in)then
    blocklength=reblock_in_blocklength
    write(6,*)'Using default from reblock.in : ',trim(i2s(blocklength))
   else
    read(5,*)blocklength
   endif
   if(blocklength==0)stop
   if(blocklength>2**(jmax-1))then
    write(6,*)
    write(6,*)'That''s too long. Try again.'
   elseif(blocklength<1)then
    write(6,*)
    write(6,*)'Oh for God''s sake.. Pay attention.'
   else
    exit
   endif
  enddo

! Corresponding number of blocks.
  if(mod(nave,blocklength)==0)then
   no_blocks=nave/blocklength
  else
   no_blocks=nave/blocklength+1
  endif

! Number of blocks, including last block as a fraction.
  eff_no_blocks=dble(nave)/dble(blocklength)

 ! Compute stats for all components of energy with requested block length

  allocate(ave(nitems),err(nitems),stat=ialloc)
  if(ialloc/=0)then
   write(6,*)'Sorry, can''t allocate arrays for averages and std. errors.'
   stop
  endif

  do i=1,nitems

! Compute block averages for all but the last block. Weights are 1.
   do ib=1,no_blocks-1
    blockstart=(ib-1)*blocklength+1
    blockstop=ib*blocklength
    call compute_av_var_unweighted(blocklength, &
     &items(i,blockstart:blockstop),.false.,block_energy(ib),var)
   enddo

! Compute block average for last (possibly incomplete) block.
! Weight is given by ratio of its length to that of the other blocks.
   blockstart=(no_blocks-1)*blocklength+1
   blockstop=nave
   call compute_av_var_unweighted(blockstop-blockstart+1, &
    &items(i,blockstart:blockstop),.false.,block_energy(no_blocks),var)
   last_block_wt=dble(blockstop-blockstart+1)/dble(blocklength)

! Evaluate average and variance of block averages.
   call compute_av_var_one_wt(no_blocks,block_energy(1:no_blocks), &
    &last_block_wt,.true.,avetot,vartot)

   ave(i)=avetot
   err(i)=sqrt(vartot)/sqrt(eff_no_blocks)

  enddo

  if(periodicity==0)then
   isperiodic=.false.
  else
   isperiodic=.true.
  endif

  call write_out_vmc(nitems,no_blocks,blocklength,ave,err,isperiodic,e_units)

  close(io_out)
  deallocate(block_energy,ave,err)

 END SUBROUTINE analyse_vmc


 SUBROUTINE analyse_dmc(e_units,got_eionion)
!------------------------------------------------!
! Perform reblocking analysis on the DMC data.   !
!------------------------------------------------!
  USE dmc_records
  USE qmcdata
  USE stats_calcs
  USE try_to_find_things
  USE io_module, ONLY : write_out_dmc,i2s
  IMPLICIT NONE
  CHARACTER(9),INTENT(in) :: e_units
  LOGICAL,INTENT(in) :: got_eionion
  INTEGER,PARAMETER :: io_out=23
  INTEGER i,j,jmax,nequil,nstart,nstop,nave,&
   &blocklength,no_blocks,ib,ialloc,blockstart,blockstop
  DOUBLE PRECISION avetot,vartot,averef,varref,kurt,skew,&
   &minval,maxval,stderr,stderr2,ave_e_other_coul,stderr_e_other_coul,temp, &
   &eff_no_blocks,corr_tau,corr_tau_err
  DOUBLE PRECISION,ALLOCATABLE :: block_energy(:),ave(:),err(:), &
   &reblocked_av_etot(:),reblocked_stderr_etot(:),block_weight(:)
  LOGICAL isperiodic,have_dteff_info

! Sort out Ewald/MPC stuff
  if(dmc2)then
   if(items2(hist2_pote,1)==0.d0.or.items2(hist2_smpc,1)==0.d0)then
    if(items2(hist2_pote,1)==0.d0.and.items2(hist2_smpc,1)/=0.d0)iterac=2
    if(items2(hist2_smpc,1)==0.d0.and.items2(hist2_pote,1)/=0.d0)iterac=1
    if(nelec>1)then
     if(items2(hist2_smpc,1)==0.d0.and.items2(hist2_pote,1)==0.d0)then
      write(6,*)'The e-e energies with both MPC and Ewald = 0.d0 in dmc.hist2 &
       &file. Problem.'
      stop
     endif
    endif
   else
    write(6,*)
    write(6,*)'Both Ewald and MPC energies present in dmc.hist2 file.'
    call try_to_find_iterac
    write(6,*)
   endif
  endif

  if(nitems>=13)then
   write(6,*)
   write(6,*)'Total weights at each iteration will be used to weight'
   write(6,*)'the energies in the calculation of averages.'
   write(6,*)
   have_total_weights=.true.
  else
   write(6,*)
   write(6,*)'Total weights at each iteration are unavailable.'
   write(6,*)'Will simply take unweighted averages of energies.'
   write(6,*)
   have_total_weights=.false.
  endif

  if(dmc2.and.nitems2>=17)then
   have_dteff_info=.true.
  else
   have_dteff_info=.false.
  endif

! Calculate number of equilibration steps.. (can be tricky)
  call compute_dmc_nequil(nequil,6)

  do
   write(6,*)'Compute stats from which line? ( <0 to count back from end ).'
   if(reblock_in)then
    nstart=reblock_in_nstart
    write(6,*)'Using default from reblock.in file : ',trim(i2s(nstart))
   else
    read(5,*)nstart
   endif
   if(nstart==0)then
    write(6,*)'Line zero does not exist.'
   elseif(abs(nstart)>nlines)then
    write(6,*)'That''s outside the data range.'
   else
    exit
   endif
   write(6,*)
  enddo
  if(nstart<0)nstart=nlines+1+nstart

  if(nstart<=nequil)then
   write(6,*)'WARNING : Equilibration included in averaged data.'
   write(6,*)'(Statistics accumulation starts at line ',nequil+1,'.)'
   write(6,*)
  endif

  write(6,*)'Average over how many lines ? (''0'' to include all &
   &remaining data.)'
  if(reblock_in)then
   nave=reblock_in_nave
   write(6,*)'Using default from reblock.in file : ',trim(i2s(nave))
  else
   read(5,*)nave
  endif
  if(nave>0)then
   nstop=nstart+nave
   if(nstop>nlines)then
    nstop=nlines
    nave=nstop-nstart+1
    write(6,*)'Too many lines! Using',nave,'lines instead.'
   endif
  else
   nave=nlines-nstart+1
   nstop=nlines
  endif

  write(6,'(7a)')'Averaging ',trim(i2s(nave)),' lines (',trim(i2s(nstart)), &
   &' -> ',trim(i2s(nstop)),')'

! Basic statistical analysis of data from nstart -> nstop

  call compute_stats_unweighted(nave,items(hist_etot,nstart:nstop), &
   &avetot,vartot,skew,kurt,maxval,minval)
  call correlation_time(items(hist_etot,nstart:nstop),corr_tau,corr_tau_err)
  write(6,*)
  write(6,*)'Analysis of total energy'
  write(6,*)'------------------------'
  write(6,*)'Minimum energy (',trim(e_units),') = ',minval
  write(6,*)'   Mean energy (',trim(e_units),') = ',avetot
  write(6,*)'Maximum energy (',trim(e_units),') = ',maxval
  write(6,*)'Std. Dev.   s2 (',trim(e_units),') = ',sqrt(vartot)
  write(6,*)'Skewness    s3 (',trim(e_units),') = ',skew
  write(6,*)'Kurtosis    s4 (',trim(e_units),') = ',kurt
  write(6,*)
  write(6,*)'Correlation time            (steps) = ',corr_tau
  write(6,*)'Error in correlation time       +/- = ',corr_tau_err

  write(6,*)
  write(6,*)'Best estimate of energy (mixed estimator)'
  write(6,*)'-----------------------------------------'
  write(6,*)'Best estimate of energy (',trim(e_units),') at end of run = ', &
   &items(hist_ebest,nstop)

  if(dmc_algorithm==4)then
   if(any(items(hist_egrowth,:)/=0.d0))then
    write(6,*)
    write(6,*)'Growth estimator of energy'
    write(6,*)'--------------------------'
    write(6,*)'Growth estimator of energy (',trim(e_units), &
     &') at end of run = ',items(hist_egrowth,nstop)
   endif
  endif

  call compute_av_var_unweighted(nave,items(hist_eref,nstart:nstop), &
   &.true.,averef,varref)
  write(6,*)
  write(6,*)'Reference energy'
  write(6,*)'----------------'
  write(6,*)'   Mean energy (',trim(e_units),') = ',averef
  write(6,*)'Std. dev.   s2 (',trim(e_units),') = ',sqrt(varref)

! Analyse data by thirds - catch failed, short and poorly equilibrated runs
  call an_by_thirds(nave,hist_etot,nstart)

! Reblock data to estimate standard error of mean

! First work out how many block transformations are possible
  do j=1,nave
   if(2**j>=nave)then
    jmax=j-1
    exit
   endif
  enddo

  allocate(block_energy(nave),block_weight(nave),reblocked_av_etot(0:jmax-1), &
   &reblocked_stderr_etot(0:jmax-1),stat=ialloc)
  if(ialloc/=0)then
   write(6,*)'Sorry, can''t allocate arrays for reblocking.'
   stop
  endif

  write (6,6)
6 format(/t2,'Reblocking of total energy',&
   &/'--------------------------')
  write(6,*)'Performing analysis on set of records from nstart to nstop'
  write(6,*)'which will accomodate exactly blocks with lengths corresponding'
  write(6,*)'to powers of 2 ranging from 0 to maxpower, where'
  write(6,'(1x,7a)')'nstart=',trim(i2s(nstart)),',  nstop=',trim(i2s(nstop)), &
   &' and maxpower=',trim(i2s(jmax-1)),'.'
  write(6,*)
  write(6,*)'#Bt  #Blocks  Length  Var block means       Std err of mean&
   &    Err in seom'

! File for writing reblocking info for plotting
  open(io_out,file='reblock.plot',status='unknown')
  write(io_out,fmt="('#BT  Std err of mean   Err in seom')")

! Loop over reblock transformation number.
  do j=0,jmax-1

! Block length for this reblock transformation.
   blocklength=2**j

! Number of such blocks (including partially completed block).
   if(mod(nave,blocklength)==0)then
    no_blocks=nave/blocklength
   else
    no_blocks=nave/blocklength+1
   endif

! Number of blocks, including last block as a fraction.
   eff_no_blocks=dble(nave)/dble(blocklength)

! Compute block average energies and weights for all but the last block.
   do ib=1,no_blocks-1
    blockstart=(ib-1)*blocklength+nstart
    blockstop=ib*blocklength+nstart-1
    if(have_total_weights)then
     call compute_av_var_weighted(blocklength, &
      &items(hist_etot,blockstart:blockstop), &
      &items(hist_totmult,blockstart:blockstop),.false.,block_energy(ib),temp)
     block_weight(ib)=sum(items(hist_totmult,blockstart:blockstop))
    else
     call compute_av_var_unweighted(blocklength, &
      &items(hist_etot,blockstart:blockstop),.false.,block_energy(ib),temp)
     block_weight(ib)=dble(blocklength)
    endif
   enddo

! Compute block average energy and weight for last (possibly incomplete) block.
   blockstart=(no_blocks-1)*blocklength+nstart
   blockstop=nstop
   if(have_total_weights)then
    call compute_av_var_weighted(blockstop-blockstart+1, &
     &items(hist_etot,blockstart:blockstop), &
     &items(hist_totmult,blockstart:blockstop), &
     &.false.,block_energy(no_blocks),temp)
    block_weight(no_blocks)=sum(items(hist_totmult,blockstart:blockstop))
   else
    call compute_av_var_unweighted(blockstop-blockstart+1, &
     &items(hist_etot,blockstart:blockstop),.false.,block_energy(no_blocks), &
     &temp)
    block_weight(no_blocks)=dble(blockstop-blockstart+1)
   endif

! Compute average and variance of block averages.
   call compute_av_var_weighted(no_blocks,block_energy(1:no_blocks), &
    &block_weight(1:no_blocks),.true.,avetot,vartot)

   stderr=sqrt(vartot)/sqrt(eff_no_blocks)
   stderr2=stderr/sqrt(2.d0*(eff_no_blocks-1.d0))

   reblocked_av_etot(j)=avetot
   reblocked_stderr_etot(j)=stderr

   write(6,7)j,no_blocks,blocklength,vartot,stderr,stderr2
7  format(i3,1x,i7,1x,i7,1x,f18.10,1x,f18.10,1x,f18.10)
   write(io_out,*)j,stderr,stderr2

  enddo

  write(6,*)
  write(6,*)'Reblocking info written to reblock.plot file.'
  write(6,*)'(Standard error against reblock transformation number.)'
  write(6,*)'View with plot_reblock utility.'

  close(io_out)

! Now choose block length for computing averages of all quantities.
  do
   write(6,*)
   write(6,*)'Choose block length for computing final statistics: (0 to exit)'
   if(reblock_in)then
    blocklength=reblock_in_blocklength
    write(6,*)'Using default from reblock.in file : ',trim(i2s(blocklength))
   else
    read(5,*)blocklength
   endif
   if(blocklength==0)stop
   if(blocklength>2**(jmax-1))then
    write(6,*)'That''s too long. Try again.'
   elseif(blocklength<1)then
    write(6,*)'Oh for God''s sake.. Pay attention.'
   else
    write(6,*)
    exit
   endif
  enddo

! Corresponding number of blocks.
  if(mod(nave,blocklength)==0)then
   no_blocks=nave/blocklength
  else
   no_blocks=nave/blocklength+1
  endif

! Number of blocks, including last block as a fraction.
  eff_no_blocks=dble(nave)/dble(blocklength)

 ! Corresponding reblocking transformation number
  j=nint(log(dble(blocklength))/log(2.d0))

  allocate(ave(nitems2),err(nitems2),stat=ialloc)
  if(ialloc/=0)then
   write(6,*)'Sorry, can''t allocate arrays for averages and std. errors.'
   stop
  endif

 ! Compute stats for all components of energy with requested block length
  if(dmc2)then

! Compute weights for each block; either sum of move weights in block
! (if available) or simply number of data points in block (otherwise).
   if(have_total_weights)then
    do ib=1,no_blocks-1
     blockstart=(ib-1)*blocklength+nstart
     blockstop=ib*blocklength+nstart-1
     block_weight(ib)=sum(items(hist_totmult,blockstart:blockstop))
    enddo
    blockstart=(no_blocks-1)*blocklength+nstart
    blockstop=nstop
    block_weight(no_blocks)=sum(items(hist_totmult,blockstart:blockstop))
   else
    do ib=1,no_blocks-1
     block_weight(ib)=dble(blocklength)
    enddo
    blockstart=(no_blocks-1)*blocklength+nstart
    blockstop=nstop
    block_weight(no_blocks)=dble(blockstop-blockstart+1)
   endif

   do i=1,nitems2 ! Want items2(1-11 and 17)

    if(i<=hist2_vcppee.or.i>=hist2_dteff)then
! Compute block averages for all but the last block.
     do ib=1,no_blocks-1
      blockstart=(ib-1)*blocklength+nstart
      blockstop=ib*blocklength+nstart-1
      if(have_total_weights)then
       call compute_av_var_weighted(blocklength, &
        &items2(i,blockstart:blockstop), &
        &items(hist_totmult,blockstart:blockstop),.false.,block_energy(ib),temp)
      else
       call compute_av_var_unweighted(blocklength, &
        &items2(i,blockstart:blockstop),.false.,block_energy(ib),temp)
      endif
     enddo

! Compute block average for last (possibly incomplete) block.
     blockstart=(no_blocks-1)*blocklength+nstart
     blockstop=nstop
     if(have_total_weights)then
      call compute_av_var_weighted(blockstop-blockstart+1, &
       &items2(i,blockstart:blockstop), &
       &items(hist_totmult,blockstart:blockstop),.false., &
       &block_energy(no_blocks),temp)
     else
      call compute_av_var_unweighted(blockstop-blockstart+1, &
       &items2(i,blockstart:blockstop),.false.,block_energy(no_blocks),temp)
     endif

! Evaluate average and variance of block averages.
     call compute_av_var_weighted(no_blocks,block_energy(1:no_blocks), &
      &block_weight(1:no_blocks),.true.,avetot,vartot)

     ave(i)=avetot
     err(i)=sqrt(vartot)/sqrt(eff_no_blocks)

    endif ! i==label of required item.

   enddo

! If both MPC and Ewald were accumulated, figure out the total energy
! with the component that was not used in DMC propagation and therefore
! accumulated in the dmc.hist file.

   ave_e_other_coul=0.d0 ; stderr_e_other_coul=0.d0
   if(iterac==3.or.iterac==4)then

    select case (iterac)
    case(3) ! Ewald done - compute total energy with MPC interaction
     items(hist_etot,:)=items(hist_etot,:)-items2(hist2_pote,:)+&
      &items2(hist2_smpc,:)+items2(hist2_lmpc,:)
    case(4) ! MPC done - compute total energy with Ewald interaction
     items(hist_etot,:)=items(hist_etot,:)+items2(hist2_pote,:)-&
      &items2(hist2_smpc,:)-items2(hist2_lmpc,:)
    end select

! Compute block averages for all but the last block.
    do ib=1,no_blocks-1
     blockstart=(ib-1)*blocklength+nstart
     blockstop=ib*blocklength+nstart-1
     call compute_av_var_weighted(blocklength, &
      &items(hist_etot,blockstart:blockstop), &
      &items(hist_totmult,blockstart:blockstop),.false.,block_energy(ib),temp)
    enddo

! Compute block average for last (possibly incomplete) block.
    blockstart=(no_blocks-1)*blocklength+nstart
    blockstop=nstop
    call compute_av_var_weighted(blockstop-blockstart+1, &
     &items(hist_etot,blockstart:blockstop), &
     &items(hist_totmult,blockstart:blockstop), &
     &.false.,block_energy(no_blocks),temp)

! Evaluate average and variance of block averages.
    call compute_av_var_weighted(no_blocks,block_energy(1:no_blocks), &
     &block_weight(1:no_blocks),.true.,avetot,vartot)

    ave_e_other_coul=avetot
    stderr_e_other_coul=sqrt(vartot)/sqrt(eff_no_blocks)

   endif

  endif ! dmc2

  if(periodicity==0)then
   isperiodic=.false.
  else
   isperiodic=.true.
  endif

  call write_out_dmc(nitems2,no_blocks,blocklength,reblocked_av_etot(j), &
   &reblocked_stderr_etot(j),ave_e_other_coul,stderr_e_other_coul,ave,err, &
   &isperiodic,iterac,got_eionion,e_units)

  deallocate(reblocked_av_etot,reblocked_stderr_etot)

! Analyse run related data (populations, efficiency)

  call compute_stats_unweighted(nave,items(hist_pop,nstart:nstop),avetot, &
   &vartot,skew,kurt,maxval,minval)

  write(6,*)
  write(6,*)'Analysis of population'
  write(6,*)'----------------------'
  write(6,*)' Minimum pop : ',minval
  write(6,*)'    Mean pop : ',avetot
  write(6,*)' Maximum pop : ',maxval
  write(6,*)'Std. Dev. s2 : ',sqrt(vartot)
  write(6,*)'Kurtosis  s3 : ',kurt
  write(6,*)'Skewness  s4 : ',skew
  if(avetot-minval>0.25d0*avetot.or.maxval-avetot>0.25d0*avetot)write(6,8)
8 format('WARNING : Population fluctuated by more than 25% of mean.')

  call compute_stats_unweighted(nave,items(hist_acceptr,nstart:nstop),avetot, &
   &vartot,skew,kurt,maxval,minval)

  write (6,*)
  write (6,*)'Analysis of acceptance ratio'
  write (6,*)'----------------------------'
  write (6,*)' Minimum acc : ',minval
  write (6,*)'    Mean acc : ',avetot
  write (6,*)' Maximum acc : ',maxval

  if(have_dteff_info)then
   call compute_stats_unweighted(nave,items2(hist2_dteff,nstart:nstop), &
    &avetot,vartot,skew,kurt,maxval,minval)
   write(6,*)
   write(6,*)'Analysis of effective timestep'
   write(6,*)'------------------------------'
   write(6,*)'   Minimum value : ',minval
   write(6,*)'      Mean value : ',avetot
   write(6,*)' Mixed estimator : ',ave(hist2_dteff),' +/-',err(hist2_dteff)
   write(6,*)'     Maximum acc : ',maxval
  endif

  deallocate(block_energy,block_weight,ave,err)

 END SUBROUTINE analyse_dmc


 SUBROUTINE compute_dmc_nequil(nequil,io)
!------------------------------------------------------------------!
! Work out number of DMC equilibration moves.                      !
! An easier way of doing this (for algorithm 4 at least) would be  !
! to look at where the DMC move number is reset in dmc.hist.       !
!------------------------------------------------------------------!
  USE dmc_records
  USE qmcdata
  USE io_module, ONLY : i2s
  IMPLICIT NONE
  INTEGER,INTENT(in) :: io
  INTEGER,INTENT(out) :: nequil
  INTEGER i,j,k,n,blocksize,numblocks,bl_eq
  DOUBLE PRECISION ebest,ebest0,ebest1

! Compute number of DMC equilibration steps in data

  if(dmc_algorithm==1)then ! OLD ALGORITHM

! (etrial/=eref during equilibration - usually!!)
   do i=nlines,1,-1
    if(items(hist_eref,i)/=items(hist_ebest,i))exit
   enddo
   nequil=i

   write(io,*)
   write(io,'(1x,2a)')trim(i2s(nlines)),' lines of data in total:'
   write(io,'(1x,2a)')trim(i2s(nequil)),' lines of equilibration (probably), &
    &followed by'
   write(io,'(1x,2a)')trim(i2s(nlines-nequil)),' lines of accumulation.'

  elseif(dmc_algorithm==2.or.dmc_algorithm==3)then ! NEW ALGORITHM (MDT 3.2001)

! Only possible test is different no of moves in blocks for equil and accum(?)
! In any given block ebest is constant.

   ebest=items(hist_ebest,nlines)
   blocksize=1
   do i=nlines-1,1,-1
    if(ebest/=items(hist_ebest,i))exit
    blocksize=blocksize+1
   enddo

   numblocks=nlines/blocksize

   n=nlines
s: do i=1,numblocks
    ebest=items(hist_ebest,n)
    k=1
    do j=n-1,n-blocksize+1,-1
     if(items(hist_ebest,j)/=ebest)then
      nequil=j+k
      exit s
     endif
     k=k+1
    enddo
    n=n-blocksize
    nequil=nlines
   enddo s

   if(nequil==nlines)then
    write(io,*)
    write(io,*)'Couldn''t figure out number of equilibration lines in the DMC &
     &file - all block sizes the same.'
    write(io,*)
   else
    write(io,*)
    write(io,'(1x,2a)')trim(i2s(nlines)),' lines of data in total:'
    write(io,'(1x,2a)')trim(i2s(nequil)),' lines of equilibration (probably).'
    write(io,'(1x,2a)')trim(i2s(nlines-nequil)),' lines of accumulation.'
   endif

  elseif(dmc_algorithm==4)then

! Calculate block length in equilibration.
   ebest0=items(hist_ebest,1)
   ebest1=ebest0
   bl_eq=1
   do
    if(bl_eq+1>nlines)exit
    ebest1=items(hist_ebest,bl_eq+1)
    if(ebest0/=ebest1)exit
    bl_eq=bl_eq+1
   enddo

   if(bl_eq==1)then
    write(io,*)
    write(io,*)'Block length in equilibration appears to be one. => Can''t &
     &tell when statistics accumulation starts.'
    write(io,*)
    nequil=nlines
   else
! Look for where best estimate of energy is updated every move -- start of
! stats accumulation.
    i=0
    do
     if(i*bl_eq+2>nlines)then
      nequil=nlines
      exit
     endif
     if(items(hist_ebest,i*bl_eq+1)/=items(hist_ebest,i*bl_eq+2))then
      nequil=i*bl_eq
      exit
     endif
     i=i+1
    enddo

    write(io,*)
    write(io,'(1x,2a)')trim(i2s(nlines)),' lines of data in total:'
    write(io,'(1x,2a)')trim(i2s(nequil)),' lines of equilibration (probably).'
    write(io,'(1x,2a)')trim(i2s(nlines-nequil)),' lines of accumulation.'
   endif

  else
   write(io,*)'DMC algorithm confusion.'
   stop
  endif

 END SUBROUTINE compute_dmc_nequil


 SUBROUTINE an_by_thirds(nave,hist_etot,nstart)
!-----------------------------------------------------------------------!
! Divide data set into three and compare stats calculated over these    !
! three subsets.                                                        !
!-----------------------------------------------------------------------!
  USE qmcdata
  USE stats_calcs
  USE dmc_records, ONLY : hist_totmult
  IMPLICIT NONE
  INTEGER,INTENT(in) :: nave,hist_etot,nstart
  INTEGER nthird,i,nthirdstart
  DOUBLE PRECISION avet(3),vart(3),varsum

  write(6,3)
  if(nave<3)then
   write(6,*)'Not enough data for analysis by thirds.'
   stop
  endif
  nthird=nave/3
  do i=1,3
   nthirdstart=(i-1)*nthird+nstart
   if(have_total_weights)then
    call compute_av_var_weighted(nthird, &
     &items(hist_etot,nthirdstart:nthirdstart+nthird-1), &
     &items(hist_totmult,nthirdstart:nthirdstart+nthird-1),.true.,avet(i), &
     &vart(i))
   else
    call compute_av_var_unweighted(nthird, &
     &items(hist_etot,nthirdstart:nthirdstart+nthird-1),.true.,avet(i),vart(i))
   endif
   write(6,4)nthirdstart,nthirdstart+nthird-1,avet(i),vart(i),sqrt(vart(i))
  enddo
3 format(/t2,'Total energy analysed by thirds', &
   &/t2,'-------------------------------', &
   &/t2,'   From        To          Ave            Var          Std.dev')
4 format(1x,i8,'-',i8,1x,e14.6,1x,e14.6,1x,e14.6)
  varsum=sum(vart)
  do i=1,3
   if(vart(i)>2*(varsum-vart(i)))then
    select case(i)
     case(1)
     write(6,5)'1st'
    case(2)
     write(6,5)'2nd'
    case(3)
     write(6,5)'3rd'
    end select
5  format(/1x,'WARNING: Variance of ',a3,' third very large. Possible run &
    &failure.')
    exit
   endif
  enddo

 END SUBROUTINE an_by_thirds


 SUBROUTINE correlation_time(corr_ener,tau,err_tau)
!------------------------------------------------------------------------!
! Obtain correlation time from a set of data                             !
!------------------------------------------------------------------------!
 IMPLICIT NONE
 DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: corr_ener
 DOUBLE PRECISION,               INTENT(out):: tau,err_tau
 DOUBLE PRECISION autocorr,corr_aver,corr_var,delta_t,delta_0,tol
 INTEGER corr_ndata,i,j

 tol=1.d-15
 corr_ndata=size(corr_ener,1)

 if(corr_ndata<10)then
  tau=-1.d0
  err_tau=0.d0
  return
 endif

! Average

 corr_aver=0.d0
 do i=1,corr_ndata
  corr_aver=corr_aver+corr_ener(i)
 enddo
 corr_aver=corr_aver/dble(corr_ndata)

! Variance

 corr_var=0.d0
 do i=1,corr_ndata
  corr_var=corr_var+corr_ener(i)*corr_ener(i)
 enddo
 corr_var=corr_var/dble(corr_ndata)-corr_aver*corr_aver

 if(corr_var<tol)then
  tau=-1.d0
  err_tau=-1.d0
  return
 endif

! Autocorrelation for i <= cut-off. Tau

 tau=1.d0
 do i=1,corr_ndata-1
  autocorr=0.d0
  do j=1,corr_ndata-i
   delta_0=corr_ener(j)-corr_aver
   delta_t=corr_ener(i+j)-corr_aver
   autocorr=autocorr+delta_0*delta_t
  enddo
  autocorr=autocorr/(dble(corr_ndata-i)*corr_var)
  tau=tau+2.d0*autocorr
  if(i>=3.d0*tau)exit
 enddo

 err_tau=tau*sqrt((4.d0*dble(i)+2.d0)/dble(corr_ndata))

 END SUBROUTINE correlation_time


END MODULE stats_analysis



PROGRAM reblock
!----------------------------!
! Main program starts here.  !
!----------------------------!
 USE dmc_records
 USE io_module
 USE qmcdata
 USE try_to_find_things
 USE stats_analysis
 IMPLICIT NONE
 INTEGER ierr,nscale
 DOUBLE PRECISION escale,exscale,reduced
 CHARACTER(3) vmc_or_dmc
 CHARACTER(9) e_units
 LOGICAL got_eionion,inp

 write(6,*)'REBLOCK'
 write(6,*)'======='

 inquire(file='vmc.hist',exist=vmc)
 inquire(file='dmc.hist',exist=dmc)
 inquire(file='dmc.hist2',exist=dmc2)
 inquire(file='reblock.in',exist=reblock_in)
 inquire(file='input',exist=inp)

 if(reblock_in.and..not.inp)then
  write(6,*)'Need input file to work in non-interactive mode.'
  write(6,*)'Delete reblock.in file or find an input file.'
  write(6,*)
  stop
 endif

 do while(vmc.and.dmc)
  write(6,*)
  if(reblock_in)then
   write(6,*)'Both dmc.hist and vmc.hist exist in this directory. Quitting.'
   stop
  else
   write(6,*)'Both dmc.hist and vmc.hist exist in this directory.'
   write(6,*)'Which would you rather analyse? "VMC" or "DMC".'
  endif
  read(5,*)vmc_or_dmc
  if(vmc_or_dmc=='VMC'.or.vmc_or_dmc=='vmc')dmc=.false.
  if(vmc_or_dmc=='DMC'.or.vmc_or_dmc=='dmc')vmc=.false.
 enddo

 if(.not.vmc.and..not.dmc)then
  write(6,*)
  write(6,*)'REBLOCK cannot find a dmc.hist or a vmc.hist file to analyse.'
  write(6,*)'Quitting.'
  write(6,*)
  stop
 endif

 ierr=0
 if(vmc)then

  write(6,*)'Found vmc.hist file.'
  open(8,file='vmc.hist',status='old',iostat=ierr)
  if(ierr/=0)then
   write(6,*)'Error opening vmc.hist file.'
   stop
  endif
  write(6,*)'Reading data.'
  call read_vmc_data(8,6)
  close(8)

  if(reblock_in)then
   write(6,*)'Found reblock.in file.'
   open(8,file='reblock.in',status='old',iostat=ierr)
   if(ierr/=0)then
    write(6,*)'Error opening reblock.in file.'
    stop
   endif
   write(6,*)'Reading reblock.in file.'
   read(8,*,iostat=ierr)reblock_in_units,reblock_in_blocklength, &
    &reblock_in_iterac
   if(ierr/=0)then
    write(6,*)'Error reading reblock.in file.'
    stop
   endif
   if(reblock_in_units<1.or.reblock_in_units>5)then
    write(6,*)'Units in reblock.in out of range (1-5).'
    stop
   endif
   if(reblock_in_blocklength<1.or.reblock_in_blocklength>nlines)then
    write(6,*)'Block length invalid in reblock.in file.'
    stop
   endif
   if(reblock_in_iterac/=1.and.reblock_in_iterac/=2)then
    write(6,*)'Interaction type invalid in reblock.in file.'
    stop
   endif
   close(8)
  endif

 else

  write(6,*)'Found dmc.hist file.'
  open(8,file='dmc.hist',status='old',iostat=ierr)
  if(ierr/=0)then
   write(6,*)'Error opening dmc.hist file.'
   stop
  endif
  write(6,*)'Reading data.'
  call read_dmc_data(8,6)
  close(8)
  nbasis=int(items(hist_nbasis,1))
  periodicity=items(hist_periodicity,1)
  if(dmc2)then
   write(6,*)
   write(6,*)'Found dmc.hist2 file.'
   open(8,file='dmc.hist2',status='old',iostat=ierr)
   if(ierr/=0)then
    write(6,*)'Error opening dmc.hist2 file.'
    stop
   endif
   write(6,*)'Reading data.'
   call read_dmc2_data(8,6)
   if(nlines2/=nlines)then
    write(6,*)'The dmc.hist and dmc.hist2 files are required to have the &
     &same number of lines. They don''t.'
    stop
   endif
   close(8)
   eionion=0.d0
   if(nbasis>1)call try_to_find_eionion(got_eionion)
  endif
  if(nbasis==0)call try_to_find_number_of_holes
  call try_to_find_nelec
  call try_to_find_dmc_algorithm

  if(reblock_in)then
   write(6,*)'Found reblock.in file.'
   open(8,file='reblock.in',status='old',iostat=ierr)
   if(ierr/=0)then
    write(6,*)'Error opening reblock.in file.'
    stop
   endif
   write(6,*)'Reading reblock.in file.'
   read(8,*,iostat=ierr)reblock_in_units,reblock_in_nstart,reblock_in_nave,&
    &reblock_in_blocklength
   if(ierr/=0)then
    write(6,*)'Error reading reblock.in file.'
    stop
   endif
   if(reblock_in_units<1.or.reblock_in_units>5)then
    write(6,*)'Units in reblock.in out of range (1-5).'
    stop
   endif
   if(reblock_in_nstart<0.or.reblock_in_nstart>nlines)then
    write(6,*)'Start line invalid in reblock.in file.'
    stop
   endif
   if(reblock_in_nave<0.or.reblock_in_nave>(nlines-reblock_in_nstart))then
    write(6,*)'Number of lines to average over invalid in reblock.in file.'
    stop
   endif
   if(reblock_in_blocklength<1.or.reblock_in_blocklength>nlines)then
    write(6,*)'Block length invalid in reblock.in file.'
    stop
   endif
   close(8)
  endif


 endif ! vmc or dmc

! Assume input in au/prim.cell (or equivalently, au/particle for E- and EH-gas)
! Get scaling factors for conversion to other units.
 do
  write(6,*)
  write(6,*)'Your choice of output units:'
  if(nbasis==0)then               ! vmc and dmc for electron-(hole) gases
   write(6,*)'(1) Ha/electron'
   write(6,*)'(2) eV/electron'
   if(nele(3)>0.or.nele(4)>0)then ! i.e. if electron-hole gas
    write(6,*)'(3) Ex.u/exciton'
   endif
  endif
  if(nbasis/=0)then               ! vmc and dmc for ion-containing systems
   if(periodicity==0)then
    write(6,*)'(1) Ha/molecule'
    write(6,*)'(2) eV/molecule'
    write(6,*)'(3) Ha/atom'
    write(6,*)'(4) eV/atom'
    write(6,*)'(5) kcal/mol'
   else
    write(6,*)'(1) Ha/primitive cell'
    write(6,*)'(2) eV/primitive cell'
    write(6,*)'(3) Ha/atom'
    write(6,*)'(4) eV/atom'
   endif
  endif
  write(6,*)'?'
  write(6,*)
  if(reblock_in)then
   nscale=reblock_in_units
   write(6,*)'Using default from reblock.in file : ',trim(i2s(nscale))
   if(nbasis==0)then
    if(nele(3)>0.or.nele(4)>0)then ! electron-hole gas
     if(nscale/=1.and.nscale/=2.and.nscale/=3)then
      write(6,*)'Invalid for this system type.'
      stop
     endif
    else                           ! electron gas
     if(nscale/=1.and.nscale/=2)then
      write(6,*)'Invalid for this system type.'
      stop
     endif
    endif
   else
    if(periodicity==0)then
     if(nscale/=1.and.nscale/=2.and.nscale/=3.and.nscale/=4.and.nscale/=5)then
      write(6,*)'Invalid for this system type.'
      stop
     endif
    else
     if(nscale/=1.and.nscale/=2.and.nscale/=3.and.nscale/=4)then
      write(6,*)'Invalid for this system type.'
      stop
     endif
    endif
   endif
   exit
  else
   read(5,*)nscale
   if(nbasis==0)then
    if(nele(3)>0.or.nele(4)>0)then ! electron-hole gas
     if(nscale==1.or.nscale==2.or.nscale==3)exit
    else                           ! electron gas
     if(nscale==1.or.nscale==2)exit
    endif
   else
    if(periodicity==0)then
     if(nscale==1.or.nscale==2.or.nscale==3.or.nscale==4.or.nscale==5)exit
    else
     if(nscale==1.or.nscale==2.or.nscale==3.or.nscale==4)exit
    endif
   endif
  endif

 enddo

 select case (nscale)
  case(1)
   escale=1.d0
   e_units='Ha'
  case(2)
   escale=htoev
   e_units='eV'
  case(3)
   if(nbasis/=0)then ! something with ions in it
    escale=1.d0/dble(nbasis)
    e_units='Ha'
   else ! electron-hole gas (hopefully)
    if(vmc.and.(nele(3)==0.and.nele(4)==0))then
! (Note we don't have the nele array in DMC calcs)
     write(6,*)'Exciton units requested but no holes in electron-hole system?'
     stop
    endif
    if(dmc)then ! get the masses from the heg.data file, if possible.
     call try_to_find_masses
    endif
    write(6,*)
    write(6,*)'Mass(electron)    : ',mass_e
    write(6,*)'Mass(hole)        : ',mass_h
    reduced=mass_e*mass_h/(mass_e+mass_h)
    write(6,*)'Reduced mass      : ',reduced
    write(6,*)'Periodicity       : ',periodicity
    if(periodicity==3)then
     exscale=4.d0/reduced
    elseif(periodicity==2)then
     exscale=1.d0/reduced
    else
     write(6,*)'Only 2D and 3D electron-hole systems currently recognised.'
     stop
    endif
    escale=exscale
    write(6,*)'Scaling energy by : ',exscale
    write(6,*)
    e_units='Ex. units'
   endif
  case(4)
   escale=htoev/dble(nbasis)
   e_units='eV'
  case(5)
   escale=htokcal
   e_units='kcal/mol'
 end select

 if(vmc)then
  items=items*escale
  eionion=eionion*escale
  call analyse_vmc(e_units)
  deallocate(items)
 else ! dmc
  items(hist_etot:hist_ebest,:)=items(hist_etot:hist_ebest,:)*escale
  items(hist_egrowth,:)=items(hist_egrowth,:)*escale
  if(dmc2)then
   items2(hist2_ti:hist2_vcppee,:)=items2(hist2_ti:hist2_vcppee,:)*escale
   eionion=eionion*escale
  endif
  call analyse_dmc(e_units,got_eionion)
  deallocate(items)
  if(dmc2)deallocate(items2)
 endif

END PROGRAM reblock

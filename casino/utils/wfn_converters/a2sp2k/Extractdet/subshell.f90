module subshell_mod

 use anglib 
 use error_mod
 common/mt/mt(40)
 common/mt67/m76(238)
 common/skmt2/mtf(8),mt3(90)

contains

!----------------------------------------------------------------------
!---- NOTES....
!---- Do some error catching for overflowing array sizes.
!----             CURRENTLY - see error.f90 code
!---- Sort order of arrays! elements of A(i,j,k) are stored in order
!----  A(1,1,1),A(2,1,1), so make sure code takes this into account
!---- in nested loops and in access order.
!----             CURRENTLY - left it to the compiler.
!---- f-subshell? CURRENTLY - does up to 2(4?)-electrons in f-shell. 
!----                      could do more by replacing v with Nr 
!----                      (see G. Gaigalas paper on subroutine library)
!---- g-subshell? CURRENTLY - does up to 2(4?)-electrons. 
!----                      seniority is enough for these.
!---- ne=1        CURRENTLY - caught as special case at beginning
!---- ne=4*l+2    CURRENTLY - caught as special case at beginning
!---- Clumsy loops to determinant?
!----             CURRENTLY - still clumsy
!---- Need to derive the (-1)^l behaviour of the prefactor
!----   for closed shells. It appears to be correct for s,p,d
!----   but this is only by observation so far...
!----------------------------------------------------------------------
 subroutine subshell(l_shell,ne,l_tot,ml_tot,s_tot,ms_tot,       &
 &   v_tot,ml_subshell,ms_subshell,coef_subshell,num)
 implicit none
 external glcons,termls,trmf
 integer                      :: jthn
 integer                      :: maxic,max_targ,max_det
 parameter(maxic=100,max_targ=50,max_det=1024)
 integer,intent(in)           :: l_shell,ne
 integer,intent(in)           :: l_tot,s_tot,ml_tot,ms_tot,v_tot
 integer,intent(out)          :: num
 integer,intent(out)          :: ml_subshell(:,:)
 integer,intent(out)          :: ms_subshell(:,:)
 real(kind(1.d0)),intent(out) :: coef_subshell(:)
 integer            :: count  (maxic)
 real(kind(1.d0))   :: coef   (maxic,max_det)
 integer            :: dic    (maxic,max_det)
 integer            :: dl   (2,maxic,max_det)
 integer            :: ds   (2,maxic,max_det)
 integer            :: dv   (2,maxic,max_det)
 integer            :: dm   (2,maxic,max_det)
 integer            :: dms  (2,maxic,max_det)
 integer            :: dll  (0:maxic)
 integer            :: dss  (0:maxic)
 integer            :: dvv  (0:maxic)
 integer            :: dmm    (maxic,max_det)
 integer            :: dmsms  (maxic,max_det)
 integer            :: start(0:maxic)
 integer            :: stop (0:maxic)
 integer            :: ic,icc,n_up,n_down
 integer            :: l1,l2,l12,m1 ,m2 ,m12
 integer            :: s1,s2,s12,ms1,ms2,ms12
 integer            :: v1,v2,v12
 real(kind(1.d0))   :: s1_,s2_,s12_,ms1_,ms2_,ms12_
 character          :: spec(0:9),spec_lower(0:9)
 character(2)       :: dig(-9:9)
 integer            :: l_c(10),m_c(10),s_c(10),ms_c(10)
 integer            :: l(10),m(10),s(10),ms(10)
 integer            :: ml_targ(10,max_targ),ms_targ(10,max_targ)
 real(kind(1.d0))   :: coef_targ(max_targ)
 integer            :: itarg,n_targ
 real(kind(1.d0))   :: s_(10),ms_(10)
 integer            :: ic_c(11)
 integer            :: iflag
 integer            :: i,j,n,nn,i1,i2,n12,jb,jk,lll,nolist,kt,lq
 integer            :: index,num12,mxi,kkk,ne2
 integer            :: n10,n9,n8,n7,n6,n5,n4,n3,n2,n1
 real(kind(1.d0))   :: cfp,x,z_ls,fact

 data spec/'S','P','D','F','G','H','I','K','L','M'/
 data spec_lower/'s','p','d','f','g','h','i','k','l','m'/
 data dig/'-9','-8','-7','-6','-5','-4','-3','-2','-1',         &
 &         ' 0',' 1',' 2',' 3',' 4',' 5',' 6',' 7',' 8',' 9'/



!---- l_shell is l for the subshell being coupled
!---- ne electrons,
!---- ne-1 couplings
!---- for ne=1 main part of code does not do anything
!----  meaningfull so this is a special case.

 ml_subshell(:,:)=0.d0
 ms_subshell(:,:)=0.d0
 coef_subshell(:)=0.d0

 n_up  =(ne+ms_tot)/2
 n_down=(ne-ms_tot)/2



!---- Catch special cases so we can do them quickly, and do them for 
!----  l=g and above
!----  1 electron in subshell
 if (ne.eq.1) then
  num=1
  ml_subshell(1,1)=ml_tot
  ms_subshell(1,1)=ms_tot
  coef_subshell(1)=1.d0
  return
 endif

!----  Closed shell
 if (ne.eq.(4*l_shell+2)) then
  num=1
  do j=1,ne/2
   ml_subshell(j,1)     =l_shell-(j-1)
   ml_subshell(j+ne/2,1)=ml_subshell(j,1)
   ms_subshell(j,1)     = 1
   ms_subshell(j+ne/2,1)=-1
  enddo
!--- the phase in the next two expressions needs to be proven.
  coef_subshell(1)=(-1.0)**l_shell/sqrt(factorial(ne))    ! for order of decreasing m_l
!-coef_subshell(1)=  1.0          /sqrt(factorial(ne))    ! for order of increasing m_l
  return
 endif
     





!-----------------------------------------------------------------------------
!---- This section finds all determinants for this ne that give correct Ml,Ms.
!---- REMEMBER: INCORRECT L,M MAY BE INCLUDED IN THE LIST,
!---- SO MAY NOT CONTRIBUTE TO THE COUPLING LISTS
!---- also this list will not contain all of the TERMs, but
!---- only those that allow the Ml,Ms we have specified.
!---- Done for efficiency.

 nn =4*l_shell+2
 mxi=(nn**ne)-1

 ml_targ(:,:)=0;ms_targ(:,:)=0;coef_targ(:)=0.d0

 n_targ=1
 do i=0,2**nn-1
  call signature2det(i,ml_targ(:,n_targ),ms_targ(:,n_targ),            &
  &                     ne,n_up,n_down,l_shell,ne2)  

  iflag=1
!--- check that this is the correct number of electrons
  if (ne2.ne.ne) iflag=0
!--- check that this has correct Mz,Ml
  if (sum(ml_targ(1:ne,n_targ)).ne.ml_tot) iflag=0
  if (sum(ms_targ(1:ne,n_targ)).ne.ms_tot) iflag=0

  if (iflag.eq.1) then
   n_targ=n_targ+1
  endif
 enddo

 n_targ=n_targ-1
 if (n_targ.gt.max_targ) call error(6,n_targ,max_targ)

!do i=1,n_targ
!  write(*,'(i6,1x,50(i3,1x))') i,(ml_targ(j,i),j=1,ne),&
!  &                                  (ms_targ(j,i),j=1,ne)
!  write(*,*) i,det2signature(ml_targ(:,i),ms_targ(:,i),&
!  &                              ne,l_shell)
!enddo
!write(*,'(a,a1,a)') 'Doing an ',spec_lower(l_shell),' subshell'

!--- ml_targ(1:ne,1:n_targ), ms_targ(1:ne,1:n_targ) now contains
!--- representation of products that can make subshell determinants,
!--- and that have correct total Ml & Ms.
!--- These are ordered so spin up then spin down electrons,
!---  and decreasing Ml for spin up and spin down seperately.
!--- (used to make lists of coupling in next stage)







!-------------- STARTING BIG LOOP OVER TARGETS
 do itarg=1,n_targ
! zero some arrays
 count(:)=0
 dic(:,:)=0
 dl (:,:,:)=0
 ds (:,:,:)=0
 dv (:,:,:)=0
 dm (:,:,:)=0
 dms(:,:,:)=0
 dll  (:)=0
 dss  (:)=0
 dvv  (:)=0
 dmm  (:,:)=0
 dmsms(:,:)=0
 coef(:,:) =0.d0
 l_c (:)=0
 m_c (:)=0
 s_c (:)=0
 ms_c(:)=0





!-----------------------------------------------------------------------------
!--- This section creates the coupling list for coupling left-right.
!--- Couples only the Ml,Ms combinations on the list generated above,
!--- so does require us to create the sum of products that is the 
!--- full determinant, just the ordered product representation in the 
!--- target list. Note that if the target list does not 
!--- contain a state with the correct L M, an error will occur.
!--- This should not happen.


!--- set some data for 0 couplings, so j=1 in next section works correctly
 start(0)=0;stop(0)=0
 dll(0)=l_shell
 dss(0)=1
 dvv(0)=1


!---- variables for counting in coupling lists
 ic=0;icc=0


!----   1-2,12-3,123-4 ... 1..(ne-1)-ne  coupling
 do j=1,ne-1                                     ! j-loop (over couplings)
!write(*,*) 'doing ',j+1,'electron coupling'

 start(j)=ic+1
 do l12=0,(j+1)*l_shell                          ! l12-loop over children
 i1=mod(j+1,2);i2=i1+j+1
 do s12=i1,i2,2                                  ! s12-loop over children

!---- get number in the list that match l12,s12 - num12 of them
  n12=1;call get_index(index,n12,num12,l12,s12,l_shell)
!---- loop through states that match l12,s12 in list
!----  each has different seniority (and Nr)
  if (num12.gt.max_det) call error(7,num12,max_det)
  do n12=1,num12                                 ! n12 - loop over children
   ic=ic+1;n=0
   call get_index(index,n12,num12,l12,s12,l_shell)
   jb=index

   do icc=start(j-1),stop(j-1)                   ! icc-loop (over parents)

    l1=dll(icc);s1=dss(icc);v1=dvv(icc) ! parent state
    l2=l_shell;s2=1;v2=1                ! coupled to one l=l_shell electron

!------ get index in the list that matches l1,s1,v1 (parent state)
    call get_index_lsv(index,l1,s1,v1,l_shell)
    jk=index

    lll=j+1
    call mk_cfp(cfp,l_shell,lll,jb,jk)

!----- NOTE: For l_shell=3, code is only valid for ne<=4 
!-----       it will work for ne => 10 if ml_target arrays are bigger
!-----       and for 4<ne<10 if we use Nr instead of seniority, v, as index
!-----       For 4<l_shell<9 code will do ne<=2 only.
!-----       see Gaigalas papers on library, and source of library.
    if (abs(cfp).gt.1d-6) then
     if(l_shell.lt.3)then
       kt=mt(jb)
     elseif(l_shell.eq.3)then 
       kt=m76(jb)
     else
       kt=mt3(jb)
     endif
     lq=jthn(kt,3,100)
     v12=2*l_shell+2-(lq+1)

     call cppl(cfp,ic,n,l1,l2,l12,s1,s2,s12,count,coef,dl,ds,dm,   &
     &            dms,dmm,dmsms,v1,v2,dv,dic,icc,                    &
     &            lll,ml_targ,ms_targ,itarg)

!------- write out data on parent and child states
!    write(101,*) ic,count(ic),j+1
!    write(101,'(a7,1x,i2,a1,i1)') 'Parent:',s1+1 ,spec(l1),v1
!    write(101,'(a7,1x,i2,a1,i1)') 'Child :',s12+1,spec(l12),v12
!    write(101,'(a7,1x,f14.7)')    'cfp   :', cfp
!    write(101,'(a7,1x,f14.7)')    'coef  :', coef(ic,n)

     dll(ic)=l12
     dss(ic)=s12
     dvv(ic)=v12
!--------record which on the list the parent is

    endif

   enddo                                         ! icc-loop

   if (count(ic).eq.0) ic=ic-1

  enddo                                          ! n12-loop
 enddo                                           ! s12-loop
 enddo                                           ! l12-loop
 stop(j)=ic

 enddo                                           ! j-loop

 nolist=ic
 if (nolist.gt.maxic) call error(5,nolist,maxic)

!--- Now we have the coupling lists. This is contained in a number of arrays:
!--- nolist is the number of couplings that have been calculated.
!--- ic=stop(j-1),start(j-1) is the range of ic that result from
!--- coupling of j electrons, so ic is an index for all couplings.
!--- count(ic)           number of states that result from coupling ic
!---  dic(ic,n)          gives ic value of the coupling that produced
!---                       the many electron states coupled to a one electron
!---                       state in coupling ic. IMPORTANT!
!--- coupling is from coupling from a number of electrons to one electron,
!--- from left to right i=1,many electron state, i=2, single electrons state.
!---  dl (i,ic,n)        l  values of states coupled (for i=1,2)
!---  ds (i,ic,n)        s  values of states coupled (for i=1,2)
!---  dv (i,ic,n)        v  values of states coupled (for i=1,2)
!---  dm (i,ic,n)        m  values of states coupled (for i=1,2)
!---  dms(i,ic,n)        ms values of states coupled (for i=1,2)
!---  dll  (ic)          ll   value of resulting state
!---  dss  (ic)          ss   value of resulting state
!---  dvv  (ic)          vv   value of resulting state
!---  dmm  (ic,n)        mm   value of resulting state
!---  dmsms(ic,n)        msms value of resulting state
!---  coef (ic,n) =0.d0  coefficient associated with coupling, ie 
!---                     clebsch-gordon coef.* coef. of fractional parentage.






!--- List possible LS values for final subshell state
!write(*,*)
!write(*,*) 'Possible TERMS for final subshell: (2S+1)L(V)'
!do i=start(ne-1),stop(ne-1)
!  write(*,*) 'TERM: ',dig(dss(i)+1),spec(dll(i)),dig(dvv(i))(2:2),count(i)
!enddo
!write(*,*) stop(ne-1)-start(ne-1)+1,' of them'
!write(*,*)
!write(*,*) 'Doing itarg=',itarg
!write(201,*) 'Signature:',det2signature(ml_targ(:,itarg),      &
! &                                        ms_targ(:,itarg),      &
! &                                        ne,l_shell)
!write(*,'(a,i2,a1,i1)') 'Looking for :',s_tot+1,spec(l_tot),v_tot

 iflag=0
 do i=start(ne-1),stop(ne-1)
  if(l_tot.eq.dll(i).and.s_tot.eq.dss(i).and.v_tot.eq.dvv(i))then
   ic=i
   iflag=1
  endif
 enddo
 if (iflag.eq.0) then
  write(*,'(a,i2,a1,i1)') 'Could not find :',s_tot+1,spec(l_tot),v_tot
  cycle
 endif


!write(*,'(a,1x,i3,1x,f4.1)') 'Doing Ml,Ms,=',ml_tot,0.5d0*dble(ms_tot)
!write(*,'(a10,1x,i3)') 'Spin up  :',n_up
!write(*,'(a10,1x,i3)') 'Spin down:',n_down
!write(*,*)

      






!-----------------------------------------------------------------------------
!--- This section creates the sum of determinants (with coefficients)
!--- from the coupling lists calculated above. In other words it converts
!--- data on coupling 
!--- 1 elec. state to 1. elec. -> 2 electron states  (ic=start(1),stop(1))
!--- 2 elec. state to 1. elec. -> 3 electron states  (ic=start(2),stop(2))
!--- 3 elec. state to 1. elec. -> 4 electron states  (ic=start(3),stop(3))
!---              etc.
!--- That have available so to a sum of determinants (with coefficients).



 icc=ic                        ! icc is the ic index of the subshell-CSF

! all of the values defined below are physically meaningless
! these are set to ensure outer loops for electrons that do not exist
! run around once and don't do anything
! (eg by pointing ic_c(n) to maxic element of count for n>ne, 
!  which has value 1)
 ic_c(1:10)    =maxic
 count(maxic)  =1
 dic(maxic,1)  =maxic

 l_c (ne:10)   =l_tot
 m_c (ne:10)   =ml_tot
 s_c (ne:10)   =s_tot
 ms_c(ne:10)   =ms_tot

 dl (1,maxic,1)=l_tot
 ds (1,maxic,1)=s_tot
 dm (1,maxic,1)=ml_tot
 dms(1,maxic,1)=ms_tot

 dll(  maxic  )=l_tot
 dmm(  maxic,1)=ml_tot
 dss(  maxic  )=s_tot
 dmsms(maxic,1)=ms_tot






 


!-- left-right coupling, and there are ne-1 nested loops!

 if (ne.eq.10) then
   ic_c(10)=icc;l_c(10)=l_tot;s_c(10)=s_tot;m_c(10)=ml_tot;ms_c(10)=ms_tot
 endif
 do n10=1,count(ic_c(10))
 ic_c(9)=dic(  ic_c(10),n10)
 l_c (9)=dl (1,ic_c(10),n10);m_c (9)=dm (1,ic_c(10),n10)
 s_c (9)=ds (1,ic_c(10),n10);ms_c(9)=dms(1,ic_c(10),n10)
 l  (10)=dl (2,ic_c(10),n10);m  (10)=dm (2,ic_c(10),n10)
 s  (10)=ds (2,ic_c(10),n10);ms (10)=dms(2,ic_c(10),n10)
 if((l_c(10).eq.dll(ic_c(10))).and.(s_c (10).eq.dss(ic_c(10))) ) then
 if((m_c(10).eq.dmm(ic_c(10),n10)).and.(ms_c(10).eq.dmsms(ic_c(10),n10)))then

 if (ne.eq.9) then
   ic_c(9)=icc;l_c(9)=l_tot;s_c(9)=s_tot;m_c(9)=ml_tot;ms_c(9)=ms_tot
 endif
 do n9=1,count(ic_c(9))
 ic_c(8)=dic(  ic_c(9),n9)
 l_c (8)=dl (1,ic_c(9),n9);m_c (8)=dm (1,ic_c(9),n9)
 s_c (8)=ds (1,ic_c(9),n9);ms_c(8)=dms(1,ic_c(9),n9)
 l   (9)=dl (2,ic_c(9),n9);m   (9)=dm (2,ic_c(9),n9)
 s   (9)=ds (2,ic_c(9),n9);ms  (9)=dms(2,ic_c(9),n9)
 if ((l_c(9).eq.dll(ic_c(9)))   .and.(s_c (9).eq.dss  (ic_c(9))) )  then
 if ((m_c(9).eq.dmm(ic_c(9),n9)).and.(ms_c(9).eq.dmsms(ic_c(9),n9)))then

 if (ne.eq.8) then
   ic_c(8)=icc;l_c(8)=l_tot;s_c(8)=s_tot;m_c(8)=ml_tot;ms_c(8)=ms_tot
 endif
 do n8=1,count(ic_c(8))
 ic_c(7)=dic(  ic_c(8),n8)
 l_c (7)=dl (1,ic_c(8),n8);m_c (7)=dm (1,ic_c(8),n8)
 s_c (7)=ds (1,ic_c(8),n8);ms_c(7)=dms(1,ic_c(8),n8)
 l   (8)=dl (2,ic_c(8),n8);m   (8)=dm (2,ic_c(8),n8)
 s   (8)=ds (2,ic_c(8),n8);ms  (8)=dms(2,ic_c(8),n8)
 if ((l_c(8).eq.dll(ic_c(8)))   .and.(s_c (8).eq.dss  (ic_c(8))) )  then
 if ((m_c(8).eq.dmm(ic_c(8),n8)).and.(ms_c(8).eq.dmsms(ic_c(8),n8)))then

 if (ne.eq.7) then
   ic_c(7)=icc;l_c(7)=l_tot;s_c(7)=s_tot;m_c(7)=ml_tot;ms_c(7)=ms_tot
 endif
 do n7=1,count(ic_c(7))
 ic_c(6)=dic(  ic_c(7),n7)
 l_c (6)=dl (1,ic_c(7),n7);m_c (6)=dm (1,ic_c(7),n7)
 s_c (6)=ds (1,ic_c(7),n7);ms_c(6)=dms(1,ic_c(7),n7)
 l   (7)=dl (2,ic_c(7),n7);m   (7)=dm (2,ic_c(7),n7)
 s   (7)=ds (2,ic_c(7),n7);ms  (7)=dms(2,ic_c(7),n7)
 if ((l_c(7).eq.dll(ic_c(7)))   .and.(s_c (7).eq.dss  (ic_c(7))) )  then
 if ((m_c(7).eq.dmm(ic_c(7),n7)).and.(ms_c(7).eq.dmsms(ic_c(7),n7)))then

 if (ne.eq.6) then
   ic_c(6)=icc;l_c(6)=l_tot;s_c(6)=s_tot;m_c(6)=ml_tot;ms_c(6)=ms_tot
 endif
 do n6=1,count(ic_c(6))
 ic_c(5)=dic(  ic_c(6),n6)
 l_c (5)=dl (1,ic_c(6),n6);m_c (5)=dm (1,ic_c(6),n6)
 s_c (5)=ds (1,ic_c(6),n6);ms_c(5)=dms(1,ic_c(6),n6)
 l   (6)=dl (2,ic_c(6),n6);m   (6)=dm (2,ic_c(6),n6)
 s   (6)=ds (2,ic_c(6),n6);ms  (6)=dms(2,ic_c(6),n6)
 if ((l_c(6).eq.dll(ic_c(6)))   .and.(s_c (6).eq.dss  (ic_c(6))) )  then
 if ((m_c(6).eq.dmm(ic_c(6),n6)).and.(ms_c(6).eq.dmsms(ic_c(6),n6)))then

 if (ne.eq.5) then
   ic_c(5)=icc;l_c(5)=l_tot;s_c(5)=s_tot;m_c(5)=ml_tot;ms_c(5)=ms_tot
 endif
 do n5=1,count(ic_c(5))
 ic_c(4)=dic(  ic_c(5),n5)
 l_c (4)=dl (1,ic_c(5),n5);m_c (4)=dm (1,ic_c(5),n5)
 s_c (4)=ds (1,ic_c(5),n5);ms_c(4)=dms(1,ic_c(5),n5)
 l   (5)=dl (2,ic_c(5),n5);m   (5)=dm (2,ic_c(5),n5)
 s   (5)=ds (2,ic_c(5),n5);ms  (5)=dms(2,ic_c(5),n5)
 if ((l_c(5).eq.dll(ic_c(5)))   .and.(s_c (5).eq.dss  (ic_c(5))) )  then
 if ((m_c(5).eq.dmm(ic_c(5),n5)).and.(ms_c(5).eq.dmsms(ic_c(5),n5)))then

 if (ne.eq.4) then
   ic_c(4)=icc;l_c(4)=l_tot;s_c(4)=s_tot;m_c(4)=ml_tot;ms_c(4)=ms_tot
 endif
 do n4=1,count(ic_c(4))
 ic_c(3)=dic(  ic_c(4),n4)
 l_c (3)=dl (1,ic_c(4),n4);m_c (3)=dm (1,ic_c(4),n4)
 s_c (3)=ds (1,ic_c(4),n4);ms_c(3)=dms(1,ic_c(4),n4)
 l   (4)=dl (2,ic_c(4),n4);m   (4)=dm (2,ic_c(4),n4)
 s   (4)=ds (2,ic_c(4),n4);ms  (4)=dms(2,ic_c(4),n4)
 if ((l_c(4).eq.dll(ic_c(4)))   .and.(s_c (4).eq.dss  (ic_c(4))) )  then
 if ((m_c(4).eq.dmm(ic_c(4),n4)).and.(ms_c(4).eq.dmsms(ic_c(4),n4)))then

 if (ne.eq.3) then
   ic_c(3)=icc;l_c(3)=l_tot;s_c(3)=s_tot;m_c(3)=ml_tot;ms_c(3)=ms_tot
 endif
 do n3=1,count(ic_c(3))
 ic_c(2)=dic(  ic_c(3),n3)
 l_c (2)=dl (1,ic_c(3),n3);m_c (2)=dm (1,ic_c(3),n3)
 s_c (2)=ds (1,ic_c(3),n3);ms_c(2)=dms(1,ic_c(3),n3)
 l   (3)=dl (2,ic_c(3),n3);m   (3)=dm (2,ic_c(3),n3)
 s   (3)=ds (2,ic_c(3),n3);ms  (3)=dms(2,ic_c(3),n3)
 if ((l_c(3).eq.dll(ic_c(3)))   .and.(s_c (3).eq.dss  (ic_c(3))) )  then
 if ((m_c(3).eq.dmm(ic_c(3),n3)).and.(ms_c(3).eq.dmsms(ic_c(3),n3)))then

 if (ne.eq.2) then
   ic_c(2)=icc;l_c(2)=l_tot;s_c(2)=s_tot;m_c(2)=ml_tot;ms_c(2)=ms_tot
 endif
 do n2=1,count(ic_c(2))
 l   (1)=dl (1,ic_c(2),n2);m   (1)=dm (1,ic_c(2),n2)! 1st electron data 
 s   (1)=ds (1,ic_c(2),n2);ms  (1)=dms(1,ic_c(2),n2)! from 2 elec.coup.
 l   (2)=dl (2,ic_c(2),n2);m   (2)=dm (2,ic_c(2),n2)
 s   (2)=ds (2,ic_c(2),n2);ms  (2)=dms(2,ic_c(2),n2)
 if ((l_c(2).eq.dll(ic_c(2)))   .and.(s_c (2).eq.  dss(ic_c(2))) )  then
 if ((m_c(2).eq.dmm(ic_c(2),n2)).and.(ms_c(2).eq.dmsms(ic_c(2),n2)))then

!--- check that we are creating the correct state. If this breaks big prob.
  if (sum(m(1:ne) )-ml_tot.ne.0) write(*,*) 'Problem !!!!!!!!!'
  if (sum(ms(1:ne))-ms_tot.ne.0) write(*,*) 'Problem !!!!!!!!!'

!--- check for violation of exclusion - these terms should sum to zero
  kkk=201
  do i=1,ne
  do j=i+1,ne
   if((m(i).eq.m(j)).and.(ms(i).eq.ms(j))) kkk=202
  enddo
  enddo
        
  ms_=0.5d0*dble(ms)

  x=0.d0
  if (ne.eq.10) then
   x=coef(ic_c(2),n2)*coef(ic_c(3),n3)*coef(ic_c(4),n4)* &
 &   coef(ic_c(5),n5)*coef(ic_c(6),n6)*coef(ic_c(7),n7)* &
 &   coef(ic_c(8),n8)*coef(ic_c(9),n9)*coef(ic_c(10),n10)
  elseif (ne.eq.9) then
   x=coef(ic_c(2),n2)*coef(ic_c(3),n3)*coef(ic_c(4),n4)* &
 &   coef(ic_c(5),n5)*coef(ic_c(6),n6)*coef(ic_c(7),n7)* &
 &   coef(ic_c(8),n8)*coef(ic_c(9),n9)
  elseif (ne.eq.8) then
   x=coef(ic_c(2),n2)*coef(ic_c(3),n3)*coef(ic_c(4),n4)* &
 &   coef(ic_c(5),n5)*coef(ic_c(6),n6)*coef(ic_c(7),n7)* &
 &   coef(ic_c(8),n8)
  elseif (ne.eq.7) then
   x=coef(ic_c(2),n2)*coef(ic_c(3),n3)*coef(ic_c(4),n4)* &
 &   coef(ic_c(5),n5)*coef(ic_c(6),n6)*coef(ic_c(7),n7)
  elseif (ne.eq.6) then
   x=coef(ic_c(2),n2)*coef(ic_c(3),n3)*coef(ic_c(4),n4)* &
 &   coef(ic_c(5),n5)*coef(ic_c(6),n6)
  elseif (ne.eq.5) then
   x=coef(ic_c(2),n2)*coef(ic_c(3),n3)*coef(ic_c(4),n4)* &
 &   coef(ic_c(5),n5)
  elseif (ne.eq.4) then
   x=coef(ic_c(2),n2)*coef(ic_c(3),n3)*coef(ic_c(4),n4)
  elseif (ne.eq.3) then
   x=coef(ic_c(2),n2)*coef(ic_c(3),n3)
  elseif (ne.eq.2) then
   x=coef(ic_c(2),n2)
  endif


  if(kkk.eq.201)coef_targ(itarg)=coef_targ(itarg)+x



 endif; endif; enddo

 endif; endif; enddo

 endif; endif; enddo

 endif; endif; enddo

 endif; endif; enddo

 endif; endif; enddo

 endif; endif; enddo

 endif; endif; enddo

 endif; endif; enddo






!-----------------------------------------------------------------------------
!--- This section writes out the coupling lists and coefficients
!--- used to create the determinant. Probably only useful for finding
!--- any mistakes in the code.

!--- now have coeffs. for 1-2,12-3,123-4 etc. coupling ---
!--- use these to expand out in terms of Y_lm.X(s,ms) products
!do ic=start(ne-1),stop(ne-1)
!do ic=1,stop(ne-1)
! write(*,*) 'ic=',ic
! write(*,'(i2,a1,i1)') dss(ic)+1,spec(dll(ic)),dvv(ic)
! write(*,*) 'No. products in sum',count(ic)
! write(*,*) 'L   S   V  M  Ms  l1 l2 v1 v2  m1 m2  ms1  ms2   s1   s2  '
! do n=1,count(ic)
!   l1 =dl(1,ic,n); l2 =dl(2,ic,n)
!   s1 =ds(1,ic,n); s2 =ds(2,ic,n)
!   v1 =dv(1,ic,n); v2 =dv(2,ic,n)
!   m1 =dm(1,ic,n); m2 =dm(2,ic,n)
!   ms1=dms(1,ic,n);ms2=dms(2,ic,n)
!   l12 =dll  (ic)
!   s12 =dss  (ic)
!   v12 =dvv  (ic)
!   m12 =dmm  (ic,n)
!   ms12=dmsms(ic,n)
!   s1_  =0.5d0*dble(s1)
!   s2_  =0.5d0*dble(s2)
!   ms1_ =0.5d0*dble(ms1)
!   ms2_ =0.5d0*dble(ms2)
!   s12_ =0.5d0*dble(s12)
!   ms12_=0.5d0*dble(ms12)
!   z_ls=coef(ic,n)
!   write(*,102)  l12,s12_,v12,m12,ms12_,l1,l2,v1,v2,        &
!&        m1,m2,ms1_,ms2_,s1_,s2_,z_ls,dic(ic,n)
! enddo
! write(*,*) 
!enddo






!------ END BIG LOOP OVER TARGETS
 enddo








 fact=1.d0
 do j=1,ne
  fact=fact*dble(j)
 enddo

!--- provide only those determinants in sum which have non-zero coefficients
 i=0
 do itarg=1,n_targ
  if (abs(coef_targ(itarg)).gt.1d-6) then
   i=i+1
   coef_subshell(i)=coef_targ(itarg)
   do j=1,ne
    ml_subshell(j,i)=ml_targ(j,itarg)
    ms_subshell(j,i)=ms_targ(j,itarg)
   enddo
  endif
 enddo
 num=i







101   format(2(i3,1x,f4.1),1x,4(i2,1x),4(f4.1,1x),1x,f8.5)
102   format(i2,1x,f4.1,1x,i2,1x,i2,1x,f4.1,1x,6(i2,1x),4(f4.1,1x), &
    &        1x,f8.5,1x,i3)

 end subroutine subshell









 subroutine signature2det(sig,m,ms,ne,n_up,n_down,l_shell,ne2)  
 integer, intent(in)  :: sig,ne,n_up,n_down,l_shell
 integer, intent(out) :: ne2
 integer              :: i
 integer              :: m(:),ms(:)
 integer              :: iflag

!-- call this routine with signature of determinant
!-- ( range 0 - 4*l_shell+1, see det2signature function)
!--  and returns the determinant as the product that is ordered
!-- in decreasing m_l and m_s.

 nn=4*l_shell+2

 i=sig
 iii=0
 ne2=0
 m(:)=0;ms(:)=0
 do j=nn-1,0,-1
   if ( (i/2**j).eq.1) then
    ne2=ne2+1
    iii=iii+1
!-- Avoid going outside array bounds for l=>3.
    if (ne2.le.ne) then
     m (iii)=j/2-l_shell
     ms(iii)=2*mod(j,2)-1
    endif
   endif
!------ write(*,*) j,i/2**j, j/2-l_shell, 2*mod(j,2)-1
   i=mod(i,2**j)
 enddo

!-------- Note sorting only works correctly if ne=ne2
 do
  iflag=0
  do i=1,ne-1
   if (ms(i).lt.ms(i+1)) then
     i1=ms(i);i2=ms(i+1)
     ms(i)=i2;ms(i+1)=i1
     i1=m (i);i2=m (i+1)
     m (i)=i2;m (i+1)=i1
     iflag=1
   endif
  enddo
  if (iflag.eq.0) exit
 enddo

 do
  iflag=0
  do i=1,n_up-1
   if (m (i).lt.m (i+1)) then
     i1=ms(i);i2=ms(i+1)
     ms(i)=i2;ms(i+1)=i1
     i1=m (i);i2=m (i+1)
     m (i)=i2;m (i+1)=i1
     iflag=1
   endif
  enddo
  if (iflag.eq.0) exit
 enddo

 do
  iflag=0
  do i=n_up+1,ne-1
   if (m (i).lt.m (i+1)) then
     i1=ms(i);i2=ms(i+1)
     ms(i)=i2;ms(i+1)=i1
     i1=m (i);i2=m (i+1)
     m (i)=i2;m (i+1)=i1
     iflag=1
   endif
  enddo
  if (iflag.eq.0) exit
 enddo

     
 end subroutine signature2det









 integer function det2signature(m,ms,ne,l_shell)
 integer, intent(in)  :: ne,l_shell
 integer, intent(in)  :: m(:),ms(:)
 integer              :: j,jj

!--- calc. number that is the same for all permutations 
!---  that characterises the determinant

 det2signature=0

!---  In what follows, jj has range 0..4*l_shell+1, and describes 
!---   one electron (ml,ms) pair.
!---  FOR A DETERMINANT this occurs once or not at all,
!---  (0 or 1)  so if we make a 4*l_shell+2 digit binary number 
!---  from this, it is unique to the determinat

 do j=1,ne
  jj=2*(m(j)+l_shell) + (ms(j)+1)/2
  det2signature=det2signature+2**jj
 enddo

 end function det2signature







 Subroutine get_index(index,n,num,li,si,l_shell)
 implicit doubleprecision (a-h,o-z)
 external glcons,termls,trmf
 integer    :: n,num,li,si
 integer    :: qq,ss,ll,l_shell

!-- returns index of nth in the list that is (2S+1)L state

 select case (l_shell)
 case (0)
  i1= 1;i2= 2
 case (1)
  i1= 3;i2= 8
 case (2)
  i1= 9;i2=40
 case (3)
  i1= 1;i2=238
 case (4)
  i1= 1;i2=10
 case (5)
  i1=11;i2=22
 case (6)
  i1=23;i2=36
 case (7)
  i1=37;i2=52
 case (8)
  i1=53;i2=70
 case (9)
  i1=71;i2=90
 end select

 num=0;index=0
 do i=i1,i2
   if(l_shell.lt.3)then
     kt=mt(i)
   elseif(l_shell.eq.3)then
     kt=m76(i)
   else
     kt=mt3(i)
   endif
   qq=jthn(kt,3,100)
   ss=jthn(kt,2,100)
   ll=jthn(kt,1,100)
   if ((ss.eq.si).and.(ll.eq.2*li)) then
     num=num+1
     if (num.eq.n) index=i
   endif
 enddo

 end subroutine get_index






 Subroutine get_index_lsv(index,li,si,v,l_shell)
 implicit doubleprecision (a-h,o-z)
 external glcons,termls,trmf
 integer    :: li,si,v
 integer    :: qq,ll,ss,vv,l_shell

!-- returns index of TERM in list that is (2S+1)Lv states
!-- if there are several (eg f) it returns the last one found.
!-- for ne <=2 (probably <=4 and >=10) this is not a problem as
!-- TERMS are unique - only one (2S+1)Lv exists.
 select case (l_shell)
 case (0)
  i1= 1;i2= 2
 case (1)
  i1= 3;i2= 8
 case (2)
  i1= 9;i2=40
 case (3)
  i1= 1;i2=238
 case (4)
  i1= 1;i2=10
 case (5)
  i1=11;i2=22
 case (6)
  i1=23;i2=36
 case (7)
  i1=37;i2=52
 case (8)
  i1=53;i2=70
 case (9)
  i1=71;i2=90
 end select
 
 num=0;index=0
 do i=i1,i2
   if(l_shell.lt.3)then
     kt=mt(i)
   elseif(l_shell.eq.3)then
     kt=m76(i)
   else
     kt=mt3(i)
   endif
   qq=jthn(kt,3,100)
   ss=jthn(kt,2,100)
   ll=jthn(kt,1,100)
   vv=2*l_shell+2-(qq+1)
   if ((ss.eq.si).and.(ll.eq.2*li).and.(vv.eq.v)) then
     index=i
   endif
 enddo

 end subroutine get_index_lsv








 Subroutine mk_cfp(cfp,li,n,jb,jk)
 implicit doubleprecision (a-h,o-z)
 integer li
 external glcons,termls,trmf

 if(li.ne.3) then
   call rumt(jb,li,iqb,isb,ilb)
   call rumt(jk,li,iqk,isk,ilk)
 else
   call rumt67(jb,nrb,iqb,isb,ilb)
   call rumt67(jk,nrk,iqk,isk,ilk)
 endif
 call sls(li,jb,iqb,ilb,isb,jk,iqk,ilk,isk,sss)
 qk=dble(iqk)*0.5
 qb=dble(iqb)*0.5
 qmk=dble(n-1-2*li-1)*0.5
 qmb=dble(n-2*li-1)*0.5
 qm=0.5
 call c0t5s(qk,qmk,qm,qb,qmb,a)
 cfp=sss*a/dsqrt(dble(n*(iqb+1)*(ilb+1)*(isb+1)))
 if(((n-1)/2)*2.ne.n-1) cfp=-cfp

 end subroutine mk_cfp









 Subroutine cppl(fact,ic,n,l1,l2,l12,s1,s2,s12,                  &
 &               count,coef,dl,ds,dm,                            &
 &               dms,dmm,dmsms,v1,v2,dv,dic,icc,                 &
 &               lll,m,ms,itarg)
 implicit doubleprecision (a-h,o-z)
 integer          :: ic,n,v1,v2
 integer,intent(in):: lll
 integer          :: l1,l2,l12,m1 ,m2 ,m12
 integer          :: s1,s2,s12,ms1,ms2,ms12
 integer          :: i,itarg
 integer          :: count(:)
 real(kind(1.d0)) :: coef(:,:)
 integer          :: dl (:,:,:)
 integer          :: ds (:,:,:)
 integer          :: dm (:,:,:)
 integer          :: dv (:,:,:)
 integer          :: dms(:,:,:)
 integer          :: dmm  (:,:)
 integer          :: dmsms(:,:)
 integer          :: dic  (:,:)
 integer,intent(in):: m(:,:),ms(:,:)

!-- <l1m1l2m2|LMmlms> stored in coef, indexed by d..() arrays
!-- coef is coef of nth term in expansion

 do m12=-l12,l12
 do ms12=-s12,s12,2

 if (lll.ne.2) then

  do ms1=-s1,s1,2
  do m1=-l1,l1
    ms2=ms(lll,itarg)
    m2 =m (lll,itarg)

    x_l=cleb(2*l1,2*m1,2*l2,2*m2,2*l12,2*m12)
    y_s=cleb(s1,ms1,s2,ms2,s12,ms12)
    z_ls=x_l*y_s

    if (abs(z_ls).gt.1d-6) then
      n=n+1
      coef(ic,n)=fact*z_ls
      dl (1,ic,n)=l1 ;dl (2,ic,n)=l2
      ds (1,ic,n)=s1 ;ds (2,ic,n)=s2
      dv (1,ic,n)=v1 ;dv (2,ic,n)=v2
      dm (1,ic,n)=m1 ;dm (2,ic,n)=m2
      dms(1,ic,n)=ms1;dms(2,ic,n)=ms2
      dmm  (ic,n)=m12
      dmsms(ic,n)=ms12
      dic  (ic,n)=icc

    endif


  enddo
  enddo

 else

    ms1=ms(lll-1,itarg)
    m1 =m (lll-1,itarg)
    ms2=ms(lll,itarg)
    m2 =m (lll,itarg)

    x_l=cleb(2*l1,2*m1,2*l2,2*m2,2*l12,2*m12)
    y_s=cleb(s1,ms1,s2,ms2,s12,ms12)
    z_ls=x_l*y_s

    if (abs(z_ls).gt.1d-6) then
      n=n+1
      coef(ic,n)=fact*z_ls
      dl (1,ic,n)=l1 ;dl (2,ic,n)=l2
      ds (1,ic,n)=s1 ;ds (2,ic,n)=s2
      dv (1,ic,n)=v1 ;dv (2,ic,n)=v2
      dm (1,ic,n)=m1 ;dm (2,ic,n)=m2
      dms(1,ic,n)=ms1;dms(2,ic,n)=ms2
      dmm  (ic,n)=m12
      dmsms(ic,n)=ms12
      dic  (ic,n)=icc
    endif

 endif

 enddo
 enddo


 count(ic)=n
!write(100,*) ic,n

 end subroutine cppl






end module subshell_mod

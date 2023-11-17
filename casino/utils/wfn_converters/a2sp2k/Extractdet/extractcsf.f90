module extractcsf_mod

 use subshell_mod
 use error_mod

contains

 subroutine extractcsf(no_csf,iuf,ouf1,ouf2,ouf3,l_tot,s_tot,ml_tot,   &
      &                ms_tot,no_up,no_down,no_det)
 implicit none
 integer,intent(in):: no_csf,iuf,ouf1,ouf2,ouf3
 integer          :: max_targ,max_subshell,max_csf,max_wk
 parameter (max_targ=256,max_subshell=50,max_csf=2048,max_wk=65536)
 integer          :: iii,i,j,jj,k,s,l,v,n
 integer          :: ml_targ(16,max_targ),ms_targ(16,max_targ)
 real(kind(1.d0)) :: coef(max_targ)
 integer          :: i1,i2,ii1,ii2,no_det,no_electron,no_up,no_down
 integer          :: no_wk
 integer          :: n_shell,l_shell,ne,ic,icc,index,iflag
 integer          :: l_tot,s_tot,ml_tot,ms_tot,v_tot
 integer          :: l1,l2,m1,m2,s1,s2,ms1,ms2,v1,v2
 integer          :: l12,m12,s12,ms12,v12,p_num
 integer          :: itarg,itarg_1,itarg_2
 real(kind(1.d0)) :: t1,x,fact,x_l,y_s,z_ls
 integer          :: num(16)
 integer          :: p_max(16),pp(0:16),ppp
 integer          :: n_sub(16)
 integer          :: l_sub(16)
 integer          :: s_sub_tot(16)
 integer          :: l_sub_tot(16)
 integer          :: v_sub_tot(16)
 integer          :: dss    (0:16)
 integer          :: dll    (0:16)
 integer          :: ml_subshell(20,max_subshell)
 integer          :: ms_subshell(20,max_subshell)
 real(kind(1.d0)) :: coef_subshell (max_subshell)
 integer ::    n_csf(max_subshell,max_csf),  n_wk (max_subshell,max_wk)
 integer ::    l_csf(max_subshell,max_csf),  l_wk (max_subshell,max_wk)
 integer ::   ml_csf(max_subshell,max_csf), ml_wk (max_subshell,max_wk)
 integer ::    s_csf(max_subshell,max_csf),  s_wk (max_subshell,max_wk)
 integer ::   ms_csf(max_subshell,max_csf), ms_wk (max_subshell,max_wk)
 integer          :: mask_clsd(max_subshell)
 real(kind(1.d0)) :: coef_csf(max_csf)
 real(kind(1.d0)) :: coef_expansion
 complex(kind(1.d0)) :: coef_wk(max_wk)
 complex(kind(1.d0)) :: rm1,c(2)
 real(kind(1.d0)) :: root2
 integer          :: occ(16),nocc,nocc_clsd,nocc_opn,n_targ
 character(3)     :: e_clsd(8),elc(8),couple(16)
 character        :: spec(0:9),spec_lower(0:9)
 character(2)     :: dig(-9:9)
 data spec/'S','P','D','F','G','H','I','K','L','M'/
 data spec_lower/'s','p','d','f','g','h','i','k','l','m'/
 data dig/'-9','-8','-7','-6','-5','-4','-3','-2','-1',         &
     &         ' 0',' 1',' 2',' 3',' 4',' 5',' 6',' 7',' 8',' 9'/

 rm1=(0.d0,1.d0)
 root2=dsqrt(2.d0)


!--- read in closed shell data part of CSF from cfg.inp file.
 open(iuf,file ='cfg.inp')
 read(iuf,*)
 read(iuf,'(16(1x,a3))') (e_clsd(j),j=1,8)
 close(iuf)

 nocc_clsd = 0
 do while(e_clsd(nocc_clsd+1).ne. '   ')
  nocc_clsd = nocc_clsd + 1
 enddo

 do j=1,nocc_clsd
  call spec2nl(e_clsd(j),n,l)
  n_sub(j)=n; l_sub(j)=l
  occ(j)=4*l+2                                            ! closed shell
  s_sub_tot(j)=0;l_sub_tot(j)=0;v_sub_tot(j)=0            ! closed shell
 enddo
! coupling closed-closed subshells
 do j=1,nocc_clsd-1
  dss(j)=0
  dll(j)=0
 enddo

!-- this writes out closed shells only.
!do j=1,nocc_clsd
!  write(*,'(i2,a1,a1,i2,a1,2x,i2,a1,i1)')                              &
!& n_sub(j),spec_lower(l_sub(j)),'(',occ(j),')  ',                      &
!& s_sub_tot(j)+1,spec(l_sub_tot(j)),v_sub_tot(j)
!enddo
!write(*,'(15(1x,a3))') (dig(dss(j)+1)//spec(dll(j)),j=1,nocc_clsd)
!write(*,*)

!--------------nocc_clsd=0

!--- read in open shell data part of CSF from cfg.inp file.
 open(iuf,file ='cfg.inp')
 read(iuf,*);read(iuf,*)
! step forward to CSF no. no_csf in list
 do i=1,no_csf-1
  read(iuf,*);read(iuf,*)
 enddo
! CSF open shells
 read(iuf,'(16(1x,a3,1x,i2,1x))') (elc(j),occ(nocc_clsd+j),j=1,8)
! read subshell TERMS,then TERMS resulting from coupling
! 1st (1:nocc_opn) are subshell,  (nocc_opn+1:2*nocc_opn-1) are coupling
 read(iuf,'(16(1x,a3))') (couple(j),j=1,16)
 close(iuf)

!--- read the expansion coefficient for this CSF from fort.70
!open(iuf,file ='fort.70')
! read(iuf,*);read(iuf,*)
! do i=1,no_csf-1
!   read(iuf,*)
! enddo
! read(iuf,*) coef_expansion
!close(iuf)
!--- Use standard fort.60 file instead.
 open(iuf,file ='fort.60')
  read(iuf,*);read(iuf,*)
  read(iuf,'(7F11.7)') (t1, i=1,no_csf)
 close(iuf)
 coef_expansion=t1

 nocc_opn = 0
 do while(elc(nocc_opn+1).ne. '   ')
  nocc_opn = nocc_opn + 1
 enddo


!--- elc - spec notation of shell nl,occ - is occupancy 
!--- couple is term of coupling 
!--- 1st TERM of 1st subshell, then 2nd, 3rd   etc...     then
!--- TERM of 1-2 coupling, then 12-3 coupling  etc...
!--- left-right coupling

 do j=1,nocc_opn
  call spec2nl(elc(j),n,l)
  n_sub(nocc_clsd+j)=n; l_sub(nocc_clsd+j)=l
  call spec2slv(couple(j),s,l,v)
  s_sub_tot(nocc_clsd+j)=s
  l_sub_tot(nocc_clsd+j)=l
  v_sub_tot(nocc_clsd+j)=v
 enddo
! result of coupling open subshells is returned by spec2slv
 do j=1,nocc_opn-1
  call spec2slv(couple(j+nocc_opn),s,l,v)
  dss(nocc_clsd+j)=s
  dll(nocc_clsd+j)=l
 enddo

!--- coupling supplied in input does not include coupling of 1st 
!--- subshell to nothing, or the coupling of the 1st open shell to 
!--- the closed shells. This is put in in the next two lines

 dss(0)=s_sub_tot(1)                     ! coupling of 1st subshell to nothing
 dll(0)=l_sub_tot(1)                     ! just TERM of 1st subshell

 dss(nocc_clsd)=s_sub_tot(nocc_clsd+1)   ! coupling of 1st open shell to closed
                                         !   shells
 dll(nocc_clsd)=l_sub_tot(nocc_clsd+1)   ! just TERM of 1st open subshell
                                         !  (closed shell are all 1S)


 nocc=nocc_clsd+nocc_opn


!--- CSF information has been put into:
!--- n_sub(j) - n no. of subshell j
!--- l_sub(j) - l no. of subshell j
!--- dll(j)   - l no. of subshell j-1 coupled to subshell j
!--- dss(j)   - s no. of subshell j-1 coupled to subshell j
!------ write out configuration and couplings to standard out
 do j=1,nocc
   write(*,'(i2,a1,a1,i2,a1,2x,i2,a1,i1)')                              &
 & n_sub(j),spec_lower(l_sub(j)),'(',occ(j),')  ',                      &
 & s_sub_tot(j)+1,spec(l_sub_tot(j)),v_sub_tot(j)
 enddo
 write(*,'(16(1x,a3))') (dig(dss(j)+1)//spec(dll(j)),j=0,nocc-1)

!----------------------------------------------------------------------
!---- SET Z ANGULAR MOMENTUM OF csf HERE, and number of electrons -----
 if (l_tot.ne.dll(nocc-1)) then
  write(*,*) 'mismatch of L';stop
 endif
 if (s_tot.ne.dss(nocc-1)) then
  write(*,*) 'mismatch of S';stop
 endif
 l_tot=dll(nocc-1)
 s_tot=dss(nocc-1)
 no_electron=sum(occ(1:nocc))
 no_up  =(no_electron+ms_tot)/2
 no_down=(no_electron-ms_tot)/2
 write(*,'(a,1x,a)')          'Doing TERM  =',dig(s_tot+1)//spec(l_tot)
 write(*,'(a,1x,i3,1x,f4.1)') '      Ml,Ms =',ml_tot,0.5d0*dble(ms_tot)
 write(*,'(a10,1x,i3)') 'Spin up  :',no_up
 write(*,'(a10,1x,i3)') 'Spin down:',no_down
 write(*,*)
!----------------------------------------------------------------------

!--- Set array mask_clsd. Set to one normally, but set to zero if 
!---  electron is a member of a closed shell.
!---  used later to make multiplication into chemical orbitals Xlm 
!---  more efficient.
 i=1
 do j=1,nocc
  if (occ(j).eq.4*l_sub(j)+2) then
   mask_clsd(i:i+occ(j))=0
  else
   mask_clsd(i:i+occ(j))=1
  endif
  i=i+occ(j)
 enddo
!write(*,'(50(1x,i2))') mask_clsd(1:no_electron)



!--- Find the ml/ms values for each subshell that add up to the total 
!---  values, ml_tot and ms_tot . These are the targets for coupling 
!---  subshells.
 write(*,*)
 pp(:)=0;p_max(:)=0;ppp=1
 do j=1,nocc
  l=l_sub_tot(j);s=s_sub_tot(j)
  p_max(j)=(s+1)*(2*l+1)
  ppp=ppp*p_max(j)
  pp(j) = ppp
  write(*,*) j,p_max(j),pp(j)
 enddo
 pp(0)=1

 n_targ=1
 do i=0,ppp-1
  j=i

  do k=nocc,1,-1
   l=l_sub_tot(k);s=s_sub_tot(k)
   iii=  j/pp(k-1)
   j=j-iii*pp(k-1)
   ml_targ(k,n_targ)=iii/(s+1)-l
   ms_targ(k,n_targ)=2*mod(iii,s+1)-s
  enddo

  iflag=1
  if (sum(ml_targ(1:nocc,n_targ)).ne.ml_tot) iflag=0
  if (sum(ms_targ(1:nocc,n_targ)).ne.ms_tot) iflag=0

  if (iflag.eq.1) then
   n_targ=n_targ+1
  endif

 enddo
 n_targ=n_targ-1
 if (n_targ.gt.max_targ) call error(1,n_targ,max_targ)

!--- Candidates are now stored in 
!---  ml_targ(1:nocc,itarg),ms_targ(1:nocc,itarg)
 




 coef(:)=1.d0

 do itarg=1,n_targ                        ! which target do we want to do?
  do k=1,nocc-1                           !j-loop (over couplings)

   l1 =dll(k-1)                           !parent 1
   s1 =dss(k-1)
   m1 =sum(ml_targ(1:k,itarg))
   ms1=sum(ms_targ(1:k,itarg))

   l2 =l_sub_tot(k+1)                     !parent 2
   s2 =s_sub_tot(k+1)
   m2 =ml_targ(k+1,itarg)
   ms2=ms_targ(k+1,itarg)

   l12 =dll(k)                            !child
   s12 =dss(k)
   m12 =m1+m2
   ms12=ms1+ms2

   x_l=cleb(2*l1,2*m1,2*l2,2*m2,2*l12,2*m12)
   y_s=cleb(s1,ms1,s2,ms2,s12,ms12)
   z_ls=x_l*y_s

!  write(*,*) spec(l1),m1,spec(l2),m2
!  write(*,*) dig(s12+1)//spec(l12),x_l,y_s
   if (abs(z_ls).gt.1d-6) then
    coef(itarg)=coef(itarg)*z_ls
   else
    coef(itarg)=0.d0
   endif
!  write(*,*) 'k=',k
!  write(*,*) s1+1, spec(l1), s2+1,spec(l2)
!  write(*,*) s12+1,spec(l12)
!  write(*,*) coef(itarg)

  enddo
 enddo

!--- remove zero coefficient ml_targ's from the list
 i=1
 do j=1,n_targ
  if(abs(coef(j)).gt.1d-6)then
    if(i.ne.j) then
     ml_targ(1:nocc,i)=ml_targ(1:nocc,j)
     ms_targ(1:nocc,i)=ms_targ(1:nocc,j)
     coef(i )=coef(j )
    endif
    i=i+1
  endif
 enddo
 n_targ=i-1    ! n_targ is now number in list
!--- products of CB coefficients now in coef(itarg), all non-zero




!--- Next expand out the `products' of sums
!--- by `product' i mean a representation of a determinant.
!--- This means that it is antisymetric by definition.

 n_csf(:,:)=0
 l_csf(:,:)=0;ml_csf(:,:)=0
 s_csf(:,:)=0;ms_csf(:,:)=0
 coef_csf(:)=0.d0
 ic=0


 do itarg_1=1,n_targ

  do j=nocc,1,-1
    n_shell=n_sub(j)              ! set n_shell of subshell j
    l_shell=l_sub(j)              ! set l_shell of subshell j
    l1 =l_sub_tot(j)              ! total L of subshell j
    s1 =s_sub_tot(j)              ! total S of subshell j
    m1 =ml_targ(j,itarg_1)        ! total Ml of subshell j
    ms1=ms_targ(j,itarg_1)        ! total Ms of subshell j
    v1 =v_sub_tot(j)              ! seniority of subshell j
    ne =occ(j)                    ! no electrons in subshell j

    call subshell(l_shell,ne,l1,m1,s1,ms1,v1,                    &
                ml_subshell,ms_subshell,coef_subshell,num(j))
    if (num(j).gt.max_subshell) call error(2,num(j),max_subshell)
    i1=sum(occ(1:j-1))+1;i2=sum(occ(1:j))

    if (j.eq.nocc) then
     do itarg_2=1,num(j)
       index=ic+itarg_2
       ml_csf(i1:i2,index)=ml_subshell(1:ne,itarg_2)
       ms_csf(i1:i2,index)=ms_subshell(1:ne,itarg_2)
       coef_csf(index)=coef_subshell(itarg_2)*coef(itarg_1)
     enddo
     p_num=num(j)
     ii1=i1
    else
     do itarg_2=2,num(j)
     do i=1,p_num
       index=ic+p_num*(itarg_2-1)+i
       ml_csf(ii1:no_electron,index)=ml_csf(ii1:no_electron,ic+i)
       ms_csf(ii1:no_electron,index)=ms_csf(ii1:no_electron,ic+i)
       coef_csf(index)=coef_csf(ic+i)
     enddo
     enddo
     ii1=i1

     do itarg_2=1,num(j)
     do i=1,p_num
       index=ic+p_num*(itarg_2-1)+i
       ml_csf(i1:i2,index)=ml_subshell(1:ne,itarg_2)
       ms_csf(i1:i2,index)=ms_subshell(1:ne,itarg_2)
       coef_csf(index)=coef_subshell(itarg_2)*coef_csf(index)
     enddo
     enddo
     p_num=p_num*num(j)
    endif
    
  enddo

!--- if target config itarg_1 has no entries on list for any j, then
!---  this does not include it on the list. Critical bit is the 
!---  condition on ic=ic+p_num part at the end.
  if (all( num(1:nocc).ne.0 )) then
   do i=1,p_num
   do j=nocc,1,-1
    i1=sum(occ(1:j-1))+1;i2=sum(occ(1:j))
    n_csf (i1:i2,ic+i)=n_sub(j)
    l_csf (i1:i2,ic+i)=l_sub(j)
    s_csf (i1:i2,ic+i)=1
   enddo
   enddo
   ic=ic+p_num
  endif

 enddo

!------ check that we have found something for this CSF
 no_det=ic
 if (no_det.eq.0) then
   write(*,*) 'Could not find this CSF !   Problem.'
   stop
 endif
 if (no_det.gt.max_csf) call error(3,no_det,max_csf)

 fact=1.d0
 do j=1,nocc
  fact=fact*dsqrt( factorial(occ(j)) )
 enddo
!---- next line is proper normalisation, so that the full many body WF
!----  is normalised.
!---- If it is commented out then each each permutation 
!----  is normalised seperately, so norm of full many body WF is N!.
!fact=fact/dsqrt( factorial(no_electron) )

 coef_csf(1:no_det)=coef_csf(1:no_det)*fact








!----------- apply a l+ ladder operator to sum of determinants ----
!write(*,*) 'Apply l+ ladder operator'
!icc=1
!do ic=1,no_det
! do j=1,no_electron
!  fact=dble(l_csf(j,ic)*(l_csf(j,ic)+1) - ml_csf(j,ic)*(ml_csf(j,ic)+1))
!  fact=sqrt(fact)
!  ml_wk(1:no_electron,icc)=ml_csf(1:no_electron,ic)
!  ml_wk(j,            icc)=ml_csf(j,            ic)+1
!  ms_wk(1:no_electron,icc)=ms_csf(1:no_electron,ic)
!   n_wk(1:no_electron,icc)= n_csf(1:no_electron,ic)
!   l_wk(1:no_electron,icc)= l_csf(1:no_electron,ic)
!  do jj=1,no_electron
!   if(jj.ne.j) then
!    if ( (ml_wk(j,icc).eq.ml_wk(jj,icc)).and.                        &
!         (ms_wk(j,icc).eq.ms_wk(jj,icc)).and.                        &
!         ( n_wk(j,icc).eq. n_wk(jj,icc)).and.                        &
!         ( l_wk(j,icc).eq. l_wk(jj,icc)) ) then
!       fact=0.d0
!    endif
!   endif
!  enddo
!  if (fact.ne.0.d0) then
!   coef_wk(icc)=coef_csf(ic)*fact
!   icc=icc+1
!  endif
! enddo
!enddo
!no_wk=icc-1
!do ic=1,no_wk
!do icc=ic+1,no_wk                       
! if (all( (ms_wk(1:no_electron,ic)-ms_wk(1:no_electron,icc)).eq.0 )) then
! if (all( (ml_wk(1:no_electron,ic)-ml_wk(1:no_electron,icc)).eq.0 )) then
!   coef_wk(ic)=coef_wk(ic)+coef_wk(icc)
!   coef_wk(icc)=0.d0
! endif
! endif
!enddo
!enddo
!icc=1
!do ic=1,no_wk
! if(abs(coef_wk(ic)).gt.1d-6)then
!   if(ic.ne.icc) then
!    n_wk (1:no_electron,icc)=n_wk (1:no_electron,ic)
!    l_wk (1:no_electron,icc)=l_wk (1:no_electron,ic)
!    s_wk (1:no_electron,icc)=s_wk (1:no_electron,ic)
!    ms_wk(1:no_electron,icc)=ms_wk(1:no_electron,ic)
!    ml_wk(1:no_electron,icc)=ml_wk(1:no_electron,ic)
!    coef_wk(icc )=coef_wk(ic)
!    coef_wk(ic  )=0.d0
!   endif
!   icc=icc+1
! endif
!enddo
!no_wk=icc-1
!write(*,*) 'No_det L+ phi:',no_wk
!do ic=1,no_wk
!   write(*,'(10(i3,1x))',advance='no') (ml_wk(jj,ic),jj=1,no_electron)
!   write(*,'(10(i3,1x))',advance='no') (ms_wk(jj,ic),jj=1,no_electron)
!   write(*,'(f10.7,1x)') coef_wk(ic)
!enddo
!write(*,*)
!----------- apply a s+ ladder operator to sum of determinants ----
!write(*,*) 'Apply s+ ladder operator'
!do ic=1,no_det
! do j=1,no_electron
!  fact=dble(s_csf(j,ic)*(s_csf(j,ic)+2) - ms_csf(j,ic)*(ms_csf(j,ic)+2))
!  fact=0.5d0*sqrt(fact)
!  ms_wk(1:no_electron,ic)=ms_csf(1:no_electron,ic)
!  ms_wk(j,ic)=ms_wk(j,ic)+2
!  do jj=1,no_electron
!   if(jj.ne.j) then
!    if ( (ml_csf(j,ic).eq.ml_csf(jj,ic)).and.                        &
!         (ms_wk (j,ic).eq.ms_csf(jj,ic)).and.                        &
!         ( n_csf(j,ic).eq. n_csf(jj,ic)).and.                        &
!         ( l_csf(j,ic).eq. l_csf(jj,ic)) ) then
!       fact=0.d0
!    endif
!   endif
!  enddo
!  if (fact.ne.0.d0) then
!   write(*,'(10(i3,1x))',advance='no') (ml_csf(jj,ic),jj=1,no_electron)
!   write(*,'(10(i3,1x))',advance='no') (ms_wk (jj,ic),jj=1,no_electron)
!   write(*,'(f10.7,1x)') coef_csf(ic)*fact
!  endif
! enddo
!enddo
!write(*,*)
!------ for cancelation to occur these must be re-ordered to same ml order
!------ does not do this at the moment, so cancelation is not obvious.


!------ sorting csf arrays.
!--- sort by spin quantum number 
!--- sort by nl quantum   number              (not done yet)
!--- convert to chemical orbitals for casino






!--- write out the CSF

 x=0.d0
 do ic=1,no_det
! write(*,'(i2,1x,a)') ic,'----------------------------------------------' 
! do j=1,no_electron
!  write(*,'(i3,1x)',advance='no') n_csf(j,ic)
! enddo;write(*,*)
! do j=1,no_electron
!  write(*,'(i3,1x)',advance='no') l_csf(j,ic)
! enddo;write(*,*)
! do j=1,no_electron
!  write(*,'(i3,1x)',advance='no') s_csf(j,ic)
! enddo;write(*,*)
  do j=1,no_electron
   write(*,'(i3,1x)',advance='no') ml_csf(j,ic)
  enddo
  do j=1,no_electron
   write(*,'(i2,1x)',advance='no') ms_csf(j,ic)
  enddo
  write(*,'(f10.7,1x)') coef_csf(ic)
  x=x+coef_csf(ic)**2
 enddo
 write(*,*) 'Sum ^2 :',x
 write(*,*) 'No_det Ylm:',no_det
 if (abs(1.d0-x).gt.1d-6) call error(21,0,0)





!--- Next we convert the sum of determinants specified 
!--- in the *_csf arrays into expansion in `chemical spherical harmonics' 
!--- instead of the normal spherical harmonics.
!--- In general result will be a sum of real determinants with 
!--- complex coefficients, instead of complex determinants with real
!--- coefficients.

 n_wk(:,:)=0;l_wk(:,:)=0;ml_wk(:,:)=0;s_wk(:,:)=0;ml_wk(:,:)=0
 coef_wk(:)=(0.d0,0.d0)

 ic=0

 do itarg_1=1,no_det

!--- 1.
!--- Multiply out the representation of the determinants
!---  as if they were just products.
!---  Integer array mask_clsd is 1 normally, but zero if electron is a member
!---  of a closed shell. For these electrons  a value of zero means
!---  we do not calculate terms in the expansion that we know will cancel.
!---  it does this by preventing the use of the c(2) parameter in the 
!---  expansion for these cases. (the coefficient of the flipped ml term)
  do j=no_electron,1,-1
    if(ml_csf(j,itarg_1).eq.0)then
      k=1;c(1)=1.d0;c(2)=0.d0
    else if(ml_csf(j,itarg_1).gt.0)then
      fact=dble( (-1)**ml_csf(j,itarg_1) )/root2
      if (mask_clsd(j).eq.1) then
       k=2;c(1)=fact;c(2)=rm1*fact
      else
       k=1;c(1)=root2*fact
      endif
    else if(ml_csf(j,itarg_1).lt.0)then
      fact=1.d0/root2
      if (mask_clsd(j).eq.1) then
       k=2;c(1)=-rm1*fact;c(2)=fact
      else
       k=1;c(1)=-root2*rm1*fact
      endif
    endif

    if (j.eq.no_electron) then
     do itarg_2=1,k
       index=ic+itarg_2
       ml_wk(j,index)=ml_csf(j,itarg_1)*(3-2*itarg_2) ! flip sign for itarg_2=2
       ms_wk(j,index)=ms_csf(j,itarg_1)
       coef_wk(index)=dcmplx(coef_csf(itarg_1))*c(itarg_2)
     enddo
     p_num=k
     ii1=j
    else
     do itarg_2=2,k
     do i=1,p_num
       index=ic+p_num*(itarg_2-1)+i
       ml_wk(ii1:no_electron,index)=ml_wk(ii1:no_electron,ic+i)
       ms_wk(ii1:no_electron,index)=ms_wk(ii1:no_electron,ic+i)
       coef_wk(index)=coef_wk(ic+i)
     enddo
     enddo
     ii1=i1

     do itarg_2=1,k
     do i=1,p_num
       index=ic+p_num*(itarg_2-1)+i
       ml_wk(j,index)=ml_csf(j,itarg_1)*(3-2*itarg_2) ! flip sign for itarg_2=2
       ms_wk(j,index)=ms_csf(j,itarg_1)
       coef_wk(index)=coef_wk(index)*c(itarg_2)
     enddo
     enddo
     p_num=p_num*k
    endif
    
  enddo

 if (ic+p_num.gt.max_wk) call error(4,ic+p_num,max_wk)

  do i=1,p_num
  do j=no_electron,1,-1
   n_wk (j,ic+i)=n_csf(j,itarg_1)
   l_wk (j,ic+i)=l_csf(j,itarg_1)
   s_wk (j,ic+i)=s_csf(j,itarg_1)
  enddo
  enddo

!--- 2.
!--- Remove terms that are symmetric on permutation - the ones that do not 
!---  contribute to determinants.
  do k=ic+1,ic+p_num                         ! set coefs of these to zero
   iflag=0
   do i=1,no_electron-1
   do j=i+1,no_electron
    if(n_wk (i,k).eq.n_wk (j,k))then
    if(l_wk (i,k).eq.l_wk (j,k))then
    if(ms_wk(i,k).eq.ms_wk(j,k))then
    if(ml_wk(i,k).eq.ml_wk(j,k))then
      iflag=1
    endif
    endif
    endif
    endif
   enddo
   enddo
   if (iflag.eq.1)coef_wk(k)=0.d0
  enddo

  icc=ic+1                                     ! removes zeroes from the list
  do k=ic+1,ic+p_num
   if(abs(coef_wk(k)).gt.1d-6)then
     if(k.ne.icc) then
      n_wk (1:no_electron,icc)=n_wk (1:no_electron,k)
      l_wk (1:no_electron,icc)=l_wk (1:no_electron,k)
      s_wk (1:no_electron,icc)=s_wk (1:no_electron,k)
      ms_wk(1:no_electron,icc)=ms_wk(1:no_electron,k)
      ml_wk(1:no_electron,icc)=ml_wk(1:no_electron,k)
      coef_wk(icc )=coef_wk(k)
      coef_wk(k   )=0.d0
     endif
     icc=icc+1
   endif
  enddo
  p_num=icc-(ic+1)
!-- p_num is now number of non-zero coeff that make up list.


  ic=ic+p_num
 enddo
 no_det=ic



!--- 3.
!--- Sort so each subshell is ordered in decreasing ml, spin up then spin down
 do ic=1,no_det
 fact=1.d0
 do 
  iflag=0
  do j=1,no_electron-1
   if(n_wk (j,ic).eq.n_wk (j+1,ic))then
   if(l_wk (j,ic).eq.l_wk (j+1,ic))then
   if(ms_wk(j,ic).eq.ms_wk(j+1,ic))then
   if(ml_wk(j,ic).lt.ml_wk(j+1,ic)) then
    i1= n_wk(j,ic);i2= n_wk(j+1,ic); n_wk(j,ic)=i2; n_wk(j+1,ic)=i1
    i1= l_wk(j,ic);i2= l_wk(j+1,ic); l_wk(j,ic)=i2; l_wk(j+1,ic)=i1
    i1= s_wk(j,ic);i2= s_wk(j+1,ic); s_wk(j,ic)=i2; s_wk(j+1,ic)=i1
    i1=ml_wk(j,ic);i2=ml_wk(j+1,ic);ml_wk(j,ic)=i2;ml_wk(j+1,ic)=i1
    i1=ms_wk(j,ic);i2=ms_wk(j+1,ic);ms_wk(j,ic)=i2;ms_wk(j+1,ic)=i1
    fact=-fact            !---- IS THIS CORRECT? - check...
    iflag=1
   endif
   endif
   endif
   endif
  enddo
 if(iflag.eq.0) exit
 enddo
 coef_wk(ic)=coef_wk(ic)*fact
 enddo

!--- 4.
!--- Add together coefficients of equal determinants.
!---  add coefs into 1st occurence, setting rest to zero.
 do ic=1,no_det
 do icc=ic+1,no_det                        
  if (all( (ms_wk(1:no_electron,ic)-ms_wk(1:no_electron,icc)).eq.0 )) then
  if (all( (ml_wk(1:no_electron,ic)-ml_wk(1:no_electron,icc)).eq.0 )) then
    coef_wk(ic)=coef_wk(ic)+coef_wk(icc)
    coef_wk(icc)=0.d0
  endif
  endif
 enddo
 enddo
!---  remove zero coefficient terms from list.
 icc=1
 do ic=1,no_det
  if(abs(coef_wk(ic)).gt.1d-6)then
    if(ic.ne.icc) then
     n_wk (1:no_electron,icc)=n_wk (1:no_electron,ic)
     l_wk (1:no_electron,icc)=l_wk (1:no_electron,ic)
     s_wk (1:no_electron,icc)=s_wk (1:no_electron,ic)
     ms_wk(1:no_electron,icc)=ms_wk(1:no_electron,ic)
     ml_wk(1:no_electron,icc)=ml_wk(1:no_electron,ic)
     coef_wk(icc )=coef_wk(ic)
     coef_wk(ic  )=0.d0
    endif
    icc=icc+1
  endif
 enddo
 no_det=icc-1


!---------- check up on lists?
!write(*,*) ml_tot
!x=0.d0
!do ic=1,no_det
!do j=1,no_electron
!  write(*,'(i3,1x)',advance='no') ml_wk(j,ic)
!enddo
!do j=1,no_electron
!  write(*,'(i2,1x)',advance='no') ms_wk(j,ic)
!enddo
!x=x+abs(coef_wk(ic)**2)
!write(*,'(2(f10.7,1x))') (coef_wk(ic))
!enddo
!write(*,*) 'Sum ^2 :',x
!write(*,*) 'No_det Xlm:',no_det


!--- 5.
!--- Take appropriate sum/differences of Ml, -Ml CSF's so we have 
!---  a real wavefunction that is consistent for all CSF in MD expansion.
!--- Put these into CSF arrays.
 if ( ml_tot.eq.0) then                             ! ml_tot=0
  if(no_det.gt.max_csf) call error(3,no_det,max_csf)
  icc=1
  do ic=1,no_det
   if ( abs(dble (coef_wk(ic))).gt.1d-6 ) then      ! Case of real coeff.
    n_csf (1:no_electron,icc)=n_wk (1:no_electron,ic)
    l_csf (1:no_electron,icc)=l_wk (1:no_electron,ic)
    s_csf (1:no_electron,icc)=s_wk (1:no_electron,ic)
    ms_csf(1:no_electron,icc)=ms_wk(1:no_electron,ic)
    ml_csf(1:no_electron,icc)=ml_wk(1:no_electron,ic)
    coef_csf(icc)=dble (coef_wk(ic))                ! since |+Ml> = |-Ml>
    icc=icc+1
   endif
   if ( abs(dimag(coef_wk(ic))).gt.1d-6 ) then      ! Case of imag coeff.
    n_csf (1:no_electron,icc)=n_wk (1:no_electron,ic)
    l_csf (1:no_electron,icc)=l_wk (1:no_electron,ic)
    s_csf (1:no_electron,icc)=s_wk (1:no_electron,ic)
    ms_csf(1:no_electron,icc)=ms_wk(1:no_electron,ic)
    ml_csf(1:no_electron,icc)=ml_wk(1:no_electron,ic)
    coef_csf(icc)=dimag(coef_wk(ic))                ! since |+Ml> = |-Ml>
    icc=icc+1
   endif
  enddo
  no_det=icc-1
 else if ( mod(ml_tot,2).eq.0) then                 ! even ml_tot - take real parts
  if(no_det/2.gt.max_csf) call error(3,no_det/2,max_csf)
  icc=1
  do ic=1,no_det
   if ( abs(dble (coef_wk(ic))).gt.1d-6 ) then      ! (|+Ml> + |-Ml>)/sqrt(+2)
!--if ( abs(dimag(coef_wk(ic))).gt.1d-6 ) then      ! (|+Ml> - |-Ml>)/sqrt(-2)
    n_csf (1:no_electron,icc)=n_wk (1:no_electron,ic)
    l_csf (1:no_electron,icc)=l_wk (1:no_electron,ic)
    s_csf (1:no_electron,icc)=s_wk (1:no_electron,ic)
    ms_csf(1:no_electron,icc)=ms_wk(1:no_electron,ic)
    ml_csf(1:no_electron,icc)=ml_wk(1:no_electron,ic)
    coef_csf(icc)=root2*dble (coef_wk(ic))          ! (|+Ml> + |-Ml>)/sqrt(+2)
!-- coef_csf(icc)=root2*dimag(coef_wk(ic))          ! (|+Ml> - |-Ml>)/sqrt(-2)
    icc=icc+1
   endif
  enddo
  no_det=icc-1
 else                                               ! odd ml_tot - take imag. parts
  if(no_det/2.gt.max_csf) call error(3,no_det/2,max_csf)
  icc=1
  do ic=1,no_det
   if ( abs(dimag(coef_wk(ic))).gt.1d-6 ) then       ! (|+Ml> + |-Ml>)/sqrt(-2)
!--if ( abs(dble (coef_wk(ic))).gt.1d-6 ) then       ! (|+Ml> - |-Ml>)/sqrt(+2)
    n_csf (1:no_electron,icc)=n_wk (1:no_electron,ic)
    l_csf (1:no_electron,icc)=l_wk (1:no_electron,ic)
    s_csf (1:no_electron,icc)=s_wk (1:no_electron,ic)
    ms_csf(1:no_electron,icc)=ms_wk(1:no_electron,ic)
    ml_csf(1:no_electron,icc)=ml_wk(1:no_electron,ic)
    coef_csf(icc)=root2*dimag(coef_wk(ic))           ! (|+Ml> + |-Ml>)/sqrt(-2)
!-- coef_csf(icc)=root2*dble (coef_wk(ic))           ! (|+Ml> - |-Ml>)/sqrt(+2)
    icc=icc+1
   endif
  enddo
  no_det=icc-1
 endif

!--- Reorder so it is all spin up, then spin down
 do ic=1,no_det
  fact=1.d0
  do 
   iflag=0
   do i=1,no_electron-1
    if(ms_csf(i,ic).lt.ms_csf(i+1,ic)) then
     i1= n_csf(i,ic);i2= n_csf(i+1,ic); n_csf(i,ic)=i2; n_csf(i+1,ic)=i1
     i1= l_csf(i,ic);i2= l_csf(i+1,ic); l_csf(i,ic)=i2; l_csf(i+1,ic)=i1
     i1= s_csf(i,ic);i2= s_csf(i+1,ic); s_csf(i,ic)=i2; s_csf(i+1,ic)=i1
     i1=ml_csf(i,ic);i2=ml_csf(i+1,ic);ml_csf(i,ic)=i2;ml_csf(i+1,ic)=i1
     i1=ms_csf(i,ic);i2=ms_csf(i+1,ic);ms_csf(i,ic)=i2;ms_csf(i+1,ic)=i1
     fact=-fact            !---- IS THIS CORRECT? - check...
     iflag=1
    endif
   enddo
  if(iflag.eq.0) exit
  enddo
  coef_csf(ic)=coef_csf(ic)*fact
 enddo




 write(*,*)
 x=0.d0
 do ic=1,no_det
!do j=1,no_electron
!  write(*,'(i3,1x)',advance='no') n_csf(j,ic)
!enddo;write(*,*)
!do j=1,no_electron
!  write(*,'(i3,1x)',advance='no') l_csf(j,ic)
!enddo;write(*,*)
!do j=1,no_electron
!  write(*,'(i3,1x)',advance='no') ml_csf(j,ic)
!enddo
!do j=1,no_electron
!  write(*,'(i2,1x)',advance='no') ms_csf(j,ic)
!enddo
 x=x+abs(coef_csf(ic)**2)
!write(*,'(f10.7,1x)') (coef_csf(ic))
 enddo
 write(*,*) 'Sum ^2 :',x
 write(*,*) 'No_det Xlm:',no_det
 if (abs(1.d0-x).gt.1d-6) call error(22,0,0)


!--- write out in awfn.data format
 do ic=1,no_det
  do j=1,no_up
   write(ouf1,'(4(1x,i2))') j      ,n_csf(j,ic),l_csf(j,ic),ml_csf(j,ic)
  enddo
  do j=no_up+1,no_electron
   write(ouf1,'(4(1x,i2))') j-no_up,n_csf(j,ic),l_csf(j,ic),ml_csf(j,ic)
  enddo
 enddo

!--- multipy coefficient by expansion coefficeient of CSF in MC expansion
 coef_csf(1:no_det)=coef_csf(1:no_det)*coef_expansion

!--- write out file portion to be added to input
 if(no_csf.eq.1)then
  write(ouf3,'(2x,e22.15,i6,1x,i3)') (coef_csf(ic),no_csf,0,ic=1,no_det)
 else
  write(ouf3,'(2x,e22.15,i6,1x,i3)') (coef_csf(ic),no_csf,1,ic=1,no_det)
 endif






 write(*,*) 'At end......'

 end subroutine extractcsf





 subroutine spec2slv(work,s,l,v)
 implicit none
 integer            :: i,j,k,s,l,v,n
 character(3)       :: work
 character          :: spec(0:9),spec_lower(0:9)
 character(2)       :: dig(-9:9)
 data spec/'S','P','D','F','G','H','I','K','L','M'/
 data spec_lower/'s','p','d','f','g','h','i','k','l','m'/
 data dig/'-9','-8','-7','-6','-5','-4','-3','-2','-1',         &
     &         ' 0',' 1',' 2',' 3',' 4',' 5',' 6',' 7',' 8',' 9'/

 do i=1,3
  do j=0,9
   if(work(i:i)==spec(j).or.work(i:i)==spec_lower(j)) then
    l=j
    k=i
   endif
  enddo
 enddo

 do i=0,9
   if(dig(i)(2:2)==work(k-1:k-1))s=i-1
   if(dig(i)(2:2)==work(k+1:k+1))v=i
 enddo
      
 end subroutine spec2slv


 subroutine spec2nl(work,n,l)
 implicit none
 integer            :: i,j,k,s,l,v,n
 character(3)       :: work
 character          :: spec(0:9),spec_lower(0:9)
 character(2)       :: dig(-9:9)
 data spec/'S','P','D','F','G','H','I','K','L','M'/
 data spec_lower/'s','p','d','f','g','h','i','k','l','m'/
 data dig/'-9','-8','-7','-6','-5','-4','-3','-2','-1',         &
     &         ' 0',' 1',' 2',' 3',' 4',' 5',' 6',' 7',' 8',' 9'/

 do i=1,3
  do j=0,9
   if(work(i:i)==spec_lower(j)) then
    l=j
    k=i
   endif
  enddo
 enddo

 do i=0,9
   if(dig(i)(2:2)==work(k-1:k-1))n=i
 enddo
      
 end subroutine spec2nl




end module extractcsf_mod

Program extractdet
!----------------------------------------------------------------------!
! ATSP2K to CASINO interface.                                          !
! Code to construct the multideterminat expansion of the               !
! configuration-state-function expansion provided by ATSP2K,           !
! implementing all the signs/ordering conventions correctly and        !
! consistently so signs of determinants within each CSF are consistent,!
! and that signs of CSF expansion coefficients are correct.            !
! Output is formatted for casino.                                      !
! J.R. Trail, TCM group, Cavendish Laboratory, Cambridge, UK 2013      !
! usage ~/bin/extractdet                                               !
! files required: wfn.out fort.60 cfg.inp                              !
! files provided: awfn.data correlation.data_foot                      !
!----------------------------------------------------------------------!
 use extractcsf_mod
 implicit none
 integer             :: nrmax,nod
 parameter (nrmax=10000,nod=220)
 integer             :: i,j,k,m,n,ii,jj,kk,nn
 integer             :: i1,i2,i3,i4,i5
 integer             :: i_csf,norbtot,niclosed,ncsf
 complex(kind(1.d0)) :: rm1
 real(kind(1.d0))    :: root2
 character           :: spec(0:9),spec_lower(0:9)
 character(2)        :: dig(-9:9)
 integer             :: nr,lo(30),no(30)
 real(kind(1.d0))    :: rmax,znuc
 character(6)        :: at,tt,atom,term
 character(3)        :: el1,new,clsd_shell(8),elc(8),couple(16),cterm
 character(24)       :: input,output1,output2,output3
 character(10)       :: lsym,l2sym
 character(11)       :: nsym
 character(60)       :: work
 integer             :: q(8)
 integer             :: n_dat_up(100,5),n_dat_down(100,5)
 integer             :: ouf1,ouf2,ouf3,iuf
 real(kind(1.d0))    :: p(nod,30),el(30),r(nod),r2(nod)
 real(kind(1.d0))    :: wfout(nrmax,30),rab(nrmax),rout(nrmax)
 integer             :: mm(30)
 integer             :: no_up,no_down,no_det
 integer             :: l_tot,s_tot,ml_tot,ms_tot
 data spec/'S','P','D','F','G','H','I','K','L','M'/
 data spec_lower/'s','p','d','f','g','h','i','k','l','m'/
 data dig/'-9','-8','-7','-6','-5','-4','-3','-2','-1',         &
     &         ' 0',' 1',' 2',' 3',' 4',' 5',' 6',' 7',' 8',' 9'/

 rm1=(0.d0,1.d0)
 root2=dsqrt(2.d0)
 nsym='0123456789 '
 lsym ='spdfghiklm'
 l2sym='SPDFGHIKLM'

 iuf=30
 ouf1=31
 ouf2=32
 ouf3=33





 input='wfn.out'
 call read_orbs(iuf,input,atom,r,r2,p,mm,no,lo,znuc,norbtot)
!--- orbitals now available on coarse grid of up to 220 points. not very fine.
 call regrid(nr,rout,rab,wfout,r,r2,p,norbtot,mm,no,lo,znuc)
!--- now on fine grid, rout, in array wfout.

!--- orbital n is in wfout(*,n),has n=no(n),l=lo(n), and there are norbtot
!---  of them.

!--- output to check them
!     do n=1,norbtot
!      do i=1,mm(n); write(100+n,*) r(i),p(i,n)*r2(i); enddo
!      do i=1,nr; write(200+n,*) rout(i),wfout(i,n); enddo
!---  enddo

!------------------------------------------------------------------------------
!---- All orbitals now available on fine grid
!------------------------------------------------------------------------------

!--- get TERM value from input file, and set ml_tot and ms_tot
 input='cfg.inp'
 open(unit=iuf,file=input,status='OLD',form='FORMATTED')
 read(iuf,*); read(iuf,*); read(iuf,*)
 read(iuf,'(16(1x,a3))')         (couple(j),j=1,16)
 close(unit=iuf)
 do j=1,16
  if(couple(j).ne.'   ')cterm=couple(j)
 enddo
 l_tot=scan(l2sym,cterm(2:2))-1
 s_tot=char_to_int(' '//cterm(1:1))-1

!--- Set Ml and Ms values here - these are used for all CSF
 ml_tot=l_tot
!write(*,*) 'L_tot',l_tot		!TMP
!write(*,*) 'ml_tot ?'			!TMP
!read(*,*)  ml_tot			!TMP
!--- Use ms_tot=s_tot for 1/4 in cusp condition (most of the time!)
!---   see Umrigar & Huang on cusps, about value of l_0 in their notation.
!---   IF ONE of the CSF is single det then l_0=0 definitely,
!---    and so 1/4 in cusp condition is correct.
!---   Choice of ms_tot just decides number of spin up/down electrons,
!---    within limits.
!---   eg C 3P, S=1, S_z   N_up   N_down
!---                  -1    2      4
!---                   0    3      3
!---                   1    4      2
!---   ms_tot=s_tot chosen so it looks like Hunds rule
 ms_tot=s_tot




!--- find out how many CSF there are
 open(unit=iuf,file=input,status='OLD',form='FORMATTED')
 i=0
 do
  i=i+1
  read(iuf,'(a)') work
  if (work(1:1).eq.'*') exit
 enddo
 close(unit=iuf)
 ncsf=(i-3)/2
!write(*,*) 'number CSFs :',ncsf

 output1='awfn.data'
 output3='correlation.data_foot'
 open(unit=ouf1,file=output1,form='FORMATTED')
 open(unit=ouf3,file=output3,form='FORMATTED')


 write(ouf1,'(3a)') 'Atomic ',atom(1:2),' wave function in real space'
 write(ouf1,'(a)') 'Atomic number'
 write(ouf1,'(1(i3))') nint(znuc)
 write(ouf1,'(a)') 'Total number of orbitals'
 write(ouf1,'(1(i3))') norbtot
 if (cterm(3:3).eq.' ') then
   write(ouf1,'(3a)') 'The ','['//cterm(1:2)//']',' state electronic configuration'
 else
   write(ouf1,'(3a)') 'The ','['//cterm(1:3)//']',' state electronic configuration'
 endif
 write(ouf1,'(a)') 'Number of up, downspin electrons, determinants '
 write(ouf1,'(3(1x,i2,1x,i2,1x,i5))') 0,0,0 ! These must be rewritten later
 write(ouf1,'(a)') 'States'


 write(ouf3,'(a)') ' START MDET'
 write(ouf3,'(a)') '  Title'
 write(ouf3,'(4a)') '  ',atom(1:2),' atom (numerical orbitals)'
 write(ouf3,'(a)') '  Type'
 write(ouf3,'(a)') '   MD'
 write(ouf3,'(a)') '   00'                  ! rewritten file later



!----- write out the determinants
 no_det=0
 do i_csf=1,ncsf
  write(*,*) 'Doing CSF no:',i_csf
  call extractcsf(i_csf,iuf,ouf1,ouf2,ouf3,l_tot,s_tot,ml_tot,  &
          &       ms_tot,no_up,no_down,nn)
  no_det=no_det+nn
 enddo

 write(*,*) 'Total determinants:',no_det


!----- write out the tabulated orbitals
 write(ouf1,'(a)') 'Radial grid (a.u.)'
 write(ouf1,'(i5)') nr
 do j=1,nr
   write(ouf1,'(e25.15)') rout(j)
 enddo
 do i=1,norbtot
  write(ouf1,'(a,i2,a,i2,a,a)')                                &
 &    'Orbital # ',i,' [',no(i),lsym(lo(i)+1:lo(i)+1),']'
  write(ouf1,'(3(1x,i2))') 0,no(i),lo(i)
  do j=1,nr
    write(ouf1,'(e25.15)') wfout(j,i)
  enddo
 enddo





 write(ouf3,'(a)') ' END MDET'




 close(unit=ouf1)
 close(unit=ouf3)


!--- Put correct no_up,no_down,no_det values in ouf1
 open(unit=ouf1,file=output1,access='direct',recl=1,form='FORMATTED')
 i1=(no_up)/10; i2=(no_up-i1*10)
 if(no_up.ge.10)write(ouf1,'(a1)',rec=175)  nsym(i1+1:i1+1)
                write(ouf1,'(a1)',rec=176)  nsym(i2+1:i2+1)

 i1=(no_down)/10; i2=(no_down-i1*10)
 if(no_up.ge.10)write(ouf1,'(a1)',rec=178)  nsym(i1+1:i1+1)
                write(ouf1,'(a1)',rec=179)  nsym(i2+1:i2+1)

 j =no_det
 i1=j/10000         ;j=j-i1*10000
 i2=j/1000          ;j=j-i2*1000
 i3=j/100           ;j=j-i3*100
 i4=j/10            ;j=j-i4*10
 i5=j
 if(no_det.ge.10000)write(ouf1,'(a1)',rec=181)  nsym(i1+1:i1+1)
 if(no_det.ge.1000 )write(ouf1,'(a1)',rec=182)  nsym(i2+1:i2+1)
 if(no_det.ge.100  )write(ouf1,'(a1)',rec=183)  nsym(i3+1:i3+1)
 if(no_det.ge.10   )write(ouf1,'(a1)',rec=184)  nsym(i4+1:i4+1)
                    write(ouf1,'(a1)',rec=185)  nsym(i5+1:i5+1)
 close(unit=ouf1)


!--- Put correct no_det value in ouf3
 open(unit=ouf3,file=output3,access='direct',recl=1,form='FORMATTED')

 j =no_det
 i1=j/10000         ;j=j-i1*10000
 i2=j/1000          ;j=j-i2*1000
 i3=j/100           ;j=j-i3*100
 i4=j/10            ;j=j-i4*10
 i5=j
 if(no_det.ge.10000)write(ouf3,'(a1)',rec=65)  nsym(i1+1:i1+1)
 if(no_det.ge.1000 )write(ouf3,'(a1)',rec=66)  nsym(i2+1:i2+1)
 if(no_det.ge.100  )write(ouf3,'(a1)',rec=67)  nsym(i3+1:i3+1)
 if(no_det.ge.10   )write(ouf3,'(a1)',rec=68)  nsym(i4+1:i4+1)
                    write(ouf3,'(a1)',rec=69)  nsym(i5+1:i5+1)
 close(unit=ouf3)

contains







!--- Read in orbitals and data from file input on iuf
 subroutine read_orbs(iuf,input,atom,r,r2,p,mm,no,lo,znuc,norb)
 implicit none
 integer           :: nrmax,nod
 parameter (nrmax=10000,nod=220)
 integer           :: iuf
 character(24)     :: input
 character(6)      :: atom
 real(kind(1.d0))  :: r(:),r2(:),p(:,:)
 integer           :: mm(:),no(:),lo(:)
 real(kind(1.d0))  :: znuc
 integer           :: norb
 character(6)      :: at,tt,term
 character(3)      :: el1
 integer           :: m,nr,j
 real(kind(1.d0))  :: zt,eti,eki,azi,rho,h

 open(unit=iuf,file=input,status='OLD',form='UNFORMATTED')

 norb=1

 do
 p(1,norb) = 0.d0
  read(iuf,END=5) at,tt,el1,m,zt,eti,eki,azi,(p(j,norb),j=2,m+1)
  mm(norb)=m
  write(*,*) at                    ! name of atom
  write(*,*) tt                    ! TERM
  write(*,*) el1                   ! nl
  write(*,*) m                     ! no. points recorded
  write(*,*) zt                    ! atomic no.
  write(*,*) eti                   ! energy
  write(*,*) eki                   ! kinetic energy
  write(*,*) azi                   ! AZ (see output for above for MCHF)
  write(*,*)
!- jrt above report not needed

  no(norb)=char_to_int(el1(1:2))         ! store n value of orbital
  lo(norb)=spec_to_l  (el1(3:3))         ! store l value of orbital
  znuc=zt                                ! atomic number

  norb=norb+1
 enddo

5 close(unit=iuf)

 norb=norb-1
 atom=at

!----- define grid used in Atsp2K progs.
 rho  = -4.d0
 h    = 1.d0/16.d0
 r(1) = 0.d0
 r2(1)= 0.d0
 do j=2,220
  r(j) =exp(rho)/znuc
  r2(j)=sqrt(r(j))
  rho  =rho + h
 enddo

 end subroutine read_orbs








!------ Regrid orbitals and put into array wfout on rout grid.
 subroutine regrid(nr,rout,rab,wfout,r,r2,p,              &
 &                  norbtot,mm,no,lo,znuc)
 implicit none
 integer           :: nrmax,nod
 parameter (nrmax=10000,nod=220)
 integer           :: nr,norbtot,mm(:),lo(:),no(:)
 real(kind(1.d0))  :: rout(:),rab(:),wfout(:,:)
 real(kind(1.d0))  :: r(nod),r2(nod),p(:,:),znuc
 real(kind(1.d0))  :: rmax
 real(kind(1.d0))  :: wf(nrmax),wfd(nrmax),wfdd(nrmax)
 real(kind(1.d0))  :: aa,bb,a,b
 integer           :: i,n

!--- define the new grid
 aa  =  8.d0
 bb  = 75.d0
 rmax=100.d0
 if (aa.eq.0.D0) aa= 5.d0
 if (bb.eq.0.D0) bb=40.d0
 if (rmax.le.0.d0) rmax=80.d0
 a = bb*exp(-aa*log(10.d0))/znuc
 b = 1/bb

 do i=1,nrmax
   rout(i) = a * (exp(b*(i-1)) - 1)
   rab(i) = (rout(i) + a) * b                 ! dr/di=b(r+a)
   nr = i
   if (rout(i).ge.rmax) exit
 enddo

!--- write(*,*) znuc,aa,bb,nrmax,r(nrmax)

 do n=1,norbtot

  do i=1,mm(n)
    wf(i)=p(i,n)*r2(i)
  enddo
  do i=mm(n)+1,nod
    wf(i)=0.d0
  enddo
  call newgrid_wf(r,wf,nod,rout,wfout(:,n),wfd,wfdd,nr,lo(n),znuc)

 enddo

 end subroutine regrid






 subroutine newgrid_wf(rold,vold,nold,rnew,vnew,vd,vdd,nnew,l,znuc)
 implicit none
 integer           :: i,igr,is,it,iel,nr
 integer           :: nold,nnew,npoly,mpoly,igrid
 parameter (mpoly=15)
 real(kind(1.d0))  :: dy,y,yd,ydd,x
 real(kind(1.d0))  :: a0,a1,znuc
 real(kind(1.d0))  :: rnew(*),rold(*),vold(*),vnew(*)
 real(kind(1.d0))  :: xa(mpoly),ya(mpoly)
 real(kind(1.d0))  :: vd(*),vdd(*),r,ldp
 integer           :: l,p
!---  Evaluate a wavefunction  on a new grid of points
!---  The old grid rold contains nold points and the old potential is vold
!---  The new grid rnew contains nnew points and the new potential is vnew


 npoly=8+max(0,l-2)
 ldp=dble(l)

!--- wf is described as  1   .. r**(npoly-1) for r>r(1) - npoly     points
!--- wf is described as r**l .. r**(npoly-1) for r<r(1) - npoly-l   points
!---  with r**(l+1) term chosen to cancel out coulomb potential.
!---  NOT CORRECT FOR PSEUDO-ORBITALS!!!!!!!!!!!!!!!!!!!!!!
!---      npoly => 4+max(0,l-2)
 do igrid=1,nnew
   r  = rnew(igrid)
   call locate(rold,nold,r,igr)
   igr=min(max(igr-(npoly-1)/2,1),nold+1-npoly)
   if (igr.ne.1) then
     xa(1:npoly) = rold(igr:igr+npoly-1)
     ya(1:npoly) = vold(igr:igr+npoly-1)/rold(igr:igr+npoly-1)
     call polintd(xa(1),ya(1),npoly  ,r,y,yd,ydd)
     vnew(igrid) = r*y
     vd  (igrid) = r*yd
     vdd (igrid) = r*ydd
   else
     xa(2:npoly-l+1) = rold(2:npoly-l+1)
     ya(2:npoly-l+1) = vold(2:npoly-l+1)/rold(2:npoly-l+1)**(l+1)
     call polintd(xa(2),ya(2),npoly-l,0.d0,y,yd,ydd)
     a0=y;a1=-a0*znuc/dble(l+1)
     ya(2:npoly-l-1)=ya(2:npoly-l-1)-a0-a1*rold(2:npoly-l-1)
     ya(2:npoly-l-1)=ya(2:npoly-l-1)/rold(2:npoly-l-1)**2
     call polintd(xa(2),ya(2),npoly-l-2,r,y,yd,ydd)
     call polintd(xa(2),ya(2),npoly-l-2,r,y,yd,ydd)
     ydd= 2.d0*y+4.d0*r*yd+r**2*ydd
     yd = 2.d0*r*y+r**2*yd          +a1
     y  = r**2*y                    +a1*r+a0
!----- v this would be appropriate instead of above if pseudo-orbitals
!-----    call polintd(xa(2),ya(2),npoly-l,r,y,yd,ydd)
!----- ^ this would be appropriate instead of above if pseudo-orbitals
     if (l.eq.0) then
      vnew(igrid) = r*y
      vd  (igrid) = r*yd
      vdd (igrid) = r*ydd
     else if (l.eq.1) then
      vnew(igrid) = r**2*y
      vd  (igrid) = r**2*yd+r*y
      vdd (igrid) = r**2*ydd+2.d0*r*yd
     else if (l.ge.2) then
      vnew(igrid) = r**(l+1)*y
      vd  (igrid) = r**(l+1)*yd+ldp*r**l*y
      vdd (igrid) = r**(l+1)*ydd+2.d0*ldp*r**l*yd+             &
 &                   ldp*(ldp-1.d0)*r**(l-1)*y
     endif
   endif
 enddo

 end subroutine newgrid_wf





!---  *****************************************************************
!---  *****************************************************************
!---  **                                                             **
!---  ** From Numerical Recipies (2nd ed.), page 111. This routine   **
!---  ** works as follows. Given an array XX(1:N) and given a value  **
!---  ** X, it returns a value J such that X is between XX(J) and    **
!---  ** XX(J+1). XX(1:N) must be monotonic increasing or decreasing **
!---  ** and J=0 or J=N is returned to indicate that X is out of     **
!---  ** range.                                                      **
!---  **                                                             **
!---  *****************************************************************
!---  *****************************************************************
 subroutine locate(xx,n,x,j)
 implicit none
 real(kind(1.d0))  :: x
 integer           :: j,n
 real(kind(1.d0))  :: xx(n)
 integer           :: jl,jm,ju
!---------------------------------------------------------------------

!---  Initialize lower and upper limits.
 jl = 0
 ju = n + 1

!---  If not yet done, compute a midpoint and replace either the lower
!---  or upper limit as appropriate.
 do
  if (ju-jl.le.1) exit
  jm = (ju+jl)/2
  if ((xx(n).gt.xx(1)) .eqv. (x.gt.xx(jm))) then
    jl = jm
  else
    ju = jm
  end if
!--- Repeat until test condition is satisfied.
 enddo
 j = jl
 
 end subroutine locate






!---------------------------------------------------------------------!
! This routine is a generalisation of Numerical Recipies              !
! (2nd Ed.) page 103 written by JRT (2005) which uses                 !
! Neville's algorithm to accurately make a polynomial                 !
! interpolation through N points.                                     !
! The generalisation is that this routine also provides the           !
! first and second derivatives of the interpolating                   !
! polynomial without calculating the coefficients of the              !
! polynomial. Prevents numerical errors, as the original              !
! Neville's did for interpolation alone.                              !
!---------------------------------------------------------------------!
 subroutine polintd(xa,ya,n,x,y,y_d,y_dd)
 implicit none
 integer nmax
 parameter (nmax=10)
 real(kind(1.d0))  :: dy,dy_d,dy_dd,x,y,y_d,y_dd
 integer n
 real(kind(1.d0))  :: xa(n),ya(n)
 real(kind(1.d0))  :: den,dif,dift,ho,hp,w,w_d,w_dd
 integer i,m,ns
 real(kind(1.d0))  :: c   (nmax),d   (nmax)
 real(kind(1.d0))  :: c_d (nmax),d_d (nmax)
 real(kind(1.d0))  :: c_dd(nmax),d_dd(nmax)
 intrinsic abs
!----------------------------------------------------------------------

 ns = 1
 dif = abs(x-xa(1))

!--- Here we find the index ns, of the closest table entry.

 do i = 1,n

   dift = abs(x-xa(I))
   if (dift.lt.dif) then
     ns  = i
     dif = dift
   endif
!--- Initializing the tableau of entries of c's and d's.
   c   (i) = ya(i)
   d   (i) = ya(i)
   c_d (i) = 0.d0
   d_d (i) = 0.d0
   c_dd(i) = 0.d0
   d_dd(i) = 0.d0

 enddo


!---   This is the initial approximation to y.
 y   = ya(ns)
 y_d = 0.d0
 y_dd= 0.d0
 ns  = ns - 1

!---  For each column of the tableu, loop over the current c's and d's
!---  and update them.

 do m = 1,n - 1

   do i = 1,n - m
     ho   = xa(i) - x
     hp   = xa(i+m) - x
     w    = c   (i+1) - d   (i)
     w_d  = c_d (i+1) - d_d (i)
     w_dd = c_dd(i+1) - d_dd(i)
     den  = ho - hp
!---  Here the c's and d's are updated.
     d   (i) = hp*w/den
     c   (i) = ho*w/den
     d_d (i) = (-w + hp*w_d)/den
     c_d (i) = (-w + ho*w_d)/den
     d_dd(i) = (-2.d0*w_d +hp* w_dd)/den
     c_dd(i) = (-2.d0*w_d +ho* w_dd)/den
   enddo

!---  After each column in the tableau is completed, decide which
!---  correction, c or d needs to be added to the accumulating value
!---  of y i.e. which path to take through the tableau-forking  up or
!---  down. This is done in such a way as to take the most "stright-line"
!---  route through the tableau to its apex, updating ns accordingly, to
!---  keep track of where we are. This route keeps the partial approx
!---  centered on the target x. The last dy added is the error indication.

   if (2*ns.lt.n-m) then
     dy   = c   (ns+1)
     dy_d = c_d (ns+1)
     dy_dd= c_dd(ns+1)
   else
     dy    = d   (ns)
     dy_d  = d_d (ns)
     dy_dd = d_dd(ns)
     ns = ns - 1
   endif

   y    = y    + dy
   y_d  = y_d  + dy_d
   y_dd = y_dd + dy_dd

 enddo

 end subroutine polintd





!--- Convert 1 character string of spd to a l value
 integer function spec_to_l(string)
 implicit none
 character(10)    :: lsym
 character(1)     :: string

 lsym='spdfghiklm'
 spec_to_l=scan(lsym,string)-1

 end function spec_to_l


!--- Convert 2 character string of digits to an integer
 integer function char_to_int(string)
 implicit none
 character(11)    :: nsym
 character(2)     :: string
 integer          :: i1,i2

 nsym='0123456789 '
 i2=scan(nsym,string(1:1))-1
 i1=scan(nsym,string(2:2))-1
 if(i2.eq.10)i2=0
 char_to_int=10*i2+i1

 end function char_to_int






end program extractdet

!**************************************************************************
!                                                                         *
! PTM                                                                     *
! Perform various transformations on QMC pseudopotentials.                *
! Mike Towler Sep 1999                                                    *
!                                                                         *
! PTM will read in the QMC x_pp.data format, or Gaussian expansions       *
! of such functions in the format documented in the accompanying          *
! README file.                                                            *
!                                                                         *
! PTM cam perform the following transformations on x_pp.data files:       *
! (1) Multiply or divide data by a power of r                             *
! (2) Calculate differences between different components                  *
! (3) Convert to xmgr format                                              *
!                                                                         *
! and the following transformations on the fit file                       *
! (1) perform the Gaussian expansion and plot results on radial grid      *
! with options to subtract Z/r and to convert between Hartree and Rydbergs*
! The result can be plotted as x_pp.data or xmgr plot.                    *
!                                                                         *
! NB: Currently doesn't handle core polarization potentials which can     *
! live at the bottom of the x_pp.data file.                               *
!                                                                         *
!**************************************************************************
      MODULE transform_defs
      INTEGER ialloc,rpower,rnpower
      CHARACTER(72) filename
      END MODULE transform_defs

      PROGRAM pp_transform
      USE transform_defs
      IMPLICIT NONE
      INTEGER i
      LOGICAL ex
      INTERFACE
       SUBROUTINE pp1
       END SUBROUTINE pp1
       SUBROUTINE pp2
       END SUBROUTINE pp2
      END INTERFACE
      write(6,*)'PTM'
      write(6,*)'==='
10    write(6,*)'choose input type:'
      write(6,*)&
     &'(1) Up to 3 pp components on a radial grid in x_pp.data format'
      write(6,*)'(2) gaussian exponents, coefficients and powers of r'
      read(5,*)i
      if(i/=1.and.i/=2)goto 10
      write(6,*)'Name of input file?'
      read(5,'(a)')filename
      inquire(file=filename,exist=ex)
      if(.not.ex)then
       write(6,*)'This file does not exist'
       stop
      endif
      if(i==1)then
       call pp1 ! manipulate pp on grid
      else
       call pp2 ! manipulate gaussian expansion of pp
      endif
      END PROGRAM pp_transform

      SUBROUTINE pp1 ! grid
      USE transform_defs
      IMPLICIT NONE
      INTEGER i,n,ndata,nx,ny,ncomp
      LOGICAL xmgr,again
      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: r,ys,yp,yd,yf,yg
      CHARACTER(1) a
      CHARACTER(4) rlabel(2)
      CHARACTER(20) label1,label2,label3,label4,label5
      CHARACTER(72) line(20)
      open(20,file=trim(filename),status='old')
13    write(6,*)&
     & 'What power of r is this data already multiplied by (0,1,2)?'
      read(5,*)rpower
      if(rpower/=0.and.rpower/=1.and.rpower/=2)goto 13
      rlabel(1)='  r*'
      rlabel(2)='r^2*'
      xmgr=.false.
      again=.true.

14    write(6,*)'How many pp components are there (1,2,3,4,5)?'
      read(5,*)ncomp
      if(ncomp/=1.and.ncomp/=2.and.ncomp/=3.and.ncomp/=4.and.ncomp/=5)goto 14

! Read data file
      read(20,'(a)',err=15,end=15)line(1)
      read(20,'(a)',err=15,end=15)line(2)
      read(20,'(a)',err=15,end=15)line(3)
      read(20,'(a)',err=15,end=15)line(4)
      read(20,'(a)',err=15,end=15)line(5)
      read(20,'(a)',err=15,end=15)line(6)
      read(20,'(a)',err=15,end=15)line(7)
      read(20,'(a)',err=15,end=15)line(8)
      read(20,'(a)',err=15,end=15)line(9)
      read(20,'(a)',err=15,end=15)line(10)
      read(20,*,err=15,end=15)ndata
      ialloc=0
      allocate(r(ndata),ys(ndata),yp(ndata),yd(ndata),yf(ndata),yg(ndata),&
       &stat=ialloc)
      if(ialloc/=0)then
       write(6,*)'Allocation error.'
       stop
      endif ! ialloc/=0
      read(20,'(a)',err=15,end=15)line(11)
      do i=1,ndata
       read(20,*,err=15,end=15)r(i)
      enddo
      read(20,*,err=15,end=15)line(12)
      do i=1,ndata
       read(20,*,err=15,end=15)ys(i)
      enddo
      if(ncomp>1)then
       read(20,*,err=15,end=15)line(13)
       do i=1,ndata
        read(20,*,err=15,end=15)yp(i)
       enddo
      endif
      if(ncomp>2)then
       read(20,*,err=15,end=15)line(14)
       do i=1,ndata
        read(20,*,err=15)yd(i)
       enddo
      endif
      if(ncomp>3)then
       read(20,*,err=15,end=15)line(15)
       do i=1,ndata
        read(20,*,err=15)yf(i)
       enddo
      endif
      if(ncomp>4)then
       read(20,*,err=15,end=15)line(15)
       do i=1,ndata
        read(20,*,err=15)yg(i)
       enddo
      endif      
      goto 30
15    write(6,*)'Problem reading file.'
      stop
! What do you want me to do with this?

30    write(6,*)'File read in successfully.'
40    write(6,*)'What do you want me to do with it?'
      write(6,*)'(1) Multiply or divide data by a power of r'
      write(6,*)'(2) Calculate differences between different components'
      write(6,*)'(3) Write it out as an xmgr plot'
      read(5,*)i
      label1='potential (L=0) in Ry    '
      label2='potential (L=1) in Ry    '
      label3='potential (L=2) in Ry    '
      label4='potential (L=3) in Ry    '
      label5='potential (L=4) in Ry    '

      select case (i)

      case(1) ! powers of r
       write(6,*)'Power (+ve or -ve) to multiply by:'
       read(5,*)n
       if(n>2.or.n<-2)then
        write(6,*)'Power out of range -2 to +2'
        stop
       endif
       if(rpower+n<0)then
        write(6,*)'This power not compatible with existing rpower'
        write(6,*)'(See program if you don''t understand this.)'
        stop
       endif
       if(n<0)then
        ys(2:ndata)=ys(2:ndata)*r(2:ndata)**n
        if(ncomp>1)yp(2:ndata)=yp(2:ndata)*r(2:ndata)**n
        if(ncomp>2)yd(2:ndata)=yd(2:ndata)*r(2:ndata)**n
        if(ncomp>3)yf(2:ndata)=yf(2:ndata)*r(2:ndata)**n
        if(ncomp>4)yg(2:ndata)=yg(2:ndata)*r(2:ndata)**n
        ys(1)=ys(2)
        if(ncomp>1)yp(1)=yp(2)
        if(ncomp>2)yd(1)=yd(2)
        if(ncomp>3)yf(1)=yf(2)
        if(ncomp>4)yg(1)=yg(2)
       else
        ys=ys*r**n
        if(ncomp>1)yp=yp*r**n
        if(ncomp>2)yd=yd*r**n
        if(ncomp>3)yf=yf*r**n
        if(ncomp>4)yg=yg*r**n
       endif
45     write(6,*)'Write output in:'
       write(6,*)'(1) x_pp.data format for QMC'
       write(6,*)'(2) xmgr format for plotting'
       read(5,*)i
       if(i/=1.and.i/=2)goto 45
       if(i==2)xmgr=.true.

      case(2) ! difference
       if(ncomp<2)then
        write(6,*)&
     &  'I need at least 2 components to calculate differences.'
        stop
       endif
47     write(6,*)'I will calculate Vx-Vy'
50     write(6,*)'x corresponds to (1,2,3)?'
       read(5,*)nx
       if(nx/=1.and.nx/=2.and.nx/=3)goto 50
60     write(6,*)'y corresponds to (1,2,3)?'
       read(5,*)ny
       if(ny/=1.and.ny/=2.and.ny/=3)goto 60
       if(nx==ny)then
        write(6,*)'nx cannot equal ny'
        goto 50
       endif
       if(nx==3.or.ny==3)then
        if(ncomp==2)&
     &  write(6,*)'3rd component requested. I only have two.'
       endif
       select case(nx)
       case(1)
        select case(ny)
        case(2)
         ys=ys-yp
         label1='Vs-Vp'
        case(3)
         ys=ys-yd
         label1='Vs-Vd'
        end select
       case(2)
        select case(ny)
        case(1)
         yp=yp-ys
         label2='Vp-Vs'
        case(3)
         yp=yp-yd
         label2='Vp-Vd'
        end select
       case(3)
        select case(ny)
        case(1)
         yd=yd-ys
         label3='Vd-Vs'
        case(2)
         yd=yd-yp
         label3='Vd-Vp'
        end select
       end select
       if(again)then
        write(6,*)'Do you wish to calculate another difference?'
        read(5,'(a)')a
        if(a=='y')then
         again=.false.
         goto 47
        endif
       endif
70     write(6,*)'Write output in:'
       write(6,*)'(1) x_pp.data format for QMC'
       write(6,*)'(2) xmgr format for plotting'
       read(5,*)i
       if(i/=1.and.i/=2)goto 70
       if(i==2)xmgr=.true.
      case(3) ! xmgr plot
       xmgr=.true.
       goto 100
      case default
       write(6,*)
       goto 40
      end select


!.....print out the transformed data
100    open(21,file='PP_out',status='unknown')
       if(xmgr)then ! write out in xmgr format
       do i=1,ndata
        write(21,*)r(i),ys(i)
       enddo
       if(ncomp>1)then
        write(21,*)'&'
        do i=1,ndata
         write(21,*)r(i),yp(i)
        enddo
       endif
       if(ncomp>2)then
        write(21,*)'&'
        do i=1,ndata
         write(21,*)r(i),yd(i)
        enddo
       endif
       if(ncomp>3)then
        write(21,*)'&'
        do i=1,ndata
         write(21,*)r(i),yf(i)
        enddo
       endif
       if(ncomp>4)then
        write(21,*)'&'
        do i=1,ndata
         write(21,*)r(i),yg(i)
        enddo
       endif        
      else ! write out in QMC x_pp.data format
       rnpower=rpower+n
       write(21,*)line(1)
       write(21,*)line(2)
       write(21,*)line(3)
       write(21,*)line(4)
       write(21,*)line(5)
       write(21,*)line(6)
       write(21,*)line(7)
       write(21,*)line(8)
       write(21,*)line(9)
       write(21,*)line(10)
       write(21,*)ndata
       write(21,*)line(11)
       do i=1,ndata
        write(21,*)r(i)
       enddo
       if(rnpower>0)then
        write(21,*)rlabel(rnpower),label1
       else
        write(21,*)label1
       endif
       do i=1,ndata
        write(21,*)ys(i)
       enddo
       if(ncomp>1)then
        if(rnpower>0)then
         write(21,*)rlabel(rnpower),label2
        else
         write(21,*)label2
        endif
        do i=1,ndata
         write(21,*)yp(i)
        enddo
       endif
       if(ncomp>2)then
        if(rnpower>0)then
         write(21,*)rlabel(rnpower),label3
        else
         write(21,*)label3
        endif
        do i=1,ndata
         write(21,*)yd(i)
        enddo
       endif
      if(ncomp>3)then
        if(rnpower>0)then
         write(21,*)rlabel(rnpower),label4
        else
         write(21,*)label4
        endif
        do i=1,ndata
         write(21,*)yf(i)
        enddo
       endif
      if(ncomp>4)then
        if(rnpower>0)then
         write(21,*)rlabel(rnpower),label5
        else
         write(21,*)label5
        endif
        do i=1,ndata
         write(21,*)yg(i)
        enddo
       endif
      endif
      write(6,*)'Transformed data written to PP_out'
      stop
      END SUBROUTINE pp1


      SUBROUTINE pp2 ! gaussian expansion
! NB: It is implicitly assumed that you feed this routine the pure pseudo
! (i.e. not multiplied by r)
      USE transform_defs
      IMPLICIT NONE
      INTEGER maxcomp,ngridpoints,maxg
      PARAMETER(maxcomp=5,ngridpoints=3000,maxg=100)
      INTEGER i,j,n,ncomp,np(maxg),nk(maxcomp,maxg),nr,coeffs
      DOUBLE PRECISION e(maxcomp,maxg),c(maxcomp,maxg),znuc,&
     & rmax,aa,bb,a,b,r(ngridpoints),pp(maxcomp,ngridpoints),&
     & atnum,conv,temp1(maxcomp,maxg),temp2(maxcomp,maxg)
      CHARACTER(1) yes
      LOGICAL xmgr,coulomb(5),rydberg
      xmgr=.false.
      coulomb=.false.

! Read in the pseudo with following format:

!     Nickel Durand            ! Title
!     28.0 10.0                ! atomic number, effective nuclear charge
!     3                        ! no of components of the pseudo (s,p,d etc.)
!     2 1 4                    ! (no of Gaussians describing each component)
!     P0                       ! label
!     1.02153800     10.55013200  -1  !  (exponent, coefficient, power of r)
!     1.02153800     -1.01208800   2
!     P1
!     1.40000000     16.41825700   0
!     P2
!     2.73853000     -0.58105000  -2
!     2.73853000     -6.35383400  -1
!     2.73853000      6.65505900   0
!     2.73853000    -11.00888000   2

      open(20,file=trim(filename),status='old')

      read(20,*)
! Atomic number, effective nuclear charge
      read(20,*)atnum,znuc
! No of components
      read(20,*)ncomp
! No of gaussians in each component
      read(20,*)(np(i),i=1,ncomp)
! Read the Gaussians, exponents and powers of r.
      do i=1,ncomp
       read(20,*) ! label
       do j=1,np(i)
! Exponents, coefficients, power of r
         read (20,*)temp1(i,j),temp2(i,j),nk(i,j)
       enddo
      enddo

      coeffs=0
      if((any(temp1(:,:)<0.d0)).and.(any(temp2(:,:)<0.d0)))then
       write(6,*)'Negative exponents and coefficients?'
       write(6,*)'Don''t think so.'
       stop
      endif
      if(any(temp1(:,:)<0.d0))coeffs=2 ! Gaussian format
      if(any(temp2(:,:)<0.d0))coeffs=1 ! CRYSTAL format
      if(coeffs==0)then
1      write(6,*)'Data is in :'
       write(6,*)&
     & '(1) CRYSTAL95 format (alpha, coefficient, power of r)?'
       write(6,*)'(2) GAUSSIAN format (coefficient, alpha, power of r)?'
       read(5,*)coeffs
       if(coeffs/=1.and.coeffs/=2)goto 1
      endif

      if(coeffs==1)then ! CRYSTAL
       e=temp1
       c=temp2
      else ! GAUSSIAN
       e=temp2
       c=temp1
      endif

! Construct the grid.
      rmax=50.d0
      aa=8.d0
      bb=80.d0
      a=bb*exp(-aa*log(10.d0))/znuc
      b=1/bb
      do i=1,ngridpoints
       r(i)=a*(exp(b*(i-1))-1)
       nr=i
       if(r(i)>rmax)goto 200
      enddo

! Construct the pseudopotential components at each point on the grid.
200   do n=2,nr
       do i=1,ncomp
        pp(i,n)=0.d0
        do j=1,np(i)
         pp(i,n)=pp(i,n)+&
     & r(n)**dble(nk(i,j))*c(i,j)*dexp(-(e(i,j)*r(n)**2.d0))
        enddo
       enddo
      enddo
      do i=1,ncomp
       pp(i,1)=pp(i,2)
      enddo

      do i=1,ncomp
       write(6,*)&
     &  'Do you want to subtract Z/r off the result of the expansion'
       write(6,*)'for component ',i,'? (you do if it''s local)'
       read(5,'(a)')yes
       if(yes=='y')coulomb(i)=.true.
      enddo

300   write(6,*)'Write output in:'
      write(6,*)'(1) x_pp.data format for QMC'
      write(6,*)'(2) xmgr format for plotting'
      read(5,*)i
      if(i/=1.and.i/=2)goto 300
      if(i==2)xmgr=.true.

400   write(6,*)'Want output in:'
      write(6,*)'(1) Hartree (i.e. leave it alone)'
      write(6,*)'(2) Rydberg (i.e. multiply current data by 2)'
      read(5,*)i
      if(i/=1.and.i/=2)goto 400
      if(i==1)then
       conv=1.d0
       rydberg=.false.
      endif
      if(i==2)then
       conv=2.d0
       rydberg=.true.
      endif

      open(21,file='PP_out',status='unknown')
      if(xmgr)then

       do i=1,nr
        if(coulomb(1))then
         if(i/=1)then
          write(21,*)r(i),conv*(pp(1,i)-znuc*r(i)**(-1.d0))
         else
          write(21,*)r(1),conv*(pp(1,2)-znuc*r(2)**(-1.d0))
         endif
        else
         write(21,*)r(i),conv*pp(1,i)
        endif
       enddo
       if(ncomp>1)then
        write(21,*)'&'
        do i=1,nr
         if(coulomb(2))then
          if(i/=1)then
           write(21,*)r(i),conv*(pp(2,i)-znuc*r(i)**(-1.d0))
          else
           write(21,*)r(1),conv*(pp(2,2)-znuc*r(2)**(-1.d0))
          endif
         else
          write(21,*)r(i),conv*pp(2,i)
         endif
        enddo
       endif
       if(ncomp>2)then
        write(21,*)'&'
        do i=1,nr
         if(coulomb(3))then
          if(i/=1)then
           write(21,*)r(i),conv*(pp(3,i)-znuc*r(i)**(-1.d0))
          else
           write(21,*)r(1),conv*(pp(3,2)-znuc*r(2)**(-1.d0))
          endif
         else
          write(21,*)r(i),conv*pp(3,i)
         endif
        enddo
       endif
      if(ncomp>3)then
        write(21,*)'&'
        do i=1,nr
         if(coulomb(4))then
          if(i/=1)then
           write(21,*)r(i),conv*(pp(4,i)-znuc*r(i)**(-1.d0))
          else
           write(21,*)r(1),conv*(pp(4,2)-znuc*r(2)**(-1.d0))
          endif
         else
          write(21,*)r(i),conv*pp(4,i)
         endif
        enddo
       endif
      if(ncomp>4)then
        write(21,*)'&'
        do i=1,nr
         if(coulomb(5))then
          if(i/=1)then
           write(21,*)r(i),conv*(pp(5,i)-znuc*r(i)**(-1.d0))
          else
           write(21,*)r(1),conv*(pp(5,2)-znuc*r(2)**(-1.d0))
          endif
         else
          write(21,*)r(i),conv*pp(5,i)
         endif
        enddo
       endif

      else
       write(21,*)'Pseudopotential Title goes here'
       write(21,*)'Atomic number and pseudo-charge'
       write(21,'(f5.1,1x,f5.1)')atnum,znuc
       write(21,*)'Energy units (rydberg/hartree/ev):'
       if(rydberg)then
        write(21,*)'rydberg'
       else
        write(21,*)'hartree'
       endif
       write(21,*)'Angular momentum of local component (0=s,1=p,2=d..) for &
        &DFT and QMC- EDIT IF NECESSARY'
       write(21,*)'0 2'
       write(21,*)&
     &'NLRULE override (1) VMC/DMC (2) config gen ',&
     &'(0 ==> input/default value)'
       write(21,*)'0 0'
       write(21,*)'Number of grid points'
       write(21,*)nr
       write(21,*)'R(i) in atomic units'
       do n=1,nr
        write(21,*)r(n)
       enddo
       write(21,*)'r*potential (L=0) in Ry '
       write(21,*)'0.000 0.000'
       do n=2,nr
        if(coulomb(1))then
         if(n/=1)then
          write(21,*)conv*(pp(1,n)-znuc*r(n)**(-1.d0))
         else
          write(21,*)conv*(pp(1,2)-znuc*r(2)**(-1.d0))
         endif
        else
         write(21,*)conv*pp(1,n)
        endif
       enddo
       if(ncomp>1)then
        write(21,*)'r*potential (L=1) in Ry '
        do n=1,nr
         if(coulomb(2))then
          if(n/=1)then
           write(21,*)conv*(pp(2,n)-znuc*r(n)**(-1.d0))
          else
           write(21,*)conv*(pp(2,2)-znuc*r(2)**(-1.d0))
          endif
         else
          write(21,*)conv*pp(2,n)
         endif
        enddo
       endif
       if(ncomp>2)then
        write(21,*)'r*potential (L=2) in Ry '
        do n=1,nr
         if(coulomb(3))then
          if(n/=1)then
           write(21,*)conv*(pp(3,n)-znuc*r(n)**(-1.d0))
          else
           write(21,*)conv*(pp(3,2)-znuc*r(2)**(-1.d0))
          endif
         else
          write(21,*)conv*pp(3,n)
         endif
        enddo
       endif
       if(ncomp>3)then
        write(21,*)'r*potential (L=3) in Ry '
        do n=1,nr
         if(coulomb(4))then
          if(n/=1)then
           write(21,*)conv*(pp(4,n)-znuc*r(n)**(-1.d0))
          else
           write(21,*)conv*(pp(4,2)-znuc*r(2)**(-1.d0))
          endif
         else
          write(21,*)conv*pp(4,n)
         endif  
        enddo
       endif
       if(ncomp>4)then
        write(21,*)'r*potential (L=4) in Ry '
        do n=1,nr
         if(coulomb(5))then
          if(n/=1)then
           write(21,*)conv*(pp(5,n)-znuc*r(n)**(-1.d0))
          else
           write(21,*)conv*(pp(5,2)-znuc*r(2)**(-1.d0))
          endif
         else
          write(21,*)conv*pp(5,n)
         endif  
        enddo
       endif       
      endif

      write(6,*)'Transformed data written to PP_out'
      write(6,*)
      write(6,*)'Note that the CASINO x_pp.data format assumes that the'
      write(6,*)'pseudopotential components are multiplied by r. If you intend'
      write(6,*)'to use the data just written to PP_out in CASINO, you should' 
      write(6,*)'now use ptm again to do the r multiplication.'
      write(6,*)
      stop
      END SUBROUTINE pp2

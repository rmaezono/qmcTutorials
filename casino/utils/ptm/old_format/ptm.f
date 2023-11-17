!*************************************************************************
!                                                                        *
! PTM                                                                    *
! Perform various transformations on QMC pseudopotentials.               *
! Mike Towler Sep 1999                                                   *
!                                                                        *
! PTM will read in the QMC x_pp.data format, or Gaussian expansions      *
! of such functions in the format documented in the accompanying         *
! README file.                                                           *
!                                                                        *
! PTM cam perform the following transformations on x_pp.data files:      *
! (1) Multiply or divide data by a power of r                            *
! (2) Calculate differences between different components                 *
! (3) Convert to xmgr format                                             *
!                                                                        *
! and the following transformations on the fit file                      *
! (1) perform the Gaussian expansion and plot results on radial grid     *
! wth options to subtract Z/r and to convert between Hartree and Rydbergs*
! The result can be plotted as x_pp.data or xmgr plot.                   *
!                                                                        *
!                                                                        *
!*************************************************************************
      MODULE transform_defs
      INTEGER ialloc,rpower,rnpower
      CHARACTER*72 filename
      END MODULE transform_defs

      PROGRAM pp_transform
      USE transform_defs
      IMPLICIT NONE
      INTEGER i
      LOGICAL ex
      write(6,*)'PTM'
      write(6,*)'==='
10    write(6,*)'choose input type:'
      write(6,*)
     & '(1) Up to 3 pp components on a radial grid in x_pp.data format'
      write(6,*)'(2) gaussian exponents, coefficients and powers of r'
      read(5,*)i
      if(i.ne.1.and.i.ne.2)goto 10
      write(6,*)'Name of input file?'
      read(5,'(a)')filename
      inquire(file=filename,exist=ex)
      if(.not.ex)then
       write(6,*)'This file does not exist'
       stop
      endif
      if(i.eq.1)then
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
      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: r,ys,yp,yd,ys_temp,
     & yp_temp,yd_temp,ydiff
      character*1 a
      CHARACTER*4 rlabel(4)
      CHARACTER*5 label1,label2,label3
      CHARACTER*72 line1,line2,line3
      open(20,file=trim(filename),status='old')
13    write(6,*)
     & 'What power of r is this data already multiplied by (0,1,2)?'
      read(5,*)rpower
      if(rpower.ne.0.and.rpower.ne.1.and.rpower.ne.2)goto 13
      rlabel(1)='  r*'
      rlabel(2)='r^2*'
      rlabel(3)='r^3*'
      rlabel(4)='r^4*'
      xmgr=.false.
      again=.true.

14    write(6,*)'How many pp components are there (1,2,3)?'
      read(5,*)ncomp
      if(ncomp.ne.1.and.ncomp.ne.2.and.ncomp.ne.3)goto 14

!.....read data file
      read(20,'(a)',err=15,end=15)line1
      read(20,'(a)',err=15,end=15)line2
      read(20,*,err=15,end=15)ndata
      ialloc=0
      allocate(r(ndata),ys(ndata),yp(ndata),yd(ndata),ydiff(ndata),
     & stat=ialloc)
      read(20,'(a)',err=15,end=15)line3
      do i=1,ndata
       read(20,*,err=15,end=15)r(i)
      enddo
      read(20,*,err=15,end=15)
      do i=1,ndata
       read(20,*,err=15,end=15)ys(i)
      enddo
      if(ncomp.gt.1)then
       read(20,*,err=15,end=15)
       do i=1,ndata
        read(20,*,err=15,end=15)yp(i)
       enddo
      endif
      if(ncomp.gt.2)then
       read(20,*,err=15,end=15)
       do i=1,ndata
        read(20,*,err=15)yd(i)
       enddo
      endif
      goto 30
15    write(6,*)'Problem reading file.'
      stop
c.....what do you want me to do with this?

30    write(6,*)'File read in successfully.'
40    write(6,*)'What do you want me to do with it?'
      write(6,*)'(1) Multiply or divide data by a power of r'
      write(6,*)'(2) Calculate differences between different components'
      write(6,*)'(3) Write it out as an xmgr plot'
      read(5,*)i
      label1='Vs   '
      label2='Vp   '
      label3='Vd   '
      select case (i)

      case(1) ! powers of r
       write(6,*)'Power (+ve or -ve) to multiply by:'
       read(5,*)n
       if(n.gt.2.or.n.lt.-2)then
        write(6,*)'Power out of range -2 to +2'
        stop
       endif
       if(rpower+n.lt.0)then
        write(6,*)'This power not compatible with existing rpower'
        write(6,*)'(See program if you don''t understand this.)'
        stop
       endif
       if(n.lt.0)then
        ys(2:ndata)=ys(2:ndata)*r(2:ndata)**n
        if(ncomp.gt.1)yp(2:ndata)=yp(2:ndata)*r(2:ndata)**n
        if(ncomp.gt.2)yd(2:ndata)=yd(2:ndata)*r(2:ndata)**n
        ys(1)=ys(2)
        if(ncomp.gt.1)yp(1)=yp(2)
        if(ncomp.gt.2)yd(1)=yd(2)
       else
        ys=ys*r**n
        if(ncomp.gt.1)yp=yp*r**n
        if(ncomp.gt.2)yd=yd*r**n
       endif
45     write(6,*)'Write output in:'
       write(6,*)'(1) x_pp.data format for QMC'
       write(6,*)'(2) xmgr format for plotting'
       read(5,*)i
       if(i.ne.1.and.i.ne.2)goto 45
       if(i.eq.2)xmgr=.true.

      case(2) ! difference
       if(ncomp.lt.2)then
        write(6,*)
     &  'I need at least 2 components to calculate differences.'
        stop
       endif
47     write(6,*)'I will calculate Vx-Vy'
50     write(6,*)'x corresponds to (1,2,3)?'
       read(5,*)nx
       if(nx.ne.1.and.nx.ne.2.and.nx.ne.3)goto 50
60     write(6,*)'y corresponds to (1,2,3)?'
       read(5,*)ny
       if(ny.ne.1.and.ny.ne.2.and.ny.ne.3)goto 60
       if(nx.eq.ny)then
        write(6,*)'nx cannot equal ny'
        goto 50
       endif
       if(nx.eq.3.or.ny.eq.3)then
        if(ncomp.eq.2)
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
        if(a.eq.'y')then
         again=.false.
         goto 47
        endif
       endif
70     write(6,*)'Write output in:'
       write(6,*)'(1) x_pp.data format for QMC'
       write(6,*)'(2) xmgr format for plotting'
       read(5,*)i
       if(i.ne.1.and.i.ne.2)goto 70
       if(i.eq.2)xmgr=.true.
      case(3) ! xmgr plot
       xmgr=.true.
       goto 100
      case default
       write(6,*)
       goto 40
      end select


c.....print out the transformed data
100    open(21,file='PP_out',status='unknown')
       if(xmgr)then ! write out in xmgr format
       do i=1,ndata
        write(21,*)r(i),ys(i)
       enddo
       if(ncomp.gt.1)then
        write(21,*)'&'
        do i=1,ndata
         write(21,*)r(i),yp(i)
        enddo
       endif
       if(ncomp.gt.2)then
        write(21,*)'&'
        do i=1,ndata
         write(21,*)r(i),yd(i)
        enddo
       endif
      else ! write out in QMC x_pp.data format
       rnpower=rpower+n
       write(21,*)line1
       write(21,*)line2
       write(21,*)ndata
       write(21,*)line3
       do i=1,ndata
        write(21,*)r(i)
       enddo
       if(rnpower.gt.0)then
        write(21,*)rlabel(rnpower),label1
       else
        write(21,*)label1
       endif
       do i=1,ndata
        write(21,*)ys(i)
       enddo
       if(ncomp.gt.1)then
        if(rnpower.gt.0)then
         write(21,*)rlabel(rnpower),label2
        else
         write(21,*)label2
        endif
        do i=1,ndata
         write(21,*)yp(i)
        enddo
       endif
       if(ncomp.gt.2)then
        if(rnpower.gt.0)then
         write(21,*)rlabel(rnpower),label3
        else
         write(21,*)label3
        endif
        do i=1,ndata
         write(21,*)yd(i)
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
      PARAMETER(maxcomp=3,ngridpoints=3000,maxg=100)
      INTEGER i,j,n,ncomp,np(maxg),nk(maxcomp,maxg),nr,coeffs
      DOUBLE PRECISION e(maxcomp,maxg),c(maxcomp,maxg),znuc,
     & rmax,aa,bb,a,b,r(ngridpoints),pp(maxcomp,ngridpoints),
     & atnum,conv,temp1(maxcomp,maxg),temp2(maxcomp,maxg)
      CHARACTER*1 yes
      LOGICAL xmgr,coulomb(3)
      xmgr=.false.
      coulomb=.false.

! read in the pseudo with following format:

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
! atomic number, effective nuclear charge
      read(20,*)atnum,znuc
! no of components
      read(20,*)ncomp
! no of gaussians in each component
      read(20,*)(np(i),i=1,ncomp)
! read the Gaussians, exponents and powers of r
      do i=1,ncomp
       read(20,*) ! label
       do j=1,np(i)
! exponents, coefficients, power of r
         read (20,*)temp1(i,j),temp2(i,j),nk(i,j)
       enddo
      enddo

      coeffs=0
      if((any(temp1(:,:).lt.0.d0)).and.(any(temp2(:,:).lt.0.d0)))then
       write(6,*)'Negative exponents and coefficients?'
       write(6,*)'Don''t think so.'
       stop
      endif
      if(any(temp1(:,:).lt.0.d0))coeffs=2 ! Gaussian format
      if(any(temp2(:,:).lt.0.d0))coeffs=1 ! CRYSTAL format
      if(coeffs.eq.0)then
1      write(6,*)'Data is in :'
       write(6,*)
     & '(1) CRYSTAL95 format (alpha, coefficient, power of r)?'
       write(6,*)'(2) GAUSSIAN format (coefficient, alpha, power of r)?'
       read(5,*)coeffs
       if(coeffs.ne.1.and.coeffs.ne.2)goto 1
      endif

      if(coeffs.eq.1)then ! CRYSTAL
       e=temp1
       c=temp2
      else ! GAUSSIAN
       e=temp2
       c=temp1
      endif

! construct the grid
      rmax=50.d0
      aa=8.d0
      bb=80.d0
      a=bb*exp(-aa*log(10.d0))/znuc
      b=1/bb
      do i=1,ngridpoints
       r(i)=a*(exp(b*(i-1))-1)
       nr=i
       if(r(i).gt.rmax)goto 200
      enddo

! construct the pseudopotential components at each point on the grid
200   do n=2,nr
       do i=1,ncomp
        pp(i,n)=0.d0
        do j=1,np(i)
         pp(i,n)=pp(i,n)+
     & r(n)**dfloat(nk(i,j))*c(i,j)*dexp(-(e(i,j)*r(n)**2.d0))
        enddo
       enddo
      enddo
      do i=1,ncomp
       pp(i,1)=pp(i,2)
      enddo

      do i=1,ncomp
       write(6,*)
     &  'Do you want to subtract Z/r off the result of the expansion'
       write(6,*)'for component ',i,'? (you do if it''s local)'
       read(5,'(a)')yes
       if(yes.eq.'y')coulomb(i)=.true.
      enddo

300   write(6,*)'Write output in:'
      write(6,*)'(1) x_pp.data format for QMC'
      write(6,*)'(2) xmgr format for plotting'
      read(5,*)i
      if(i.ne.1.and.i.ne.2)goto 300
      if(i.eq.2)xmgr=.true.

400   write(6,*)'Want output in:'
      write(6,*)'(1) Hartree (i.e. leave it alone)'
      write(6,*)'(2) Rydberg (i.e. multiply current data by 2)'
      read(5,*)i
      if(i.ne.1.and.i.ne.2)goto 400
      conv=1.d0
      if(i.eq.2)conv=2.d0

      open(21,file='PP_out',status='unknown')
      if(xmgr)then

       do i=1,nr
        if(coulomb(1))then
         if(i.ne.1)then
          write(21,*)r(i),conv*(pp(1,i)-znuc*r(i)**(-1.d0))
         else
          write(21,*)r(1),conv*(pp(1,2)-znuc*r(2)**(-1.d0))
         endif
        else
         write(21,*)r(i),conv*pp(1,i)
        endif
       enddo
       if(ncomp.gt.1)then
        write(21,*)'&'
        do i=1,nr
         if(coulomb(2))then
          if(i.ne.1)then
           write(21,*)r(i),conv*(pp(2,i)-znuc*r(i)**(-1.d0))
          else
           write(21,*)r(1),conv*(pp(2,2)-znuc*r(2)**(-1.d0))
          endif
         else
          write(21,*)r(i),conv*pp(2,i)
         endif
        enddo
       endif
       if(ncomp.gt.2)then
        write(21,*)'&'
        do i=1,nr
         if(coulomb(3))then
          if(i.ne.1)then
           write(21,*)r(i),conv*(pp(3,i)-znuc*r(i)**(-1.d0))
          else
           write(21,*)r(1),conv*(pp(3,2)-znuc*r(2)**(-1.d0))
          endif
         else
          write(21,*)r(i),conv*pp(3,i)
         endif
        enddo
       endif

      else
       write(21,*)'Pseudopotential in real space'
       write(21,*)'Number of grid points'
       write(21,*)nr
       write(21,*)'R(i) in atomic units'
       do n=1,nr
        write(21,*)r(n)
       enddo
       write(21,*)'r*s pseudopotential in Ry'
       write(21,*)'0.000 0.000'
       do n=2,nr
        if(coulomb(1))then
         if(n.ne.1)then
          write(21,*)conv*(pp(1,n)-znuc*r(n)**(-1.d0))
         else
          write(21,*)conv*(pp(1,2)-znuc*r(2)**(-1.d0))
         endif
        else
         write(21,*)conv*pp(1,n)
        endif
       enddo
       if(ncomp.gt.1)then
        write(21,*)'r*p pseudopotential in Ry'
        do n=1,nr
         if(coulomb(2))then
          if(n.ne.1)then
           write(21,*)conv*(pp(2,n)-znuc*r(n)**(-1.d0))
          else
           write(21,*)conv*(pp(2,2)-znuc*r(2)**(-1.d0))
          endif
         else
          write(21,*)conv*pp(2,n)
         endif
        enddo
       endif
       if(ncomp.gt.2)then
        write(21,*)'r*d pseudopotential in Ry'
        do n=1,nr
         if(coulomb(3))then
          if(n.ne.1)then
           write(21,*)conv*(pp(3,n)-znuc*r(n)**(-1.d0))
          else
           write(21,*)conv*(pp(3,2)-znuc*r(2)**(-1.d0))
          endif
         else
          write(21,*)conv*pp(3,n)
         endif
        enddo
       endif
       if(ncomp.gt.2)then
        write(21,*)'THE CUTOFF RADIUS FOR THE LOCAL PART OF PP'
     * ,'(10-6 tolerance):'
        write(21,*)
        write(21,*)'THE CUTOFF RADIUS FOR THE NON-LOCAL PART OF PP',
     * ' (10-6 tolerance):'
        write(21,*)
        write(21,*)'THE VALUE OF NSAMP : *NSAMP set in input file,',
     *  ' number put here is ignored'
        write(21,*)'1'
        write(21,*)'THE VALUE OF LEXACT :'
        write(21,*)
        write(21,*)'THE VALUE OF NLANG,L_LOCAL :  *Andrew ',
     *  'added L_Local to make it work.'
        write(21,*)
        write(21,*)'ANG MOM OF NON_LOCAL PP :'
        write(21,*)
        write(21,*)'Atomic number and pseudo-charge'
        write(21,*)atnum,znuc
       endif
      endif

      write(6,*)'Transformed data written to PP_out'
      stop
      END SUBROUTINE pp2

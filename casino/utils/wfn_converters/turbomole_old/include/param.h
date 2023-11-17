! CVS ------------------------------------------------------------------
! $Id: param.h,v 2.1 2000/03/28 09:20:14 babs Exp $
! CVS ------------------------------------------------------------------
! ----------------------------------------------------------------------
!
!     parameter statements for turbomodules with DYNAMIC CORE ALLOCATION
!
! ----------------------------------------------------------------------
      integer ndi3,ndi4,ndi9,ndi10,ndi11,ndi13,ndi14,ndi15,ndi17,ndi18
      integer ndi20,ndi21,ndi22,ndi23,ndi24,ndi25,ndi27,ndi30,ndi31
      integer ndirec,ndilmx,ndirt,ndiddd,ndifff,ndiggg,ndihhh,ndiiii
      integer ndibas,ndistr,lmcomm,ndi23r,ndi24r,mxincor,mxsub,mxbond

      integer npmx,npqmx,mxgout
      integer mxtsk,mxrita,mxprc

      parameter (ndi3  =   4000)
      parameter (ndi4  =   8000)
      parameter (ndi9  =   4000)
      parameter (ndi10 =    500)
      parameter (ndi11 =    150)
      parameter (ndi13 =    900)
      parameter (ndi14 =    120)
      parameter (ndi15 =   2000)
      parameter (ndi17 =     16)
      parameter (ndi18 =     80)
      parameter (ndi20 =     20)
      parameter (ndi21 =     40)
      parameter (ndi22 = ndi10*ndi10)
      parameter (ndi23 = 3*ndi10)
      parameter (ndi24 =    242)
      parameter (ndi25 = ndi14*ndi24)
      parameter (ndi27 =    400)
      parameter (ndi30 = 30*ndi10)
      parameter (ndi31 =     20)
      parameter (ndirec= 500000)
!    ------------------------------------------------------------------
!     the following dimensions are are related with maximum l quantum
!     number ndilmx
      parameter (ndilmx=      7)
      parameter (ndirt = 2*ndilmx-1)
      parameter (ndiddd= ndilmx/3-ndilmx/6-ndilmx/9)
      parameter (ndifff= ndilmx/4-ndilmx/8)
      parameter (ndiggg= ndilmx/5-ndilmx/10)
      parameter (ndihhh= ndilmx/6)
      parameter (ndiiii= ndilmx/7)
!    ------------------------------------------------------------------
!     ndiddd, ndifff, ndiggg, ...  indicate if d, f, g, ... functions
!     are allowed (their values are 1 or 0)
!    ------------------------------------------------------------------
      parameter (npmx  = (ndilmx * (ndilmx + 1) / 2))
!     npmx : number of reducible functions for a given l quantum number
      parameter (npqmx = npmx * npmx)
!     npqmx: number of reduc. funct. for p,q-pair
      parameter (mxgout= npqmx * npqmx)
!     mxgout: number of reduc. funct. for p,q,r,s-pair
!-----------------------------------------------------------------------
!     mxtsk: maximal number of parallel tasks in dft-calculations
      parameter (mxtsk=50000)
!     mxrita: maximal number of parallel tasks in RI-Coulomb Part
      parameter (mxrita=2000)
!     mxprc: maximal number of processors
!     needed in parallel RI-Coulomb Part
      parameter (mxprc=512)
!     array with ecp-formulas and some triangular matrices
!     used in new ECP-routines (W. Thiel and J. Breidung Jan 1997)
      parameter (lmcomm=345853)
!-----------------------------------------------------------------------
!     ndi3      max. number of sao's in one irrep (and one column)
!     ndi4      total number of basis functions
!     ndi9      maximum number of shells
!     ndi10     max.number of atoms
!     ndi11     max.number of shell types
!     ndi13     max.number of types of primitive shells
!     ndi14     max.number of symmetry operations in symmetry group
!     ndi15     max. number of occupied mos
!     ndi17     max.number of real representations in symmetry group
!     ndi18     max. number of ao's in each set of symmetryequivalent
!               shells
!     ndi20     max.number of linearly combined density matrices
!     ndi21     max.number of primitives in a given shell
!     ndi22     max.number of atomic pairs, which suffice certain
!               distance criteria
!     ndi23     max.number of "internal" coordinates
!               (= 3*number of atoms to enable optimization in
!                cartesian space !)
!     ndi24     max.number of linearly combined internal coordinates
!               contributing to a given internal coordinate
!               or maximum number of geometries to be used for GDIIS or
!               SCHLEGEL update in geometry optimization tasks
!     ndi25     same as ndi24, but also including symmetry equivalent
!               internal coordinates
!     ndi27     maximum number of parameters for ecps
!     ndi30     maximum number of atomic neighbors in bond analysis
!     ndi31     maximum number of bond partners of a given atom
!     ndilmx    maximum l quantum number+1
!     ndirt     maximum number of roots for numerical evalution
!               of two electron integrals using RYS polynoms
!    ------------------------------------------------------------------

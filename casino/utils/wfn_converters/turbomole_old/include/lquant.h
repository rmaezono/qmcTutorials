! CVS ------------------------------------------------------------------
! $Id: lquant.h,v 2.1 2000/03/17 17:02:05 marcok Exp $
! CVS ------------------------------------------------------------------
!-----------------------------------------------------------------------
!   SPECIAL parameters needed in gradient calculations
!   they are adjusted to enable the program to evaluate
!   -  cartesian gradients for gaussians up to ndilmx
!   -  exponent  gradients for gaussians up to ndilmx-1
!-----------------------------------------------------------------------
!
!   ndilmx                   max. angular momentum number+1
!   nfired                   number of irreducible functions in {ndilmx}
!   nft                      number of reducible functions for a given l
!                               quantum number
!   nftij                    number of reduc. funct. for i,j-pair
!   nftcrt                   number of reducible functions in
!                               {ndilmx-1} + {ndilmx+1}
!                               (needed for cartesian gradients)
!   nftexp                   number of reducible functions in {(ndilmx-1)+2}
!                               (needed for exponent gradients)
!   nftmax                   total number of reducible gaussians up to
!                               {ndilmx+1}
!                               (l+1)*(l+2)*(l+3)/6 { = sum[j(j+1)/2] ,
!                               j=1,l+1 }
!   nfmax                    total number of reducible gaussians up to
!                               {ndilmx} = sum[l(l+1)/2] , l=1,ndilmx
!   ldredn                   nft**n
!-----------------------------------------------------------------------
      integer nfired, nft, nftm1, nftp1, nftij, nftcrt, nftexp
      integer nftmax, nfmax, ldred3, ldred4

      parameter (nfired = 2*ndilmx-1)
      parameter (nft    = ndilmx*(ndilmx+1)/2)
      parameter (nftm1  = (ndilmx-1)*ndilmx/2)
      parameter (nftp1  = (ndilmx+1)*(ndilmx+2)/2)
      parameter (nftij  = nft*nft)
      parameter (nftcrt = nftm1+nftp1)
      parameter (nftexp = nftp1)
      parameter (nftmax = (ndilmx+1)*(ndilmx+2)*(ndilmx+3)/6)
      parameter (nfmax  = ndilmx*(ndilmx+1)*(ndilmx+2)/6)
      parameter (ldred3 = nft*nft*nft)
      parameter (ldred4 = ldred3*nft)
!----------------------------------------------------------------------------

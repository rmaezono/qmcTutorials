 MODULE sphfac_module
  USE numbers
  USE lmaxxx
  REAL(float),DIMENSION(lmax_dft200) :: fsph
  INTEGER,DIMENSION(0:lmax_dft7) :: nu3p
  INTEGER,DIMENSION(lmax_dft49+1) :: ietap
  INTEGER,DIMENSION(lmax_dft200) :: net1
  INTEGER,DIMENSION(lmax_dft200) :: net2
  INTEGER,DIMENSION(lmax_dft200) :: net3
 END MODULE sphfac_module

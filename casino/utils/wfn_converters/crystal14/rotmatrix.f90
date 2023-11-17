 MODULE rotmatrix
  USE numbers
  USE lmaxxx
! Allocatable vectors for rotation matrices (dimension limrot1)
  INTEGER limrot1
  REAL(float),ALLOCATABLE :: tto(:),zo34(:)
  INTEGER,ALLOCATABLE :: mmo(:),mmom(:)
! Allocatable vectors for rotation matrices product (dimension limrot2)
  INTEGER limrot2
  INTEGER,ALLOCATABLE :: no34(:),nom34(:),ko34(:)
! Other related vectors & matrices
  INTEGER minz(48*(lmax_dft7+1)+1)
  INTEGER novf34(0:48,0:lmax_dft7,0:lmax_dft7)
  INTEGER lovf34(0:48,0:lmax_dft7)
  INTEGER,ALLOCATABLE :: mom12(:),no12(:),mo34(:)
 END MODULE rotmatrix

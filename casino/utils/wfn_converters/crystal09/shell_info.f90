 MODULE shell_info
! Variables related to shell properties
  USE numbers
  USE lmaxxx
! Name, L, L+1, No.AO, 2*No.AO
  CHARACTER,PARAMETER :: shtxt(0:lmax_dft7)*2 = &
   &(/'S ','SP','P ','D ','F ','G ','H ','I '/)
  INTEGER,PARAMETER :: shmxl(0:lmax_dft7)=&
   &(/  0 ,  1 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 /)
  INTEGER,PARAMETER :: shnao(0:lmax_dft7)=&
   &(/  1 ,  4 ,  3 ,  5 ,  7 ,  9 , 11 , 13 /)
  INTEGER,PARAMETER :: nshchg(0:lmax_dft7)=&
   &(/  2 ,  8 ,  6 , 10 , 14 , 18 , 22 , 26 /)
  INTEGER,PARAMETER :: shnaosq(0:lmax_dft7)=&
   &(/  1 , 17,  26 , 51 ,100 ,181 ,302 ,471 /)
! Name of AO's
  CHARACTER,PARAMETER :: aotxt(0:lmax_dft6,lmax_dft13)*3=reshape(&
   &(/'  S','---','---','---','---','---','---','---','---',   &
   &'---','---','---','---',                           &
   &' PX',' PY',' PZ','---','---','---','---','---','---',   &
   &'---','---','---','---',                           &
   &'D00','D+1','D-1','D+2','D-2','---','---','---','---',   &
   &'---','---','---','---',                           &
   &'F00','F+1','F-1','F+2','F-2','F+3','F-3','---','---',   &
   &'---','---','---','---',                           &
   &'G00','G+1','G-1','G+2','G-2','G+3','G-3','G+4','G-4',   &
   &'---','---','---','---',                           &
   &'H00','H+1','H-1','H+2','H-2','H+3','H-3','H+4','H-4',   &
   &'H+5','H-5','---','---',                           &
   &'I00','I+1','I-1','I+2','I-2','I+3','I-3','I+4','I-4',   &
   &'I+5','I-5','I+6','I-6' /),                        &
   &(/lmax_dft7,lmax_dft13/) )
! No. of different elements in a square symmetric matrix
  INTEGER,PARAMETER :: iky(0:lmax_dft13)=&
   &(/  0,  1,  3,  6, 10, 15, 21, 28, 36, 45, 55, 66, 78, 91/)
 END MODULE shell_info

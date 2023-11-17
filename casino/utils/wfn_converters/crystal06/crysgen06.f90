PROGRAM crystal2casino
 USE casino_interface
 USE xyvdim_module

 write(6,*)
 write(6,'(1x,a)')'CRYSTAL06-CASINO INTERFACE FOR QUANTUM MONTE CARLO &
  &CALCULATIONS'
 write(6,'(1x,a)')'---------------------------------------------------&
  &------------'
 write(6,*)

 write(6,*)'Reading data...'
 write(6,*)

 call cryread
 call kredin
 call vrslat
 call condft
 call gsym11(xyv)
 call qmc_main

END PROGRAM crystal2casino

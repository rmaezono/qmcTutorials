 MODULE lmaxxx
  IMPLICIT INTEGER (i-n)
  PARAMETER (max_gvec=48)
  PARAMETER (max_buff=1000)
! l max values for mono
  PARAMETER (lmax_mul6=6)
  PARAMETER (lmax_mul7=lmax_mul6+1)
  PARAMETER (lmax_mul8=lmax_mul6+2)
  PARAMETER (lmax_mul12=lmax_mul6+lmax_mul6)
  PARAMETER (lmax_mul13=lmax_mul6+lmax_mul7)
  PARAMETER (lmax_mul49=lmax_mul7*lmax_mul7)
  PARAMETER (lmax_mul84=(lmax_mul8+1)*lmax_mul8*lmax_mul7/6)
  PARAMETER (lmax_mul120=lmax_mul8*(lmax_mul8+1)*(lmax_mul8+2)/6)
  PARAMETER (lmax_mul196=lmax_mul49*4)
  PARAMETER (lmax_mul36=lmax_mul8*(lmax_mul8+1)/2)
  PARAMETER (lmax_mul28=lmax_mul8*lmax_mul7/2)
  PARAMETER (lmax_mul200=lmax_mul7*lmax_mul6*(lmax_mul6-1)) ! upper bound if
                                                            ! lmax_mul6>=4
! l max values for DFT fitting basis
  PARAMETER (lmax_dft6=6)
  PARAMETER (lmax_dft7=lmax_dft6+1)
  PARAMETER (lmax_dft28=lmax_dft7*(lmax_dft7+1)/2)
  PARAMETER (lmax_dft13=lmax_dft6+lmax_dft7)
  PARAMETER (lmax_dft49=lmax_dft7*lmax_dft7)
  PARAMETER (lmax_dft91=lmax_dft13*(lmax_dft13+1)/2)
  PARAMETER (lmax_dft200=lmax_dft7*lmax_dft6*(lmax_dft6-1)) ! upper bound if
                                                            ! lmax_dft6>=4
! l max values for basis set
  PARAMETER (lmaxx2=3)
  PARAMETER (lmaxx3=lmaxx2+1)
  PARAMETER (lmaxxfour=lmaxx3+1)
  PARAMETER (lmaxx4=lmaxx2+lmaxx2)
  PARAMETER (lmaxx8=lmaxx4+lmaxx4)
  PARAMETER (lmaxx9=lmaxx8+1)
  PARAMETER (lmaxx5=lmaxx4+1)
  PARAMETER (lmaxx6=lmaxx4+2)
  PARAMETER (lmaxx33=lmaxx3*lmaxx3)
  PARAMETER (lmaxx26=lmaxx3*lmaxx33-1)
  PARAMETER (lmaxxsix=lmaxx3*(lmaxx3+1)/2)
  PARAMETER (lmaxx200=lmaxx3*lmaxx2*lmaxx2) ! upper bound if lmaxx2>=2
  PARAMETER (lmaxx25=lmaxx5*lmaxx5)
  PARAMETER (lmaxx35=(lmaxx5+2)*lmaxx6*lmaxx5/6)
  PARAMETER (lmaxx56=(lmaxx5+3)*(lmaxx5+2)*lmaxx6/6)
  PARAMETER (lmaxx84=(lmaxx5+4)*(lmaxx5+3)*(lmaxx5+2)/6)
  PARAMETER (lmaxx165=(lmaxx9+2)*(lmaxx9+1)*lmaxx9/6)
  PARAMETER (lmaxx220=(lmaxx9+3)*(lmaxx9+2)*(lmaxx9+1)/6)
  PARAMETER (lmaxx3525=lmaxx35*lmaxx25)
  PARAMETER (lmaxx5625=lmaxx56*lmaxx25)
  PARAMETER (lmaxx5635=lmaxx56*lmaxx35)
  PARAMETER (lmaxx2525=lmaxx25*lmaxx25)
  PARAMETER (lmaxx2525482=lmaxx2525*max_gvec*2)
  PARAMETER (lmaxx2525486=lmaxx2525*max_gvec*6)
  PARAMETER (lmaxx25254812=lmaxx2525*max_gvec*12)
  PARAMETER (lmaxx3515=lmaxx5*(lmaxx2+3)*(lmaxx2+4)*(lmaxx2+5)/2)
  PARAMETER (lmaxx3535=lmaxx35*lmaxx35)
  PARAMETER (lmaxxee=lmaxx6+lmax_mul6)
  PARAMETER (lmaxx364=lmaxxee*(lmaxxee+1)*(lmaxxee+2)/6)
  PARAMETER (lmaxx286=(lmaxxee-1)*lmaxxee*(lmaxxee+1)/6)
  PARAMETER (lmaxx78=lmaxxee*(lmaxxee+1)/2)
  PARAMETER (lmaxx66=lmaxxee*(lmaxxee-1)/2)
  PARAMETER (lmaxx11=lmaxx5+lmax_mul6) ! USE if lmaxx4 <= lmaxmul6
!  PARAMETER (lmaxx11=lmaxx5+lmaxx4) ! USE if lmaxx4 > lmaxmul6
  PARAMETER (lmaxx12=lmaxx11+1)
  PARAMETER (lmaxxfive=(lmaxx12-2)/2)
  PARAMETER (lmaxx14=lmaxx11+3)
  PARAMETER (lmaxx14014=(lmaxx14+1)*1001)
  PARAMETER (lmaxx560=lmaxx12*(lmaxx12+1)*lmaxx14/6)
  PARAMETER (lmaxx2549=lmaxx25*lmax_mul49)
  PARAMETER (lmaxx254948=lmaxx2549*max_gvec)
! l max values for libxg (derivatives of charge density)
  PARAMETER (max_grad4=4)
  PARAMETER (max_grad8=max_grad4+lmaxx4)
  PARAMETER (max_grad35=(max_grad4+1)*(max_grad4+2)*(max_grad4+3)/6)
  PARAMETER (max_grad45=(max_grad8+1)*(max_grad8+2)/2)
! l max values for bipolar coulomb
  PARAMETER (lmax_coul4=4)
  PARAMETER (lmax_coul5=lmax_coul4+1)
  PARAMETER (lmax_coul9=lmax_coul4+lmax_coul5)
  PARAMETER (lmax_coul25=lmax_coul5*lmax_coul5)
! l max values for bipolar exchange
  PARAMETER (lmax_exch2=2)
  PARAMETER (lmax_exch3=lmax_exch2+1)
  PARAMETER (lmax_exch5=lmax_exch2+lmax_exch3)
  PARAMETER (lmax_exch9=lmax_exch3*lmax_exch3)
! Max number of primitives per contraction
  PARAMETER (max_prim10=10)
  PARAMETER (max_prim100=max_prim10*max_prim10)
! l max values for pseudopotentials
  PARAMETER (lmax_psproj=4,lmax_psnpot=2)
  PARAMETER (lmax_pse6=lmaxx2+lmax_psproj) ! USE if lmaxx2<=lmax_psproj
!  PARAMETER (lmax_pse6=lmaxx4)  ! USE if lmaxx2>lmax_psproj
  PARAMETER (lmax_pse7=lmax_pse6+1)
  PARAMETER (lmax_pse8=lmax_pse6+2)
  PARAMETER (lmax_pse13=lmax_pse6+lmax_pse7)
  PARAMETER (lmax_pse28=lmax_pse8*lmax_pse7/2)
  PARAMETER (lmax_pse200=lmax_pse7*lmax_pse6*(lmax_pse6-1)) ! upper bound if
                                                            ! lmax_pse6.ge.4
  PARAMETER (lmax_pse12=lmax_pse6+lmax_pse6)
  PARAMETER (lmax_pse49=lmax_pse7*lmax_pse7)
  PARAMETER (lmax_ps2tot=lmaxx2+lmax_psproj)
  PARAMETER (lmaxps=lmax_ps2tot) ! USE if lmax_psproj>=lmaxx2
! PARAMETER (lmaxps=lmaxx4) ! USE if lmax_psproj<lmaxx2
  PARAMETER (lmaxpsdub=lmaxps+lmaxps)
  PARAMETER (lmax_ps2tot2=lmax_ps2tot+lmax_ps2tot)
  PARAMETER (lmax_ps4tot=lmaxx4+lmax_psnpot)
  PARAMETER (lmax_ps4vrs=lmaxx4+lmax_psnpot)  ! USE if lmax_psnpot>=e.2
!  PARAMETER (lmax_ps4vrs=lmaxx6)  ! USE if lmax_psnpot<2
  PARAMETER (lmax_ps7tot=lmaxx8+lmax_psnpot-1)
  PARAMETER (lmax_psqlng=(lmax_ps4tot+3)*(lmax_ps4tot+4)/2)
  PARAMETER (lmax_psqmmdg=(lmaxx4+3)*(lmaxx4+4)/2)
 END MODULE lmaxxx

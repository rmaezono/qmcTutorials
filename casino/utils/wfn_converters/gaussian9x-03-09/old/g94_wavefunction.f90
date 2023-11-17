MODULE g94_wavefunction
  USE paramfile

  CHARACTER*71, SAVE  :: title_txt ! Holds the G94 job title
  CHARACTER*71, SAVE  :: job_txt   ! Job details - method used etc.
  CHARACTER*11, SAVE  :: code_used ! Holds the version (94 or 98) of
                                   ! Gaussian used to do the calc.
  CHARACTER*30, SAVE  :: g94_file  ! The root of the Gaussian job name.
                                   ! The formatted checkpoint file is
                                   ! assumed to have .Fchk appended to
                                   ! this whilst the G94 output has the
                                   ! .out suffix

  INTEGER, SAVE :: ispin_lim=2 ! Limit for loop over spins - is reset
                               ! to 1 if calculation is not spin
                               ! polarised (i.e. we have only alpha
                               ! electrons)
  INTEGER, SAVE :: Natom   ! No. of atoms in the system
  INTEGER, SAVE, ALLOCATABLE :: Natomic_no(:) ! Atomic no. of each atom

  INTEGER, SAVE :: Nelec   ! No. of electrons in the system

  ! No. of coefficients defining a single MO.
  INTEGER, SAVE :: Ncoeffs

  ! The total number of alpha (or beta) MOs.
  INTEGER, SAVE :: Nmo

  ! Max. AM quantum no. of any shell in the basis
  INTEGER, SAVE :: max_AM

  ! Max. degree of contraction of any shell
  INTEGER, SAVE :: Max_contract

  ! No. of contracted shells (shell = group of f'ns with same
  ! angular dependence)
  INTEGER, SAVE :: Nshells

  ! The number of SP shells in the basis set
  INTEGER, SAVE :: Num_sp=0

  ! No. of primitive gaussian type functions (e.g. p-type gaussian f'n etc.)
  INTEGER, SAVE :: Nprimgtf

  ! No. of gaussian primitives in each shell contraction
  INTEGER, SAVE, ALLOCATABLE :: Nprim(:)

  INTEGER, SAVE :: Nspin(2) ! No. of alpha (1) and beta (2) electrons

  REAL (KIND=DP), SAVE :: Eionion ! Nuclear repulsion energy (au/atom)

! atom_posns(i,j) holds x_{i} for atom j
  REAL (KIND=DP), ALLOCATABLE, SAVE :: atom_posn(:,:)

! atomic_chrg holds the actual valence charge on each atom (to allow for
! charged systems and pseudopotentials etc.)
  REAL (KIND=DP), ALLOCATABLE, SAVE :: atomic_chrg(:)

! Shell to atom map - ish_map(i) holds which atom shell i is centred on
  INTEGER, ALLOCATABLE, SAVE        :: ish_map(:)

! Total no. of orbitals to allow room for when storing eigenvector
! coefficients.  It must be at least Nmo plus the number of electrons
! in the system (for CIS resumming).  Is set in read_fchk.
  INTEGER, SAVE                     :: Nmo_total

! shll_posns(i,j) holds x_{i} for shell j
  REAL (KIND=DP), ALLOCATABLE, SAVE :: shll_posn(:,:)

! con_coeff(primitive index, shell index) - Shell contraction coeff.'s
  REAL (KIND=DP), ALLOCATABLE, SAVE :: con_coeff(:,:)
! Same info. indexed by primitive shell so suitable for output to QMC
! input file
  REAL (KIND=DP), ALLOCATABLE, SAVE :: c_prim(:)

! As for con_coeff but holds the contraction coeff.'s for the p functions
! in a contracted SP shell.  (Those for the s are stored in con_coeff.)
  REAL (KIND=DP), ALLOCATABLE, SAVE :: sp_coeff(:,:)
! Same info. indexed by primitive shell so suitable for output to QMC
! input file
  REAL (KIND=DP), ALLOCATABLE, SAVE :: c2_prim(:)

! evcoeff1(ifun,imo,ispin) - MO expansion coefficients, ifun loops
! over basis functions, imo indexes the imo'th MO of spin ispin(=1 for alpha
! and 2 for beta). Allocated in read_fchk.
  REAL (KIND=DP), ALLOCATABLE, SAVE :: evcoeff1(:,:,:)

! More efficient storage whereby the AM quantum numbers and exponents are
! stored by shell index.  How many primitives have this index must then be
! worked out from the AM of the shell.
  INTEGER, ALLOCATABLE :: Lshell(:)
! Exponent of each primitive shell
  REAL (KIND=DP), ALLOCATABLE, SAVE :: shexpnt(:)

! hf_spectrum(i,ispin) holds the HF eigenvalue for (MO i of spin ispin)
! read from the output file (for the purposes
! of analysing a CIS excited state - only read if we have a CIS
! calculation)
  REAL (KIND=DP), ALLOCATABLE, SAVE :: hf_spectrum(:,:)

! Switch to identify whether there are CIS wavefunctions present in the
! Gaussian output file (CIS=.TRUE.) or not (.FALSE.).
! Runtime value set in read_G9xout.f90.
  LOGICAL, SAVE :: CIS =.true.

! Switch to flag whether we are dealing with a CASSCF/MCSCF wavefunction
! Runtime value set in read_G9xout.f90.
  LOGICAL, SAVE :: CAS=.true.

! Switch for spin-unrestricted (=.true.) or spin-restricted (=.false.)
! calculation
  LOGICAL, SAVE :: SPIN=.true.

END MODULE g94_wavefunction

MODULE cis_data
!------------------------------------------------------------------------!
! Contains the description of either the CIS or CASSCF wave function.    !
!------------------------------------------------------------------------!
 USE paramfile

!------------------------------------------------------------------------!
! Data relating to the CIS wave function - currently assumes a maximum   !
! of 10000 configurations per state.                                     !
!------------------------------------------------------------------------!
 INTEGER,PARAMETER :: Ncon_max=10000

! Which CIS state to output.  If this is set to zero then we output the
! GS, even if there were CIS states in the Gaussian output file.
 INTEGER,SAVE :: icis_out

! No. of excited states
 INTEGER,SAVE :: Nexcite

! No. of alpha/beta (index=1/2) configurations in the expansion of the CIS
! excited state selected for output
 INTEGER,SAVE :: Nconfig(2)

! No. of distinct occupied MOs in GS wave function
 INTEGER,SAVE :: Nocc_mo

! No. of MO eigenvalues for each spin successfully read from output file.
! This will be less than Nmo because of formatting errors in Gaussian's
! output of large eigenvalues.
 INTEGER,SAVE :: Nmo_read(2)

! Orbitals(1,i,ispin) - The index of the ispin (=1 for alpha,2 for beta)
! occ. orbital excited from for configuration i in the CIS state selected
! for output
! Orbitals(2,i,ispin) - The index of the alpha virtual orbital excited
! into for configuration i in the CIS state selected for output
 INTEGER,SAVE :: Orbitals(2,Ncon_max,2)

! Npromote(imo,ispin) - no. of promotions from each occupied (alpha/beta)
! MO that make up the configurations in the expansion of the icis_out
! excited state being resummed. Set in re_sum.
 INTEGER,ALLOCATABLE,SAVE :: Npromote(:,:)

 INTEGER,SAVE :: Ndet(2) ! Counts the no. of resummed determinants
                         ! that we'll need for the alpha (1) and beta (2)
                         ! excitations

! ci_coeff(i,ispin) - the coefficient of  configuration i in the CIS state
 REAL(KIND=dp),SAVE :: ci_coeff(Ncon_max,2)

! The sum of squares of the CIS/TD-DFT coefficients.  Used in calculating
! the percentage break-down of an excited state wave function.
 REAL(KIND=dp),SAVE :: ci_norm

! True if excited state selected for output is a singlet, false otherwise.
! Only referenced for spin-restricted calculations.
 LOGICAL,SAVE :: singlet=.false.

! Flag to identify whether CIS expansion has been resummed (.true.) or not.
 LOGICAL,SAVE :: resummed=.false.

! Whether to analyze the percentage breakdown of the CIS excited state
 LOGICAL,SAVE :: analyze_cis=.true.

!--------------------------------------------------!
! Data relating to the CASSCF/MCSCF wave function. !
!--------------------------------------------------!

! The maximum no. of configurations that G98/4 prints for the CAS
! eigenvector.  Currently this is just 50 but maybe it will increase in the
! future
 INTEGER,PARAMETER :: Ncas_config_max=50

! Which state the CASSCF calculation was converged on - this is the one we
! pull from the output file. Initialise to recognisably dummy value.
 INTEGER,SAVE :: Nroot=999

! Number of frozen, occupied MOs in CASSCF
 INTEGER,SAVE :: Nfrozen

! Multiplicity of the CAS wave function
 INTEGER,SAVE :: Multiplicity=0

! Number of orbitals in the active space
 INTEGER,SAVE :: Nact_orbs

! Number of electrons in the active space
 INTEGER,SAVE :: Nact_elecs

! Number of alpha and beta electrons in the active space
 INTEGER,SAVE :: Nspin_cas(2)

! No. of resummed/effective MOs created when resumming a CAS
! wave function Set in resum_cas and used in qmc_write.
 INTEGER,SAVE :: Nmo_resummed

! No. of MOs available for storing resummed MOs.  If this is less than
! Nmo_resummed then we cannot use a resummed CAS wave function without
! modifying CASINO (because CASINO assumes that Nmo=Ncoeffs).
  INTEGER,SAVE :: Nfree_mo

! Stores the makeup of the reference state from which all CAS
! configurations are generated.  Spin-up MOs are listed first
! followed by spin-down.
 INTEGER,ALLOCATABLE,SAVE :: icas_ref(:)

! Stores the makeup of all possible configurations as listed by Gaussian
! icasdet(i,iconfig) is the orbital occupied by active electron i in
! the iconfig'th config.  i runs over all active electrons beginning with
! those that are spin up.
 INTEGER,ALLOCATABLE,SAVE :: icasdet(:,:)

! Size of CAS basis set - i.e. No. of determinants in the CASSCF procedure
 INTEGER,SAVE :: Ncas_basis

! Stores the labels of the configurations involved in the CAS expansion of
! the state of interest...
 INTEGER,ALLOCATABLE,SAVE :: icas_config(:)

! Lists MOs that do not occur in any of the CAS configurations and are
! hence free to store resummed MOs.  Dimensioned to Nmo in read_fchk.
 INTEGER,ALLOCATABLE,SAVE :: list_of_free_mos(:)

! ...corresponding expansion coefficients
 REAL(KIND=dp),ALLOCATABLE,SAVE :: cas_coeffs(:)

! Flag is T if CAS wave function expanded in Slater determinants and
! F if it is expanded in spin configurations.
 LOGICAL,SAVE :: slater

! True if this determinant has been involved in a resummation
! Initialised and set in resum_cas
 LOGICAL,ALLOCATABLE,SAVE :: resummed_det(:)
! Only true if determinant has been involved in a resummation and
! subsequently deleted.  Initialised and set in resum_cas
 LOGICAL,ALLOCATABLE,SAVE :: deleted_det(:)

! No. of determinants deleted from expansion through resumming
 INTEGER,SAVE :: num_det_deleted

END MODULE cis_data

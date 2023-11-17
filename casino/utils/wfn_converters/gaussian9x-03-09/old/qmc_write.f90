SUBROUTINE qmc_write
!=============================================================
! Output the wavefunction obtained from the .Fchk and .out
! files in a form suitable for CASINO to take as input.
!=============================================================
  USE g94_wavefunction
  USE cis_data
  USE awk_like
  IMPLICIT none
  INTEGER, PARAMETER      :: Num_dft=10

  ! FORMAT statements
  CHARACTER*13, PARAMETER :: reals3="(3(1PE20.13))"
  CHARACTER*13, PARAMETER :: reals4="(4(1PE20.13))"
  CHARACTER*13, PARAMETER :: int8="(8I10)"
  ! dft holds the strings that must occur in the job description
  ! if a DFT job is being performed
  CHARACTER*5             :: dft(Num_dft)
  CHARACTER*130           :: line_input ! Use for copying lines
                                        ! of external files
  CHARACTER*34            :: g94_outfile ! G94 output filename
  CHARACTER*34            :: filenames(2) ! G94 job and output
                                          ! filenames
  CHARACTER*15            :: mthd ! Holds the bit of string that
                                  ! identifies the type of calc.
  ! ishll_type() holds the identity of each shell in the form
  ! used in CASINO:
  ! s=1, sp=2, p=3, d= 4, f= 5 etc.. (harmonic representation)
  !                 d=-4, f=-5       (cartesian representation)
  INTEGER :: ishll_type(Nshells)
  INTEGER :: Num_centres            ! No. of unique shell centres
  INTEGER :: ifirst_shll(Nshells+1) ! Index of the first shell at
                                    ! each unique centre
  INTEGER :: i, ishll, iband, ifun, ifun0, ifun1, nfun, icoeff
  INTEGER :: num_loop, iocc, ispin, iname, idet
  INTEGER :: ilook, iscan,ifunc ! For parsing job/calculation type
  INTEGER :: icheck ! Flag for error checking on IO operations
  INTEGER :: ivirt  ! Indexes the resummed excited state det.'s
                    ! and thus the resummed MO used to create each
  INTEGER :: ndet_out ! No. of det.'s in CIS expansion to output

  REAL (KIND=DP) :: ck(100) ! Holds the eigenvector coeffs in the
                          ! order req'd by QMC (100 at a time to
                          ! save memory).
  LOGICAL :: scan ! Use to control scanning of external files
  LOGICAL :: dftcalc ! .true. if Gaussian calc. was a DFT calc.

  WRITE(6,FMT="(/'Writing QMC wavefunction to ',A)") &
                                            &TRIM(g94_file)//'.qmc'

  OPEN(unit=22,file=TRIM(g94_file)//'.qmc',status='unknown')

  dft=(/"LSDA ","LDA  ","VWN  ","LYP  ","PL   ","P86  ","PW91 ", &
       &"BECKE","B3P  ","SP81 "/)

!=================Title Section==================================
  WRITE(22,FMT=*) TRIM(title_txt)

!=================Basic Information==============================

  WRITE(22,FMT="(/'BASIC_INFO'/'----------')")

  WRITE(22,FMT="('Generated by:'/A11)") code_used

 ! Work out whether we are doing a density functional calc.
  CALL scan_string(job_txt,fields,NF,NF_max," ")
  mthd=TRIM(fields(2))

  call capitalise(mthd,LEN(mthd))

  dftcalc=.false.
  DO ilook=1,Num_dft,1
    iscan=INDEX(STRING=mthd, SUBSTRING=TRIM(dft(ilook)))
    IF(iscan .ne. 0)THEN
      dftcalc=.true.
      ifunc=iscan
      EXIT
    END IF
  END DO

  WRITE(22,FMT="('Method:')")
  WRITE(22,FMT="(A11)") mthd
  WRITE(22,FMT="('DFT Functional:')")
  IF(dftcalc)THEN
    ! mthd will not just hold the functional at the moment - it may
    ! also hold more stuff about the calculation so try and
    ! remove it - assumes such bits are separated by hyphens.
    CALL scan_string(mthd(ifunc:LEN(mthd)),fields,NF,NF_max,"-")
    WRITE(22,*) TRIM(fields(1))
  ELSE
    WRITE(22,FMT="('none')")
  END IF

  WRITE(22,FMT="('Periodicity:'/'0')")

  WRITE(22,FMT="('Spin unrestricted:')")
  IF(SPIN)THEN
    WRITE(22,FMT="('.true.')")
  ELSE
    WRITE(22,FMT="('.false.')")
  END IF

  WRITE(22,FMT="('nuclear-nuclear repulsion energy (au/atom):')")
  WRITE(22,*) Eionion

  WRITE(22,FMT="('Number of electrons per primitive cell:')")
  WRITE(22,*) Nelec

!=================Geometry input=======================================

  WRITE(22,FMT="(/'GEOMETRY'/'--------')")

  WRITE(22,FMT="('Number of atoms:')")
  WRITE(22,*) Natom

  WRITE(22,FMT="('Atomic positions (au):')")
  WRITE(22,FMT=reals3) atom_posn

  WRITE(22,FMT="('Atomic numbers for each atom:')")
  WRITE(22,FMT=int8) Natomic_no

  WRITE(22,FMT="('Valence charges for each atom:')")
  WRITE(22,FMT=reals4) atomic_chrg

!=================Basis set==============================================

  WRITE(22,FMT="(/'BASIS SET'/'---------')")

! Identify how many unique shell centres we have and the index of
! the first shell at each such centre
  CALL shell_centres(Num_centres,ifirst_shll)

  WRITE(22,FMT="('Number of Gaussian centres')")
  WRITE(22,*) Num_centres

  WRITE(22,FMT="('Number of shells per primitive cell')")
  WRITE(22,*) Nshells

! If Ncoeffs /= Nmo then we create some fake MOs when we write out the
! eigenvectors in order to keep Casino happy.
  WRITE(22,FMT="('Number of basis functions (''AO'') per primitive&
                 & cell')")
  WRITE(22,*) Ncoeffs

  WRITE(22,FMT="('Number of Gaussian primitives per primitive cell')")
  WRITE(22,*) Nprimgtf

! Have to shift Gaussian94's value for this since L=0,1,2 for s,p,d in
! that program and QMC expects L=1 for s etc. etc.
  WRITE(22,FMT="('Highest shell angular momentum (s/p/d/f... 1/2/3/4...)')")
  WRITE(22,*) (max_AM+1)

! Convert from Gaussian94 shell types to QMC notation:
  DO i=1,Nshells,1

    shell_type: SELECT CASE( ABS(Lshell(i)) )

    CASE(0)
    ! S shell
      ishll_type(i)=1
    CASE(1)
    ! SP or P shell
      IF(Lshell(i) .eq. -1)THEN
        ishll_type(i)=2 ! SP shell
      ELSE
        ishll_type(i)=3 ! P shell
      END IF
    CASE(2)
    ! D shell - harmonic if -ve, cartesian if +ve but needs to  be
    ! the other way around for CASINO
      ishll_type(i)=-SIGN(4,Lshell(i))
    CASE(3)
    ! F shell
      ishll_type(i)=-SIGN(5,Lshell(i))
    CASE(4)
    ! G shell
      ishll_type(i)=-SIGN(6,Lshell(i))
    CASE(5)
    ! H shell
      ishll_type(i)=-SIGN(7,Lshell(i))
    CASE(6)
    ! I shell
      ishll_type(i)=-SIGN(8,Lshell(i))
    CASE DEFAULT
      WRITE(*,FMT="('Error in qmc_write: shell type is not supported, &
                    &Lshell = ',I2)") Lshell(i)

    END SELECT shell_type

  END DO

  WRITE(22,FMT="('Code for shell types (s/sp/p/d/f... 1/2/3/4/5...)')")
  WRITE(22,FMT=int8) ishll_type

  WRITE(22,FMT="('Number of primitive Gaussians in each shell')")
  WRITE(22,FMT=int8) Nprim

  WRITE(22,FMT="('Sequence number of first shell on each centre')")
  WRITE(22,FMT=int8) ifirst_shll(1:Num_centres+1)

  WRITE(22,FMT="('Exponents of Gaussian primitives')")
  WRITE(22,FMT=reals4) shexpnt

  WRITE(22,FMT="('Normalised contraction coefficients')")
  WRITE(22,FMT=reals4) c_prim

  IF(num_sp .gt. 0)THEN
    WRITE(22,FMT="('2nd contraction coefficients (p coeff. for sp shells,&
                   & 0 otherwise)')")
    WRITE(22,FMT=reals4) c2_prim
  END IF

  WRITE(22,FMT="('Position of each shell (au)')")
  WRITE(22,FMT=reals3) shll_posn

!=================Multideterminant information============================

  WRITE(22,FMT="(/'MULTIDETERMINANT INFORMATION'/&
                 &'----------------------------')")
  IF(CAS)THEN
    CALL cas_write(22)
  ELSE IF(icis_out .eq. 0)THEN
    WRITE(22,FMT="('GS')")
  ELSE
    IF(resummed)THEN
      ! RESUMMED CIS/TD-DFT wavefunction...

      ! No. of determinants in expansion
      Ndet_out=SUM(Ndet(1:ispin_lim))
      IF(.not. SPIN) Ndet_out = 2*Ndet_out
      IF(Ndet_out .gt. 1)THEN
        WRITE(22,FMT="('MD')")
        WRITE(22,*) Ndet_out
      ELSE
        WRITE(22,FMT="('SD')")
      END IF

      ! If we have a spin-restricted wavefunction then we need two
      ! determinants per excitation in order to construct a true spin
      ! state.  We therefore have twice as many determinants and hence
      ! twice as many coefficients
      IF(SPIN)THEN ! Spin-unrestricted case

        IF(ndet_out .gt. 1)THEN
          ! Coefficients of the determinants -  the CIS expansion
          ! coefficients are all included in the eigenvector coefficients
          ! due to the resumming process so we just need some '1's here
          DO ispin=1,ispin_lim,1
            DO i=1,Ndet(ispin),1
              WRITE(22,*) 1.0d0
            END DO
          END DO
        END IF

        ! Corresponding alpha then beta excited state determinants
        idet = 0  ! Counts determinants
        DO ispin=1,ispin_lim,1
          ivirt = 0 ! Count of the effective (resummed) virtual orbitals
                    ! must be reset for different spins
          DO iocc=1,Nspin(ispin),1

            IF(Npromote(iocc,ispin) .gt. 0)THEN
              ivirt = ivirt + 1 ! Next resummed orbital of spin ispin
              idet  = idet  + 1 ! Next excited determinant
              ! Promote electron in band iocc, k-point 1, to band
              ! ivirt+Nspin(ispin) (ivirt is shift above occupied bands),
              ! kpoint 1, in determinant ivirt, spin ispin.  ivirt essentially
              ! counts the no. of (resummed) determinants
              WRITE(22,FMT="('DET ',I3,1x,I1,1x,'PR',2(1x,I3,1x,I1))") &
                               &idet,ispin,iocc,1,(ivirt+Nspin(ispin)),1
            END IF
          END DO
        END DO

      ELSE ! Spin-restricted case

        IF(ndet_out .gt. 1)THEN
          IF(singlet)THEN
            ! Coefficients of the determinants - again, all unity because
            ! of resumming.
            DO i=1,Ndet_out,1
              WRITE(22,*) 1.0d0
            END DO
          ELSE
            DO i=1,(Ndet_out/2),1
              ! Triplet is spin-down excitation minus spin-up excitation
              WRITE(22,*) 1.0d0
              WRITE(22,*) -1.0d0
            END DO
          END IF
        END IF

        ! This is a spin-restricted calculation and therefore we have to
        ! construct pure spin states
        ivirt = 0
        idet = 0  ! Counts determinants
        DO iocc=1,Nspin(1),1

          IF(Npromote(iocc,1) .gt. 0)THEN
            ivirt = ivirt + 1
            IF(singlet)THEN
              ! Construct a spin-singlet excitation
              idet = idet + 1
              WRITE(22,FMT="('DET ',I3,1x,I1,1x,'PR',2(1x,I3,1x,I1))") &
                               &idet,1,iocc,1,(ivirt+Nspin(1)),1
              idet = idet + 1
              WRITE(22,FMT="('DET ',I3,1x,I1,1x,'PR',2(1x,I3,1x,I1))") &
                               &idet,2,iocc,1,(ivirt+Nspin(1)),1
            ELSE
              ! Construct a spin-triplet excitation
! Sz=0 (two-determinant form) triplet - legal for constructing an
! excitation from a ground state with Sz=0
              idet = idet + 1
              WRITE(22,FMT="('DET ',I3,1x,I1,1x,'PR',2(1x,I3,1x,I1))") &
                               &idet,2,iocc,1,(ivirt+Nspin(1)),1
              idet = idet + 1
              WRITE(22,FMT="('DET ',I3,1x,I1,1x,'PR',2(1x,I3,1x,I1))") &
                               &idet,1,iocc,1,(ivirt+Nspin(1)),1
! This bit for an Sz=-1 triplet wavefunction (which has half as many
! determinants as the Sz=0 form).
!              WRITE(22,FMT="('DET ',I3,1x,I1,1x,'MI',1x,I3,1x,I1)") &
!                                              ivirt,1,iocc,1
!              WRITE(22,FMT="('DET ',I3,1x,I1,1x,'PL',1x,I3,1x,I1)") &
!                                              ivirt,2,(ivirt+Nspin(1)),1
            END IF
          END IF
        END DO
      END IF

    ELSE  ! Excited state expansion has NOT been resummed

      ! If we have a spin-restricted wavefunction then can construct
      ! a spin eigenstate and thus have twice as many determinants
      ! and hence twice as many coefficients
      IF(.not. SPIN)THEN
        Ndet_out=2*SUM(Nconfig(1:ispin_lim))
      ELSE
        Ndet_out=SUM(Nconfig(1:ispin_lim))
      END IF

      IF(Ndet_out .gt. 1)THEN
        WRITE(22,FMT="('MD')")
        WRITE(22,*) Ndet_out

        ! CIS expansion coefficients for the alpha (ispin=1) and beta
        ! (ispin=2) excitations...
        DO ispin=1,ispin_lim,1
          DO i=1,Nconfig(ispin),1
            WRITE(22,*) ci_coeff(i,ispin)
            IF(.not. SPIN)THEN
              IF(singlet)THEN
                ! Each singlet is constructed from an alpha and a beta
                ! excitation, each share the same expansion coefficient
                WRITE(22,*) ci_coeff(i,ispin)
              ELSE
                ! Similarly for the triplet except that it is formed
                ! from a spin-down excitation MINUS a spin-up
                WRITE(22,*) -ci_coeff(i,ispin)
              END IF
            END IF
          END DO
        END DO

      ELSE
        WRITE(22,FMT="('SD')")
      END IF

      IF(SPIN)THEN
        ! Corresponding definitions of the alpha & beta excited states
        idet=0
        DO ispin=1,ispin_lim,1
          DO i=1,Nconfig(ispin),1
            idet=idet + 1
            WRITE(22,FMT="('DET ',I3,1x,I1,1x,'PR',2(1x,I3,1x,I1))") idet,ispin, &
                &Orbitals(1,i,ispin),1,Orbitals(2,i,ispin),1
          END DO
        END DO
      ELSE
        ! This is a spin-restricted calculation and therefore we have to
        ! construct pure spin states.  A spin restricted calc. is stored
        ! as an alpha state here so ispin has been set to 1 below.

        idet = 0
        DO i=1,Nconfig(1),1
          iocc=Orbitals(1,i,1)
          ivirt=Orbitals(2,i,1)
          IF(singlet)THEN
            ! Construct a singlet state
            idet=idet + 1
            WRITE(22,FMT="('DET ',I3,1x,I1,1x,'PR',2(1x,I3,1x,I1))") idet,1, &
                &iocc,1,ivirt,1
            idet=idet + 1
            WRITE(22,FMT="('DET ',I3,1x,I1,1x,'PR',2(1x,I3,1x,I1))") idet,2, &
                &iocc,1,ivirt,1
          ELSE
          ! Construct a triplet state
! Sz=0 triplet state
            idet=idet + 1
            WRITE(22,FMT="('DET ',I3,1x,I1,1x,'PR',2(1x,I3,1x,I1))") &
                           &idet,2,iocc,1,ivirt,1
            idet=idet + 1
            WRITE(22,FMT="('DET ',I3,1x,I1,1x,'PR',2(1x,I3,1x,I1))") &
                           &idet,1,iocc,1,ivirt,1
! Sz=-1 triplet state
!            WRITE(22,FMT="('DET ',I3,1x,I1,1x,'MI',1x,I3,1x,I1)") i,1,iocc,1
!            WRITE(22,FMT="('DET ',I3,1x,I1,1x,'PL',1x,I3,1x,I1)") i,2,ivirt,1
          END IF
        END DO
      END IF

    END IF
! The no. of up- and down-spin electrons specified in 'input' will have
! to be modified if we're constructing an Sz!=0 triplet state
!    IF((.not. SPIN) .and. (.not. singlet))THEN
!      WRITE(6,FMT="(/'******************************************************&
!                    &********************'/&
!                    &'WARNING: you must set NEU=',I3,' and NED=',I3,' in the&
!                    & xqmc ''input'' file in order'/'to run this calculation&
!                    & correctly.'/&
!                    &'******************************************************&
!                    &********************')") (Nspin(1)-1),(Nspin(2)+1)
!    END IF
  END IF

!=================Eigenvector components==================================

  WRITE(22,FMT="(/'EIGENVECTOR COEFFICIENTS'/'------------------------')")

  IF(Nmo /= Ncoeffs)THEN
! In G98 Nmo is not always equal to Ncoeffs.  Since CASINO assumes
! that Nmo=Ncoeffs, we have to create some fake MOs to make up the
! difference. To do this, we copy the highest MO (no. Nmo) into the
! (Ncoeffs-Nmo) fake MOs.
    DO ispin=1,ispin_lim,1
      DO iband=(Nmo + 1),Ncoeffs,1
        evcoeff1(1:Ncoeffs,iband,ispin)=evcoeff1(1:Ncoeffs,Nmo,ispin)
      END DO
    END DO
  END IF

  icoeff=0
  DO ispin=1,ispin_lim,1
    DO iband=1,Ncoeffs,1
      DO ifun=1,Ncoeffs,1
        icoeff=icoeff+1
        ck(icoeff)=evcoeff1(ifun,iband,ispin)
! Write the eigenvector coefficients to the file, 100 at a time.
        IF(icoeff .eq. 100)THEN
          icoeff=0
          WRITE(22,FMT=reals4) ck
        END IF
      END DO
    END DO
  END DO

! Allow for the fact that the total number of eigenvector coeff.'s will
! not in general be exactly divisible by 100 so we have some left...
  IF(icoeff .ne. 0)THEN
    WRITE(22,FMT=reals4) ck(1:icoeff)
  END IF

!====================Append any other info. to the input file===========

  WRITE(22,FMT="(/'Job and output files for this calculation (not read &
                  &by QMC program)')")
  WRITE(22,FMT="('=====================================================&
                 &===========================')")

  ! Files to append to the QMC input file
  filenames(1)=TRIM(g94_file)
  filenames(2)=TRIM(g94_file)//".out"

  DO iname=1,2,1

  ! Open the file
    OPEN(unit=23,file=filenames(iname),status='old', IOSTAT=icheck)

    IF(icheck .GT. 0)THEN
      WRITE(*,FMT="('Cannot find the file: ',A30)") filenames(iname)
      WRITE(22,FMT="('Could not find the file, ',A30)") filenames(iname)
    ELSE
      WRITE(22,*)
      scan=.true.
      DO WHILE(scan)

        READ(23,FMT="(A130)",IOSTAT=icheck) line_input

        IF(icheck .eq. 0)THEN
          WRITE(22,FMT=*) TRIM(line_input)
        ELSE IF(icheck .lt. 0)THEN
          scan=.false.
        ELSE
          WRITE(*,FMT="('Error reading ',A30,' ...skipping.')") filenames(iname)
          scan=.false.
        END IF

      END DO
    END IF

    CLOSE(unit=23)

  END DO

  CLOSE(unit=22)

END SUBROUTINE qmc_write

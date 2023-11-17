!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This routine is inspired by the former routine pw2casino of
! Norbert Nemec
! (C) 2010 by Norbert Nemec <Norbert@Nemec-online.de>
!----------------------------------------------------------------------------
SUBROUTINE pw2casino()
  !----------------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  !
  USE mp_global,     ONLY : npool, nimage
  !
  USE control_flags, ONLY : istep, nstep
  !
  USE io_files, ONLY : tmp_dir
  !
  USE plugin_flags, ONLY : use_pw2casino
  !
  IMPLICIT NONE
  !
  CHARACTER(len=4) :: postfix
  !
  CHARACTER(len=6), EXTERNAL :: int_to_char
  !
  INTEGER, EXTERNAL :: find_free_unit
  !
  INTEGER :: tmp_unit
  !
  INTEGER  :: ios
  LOGICAL  :: casino_gather = .true.
  LOGICAL  :: blip_convert = .true.
  LOGICAL  :: blip_binary = .true.
  LOGICAL  :: blip_single_prec = .false.
  REAL(dp) :: blip_multiplicity = 1.d0
  INTEGER  :: n_points_for_test = 0
  !
  NAMELIST / inputpp / &
   blip_convert, &
   blip_multiplicity, &
   blip_binary, &
   blip_single_prec, &
   n_points_for_test
  !
  !
  IF ( use_pw2casino ) THEN
    !
    IF ( npool > 1 .or. nimage > 1) THEN
      CALL errore('pw2casino', 'pool or image parallelization not (yet) implemented',1)
    ENDIF
    !
    tmp_unit = find_free_unit()
    OPEN(unit=tmp_unit,file = trim(tmp_dir)//'/'//'pw2casino.dat',status='old',err=20)
    READ(tmp_unit,inputpp,iostat=ios)
    CLOSE(tmp_unit)

20 CONTINUE

   IF ( .not. blip_convert ) blip_binary = .false.

    IF ( nstep == 1 ) THEN
      write(postfix,*) ''
      CALL write_casino_wfn( &
               casino_gather, & ! gather
               blip_convert,  & ! blip
               blip_multiplicity,   & ! multiplicity
               blip_binary,  & ! binwrite
               blip_single_prec,  & ! single_precision_blips
               n_points_for_test,       & ! n_points_for_test
               postfix)   ! postfix

    ELSE
!      write(postfix,'(i4.4)') istep
      postfix=trim(int_to_char(istep))
      !
      CALL write_casino_wfn( &
               casino_gather, & ! gather
               blip_convert,  & ! blip
               blip_multiplicity,   & ! multiplicity
               blip_binary,  & ! binwrite
               blip_single_prec,  & ! single_precision_blips
               n_points_for_test,       & ! n_points_for_test
               '.'//postfix)   ! postfix
    ENDIF
  ENDIF
  !
  !
END SUBROUTINE pw2casino

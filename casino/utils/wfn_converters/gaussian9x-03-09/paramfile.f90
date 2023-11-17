MODULE paramfile
 IMPLICIT NONE
 INTEGER,PARAMETER :: dp=kind(1.d0)
 REAL(KIND=dp),PARAMETER :: Pi=3.1415926535897932384626433832795d0
 REAL(KIND=dp),PARAMETER :: oneover_pi=1.0d0/Pi
 REAL(KIND=dp),PARAMETER :: oneover_pi_cubed=oneover_pi**3
 REAL(KIND=dp),PARAMETER :: two_rtpi=3.5449077018110320545963349666822d0
! Following two constants taken from Cohen and Taylor, Physics Today, Aug 1994.
 REAL(KIND=dp),PARAMETER :: eV=27.2113962d0
 REAL(KIND=dp),PARAMETER :: Bohr=0.529177249d0
END MODULE paramfile

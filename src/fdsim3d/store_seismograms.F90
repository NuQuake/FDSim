!=======================================================================
! Copyright (c) 2014 Jozef Kristek, Peter Moczo
!
! Licensed under the EUPL (the "Licence")
! You may not use this work except in compliance with the Licence.
! You may obtain a copy of the Licence at http://ec.europa.eu/idabc/eupl
! Unless required by applicable law or agreed to in writing, software
! distributed under the Licence is distributed on an "AS IS" basis,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and
! limitations under the Licence.
!
! REVISION HISTORY:
! Version 1.0, 2014-03-19: Jozef Kristek, Peter Moczo
!=======================================================================

SUBROUTINE  STORE_SEISMOGRAMS(T)

  USE PRECISION   , ONLY: PP

  USE CONTROL_DATA, ONLY: MR, IREC, KREC, LREC, LPAS, REC_NAME, DT

  USE AUXIL       , ONLY: X_MIN ,X_MAX , Y_MIN ,Y_MAX , Z_MIN ,Z_MAX , &
                          X_MINF,X_MAXF, Y_MINF,Y_MAXF, Z_MINF,Z_MAXF, &
                          ORDER

  USE WAVEFIELD   , ONLY: UM, VM, WM, UMF, VMF, WMF, SEISU, SEISV, SEISW

#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: T

  INTEGER :: J, IR, KR, LR

  REAL(PP), DIMENSION (1:MR) :: SU, SV, SW

!-----------------------------------------------------------------------

  INTERFACE
    FUNCTION     INTERPOL(X,Y,Z,A,ORDER)
      USE PRECISION   , ONLY: PP
      REAL(PP),                          INTENT(IN) :: X,Y,Z
      REAL(PP), DIMENSION (1:2,1:2,1:2), INTENT(IN) :: A
      REAL(PP)                                      :: INTERPOL
      INTEGER,                           INTENT(IN) :: ORDER
    END FUNCTION INTERPOL
  END INTERFACE

!-----------------------------------------------------------------------

  SU = 0.
  SV = 0.
  SW = 0.

  DO J = 1, MR

    IR = INT(IREC(J)+0.5)
    KR = INT(KREC(J)    )
    LR = INT(LREC(J)-0.5)
    IF ( LREC(J) < LPAS  ) THEN
      IF ( ( X_MINF <= IR ) .AND. ( IR <= X_MAXF) .AND.                &
           ( Y_MINF <= KR ) .AND. ( KR <= Y_MAXF) .AND.                &
           ( Z_MINF <= LR ) .AND. ( LR <= Z_MAXF)       ) THEN
        SU(J) = INTERPOL( IREC(J)+0.5-IR, KREC(J)-KR, LREC(J)-0.5-LR,  &
                          UMF ( IR:IR+1, KR:KR+1, LR:LR+1 ) ,ORDER )
      END IF
    ELSE
      IF ( ( X_MIN <= IR ) .AND. ( IR <= X_MAX) .AND.                  &
           ( Y_MIN <= KR ) .AND. ( KR <= Y_MAX) .AND.                  &
           ( Z_MIN <= LR ) .AND. ( LR <= Z_MAX)       ) THEN
        SU(J) = INTERPOL( IREC(J)+0.5-IR, KREC(J)-KR, LREC(J)-0.5-LR,  &
                          UM  ( IR:IR+1, KR:KR+1, LR:LR+1 ) ,ORDER )
      END IF
    END IF

    IR = INT(IREC(J)    )
    KR = INT(KREC(J)+0.5)
    LR = INT(LREC(J)-0.5)
    IF ( LREC(J) < LPAS  ) THEN
      IF ( ( X_MINF <= IR ) .AND. ( IR <= X_MAXF) .AND.                &
           ( Y_MINF <= KR ) .AND. ( KR <= Y_MAXF) .AND.                &
           ( Z_MINF <= LR ) .AND. ( LR <= Z_MAXF)       ) THEN
        SV(J) = INTERPOL( IREC(J)-IR, KREC(J)+0.5-KR, LREC(J)-0.5-LR,  &
                          VMF ( IR:IR+1, KR:KR+1, LR:LR+1 ) ,ORDER )
      END IF
    ELSE
      IF ( ( X_MIN <= IR ) .AND. ( IR <= X_MAX) .AND.                  &
           ( Y_MIN <= KR ) .AND. ( KR <= Y_MAX) .AND.                  &
           ( Z_MIN <= LR ) .AND. ( LR <= Z_MAX)       ) THEN
        SV(J) = INTERPOL( IREC(J)-IR, KREC(J)+0.5-KR, LREC(J)-0.5-LR,  &
                          VM  ( IR:IR+1, KR:KR+1, LR:LR+1 ) ,ORDER )
      END IF
    END IF

    IR = INT(IREC(J)    )
    KR = INT(KREC(J)    )
    LR = INT(LREC(J)    )
    IF ( LREC(J) < LPAS  ) THEN
      IF ( ( X_MINF <= IR ) .AND. ( IR <= X_MAXF) .AND.                &
           ( Y_MINF <= KR ) .AND. ( KR <= Y_MAXF) .AND.                &
           ( Z_MINF <= LR ) .AND. ( LR <= Z_MAXF)       ) THEN
        SW(J) = INTERPOL( IREC(J)-IR, KREC(J)-KR, LREC(J)-LR,          &
                          WMF ( IR:IR+1, KR:KR+1, LR:LR+1 ) ,ORDER )
      END IF
    ELSE
      IF ( ( X_MIN <= IR ) .AND. ( IR <= X_MAX) .AND.                  &
           ( Y_MIN <= KR ) .AND. ( KR <= Y_MAX) .AND.                  &
           ( Z_MIN <= LR ) .AND. ( LR <= Z_MAX)       ) THEN
        SW(J) = INTERPOL( IREC(J)-IR, KREC(J)-KR, LREC(J)-LR,          &
                          WM  ( IR:IR+1, KR:KR+1, LR:LR+1 ) ,ORDER )
      END IF
    END IF

  END DO

#if(USE_MPI)
  CALL MPI_REDUCE                                                      &
           (SU, SEISU, MR, MPI_REAL, MPI_SUM, root, comm_cart, errcode)
  CALL MPI_REDUCE                                                      &
           (SV, SEISV, MR, MPI_REAL, MPI_SUM, root, comm_cart, errcode)
  CALL MPI_REDUCE                                                      &
           (SW, SEISW, MR, MPI_REAL, MPI_SUM, root, comm_cart, errcode)


  IF ( myid == root ) THEN
    DO J = 1, MR
      OPEN ( 21, FILE=REC_NAME(J)//'.asc', STATUS='UNKNOWN',           &
                                                     POSITION='APPEND')
      WRITE (21, '(4(1X,E13.6))') DT*T, SEISU (J), SEISV (J), -SEISW (J)
      CLOSE (21)
    END DO
  END IF
#else
  DO J = 1, MR
    OPEN ( 21, FILE=REC_NAME(J)//'.asc', STATUS='UNKNOWN',           &
                                                   POSITION='APPEND')
    WRITE (21, '(4(1X,E13.6))') DT*T, SU (J), SV (J), -SW (J)
    CLOSE (21)
  END DO
#endif


END SUBROUTINE  STORE_SEISMOGRAMS

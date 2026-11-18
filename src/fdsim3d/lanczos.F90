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

SUBROUTINE  LANCZOS

  USE PRECISION   , ONLY: PP

  USE CONTROL_DATA, ONLY: REFINE

  USE AUXIL       , ONLY: LANC, REF2, PI, LLAN,                        \
                           X_MIN_DF, X_MAX_DF,                         &
                           Y_MIN_DF, Y_MAX_DF,                         &
                           VARF1

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER  :: I, K

  REAL(PP) :: X, Y, R, SUMA, SINC_X, SINC_Y

!-----------------------------------------------------------------------

  ALLOCATE ( VARF1 ( X_MIN_DF-LLAN*REFINE:X_MAX_DF+LLAN*REFINE,     &
                     Y_MIN_DF-LLAN*REFINE:Y_MAX_DF+LLAN*REFINE ) )

  REF2 = LLAN*REFINE

  ALLOCATE ( LANC ( -REF2:REF2, -REF2:REF2 ) )

  SUMA = 0.

  DO I = -REF2, REF2
    DO K = -REF2, REF2

      X = REAL(I)/REFINE
      Y = REAL(K)/REFINE
      R = SQRT(X*X+Y*Y)

      IF ( R < REAL(LLAN) ) THEN
        IF ( R > 0. ) THEN
          IF ( I == 0 ) THEN
            SINC_X = 1.
          ELSE
            SINC_X = SIN (PI*X) / PI / X
          END IF
          IF ( K == 0 ) THEN
            SINC_Y = 1.
          ELSE
            SINC_Y = SIN (PI*Y) / PI / Y
          END IF

          LANC (I,K) = SINC_X * SINC_Y * SIN (PI*R) / (PI*R/REAL(LLAN))
        ELSE
          LANC (I,K) = 1.
        END IF
      ELSE
        LANC (I,K) = 0.
      END IF

      SUMA = SUMA + LANC(I,K)

    END DO
  END DO

  LANC = LANC / SUMA

END SUBROUTINE  LANCZOS

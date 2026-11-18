!=======================================================================
! Copyright (c) 2014 Jozef Kristek, Peter Moczo, Martin Galis
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
! Version 1.0, 2014-03-19: Jozef Kristek, Peter Moczo, Martin Galis
!=======================================================================
MODULE MOD_NR

USE NRTYPE

! POLINT: piecewise-linear interpolation in a monotonic table
! (own implementation; only the routine name follows Numerical Recipes)

IMPLICIT NONE

!-------------------------------------------------------- implementation
CONTAINS

!---------------------------------------------------------------- POLINT
  SUBROUTINE POLINT(XA,YA,X,Y,N)

  IMPLICIT NONE

  INTEGER            :: N, I

  REAL(DP), DIMENSION(:) :: XA, YA
  REAL(DP)               :: X , Y

  Y=0

  IF ( XA(1)<XA(N) ) THEN

    DO I = 1, N-1

      IF      ( X == XA(I  ) ) THEN
        Y = YA(I)
        RETURN
      ELSE IF ( X == XA(I+1) ) THEN
        Y = YA(I+1)
        RETURN
      ELSE IF ( X <  XA(I+1) ) THEN
        Y = (YA(I+1)-YA(I))/(XA(I+1)-XA(I))*(X-XA(I))+YA(I)
        RETURN
      END IF

    END DO
  ELSE
    DO I = 1, N-1

      IF      ( X == XA(I  ) ) THEN
        Y = YA(I)
        RETURN
      ELSE IF ( X == XA(I+1) ) THEN
        Y = YA(I+1)
        RETURN
      ELSE IF ( X >  XA(I+1) ) THEN
        Y = (YA(I+1)-YA(I))/(XA(I+1)-XA(I))*(X-XA(I))+YA(I)
        RETURN
      END IF

    END DO
  END IF

  END SUBROUTINE POLINT

END MODULE MOD_NR

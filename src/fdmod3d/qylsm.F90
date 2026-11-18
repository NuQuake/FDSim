!=======================================================================
! Copyright (c) 2014 Jozef Kristek, Peter Moczo, Martin Galis
!
! Licensed under the EUPL (the "Licence")
! You may not use this work except in compliance with the Licence.
! You may obtain a copy of the Licence at
! https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12
! Unless required by applicable law or agreed to in writing, software
! distributed under the Licence is distributed on an "AS IS" basis,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and
! limitations under the Licence.
!
! REVISION HISTORY:
! Version 1.0, 2014-03-19: Jozef Kristek, Peter Moczo, Martin Galis
!=======================================================================
SUBROUTINE QYLSM(Q,Y)

  USE NRTYPE   , ONLY: WP
  USE MOD_CONST, ONLY: N_FREQ, I2N_FREQ1
  USE MOD_ATTEN, ONLY: FJJ, FJK, A, B, P

!----------------------------------------------------------------------

  IMPLICIT NONE

  REAL(WP), DIMENSION (  N_FREQ ), INTENT(INOUT) :: Y
  REAL(WP), DIMENSION (I2N_FREQ1), INTENT(INOUT) :: Q

  INTEGER  :: I, J, K
  REAL(WP) :: FJ, FI, SUMN

!----------------------------------------------------------------------

  IF ( MINVAL(Q) > 1000. ) THEN
     Y = 0.
  ELSE

    DO J = 1, N_FREQ
      B(J) = 0.
      A(J,:) = 0.
      DO K = 1, I2N_FREQ1
        FJ = FJK(J,K) + FJJ(J,K)/Q(K)
        B(J) = B(J)+ FJ/Q(K)
        DO I = 1, N_FREQ
          FI = FJK(I,K) + FJJ(I,K)/Q(K)
          A(J,I) = A(J,I) + FI*FJ
        END DO
      END DO
    END DO

    DO I = 1, N_FREQ
      SUMN = A(I,I)-DOT_PRODUCT(A(I,1:I-1),A(I,1:I-1))
      IF (SUMN <= 0.0) THEN
        PRINT *, 'CHOLECKY DECOMPOSITION FAILED'
        STOP
      END IF
      P(I)=SQRT(SUMN)
      A(I+1:N_FREQ,I)= ( A(I,I+1:N_FREQ)                                 &
                       - MATMUL(A(I+1:N_FREQ,1:I-1),A(I,1:I-1)) ) / P(I)
    END DO

    DO I = 1, N_FREQ
      Y(I)=(B(I)-DOT_PRODUCT(A(I,1:I-1),Y(1:I-1)))/P(I)
    END DO
    DO I = N_FREQ, 1, -1
      Y(I)=(Y(I)-DOT_PRODUCT(A(I+1:N_FREQ,I),Y(I+1:N_FREQ)))/P(I)
    END DO

  END IF

END SUBROUTINE QYLSM

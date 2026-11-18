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

SUBROUTINE  COARSE2FINE_TXZ

  USE PRECISION   , ONLY: PP

  USE CONTROL_DATA, ONLY: LPAS, REFINE

  USE GRID_MEDIUM , ONLY: H

  USE STRESSFIELD , ONLY: TXZF, TXZ

  USE AUXIL       , ONLY:  X_MIN_D , X_MAX_D , X_MIN , X_MAX ,         &
                           Y_MIN_D , Y_MAX_D , Y_MIN , Y_MAX , PMLL

  USE INTERFACES  , ONLY: BCUCOF

!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL(PP), DIMENSION (4)   :: Y, Y1, Y2, Y12
  REAL(PP), DIMENSION (4,4) :: C

  REAL(PP), DIMENSION (X_MIN-1:X_MAX+1,Y_MIN-1:Y_MAX) :: Y1A, Y2A, Y12A


  INTEGER :: I, K, J, JI, JK, IL, I0, IR, KL, KR, JI_MIN, JI_MAX,      &
             JK_MIN, ISTART, KSTART
  REAL(PP):: U, T, ANSY, DX, DY, TXZ_IR_KL, TXZ_IR_K, TXZ_IR_KR

!-----------------------------------------------------------------------

  ISTART = X_MIN
  KSTART = Y_MIN

  IF ( X_MIN == X_MIN_D ) THEN
    ISTART = ISTART+1
  END IF
  IF ( Y_MIN == Y_MIN_D ) THEN
    KSTART = KSTART+1
  END IF

  DO K = KSTART-1, Y_MAX
    DO I = ISTART-1, X_MAX+1

      IR = I + 1
      IL = I - 1
      KR = K + 1
      KL = K - 1
      DX = 2.*H
      DY = 2.*H
      IF ( IL < X_MIN_D ) THEN
        IL = I
        DX = H
      END IF
      IF ( KR > Y_MAX_D ) THEN
        KR = K
        DY = H
      END IF
      IF ( KL < Y_MIN_D ) THEN
        KL = K
        DY = H
      END IF
      IF ( IR > X_MAX_D ) THEN
        I0 = I
        IF ( IR > X_MAX_D+1 ) THEN
          I0 = I -1
          IL = IL-1
        END IF
!        TXZ_IR_KL= 2.*TXZ(X_MAX,KL,LPAS) - TXZ (X_MAX-1,KL,LPAS)
!        TXZ_IR_K = 2.*TXZ(X_MAX,K ,LPAS) - TXZ (X_MAX-1,K ,LPAS)
!        TXZ_IR_KR= 2.*TXZ(X_MAX,KR,LPAS) - TXZ (X_MAX-1,KR,LPAS)
        TXZ_IR_KL= 0.
        TXZ_IR_K = 0.
        TXZ_IR_KR= 0.
        Y1A (I,K) = ( TXZ_IR_K        - TXZ(IL,K ,LPAS) ) /  DX
        Y2A (I,K) = ( TXZ(I0,KR,LPAS) - TXZ(I0,KL,LPAS) ) /  DY
        Y12A(I,K) = ( TXZ_IR_KR       - TXZ_IR_KL                      &
                    - TXZ(IL,KR,LPAS) + TXZ(IL,KL,LPAS) ) / (DX*DY)
      ELSE
        Y1A (I,K) = ( TXZ(IR,K ,LPAS) - TXZ(IL,K ,LPAS) ) /  DX
        Y2A (I,K) = ( TXZ(I ,KR,LPAS) - TXZ(I ,KL,LPAS) ) /  DY
        Y12A(I,K) = ( TXZ(IR,KR,LPAS) - TXZ(IR,KL,LPAS)                &
                    - TXZ(IL,KR,LPAS) + TXZ(IL,KL,LPAS) ) / (DX*DY)
      END IF

      IF ( I == X_MIN_D) THEN
        IF ( (K == Y_MIN_D ) .OR. ( K == Y_MAX_D ) ) THEN
          Y1A (I,K) = 0.
          Y2A (I,K) = 0.
          Y12A(I,K) = 0.
        END IF
      END IF
      IF ( I >= X_MAX_D ) THEN
        IF ( (K == Y_MIN_D ) .OR. ( K == Y_MAX_D ) ) THEN
          Y1A (I,K) = 0.
          Y2A (I,K) = 0.
          Y12A(I,K) = 0.
        END IF
      END IF

    END DO
  END DO

  DO K = KSTART, Y_MAX
    DO I = ISTART, X_MAX+1

      Y  (1) = TXZ  (I-1, K-1, LPAS)
      Y1 (1) = Y1A  (I-1, K-1)
      Y2 (1) = Y2A  (I-1, K-1)
      Y12(1) = Y12A (I-1, K-1)

      IF ( I == X_MAX_D+1 ) THEN
!        Y(2) = 2.*TXZ(X_MAX,K-1,LPAS) - TXZ (X_MAX-1,K-1,LPAS)
        Y(2) = 0.
      ELSE
        Y(2) = TXZ  (I  , K-1, LPAS)
      END IF
      Y1 (2) = Y1A  (I  , K-1)
      Y2 (2) = Y2A  (I  , K-1)
      Y12(2) = Y12A (I  , K-1)

      IF ( I == X_MAX_D+1 ) THEN
!        Y(3) = 2.*TXZ(X_MAX,K  ,LPAS) - TXZ (X_MAX-1,K  ,LPAS)
        Y(3) = 0.
      ELSE
        Y(3) = TXZ  (I  , K  , LPAS)
      END IF
      Y1 (3) = Y1A  (I  , K  )
      Y2 (3) = Y2A  (I  , K  )
      Y12(3) = Y12A (I  , K  )

      Y  (4) = TXZ  (I-1, K  , LPAS)
      Y1 (4) = Y1A  (I-1, K  )
      Y2 (4) = Y2A  (I-1, K  )
      Y12(4) = Y12A (I-1, K  )

      CALL BCUCOF (Y,Y1,Y2,Y12, H, H, C)

      IF ( ( K == KSTART ) .AND. ( Y_MIN == Y_MIN_D ) ) THEN
        JK_MIN = 0
      ELSE
        JK_MIN = 1
      END IF

      DO JK = JK_MIN, REFINE

        IF ( I == ISTART ) THEN
            JI_MIN = (REFINE+1)/2
        ELSE
            JI_MIN = 1
        END IF

        IF ( I == X_MAX+1 ) THEN
          JI_MAX = (REFINE+1)/2 - 1
        ELSE
          JI_MAX = REFINE
        END IF

        DO JI = JI_MIN, JI_MAX

          T = REAL(JI)/REAL(REFINE)
          U = REAL(JK)/REAL(REFINE)
          ANSY = 0.
          DO J = 4, 1, -1
            ANSY=T*ANSY+((C(J,4)*U+C(J,3))*U+C(J,2))*U+C(J,1)
          END DO

          TXZF( ((2*I-5)*REFINE+3)/2+JI, (K-2)*REFINE+1+JK, LPAS) = ANSY

        END DO
      END DO

    END DO
  END DO


END SUBROUTINE  COARSE2FINE_TXZ

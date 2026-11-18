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

SUBROUTINE  STRESS_XY (L    , MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,  &
                                          X_MIN, X_MAX, Y_MIN, Y_MAX,  &
                                          TPML, TXYT, EXYT, XXYT, JMT, &
                           XXYK2T, XXYKMYT, XXYI2T, XXYIMXT )

  USE PRECISION   , ONLY: PP, IK
  USE CONTROL_DATA, ONLY: DT, TAU_EPS
  USE GRID_MEDIUM , ONLY: MUXY, YMUXY
  USE AUXIL       , ONLY: OTA, OTB
  USE INTERFACES  , ONLY: IND

!-----------------------------------------------------------------------

  IMPLICIT NONE


  INTEGER,                                           INTENT (IN)    :: &
                       MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH, L, TPML,&
                       X_MIN, X_MAX, Y_MIN, Y_MAX

  INTEGER(IK), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                   JMT

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EXYT

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XXYT

  REAL   (PP), DIMENSION (MXL:MXH,       4, MZL:MZH), INTENT(INOUT) :: &
                                                      XXYK2T, XXYKMYT

  REAL   (PP), DIMENSION (      MYL:MYH, 4, MZL:MZH), INTENT(INOUT) :: &
                                                      XXYI2T, XXYIMXT

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  TXYT

!-----------------------------------------------------------------------

  INTEGER       :: I, K,           JM1,  IO,                           &
                   IN_I_K_, IN_I1K_, IN_I_K1, IN_I1K1

  REAL (KIND=PP):: XY_I_K_, XY_I1K_, XY_I_K1, XY_I1K1, TXY, XXY

!-----------------------------------------------------------------------

!========================= TXYT ===========================
!==== 2 ====
K = 2 - TPML
IF ( ( Y_MIN <= K ) .AND. ( K <= Y_MAX ) ) THEN
  DO  I = MAX(2 - TPML,X_MIN), MIN(MXT + TPML,X_MAX)

    JM1     = JMT (I,K, L )

    TXY =  MUXY(JM1)         * EXYT(I,K)

    DO IO = 1, 4
      XXY = OTA(IO)*EXYT(I,K) + OTB(IO)*XXYK2T(I,IO, L )
      TXY = TXY - YMUXY(JM1,IO) *  XXY
      XXYK2T(I,IO, L ) = XXY
    END DO

    TXYT(I,K, L) = TXYT(I,K, L) + ( DT + TAU_EPS ) * TXY

  END DO
END IF


!==== 3 - MYT-1 ====
DO  K = MAX(3-TPML,Y_MIN), MIN(MYT-1+TPML,Y_MAX)
  I = 2 - TPML
  IF ( ( X_MIN <= I ) .AND. ( I <= X_MAX) ) THEN
    JM1     = JMT (I,K, L )

    TXY =  MUXY(JM1)         * EXYT(I,K)

    DO IO = 1, 4
      XXY = OTA(IO)*EXYT(I,K) + OTB(IO)*XXYI2T(K,IO, L )
      TXY = TXY - YMUXY(JM1,IO) *  XXY
      XXYI2T(K,IO, L ) = XXY
    END DO

    TXYT(I,K, L) = TXYT(I,K, L) + ( DT + TAU_EPS ) * TXY
  END IF

  DO  I = MAX(3-TPML,X_MIN), MIN(MXT-1+TPML,X_MAX)

    JM1     = JMT (I,K, L )

    IN_I_K_ = IND (I  ,K  )
    IN_I1K_ = IND (I-1,K  )
    IN_I_K1 = IND (I  ,K-1)
    IN_I1K1 = IND (I-1,K-1)

    XY_I_K_ =   XXYT(I  ,K  , L )
    XY_I1K_ = ( XXYT(I-1,K  , L ) + XXYT(I+1,K  , L ) ) * 0.50_PP
    XY_I_K1 = ( XXYT(I  ,K-1, L ) + XXYT(I  ,K+1, L ) ) * 0.50_PP
    XY_I1K1 = ( XXYT(I-1,K-1, L ) + XXYT(I+1,K-1, L )                  &
            +   XXYT(I-1,K+1, L ) + XXYT(I+1,K+1, L ) ) * 0.25_PP

    TXYT(I,K, L)    = TXYT(I,K, L) + ( DT + TAU_EPS ) * (                            &
                       MUXY(JM1)         * EXYT(I,K   )                &
                    - YMUXY(JM1,IN_I_K_) *  XY_I_K_                    &
                    - YMUXY(JM1,IN_I1K_) *  XY_I1K_                    &
                    - YMUXY(JM1,IN_I_K1) *  XY_I_K1                    &
                    - YMUXY(JM1,IN_I1K1) *  XY_I1K1      )
  END DO

  I = MXT + TPML
  IF ( ( X_MIN <= I ) .AND. ( I <= X_MAX) ) THEN
    JM1     = JMT (I,K, L )

    TXY =  MUXY(JM1)         * EXYT(I,K)

    DO IO = 1, 4
      XXY = OTA(IO)*EXYT(I,K) + OTB(IO)*XXYIMXT(K,IO, L )
      TXY = TXY - YMUXY(JM1,IO) *  XXY
      XXYIMXT(K,IO, L ) = XXY
    END DO

    TXYT(I,K, L) = TXYT(I,K, L) + ( DT + TAU_EPS ) * TXY
  END IF

END DO


!==== MYT ====
K = MYT + TPML
IF ( ( Y_MIN <= K ) .AND. ( K <= Y_MAX) ) THEN
  DO  I = MAX(2-TPML,X_MIN), MIN(MXT+TPML,X_MAX)

    JM1     = JMT (I,K, L )

    TXY =  MUXY(JM1)         * EXYT(I,K)

    DO IO = 1, 4
      XXY = OTA(IO)*EXYT(I,K) + OTB(IO)*XXYKMYT(I,IO, L )
      TXY = TXY - YMUXY(JM1,IO) *  XXY
      XXYKMYT(I,IO, L ) = XXY
    END DO

    TXYT(I,K, L) = TXYT(I,K, L) + ( DT + TAU_EPS ) * TXY

  END DO
END IF

END SUBROUTINE  STRESS_XY

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

SUBROUTINE  STRESS_NN (L    , MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,  &
                     X_MIN, X_MAX, Y_MIN, Y_MAX,                       &
                     TPML, TXXT,       TYYT,       TZZT,               &
                           XXXT, EXXT, XYYT, EYYT, XZZT, EZZT, JMT,    &
                           XXXK2T, XXXI2T, XXXKMY1T, XXXIMX1T,         &
                           XYYK2T, XYYI2T, XYYKMY1T, XYYIMX1T,         &
                           XZZK2T, XZZI2T, XZZKMY1T, XZZIMX1T )

  USE PRECISION   , ONLY: PP, IK
  USE CONTROL_DATA, ONLY: DT, TAU_EPS
  USE GRID_MEDIUM , ONLY: L2MX , L2MY , L2MZ , LAXY , LAXZ , LAYZ ,    &
                          YL2MX , YL2MY , YL2MZ ,                      &
                          YLAXY , YLAXZ , YLAYZ
  USE AUXIL       , ONLY: OTA, OTB
  USE INTERFACES  , ONLY: IND

!-----------------------------------------------------------------------

  IMPLICIT NONE


  INTEGER,                                            INTENT(IN)    :: &
                       MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH, L, TPML,&
                       X_MIN, X_MAX, Y_MIN, Y_MAX

  INTEGER(IK), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                   JMT

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                      EXXT, EYYT, EZZT

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                      XXXT, XYYT, XZZT

  REAL   (PP), DIMENSION (MXL:MXH,       4, MZL:MZH), INTENT(INOUT) :: &
                  XXXK2T, XXXKMY1T, XYYK2T, XYYKMY1T, XZZK2T, XZZKMY1T

  REAL   (PP), DIMENSION (      MYL:MYH, 4, MZL:MZH), INTENT(INOUT) :: &
                  XXXI2T, XXXIMX1T, XYYI2T, XYYIMX1T, XZZI2T, XZZIMX1T

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                      TXXT, TYYT, TZZT

!-----------------------------------------------------------------------

  INTEGER       :: I, K,               IO,                             &
                   JM1, IN_I_K_, IN_I1K_, IN_I_K1, IN_I1K1

  REAL (KIND=PP):: XX_I_K_, XX_I1K_, XX_I_K1, XX_I1K1,                 &
                   YY_I_K_, YY_I1K_, YY_I_K1, YY_I1K1,                 &
                   ZZ_I_K_, ZZ_I1K_, ZZ_I_K1, ZZ_I1K1,                 &
                   XXX, XYY, XZZ, TXX, TYY, TZZ

!-----------------------------------------------------------------------

!========================= TXXT, TYYT, TZZT ===========================
!==== 2 ====
K = 2 - TPML
IF ( ( Y_MIN <= K ) .AND. ( K <= Y_MAX ) ) THEN
  DO  I = MAX(2 - TPML,X_MIN), MIN(MXT-1 + TPML,X_MAX)

    JM1     = JMT (I,K, L)

    TXX             =  L2MX(JM1)         * EXXT(I,K   )                &
                    +  LAXY(JM1)         * EYYT(I,K   )                &
                    +  LAXZ(JM1)         * EZZT(I,K   )

    TYY             =  LAXY(JM1)         * EXXT(I,K   )                &
                    +  L2MY(JM1)         * EYYT(I,K   )                &
                    +  LAYZ(JM1)         * EZZT(I,K   )

    TZZ             =  LAXZ(JM1)         * EXXT(I,K   )                &
                    +  LAYZ(JM1)         * EYYT(I,K   )                &
                    +  L2MZ(JM1)         * EZZT(I,K   )

    DO IO = 1, 4
      XXX = OTA(IO)*EXXT(I,K) + OTB(IO)*XXXK2T(I,IO, L)
      XYY = OTA(IO)*EYYT(I,K) + OTB(IO)*XYYK2T(I,IO, L)
      XZZ = OTA(IO)*EZZT(I,K) + OTB(IO)*XZZK2T(I,IO, L)
      TXX = TXX - YL2MX(JM1,IO) *  XXX                                 &
                - YLAXY(JM1,IO) *  XYY                                 &
                - YLAXZ(JM1,IO) *  XZZ
      TYY = TYY - YLAXY(JM1,IO) *  XXX                                 &
                - YL2MY(JM1,IO) *  XYY                                 &
                - YLAYZ(JM1,IO) *  XZZ
      TZZ = TZZ - YLAXZ(JM1,IO) *  XXX                                 &
                - YLAYZ(JM1,IO) *  XYY                                 &
                - YL2MZ(JM1,IO) *  XZZ
      XXXK2T(I,IO, L) = XXX
      XYYK2T(I,IO, L) = XYY
      XZZK2T(I,IO, L) = XZZ
    END DO

    TXXT(I,K, L) = TXXT(I,K, L) + ( DT + TAU_EPS) * TXX
    TYYT(I,K, L) = TYYT(I,K, L) + ( DT + TAU_EPS) * TYY
    TZZT(I,K, L) = TZZT(I,K, L) + ( DT + TAU_EPS) * TZZ

  END DO
END IF


!==== 3 - MYT-2 ====
DO  K = MAX(3 - TPML,Y_MIN), MIN(MYT-2 + TPML,Y_MAX)

  I = 2 - TPML
  IF ( ( X_MIN <= I ) .AND. ( I <= X_MAX) ) THEN

    JM1     = JMT (I,K, L)

    TXX             =  L2MX(JM1)         * EXXT(I,K   )                &
                    +  LAXY(JM1)         * EYYT(I,K   )                &
                    +  LAXZ(JM1)         * EZZT(I,K   )

    TYY             =  LAXY(JM1)         * EXXT(I,K   )                &
                    +  L2MY(JM1)         * EYYT(I,K   )                &
                    +  LAYZ(JM1)         * EZZT(I,K   )

    TZZ             =  LAXZ(JM1)         * EXXT(I,K   )                &
                    +  LAYZ(JM1)         * EYYT(I,K   )                &
                    +  L2MZ(JM1)         * EZZT(I,K   )

    DO IO = 1, 4
      XXX = OTA(IO)*EXXT(I,K) + OTB(IO)*XXXI2T(K,IO, L)
      XYY = OTA(IO)*EYYT(I,K) + OTB(IO)*XYYI2T(K,IO, L)
      XZZ = OTA(IO)*EZZT(I,K) + OTB(IO)*XZZI2T(K,IO, L)
      TXX = TXX - YL2MX(JM1,IO) *  XXX                                 &
                - YLAXY(JM1,IO) *  XYY                                 &
                - YLAXZ(JM1,IO) *  XZZ
      TYY = TYY - YLAXY(JM1,IO) *  XXX                                 &
                - YL2MY(JM1,IO) *  XYY                                 &
                - YLAYZ(JM1,IO) *  XZZ
      TZZ = TZZ - YLAXZ(JM1,IO) *  XXX                                 &
                - YLAYZ(JM1,IO) *  XYY                                 &
                - YL2MZ(JM1,IO) *  XZZ
      XXXI2T(K,IO, L) = XXX
      XYYI2T(K,IO, L) = XYY
      XZZI2T(K,IO, L) = XZZ
    END DO

    TXXT(I,K, L) = TXXT(I,K, L) + ( DT + TAU_EPS) * TXX
    TYYT(I,K, L) = TYYT(I,K, L) + ( DT + TAU_EPS) * TYY
    TZZT(I,K, L) = TZZT(I,K, L) + ( DT + TAU_EPS) * TZZ
  END IF

  DO  I = MAX(3 - TPML,X_MIN), MIN(MXT-2 + TPML,X_MAX)

    JM1     = JMT (I,K, L)

    IN_I_K_ = IND (I  ,K  )
    IN_I1K_ = IND (I-1,K  )
    IN_I_K1 = IND (I  ,K-1)
    IN_I1K1 = IND (I-1,K-1)

    XX_I_K_ =   XXXT(I  ,K  , L)
    XX_I1K_ = ( XXXT(I-1,K  , L) + XXXT(I+1,K  , L) ) * 0.50_PP
    XX_I_K1 = ( XXXT(I  ,K-1, L) + XXXT(I  ,K+1, L) ) * 0.50_PP
    XX_I1K1 = ( XXXT(I-1,K-1, L) + XXXT(I+1,K-1, L)                    &
            +   XXXT(I-1,K+1, L) + XXXT(I+1,K+1, L) ) * 0.25_PP

    YY_I_K_ =   XYYT(I  ,K  , L)
    YY_I1K_ = ( XYYT(I-1,K  , L) + XYYT(I+1,K  , L) ) * 0.50_PP
    YY_I_K1 = ( XYYT(I  ,K-1, L) + XYYT(I  ,K+1, L) ) * 0.50_PP
    YY_I1K1 = ( XYYT(I-1,K-1, L) + XYYT(I+1,K-1, L)                    &
            +   XYYT(I-1,K+1, L) + XYYT(I+1,K+1, L) ) * 0.25_PP

    ZZ_I_K_ =   XZZT(I  ,K  , L)
    ZZ_I1K_ = ( XZZT(I-1,K  , L) + XZZT(I+1,K  , L) ) * 0.50_PP
    ZZ_I_K1 = ( XZZT(I  ,K-1, L) + XZZT(I  ,K+1, L) ) * 0.50_PP
    ZZ_I1K1 = ( XZZT(I-1,K-1, L) + XZZT(I+1,K-1, L)                    &
            +   XZZT(I-1,K+1, L) + XZZT(I+1,K+1, L) ) * 0.25_PP

    TXXT(I,K, L) = TXXT(I,K, L) + ( DT + TAU_EPS) * (                               &
                       L2MX(JM1)         * EXXT(I,K   )                &
                    - YL2MX(JM1,IN_I_K_) *  XX_I_K_                    &
                    - YL2MX(JM1,IN_I1K_) *  XX_I1K_                    &
                    - YL2MX(JM1,IN_I_K1) *  XX_I_K1                    &
                    - YL2MX(JM1,IN_I1K1) *  XX_I1K1                    &
                    +  LAXY(JM1)         * EYYT(I,K   )                &
                    - YLAXY(JM1,IN_I_K_) *  YY_I_K_                    &
                    - YLAXY(JM1,IN_I1K_) *  YY_I1K_                    &
                    - YLAXY(JM1,IN_I_K1) *  YY_I_K1                    &
                    - YLAXY(JM1,IN_I1K1) *  YY_I1K1                    &
                    +  LAXZ(JM1)         * EZZT(I,K   )                &
                    - YLAXZ(JM1,IN_I_K_) *  ZZ_I_K_                    &
                    - YLAXZ(JM1,IN_I1K_) *  ZZ_I1K_                    &
                    - YLAXZ(JM1,IN_I_K1) *  ZZ_I_K1                    &
                    - YLAXZ(JM1,IN_I1K1) *  ZZ_I1K1      )

    TYYT(I,K, L) = TYYT(I,K, L) + ( DT + TAU_EPS) * (                               &
                       LAXY(JM1)         * EXXT(I,K   )                &
                    - YLAXY(JM1,IN_I_K_) *  XX_I_K_                    &
                    - YLAXY(JM1,IN_I1K_) *  XX_I1K_                    &
                    - YLAXY(JM1,IN_I_K1) *  XX_I_K1                    &
                    - YLAXY(JM1,IN_I1K1) *  XX_I1K1                    &
                    +  L2MY(JM1)         * EYYT(I,K   )                &
                    - YL2MY(JM1,IN_I_K_) *  YY_I_K_                    &
                    - YL2MY(JM1,IN_I1K_) *  YY_I1K_                    &
                    - YL2MY(JM1,IN_I_K1) *  YY_I_K1                    &
                    - YL2MY(JM1,IN_I1K1) *  YY_I1K1                    &
                    +  LAYZ(JM1)         * EZZT(I,K   )                &
                    - YLAYZ(JM1,IN_I_K_) *  ZZ_I_K_                    &
                    - YLAYZ(JM1,IN_I1K_) *  ZZ_I1K_                    &
                    - YLAYZ(JM1,IN_I_K1) *  ZZ_I_K1                    &
                    - YLAYZ(JM1,IN_I1K1) *  ZZ_I1K1      )

    TZZT(I,K, L) = TZZT(I,K, L) + ( DT + TAU_EPS) * (                               &
                       LAXZ(JM1)         * EXXT(I,K   )                &
                    - YLAXZ(JM1,IN_I_K_) *  XX_I_K_                    &
                    - YLAXZ(JM1,IN_I1K_) *  XX_I1K_                    &
                    - YLAXZ(JM1,IN_I_K1) *  XX_I_K1                    &
                    - YLAXZ(JM1,IN_I1K1) *  XX_I1K1                    &
                    +  LAYZ(JM1)         * EYYT(I,K   )                &
                    - YLAYZ(JM1,IN_I_K_) *  YY_I_K_                    &
                    - YLAYZ(JM1,IN_I1K_) *  YY_I1K_                    &
                    - YLAYZ(JM1,IN_I_K1) *  YY_I_K1                    &
                    - YLAYZ(JM1,IN_I1K1) *  YY_I1K1                    &
                    +  L2MZ(JM1)         * EZZT(I,K   )                &
                    - YL2MZ(JM1,IN_I_K_) *  ZZ_I_K_                    &
                    - YL2MZ(JM1,IN_I1K_) *  ZZ_I1K_                    &
                    - YL2MZ(JM1,IN_I_K1) *  ZZ_I_K1                    &
                    - YL2MZ(JM1,IN_I1K1) *  ZZ_I1K1      )

  END DO

  I = MXT-1 + TPML
  IF ( ( X_MIN <= I ) .AND. ( I <= X_MAX) ) THEN

    JM1     = JMT (I,K, L)

    TXX             =  L2MX(JM1)         * EXXT(I,K   )                &
                    +  LAXY(JM1)         * EYYT(I,K   )                &
                    +  LAXZ(JM1)         * EZZT(I,K   )

    TYY             =  LAXY(JM1)         * EXXT(I,K   )                &
                    +  L2MY(JM1)         * EYYT(I,K   )                &
                    +  LAYZ(JM1)         * EZZT(I,K   )

    TZZ             =  LAXZ(JM1)         * EXXT(I,K   )                &
                    +  LAYZ(JM1)         * EYYT(I,K   )                &
                    +  L2MZ(JM1)         * EZZT(I,K   )

    DO IO = 1, 4
      XXX = OTA(IO)*EXXT(I,K) + OTB(IO)*XXXIMX1T(K,IO, L)
      XYY = OTA(IO)*EYYT(I,K) + OTB(IO)*XYYIMX1T(K,IO, L)
      XZZ = OTA(IO)*EZZT(I,K) + OTB(IO)*XZZIMX1T(K,IO, L)
      TXX = TXX - YL2MX(JM1,IO) *  XXX                                 &
                - YLAXY(JM1,IO) *  XYY                                 &
                - YLAXZ(JM1,IO) *  XZZ
      TYY = TYY - YLAXY(JM1,IO) *  XXX                                 &
                - YL2MY(JM1,IO) *  XYY                                 &
                - YLAYZ(JM1,IO) *  XZZ
      TZZ = TZZ - YLAXZ(JM1,IO) *  XXX                                 &
                - YLAYZ(JM1,IO) *  XYY                                 &
                - YL2MZ(JM1,IO) *  XZZ
      XXXIMX1T(K,IO, L) = XXX
      XYYIMX1T(K,IO, L) = XYY
      XZZIMX1T(K,IO, L) = XZZ
    END DO

    TXXT(I,K, L) = TXXT(I,K, L) + ( DT + TAU_EPS) * TXX
    TYYT(I,K, L) = TYYT(I,K, L) + ( DT + TAU_EPS) * TYY
    TZZT(I,K, L) = TZZT(I,K, L) + ( DT + TAU_EPS) * TZZ
  END IF

END DO


!==== MYT-1 ====
K = MYT-1 + TPML
IF ( ( Y_MIN <= K ) .AND. ( K <= Y_MAX ) ) THEN
  DO  I = MAX(2 - TPML,X_MIN), MIN(MXT-1 + TPML,X_MAX)

    JM1     = JMT (I,K, L)

    TXX             =  L2MX(JM1)         * EXXT(I,K   )                &
                    +  LAXY(JM1)         * EYYT(I,K   )                &
                    +  LAXZ(JM1)         * EZZT(I,K   )

    TYY             =  LAXY(JM1)         * EXXT(I,K   )                &
                    +  L2MY(JM1)         * EYYT(I,K   )                &
                    +  LAYZ(JM1)         * EZZT(I,K   )

    TZZ             =  LAXZ(JM1)         * EXXT(I,K   )                &
                    +  LAYZ(JM1)         * EYYT(I,K   )                &
                    +  L2MZ(JM1)         * EZZT(I,K   )

    DO IO = 1, 4
      XXX = OTA(IO)*EXXT(I,K) + OTB(IO)*XXXKMY1T(I,IO, L)
      XYY = OTA(IO)*EYYT(I,K) + OTB(IO)*XYYKMY1T(I,IO, L)
      XZZ = OTA(IO)*EZZT(I,K) + OTB(IO)*XZZKMY1T(I,IO, L)
      TXX = TXX - YL2MX(JM1,IO) *  XXX                                 &
                - YLAXY(JM1,IO) *  XYY                                 &
                - YLAXZ(JM1,IO) *  XZZ
      TYY = TYY - YLAXY(JM1,IO) *  XXX                                 &
                - YL2MY(JM1,IO) *  XYY                                 &
                - YLAYZ(JM1,IO) *  XZZ
      TZZ = TZZ - YLAXZ(JM1,IO) *  XXX                                 &
                - YLAYZ(JM1,IO) *  XYY                                 &
                - YL2MZ(JM1,IO) *  XZZ
      XXXKMY1T(I,IO, L) = XXX
      XYYKMY1T(I,IO, L) = XYY
      XZZKMY1T(I,IO, L) = XZZ
    END DO

    TXXT(I,K, L) = TXXT(I,K, L) + ( DT + TAU_EPS) * TXX
    TYYT(I,K, L) = TYYT(I,K, L) + ( DT + TAU_EPS) * TYY
    TZZT(I,K, L) = TZZT(I,K, L) + ( DT + TAU_EPS) * TZZ

  END DO
END IF


END SUBROUTINE  STRESS_NN

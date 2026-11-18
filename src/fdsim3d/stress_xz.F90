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

SUBROUTINE  STRESS_XZ (L    , MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,  &
                                          X_MIN, X_MAX, Y_MIN, Y_MAX,  &
                                          TPML, TXZT, EXZT, XXZT, JMT, &
                           XXZK2T, XXZKMY1T, XXZI2T, XXZIMXT )

  USE PRECISION   , ONLY: PP, IK
  USE CONTROL_DATA, ONLY: DT, TAU_EPS
  USE GRID_MEDIUM , ONLY: MUXZ, YMUXZ
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
                                                                  EXZT

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XXZT

  REAL   (PP), DIMENSION (MXL:MXH,       4, MZL:MZH), INTENT(INOUT) :: &
                                                      XXZK2T, XXZKMY1T

  REAL   (PP), DIMENSION (      MYL:MYH, 4, MZL:MZH), INTENT(INOUT) :: &
                                                      XXZI2T, XXZIMXT

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  TXZT

!-----------------------------------------------------------------------

  INTEGER       :: I, K,          JM1,      IO,                        &
                   IN_I_K_, IN_I1K_, IN_I_K1, IN_I1K1

  REAL (KIND=PP):: XZ_I_K_, XZ_I1K_, XZ_I_K1, XZ_I1K1, TXZ, XXZ

!-----------------------------------------------------------------------

!========================= TXZT ===========================
!==== 2 ====
K = 2 - TPML
IF ( ( Y_MIN <= K ) .AND. ( K <= Y_MAX ) ) THEN
  DO  I = MAX(2 - TPML,X_MIN), MIN(MXT + TPML,X_MAX)

    JM1     = JMT (I,K, L )

    TXZ =  MUXZ(JM1)         * EXZT(I,K)

    DO IO = 1, 4
      XXZ = OTA(IO)*EXZT(I,K) + OTB(IO)*XXZK2T(I,IO, L )
      TXZ = TXZ - YMUXZ(JM1,IO) *  XXZ
      XXZK2T(I,IO, L ) = XXZ
    END DO

    TXZT(I,K, L ) = TXZT(I,K, L ) + ( DT + TAU_EPS ) * TXZ

  END DO
END IF


!==== 3 - MYT-2 ====
DO  K = MAX(3-TPML,Y_MIN), MIN(MYT-2 + TPML,Y_MAX)
  I = 2 - TPML
  IF ( ( X_MIN <= I ) .AND. ( I <= X_MAX) ) THEN
    JM1     = JMT (I,K, L )

    TXZ =  MUXZ(JM1)         * EXZT(I,K)

    DO IO = 1, 4
      XXZ = OTA(IO)*EXZT(I,K) + OTB(IO)*XXZI2T(K,IO, L )
      TXZ = TXZ - YMUXZ(JM1,IO) *  XXZ
      XXZI2T(K,IO, L ) = XXZ
    END DO

    TXZT(I,K, L ) = TXZT(I,K, L ) + ( DT + TAU_EPS ) * TXZ
  END IF

  DO  I = MAX(3-TPML,X_MIN), MIN(MXT-1+TPML,X_MAX)

    JM1     = JMT (I,K, L )

    IN_I_K_ = IND (I  ,K  )
    IN_I1K_ = IND (I-1,K  )
    IN_I_K1 = IND (I  ,K-1)
    IN_I1K1 = IND (I-1,K-1)

    XZ_I_K_ =   XXZT(I  ,K  , L )
    XZ_I1K_ = ( XXZT(I-1,K  , L ) + XXZT(I+1,K  , L ) ) * 0.50_PP
    XZ_I_K1 = ( XXZT(I  ,K-1, L ) + XXZT(I  ,K+1, L ) ) * 0.50_PP
    XZ_I1K1 = ( XXZT(I-1,K-1, L ) + XXZT(I+1,K-1, L )                  &
            +   XXZT(I-1,K+1, L ) + XXZT(I+1,K+1, L ) ) * 0.25_PP

    TXZT(I,K, L ) = TXZT(I,K, L ) + ( DT + TAU_EPS ) * (                             &
                       MUXZ(JM1)         * EXZT(I,K   )                &
                    - YMUXZ(JM1,IN_I_K_) *  XZ_I_K_                    &
                    - YMUXZ(JM1,IN_I1K_) *  XZ_I1K_                    &
                    - YMUXZ(JM1,IN_I_K1) *  XZ_I_K1                    &
                    - YMUXZ(JM1,IN_I1K1) *  XZ_I1K1     )
  END DO

  I = MXT + TPML
  IF ( ( X_MIN <= I ) .AND. ( I <= X_MAX) ) THEN
    JM1     = JMT (I,K, L )

    TXZ =  MUXZ(JM1)         * EXZT(I,K)

    DO IO = 1, 4
      XXZ = OTA(IO)*EXZT(I,K) + OTB(IO)*XXZIMXT(K,IO, L )
      TXZ = TXZ - YMUXZ(JM1,IO) *  XXZ
      XXZIMXT(K,IO, L ) = XXZ
    END DO

    TXZT(I,K, L ) = TXZT(I,K, L ) + ( DT + TAU_EPS ) * TXZ
  END IF

END DO


!==== MYT-1 ====
K = MYT-1 + TPML
IF ( ( Y_MIN <= K ) .AND. ( K <= Y_MAX ) ) THEN
  DO  I = MAX(2 - TPML,X_MIN), MIN(MXT + TPML,X_MAX)

    JM1     = JMT (I,K, L )

    TXZ =  MUXZ(JM1)         * EXZT(I,K)

    DO IO = 1, 4
      XXZ = OTA(IO)*EXZT(I,K) + OTB(IO)*XXZKMY1T(I,IO, L )
      TXZ = TXZ - YMUXZ(JM1,IO) *  XXZ
      XXZKMY1T(I,IO, L ) = XXZ
    END DO

    TXZT(I,K, L ) = TXZT(I,K, L ) + ( DT + TAU_EPS ) * TXZ

  END DO
END IF


END SUBROUTINE  STRESS_XZ

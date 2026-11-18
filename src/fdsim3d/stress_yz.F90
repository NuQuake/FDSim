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

SUBROUTINE  STRESS_YZ (L    , MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,  &
                                          X_MIN, X_MAX, Y_MIN, Y_MAX,  &
                                          TPML, TYZT, EYZT, XYZT, JMT, &
                           XYZK2T, XYZKMYT, XYZI2T, XYZIMX1T )

  USE PRECISION   , ONLY: PP, IK
  USE CONTROL_DATA, ONLY: DT, TAU_EPS
  USE GRID_MEDIUM , ONLY: MUYZ, YMUYZ
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
                                                                  EYZT

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XYZT

  REAL   (PP), DIMENSION (MXL:MXH,       4, MZL:MZH), INTENT(INOUT) :: &
                                                      XYZK2T, XYZKMYT

  REAL   (PP), DIMENSION (      MYL:MYH, 4, MZL:MZH), INTENT(INOUT) :: &
                                                      XYZI2T, XYZIMX1T

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  TYZT

!-----------------------------------------------------------------------

  INTEGER       :: I, K,          JM1,      IO,                        &
                   IN_I_K_, IN_I1K_, IN_I_K1, IN_I1K1

  REAL (KIND=PP):: YZ_I_K_, YZ_I1K_, YZ_I_K1, YZ_I1K1, TYZ, XYZ

!-----------------------------------------------------------------------

!========================= TYZT ===========================
!==== 2 ====
K = 2 - TPML
IF ( ( Y_MIN <= K ) .AND. ( K <= Y_MAX ) ) THEN
  DO  I = MAX(2 - TPML,X_MIN), MIN(MXT-1 + TPML,X_MAX)

    JM1     = JMT (I,K, L )

    TYZ =  MUYZ(JM1)         * EYZT(I,K)

    DO IO = 1, 4
      XYZ = OTA(IO)*EYZT(I,K) + OTB(IO)*XYZK2T(I,IO, L )
      TYZ = TYZ - YMUYZ(JM1,IO) *  XYZ
      XYZK2T(I,IO, L ) = XYZ
    END DO

    TYZT(I,K, L ) = TYZT(I,K, L ) + ( DT + TAU_EPS ) * TYZ

  END DO
END IF


!==== 3 - MYT-1 ====
DO  K = MAX(3-TPML,Y_MIN), MIN(MYT-1 + TPML,Y_MAX)
  I = 2 - TPML
  IF ( ( X_MIN <= I ) .AND. ( I <= X_MAX) ) THEN
    JM1     = JMT (I,K, L )

    TYZ =  MUYZ(JM1)         * EYZT(I,K)

    DO IO = 1, 4
      XYZ = OTA(IO)*EYZT(I,K) + OTB(IO)*XYZI2T(K,IO, L )
      TYZ = TYZ - YMUYZ(JM1,IO) *  XYZ
      XYZI2T(K,IO, L ) = XYZ
    END DO

    TYZT(I,K, L ) = TYZT(I,K, L ) + ( DT + TAU_EPS ) * TYZ
  END IF

  DO  I = MAX(3-TPML,X_MIN), MIN(MXT-2+TPML,X_MAX)

    JM1     = JMT (I,K, L )

    IN_I_K_ = IND (I  ,K  )
    IN_I1K_ = IND (I-1,K  )
    IN_I_K1 = IND (I  ,K-1)
    IN_I1K1 = IND (I-1,K-1)

    YZ_I_K_ =   XYZT(I  ,K  , L )
    YZ_I1K_ = ( XYZT(I-1,K  , L ) + XYZT(I+1,K  , L ) ) * 0.50_PP
    YZ_I_K1 = ( XYZT(I  ,K-1, L ) + XYZT(I  ,K+1, L ) ) * 0.50_PP
    YZ_I1K1 = ( XYZT(I-1,K-1, L ) + XYZT(I+1,K-1, L )                  &
            +   XYZT(I-1,K+1, L ) + XYZT(I+1,K+1, L ) ) * 0.25_PP

    TYZT(I,K, L ) = TYZT(I,K, L ) + ( DT + TAU_EPS ) * (                             &
                       MUYZ(JM1)         * EYZT(I,K   )                &
                    - YMUYZ(JM1,IN_I_K_) *  YZ_I_K_                    &
                    - YMUYZ(JM1,IN_I1K_) *  YZ_I1K_                    &
                    - YMUYZ(JM1,IN_I_K1) *  YZ_I_K1                    &
                    - YMUYZ(JM1,IN_I1K1) *  YZ_I1K1      )
  END DO

  I = MXT-1+TPML
  IF ( ( X_MIN <= I ) .AND. ( I <= X_MAX) ) THEN
    JM1     = JMT (I,K, L )

    TYZ =  MUYZ(JM1)         * EYZT(I,K)

    DO IO = 1, 4
      XYZ = OTA(IO)*EYZT(I,K) + OTB(IO)*XYZIMX1T(K,IO, L )
      TYZ = TYZ - YMUYZ(JM1,IO) *  XYZ
      XYZIMX1T(K,IO, L ) = XYZ
    END DO

    TYZT(I,K, L ) = TYZT(I,K, L ) + ( DT + TAU_EPS ) * TYZ
  END IF

END DO


!==== MYT ====
K = MYT+TPML
IF ( ( Y_MIN <= K ) .AND. ( K <= Y_MAX ) ) THEN
  DO  I = MAX(2-TPML,X_MIN), MIN(MXT-1+TPML,X_MAX)

    JM1     = JMT (I,K, L )

    TYZ =  MUYZ(JM1)         * EYZT(I,K)

    DO IO = 1, 4
      XYZ = OTA(IO)*EYZT(I,K) + OTB(IO)*XYZKMYT(I,IO, L )
      TYZ = TYZ - YMUYZ(JM1,IO) *  XYZ
      XYZKMYT(I,IO, L ) = XYZ
    END DO

    TYZT(I,K, L ) = TYZT(I,K, L ) + ( DT + TAU_EPS ) * TYZ

  END DO
END IF


END SUBROUTINE  STRESS_YZ

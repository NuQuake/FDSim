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

SUBROUTINE  INH_V_V_2ND(L, MXL, MXH, PL, MXT, PXH, MX1,                &
                           MYL, MYH,     MYT, PYH, MY1,                &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX, Y_MIN, Y_MAX,                 &
                           JMT, DENVT, VT,                             &
                           TXYT, TYYT, TYZT,                           &
                                           RYYLT, RYYHT,               &
                                           RYXLT, RYXHT,               &
      ODA_XL, ODB_XL, ODC_XL, ODD_XL,  ODA_XH, ODB_XH, ODC_XH, ODD_XH, &
     ODA2_YL,ODB2_YL,ODC2_YL,ODD2_YL, ODA2_YH,ODB2_YH,ODC2_YH,ODD2_YH )

  USE PRECISION   , ONLY:                                              &
                          PP, IK
  USE GRID_MEDIUM , ONLY:                                              &
                          JMNUM

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER    ,                                        INTENT(IN)    :: &
                                        L, MXL, MXH, PL, MXT, MX1, PXH,&
                                           MYL, MYH,     MYT, MY1, PYH,&
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX, Y_MIN, Y_MAX
  INTEGER(IK), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
  REAL   (PP), DIMENSION (1:JMNUM                  ), INTENT(INOUT) :: &
                                                                  DENVT
  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                       TXYT, TYYT, TYZT
  REAL   (PP), DIMENSION (PL :2  , MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  RYXLT
  REAL   (PP), DIMENSION (MX1:PXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  RYXHT
  REAL   (PP), DIMENSION (MXL:MXH, PL :3  , MZL:MZH), INTENT(INOUT) :: &
                                                                  RYYLT
  REAL   (PP), DIMENSION (MXL:MXH, MY1:PYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  RYYHT
  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                     VT
  REAL   (PP), DIMENSION ( 0:TPML, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
  ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,  ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH
  REAL   (PP), DIMENSION(0:TPML+1, MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
 ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL, ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH

  INTEGER   :: I, K, IN2, IN3, JM1, MY2
  REAL (PP) :: DYYYA, DXYXA

!----------------------------------------------------------------------

IN2  = L+0
IN3  = L+1

MY2  = MYT - 2

!_____________________________________________________________ REAR ____
DO K = MAX(3-TPML,Y_MIN), MIN(3,Y_MAX)
  DO  I = MAX(2-TPML,X_MIN), MIN(1,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DYYYA         = TYYT (I  , K, L    ) - TYYT (I, K-1, L    )
    DXYXA         = TXYT (I+1, K, L    ) - TXYT (I, K  , L    )

    VT (I,K,L ) = VT (I,K,L ) + DENVT(JM1)*                            &
                  ( ODA2_YL(3-K,I,L) * ( ODD2_YL(3-K,I,L)*DYYYA        &
                                       + RYYLT(I,K,L)             )    &
                  + ODA_XL (2-I,K,L) * ( ODD_XL (2-I,K,L)*DXYXA        &
                                       + RYXLT(I,K,L)             )    &
                  + TYZT (I  , K, IN3  ) - TYZT (I, K  , IN2  ) )

    RYYLT (I,K,L )= ODB2_YL(3-K,I,L) * RYYLT(I,K,L)                    &
                  + ODC2_YL(3-K,I,L) * DYYYA

    RYXLT (I,K,L )= ODB_XL (2-I,K,L) * RYXLT(I,K,L)                    &
                  + ODC_XL (2-I,K,L) * DXYXA

  END DO

  DO  I = MAX(2,X_MIN), MIN(MX1 ,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DYYYA         = TYYT (I  , K, L    ) - TYYT (I, K-1, L    )
    DXYXA         = TXYT (I+1, K, L    ) - TXYT (I, K  , L    )

    VT (I,K,L ) = VT (I,K,L ) + DENVT(JM1)*                            &
                  ( ODA2_YL(3-K,I,L) * ( ODD2_YL(3-K,I,L)*DYYYA        &
                                       + RYYLT(I,K,L)             )    &
                  +                      DXYXA                         &
                  + TYZT (I  , K, IN3  ) - TYZT (I, K  , IN2  ) )

    RYYLT (I,K,L )= ODB2_YL(3-K,I,L) * RYYLT(I,K,L)                    &
                  + ODC2_YL(3-K,I,L) * DYYYA

  END DO

  DO  I = MAX(MXT,X_MIN), MIN(MX1  + TPML,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DYYYA         = TYYT (I  , K, L    ) - TYYT (I, K-1, L    )
    DXYXA         = TXYT (I+1, K, L    ) - TXYT (I, K  , L    )

    VT (I,K,L ) = VT (I,K,L ) + DENVT(JM1)*                            &
                  ( ODA2_YL(3-K  ,I,L) * ( ODD2_YL(3-K  ,I,L)*DYYYA    &
                                         + RYYLT(I,K,L)             )  &
                  + ODA_XH (I-MX1,K,L) * ( ODD_XH (I-MX1,K,L)*DXYXA    &
                                         + RYXHT(I,K,L)             )  &
                  + TYZT (I  , K, IN3  ) - TYZT (I, K  , IN2  ) )

    RYYLT (I,K,L )= ODB2_YL(3-K  ,I,L) * RYYLT(I,K,L)                  &
                  + ODC2_YL(3-K  ,I,L) * DYYYA

    RYXHT (I,K,L )= ODB_XH (I-MX1,K,L) * RYXHT(I,K,L)                  &
                  + ODC_XH (I-MX1,K,L) * DXYXA

  END DO
END DO
!____________________________________________________________ FRONT ____
DO K = MAX(MY1,Y_MIN), MIN(MY1  + TPML,Y_MAX)
  DO I = MAX(2 - TPML,X_MIN), MIN(1,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DYYYA         = TYYT (I  , K, L    ) - TYYT (I, K-1, L    )
    DXYXA         = TXYT (I+1, K, L    ) - TXYT (I, K  , L    )

    VT (I,K,L ) = VT (I,K,L ) + DENVT(JM1)*                            &
                  ( ODA2_YH(K-MY1,I,L) * ( ODD2_YH(K-MY1,I,L)*DYYYA    &
                                         + RYYHT(I,K,L)             )  &
                  + ODA_XL (2-I  ,K,L) * ( ODD_XL (2-I  ,K,L)*DXYXA    &
                                         + RYXLT(I,K,L)             )  &
                  + TYZT (I  , K, IN3  ) - TYZT (I, K  , IN2  )   )

    RYYHT (I,K,L )= ODB2_YH(K-MY1,I,L) * RYYHT(I,K,L)                  &
                  + ODC2_YH(K-MY1,I,L) * DYYYA

    RYXLT (I,K,L )= ODB_XL (2-I  ,K,L) * RYXLT(I,K,L)                  &
                  + ODC_XL (2-I  ,K,L) * DXYXA

  END DO

  DO  I = MAX(2,X_MIN), MIN(MX1,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DYYYA         = TYYT (I  , K, L    ) - TYYT (I, K-1, L    )
    DXYXA         = TXYT (I+1, K, L    ) - TXYT (I, K  , L    )

    VT (I,K,L ) = VT (I,K,L ) + DENVT(JM1)*                            &
                  ( ODA2_YH(K-MY1,I,L) * ( ODD2_YH(K-MY1,I,L)*DYYYA    &
                                         + RYYHT(I,K,L)             )  &
                  +                        DXYXA                       &
                  + TYZT (I  , K, IN3  ) - TYZT (I, K  , IN2  )   )

    RYYHT (I,K,L )= ODB2_YH(K-MY1,I,L) * RYYHT(I,K,L)                  &
                  + ODC2_YH(K-MY1,I,L) * DYYYA

  END DO

  DO  I = MAX(MXT,X_MIN), MIN(MX1  + TPML,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DYYYA         = TYYT (I  , K, L    ) - TYYT (I, K-1, L    )
    DXYXA         = TXYT (I+1, K, L    ) - TXYT (I, K  , L    )

    VT (I,K,L ) = VT (I,K,L ) + DENVT(JM1)*                            &
                  ( ODA2_YH(K-MY1,I,L) * ( ODD2_YH(K-MY1,I,L)*DYYYA    &
                                         + RYYHT(I,K,L)             )  &
                  + ODA_XH (I-MX1,K,L) * ( ODD_XH (I-MX1,K,L)*DXYXA    &
                                         + RYXHT(I,K,L)             )  &
                  + TYZT (I  , K, IN3  ) - TYZT (I, K  , IN2  )   )

    RYYHT (I,K,L )= ODB2_YH(K-MY1,I,L) * RYYHT(I,K,L)                  &
                  + ODC2_YH(K-MY1,I,L) * DYYYA

    RYXHT (I,K,L )= ODB_XH (I-MX1,K,L) * RYXHT(I,K,L)                  &
                  + ODC_XH (I-MX1,K,L) * DXYXA

  END DO
END DO

!_____________________________________________________________ LEFT ____
DO  K = MAX(4,Y_MIN), MIN(MY2,Y_MAX)
  DO  I = MAX(2 - TPML,X_MIN), MIN(2,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DYYYA         = TYYT (I  , K, L    ) - TYYT (I, K-1, L    )
    DXYXA         = TXYT (I+1, K, L    ) - TXYT (I, K  , L    )

    VT (I,K,L ) = VT (I,K,L ) + DENVT(JM1)*                            &
                  (                        DYYYA                       &
                  + ODA_XL (2-I  ,K,L) * ( ODD_XL (2-I  ,K,L)*DXYXA    &
                                         + RYXLT(I,K,L)             )  &
                  + TYZT (I  , K, IN3  ) - TYZT (I, K  , IN2  )   )

    RYXLT (I,K,L )= ODB_XL (2-I  ,K,L) * RYXLT(I,K,L)                  &
                  + ODC_XL (2-I  ,K,L) * DXYXA

  END DO
END DO

!____________________________________________________________ RIGHT ____
DO  K = MAX(4,Y_MIN), MIN(MY2,Y_MAX)
  DO  I = MAX(MX1,X_MIN), MIN(MX1  + TPML,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DYYYA         = TYYT (I  , K, L    ) - TYYT (I, K-1, L    )
    DXYXA         = TXYT (I+1, K, L    ) - TXYT (I, K  , L    )

    VT (I,K,L ) = VT (I,K,L ) + DENVT(JM1)*                            &
                  (                        DYYYA                       &
                  + ODA_XH (I-MX1,K,L) * ( ODD_XH (I-MX1,K,L)* DXYXA   &
                                         + RYXHT(I,K,L)              ) &
                  + TYZT (I  , K, IN3  ) - TYZT (I, K  , IN2  )   )

    RYXHT (I,K,L )= ODB_XH (I-MX1,K,L) * RYXHT(I,K,L)                  &
                  + ODC_XH (I-MX1,K,L) * DXYXA

  END DO
END DO


END SUBROUTINE  INH_V_V_2ND

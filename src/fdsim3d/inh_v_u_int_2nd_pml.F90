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

SUBROUTINE  INH_V_U_INT_2ND_PML                                        &
                  ( L, LL, MXL, MXH, PL, MXT, PXH, MX1,                &
                           MYL, MYH,     MYT, PYH, MY1,                &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX, Y_MIN, Y_MAX,                 &
                           JMT, DENUT, UT,                             &
                           TXXT, TXYT, TXZT,                           &
                                           RXXLT, RXXHT,               &
                                           RXYLT, RXYHT,               &
                                                  RXZHT,               &
     ODA2_XL,ODB2_XL,ODC2_XL,ODD2_XL, ODA2_XH,ODB2_XH,ODC2_XH,ODD2_XH, &
     ODA_YL ,ODB_YL ,ODC_YL ,ODD_YL , ODA_YH ,ODB_YH ,ODC_YH ,ODD_YH  ,&
                                      ODA_ZH ,ODB_ZH ,ODC_ZH ,ODD_ZH  )

  USE PRECISION   , ONLY:                                              &
                          PP, IK
  USE GRID_MEDIUM , ONLY:                                              &
                          JMNUM

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER    ,                                        INTENT(IN)    :: &
                                    L, LL, MXL, MXH, PL, MXT, MX1, PXH,&
                                           MYL, MYH,     MYT, MY1, PYH,&
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX, Y_MIN, Y_MAX
  INTEGER(IK), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
  REAL   (PP), DIMENSION (1:JMNUM                  ), INTENT(INOUT) :: &
                                                                  DENUT
  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                       TXXT, TXYT, TXZT
  REAL   (PP), DIMENSION (PL :3  , MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  RXXLT
  REAL   (PP), DIMENSION (MX1:PXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  RXXHT
  REAL   (PP), DIMENSION (MXL:MXH, PL :2  , MZL:MZH), INTENT(INOUT) :: &
                                                                  RXYLT
  REAL   (PP), DIMENSION (MXL:MXH, MY1:PYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  RXYHT
  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH,  0:TPML), INTENT(INOUT) :: &
                                                                  RXZHT
  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                     UT
  REAL   (PP), DIMENSION ( 0:TPML, MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
  ODA_YL,  ODB_YL,  ODC_YL,  ODD_YL,  ODA_YH,  ODB_YH,  ODC_YH,  ODD_YH
  REAL   (PP), DIMENSION(0:TPML+1, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
  ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL, ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH
  REAL   (PP), DIMENSION ( 0:TPML, MXL:MXH, MYL:MYH), INTENT(INOUT) :: &
                                      ODA_ZH,  ODB_ZH,  ODC_ZH,  ODD_ZH

  INTEGER   :: I, K, IN2, IN3, JM1, MY2
  REAL (PP) :: DXXXA, DXYYA, DXZZA

!----------------------------------------------------------------------

IN2  = L+0
IN3  = L+1

MY2  = MYT - 2

!_____________________________________________________________ REAR ____
DO K = MAX(2 - TPML,Y_MIN), MIN(2,Y_MAX)
  DO  I = MAX(3 - TPML,X_MIN), MIN(2,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DXXXA         = TXXT (I, K  ,L    ) - TXXT (I-1, K, L    )
    DXYYA         = TXYT (I, K+1,L    ) - TXYT (I  , K, L    )
    DXZZA         = TXZT (I, K  ,IN3  ) - TXZT (I  , K, IN2  )

    UT (I,K,L ) = UT (I,K,L ) + DENUT(JM1)*                            &
                  ( ODA2_XL(3-I,K,L) * ( ODD2_XL(3-I,K,L)*DXXXA        &
                                       + RXXLT(I,K,L)             )    &
                  + ODA_YL (2-K,I,L) * ( ODD_YL (2-K,I,L)*DXYYA        &
                                       + RXYLT(I,K,L)             )    &
                  + ODA_ZH (LL ,I,K) * ( ODD_ZH (LL ,I,K)*DXZZA        &
                                       + RXZHT(I,K,LL)            ) )

    RXXLT (I,K,L )= ODB2_XL(3-I,K,L) * RXXLT(I,K,L )                   &
                  + ODC2_XL(3-I,K,L) * DXXXA

    RXYLT (I,K,L )= ODB_YL (2-K,I,L) * RXYLT(I,K,L )                   &
                  + ODC_YL (2-K,I,L) * DXYYA

    RXZHT (I,K,LL)= ODB_ZH (LL ,I,K) * RXZHT(I,K,LL)                   &
                  + ODC_ZH (LL ,I,K) * DXZZA

  END DO

  DO  I = MAX(3,X_MIN), MIN(MX1,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DXXXA         = TXXT (I, K  ,L    ) - TXXT (I-1, K, L    )
    DXYYA         = TXYT (I, K+1,L    ) - TXYT (I  , K, L    )
    DXZZA         = TXZT (I, K  ,IN3  ) - TXZT (I  , K, IN2  )

    UT (I,K,L ) = UT (I,K,L ) + DENUT(JM1)*                            &
                  ( DXXXA                                              &
                  + ODA_YL (2-K,I,L) * ( ODD_YL (2-K,I,L)*DXYYA        &
                                       + RXYLT(I,K,L)             )    &
                  + ODA_ZH (LL ,I,K) * ( ODD_ZH (LL ,I,K)*DXZZA        &
                                       + RXZHT(I,K,LL)            ) )

    RXYLT (I,K,L )= ODB_YL (2-K,I,L) * RXYLT(I,K,L )                   &
                  + ODC_YL (2-K,I,L) * DXYYA

    RXZHT (I,K,LL)= ODB_ZH (LL ,I,K) * RXZHT(I,K,LL)                   &
                  + ODC_ZH (LL ,I,K) * DXZZA

  END DO

  DO  I = MAX(MXT,X_MIN), MIN(MX1  + TPML,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DXXXA         = TXXT (I, K  ,L    ) - TXXT (I-1, K, L    )
    DXYYA         = TXYT (I, K+1,L    ) - TXYT (I  , K, L    )
    DXZZA         = TXZT (I, K  ,IN3  ) - TXZT (I  , K, IN2  )

    UT (I,K,L ) = UT (I,K,L ) + DENUT(JM1)*                            &
                  ( ODA2_XH(I-MX1,K,L) * ( ODD2_XH(I-MX1,K,L)*DXXXA    &
                                         + RXXHT(I,K,L)             )  &
                  + ODA_YL (2-K  ,I,L) * ( ODD_YL (2-K  ,I,L)*DXYYA    &
                                         + RXYLT(I,K,L)             )  &
                  + ODA_ZH (LL   ,I,K) * ( ODD_ZH (LL   ,I,K)*DXZZA    &
                                         + RXZHT(I,K,LL)            ) )

    RXXHT (I,K,L )= ODB2_XH(I-MX1,K,L) * RXXHT(I,K,L)                  &
                  + ODC2_XH(I-MX1,K,L) * DXXXA

    RXYLT (I,K,L )= ODB_YL (2-K  ,I,L) * RXYLT(I,K,L)                  &
                  + ODC_YL (2-K  ,I,L) * DXYYA

    RXZHT (I,K,LL)= ODB_ZH (LL ,I,K) * RXZHT(I,K,LL)                   &
                  + ODC_ZH (LL ,I,K) * DXZZA

  END DO
END DO
!____________________________________________________________ FRONT ____
DO K = MAX(MY1 ,Y_MIN), MIN(MY1  + TPML,Y_MAX)
  DO  I = MAX(3 - TPML,X_MIN), MIN(2,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DXXXA         = TXXT (I, K  ,L    ) - TXXT (I-1, K, L    )
    DXYYA         = TXYT (I, K+1,L    ) - TXYT (I  , K, L    )
    DXZZA         = TXZT (I, K  ,IN3  ) - TXZT (I  , K, IN2  )

    UT (I,K,L ) = UT (I,K,L ) + DENUT(JM1)*                            &
                  ( ODA2_XL(3-I  ,K,L) * ( ODD2_XL(3-I  ,K,L)*DXXXA    &
                                         + RXXLT(I,K,L)             )  &
                  + ODA_YH (K-MY1,I,L) * ( ODD_YH (K-MY1,I,L)*DXYYA    &
                                         + RXYHT(I,K,L)             )  &
                  + ODA_ZH (LL   ,I,K) * ( ODD_ZH (LL   ,I,K)*DXZZA    &
                                         + RXZHT(I,K,LL)            ) )

    RXXLT (I,K,L )= ODB2_XL(3-I  ,K,L) * RXXLT(I,K,L )                 &
                  + ODC2_XL(3-I  ,K,L) * DXXXA

    RXYHT (I,K,L )= ODB_YH (K-MY1,I,L) * RXYHT(I,K,L )                 &
                  + ODC_YH (K-MY1,I,L) * DXYYA

    RXZHT (I,K,LL)= ODB_ZH (LL   ,I,K) * RXZHT(I,K,LL)                 &
                  + ODC_ZH (LL   ,I,K) * DXZZA

  END DO

  DO  I = MAX(3,X_MIN), MIN(MX1,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DXXXA         = TXXT (I, K  ,L    ) - TXXT (I-1, K, L    )
    DXYYA         = TXYT (I, K+1,L    ) - TXYT (I  , K, L    )
    DXZZA         = TXZT (I, K  ,IN3  ) - TXZT (I  , K, IN2  )

    UT (I,K,L ) = UT (I,K,L ) + DENUT(JM1)*                            &
                  ( DXXXA                                              &
                  + ODA_YH (K-MY1,I,L) * ( ODD_YH (K-MY1,I,L)*DXYYA    &
                                         + RXYHT(I,K,L)             )  &
                  + ODA_ZH (LL   ,I,K) * ( ODD_ZH (LL   ,I,K)*DXZZA    &
                                         + RXZHT(I,K,LL)            ) )

    RXYHT (I,K,L )= ODB_YH (K-MY1,I,L) * RXYHT(I,K,L )                 &
                  + ODC_YH (K-MY1,I,L) * DXYYA

    RXZHT (I,K,LL)= ODB_ZH (LL   ,I,K) * RXZHT(I,K,LL)                 &
                  + ODC_ZH (LL   ,I,K) * DXZZA

  END DO

  DO  I = MAX(MXT,X_MIN), MIN(MX1  + TPML,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DXXXA         = TXXT (I, K  ,L    ) - TXXT (I-1, K, L    )
    DXYYA         = TXYT (I, K+1,L    ) - TXYT (I  , K, L    )
    DXZZA         = TXZT (I, K  ,IN3  ) - TXZT (I  , K, IN2  )

    UT (I,K,L ) = UT (I,K,L ) + DENUT(JM1)*                            &
                  ( ODA2_XH(I-MX1,K,L) * ( ODD2_XH(I-MX1,K,L)*DXXXA    &
                                         + RXXHT(I,K,L)             )  &
                  + ODA_YH (K-MY1,I,L) * ( ODD_YH (K-MY1,I,L)*DXYYA    &
                                         + RXYHT(I,K,L)             )  &
                  + ODA_ZH (LL   ,I,K) * ( ODD_ZH (LL   ,I,K)*DXZZA    &
                                         + RXZHT(I,K,LL) ) )

    RXXHT (I,K,L )= ODB2_XH(I-MX1,K,L) * RXXHT(I,K,L )                 &
                  + ODC2_XH(I-MX1,K,L) * DXXXA

    RXYHT (I,K,L )= ODB_YH (K-MY1,I,L) * RXYHT(I,K,L )                 &
                  + ODC_YH (K-MY1,I,L) * DXYYA

    RXZHT (I,K,LL)= ODB_ZH (LL   ,I,K) * RXZHT(I,K,LL)                 &
                  + ODC_ZH (LL   ,I,K) * DXZZA

  END DO
END DO

!_____________________________________________________________ LEFT ____
DO  K = MAX(3,Y_MIN), MIN(MY2,Y_MAX)
  DO  I = MAX(3 - TPML,X_MIN), MIN(3,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DXXXA         = TXXT (I, K  ,L    ) - TXXT (I-1, K, L    )
    DXYYA         = TXYT (I, K+1,L    ) - TXYT (I  , K, L    )
    DXZZA         = TXZT (I, K  ,IN3  ) - TXZT (I  , K, IN2  )

    UT (I,K,L ) = UT (I,K,L ) + DENUT(JM1)*                            &
                  ( ODA2_XL(3-I  ,K,L) * ( ODD2_XL(3-I  ,K,L)*DXXXA    &
                                         + RXXLT(I,K,L)             )  &
                  + DXYYA                                              &
                  + ODA_ZH (LL   ,I,K) * ( ODD_ZH (LL   ,I,K)*DXZZA    &
                                         + RXZHT(I,K,LL)            ) )

    RXXLT (I,K,L )= ODB2_XL(3-I  ,K,L) * RXXLT(I,K,L )                 &
                  + ODC2_XL(3-I  ,K,L) * DXXXA

    RXZHT (I,K,LL)= ODB_ZH (LL   ,I,K) * RXZHT(I,K,LL)                 &
                  + ODC_ZH (LL   ,I,K) * DXZZA

  END DO

  DO  I = MAX(4,X_MIN), MIN(MX1 -1,X_MAX)
    JM1 = JMT (I  ,K  ,L   )

    DXXXA         = TXXT (I, K  ,L    ) - TXXT (I-1, K, L    )
    DXYYA         = TXYT (I, K+1,L    ) - TXYT (I  , K, L    )
    DXZZA         = TXZT (I, K  ,IN3  ) - TXZT (I  , K, IN2  )

    UT (I,K,L ) = UT (I,K,L ) + DENUT(JM1)*                            &
                  ( DXXXA                                              &
                  + DXYYA                                              &
                  + ODA_ZH (LL   ,I,K) * ( ODD_ZH (LL   ,I,K)*DXZZA    &
                                         + RXZHT(I,K,LL)            ) )

    RXZHT (I,K,LL)= ODB_ZH (LL   ,I,K) * RXZHT(I,K,LL)                 &
                  + ODC_ZH (LL   ,I,K) * DXZZA

  END DO

!____________________________________________________________ RIGHT ____

  DO  I = MAX(MX1 ,X_MIN), MIN(MX1  + TPML,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DXXXA         = TXXT (I, K  ,L    ) - TXXT (I-1, K, L    )
    DXYYA         = TXYT (I, K+1,L    ) - TXYT (I  , K, L    )
    DXZZA         = TXZT (I, K  ,IN3  ) - TXZT (I  , K, IN2  )

    UT (I,K,L ) = UT (I,K,L ) + DENUT(JM1)*                            &
                  ( ODA2_XH(I-MX1,K,L) * ( ODD2_XH(I-MX1,K,L)* DXXXA   &
                                         + RXXHT(I,K,L)              ) &
                  + DXYYA                                              &
                  + ODA_ZH (LL   ,I,K) * ( ODD_ZH (LL   ,I,K)* DXZZA   &
                                         + RXZHT(I,K,LL)             ) )

    RXXHT (I,K,L )= ODB2_XH(I-MX1,K,L) * RXXHT(I,K,L )                 &
                  + ODC2_XH(I-MX1,K,L) * DXXXA

    RXZHT (I,K,LL)= ODB_ZH (LL   ,I,K) * RXZHT(I,K,LL)                 &
                  + ODC_ZH (LL   ,I,K) * DXZZA

  END DO
END DO


END SUBROUTINE  INH_V_U_INT_2ND_PML

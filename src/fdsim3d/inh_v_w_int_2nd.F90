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

SUBROUTINE  INH_V_W_INT_2ND                                            &
                  ( L,     MXL, MXH, PL, MXT, PXH, MX1,                &
                           MYL, MYH,     MYT, PYH, MY1,                &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX, Y_MIN, Y_MAX,                 &
                           JMT, DENWT, WT,                             &
                           TZZT, TXZT, TYZT,                           &
                                           RZXLT, RZXHT,               &
                                           RZYLT, RZYHT,               &
ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,  ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH, &
ODA_YL,  ODB_YL,  ODC_YL,  ODD_YL,  ODA_YH,  ODB_YH,  ODC_YH,  ODD_YH )

  USE PRECISION   , ONLY:                                              &
                          PP, IK
  USE GRID_MEDIUM , ONLY:                                              &
                          JMNUM

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER    ,                                        INTENT(IN)    :: &
                                             L, MXL, MXH, PL, MXT, PXH,&
                                                MYL, MYH,     MYT, PYH,&
                                                MZL, MZH,  TPML,       &
                                                MX1, MY1,              &
                                            X_MIN, X_MAX, Y_MIN, Y_MAX
  INTEGER(IK), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
  REAL   (PP), DIMENSION (1:JMNUM                  ), INTENT(INOUT) :: &
                                                                  DENWT
  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                       TXZT, TYZT, TZZT
  REAL   (PP), DIMENSION (PL :2  , MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  RZXLT
  REAL   (PP), DIMENSION (MX1:PXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  RZXHT
  REAL   (PP), DIMENSION (MXL:MXH, PL :2  , MZL:MZH), INTENT(INOUT) :: &
                                                                  RZYLT
  REAL   (PP), DIMENSION (MXL:MXH, MY1:PYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  RZYHT
  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                     WT
  REAL   (PP), DIMENSION ( 0:TPML, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
  ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,  ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH
  REAL   (PP), DIMENSION ( 0:TPML, MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
  ODA_YL,  ODB_YL,  ODC_YL,  ODD_YL,  ODA_YH,  ODB_YH,  ODC_YH,  ODD_YH

  INTEGER   :: I, K, IN2, IN3, JM1, MY2
  REAL (PP) :: DXZXA, DYZYA

!-----------------------------------------------------------------------

IN2  = L-1
IN3  = L

MY2  = MYT - 2

!_____________________________________________________________ REAR ____
DO K = MAX(2 - TPML,Y_MIN), MIN(2,Y_MAX)
  DO  I = MAX(2 - TPML,X_MIN), MIN(1,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DXZXA         = TXZT (I+1,K  ,L    ) - TXZT (I  ,K  ,L    )
    DYZYA         = TYZT (I  ,K+1,L    ) - TYZT (I  ,K  ,L    )

    WT (I,K,L ) = WT (I,K,L ) + DENWT(JM1)*                            &
                 (     TZZT (I  ,K  ,IN3  ) - TZZT (I  ,K  ,IN2  )     &
                 + ODA_XL (2-I,K,L) * ( ODD_XL (2-I,K,L)*DXZXA         &
                                      + RZXLT(I,K,L)             )     &
                 + ODA_YL (2-K,I,L) * ( ODD_YL (2-K,I,L)*DYZYA         &
                                      + RZYLT(I,K,L)             ) )

    RZXLT (I,K,L )= ODB_XL (2-I,K,L) * RZXLT(I,K,L)                    &
                  + ODC_XL (2-I,K,L) * DXZXA

    RZYLT (I,K,L )= ODB_YL (2-K,I,L) * RZYLT(I,K,L)                    &
                  + ODC_YL (2-K,I,L) * DYZYA

  END DO

  DO  I = MAX(2,X_MIN), MIN(MX1 ,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DXZXA         = TXZT (I+1,K  ,L    ) - TXZT (I  ,K  ,L    )
    DYZYA         = TYZT (I  ,K+1,L    ) - TYZT (I  ,K  ,L    )

    WT (I,K,L ) = WT (I,K,L ) + DENWT(JM1)*                            &
                 (     TZZT (I  ,K  ,IN3  ) - TZZT (I  ,K  ,IN2  )     &
                 +                      DXZXA                          &
                 + ODA_YL (2-K,I,L) * ( ODD_YL (2-K,I,L)*DYZYA         &
                                      + RZYLT(I,K,L)           ) )

    RZYLT (I,K,L )= ODB_YL (2-K,I,L) * RZYLT(I,K,L)                    &
                  + ODC_YL (2-K,I,L) * DYZYA

  END DO

  DO  I = MAX(MXT,X_MIN), MIN(MX1  + TPML,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DXZXA         = TXZT (I+1,K  ,L    ) - TXZT (I  ,K  ,L    )
    DYZYA         = TYZT (I  ,K+1,L    ) - TYZT (I  ,K  ,L    )

    WT (I,K,L ) = WT (I,K,L ) + DENWT(JM1)*                            &
                 (     TZZT (I  ,K  ,IN3  ) - TZZT (I  ,K  ,IN2  )     &
                 + ODA_XH (I-MX1,K,L) * ( ODD_XH (I-MX1,K,L)*DXZXA     &
                                        + RZXHT(I,K,L)             )   &
                 + ODA_YL (2-K  ,I,L) * ( ODD_YL (2-K  ,I,L)*DYZYA     &
                                        + RZYLT(I,K,L)             ) )

    RZXHT (I,K,L )= ODB_XH (I-MX1,K,L) * RZXHT(I,K,L)                  &
                  + ODC_XH (I-MX1,K,L) * DXZXA

    RZYLT (I,K,L )= ODB_YL (2-K,I,L) * RZYLT(I,K,L)                    &
                  + ODC_YL (2-K,I,L) * DYZYA

  END DO
END DO
!____________________________________________________________ FRONT ____
DO K = MAX(MY1 ,Y_MIN), MIN(MY1  + TPML,Y_MAX)
  DO  I = MAX(2 - TPML,X_MIN), MIN(1,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DXZXA         = TXZT (I+1,K  ,L    ) - TXZT (I  ,K  ,L    )
    DYZYA         = TYZT (I  ,K+1,L    ) - TYZT (I  ,K  ,L    )

    WT (I,K,L ) = WT (I,K,L ) + DENWT(JM1)*                            &
                 (     TZZT (I  ,K  ,IN3  ) - TZZT (I  ,K  ,IN2  )     &
                 + ODA_XL (2-I  ,K,L) * ( ODD_XL (2-I  ,K,L)*DXZXA     &
                                        + RZXLT(I,K,L)             )   &
                 + ODA_YH (K-MY1,I,L) * ( ODD_YH (K-MY1,I,L)*DYZYA     &
                                        + RZYHT(I,K,L)             ) )

    RZXLT (I,K,L )= ODB_XL (2-I,K,L) * RZXLT(I,K,L)                    &
                  + ODC_XL (2-I,K,L) * DXZXA

    RZYHT (I,K,L )= ODB_YH (K-MY1,I,L) * RZYHT(I,K,L)                  &
                  + ODC_YH (K-MY1,I,L) * DYZYA

  END DO

  DO  I = MAX(2,X_MIN), MIN(MX1 ,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DXZXA         = TXZT (I+1,K  ,L    ) - TXZT (I  ,K  ,L    )
    DYZYA         = TYZT (I  ,K+1,L    ) - TYZT (I  ,K  ,L    )

    WT (I,K,L ) = WT (I,K,L ) + DENWT(JM1)*                            &
                 (     TZZT (I  ,K  ,IN3  ) - TZZT (I  ,K  ,IN2  )     &
                     + DXZXA                                           &
                 + ODA_YH (K-MY1,I,L) * ( ODD_YH (K-MY1,I,L)*DYZYA     &
                                        + RZYHT(I,K,L)             ) )

    RZYHT (I,K,L )= ODB_YH (K-MY1,I,L) * RZYHT(I,K,L)                  &
                  + ODC_YH (K-MY1,I,L) * DYZYA

  END DO

  DO  I = MAX(MXT,X_MIN), MIN(MX1  + TPML,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DXZXA         = TXZT (I+1,K  ,L    ) - TXZT (I  ,K  ,L    )
    DYZYA         = TYZT (I  ,K+1,L    ) - TYZT (I  ,K  ,L    )

    WT (I,K,L ) = WT (I,K,L ) + DENWT(JM1)*                            &
                 (     TZZT (I  ,K  ,IN3  ) - TZZT (I  ,K  ,IN2  )     &
                 + ODA_XH (I-MX1,K,L) * ( ODD_XH (I-MX1,K,L)*DXZXA     &
                                        + RZXHT(I,K,L)             )   &
                 + ODA_YH (K-MY1,I,L) * ( ODD_YH (K-MY1,I,L)*DYZYA     &
                                        + RZYHT(I,K,L)             ) )

    RZXHT (I,K,L )= ODB_XH (I-MX1,K,L) * RZXHT(I,K,L)                  &
                  + ODC_XH (I-MX1,K,L) * DXZXA

    RZYHT (I,K,L )= ODB_YH (K-MY1,I,L) * RZYHT(I,K,L)                  &
                  + ODC_YH (K-MY1,I,L) * DYZYA

  END DO
END DO

!_____________________________________________________________ LEFT ____
DO  K = MAX(3,Y_MIN), MIN(MY2,Y_MAX)
  DO  I = MAX(2 - TPML,X_MIN), MIN(2,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DXZXA         = TXZT (I+1,K  ,L    ) - TXZT (I  ,K  ,L    )
    DYZYA         = TYZT (I  ,K+1,L    ) - TYZT (I  ,K  ,L    )

    WT (I,K,L ) = WT (I,K,L ) + DENWT(JM1)*                            &
                 (     TZZT (I  ,K  ,IN3  ) - TZZT (I  ,K  ,IN2  )     &
                 + ODA_XL (2-I  ,K,L) * ( ODD_XL (2-I  ,K,L)*DXZXA     &
                                        + RZXLT(I,K,L)             )   &
                     + DYZYA )

    RZXLT (I,K,L )= ODB_XL (2-I,K,L) * RZXLT(I,K,L)                    &
                  + ODC_XL (2-I,K,L) * DXZXA

  END DO

  DO  I = MAX(3,X_MIN), MIN(MX1 -1,X_MAX)
    JM1 = JMT (I  ,K  ,L   )

    DXZXA         = TXZT (I+1,K  ,L    ) - TXZT (I  ,K  ,L    )
    DYZYA         = TYZT (I  ,K+1,L    ) - TYZT (I  ,K  ,L    )

    WT (I,K,L ) = WT (I,K,L ) + DENWT(JM1)*                            &
                 (     TZZT (I  ,K  ,IN3  ) - TZZT (I  ,K  ,IN2  )     &
                     + DXZXA                                           &
                     + DYZYA )
  END DO

!____________________________________________________________ RIGHT ____

  DO  I = MAX(MX1 ,X_MIN), MIN(MX1  + TPML,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    DXZXA         = TXZT (I+1,K  ,L    ) - TXZT (I  ,K  ,L    )
    DYZYA         = TYZT (I  ,K+1,L    ) - TYZT (I  ,K  ,L    )

    WT (I,K,L ) = WT (I,K,L ) + DENWT(JM1)*                            &
                 (     TZZT (I  ,K  ,IN3  ) - TZZT (I  ,K  ,IN2  )     &
                 + ODA_XH (I-MX1,K,L) * ( ODD_XH (I-MX1,K,L)*DXZXA     &
                                        + RZXHT(I,K,L)             )   &
                     + DYZYA )

    RZXHT (I,K,L )= ODB_XH (I-MX1,K,L) * RZXHT(I,K,L)                  &
                  + ODC_XH (I-MX1,K,L) * DXZXA

  END DO
END DO


END SUBROUTINE  INH_V_W_INT_2ND

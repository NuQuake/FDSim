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

SUBROUTINE  INH_V_UV_INT(L, MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,    &
                            X_MIN, X_MAX, Y_MIN, Y_MAX,                &
                            JMT, DENUT, DENVT, UT, VT,                 &
                            TXXT, TXYT, TXZT, TYYT,       TYZT )

  USE PRECISION   , ONLY:                                              &
                          PP, IK
  USE GRID_MEDIUM , ONLY:                                              &
                          JMNUM
  USE AUXIL       , ONLY:                                              &
                          A, B

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,                                            INTENT(IN)    :: &
                             L, MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,&
                             X_MIN, X_MAX, Y_MIN, Y_MAX
  INTEGER(IK), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
  REAL   (PP), DIMENSION (1:JMNUM                  ), INTENT(INOUT) :: &
                                                           DENUT, DENVT
  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                           TXXT, TXYT, TYYT, TXZT, TYZT
  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                               UT  , VT

  INTEGER              :: I, K, IN1, IN2, IN3, IN4, JM1

!-----------------------------------------------------------------------

IN1  = L-1
IN2  = L+0
IN3  = L+1
IN4  = L+2

DO  K = MAX(3,Y_MIN), MIN(MYT-2,Y_MAX)
  DO  I = MAX(4,X_MIN), MIN(MXT-2,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    UT (I,K,L ) = UT (I,K,L ) + DENUT(JM1)*                            &
                 ( A*( TXXT (I+1,K  ,L    ) - TXXT (I-2,K  ,L    )     &
                     + TXYT (I  ,K+2,L    ) - TXYT (I  ,K-1,L    )     &
                     + TXZT (I  ,K  ,IN4  ) - TXZT (I  ,K  ,IN1  ) )   &
                 + B*( TXXT (I  ,K  ,L    ) - TXXT (I-1,K  ,L    )     &
                     + TXYT (I  ,K+1,L    ) - TXYT (I  ,K  ,L    )     &
                     + TXZT (I  ,K  ,IN3  ) - TXZT (I  ,K  ,IN2  ) ) )

  END DO
END DO

DO  K = MAX(4,Y_MIN), MIN(MYT-2,Y_MAX)
  DO  I = MAX(3,X_MIN), MIN(MXT-2,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

    VT (I,K,L ) = VT (I,K,L ) + DENVT(JM1)*                            &
                 ( A*( TYYT (I  ,K+1,L    ) - TYYT (I  ,K-2, L    )    &
                     + TXYT (I+2,K  ,L    ) - TXYT (I-1,K  , L    )    &
                     + TYZT (I  ,K  ,IN4  ) - TYZT (I  ,K  , IN1  ) )  &
                 + B*( TYYT (I  ,K  ,L    ) - TYYT (I  ,K-1, L    )    &
                     + TXYT (I+1,K  ,L    ) - TXYT (I  ,K  , L    )    &
                     + TYZT (I  ,K  ,IN3  ) - TYZT (I  ,K  , IN2  ) ) )

  END DO
END DO


END SUBROUTINE  INH_V_UV_INT

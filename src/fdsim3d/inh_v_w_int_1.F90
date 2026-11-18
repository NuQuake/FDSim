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

SUBROUTINE  INH_V_W_INT_1 ( MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,    &
                            X_MIN, X_MAX, Y_MIN, Y_MAX,                &
                            JMT, DENWT, WT, TZZT, TXZT, TYZT )

  USE PRECISION   , ONLY:                                              &
                          PP, IK
  USE GRID_MEDIUM , ONLY:                                              &
                          JMNUM
  USE AUXIL       , ONLY:                                              &
                          A, B

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER    ,                                        INTENT(IN)    :: &
                               MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH, &
                               X_MIN, X_MAX, Y_MIN, Y_MAX
  INTEGER(IK), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
  REAL   (PP), DIMENSION (1:JMNUM                  ), INTENT(INOUT) :: &
                                                                  DENWT
  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                       TXZT, TYZT, TZZT
  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                     WT

  INTEGER   :: I, K, L, IN1, IN2, IN3, IN4, JM1

!-----------------------------------------------------------------------

L    = 1

IN2  = L-1
IN3  = L
IN4  = L+1
IN1  = L+2


DO  K = MAX(3,Y_MIN), MIN(MYT-2,Y_MAX)
  DO  I = MAX(3,X_MIN), MIN(MXT-2,X_MAX)

    JM1 = JMT (I  ,K  ,L   )

#if(SEDOND)
    WT (I,K,L ) = WT (I,K,L ) + DENWT(JM1)*                            &
                 (     TXZT (I+1,K  ,L    ) - TXZT (I  ,K  ,L    )     &
                 +     TYZT (I  ,K+1,L    ) - TYZT (I  ,K  ,L    )     &
                 +     TZZT (I  ,K   ,IN3 ) - TZZT (I  ,K  ,IN2  ) )
#else
    WT (I,K,L ) = WT (I,K,L ) + DENWT(JM1)*                            &
                 ( A*( TXZT (I+2,K  ,L    ) - TXZT (I-1,K  ,L    )     &
                     + TYZT (I  ,K+2,L    ) - TYZT (I  ,K-1,L    ) )   &
                 + B*( TXZT (I+1,K  ,L    ) - TXZT (I  ,K  ,L    )     &
                     + TYZT (I  ,K+1,L    ) - TYZT (I  ,K  ,L    ) )   &
                 - 31._PP/ 24._PP * TZZT (I  ,K  ,IN2  )               &
                 + 29._PP/ 24._PP * TZZT (I  ,K  ,IN3  )               &
                 -  3._PP/ 40._PP * TZZT (I  ,K  ,IN4  )               &
                 +  1._PP/168._PP * TZZT (I  ,K  ,IN1  ) )
#endif
  END DO
END DO


END SUBROUTINE  INH_V_W_INT_1

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

SUBROUTINE  INH_V_W_2ND_0( MXL, MXH,     MXT,                          &
                           MYL, MYH,     MYT,                          &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX, Y_MIN, Y_MAX,                 &
                           JMT, DENWT, WT, TZZT )

  USE PRECISION   , ONLY:                                              &
                          PP, IK
  USE GRID_MEDIUM , ONLY:                                              &
                          JMNUM

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER    ,                                        INTENT(IN)    :: &
                                    MXL, MXH, MYL, MYH, MXT, MYT, TPML,&
                                    MZL, MZH,                          &
                                    X_MIN, X_MAX, Y_MIN, Y_MAX
  INTEGER(IK), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
  REAL   (PP), DIMENSION (1:JMNUM                  ), INTENT(INOUT) :: &
                                                                  DENWT
  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                   TZZT
  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                     WT

  INTEGER              :: I, K, L, JM1, MXT1, MYT1, MYT2

!-----------------------------------------------------------------------
L    = 0

MXT1 = MXT - 1
MYT1 = MYT - 1
MYT2 = MYT - 2

!_____________________________________________________________ REAR ____
DO K = MAX(2 - TPML,Y_MIN), MIN(2,Y_MAX)
  DO  I = MAX(2 - TPML,X_MIN), MIN(MXT1 + TPML,X_MAX)

    JM1 = JMT (I  ,K  , L  )

    WT (I,K, L) = WT (I,K, L) + DENWT(JM1)*                            &
                 ( 3._PP         * TZZT (I  ,K  ,  0  )                &
                 - 1._PP/3._PP   * TZZT (I  ,K  ,  1  ) )
!                 ( 35._PP/ 8._PP * TZZT (I  ,K  ,  0  )                &
!                 - 35._PP/24._PP * TZZT (I  ,K  ,  1  )                &
!                 + 21._PP/40._PP * TZZT (I  ,K  ,  2  )                &
!                 -  5._PP/56._PP * TZZT (I  ,K  ,  3  ) )

  END DO
END DO

!____________________________________________________________ FRONT ____
DO K = MAX(MYT1,Y_MIN), MIN(MYT1 + TPML,Y_MAX)
  DO I = MAX(2 - TPML,X_MIN), MIN(MXT1 + TPML,X_MAX)

    JM1 = JMT (I  ,K  , L  )

    WT (I,K, L) = WT (I,K, L) + DENWT(JM1)*                            &
                 ( 3._PP         * TZZT (I  ,K  ,  0  )                &
                 - 1._PP/3._PP   * TZZT (I  ,K  ,  1  ) )
!                 ( 35._PP/ 8._PP * TZZT (I  ,K  ,  0  )                &
!                 - 35._PP/24._PP * TZZT (I  ,K  ,  1  )                &
!                 + 21._PP/40._PP * TZZT (I  ,K  ,  2  )                &
!                 -  5._PP/56._PP * TZZT (I  ,K  ,  3  ) )

  END DO
END DO

!_____________________________________________________________ LEFT ____
DO I = MAX(2 - TPML,X_MIN), MIN(2,X_MAX)
  DO K = MAX(3,Y_MIN), MIN(MYT2,Y_MAX)

    JM1 = JMT (I  ,K  , L  )

    WT (I,K, L) = WT (I,K, L) + DENWT(JM1)*                            &
                 ( 3._PP         * TZZT (I  ,K  ,  0  )                &
                 - 1._PP/3._PP   * TZZT (I  ,K  ,  1  ) )
!                 ( 35._PP/ 8._PP * TZZT (I  ,K  ,  0  )                &
!                 - 35._PP/24._PP * TZZT (I  ,K  ,  1  )                &
!                 + 21._PP/40._PP * TZZT (I  ,K  ,  2  )                &
!                 -  5._PP/56._PP * TZZT (I  ,K  ,  3  ) )

  END DO
END DO

!____________________________________________________________ RIGHT ____
DO I = MAX(MXT1,X_MIN), MIN(MXT1 + TPML,X_MAX)
  DO K = MAX(3,Y_MIN), MIN(MYT2,Y_MAX)

    JM1 = JMT (I  ,K  , L  )

    WT (I,K, L) = WT (I,K, L) + DENWT(JM1)*                            &
                 ( 3._PP         * TZZT (I  ,K  ,  0  )                &
                 - 1._PP/3._PP   * TZZT (I  ,K  ,  1  ) )
!                 ( 35._PP/ 8._PP * TZZT (I  ,K  ,  0  )                &
!                 - 35._PP/24._PP * TZZT (I  ,K  ,  1  )                &
!                 + 21._PP/40._PP * TZZT (I  ,K  ,  2  )                &
!                 -  5._PP/56._PP * TZZT (I  ,K  ,  3  ) )

  END DO
END DO


END SUBROUTINE  INH_V_W_2ND_0

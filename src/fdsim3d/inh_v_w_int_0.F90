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

SUBROUTINE  INH_V_W_INT_0 ( MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,    &
                            X_MIN, X_MAX, Y_MIN, Y_MAX,                &
                            JMT, DENWT, WT, TZZT )

  USE PRECISION   , ONLY:                                              &
                          PP, IK
  USE GRID_MEDIUM , ONLY:                                              &
                          JMNUM

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,                                            INTENT(IN)    :: &
                                MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,&
                                X_MIN, X_MAX, Y_MIN, Y_MAX
  INTEGER(IK), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
  REAL   (PP), DIMENSION (1:JMNUM                  ), INTENT(INOUT) :: &
                                                                  DENWT
  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                   TZZT
  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                     WT

  INTEGER :: I, K, L, JM1

!-----------------------------------------------------------------------

L    = 0

DO  K = MAX(3,Y_MIN), MIN(MYT-2,Y_MAX)
  DO  I = MAX(3,X_MIN), MIN(MXT-2,X_MAX)

    JM1 = JMT (I  ,K  , L  )

    WT (I,K, L) = WT (I,K, L) + DENWT(JM1)*                            &
                 ( 35._PP/ 8._PP * TZZT (I  ,K  ,  0  )                &
                 - 35._PP/24._PP * TZZT (I  ,K  ,  1  )                &
                 + 21._PP/40._PP * TZZT (I  ,K  ,  2  )                &
                 -  5._PP/56._PP * TZZT (I  ,K  ,  3  ) )
  END DO
END DO

END SUBROUTINE  INH_V_W_INT_0

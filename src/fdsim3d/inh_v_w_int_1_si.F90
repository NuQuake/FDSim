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

SUBROUTINE  INH_V_W_INT_1_SI ( MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH, &
                            X_MIN, X_MAX, Y_MIN, Y_MAX, TPML,          &
                            JMT, UT, VT, WT )

  USE PRECISION   , ONLY:                                              &
                          PP, IK
  USE GRID_MEDIUM , ONLY:                                              &
                          L2MZ, LAXZ, LAYZ

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER    ,                                        INTENT(IN)    :: &
                               MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH, &
                               X_MIN, X_MAX, Y_MIN, Y_MAX, TPML
  INTEGER(IK), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                             UT, VT, WT

  INTEGER   :: I, K, L, JM1

!-----------------------------------------------------------------------

L    = 1

DO  K = MAX(3-TPML,Y_MIN), MIN(MYT-2+TPML,Y_MAX)
  DO  I = MAX(3-TPML,X_MIN), MIN(MXT-2+TPML,X_MAX)

    JM1 = JMT (I  ,K  ,0  )

    WT (I,K,L) = WT  (I,K,3)                                           &
               + LAXZ(JM1)*( UT (I+1,K  ,1) - UT (I  ,K  ,1)           &
                           + UT (I+1,K  ,2) - UT (I  ,K  ,2))/L2MZ(JM1)&
               + LAYZ(JM1)*( VT (I  ,K+1,1) - VT (I  ,K  ,1)           &
                           + VT (I  ,K+1,2) - VT (I  ,K  ,2))/L2MZ(JM1)

  END DO
END DO


END SUBROUTINE  INH_V_W_INT_1_SI

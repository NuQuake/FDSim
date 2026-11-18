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

SUBROUTINE  INH_V_UV_INT_1_SI                                          &
                          ( MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,    &
                            X_MIN, X_MAX, Y_MIN, Y_MAX, TPML,          &
                            UT, VT, WT )

  USE PRECISION   , ONLY: PP

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,                                            INTENT(IN)    :: &
                               MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH, &
                               X_MIN, X_MAX, Y_MIN, Y_MAX, TPML
  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                            UT, VT, WT

  INTEGER              :: I, K, L

!-----------------------------------------------------------------------

L    = 1

DO  K = MAX(2-TPML,Y_MIN), MIN(MYT-1+TPML,Y_MAX)
  DO  I = MAX(3-TPML,X_MIN), MIN(MXT-1+TPML,X_MAX)

    UT (I,K, L) = UT (I,K,2) + WT (I,K,2) - WT (I-1,K,2)

  END DO
END DO

DO  K = MAX(3-TPML,Y_MIN), MIN(MYT-1+TPML,Y_MAX)
  DO  I = MAX(2-TPML,X_MIN), MIN(MXT-1+TPML,X_MAX)

    VT (I,K, L) = VT (I,K,2) + WT (I,K,2) - WT (I,K-1,2)

  END DO
END DO

END SUBROUTINE  INH_V_UV_INT_1_SI

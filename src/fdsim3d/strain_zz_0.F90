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

SUBROUTINE  STRAIN_ZZ_0 ( MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,      &
                          X_MIN, X_MAX, Y_MIN, Y_MAX, TPML, C,         &
                          WTM, XZZT, EZZT )

  USE PRECISION   , ONLY: PP
  USE AUXIL       , ONLY: OTA, OTB
  USE INTERFACES  , ONLY: IND

!-----------------------------------------------------------------------

  IMPLICIT NONE


  INTEGER,                                           INTENT (IN)    :: &
                          MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH, TPML,&
                          X_MIN, X_MAX, Y_MIN, Y_MAX

  REAL   (PP),                                       INTENT (IN)    :: &
                                                                     C

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EZZT

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                   WTM

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XZZT

!-----------------------------------------------------------------------

  INTEGER       :: I, K, L, IO, LP1, MXT1, MYT1

!-----------------------------------------------------------------------

L    = 0
LP1  = 1

MXT1  = MXT - 1
MYT1  = MYT - 1

!================================== EZZT ===============================
!----- ALL POINTS COMPUTED BY 2ND ORDER FD APPROX.
DO  K = MAX(2-TPML,Y_MIN), MIN(MYT1+TPML,Y_MAX)
  DO  I = MAX(2-TPML,X_MIN), MIN(MXT1+TPML,X_MAX)

    EZZT(I,K) = ( C*( WTM (I  ,K  , LP1) - WTM (I  ,K  , L  ) ) )

  END DO
END DO

!----- INTERIOR POINTS COMPUTED BY 4TH ORDER FD APPROX.
!DO  K = MAX(2,Y_MIN), MIN(MYT1,Y_MAX)
!  DO  I = MAX(2,X_MIN), MIN(MXT1,X_MAX)

!    EZZT(I,K) = C*( - 11._PP/12._PP * WTM (I,K,0)                      &
!                    + 17._PP/24._PP * WTM (I,K,1)                      &
!                    +  3._PP/ 8._PP * WTM (I,K,2)                      &
!                    -  5._PP/24._PP * WTM (I,K,3)                      &
!                    +  1._PP/24._PP * WTM (I,K,4) )

!  END DO
!END DO

!----- COMPUTATION XZZT AT TIME LEVEL M+1
DO  K = MAX(2-TPML,Y_MIN), MIN(MYT1+TPML,Y_MAX)
  DO  I = MAX(2-TPML,X_MIN), MIN(MXT1+TPML,X_MAX)

    IO  = IND (I, K)
    XZZT(I,K, L ) = OTA(IO)*EZZT(I,K) + OTB(IO)*XZZT(I,K, L )

  END DO
END DO

END SUBROUTINE  STRAIN_ZZ_0

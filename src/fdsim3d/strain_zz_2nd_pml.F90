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

SUBROUTINE STRAIN_ZZ_2ND_PML ( L, LL, C, TPML,                         &
                          MXL, MXH,     MXT,                           &
                          MYL, MYH,     MYT,                           &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX, Y_MIN, Y_MAX,                  &
                          WTM, XZZT, EZZT, PZZHT,                      &
                          ODA_ZH, ODB_ZH, ODC_ZH, ODD_ZH )

  USE PRECISION   , ONLY: PP
  USE AUXIL       , ONLY: OTA, OTB
  USE INTERFACES  , ONLY: IND

!-----------------------------------------------------------------------

  IMPLICIT NONE


  INTEGER,                                            INTENT(IN)    :: &
                                               MXL, MXH,     MXT,      &
                                               MYL, MYH,     MYT,      &
                                               MZL, MZH, L, LL, TPML,  &
                                            X_MIN, X_MAX, Y_MIN, Y_MAX

  REAL   (PP),                                        INTENT(IN)    :: &
                                                                     C

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EZZT

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                   WTM

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XZZT

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH,  0:TPML), INTENT(INOUT) :: &
                                                                 PZZHT

  REAL   (PP), DIMENSION ( 0:TPML, MXL:MXH, MYL:MYH), INTENT(INOUT) :: &
                                        ODA_ZH, ODB_ZH, ODC_ZH, ODD_ZH

!-----------------------------------------------------------------------

  INTEGER       :: I, K, LP1, IO, MXT1, MYT1

  REAL (PP)     :: DZZA

!-----------------------------------------------------------------------

LP1    = L + 1

MXT1  = MXT - 1
MYT1  = MYT - 1

!================================== EZZT ===============================
!==== 2 - MYT1 ====
  DO  K = MAX(2-TPML,Y_MIN), MIN(MYT1+TPML,Y_MAX)
    DO  I = MAX(2-TPML,X_MIN), MIN(MXT1+TPML,X_MAX)

      DZZA          = C*(   WTM (I  ,K  , LP1 ) - WTM (I  ,K  , L ) )

      EZZT(I,K)     = ODA_ZH(LL,I,K) * ( ODD_ZH(LL,I,K)*DZZA           &
                                       + PZZHT(I,K,LL)       )

      PZZHT(I,K,LL) = ODB_ZH(LL,I,K) * PZZHT(I,K,LL)                   &
                    + ODC_ZH(LL,I,K) * DZZA

      IO  = IND (I, K)
      XZZT(I,K, L) = OTA(IO)*EZZT(I,K) + OTB(IO)*XZZT(I,K, L)

    END DO
  END DO


END SUBROUTINE  STRAIN_ZZ_2ND_PML

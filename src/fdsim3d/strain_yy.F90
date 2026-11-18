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

SUBROUTINE  STRAIN_YY (L, MXL, MXH,     MXT,                           &
                          MYL, MYH, PL, MYT, PYH,                      &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX, Y_MIN, Y_MAX,                  &
                          A, B, C, TPML,                               &
                          VTM, XYYT, EYYT, PYYLT, PYYHT,               &
                          ODA_YL, ODB_YL, ODC_YL, ODD_YL,              &
                          ODA_YH, ODB_YH, ODC_YH, ODD_YH )

  USE PRECISION   , ONLY: PP
  USE AUXIL       , ONLY: OTA, OTB
  USE INTERFACES  , ONLY: IND

!-----------------------------------------------------------------------

  IMPLICIT NONE


  INTEGER,                                           INTENT (IN)    :: &
                                               MXL, MXH,     MXT,      &
                                               MYL, MYH, PL, MYT, PYH, &
                                               MZL, MZH,  TPML,   L,   &
                                            X_MIN, X_MAX, Y_MIN, Y_MAX

  REAL   (PP),                                       INTENT (IN)    :: &
                                                               A, B, C

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EYYT

  REAL   (PP), DIMENSION (MXL:MXH, PL :  1, MZL:MZH), INTENT(INOUT) :: &
                                                                 PYYLT

  REAL   (PP), DIMENSION (MXL:MXH, MYT:PYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PYYHT

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                   VTM

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XYYT

  REAL   (PP), DIMENSION ( 0:TPML, MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
        ODA_YL, ODB_YL, ODC_YL, ODD_YL, ODA_YH, ODB_YH, ODC_YH, ODD_YH

!-----------------------------------------------------------------------

  INTEGER       :: I, K, IO, MXT1, MYT1

  REAL (PP)     :: DYYA

!-----------------------------------------------------------------------

MXT1  = MXT - 1
MYT1  = MYT - 1

!================================== EYYT ===============================
!----- FIRST ROW OF PML COMPUTED BY 2ND ORDER FD APPROX.
DO  K = MAX(2-TPML,Y_MIN), MIN(1,Y_MAX)
  DO  I = MAX(2-TPML,X_MIN), MIN(MXT1+TPML,X_MAX)

    DYYA          = C*( VTM (I  ,K+1, L  ) - VTM (I  ,K  , L ) )

    EYYT(I,K)     = ODA_YL(2-K,I,L) * ( ODD_YL(2-K,I,L)*DYYA           &
                                      + PYYLT(I,K,L)         )

    PYYLT(I,K, L )= ODB_YL(2-K,I,L) * PYYLT(I,K,L)                     &
                  + ODC_YL(2-K,I,L) * DYYA

  END DO
END DO

!----- FIRST ROW COMPUTED BY 2ND ORDER FD APPROX.
IF ( ( Y_MIN <= 2 ) .AND. ( 2 <= Y_MAX ) ) THEN
  K = 2
  DO  I = MAX(2-TPML,X_MIN), MIN(MXT1+TPML,X_MAX)

    EYYT(I,K)     =    C*( VTM (I  ,K+1, L  ) - VTM (I  ,K  , L ) )

  END DO
END IF

!----- INTERIOR POINTS COMPUTED BY 4TH ORDER FD APPROX.
DO  K = MAX(3,Y_MIN), MIN(MYT-2,Y_MAX)
  DO  I = MAX(2-TPML,X_MIN), MIN(1,X_MAX)

    EYYT(I,K) = C*( VTM (I  ,K+1, L  ) - VTM (I  ,K  , L  ) )

  END DO
  DO  I = MAX(2,X_MIN), MIN(MXT1,X_MAX)

    EYYT(I,K) = A*( VTM (I  ,K+2, L  ) - VTM (I  ,K-1, L  ) )          &
              + B*( VTM (I  ,K+1, L  ) - VTM (I  ,K  , L  ) )

  END DO
  DO  I = MAX(MXT,X_MIN), MIN(MXT1+TPML,X_MAX)

    EYYT(I,K) = C*( VTM (I  ,K+1, L  ) - VTM (I  ,K  , L  ) )

  END DO
END DO

!----- LAST ROW COMPUTED BY 2ND ORDER FD APPROX.
IF ( ( Y_MIN <= MYT1 ) .AND. ( MYT1 <= Y_MAX ) ) THEN
  K = MYT1
  DO  I = MAX(2-TPML,X_MIN), MIN(MXT1+TPML,X_MAX)

    EYYT(I,K)     =    C*( VTM (I  ,K+1, L  ) - VTM (I  ,K  , L ) )

  END DO
END IF

!----- LAST ROW OF PML COMPUTED BY 2ND ORDER FD APPROX.
DO  K = MAX(MYT,Y_MIN), MIN(MYT1+TPML,Y_MAX)
  DO  I = MAX(2-TPML,X_MIN), MIN(MXT1+TPML,X_MAX)

    DYYA          = C*( VTM (I  ,K+1, L  ) - VTM (I  ,K  , L ) )

    EYYT(I,K)     = ODA_YH(K-MYT1,I,L) * ( ODD_YH(K-MYT1,I,L)*DYYA        &
                                         + PYYHT(I,K,L)            )

    PYYHT(I,K, L )= ODB_YH(K-MYT1,I,L) * PYYHT(I,K,L)                     &
                  + ODC_YH(K-MYT1,I,L) * DYYA

  END DO
END DO

!----- COMPUTATION XYYT AT TIME LEVEL M+1
DO K = MAX(2-TPML,Y_MIN), MIN(MYT1+TPML,Y_MAX)
  DO  I = MAX(2-TPML,X_MIN), MIN(MXT1+TPML,X_MAX)

    IO  = IND (I, K)
    XYYT(I,K, L ) = OTA(IO)*EYYT(I,K) + OTB(IO)*XYYT(I,K, L )

  END DO
END DO

END SUBROUTINE  STRAIN_YY

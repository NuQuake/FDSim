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

SUBROUTINE  STRAIN_XZ_1 ( MXL, MXH, PL, MXT, PXH,                      &
                          MYL, MYH,     MYT,                           &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX, Y_MIN, Y_MAX,                  &
                          A, B, C, TPML,                               &
                          UTM, WTM, XXZT, EXZT,                        &
                          PZXLT, PZXHT,                                &
                          ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,          &
                          ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH )

  USE PRECISION   , ONLY: PP
  USE AUXIL       , ONLY: OTA, OTB
  USE INTERFACES  , ONLY: IND

!-----------------------------------------------------------------------

  IMPLICIT NONE


  INTEGER,                                            INTENT(IN)    :: &
                                               MXL, MXH, PL, MXT, PXH, &
                                               MYL, MYH,     MYT,      &
                                               MZL, MZH, TPML,         &
                                           X_MIN, X_MAX, Y_MIN, Y_MAX

  REAL   (PP),                                        INTENT(IN)    :: &
                                                               A, B, C

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EXZT

  REAL   (PP), DIMENSION (PL :  2, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXLT

  REAL   (PP), DIMENSION (MXT:PXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXHT

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                              UTM, WTM

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XXZT

  REAL   (PP), DIMENSION(0:TPML+1, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
 ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL, ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH

!-----------------------------------------------------------------------

  INTEGER       :: I, K, L1, L , IO, MXT1, MYT1

  REAL (KIND=PP):: DUIKM1, DZXA

!-----------------------------------------------------------------------

L     = 1
L1    = L - 1

MXT1  = MXT - 1
MYT1  = MYT - 1

!================================== EXZT ===============================
!----- FIRST ROW OF PML COMPUTED BY 2ND ORDER FD APPROX.
DO  K = MAX(PL+1,Y_MIN), MIN(1,Y_MAX)
  DO  I = MAX(PL+1,X_MIN), MIN(2,X_MAX)

    DZXA          = C*(   WTM (I  ,K  , L  ) - WTM (I-1,K  , L ) )

    EXZT(I,K)     = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K  , L1) )     &
                  + ODA2_XL(3-I,K,L) * ( ODD2_XL(3-I,K,L)*DZXA         &
                                       + PZXLT (I,K,L)         )

    PZXLT(I,K, L )= ODB2_XL(3-I,K,L) * PZXLT(I,K,L)                    &
                  + ODC2_XL(3-I,K,L) * DZXA

  END DO
  DO  I = MAX(3,X_MIN), MIN(MXT-1,X_MAX)

    EXZT(I,K)   =  C*(   UTM (I  ,K  , L  ) - UTM (I  ,K  , L1)        &
                       + WTM (I  ,K  , L  ) - WTM (I-1,K  , L ) )

  END DO
  DO  I = MAX(MXT,X_MIN), MIN(PXH,X_MAX)

    DZXA          = C*(   WTM (I  ,K  , L  ) - WTM (I-1,K  , L ) )

    EXZT(I,K)     = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K  , L1) )     &
                  + ODA2_XH(I-MXT1,K,L) * ( ODD2_XH(I-MXT1,K,L)*DZXA   &
                                          + PZXHT(I,K,L)             )

    PZXHT(I,K, L )= ODB2_XH(I-MXT1,K,L) * PZXHT(I,K,L)                 &
                  + ODC2_XH(I-MXT1,K,L) * DZXA

  END DO

END DO

DO  K = MAX(2,Y_MIN), MIN(MYT1,Y_MAX)
! FIRST COLUMN OF PML COMPUTED BY 2ND ORDER FD APPROX.
  DO  I = MAX(PL+1,X_MIN), MIN(2,X_MAX)

    DZXA          = C*(   WTM (I  ,K  , L  ) - WTM (I-1,K  , L ) )

    EXZT(I,K)     = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K  , L1) )     &
                  + ODA2_XL(3-I,K,L) * ( ODD2_XL(3-I,K,L)*DZXA         &
                                       + PZXLT (I,K,L)         )

    PZXLT(I,K, L )= ODB2_XL(3-I,K,L) * PZXLT(I,K,L)                    &
                  + ODC2_XL(3-I,K,L) * DZXA

  END DO

! INTERIOR POINTS COMPUTED BY 4TH ORDER FD APPROX.
  DO  I = MAX(3,X_MIN), MIN(MXT-1,X_MAX)

    DUIKM1    = - A*( WTM (I+1,K,0) - WTM (I-2,K,0) )                  &
                - B*( WTM (I  ,K,0) - WTM (I-1,K,0) )

    EXZT(I,K)     =    A*( WTM (I+1,K  , L   ) - WTM (I-2,K  , L   ) ) &
                     + B*( WTM (I  ,K  , L   ) - WTM (I-1,K  , L   ) ) &
                     + C*( UTM (I  ,K  , 1   ) - UTM (I  ,K  , 0   ) )
!                     -       1._PP/ 22._PP * DUIKM1                    &
!                     + C*(                                             &
!                         - 577._PP/528._PP * UTM (I,K,0)               &
!                         + 201._PP/176._PP * UTM (I,K,1)               &
!                         -   9._PP/176._PP * UTM (I,K,2)               &
!                         +   1._PP/528._PP * UTM (I,K,3) )

  END DO

! LAST COLUMN OF PML COMPUTED BY 2ND ORDER FD APPROX.
  DO  I = MAX(MXT,X_MIN), MIN(PXH,X_MAX)

    DZXA          = C*(   WTM (I  ,K  , L  ) - WTM (I-1,K  , L ) )

    EXZT(I,K)     = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K  , L1) )     &
                  + ODA2_XH(I-MXT1,K,L) * ( ODD2_XH(I-MXT1,K,L)*DZXA   &
                                          + PZXHT(I,K,L)             )

    PZXHT(I,K, L )= ODB2_XH(I-MXT1,K,L) * PZXHT(I,K,L)                 &
                  + ODC2_XH(I-MXT1,K,L) * DZXA

  END DO

END DO

!----- LAST COLUMN OF PML COMPUTED BY 2ND ORDER FD APPROX.
DO  K = MAX(MYT,Y_MIN), MIN(MYT1+TPML,Y_MAX)
  DO  I = MAX(PL+1,X_MIN), MIN(2,X_MAX)

    DZXA          = C*(   WTM (I  ,K  , L  ) - WTM (I-1,K  , L ) )

    EXZT(I,K)     = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K  , L1) )     &
                  + ODA2_XL(3-I,K,L) * ( ODD2_XL(3-I,K,L)*DZXA         &
                                       + PZXLT (I,K,L)         )

    PZXLT(I,K, L )= ODB2_XL(3-I,K,L) * PZXLT(I,K,L)                    &
                  + ODC2_XL(3-I,K,L) * DZXA

  END DO
  DO  I = MAX(3,X_MIN), MIN(MXT-1,X_MAX)

    EXZT(I,K)   =  C*(   UTM (I  ,K  , L  ) - UTM (I  ,K  , L1)        &
                       + WTM (I  ,K  , L  ) - WTM (I-1,K  , L ) )

  END DO
  DO  I = MAX(MXT,X_MIN), MIN(PXH,X_MAX)

    DZXA          = C*(   WTM (I  ,K  , L  ) - WTM (I-1,K  , L ) )

    EXZT(I,K)     = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K  , L1) )     &
                  + ODA2_XH(I-MXT1,K,L) * ( ODD2_XH(I-MXT1,K,L)*DZXA   &
                                          + PZXHT(I,K,L)             )

    PZXHT(I,K, L )= ODB2_XH(I-MXT1,K,L) * PZXHT(I,K,L)                 &
                  + ODC2_XH(I-MXT1,K,L) * DZXA

  END DO

END DO

!----- COMPUTATION XXZT AT TIME LEVEL M+1
DO K = MAX(PL+1,Y_MIN), MIN(MYT1+TPML,Y_MAX)
  DO  I = MAX(PL+1,X_MIN), MIN(PXH,X_MAX)

    IO  = IND (I, K)
    XXZT(I,K, L ) = OTA(IO)*EXZT(I,K) + OTB(IO)*XXZT(I,K, L )

  END DO
END DO

END SUBROUTINE  STRAIN_XZ_1

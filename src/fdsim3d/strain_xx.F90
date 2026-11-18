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

SUBROUTINE  STRAIN_XX (L, MXL, MXH, PL, MXT, PXH,                      &
                          MYL, MYH,     MYT,                           &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX, Y_MIN, Y_MAX,                  &
                          A, B, C, TPML,                               &
                          UTM, XXXT, EXXT, PXXLT, PXXHT,               &
                          ODA_XL, ODB_XL, ODC_XL, ODD_XL,              &
                          ODA_XH, ODB_XH, ODC_XH, ODD_XH )

  USE PRECISION   , ONLY: PP
  USE AUXIL       , ONLY: OTA, OTB
  USE INTERFACES  , ONLY: IND

!-----------------------------------------------------------------------

  IMPLICIT NONE


  INTEGER,                                            INTENT(IN)    :: &
                                               MXL, MXH, PL, MXT, PXH, &
                                               MYL, MYH,     MYT,      &
                                               MZL, MZH, L, TPML,      &
                                           X_MIN, X_MAX, Y_MIN, Y_MAX

  REAL   (PP),                                       INTENT (IN)    :: &
                                                               A, B, C

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EXXT

  REAL   (PP), DIMENSION (PL :  1, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PXXLT

  REAL   (PP), DIMENSION (MXT:PXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PXXHT

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                   UTM

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XXXT

  REAL   (PP), DIMENSION ( 0:TPML, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
        ODA_XL, ODB_XL, ODC_XL, ODD_XL, ODA_XH, ODB_XH, ODC_XH, ODD_XH

!-----------------------------------------------------------------------

  INTEGER       :: I, K, IO, MXT1, MYT1

  REAL (PP)     :: DXXA

!-----------------------------------------------------------------------

MXT1  = MXT - 1
MYT1  = MYT - 1

!================================== EXXT ===============================
!----- FIRST ROW OF PML COMPUTED BY 2ND ORDER FD APPROX.
DO  K = MAX(2-TPML,Y_MIN), MIN(1,Y_MAX)
  DO  I = MAX(2-TPML,X_MIN), MIN(1,X_MAX)

    DXXA          = C*( UTM (I+1,K  , L  ) - UTM (I  ,K  , L ) )

    EXXT(I,K)     = ODA_XL(2-I,K,L) * ( ODD_XL(2-I,K,L)*DXXA           &
                                      + PXXLT(I,K,L)         )

    PXXLT(I,K, L )= ODB_XL(2-I,K,L) * PXXLT(I,K,L)                     &
                  + ODC_XL(2-I,K,L) * DXXA

  END DO

  DO  I = MAX(2,X_MIN), MIN(MXT1,X_MAX)

    EXXT(I,K) = C*( UTM (I+1,K  , L  ) - UTM (I  ,K  , L  ) )

  END DO

  DO  I = MAX(MXT,X_MIN), MIN(MXT1+TPML,X_MAX)

    DXXA          = C*( UTM (I+1,K  , L  ) - UTM (I  ,K  , L ) )

    EXXT(I,K)     = ODA_XH(I-MXT1,K,L) * ( ODD_XH(I-MXT1,K,L)*DXXA     &
                                         + PXXHT(I,K,L)            )

    PXXHT(I,K, L )= ODB_XH(I-MXT1,K,L) * PXXHT(I,K,L)                  &
                  + ODC_XH(I-MXT1,K,L) * DXXA

  END DO

END DO

DO  K = MAX(2,Y_MIN), MIN(MYT1,Y_MAX)
! FIRST COLUMN OF PML COMPUTED BY 2ND ORDER FD APPROX.
  DO  I = MAX(2-TPML,X_MIN), MIN(1,X_MAX)

    DXXA          = C*( UTM (I+1,K  , L  ) - UTM (I  ,K  , L ) )

    EXXT(I,K)     = ODA_XL(2-I,K,L) * ( ODD_XL(2-I,K,L)*DXXA           &
                                      + PXXLT(I,K,L)         )

    PXXLT(I,K, L )= ODB_XL(2-I,K,L) * PXXLT(I,K,L)                     &
                  + ODC_XL(2-I,K,L) * DXXA


  END DO

!----- FIRST COLUMN COMPUTED BY 2ND ORDER FD APPROX.
  IF ( ( X_MIN <= 2 ) .AND. ( 2 <= X_MAX ) ) THEN
    I = 2
    EXXT(I,K)     =    C*( UTM (I+1,K  , L  ) - UTM (I  ,K  , L ) )
  END IF

!----- INTERIOR POINTS COMPUTED BY 4TH ORDER FD APPROX.
  DO  I = MAX(3,X_MIN), MIN(MXT-2,X_MAX)

    EXXT(I,K) = A*( UTM (I+2,K  , L  ) - UTM (I-1,K  , L  ) )          &
              + B*( UTM (I+1,K  , L  ) - UTM (I  ,K  , L  ) )

  END DO

!----- LAST COLUMN COMPUTED BY 2ND ORDER FD APPROX.
  IF ( ( X_MIN <= MXT1 ) .AND. ( MXT1 <= X_MAX) ) THEN
    I = MXT1
    EXXT(I,K)     =    C*( UTM (I+1,K  , L  ) - UTM (I  ,K  , L ) )
  END IF

! LAST COLUMN OF PML COMPUTED BY 2ND ORDER FD APPROX.
  DO  I = MAX(MXT,X_MIN), MIN(MXT1+TPML,X_MAX)

    DXXA          = C*( UTM (I+1,K  , L  ) - UTM (I  ,K  , L ) )

    EXXT(I,K)     = ODA_XH(I-MXT1,K,L) * ( ODD_XH(I-MXT1,K,L)*DXXA     &
                                         + PXXHT(I,K,L)            )

    PXXHT(I,K, L )= ODB_XH(I-MXT1,K,L) * PXXHT(I,K,L)                  &
                  + ODC_XH(I-MXT1,K,L) * DXXA


  END DO

END DO

!----- LAST COLUMN OF PML COMPUTED BY 2ND ORDER FD APPROX.
DO  K = MAX(MYT,Y_MIN), MIN(MYT1+TPML,Y_MAX)
  DO  I = MAX(2-TPML,X_MIN), MIN(1,X_MAX)

    DXXA          = C*( UTM (I+1,K  , L  ) - UTM (I  ,K  , L ) )

    EXXT(I,K)     = ODA_XL(2-I,K,L) * ( ODD_XL(2-I,K,L)*DXXA           &
                                      + PXXLT(I,K,L)         )

    PXXLT(I,K, L )= ODB_XL(2-I,K,L) * PXXLT(I,K,L)                     &
                  + ODC_XL(2-I,K,L) * DXXA


  END DO

  DO  I = MAX(2,X_MIN), MIN(MXT1,X_MAX)

    EXXT(I,K) = C*( UTM (I+1,K  , L  ) - UTM (I  ,K  , L  ) )

  END DO

  DO  I = MAX(MXT,X_MIN), MIN(MXT1+TPML,X_MAX)

    DXXA          = C*( UTM (I+1,K  , L  ) - UTM (I  ,K  , L ) )

    EXXT(I,K)     = ODA_XH(I-MXT1,K,L) * ( ODD_XH(I-MXT1,K,L)*DXXA     &
                                         + PXXHT(I,K,L)            )

    PXXHT(I,K, L )= ODB_XH(I-MXT1,K,L) * PXXHT(I,K,L)                  &
                  + ODC_XH(I-MXT1,K,L) * DXXA


  END DO

END DO

!----- COMPUTATION XXXT AT TIME LEVEL M+1
DO K = MAX(2-TPML,Y_MIN), MIN(MYT1+TPML,Y_MAX)
  DO  I = MAX(2-TPML,X_MIN), MIN(MXT1+TPML,X_MAX)

    IO  = IND (I, K)
    XXXT(I,K, L ) = OTA(IO)*EXXT(I,K) + OTB(IO)*XXXT(I,K, L )

  END DO
END DO

END SUBROUTINE  STRAIN_XX

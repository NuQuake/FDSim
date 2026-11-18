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

SUBROUTINE  STRAIN_XY (L, MXL, MXH,     MXT, PXH,                      &
                          MYL, MYH, PL, MYT, PYH,                      &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX, Y_MIN, Y_MAX,                  &
                          A, B, C, TPML,                               &
                          UTM, VTM, XXYT, EXYT,                        &
                          PYXLT, PYXHT,                                &
                          PXYLT, PXYHT,                                &
                          ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,          &
                          ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH,          &
                          ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL,          &
                          ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH )

  USE PRECISION   , ONLY: PP
  USE AUXIL       , ONLY: OTA, OTB
  USE INTERFACES  , ONLY: IND

!-----------------------------------------------------------------------

  IMPLICIT NONE


  INTEGER,                                           INTENT (IN)    :: &
                                               MXL, MXH,     MXT, PXH, &
                                               MYL, MYH, PL, MYT, PYH, &
                                               MZL, MZH, L, TPML,      &
                                            X_MIN, X_MAX, Y_MIN, Y_MAX

  REAL   (PP),                                       INTENT (IN)    :: &
                                                               A, B, C

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EXYT

  REAL   (PP), DIMENSION (PL :  2, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PYXLT

  REAL   (PP), DIMENSION (MXT:PXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PYXHT

  REAL   (PP), DIMENSION (MXL:MXH, MYL:  2, MZL:MZH), INTENT(INOUT) :: &
                                                                 PXYLT

  REAL   (PP), DIMENSION (MXL:MXH, MYT:PYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PXYHT

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                              UTM, VTM

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XXYT

  REAL   (PP), DIMENSION(0:TPML+1, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
 ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL, ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH

  REAL   (PP), DIMENSION(0:TPML+1, MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
 ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL, ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH

!-----------------------------------------------------------------------

  INTEGER       :: I, K, IO, MXT1, MYT1

  REAL (PP)     :: DYXA, DXYA

!-----------------------------------------------------------------------

MXT1  = MXT - 1
MYT1  = MYT - 1

!================================== EXYT ===============================
!----- FIRST ROW OF PML COMPUTED BY 2ND ORDER FD APPROX.
IF ( ( Y_MIN <= 2-TPML ) .AND. ( 2-TPML <= Y_MAX ) ) THEN
  K = 2-TPML
  DO  I = MAX(3-TPML,X_MIN), MIN(2,X_MAX)

    DYXA          = C*(   VTM (I  ,K  , L  ) - VTM (I-1,K  , L  ) )
    DXYA          = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K-1, L  ) )

    EXYT(I,K)     = ODA2_XL(3-I,K,L) * ( ODD2_XL(3-I,K,L)*DYXA         &
                                       + PYXLT(I,K,L)          )       &
                  + ODA2_YL(3-K,I,L) * ( ODD2_YL(3-K,I,L)*DXYA         &
                                       + PXYLT(I,K,L)          )

    PYXLT(I,K, L )= ODB2_XL(3-I,K,L) * PYXLT(I,K,L)                    &
                  + ODC2_XL(3-I,K,L) * DYXA

    PXYLT(I,K, L )= ODB2_YL(3-K,I,L) * PXYLT(I,K,L )                   &
                  + ODC2_YL(3-K,I,L) * DXYA

  END DO

  DO  I = MAX(3,X_MIN), MIN(MXT1,X_MAX)

    DYXA          = C*(   VTM (I  ,K  , L  ) - VTM (I-1,K  , L  ) )
    DXYA          = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K-1, L  ) )

    EXYT(I,K)     =                      DYXA                          &
                  + ODA2_YL(3-K,I,L) * ( ODD2_YL(3-K,I,L)*DXYA         &
                                       + PXYLT(I,K,L)          )

    PXYLT(I,K, L )= ODB2_YL(3-K,I,L) * PXYLT(I,K,L )                   &
                  + ODC2_YL(3-K,I,L) * DXYA


  END DO

  DO  I = MAX(MXT,X_MIN), MIN(MXT1+TPML,X_MAX)

    DYXA          = C*(   VTM (I  ,K  , L  ) - VTM (I-1,K  , L  ) )
    DXYA          = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K-1, L  ) )

    EXYT(I,K)     = ODA2_XH(I-MXT1,K,L) * ( ODD2_XH(I-MXT1,K,L)*DYXA   &
                                          + PYXHT(I,K,L)             ) &
                  + ODA2_YL(3-K   ,I,L) * ( ODD2_YL(3-K   ,I,L)*DXYA   &
                                          + PXYLT(I,K,L)             )

    PYXHT(I,K, L )= ODB2_XH(I-MXT1,K,L) * PYXHT(I,K,L)                 &
                  + ODC2_XH(I-MXT1,K,L) * DYXA

    PXYLT(I,K, L )= ODB2_YL(3-K,I,L) * PXYLT(I,K,L )                   &
                  + ODC2_YL(3-K,I,L) * DXYA

  END DO

  DO  I = MAX(3-TPML,X_MIN), MIN(MXT1+TPML,X_MAX)

    IO = IND (I  ,K)
    XXYT(I  ,K, L ) =         OTA(IO)*EXYT(I,K) + OTB(IO)*XXYT(I,K, L )

  END DO

END IF

!=== K = 3-TPML, 2
DO  K = MAX(3-TPML,Y_MIN), MIN(2,Y_MAX)
  DO  I = MAX(2-TPML,X_MIN), MIN(2,X_MAX)

    DYXA          = C*(   VTM (I  ,K  , L  ) - VTM (I-1,K  , L  ) )
    DXYA          = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K-1, L  ) )

    EXYT(I,K)     = ODA2_XL(3-I,K,L) * ( ODD2_XL(3-I,K,L)*DYXA         &
                                       + PYXLT(I,K,L)            )     &
                  + ODA2_YL(3-K,I,L) * ( ODD2_YL(3-K,I,L)*DXYA         &
                                       + PXYLT(I,K,L) )

    PYXLT(I,K, L )= ODB2_XL(3-I,K,L) * PYXLT(I,K,L)                    &
                  + ODC2_XL(3-I,K,L) * DYXA

    PXYLT(I,K, L )= ODB2_YL(3-K,I,L) * PXYLT(I,K,L )                   &
                  + ODC2_YL(3-K,I,L) * DXYA

  END DO

  DO  I = MAX(3,X_MIN), MIN(MXT1,X_MAX)

    DYXA          = C*(   VTM (I  ,K  , L  ) - VTM (I-1,K  , L  ) )
    DXYA          = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K-1, L  ) )

    EXYT(I,K)     =                      DYXA                          &
                  + ODA2_YL(3-K,I,L) * ( ODD2_YL(3-K,I,L)*DXYA         &
                                       + PXYLT(I,K,L)          )

    PXYLT(I,K, L )= ODB2_YL(3-K,I,L) * PXYLT(I,K,L )                   &
                  + ODC2_YL(3-K,I,L) * DXYA

  END DO

  DO  I = MAX(MXT,X_MIN), MIN(MXT+TPML,X_MAX)

    DYXA          = C*(   VTM (I  ,K  , L  ) - VTM (I-1,K  , L  ) )
    DXYA          = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K-1, L  ) )

    EXYT(I,K)     = ODA2_XH(I-MXT1,K,L) * ( ODD2_XH(I-MXT1,K,L)*DYXA   &
                                          + PYXHT(I,K,L)             ) &
                  + ODA2_YL(3-K   ,I,L) * ( ODD2_YL(3-K   ,I,L)*DXYA   &
                                          + PXYLT(I,K,L)             )

    PYXHT(I,K, L )= ODB2_XH(I-MXT1,K,L) * PYXHT(I,K,L)                 &
                  + ODC2_XH(I-MXT1,K,L) * DYXA

    PXYLT(I,K, L )= ODB2_YL(3-K,I,L) * PXYLT(I,K,L )                   &
                  + ODC2_YL(3-K,I,L) * DXYA

  END DO

  DO  I = MAX(2-TPML,X_MIN), MIN(MXT+TPML,X_MAX)

    IO = IND (I  ,K)
    XXYT(I  ,K, L ) =         OTA(IO)*EXYT(I,K) + OTB(IO)*XXYT(I,K, L )

  END DO
END DO

!=== K = 3, MYT1
DO  K = MAX(3,Y_MIN), MIN(MYT1,Y_MAX)
  DO  I = MAX(2-TPML,X_MIN), MIN(2,X_MAX)

    DYXA          = C*(   VTM (I  ,K  , L  ) - VTM (I-1,K  , L  ) )
    DXYA          = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K-1, L  ) )

    EXYT(I,K)     = ODA2_XL(3-I,K,L) * ( ODD2_XL(3-I,K,L)*DYXA         &
                                       + PYXLT(I,K,L)            )     &
                  +                      DXYA

    PYXLT(I,K, L )= ODB2_XL(3-I,K,L) * PYXLT(I,K,L)                    &
                  + ODC2_XL(3-I,K,L) * DYXA

  END DO

!----- INTERIOR POINTS COMPUTED BY 4TH ORDER FD APPROX.
  DO  I = MAX(3,X_MIN), MIN(MXT1,X_MAX)

    EXYT(I,K) =          A*( UTM (I  ,K+1, L  ) - UTM (I  ,K-2, L  )   &
                           + VTM (I+1,K  , L  ) - VTM (I-2,K  , L  ) ) &
                       + B*( UTM (I  ,K  , L  ) - UTM (I  ,K-1, L  )   &
                           + VTM (I  ,K  , L  ) - VTM (I-1,K  , L  ) )

  END DO

  DO  I = MAX(MXT,X_MIN), MIN(MXT+TPML,X_MAX)

    DYXA          = C*(   VTM (I  ,K  , L  ) - VTM (I-1,K  , L  ) )
    DXYA          = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K-1, L  ) )

    EXYT(I,K)     = ODA2_XH(I-MXT1,K,L) * ( ODD2_XH(I-MXT1,K,L)*DYXA   &
                                          + PYXHT(I,K,L)             ) &
                  +                         DXYA

    PYXHT(I,K, L )= ODB2_XH(I-MXT1,K,L) * PYXHT(I,K,L)                 &
                  + ODC2_XH(I-MXT1,K,L) * DYXA

  END DO

  DO  I = MAX(2-TPML,X_MIN), MIN(MXT+TPML,X_MAX)

    IO  = IND (I, K)
    XXYT(I,K, L ) = OTA(IO)*EXYT(I,K) + OTB(IO)*XXYT(I,K, L )

  END DO

END DO

!----- LAST ROW COMPUTED BY 2ND ORDER FD APPROX.
!=== K = MYT, MYT1+TPML
DO  K = MAX(MYT,Y_MIN), MIN(MYT1+TPML,Y_MAX)
  DO  I = MAX(2-TPML,X_MIN), MIN(2,X_MAX)

    DYXA          = C*(   VTM (I  ,K  , L  ) - VTM (I-1,K  , L  ) )
    DXYA          = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K-1, L  ) )

    EXYT(I,K)     = ODA2_XL(3-I   ,K,L) * ( ODD2_XL(3-I   ,K,L)*DYXA   &
                                          + PYXLT(I,K,L)             ) &
                  + ODA2_YH(K-MYT1,I,L) * ( ODD2_YH(K-MYT1,I,L)*DXYA   &
                                          + PXYHT(I,K,L)             )

    PYXLT(I,K, L )= ODB2_XL(3-I   ,K,L) * PYXLT(I,K,L)                 &
                  + ODC2_XL(3-I   ,K,L) * DYXA

    PXYHT(I,K, L )= ODB2_YH(K-MYT1,I,L) * PXYHT(I,K,L )                &
                  + ODC2_YH(K-MYT1,I,L) * DXYA

  END DO

  DO  I = MAX(3,X_MIN), MIN(MXT1,X_MAX)

    DYXA          = C*(   VTM (I  ,K  , L  ) - VTM (I-1,K  , L  ) )
    DXYA          = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K-1, L  ) )

    EXYT(I,K)     =                         DYXA                       &
                  + ODA2_YH(K-MYT1,I,L) * ( ODD2_YH(K-MYT1,I,L)*DXYA   &
                                          + PXYHT(I,K,L)             )

    PXYHT(I,K, L )= ODB2_YH(K-MYT1,I,L) * PXYHT(I,K,L )                &
                  + ODC2_YH(K-MYT1,I,L) * DXYA

  END DO

  DO  I = MAX(MXT,X_MIN), MIN(MXT+TPML,X_MAX)

    DYXA          = C*(   VTM (I  ,K  , L  ) - VTM (I-1,K  , L  ) )
    DXYA          = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K-1, L  ) )

    EXYT(I,K)     = ODA2_XH(I-MXT1,K,L) * ( ODD2_XH(I-MXT1,K,L)*DYXA   &
                                          + PYXHT(I,K,L)             ) &
                  + ODA2_YH(K-MYT1,I,L) * ( ODD2_YH(K-MYT1,I,L)*DXYA   &
                                          + PXYHT(I,K,L)             )

    PYXHT(I,K, L )= ODB2_XH(I-MXT1,K,L) * PYXHT(I,K,L)                 &
                  + ODC2_XH(I-MXT1,K,L) * DYXA

    PXYHT(I,K, L )= ODB2_YH(K-MYT1,I,L) * PXYHT(I,K,L )                &
                  + ODC2_YH(K-MYT1,I,L) * DXYA

  END DO

  DO  I = MAX(2-TPML,X_MIN), MIN(MXT+TPML,X_MAX)

    IO = IND (I  ,K)
    XXYT(I  ,K, L ) =         OTA(IO)*EXYT(I,K) + OTB(IO)*XXYT(I,K, L )

  END DO
END DO

!=== K = MYT+TPML
IF ( ( Y_MIN <= MYT+TPML ) .AND. ( MYT+TPML <= Y_MAX ) ) THEN
  K = MYT+TPML
  DO  I = MAX(3-TPML,X_MIN), MIN(2,X_MAX)

    DYXA          = C*(   VTM (I  ,K  , L  ) - VTM (I-1,K  , L  ) )
    DXYA          = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K-1, L  ) )

    EXYT(I,K)     = ODA2_XL(3-I   ,K,L) * ( ODD2_XL(3-I   ,K,L)*DYXA   &
                                          + PYXLT(I,K,L)             ) &
                  + ODA2_YH(K-MYT1,I,L) * ( ODD2_YH(K-MYT1,I,L)*DXYA   &
                                          + PXYHT(I,K,L)             )

    PYXLT(I,K, L )= ODB2_XL(3-I   ,K,L) * PYXLT(I,K,L)                 &
                  + ODC2_XL(3-I   ,K,L) * DYXA

    PXYHT(I,K, L )= ODB2_YH(K-MYT1,I,L) * PXYHT(I,K,L )                &
                  + ODC2_YH(K-MYT1,I,L) * DXYA

  END DO

  DO  I = MAX(3,X_MIN), MIN(MXT1,X_MAX)

    DYXA          = C*(   VTM (I  ,K  , L  ) - VTM (I-1,K  , L  ) )
    DXYA          = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K-1, L  ) )

    EXYT(I,K)     =                         DYXA                       &
                  + ODA2_YH(K-MYT1,I,L) * ( ODD2_YH(K-MYT1,I,L)*DXYA   &
                                          + PXYHT(I,K,L)             )

    PXYHT(I,K, L )= ODB2_YH(K-MYT1,I,L) * PXYHT(I,K,L )                &
                  + ODC2_YH(K-MYT1,I,L) * DXYA

  END DO

  DO  I = MAX(MXT,X_MIN), MIN(MXT1+TPML,X_MAX)

    DYXA          = C*(   VTM (I  ,K  , L  ) - VTM (I-1,K  , L  ) )
    DXYA          = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K-1, L  ) )

    EXYT(I,K)     = ODA2_XH(I-MXT1,K,L) * ( ODD2_XH(I-MXT1,K,L)*DYXA   &
                                          + PYXHT(I,K,L)             ) &
                  + ODA2_YH(K-MYT1,I,L) * ( ODD2_YH(K-MYT1,I,L)*DXYA   &
                                          + PXYHT(I,K,L)             )

    PYXHT(I,K, L )= ODB2_XH(I-MXT1,K,L) * PYXHT(I,K,L)                 &
                  + ODC2_XH(I-MXT1,K,L) * DYXA

    PXYHT(I,K, L )= ODB2_YH(K-MYT1,I,L) * PXYHT(I,K,L )                &
                  + ODC2_YH(K-MYT1,I,L) * DXYA

  END DO

  DO  I = MAX(3-TPML,X_MIN), MIN(MXT1+TPML,X_MAX)

    IO = IND (I  ,K)
    XXYT(I  ,K, L ) =         OTA(IO)*EXYT(I,K) + OTB(IO)*XXYT(I,K, L )

  END DO
END IF


END SUBROUTINE  STRAIN_XY

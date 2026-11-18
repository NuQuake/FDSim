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

SUBROUTINE STRAIN_YZ_2ND_PML ( L, LL, C, TPML,                         &
                          MXL, MXH,     MXT,                           &
                          MYL, MYH, PL, MYT, PYH,                      &
                          MZL, MZH,                                    &
                       X_MIN, X_MAX, Y_MIN, Y_MAX,                     &
                       VTM, WTM, XYZT, EYZT,                           &
                          PZYLT, PZYHT,                                &
                                 PYZHT,                                &
                          ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL,          &
                          ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH,          &
                          ODA2_ZH, ODB2_ZH, ODC2_ZH, ODD2_ZH )

  USE PRECISION   , ONLY: PP
  USE AUXIL       , ONLY: OTA, OTB
  USE INTERFACES  , ONLY: IND

!-----------------------------------------------------------------------

  IMPLICIT NONE


  INTEGER,                                            INTENT(IN)    :: &
                                               MXL, MXH, PL, MXT, PYH, &
                                               MYL, MYH,     MYT,      &
                                               MZL, MZH, L, LL, TPML,  &
                                           X_MIN, X_MAX, Y_MIN, Y_MAX

  REAL   (PP),                                        INTENT(IN)    :: &
                                                               C

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EYZT

  REAL   (PP), DIMENSION (MXL:MXH, PL :  2, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZYLT

  REAL   (PP), DIMENSION (MXL:MXH, MYT:PYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZYHT

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH,0:TPML+1), INTENT(INOUT) :: &
                                                                 PYZHT

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                              VTM, WTM

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XYZT

  REAL   (PP), DIMENSION(0:TPML+1, MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
 ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL, ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH

  REAL   (PP), DIMENSION(0:TPML+1, MXL:MXH, MYL:MYH), INTENT(INOUT) :: &
                                     ODA2_ZH, ODB2_ZH, ODC2_ZH, ODD2_ZH

!-----------------------------------------------------------------------

  INTEGER       :: I, K, IO, LP1 , L1  , L2  , MXT1, MYT1

  REAL(PP)      :: DZYA, DYZA

!-----------------------------------------------------------------------

LP1  = L  + 1
L1   = L  - 1
L2   = L  - 2

MXT1 = MXT - 1
MYT1 = MYT - 1

!================================== EYZT ===============================
!----- FIRST ROW OF PML COMPUTED BY 2ND ORDER FD APPROX.
DO  K = MAX(2-TPML,Y_MIN), MIN(2,Y_MAX)
  DO  I = MAX(2-TPML,X_MIN), MIN(MXT1+TPML,X_MAX)

    DZYA          = C*(   WTM (I  ,K  , L ) - WTM (I  ,K-1, L ) )
    DYZA          = C*(   VTM (I  ,K  , L ) - VTM (I  ,K  , L1) )

    EYZT(I,K)     = ODA2_ZH(LL ,I,K) * ( ODD2_ZH(LL ,I,K)*DYZA         &
                                       + PYZHT (I,K,LL)         )      &
                  + ODA2_YL(3-K,I,L) * ( ODD2_YL(3-K,I,L)*DZYA         &
                                       + PZYLT (I,K,L )         )

    PZYLT(I,K, L )= ODB2_YL(3-K,I,L) * PZYLT(I,K,L)                    &
                  + ODC2_YL(3-K,I,L) * DZYA

    PYZHT(I,K,LL) = ODB2_ZH(LL,I,K) * PYZHT(I,K,LL)                    &
                  + ODC2_ZH(LL,I,K) * DYZA

  END DO
END DO

!----- INTERIOR POINTS COMPUTED BY 4TH ORDER FD APPROX.
DO  K = MAX(3,Y_MIN), MIN(MYT1,Y_MAX)
  DO  I = MAX(2-TPML,X_MIN), MIN(MXT1+TPML,X_MAX)

    DZYA          = C*(   WTM (I  ,K  , L ) - WTM (I  ,K-1, L ) )
    DYZA          = C*(   VTM (I  ,K  , L ) - VTM (I  ,K  , L1) )

    EYZT(I,K)     = ODA2_ZH(LL ,I,K) * ( ODD2_ZH(LL ,I,K)*DYZA         &
                                       + PYZHT (I,K,LL)            )   &
                  +                      DZYA

    PYZHT(I,K,LL) = ODB2_ZH(LL,I,K) * PYZHT(I,K,LL)                    &
                  + ODC2_ZH(LL,I,K) * DYZA

  END DO
END DO

!----- LAST ROW OF PML COMPUTED BY 2ND ORDER FD APPROX.
DO  K = MAX(MYT,Y_MIN), MIN(MYT+TPML,Y_MAX)
  DO  I = MAX(2-TPML,X_MIN), MIN(MXT1+TPML,X_MAX)

    DZYA          = C*(   WTM (I  ,K  , L ) - WTM (I  ,K-1, L ) )
    DYZA          = C*(   VTM (I  ,K  , L ) - VTM (I  ,K  , L1) )

    EYZT(I,K)     = ODA2_ZH(LL    ,I,K) * ( ODD2_ZH(LL    ,I,K)*DYZA   &
                                          + PYZHT (I,K,LL)           ) &
                  + ODA2_YH(K-MYT1,I,L) * ( ODD2_YH(K-MYT1,I,L)*DZYA   &
                                          + PZYHT (I,K,L )           )

    PZYHT(I,K, L )= ODB2_YH(K-MYT1,I,L) * PZYHT(I,K,L)                 &
                  + ODC2_YH(K-MYT1,I,L) * DZYA

    PYZHT(I,K,LL) = ODB2_ZH(LL,I,K) * PYZHT(I,K,LL)                    &
                  + ODC2_ZH(LL,I,K) * DYZA

  END DO
END DO

!----- COMPUTATION XYZT AT TIME LEVEL M+1
DO K = MAX(2-TPML,Y_MIN), MIN(MYT+TPML,Y_MAX)
  DO  I = MAX(2-TPML,X_MIN), MIN(MXT1+TPML,X_MAX)

    IO  = IND (I, K)
    XYZT(I,K, L ) = OTA(IO)*EYZT(I,K) + OTB(IO)*XYZT(I,K, L )

  END DO
END DO

END SUBROUTINE  STRAIN_YZ_2ND_PML

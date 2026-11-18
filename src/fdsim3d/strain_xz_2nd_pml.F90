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

SUBROUTINE STRAIN_XZ_2ND_PML ( L, LL, C, TPML,                         &
                          MXL, MXH, PL, MXT, PXH,                      &
                          MYL, MYH,     MYT,                           &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX, Y_MIN, Y_MAX,                  &
                           UTM, WTM, XXZT, EXZT,                       &
                          PZXLT, PZXHT, PXZHT,                         &
                          ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,          &
                          ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH,          &
                          ODA2_ZH, ODB2_ZH, ODC2_ZH, ODD2_ZH )

  USE PRECISION   , ONLY: PP
  USE AUXIL       , ONLY: OTA, OTB
  USE INTERFACES  , ONLY: IND

!-----------------------------------------------------------------------

  IMPLICIT NONE


  INTEGER,                                            INTENT(IN)    :: &
                                               MXL, MXH, PL, MXT, PXH, &
                                               MYL, MYH,     MYT,      &
                                               MZL, MZH, L, LL, TPML,  &
                                           X_MIN, X_MAX, Y_MIN, Y_MAX

  REAL   (PP),                                        INTENT(IN)    :: &
                                                               C

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EXZT

  REAL   (PP), DIMENSION (PL :  2, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXLT

  REAL   (PP), DIMENSION (MXT:PXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXHT

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH,0:TPML+1), INTENT(INOUT) :: &
                                                                 PXZHT

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                              UTM, WTM

  REAL   (PP), DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XXZT

  REAL   (PP), DIMENSION(0:TPML+1, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
 ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL, ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH

  REAL   (PP), DIMENSION(0:TPML+1, MXL:MXH, MYL:MYH), INTENT(INOUT) :: &
                                     ODA2_ZH, ODB2_ZH, ODC2_ZH, ODD2_ZH

!-----------------------------------------------------------------------

  INTEGER       :: I, K, IO, LP1 , L1  , L2  , MXT1, MYT1

  REAL(PP)      :: DZXA, DXZA

!-----------------------------------------------------------------------

LP1  = L  + 1
L1   = L  - 1
L2   = L  - 2

MXT1  = MXT - 1
MYT1  = MYT - 1

!================================== EXZT ===============================
!----- FIRST ROW OF PML COMPUTED BY 2ND ORDER FD APPROX.
DO  K = MAX(2-TPML,Y_MIN), MIN(MYT1+TPML,Y_MAX)
  DO  I = MAX(2-TPML,X_MIN), MIN(2,X_MAX)

    DZXA          = C*(   WTM (I  ,K  , L  ) - WTM (I-1,K  , L ) )
    DXZA          = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K  , L1) )

    EXZT(I,K)     = ODA2_ZH(LL ,I,K) * ( ODD2_ZH(LL ,I,K)*DXZA         &
                                       + PXZHT (I,K,LL)            )   &
                  + ODA2_XL(3-I,K,L) * ( ODD2_XL(3-I,K,L)*DZXA         &
                                       + PZXLT (I,K,L )            )

    PXZHT(I,K, LL)= ODB2_ZH(LL ,I,K) * PXZHT(I,K,LL)                   &
                  + ODC2_ZH(LL ,I,K) * DXZA

    PZXLT(I,K, L )= ODB2_XL(3-I,K,L) * PZXLT(I,K,L)                    &
                  + ODC2_XL(3-I,K,L) * DZXA

  END DO
  DO  I = MAX(3,X_MIN), MIN(MXT1,X_MAX)

    DZXA          = C*(   WTM (I  ,K  , L  ) - WTM (I-1,K  , L ) )
    DXZA          = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K  , L1) )

    EXZT(I,K)     = ODA2_ZH(LL ,I,K) * ( ODD2_ZH(LL ,I,K)*DXZA         &
                                       + PXZHT (I,K,LL)            )   &
                  +                      DZXA

    PXZHT(I,K, LL)= ODB2_ZH(LL ,I,K) * PXZHT(I,K,LL)                   &
                  + ODC2_ZH(LL ,I,K) * DXZA

  END DO
  DO  I = MAX(MXT,X_MIN), MIN(MXT+TPML,X_MAX)

    DZXA          = C*(   WTM (I  ,K  , L  ) - WTM (I-1,K  , L ) )
    DXZA          = C*(   UTM (I  ,K  , L  ) - UTM (I  ,K  , L1) )

    EXZT(I,K)     = ODA2_ZH(LL    ,I,K) * ( ODD2_ZH(LL    ,I,K)*DXZA   &
                                          + PXZHT (I,K,LL)           ) &
                  + ODA2_XH(I-MXT1,K,L) * ( ODD2_XH(I-MXT1,K,L)*DZXA   &
                                          + PZXHT (I,K,L )           )

    PXZHT(I,K, LL)= ODB2_ZH(LL ,I,K) * PXZHT(I,K,LL)                   &
                  + ODC2_ZH(LL ,I,K) * DXZA

    PZXHT(I,K, L )= ODB2_XH(I-MXT1,K,L) * PZXHT(I,K,L)                 &
                  + ODC2_XH(I-MXT1,K,L) * DZXA

  END DO

END DO


!----- COMPUTATION XXZT AT TIME LEVEL M+1
DO K = MAX(2-TPML,Y_MIN), MIN(MYT1+TPML,Y_MAX)
  DO  I = MAX(2-TPML,X_MIN), MIN(MXT+TPML,X_MAX)

    IO  = IND (I, K)
    XXZT(I,K, L ) = OTA(IO)*EXZT(I,K) + OTB(IO)*XXZT(I,K, L )

  END DO
END DO

END SUBROUTINE  STRAIN_XZ_2ND_PML

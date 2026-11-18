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

SUBROUTINE  HPLANE_1_VELOCITY

  USE PRECISION   , ONLY:                                              &
                          PP
  USE GRID_MEDIUM , ONLY:                                              &
                          JM ,                                         &
                          DENU , DENV , DENW
  USE CONTROL_DATA, ONLY:                                              &
                          MX, MY, TPML, STRESS_IMAGING
  USE AUXIL       , ONLY:                                              &
                          MX1 , MY1 ,                                  &
                          PMLL, PMLXH, PMLYH,                          &
                           X_MIN  , X_MAX  ,                           &
                           Y_MIN  , Y_MAX  ,                           &
                           X_MIN_D, X_MAX_D,                           &
                           Y_MIN_D, Y_MAX_D,                           &
                           Z_MIN_D, Z_MAX_D
  USE WAVEFIELD   , ONLY:                                              &
                          UM, VM, WM
  USE STRESSFIELD , ONLY:                                              &
                          TXX , TYY , TZZ ,                            &
                          TXY ,       TXZ ,       TYZ
  USE NONREF_BOUND, ONLY:                                              &
                          AXU , AYU , AXV , AYV , AXW , AYW ,          &
                                    SX_U  , SX_UP  , SY_U  , SY_UP  ,  &
                                    SX_V  , SX_VP  , SY_V  , SY_VP  ,  &
                                    SX_W  , SX_WP  , SY_W  , SY_WP
  USE PML         , ONLY:                                              &
            ODA_XL ,  ODB_XL ,  ODC_XL ,  ODA_XH ,  ODB_XH ,  ODC_XH , &
            ODD_XL ,                      ODD_XH ,                     &
            ODA_YL ,  ODB_YL ,  ODC_YL ,  ODA_YH ,  ODB_YH ,  ODC_YH , &
            ODD_YL ,                      ODD_YH ,                     &
           ODA2_XL , ODB2_XL , ODC2_XL , ODA2_XH , ODB2_XH , ODC2_XH , &
           ODD2_XL ,                     ODD2_XH ,                     &
           ODA2_YL , ODB2_YL , ODC2_YL , ODA2_YH , ODB2_YH , ODC2_YH , &
           ODD2_YL ,                     ODD2_YH ,                     &
                                            RXXL  , RXXH  ,            &
                                            RXYL  , RXYH  ,            &
                                            RYYL  , RYYH  ,            &
                                            RYXL  , RYXH  ,            &
                                            RZXL  , RZXH  ,            &
                                            RZYL  , RZYH
  USE INTERFACES  , ONLY:                                              &
                          INH_V_UV_INT     , INH_V_W_INT_1,            &
                          INH_V_UV_INT_1_SI, INH_V_W_INT_1_SI,         &
                          INH_V_U_2ND ,INH_V_V_2ND, INH_V_W_2ND,       &
                          STO_BOR, INH_NR

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER             :: L

!-----------------------------------------------------------------------

  L   = 1

IF (STRESS_IMAGING) THEN
  CALL INH_V_UV_INT_1_SI                                               &
                     ( MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,        &
                       MYL = Y_MIN_D, MYH = Y_MAX_D, MYT = MY ,        &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX , TPML = TPML,   &
                       UT = UM, VT = VM, WT = WM )

  CALL INH_V_W_INT_1_SI                                                &
                     ( MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,        &
                       MYL = Y_MIN_D, MYH = Y_MAX_D, MYT = MY ,        &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX , TPML = TPML,   &
                       JMT = JM , UT = UM, VT = VM, WT = WM )
ELSE

  CALL  STO_BOR    ( MXL = X_MIN_D, MXH = X_MAX_D,                     &
                     MYL = Y_MIN_D, MYH = Y_MAX_D,                     &
                     PXH = PMLXH , PYH = PMLYH ,                       &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX ,                &
                     SX_T = SX_U  (:,L ), SY_T = SY_U  (:,L ),         &
                     SX_TP= SX_UP (:   ), SY_TP= SY_UP (:   ),         &
                     T = UM (:,:,L ), I1 = 2, K1 = 1, TPML = TPML )

  CALL  STO_BOR    ( MXL = X_MIN_D, MXH = X_MAX_D,                     &
                     MYL = Y_MIN_D, MYH = Y_MAX_D,                     &
                     PXH = PMLXH , PYH = PMLYH ,                       &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX ,                &
                     SX_T = SX_V  (:,L ), SY_T = SY_V  (:,L ),         &
                     SX_TP= SX_VP (:   ), SY_TP= SY_VP (:   ),         &
                     T = VM (:,:,L ), I1 = 1, K1 = 2, TPML = TPML )

  CALL  STO_BOR    ( MXL = X_MIN_D, MXH = X_MAX_D,                     &
                     MYL = Y_MIN_D, MYH = Y_MAX_D,                     &
                     PXH = PMLXH , PYH = PMLYH ,                       &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX ,                &
                     SX_T = SX_W  (:,L ), SY_T = SY_W  (:,L ),         &
                     SX_TP= SX_WP (:   ), SY_TP= SY_WP (:   ),         &
                     T = WM (:,:,L ), I1 = 1, K1 = 1, TPML = TPML )


  CALL  INH_V_UV_INT                                                   &
              ( L = L, MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,        &
                       MYL = Y_MIN_D, MYH = Y_MAX_D, MYT = MY ,        &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX ,                &
                       JMT = JM , DENUT = DENU , DENVT = DENV ,        &
                       UT = UM, VT = VM,                               &
                       TXXT = TXX , TXYT = TXY , TXZT = TXZ ,          &
                       TYYT = TYY ,              TYZT = TYZ  )

  CALL  INH_V_W_INT_1                                                  &
                   (   MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,        &
                       MYL = Y_MIN_D, MYH = Y_MAX_D, MYT = MY ,        &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX ,                &
                       JMT = JM , DENWT = DENW ,                       &
                       WT = WM,                                        &
                       TZZT = TZZ , TXZT = TXZ , TYZT = TYZ  )

  CALL  INH_V_U_2ND                                                    &
              ( L = L, MXL = X_MIN_D, MXH = X_MAX_D, PL  = PMLL,       &
                       MXT = MX, MX1 = MX1,          PXH = PMLXH,      &
                       MYL = Y_MIN_D, MYH = Y_MAX_D,                   &
                       MYT = MY, MY1 = MY1,          PYH = PMLYH,      &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,  TPML = TPML,     &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX ,                &
                       JMT = JM , DENUT = DENU ,                       &
                       UT = UM,                                        &
                       TXXT = TXX ,              TXYT = TXY ,          &
                       TXZT = TXZ ,                                    &
                       RXXLT = RXXL , RXXHT = RXXH ,                   &
                       RXYLT = RXYL , RXYHT = RXYH ,                   &
           ODA2_XL = ODA2_XL , ODB2_XL = ODB2_XL ,                     &
           ODC2_XL = ODC2_XL , ODD2_XL = ODD2_XL ,                     &
           ODA2_XH = ODA2_XH , ODB2_XH = ODB2_XH ,                     &
           ODC2_XH = ODC2_XH , ODD2_XH = ODD2_XH ,                     &
            ODA_YL = ODA_YL  ,  ODB_YL = ODB_YL  ,                     &
            ODC_YL = ODC_YL  ,  ODD_YL = ODD_YL  ,                     &
            ODA_YH = ODA_YH  ,  ODB_YH = ODB_YH  ,                     &
            ODC_YH = ODC_YH  ,  ODD_YH = ODD_YH  )

  CALL  INH_V_V_2ND                                                    &
              ( L = L, MXL = X_MIN_D, MXH = X_MAX_D, PL  = PMLL,       &
                       MXT = MX, MX1 = MX1,          PXH = PMLXH,      &
                       MYL = Y_MIN_D, MYH = Y_MAX_D,                   &
                       MYT = MY, MY1 = MY1,          PYH = PMLYH,      &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,  TPML = TPML,     &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX ,                &
                       JMT = JM ,                DENVT = DENV ,        &
                                VT = VM,                               &
                                    TYYT = TYY , TXYT = TXY ,          &
                                    TYZT = TYZ ,                       &
                       RYYLT = RYYL , RYYHT = RYYH ,                   &
                       RYXLT = RYXL , RYXHT = RYXH ,                   &
            ODA_XL = ODA_XL  ,  ODB_XL = ODB_XL  ,                     &
            ODC_XL = ODC_XL  ,  ODD_XL = ODD_XL  ,                     &
            ODA_XH = ODA_XH  ,  ODB_XH = ODB_XH  ,                     &
            ODC_XH = ODC_XH  ,  ODD_XH = ODD_XH  ,                     &
           ODA2_YL = ODA2_YL , ODB2_YL = ODB2_YL ,                     &
           ODC2_YL = ODC2_YL , ODD2_YL = ODD2_YL ,                     &
           ODA2_YH = ODA2_YH , ODB2_YH = ODB2_YH ,                     &
           ODC2_YH = ODC2_YH , ODD2_YH = ODD2_YH )

  CALL  INH_V_W_2ND                                                    &
              ( L = L, MXL = X_MIN_D, MXH = X_MAX_D, PL  = PMLL,       &
                       MXT = MX, MX1 = MX1,          PXH = PMLXH,      &
                       MYL = Y_MIN_D, MYH = Y_MAX_D,                   &
                       MYT = MY, MY1 = MY1,          PYH = PMLYH,      &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,  TPML = TPML,     &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX ,                &
                       JMT = JM , DENWT = DENW , WT = WM,              &
                       TZZT = TZZ , TXZT = TXZ , TYZT = TYZ ,          &
                       RZXLT = RZXL , RZXHT = RZXH ,                   &
                       RZYLT = RZYL , RZYHT = RZYH ,                   &
            ODA_XL = ODA_XL  ,  ODB_XL = ODB_XL  ,                     &
            ODC_XL = ODC_XL  ,  ODD_XL = ODD_XL  ,                     &
            ODA_XH = ODA_XH  ,  ODB_XH = ODB_XH  ,                     &
            ODC_XH = ODC_XH  ,  ODD_XH = ODD_XH  ,                     &
            ODA_YL = ODA_YL  ,  ODB_YL = ODB_YL  ,                     &
            ODC_YL = ODC_YL  ,  ODD_YL = ODD_YL  ,                     &
            ODA_YH = ODA_YH  ,  ODB_YH = ODB_YH  ,                     &
            ODC_YH = ODC_YH  ,  ODD_YH = ODD_YH  )


  CALL  INH_NR     ( MXL = X_MIN_D, MXH = X_MAX_D, PXH = PMLXH,        &
                     MYL = Y_MIN_D, MYH = Y_MAX_D, PYH = PMLYH,        &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX ,                &
                     AXT = AXU   (:,L ), AYT = AYU    (:,L ),          &
                     SX  = SX_U  (:,L ), SXP = SX_UP  (:   ),          &
                     SY  = SY_U  (:,L ), SYP = SY_UP  (:   ),          &
                     T = UM (:,:,L ), I1 = 2, K1 = 1, TPML = TPML)

  CALL  INH_NR     ( MXL = X_MIN_D, MXH = X_MAX_D, PXH = PMLXH,        &
                     MYL = Y_MIN_D, MYH = Y_MAX_D, PYH = PMLYH,        &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX ,                &
                     AXT = AXV   (:,L ), AYT = AYV    (:,L ),          &
                     SX  = SX_V  (:,L ), SXP = SX_VP  (:   ),          &
                     SY  = SY_V  (:,L ), SYP = SY_VP  (:   ),          &
                     T = VM (:,:,L ), I1 = 1, K1 = 2, TPML = TPML)

  CALL  INH_NR     ( MXL = X_MIN_D, MXH = X_MAX_D, PXH = PMLXH,        &
                     MYL = Y_MIN_D, MYH = Y_MAX_D, PYH = PMLYH,        &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX ,                &
                     AXT = AXW   (:,L ), AYT = AYW    (:,L ),          &
                     SX  = SX_W  (:,L ), SXP = SX_WP  (:   ),          &
                     SY  = SY_W  (:,L ), SYP = SY_WP  (:   ),          &
                     T = WM (:,:,L ), I1 = 1, K1 = 1, TPML = TPML)
END IF


END SUBROUTINE HPLANE_1_VELOCITY

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

SUBROUTINE  HPLANE_IF_VELOCITY ( L )

  USE PRECISION   , ONLY:                                              &
                          PP
  USE CONTROL_DATA, ONLY:                                              &
                          LPAS, REFINE
  USE GRID_MEDIUM , ONLY:                                              &
                          JMF,                                         &
                          DENUF, DENVF, DENWF
  USE AUXIL       , ONLY:                                              &
                          MXF, MYF, TPMLF, MXF1 , MYF1 ,               &
                          PMLLF, PMLXHF, PMLYHF,                       &
                           X_MINF  , X_MAXF  ,                         &
                           Y_MINF  , Y_MAXF  ,                         &
                           X_MIN_DF, X_MAX_DF,                         &
                           Y_MIN_DF, Y_MAX_DF,                         &
                           Z_MIN_DF, Z_MAX_DF
  USE WAVEFIELD   , ONLY:                                              &
                          UMF, VMF, WMF
  USE STRESSFIELD , ONLY:                                              &
                          TXXF , TYYF , TZZF ,                         &
                          TXYF ,        TXZF ,       TYZF
  USE NONREF_BOUND, ONLY:                                              &
                          AXUF, AYUF, AXVF, AYVF, AXWF, AYWF,          &
                                    SX_UF , SX_UPF , SY_UF , SY_UPF ,  &
                                    SX_VF , SX_VPF , SY_VF , SY_VPF ,  &
                                    SX_WF , SX_WPF , SY_WF , SY_WPF
  USE PML         , ONLY:                                              &
            ODA_XLF,  ODB_XLF,  ODC_XLF,  ODA_XHF,  ODB_XHF,  ODC_XHF, &
            ODD_XLF,                      ODD_XHF,                     &
            ODA_YLF,  ODB_YLF,  ODC_YLF,  ODA_YHF,  ODB_YHF,  ODC_YHF, &
            ODD_YLF,                      ODD_YHF,                     &
           ODA2_XLF, ODB2_XLF, ODC2_XLF, ODA2_XHF, ODB2_XHF, ODC2_XHF, &
           ODD2_XLF,                     ODD2_XHF,                     &
           ODA2_YLF, ODB2_YLF, ODC2_YLF, ODA2_YHF, ODB2_YHF, ODC2_YHF, &
           ODD2_YLF,                     ODD2_YHF,                     &
                                            RXXLF , RXXHF ,            &
                                            RXYLF , RXYHF ,            &
                                            RYYLF , RYYHF ,            &
                                            RYXLF , RYXHF ,            &
                                            RZXLF , RZXHF ,            &
                                            RZYLF , RZYHF
  USE INTERFACES  , ONLY:                                              &
                          INH_V_UV_INT  , INH_V_W_INT  ,               &
                          INH_V_U_2ND ,INH_V_V_2ND, INH_V_W_2ND,       &
                          STO_BOR, INH_NR

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: L

!-----------------------------------------------------------------------

  CALL  STO_BOR    ( MXL = X_MIN_DF, MXH = X_MAX_DF,                   &
                     MYL = Y_MIN_DF, MYH = Y_MAX_DF,                   &
                     PXH = PMLXHF, PYH = PMLYHF,                       &
                        X_MIN = X_MINF, X_MAX = X_MAXF,                &
                        Y_MIN = Y_MINF, Y_MAX = Y_MAXF,                &
                     SX_T = SX_UF (:,L ), SY_T = SY_UF (:,L ),         &
                     SX_TP= SX_UPF(:   ), SY_TP= SY_UPF(:   ),         &
                     T = UMF(:,:,L ), I1 = 2, K1 = 1, TPML = TPMLF)

  CALL  STO_BOR    ( MXL = X_MIN_DF, MXH = X_MAX_DF,                   &
                     MYL = Y_MIN_DF, MYH = Y_MAX_DF,                   &
                     PXH = PMLXHF, PYH = PMLYHF,                       &
                        X_MIN = X_MINF, X_MAX = X_MAXF,                &
                        Y_MIN = Y_MINF, Y_MAX = Y_MAXF,                &
                     SX_T = SX_VF (:,L ), SY_T = SY_VF (:,L ),         &
                     SX_TP= SX_VPF(:   ), SY_TP= SY_VPF(:   ),         &
                     T = VMF(:,:,L ), I1 = 1, K1 = 2, TPML = TPMLF)

  CALL  STO_BOR    ( MXL = X_MIN_DF, MXH = X_MAX_DF,                   &
                     MYL = Y_MIN_DF, MYH = Y_MAX_DF,                   &
                     PXH = PMLXHF, PYH = PMLYHF,                       &
                        X_MIN = X_MINF, X_MAX = X_MAXF,                &
                        Y_MIN = Y_MINF, Y_MAX = Y_MAXF,                &
                     SX_T = SX_WF (:,L ), SY_T = SY_WF (:,L ),         &
                     SX_TP= SX_WPF(:   ), SY_TP= SY_WPF(:   ),         &
                     T = WMF(:,:,L ), I1 = 1, K1 = 1, TPML = TPMLF)


  CALL  INH_V_UV_INT                                                   &
              ( L = L, MXL = X_MIN_DF, MXH = X_MAX_DF, MXT = MXF ,     &
                       MYL = Y_MIN_DF, MYH = Y_MAX_DF, MYT = MYF ,     &
                       MZL = Z_MIN_DF, MZH = Z_MAX_DF,                 &
                        X_MIN = X_MINF, X_MAX = X_MAXF,                &
                        Y_MIN = Y_MINF, Y_MAX = Y_MAXF,                &
                       JMT = JMF, DENUT = DENUF, DENVT = DENVF,        &
                       UT = UMF, VT = VMF,                             &
                       TXXT = TXXF, TXYT = TXYF, TXZT = TXZF,          &
                       TYYT = TYYF,              TYZT = TYZF )

  CALL  INH_V_W_INT                                                    &
              ( L = L, MXL = X_MIN_DF, MXH = X_MAX_DF, MXT = MXF,      &
                       MYL = Y_MIN_DF, MYH = Y_MAX_DF, MYT = MYF,      &
                       MZL = Z_MIN_DF, MZH = Z_MAX_DF,                 &
                        X_MIN = X_MINF, X_MAX = X_MAXF,                &
                        Y_MIN = Y_MINF, Y_MAX = Y_MAXF,                &
                       JMT = JMF, DENWT = DENWF,                       &
                       WT = WMF,                                       &
                       TZZT = TZZF, TXZT = TXZF, TYZT = TYZF )

  CALL  INH_V_U_2ND                                                    &
              ( L = L, MXL = X_MIN_DF, MXH = X_MAX_DF, PL  = PMLLF,    &
                       MXT = MXF, MX1 = MXF1,          PXH = PMLXHF,   &
                       MYL = Y_MIN_DF, MYH = Y_MAX_DF,                 &
                       MYT = MYF, MY1 = MYF1,          PYH = PMLYHF,   &
                       MZL = Z_MIN_DF, MZH = Z_MAX_DF,  TPML = TPMLF,  &
                        X_MIN = X_MINF, X_MAX = X_MAXF,                &
                        Y_MIN = Y_MINF, Y_MAX = Y_MAXF,                &
                       JMT = JMF, DENUT = DENUF,                       &
                       UT = UMF,                                       &
                       TXXT = TXXF,              TXYT = TXYF,          &
                       TXZT = TXZF,                                    &
                       RXXLT = RXXLF, RXXHT = RXXHF,                   &
                       RXYLT = RXYLF, RXYHT = RXYHF,                   &
           ODA2_XL = ODA2_XLF, ODB2_XL = ODB2_XLF,                     &
           ODC2_XL = ODC2_XLF, ODD2_XL = ODD2_XLF,                     &
           ODA2_XH = ODA2_XHF, ODB2_XH = ODB2_XHF,                     &
           ODC2_XH = ODC2_XHF, ODD2_XH = ODD2_XHF,                     &
            ODA_YL = ODA_YLF ,  ODB_YL = ODB_YLF ,                     &
            ODC_YL = ODC_YLF ,  ODD_YL = ODD_YLF ,                     &
            ODA_YH = ODA_YHF ,  ODB_YH = ODB_YHF ,                     &
            ODC_YH = ODC_YHF ,  ODD_YH = ODD_YHF )

  CALL  INH_V_V_2ND                                                    &
              ( L = L, MXL = X_MIN_DF, MXH = X_MAX_DF, PL  = PMLLF,    &
                       MXT = MXF, MX1 = MXF1,          PXH = PMLXHF,   &
                       MYL = Y_MIN_DF, MYH = Y_MAX_DF,                 &
                       MYT = MYF, MY1 = MYF1,          PYH = PMLYHF,   &
                       MZL = Z_MIN_DF, MZH = Z_MAX_DF,  TPML = TPMLF,  &
                        X_MIN = X_MINF, X_MAX = X_MAXF,                &
                        Y_MIN = Y_MINF, Y_MAX = Y_MAXF,                &
                       JMT = JMF,                DENVT = DENVF,        &
                                VT = VMF,                              &
                                    TYYT = TYYF, TXYT = TXYF,          &
                                    TYZT = TYZF,                       &
                       RYYLT = RYYLF, RYYHT = RYYHF,                   &
                       RYXLT = RYXLF, RYXHT = RYXHF,                   &
            ODA_XL = ODA_XLF ,  ODB_XL = ODB_XLF ,                     &
            ODC_XL = ODC_XLF ,  ODD_XL = ODD_XLF ,                     &
            ODA_XH = ODA_XHF ,  ODB_XH = ODB_XHF ,                     &
            ODC_XH = ODC_XHF ,  ODD_XH = ODD_XHF ,                     &
           ODA2_YL = ODA2_YLF, ODB2_YL = ODB2_YLF,                     &
           ODC2_YL = ODC2_YLF, ODD2_YL = ODD2_YLF,                     &
           ODA2_YH = ODA2_YHF, ODB2_YH = ODB2_YHF,                     &
           ODC2_YH = ODC2_YHF, ODD2_YH = ODD2_YHF )

  CALL  INH_V_W_2ND                                                    &
              ( L = L, MXL = X_MIN_DF, MXH = X_MAX_DF, PL  = PMLLF,    &
                       MXT = MXF, MX1 = MXF1,          PXH = PMLXHF,   &
                       MYL = Y_MIN_DF, MYH = Y_MAX_DF,                 &
                       MYT = MYF, MY1 = MYF1,          PYH = PMLYHF,   &
                       MZL = Z_MIN_DF, MZH = Z_MAX_DF,  TPML = TPMLF,  &
                        X_MIN = X_MINF, X_MAX = X_MAXF,                &
                        Y_MIN = Y_MINF, Y_MAX = Y_MAXF,                &
                       JMT = JMF, DENWT = DENWF, WT = WMF,             &
                       TZZT = TZZF, TXZT = TXZF, TYZT = TYZF,          &
                       RZXLT = RZXLF, RZXHT = RZXHF,                   &
                       RZYLT = RZYLF, RZYHT = RZYHF,                   &
            ODA_XL = ODA_XLF ,  ODB_XL = ODB_XLF ,                     &
            ODC_XL = ODC_XLF ,  ODD_XL = ODD_XLF ,                     &
            ODA_XH = ODA_XHF ,  ODB_XH = ODB_XHF ,                     &
            ODC_XH = ODC_XHF ,  ODD_XH = ODD_XHF ,                     &
            ODA_YL = ODA_YLF ,  ODB_YL = ODB_YLF ,                     &
            ODC_YL = ODC_YLF ,  ODD_YL = ODD_YLF ,                     &
            ODA_YH = ODA_YHF ,  ODB_YH = ODB_YHF ,                     &
            ODC_YH = ODC_YHF ,  ODD_YH = ODD_YHF )


  CALL  INH_NR     ( MXL = X_MIN_DF, MXH = X_MAX_DF, PXH = PMLXHF,     &
                     MYL = Y_MIN_DF, MYH = Y_MAX_DF, PYH = PMLYHF,     &
                        X_MIN = X_MINF, X_MAX = X_MAXF,                &
                        Y_MIN = Y_MINF, Y_MAX = Y_MAXF,                &
                     AXT = AXUF  (:,L ), AYT = AYUF   (:,L ),          &
                     SX  = SX_UF (:,L ), SXP = SX_UPF (:   ),          &
                     SY  = SY_UF (:,L ), SYP = SY_UPF (:   ),          &
                     T = UMF(:,:,L ), I1 = 2, K1 = 1, TPML = TPMLF)

  CALL  INH_NR     ( MXL = X_MIN_DF, MXH = X_MAX_DF, PXH = PMLXHF,     &
                     MYL = Y_MIN_DF, MYH = Y_MAX_DF, PYH = PMLYHF,     &
                        X_MIN = X_MINF, X_MAX = X_MAXF,                &
                        Y_MIN = Y_MINF, Y_MAX = Y_MAXF,                &
                     AXT = AXVF  (:,L ), AYT = AYVF   (:,L ),          &
                     SX  = SX_VF (:,L ), SXP = SX_VPF (:   ),          &
                     SY  = SY_VF (:,L ), SYP = SY_VPF (:   ),          &
                     T = VMF(:,:,L ), I1 = 1, K1 = 2, TPML = TPMLF)

  CALL  INH_NR     ( MXL = X_MIN_DF, MXH = X_MAX_DF, PXH = PMLXHF,     &
                     MYL = Y_MIN_DF, MYH = Y_MAX_DF, PYH = PMLYHF,     &
                        X_MIN = X_MINF, X_MAX = X_MAXF,                &
                        Y_MIN = Y_MINF, Y_MAX = Y_MAXF,                &
                     AXT = AXWF  (:,L ), AYT = AYWF   (:,L ),          &
                     SX  = SX_WF (:,L ), SXP = SX_WPF (:   ),          &
                     SY  = SY_WF (:,L ), SYP = SY_WPF (:   ),          &
                     T = WMF(:,:,L ), I1 = 1, K1 = 1, TPML = TPMLF)


END SUBROUTINE HPLANE_IF_VELOCITY

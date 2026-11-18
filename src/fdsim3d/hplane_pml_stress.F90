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

SUBROUTINE  HPLANE_PML_STRESS ( L, LL )

  USE PRECISION   , ONLY:                                              &
                          PP
  USE GRID_MEDIUM , ONLY:                                              &
                          JM
  USE CONTROL_DATA, ONLY:                                              &
                          MX, MY, TPML
  USE AUXIL       , ONLY:                                              &
                          AH , BH , CH ,                               &
                               MX1 ,            MY1 ,                  &
                          PMLL, PMLXH, PMLYH,                          &
                           X_MIN  , X_MAX  ,                           &
                           Y_MIN  , Y_MAX  ,                           &
                           X_MIN_D, X_MAX_D,                           &
                           Y_MIN_D, Y_MAX_D,                           &
                           Z_MIN_D, Z_MAX_D
  USE WAVEFIELD   , ONLY:                                              &
                          UM , VM , WM ,                               &
                          XXX , XYY , XZZ , XXY , XXZ , XYZ
  USE STRESSFIELD , ONLY:                                              &
                          TXX , TYY , TZZ ,                            &
                          TXY ,       TXZ ,       TYZ
  USE STRAINFIELD , ONLY:                                              &
                                        EXX , EYY , EXY ,              &
                                        EXZ , EYZ , EZZ
  USE NONREF_BOUND, ONLY:                                              &
                          XXYI2 , XXYIMX  , XXYK2   , XXYKMY  ,        &
                          XXZK2 , XXZKMY1 , XXZI2   , XXZIMX  ,        &
                          XYZK2 , XYZKMY  , XYZI2   , XYZIMX1 ,        &
                          XXXK2 , XXXI2   , XXXKMY1 , XXXIMX1 ,        &
                          XYYK2 , XYYI2   , XYYKMY1 , XYYIMX1 ,        &
                          XZZK2 , XZZI2   , XZZKMY1 , XZZIMX1
  USE PML         , ONLY:                                              &
            ODA_XL ,  ODB_XL ,  ODC_XL ,  ODA_XH ,  ODB_XH ,  ODC_XH , &
            ODD_XL ,                      ODD_XH ,                     &
            ODA_YL ,  ODB_YL ,  ODC_YL ,  ODA_YH ,  ODB_YH ,  ODC_YH , &
            ODD_YL ,                      ODD_YH ,                     &
                                          ODA_ZH ,  ODB_ZH ,  ODC_ZH , &
                                          ODD_ZH ,                     &
           ODA2_XL , ODB2_XL , ODC2_XL , ODA2_XH , ODB2_XH , ODC2_XH , &
           ODD2_XL ,                     ODD2_XH ,                     &
           ODA2_YL , ODB2_YL , ODC2_YL , ODA2_YH , ODB2_YH , ODC2_YH , &
           ODD2_YL ,                     ODD2_YH ,                     &
                                         ODA2_ZH , ODB2_ZH , ODC2_ZH , &
                     ODD2_ZH ,                     &
                                          PZXL  , PZXH  ,              &
                                          PZYL  , PZYH  ,              &
                                          PXXL  , PXXH  ,              &
                                          PYYL  , PYYH  ,              &
                                          PYXL  , PYXH  ,              &
                                          PXYL  , PXYH  ,              &
                                                  PXZH ,               &
                                                  PYZH ,               &
                                                  PZZH
  USE INTERFACES  , ONLY:                                              &
                          STRAIN_XZ_2ND_PML, STRAIN_YZ_2ND_PML,        &
                          STRAIN_ZZ_2ND_PML,                           &
                          STRAIN_XX, STRAIN_YY, STRAIN_XY,             &
                          STRESS_XZ, STRESS_YZ, STRESS_XY,             &
                          STRESS_NN

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: L, LL

!-----------------------------------------------------------------------

  CALL STRAIN_XZ_2ND_PML  ( L = L , LL = LL  ,                         &
       MXL = X_MIN_D, MXH = X_MAX_D, PL = PMLL, MXT = MX, PXH = PMLXH, &
       MYL = Y_MIN_D, MYH = Y_MAX_D,            MYT = MY,              &
       MZL = Z_MIN_D, MZH = Z_MAX_D,                                   &
                        TPML = TPML,                   C = CH ,        &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX ,                &
                        UTM = UM , WTM = WM , XXZT = XXZ , EXZT = EXZ ,&
                        PZXLT = PZXL , PZXHT = PZXH ,                  &
                                       PXZHT = PXZH ,                  &
       ODA2_XL = ODA2_XL , ODB2_XL = ODB2_XL ,                         &
       ODC2_XL = ODC2_XL , ODD2_XL = ODD2_XL ,                         &
       ODA2_XH = ODA2_XH , ODB2_XH = ODB2_XH ,                         &
       ODC2_XH = ODC2_XH , ODD2_XH = ODD2_XH ,                         &
       ODA2_ZH = ODA2_ZH , ODB2_ZH = ODB2_ZH ,                         &
     ODC2_ZH = ODC2_ZH , ODD2_ZH = ODD2_ZH )

  CALL STRAIN_YZ_2ND_PML  ( L = L , LL = LL  ,                         &
       MXL = X_MIN_D, MXH = X_MAX_D,            MXT = MX,              &
       MYL = Y_MIN_D, MYH = Y_MAX_D, PL = PMLL, MYT = MY, PYH = PMLYH, &
       MZL = Z_MIN_D, MZH = Z_MAX_D,                                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX ,                &
                        TPML = TPML,                   C = CH ,        &
                        VTM = VM , WTM = WM , XYZT = XYZ , EYZT = EYZ ,&
                        PZYLT = PZYL , PZYHT = PZYH ,                  &
                                       PYZHT = PYZH ,                  &
       ODA2_YL = ODA2_YL , ODB2_YL = ODB2_YL ,                         &
       ODC2_YL = ODC2_YL , ODD2_YL = ODD2_YL ,                         &
       ODA2_YH = ODA2_YH , ODB2_YH = ODB2_YH ,                         &
       ODC2_YH = ODC2_YH , ODD2_YH = ODD2_YH ,                         &
       ODA2_ZH = ODA2_ZH , ODB2_ZH = ODB2_ZH ,                         &
     ODC2_ZH = ODC2_ZH , ODD2_ZH = ODD2_ZH )

  CALL STRAIN_ZZ_2ND_PML  ( L = L , LL = LL  ,                         &
       MXL = X_MIN_D, MXH = X_MAX_D,            MXT = MX,              &
       MYL = Y_MIN_D, MYH = Y_MAX_D,            MYT = MY,              &
       MZL = Z_MIN_D, MZH = Z_MAX_D,                                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX ,                &
                        TPML = TPML,                   C = CH ,        &
                        WTM = WM ,            XZZT = XZZ , EZZT = EZZ ,&
                                       PZZHT = PZZH ,                  &
                ODA_ZH = ODA_ZH, ODB_ZH = ODB_ZH,                      &
                ODC_ZH = ODC_ZH, ODD_ZH = ODD_ZH )

  CALL STRAIN_XX(L = L ,                                               &
       MXL = X_MIN_D, MXH = X_MAX_D, PL = PMLL, MXT = MX, PXH = PMLXH, &
       MYL = Y_MIN_D, MYH = Y_MAX_D,            MYT = MY,              &
       MZL = Z_MIN_D, MZH = Z_MAX_D,                                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX ,                &
       A = AH, B = BH, C = CH, TPML = TPML,                            &
                        UTM = UM ,            XXXT = XXX , EXXT = EXX ,&
                        PXXLT = PXXL , PXXHT = PXXH ,                  &
             ODA_XL = ODA_XL , ODB_XL = ODB_XL ,                       &
             ODC_XL = ODC_XL , ODD_XL = ODD_XL ,                       &
             ODA_XH = ODA_XH , ODB_XH = ODB_XH ,                       &
             ODC_XH = ODC_XH , ODD_XH = ODD_XH )

  CALL STRAIN_YY(L = L ,                                               &
       MXL = X_MIN_D, MXH = X_MAX_D,            MXT = MX,              &
       MYL = Y_MIN_D, MYH = Y_MAX_D, PL = PMLL, MYT = MY, PYH = PMLYH, &
       MZL = Z_MIN_D, MZH = Z_MAX_D,                                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX ,                &
       A = AH, B = BH, C = CH, TPML = TPML,                            &
                        VTM = VM ,            XYYT = XYY , EYYT = EYY ,&
                        PYYLT = PYYL , PYYHT = PYYH ,                  &
             ODA_YL = ODA_YL , ODB_YL = ODB_YL ,                       &
             ODC_YL = ODC_YL , ODD_YL = ODD_YL ,                       &
             ODA_YH = ODA_YH , ODB_YH = ODB_YH ,                       &
             ODC_YH = ODC_YH , ODD_YH = ODD_YH )

  CALL STRAIN_XY(L = L ,                                               &
       MXL = X_MIN_D, MXH = X_MAX_D,            MXT = MX, PXH = PMLXH, &
       MYL = Y_MIN_D, MYH = Y_MAX_D, PL = PMLL, MYT = MY, PYH = PMLYH, &
       MZL = Z_MIN_D, MZH = Z_MAX_D,                                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX ,                &
       A = AH, B = BH, C = CH, TPML = TPML,                            &
                        UTM = UM , VTM = VM , XXYT = XXY , EXYT = EXY ,&
                        PYXLT = PYXL , PYXHT = PYXH ,                  &
                        PXYLT = PXYL , PXYHT = PXYH ,                  &
       ODA2_XL = ODA2_XL , ODB2_XL = ODB2_XL ,                         &
       ODC2_XL = ODC2_XL , ODD2_XL = ODD2_XL ,                         &
       ODA2_XH = ODA2_XH , ODB2_XH = ODB2_XH ,                         &
       ODC2_XH = ODC2_XH , ODD2_XH = ODD2_XH ,                         &
       ODA2_YL = ODA2_YL , ODB2_YL = ODB2_YL ,                         &
       ODC2_YL = ODC2_YL , ODD2_YL = ODD2_YL ,                         &
       ODA2_YH = ODA2_YH , ODB2_YH = ODB2_YH ,                         &
       ODC2_YH = ODC2_YH , ODD2_YH = ODD2_YH )

  CALL STRESS_XY (L=L ,MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,        &
                       MYL = Y_MIN_D, MYH = Y_MAX_D, MYT = MY ,        &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX ,                &
                       TXYT = TXY , EXYT = EXY , XXYT = XXY ,          &
                       JMT = JM , TPML = TPML,                         &
                       XXYK2T = XXYK2 , XXYKMYT  = XXYKMY  ,           &
                       XXYI2T = XXYI2 , XXYIMXT  = XXYIMX  )

  CALL STRESS_XZ (L=L ,MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,        &
                       MYL = Y_MIN_D, MYH = Y_MAX_D, MYT = MY ,        &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX ,                &
                       TXZT = TXZ , EXZT = EXZ , XXZT = XXZ ,          &
                       JMT = JM , TPML = TPML,                         &
                       XXZK2T = XXZK2 , XXZKMY1T = XXZKMY1 ,           &
                       XXZI2T = XXZI2 , XXZIMXT  = XXZIMX  )

  CALL STRESS_YZ (L=L ,MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,        &
                       MYL = Y_MIN_D, MYH = Y_MAX_D, MYT = MY ,        &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX ,                &
                       TYZT = TYZ , EYZT = EYZ , XYZT = XYZ ,          &
                       JMT = JM , TPML = TPML,                         &
                       XYZK2T = XYZK2 , XYZKMYT  = XYZKMY  ,           &
                       XYZI2T = XYZI2 , XYZIMX1T = XYZIMX1  )


  CALL STRESS_NN (L=L ,MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,        &
                       MYL = Y_MIN_D, MYH = Y_MAX_D, MYT = MY ,        &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        Y_MIN = Y_MIN , Y_MAX = Y_MAX ,                &
                       TXXT = TXX , TYYT = TYY , TZZT = TZZ ,          &
                       EXXT = EXX , EYYT = EYY , EZZT = EZZ ,          &
                       XXXT = XXX , XYYT = XYY , XZZT = XZZ ,          &
                       JMT = JM , TPML = TPML,                         &
                       XXXK2T   = XXXK2   , XXXI2T   = XXXI2   ,       &
                       XXXKMY1T = XXXKMY1 , XXXIMX1T = XXXIMX1 ,       &
                       XYYK2T   = XYYK2   , XYYI2T   = XYYI2   ,       &
                       XYYKMY1T = XYYKMY1 , XYYIMX1T = XYYIMX1 ,       &
                       XZZK2T   = XZZK2   , XZZI2T   = XZZI2   ,       &
                       XZZKMY1T = XZZKMY1 , XZZIMX1T = XZZIMX1  )


END SUBROUTINE HPLANE_PML_STRESS

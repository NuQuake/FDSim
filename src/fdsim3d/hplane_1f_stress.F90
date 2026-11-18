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

SUBROUTINE  HPLANE_1F_STRESS

  USE PRECISION   , ONLY:                                              &
                          PP
  USE GRID_MEDIUM , ONLY:                                              &
                          JMF
  USE CONTROL_DATA, ONLY:                                              &
                          STRESS_IMAGING
  USE AUXIL       , ONLY:                                              &
                          MXF, MYF, TPMLF,                             &
                          AHF, BHF, CHF,                               &
                          PMLLF, PMLXHF, PMLYHF,                       &
                           X_MINF  , X_MAXF  ,                         &
                           Y_MINF  , Y_MAXF  ,                         &
                           X_MIN_DF, X_MAX_DF,                         &
                           Y_MIN_DF, Y_MAX_DF,                         &
                           Z_MIN_DF, Z_MAX_DF
  USE WAVEFIELD   , ONLY:                                              &
                          UMF, VMF, WMF,                               &
                          XXXF, XYYF, XZZF, XXYF, XXZF, XYZF
  USE STRESSFIELD , ONLY:                                              &
                          TXXF, TYYF, TZZF,                            &
                          TXYF,       TXZF,       TYZF
  USE STRAINFIELD , ONLY:                                              &
                                        EXXF, EYYF, EXYF,              &
                                        EXZF, EYZF, EZZF
  USE NONREF_BOUND, ONLY:                                              &
                          XXYI2F, XXYIMXF , XXYK2F  , XXYKMYF ,        &
                          XXZK2F, XXZKMY1F, XXZI2F  , XXZIMXF ,        &
                          XYZK2F, XYZKMYF , XYZI2F  , XYZIMX1F,        &
                          XXXK2F, XXXI2F  , XXXKMY1F, XXXIMX1F,        &
                          XYYK2F, XYYI2F  , XYYKMY1F, XYYIMX1F,        &
                          XZZK2F, XZZI2F  , XZZKMY1F, XZZIMX1F
  USE PML         , ONLY:                                              &
            ODA_XLF,  ODB_XLF,  ODC_XLF,  ODA_XHF,  ODB_XHF,  ODC_XHF, &
            ODD_XLF,                      ODD_XHF,                     &
            ODA_YLF,  ODB_YLF,  ODC_YLF,  ODA_YHF,  ODB_YHF,  ODC_YHF, &
            ODD_YLF,                      ODD_YHF,                     &
           ODA2_XLF, ODB2_XLF, ODC2_XLF, ODA2_XHF, ODB2_XHF, ODC2_XHF, &
           ODD2_XLF,                     ODD2_XHF,                     &
           ODA2_YLF, ODB2_YLF, ODC2_YLF, ODA2_YHF, ODB2_YHF, ODC2_YHF, &
           ODD2_YLF,                     ODD2_YHF,                     &
                                          PZXLF , PZXHF ,              &
                                          PZYLF , PZYHF ,              &
                                          PXXLF , PXXHF ,              &
                                          PYYLF , PYYHF ,              &
                                          PYXLF , PYXHF ,              &
                                          PXYLF , PXYHF
  USE INTERFACES  , ONLY:                                              &
                          STRAIN_XX, STRAIN_YY, STRAIN_ZZ,             &
                                                STRAIN_XY,             &
                          STRESS_XZ, STRESS_YZ, STRESS_XY,             &
                          STRESS_NN,                                   &
                          STRAIN_XZ_1, STRAIN_YZ_1

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: I,K

!-----------------------------------------------------------------------

IF (STRESS_IMAGING) THEN
  TXZF (:,:,1) = -TXZF (:,:,3)
  TYZF (:,:,1) = -TYZF (:,:,3)
  TZZF (:,:,1) = -TZZF (:,:,2)
  DO  K = MAX(3,Y_MINF), MIN(MYF-2,Y_MAXF)
    DO  I = MAX(3,X_MINF), MIN(MXF-2,X_MAXF)
      TZZF (I,K,1) =  TZZF (I,K,1)-TZZF (I,K,2)+1./3.*TZZF (I,K,3)
    END DO
  END DO

  TXZF (:,:,2) = 0.
  TYZF (:,:,2) = 0.
ELSE

  CALL STRAIN_XZ_1                                                     &
 ( MXL = X_MIN_DF, MXH = X_MAX_DF, PL = PMLLF, MXT = MXF, PXH = PMLXHF,&
   MYL = Y_MIN_DF, MYH = Y_MAX_DF,             MYT = MYF,              &
   MZL = Z_MIN_DF, MZH = Z_MAX_DF,                                     &
                        X_MIN = X_MINF, X_MAX = X_MAXF,                &
                        Y_MIN = Y_MINF, Y_MAX = Y_MAXF,                &
       A = AHF, B = BHF, C = CHF, TPML = TPMLF,                        &
                        UTM = UMF, WTM = WMF, XXZT = XXZF, EXZT = EXZF,&
                        PZXLT = PZXLF, PZXHT = PZXHF,                  &
       ODA2_XL = ODA2_XLF, ODB2_XL = ODB2_XLF,                         &
       ODC2_XL = ODC2_XLF, ODD2_XL = ODD2_XLF,                         &
       ODA2_XH = ODA2_XHF, ODB2_XH = ODB2_XHF,                         &
       ODC2_XH = ODC2_XHF, ODD2_XH = ODD2_XHF )

  CALL STRAIN_YZ_1                                                     &
 ( MXL = X_MIN_DF, MXH = X_MAX_DF,             MXT = MXF,              &
   MYL = Y_MIN_DF, MYH = Y_MAX_DF, PL = PMLLF, MYT = MYF, PYH = PMLYHF,&
   MZL = Z_MIN_DF, MZH = Z_MAX_DF,                                     &
                        X_MIN = X_MINF, X_MAX = X_MAXF,                &
                        Y_MIN = Y_MINF, Y_MAX = Y_MAXF,                &
       A = AHF, B = BHF, C = CHF, TPML = TPMLF,                        &
                        VTM = VMF, WTM = WMF, XYZT = XYZF, EYZT = EYZF,&
                        PZYLT = PZYLF, PZYHT = PZYHF,                  &
       ODA2_YL = ODA2_YLF, ODB2_YL = ODB2_YLF,                         &
       ODC2_YL = ODC2_YLF, ODD2_YL = ODD2_YLF,                         &
       ODA2_YH = ODA2_YHF, ODB2_YH = ODB2_YHF,                         &
       ODC2_YH = ODC2_YHF, ODD2_YH = ODD2_YHF )

  CALL STRAIN_ZZ(L = 1, MXL = X_MIN_DF, MXH = X_MAX_DF, MXT = MXF ,    &
                        MYL = Y_MIN_DF, MYH = Y_MAX_DF, MYT = MYF ,    &
                        MZL = Z_MIN_DF, MZH = Z_MAX_DF,                &
                        X_MIN = X_MINF, X_MAX = X_MAXF,                &
                        Y_MIN = Y_MINF, Y_MAX = Y_MAXF,                &
                        TPML = TPMLF, A = AHF, B = BHF, C = CHF,       &
                        WTM = WMF,            XZZT = XZZF, EZZT = EZZF )

  CALL STRAIN_XX(L = 1 ,                                               &
   MXL = X_MIN_DF, MXH = X_MAX_DF, PL = PMLLF, MXT = MXF, PXH = PMLXHF,&
   MYL = Y_MIN_DF, MYH = Y_MAX_DF,             MYT = MYF,              &
   MZL = Z_MIN_DF, MZH = Z_MAX_DF,                                     &
                        X_MIN = X_MINF, X_MAX = X_MAXF,                &
                        Y_MIN = Y_MINF, Y_MAX = Y_MAXF,                &
       A = AHF, B = BHF, C = CHF, TPML = TPMLF,                        &
                        UTM = UMF,            XXXT = XXXF, EXXT = EXXF,&
                        PXXLT = PXXLF, PXXHT = PXXHF,                  &
             ODA_XL = ODA_XLF, ODB_XL = ODB_XLF,                       &
             ODC_XL = ODC_XLF, ODD_XL = ODD_XLF,                       &
             ODA_XH = ODA_XHF, ODB_XH = ODB_XHF,                       &
             ODC_XH = ODC_XHF, ODD_XH = ODD_XHF )

  CALL STRAIN_YY(L = 1 ,                                               &
   MXL = X_MIN_DF, MXH = X_MAX_DF,             MXT = MXF,              &
   MYL = Y_MIN_DF, MYH = Y_MAX_DF, PL = PMLLF, MYT = MYF, PYH = PMLYHF,&
   MZL = Z_MIN_DF, MZH = Z_MAX_DF,                                     &
                        X_MIN = X_MINF, X_MAX = X_MAXF,                &
                        Y_MIN = Y_MINF, Y_MAX = Y_MAXF,                &
       A = AHF, B = BHF, C = CHF, TPML = TPMLF,                        &
                        VTM = VMF,            XYYT = XYYF, EYYT = EYYF,&
                        PYYLT = PYYLF, PYYHT = PYYHF,                  &
             ODA_YL = ODA_YLF, ODB_YL = ODB_YLF,                       &
             ODC_YL = ODC_YLF, ODD_YL = ODD_YLF,                       &
             ODA_YH = ODA_YHF, ODB_YH = ODB_YHF,                       &
             ODC_YH = ODC_YHF, ODD_YH = ODD_YHF )

  CALL STRAIN_XY(L = 1 ,                                               &
   MXL = X_MIN_DF, MXH = X_MAX_DF,             MXT = MXF, PXH = PMLXHF,&
   MYL = Y_MIN_DF, MYH = Y_MAX_DF, PL = PMLLF, MYT = MYF, PYH = PMLYHF,&
   MZL = Z_MIN_DF, MZH = Z_MAX_DF,                                     &
                        X_MIN = X_MINF, X_MAX = X_MAXF,                &
                        Y_MIN = Y_MINF, Y_MAX = Y_MAXF,                &
       A = AHF, B = BHF, C = CHF, TPML = TPMLF,                        &
                        UTM = UMF, VTM = VMF, XXYT = XXYF, EXYT = EXYF,&
                        PYXLT = PYXLF, PYXHT = PYXHF,                  &
                        PXYLT = PXYLF, PXYHT = PXYHF,                  &
       ODA2_XL = ODA2_XLF, ODB2_XL = ODB2_XLF,                         &
       ODC2_XL = ODC2_XLF, ODD2_XL = ODD2_XLF,                         &
       ODA2_XH = ODA2_XHF, ODB2_XH = ODB2_XHF,                         &
       ODC2_XH = ODC2_XHF, ODD2_XH = ODD2_XHF,                         &
       ODA2_YL = ODA2_YLF, ODB2_YL = ODB2_YLF,                         &
       ODC2_YL = ODC2_YLF, ODD2_YL = ODD2_YLF,                         &
       ODA2_YH = ODA2_YHF, ODB2_YH = ODB2_YHF,                         &
       ODC2_YH = ODC2_YHF, ODD2_YH = ODD2_YHF )

  CALL STRESS_XY (L=1, MXL = X_MIN_DF, MXH = X_MAX_DF, MXT = MXF ,     &
                       MYL = Y_MIN_DF, MYH = Y_MAX_DF, MYT = MYF ,     &
                       MZL = Z_MIN_DF, MZH = Z_MAX_DF,                 &
                        X_MIN = X_MINF, X_MAX = X_MAXF,                &
                        Y_MIN = Y_MINF, Y_MAX = Y_MAXF,                &
                       TXYT = TXYF, EXYT = EXYF, XXYT = XXYF,          &
                       JMT = JMF , TPML = TPMLF,                       &
                       XXYK2T = XXYK2F, XXYKMYT  = XXYKMYF ,           &
                       XXYI2T = XXYI2F, XXYIMXT  = XXYIMXF )

  CALL STRESS_XZ (L=1 ,MXL = X_MIN_DF, MXH = X_MAX_DF, MXT = MXF,      &
                       MYL = Y_MIN_DF, MYH = Y_MAX_DF, MYT = MYF,      &
                       MZL = Z_MIN_DF, MZH = Z_MAX_DF,                 &
                        X_MIN = X_MINF, X_MAX = X_MAXF,                &
                        Y_MIN = Y_MINF, Y_MAX = Y_MAXF,                &
                       TXZT = TXZF, EXZT = EXZF, XXZT = XXZF,          &
                       JMT = JMF, TPML = TPMLF,                        &
                       XXZK2T = XXZK2F, XXZKMY1T = XXZKMY1F,           &
                       XXZI2T = XXZI2F, XXZIMXT  = XXZIMXF )

  CALL STRESS_YZ (L=1 ,MXL = X_MIN_DF, MXH = X_MAX_DF, MXT = MXF ,     &
                       MYL = Y_MIN_DF, MYH = Y_MAX_DF, MYT = MYF ,     &
                       MZL = Z_MIN_DF, MZH = Z_MAX_DF,                 &
                        X_MIN = X_MINF, X_MAX = X_MAXF,                &
                        Y_MIN = Y_MINF, Y_MAX = Y_MAXF,                &
                       TYZT = TYZF, EYZT = EYZF, XYZT = XYZF,          &
                       JMT = JMF, TPML = TPMLF,                        &
                       XYZK2T = XYZK2F, XYZKMYT  = XYZKMYF ,           &
                       XYZI2T = XYZI2F, XYZIMX1T = XYZIMX1F )


  CALL STRESS_NN (L=1 ,MXL = X_MIN_DF, MXH = X_MAX_DF, MXT = MXF ,     &
                       MYL = Y_MIN_DF, MYH = Y_MAX_DF, MYT = MYF ,     &
                       MZL = Z_MIN_DF, MZH = Z_MAX_DF,                 &
                        X_MIN = X_MINF, X_MAX = X_MAXF,                &
                        Y_MIN = Y_MINF, Y_MAX = Y_MAXF,                &
                       TXXT = TXXF, TYYT = TYYF, TZZT = TZZF,          &
                       EXXT = EXXF, EYYT = EYYF, EZZT = EZZF,          &
                       XXXT = XXXF, XYYT = XYYF, XZZT = XZZF,          &
                       JMT = JMF, TPML = TPMLF,                        &
                       XXXK2T   = XXXK2F  , XXXI2T   = XXXI2F  ,       &
                       XXXKMY1T = XXXKMY1F, XXXIMX1T = XXXIMX1F,       &
                       XYYK2T   = XYYK2F  , XYYI2T   = XYYI2F  ,       &
                       XYYKMY1T = XYYKMY1F, XYYIMX1T = XYYIMX1F,       &
                       XZZK2T   = XZZK2F  , XZZI2T   = XZZI2F  ,       &
                       XZZKMY1T = XZZKMY1F, XZZIMX1T = XZZIMX1F )
END IF

END SUBROUTINE HPLANE_1F_STRESS

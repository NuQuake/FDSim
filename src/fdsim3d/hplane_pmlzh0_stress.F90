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

SUBROUTINE  HPLANE_PMLZH0_STRESS

  USE PRECISION   , ONLY:                                              &
                          PP
  USE GRID_MEDIUM , ONLY:                                              &
                          JM
  USE CONTROL_DATA, ONLY:                                              &
                          MX, MY, MZ, TPML
  USE AUXIL       , ONLY:                                              &
                          AH , BH , CH ,                               &
                          PMLL, PMLXH, PMLYH,                          &
                           X_MIN  , X_MAX  ,                           &
                           Y_MIN  , Y_MAX  ,                           &
                           X_MIN_D, X_MAX_D,                           &
                           Y_MIN_D, Y_MAX_D,                           &
                           Z_MIN_D, Z_MAX_D
  USE WAVEFIELD   , ONLY:                                              &
                          UM , VM , WM ,                               &
                          XXZ , XYZ
  USE STRESSFIELD , ONLY:                                              &
                                TXZ ,       TYZ
  USE STRAINFIELD , ONLY:                                              &
                                        EXZ , EYZ
  USE NONREF_BOUND, ONLY:                                              &
                          XXZK2 , XXZKMY1 , XXZI2   , XXZIMX  ,        &
                          XYZK2 , XYZKMY  , XYZI2   , XYZIMX1
  USE PML         , ONLY:                                              &
           ODA2_XL , ODB2_XL , ODC2_XL , ODA2_XH , ODB2_XH , ODC2_XH , &
           ODD2_XL ,                     ODD2_XH ,                     &
           ODA2_YL , ODB2_YL , ODC2_YL , ODA2_YH , ODB2_YH , ODC2_YH , &
           ODD2_YL ,                     ODD2_YH ,                     &
                                         ODA2_ZH , ODB2_ZH , ODC2_ZH , &
                     ODD2_ZH ,                     &
                                          PZXL  , PZXH  ,              &
                                          PZYL  , PZYH  ,              &
                                                  PXZH ,               &
                                                  PYZH
  USE INTERFACES  , ONLY:                                              &
                          STRAIN_XZ_2ND_PML, STRAIN_YZ_2ND_PML,        &
                          STRESS_XZ    , STRESS_YZ

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER             :: L, LL

!-----------------------------------------------------------------------

  L  = MZ+TPML
  LL = TPML+1

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


END SUBROUTINE HPLANE_PMLZH0_STRESS

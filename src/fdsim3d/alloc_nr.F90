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

SUBROUTINE ALLOC_NR

  USE PRECISION   , ONLY:                                              &
                          PP
  USE CONTROL_DATA, ONLY:                                              &
                          MX, MY, MZ, LPAS

  USE NONREF_BOUND, ONLY:                                              &
                          NR_DATA_X, NR_DATA_Y, NR_DATA_Z,             &
                                AXU, AXUF, AYU, AYUF,                  &
                                AXV, AXVF, AYV, AYVF,                  &
                                AXW, AXWF, AYW, AYWF,                  &
                                         SX_U , SX_UP ,                &
                                         SX_V , SX_VP ,                &
                                         SX_W , SX_WP ,                &
                                         SY_U , SY_UP ,                &
                                         SY_V , SY_VP ,                &
                                         SY_W , SY_WP ,                &
                                         SX_UF, SX_UPF,                &
                                         SX_VF, SX_VPF,                &
                                         SX_WF, SX_WPF,                &
                                         SY_UF, SY_UPF,                &
                                         SY_VF, SY_VPF,                &
                                         SY_WF, SY_WPF,                &
                                   SZ_U, SZ_UP,                        &
                                   SZ_V, SZ_VP,                        &
                                   SZ_W, SZ_WP,                        &
                          XXYI2 , XXYIMX  , XXYK2   , XXYKMY  ,        &
                          XXZK2 , XXZKMY1 , XXZI2   , XXZIMX  ,        &
                          XYZK2 , XYZKMY  , XYZI2   , XYZIMX1 ,        &
                          XXXK2 , XXXI2   , XXXKMY1 , XXXIMX1 ,        &
                          XYYK2 , XYYI2   , XYYKMY1 , XYYIMX1 ,        &
                          XZZK2 , XZZI2   , XZZKMY1 , XZZIMX1 ,        &
                          XXYI2F, XXYIMXF , XXYK2F  , XXYKMYF ,        &
                          XXZK2F, XXZKMY1F, XXZI2F  , XXZIMXF ,        &
                          XYZK2F, XYZKMYF , XYZI2F  , XYZIMX1F,        &
                          XXXK2F, XXXI2F  , XXXKMY1F, XXXIMX1F,        &
                          XYYK2F, XYYI2F  , XYYKMY1F, XYYIMX1F,        &
                          XZZK2F, XZZI2F  , XZZKMY1F, XZZIMX1F

  USE AUXIL       , ONLY:                                              &
                          X_MIN_D, X_MAX_D, X_MIN_DF, X_MAX_DF,        &
                          Y_MIN_D, Y_MAX_D, Y_MIN_DF, Y_MAX_DF,        &
                          Z_MIN_D, Z_MAX_D, Z_MIN_DF, Z_MAX_DF

!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------

  IF ( LPAS /= 0 ) THEN
    ALLOCATE (  AXUF  (X_MIN_DF:X_MAX_DF , Z_MIN_DF:Z_MAX_DF ),        &
                AXVF  (X_MIN_DF:X_MAX_DF , Z_MIN_DF:Z_MAX_DF ),        &
                AXWF  (X_MIN_DF:X_MAX_DF , Z_MIN_DF:Z_MAX_DF ),        &
                AYUF  (Y_MIN_DF:Y_MAX_DF , Z_MIN_DF:Z_MAX_DF ),        &
                AYVF  (Y_MIN_DF:Y_MAX_DF , Z_MIN_DF:Z_MAX_DF ),        &
                AYWF  (Y_MIN_DF:Y_MAX_DF , Z_MIN_DF:Z_MAX_DF )   )
    ALLOCATE (  SX_UF (X_MIN_DF:X_MAX_DF , Z_MIN_DF:Z_MAX_DF ),        &
                SX_VF (X_MIN_DF:X_MAX_DF , Z_MIN_DF:Z_MAX_DF ),        &
                SX_WF (X_MIN_DF:X_MAX_DF , Z_MIN_DF:Z_MAX_DF ),        &
                SY_UF (Y_MIN_DF:Y_MAX_DF , Z_MIN_DF:Z_MAX_DF ),        &
                SY_VF (Y_MIN_DF:Y_MAX_DF , Z_MIN_DF:Z_MAX_DF ),        &
                SY_WF (Y_MIN_DF:Y_MAX_DF , Z_MIN_DF:Z_MAX_DF ),        &
                SX_UPF(X_MIN_DF:X_MAX_DF),                             &
                SX_VPF(X_MIN_DF:X_MAX_DF),                             &
                SX_WPF(X_MIN_DF:X_MAX_DF),                             &
                SY_UPF(Y_MIN_DF:Y_MAX_DF),                             &
                SY_VPF(Y_MIN_DF:Y_MAX_DF),                             &
                SY_WPF(Y_MIN_DF:Y_MAX_DF) )
    ALLOCATE (  XXYI2F  (Y_MIN_DF:Y_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XXYIMXF (Y_MIN_DF:Y_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XXYK2F  (X_MIN_DF:X_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XXYKMYF (X_MIN_DF:X_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XXZK2F  (X_MIN_DF:X_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XXZKMY1F(X_MIN_DF:X_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XXZI2F  (Y_MIN_DF:Y_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XXZIMXF (Y_MIN_DF:Y_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XYZK2F  (X_MIN_DF:X_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XYZKMYF (X_MIN_DF:X_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XYZI2F  (Y_MIN_DF:Y_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XYZIMX1F(Y_MIN_DF:Y_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XXXK2F  (X_MIN_DF:X_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XXXKMY1F(X_MIN_DF:X_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XXXI2F  (Y_MIN_DF:Y_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XXXIMX1F(Y_MIN_DF:Y_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XYYK2F  (X_MIN_DF:X_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XYYKMY1F(X_MIN_DF:X_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XYYI2F  (Y_MIN_DF:Y_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XYYIMX1F(Y_MIN_DF:Y_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XZZK2F  (X_MIN_DF:X_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XZZKMY1F(X_MIN_DF:X_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XZZI2F  (Y_MIN_DF:Y_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF ),   &
                XZZIMX1F(Y_MIN_DF:Y_MAX_DF , 4, Z_MIN_DF:Z_MAX_DF )   )

    SX_UF = NR_DATA_X ( 0._PP, 0._PP, 0._PP, 0._PP, 0._PP, 0._PP )
    SX_VF = NR_DATA_X ( 0._PP, 0._PP, 0._PP, 0._PP, 0._PP, 0._PP )
    SX_WF = NR_DATA_X ( 0._PP, 0._PP, 0._PP, 0._PP, 0._PP, 0._PP )
    SY_UF = NR_DATA_Y ( 0._PP, 0._PP, 0._PP, 0._PP, 0._PP, 0._PP )
    SY_VF = NR_DATA_Y ( 0._PP, 0._PP, 0._PP, 0._PP, 0._PP, 0._PP )
    SY_WF = NR_DATA_Y ( 0._PP, 0._PP, 0._PP, 0._PP, 0._PP, 0._PP )

    XXYI2F  = 0.
    XXYIMXF = 0.
    XXYK2F  = 0.
    XXYKMYF = 0.
    XXZK2F  = 0.
    XXZKMY1F= 0.
    XXZI2F  = 0.
    XXZIMXF = 0.
    XYZK2F  = 0.
    XYZKMYF = 0.
    XYZI2F  = 0.
    XYZIMX1F= 0.
    XXXK2F  = 0.
    XXXI2F  = 0.
    XXXKMY1F= 0.
    XXXIMX1F= 0.
    XYYK2F  = 0.
    XYYI2F  = 0.
    XYYKMY1F= 0.
    XYYIMX1F= 0.
    XZZK2F  = 0.
    XZZI2F  = 0.
    XZZKMY1F= 0.
    XZZIMX1F= 0.
  END IF

  ALLOCATE (  AXU   (X_MIN_D:X_MAX_D , Z_MIN_D:Z_MAX_D),               &
              AXV   (X_MIN_D:X_MAX_D , Z_MIN_D:Z_MAX_D),               &
              AXW   (X_MIN_D:X_MAX_D , Z_MIN_D:Z_MAX_D),               &
              AYU   (Y_MIN_D:Y_MAX_D , Z_MIN_D:Z_MAX_D),               &
              AYV   (Y_MIN_D:Y_MAX_D , Z_MIN_D:Z_MAX_D),               &
              AYW   (Y_MIN_D:Y_MAX_D , Z_MIN_D:Z_MAX_D) )
  ALLOCATE (  SX_U  (X_MIN_D:X_MAX_D , Z_MIN_D:Z_MAX_D),               &
              SX_V  (X_MIN_D:X_MAX_D , Z_MIN_D:Z_MAX_D),               &
              SX_W  (X_MIN_D:X_MAX_D , Z_MIN_D:Z_MAX_D),               &
              SY_U  (Y_MIN_D:Y_MAX_D , Z_MIN_D:Z_MAX_D),               &
              SY_V  (Y_MIN_D:Y_MAX_D , Z_MIN_D:Z_MAX_D),               &
              SY_W  (Y_MIN_D:Y_MAX_D , Z_MIN_D:Z_MAX_D),               &
              SX_UP (X_MIN_D:X_MAX_D),                                 &
              SX_VP (X_MIN_D:X_MAX_D),                                 &
              SX_WP (X_MIN_D:X_MAX_D),                                 &
              SY_UP (Y_MIN_D:Y_MAX_D),                                 &
              SY_VP (Y_MIN_D:Y_MAX_D),                                 &
              SY_WP (Y_MIN_D:Y_MAX_D),                                 &
              SZ_U  (X_MIN_D:X_MAX_D , Y_MIN_D:Y_MAX_D),               &
              SZ_V  (X_MIN_D:X_MAX_D , Y_MIN_D:Y_MAX_D),               &
              SZ_W  (X_MIN_D:X_MAX_D , Y_MIN_D:Y_MAX_D),               &
              SZ_UP (X_MIN_D:X_MAX_D , Y_MIN_D:Y_MAX_D),               &
              SZ_VP (X_MIN_D:X_MAX_D , Y_MIN_D:Y_MAX_D),               &
              SZ_WP (X_MIN_D:X_MAX_D , Y_MIN_D:Y_MAX_D)  )
  ALLOCATE (  XXYI2   (Y_MIN_D:Y_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XXYIMX  (Y_MIN_D:Y_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XXYK2   (X_MIN_D:X_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XXYKMY  (X_MIN_D:X_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XXZK2   (X_MIN_D:X_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XXZKMY1 (X_MIN_D:X_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XXZI2   (Y_MIN_D:Y_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XXZIMX  (Y_MIN_D:Y_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XYZK2   (X_MIN_D:X_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XYZKMY  (X_MIN_D:X_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XYZI2   (Y_MIN_D:Y_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XYZIMX1 (Y_MIN_D:Y_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XXXK2   (X_MIN_D:X_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XXXKMY1 (X_MIN_D:X_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XXXI2   (Y_MIN_D:Y_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XXXIMX1 (Y_MIN_D:Y_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XYYK2   (X_MIN_D:X_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XYYKMY1 (X_MIN_D:X_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XYYI2   (Y_MIN_D:Y_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XYYIMX1 (Y_MIN_D:Y_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XZZK2   (X_MIN_D:X_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XZZKMY1 (X_MIN_D:X_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XZZI2   (Y_MIN_D:Y_MAX_D , 4, Z_MIN_D:Z_MAX_D ),         &
              XZZIMX1 (Y_MIN_D:Y_MAX_D , 4, Z_MIN_D:Z_MAX_D )   )

  SX_U  = NR_DATA_X ( 0._PP, 0._PP, 0._PP, 0._PP, 0._PP, 0._PP )
  SX_V  = NR_DATA_X ( 0._PP, 0._PP, 0._PP, 0._PP, 0._PP, 0._PP )
  SX_W  = NR_DATA_X ( 0._PP, 0._PP, 0._PP, 0._PP, 0._PP, 0._PP )
  SY_U  = NR_DATA_Y ( 0._PP, 0._PP, 0._PP, 0._PP, 0._PP, 0._PP )
  SY_V  = NR_DATA_Y ( 0._PP, 0._PP, 0._PP, 0._PP, 0._PP, 0._PP )
  SY_W  = NR_DATA_Y ( 0._PP, 0._PP, 0._PP, 0._PP, 0._PP, 0._PP )

  SZ_U  = NR_DATA_Z ( 0._PP, 0._PP, 0._PP )
  SZ_V  = NR_DATA_Z ( 0._PP, 0._PP, 0._PP )
  SZ_W  = NR_DATA_Z ( 0._PP, 0._PP, 0._PP )

  XXYI2   = 0.
  XXYIMX  = 0.
  XXYK2   = 0.
  XXYKMY  = 0.
  XXZK2   = 0.
  XXZKMY1 = 0.
  XXZI2   = 0.
  XXZIMX  = 0.
  XYZK2   = 0.
  XYZKMY  = 0.
  XYZI2   = 0.
  XYZIMX1 = 0.
  XXXK2   = 0.
  XXXI2   = 0.
  XXXKMY1 = 0.
  XXXIMX1 = 0.
  XYYK2   = 0.
  XYYI2   = 0.
  XYYKMY1 = 0.
  XYYIMX1 = 0.
  XZZK2   = 0.
  XZZI2   = 0.
  XZZKMY1 = 0.
  XZZIMX1 = 0.

END SUBROUTINE ALLOC_NR

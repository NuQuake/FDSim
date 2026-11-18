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

MODULE NONREF_BOUND

  USE PRECISION, ONLY: PP

!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------

  TYPE NR_PAR_A_Y
    REAL (KIND=PP) ::                                                  &
                    ALE01,   ALE02,   ALE10,   ALE11,                  &
                    ALE12,   ALE20,   ALE21,   ALE22,                  &
                    ARI01,   ARI02,   ARI10,   ARI11,                  &
                    ARI12,   ARI20,   ARI21,   ARI22
  END TYPE NR_PAR_A_Y

  TYPE NR_PAR_A_X
    REAL (KIND=PP) ::                                                  &
                    ARE01,   ARE02,   ARE10,   ARE11,                  &
                    ARE12,   ARE20,   ARE21,   ARE22,                  &
                    AFR01,   AFR02,   AFR10,   AFR11,                  &
                    AFR12,   AFR20,   AFR21,   AFR22
  END TYPE NR_PAR_A_X

  TYPE NR_DATA_Z
    REAL (KIND=PP) ::                                                  &
                    TIKMZ0  , TIKMZ1  , TIKMZ2
  END TYPE NR_DATA_Z

  TYPE NR_DATA_Y
    REAL (KIND=PP) ::                                                  &
                    T1KL    , T2KL    , T3KL    ,                      &
                    TMX0KL  , TMX1KL  , TMX2KL
  END TYPE NR_DATA_Y

  TYPE NR_DATA_X
    REAL (KIND=PP) ::                                                  &
                    TI1L    , TI2L    , TI3L    ,                      &
                    TIMY0L  , TIMY1L  , TIMY2L
  END TYPE NR_DATA_X

!-----------------------------------------------------------------------

  TYPE (NR_PAR_A_X), DIMENSION (:,:), ALLOCATABLE ::                   &
                    AXU, AXV, AXW, AXUF, AXVF, AXWF

  TYPE (NR_PAR_A_Y), DIMENSION (:,:), ALLOCATABLE ::                   &
                    AYU, AYV, AYW, AYUF, AYVF, AYWF

  TYPE (NR_DATA_X), DIMENSION (:,:), ALLOCATABLE ::                    &
                             SX_U, SX_UF,                              &
                             SX_V, SX_VF,                              &
                             SX_W, SX_WF

  TYPE (NR_DATA_X), DIMENSION (:)  , ALLOCATABLE ::                    &
                            SX_UP, SX_UPF,                             &
                            SX_VP, SX_VPF,                             &
                            SX_WP, SX_WPF

  TYPE (NR_DATA_Y), DIMENSION (:,:), ALLOCATABLE ::                    &
                             SY_U, SY_UF,                              &
                             SY_V, SY_VF,                              &
                             SY_W, SY_WF

  TYPE (NR_DATA_Y), DIMENSION (:),   ALLOCATABLE ::                    &
                            SY_UP, SY_UPF,                             &
                            SY_VP, SY_VPF,                             &
                            SY_WP, SY_WPF

  TYPE (NR_DATA_Z), DIMENSION (:,:), ALLOCATABLE ::                    &
                                   SZ_U, SZ_UP,                        &
                                   SZ_V, SZ_VP,                        &
                                   SZ_W, SZ_WP

  REAL (KIND=PP)   ::                                                  &
                    ABW01,   ABW02,   ABW10,   ABW11,                  &
                    ABW12,   ABW20,   ABW21,   ABW22,                  &
                    ABH01,   ABH02,   ABH10,   ABH11,                  &
                    ABH12,   ABH20,   ABH21,   ABH22

  REAL (KIND=PP), DIMENSION(:,:,:), ALLOCATABLE ::                     &
                    XXYI2 , XXYIMX  , XXYK2   , XXYKMY  ,              &
                    XXZK2 , XXZKMY1 , XXZI2   , XXZIMX  ,              &
                    XYZK2 , XYZKMY  , XYZI2   , XYZIMX1 ,              &
                    XXXK2 , XXXI2   , XXXKMY1 , XXXIMX1 ,              &
                    XYYK2 , XYYI2   , XYYKMY1 , XYYIMX1 ,              &
                    XZZK2 , XZZI2   , XZZKMY1 , XZZIMX1 ,              &
                    XXYI2F, XXYIMXF , XXYK2F  , XXYKMYF ,              &
                    XXZK2F, XXZKMY1F, XXZI2F  , XXZIMXF ,              &
                    XYZK2F, XYZKMYF , XYZI2F  , XYZIMX1F,              &
                    XXXK2F, XXXI2F  , XXXKMY1F, XXXIMX1F,              &
                    XYYK2F, XYYI2F  , XYYKMY1F, XYYIMX1F,              &
                    XZZK2F, XZZI2F  , XZZKMY1F, XZZIMX1F

CONTAINS

  FUNCTION     NR_DATA_X_PLUS (A, B)

    TYPE(NR_DATA_X), INTENT(IN) :: A, B
    TYPE(NR_DATA_X)             :: NR_DATA_X_PLUS

    NR_DATA_X_PLUS%TI1L   = A%TI1L   + B%TI1L
    NR_DATA_X_PLUS%TI2L   = A%TI2L   + B%TI2L
    NR_DATA_X_PLUS%TI3L   = A%TI3L   + B%TI3L
    NR_DATA_X_PLUS%TIMY0L = A%TIMY0L + B%TIMY0L
    NR_DATA_X_PLUS%TIMY1L = A%TIMY1L + B%TIMY1L
    NR_DATA_X_PLUS%TIMY2L = A%TIMY2L + B%TIMY2L

  END FUNCTION NR_DATA_X_PLUS

  FUNCTION     NR_DATA_Y_PLUS (A, B)

    TYPE(NR_DATA_Y), INTENT(IN) :: A, B
    TYPE(NR_DATA_Y)             :: NR_DATA_Y_PLUS

    NR_DATA_Y_PLUS%T1KL   = A%T1KL   + B%T1KL
    NR_DATA_Y_PLUS%T2KL   = A%T2KL   + B%T2KL
    NR_DATA_Y_PLUS%T3KL   = A%T3KL   + B%T3KL
    NR_DATA_Y_PLUS%TMX0KL = A%TMX0KL + B%TMX0KL
    NR_DATA_Y_PLUS%TMX1KL = A%TMX1KL + B%TMX1KL
    NR_DATA_Y_PLUS%TMX2KL = A%TMX2KL + B%TMX2KL

  END FUNCTION NR_DATA_Y_PLUS

  FUNCTION     NR_DATA_Z_PLUS (A, B)

    TYPE(NR_DATA_Z), INTENT(IN) :: A, B
    TYPE(NR_DATA_Z)             :: NR_DATA_Z_PLUS

    NR_DATA_Z_PLUS%TIKMZ0 = A%TIKMZ0 + B%TIKMZ0
    NR_DATA_Z_PLUS%TIKMZ1 = A%TIKMZ1 + B%TIKMZ1
    NR_DATA_Z_PLUS%TIKMZ2 = A%TIKMZ2 + B%TIKMZ2

  END FUNCTION NR_DATA_Z_PLUS


  FUNCTION     NR_DATA_X_TIMES(A, B)

    REAL (PP)      , INTENT(IN) :: A
    TYPE(NR_DATA_X), INTENT(IN) :: B
    TYPE(NR_DATA_X)             :: NR_DATA_X_TIMES

    NR_DATA_X_TIMES%TI1L   = A * B%TI1L
    NR_DATA_X_TIMES%TI2L   = A * B%TI2L
    NR_DATA_X_TIMES%TI3L   = A * B%TI3L
    NR_DATA_X_TIMES%TIMY0L = A * B%TIMY0L
    NR_DATA_X_TIMES%TIMY1L = A * B%TIMY1L
    NR_DATA_X_TIMES%TIMY2L = A * B%TIMY2L

  END FUNCTION NR_DATA_X_TIMES

  FUNCTION     NR_DATA_Y_TIMES(A, B)

    REAL (PP)      , INTENT(IN) :: A
    TYPE(NR_DATA_Y), INTENT(IN) :: B
    TYPE(NR_DATA_Y)             :: NR_DATA_Y_TIMES

    NR_DATA_Y_TIMES%T1KL   = A * B%T1KL
    NR_DATA_Y_TIMES%T2KL   = A * B%T2KL
    NR_DATA_Y_TIMES%T3KL   = A * B%T3KL
    NR_DATA_Y_TIMES%TMX0KL = A * B%TMX0KL
    NR_DATA_Y_TIMES%TMX1KL = A * B%TMX1KL
    NR_DATA_Y_TIMES%TMX2KL = A * B%TMX2KL

  END FUNCTION NR_DATA_Y_TIMES

  FUNCTION     NR_DATA_Z_TIMES(A, B)

    REAL (PP)      , INTENT(IN) :: A
    TYPE(NR_DATA_Z), INTENT(IN) :: B
    TYPE(NR_DATA_Z)             :: NR_DATA_Z_TIMES

    NR_DATA_Z_TIMES%TIKMZ0 = A * B%TIKMZ0
    NR_DATA_Z_TIMES%TIKMZ1 = A * B%TIKMZ1
    NR_DATA_Z_TIMES%TIKMZ2 = A * B%TIKMZ2

  END FUNCTION NR_DATA_Z_TIMES

END MODULE NONREF_BOUND

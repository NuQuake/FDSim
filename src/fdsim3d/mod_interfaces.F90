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

MODULE  INTERFACES

  INTERFACE

    FUNCTION             IND (I,K)
      INTEGER, INTENT(IN ) :: I,K
      INTEGER              :: IND
    END FUNCTION         IND


    SUBROUTINE      ALLOC_NR
    END SUBROUTINE  ALLOC_NR

    SUBROUTINE      ALLOC_T_SEIS
    END SUBROUTINE  ALLOC_T_SEIS

    SUBROUTINE      AUX1
    END SUBROUTINE  AUX1


    SUBROUTINE      BCUCOF(Y,Y1,Y2,Y12,D1,D2,C)
      USE PRECISION, ONLY: PP
      IMPLICIT NONE
      REAL(PP),                 INTENT(IN) :: D1,D2
      REAL(PP), DIMENSION(4),   INTENT(IN) :: Y,Y1,Y2,Y12
      REAL(PP), DIMENSION(4,4), INTENT(OUT) :: C
    END SUBROUTINE  BCUCOF

    SUBROUTINE      COARSE2FINE_TXZ
    END SUBROUTINE  COARSE2FINE_TXZ

    SUBROUTINE      COARSE2FINE_TYZ
    END SUBROUTINE  COARSE2FINE_TYZ

    SUBROUTINE      COARSE2FINE_W
    END SUBROUTINE  COARSE2FINE_W

    SUBROUTINE      COMMUNICATE_STRESSFIELD
    END SUBROUTINE  COMMUNICATE_STRESSFIELD

    SUBROUTINE      COMMUNICATE_WAVEFIELD
    END SUBROUTINE  COMMUNICATE_WAVEFIELD

    SUBROUTINE      DISTRIBUTE
    END SUBROUTINE  DISTRIBUTE

    SUBROUTINE      DOMAIN_DECOMPOSITION
    END SUBROUTINE  DOMAIN_DECOMPOSITION

    SUBROUTINE      FINE2COARSE (L, L_F, VAR, VARF, I_SHIFT, K_SHIFT)
      USE PRECISION   , ONLY: PP
      USE AUXIL       , ONLY:  X_MIN_D , X_MAX_D ,                     &
                               Y_MIN_D , Y_MAX_D ,                     &
                               Z_MIN_D , Z_MAX_D ,                     &
                               X_MIN_DF, X_MAX_DF,                     &
                               Y_MIN_DF, Y_MAX_DF,                     &
                               Z_MIN_DF, Z_MAX_DF
      IMPLICIT NONE
      LOGICAL , INTENT (IN)    :: I_SHIFT, K_SHIFT
      INTEGER , INTENT (IN)    :: L, L_F
      REAL(PP), INTENT (INOUT), DIMENSION                              &
         (X_MIN_D :X_MAX_D ,Y_MIN_D :Y_MAX_D ,Z_MIN_D :Z_MAX_D ) :: VAR
      REAL(PP), INTENT (INOUT), DIMENSION                              &
         (X_MIN_DF:X_MAX_DF,Y_MIN_DF:Y_MAX_DF,Z_MIN_DF:Z_MAX_DF) :: VARF
    END SUBROUTINE  FINE2COARSE


    SUBROUTINE      HPLANE_0_STRESS
    END SUBROUTINE  HPLANE_0_STRESS

    SUBROUTINE      HPLANE_0_VELOCITY
    END SUBROUTINE  HPLANE_0_VELOCITY

    SUBROUTINE      HPLANE_0F_STRESS
    END SUBROUTINE  HPLANE_0F_STRESS

    SUBROUTINE      HPLANE_0F_VELOCITY
    END SUBROUTINE  HPLANE_0F_VELOCITY

    SUBROUTINE      HPLANE_1_STRESS
    END SUBROUTINE  HPLANE_1_STRESS

    SUBROUTINE      HPLANE_1_VELOCITY
    END SUBROUTINE  HPLANE_1_VELOCITY

    SUBROUTINE      HPLANE_1F_STRESS
    END SUBROUTINE  HPLANE_1F_STRESS

    SUBROUTINE      HPLANE_1F_VELOCITY
    END SUBROUTINE  HPLANE_1F_VELOCITY

    SUBROUTINE      HPLANE_IC_STRESS ( L )
      INTEGER, INTENT(IN) ::           L
    END SUBROUTINE  HPLANE_IC_STRESS

    SUBROUTINE      HPLANE_IC_VELOCITY ( L )
      INTEGER, INTENT(IN) ::             L
    END SUBROUTINE  HPLANE_IC_VELOCITY

    SUBROUTINE      HPLANE_IF_STRESS ( L )
      INTEGER, INTENT(IN) ::           L
    END SUBROUTINE  HPLANE_IF_STRESS

    SUBROUTINE      HPLANE_IF_VELOCITY ( L )
      INTEGER, INTENT(IN) ::             L
    END SUBROUTINE  HPLANE_IF_VELOCITY

    SUBROUTINE      HPLANE_LIN_STRESS (SOURTF_LIN, SOURTF_LOU)
      USE PRECISION, ONLY: PP
      REAL(PP), INTENT(IN) ::          SOURTF_LIN, SOURTF_LOU
    END SUBROUTINE  HPLANE_LIN_STRESS

    SUBROUTINE      HPLANE_LIN_W_STRESS (SOURTF_LIN, SOURTF_LOU)
      USE PRECISION, ONLY: PP
      REAL(PP), INTENT(IN) ::          SOURTF_LIN, SOURTF_LOU
    END SUBROUTINE  HPLANE_LIN_W_STRESS

    SUBROUTINE      HPLANE_LIN_VELOCITY
    END SUBROUTINE  HPLANE_LIN_VELOCITY

    SUBROUTINE      HPLANE_LIN_W_VELOCITY
    END SUBROUTINE  HPLANE_LIN_W_VELOCITY

    SUBROUTINE      HPLANE_LINPM1_STRESS ( L )
      INTEGER, INTENT(IN) ::               L
    END SUBROUTINE  HPLANE_LINPM1_STRESS

    SUBROUTINE      HPLANE_LINPM1_VELOCITY ( L )
      INTEGER, INTENT(IN) ::                 L
    END SUBROUTINE  HPLANE_LINPM1_VELOCITY

    SUBROUTINE      HPLANE_LPASF_STRESS
    END SUBROUTINE  HPLANE_LPASF_STRESS

    SUBROUTINE      HPLANE_LPASM1F_STRESS
    END SUBROUTINE  HPLANE_LPASM1F_STRESS

    SUBROUTINE      HPLANE_LPASM1F_VELOCITY
    END SUBROUTINE  HPLANE_LPASM1F_VELOCITY

    SUBROUTINE      HPLANE_LPASM1M2_STRESS
    END SUBROUTINE  HPLANE_LPASM1M2_STRESS

    SUBROUTINE      HPLANE_LPASM1M2_VELOCITY
    END SUBROUTINE  HPLANE_LPASM1M2_VELOCITY

    SUBROUTINE      HPLANE_MZ1_STRESS
    END SUBROUTINE  HPLANE_MZ1_STRESS

    SUBROUTINE      HPLANE_PML_STRESS ( L, LL )
      INTEGER, INTENT(IN) ::            L, LL
    END SUBROUTINE  HPLANE_PML_STRESS

    SUBROUTINE      HPLANE_PML_VELOCITY ( L, LL )
      INTEGER, INTENT(IN) ::              L, LL
    END SUBROUTINE  HPLANE_PML_VELOCITY

    SUBROUTINE      HPLANE_PMLZH0_STRESS
    END SUBROUTINE  HPLANE_PMLZH0_STRESS

    SUBROUTINE      HPLANE_PMLZH0_VELOCITY
    END SUBROUTINE  HPLANE_PMLZH0_VELOCITY

    SUBROUTINE      INH_NR                                             &
                    ( MXL, MXH, PXH, MYL, MYH, PYH,                    &
                      X_MIN, X_MAX, Y_MIN, Y_MAX,                      &
                      AXT, AYT, SX, SXP, SY, SYP, T, K1, I1, TPML)
      USE PRECISION   , ONLY: PP
      USE NONREF_BOUND, ONLY: NR_PAR_A_X, NR_DATA_X,                   &
                              NR_PAR_A_Y, NR_DATA_Y
      IMPLICIT NONE
      INTEGER          ,                           INTENT (IN)    ::   &
                      X_MIN, X_MAX, Y_MIN, Y_MAX,                      &
                             K1, I1, MXL, MXH, PXH, MYL, MYH, PYH, TPML
      TYPE (NR_PAR_A_X), DIMENSION (MXL:MXH), INTENT (INOUT) :: AXT
      TYPE (NR_PAR_A_Y), DIMENSION (MYL:MYH), INTENT (INOUT) :: AYT
      TYPE (NR_DATA_X) , DIMENSION (MXL:MXH), INTENT (INOUT) :: SX, SXP
      TYPE (NR_DATA_Y) , DIMENSION (MYL:MYH), INTENT (INOUT) :: SY, SYP
      REAL(PP), DIMENSION (MXL:MXH, MYL:MYH), INTENT (INOUT) :: T
    END SUBROUTINE  INH_NR

    SUBROUTINE      INH_V_U_2ND                                        &
                       (L, MXL, MXH, PL, MXT, PXH, MX1,                &
                           MYL, MYH,     MYT, PYH, MY1,                &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX, Y_MIN, Y_MAX,                 &
                           JMT, DENUT, UT,                             &
                           TXXT, TXYT, TXZT,                           &
                                           RXXLT, RXXHT,               &
                                           RXYLT, RXYHT,               &
                 ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,                   &
                 ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH,                   &
                  ODA_YL,  ODB_YL,  ODC_YL,  ODD_YL,                   &
                  ODA_YH,  ODB_YH,  ODC_YH,  ODD_YH )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                        L, MXL, MXH, PL, MXT, MX1, PXH,&
                                           MYL, MYH,     MYT, MY1, PYH,&
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX, Y_MIN, Y_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,MYL:MYH,MZL:MZH), INTENT(INOUT) :: &
                                                                   JMT
      REAL   (PP),DIMENSION(1:JMNUM                ), INTENT(INOUT) :: &
                                                                  DENUT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH,MZL:MZH), INTENT(INOUT) :: &
                                                       TXXT, TXYT, TXZT
      REAL   (PP),DIMENSION(PL :3  ,MYL:MYH,MZL:MZH), INTENT(INOUT) :: &
                                                                  RXXLT
      REAL   (PP),DIMENSION(MX1:PXH,MYL:MYH,MZL:MZH), INTENT(INOUT) :: &
                                                                  RXXHT
      REAL   (PP),DIMENSION(MXL:MXH,PL :2  ,MZL:MZH), INTENT(INOUT) :: &
                                                                  RXYLT
      REAL   (PP),DIMENSION(MXL:MXH,MY1:PYH,MZL:MZH), INTENT(INOUT) :: &
                                                                  RXYHT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH,MZL:MZH), INTENT(INOUT) :: &
                                                                     UT
      REAL   (PP),DIMENSION( 0:TPML,MXL:MXH,MZL:MZH), INTENT(INOUT) :: &
                                    ODA_YL,  ODB_YL,  ODC_YL,  ODD_YL, &
                                    ODA_YH,  ODB_YH,  ODC_YH,  ODD_YH
      REAL   (PP),DIMENSION(0:TPML+1,MYL:MYH,MZL:MZH),INTENT(INOUT) :: &
                                   ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL, &
                                   ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH
    END SUBROUTINE  INH_V_U_2ND

    SUBROUTINE      INH_V_U_2ND_0                                      &
                       (   MXL, MXH, PL, MXT, PXH, MX1,                &
                           MYL, MYH,     MYT, PYH, MY1,                &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX, Y_MIN, Y_MAX,                 &
                           JMT, DENUT, UT,                             &
                           TXXT, TXYT, TXZT,                           &
                                           RXXLT, RXXHT,               &
                                           RXYLT, RXYHT,               &
                 ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,                   &
                 ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH,                   &
                  ODA_YL,  ODB_YL,  ODC_YL,  ODD_YL,                   &
                  ODA_YH,  ODB_YH,  ODC_YH,  ODD_YH )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                           MXL, MXH, PL, MXT, MX1, PXH,&
                                           MYL, MYH,     MYT, MY1, PYH,&
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX, Y_MIN, Y_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,MYL:MYH,MZL:MZH), INTENT(INOUT) :: &
                                                                   JMT
      REAL   (PP),DIMENSION(1:JMNUM                ), INTENT(INOUT) :: &
                                                                  DENUT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH,MZL:MZH), INTENT(INOUT) :: &
                                                       TXXT, TXYT, TXZT
      REAL   (PP),DIMENSION(PL :3  ,MYL:MYH,MZL:MZH), INTENT(INOUT) :: &
                                                                  RXXLT
      REAL   (PP),DIMENSION(MX1:PXH,MYL:MYH,MZL:MZH), INTENT(INOUT) :: &
                                                                  RXXHT
      REAL   (PP),DIMENSION(MXL:MXH,PL :2  ,MZL:MZH), INTENT(INOUT) :: &
                                                                  RXYLT
      REAL   (PP),DIMENSION(MXL:MXH,MY1:PYH,MZL:MZH), INTENT(INOUT) :: &
                                                                  RXYHT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH,MZL:MZH), INTENT(INOUT) :: &
                                                                     UT
      REAL   (PP),DIMENSION( 0:TPML,MXL:MXH,MZL:MZH), INTENT(INOUT) :: &
                                    ODA_YL,  ODB_YL,  ODC_YL,  ODD_YL, &
                                    ODA_YH,  ODB_YH,  ODC_YH,  ODD_YH
      REAL   (PP),DIMENSION(0:TPML+1,MYL:MYH,MZL:MZH),INTENT(INOUT) :: &
                                   ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL, &
                                   ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH
    END SUBROUTINE  INH_V_U_2ND_0

    SUBROUTINE      INH_V_U_INT_2ND                                    &
                   (L,     MXL, MXH, PL, MXT, PXH, MX1,                &
                           MYL, MYH,     MYT, PYH, MY1,                &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX, Y_MIN, Y_MAX,                 &
                           JMT, DENUT, UT,                             &
                           TXXT, TXYT, TXZT,                           &
                                           RXXLT, RXXHT,               &
                                           RXYLT, RXYHT,               &
                 ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,                   &
                 ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH,                   &
                  ODA_YL,  ODB_YL,  ODC_YL,  ODD_YL,                   &
                  ODA_YH,  ODB_YH,  ODC_YH,  ODD_YH )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                    L,     MXL, MXH, PL, MXT, MX1, PXH,&
                                           MYL, MYH,     MYT, MY1, PYH,&
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX, Y_MIN, Y_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,MYL:MYH,MZL:MZH), INTENT(INOUT) :: &
                                                                   JMT
      REAL   (PP),DIMENSION(1:JMNUM                ), INTENT(INOUT) :: &
                                                                  DENUT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH,MZL:MZH), INTENT(INOUT) :: &
                                                       TXXT, TXYT, TXZT
      REAL   (PP),DIMENSION(PL :3  ,MYL:MYH,MZL:MZH), INTENT(INOUT) :: &
                                                                  RXXLT
      REAL   (PP),DIMENSION(MX1:PXH,MYL:MYH,MZL:MZH), INTENT(INOUT) :: &
                                                                  RXXHT
      REAL   (PP),DIMENSION(MXL:MXH,PL :2  ,MZL:MZH), INTENT(INOUT) :: &
                                                                  RXYLT
      REAL   (PP),DIMENSION(MXL:MXH,MY1:PYH,MZL:MZH), INTENT(INOUT) :: &
                                                                  RXYHT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH,MZL:MZH), INTENT(INOUT) :: &
                                                                     UT
      REAL   (PP),DIMENSION( 0:TPML,MXL:MXH,MZL:MZH), INTENT(INOUT) :: &
                                    ODA_YL,  ODB_YL,  ODC_YL,  ODD_YL, &
                                    ODA_YH,  ODB_YH,  ODC_YH,  ODD_YH
      REAL   (PP),DIMENSION(0:TPML+1,MYL:MYH,MZL:MZH),INTENT(INOUT) :: &
                                   ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL, &
                                   ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH
    END SUBROUTINE  INH_V_U_INT_2ND

    SUBROUTINE      INH_V_U_INT_2ND_PML                                &
                  ( L, LL, MXL, MXH, PL, MXT, PXH, MX1,                &
                           MYL, MYH,     MYT, PYH, MY1,                &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX, Y_MIN, Y_MAX,                 &
                           JMT, DENUT, UT,                             &
                           TXXT, TXYT, TXZT,                           &
                                           RXXLT, RXXHT,               &
                                           RXYLT, RXYHT,               &
                                                  RXZHT,               &
                 ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,                   &
                 ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH,                   &
                  ODA_YL,  ODB_YL,  ODC_YL,  ODD_YL,                   &
                  ODA_YH,  ODB_YH,  ODC_YH,  ODD_YH,                   &
                  ODA_ZH,  ODB_ZH,  ODC_ZH,  ODD_ZH )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                    L, LL, MXL, MXH, PL, MXT, MX1, PXH,&
                                           MYL, MYH,     MYT, MY1, PYH,&
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX, Y_MIN, Y_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                                                  DENUT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                       TXXT, TXYT, TXZT
      REAL   (PP),DIMENSION(PL :3  ,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RXXLT
      REAL   (PP),DIMENSION(MX1:PXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RXXHT
      REAL   (PP),DIMENSION(MXL:MXH,PL :2  , MZL:MZH),INTENT(INOUT) :: &
                                                                  RXYLT
      REAL   (PP),DIMENSION(MXL:MXH,MY1:PYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RXYHT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH,  0:TPML),INTENT(INOUT) :: &
                                                                  RXZHT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                     UT
      REAL   (PP),DIMENSION( 0:TPML,MXL:MXH,MZL:MZH), INTENT(INOUT) :: &
                                    ODA_YL,  ODB_YL,  ODC_YL,  ODD_YL, &
                                    ODA_YH,  ODB_YH,  ODC_YH,  ODD_YH
      REAL   (PP),DIMENSION(0:TPML+1,MYL:MYH,MZL:MZH),INTENT(INOUT) :: &
                                   ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL, &
                                   ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH
      REAL   (PP),DIMENSION( 0:TPML,MXL:MXH, MYL:MYH),INTENT(INOUT) :: &
                                    ODA_ZH,  ODB_ZH,  ODC_ZH,  ODD_ZH
    END SUBROUTINE  INH_V_U_INT_2ND_PML

    SUBROUTINE      INH_V_UV_INT                                       &
                        (L, MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,    &
                            X_MIN, X_MAX, Y_MIN, Y_MAX,                &
                            JMT, DENUT, DENVT, UT, VT,                 &
                            TXXT, TXYT, TXZT, TYYT,       TYZT )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                             L, MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,&
                             X_MIN, X_MAX, Y_MIN, Y_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                                           DENUT, DENVT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                           TXXT, TXYT, TYYT, TXZT, TYZT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                               UT  , VT
    END SUBROUTINE  INH_V_UV_INT

    SUBROUTINE      INH_V_UV_INT_0                                     &
                          ( MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,    &
                            X_MIN, X_MAX, Y_MIN, Y_MAX,                &
                            JMT, DENUT, DENVT, UT, VT,                 &
                            TXXT, TXYT, TXZT, TYYT,       TYZT )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                               MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH, &
                               X_MIN, X_MAX, Y_MIN, Y_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                                           DENUT, DENVT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                           TXXT, TXYT, TYYT, TXZT, TYZT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                               UT  , VT
    END SUBROUTINE  INH_V_UV_INT_0

    SUBROUTINE      INH_V_UV_INT_1_SI                                  &
                          ( MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,    &
                            X_MIN, X_MAX, Y_MIN, Y_MAX, TPML,          &
                            UT, VT, WT )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                               MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH, &
                               X_MIN, X_MAX, Y_MIN, Y_MAX, TPML
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH,MZL:MZH), INTENT(INOUT) :: &
                                                          UT, VT, WT
    END SUBROUTINE  INH_V_UV_INT_1_SI

    SUBROUTINE      INH_V_V_2ND                                        &
                       (L, MXL, MXH, PL, MXT, PXH, MX1,                &
                           MYL, MYH,     MYT, PYH, MY1,                &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX, Y_MIN, Y_MAX,                 &
                           JMT, DENVT, VT,                             &
                           TXYT, TYYT, TYZT,                           &
                                           RYYLT, RYYHT,               &
                                           RYXLT, RYXHT,               &
                  ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,                   &
                  ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH,                   &
                 ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL,                   &
                 ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                        L, MXL, MXH, PL, MXT, MX1, PXH,&
                                           MYL, MYH,     MYT, MY1, PYH,&
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX, Y_MIN, Y_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                                                  DENVT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                       TXYT, TYYT, TYZT
      REAL   (PP),DIMENSION(PL :2  ,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RYXLT
      REAL   (PP),DIMENSION(MX1:PXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RYXHT
      REAL   (PP),DIMENSION(MXL:MXH,PL :3  , MZL:MZH),INTENT(INOUT) :: &
                                                                  RYYLT
      REAL   (PP),DIMENSION(MXL:MXH,MY1:PYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RYYHT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                     VT
      REAL   (PP),DIMENSION( 0:TPML,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                    ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,                 &
                    ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH
      REAL   (PP),DIMENSION(0:TPML+1,MXL:MXH, MZL:MZH),INTENT(INOUT):: &
                   ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL,                 &
                   ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH
    END SUBROUTINE  INH_V_V_2ND

    SUBROUTINE      INH_V_V_2ND_0                                      &
                       (   MXL, MXH, PL, MXT, PXH, MX1,                &
                           MYL, MYH,     MYT, PYH, MY1,                &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX, Y_MIN, Y_MAX,                 &
                           JMT, DENVT, VT,                             &
                           TXYT, TYYT, TYZT,                           &
                                           RYYLT, RYYHT,               &
                                           RYXLT, RYXHT,               &
                  ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,                   &
                  ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH,                   &
                 ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL,                   &
                 ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                           MXL, MXH, PL, MXT, MX1, PXH,&
                                           MYL, MYH,     MYT, MY1, PYH,&
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX, Y_MIN, Y_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                                                  DENVT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                       TXYT, TYYT, TYZT
      REAL   (PP),DIMENSION(PL :2  ,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RYXLT
      REAL   (PP),DIMENSION(MX1:PXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RYXHT
      REAL   (PP),DIMENSION(MXL:MXH,PL :3  , MZL:MZH),INTENT(INOUT) :: &
                                                                  RYYLT
      REAL   (PP),DIMENSION(MXL:MXH,MY1:PYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RYYHT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                     VT
      REAL   (PP),DIMENSION( 0:TPML,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                    ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,                 &
                    ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH
      REAL   (PP),DIMENSION(0:TPML+1,MXL:MXH, MZL:MZH),INTENT(INOUT):: &
                   ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL,                 &
                   ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH
    END SUBROUTINE  INH_V_V_2ND_0

    SUBROUTINE      INH_V_V_INT_2ND                                    &
                   (L,     MXL, MXH, PL, MXT, PXH, MX1,                &
                           MYL, MYH,     MYT, PYH, MY1,                &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX, Y_MIN, Y_MAX,                 &
                           JMT, DENVT, VT,                             &
                           TXYT, TYYT, TYZT,                           &
                                           RYYLT, RYYHT,               &
                                           RYXLT, RYXHT,               &
                  ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,                   &
                  ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH,                   &
                 ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL,                   &
                 ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                    L,     MXL, MXH, PL, MXT, MX1, PXH,&
                                           MYL, MYH,     MYT, MY1, PYH,&
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX, Y_MIN, Y_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                                                  DENVT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                       TXYT, TYYT, TYZT
      REAL   (PP),DIMENSION(PL :2  ,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RYXLT
      REAL   (PP),DIMENSION(MX1:PXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RYXHT
      REAL   (PP),DIMENSION(MXL:MXH,PL :3  , MZL:MZH),INTENT(INOUT) :: &
                                                                  RYYLT
      REAL   (PP),DIMENSION(MXL:MXH,MY1:PYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RYYHT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                     VT
      REAL   (PP),DIMENSION( 0:TPML,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                    ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,                 &
                    ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH
      REAL   (PP),DIMENSION(0:TPML+1,MXL:MXH, MZL:MZH),INTENT(INOUT):: &
                   ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL,                 &
                   ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH
    END SUBROUTINE  INH_V_V_INT_2ND

    SUBROUTINE      INH_V_V_INT_2ND_PML                                &
                  ( L, LL, MXL, MXH, PL, MXT, PXH, MX1,                &
                           MYL, MYH,     MYT, PYH, MY1,                &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX, Y_MIN, Y_MAX,                 &
                           JMT, DENVT, VT,                             &
                           TXYT, TYYT, TYZT,                           &
                                           RYYLT, RYYHT,               &
                                           RYXLT, RYXHT,               &
                                                  RYZHT,               &
                  ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,                   &
                  ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH,                   &
                 ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL,                   &
                 ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH,                   &
                  ODA_ZH,  ODB_ZH,  ODC_ZH,  ODD_ZH )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                    L, LL, MXL, MXH, PL, MXT, MX1, PXH,&
                                           MYL, MYH,     MYT, MY1, PYH,&
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX, Y_MIN, Y_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                                                  DENVT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                       TXYT, TYYT, TYZT
      REAL   (PP),DIMENSION(PL :2  ,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RYXLT
      REAL   (PP),DIMENSION(MX1:PXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RYXHT
      REAL   (PP),DIMENSION(MXL:MXH,PL :3  , MZL:MZH),INTENT(INOUT) :: &
                                                                  RYYLT
      REAL   (PP),DIMENSION(MXL:MXH,MY1:PYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RYYHT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH,  0:TPML),INTENT(INOUT) :: &
                                                                  RYZHT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                     VT
      REAL   (PP),DIMENSION( 0:TPML,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                    ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,                 &
                    ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH
      REAL   (PP),DIMENSION(0:TPML+1,MXL:MXH, MZL:MZH),INTENT(INOUT):: &
                   ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL,                 &
                   ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH
      REAL   (PP),DIMENSION( 0:TPML,MXL:MXH, MYL:MYH),INTENT(INOUT) :: &
                    ODA_ZH,  ODB_ZH,  ODC_ZH,  ODD_ZH
    END SUBROUTINE  INH_V_V_INT_2ND_PML

    SUBROUTINE      INH_V_W_2ND                                        &
                       (L, MXL, MXH, PL, MXT, PXH, MX1,                &
                           MYL, MYH,     MYT, PYH, MY1,                &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX, Y_MIN, Y_MAX,                 &
                           JMT, DENWT, WT,                             &
                           TZZT, TXZT, TYZT,                           &
                                           RZXLT, RZXHT,               &
                                           RZYLT, RZYHT,               &
                  ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,                   &
                  ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH,                   &
                  ODA_YL,  ODB_YL,  ODC_YL,  ODD_YL,                   &
                  ODA_YH,  ODB_YH,  ODC_YH,  ODD_YH )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                             L, MXL, MXH, PL, MXT, PXH,&
                                                MYL, MYH,     MYT, PYH,&
                                                MZL, MZH,  TPML,       &
                                                MX1, MY1,              &
                                             X_MIN, X_MAX, Y_MIN, Y_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                                                  DENWT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                       TXZT, TYZT, TZZT
      REAL   (PP),DIMENSION(PL :2  ,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RZXLT
      REAL   (PP),DIMENSION(MX1:PXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RZXHT
      REAL   (PP),DIMENSION(MXL:MXH,PL :2  , MZL:MZH),INTENT(INOUT) :: &
                                                                  RZYLT
      REAL   (PP),DIMENSION(MXL:MXH,MY1:PYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RZYHT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                     WT
      REAL   (PP),DIMENSION( 0:TPML,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                    ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,                 &
                    ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH
      REAL   (PP),DIMENSION( 0:TPML,MXL:MXH, MZL:MZH),INTENT(INOUT) :: &
                    ODA_YL,  ODB_YL,  ODC_YL,  ODD_YL,                 &
                    ODA_YH,  ODB_YH,  ODC_YH,  ODD_YH
    END SUBROUTINE  INH_V_W_2ND

    SUBROUTINE      INH_V_W_2ND_0                                      &
                         ( MXL, MXH,     MXT,                          &
                           MYL, MYH,     MYT,                          &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX, Y_MIN, Y_MAX,                 &
                           JMT, DENWT, WT, TZZT )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                    MXL, MXH, MYL, MYH, MXT, MYT, TPML,&
                                    MZL, MZH,                          &
                                    X_MIN, X_MAX, Y_MIN, Y_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                                                  DENWT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                   TZZT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                     WT
    END SUBROUTINE  INH_V_W_2ND_0

    SUBROUTINE      INH_V_W_INT                                        &
                       ( L, MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,    &
                            X_MIN, X_MAX, Y_MIN, Y_MAX,                &
                            JMT, DENWT, WT, TZZT, TXZT, TYZT )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                            L, MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH, &
                            X_MIN, X_MAX, Y_MIN, Y_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                                                  DENWT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                       TXZT, TYZT, TZZT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                     WT
    END SUBROUTINE  INH_V_W_INT

    SUBROUTINE      INH_V_W_INT_0                                      &
                       (    MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,    &
                            X_MIN, X_MAX, Y_MIN, Y_MAX,                &
                            JMT, DENWT, WT, TZZT )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                               MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH, &
                            X_MIN, X_MAX, Y_MIN, Y_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                                                  DENWT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  TZZT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                     WT
    END SUBROUTINE  INH_V_W_INT_0

    SUBROUTINE      INH_V_W_INT_1                                      &
                       (    MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,    &
                            X_MIN, X_MAX, Y_MIN, Y_MAX,                &
                            JMT, DENWT, WT, TZZT, TXZT, TYZT )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                               MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH, &
                            X_MIN, X_MAX, Y_MIN, Y_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                                                  DENWT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                       TXZT, TYZT, TZZT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                     WT
    END SUBROUTINE  INH_V_W_INT_1

    SUBROUTINE      INH_V_W_INT_1_SI                                   &
                          ( MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,    &
                            X_MIN, X_MAX, Y_MIN, Y_MAX, TPML,          &
                            JMT, UT, VT, WT )
      USE PRECISION   , ONLY: PP, IK
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                               MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH, &
                               X_MIN, X_MAX, Y_MIN, Y_MAX, TPML
      INTEGER(IK),DIMENSION(MXL:MXH,MYL:MYH,MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH,MZL:MZH), INTENT(INOUT) :: &
                                                             UT, VT, WT
    END SUBROUTINE  INH_V_W_INT_1_SI

    SUBROUTINE      INH_V_W_INT_2ND                                    &
                       (L, MXL, MXH, PL, MXT, PXH, MX1,                &
                           MYL, MYH,     MYT, PYH, MY1,                &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX, Y_MIN, Y_MAX,                 &
                           JMT, DENWT, WT,                             &
                           TZZT, TXZT, TYZT,                           &
                                           RZXLT, RZXHT,               &
                                           RZYLT, RZYHT,               &
                  ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,                   &
                  ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH,                   &
                  ODA_YL,  ODB_YL,  ODC_YL,  ODD_YL,                   &
                  ODA_YH,  ODB_YH,  ODC_YH,  ODD_YH )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                             L, MXL, MXH, PL, MXT, PXH,&
                                                MYL, MYH,     MYT, PYH,&
                                                MZL, MZH,  TPML,       &
                                                MX1, MY1,              &
                                             X_MIN, X_MAX, Y_MIN, Y_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                                                  DENWT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                       TXZT, TYZT, TZZT
      REAL   (PP),DIMENSION(PL :2  ,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RZXLT
      REAL   (PP),DIMENSION(MX1:PXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RZXHT
      REAL   (PP),DIMENSION(MXL:MXH,PL :2  , MZL:MZH),INTENT(INOUT) :: &
                                                                  RZYLT
      REAL   (PP),DIMENSION(MXL:MXH,MY1:PYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RZYHT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                     WT
      REAL   (PP),DIMENSION( 0:TPML,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                    ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,                 &
                    ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH
      REAL   (PP),DIMENSION( 0:TPML,MXL:MXH, MZL:MZH),INTENT(INOUT) :: &
                    ODA_YL,  ODB_YL,  ODC_YL,  ODD_YL,                 &
                    ODA_YH,  ODB_YH,  ODC_YH,  ODD_YH
    END SUBROUTINE  INH_V_W_INT_2ND

    SUBROUTINE      INH_V_W_INT_2ND_PML                                &
                  ( L, LL, MXL, MXH, PL, MXT, PXH, MX1,                &
                           MYL, MYH,     MYT, PYH, MY1,                &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX, Y_MIN, Y_MAX,                 &
                           JMT, DENWT, WT,                             &
                           TZZT, TXZT, TYZT,                           &
                                           RZXLT, RZXHT,               &
                                           RZYLT, RZYHT,               &
                                                  RZZHT,               &
                  ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,                   &
                  ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH,                   &
                  ODA_YL,  ODB_YL,  ODC_YL,  ODD_YL,                   &
                  ODA_YH,  ODB_YH,  ODC_YH,  ODD_YH,                   &
                 ODA2_ZH, ODB2_ZH, ODC2_ZH, ODD2_ZH )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                         L, LL, MXL, MXH, PL, MXT, PXH,&
                                                MYL, MYH,     MYT, PYH,&
                                                MZL, MZH,  TPML,       &
                                                MX1, MY1,              &
                                            X_MIN, X_MAX, Y_MIN, Y_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                                                  DENWT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                       TXZT, TYZT, TZZT
      REAL   (PP),DIMENSION(PL :2  ,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RZXLT
      REAL   (PP),DIMENSION(MX1:PXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RZXHT
      REAL   (PP),DIMENSION(MXL:MXH,PL :2  , MZL:MZH),INTENT(INOUT) :: &
                                                                  RZYLT
      REAL   (PP),DIMENSION(MXL:MXH,MY1:PYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  RZYHT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                     WT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH,  0:TPML),INTENT(INOUT) :: &
                                                                  RZZHT
      REAL   (PP),DIMENSION( 0:TPML,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                    ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,                 &
                    ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH
      REAL   (PP),DIMENSION( 0:TPML,MXL:MXH, MZL:MZH),INTENT(INOUT) :: &
                    ODA_YL,  ODB_YL,  ODC_YL,  ODD_YL,                 &
                    ODA_YH,  ODB_YH,  ODC_YH,  ODD_YH
      REAL   (PP),DIMENSION(0:TPML+1,MXL:MXH, MYL:MYH),INTENT(INOUT):: &
                   ODA2_ZH, ODB2_ZH, ODC2_ZH, ODD2_ZH
    END SUBROUTINE  INH_V_W_INT_2ND_PML

    SUBROUTINE      LANCZOS
    END SUBROUTINE  LANCZOS

    SUBROUTINE      OPEN_SR_FILE
    END SUBROUTINE  OPEN_SR_FILE

    SUBROUTINE      PPML
    END SUBROUTINE  PPML

    SUBROUTINE      PPMLF
    END SUBROUTINE  PPMLF

    SUBROUTINE      PRE_NR_COEF1
    END SUBROUTINE  PRE_NR_COEF1

    SUBROUTINE      PRE_NR_COEF2
    END SUBROUTINE  PRE_NR_COEF2

    SUBROUTINE      READ_MODEL
    END SUBROUTINE  READ_MODEL

    SUBROUTINE      READ_MODEL_PAR
    END SUBROUTINE  READ_MODEL_PAR

    SUBROUTINE      READ_OPEN         ( WORK_NAME )
      IMPLICIT NONE
      CHARACTER (LEN=*), INTENT (IN) :: WORK_NAME
    END SUBROUTINE  READ_OPEN

    SUBROUTINE      SRC_ADD ( TILE )
      INTEGER , INTENT(IN) :: TILE
    END SUBROUTINE  SRC_ADD

    SUBROUTINE      STO_BOR                                            &
                    ( MXL, MXH, MYL, MYH, PXH, PYH, I1, K1,            &
                      X_MIN, X_MAX, Y_MIN, Y_MAX, T,                   &
                      SX_T, SY_T, SX_TP, SY_TP, TPML )
      USE PRECISION   , ONLY: PP
      USE NONREF_BOUND, ONLY: NR_DATA_X, NR_DATA_Y
      IMPLICIT NONE
      INTEGER                                  , INTENT (IN)    ::     &
                                     MXL, MXH, MYL, MYH, I1, K1, TPML, &
                                     PXH, PYH,                         &
                                     X_MIN, X_MAX, Y_MIN, Y_MAX
      TYPE(NR_DATA_X),DIMENSION(MXL:MXH)        ,INTENT (INOUT) ::     &
                                                            SX_T, SX_TP
      TYPE(NR_DATA_Y),DIMENSION(MYL:MYH)        ,INTENT (INOUT) ::     &
                                                            SY_T, SY_TP
      REAL(KIND=PP)  ,DIMENSION(MXL:MXH,MYL:MYH),INTENT (INOUT) :: T
    END SUBROUTINE  STO_BOR

    SUBROUTINE      STORE_SEISMOGRAMS (T)
      INTEGER, INTENT (IN) ::          T
    END SUBROUTINE  STORE_SEISMOGRAMS

    SUBROUTINE      STORE_SNAPSHOTS (ITIST)
      INTEGER, INTENT (IN) ::        ITIST
    END SUBROUTINE  STORE_SNAPSHOTS

    SUBROUTINE      STORE_SNAPSHOTS_F (ITIST)
      INTEGER, INTENT (IN) ::          ITIST
    END SUBROUTINE  STORE_SNAPSHOTS_F

    SUBROUTINE      STRAIN_XX                                          &
                      (L, MXL, MXH, PL, MXT, PXH,                      &
                          MYL, MYH,     MYT,                           &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX, Y_MIN, Y_MAX,                  &
                          A, B, C, TPML,                               &
                          UTM, XXXT, EXXT, PXXLT, PXXHT,               &
                          ODA_XL, ODB_XL, ODC_XL, ODD_XL,              &
                          ODA_XH, ODB_XH, ODC_XH, ODD_XH )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                                               MXL, MXH, PL, MXT, PXH, &
                                               MYL, MYH,     MYT,      &
                                               MZL, MZH, L, TPML,      &
                                           X_MIN, X_MAX, Y_MIN, Y_MAX
      REAL   (PP),                                    INTENT (IN)   :: &
                                                               A, B, C
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH         ),INTENT(INOUT) :: &
                                                                  EXXT
      REAL   (PP),DIMENSION(PL :  1,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                 PXXLT
      REAL   (PP),DIMENSION(MXT:PXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                 PXXHT
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                   UTM
      REAL   (PP),DIMENSION(MXL:MXH,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
                                                                  XXXT
      REAL   (PP),DIMENSION( 0:TPML,MYL:MYH, MZL:MZH),INTENT(INOUT) :: &
        ODA_XL, ODB_XL, ODC_XL, ODD_XL, ODA_XH, ODB_XH, ODC_XH, ODD_XH
    END SUBROUTINE  STRAIN_XX

    SUBROUTINE      STRAIN_XY                                          &
                      (L, MXL, MXH,     MXT, PXH,                      &
                          MYL, MYH, PL, MYT, PYH,                      &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX, Y_MIN, Y_MAX,                  &
                          A, B, C, TPML,                               &
                       UTM, VTM, XXYT, EXYT,                           &
                          PYXLT, PYXHT,                                &
                          PXYLT, PXYHT,                                &
                          ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,          &
                          ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH,          &
                          ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL,          &
                          ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                       INTENT (IN)    :: &
                                               MXL, MXH,     MXT, PXH, &
                                               MYL, MYH, PL, MYT, PYH, &
                                               MZL, MZH, L, TPML,      &
                                            X_MIN, X_MAX, Y_MIN, Y_MAX
      REAL(PP),                                      INTENT (IN)    :: &
                                                               A, B, C
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EXYT
      REAL(PP),DIMENSION (PL :  2, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PYXLT
      REAL(PP),DIMENSION (MXT:PXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PYXHT
      REAL(PP),DIMENSION (MXL:MXH, MYL:  2, MZL:MZH), INTENT(INOUT) :: &
                                                                 PXYLT
      REAL(PP),DIMENSION (MXL:MXH, MYT:PYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PXYHT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                              UTM, VTM
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XXYT
      REAL(PP),DIMENSION(0:TPML+1, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                  ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,                  &
                  ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH
      REAL(PP),DIMENSION(0:TPML+1, MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
                  ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL,                  &
                  ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH
    END SUBROUTINE  STRAIN_XY

    SUBROUTINE      STRAIN_XZ                                          &
                     ( L, MXL, MXH, PL, MXT, PXH,                      &
                          MYL, MYH,     MYT,                           &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX, Y_MIN, Y_MAX,                  &
                          A, B, C, TPML,                               &
                          UTM, WTM, XXZT, EXZT,                        &
                          PZXLT, PZXHT,                                &
                          ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,          &
                          ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                                               MXL, MXH, PL, MXT, PXH, &
                                               MYL, MYH,     MYT,      &
                                               MZL, MZH, L, TPML,      &
                                            X_MIN, X_MAX, Y_MIN, Y_MAX
      REAL(PP),                                       INTENT(IN)    :: &
                                                               A, B, C
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EXZT
      REAL(PP),DIMENSION (PL :  2, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXLT
      REAL(PP),DIMENSION (MXT:PXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXHT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                              UTM, WTM
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XXZT
      REAL(PP),DIMENSION(0:TPML+1, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                  ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,                  &
                  ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH
    END SUBROUTINE  STRAIN_XZ

    SUBROUTINE      STRAIN_XZ_1                                        &
                     (    MXL, MXH, PL, MXT, PXH,                      &
                          MYL, MYH,     MYT,                           &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX, Y_MIN, Y_MAX,                  &
                          A, B, C, TPML,                               &
                          UTM, WTM, XXZT, EXZT,                        &
                          PZXLT, PZXHT,                                &
                          ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,          &
                          ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                                               MXL, MXH, PL, MXT, PXH, &
                                               MYL, MYH,     MYT,      &
                                               MZL, MZH,    TPML,      &
                                            X_MIN, X_MAX, Y_MIN, Y_MAX
      REAL(PP),                                       INTENT(IN)    :: &
                                                               A, B, C
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EXZT
      REAL(PP),DIMENSION (PL :  2, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXLT
      REAL(PP),DIMENSION (MXT:PXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXHT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                              UTM, WTM
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XXZT
      REAL(PP),DIMENSION(0:TPML+1, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                  ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,                  &
                  ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH
    END SUBROUTINE  STRAIN_XZ_1

    SUBROUTINE      STRAIN_XZ_2ND                                      &
                         ( L, C, TPML, X_MIN, X_MAX, Y_MIN, Y_MAX,     &
                          MXL, MXH, PL, MXT, PXH,                      &
                          MYL, MYH,     MYT,                           &
                          MZL, MZH,                                    &
                           UTM, WTM, XXZT, EXZT,                       &
                          PZXLT, PZXHT,                                &
                          ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,          &
                          ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                                               MXL, MXH, PL, MXT, PXH, &
                                               MYL, MYH,     MYT,      &
                                               MZL, MZH, L, TPML,      &
                                            X_MIN, X_MAX, Y_MIN, Y_MAX
      REAL(PP),                                       INTENT(IN)    :: &
                                                               C
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EXZT
      REAL(PP),DIMENSION (PL :  2, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXLT
      REAL(PP),DIMENSION (MXT:PXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXHT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                              UTM, WTM
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XXZT
      REAL(PP),DIMENSION(0:TPML+1, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                  ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,                  &
                  ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH
    END SUBROUTINE  STRAIN_XZ_2ND

    SUBROUTINE      STRAIN_XZ_2ND_PML                                  &
                             ( L, LL, C, TPML,                         &
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
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                                               MXL, MXH, PL, MXT, PXH, &
                                               MYL, MYH,     MYT,      &
                                               MZL, MZH, L, LL, TPML,  &
                                           X_MIN, X_MAX, Y_MIN, Y_MAX
      REAL(PP),                                       INTENT(IN)    :: &
                                                               C
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EXZT
      REAL(PP),DIMENSION (PL :  2, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXLT
      REAL(PP),DIMENSION (MXT:PXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXHT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH,0:TPML+1), INTENT(INOUT) :: &
                                                                 PXZHT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                              UTM, WTM
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XXZT
      REAL(PP),DIMENSION(0:TPML+1, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                  ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,                  &
                  ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH
      REAL(PP),DIMENSION(0:TPML+1, MXL:MXH, MYL:MYH), INTENT(INOUT) :: &
                  ODA2_ZH, ODB2_ZH, ODC2_ZH, ODD2_ZH
    END SUBROUTINE  STRAIN_XZ_2ND_PML

    SUBROUTINE      STRAIN_YY                                          &
                      (L, MXL, MXH,     MXT,                           &
                          MYL, MYH, PL, MYT, PYH,                      &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX, Y_MIN, Y_MAX,                  &
                          A, B, C, TPML,                               &
                          VTM, XYYT, EYYT, PYYLT, PYYHT,               &
                          ODA_YL, ODB_YL, ODC_YL, ODD_YL,              &
                          ODA_YH, ODB_YH, ODC_YH, ODD_YH )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                       INTENT (IN)    :: &
                                               MXL, MXH,     MXT,      &
                                               MYL, MYH, PL, MYT, PYH, &
                                               MZL, MZH,  TPML,   L,   &
                                            X_MIN, X_MAX, Y_MIN, Y_MAX
      REAL(PP),                                      INTENT (IN)    :: &
                                                               A, B, C
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EYYT
      REAL(PP),DIMENSION (MXL:MXH, PL :  1, MZL:MZH), INTENT(INOUT) :: &
                                                                 PYYLT
      REAL(PP),DIMENSION (MXL:MXH, MYT:PYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PYYHT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                   VTM
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XYYT
      REAL(PP),DIMENSION ( 0:TPML, MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
        ODA_YL, ODB_YL, ODC_YL, ODD_YL, ODA_YH, ODB_YH, ODC_YH, ODD_YH
    END SUBROUTINE  STRAIN_YY

    SUBROUTINE      STRAIN_YZ                                          &
                     ( L, MXL, MXH,     MXT,                           &
                          MYL, MYH, PL, MYT, PYH,                      &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX, Y_MIN, Y_MAX,                  &
                          A, B, C, TPML,                               &
                          VTM, WTM, XYZT, EYZT,                        &
                          PZYLT, PZYHT,                                &
                          ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL,          &
                          ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                       INTENT (IN)    :: &
                                               MXL, MXH,     MXT,      &
                                               MYL, MYH, PL, MYT, PYH, &
                                               MZL, MZH, L, TPML,      &
                                            X_MIN, X_MAX, Y_MIN, Y_MAX
      REAL(PP),                                      INTENT (IN)    :: &
                                                               A, B, C
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EYZT
      REAL(PP),DIMENSION (MXL:MXH, PL :  2, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZYLT
      REAL(PP),DIMENSION (MXL:MXH, MYT:PYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZYHT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                              VTM, WTM
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XYZT
      REAL(PP),DIMENSION(0:TPML+1, MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
                  ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL,                  &
                  ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH
    END SUBROUTINE  STRAIN_YZ

    SUBROUTINE      STRAIN_YZ_1                                        &
                     (    MXL, MXH,     MXT,                           &
                          MYL, MYH, PL, MYT, PYH,                      &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX, Y_MIN, Y_MAX,                  &
                          A, B, C, TPML,                               &
                          VTM, WTM, XYZT, EYZT,                        &
                          PZYLT, PZYHT,                                &
                          ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL,          &
                          ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                       INTENT (IN)    :: &
                                               MXL, MXH,     MXT,      &
                                               MYL, MYH, PL, MYT, PYH, &
                                               MZL, MZH,    TPML,      &
                                            X_MIN, X_MAX, Y_MIN, Y_MAX
      REAL(PP),                                      INTENT (IN)    :: &
                                                               A, B, C
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EYZT
      REAL(PP),DIMENSION (MXL:MXH, PL :  2, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZYLT
      REAL(PP),DIMENSION (MXL:MXH, MYT:PYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZYHT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                              VTM, WTM
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XYZT
      REAL(PP),DIMENSION(0:TPML+1, MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
                  ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL,                  &
                  ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH
    END SUBROUTINE  STRAIN_YZ_1

    SUBROUTINE      STRAIN_YZ_2ND                                      &
                         ( L, C, TPML, X_MIN, X_MAX, Y_MIN, Y_MAX,     &
                          MXL, MXH,     MXT,                           &
                          MYL, MYH, PL, MYT, PYH,                      &
                          MZL, MZH,                                    &
                       VTM, WTM, XYZT, EYZT,                           &
                          PZYLT, PZYHT,                                &
                          ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL,          &
                          ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                                               MXL, MXH, PL, MXT, PYH, &
                                               MYL, MYH,     MYT,      &
                                               MZL, MZH, L, TPML,      &
                                           X_MIN, X_MAX, Y_MIN, Y_MAX
      REAL(PP),                                       INTENT(IN)    :: &
                                                               C
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EYZT
      REAL(PP),DIMENSION (MXL:MXH, PL :  2, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZYLT
      REAL(PP),DIMENSION (MXL:MXH, MYT:PYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZYHT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                              VTM, WTM
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XYZT
      REAL(PP),DIMENSION(0:TPML+1, MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
                  ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL,                  &
                  ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH
    END SUBROUTINE  STRAIN_YZ_2ND

    SUBROUTINE      STRAIN_YZ_2ND_PML                                  &
                             ( L, LL, C, TPML,                         &
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
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                                               MXL, MXH, PL, MXT, PYH, &
                                               MYL, MYH,     MYT,      &
                                               MZL, MZH, L, LL, TPML,  &
                                           X_MIN, X_MAX, Y_MIN, Y_MAX
      REAL(PP),                                       INTENT(IN)    :: &
                                                               C
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EYZT
      REAL(PP),DIMENSION (MXL:MXH, PL :  2, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZYLT
      REAL(PP),DIMENSION (MXL:MXH, MYT:PYH, MZL:MZH), INTENT(INOUT) :: &
                                                                 PZYHT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH,0:TPML+1), INTENT(INOUT) :: &
                                                                 PYZHT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                              VTM, WTM
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XYZT
      REAL(PP),DIMENSION(0:TPML+1, MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
                  ODA2_YL, ODB2_YL, ODC2_YL, ODD2_YL,                  &
                  ODA2_YH, ODB2_YH, ODC2_YH, ODD2_YH
      REAL(PP),DIMENSION(0:TPML+1, MXL:MXH, MYL:MYH), INTENT(INOUT) :: &
                  ODA2_ZH, ODB2_ZH, ODC2_ZH, ODD2_ZH
    END SUBROUTINE  STRAIN_YZ_2ND_PML

    SUBROUTINE      STRAIN_ZZ                                          &
                     ( L, MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,      &
                          X_MIN, X_MAX, Y_MIN, Y_MAX,                  &
                          TPML, A, B, C,                               &
                          WTM, XZZT, EZZT )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                       INTENT (IN)    :: &
                                 MXL, MXH, MXT, MYL, MYH, MYT, TPML, L,&
                                 MZL, MZH,                             &
                                 X_MIN, X_MAX, Y_MIN, Y_MAX
      REAL(PP),                                      INTENT (IN)    :: &
                                                               A, B, C
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EZZT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                   WTM
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XZZT
    END SUBROUTINE  STRAIN_ZZ

    SUBROUTINE      STRAIN_ZZ_0                                        &
                        ( MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,      &
                          X_MIN, X_MAX, Y_MIN, Y_MAX, TPML, C,         &
                          WTM, XZZT, EZZT )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                       INTENT (IN)    :: &
                          MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH, TPML,&
                          X_MIN, X_MAX, Y_MIN, Y_MAX
      REAL(PP),                                      INTENT (IN)    :: &
                                                                     C
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EZZT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                   WTM
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XZZT
    END SUBROUTINE  STRAIN_ZZ_0

    SUBROUTINE      STRAIN_ZZ_2ND                                      &
                         ( L, C, TPML, X_MIN, X_MAX, Y_MIN, Y_MAX,     &
                          MXL, MXH,     MXT,                           &
                          MYL, MYH,     MYT,                           &
                          MZL, MZH,                                    &
                          WTM, XZZT, EZZT )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                                               MXL, MXH,     MXT,      &
                                               MYL, MYH,     MYT,      &
                                               MZL, MZH, L,     TPML,  &
                                           X_MIN, X_MAX, Y_MIN, Y_MAX
      REAL(PP),                                       INTENT(IN)    :: &
                                                               C
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EZZT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                   WTM
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XZZT
    END SUBROUTINE  STRAIN_ZZ_2ND

    SUBROUTINE      STRAIN_ZZ_2ND_PML                                  &
                             ( L, LL, C, TPML,                         &
                          MXL, MXH,     MXT,                           &
                          MYL, MYH,     MYT,                           &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX, Y_MIN, Y_MAX,                  &
                          WTM, XZZT, EZZT, PZZHT,                      &
                          ODA_ZH, ODB_ZH, ODC_ZH, ODD_ZH )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                                               MXL, MXH,     MXT,      &
                                               MYL, MYH,     MYT,      &
                                               MZL, MZH, L, LL, TPML,  &
                                            X_MIN, X_MAX, Y_MIN, Y_MAX
      REAL(PP),                                       INTENT(IN)    :: &
                                                                     C
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EZZT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                   WTM
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XZZT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH,  0:TPML), INTENT(INOUT) :: &
                                                                 PZZHT
      REAL(PP),DIMENSION ( 0:TPML, MXL:MXH, MYL:MYH), INTENT(INOUT) :: &
                                        ODA_ZH, ODB_ZH, ODC_ZH, ODD_ZH
    END SUBROUTINE  STRAIN_ZZ_2ND_PML

    SUBROUTINE      STRESS_NN                                          &
                      (L    , MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,  &
                     X_MIN, X_MAX, Y_MIN, Y_MAX,                       &
                     TPML, TXXT,       TYYT,       TZZT,               &
                           XXXT, EXXT, XYYT, EYYT, XZZT, EZZT, JMT,    &
                           XXXK2T, XXXI2T, XXXKMY1T, XXXIMX1T,         &
                           XYYK2T, XYYI2T, XYYKMY1T, XYYIMX1T,         &
                           XZZK2T, XZZI2T, XZZKMY1T, XZZIMX1T )
      USE PRECISION   , ONLY: PP, IK
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                       MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH, L, TPML,&
                       X_MIN, X_MAX, Y_MIN, Y_MAX
      INTEGER(IK),DIMENSION(MXL:MXH, MYL:MYH, MZL:MZH),INTENT(INOUT):: &
                                                                   JMT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                      EXXT, EYYT, EZZT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                      XXXT, XYYT, XZZT
      REAL(PP),DIMENSION (MXL:MXH,       4, MZL:MZH), INTENT(INOUT) :: &
                  XXXK2T, XXXKMY1T, XYYK2T, XYYKMY1T, XZZK2T, XZZKMY1T
      REAL(PP),DIMENSION (      MYL:MYH, 4, MZL:MZH), INTENT(INOUT) :: &
                  XXXI2T, XXXIMX1T, XYYI2T, XYYIMX1T, XZZI2T, XZZIMX1T
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                      TXXT, TYYT, TZZT
    END SUBROUTINE  STRESS_NN

    SUBROUTINE      STRESS_XY                                          &
                      (L    , MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,  &
                                          X_MIN, X_MAX, Y_MIN, Y_MAX,  &
                                          TPML, TXYT, EXYT, XXYT, JMT, &
                           XXYK2T, XXYKMYT, XXYI2T, XXYIMXT )
      USE PRECISION   , ONLY: PP, IK
      IMPLICIT NONE
      INTEGER,                                       INTENT (IN)    :: &
                       MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH, L, TPML,&
                       X_MIN, X_MAX, Y_MIN, Y_MAX
      INTEGER(IK),DIMENSION(MXL:MXH, MYL:MYH, MZL:MZH),INTENT(INOUT):: &
                                                                   JMT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EXYT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XXYT
      REAL(PP),DIMENSION (MXL:MXH,       4, MZL:MZH), INTENT(INOUT) :: &
                                                      XXYK2T, XXYKMYT
      REAL(PP),DIMENSION (      MYL:MYH, 4, MZL:MZH), INTENT(INOUT) :: &
                                                      XXYI2T, XXYIMXT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  TXYT
    END SUBROUTINE  STRESS_XY

    SUBROUTINE      STRESS_XZ                                          &
                      (L    , MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,  &
                                          X_MIN, X_MAX, Y_MIN, Y_MAX,  &
                                          TPML, TXZT, EXZT, XXZT, JMT, &
                           XXZK2T, XXZKMY1T, XXZI2T, XXZIMXT )
      USE PRECISION   , ONLY: PP, IK
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                       MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH, L, TPML,&
                       X_MIN, X_MAX, Y_MIN, Y_MAX
      INTEGER(IK),DIMENSION(MXL:MXH, MYL:MYH, MZL:MZH),INTENT(INOUT):: &
                                                                   JMT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EXZT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XXZT
      REAL(PP),DIMENSION (MXL:MXH,       4, MZL:MZH), INTENT(INOUT) :: &
                                                      XXZK2T, XXZKMY1T
      REAL(PP),DIMENSION (      MYL:MYH, 4, MZL:MZH), INTENT(INOUT) :: &
                                                      XXZI2T, XXZIMXT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  TXZT
    END SUBROUTINE  STRESS_XZ

    SUBROUTINE      STRESS_YZ                                          &
                      (L    , MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH,  &
                                          X_MIN, X_MAX, Y_MIN, Y_MAX,  &
                                          TPML, TYZT, EYZT, XYZT, JMT, &
                           XYZK2T, XYZKMYT, XYZI2T, XYZIMX1T )
      USE PRECISION   , ONLY: PP, IK
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                       MXL, MXH, MXT, MYL, MYH, MYT, MZL, MZH, L, TPML,&
                       X_MIN, X_MAX, Y_MIN, Y_MAX
      INTEGER(IK),DIMENSION(MXL:MXH, MYL:MYH, MZL:MZH),INTENT(INOUT):: &
                                                                   JMT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH         ), INTENT(INOUT) :: &
                                                                  EYZT
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  XYZT
      REAL(PP),DIMENSION (MXL:MXH,       4, MZL:MZH), INTENT(INOUT) :: &
                                                      XYZK2T, XYZKMYT
      REAL(PP),DIMENSION (      MYL:MYH, 4, MZL:MZH), INTENT(INOUT) :: &
                                                      XYZI2T, XYZIMX1T
      REAL(PP),DIMENSION (MXL:MXH, MYL:MYH, MZL:MZH), INTENT(INOUT) :: &
                                                                  TYZT
    END SUBROUTINE  STRESS_YZ

    SUBROUTINE      TIME_LOOP
    END SUBROUTINE  TIME_LOOP

  END INTERFACE

END MODULE  INTERFACES

!=======================================================================
! Copyright (c) 2014 Jozef Kristek, Peter Moczo, Martin Galis
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
! Version 1.0, 2014-03-19: Jozef Kristek, Peter Moczo, Martin Galis
!=======================================================================
!-----------------------------------------------------------------------
! function returns values of YMU and YL2M as function of coordinates x,y,z
!-----------------------------------------------------------------------
FUNCTION FUNC_YMU_YL2M (XA,YA,ZA)

  USE NRTYPE    , ONLY: WP
  USE MOD_CONST , ONLY: N_FREQ
  USE MOD_MODEL , ONLY: LAYER_1, LAYER_2
  USE MOD_ATTEN , ONLY: PLO, Q, Y
  USE MOD_PARAMETERIZATION                                             &
                , ONLY: YMU_HOM, YL2M_HOM, QS_OLD, QP_OLD, YMU_YL2M_OLD
  USE MOD_INTERFACES, ONLY: GET_LAYER
  USE MOD_FUNC  , ONLY: FUNCQS, FUNCQS_1D, FUNCQP, FUNCQP_1D

!----------------------------------------------------------------------

  IMPLICIT NONE

  REAL(WP), INTENT(IN)          :: XA,YA,ZA
  REAL(WP), DIMENSION(N_FREQ,2) :: FUNC_YMU_YL2M

  INTEGER                         :: KF, LAYER
  REAL(WP)                        :: FQS, FQP
  REAL(WP), DIMENSION(2*N_FREQ-1) :: HEQSW, HEQPW

!----------------------------------------------------------------------

  LAYER = GET_LAYER(XA,YA,ZA)

  SELECT CASE (LAYER)
  CASE(LAYER_1)  ! in sediments (1st layer)

    FQS = FUNCQS(XA,YA,ZA)
    FQP = FUNCQP(XA,YA,ZA)

    IF  ( FQS == QS_OLD ) THEN
      IF ( FQP == QP_OLD ) THEN
        FUNC_YMU_YL2M = YMU_YL2M_OLD
        RETURN
      END IF
    END IF

    DO KF = 1, 2*N_FREQ-1
      HEQSW(KF) = FQS - PLO(KF)
      HEQPW(KF) = FQP - PLO(KF)
    END DO

    Q   = HEQSW
    CALL QYLSM(Q,Y)

    FUNC_YMU_YL2M (:,1) = Y
    QS_OLD = FQS
    YMU_YL2M_OLD  (:,1) = Y

    Q   = HEQPW
    CALL QYLSM(Q,Y)

    FUNC_YMU_YL2M (:,2) = Y
    QP_OLD = FQP
    YMU_YL2M_OLD  (:,2) = Y

    RETURN
  CASE(LAYER_2)  ! in bedrock   (2nd layer)

    FQS = FUNCQS_1D(ZA)
    FQP = FUNCQP_1D(ZA)

    IF  ( FQS == QS_OLD ) THEN
      IF ( FQP == QP_OLD ) THEN
        FUNC_YMU_YL2M = YMU_YL2M_OLD
        RETURN
      END IF
    END IF

    DO KF = 1, 2*N_FREQ-1
      HEQSW(KF) = FQS - PLO(KF)
      HEQPW(KF) = FQP - PLO(KF)
    END DO

    Q   = HEQSW
    CALL QYLSM(Q,Y)

    FUNC_YMU_YL2M (:,1) = Y
    QS_OLD = FQS
    YMU_YL2M_OLD  (:,1) = Y

    Q   = HEQPW
    CALL QYLSM(Q,Y)

    FUNC_YMU_YL2M (:,2) = Y
    QP_OLD = FQP
    YMU_YL2M_OLD  (:,2) = Y

    RETURN
  END SELECT

END FUNCTION FUNC_YMU_YL2M

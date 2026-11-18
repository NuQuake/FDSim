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
! this subroutine calculates Y coefficients in every layer in model
!-----------------------------------------------------------------------
SUBROUTINE PREPARE_Y

  USE NRTYPE              , ONLY: WP
  USE MOD_FILES
  USE MOD_CONST           , ONLY: N_FREQ, I2N_FREQ1
  USE MOD_ATTEN           , ONLY: Q, Y, FJK, FJJ, QYLSM, PLO, FRJ, FRK,&
                                  CF, PFP21, PF21
  USE MOD_PARAMETERIZATION, ONLY: QS_HOM, QP_HOM, YMU_HOM, YL2M_HOM,   &
                                  FRANGE, FRJMAX, FREF,                &
                                  QS_OLD, QP_OLD, YMU_YL2M_OLD
!----------------------------------------------------------------------

  IMPLICIT NONE

  REAL(WP), PARAMETER  :: PI=3.141592654

  INTEGER    :: NL, I, M, KF, J, K
  REAL(WP)   :: L2M, FQS, FQP, PF
  REAL(WP)   :: RN_FREQ, R2N_FREQ
  REAL(WP), DIMENSION(  N_FREQ ) :: YL2M
  REAL(WP), DIMENSION(I2N_FREQ1) :: HEQSW, HEQPW

!----------------------------------------------------------------------

!  WRITE(*,'(A)')'Preparing Y coefficients...'

  !-------------------------------------- code for calculating Y coef.
  RN_FREQ   = REAL(N_FREQ)
  R2N_FREQ  = 2*RN_FREQ

  FREF = FREF*2.*PI

  DO I = 1, N_FREQ
    FRJ(I)=FRJMAX/(10.**(FRANGE*REAL(N_FREQ-I)/REAL(N_FREQ-1)))*2.*PI
  END DO

  DO I = 1, I2N_FREQ1
    FRK(I)=FRJMAX/(10.**(FRANGE*REAL(I2N_FREQ1-I)/REAL(I2N_FREQ1)))*2.*PI
  END DO

  DO KF = 1, I2N_FREQ1
!    PLO(KF)= LOG (FRK(KF)/FREF) / PI !FUTTERMAN Q(w) LAW
    PLO(KF)= 0.                       !CONSTANT Q(w) LAW
  END DO

  DO J = 1, N_FREQ
    PF      = FREF / FRJ(J)
    PF21 (J) = PF*PF + 1._WP
    PFP21(J) = PF / PF21(J)
    DO K = 1, I2N_FREQ1
      FJK(J,K) = FRJ(J)*FRK(K) / ( FRJ(J)*FRJ(J) + FRK(K)*FRK(K) )
      FJJ(J,K) = FRJ(J)*FRJ(J) / ( FRJ(J)*FRJ(J) + FRK(K)*FRK(K) )
      CF (J,K) = FRJ(J) / ( CMPLX(0.,FRK(K)) + FRJ(J) )
    END DO
  END DO

  !---------------------------------------------------- in every layer
  FQS = QS_HOM
  FQP = QP_HOM

  DO KF = 1, I2N_FREQ1
    HEQSW(KF) = FQS - PLO(KF)
    HEQPW(KF) = FQP - PLO(KF)
  END DO

  Q   = HEQSW

  CALL QYLSM (Q,Y)
  YMU_HOM(:) = Y

  Q   = HEQPW
  CALL QYLSM (Q,Y)
  YL2M_HOM(:) = Y

  QS_OLD = QS_HOM
  QP_OLD = QP_HOM
  YMU_YL2M_OLD(:,1) =  YMU_HOM
  YMU_YL2M_OLD(:,2) = YL2M_HOM


END SUBROUTINE PREPARE_Y


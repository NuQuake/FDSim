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

MODULE MOD_EFFMATPAR_Y

CONTAINS

FUNCTION MUXZ_Y (II,KI,LI,YMU)

  USE NRTYPE         , ONLY: WP
  USE MOD_INTEGRATION, ONLY: PTS_1, PTS_2, PTS2
  USE MOD_CONST      , ONLY: I2N_FREQ1

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: II, KI, LI
  REAL(WP), DIMENSION(I2N_FREQ1) :: MUXZ_Y

  COMPLEX(WP), DIMENSION(:,:,:,:), INTENT(INOUT) :: YMU

  INTEGER              :: I, K, L
  COMPLEX(WP), DIMENSION(I2N_FREQ1) :: SUM1, SUM2

!----------------------------------------------------------------------

  SUM2 = CMPLX(0.)

  DO L = -PTS2+1, PTS2
    DO I = -PTS2+1, PTS2
      SUM1 = CMPLX(0.)
      DO K = -PTS2+1, PTS2
        SUM1(:) = SUM1(:) + YMU(II+I, KI+K, LI+L,:)
      END DO
      SUM2 = SUM2 + PTS_1 / SUM1
    END DO
  END DO

  SUM2 = PTS_2 / SUM2

  MUXZ_Y = REAL(SUM2)/AIMAG(SUM2)

END FUNCTION MUXZ_Y

!======================================================================

FUNCTION MUYZ_Y (II,KI,LI,YMU)

  USE NRTYPE         , ONLY: WP
  USE MOD_INTEGRATION, ONLY: PTS_1, PTS_2, PTS2
  USE MOD_CONST      , ONLY: I2N_FREQ1

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: II, KI, LI
  REAL(WP), DIMENSION(I2N_FREQ1) :: MUYZ_Y

  COMPLEX(WP), DIMENSION(:,:,:,:), INTENT(INOUT) :: YMU

  INTEGER              :: I, K, L
  COMPLEX(WP), DIMENSION(I2N_FREQ1) :: SUM1, SUM2

!----------------------------------------------------------------------

  SUM2 = CMPLX(0.)

  DO L = -PTS2+1, PTS2
    DO K = -PTS2+1, PTS2
      SUM1 = CMPLX(0.)
      DO I = -PTS2+1, PTS2
        SUM1(:) = SUM1(:) + YMU(II+I, KI+K, LI+L,:)
      END DO
      SUM2 = SUM2 + PTS_1 / SUM1
    END DO
  END DO

  SUM2 = PTS_2 / SUM2

  MUYZ_Y = REAL(SUM2)/AIMAG(SUM2)

END FUNCTION MUYZ_Y

!======================================================================

FUNCTION MUXY_Y (II,KI,LI,YMU)

  USE NRTYPE         , ONLY: WP
  USE MOD_INTEGRATION, ONLY: PTS_1, PTS_2, PTS2
  USE MOD_CONST      , ONLY: I2N_FREQ1

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: II, KI, LI
  REAL(WP), DIMENSION(I2N_FREQ1) :: MUXY_Y

  COMPLEX(WP), DIMENSION(:,:,:,:), INTENT(INOUT) :: YMU

  INTEGER              :: I, K, L
  COMPLEX(WP), DIMENSION(I2N_FREQ1) :: SUM1, SUM2

!----------------------------------------------------------------------

  SUM2 = CMPLX(0.)

  DO K = -PTS2+1, PTS2
    DO I = -PTS2+1, PTS2
      SUM1 = CMPLX(0.)
      DO L = -PTS2+1, PTS2
        SUM1(:) = SUM1(:) + YMU(II+I, KI+K, LI+L,:)
      END DO
      SUM2 = SUM2 + PTS_1 / SUM1
    END DO
  END DO

  SUM2 = PTS_2 / SUM2

  MUXY_Y = REAL(SUM2)/AIMAG(SUM2)

END FUNCTION MUXY_Y

!======================================================================

FUNCTION L2MX_Y (II,KI,LI, YLAM, YMU)

  USE NRTYPE         , ONLY: WP
  USE MOD_INTEGRATION, ONLY: PTS_1, PTS_2, PTS2
  USE MOD_CONST      , ONLY: I2N_FREQ1

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: II, KI, LI
  REAL(WP), DIMENSION(I2N_FREQ1) :: L2MX_Y

  COMPLEX(WP), DIMENSION(:,:,:,:), INTENT(INOUT) :: YMU, YLAM

  INTEGER              :: I
  COMPLEX(WP), DIMENSION(I2N_FREQ1) :: SUM1

!----------------------------------------------------------------------

  SUM1 = CMPLX(0.)
  DO I = -PTS2+1, PTS2
    SUM1 = SUM1 + CMPLX(1.)/L2M_YZ()
  END DO

  SUM1 = PTS_1/SUM1

  L2MX_Y = REAL(SUM1)/AIMAG(SUM1)

CONTAINS

  FUNCTION L2M_YZ()

    COMPLEX(WP), DIMENSION(I2N_FREQ1) :: L2M_YZ

    INTEGER              :: K, L
    COMPLEX(WP), DIMENSION(I2N_FREQ1) :: S1, S2, S3, LAM2MU, LAMBDA

  !----------------------------------------------------------------------

    S1 = CMPLX(0.)
    S2 = CMPLX(0.)
    S3 = CMPLX(0.)
    DO L = -PTS2+1, PTS2
      DO K = -PTS2+1, PTS2
        LAMBDA = YLAM(II+I, KI+K, LI+L,:)
        LAM2MU = LAMBDA + 2. * YMU(II+I, KI+K, LI+L,:)
        S1 = S1 + LAMBDA / LAM2MU
        S2 = S2 +      1./ LAM2MU
        S3 = S3 + LAM2MU - LAMBDA*LAMBDA/LAM2MU
      END DO
    END DO

    S1 = S1/PTS_2
    S2 = S2/PTS_2
    S3 = S3/PTS_2

    L2M_YZ = S1*S1/S2+S3

  END FUNCTION L2M_YZ


END FUNCTION L2MX_Y

!======================================================================

FUNCTION L2MY_Y (II,KI,LI, YLAM, YMU)

  USE NRTYPE         , ONLY: WP
  USE MOD_INTEGRATION, ONLY: PTS_1, PTS_2, PTS2
  USE MOD_CONST      , ONLY: I2N_FREQ1

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: II, KI, LI
  REAL(WP), DIMENSION(I2N_FREQ1) :: L2MY_Y

  COMPLEX(WP), DIMENSION(:,:,:,:), INTENT(INOUT) :: YMU, YLAM

  INTEGER              :: K
  COMPLEX(WP), DIMENSION(I2N_FREQ1) :: SUM1

!----------------------------------------------------------------------

  SUM1 = CMPLX(0.)
  DO K = -PTS2+1, PTS2
    SUM1 = SUM1 + CMPLX(1.)/L2M_XZ()
  END DO

  SUM1 = PTS_1/SUM1

  L2MY_Y = REAL(SUM1)/AIMAG(SUM1)

CONTAINS

  FUNCTION L2M_XZ()

    COMPLEX(WP), DIMENSION(I2N_FREQ1) :: L2M_XZ

    INTEGER              :: I, L
    COMPLEX(WP), DIMENSION(I2N_FREQ1) :: S1, S2, S3, LAM2MU, LAMBDA

  !----------------------------------------------------------------------

    S1 = CMPLX(0.)
    S2 = CMPLX(0.)
    S3 = CMPLX(0.)
    DO L = -PTS2+1, PTS2
      DO I = -PTS2+1, PTS2
        LAMBDA = YLAM(II+I, KI+K, LI+L,:)
        LAM2MU = LAMBDA + 2. * YMU(II+I, KI+K, LI+L,:)
        S1 = S1 + LAMBDA / LAM2MU
        S2 = S2 +      1./ LAM2MU
        S3 = S3 + LAM2MU - LAMBDA*LAMBDA/LAM2MU
      END DO
    END DO

    S1 = S1/PTS_2
    S2 = S2/PTS_2
    S3 = S3/PTS_2

    L2M_XZ = S1*S1/S2+S3

  END FUNCTION L2M_XZ


END FUNCTION L2MY_Y

!======================================================================

FUNCTION L2MZ_Y (II,KI,LI, YLAM, YMU)

  USE NRTYPE         , ONLY: WP
  USE MOD_INTEGRATION, ONLY: PTS_1, PTS_2, PTS2
  USE MOD_CONST      , ONLY: I2N_FREQ1

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: II, KI, LI
  REAL(WP), DIMENSION(I2N_FREQ1) :: L2MZ_Y

  COMPLEX(WP), DIMENSION(:,:,:,:), INTENT(INOUT) :: YMU, YLAM

  INTEGER              :: L
  COMPLEX(WP), DIMENSION(I2N_FREQ1) :: SUM1

!----------------------------------------------------------------------

  SUM1 = CMPLX(0.)
  DO L = -PTS2+1, PTS2
    SUM1 = SUM1 + CMPLX(1.)/L2M_XY()
  END DO

  SUM1 = PTS_1/SUM1

  L2MZ_Y = REAL(SUM1)/AIMAG(SUM1)

CONTAINS

  FUNCTION L2M_XY()

    COMPLEX(WP), DIMENSION(I2N_FREQ1) :: L2M_XY

    INTEGER              :: I, K
    COMPLEX(WP), DIMENSION(I2N_FREQ1) :: S1, S2, S3, LAM2MU, LAMBDA

  !----------------------------------------------------------------------

    S1 = CMPLX(0.)
    S2 = CMPLX(0.)
    S3 = CMPLX(0.)
    DO K = -PTS2+1, PTS2
      DO I = -PTS2+1, PTS2
        LAMBDA = YLAM(II+I, KI+K, LI+L,:)
        LAM2MU = LAMBDA + 2. * YMU(II+I, KI+K, LI+L,:)
        S1 = S1 + LAMBDA / LAM2MU
        S2 = S2 +      1./ LAM2MU
        S3 = S3 + LAM2MU - LAMBDA*LAMBDA/LAM2MU
      END DO
    END DO

    S1 = S1/PTS_2
    S2 = S2/PTS_2
    S3 = S3/PTS_2

    L2M_XY = S1*S1/S2+S3

  END FUNCTION L2M_XY


END FUNCTION L2MZ_Y

!======================================================================

FUNCTION LAXZ_Y (II,KI,LI, YLAM, YMU)

  USE NRTYPE         , ONLY: WP
  USE MOD_INTEGRATION, ONLY: PTS_1, PTS2
  USE MOD_CONST      , ONLY: I2N_FREQ1

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: II, KI, LI
  REAL(WP), DIMENSION(I2N_FREQ1) :: LAXZ_Y

  COMPLEX(WP), DIMENSION(:,:,:,:), INTENT(INOUT) :: YMU, YLAM

  INTEGER              :: I, L
  COMPLEX(WP), DIMENSION(I2N_FREQ1) :: SUM1, SUM2, L2M

!----------------------------------------------------------------------

  SUM1 = CMPLX(0.)
  SUM2 = CMPLX(0.)

  DO L = -PTS2+1, PTS2
    DO I = -PTS2+1, PTS2
      L2M  = L2M_Y()
      SUM1 = SUM1 + (L2M-2.*MU_Y())/L2M
      SUM2 = SUM2 +              1./L2M
    END DO
  END DO

  SUM1 = SUM1/SUM2

  LAXZ_Y = REAL(SUM1)/AIMAG(SUM1)

CONTAINS

  FUNCTION L2M_Y()

    COMPLEX(WP), DIMENSION(I2N_FREQ1) :: L2M_Y
    INTEGER                           :: K
    COMPLEX(WP), DIMENSION(I2N_FREQ1) :: S1, S2, S3, LAM2MU, LAMBDA

    S1 = CMPLX(0.)
    S2 = CMPLX(0.)
    S3 = CMPLX(0.)
    DO K = -PTS2+1, PTS2
      LAMBDA = YLAM(II+I, KI+K, LI+L,:)
      LAM2MU = LAMBDA + 2. * YMU(II+I, KI+K, LI+L,:)
      S1 = S1 + LAMBDA / LAM2MU
      S2 = S2 +      1./ LAM2MU
      S3 = S3 + LAM2MU - LAMBDA*LAMBDA/LAM2MU
    END DO

    S1 = S1/PTS_1
    S2 = S2/PTS_1
    S3 = S3/PTS_1

    L2M_Y = S1*S1/S2+S3

  END FUNCTION L2M_Y

  FUNCTION MU_Y()

    COMPLEX(WP), DIMENSION(I2N_FREQ1) :: MU_Y
    INTEGER                           :: K
    COMPLEX(WP), DIMENSION(I2N_FREQ1) :: S1

    S1 = CMPLX(0.)
    DO K = -PTS2+1, PTS2
      S1 = S1 + YMU(II+I, KI+K, LI+L,:)
    END DO

    MU_Y = S1/PTS_1

  END FUNCTION MU_Y

END FUNCTION LAXZ_Y

!======================================================================

FUNCTION LAYZ_Y (II,KI,LI, YLAM, YMU)

  USE NRTYPE         , ONLY: WP
  USE MOD_INTEGRATION, ONLY: PTS_1, PTS2
  USE MOD_CONST      , ONLY: I2N_FREQ1

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: II, KI, LI
  REAL(WP), DIMENSION(I2N_FREQ1) :: LAYZ_Y

  COMPLEX(WP), DIMENSION(:,:,:,:), INTENT(INOUT) :: YMU, YLAM

  INTEGER              :: K, L
  COMPLEX(WP), DIMENSION(I2N_FREQ1) :: SUM1, SUM2, L2M

!----------------------------------------------------------------------

  SUM1 = CMPLX(0.)
  SUM2 = CMPLX(0.)

  DO L = -PTS2+1, PTS2
    DO K = -PTS2+1, PTS2
      L2M  = L2M_X()
      SUM1 = SUM1 + (L2M-2.*MU_X())/L2M
      SUM2 = SUM2 +              1./L2M
    END DO
  END DO

  SUM1 = SUM1/SUM2

  LAYZ_Y = REAL(SUM1)/AIMAG(SUM1)

CONTAINS

  FUNCTION L2M_X()

    COMPLEX(WP), DIMENSION(I2N_FREQ1) :: L2M_X
    INTEGER                           :: I
    COMPLEX(WP), DIMENSION(I2N_FREQ1) :: S1, S2, S3, LAM2MU, LAMBDA

    S1 = CMPLX(0.)
    S2 = CMPLX(0.)
    S3 = CMPLX(0.)
    DO I = -PTS2+1, PTS2
      LAMBDA = YLAM(II+I, KI+K, LI+L,:)
      LAM2MU = LAMBDA + 2. * YMU(II+I, KI+K, LI+L,:)
      S1 = S1 + LAMBDA / LAM2MU
      S2 = S2 +      1./ LAM2MU
      S3 = S3 + LAM2MU - LAMBDA*LAMBDA/LAM2MU
    END DO

    S1 = S1/PTS_1
    S2 = S2/PTS_1
    S3 = S3/PTS_1

    L2M_X = S1*S1/S2+S3

  END FUNCTION L2M_X

  FUNCTION MU_X()

    COMPLEX(WP), DIMENSION(I2N_FREQ1) :: MU_X
    INTEGER                           :: I
    COMPLEX(WP), DIMENSION(I2N_FREQ1) :: S1

    S1 = CMPLX(0.)
    DO I = -PTS2+1, PTS2
      S1 = S1 + YMU(II+I, KI+K, LI+L,:)
    END DO

    MU_X = S1/PTS_1

  END FUNCTION MU_X

END FUNCTION LAYZ_Y

!======================================================================

FUNCTION LAXY_Y (II,KI,LI, YLAM, YMU)

  USE NRTYPE         , ONLY: WP
  USE MOD_INTEGRATION, ONLY: PTS_1, PTS2
  USE MOD_CONST      , ONLY: I2N_FREQ1

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: II, KI, LI
  REAL(WP), DIMENSION(I2N_FREQ1) ::LAXY_Y

  COMPLEX(WP), DIMENSION(:,:,:,:), INTENT(INOUT) :: YMU, YLAM

  INTEGER              :: I, K
  COMPLEX(WP), DIMENSION(I2N_FREQ1) :: SUM1, SUM2, L2M

!----------------------------------------------------------------------

  SUM1 = CMPLX(0.)
  SUM2 = CMPLX(0.)

  DO K = -PTS2+1, PTS2
    DO I = -PTS2+1, PTS2
      L2M  = L2M_Z()
      SUM1 = SUM1 + (L2M-2.*MU_Z())/L2M
      SUM2 = SUM2 +              1./L2M
    END DO
  END DO

  SUM1 = SUM1/SUM2

  LAXY_Y = REAL(SUM1)/AIMAG(SUM1)

CONTAINS

  FUNCTION L2M_Z()

    COMPLEX(WP), DIMENSION(I2N_FREQ1) :: L2M_Z
    INTEGER                           :: L
    COMPLEX(WP), DIMENSION(I2N_FREQ1) :: S1, S2, S3, LAM2MU, LAMBDA

    S1 = CMPLX(0.)
    S2 = CMPLX(0.)
    S3 = CMPLX(0.)
    DO L = -PTS2+1, PTS2
      LAMBDA = YLAM(II+I, KI+K, LI+L,:)
      LAM2MU = LAMBDA + 2. * YMU(II+I, KI+K, LI+L,:)
      S1 = S1 + LAMBDA / LAM2MU
      S2 = S2 +      1./ LAM2MU
      S3 = S3 + LAM2MU - LAMBDA*LAMBDA/LAM2MU
    END DO

    S1 = S1/PTS_1
    S2 = S2/PTS_1
    S3 = S3/PTS_1

    L2M_Z = S1*S1/S2+S3

  END FUNCTION L2M_Z

  FUNCTION MU_Z()

    COMPLEX(WP), DIMENSION(I2N_FREQ1) :: MU_Z
    INTEGER                           :: L
    COMPLEX(WP), DIMENSION(I2N_FREQ1) :: S1

    S1 = 0.
    DO L = -PTS2+1, PTS2
      S1 = S1 + YMU(II+I, KI+K, LI+L,:)
    END DO

    MU_Z = S1/PTS_1

  END FUNCTION MU_Z

END FUNCTION LAXY_Y


END MODULE MOD_EFFMATPAR_Y

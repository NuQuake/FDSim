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

MODULE MOD_EFFMATPAR

CONTAINS

FUNCTION MUXZ (II,KI,LI,MU)

  USE NRTYPE         , ONLY: WP
  USE MOD_INTEGRATION, ONLY: PTS_1, PTS_2, PTS2

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: II, KI, LI
  REAL(WP)             :: MUXZ

  REAL(WP), DIMENSION(:,:,:), INTENT(INOUT) :: MU

  INTEGER              :: I, K, L
  REAL(WP)             :: SUM1, SUM2

!----------------------------------------------------------------------

  SUM2 = 0.

  DO L = -PTS2+1, PTS2
    DO I = -PTS2+1, PTS2
      SUM1 = 0.
      DO K = -PTS2+1, PTS2
        SUM1 = SUM1 + MU(II+I, KI+K, LI+L)
      END DO
      SUM2 = SUM2 + PTS_1 / SUM1
    END DO
  END DO

  MUXZ = PTS_2/SUM2

END FUNCTION MUXZ

!======================================================================

FUNCTION MUYZ (II,KI,LI,MU)

  USE NRTYPE         , ONLY: WP
  USE MOD_INTEGRATION, ONLY: PTS_1, PTS_2, PTS2

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: II, KI, LI
  REAL(WP)             :: MUYZ

  REAL(WP), DIMENSION(:,:,:), INTENT(INOUT) :: MU

  INTEGER              :: I, K, L
  REAL(WP)             :: SUM1, SUM2

!----------------------------------------------------------------------

  SUM2 = 0.

  DO L = -PTS2+1, PTS2
    DO K = -PTS2+1, PTS2
      SUM1 = 0.
      DO I = -PTS2+1, PTS2
        SUM1 = SUM1 + MU(II+I, KI+K, LI+L)
      END DO
      SUM2 = SUM2 + PTS_1 / SUM1
    END DO
  END DO

  MUYZ = PTS_2/SUM2

END FUNCTION MUYZ

!======================================================================

FUNCTION MUXY (II,KI,LI,MU)

  USE NRTYPE         , ONLY: WP
  USE MOD_INTEGRATION, ONLY: PTS_1, PTS_2, PTS2

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: II, KI, LI
  REAL(WP)             :: MUXY

  REAL(WP), DIMENSION(:,:,:), INTENT(INOUT) :: MU

  INTEGER              :: I, K, L
  REAL(WP)             :: SUM1, SUM2

!----------------------------------------------------------------------

  SUM2 = 0.

  DO K = -PTS2+1, PTS2
    DO I = -PTS2+1, PTS2
      SUM1 = 0.
      DO L = -PTS2+1, PTS2
        SUM1 = SUM1 + MU(II+I, KI+K, LI+L)
      END DO
      SUM2 = SUM2 + PTS_1 / SUM1
    END DO
  END DO

  MUXY = PTS_2/SUM2

END FUNCTION MUXY

!======================================================================

FUNCTION L2MX (II,KI,LI, LAM, MU)

  USE NRTYPE         , ONLY: WP
  USE MOD_INTEGRATION, ONLY: PTS_1, PTS_2, PTS2

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: II, KI, LI
  REAL(WP)             :: L2MX

  REAL(WP), DIMENSION(:,:,:), INTENT(INOUT) :: MU, LAM

  INTEGER              :: I
  REAL(WP)             :: SUM1

!----------------------------------------------------------------------

  SUM1 = 0.
  DO I = -PTS2+1, PTS2
    SUM1 = SUM1 + 1./L2M_YZ()
  END DO

  L2MX = PTS_1/SUM1

CONTAINS

  FUNCTION L2M_YZ()

    REAL(WP)             :: L2M_YZ

    INTEGER              :: K, L
    REAL(WP)             :: S1, S2, S3, LAM2MU, LAMBDA

  !----------------------------------------------------------------------

    S1 = 0.
    S2 = 0.
    S3 = 0.
    DO L = -PTS2+1, PTS2
      DO K = -PTS2+1, PTS2
        LAMBDA = LAM(II+I, KI+K, LI+L)
        LAM2MU = LAMBDA + 2. * MU(II+I, KI+K, LI+L)
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


END FUNCTION L2MX

!======================================================================

FUNCTION L2MY (II,KI,LI, LAM, MU)

  USE NRTYPE         , ONLY: WP
  USE MOD_INTEGRATION, ONLY: PTS_1, PTS_2, PTS2

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: II, KI, LI
  REAL(WP)             :: L2MY

  REAL(WP), DIMENSION(:,:,:), INTENT(INOUT) :: MU, LAM

  INTEGER              :: K
  REAL(WP)             :: SUM1

!----------------------------------------------------------------------

  SUM1 = 0.
  DO K = -PTS2+1, PTS2
    SUM1 = SUM1 + 1./L2M_XZ()
  END DO

  L2MY = PTS_1/SUM1

CONTAINS

  FUNCTION L2M_XZ()

    REAL(WP)             :: L2M_XZ

    INTEGER              :: I, L
    REAL(WP)             :: S1, S2, S3, LAM2MU, LAMBDA

  !----------------------------------------------------------------------

    S1 = 0.
    S2 = 0.
    S3 = 0.
    DO L = -PTS2+1, PTS2
      DO I = -PTS2+1, PTS2
        LAMBDA = LAM(II+I, KI+K, LI+L)
        LAM2MU = LAMBDA + 2. * MU(II+I, KI+K, LI+L)
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


END FUNCTION L2MY

!======================================================================

FUNCTION L2MZ (II,KI,LI, LAM, MU)

  USE NRTYPE         , ONLY: WP
  USE MOD_INTEGRATION, ONLY: PTS_1, PTS_2, PTS2

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: II, KI, LI
  REAL(WP)             :: L2MZ

  REAL(WP), DIMENSION(:,:,:), INTENT(INOUT) :: MU, LAM

  INTEGER              :: L
  REAL(WP)             :: SUM1

!----------------------------------------------------------------------

  SUM1 = 0.
  DO L = -PTS2+1, PTS2
    SUM1 = SUM1 + 1./L2M_XY()
  END DO

  L2MZ = PTS_1/SUM1

CONTAINS

  FUNCTION L2M_XY()

    REAL(WP)             :: L2M_XY

    INTEGER              :: I, K
    REAL(WP)             :: S1, S2, S3, LAM2MU, LAMBDA

  !----------------------------------------------------------------------

    S1 = 0.
    S2 = 0.
    S3 = 0.
    DO K = -PTS2+1, PTS2
      DO I = -PTS2+1, PTS2
        LAMBDA = LAM(II+I, KI+K, LI+L)
        LAM2MU = LAMBDA + 2. * MU(II+I, KI+K, LI+L)
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


END FUNCTION L2MZ

!======================================================================

FUNCTION LAXZ (II,KI,LI, LAM, MU)

  USE NRTYPE         , ONLY: WP
  USE MOD_INTEGRATION, ONLY: PTS_1, PTS2

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: II, KI, LI
  REAL(WP)             :: LAXZ

  REAL(WP), DIMENSION(:,:,:), INTENT(INOUT) :: MU, LAM

  INTEGER              :: I, L
  REAL(WP)             :: SUM1, SUM2, L2M

!----------------------------------------------------------------------

  SUM1 = 0.
  SUM2 = 0.

  DO L = -PTS2+1, PTS2
    DO I = -PTS2+1, PTS2
      L2M  = L2M_Y()
      SUM1 = SUM1 + (L2M-2.*MU_Y())/L2M
      SUM2 = SUM2 +              1./L2M
    END DO
  END DO

  LAXZ = SUM1/SUM2

CONTAINS

  FUNCTION L2M_Y()

    REAL(WP)             :: L2M_Y
    INTEGER              :: K
    REAL(WP)             :: S1, S2, S3, LAM2MU, LAMBDA

    S1 = 0.
    S2 = 0.
    S3 = 0.
    DO K = -PTS2+1, PTS2
      LAMBDA = LAM(II+I, KI+K, LI+L)
      LAM2MU = LAMBDA + 2. * MU(II+I, KI+K, LI+L)
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

    REAL(WP)             :: MU_Y
    INTEGER              :: K
    REAL(WP)             :: S1

    S1 = 0.
    DO K = -PTS2+1, PTS2
      S1 = S1 + MU(II+I, KI+K, LI+L)
    END DO

    MU_Y = S1/PTS_1

  END FUNCTION MU_Y

END FUNCTION LAXZ

!======================================================================

FUNCTION LAYZ (II,KI,LI, LAM, MU)

  USE NRTYPE         , ONLY: WP
  USE MOD_INTEGRATION, ONLY: PTS_1, PTS2

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: II, KI, LI
  REAL(WP)             :: LAYZ

  REAL(WP), DIMENSION(:,:,:), INTENT(INOUT) :: MU, LAM

  INTEGER              :: K, L
  REAL(WP)             :: SUM1, SUM2, L2M

!----------------------------------------------------------------------

  SUM1 = 0.
  SUM2 = 0.

  DO L = -PTS2+1, PTS2
    DO K = -PTS2+1, PTS2
      L2M  = L2M_X()
      SUM1 = SUM1 + (L2M-2.*MU_X())/L2M
      SUM2 = SUM2 +              1./L2M
    END DO
  END DO

  LAYZ = SUM1/SUM2

CONTAINS

  FUNCTION L2M_X()

    REAL(WP)             :: L2M_X
    INTEGER              :: I
    REAL(WP)             :: S1, S2, S3, LAM2MU, LAMBDA

    S1 = 0.
    S2 = 0.
    S3 = 0.
    DO I = -PTS2+1, PTS2
      LAMBDA = LAM(II+I, KI+K, LI+L)
      LAM2MU = LAMBDA + 2. * MU(II+I, KI+K, LI+L)
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

    REAL(WP)             :: MU_X
    INTEGER              :: I
    REAL(WP)             :: S1

    S1 = 0.
    DO I = -PTS2+1, PTS2
      S1 = S1 + MU(II+I, KI+K, LI+L)
    END DO

    MU_X = S1/PTS_1

  END FUNCTION MU_X

END FUNCTION LAYZ

!======================================================================

FUNCTION LAXY (II,KI,LI, LAM, MU)

  USE NRTYPE         , ONLY: WP
  USE MOD_INTEGRATION, ONLY: PTS_1, PTS2

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: II, KI, LI
  REAL(WP)             :: LAXY

  REAL(WP), DIMENSION(:,:,:), INTENT(INOUT) :: MU, LAM

  INTEGER              :: I, K
  REAL(WP)             :: SUM1, SUM2, L2M

!----------------------------------------------------------------------

  SUM1 = 0.
  SUM2 = 0.

  DO K = -PTS2+1, PTS2
    DO I = -PTS2+1, PTS2
      L2M  = L2M_Z()
      SUM1 = SUM1 + (L2M-2.*MU_Z())/L2M
      SUM2 = SUM2 +              1./L2M
    END DO
  END DO

  LAXY = SUM1/SUM2

CONTAINS

  FUNCTION L2M_Z()

    REAL(WP)             :: L2M_Z
    INTEGER              :: L
    REAL(WP)             :: S1, S2, S3, LAM2MU, LAMBDA

    S1 = 0.
    S2 = 0.
    S3 = 0.
    DO L = -PTS2+1, PTS2
      LAMBDA = LAM(II+I, KI+K, LI+L)
      LAM2MU = LAMBDA + 2. * MU(II+I, KI+K, LI+L)
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

    REAL(WP)             :: MU_Z
    INTEGER              :: L
    REAL(WP)             :: S1

    S1 = 0.
    DO L = -PTS2+1, PTS2
      S1 = S1 + MU(II+I, KI+K, LI+L)
    END DO

    MU_Z = S1/PTS_1

  END FUNCTION MU_Z

END FUNCTION LAXY


END MODULE MOD_EFFMATPAR

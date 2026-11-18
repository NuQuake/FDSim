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
! function returns aproximated value of volume integral of field
! FIELD, where volume is cube with side H/PTS and center with
! coordinates (II,KI,LI)
! note: integral is aproximated by PTS points per each side of cube
!       (that means by PTS^3 points in whole volume)
!-----------------------------------------------------------------------
FUNCTION INT3D (II,KI,LI,FIELD)

  USE NRTYPE         , ONLY: WP
  USE MOD_INTEGRATION, ONLY: PTS, PTS2, PTS_3

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: II, KI, LI
  REAL(WP)             :: INT3D

  REAL(WP), DIMENSION(:,:,:), INTENT(INOUT) :: FIELD

  INTEGER              :: I, K, L
  REAL(WP)             :: SUM

!----------------------------------------------------------------------

  SUM = 0._WP

  DO L = -PTS2+1, PTS2

    DO K = -PTS2+1, PTS2

      DO I = -PTS2+1, PTS2

        SUM = SUM + FIELD(II+I, KI+K, LI+L)

      END DO
    END DO
  END DO

  INT3D = SUM/PTS_3

END FUNCTION INT3D

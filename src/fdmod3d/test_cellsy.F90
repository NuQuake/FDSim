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
! tests 2 FD cells
!-----------------------------------------------------------------------
LOGICAL FUNCTION TEST_CELLSY( LAST_CELL, CELLY )

  USE NRTYPE              , ONLY: SP
  USE MOD_CONST           , ONLY: N_Y_PAR, N_FREQ
  USE MOD_PARAMETERIZATION, ONLY: CELLSY

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER ,                             INTENT(IN) :: LAST_CELL
  REAL(SP), DIMENSION(N_Y_PAR, N_FREQ), INTENT(IN) :: CELLY

  INTEGER :: I

!----------------------------------------------------------------------

  DO I = 1, N_Y_PAR
    IF ( CELLSY(LAST_CELL,I,1)/=CELLY(I,1) ) THEN
      TEST_CELLSY = .FALSE.
      RETURN
    END IF
  END DO

  TEST_CELLSY = .TRUE.

END FUNCTION TEST_CELLSY

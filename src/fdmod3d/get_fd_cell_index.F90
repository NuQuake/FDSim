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
! function returns index if FD cell with the same parameters already
! have one, othewise it return 0 (that means: new FD cell type found)
!-----------------------------------------------------------------------
FUNCTION GET_FD_CELL_INDEX()

  USE MOD_MODEL           , ONLY: KEY_Q
  USE MOD_PARAMETERIZATION, ONLY: CELL, CELLY, LAST_CELL, N_CELLS
  USE MOD_INTERFACES      , ONLY: TEST_CELLS, TEST_CELLSY,             &
                                  ADD_NEW_CELL_TYPE

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: GET_FD_CELL_INDEX
  INTEGER :: J

!----------------------------------------------------------------------

                                        !try last FD cell first
  IF ( N_CELLS > 0 ) THEN
    IF (TEST_CELLS( LAST_CELL, CELL)) THEN
      IF (KEY_Q) THEN
        IF (TEST_CELLSY( LAST_CELL, CELLY)) THEN
          GET_FD_CELL_INDEX = LAST_CELL
          RETURN
        END IF
      ELSE
        GET_FD_CELL_INDEX = LAST_CELL
        RETURN
      END IF
    END IF                               !try all other FD cells

    DO J = 1, N_CELLS
      IF (TEST_CELLS( J, CELL )) THEN
        IF (KEY_Q) THEN

          IF (TEST_CELLSY( J, CELLY)) THEN
            GET_FD_CELL_INDEX = J
            LAST_CELL         = J
            RETURN
          END IF

        ELSE

          GET_FD_CELL_INDEX = J
          LAST_CELL         = J
          RETURN

        END IF
      END IF
    END DO

  END IF

  CALL ADD_NEW_CELL_TYPE

  LAST_CELL         = N_CELLS
  GET_FD_CELL_INDEX = N_CELLS

END FUNCTION GET_FD_CELL_INDEX


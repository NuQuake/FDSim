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
! adds new cell type
!-----------------------------------------------------------------------
SUBROUTINE ADD_NEW_CELL_TYPE

  USE MOD_MODEL           , ONLY: KEY_Q
  USE MOD_PARAMETERIZATION, ONLY: N_CELLS, MAXCELLS,                   &
                                  CELLS, CELLSY, CELL, CELLY

!----------------------------------------------------------------------

  IMPLICIT NONE

!----------------------------------------------------------------------

  N_CELLS = N_CELLS + 1

  IF ( N_CELLS > MAXCELLS ) THEN
    CALL ENLARGE
  END IF

  CELLS(N_CELLS,:) = CELL

  IF(KEY_Q) CELLSY(N_CELLS,:,:) = CELLY

END SUBROUTINE ADD_NEW_CELL_TYPE


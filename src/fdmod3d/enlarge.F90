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
! enlarges array CELLS
!-----------------------------------------------------------------------
SUBROUTINE ENLARGE

  USE NRTYPE              , ONLY: SP
  USE MOD_CONST           , ONLY: N_E_PAR, N_Y_PAR, N_FREQ,            &
                                  NUM_OF_CELLS_TO_ENLARGE
  USE MOD_MODEL           , ONLY: KEY_Q
  USE MOD_PARAMETERIZATION, ONLY: CELLS, CELLSY, MAXCELLS

!----------------------------------------------------------------------

  IMPLICIT NONE

  REAL(SP), DIMENSION(:,:  ), ALLOCATABLE :: TMP
  REAL(SP), DIMENSION(:,:,:), ALLOCATABLE :: TMPY

!----------------------------------------------------------------------

  ALLOCATE ( TMP(MAXCELLS, N_E_PAR) )
  TMP = CELLS
  DEALLOCATE ( CELLS )

  IF ( KEY_Q ) THEN
    ALLOCATE ( TMPY(MAXCELLS, N_Y_PAR, N_FREQ) )
    TMPY = CELLSY
    DEALLOCATE ( CELLSY )
  END IF

  ALLOCATE( CELLS (MAXCELLS + NUM_OF_CELLS_TO_ENLARGE, N_E_PAR) )
  IF (KEY_Q) THEN
    ALLOCATE( CELLSY(MAXCELLS+NUM_OF_CELLS_TO_ENLARGE, N_Y_PAR,N_FREQ) )
  END IF

  CELLS (1:MAXCELLS,: ) = TMP
  DEALLOCATE ( TMP )

  IF ( KEY_Q ) THEN
    CELLSY (1:MAXCELLS,:,:) = TMPY
    DEALLOCATE ( TMPY )
  END IF

  MAXCELLS = MAXCELLS + NUM_OF_CELLS_TO_ENLARGE

END SUBROUTINE ENLARGE


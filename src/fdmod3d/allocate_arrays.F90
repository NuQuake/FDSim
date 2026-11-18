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
!***********************************************************************
! allocates arrays for physical parameters
!***********************************************************************
SUBROUTINE ALLOCATE_ARRAYS (WHICH_GRID)

  USE MOD_CONST           , ONLY: N_FREQ, N_E_PAR, N_Y_PAR
  USE MOD_MODEL           , ONLY: KEY_Q, MXF, MYF, MX, MY, LPAS, MZ
  USE MOD_PARAMETERIZATION, ONLY: FINER_GRID, COARSER_GRID,            &
                                  JMH, JMH3, CELLS, CELLSY, MAXCELLS,  &
                                  N_CELLS

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: WHICH_GRID

!----------------------------------------------------------------------

  SELECT CASE (WHICH_GRID)
    CASE(FINER_GRID)
      ALLOCATE( JMH3  (MYF, MXF, 0:LPAS) )
      ALLOCATE( CELLS (MXF*MYF, N_E_PAR) )
      IF (KEY_Q) THEN
        ALLOCATE( CELLSY(MXF*MYF, N_Y_PAR, N_FREQ) )
      ELSE
        ALLOCATE( CELLSY(      1, N_Y_PAR, N_FREQ) )
      END IF
      N_CELLS = 0
      MAXCELLS = MXF*MYF
    CASE(COARSER_GRID)
      ALLOCATE( JMH   (MY , MX , LPAS:MZ) )
  END SELECT

END SUBROUTINE ALLOCATE_ARRAYS


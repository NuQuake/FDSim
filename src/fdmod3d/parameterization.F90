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
! main procedure of this module
!-----------------------------------------------------------------------
SUBROUTINE PARAMETERIZATION

  USE MOD_MODEL           , ONLY: MZ, LPAS
  USE MOD_PARAMETERIZATION, ONLY: LSTART,                              &
                                  N_PLANES, N_PLANESF,                 &
                                  FINER_GRID, COARSER_GRID
  USE MOD_INTERFACES      , ONLY: ALLOCATE_ARRAYS,                     &
                                  DEALLOCATE_ARRAYS
#if(USE_MPI)
  USE MPI_SHARED
#endif
!----------------------------------------------------------------------

  IMPLICIT NONE

!----------------------------------------------------------------------

  INTERFACE
    SUBROUTINE     PARAMETRIZE_FINER_GRID
    END SUBROUTINE PARAMETRIZE_FINER_GRID

    SUBROUTINE     PARAMETRIZE_COARSER_GRID
    END SUBROUTINE PARAMETRIZE_COARSER_GRID

    SUBROUTINE     PREPARE_INPUT
    END SUBROUTINE PREPARE_INPUT

    SUBROUTINE     PREPARE_INT
    END SUBROUTINE PREPARE_INT
  END INTERFACE

!----------------------------------------------------------------------

#if(USE_MPI)
  IF ( myid == root ) THEN
#endif
    WRITE(*,'(A)')'Parameterize medium using avaraging of physical &
                  &parameters'
#if(USE_MPI)
  END IF
#endif

  CALL PREPARE_INPUT

  CALL PREPARE_INT

  CALL ALLOCATE_ARRAYS(FINER_GRID)

  LSTART = 1

  N_PLANESF = LPAS+1

  CALL PARAMETERIZE_FINER_GRID

  CALL ALLOCATE_ARRAYS(COARSER_GRID)

  LSTART = 1

  N_PLANES = MZ - LPAS + 1

  CALL PARAMETERIZE_COARSER_GRID

END SUBROUTINE PARAMETERIZATION

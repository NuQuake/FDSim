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
! subroutine reads input files with elastic parameters and prepares
! needed variables for anelastic model
!-----------------------------------------------------------------------
SUBROUTINE PREPARE_INPUT

  USE NRTYPE              , ONLY: WP
  USE MOD_FILES
  USE MOD_MODEL           , ONLY: KEY_Q
  USE MOD_PARAMETERIZATION, ONLY: MU_HOM,LA_HOM,RHO_HOM,VS_HOM,VP_HOM

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER  :: LEN, K
  REAL(WP) :: VS,VP

!----------------------------------------------------------------------

  INTERFACE
    SUBROUTINE     PREPARE_Y
    END SUBROUTINE PREPARE_Y
  END INTERFACE

!----------------------------------------------------------------------

  MU_HOM = RHO_HOM * VS_HOM * VS_HOM
  LA_HOM = RHO_HOM * VP_HOM * VP_HOM - 2._WP*MU_HOM

  IF (KEY_Q) THEN
    CALL PREPARE_Y
  END IF

END SUBROUTINE PREPARE_INPUT

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

SUBROUTINE  HPLANE_LPASM1M2_STRESS

  USE CONTROL_DATA, ONLY: LPAS, REFINE

  USE STRESSFIELD , ONLY: TXZ, TYZ, TZZ, TXZF, TYZF, TZZF

  USE INTERFACES  , ONLY: FINE2COARSE

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: L, L_F

!-----------------------------------------------------------------------

  L   = LPAS - 1
  L_F = LPAS - (REFINE+1)/2

  CALL FINE2COARSE ( L = L, L_F = L_F, VAR = TZZ, VARF = TZZF,         &
                     I_SHIFT = .FALSE., K_SHIFT = .FALSE. )


  L   = LPAS - 1
  L_F = LPAS - REFINE

  CALL FINE2COARSE ( L = L, L_F = L_F, VAR = TXZ, VARF = TXZF,         &
                     I_SHIFT = .TRUE. , K_SHIFT = .FALSE. )
  CALL FINE2COARSE ( L = L, L_F = L_F, VAR = TYZ, VARF = TYZF,         &
                     I_SHIFT = .FALSE., K_SHIFT = .TRUE.  )


  L   = LPAS - 2
  L_F = LPAS - (3*REFINE+1)/2

  CALL FINE2COARSE ( L = L, L_F = L_F, VAR = TZZ, VARF = TZZF,         &
                     I_SHIFT = .FALSE., K_SHIFT = .FALSE. )


END SUBROUTINE  HPLANE_LPASM1M2_STRESS

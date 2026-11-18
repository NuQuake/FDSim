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
MODULE MOD_CONST

! module contains all build-time constants
  USE NRTYPE, ONLY: WP

  IMPLICIT NONE

  INTEGER, PARAMETER  :: N_FREQ  = 4 ,                                 &
                         N_E_PAR = 12,                                 &
                         N_Y_PAR = 9 ,                                 &
                         I2N_FREQ1 = 7,                                &
                         NUM_OF_CELLS_TO_ENLARGE = 50000

END MODULE MOD_CONST

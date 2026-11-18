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
MODULE PAR_HA

! This module contains procedures and variables special for
! parameterization of media using volume arithemtic avaraging to
! calculate grid densities and volume harmonic avarage to calculate
! grid parameters bulk & shear moduli and dampnig coeficients Y

! This module uses input coordinate system, but positions in FD cell
! are considered as in output coordinate system.

USE MOD_PARAMETERIZATION
USE MOD_CONST
USE MOD_MODEL
USE MOD_FILES
!USE MOD_STOP

!------------------------------------------------------------- variables
IMPLICIT NONE

INTEGER,  SAVE :: WHICH_GRID
LOGICAL,  SAVE :: RESTART

END MODULE PAR_HA

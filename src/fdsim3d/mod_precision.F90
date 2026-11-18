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

MODULE PRECISION
  IMPLICIT NONE
  INTEGER, PARAMETER :: SP = SELECTED_REAL_KIND (6,30)
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND (13,60)
#if(USE_MPI)
  INTEGER, PARAMETER :: IK = SELECTED_INT_KIND  (8)  !!! DO NOT CHANGE !!!
  INTEGER, PARAMETER :: PP = SP                      !!! DO NOT CHANGE !!!
#else
  INTEGER, PARAMETER :: IK = SELECTED_INT_KIND  (8)
  INTEGER, PARAMETER :: PP = SP
#endif

END MODULE PRECISION

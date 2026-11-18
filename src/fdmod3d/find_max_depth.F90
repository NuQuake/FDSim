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
! function returns the maximum depth of sediments/bedrock interface
! (i.e. interface NI)
!-----------------------------------------------------------------------
SUBROUTINE FIND_MAX_DEPTH

  USE NRTYPE   , ONLY: WP
  USE MOD_MODEL, ONLY: MAXDEPTH, ITF1, NI

!----------------------------------------------------------------------

  IMPLICIT NONE

!----------------------------------------------------------------------

  MAXDEPTH = REAL(MINVAL(ITF1(:,:,NI)),WP)

END SUBROUTINE FIND_MAX_DEPTH

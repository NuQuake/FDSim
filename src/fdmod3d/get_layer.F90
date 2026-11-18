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
! function returns layer, into which point x,y,z belongs
!-----------------------------------------------------------------------
FUNCTION GET_LAYER (X,Y,Z)

  USE NRTYPE   , ONLY: WP
  USE MOD_MODEL, ONLY: LAYER_1, LAYER_2, NI
  USE MOD_INTERFACES, ONLY: GETDEPTH

!----------------------------------------------------------------------

  IMPLICIT NONE

  REAL(WP), INTENT(IN) :: X, Y, Z
  INTEGER              :: GET_LAYER

  REAL(WP) :: ZS

!----------------------------------------------------------------------

  ZS = GETDEPTH(X,Y,NI)

  IF ( Z > ZS ) THEN
                          ! in sediments (1st layer)
    GET_LAYER = LAYER_1
    RETURN
  ELSE
                          ! in bedrock   (2nd layer)
    GET_LAYER = LAYER_2
    RETURN
  END IF

END FUNCTION GET_LAYER

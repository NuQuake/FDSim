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
! function returns value of MU, LAMBDA and RHO as function
! of coordinates x,y,z
!-----------------------------------------------------------------------
FUNCTION FUNC_MU_LA_RHO(X,Y,Z)

  USE NRTYPE              , ONLY: WP
  USE MOD_MODEL           , ONLY: LAYER_1, LAYER_2
  USE MOD_PARAMETERIZATION, ONLY: MU_HOM, LA_HOM, RHO_HOM
  USE MOD_INTERFACES      , ONLY: GET_LAYER
  USE MOD_FUNC            , ONLY: FUNCVS   , FUNCVP   , FUNCRHO    ,   &
                                  FUNCVS_1D, FUNCVP_1D, FUNCRHO_1D

!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL(WP), INTENT(IN) :: X,Y,Z
  REAL(WP)             :: FUNC_MU_LA_RHO(3)

  REAL(WP) :: VS, VP, RHO, MU
  INTEGER  :: LAYER

!-----------------------------------------------------------------------

  LAYER = GET_LAYER(X,Y,Z)

  SELECT CASE (LAYER)
  CASE(LAYER_1)  ! in sediments (1st layer)
    VP     = FUNCVP ( X,Y,Z )
    VS     = FUNCVS ( X,Y,Z )
    RHO    = FUNCRHO( X,Y,Z )
    MU     = RHO * VS * VS
    FUNC_MU_LA_RHO (1) = MU
    FUNC_MU_LA_RHO (2) = RHO*VP*VP - 2.*MU
    FUNC_MU_LA_RHO (3) = RHO
    RETURN
  CASE(LAYER_2)  ! in bedrock   (2nd layer)
    VP     = FUNCVP_1D ( Z )
    VS     = FUNCVS_1D ( Z )
    RHO    = FUNCRHO_1D( Z )
    MU     = RHO * VS * VS
    FUNC_MU_LA_RHO (1) = MU
    FUNC_MU_LA_RHO (2) = RHO*VP*VP - 2.*MU
    FUNC_MU_LA_RHO (3) = RHO
    RETURN
  END SELECT


END FUNCTION FUNC_MU_LA_RHO

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
! function returns true if whole volume over which avarage of any
! parameter should be calculated is inside 1D-inhomogeneous medium
! (lies inside 1D-IH bedrock)
! (X,Y,Z) - position of "free point" in FD cell
!***********************************************************************
FUNCTION IS_IN_1D_MEDIUM(X,Y,Z,WHICH_GRID)

  USE NRTYPE        , ONLY: WP
  USE MOD_MODEL     , ONLY: H, HF, H2, HF2, ni
  USE MOD_PARAMETERIZATION, ONLY: FINER_GRID, COARSER_GRID
  USE MOD_INTERFACES, ONLY: GETDEPTH

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,  INTENT(IN) :: WHICH_GRID
  REAL(WP), INTENT(IN) :: X,Y,Z
  LOGICAL              :: IS_IN_1D_MEDIUM

  REAL(WP) :: X1,X2,Y1,Y2,GH,GH2,Z0

!----------------------------------------------------------------------

  SELECT CASE(WHICH_GRID)
  CASE(FINER_GRID)
    X1 = X - HF2
    X2 = X + HF
    Y1 = Y - HF2
    Y2 = Y + HF
    GH  = HF
    GH2 = HF2
  CASE(COARSER_GRID)
    X1 = X - H2
    X2 = X + H
    Y1 = Y - H2
    Y2 = Y + H
    GH  = H
    GH2 = H2
  END SELECT

  Z0 = GETDEPTH(X1,Y1,ni) - GH2
  IF( Z0 > Z ) THEN
    Z0 = GETDEPTH(X2,Y1,ni) - GH2
    IF( Z0 > Z ) THEN
      Z0 = GETDEPTH(X1,Y2,ni) - GH2
      IF( Z0 > Z ) THEN
        Z0 = GETDEPTH(X2,Y2,ni) - GH2
        IF( Z0 > Z ) THEN
          IS_IN_1D_MEDIUM = .TRUE.
        ELSE
          IS_IN_1D_MEDIUM = .FALSE.
        END IF
      ELSE
        IS_IN_1D_MEDIUM = .FALSE.
      END IF
    ELSE
      IS_IN_1D_MEDIUM = .FALSE.
    END IF
  ELSE
    IS_IN_1D_MEDIUM = .FALSE.
  END IF

END FUNCTION IS_IN_1D_MEDIUM

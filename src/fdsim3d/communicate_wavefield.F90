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

SUBROUTINE COMMUNICATE_WAVEFIELD

  USE CONTROL_DATA, ONLY: LPAS
  USE WAVEFIELD   , ONLY: UM, VM, WM, UMF, VMF, WMF

!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------

  INTERFACE
    SUBROUTINE     COMMUNICATE (VARIABLE)
    USE PRECISION   , ONLY : PP
    USE AUXIL       , ONLY : X_MIN_D, X_MAX_D,                         &
                             Y_MIN_D, Y_MAX_D,                         &
                             Z_MIN_D, Z_MAX_D
    IMPLICIT NONE
    REAL(PP), INTENT(INOUT),                                           &
      DIMENSION(X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) ::  &
                                                              VARIABLE
    END SUBROUTINE COMMUNICATE

    SUBROUTINE     COMMUNICATE_F (VARIABLE)
    USE PRECISION   , ONLY : PP
    USE AUXIL       , ONLY : X_MIN_DF, X_MAX_DF,                       &
                             Y_MIN_DF, Y_MAX_DF,                       &
                             Z_MIN_DF, Z_MAX_DF
    IMPLICIT NONE
    REAL(PP), INTENT(INOUT), DIMENSION                                 &
         (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) ::  &
                                                              VARIABLE
    END SUBROUTINE COMMUNICATE_F

  END INTERFACE

!-----------------------------------------------------------------------

  CALL COMMUNICATE   ( UM  )
  CALL COMMUNICATE   ( VM  )
  CALL COMMUNICATE   ( WM  )

  IF ( LPAS /= 0 ) THEN
    CALL COMMUNICATE_F ( UMF )
    CALL COMMUNICATE_F ( VMF )
    CALL COMMUNICATE_F ( WMF )
  END IF


END SUBROUTINE COMMUNICATE_WAVEFIELD

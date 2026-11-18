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

SUBROUTINE COMMUNICATE_STRESSFIELD

  USE CONTROL_DATA, ONLY: LPAS, PWINC
  USE WAVEFIELD   , ONLY: XXX,   XYY,   XZZ,   XXY,   XXZ,   XYZ ,     &
                          XXXF,  XYYF,  XZZF,  XXYF,  XXZF,  XYZF
  USE STRESSFIELD , ONLY: TXX,   TYY,   TZZ,   TXY,   TXZ,   TYZ ,     &
                          TXXF,  TYYF,  TZZF,  TXYF,  TXZF,  TYZF
  USE AUXIL       , ONLY: TXZ_RES, TYZ_RES, TXX_TOT, TYY_TOT, TZZ_TOT

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

    SUBROUTINE     COMMUNICATE_PLANE (VARIABLE)
    USE PRECISION   , ONLY : PP
    USE AUXIL       , ONLY : X_MIN_D, X_MAX_D,                         &
                             Y_MIN_D, Y_MAX_D
    IMPLICIT NONE
    REAL(PP), INTENT(INOUT),                                           &
               DIMENSION(X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D) :: VARIABLE
    END SUBROUTINE COMMUNICATE_PLANE

    SUBROUTINE     COMMUNICATE_F (VARIABLE)
    USE PRECISION   , ONLY : PP
    USE AUXIL       , ONLY : X_MIN_DF, X_MAX_DF,                       &
                             Y_MIN_DF, Y_MAX_DF,                       &
                             Z_MIN_DF, Z_MAX_DF
    IMPLICIT NONE
    REAL(PP), INTENT(INOUT), DIMENSION                                 &
            (X_MIN_DF:X_MAX_DF,Y_MIN_DF:Y_MAX_DF,Z_MIN_DF:Z_MAX_DF)::  &
                                                              VARIABLE
    END SUBROUTINE COMMUNICATE_F

  END INTERFACE

!-----------------------------------------------------------------------

  CALL COMMUNICATE   ( TXX  )
  CALL COMMUNICATE   ( TYY  )
  CALL COMMUNICATE   ( TZZ  )
  CALL COMMUNICATE   ( TXY  )
  CALL COMMUNICATE   ( TXZ  )
  CALL COMMUNICATE   ( TYZ  )
  CALL COMMUNICATE   ( XXX  )
  CALL COMMUNICATE   ( XYY  )
  CALL COMMUNICATE   ( XZZ  )
  CALL COMMUNICATE   ( XXY  )
  CALL COMMUNICATE   ( XXZ  )
  CALL COMMUNICATE   ( XYZ  )

  IF ( PWINC == 1 ) THEN
    CALL COMMUNICATE_PLANE   ( TXZ_RES  )
    CALL COMMUNICATE_PLANE   ( TYZ_RES  )
  END IF

  IF ( PWINC == 2 ) THEN
    CALL COMMUNICATE_PLANE   ( TXX_TOT  )
    CALL COMMUNICATE_PLANE   ( TYY_TOT  )
    CALL COMMUNICATE_PLANE   ( TZZ_TOT  )
  END IF

  IF ( LPAS /= 0 ) THEN
    CALL COMMUNICATE_F ( TXXF )
    CALL COMMUNICATE_F ( TYYF )
    CALL COMMUNICATE_F ( TZZF )
    CALL COMMUNICATE_F ( TXYF )
    CALL COMMUNICATE_F ( TXZF )
    CALL COMMUNICATE_F ( TYZF )
    CALL COMMUNICATE_F ( XXXF )
    CALL COMMUNICATE_F ( XYYF )
    CALL COMMUNICATE_F ( XZZF )
    CALL COMMUNICATE_F ( XXYF )
    CALL COMMUNICATE_F ( XXZF )
    CALL COMMUNICATE_F ( XYZF )
  END IF


END SUBROUTINE COMMUNICATE_STRESSFIELD

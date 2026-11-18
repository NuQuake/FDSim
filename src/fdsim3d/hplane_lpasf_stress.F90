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

SUBROUTINE  HPLANE_LPASF_STRESS

  USE CONTROL_DATA, ONLY: LPAS

  USE AUXIL       , ONLY:  X_MIN_D , X_MAX_D , X_MIN , X_MAX ,         &
                           Y_MIN_D , Y_MAX_D , Y_MIN , Y_MAX

  USE STRESSFIELD , ONLY: TXZ, TYZ

  USE INTERFACES  , ONLY: COARSE2FINE_TXZ, COARSE2FINE_TYZ

!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------

  IF ( X_MIN_D == X_MIN ) THEN
    TXZ (X_MIN,:,LPAS) = 0.
    TYZ (X_MIN,:,LPAS) = 0.
!    TXZ (X_MIN,:,LPAS) = 2.*TXZ(X_MIN+1,:,LPAS) - TXZ (X_MIN+2,:,LPAS)
!    TYZ (X_MIN,:,LPAS) = 2.*TYZ(X_MIN+1,:,LPAS) - TYZ (X_MIN+2,:,LPAS)
  END IF

  IF ( Y_MIN_D == Y_MIN ) THEN
    TXZ (:,Y_MIN,LPAS) = 0.
    TYZ (:,Y_MIN,LPAS) = 0.
!    TXZ (:,Y_MIN,LPAS) = 2.*TXZ(:,Y_MIN+1,LPAS) - TXZ (:,Y_MIN+2,LPAS)
!    TYZ (:,Y_MIN,LPAS) = 2.*TYZ(:,Y_MIN+1,LPAS) - TYZ (:,Y_MIN+2,LPAS)
  END IF

  IF ( X_MAX_D == X_MAX ) THEN
!    TYZ (X_MAX,:,LPAS) = 2.*TYZ(X_MAX-1,:,LPAS) - TYZ (X_MAX-2,:,LPAS)
    TYZ (X_MAX,:,LPAS) = 0.
  END IF

  IF ( Y_MAX_D == Y_MAX ) THEN
!    TXZ (:,Y_MAX,LPAS) = 2.*TXZ(:,Y_MAX-1,LPAS) - TXZ (:,Y_MAX-2,LPAS)
    TXZ (:,Y_MAX,LPAS) = 0.
  END IF

  CALL COARSE2FINE_TXZ
  CALL COARSE2FINE_TYZ


END SUBROUTINE  HPLANE_LPASF_STRESS

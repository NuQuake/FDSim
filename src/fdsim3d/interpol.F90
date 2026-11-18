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

FUNCTION INTERPOL (X,Y,Z, A, ORDER)

  USE PRECISION   , ONLY: PP

  IMPLICIT NONE

  REAL(PP),                                      INTENT(IN) :: X,Y,Z
  REAL(PP), DIMENSION (1:2,1:2,1:2), INTENT(IN) :: A
  REAL(PP)                                                  :: INTERPOL
  INTEGER ,                                      INTENT(IN) :: ORDER

  REAL(PP), DIMENSION (1:2)                             :: XX, YY, ZZ
  REAL(PP)                                                  :: RES
  INTEGER                                                   :: I, K, L

  IF ( ORDER == 2 ) THEN

    XX(1) = 1. - X
    XX(2) =      X
    YY(1) = 1. - Y
    YY(2) =      Y
    ZZ(1) = 1. - Z
    ZZ(2) =      Z

  ELSE IF ( ORDER == 4 ) THEN  ! STORE_SEISMOGRAMS HAS TO BE CHANGED !!!

!    XX(1) = 1./6.*(    -2*X + 3*X*X - X*X*X)
!    XX(2) = 1./2.*( 2. -  X - 2*X*X + X*X*X)
!    XX(3) = 1./2.*(     2*X +   X*X - X*X*X)
!    XX(4) = 1./6.*(    -  X +         X*X*X)
!
!    YY(1) = 1./6.*(    -2*Y + 3*Y*Y - Y*Y*Y)
!    YY(2) = 1./2.*( 2. -  Y - 2*Y*Y + Y*Y*Y)
!    YY(3) = 1./2.*(     2*Y +   Y*Y - Y*Y*Y)
!    YY(4) = 1./6.*(    -  Y +         Y*Y*Y)
!
!    ZZ(1) = 1./6.*(    -2*Z + 3*Z*Z - Z*Z*Z)
!    ZZ(2) = 1./2.*( 2. -  Z - 2*Z*Z + Z*Z*Z)
!    ZZ(3) = 1./2.*(     2*Z +   Z*Z - Z*Z*Z)
!    ZZ(4) = 1./6.*(    -  Z +         Z*Z*Z)

  END IF

  RES = 0.
  DO I = 1, ORDER
    DO K = 1, ORDER
      DO L = 1, ORDER
        RES = RES + A(I,K,L)*XX(I)*YY(K)*ZZ(L)
      END DO
    END DO
  END DO

  INTERPOL = RES

  !INTERPOL = A(1,1,1) * (1.-X) * (1.-Y) * (1.-Z)                       &
  !         + A(2,1,1) *     X  * (1.-Y) * (1.-Z)                       &
  !         + A(1,2,1) * (1.-X) *     Y  * (1.-Z)                       &
  !         + A(2,2,1) *     X  *     Y  * (1.-Z)                       &
  !         + A(1,1,2) * (1.-X) * (1.-Y) *     Z                        &
  !         + A(2,1,2) *     X  * (1.-Y) *     Z                        &
  !         + A(1,2,2) * (1.-X) *     Y  *     Z                        &
  !         + A(2,2,2) *     X  *     Y  *     Z

END FUNCTION

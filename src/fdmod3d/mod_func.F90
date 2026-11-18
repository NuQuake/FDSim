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
MODULE MOD_FUNC

! module contains functions for calculation of
! P and S-wave velocities, density and quality factors
! in 3D sediments (1st layer) and in 1D bedrock (2nd layer)

! Note: Z IS LESS OR EQUAL TO 0

use mod_model,      only: ni
use mod_interfaces, only: getdepth

IMPLICIT NONE

CONTAINS
!==================== 3D INHOMOGENEOUS STRUCTURE =======================
!------------------------------------------------------- S-wave velocity
  FUNCTION FUNCVS(X,Y,Z)

    USE NRTYPE, ONLY: WP

    REAL(WP), INTENT(IN) :: X,Y,Z
    REAL(WP)             :: FUNCVS
    REAL(WP)             :: H1, H2, H3, Z1, Z2, Z3

    H1 = - 17.3
    H2 = - 72.5
    H3 = -115.6

    IF ( ( Y < 2500. ) .AND. ( Y > -2500. ) ) THEN
      Z1 =      H1 * MIN ( 1., (2500.-Y) / 1500. )
      Z2 = Z1 + H2 * MIN ( 1., (2500.-Y) / 1500. )
      Z3 = Z2 + H3 * MIN ( 1., (2500.-Y) / 1500. )
    ELSE
      Z1 = 0.
      Z2 = 0.
      Z3 = 0.
    END IF

    IF ( Z > Z1 ) THEN
      FUNCVS = 200.
    ELSE IF ( Z > Z2 ) THEN
      FUNCVS = 350.
    ELSE IF ( Z > Z3 ) THEN
      FUNCVS = 650.
    ELSE
      FUNCVS = 2600.
    END IF

  END FUNCTION FUNCVS

!------------------------------------------------------- P-wave velocity
  FUNCTION FUNCVP(X,Y,Z)

    USE NRTYPE, ONLY: WP

    REAL(WP), INTENT(IN) :: X,Y,Z
    REAL(WP)             :: FUNCVP
    REAL(WP)             :: H1, H2, H3, Z1, Z2, Z3

    H1 = - 17.3
    H2 = - 72.5
    H3 = -115.6

    IF ( ( Y < 2500. ) .AND. ( Y > -2500. ) ) THEN
      Z1 =      H1 * MIN ( 1., (2500.-Y) / 1500. )
      Z2 = Z1 + H2 * MIN ( 1., (2500.-Y) / 1500. )
      Z3 = Z2 + H3 * MIN ( 1., (2500.-Y) / 1500. )
    ELSE
      Z1 = 0.
      Z2 = 0.
      Z3 = 0.
    END IF

    IF ( Z > Z1 ) THEN
      FUNCVP = 1500.
    ELSE IF ( Z > Z2 ) THEN
      FUNCVP = 1800.
    ELSE IF ( Z > Z3 ) THEN
      FUNCVP = 2500.
    ELSE
      FUNCVP = 4500.
    END IF

  END FUNCTION FUNCVP

!--------------------------------------------------------------- Density
  FUNCTION FUNCRHO(X,Y,Z)

    USE NRTYPE, ONLY: WP

    REAL(WP), INTENT(IN) :: X,Y,Z
    REAL(WP)             :: FUNCRHO
    REAL(WP)             :: H1, H2, H3, Z1, Z2, Z3

    H1 = - 17.3
    H2 = - 72.5
    H3 = -115.6

    IF ( ( Y < 2500. ) .AND. ( Y > -2500. ) ) THEN
      Z1 =      H1 * MIN ( 1., (2500.-Y) / 1500. )
      Z2 = Z1 + H2 * MIN ( 1., (2500.-Y) / 1500. )
      Z3 = Z2 + H3 * MIN ( 1., (2500.-Y) / 1500. )
    ELSE
      Z1 = 0.
      Z2 = 0.
      Z3 = 0.
    END IF

    IF ( Z > Z1 ) THEN
      FUNCRHO = 2100.
    ELSE IF ( Z > Z2 ) THEN
      FUNCRHO = 2100.
    ELSE IF ( Z > Z3 ) THEN
      FUNCRHO = 2200.
    ELSE
      FUNCRHO = 2600.
    END IF

  END FUNCTION FUNCRHO

!------------------------------------------------- S-wave quality factor
  FUNCTION FUNCQS(X,Y,Z)

    USE NRTYPE, ONLY: WP

    REAL(WP), INTENT(IN) :: X,Y,Z
    REAL(WP)             :: FUNCQS
    REAL(WP)             :: Z1, Z2, Z3

    FUNCQS = 9999.

  END FUNCTION FUNCQS

!------------------------------------------------- P-wave quality factor
  FUNCTION FUNCQP(X,Y,Z)

    USE NRTYPE, ONLY: WP

    REAL(WP), INTENT(IN) :: X,Y,Z
    REAL(WP)             :: FUNCQP

    FUNCQP = 9999.

  END FUNCTION FUNCQP

!==================== 1D INHOMOGENEOUS BEDROCK =========================
!------------------------------------------------------- S-wave velocity
  FUNCTION FUNCVS_1D(Z)

    USE NRTYPE, ONLY: WP

    REAL(WP), INTENT(IN)   :: Z
    REAL(WP)               :: FUNCVS_1D

    FUNCVS_1D = 2600.

  END FUNCTION FUNCVS_1D

!------------------------------------------------------- P-wave velocity
  FUNCTION FUNCVP_1D(Z)

    USE NRTYPE, ONLY: WP

    REAL(WP), INTENT(IN) :: Z
    REAL(WP)             :: FUNCVP_1D

    FUNCVP_1D = 4500.

  END FUNCTION FUNCVP_1D


!--------------------------------------------------------------- Density
  FUNCTION FUNCRHO_1D(Z)

    USE NRTYPE, ONLY: WP

    REAL(WP), INTENT(IN) :: Z
    REAL(WP)             :: FUNCRHO_1D

    FUNCRHO_1D = 2600.

  END FUNCTION FUNCRHO_1D

!------------------------------------------------- S-wave quality factor
  FUNCTION FUNCQS_1D(Z)

    USE NRTYPE, ONLY: WP

    REAL(WP), INTENT(IN) :: Z
    REAL(WP)             :: FUNCQS_1D

    FUNCQS_1D = 9999.

  END FUNCTION FUNCQS_1D

!------------------------------------------------- P-wave quality factor
  FUNCTION FUNCQP_1D(Z)

    USE NRTYPE, ONLY: WP

    REAL(WP), INTENT(IN) :: Z
    REAL(WP)             :: FUNCQP_1D

    FUNCQP_1D = 9999.

  END FUNCTION

END MODULE

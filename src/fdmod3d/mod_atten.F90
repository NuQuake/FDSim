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
MODULE MOD_ATTEN

  USE NRTYPE    , ONLY: WP
  USE MOD_CONST , ONLY: N_FREQ, I2N_FREQ1

  IMPLICIT NONE

  REAL   (WP), DIMENSION(N_FREQ,I2N_FREQ1) :: FJJ, FJK
  COMPLEX(WP), DIMENSION(N_FREQ,I2N_FREQ1) :: CF
  REAL   (WP), DIMENSION(  N_FREQ)         :: Y, B, P, FRJ, PFP21, PF21
  REAL   (WP), DIMENSION(  N_FREQ, N_FREQ) :: A
  REAL   (WP), DIMENSION(I2N_FREQ1)        :: Q, PLO, FRK

  INTERFACE

    SUBROUTINE     QYLSM(QQ,YY)
      USE NRTYPE   , ONLY: WP
      USE MOD_CONST, ONLY: N_FREQ, I2N_FREQ1
      IMPLICIT NONE
      REAL(WP), DIMENSION (  N_FREQ ), INTENT(INOUT) :: YY
      REAL(WP), DIMENSION (I2N_FREQ1), INTENT(INOUT) :: QQ
    END SUBROUTINE QYLSM

  END INTERFACE

END MODULE MOD_ATTEN

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

SUBROUTINE HPLANE_PMLZH0_VELOCITY

  USE PRECISION   , ONLY:                                              &
                          PP
  USE CONTROL_DATA, ONLY:                                              &
                          MX, MY, MZ, TPML
  USE WAVEFIELD   , ONLY:                                              &
                          UM, VM, WM
  USE AUXIL       , ONLY:                                              &
                           X_MIN  , X_MAX  ,                           &
                           Y_MIN  , Y_MAX
  USE NONREF_BOUND, ONLY:                                              &
                          NR_DATA_Z,                                   &
                          SZ_U, SZ_V, SZ_W, SZ_UP, SZ_VP, SZ_WP,       &
                          ABW01,   ABW02,   ABW10,   ABW11,            &
                          ABW12,   ABW20,   ABW21,   ABW22,            &
                          ABH01,   ABH02,   ABH10,   ABH11,            &
                          ABH12,   ABH20,   ABH21,   ABH22

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER             :: I, K, LST, LST1, LST2

!-----------------------------------------------------------------------

  LST  = MZ + TPML
  LST1 = LST-1
  LST2 = LST-2

  DO  K = MAX(1-TPML,Y_MIN), MIN(MY+TPML,Y_MAX)
    DO  I = MAX(1-TPML,X_MIN), MIN(MX+TPML,X_MAX)

      UM(I,K,LST) =                  ABH01 *    UM(I,K,   LST1)        &
                                   + ABH02 *    UM(I,K,   LST2)        &
       + ABH10 * SZ_U (I,K)%TIKMZ0 + ABH11 * SZ_U (I,K)%TIKMZ1         &
                                   + ABH12 * SZ_U (I,K)%TIKMZ2         &
       + ABH20 * SZ_UP(I,K)%TIKMZ0 + ABH21 * SZ_UP(I,K)%TIKMZ1         &
                                   + ABH22 * SZ_UP(I,K)%TIKMZ2

      VM(I,K,LST) =                  ABH01 *    VM(I,K,   LST1)        &
                                   + ABH02 *    VM(I,K,   LST2)        &
       + ABH10 * SZ_V (I,K)%TIKMZ0 + ABH11 * SZ_V (I,K)%TIKMZ1         &
                                   + ABH12 * SZ_V (I,K)%TIKMZ2         &
       + ABH20 * SZ_VP(I,K)%TIKMZ0 + ABH21 * SZ_VP(I,K)%TIKMZ1         &
                                   + ABH22 * SZ_VP(I,K)%TIKMZ2

      WM(I,K,LST) =                  ABW01 *    WM(I,K,   LST1)        &
                                   + ABW02 *    WM(I,K,   LST2)        &
       + ABW10 * SZ_W (I,K)%TIKMZ0 + ABW11 * SZ_W (I,K)%TIKMZ1         &
                                   + ABW12 * SZ_W (I,K)%TIKMZ2         &
       + ABW20 * SZ_WP(I,K)%TIKMZ0 + ABW21 * SZ_WP(I,K)%TIKMZ1         &
                                   + ABW22 * SZ_WP(I,K)%TIKMZ2

    END DO
  END DO


END SUBROUTINE HPLANE_PMLZH0_VELOCITY

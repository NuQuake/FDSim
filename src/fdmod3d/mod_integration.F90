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
MODULE MOD_INTEGRATION

  USE NRTYPE   , ONLY: WP

  IMPLICIT NONE

  INTEGER  :: PTS, PTS2
  REAL(WP) :: PTS_1, PTS_2, PTS_3

  REAL(WP)   , DIMENSION (:,:,:)  , ALLOCATABLE ::                     &
               MU_FINE, LA_FINE, RHO_FINE

  COMPLEX(WP), DIMENSION (:,:,:,:), ALLOCATABLE ::                     &
              YMU_FINE, YLA_FINE

END MODULE MOD_INTEGRATION

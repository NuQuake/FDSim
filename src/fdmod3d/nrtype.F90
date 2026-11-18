!=======================================================================
! Copyright (c) 2014-2026 Jozef Kristek, Peter Moczo, Martin Galis
!
! Licensed under the EUPL (the "Licence")
! You may not use this work except in compliance with the Licence.
! You may obtain a copy of the Licence at
! https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12
! Unless required by applicable law or agreed to in writing, software
! distributed under the Licence is distributed on an "AS IS" basis,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and
! limitations under the Licence.
!
! REVISION HISTORY:
! Version 1.0, 2014-03-19: Jozef Kristek, Peter Moczo, Martin Galis
! 2026-11-18: rewritten; contains only the kinds used by FDMod3D
!=======================================================================
!
! Real kind parameters used throughout FDMod3D.
!   SP - single precision
!   DP - double precision
!   WP - working precision of the model preparation
!-----------------------------------------------------------------------
MODULE NRTYPE

  IMPLICIT NONE

  INTEGER, PARAMETER :: SP = KIND(1.0)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER, PARAMETER :: WP = DP

END MODULE NRTYPE

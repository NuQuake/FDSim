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

MODULE STRESSFIELD
  USE PRECISION
  IMPLICIT NONE

  REAL    (KIND=PP), DIMENSION (:,:,:), ALLOCATABLE ::                 &
                                        TXX , TYY , TXY ,              &
                                        TXZ , TYZ , TZZ ,              &
                                        TXXF, TYYF, TXYF,              &
                                        TXZF, TYZF, TZZF
!$SGI DISTRIBUTE TXZ  (*,BLOCK,*,*),TYZ  (*,BLOCK,*,*)
!$SGI DISTRIBUTE TXX  (*,BLOCK,*,*),TYY  (*,BLOCK,*,*)
!$SGI DISTRIBUTE TZZ  (*,BLOCK,*,*),TXY  (*,BLOCK,*,*)

END MODULE STRESSFIELD

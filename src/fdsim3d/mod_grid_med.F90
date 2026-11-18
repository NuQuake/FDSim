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

MODULE GRID_MEDIUM
  USE PRECISION, ONLY: PP, IK
  IMPLICIT NONE

  CHARACTER(LEN=40) :: TEXT

  INTEGER           :: JMNUM

  INTEGER(IK), DIMENSION (:,:,:), ALLOCATABLE :: JM, JMF
!$SGI DISTRIBUTE JM(*,BLOCK,*)

  REAL (KIND=PP) :: H

  REAL (KIND=PP), DIMENSION (:), ALLOCATABLE ::                        &
                           DENU , DENV , DENW ,                        &
                           DENUF, DENVF, DENWF,                        &
                           L2MX , L2MY , L2MZ , LAXY , LAXZ , LAYZ,    &
                                                MUXY , MUXZ , MUYZ,    &
                          CL2MX ,CL2MY ,CL2MZ ,CLAXY ,CLAXZ ,CLAYZ,    &
                          CMUXY ,CMUXZ ,CMUYZ
  REAL (KIND=PP), DIMENSION (:,:), ALLOCATABLE ::                      &
                          YL2MX , YL2MY , YL2MZ ,                      &
                          YLAXY , YLAXZ , YLAYZ ,                      &
                          YMUXY , YMUXZ , YMUYZ

END MODULE GRID_MEDIUM

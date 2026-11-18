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

SUBROUTINE  COMPLIANCES

  USE PRECISION   , ONLY:                                              &
                          PP
  USE GRID_MEDIUM , ONLY:                                              &
                          JMNUM,                                       &
                          L2MX , L2MY , L2MZ , LAXY , LAXZ , LAYZ ,    &
                          MUXY , MUXZ , MUYZ ,                         &
                         CL2MX ,CL2MY ,CL2MZ ,CLAXY ,CLAXZ ,CLAYZ ,    &
                         CMUXY ,CMUXZ ,CMUYZ

!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL(PP),DIMENSION(JMNUM) :: DENOM

!-----------------------------------------------------------------------

  CMUXY = 1./MUXY
  CMUXZ = 1./MUXZ
  CMUYZ = 1./MUYZ

  DENOM = L2MX*L2MY*L2MZ + 2*LAXY*LAXZ*LAYZ                            &
        - L2MY*LAXZ**2 - L2MX*LAYZ**2 - L2MZ*LAXY**2

  CL2MX = (L2MY*L2MZ - LAYZ**2) / DENOM
  CL2MY = (L2MX*L2MZ - LAXZ**2) / DENOM
  CL2MZ = (L2MX*L2MY - LAXY**2) / DENOM

  CLAXY = ( -L2MZ*LAXY + LAXZ*LAYZ ) / DENOM
  CLAXZ = ( -L2MY*LAXZ + LAXY*LAYZ ) / DENOM
  CLAYZ = ( -L2MX*LAYZ + LAXY*LAXZ ) / DENOM


END SUBROUTINE  COMPLIANCES

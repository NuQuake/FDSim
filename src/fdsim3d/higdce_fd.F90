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

SUBROUTINE  HIGCE_FD  (WB, ALFA, BETA, DT, H,                          &
                       A00, A01, A02, A10, A11, A12, A20, A21, A22)
  USE PRECISION

  IMPLICIT NONE

  REAL (PP), INTENT (IN)  :: WB, ALFA, BETA, DT, H

  REAL (DP), INTENT (OUT) :: A00, A01, A02, A10, A11, A12, A20, A21, A22

  REAL (DP)               :: H1X, H1T, HXT, C1T, CXT, GP, GS, GS1, HC


  GP = ALFA*DT/H
  GS = BETA*DT/H

  GS1 = 1._DP + GS

  H1X = ( GS    - WB*GS1 ) / ( GS1*(1-WB) )
  H1T = ( 1._DP - WB*GS1 ) / ( GS1*(1-WB) )
  HXT = WB / ( 1._DP - WB )

  C1T = 1._DP - GP
  CXT = GP

  HC = H1T * C1T


  A00 =     0._DP

  A01 =     H1X
  A02 =     0._DP
  A10 =         C1T + H1T
  A11 =         CXT + HXT    - H1X*C1T
  A12 =   - H1X*CXT
  A20 =   - H1T*C1T
  A21 =   - HXT*C1T - H1T*CXT
  A22 =   - HXT*CXT

END SUBROUTINE  HIGCE_FD

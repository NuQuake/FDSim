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

SUBROUTINE  HIGCE2_FD  (WB, VELH, VELCE, DT, H,                        &
                        A00, A01, A02, A10, A11, A12, A20, A21, A22)
  USE PRECISION

  IMPLICIT NONE

  REAL (KIND=PP), INTENT (IN)  ::                                      &
                    WB, VELH, VELCE, DT, H

  REAL (KIND=DP), INTENT (OUT) ::                                      &
                    A00, A01, A02, A10, A11, A12, A20, A21, A22

  REAL (KIND=DP)              :: QX, QT, QXT, RX, RT, RXT, NI, G


  NI = REAL(VELH ,DP)*REAL(DT,DP)/REAL(H,DP)
  G  = REAL(VELCE,DP)*REAL(DT,DP)/REAL(H,DP)

  QX  = ( WB*(1._DP + NI) -    NI )/(1._DP + NI)/(1._DP - WB)
  RX  = ( 1._DP - G )/( 1._DP + G )

  QT  = ( WB*(1._DP + NI) - 1._DP )/(1._DP + NI)/(1._DP - WB)
  RT  =-( 1._DP - G )/( 1._DP + G )

  QXT =   WB/(WB - 1._DP)
  RXT =  -1._DP

  A00 =   0._DP

  A01 = -(QX + RX)
  A10 = -(QT + RT)

  A02 = - QX * RX
  A20 = - QT * RT

  A11 = -(QX*RT  + QT*RX + QXT + RXT)

  A12 = -(QX*RXT + RX*QXT)
  A21 = -(QT*RXT + RT*QXT)

  A22 = - QXT* RXT

END SUBROUTINE  HIGCE2_FD

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
!-----------------------------------------------------------------------
! subroutine calculates some expresions needed in numerical integration
!-----------------------------------------------------------------------
SUBROUTINE PREPARE_INT
  USE NRTYPE         , ONLY: WP
  USE MOD_MODEL      , ONLY: H, HF, MXF, MYF, NK
  USE MOD_CONST      , ONLY: I2N_FREQ1
  USE MOD_INTEGRATION, ONLY: PTS, PTS2, PTS_3, PTS_2, PTS_1,           &
                              MU_FINE,  LA_FINE, RHO_FINE,             &
                             YMU_FINE, YLA_FINE


!----------------------------------------------------------------------

  IMPLICIT NONE

!----------------------------------------------------------------------

  PTS2 = PTS/2
  PTS_1= REAL(PTS,WP)
  PTS_3= PTS_1*PTS_1*PTS_1
  PTS_2= PTS_1*PTS_1

  ALLOCATE (                                                           &
         MU_FINE ( MXF*PTS+PTS2, NK*PTS+PTS2, PTS+PTS2            ),   &
         LA_FINE ( MXF*PTS+PTS2, NK*PTS+PTS2, PTS+PTS2            ),   &
        RHO_FINE ( MXF*PTS+PTS2, NK*PTS+PTS2, PTS+PTS2            ),   &
        YMU_FINE ( MXF*PTS+PTS2, NK*PTS+PTS2, PTS+PTS2, I2N_FREQ1 ),   &
        YLA_FINE ( MXF*PTS+PTS2, NK*PTS+PTS2, PTS+PTS2, I2N_FREQ1 ) )

END SUBROUTINE PREPARE_INT

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

SUBROUTINE  FDL2MT ( SFIS,SDEL,SLAM,M0, SXX,SXY,SXZ,SYY,SYZ,SZZ )
  USE PRECISION   , ONLY:                                              &
                          DP, PP
  USE GRID_MEDIUM , ONLY:                                              &
                          H
  USE AUXIL       , ONLY:                                              &
                          PI
  USE CONTROL_DATA, ONLY:                                              &
                          DT

  IMPLICIT NONE


  REAL (KIND=PP), INTENT(IN)  :: SFIS, SDEL, SLAM, M0
  REAL (KIND=PP), INTENT(OUT) :: SXX, SXY, SXZ, SYY, SYZ, SZZ

  REAL (KIND=DP):: SF, SD, SL, S2F, S2D, HF1,  MXX, MXY, MXZ,          &
                   CF, CD, CL, C2F, C2D, HF4,  MYY, MYZ, MZZ,          &
                   SFISR, SDELR, SLAMR


  IF ( ( SFIS < 0._PP ) .AND.                                          &
       ( SDEL < 0._PP ) .AND.                                          &
       ( SLAM < 0._PP ) ) THEN

    MXX = 1._PP/3._PP
    MYY = 1._PP/3._PP
    MZZ = 1._PP/3._PP
    MXY = 0._PP
    MXZ = 0._PP
    MYZ = 0._PP

  ELSE

    SFISR = PI*SFIS/180._DP
    SDELR = PI*SDEL/180._DP
    SLAMR = PI*SLAM/180._DP

    SF = SIN (SFISR)
    SD = SIN (SDELR)
    SL = SIN (SLAMR)
    CF = COS (SFISR)
    CD = COS (SDELR)
    CL = COS (SLAMR)

    S2F = 2._DP*SF*CF
    S2D = 2._DP*SD*CD
    C2F = CF*CF - SF*SF
    C2D = CD*CD - SD*SD

    MXX = -( SD*CL*S2F +       S2D*SL*SF*SF )

    MXY =  ( SD*CL*C2F + .5_DP*S2D*SL*S2F   )

    MXZ = -( CD*CL*CF  +       C2D*SL*SF    )

    MYY =  ( SD*CL*S2F -       S2D*SL*CF*CF )

    MYZ = -( CD*CL*SF  -       C2D*SL*CF    )

    MZZ =  S2D*SL

  END IF

  HF1 = M0/(      H*H*H)
  HF4 = M0/(4._DP*H*H*H)

  SXX = HF1*MXX
  SXY = HF4*MXY
  SXZ = HF4*MXZ
  SYY = HF1*MYY
  SYZ = HF4*MYZ
  SZZ = HF1*MZZ

END SUBROUTINE  FDL2MT
!=======================================================================

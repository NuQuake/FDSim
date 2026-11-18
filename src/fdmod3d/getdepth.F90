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
! function returns depth of sediments/bedrock interface
! (i.e. interface NI) in point (x,y)
!-----------------------------------------------------------------------
FUNCTION GETDEPTH (X,Y,NINF)

  USE NRTYPE   , ONLY: WP
  USE MOD_MODEL, ONLY: XIMAX, YIMAX, XIMIN, YIMIN, ITF1, H10, NI,      &
                       ZHF, HF, H

!----------------------------------------------------------------------

  IMPLICIT NONE

  REAL(WP),INTENT(IN) :: X,Y
  INTEGER, intent(in) :: ninf
  REAL(WP)            :: GETDEPTH
  REAL(WP)            :: XA, YA
  INTEGER             :: XIND, YIND

!----------------------------------------------------------------------

  XA = X
  YA = Y

  IF ( XA > XIMAX ) THEN
    XA = XIMAX
  END IF

  IF ( YA > YIMAX ) THEN
    YA = YIMAX
  END IF

  IF ( XA < XIMIN ) THEN
    XA = XIMIN
  END IF

  IF ( YA < YIMIN ) THEN
    YA = YIMIN
  END IF

  XIND = INT( (XA-XIMIN) / H10 + 1.5_WP )
  YIND = INT( (YA-YIMIN) / H10 + 1.5_WP )

  GETDEPTH=REAL(ITF1(XIND,YIND,NInf),WP)

  IF ( ZHF == 0._WP ) THEN
!!!    IF ( GETDEPTH > -2.*H  ) GETDEPTH = 0._WP
  ELSE
!!!    IF ( GETDEPTH > -2.*HF ) GETDEPTH = 0._WP
  END IF

END FUNCTION GETDEPTH

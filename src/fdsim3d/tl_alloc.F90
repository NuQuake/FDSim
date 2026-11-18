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

SUBROUTINE  TL_ALLOC

  USE GRID_MEDIUM , ONLY:                                              &
                          JM
  USE CONTROL_DATA, ONLY:                                              &
                          PWINC
  USE AUXIL       , ONLY:                                              &
                          X_MIN_D, X_MAX_D,                            &
                          Y_MIN_D, Y_MAX_D,                            &
                          Z_MIN_D, Z_MAX_D, TXZ_RES, TYZ_RES,          &
                          TXX_TOT, TYY_TOT, TZZ_TOT
  USE STRESSFIELD , ONLY:                                              &
                  TXX , TYY , TZZ , TXY ,       TXZ ,       TYZ
  USE WAVEFIELD   , ONLY:                                              &
                    UM, VM, WM,                                        &
                           XXX,   XYY,   XZZ,   XXZ,   XYZ,   XXY

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: ALLOSTAT

!-----------------------------------------------------------------------

  IF ( PWINC == 1 ) THEN
    ALLOCATE ( TXZ_RES (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D),             &
               TYZ_RES (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D) )
    TXZ_RES = 0.
    TYZ_RES = 0.
  END IF

  IF ( PWINC == 2 ) THEN
    ALLOCATE ( TXX_TOT (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D),             &
               TYY_TOT (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D),             &
               TZZ_TOT (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D) )
    TXX_TOT = 0.
    TYY_TOT = 0.
    TZZ_TOT = 0.
  END IF

  ALLOCATE (  UM  (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D), &
              VM  (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D), &
              WM  (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D), &
              XXX (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D), &
              XYY (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D), &
              XZZ (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D), &
              XXZ (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D), &
              XYZ (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D), &
              XXY (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D), &
              TXX (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D), &
              TYY (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D), &
              TZZ (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D), &
              TXZ (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D), &
              TYZ (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D), &
              TXY (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D), &
                                                       STAT=ALLOSTAT )
    IF  ( ALLOSTAT > 0 )  THEN
      WRITE (11,*) '  ALLOCATION ERROR T_ALLOC_02, STAT= ', ALLOSTAT
      STOP
    END IF

  UM  (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  VM  (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  WM  (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  XXX (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  XYY (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  XZZ (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  XXZ (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  XYZ (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  XXY (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  TXX (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  TYY (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  TZZ (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  TXZ (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  TYZ (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  TXY (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0.


END  SUBROUTINE  TL_ALLOC

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

SUBROUTINE  TL_ALLOCF

  USE AUXIL       , ONLY:                                              &
                          X_MIN_DF, X_MAX_DF,                          &
                          Y_MIN_DF, Y_MAX_DF,                          &
                          Z_MIN_DF, Z_MAX_DF
  USE STRESSFIELD , ONLY:                                              &
                  TXXF, TYYF, TZZF, TXYF,       TXZF,       TYZF
  USE WAVEFIELD   , ONLY:                                              &
                    UMF, VMF, WMF,                                     &
                           XXXF,  XYYF,  XZZF,  XXZF,  XYZF,  XXYF

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: ALLOSTAT

!-----------------------------------------------------------------------

  ALLOCATE                                                             &
    (  UMF  (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF), &
       VMF  (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF), &
       WMF  (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF), &
       XXXF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF), &
       XYYF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF), &
       XZZF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF), &
       XXZF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF), &
       XYZF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF), &
       XXYF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF), &
       TXXF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF), &
       TYYF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF), &
       TZZF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF), &
       TXZF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF), &
       TYZF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF), &
       TXYF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF), &
                                                       STAT=ALLOSTAT )
    IF  ( ALLOSTAT > 0 )  THEN
      WRITE (11,*) '  ALLOCATION ERROR T_ALLOC_02F, STAT= ', ALLOSTAT
      STOP
    END IF

  UMF  (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0.
  VMF  (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0.
  WMF  (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0.
  XXXF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0.
  XYYF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0.
  XZZF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0.
  XXZF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0.
  XYZF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0.
  XXYF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0.
  TXXF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0.
  TYYF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0.
  TZZF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0.
  TXZF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0.
  TYZF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0.
  TXYF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0.


END  SUBROUTINE  TL_ALLOCF

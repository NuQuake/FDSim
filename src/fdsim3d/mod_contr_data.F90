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

MODULE CONTROL_DATA
  USE PRECISION, ONLY: PP

  IMPLICIT NONE

  CHARACTER(LEN=20) :: SR_FILE_NAME,                                   &
                      MO_FILE_NAME, JMH_FILE_NAME, Q_FILE_NAME,        &
                                   JMHF_FILE_NAME

  CHARACTER(LEN=6), DIMENSION(:), ALLOCATABLE :: REC_NAME

  LOGICAL       ::  KEY_TLV, KEY_SNV,                                  &
                    KEY_ONE_STF, KEY_EN, STRESS_IMAGING

  INTEGER       ::  MX, MY, MZ,                                        &
                    MT1, MT2, IPAS2, MR, NPS, MSNP,                    &
                    TPML, PWINC,                                       &
                    LPAS, REFINE, LIN, N_FREQ

  INTEGER       ::  KTLE, KTRI, KTRE, KTFR, KTBO

  REAL (KIND=PP)::  DT, FRJMAX, FRANGE, TAU_EPS, XBMIN, YBMIN

  REAL (KIND=PP)::  SFIS, SDEL, SLAM,  M0, CU, CV

  REAL (KIND=PP)::  R, TAU, POW, SF, WC

  REAL (KIND=PP)::  OMG, WB,                                           &
                    THPPLE, THPPRI, THPPRE, THPPFR, THPPBO,            &
                    THSSLE, THSSRI, THSSRE, THSSFR, THSSBO

  REAL (KIND=PP), DIMENSION (:) , ALLOCATABLE ::                       &
                    IREC, KREC, LREC

  INTEGER, DIMENSION (:) , ALLOCATABLE ::                              &
                    IS, KS, LS, ITSB, ITSE

END MODULE CONTROL_DATA

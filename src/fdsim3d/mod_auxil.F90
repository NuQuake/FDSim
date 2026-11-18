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

MODULE AUXIL
  USE PRECISION, ONLY: PP, DP
  USE NONREF_BOUND, ONLY: NR_DATA_X, NR_DATA_Y

!----------------------------------------------------------------------

  IMPLICIT NONE

!----------------------------------------------------------------------

  INTEGER       , PARAMETER :: LLAN = 2, ORDER = 2 ! FOR VALUE 4: STORE_SEISMOGRAMS HAS TO BE CHANGED !!!

#if(SECOND)
  REAL (KIND=PP), PARAMETER :: PI=3.141592653589793,                   &
                               A = 0., B = 1.
#else
  REAL (KIND=PP), PARAMETER :: PI=3.141592653589793,                   &
                               A = -1./24., B = 9./8.
#endif

  INTEGER       ::  MXF, MYF

  INTEGER       ::  NPST, REF2

  INTEGER       ::  PMLL, PMLXH, PMLYH, PMLZH, PMLZH1, PMLZH2,         &
                           X_MIN_D , X_MAX_D , X_MIN , X_MAX ,         &
                           Y_MIN_D , Y_MAX_D , Y_MIN , Y_MAX ,         &
                           Z_MIN_D , Z_MAX_D , Z_MIN , Z_MAX ,         &
                           X_MIN_DF, X_MAX_DF, X_MINF, X_MAXF,         &
                           Y_MIN_DF, Y_MAX_DF, Y_MINF, Y_MAXF,         &
                           Z_MIN_DF, Z_MAX_DF, Z_MINF, Z_MAXF,         &
                    MX1, MY1,                                          &
                    PMLLF, PMLXHF, PMLYHF, MXF1, MYF1, TPMLF,          &
                           source_id0A , dest_id0A ,                   &
                           source_id0B , dest_id0B ,                   &
                           source_id1A , dest_id1A ,                   &
                           source_id1B , dest_id1B ,                   &
                           source_id2A , dest_id2A ,                   &
                           source_id2B , dest_id2B ,                   &
                           source_id2AF, dest_id2AF,                   &
                           source_id2BF, dest_id2BF

  REAL (KIND=PP)::  AH, BH, CH, HF, AHF, BHF, CHF, H32, H2, HF3, H3

  REAL (KIND=PP)::  SXX, SXY, SXZ, SYY, SYZ, SZZ, SOURTF

  REAL (KIND=PP), DIMENSION (:), ALLOCATABLE ::                        &
                    SXXT, SXYT, SXZT, SYYT, SYZT, SZZT,                &
                    BUFSF, BUFS

  LOGICAL , DIMENSION (:), ALLOCATABLE ::  LOGSNP

  INTEGER , DIMENSION (:), ALLOCATABLE ::  SRC_P

  REAL (KIND=PP), DIMENSION (4) ::                                     &
                    OTA, OTB, OTC, OTD

  REAL (KIND=PP), DIMENSION (:  ), ALLOCATABLE :: CUS, CVS, CWS,       &
                                                  SRCDELI

  REAL (KIND=PP), DIMENSION (:,:), ALLOCATABLE :: SOURTF_A,            &
                                                  TMP, LANC, VARF1,    &
                                                  TXZ_RES, TYZ_RES,    &
                                         TXX_TOT, TYY_TOT, TZZ_TOT
  REAL (KIND=PP), DIMENSION (:,:), ALLOCATABLE :: BUFXS_P, BUFXR_P,    &
                                                  BUFYS_P, BUFYR_P
  REAL (KIND=PP), DIMENSION (:,:,:), ALLOCATABLE::BUFXS ,BUFYS ,BUFZS ,&
                                                  BUFXSF,BUFYSF,BUFZSF,&
                                                  BUFXR ,BUFYR ,BUFZR ,&
                                                  BUFXRF,BUFYRF,BUFZRF

END MODULE AUXIL

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

SUBROUTINE  SRC_V ( IC, KC, LC, MXL, MXH, MYL, MYH, MZL, MZH,          &
                    UT, VT, WT, JMT, DENUT, DENVT, DENWT, REFINE )

  USE PRECISION   , ONLY :                                             &
                           PP, IK
  USE GRID_MEDIUM , ONLY :                                             &
                           JMNUM
  USE AUXIL       , ONLY :                                             &
                           SXX, SXY, SXZ, SYY, SYZ, SZZ, SOURTF

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,                                         INTENT(INOUT) ::    &
                                                     IC, KC, LC
  INTEGER,                                         INTENT(IN)    ::    &
                                                     MXL, MXH,         &
                                                     MYL, MYH,         &
                                                     MZL, MZH, REFINE
  INTEGER(IK), DIMENSION(MXL:MXH,MYL:MYH,MZL:MZH), INTENT(INOUT) ::    &
                                                                    JMT
  REAL   (PP), DIMENSION(1:JMNUM)                , INTENT(INOUT) ::    &
                                                    DENUT, DENVT, DENWT
  REAL   (PP), DIMENSION(MXL:MXH,MYL:MYH,MZL:MZH), INTENT(INOUT) ::    &
                                                             UT, VT, WT

  INTEGER :: JM1, JM2, JM3, JM4, JM5, JM6, R, S, T,                    &
             ICPS, ICMT, ICPR, ICMR,                                   &
             KCPS, KCMT, KCPR, KCMR,                                   &
             LCPS, LCMT, LCPR, LCMR

  REAL(PP):: CICPS, CICMT, CICPR, CICMR,                               &
             CKCPS, CKCMT, CKCPR, CKCMR,                               &
             CLCPS, CLCMT, CLCPR, CLCMR, CIC, CKC, CLC

!----------------------------------------------------------------------

  R =  REFINE
  S = (REFINE+1)/2
  T = (REFINE-1)/2

  ICPS  = IC + S
  ICMT  = IC - T
  ICPR  = IC + R
  ICMR  = IC - R
  CIC   = 1.
  CICPS = 1.
  CICMT = 1.
  CICPR = 1.
  CICMR = 1.
  IF ( IC > MXH ) THEN
    IC  = MXH
    CIC = 0.
  END IF
  IF ( IC < MXL ) THEN
    IC = MXL
    CIC = 0.
  END IF
  IF ( ICPS > MXH ) THEN
    ICPS = IC
    CICPS = 0.
  END IF
  IF ( ICMT < MXL ) THEN
    ICMT = IC
    CICMT = 0.
  END IF
  IF ( ICPR > MXH ) THEN
    ICPR = IC
    CICPR = 0.
  END IF
  IF ( ICMR < MXL ) THEN
    ICMR = IC
    CICMR = 0.
  END IF

  KCPS  = KC + S
  KCMT  = KC - T
  KCPR  = KC + R
  KCMR  = KC - R
  CKC   = 1.
  CKCPS = 1.
  CKCMT = 1.
  CKCPR = 1.
  CKCMR = 1.
  IF ( KC > MYH ) THEN
    KC  = MYH
    CKC = 0.
  END IF
  IF ( KC < MYL ) THEN
    KC = MYL
    CKC = 0.
  END IF
  IF ( KCPS > MYH ) THEN
    KCPS = KC
    CKCPS = 0.
  END IF
  IF ( KCMT < MYL ) THEN
    KCMT = KC
    CKCMT = 0.
  END IF
  IF ( KCPR > MYH ) THEN
    KCPR = KC
    CKCPR = 0.
  END IF
  IF ( KCMR < MYL ) THEN
    KCMR = KC
    CKCMR = 0.
  END IF

  LCPS  = LC + S
  LCMT  = LC - T
  LCPR  = LC + R
  LCMR  = LC - R
  CLC   = 1.
  CLCPS = 1.
  CLCMT = 1.
  CLCPR = 1.
  CLCMR = 1.
  IF ( LC > MZH ) THEN
    LC  = MZH
    CLC = 0.
  END IF
  IF ( LC < MZL ) THEN
    LC = MZL
    CLC = 0.
  END IF
  IF ( LCPS > MZH ) THEN
    LCPS = LC
    CLCPS = 0.
  END IF
  IF ( LCMT < MZL ) THEN
    LCMT = LC
    CLCMT = 0.
  END IF
  IF ( LCPR > MZH ) THEN
    LCPR = LC
    CLCPR = 0.
  END IF
  IF ( LCMR < MZL ) THEN
    LCMR = LC
    CLCMR = 0.
  END IF

  JM1   =                 JMT ( ICPS,KC  , LC   )
  UT ( ICPS,KC  , LC   ) = UT ( ICPS,KC  , LC   )                      &
                         + SXX * SOURTF * DENUT (JM1)*CICPS*CKC*CLC
  JM2   =                 JMT ( ICPS,KCPR, LC   )
  UT ( ICPS,KCPR, LC   ) = UT ( ICPS,KCPR, LC   )                      &
                         + SXY * SOURTF * DENUT (JM2)*CICPS*CKCPR*CLC
  JM3   =                 JMT ( ICMT,KCPR, LC   )
  UT ( ICMT,KCPR, LC   ) = UT ( ICMT,KCPR, LC   )                      &
                         + SXY * SOURTF * DENUT (JM3)*CICMT*CKCPR*CLC

  JM4   =                 JMT ( ICMT,KC  , LC   )
  UT ( ICMT,KC  , LC   ) = UT ( ICMT,KC  , LC   )                      &
                         - SXX * SOURTF * DENUT (JM4)*CICMT*CKC*CLC
  JM5   =                 JMT ( ICPS,KCMR, LC   )
  UT ( ICPS,KCMR, LC   ) = UT ( ICPS,KCMR, LC   )                      &
                         - SXY * SOURTF * DENUT (JM5)*CICPS*CKCMR*CLC
  JM6   =                 JMT ( ICMT,KCMR, LC   )
  UT ( ICMT,KCMR, LC   ) = UT ( ICMT,KCMR, LC   )                      &
                         - SXY * SOURTF * DENUT (JM6)*CICMT*CKCMR*CLC


  JM1   =  JMT ( IC  , KCPS,  LC )
  JM2   =  JMT ( ICPR, KCPS,  LC )
  JM3   =  JMT ( ICPR, KCMT,  LC )

  JM4   =  JMT ( IC  , KCMT,  LC )
  JM5   =  JMT ( ICMR, KCPS,  LC )
  JM6   =  JMT ( ICMR, KCMT,  LC )

  VT ( IC  ,KCPS, LC   ) = VT ( IC  ,KCPS, LC   )                      &
                         + SYY * SOURTF * DENVT (JM1)*CKCPS*CIC*CLC
  VT ( ICPR,KCPS, LC   ) = VT ( ICPR,KCPS, LC   )                      &
                         + SXY * SOURTF * DENVT (JM2)*CICPR*CKCPS*CLC
  VT ( ICPR,KCMT, LC   ) = VT ( ICPR,KCMT, LC   )                      &
                         + SXY * SOURTF * DENVT (JM3)*CICPR*CKCMT*CLC

  VT ( IC  ,KCMT, LC   ) = VT ( IC  ,KCMT, LC   )                      &
                         - SYY * SOURTF * DENVT (JM4)*CKCMT*CIC*CLC
  VT ( ICMR,KCPS, LC   ) = VT ( ICMR,KCPS, LC   )                      &
                         - SXY * SOURTF * DENVT (JM5)*CICMR*CKCPS*CLC
  VT ( ICMR,KCMT, LC   ) = VT ( ICMR,KCMT, LC   )                      &
                         - SXY * SOURTF * DENVT (JM6)*CICMR*CKCMT*CLC


  JM1   =  JMT ( ICPR, KC  ,  LCMT )
  JM3   =  JMT ( IC  , KCPR,  LCMT )

  JM4   =  JMT ( IC  , KC  ,  LCMT )
  JM5   =  JMT ( ICMR, KC  ,  LCMT )
  JM6   =  JMT ( IC  , KCMR,  LCMT )

  WT ( ICPR,KC  , LCMT ) = WT ( ICPR,KC  , LCMT )                      &
                         + SXZ * SOURTF * DENWT (JM1)*CICPR*CLCMT*CKC
  WT ( IC  ,KCPR, LCMT ) = WT ( IC  ,KCPR, LCMT )                      &
                         + SYZ * SOURTF * DENWT (JM3)*CKCPR*CLCMT*CIC

  WT ( IC  ,KC  , LCMT ) = WT ( IC  ,KC  , LCMT )                      &
                         - SZZ * SOURTF * DENWT (JM4)*CLCMT*CIC*CKC
  WT ( ICMR,KC  , LCMT ) = WT ( ICMR,KC  , LCMT )                      &
                         - SXZ * SOURTF * DENWT (JM5)*CICMR*CLCMT*CKC
  WT ( IC  ,KCMR, LCMT ) = WT ( IC  ,KCMR, LCMT )                      &
                         - SYZ * SOURTF * DENWT (JM6)*CKCMR*CLCMT*CIC


  JM1   =  JMT ( ICPS, KC  ,  LCMR )
  JM3   =  JMT ( ICMT, KC  ,  LCMR )

  JM4   =  JMT ( IC  , KCPS,  LCMR )
  JM6   =  JMT ( IC  , KCMT,  LCMR )

  UT ( ICPS,KC  , LCMR ) = UT ( ICPS,KC  , LCMR )                      &
                         - SXZ * SOURTF * DENUT (JM1)*CICPS*CLCMR*CKC
  UT ( ICMT,KC  , LCMR ) = UT ( ICMT,KC  , LCMR )                      &
                         - SXZ * SOURTF * DENUT (JM3)*CICMT*CLCMR*CKC

  VT ( IC  ,KCPS, LCMR ) = VT ( IC  ,KCPS, LCMR )                      &
                         - SYZ * SOURTF * DENVT (JM4)*CKCPS*CLCMR*CIC
  VT ( IC  ,KCMT, LCMR ) = VT ( IC  ,KCMT, LCMR )                      &
                         - SYZ * SOURTF * DENVT (JM6)*CKCMT*CLCMR*CIC


  JM1   =  JMT ( ICPS, KC  ,  LCPR )
  JM3   =  JMT ( ICMT, KC  ,  LCPR )

  JM4   =  JMT ( IC  , KCPS,  LCPR )
  JM6   =  JMT ( IC  , KCMT,  LCPR )

  UT ( ICPS,KC  , LCPR ) = UT ( ICPS,KC  , LCPR )                      &
                         + SXZ * SOURTF * DENUT (JM1)*CICPS*CLCPR*CKC
  UT ( ICMT,KC  , LCPR ) = UT ( ICMT,KC  , LCPR )                      &
                         + SXZ * SOURTF * DENUT (JM3)*CICMT*CLCPR*CKC

  VT ( IC  ,KCPS, LCPR ) = VT ( IC  ,KCPS, LCPR )                      &
                         + SYZ * SOURTF * DENVT (JM4)*CKCPS*CLCPR*CIC
  VT ( IC  ,KCMT, LCPR ) = VT ( IC  ,KCMT, LCPR )                      &
                         + SYZ * SOURTF * DENVT (JM6)*CKCMT*CLCPR*CIC


  JM1   =  JMT ( IC  , KC  ,  LCPS )
  JM2   =  JMT ( ICPR, KC  ,  LCPS )
  JM3   =  JMT ( IC  , KCPR,  LCPS )

  JM4   =  JMT ( ICMR, KC  ,  LCPS )
  JM6   =  JMT ( IC  , KCMR,  LCPS )

  WT ( IC  ,KC  , LCPS ) = WT ( IC  ,KC  , LCPS )                      &
                         + SZZ * SOURTF * DENWT (JM1)*CLCPS*CIC*CKC
  WT ( ICPR,KC  , LCPS ) = WT ( ICPR,KC  , LCPS )                      &
                         + SXZ * SOURTF * DENWT (JM2)*CICPR*CLCPS*CKC
  WT ( IC  ,KCPR, LCPS ) = WT ( IC  ,KCPR, LCPS )                      &
                         + SYZ * SOURTF * DENWT (JM3)*CKCPR*CLCPS*CIC
  WT ( ICMR,KC  , LCPS ) = WT ( ICMR,KC  , LCPS )                      &
                         - SXZ * SOURTF * DENWT (JM4)*CICMR*CLCPS*CKC
  WT ( IC  ,KCMR, LCPS ) = WT ( IC  ,KCMR, LCPS )                      &
                         - SYZ * SOURTF * DENWT (JM6)*CKCMR*CLCPS*CIC


END SUBROUTINE  SRC_V

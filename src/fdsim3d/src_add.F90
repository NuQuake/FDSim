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

SUBROUTINE  SRC_ADD ( TILE )

  USE PRECISION   , ONLY:                                              &
                          PP, IK
  USE CONTROL_DATA, ONLY:                                              &
                          MZ, ITSB, ITSE, KEY_ONE_STF, IS, KS, LS,     &
                          TPML, PWINC, NPS, LPAS, REFINE, MT1
  USE GRID_MEDIUM , ONLY:                                              &
                          DENU , DENV , DENW , JM ,                    &
                          DENUF, DENVF, DENWF, JMF
  USE WAVEFIELD   , ONLY:                                              &
                          UM, VM, WM, UMF, VMF, WMF
  USE AUXIL       , ONLY:                                              &
                          SXX , SXY , SXZ , SYY , SYZ , SZZ ,          &
                          SXXT, SXYT, SXZT, SYYT, SYZT, SZZT,          &
                          SOURTF, SOURTF_A,                            &
                          X_MIN_DF, X_MAX_DF, X_MINF, X_MAXF,          &
                          Y_MIN_DF, Y_MAX_DF, Y_MINF, Y_MAXF,          &
                          Z_MIN_DF, Z_MAX_DF, Z_MINF, Z_MAXF,          &
                          X_MIN_D , X_MAX_D , X_MIN , X_MAX ,          &
                          Y_MIN_D , Y_MAX_D , Y_MIN , Y_MAX ,          &
                          Z_MIN_D , Z_MAX_D , Z_MIN , Z_MAX ,          &
                          CUS, CVS, CWS, SRCDELI, NPST, SRC_P

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: TILE

  INTEGER             :: ICS, KCS, LCS, SRC, SRCT, TILAU, I

#if(NOISE)
  REAL(KIND=PP)        :: CU, CV, CW
#endif
!-----------------------------------------------------------------------

  INTERFACE

    SUBROUTINE    SRC_V ( IC, KC, LC, MXL, MXH, MYL, MYH, MZL, MZH,    &
                          UT, VT, WT, JMT, DENUT, DENVT, DENWT, REFINE )
      USE PRECISION  , ONLY: PP, IK
      USE GRID_MEDIUM, ONLY: JMNUM
      IMPLICIT NONE
      INTEGER,                                         INTENT(INOUT) ::&
                                                         IC, KC, LC
      INTEGER,                                         INTENT(IN)    ::&
                                                        MXL, MXH,      &
                                                        MYL, MYH,      &
                                                        MZL, MZH, REFINE
      INTEGER(IK), DIMENSION(MXL:MXH,MYL:MYH,MZL:MZH), INTENT(INOUT) ::&
                                                                    JMT
      REAL   (PP), DIMENSION(1:JMNUM)                , INTENT(INOUT) ::&
                                                    DENUT, DENVT, DENWT
      REAL   (PP), DIMENSION(MXL:MXH,MYL:MYH,MZL:MZH), INTENT(INOUT) ::&
                                                             UT, VT, WT
    END SUBROUTINE  SRC_V

  END INTERFACE

!-----------------------------------------------------------------------

  DO SRC = 1, NPS
    IF( (ITSB(SRC) <= TILE) .AND. (TILE <= ITSE(SRC)) .AND.            &
        (ITSB(SRC) >= MT1 )  ) THEN
      TILAU = TILE - ITSB(SRC) + 1
      ICS = IS  (SRC)
      KCS = KS  (SRC)
      LCS = LS  (SRC)
#if(NOISE)
      IF (ITSE(SRC) == 999999 ) THEN
        SOURTF = SRCDELI(TILAU)
      ELSE
        DO I = 1, NPST
          IF (SRC_P(I) == SRC) THEN
            SOURTF = SOURTF_A(I,TILAU)
          END IF
        END DO
      END IF

      CU  = CUS (SRC)
      CV  = CVS (SRC)
      CW  = CWS (SRC)

      IF ( ( X_MIN  -1 <= ICS ) .AND. ( ICS <= X_MAX +1) .AND.         &
           ( Y_MIN  -1 <= KCS ) .AND. ( KCS <= Y_MAX +1) .AND.         &
           ( Z_MIN  -1 <= LCS ) .AND. ( LCS <= Z_MAX +1) ) THEN
        UM ( ICS+1,KCS  ,LCS   ) = UM ( ICS+1,KCS  ,LCS   )            &
                         + CU  * SOURTF
        UM ( ICS  ,KCS  ,LCS   ) = UM ( ICS  ,KCS  ,LCS   )            &
                         + CU  * SOURTF
        VM ( ICS  ,KCS+1,LCS   ) = VM ( ICS  ,KCS+1,LCS   )            &
                         + CV  * SOURTF
        VM ( ICS  ,KCS  ,LCS   ) = VM ( ICS  ,KCS  ,LCS   )            &
                         + CV  * SOURTF
        WM ( ICS  ,KCS  ,LCS+1 ) = WM ( ICS  ,KCS  ,LCS+1 )            &
                         + CW  * SOURTF
        WM ( ICS  ,KCS  ,LCS   ) = WM ( ICS  ,KCS  ,LCS   )            &
                         + CW  * SOURTF
      END IF

      IF ( LPAS /= 0 ) THEN
        IF ( ( X_MINF -1 <= ICS ) .AND. ( ICS <= X_MAXF +1) .AND.      &
             ( Y_MINF -1 <= KCS ) .AND. ( KCS <= Y_MAXF +1) .AND.      &
             ( Z_MINF -1 <= LCS ) .AND. ( LCS <= Z_MAXF +1) .AND.      &
             ( LCS <= LPAS ) ) THEN
          UMF( ICS+1,KCS  ,LCS   ) = UMF( ICS+1,KCS  ,LCS   )          &
                         + CU  * SOURTF
          UMF( ICS  ,KCS  ,LCS   ) = UMF( ICS  ,KCS  ,LCS   )          &
                         + CU  * SOURTF
          VMF( ICS  ,KCS+1,LCS   ) = VMF( ICS  ,KCS+1,LCS   )          &
                         + CV  * SOURTF
          VMF( ICS  ,KCS  ,LCS   ) = VMF( ICS  ,KCS  ,LCS   )          &
                         + CV  * SOURTF
          WMF( ICS  ,KCS  ,LCS+1 ) = WMF( ICS  ,KCS  ,LCS+1 )          &
                         + CW  * SOURTF
          WMF( ICS  ,KCS  ,LCS   ) = WMF( ICS  ,KCS  ,LCS   )          &
                         + CW  * SOURTF
        END IF
      END IF
#else
      IF ( KEY_ONE_STF ) THEN
        SRCT = 1
      ELSE
        SRCT = SRC
      END IF

      SOURTF = SOURTF_A(SRCT,TILAU)
      SXX = SXXT(SRC)
      SXY = SXYT(SRC)
      SXZ = SXZT(SRC)
      SYY = SYYT(SRC)
      SYZ = SYZT(SRC)
      SZZ = SZZT(SRC)

      IF      ( PWINC == 0 ) THEN

        IF ( LCS < LPAS ) THEN
          LCS = NINT((REAL(LCS-LPAS)+0.5)/REFINE+REAL(LPAS)-0.5)
          ICS = (ICS-1)/REFINE+1
          KCS = (KCS-1)/REFINE+1
        END IF

        IF ( ( X_MIN  -1 <= ICS ) .AND. ( ICS <= X_MAX +1) .AND.       &
             ( Y_MIN  -1 <= KCS ) .AND. ( KCS <= Y_MAX +1) .AND.       &
             ( Z_MIN  -1 <= LCS ) .AND. ( LCS <= Z_MAX +1) ) THEN
          CALL SRC_V ( IC = ICS, KC = KCS, LC = LCS,                   &
                       MXL = X_MIN_D, MXH = X_MAX_D,                   &
                       MYL = Y_MIN_D, MYH = Y_MAX_D,                   &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,                   &
                       UT = UM, VT = VM, WT = WM, JMT = JM,            &
                       DENUT = DENU, DENVT = DENV, DENWT = DENW,       &
                       REFINE = 1 )
        END IF

        ICS = IS  (SRC)
        KCS = KS  (SRC)
        LCS = LS  (SRC)

        IF ( LCS == LPAS ) THEN
          ICS = (ICS-1)*REFINE+1
          KCS = (KCS-1)*REFINE+1
        END IF

        IF ( LPAS /= 0 ) THEN
          IF ( ( X_MINF -REFINE <= ICS ) .AND. ( ICS <= X_MAXF +REFINE) .AND.    &
               ( Y_MINF -REFINE <= KCS ) .AND. ( KCS <= Y_MAXF +REFINE) .AND.    &
               ( Z_MINF -REFINE <= LCS ) .AND. ( LCS <= Z_MAXF +REFINE) .AND.    &
               ( LCS <= LPAS ) ) THEN
            CALL SRC_V ( IC = ICS, KC = KCS, LC = LCS,                 &
                         MXL = X_MIN_DF, MXH = X_MAX_DF,               &
                         MYL = Y_MIN_DF, MYH = Y_MAX_DF,               &
                         MZL = Z_MIN_DF, MZH = Z_MAX_DF,               &
                         UT = UMF, VT = VMF, WT = WMF, JMT = JMF,      &
                         DENUT = DENU , DENVT = DENV , DENWT = DENW ,  &
                         REFINE = 1 )
          END IF
        END IF

!      ELSE IF ( PWINC == 1 ) THEN

      ELSE IF ( PWINC == 2 ) THEN

        IF ( ( Z_MIN <= LCS ) .AND. ( LCS <= Z_MAX ) ) THEN
          VM(:,:,LCS) = VM(:,:,LCS) + SOURTF
        END IF
        IF ( ( Z_MINF <= LCS ) .AND. ( LCS <= Z_MAXF ) ) THEN
          VMF(:,:,LCS) = VMF(:,:,LCS) + SOURTF
        END IF

      ELSE IF ( PWINC == 3 ) THEN

        IF ( ( Z_MIN <= LCS ) .AND. ( LCS <= Z_MAX ) ) THEN
          WM(:,:,LCS) = WM(:,:,LCS) + SOURTF
        END IF
        IF ( ( Z_MINF <= LCS ) .AND. ( LCS <= Z_MAXF ) ) THEN
          WMF(:,:,LCS) = WMF(:,:,LCS) + SOURTF
        END IF

      END IF
#endif

    END IF
  END DO


END SUBROUTINE  SRC_ADD

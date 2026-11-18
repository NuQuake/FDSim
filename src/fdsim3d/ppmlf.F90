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

SUBROUTINE PPMLF

  USE PRECISION   , ONLY:                                              &
                          PP
  USE CONTROL_DATA, ONLY:                                              &
                          DT, R, TAU, POW, SF, WC
  USE GRID_MEDIUM , ONLY:                                              &
                          JMF, L2MX, L2MY, L2MZ,                       &
                          MUXY, MUXZ,                                  &
                          DENUF, DENVF, DENWF

  USE PML         , ONLY:                                              &
            ODA_XLF,  ODB_XLF,  ODC_XLF,  ODA_XHF,  ODB_XHF,  ODC_XHF, &
            ODD_XLF,                      ODD_XHF,                     &
            ODA_YLF,  ODB_YLF,  ODC_YLF,  ODA_YHF,  ODB_YHF,  ODC_YHF, &
            ODD_YLF,                      ODD_YHF,                     &
           ODA2_XLF, ODB2_XLF, ODC2_XLF, ODA2_XHF, ODB2_XHF, ODC2_XHF, &
           ODD2_XLF,                     ODD2_XHF,                     &
           ODA2_YLF, ODB2_YLF, ODC2_YLF, ODA2_YHF, ODB2_YHF, ODC2_YHF, &
           ODD2_YLF,                     ODD2_YHF,                     &
                                        PZXLF, PZXHF,                  &
                                        PZYLF, PZYHF,                  &
                                        PXXLF, PXXHF,                  &
                                        PYYLF, PYYHF,                  &
                                        PYXLF, PYXHF,                  &
                                        PXYLF, PXYHF,                  &
                                          RXXLF, RXXHF,                &
                                          RXYLF, RXYHF,                &
                                          RYYLF, RYYHF,                &
                                          RYXLF, RYXHF,                &
                                          RZXLF, RZXHF,                &
                                          RZYLF, RZYHF

  USE AUXIL       , ONLY:                                              &
                           PMLLF, PMLXHF, PMLYHF,                      &
                           X_MIN_DF, X_MAX_DF, X_MINF, X_MAXF,         &
                           Y_MIN_DF, Y_MAX_DF, Y_MINF, Y_MAXF,         &
                           Z_MIN_DF, Z_MAX_DF, Z_MINF, Z_MAXF,         &
                           MXF, MYF, HF, TPMLF

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER        :: I, L, K, J, JM1

  REAL (KIND=PP) :: VEL, VELS, VELP, OMG, OMGPDT, AP, DENOM

!-----------------------------------------------------------------------

IF ( Z_MIN_DF <= Z_MAX_DF ) THEN
  ALLOCATE( ODA_XLF  (0:TPMLF  , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODB_XLF  (0:TPMLF  , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODC_XLF  (0:TPMLF  , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODD_XLF  (0:TPMLF  , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODA_XHF  (0:TPMLF  , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODB_XHF  (0:TPMLF  , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODC_XHF  (0:TPMLF  , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODD_XHF  (0:TPMLF  , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODA_YLF  (0:TPMLF  , X_MIN_DF:X_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODB_YLF  (0:TPMLF  , X_MIN_DF:X_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODC_YLF  (0:TPMLF  , X_MIN_DF:X_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODD_YLF  (0:TPMLF  , X_MIN_DF:X_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODA_YHF  (0:TPMLF  , X_MIN_DF:X_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODB_YHF  (0:TPMLF  , X_MIN_DF:X_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODC_YHF  (0:TPMLF  , X_MIN_DF:X_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODD_YHF  (0:TPMLF  , X_MIN_DF:X_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODA2_XLF (0:TPMLF+1, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODB2_XLF (0:TPMLF+1, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODC2_XLF (0:TPMLF+1, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODD2_XLF (0:TPMLF+1, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODA2_XHF (0:TPMLF+1, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODB2_XHF (0:TPMLF+1, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODC2_XHF (0:TPMLF+1, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODD2_XHF (0:TPMLF+1, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODA2_YLF (0:TPMLF+1, X_MIN_DF:X_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODB2_YLF (0:TPMLF+1, X_MIN_DF:X_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODC2_YLF (0:TPMLF+1, X_MIN_DF:X_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODD2_YLF (0:TPMLF+1, X_MIN_DF:X_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODA2_YHF (0:TPMLF+1, X_MIN_DF:X_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODB2_YHF (0:TPMLF+1, X_MIN_DF:X_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODC2_YHF (0:TPMLF+1, X_MIN_DF:X_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
            ODD2_YHF (0:TPMLF+1, X_MIN_DF:X_MAX_DF, Z_MIN_DF:Z_MAX_DF) )

  ALLOCATE (                                                           &
       PZXLF (PMLLF   :2       , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
       PZXHF (MXF     :PMLXHF  , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
       PZYLF (X_MIN_DF:X_MAX_DF, PMLLF   :2       , Z_MIN_DF:Z_MAX_DF),&
       PZYHF (X_MIN_DF:X_MAX_DF, MYF     :PMLYHF  , Z_MIN_DF:Z_MAX_DF),&
       PXXLF (PMLLF   :1       , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
       PXXHF (MXF     :PMLXHF  , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
       PYYLF (X_MIN_DF:X_MAX_DF, PMLLF   :1       , Z_MIN_DF:Z_MAX_DF),&
       PYYHF (X_MIN_DF:X_MAX_DF, MYF     :PMLYHF  , Z_MIN_DF:Z_MAX_DF),&
       PYXLF (PMLLF   :2       , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
       PYXHF (MXF     :PMLXHF  , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
       PXYLF (X_MIN_DF:X_MAX_DF, PMLLF   :2       , Z_MIN_DF:Z_MAX_DF),&
       PXYHF (X_MIN_DF:X_MAX_DF, MYF     :PMLYHF  , Z_MIN_DF:Z_MAX_DF) )

  ALLOCATE (                                                           &
       RXXLF (PMLLF   :3       , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
       RXXHF (MXF-1   :PMLXHF  , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
       RXYLF (X_MIN_DF:X_MAX_DF, PMLLF   :2       , Z_MIN_DF:Z_MAX_DF),&
       RXYHF (X_MIN_DF:X_MAX_DF, MYF-1   :PMLYHF  , Z_MIN_DF:Z_MAX_DF),&
       RYYLF (X_MIN_DF:X_MAX_DF, PMLLF   :3       , Z_MIN_DF:Z_MAX_DF),&
       RYYHF (X_MIN_DF:X_MAX_DF, MYF-1   :PMLYHF  , Z_MIN_DF:Z_MAX_DF),&
       RYXLF (PMLLF   :2       , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
       RYXHF (MXF-1   :PMLXHF  , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
       RZXLF (PMLLF   :2       , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
       RZXHF (MXF-1   :PMLXHF  , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF),&
       RZYLF (X_MIN_DF:X_MAX_DF, PMLLF   :2       , Z_MIN_DF:Z_MAX_DF),&
       RZYHF (X_MIN_DF:X_MAX_DF, MYF-1   :PMLYHF  , Z_MIN_DF:Z_MAX_DF) )

  PZXLF (PMLLF   :2       , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0._PP
  PZXHF (MXF     :PMLXHF  , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0._PP
  PZYLF (X_MIN_DF:X_MAX_DF, PMLLF   :2       , Z_MIN_DF:Z_MAX_DF) = 0._PP
  PZYHF (X_MIN_DF:X_MAX_DF, MYF     :PMLYHF  , Z_MIN_DF:Z_MAX_DF) = 0._PP
  PXXLF (PMLLF   :1       , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0._PP
  PXXHF (MXF     :PMLXHF  , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0._PP
  PYYLF (X_MIN_DF:X_MAX_DF, PMLLF   :1       , Z_MIN_DF:Z_MAX_DF) = 0._PP
  PYYHF (X_MIN_DF:X_MAX_DF, MYF     :PMLYHF  , Z_MIN_DF:Z_MAX_DF) = 0._PP
  PYXLF (PMLLF   :2       , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0._PP
  PYXHF (MXF     :PMLXHF  , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0._PP
  PXYLF (X_MIN_DF:X_MAX_DF, PMLLF   :2       , Z_MIN_DF:Z_MAX_DF) = 0._PP
  PXYHF (X_MIN_DF:X_MAX_DF, MYF     :PMLYHF  , Z_MIN_DF:Z_MAX_DF) = 0._PP

  RXXLF (PMLLF   :3       , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0._PP
  RXXHF (MXF-1   :PMLXHF  , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0._PP
  RXYLF (X_MIN_DF:X_MAX_DF, PMLLF   :2       , Z_MIN_DF:Z_MAX_DF) = 0._PP
  RXYHF (X_MIN_DF:X_MAX_DF, MYF-1   :PMLYHF  , Z_MIN_DF:Z_MAX_DF) = 0._PP
  RYYLF (X_MIN_DF:X_MAX_DF, PMLLF   :3       , Z_MIN_DF:Z_MAX_DF) = 0._PP
  RYYHF (X_MIN_DF:X_MAX_DF, MYF-1   :PMLYHF  , Z_MIN_DF:Z_MAX_DF) = 0._PP
  RYXLF (PMLLF   :2       , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0._PP
  RYXHF (MXF-1   :PMLXHF  , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0._PP
  RZXLF (PMLLF   :2       , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0._PP
  RZXHF (MXF-1   :PMLXHF  , Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) = 0._PP
  RZYLF (X_MIN_DF:X_MAX_DF, PMLLF   :2       , Z_MIN_DF:Z_MAX_DF) = 0._PP
  RZYHF (X_MIN_DF:X_MAX_DF, MYF-1   :PMLYHF  , Z_MIN_DF:Z_MAX_DF) = 0._PP

!                                                   HETEROGENEOUS MEDIUM
  IF ( TPMLF == 0 ) THEN
    OMG = 0._PP
  ELSE
    OMG = -TAU*LOG(R)*DT/REAL(TPMLF)/(REAL(TPMLF)**POW)/HF
  END IF

  DO  L = Z_MINF, Z_MAXF
!______________________________________________________ REAR TRANSPARENT

    IF ( Y_MINF <= 1 ) THEN
      DO  I = X_MINF, X_MAXF
        JM1 = JMF ( I  , Y_MINF, L )

        VELS = SQRT( MUXY(JM1) *HF*DENUF(JM1)/DT )
        VELP = SQRT( L2MY(JM1) *HF*DENVF(JM1)/DT )
        VEL  = VELP

        ODA_YLF  (:,I,L) = 1._PP
        ODA2_YLF (:,I,L) = 1._PP

        ODB_YLF  (:,I,L) = 1._PP
        ODB2_YLF (:,I,L) = 1._PP

        ODC_YLF  (:,I,L) = 0._PP
        ODC2_YLF (:,I,L) = 0._PP

#if(SIMPLE_PML)
        DO J = 1, TPMLF
          OMGPDT = OMG*VEL* REAL(J)     **POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA_YLF (J,I,L) =  2._PP                         / DENOM
          ODB_YLF (J,I,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC_YLF (J,I,L) = -2._PP*     OMGPDT             / DENOM

          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA2_YLF(J,I,L) =  2._PP                         / DENOM
          ODB2_YLF(J,I,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC2_YLF(J,I,L) = -2._PP*     OMGPDT             / DENOM

        END DO

        J = TPMLF+1
          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA2_YLF(J,I,L) =  2._PP                         / DENOM
          ODB2_YLF(J,I,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC2_YLF(J,I,L) = -2._PP*     OMGPDT             / DENOM

#else
        DO J = 1, TPMLF
          OMGPDT = OMG*VEL* REAL(J)     **POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA_YLF (J,I,L) = (DENOM+1.)/2./AP
          ODB_YLF (J,I,L) =  DENOM
          ODC_YLF (J,I,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA2_YLF(J,I,L) = (DENOM+1.)/2./AP
          ODB2_YLF(J,I,L) =  DENOM
          ODC2_YLF(J,I,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

        END DO

        J = TPMLF+1
          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA2_YLF(J,I,L) = (DENOM+1.)/2./AP
          ODB2_YLF(J,I,L) =  DENOM
          ODC2_YLF(J,I,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)
#endif

        ODD_YLF  (:,I,L) = (ODC_YLF (:,I,L)+2.)/(ODB_YLF (:,I,L)+1.)
        ODD2_YLF (:,I,L) = (ODC2_YLF(:,I,L)+2.)/(ODB2_YLF(:,I,L)+1.)

      END DO
    END IF

!_____________________________________________________ FRONT TRANSPARENT

    IF ( Y_MAXF >= MYF ) THEN
      DO  I = X_MINF, X_MAXF
        JM1 = JMF ( I  , Y_MAXF, L )

        VELS = SQRT( MUXY(JM1) *HF*DENUF(JM1)/DT )
        VELP = SQRT( L2MY(JM1) *HF*DENVF(JM1)/DT )
        VEL  = VELP

        ODA_YHF  (:,I,L) = 1._PP
        ODA2_YHF (:,I,L) = 1._PP

        ODB_YHF  (:,I,L) = 1._PP
        ODB2_YHF (:,I,L) = 1._PP

        ODC_YHF  (:,I,L) = 0._PP
        ODC2_YHF (:,I,L) = 0._PP

#if(SIMPLE_PML)
        DO J = 1, TPMLF
          OMGPDT = OMG*VEL* REAL(J)     **POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA_YHF (J,I,L) =  2._PP                         / DENOM
          ODB_YHF (J,I,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC_YHF (J,I,L) = -2._PP*     OMGPDT             / DENOM

          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA2_YHF(J,I,L) =  2._PP                         / DENOM
          ODB2_YHF(J,I,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC2_YHF(J,I,L) = -2._PP*     OMGPDT             / DENOM

        END DO

        J = TPMLF+1
          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA2_YHF(J,I,L) =  2._PP                         / DENOM
          ODB2_YHF(J,I,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC2_YHF(J,I,L) = -2._PP*     OMGPDT             / DENOM

#else
        DO J = 1, TPMLF
          OMGPDT = OMG*VEL* REAL(J)     **POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA_YHF (J,I,L) = (DENOM+1.)/2./AP
          ODB_YHF (J,I,L) =  DENOM
          ODC_YHF (J,I,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA2_YHF(J,I,L) = (DENOM+1.)/2./AP
          ODB2_YHF(J,I,L) =  DENOM
          ODC2_YHF(J,I,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

        END DO

        J = TPMLF+1
          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA2_YHF(J,I,L) = (DENOM+1.)/2./AP
          ODB2_YHF(J,I,L) =  DENOM
          ODC2_YHF(J,I,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)
#endif

        ODD_YHF  (:,I,L) = (ODC_YHF (:,I,L)+2.)/(ODB_YHF (:,I,L)+1.)
        ODD2_YHF (:,I,L) = (ODC2_YHF(:,I,L)+2.)/(ODB2_YHF(:,I,L)+1.)

      END DO
    END IF

!______________________________________________________ LEFT TRANSPARENT

    IF ( X_MINF <= 1 ) THEN
      DO  K = Y_MINF, Y_MAXF
        JM1 = JMF ( X_MINF, K, L )
        VELP = SQRT( L2MX(JM1) *HF*DENUF(JM1)/DT )
        VELS = SQRT( MUXY(JM1) *HF*DENVF(JM1)/DT )
        VEL  = VELP

        ODA_XLF  (:,K,L) = 1._PP
        ODA2_XLF (:,K,L) = 1._PP

        ODB_XLF  (:,K,L) = 1._PP
        ODB2_XLF (:,K,L) = 1._PP

        ODC_XLF  (:,K,L) = 0._PP
        ODC2_XLF (:,K,L) = 0._PP

#if(SIMPLE_PML)
        DO J = 1, TPMLF
          OMGPDT = OMG*VEL* REAL(J)     **POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA_XLF (J,K,L) =  2._PP                         / DENOM
          ODB_XLF (J,K,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC_XLF (J,K,L) = -2._PP*     OMGPDT             / DENOM

          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA2_XLF(J,K,L) =  2._PP                         / DENOM
          ODB2_XLF(J,K,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC2_XLF(J,K,L) = -2._PP*     OMGPDT             / DENOM

        END DO

        J = TPMLF+1
          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA2_XLF(J,K,L) =  2._PP                         / DENOM
          ODB2_XLF(J,K,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC2_XLF(J,K,L) = -2._PP*     OMGPDT             / DENOM

#else
        DO J = 1, TPMLF
          OMGPDT = OMG*VEL* REAL(J)     **POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA_XLF (J,K,L) = (DENOM+1.)/2./AP
          ODB_XLF (J,K,L) =  DENOM
          ODC_XLF (J,K,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA2_XLF(J,K,L) = (DENOM+1.)/2./AP
          ODB2_XLF(J,K,L) =  DENOM
          ODC2_XLF(J,K,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

        END DO

        J = TPMLF+1
          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA2_XLF(J,K,L) = (DENOM+1.)/2./AP
          ODB2_XLF(J,K,L) =  DENOM
          ODC2_XLF(J,K,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)
#endif

        ODD_XLF  (:,K,L) = (ODC_XLF (:,K,L)+2.)/(ODB_XLF (:,K,L)+1.)
        ODD2_XLF (:,K,L) = (ODC2_XLF(:,K,L)+2.)/(ODB2_XLF(:,K,L)+1.)

     END DO
    END IF

!_____________________________________________________ RIGHT TRANSPARENT

    IF ( X_MAXF >= MXF ) THEN
      DO  K = Y_MINF, Y_MAXF
        JM1 = JMF (X_MAXF ,K,L)
        VELP = SQRT( L2MX(JM1) *HF*DENUF(JM1)/DT )
        VELS = SQRT( MUXY(JM1) *HF*DENVF(JM1)/DT )
        VEL  = VELP

        ODA_XHF  (:,K,L) = 1._PP
        ODA2_XHF (:,K,L) = 1._PP

        ODB_XHF  (:,K,L) = 1._PP
        ODB2_XHF (:,K,L) = 1._PP

        ODC_XHF  (:,K,L) = 0._PP
        ODC2_XHF (:,K,L) = 0._PP

#if(SIMPLE_PML)
       DO J = 1, TPMLF
          OMGPDT = OMG*VEL* REAL(J)     **POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA_XHF (J,K,L) =  2._PP                         / DENOM
          ODB_XHF (J,K,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC_XHF (J,K,L) = -2._PP*     OMGPDT             / DENOM

          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA2_XHF(J,K,L) =  2._PP                         / DENOM
          ODB2_XHF(J,K,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC2_XHF(J,K,L) = -2._PP*     OMGPDT             / DENOM

        END DO

        J = TPMLF+1
          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA2_XHF(J,K,L) =  2._PP                         / DENOM
          ODB2_XHF(J,K,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC2_XHF(J,K,L) = -2._PP*     OMGPDT             / DENOM

#else
        DO J = 1, TPMLF
          OMGPDT = OMG*VEL* REAL(J)     **POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA_XHF (J,K,L) = (DENOM+1.)/2./AP
          ODB_XHF (J,K,L) =  DENOM
          ODC_XHF (J,K,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA2_XHF(J,K,L) = (DENOM+1.)/2./AP
          ODB2_XHF(J,K,L) =  DENOM
          ODC2_XHF(J,K,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

        END DO

        J = TPMLF+1
          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA2_XHF(J,K,L) = (DENOM+1.)/2./AP
          ODB2_XHF(J,K,L) =  DENOM
          ODC2_XHF(J,K,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)
#endif

        ODD_XHF  (:,K,L) = (ODC_XHF (:,K,L)+2.)/(ODB_XHF (:,K,L)+1.)
        ODD2_XHF (:,K,L) = (ODC2_XHF(:,K,L)+2.)/(ODB2_XHF(:,K,L)+1.)

      END DO
    END IF

  END DO

END IF


END SUBROUTINE PPMLF

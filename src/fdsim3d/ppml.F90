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

SUBROUTINE PPML

  USE PRECISION   , ONLY:                                              &
                          PP
  USE CONTROL_DATA, ONLY:                                              &
                          MX, MY, MZ,  DT,                             &
                          TPML, R, TAU, POW, SF, WC
  USE GRID_MEDIUM , ONLY:                                              &
                          JM, L2MX, L2MY, L2MZ,                        &
                          MUXY, MUXZ,                                  &
                          DENU, DENV, DENW, H

  USE PML         , ONLY:                                              &
            ODA_XL ,  ODB_XL ,  ODC_XL ,  ODA_XH ,  ODB_XH ,  ODC_XH , &
            ODD_XL ,                      ODD_XH ,                     &
            ODA_YL ,  ODB_YL ,  ODC_YL ,  ODA_YH ,  ODB_YH ,  ODC_YH , &
            ODD_YL ,                      ODD_YH ,                     &
                                          ODA_ZH ,  ODB_ZH ,  ODC_ZH , &
                                          ODD_ZH ,                     &
           ODA2_XL , ODB2_XL , ODC2_XL , ODA2_XH , ODB2_XH , ODC2_XH , &
           ODD2_XL ,                     ODD2_XH ,                     &
           ODA2_YL , ODB2_YL , ODC2_YL , ODA2_YH , ODB2_YH , ODC2_YH , &
           ODD2_YL ,                     ODD2_YH ,                     &
                                         ODA2_ZH , ODB2_ZH , ODC2_ZH , &
                                         ODD2_ZH ,                     &
                                        PZXL , PZXH ,                  &
                                        PZYL , PZYH ,                  &
                                        PXXL , PXXH ,                  &
                                        PYYL , PYYH ,                  &
                                        PYXL , PYXH ,                  &
                                        PXYL , PXYH ,                  &
                                               PXZH ,                  &
                                               PYZH ,                  &
                                               PZZH ,                  &
                                          RXXL , RXXH ,                &
                                          RXYL , RXYH ,                &
                                          RYYL , RYYH ,                &
                                          RYXL , RYXH ,                &
                                          RZXL , RZXH ,                &
                                          RZYL , RZYH ,                &
                                  RXZH ,                               &
                                  RYZH ,                               &
                                  RZZH

  USE AUXIL       , ONLY:                                              &
                          PMLL, PMLXH, PMLYH, PMLZH,                   &
                           X_MIN_D, X_MAX_D, X_MIN, X_MAX,             &
                           Y_MIN_D, Y_MAX_D, Y_MIN, Y_MAX,             &
                           Z_MIN_D, Z_MAX_D, Z_MIN, Z_MAX

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER        :: I, L, K, J, JM1

  REAL (KIND=PP) :: VEL, VELS, VELP, OMG, OMGPDT, AP, DENOM

!-----------------------------------------------------------------------

IF ( Z_MIN_D <= Z_MAX_D ) THEN
  ALLOCATE (ODA_XL   (0:TPML   , Y_MIN_D :Y_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODB_XL   (0:TPML   , Y_MIN_D :Y_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODC_XL   (0:TPML   , Y_MIN_D :Y_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODD_XL   (0:TPML   , Y_MIN_D :Y_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODA_XH   (0:TPML   , Y_MIN_D :Y_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODB_XH   (0:TPML   , Y_MIN_D :Y_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODC_XH   (0:TPML   , Y_MIN_D :Y_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODD_XH   (0:TPML   , Y_MIN_D :Y_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODA_YL   (0:TPML   , X_MIN_D :X_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODB_YL   (0:TPML   , X_MIN_D :X_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODC_YL   (0:TPML   , X_MIN_D :X_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODD_YL   (0:TPML   , X_MIN_D :X_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODA_YH   (0:TPML   , X_MIN_D :X_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODB_YH   (0:TPML   , X_MIN_D :X_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODC_YH   (0:TPML   , X_MIN_D :X_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODD_YH   (0:TPML   , X_MIN_D :X_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODA_ZH   (0:TPML   , X_MIN_D :X_MAX_D , Y_MIN_D :Y_MAX_D ),&
            ODB_ZH   (0:TPML   , X_MIN_D :X_MAX_D , Y_MIN_D :Y_MAX_D ),&
            ODC_ZH   (0:TPML   , X_MIN_D :X_MAX_D , Y_MIN_D :Y_MAX_D ),&
            ODD_ZH   (0:TPML   , X_MIN_D :X_MAX_D , Y_MIN_D :Y_MAX_D ),&
            ODA2_XL  (0:TPML +1, Y_MIN_D :Y_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODB2_XL  (0:TPML +1, Y_MIN_D :Y_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODC2_XL  (0:TPML +1, Y_MIN_D :Y_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODD2_XL  (0:TPML +1, Y_MIN_D :Y_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODA2_XH  (0:TPML +1, Y_MIN_D :Y_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODB2_XH  (0:TPML +1, Y_MIN_D :Y_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODC2_XH  (0:TPML +1, Y_MIN_D :Y_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODD2_XH  (0:TPML +1, Y_MIN_D :Y_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODA2_YL  (0:TPML +1, X_MIN_D :X_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODB2_YL  (0:TPML +1, X_MIN_D :X_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODC2_YL  (0:TPML +1, X_MIN_D :X_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODD2_YL  (0:TPML +1, X_MIN_D :X_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODA2_YH  (0:TPML +1, X_MIN_D :X_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODB2_YH  (0:TPML +1, X_MIN_D :X_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODC2_YH  (0:TPML +1, X_MIN_D :X_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODD2_YH  (0:TPML +1, X_MIN_D :X_MAX_D , Z_MIN_D :Z_MAX_D ),&
            ODA2_ZH  (0:TPML +1, X_MIN_D :X_MAX_D , Y_MIN_D :Y_MAX_D ),&
            ODB2_ZH  (0:TPML +1, X_MIN_D :X_MAX_D , Y_MIN_D :Y_MAX_D ),&
            ODC2_ZH  (0:TPML +1, X_MIN_D :X_MAX_D , Y_MIN_D :Y_MAX_D ),&
            ODD2_ZH  (0:TPML +1, X_MIN_D :X_MAX_D , Y_MIN_D :Y_MAX_D ) )

  ALLOCATE (  PZXL (PMLL   :2      , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D),&
              PZXH (MX     :PMLXH  , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D),&
              PZYL (X_MIN_D:X_MAX_D, PMLL   :2      , Z_MIN_D:Z_MAX_D),&
              PZYH (X_MIN_D:X_MAX_D, MY     :PMLYH  , Z_MIN_D:Z_MAX_D),&
              PXXL (PMLL   :1      , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D),&
              PXXH (MX     :PMLXH  , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D),&
              PYYL (X_MIN_D:X_MAX_D, PMLL   :1      , Z_MIN_D:Z_MAX_D),&
              PYYH (X_MIN_D:X_MAX_D, MY     :PMLYH  , Z_MIN_D:Z_MAX_D),&
              PYXL (PMLL   :2      , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D),&
              PYXH (MX     :PMLXH  , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D),&
              PXYL (X_MIN_D:X_MAX_D, PMLL   :2      , Z_MIN_D:Z_MAX_D),&
              PXYH (X_MIN_D:X_MAX_D, MY     :PMLYH  , Z_MIN_D:Z_MAX_D),&
              PXZH (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D,       0:TPML+1 ),&
              PYZH (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D,       0:TPML+1 ),&
              PZZH (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D,       0:TPML   ) )

  ALLOCATE (   RXXL(PMLL   :3      , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D),&
               RXXH(MX-1   :PMLXH  , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D),&
               RXYL(X_MIN_D:X_MAX_D, PMLL   :2      , Z_MIN_D:Z_MAX_D),&
               RXYH(X_MIN_D:X_MAX_D, MY-1   :PMLYH  , Z_MIN_D:Z_MAX_D),&
               RYYL(X_MIN_D:X_MAX_D, PMLL   :3      , Z_MIN_D:Z_MAX_D),&
               RYYH(X_MIN_D:X_MAX_D, MY-1   :PMLYH  , Z_MIN_D:Z_MAX_D),&
               RYXL(PMLL   :2      , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D),&
               RYXH(MX-1   :PMLXH  , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D),&
               RZXL(PMLL   :2      , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D),&
               RZXH(MX-1   :PMLXH  , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D),&
               RZYL(X_MIN_D:X_MAX_D, PMLL   :2      , Z_MIN_D:Z_MAX_D),&
               RZYH(X_MIN_D:X_MAX_D, MY-1   :PMLYH  , Z_MIN_D:Z_MAX_D),&
               RXZH(X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, 0      :TPML   ),&
               RYZH(X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, 0      :TPML   ),&
               RZZH(X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, 0      :TPML   ) )

  PZXL (PMLL   :2      , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0._PP
  PZXH (MX     :PMLXH  , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0._PP
  PZYL (X_MIN_D:X_MAX_D, PMLL   :2      , Z_MIN_D:Z_MAX_D) = 0._PP
  PZYH (X_MIN_D:X_MAX_D, MY     :PMLYH  , Z_MIN_D:Z_MAX_D) = 0._PP
  PXXL (PMLL   :1      , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0._PP
  PXXH (MX     :PMLXH  , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0._PP
  PYYL (X_MIN_D:X_MAX_D, PMLL   :1      , Z_MIN_D:Z_MAX_D) = 0._PP
  PYYH (X_MIN_D:X_MAX_D, MY     :PMLYH  , Z_MIN_D:Z_MAX_D) = 0._PP
  PYXL (PMLL   :2      , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0._PP
  PYXH (MX     :PMLXH  , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0._PP
  PXYL (X_MIN_D:X_MAX_D, PMLL   :2      , Z_MIN_D:Z_MAX_D) = 0._PP
  PXYH (X_MIN_D:X_MAX_D, MY     :PMLYH  , Z_MIN_D:Z_MAX_D) = 0._PP
  PXZH (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, 0      :TPML+1 ) = 0._PP
  PYZH (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, 0      :TPML+1 ) = 0._PP
  PZZH (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, 0      :TPML   ) = 0._PP

  RXXL (PMLL   :3      , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0._PP
  RXXH (MX-1   :PMLXH  , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0._PP
  RXYL (X_MIN_D:X_MAX_D, PMLL   :2      , Z_MIN_D:Z_MAX_D) = 0._PP
  RXYH (X_MIN_D:X_MAX_D, MY-1   :PMLYH  , Z_MIN_D:Z_MAX_D) = 0._PP
  RYYL (X_MIN_D:X_MAX_D, PMLL   :3      , Z_MIN_D:Z_MAX_D) = 0._PP
  RYYH (X_MIN_D:X_MAX_D, MY-1   :PMLYH  , Z_MIN_D:Z_MAX_D) = 0._PP
  RYXL (PMLL   :2      , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0._PP
  RYXH (MX-1   :PMLXH  , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0._PP
  RZXL (PMLL   :2      , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0._PP
  RZXH (MX-1   :PMLXH  , Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) = 0._PP
  RZYL (X_MIN_D:X_MAX_D, PMLL   :2      , Z_MIN_D:Z_MAX_D) = 0._PP
  RZYH (X_MIN_D:X_MAX_D, MY-1   :PMLYH  , Z_MIN_D:Z_MAX_D) = 0._PP
  RXZH (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, 0      :TPML   ) = 0._PP
  RYZH (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, 0      :TPML   ) = 0._PP
  RZZH (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, 0      :TPML   ) = 0._PP

!                                                   HETEROGENEOUS MEDIUM
  IF ( TPML == 0 ) THEN
    OMG = 0._PP
  ELSE
    OMG = -TAU*LOG(R)*DT/REAL(TPML)/(REAL(TPML)**POW)/H
  END IF

  DO  L = Z_MIN, Z_MAX
!______________________________________________________ REAR TRANSPARENT

    IF ( Y_MIN <= 1 ) THEN
      DO  I = X_MIN, X_MAX
        JM1 = JM ( I  , Y_MIN, L )

        VELS = SQRT( MUXY(JM1) *H*DENU(JM1)/DT )
        VELP = SQRT( L2MY(JM1) *H*DENV(JM1)/DT )
        VEL  = VELP

        ODA_YL   (:,I,L) = 1._PP
        ODA2_YL  (:,I,L) = 1._PP

        ODB_YL   (:,I,L) = 1._PP
        ODB2_YL  (:,I,L) = 1._PP

        ODC_YL   (:,I,L) = 0._PP
        ODC2_YL  (:,I,L) = 0._PP

#if(SIMPLE_PML)
        DO J = 1, TPML
          OMGPDT = OMG*VEL* REAL(J)     **POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA_YL  (J,I,L) =  2._PP                         / DENOM
          ODB_YL  (J,I,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC_YL  (J,I,L) = -2._PP*     OMGPDT             / DENOM

          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA2_YL (J,I,L) =  2._PP                         / DENOM
          ODB2_YL (J,I,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC2_YL (J,I,L) = -2._PP*     OMGPDT             / DENOM

        END DO

        J = TPML +1
          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA2_YL (J,I,L) =  2._PP                         / DENOM
          ODB2_YL (J,I,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC2_YL (J,I,L) = -2._PP*     OMGPDT             / DENOM
#else
        DO J = 1, TPML
          OMGPDT = OMG*VEL* REAL(J)     **POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA_YL  (J,I,L) = (DENOM+1.)/2./AP
          ODB_YL  (J,I,L) =  DENOM
          ODC_YL  (J,I,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA2_YL (J,I,L) = (DENOM+1.)/2./AP
          ODB2_YL (J,I,L) =  DENOM
          ODC2_YL (J,I,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

        END DO

        J = TPML +1
          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA2_YL (J,I,L) = (DENOM+1.)/2./AP
          ODB2_YL (J,I,L) =  DENOM
          ODC2_YL (J,I,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)
#endif

        ODD_YL   (:,I,L) = (ODC_YL  (:,I,L)+2.)/(ODB_YL  (:,I,L)+1.)
        ODD2_YL  (:,I,L) = (ODC2_YL (:,I,L)+2.)/(ODB2_YL (:,I,L)+1.)

      END DO
    END IF

!_____________________________________________________ FRONT TRANSPARENT

    IF ( Y_MAX >= MY ) THEN
      DO  I = X_MIN, X_MAX
        JM1 = JM ( I  , Y_MAX, L )

        VELS = SQRT( MUXY(JM1) *H*DENU(JM1)/DT )
        VELP = SQRT( L2MY(JM1) *H*DENV(JM1)/DT )
        VEL  = VELP

        ODA_YH   (:,I,L) = 1._PP
        ODA2_YH  (:,I,L) = 1._PP

        ODB_YH   (:,I,L) = 1._PP
        ODB2_YH  (:,I,L) = 1._PP

        ODC_YH   (:,I,L) = 0._PP
        ODC2_YH  (:,I,L) = 0._PP

#if(SIMPLE_PML)
        DO J = 1, TPML
          OMGPDT = OMG*VEL* REAL(J)     **POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA_YH  (J,I,L) =  2._PP                         / DENOM
          ODB_YH  (J,I,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC_YH  (J,I,L) = -2._PP*     OMGPDT             / DENOM

          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA2_YH (J,I,L) =  2._PP                         / DENOM
          ODB2_YH (J,I,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC2_YH (J,I,L) = -2._PP*     OMGPDT             / DENOM

        END DO

        J = TPML +1
          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA2_YH (J,I,L) =  2._PP                         / DENOM
          ODB2_YH (J,I,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC2_YH (J,I,L) = -2._PP*     OMGPDT             / DENOM
#else
        DO J = 1, TPML
          OMGPDT = OMG*VEL* REAL(J)     **POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA_YH  (J,I,L) = (DENOM+1.)/2./AP
          ODB_YH  (J,I,L) =  DENOM
          ODC_YH  (J,I,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA2_YH (J,I,L) = (DENOM+1.)/2./AP
          ODB2_YH (J,I,L) =  DENOM
          ODC2_YH (J,I,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

        END DO

        J = TPML +1
          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA2_YH (J,I,L) = (DENOM+1.)/2./AP
          ODB2_YH (J,I,L) =  DENOM
          ODC2_YH (J,I,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)
#endif

        ODD_YH   (:,I,L) = (ODC_YH  (:,I,L)+2.)/(ODB_YH  (:,I,L)+1.)
        ODD2_YH  (:,I,L) = (ODC2_YH (:,I,L)+2.)/(ODB2_YH (:,I,L)+1.)

      END DO
    END IF

!______________________________________________________ LEFT TRANSPARENT

    IF ( X_MIN <= 1 ) THEN
      DO  K = Y_MIN, Y_MAX
        JM1 = JM ( X_MIN, K, L )
        VELP = SQRT( L2MX(JM1) *H*DENU(JM1)/DT )
        VELS = SQRT( MUXY(JM1) *H*DENV(JM1)/DT )
        VEL  = VELP

        ODA_XL   (:,K,L) = 1._PP
        ODA2_XL  (:,K,L) = 1._PP

        ODB_XL   (:,K,L) = 1._PP
        ODB2_XL  (:,K,L) = 1._PP

        ODC_XL   (:,K,L) = 0._PP
        ODC2_XL  (:,K,L) = 0._PP

#if(SIMPLE_PML)
        DO J = 1, TPML
          OMGPDT = OMG*VEL* REAL(J)     **POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA_XL  (J,K,L) =  2._PP                         / DENOM
          ODB_XL  (J,K,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC_XL  (J,K,L) = -2._PP*     OMGPDT             / DENOM

          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA2_XL (J,K,L) =  2._PP                         / DENOM
          ODB2_XL (J,K,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC2_XL (J,K,L) = -2._PP*     OMGPDT             / DENOM

        END DO

        J = TPML +1
          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA2_XL (J,K,L) =  2._PP                         / DENOM
          ODB2_XL (J,K,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC2_XL (J,K,L) = -2._PP*     OMGPDT             / DENOM

#else
        DO J = 1, TPML
          OMGPDT = OMG*VEL* REAL(J)     **POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA_XL  (J,K,L) = (DENOM+1.)/2./AP
          ODB_XL  (J,K,L) =  DENOM
          ODC_XL  (J,K,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA2_XL (J,K,L) = (DENOM+1.)/2./AP
          ODB2_XL (J,K,L) =  DENOM
          ODC2_XL (J,K,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

        END DO

        J = TPML +1
          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA2_XL (J,K,L) = (DENOM+1.)/2./AP
          ODB2_XL (J,K,L) =  DENOM
          ODC2_XL (J,K,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)
#endif

        ODD_XL   (:,K,L) = (ODC_XL  (:,K,L)+2.)/(ODB_XL  (:,K,L)+1.)
        ODD2_XL  (:,K,L) = (ODC2_XL (:,K,L)+2.)/(ODB2_XL (:,K,L)+1.)

      END DO
    END IF

!_____________________________________________________ RIGHT TRANSPARENT

    IF ( X_MAX >= MX ) THEN
      DO  K = Y_MIN, Y_MAX
        JM1 = JM(X_MAX ,K,L)
        VELP = SQRT( L2MX(JM1) *H*DENU(JM1)/DT )
        VELS = SQRT( MUXY(JM1) *H*DENV(JM1)/DT )
        VEL  = VELP

        ODA_XH   (:,K,L) = 1._PP
        ODA2_XH  (:,K,L) = 1._PP

        ODB_XH   (:,K,L) = 1._PP
        ODB2_XH  (:,K,L) = 1._PP

        ODC_XH   (:,K,L) = 0._PP
        ODC2_XH  (:,K,L) = 0._PP

#if(SIMPLE_PML)
        DO J = 1, TPML
          OMGPDT = OMG*VEL* REAL(J)     **POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA_XH  (J,K,L) =  2._PP                         / DENOM
          ODB_XH  (J,K,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC_XH  (J,K,L) = -2._PP*     OMGPDT             / DENOM

          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA2_XH (J,K,L) =  2._PP                         / DENOM
          ODB2_XH (J,K,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC2_XH (J,K,L) = -2._PP*     OMGPDT             / DENOM

        END DO

        J = TPML +1
          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
          ODA2_XH (J,K,L) =  2._PP                         / DENOM
          ODB2_XH (J,K,L) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
          ODC2_XH (J,K,L) = -2._PP*     OMGPDT             / DENOM

#else
        DO J = 1, TPML
          OMGPDT = OMG*VEL* REAL(J)     **POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA_XH  (J,K,L) = (DENOM+1.)/2./AP
          ODB_XH  (J,K,L) =  DENOM
          ODC_XH  (J,K,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA2_XH (J,K,L) = (DENOM+1.)/2./AP
          ODB2_XH (J,K,L) =  DENOM
          ODC2_XH (J,K,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

        END DO

        J = TPML +1
          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA2_XH (J,K,L) = (DENOM+1.)/2./AP
          ODB2_XH (J,K,L) =  DENOM
          ODC2_XH (J,K,L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)
#endif

        ODD_XH   (:,K,L) = (ODC_XH  (:,K,L)+2.)/(ODB_XH  (:,K,L)+1.)
        ODD2_XH  (:,K,L) = (ODC2_XH (:,K,L)+2.)/(ODB2_XH (:,K,L)+1.)

      END DO
    END IF

  END DO

!____________________________________________________ BOTTOM TRANSPARENT

  IF ( Z_MAX >= MZ ) THEN
    JM1 = JM ( X_MIN, Y_MIN, Z_MAX )
    VELS = SQRT( MUXZ(JM1) *H*DENU(JM1)/DT )
    VELP = SQRT( L2MZ(JM1) *H*DENW(JM1)/DT )
    VEL  = VELP

    ODA_ZH   (:,:,:) = 1._PP
    ODA2_ZH  (:,:,:) = 1._PP

    ODB_ZH   (:,:,:) = 1._PP
    ODB2_ZH  (:,:,:) = 1._PP

    ODC_ZH   (:,:,:) = 0._PP
    ODC2_ZH  (:,:,:) = 0._PP

#if(SIMPLE_PML)
    DO J = 1, TPML
      OMGPDT = OMG*VEL* REAL(J)     **POW
      AP     = 1.+SF*OMGPDT*0.075/0.91
      DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
      ODA_ZH  (J,:,:) =  2._PP                         / DENOM
      ODB_ZH  (J,:,:) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
      ODC_ZH  (J,:,:) = -2._PP*     OMGPDT             / DENOM

      OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
      AP     = 1.+SF*OMGPDT*0.075/0.91
      DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
      ODA2_ZH (J,:,:) =  2._PP                         / DENOM
      ODB2_ZH (J,:,:) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
      ODC2_ZH (J,:,:) = -2._PP*     OMGPDT             / DENOM

    END DO

    J = TPML +1
      OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
      AP     = 1.+SF*OMGPDT*0.075/0.91
      DENOM  = 2._PP*AP + OMGPDT + WC*AP*DT
      ODA2_ZH (J,:,:) =  2._PP                         / DENOM
      ODB2_ZH (J,:,:) = (2._PP*AP - OMGPDT - WC*AP*DT) / DENOM
      ODC2_ZH (J,:,:) = -2._PP*     OMGPDT             / DENOM

#else
    DO J = 1, TPML
      OMGPDT = OMG*VEL* REAL(J)     **POW
      AP     = 1.+SF*OMGPDT*0.075/0.91
      DENOM = EXP(-OMGPDT/AP-WC*DT)
      ODA_ZH  (J,:,:) = (DENOM+1.)/2./AP
      ODB_ZH  (J,:,:) =  DENOM
      ODC_ZH  (J,:,:) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

      OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
      AP     = 1.+SF*OMGPDT*0.075/0.91
      DENOM = EXP(-OMGPDT/AP-WC*DT)
      ODA2_ZH (J,:,:) = (DENOM+1.)/2./AP
      ODB2_ZH (J,:,:) =  DENOM
      ODC2_ZH (J,:,:) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

    END DO

    J = TPML +1
      OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
      AP     = 1.+SF*OMGPDT*0.075/0.91
      DENOM  = EXP(-OMGPDT/AP-WC*DT)
      ODA2_ZH (J,:,:) = (DENOM+1.)/2./AP
      ODB2_ZH (J,:,:) =  DENOM
      ODC2_ZH (J,:,:) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)
#endif

    ODD_ZH   (:,:,:) = (ODC_ZH  (:,:,:)+2.)/(ODB_ZH  (:,:,:)+1.)
    ODD2_ZH  (:,:,:) = (ODC2_ZH (:,:,:)+2.)/(ODB2_ZH (:,:,:)+1.)

  END IF

END IF


END SUBROUTINE PPML

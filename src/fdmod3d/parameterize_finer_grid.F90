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
!***********************************************************************
! SUBROUTINE PARAMETERIZES FINER GRID
!***********************************************************************
SUBROUTINE PARAMETERIZE_FINER_GRID

  USE NRTYPE   , ONLY: WP, SP
  USE MOD_CONST, ONLY: N_FREQ
  USE MOD_ATTEN, ONLY: Q, Y, QYLSM, FRJ
  USE MOD_MODEL, ONLY: XBMIN, YBMIN, HF, HF2, ZHF, MXF, MYF,           &
                       KEY_Q, NK, MAXDEPTH
  USE MOD_INTEGRATION     , ONLY: PTS, PTS2,                           &
                                  RHO_FINE, MU_FINE, LA_FINE,          &
                                  YMU_FINE, YLA_FINE
  USE MOD_PARAMETERIZATION, ONLY: LSTART, LSTOP, N_PLANESF, FREF,      &
                                  CELL, CELLY, FINER_GRID,             &
                                  RHO_HOM, MU_HOM, LA_HOM,             &
                                  YMU_HOM, YL2M_HOM,                   &
                                  JMH3, N_CELLS
  USE MOD_EFFMATPAR       , ONLY: MUXZ, MUYZ, MUXY,                    &
                                  LAXZ, LAYZ, LAXY,                    &
                                  L2MX, L2MY, L2MZ
  USE MOD_EFFMATPAR_Y     , ONLY: MUXZ_Y, MUYZ_Y, MUXY_Y,              &
                                  LAXZ_Y, LAYZ_Y, LAXY_Y,              &
                                  L2MX_Y, L2MY_Y, L2MZ_Y
  USE MOD_INTERFACES      , ONLY: IS_IN_1D_MEDIUM,                     &
                                  PREPARE_FINE_PAR, PREPARE_FINE_PAR_Y,&
                                  INT3D, GET_FD_CELL_INDEX
#if(USE_MPI)
  USE MPI_SHARED
#endif

!----------------------------------------------------------------------

  IMPLICIT NONE

  LOGICAL  :: FIRST1D
  INTEGER  :: I, K, L, J, NKT, II, KI, K1, JM
  REAL(WP) :: XSTART, YSTART, XA,YA,ZA,                                &
              R1, R2, THE11, THE12, THE21, THE22, POFRJ, POFRJ2P1,     &
              MU, L2M

  REAL(SP) :: RHO_1D_UV,RHO_1D_W,                                      &
               MU_1D_XY, MU_1D_XZYZ,                                   &
               LA_1D_XY, LA_1D_XZYZ,                                   &
               L2M_1D_X_Y, L2M_1D_Z
  REAL(SP), DIMENSION(N_FREQ) :: YMU_1D_XZYZ, YMU_1D_XY,               &
                                 YLA_1D_XZYZ, YLA_1D_XY,               &
                                 YL2M_1D_X_Y, YL2M_1D_Z

!----------------------------------------------------------------------

  YSTART = YBMIN-HF2
  XSTART = XBMIN-HF2
  LSTOP  = ABS(ZHF/HF) + 1

  DO L = LSTART, LSTOP

    ZA = - REAL(L-1,WP) * HF
    FIRST1D = .TRUE.

      IF ( ZA < (MAXDEPTH-HF2) ) THEN

        CALL   PREPARE_FINE_PAR   (1, 1, 1, L, HF, 0.)

        IF (KEY_Q) THEN
          CALL PREPARE_FINE_PAR_Y (1, 1, 1, L, HF, 0.)
        END IF

        II = PTS
        KI = PTS
                                                          !RHO U
        CELL(10) = INT3D (II     , KI-PTS2, PTS , RHO_FINE )
                                                          !RHO V
        CELL(11) = INT3D (II-PTS2, KI     , PTS , RHO_FINE )
                                                          !RHO W
        CELL(12) = INT3D (II     , KI     , PTS2, RHO_FINE )
                                                          !MUXZ
        CELL(1) = MUXZ (II     , KI-PTS2, PTS2, MU_FINE )
                                                          !MUYZ
        CELL(2) = MUYZ (II-PTS2, KI     , PTS2, MU_FINE )
                                                          !MUXY
        CELL(3) = MUXY (II-PTS2, KI-PTS2, PTS , MU_FINE )
                                                          !LAXZ
        CELL(4) = LAXZ (II     , KI     , PTS , LA_FINE, MU_FINE )
                                                          !LAYZ
        CELL(5) = LAYZ (II     , KI     , PTS , LA_FINE, MU_FINE )
                                                          !LAXY
        CELL(6) = LAXY (II     , KI     , PTS , LA_FINE, MU_FINE )
                                                          !L2MX
        CELL(7) = L2MX (II     , KI     , PTS , LA_FINE, MU_FINE )
                                                          !L2MY
        CELL(8) = L2MY (II     , KI     , PTS , LA_FINE, MU_FINE )
                                                          !L2MZ
        CELL(9) = L2MZ (II     , KI     , PTS , LA_FINE, MU_FINE )

        IF (KEY_Q) THEN
          Q = MUXZ_Y (II     , KI-PTS2, PTS2, YMU_FINE )
          CALL QYLSM(Q,Y)
          CELLY(1,:) = REAL(Y,SP)

          Q = MUYZ_Y (II-PTS2, KI     , PTS2, YMU_FINE )
          CALL QYLSM(Q,Y)
          CELLY(2,:) = REAL(Y,SP)

          Q = MUXY_Y (II-PTS2, KI-PTS2, PTS , YMU_FINE )
          CALL QYLSM(Q,Y)
          CELLY(3,:) = REAL(Y,SP)

          Q = LAXZ_Y (II     , KI     , PTS , YLA_FINE, YMU_FINE )
          CALL QYLSM(Q,Y)
          CELLY(4,:) = REAL(Y,SP)

          Q = LAYZ_Y (II     , KI     , PTS , YLA_FINE, YMU_FINE )
          CALL QYLSM(Q,Y)
          CELLY(5,:) = REAL(Y,SP)

          Q = LAXY_Y (II     , KI     , PTS , YLA_FINE, YMU_FINE )
          CALL QYLSM(Q,Y)
          CELLY(6,:) = REAL(Y,SP)

          Q = L2MX_Y (II     , KI     , PTS , YLA_FINE, YMU_FINE )
          CALL QYLSM(Q,Y)
          CELLY(7,:) = REAL(Y,SP)

          Q = L2MY_Y (II     , KI     , PTS , YLA_FINE, YMU_FINE )
          CALL QYLSM(Q,Y)
          CELLY(8,:) = REAL(Y,SP)

          Q = L2MZ_Y (II     , KI     , PTS , YLA_FINE, YMU_FINE )
          CALL QYLSM(Q,Y)
          CELLY(9,:) = REAL(Y,SP)
        END IF

        JM = GET_FD_CELL_INDEX()
        JMH3 (:,:,L-1) = JM

      ELSE

        DO K1 = 1, MYF, NK
          NKT = MIN(NK, MYF-K1+1)
          CALL   PREPARE_FINE_PAR   (K1, MXF, NKT, L, HF, 0.)

          IF (KEY_Q) THEN
            CALL PREPARE_FINE_PAR_Y (K1, MXF, NKT, L, HF, 0.)
          END IF

          DO K = K1, K1+NKT-1
            YA = YSTART + REAL(K-1,WP)*HF
            KI = (K-K1+1)*PTS
            DO I = 1, MXF
              XA = XSTART + REAL(I-1,WP)*HF
              II = I*PTS

              IF (IS_IN_1D_MEDIUM(XA,YA,ZA,FINER_GRID)) THEN
          !------------------------------------------------1D-IH BEDROCK
                IF (FIRST1D) THEN
                                                                    !RHO U
                  RHO_1D_UV = INT3D (II     , KI-PTS2, PTS , RHO_FINE )
                                                                    !RHO W
                  RHO_1D_W  = INT3D (II     , KI     , PTS2, RHO_FINE )
                                                                    !MUXZ
                  MU_1D_XZYZ = MUXZ (II     , KI-PTS2, PTS2, MU_FINE )
                                                                    !MUXY
                  MU_1D_XY =   MUXY (II-PTS2, KI-PTS2, PTS , MU_FINE )
                                                                    !LAXZ
                  LA_1D_XZYZ = LAXZ (II     , KI     , PTS , LA_FINE, MU_FINE )
                                                                    !LAXY
                  LA_1D_XY   = LAXY (II     , KI     , PTS , LA_FINE, MU_FINE )
                                                                    !L2MX
                  L2M_1D_X_Y = L2MX (II     , KI     , PTS , LA_FINE, MU_FINE )
                                                                    !L2MZ
                  L2M_1D_Z   = L2MZ (II     , KI     , PTS , LA_FINE, MU_FINE )

                  IF (KEY_Q) THEN
                    Q = MUXZ_Y (II     , KI-PTS2, PTS2, YMU_FINE )
                    CALL QYLSM(Q,Y)
                    YMU_1D_XZYZ = REAL(Y,SP)

                    Q = MUXY_Y (II-PTS2, KI-PTS2, PTS , YMU_FINE )
                    CALL QYLSM(Q,Y)
                    YMU_1D_XY = REAL(Y,SP)

                    Q = LAXZ_Y (II     , KI     , PTS , YLA_FINE, YMU_FINE )
                    CALL QYLSM(Q,Y)
                    YLA_1D_XZYZ = REAL( Y, SP)

                    Q = LAXY_Y (II     , KI     , PTS , YLA_FINE, YMU_FINE )
                    CALL QYLSM(Q,Y)
                    YLA_1D_XY = REAL(Y,SP)

                    Q = L2MX_Y (II     , KI     , PTS , YLA_FINE, YMU_FINE )
                    CALL QYLSM(Q,Y)
                    YL2M_1D_X_Y = REAL( Y, SP)

                    Q = L2MZ_Y (II     , KI     , PTS , YLA_FINE, YMU_FINE )
                    CALL QYLSM(Q,Y)
                    YL2M_1D_Z = REAL(Y,SP)

                  END IF

                  FIRST1D = .FALSE.

                END IF
                                                                  !RHO U
                CELL(10) = RHO_1D_UV
                                                                  !RHO V
                CELL(11) = RHO_1D_UV
                                                                  !RHO W
                CELL(12) = RHO_1D_W
                                                                  !MUXZ
                CELL(1) = MU_1D_XZYZ
                                                                  !MUYZ
                CELL(2) = MU_1D_XZYZ
                                                                  !MUXY
                CELL(3) = MU_1D_XY
                                                                  !LAXZ
                CELL(4) = LA_1D_XZYZ
                                                                  !LAXZ
                CELL(5) = LA_1D_XZYZ
                                                                  !LAXY
                CELL(6) = LA_1D_XY
                                                                  !L2MX
                CELL(7) = L2M_1D_X_Y
                                                                  !L2MY
                CELL(8) = L2M_1D_X_Y
                                                                  !L2MZ
                CELL(9) = L2M_1D_Z

                IF ( KEY_Q ) THEN
                  CELLY(1,:) = YMU_1D_XZYZ(:)
                  CELLY(2,:) = YMU_1D_XZYZ(:)
                  CELLY(3,:) = YMU_1D_XY  (:)
                  CELLY(4,:) = YLA_1D_XZYZ(:)
                  CELLY(5,:) = YLA_1D_XZYZ(:)
                  CELLY(6,:) = YLA_1D_XY  (:)
                  CELLY(7,:) = YL2M_1D_X_Y(:)
                  CELLY(8,:) = YL2M_1D_X_Y(:)
                  CELLY(9,:) = YL2M_1D_Z  (:)
                END IF

                JM = GET_FD_CELL_INDEX()
                JMH3 (K,I,L-1) = JM

              ELSE
          !----------------------------------------------3D-IH STRUCTURE

                                                                  !RHO U
                CELL(10) = INT3D (II     , KI-PTS2, PTS , RHO_FINE )
                                                                  !RHO V
                CELL(11) = INT3D (II-PTS2, KI     , PTS , RHO_FINE )
                                                                  !RHO W
                CELL(12) = INT3D (II     , KI     , PTS2, RHO_FINE )
                                                                  !MUXZ
                CELL(1) = MUXZ (II     , KI-PTS2, PTS2, MU_FINE )
                                                                  !MUYZ
                CELL(2) = MUYZ (II-PTS2, KI     , PTS2, MU_FINE )
                                                                  !MUXY
                CELL(3) = MUXY (II-PTS2, KI-PTS2, PTS , MU_FINE )
                                                                  !LAXZ
                CELL(4) = LAXZ (II     , KI     , PTS , LA_FINE, MU_FINE )
                                                                  !LAYZ
                CELL(5) = LAYZ (II     , KI     , PTS , LA_FINE, MU_FINE )
                                                                  !LAXY
                CELL(6) = LAXY (II     , KI     , PTS , LA_FINE, MU_FINE )
                                                                  !L2MX
                CELL(7) = L2MX (II     , KI     , PTS , LA_FINE, MU_FINE )
                                                                  !L2MY
                CELL(8) = L2MY (II     , KI     , PTS , LA_FINE, MU_FINE )
                                                                  !L2MZ
                CELL(9) = L2MZ (II     , KI     , PTS , LA_FINE, MU_FINE )

                IF (KEY_Q) THEN
                  Q = MUXZ_Y (II     , KI-PTS2, PTS2, YMU_FINE )
                  CALL QYLSM(Q,Y)
                  CELLY(1,:) = REAL(Y,SP)

                  Q = MUYZ_Y (II-PTS2, KI     , PTS2, YMU_FINE )
                  CALL QYLSM(Q,Y)
                  CELLY(2,:) = REAL(Y,SP)

                  Q = MUXY_Y (II-PTS2, KI-PTS2, PTS , YMU_FINE )
                  CALL QYLSM(Q,Y)
                  CELLY(3,:) = REAL(Y,SP)

                  Q = LAXZ_Y (II     , KI     , PTS , YLA_FINE, YMU_FINE )
                  CALL QYLSM(Q,Y)
                  CELLY(4,:) = REAL(Y,SP)

                  Q = LAYZ_Y (II     , KI     , PTS , YLA_FINE, YMU_FINE )
                  CALL QYLSM(Q,Y)
                  CELLY(5,:) = REAL(Y,SP)

                  Q = LAXY_Y (II     , KI     , PTS , YLA_FINE, YMU_FINE )
                  CALL QYLSM(Q,Y)
                  CELLY(6,:) = REAL(Y,SP)

                  Q = L2MX_Y (II     , KI     , PTS , YLA_FINE, YMU_FINE )
                  CALL QYLSM(Q,Y)
                  CELLY(7,:) = REAL(Y,SP)

                  Q = L2MY_Y (II     , KI     , PTS , YLA_FINE, YMU_FINE )
                  CALL QYLSM(Q,Y)
                  CELLY(8,:) = REAL(Y,SP)

                  Q = L2MZ_Y (II     , KI     , PTS , YLA_FINE, YMU_FINE )
                  CALL QYLSM(Q,Y)
                  CELLY(9,:) = REAL(Y,SP)
                END IF

                JM = GET_FD_CELL_INDEX()
                JMH3 (K,I,L-1) = JM

              END IF

            END DO
          END DO
        END DO

      END IF

  END DO

END SUBROUTINE PARAMETERIZE_FINER_GRID

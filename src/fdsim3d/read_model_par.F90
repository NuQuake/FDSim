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

SUBROUTINE  READ_MODEL_PAR

  USE PRECISION   , ONLY:                                              &
                          PP
  USE GRID_MEDIUM , ONLY:                                              &
                          JMNUM,                                       &
                          H,                                           &
                          DENU , DENV , DENW ,                         &
                          DENUF, DENVF, DENWF,                         &
                          L2MX , L2MY , L2MZ , LAXY , LAXZ , LAYZ ,    &
                          MUXY , MUXZ , MUYZ ,                         &
                         CL2MX ,CL2MY ,CL2MZ ,CLAXY ,CLAXZ ,CLAYZ ,    &
                         CMUXY ,CMUXZ ,CMUYZ ,                         &
                          YL2MX , YL2MY , YL2MZ ,                      &
                          YLAXY , YLAXZ , YLAYZ ,                      &
                          YMUXY , YMUXZ , YMUYZ
  USE CONTROL_DATA, ONLY:                                              &
                          DT, LPAS, KEY_EN
  USE AUXIL       , ONLY:                                              &
                          OTC, OTD, HF
#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: J, JM1, ALLOSTAT, JMNUMR

!-----------------------------------------------------------------------

    READ ( 14 ) JMNUMR

    JMNUM = ABS(JMNUMR)

  ALLOCATE ( DENU (JMNUM), DENV (JMNUM), DENW (JMNUM),                 &
             L2MX (JMNUM), L2MY (JMNUM), L2MZ (JMNUM),                 &
             LAXY (JMNUM), LAXZ (JMNUM), LAYZ (JMNUM),                 &
             MUXY (JMNUM), MUXZ (JMNUM), MUYZ (JMNUM) )

  IF ( KEY_EN ) THEN
    ALLOCATE ( CL2MX (JMNUM), CL2MY (JMNUM), CL2MZ (JMNUM),            &
               CLAXY (JMNUM), CLAXZ (JMNUM), CLAYZ (JMNUM),            &
               CMUXY (JMNUM), CMUXZ (JMNUM), CMUYZ (JMNUM) )
  END IF

  IF ( LPAS /= 0 ) THEN
    ALLOCATE ( DENUF (JMNUM), DENVF (JMNUM), DENWF (JMNUM) )
  END IF

    IF ( JMNUMR < 0 ) THEN   ! old format (DVS5.0)

      READ  ( 14 ) ( DENU(JM1), DENV(JM1), DENW(JM1),                  &
                     L2MX(JM1), LAXY(JM1),                             & !L2MX = LAM, LAXY = MU
                     MUXY(JM1), MUXZ(JM1), MUYZ(JM1), JM1 = 1, JMNUM)

      L2MY = L2MX + 2._PP * LAXY
      LAXZ = L2MX
      LAYZ = L2MX
      LAXY = L2MX
      L2MX = L2MY
      L2MZ = L2MY

    ELSE

      READ  ( 14 )                                                     &
           ( DENU (JM1), DENV (JM1), DENW (JM1),                       &
             L2MX (JM1), L2MY (JM1), L2MZ (JM1),                       &
             LAXY (JM1), LAXZ (JM1), LAYZ (JM1),                       &
             MUXY (JM1), MUXZ (JM1), MUYZ (JM1), JM1 = 1, JMNUM)

    END IF

  IF ( KEY_EN ) CALL COMPLIANCES

  IF ( LPAS /= 0 ) THEN
    DENUF = DT/HF/DENU
    DENVF = DT/HF/DENV
    DENWF = DT/HF/DENW
  END IF

  DENU  = DT/H/DENU
  DENV  = DT/H/DENV
  DENW  = DT/H/DENW

  ALLOCATE ( YL2MX (JMNUM, 4), YL2MY (JMNUM, 4), YL2MZ (JMNUM, 4),     &
             YLAXY (JMNUM, 4), YLAXZ (JMNUM, 4), YLAYZ (JMNUM, 4),     &
             YMUXY (JMNUM, 4), YMUXZ (JMNUM, 4), YMUYZ (JMNUM, 4),     &
                                                     STAT=ALLOSTAT )

    IF  ( ALLOSTAT > 0 )  THEN
      WRITE (11,*) '  ALLOCATION ERROR R_M_PAR_02, STAT= ', ALLOSTAT
      STOP
    END IF

    IF ( JMNUMR < 0 ) THEN

      READ  ( 15 ) ( YL2MX(JM1,:), YLAXY(JM1,:),                       & !YL2MX = YLAM, YLAXY = YMU
                     YMUXY(JM1,:), YMUXZ(JM1,:), YMUYZ(JM1,:),         &
                                                       JM1 = 1, JMNUM)

      DO J = 1, 4
        YL2MY(:,J)=(YL2MX(:,J)*LAXY(:) + YLAXY(:,J)*(L2MX(:)-LAXY(:)) )&
                  /   L2MX(:)
      END DO
      YLAXY = YL2MX
      YLAXZ = YL2MX
      YLAYZ = YL2MX
      YL2MX = YL2MY
      YL2MZ = YL2MY

    ELSE

      READ  ( 15 )                                                     &
             ( YL2MX (JM1, :), YL2MY (JM1, :), YL2MZ (JM1, :),         &
               YLAXY (JM1, :), YLAXZ (JM1, :), YLAYZ (JM1, :),         &
               YMUXY (JM1, :), YMUXZ (JM1, :), YMUYZ (JM1, :),         &
                                                       JM1 = 1, JMNUM)

    END IF

  DO J = 1, 4
    YL2MX(:,J) = OTC(J) * YL2MX(:,J) * L2MX(:)
    YL2MY(:,J) = OTC(J) * YL2MY(:,J) * L2MY(:)
    YL2MZ(:,J) = OTC(J) * YL2MZ(:,J) * L2MZ(:)
    YLAXY(:,J) = OTC(J) * YLAXY(:,J) * LAXY(:)
    YLAXZ(:,J) = OTC(J) * YLAXZ(:,J) * LAXZ(:)
    YLAYZ(:,J) = OTC(J) * YLAYZ(:,J) * LAYZ(:)
    YMUXY(:,J) = OTC(J) * YMUXY(:,J) * MUXY(:)
    YMUXZ(:,J) = OTC(J) * YMUXZ(:,J) * MUXZ(:)
    YMUYZ(:,J) = OTC(J) * YMUYZ(:,J) * MUYZ(:)
  END DO

  DO J = 1, 4
    L2MX(:) = L2MX(:) - OTD(J)/OTC(J)*YL2MX(:,J)
    L2MY(:) = L2MY(:) - OTD(J)/OTC(J)*YL2MY(:,J)
    L2MZ(:) = L2MZ(:) - OTD(J)/OTC(J)*YL2MZ(:,J)
    LAXY(:) = LAXY(:) - OTD(J)/OTC(J)*YLAXY(:,J)
    LAXZ(:) = LAXZ(:) - OTD(J)/OTC(J)*YLAXZ(:,J)
    LAYZ(:) = LAYZ(:) - OTD(J)/OTC(J)*YLAYZ(:,J)
    MUXY(:) = MUXY(:) - OTD(J)/OTC(J)*YMUXY(:,J)
    MUXZ(:) = MUXZ(:) - OTD(J)/OTC(J)*YMUXZ(:,J)
    MUYZ(:) = MUYZ(:) - OTD(J)/OTC(J)*YMUYZ(:,J)
  END DO

END SUBROUTINE  READ_MODEL_PAR

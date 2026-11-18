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
!-----------------------------------------------------------------------
! subroutine reads program configuration, set some run-time constants
! and check if some parameters fulfill required conditions
!-----------------------------------------------------------------------
SUBROUTINE SETPARAMS

  USE NRTYPE   , ONLY: WP
  USE MOD_MODEL, ONLY: H, H10, KEY_Q, KEY_FILTER,          &
                       NPX, NPY, FLTR_RANGE, HF,             &
                       XBMIN, XBMAX, YBMIN, YBMAX, ZBMIN, ZHF,   &
                       XIMIN, XIMAX, YIMIN, YIMAX,                     &
                       XBMIN_G, XBMAX_G, XIMIN_G, XIMAX_G, MX_G, MXF_G,&
                       MX, MY, MZ, LPAS, NX10, NX10_G, NY10, HF2, H2,  &
                       MXF, MYF, NI, NK, KEY_QK_INF, KEY_GBM3D
  USE MOD_CONST, ONLY: N_F => N_FREQ
  USE MOD_FILES, ONLY: JMH_FILE_NAME, JMHF_FILE_NAME,            &
                             MO_FILE_NAME, Q_FILE_NAME, FSED
  USE MOD_INTEGRATION, ONLY: PTS
  USE MOD_INTERFACES , ONLY: PRINT_INFO, TEST_PARAMETERS

  USE MOD_PARAMETERIZATION, ONLY:                                      &
                             FRJMAX, FRANGE, FREF, REFINE
#if(USE_MPI)
  USE MPI_SHARED
#endif
!----------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER(LEN=17) :: JOB_NAME
  CHARACTER(LEN=20) :: MD_FILE_NAME

  INTEGER :: I, N_FREQ

  REAL(WP):: FRJMIN, X_STEP
  REAL(WP):: DT, V_MAX

  NAMELIST /CONTROL/ H, PTS, KEY_Q, NK,                                &
                     KEY_FILTER, NPX, NPY, NI,                         &
                     FLTR_RANGE, REFINE, KEY_QK_INF, KEY_GBM3D

  NAMELIST /CONTROLDATA2/ MX, MY, MZ, LPAS, H, XBMIN, YBMIN, REFINE

  NAMELIST /FILENAMES/    MO_FILE_NAME,  JMH_FILE_NAME, Q_FILE_NAME,   &
                                        JMHF_FILE_NAME

  NAMELIST /BOUNDARIES/ XBMIN, XBMAX, YBMIN, YBMAX, ZBMIN, ZHF

  NAMELIST /ATTEN/  FRJMIN, FRJMAX, FRANGE, FREF, N_FREQ
!----------------------------------------------------------------------

#if(USE_MPI)
  IF ( myid == root ) THEN
#endif
    NI         = 1
    NK         = 1
    KEY_Q      = .TRUE.
    KEY_FILTER = .FALSE.
    FRJMIN     = 0.
    KEY_QK_INF = .FALSE.
    KEY_GBM3D  = .FALSE.
    REFINE     = 1
    PTS        = 8
    FLTR_RANGE = -1.

    OPEN  (10, FILE='HF_FDSim3D', STATUS='OLD' )
    READ  (10,'(A)') JOB_NAME
    CLOSE (10)

   JMHF_FILE_NAME = TRIM(JOB_NAME)//'.JMHF'       !Fine-grid indices
    JMH_FILE_NAME = TRIM(JOB_NAME)//'.JMH'        !Coarse-grid indices
     MO_FILE_NAME = TRIM(JOB_NAME)//'.MO'         !Elastic model parameters
      Q_FILE_NAME = TRIM(JOB_NAME)//'.Q'          !Anelastic model parameters
     MD_FILE_NAME = TRIM(JOB_NAME)//'.MD'         !Model deffinitin

    OPEN(10, FILE= MD_FILE_NAME, STATUS='OLD' )
    READ(10, NML = CONTROL   )
    READ(10, NML = BOUNDARIES)
    READ(10, NML = ATTEN     )
    H10        = H/10./REFINE
    IF ( FLTR_RANGE < 0. ) FLTR_RANGE = 3.*H
    IF ( FRJMIN     > 0. ) FRANGE     = LOG10(FRJMAX/FRJMIN)

    ALLOCATE(FSED(NI))
    DO I = 1, NI
      READ(10,*)FSED(I)
    END DO
    CLOSE (10)

    XBMIN_G = XBMIN
    XBMAX_G = XBMAX
#if(USE_MPI)
  END IF

  CALL MPI_Bcast(H,         1, MPI_REAL     , root, comm, errcode)
  CALL MPI_Bcast(H10,       1, MPI_REAL     , root, comm, errcode)
  CALL MPI_Bcast(PTS,       1, MPI_INTEGER  , root, comm, errcode)
  CALL MPI_Bcast(KEY_Q,     1, MPI_LOGICAL  , root, comm, errcode)
  CALL MPI_Bcast(KEY_QK_INF,1, MPI_LOGICAL  , root, comm, errcode)
  CALL MPI_Bcast(KEY_GBM3D, 1, MPI_LOGICAL  , root, comm, errcode)
  CALL MPI_Bcast(NK,        1, MPI_INTEGER  , root, comm, errcode)
  CALL MPI_Bcast(KEY_FILTER,1, MPI_LOGICAL  , root, comm, errcode)
  CALL MPI_Bcast(NPX,       1, MPI_INTEGER  , root, comm, errcode)
  CALL MPI_Bcast(NPY,       1, MPI_INTEGER  , root, comm, errcode)
  CALL MPI_Bcast(NI,        1, MPI_INTEGER  , root, comm, errcode)
  CALL MPI_Bcast(FLTR_RANGE,1, MPI_REAL     , root, comm, errcode)
  CALL MPI_Bcast(REFINE,    1, MPI_INTEGER  , root, comm, errcode)

  CALL MPI_Bcast(YBMIN,     1, MPI_REAL     , root, comm, errcode)
  CALL MPI_Bcast(YBMAX,     1, MPI_REAL     , root, comm, errcode)
  CALL MPI_Bcast(ZBMIN,     1, MPI_REAL     , root, comm, errcode)
  CALL MPI_Bcast(ZHF,       1, MPI_REAL     , root, comm, errcode)

  CALL MPI_Bcast(FRJMAX,    1, MPI_REAL     , root, comm, errcode)
  CALL MPI_Bcast(FRANGE,    1, MPI_REAL     , root, comm, errcode)
  CALL MPI_Bcast(FREF,      1, MPI_REAL     , root, comm, errcode)

  CALL MPI_Bcast(XBMIN_G,   1, MPI_REAL     , root, comm, errcode)
  CALL MPI_Bcast(XBMAX_G,   1, MPI_REAL     , root, comm, errcode)

  X_STEP = INT((XBMAX_G-XBMIN_G)/H/numprocs)*H
  numprocsplus = (XBMAX_G-XBMIN_G)/H - INT((XBMAX_G-XBMIN_G)/H/numprocs)*numprocs

  IF ( myid == root ) THEN
    DO I = 0, numprocsplus-1
      XBMIN = XBMIN_G + I*(X_STEP+1)
      XBMAX = XBMIN   + X_STEP + 1

      IF ( I == numprocs-1 ) THEN
        XBMAX = XBMAX_G
      END IF

      XIMIN = XBMIN-H
      XIMAX = XBMAX+H

      MX_A   (I) = (XBMAX-XBMIN)/H   + 1
      MXF_A  (I) = ( MX_A(I) - 1 ) * REFINE + 1

      NX10_A (I) = (XIMAX-XIMIN)/H10 + 1
    END DO
    DO I = numprocsplus, numprocs-1
      XBMIN = XBMIN_G + I*X_STEP
      XBMAX = XBMIN   + X_STEP

      IF ( I == numprocs-1 ) THEN
        XBMAX = XBMAX_G
      END IF

      XIMIN = XBMIN-H
      XIMAX = XBMAX+H

      MX_A   (I) = (XBMAX-XBMIN)/H   + 1
      MXF_A  (I) = ( MX_A(I) - 1 ) * REFINE + 1

      NX10_A (I) = (XIMAX-XIMIN)/H10 + 1
    END DO
  END IF

  IF ( myid < numprocsplus ) THEN
    XBMIN = XBMIN_G + myid*(X_STEP+1)
  ELSE
    XBMIN = XBMIN_G + myid* X_STEP
  END IF
  XBMAX = XBMIN   + X_STEP

  IF ( myid == numprocs-1 ) THEN
    XBMAX = XBMAX_G
  END IF

#endif

!----------------------------------------- set some run-time constants
  HF = H/REFINE     ! spacing in finer grid

  XIMIN = XBMIN-H
  XIMAX = XBMAX+H
  YIMIN = YBMIN-H
  YIMAX = YBMAX+H

  XIMIN_G = XBMIN_G-H
  XIMAX_G = XBMAX_G+H

                          !-------------------------- model parameters
  MX   = (XBMAX  -XBMIN  )/H  + 1
  MX_G = (XBMAX_G-XBMIN_G)/H  + 1
  MY   = (YBMAX  -YBMIN  )/H  + 1

  LPAS = ABS(ZHF)/HF
  MZ   = ABS(ZHF-ZBMIN)/H + LPAS

  MXF  = ( MX   - 1 ) * REFINE + 1
  MXF_G= ( MX_G - 1 ) * REFINE + 1
  MYF  = ( MY   - 1 ) * REFINE + 1

                          !--number of points in aproximated interface
  NX10 = (XIMAX-XIMIN)/H10 + 1
  NY10 = (YIMAX-YIMIN)/H10 + 1

  NX10_G = (XIMAX_G-XIMIN_G)/H10 + 1

  HF2    = HF/2._WP
  H2     = H /2._WP

#if(USE_MPI)
  IF ( myid == root ) THEN
#endif
  !------------------------------------ tests of some input parameters
    CALL TEST_PARAMETERS

  !---------------------------------- print info about model on screen
    CALL PRINT_INFO

    N_FREQ = N_F

    MX   = (XBMAX_G-XBMIN_G)/H  + 1

    OPEN  ( 10, FILE=TRIM(JOB_NAME)//'.INM', STATUS = 'NEW' )
    WRITE ( 10, NML = FILENAMES )
    WRITE ( 10, NML = CONTROLDATA2 )
    WRITE ( 10, NML = ATTEN )
    CLOSE ( 10 )
    MX   = (XBMAX  -XBMIN  )/H  + 1

#if(USE_MPI)
  END IF
#endif

END SUBROUTINE SETPARAMS


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

SUBROUTINE  READ_OPEN  ( WORK_NAME )

  USE PRECISION   , ONLY:                                              &
                          PP

  USE GRID_MEDIUM , ONLY:                                              &
                          TEXT, JM, H

  USE CONTROL_DATA, ONLY:                                              &
                          SR_FILE_NAME, REC_NAME,                      &
                          KEY_TLV, KEY_SNV, KEY_EN, XBMIN, YBMIN,      &
                          MX, MY, MZ, LPAS, REFINE,                    &
                          MT1, MT2, IPAS2, MR, MSNP,                   &
                          KTLE, KTRI, KTRE, KTFR, KTBO,                &
                          DT,  REFINE, LPAS, TAU_EPS,                  &
                          OMG, WB,                                     &
                          THPPLE, THPPRI, THPPRE, THPPFR, THPPBO,      &
                          THSSLE, THSSRI, THSSRE, THSSFR, THSSBO,      &
                          IREC, KREC, LREC,                            &
                          FRJMAX, FRANGE, N_FREQ,                      &
                          TPML, R, TAU, POW, SF, WC,                   &
                          MO_FILE_NAME, JMH_FILE_NAME, Q_FILE_NAME,    &
                                       JMHF_FILE_NAME,                 &
                          STRESS_IMAGING

  USE AUXIL       , ONLY:                                              &
                          LOGSNP, MXF, MYF

#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif

!-----------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER (LEN=*), INTENT (IN) :: WORK_NAME

  CHARACTER (LEN=20) LOG_FILE_NAME,  IN_FILE_NAME, INM_FILE_NAME,      &
                     REC_FILE_NAME

  LOGICAL :: EX

  INTEGER :: IOS, I, IOLENINT, ALLOSTAT, LSNP, PROCX, PROCY, PROCZ

  REAL    :: FRJMIN, FREF, TW

  NAMELIST /FILENAMES/    MO_FILE_NAME,  JMH_FILE_NAME, Q_FILE_NAME,   &
                          JMHF_FILE_NAME

  NAMELIST /CONTROLDATA2/ MX, MY, MZ, LPAS, H, XBMIN, YBMIN, REFINE

  NAMELIST /ATTEN/        FRJMIN, FRJMAX, FRANGE, FREF, N_FREQ

  NAMELIST /CONTROLDATA/  TPML, PROCX, PROCY, PROCZ, STRESS_IMAGING

  NAMELIST /CONTROLDATA1/ TW, DT, R, TAU , POW, SF, WC, TAU_EPS

  NAMELIST /NAMES_SR_REC/ SR_FILE_NAME, REC_FILE_NAME

  NAMELIST /KEYS/         KEY_TLV, KEY_SNV, KEY_EN

  NAMELIST /NONREF/       OMG, WB, KTLE, KTRI, KTRE, KTFR, KTBO

  NAMELIST /TXT/          TEXT

  NAMELIST /SNAP/         IPAS2, MSNP

  NAMELIST /REC/          MR

#if(USE_MPI)
  IF ( myid == root ) THEN
#endif
!______________________________________________ LOG FILE __________ 11 _

    LOG_FILE_NAME = WORK_NAME//'.LOG'

    OPEN  (11, FILE = LOG_FILE_NAME, STATUS='UNKNOWN' )
    WRITE (11,*) 'LOG_FILE: JOB_NAME =', WORK_NAME

!______________________________________________ CONTROL DATA FILE _ 12 _

    IN_FILE_NAME = WORK_NAME//'.IN'

    INQUIRE ( FILE = IN_FILE_NAME, EXIST=EX, IOSTAT=IOS )
      IF  ( IOS > 0 )  THEN
        WRITE (11,*) 'ERROR HAS OCCURED DURING THE EXECUTION           &
                     &OF THE INQUIRE STATEMENT.',                      &
                     'IOSTAT (IN_FILE_NAME)=',IOS
        STOP
      END IF
      IF  ( .NOT.EX )  THEN
        WRITE (11,*) IN_FILE_NAME,' NOT FOUND'
        STOP
      END IF
    OPEN  (12, FILE = IN_FILE_NAME, STATUS='OLD' )

!_____________________________________________ CONTROL DATA FILE 2_ 10 _

    INM_FILE_NAME = WORK_NAME//'.INM'

    INQUIRE ( FILE = INM_FILE_NAME, EXIST=EX, IOSTAT=IOS )
      IF  ( IOS > 0 )  THEN
        WRITE (11,*) 'ERROR HAS OCCURED DURING THE EXECUTION           &
                     &OF THE INQUIRE STATEMENT.',                      &
                     'IOSTAT (IN_FILE_NAME)=',IOS
        STOP
      END IF
      IF  ( .NOT.EX )  THEN
        WRITE (11,*) INM_FILE_NAME,' NOT FOUND'
        STOP
      END IF
    OPEN  (10, FILE = INM_FILE_NAME, STATUS='OLD' )

!______________________________________________ READ CONTROL DATA ______

    READ  (10, NML=FILENAMES)
    WRITE (11, NML=FILENAMES)

    REFINE       = 1
    LPAS         = 0
    XBMIN        = 0.
    YBMIN        = 0.
    READ  (10, NML=CONTROLDATA2)
    WRITE (11, NML=CONTROLDATA2)

    FRJMIN  = 0.
    N_FREQ  = 4
    READ  (10, NML=ATTEN)
    WRITE (11, NML=ATTEN)
    IF ( FRJMIN >  0. ) FRANGE = LOG10(FRJMAX/FRJMIN)

    TPML         = 0
    PROCX        = 0
    PROCY        = 0
    PROCZ        = 0
    STRESS_IMAGING = .FALSE.
    READ  (12, NML=CONTROLDATA)
    WRITE (11, NML=CONTROLDATA)

    IF (STRESS_IMAGING) THEN
      IF ( LPAS > 0 ) THEN
        LPAS = LPAS + 2
      END IF
      MZ = MZ+2
    END IF

#if(USE_MPI)
    dims(1) = PROCY
    dims(2) = PROCX
    dims(3) = PROCZ
#endif

    R       = 0.001
    TAU     = 1.5
    POW     = 2
    SF      = 0.
    WC      = 0.
    TAU_EPS = 0.
    READ  (12, NML=CONTROLDATA1)
    WRITE (11, NML=CONTROLDATA1)
    MT1 = 1
    MT2 = INT(TW/DT) + 1

    READ  (12, NML=NAMES_SR_REC)
    WRITE (11, NML=NAMES_SR_REC)

    KEY_EN    = .FALSE.
    KEY_TLV   = .FALSE.
    KEY_SNV   = .FALSE.
    READ  (12, NML=KEYS)
    WRITE (11, NML=KEYS)

    IF  ( .NOT. ( KEY_TLV .OR. KEY_SNV ) ) THEN
      WRITE (11,*) ' NONE OF TL_ AND SN_ OUTPUT FILES IS REQUIRED'
      STOP
    END IF

    WB     = 0.4_PP
    OMG    = FREF
    THPPLE = 0.19757_PP
    THPPRI = 0.19757_PP
    THPPRE = 0.19757_PP
    THPPFR = 0.19757_PP
    THPPBO = 0.19757_PP
    THSSLE = 0.70511_PP
    THSSRI = 0.70511_PP
    THSSRE = 0.70511_PP
    THSSFR = 0.70511_PP
    THSSBO = 0.70511_PP

    READ  (12, NML=NONREF)
    WRITE (11, NML=NONREF)
    OMG = OMG/6.283185307179586476925286766559

    READ  (12, NML=TXT)
    WRITE (11, NML=TXT)

    IPAS2  = 1
#if(USE_MPI)
  END IF

  CALL MPI_Bcast( JMH_FILE_NAME, 20, MPI_CHARACTER, root, comm, errcode)
  CALL MPI_Bcast(JMHF_FILE_NAME, 20, MPI_CHARACTER, root, comm, errcode)
  CALL MPI_Bcast(  MO_FILE_NAME, 20, MPI_CHARACTER, root, comm, errcode)
  CALL MPI_Bcast(   Q_FILE_NAME, 20, MPI_CHARACTER, root, comm, errcode)

  CALL MPI_Bcast(TPML,      1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(REFINE,    1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(MX,        1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(MY,        1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(MZ,        1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(KEY_TLV  , 1, MPI_LOGICAL, root, comm, errcode)
  CALL MPI_Bcast(KEY_EN   , 1, MPI_LOGICAL, root, comm, errcode)
  CALL MPI_Bcast(dims,  ndims, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(STRESS_IMAGING, 1, MPI_LOGICAL, root, comm, errcode)
  CALL MPI_Bcast(TAU_EPS  , 1, MPI_REAL   , root, comm, errcode)

#endif

  MXF = (MX-1)*REFINE + 1
  MYF = (MY-1)*REFINE + 1

  ALLOCATE ( LOGSNP (0:MZ+TPML), STAT=ALLOSTAT )
    IF  ( ALLOSTAT > 0 )  THEN
      WRITE (11,*) '  ALLOCATION ERROR R_O_01, STAT=', ALLOSTAT
      STOP
    END IF

#if(USE_MPI)
  IF ( myid == root ) THEN
#endif

    LOGSNP = .FALSE.
    IF  ( KEY_SNV )  THEN
      READ  (12, NML=SNAP)
      WRITE (11, NML=SNAP)
    END IF

#if(USE_MPI)
  END IF
#endif

  IF  ( KEY_TLV )  THEN

#if(USE_MPI)
    IF ( myid == root ) THEN
#endif
      READ  (12, NML=REC)
      WRITE (11, NML=REC)
#if(USE_MPI)
    END IF

    CALL MPI_Bcast(MR, 1, MPI_INTEGER, root, comm, errcode)
#endif

    ALLOCATE ( IREC (1:MR), KREC (1:MR) , LREC (1:MR), REC_NAME(1:MR), &
                                                         STAT=ALLOSTAT )
      IF  ( ALLOSTAT > 0 )  THEN
        WRITE (11,*) '  ALLOCATION ERROR R_O_02, STAT=', ALLOSTAT
        STOP
      END IF

#if(USE_MPI)
    IF ( myid == root ) THEN
#endif
      OPEN ( 22, FILE = REC_FILE_NAME, STATUS = 'OLD' )
      READ  (22,*)  ( REC_NAME(I), KREC(I), IREC(I), LREC(I), I=1,MR )
      CLOSE( 22 )

      LREC = ABS(LREC)

      IF (STRESS_IMAGING) THEN
        IF ( LPAS > 0 ) THEN
          LREC = LREC + 2*H/REFINE
        ELSE
          LREC = LREC + 2*H
        END IF
      END IF

      DO I = 1, MR
        IF ( LREC(I) >= LPAS*H/REFINE ) THEN
          KREC(I) = (KREC(I) - XBMIN) / H + 1
          IREC(I) = (IREC(I) - YBMIN) / H + 1
          LREC(I) = (LREC(I)-LPAS*H/REFINE) / H + LPAS
        ELSE
          KREC(I) = (KREC(I) - XBMIN) / H * REFINE + 1
          IREC(I) = (IREC(I) - YBMIN) / H * REFINE + 1
          LREC(I) =  LREC(I)          / H * REFINE
        END IF
      END DO

#if(USE_MPI)
    END IF
#endif
  END IF

#if(USE_MPI)
  IF ( myid == root ) THEN
#endif
    IF  ( KEY_SNV )  THEN
      DO I = 1, MSNP
        READ  (12,*)  LSNP
        IF (STRESS_IMAGING) THEN
          LOGSNP ( LSNP+2 ) = .TRUE.
        ELSE
          LOGSNP ( LSNP   ) = .TRUE.
        END IF
      END DO
    END IF

    CLOSE (12)

#if(USE_MPI)
#else
!___________________________________________ MODEL FILES _ 13,14,15,16 _

    INQUIRE ( FILE = MO_FILE_NAME, EXIST=EX, IOSTAT=IOS )
      IF  ( IOS > 0 )  THEN
        WRITE (11,*) 'ERROR HAS OCCURED DURING THE EXECUTION           &
                     &OF THE INQUIRE STATEMENT.',                      &
                     'IOSTAT (MO_FILE_NAME)=',IOS
        STOP
      END IF
      IF  ( .NOT.EX )  THEN
        WRITE (11,*) MO_FILE_NAME,' NOT FOUND'
        STOP
      END IF

    INQUIRE ( FILE = JMH_FILE_NAME, EXIST=EX, IOSTAT=IOS )
      IF  ( IOS > 0 )  THEN
        WRITE (11,*) 'ERROR HAS OCCURED DURING THE EXECUTION           &
                     &OF THE INQUIRE STATEMENT.',                      &
                     'IOSTAT (JMH_FILE_NAME)=',IOS
        STOP
      END IF
      IF  ( .NOT.EX )  THEN
        WRITE (11,*) JMH_FILE_NAME,' NOT FOUND'
        STOP
      END IF

    INQUIRE ( FILE =  Q_FILE_NAME, EXIST=EX, IOSTAT=IOS )
      IF  ( IOS > 0 )  THEN
        WRITE (11,*) 'ERROR HAS OCCURED DURING THE EXECUTION           &
                     &OF THE INQUIRE STATEMENT.',                      &
                     'IOSTAT (DE_FILE_NAME)=',IOS
        STOP
      END IF
      IF  ( .NOT.EX )  THEN
        WRITE (11,*)  Q_FILE_NAME,' NOT FOUND'
        STOP
      END IF

    OPEN ( 13, FILE   = JMH_FILE_NAME, FORM = 'UNFORMATTED',           &
                                                          STATUS='OLD' )

    OPEN ( 14, FILE = MO_FILE_NAME, FORM = 'UNFORMATTED', STATUS='OLD' )

    OPEN ( 15, FILE =  Q_FILE_NAME, FORM = 'UNFORMATTED', STATUS='OLD' )

    IF ( LPAS /= 0 ) THEN

      OPEN ( 16, FILE   = JMHF_FILE_NAME, FORM = 'UNFORMATTED',        &
                                                          STATUS='OLD' )

    END IF
#endif

#if(USE_MPI)
  END IF
#endif


END SUBROUTINE  READ_OPEN

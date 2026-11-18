!=======================================================================
! Copyright (c) 2014 Jozef Kristek, Peter Moczo
!
! Licensed under the EUPL (the "Licence")
! You may not use this work except in compliance with the Licence.
! You may obtain a copy of the Licence at
! https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12
! Unless required by applicable law or agreed to in writing, software
! distributed under the Licence is distributed on an "AS IS" basis,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and
! limitations under the Licence.
!
! REVISION HISTORY:
! Version 1.0, 2014-03-19: Jozef Kristek, Peter Moczo
!=======================================================================
!=======================================================================

PROGRAM FDModPar3D

!-----------------------------------------------------------------------

  IMPLICIT NONE

  LOGICAL :: STRESS_IMAGING, EX

  CHARACTER (LEN=17) :: JOB_NAME
  
  CHARACTER (LEN=20) ::  MO_FILE_NAME, JMH_FILE_NAME, Q_FILE_NAME,     &
                                      JMHF_FILE_NAME
  CHARACTER (LEN=3 )  CH_N

  INTEGER, PARAMETER :: IK = SELECTED_INT_KIND  (8)

  INTEGER :: I,K,L, J, JM1, JMNUM, MX, MY, MZ, LPAS, REFINE, MXF, MYF, &
             JMNUMN, TPML, PROCX, PROCY, PROCZ, N, JMNUMR

  REAL    :: H, DT, XBMIN, YBMIN
  
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: JM, JMF, JMN, JMFN
  INTEGER, ALLOCATABLE, DIMENSION(:)     :: COD

  REAL       , ALLOCATABLE, DIMENSION(:,:)   :: PAR , PARN
  REAL       , ALLOCATABLE, DIMENSION(:,:,:) :: PARY, PARYN

  INTEGER, ALLOCATABLE, DIMENSION(:) :: X_MIN_DA , X_MAX_DA ,          &
                                        Y_MIN_DA , Y_MAX_DA ,          &
                                        Z_MIN_DA , Z_MAX_DA ,          &
                                        X_MIN_DAF, X_MAX_DAF,          &
                                        Y_MIN_DAF, Y_MAX_DAF,          &
                                        Z_MIN_DAF, Z_MAX_DAF,          &
                                        X_MIN_A  , X_MAX_A  ,          &
                                        Y_MIN_A  , Y_MAX_A  ,          &
                                        Z_MIN_A  , Z_MAX_A  ,          &
                                        X_MIN_AF , X_MAX_AF ,          &
                                        Y_MIN_AF , Y_MAX_AF ,          &
                                        Z_MIN_AF , Z_MAX_AF

  INTEGER :: IOLENINT, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, X,Y,Z,      &
             PMLL, PMLXH, PMLYH, PMLZH, PMLLF, PMLXHF, PMLYHF, TPMLF

!-----------------------------------------------------------------------

  NAMELIST /FILENAMES/    MO_FILE_NAME,  JMH_FILE_NAME, Q_FILE_NAME,   &
                          JMHF_FILE_NAME

  NAMELIST /CONTROLDATA/  TPML, PROCX, PROCY, PROCZ, STRESS_IMAGING
  
  NAMELIST /CONTROLDATA2/ MX, MY, MZ, LPAS, H, DT, XBMIN, YBMIN, REFINE

!-----------------------------------------------------------------------

  OPEN  (10, FILE='HF_FDSim3D', STATUS='OLD' )
  READ  (10,'(A)') JOB_NAME
  CLOSE (10)
  
  OPEN  (12, FILE = TRIM(JOB_NAME)//'.INM', STATUS='OLD' )
  READ  (12, NML=FILENAMES)
  REFINE       = 1
  LPAS         = 0
  XBMIN        = 0.
  YBMIN        = 0.
  READ  (12, NML=CONTROLDATA2)
  CLOSE (12)
  
  OPEN  (12, FILE = TRIM(JOB_NAME)//'.IN', STATUS='OLD' )
  TPML         = 0
  STRESS_IMAGING = .FALSE.
  READ  (12, NML=CONTROLDATA)
  CLOSE (12)

  IF (STRESS_IMAGING) THEN
    IF ( LPAS /= 0 ) LPAS = LPAS + 2
    MZ   = MZ   + 2
  END IF

  PMLL   =   1 - TPML
  PMLXH  = MX  + TPML
  PMLYH  = MY  + TPML
  PMLZH  = MZ  + TPML

  IF ( LPAS /= 0 ) THEN
    TPMLF= TPML * REFINE

    MXF  = (MX-1)*REFINE+1
    MYF  = (MY-1)*REFINE+1


    PMLLF   =   1 - TPMLF
    PMLXHF  = MXF + TPMLF
    PMLYHF  = MYF + TPMLF
  END IF

  IF ( LPAS /= 0 ) THEN
    ALLOCATE ( JMF(PMLLF:PMLXHF,PMLLF:PMLYHF,0:LPAS) )

    OPEN ( 16, FILE   = JMHF_FILE_NAME, FORM = 'UNFORMATTED',          &
                                                          STATUS='OLD' )

    IF (STRESS_IMAGING) THEN
      DO L= 2, LPAS
        DO K = 1, MYF
          READ  ( 16 )  JMF( 1:MXF, K, L)
        END DO
      END DO
      JMF (:,:,1) = JMF (:,:,2)
      JMF (:,:,0) = JMF (:,:,3)
    ELSE
      DO L= 0, LPAS
        DO K = 1, MYF
          READ  ( 16 )  JMF( 1:MXF, K, L)
        END DO
      END DO
    END IF

    DO L= 0, LPAS
      DO J = 0, TPMLF
        JMF (1  -J, :    , L) = JMF (  1,   :, L)
        JMF (MXF+J, :    , L) = JMF (MXF,   :, L)
        JMF (:    , 1  -J, L) = JMF (  :,   1, L)
        JMF (:    , MYF+J, L) = JMF (  :, MYF, L)
      END DO
    END DO

    ALLOCATE ( JM (PMLL:PMLXH,PMLL:PMLYH,LPAS-2:PMLZH) )

    OPEN ( 13, FILE   = JMH_FILE_NAME, FORM = 'UNFORMATTED',           &
                                                          STATUS='OLD' )

    DO L = LPAS, MZ
      DO K = 1, MY 
        READ  ( 13 )  JM ( 1:MX , K, L)
      END DO
      DO J = 0, TPML
        JM (1 -J, :   , L) = JM ( 1,  :, L)
        JM (MX+J, :   , L) = JM (MX,  :, L)
        JM (:   , 1 -J, L) = JM ( :,  1, L)
        JM (:   , MY+J, L) = JM ( :, MY, L)
      END DO
    END DO

    DO L= MZ+1, PMLZH
      JM (:,:,L) = JM (:,:,MZ)
    END DO
    DO L= LPAS-2, LPAS-1
      JM (:,:,L) = JM (:,:,LPAS)
    END DO

  ELSE

    ALLOCATE ( JM (PMLL:PMLXH,PMLL:PMLYH,0:PMLZH) )

    OPEN ( 13, FILE   = JMH_FILE_NAME, FORM = 'UNFORMATTED',           &
                                                          STATUS='OLD' )


    IF (STRESS_IMAGING) THEN
      DO L= 2, MZ
        DO K = 1, MY
          READ  ( 13 )  JM ( 1:MX , K, L)
        END DO
      END DO
      JM (:,:,1) = JM (:,:,2)
      JM (:,:,0) = JM (:,:,3)
    ELSE
      DO L= 0, MZ
        DO K = 1, MY 
          READ  ( 13 )  JM ( 1:MX , K, L)
        END DO
      END DO
    END IF

    DO L = 0, MZ
      DO J = 0, TPML
        JM (1 -J, :   , L) = JM ( 1,  :, L)
        JM (MX+J, :   , L) = JM (MX,  :, L)
        JM (:   , 1 -J, L) = JM ( :,  1, L)
        JM (:   , MY+J, L) = JM ( :, MY, L)
      END DO
    END DO

    DO L= MZ+1, PMLZH
      JM (:,:,L) = JM (:,:,MZ)
    END DO

  END IF

!___________________________________________ MODEL FILES _ 13,14,15,16 _

  OPEN ( 14, FILE = MO_FILE_NAME, FORM = 'UNFORMATTED', STATUS='OLD' )
  OPEN ( 15, FILE =  Q_FILE_NAME, FORM = 'UNFORMATTED', STATUS='OLD' )

!___________________________________________ READ MO FILE

  READ ( 14 ) JMNUMR
  JMNUM = ABS(JMNUMR)
  
  IF ( JMNUMR > 0 ) THEN
    ALLOCATE ( PAR (JMNUM, 12), PARN (JMNUM, 12), COD(JMNUM) )
    READ  ( 14 ) ( PAR (JM1, 1), PAR (JM1, 2), PAR (JM1, 3),           &
                   PAR (JM1, 4), PAR (JM1, 5), PAR (JM1, 6),           &
                   PAR (JM1, 7), PAR (JM1, 8), PAR (JM1, 9),           &
                   PAR (JM1,10), PAR (JM1,11), PAR (JM1,12),           &
                   JM1 = 1, JMNUM)
  ELSE
    ALLOCATE ( PAR (JMNUM,  8), PARN (JMNUM,  8), COD(JMNUM) )
    READ  ( 14 ) ( PAR (JM1,1), PAR (JM1,2), PAR (JM1,3),              &
                   PAR (JM1,4), PAR (JM1,5), PAR (JM1,6),              &
                   PAR (JM1,7), PAR (JM1,8), JM1 = 1, JMNUM)
  END IF

!___________________________________________ READ Q FILE

  IF ( JMNUMR > 0 ) THEN
    ALLOCATE ( PARY (JMNUM, 9, 4), PARYN (JMNUM, 9, 4) )
    READ  ( 15 ) ( PARY(JM1,1,1:4), PARY(JM1,2,1:4), PARY(JM1,3,1:4),  &
                   PARY(JM1,4,1:4), PARY(JM1,5,1:4), PARY(JM1,6,1:4),  &
                   PARY(JM1,7,1:4), PARY(JM1,8,1:4), PARY(JM1,9,1:4),  &
                   JM1 = 1, JMNUM)
  ELSE
    ALLOCATE ( PARY (JMNUM, 5, 4), PARYN (JMNUM, 5, 4) )
    READ  ( 15 ) ( PARY(JM1,1,1:4), PARY(JM1,2,1:4), PARY(JM1,3,1:4),  &
                   PARY(JM1,4,1:4), PARY(JM1,5,1:4), JM1 = 1, JMNUM)
  END IF
  CLOSE ( 13 )
  CLOSE ( 14 )
  CLOSE ( 15 )

!_______________________________________________________________________

  CALL DOMAIN_DECOMPOSITION

  N = 0

  DO Z = 0, PROCZ-1
    DO Y = 0, PROCY-1
      DO X = 0, PROCX-1

        XMIN = X_MIN_DA(X)
        XMAX = X_MAX_DA(X)
        YMIN = Y_MIN_DA(Y)
        YMAX = Y_MAX_DA(Y)
        ZMIN = Z_MIN_DA(Z)
        ZMAX = Z_MAX_DA(Z)

        JMNUMN = 0
        COD    = 0

        ALLOCATE ( JMN( XMIN:XMAX, YMIN:YMAX, ZMIN:ZMAX ) )

        WRITE ( CH_N, '(I3.3)') N

        OPEN ( 23, FILE = TRIM(JMH_FILE_NAME)//CH_N, FORM = 'UNFORMATTED',   &
                                                         STATUS='NEW' )
        OPEN ( 24, FILE = TRIM( MO_FILE_NAME)//CH_N, FORM = 'UNFORMATTED',   &
                                                         STATUS='NEW' )
        OPEN ( 25, FILE = TRIM(  Q_FILE_NAME)//CH_N, FORM = 'UNFORMATTED',   &
                                                         STATUS='NEW' )

        DO L = ZMIN, ZMAX
          DO K = YMIN, YMAX
            DO I = XMIN, XMAX

              JM1 = JM (I,K,L)

              IF ( COD(JM1) > 0 ) THEN
                JMN (I,K,L) = COD(JM1)
              ELSE
                JMNUMN = JMNUMN + 1

                JMN (I,K,L) = JMNUMN

                PARN (JMNUMN, :   ) = PAR (JM1, :   )
                PARYN(JMNUMN, :, :) = PARY(JM1, :, :)

                COD(JM1) = JMNUMN
              END IF

            END DO
            WRITE ( 23 )  JMN ( :, K, L)
          END DO
        END DO

        CLOSE ( 23 )

        IF ( LPAS /= 0 ) THEN

          XMIN = X_MIN_DAF(X)
          XMAX = X_MAX_DAF(X)
          YMIN = Y_MIN_DAF(Y)
          YMAX = Y_MAX_DAF(Y)
          ZMIN = Z_MIN_DAF(Z)
          ZMAX = Z_MAX_DAF(Z)

          ALLOCATE ( JMFN( XMIN:XMAX, YMIN:YMAX, ZMIN:ZMAX ) )

          OPEN ( 26, FILE = TRIM(JMHF_FILE_NAME)//CH_N, FORM = 'UNFORMATTED',&
                                                          STATUS='NEW' )


          DO L = ZMIN, ZMAX
            DO K = YMIN, YMAX
              DO I = XMIN, XMAX

                JM1 = JMF (I,K,L)

                IF ( COD(JM1) > 0 ) THEN
                  JMFN (I,K,L) = COD(JM1)
                ELSE
                  JMNUMN = JMNUMN + 1

                  JMFN (I,K,L) = JMNUMN

                  PARN (JMNUMN, :   ) = PAR (JM1, :   )
                  PARYN(JMNUMN, :, :) = PARY(JM1, :, :)

                  COD(JM1) = JMNUMN
                END IF

              END DO
              WRITE ( 26 )  JMFN( :, K, L)
            END DO
          END DO

          CLOSE ( 26 )

          DEALLOCATE ( JMFN )

        END IF

        IF ( JMNUMR > 0 ) THEN
          WRITE ( 24 ) JMNUMN

          WRITE ( 24 ) ( PARN (JM1, 1), PARN (JM1, 2), PARN (JM1, 3),  &
                         PARN (JM1, 4), PARN (JM1, 5), PARN (JM1, 6),  &
                         PARN (JM1, 7), PARN (JM1, 8), PARN (JM1, 9),  &
                         PARN (JM1,10), PARN (JM1,11), PARN (JM1,12),  &
                         JM1 = 1, JMNUMN)

          WRITE ( 25 ) ( PARYN(JM1,1,1:4), PARYN(JM1,2,1:4), PARYN(JM1,3,1:4), &
                         PARYN(JM1,4,1:4), PARYN(JM1,5,1:4), PARYN(JM1,6,1:4), &
                         PARYN(JM1,7,1:4), PARYN(JM1,8,1:4), PARYN(JM1,9,1:4), &
                         JM1 = 1, JMNUMN)
        ELSE
          WRITE ( 24 ) -JMNUMN

          WRITE ( 24 ) ( PARN (JM1,1), PARN (JM1,2), PARN (JM1,3),       &
                         PARN (JM1,4), PARN (JM1,5), PARN (JM1,6),       &
                         PARN (JM1,7), PARN (JM1,8), JM1 = 1, JMNUMN)

          WRITE ( 25 ) ( PARYN(JM1,1,1:4), PARYN(JM1,2,1:4), PARYN(JM1,3,1:4), &
                         PARYN(JM1,4,1:4), PARYN(JM1,5,1:4), JM1 = 1, JMNUMN)
        END IF

        CLOSE ( 24 )
        CLOSE ( 25 )

        DEALLOCATE ( JMN )

        N = N + 1

      END DO
    END DO
  END DO

CONTAINS

  SUBROUTINE  DOMAIN_DECOMPOSITION

    INTEGER :: J, OVERLAP
    REAL    :: WL

!-----------------------------------------------------------------------

    OVERLAP = 2

    ALLOCATE ( X_MIN_A (0:PROCX-1), X_MAX_A (0:PROCX-1),               &
               Y_MIN_A (0:PROCY-1), Y_MAX_A (0:PROCY-1),               &
               Z_MIN_A (0:PROCZ-1), Z_MAX_A (0:PROCZ-1) )

    X_MIN_A(0      ) = PMLL
    X_MAX_A(PROCX-1) = PMLXH
    DO J = 0, PROCX-2
      X_MAX_A(J  ) = (J+1)*(PMLXH-PMLL+1)/PROCX + PMLL
      X_MIN_A(J+1) = (J+1)*(PMLXH-PMLL+1)/PROCX + PMLL + 1
    END DO

    Y_MIN_A(0      ) = PMLL
    Y_MAX_A(PROCY-1) = PMLYH
    DO J = 0, PROCY-2
      Y_MAX_A(J  ) = (J+1)*(PMLYH-PMLL+1)/PROCY + PMLL
      Y_MIN_A(J+1) = (J+1)*(PMLYH-PMLL+1)/PROCY + PMLL + 1
    END DO

    IF ( LPAS /= 0 ) THEN
      ALLOCATE ( X_MIN_AF (0:PROCX-1), X_MAX_AF (0:PROCX-1),           &
                 Y_MIN_AF (0:PROCY-1), Y_MAX_AF (0:PROCY-1),           &
                 Z_MIN_AF (0:PROCZ-1), Z_MAX_AF (0:PROCZ-1) )

      X_MIN_AF(0      ) = PMLLF
      X_MAX_AF(PROCX-1) = PMLXHF
      DO J = 0, PROCX-2
        X_MAX_AF(J  ) = (X_MAX_A (J)-1)*REFINE+1
        X_MIN_AF(J+1) =  X_MAX_AF(J) + 1
      END DO

      Y_MIN_AF(0      ) = PMLLF
      Y_MAX_AF(PROCY-1) = PMLYHF
      DO J = 0, PROCY-2
        Y_MAX_AF(J  ) = (Y_MAX_A (J)-1)*REFINE+1
        Y_MIN_AF(J+1) =  Y_MAX_AF(J) + 1
      END DO

      Z_MIN_AF(0      ) = 0
      Z_MAX_AF(PROCZ-1) = LPAS
      Z_MIN_A (0      ) = LPAS-2
      Z_MAX_A (PROCZ-1) = PMLZH
      WL = ((PMLZH-LPAS+1)/REFINE+(LPAS+1)*REFINE)/PROCZ
      DO J = 0, PROCZ-2
        Z_MAX_AF(J  ) = MIN(NINT((J+1)*WL/REFINE),LPAS)
        Z_MIN_AF(J+1) = MIN(NINT((J+1)*WL/REFINE),LPAS) + 1
      END DO
      DO J = PROCZ-1,1,-1
        Z_MIN_A (J  ) = MAX(PMLZH-NINT((PROCZ-J)*WL*REFINE),LPAS-2)
        Z_MAX_A (J-1) = MAX(PMLZH-NINT((PROCZ-J)*WL*REFINE),LPAS-2)-1
      END DO

    ELSE
      Z_MIN_A(0      ) = 0
      Z_MAX_A(PROCZ-1) = PMLZH
      DO J = 0, PROCZ-2
        Z_MAX_A(J  ) = (J+1)*(PMLZH     +1)/PROCZ
        Z_MIN_A(J+1) = (J+1)*(PMLZH     +1)/PROCZ + 1
      END DO
    END IF

    ALLOCATE ( X_MIN_DA (0:PROCX-1), X_MAX_DA (0:PROCX-1),             &
               Y_MIN_DA (0:PROCY-1), Y_MAX_DA (0:PROCY-1),             &
               Z_MIN_DA (0:PROCZ-1), Z_MAX_DA (0:PROCZ-1) )

    X_MIN_DA(0      ) = PMLL
    X_MAX_DA(PROCX-1) = PMLXH
    DO J = 0, PROCX-2
      X_MAX_DA(J  ) = X_MAX_A(J  ) + OVERLAP
      X_MIN_DA(J+1) = X_MIN_A(J+1) - OVERLAP
    END DO

    Y_MIN_DA(0        ) = PMLL
    Y_MAX_DA(PROCY-1) = PMLYH
    DO J = 0, PROCY-2
      Y_MAX_DA(J  ) = Y_MAX_A(J  ) + OVERLAP
      Y_MIN_DA(J+1) = Y_MIN_A(J+1) - OVERLAP
    END DO

    IF ( LPAS /= 0 ) THEN
      ALLOCATE ( X_MIN_DAF (0:PROCX-1), X_MAX_DAF (0:PROCX-1),         &
                 Y_MIN_DAF (0:PROCY-1), Y_MAX_DAF (0:PROCY-1),         &
                 Z_MIN_DAF (0:PROCZ-1), Z_MAX_DAF (0:PROCZ-1) )

      X_MIN_DAF(0      ) = PMLLF
      X_MAX_DAF(PROCX-1) = PMLXHF
      DO J = 0, PROCX-2
        X_MAX_DAF(J  ) = X_MAX_AF(J  ) + OVERLAP
        X_MIN_DAF(J+1) = X_MIN_AF(J+1) - OVERLAP
      END DO

      Y_MIN_DAF(0      ) = PMLLF
      Y_MAX_DAF(PROCY-1) = PMLYHF
      DO J = 0, PROCY-2
        Y_MAX_DAF(J  ) = Y_MAX_AF(J  ) + OVERLAP
        Y_MIN_DAF(J+1) = Y_MIN_AF(J+1) - OVERLAP
      END DO

      Z_MIN_DAF(0      ) = 0
      Z_MAX_DAF(PROCZ-1) = LPAS
      Z_MIN_DA (0      ) = LPAS-2
      Z_MAX_DA (PROCZ-1) = PMLZH
      WL = ((PMLZH-LPAS+1)/REFINE+(LPAS+1)*REFINE)/PROCZ
      DO J = 0, PROCZ-2
        Z_MAX_DAF(J  ) = MIN(NINT((J+1)*WL/REFINE) + OVERLAP, LPAS)
        Z_MIN_DAF(J+1) = MIN(NINT((J+1)*WL/REFINE) - OVERLAP, LPAS) + 1
      END DO
      DO J = PROCZ-1,1,-1
        Z_MIN_DA (J  ) = MAX(PMLZH-NINT((PROCZ-J)*WL*REFINE)+1-OVERLAP,LPAS-2)
        Z_MAX_DA (J-1) = MAX(PMLZH-NINT((PROCZ-J)*WL*REFINE)+1+OVERLAP,LPAS-2)-1
      END DO
    ELSE
      Z_MIN_DA(0      ) = 0
      Z_MAX_DA(PROCZ-1) = PMLZH
      DO J = 0, PROCZ-2
        Z_MAX_DA(J  ) = (J+1)*(PMLZH     +1)/PROCZ     + OVERLAP
        Z_MIN_DA(J+1) = (J+1)*(PMLZH     +1)/PROCZ + 1 - OVERLAP
      END DO
    END IF

  END SUBROUTINE  DOMAIN_DECOMPOSITION

END PROGRAM

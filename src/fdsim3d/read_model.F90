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

SUBROUTINE READ_MODEL

  USE PRECISION   , ONLY:                                              &
                          IK
  USE CONTROL_DATA, ONLY:                                              &
                          MX, MY, MZ, TPML, LPAS,                      &
                          MO_FILE_NAME, JMH_FILE_NAME, Q_FILE_NAME,    &
                                       JMHF_FILE_NAME,                 &
                          STRESS_IMAGING
  USE GRID_MEDIUM , ONLY:                                              &
                          JM, JMF
  USE AUXIL       , ONLY :                                             &
                          PMLL , PMLXH , PMLYH , PMLZH,                &
                          PMLLF, PMLXHF, PMLYHF, TPMLF,                &
                          X_MIN_D , X_MAX_D,                           &
                          Y_MIN_D , Y_MAX_D,                           &
                          Z_MIN_D , Z_MAX_D,                           &
                          X_MIN_DF, X_MAX_DF,                          &
                          Y_MIN_DF, Y_MAX_DF,                          &
                          Z_MIN_DF, Z_MAX_DF, MXF, MYF
#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif

!-----------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER (LEN=3) :: my_N_CH

  INTEGER :: K, L, J, SIZE
  INTEGER(IK), ALLOCATABLE, DIMENSION(:,:,:) :: JMT, JMTF

!-----------------------------------------------------------------------

  IF ( LPAS /= 0 ) THEN
    ALLOCATE ( JMF( X_MIN_DF:X_MAX_DF,                                 &
                    Y_MIN_DF:Y_MAX_DF,                                 &
                    Z_MIN_DF:Z_MAX_DF) )
  END IF

  ALLOCATE ( JM (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) )

#if(USE_MPI)
  my_N = ( my_coord(3)*dims(2) + my_coord(2) )*dims(1) + my_coord(1)

  WRITE (my_N_CH, '(I3.3)') my_N

  OPEN (13,FILE = TRIM(JMH_FILE_NAME)//my_N_CH, FORM = 'UNFORMATTED',&
                                                      STATUS='OLD' )

  DO L = Z_MIN_D, Z_MAX_D
    DO K = Y_MIN_D, Y_MAX_D
      READ ( 13 )  JM ( :, K, L)
    END DO
  END DO

  OPEN ( 14, FILE =  TRIM(MO_FILE_NAME)//my_N_CH, FORM = 'UNFORMATTED',&
                                                        STATUS='OLD' )

  OPEN ( 15, FILE =   TRIM(Q_FILE_NAME)//my_N_CH, FORM = 'UNFORMATTED',&
                                                        STATUS='OLD' )

  IF ( LPAS /= 0 ) THEN

    OPEN(16,FILE=TRIM(JMHF_FILE_NAME)//my_N_CH, FORM = 'UNFORMATTED',&
                                                      STATUS='OLD' )

    DO L = Z_MIN_DF, Z_MAX_DF
      DO K = Y_MIN_DF, Y_MAX_DF
        READ ( 16 )  JMF( :, K, L)
      END DO
    END DO

  END IF


#else

  IF ( LPAS /= 0 ) THEN
    ALLOCATE ( JMTF(PMLLF:PMLXHF,PMLLF:PMLYHF,0:LPAS) )

    IF (STRESS_IMAGING) THEN
      DO L= 2, LPAS

        DO K = 1, MYF
          READ  ( 16 )  JMTF( 1:MXF, K, L)
        END DO

        JMTF (:,:,1) = JMTF (:,:,2)
        JMTF (:,:,0) = JMTF (:,:,3)

        DO J = 0, TPMLF
          JMTF (1  -J, :    , L) = JMTF (  1,   :, L)
          JMTF (MXF+J, :    , L) = JMTF (MXF,   :, L)
          JMTF (:    , 1  -J, L) = JMTF (  :,   1, L)
          JMTF (:    , MYF+J, L) = JMTF (  :, MYF, L)
        END DO

      END DO
    ELSE
      DO L= 0, LPAS

        DO K = 1, MYF
          READ  ( 16 )  JMTF( 1:MXF, K, L)
        END DO

        DO J = 0, TPMLF
          JMTF (1  -J, :    , L) = JMTF (  1,   :, L)
          JMTF (MXF+J, :    , L) = JMTF (MXF,   :, L)
          JMTF (:    , 1  -J, L) = JMTF (  :,   1, L)
          JMTF (:    , MYF+J, L) = JMTF (  :, MYF, L)
        END DO

      END DO
    END IF

    ALLOCATE ( JMT (PMLL:PMLXH,PMLL:PMLYH,LPAS-2:PMLZH) )

    DO L=LPAS, MZ

      DO K = 1, MY
        READ ( 13)  JMT ( 1:MX, K, L)
      END DO

      DO J = 0, TPML
        JMT (1 -J, :   , L) = JMT ( 1,  :, L)
        JMT (MX+J, :   , L) = JMT (MX,  :, L)
        JMT (:   , 1 -J, L) = JMT ( :,  1, L)
        JMT (:   , MY+J, L) = JMT ( :, MY, L)
      END DO

    END DO

    DO L= MZ+1, PMLZH
      JMT (:,:,L) = JMT (:,:,MZ)
    END DO
    DO L= LPAS-2, LPAS-1
      JMT (:,:,L) = JMT (:,:,LPAS)
    END DO

  ELSE

    ALLOCATE ( JMT (PMLL:PMLXH,PMLL:PMLYH,0:PMLZH) )

    IF (STRESS_IMAGING) THEN
      DO L= 2, MZ

        DO K = 1, MY
          READ  ( 13 )  JMT ( 1:MX, K, L)
        END DO

        JMT (:,:,1) = JMT (:,:,2)
        JMT (:,:,0) = JMT (:,:,3)
        DO J = 0, TPML
          JMT (1 -J, :   , L) = JMT ( 1,  :, L)
          JMT (MX+J, :   , L) = JMT (MX,  :, L)
          JMT (:   , 1 -J, L) = JMT ( :,  1, L)
          JMT (:   , MY+J, L) = JMT ( :, MY, L)
        END DO

      END DO
    ELSE
      DO L= 0, MZ

        DO K = 1, MY
          READ ( 13 )  JMT ( 1:MX, K, L)
        END DO

        DO J = 0, TPML
          JMT (1 -J, :   , L) = JMT ( 1,  :, L)
          JMT (MX+J, :   , L) = JMT (MX,  :, L)
          JMT (:   , 1 -J, L) = JMT ( :,  1, L)
          JMT (:   , MY+J, L) = JMT ( :, MY, L)
        END DO

      END DO
    END IF

    DO L= MZ+1, PMLZH
      JMT (:,:,L) = JMT (:,:,MZ)
    END DO

  END IF

  JM = JMT

  IF ( LPAS /= 0 ) THEN
    JMF = JMTF
  END IF

  DEALLOCATE(JMT)

  IF ( LPAS /= 0 ) THEN
    DEALLOCATE(JMTF)
  END IF
#endif


END SUBROUTINE READ_MODEL

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

SUBROUTINE COMMUNICATE_PLANE (VARIABLE)

  USE PRECISION   , ONLY : PP
  USE AUXIL       , ONLY : X_MIN_D, X_MAX_D, X_MIN, X_MAX,             &
                           Y_MIN_D, Y_MAX_D, Y_MIN, Y_MAX,             &
                           BUFXS_P, BUFYS_P,                           &
                           BUFXR_P, BUFYR_P,                           &
                           source_id0A, dest_id0A,                     &
                           source_id0B, dest_id0B,                     &
                           source_id1A, dest_id1A,                     &
                           source_id1B, dest_id1B

#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif
!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL(PP), INTENT(INOUT),                                             &
                DIMENSION(X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D) :: VARIABLE


  INTEGER       :: stag, rtag, SIZE

!-----------------------------------------------------------------------

#if(USE_MPI)
!!  CALL MPI_BARRIER (comm_cart,errcode)

! --- communicate in x direction

  SIZE = OVERLAP*(Y_MAX-Y_MIN+1)

  IF ( SIZE > 0 )  THEN

    stag = 0
    rtag = 0
    IF ( dest_id0A /= MPI_PROC_NULL ) THEN
      BUFXS_P (1:OVERLAP,1:Y_MAX-Y_MIN+1) =                            &
                            VARIABLE (X_MAX+1-OVERLAP:X_MAX,Y_MIN:Y_MAX)
      IF ( source_id0A == MPI_PROC_NULL ) THEN
        CALL MPI_SEND            (BUFXS_P, SIZE, MPI_REAL,             &
                                   dest_id0A, stag,                    &
                                   comm_cart, errcode)
      ELSE
        CALL MPI_ISEND           (BUFXS_P, SIZE, MPI_REAL,             &
                                   dest_id0A, stag,                    &
                                   comm_cart, request, errcode)
        CALL MPI_RECV            (BUFXR_P, SIZE, MPI_REAL,             &
                                                    source_id0A, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE (X_MIN_D:X_MIN_D+OVERLAP-1,Y_MIN:Y_MAX) =             &
                                    BUFXR_P (1:OVERLAP,1:Y_MAX-Y_MIN+1)
        CALL MPI_WAIT (request, status, errcode)
      END IF
    ELSE
      IF ( source_id0A /= MPI_PROC_NULL ) THEN
        CALL MPI_RECV            (BUFXR_P, SIZE, MPI_REAL,             &
                                                    source_id0A, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE (X_MIN_D:X_MIN_D+OVERLAP-1,Y_MIN:Y_MAX) =             &
                                    BUFXR_P (1:OVERLAP,1:Y_MAX-Y_MIN+1)
      END IF
    END IF

    stag = 1
    rtag = 1
    IF ( dest_id0B /= MPI_PROC_NULL ) THEN
      BUFXS_P = VARIABLE (X_MIN:X_MIN+OVERLAP-1,Y_MIN:Y_MAX)
      IF ( source_id0B == MPI_PROC_NULL ) THEN
        CALL MPI_SEND            (BUFXS_P, SIZE, MPI_REAL,             &
                                   dest_id0B, stag,                    &
                                   comm_cart, errcode)
      ELSE
        CALL MPI_ISEND           (BUFXS_P, SIZE, MPI_REAL,             &
                                   dest_id0B, stag,                    &
                                   comm_cart, request, errcode)
        CALL MPI_RECV            (BUFXR_P, SIZE, MPI_REAL,             &
                                                    source_id0B, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE (X_MAX_D+1-OVERLAP:X_MAX_D,Y_MIN:Y_MAX) = BUFXR_P
        CALL MPI_WAIT (request, status, errcode)
      END IF
    ELSE
      IF ( source_id0B /= MPI_PROC_NULL ) THEN
        CALL MPI_RECV            (BUFXR_P, SIZE, MPI_REAL,             &
                                                    source_id0B, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE (X_MAX_D+1-OVERLAP:X_MAX_D,Y_MIN:Y_MAX) = BUFXR_P
      END IF
    END IF

  END IF


! --- communicate in y direction

  SIZE = (X_MAX_D-X_MIN_D+1)*OVERLAP

  IF ( SIZE > 0 )  THEN

    stag = 2
    rtag = 2
    IF ( dest_id1A /= MPI_PROC_NULL ) THEN
      BUFYS_P = VARIABLE (X_MIN_D:X_MAX_D,Y_MAX+1-OVERLAP:Y_MAX)
      IF ( source_id1A == MPI_PROC_NULL ) THEN
        CALL MPI_SEND            (BUFYS_P, SIZE, MPI_REAL,             &
                                   dest_id1A, stag,                    &
                                   comm_cart, errcode)
      ELSE
        CALL MPI_ISEND           (BUFYS_P, SIZE, MPI_REAL,             &
                                   dest_id1A, stag,                    &
                                   comm_cart, request, errcode)
        CALL MPI_RECV            (BUFYR_P, SIZE, MPI_REAL,             &
                                                    source_id1A, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE(X_MIN_D:X_MAX_D,Y_MIN_D:Y_MIN_D+OVERLAP-1) = BUFYR_P
        CALL MPI_WAIT (request, status, errcode)
      END IF
    ELSE
      IF ( source_id1A /= MPI_PROC_NULL ) THEN
        CALL MPI_RECV            (BUFYR_P, SIZE, MPI_REAL,             &
                                                    source_id1A, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE(X_MIN_D:X_MAX_D,Y_MIN_D:Y_MIN_D+OVERLAP-1) = BUFYR_P
      END IF
    END IF


    stag = 3
    rtag = 3
    IF ( dest_id1B /= MPI_PROC_NULL ) THEN
      BUFYS_P = VARIABLE (X_MIN_D:X_MAX_D,Y_MIN:Y_MIN+OVERLAP-1)
      IF ( source_id1B == MPI_PROC_NULL ) THEN
        CALL MPI_SEND            (BUFYS_P, SIZE, MPI_REAL,             &
                                   dest_id1B, stag,                    &
                                   comm_cart, errcode)
      ELSE
        CALL MPI_ISEND           (BUFYS_P, SIZE, MPI_REAL,             &
                                   dest_id1B, stag,                    &
                                   comm_cart, request, errcode)
        CALL MPI_RECV            (BUFYR_P, SIZE, MPI_REAL,             &
                                                    source_id1B, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE(X_MIN_D:X_MAX_D,Y_MAX_D+1-OVERLAP:Y_MAX_D) = BUFYR_P
        CALL MPI_WAIT (request, status, errcode)
      END IF
    ELSE
      IF ( source_id1B /= MPI_PROC_NULL ) THEN
        CALL MPI_RECV            (BUFYR_P, SIZE, MPI_REAL,             &
                                                    source_id1B, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE(X_MIN_D:X_MAX_D,Y_MAX_D+1-OVERLAP:Y_MAX_D) = BUFYR_P
      END IF
    END IF

  END IF

#endif


END SUBROUTINE COMMUNICATE_PLANE

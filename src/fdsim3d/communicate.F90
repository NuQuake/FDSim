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

SUBROUTINE COMMUNICATE (VARIABLE)

  USE PRECISION   , ONLY : PP
  USE CONTROL_DATA, ONLY : LPAS
  USE AUXIL       , ONLY : X_MIN_D, X_MAX_D, X_MIN, X_MAX,             &
                           Y_MIN_D, Y_MAX_D, Y_MIN, Y_MAX,             &
                           Z_MIN_D, Z_MAX_D, Z_MIN, Z_MAX,             &
                           BUFXS, BUFYS, BUFZS,                        &
                           BUFXR, BUFYR, BUFZR,                        &
                           source_id0A, dest_id0A,                     &
                           source_id0B, dest_id0B,                     &
                           source_id1A, dest_id1A,                     &
                           source_id1B, dest_id1B,                     &
                           source_id2A, dest_id2A,                     &
                           source_id2B, dest_id2B

#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif
!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL(PP), INTENT(INOUT),                                             &
    DIMENSION(X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D, Z_MIN_D:Z_MAX_D) ::    &
                                                              VARIABLE


  INTEGER       :: stag, rtag, SIZE

!-----------------------------------------------------------------------

#if(USE_MPI)
!!  CALL MPI_BARRIER (comm_cart,errcode)

! --- communicate in x direction

  SIZE = OVERLAP*(Y_MAX-Y_MIN+1)*(Z_MAX-Z_MIN+1)

  IF ( SIZE > 0 )  THEN

    stag = 0
    rtag = 0
    IF ( dest_id0A /= MPI_PROC_NULL ) THEN
      BUFXS(1:OVERLAP,1:Y_MAX-Y_MIN+1,1:Z_MAX-Z_MIN+1) =               &
        VARIABLE (X_MAX+1-OVERLAP:X_MAX,Y_MIN:Y_MAX,Z_MIN:Z_MAX)
      IF ( source_id0A == MPI_PROC_NULL ) THEN
        CALL MPI_SEND            (BUFXS, SIZE, MPI_REAL,               &
                                   dest_id0A, stag,                    &
                                   comm_cart, errcode)
      ELSE
!        CALL MPI_SENDRECV (BUFXS, SIZE, MPI_REAL,   dest_id0A, stag,   &
!                           BUFXR, SIZE, MPI_REAL, source_id0A, rtag,   &
!                                   comm_cart, status, errcode)
        CALL MPI_ISEND           (BUFXS, SIZE, MPI_REAL,               &
                                   dest_id0A, stag,                    &
                                   comm_cart, request, errcode)
        CALL MPI_RECV            (BUFXR, SIZE, MPI_REAL,               &
                                                    source_id0A, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE (X_MIN_D:X_MIN_D+OVERLAP-1,Y_MIN:Y_MAX,Z_MIN:Z_MAX) = &
          BUFXR(1:OVERLAP,1:Y_MAX-Y_MIN+1,1:Z_MAX-Z_MIN+1)
        CALL MPI_WAIT (request, status, errcode)
      END IF
    ELSE
      IF ( source_id0A /= MPI_PROC_NULL ) THEN
        CALL MPI_RECV            (BUFXR, SIZE, MPI_REAL,               &
                                                    source_id0A, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE (X_MIN_D:X_MIN_D+OVERLAP-1,Y_MIN:Y_MAX,Z_MIN:Z_MAX) = &
          BUFXR(1:OVERLAP,1:Y_MAX-Y_MIN+1,1:Z_MAX-Z_MIN+1)
      END IF
    END IF

    stag = 1
    rtag = 1
    IF ( dest_id0B /= MPI_PROC_NULL ) THEN
      BUFXS = VARIABLE (X_MIN:X_MIN+OVERLAP-1,Y_MIN:Y_MAX,Z_MIN:Z_MAX)
      IF ( source_id0B == MPI_PROC_NULL ) THEN
        CALL MPI_SEND            (BUFXS, SIZE, MPI_REAL,               &
                                   dest_id0B, stag,                    &
                                   comm_cart, errcode)
      ELSE
!        CALL MPI_SENDRECV (BUFXS, SIZE, MPI_REAL,   dest_id0B, stag,   &
!                           BUFXR, SIZE, MPI_REAL, source_id0B, rtag,   &
!                                   comm_cart, status, errcode)
        CALL MPI_ISEND           (BUFXS, SIZE, MPI_REAL,               &
                                   dest_id0B, stag,                    &
                                   comm_cart, request, errcode)
        CALL MPI_RECV            (BUFXR, SIZE, MPI_REAL,               &
                                                    source_id0B, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE (X_MAX_D+1-OVERLAP:X_MAX_D,Y_MIN:Y_MAX,Z_MIN:Z_MAX) = &
          BUFXR
        CALL MPI_WAIT (request, status, errcode)
      END IF
    ELSE
      IF ( source_id0B /= MPI_PROC_NULL ) THEN
        CALL MPI_RECV            (BUFXR, SIZE, MPI_REAL,               &
                                                    source_id0B, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE (X_MAX_D+1-OVERLAP:X_MAX_D,Y_MIN:Y_MAX,Z_MIN:Z_MAX) = &
          BUFXR
      END IF
    END IF

  END IF


! --- communicate in y direction

  SIZE = (X_MAX_D-X_MIN_D+1)*OVERLAP*(Z_MAX-Z_MIN+1)

  IF ( SIZE > 0 )  THEN

    stag = 2
    rtag = 2
    IF ( dest_id1A /= MPI_PROC_NULL ) THEN
      BUFYS = VARIABLE (X_MIN_D:X_MAX_D,Y_MAX+1-OVERLAP:Y_MAX,Z_MIN:Z_MAX)
      IF ( source_id1A == MPI_PROC_NULL ) THEN
        CALL MPI_SEND            (BUFYS, SIZE, MPI_REAL,               &
                                   dest_id1A, stag,                    &
                                   comm_cart, errcode)
      ELSE
!        CALL MPI_SENDRECV (BUFYS, SIZE, MPI_REAL,   dest_id1A, stag,   &
!                           BUFYR, SIZE, MPI_REAL, source_id1A, rtag,   &
!                                   comm_cart, status, errcode)
        CALL MPI_ISEND           (BUFYS, SIZE, MPI_REAL,               &
                                   dest_id1A, stag,                    &
                                   comm_cart, request, errcode)
        CALL MPI_RECV            (BUFYR, SIZE, MPI_REAL,               &
                                                    source_id1A, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE(X_MIN_D:X_MAX_D,Y_MIN_D:Y_MIN_D+OVERLAP-1,Z_MIN:Z_MAX)&
          = BUFYR
        CALL MPI_WAIT (request, status, errcode)
      END IF
    ELSE
      IF ( source_id1A /= MPI_PROC_NULL ) THEN
        CALL MPI_RECV            (BUFYR, SIZE, MPI_REAL,               &
                                                    source_id1A, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE(X_MIN_D:X_MAX_D,Y_MIN_D:Y_MIN_D+OVERLAP-1,Z_MIN:Z_MAX)&
          = BUFYR
      END IF
    END IF


    stag = 3
    rtag = 3
    IF ( dest_id1B /= MPI_PROC_NULL ) THEN
      BUFYS = VARIABLE (X_MIN_D:X_MAX_D,Y_MIN:Y_MIN+OVERLAP-1,Z_MIN:Z_MAX)
      IF ( source_id1B == MPI_PROC_NULL ) THEN
        CALL MPI_SEND            (BUFYS, SIZE, MPI_REAL,               &
                                   dest_id1B, stag,                    &
                                   comm_cart, errcode)
      ELSE
!        CALL MPI_SENDRECV (BUFYS, SIZE, MPI_REAL,   dest_id1B, stag,   &
!                           BUFYR, SIZE, MPI_REAL, source_id1B, rtag,   &
!                                   comm_cart, status, errcode)
        CALL MPI_ISEND           (BUFYS, SIZE, MPI_REAL,               &
                                   dest_id1B, stag,                    &
                                   comm_cart, request, errcode)
        CALL MPI_RECV            (BUFYR, SIZE, MPI_REAL,               &
                                                    source_id1B, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE(X_MIN_D:X_MAX_D,Y_MAX_D+1-OVERLAP:Y_MAX_D,Z_MIN:Z_MAX)&
          = BUFYR
        CALL MPI_WAIT (request, status, errcode)
      END IF
    ELSE
      IF ( source_id1B /= MPI_PROC_NULL ) THEN
        CALL MPI_RECV            (BUFYR, SIZE, MPI_REAL,               &
                                                    source_id1B, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE(X_MIN_D:X_MAX_D,Y_MAX_D+1-OVERLAP:Y_MAX_D,Z_MIN:Z_MAX)&
          = BUFYR
      END IF
    END IF

  END IF

! --- communicate in z direction

  SIZE = (X_MAX_D-X_MIN_D+1)*(Y_MAX_D-Y_MIN_D+1)*OVERLAP

  IF ( Z_MIN <= Z_MAX ) THEN

    stag = 4
    rtag = 4
    IF ( dest_id2A /= MPI_PROC_NULL ) THEN
      BUFZS = VARIABLE                                                 &
                (X_MIN_D:X_MAX_D,Y_MIN_D:Y_MAX_D,Z_MAX+1-OVERLAP:Z_MAX)
      IF ( source_id2A == MPI_PROC_NULL ) THEN
        CALL MPI_SEND            (BUFZS, SIZE, MPI_REAL,               &
                                   dest_id2A, stag,                    &
                                   comm_cart, errcode)
      ELSE
!        CALL MPI_SENDRECV (BUFZS, SIZE, MPI_REAL,   dest_id2A, stag,   &
!                           BUFZR, SIZE, MPI_REAL, source_id2A, rtag,   &
!                                   comm_cart, status, errcode)
        CALL MPI_ISEND           (BUFZS, SIZE, MPI_REAL,               &
                                   dest_id2A, stag,                    &
                                   comm_cart, request, errcode)
        CALL MPI_RECV            (BUFZR, SIZE, MPI_REAL,               &
                                                    source_id2A, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE                                                       &
         (X_MIN_D:X_MAX_D,Y_MIN_D:Y_MAX_D,Z_MIN_D:Z_MIN_D+OVERLAP-1) = &
         BUFZR
        CALL MPI_WAIT (request, status, errcode)
      END IF
    ELSE
      IF ( source_id2A /= MPI_PROC_NULL ) THEN
        CALL MPI_RECV            (BUFZR, SIZE, MPI_REAL,               &
                                                    source_id2A, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE                                                       &
         (X_MIN_D:X_MAX_D,Y_MIN_D:Y_MAX_D,Z_MIN_D:Z_MIN_D+OVERLAP-1) = &
         BUFZR
      END IF
    END IF


    stag = 5
    rtag = 5
    IF ( dest_id2B /= MPI_PROC_NULL ) THEN
      BUFZS = VARIABLE                                                 &
                (X_MIN_D:X_MAX_D,Y_MIN_D:Y_MAX_D,Z_MIN:Z_MIN+OVERLAP-1)
      IF ( source_id2B == MPI_PROC_NULL ) THEN
        CALL MPI_SEND            (BUFZS, SIZE, MPI_REAL,               &
                                   dest_id2B, stag,                    &
                                   comm_cart, errcode)
      ELSE
!        CALL MPI_SENDRECV (BUFZS, SIZE, MPI_REAL,   dest_id2B, stag,   &
!                           BUFZR, SIZE, MPI_REAL, source_id2B, rtag,   &
!                                   comm_cart, status, errcode)
        CALL MPI_ISEND           (BUFZS, SIZE, MPI_REAL,               &
                                   dest_id2B, stag,                    &
                                   comm_cart, request, errcode)
        CALL MPI_RECV            (BUFZR, SIZE, MPI_REAL,               &
                                                    source_id2B, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE                                                       &
           (X_MIN_D:X_MAX_D,Y_MIN_D:Y_MAX_D,Z_MAX_D+1-OVERLAP:Z_MAX_D) &
           = BUFZR
        CALL MPI_WAIT (request, status, errcode)
      END IF
    ELSE
      IF ( source_id2B /= MPI_PROC_NULL ) THEN
        CALL MPI_RECV            (BUFZR, SIZE, MPI_REAL,               &
                                                    source_id2B, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE                                                       &
           (X_MIN_D:X_MAX_D,Y_MIN_D:Y_MAX_D,Z_MAX_D+1-OVERLAP:Z_MAX_D) &
           = BUFZR
      END IF
    END IF

  END IF
#endif


END SUBROUTINE COMMUNICATE

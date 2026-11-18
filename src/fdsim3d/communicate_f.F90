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

SUBROUTINE COMMUNICATE_F (VARIABLE)

  USE PRECISION   , ONLY : PP
  USE CONTROL_DATA, ONLY : LPAS
  USE AUXIL       , ONLY : X_MIN_DF, X_MAX_DF, X_MINF, X_MAXF,         &
                           Y_MIN_DF, Y_MAX_DF, Y_MINF, Y_MAXF,         &
                           Z_MIN_DF, Z_MAX_DF, Z_MINF, Z_MAXF,         &
                           BUFXSF, BUFYSF, BUFZSF,                     &
                           BUFXRF, BUFYRF, BUFZRF,                     &
                           source_id0A , dest_id0A ,                   &
                           source_id0B , dest_id0B ,                   &
                           source_id1A , dest_id1A ,                   &
                           source_id1B , dest_id1B ,                   &
                           source_id2AF, dest_id2AF,                   &
                           source_id2BF, dest_id2BF
#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif
!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL(PP), INTENT(INOUT), DIMENSION                                   &
       (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF, Z_MIN_DF:Z_MAX_DF) ::    &
                                                              VARIABLE


  INTEGER       :: stag, rtag, SIZE

!-----------------------------------------------------------------------

#if(USE_MPI)
!!  CALL MPI_BARRIER (comm_cart,errcode)

! --- communicate in x direction

  SIZE = OVERLAP*(Y_MAXF-Y_MINF+1)*(Z_MAXF-Z_MINF+1)

  IF ( SIZE > 0 )  THEN

    stag = 0
    rtag = 0
    IF ( dest_id0A /= MPI_PROC_NULL ) THEN
      BUFXSF(1:OVERLAP,1:Y_MAXF-Y_MINF+1,1:Z_MAXF-Z_MINF+1) =          &
        VARIABLE (X_MAXF+1-OVERLAP:X_MAXF,Y_MINF:Y_MAXF,Z_MINF:Z_MAXF)
      IF ( source_id0A == MPI_PROC_NULL ) THEN
        CALL MPI_SEND            (BUFXSF, SIZE, MPI_REAL,              &
                                   dest_id0A, stag,                    &
                               comm_cart, errcode)
      ELSE
!        CALL MPI_SENDRECV (BUFXSF, SIZE, MPI_REAL,   dest_id0A, stag,  &
!                           BUFXRF, SIZE, MPI_REAL, source_id0A, rtag,  &
!                                   comm_cart, status, errcode)
        CALL MPI_ISEND           (BUFXSF, SIZE, MPI_REAL,              &
                                   dest_id0A, stag,                    &
                               comm_cart, request, errcode)
        CALL MPI_RECV            (BUFXRF, SIZE, MPI_REAL,              &
                                                    source_id0A, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE (X_MIN_DF:X_MIN_DF+OVERLAP-1,Y_MINF:Y_MAXF,Z_MINF:Z_MAXF) = &
          BUFXRF(1:OVERLAP,1:Y_MAXF-Y_MINF+1,1:Z_MAXF-Z_MINF+1)
        CALL MPI_WAIT (request, status, errcode)
      END IF
    ELSE
      IF ( source_id0A /= MPI_PROC_NULL ) THEN
        CALL MPI_RECV            (BUFXRF, SIZE, MPI_REAL,              &
                                                    source_id0A, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE (X_MIN_DF:X_MIN_DF+OVERLAP-1,Y_MINF:Y_MAXF,Z_MINF:Z_MAXF) = &
          BUFXRF(1:OVERLAP,1:Y_MAXF-Y_MINF+1,1:Z_MAXF-Z_MINF+1)
      END IF
    END IF


    stag = 1
    rtag = 1
    IF ( dest_id0B /= MPI_PROC_NULL ) THEN
      BUFXSF = VARIABLE (X_MINF:X_MINF+OVERLAP-1,Y_MINF:Y_MAXF,Z_MINF:Z_MAXF)
      IF ( source_id0B == MPI_PROC_NULL ) THEN
        CALL MPI_SEND            (BUFXSF, SIZE, MPI_REAL,              &
                                   dest_id0B, stag,                    &
                               comm_cart, errcode)
      ELSE
!        CALL MPI_SENDRECV (BUFXSF, SIZE, MPI_REAL,   dest_id0B, stag,  &
!                           BUFXRF, SIZE, MPI_REAL, source_id0B, rtag,  &
!                                   comm_cart, status, errcode)
        CALL MPI_ISEND           (BUFXSF, SIZE, MPI_REAL,              &
                                   dest_id0B, stag,                    &
                               comm_cart, request, errcode)
        CALL MPI_RECV            (BUFXRF, SIZE, MPI_REAL,              &
                                                    source_id0B, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE (X_MAX_DF+1-OVERLAP:X_MAX_DF,Y_MINF:Y_MAXF,Z_MINF:Z_MAXF) = &
          BUFXRF
        CALL MPI_WAIT (request, status, errcode)
      END IF
    ELSE
      IF ( source_id0B /= MPI_PROC_NULL ) THEN
        CALL MPI_RECV            (BUFXRF, SIZE, MPI_REAL,              &
                                                    source_id0B, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE (X_MAX_DF+1-OVERLAP:X_MAX_DF,Y_MINF:Y_MAXF,Z_MINF:Z_MAXF) = &
          BUFXRF
      END IF
    END IF

  END IF


! --- communicate in y direction

  SIZE = (X_MAX_DF-X_MIN_DF+1)*OVERLAP*(Z_MAXF-Z_MINF+1)

  IF ( SIZE > 0 )  THEN

    stag = 2
    rtag = 2
    IF ( dest_id1A /= MPI_PROC_NULL ) THEN
      BUFYSF = VARIABLE (X_MIN_DF:X_MAX_DF,Y_MAXF+1-OVERLAP:Y_MAXF,Z_MINF:Z_MAXF)
      IF ( source_id1A == MPI_PROC_NULL ) THEN
        CALL MPI_SEND            (BUFYSF, SIZE, MPI_REAL,               &
                                   dest_id1A, stag,                    &
                               comm_cart, errcode)
      ELSE
!        CALL MPI_SENDRECV (BUFYSF, SIZE, MPI_REAL,   dest_id1A, stag,  &
!                           BUFYRF, SIZE, MPI_REAL, source_id1A, rtag,  &
!                                   comm_cart, status, errcode)
        CALL MPI_ISEND           (BUFYSF, SIZE, MPI_REAL,               &
                                   dest_id1A, stag,                    &
                               comm_cart, request, errcode)
        CALL MPI_RECV            (BUFYRF, SIZE, MPI_REAL,              &
                                                    source_id1A, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE(X_MIN_DF:X_MAX_DF,Y_MIN_DF:Y_MIN_DF+OVERLAP-1,Z_MINF:Z_MAXF)&
          = BUFYRF
        CALL MPI_WAIT (request, status, errcode)
      END IF
    ELSE
      IF ( source_id1A /= MPI_PROC_NULL ) THEN
        CALL MPI_RECV            (BUFYRF, SIZE, MPI_REAL,              &
                                                    source_id1A, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE(X_MIN_DF:X_MAX_DF,Y_MIN_DF:Y_MIN_DF+OVERLAP-1,Z_MINF:Z_MAXF)&
          = BUFYRF
      END IF
    END IF


    stag = 3
    rtag = 3
    IF ( dest_id1B /= MPI_PROC_NULL ) THEN
      BUFYSF = VARIABLE (X_MIN_DF:X_MAX_DF,Y_MINF:Y_MINF+OVERLAP-1,Z_MINF:Z_MAXF)
      IF ( source_id1B == MPI_PROC_NULL ) THEN
        CALL MPI_SEND            (BUFYSF, SIZE, MPI_REAL,               &
                                   dest_id1B, stag,                    &
                               comm_cart, errcode)
      ELSE
!        CALL MPI_SENDRECV (BUFYSF, SIZE, MPI_REAL,   dest_id1B, stag,  &
!                           BUFYRF, SIZE, MPI_REAL, source_id1B, rtag,  &
!                                   comm_cart, status, errcode)
        CALL MPI_ISEND           (BUFYSF, SIZE, MPI_REAL,               &
                                   dest_id1B, stag,                    &
                               comm_cart, request, errcode)
        CALL MPI_RECV            (BUFYRF, SIZE, MPI_REAL,              &
                                                    source_id1B, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE(X_MIN_DF:X_MAX_DF,Y_MAX_DF+1-OVERLAP:Y_MAX_DF,Z_MINF:Z_MAXF)&
          = BUFYRF
        CALL MPI_WAIT (request, status, errcode)
      END IF
    ELSE
      IF ( source_id1B /= MPI_PROC_NULL ) THEN
        CALL MPI_RECV            (BUFYRF, SIZE, MPI_REAL,              &
                                                    source_id1B, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE(X_MIN_DF:X_MAX_DF,Y_MAX_DF+1-OVERLAP:Y_MAX_DF,Z_MINF:Z_MAXF)&
          = BUFYRF
      END IF
    END IF

  END IF

! --- communicate in z direction

  SIZE = (X_MAX_DF-X_MIN_DF+1)*(Y_MAX_DF-Y_MIN_DF+1)*OVERLAP

  IF ( Z_MINF <= Z_MAXF ) THEN

    stag = 4
    rtag = 4
    IF ( dest_id2AF /= MPI_PROC_NULL ) THEN
      BUFZSF = VARIABLE                                                &
                (X_MIN_DF:X_MAX_DF,Y_MIN_DF:Y_MAX_DF,Z_MAXF+1-OVERLAP:Z_MAXF)
      IF ( source_id2AF == MPI_PROC_NULL ) THEN
        CALL MPI_SEND            (BUFZSF, SIZE, MPI_REAL,              &
                                   dest_id2AF, stag,                   &
                               comm_cart, errcode)
      ELSE
!        CALL MPI_SENDRECV (BUFZSF, SIZE, MPI_REAL,   dest_id2AF, stag, &
!                           BUFZRF, SIZE, MPI_REAL, source_id2AF, rtag, &
!                                   comm_cart, status, errcode)
        CALL MPI_ISEND           (BUFZSF, SIZE, MPI_REAL,              &
                                   dest_id2AF, stag,                   &
                               comm_cart, request, errcode)
        CALL MPI_RECV            (BUFZRF, SIZE, MPI_REAL,              &
                                                   source_id2AF, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE                                                       &
       (X_MIN_DF:X_MAX_DF,Y_MIN_DF:Y_MAX_DF,Z_MIN_DF:Z_MIN_DF+OVERLAP-1)&
       = BUFZRF
        CALL MPI_WAIT (request, status, errcode)
      END IF
    ELSE
      IF ( source_id2AF /= MPI_PROC_NULL ) THEN
        CALL MPI_RECV            (BUFZRF, SIZE, MPI_REAL,              &
                                                    source_id2AF, rtag,&
                                   comm_cart, status, errcode)
        VARIABLE                                                       &
         (X_MIN_DF:X_MAX_DF,Y_MIN_DF:Y_MAX_DF,Z_MIN_DF:Z_MIN_DF+OVERLAP-1) = &
         BUFZRF
      END IF
    END IF


    stag = 5
    rtag = 5
    IF ( dest_id2BF /= MPI_PROC_NULL ) THEN
      BUFZSF = VARIABLE                                                &
                (X_MIN_DF:X_MAX_DF,Y_MIN_DF:Y_MAX_DF,Z_MINF:Z_MINF+OVERLAP-1)
      IF ( source_id2BF == MPI_PROC_NULL ) THEN
        CALL MPI_SEND            (BUFZSF, SIZE, MPI_REAL,              &
                                   dest_id2BF, stag,                   &
                               comm_cart, errcode)
      ELSE
!        CALL MPI_SENDRECV (BUFZSF, SIZE, MPI_REAL,   dest_id2BF, stag, &
!                           BUFZRF, SIZE, MPI_REAL, source_id2BF, rtag, &
!                                   comm_cart, status, errcode)
        CALL MPI_ISEND           (BUFZSF, SIZE, MPI_REAL,              &
                                   dest_id2BF, stag,                   &
                               comm_cart, request, errcode)
        CALL MPI_RECV            (BUFZRF, SIZE, MPI_REAL,              &
                                                   source_id2BF, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE                                                       &
           (X_MIN_DF:X_MAX_DF,Y_MIN_DF:Y_MAX_DF,Z_MAX_DF+1-OVERLAP:Z_MAX_DF) &
           = BUFZRF
        CALL MPI_WAIT (request, status, errcode)
      END IF
    ELSE
      IF ( source_id2BF /= MPI_PROC_NULL ) THEN
        CALL MPI_RECV            (BUFZRF, SIZE, MPI_REAL,              &
                                                   source_id2BF, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE                                                       &
           (X_MIN_DF:X_MAX_DF,Y_MIN_DF:Y_MAX_DF,Z_MAX_DF+1-OVERLAP:Z_MAX_DF) &
           = BUFZRF
      END IF
    END IF

  END IF
#endif


END SUBROUTINE COMMUNICATE_F

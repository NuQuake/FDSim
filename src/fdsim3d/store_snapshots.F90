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

SUBROUTINE  STORE_SNAPSHOTS(ITIST)

  USE PRECISION   , ONLY: PP

  USE CONTROL_DATA, ONLY: MZ, LPAS, STRESS_IMAGING

  USE AUXIL       , ONLY: X_MIN, X_MAX, Y_MIN, Y_MAX, Z_MIN, Z_MAX,    &
                          LOGSNP, BUFS

  USE WAVEFIELD   , ONLY: UM, VM, WM, SNAPU, SNAPV, SNAPW
#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: ITIST

  CHARACTER (LEN = 14) :: TLBUF

  INTEGER :: J, JJ, I, K, L, SIZE

!-----------------------------------------------------------------------

#if(USE_MPI)
  DO L = LPAS, MZ

    IF ( LOGSNP(L) ) THEN

      IF ( myid == root ) THEN

        DO J = 0, numprocs-1

          CALL MPI_CART_COORDS (comm_cart, J, ndims, coord, errcode)

          IF ( ( Z_MIN_A(coord(3)) <= L ) .AND.                        &
                                    ( L <= Z_MAX_A(coord(3)) ) ) THEN

            SIZE = ( X_MAX_A(coord(1)) - X_MIN_A(coord(1)) + 1 )       &
                 * ( Y_MAX_A(coord(2)) - Y_MIN_A(coord(2)) + 1 )

            IF ( J /= root ) THEN

              tag  = 9

              CALL MPI_RECV (BUFS , SIZE, MPI_REAL, J, tag, comm_cart,   &
                                                       status, errcode)
              JJ=1
              DO K = Y_MIN_A (coord(2)), Y_MAX_A (coord(2))
                DO I = X_MIN_A (coord(1)), X_MAX_A (coord(1))
                  SNAPU  ( I, K )=BUFS (JJ)
                  JJ = JJ + 1
                END DO
              END DO

              CALL MPI_RECV (BUFS , SIZE, MPI_REAL, J, tag, comm_cart,   &
                                                       status, errcode)
              JJ=1
              DO K = Y_MIN_A (coord(2)), Y_MAX_A (coord(2))
                DO I = X_MIN_A (coord(1)), X_MAX_A (coord(1))
                  SNAPV  ( I, K )=BUFS (JJ)
                  JJ = JJ + 1
                END DO
              END DO

              CALL MPI_RECV (BUFS , SIZE, MPI_REAL, J, tag, comm_cart,   &
                                                       status, errcode)
              JJ=1
              DO K = Y_MIN_A (coord(2)), Y_MAX_A (coord(2))
                DO I = X_MIN_A (coord(1)), X_MAX_A (coord(1))
                  SNAPW  ( I, K )=BUFS (JJ)
                  JJ = JJ + 1
                END DO
              END DO

            ELSE
              SNAPU( X_MIN : X_MAX, Y_MIN : Y_MAX ) =                  &
                UM ( X_MIN : X_MAX, Y_MIN : Y_MAX, L )
              SNAPV( X_MIN : X_MAX, Y_MIN : Y_MAX ) =                  &
                VM ( X_MIN : X_MAX, Y_MIN : Y_MAX, L )
              SNAPW( X_MIN : X_MAX, Y_MIN : Y_MAX ) =                  &
                WM ( X_MIN : X_MAX, Y_MIN : Y_MAX, L )
            END IF
          END IF

        END DO

        IF (STRESS_IMAGING) THEN
          WRITE ( TLBUF,'(A4, I3.3, A2, I5.5)' )  'SNAP', L-2, '.V', ITIST
        ELSE
          WRITE ( TLBUF,'(A4, I3.3, A2, I5.5)' )  'SNAP', L  , '.V', ITIST
        ENDIF
        OPEN  ( 18, FILE   = TLBUF(1:14),  FORM   ='UNFORMATTED',      &
                STATUS = 'REPLACE')

        WRITE ( 18 ) SNAPU, SNAPV, -SNAPW
        CLOSE ( 18 )

      ELSE

        IF ( ( Z_MIN  <= L ) .AND. ( L <= Z_MAX ) ) THEN

          SIZE = ( X_MAX  - X_MIN  + 1 ) * ( Y_MAX  - Y_MIN  + 1 )

          tag=9

          JJ=1
          DO K = Y_MIN , Y_MAX
            DO I = X_MIN , X_MAX
              BUFS (JJ) = UM  ( I, K, L )
              JJ = JJ + 1
            END DO
          END DO
          CALL MPI_SEND(BUFS ,SIZE,MPI_REAL,root,tag,comm_cart,errcode)

          JJ=1
          DO K = Y_MIN , Y_MAX
            DO I = X_MIN , X_MAX
              BUFS (JJ) = VM  ( I, K, L )
              JJ = JJ + 1
            END DO
          END DO
          CALL MPI_SEND(BUFS ,SIZE,MPI_REAL,root,tag,comm_cart,errcode)

          JJ=1
          DO K = Y_MIN , Y_MAX
            DO I = X_MIN , X_MAX
              BUFS (JJ) = WM  ( I, K, L )
              JJ = JJ + 1
            END DO
          END DO
          CALL MPI_SEND(BUFS ,SIZE,MPI_REAL,root,tag,comm_cart,errcode)

        END IF

      END IF

    END IF

  END DO

#else
  DO L = LPAS, MZ
    IF ( LOGSNP(L) ) THEN
      IF(STRESS_IMAGING) THEN
        WRITE ( TLBUF,'(A4, I3.3, A2, I5.5)' )  'SNAP', L-2, '.V', ITIST
      ELSE
        WRITE ( TLBUF,'(A4, I3.3, A2, I5.5)' )  'SNAP', L  , '.V', ITIST
      END IF
      OPEN  ( 18, FILE   = TLBUF(1:14),  FORM   ='UNFORMATTED',        &
                  STATUS = 'REPLACE')

      WRITE ( 18 ) UM (: ,: ,L), VM (: ,: ,L), -WM (: ,: ,L)

      CLOSE ( 18 )
    END IF
  END DO
#endif


END SUBROUTINE  STORE_SNAPSHOTS

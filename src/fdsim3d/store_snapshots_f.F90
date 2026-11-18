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

SUBROUTINE  STORE_SNAPSHOTS_F(ITIST)

  USE PRECISION   , ONLY: PP

  USE CONTROL_DATA, ONLY: LPAS, STRESS_IMAGING

  USE AUXIL       , ONLY: X_MINF,X_MAXF, Y_MINF,Y_MAXF, Z_MINF,Z_MAXF, &
                          LOGSNP, BUFSF

  USE WAVEFIELD   , ONLY: UMF, VMF, WMF, SNAPUF, SNAPVF, SNAPWF
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
  DO L = 0, LPAS-1

    IF ( LOGSNP(L) ) THEN

      IF ( myid == root ) THEN

        DO J = 0, numprocs-1

          CALL MPI_CART_COORDS (comm_cart, J, ndims, coord, errcode)

          IF ( ( Z_MIN_AF(coord(3)) <= L ) .AND.                       &
                                     ( L <= Z_MAX_AF(coord(3)) ) ) THEN

            SIZE = ( X_MAX_AF(coord(1)) - X_MIN_AF(coord(1)) + 1 )     &
                 * ( Y_MAX_AF(coord(2)) - Y_MIN_AF(coord(2)) + 1 )

            IF ( J /= root ) THEN

              tag  = 9

              CALL MPI_RECV (BUFSF, SIZE, MPI_REAL, J, tag, comm_cart,   &
                                                       status, errcode)
              JJ=1
              DO K = Y_MIN_AF(coord(2)), Y_MAX_AF(coord(2))
                DO I = X_MIN_AF(coord(1)), X_MAX_AF(coord(1))
                  SNAPUF ( I, K )=BUFSF(JJ)
                  JJ = JJ + 1
                END DO
              END DO

              CALL MPI_RECV (BUFSF, SIZE, MPI_REAL, J, tag, comm_cart,   &
                                                       status, errcode)
              JJ=1
              DO K = Y_MIN_AF(coord(2)), Y_MAX_AF(coord(2))
                DO I = X_MIN_AF(coord(1)), X_MAX_AF(coord(1))
                  SNAPVF ( I, K )=BUFSF(JJ)
                  JJ = JJ + 1
                END DO
              END DO

              CALL MPI_RECV (BUFSF, SIZE, MPI_REAL, J, tag, comm_cart,   &
                                                       status, errcode)
              JJ=1
              DO K = Y_MIN_AF(coord(2)), Y_MAX_AF(coord(2))
                DO I = X_MIN_AF(coord(1)), X_MAX_AF(coord(1))
                  SNAPWF ( I, K )=BUFSF(JJ)
                  JJ = JJ + 1
                END DO
              END DO

            ELSE
              SNAPUF( X_MINF : X_MAXF, Y_MINF : Y_MAXF ) =             &
                UMF ( X_MINF : X_MAXF, Y_MINF : Y_MAXF, L )
              SNAPVF( X_MINF : X_MAXF, Y_MINF : Y_MAXF ) =             &
                VMF ( X_MINF : X_MAXF, Y_MINF : Y_MAXF, L )
              SNAPWF( X_MINF : X_MAXF, Y_MINF : Y_MAXF ) =             &
                WMF ( X_MINF : X_MAXF, Y_MINF : Y_MAXF, L )
            END IF
          END IF

        END DO

        IF (STRESS_IMAGING) THEN
          WRITE ( TLBUF,'(A4, I3.3, A2, I5.5)' )  'SNAP', L-2, '.V', ITIST
        ELSE
          WRITE ( TLBUF,'(A4, I3.3, A2, I5.5)' )  'SNAP', L  , '.V', ITIST
        END IF
        OPEN  ( 18, FILE   = TLBUF(1:14),  FORM   ='UNFORMATTED',      &
                STATUS = 'REPLACE')

        WRITE ( 18 ) SNAPUF, SNAPVF, -SNAPWF
        CLOSE ( 18 )

      ELSE

        IF ( ( Z_MINF <= L ) .AND. ( L <= Z_MAXF) ) THEN

          SIZE = ( X_MAXF - X_MINF + 1 ) * ( Y_MAXF - Y_MINF + 1 )

          tag=9

          JJ=1
          DO K = Y_MINF, Y_MAXF
            DO I = X_MINF, X_MAXF
              BUFSF(JJ) = UMF ( I, K, L )
              JJ = JJ + 1
            END DO
          END DO
          CALL MPI_SEND(BUFSF,SIZE,MPI_REAL,root,tag,comm_cart,errcode)

          JJ=1
          DO K = Y_MINF, Y_MAXF
            DO I = X_MINF, X_MAXF
              BUFSF(JJ) = VMF ( I, K, L )
              JJ = JJ + 1
            END DO
          END DO
          CALL MPI_SEND(BUFSF,SIZE,MPI_REAL,root,tag,comm_cart,errcode)

          JJ=1
          DO K = Y_MINF, Y_MAXF
            DO I = X_MINF, X_MAXF
              BUFSF(JJ) = WMF ( I, K, L )
              JJ = JJ + 1
            END DO
          END DO
          CALL MPI_SEND(BUFSF,SIZE,MPI_REAL,root,tag,comm_cart,errcode)

        END IF

      END IF

    END IF

  END DO

#else
  DO L = 0, LPAS-1
    IF ( LOGSNP(L) ) THEN
      IF (STRESS_IMAGING) THEN
        WRITE ( TLBUF,'(A4, I3.3, A2, I5.5)' )  'SNAP', L-2, '.V', ITIST
      ELSE
        WRITE ( TLBUF,'(A4, I3.3, A2, I5.5)' )  'SNAP', L  , '.V', ITIST
      END IF
      OPEN  ( 18, FILE   = TLBUF(1:14),  FORM   ='UNFORMATTED',        &
                  STATUS = 'REPLACE')

      WRITE ( 18 ) UMF (: ,: ,L), VMF (: ,: ,L), -WMF (: ,: ,L)

      CLOSE ( 18 )
    END IF
  END DO
#endif


END SUBROUTINE  STORE_SNAPSHOTS_F

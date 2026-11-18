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
SUBROUTINE WRITE_PARAM

  USE NRTYPE              , ONLY: SP
  USE MOD_CONST           , ONLY: N_FREQ, N_E_PAR, N_Y_PAR
  USE MOD_PARAMETERIZATION, ONLY: N_CELLS, JMH, JMH3, CELLS, CELLSY,   &
                                  N_CELLS_G, JMH_G, JMH3_G, CELLS_G,   &
                                  CELLSY_G
  USE MOD_MODEL           , ONLY: MX, MY, MXF, MYF, LPAS, MZ, KEY_Q,   &
                                  KEY_QK_INF, MX_G, MXF_G
  USE MOD_FILES           , ONLY: JMH_FILE_NAME, JMHF_FILE_NAME,       &
                                  MO_FILE_NAME, Q_FILE_NAME
#if(USE_MPI)
  USE MPI_SHARED
#endif

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: JM1, L, TMP, NC, I

  INTEGER , ALLOCATABLE :: BUF_I (:,:,:)

  REAL(SP), ALLOCATABLE :: BUF(:,:), BUFY(:,:,:)

!----------------------------------------------------------------------

#if(USE_MPI)
  call MPI_BARRIER ( comm, errcode )

  IF ( myid /= root ) THEN

    tag = 1
    CALL MPI_SEND ( N_CELLS, 1, MPI_INTEGER, root, tag, comm, errcode )

    ALLOCATE ( BUF (N_CELLS,N_E_PAR) )
    BUF = CELLS(1:N_CELLS,1:N_E_PAR)

    tag = 2
    CALL MPI_SEND ( BUF, N_CELLS*N_E_PAR, MPI_REAL, root, tag, comm, errcode )
    DEALLOCATE ( BUF )

    IF (KEY_Q) THEN
      ALLOCATE ( BUFY (N_CELLS,N_Y_PAR,N_FREQ) )
      BUFY = CELLSY(1:N_CELLS,1:N_Y_PAR,1:N_FREQ)
      tag = 22
      CALL MPI_SEND ( BUFY, N_CELLS*N_Y_PAR*N_FREQ, MPI_REAL, root, tag, comm, errcode )
      DEALLOCATE ( BUFY )
    END IF

    tag = 3
    CALL MPI_SEND ( JMH, MX*MY*(MZ-LPAS+1), MPI_INTEGER, root, tag,    &
                    comm, errcode )
    tag = 4
    CALL MPI_SEND ( JMH3, MXF*MYF*(LPAS+1), MPI_INTEGER, root, tag,    &
                    comm, errcode )

  ELSE

    ALLOCATE ( JMH_G  (1:MY , 1:MX_G , LPAS:MZ   ),                    &
               JMH3_G (1:MYF, 1:MXF_G, 0   :LPAS )  )
    ALLOCATE ( N_CELLS_A (0:numprocs-1) )

    DO PROC = 0, numprocs-1
      IF ( PROC /= root ) THEN
        tag = 1
        CALL MPI_RECV ( TMP, 1, MPI_INTEGER, PROC, tag, comm, status, errcode )
        N_CELLS_A (PROC) = TMP
      ELSE
        N_CELLS_A (PROC) = N_CELLS
      END IF
    END DO

    N_CELLS_G = SUM(N_CELLS_A)

    ALLOCATE ( CELLS_G(1:N_CELLS_G,1:N_E_PAR) )
    JM1 = 1
    DO PROC = 0, numprocs-1
      IF ( PROC /= root ) THEN
        tag = 2
        ALLOCATE ( BUF (N_CELLS_A(PROC),N_E_PAR) )
        CALL MPI_RECV ( BUF, N_CELLS_A(PROC)*N_E_PAR, MPI_REAL, PROC, tag,   &
                        comm, status, errcode )
        CELLS_G ( JM1:JM1+N_CELLS_A(PROC)-1, 1:N_E_PAR ) = BUF
        DEALLOCATE ( BUF )
      ELSE
        CELLS_G ( JM1:JM1+N_CELLS_A(PROC)-1, 1:N_E_PAR ) = CELLS(1:N_CELLS,1:N_E_PAR)
      END IF
      JM1 = JM1 + N_CELLS_A(PROC)
    END DO

    ALLOCATE ( CELLSY_G(1:N_CELLS_G,1:N_Y_PAR,1:N_FREQ) )
    IF (KEY_Q) THEN
      JM1 = 1
      DO PROC = 0, numprocs-1
        IF ( PROC /= root ) THEN
          tag = 22
          ALLOCATE ( BUFY (N_CELLS_A(PROC),N_Y_PAR,N_FREQ) )
          CALL MPI_RECV ( BUFY, N_CELLS_A(PROC)*N_Y_PAR*N_FREQ, MPI_REAL, PROC, tag,   &
                          comm, status, errcode )
          CELLSY_G ( JM1:JM1+N_CELLS_A(PROC)-1, 1:N_Y_PAR, 1:N_FREQ ) = BUFY
          DEALLOCATE ( BUFY )
        ELSE
          CELLSY_G ( JM1:JM1+N_CELLS_A(PROC)-1, 1:N_Y_PAR, 1:N_FREQ ) = CELLSY(1:N_CELLS,:,:)
        END IF
        JM1 = JM1 + N_CELLS_A(PROC)
      END DO
    END IF

    JM1 = 1
    NC  = 0
    DO PROC = 0, numprocs-1
      IF ( PROC /= root ) THEN
        tag = 3
        ALLOCATE ( BUF_I (MY,MX_A(PROC),LPAS:MZ ) )
        CALL MPI_RECV ( BUF_I, MX_A(PROC)*MY*(MZ-LPAS+1), MPI_INTEGER, &
                        PROC, tag, comm, status, errcode )
        JMH_G ( :,JM1:JM1+MX_A(PROC)-1,: ) = BUF_I + NC
        DEALLOCATE ( BUF_I )
      ELSE
        JMH_G ( :,JM1:JM1+MX_A(PROC)-1,: ) = JMH (:,1:MX,:) + NC
      END IF
      JM1 = JM1 + MX_A(PROC) - 1
      NC  = NC  + N_CELLS_A(PROC)
    END DO

    JM1 = 1
    NC  = 0
    DO PROC = 0, numprocs-1
      IF ( PROC /= root ) THEN
        tag = 4
        ALLOCATE ( BUF_I (MYF,MXF_A(PROC),0:LPAS) )
        CALL MPI_RECV ( BUF_I, MXF_A(PROC)*MYF*(LPAS+1), MPI_INTEGER,  &
                        PROC, tag, comm, status, errcode )
              JMH3_G ( : ,JM1:JM1+MXF_A(PROC)-1, : ) = BUF_I + NC
        DEALLOCATE ( BUF_I )
      ELSE
        JMH3_G ( :,JM1:JM1+MXF_A(PROC)-1,: ) = JMH3 (:,1:MXF,:) + NC
      END IF

      JM1 = JM1 + MXF_A(PROC) - 1
      NC  = NC  + N_CELLS_A(PROC)
    END DO

#else
  N_CELLS_G = N_CELLS
#endif

    OPEN ( 13, FILE   = JMH_FILE_NAME, FORM   = 'UNFORMATTED',         &
                                       STATUS ='NEW' )

    DO L = LPAS, MZ
      DO I = 1, MX_G
#if(USE_MPI)
        WRITE ( 13 )  JMH_G ( 1:MY, I, L)
#else
        WRITE ( 13 )  JMH   ( 1:MY, I, L)
#endif
      END DO
    END DO

    CLOSE( 13 )

    IF ( LPAS > 0 ) THEN

      OPEN ( 33, FILE   = JMHF_FILE_NAME, FORM   = 'UNFORMATTED',      &
                                          STATUS ='NEW' )

      DO L = 0, LPAS
        DO I = 1, MXF_G
#if(USE_MPI)
          WRITE( 33 )  JMH3_G( 1:MYF, I, L)
#else
          WRITE( 33 )  JMH3  ( 1:MYF, I, L)
#endif
        END DO
      END DO

      CLOSE( 33 )

    END IF

    OPEN ( 14, FILE = MO_FILE_NAME, FORM = 'UNFORMATTED', STATUS='NEW' )

    OPEN ( 15, FILE =  Q_FILE_NAME, FORM = 'UNFORMATTED', STATUS='NEW' )


    WRITE( 14 ) N_CELLS_G

#if(USE_MPI)
    WRITE( 14 ) ( CELLS_G(JM1,10), CELLS_G(JM1,11), CELLS_G(JM1,12),   &
                  CELLS_G(JM1, 7), CELLS_G(JM1, 8), CELLS_G(JM1, 9),   &
                  CELLS_G(JM1, 6), CELLS_G(JM1, 4), CELLS_G(JM1, 5),   &
                  CELLS_G(JM1, 3), CELLS_G(JM1, 1), CELLS_G(JM1, 2),   &
                                                    JM1 = 1, N_CELLS_G )
#else
    WRITE( 14 ) ( CELLS  (JM1,10), CELLS  (JM1,11), CELLS  (JM1,12),   &
                  CELLS  (JM1, 7), CELLS  (JM1, 8), CELLS  (JM1, 9),   &
                  CELLS  (JM1, 6), CELLS  (JM1, 4), CELLS  (JM1, 5),   &
                  CELLS  (JM1, 3), CELLS  (JM1, 1), CELLS  (JM1, 2),   &
                                                    JM1 = 1, N_CELLS   )
#endif

    CLOSE(14)

    IF ( KEY_Q ) THEN

      IF ( KEY_QK_INF ) THEN

        DO L = 1, N_CELLS
#if(USE_MPI)
          CELLSY_G(L,7,:) =  4./3. * CELLSY_G(L,2,:) * CELLS_G(L,2)/CELLS_G(L,7)
          CELLSY_G(L,8,:) =  4./3. * CELLSY_G(L,1,:) * CELLS_G(L,1)/CELLS_G(L,8)
          CELLSY_G(L,9,:) =  4./3. * CELLSY_G(L,3,:) * CELLS_G(L,3)/CELLS_G(L,9)
          CELLSY_G(L,6,:) = -2./3. * CELLSY_G(L,3,:) * CELLS_G(L,3)/CELLS_G(L,6)
          CELLSY_G(L,4,:) = -2./3. * CELLSY_G(L,1,:) * CELLS_G(L,1)/CELLS_G(L,4)
          CELLSY_G(L,5,:) = -2./3. * CELLSY_G(L,2,:) * CELLS_G(L,2)/CELLS_G(L,5)
#else
          CELLSY  (L,7,:) =  4./3. * CELLSY  (L,2,:) * CELLS  (L,2)/CELLS  (L,7)
          CELLSY  (L,8,:) =  4./3. * CELLSY  (L,1,:) * CELLS  (L,1)/CELLS  (L,8)
          CELLSY  (L,9,:) =  4./3. * CELLSY  (L,3,:) * CELLS  (L,3)/CELLS  (L,9)
          CELLSY  (L,6,:) = -2./3. * CELLSY  (L,3,:) * CELLS  (L,3)/CELLS  (L,6)
          CELLSY  (L,4,:) = -2./3. * CELLSY  (L,1,:) * CELLS  (L,1)/CELLS  (L,4)
          CELLSY  (L,5,:) = -2./3. * CELLSY  (L,2,:) * CELLS  (L,2)/CELLS  (L,5)
#endif
        END DO
      END IF

    ELSE
#if(USE_MPI)
      CELLSY_G = 0.
    END IF

    WRITE(15) (CELLSY_G(JM1,7,:), CELLSY_G(JM1,8,:), CELLSY_G(JM1,9,:),&
               CELLSY_G(JM1,6,:), CELLSY_G(JM1,4,:), CELLSY_G(JM1,5,:),&
               CELLSY_G(JM1,3,:), CELLSY_G(JM1,1,:), CELLSY_G(JM1,2,:),&
                                                    JM1 = 1, N_CELLS_G )
#else
      CELLSY   = 0.
    END IF

    WRITE ( 15 ) ( CELLSY(JM1,7,:), CELLSY(JM1,8,:), CELLSY(JM1,9,:),  &
                   CELLSY(JM1,6,:), CELLSY(JM1,4,:), CELLSY(JM1,5,:),  &
                   CELLSY(JM1,3,:), CELLSY(JM1,1,:), CELLSY(JM1,2,:),  &
                                                      JM1 = 1, N_CELLS )
#endif

    CLOSE(15)

#if(USE_MPI)
  END IF

  CALL MPI_BARRIER( comm, errcode )

#endif

END SUBROUTINE WRITE_PARAM

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
!
! This module linearly interpolate NI input grids.
! Input grids are characterized by number of points
! in X and Y directions NPX and NPY.
! Steps in X and Y directions can variate.
!
! Example of GOOD input mesh:
!
! y ^
!   |  *--*--*--*--*
!   |  |           |
!   |  *--*---*--*-*
!   |  |           |
!   |  |           |
!   |  *--*--*-*---*
!   +----------------> x
!
! Example of BAD input mesh:
!
! y ^
!   |  *--*-*------*
!   |  |  | |      |
!   |  *  | *      *
!   |  |  * |      |
!   |  |  | |      |
!   |  *--*-*------*
!   +----------------> x
!
! Output is spatial grid with step H10.
! Whole output region have to be inside the input mesh.
!
!-----------------------------------------------------------------------

SUBROUTINE MODEL_LIN

  USE NRTYPE   , ONLY: WP, SP
  USE MOD_MODEL, ONLY: NPX, NPY, NX10, NY10, H, H10, NX10_G,           &
                       XIMIN, XIMAX, YIMIN, YIMAX, ITF1, NI,           &
                       XIMIN_G, XIMAX_G
  USE MOD_FILES, ONLY: FSED
#if(USE_MPI)
  USE MPI_SHARED
#endif

  !----------------------------------------------------------- variables

  IMPLICIT NONE

  REAL(SP), ALLOCATABLE :: ITF1_R(:,:), BUF_S (:,:), BUF(:,:)
  REAL(WP), ALLOCATABLE :: FILE_CONTENT(:,:), MESH(:,:,:), OUTY(:,:,:)

  INTEGER :: P, NP, I, J, K , EX, LI, IMIN, IMAX, I_OVERLAP
  REAL(WP) :: XD, YD, ZD, YLAST,X1, X2, Z1, Z2, Y1,Y2, XA, YA,          &
             XMIN, XMAX, YMIN, YMAX, MAXNUM

  !---------------------------------------------------- subroutines body

  ALLOCATE( ITF1( NX10, NY10, NI ), BUF (NX10, NY10) )

#if(USE_MPI)
  IF ( myid == root ) THEN
#endif

    ALLOCATE ( ITF1_R (NX10_G, NY10) )
    WRITE(*,'(A)')'Linear model of sediments/bedrock interface...'

    DO LI = 1, NI

    !----------------- read input file
      OPEN( 10, FILE=FSED(LI), STATUS='OLD')

      NP=0

    !--------------- count number of points in mesh
      DO
        READ(10,FMT=*,IOSTAT=EX)XD,YD,ZD
        IF (EX/=0 ) EXIT
        NP = NP+1
      END DO

      REWIND(10)

    !--------------- test if correct number of input planes was set
      IF (NPX*NPY/=NP) THEN
        WRITE(*,*)
        WRITE(*,'(A,A)')'NUMBER OF POINTS IN FILE ',FSED(LI)
        WRITE(*,'(A)')'DIFFER FROM NUMBER CALCULATED FROM NUMBERS OF &
                      &INPUT PLANES:'
        WRITE(*,'(A,I8,A,I8)')'POINTS IN FILE: ',NP,'NPX*NPY = ',NPX*NPY
        STOP
      END IF

      ALLOCATE( FILE_CONTENT(NP,3) )

      DO I = 1, NP
        READ(10,*) FILE_CONTENT(I,:)
      END DO

    !--------------- test if area that have to be aproximated lie inside
    !--------------- input mesh
      XMIN = MINVAL(FILE_CONTENT(:,1))
      XMAX = MAXVAL(FILE_CONTENT(:,1))
      YMIN = MINVAL(FILE_CONTENT(:,2))
      YMAX = MAXVAL(FILE_CONTENT(:,2))

      IF ( (XMIN>XIMIN_G).OR. &
           (XMAX<XIMAX_G).OR. &
           (YMIN>YIMIN  ).OR. &
           (YMAX<YIMAX  )     ) THEN
        WRITE(*,*)' '
        WRITE(*,'(A)')'AREA TO BE APROXIMATED DOES NOT LIE &
                      &INSIDE KNOWN INTERFACE BOUNDARIES.'
        WRITE(*,'(A,F23.10,A,F23.10,A)')'KNOWN INTERFACE:        BOX &
                &[ ',XMIN ,', ',YMIN ,']'
        WRITE(*,'(A,F23.10,A,F23.10,A)')'                            &
                &[ ',XMAX ,', ',YMAX ,']'
        WRITE(*,'(A,F23.10,A,F23.10,A)')'AREA TO BE APROXIMATED: BOX &
                &[ ',XIMIN_G,', ',YIMIN,']'
        WRITE(*,'(A,F23.10,A,F23.10,A)')'                            &
                &[ ',XIMAX_G,', ',YIMAX,']'
        STOP
      END IF

      CLOSE(10)

    !----------------- change format of input data
      ALLOCATE(MESH(NPX,NPY,3))

      MESH  = 1
      YLAST = FILE_CONTENT(1,2)
      MAXNUM = MAXVAL(FILE_CONTENT(:,2))+10.

      DO J = 1, NPY
        YLAST = MINVAL(FILE_CONTENT(:,2))
        I     = 1
        DO K = 1, NP

          IF ( FILE_CONTENT(K,2)==YLAST ) THEN
            MESH(I,J,:) = FILE_CONTENT(K,:)
            FILE_CONTENT(K,2) = MAXNUM
            I=I+1
          END IF

        END DO

        IF ( (I-1) /= NPX ) THEN
          PRINT *, 'ERROR IN INPUT DATA !!! '
          STOP
        END IF
      END DO

      DEALLOCATE(FILE_CONTENT)

    !============================================== aproximation of mesh

    !------------------ aproximating y slices
      ALLOCATE( OUTY(NX10_G,NPY,3) )

      DO J = 1, NPY

      !----- init
        XA = XIMIN_G
        P=1
        DO
          X1 = MESH(P  ,J,1)
          X2 = MESH(P+1,J,1)
          IF ((X1<=XA).AND.(XA<X2)) EXIT
          P = P + 1
        END DO

        Z1 = MESH(P  ,J,3)
        Z2 = MESH(P+1,J,3)

      !----- main cycle of aproximation of y slice
        DO I = 1, NX10_G
          XA = XIMIN_G + DBLE(I-1)*H10
          IF ( XA > X2 ) THEN
            P = P+1

            X1 = X2
            Z1 = Z2

            X2 = MESH(P+1,J,1)
            Z2 = MESH(P+1,J,3)

          END IF

          OUTY( I, J, 1 ) = XA
          OUTY( I, J, 2 ) = MESH(1,J,2)
          OUTY( I, J, 3 ) = ( (Z2-Z1)*XA + (Z1*X2-Z2*X1) ) / ( X2-X1 )

        END DO

      END DO

      DEALLOCATE(MESH)

    !------------------ aproximating x slices

      DO I = 1, NX10_G

      !----- init
        YA = YIMIN
        P=1
        DO
          Y1 = OUTY(I,P  ,2)
          Y2 = OUTY(I,P+1,2)
          IF ((Y1<=YA).AND.(YA<Y2)) EXIT
          P = P + 1
        END DO

        Z1 = OUTY( I, P  , 3)
        Z2 = OUTY( I, P+1, 3)

      !----- main cycle of aproximation of x slices
        DO J = 1, NY10
          YA = YIMIN + DBLE(J-1)*H10
          IF ( YA > Y2 ) THEN

            P = P+1

            Y1 = Y2
            Z1 = Z2

            Y2 = OUTY(I,P+1,2)
            Z2 = OUTY(I,P+1,3)

          END IF

          ITF1_R(I,J) = REAL( ( (Z2-Z1)*YA + (Z1*Y2-Z2*Y1) ) /         &
                              ( Y2-Y1 ), SP)

        END DO

      END DO

      DEALLOCATE(OUTY)

#if(USE_MPI)
      I_OVERLAP = 2*NINT(H/H10)+1
      IMAX = I_OVERLAP
      DO PROC = 0, numprocs-1
        IMIN = IMAX + 1 - I_OVERLAP
        IMAX = IMIN - 1 + NX10_A(PROC)
        SIZE = (IMAX - IMIN + 1)*NY10
        PRINT *, PROC, IMIN, IMAX, NX10_A(PROC), NX10_G
        IF ( PROC /= root ) THEN
          tag = 1
          ALLOCATE ( BUF_S(1:NX10_A(PROC),NY10) )
          BUF_S = ITF1_R(IMIN:IMAX,1:NY10)
          CALL MPI_SEND (BUF_S, SIZE, MPI_REAL, PROC, tag, comm, errcode)
          DEALLOCATE (BUF_S)
        ELSE
          ITF1 (1:NX10,1:NY10,LI) = ITF1_R (IMIN:IMAX,1:NY10)
        END IF
      END DO
#else
      ITF1 (:,:,LI) = ITF1_R (:,:)
#endif

      WRITE(*,'(A,I2,A)')'Interface #',LI,'OK'

    END DO

#if(USE_MPI)
  ELSE

    DO LI = 1, NI
      SIZE = NX10*NY10
      tag  = 1
      CALL MPI_RECV (BUF, SIZE, MPI_REAL, root, tag, comm, status, errcode)
      ITF1 (1:NX10,1:NY10,LI) = BUF(1:NX10,1:NY10)
    END DO

  END IF

#endif

END SUBROUTINE MODEL_LIN

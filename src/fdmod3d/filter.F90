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
! subroutine implements running avarage on aproximated sediments/bedrock
! interface with window width = RANGE
!-----------------------------------------------------------------------
SUBROUTINE FILTER

  USE NRTYPE   , ONLY: SP
  USE MOD_MODEL, ONLY: NX10, NY10, FLTR_RANGE, H10, ITF1, ni

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER  :: I, J, IM, JM, IMP, JMP, RACNT, RACNT_START, RACNT_STOPX, &
              RACNT_STOPY, NORM, k
  REAL(SP) :: SUM

  REAL(SP), DIMENSION(:,:), ALLOCATABLE :: OUT

!----------------------------------------------------------------------

  DO K = 1, NI
!    WRITE(*,'(A)')'Filtering interface(-s)...'

    ALLOCATE(OUT(NX10,NY10))

    OUT = 100. !bad value to indicate bad filtering

    RACNT = INT(FLTR_RANGE/2./H10+1.)

    RACNT_START = RACNT + 1
    RACNT_STOPX = NX10-RACNT
    RACNT_STOPY = NY10-RACNT

    ! Inner points
    NORM = (2*RACNT+1)*(2*RACNT+1)
    !(i,j) is position of point in which calculating running avarage
                                         !cycle through all inner points
    DO J = RACNT_START, RACNT_STOPY
      DO I = RACNT_START, RACNT_STOPX
        SUM = 0.
                                      !cycle through all adjacent points
        DO JM = -RACNT,RACNT
          DO IM = -RACNT,RACNT
            SUM = SUM + ITF1(I+IM,J+JM,k)
          END DO
        END DO

        OUT(I,J) = SUM/REAL(NORM,SP)

      END DO
    END DO

    ! Top margin (y around yimax)
    JMP = RACNT-1
    DO J = RACNT_STOPY+1, NY10
      NORM = (RACNT+1+JMP)*(2*RACNT+1)
      DO I = RACNT_START, RACNT_STOPX
        SUM  = 0.
        DO JM = -RACNT,JMP
          DO IM = -RACNT,RACNT
            SUM = SUM + ITF1(I+IM,J+JM,k)
          END DO
        END DO
        OUT(I,J) = SUM/REAL(NORM,SP)
      END DO
      JMP = JMP-1
    END DO

    ! Right margin (x around ximax)
    DO J = RACNT_START, RACNT_STOPY
      IMP = RACNT-1
      DO I = RACNT_STOPX+1, NX10
        NORM = (RACNT+1+IMP)*(2*RACNT+1)
        SUM  = 0.
        DO JM = -RACNT,RACNT
          DO IM = -RACNT,IMP
            SUM = SUM + ITF1(I+IM,J+JM,k)
          END DO
        END DO
        OUT(I,J) = SUM/REAL(NORM,SP)
        IMP = IMP-1
      END DO
    END DO

    ! Left margin (x around ximin)
    DO J = RACNT_START, RACNT_STOPY
      IMP = 0
      DO I = 1, RACNT
        NORM = (RACNT+1+IMP)*(2*RACNT+1)
        SUM  = 0.
        DO JM = -RACNT,RACNT
          DO IM = -IMP,RACNT
            SUM = SUM + ITF1(I+IM,J+JM,k)
          END DO
        END DO
        OUT(I,J) = SUM/REAL(NORM,SP)
        IMP = IMP+1
      END DO
    END DO

    ! Bottom margin (y around yimin)
    JMP = 0
    DO J = 1, RACNT
      NORM = (JMP+1+RACNT)*(2*RACNT+1)
      DO I = RACNT_START, RACNT_STOPX
        SUM  = 0.

        DO JM = -JMP,RACNT
          DO IM = -RACNT,RACNT
            SUM = SUM + ITF1(I+IM,J+JM,k)
          END DO
        END DO

        OUT(I,J) = SUM/REAL(NORM,SP)

      END DO
      JMP=JMP+1
    END DO

    ! Right-Top corner
    JMP = RACNT-1
    DO J = RACNT_STOPY+1, NY10
      IMP = RACNT-1
      DO I = RACNT_STOPX+1, NX10
        NORM = (RACNT+1+IMP)*(RACNT+1+JMP)
        SUM  = 0.
        DO JM = -RACNT, JMP
          DO IM = -RACNT, IMP
            SUM = SUM + ITF1(I+IM,J+JM,k)
          END DO
        END DO
        OUT(I,J) = SUM/REAL(NORM,SP)
        IMP = IMP - 1
      END DO
      JMP = JMP - 1
    END DO

    ! Left-Top corner
    JMP = RACNT-1
    DO J = RACNT_STOPY+1, NY10
      IMP = 0
      DO I = 1, RACNT
        NORM = (RACNT+1+IMP)*(RACNT+1+JMP)
        SUM  = 0.
        DO JM = -RACNT, JMP
          DO IM = -IMP, RACNT
            SUM = SUM + ITF1(I+IM,J+JM,k)
          END DO
        END DO
        OUT(I,J) = SUM/REAL(NORM,SP)
        IMP = IMP + 1
      END DO
      JMP = JMP - 1
    END DO

    ! Left-Bottom corner
    JMP = 0
    DO J = 1, RACNT
      IMP = 0
      DO I = 1, RACNT
        NORM = (RACNT+1+IMP)*(RACNT+1+JMP)
        SUM  = 0.
        DO JM = -JMP, RACNT
          DO IM = -IMP, RACNT
            SUM = SUM + ITF1(I+IM,J+JM,k)
          END DO
        END DO
        OUT(I,J) = SUM/REAL(NORM,SP)
        IMP = IMP + 1
      END DO
      JMP = JMP + 1
    END DO

    ! Right-Bottom corner
    JMP = 0
    DO J = 1, RACNT
      IMP = RACNT-1
      DO I = RACNT_STOPX+1, NX10
        NORM = (RACNT+1+IMP)*(RACNT+1+JMP)
        SUM  = 0.
        DO JM = -JMP, RACNT
          DO IM = -RACNT, IMP
            SUM = SUM + ITF1(I+IM,J+JM,k)
          END DO
        END DO
        OUT(I,J) = SUM/REAL(NORM,SP)
        IMP = IMP - 1
      END DO
      JMP = JMP + 1
    END DO

    ITF1(:,:,k) = OUT

    DEALLOCATE(OUT)

!    WRITE(*,'(A)')'OK'

  end do

END SUBROUTINE FILTER


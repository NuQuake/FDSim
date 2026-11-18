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
PROGRAM FDMod3D

USE MOD_MODEL     , ONLY: KEY_FILTER, KEY_SRFOUT, KEY_GBM3D
USE MOD_INTERFACES, ONLY: SETPARAMS, PRINT_INFO, WRITE_GRD,            &
                          MODEL_LIN, PARAMETERIZATION, FIND_MAX_DEPTH, &
                          WRITE_PARAM

#if(USE_MPI)
  USE MPI_SHARED
#endif

!-----------------------------------------------------------------------

IMPLICIT NONE

REAL :: TIME, TIMEB

!-----------------------------------------------------------------------

#if(USE_MPI)
  CALL MPI_INIT (errcode)
  comm = MPI_COMM_WORLD
  CALL MPI_COMM_RANK (comm, myid,     errcode)
  CALL MPI_Get_ProCessor_name   (nam_proc, len, errcode)
  CALL MPI_COMM_SIZE (comm, numprocs, errcode)
  WRITE(0,*) "Process ", myid, "of", numprocs, "is alive : ",nam_proc
  CALL MPI_Barrier   (comm,           errcode)
  root = 0
  ALLOCATE ( MX_A (0:numprocs-1), MXF_A (0:numprocs-1), NX10_A (0:numprocs-1) )
#endif

#if(USE_MPI)
  IF ( myid == root ) THEN
#endif
    CALL CPU_TIME (TIMEB)
#if(USE_MPI)
  END IF
#endif

CALL SETPARAMS

CALL MODEL_LIN

IF (KEY_GBM3D)  CALL READ_GBM3D

!----- filter aproximated interface

IF (KEY_FILTER) CALL FILTER

CALL FIND_MAX_DEPTH

CALL PARAMETERIZATION

CALL WRITE_PARAM

#if(USE_MPI)
IF ( myid == root ) THEN
#endif
CALL PRINT_INFO

CALL LOG_INFO

CALL CPU_TIME(TIME)

TIME = TIME-TIMEB
PRINT *,'CPU TIME: ',TIME,' [s]'
#if(USE_MPI)
END IF

CALL MPI_FINALIZE ( errcode )

#endif

END PROGRAM FDMod3D

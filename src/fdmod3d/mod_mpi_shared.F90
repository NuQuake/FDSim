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
!_______________________________________________________________________
MODULE MPI_SHARED

#if(USE_MPI)

#if(MPI2)
  USE MPI
#endif
!----------------------------------------------------------------------

  IMPLICIT NONE
#if(MPI2)
#else
  INCLUDE 'mpif.h'
#endif

!----------------------------------------------------------------------

  CHARACTER*40       :: nam_proc

  INTEGER            :: comm, errcode, myid, root,                     &
                        numprocs, numprocsplus,                        &
                        status(MPI_STATUS_SIZE)

  INTEGER            :: tag, len, PROC, SIZE

  INTEGER, ALLOCATABLE, DIMENSION(:) :: MX_A, MXF_A, NX10_A, N_CELLS_A
#endif

END MODULE MPI_SHARED

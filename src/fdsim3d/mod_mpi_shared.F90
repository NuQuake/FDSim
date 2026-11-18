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

MODULE MPI_SHARED

#if(USE_MPI)

#if(MPI2)
  USE MPI
#endif
!----------------------------------------------------------------------

  IMPLICIT NONE

!----------------------------------------------------------------------
#if(MPI2)
#else
  INCLUDE 'mpif.h'
#endif

  INTEGER, PARAMETER :: ndims = 3

  CHARACTER(LEN=40)  :: nam_proc

  LOGICAL        :: iamroot, periods(ndims)
  INTEGER        :: my_N
  INTEGER        :: comm, comm_cart, errcode, myid, root, numprocs,    &
                    status(MPI_STATUS_SIZE), request

  INTEGER        :: len, tag, OVERLAP

  INTEGER, DIMENSION(ndims)     :: dims, my_coord, coord

  INTEGER, ALLOCATABLE, DIMENSION(:) :: X_MIN_DA , X_MAX_DA ,          &
                                        Y_MIN_DA , Y_MAX_DA ,          &
                                        Z_MIN_DA , Z_MAX_DA ,          &
                                        X_MIN_A  , X_MAX_A  ,          &
                                        Y_MIN_A  , Y_MAX_A  ,          &
                                        Z_MIN_A  , Z_MAX_A  ,          &
                                        X_MIN_DAF, X_MAX_DAF,          &
                                        Y_MIN_DAF, Y_MAX_DAF,          &
                                        Z_MIN_DAF, Z_MAX_DAF,          &
                                        X_MIN_AF , X_MAX_AF ,          &
                                        Y_MIN_AF , Y_MAX_AF ,          &
                                        Z_MIN_AF , Z_MAX_AF


#endif

END MODULE MPI_SHARED

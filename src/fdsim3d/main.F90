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
! REVISION HISTORY (for details see CHANGELOG):
! Version 1.0, 2014-03-19: Jozef Kristek, Peter Moczo
!=======================================================================

PROGRAM FDSim3D

  USE PRECISION   , ONLY:                                              &
                          PP
  USE GRID_MEDIUM , ONLY:                                              &
                          H
  USE CONTROL_DATA, ONLY:                                              &
                          LPAS
  USE INTERFACES  , ONLY:                                              &
                          ALLOC_NR, ALLOC_T_SEIS, AUX1, DISTRIBUTE,    &
                          DOMAIN_DECOMPOSITION, OPEN_SR_FILE,          &
                          PPML, PPMLF, PRE_NR_COEF1, PRE_NR_COEF2,     &
                          READ_MODEL, READ_MODEL_PAR, READ_OPEN,       &
                          TIME_LOOP, LANCZOS
  USE TIME_FD     , ONLY:                                              &
                          time_begin, time_end, DATE_TIME, BIG_BEN

#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif

!-----------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER(LEN=17) :: JOB_NAME

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
#endif

#if(USE_MPI)
  IF ( myid == root ) THEN
#endif
    CALL CPU_TIME ( time_begin )
    CALL DATE_AND_TIME ( BIG_BEN(1), BIG_BEN(2), BIG_BEN(3), DATE_TIME )
    WRITE (*,'(''Current time: '',I2,'':'',I2,'':'',F6.3)')            &
      DATE_TIME(5), DATE_TIME(6), DATE_TIME(7)+REAL(DATE_TIME(8)/1000.)

    OPEN  (10, FILE = 'HF_FDSim3D', STATUS='OLD')
    READ  (10,'(A)') JOB_NAME
    CLOSE (10)

#if(USE_MPI)
  END IF
#endif

  CALL  READ_OPEN ( JOB_NAME (1:LEN_TRIM(JOB_NAME)) )
                         print *,  'read_open'
  CALL  DISTRIBUTE
                         print *,  'distribute'
  CALL  AUX1
                         print *,  'aux1'
  CALL  DOMAIN_DECOMPOSITION
                         print *,  'domain_decomposition'
  CALL  ALLOC_NR
                         print *,  'alloc_nr'
  CALL  ALLOC_T_SEIS
                         print *,  'alloc_t_seis'
  CALL  READ_MODEL
                         print *,  'read_model'
  CALL  READ_MODEL_PAR
                         print *,  'read_model_par'
  CALL  OPEN_SR_FILE
                         print *,  'open_sr_file'
  CALL  PRE_NR_COEF1
                         print *,  'pre_nr_coef1'
  IF ( LPAS /= 0 ) THEN
    CALL  PRE_NR_COEF2
                         print *,  'pre_nr_coef2'
    CALL  LANCZOS
                         print *,  'lanczos'
  END IF
  CALL  PPML
                         print *,  'ppml'
  IF ( LPAS /= 0 ) THEN
    CALL  PPMLF
                         print *,  'ppmlf'
  END IF
  CALL  TIME_LOOP
                         print *,  'time_loop'

#if(USE_MPI)
  IF ( myid == root ) THEN
#endif
    CALL CPU_TIME ( time_end )
    CALL DATE_AND_TIME ( BIG_BEN(1), BIG_BEN(2), BIG_BEN(3), DATE_TIME )
    PRINT     *, 'CPU time: ', time_end - time_begin, ' seconds'
    WRITE (*,'(''Current time: '',I2,'':'',I2,'':'',F6.3)')            &
      DATE_TIME(5), DATE_TIME(6), DATE_TIME(7)+REAL(DATE_TIME(8)/1000.)
    WRITE (11,'(F10.2,1X,I4,I3,I3,I5,3I3,I4)')                         &
             time_end - time_begin, DATE_TIME
#if(USE_MPI)
  END IF
  CALL MPI_Barrier   (comm,           errcode)
  CALL MPI_FINALIZE ( errcode )
#endif


END PROGRAM FDSim3D

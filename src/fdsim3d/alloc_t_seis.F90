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

SUBROUTINE ALLOC_T_SEIS

  USE PRECISION   , ONLY:                                              &
                          PP

  USE CONTROL_DATA, ONLY:                                              &
                          MR, LPAS

  USE AUXIL       , ONLY:                                              &
                          PMLL , PMLXH , PMLYH ,                       &
                          PMLLF, PMLXHF, PMLYHF,                       &
                          X_MIN_D , X_MAX_D ,                          &
                          Y_MIN_D , Y_MAX_D ,                          &
                          X_MIN_DF, X_MAX_DF,                          &
                          Y_MIN_DF, Y_MAX_DF

  USE STRAINFIELD , ONLY:                                              &
                                        EXX , EYY , EXY ,              &
                                        EXZ , EYZ , EZZ ,              &
                                        EXXF, EYYF, EXYF,              &
                                        EXZF, EYZF, EZZF

  USE WAVEFIELD   , ONLY:                                              &
                          SEISU , SEISV , SEISW ,                      &
                          SNAPU , SNAPV , SNAPW ,                      &
                          SNAPUF, SNAPVF, SNAPWF

#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: ALLOSTAT

!-----------------------------------------------------------------------

  IF ( LPAS /= 0 ) THEN
    ALLOCATE ( EXXF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF),            &
               EYYF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF),            &
               EXYF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF),            &
               EZZF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF),            &
               EXZF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF),            &
               EYZF (X_MIN_DF:X_MAX_DF, Y_MIN_DF:Y_MAX_DF) )
  END IF

  ALLOCATE ( EXX (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D),                   &
             EYY (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D),                   &
             EXY (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D),                   &
             EZZ (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D),                   &
             EXZ (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D),                   &
             EYZ (X_MIN_D:X_MAX_D, Y_MIN_D:Y_MAX_D),                   &
                                                        STAT=ALLOSTAT )
    IF  ( ALLOSTAT > 0 )  THEN
      WRITE (11,*) '  ALLOCATION ERROR ALLOC_T_SEIS_2, STAT= ', ALLOSTAT
      STOP
    END IF

  ALLOCATE (  SEISU (1:MR),  SEISV (1:MR),  SEISW (1:MR),              &
                                                        STAT=ALLOSTAT )
      IF  ( ALLOSTAT > 0 )  THEN
        WRITE (11,*) 'ALLOCATION ERROR ALLOC_T_SEIS_5, STAT= ', ALLOSTAT
        STOP
      END IF

#if(USE_MPI)
  IF ( myid == root ) THEN
    ALLOCATE ( SNAPU(PMLL:PMLXH,PMLL:PMLYH),                           &
               SNAPV(PMLL:PMLXH,PMLL:PMLYH),                           &
               SNAPW(PMLL:PMLXH,PMLL:PMLYH)   )

    IF ( LPAS /= 0 ) THEN
      ALLOCATE ( SNAPUF(PMLLF:PMLXHF,PMLLF:PMLYHF),                    &
                 SNAPVF(PMLLF:PMLXHF,PMLLF:PMLYHF),                    &
                 SNAPWF(PMLLF:PMLXHF,PMLLF:PMLYHF)   )
    END IF

  END IF
#endif

END SUBROUTINE ALLOC_T_SEIS

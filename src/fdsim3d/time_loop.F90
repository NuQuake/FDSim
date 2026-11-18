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

SUBROUTINE  TIME_LOOP

  USE CONTROL_DATA, ONLY:                                              &
                                     KEY_SNV, KEY_TLV, KEY_EN,         &
                          MZ, MT1, MT2, DT,                            &
                          REFINE, IPAS2,        TPML, LPAS,            &
                          PWINC, LIN, STRESS_IMAGING
  USE GRID_MEDIUM , ONLY:                                              &
                          H
  USE T_LOOP      , ONLY:                                              &
                          ITILE
  USE AUXIL       , ONLY:                                              &
                          PMLZH, PMLZH1, PMLZH2,                       &
                          Z_MIN, Z_MAX, Z_MINF, Z_MAXF, SOURTF_A
  USE WAVEFIELD   , ONLY:                                              &
                          UM, VM, WM, TOTAL_ENERGY
  USE TIME_FD     , ONLY:                                              &
                          time_begin, time_end, DATE_TIME, BIG_BEN
  USE NONREF_BOUND, ONLY:                                              &
                          SZ_UP, SZ_VP, SZ_WP, SZ_U, SZ_V, SZ_W
  USE INTERFACES  , ONLY:                                              &
                    HPLANE_0_STRESS       , HPLANE_0_VELOCITY ,        &
                    HPLANE_0F_STRESS      , HPLANE_0F_VELOCITY,        &
                    HPLANE_1_STRESS       , HPLANE_1_VELOCITY ,        &
                    HPLANE_1F_STRESS      , HPLANE_1F_VELOCITY,        &
                    HPLANE_IC_STRESS      , HPLANE_IC_VELOCITY,        &
                    HPLANE_IF_STRESS      , HPLANE_IF_VELOCITY,        &
                    HPLANE_LINPM1_STRESS  , HPLANE_LINPM1_VELOCITY,    &
                    HPLANE_LIN_STRESS     , HPLANE_LIN_VELOCITY,       &
                    HPLANE_LIN_W_STRESS   , HPLANE_LIN_W_VELOCITY,     &
                    HPLANE_LPASF_STRESS   ,                            &
                    HPLANE_LPASM1F_STRESS , HPLANE_LPASM1F_VELOCITY,   &
                    HPLANE_LPASM1M2_STRESS, HPLANE_LPASM1M2_VELOCITY,  &
                    HPLANE_MZ1_STRESS     ,                            &
                    HPLANE_PML_STRESS     , HPLANE_PML_VELOCITY,       &
                    HPLANE_PMLZH0_STRESS  , HPLANE_PMLZH0_VELOCITY,    &
                    SRC_ADD, STORE_SEISMOGRAMS,                        &
                    STORE_SNAPSHOTS, STORE_SNAPSHOTS_F
#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif

!-----------------------------------------------------------------------

  IMPLICIT NONE

  LOGICAL :: EX

  INTEGER :: L

!-----------------------------------------------------------------------

  INTERFACE

    SUBROUTINE     TL_ALLOC
    END SUBROUTINE TL_ALLOC

    SUBROUTINE     TL_ALLOCF
    END SUBROUTINE TL_ALLOCF

    SUBROUTINE     ENERGY   (ITILE)
      INTEGER, INTENT(IN) :: ITILE
    END SUBROUTINE ENERGY

  END INTERFACE

!-----------------------------------------------------------------------

  CALL TL_ALLOC

  IF ( LPAS /= 0 ) THEN
    CALL TL_ALLOCF
  END IF

  DO ITILE = MT1, MT2

    TOTAL_ENERGY = 0.

#if(USE_MPI)
    IF ( myid == root ) THEN
#endif
      CALL CPU_TIME ( time_end )
      CALL DATE_AND_TIME ( BIG_BEN(1),BIG_BEN(2),BIG_BEN(3), DATE_TIME )
      PRINT     *, 'CPU time: ', time_end - time_begin, ' seconds'
      WRITE (*,'(''Current time: '',I2,'':'',I2,'':'',F6.3)')            &
      DATE_TIME(5), DATE_TIME(6), DATE_TIME(7)+REAL(DATE_TIME(8)/1000.)
      PRINT     *, 'time level:', ITILE
      WRITE (11,'(I6,1X,F10.2,1X,I4,I3,I3,I5,3I3,I4)')                   &
                ITILE, time_end - time_begin, DATE_TIME

      INQUIRE ( FILE = 'STOP', EXIST = EX )
      IF ( EX ) THEN
        MT2 = ITILE-1
        OPEN  ( 35, FILE='STOP', STATUS='UNKNOWN')
        WRITE ( 35,* ) ITILE
        CLOSE ( 35 )
        STOP
      END IF

#if(USE_MPI)
    END IF
#endif

!-- UPDATE STRESS FIELD

    IF ( LPAS /= 0 ) THEN

      IF ( Z_MINF == 0 ) THEN
        DO L =     2, MIN (Z_MAXF, LPAS-2)
          CALL  HPLANE_IF_STRESS ( L )
        END DO
        CALL  HPLANE_0F_STRESS
        CALL  HPLANE_1F_STRESS
      ELSE
        DO L = Z_MINF, MIN (Z_MAXF, LPAS-2)
          CALL  HPLANE_IF_STRESS ( L )
        END DO
      END IF
      IF ( ( Z_MINF <= LPAS-1 ) .AND. ( LPAS-1 <= Z_MAXF ) ) THEN
        CALL  HPLANE_LPASM1F_STRESS
      END IF

      DO L = MAX(Z_MIN, LPAS), MIN (Z_MAX, MZ-2)
        CALL  HPLANE_IC_STRESS ( L )
      END DO

    ELSE

      IF      ( PWINC == 1 ) THEN
        DO L = MAX(Z_MIN, 2), MIN (Z_MAX, LIN-2)
          CALL  HPLANE_IC_STRESS ( L )
        END DO
        IF ( ( Z_MIN <= LIN-1 ) .AND. ( LIN-1 <= Z_MAX ) ) THEN
          CALL  HPLANE_LINPM1_STRESS (LIN-1)
        END IF
        IF ( ( Z_MIN <= LIN   ) .AND. ( LIN   <= Z_MAX ) ) THEN
           CALL  HPLANE_LIN_STRESS ( SOURTF_LIN = SOURTF_A(1,ITILE), &
                                     SOURTF_LOU = SOURTF_A(2,ITILE) )
        END IF
        IF ( ( Z_MIN <= LIN+1 ) .AND. ( LIN+1 <= Z_MAX ) ) THEN
          CALL  HPLANE_LINPM1_STRESS (LIN+1)
        END IF
        DO L = MAX(Z_MIN, LIN+2), MIN (Z_MAX, MZ-2)
          CALL  HPLANE_IC_STRESS ( L )
        END DO
      ELSE IF ( PWINC == 2 ) THEN
        DO L = MAX(Z_MIN, 2), MIN (Z_MAX, LIN-2)
          CALL  HPLANE_IC_STRESS ( L )
        END DO
        IF ( ( Z_MIN <= LIN-1 ) .AND. ( LIN-1 <= Z_MAX ) ) THEN
          CALL  HPLANE_LINPM1_STRESS (LIN-1)
        END IF
        IF ( ( Z_MIN <= LIN   ) .AND. ( LIN   <= Z_MAX ) ) THEN
           CALL  HPLANE_LIN_W_STRESS(SOURTF_LIN = SOURTF_A(1,ITILE), &
                                     SOURTF_LOU = SOURTF_A(2,ITILE) )
        END IF
        IF ( ( Z_MIN <= LIN+1 ) .AND. ( LIN+1 <= Z_MAX ) ) THEN
          CALL  HPLANE_LINPM1_STRESS (LIN+1)
        END IF
        DO L = MAX(Z_MIN, LIN+2), MIN (Z_MAX, MZ-2)
          CALL  HPLANE_IC_STRESS ( L )
        END DO
      ELSE
        DO L = MAX(Z_MIN, 2), MIN (Z_MAX, MZ-2)
          CALL  HPLANE_IC_STRESS ( L )
        END DO
      END IF

      IF ( Z_MIN == 0 ) THEN
        CALL  HPLANE_0_STRESS
        CALL  HPLANE_1_STRESS
      END IF

    END IF

    IF ( ( Z_MIN <= MZ-1 ) .AND. ( MZ-1 <= Z_MAX ) ) THEN
      CALL  HPLANE_MZ1_STRESS
    END IF

    DO L = MAX(Z_MIN, MZ), MIN(Z_MAX, PMLZH1)
      CALL  HPLANE_PML_STRESS ( L, L-MZ+1 )
    END DO

    IF ( Z_MAX == PMLZH ) THEN
      CALL  HPLANE_PMLZH0_STRESS
    END IF

!--- COMMUNICATE STRESSFIELD

#if(USE_MPI)
    CALL COMMUNICATE_STRESSFIELD
#endif

!--- INTERPOLATE STRESSFIELD ON LPAS PLANE

    IF ( LPAS /= 0 ) THEN
      IF ( ( Z_MINF <= LPAS ) .AND. ( LPAS <= Z_MAXF ) ) THEN
        CALL  HPLANE_LPASF_STRESS
      END IF
    END IF

!--- COPY STRESSFIELD TO PLANES LPAS-1 AND LPAS-2

    IF ( LPAS /= 0 ) THEN
      IF ( ( Z_MINF <= LPAS-(REFINE*3+1)/2 ) .AND. ( LPAS <= Z_MAXF ) ) THEN
        CALL  HPLANE_LPASM1M2_STRESS
      END IF
#if(USE_MPI)
      CALL MPI_BARRIER (comm_cart,errcode)
#endif
    END IF

!--- STORE BOTTOM BOUNDARY

    IF ( Z_MAX == PMLZH ) THEN
      SZ_UP(:,:)         = SZ_U (:,:)
      SZ_VP(:,:)         = SZ_V (:,:)
      SZ_WP(:,:)         = SZ_W (:,:)

      SZ_U (:,:)%TIKMZ2  = UM( :, :, PMLZH2 )
      SZ_V (:,:)%TIKMZ2  = VM( :, :, PMLZH2 )
      SZ_W (:,:)%TIKMZ2  = WM( :, :, PMLZH2 )

      SZ_U (:,:)%TIKMZ1  = UM( :, :, PMLZH1 )
      SZ_V (:,:)%TIKMZ1  = VM( :, :, PMLZH1 )
      SZ_W (:,:)%TIKMZ1  = WM( :, :, PMLZH1 )

      SZ_U (:,:)%TIKMZ0  = UM( :, :, PMLZH  )
      SZ_V (:,:)%TIKMZ0  = VM( :, :, PMLZH  )
      SZ_W (:,:)%TIKMZ0  = WM( :, :, PMLZH  )
    END IF

!--------------------------

    IF ( LPAS /= 0 ) THEN

      IF ( Z_MINF == 0 ) THEN
        DO L =     2, MIN (Z_MAXF, LPAS-2)
          CALL  HPLANE_IF_VELOCITY ( L )
        END DO
        IF (STRESS_IMAGING) THEN
          CALL  HPLANE_1F_VELOCITY
        ELSE
          CALL  HPLANE_0F_VELOCITY
          CALL  HPLANE_1F_VELOCITY
        END IF
      ELSE
        DO L = Z_MINF, MIN (Z_MAXF, LPAS-2)
          CALL  HPLANE_IF_VELOCITY ( L )
        END DO
      END IF
      IF ( ( Z_MINF <= LPAS-1 ) .AND. ( LPAS-1 <= Z_MAXF ) ) THEN
        CALL  HPLANE_LPASM1F_VELOCITY
      END IF

      DO L = MAX(Z_MIN, LPAS), MIN(Z_MAX, MZ-2)
        CALL  HPLANE_IC_VELOCITY ( L )
      END DO

    ELSE

      IF ( PWINC == 1 ) THEN
        DO L = MAX(Z_MIN, 2), MIN (Z_MAX, LIN-2)
          CALL  HPLANE_IC_VELOCITY ( L )
        END DO
        IF ( ( Z_MIN <= LIN-1 ) .AND. ( LIN-1 <= Z_MAX ) ) THEN
          CALL  HPLANE_LINPM1_VELOCITY (LIN-1)
        END IF
        IF ( ( Z_MIN <= LIN   ) .AND. ( LIN   <= Z_MAX ) ) THEN
           CALL  HPLANE_LIN_VELOCITY
        END IF
        IF ( ( Z_MIN <= LIN+1 ) .AND. ( LIN+1 <= Z_MAX ) ) THEN
          CALL  HPLANE_LINPM1_VELOCITY (LIN+1)
        END IF
        DO L = MAX(Z_MIN, LIN+2), MIN (Z_MAX, MZ-2)
          CALL  HPLANE_IC_VELOCITY ( L )
        END DO
      ELSE IF ( PWINC == 2 ) THEN
        DO L = MAX(Z_MIN, 2), MIN (Z_MAX, LIN-2)
          CALL  HPLANE_IC_VELOCITY ( L )
        END DO
        IF ( ( Z_MIN <= LIN-1 ) .AND. ( LIN-1 <= Z_MAX ) ) THEN
          CALL  HPLANE_LINPM1_VELOCITY (LIN-1)
        END IF
        IF ( ( Z_MIN <= LIN   ) .AND. ( LIN   <= Z_MAX ) ) THEN
           CALL  HPLANE_LIN_W_VELOCITY
        END IF
        IF ( ( Z_MIN <= LIN+1 ) .AND. ( LIN+1 <= Z_MAX ) ) THEN
          CALL  HPLANE_LINPM1_VELOCITY (LIN+1)
        END IF
        DO L = MAX(Z_MIN, LIN+2), MIN (Z_MAX, MZ-2)
          CALL  HPLANE_IC_VELOCITY ( L )
        END DO
      ELSE
        DO L = MAX(Z_MIN,2), MIN (Z_MAX, MZ-2)
          CALL  HPLANE_IC_VELOCITY ( L )
        END DO
      END IF

      IF ( Z_MIN == 0 ) THEN
        IF (STRESS_IMAGING) THEN
          CALL  HPLANE_1_VELOCITY
        ELSE
          CALL  HPLANE_0_VELOCITY
          CALL  HPLANE_1_VELOCITY
        END IF
      END IF
    END IF

    DO L = MAX(Z_MIN, MZ-1), MIN(Z_MAX, PMLZH1)
      CALL  HPLANE_PML_VELOCITY ( L, L-MZ+1 )
    END DO

    IF ( Z_MAX == PMLZH ) THEN
      CALL  HPLANE_PMLZH0_VELOCITY
    END IF

!--- COMMUNICATE WAVEFIELD

#if(USE_MPI)
    CALL COMMUNICATE_WAVEFIELD
#endif

!--- INTERPOLATE WAVEFIELD ON LPAS PLANE

    IF ( LPAS /= 0 ) THEN
      IF ( ( Z_MINF <= LPAS ) .AND. ( LPAS <= Z_MAXF ) ) THEN
        CALL  COARSE2FINE_W
      END IF
    END IF

!--- COPY WAVEFIELD TO PLANES LPAS-1 AND LPAS-2

    IF ( LPAS /= 0 ) THEN
      IF ( ( Z_MINF <= LPAS-(REFINE*3+1)/2 ) .AND. ( LPAS <= Z_MAXF ) ) THEN
        CALL  HPLANE_LPASM1M2_VELOCITY
      END IF
    END IF

!--- SOURCE

    CALL SRC_ADD ( ITILE )

!--- SYNCHRONIZATION

#if(USE_MPI)
    CALL MPI_BARRIER (comm_cart,errcode)
#endif

!--- STORE RESULTS

    IF ( KEY_TLV ) THEN
      CALL STORE_SEISMOGRAMS (ITILE)
    END IF

    IF ( KEY_SNV .AND. (MOD(ITILE,IPAS2)==0) ) THEN
      CALL STORE_SNAPSHOTS(ITILE)
      IF ( LPAS /= 0 ) THEN
        CALL STORE_SNAPSHOTS_F(ITILE)
      END IF
    END IF

    IF ( KEY_EN ) CALL ENERGY(ITILE)

  END DO


END  SUBROUTINE  TIME_LOOP

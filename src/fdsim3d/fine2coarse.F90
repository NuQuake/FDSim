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

SUBROUTINE  FINE2COARSE (L, L_F, VAR, VARF, I_SHIFT, K_SHIFT)

  USE PRECISION   , ONLY: PP

  USE CONTROL_DATA, ONLY: REFINE

  USE AUXIL       , ONLY:  X_MIN_D , X_MAX_D , X_MIN , X_MAX ,         &
                           Y_MIN_D , Y_MAX_D , Y_MIN , Y_MAX ,         &
                           Z_MIN_D , Z_MAX_D , Z_MIN , Z_MAX ,         &
                           X_MIN_DF, X_MAX_DF,                         &
                           Y_MIN_DF, Y_MAX_DF,                         &
                           Z_MIN_DF, Z_MAX_DF,  PMLLF, TMP,            &
                           Z_MINF  , Z_MAXF ,                          &
                           LANC, REF2, LLAN, VARF1

#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif
!-----------------------------------------------------------------------

  IMPLICIT NONE

  LOGICAL , INTENT (IN)    :: I_SHIFT, K_SHIFT
  INTEGER , INTENT (IN)    :: L, L_F

  REAL(PP), INTENT (INOUT), DIMENSION                                  &
       (X_MIN_D :X_MAX_D ,Y_MIN_D :Y_MAX_D ,Z_MIN_D :Z_MAX_D ) :: VAR
  REAL(PP), INTENT (INOUT), DIMENSION                                  &
       (X_MIN_DF:X_MAX_DF,Y_MIN_DF:Y_MAX_DF,Z_MIN_DF:Z_MAX_DF) :: VARF

  INTEGER                  :: I, K, J, I_F, K_F, SIZE1, I1, K1

!-----------------------------------------------------------------------

  IF ( ( Z_MINF <= L_F ) .AND. ( L_F <= Z_MAXF ) ) THEN

    VARF1  ( X_MIN_DF:X_MAX_DF,Y_MIN_DF:Y_MAX_DF      ) =              &
      VARF ( X_MIN_DF:X_MAX_DF,Y_MIN_DF:Y_MAX_DF, L_F )

    DO J = 1, LLAN*REFINE
      VARF1 ( X_MIN_DF-J, Y_MIN_DF : Y_MAX_DF      ) =                 &
       VARF1( X_MIN_DF+J, Y_MIN_DF : Y_MAX_DF      )
      VARF1 ( X_MAX_DF+J, Y_MIN_DF : Y_MAX_DF      ) =                 &
       VARF1( X_MAX_DF-J, Y_MIN_DF : Y_MAX_DF      )
    END DO

    DO J = 1, LLAN*REFINE
      VARF1 ( :, Y_MIN_DF-J ) = VARF1 ( :, Y_MIN_DF+J )
      VARF1 ( :, Y_MAX_DF+J ) = VARF1 ( :, Y_MAX_DF-J )
    END DO

    DO K = Y_MIN, Y_MAX
      DO I = X_MIN, X_MAX

        I_F =  (  I-1)*REFINE+1
        K_F =  (  K-1)*REFINE+1
        IF ( I_SHIFT ) THEN
          I_F = ((2*I-3)*REFINE+3)/2
          IF ( I_F < PMLLF ) I_F = PMLLF
        END IF
        IF ( K_SHIFT ) THEN
          K_F = ((2*K-3)*REFINE+3)/2
          IF ( K_F < PMLLF ) K_F = PMLLF
        END IF

        TMP (I,K) = 0.

        DO K1 = -REF2, REF2
          DO I1 = -REF2, REF2
            TMP(I,K) = TMP(I,K) + LANC(I1,K1) * VARF1 (I_F+I1,K_F+K1)
          END DO
        END DO

      END DO
    END DO

    IF ( ( Z_MIN <= L ) .AND. ( L <= Z_MAX ) ) THEN

      VAR (X_MIN:X_MAX,Y_MIN:Y_MAX,L) = TMP (X_MIN:X_MAX,Y_MIN:Y_MAX)

#if(USE_MPI)
    ELSE

      DO J = 0, numprocs-1

        CALL MPI_CART_COORDS (comm_cart, J, ndims, coord,errcode)

        IF ( ( Z_MIN_A(coord(3)) <= L ) .AND.                          &
                                  ( L <= Z_MAX_A(coord(3)) ) .AND.     &
             ( X_MIN == X_MIN_A(coord(1)) ) .AND.                      &
             ( X_MAX == X_MAX_A(coord(1)) ) .AND.                      &
             ( Y_MIN == Y_MIN_A(coord(2)) ) .AND.                      &
             ( Y_MAX == Y_MAX_A(coord(2)) )  ) THEN

          SIZE1 = ( X_MAX - X_MIN + 1 ) * ( Y_MAX - Y_MIN + 1 )

          tag  = 0

          CALL MPI_SEND (TMP, SIZE1, MPI_REAL, J, tag, comm_cart,      &
                                                       errcode)

        END IF

      END DO
#endif
    END IF
#if(USE_MPI)
  ELSE

    IF ( ( Z_MIN <= L ) .AND. ( L <= Z_MAX ) ) THEN

      DO J = 0, numprocs-1

        CALL MPI_CART_COORDS (comm_cart, J, ndims, coord,errcode)

        IF ( ( Z_MIN_AF(coord(3)) <= L_F ) .AND.                       &
                                   ( L_F <= Z_MAX_AF(coord(3)) ) .AND. &
             ( X_MIN == X_MIN_A(coord(1)) ) .AND.                      &
             ( X_MAX == X_MAX_A(coord(1)) ) .AND.                      &
             ( Y_MIN == Y_MIN_A(coord(2)) ) .AND.                      &
             ( Y_MAX == Y_MAX_A(coord(2)) )  ) THEN
          SIZE1 = ( X_MAX - X_MIN + 1 ) * ( Y_MAX - Y_MIN + 1 )

          tag  = 0

          CALL MPI_RECV (TMP, SIZE1, MPI_REAL, J, tag, comm_cart,      &
                                                       status, errcode)

          VAR (X_MIN:X_MAX,Y_MIN:Y_MAX,L) = TMP (X_MIN:X_MAX,Y_MIN:Y_MAX)

        END IF

      END DO

    END IF
#endif
  END IF


END SUBROUTINE  FINE2COARSE

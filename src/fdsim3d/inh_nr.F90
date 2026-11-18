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

SUBROUTINE INH_NR  ( MXL, MXH, PXH, MYL, MYH, PYH,                     &
                          X_MIN, X_MAX, Y_MIN, Y_MAX,                  &
                     AXT, AYT, SX, SXP, SY, SYP, T, K1, I1, TPML)

  USE PRECISION   , ONLY:                                              &
                          PP
  USE CONTROL_DATA, ONLY:                                              &
                          LS, PWINC, REFINE
  USE AUXIL       , ONLY:                                              &
                          X_MAX_D
  USE NONREF_BOUND, ONLY:                                              &
                          NR_PAR_A_X, NR_DATA_X,                       &
                          NR_PAR_A_Y, NR_DATA_Y

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER          ,                               INTENT (IN)    ::   &
                            MXL, MXH, MYL, MYH, PXH, PYH, I1, K1, TPML,&
                           X_MIN, X_MAX, Y_MIN, Y_MAX
  TYPE (NR_PAR_A_X), DIMENSION (MXL:MXH         ), INTENT (INOUT) ::   &
                                                                    AXT
  TYPE (NR_PAR_A_Y), DIMENSION (MYL:MYH         ), INTENT (INOUT) ::   &
                                                                    AYT
  TYPE (NR_DATA_X) , DIMENSION (MXL:MXH         ), INTENT (INOUT) ::   &
                                                                SX, SXP
  TYPE (NR_DATA_Y) , DIMENSION (MYL:MYH         ), INTENT (INOUT) ::   &
                                                                SY, SYP
  REAL (KIND=PP)   , DIMENSION (MXL:MXH, MYL:MYH), INTENT (INOUT) ::   &
                                                                      T

  INTEGER           :: I, K, PXH1, PXH2, PYH1, PYH2, IT1, IT2, IT3,    &
                                                     KT1, KT2, KT3, RR
  TYPE (NR_PAR_A_X) :: AI
  TYPE (NR_PAR_A_Y) :: AK
  TYPE (NR_DATA_X)  :: SXI, SXPI
  TYPE (NR_DATA_Y)  :: SYK, SYPK

!-----------------------------------------------------------------------

  IT1  = I1  - TPML
  KT1  = K1  - TPML
  PXH1 = PXH - 1
  PXH2 = PXH - 2
  PYH1 = PYH - 1
  PYH2 = PYH - 2
  IT2  = IT1 + 1
  IT3  = IT1 + 2
  KT2  = KT1 + 1
  KT3  = KT1 + 2

  IF ( MXH == X_MAX_D ) THEN
    RR = 1
  ELSE
    RR = REFINE
  END IF

  IF ( PWINC > 0 ) THEN
    IF ( PWINC == 2 ) THEN
!    P wave
      IF ( I1 + K1 == 3 ) THEN
!     U or V component
!       IF ( IT1 == 2 ) THEN
!       U component
          IF ( ( X_MIN <= IT1 ) .AND. ( IT1 <= X_MAX ) ) THEN
            T (IT1  ,  :  ) =  0.
            T (IT2  ,  :  ) =  0.
            T (IT3  ,  :  ) =  0.
          END IF
          IF ( PXH == X_MAX ) THEN
            T (PXH  ,  :  ) =  0.
            T (PXH1 ,  :  ) =  0.
            T (PXH2 ,  :  ) =  0.
          END IF
          IF ( ( Y_MIN <= KT1 ) .AND. ( KT1 <= Y_MAX ) ) THEN
            T (   :, KT1  ) =  0.
            T (   :, KT2  ) =  0.
            T (   :, KT3  ) =  0.
          END IF
          IF ( PYH == Y_MAX ) THEN
            T (   :, PYH  ) =  0.
            T (   :, PYH1 ) =  0.
            T (   :, PYH2 ) =  0.
          END IF
!       END IF
      ELSE
!      W component
        IF ( ( X_MIN <= IT1 ) .AND. ( IT1 <= X_MAX ) ) THEN
          T (IT1  ,  :  ) =  T (    7-TPML, :    )
          T (IT2  ,  :  ) =  T (    6-TPML, :    )
          T (IT3  ,  :  ) =  T (    5-TPML, :    )
        END IF
        IF ( PXH == X_MAX ) THEN
          T (PXH  ,  :  ) =  T (PXH-6     , :    )
          T (PXH1 ,  :  ) =  T (PXH-5     , :    )
          T (PXH2 ,  :  ) =  T (PXH-4     , :    )
        END IF
        IF ( ( Y_MIN <= KT1 ) .AND. ( KT1 <= Y_MAX ) ) THEN
          T (   :, KT1  ) =  T (    :, 7-TPML    )
          T (   :, KT2  ) =  T (    :, 6-TPML    )
          T (   :, KT3  ) =  T (    :, 5-TPML    )
        END IF
        IF ( PYH == Y_MAX ) THEN
          T (   :, PYH  ) =  T (    :, PYH-6     )
          T (   :, PYH1 ) =  T (    :, PYH-5     )
          T (   :, PYH2 ) =  T (    :, PYH-4     )
        END IF
      END IF
    ELSE IF ( PWINC == 1 ) THEN
!    S wave
      IF ( I1 + K1 == 3 ) THEN
!      U or V component
        IF ( I1 == 2 ) THEN
!      U component
          IF ( ( X_MIN <= IT1 ) .AND. ( IT1 <= X_MAX ) ) THEN
            DO I = 1, RR
              T (IT1+I-1,  :  ) =  T (2*RR+1-TPML-I+1, :    )
            END DO
          END IF
          IF ( PXH == X_MAX ) THEN
            DO I = 1, RR
              T (PXH-I+1,  :  ) =  T (PXH-2*RR+1+I-1 , :    )
            END DO
          END IF
          IF ( ( Y_MIN <= KT1 ) .AND. ( KT1 <= Y_MAX ) ) THEN
            DO K = 1, RR
              T (   :, KT1+K-1) =  T (    :, 2*RR+1-TPML-K+1)
            END DO
          END IF
          IF ( PYH == Y_MAX ) THEN
            DO K = 1, RR
              T (   :, PYH-K+1) =  T (    :, PYH-2*RR+  K-1 )
            END DO
          END IF
        ELSE
!      V component
          IF ( ( X_MIN <= IT1 ) .AND. ( IT1 <= X_MAX ) ) THEN
            DO I = 1, RR
              T (IT1+I-1,  :  ) =  T (2*RR+1-TPML-I+1, :    )
            END DO
          END IF
          IF ( PXH == X_MAX ) THEN
            DO I = 1, RR
              T (PXH-I+1,  :  ) =  T (PXH-2*RR  +I-1 , :    )
            END DO
          END IF
          IF ( ( Y_MIN <= KT1 ) .AND. ( KT1 <= Y_MAX ) ) THEN
            DO K = 1, RR
              T (   :, KT1+K-1) =  T (    :, 2*RR+1-TPML-K+1)
            END DO
          END IF
          IF ( PYH == Y_MAX ) THEN
            DO K = 1, RR
              T (   :, PYH-K+1) =  T (    :, PYH-2*RR+1+K-1 )
            END DO
          END IF
        END IF
      ELSE
!      W component
        IF ( ( X_MIN <= IT1 ) .AND. ( IT1 <= X_MAX ) ) THEN
          DO I = 1, RR
            T (IT1+I-1,  :  ) =  0
          END DO
        END IF
        IF ( PXH == X_MAX ) THEN
          DO I = 1, RR
            T (PXH-I+1,  :  ) =  0
          END DO
        END IF
        IF ( ( Y_MIN <= KT1 ) .AND. ( KT1 <= Y_MAX ) ) THEN
          DO K = 1, RR
            T (   :, KT1+K-1) =  0
          END DO
        END IF
        IF ( PYH == Y_MAX ) THEN
          DO K = 1, RR
            T (   :, PYH-K+1) =  0
          END DO
        END IF
      END IF
    END IF

  ELSE
!-----------------------------------------------------------------------
    DO  K = MAX(KT2,Y_MIN), MIN(PYH1,Y_MAX)

      AK   = AYT (K)
      SYK  = SY  (K)
      SYPK = SYP (K)
!_____________________________________________________________ LEFT ____

      IF ( ( X_MIN <= IT1 ) .AND. ( IT1 <= X_MAX ) ) THEN

        T (IT1,K) =                AK%ALE01 *    T(IT2,K)              &
                                 + AK%ALE02 *    T(IT3,K)              &
          + AK%ALE10 *  SYK%T1KL + AK%ALE11 *  SYK%T2KL                &
                                 + AK%ALE12 *  SYK%T3KL                &
          + AK%ALE20 * SYPK%T1KL + AK%ALE21 * SYPK%T2KL                &
                                 + AK%ALE22 * SYPK%T3KL

      END IF
!____________________________________________________________ RIGHT ____

      IF ( PXH == X_MAX ) THEN
        T(PXH,K) =                   AK%ARI01 *     T(PXH1,K)          &
                                   + AK%ARI02 *     T(PXH2,K)          &
          + AK%ARI10 *  SYK%TMX0KL + AK%ARI11 *  SYK%TMX1KL            &
                                   + AK%ARI12 *  SYK%TMX2KL            &
          + AK%ARI20 * SYPK%TMX0KL + AK%ARI21 * SYPK%TMX1KL            &
                                   + AK%ARI22 * SYPK%TMX2KL
      END IF
    END DO

!-----------------------------------------------------------------------
    DO  I = MAX(IT2,X_MIN), MIN(PXH1,X_MAX)

      AI   = AXT (I)
      SXI  = SX  (I)
      SXPI = SXP (I)

!_____________________________________________________________ REAR ____

      IF ( ( Y_MIN <= KT1 ) .AND. ( KT1 <= Y_MAX ) ) THEN
        T (I,KT1) =                AI%ARE01 *    T(I,KT2)              &
                                 + AI%ARE02 *    T(I,KT3)              &
          + AI%ARE10 *  SXI%TI1L + AI%ARE11 *  SXI%TI2L                &
                                 + AI%ARE12 *  SXI%TI3L                &
          + AI%ARE20 * SXPI%TI1L + AI%ARE21 * SXPI%TI2L                &
                                 + AI%ARE22 * SXPI%TI3L
      END IF
!____________________________________________________________ FRONT ____

      IF ( PYH == Y_MAX ) THEN
        T(I,PYH) =                   AI%AFR01 *     T(I,PYH1)          &
                                   + AI%AFR02 *     T(I,PYH2)          &
          + AI%AFR10 *  SXI%TIMY0L + AI%AFR11 *  SXI%TIMY1L            &
                                   + AI%AFR12 *  SXI%TIMY2L            &
          + AI%AFR20 * SXPI%TIMY0L + AI%AFR21 * SXPI%TIMY1L            &
                                   + AI%AFR22 * SXPI%TIMY2L
      END IF
    END DO

!    IF ( ( X_MIN == MXL ) .AND. ( Y_MIN == MYL ) ) THEN
!      T (IT1, KT1) = ( T (IT2 ,KT1) + T (IT1,KT2 ) ) / 2.
!    END IF
!
!    IF ( ( X_MIN == MXL ) .AND. ( Y_MAX == MYH ) ) THEN
!      T (IT1, PYH) = ( T (IT2 ,PYH) + T (IT1,PYH1) ) / 2.
!    END IF
!
!    IF ( ( X_MAX == MXH ) .AND. ( Y_MIN == MYL ) ) THEN
!      T (PXH, KT1) = ( T (PXH1,KT1) + T (PXH,KT2 ) ) / 2.
!    END IF
!
!    IF ( ( X_MAX == MXH ) .AND. ( Y_MAX == MYH ) ) THEN
!      T (PXH, PYH) = ( T (PXH1,PYH) + T (PXH,PYH1) ) / 2.
!  END IF

    IF ( ( X_MIN == MXL ) .AND. ( Y_MIN == MYL ) ) THEN
      T (IT1, KT1) = 0.
    END IF

    IF ( ( X_MIN == MXL ) .AND. ( Y_MAX == MYH ) ) THEN
      T (IT1, PYH) = 0.
    END IF

    IF ( ( X_MAX == MXH ) .AND. ( Y_MIN == MYL ) ) THEN
      T (PXH, KT1) = 0.
    END IF

    IF ( ( X_MAX == MXH ) .AND. ( Y_MAX == MYH ) ) THEN
      T (PXH, PYH) = 0.
    END IF

  END IF

END SUBROUTINE INH_NR

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

SUBROUTINE PRE_NR_COEF2

  USE PRECISION   , ONLY:                                              &
                          PP
  USE CONTROL_DATA, ONLY:                                              &
                          DT, KTLE, KTRI, KTRE, KTFR, LPAS
  USE AUXIL       , ONLY:                                              &
                           PMLLF, PMLXHF , PMLYHF ,                    &
                                             X_MINF, X_MAXF,           &
                                             Y_MINF, Y_MAXF,           &
                                             Z_MINF, Z_MAXF, HF
  USE GRID_MEDIUM , ONLY:                                              &
                          JMF, L2MX, L2MY, L2MZ,                       &
                          MUXY, MUYZ,       MUXZ,                      &
                          DENUF, DENVF, DENWF
  USE NONREF_BOUND, ONLY:                                              &
                          AXUF, AYUF, AXVF, AYVF, AXWF, AYWF

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER        :: I, L, K, JM1, JM2, JM3, JM4

  REAL (KIND=PP) :: VEL, VELS, VELP, KT, K1, K2

!-----------------------------------------------------------------------

  INTERFACE

    SUBROUTINE     NONREF (H, KEY, VEL, VELS, VELP,                    &
                           AX01,AX02,AX10,AX11,AX12,AX20,AX21,AX22)
      USE PRECISION, ONLY: PP
      INTEGER,        INTENT (IN)  :: KEY
      REAL (KIND=PP), INTENT (IN)  :: H, VEL, VELS, VELP
      REAL (KIND=PP), INTENT (OUT) ::                                  &
                           AX01,AX02,AX10,AX11,AX12,AX20,AX21,AX22
    END SUBROUTINE NONREF

  END INTERFACE

!-----------------------------------------------------------------------
!                                                   HETEROGENEOUS MEDIUM
  DO  L = Z_MINF, MIN(Z_MAXF,LPAS-1)

!______________________________________________________ REAR TRANSPARENT

   IF ( Y_MINF == PMLLF ) THEN
    DO  I = X_MINF, X_MAXF
      JM1 = JMF ( I  , PMLLF  , L )
      JM2 = JMF ( I  , PMLLF+1, L )

      IF ( I > PMLLF ) THEN
        JM3 = JMF ( I-1, PMLLF, L )
      ELSE
        JM3 = JM1
      END IF

      K1  = L2MX(JM1)
      K2  = L2MX(JM3)
      KT  = 2._PP*K1*K2/(K1+K2)

      VELS = SQRT( MUXY(JM2) *HF*DENUF(JM1)/DT )
      VELP = SQRT( KT        *HF*DENUF(JM1)/DT )
      VEL  = VELS
      CALL NONREF ( HF, KTRE, VEL, VELS, VELP,                         &
                              AXUF(I,L)%ARE01, AXUF(I,L)%ARE02,        &
             AXUF(I,L)%ARE10, AXUF(I,L)%ARE11, AXUF(I,L)%ARE12,        &
             AXUF(I,L)%ARE20, AXUF(I,L)%ARE21, AXUF(I,L)%ARE22   )

      IF ( I < PMLXHF-1 ) THEN
        JM3 = JMF ( I+1, PMLLF+1, L )
      ELSE
        JM3 = JM2
      END IF
      JM4 = JMF ( I, PMLLF+1, L+1 )

      KT  = 4._PP / ( 1._PP/MUXY(JM2) + 1._PP/MUXY(JM3)                &
                    + 1._PP/MUYZ(JM2) + 1._PP/MUYZ(JM4) )

      VELS = SQRT( KT        *HF*DENVF(JM2)/DT )
      VELP = SQRT( L2MY(JM2) *HF*DENVF(JM2)/DT )
      VEL  = VELP
      CALL NONREF ( HF, KTRE, VEL, VELS, VELP,                         &
                              AXVF(I,L)%ARE01, AXVF(I,L)%ARE02,        &
             AXVF(I,L)%ARE10, AXVF(I,L)%ARE11, AXVF(I,L)%ARE12,        &
             AXVF(I,L)%ARE20, AXVF(I,L)%ARE21, AXVF(I,L)%ARE22   )

      IF ( L > 0 ) THEN
        JM3 = JMF ( I, PMLLF, L-1 )
      ELSE
        JM3 = JM1
      END IF

      K1  = L2MZ(JM1)
      K2  = L2MZ(JM3)
      KT  = 2._PP*K1*K2/(K1+K2)

      VELS = SQRT( MUYZ(JM2) *HF*DENWF(JM1)/DT )
      VELP = SQRT( KT        *HF*DENWF(JM1)/DT )
      VEL  = VELS
      CALL NONREF ( HF, KTRE, VEL, VELS, VELP,                         &
                              AXWF(I,L)%ARE01, AXWF(I,L)%ARE02,        &
             AXWF(I,L)%ARE10, AXWF(I,L)%ARE11, AXWF(I,L)%ARE12,        &
             AXWF(I,L)%ARE20, AXWF(I,L)%ARE21, AXWF(I,L)%ARE22   )
    END DO
   END IF

!_____________________________________________________ FRONT TRANSPARENT

   IF ( Y_MAXF == PMLYHF ) THEN
    DO  I = X_MINF, X_MAXF
      JM1 = JMF ( I  , PMLYHF  , L )
      JM2 = JMF ( I  , PMLYHF-1, L )

      IF ( I > PMLLF ) THEN
        JM3 = JMF ( I-1, PMLYHF, L )
      ELSE
        JM3 = JM1
      END IF

      K1  = L2MX(JM1)
      K2  = L2MX(JM3)
      KT  = 2._PP*K1*K2/(K1+K2)

      VELS = SQRT( MUXY(JM1) *HF*DENUF(JM1)/DT )
      VELP = SQRT( KT        *HF*DENUF(JM1)/DT )
      VEL  = VELS
      CALL NONREF ( HF, KTFR, VEL, VELS, VELP,                         &
                              AXUF(I,L)%AFR01, AXUF(I,L)%AFR02,        &
             AXUF(I,L)%AFR10, AXUF(I,L)%AFR11, AXUF(I,L)%AFR12,        &
             AXUF(I,L)%AFR20, AXUF(I,L)%AFR21, AXUF(I,L)%AFR22   )

      IF ( I < PMLXHF-1 ) THEN
        JM3 = JMF ( I+1, PMLYHF, L )
      ELSE
        JM3 = JM1
      END IF
      JM4 = JMF ( I, PMLYHF, L+1 )

      KT  = 4._PP / ( 1._PP/MUXY(JM1) + 1._PP/MUXY(JM3)                &
                    + 1._PP/MUYZ(JM1) + 1._PP/MUYZ(JM4) )

      VELS = SQRT( KT        *HF*DENVF(JM1)/DT )
      VELP = SQRT( L2MY(JM2) *HF*DENVF(JM1)/DT )
      VEL  = VELP
      CALL NONREF ( HF, KTFR, VEL, VELS, VELP,                         &
                              AXVF(I,L)%AFR01, AXVF(I,L)%AFR02,        &
             AXVF(I,L)%AFR10, AXVF(I,L)%AFR11, AXVF(I,L)%AFR12,        &
             AXVF(I,L)%AFR20, AXVF(I,L)%AFR21, AXVF(I,L)%AFR22   )

      IF ( L > 0 ) THEN
        JM3 = JMF ( I, PMLYHF, L-1 )
      ELSE
        JM3 = JM1
      END IF

      K1  = L2MZ(JM1)
      K2  = L2MZ(JM3)
      KT  = 2._PP*K1*K2/(K1+K2)

      VELS = SQRT( MUYZ(JM1) *HF*DENWF(JM1)/DT )
      VELP = SQRT( KT        *HF*DENWF(JM1)/DT )
      VEL  = VELS
      CALL NONREF ( HF, KTFR, VEL, VELS, VELP,                         &
                              AXWF(I,L)%AFR01, AXWF(I,L)%AFR02,        &
             AXWF(I,L)%AFR10, AXWF(I,L)%AFR11, AXWF(I,L)%AFR12,        &
             AXWF(I,L)%AFR20, AXWF(I,L)%AFR21, AXWF(I,L)%AFR22   )
    END DO
   END IF
!______________________________________________________ LEFT TRANSPARENT

   IF ( X_MINF == PMLLF ) THEN
    DO  K = Y_MINF, Y_MAXF

     IF ( KTLE >= 0 ) THEN
      JM1 = JMF ( PMLLF  , K, L )
      JM2 = JMF ( PMLLF+1, K, L )

      IF ( K < PMLYHF-1 ) THEN
        JM3 = JMF ( PMLLF+1, K+1, L )
      ELSE
        JM3 = JM2
      END IF
      JM4 = JMF ( PMLLF+1, K, L+1 )

      KT  = 4._PP / ( 1._PP/MUXY(JM2) + 1._PP/MUXY(JM3)                &
                    + 1._PP/MUXZ(JM2) + 1._PP/MUXZ(JM4) )
      VELS = SQRT( KT        *HF*DENUF(JM2)/DT )
      VELP = SQRT( L2MX(JM2) *HF*DENUF(JM2)/DT )
      VEL  = VELP
      CALL NONREF ( HF, KTLE, VEL, VELS, VELP,                         &
                              AYUF(K,L)%ALE01, AYUF(K,L)%ALE02,        &
             AYUF(K,L)%ALE10, AYUF(K,L)%ALE11, AYUF(K,L)%ALE12,        &
             AYUF(K,L)%ALE20, AYUF(K,L)%ALE21, AYUF(K,L)%ALE22   )

      IF ( K > PMLLF ) THEN
        JM3 = JMF ( PMLLF, K-1, L )
      ELSE
        JM3 = JM1
      END IF

      K1  = L2MY(JM1)
      K2  = L2MY(JM3)
      KT  = 2._PP*K1*K2/(K1+K2)

      VELS = SQRT( MUXY(JM2) *HF*DENVF(JM1)/DT )
      VELP = SQRT( KT        *HF*DENVF(JM1)/DT )
      VEL  = VELS
      CALL NONREF ( HF, KTLE, VEL, VELS, VELP,                         &
                              AYVF(K,L)%ALE01, AYVF(K,L)%ALE02,        &
             AYVF(K,L)%ALE10, AYVF(K,L)%ALE11, AYVF(K,L)%ALE12,        &
             AYVF(K,L)%ALE20, AYVF(K,L)%ALE21, AYVF(K,L)%ALE22   )

      IF ( L > 0 ) THEN
        JM3 = JMF ( PMLLF, K , L-1 )
      ELSE
        JM3 = JM1
      END IF

      K1  = L2MZ(JM1)
      K2  = L2MZ(JM3)
      KT  = 2._PP*K1*K2/(K1+K2)

      VELS = SQRT( MUXZ(JM2) *HF*DENWF(JM1)/DT )
      VELP = SQRT( KT        *HF*DENWF(JM1)/DT )
      VEL  = VELS
      CALL NONREF ( HF, KTLE, VEL, VELS, VELP,                         &
                              AYWF(K,L)%ALE01, AYWF(K,L)%ALE02,        &
             AYWF(K,L)%ALE10, AYWF(K,L)%ALE11, AYWF(K,L)%ALE12,        &
             AYWF(K,L)%ALE20, AYWF(K,L)%ALE21, AYWF(K,L)%ALE22   )

     ELSE

      AYUF(K,L)%ALE01 = -1.
      AYUF(K,L)%ALE02 =  0.
      AYUF(K,L)%ALE10 =  0.
      AYUF(K,L)%ALE11 =  0.
      AYUF(K,L)%ALE12 =  0.
      AYUF(K,L)%ALE20 =  0.
      AYUF(K,L)%ALE21 =  0.
      AYUF(K,L)%ALE22 =  0.

      AYVF(K,L)%ALE01 =  0.
      AYVF(K,L)%ALE02 =  1.
      AYVF(K,L)%ALE10 =  0.
      AYVF(K,L)%ALE11 =  0.
      AYVF(K,L)%ALE12 =  0.
      AYVF(K,L)%ALE20 =  0.
      AYVF(K,L)%ALE21 =  0.
      AYVF(K,L)%ALE22 =  0.

      AYWF(K,L)%ALE01 =  0.
      AYWF(K,L)%ALE02 =  1.
      AYWF(K,L)%ALE10 =  0.
      AYWF(K,L)%ALE11 =  0.
      AYWF(K,L)%ALE12 =  0.
      AYWF(K,L)%ALE20 =  0.
      AYWF(K,L)%ALE21 =  0.
      AYWF(K,L)%ALE22 =  0.

     END IF

    END DO
   END IF

!_____________________________________________________ RIGHT TRANSPARENT

   IF ( X_MAXF == PMLXHF ) THEN
    DO  K = Y_MINF, Y_MAXF

     IF ( KTRI >= 0 ) THEN
      JM1 = JMF(PMLXHF  ,K,L)
      JM2 = JMF(PMLXHF-1,K,L)

      IF ( K < PMLYHF-1 ) THEN
        JM3 = JMF ( PMLXHF, K+1, L )
      ELSE
        JM3 = JM2
      END IF
      JM4 = JMF ( PMLXHF, K, L+1 )

      KT  = 4._PP / ( 1._PP/MUXY(JM1) + 1._PP/MUXY(JM3)                &
                    + 1._PP/MUXZ(JM1) + 1._PP/MUXZ(JM4) )
      VELS = SQRT( KT        *HF*DENUF(JM1)/DT )
      VELP = SQRT( L2MX(JM2) *HF*DENUF(JM1)/DT )
      VEL  = VELP
      CALL NONREF ( HF, KTRI, VEL, VELS, VELP,                         &
                              AYUF(K,L)%ARI01, AYUF(K,L)%ARI02,        &
             AYUF(K,L)%ARI10, AYUF(K,L)%ARI11, AYUF(K,L)%ARI12,        &
             AYUF(K,L)%ARI20, AYUF(K,L)%ARI21, AYUF(K,L)%ARI22   )

      IF ( K > PMLLF ) THEN
        JM3 = JMF ( PMLXHF, K-1, L )
      ELSE
        JM3 = JM1
      END IF

      K1  = L2MY(JM1)
      K2  = L2MY(JM3)
      KT  = 2._PP*K1*K2/(K1+K2)

      VELS = SQRT( MUXY(JM1) *HF*DENVF(JM1)/DT )
      VELP = SQRT( KT        *HF*DENVF(JM1)/DT )
      VEL  = VELS
      CALL NONREF ( HF, KTRI, VEL, VELS, VELP,                         &
                              AYVF(K,L)%ARI01, AYVF(K,L)%ARI02,        &
             AYVF(K,L)%ARI10, AYVF(K,L)%ARI11, AYVF(K,L)%ARI12,        &
             AYVF(K,L)%ARI20, AYVF(K,L)%ARI21, AYVF(K,L)%ARI22   )

      IF ( L > 0 ) THEN
        JM3 = JMF ( PMLXHF , K, L-1 )
      ELSE
        JM3 = JM1
      END IF

      K1  = L2MZ(JM1)
      K2  = L2MZ(JM3)
      KT  = 2._PP*K1*K2/(K1+K2)

      VELS = SQRT( MUXZ(JM1) *HF*DENWF(JM1)/DT )
      VELP = SQRT( KT        *HF*DENWF(JM1)/DT )
      VEL  = VELS
      CALL NONREF ( HF, KTRI, VEL, VELS, VELP,                         &
                              AYWF(K,L)%ARI01, AYWF(K,L)%ARI02,        &
             AYWF(K,L)%ARI10, AYWF(K,L)%ARI11, AYWF(K,L)%ARI12,        &
             AYWF(K,L)%ARI20, AYWF(K,L)%ARI21, AYWF(K,L)%ARI22   )

     ELSE

      AYUF(K,L)%ARI01 = -1.
      AYUF(K,L)%ARI02 =  0.
      AYUF(K,L)%ARI10 =  0.
      AYUF(K,L)%ARI11 =  0.
      AYUF(K,L)%ARI12 =  0.
      AYUF(K,L)%ARI20 =  0.
      AYUF(K,L)%ARI21 =  0.
      AYUF(K,L)%ARI22 =  0.

      AYVF(K,L)%ARI01 =  0.
      AYVF(K,L)%ARI02 =  1.
      AYVF(K,L)%ARI10 =  0.
      AYVF(K,L)%ARI11 =  0.
      AYVF(K,L)%ARI12 =  0.
      AYVF(K,L)%ARI20 =  0.
      AYVF(K,L)%ARI21 =  0.
      AYVF(K,L)%ARI22 =  0.

      AYWF(K,L)%ARI01 =  0.
      AYWF(K,L)%ARI02 =  1.
      AYWF(K,L)%ARI10 =  0.
      AYWF(K,L)%ARI11 =  0.
      AYWF(K,L)%ARI12 =  0.
      AYWF(K,L)%ARI20 =  0.
      AYWF(K,L)%ARI21 =  0.
      AYWF(K,L)%ARI22 =  0.

     END IF

    END DO
   END IF
  END DO


END SUBROUTINE PRE_NR_COEF2

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

SUBROUTINE PRE_NR_COEF1

  USE PRECISION   , ONLY:                                              &
                          PP
  USE CONTROL_DATA, ONLY:                                              &
                          MZ,  DT,                                     &
                          KTLE, KTRI, KTRE, KTFR, KTBO, LPAS
  USE AUXIL       , ONLY:                                              &
                           PMLL, PMLXH , PMLYH , PMLZH,                &
                                             X_MIN, X_MAX,             &
                                             Y_MIN, Y_MAX,             &
                                             Z_MIN, Z_MAX
  USE GRID_MEDIUM , ONLY:                                              &
                          JM, L2MX, L2MY, L2MZ,                        &
                          MUXY, MUYZ,       MUXZ,                      &
                          DENU, DENV, DENW, H
  USE NONREF_BOUND, ONLY:                                              &
                          AXU, AYU, AXV, AYV, AXW, AYW,                &
                          ABW01,   ABW02,   ABW10,   ABW11,            &
                          ABW12,   ABW20,   ABW21,   ABW22,            &
                          ABH01,   ABH02,   ABH10,   ABH11,            &
                          ABH12,   ABH20,   ABH21,   ABH22

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
  DO  L = MAX(Z_MIN,LPAS), MIN(Z_MAX,MZ-1)

!______________________________________________________ REAR TRANSPARENT

   IF ( Y_MIN == PMLL ) THEN
    DO  I = X_MIN, X_MAX
      JM1 = JM ( I  , PMLL  , L )
      JM2 = JM ( I  , PMLL+1, L )

      IF ( I > PMLL ) THEN
        JM3 = JM ( I-1, PMLL, L )
      ELSE
        JM3 = JM1
      END IF

      K1  = L2MX(JM1)
      K2  = L2MX(JM3)
      KT  = 2._PP*K1*K2/(K1+K2)

      VELS = SQRT( MUXY(JM2) *H*DENU(JM1)/DT )
      VELP = SQRT( KT        *H*DENU(JM1)/DT )
      VEL  = VELS
      CALL NONREF ( H, KTRE, VEL, VELS, VELP,                          &
                                AXU(I,L)%ARE01, AXU(I,L)%ARE02,        &
                AXU(I,L)%ARE10, AXU(I,L)%ARE11, AXU(I,L)%ARE12,        &
                AXU(I,L)%ARE20, AXU(I,L)%ARE21, AXU(I,L)%ARE22   )

      IF ( I < PMLXH-1 ) THEN
        JM3 = JM ( I+1, PMLL+1, L )
      ELSE
        JM3 = JM2
      END IF
      JM4 = JM ( I, PMLL+1, L+1 )

      KT  = 4._PP / ( 1._PP/MUXY(JM2) + 1._PP/MUXY(JM3)                &
                    + 1._PP/MUYZ(JM2) + 1._PP/MUYZ(JM4) )

      VELS = SQRT( KT        *H*DENV(JM2)/DT )
      VELP = SQRT( L2MY(JM2) *H*DENV(JM2)/DT )
      VEL  = VELP
      CALL NONREF ( H, KTRE, VEL, VELS, VELP,                          &
                                AXV(I,L)%ARE01, AXV(I,L)%ARE02,        &
                AXV(I,L)%ARE10, AXV(I,L)%ARE11, AXV(I,L)%ARE12,        &
                AXV(I,L)%ARE20, AXV(I,L)%ARE21, AXV(I,L)%ARE22   )

      IF ( L > 0 ) THEN
        JM3 = JM ( I, PMLL, L-1 )
      ELSE
        JM3 = JM1
      END IF

      K1  = L2MZ(JM1)
      K2  = L2MZ(JM3)
      KT  = 2._PP*K1*K2/(K1+K2)

      VELS = SQRT( MUYZ(JM2) *H*DENW(JM1)/DT )
      VELP = SQRT( KT        *H*DENW(JM1)/DT )
      VEL  = VELS
      CALL NONREF ( H, KTRE, VEL, VELS, VELP,                          &
                                AXW(I,L)%ARE01, AXW(I,L)%ARE02,        &
                AXW(I,L)%ARE10, AXW(I,L)%ARE11, AXW(I,L)%ARE12,        &
                AXW(I,L)%ARE20, AXW(I,L)%ARE21, AXW(I,L)%ARE22   )
    END DO
    IF ( Z_MAX >= MZ ) THEN
      DO I = X_MIN, X_MAX
        AXU(I,MZ:Z_MAX) = AXU(I,MZ-1)
        AXV(I,MZ:Z_MAX) = AXV(I,MZ-1)
        AXW(I,MZ:Z_MAX) = AXW(I,MZ-1)
      END DO
    END IF
   END IF

!_____________________________________________________ FRONT TRANSPARENT

   IF ( Y_MAX == PMLYH ) THEN
    DO  I = X_MIN, X_MAX
      JM1 = JM ( I  , PMLYH  , L )
      JM2 = JM ( I  , PMLYH-1, L )

      IF ( I > PMLL ) THEN
        JM3 = JM ( I-1, PMLYH, L )
      ELSE
        JM3 = JM1
      END IF

      K1  = L2MX(JM1)
      K2  = L2MX(JM3)
      KT  = 2._PP*K1*K2/(K1+K2)

      VELS = SQRT( MUXY(JM1) *H*DENU(JM1)/DT )
      VELP = SQRT( KT        *H*DENU(JM1)/DT )
      VEL  = VELS
      CALL NONREF ( H, KTFR, VEL, VELS, VELP,                          &
                                AXU(I,L)%AFR01, AXU(I,L)%AFR02,        &
                AXU(I,L)%AFR10, AXU(I,L)%AFR11, AXU(I,L)%AFR12,        &
                AXU(I,L)%AFR20, AXU(I,L)%AFR21, AXU(I,L)%AFR22   )

      IF ( I < PMLXH-1 ) THEN
        JM3 = JM ( I+1, PMLYH, L )
      ELSE
        JM3 = JM1
      END IF
      JM4 = JM ( I, PMLYH, L+1 )

      KT  = 4._PP / ( 1._PP/MUXY(JM1) + 1._PP/MUXY(JM3)                &
                    + 1._PP/MUYZ(JM1) + 1._PP/MUYZ(JM4) )

      VELS = SQRT( KT        *H*DENV(JM1)/DT )
      VELP = SQRT( L2MY(JM2) *H*DENV(JM1)/DT )
      VEL  = VELP
      CALL NONREF ( H, KTFR, VEL, VELS, VELP,                          &
                                AXV(I,L)%AFR01, AXV(I,L)%AFR02,        &
                AXV(I,L)%AFR10, AXV(I,L)%AFR11, AXV(I,L)%AFR12,        &
                AXV(I,L)%AFR20, AXV(I,L)%AFR21, AXV(I,L)%AFR22   )

      IF ( L > 0 ) THEN
        JM3 = JM ( I, PMLYH, L-1 )
      ELSE
        JM3 = JM1
      END IF

      K1  = L2MZ(JM1)
      K2  = L2MZ(JM3)
      KT  = 2._PP*K1*K2/(K1+K2)

      VELS = SQRT( MUYZ(JM1) *H*DENW(JM1)/DT )
      VELP = SQRT( KT        *H*DENW(JM1)/DT )
      VEL  = VELS
      CALL NONREF ( H, KTFR, VEL, VELS, VELP,                          &
                                AXW(I,L)%AFR01, AXW(I,L)%AFR02,        &
                AXW(I,L)%AFR10, AXW(I,L)%AFR11, AXW(I,L)%AFR12,        &
                AXW(I,L)%AFR20, AXW(I,L)%AFR21, AXW(I,L)%AFR22   )
    END DO
    IF ( Z_MAX >= MZ ) THEN
      DO I = X_MIN, X_MAX
        AXU(I,MZ:Z_MAX) = AXU(I,MZ-1)
        AXV(I,MZ:Z_MAX) = AXV(I,MZ-1)
        AXW(I,MZ:Z_MAX) = AXW(I,MZ-1)
      END DO
    END IF
   END IF
!______________________________________________________ LEFT TRANSPARENT

   IF ( X_MIN == PMLL ) THEN
    DO  K = Y_MIN, Y_MAX

     IF ( KTLE >=0 ) THEN

      JM1 = JM ( PMLL  , K, L )
      JM2 = JM ( PMLL+1, K, L )

      IF ( K < PMLYH-1 ) THEN
        JM3 = JM ( PMLL+1, K+1, L )
      ELSE
        JM3 = JM2
      END IF
      JM4 = JM ( PMLL+1, K, L+1 )

      KT  = 4._PP / ( 1._PP/MUXY(JM2) + 1._PP/MUXY(JM3)                &
                    + 1._PP/MUXZ(JM2) + 1._PP/MUXZ(JM4) )
      VELS = SQRT( KT        *H*DENU(JM2)/DT )
      VELP = SQRT( L2MX(JM2) *H*DENU(JM2)/DT )
      VEL  = VELP
      CALL NONREF ( H, KTLE, VEL, VELS, VELP,                          &
                                AYU(K,L)%ALE01, AYU(K,L)%ALE02,        &
                AYU(K,L)%ALE10, AYU(K,L)%ALE11, AYU(K,L)%ALE12,        &
                AYU(K,L)%ALE20, AYU(K,L)%ALE21, AYU(K,L)%ALE22   )

      IF ( K > PMLL ) THEN
        JM3 = JM ( PMLL, K-1, L )
      ELSE
        JM3 = JM1
      END IF

      K1  = L2MY(JM1)
      K2  = L2MY(JM3)
      KT  = 2._PP*K1*K2/(K1+K2)

      VELS = SQRT( MUXY(JM2) *H*DENV(JM1)/DT )
      VELP = SQRT( KT        *H*DENV(JM1)/DT )
      VEL  = VELS
      CALL NONREF ( H, KTLE, VEL, VELS, VELP,                          &
                                AYV(K,L)%ALE01, AYV(K,L)%ALE02,        &
                AYV(K,L)%ALE10, AYV(K,L)%ALE11, AYV(K,L)%ALE12,        &
                AYV(K,L)%ALE20, AYV(K,L)%ALE21, AYV(K,L)%ALE22   )

      IF ( L > 0 ) THEN
        JM3 = JM ( PMLL, K , L-1 )
      ELSE
        JM3 = JM1
      END IF

      K1  = L2MZ(JM1)
      K2  = L2MZ(JM3)
      KT  = 2._PP*K1*K2/(K1+K2)

      VELS = SQRT( MUXZ(JM2) *H*DENW(JM1)/DT )
      VELP = SQRT( KT        *H*DENW(JM1)/DT )
      VEL  = VELS
      CALL NONREF ( H, KTLE, VEL, VELS, VELP,                          &
                                AYW(K,L)%ALE01, AYW(K,L)%ALE02,        &
                AYW(K,L)%ALE10, AYW(K,L)%ALE11, AYW(K,L)%ALE12,        &
                AYW(K,L)%ALE20, AYW(K,L)%ALE21, AYW(K,L)%ALE22   )

     ELSE

      AYU (K,L)%ALE01 = -1.
      AYU (K,L)%ALE02 =  0.
      AYU (K,L)%ALE10 =  0.
      AYU (K,L)%ALE11 =  0.
      AYU (K,L)%ALE12 =  0.
      AYU (K,L)%ALE20 =  0.
      AYU (K,L)%ALE21 =  0.
      AYU (K,L)%ALE22 =  0.

      AYV (K,L)%ALE01 =  0.
      AYV (K,L)%ALE02 =  1.
      AYV (K,L)%ALE10 =  0.
      AYV (K,L)%ALE11 =  0.
      AYV (K,L)%ALE12 =  0.
      AYV (K,L)%ALE20 =  0.
      AYV (K,L)%ALE21 =  0.
      AYV (K,L)%ALE22 =  0.

      AYW (K,L)%ALE01 =  0.
      AYW (K,L)%ALE02 =  1.
      AYW (K,L)%ALE10 =  0.
      AYW (K,L)%ALE11 =  0.
      AYW (K,L)%ALE12 =  0.
      AYW (K,L)%ALE20 =  0.
      AYW (K,L)%ALE21 =  0.
      AYW (K,L)%ALE22 =  0.

     END IF

    END DO
    IF ( Z_MAX >= MZ ) THEN
      DO K = Y_MIN, Y_MAX
        AYU(K,MZ:Z_MAX) = AYU(K,MZ-1)
        AYV(K,MZ:Z_MAX) = AYV(K,MZ-1)
        AYW(K,MZ:Z_MAX) = AYW(K,MZ-1)
      END DO
    END IF
   END IF

!_____________________________________________________ RIGHT TRANSPARENT

   IF ( X_MAX == PMLXH ) THEN
    DO  K = Y_MIN, Y_MAX

     IF ( KTRI >=0 ) THEN

      JM1 = JM(PMLXH  ,K,L)
      JM2 = JM(PMLXH-1,K,L)

      IF ( K < PMLYH-1 ) THEN
        JM3 = JM ( PMLXH, K+1, L )
      ELSE
        JM3 = JM2
      END IF
      JM4 = JM ( PMLXH, K, L+1 )

      KT  = 4._PP / ( 1._PP/MUXY(JM1) + 1._PP/MUXY(JM3)                &
                    + 1._PP/MUXZ(JM1) + 1._PP/MUXZ(JM4) )
      VELS = SQRT( KT        *H*DENU(JM1)/DT )
      VELP = SQRT( L2MX(JM2) *H*DENU(JM1)/DT )
      VEL  = VELP
      CALL NONREF ( H, KTRI, VEL, VELS, VELP,                          &
                                AYU(K,L)%ARI01, AYU(K,L)%ARI02,        &
                AYU(K,L)%ARI10, AYU(K,L)%ARI11, AYU(K,L)%ARI12,        &
                AYU(K,L)%ARI20, AYU(K,L)%ARI21, AYU(K,L)%ARI22   )

      IF ( K > PMLL ) THEN
        JM3 = JM ( PMLXH, K-1, L )
      ELSE
        JM3 = JM1
      END IF

      K1  = L2MY(JM1)
      K2  = L2MY(JM3)
      KT  = 2._PP*K1*K2/(K1+K2)

      VELS = SQRT( MUXY(JM1) *H*DENV(JM1)/DT )
      VELP = SQRT( KT        *H*DENV(JM1)/DT )
      VEL  = VELS
      CALL NONREF ( H, KTRI, VEL, VELS, VELP,                          &
                                AYV(K,L)%ARI01, AYV(K,L)%ARI02,        &
                AYV(K,L)%ARI10, AYV(K,L)%ARI11, AYV(K,L)%ARI12,        &
                AYV(K,L)%ARI20, AYV(K,L)%ARI21, AYV(K,L)%ARI22   )

      IF ( L > 0 ) THEN
        JM3 = JM ( PMLXH , K, L-1 )
      ELSE
        JM3 = JM1
      END IF

      K1  = L2MZ(JM1)
      K2  = L2MZ(JM3)
      KT  = 2._PP*K1*K2/(K1+K2)

      VELS = SQRT( MUXZ(JM1) *H*DENW(JM1)/DT )
      VELP = SQRT( KT        *H*DENW(JM1)/DT )
      VEL  = VELS
      CALL NONREF ( H, KTRI, VEL, VELS, VELP,                          &
                                AYW(K,L)%ARI01, AYW(K,L)%ARI02,        &
                AYW(K,L)%ARI10, AYW(K,L)%ARI11, AYW(K,L)%ARI12,        &
                AYW(K,L)%ARI20, AYW(K,L)%ARI21, AYW(K,L)%ARI22   )

     ELSE

      AYU (K,L)%ARI01 = -1.
      AYU (K,L)%ARI02 =  0.
      AYU (K,L)%ARI10 =  0.
      AYU (K,L)%ARI11 =  0.
      AYU (K,L)%ARI12 =  0.
      AYU (K,L)%ARI20 =  0.
      AYU (K,L)%ARI21 =  0.
      AYU (K,L)%ARI22 =  0.

      AYV (K,L)%ARI01 =  0.
      AYV (K,L)%ARI02 =  1.
      AYV (K,L)%ARI10 =  0.
      AYV (K,L)%ARI11 =  0.
      AYV (K,L)%ARI12 =  0.
      AYV (K,L)%ARI20 =  0.
      AYV (K,L)%ARI21 =  0.
      AYV (K,L)%ARI22 =  0.

      AYW (K,L)%ARI01 =  0.
      AYW (K,L)%ARI02 =  1.
      AYW (K,L)%ARI10 =  0.
      AYW (K,L)%ARI11 =  0.
      AYW (K,L)%ARI12 =  0.
      AYW (K,L)%ARI20 =  0.
      AYW (K,L)%ARI21 =  0.
      AYW (K,L)%ARI22 =  0.

     END IF

    END DO
    IF ( Z_MAX >= MZ ) THEN
      DO K = Y_MIN, Y_MAX
        AYU(K,MZ:Z_MAX) = AYU(K,MZ-1)
        AYV(K,MZ:Z_MAX) = AYV(K,MZ-1)
        AYW(K,MZ:Z_MAX) = AYW(K,MZ-1)
      END DO
    END IF
   END IF
  END DO

!____________________________________________________ BOTTOM TRANSPARENT

  IF ( Z_MAX >= MZ ) THEN
    JM1 = JM ( X_MIN, Y_MIN, MZ   )
    JM2 = JM ( X_MIN, Y_MIN, MZ-1 )
    VELS = SQRT( MUXZ(JM1) *H*DENU(JM1)/DT )
    VELP = SQRT( L2MZ(JM2) *H*DENW(JM1)/DT )
    VEL  = VELP
    CALL NONREF ( H,KTBO, VEL, VELS, VELP,                             &
                          ABW01, ABW02,                                &
                   ABW10, ABW11, ABW12,                                &
                   ABW20, ABW21, ABW22  )

    VEL  = VELS
    CALL NONREF ( H,KTBO, VEL, VELS, VELP,                             &
                          ABH01, ABH02,                                &
                   ABH10, ABH11, ABH12,                                &
                   ABH20, ABH21, ABH22  )
  END IF


END SUBROUTINE PRE_NR_COEF1

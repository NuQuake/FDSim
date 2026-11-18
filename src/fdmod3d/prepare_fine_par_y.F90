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
SUBROUTINE PREPARE_FINE_PAR_Y (K1, MXT, NK, L, H, Z0)

  USE NRTYPE              , ONLY: SP, WP
  USE MOD_CONST           , ONLY: N_FREQ, I2N_FREQ1
  USE MOD_MODEL           , ONLY: XBMIN, YBMIN
  USE MOD_ATTEN           , ONLY: FRK, FRJ, PFP21, PF21, CF
  USE MOD_INTEGRATION     , ONLY: PTS, PTS2,                           &
                                  MU_FINE, LA_FINE, YMU_FINE, YLA_FINE
  USE MOD_PARAMETERIZATION, ONLY: FREF
  USE MOD_INTERFACES      , ONLY: FUNC_YMU_YL2M

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER , INTENT(IN) :: K1, MXT, NK, L
  REAL(SP), INTENT(IN) :: H, Z0

  INTEGER              :: II, KI, LI, J

  REAL(WP)             :: HPTS_H, H2_H, XSTART, YSTART, ZSTART,        &
                          POFRJ, POFRJ2P1,                             &
                          THE11, THE12, R1,                            &
                          THE21, THE22, R2,                            &
                          MU, L2M, X, Y, Z

  REAL(WP)   , DIMENSION(  N_FREQ, 2 ) :: YY

  COMPLEX(WP), DIMENSION(I2N_FREQ1)    :: SUM1, SUM2

!-----------------------------------------------------------------------


  HPTS_H  =  H/REAL(PTS,WP)
  H2_H    = (H - HPTS_H)/2._WP

  ZSTART = Z0              - REAL(L -1,WP)*H + H2_H
  YSTART = YBMIN - H/2._WP + REAL(K1-1,WP)*H - H2_H
  XSTART = XBMIN - H/2._WP                   - H2_H

  DO LI = PTS+PTS2, PTS2+1, -1
    Z = ZSTART - REAL(LI-1,WP)*HPTS_H

    DO KI = 1, NK*PTS+PTS2
      Y = YSTART + REAL(KI-1,WP)*HPTS_H

      DO II = 1, MXT*PTS+PTS2
        X = XSTART + REAL(II-1,WP)*HPTS_H

        YY  = FUNC_YMU_YL2M (X,Y,Z)

        SUM1 = CMPLX(0.,0.)
        SUM2 = CMPLX(0.,0.)
        THE11 = 0._WP
        THE12 = 0._WP
        THE21 = 0._WP
        THE22 = 0._WP
        DO J = 1, N_FREQ

          POFRJ    = FREF / FRJ(J)
          POFRJ2P1 = POFRJ * POFRJ + 1._WP

          SUM1(:)  = SUM1(:) + YY(J,1)*CF(J,:)

          SUM2(:)  = SUM2(:) + YY(J,2)*CF(J,:)

          THE11    = THE11 + YY(J,1) / PF21 (J)
          THE12    = THE12 + YY(J,1) * PFP21(J)
          THE21    = THE21 + YY(J,2) / PF21 (J)
          THE22    = THE22 + YY(J,2) * PFP21(J)

        END DO

        THE11 = 1._WP - THE11
        R1    = SQRT( THE11*THE11 + THE12*THE12 )

        THE21 = 1._WP - THE21
        R2    = SQRT( THE21*THE21 + THE22*THE22 )

        MU  = MU_FINE(II,KI,LI)                                        &
            * ((R1+THE11)/(2._WP*R1*R1))

        L2M = ( LA_FINE(II,KI,LI) + 2.*MU_FINE(II,KI,LI) )             &
            * ((R2+THE21)/(2._WP*R2*R2))

        YMU_FINE (II,KI,LI,:) =   MU * (1. - SUM1(:))

        YLA_FINE (II,KI,LI,:) =   L2M* (1. - SUM2(:))                  &
                              -2.* MU* (1. - SUM1(:))

         MU_FINE (II,KI,LI  ) = MU
         LA_FINE (II,KI,LI  ) = L2M - 2.*MU

      END DO
    END DO
  END DO

  IF ( L == 1 ) THEN

    DO LI = 1, PTS2
      YMU_FINE (:,:,LI,:) = YMU_FINE (:,:,PTS+1-LI,:)
      YLA_FINE (:,:,LI,:) = YLA_FINE (:,:,PTS+1-LI,:)
       MU_FINE (:,:,LI  ) =  MU_FINE (:,:,PTS+1-LI  )
       LA_FINE (:,:,LI  ) =  LA_FINE (:,:,PTS+1-LI  )
    END DO

  ELSE

    DO LI = 1, PTS2
      Z = ZSTART - REAL(LI-1,WP)*HPTS_H

      DO KI = 1, NK*PTS+PTS2
        Y = YSTART + REAL(KI-1,WP)*HPTS_H

        DO II = 1, MXT*PTS+PTS2
          X = XSTART + REAL(II-1,WP)*HPTS_H

          YY  = FUNC_YMU_YL2M (X,Y,Z)

          SUM1 = CMPLX(0.,0.)
          SUM2 = CMPLX(0.,0.)
          THE11 = 0._WP
          THE12 = 0._WP
          THE21 = 0._WP
          THE22 = 0._WP
          DO J = 1, N_FREQ

            SUM1(:)  = SUM1(:) + YY(J,1)*CF(J,:)

            SUM2(:)  = SUM2(:) + YY(J,2)*CF(J,:)

            THE11    = THE11 + YY(J,1) / PF21 (J)
            THE12    = THE12 + YY(J,1) * PFP21(J)
            THE21    = THE21 + YY(J,2) / PF21 (J)
            THE22    = THE22 + YY(J,2) * PFP21(J)

          END DO

          THE11 = 1._WP - THE11
          R1    = SQRT( THE11*THE11 + THE12*THE12 )

          THE21 = 1._WP - THE21
          R2    = SQRT( THE21*THE21 + THE22*THE22 )

          MU  = MU_FINE(II,KI,LI)                                      &
              * ((R1+THE11)/(2._WP*R1*R1))

          L2M = ( LA_FINE(II,KI,LI) + 2.*MU_FINE(II,KI,LI) )           &
              * ((R2+THE21)/(2._WP*R2*R2))

          YMU_FINE (II,KI,LI,:) =     MU * (1. - SUM1(:))

          YLA_FINE (II,KI,LI,:) =    L2M * (1. - SUM2(:))              &
                                - 2.* MU * (1. - SUM1(:))

           MU_FINE (II,KI,LI  ) = MU
           LA_FINE (II,KI,LI  ) = L2M - 2.*MU

        END DO
      END DO
    END DO

  END IF

END SUBROUTINE PREPARE_FINE_PAR_Y

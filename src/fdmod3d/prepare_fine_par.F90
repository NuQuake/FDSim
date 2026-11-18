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
SUBROUTINE PREPARE_FINE_PAR (K1, MXT, NK, L, H, Z0)

  USE NRTYPE         , ONLY: WP, SP
  USE MOD_MODEL      , ONLY: XBMIN, YBMIN
  USE MOD_INTEGRATION, ONLY: PTS, PTS2, MU_FINE, LA_FINE, RHO_FINE
  USE MOD_INTERFACES , ONLY: FUNC_MU_LA_RHO
!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER , INTENT(IN) :: K1, MXT, NK, L
  REAL(SP), INTENT(IN) :: H, Z0

  INTEGER              :: II, KI, LI

  REAL(WP)             :: HPTS_H, H2_H, XSTART, YSTART, ZSTART,        &
                          PAR(3), X, Y, Z

!----------------------------------------------------------------------

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

        PAR = FUNC_MU_LA_RHO (X,Y,Z)
         MU_FINE (II,KI,LI) = PAR(1)
         LA_FINE (II,KI,LI) = PAR(2)
        RHO_FINE (II,KI,LI) = PAR(3)

      END DO
    END DO
  END DO

  IF ( L == 1 ) THEN

    DO LI = 1, PTS2
       MU_FINE (:,:,LI) =  MU_FINE (:,:,PTS+1-LI)
       LA_FINE (:,:,LI) =  LA_FINE (:,:,PTS+1-LI)
      RHO_FINE (:,:,LI) = RHO_FINE (:,:,PTS+1-LI)
    END DO

  ELSE

    DO LI = 1, PTS2
      Z = ZSTART - REAL(LI-1,WP)*HPTS_H

      DO KI = 1, NK*PTS+PTS2
        Y = YSTART + REAL(KI-1,WP)*HPTS_H

        DO II = 1, MXT*PTS+PTS2
          X = XSTART + REAL(II-1,WP)*HPTS_H

          PAR = FUNC_MU_LA_RHO (X,Y,Z)
           MU_FINE (II,KI,LI) = PAR(1)
           LA_FINE (II,KI,LI) = PAR(2)
          RHO_FINE (II,KI,LI) = PAR(3)

        END DO
      END DO
    END DO

  END IF

END SUBROUTINE PREPARE_FINE_PAR

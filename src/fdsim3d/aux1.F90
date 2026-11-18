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

SUBROUTINE  AUX1
  USE PRECISION   , ONLY :                                             &
                           PP

  USE GRID_MEDIUM , ONLY :                                             &
                           H

  USE CONTROL_DATA, ONLY :                                             &
                           MX, MY, MZ, DT, FRJMAX, FRANGE, TPML, LPAS, &
                           REFINE, N_FREQ

  USE AUXIL       , ONLY :                                             &
                           PI,                                         &
                           A, B, AH, BH, CH, H32, H2,                  &
                           OTA, OTB, OTC, OTD,                         &
                           PMLL, PMLXH , PMLYH , PMLZH, PMLZH1, PMLZH2,&
                           MX1, MY1,                                   &
                           AHF, BHF, CHF, MXF, MYF, MXF1, MYF1,  HF,   &
                           PMLLF, PMLXHF, PMLYHF, TPMLF, HF3, H3

#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER        :: J

  REAL (KIND=PP) :: OH

!-----------------------------------------------------------------------

  MX1 = MX - 1
  MY1 = MY - 1

  PMLL   =   1 - TPML
  PMLXH  = MX  + TPML
  PMLYH  = MY  + TPML
  PMLZH  = MZ  + TPML
  PMLZH1 = MZ  + TPML - 1
  PMLZH2 = MZ  + TPML - 2

  AH  = A    /H
  BH  = B    /H
  CH  = 1._PP/H

  H32 = H*H*H/2._PP
  H2  = H*H

  HF  = H

  IF ( LPAS /= 0 ) THEN
    TPMLF= TPML * REFINE
    HF   = H    / REFINE

    MXF  = (MX-1)*REFINE+1
    MYF  = (MY-1)*REFINE+1

    MXF1 = MXF - 1
    MYF1 = MYF - 1

    PMLLF   =   1 - TPMLF
    PMLXHF  = MXF + TPMLF
    PMLYHF  = MYF + TPMLF

    AHF  = A    /HF
    BHF  = B    /HF
    CHF  = 1._PP/HF

  END IF

  H3  = H *H *H
  HF3 = HF*HF*HF

  OTA = 0.
  OTB = 1.
  OTC = 1.
  OTD = 0.

  DO J = 1, N_FREQ
    OH = FRJMAX/(10._PP**(FRANGE*REAL(N_FREQ-J,PP)/REAL(N_FREQ-1)))    &
       * DT *2._PP*PI
    OTA (J)= 2._PP*OH /(2._PP+OH)
    OTB (J)=(2._PP-OH)/(2._PP+OH)
    OTC (J)= 2._PP    /(2._PP-OH)
    OTD (J)=-      OH /(2._PP-OH)
  END DO


END SUBROUTINE  AUX1

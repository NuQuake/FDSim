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

SUBROUTINE  ENERGY(ITILE)

  USE PRECISION   , ONLY: PP

  USE CONTROL_DATA, ONLY: MX, MY, MZ, LPAS, DT

  USE GRID_MEDIUM , ONLY:                                              &
                          JM , JMF, H,                                 &
                          DENUF, DENVF, DENWF,                         &
                          DENU , DENV , DENW ,                         &
                         CL2MX ,CL2MY ,CL2MZ ,CLAXY ,CLAXZ ,CLAYZ ,    &
                         CMUXY ,CMUXZ ,CMUYZ
  USE AUXIL       , ONLY: X_MIN ,X_MAX , Y_MIN ,Y_MAX , Z_MIN ,Z_MAX , &
                          X_MINF,X_MAXF, Y_MINF,Y_MAXF, Z_MINF,Z_MAXF, &
                          MXF, MYF, HF, HF3, H3

  USE WAVEFIELD   , ONLY: UM, VM, WM, UMF, VMF, WMF, TOTAL_ENERGY

  USE STRESSFIELD , ONLY:               TXX , TYY , TXY ,              &
                                        TXZ , TYZ , TZZ ,              &
                                        TXXF, TYYF, TXYF,              &
                                        TXZF, TYZF, TZZF


#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ITILE

  INTEGER  :: I, K, L, JM1

  REAL(PP) :: T_E, DEN

!-----------------------------------------------------------------------

  IF ( LPAS /= 0 ) THEN

    DO L = MAX(0,Z_MINF), MIN(LPAS-1,Z_MAXF)
      DO K = MAX(3,Y_MINF), MIN(MYF-1,Y_MAXF)
        DO I = MAX(3,X_MINF), MIN(MXF-1,X_MAXF)

          JM1 = JMF (I,K,L)

          DEN = DT/HF/DENUF(JM1)
          TOTAL_ENERGY = TOTAL_ENERGY + DEN*UMF(I,K,L)*UMF(I,K,L)*HF3

          DEN = DT/HF/DENVF(JM1)
          TOTAL_ENERGY = TOTAL_ENERGY + DEN*VMF(I,K,L)*VMF(I,K,L)*HF3

          DEN = DT/HF/DENWF(JM1)
          TOTAL_ENERGY = TOTAL_ENERGY + DEN*WMF(I,K,L)*WMF(I,K,L)*HF3

          TOTAL_ENERGY = TOTAL_ENERGY                                  &
                       + TXYF(I,K,L)*TXYF(I,K,L)*CMUXY(JM1)*HF3

          TOTAL_ENERGY = TOTAL_ENERGY                                  &
                       + TXZF(I,K,L)*TXZF(I,K,L)*CMUXZ(JM1)*HF3

          TOTAL_ENERGY = TOTAL_ENERGY                                  &
                       + TYZF(I,K,L)*TYZF(I,K,L)*CMUYZ(JM1)*HF3

          TOTAL_ENERGY = TOTAL_ENERGY                                  &
                       + TXXF(I,K,L)*TXXF(I,K,L)*CL2MX(JM1)*HF3        &
                       + TYYF(I,K,L)*TYYF(I,K,L)*CLAXY(JM1)*HF3        &
                       + TZZF(I,K,L)*TZZF(I,K,L)*CLAXZ(JM1)*HF3

          TOTAL_ENERGY = TOTAL_ENERGY                                  &
                       + TXXF(I,K,L)*TXXF(I,K,L)*CLAXY(JM1)*HF3        &
                       + TYYF(I,K,L)*TYYF(I,K,L)*CL2MY(JM1)*HF3        &
                       + TZZF(I,K,L)*TZZF(I,K,L)*CLAYZ(JM1)*HF3

          TOTAL_ENERGY = TOTAL_ENERGY                                  &
                       + TXXF(I,K,L)*TXXF(I,K,L)*CLAXZ(JM1)*HF3        &
                       + TYYF(I,K,L)*TYYF(I,K,L)*CLAYZ(JM1)*HF3        &
                       + TZZF(I,K,L)*TZZF(I,K,L)*CL2MZ(JM1)*HF3


        END DO
      END DO
    END DO

  ENDIF

  DO L = MAX(LPAS,Z_MIN ), MIN(MZ -1,Z_MAX )
    DO K = MAX(   3,Y_MIN ), MIN(MY -1,Y_MAX )
      DO I = MAX(   3,X_MIN ), MIN(MX -1,X_MAX )

        JM1 = JM (I,K,L)

        DEN = DT/H /DENU (JM1)
        TOTAL_ENERGY = TOTAL_ENERGY + DEN*UM (I,K,L)*UM (I,K,L)*H3

        DEN = DT/H /DENV (JM1)
        TOTAL_ENERGY = TOTAL_ENERGY + DEN*VM (I,K,L)*VM (I,K,L)*H3

        DEN = DT/H /DENW (JM1)
        TOTAL_ENERGY = TOTAL_ENERGY + DEN*WM (I,K,L)*WM (I,K,L)*H3

        TOTAL_ENERGY = TOTAL_ENERGY                                    &
                     + TXY (I,K,L)*TXY (I,K,L)*CMUXY(JM1)*H3

        TOTAL_ENERGY = TOTAL_ENERGY                                    &
                     + TXZ (I,K,L)*TXZ (I,K,L)*CMUXZ(JM1)*H3

        TOTAL_ENERGY = TOTAL_ENERGY                                    &
                     + TYZ (I,K,L)*TYZ (I,K,L)*CMUYZ(JM1)*H3

        TOTAL_ENERGY = TOTAL_ENERGY                                    &
                     + TXX (I,K,L)*TXX (I,K,L)*CL2MX(JM1)*H3           &
                     + TYY (I,K,L)*TYY (I,K,L)*CLAXY(JM1)*H3           &
                     + TZZ (I,K,L)*TZZ (I,K,L)*CLAXZ(JM1)*H3

        TOTAL_ENERGY = TOTAL_ENERGY                                    &
                     + TXX (I,K,L)*TXX (I,K,L)*CLAXY(JM1)*H3           &
                     + TYY (I,K,L)*TYY (I,K,L)*CL2MY(JM1)*H3           &
                     + TZZ (I,K,L)*TZZ (I,K,L)*CLAYZ(JM1)*H3

        TOTAL_ENERGY = TOTAL_ENERGY                                    &
                     + TXX (I,K,L)*TXX (I,K,L)*CLAXZ(JM1)*H3           &
                     + TYY (I,K,L)*TYY (I,K,L)*CLAYZ(JM1)*H3           &
                     + TZZ (I,K,L)*TZZ (I,K,L)*CL2MZ(JM1)*H3

      END DO
    END DO
  END DO


#if(USE_MPI)
  CALL MPI_REDUCE                                                      &
     (TOTAL_ENERGY, T_E, 1, MPI_REAL, MPI_SUM, root, comm_cart, errcode)

  IF ( myid == root ) THEN
    OPEN  ( 21, FILE='ENERGY.DAT', STATUS='UNKNOWN', POSITION='APPEND')
    WRITE ( 21, '(2(1X,E13.6))') (ITILE-1)*DT, T_E/2.
    CLOSE ( 21)
  END IF
#else
  OPEN  ( 21, FILE='ENERGY.DAT', STATUS='UNKNOWN', POSITION='APPEND')
  WRITE ( 21, '(2(1X,E13.6))')  (ITILE-1)*DT, TOTAL_ENERGY/2.
  CLOSE ( 21)
#endif


END SUBROUTINE  ENERGY

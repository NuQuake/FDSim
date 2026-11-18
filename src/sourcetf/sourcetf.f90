!=======================================================================
! Copyright (c) 2014 Jozef Kristek, Peter Moczo
!
! Licensed under the EUPL (the "Licence")
! You may not use this work except in compliance with the Licence.
! You may obtain a copy of the Licence at
! https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12
! Unless required by applicable law or agreed to in writing, software
! distributed under the Licence is distributed on an "AS IS" basis,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and
! limitations under the Licence.
!
! REVISION HISTORY:
! Version 1.0, 2014-03-19: Jozef Kristek, Peter Moczo
!=======================================================================
PROGRAM SOURCETF

  USE PRECISION, ONLY : WP

!----------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER*20 :: SRCBUF

  INTEGER :: I, J, K, NPS, NSIG, NT, ORDER, KN, NF, NF1, JSTART

  REAL (KIND=WP)            :: TP, TS, GAMA, FP, PSI, ZETA, TRIAN,     &
                               A1, A2, A3, DT, T, AA, BB, CC, T00,     &
                               FMIN, FMAX, BZERO, PHASE, DF, TMP, MESTF

  REAL, DIMENSION(:), ALLOCATABLE :: FA1, FA2, STF, ESTF

  REAL (KIND=WP), PARAMETER :: PI=3.141592654

  COMPLEX, DIMENSION(131072):: SSTF, CSTF, CESTF


!----------------------------------------------------------------------


  NAMELIST /INPUT/        NSIG, DT
  NAMELIST /SIGNAL_1/     TP
  NAMELIST /SIGNAL_2/     TP, TS
  NAMELIST /SIGNAL_3/     GAMA, FP, PSI , TS
  NAMELIST /SIGNAL_4/     GAMA, FP, ZETA, TS

!----------------------------------------------------------------------

  OPEN (10, FILE='SOURCETF.IN' , STATUS='OLD')
  OPEN (11, FILE='SIGNAL.DAT'  , STATUS='NEW')
  OPEN (12, FILE='ENVELOPE.DAT', STATUS='NEW')
  OPEN (13, FILE='SPECTRUM.DAT', STATUS='NEW')

  KN = 14
  READ ( 10, NML = INPUT )

  SELECT CASE ( NSIG )
  CASE ( 1 )
    READ  (10, NML=SIGNAL_1)
  CASE ( 2 )
    READ  (10, NML=SIGNAL_2)
  CASE ( 3 )
    READ  (10, NML=SIGNAL_3)
  CASE ( 4 )
    READ  (10, NML=SIGNAL_4)
  END SELECT

  SELECT CASE ( NSIG )
  !_____________________________________________ KUPPER ________
    CASE ( 1 )
      A1 = 2._WP*PI/TP
      A2 = 2._WP*A1
      NT = TP/DT + 1._WP
  !_____________________________________________ RICKER ________
    CASE ( 2 )
      A1 = TP*SQRT(6._WP)/PI
      A1 = A1*A1
      A2 = SQRT(PI)/2._WP
      IF  ( TS < .0001_WP ) TS=1.1_WP*TP
  !_____________________________________________ GABOR _________
    CASE ( 3 )
      A1 = 2._WP*PI*FP
      A2 = A1/GAMA
      A2 = A2*A2
      A3 = PSI
      IF  ( TS < .0001_WP ) TS=0.45_WP*GAMA/FP
  !_____________________________________________ BERLAGE _______
    CASE ( 4 )
      A1 = 2._WP*PI*FP
      A2 = A1/GAMA
      A3 = ZETA
  
    END SELECT

    IF  ( NSIG /= 1 )  THEN
      NT = TS/DT +.5_WP
      TS = DT*NT
      NT = NT+1
      IF  ( NSIG /= 4 )  THEN
        NT  = 2*NT -1
      END IF
    END IF

    ALLOCATE ( STF (NT), ESTF (NT) )

    DO J = 1, NT
      T = (J-1)*DT

      SELECT CASE ( NSIG )
!_______________________________________ KUPPER _________
      CASE ( 1 )
        AA = T*A1
        BB = T*A2
        STF(J) = SIN(AA) - .5_WP*SIN(BB)
!_______________________________________ RICKER _________
      CASE ( 2 )
        AA = T - TS
        AA = (6._WP*AA*AA)/A1
        STF(J) = A2*( AA - .5_WP )*EXP(-AA)
!_______________________________________ GABOR __________
      CASE ( 3 )
        AA = T - TS
        STF(J) = EXP(-A2*AA*AA)*COS(A1*AA + A3)
!_______________________________________ BERLAGE ________
      CASE ( 4 )
        STF(J) =(T**A3)*EXP(-A2*T)*SIN(A1*T)
    
      END SELECT

    END DO

    NF  = 2**KN
    DF  = 1./(DT*REAL(NF,WP))
    NF1 = NF/2+1

    DO J = 1, NT
      SSTF(J) = CMPLX(STF(J),0.)
    END DO
    DO J = NT+1,NF
      SSTF(J)=(0.,0.)
    END DO

    CALL FCOOLR(KN,SSTF,-1.)

    SSTF = SSTF*DT

    DO J = NF1+1, NF
      SSTF(J)=CONJG (SSTF(NF-J+2))
    END DO

    CSTF = SSTF

    CALL FCOOLR(KN,CSTF, 1.)

    DO J = 1, NF1
      CESTF (J) =   SSTF (J) * CMPLX (0.,1.)
    END DO

    DO J = NF1 + 1, NF
      CESTF (J) = - SSTF (J) * CMPLX (0.,1.)
    END DO

    CALL FCOOLR ( KN, CESTF, 1. )

    DO J = 1, NT
      TMP      = REAL ( CESTF(J) ) * DF
      ESTF (J) = SQRT ( STF(J)*STF(J) + TMP*TMP )
    END DO

    DO J = 1, NT

      WRITE ( 11, * ) (J-1)*DT,STF(J)
      WRITE ( 12, * ) (J-1)*DT,ESTF(J)

    END DO

    CLOSE ( 12 )

    DO J = 1, NF1
      IF ( ABS(SSTF(J)) > 1.e-30 ) THEN
        PHASE = ATAN(IMAG(SSTF(J))/REAL(SSTF(J)))
      ELSE
        PHASE = 0.
      END IF
      WRITE ( 13, '(4(E15.4,2X))') (J-1)*DF, ABS(SSTF(J)*SSTF(J)),     &
                            ABS(SSTF(J)), PHASE
    END DO


  CLOSE ( 10 )
  CLOSE ( 11 )
  CLOSE ( 13 )

END PROGRAM SOURCETF

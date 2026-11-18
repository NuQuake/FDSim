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
! 2026-11-18: spline interpolation of source time functions reimplemented
!=======================================================================

SUBROUTINE  OPEN_SR_FILE

  USE PRECISION   , ONLY:                                              &
                          PP, DP

  USE CONTROL_DATA, ONLY:                                              &
                          SR_FILE_NAME, KEY_ONE_STF,                   &
                          MX, MY,       NPS, MT1, MT2,                 &
                          SFIS, SDEL, SLAM,  M0, DT,                   &
                          IS, KS, LS, ITSB, ITSE, PWINC,               &
                          LPAS, REFINE, LIN, CU, CV,                   &
                          STRESS_IMAGING, XBMIN, YBMIN

  USE GRID_MEDIUM , ONLY: JM, MUXZ, L2MZ, DENU, DENW, H

  USE AUXIL       , ONLY:                                              &
                          SXXT, SXYT, SXZT, SYYT, SYZT, SZZT,          &
                          PI, SOURTF_A,                                &
                          CUS, CVS, CWS, SRCDELI, SRC_P, NPST,         &
                          X_MIN, Y_MIN, Z_MIN, Z_MAX

#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif

!-----------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER (LEN = 12) :: SRCBUF

  LOGICAL       ::  EX

  INTEGER       ::  I, ALLOSTAT, NSTFS, SRC, IOLEN, IOS, NSTFSP, JM1

  REAL (KIND=PP)::  SXX, SXY, SXZ, SYY, SYZ, SZZ, SOURTF, VEL, VELX, T,&
                    TMP, XS, YS, ZS, TB, TE, ANGLE, DPW

  REAL (KIND=PP), DIMENSION(:), ALLOCATABLE ::  TIME, STF, STF2

#if(NOISE)
  NAMELIST /SOURCE/       NPS
#else
  NAMELIST /SOURCE/       NPS, KEY_ONE_STF, PWINC
#endif

!-----------------------------------------------------------------------

  INTERFACE

    SUBROUTINE  FDL2MT (SFIS,SDEL,SLAM,M0,SXX,SXY,SXZ,SYY,SYZ,SZZ)
      USE PRECISION
      REAL (KIND=PP), INTENT(IN)  :: SFIS, SDEL, SLAM, M0
      REAL (KIND=PP), INTENT(OUT) :: SXX, SXY, SXZ, SYY, SYZ, SZZ
    END SUBROUTINE FDL2MT

  END INTERFACE

!______________________________________________ SOURCE FILE _ 10 ______

#if(USE_MPI)
  IF ( myid == root ) THEN
#endif

    OPEN  (10, FILE = SR_FILE_NAME, STATUS='OLD' )

    PWINC = 0
    KEY_ONE_STF = .FALSE.
    READ  (10, NML=SOURCE)
    WRITE (11, NML=SOURCE)

    IF ( PWINC > 0 ) THEN
      NPS = 2
      KEY_ONE_STF = .FALSE.
      IF ( LPAS /= 0 ) THEN
        WRITE (11,*) ' PLANE WAVE INCIDENCE IS NOT IMPLEMENTED IN DISCONTINUOUS GRID !'
        STOP
      END IF
    END IF

#if(USE_MPI)
  END IF

  CALL MPI_Barrier(comm_cart,errcode)

  CALL MPI_Bcast  (NPS        , 1, MPI_INTEGER, root, comm_cart, errcode)
  CALL MPI_Bcast  (PWINC      , 1, MPI_INTEGER, root, comm_cart, errcode)
  CALL MPI_Bcast  (KEY_ONE_STF, 1, MPI_LOGICAL, root, comm_cart, errcode)
  CALL MPI_Barrier(comm_cart,errcode)
#endif

#if(NOISE)
  ALLOCATE ( IS  (1:NPS ), KS (1:NPS ), LS (1:NPS ),                   &
             ITSB(1:NPS ),ITSE(1:NPS ),                                &
             CUS (1:NPS ),CVS (1:NPS ),CWS (1:NPS ),   STAT = ALLOSTAT )
  ALLOCATE (SRCDELI(MT2))
#else
  ALLOCATE ( IS  (1:NPS ), KS (1:NPS ), LS (1:NPS ),                   &
             ITSB(1:NPS ),ITSE(1:NPS ),                                &
             SXXT(1:NPS ),SXYT(1:NPS ),SXZT(1:NPS ),                   &
             SYYT(1:NPS ),SYZT(1:NPS ),SZZT(1:NPS ),   STAT = ALLOSTAT )
#endif

    IF  ( ALLOSTAT > 0 )  THEN
      WRITE (11,*) '  ALLOCATION ERROR O_S_F_01, STAT=', ALLOSTAT
      STOP
    END IF

#if(USE_MPI)
  IF ( myid == root ) THEN
#endif

#if(NOISE)
    DO I = 1, NPS
      READ (10,*) IS(I), KS(I), LS(I), ITSB(I), ITSE(I)
      READ (10,*) CUS(I), CVS(I), CWS(I)
    END DO
#else
    DO I = 1, NPS
      IF ( PWINC > 0 ) THEN
        IF ( I == 1 ) THEN
          READ (10,*) DPW
          LIN = NINT( ABS(DPW) / H )
          IF ( PWINC == 1 ) THEN
            READ (10,*) ANGLE
            CU = COS(ANGLE*PI/180.)
            CV = SIN(ANGLE*PI/180.)
          END IF
        END IF
      ELSE
        READ (10,*) XS, YS, ZS, TB, TE
        READ (10,*) SFIS, SDEL, SLAM, M0

        ITSB(I) = TB/DT + 1
        ITSE(I) = TE/DT + 1

        ZS = ABS(ZS)

        IF ( LPAS > 0 ) THEN
          IF ( ZS < 1.5*H/REFINE ) ZS = 1.5*H/REFINE
        ELSE
          IF ( ZS < 1.5*H        ) ZS = 1.5*H
        END IF

        IF (STRESS_IMAGING) THEN
          IF ( LPAS > 0 ) THEN
            ZS = ZS + 2*H/REFINE
          ELSE
            ZS = ZS + 2*H
          END IF
        END IF

        IF ( ZS >= LPAS*H/REFINE ) THEN
          KS(I) = NINT((XS - XBMIN) / H + 1 )
          IS(I) = NINT((YS - YBMIN) / H + 1 )
          LS(I) = NINT((ZS-LPAS*H/REFINE) / H + LPAS - 0.5 )
        ELSE
          KS(I) = NINT((XS - XBMIN) / H * REFINE + 1 )
          IS(I) = NINT((YS - YBMIN) / H * REFINE + 1 )
          LS(I) = NINT( ZS          / H * REFINE - 0.5 )
        END IF

        CALL FDL2MT(SFIS,SDEL,SLAM,M0, SXX,SXY,SXZ,SYY,SYZ,SZZ)
          SXXT(I) = SXX
          SXYT(I) = SXY
          SXZT(I) = SXZ
          SYYT(I) = SYY
          SYZT(I) = SYZ
          SZZT(I) = SZZ
      ENDIF
    END DO
#endif


    IF (STRESS_IMAGING) THEN
      LS  = LS  + 2
      LIN = LIN + 2
    END IF

    CLOSE(10)

#if(USE_MPI)
  END IF

  CALL MPI_Barrier(comm_cart,errcode)
  CALL MPI_Bcast  (IS  , NPS, MPI_INTEGER, root, comm_cart, errcode)
  CALL MPI_Bcast  (KS  , NPS, MPI_INTEGER, root, comm_cart, errcode)
  CALL MPI_Bcast  (LS  , NPS, MPI_INTEGER, root, comm_cart, errcode)
  CALL MPI_Bcast  (ITSB, NPS, MPI_INTEGER, root, comm_cart, errcode)
  CALL MPI_Bcast  (ITSE, NPS, MPI_INTEGER, root, comm_cart, errcode)
#if(NOISE)
  CALL MPI_Bcast  (CUS , NPS, MPI_REAL   , root, comm_cart, errcode)
  CALL MPI_Bcast  (CVS , NPS, MPI_REAL   , root, comm_cart, errcode)
  CALL MPI_Bcast  (CWS , NPS, MPI_REAL   , root, comm_cart, errcode)
#else
  CALL MPI_Bcast  (CU  , 1  , MPI_REAL   , root, comm_cart, errcode)
  CALL MPI_Bcast  (CV  , 1  , MPI_REAL   , root, comm_cart, errcode)
  CALL MPI_Bcast  (LIN , 1  , MPI_INTEGER, root, comm_cart, errcode)
  CALL MPI_Bcast  (SXXT, NPS, MPI_REAL   , root, comm_cart, errcode)
  CALL MPI_Bcast  (SXYT, NPS, MPI_REAL   , root, comm_cart, errcode)
  CALL MPI_Bcast  (SXZT, NPS, MPI_REAL   , root, comm_cart, errcode)
  CALL MPI_Bcast  (SYYT, NPS, MPI_REAL   , root, comm_cart, errcode)
  CALL MPI_Bcast  (SYZT, NPS, MPI_REAL   , root, comm_cart, errcode)
  CALL MPI_Bcast  (SZZT, NPS, MPI_REAL   , root, comm_cart, errcode)
#endif
  CALL MPI_Barrier(comm_cart,errcode)
#endif

#if(NOISE)
  NSTFS = 0
  NPST  = 0
  DO SRC = 1, NPS
    IF ( ITSE(SRC) < 999999 ) THEN
      NSTFS = MAX(NSTFS,ITSE(SRC)-ITSB(SRC)+1)
      NPST  = NPST + 1
    END IF
  END DO
  ALLOCATE ( SRC_P(NPST) )
  I = 0
  DO SRC = 1, NPS
    IF ( ITSE(SRC) < 999999 ) THEN
      I = I + 1
      SRC_P(I) = SRC
    END IF
  END DO
#else
  NSTFS = MAXVAL(ITSE(:)-ITSB(:)+1)
  IF ( PWINC > 0) NSTFS = MT2

  IF ( KEY_ONE_STF ) THEN
    NPST = 1
  ELSE
    NPST = NPS
  END IF
#endif
  ALLOCATE( SOURTF_A(NPST,NSTFS) )
  SOURTF_A = 0.

  IF ( PWINC > 0 ) THEN
    IF ( ( Z_MIN <=LIN ) .AND. ( LIN <= Z_MAX ) ) THEN
        JM1 = JM (X_MIN,Y_MIN, LIN)
      IF ( PWINC == 1 ) THEN
        VELX = SQRT(MUXZ(JM1)*DENU(JM1)/DT*H)
      ELSE
        VELX = SQRT(L2MZ(JM1)*DENW(JM1)/DT*H)
      END IF
    ELSE
      VELX = 0.
    END IF

#if(USE_MPI)
    CALL MPI_REDUCE                                                      &
             (VELX ,VEL, 1, MPI_REAL, MPI_MAX, root, comm_cart, errcode)
#else
    VEL = VELX
#endif
  END IF

#if(USE_MPI)
  IF ( myid == root ) THEN
#endif

#if(NOISE)
    INQUIRE ( FILE = 'SRCDELI.DAT', EXIST=EX, IOSTAT=IOS )
    IF  ( EX )  THEN
      INQUIRE ( IOLENGTH = IOLEN ) SOURTF
      OPEN  ( 10, FILE = 'SRCDELI.DAT', FORM ='UNFORMATTED',             &
                ACCESS = 'DIRECT',   RECL = IOLEN, STATUS='OLD' )
      DO I = 1, MT2
        READ  ( 10, REC=I ) SRCDELI(I)
      END DO
      CLOSE(10)
    ELSE
      SRCDELI(:) = 0.
    END IF

    DO SRC = 1, NPST

      WRITE ( SRCBUF,'(A3, I5.5, A4)' )  'SRC_', SRC_P(SRC), '.DAT'
      INQUIRE ( IOLENGTH = IOLEN ) SOURTF
      OPEN  ( 10, FILE = SRCBUF(1:12), FORM ='UNFORMATTED',            &
                ACCESS = 'DIRECT',   RECL = IOLEN, STATUS='OLD' )
      DO I = 1, ITSE(SRC_P(SRC)) - ITSB(SRC_P(SRC))
        READ  ( 10, REC=I ) SOURTF_A(SRC,I)
      END DO
      CLOSE ( 10 )

    END DO
#else
    IF ( PWINC < 1 ) THEN

      DO SRC = 1, NPST

        WRITE ( SRCBUF,'(A4, I3.3, A4)' )  'SRC_', SRC, '.DAT'

        OPEN  ( 10, FILE = SRCBUF, STATUS='UNKNOWN' )

        I = 0
        DO
          READ(10,*,IOSTAT=IOS) TMP, TMP
          IF ( IOS < 0 ) EXIT
          I = I+1
        END DO
        REWIND (10)
        NSTFSP = I+2
        ALLOCATE ( TIME(NSTFSP), STF(NSTFSP), STF2(NSTFSP) )

        DO I = 2, NSTFSP-1
          READ(10,*) TIME(I), STF(I)
        END DO
        CLOSE(10)
        TIME(1    ) = -DT
        STF (1    ) = STF(2)
        IF ( TIME(NSTFSP-1) < (ITSE(SRC)-ITSB(SRC))*DT ) THEN
          TIME(NSTFSP) = (ITSE(SRC)-ITSB(SRC))*DT
        ELSE
          TIME(NSTFSP) = TIME(NSTFSP-1)+DT
        END IF
        STF (NSTFSP) = STF(NSTFSP-1)

        CALL SPLINE_CLAMPED_SETUP ( TIME, STF, 0._PP, 0._PP, STF2 )

        DO I = 1, ITSE(SRC)-ITSB(SRC)+1
          T = (I-1)*DT
          SOURTF_A (SRC,I) = SPLINE_EVAL (TIME,STF,STF2,T)
        END DO
        CLOSE(10)

        DEALLOCATE ( TIME, STF, STF2 )

      END DO

    ELSE
      OPEN  ( 10, FILE = 'PW_STF.DAT', STATUS='OLD' )

      I = 0
      DO
        READ(10,*,IOSTAT=IOS) TMP, TMP
        IF ( IOS < 0 ) EXIT
        I = I+1
      END DO
      REWIND (10)
      NSTFSP = I+2
      ALLOCATE ( TIME(NSTFSP), STF(NSTFSP), STF2(NSTFSP) )

      DO I = 2, NSTFSP-1
        READ(10,*) TIME(I), STF(I)
      END DO
      CLOSE(10)
      TIME(1    ) = -H/VEL
      STF (1    ) = STF(2)
      IF ( TIME(NSTFSP-1) < MT2*DT ) THEN
        TIME(NSTFSP) = MT2*DT
      ELSE
        TIME(NSTFSP) = TIME(NSTFSP-1)+H/VEL
      END IF
      STF (NSTFSP) = STF(NSTFSP-1)

      CALL SPLINE_CLAMPED_SETUP ( TIME, STF, 0._PP, 0._PP, STF2 )

      DO I = MT1, MT2
        T = (I-1)*DT
        SOURTF_A (1,I) = SPLINE_EVAL (TIME,STF,STF2,T)
        T = (I-1)*DT - H/VEL
        SOURTF_A (2,I) = SPLINE_EVAL (TIME,STF,STF2,T)
      END DO
    END IF

#endif

#if(USE_MPI)
  END IF

#if(NOISE)
  CALL MPI_Bcast  (SRCDELI ,        MT2, MPI_REAL, root, comm_cart,    &
                                                               errcode)
#endif
  CALL MPI_Bcast  (SOURTF_A, NSTFS*NPST, MPI_REAL, root, comm_cart,    &
                                                               errcode)
  CALL MPI_Barrier(comm_cart,errcode)
#endif

!=======================================================================

CONTAINS

!-----------------------------------------------------------------------
! Cubic spline with prescribed first derivatives at both ends.
!   X(1:N) strictly increasing abscissae, F(1:N) values,
!   DF1, DFN  first derivatives at X(1) and X(N).
! Returns S2(1:N), the second derivatives of the spline at the knots,
! from the continuity of the first derivative at interior knots:
!   h(i-1) S2(i-1) + 2 (h(i-1)+h(i)) S2(i) + h(i) S2(i+1)
!        = 6 [ (F(i+1)-F(i))/h(i) - (F(i)-F(i-1))/h(i-1) ]
! and the end conditions
!   2 h(1)   S2(1) + h(1)   S2(2) = 6 [ (F(2)-F(1))/h(1) - DF1 ]
!   h(N-1) S2(N-1) + 2 h(N-1) S2(N) = 6 [ DFN - (F(N)-F(N-1))/h(N-1) ]
! The tridiagonal system is solved by forward elimination and
! back substitution.
!-----------------------------------------------------------------------
  SUBROUTINE SPLINE_CLAMPED_SETUP ( X, F, DF1, DFN, S2 )
    USE PRECISION, ONLY: PP
    IMPLICIT NONE
    REAL(PP), DIMENSION(:), INTENT(IN)  :: X, F
    REAL(PP),               INTENT(IN)  :: DF1, DFN
    REAL(PP), DIMENSION(:), INTENT(OUT) :: S2
    REAL(PP), DIMENSION(SIZE(X)) :: SUB, DIAG, SUP, RHS
    REAL(PP) :: H_L, H_R, W
    INTEGER  :: N, I

    N = SIZE(X)
    IF ( N < 2 ) STOP 'SPLINE_CLAMPED_SETUP: at least 2 knots required'

    H_R      = X(2) - X(1)
    IF ( H_R <= 0. ) STOP 'SPLINE_CLAMPED_SETUP: abscissae not increasing'
    SUB (1)  = 0.
    DIAG(1)  = 2.*H_R
    SUP (1)  = H_R
    RHS (1)  = 6.*( (F(2)-F(1))/H_R - DF1 )

    DO I = 2, N-1
      H_L = X(I  ) - X(I-1)
      H_R = X(I+1) - X(I  )
      IF ( H_R <= 0. ) STOP 'SPLINE_CLAMPED_SETUP: abscissae not increasing'
      SUB (I) = H_L
      DIAG(I) = 2.*(H_L+H_R)
      SUP (I) = H_R
      RHS (I) = 6.*( (F(I+1)-F(I))/H_R - (F(I)-F(I-1))/H_L )
    END DO

    H_L      = X(N) - X(N-1)
    SUB (N)  = H_L
    DIAG(N)  = 2.*H_L
    SUP (N)  = 0.
    RHS (N)  = 6.*( DFN - (F(N)-F(N-1))/H_L )

    ! forward elimination
    DO I = 2, N
      W       = SUB(I) / DIAG(I-1)
      DIAG(I) = DIAG(I) - W*SUP(I-1)
      RHS (I) = RHS (I) - W*RHS(I-1)
    END DO

    ! back substitution
    S2(N) = RHS(N) / DIAG(N)
    DO I = N-1, 1, -1
      S2(I) = ( RHS(I) - SUP(I)*S2(I+1) ) / DIAG(I)
    END DO

  END SUBROUTINE SPLINE_CLAMPED_SETUP

!-----------------------------------------------------------------------
! Value at XI of the cubic spline defined by knots X, values F and
! second derivatives S2 (from SPLINE_CLAMPED_SETUP). Outside [X(1),X(N)]
! the end polynomial pieces are extended.
!-----------------------------------------------------------------------
  FUNCTION SPLINE_EVAL ( X, F, S2, XI )
    USE PRECISION, ONLY: PP
    IMPLICIT NONE
    REAL(PP), DIMENSION(:), INTENT(IN) :: X, F, S2
    REAL(PP),               INTENT(IN) :: XI
    REAL(PP) :: SPLINE_EVAL
    REAL(PP) :: H, A, B
    INTEGER  :: N, LO, HI, MID

    N = SIZE(X)

    ! bisection for the interval X(LO) <= XI < X(LO+1), LO in [1,N-1]
    LO = 1
    HI = N
    DO WHILE ( HI - LO > 1 )
      MID = (LO + HI) / 2
      IF ( XI >= X(MID) ) THEN
        LO = MID
      ELSE
        HI = MID
      END IF
    END DO
    HI = LO + 1

    H = X(HI) - X(LO)
    A = ( X(HI) - XI ) / H
    B = ( XI - X(LO) ) / H

    SPLINE_EVAL = A*F(LO) + B*F(HI)                                    &
                + ( (A*A*A-A)*S2(LO) + (B*B*B-B)*S2(HI) ) * H*H/6.

  END FUNCTION SPLINE_EVAL


END SUBROUTINE  OPEN_SR_FILE

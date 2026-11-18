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
SUBROUTINE READ_GBM3D

  USE MOD_MODEL, ONLY: X_MIN, Y_MIN, XIMIN, XIMAX, YIMIN, YIMAX,       &
                       RHO, VS, VP, QS

#if(USE_MPI)
  USE MPI_SHARED
#endif

  IMPLICIT NONE

  CHARACTER (LEN=2500) :: TMP_CHAR

  INTEGER :: IOS, I_MAX, K_MAX, L_MAX, NL_MAX, X_MAX, Y_MAX,           &
             X_MIN_M, X_MAX_M, Y_MIN_M, Y_MAX_M, XT, YT, I, K, L, NL

  REAL    :: ZT, ZB, RTMP


!=======================================================================

  OPEN ( 71, FILE = 'GBM3D_RHO.txt', STATUS = 'OLD' )
  OPEN ( 72, FILE = 'GBM3D_Vp.txt' , STATUS = 'OLD' )
  OPEN ( 73, FILE = 'GBM3D_Vs.txt' , STATUS = 'OLD' )
  OPEN ( 74, FILE = 'GBM3D_Qs.txt' , STATUS = 'OLD' )

  READ ( 71, '(A)' ) TMP_CHAR

   X_MIN_M = INT(XIMIN/10.       )*10
   X_MAX_M = INT(XIMAX/10.+0.9999)*10
   Y_MIN_M = INT(YIMIN/10.       )*10
   Y_MAX_M = INT(YIMAX/10.+0.9999)*10
   X_MIN   = X_MAX_M
   X_MAX   = X_MIN_M
   Y_MIN   = Y_MAX_M
   Y_MAX   = Y_MIN_M
  NL_MAX   = 0

  DO
    TMP_CHAR = ''
    READ ( 71, '(A)', IOSTAT = IOS ) TMP_CHAR

    IF ( IOS < 0 ) EXIT

    READ ( TMP_CHAR, * ) XT, YT, ZT, RTMP, RTMP, ZB
    NL = INT ( ZT - ZB )

    IF ( ( X_MIN_M <= XT ) .AND. ( XT <= X_MAX_M ) .AND.               &
         ( Y_MIN_M <= YT ) .AND. ( YT <= Y_MAX_M ) ) THEN

      IF ( XT <  X_MIN )  X_MIN = XT
      IF ( XT >  X_MAX )  X_MAX = XT
      IF ( YT <  Y_MIN )  Y_MIN = YT
      IF ( YT >  Y_MAX )  Y_MAX = YT
      IF ( NL > NL_MAX ) NL_MAX = NL

    END IF
  END DO

  REWIND ( 71 )

  I_MAX = (X_MAX - X_MIN) / 10 + 2
  K_MAX = (Y_MAX - Y_MIN) / 10 + 2
  L_MAX = NL_MAX + 2

!  print *, 'myid:', myid, I_MAX, K_MAX, L_MAX, X_MIN, Y_MIN
  ALLOCATE ( RHO (I_MAX, K_MAX, L_MAX), VP  (I_MAX, K_MAX, L_MAX),     &
             QS  (I_MAX, K_MAX, L_MAX), VS  (I_MAX, K_MAX, L_MAX) )

  RHO = -1.
  VP  = -1.
  VS  = -1.
  QS  = -1.


  READ ( 71, '(A)' ) TMP_CHAR
  READ ( 72, '(A)' ) TMP_CHAR
  READ ( 73, '(A)' ) TMP_CHAR
  READ ( 74, '(A)' ) TMP_CHAR

  DO
    TMP_CHAR = ''
    READ ( 71, '(A)', IOSTAT = IOS ) TMP_CHAR

    IF ( IOS < 0 ) EXIT

    READ ( TMP_CHAR, * ) XT, YT, ZT, RTMP, RTMP, ZB
    NL = INT ( ZT - ZB )

    IF ( ( X_MIN_M <= XT ) .AND. ( XT <= X_MAX_M ) .AND.               &
         ( Y_MIN_M <= YT ) .AND. ( YT <= Y_MAX_M ) ) THEN

      I = ( XT - X_MIN ) / 10 + 1
      K = ( YT - Y_MIN ) / 10 + 1

      READ ( TMP_CHAR, * ) XT, YT, ZT, RTMP, RTMP, ZB,                 &
                           ( RHO (I,K,L), L = 1, NL + 1 )

      TMP_CHAR = ''
      READ ( 72, '(A)' ) TMP_CHAR
      READ ( TMP_CHAR, * ) XT, YT, ZT, RTMP, RTMP, ZB,                 &
                           ( VP  (I,K,L), L = 1, NL + 1 )

      TMP_CHAR = ''
      READ ( 73, '(A)' ) TMP_CHAR
      READ ( TMP_CHAR, * ) XT, YT, ZT, RTMP, RTMP, ZB,                 &
                           ( VS  (I,K,L), L = 1, NL + 1 )

      TMP_CHAR = ''
      READ ( 74, '(A)' ) TMP_CHAR
      READ ( TMP_CHAR, * ) XT, YT, ZT, RTMP, RTMP, ZB,                 &
                           ( QS  (I,K,L), L = 1, NL + 1 )

    ELSE
      READ ( 72, '(A)' ) TMP_CHAR
      READ ( 73, '(A)' ) TMP_CHAR
      READ ( 74, '(A)' ) TMP_CHAR
    END IF

  END DO

  CLOSE ( 71 )
  CLOSE ( 72 )
  CLOSE ( 73 )
  CLOSE ( 74 )

END SUBROUTINE READ_GBM3D

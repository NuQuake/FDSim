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

SUBROUTINE  STO_BOR ( MXL, MXH, MYL, MYH, PXH, PYH, I1, K1,            &
                      X_MIN, X_MAX, Y_MIN, Y_MAX, T,                   &
                      SX_T, SY_T, SX_TP, SY_TP, TPML )

  USE PRECISION   , ONLY: PP
  USE NONREF_BOUND, ONLY: NR_DATA_X, NR_DATA_Y

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER                                      , INTENT (IN)    ::     &
                                     MXL, MXH, MYL, MYH, I1, K1, TPML, &
                                     PXH, PYH,                         &
                                     X_MIN, X_MAX, Y_MIN, Y_MAX
  TYPE (NR_DATA_X), DIMENSION (MXL:MXH)        , INTENT (INOUT) ::     &
                                                            SX_T, SX_TP
  TYPE (NR_DATA_Y), DIMENSION (MYL:MYH)        , INTENT (INOUT) ::     &
                                                            SY_T, SY_TP
  REAL (KIND=PP)  , DIMENSION (MXL:MXH,MYL:MYH), INTENT (INOUT) :: T


  INTEGER :: PXH1, PXH2, PYH1, PYH2, IT1, IT2, IT3, KT1, KT2, KT3

!-----------------------------------------------------------------------

  PXH1 = PXH - 1
  PXH2 = PXH - 2
  PYH1 = PYH - 1
  PYH2 = PYH - 2
  IT1  = I1  - TPML
  KT1  = K1  - TPML
  IT2  = IT1 + 1
  IT3  = IT1 + 2
  KT2  = KT1 + 1
  KT3  = KT1 + 2

  SX_TP = SX_T
  SY_TP = SY_T

  IF ( ( Y_MIN <= KT1 ) .AND. ( KT1 <= Y_MAX ) ) THEN
    SX_T (:)%TI1L   = T ( :   , KT1 )
    SX_T (:)%TI2L   = T ( :   , KT2 )
    SX_T (:)%TI3L   = T ( :   , KT3 )
  END IF

  IF ( Y_MAX == PYH ) THEN
    SX_T (:)%TIMY0L = T ( :   , PYH )
    SX_T (:)%TIMY1L = T ( :   , PYH1)
    SX_T (:)%TIMY2L = T ( :   , PYH2)
  END IF

  IF ( ( X_MIN <= IT1 ) .AND. ( IT1 <= X_MAX ) ) THEN
    SY_T (:)%T1KL   = T ( IT1 , :   )
    SY_T (:)%T2KL   = T ( IT2 , :   )
    SY_T (:)%T3KL   = T ( IT3 , :   )
  END IF

  IF ( X_MAX == PXH ) THEN
    SY_T (:)%TMX0KL = T ( PXH , :   )
    SY_T (:)%TMX1KL = T ( PXH1, :   )
    SY_T (:)%TMX2KL = T ( PXH2, :   )
  END IF

END SUBROUTINE STO_BOR

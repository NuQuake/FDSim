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
MODULE MOD_MODEL

  USE NRTYPE, ONLY: WP, SP

! module contains procedures, function and variables for
! creating and working with aproximation of sediments/bedrock interface

!------------------------------------------------------------- constants
  IMPLICIT NONE

  INTEGER, PARAMETER :: LAYER_1 = 1
  INTEGER, PARAMETER :: LAYER_2 = 2

  INTEGER, PARAMETER :: POSITION_TOP    = 0
  INTEGER, PARAMETER :: POSITION_BOTTOM = 1

!------------------------------------------------------ public variables

  INTEGER      :: NPX,NPY,NX10,NY10,MX,MXF,MY,MYF,MZ,LPAS,GRD_STEP,ni, &
                  NK, MX_G, MXF_G, NX10_G,                             &
                  X_MIN, Y_MIN
  REAL(SP)     :: H, H10, HF
  REAL(WP)     :: H2, HF2
  REAL(SP)     :: YBMIN, YBMAX, ZBMIN, ZHF,                            &
                  XBMIN_G, XBMAX_G
  REAL(WP)     :: XBMIN, XBMAX, XIMIN_G, XIMAX_G
  REAL(WP)     :: XIMIN, XIMAX, YIMIN, YIMAX, MAXDEPTH
  REAL(SP)     :: FLTR_RANGE
  LOGICAL      :: KEY_Q, KEY_FILTER, KEY_SRFOUT, KEY_QK_INF, KEY_GBM3D

  REAL(SP), ALLOCATABLE, DIMENSION (:,:,:) :: ITF1,                    &
                                              RHO, VP, VS, QS

END MODULE MOD_MODEL

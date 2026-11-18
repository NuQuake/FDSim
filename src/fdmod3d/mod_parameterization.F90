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
MODULE MOD_PARAMETERIZATION

USE NRTYPE,    ONLY: WP, SP
USE MOD_CONST, ONLY: N_FREQ, N_E_PAR, N_Y_PAR

!------------------------------------------------------------- constants
IMPLICIT NONE

INTEGER, PARAMETER :: FINER_GRID   = 0
INTEGER, PARAMETER :: COARSER_GRID = 1

INTEGER, PARAMETER :: STATUS_NEW     = 0
INTEGER, PARAMETER :: STATUS_OLD     = 1
INTEGER, PARAMETER :: STATUS_UNKNOWN = 2

!------------------------------------------------------------- variables
REAL(SP) :: RHO_HOM, MU_HOM, LA_HOM, QS_HOM, QP_HOM, VP_HOM, VS_HOM,   &
            FRJMAX, FRANGE, FREF, QP_OLD, QS_OLD
REAL(WP), DIMENSION(N_FREQ)   :: YMU_HOM, YL2M_HOM
REAL(WP), DIMENSION(N_FREQ,2) :: YMU_YL2M_OLD

INTEGER  :: N_PLANESF, N_PLANES, LSTART, LSTOP,                        &
            MAXCELLS, N_CELLS, LAST_CELL, REFINE, N_CELLS_G

INTEGER , DIMENSION(:,:,:), ALLOCATABLE :: JMH, JMH3, JMH_G, JMH3_G

REAL(SP), DIMENSION(N_E_PAR       ) :: CELL
REAL(SP), DIMENSION(N_Y_PAR,N_FREQ) :: CELLY

REAL(SP), DIMENSION(:,:  ), ALLOCATABLE :: CELLS,  CELLS_G
REAL(SP), DIMENSION(:,:,:), ALLOCATABLE :: CELLSY, CELLSY_G

END MODULE MOD_PARAMETERIZATION

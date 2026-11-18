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
!-----------------------------------------------------------------------
! print info about model on screen
!-----------------------------------------------------------------------
SUBROUTINE PRINT_INFO

  USE MOD_MODEL, ONLY: H, HF, H10, KEY_FILTER, ZHF, MX_G,MY, MXF_G,MYF,&
                       XIMIN, YIMIN, XIMAX, YIMAX, LPAS, MZ,           &
                       XBMIN_G, YBMIN, XBMAX_G, YBMAX, ZBMIN, XIMAX_G

  IMPLICIT NONE

  WRITE(*,*)'Coordinate system: x -> east'
  WRITE(*,*)'                   y -> north'
  WRITE(*,*)'                   z -> up'
  WRITE(*,*)'Spacing in coarser grid: h   = ',H
  WRITE(*,*)'Spacing in finer grid:   hf  = ',HF
  WRITE(*,*)'Spacing for integration: h10 = ',H10
  IF (KEY_FILTER) THEN
    WRITE(*,*)'Smooth interface sediments/bedrock is ON'
  ELSE
    WRITE(*,*)'Smooth interface sediments/bedrock is OFF'
  END IF
  WRITE(*,*)'Bounderies for interface sediments/bedrock aproximation:'
  WRITE(*,*)'  ximin = ',XIMIN,' ximax = ',XIMAX_G
  WRITE(*,*)'  yimin = ',YIMIN,' yimax = ',YIMAX
  WRITE(*,*)'Calculation model boundaries: '
  WRITE(*,*)'  xbmin = ',XBMIN_G,' xbmax = ',XBMAX_G
  WRITE(*,*)'  ybmin = ',YBMIN,' ybmax = ',YBMAX
  WRITE(*,*)'  zbmin = ',ZBMIN,' zhf  = ',ZHF
  WRITE(*,*)'Number of grid cells in '
  WRITE(*,*)'  finer   grid: ',MXF_G, MYF, LPAS+1
  WRITE(*,*)'  coarser grid: ',MX_G , MY , MZ-LPAS+1
  WRITE(*,*)'Model parameters: MX   = ',MY
  WRITE(*,*)'                  MY   = ',MX_G
  WRITE(*,*)'                  MZ   = ',MZ
  WRITE(*,*)'                  LPAS = ',LPAS
  WRITE(*,*)' '

END SUBROUTINE PRINT_INFO

!-----------------------------------------------------------------------
! logs info about model into a file
!-----------------------------------------------------------------------
SUBROUTINE LOG_INFO

  USE MOD_PARAMETERIZATION, ONLY: N_CELLS_G
  USE MOD_MODEL, ONLY: H, HF, H10, KEY_FILTER, ZHF, MX_G,MY, MXF_G,MYF,&
                       XIMIN, YIMIN, XIMAX, YIMAX, LPAS, MZ,           &
                       XBMIN_G, YBMIN, XBMAX_G, YBMAX, ZBMIN, XIMAX_G

  IMPLICIT NONE

  OPEN(99,file='MODEL.LOG')
  WRITE(99,*)'Coordinate system: x -> east'
  WRITE(99,*)'                   y -> north'
  WRITE(99,*)'                   z -> up'
  WRITE(99,*)'Spacing in coarser grid: h   = ',H
  WRITE(99,*)'Spacing in finer grid:   hf  = ',HF
  WRITE(99,*)'Spacing for integration: h10 = ',H10
  IF (KEY_FILTER) THEN
    WRITE(99,*)'Smooth interface sediments/bedrock is ON'
  ELSE
    WRITE(99,*)'Smooth interface sediments/bedrock is OFF'
  END IF
  WRITE(99,*)'Bounderies for interface sediments/bedrock aproximation:'
  WRITE(99,*)'  ximin = ',XIMIN,' ximax = ',XIMAX_G
  WRITE(99,*)'  yimin = ',YIMIN,' yimax = ',YIMAX
  WRITE(99,*)'Calculation model boundaries: '
  WRITE(99,*)'  xbmin = ',XBMIN_G,' xbmax = ',XBMAX_G
  WRITE(99,*)'  ybmin = ',YBMIN,' ybmax = ',YBMAX
  WRITE(99,*)'  zbmin = ',ZBMIN,' zhf  = ',ZHF
  WRITE(99,*)'Number of grid cells in '
  WRITE(99,*)'  finer   grid: ',MXF_G, MYF, LPAS+1
  WRITE(99,*)'  coarser grid: ',MX_G , MY , MZ-LPAS+1
  WRITE(99,*)'Model parameters: MX   = ',MY
  WRITE(99,*)'                  MY   = ',MX_G
  WRITE(99,*)'                  MZ   = ',MZ
  WRITE(99,*)'                  LPAS = ',LPAS
  WRITE(99,*)' '
  WRITE(99,*)'TYPES OF MAT. CELLS: ',N_CELLS_G
  WRITE(99,*)' '

END SUBROUTINE LOG_INFO

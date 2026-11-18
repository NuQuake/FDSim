!=======================================================================
! Copyright (c) 2014-2026 Jozef Kristek, Peter Moczo
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
! 2026-11-18: independent implementation (Hermite form), same interface
!=======================================================================
!
! Coefficients of the bicubic interpolant on one rectangular cell.
!
! Corners are numbered counter-clockwise from the lower left:
!   1 = (0,0), 2 = (1,0), 3 = (1,1), 4 = (0,1)   in local coordinates
!   t = (x1 - x1_corner1)/D1,  u = (x2 - x2_corner1)/D2.
!
! Input at each corner: value Y, derivatives Y1 = df/dx1, Y2 = df/dx2,
! cross derivative Y12 = d2f/dx1dx2; D1, D2 are the cell sizes.
!
! Output C(i,j), i,j = 1..4, such that
!   f(t,u) = SUM_{i,j} C(i,j) * t**(i-1) * u**(j-1).
!
! With the cubic Hermite matrix
!   M = | 1  0  0  0 |
!       | 0  0  1  0 |
!       |-3  3 -2 -1 |
!       | 2 -2  1  1 |
! and the corner-data matrix
!   G = | f(0,0)   f(0,1)   fu(0,0)   fu(0,1)  |
!       | f(1,0)   f(1,1)   fu(1,0)   fu(1,1)  |
!       | ft(0,0)  ft(0,1)  ftu(0,0)  ftu(0,1) |
!       | ft(1,0)  ft(1,1)  ftu(1,0)  ftu(1,1) |
! (derivatives with respect to the local coordinates t, u)
! the coefficients are  C = M G M^T.
!-----------------------------------------------------------------------
SUBROUTINE BCUCOF(Y,Y1,Y2,Y12,D1,D2,C)

  USE PRECISION, ONLY: PP

  IMPLICIT NONE

  REAL(PP),                 INTENT(IN)  :: D1,D2
  REAL(PP), DIMENSION(4),   INTENT(IN)  :: Y,Y1,Y2,Y12
  REAL(PP), DIMENSION(4,4), INTENT(OUT) :: C

  REAL(PP), DIMENSION(4,4), PARAMETER :: M = RESHAPE (                 &
       (/  1._PP, 0._PP,-3._PP, 2._PP,                                 &
           0._PP, 0._PP, 3._PP,-2._PP,                                 &
           0._PP, 1._PP,-2._PP, 1._PP,                                 &
           0._PP, 0._PP,-1._PP, 1._PP /), (/4,4/) )

  REAL(PP), DIMENSION(4,4) :: G

!-----------------------------------------------------------------------
! corner data in local coordinates; corner (t,u): 1=(0,0) 2=(1,0) 3=(1,1) 4=(0,1)

  G(1,1) = Y(1)
  G(1,2) = Y(4)
  G(2,1) = Y(2)
  G(2,2) = Y(3)

  G(1,3) = Y2(1)*D2
  G(1,4) = Y2(4)*D2
  G(2,3) = Y2(2)*D2
  G(2,4) = Y2(3)*D2

  G(3,1) = Y1(1)*D1
  G(3,2) = Y1(4)*D1
  G(4,1) = Y1(2)*D1
  G(4,2) = Y1(3)*D1

  G(3,3) = Y12(1)*D1*D2
  G(3,4) = Y12(4)*D1*D2
  G(4,3) = Y12(2)*D1*D2
  G(4,4) = Y12(3)*D1*D2

  C = MATMUL ( MATMUL ( M, G ), TRANSPOSE (M) )

END SUBROUTINE BCUCOF

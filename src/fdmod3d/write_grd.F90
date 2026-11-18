!=======================================================================
! Copyright (c) 2014 Jozef Kristek, Peter Moczo, Martin Galis
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
! Version 1.0, 2014-03-19: Jozef Kristek, Peter Moczo, Martin Galis
!=======================================================================
!-----------------------------------------------------------------------
! write aproximation of sediments/bedrock interface to disk in Surfer 
! grid format
!   GRD_STEP = 1 .. output mesh is with spatial step H10
!   GRD_STEP = 2 .. output mesh is with spatial step HF
!   GRD_STEP = 3 .. output mesh is with spatial step H
!-----------------------------------------------------------------------
SUBROUTINE WRITE_GRD

  USE NRTYPE   , ONLY: WP
  USE MOD_MODEL, ONLY: NX10, NY10, H10, XIMAX, XIMIN, YIMAX, YIMIN, H, &
                       ITF1, GRD_STEP, HF, ni
  USE MOD_FILES, ONLY: SRFNAME

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER      :: I,J,K,NX,NY,DI
  REAL(WP)     :: DX
  character(30):: txtbuf
!----------------------------------------------------------------------

  WRITE(*,'(A)')'Saving grd files...'

  SELECT CASE (GRD_STEP)
    CASE(1)
      !--------- output with step h10
      NX = NX10
      NY = NY10
      DX = H10
      DI = 1
    CASE(2)
      !--------- output with step hf
      NX = (XIMAX-XIMIN)/HF + 1
      NY = (YIMAX-YIMIN)/HF + 1
      DX = HF
      DI = HF/H10
    CASE(3)
      !--------- output with step h
      NX = (XIMAX-XIMIN)/H + 1
      NY = (YIMAX-YIMIN)/H + 1
      DX = H
      DI = H/H10
  END SELECT

  do k = 1, ni

    write(txtbuf,'(i2.2,a)')k,srfname
    OPEN ( 30, FILE=txtbuf )

    WRITE( 30,'(A)' )'DSAA'
    WRITE( 30,* )NY,NX
    WRITE( 30,* )YIMIN, YIMIN+DX*DBLE(NY-1)
    WRITE( 30,* )XIMIN, XIMIN+DX*DBLE(NX-1)
    WRITE( 30,* )MINVAL(ITF1(:,:,k)),MAXVAL(ITF1(:,:,k))

    DO I = 1, NX10, DI
      WRITE(30,'(10(E15.6,1X))' ) (ITF1(I,J,k), J=1,NY10,DI)
    END DO

    CLOSE(30)
  
  end do

  WRITE(*,'(A)')'OK'

END SUBROUTINE WRITE_GRD

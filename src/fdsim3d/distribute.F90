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

SUBROUTINE  DISTRIBUTE

  USE GRID_MEDIUM , ONLY:                                              &
                          H

  USE CONTROL_DATA, ONLY:                                              &
                          KEY_TLV, KEY_SNV, REC_NAME,                  &
                          MT1, MT2, IPAS2, MSNP,                       &
                          FRJMAX, FRANGE, N_FREQ,                      &
                          DT, LPAS, MR, MZ, TPML,                      &
                          IREC, KREC, LREC,                            &
                          TAU, R, POW, SF, WC,                         &
                          KTLE, KTRI, KTRE, KTFR, KTBO,                &
                          OMG, WB,                                     &
                    THPPLE, THPPRI, THPPRE, THPPFR, THPPBO,            &
                    THSSLE, THSSRI, THSSRE, THSSFR, THSSBO

  USE AUXIL       , ONLY:                                              &
                          LOGSNP
#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED

!-----------------------------------------------------------------------

  IMPLICIT NONE
  CHARACTER(LEN=6) :: REC_NAME_T
  INTEGER          :: J

!-----------------------------------------------------------------------

  CALL MPI_Bcast(H,         1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(KEY_SNV  , 1, MPI_LOGICAL, root, comm, errcode)
  CALL MPI_Bcast(LPAS,      1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(MT1,       1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(MT2,       1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(IPAS2,     1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(MSNP,      1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(DT,        1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(FRJMAX,    1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(FRANGE,    1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(N_FREQ,    1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(TAU,       1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(R,         1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(POW,       1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(SF,        1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(WC,        1, MPI_REAL   , root, comm, errcode)

  CALL MPI_Bcast(KTLE,      1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(KTRI,      1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(KTRE,      1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(KTFR,      1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(KTBO,      1, MPI_INTEGER, root, comm, errcode)

  CALL MPI_Bcast(OMG,       1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(WB ,       1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(THPPLE,    1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(THPPRI,    1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(THPPRE,    1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(THPPFR,    1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(THPPBO,    1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(THSSLE,    1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(THSSRI,    1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(THSSRE,    1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(THSSFR,    1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(THSSBO,    1, MPI_REAL   , root, comm, errcode)

  IF (KEY_TLV) THEN
    DO J = 1, MR
      REC_NAME_T = REC_NAME (J)
      CALL MPI_Bcast(REC_NAME_T, 6, MPI_CHARACTER, root, comm, errcode)
      REC_NAME(J) = REC_NAME_T
    END DO
    CALL MPI_Bcast(IREC,   MR, MPI_INTEGER, root, comm, errcode)
    CALL MPI_Bcast(KREC,   MR, MPI_INTEGER, root, comm, errcode)
    CALL MPI_Bcast(LREC,   MR, MPI_INTEGER, root, comm, errcode)
  END IF

  CALL MPI_Bcast(LOGSNP, MZ+1+TPML,                                    &
                               MPI_LOGICAL, root, comm, errcode)

#endif

END SUBROUTINE  DISTRIBUTE

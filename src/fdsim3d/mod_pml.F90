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

MODULE PML
  USE PRECISION, ONLY: PP
  IMPLICIT NONE

  REAL (KIND=PP), DIMENSION (:,:,:), ALLOCATABLE ::                    &
            ODA_XL ,  ODB_XL ,  ODC_XL ,  ODA_XH ,  ODB_XH ,  ODC_XH , &
            ODD_XL ,                      ODD_XH ,                     &
            ODA_YL ,  ODB_YL ,  ODC_YL ,  ODA_YH ,  ODB_YH ,  ODC_YH , &
            ODD_YL ,                      ODD_YH ,                     &
                                          ODA_ZH ,  ODB_ZH ,  ODC_ZH , &
                                          ODD_ZH ,                     &
           ODA2_XL , ODB2_XL , ODC2_XL , ODA2_XH , ODB2_XH , ODC2_XH , &
           ODD2_XL ,                     ODD2_XH ,                     &
           ODA2_YL , ODB2_YL , ODC2_YL , ODA2_YH , ODB2_YH , ODC2_YH , &
           ODD2_YL ,                     ODD2_YH ,                     &
                                         ODA2_ZH , ODB2_ZH , ODC2_ZH , &
                                         ODD2_ZH ,                     &
            ODA_XLF,  ODB_XLF,  ODC_XLF,  ODA_XHF,  ODB_XHF,  ODC_XHF, &
            ODD_XLF,                      ODD_XHF,                     &
            ODA_YLF,  ODB_YLF,  ODC_YLF,  ODA_YHF,  ODB_YHF,  ODC_YHF, &
            ODD_YLF,                      ODD_YHF,                     &
           ODA2_XLF, ODB2_XLF, ODC2_XLF, ODA2_XHF, ODB2_XHF, ODC2_XHF, &
           ODD2_XLF,                     ODD2_XHF,                     &
           ODA2_YLF, ODB2_YLF, ODC2_YLF, ODA2_YHF, ODB2_YHF, ODC2_YHF, &
           ODD2_YLF,                     ODD2_YHF

  REAL (KIND=PP), DIMENSION (:,:,:), ALLOCATABLE ::                    &
                                        PZXL , PZXH ,                  &
                                        PZYL , PZYH ,                  &
                                        PXXL , PXXH ,                  &
                                        PYYL , PYYH ,                  &
                                        PYXL , PYXH ,                  &
                                        PXYL , PXYH ,                  &
                                               PXZH ,                  &
                                               PYZH ,                  &
                                               PZZH ,                  &
                                          RXXL , RXXH ,                &
                                          RXYL , RXYH ,                &
                                          RYYL , RYYH ,                &
                                          RYXL , RYXH ,                &
                                          RZXL , RZXH ,                &
                                          RZYL , RZYH ,                &
                                  RXZH ,                               &
                                  RYZH ,                               &
                                  RZZH ,                               &
                                        PZXLF, PZXHF,                  &
                                        PZYLF, PZYHF,                  &
                                        PXXLF, PXXHF,                  &
                                        PYYLF, PYYHF,                  &
                                        PYXLF, PYXHF,                  &
                                        PXYLF, PXYHF,                  &
                                          RXXLF, RXXHF,                &
                                          RXYLF, RXYHF,                &
                                          RYYLF, RYYHF,                &
                                          RYXLF, RYXHF,                &
                                          RZXLF, RZXHF,                &
                                          RZYLF, RZYHF


END MODULE PML

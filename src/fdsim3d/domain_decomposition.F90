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

SUBROUTINE  DOMAIN_DECOMPOSITION

  USE PRECISION   , ONLY:                                              &
                           PP
  USE CONTROL_DATA, ONLY:                                              &
                           LPAS, REFINE
  USE AUXIL       , ONLY:                                              &
                           PMLL , PMLXH  , PMLYH  , PMLZH ,            &
                           X_MIN_D , X_MAX_D , X_MIN , X_MAX ,         &
                           Y_MIN_D , Y_MAX_D , Y_MIN , Y_MAX ,         &
                           Z_MIN_D , Z_MAX_D , Z_MIN , Z_MAX ,         &
                           PMLLF, PMLXHF , PMLYHF ,                    &
                           X_MIN_DF, X_MAX_DF, X_MINF, X_MAXF,         &
                           Y_MIN_DF, Y_MAX_DF, Y_MINF, Y_MAXF,         &
                           Z_MIN_DF, Z_MAX_DF, Z_MINF, Z_MAXF,         &
                           TMP, BUFSF, BUFS,                           &
                        BUFXS , BUFYS , BUFZS , BUFXR , BUFYR , BUFZR ,&
                        BUFXS_P,BUFYS_P,        BUFXR_P,BUFYR_P,       &
                        BUFXSF, BUFYSF, BUFZSF, BUFXRF, BUFYRF, BUFZRF,&
                           source_id0A , dest_id0A ,                   &
                           source_id0B , dest_id0B ,                   &
                           source_id1A , dest_id1A ,                   &
                           source_id1B , dest_id1B ,                   &
                           source_id2A , dest_id2A ,                   &
                           source_id2B , dest_id2B ,                   &
                           source_id2AF, dest_id2AF,                   &
                           source_id2BF, dest_id2BF

#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif

!-----------------------------------------------------------------------

  IMPLICIT NONE

#if(USE_MPI)
  LOGICAL :: REORDER

  INTEGER :: PROCX, PROCY, PROCZ, J, XS, YS, root1
  REAL(PP):: WL
#endif

!-----------------------------------------------------------------------

#if(USE_MPI)

  OVERLAP = 2

  periods = .FALSE.
  REORDER = .TRUE.

  IF ( myid == 0 ) THEN
    iamroot = .TRUE.
  ELSE
    iamroot = .FALSE.
  END IF

  CALL MPI_DIMS_CREATE (numprocs,ndims,dims, errcode)
  PROCX = dims(1)
  PROCY = dims(2)
  PROCZ = dims(3)

  print *,ndims,dims

  CALL MPI_CART_CREATE                                                 &
              (comm, ndims, dims, periods, REORDER, comm_cart, errcode)
  CALL MPI_CART_GET                                                    &
         (comm_cart, ndims, dims, periods, my_coord, errcode)
  CALL MPI_COMM_RANK (comm_cart, myid, errcode)

  IF ( iamroot) root1 = myid
  CALL MPI_Bcast(root1,      1, MPI_INTEGER, root, comm, errcode)
  root = root1

  ALLOCATE ( X_MIN_A (0:PROCX-1), X_MAX_A (0:PROCX-1),                 &
             Y_MIN_A (0:PROCY-1), Y_MAX_A (0:PROCY-1),                 &
             Z_MIN_A (0:PROCZ-1), Z_MAX_A (0:PROCZ-1) )

  X_MIN_A(0        ) = PMLL
  X_MAX_A(dims(1)-1) = PMLXH
  DO J = 0, dims(1)-2
    X_MAX_A(J  ) = (J+1)*(PMLXH-PMLL+1)/PROCX + PMLL
    X_MIN_A(J+1) = (J+1)*(PMLXH-PMLL+1)/PROCX + PMLL + 1
  END DO

  Y_MIN_A(0        ) = PMLL
  Y_MAX_A(dims(2)-1) = PMLYH
  DO J = 0, dims(2)-2
    Y_MAX_A(J  ) = (J+1)*(PMLYH-PMLL+1)/PROCY + PMLL
    Y_MIN_A(J+1) = (J+1)*(PMLYH-PMLL+1)/PROCY + PMLL + 1
  END DO

  IF ( LPAS /= 0 ) THEN
    ALLOCATE ( X_MIN_AF (0:PROCX-1), X_MAX_AF (0:PROCX-1),             &
               Y_MIN_AF (0:PROCY-1), Y_MAX_AF (0:PROCY-1),             &
               Z_MIN_AF (0:PROCZ-1), Z_MAX_AF (0:PROCZ-1) )

    X_MIN_AF(0        ) = PMLLF
    X_MAX_AF(dims(1)-1) = PMLXHF
    DO J = 0, dims(1)-2
      X_MAX_AF(J  ) = (X_MAX_A (J)-1)*REFINE+1
      X_MIN_AF(J+1) =  X_MAX_AF(J) + 1
    END DO

    Y_MIN_AF(0        ) = PMLLF
    Y_MAX_AF(dims(2)-1) = PMLYHF
    DO J = 0, dims(2)-2
      Y_MAX_AF(J  ) = (Y_MAX_A (J)-1)*REFINE+1
      Y_MIN_AF(J+1) =  Y_MAX_AF(J) + 1
    END DO

    Z_MIN_AF(0        ) = 0
    Z_MAX_AF(dims(3)-1) = LPAS
    Z_MIN_A (0        ) = LPAS-2
    Z_MAX_A (dims(3)-1) = PMLZH
    WL = ((PMLZH-LPAS+1)/REFINE+(LPAS+1)*REFINE)/PROCZ
    DO J = 0, dims(3)-2
      Z_MAX_AF(J  ) = MIN(NINT((J+1)*WL/REFINE),LPAS)
      Z_MIN_AF(J+1) = MIN(NINT((J+1)*WL/REFINE),LPAS) + 1
    END DO
    DO J = dims(3)-1,1,-1
      Z_MIN_A (J  ) = MAX(PMLZH-NINT((PROCZ-J)*WL*REFINE),LPAS-2)
      Z_MAX_A (J-1) = MAX(PMLZH-NINT((PROCZ-J)*WL*REFINE),LPAS-2)-1
    END DO

    X_MINF = X_MIN_AF ( my_coord(1) )
    X_MAXF = X_MAX_AF ( my_coord(1) )
    Y_MINF = Y_MIN_AF ( my_coord(2) )
    Y_MAXF = Y_MAX_AF ( my_coord(2) )
    Z_MINF = Z_MIN_AF ( my_coord(3) )
    Z_MAXF = Z_MAX_AF ( my_coord(3) )

  ELSE
    Z_MIN_A(0        ) = 0
    Z_MAX_A(dims(3)-1) = PMLZH
    DO J = 0, dims(3)-2
      Z_MAX_A(J  ) = (J+1)*(PMLZH     +1)/PROCZ
      Z_MIN_A(J+1) = (J+1)*(PMLZH     +1)/PROCZ + 1
    END DO
  END IF

  X_MIN = X_MIN_A ( my_coord(1) )
  X_MAX = X_MAX_A ( my_coord(1) )
  Y_MIN = Y_MIN_A ( my_coord(2) )
  Y_MAX = Y_MAX_A ( my_coord(2) )
  Z_MIN = Z_MIN_A ( my_coord(3) )
  Z_MAX = Z_MAX_A ( my_coord(3) )


  ALLOCATE ( X_MIN_DA (0:PROCX-1), X_MAX_DA (0:PROCX-1),               &
             Y_MIN_DA (0:PROCY-1), Y_MAX_DA (0:PROCY-1),               &
             Z_MIN_DA (0:PROCZ-1), Z_MAX_DA (0:PROCZ-1) )

  X_MIN_DA(0        ) = PMLL
  X_MAX_DA(dims(1)-1) = PMLXH
  DO J = 0, dims(1)-2
    X_MAX_DA(J  ) = X_MAX_A(J  ) + OVERLAP
    X_MIN_DA(J+1) = X_MIN_A(J+1) - OVERLAP
  END DO

  Y_MIN_DA(0        ) = PMLL
  Y_MAX_DA(dims(2)-1) = PMLYH
  DO J = 0, dims(2)-2
    Y_MAX_DA(J  ) = Y_MAX_A(J  ) + OVERLAP
    Y_MIN_DA(J+1) = Y_MIN_A(J+1) - OVERLAP
  END DO

  IF ( LPAS /= 0 ) THEN
    ALLOCATE ( X_MIN_DAF (0:PROCX-1), X_MAX_DAF (0:PROCX-1),           &
               Y_MIN_DAF (0:PROCY-1), Y_MAX_DAF (0:PROCY-1),           &
               Z_MIN_DAF (0:PROCZ-1), Z_MAX_DAF (0:PROCZ-1) )

    X_MIN_DAF(0        ) = PMLLF
    X_MAX_DAF(dims(1)-1) = PMLXHF
    DO J = 0, dims(1)-2
      X_MAX_DAF(J  ) = X_MAX_AF(J  ) + OVERLAP
      X_MIN_DAF(J+1) = X_MIN_AF(J+1) - OVERLAP
    END DO

    Y_MIN_DAF(0        ) = PMLLF
    Y_MAX_DAF(dims(2)-1) = PMLYHF
    DO J = 0, dims(2)-2
      Y_MAX_DAF(J  ) = Y_MAX_AF(J  ) + OVERLAP
      Y_MIN_DAF(J+1) = Y_MIN_AF(J+1) - OVERLAP
    END DO

    Z_MIN_DAF(0        ) = 0
    Z_MAX_DAF(dims(3)-1) = LPAS
    Z_MIN_DA (0        ) = LPAS-2
    Z_MAX_DA (dims(3)-1) = PMLZH
    WL = ((PMLZH-LPAS+1)/REFINE+(LPAS+1)*REFINE)/PROCZ
    DO J = 0, dims(3)-2
      Z_MAX_DAF(J  ) = MIN(NINT((J+1)*WL/REFINE) + OVERLAP, LPAS)
      Z_MIN_DAF(J+1) = MIN(NINT((J+1)*WL/REFINE) - OVERLAP, LPAS) + 1
    END DO
    DO J = dims(3)-1,1,-1
      Z_MIN_DA (J  ) = MAX(PMLZH-NINT((PROCZ-J)*WL*REFINE)+1-OVERLAP,LPAS-2)
      Z_MAX_DA (J-1) = MAX(PMLZH-NINT((PROCZ-J)*WL*REFINE)+1+OVERLAP,LPAS-2)-1
    END DO
  ELSE
    Z_MIN_DA(0        ) = 0
    Z_MAX_DA(dims(3)-1) = PMLZH
    DO J = 0, dims(3)-2
      Z_MAX_DA(J  ) = (J+1)*(PMLZH     +1)/PROCZ     + OVERLAP
      Z_MIN_DA(J+1) = (J+1)*(PMLZH     +1)/PROCZ + 1 - OVERLAP
    END DO
  END IF

  IF ( LPAS /= 0 ) THEN
    X_MIN_DF = X_MIN_DAF ( my_coord(1) )
    X_MAX_DF = X_MAX_DAF ( my_coord(1) )
    Y_MIN_DF = Y_MIN_DAF ( my_coord(2) )
    Y_MAX_DF = Y_MAX_DAF ( my_coord(2) )
    Z_MIN_DF = Z_MIN_DAF ( my_coord(3) )
    Z_MAX_DF = Z_MAX_DAF ( my_coord(3) )
  END IF

  X_MIN_D = X_MIN_DA ( my_coord(1) )
  X_MAX_D = X_MAX_DA ( my_coord(1) )
  Y_MIN_D = Y_MIN_DA ( my_coord(2) )
  Y_MAX_D = Y_MAX_DA ( my_coord(2) )
  Z_MIN_D = Z_MIN_DA ( my_coord(3) )
  Z_MAX_D = Z_MAX_DA ( my_coord(3) )



  ALLOCATE   ( BUFXS(1:OVERLAP,1:Y_MAX-Y_MIN+1,1:Z_MAX-Z_MIN+1),       &
               BUFYS(X_MIN_D:X_MAX_D,1:OVERLAP,1:Z_MAX-Z_MIN+1),       &
               BUFZS(X_MIN_D:X_MAX_D,Y_MIN_D:Y_MAX_D,1:OVERLAP) )
  ALLOCATE   ( BUFXR(1:OVERLAP,1:Y_MAX-Y_MIN+1,1:Z_MAX-Z_MIN+1),       &
               BUFYR(X_MIN_D:X_MAX_D,1:OVERLAP,1:Z_MAX-Z_MIN+1),       &
               BUFZR(X_MIN_D:X_MAX_D,Y_MIN_D:Y_MAX_D,1:OVERLAP) )

  ALLOCATE ( BUFXS_P (1:OVERLAP,1:Y_MAX-Y_MIN+1),                      &
             BUFYS_P (X_MIN_D:X_MAX_D,1:OVERLAP) )
  ALLOCATE ( BUFXR_P (1:OVERLAP,1:Y_MAX-Y_MIN+1),                      &
             BUFYR_P (X_MIN_D:X_MAX_D,1:OVERLAP) )

  CALL MPI_CART_SHIFT (comm_cart, 0, 1, source_id0A, dest_id0A, errcode)
  CALL MPI_CART_SHIFT (comm_cart, 0,-1, source_id0B, dest_id0B, errcode)
  CALL MPI_CART_SHIFT (comm_cart, 1, 1, source_id1A, dest_id1A, errcode)
  CALL MPI_CART_SHIFT (comm_cart, 1,-1, source_id1B, dest_id1B, errcode)
  CALL MPI_CART_SHIFT (comm_cart, 2, 1, source_id2A, dest_id2A, errcode)
  CALL MPI_CART_SHIFT (comm_cart, 2,-1, source_id2B, dest_id2B, errcode)

  IF ( LPAS /= 0 ) THEN

    IF      ( Z_MAX < Z_MIN ) THEN
      source_id2AF = source_id2A
        dest_id2AF =   dest_id2A
      source_id2BF = source_id2B
        dest_id2BF =   dest_id2B
      source_id2A  = MPI_PROC_NULL
        dest_id2A  = MPI_PROC_NULL
      source_id2B  = MPI_PROC_NULL
        dest_id2B  = MPI_PROC_NULL
    ELSE IF ( Z_MAXF < Z_MINF ) THEN
      source_id2AF = MPI_PROC_NULL
        dest_id2AF = MPI_PROC_NULL
      source_id2BF = MPI_PROC_NULL
        dest_id2BF = MPI_PROC_NULL
    ELSE
      source_id2AF = source_id2A
      source_id2A  = MPI_PROC_NULL
        dest_id2AF = MPI_PROC_NULL
      source_id2BF = MPI_PROC_NULL
        dest_id2BF =   dest_id2B
        dest_id2B  = MPI_PROC_NULL
    END IF

    ALLOCATE ( TMP (X_MIN:X_MAX,Y_MIN:Y_MAX) )
    ALLOCATE ( BUFXSF(1:OVERLAP,1:Y_MAXF-Y_MINF+1,1:Z_MAXF-Z_MINF+1),  &
               BUFYSF(X_MIN_DF:X_MAX_DF,1:OVERLAP,1:Z_MAXF-Z_MINF+1),  &
               BUFZSF(X_MIN_DF:X_MAX_DF,Y_MIN_DF:Y_MAX_DF,1:OVERLAP) )
    ALLOCATE ( BUFXRF(1:OVERLAP,1:Y_MAXF-Y_MINF+1,1:Z_MAXF-Z_MINF+1),  &
               BUFYRF(X_MIN_DF:X_MAX_DF,1:OVERLAP,1:Z_MAXF-Z_MINF+1),  &
               BUFZRF(X_MIN_DF:X_MAX_DF,Y_MIN_DF:Y_MAX_DF,1:OVERLAP) )

    XS = MAXVAL(X_MAX_AF-X_MIN_AF+1)
    YS = MAXVAL(Y_MAX_AF-Y_MIN_AF+1)
    ALLOCATE ( BUFSF (1:XS*YS) )
  END IF

  XS = MAXVAL(X_MAX_A-X_MIN_A+1)
  YS = MAXVAL(Y_MAX_A-Y_MIN_A+1)
  ALLOCATE ( BUFS (1:XS*YS) )

#else

  IF ( LPAS /= 0 ) THEN

    X_MIN_DF = PMLLF
    X_MAX_DF = PMLXHF
    Y_MIN_DF = PMLLF
    Y_MAX_DF = PMLYHF
    Z_MIN_DF = 0
    Z_MAX_DF = LPAS

    X_MINF = PMLLF
    X_MAXF = PMLXHF
    Y_MINF = PMLLF
    Y_MAXF = PMLYHF
    Z_MINF = 0
    Z_MAXF = LPAS

    X_MIN_D = PMLL
    X_MAX_D = PMLXH
    Y_MIN_D = PMLL
    Y_MAX_D = PMLYH
    Z_MIN_D = LPAS-2
    Z_MAX_D = PMLZH

    X_MIN = PMLL
    X_MAX = PMLXH
    Y_MIN = PMLL
    Y_MAX = PMLYH
    Z_MIN = LPAS-2
    Z_MAX = PMLZH

    ALLOCATE ( TMP (X_MIN:X_MAX,Y_MIN:Y_MAX) )

  ELSE
    X_MIN_D = PMLL
    X_MAX_D = PMLXH
    Y_MIN_D = PMLL
    Y_MAX_D = PMLYH
    Z_MIN_D = 0
    Z_MAX_D = PMLZH

    X_MIN = PMLL
    X_MAX = PMLXH
    Y_MIN = PMLL
    Y_MAX = PMLYH
    Z_MIN = 0
    Z_MAX = PMLZH
  END IF

#endif

END SUBROUTINE  DOMAIN_DECOMPOSITION

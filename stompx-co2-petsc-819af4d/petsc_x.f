!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE STOMP_PETSC_CREATE( KSPX,MATX,RHS_VECX,SOL_VECX,
     &  SOL_VEC_SX,SCATTERX,INDX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Create solution vector, problem vector, matrix, solver, and
!     solution-vector scatter scheme for the PETSc solver
!
!     INDX = 0 : Coupled flow and transport
!     INDX = 1 : Solute or species component transport
!     INDX = 2 : Geomechanics
!
!     Determine the number of nonzero entries in the solute/species
!     transport matrix
!
!     ID_NZ - number of nonzeros per row in DIAGONAL portion of 
!     local submatrix (same value is used for all local rows)
!
!     ID_NNZ - array containing the number of nonzeros in the various 
!     rows of the DIAGONAL portion of the local submatrix 
!     (possibly different for each row) or NULL, if d_nz is used to 
!     specify the nonzero structure. The size of this array is equal 
!     to the number of local rows, i.e ‘m’.
!
!     IO_NZ - number of nonzeros per row in the OFF-DIAGONAL portion 
!     of local submatrix (same value is used for all local rows).
!
!     IO_NNZ - array containing the number of nonzeros in the various 
!     rows of the OFF-DIAGONAL portion of the local submatrix 
!     (possibly different for each row) or NULL, if o_nz is used to 
!     specify the nonzero structure. The size of this array is equal 
!     to the number of local rows, i.e ‘m’.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 1 December 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE JACOB
      USE GRID
      USE GEO_MECH
      USE COUP_WELL
      USE PETSCKSP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

!
!
!  Include file for Fortran use of the Mat package in PETSc
!  Portions of this code are under:
!  Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
!




!
!
!  Include file for Fortran use of the Vec package in PETSc
!




!
!
!  Include file for Fortran use of the AO (application ordering) package in PETSc
!




!
!
!  Include file for Fortran use of the IS (index set) package in PETSc
!




!
!
!  Part of the base include file for Fortran use of PETSc.
!  Note: This file should contain only define statements and
!  not the declaration of variables.

! No spaces for #defines as some compilers (PGI) also adds
! those additional spaces during preprocessing - bad for fixed format
!
























































!
!  Include file for Fortran use of the PetscViewer package in PETSc
!













!
!  No includes needed for logging

!
!  Include file for Fortran use of the Bag package in PETSc
!







!
! The real*8,complex*16 notatiton is used so that the
! PETSc double/complex variables are not affected by
! compiler options like -r4,-r8, sometimes invoked
! by the user. NAG compiler does not like integer*4,real*8
























!
! Fortran does not support unsigned, though ISO_C_BINDING
! supports INTEGER(KIND=C_SIZE_T). We don't use that here
! only to avoid importing the module.





!

!






!

!


!



!
!     Macro for templating between real and complex
!










!
!    Allows the matrix Fortran Kernels to work with single precision
!    matrix data structures
!

!
!     PetscLogDouble variables are used to contain double precision numbers
!     that are not used in the numerical computations, but rather in logging,
!     timing etc.
!


!
!     Macros for error checking
!




















!
!  Include file for Fortran use of the type(tPetscViewer) package in PETSc
!




































































!
!  Matrix types
!


!
! MatMFFDType values
!



!
! MatSolverTypes
!


!
! GPU Storage Formats for CUSPARSE
!



!
! GPU Storage Formats for HIPSPARSE
!



!
! sparsity reducing ordering for STRUMPACK
!



!
!
!  Include file for Fortran use of the type(tVec) package in PETSc
!

!
!
!  Include file for Fortran use of the KSP package in PETSc
!




!
!
!  Include file for Fortran use of the PC (preconditioner) package in PETSc
!




!
!
!  Include file for Fortran use of the type(tMat) package in PETSc
!  Portions of this code are under:
!  Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
!


!
!  Include file for Fortran use of the DM package in PETSc
!




!
!
!  Include file for Fortran use of the type(tIS) (index set) package in PETSc
!

!
!
!  Include file for Fortran use of the type(tVec) package in PETSc
!

!
!
!  Include file for Fortran use of the type(tMat) package in PETSc
!  Portions of this code are under:
!  Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
!



















!
! GAMG types
!



!
! GAMG classical types
!



!
! Various preconditioners
!









!
!  Various Krylov subspace methods
!

!
!  Various Initial guesses for Krylov subspace methods
!

!
!----------------------Type Declarations-------------------------------!
!
      integer(kind=selected_int_kind(5)) :: NLX,NGX,ISTART,IEND,NC
      INTEGER, DIMENSION(:), ALLOCATABLE :: IDX_FROM,IDX_TO
      INTEGER, DIMENSION(:), ALLOCATABLE :: ID_NNZ,IO_NNZ
      integer(kind=selected_int_kind(5)) :: IERR
      type(tMat) :: MATX
      type(tVec) :: SOL_VECX,SOL_VEC_SX
      type(tVec) :: RHS_VECX
      type(tVecScatter) :: SCATTERX
      type(tIS) :: FROMX,TOX
      type(tKSP) :: KSPX
      character*(80) :: KTYPEX
      type(tPC) :: PCX
      character*(80) :: PTYPEX
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/STOMP_PETSC_CREATE'
!
!---  Coupled flow and transport linear system  ---
!
      IF( INDX.EQ.0 ) THEN
        NLX = NUKFL(ID+1)
        NGX = NUKFG
        ALLOCATE( ID_NNZ(1:NUKFL(ID+1)),STAT=ISTAT )
        CHMSG = 'ID_NNZ'
        CALL ALLOC_ERROR( CHMSG,ISTAT )
        ALLOCATE( IO_NNZ(1:NUKFL(ID+1)),STAT=ISTAT )
        CHMSG = 'IO_NNZ'
        CALL ALLOC_ERROR( CHMSG,ISTAT )
        ID_NZ = 0
        IO_NZ = 0
        NC = 0
        IROW1 = 1 + NUKFO(ID+1) - 1
        IROW2 = NUKFL(ID+1) + NUKFO(ID+1) - 1
        DO NEQ = 1,NUKFL(ID+1)
          ID_NNZ(NEQ) = 0
          IO_NNZ(NEQ) = 0
          DO NCOL = NLU(NEQ),(NLU(NEQ+1)-1)
            NC = NC + 1
            IF( MLU(NC).GE.IROW1 .AND. MLU(NC).LE.IROW2 ) THEN
              ID_NNZ(NEQ) = ID_NNZ(NEQ) + 1
            ELSE
              IO_NNZ(NEQ) = IO_NNZ(NEQ) + 1
            ENDIF
          ENDDO
          ID_NZ = MAX( ID_NZ,ID_NNZ(NEQ))
          IO_NZ = MAX( IO_NZ,IO_NNZ(NEQ))
        ENDDO
!
!---  Solute/species transport linear system  ---
!
      ELSEIF( INDX.EQ.1 ) THEN
        NLX = NUKTL(ID+1)
        NGX = NUKTG
        ALLOCATE( ID_NNZ(1:NUKTL(ID+1)),STAT=ISTAT )
        CHMSG = 'ID_NNZ'
        CALL ALLOC_ERROR( CHMSG,ISTAT )
        ALLOCATE( IO_NNZ(1:NUKTL(ID+1)),STAT=ISTAT )
        CHMSG = 'IO_NNZ'
        CALL ALLOC_ERROR( CHMSG,ISTAT )
        ID_NZ = 0
        IO_NZ = 0
        NC = 0
        IROW1 = 1 + NUKTO(ID+1) - 1
        IROW2 = NUKTL(ID+1) + NUKTO(ID+1) - 1
        DO NEQ = 1,NUKTL(ID+1)
          ID_NNZ(NEQ) = 0
          IO_NNZ(NEQ) = 0
          DO NCOL = NLUC(NEQ),(NLUC(NEQ+1)-1)
            NC = NC + 1
            IF( MLUC(NC).GE.IROW1 .AND. MLUC(NC).LE.IROW2 ) THEN
              ID_NNZ(NEQ) = ID_NNZ(NEQ) + 1
            ELSE
              IO_NNZ(NEQ) = IO_NNZ(NEQ) + 1
            ENDIF
          ENDDO
          ID_NZ = MAX( ID_NZ,ID_NNZ(NEQ))
          IO_NZ = MAX( IO_NZ,IO_NNZ(NEQ))
        ENDDO
!
!---  Geomechanics linear system  ---
!
      ELSEIF( INDX.EQ.2 ) THEN
        NLX = NUKGL(ID+1)
        NGX = NUKGG
        ALLOCATE( ID_NNZ(1:NUKGL(ID+1)),STAT=ISTAT )
        CHMSG = 'ID_NNZ'
        CALL ALLOC_ERROR( CHMSG,ISTAT )
        ALLOCATE( IO_NNZ(1:NUKGL(ID+1)),STAT=ISTAT )
        CHMSG = 'IO_NNZ'
        CALL ALLOC_ERROR( CHMSG,ISTAT )
        ID_NZ = 0
        IO_NZ = 0
        NC = 0
        IROW1 = 1 + NUKGO(ID+1) - 1
        IROW2 = NUKGL(ID+1) + NUKGO(ID+1) - 1
        DO NEQ = 1,NUKGL(ID+1)
          ID_NNZ(NEQ) = 0
          IO_NNZ(NEQ) = 0
          DO NCOL = NLU_GM(NEQ),(NLU_GM(NEQ+1)-1)
            NC = NC + 1
            IF( MLU_GM(NC).GE.IROW1 .AND. MLU_GM(NC).LE.IROW2 ) THEN
              ID_NNZ(NEQ) = ID_NNZ(NEQ) + 1
            ELSE
              IO_NNZ(NEQ) = IO_NNZ(NEQ) + 1
            ENDIF
          ENDDO
          ID_NZ = MAX( ID_NZ,ID_NNZ(NEQ))
          IO_NZ = MAX( IO_NZ,IO_NNZ(NEQ))
        ENDDO
      ENDIF
!
!---  Creates a sparse parallel matrix in 'aij' format (the default 
!     parallel PETSc format). For good matrix assembly performance 
!     the user should preallocate the matrix storage by setting the 
!     parameters d_nz (or d_nnz) and o_nz (or o_nnz). By setting these 
!     parameters accurately, performance can be increased by more 
!     than a factor of 50.
!
      CALL MatCreateAIJ( MPI_COMM_WORLD,NLX,NLX,NGX,NGX,
     &  ID_NZ,ID_NNZ,IO_NZ,IO_NNZ,MATX,IERR )
!      PRINT *,'NLX = ',NLX,'NGX = ',NGX,'ID = ',ID
!
!---  Creates a matrix where the type is determined from the options 
!     database. Generates a parallel MPI matrix if the communicator 
!     has more than one processor. The default matrix type is 'aij', 
!     using the routines MatCreateSeqAIJ() and MatCreateAIJ() if you do 
!     not select a type in the options database.
!
      CALL MatSetFromOptions( MATX,IERR )
!
!---  Sets up the internal matrix data structures for later use.
!
      CALL MatSetUp( MATX,IERR )
!!
!!---  Currently, all PETSc parallel matrix formats are partitioned by
!!     contiguous chunks of rows across the processors.  Determine which
!!     rows of the matrix are locally owned.
!!
!      CALL MatGetOwnershipRange( MATX,ISMX,IEMX,IERR )
!      PRINT *,'ISMX = ',ISMX,'IEMX = ',IEMX,'NLX = ',NLX,'NGX = ',NGX,
!     &  'ID = ',ID
!
!---  Create solution (SOL) vector  ---
!
      CALL VecCreate( MPI_COMM_WORLD,SOL_VECX,IERR )
      CALL VecSetSizes( SOL_VECX,NLX,NGX,IERR )
      CALL VecSetFromOptions( SOL_VECX,IERR )
!
!---  Create problem (RHS) vector  ---
!
      CALL VecCreate( MPI_COMM_WORLD,RHS_VECX,IERR )
      CALL VecSetSizes( RHS_VECX,NLX,NGX,IERR )
      CALL VecSetFromOptions( RHS_VECX,IERR )
!
!---  Create a vector scatter object to communicate the solution
!     vector to the field cells, ghost cells, and wells for the
!     coupled flow and transport solution  ---
!
      IF( INDX.EQ.0 ) THEN
        NC = 0
        DO N = 1,NFCGC(ID+1)
          IF( IXP(N).EQ.0 ) CYCLE
          DO M = 1,ISVC
            NC = NC + 1
          ENDDO
        ENDDO
        IF( ID+1.EQ.NP ) THEN
          DO NCW = 1,N_CW
            NC = NC + 1
          ENDDO
        ENDIF
        NUKFP(ID+1) = NC
        ALLOCATE( IDX_FROM(1:NUKFP(ID+1)),STAT=ISTAT )
        CHMSG = 'IDX_FROM'
        CALL ALLOC_ERROR( CHMSG,ISTAT )
        ALLOCATE( IDX_TO(1:NUKFP(ID+1)),STAT=ISTAT )
        CHMSG = 'IDX_TO'
        CALL ALLOC_ERROR( CHMSG,ISTAT )
!        NUKX = NFCGC(ID+1)*ISVC
!        IF( ID+1.EQ.NP ) NUKX = NUKX + N_CW
!        ALLOCATE( IDX_FROM(1:NUKX),STAT=ISTAT )
!        CHMSG = 'IDX_FROM'
!        CALL ALLOC_ERROR( CHMSG,ISTAT )
!        ALLOCATE( IDX_TO(1:NUKX),STAT=ISTAT )
!        CHMSG = 'IDX_TO'
!        CALL ALLOC_ERROR( CHMSG,ISTAT )
        NC = 0
        DO N = 1,NFCGC(ID+1)
          IF( IXP(N).EQ.0 ) CYCLE
          DO M = 1,ISVC
            NC = NC + 1
            IDX_FROM(NC) = (IXP(N)-1)*ISVC + M - 1
            IDX_TO(NC) = NC - 1
          ENDDO
        ENDDO
        IF( ID+1.EQ.NP ) THEN
          DO NCW = 1,N_CW
            NC = NC + 1
            IDX_FROM(NC) = JM_CW(NCW) - 1
            IDX_TO(NC) = NC - 1
          ENDDO
        ENDIF
        CALL VecCreateSeq( PETSC_COMM_SELF,NC,SOL_VEC_SX,IERR )
        CALL ISCreateGeneral( PETSC_COMM_SELF,NC,IDX_FROM,
     &    PETSC_COPY_VALUES,FROMX,IERR )
        CALL ISCreateGeneral( PETSC_COMM_SELF,NC,IDX_TO,
     &    PETSC_COPY_VALUES,TOX,IERR )
        CALL VecScatterCreate( SOL_VECX,FROMX,SOL_VEC_SX,TOX,
     &    SCATTERX,IERR)
      ENDIF
!
!---  Create solver  ---
!
      CALL KSPCreate( MPI_COMM_WORLD,KSPX,IERR )
!
!---  Set the options for the solver  ---
!
      CALL KSPSetOperators( KSPX,MATX,MATX,IERR )
      KTYPEX = 'bicg'
      CALL KSPSetType( KSPX,KTYPEX,IERR )
!      PRINT *,'RTOL = ',RTOL_PETSC,'ATOL = ',ATOL_PETSC,
!     &  'DTOL = ',DTOL_PETSC,'MAXIT = ',MAXITS_PETSC,'ID = ',ID
      CALL KSPSetTolerances( KSPX,RTOL_PETSC,ATOL_PETSC,DTOL_PETSC,
     &  MAXITS_PETSC,IERR )
!
!---  Deallocate temporary memory  ---
!
      DEALLOCATE( ID_NNZ,STAT=ISTAT )
      CHMSG = 'ID_NNZ'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
      DEALLOCATE( IO_NNZ,STAT=ISTAT )
      CHMSG = 'IO_NNZ'
      CALL DEALLOC_ERROR( CHMSG,ISTAT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of STOMP_PETSC_CREATE group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE STOMP_PETSC_DESTROY( KSPX,MATX,RHS_VECX,SOL_VECX,
     &  SOL_VEC_SX,SCATTERX,INDX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Destroy solution vector, problem vector, matrix, solver, and
!     solution-vector scatter scheme for the PETSc solver
!
!     INDX = 0 : Coupled flow and transport
!     INDX = 1 : Solute or species component transport
!     INDX = 2 : Geomechanics
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 7 December 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE GRID
      USE COUP_WELL
      USE PETSCKSP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

!
!
!  Include file for Fortran use of the type(tMat) package in PETSc
!  Portions of this code are under:
!  Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
!

!
!
!  Include file for Fortran use of the type(tVec) package in PETSc
!

!
!
!  Include file for Fortran use of the type(tKSP) package in PETSc
!

!
!----------------------Type Declarations-------------------------------!
!
      integer(kind=selected_int_kind(5)) :: IERR
      type(tMat) :: MATX
      type(tVec) :: SOL_VECX,SOL_VEC_SX
      type(tVec) :: RHS_VECX
      type(tVecScatter) :: SCATTERX
      type(tKSP) :: KSPX
      CHARACTER(32) :: CHMSG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/STOMP_PETSC_DESTROY'
!
!---  Destroy parallel matrix
!
      CALL MatDestroy( MATX,IERR )
!
!---  Destroy solution (SOL) vector  ---
!
      CALL VecDestroy( SOL_VECX,IERR )
!
!---  Destroy problem (RHS) vector  ---
!
      CALL VecDestroy( RHS_VECX,IERR )
!
!---  Destroy (SOL) scatter vector  ---
!
      CALL VecDestroy( SOL_VEC_SX,IERR )
!
!---  Destroy the (SOL) vector scatter index set  ---
!
      CALL VecScatterDestroy( SCATTERX,IERR )
!
!---  Destroy type(tKSP) context  ---
!
      CALL KSPDestroy( KSPX,IERR )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of STOMP_PETSC_DESTROY group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE STOMP_PETSC_SOLVE( KSPX,MATX,RHS_VECX,SOL_VECX,
     &  SOL_VEC_SX,SCATTERX,INDX )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Solve the linear system A x = b with PETSc
!
!     INDX = 0 : Coupled flow and transport
!     INDX = 1 : Solute or species component transport
!     INDX = 2 : Geomechanics
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 December 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE PETSCKSP
      USE JACOB
      USE GRID
      USE GEO_MECH
      USE COUP_WELL
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

!
!
!  Include file for Fortran use of the type(tMat) package in PETSc
!  Portions of this code are under:
!  Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
!

!
!
!  Include file for Fortran use of the type(tVec) package in PETSc
!

!
!
!  Include file for Fortran use of the type(tKSP) package in PETSc
!

!
!----------------------Type Declarations-------------------------------!
!
      integer(kind=selected_int_kind(5)) :: IERR
      real(kind=selected_real_kind(10)), POINTER :: SOL_VEC_ARRAY(:)
      type(tMat) :: MATX
      type(tVec) :: SOL_VECX,SOL_VEC_SX
      type(tVec) :: RHS_VECX
      type(tVecScatter) :: SCATTERX
      type(tKSP) :: KSPX
      character*(80) :: KTYPEX
      type(tPC) :: PCX
      character*(80) :: PTYPEX
      integer(kind=selected_int_kind(10)) :: IOFFSETX
      type(tPetscViewer) :: VIEWERX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/STOMP_PETSC_SOLVE'
!
!---  Assemble the matrix, using the 2-step process
!
      CALL MatAssemblyBegin( MATX,MAT_FINAL_ASSEMBLY,IERR )
      CALL MatAssemblyEnd( MATX,MAT_FINAL_ASSEMBLY,IERR )
!      PRINT *,'MatAssemblyEnd: ID = ',ID
!
!---  Assemble the problem vector, using the 2-step process
!
      CALL VecAssemblyBegin( RHS_VECX,IERR )
      CALL VecAssemblyEnd( RHS_VECX,IERR )
!      PRINT *,'VecAssemblyEnd: RHS_VECX: ID = ',ID
!
!---  Output coupled-flow problem vector  ---
!
      IF( ISLC(34).EQ.1 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER ) THEN
        CALL PetscViewerASCIIOpen( MPI_COMM_WORLD,"frhs_vecx.dat",
     &    VIEWERX,IERR )
        CALL VecView( RHS_VECX,VIEWERX,IERR )
        ICNV = 4
        IF( ID.EQ.0 ) PRINT *,'Output: frhs_vecx.dat'
        RETURN
      ENDIF
!
!---  Output coupled-flow matrix  ---
!
      IF( ISLC(34).EQ.2 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER ) THEN
        CALL PetscViewerASCIIOpen( MPI_COMM_WORLD,"fmatx.dat",
     &    VIEWERX,IERR )
        CALL MatView( MATX,VIEWERX,IERR )
        ICNV = 4
        IF( ID.EQ.0 ) PRINT *,'Output: fmatx.dat'
        RETURN
      ENDIF
!
!---  Output transport problem vector  ---
!
      IF( ISLC(34).EQ.11 .AND. INDX.EQ.1 .AND. ISLC(83).LE.NSTEP ) THEN
        CALL PetscViewerASCIIOpen( MPI_COMM_WORLD,"crhs_vecx.dat",
     &    VIEWERX,IERR )
        CALL VecView( RHS_VECX,VIEWERX,IERR )
        ICNV = 4
        IF( ID.EQ.0 ) PRINT *,'Output: crhs_vecx.dat'
        RETURN
      ENDIF
!
!---  Output transport matrix  ---
!
      IF( ISLC(34).EQ.12 .AND. INDX.EQ.1 .AND. ISLC(83).LE.NSTEP ) THEN
        CALL PetscViewerASCIIOpen( MPI_COMM_WORLD,"cmatx.dat",
     &    VIEWERX,IERR )
        CALL MatView( MATX,VIEWERX,IERR )
        ICNV = 4
        IF( ID.EQ.0 ) PRINT *,'Output: cmatx.dat'
        RETURN
      ENDIF
!
!---  Output geomechanics problem vector  ---
!
      IF( ISLC(34).EQ.21 .AND. INDX.EQ.2 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER .AND. ISLC(85).EQ.IREF_GM ) THEN
        CALL PetscViewerASCIIOpen( MPI_COMM_WORLD,"grhs_vecx.dat",
     &    VIEWERX,IERR )
        CALL VecView( RHS_VECX,VIEWERX,IERR )
        ICNV = 4
        IF( ID.EQ.0 ) PRINT *,'Output: grhs_vecx.dat'
        RETURN
      ENDIF
!
!---  Output geomechanics matrix  ---
!
      IF( ISLC(34).EQ.22 .AND. INDX.EQ.2 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER .AND. ISLC(85).EQ.IREF_GM ) THEN
        CALL PetscViewerASCIIOpen( MPI_COMM_WORLD,"gmatx.dat",
     &    VIEWERX,IERR )
        CALL MatView( MATX,VIEWERX,IERR )
        ICNV = 4
        IF( ID.EQ.0 ) PRINT *,'Output: gmatx.dat'
        RETURN
      ENDIF
!
!---  Solve linear system A x = b
!
      CALL KSPSolve( KSPX,RHS_VECX,SOL_VECX,IERR )
      CALL KSPGetIterationNumber( KSPX,ITSX,IERR )
      IF( ITSX.GE.MAXITS_PETSC ) ICNV = 1
!      IF( ID.EQ.0 ) PRINT *,'ITSX = ',ITSX
!
!---  Output coupled-flow solution vector  ---
!
      IF( ISLC(34).EQ.3 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER ) THEN
        CALL PetscViewerASCIIOpen( MPI_COMM_WORLD,"fsol_vecx.dat",
     &    VIEWERX,IERR )
        CALL VecView( SOL_VECX,VIEWERX,IERR )
        ICNV = 4
        IF( ID.EQ.0 ) PRINT *,'Output: fsol_vecx.dat'
        RETURN
      ENDIF
!
!---  Output transport solution vector  ---
!
      IF( ISLC(34).EQ.13 .AND. INDX.EQ.1 .AND. ISLC(83).LE.NSTEP ) THEN
        CALL PetscViewerASCIIOpen( MPI_COMM_WORLD,"csol_vecx.dat",
     &    VIEWERX,IERR )
        CALL VecView( SOL_VECX,VIEWERX,IERR )
        ICNV = 4
        IF( ID.EQ.0 ) PRINT *,'Output: csol_vecx.dat'
        RETURN
      ENDIF
!
!---  Output geomechanics solution vector  ---
!
      IF( ISLC(34).EQ.23 .AND. INDX.EQ.2 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER .AND. ISLC(85).EQ.IREF_GM ) THEN
        CALL PetscViewerASCIIOpen( MPI_COMM_WORLD,"gsol_vecx.dat",
     &    VIEWERX,IERR )
        CALL VecView( SOL_VECX,VIEWERX,IERR )
        ICNV = 4
        IF( ID.EQ.0 ) PRINT *,'Output: gsol_vecx.dat'
        RETURN
      ENDIF
!
!---  Scatter solution vector
!
!      PRINT *,'SCATTERX = ',SCATTERX,'ID = ',ID
      CALL VecScatterBegin( SCATTERX,SOL_VECX,SOL_VEC_SX,INSERT_VALUES,
     &  SCATTER_FORWARD,IERR )
      CALL VecScatterEnd( SCATTERX,SOL_VECX,SOL_VEC_SX,INSERT_VALUES,
     &  SCATTER_FORWARD,IERR )
      CALL VecGetArrayReadF90( SOL_VEC_SX,SOL_VEC_ARRAY,IERR )
!
!---  Coupled flow and heat transport solution
!
      IF( INDX.EQ.0 ) THEN
!        NUKX = NFCGC(ID+1)*ISVC
!        IF( ID+1.EQ.NP ) NUKX = NUKX + N_CW
!        DO M = 1,NUKX
!          BLU(M) = SOL_VEC_ARRAY(M)
!        ENDDO
        DO M = 1,NUKFP(ID+1)
          BLU(M) = SOL_VEC_ARRAY(M)
        ENDDO
!
!---  Solute/species transport solution
!
      ELSEIF( INDX.EQ.1 ) THEN
!        NUKX = NFCGC(ID+1)
!        DO M = 1,NUKX
!          BLU(M) = SOL_VEC_ARRAY(M)
!        ENDDO
        DO M = 1,NUKTP(ID+1)
          BLU(M) = SOL_VEC_ARRAY(M)
        ENDDO
!
!---  Geomechanics solution
!
      ELSEIF( INDX.EQ.2 ) THEN
!        NUKX = NFNGN(ID+1)
!        DO M = 1,NUKX
!          BLU(M) = SOL_VEC_ARRAY(M)
!        ENDDO
        DO M = 1,NUKGP(ID+1)
          BLU(M) = SOL_VEC_ARRAY(M)
        ENDDO
      ENDIF
      CALL VecRestoreArrayReadF90( SOL_VEC_SX,SOL_VEC_ARRAY,IERR )
!
!---  Initialize problem vector for next solve  ---
!
      CALL VecSet( RHS_VECX,0.D+0,IERR )
      IF( IERR.NE.0 ) THEN
        I_ERR(1) = 0
        I_ERR(2) = 0
        I_ERR(3) = 0
        I_ERR(4) = ID
        IF( INDX.EQ.0 ) THEN
          M_ERR(1) = 'Coupled Flow and Transport PETSc VecSet'
        ELSEIF( INDX.EQ.1 ) THEN
          M_ERR(1) = 'Solute Transport  PETSc VecSet'
        ELSEIF( INDX.EQ.2 ) THEN
          M_ERR(1) = 'Geomechanics  PETSc VecSet'
        ELSE
          M_ERR(1) = ' PETSc VecSet'
        ENDIF
        RETURN
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of STOMP_PETSC_SOLVE group
!
      RETURN
      END



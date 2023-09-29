!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE STOMP_PETSC_CREATE( KSPX,MATX,RHS_VECX,SOL_VECX,INDX )
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
!     Create solution vector, problem vector, matrix, and solver
!     for the PETSc solver
!
!     INDX = 0 : Coupled flow and transport
!     INDX = 1 : Solute or species component transport
!     INDX = 2 : Geomechanics
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 December 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_CL
      USE SOLTN
      USE LEAK_WELL
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
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
      integer(kind=selected_int_kind(5)) :: NGX
      INTEGER, DIMENSION(:), ALLOCATABLE :: NNZ
      integer(kind=selected_int_kind(5)) :: IERR
      type(tMat) :: MATX
      type(tVec) :: SOL_VECX
      type(tVec) :: RHS_VECX
      type(tKSP) :: KSPX
      character*(80) :: KTYPEX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/STOMP_PETSC_CREATE'
!
!---  Number of unknowns for coupled flow and heat transport  ---
!
      IF( INDX.EQ.0 ) THEN
!
!---  Number of unknowns equals number of number of coupled equations x
!     (active nodes + fracture triangles) + number of coupled wells  ---
!
        IF( N_CW.GT.0 ) THEN
          NGX = ISVC*(NFBN - NRFN - NXP + NWN_LW + NT_FRC - NXP_FRC) + 
     &      N_CW
!
!---    Number of unknowns equals number of number of coupled equations
!       x (number of active nodes + number of wells +
!       number of surface nodes)
!
        ELSEIF( LSPILL.GT.0 ) THEN
          NGX = ISVC*(NFBN - NRFN - NXP + NWLN + NWN_LW + IJFLD)
!
!---    Number of unknowns equals number of number of coupled equations
!       x (number of active nodes + number of wells + fracture triangles
!       + number of borehole nodes)
!
        ELSE
          NGX = ISVC*(NFBN - NRFN - NXP + NWLN + NWN_LW + NT_FRC - 
     &      NXP_FRC + NBN_BH - NXP_BH)
        ENDIF
!
!---    Dual-porosity option  ---
!
        IF( ISLC(11).EQ.1 ) NGX = NGX + ISVC*(NFBN-NRFN-NXP)
!
!---  Number of unknowns equals number of active nodes for
!     transport equations ---
!
      ELSEIF( INDX.EQ.1 ) THEN
        NGX = NFBN - NRFN - NXP + NWN_LW + NT_FRC - NXP_FRC + 
     &    NBN_BH - NXP_BH
!
!---  Geomechanics - number of unknowns equals number of active
!     finite element nodes x 3
!
      ELSEIF( INDX.EQ.2 ) THEN
        NGX = 3*NFEN_GM
      ENDIF 
!
!---  Array containing the number of non-zeros per equation row  ---
!
      ALLOCATE( NNZ(1:NGX),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: NNZ'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Number of non-zeros for coupled flow and heat transport  ---
!
      IF( INDX.EQ.0 ) THEN
        NNZFR = 0
        DO N = 1,NGX
          NNZ(N) = NLU(N+1)-NLU(N)
          NNZFR = MAX( NNZ(N),NNZFR )
        ENDDO
!
!---  Number of non-zeros for solute/species transport  ---
!
      ELSEIF( INDX.EQ.1 ) THEN
        NNZTR = 0
        DO N = 1,NGX
          NNZ(N) = NLUC(N+1)-NLUC(N)
          NNZTR = MAX( NNZ(N),NNZTR )
        ENDDO
!
!---  Number of non-zeros for geomechanics  ---
!
      ELSEIF( INDX.EQ.2 ) THEN
        NNZR_GM = 0
        DO N = 1,NGX
          NNZ(N) = NLU_GM(N+1)-NLU_GM(N)
          NNZR_GM = MAX( NNZ(N),NNZR_GM )
        ENDDO
      ENDIF
!
!---  Create a sparse matrix in compressed row format (the default
!     parallel PETSc format
!     
      CALL MatCreateSeqAIJ( PETSC_COMM_SELF,NGX,NGX,
     &  PETSC_DEFAULT_INTEGER,NNZ,MATX,IERR )
      CALL MatSetFromOptions( MATX,IERR )
      CALL MatSetUp( MATX,IERR )
!
!---  Initialize matrix
!
      CALL MatZeroEntries( MATX,IERR ) 
!
!---  Create solution (SOL) vector  ---
!
      CALL VecCreateSeq( PETSC_COMM_SELF,NGX,SOL_VECX,IERR )
!
!---  Initialize solution (SOL) vector  ---
!
      CALL VecZeroEntries( SOL_VECX,IERR ) 
!
!---  Create problem (RHS) vector  ---
!
      CALL VecCreateSeq( PETSC_COMM_SELF,NGX,RHS_VECX,IERR )
!
!---  Initialize problem (RHS) vector  ---
!
      CALL VecZeroEntries( RHS_VECX,IERR ) 
!
!---  Create solver  ---
!
      CALL KSPCreate( PETSC_COMM_SELF,KSPX,IERR )
!
!---  Set the options for the solver  ---
!
      CALL KSPSetOperators( KSPX,MATX,MATX,IERR )
      KTYPEX = 'ibcgs'
      CALL KSPSetType( KSPX,KTYPEX,IERR )
      CALL KSPSetTolerances( KSPX,RTOL_PETSC,ATOL_PETSC,DTOL_PETSC,
     &  MAXITS_PETSC,IERR )
!
!---  Deallocation number of non-zeros array  ---
!
      DEALLOCATE( NNZ,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: NNZ'
        CALL WRMSGS( INDX )
      ENDIF
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
      SUBROUTINE STOMP_PETSC_DESTROY( KSPX,MATX,RHS_VECX,SOL_VECX )
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
!     Destroy solution vector, problem vector, matrix, and solver
!     for the PETSc solver
!
!     INDX = 0 : Coupled flow and transport
!     INDX = 1 : Solute or species component transport
!     INDX = 2 : Geomechanics
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 December 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE SOLTN
      USE GRID
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
      type(tVec) :: SOL_VECX
      type(tVec) :: RHS_VECX
      type(tKSP) :: KSPX
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
      SUBROUTINE STOMP_PETSC_SOLVE( KSPX,MATX,RHS_VECX,SOL_VECX,INDX )
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
!     Written by M.D. White, PNNL, 8 December 2022
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_CL
      USE SOLTN
      USE LEAK_WELL
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE GEO_MECH
      USE FILES
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
!
!  Part of the base include file for Fortran use of PETSc.
!  Note: This file should contain only define statements and
!  not the declaration of variables.

! No spaces for #defines as some compilers (PGI) also adds
! those additional spaces during preprocessing - bad for fixed format
!

!
!----------------------Type Declarations-------------------------------!
!
      integer(kind=selected_int_kind(5)) :: IERR
      real(kind=selected_real_kind(10)) :: SOL_VEC_ARRAY(1)
      type(tMat) :: MATX
      type(tVec) :: SOL_VECX
      type(tVec) :: RHS_VECX
      type(tKSP) :: KSPX
      integer(kind=selected_int_kind(10)) :: IOFFSETX
      type(tPetscViewer) :: VIEWERX
      REAL*8, DIMENSION(:), ALLOCATABLE :: BUFFER
      INTEGER, DIMENSION(:), ALLOCATABLE :: ICOL
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/STOMP_PETSC_SOLVE'
!
!---  Number of unknowns for coupled flow and heat transport  ---
!
      IF( INDX.EQ.0 ) THEN
        ALLOCATE( BUFFER(1:NNZFR),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: BUFFER'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( ICOL(1:NNZFR),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: ICOL'
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Number of unknowns equals number of number of coupled 
!       equations x (active nodes + fracture triangles) + number of 
!       coupled wells  ---
!
        IF( N_CW.GT.0 ) THEN
          NGX = ISVC*(NFBN - NRFN - NXP + NWN_LW + NT_FRC - NXP_FRC) + 
     &      N_CW
!
!---    Number of unknowns equals number of number of coupled equations
!       x (number of active nodes + number of wells +
!       number of surface nodes)
!
        ELSEIF( LSPILL.GT.0 ) THEN
          NGX = ISVC*(NFBN - NRFN - NXP + NWLN + NWN_LW + IJFLD)
!
!---    Number of unknowns equals number of number of coupled equations
!       x (number of active nodes + number of wells + fracture triangles
!       + number of borehole nodes)
!
        ELSE
          NGX = ISVC*(NFBN - NRFN - NXP + NWLN + NWN_LW + NT_FRC - 
     &      NXP_FRC + NBN_BH - NXP_BH)
        ENDIF
!
!---    Dual-porosity option  ---
!
        IF( ISLC(11).EQ.1 ) NGX = NGX + ISVC*(NFBN-NRFN-NXP)
!
!---  Number of unknowns equals number of active nodes for
!     transport equations ---
!
      ELSEIF( INDX.EQ.1 ) THEN
        ALLOCATE( BUFFER(1:NNZTR),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: BUFFER'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( ICOL(1:NNZTR),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: ICOL'
          CALL WRMSGS( INDX )
        ENDIF
        NGX = NFBN - NRFN - NXP + NWN_LW + NT_FRC - NXP_FRC + 
     &    NBN_BH - NXP_BH
!
!---  Geomechanics - number of unknowns equals number of active
!     finite element nodes x 3
!
      ELSEIF( INDX.EQ.2 ) THEN
        ALLOCATE( BUFFER(1:NNZR_GM),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: BUFFER'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( ICOL(1:NNZR_GM),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: ICOL'
          CALL WRMSGS( INDX )
        ENDIF
        NGX = 3*NFEN_GM
      ENDIF 
!
!---  Load PETSc problem vector  ---
!
      CALL VecSetValues( RHS_VECX,NGX,ILU,BLU,INSERT_VALUES,IERR )
!
!---  Load PETSc matrix for coupled flow and heat transport  ---
!
      IF( INDX.EQ.0 ) THEN
        NC = 0
        DO NEQ = 1,NGX
          MC = 0
          DO NCOL = NLU(NEQ),(NLU(NEQ+1)-1)
            NC = NC + 1
            MC = MC + 1
            BUFFER(MC) = DLU(NC)
            ICOL(MC) = MLU(NC)
          ENDDO
          IROW = NEQ - 1
          CALL MatSetValues( MATX,1,IROW,MC,ICOL,BUFFER,
     &      ADD_VALUES,IERR )
        ENDDO
!
!---  Load PETSc matrix for solute/species transport  ---
!
      ELSEIF( INDX.EQ.1 ) THEN
        NC = 0
        DO NEQ = 1,NGX
          MC = 0
          DO NCOL = NLUC(NEQ),(NLUC(NEQ+1)-1)
            NC = NC + 1
            MC = MC + 1
            BUFFER(MC) = DLU(NC)
            ICOL(MC) = MLUC(NC)
          ENDDO
          IROW = NEQ - 1
          CALL MatSetValues( MATX,1,IROW,MC,ICOL,BUFFER,
     &      ADD_VALUES,IERR )
        ENDDO
!
!---  Load PETSc matrix for geomechanics  ---
!
      ELSEIF( INDX.EQ.1 ) THEN
        NC = 0
        DO NEQ = 1,NGX
          MC = 0
          DO NCOL = NLU_GM(NEQ),(NLU_GM(NEQ+1)-1)
            NC = NC + 1
            MC = MC + 1
            BUFFER(MC) = DLU(NC)
            ICOL(MC) = MLU_GM(NC)
          ENDDO
          IROW = NEQ - 1
          CALL MatSetValues( MATX,1,IROW,MC,ICOL,BUFFER,
     &      ADD_VALUES,IERR )
        ENDDO
      ENDIF
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
        CALL PetscViewerASCIIOpen( PETSC_COMM_SELF,"frhs_vec.dat",
     &    VIEWERX,IERR )
        CALL PetscViewerPushFormat( VIEWERX,PETSC_VIEWER_ASCII_DETAIL,
     &    IERR )
        CALL VecView( RHS_VECX,VIEWERX,IERR )
        WRITE(ISC,'(/,A)') ' Output: frhs_vec.dat'
        STOP
      ENDIF
!
!---  Output coupled-flow matrix  ---
!
      IF( ISLC(34).EQ.2 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER ) THEN
        CALL PetscViewerASCIIOpen( PETSC_COMM_SELF,"fmat.dat",
     &    VIEWERX,IERR )
        CALL PetscViewerPushFormat( VIEWERX,PETSC_VIEWER_ASCII_DETAIL,
     &    IERR )
        CALL MatView( MATX,VIEWERX,IERR )
        WRITE(ISC,'(/,A)') ' Output: fmat.dat'
        STOP
      ENDIF
!
!---  Output transport problem vector  ---
!
      IF( ISLC(34).EQ.11 .AND. INDX.EQ.1 .AND. ISLC(83).LE.NSTEP ) THEN
        CALL PetscViewerASCIIOpen( PETSC_COMM_SELF,"crhs_vec.dat",
     &    VIEWERX,IERR )
        CALL PetscViewerPushFormat( VIEWERX,PETSC_VIEWER_ASCII_DETAIL,
     &    IERR )
        CALL VecView( RHS_VECX,VIEWERX,IERR )
        WRITE(ISC,'(/,A)') ' Output: crhs_vec.dat'
        STOP
      ENDIF
!
!---  Output transport matrix  ---
!
      IF( ISLC(34).EQ.12 .AND. INDX.EQ.1 .AND. ISLC(83).LE.NSTEP ) THEN
        CALL PetscViewerASCIIOpen( PETSC_COMM_SELF,"cmat.dat",
     &    VIEWERX,IERR )
        CALL PetscViewerPushFormat( VIEWERX,PETSC_VIEWER_ASCII_DETAIL,
     &    IERR )
        CALL MatView( MATX,VIEWERX,IERR )
        WRITE(ISC,'(/,A)') ' Output: cmat.dat'
        STOP
      ENDIF
!
!---  Output geomechanics problem vector  ---
!
      IF( ISLC(34).EQ.21 .AND. INDX.EQ.2 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER .AND. ISLC(85).EQ.IREF_GM ) THEN
        CALL PetscViewerASCIIOpen( PETSC_COMM_SELF,"grhs_vec.dat",
     &    VIEWERX,IERR )
        CALL PetscViewerPushFormat( VIEWERX,PETSC_VIEWER_ASCII_DETAIL,
     &    IERR )
        CALL VecView( RHS_VECX,VIEWERX,IERR )
        WRITE(ISC,'(/,A)') ' Output: grhs_vec.dat'
        STOP
      ENDIF
!
!---  Output geomechanics matrix  ---
!
      IF( ISLC(34).EQ.22 .AND. INDX.EQ.2 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER .AND. ISLC(85).EQ.IREF_GM ) THEN
        CALL PetscViewerASCIIOpen( PETSC_COMM_SELF,"gmat.dat",
     &    VIEWERX,IERR )
        CALL PetscViewerPushFormat( VIEWERX,PETSC_VIEWER_ASCII_DETAIL,
     &    IERR )
        CALL MatView( MATX,VIEWERX,IERR )
        WRITE(ISC,'(/,A)') ' Output: gmat.dat'
        STOP
      ENDIF
!
!---  Solve linear system A x = b
!
      CALL KSPSolve( KSPX,RHS_VECX,SOL_VECX,IERR )
!      CALL KSPGetIterationNumber( KSPX,ITSX,IERR )
!      PRINT *,'ITSX = ',ITSX
!
!---  Output coupled-flow solution vector  ---
!
      IF( ISLC(34).EQ.3 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER ) THEN
        CALL PetscViewerASCIIOpen( PETSC_COMM_SELF,"fsol_vec.dat",
     &    VIEWERX,IERR )
        CALL PetscViewerPushFormat( VIEWERX,PETSC_VIEWER_ASCII_DETAIL,
     &    IERR )
        CALL VecView( SOL_VECX,VIEWERX,IERR )
        WRITE(ISC,'(/,A)') ' Output: fsol_vec.dat'
        STOP
      ENDIF
!
!---  Output transport solution vector  ---
!
      IF( ISLC(34).EQ.13 .AND. INDX.EQ.1 .AND. ISLC(83).LE.NSTEP ) THEN
        CALL PetscViewerASCIIOpen( PETSC_COMM_SELF,"csol_vec.dat",
     &    VIEWERX,IERR )
        CALL PetscViewerPushFormat( VIEWERX,PETSC_VIEWER_ASCII_DETAIL,
     &    IERR )
        CALL VecView( SOL_VECX,VIEWERX,IERR )
        WRITE(ISC,'(/,A)') ' Output: csol_vec.dat'
        STOP
      ENDIF
!
!---  Output geomechanics solution vector  ---
!
      IF( ISLC(34).EQ.23 .AND. INDX.EQ.2 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER .AND. ISLC(85).EQ.IREF_GM ) THEN
        CALL PetscViewerASCIIOpen( PETSC_COMM_SELF,"gsol_vec.dat",
     &    VIEWERX,IERR )
        CALL PetscViewerPushFormat( VIEWERX,PETSC_VIEWER_ASCII_DETAIL,
     &    IERR )
        CALL VecView( SOL_VECX,VIEWERX,IERR )
        WRITE(ISC,'(/,A)') ' Output: gsol_vec.dat'
        STOP
      ENDIF
!
!---  Get the solution vector values
!
      CALL VecGetArray( SOL_VECX,SOL_VEC_ARRAY,IOFFSETX,IERR )
!
!---  Load the BLU array with the solution values
!
      DO M = 1,NGX
        BLU(M) = SOL_VEC_ARRAY(IOFFSETX+M)
      ENDDO
!
!---  Restore the solution vector values
!
      CALL VecRestoreArray( SOL_VECX,SOL_VEC_ARRAY,IOFFSETX,IERR )
!
!---  Initialize problem vector for next solve  ---
!
      CALL MatZeroEntries( MATX,IERR ) 
      CALL VecZeroEntries( RHS_VECX,IERR ) 
!
!---  Deallocate temporary matrix row arrays  ---
!
      IF( ALLOCATED(BUFFER) ) THEN
        DEALLOCATE( BUFFER,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: BUFFER'
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
      IF( ALLOCATED(ICOL) ) THEN
        DEALLOCATE( ICOL,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: ICOL'
          CALL WRMSGS( INDX )
        ENDIF
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



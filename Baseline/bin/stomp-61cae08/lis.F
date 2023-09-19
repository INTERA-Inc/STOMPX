!---------------------Fortran 90 Module--------------------------------!
!
      MODULE STOMP_LIS_MODULE
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
!     Lis solver module.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 11 December 2015.
!
!---------------------Type Declarations--------------------------!
!
#include "include/lisf.h"

      LIS_MATRIX, PUBLIC :: F_MAT
      LIS_MATRIX, PUBLIC :: T_MAT
      LIS_MATRIX, PUBLIC :: G_MAT
      LIS_VECTOR, PUBLIC :: F_SOL_VEC
      LIS_VECTOR, PUBLIC :: T_SOL_VEC
      LIS_VECTOR, PUBLIC :: G_SOL_VEC
      LIS_VECTOR, PUBLIC :: F_RHS_VEC
      LIS_VECTOR, PUBLIC :: T_RHS_VEC
      LIS_VECTOR, PUBLIC :: G_RHS_VEC
      LIS_SOLVER, PUBLIC :: F_KSP
      LIS_SOLVER, PUBLIC :: T_KSP
      LIS_SOLVER, PUBLIC :: G_KSP
      
      PUBLIC :: STOMP_LIS_CREATE,STOMP_LIS_DESTROY,STOMP_LIS_INIT,
     &          STOMP_LIS_SOLVE


      CONTAINS

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE STOMP_LIS_CREATE( ISLV,KSP,MAT,RHS_VEC,SOL_VEC,INDX )
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
!     Create solution vector, problem vector and matrix for Lis solver
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_CL
      USE SOLTN
      USE LEAK_WELL
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE GEO_MECH
      USE COUP_WELL
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      LIS_INTEGER :: NGX,NLX,IERR
      LIS_MATRIX :: MAT
      LIS_VECTOR :: SOL_VEC
      LIS_VECTOR :: RHS_VEC
      LIS_SOLVER :: KSP
!
!----------------------Executable Lines--------------------------------!
!
      NLX = 0
!
!---  Number of unknowns equals number of active nodes for
!     transport equations ---
!
      IF( INDX.EQ.1 ) THEN
        NGX = NFBN - NRFN - NXP + NWN_LW + NT_FRC - NXP_FRC + 
     &    NBN_BH - NXP_BH
!
!---  Geomechanics - number of unknowns equals number of active
!     finite element nodes x 3
!
      ELSEIF( INDX.EQ.2 ) THEN
        NGX = 3*NFEN_GM
!
!---  Number of unknowns equals number of number of coupled equations x
!     (active nodes + fracture triangles) + number of coupled wells  ---
!
      ELSEIF( N_CW.GT.0 ) THEN
        NGX = ISLV*(NFBN - NRFN - NXP + NWN_LW + NT_FRC - NXP_FRC) + 
     &    N_CW
!
!---  Number of unknowns equals number of number of coupled equations
!     x (number of active nodes + number of wells +
!     number of surface nodes)
!
      ELSEIF( LSPILL.GT.0 ) THEN
        NGX = ISLV*(NFBN - NRFN - NXP + NWLN + NWN_LW + IJFLD)
!
!---  Number of unknowns equals number of number of coupled equations
!     x (number of active nodes + number of wells + fracture triangles
!     + number of borehole nodes)
!
      ELSE
        NGX = ISLV*(NFBN - NRFN - NXP + NWLN + NWN_LW + NT_FRC - 
     &    NXP_FRC + NBN_BH - NXP_BH)
      ENDIF
!
!---  Dual-porosity option  ---
!
      IF( ISLC(11).EQ.1 ) NGX = NGX + ISLV*(NFBN-NRFN-NXP)
!
!---  Create Jacobian matrix  ---
!
      CALL lis_matrix_create(LIS_COMM_WORLD,MAT,IERR)
      CALL lis_matrix_set_size(MAT,NLX,NGX,IERR)
!
!---  Create solution vector  ---
!
      CALL lis_vector_duplicate(MAT,SOL_VEC,IERR)
!
!---  Create problem (RHS) vector  ---
!
      CALL lis_vector_duplicate(MAT,RHS_VEC,IERR)
!
!---  Create solver  ---
!
      CALL lis_solver_create(KSP,IERR)
!
!---  Set the options for the solver  ---
!
      CALL lis_solver_set_option('-i bicgstab',KSP,IERR)
      CALL lis_solver_set_option('-p ilu -ilu_fill 0',KSP,IERR)
      IF( INDX.EQ.2 ) 
     &  CALL lis_solver_set_option('-tol 1.0e-7',KSP,IERR)
      CALL lis_solver_set_option('-print none -adds true',KSP,IERR)
      CALL lis_solver_set_optionC(KSP,IERR)
!
!---  End of STOMP_LIS_CREATE group
!
      RETURN
      END SUBROUTINE STOMP_LIS_CREATE

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE STOMP_LIS_DESTROY( KSP,MAT,RHS_VEC,SOL_VEC )
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
!     Destroy the Lis solver, matrix, and solution and problem vectors
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      LIS_INTEGER :: IERR
      LIS_MATRIX :: MAT
      LIS_VECTOR :: SOL_VEC
      LIS_VECTOR :: RHS_VEC
      LIS_SOLVER :: KSP
!
!----------------------Executable Lines--------------------------------!
!
!
!---  Unset Jacobian matrix from DLU array  ---
!
      CALL lis_matrix_unset(MAT,IERR)
!
!---  Destroy Jacobian matrix  ---
!
      CALL lis_matrix_destroy(MAT,IERR)
!
!---  Destroy solution vector  ---
!
      CALL lis_vector_destroy(SOL_VEC,IERR)
!
!---  Destroy problem vector  ---
!
      CALL lis_vector_destroy(RHS_VEC,IERR)
!
!---  Destroy solver  ---
!
      CALL lis_solver_destroy(KSP,IERR)
!
!---  End of STOMP_LIS_DESTROY group
!
      RETURN
      END SUBROUTINE STOMP_LIS_DESTROY

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE STOMP_LIS_SOLVE( ISLV,KSP,MAT,RHS_VEC,SOL_VEC,INDX )
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
!     Solve the linear system
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
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#ifdef lis
#include "include/lisf.h"
#endif
!
!----------------------Type Declarations-------------------------------!
!
      LIS_INTEGER :: ISX,ICX,NGX,IERR,NNZ
      LIS_MATRIX :: MAT
      LIS_VECTOR :: SOL_VEC
      LIS_VECTOR :: RHS_VEC
      LIS_SOLVER :: KSP
      INTEGER, EXTERNAL :: VALUERSS
!
!----------------------Executable Lines--------------------------------!
!
      IERR = 0
      NLX = 0
!
!---  Number of unknowns equals number of active nodes for
!     transport equations ---
!
      IF( INDX.EQ.1 ) THEN
        NGX = NFBN - NRFN - NXP + NWN_LW + NT_FRC - NXP_FRC + 
     &    NBN_BH - NXP_BH
!
!---  Geomechanics - number of unknowns equals number of active
!     finite element nodes x 3
!
      ELSEIF( INDX.EQ.2 ) THEN
        NGX = 3*NFEN_GM
!
!---  Number of unknowns equals number of number of coupled equations x
!     (active nodes + fracture triangles) + number of coupled wells  ---
!
      ELSEIF( N_CW.GT.0 ) THEN
        NGX = ISLV*(NFBN - NRFN - NXP + NWN_LW + NT_FRC - NXP_FRC) + 
     &    N_CW
!
!---  Number of unknowns equals number of number of coupled equations
!     x (number of active nodes + number of wells +
!     number of surface nodes)
!
      ELSEIF( LSPILL.GT.0 ) THEN
        NGX = ISLV*(NFBN - NRFN - NXP + NWLN + NWN_LW + IJFLD)
!
!---  Number of unknowns equals number of number of coupled equations
!     x (number of active nodes + number of wells + fracture triangles
!     + number of borehole nodes)
!
      ELSE
        NGX = ISLV*(NFBN - NRFN - NXP + NWLN + NWN_LW + NT_FRC - 
     &    NXP_FRC + NBN_BH - NXP_BH)
      ENDIF
!
!---  Dual-porosity option  ---
!
      IF( ISLC(11).EQ.1 ) NGX = NGX + ISLV*(NFBN-NRFN-NXP)
!
!---  Load Lis problem vector  ---
!
      CALL lis_vector_set_values(LIS_INS_VALUE,NGX,ILU,BLU,RHS_VEC,IERR)
      IF( ISLC(79).EQ.1 ) THEN
        INDX = VALUERSS()
        PRINT *,'Post Lis Vector Set Values: ',INDX,'IERR = ',IERR
      ENDIF
!
!---  Load nonzero elements of the Lis matrix for transport  ---
!
      IF( ISLV.EQ.-1 ) THEN
        NNZ = NLUC(NGX+1)
        CALL lis_matrix_set_csr(NNZ,NLUC,MLUC,DLU,MAT,IERR)
!
!---  Load nonzero elements of the Lis matrix for geomechanics  ---
!
      ELSEIF( INDX.EQ.2 ) THEN
        NNZ = NLU_GM(NGX+1)
        CALL lis_matrix_set_csr(NNZ,NLU_GM,MLU_GM,DLU,MAT,IERR)
!
!---  Load nonzero elements of the Lis matrix for coupled-flow  ---
!
      ELSE
        NNZ = NLU(NGX+1)
        CALL lis_matrix_set_csr(NNZ,NLU,MLU,DLU,MAT,IERR)
      ENDIF
      IF( ISLC(79).EQ.1 ) THEN
        INDX = VALUERSS()
        PRINT *,'Post Lis Matrix Set CSR: ',INDX,'IERR = ',IERR
      ENDIF
      CALL lis_matrix_assemble(MAT,IERR)
      IF( ISLC(79).EQ.1 ) THEN
        INDX = VALUERSS()
        PRINT *,'Post Lis Matrix Assemble: ',INDX,'IERR = ',IERR
      ENDIF
!
!---  Output coupled-flow problem vector  ---
!
      IF( ISLC(34).EQ.1 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER ) THEN
        CALL lis_output_vector(RHS_VEC,LIS_FMT_MM,"frhs_vec.dat",IERR)
        WRITE(ISC,'(/,A)') ' Output: frhs_vec.dat'
        STOP
      ENDIF
!
!---  Output coupled-flow matrix  ---
!
      IF( ISLC(34).EQ.2 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER ) THEN
        CALL lis_output_matrix(MAT,LIS_FMT_MM,"fmat.dat",IERR)
        WRITE(ISC,'(/,A)') ' Output: fmat.dat'
        STOP
      ENDIF
!
!---  Output transport problem vector  ---
!
      IF( ISLC(34).EQ.11 .AND. ISLV.EQ.-1 .AND. ISLC(83).LE.NSTEP ) THEN
        CALL lis_output_vector(RHS_VEC,LIS_FMT_MM,"crhs_vec.dat",IERR)
        WRITE(ISC,'(/,A)') ' Output: crhs_vec.dat'
        STOP
      ENDIF
!
!---  Output transport matrix  ---
!
      IF( ISLC(34).EQ.12 .AND. ISLV.EQ.-1 .AND. ISLC(83).LE.NSTEP ) THEN
        CALL lis_output_matrix(MAT,LIS_FMT_MM,"cmat.dat",IERR)
        WRITE(ISC,'(/,A)') ' Output: cmat.dat'
        STOP
      ENDIF
!
!---  Output geomechanics problem vector  ---
!
      IF( ISLC(34).EQ.21 .AND. INDX.EQ.2 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER .AND. ISLC(85).EQ.IREF_GM ) THEN
        CALL lis_output_vector(RHS_VEC,LIS_FMT_MM,"grhs_vec.dat",IERR)
        WRITE(ISC,'(/,A)') ' Output: grhs_vec.dat'
        STOP
      ENDIF
!
!---  Output geomechanics matrix  ---
!
      IF( ISLC(34).EQ.22 .AND. INDX.EQ.2 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER .AND. ISLC(85).EQ.IREF_GM ) THEN
        CALL lis_output_matrix(MAT,LIS_FMT_MM,"gmat.dat",IERR)
        WRITE(ISC,'(/,A)') ' Output: gmat.dat'
        STOP
      ENDIF
!
!---  Solve  ---
!
      CALL lis_solve(MAT,RHS_VEC,SOL_VEC,KSP,IERR)
      IF( IERR.NE.0 ) ICNV = 1
      IF( ISLC(79).EQ.1 ) THEN
        INDX = VALUERSS()
        PRINT *,'Post Lis Solve: ',INDX,'IERR = ',IERR
      ENDIF
!      CALL LIS_SOLVER_GET_STATUS( KSP,IVARX,IERR )
!      PRINT *,'Lis Solver Status = ',IVARX
!      CALL LIS_SOLVER_GET_ITER( KSP,ITER,IERR )
!      PRINT *,'Lis Solver Iterations = ',ITER
!      CALL LIS_SOLVER_GET_RESIDUALNORM( KSP,VARX,IERR )
!      PRINT *,'Lis Solver Residual Norm = ',VARX
!
!---  Output coupled-flow solution vector  ---
!
      IF( ISLC(34).EQ.3 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER ) THEN
        CALL lis_output_vector(SOL_VEC,LIS_FMT_MM,"fsol_vec.dat",IERR)
        WRITE(ISC,'(/,A)') ' Output: fsol_vec.dat'
        STOP
      ENDIF
!
!---  Output transport solution vector  ---
!
      IF( ISLC(34).EQ.13 .AND. ISLV.EQ.-1 .AND. ISLC(83).LE.NSTEP ) THEN
        CALL lis_output_vector(SOL_VEC,LIS_FMT_MM,"csol_vec.dat",IERR)
        WRITE(ISC,'(/,A)') ' Output: csol_vec.dat'
        STOP
      ENDIF
!
!---  Output geomechanics solution vector  ---
!
      IF( ISLC(34).EQ.23 .AND. INDX.EQ.2 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER .AND. ISLC(85).EQ.IREF_GM ) THEN
        CALL lis_output_vector(SOL_VEC,LIS_FMT_MM,"gsol_vec.dat",IERR)
        WRITE(ISC,'(/,A)') ' Output: gsol_vec.dat'
        STOP
      ENDIF
!
!---  Get solution  ---
!
      ISX = 1
      ICX = NGX
      CALL lis_vector_get_values(SOL_VEC,ISX,ICX,BLU,IERR)
      IF( ISLC(79).EQ.1 ) THEN
        INDX = VALUERSS()
        PRINT *,'Post Lis Get Values: ',INDX,'IERR = ',IERR
      ENDIF
      IF( IERR.NE.0 ) ICNV = 1
      DO M = 1,NGX
        IF( ABS(BLU(M))/EPSL.LT.EPSL ) BLU(M) = 0.D+0
      ENDDO
!!
!!---  Unset Jacobian matrix from DLU array  ---
!!
!      CALL lis_matrix_unset(MAT,IERR)
!      IF( ISLC(79).EQ.1 ) THEN
!        INDX = VALUERSS()
!        PRINT *,'Post Lis Matrix Unset: ',INDX,'IERR = ',IERR
!      ENDIF
!!
!!---  Destroy Jacobian matrix  ---
!!
!      CALL lis_matrix_destroy(MAT,IERR)
!      IF( ISLC(79).EQ.1 ) THEN
!        INDX = VALUERSS()
!        PRINT *,'Post Lis Matrix Destroy: ',INDX,'IERR = ',IERR
!      ENDIF
!!
!!---  Destroy solution vector  ---
!!
!      CALL lis_vector_destroy(SOL_VEC,IERR)
!      IF( ISLC(79).EQ.1 ) THEN
!        INDX = VALUERSS()
!        PRINT *,'Post Lis Solution Vector Destroy: ',INDX,'IERR = ',IERR
!      ENDIF
!!
!!---  Destroy problem vector  ---
!!
!      CALL lis_vector_destroy(RHS_VEC,IERR)
!      IF( ISLC(79).EQ.1 ) THEN
!        INDX = VALUERSS()
!        PRINT *,'Post Lis Problem Vector Destroy: ',INDX,'IERR = ',IERR
!      ENDIF
!!
!!---  Destroy solver  ---
!!
!      CALL lis_solver_destroy(KSP,IERR)
!      IF( ISLC(79).EQ.1 ) THEN
!        INDX = VALUERSS()
!        PRINT *,'Post Lis Solver Destroy: ',INDX,'IERR = ',IERR
!      ENDIF
!!
!!---  Create Jacobian matrix  ---
!!
!      CALL lis_matrix_create(LIS_COMM_WORLD,MAT,IERR)
!      IF( ISLC(79).EQ.1 ) THEN
!        INDX = VALUERSS()
!        PRINT *,'Post Lis Matrix Create: ',INDX,'IERR = ',IERR
!      ENDIF
!      CALL lis_matrix_set_size(MAT,NLX,NGX,IERR)
!      IF( ISLC(79).EQ.1 ) THEN
!        INDX = VALUERSS()
!        PRINT *,'Post Lis Matrix Set Size: ',INDX,'IERR = ',IERR
!      ENDIF
!!
!!---  Create solution vector  ---
!!
!      CALL lis_vector_duplicate(MAT,SOL_VEC,IERR)
!      IF( ISLC(79).EQ.1 ) THEN
!        INDX = VALUERSS()
!        PRINT *,'Post Lis Solution Vector Duplicate: ',INDX,
!     &    'IERR = ',IERR
!      ENDIF
!!
!!---  Create problem (RHS) vector  ---
!!
!      CALL lis_vector_duplicate(MAT,RHS_VEC,IERR)
!      IF( ISLC(79).EQ.1 ) THEN
!        INDX = VALUERSS()
!        PRINT *,'Post Lis Problem Vector Duplicate: ',INDX,
!     &    'IERR = ',IERR
!      ENDIF
!!
!!---  Create solver  ---
!!
!      CALL lis_solver_create(KSP,IERR)
!      IF( ISLC(79).EQ.1 ) THEN
!        INDX = VALUERSS()
!        PRINT *,'Post Lis Solver Create: ',INDX,'IERR = ',IERR
!      ENDIF
!!
!!---  Set the options for the solver  ---
!!
!      CALL lis_solver_set_option('-i bicgstab',KSP,IERR)
!      CALL lis_solver_set_option('-p ilu -ilu_fill 0',KSP,IERR)
!      CALL lis_solver_set_option('-print none -adds true',KSP,IERR)
!      CALL lis_solver_set_optionC(KSP,IERR)
!      IF( ISLC(79).EQ.1 ) THEN
!        INDX = VALUERSS()
!        PRINT *,'Post Lis Solver Set Options: ',INDX,'IERR = ',IERR
!      ENDIF
!
!---  End of STOMP_LIS_SOLVE group
!
      RETURN
      END SUBROUTINE STOMP_LIS_SOLVE
!
!---  End of module  ---
!
      END MODULE STOMP_LIS_MODULE

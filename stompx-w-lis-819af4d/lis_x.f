!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE STOMP_LIS_CREATE( KSP,MAT,RHS_VEC,SOL_VEC,NLX )
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
      USE MPI
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)











































	









































	



















































































































































































































































































































































































































      






























































































!
!----------------------Type Declarations-------------------------------!
!
      integer :: NLX,IERR
      integer*8 :: MAT
      integer*8 :: SOL_VEC
      integer*8 :: RHS_VEC
      integer*8 :: KSP
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/STOMP_LIS_CREATE'
!
!---  Create Jacobian matrix  ---
!
      CALL lis_matrix_create_f( MPI_COMM_WORLD,MAT,IERR )
      CALL lis_matrix_set_size_f( MAT,NLX,0,IERR )
      CALL lis_matrix_set_type_f( MAT,1,IERR )
!
!---  Create solution vector  ---
!
      CALL lis_vector_create_f( MPI_COMM_WORLD,SOL_VEC,IERR )
      CALL lis_vector_set_size_f( SOL_VEC,NLX,0,IERR )
!
!---  Create problem (RHS) vector  ---
!
      CALL lis_vector_create_f( MPI_COMM_WORLD,RHS_VEC,IERR )
      CALL lis_vector_set_size_f( RHS_VEC,NLX,0,IERR )
      CALL lis_vector_set_all_f( 0.D+0,RHS_VEC,IERR )
!
!---  Create solver  ---
!
      CALL lis_solver_create_f( KSP,IERR )
!
!---  Set the options for the solver  ---
!
      CALL lis_solver_set_option_f('-i bicgstab',KSP,IERR)
      CALL lis_solver_set_option_f('-p ilu -ilu_fill 0',KSP,IERR)
      CALL lis_solver_set_option_f('-print none -adds true',KSP,IERR)
      CALL lis_solver_set_optionc_f( KSP,IERR )
!
!---  End of STOMP_LIS_CREATE group
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

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
!     Destroy solution vector, problem vector and matrix for Lis solver
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

!
!----------------------Type Declarations-------------------------------!
!
      integer :: NLX,IERR
      integer*8 :: MAT
      integer*8 :: SOL_VEC
      integer*8 :: RHS_VEC
      integer*8 :: KSP
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/STOMP_LIS_DESTORY'
!
!---  Destroy Jacobian matrix  ---
!
      CALL lis_matrix_destroy_f( MAT,IERR )
!
!---  Destroy solution vector  ---
!
      CALL lis_vector_destroy_f( SOL_VEC,IERR )
!
!---  Destroy problem (RHS) vector  ---
!
      CALL lis_vector_destroy_f( RHS_VEC,IERR )
!
!---  Destroy solver  ---
!
      CALL lis_solver_destroy_f( KSP,IERR )
!
!---  End of STOMP_LIS_DESTROY group
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE STOMP_LIS_SOLVE( KSP,MAT,RHS_VEC,SOL_VEC,
     &  NUKOX,NUKLX,INDX )
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
!     Solve the linear system A x = b for Lis solver
!
!     INDX = 0 : Coupled flow and transport
!     INDX = 1 : Solute or species component transport
!     INDX = 2 : Geomechanics
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE SOLTN
      USE JACOB
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

!
!----------------------Type Declarations-------------------------------!
!
      integer :: IERR,ISX,ICX
      integer*8 :: MAT
      integer*8 :: SOL_VEC
      integer*8 :: RHS_VEC
      integer*8 :: KSP
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/STOMP_LIS_SOLVE'
!
!---  Assemble matrix  ---
!
      CALL lis_matrix_assemble_f( MAT,IERR )
      IF( IERR.NE.0 ) THEN
        I_ERR(1) = 0
        I_ERR(2) = 0
        I_ERR(3) = 0
        I_ERR(4) = ID
        IF( INDX.EQ.0 ) THEN
          M_ERR(1) = 'Coupled Flow and Transport Lis Matrix Assemble'
        ELSEIF( INDX.EQ.1 ) THEN
          M_ERR(1) = 'Solute Transport Lis Matrix Assemble'
        ELSEIF( INDX.EQ.2 ) THEN
          M_ERR(1) = 'Geomechanics Lis Matrix Assemble'
        ELSE
          M_ERR(1) = 'Lis Matrix Assemble'
        ENDIF
        RETURN
      ENDIF
!
!---  Output coupled-flow problem vector  ---
!
      IF( ISLC(34).EQ.1 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER .AND. INDX.EQ.0 ) THEN
        CALL lis_output_vector_f(RHS_VEC,2,"frhs_vecx.dat",IERR)
        IF( IERR.NE.0 ) THEN
          I_ERR(1) = 0
          I_ERR(2) = 0
          I_ERR(3) = 0
          I_ERR(4) = ID
          M_ERR(1) = 'Coupled-Flow Lis Output Vector'
          RETURN
        ENDIF
        ICNV = 4
        IF( ID.EQ.0 ) PRINT *,'Output: frhs_vecx.dat'
        RETURN
      ENDIF
!
!---  Output coupled-flow matrix  ---
!
      IF( ISLC(34).EQ.2 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER .AND. INDX.EQ.0 ) THEN
        CALL lis_output_matrix_f( MAT,2,'fmatx.dat',IERR )
        IF( IERR.NE.0 ) THEN
          I_ERR(1) = 0
          I_ERR(2) = 0
          I_ERR(3) = 0
          I_ERR(4) = ID
          M_ERR(1) = 'Coupled-Flow Lis Output Matrix'
          RETURN
        ENDIF
        ICNV = 4
        IF( ID.EQ.0 ) PRINT *,'Output: fmatx.dat'
        RETURN
      ENDIF
!
!---  Output transport problem vector  ---
!
      IF( ISLC(34).EQ.11 .AND. ISLC(83).LE.NSTEP .AND. INDX.EQ.1 ) THEN
        CALL lis_output_vector_f(RHS_VEC,2,"crhs_vecx.dat",IERR)
        IF( IERR.NE.0 ) THEN
          I_ERR(1) = 0
          I_ERR(2) = 0
          I_ERR(3) = 0
          I_ERR(4) = ID
          M_ERR(1) = 'Transport Lis Output Vector'
          RETURN
        ENDIF
        ICNV = 4
        IF( ID.EQ.0 ) PRINT *,'Output: crhs_vecx.dat'
        RETURN
      ENDIF
!
!---  Output transport matrix  ---
!
      IF( ISLC(34).EQ.12 .AND. ISLC(83).LE.NSTEP .AND. INDX.EQ.1 ) THEN
        CALL lis_output_matrix_f( MAT,2,'cmatx.dat',IERR )
        IF( IERR.NE.0 ) THEN
          I_ERR(1) = 0
          I_ERR(2) = 0
          I_ERR(3) = 0
          I_ERR(4) = ID
          M_ERR(1) = 'Transport Lis Output Matrix'
          RETURN
        ENDIF
        ICNV = 4
        IF( ID.EQ.0 ) PRINT *,'Output: cmatx.dat'
        RETURN
      ENDIF
!
!---  Output geomechanics problem vector  ---
!
      IF( ISLC(34).EQ.21 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER .AND. INDX.EQ.2 ) THEN
        CALL lis_output_vector_f(RHS_VEC,2,"grhs_vecx.dat",IERR)
        IF( IERR.NE.0 ) THEN
          I_ERR(1) = 0
          I_ERR(2) = 0
          I_ERR(3) = 0
          I_ERR(4) = ID
          M_ERR(1) = 'Geomechanics Lis Output Vector'
          RETURN
        ENDIF
        ICNV = 4
        IF( ID.EQ.0 ) PRINT *,'Output: grhs_vecx.dat'
        RETURN
      ENDIF
!
!---  Output geomechanics matrix  ---
!
      IF( ISLC(34).EQ.22 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER .AND. INDX.EQ.2 ) THEN
        CALL lis_output_matrix_f( MAT,2,'gmatx.dat',IERR )
        IF( IERR.NE.0 ) THEN
          I_ERR(1) = 0
          I_ERR(2) = 0
          I_ERR(3) = 0
          I_ERR(4) = ID
          M_ERR(1) = 'Geomechanics Lis Output Matrix'
          RETURN
        ENDIF
        ICNV = 4
        IF( ID.EQ.0 ) PRINT *,'Output: gmatx.dat'
        RETURN
      ENDIF
!
!---  Solve linear system  ---
!
      CALL lis_solve_f( MAT,RHS_VEC,SOL_VEC,KSP,IERR )
      IF( IERR.NE.0 ) THEN
        I_ERR(1) = 0
        I_ERR(2) = 0
        I_ERR(3) = 0
        I_ERR(4) = ID
        IF( INDX.EQ.0 ) THEN
          M_ERR(1) = 'Coupled Flow and Transport Lis Solve'
        ELSEIF( INDX.EQ.1 ) THEN
          M_ERR(1) = 'Solute Transport Lis Solve'
        ELSEIF( INDX.EQ.2 ) THEN
          M_ERR(1) = 'Geomechanics Lis Solve'
        ELSE
          M_ERR(1) = 'Lis Solve'
        ENDIF
        RETURN
      ENDIF
!     CALL lis_solver_get_iter_f( KSP,ITER,IERR )
!
!---  Output coupled-flow solution vector  ---
!
      IF( ISLC(34).EQ.3 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER .AND. INDX.EQ.0 ) THEN
        CALL lis_output_vector_f( SOL_VEC,2,'fsol_vecx.dat',
     &    IERR )
        IF( IERR.NE.0 ) THEN
          I_ERR(1) = 0
          I_ERR(2) = 0
          I_ERR(3) = 0
          I_ERR(4) = ID
          M_ERR(1) = 'Coupled-Flow Lis Output Vector'
          RETURN
        ENDIF
        ICNV = 4
        IF( ID.EQ.0 ) PRINT *,'Output: fsol_vecx.dat'
        RETURN
      ENDIF
!
!---  Output transport solution vector  ---
!
      IF( ISLC(34).EQ.13 .AND. ISLC(83).LE.NSTEP .AND. INDX.EQ.1 ) THEN
        CALL lis_output_vector_f( SOL_VEC,2,'csol_vecx.dat',
     &    IERR )
        IF( IERR.NE.0 ) THEN
          I_ERR(1) = 0
          I_ERR(2) = 0
          I_ERR(3) = 0
          I_ERR(4) = ID
          M_ERR(1) = 'Transport Lis Output Vector'
          RETURN
        ENDIF
        ICNV = 4
        IF( ID.EQ.0 ) PRINT *,'Output: csol_vecx.dat'
        RETURN
      ENDIF
!
!---  Output geomechanics solution vector  ---
!
      IF( ISLC(34).EQ.23 .AND. ISLC(83).LE.NSTEP .AND. 
     &  ISLC(84).LE.NITER .AND. INDX.EQ.2 ) THEN
        CALL lis_output_vector_f( SOL_VEC,2,'gsol_vecx.dat',
     &    IERR )
        IF( IERR.NE.0 ) THEN
          I_ERR(1) = 0
          I_ERR(2) = 0
          I_ERR(3) = 0
          I_ERR(4) = ID
          M_ERR(1) = 'Geomechanics Lis Output Vector'
          RETURN
        ENDIF
        ICNV = 4
        IF( ID.EQ.0 ) PRINT *,'Output: gsol_vecx.dat'
        RETURN
      ENDIF
!
!---  Get linear solution  ---
!
      ISX = NUKOX + 1
      ICX = NUKLX
      CALL lis_vector_get_values_f( SOL_VEC,ISX,ICX,BLU,IERR )
      IF( IERR.NE.0 ) THEN
        I_ERR(1) = 0
        I_ERR(2) = 0
        I_ERR(3) = 0
        I_ERR(4) = ID
        IF( INDX.EQ.0 ) THEN
          M_ERR(1) = 'Coupled Flow and Transport Lis Vector Get Values'
        ELSEIF( INDX.EQ.1 ) THEN
          M_ERR(1) = 'Solute Transport Lis Vector Get Values'
        ELSEIF( INDX.EQ.2 ) THEN
          M_ERR(1) = 'Geomechanics Lis Vector Get Values'
        ELSE
          M_ERR(1) = 'Lis Vector Get Values'
        ENDIF
        RETURN
      ENDIF
!
!---  Set underflows to zero  ---
!
      DO M = 1,NUKLX
        IF( ABS(BLU(M))/EPSL.LT.EPSL ) BLU(M) = 0.D+0
      ENDDO
!!
!!---  Initialize matrix for next solve  ---
!!
!      CALL lis_matrix_destroy_f( MAT,IERR )
!      IF( IERR.NE.0 ) THEN
!        I_ERR(1) = 0
!        I_ERR(2) = 0
!        I_ERR(3) = 0
!        I_ERR(4) = ID
!        IF( INDX.EQ.0 ) THEN
!          M_ERR(1) = 'Coupled Flow and Transport Lis Matrix Destroy'
!        ELSEIF( INDX.EQ.1 ) THEN
!          M_ERR(1) = 'Solute Transport Lis Matrix Destroy'
!        ELSEIF( INDX.EQ.2 ) THEN
!          M_ERR(1) = 'Geomechanics Lis Matrix Destroy'
!        ELSE
!          M_ERR(1) = 'Lis Matrix Destroy'
!        ENDIF
!        RETURN
!      ENDIF
!      CALL lis_matrix_create_f( MPI_COMM_WORLD,MAT,IERR )
!      IF( IERR.NE.0 ) THEN
!        I_ERR(1) = 0
!        I_ERR(2) = 0
!        I_ERR(3) = 0
!        I_ERR(4) = ID
!        IF( INDX.EQ.0 ) THEN
!          M_ERR(1) = 'Coupled Flow and Transport Lis Matrix Create'
!        ELSEIF( INDX.EQ.1 ) THEN
!          M_ERR(1) = 'Solute Transport Lis Matrix Create'
!        ELSEIF( INDX.EQ.2 ) THEN
!          M_ERR(1) = 'Geomechanics Lis Matrix Create'
!        ELSE
!          M_ERR(1) = 'Lis Matrix Create'
!        ENDIF
!        RETURN
!      ENDIF
!      CALL lis_matrix_set_size_f( MAT,NUKLX,0,IERR )
!      IF( IERR.NE.0 ) THEN
!        I_ERR(1) = 0
!        I_ERR(2) = 0
!        I_ERR(3) = 0
!        I_ERR(4) = ID
!        IF( INDX.EQ.0 ) THEN
!          M_ERR(1) = 'Coupled Flow and Transport Lis Matrix Set Size'
!        ELSEIF( INDX.EQ.1 ) THEN
!          M_ERR(1) = 'Solute Transport Lis Matrix Set Size'
!        ELSEIF( INDX.EQ.2 ) THEN
!          M_ERR(1) = 'Geomechanics Lis Matrix Set Size'
!        ELSE
!          M_ERR(1) = 'Lis Matrix Set Size'
!        ENDIF
!        RETURN
!      ENDIF
!      CALL lis_matrix_set_type_f( MAT,1,IERR )
!      IF( IERR.NE.0 ) THEN
!        I_ERR(1) = 0
!        I_ERR(2) = 0
!        I_ERR(3) = 0
!        I_ERR(4) = ID
!        IF( INDX.EQ.0 ) THEN
!          M_ERR(1) = 'Coupled Flow and Transport Lis Matrix Set Type'
!        ELSEIF( INDX.EQ.1 ) THEN
!          M_ERR(1) = 'Solute Transport Lis Matrix Set Type'
!        ELSEIF( INDX.EQ.2 ) THEN
!          M_ERR(1) = 'Geomechanics Lis Matrix Set Type'
!        ELSE
!          M_ERR(1) = 'Lis Matrix Set Type'
!        ENDIF
!        RETURN
!      ENDIF
!
!---  Initialize problem vector for next solve  ---
!
      CALL lis_vector_set_all_f( 0.D+0,RHS_VEC,IERR )
      IF( IERR.NE.0 ) THEN
        I_ERR(1) = 0
        I_ERR(2) = 0
        I_ERR(3) = 0
        I_ERR(4) = ID
        IF( INDX.EQ.0 ) THEN
          M_ERR(1) = 'Coupled Flow and Transport Lis Vector Set All'
        ELSEIF( INDX.EQ.1 ) THEN
          M_ERR(1) = 'Solute Transport Lis Vector Set All'
        ELSEIF( INDX.EQ.2 ) THEN
          M_ERR(1) = 'Geomechanics Lis Vector Set All'
        ELSE
          M_ERR(1) = 'Lis Vector Set All'
        ENDIF
        RETURN
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of STOMP_LIS_SOLVE group
!
      RETURN
      END



!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PUMFPK( ISLV,IRF )
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
!     Calling routine for the UMFPACK Version 2.0, Unsymmetric-pattern
!     Multifrontal Package for solving systems of sparse linear systems,
!     Ax=b; where A is sparse and can be unsymmetric.
!
!     Timothy A. Davis and Iain S. Duff, Copyright (C) 1995.
!     September 13, 1995.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1996.
!     Last Modified by MD White on August 1, 1996.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
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
      LOGICAL TRANSA
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PUMFPK'
#ifdef um
!
!---  Initialize control parameters for the UMFPACK routines  ---
!
      CALL UMD2IN( ICLU,CNLU,KPLU )
!
!---  Pivot from the diagonal of the Jacobian.  The nonzero pattern
!     of the Jacobian matrix is basically symmetric, therefore preserve
!     the pivots on the diagonal  ---
!
      ICLU(6) = 1
!
!---  Factorize Jacobian matrix ---
!
      N = (NFLD-NXP)*ISLV
      NE = NJME*ISLV*ISLV
      JOB = 0
      TRANSA = .FALSE.
      LVALUE = LJD
      LINDEX = LJI
      IF( IRF.EQ.1 ) THEN
        CALL UMD2FA( N,NE,JOB,TRANSA,LVALUE,LINDEX,ALU,ILU,KPLU,
     &    CNLU,ICLU,IFLU,RILU )
      ELSE
        CALL UMD2RF( N,NE,JOB,TRANSA,LVALUE,LINDEX,ALU,ILU,KPLU,
     &    CNLU,ICLU,IFLU,RILU )
      ENDIF
!
!---  Solve the linear system  ---
!
      CALL UMD2SO( N,JOB,TRANSA,LVALUE,LINDEX,ALU,ILU,KPLU,BLU,
     &  CLU,WLU,CNLU,ICLU,IFLU,RILU )
!
!---  Fill solution vector with current solution  ---
!
      DO 100 M = 1,N
        IF( ABS(CLU(M))/EPSL.GT.EPSL ) BLU(M) = CLU(M)
  100 CONTINUE
#endif
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PUMFPK group
!
      RETURN
      END


!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PSPLIB( ISLV,INDX )
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
!  Calling routine for the SPLIB: A Library of Iterative Methods for
!  Sparse Linear Systems.
!
!  Randall Bramley and Xiaoge Wang
!  Department of Computer Science
!  Indiana University - Bloomington
!  December 18, 1995
!
!  SOLMETH specifies the iterative solver following the schema:
!
!  SOLMETH = 1 --> Bi-Conjugate Gradient
!  SOLMETH = 2 --> Conjugate Gradient for system AA'y = b, x = A'y
!  SOLMETH = 3 --> Conjugate Gradient for system A'Ax = A'b
!  SOLMETH = 4 --> Conjugate Gradient Squared
!  SOLMETH = 5 --> Conjugate Gradients Stabilized
!  SOLMETH = 6 --> GMRES(k)
!  SOLMETH = 7 --> Transpose-free QMR
!  SOLMETH = 8 --> Templates version of Conjugate Gradient Stablized
!  SOLMETH = 9 --> Templates version of GMRES
!  SOLMETH = 10 --> Jacobi Iteration
!  SOLMETH = 11 --> Gauss-Siedel
!  SOLMETH = 12 --> Successive Over-relaxation (SOR)
!  SOLMETH = 13 --> Orthomin(m)
!
!  PCMETH specifies the preconditioner following the schema:
!
!  PCMETH = 0 --> no preconditioner
!  PCMETH = 1 --> ILU(s)
!  PCMETH = 2 --> MILU(s,rpcprm)
!  PCMETH = 3 --> ILUT(s,rpcprm)
!  PCMETH = 4 --> SSOR(rpcprm)
!  PCMETH = 5 --> TRID(s), where s is the block size
!  PCMETH = 6 --> ILU0, the space saver version of ILU(O)
!  PCMETH = 7 --> ECIMGS(rpcprm)
!
!  ERRFLG  : Integer variable indicating error conditions.  For the
!            solvers, the interpretation is:
!
!  ERRFLG = -m --> Not enough workspace provided; m more double
!                  precision words needed in array work().
!  ERRFLG = 0  --> Successful return; tolerance achieved.
!  ERRFLG = 1  --> The initial guess satisfies the stopping test.
!                  This is often a clue that rhs not set properly.
!  ERRFLG = 2  --> Tolerance reset 8 times in attempt to
!                  get true residual small enough.
!  ERRFLG = 3  --> Convergence not achieved in maxits iterations.
!  ERRFLG = 4  --> Preconditioned residual is too large relative
!                  to the unpreconditioned residual.  Indicates
!                  an unstable preconditioner.
!  4 < ERRFLG < 9 --> Breakdown condition, usually from an attempt to
!                     divide by a near-zero number in algorithm.
!  ERRFLG = 10+k  --> Error condition k occured in preconditioner.
!
!  Note that any of conditions 0, 1, 2 may mean success.
!
!  IREDO   : Recompute preconditioner if redo = 0, reuse if redo = 1.
!            This is useful if you wish to reuse a preconditioner from
!            an earlier call to SPLIB.  Converted internally to logical
!            variable redo
!
!  IUNIT   : vector of I/O unit numbers to use for output.
!            The I/O unit numbers are provided in the vector iunit(16).
!            The units are optional, depending on the error norm(s) 
!            you wish toplot/use.  However, if you specify instlvl >= 0,
!            the following units will be opened for writing:
!
!  IUNIT(1):    File for output of summary data
!  IUNIT(2-3):  Not used here
!  IUNIT(4):    || Ax - b ||  vs. iter
!  IUNIT(5):    || inv(M) (Ax-b) || vs. iter
!  IUNIT(6):    || Ax - b || / || r_0 || vs. iter
!  IUNIT(7):    || inv(M) (Ax-b) || / || inv(M) r_0 || vs. iter
!  IUNIT(8):    || Ax - b ||_inf / ( ||A||_inf ||x||_1 + ||b||_inf) vs. iter
!  IUNIT(9):    || x - xstar || vs. iter
!  IUNIT(10):   Not used.
!  IUNIT(11):   || Ax - b || vs. time
!  IUNIT(12):   || inv(M) (Ax-b) || vs. time
!  IUNIT(13):   || Ax - b || / || r_0 || vs. time
!  IUNIT(14):   || inv(M) (Ax-b) || / || inv(M) r_0 || vs. time
!  IUNIT(15):   || Ax - b ||_inf / ( ||A||_inf ||x||_1 + ||b||_inf) vs. time
!  IUNIT(16):   || x - xstar || vs. time
!
!  WORK    : Array of work storage, for preconditioners, solvers, and.
!            instrumentation routines.  The matrices and vectors used
!            are all indexed from this integer array, using pointers.
!
!  SPACE   : maximum number of integer words allowed for all problem
!            dependent arrays.  This includes the preconditioners
!            and the additional work vectors needed in the solvers.
!            On return, space holds the actual amount of storage
!            required by SPLIB.
!
!  INSTLVL : Level of instrumentation.  The level is cummulative, so specifying
!            output of Oetli/Prager norms implies residual norms, etc.
!
!  INSTLVL < 0  ---> no output from splib
!  INSTLVL = 0  ---> only summary data
!  INSTLVL = 1  ---> residual norm data
!  INSTLVL = 2  ---> preconditioned residual norm
!  INSTLVL = 3  ---> relative residual norm
!  INSTLVL = 4  ---> relative preconditioned residual norm
!  INSTLVL = 5  ---> Oetli/Prager norms
!  INSTLVL = 6  ---> error norms (if true solution vector is provided
!                    via a common block)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, November, 1998.
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
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER(LSPACE=600)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER N
#ifdef sp
      INTEGER IPARMS(2),ERRFLG,IUNIT(16),INSTLVL
      INTEGER PCMETH,SOLMETH,SPACE,RELTV,MAXITS,SCALE,IREDO
      REAL*8 TOL,RPARMS(2)
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE ::  DL,DR
      INTEGER, DIMENSION(:), ALLOCATABLE ::  WORK
      CHARACTER*64 FDUM
!
!----------------------Data Statements---------------------------------!
!
      DATA IUNIT / 41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56 /
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PSPLIB'
      ERRFLG = 0
!
!---  Number of unknowns equals number of active nodes + 
!     number of fracture triangles + number of borehole nodes for
!     transport equations ---
!
      IF( INDX.EQ.1 ) THEN
        N = NFBN - NRFN - NXP + NWN_LW + NT_FRC - NXP_FRC + NBN_BH - 
     &    NXP_BH
!
!---  Geomechanics - number of unknowns equals number of active
!     finite element nodes x 3
!
      ELSEIF( INDX.EQ.2 ) THEN
        N = 3*NFEN_GM
!
!---  Number of unknowns equals number of number of coupled equations x
!     (active nodes + fracture triangles) + number of coupled wells  ---
!
      ELSEIF( N_CW.GT.0 ) THEN
        N = ISLV*(NFBN - NRFN - NXP + NWN_LW + NT_FRC - NXP_FRC) + N_CW
!
!---  Number of unknowns equals number of number of coupled equations
!     x (number of active nodes + number of wells +
!     number of surface nodes)
!
      ELSEIF( LSPILL.GT.0 ) THEN
        N = ISLV*(NFBN - NRFN - NXP + NWLN + NWN_LW + IJFLD)
!
!---  Number of unknowns equals number of number of coupled equations
!     x (number of active nodes + number of wells + fracture triangles
!     + number of borehole nodes)
!
      ELSE
        N = ISLV*(NFBN - NRFN - NXP + NWLN + NWN_LW + NT_FRC - NXP_FRC
     &    + NBN_BH - NXP_BH)
      ENDIF
!
!---  Dual-porosity option  ---
!
      IF( ISLC(11).EQ.1 .AND. INDX.EQ.1 ) 
     &  N = N + ISLV*(NFBN-NRFN-NXP+NWN_LW)
#ifdef sp
      IF( INDX.EQ.-1 ) THEN
        IF( MOD(ISLC(14),10).EQ.1 ) THEN
          OPEN(UNIT=IUNIT(1), FILE='splib_out', STATUS='UNKNOWN',
     &    FORM='FORMATTED')
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      ALLOCATE( DL(1:LJF),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: DL'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( DR(1:LJF),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: DR'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( WORK(1:LSPACE*LJF),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: WORK'
        CALL WRMSGS( INDX )
      ENDIF
      TOL = 1.D-6
      RELTV = 1
      MAXITS = 200
      PCMETH = 1
      SOLMETH = 5
      SPACE = LSPACE*LJF
      SCALE = 0
      IF( MOD(ISLC(14),100)/10.GT.0 ) SCALE = 1
      IPARMS(1) = 1
      IPARMS(2) = 0
      RPARMS(1) = 0.D+0
      RPARMS(2) = 0.D+0
      IREDO = 0
      INSTLVL = -1
      IF( MOD(ISLC(14),10).EQ.1 ) INSTLVL = 0
!
!---  Initialize solution vector  ---
!
      DO 100 M = 1,N
        CLU(M) = 0.D+0
  100 CONTINUE
!
!---  Call SPLIB for solute transport,
!     skip solution for zero solute concentration throughout
!     active domain ---
!
      IF( INDX.EQ.1 ) THEN
        DO 110 M = 1,N
          IF( ABS(BLU(M))/EPSL.GT.EPSL ) GOTO 112
          BLU(M) = 0.D+0
  110   CONTINUE
        GOTO 114
  112   CONTINUE
!
!---    Output transport problem vector  ---
!
        IF( ISLC(34).EQ.11 ) THEN
          FDUM = 'rhs_vec_c.dat'
          CALL SPLIB_OUTPUT_VECTOR( N,FDUM )
          STOP
        ENDIF
!
!---    Output transport matrix  ---
!
        IF( ISLC(34).EQ.12 ) THEN
          FDUM = 'mat_c.dat'
          CALL SPLIB_OUTPUT_MATRIX( MLUC,NLUC,N,FDUM )
          STOP
        ENDIF
        IF( N.EQ.1 ) THEN
          CLU(1) = BLU(1)/DLU(1)
        ELSE
          CALL SPLIB( DLU,MLUC,NLUC,N,CLU,BLU,WORK,SPACE,TOL,RELTV,
     &      MAXITS,PCMETH,SOLMETH,SCALE,DL,DR,IPARMS,RPARMS,IREDO,
     &      ERRFLG,IUNIT,INSTLVL )
        ENDIF
  114   CONTINUE
!
!---  Call SPLIB for geomechanics ---
!
      ELSEIF( INDX.EQ.2 ) THEN
!
!---    Output geomechanics problem vector  ---
!
        IF( ISLC(34).EQ.21 ) THEN
          FDUM = 'rhs_vec_g.dat'
          CALL SPLIB_OUTPUT_VECTOR( N,FDUM )
          STOP
        ENDIF
!
!---    Output geomechanics matrix  ---
!
        IF( ISLC(34).EQ.22 ) THEN
          FDUM = 'mat_g.dat'
          CALL SPLIB_OUTPUT_MATRIX( MLU_GM,NLU_GM,N,FDUM )
          STOP
        ENDIF
        CALL SPLIB( DLU,MLU_GM,NLU_GM,N,CLU,BLU,WORK,SPACE,TOL,RELTV,
     &    MAXITS,PCMETH,SOLMETH,SCALE,DL,DR,IPARMS,RPARMS,IREDO,
     &    ERRFLG,IUNIT,INSTLVL )
!
!---  Call SPLIB for coupled flow ---
!
      ELSE
!
!---    Output coupled-flow problem vector  ---
!
        IF( ISLC(34).EQ.1 ) THEN
          FDUM = 'rhs_vec.dat'
          CALL SPLIB_OUTPUT_VECTOR( N,FDUM )
          STOP
        ENDIF
!
!---    Output coupled-flow matrix  ---
!
        IF( ISLC(34).EQ.2 ) THEN
          FDUM = 'mat.dat'
          CALL SPLIB_OUTPUT_MATRIX( MLU,NLU,N,FDUM )
          STOP
        ENDIF
        CALL SPLIB( DLU,MLU,NLU,N,CLU,BLU,WORK,SPACE,TOL,RELTV,
     &    MAXITS,PCMETH,SOLMETH,SCALE,DL,DR,IPARMS,RPARMS,IREDO,
     &    ERRFLG,IUNIT,INSTLVL )
      ENDIF
!
!---  SPLIB error flags  ---
!
!      IF( ERRFLG.EQ.4 ) THEN
!        INDX = 26
!        CHMSG = 'SPLIB Error: Unstable Preconditioner: '
!        CALL WRMSGS( INDX )
!      ELSEIF( ERRFLG.GT.1 .AND. N.GT.1 ) THEN
      IF( ERRFLG.GT.1 .AND. N.GT.1 ) THEN
        INDX = 20
        IMSG = ERRFLG
        IF( ERRFLG.EQ.1 ) THEN
          CHMSG = 'SPLIB: Initial Guess Statisfies Stopping Test: '
        ELSEIF( ERRFLG.EQ.2 ) THEN
          CHMSG = 'SPLIB: Tolerance Reset 8 Times: '
        ELSEIF( ERRFLG.EQ.3 ) THEN
          CHMSG = 'SPLIB: No Convergence after Maximum Iterations: '
        ELSEIF( ERRFLG.GT.4 ) THEN
          CHMSG = 'SPLIB: Unstable Preconditioner: '
        ELSEIF( ERRFLG.GT.5 .AND. ERRFLG.LT.9 ) THEN
          CHMSG = 'SPLIB: Divide by Near-Zero Number: '
        ELSE
          CHMSG = 'SPLIB Failure: (Non-Zero Error Flag): '
        ENDIF
        CALL WRMSGS( INDX )
        DO 120 M = 1,N
          BLU(M) = 0.D+0
  120   CONTINUE
        ICNV = 1
        IF( ALLOCATED(DL) ) THEN
        DEALLOCATE( DL,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Deallocation Error: DL'
          CALL WRMSGS( INDX )
        ENDIF
        ENDIF
        IF( ALLOCATED(DR) ) THEN
        DEALLOCATE( DR,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Deallocation Error: DR'
          CALL WRMSGS( INDX )
        ENDIF
        ENDIF
        IF( ALLOCATED(WORK) ) THEN
        DEALLOCATE( WORK,STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Deallocation Error: WORK'
          CALL WRMSGS( INDX )
        ENDIF
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Fill Solution Vector with Current Solution
!
      DO 200 M = 1,N
        IF( ABS(CLU(M))/EPSL.GT.EPSL ) BLU(M) = CLU(M)
  200 CONTINUE
!
!---  Output transport solution vector  ---
!
      IF( ISLC(34).EQ.13 .AND. INDX.EQ.1 ) THEN
        FDUM = 'fsol_vec_c.dat'
        CALL SPLIB_OUTPUT_VECTOR( N,FDUM )
        STOP
      ENDIF
!
!---  Output geomechanics solution vector  ---
!
      IF( ISLC(34).EQ.23 .AND. INDX.EQ.2 ) THEN
        FDUM = 'fsol_vec_g.dat'
        CALL SPLIB_OUTPUT_VECTOR( N,FDUM )
        STOP
      ENDIF
!
!---  Output fracture flow solution vector  ---
!
      IF( ISLC(34).EQ.33 .AND. INDX.EQ.3 ) THEN
        FDUM = 'fsol_vec_frc.dat'
        CALL SPLIB_OUTPUT_VECTOR( N,FDUM )
        STOP
      ENDIF
!
!---  Output coupled-flow solution vector  ---
!
      IF( ISLC(34).EQ.3 ) THEN
        FDUM = 'fsol_vec.dat'
        CALL SPLIB_OUTPUT_VECTOR( N,FDUM )
        STOP
      ENDIF
      IF( ALLOCATED(DL) ) THEN
      DEALLOCATE( DL,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: DL'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
      IF( ALLOCATED(DR) ) THEN
      DEALLOCATE( DR,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: DR'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
      IF( ALLOCATED(WORK) ) THEN
      DEALLOCATE( WORK,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: WORK'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
#endif
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PSPLIB group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SPLIB_OUTPUT_MATRIX( MLUX,NLUX,NUK,FDUM )
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
!    Output the Jacobian or transport matrix
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 03 May 2016
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER MLUX(*),NLUX(*)
      CHARACTER*64 ADUM,FDUM
      CHARACTER*16 FORM1
      CHARACTER*22 FORM2
      CHARACTER*21 FORM3
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM1 /'(I1,1X,I1,1X,I1)'/
      DATA FORM2 /'(I1,1X,I1,2X,1PE27.20)'/
      DATA FORM3 /'(I1,1X,I1,2X,1PE12.5)'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SPLIB_OUTPUT_MATRIX'
      NCH = LEN_TRIM(FDUM)
      OPEN(UNIT=26,FILE=FDUM(1:NCH),STATUS='UNKNOWN',FORM='FORMATTED')
      CLOSE(UNIT=26,STATUS='DELETE')
      OPEN(UNIT=26, FILE=FDUM(1:NCH), STATUS='NEW', FORM='FORMATTED')
      WRITE(26,'(A)') '%%MatrixMarket matrix coordinate real general'
      ICX = ICOUNT(NUK)
      WRITE(FORM1(3:3),'(I1)') ICX
      WRITE(FORM1(9:9),'(I1)') ICX
      NNZ = NLUX(NUK+1)-1
      ICX = ICOUNT(NNZ)
      WRITE(FORM1(15:15),'(I1)') ICX
      WRITE(26,FORM1) NUK,NUK,NNZ
      DO N = 1,NUK
        DO M = NLUX(N),(NLUX(N+1)-1)
          ICX = ICOUNT(N)
          WRITE(FORM3(3:3),'(I1)') ICX
          ICX = ICOUNT(MLUX(M))
          WRITE(FORM3(9:9),'(I1)') ICX
          WRITE(ADUM,FORM3) N,MLUX(M),DLU(M)
          CALL LCASE( ADUM )
          NCH = LEN_TRIM(ADUM)
          WRITE(26,'(A)') ADUM(1:NCH)
        ENDDO
      ENDDO
      CLOSE(UNIT=26)
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SPLIB_OUTPUT_MATRIX group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SPLIB_OUTPUT_VECTOR( NUK,FDUM )
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
!    Output the problem or solution vector
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 03 May 2016
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM,FDUM
      CHARACTER*6 FORM1
      CHARACTER*16 FORM2
      CHARACTER*15 FORM3
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM1 /'(I1)'/
      DATA FORM2 /'(I1,2X,1pe27.20)'/
      DATA FORM3 /'(I1,2X,1pe12.5)'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SPLIB_OUTPUT_VECTOR'
      NCH = LEN_TRIM(FDUM)
      OPEN(UNIT=26,FILE=FDUM(1:NCH),STATUS='UNKNOWN',FORM='FORMATTED')
      CLOSE(UNIT=26,STATUS='DELETE')
      OPEN(UNIT=26, FILE=FDUM(1:NCH), STATUS='NEW', FORM='FORMATTED')
      WRITE(26,'(A)') '%%MatrixMarket vector coordinate real general'
      ICX = ICOUNT(NUK)
      WRITE(FORM1(3:3),'(I1)') ICX
      WRITE(26,FORM1) NUK
      DO N = 1,NUK
        ICX = ICOUNT(N)
        WRITE(FORM3(3:3),'(I1)') ICX
        WRITE(ADUM,FORM3) N,BLU(N)
        CALL LCASE( ADUM )
        NCH = LEN_TRIM(ADUM)
        WRITE(26,'(A)') ADUM(1:NCH)
      ENDDO
      CLOSE(UNIT=26)
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SPLIB_OUTPUT_VECTOR group
!
      RETURN
      END


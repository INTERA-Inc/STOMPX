!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PNSPCG( ISLV )
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
!     Calling routine for the NSPCG (Nonsymmetric Preconditioned
!     Conjugate Gradient) solver for large sparse linear systems by
!     various iterative methods.
!
!     NSPCG User's Guide, Version 1.0
!     T. C. Oppe, W. D. Joubert, D. R. Kincaid
!     Center for Numerical Analysis
!     The University of Texas at Austin
!     April, 1988.
!
!     Parameter Arrays IPARM and RPARM
!
!     IPARM(1)  : (NTEST) Stopping Test Number [2]
!     IPARM(2)  : (ITMAX) Maximum Number of Iterations [100]
!     IPARM(3)  : (LEVEL) Output Level [0]
!     IPARM(4)  : (NOUT) Fortran Unit Number for Output [6]
!     IPARM(5)  : (IDGTS) Error Analysis Switch [0]
!     IPARM(6)  : (MAXADP) Adaptive Procedure Switch for EMAX [1]
!     IPARM(7)  : (MINADP) Adaptive Procedure Switch for EMIN [1]
!     IPARM(8)  : (IOMAD) Adaptive Procedure Switch for OMEGA [1]
!     IPARM(9)  : (NS1) Number of Old Vectors Saved [5]
!     IPARM(10) : (NS2) Frequency of Restarting (Large Number) [100000]
!     IPARM(11) : (NS3) Size of the Largest Hessenberg Matrix [0]
!     IPARM(12) : (NSTORE) Storage Mode [1]
!     IPARM(13) : (ISCALE) Matrix Scaling Switch [0]
!     IPARM(14) : (IPERM) Matrix Permuting Switch [0]
!     IPARM(15) : (IFACT) Factorization Switch [1]
!     IPARM(16) : (LVFILL) Level of Point or Block Fill-In [0]
!     IPARM(17) : (LTRUNC) Truncation Bandwidth [0]
!     IPARM(18) : (IPROPA) Matrix Property "A" Switch [2]
!     IPARM(19) : (KBLSZ) 1-D Block Size [-1]
!     IPARM(20) : (NBL2D) 2-D Block Size [-1]
!     IPARM(21) : (IFCTV) Scalar/Vectorization Switch [1]
!     IPARM(22) : (IQLR) Orientation of the Basic Preconditioner [1]
!     IPARM(23) : (ISYMM) Matrix Symmetry Switch [2]
!     IPARM(24) : (IELIM) Row Elimination Switch [0]
!     IPARM(25) : (NDEG) Degree of Polynomial Preconditioner [1]
!     RPARM(1)  : (ZETA) Stopping Test Value [1.E-6]
!     RPARM(2)  : (EMAX) Eigenvalue Estimate of Precond. Matrix [2.]
!     RPARM(3)  : (EMIN) Eigenvalue Estimate of Precond. Matrix [1.]
!     RPARM(4)  : (FF) Adapative Procedure Damping Factor [0.75]
!     RPARM(5)  : (FFF) Adapative Procedure Damping Factor [0.75]
!     RPARM(6)  : (TIMIT) Iteration Time in Seconds [0.]
!     RPARM(7)  : (DIGIT1) Digits Accuracy of Solution [0.]
!     RPARM(8)  : (DIGIT2) Digits Accuracy of Solution [0.]
!     RPARM(9)  : (OMEGA) Overrelaxation Parameter [1.]
!     RPARM(10) : (ALPHAB) Estimate of the Minimum Eigenvalue [0.]
!     RPARM(11) : (BETAB) Estimate of the Maximum Eigenvalue [0.]
!     RPARM(12) : (SPECR) Estimate of the Spectral Radius [0.]
!     RPARM(13) : (TIMFAC) Factorization Time in Seconds [0.]
!     RPARM(14) : (TIMTOT) Total Time for NSPCG in Seconds [0.]
!     RPARM(15) : (TOL) Tolerance for Eliminating Eqns. [500*SRELPR]
!     RPARM(16) : (AINF) Infinity Norm of the Matrix "A" [0.]
!
!     Error Messages Returned by Index IER
!
!     0  : No error
!    -1  : Nonpositive matrix size
!    -2  : Insufficient real work space
!    -3  : Insufficient integer work space
!    -4  : Nonpositive diagonal element
!    -5  : Nonexistent diagonal element
!    -6  : A is not positive definite
!    -7  : Q is not positive definite
!    -8  : Cannot permute matrix as requested
!    -9  : MDIM is not large enough to allow expansion of matrix
!   -10  : Inadmissible parameter encountered
!   -11  : Incorrect storage method for block method
!   -12  : Zero pivot encountered in factorization
!   -13  : Breakdown in direction vector calculation
!   -14  : Breakdown in attempt to perform rotation
!   -15  : Breakdown in iterate calculation
!   -16  : Unimplemented combination of parameters
!   -17  : Unable to perform eigenvalue stimation
!     1  : Failure to converg in itmax iterations
!     2  : ZETA too small - reset to 500*SPELPR
!     3  : ZBRENT failed to converg in MAXFN iterations
!          (signifies difficulty in eigenvalue estimation)
!     4  : In ZBRENT, F(A) and F(B) have the same sign
!          (signifies difficulty in eigenvalue estimation)
!     5  : Negative pivot encountered in factorization
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1994.
!     Last Modified by MD White, Battelle, PNL, April 14, 1994.
!     pnspcg.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE FILES
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!---  Generalized Minimal Residual Method (GMRES) Accelerator  ---
!
      PARAMETER(LNV=100*LCG,LNH=102*LCG)
      PARAMETER(LNWA=LCG*((LNH*(LNV+3))+(LJA*(LNV+3))+(2*(LNV+2)**2)+
     &  (7*LNV)+17+(LNH**2)+(2*LNH)+LJA*(LNV+3)+LNV+1))
      PARAMETER(LNWS=3*LJA)
!
!---  Incomplete Cholesky (IC) Preconditioner  ---
!
      PARAMETER(LNWP=LJA*LCG,LINWP=0)
      PARAMETER(LNWF=4*LJA*LJE*LCG,LINWF=4*LJA*LJE*LCG)
!
!---  Conjugate Gradient Parameters  ---
!
      PARAMETER(LNW=LNWA+LNWS+LNWP+LNWF+LBD)
      PARAMETER(LINW=LINWP+LINWF+LBD)
#ifdef cg
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER IP(1),JP(1),IWKSP(LINW),IPARM(25)
      REAL*8 DLU(1),WKSP(LNW),RPARM(16)
      EXTERNAL IC1,GMRES
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PNSPCG'
      IER = 0
      N = ISLV*(NFLD-NXP)
#ifdef cg
      NROW = LJD
      NCOL = LJE
      MAXNZ = ISLV*(2*NDIM+1)
      NW = LNW
      INW = LINW
      IF( N .GT. NROW ) THEN
        WRITE(ISC,'(2X,A)') 'ERROR: Jacobian Matrix Row Index'
        WRITE(IWR,'(2X,A)') 'ERROR: Jacobian Matrix Row Index'
        STOP
      ELSEIF( MAXNZ .GT. NCOL ) THEN
        WRITE(ISC,'(2X,A)') 'ERROR: Jacobian Matrix Column Index'
        WRITE(IWR,'(2X,A)') 'ERROR: Jacobian Matrix Column Index'
        STOP
      ENDIF
!
!---  Assign default values for the IPARM and RPARM arrays  ---
!
      CALL DFAULT ( IPARM,RPARM )
!
!---  Customized values for the IPARM and RPARM arrays  ---
!
      IPARM(4) = ISC
      IPARM(12) = 1
!
!---  Call NSPCG Routines  ---
!
      CALL NSPCG( IC1,GMRES,NROW,NCOL,N,MAXNZ,ALU,KLU,IP,JP,CLU,DLU,
     &  BLU,WKSP,IWKSP,NW,INW,IPARM,RPARM,IER )
!
!---  Fill Solution Vector with Current Solution
!
      DO 100 M = 1,N
        IF( ABS(CLU(M))/EPSL.GT.EPSL ) BLU(M) = CLU(M)
  100 CONTINUE
!
!---  Handle NSPCG Error Messages  ---
!
      IF( IER .EQ. -1 ) THEN
        WRITE(ISC,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(IWR,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(ISC,'(2X,A)') 'Nonpositive Matrix Size'
        WRITE(IWR,'(2X,A)') 'Nonpositive Matrix Size'
        STOP
      ELSEIF( IER .EQ. -2 ) THEN
        WRITE(ISC,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(IWR,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(ISC,'(2X,A)') 'Insufficient Real Work Space'
        WRITE(IWR,'(2X,A)') 'Insufficient Real Work Space'
        WRITE(ISC,'(2X,A,I4)') 'NW = ',NW
        WRITE(IWR,'(2X,A,I4)') 'NW = ',NW
        STOP
      ELSEIF( IER .EQ. -3 ) THEN
        WRITE(ISC,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(IWR,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(ISC,'(2X,A)') 'Insufficient Integer Work Space'
        WRITE(IWR,'(2X,A)') 'Insufficient Integer Work Space'
        WRITE(ISC,'(2X,A,I4)') 'INW = ',INW
        WRITE(IWR,'(2X,A,I4)') 'INW = ',INW
        STOP
      ELSEIF( IER .EQ. -4 ) THEN
        WRITE(ISC,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(IWR,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(ISC,'(2X,A)') 'Nonpositive Diagonal Element'
        WRITE(IWR,'(2X,A)') 'Nonpositive Diagonal Element'
        STOP
      ELSEIF( IER .EQ. -5 ) THEN
        WRITE(ISC,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(IWR,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(ISC,'(2X,A)') 'Nonexisent Diagonal Element'
        WRITE(IWR,'(2X,A)') 'Nonexisent Diagonal Element'
        STOP
      ELSEIF( IER .EQ. -6 ) THEN
        WRITE(ISC,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(IWR,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(ISC,'(2X,A)') 'A is Not Positive Definite'
        WRITE(IWR,'(2X,A)') 'A is Not Positive Definite'
        STOP
      ELSEIF( IER .EQ. -7 ) THEN
        WRITE(ISC,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(IWR,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(ISC,'(2X,A)') 'Q is Not Positive Definite'
        WRITE(IWR,'(2X,A)') 'Q is Not Positive Definite'
        STOP
      ELSEIF( IER .EQ. -8 ) THEN
        WRITE(ISC,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(IWR,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(ISC,'(2X,A)') 'Cannot Permute Matrix as Requested'
        WRITE(IWR,'(2X,A)') 'Cannot Permute Matrix as Requested'
        STOP
      ELSEIF( IER .EQ. -9 ) THEN
        WRITE(ISC,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(IWR,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(ISC,'(2X,A)') 'MDIM is Not Large Enough to Allow Expansion
     & of Matrix'
        WRITE(IWR,'(2X,A)') 'MDIM is Not Large Enough to Allow Expansion
     & of Matrix'
        STOP
      ELSEIF( IER .EQ. -10 ) THEN
        WRITE(ISC,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(IWR,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(ISC,'(2X,A)') 'Inadmissible Parameter Encountered'
        WRITE(IWR,'(2X,A)') 'Inadmissible Parameter Encountered'
        STOP
      ELSEIF( IER .EQ. -11 ) THEN
        WRITE(ISC,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(IWR,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(ISC,'(2X,A)') 'Incorrect Storage Method for Block Method'
        WRITE(IWR,'(2X,A)') 'Incorrect Storage Method for Block Method'
        STOP
      ELSEIF( IER .EQ. -12 ) THEN
        WRITE(ISC,'(2X,A)') 'WARNING: NSPCG Convergence Failure'
        WRITE(IWR,'(2X,A)') 'WARNING: NSPCG Convergence Failure'
        WRITE(ISC,'(2X,A)') 'Zero Pivot Encountered in Factorization'
        WRITE(IWR,'(2X,A)') 'Zero Pivot Encountered in Factorization'
        ICNV = 1
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ELSEIF( IER .EQ. -13 ) THEN
        WRITE(ISC,'(2X,A)') 'WARNING: NSPCG Convergence Failure'
        WRITE(IWR,'(2X,A)') 'WARNING: NSPCG Convergence Failure'
        WRITE(ISC,'(2X,A)') 'Breakdown in Direction Vector Calculation'
        WRITE(IWR,'(2X,A)') 'Breakdown in Direction Vector Calculation'
        ICNV = 1
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ELSEIF( IER .EQ. -14 ) THEN
        WRITE(ISC,'(2X,A)') 'WARNING: NSPCG Convergence Failure'
        WRITE(IWR,'(2X,A)') 'WARNING: NSPCG Convergence Failure'
        WRITE(ISC,'(2X,A)') 'Breakdown in Attempt to Perform Rotation'
        WRITE(IWR,'(2X,A)') 'Breakdown in Attempt to Perform Rotation'
        ICNV = 1
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ELSEIF( IER .EQ. -15 ) THEN
        WRITE(ISC,'(2X,A)') 'WARNING: NSPCG Convergence Failure'
        WRITE(IWR,'(2X,A)') 'WARNING: NSPCG Convergence Failure'
        WRITE(ISC,'(2X,A)') 'Breakdown in Iterate Calculation'
        WRITE(IWR,'(2X,A)') 'Breakdown in Iterate Calculation'
        ICNV = 1
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ELSEIF( IER .EQ. -16 ) THEN
        WRITE(ISC,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(IWR,'(2X,A)') 'ERROR: NSPCG Convergence Failure'
        WRITE(ISC,'(2X,A)') 'Unimplemented Combination of Parameters'
        WRITE(IWR,'(2X,A)') 'Unimplemented Combination of Parameters'
        STOP
      ELSEIF( IER .EQ. -17 ) THEN
        WRITE(ISC,'(2X,A)') 'WARNING: Unable to Perform Eigenvalue Estim
     &ation'
        WRITE(IWR,'(2X,A)') 'WARNING: Unable to Perform Eigenvalue Estim
     &ation'
        ICNV = 1
        ISUB_LOG = ISUB_LOG-1
        RETURN
c      ELSEIF( IER .EQ. 1 ) THEN
c        WRITE(ISC,'(2X,A)') 'WARNING: Failure to Converge in ITMAX Itera
c     &tions'
c        WRITE(IWR,'(2X,A)') 'WARNING: Failure to Converge in ITMAX Itera
c     &tions'
c        ICNV = 1
c        ISUB_LOG = ISUB_LOG-1
c        RETURN
c      ELSEIF( IER .EQ. 2 ) THEN
c        WRITE(ISC,'(2X,A)') 'WARNING: ZETA Too Small - Reset to 500*SREL
c     &PR'
c        WRITE(IWR,'(2X,A)') 'WARNING: ZETA Too Small - Reset to 500*SREL
c     &PR'
c        STOP
c      ELSEIF( IER .EQ. 3 ) THEN
c        WRITE(ISC,'(2X,A)') 'WARNING: ZBRENT Difficulty in Eigenvalue Es
c     &timation'
c        WRITE(IWR,'(2X,A)') 'WARNING: ZBRENT Difficulty in Eigenvalue Es
c     &timation'
c        ICNV = 1
c        ISUB_LOG = ISUB_LOG-1
c        RETURN
c      ELSEIF( IER .EQ. 4 ) THEN
c        WRITE(ISC,'(2X,A)') 'WARNING: ZBRENT Difficulty in Eigenvalue Es
c     &timation'
c        WRITE(IWR,'(2X,A)') 'WARNING: ZBRENT Difficulty in Eigenvalue Es
c     &timation'
c        ICNV = 1
c        ISUB_LOG = ISUB_LOG-1
c        RETURN
c      ELSEIF( IER .EQ. 5 ) THEN
c        WRITE(ISC,'(2X,A)')'WARNING: Negative Pivot Encountered in Facto
c     &rization'
c        WRITE(IWR,'(2X,A)')'WARNING: Negative Pivot Encountered in Facto
c     &rization'
c        ISUB_LOG = ISUB_LOG-1
c        RETURN
      ENDIF
#endif
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PNSPCG group
!
      RETURN
      END


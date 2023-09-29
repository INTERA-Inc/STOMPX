!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBRS
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
!     Row scale the Jacobian matrix elements and solution vector.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, March, 1999.
!     Last Modified by Mark White, Battelle, March 4, 1999.
!     jcbrs.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBRS'
!
!---  Direct banded matrix solver  ---
!
      IF( ILES.EQ.1 ) THEN
        DO 200 K = 1,KFLD
        DO 200 J = 1,JFLD
        DO 200 I = 1,IFLD
          N = ND(I,J,K)
          NMP = IXP(N)
          IF( NMP.LE.0 ) GOTO 200
          DO 190 M = 1,ISVC
            MP = IM(M,NMP)
            DIAG = ALU(MDC,MP)
            BLU(MP) = BLU(MP)/DIAG
            DO 10 L = 1,ISVC
              MCOL = IM(L,NMP)
              MROW = MP-MCOL+MDC
              ALU(MROW,MCOL) = ALU(MROW,MCOL)/DIAG
   10       CONTINUE
            IF( K.NE.1 ) THEN
              NB = N-IJFLD
              IF( IXP(NB).EQ.0 .OR. INBS(1,N).GT.0 ) GOTO 30
              NMB = IXP(NB)
              DO 20 L = 1,ISVC
                MCOL = IM(L,NMB)
                MROW = MP-MCOL+MDC
                ALU(MROW,MCOL) = ALU(MROW,MCOL)/DIAG
   20         CONTINUE
            ENDIF
   30       CONTINUE
            IF( J.NE.1 ) THEN
              NS = N-IFLD
              IF( IXP(NS).EQ.0 .OR. INBS(2,N).GT.0 ) GOTO 50
              NMS = IXP(NS)
              DO 40 L = 1,ISVC
                MCOL = IM(L,NMS)
                MROW = MP-MCOL+MDC
                ALU(MROW,MCOL) = ALU(MROW,MCOL)/DIAG
   40         CONTINUE
            ENDIF
   50       CONTINUE
            IF( I.NE.1 ) THEN
              NW = N-1
              IF( IXP(NW).EQ.0 .OR. INBS(3,N).GT.0 ) GOTO 70
              NMW = IXP(NW)
              DO 60 L = 1,ISVC
                MCOL = IM(L,NMW)
                MROW = MP-MCOL+MDC
                ALU(MROW,MCOL) = ALU(MROW,MCOL)/DIAG
   60         CONTINUE
            ENDIF
   70       CONTINUE
            IF( I.NE.IFLD ) THEN
              NE = N+1
              IF( IXP(NE).EQ.0 .OR. INBS(4,N).GT.0 ) GOTO 90
              NME = IXP(NE)
              DO 80 L = 1,ISVC
                MCOL = IM(L,NME)
                MROW = MP-MCOL+MDC
                ALU(MROW,MCOL) = ALU(MROW,MCOL)/DIAG
   80         CONTINUE
            ENDIF
   90       CONTINUE
            IF( J.NE.JFLD ) THEN
              NN = N+IFLD
              IF( IXP(NN).EQ.0 .OR. INBS(5,N).GT.0 ) GOTO 110
              NMN = IXP(NN)
              DO 100 L = 1,ISVC
                MCOL = IM(L,NMN)
                MROW = MP-MCOL+MDC
                ALU(MROW,MCOL) = ALU(MROW,MCOL)/DIAG
  100         CONTINUE
            ENDIF
  110       CONTINUE
            IF( K.NE.KFLD ) THEN
              NT = N+IJFLD
              IF( IXP(NT).EQ.0 .OR. INBS(6,N).GT.0 ) GOTO 130
              NMT = IXP(NT)
              DO 120 L = 1,ISVC
                MCOL = IM(L,NMT)
                MROW = MP-MCOL+MDC
                ALU(MROW,MCOL) = ALU(MROW,MCOL)/DIAG
  120         CONTINUE
            ENDIF
  130       CONTINUE
  190     CONTINUE
  200   CONTINUE
!
!---  SPLIB Solver  ---
!
      ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
        MA = 0
        DO 400 K = 1,KFLD
        DO 400 J = 1,JFLD
        DO 400 I = 1,IFLD
          N = ND(I,J,K)
          NMP = IXP(N)
          IF( NMP.LE.0 ) GOTO 400
          DO 390 M = 1,ISVC
            MP = IM(M,NMP)
            MCOL = KLU(MP,M+MA)
            DIAG = DLU(MCOL)
            BLU(MP) = BLU(MP)/DIAG
  390     CONTINUE
  400   CONTINUE
      ELSE
        WRITE(ISC,'(A)') 'ERROR: Unknown Linear Equation Solver'
        WRITE(IWR,'(A)') 'ERROR: Unknown Linear Equation Solver'
        STOP
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBRS group  ---
!
      RETURN
      END

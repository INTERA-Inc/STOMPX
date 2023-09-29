!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBT33
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
!     Load the Jacobian matrix for the energy equation with
!     aqueous-phase and gas-phase contributions
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, PNNL, 28 May 2002.
!     Last Modified by MD White, PNNL, 28 May 2002.
!     jcbt_co2e.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVS
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STT(LUK+1),RTP(LUK),RTA(LUK,6),FT(LSFV,6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBT33'
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO 1000 K = 1,KFLD
        DO 1000 J = 1,JFLD
          DO 1000 I = 1,IFLD
        N = ND(I,J,K)
        IF( IXP(N).EQ.0 ) GOTO 1000
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NPX + 1
        NQY = NPY + IFLD
        NQZ = NPZ + IJFLD
        IZN = IZ(N)
!
!---    First-order, forward-difference, time differential
!       assume rock and entrapped water density increases
!       porportionally with increasing porosity  ---
!
        STT1 = (1.D+0-PORT(1,N))*RHOS(IZN)*CPS(IZN)*T(1,N) +
     &    PORD(1,N)*SN(1,N)*RHOSP(1,N)*HSP(1,N) +
     &    (PORT(1,N)-PORD(1,N))*RHOL(1,N)*HL(1,N) +
     &    PORD(1,N)*(SL(1,N)*RHOL(1,N)*HL(1,N) +
     &    SG(1,N)*RHOG(1,N)*UEG(1,N))
        DO 100 M = 1,ISVC+1
          MP = M + 1
          STT0 = (1.D+0-PORT(MP,N))*RHOS(IZN)*CPS(IZN)*T(MP,N) +
     &      PORD(MP,N)*SN(MP,N)*RHOSP(MP,N)*HSP(MP,N) +
     &      (PORT(MP,N)-PORD(MP,N))*RHOL(MP,N)*HL(MP,N) +
     &      PORD(MP,N)*(SL(MP,N)*RHOL(MP,N)*HL(MP,N) +
     &      SG(MP,N)*RHOG(MP,N)*UEG(MP,N))
          STT(M) = (STT0-STT1)*DTI*VOL(N)
  100   CONTINUE
!
!---  Initialize surface fluxes  ---
!
        DO 210 MD = 1,6
          DO 200 M = 1,ISVF
            FT(M,MD) = 0.D+0
  200     CONTINUE
  210   CONTINUE
!
!---  Compute surface fluxes  ---
!
!---  Bottom ---
!
        IF( K.NE.1 ) THEN
          NB = N-IJFLD
          IF( IXP(NB).EQ.0 .OR. INBS(1,N).GT.0 ) GOTO 310
          DO 300 M = 1,ISVF
            FT(M,1) = -AFZ(NPZ)*WQ(M,NPZ)
  300     CONTINUE
        ENDIF
  310   CONTINUE
!
!---  South ---
!
        IF( J.NE.1 ) THEN
          NS = N-IFLD
          IF( IXP(NS).EQ.0 .OR. INBS(2,N).GT.0 ) GOTO 410
          DO 400 M = 1,ISVF
            FT(M,2) = -AFY(NPY)*VQ(M,NPY)
  400     CONTINUE
        ENDIF
  410   CONTINUE
!
!---  West ---
!
        IF( I.NE.1 ) THEN
          NW = N-1
          IF( IXP(NW).EQ.0 .OR. INBS(3,N).GT.0 ) GOTO 510
          DO 500 M = 1,ISVF
            FT(M,3) = -AFX(NPX)*UQ(M,NPX)
  500     CONTINUE
        ENDIF
  510   CONTINUE
!
!---  East ---
!
        IF( I.NE.IFLD ) THEN
          NE = N+1
          IF( IXP(NE).EQ.0 .OR. INBS(4,N).GT.0 ) GOTO 610
          DO 600 M = 1,ISVF
            MF = MFLX(M)
            FT(M,4) = AFX(NQX)*UQ(MF,NQX)
  600     CONTINUE
        ENDIF
  610   CONTINUE
!
!---  North ---
!
        IF( J.NE.JFLD ) THEN
          NN = N+IFLD
          IF( IXP(NN).EQ.0 .OR. INBS(5,N).GT.0 ) GOTO 710
          DO 700 M = 1,ISVF
            MF = MFLX(M)
            FT(M,5) = AFY(NQY)*VQ(MF,NQY)
  700     CONTINUE
        ENDIF
  710   CONTINUE
!
!---  Top ---
!
        IF( K.NE.KFLD ) THEN
          NT = N+IJFLD
          IF( IXP(NT).EQ.0 .OR. INBS(6,N).GT.0 ) GOTO 810
          DO 800 M = 1,ISVF
            MF = MFLX(M)
            FT(M,6) = AFZ(NQZ)*WQ(MF,NQZ)
  800     CONTINUE
        ENDIF
  810   CONTINUE
!
!---  Compute air equation residuals  ---
!
        RTS = STT(1) - SRCT(2,N)
        DO 900 MD = 1,6
          RTS = RTS + FT(1,MD)
  900   CONTINUE
        DO 920 M = 1,ISVC
          RTP(M) = STT(M+1) - SRCT(M+2,N)
          MM = 2*M
          DO 910 MD = 1,6
            RTP(M) = RTP(M) + FT(MM,MD)
  910     CONTINUE
  920   CONTINUE
        DO 940 M = 1,ISVC
          MM = 2*M + 1
          DO 930 MD = 1,6
            RTA(M,MD) = RTS - FT(1,MD) + FT(MM,MD)
  930     CONTINUE
  940   CONTINUE
!
!---  Load Jacobian Matrix  ---
!
        CALL JCBBM( RTS,RTP,RTA,N,I,J,K,IEQT )
!
!---  Continue to Next Node  ---
!
 1000 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBT33 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBC_COMP( IGC )
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
!     Load the Jacobian matrix for the air equation with
!     aqueous-phase and gas-phase contributions
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STC(LUK+1),RAP(LUK),RAA(LUK,6),FC(LSFV,6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBC_COMP'
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO 1000 K = 1,KFLD
        DO 1000 J = 1,JFLD
          DO 1000 I = 1,IFLD
        N = ND(I,J,K)
        IF( IXP(N).LE.0 ) GOTO 1000
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NPX + 1
        NQY = NPY + IFLD
        NQZ = NPZ + IJFLD
        NB = N-IJFLD
        NS = N-IFLD
        NW = N-1
        NE = N+1
        NN = N+IFLD
        NT = N+IJFLD
!
!---    First-order, forward-difference, time differential  ---
!
        STC1 = PORD(1,N)*(XLC(IGC,1,N)*RHOL(1,N)*SL(1,N) +
     &      XGC(IGC,1,N)*RHOG(1,N)*SG(1,N))
        DO 100 M = 1,ISVC+1
          MP = M + 1
          STC0 = PORD(MP,N)*(XLC(IGC,MP,N)*RHOL(MP,N)*SL(MP,N) +
     &      XGC(IGC,MP,N)*RHOG(MP,N)*SG(MP,N))
          STC(M) = (STC0-STC1)*DTI*VOL(N)
  100   CONTINUE
!
!---  Initialize surface fluxes  ---
!
        DO 210 MD = 1,6
          DO 200 M = 1,ISVF
            FC(M,MD) = 0.D+0
  200     CONTINUE
  210   CONTINUE
!
!---  Compute surface fluxes  ---
!
!---  Bottom ---
!
        IF( K.NE.1 ) THEN
          NB = N-IJFLD
          IF( IXP(NB).EQ.0 ) GOTO 310
          DO 300 M = 1,ISVF
            MB = MADJ(M)
            MP = MNOD(M)
            FLCB = XLC(IGC,MB,NB)*RHOL(MB,NB)
            FLCP = XLC(IGC,MP,N)*RHOL(MP,N)
            INDX = 2
            FLC = DIFMN( FLCB,FLCP,DZGF(NB),DZGF(N),WL(1,NPZ),
     &            INDX )
            FGCB = XGC(IGC,MB,NB)*RHOG(MB,NB)
            FGCP = XGC(IGC,MP,N)*RHOG(MP,N)
            INDX = 3
            FGC = DIFMN( FGCB,FGCP,DZGF(NB),DZGF(N),WG(1,NPZ),
     &             INDX )
            FC(M,1) = -AFZ(NPZ)*(WL(M,NPZ)*FLC + WG(M,NPZ)*FGC
     &        + GCPP(1,IGC)*WDLC(IGC,M,NPZ))
            IF( IGC.EQ.1 ) THEN
              FC(M,1) = FC(M,1) -AFZ(NPZ)*(-GCPP(1,IGC)*WDGW(M,NPZ))
              DO JGC = 2,NGC
                FC(M,1) = FC(M,1) -AFZ(NPZ)*(-GCPP(1,JGC)
     &            *WDGC(JGC,M,NPZ))
              END DO
            ELSE
              FC(M,1) = FC(M,1) -AFZ(NPZ)*(GCPP(1,IGC)*WDGC(IGC,M,NPZ))
            ENDIF
  300     CONTINUE
        ENDIF
  310   CONTINUE
!
!---  South ---
!
        IF( J.NE.1 ) THEN
          NS = N-IFLD
          IF( IXP(NS).EQ.0 ) GOTO 410
          DO 400 M = 1,ISVF
            MS = MADJ(M)
            MP = MNOD(M)
            FLCS = XLC(IGC,MS,NS)*RHOL(MS,NS)
            FLCP = XLC(IGC,MP,N)*RHOL(MP,N)
            INDX = 2
            FLC = DIFMN( FLCS,FLCP,DYGF(NS),DYGF(N),VL(1,NPY),
     &          INDX )
            FGCS = XGC(IGC,MS,NS)*RHOG(MS,NS)
            FGCP = XGC(IGC,MP,N)*RHOG(MP,N)
            INDX = 3
            FGC = DIFMN( FGCS,FGCP,DYGF(NS),DYGF(N),VG(1,NPY),
     &          INDX )
            FC(M,2) = -AFY(NPY)*(VL(M,NPY)*FLC + VG(M,NPY)*FGC
     &          + GCPP(1,IGC)*VDLC(IGC,M,NPY))
            IF( IGC.EQ.1 ) THEN
              FC(M,2) = FC(M,2) -AFY(NPY)*(-GCPP(1,IGC)*VDGW(M,NPY))
              DO JGC = 2,NGC
                FC(M,2) = FC(M,2) -AFY(NPY)*(-GCPP(1,JGC)
     &            *VDGC(JGC,M,NPY))
              END DO
            ELSE
              FC(M,2) = FC(M,2) -AFY(NPY)*(GCPP(1,IGC)*VDGC(IGC,M,NPY))
            ENDIF
  400     CONTINUE
        ENDIF
  410   CONTINUE
!
!---  West ---
!
        IF( I.NE.1 ) THEN
          NW = N-1
          IF( IXP(NW).EQ.0 ) GOTO 510
          DO 500 M = 1,ISVF
            MW = MADJ(M)
            MP = MNOD(M)
            FLCW = XLC(IGC,MW,NW)*RHOL(MW,NW)
            FLCP = XLC(IGC,MP,N)*RHOL(MP,N)
            INDX = 2
            FLC = DIFMN( FLCW,FLCP,DXGF(NW),DXGF(N),UL(1,NPX),
     &          INDX )
            FGCW = XGC(IGC,MW,NW)*RHOG(MW,NW)
            FGCP = XGC(IGC,MP,N)*RHOG(MP,N)
            INDX = 3
            FGC = DIFMN( FGCW,FGCP,DXGF(NW),DXGF(N),UG(1,NPX),
     &         INDX )
            FC(M,3) = -AFX(NPX)*(UL(M,NPX)*FLC + UG(M,NPX)*FGC
     &         + GCPP(1,IGC)*UDLC(IGC,M,NPX))
            IF( IGC.EQ.1 ) THEN
              FC(M,3) = FC(M,3) -AFX(NPX)*(-GCPP(1,IGC)*UDGW(M,NPX))
              DO JGC = 2,NGC
                FC(M,3) = FC(M,3) -AFX(NPX)*(-GCPP(1,JGC)
     &            *UDGC(JGC,M,NPX))
              END DO
            ELSE
              FC(M,3) = FC(M,3) -AFX(NPX)*(GCPP(1,IGC)*UDGC(IGC,M,NPX))
            ENDIF
  500     CONTINUE
        ENDIF
  510   CONTINUE
!
!---  East ---
!
        IF( I.NE.IFLD ) THEN
          NE = N+1
          IF( IXP(NE).EQ.0 ) GOTO 610
          DO 600 M = 1,ISVF
            ME = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FLCE = XLC(IGC,ME,NE)*RHOL(ME,NE)
            FLCP = XLC(IGC,MP,N)*RHOL(MP,N)
            INDX = 2
            FLC = DIFMN( FLCP,FLCE,DXGF(N),DXGF(NE),UL(1,NQX),
     &        INDX )
            FGCE = XGC(IGC,ME,NE)*RHOG(ME,NE)
            FGCP = XGC(IGC,MP,N)*RHOG(MP,N)
            INDX = 3
            FGC = DIFMN( FGCP,FGCE,DXGF(N),DXGF(NE),UG(1,NQX),
     &        INDX )
            FC(M,4) = AFX(NQX)*(UL(MF,NQX)*FLC + UG(MF,NQX)*FGC
     &        + GCPP(1,IGC)*UDLC(IGC,MF,NQX))
            IF( IGC.EQ.1 ) THEN
              FC(M,4) = FC(M,4) + AFX(NQX)*(-GCPP(1,IGC)*UDGW(MF,NQX))
              DO JGC = 2,NGC
                FC(M,4) = FC(M,4) + AFX(NQX)*(-GCPP(1,JGC)
     &            *UDGC(JGC,MF,NQX))
              END DO
            ELSE
              FC(M,4) = FC(M,4) + AFX(NQX)*(GCPP(1,IGC)
     &          *UDGC(IGC,MF,NQX))
            ENDIF
  600     CONTINUE
        ENDIF
  610   CONTINUE
!
!---  North ---
!
        IF( J.NE.JFLD ) THEN
          NN = N+IFLD
          IF( IXP(NN).EQ.0 ) GOTO 710
          DO 700 M = 1,ISVF
            MN = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FLCN = XLC(IGC,MN,NN)*RHOL(MN,NN)
            FLCP = XLC(IGC,MP,N)*RHOL(MP,N)
            INDX = 2
            FLC = DIFMN( FLCP,FLCN,DYGF(N),DYGF(NN),VL(1,NQY),
     &        INDX )
            FGCN = XGC(IGC,MN,NN)*RHOG(MN,NN)
            FGCP = XGC(IGC,MP,N)*RHOG(MP,N)
            INDX = 3
            FGC = DIFMN( FGCP,FGCN,DYGF(N),DYGF(NN),VG(1,NQY),
     &        INDX )
            FC(M,5) = AFY(NQY)*(VL(MF,NQY)*FLC + VG(MF,NQY)*FGC
     &         + GCPP(1,IGC)*VDLC(IGC,MF,NQY))
            IF( IGC.EQ.1 ) THEN
              FC(M,5) = FC(M,5) + AFY(NQY)*(-GCPP(1,IGC)*VDGW(MF,NQY))
              DO JGC = 2,NGC
                FC(M,5) = FC(M,5) + AFY(NQY)*(-GCPP(1,JGC)
     &            *VDGC(JGC,MF,NQY))
              END DO
            ELSE
              FC(M,5) = FC(M,5) + AFY(NQY)*(GCPP(1,IGC)*
     &          VDGC(IGC,MF,NQY))
            ENDIF
  700     CONTINUE
        ENDIF
  710   CONTINUE
!
!---  Top ---
!
        IF( K.NE.KFLD ) THEN
          NT = N+IJFLD
          IF( IXP(NT).EQ.0 ) GOTO 810
          DO 800 M = 1,ISVF
            MT = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FLCT = XLC(IGC,MT,NT)*RHOL(MT,NT)
            FLCP = XLC(IGC,MP,N)*RHOL(MP,N)
            INDX = 2
            FLC = DIFMN( FLCP,FLCT,DZGF(N),DZGF(NT),WL(1,NQZ),
     &         INDX )
            FGCT = XGC(IGC,MT,NT)*RHOG(MT,NT)
            FGCP = XGC(IGC,MP,N)*RHOG(MP,N)
            INDX = 3
            FGC = DIFMN( FGCP,FGCT,DZGF(N),DZGF(NT),WG(1,NQZ),
     &         INDX )
            FC(M,6) = AFZ(NQZ)*(WL(MF,NQZ)*FLC + WG(MF,NQZ)*FGC
     &         + GCPP(1,IGC)*WDLC(IGC,MF,NQZ))
            IF( IGC.EQ.1 ) THEN
              FC(M,6) = FC(M,6) + AFZ(NQZ)*(-GCPP(1,IGC)*WDGW(MF,NQZ))
              DO JGC = 2,NGC
                FC(M,6) = FC(M,6) + AFZ(NQZ)*(-GCPP(1,JGC)
     &            *WDGC(JGC,MF,NQZ))
              END DO
            ELSE
              FC(M,6) = FC(M,6) + AFZ(NQZ)*(GCPP(1,IGC)
     &          *WDGC(IGC,MF,NQZ))
            ENDIF
  800     CONTINUE
        ENDIF
  810   CONTINUE
!
!---  Compute air equation residuals  ---
!
        RAS = STC(1) - SRCGC(IGC,2,N)
        DO 900 MD = 1,6
          RAS = RAS + FC(1,MD)
  900   CONTINUE
        DO 920 M = 1,ISVC
          RAP(M) = STC(M+1) - SRCGC(IGC,M+2,N)
          MM = 2*M
          DO 910 MD = 1,6
            RAP(M) = RAP(M) + FC(MM,MD)
  910     CONTINUE
  920   CONTINUE
        DO 940 M = 1,ISVC
          MM = 2*M + 1
          DO 930 MD = 1,6
            RAA(M,MD) = RAS - FC(1,MD) + FC(MM,MD)
  930     CONTINUE
  940   CONTINUE
!
!---  Load Jacobian Matrix  ---
!
        CALL JCBBM( RAS,RAP,RAA,N,I,J,K,IEQGC(IGC) )
!
!---  Continue to Next Node  ---
!
 1000 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBC_COMP group
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBS_COMP
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
!     Load the Jacobian matrix for the salt equation aqueous-phase
!     contributions (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STS(LUK+1),RSP(LUK),RSA(LUK,6),FS(LSFV,6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBS_COMP'
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO 1000 K = 1,KFLD
        DO 1000 J = 1,JFLD
          DO 1000 I = 1,IFLD
        N = ND(I,J,K)
        IF( IXP(N).LE.0 ) GOTO 1000
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NPX + 1
        NQY = NPY + IFLD
        NQZ = NPZ + IJFLD
        DO 10 M = 1,ISVC+1
          MP = M + 1
          STS(M) = (TMS(MP,N)-TMS(1,N))*VOL(N)*DTI
   10   CONTINUE
!
!---  Initialize surface fluxes  ---
!
        DO 30 MD = 1,6
          DO 20 M = 1,ISVF
            FS(M,MD) = 0.D+0
   20     CONTINUE
   30   CONTINUE
!
!---  Bottom surface flux component  ---
!
        IF( K.NE.1 ) THEN
          NB = N-IJFLD
          IF( IXP(NB).EQ.0 ) GOTO 110
          DO 100 M = 1,ISVF
            FS(M,1) = -AFZ(NPZ)*WS(M,NPZ)
  100     CONTINUE
        ENDIF
  110   CONTINUE
!
!---  South surface flux component ---
!
        IF( J.NE.1 ) THEN
          NS = N-IFLD
          IF( IXP(NS).EQ.0 ) GOTO 210
          DO 200 M = 1,ISVF
            FS(M,2) = -AFY(NPY)*VS(M,NPY)
  200     CONTINUE
        ENDIF
  210   CONTINUE
!
!---  West surface flux component ---
!
        IF( I.NE.1 ) THEN
          NW = N-1
          IF( IXP(NW).EQ.0 ) GOTO 310
          DO 300 M = 1,ISVF
            FS(M,3) = -AFX(NPX)*US(M,NPX)
  300     CONTINUE
        ENDIF
  310   CONTINUE
!
!---  East surface flux component  ---
!
        IF( I.NE.IFLD ) THEN
          NE = N+1
          IF( IXP(NE).EQ.0 ) GOTO 410
          DO 400 M = 1,ISVF
            MF = MFLX(M)
            FS(M,4) = AFX(NQX)*US(MF,NQX)
  400     CONTINUE
        ENDIF
  410   CONTINUE
!
!---  North surface flux component  ---
!
        IF( J.NE.JFLD ) THEN
          NN = N+IFLD
          IF( IXP(NN).EQ.0 ) GOTO 510
          DO 500 M = 1,ISVF
            MF = MFLX(M)
            FS(M,5) = AFY(NQY)*VS(MF,NQY)
  500     CONTINUE
        ENDIF
  510   CONTINUE
!
!---  Top surface flux component  ---
!
        IF( K.NE.KFLD ) THEN
          NT = N+IJFLD
          IF( IXP(NT).EQ.0 ) GOTO 610
          DO 600 M = 1,ISVF
            MF = MFLX(M)
            FS(M,6) = AFZ(NQZ)*WS(MF,NQZ)
  600     CONTINUE
        ENDIF
  610   CONTINUE
!
!---  Salt equation residuals  ---
!
        RSS = STS(1) - SRCS(2,N)
        DO 900 MD = 1,6
          RSS = RSS + FS(1,MD)
  900   CONTINUE
        DO 920 M = 1,ISVC
          RSP(M) = STS(M+1) - SRCS(M+2,N)
          MM = 2*M
          DO 910 MD = 1,6
            RSP(M) = RSP(M) + FS(MM,MD)
  910     CONTINUE
  920   CONTINUE
        DO 940 M = 1,ISVC
          MM = 2*M + 1
          DO 930 MD = 1,6
            RSA(M,MD) = RSS - FS(1,MD) + FS(MM,MD)
  930     CONTINUE
  940   CONTINUE
!
!---  Jacobian matrix loader  --
!
        CALL JCBBM( RSS,RSP,RSA,N,I,J,K,IEQS )
!
!---  Continue to Next Node  ---
!
 1000 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBS_COMP group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBT_COMP
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
!     Written by DH Bacon, PNNL, 31 October 2012!
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
      REAL*8 STTX(LUK+1),RTP(LUK),RTA(LUK,6),FT(LSFV,6)
!     REAL*8 UEL(LSV)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBT_COMP'
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO 1000 K = 1,KFLD
        DO 1000 J = 1,JFLD
          DO 1000 I = 1,IFLD
        N = ND(I,J,K)
        IF( IXP(N).LE.0 ) GOTO 1000
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
!       UEL(1) = HL(1,N) - (PG(1,N)+PATM)/RHOL(1,N)
        STT1 = (1.D+0-PORT(1,N))*RHOS(IZN)*CPS(IZN)*T(1,N) +
     &    PORD(1,N)*SS(1,N)*RHOSP(1,N)*HSP(1,N) +
     &    (PORT(1,N)-PORD(1,N))*RHOL(1,N)*HL(1,N) +
     &    PORD(1,N)*(SL(1,N)*RHOL(1,N)*HL(1,N) +
     &    SG(1,N)*RHOG(1,N)*UEG(1,N))
        DO 100 M = 1,ISVC+1
          MP = M + 1
!         UEL(MP) = HL(MP,N) - (PG(MP,N)+PATM)/RHOL(MP,N)
          STT0 = (1.D+0-PORT(MP,N))*RHOS(IZN)*CPS(IZN)*T(MP,N) +
     &      PORD(MP,N)*SS(MP,N)*RHOSP(MP,N)*HSP(MP,N) +
     &      (PORT(MP,N)-PORD(MP,N))*RHOL(MP,N)*HL(MP,N) +
     &      PORD(MP,N)*(SL(MP,N)*RHOL(MP,N)*HL(MP,N) +
     &      SG(MP,N)*RHOG(MP,N)*UEG(MP,N))
          STTX(M) = (STT0-STT1)*DTI*VOL(N)
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
          IF( IXP(NB).EQ.0 ) GOTO 310
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
          IF( IXP(NS).EQ.0 ) GOTO 410
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
          IF( IXP(NW).EQ.0 ) GOTO 510
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
          IF( IXP(NE).EQ.0 ) GOTO 610
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
          IF( IXP(NN).EQ.0 ) GOTO 710
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
          IF( IXP(NT).EQ.0 ) GOTO 810
          DO 800 M = 1,ISVF
            MF = MFLX(M)
            FT(M,6) = AFZ(NQZ)*WQ(MF,NQZ)
  800     CONTINUE
        ENDIF
  810   CONTINUE
!
!---  Compute air equation residuals  ---
!
        RTS = STTX(1) - SRCT(2,N)
        DO 900 MD = 1,6
          RTS = RTS + FT(1,MD)
  900   CONTINUE
        DO 920 M = 1,ISVC
          RTP(M) = STTX(M+1) - SRCT(M+2,N)
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
!---  End of JCBT_COMP group
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBW_COMP
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
!     Load the Jacobian matrix for the water equation with
!     aqueous-phase, gas-phase and salt contributions
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE FDVI
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STWX(LUK+1),RWP(LUK),RWA(LUK,6),FW(LSFV,6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBW_COMP'
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO 1000 K = 1,KFLD
        DO 1000 J = 1,JFLD
          DO 1000 I = 1,IFLD
        N = ND(I,J,K)
        IF( IXP(N).LE.0 ) GOTO 1000
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NPX + 1
        NQY = NPY + IFLD
        NQZ = NPZ + IJFLD
!
!---    First-order, forward-difference, time differential  ---
!
        STW1 = PORD(1,N)*(XLW(1,N)*RHOL(1,N)*SL(1,N) +
     &    XGW(1,N)*RHOG(1,N)*SG(1,N))
        DO 100 M = 1,ISVC+1
          MP = M + 1
          STW0 = PORD(MP,N)*(XLW(MP,N)*RHOL(MP,N)*SL(MP,N) +
     &      XGW(MP,N)*RHOG(MP,N)*SG(MP,N))
          STWX(M) = (STW0-STW1)*DTI*VOL(N)
  100   CONTINUE
!
!---  Freezing conditions  ---
!
!       IF( ISLC(5).EQ.1 ) THEN
!
!---    First-order, forward-difference, time differential  ---
!
         DO 104 M = 1,ISVC+1
           MP = M + 1
           STWX(M) = STWX(M) + (PORD(MP,N)*RHOI(MP,N)*SI(MP,N) -
     &       PORD(1,N)*RHOI(1,N)*SI(1,N))*DTI*VOL(N)
  104    CONTINUE
!
!---  Initialize surface fluxes  ---
!
        DO 210 MD = 1,6
          DO 200 M = 1,ISVF
            FW(M,MD) = 0.D+0
  200     CONTINUE
  210   CONTINUE
!
!---  Compute surface fluxes  ---
!
!---  Bottom ---
!
        IF( K.NE.1 ) THEN
          NB = N-IJFLD
          IF( IXP(NB).EQ.0 ) GOTO 310
          DO 300 M = 1,ISVF
            MB = MADJ(M)
            MP = MNOD(M)
            FLWB = XLW(MB,NB)*RHOL(MB,NB)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWB,FLWP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
            FGWB = XGW(MB,NB)*RHOG(MB,NB)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            INDX = 3
            FGW = DIFMN( FGWB,FGWP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
            FW(M,1) = -AFZ(NPZ)*(WL(M,NPZ)*FLW + WG(M,NPZ)*FGW
     &        + WTMW*WDGW(M,NPZ) - WTMW*WDS(M,NPZ)/WTMS)
            DO IGC = 1,NGC
              FW(M,1) = FW(M,1) - AFZ(NPZ)*(-WTMW*WDLC(IGC,M,NPZ))
            ENDDO
  300     CONTINUE
        ENDIF
  310   CONTINUE
!
!---  South ---
!
        IF( J.NE.1 ) THEN
          NS = N-IFLD
          IF( IXP(NS).EQ.0 ) GOTO 410
          DO 400 M = 1,ISVF
            MS = MADJ(M)
            MP = MNOD(M)
            FLWS = XLW(MS,NS)*RHOL(MS,NS)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWS,FLWP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
            FGWS = XGW(MS,NS)*RHOG(MS,NS)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            INDX = 3
            FGW = DIFMN( FGWS,FGWP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
            FW(M,2) = -AFY(NPY)*(VL(M,NPY)*FLW + VG(M,NPY)*FGW
     &        + WTMW*VDGW(M,NPY) - WTMW*VDS(M,NPY)/WTMS)
            DO IGC = 1,NGC
              FW(M,2) = FW(M,2) - AFY(NPY)*(-WTMW*VDLC(IGC,M,NPY))
            ENDDO
  400     CONTINUE
        ENDIF
  410   CONTINUE
!
!---  West ---
!
        IF( I.NE.1 ) THEN
          NW = N-1
          IF( IXP(NW).EQ.0 ) GOTO 510
          DO 500 M = 1,ISVF
            MW = MADJ(M)
            MP = MNOD(M)
            FLWW = XLW(MW,NW)*RHOL(MW,NW)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWW,FLWP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            FGWW = XGW(MW,NW)*RHOG(MW,NW)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            INDX = 3
            FGW = DIFMN( FGWW,FGWP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            FW(M,3) = -AFX(NPX)*(UL(M,NPX)*FLW + UG(M,NPX)*FGW
     &        + WTMW*UDGW(M,NPX) - WTMW*UDS(M,NPX)/WTMS)
            DO IGC = 1,NGC
              FW(M,3) = FW(M,3) - AFX(NPX)*(-WTMW*UDLC(IGC,M,NPX))
            ENDDO
  500     CONTINUE
        ENDIF
  510   CONTINUE
!
!---  East ---
!
        IF( I.NE.IFLD ) THEN
          NE = N+1
          IF( IXP(NE).EQ.0 ) GOTO 610
          DO 600 M = 1,ISVF
            ME = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FLWE = XLW(ME,NE)*RHOL(ME,NE)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWP,FLWE,DXGF(N),DXGF(NE),UL(1,NQX),INDX )
            FGWE = XGW(ME,NE)*RHOG(ME,NE)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            INDX = 3
            FGW = DIFMN( FGWP,FGWE,DXGF(N),DXGF(NE),UG(1,NQX),INDX )
            FW(M,4) = AFX(NQX)*(UL(MF,NQX)*FLW + UG(MF,NQX)*FGW
     &        + WTMW*UDGW(MF,NQX) - WTMW*UDS(MF,NQX)/WTMS)
            DO IGC = 1,NGC
              FW(M,4) = FW(M,4) + AFX(NQX)*(-WTMW*UDLC(IGC,MF,NQX))
            ENDDO
  600     CONTINUE
        ENDIF
  610   CONTINUE
!
!---  North ---
!
        IF( J.NE.JFLD ) THEN
          NN = N+IFLD
          IF( IXP(NN).EQ.0 ) GOTO 710
          DO 700 M = 1,ISVF
            MN = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FLWN = XLW(MN,NN)*RHOL(MN,NN)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWP,FLWN,DYGF(N),DYGF(NN),VL(1,NQY),INDX )
            FGWN = XGW(MN,NN)*RHOG(MN,NN)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            INDX = 3
            FGW = DIFMN( FGWP,FGWN,DYGF(N),DYGF(NN),VG(1,NQY),INDX )
            FW(M,5) = AFY(NQY)*(VL(MF,NQY)*FLW + VG(MF,NQY)*FGW
     &        + WTMW*VDGW(MF,NQY) - WTMW*VDS(MF,NQY)/WTMS)
            DO IGC = 1,NGC
              FW(M,5) = FW(M,5) + AFY(NQY)*(-WTMW*VDLC(IGC,MF,NQY))
            ENDDO
  700     CONTINUE
        ENDIF
  710   CONTINUE
!
!---  Top ---
!
        IF( K.NE.KFLD ) THEN
          NT = N+IJFLD
          IF( IXP(NT).EQ.0 ) GOTO 810
          DO 800 M = 1,ISVF
            MT = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FLWT = XLW(MT,NT)*RHOL(MT,NT)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWP,FLWT,DZGF(N),DZGF(NT),WL(1,NQZ),INDX )
            FGWT = XGW(MT,NT)*RHOG(MT,NT)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            INDX = 3
            FGW = DIFMN( FGWP,FGWT,DZGF(N),DZGF(NT),WG(1,NQZ),INDX )
            FW(M,6) = AFZ(NQZ)*(WL(MF,NQZ)*FLW + WG(MF,NQZ)*FGW
     &        + WTMW*WDGW(MF,NQZ) - WTMW*WDS(MF,NQZ)/WTMS)
            DO IGC = 1,NGC
              FW(M,6) = FW(M,6) + AFZ(NQZ)*(-WTMW*WDLC(IGC,MF,NQZ))
            ENDDO
  800     CONTINUE
        ENDIF
  810   CONTINUE
!
!---  Compute water equation residuals  ---
!
        RWS = STWX(1) - SRCW(2,N)
        DO 900 MD = 1,6
          RWS = RWS + FW(1,MD)
  900   CONTINUE
        DO 920 M = 1,ISVC
          RWP(M) = STWX(M+1) - SRCW(M+2,N)
          MM = 2*M
          DO 910 MD = 1,6
            RWP(M) = RWP(M) + FW(MM,MD)
  910     CONTINUE
  920   CONTINUE
        DO 940 M = 1,ISVC
          MM = 2*M + 1
          DO 930 MD = 1,6
            RWA(M,MD) = RWS - FW(1,MD) + FW(MM,MD)
  930     CONTINUE
  940   CONTINUE
!
!---  Load Jacobian Matrix  ---
!
        CALL JCBBM( RWS,RWP,RWA,N,I,J,K,IEQW )
!
!---  Continue to Next Node  ---
!
 1000 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBW_COMP group
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGCB_COMP( IGC,N,NB,NPZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, air equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE CCP
      USE BCVP
      USE BCVGC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGCB_COMP'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLCB = XGCB(IGC,MB,NB)*RHOGB(MB,NB)
        FLCP = XGC(IGC,MB,N)*RHOG(MB,N)
        INDX = 3
        FLC = DIFMN( FLCB,FLCP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
        RS(M) = -AFZ(NPZ)*WG(MP,NPZ)*FLC
        IF( IGC.EQ.1 ) THEN
          RS(M) = RS(M) -AFZ(NPZ)*(-GCPP(1,IGC)*WDGW(MP,NPZ))
          DO JGC = 2,NGC
            RS(M) = RS(M) -AFZ(NPZ)*(-GCPP(1,JGC)*WDGC(JGC,MP,NPZ))
          END DO
        ELSE
          RS(M) = RS(M) -AFZ(NPZ)*(GCPP(1,IGC)*WDGC(IGC,MP,NPZ))
        ENDIF
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQGC(IGC) )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGCB_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGCS_COMP( IGC,N,NB,NPY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, air equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE CCP
      USE BCVP
      USE BCVGC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGCS_COMP'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLCB = XGCB(IGC,MB,NB)*RHOGB(MB,NB)
        FLCP = XGC(IGC,MB,N)*RHOG(MB,N)
        INDX = 3
        FLC = DIFMN( FLCB,FLCP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
        RS(M) = -AFY(NPY)*VG(MP,NPY)*FLC
        IF( IGC.EQ.1 ) THEN
          RS(M) = RS(M) -AFY(NPY)*(-GCPP(1,IGC)*VDGW(MP,NPY))
          DO JGC = 2,NGC
            RS(M) = RS(M) -AFY(NPY)*(-GCPP(1,JGC)*VDGC(JGC,MP,NPY))
          END DO
        ELSE
          RS(M) = RS(M) -AFY(NPY)*(GCPP(1,IGC)*VDGC(IGC,MP,NPY))
        ENDIF
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQGC(IGC) )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGCS_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGCW_COMP( IGC,N,NB,NPX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, air equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE CCP
      USE BCVP
      USE BCVGC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGCW_COMP'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLCB = XGCB(IGC,MB,NB)*RHOGB(MB,NB)
        FLCP = XGC(IGC,MB,N)*RHOG(MB,N)
        INDX = 3
        FLC = DIFMN( FLCB,FLCP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
        RS(M) = -AFX(NPX)*UG(MP,NPX)*FLC
        IF( IGC.EQ.1 ) THEN
          RS(M) = RS(M) -AFX(NPX)*(-GCPP(1,IGC)*UDGW(MP,NPX))
          DO JGC = 2,NGC
            RS(M) = RS(M) -AFX(NPX)*(-GCPP(1,JGC)*UDGC(JGC,MP,NPX))
          END DO
        ELSE
          RS(M) = RS(M) -AFX(NPX)*(GCPP(1,IGC)*UDGC(IGC,MP,NPX))
        ENDIF
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQGC(IGC) )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGCW_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGCE_COMP( IGC,N,NB,NQX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, air equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE CCP
      USE BCVP
      USE BCVGC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGCE_COMP'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLCB = XGCB(IGC,MB,NB)*RHOGB(MB,NB)
        FLCP = XGC(IGC,MB,N)*RHOG(MB,N)
        INDX = 3
        FLC = DIFMN( FLCP,FLCB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
        RS(M) = AFX(NQX)*UG(MN,NQX)*FLC 
        IF( IGC.EQ.1 ) THEN
          RS(M) = RS(M) + AFX(NQX)*(-GCPP(1,IGC)*UDGW(MN,NQX))
          DO JGC = 2,NGC
            RS(M) = RS(M) + AFX(NQX)*(-GCPP(1,JGC)*UDGC(JGC,MN,NQX))
          END DO
        ELSE
          RS(M) = RS(M) + AFX(NQX)*(GCPP(1,IGC)*UDGC(IGC,MN,NQX))
        ENDIF
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQGC(IGC) )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGCE_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGCN_COMP( IGC,N,NB,NQY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, air equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE CCP
      USE BCVP
      USE BCVGC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGCN_COMP'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLCB = XGCB(IGC,MB,NB)*RHOGB(MB,NB)
        FLCP = XGC(IGC,MB,N)*RHOG(MB,N)
        INDX = 3
        FLC = DIFMN( FLCP,FLCB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
        RS(M) = AFY(NQY)*VG(MN,NQY)*FLC
        IF( IGC.EQ.1 ) THEN
          RS(M) = RS(M) + AFY(NQY)*(-GCPP(1,IGC)*VDGW(MN,NQY))
          DO JGC = 2,NGC
            RS(M) = RS(M) + AFY(NQY)*(-GCPP(1,JGC)*VDGC(JGC,MN,NQY))
          END DO
        ELSE
          RS(M) = RS(M) + AFY(NQY)*(GCPP(1,IGC)*VDGC(IGC,MN,NQY))
        ENDIF
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQGC(IGC) )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGCN_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGCT_COMP( IGC,N,NB,NQZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, air equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE CCP
      USE BCVP
      USE BCVGC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGCT_COMP'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLCB = XGCB(IGC,MB,NB)*RHOGB(MB,NB)
        FLCP = XGC(IGC,MB,N)*RHOG(MB,N)
        INDX = 3
        FLC = DIFMN( FLCP,FLCB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
        RS(M) = AFZ(NQZ)*WG(MN,NQZ)*FLC
        IF( IGC.EQ.1 ) THEN
          RS(M) = RS(M) + AFZ(NQZ)*(-GCPP(1,IGC)*WDGW(MN,NQZ))
          DO JGC = 2,NGC
            RS(M) = RS(M) + AFZ(NQZ)*(-GCPP(1,JGC)*WDGC(JGC,MN,NQZ))
          END DO
        ELSE
          RS(M) = RS(M) + AFZ(NQZ)*(GCPP(1,IGC)*WDGC(IGC,MN,NQZ))
        ENDIF
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQGC(IGC) )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGCT_COMP group  ---
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGWB_COMP( N,NB,NPZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, water equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGWB_COMP'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLWB = XGWB(MB,NB)*RHOGB(MB,NB)
        FLWP = XGW(MB,N)*RHOG(MB,N)
        INDX = 3
        FLW = DIFMN( FLWB,FLWP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
        RS(M) = -AFZ(NPZ)*(WG(MP,NPZ)*FLW + WTMW*WDGW(MP,NPZ))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGWB_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGWS_COMP( N,NB,NPY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, water equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGWS_COMP'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLWB = XGWB(MB,NB)*RHOGB(MB,NB)
        FLWP = XGW(MB,N)*RHOG(MB,N)
        INDX = 3
        FLW = DIFMN( FLWB,FLWP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
        RS(M) = -AFY(NPY)*(VG(MP,NPY)*FLW + WTMW*VDGW(MP,NPY))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGWS_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGWW_COMP( N,NB,NPX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, water equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGWW_COMP'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLWB = XGWB(MB,NB)*RHOGB(MB,NB)
        FLWP = XGW(MB,N)*RHOG(MB,N)
        INDX = 3
        FLW = DIFMN( FLWB,FLWP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
        RS(M) = -AFX(NPX)*(UG(MP,NPX)*FLW + WTMW*UDGW(MP,NPX))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGWW_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGWE_COMP( N,NB,NQX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, water equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGWE_COMP'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLWB = XGWB(MB,NB)*RHOGB(MB,NB)
        FLWP = XGW(MB,N)*RHOG(MB,N)
        INDX = 3
        FLW = DIFMN( FLWP,FLWB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
        RS(M) = AFX(NQX)*(UG(MN,NQX)*FLW + WTMW*UDGW(MN,NQX))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGWE_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGWN_COMP( N,NB,NQY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, water equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGWN_COMP'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLWB = XGWB(MB,NB)*RHOGB(MB,NB)
        FLWP = XGW(MB,N)*RHOG(MB,N)
        INDX = 3
        FLW = DIFMN( FLWP,FLWB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
        RS(M) = AFY(NQY)*(VG(MN,NQY)*FLW + WTMW*VDGW(MN,NQY))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGWN_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGWT_COMP( N,NB,NQZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, water equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGWT_COMP'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLWB = XGWB(MB,NB)*RHOGB(MB,NB)
        FLWP = XGW(MB,N)*RHOG(MB,N)
        INDX = 3
        FLW = DIFMN( FLWP,FLWB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
        RS(M) = AFZ(NQZ)*(WG(MN,NQZ)*FLW + WTMW*WDGW(MN,NQZ))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGWT_COMP group  ---
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLCB_COMP( IGC,N,NB,NPZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, air equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE CCP
      USE BCVP
      USE BCVGC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLCB_COMP'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLCB = XLCB(IGC,MB,NB)*RHOLB(MB,NB)
        FLCP = XLC(IGC,MB,N)*RHOL(MB,N)
        INDX = 2
        FLC = DIFMN( FLCB,FLCP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
        RS(M) = -AFZ(NPZ)*(WL(MP,NPZ)*FLC+WDLC(IGC,MP,NPZ)*GCPP(1,IGC))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQGC(IGC) )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLCB_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLCS_COMP( IGC,N,NB,NPY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, air equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE CCP
      USE BCVP
      USE BCVGC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLCS_COMP'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLCB = XLCB(IGC,MB,NB)*RHOLB(MB,NB)
        FLCP = XLC(IGC,MB,N)*RHOL(MB,N)
        INDX = 2
        FLC = DIFMN( FLCB,FLCP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
        RS(M) = -AFY(NPY)*(VL(MP,NPY)*FLC+VDLC(IGC,MP,NPY)*GCPP(1,IGC))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQGC(IGC) )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLCS_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLCW_COMP( IGC,N,NB,NPX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, air equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE CCP
      USE BCVP
      USE BCVGC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLCW_COMP'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLCB = XLCB(IGC,MB,NB)*RHOLB(MB,NB)
        FLCP = XLC(IGC,MB,N)*RHOL(MB,N)
        INDX = 2
        FLC = DIFMN( FLCB,FLCP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
        RS(M) = -AFX(NPX)*(UL(MP,NPX)*FLC+UDLC(IGC,MP,NPX)*GCPP(1,IGC))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQGC(IGC) )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLCW_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLCE_COMP( IGC,N,NB,NQX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, air equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE CCP
      USE BCVP
      USE BCVGC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLCE_COMP'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLCB = XLCB(IGC,MB,NB)*RHOLB(MB,NB)
        FLCP = XLC(IGC,MB,N)*RHOL(MB,N)
        INDX = 2
        FLC = DIFMN( FLCP,FLCB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
        RS(M) = AFX(NQX)*(UL(MN,NQX)*FLC+UDLC(IGC,MN,NQX)*GCPP(1,IGC))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQGC(IGC) )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLCE_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLCN_COMP( IGC,N,NB,NQY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, air equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE CCP
      USE BCVP
      USE BCVGC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLCN_COMP'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLCB = XLCB(IGC,MB,NB)*RHOLB(MB,NB)
        FLCP = XLC(IGC,MB,N)*RHOL(MB,N)
        INDX = 2
        FLC = DIFMN( FLCP,FLCB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
        RS(M) = AFY(NQY)*(VL(MN,NQY)*FLC+VDLC(IGC,MN,NQY)*GCPP(1,IGC))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQGC(IGC) )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLCN_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLCT_COMP( IGC,N,NB,NQZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, air equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE CCP
      USE BCVP
      USE BCVGC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLCT_COMP'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLCB = XLCB(IGC,MB,NB)*RHOLB(MB,NB)
        FLCP = XLC(IGC,MB,N)*RHOL(MB,N)
        INDX = 2
        FLC = DIFMN( FLCP,FLCB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
        RS(M) = AFZ(NQZ)*(WL(MN,NQZ)*FLC+WDLC(IGC,MN,NQZ)*GCPP(1,IGC))
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQGC(IGC) )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLCT_COMP group  ---
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWB_COMP( N,NB,NPZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, water equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLWB_COMP'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLWB = XLWB(MB,NB)*RHOLB(MB,NB)
        FLWP = XLW(MB,N)*RHOL(MB,N)
        INDX = 2
        FLW = DIFMN( FLWB,FLWP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
        WDLA(MP,NPZ) = 0.D+0
        DO IGC = 1, NGC
         WDLA(MP,NPZ) = WDLA(MP,NPZ) + WDLC(IGC,MP,NPZ)
        ENDDO
        RS(M) = -AFZ(NPZ)*(WL(MP,NPZ)*FLW - WTMW*WDLA(MP,NPZ)
     &    - WTMW*WDS(MP,NPZ)/WTMS)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWB_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWS_COMP( N,NB,NPY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, water equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLWS_COMP'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLWB = XLWB(MB,NB)*RHOLB(MB,NB)
        FLWP = XLW(MB,N)*RHOL(MB,N)
        INDX = 2
        FLW = DIFMN( FLWB,FLWP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
        VDLA(MP,NPY) = 0.D+0
        DO IGC = 1, NGC
         VDLA(MP,NPY) = VDLA(MP,NPY) + VDLC(IGC,MP,NPY)
        ENDDO
        RS(M) = -AFY(NPY)*(VL(MP,NPY)*FLW - WTMW*VDLA(MP,NPY)
     &    - WTMW*VDS(MP,NPY)/WTMS)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWS_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWW_COMP( N,NB,NPX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, water equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLWW_COMP'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLWB = XLWB(MB,NB)*RHOLB(MB,NB)
        FLWP = XLW(MB,N)*RHOL(MB,N)
        INDX = 2
        FLW = DIFMN( FLWB,FLWP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
        UDLA(MP,NPX) = 0.D+0
        DO IGC = 1, NGC
         UDLA(MP,NPX) = UDLA(MP,NPX) + UDLC(IGC,MP,NPX)
        ENDDO
        RS(M) = -AFX(NPX)*(UL(MP,NPX)*FLW - WTMW*UDLA(MP,NPX)
     &    - WTMW*UDS(MP,NPX)/WTMS)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWW_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWE_COMP( N,NB,NQX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, water equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLWE_COMP'
      I = ID(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLWB = XLWB(MB,NB)*RHOLB(MB,NB)
        FLWP = XLW(MB,N)*RHOL(MB,N)
        INDX = 2
        FLW = DIFMN( FLWP,FLWB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
        UDLA(MN,NQX) = 0.D+0
        DO IGC = 1, NGC
         UDLA(MN,NQX) = UDLA(MN,NQX) + UDLC(IGC,MN,NQX)
        ENDDO
        RS(M) = AFX(NQX)*(UL(MN,NQX)*FLW - WTMW*UDLA(MN,NQX)
     &    - WTMW*UDS(MN,NQX)/WTMS)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWE_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWN_COMP( N,NB,NQY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, water equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLWN_COMP'
      J = JD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLWB = XLWB(MB,NB)*RHOLB(MB,NB)
        FLWP = XLW(MB,N)*RHOL(MB,N)
        INDX = 2
        FLW = DIFMN( FLWP,FLWB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
        VDLA(MN,NQY) = 0.D+0
        DO IGC = 1, NGC
         VDLA(MN,NQY) = VDLA(MN,NQY) + VDLC(IGC,MN,NQY)
        ENDDO
        RS(M) = AFY(NQY)*(VL(MN,NQY)*FLW - WTMW*VDLA(MN,NQY)
     &    - WTMW*VDS(MN,NQY)/WTMS)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWN_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWT_COMP( N,NB,NQZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, water equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLWT_COMP'
      K = KD(N)
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLWB = XLWB(MB,NB)*RHOLB(MB,NB)
        FLWP = XLW(MB,N)*RHOL(MB,N)
        INDX = 2
        FLW = DIFMN( FLWP,FLWB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
        WDLA(MN,NQZ) = 0.D+0
        DO IGC = 1, NGC
         WDLA(MN,NQZ) = WDLA(MN,NQZ) + WDLC(IGC,MN,NQZ)
        ENDDO
        RS(M) = AFZ(NQZ)*(WL(MN,NQZ)*FLW - WTMW*WDLA(MN,NQZ)
     &    - WTMW*WDS(MN,NQZ)/WTMS)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQW )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWT_COMP group  ---
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBSB_COMP( N,NB,NPZ )
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
!     Modify the surfactant equations of the Jacobian matrix for
!     boundary conditions on the bottom surface in the aqueous phase.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBSB_COMP'
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFZ(NPZ)*WS(MP,NPZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQS )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBSB_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBSS_COMP( N,NB,NPY )
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
!     Modify the surfactant equations of the Jacobian matrix for
!     boundary conditions on the south surface in the aqueous phase.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBSS_COMP'
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFY(NPY)*VS(MP,NPY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQS )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBSS_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBSW_COMP( N,NB,NPX )
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
!     Modify the surfactant equations of the Jacobian matrix for
!     boundary conditions on the west surface in the aqueous phase.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBSW_COMP'
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFX(NPX)*US(MP,NPX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQS )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBSW_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBSE_COMP( N,NB,NQX )
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
!     Modify the surfactant equations of the Jacobian matrix for
!     boundary conditions on the east surface in the aqueous phase.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBSE_COMP'
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFX(NQX)*US(MN,NQX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQS )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBSE_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBSN_COMP( N,NB,NQY )
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
!     Modify the surfactant equations of the Jacobian matrix for
!     boundary conditions on the north surface in the aqueous phase.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBSN_COMP'
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFY(NQY)*VS(MN,NQY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQS )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBSN_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBST_COMP( N,NB,NQZ )
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
!     Modify the surfactant equations of the Jacobian matrix for
!     boundary conditions on the top surface in the aqueous phase.
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBST_COMP'
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFZ(NQZ)*WS(MN,NQZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQS )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBST_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTB_COMP( N,NB,NPZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (thermal boundary, energy equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBTB_COMP'
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFZ(NPZ)*WQ(MP,NPZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQT )
!
!---  End of JCBTB_COMP group  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTS_COMP( N,NB,NPY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (thermal boundary, energy equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBTS_COMP'
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFY(NPY)*VQ(MP,NPY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQT )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTS_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTW_COMP( N,NB,NPX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (thermal boundary, energy equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBTW_COMP'
      DO 100 M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFX(NPX)*UQ(MP,NPX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQT )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTW_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTE_COMP( N,NB,NQX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (thermal boundary, energy equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBTE_COMP'
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFX(NQX)*UQ(MN,NQX)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQT )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTE_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTN_COMP( N,NB,NQY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (thermal boundary, energy equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBTN_COMP'
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFY(NQY)*VQ(MN,NQY)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQT )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTN_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTT_COMP( N,NB,NQZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (thermal boundary, energy equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBTT_COMP'
      DO 100 M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFZ(NQZ)*WQ(MN,NQZ)
  100 CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBMB( RS,N,IEQT )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTT_COMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBOLGN
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
!     Load the Jacobian matrix for the oil equation with
!     aqueous-phase, gas-phase, and NAPL contributions
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, October 19, 1999.
!     Last Modified by M Oostrom, PNNL, 20 March 2003.
!     Last Modified by MD White, PNNL, 09 January 2004.
!     jcbolgn.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE PORMED
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXN
      USE FLUXD
      USE FDVP
      USE FDVN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STO(LUK+1),ROP(LUK),ROA(LUK,6),FO(LSFV,6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBOLGN'
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO 1000 K = 1,KFLD
        DO 1000 J = 1,JFLD
          DO 1000 I = 1,IFLD
        N = ND(I,J,K)
        IF( IXP(N).EQ.0 ) GOTO 1000
        IZN = IZ(N)
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NPX + 1
        NQY = NPY + IFLD
        NQZ = NPZ + IJFLD
!
!---    First-order, forward-difference, time differential  ---
!
        STO1 = PORD(1,N)*(XLO(1,N)*RHOL(1,N)*SL(1,N) +
     &      XGO(1,N)*RHOG(1,N)*SG(1,N) + RHON(1,N)*SN(1,N)) +
     &      XSO(1,N)*RHOS(IZN)*(1.D+0-PORT(1,N))
        DO 100 M = 1,ISVC+1
          MP = M + 1
          STO0 = PORD(MP,N)*(XLO(MP,N)*RHOL(MP,N)*SL(MP,N) +
     &      XGO(MP,N)*RHOG(MP,N)*SG(MP,N) + RHON(MP,N)*SN(MP,N)) +
     &      XSO(MP,N)*RHOS(IZN)*(1.D+0-PORT(MP,N))
          STO(M) = (STO0-STO1)*DTI*VOL(N)
  100   CONTINUE
!
!---  Initialize surface fluxes  ---
!
        DO 210 MD = 1,6
          DO 200 M = 1,ISVF
            FO(M,MD) = 0.D+0
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
            MB = MADJ(M)
            MP = MNOD(M)
            INDX = 4
            FNOB = RHON(MB,NB)
            FNOP = RHON(MP,N)
            FNO = DIFMN( FNOB,FNOP,DZGF(NB),DZGF(N),WN(1,NPZ),INDX )
            IF( ISLC(12).EQ.1 ) THEN
              FO(M,1) = -AFZ(NPZ)*(WGO(M,NPZ) +
     &          WN(M,NPZ)*FNO + WLO(M,NPZ))
            ELSE
              FLOB = XLO(MB,NB)*RHOL(MB,NB)
              FLOP = XLO(MP,N)*RHOL(MP,N)
              INDX = 2
              FLO = DIFMN( FLOB,FLOP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
              FGOB = XGO(MB,NB)*RHOG(MB,NB)
              FGOP = XGO(MP,N)*RHOG(MP,N)
              INDX = 3
              FGO = DIFMN( FGOB,FGOP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
              FO(M,1) = -AFZ(NPZ)*(WL(M,NPZ)*FLO + WG(M,NPZ)*FGO +
     &          WN(M,NPZ)*FNO + WTMO*(WDGO(M,NPZ)+WDLO(M,NPZ)))
              WLO(M,NPZ) = WL(M,NPZ)*FLO + WTMO*WDLO(M,NPZ)
              WGO(M,NPZ) = WG(M,NPZ)*FGO + WTMO*WDGO(M,NPZ)
            ENDIF
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
            MS = MADJ(M)
            MP = MNOD(M)
            INDX = 4
            FNOS = RHON(MS,NS)
            FNOP = RHON(MP,N)
            FNO = DIFMN( FNOS,FNOP,DYGF(NS),DYGF(N),VN(1,NPY),INDX )
            IF( ISLC(12).EQ.1 ) THEN
              FO(M,2) = -AFY(NPY)*(VGO(M,NPY) +
     &          VN(M,NPY)*FNO + VLO(M,NPY))
            ELSE
              FLOS = XLO(MS,NS)*RHOL(MS,NS)
              FLOP = XLO(MP,N)*RHOL(MP,N)
              INDX = 2
              FLO = DIFMN( FLOS,FLOP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
              FGOS = XGO(MS,NS)*RHOG(MS,NS)
              FGOP = XGO(MP,N)*RHOG(MP,N)
              INDX = 3
              FGO = DIFMN( FGOS,FGOP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
              FO(M,2) = -AFY(NPY)*(VL(M,NPY)*FLO + VG(M,NPY)*FGO +
     &          VN(M,NPY)*FNO + WTMO*(VDGO(M,NPY)+VDLO(M,NPY)))
              VLO(M,NPY) = VL(M,NPY)*FLO + WTMO*VDLO(M,NPY)
              VGO(M,NPY) = VG(M,NPY)*FGO + WTMO*VDGO(M,NPY)
            ENDIF
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
            MW = MADJ(M)
            MP = MNOD(M)
            INDX = 4
            FNOW = RHON(MW,NW)
            FNOP = RHON(MP,N)
            FNO = DIFMN( FNOW,FNOP,DXGF(NW),DXGF(N),UN(1,NPX),INDX )
            IF( ISLC(12).EQ.1 ) THEN
              FO(M,3) = -AFX(NPX)*(UGO(M,NPX) +
     &          UN(M,NPX)*FNO + ULO(M,NPX))
            ELSE
              FLOW = XLO(MW,NW)*RHOL(MW,NW)
              FLOP = XLO(MP,N)*RHOL(MP,N)
              INDX = 2
              FLO = DIFMN( FLOW,FLOP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
              FGOW = XGO(MW,NW)*RHOG(MW,NW)
              FGOP = XGO(MP,N)*RHOG(MP,N)
              INDX = 3
              FGO = DIFMN( FGOW,FGOP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
              FO(M,3) = -AFX(NPX)*(UL(M,NPX)*FLO + UG(M,NPX)*FGO +
     &          UN(M,NPX)*FNO + WTMO*(UDGO(M,NPX)+UDLO(M,NPX)))
              ULO(M,NPX) = UL(M,NPX)*FLO + WTMO*UDLO(M,NPX)
              UGO(M,NPX) = UG(M,NPX)*FGO + WTMO*UDGO(M,NPX)
            ENDIF
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
            ME = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            INDX = 4
            FNOP = RHON(MP,N)
            FNOE = RHON(ME,NE)
            FNO = DIFMN( FNOP,FNOE,DXGF(N),DXGF(NE),UN(1,NQX),INDX )
            IF( ISLC(12).EQ.1 ) THEN
              FO(M,4) = AFX(NQX)*(UGO(MF,NQX) +
     &          UN(MF,NQX)*FNO + ULO(MF,NQX))
            ELSE
              FLOE = XLO(ME,NE)*RHOL(ME,NE)
              FLOP = XLO(MP,N)*RHOL(MP,N)
              INDX = 2
              FLO = DIFMN( FLOP,FLOE,DXGF(N),DXGF(NE),UL(1,NQX),INDX )
              FGOE = XGO(ME,NE)*RHOG(ME,NE)
              FGOP = XGO(MP,N)*RHOG(MP,N)
              INDX = 3
              FGO = DIFMN( FGOP,FGOE,DXGF(N),DXGF(NE),UG(1,NQX),INDX )
              FO(M,4) = AFX(NQX)*(UL(MF,NQX)*FLO + UG(MF,NQX)*FGO +
     &          UN(MF,NQX)*FNO + WTMO*(UDGO(MF,NQX)+UDLO(MF,NQX)))
              ULO(MF,NQX) = UL(MF,NQX)*FLO + WTMO*UDLO(MF,NQX)
              UGO(MF,NQX) = UG(MF,NQX)*FLO + WTMO*UDGO(MF,NQX)
            ENDIF
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
            MN = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            INDX = 4
            FNOP = RHON(MP,N)
            FNON = RHON(MN,NN)
            FNO = DIFMN( FNOP,FNON,DYGF(N),DYGF(NN),VN(1,NQY),INDX )
            IF( ISLC(12).EQ.1 ) THEN
              FO(M,5) = AFY(NQY)*(VGO(MF,NQY) +
     &          VN(MF,NQY)*FNO + VLO(MF,NQY))
            ELSE
              FLON = XLO(MN,NN)*RHOL(MN,NN)
              FLOP = XLO(MP,N)*RHOL(MP,N)
              INDX = 2
              FLO = DIFMN( FLOP,FLON,DYGF(N),DYGF(NN),VL(1,NQY),INDX )
              FGON = XGO(MN,NN)*RHOG(MN,NN)
              FGOP = XGO(MP,N)*RHOG(MP,N)
              INDX = 3
              FGO = DIFMN( FGOP,FGON,DYGF(N),DYGF(NN),VG(1,NQY),INDX )
              FO(M,5) = AFY(NQY)*(VL(MF,NQY)*FLO + VG(MF,NQY)*FGO +
     &          VN(MF,NQY)*FNO + WTMO*(VDGO(MF,NQY)+VDLO(MF,NQY)))
              VLO(MF,NQY) = VL(MF,NQY)*FLO + WTMO*VDLO(MF,NQY)
              VGO(MF,NQY) = VG(MF,NQY)*FGO + WTMO*VDGO(MF,NQY)
            ENDIF
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
            MT = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            INDX = 4
            FNOP = RHON(MP,N)
            FNOT = RHON(MT,NT)
            FNO = DIFMN( FNOP,FNOT,DZGF(N),DZGF(NT),WN(1,NQZ),INDX )
            IF( ISLC(12).EQ.1 ) THEN
              FO(M,6) = AFZ(NQZ)*(WGO(MF,NQZ) +
     &          WN(MF,NQZ)*FNO + WLO(MF,NQZ))
            ELSE
              FLOT = XLO(MT,NT)*RHOL(MT,NT)
              FLOP = XLO(MP,N)*RHOL(MP,N)
              INDX = 2
              FLO = DIFMN( FLOP,FLOT,DZGF(N),DZGF(NT),WL(1,NQZ),INDX )
              FGOT = XGO(MT,NT)*RHOG(MT,NT)
              FGOP = XGO(MP,N)*RHOG(MP,N)
              INDX = 3
              FGO = DIFMN( FGOP,FGOT,DZGF(N),DZGF(NT),WG(1,NQZ),INDX )
              FO(M,6) = AFZ(NQZ)*(WL(MF,NQZ)*FLO + WG(MF,NQZ)*FGO +
     &          WN(MF,NQZ)*FNO + WTMO*(WDGO(MF,NQZ)+WDLO(MF,NQZ)))
              WLO(MF,NQZ) = WL(MF,NQZ)*FLO + WTMO*WDLO(MF,NQZ)
              WGO(MF,NQZ) = WG(MF,NQZ)*FLO + WTMO*WDGO(MF,NQZ)
            ENDIF
  800     CONTINUE
        ENDIF
  810   CONTINUE
!
!---  Compute VOC equation residuals  ---
!
        ROS = STO(1) - SRCO(2,N)
        DO 900 MD = 1,6
          ROS = ROS + FO(1,MD)
  900   CONTINUE
        DO 920 M = 1,ISVC
          ROP(M) = STO(M+1) - SRCO(M+2,N)
          MM = 2*M
          DO 910 MD = 1,6
            ROP(M) = ROP(M) + FO(MM,MD)
  910     CONTINUE
  920   CONTINUE
        DO 940 M = 1,ISVC
          MM = 2*M + 1
          DO 930 MD = 1,6
            ROA(M,MD) = ROS - FO(1,MD) + FO(MM,MD)
  930     CONTINUE
  940   CONTINUE
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBBM( ROS,ROP,ROA,N,I,J,K,IEQO )
!
!---  Continue to Next Node  ---
!
 1000 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBOLGN group
!
      RETURN
      END

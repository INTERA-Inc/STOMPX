!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBWTM
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
!     Load the Jacobian matrix for the water equation for
!     an aqueous-phase/NAPL ternary mixture system
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, October 19, 1999.
!     jcbwtm.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXN
      USE FDVP
      USE FDVA
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
      SUB_LOG(ISUB_LOG) = '/JCBWTM'
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
!
!---    First-order, forward-difference, time differential  ---
!
        STW1 = PORD(1,N)*(XLW(1,N)*RHOL(1,N)*SL(1,N) +
     &    XNW(1,N)*RHON(1,N)*SN(1,N))
        DO 100 M = 1,ISVC+1
          MP = M + 1
          STW0 = PORD(MP,N)*(XLW(MP,N)*RHOL(MP,N)*SL(MP,N) +
     &      XNW(MP,N)*RHON(MP,N)*SN(MP,N))
          STWX(M) = (STW0-STW1)*DTI*VOL(N)
  100   CONTINUE
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
          IF( IXP(NB).EQ.0 .OR. INBS(1,N).GT.0 ) GOTO 310
          DO 300 M = 1,ISVF
            MB = MADJ(M)
            MP = MNOD(M)
            FLWB = XLW(MB,NB)*RHOL(MB,NB)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWB,FLWP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
            FNWB = XNW(MB,NB)*RHON(MB,NB)
            FNWP = XNW(MP,N)*RHON(MP,N)
            INDX = 4
            FNW = DIFMN( FNWB,FNWP,DZGF(NB),DZGF(N),WN(1,NPZ),INDX )
            FW(M,1) = -AFZ(NPZ)*(WL(M,NPZ)*FLW + WN(M,NPZ)*FNW +
     &        (WDGW(M,NPZ)-WDLO(M,NPZ)-WDLA(M,NPZ))*WTMW)
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
            FLWS = XLW(MS,NS)*RHOL(MS,NS)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWS,FLWP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
            FNWS = XNW(MS,NS)*RHON(MS,NS)
            FNWP = XNW(MP,N)*RHON(MP,N)
            INDX = 4
            FNW = DIFMN( FNWS,FNWP,DYGF(NS),DYGF(N),VN(1,NPY),INDX )
            FW(M,2) = -AFY(NPY)*(VL(M,NPY)*FLW + VN(M,NPY)*FNW +
     &        (VDGW(M,NPY)-VDLO(M,NPY)-VDLA(M,NPY))*WTMW)
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
            FLWW = XLW(MW,NW)*RHOL(MW,NW)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWW,FLWP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            FNWW = XNW(MW,NW)*RHON(MW,NW)
            FNWP = XNW(MP,N)*RHON(MP,N)
            INDX = 4
            FNW = DIFMN( FNWW,FNWP,DXGF(NW),DXGF(N),UN(1,NPX),INDX )
            FW(M,3) = -AFX(NPX)*(UL(M,NPX)*FLW + UN(M,NPX)*FNW +
     &        (UDGW(M,NPX)-UDLO(M,NPX)-UDLA(M,NPX))*WTMW)
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
            FLWE = XLW(ME,NE)*RHOL(ME,NE)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWP,FLWE,DXGF(N),DXGF(NE),UL(1,NQX),INDX )
            FNWE = XNW(ME,NE)*RHON(ME,NE)
            FNWP = XNW(MP,N)*RHON(MP,N)
            INDX = 4
            FNW = DIFMN( FNWP,FNWE,DXGF(N),DXGF(NE),UN(1,NQX),INDX )
            FW(M,4) = AFX(NQX)*(UL(MF,NQX)*FLW + UN(MF,NQX)*FNW +
     &        (UDGW(MF,NQX)-UDLO(MF,NQX)-UDLA(MF,NQX))*WTMW)
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
            FLWN = XLW(MN,NN)*RHOL(MN,NN)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWP,FLWN,DYGF(N),DYGF(NN),VL(1,NQY),INDX )
            FNWN = XNW(MN,NN)*RHON(MN,NN)
            FNWP = XNW(MP,N)*RHON(MP,N)
            INDX = 4
            FNW = DIFMN( FNWP,FNWN,DYGF(N),DYGF(NN),VN(1,NQY),INDX )
            FW(M,5) = AFY(NQY)*(VL(MF,NQY)*FLW + VN(MF,NQY)*FNW +
     &        (VDGW(MF,NQY)-VDLO(MF,NQY)-VDLA(MF,NQY))*WTMW)
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
            FLWT = XLW(MT,NT)*RHOL(MT,NT)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWP,FLWT,DZGF(N),DZGF(NT),WL(1,NQZ),INDX )
            FW(M,6) = AFZ(NQZ)*(WL(MF,NQZ)*FLW - WTMW*WDLO(MF,NQZ))
            FNWT = XNW(MT,NT)*RHON(MT,NT)
            FNWP = XNW(MP,N)*RHON(MP,N)
            INDX = 4
            FNW = DIFMN( FNWP,FNWT,DZGF(N),DZGF(NT),WN(1,NQZ),INDX )
            FW(M,6) = AFZ(NQZ)*(WL(MF,NQZ)*FLW + WN(MF,NQZ)*FNW +
     &        (WDGW(MF,NQZ)-WDLO(MF,NQZ)-WDLA(MF,NQZ))*WTMW)
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
!---  End of JCBWTM group
!
      RETURN
      END

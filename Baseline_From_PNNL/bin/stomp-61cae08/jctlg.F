!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCTLG
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
!     Load the Jacobian matrix for the VOC transport equation
!     aqueous-phase contributions, gas-phase contributions,
!     no NAPL contributions (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1996.
!     Last Modified by MD White on September 5, 1996.
!     jctlg.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVG
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCTLG'
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO 1000 K = 1,KFLD
        DO 1000 J = 1,JFLD
          DO 1000 I = 1,IFLD
        N = ND(I,J,K)
        IZN = IZ(N)
        IF( IXP(N).EQ.0 ) GOTO 1000
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NPX + 1
        NQY = NPY + IFLD
        NQZ = NPZ + IJFLD
!
!---  First-order, forward-difference, time differential  ---
!
        SC = VOL(N)*DTI
        MP = IXP(N)
        MA = 0
        IF( ILES.EQ.1 ) THEN
          MCD = MP
          MRD = MDT
          ALU(MRD,MCD) = SC
        ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
          MCD = KLU(MP,MA+1)
          DLU(MCD) = SC
        ENDIF
        MA = MA + 1
!
!---  Compute diffusion coefficients  ---
!
        DGP = TORG(2,N)*PORD(2,N)*SG(2,N)*DFGO(2,N)
        DLP = TORL(2,N)*PORD(2,N)*SL(2,N)*DFLO(2,N)
!
!---  Compute phase fractions  ---
!
        FCL = 0.D+0
        FCG = 0.D+0
        IF( SL(2,N).GT.SMALL ) FCL = YL(N,1)/(SL(2,N)*PORD(2,N))
        IF( SG(2,N).GT.SMALL ) FCG = YG(N,1)/(SG(2,N)*PORD(2,N))
!
!---  Bottom ---
!
        IF( K.NE.1 ) THEN
          M = 1
          NB = N-IJFLD
          IZNB = IZ(NB)
          IF( IXP(NB).EQ.0 .OR. INBS(1,N).GT.0 ) GOTO 310
          FCLB = 0.D+0
          FCGB = 0.D+0
          IF( SL(2,NB).GT.SMALL )FCLB = YL(NB,1)/(SL(2,NB)*PORD(2,NB))
          IF( SG(2,NB).GT.SMALL )FCGB = YG(NB,1)/(SG(2,NB)*PORD(2,NB))
          DGB = TORG(2,NB)*PORD(2,NB)*SG(2,NB)*DFGO(2,NB)
          IF( IDISP.EQ.1 ) THEN
            CALL ADVB( PORD,SG,UG,VG,WG,UGBX,VGBX,WGBX,N,M )
            INDX = 17
            DSPLB = DIFMN(DISPL(IZNB),DISPL(IZN),DZGF(NB),DZGF(N),
     &        WGBX,INDX)
            DSPTB = DIFMN(DISPT(IZNB),DISPT(IZN),DZGF(NB),DZGF(N),
     &        WGBX,INDX)
            CALL SHDP( WGBX,UGBX,VGBX,DSPLB,DSPTB,DPGB )
          ELSE
            DPGB = 0.D+0
          ENDIF
          DCG = C(NB,1)*FCGB - C(N,1)*FCG
          INDX = 16
          DGZ = (DIFMN( DGB,DGP,DZGF(NB),DZGF(N),DCG,INDX)+DPGB)/
     &      DZGP(NPZ)
          AG = MAX( WG(1,NPZ),ZERO ) +
     &      DGZ*MAX((ONE-(TENTH*ABS(WG(1,NPZ))/(DGZ+SMALL)))**5,ZERO )
          DLB = TORL(2,NB)*PORD(2,NB)*SL(2,NB)*DFLO(2,NB)
          IF( IDISP.EQ.1 ) THEN
            CALL ADVB( PORD,SL,UL,VL,WL,ULBX,VLBX,WLBX,N,M )
            INDX = 17
            DSPLB = DIFMN(DISPL(IZNB),DISPL(IZN),DZGF(NB),DZGF(N),
     &        WLBX,INDX)
            DSPTB = DIFMN(DISPT(IZNB),DISPT(IZN),DZGF(NB),DZGF(N),
     &        WLBX,INDX)
            CALL SHDP( WLBX,ULBX,VLBX,DSPLB,DSPTB,DPLB )
          ELSE
            DPLB = 0.D+0
          ENDIF
          DCL = C(NB,1)*FCLB - C(N,1)*FCL
          INDX = 16
          DLZ = (DIFMN( DLB,DLP,DZGF(NB),DZGF(N),DCL,INDX)+DPLB)/
     &      DZGP(NPZ)
          AL = MAX( WL(1,NPZ),ZERO ) +
     &      DLZ*MAX((ONE-(TENTH*ABS(WL(1,NPZ))/(DLZ+SMALL)))**5,ZERO )
          AP = AFZ(NPZ)*(((AL-WL(1,NPZ))*FCL + (AG-WG(1,NPZ))*FCG)
          AB = AFZ(NPZ)*(AL*FCLB + AG*FCGB)
          MB = IXP(NB)
          IF( ILES.EQ.1 ) THEN
            ALU(MRD,MCD) = ALU(MRD,MCD) + AP
            MCOL = MB
            MROW = MP-MB+MDT
            ALU(MROW,MCOL) = -AB
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MCOL = KLU(MP,MA+1)
            DLU(MCOL) = -AB
            MA = MA + 1
          ENDIF
        ENDIF
  310   CONTINUE
!
!---  South ---
!
        IF( J.NE.1 ) THEN
          M = 1
          NS = N-IFLD
          IZNS = IZ(NS)
          IF( IXP(NS).EQ.0 .OR. INBS(2,N).GT.0 ) GOTO 410
          FCLS = 0.D+0
          FCGS = 0.D+0
          IF( SL(2,NS).GT.SMALL )FCLS = YL(NS,1)/(SL(2,NS)*PORD(2,NS))
          IF( SG(2,NS).GT.SMALL )FCGS = YG(NS,1)/(SG(2,NS)*PORD(2,NS))
          DGS = TORG(2,NS)*PORD(2,NS)*SG(2,NS)*DFGO(2,NS)
          IF( IDISP.EQ.1 ) THEN
            CALL ADVS( PORD,SG,UG,VG,WG,UGSX,VGSX,WGSX,N,M )
            INDX = 17
            DSPLS = DIFMN(DISPL(IZNS),DISPL(IZN),DYGF(NS),DYGF(N),
     &        VGSX,INDX)
            DSPTS = DIFMN(DISPT(IZNS),DISPT(IZN),DYGF(NS),DYGF(N),
     &        VGSX,INDX)
            CALL SHDP( VGSX,WGSX,UGSX,DSPLS,DSPTS,DPGS )
          ELSE
            DPGS = 0.D+0
          ENDIF
          DCG = C(NS,1)*FCGS - C(N,1)*FCG
          INDX = 16
          DGY = (DIFMN( DGS,DGP,DYGF(NS),DYGF(N),DCG,INDX)+DPGS)
     &      DYGP(NPY)/RP(I)
          AG = MAX( VG(1,NPY),ZERO ) +
     &      DGY*MAX((ONE-(TENTH*ABS(VG(1,NPY))/(DGY+SMALL)))**5,ZERO )
          DLS = TORL(2,NS)*PORD(2,NS)*SL(2,NS)*DFLO(2,NS)
          IF( IDISP.EQ.1 ) THEN
            CALL ADVS( PORD,SL,UL,VL,WL,ULSX,VLSX,WLSX,N,M )
            INDX = 17
            DSPLS = DIFMN(DISPL(IZNS),DISPL(IZN),DYGF(NS),DYGF(N),
     &        VLSX,INDX)
            DSPTS = DIFMN(DISPT(IZNS),DISPT(IZN),DYGF(NS),DYGF(N),
     &        VLSX,INDX)
            CALL SHDP( VLSX,WLSX,ULSX,DSPLS,DSPTS,DPLS )
          ELSE
            DPLS = 0.D+0
          ENDIF
          DCL = C(NS,1)*FCLS - C(N,1)*FCL
          INDX = 16
          DLY = (DIFMN( DLS,DLP,DYGF(NS),DYGF(N),DCL,INDX)+DPLS)
     &      DYGP(NPY)/RP(I)
          AL = MAX( VL(1,NPY),ZERO ) +
     &      DLY*MAX((ONE-(TENTH*ABS(VL(1,NPY))/(DLY+SMALL)))**5,ZERO )
          AP = AFY(NPY)*((AL-VL(1,NPY))*FCL + (AG-VG(1,NPY))*FCG)
          AS = AFY(NPY)*(AL*FCLS + AG*FCGS)
          MS = IXP(NS)
          IF( ILES.EQ.1 ) THEN
            ALU(MRD,MCD) = ALU(MRD,MCD) + AP
            MCOL = MS
            MROW = MP-MS+MDT
            ALU(MROW,MCOL) = -AS
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MCOL = KLU(MP,MA+1)
            DLU(MCOL) = -AS
            MA = MA + 1
          ENDIF
        ENDIF
  410   CONTINUE
!
!---  West ---
!
        IF( I.NE.1 ) THEN
          M = 1
          NW = N-1
          IZNW = IZ(NW)
          IF( IXP(NW).EQ.0 .OR. INBS(3,N).GT.0 ) GOTO 510
          FCLW = 0.D+0
          FCGW = 0.D+0
          IF( SL(2,NW).GT.SMALL )FCLW = YL(NW,1)/(SL(2,NW)*PORD(2,NW))
          IF( SG(2,NW).GT.SMALL )FCGW = YG(NW,1)/(SG(2,NW)*PORD(2,NW))
          DGW = TORG(2,NW)*PORD(2,NW)*SG(2,NW)*DFGO(2,NW)
          IF( IDISP.EQ.1 ) THEN
            CALL ADVW( PORD,SG,UG,VG,WG,UGX,VGX,WGX,N,M )
            INDX = 17
            DSPLW = DIFMN(DISPL(IZNW),DISPL(IZN),DXGF(NW),DXGF(N),
     &        UGX,INDX)
            DSPTW = DIFMN(DISPT(IZNW),DISPT(IZN),DXGF(NW),DXGF(N),
     &        UGX,INDX)
            CALL SHDP( UGX,VGX,WGX,DSPLW,DSPTW,DPGW )
          ELSE
            DPGW = 0.D+0
          ENDIF
          DCG = C(NW,1)*FCGW - C(N,1)*FCG
          INDX = 16
          DGX = (DIFMN( DGW,DGP,DXGF(NW),DXGF(N),DCG,INDX)+DPGW)/
     &      DXGP(NPX)
          AG = MAX( UG(1,NPX),ZERO ) +
     &      DGX*MAX((ONE-(TENTH*ABS(UG(1,NPX))/(DGX+SMALL)))**5,ZERO )
          DLW = TORL(2,NW)*PORD(2,NW)*SL(2,NW)*DFLO(2,NW)
          IF( IDISP.EQ.1 ) THEN
            CALL ADVW( PORD,SL,UL,VL,WL,ULX,VLX,WLX,N,M )
            INDX = 17
            DSPLW = DIFMN(DISPL(IZNW),DISPL(IZN),DXGF(NW),DXGF(N),
     &        ULX,INDX)
            DSPTW = DIFMN(DISPT(IZNW),DISPT(IZN),DXGF(NW),DXGF(N),
     &        ULX,INDX)
            CALL SHDP( ULX,VLX,WLX,DSPLW,DSPTW,DPLW )
          ELSE
            DPLW = 0.D+0
          ENDIF
          DCL = C(NW,1)*FCLW - C(N,1)*FCL
          INDX = 16
          DLX = (DIFMN( DLW,DLP,DXGF(NW),DXGF(N),DCL,INDX)+DPLW)/
     &      DXGP(NPX)
          AL = MAX( UL(1,NPX),ZERO ) +
     &      DLX*MAX((ONE-(TENTH*ABS(UL(1,NPX))/(DLX+SMALL)))**5,ZERO )
          AP = AFX(NPX)*((AL-UL(1,NPX))*FCL  + (AG-UG(1,NPX))*FCG)
          AW = AFX(NPX)*(AL*FCLW + AG*FCGW)
          MW = IXP(NW)
          IF( ILES.EQ.1 ) THEN
            ALU(MRD,MCD) = ALU(MRD,MCD) + AP
            MCOL = MW
            MROW = MP-MW+MDT
            ALU(MROW,MCOL) = -AW
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MCOL = KLU(MP,MA+1)
            DLU(MCOL) = -AW
            MA = MA + 1
          ENDIF
        ENDIF
  510   CONTINUE
!
!---  East ---
!
        IF( I.NE.IFLD ) THEN
          M = 1
          NE = N+1
          IZNE = IZ(NE)
          IF( IXP(NE).EQ.0 .OR. INBS(4,N).GT.0 ) GOTO 610
          FCLE = 0.D+0
          FCGE = 0.D+0
          IF( SL(2,NE).GT.SMALL )FCLE = YL(NE,1)/(SL(2,NE)*PORD(2,NE))
          IF( SG(2,NE).GT.SMALL )FCGE = YG(NE,1)/(SG(2,NE)*PORD(2,NE))
          DGE = TORG(2,NE)*PORD(2,NE)*SG(2,NE)*DFGO(2,NE)
          IF( IDISP.EQ.1 ) THEN
            CALL ADVE( PORD,SG,UG,VG,WG,UGEX,VGEX,WGEX,N,M )
            INDX = 17
            DSPLE = DIFMN(DISPL(IZN),DISPL(IZNE),DXGF(N),DXGF(NE),
     &        UGEX,INDX)
            DSPTE = DIFMN(DISPT(IZN),DISPT(IZNE),DXGF(N),DXGF(NE),
     &        UGEX,INDX)
            CALL SHDP( UGEX,VGEX,WGEX,DSPLE,DSPTE,DPGE )
          ELSE
            DPGE = 0.D+0
          ENDIF
          DCG = C(N,1)*FCG - C(NE,1)*FCGE
          INDX = 16
          DGX = (DIFMN(DGP,DGE,DXGF(N),DXGF(NE),DCG,INDX)
     &      +DPGE)/DXGP(NQX)
          AG = MAX( -UG(1,NQX),ZERO ) +
     &      DGX*MAX((ONE-(TENTH*ABS(UG(1,NQX))/(DGX+SMALL)))**5,ZERO )
          DLE = TORL(2,NE)*PORD(2,NE)*SL(2,NE)*DFLO(2,NE)
          IF( IDISP.EQ.1 ) THEN
            CALL ADVE( PORD,SL,UL,VL,WL,ULEX,VLEX,WLEX,N,M )
            INDX = 17
            DSPLE = DIFMN(DISPL(IZN),DISPL(IZNE),DXGF(N),DXGF(NE),
     &        ULEX,INDX)
            DSPTE = DIFMN(DISPT(IZN),DISPT(IZNE),DXGF(N),DXGF(NE),
     &        ULEX,INDX)
            CALL SHDP( ULEX,VLEX,WLEX,DSPLE,DSPTE,DPLE )
          ELSE
            DPLE = 0.D+0
          ENDIF
          DCL = C(N,1)*FCL - C(NE,1)*FCLE
          INDX = 16
          DLX = (DIFMN(DLP,DLE,DXGF(N),DXGF(NE),DCL,INDX)
     &      +DPLE)/DXGP(NQX)
          AL = MAX( -UL(1,NQX),ZERO ) +
     &      DLX*MAX((ONE-(TENTH*ABS(UL(1,NQX))/(DLX+SMALL)))**5,ZERO )
          AP = AFX(NQX)*((AL+UL(1,NQX))*FCL + (AG+UG(1,NQX))*FCG)
          AE = AFX(NQX)*(AL*FCLE + AG*FCGE)
          ME = IXP(NE)
          IF( ILES.EQ.1 ) THEN
            ALU(MRD,MCD) = ALU(MRD,MCD) + AP
            MCOL = ME
            MROW = MP-ME+MDT
            ALU(MROW,MCOL) = -AE
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MCOL = KLU(MP,MA+1)
            DLU(MCOL) = -AE
            MA = MA + 1
          ENDIF
        ENDIF
  610   CONTINUE
!
!---  North ---
!
        IF( J.NE.JFLD ) THEN
          M = 1
          NN = N+IFLD
          IZNN = IZ(NN)
          IF( IXP(NN).EQ.0 .OR. INBS(5,N).GT.0 ) GOTO 710
          FCLN = 0.D+0
          FCGN = 0.D+0
          IF( SL(2,NN).GT.SMALL )FCLN = YL(NN,1)/(SL(2,NN)*PORD(2,NN))
          IF( SG(2,NN).GT.SMALL )FCGN = YG(NN,1)/(SG(2,NN)*PORD(2,NN))
          DGN = TORG(2,NN)*PORD(2,NN)*SG(2,NN)*DFGO(2,NN)
          IF( IDISP.EQ.1 ) THEN
            CALL ADVN( PORD,SG,UG,VG,WG,UGNX,VGNX,WGNX,N,M )
            INDX = 17
            DSPLN = DIFMN(DISPL(IZN),DISPL(IZNN),DYGF(N),DYGF(NN),
     &        VGNX,INDX)
            DSPTN = DIFMN(DISPT(IZN),DISPT(IZNN),DYGF(N),DYGF(NN),
     &        VGNX,INDX)
            CALL SHDP( VGNX,WGNX,UGNX,DSPLN,DSPTN,DPGN )
          ELSE
            DPGN = 0.D+0
          ENDIF
          DCG = C(N,1)*FCG - C(NN,1)*FCGN
          INDX = 16
          DGY = (DIFMN(DGP,DGN,DYGF(N),DYGF(NN),DCG,INDX)
     &      +DPGN)/DYGP(NQY)/RP(I)
          AG = MAX( -VG(1,NQY),ZERO ) +
     &      DGY*MAX((ONE-(TENTH*ABS(VG(1,NQY))/(DGY+SMALL)))**5,ZERO )
          DLN = TORL(2,NN)*PORD(2,NN)*SL(2,NN)*DFLO(2,NN)
          IF( IDISP.EQ.1 ) THEN
            CALL ADVN( PORD,SL,UL,VL,WL,ULNX,VLNX,WLNX,N,M )
            INDX = 17
            DSPLN = DIFMN(DISPL(IZN),DISPL(IZNN),DYGF(N),DYGF(NN),
     &        VLNX,INDX)
            DSPTN = DIFMN(DISPT(IZN),DISPT(IZNN),DYGF(N),DYGF(NN),
     &        VLNX,INDX)
            CALL SHDP( VLNX,WLNX,ULNX,DSPLN,DSPTN,DPLN )
          ELSE
            DPLN = 0.D+0
          ENDIF
          DCL = C(N,1)*FCL - C(NN,1)*FCLN
          INDX = 16
          DLY = (DIFMN(DLP,DLN,DYGF(N),DYGF(NN),DCL,INDX)
     &      +DPLN)/DYGP(NQY)/RP(I)
          AL = MAX( -VL(1,NQY),ZERO ) +
     &      DLY*MAX((ONE-(TENTH*ABS(VL(1,NQY))/(DLY+SMALL)))**5,ZERO )
          AP = AFY(NQY)*((AL+VL(1,NQY))*FCL + (AG+VG(1,NQY))*FCG)
          AN = AFY(NQY)*(AL*FCLN + AG*FCGN)
          MN = IXP(NN)
          IF( ILES.EQ.1 ) THEN
            ALU(MRD,MCD) = ALU(MRD,MCD) + AP
            MCOL = MN
            MROW = MP-MN+MDT
            ALU(MROW,MCOL) = -AN
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MCOL = KLU(MP,MA+1)
            DLU(MCOL) = -AN
            MA = MA + 1
          ENDIF
        ENDIF
  710   CONTINUE
!
!---  Top ---
!
        IF( K.NE.KFLD ) THEN
          M = 1
          NT = N+IJFLD
          IZNT = IZ(NT)
          IF( IXP(NT).EQ.0 .OR. INBS(6,N).GT.0 ) GOTO 810
          FCLT = 0.D+0
          FCGT = 0.D+0
          IF( SL(2,NT).GT.SMALL )FCLT = YL(NT,1)/(SL(2,NT)*PORD(2,NT))
          IF( SG(2,NT).GT.SMALL )FCGT = YG(NT,1)/(SG(2,NT)*PORD(2,NT))
          DGT = TORG(2,NT)*PORD(2,NT)*SG(2,NT)*DFGO(2,NT)
          IF( IDISP.EQ.1 ) THEN
            CALL ADVT( PORD,SG,UG,VG,WG,UGTX,VGTX,WGTX,N,M )
            INDX = 17
            DSPLT = DIFMN(DISPL(IZN),DISPL(IZNT),DZGF(N),DZGF(NT),
     &        WGTX,INDX)
            DSPTT = DIFMN(DISPT(IZN),DISPT(IZNT),DZGF(N),DZGF(NT),
     &        WGTX,INDX)
            CALL SHDP( WGTX,UGTX,VGTX,DSPLT,DSPTT,DPGT )
          ELSE
            DPGT = 0.D+0
          ENDIF
          DCG = C(N,1)*FCG - C(NT,1)*FCGT
          INDX = 16
          DGZ = (DIFMN( DGP,DGT,DZGF(N),DZGF(NT),DCG,INDX)
     &      +DPGT)/DZGP(NQZ)
          AG = MAX( -WG(1,NQZ),ZERO ) +
     &      DGZ*MAX((ONE-(TENTH*ABS(WG(1,NQZ))/(DGZ+SMALL)))**5,ZERO )
          DLT = TORL(2,NT)*PORD(2,NT)*SL(2,NT)*DFLO(2,NT)
          IF( IDISP.EQ.1 ) THEN
            CALL ADVT( PORD,SL,UL,VL,WL,ULTX,VLTX,WLTX,N,M )
            INDX = 17
            DSPLT = DIFMN(DISPL(IZN),DISPL(IZNT),DZGF(N),DZGF(NT),
     &        WLTX,INDX)
            DSPTT = DIFMN(DISPT(IZN),DISPT(IZNT),DZGF(N),DZGF(NT),
     &        WLTX,INDX)
            CALL SHDP( WLTX,ULTX,VLTX,DSPLT,DSPTT,DPLT )
          ELSE
            DPLT = 0.D+0
          ENDIF
          DCL = C(N,1)*FCL - C(NT,1)*FCLT
          INDX = 16
          DLZ = (DIFMN( DLP,DLT,DZGF(N),DZGF(NT),DCL,INDX)
     &      +DPLT)/DZGP(NQZ)
          AL = MAX( -WL(1,NQZ),ZERO ) +
     &      DLZ*MAX((ONE-(TENTH*ABS(WL(1,NQZ))/(DLZ+SMALL)))**5,ZERO )
          AP = AFZ(NQZ)*((AL+WL(1,NQZ))*FCL  + (AG+WG(1,NQZ))*FCG)
          AT = AFZ(NQZ)*(AL*FCLT + AG*FCGT)
          MT = IXP(NT)
          IF( ILES.EQ.1 ) THEN
            ALU(MRD,MCD) = ALU(MRD,MCD) + AP
            MCOL = MT
            MROW = MP-MT+MDT
            ALU(MROW,MCOL) = -AT
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MCOL = KLU(MP,MA+1)
            DLU(MCOL) = -AT
            MA = MA + 1
          ENDIF
        ENDIF
  810   CONTINUE
!
!---  Solution Vector  ---
!
        BLU(MP) = C(N,1)*SC
 1000 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCTLG group
!
      RETURN
      END

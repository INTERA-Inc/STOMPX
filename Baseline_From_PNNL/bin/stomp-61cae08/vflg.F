!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VFLG
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
!     Compute VOC transport fluxes, gas-phase diffusion,
!     aqueous-phase advection, gas-phase advection,
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, April 8, 1994.
!     vflg.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
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
      SUB_LOG(ISUB_LOG) = '/VFLG'
      M = 1
!
!---  X-direction VOC transport flux, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
      DO 100 K = 1,KFLD
        DO 100 J = 1,JFLD
          DO 100 I = 2,IFLD
        N = ND(I,J,K)
        NW = N-1
        IZN = IZ(N)
        IZNW = IZ(NW)
        IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &    INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 100
        NPX = NSX(N)
        FCL = 0.D+0
        FCG = 0.D+0
        IF( SL(2,N).GT.SMALL ) FCL = YL(N,1)/(SL(2,N)*PORD(2,N))
        IF( SG(2,N).GT.SMALL ) FCG = YG(N,1)/(SG(2,N)*PORD(2,N))
        FCLW = 0.D+0
        FCGW = 0.D+0
        IF( SL(2,NW).GT.SMALL ) FCLW = YL(NW,1)/(SL(2,NW)*PORD(2,NW))
        IF( SG(2,NW).GT.SMALL ) FCGW = YG(NW,1)/(SG(2,NW)*PORD(2,NW))
        DGP = TORG(2,N)*PORD(2,N)*SG(2,N)*DFGO(2,N)
        DGW = TORG(2,NW)*PORD(2,NW)*SG(2,NW)*DFGO(2,NW)
        IF( IDISP.EQ.1 ) THEN
          CALL ADVW( PORD,SG,UG,VG,WG,UGX,VGX,WGX,N,M )
          INDX = 17
          DSPLW = DIFMN(DISPL(IZNW),DISPL(IZN),DXGF(NW),DXGF(N),
     &      UGX,INDX)
          DSPTW = DIFMN(DISPT(IZNW),DISPT(IZN),DXGF(NW),DXGF(N),
     &      UGX,INDX)
          CALL SHDP( UGX,VGX,WGX,DSPLW,DSPTW,DPGW )
        ELSE
          DPGW = 0.D+0
        ENDIF
        DCG = C(NW,1)*FCGW - C(N,1)*FCG
        INDX = 16
        DGX = (DIFMN( DGW,DGP,DXGF(NW),DXGF(N),DCG,INDX)+DPGW)/DXGP(NPX)
        AG = MAX( UG(1,NPX),ZERO ) +
     &    DGX*MAX( (ONE-(TENTH*ABS(UG(1,NPX))/(DGX+SMALL)))**5,ZERO )
        AGP = MAX( -UG(1,NPX),ZERO ) +
     &    DGX*MAX( (ONE-(TENTH*ABS(UG(1,NPX))/(DGX+SMALL)))**5,ZERO )
        DLP = TORL(2,N)*PORD(2,N)*SL(2,N)*DFLO(2,N)
        DLW = TORL(2,NW)*PORD(2,NW)*SL(2,NW)*DFLO(2,NW)
        IF( IDISP.EQ.1 ) THEN
          CALL ADVW( PORD,SL,UL,VL,WL,ULX,VLX,WLX,N,M )
          INDX = 17
          DSPLW = DIFMN(DISPL(IZNW),DISPL(IZN),DXGF(NW),DXGF(N),
     &      ULX,INDX)
          DSPTW = DIFMN(DISPT(IZNW),DISPT(IZN),DXGF(NW),DXGF(N),
     &      ULX,INDX)
          CALL SHDP( ULX,VLX,WLX,DSPLW,DSPTW,DPLW )
        ELSE
          DPLW = 0.D+0
        ENDIF
        DCL = C(NW,1)*FCLW - C(N,1)*FCL
        INDX = 16
        DLX = (DIFMN( DLW,DLP,DXGF(NW),DXGF(N),DCL,INDX)+DPLW)/DXGP(NPX)
        AL = MAX( UL(1,NPX),ZERO ) +
     &    DLX*MAX( (ONE-(TENTH*ABS(UL(1,NPX))/(DLX+SMALL)))**5,ZERO )
        ALP = MAX( -UL(1,NPX),ZERO ) +
     &    DLX*MAX( (ONE-(TENTH*ABS(UL(1,NPX))/(DLX+SMALL)))**5,ZERO )
        UC(NPX,1) = AFX(NPX)*
     &    ( C(NW,1)*(AL*FCLW + AG*FCGW) - C(N,1)*(ALP*FCL + AGP*FCG) )
  100 CONTINUE
      ENDIF
!
!---  Y-direction VOC transport flux, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
      DO 200 I = 1,IFLD
        DO 200 K = 1,KFLD
          DO 200 J = 2,JFLD
        N = ND(I,J,K)
        NS = N-IFLD
        IZN = IZ(N)
        IZNS = IZ(NS)
        IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &    INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 200
        NPY = NSY(N)
        FCL = 0.D+0
        FCG = 0.D+0
        IF( SL(2,N).GT.SMALL ) FCL = YL(N,1)/(SL(2,N)*PORD(2,N))
        IF( SG(2,N).GT.SMALL ) FCG = YG(N,1)/(SG(2,N)*PORD(2,N))
        FCLS = 0.D+0
        FCGS = 0.D+0
        IF( SL(2,NS).GT.SMALL ) FCLS = YL(NS,1)/(SL(2,NS)*PORD(2,NS))
        IF( SG(2,NS).GT.SMALL ) FCGS = YG(NS,1)/(SG(2,NS)*PORD(2,NS))
        DGP = TORG(2,N)*PORD(2,N)*SG(2,N)*DFGO(2,N)
        DGS = TORG(2,NS)*PORD(2,NS)*SG(2,NS)*DFGO(2,NS)
        IF( IDISP.EQ.1 ) THEN
          CALL ADVS( PORD,SG,UG,VG,WG,UGSX,VGSX,WGSX,N,M )
          INDX = 17
          DSPLS = DIFMN(DISPL(IZNS),DISPL(IZN),DYGF(NS),DYGF(N),
     &      VGSX,INDX)
          DSPTS = DIFMN(DISPT(IZNS),DISPT(IZN),DYGF(NS),DYGF(N),
     &      VGSX,INDX)
          CALL SHDP( VGSX,WGSX,UGSX,DSPLS,DSPTS,DPGS )
        ELSE
          DPGS = 0.D+0
        ENDIF
        DCG = C(NS,1)*FCGS - C(N,1)*FCG
        INDX = 16
        DGY = (DIFMN( DGS,DGP,DYGF(NS),DYGF(N),DCG,INDX)+DPGS)
     &      /DYGP(NPY)/RP(I)
        AG = MAX( VG(1,NPY),ZERO ) +
     &    DGY*MAX( (ONE-(TENTH*ABS(VG(1,NPY))/(DGY+SMALL)))**5,ZERO )
        AGP = MAX( -VG(1,NPY),ZERO ) +
     &    DGY*MAX( (ONE-(TENTH*ABS(VG(1,NPY))/(DGY+SMALL)))**5,ZERO )
        DLP = TORL(2,N)*PORD(2,N)*SL(2,N)*DFLO(2,N)
        DLS = TORL(2,NS)*PORD(2,NS)*SL(2,NS)*DFLO(2,NS)
        IF( IDISP.EQ.1 ) THEN
          CALL ADVS( PORD,SL,UL,VL,WL,ULSX,VLSX,WLSX,N,M )
          INDX = 17
          DSPLS = DIFMN(DISPL(IZNS),DISPL(IZN),DYGF(NS),DYGF(N),
     &      VLSX,INDX)
          DSPTS = DIFMN(DISPT(IZNS),DISPT(IZN),DYGF(NS),DYGF(N),
     &      VLSX,INDX)
          CALL SHDP( VLSX,WLSX,ULSX,DSPLS,DSPTS,DPLS )
        ELSE
          DPLS = 0.D+0
        ENDIF
        DCL = C(NS,1)*FCLS - C(N,1)*FCL
        INDX = 16
        DLY = (DIFMN( DLS,DLP,DYGF(NS),DYGF(N),DCL,INDX)+DPLS)
     &    /DYGP(NPY)/RP(I)
        AL = MAX( VL(1,NPY),ZERO ) +
     &    DLY*MAX( (ONE-(TENTH*ABS(VL(1,NPY))/(DLY+SMALL)))**5,ZERO )
        ALP = MAX( -VL(1,NPY),ZERO ) +
     &    DLY*MAX( (ONE-(TENTH*ABS(VL(1,NPY))/(DLY+SMALL)))**5,ZERO )
        VC(NPY,1) = AFY(NPY)*
     &    ( C(NS,1)*(AL*FCLS + AG*FCGS) - C(N,1)*(ALP*FCL + AGP*FCG) )
  200 CONTINUE
      ENDIF
!
!---  Z-direction VOC transport flux, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
      DO 300 J = 1,JFLD
        DO 300 I = 1,IFLD
          DO 300 K = 2,KFLD
        N = ND(I,J,K)
        NB = N-IJFLD
        IZN = IZ(N)
        IZNB = IZ(NB)
        IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &    INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 300
        NPZ = NSZ(N)
        FCL = 0.D+0
        FCG = 0.D+0
        IF( SL(2,N).GT.SMALL ) FCL = YL(N,1)/(SL(2,N)*PORD(2,N))
        IF( SG(2,N).GT.SMALL ) FCG = YG(N,1)/(SG(2,N)*PORD(2,N))
        FCLB = 0.D+0
        FCGB = 0.D+0
        IF( SL(2,NB).GT.SMALL ) FCLB = YL(NB,1)/(SL(2,NB)*PORD(2,NB))
        IF( SG(2,NB).GT.SMALL ) FCGB = YG(NB,1)/(SG(2,NB)*PORD(2,NB))
        DGP = TORG(2,N)*PORD(2,N)*SG(2,N)*DFGO(2,N)
        DGB = TORG(2,NB)*PORD(2,NB)*SG(2,NB)*DFGO(2,NB)
        IF( IDISP.EQ.1 ) THEN
          CALL ADVB( PORD,SG,UG,VG,WG,UGBX,VGBX,WGBX,N,M )
          INDX = 17
          DSPLB = DIFMN(DISPL(IZNB),DISPL(IZN),DZGF(NB),DZGF(N),
     &      WGBX,INDX)
          DSPTB = DIFMN(DISPT(IZNB),DISPT(IZN),DZGF(NB),DZGF(N),
     &      WGBX,INDX)
          CALL SHDP( WGBX,UGBX,VGBX,DSPLB,DSPTB,DPGB )
        ELSE
          DPGB = 0.D+0
        ENDIF
        DCG = C(NB,1)*FCGB - C(N,1)*FCG
        INDX = 16
        DGZ = (DIFMN( DGB,DGP,DZGF(NB),DZGF(N),DCG,INDX)+DPGB)/DZGP(NPZ)
        AG = MAX( WG(1,NPZ),ZERO ) +
     &    DGZ*MAX( (ONE-(TENTH*ABS(WG(1,NPZ))/(DGZ+SMALL)))**5,ZERO )
        AGP = MAX( -WG(1,NPZ),ZERO ) +
     &    DGZ*MAX( (ONE-(TENTH*ABS(WG(1,NPZ))/(DGZ+SMALL)))**5,ZERO )
        DLP = TORL(2,N)*PORD(2,N)*SL(2,N)*DFLO(2,N)
        DLB = TORL(2,NB)*PORD(2,NB)*SL(2,NB)*DFLO(2,NB)
        IF( IDISP.EQ.1 ) THEN
          CALL ADVB( PORD,SL,UL,VL,WL,ULBX,VLBX,WLBX,N,M )
          INDX = 17
          DSPLB = DIFMN(DISPL(IZNB),DISPL(IZN),DZGF(NB),DZGF(N),
     &      WLBX,INDX)
          DSPTB = DIFMN(DISPT(IZNB),DISPT(IZN),DZGF(NB),DZGF(N),
     &      WLBX,INDX)
          CALL SHDP( WLBX,ULBX,VLBX,DSPLB,DSPTB,DPLB )
        ELSE
          DPLB = 0.D+0
        ENDIF
        DCL = C(NB,1)*FCLB - C(N,1)*FCL
        INDX = 16
        DLZ = (DIFMN( DLB,DLP,DZGF(NB),DZGF(N),DCL,INDX)+DPLB)/DZGP(NPZ)
        AL = MAX( WL(1,NPZ),ZERO ) +
     &    DLZ*MAX( (ONE-(TENTH*ABS(WL(1,NPZ))/(DLZ+SMALL)))**5,ZERO )
        ALP = MAX( -WL(1,NPZ),ZERO ) +
     &    DLZ*MAX( (ONE-(TENTH*ABS(WL(1,NPZ))/(DLZ+SMALL)))**5,ZERO )
        WC(NPZ,1) = AFZ(NPZ)*
     &    ( C(NB,1)*(AL*FCLB + AG*FCGB)- C(N,1)*(ALP*FCL + AGP*FCG) )
  300 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VFLG group
!
      RETURN
      END

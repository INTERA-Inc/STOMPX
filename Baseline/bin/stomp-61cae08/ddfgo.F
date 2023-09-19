!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DDFGO
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
!     Compute advective, diffusive-dispersive oil gas fluxes on
!     interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 January 2004.
!     Last Modified by MD White, PNNL, 25 January 2004.
!     ddfgo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
      USE FLUXD
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
      SUB_LOG(ISUB_LOG) = '/DDFGO'
!
!---  X-direction  ---
!
      IF( IFLD.GT.1 ) THEN
      DO 200 K = 1,KFLD
      DO 200 J = 1,JFLD
      DO 200 I = 2,IFLD
        N = ND(I,J,K)
        NW = N-1
        IZN = IZ(N)
        IZNW = IZ(NW)
        IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &    INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
        NPX = NSX(N)
        DXMGO = XMGO(2,NW)-XMGO(2,N)
        DO 100 M = 1,ISVF
          MN = MNEG(M)
          MP = MPOS(M)
!
!---      Diffusion coefficients  ---
!
          DFSP = TORG(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &      PORD(MP,N)*DFGO(MP,N)*RHOMG(MP,N)
          DFSW = TORG(MN,NW)*(SG(MN,NW)-SGT(MN,NW))*
     &      PORD(MN,NW)*DFGO(MN,NW)*RHOMG(MN,NW)
          INDX = 13
          DFSW = DIFMN(DFSW,DFSP,DXGF(NW),DXGF(N),DXMGO,INDX)
!
!---      Hydraulic dispersion  ---
!
          IF( IDSPD.EQ.1 ) THEN
            CALL ADVW( PORD,SG,UG,VG,WG,UGX,VGX,WGX,N,M )
            UGX = UGX*UGX
            VGX = VGX*VGX
            WGX = WGX*WGX
            ZGW = SQRT(UGX+VGX+WGX)
            INDX = 17
            DPLW = DIFMN(DPLD(IZNW),DPLD(IZN),DXGF(NW),DXGF(N),
     &        UG(1,NPX),INDX)
            DPTW = DIFMN(DPTD(IZNW),DPTD(IZN),DXGF(NW),DXGF(N),
     &        UG(1,NPX),INDX)
            DPSW = (DPLW*UGX + DPTW*(VGX+WGX))/(ZGW+SMALL)
          ELSE
            DPSW = 0.D+0
          ENDIF
!
!---      Gas oil flux by advection, diffusion, and dispersion  ---
!
          DDGW = (DFSW+DPSW)/DXGP(NPX)
          VMCP = (SG(MP,N)-SGT(MP,N))*PORD(MP,N)
          VMCW = (SG(MN,NW)-SGT(MN,NW))*PORD(MN,NW)
          INDX = 13
          VMCW = DIFMN(VMCW,VMCP,DXGF(NW),DXGF(N),UG(1,NPX),INDX)
          CRGW = ABS(UG(M,NPX))*DT/(DXGP(NPX)*VMCW+SMALL)
!
!---      Leonard-TVD and Roe-Superbee schemes  ---
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN
            UDGO(M,NPX) = DDGW*(XMGO(MN,NW)-XMGO(MP,N))
            IF( UG(1,NPX).GE.0.D+0 .AND. I.GT.2 ) THEN
              NWW = NW-1
              IF( IXP(NWW).GT.0 ) THEN
                R = ((XMGO(1,NW)-XMGO(1,NWW))
     &            /(XMGO(1,N)-XMGO(1,NW)+SMALL))
     &            *((DXGF(N)+DXGF(NW))/(DXGF(NW)+DXGF(NWW)))
              ELSE
                R = 0.D+0
              ENDIF
              THETA = FLIMIT( R,CRGW,ISLC(8) )
!              DXF = DXGF(NW)/(DXGF(N)+DXGF(NW))
              DXF = 5.D-1
              UGO(M,NPX) = XGO(MP,N)*RHOG(MP,N)*UG(M,NPX)*THETA*DXF
     &          + XGO(MN,NW)*RHOG(MN,NW)*UG(M,NPX)*(1.D+0-THETA*DXF)
!
!---        Assumes a zero flux west boundary condition, see boundary
!           condition routines for other boundary conditions  ---
!
            ELSEIF( UG(1,NPX).GE.0.D+0 .AND. I.EQ.2 ) THEN
              UGO(M,NPX) = XGO(MN,NW)*RHOG(MN,NW)*UG(M,NPX)
            ELSEIF( UG(1,NPX).LT.0.D+0 .AND. I.LT.IFLD ) THEN
              NE = N+1
              IF( IXP(NE).GT.0 ) THEN
                R = ((XMGO(1,N)-XMGO(1,NE))
     &          /(XMGO(1,NW)-XMGO(1,N)+SMALL))
     &          *((DXGF(NW)+DXGF(N))/(DXGF(N)+DXGF(NE)))
              ELSE
                R = 0.D+0
              ENDIF
              THETA = FLIMIT( R,CRGW,ISLC(8) )
!              DXF = DXGF(N)/(DXGF(N)+DXGF(NW))
              DXF = 5.D-1
              UGO(M,NPX) = XGO(MN,NW)*RHOG(MN,NW)*UG(M,NPX)*THETA*DXF
     &          + XGO(MP,N)*RHOG(MP,N)*UG(M,NPX)*(1.D+0-THETA*DXF)
!
!---        Assumes a zero flux east boundary condition, see boundary
!           condition routines for other boundary conditions  ---
!
            ELSEIF( UG(1,NPX).LT.0.D+0 .AND. I.EQ.IFLD ) THEN
              UGO(M,NPX) = XGO(MP,N)*RHOG(MP,N)*UG(M,NPX)
            ENDIF
            UGO(M,NPX) = UGO(M,NPX) + UDGO(M,NPX)*WTMO
!
!---      Patankar scheme  ---
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            AG = MAX( UG(M,NPX),ZERO ) +
     &        DDGW*MAX((ONE-(TENTH*ABS(UG(M,NPX))/
     &        (DDGW+SMALL)))**5,ZERO)
            AGP = MAX( -UG(M,NPX),ZERO ) +
     &        DDGW*MAX((ONE-(TENTH*ABS(UG(M,NPX))/
     &        (DDGW+SMALL)))**5,ZERO)
            UGO(M,NPX) = WTMO*(XMGO(MN,NW)*AG -
     &        XMGO(MP,N)*AGP)
            UDGO(M,NPX) = DDGW*(XMGO(MN,NW)-XMGO(MP,N))
!
!---      Upwind scheme  ---
!
          ELSE
            UDGO(M,NPX) = DDGW*(XMGO(MN,NW)-XMGO(MP,N))
            FGOP = XGO(MP,N)*RHOG(MP,N)
            FGOW = XGO(MN,NW)*RHOG(MN,NW)
            INDX = -4
            FGO = DIFMN( FGOW,FGOP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            UGO(M,NPX) = UG(M,NPX)*FGO + WTMO*UDGO(M,NPX)
          ENDIF
  100   CONTINUE
  200 CONTINUE
      ENDIF
!
!---  Y-direction  ---
!
      IF( JFLD.GT.1 ) THEN
      DO 400 I = 1,IFLD
      DO 400 K = 1,KFLD
      DO 400 J = 2,JFLD
        N = ND(I,J,K)
        NS = N-IFLD
        IZN = IZ(N)
        IZNS = IZ(NS)
        IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &    INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
        NPY = NSY(N)
        DXMGO = XMGO(2,NS)-XMGO(2,N)
        DO 300 M = 1,ISVF
          MN = MNEG(M)
          MP = MPOS(M)
!
!---      Diffusion coefficients  ---
!
          DFSP = TORG(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &      PORD(MP,N)*DFGO(MP,N)*RHOMG(MP,N)
          DFSS = TORG(MN,NS)*(SG(MN,NS)-SGT(MN,NS))*
     &      PORD(MN,NS)*DFGO(MN,NS)*RHOMG(MN,NS)
          INDX = 13
          DFSS = DIFMN(DFSS,DFSP,DYGF(NS),DYGF(N),DXMGO,INDX)
!
!---      Hydraulic dispersion  ---
!
          IF( IDSPD.EQ.1 ) THEN
            CALL ADVS( PORD,SG,UG,VG,WG,UGX,VGX,WGX,N,M )
            UGX = UGX*UGX
            VGX = VGX*VGX
            WGX = WGX*WGX
            ZGS = SQRT(UGX+VGX+WGX)
            INDX = 17
            DPLS = DIFMN(DPLD(IZNS),DPLD(IZN),DYGF(NS),DYGF(N),
     &        VG(1,NPY),INDX)
            DPTS = DIFMN(DPTD(IZNS),DPTD(IZN),DYGF(NS),DYGF(N),
     &        VG(1,NPY),INDX)
            DPSS = (DPLS*VGX + DPTS*(WGX+UGX))/(ZGS+SMALL)
          ELSE
            DPSS = 0.D+0
          ENDIF
!
!---      Gas oil flux by advection, diffusion, and dispersion  ---
!
          DDGS = (DFSS+DPSS)/DYGP(NPY)/RP(I)
          VMCP = (SG(MP,N)-SGT(MP,N))*PORD(MP,N)
          VMCS = (SG(MN,NS)-SGT(MN,NS))*PORD(MN,NS)
          INDX = 13
          VMCS = DIFMN(VMCS,VMCP,DYGF(NS),DYGF(N),VG(1,NPY),INDX)
          CRGS = ABS(VG(M,NPY))*DT/(DYGP(NPY)*VMCS+SMALL)/RP(I)
!
!---      Leonard-TVD and Roe-Superbee schemes  ---
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN
            VDGO(M,NPY) = DDGS*(XMGO(MN,NS)-XMGO(MP,N))
            IF( VG(1,NPY).GE.0.D+0 .AND. J.GT.2 ) THEN
              NSS = NS-IFLD
              IF( IXP(NSS).GT.0 ) THEN
                R = ((XMGO(1,NS)-XMGO(1,NSS))
     &            /(XMGO(1,N)-XMGO(1,NS)+SMALL))
     &            *((DYGF(N)+DYGF(NS))/(DYGF(NS)+DYGF(NSS)))
              ELSE
                R = 0.D+0
              ENDIF
              THETA = FLIMIT( R,CRGS,ISLC(8) )
!              DYF = DYGF(NS)/(DYGF(N)+DYGF(NS))
              DYF = 5.D-1
              VGO(M,NPY) = XGO(MP,N)*RHOG(MP,N)*VG(M,NPY)*THETA*DYF
     &          + XGO(MN,NS)*RHOG(MN,NS)*VG(M,NPY)*(1.D+0-THETA*DYF)
!
!---        Assumes a zero flux south boundary condition, see boundary
!           condition routines for other boundary conditions  ---
!
            ELSEIF( VG(1,NPY).GE.0.D+0 .AND. J.EQ.2 ) THEN
              VGO(M,NPY) = XGO(MN,NS)*RHOG(MN,NS)*VG(M,NPY)
            ELSEIF( VG(1,NPY).LT.0.D+0 .AND. J.LT.JFLD ) THEN
              NN = N+IFLD
              IF( IXP(NN).GT.0 ) THEN
                R = ((XMGO(1,N)-XMGO(1,NN))
     &            /(XMGO(1,NS)-XMGO(1,N)+SMALL))
     &            *((DYGF(NS)+DYGF(N))/(DYGF(N)+DYGF(NN)))
              ELSE
                R = 0.D+0
              ENDIF
              THETA = FLIMIT( R,CRGS,ISLC(8) )
!              DYF = DYGF(N)/(DYGF(N)+DYGF(NS))
              DYF = 5.D-1
              VGO(M,NPY) = XGO(MN,NS)*RHOG(MN,NS)*VG(M,NPY)*THETA*DYF
     &          + XGO(MP,N)*RHOG(MP,N)*VG(M,NPY)*(1.D+0-THETA*DYF)
!
!---        Assumes a zero flux north boundary condition, see boundary
!           condition routines for other boundary conditions  ---
!
            ELSEIF( VG(1,NPY).LT.0.D+0 .AND. J.EQ.JFLD ) THEN
              VGO(M,NPY) = XGO(MP,N)*RHOG(MP,N)*VG(M,NPY)
            ENDIF
            VGO(M,NPY) = VGO(M,NPY) + VDGO(M,NPY)*WTMO
!
!---      Patankar scheme  ---
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            AG = MAX( VG(M,NPY),ZERO ) +
     &        DDGS*MAX((ONE-(TENTH*ABS(VG(M,NPY))/
     &        (DDGS+SMALL)))**5,ZERO)
            AGP = MAX( -VG(M,NPY),ZERO ) +
     &        DDGS*MAX((ONE-(TENTH*ABS(VG(M,NPY))/
     &        (DDGS+SMALL)))**5,ZERO)
            VGO(M,NPY) = WTMO*(XMGO(MN,NS)*AG -
     &        XMGO(MP,N)*AGP)
            VDGO(M,NPY) = DDGS*(XMGO(MN,NS)-XMGO(MP,N))
!
!---      Upwind scheme  ---
!
          ELSE
            VDGO(M,NPY) = DDGS*(XMGO(MN,NS)-XMGO(MP,N))
            FGOP = XGO(MP,N)*RHOG(MP,N)
            FGOS = XGO(MN,NS)*RHOG(MN,NS)
            INDX = -4
            FGO = DIFMN( FGOS,FGOP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
            VGO(M,NPY) = VG(M,NPY)*FGO + WTMO*VDGO(M,NPY)
          ENDIF
  300   CONTINUE
  400 CONTINUE
      ENDIF
!
!---  Z-direction  ---
!
      IF( KFLD.GT.1 ) THEN
      DO 600 J = 1,JFLD
      DO 600 I = 1,IFLD
      DO 600 K = 2,KFLD
        N = ND(I,J,K)
        NB = N-IJFLD
        IZN = IZ(N)
        IZNB = IZ(NB)
        IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &    INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
        NPZ = NSZ(N)
        DXMGO = XMGO(2,NB)-XMGO(2,N)
        DO 500 M = 1,ISVF
          MN = MNEG(M)
          MP = MPOS(M)
!
!---      Diffusion coefficients  ---
!
          DFSP = TORG(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &      PORD(MP,N)*DFGO(MP,N)*RHOMG(MP,N)
          DFSB = TORG(MN,NB)*(SG(MN,NB)-SGT(MN,NB))*
     &      PORD(MN,NB)*DFGO(MN,NB)*RHOMG(MN,NB)
          INDX = 13
          DFSB = DIFMN(DFSB,DFSP,DZGF(NB),DZGF(N),DXMGO,INDX)
!
!---      Hydraulic dispersion  ---
!
          IF( IDSPD.EQ.1 ) THEN
            CALL ADVB( PORD,SG,UG,VG,WG,UGX,VGX,WGX,N,M )
            UGX = UGX*UGX
            VGX = VGX*VGX
            WGX = WGX*WGX
            ZGBX = SQRT(UGX+VGX+WGX)
            INDX = 17
            DPLB = DIFMN(DPLD(IZNB),DPLD(IZN),DZGF(NB),DZGF(N),
     &        WG(1,NPZ),INDX)
            DPTB = DIFMN(DPTD(IZNB),DPTD(IZN),DZGF(NB),DZGF(N),
     &        WG(1,NPZ),INDX)
            DPSB = (DPLB*WGX + DPTB*(UGX+VGX))/(ZGBX+SMALL)
          ELSE
            DPSB = 0.D+0
          ENDIF
!
!---      Gas oil flux by advection, diffusion, and dispersion  ---
!
          DDGB = (DFSB+DPSB)/DZGP(NPZ)
          VMCP = (SG(MP,N)-SGT(MP,N))*PORD(MP,N)
          VMCB = (SG(MN,NB)-SGT(MN,NB))*PORD(MN,NB)
          INDX = 13
          VMCB = DIFMN(VMCB,VMCP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX)
          CRGB = ABS(WG(M,NPZ))*DT/(DZGP(NPZ)*VMCB+SMALL)
!
!---      Leonard-TVD and Roe-Superbee schemes  ---
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN
            WDGO(M,NPZ) = DDGB*(XMGO(MN,NB)-XMGO(MP,N))
            IF( WG(1,NPZ).GE.0.D+0 .AND. K.GT.2 ) THEN
              NBB = NB-IJFLD
              IF( IXP(NBB).GT.0 ) THEN
                R = ((XMGO(1,NB)-XMGO(1,NBB))
     &            /(XMGO(1,N)-XMGO(1,NB)+SMALL))
     &            *((DZGF(N)+DZGF(NB))/(DZGF(NB)+DZGF(NBB)))
              ELSE
                R = 0.D+0
              ENDIF
              THETA = FLIMIT( R,CRGB,ISLC(8) )
!              DZF = DZGF(NB)/(DZGF(N)+DZGF(NB))
              DZF = 5.D-1
              WGO(M,NPZ) = XGO(MP,N)*RHOG(MP,N)*WG(M,NPZ)*THETA*DZF
     &          + XGO(MN,NB)*RHOG(MN,NB)*WG(M,NPZ)*(1.D+0-THETA*DZF)
!
!---        Assumes a zero flux bottom boundary condition, see boundary
!           condition routines for other boundary conditions  ---
!
            ELSEIF( WG(1,NPZ).GE.0.D+0 .AND. K.EQ.2 ) THEN
              WGO(M,NPZ) = XGO(MN,NB)*RHOG(MN,NB)*WG(M,NPZ)
            ELSEIF( WG(1,NPZ).LT.0.D+0 .AND. K.LT.KFLD ) THEN
              NT = N+IJFLD
              IF( IXP(NT).GT.0 ) THEN
               R = ((XMGO(1,N)-XMGO(1,NT))
     &           /(XMGO(1,NB)-XMGO(1,N)+SMALL))
     &           *((DZGF(NB)+DZGF(N))/(DZGF(N)+DZGF(NT)))
              ELSE
                R = 0.D+0
              ENDIF
              THETA = FLIMIT( R,CRGB,ISLC(8) )
!              DZF = DZGF(N)/(DZGF(N)+DZGF(NB))
              DZF = 5.D-1
              WGO(M,NPZ) = XGO(MN,NB)*RHOG(MN,NB)*WG(M,NPZ)*THETA*DZF
     &          + XGO(MP,N)*RHOG(MP,N)*WG(M,NPZ)*(1.D+0-THETA*DZF)
!
!---        Assumes a zero flux top boundary condition, see boundary
!           condition routines for other boundary conditions  ---
!
            ELSEIF( WG(1,NPZ).LT.0.D+0 .AND. K.EQ.KFLD ) THEN
              WGO(M,NPZ) = XGO(MP,N)*RHOG(MP,N)*WG(M,NPZ)
            ENDIF
             WGO(M,NPZ) = WGO(M,NPZ) + WDGO(M,NPZ)*WTMO
!
!---      Patankar scheme  ---
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            AG = MAX( WG(M,NPZ),ZERO ) +
     &        DDGB*MAX((ONE-(TENTH*ABS(WG(M,NPZ))/
     &        (DDGB+SMALL)))**5,ZERO)
            AGP = MAX( -WG(M,NPZ),ZERO ) +
     &        DDGB*MAX((ONE-(TENTH*ABS(WG(M,NPZ))/
     &        (DDGB+SMALL)))**5,ZERO)
            WGO(M,NPZ) = WTMO*(XMGO(MN,NB)*AG -
     &        XMGO(MP,N)*AGP)
            WDGO(M,NPZ) = DDGB*(XMGO(MN,NB)-XMGO(MP,N))
!
!---      Upwind scheme  ---
!
          ELSE
            WDGO(M,NPZ) = DDGB*(XMGO(MN,NB)-XMGO(MP,N))
            FGOP = XGO(MP,N)*RHOG(MP,N)
            FGOB = XGO(MN,NB)*RHOG(MN,NB)
            INDX = -4
            FGO = DIFMN( FGOB,FGOP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
            WGO(M,NPZ) = WG(M,NPZ)*FGO + WTMO*WDGO(M,NPZ)
          ENDIF
  500   CONTINUE
  600 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DDFGO group  ---
!
      RETURN
      END


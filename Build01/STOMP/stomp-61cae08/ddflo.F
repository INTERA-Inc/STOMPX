!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DDFLO
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
!     Compute advective, diffusive, dispersive aqueous oil fluxes 
!     on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1996.
!     Last Modified by MD White on August 1, 1996.
!     ddflo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/DDFLO'
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
        DXMLO = XMLO(2,NW)-XMLO(2,N)
        DO 100 M = 1,ISVF
          MN = MNEG(M)
          MP = MPOS(M)
!
!---      Diffusion coefficients  ---
!
          DFSP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLO(MP,N)
     &      *RHOML(MP,N)
          DFSW = TORL(MN,NW)*SL(MN,NW)*PORD(MN,NW)*DFLO(MN,NW)
     &      *RHOML(MN,NW)
          INDX = 15
          DFSW = DIFMN(DFSW,DFSP,DXGF(NW),DXGF(N),DXMLO,INDX)
!
!---      Hydraulic dispersion  ---
!
          IF( IDSPD.EQ.1 ) THEN
            CALL ADVW( PORD,SL,UL,VL,WL,ULX,VLX,WLX,N,M )
            ULX = ULX*ULX
            VLX = VLX*VLX
            WLX = WLX*WLX
            ZVW = SQRT(ULX+VLX+WLX)
            INDX = 17
            DPLW = DIFMN(DPLD(IZNW),DPLD(IZN),DXGF(NW),DXGF(N),
     &        UL(1,NPX),INDX)
            DPTW = DIFMN(DPTD(IZNW),DPTD(IZN),DXGF(NW),DXGF(N),
     &        UL(1,NPX),INDX)
            DPSW = (DPLW*ULX + DPTW*(VLX+WLX))/(ZVW+SMALL)
          ELSE
            DPSW = 0.D+0
          ENDIF
!
!---      Aqueous oil flux by advection, diffusion, 
!         and dispersion  ---
!
          DDLW = (DFSW+DPSW)/DXGP(NPX)
          VMCP = SL(MP,N)*PORD(MP,N)
          VMCW = SL(MN,NW)*PORD(MN,NW)
          INDX = 15
          VMCW = DIFMN(VMCW,VMCP,DXGF(NW),DXGF(N),UL(1,NPX),INDX)
          CRLW = ABS(UL(M,NPX))*DT/(DXGP(NPX)*VMCW+SMALL)
!
!---      Leonard-TVD and Roe-Superbee schemes  ---
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN
            UDLO(M,NPX) = DDLW*(XMLO(MN,NW)-XMLO(MP,N))
            IF( UL(1,NPX).GE.0.D+0 .AND. I.GT.2 ) THEN
              NWW = NW-1
              IF( IXP(NWW).GT.0 ) THEN
                R = ((XMLO(1,NW)-XMLO(1,NWW))
     &            /(XMLO(1,N)-XMLO(1,NW)+SMALL))
     &            *((DXGF(N)+DXGF(NW))/(DXGF(NW)+DXGF(NWW)))
              ELSE
                R = 0.D+0
              ENDIF
              THETA = FLIMIT( R,CRLW,ISLC(8) )
!              DXF = DXGF(NW)/(DXGF(N)+DXGF(NW))
              DXF = 5.D-1
              ULO(M,NPX) = XLO(MP,N)*RHOL(MP,N)*UL(M,NPX)*THETA*DXF
     &          + XLO(MN,NW)*RHOL(MN,NW)*UL(M,NPX)*(1.D+0-THETA*DXF)
!
!---        Assumes a zero flux west boundary condition, see boundary
!           condition routines for other boundary conditions  ---
!
            ELSEIF( UL(1,NPX).GE.0.D+0 .AND. I.EQ.2 ) THEN
              ULO(M,NPX) = XLO(MN,NW)*RHOL(MN,NW)*UL(M,NPX)
            ELSEIF( UL(1,NPX).LT.0.D+0 .AND. I.LT.IFLD ) THEN
              NE = N+1
              IF( IXP(NE).GT.0 ) THEN
                R = ((XMLO(1,N)-XMLO(1,NE))
     &          /(XMLO(1,NW)-XMLO(1,N)+SMALL))
     &          *((DXGF(NW)+DXGF(N))/(DXGF(N)+DXGF(NE)))
              ELSE
                R = 0.D+0
              ENDIF
              THETA = FLIMIT( R,CRLW,ISLC(8) )
!              DXF = DXGF(N)/(DXGF(N)+DXGF(NW))
              DXF = 5.D-1
              ULO(M,NPX) = XLO(MN,NW)*RHOL(MN,NW)*UL(M,NPX)*THETA*DXF
     &          + XLO(MP,N)*RHOL(MP,N)*UL(M,NPX)*(1.D+0-THETA*DXF)
!
!---        Assumes a zero flux east boundary condition, see boundary
!           condition routines for other boundary conditions  ---
!
            ELSEIF( UL(1,NPX).LT.0.D+0 .AND. I.EQ.IFLD ) THEN
              ULO(M,NPX) = XLO(MP,N)*RHOL(MP,N)*UL(M,NPX)
            ENDIF
            ULO(M,NPX) = ULO(M,NPX) + UDLO(M,NPX)*WTMO
!
!---      Patankar scheme  ---
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            AL = MAX( UL(M,NPX),0.D+0 ) +
     &        DDLW*MAX((ONE-(TENTH*ABS(UL(M,NPX))/
     &        (DDLW+SMALL)))**5,0.D+0)
            ALP = MAX( -UL(M,NPX),0.D+0 ) +
     &        DDLW*MAX((ONE-(TENTH*ABS(UL(M,NPX))/
     &        (DDLW+SMALL)))**5,0.D+0)
            ULO(M,NPX) = WTMO*(XMLO(MN,NW)*AL-XMLO(MP,N)*ALP)
            UDLO(M,NPX) = DDLW*(XMLO(MN,NW)-XMLO(MP,N))
!
!---      Upwind scheme  ---
!
          ELSE
            UDLO(M,NPX) = DDLW*(XMLO(MN,NW) -
     &        XMLO(MP,N))
            FLOP = XLO(MP,N)*RHOL(MP,N)
            FLOW = XLO(MN,NW)*RHOL(MN,NW)
            INDX = -4
            FLO = DIFMN( FLOW,FLOP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            ULO(M,NPX) = UL(M,NPX)*FLO + WTMO*UDLO(M,NPX)
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
        DXMLO = (XMLO(2,NS)-XMLO(2,N))
        DO 300 M = 1,ISVF
          MN = MNEG(M)
          MP = MPOS(M)
!
!---      Diffusion coefficients  ---
!
          DFSP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLO(MP,N)
     &      *RHOML(MP,N)
          DFSS = TORL(MN,NS)*SL(MN,NS)*PORD(MN,NS)*DFLO(MN,NS)
     &      *RHOML(MN,NS)
          INDX = 15
          DFSS = DIFMN(DFSS,DFSP,DYGF(NS),DYGF(N),DXMLO,INDX)
!
!---      Hydraulic dispersion  ---
!
          IF( IDSPD.EQ.1 ) THEN
            CALL ADVS( PORD,SL,UL,VL,WL,ULX,VLX,WLX,N,M )
            ULX = ULX*ULX
            VLX = VLX*VLX
            WLX = WLX*WLX
            ZVS = SQRT(ULX+VLX+WLX)
            INDX = 17
            DPLS = DIFMN(DPLD(IZNS),DPLD(IZN),DYGF(NS),DYGF(N),
     &        VL(1,NPY),INDX)
            DPTS = DIFMN(DPTD(IZNS),DPTD(IZN),DYGF(NS),DYGF(N),
     &        VL(1,NPY),INDX)
            DPSS = (DPLS*VLX + DPTS*(WLX+ULX))/(ZVS+SMALL)
          ELSE
            DPSS = 0.D+0
          ENDIF
!
!---      Aqueous oil flux by advection, diffusion,
!         and dispersion  ---
!
          DDLS = (DFSS+DPSS)/DYGP(NPY)/RP(I)
          VMCP = SL(MP,N)*PORD(MP,N)
          VMCS = SL(MN,NS)*PORD(MN,NS)
          INDX = 15
          VMCS = DIFMN(VMCS,VMCP,DYGF(NS),DYGF(N),VL(1,NPY),INDX)
          CRLS = ABS(VL(M,NPY))*DT/(DYGP(NPY)*VMCS+SMALL)/RP(I)
!
!---      Leonard-TVD and Roe-Superbee schemes  ---
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN
            VDLO(M,NPY) = DDLS*(XMLO(MN,NS)-XMLO(MP,N))
            IF( VL(1,NPY).GE.0.D+0 .AND. J.GT.2 ) THEN
              NSS = NS-IFLD
              IF( IXP(NSS).GT.0 ) THEN
                R = ((XMLO(1,NS)-XMLO(1,NSS))
     &            /(XMLO(1,N)-XMLO(1,NS)+SMALL))
     &            *((DYGF(N)+DYGF(NS))/(DYGF(NS)+DYGF(NSS)))
              ELSE
                R = 0.D+0
              ENDIF
              THETA = FLIMIT( R,CRLS,ISLC(8) )
!              DYF = DYGF(NS)/(DYGF(N)+DYGF(NS))
              DYF = 5.D-1
              VLO(M,NPY) = XLO(MP,N)*RHOL(MP,N)*VL(M,NPY)*THETA*DYF
     &          + XLO(MN,NS)*RHOL(MN,NS)*VL(M,NPY)*(1.D+0-THETA*DYF)
!
!---        Assumes a zero flux south boundary condition, see boundary
!           condition routines for other boundary conditions  ---
!
            ELSEIF( VL(1,NPY).GE.0.D+0 .AND. J.EQ.2 ) THEN
              VLO(M,NPY) = XLO(MN,NS)*RHOL(MN,NS)*VL(M,NPY)
            ELSEIF( VL(1,NPY).LT.0.D+0 .AND. J.LT.JFLD ) THEN
              NN = N+IFLD
              IF( IXP(NN).GT.0 ) THEN
                R = ((XMLO(1,N)-XMLO(1,NN))
     &            /(XMLO(1,NS)-XMLO(1,N)+SMALL))
     &            *((DYGF(NS)+DYGF(N))/(DYGF(N)+DYGF(NN)))
              ELSE
                R = 0.D+0
              ENDIF
              THETA = FLIMIT( R,CRLS,ISLC(8) )
!              DYF = DYGF(N)/(DYGF(N)+DYGF(NS))
              DYF = 5.D-1
              VLO(M,NPY) = XLO(MN,NS)*RHOL(MN,NS)*VL(M,NPY)*THETA*DYF
     &          + XLO(MP,N)*RHOL(MP,N)*VL(M,NPY)*(1.D+0-THETA*DYF)
!
!---        Assumes a zero flux north boundary condition, see boundary
!           condition routines for other boundary conditions  ---
!
            ELSEIF( VL(1,NPY).LT.0.D+0 .AND. J.EQ.JFLD ) THEN
              VLO(M,NPY) = XLO(MP,N)*RHOL(MP,N)*VL(M,NPY)
            ENDIF
            VLO(M,NPY) = VLO(M,NPY) + VDLO(M,NPY)*WTMO
!
!---      Patankar scheme  ---
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            AL = MAX( VL(M,NPY),0.D+0 ) +
     &        DDLS*MAX((ONE-(TENTH*ABS(VL(M,NPY))/
     &        (DDLS+SMALL)))**5,0.D+0)
            ALP = MAX( -VL(M,NPY),0.D+0 ) +
     &        DDLS*MAX((ONE-(TENTH*ABS(VL(M,NPY))/
     &        (DDLS+SMALL)))**5,0.D+0)
            VLO(M,NPY) = WTMO*(XMLO(MN,NS)*AL-XMLO(MP,N)*ALP)
            VDLO(M,NPY) = DDLS*(XMLO(MN,NS)-XMLO(MP,N))
!
!---      Upwind scheme  ---
!
          ELSE
            VDLO(M,NPY) = DDLS*(XMLO(MN,NS)-XMLO(MP,N))
            FLOP = XLO(MP,N)*RHOL(MP,N)
            FLOS = XLO(MN,NS)*RHOL(MN,NS)
            INDX = -4
            FLO = DIFMN( FLOS,FLOP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
            VLO(M,NPY) = VL(M,NPY)*FLO + WTMO*VDLO(M,NPY)
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
        DXMLO = (XMLO(2,NB)-XMLO(2,N))
        DO 500 M = 1,ISVF
          MN = MNEG(M)
          MP = MPOS(M)
!
!---      Diffusion coefficients  ---
!
          DFSP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLO(MP,N)
     &      *RHOML(MP,N)
          DFSB = TORL(MN,NB)*SL(MN,NB)*PORD(MN,NB)*DFLO(MN,NB)
     &      *RHOML(MN,NB)
          INDX = 15
          DFSB = DIFMN(DFSB,DFSP,DZGF(NB),DZGF(N),DXMLO,INDX)
!
!---      Hydraulic dispersion  ---
!
          IF( IDSPD.EQ.1 ) THEN
            CALL ADVB( PORD,SL,UL,VL,WL,ULX,VLX,WLX,N,M )
            ULX = ULX*ULX
            VLX = VLX*VLX
            WLX = WLX*WLX
            ZVB = SQRT(ULX+VLX+WLX)
            INDX = 17
            DPLB = DIFMN(DPLD(IZNB),DPLD(IZN),DZGF(NB),DZGF(N),
     &        WL(1,NPZ),INDX)
            DPTB = DIFMN(DPTD(IZNB),DPTD(IZN),DZGF(NB),DZGF(N),
     &        WL(1,NPZ),INDX)
            DPSB = (DPLB*WLX + DPTB*(ULX+VLX))/(ZVB+SMALL)
          ELSE
            DPSB = 0.D+0
          ENDIF
!
!---      Aqueous oil flux by advection, diffusion,
!         and dispersion  ---
!
          DDLB = (DFSB+DPSB)/DZGP(NPZ)
          VMCP = SL(MP,N)*PORD(MP,N)
          VMCB = SL(MN,NB)*PORD(MN,NB)
          INDX = 15
          VMCB = DIFMN(VMCB,VMCP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX)
          CRLB = ABS(WL(M,NPZ))*DT/(DZGP(NPZ)*VMCB+SMALL)
!
!---      Leonard-TVD and Roe-Superbee schemes  ---
!
          IF( ISLC(8).EQ.1 .OR. ISLC(8).EQ.2 ) THEN
            WDLO(M,NPZ) = DDLB*(XMLO(MN,NB) -
     &        XMLO(MP,N))
            IF( WL(1,NPZ).GE.0.D+0 .AND. K.GT.2 ) THEN
              NBB = NB-IJFLD
              IF( IXP(NBB).GT.0 ) THEN
                R = ((XMLO(1,NB)-XMLO(1,NBB))
     &            /(XMLO(1,N)-XMLO(1,NB)+SMALL))
     &            *((DZGF(N)+DZGF(NB))/(DZGF(NB)+DZGF(NBB)))
              ELSE
                R = 0.D+0
              ENDIF
              THETA = FLIMIT( R,CRLB,ISLC(8) )
!              DZF = DZGF(NB)/(DZGF(N)+DZGF(NB))
              DZF = 5.D-1
              WLO(M,NPZ) = XLO(MP,N)*RHOL(MP,N)*WL(M,NPZ)*THETA*DZF
     &          + XLO(MN,NB)*RHOL(MN,NB)*WL(M,NPZ)*(1.D+0-THETA*DZF)
!
!---        Assumes a zero flux bottom boundary condition, see boundary
!           condition routines for other boundary conditions  ---
!
            ELSEIF( WL(1,NPZ).GE.0.D+0 .AND. K.EQ.2 ) THEN
              WLO(M,NPZ) = XLO(MN,NB)*RHOL(MN,NB)*WL(M,NPZ)
            ELSEIF( WL(1,NPZ).LT.0.D+0 .AND. K.LT.KFLD ) THEN
              NT = N+IJFLD
              IF( IXP(NT).GT.0 ) THEN
               R = ((XMLO(1,N)-XMLO(1,NT))
     &           /(XMLO(1,NB)-XMLO(1,N)+SMALL))
     &           *((DZGF(NB)+DZGF(N))/(DZGF(N)+DZGF(NT)))
              ELSE
                R = 0.D+0
              ENDIF
              THETA = FLIMIT( R,CRLB,ISLC(8) )
!              DZF = DZGF(N)/(DZGF(N)+DZGF(NB))
              DZF = 5.D-1
              WLO(M,NPZ) = XLO(MN,NB)*RHOL(MN,NB)*WL(M,NPZ)*THETA*DZF
     &          + XLO(MP,N)*RHOL(MP,N)*WL(M,NPZ)*(1.D+0-THETA*DZF)
!
!---        Assumes a zero flux top boundary condition, see boundary
!           condition routines for other boundary conditions  ---
!
            ELSEIF( WL(1,NPZ).LT.0.D+0 .AND. K.EQ.KFLD ) THEN
              WLO(M,NPZ) = XLO(MP,N)*RHOL(MP,N)*WL(M,NPZ)
            ENDIF
            WLO(M,NPZ) = WLO(M,NPZ) + WDLO(M,NPZ)*WTMO
!
!---      Patankar scheme  ---
!
          ELSEIF( ISLC(8).EQ.0 ) THEN
            AL = MAX( WL(M,NPZ),0.D+0 ) +
     &        DDLB*MAX((ONE-(TENTH*ABS(WL(M,NPZ))/
     &        (DDLB+SMALL)))**5,0.D+0)
            ALP = MAX( -WL(M,NPZ),0.D+0 ) +
     &        DDLB*MAX((ONE-(TENTH*ABS(WL(M,NPZ))/
     &        (DDLB+SMALL)))**5,0.D+0)
            WLO(M,NPZ) = WTMO*(XMLO(MN,NB)*AL-XMLO(MP,N)*ALP)
            WDLO(M,NPZ) = DDLB*(XMLO(MN,NB)-XMLO(MP,N))
!
!---      Upwind scheme  ---
!
          ELSE
            WDLO(M,NPZ) = DDLB*(XMLO(MN,NB)-XMLO(MP,N))
            FLOP = XLO(MP,N)*RHOL(MP,N)
            FLOB = XLO(MN,NB)*RHOL(MN,NB)
            INDX = -4
            FLO = DIFMN( FLOB,FLOP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
            WLO(M,NPZ) = WL(M,NPZ)*FLO + WTMO*WDLO(M,NPZ)
          ENDIF
  500   CONTINUE
  600 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DDFLO group  ---
!
      RETURN
      END


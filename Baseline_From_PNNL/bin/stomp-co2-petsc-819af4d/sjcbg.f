!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SJCBG( NSL )
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
!     Loads the matrix elements and solution vector for the
!     gas-phase convective-dispersive mass transport equation.
!
!     The Jacobian matrix is initially configured assuming zero-flux
!     boundary conditions.  The matrix is then updated for other
!     user-specified boundary conditions.
!
!     Matrix elements are stored in the array ALU.
!     Elements for the right-hand-side are stored in the array BLU.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1996.
!     Last Modified by MD White on September 5, 1996.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE JACOB
      USE HYST
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
      SUB_LOG(ISUB_LOG) = '/SJCBG'
!
!---  Fill matrix elements  ---
!
      DO 900 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 900
        I = ID(N)
        J = JD(N)
        K = KD(N)
        NW = N-1
        NE = N+1
        NS = N-IFLD
        NN = N+IFLD
        NB = N-IJFLD
        NT = N+IJFLD
        IZN = IZ(N)
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
        DYR = RP(I)*DYGF(N)
!
!---    Storage terms  ---
!
        SC = VOL(N)*DTI
        CC = 0.D+0
        MCP = IXP(N)
        MA = 1
!
!---    Banded solver  ---
!
        IF( ILES.EQ.1 ) THEN
          MCD = MCP
          MRD = MDT
          IF( IEQW.EQ.0 ) ALU(MRD,MCD) = ALU(MRD,MCD) + SC
!
!---    SPLib solver  ---
!
        ELSEIF( ILES.EQ.3 ) THEN
          MCD = KLUC(MCP,MA)
          MA = MA + 1
          IF( IEQW.EQ.0 ) DLU(MCD) = DLU(MCD) + SC

!
!---    PETSc solver  ---
!
        ELSEIF( ILES.EQ.5 ) THEN
          MCD = KLUC(MCP,MA)
          MA = MA + 1
          IF( IEQW.EQ.0 ) DLU(MCD) = DLU(MCD) + SC

        ENDIF
!
!---    Molecular diffusion coefficients at the nodes  ---
!
        TCOR = (T(2,N)+TABS)/TSPRF
        PCOR = (PG(2,N)+PATM)/PATM
        SDFGP = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
        VMCP = SG(2,N)*PORD(2,N)
        FCGP = YG(N,NSL)/(VMCP+SMALL)
        DGP = TORG(2,N)*(SG(2,N)-SGT(2,N))*PORD(2,N)*SDFGP
!
!---    Hydrodynamic dispersion coefficients at cell faces  ---
!
        CALL SHDPG( N,DPGB,DPGS,DPGW,DPGE,DPGN,DPGT )
!
!---    Bottom face diffusion and advection terms  ---
!
        IF( K.NE.1 ) THEN
          IF( IXP(NB).EQ.0 .OR. INBS(1,N).GT.0 ) GOTO 110
          TCOR = (T(2,NB)+TABS)/TSPRF
          PCOR = (PG(2,NB)+PATM)/PATM
          SDFGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCB = SG(2,NB)*PORD(2,NB)
          FCGB = YG(NB,NSL)/(VMCB+SMALL)
          DGB = TORG(2,NB)*(SG(2,NB)-SGT(2,NB))*PORD(2,NB)*SDFGB
          INDX = 16
          DGZ = DIFMN(DGB,DGP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX)
          DGZ = AFZ(NPZ)*(DGZ+DPGB)/DZGP(NPZ)
          IF( MOD(ISLC(23),100)/10.EQ.1 ) FGB = 0.D+0
          FGB = AFZ(NPZ)*WG(1,NPZ)
          VMCBX = DIFMN(VMCB,VMCP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX)
          CRGB = ABS(WG(1,NPZ))*DT/(DZGP(NPZ)*VMCBX+SMALL)
!
!---      TVD solute transport  ---
!
          IF( ISLC(1).GE.1 ) THEN
            IF( FGB.GE.ZERO .AND. K.GT.2 ) THEN
             IF( IXP(NB-IJFLD).GT.0 ) THEN
              FCGBB = YG(NB-IJFLD,NSL)/(SG(2,NB-IJFLD)
     &          *PORD(2,NB-IJFLD)+SMALL)
              R = ((CO(NB,NSL)*FCGB-CO(NB-IJFLD,NSL)*FCGBB)
     &          /(CO(N,NSL)*FCGP-CO(NB,NSL)*FCGB+EPSL))
     &          *((DZGF(N)+DZGF(NB))/(DZGF(NB)+DZGF(NB-IJFLD)))
              THETA = FLIMIT( R,CRGB,ISLC(1) )
              DZF = DZGF(NB)/(DZGF(N)+DZGF(NB))
              WCZ = CO(NB,NSL)*FGB*(1.D+0-THETA*DZF)*FCGB
     &          + CO(N,NSL)*FGB*THETA*DZF*FCGP
              WC(NPZ,NSL) = WC(NPZ,NSL) + WCZ/AFZ(NPZ)
              CC = CC + WCZ
             ELSE
              WCZ = CO(NB,NSL)*FGB*FCGB
              WC(NPZ,NSL) = WC(NPZ,NSL) + WCZ/AFZ(NPZ)
              CC = CC + WCZ
             ENDIF
            ELSEIF( FGB.GE.ZERO .AND. K.EQ.2 ) THEN
              WCZ = CO(NB,NSL)*FGB*FCGB
              WC(NPZ,NSL) = WC(NPZ,NSL) + WCZ/AFZ(NPZ)
              CC = CC + WCZ
            ELSEIF( FGB.LT.ZERO .AND. K.LT.KFLD ) THEN
             IF( IXP(NT).GT.0 ) THEN
              FCGT = YG(NT,NSL)/(SG(2,NT)*PORD(2,NT)+SMALL)
              R = ((CO(N,NSL)*FCGP-CO(NT,NSL)*FCGT)
     &          /(CO(NB,NSL)*FCGB-CO(N,NSL)*FCGP+EPSL))
     &          *((DZGF(NB)+DZGF(N))/(DZGF(N)+DZGF(NT)))
              THETA = FLIMIT( R,CRGB,ISLC(1) )
              DZF = DZGF(N)/(DZGF(N)+DZGF(NB))
              WCZ = CO(NB,NSL)*FGB*THETA*DZF*FCGB
     &          + CO(N,NSL)*FGB*(1.D+0-THETA*DZF)*FCGP
              WC(NPZ,NSL) = WC(NPZ,NSL) + WCZ/AFZ(NPZ)
              CC = CC + WCZ
             ELSE
              WCZ = CO(N,NSL)*FGB*FCGP
              WC(NPZ,NSL) = WC(NPZ,NSL) + WCZ/AFZ(NPZ)
              CC = CC + WCZ
             ENDIF
            ELSEIF( FGB.LT.ZERO .AND. K.EQ.KFLD ) THEN
              WCZ = CO(N,NSL)*FGB*FCGP
              WC(NPZ,NSL) = WC(NPZ,NSL) + WCZ/AFZ(NPZ)
              CC = CC + WCZ
            ENDIF
            AB = DGZ*FCGB
            AP = DGZ*FCGP
!
!---      Patankar solute transport  ---
!
          ELSE
            AGB = MAX(FGB,ZERO)
     &        + DGZ*MAX((ONE-(TENTH*ABS(FGB)/(DGZ+SMALL)))**5,ZERO)
            AP = (AGB-FGB)*FCGP
            AB = AGB*FCGB
          ENDIF
          MCB = IXP(NB)
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
            MCOL = MCB
            MROW = MCP-MCB+MDT
            ALU(MRD,MCD) = ALU(MRD,MCD) + AP
            ALU(MROW,MCOL) = ALU(MROW,MCOL) - AB
!
!---      SPLib solver  ---
!
          ELSEIF( ILES.EQ.3 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AB

!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AB

          ENDIF
        ENDIF
  110   CONTINUE
!
!---    South face diffusion and advection terms  ---
!
        IF( J.NE.1 ) THEN
          IF( IXP(NS).EQ.0 .OR. INBS(2,N).GT.0 ) GOTO 120
          DYSR = RP(I)*DYGF(NS)
          TCOR = (T(2,NS)+TABS)/TSPRF
          PCOR = (PG(2,NS)+PATM)/PATM
          SDFGS = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCS = SG(2,NS)*PORD(2,NS)
          FCGS = YG(NS,NSL)/(VMCS+SMALL)
          DGS = TORG(2,NS)*(SG(2,NS)-SGT(2,NS))*PORD(2,NS)*SDFGS
          INDX = 16
          DGY = DIFMN(DGS,DGP,DYSR,DYR,VG(1,NPY),INDX)
          DGY = AFY(NPY)*(DGY+DPGS)/(DYGP(NPY)*RP(I))
          FGS = AFY(NPY)*VG(1,NPY)
          IF( MOD(ISLC(23),100)/10.EQ.1 ) FGS = 0.D+0
          VMCSX = DIFMN(VMCS,VMCP,DYSR,DYR,VG(1,NPY),INDX)
          CRGS = ABS(VG(1,NPY))*DT/(DYGP(NPY)*VMCSX+SMALL)/RP(I)
!
!---      TVD solute transport  ---
!
          IF( ISLC(1).GE.1 ) THEN
            IF( FGS.GE.ZERO .AND. J.GT.2 ) THEN
             IF( IXP(NS-IFLD).GT.0 ) THEN
              FCGSS = YG(NS-IFLD,NSL)/(SG(2,NS-IFLD)
     &          *PORD(2,NS-IFLD)+SMALL)
              R = ((CO(NS,NSL)*FCGS-CO(NS-IFLD,NSL)*FCGSS)
     &          /(CO(N,NSL)*FCGP-CO(NS,NSL)*FCGS+EPSL))
     &          *((DYGF(N)+DYGF(NS))/(DYGF(NS)+DYGF(NN-IFLD)))
              THETA = FLIMIT( R,CRGS,ISLC(1) )
              DYF = DYGF(NS)/(DYGF(N)+DYGF(NS))
              VCY = CO(NS,NSL)*FGS*(1.D+0-THETA*DYF)*FCGS
     &          + CO(N,NSL)*FGS*THETA*DYF*FCGP
              VC(NPY,NSL) = VC(NPY,NSL) + VCY/AFY(NPY)
              CC = CC + VCY
             ELSE
              VCY = CO(NS,NSL)*FGS*FCGS
              VC(NPY,NSL) = VC(NPY,NSL) + VCY/AFY(NPY)
              CC = CC + VCY
             ENDIF
            ELSEIF( FGS.GE.ZERO .AND. J.EQ.2 ) THEN
              VCY = CO(NS,NSL)*FGS*FCGS
              VC(NPY,NSL) = VC(NPY,NSL) + VCY/AFY(NPY)
              CC = CC + VCY
            ELSEIF( FGS.LT.ZERO .AND. J.LT.JFLD ) THEN
             IF( IXP(NN).GT.0 ) THEN
              FCGN = YG(NN,NSL)/(SG(2,NN)*PORD(2,NN)+SMALL)
              R = ((CO(N,NSL)*FCGP-CO(NN,NSL)*FCGN)
     &          /(CO(NS,NSL)*FCGS-CO(N,NSL)*FCGP+EPSL))
     &          *((DYGF(NS)+DYGF(N))/(DYGF(N)+DYGF(NN)))
              THETA = FLIMIT( R,CRGS,ISLC(1) )
              DYF = DYGF(N)/(DYGF(N)+DYGF(NS))
              VCY = CO(NS,NSL)*FGS*THETA*DYF*FCGS
     &          + CO(N,NSL)*FGS*(1.D+0-THETA*DYF)*FCGP
              VC(NPY,NSL) = VC(NPY,NSL) + VCY/AFY(NPY)
              CC = CC + VCY
             ELSE
              VCY = CO(N,NSL)*FGS*FCGP
              VC(NPY,NSL) = VC(NPY,NSL) + VCY/AFY(NPY)
              CC = CC + VCY
             ENDIF
            ELSEIF( FGS.LT.ZERO .AND. J.EQ.JFLD ) THEN
              VCY = CO(N,NSL)*FGS*FCGP
              VC(NPY,NSL) = VC(NPY,NSL) + VCY/AFY(NPY)
              CC = CC + VCY
            ENDIF
            AS = DGY*FCGS
            AP = DGY*FCGP
!
!---      Patankar solute transport  ---
!
          ELSE
            AGS = MAX(FGS,ZERO)
     &        + DGY*MAX((ONE-(TENTH*ABS(FGS)/(DGY+SMALL)))**5,ZERO)
            AP = (AGS-FGS)*FCGP
            AS = AGS*FCGS
          ENDIF
          MCS = IXP(NS)
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
            MCOL = MCS
            MROW = MCP-MCS+MDT
            ALU(MRD,MCD) = ALU(MRD,MCD) + AP
            ALU(MROW,MCOL) = ALU(MROW,MCOL) - AS
!
!---      SPLib solver  ---
!
          ELSEIF( ILES.EQ.3 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AS

!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AS

          ENDIF
        ENDIF
  120   CONTINUE
!
!---    West face diffusion and advection terms  ---
!
        IF( I.NE.1 ) THEN
          IF( IXP(NW).EQ.0 .OR. INBS(3,N).GT.0 ) GOTO 130
          TCOR = (T(2,NW)+TABS)/TSPRF
          PCOR = (PG(2,NW)+PATM)/PATM
          SDFGW = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCW = SG(2,NW)*PORD(2,NW)
          FCGW = YG(NW,NSL)/(VMCW+SMALL)
          DGW = TORG(2,NW)*(SG(2,NW)-SGT(2,NW))*PORD(2,NW)*SDFGW
          INDX = 16
          DGX = DIFMN(DGW,DGP,DXGF(NW),DXGF(N),UG(1,NPX),INDX)
          DGX = AFX(NPX)*(DGX+DPGW)/DXGP(NPX)
          FGW = AFX(NPX)*UG(1,NPX)
          IF( MOD(ISLC(23),100)/10.EQ.1 ) FGW = 0.D+0
          VMCWX = DIFMN(VMCW,VMCP,DXGF(NW),DXGF(N),UG(1,NPX),INDX)
          CRGW = ABS(UG(1,NPX))*DT/(DXGP(NPX)*VMCWX+SMALL)
!
!---      TVD solute transport  ---
!
          IF( ISLC(1).GE.1 ) THEN
            IF( FGW.GE.ZERO .AND. I.GT.2 ) THEN
             IF( IXP(NW-1).GT.0 ) THEN
              FCGWW = YG(NW-1,NSL)/(SG(2,NW-1)*PORD(2,NW-1)+SMALL)
              R = ((CO(NW,NSL)*FCGW-CO(NW-1,NSL)*FCGWW)
     &          /(CO(N,NSL)*FCGP-CO(NW,NSL)*FCGW+EPSL))
     &          *((DXGF(N)+DXGF(NW))/(DXGF(NW)+DXGF(NW-1)))
              THETA = FLIMIT( R,CRGW,ISLC(1) )
              DXF = DXGF(NW)/(DXGF(N)+DXGF(NW))
              UCX = CO(NW,NSL)*FGW*(1.D+0-THETA*DXF)*FCGW
     &          + CO(N,NSL)*FGW*THETA*DXF*FCGP
              UC(NPX,NSL) = UC(NPX,NSL) + UCX/AFX(NPX)
              CC = CC + UCX
             ELSE
              UCX = CO(NW,NSL)*FGW*FCGW
              UC(NPX,NSL) = UC(NPX,NSL) + UCX/AFX(NPX)
              CC = CC + UCX
             ENDIF
            ELSEIF( FGW.GE.ZERO .AND. I.EQ.2 ) THEN
              UCX = CO(NW,NSL)*FGW*FCGW
              UC(NPX,NSL) = UC(NPX,NSL) + UCX/AFX(NPX)
              CC = CC + UCX
            ELSEIF( FGW.LT.ZERO .AND. I.LT.IFLD ) THEN
             IF( IXP(NE).GT.0 ) THEN
              FCGE = YG(NE,NSL)/(SG(2,NE)*PORD(2,NE)+SMALL)
              R = ((CO(N,NSL)*FCGP-CO(NE,NSL)*FCGE)
     &          /(CO(NW,NSL)*FCGW-CO(N,NSL)*FCGP+EPSL))
     &          *((DXGF(NW)+DXGF(N))/(DXGF(N)+DXGF(NE)))
              THETA = FLIMIT( R,CRGW,ISLC(1) )
              DXF = DXGF(N)/(DXGF(N)+DXGF(NW))
              UCX = CO(NW,NSL)*FGW*THETA*DXF*FCGW
     &          + CO(N,NSL)*FGW*(1.D+0-THETA*DXF)*FCGP
              UC(NPX,NSL) = UC(NPX,NSL) + UCX/AFX(NPX)
              CC = CC + UCX
             ELSE
              UCX = CO(N,NSL)*FGW*FCGP
              UC(NPX,NSL) = UC(NPX,NSL) + UCX/AFX(NPX)
              CC = CC + UCX
             ENDIF
            ELSEIF( FGW.LT.ZERO .AND. I.EQ.IFLD ) THEN
              UCX = CO(N,NSL)*FGW*FCGP
              UC(NPX,NSL) = UC(NPX,NSL) + UCX/AFX(NPX)
              CC = CC + UCX
            ENDIF
            AW = DGX*FCGW
            AP = DGX*FCGP
!
!---      Patankar solute transport  ---
!
          ELSE
            AGW = MAX(FGW,ZERO)
     &        + DGX*MAX((ONE-(TENTH*ABS(FGW)/(DGX+SMALL)))**5,ZERO)
            AP = (AGW-FGW)*FCGP
            AW = AGW*FCGW
          ENDIF
          MCW = IXP(NW)
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
            MCOL = MCW
            MROW = MCP-MCW+MDT
            ALU(MRD,MCD) = ALU(MRD,MCD) + AP
            ALU(MROW,MCOL) = ALU(MROW,MCOL) - AW
!
!---      SPLib solver  ---
!
          ELSEIF( ILES.EQ.3 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AW

!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AW

          ENDIF
        ENDIF
  130   CONTINUE
!
!---    East face diffusion and advection terms  ---
!
        IF( I.NE.IFLD ) THEN
          IF( IXP(NE).EQ.0 .OR. INBS(4,N).GT.0 ) GOTO 140
          TCOR = (T(2,NE)+TABS)/TSPRF
          PCOR = (PG(2,NE)+PATM)/PATM
          SDFGE = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCE = SG(2,NE)*PORD(2,NE)
          FCGE = YG(NE,NSL)/(VMCE+SMALL)
          DGE = TORG(2,NE)*(SG(2,NE)-SGT(2,NE))*PORD(2,NE)*SDFGE
          INDX = 16
          DGX = DIFMN(DGP,DGE,DXGF(N),DXGF(NE),UG(1,NQX),INDX)
          DGX = AFX(NQX)*(DGX+DPGE)/DXGP(NQX)
          FGE = AFX(NQX)*UG(1,NQX)
          IF( MOD(ISLC(23),100)/10.EQ.1 ) FGE = 0.D+0
          VMCEX = DIFMN(VMCP,VMCE,DXGF(N),DXGF(NE),UG(1,NQX),INDX)
          CRGE = ABS(UG(1,NQX))*DT/(DXGP(NQX)*VMCEX+SMALL)
!
!---      TVD solute transport  ---
!
          IF( ISLC(1).GE.1 ) THEN
            IF( FGE.GE.ZERO .AND. I.GT.1 ) THEN
             IF( IXP(NW).GT.0 ) THEN
              FCGW = YG(NW,NSL)/(SG(2,NW)*PORD(2,NW)+SMALL)
              R = ((CO(N,NSL)*FCGP-CO(NW,NSL)*FCGW)
     &          /(CO(NE,NSL)*FCGE-CO(N,NSL)*FCGP+EPSL))
     &          *((DXGF(NE)+DXGF(N))/(DXGF(N)+DXGF(NW)))
              THETA = FLIMIT( R,CRGE,ISLC(1) )
              DXF = DXGF(N)/(DXGF(N)+DXGF(NE))
              UCX = CO(NE,NSL)*FGE*THETA*DXF*FCGE
     &          + CO(N,NSL)*FGE*(1.D+0-THETA*DXF)*FCGP
              CC = CC - UCX
             ELSE
              UCX = CO(N,NSL)*FGE*FCGP
              CC = CC - UCX
             ENDIF
            ELSEIF( FGE.GE.ZERO .AND. I.EQ.1 ) THEN
              UCX = CO(N,NSL)*FGE*FCGP
              CC = CC - UCX
            ELSEIF( FGE.LT.ZERO .AND. I.LT.IFLD-1 ) THEN
             IF( IXP(NE+1).GT.0 ) THEN
              FCGEE = YG(NE+1,NSL)/(SG(2,NE+1)*PORD(2,NE+1)+SMALL)
              R = ((CO(NE,NSL)*FCGE-CO(NE+1,NSL)*FCGEE)
     &          /(CO(N,NSL)*FCGP-CO(NE,NSL)*FCGE+EPSL))
     &          *((DXGF(N)+DXGF(NE))/(DXGF(NE)+DXGF(NE+1)))
              THETA = FLIMIT( R,CRGE,ISLC(1) )
              DXF = DXGF(NE)/(DXGF(N)+DXGF(NE))
              UCX = CO(NE,NSL)*FGE*(1.D+0-THETA*DXF)*FCGE
     &          + CO(N,NSL)*FGE*THETA*DXF*FCGP
              CC = CC - UCX
             ELSE
              UCX = CO(NE,NSL)*FGE*FCGE
              CC = CC - UCX
             ENDIF
            ELSEIF( FGE.LT.ZERO .AND. I.EQ.IFLD-1 ) THEN
              UCX = CO(NE,NSL)*FGE*FCGE
              CC = CC - UCX
            ENDIF
            AE = DGX*FCGE
            AP = DGX*FCGP
!
!---      Patankar solute transport  ---
!
          ELSE
            AGE = MAX(-FGE,ZERO)
     &        + DGX*MAX((ONE-(TENTH*ABS(FGE)/(DGX+SMALL)))**5,ZERO)
            AP = (AGE+FGE)*FCGP
            AE = AGE*FCGE
          ENDIF
          MCE = IXP(NE)
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
            MCOL = MCE
            MROW = MCP-MCE+MDT
            ALU(MRD,MCD) = ALU(MRD,MCD) + AP
            ALU(MROW,MCOL) = ALU(MROW,MCOL) - AE
!
!---      SPLib solver  ---
!
          ELSEIF( ILES.EQ.3 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AE

!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AE

          ENDIF
        ENDIF
  140   CONTINUE
!
!---    North face diffusion and advection terms  ---
!
        IF( J.NE.JFLD ) THEN
          IF( IXP(NN).EQ.0 .OR. INBS(5,N).GT.0 ) GOTO 150
          DYNR = RP(I)*DYGF(NN)
          TCOR = (T(2,NN)+TABS)/TSPRF
          PCOR = (PG(2,NN)+PATM)/PATM
          SDFGN = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCN = SG(2,NN)*PORD(2,NN)
          FCGN = YG(NN,NSL)/(VMCN+SMALL)
          DGN = TORG(2,NN)*(SG(2,NN)-SGT(2,NN))*PORD(2,NN)*SDFGN
          INDX = 16
          DGY = DIFMN(DGP,DGN,DYNR,DYR,VG(1,NQY),INDX)
          DGY = AFY(NQY)*(DGY+DPGN)/(DYGP(NQY)*RP(I))
          FGN = AFY(NQY)*VG(1,NQY)
          IF( MOD(ISLC(23),100)/10.EQ.1 ) FGN = 0.D+0
          VMCNX = DIFMN(VMCP,VMCN,DYNR,DYR,VG(1,NQY),INDX)
          CRGN = ABS(VG(1,NQY))*DT/(DYGP(NQY)*VMCNX+SMALL)/RP(I)
          IF( ISLC(1).GE.1 ) THEN
!
!---        TVD solute transport  ---
!
            IF( FGN.GE.ZERO .AND. J.GT.1 ) THEN
             IF( IXP(NS).GT.0 ) THEN
              FCGS = YG(NS,NSL)/(SG(2,NS)*PORD(2,NS)+SMALL)
              R = ((CO(N,NSL)*FCGP-CO(NS,NSL)*FCGS)
     &          /(CO(NN,NSL)*FCGN-CO(N,NSL)*FCGP+EPSL))
     &          *((DYGF(NN)+DYGF(N))/(DYGF(N)+DYGF(NS)))
              THETA = FLIMIT( R,CRGN,ISLC(1) )
              DYF = DYGF(N)/(DYGF(N)+DYGF(NN))
              VCY = CO(NN,NSL)*FGN*THETA*DYF*FCGN
     &          + CO(N,NSL)*FGN*(1.D+0-THETA*DYF)*FCGP
              CC = CC - VCY
             ELSE
              VCY = CO(N,NSL)*FGN*FCGP
              CC = CC - VCY
             ENDIF
            ELSEIF( FGN.GE.ZERO .AND. J.EQ.1 ) THEN
              VCY = CO(N,NSL)*FGN*FCGP
              CC = CC - VCY
            ELSEIF( FGN.LT.ZERO .AND. J.LT.JFLD-1 ) THEN
             IF( IXP(NN+IFLD).GT.0 ) THEN
              FCGNN = YG(NN+IFLD,NSL)/(SG(2,NN+IFLD)
     &          *PORD(2,NN+IFLD)+SMALL)
              R = ((CO(NN,NSL)*FCGN-CO(NN+IFLD,NSL)*FCGNN)
     &          /(CO(N,NSL)*FCGP-CO(NN,NSL)*FCGN+EPSL))
     &          *((DYGF(N)+DYGF(NN))/(DYGF(NN)+DYGF(NN+IFLD)))
              THETA = FLIMIT( R,CRGN,ISLC(1) )
              DYF = DYGF(NN)/(DYGF(N)+DYGF(NN))
              VCY = CO(NN,NSL)*FGN*(1.D+0-THETA*DYF)*FCGN
     &          + CO(N,NSL)*FGN*THETA*DYF*FCGP
              CC = CC - VCY
             ELSE
              VCY = CO(NN,NSL)*FGN*FCGN
              CC = CC - VCY
             ENDIF
            ELSEIF( FGN.LT.ZERO .AND. J.EQ.JFLD-1 ) THEN
              VCY = CO(NN,NSL)*FGN*FCGN
              CC = CC - VCY
            ENDIF
            AN = DGY*FCGN
            AP = DGY*FCGP
          ELSE
!
!---        Patankar solute transport  ---
!
            AGN = MAX(-FGN,ZERO)
     &        + DGY*MAX((ONE-(TENTH*ABS(FGN)/(DGY+SMALL)))**5,ZERO)
            AP = (AGN+FGN)*FCGP
            AN = AGN*FCGN
          ENDIF
          MCN = IXP(NN)
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
            MCOL = MCN
            MROW = MCP-MCN+MDT
            ALU(MRD,MCD) = ALU(MRD,MCD) + AP
            ALU(MROW,MCOL) = ALU(MROW,MCOL) - AN
!
!---      SPLib solver  ---
!
          ELSEIF( ILES.EQ.3 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AN

!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AN

          ENDIF
        ENDIF
  150   CONTINUE
!
!---    Top face diffusion and advection terms  ---
!
        IF( K.NE.KFLD ) THEN
          IF( IXP(NT).EQ.0 .OR. INBS(6,N).GT.0 ) GOTO 160
          TCOR = (T(2,NT)+TABS)/TSPRF
          PCOR = (PG(2,NT)+PATM)/PATM
          SDFGT = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCT = SG(2,NT)*PORD(2,NT)
          FCGT = YG(NT,NSL)/(VMCT+SMALL)
          DGT = TORG(2,NT)*(SG(2,NT)-SGT(2,NT))*PORD(2,NT)*SDFGT
          INDX = 16
          DGZ = DIFMN(DGP,DGT,DZGF(N),DZGF(NT),WG(1,NQZ),INDX)
          DGZ = AFZ(NQZ)*(DGZ+DPGT)/DZGP(NQZ)
          FGT = AFZ(NQZ)*WG(1,NQZ)
          IF( MOD(ISLC(23),100)/10.EQ.1 ) FGT = 0.D+0
          VMCTX = DIFMN(VMCP,VMCT,DZGF(N),DZGF(NT),WG(1,NQZ),INDX)
          CRGT = ABS(WG(1,NQZ))*DT/(DZGP(NQZ)*VMCTX+SMALL)
          IF( ISLC(1).GE.1 ) THEN
!
!---        TVD solute transport  ---
!
            IF( FGT.GE.ZERO .AND. K.GT.1 ) THEN
             IF( IXP(NB).GT.0 ) THEN
              FCGB = YG(NB,NSL)/(SG(2,NB)*PORD(2,NB)+SMALL)
              R = ((CO(N,NSL)*FCGP-CO(NB,NSL)*FCGB)
     &          /(CO(NT,NSL)*FCGT-CO(N,NSL)*FCGP+EPSL))
     &          *((DZGF(NT)+DZGF(N))/(DZGF(N)+DZGF(NB)))
              THETA = FLIMIT( R,CRGT,ISLC(1) )
              DZF = DZGF(N)/(DZGF(N)+DZGF(NT))
              WCZ = CO(NT,NSL)*FGT*THETA*DZF*FCGT
     &          + CO(N,NSL)*FGT*(1.D+0-THETA*DZF)*FCGP
              CC = CC - WCZ
             ELSE
              WCZ = CO(N,NSL)*FGT*FCGP
              CC = CC - WCZ
             ENDIF
            ELSEIF( FGT.GE.ZERO .AND. K.EQ.1 ) THEN
              WCZ = CO(N,NSL)*FGT*FCGP
              CC = CC - WCZ
            ELSEIF( FGT.LT.ZERO .AND. K.LT.KFLD-1 ) THEN
             IF( IXP(NT+IJFLD).GT.0 ) THEN
              FCGTT = YG(NT+IJFLD,NSL)/(SG(2,NT+IJFLD)
     &          *PORD(2,NT+IJFLD)+SMALL)
              R = ((CO(NT,NSL)*FCGT-CO(NT+IJFLD,NSL)*FCGTT)
     &          /(CO(N,NSL)*FCGP-CO(NT,NSL)*FCGT+EPSL))
     &          *((DZGF(N)+DZGF(NT))/(DZGF(NT)+DZGF(NT+IJFLD)))
              THETA = FLIMIT( R,CRGT,ISLC(1) )
              DZF = DZGF(NT)/(DZGF(N)+DZGF(NT))
              WCZ = CO(NT,NSL)*FGT*(1.D+0-THETA*DZF)*FCGT
     &          + CO(N,NSL)*FGT*THETA*DZF*FCGP
              CC = CC - WCZ
             ELSE
              WCZ = CO(NT,NSL)*FGT*FCGT
              CC = CC - WCZ
             ENDIF
            ELSEIF( FGT.LT.ZERO .AND. K.EQ.KFLD-1 ) THEN
              WCZ = CO(NT,NSL)*FGT*FCGT
              CC = CC - WCZ
            ENDIF
            AT = DGZ*FCGT
            AP = DGZ*FCGP
!
!---      Patankar solute transport  ---
!
          ELSE
            AGT = MAX(-FGT,ZERO)
     &        + DGZ*MAX((ONE-(TENTH*ABS(FGT)/(DGZ+SMALL)))**5,ZERO)
            AP = (AGT+FGT)*FCGP
            AT = AGT*FCGT
          ENDIF
          MCT = IXP(NT)
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
            MCOL = MCT
            MROW = MCP-MCT+MDT
            ALU(MRD,MCD) = ALU(MRD,MCD) + AP
            ALU(MROW,MCOL) = ALU(MROW,MCOL) - AT
!
!---      SPLib solver  ---
!
          ELSEIF( ILES.EQ.3 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AT

!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AT

          ENDIF
        ENDIF
  160   CONTINUE
!
!---    Solution vector  ---
!
        IF( IEQW.EQ.0 ) THEN
!
!---      Solution vector and parent 1 decay  ---
!
          BLU(MCP) = BLU(MCP) + CO(N,NSL)*SC
!     &      - 6.931D-1*CO(N,NSL)*VOL(N)/HLF(NSL)
!!          BLU(MCP) = BLU(MCP) + CO(N,NSL)*(2.D+0**(-DT/HLF(NSL)))*SC
!!
!!---      Daughter 1 chain decay  ---
!!
!          IF( NSL.LE.NSOLU ) THEN
!            DO 170 NPSL = 1,NSOLU
!              IF( NPSL.EQ.NSL ) GOTO 170
!              BLU(MCP) = BLU(MCP) +
!     &          CHDF(NPSL,NSL)*6.931D-1*CO(N,NPSL)*VOL(N)/HLF(NPSL)
!!              BLU(MCP) = BLU(MCP) +
!!     &          CHDF(NPSL,NSL)*CO(N,NPSL)*(2.D+0**(-DT/HLF(NPSL)))*SC
!  170       CONTINUE
!          ENDIF
        ENDIF
        IF( ISLC(1).GE.1 ) BLU(MCP) = BLU(MCP) + CC
  900 CONTINUE
!
!---  End of SJCBG group  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END



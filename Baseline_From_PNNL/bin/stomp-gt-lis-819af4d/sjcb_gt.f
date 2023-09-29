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
      DO N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) CYCLE
        I = ID(N)
        IZN = IZ(N)
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
!---    SPLib  or LIS solver  ---
!
        ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
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
!---    Bottom face diffusion and advection terms  ---
!
        DO NC = 1,4
          NB = ICM(NC,1,N)
          IF( NB.EQ.0 ) EXIT
!
!---      Multiple bottom connections  ---
!
          IF( INBS(1,N).LT.0 ) THEN
            NPZ = NSSZ(NB)
!
!---      Single bottom connections  ---
!
          ELSE
            NPZ = NSZ(N)
          ENDIF
          TCOR = (T(2,NB)+TABS)/TSPRF
          PCOR = (PG(2,NB)+PATM)/PATM
          SDFGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCB = SG(2,NB)*PORD(2,NB)
          IF( IDISP.EQ.1 ) THEN
!
!---        Hydrodynamic dispersion coefficient  ---
!
            UBX = 5.D-1*(UG(1,NSX(N))+UG(1,NSX(N)+1))
            VBX = 5.D-1*(VG(1,NSY(N))+VG(1,NSY(N)+IFLD))
            WBX = WG(1,NSZ(N))
            INDX = 17
            VMCBX = DIFMN(VMCB,VMCP,DZGF(NB),DZGF(N),WBX,INDX)
!
!---        Hydrodynamic dispersion coefficient w/ flow velocity  ---
!
            IF( ISLC(52).EQ.0 ) THEN
              UBX = UBX/VMCBX
              VBX = VBX/VMCBX
              WBX = WBX/VMCBX
            ENDIF
            UBSQ = UBX*UBX
            VBSQ = VBX*VBX
            WBSQ = WBX*WBX
            ZVB = SQRT(UBSQ+VBSQ+WBSQ)
            INDX = 17
            DISPLB = DISPL(IZ(NB))*SMDEF(IZ(NB),NSL)
            DISPLP = DISPL(IZ(N))*SMDEF(IZ(N),NSL)
            DISPTB = DISPT(IZ(NB))*SMDEF(IZ(NB),NSL)
            DISPTP = DISPT(IZ(N))*SMDEF(IZ(N),NSL)
            DPLB = DIFMN(DISPLB,DISPLP,DZGF(NB),DZGF(N),WBX,INDX)
            DPTB = DIFMN(DISPTB,DISPTP,DZGF(NB),DZGF(N),WBX,INDX)
            DPB = (DPLB*WBSQ + DPTB*(VBSQ+UBSQ))/(ZVB+SMALL)
          ELSE
            DPB = 0.D+0
          ENDIF
          FCGB = YG(NB,NSL)/(VMCB+SMALL)
          DGB = TORG(2,NB)*(SG(2,NB)-SGT(2,NB))*PORD(2,NB)*SDFGB
          INDX = 16
          DGZ = DIFMN(DGB,DGP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX)
          DGZ = AFZ(NPZ)*(DGZ+DPB)/DZGP(NPZ)
          FGB = AFZ(NPZ)*WG(1,NPZ)
          IF( MOD(ISLC(23),100)/10.EQ.1 ) FGB = 0.D+0
!
!---      Patankar solute transport  ---
!
          AGB = MAX(FGB,ZERO)
     &      + DGZ*MAX((ONE-(TENTH*ABS(FGB)/(DGZ+SMALL)))**5,ZERO)
          AP = (AGB-FGB)*FCGP
          AB = AGB*FCGB
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
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCD) = DLU(MCD) + AP
!
!---        Block refinement scheme  ---
!
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AB

          ENDIF
        ENDDO
!
!---    South face diffusion and advection terms  ---
!
        DO NC = 1,4
          NS = ICM(NC,2,N)
          IF( NS.EQ.0 ) EXIT
!
!---      Multiple south connections  ---
!
          IF( INBS(2,N).LT.0 ) THEN
            NPY = NSSY(NS)
!
!---      Single south connections  ---
!
          ELSE
            NPY = NSY(N)
          ENDIF
          TCOR = (T(2,NS)+TABS)/TSPRF
          PCOR = (PG(2,NS)+PATM)/PATM
          SDFGS = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCS = SG(2,NS)*PORD(2,NS)
          IF( IDISP.EQ.1 ) THEN
!
!---        Hydrodynamic dispersion coefficient  ---
!
            USX = 5.D-1*(UG(1,NSX(N))+UG(1,NSX(N)+1))
            WSX = 5.D-1*(WG(1,NSZ(N))+WG(1,NSZ(N)+IJFLD))
            VSX = VG(1,NSY(N))
            INDX = 17
            VMCSX = DIFMN(VMCS,VMCP,DYGF(NS),DYGF(N),VSX,INDX)
!
!---        Hydrodynamic dispersion coefficient w/ flow velocity  ---
!
            IF( ISLC(52).EQ.0 ) THEN
              USX = USX/VMCSX
              VSX = VSX/VMCSX
              WSX = WSX/VMCSX
            ENDIF
            USSQ = USX*USX
            VSSQ = VSX*VSX
            WSSQ = WSX*WSX
            ZVS = SQRT(USSQ+VSSQ+WSSQ)
            INDX = 17
            DISPLS = DISPL(IZ(NS))*SMDEF(IZ(NS),NSL)
            DISPLP = DISPL(IZ(N))*SMDEF(IZ(N),NSL)
            DISPTS = DISPT(IZ(NS))*SMDEF(IZ(NS),NSL)
            DISPTP = DISPT(IZ(N))*SMDEF(IZ(N),NSL)
            DPLS = DIFMN(DISPLS,DISPLP,DYGF(NS),DYGF(N),VSX,INDX)
            DPTS = DIFMN(DISPTS,DISPTP,DYGF(NS),DYGF(N),VSX,INDX)
            DPS = (DPLS*VSSQ + DPTS*(USSQ+WSSQ))/(ZVS+SMALL)
          ELSE
            DPS = 0.D+0
          ENDIF
          DYSR = RP(I)*DYGF(NS)
          FCGS = YG(NS,NSL)/(VMCS+SMALL)
          DGS = TORG(2,NS)*(SG(2,NS)-SGT(2,NS))*PORD(2,NS)*SDFGS
          INDX = 16
          DGY = DIFMN(DGS,DGP,DYSR,DYR,VG(1,NPY),INDX)
          DGY = AFY(NPY)*(DGY+DPS)/(DYGP(NPY)*RP(I))
          FGS = AFY(NPY)*VG(1,NPY)
          IF( MOD(ISLC(23),100)/10.EQ.1 ) FGS = 0.D+0
!
!---      Patankar solute transport  ---
!
          AGS = MAX(FGS,ZERO)
     &      + DGY*MAX((ONE-(TENTH*ABS(FGS)/(DGY+SMALL)))**5,ZERO)
          AP = (AGS-FGS)*FCGP
          AS = AGS*FCGS
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
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCD) = DLU(MCD) + AP
!
!---        Block refinement scheme  ---
!
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AS

          ENDIF
        ENDDO
!
!---    West face diffusion and advection terms  ---
!
        DO NC = 1,4
          NW = ICM(NC,3,N)
          IF( NW.EQ.0 ) EXIT
!
!---      Multiple west connections  ---
!
          IF( INBS(3,N).LT.0 ) THEN
            NPX = NSSX(NW)
!
!---      Single west connections  ---
!
          ELSE
            NPX = NSX(N)
          ENDIF
          TCOR = (T(2,NW)+TABS)/TSPRF
          PCOR = (PG(2,NW)+PATM)/PATM
          SDFGW = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCW = SG(2,NW)*PORD(2,NW)
          IF( IDISP.EQ.1 ) THEN
!
!---        Hydrodynamic dispersion coefficient  ---
!
            UWX = UG(1,NSX(N))
            VWX = 5.D-1*(VG(1,NSY(N))+VG(1,NSY(N)+IFLD))
            WWX = 5.D-1*(WG(1,NSZ(N))+WG(1,NSZ(N)+IJFLD))
            INDX = 17
            VMCWX = DIFMN(VMCW,VMCP,DXGF(NW),DXGF(N),UWX,INDX)
!
!---        Hydrodynamic dispersion coefficient w/ flow velocity  ---
!
            IF( ISLC(52).EQ.0 ) THEN
              UWX = UWX/VMCWX
              VWX = VWX/VMCWX
              WWX = WWX/VMCWX
            ENDIF
            UWSQ = UWX*UWX
            VWSQ = VWX*VWX
            WWSQ = WWX*WWX
            ZVW = SQRT(UWSQ+VWSQ+WWSQ)
            INDX = 17
            DISPLW = DISPL(IZ(NW))*SMDEF(IZ(NW),NSL)
            DISPLP = DISPL(IZ(N))*SMDEF(IZ(N),NSL)
            DISPTW = DISPT(IZ(NW))*SMDEF(IZ(NW),NSL)
            DISPTP = DISPT(IZ(N))*SMDEF(IZ(N),NSL)
            DPLW = DIFMN(DISPLW,DISPLP,DXGF(NW),DXGF(N),UWX,INDX)
            DPTW = DIFMN(DISPTW,DISPTP,DXGF(NW),DXGF(N),UWX,INDX)
            DPW = (DPLW*UWSQ + DPTW*(VWSQ+WWSQ))/(ZVW+SMALL)
          ELSE
            DPW = 0.D+0
          ENDIF
          FCGW = YG(NW,NSL)/(VMCW+SMALL)
          DGW = TORG(2,NW)*(SG(2,NW)-SGT(2,NW))*PORD(2,NW)*SDFGW
          INDX = 16
          DGX = DIFMN(DGW,DGP,DXGF(NW),DXGF(N),UG(1,NPX),INDX)
          DGX = AFX(NPX)*(DGX+DPW)/DXGP(NPX)
          FGW = AFX(NPX)*UG(1,NPX)
          IF( MOD(ISLC(23),100)/10.EQ.1 ) FGW = 0.D+0
!
!---      Patankar solute transport  ---
!
          AGW = MAX(FGW,ZERO)
     &      + DGX*MAX((ONE-(TENTH*ABS(FGW)/(DGX+SMALL)))**5,ZERO)
          AP = (AGW-FGW)*FCGP
          AW = AGW*FCGW
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
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCD) = DLU(MCD) + AP
!
!---        Block refinement scheme  ---
!
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AW

          ENDIF
        ENDDO
!
!---    East face diffusion and advection terms  ---
!
        DO NC = 1,4
          NE = ICM(NC,4,N)
          IF( NE.EQ.0 ) EXIT
!
!---      Multiple west connections for east node  ---
!
          IF( INBS(3,NE).LT.0 ) THEN
            NQX = NSSX(N)
!
!---      Single west connection for east node  ---
!
          ELSE
            NQX = NSX(NE)
          ENDIF
          TCOR = (T(2,NE)+TABS)/TSPRF
          PCOR = (PG(2,NE)+PATM)/PATM
          SDFGE = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCE = SG(2,NE)*PORD(2,NE)
          IF( IDISP.EQ.1 ) THEN
!
!---        Hydrodynamic dispersion coefficient  ---
!
            UEX = UG(1,NSX(N)+1)
            VEX = 5.D-1*(VG(1,NSY(N))+VG(1,NSY(N)+IFLD))
            WEX = 5.D-1*(WG(1,NSZ(N))+WG(1,NSZ(N)+IJFLD))
            INDX = 17
            VMCEX = DIFMN(VMCP,VMCE,DXGF(N),DXGF(NE),UEX,INDX)
!
!---        Hydrodynamic dispersion coefficient w/ flow velocity  ---
!
            IF( ISLC(52).EQ.0 ) THEN
              UEX = UEX/VMCEX
              VEX = VEX/VMCEX
              WEX = WEX/VMCEX
            ENDIF
            UESQ = UEX*UEX
            VESQ = VEX*VEX
            WESQ = WEX*WEX
            ZVE = SQRT(UESQ+VESQ+WESQ)
            INDX = 17
            DISPLE = DISPL(IZ(NE))*SMDEF(IZ(NE),NSL)
            DISPLP = DISPL(IZ(N))*SMDEF(IZ(N),NSL)
            DISPTE = DISPT(IZ(NE))*SMDEF(IZ(NE),NSL)
            DISPTP = DISPT(IZ(N))*SMDEF(IZ(N),NSL)
            DPLE = DIFMN(DISPLP,DISPLE,DXGF(N),DXGF(NE),UEX,INDX)
            DPTE = DIFMN(DISPTP,DISPTE,DXGF(N),DXGF(NE),UEX,INDX)
            DPE = (DPLE*UESQ + DPTE*(VESQ+WESQ))/(ZVE+SMALL)
          ELSE
            DPE = 0.D+0
          ENDIF
          FCGE = YG(NE,NSL)/(VMCE+SMALL)
          DGE = TORG(2,NE)*(SG(2,NE)-SGT(2,NE))*PORD(2,NE)*SDFGE
          INDX = 16
          DGX = DIFMN(DGP,DGE,DXGF(N),DXGF(NE),UG(1,NQX),INDX)
          DGX = AFX(NQX)*(DGX+DPE)/DXGP(NQX)
          FGE = AFX(NQX)*UG(1,NQX)
          IF( MOD(ISLC(23),100)/10.EQ.1 ) FGE = 0.D+0
!
!---      Patankar solute transport  ---
!
          AGE = MAX(-FGE,ZERO)
     &      + DGX*MAX((ONE-(TENTH*ABS(FGE)/(DGX+SMALL)))**5,ZERO)
          AP = (AGE+FGE)*FCGP
          AE = AGE*FCGE
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
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCD) = DLU(MCD) + AP
!
!---        Block refinement scheme  ---
!
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AE

          ENDIF
        ENDDO
!
!---    North face diffusion and advection terms  ---
!
        DO NC = 1,4
          NN = ICM(NC,5,N)
          IF( NN.EQ.0 ) EXIT
!
!---      Multiple south connections for north node  ---
!
          IF( INBS(2,NN).LT.0 ) THEN
            NQY = NSSY(N)
!
!---      Single south connection for north node  ---
!
          ELSE
            NQY = NSY(NN)
          ENDIF
          TCOR = (T(2,NN)+TABS)/TSPRF
          PCOR = (PG(2,NN)+PATM)/PATM
          SDFGN = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCN = SG(2,NN)*PORD(2,NN)
          IF( IDISP.EQ.1 ) THEN
!
!---        Hydrodynamic dispersion coefficient  ---
!
            UNX = 5.D-1*(UG(1,NSX(N))+UG(1,NSX(N)+1))
            WNX = 5.D-1*(WG(1,NSZ(N))+WG(1,NSZ(N)+IJFLD))
            VNX = VG(1,NSY(N)+IFLD)
            INDX = 17
            VMCNX = DIFMN(VMCP,VMCN,DYNR,DYR,VNX,INDX)
!
!---        Hydrodynamic dispersion coefficient w/ flow velocity  ---
!
            IF( ISLC(52).EQ.0 ) THEN
              UNX = UNX/VMCNX
              VNX = VNX/VMCNX
              WNX = WNX/VMCNX
            ENDIF
            UNSQ = UNX*UNX
            VNSQ = VNX*VNX
            WNSQ = WNX*WNX
            ZVN = SQRT(UNSQ+VNSQ+WNSQ)
            INDX = 17
            DISPLN = DISPL(IZ(NN))*SMDEF(IZ(NN),NSL)
            DISPLP = DISPL(IZ(N))*SMDEF(IZ(N),NSL)
            DISPTN = DISPT(IZ(NN))*SMDEF(IZ(NN),NSL)
            DISPTP = DISPT(IZ(N))*SMDEF(IZ(N),NSL)
            DPLN = DIFMN(DISPLP,DISPLN,DYGF(N),DYGF(NN),VNX,INDX)
            DPTN = DIFMN(DISPTP,DISPTN,DYGF(N),DYGF(NN),VNX,INDX)
            DPN = (DPLN*VNSQ + DPTN*(UNSQ+WNSQ))/(ZVN+SMALL)
          ELSE
            DPN = 0.D+0
          ENDIF
          DYNR = RP(I)*DYGF(NN)
          FCGN = YG(NN,NSL)/(VMCN+SMALL)
          DGN = TORG(2,NN)*(SG(2,NN)-SGT(2,NN))*PORD(2,NN)*SDFGN
          INDX = 16
          DGY = DIFMN(DGP,DGN,DYNR,DYR,VG(1,NQY),INDX)
          DGY = AFY(NQY)*(DGY+DPN)/(DYGP(NQY)*RP(I))
          FGN = AFY(NQY)*VG(1,NQY)
          IF( MOD(ISLC(23),100)/10.EQ.1 ) FGN = 0.D+0
!
!---      Patankar solute transport  ---
!
          AGN = MAX(-FGN,ZERO)
     &      + DGY*MAX((ONE-(TENTH*ABS(FGN)/(DGY+SMALL)))**5,ZERO)
          AP = (AGN+FGN)*FCGP
          AN = AGN*FCGN
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
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCD) = DLU(MCD) + AP
!
!---        Block refinement scheme  ---
!
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AN

          ENDIF
        ENDDO
!
!---    Top face diffusion and advection terms  ---
!
        DO NC = 1,4
          NT = ICM(NC,6,N)
          IF( NT.EQ.0 ) EXIT
!
!---      Multiple bottom connections for top node  ---
!
          IF( INBS(1,NT).LT.0 ) THEN
            NQZ = NSSZ(N)
!
!---      Single bottom connection for top node  ---
!
          ELSE
            NQZ = NSZ(NT)
          ENDIF
          TCOR = (T(2,NT)+TABS)/TSPRF
          PCOR = (PG(2,NT)+PATM)/PATM
          SDFGT = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          VMCT = SG(2,NT)*PORD(2,NT)
          IF( IDISP.EQ.1 ) THEN
!
!---        Hydrodynamic dispersion coefficient  ---
!
            UTX = 5.D-1*(UG(1,NSX(N))+UG(1,NSX(N)+1))
            VTX = 5.D-1*(VG(1,NSY(N))+VG(1,NSY(N)+IFLD))
            WTX = WG(1,NSZ(N)+IJFLD)
            INDX = 17
            VMCTX = DIFMN(VMCP,VMCT,DZGF(N),DZGF(NT),WTX,INDX)
!
!---        Hydrodynamic dispersion coefficient w/ flow velocity  ---
!
            IF( ISLC(52).EQ.0 ) THEN
              UTX = UTX/VMCTX
              VTX = VTX/VMCTX
              WTX = WTX/VMCTX
            ENDIF
            UTSQ = UTX*UTX
            VTSQ = VTX*VTX
            WTSQ = WTX*WTX
            ZVT = SQRT(UTSQ+VTSQ+WTSQ)
            INDX = 17
            DISPLT = DISPL(IZ(NT))*SMDEF(IZ(NT),NSL)
            DISPLP = DISPL(IZ(N))*SMDEF(IZ(N),NSL)
            DISPTT = DISPT(IZ(NT))*SMDEF(IZ(NT),NSL)
            DISPTP = DISPT(IZ(N))*SMDEF(IZ(N),NSL)
            DPLT = DIFMN(DISPLP,DISPLT,DZGF(N),DZGF(NT),WTX,INDX)
            DPTT = DIFMN(DISPTP,DISPTT,DZGF(N),DZGF(NT),WTX,INDX)
            DPT = (DPLT*WTSQ + DPTT*(VTSQ+UTSQ))/(ZVT+SMALL)
          ELSE
            DPT = 0.D+0
          ENDIF
          FCGT = YG(NT,NSL)/(VMCT+SMALL)
          DGT = TORG(2,NT)*(SG(2,NT)-SGT(2,NT))*PORD(2,NT)*SDFGT
          INDX = 16
          DGZ = DIFMN(DGP,DGT,DZGF(N),DZGF(NT),WG(1,NQZ),INDX)
          DGZ = AFZ(NQZ)*(DGZ+DPT)/DZGP(NQZ)
          FGT = AFZ(NQZ)*WG(1,NQZ)
          IF( MOD(ISLC(23),100)/10.EQ.1 ) FGT = 0.D+0
!
!---      Patankar solute transport  ---
!
          AGT = MAX(-FGT,ZERO)
     &      + DGZ*MAX((ONE-(TENTH*ABS(FGT)/(DGZ+SMALL)))**5,ZERO)
          AP = (AGT+FGT)*FCGP
          AT = AGT*FCGT
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
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCD) = DLU(MCD) + AP
!
!---        Block refinement scheme  ---
!
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AT

          ENDIF
        ENDDO
!
!---    Solution vector  ---
!
        IF( IEQW.EQ.0 ) THEN
!!
!!---      Schroth half-life modification model  ---
!!
!          IF( PCLN(3,NSL).LT.0.D+0 ) THEN
!            IF( SL(2,N).LT.PCLN(1,NSL) ) THEN
!              HLFX = 1.D+12
!            ELSE
!              HLFX = HLF(NSL)*(PCLN(2,N) + SL(2,N) - 2.D+0*PCLN(1,N))/
!     &          (SL(2,N) - PCLN(1,N) + SMALL)
!            ENDIF
!          ELSE
!            HLFX = HLF(NSL)
!          ENDIF
!
!---      Solution vector and parent 1 decay  ---
!
          BLU(MCP) = BLU(MCP) + CO(N,NSL)*SC
!     &      - 6.931D-1*CO(N,NSL)*VOL(N)/HLFX
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
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SJCBG group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SJCBL( NSL )
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
!     aqueous-phase convective-dispersive mass transport equation.
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
      SUB_LOG(ISUB_LOG) = '/SJCBL'
!
!---  Fill matrix elements  ---
!
      DO N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) CYCLE
        I = ID(N)
        IZN = IZ(N)
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
          ALU(MRD,MCD) = ALU(MRD,MCD) + SC
!
!---    SPLib solver  ---
!
        ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
          MCD = KLUC(MCP,MA)
          MA = MA + 1
          DLU(MCD) = DLU(MCD) + SC

        ENDIF
!
!---    Molecular diffusion coefficients at the nodes  ---
!
        TCOR = (T(2,N)+TABS)/TSPRF
        SDFLP = SMDL(NSL)*TCOR*(VISRL/VISL(2,N))
        VMCP = SL(2,N)*PORD(2,N)
        FCLP = YL(N,NSL)/(VMCP+SMALL)
        IF( IEDL(NSL).EQ.2 ) THEN
          DLP = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &      EXP(VMCP*SDCL(3,IZN,NSL))
        ELSEIF( IEDL(NSL).EQ.3 ) THEN
          DLP = TORL(2,N)*VMCP*SMDL(NSL)
        ELSEIF( IEDL(NSL).EQ.4 ) THEN
          DLP = SDCL(1,IZN,NSL)*SDCL(2,IZN,NSL)*
     &      VMCP**SDCL(3,IZN,NSL)
        ELSE
          DLP = TORL(2,N)*VMCP*SDFLP
        ENDIF
!
!---    Bottom face diffusion and advection terms  ---
!
        DO NC = 1,4
          NB = ICM(NC,1,N)
          IF( NB.EQ.0 ) EXIT
!
!---      Multiple bottom connections  ---
!
          IF( INBS(1,N).LT.0 ) THEN
            NPZ = NSSZ(NB)
!
!---      Single bottom connections  ---
!
          ELSE
            NPZ = NSZ(N)
          ENDIF
          TCOR = (T(2,NB)+TABS)/TSPRF
          SDFLB = SMDL(NSL)*TCOR*(VISRL/VISL(2,NB))
          VMCB = SL(2,NB)*PORD(2,NB)
          IF( IDISP.EQ.1 ) THEN
!
!---        Hydrodynamic dispersion coefficient  ---
!
            UBX = 5.D-1*(UL(1,NSX(N))+UL(1,NSX(N)+1))
            VBX = 5.D-1*(VL(1,NSY(N))+VL(1,NSY(N)+IFLD))
            WBX = WL(1,NSZ(N))
            INDX = 17
            VMCBX = DIFMN(VMCB,VMCP,DZGF(NB),DZGF(N),WBX,INDX)
!
!---        Hydrodynamic dispersion coefficient w/ flow velocity  ---
!
            IF( ISLC(52).EQ.0 ) THEN
              UBX = UBX/VMCBX
              VBX = VBX/VMCBX
              WBX = WBX/VMCBX
            ENDIF
            UBSQ = UBX*UBX
            VBSQ = VBX*VBX
            WBSQ = WBX*WBX
            ZVB = SQRT(UBSQ+VBSQ+WBSQ)
            INDX = 17
            DISPLB = DISPL(IZ(NB))*SMDEF(IZ(NB),NSL)
            DISPLP = DISPL(IZ(N))*SMDEF(IZ(N),NSL)
            DISPTB = DISPT(IZ(NB))*SMDEF(IZ(NB),NSL)
            DISPTP = DISPT(IZ(N))*SMDEF(IZ(N),NSL)
            DPLB = DIFMN(DISPLB,DISPLP,DZGF(NB),DZGF(N),WBX,INDX)
            DPTB = DIFMN(DISPTB,DISPTP,DZGF(NB),DZGF(N),WBX,INDX)
            DPB = (DPLB*WBSQ + DPTB*(VBSQ+UBSQ))/(ZVB+SMALL)
          ELSE
            DPB = 0.D+0
          ENDIF
          FCLB = YL(NB,NSL)/(VMCB+SMALL)
          IF( IEDL(NSL).EQ.2 ) THEN
            DLB = SDCL(1,IZ(NB),NSL)*SDCL(2,IZ(NB),NSL)*
     &        EXP(VMCB*SDCL(3,IZ(NB),NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DLB = TORL(2,NB)*VMCB*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DLB = SDCL(1,IZ(NB),NSL)*SDCL(2,IZ(NB),NSL)*
     &        VMCB**SDCL(3,IZ(NB),NSL)
          ELSE
            DLB = TORL(2,NB)*VMCB*SDFLB
          ENDIF
          INDX = 16
          DLZ = DIFMN(DLB,DLP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX)
          DLZ = AFZ(NPZ)*(DLZ+DPB)/DZGP(NPZ)
          FLB = AFZ(NPZ)*WL(1,NPZ)
          IF( MOD(ISLC(23),10).EQ.1 ) FLB = 0.D+0
!
!---      Patankar solute transport  ---
!
          ALB = MAX(FLB,ZERO)
     &      + DLZ*MAX((ONE-(TENTH*ABS(FLB)/(DLZ+SMALL)))**5,ZERO)
          AP = (ALB-FLB)*FCLP
          AB = ALB*FCLB
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
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AB

          ENDIF
        ENDDO
!
!---    South face diffusion and advection terms  ---
!
        DO NC = 1,4
          NS = ICM(NC,2,N)
          IF( NS.EQ.0 ) EXIT
!
!---      Multiple south connections  ---
!
          IF( INBS(2,N).LT.0 ) THEN
            NPY = NSSY(NS)
!
!---      Single south connections  ---
!
          ELSE
            NPY = NSY(N)
          ENDIF
          TCOR = (T(2,NS)+TABS)/TSPRF
          SDFLS = SMDL(NSL)*TCOR*(VISRL/VISL(2,NS))
          VMCS = SL(2,NS)*PORD(2,NS)
          IF( IDISP.EQ.1 ) THEN
!
!---        Hydrodynamic dispersion coefficient  ---
!
            USX = 5.D-1*(UL(1,NSX(N))+UL(1,NSX(N)+1))
            WSX = 5.D-1*(WL(1,NSZ(N))+WL(1,NSZ(N)+IJFLD))
            VSX = VL(1,NSY(N))
            INDX = 17
            VMCSX = DIFMN(VMCS,VMCP,DYGF(NS),DYGF(N),VSX,INDX)
!
!---        Hydrodynamic dispersion coefficient w/ flow velocity  ---
!
            IF( ISLC(52).EQ.0 ) THEN
              USX = USX/VMCSX
              VSX = VSX/VMCSX
              WSX = WSX/VMCSX
            ENDIF
            USSQ = USX*USX
            VSSQ = VSX*VSX
            WSSQ = WSX*WSX
            ZVS = SQRT(USSQ+VSSQ+WSSQ)
            INDX = 17
            DISPLS = DISPL(IZ(NS))*SMDEF(IZ(NS),NSL)
            DISPLP = DISPL(IZ(N))*SMDEF(IZ(N),NSL)
            DISPTS = DISPT(IZ(NS))*SMDEF(IZ(NS),NSL)
            DISPTP = DISPT(IZ(N))*SMDEF(IZ(N),NSL)
            DPLS = DIFMN(DISPLS,DISPLP,DYGF(NS),DYGF(N),VSX,INDX)
            DPTS = DIFMN(DISPTS,DISPTP,DYGF(NS),DYGF(N),VSX,INDX)
            DPS = (DPLS*VSSQ + DPTS*(USSQ+WSSQ))/(ZVS+SMALL)
          ELSE
            DPS = 0.D+0
          ENDIF
          DYSR = RP(I)*DYGF(NS)
          FCLS = YL(NS,NSL)/(VMCS+SMALL)
          IF( IEDL(NSL).EQ.2 ) THEN
            DLS = SDCL(1,IZ(NS),NSL)*SDCL(2,IZ(NS),NSL)*
     &        EXP(VMCS*SDCL(3,IZ(NS),NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DLS = TORL(2,NS)*VMCS*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DLS = SDCL(1,IZ(NS),NSL)*SDCL(2,IZ(NS),NSL)*
     &        VMCS**SDCL(3,IZ(NS),NSL)
          ELSE
            DLS = TORL(2,NS)*VMCS*SDFLS
          ENDIF
          INDX = 16
          DLY = DIFMN(DLS,DLP,DYSR,DYR,VL(1,NPY),INDX)
          DLY = AFY(NPY)*(DLY+DPS)/(DYGP(NPY)*RP(I))
          FLS = AFY(NPY)*VL(1,NPY)
          IF( MOD(ISLC(23),10).EQ.1 ) FLS = 0.D+0
!
!---      Patankar solute transport  ---
!
          ALS = MAX(FLS,ZERO)
     &      + DLY*MAX((ONE-(TENTH*ABS(FLS)/(DLY+SMALL)))**5,ZERO)
          AP = (ALS-FLS)*FCLP
          AS = ALS*FCLS
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
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AS

          ENDIF
        ENDDO
!
!---    West face diffusion and advection terms  ---
!
        DO NC = 1,4
          NW = ICM(NC,3,N)
          IF( NW.EQ.0 ) EXIT
!
!---      Multiple west connections  ---
!
          IF( INBS(3,N).LT.0 ) THEN
            NPX = NSSX(NW)
!
!---      Single west connections  ---
!
          ELSE
            NPX = NSX(N)
          ENDIF
          TCOR = (T(2,NW)+TABS)/TSPRF
          SDFLW = SMDL(NSL)*TCOR*(VISRL/VISL(2,NW))
          VMCW = SL(2,NW)*PORD(2,NW)
          IF( IDISP.EQ.1 ) THEN
!
!---        Hydrodynamic dispersion coefficient  ---
!
            UWX = UL(1,NSX(N))
            VWX = 5.D-1*(VL(1,NSY(N))+VL(1,NSY(N)+IFLD))
            WWX = 5.D-1*(WL(1,NSZ(N))+WL(1,NSZ(N)+IJFLD))
            INDX = 17
            VMCWX = DIFMN(VMCW,VMCP,DXGF(NW),DXGF(N),UWX,INDX)
!
!---        Hydrodynamic dispersion coefficient w/ flow velocity  ---
!
            IF( ISLC(52).EQ.0 ) THEN
              UWX = UWX/VMCWX
              VWX = VWX/VMCWX
              WWX = WWX/VMCWX
            ENDIF
            UWSQ = UWX*UWX
            VWSQ = VWX*VWX
            WWSQ = WWX*WWX
            ZVW = SQRT(UWSQ+VWSQ+WWSQ)
            INDX = 17
            DISPLW = DISPL(IZ(NW))*SMDEF(IZ(NW),NSL)
            DISPLP = DISPL(IZ(N))*SMDEF(IZ(N),NSL)
            DISPTW = DISPT(IZ(NW))*SMDEF(IZ(NW),NSL)
            DISPTP = DISPT(IZ(N))*SMDEF(IZ(N),NSL)
            DPLW = DIFMN(DISPLW,DISPLP,DXGF(NW),DXGF(N),UWX,INDX)
            DPTW = DIFMN(DISPTW,DISPTP,DXGF(NW),DXGF(N),UWX,INDX)
            DPW = (DPLW*UWSQ + DPTW*(VWSQ+WWSQ))/(ZVW+SMALL)
          ELSE
            DPW = 0.D+0
          ENDIF
          FCLW = YL(NW,NSL)/(VMCW+SMALL)
          IF( IEDL(NSL).EQ.2 ) THEN
            DLW = SDCL(1,IZ(NW),NSL)*SDCL(2,IZ(NW),NSL)*
     &        EXP(VMCW*SDCL(3,IZ(NW),NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DLW = TORL(2,NW)*VMCW*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DLW = SDCL(1,IZ(NW),NSL)*SDCL(2,IZ(NW),NSL)*
     &        VMCW**SDCL(3,IZ(NW),NSL)
          ELSE
            DLW = TORL(2,NW)*VMCW*SDFLW
          ENDIF
          INDX = 16
          DLX = DIFMN(DLW,DLP,DXGF(NW),DXGF(N),UL(1,NPX),INDX)
          DLX = AFX(NPX)*(DLX+DPW)/DXGP(NPX)
          FLW = AFX(NPX)*UL(1,NPX)
          IF( MOD(ISLC(23),10).EQ.1 ) FLW = 0.D+0
!
!---      Patankar solute transport  ---
!
          ALW = MAX(FLW,ZERO)
     &      + DLX*MAX((ONE-(TENTH*ABS(FLW)/(DLX+SMALL)))**5,ZERO)
          AP = (ALW-FLW)*FCLP
          AW = ALW*FCLW
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
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AW

          ENDIF
        ENDDO
!
!---    East face diffusion and advection terms  ---
!
        DO NC = 1,4
          NE = ICM(NC,4,N)
          IF( NE.EQ.0 ) EXIT
!
!---      Multiple west connections for east node  ---
!
          IF( INBS(3,NE).LT.0 ) THEN
            NQX = NSSX(N)
!
!---      Single west connection for east node  ---
!
          ELSE
            NQX = NSX(NE)
          ENDIF
          TCOR = (T(2,NE)+TABS)/TSPRF
          SDFLE = SMDL(NSL)*TCOR*(VISRL/VISL(2,NE))
          VMCE = SL(2,NE)*PORD(2,NE)
          IF( IDISP.EQ.1 ) THEN
!
!---        Hydrodynamic dispersion coefficient  ---
!
            UEX = UL(1,NSX(N)+1)
            VEX = 5.D-1*(VL(1,NSY(N))+VL(1,NSY(N)+IFLD))
            WEX = 5.D-1*(WL(1,NSZ(N))+WL(1,NSZ(N)+IJFLD))
            INDX = 17
            VMCEX = DIFMN(VMCP,VMCE,DXGF(N),DXGF(NE),UEX,INDX)
!
!---        Hydrodynamic dispersion coefficient w/ flow velocity  ---
!
            IF( ISLC(52).EQ.0 ) THEN
              UEX = UEX/VMCEX
              VEX = VEX/VMCEX
              WEX = WEX/VMCEX
            ENDIF
            UESQ = UEX*UEX
            VESQ = VEX*VEX
            WESQ = WEX*WEX
            ZVE = SQRT(UESQ+VESQ+WESQ)
            INDX = 17
            DISPLE = DISPL(IZ(NE))*SMDEF(IZ(NE),NSL)
            DISPLP = DISPL(IZ(N))*SMDEF(IZ(N),NSL)
            DISPTE = DISPT(IZ(NE))*SMDEF(IZ(NE),NSL)
            DISPTP = DISPT(IZ(N))*SMDEF(IZ(N),NSL)
            DPLE = DIFMN(DISPLP,DISPLE,DXGF(N),DXGF(NE),UEX,INDX)
            DPTE = DIFMN(DISPTP,DISPTE,DXGF(N),DXGF(NE),UEX,INDX)
            DPE = (DPLE*UESQ + DPTE*(VESQ+WESQ))/(ZVE+SMALL)
         ELSE
            DPE = 0.D+0
          ENDIF
          FCLE = YL(NE,NSL)/(VMCE+SMALL)
          IF( IEDL(NSL).EQ.2 ) THEN
            DLE = SDCL(1,IZ(NE),NSL)*SDCL(2,IZ(NE),NSL)*
     &        EXP(VMCE*SDCL(3,IZ(NE),NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DLE = TORL(2,NE)*VMCE*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DLE = SDCL(1,IZ(NE),NSL)*SDCL(2,IZ(NE),NSL)*
     &        VMCE**SDCL(3,IZ(NE),NSL)
          ELSE
            DLE = TORL(2,NE)*VMCE*SDFLE
          ENDIF
          INDX = 16
          DLX = DIFMN(DLP,DLE,DXGF(N),DXGF(NE),UL(1,NQX),INDX)
          DLX = AFX(NQX)*(DLX+DPE)/DXGP(NQX)
          FLE = AFX(NQX)*UL(1,NQX)
          IF( MOD(ISLC(23),10).EQ.1 ) FLE = 0.D+0
!
!---      Patankar solute transport  ---
!
          ALE = MAX(-FLE,ZERO)
     &      + DLX*MAX((ONE-(TENTH*ABS(FLE)/(DLX+SMALL)))**5,ZERO)
          AP = (ALE+FLE)*FCLP
          AE = ALE*FCLE
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
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AE

          ENDIF
        ENDDO
!
!---    North face diffusion and advection terms  ---
!
        DO NC = 1,4
          NN = ICM(NC,5,N)
          IF( NN.EQ.0 ) EXIT
!
!---      Multiple south connections for north node  ---
!
          IF( INBS(2,NN).LT.0 ) THEN
            NQY = NSSY(N)
!
!---      Single south connection for north node  ---
!
          ELSE
            NQY = NSY(NN)
          ENDIF
          TCOR = (T(2,NN)+TABS)/TSPRF
          SDFLN = SMDL(NSL)*TCOR*(VISRL/VISL(2,NN))
          VMCN = SL(2,NN)*PORD(2,NN)
          IF( IDISP.EQ.1 ) THEN
!
!---        Hydrodynamic dispersion coefficient  ---
!
            UNX = 5.D-1*(UL(1,NSX(N))+UL(1,NSX(N)+1))
            WNX = 5.D-1*(WL(1,NSZ(N))+WL(1,NSZ(N)+IJFLD))
            VNX = VL(1,NSY(N)+IFLD)
            INDX = 17
            VMCNX = DIFMN(VMCP,VMCN,DYNR,DYR,VNX,INDX)
!
!---        Hydrodynamic dispersion coefficient w/ flow velocity  ---
!
            IF( ISLC(52).EQ.0 ) THEN
              UNX = UNX/VMCNX
              VNX = VNX/VMCNX
              WNX = WNX/VMCNX
            ENDIF
            UNSQ = UNX*UNX
            VNSQ = VNX*VNX
            WNSQ = WNX*WNX
            ZVN = SQRT(UNSQ+VNSQ+WNSQ)
            INDX = 17
            DISPLN = DISPL(IZ(NN))*SMDEF(IZ(NN),NSL)
            DISPLP = DISPL(IZ(N))*SMDEF(IZ(N),NSL)
            DISPTN = DISPT(IZ(NN))*SMDEF(IZ(NN),NSL)
            DISPTP = DISPT(IZ(N))*SMDEF(IZ(N),NSL)
            DPLN = DIFMN(DISPLP,DISPLN,DYGF(N),DYGF(NN),VNX,INDX)
            DPTN = DIFMN(DISPTP,DISPTN,DYGF(N),DYGF(NN),VNX,INDX)
            DPN = (DPLN*VNSQ + DPTN*(UNSQ+WNSQ))/(ZVN+SMALL)
          ELSE
            DPN = 0.D+0
          ENDIF
          DYNR = RP(I)*DYGF(NN)
          FCLN = YL(NN,NSL)/(VMCN+SMALL)
          IF( IEDL(NSL).EQ.2 ) THEN
            DLN = SDCL(1,IZ(NN),NSL)*SDCL(2,IZ(NN),NSL)*
     &        EXP(VMCN*SDCL(3,IZ(NN),NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DLN = TORL(2,NN)*VMCN*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DLN = SDCL(1,IZ(NN),NSL)*SDCL(2,IZ(NN),NSL)*
     &        VMCN**SDCL(3,IZ(NN),NSL)
          ELSE
            DLN = TORL(2,NN)*VMCN*SDFLN
          ENDIF
          INDX = 16
          DLY = DIFMN(DLP,DLN,DYNR,DYR,VL(1,NQY),INDX)
          DLY = AFY(NQY)*(DLY+DPN)/(DYGP(NQY)*RP(I))
          FLN = AFY(NQY)*VL(1,NQY)
          IF( MOD(ISLC(23),10).EQ.1 ) FLS = 0.D+0
          VMCNX = DIFMN(VMCP,VMCN,DYNR,DYR,VL(1,NQY),INDX)
          CRLN = ABS(VL(1,NQY))*DT/(DYGP(NQY)*VMCNX+SMALL)/RP(I)
!
!---      Patankar solute transport  ---
!
          ALN = MAX(-FLN,ZERO)
     &      + DLY*MAX((ONE-(TENTH*ABS(FLN)/(DLY+SMALL)))**5,ZERO)
          AP = (ALN+FLN)*FCLP
          AN = ALN*FCLN
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
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AN

          ENDIF
        ENDDO
!
!---    Top face diffusion and advection terms  ---
!
        DO NC = 1,4
          NT = ICM(NC,6,N)
          IF( NT.EQ.0 ) EXIT
!
!---      Multiple bottom connections for top node  ---
!
          IF( INBS(1,NT).LT.0 ) THEN
            NQZ = NSSZ(N)
!
!---      Single bottom connection for top node  ---
!
          ELSE
            NQZ = NSZ(NT)
          ENDIF
          TCOR = (T(2,NT)+TABS)/TSPRF
          SDFLT = SMDL(NSL)*TCOR*(VISRL/VISL(2,NT))
          VMCT = SL(2,NT)*PORD(2,NT)
          IF( IDISP.EQ.1 ) THEN
!
!---        Hydrodynamic dispersion coefficient  ---
!
            UTX = 5.D-1*(UL(1,NSX(N))+UL(1,NSX(N)+1))
            VTX = 5.D-1*(VL(1,NSY(N))+VL(1,NSY(N)+IFLD))
            WTX = WL(1,NSZ(N)+IJFLD)
            INDX = 17
            VMCTX = DIFMN(VMCP,VMCT,DZGF(N),DZGF(NT),WTX,INDX)
!
!---        Hydrodynamic dispersion coefficient w/ flow velocity  ---
!
            IF( ISLC(52).EQ.0 ) THEN
              UTX = UTX/VMCTX
              VTX = VTX/VMCTX
              WTX = WTX/VMCTX
            ENDIF
            UTSQ = UTX*UTX
            VTSQ = VTX*VTX
            WTSQ = WTX*WTX
            ZVT = SQRT(UTSQ+VTSQ+WTSQ)
            INDX = 17
            DISPLT = DISPL(IZ(NT))*SMDEF(IZ(NT),NSL)
            DISPLP = DISPL(IZ(N))*SMDEF(IZ(N),NSL)
            DISPTT = DISPT(IZ(NT))*SMDEF(IZ(NT),NSL)
            DISPTP = DISPT(IZ(N))*SMDEF(IZ(N),NSL)
            DPLT = DIFMN(DISPLP,DISPLT,DZGF(N),DZGF(NT),WTX,INDX)
            DPTT = DIFMN(DISPTP,DISPTT,DZGF(N),DZGF(NT),WTX,INDX)
            DPT = (DPLT*WTSQ + DPTT*(VTSQ+UTSQ))/(ZVT+SMALL)
          ELSE
            DPT = 0.D+0
          ENDIF
          FCLT = YL(NT,NSL)/(VMCT+SMALL)
          IF( IEDL(NSL).EQ.2 ) THEN
            DLT = SDCL(1,IZ(NT),NSL)*SDCL(2,IZ(NT),NSL)*
     &        EXP(VMCT*SDCL(3,IZ(NT),NSL))
          ELSEIF( IEDL(NSL).EQ.3 ) THEN
            DLT = TORL(2,NT)*VMCT*SMDL(NSL)
          ELSEIF( IEDL(NSL).EQ.4 ) THEN
            DLT = SDCL(1,IZ(NT),NSL)*SDCL(2,IZ(NT),NSL)*
     &        VMCT**SDCL(3,IZ(NT),NSL)
          ELSE
            DLT = TORL(2,NT)*VMCT*SDFLT
          ENDIF
          INDX = 16
          DLZ = DIFMN(DLP,DLT,DZGF(N),DZGF(NT),WL(1,NQZ),INDX)
          DLZ = AFZ(NQZ)*(DLZ+DPT)/DZGP(NQZ)
          FLT = AFZ(NQZ)*WL(1,NQZ)
          IF( MOD(ISLC(23),10).EQ.1 ) FLT = 0.D+0
          VMCTX = DIFMN(VMCP,VMCT,DZGF(N),DZGF(NT),WL(1,NQZ),INDX)
          CRLT = ABS(WL(1,NQZ))*DT/(DZGP(NQZ)*VMCTX+SMALL)
!
!---      Patankar solute transport  ---
!
          ALT = MAX(-FLT,ZERO)
     &      + DLZ*MAX((ONE-(TENTH*ABS(FLT)/(DLZ+SMALL)))**5,ZERO)
          AP = (ALT+FLT)*FCLP
          AT = ALT*FCLT
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
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
            DLU(MCD) = DLU(MCD) + AP
            MROW = KLUC(MCP,MA)
            MA = MA + 1
            DLU(MROW) = DLU(MROW) - AT

          ENDIF
        ENDDO
!!
!!---    Schroth half-life modification model  ---
!!
!        IF( PCLN(3,NSL).LT.0.D+0 ) THEN
!          IF( SL(2,N).LT.PCLN(1,NSL) ) THEN
!            HLFX = 1.D+12
!          ELSE
!            HLFX = HLF(NSL)*(PCLN(2,N) + SL(2,N) - 2.D+0*PCLN(1,N))/
!     &        (SL(2,N) - PCLN(1,N) + SMALL)
!          ENDIF
!        ELSE
!          HLFX = HLF(NSL)
!        ENDIF
!
!---    Solution vector and parent 1 decay  ---
!
        BLU(MCP) = BLU(MCP) + CO(N,NSL)*SC
!     &    - 6.931D-1*CO(N,NSL)*VOL(N)/HLFX
!!        BLU(MCP) = BLU(MCP) + CO(N,NSL)*(2.D+0**(-DT/HLF(NSL)))*SC
!!
!!---    Daughter 1 chain decay  ---
!!
!        IF( NSL.LE.NSOLU ) THEN
!          DO 170 NPSL = 1,NSOLU
!            IF( NPSL.EQ.NSL ) GOTO 170
!            BLU(MCP) = BLU(MCP) +
!     &        CHDF(NPSL,NSL)*6.931D-1*CO(N,NPSL)*VOL(N)/HLF(NPSL)
!!            BLU(MCP) = BLU(MCP) +
!!     &        CHDF(NPSL,NSL)*CO(N,NPSL)*(2.D+0**(-DT/HLF(NPSL)))*SC
!  170     CONTINUE
!        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SJCBL group  ---
!
      RETURN
      END



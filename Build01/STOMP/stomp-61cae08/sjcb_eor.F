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
#ifdef petsc
!
!----------------------PETSc Modules-----------------------------------!
!
      USE PETSC_STOMP
#endif
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE JACOB
      USE HYST
      USE GRID
      USE GLB_PAR
      USE FLUXP
      USE FDVP
      USE FDVG
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#ifdef petsc
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER :: IROW,ICOL
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SJCBG'
!
!---  Fill matrix elements  ---
!
      DO 900 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 900
        I = ID(N)
        IZN = IZ(N)
        DYR = RP(I)*DYGF(N)
        MCP = IXP(N)
        MA = 1
!
!---    Banded solver  ---
!
        IF( ILES.EQ.1 ) THEN
          MCD = MCP
          MRD = MDT
!
!---    SPLib solver  ---
!
        ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
          MCD = KLUC(MCP,MA)
          MA = MA + 1
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
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AP,
     &                        ADD_VALUES,IERR)
            ICOL = petsc_offset(NB)
            AB = -AB
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AB,
     &                        ADD_VALUES,IERR)
#endif
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
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AP,
     &                        ADD_VALUES,IERR)
            ICOL = petsc_offset(NS)
            AS = -AS
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AS,
     &                        ADD_VALUES,IERR)
#endif
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
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AP,
     &                        ADD_VALUES,IERR)
            ICOL = petsc_offset(NW)
            AW = -AW
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AW,
     &                        ADD_VALUES,IERR)
#endif
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
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AP,
     &                        ADD_VALUES,IERR)
            ICOL = petsc_offset(NE)
            AE = -AE
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AE,
     &                        ADD_VALUES,IERR)
#endif
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
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AP,
     &                        ADD_VALUES,IERR)
            ICOL = petsc_offset(NN)
            AN = -AN
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AN,
     &                        ADD_VALUES,IERR)
#endif
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
          FCGT = YG(NT,NSL)/(VMCT+SMALL)
          DGT = TORG(2,NT)*(SG(2,NT)-SGT(2,NT))*PORD(2,NT)*SDFGT
          INDX = 16
          DGZ = DIFMN(DGP,DGT,DZGF(N),DZGF(NT),WG(1,NQZ),INDX)
          DGZ = AFZ(NQZ)*(DGZ+DPGT)/DZGP(NQZ)
          FGT = AFZ(NQZ)*WG(1,NQZ)
          IF( MOD(ISLC(23),100)/10.EQ.1 ) FGT = 0.D+0
          VMCTX = DIFMN(VMCP,VMCT,DZGF(N),DZGF(NT),WG(1,NQZ),INDX)
          CRGT = ABS(WG(1,NQZ))*DT/(DZGP(NQZ)*VMCTX+SMALL)
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
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AP,
     &                        ADD_VALUES,IERR)
            ICOL = petsc_offset(NT)
            AT = -AT
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AT,
     &                        ADD_VALUES,IERR)
#endif
          ENDIF
        ENDDO
  900 CONTINUE
!
!---  End of SJCBG group  ---
!
      ISUB_LOG = ISUB_LOG-1
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
#ifdef petsc
!
!----------------------PETSc Modules-----------------------------------!
!
      USE PETSC_STOMP
#endif
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
      USE FDVP
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#ifdef petsc
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER :: IROW,ICOL
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SJCBL'
!
!---  Fill matrix elements  ---
!
      DO 900 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 900
        IZN = IZ(N)
        I = ID(N)
        DYR = RP(I)*DYGF(N)
!
!---    Storage terms  ---
!
        SC = VOL(N)*DTI
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
#ifdef petsc
!
!---    PETSc solver  ---
!
        ELSEIF( ILES.EQ.5 ) THEN
          IROW = petsc_offset(N)
          CALL MatSetValues(T_MAT,1,IROW,1,IROW,SC,
     &                      ADD_VALUES,IERR)
#endif
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
!---    Hydrodynamic dispersion coefficients at cell faces  ---
!
        CALL SHDPL( N,DPLB,DPLS,DPLW,DPLE,DPLN,DPLT,NSL )
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
          DLZ = AFZ(NPZ)*(DLZ+DPLB)/DZGP(NPZ)
          FLB = AFZ(NPZ)*WL(1,NPZ)
          IF( MOD(ISLC(23),10).EQ.1 ) FLB = 0.D+0
          VMCBX = DIFMN(VMCB,VMCP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX)
          CRLB = ABS(WL(1,NPZ))*DT/(DZGP(NPZ)*VMCBX+SMALL)
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
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AP,
     &                        ADD_VALUES,IERR)
            ICOL = petsc_offset(NB)
            AB = -AB
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AB,
     &                        ADD_VALUES,IERR)
#endif
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
          DYSR = RP(I)*DYGF(NS)
          TCOR = (T(2,NS)+TABS)/TSPRF
          SDFLS = SMDL(NSL)*TCOR*(VISRL/VISL(2,NS))
          VMCS = SL(2,NS)*PORD(2,NS)
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
          DLY = AFY(NPY)*(DLY+DPLS)/(DYGP(NPY)*RP(I))
          FLS = AFY(NPY)*VL(1,NPY)
          IF( MOD(ISLC(23),10).EQ.1 ) FLS = 0.D+0
          VMCSX = DIFMN(VMCS,VMCP,DYSR,DYR,VL(1,NPY),INDX)
          CRLS = ABS(VL(1,NPY))*DT/(DYGP(NPY)*VMCSX+SMALL)/RP(I)
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
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AP,
     &                        ADD_VALUES,IERR)
            ICOL = petsc_offset(NS)
            AS = -AS
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AS,
     &                        ADD_VALUES,IERR)
#endif
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
          DLX = AFX(NPX)*(DLX+DPLW)/DXGP(NPX)
          FLW = AFX(NPX)*UL(1,NPX)
          IF( MOD(ISLC(23),10).EQ.1 ) FLW = 0.D+0
          VMCWX = DIFMN(VMCW,VMCP,DXGF(NW),DXGF(N),UL(1,NPX),INDX)
          CRLW = ABS(UL(1,NPX))*DT/(DXGP(NPX)*VMCWX+SMALL)
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
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AP,
     &                        ADD_VALUES,IERR)
            ICOL = petsc_offset(NW)
            AW = -AW
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AW,
     &                        ADD_VALUES,IERR)
#endif
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
          DLX = AFX(NQX)*(DLX+DPLE)/DXGP(NQX)
          FLE = AFX(NQX)*UL(1,NQX)
          IF( MOD(ISLC(23),10).EQ.1 ) FLE = 0.D+0
          VMCEX = DIFMN(VMCP,VMCE,DXGF(N),DXGF(NE),UL(1,NQX),INDX)
          CRLE = ABS(UL(1,NQX))*DT/(DXGP(NQX)*VMCEX+SMALL)
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
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AP,
     &                        ADD_VALUES,IERR)
            ICOL = petsc_offset(NE)
            AE = -AE
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AE,
     &                        ADD_VALUES,IERR)
#endif
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
          DYNR = RP(I)*DYGF(NN)
          TCOR = (T(2,NN)+TABS)/TSPRF
          SDFLN = SMDL(NSL)*TCOR*(VISRL/VISL(2,NN))
          VMCN = SL(2,NN)*PORD(2,NN)
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
          DLY = AFY(NQY)*(DLY+DPLN)/(DYGP(NQY)*RP(I))
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
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AP,
     &                        ADD_VALUES,IERR)
            ICOL = petsc_offset(NN)
            AN = -AN
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AN,
     &                        ADD_VALUES,IERR)
#endif
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
          DLZ = AFZ(NQZ)*(DLZ+DPLT)/DZGP(NQZ)
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
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AP,
     &                        ADD_VALUES,IERR)
            ICOL = petsc_offset(NT)
            AT = -AT
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AT,
     &                        ADD_VALUES,IERR)
#endif
          ENDIF
        ENDDO
!
!---    Solution vector  ---
!
        BLU(MCP) = BLU(MCP) + CO(N,NSL)*SC
  900 CONTINUE
!
!---  End of SJCBL group  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SJCBN( NSL )
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
!---------------------!opyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Loads the matrix elements and solution vector for the
!     NAPL convective-dispersive mass transport equation.
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
#ifdef petsc
!
!----------------------PETSc Modules-----------------------------------!
!
      USE PETSC_STOMP
#endif
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE NAPL
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXN
      USE FDVP
      USE FDVN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#ifdef petsc
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER :: IROW,ICOL
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SJCBN'
!
!---  Fill matrix elements  ---
!
      DO 900 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 900
        I = ID(N)
        IZN = IZ(N)
        DYR = RP(I)*DYGF(N)
        MCP = IXP(N)
        MA = 1
!
!---    Banded solver  ---
!
        IF( ILES.EQ.1 ) THEN
          MCD = MCP
          MRD = MDT
!
!---    SPLib solver  ---
!
        ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
          MCD = KLUC(MCP,MA)
          MA = MA + 1
        ENDIF
!
!---  Molecular diffusion coefficients at the nodes  ---
!
        TCOR = (T(2,N)+TABS)/TSPRF
        SDFNP = SMDN(NSL)*TCOR*(VISRN/VISN(2,N))
        VMCP = SN(2,N)*PORD(2,N)
        FCNP = YN(N,NSL)/(VMCP+SMALL)
        DNP = TORN(2,N)*VMCP*SDFNP
!
!---  Hydrodynamic dispersion coefficients at cell faces  ---
!
        CALL SHDPN( N,DPNB,DPNS,DPNW,DPNE,DPNN,DPNT )
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
          SDFNB = SMDN(NSL)*TCOR*(VISRN/VISN(2,NB))
          VMCB = SN(2,NB)*PORD(2,NB)
          FCNB = YN(NB,NSL)/(VMCB+SMALL)
          DNB = TORN(2,NB)*VMCB*SDFNB
          INDX = 16
          DNZ = DIFMN(DNB,DNP,DZGF(NB),DZGF(N),WN(1,NPZ),INDX)
          DNZ = AFZ(NPZ)/DZGP(NPZ)*(DNZ+DPNB)
          FNB = AFZ(NPZ)*WN(1,NPZ)
          IF( MOD(ISLC(23),1000)/100.EQ.1 ) FNB = 0.D+0
          VMCBX = DIFMN(VMCB,VMCP,DZGF(NB),DZGF(N),WN(1,NPZ),INDX)
          CRNB = ABS(WN(1,NPZ))*DT/DZGP(NPZ)/(VMCBX+SMALL)
!
!---      Patankar solute transport  ---
!
          ANB = MAX(FNB,ZERO)
     &      + DNZ*MAX((ONE-(TENTH*ABS(FNB)/(DNZ+SMALL)))**5,ZERO)
          AP = (ANB-FNB)*FCNP
          AB = ANB*FCNB
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
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AP,
     &                        ADD_VALUES,IERR)
            ICOL = petsc_offset(NB)
            AB = -AB
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AB,
     &                        ADD_VALUES,IERR)
#endif
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
          DYSR = RP(I)*DYGF(NS)
          TCOR = (T(2,NS)+TABS)/TSPRF
          SDFNS = SMDN(NSL)*TCOR*(VISRN/VISN(2,NS))
          VMCS = SN(2,NS)*PORD(2,NS)
          FCNS = YN(NS,NSL)/(VMCS+SMALL)
          DNS = TORN(2,NS)*VMCS*SDFNS
          INDX = 16
          DNY = DIFMN(DNS,DNP,DYSR,DYR,VN(1,NPY),INDX)
          DNY = AFY(NPY)/DYGP(NPY)*(DNY+DPNS)/RP(I)
          FNS = AFY(NPY)*VN(1,NPY)
          IF( MOD(ISLC(23),1000)/100.EQ.1 ) FNS = 0.D+0
          VMCSX = DIFMN(VMCS,VMCP,DYSR,DYR,VN(1,NPY),INDX)
          CRNS = ABS(VN(1,NPY))*DT/DYGP(NPY)/(VMCSX+SMALL)/RP(I)
!
!---      Patankar solute transport  ---
!
          ANS = MAX(FNS,ZERO)
     &      + DNY*MAX((ONE-(TENTH*ABS(FNS)/(DNY+SMALL)))**5,ZERO)
          AP = (ANS-FNS)*FCNP
          AS = ANS*FCNS
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
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AP,
     &                        ADD_VALUES,IERR)
            ICOL = petsc_offset(NS)
            AS = -AS
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AS,
     &                        ADD_VALUES,IERR)
#endif
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
          SDFNW = SMDN(NSL)*TCOR*(VISRN/VISN(2,NW))
          VMCW = SN(2,NW)*PORD(2,NW)
          FCNW = YN(NW,NSL)/(VMCW+SMALL)
          DNW = TORN(2,NW)*VMCW*SDFNW
          INDX = 16
          DNX = DIFMN(DNW,DNP,DXGF(NW),DXGF(N),UN(1,NPX),INDX)
          DNX = AFX(NPX)*(DNX+DPNW)/DXGP(NPX)
          FNW = AFX(NPX)*UN(1,NPX)
          IF( MOD(ISLC(23),1000)/100.EQ.1 ) FNW = 0.D+0
          VMCWX = DIFMN(VMCW,VMCP,DXGF(NW),DXGF(N),UN(1,NPX),INDX)
          CRNW = ABS(UN(1,NPX))*DT/DXGP(NPX)/(VMCWX+SMALL)
!
!---      Patankar solute transport  ---
!
          ANW = MAX(FNW,ZERO)
     &      + DNX*MAX((ONE-(TENTH*ABS(FNW)/(DNX+SMALL)))**5,ZERO)
          AP = (ANW-FNW)*FCNP
          AW = ANW*FCNW
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
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AP,
     &                        ADD_VALUES,IERR)
            ICOL = petsc_offset(NW)
            AW = -AW
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AW,
     &                        ADD_VALUES,IERR)
#endif
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
          SDFNE = SMDN(NSL)*TCOR*(VISRN/VISN(2,NE))
          VMCE = SN(2,NE)*PORD(2,NE)
          FCNE = YN(NE,NSL)/(VMCE+SMALL)
          DNE = TORN(2,NE)*VMCE*SDFNE
          INDX = 16
          DNX = DIFMN(DNP,DNE,DXGF(N),DXGF(NE),UN(1,NQX),INDX)
          DNX = AFX(NQX)/DXGP(NQX)*(DNX+DPNE)
          FNE = AFX(NQX)*UN(1,NQX)
          IF( MOD(ISLC(23),1000)/100.EQ.1 ) FNE = 0.D+0
          VMCEX = DIFMN(VMCP,VMCE,DXGF(N),DXGF(NE),UN(1,NQX),INDX)
          CRNE = ABS(UN(1,NQX))*DT/DXGP(NQX)/(VMCEX+SMALL)
!
!---      Patankar solute transport  ---
!
          ANE = MAX(-FNE,ZERO)
     &      + DNX*MAX((ONE-(TENTH*ABS(FNE)/(DNX+SMALL)))**5,ZERO)
          AP = (ANE+FNE)*FCNP
          AE = ANE*FCNE
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
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AP,
     &                        ADD_VALUES,IERR)
            ICOL = petsc_offset(NE)
            AE = -AE
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AE,
     &                        ADD_VALUES,IERR)
#endif
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
          DYNR = RP(I)*DYGF(NN)
          TCOR = (T(2,NN)+TABS)/TSPRF
          SDFNN = SMDN(NSL)*TCOR*(VISRN/VISN(2,NN))
          VMCN = SN(2,NN)*PORD(2,NN)
          FCNN = YN(NN,NSL)/(VMCN+SMALL)
          DNN = TORN(2,NN)*VMCN*SDFNN
          INDX = 16
          DNY = DIFMN(DNP,DNN,DYNR,DYR,VN(1,NQY),INDX)
          DNY = AFY(NQY)/DYGP(NQY)*(DNY+DPNN)/RP(I)
          FNN = AFY(NQY)*VN(1,NQY)
          IF( MOD(ISLC(23),1000)/100.EQ.1 ) FNN = 0.D+0
          VMCNX = DIFMN(VMCP,VMCN,DYNR,DYR,VN(1,NQY),INDX)
          CRNN = ABS(VN(1,NQY))*DT/DYGP(NQY)/(VMCNX+SMALL)/RP(I)
!
!---      Patankar solute transport  ---
!
          ANN = MAX(-FNN,ZERO)
     &      + DNY*MAX((ONE-(TENTH*ABS(FNN)/(DNY+SMALL)))**5,ZERO)
          AP = (ANN+FNN)*FCNP
          AN = ANN*FCNN
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
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AP,
     &                        ADD_VALUES,IERR)
            ICOL = petsc_offset(NN)
            AN = -AN
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AN,
     &                        ADD_VALUES,IERR)
#endif
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
          SDFNT = SMDN(NSL)*TCOR*(VISRN/VISN(2,NT))
          VMCT = SN(2,NT)*PORD(2,NT)
          FCNT = YN(NT,NSL)/(VMCT+SMALL)
          DNT = TORN(2,NT)*VMCT*SDFNT
          INDX = 16
          DNZ = DIFMN(DNP,DNT,DZGF(N),DZGF(NT),WN(1,NQZ),INDX)
          DNZ = AFZ(NQZ)/DZGP(NQZ)*(DNZ+DPNT)
          FNT = AFZ(NQZ)*WN(1,NQZ)
          IF( MOD(ISLC(23),1000)/100.EQ.1 ) FNT = 0.D+0
          VMCTX = DIFMN(VMCP,VMCT,DZGF(N),DZGF(NT),WN(1,NQZ),INDX)
          CRNT = ABS(WN(1,NQZ))*DT/DZGP(NQZ)/(VMCTX+SMALL)
!
!---      Patankar solute transport  ---
!
          ANT = MAX(-FNT,ZERO)
     &      + DNZ*MAX((ONE-(TENTH*ABS(FNT)/(DNZ+SMALL)))**5,ZERO)
          AP = (ANT+FNT)*FCNP
          AT = ANT*FCNT
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
#ifdef petsc
!
!---      PETSc solver  ---
!
          ELSEIF( ILES.EQ.5 ) THEN
            CALL MatSetValues(T_MAT,1,IROW,1,IROW,AP,
     &                        ADD_VALUES,IERR)
            ICOL = petsc_offset(NT)
            AT = -AT
            CALL MatSetValues(T_MAT,1,IROW,1,ICOL,AT,
     &                        ADD_VALUES,IERR)
#endif
          ENDIF
        ENDDO
  900 CONTINUE
!
!---  End of SJCBN group  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END


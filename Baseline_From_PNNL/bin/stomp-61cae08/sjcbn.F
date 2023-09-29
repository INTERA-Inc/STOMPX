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
!     Last Modified by MD White on September 5, 1996.
!     sjcbn.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      INTEGER :: IROW, ICOL
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SJCBN'
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
        NWW = NW-1
        NEE = NE+1
        NSS = NS-IFLD
        NNN = NN+IFLD
        NBB = NB-IJFLD
        NTT = NT+IJFLD
        IZN = IZ(N)
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NSX(N) + 1
        NQY = NSY(N) + IFLD
        NQZ = NSZ(N) + IJFLD
        DYR = RP(I)*DYGF(N)
!
!---  Storage terms  ---
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
          IF( IEQW.EQ.0 .AND. IEQA.EQ.0 )
     &      ALU(MRD,MCD) = ALU(MRD,MCD) + SC
!
!---    SPLib solver  ---
!
        ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
          MCD = KLUC(MCP,MA)
          MA = MA + 1
          IF( IEQW.EQ.0 .AND. IEQA.EQ.0 )
     &      DLU(MCD) = DLU(MCD) + SC
#ifdef petsc
!
!---    PETSc solver  ---
!
        ELSEIF( ILES.EQ.5 ) THEN
          IROW = petsc_offset(N)
          IF( IEQW.EQ.0 .AND. IEQA.EQ.0 ) 
     &      CALL MatSetValues(T_MAT,1,IROW,1,IROW,SC,
     &                        ADD_VALUES,IERR)
#endif
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
!---  Bottom face diffusion and advection terms  ---
!
        IF( K.NE.1 ) THEN
          IF( IXP(NB).EQ.0 .OR. INBS(1,N).GT.0 ) GOTO 110
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
          IF( ISLC(1).GE.1 ) THEN
!
!---  TVD solute transport  ---
!
            IF( FNB.GE.ZERO .AND. K.GT.2 ) THEN
             IF( IXP(NB-IJFLD).GT.0 ) THEN
               FCNBB = YN(NB-IJFLD,NSL)
     &          /(SN(2,NB-IJFLD)*PORD(2,NB-IJFLD)+SMALL)
              R = ((CO(NB,NSL)*FCNB-CO(NB-IJFLD,NSL)*FCNBB)
     &          /(CO(N,NSL)*FCNP-CO(NB,NSL)*FCNB+EPSL))
     &          *((DZGF(N)+DZGF(NB))/(DZGF(NB)+DZGF(NBB)))
              THETA = FLIMIT( R,CRNB,ISLC(1) )
              DZF = DZGF(NB)/(DZGF(N)+DZGF(NB))
              WCZ = CO(NB,NSL)*FNB*(1.D+0-THETA*DZF)*FCNB
     &          + CO(N,NSL)*FNB*THETA*DZF*FCNP
              WC(NPZ,NSL) = WC(NPZ,NSL) + WCZ/AFZ(NPZ)
              CC = CC + WCZ
             ELSE
              WCZ = CO(NB,NSL)*FNB*FCNB
              WC(NPZ,NSL) = WC(NPZ,NSL) + WCZ/AFZ(NPZ)
              CC = CC + WCZ
             ENDIF
            ELSEIF( FNB.GE.ZERO .AND. K.EQ.2 ) THEN
              WCZ = CO(NB,NSL)*FNB*FCNB
              WC(NPZ,NSL) = WC(NPZ,NSL) + WCZ/AFZ(NPZ)
              CC = CC + WCZ
            ELSEIF( FNB.LT.ZERO .AND. K.LT.KFLD ) THEN
             IF( IXP(NT).GT.0 ) THEN
              FCNT = YN(NT,NSL)/(SN(2,NT)*PORD(2,NT)+SMALL)
              R = ((CO(N,NSL)*FCNP-CO(NT,NSL)*FCNT)
     &          /(CO(NB,NSL)*FCNB-CO(N,NSL)*FCNP+EPSL))
     &          *((DZGF(NB)+DZGF(N))/(DZGF(N)+DZGF(NT)))
              THETA = FLIMIT( R,CRNB,ISLC(1) )
              DZF = DZGF(N)/(DZGF(N)+DZGF(NB))
              WCZ = CO(NB,NSL)*FNB*THETA*DZF*FCNB
     &          + CO(N,NSL)*FNB*(1.D+0-THETA*DZF)*FCNP
              WC(NPZ,NSL) = WC(NPZ,NSL) + WCZ/AFZ(NPZ)
              CC = CC + WCZ
             ELSE
              WCZ = CO(N,NSL)*FNB*FCNP
              WC(NPZ,NSL) = WC(NPZ,NSL) + WCZ/AFZ(NPZ)
              CC = CC + WCZ
             ENDIF
            ELSEIF( FNB.LT.ZERO .AND. K.EQ.KFLD ) THEN
              WCZ = CO(N,NSL)*FNB*FCNP
              WC(NPZ,NSL) = WC(NPZ,NSL) + WCZ/AFZ(NPZ)
              CC = CC + WCZ
            ENDIF
            AB = DNZ*FCNB
            AP = DNZ*FCNP
          ELSE
!
!---  Patankar solute transport  ---
!
            ANB = MAX(FNB,ZERO)
     &        + DNZ*MAX((ONE-(TENTH*ABS(FNB)/(DNZ+SMALL)))**5,ZERO)
            AP = (ANB-FNB)*FCNP
            AB = ANB*FCNB
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
        ENDIF
  110   CONTINUE
!
!---  South face diffusion and advection terms  ---
!
        IF( J.NE.1 ) THEN
          IF( IXP(NS).EQ.0 .OR. INBS(2,N).GT.0 ) GOTO 120
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
          IF( ISLC(1).GE.1 ) THEN
!
!---  TVD solute transport  ---
!
            IF( FNS.GE.ZERO .AND. J.GT.2 ) THEN
             IF( IXP(NS-IFLD).GT.0 ) THEN
              FCNSS = YN(NS-IFLD,NSL)
     &          /(SN(2,NS-IFLD)*PORD(2,NS-IFLD)+SMALL)
              R = ((CO(NS,NSL)*FCNS-CO(NS-IFLD,NSL)*FCNSS)
     &          /(CO(N,NSL)*FCNP-CO(NS,NSL)*FCNS+EPSL))
     &          *((DYGF(N)+DYGF(NS))/(DYGF(NS)+DYGF(NSS)))
              THETA = FLIMIT( R,CRNS,ISLC(1) )
              DYF = DYGF(NS)/(DYGF(N)+DYGF(NS))
              VCY = CO(NS,NSL)*FNS*(1.D+0-THETA*DYF)*FCNS
     &          + CO(N,NSL)*FNS*THETA*DYF*FCNP
              VC(NPY,NSL) = VC(NPY,NSL) + VCY/AFY(NPY)
              CC = CC + VCY
             ELSE
              VCY = CO(NS,NSL)*FNS*FCNS
              VC(NPY,NSL) = VC(NPY,NSL) + VCY/AFY(NPY)
              CC = CC + VCY
             ENDIF
            ELSEIF( FNS.GE.ZERO .AND. J.EQ.2 ) THEN
              VCY = CO(NS,NSL)*FNS*FCNS
              VC(NPY,NSL) = VC(NPY,NSL) + VCY/AFY(NPY)
              CC = CC + VCY
            ELSEIF( FNS.LT.ZERO .AND. J.LT.JFLD ) THEN
             IF( IXP(NN).GT.0 ) THEN
              FCNN = YN(NN,NSL)/(SN(2,NN)*PORD(2,NN)+SMALL)
              R = ((CO(N,NSL)*FCNP-CO(NN,NSL)*FCNN)
     &          /(CO(NS,NSL)*FCNS-CO(N,NSL)*FCNP+EPSL))
     &          *((DYGF(NS)+DYGF(N))/(DYGF(N)+DYGF(NN)))
              THETA = FLIMIT( R,CRNS,ISLC(1) )
              DYF = DYGF(N)/(DYGF(N)+DYGF(NS))
              VCY = CO(NS,NSL)*FNS*THETA*DYF*FCNS
     &          + CO(N,NSL)*FNS*(1.D+0-THETA*DYF)*FCNP
              VC(NPY,NSL) = VC(NPY,NSL) + VCY/AFY(NPY)
              CC = CC + VCY
             ELSE
              VCY = CO(N,NSL)*FNS*FCNP
              VC(NPY,NSL) = VC(NPY,NSL) + VCY/AFY(NPY)
              CC = CC + VCY
             ENDIF
            ELSEIF( FNS.LT.ZERO .AND. J.EQ.JFLD ) THEN
              VCY = CO(N,NSL)*FNS*FCNP
              VC(NPY,NSL) = VC(NPY,NSL) + VCY/AFY(NPY)
              CC = CC + VCY
            ENDIF
            AS = DNY*FCNS
            AP = DNY*FCNP
          ELSE
!
!---  Patankar solute transport  ---
!
            ANS = MAX(FNS,ZERO)
     &        + DNY*MAX((ONE-(TENTH*ABS(FNS)/(DNY+SMALL)))**5,ZERO)
            AP = (ANS-FNS)*FCNP
            AS = ANS*FCNS
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
        ENDIF
  120   CONTINUE
!
!---  West face diffusion and advection terms  ---
!
        IF( I.NE.1 ) THEN
          IF( IXP(NW).EQ.0 .OR. INBS(3,N).GT.0 ) GOTO 130
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
          IF( ISLC(1).GE.1 ) THEN
!
!---  TVD solute transport  ---
!
            IF( FNW.GE.ZERO .AND. I.GT.2 ) THEN
             IF( IXP(NW-1).GT.0 ) THEN
              FCNWW = YN(NW-1,NSL)/(SN(2,NW-1)*PORD(2,NW-1)+SMALL)
              R = ((CO(NW,NSL)*FCNW-CO(NW-1,NSL)*FCNWW)
     &          /(CO(N,NSL)*FCNP-CO(NW,NSL)*FCNW+EPSL))
     &          *((DXGF(N)+DXGF(NW))/(DXGF(NW)+DXGF(NWW)))
              THETA = FLIMIT( R,CRNW,ISLC(1) )
              DXF = DXGF(NW)/(DXGF(N)+DXGF(NW))
              UCX = CO(NW,NSL)*FNW*(1.D+0-THETA*DXF)*FCNW
     &          + CO(N,NSL)*FNW*THETA*DXF*FCNP
              UC(NPX,NSL) = UC(NPX,NSL) + UCX/AFX(NPX)
              CC = CC + UCX
             ELSE
              UCX = CO(NW,NSL)*FNW*FCNW
              UC(NPX,NSL) = UC(NPX,NSL) + UCX/AFX(NPX)
              CC = CC + UCX
             ENDIF
            ELSEIF( FNW.GE.ZERO .AND. I.EQ.2 ) THEN
              UCX = CO(NW,NSL)*FNW*FCNW
              UC(NPX,NSL) = UC(NPX,NSL) + UCX/AFX(NPX)
              CC = CC + UCX
            ELSEIF( FNW.LT.ZERO .AND. I.LT.IFLD ) THEN
             IF( IXP(NE).GT.0 ) THEN
              FCNE = YN(NE,NSL)/(SN(2,NE)*PORD(2,NE)+SMALL)
              R = ((CO(N,NSL)*FCNP-CO(NE,NSL)*FCNE)
     &          /(CO(NW,NSL)*FCNW-CO(N,NSL)*FCNP+EPSL))
     &          *((DXGF(NW)+DXGF(N))/(DXGF(N)+DXGF(NE)))
              THETA = FLIMIT( R,CRNW,ISLC(1) )
              DXF = DXGF(N)/(DXGF(N)+DXGF(NW))
              UCX = CO(NW,NSL)*FNW*THETA*DXF*FCNW
     &          + CO(N,NSL)*FNW*(1.D+0-THETA*DXF)*FCNP
              UC(NPX,NSL) = UC(NPX,NSL) + UCX/AFX(NPX)
              CC = CC + UCX
             ELSE
              UCX = CO(N,NSL)*FNW*FCNP
              UC(NPX,NSL) = UC(NPX,NSL) + UCX/AFX(NPX)
              CC = CC + UCX
             ENDIF
            ELSEIF( FNW.LT.ZERO .AND. I.EQ.IFLD ) THEN
              UCX = CO(N,NSL)*FNW*FCNP
              UC(NPX,NSL) = UC(NPX,NSL) + UCX/AFX(NPX)
              CC = CC + UCX
            ENDIF
            AW = DNX*FCNW
            AP = DNX*FCNP
          ELSE
!
!---  Patankar solute transport  ---
!
            ANW = MAX(FNW,ZERO)
     &        + DNX*MAX((ONE-(TENTH*ABS(FNW)/(DNX+SMALL)))**5,ZERO)
            AP = (ANW-FNW)*FCNP
            AW = ANW*FCNW
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
        ENDIF
  130   CONTINUE
!
!---  East face diffusion and advection terms  ---
!
        IF( I.NE.IFLD ) THEN
          IF( IXP(NE).EQ.0 .OR. INBS(4,N).GT.0 ) GOTO 140
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
          IF( ISLC(1).GE.1 ) THEN
!
!---  TVD solute transport  ---
!
            IF( FNE.GE.ZERO .AND. I.GT.1 ) THEN
             IF( IXP(NW).GT.0 ) THEN
              FCNW = YN(NW,NSL)/(SN(2,NW)*PORD(2,NW)+SMALL)
              R = ((CO(N,NSL)*FCNP-CO(NW,NSL)*FCNW)
     &          /(CO(NE,NSL)*FCNE-CO(N,NSL)*FCNP+EPSL))
     &          *((DXGF(NE)+DXGF(N))/(DXGF(N)+DXGF(NW)))
              THETA = FLIMIT( R,CRNE,ISLC(1) )
              DXF = DXGF(N)/(DXGF(N)+DXGF(NE))
              UCX = CO(NE,NSL)*FNE*THETA*DXF*FCNE
     &          + CO(N,NSL)*FNE*(1.D+0-THETA*DXF)*FCNP
              CC = CC - UCX
             ELSE
              UCX = CO(N,NSL)*FNE*FCNP
              CC = CC - UCX
             ENDIF
            ELSEIF( FNE.GE.ZERO .AND. I.EQ.1 ) THEN
              UCX = CO(N,NSL)*FNE*FCNP
              CC = CC - UCX
            ELSEIF( FNE.LT.ZERO .AND. I.LT.IFLD-1 ) THEN
             IF( IXP(NE+1).GT.0 ) THEN
              FCNEE = YN(NE+1,NSL)/(SN(2,NE+1)*PORD(2,NE+1)+SMALL)
              R = ((CO(NE,NSL)*FCNE-CO(NE+1,NSL)*FCNEE)
     &          /(CO(N,NSL)*FCNP-CO(NE,NSL)*FCNE+EPSL))
     &          *((DXGF(N)+DXGF(NE))/(DXGF(NE)+DXGF(NEE)))
              THETA = FLIMIT( R,CRNE,ISLC(1) )
              DXF = DXGF(NE)/(DXGF(N)+DXGF(NE))
              UCX = CO(NE,NSL)*FNE*(1.D+0-THETA*DXF)*FCNE
     &          + CO(N,NSL)*FNE*THETA*DXF*FCNP
              CC = CC - UCX
             ELSE
              UCX = CO(NE,NSL)*FNE*FCNE
              CC = CC - UCX
             ENDIF
            ELSEIF( FNE.LT.ZERO .AND. I.EQ.IFLD-1 ) THEN
              UCX = CO(NE,NSL)*FNE*FCNE
              CC = CC - UCX
            ENDIF
            AE = DNX*FCNE
            AP = DNX*FCNP
          ELSE
!
!---  Patankar solute transport  ---
!
            ANE = MAX(-FNE,ZERO)
     &        + DNX*MAX((ONE-(TENTH*ABS(FNE)/(DNX+SMALL)))**5,ZERO)
            AP = (ANE+FNE)*FCNP
            AE = ANE*FCNE
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
        ENDIF
  140   CONTINUE
!
!---  North face diffusion and advection terms  ---
!
        IF( J.NE.JFLD ) THEN
          IF( IXP(NN).EQ.0 .OR. INBS(5,N).GT.0 ) GOTO 150
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
          IF( ISLC(1).GE.1 ) THEN
!
!---  TVD solute transport  ---
!
            IF( FNN.GE.ZERO .AND. J.GT.1 ) THEN
             IF( IXP(NS).GT.0 ) THEN
              FCNS = YN(NS,NSL)/(SN(2,NS)*PORD(2,NS)+SMALL)
              R = ((CO(N,NSL)*FCNP-CO(NS,NSL)*FCNS)
     &          /(CO(NN,NSL)*FCNN-CO(N,NSL)*FCNP+EPSL))
     &          *((DYGF(NN)+DYGF(N))/(DYGF(N)+DYGF(NS)))
              THETA = FLIMIT( R,CRNN,ISLC(1) )
              DYF = DYGF(N)/(DYGF(N)+DYGF(NN))
              VCY = CO(NN,NSL)*FNN*THETA*DYF*FCNN
     &          + CO(N,NSL)*FNN*(1.D+0-THETA*DYF)*FCNP
              CC = CC - VCY
             ELSE
              VCY = CO(N,NSL)*FNN*FCNP
              CC = CC - VCY
             ENDIF
            ELSEIF( FNN.GE.ZERO .AND. J.EQ.1 ) THEN
              VCY = CO(N,NSL)*FNN*FCNP
              CC = CC - VCY
            ELSEIF( FNN.LT.ZERO .AND. J.LT.JFLD-1 ) THEN
             IF( IXP(NN+IFLD).GT.0 ) THEN
              FCNNN = YN(NN+IFLD,NSL)
     &          /(SN(2,NN+IFLD)*PORD(2,NN+IFLD)+SMALL)
              R = ((CO(NN,NSL)*FCNN-CO(NN+IFLD,NSL)*FCNNN)
     &          /(CO(N,NSL)*FCNP-CO(NN,NSL)*FCNN+EPSL))
     &          *((DYGF(N)+DYGF(NN))/(DYGF(NN)+DYGF(NNN)))
              THETA = FLIMIT( R,CRNN,ISLC(1) )
              DYF = DYGF(NN)/(DYGF(N)+DYGF(NN))
              VCY = CO(NN,NSL)*FNN*(1.D+0-THETA*DYF)*FCNN
     &          + CO(N,NSL)*FNN*THETA*DYF*FCNP
              CC = CC - VCY
             ELSE
              VCY = CO(NN,NSL)*FNN*FCNN
              CC = CC - VCY
             ENDIF
            ELSEIF( FNN.LT.ZERO .AND. J.EQ.JFLD-1 ) THEN
              VCY = CO(NN,NSL)*FNN*FCNN
              CC = CC - VCY
            ENDIF
            AN = DNY*FCNN
            AP = DNY*FCNP
          ELSE
!
!---  Patankar solute transport  ---
!
            ANN = MAX(-FNN,ZERO)
     &        + DNY*MAX((ONE-(TENTH*ABS(FNN)/(DNY+SMALL)))**5,ZERO)
            AP = (ANN+FNN)*FCNP
            AN = ANN*FCNN
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
        ENDIF
  150   CONTINUE
!
!---  Top face diffusion and advection terms  ---
!
        IF( K.NE.KFLD ) THEN
          IF( IXP(NT).EQ.0 .OR. INBS(6,N).GT.0 ) GOTO 160
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
          IF( ISLC(1).GE.1 ) THEN
!
!---  TVD solute transport  ---
!
            IF( FNT.GE.ZERO .AND. K.GT.1 ) THEN
             IF( IXP(NB).GT.0 ) THEN
              FCNB = YN(NB,NSL)/(SN(2,NB)*PORD(2,NB)+SMALL)
              R = ((CO(N,NSL)*FCNP-CO(NS,NSL)*FCNB)
     &          /(CO(NT,NSL)*FCNT-CO(N,NSL)*FCNP+EPSL))
     &          *((DZGF(NT)+DZGF(N))/(DZGF(N)+DZGF(NB)))
              THETA = FLIMIT( R,CRNT,ISLC(1) )
              DZF = DZGF(N)/(DZGF(N)+DZGF(NT))
              WCZ = CO(NT,NSL)*FNT*THETA*DZF*FCNT
     &          + CO(N,NSL)*FNT*(1.D+0-THETA*DZF)*FCNP
              CC = CC - WCZ
             ELSE
              WCZ = CO(N,NSL)*FNT*FCNP
              CC = CC - WCZ
             ENDIF
            ELSEIF( FNT.GE.ZERO .AND. K.EQ.1 ) THEN
              WCZ = CO(N,NSL)*FNT*FCNP
              CC = CC - WCZ
            ELSEIF( FNT.LT.ZERO .AND. K.LT.KFLD-1 ) THEN
             IF( IXP(NT+IJFLD).GT.0 ) THEN
              FCNTT = YN(NT+IJFLD,NSL)
     &          /(SN(2,NT+IJFLD)*PORD(2,NT+IJFLD)+SMALL)
              R = ((CO(NT,NSL)*FCNT-CO(NT+IJFLD,NSL)*FCNTT)
     &          /(CO(N,NSL)*FCNP-CO(NT,NSL)*FCNT+EPSL))
     &          *((DZGF(N)+DZGF(NT))/(DZGF(NT)+DZGF(NTT)))
              THETA = FLIMIT( R,CRNT,ISLC(1) )
              DZF = DZGF(NT)/(DZGF(N)+DZGF(NT))
              WCZ = CO(NT,NSL)*FNT*(1.D+0-THETA*DZF)*FCNT
     &          + CO(N,NSL)*FNT*THETA*DZF*FCNP
              CC = CC - WCZ
             ELSE
              WCZ = CO(NT,NSL)*FNT*FCNT
              CC = CC - WCZ
             ENDIF
            ELSEIF( FNT.LT.ZERO .AND. K.EQ.KFLD-1 ) THEN
              WCZ = CO(NT,NSL)*FNT*FCNT
              CC = CC - WCZ
            ENDIF
            AT = DNZ*FCNT
            AP = DNZ*FCNP
          ELSE
!
!---  Patankar solute transport  ---
!
            ANT = MAX(-FNT,ZERO)
     &        + DNZ*MAX((ONE-(TENTH*ABS(FNT)/(DNZ+SMALL)))**5,ZERO)
            AP = (ANT+FNT)*FCNP
            AT = ANT*FCNT
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
        ENDIF
  160   CONTINUE
!
!---  Solution vector and parent radioactive decay  ---
!
        IF( IEQW.EQ.0 .AND. IEQA.EQ.0 ) THEN
!
!---      Solution vector and parent radioactive decay  ---
!
          BLU(MCP) = BLU(MCP) + CO(N,NSL)*SC
     &      - 6.931D-1*CO(N,NSL)*VOL(N)/HLF(NSL)
!          BLU(MCP) = BLU(MCP) + CO(N,NSL)*(2.D+0**(-DT/HLF(NSL)))*SC

!---      Daughter radioactive chain decay  ---
!
          IF( NSL.LE.NSOLU ) THEN
            DO 170 NPSL = 1,NSOLU
              IF( NPSL.EQ.NSL ) GOTO 170
              BLU(MCP) = BLU(MCP) +
     &          CHDF(NPSL,NSL)*6.931D-1*CO(N,NPSL)*VOL(N)/HLF(NPSL)
!              BLU(MCP) = BLU(MCP) +
!     &          CHDF(NPSL,NSL)*CO(N,NPSL)*(2.D+0**(-DT/HLF(NPSL)))*SC
  170       CONTINUE
          ENDIF
        ENDIF
        IF( ISLC(1).GE.1 ) BLU(MCP) = BLU(MCP) + CC
  900 CONTINUE
!
!---  End of SJCBN group  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END


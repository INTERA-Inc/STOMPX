!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBA_GT
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
!     Load the Jacobian matrix for the air equation with
!     aqueous-phase and gas-phase contributions
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
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
      USE FDVP
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STAX(LUK+1),RAP(LUK),RAA(LUK,6),FA(LSFV,6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBA_GT'
      IEQAX = IEQA
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO N = 1,NFLD
        I = ID(N)
        J = JD(N)
        K = KD(N)
        IF( IXP(N).EQ.0 ) CYCLE
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
!
!---    First-order, forward-difference, time differential  ---
!
        STA1 = PORD(1,N)*(XLA(1,N)*RHOL(1,N)*SL(1,N) +
     &      XGA(1,N)*RHOG(1,N)*SG(1,N))
        DO M = 1,ISVC+1
          MP = M + 1
          STA0 = PORD(MP,N)*(XLA(MP,N)*RHOL(MP,N)*SL(MP,N) +
     &      XGA(MP,N)*RHOG(MP,N)*SG(MP,N))
          STAX(M) = (STA0-STA1)*DTI*VOL(N)
        ENDDO
!
!---    Initialize surface fluxes  ---
!
        DO MD = 1,6
          DO M = 1,ISVF
            FA(M,MD) = 0.D+0
          ENDDO
        ENDDO
!
!---    Compute surface fluxes  ---
!
!---    Bottom ---
!
        IF( K.NE.1 ) THEN
          NB = N-IJFLD
          IF( IXP(NB).EQ.0 .OR. INBS(1,N).GT.0 ) GOTO 310
          DO M = 1,ISVF
            MB = MADJ(M)
            MP = MNOD(M)
            FLAB = XLA(MB,NB)*RHOL(MB,NB)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            INDX = 2
            FLA = DIFMN( FLAB,FLAP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
            FGAB = XGA(MB,NB)*RHOG(MB,NB)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            INDX = 3
            FGA = DIFMN( FGAB,FGAP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
            FA(M,1) = -AFZ(NPZ)*(WL(M,NPZ)*FLA + WG(M,NPZ)*FGA -
     &        WTMA*WDGW(M,NPZ) + WTMA*WDLA(M,NPZ))
          ENDDO
        ENDIF
  310   CONTINUE
!
!---    South ---
!
        IF( J.NE.1 ) THEN
          NS = N-IFLD
          IF( IXP(NS).EQ.0 .OR. INBS(2,N).GT.0 ) GOTO 410
          DO M = 1,ISVF
            MS = MADJ(M)
            MP = MNOD(M)
            FLAS = XLA(MS,NS)*RHOL(MS,NS)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            INDX = 2
            FLA = DIFMN( FLAS,FLAP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
            FGAS = XGA(MS,NS)*RHOG(MS,NS)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            INDX = 3
            FGA = DIFMN( FGAS,FGAP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
            FA(M,2) = -AFY(NPY)*(VL(M,NPY)*FLA + VG(M,NPY)*FGA -
     &        WTMA*VDGW(M,NPY) + WTMA*VDLA(M,NPY))
          ENDDO
        ENDIF
  410   CONTINUE
!
!---    West ---
!
        IF( I.NE.1 ) THEN
          NW = N-1
          IF( IXP(NW).EQ.0 .OR. INBS(3,N).GT.0 ) GOTO 510
          DO M = 1,ISVF
            MW = MADJ(M)
            MP = MNOD(M)
            FLAW = XLA(MW,NW)*RHOL(MW,NW)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            INDX = 2
            FLA = DIFMN( FLAW,FLAP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            FGAW = XGA(MW,NW)*RHOG(MW,NW)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            INDX = 3
            FGA = DIFMN( FGAW,FGAP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            FA(M,3) = -AFX(NPX)*(UL(M,NPX)*FLA + UG(M,NPX)*FGA -
     &        WTMA*UDGW(M,NPX) + WTMA*UDLA(M,NPX))
          ENDDO
        ENDIF
  510   CONTINUE
!
!---    East ---
!
        IF( I.NE.IFLD ) THEN
          NE = N+1
          IF( IXP(NE).EQ.0 .OR. INBS(4,N).GT.0 ) GOTO 610
          DO M = 1,ISVF
            ME = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FLAE = XLA(ME,NE)*RHOL(ME,NE)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            INDX = 2
            FLA = DIFMN( FLAP,FLAE,DXGF(N),DXGF(NE),UL(1,NQX),INDX )
            FGAE = XGA(ME,NE)*RHOG(ME,NE)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            INDX = 3
            FGA = DIFMN( FGAP,FGAE,DXGF(N),DXGF(NE),UG(1,NQX),INDX )
            FA(M,4) = AFX(NQX)*(UL(MF,NQX)*FLA + UG(MF,NQX)*FGA -
     &        WTMA*UDGW(MF,NQX) + WTMA*UDLA(MF,NQX))
          ENDDO
        ENDIF
  610   CONTINUE
!
!---    North ---
!
        IF( J.NE.JFLD ) THEN
          NN = N+IFLD
          IF( IXP(NN).EQ.0 .OR. INBS(5,N).GT.0 ) GOTO 710
          DO M = 1,ISVF
            MN = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FLAN = XLA(MN,NN)*RHOL(MN,NN)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            INDX = 2
            FLA = DIFMN( FLAP,FLAN,DYGF(N),DYGF(NN),VL(1,NQY),INDX )
            FGAN = XGA(MN,NN)*RHOG(MN,NN)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            INDX = 3
            FGA = DIFMN( FGAP,FGAN,DYGF(N),DYGF(NN),VG(1,NQY),INDX )
            FA(M,5) = AFY(NQY)*(VL(MF,NQY)*FLA + VG(MF,NQY)*FGA -
     &        WTMA*VDGW(MF,NQY) + WTMA*VDLA(MF,NQY))
          ENDDO
        ENDIF
  710   CONTINUE
!
!---    Top ---
!
        IF( K.NE.KFLD ) THEN
          NT = N+IJFLD
          IF( IXP(NT).EQ.0 .OR. INBS(6,N).GT.0 ) GOTO 810
          DO M = 1,ISVF
            MT = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FLAT = XLA(MT,NT)*RHOL(MT,NT)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            INDX = 2
            FLA = DIFMN( FLAP,FLAT,DZGF(N),DZGF(NT),WL(1,NQZ),INDX )
            FGAT = XGA(MT,NT)*RHOG(MT,NT)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            INDX = 3
            FGA = DIFMN( FGAP,FGAT,DZGF(N),DZGF(NT),WG(1,NQZ),INDX )
            FA(M,6) = AFZ(NQZ)*(WL(MF,NQZ)*FLA + WG(MF,NQZ)*FGA -
     &        WTMA*WDGW(MF,NQZ) + WTMA*WDLA(MF,NQZ))
          ENDDO
        ENDIF
  810   CONTINUE
!
!---    Compute air equation residuals  ---
!
        RAS = STAX(1) - SRCA(2,N)
        DO MD = 1,6
          RAS = RAS + FA(1,MD)
        ENDDO
        DO M = 1,ISVC
          RAP(M) = STAX(M+1) - SRCA(M+2,N)
          MM = 2*M
          DO MD = 1,6
            RAP(M) = RAP(M) + FA(MM,MD)
          ENDDO
        ENDDO
        DO M = 1,ISVC
          MM = 2*M + 1
          DO MD = 1,6
            RAA(M,MD) = RAS - FA(1,MD) + FA(MM,MD)
          ENDDO
        ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        RAS = 0.D+0
        DO M = 1,ISVC
          RAP(M) = 0.D+0
          DO MD = 1,6
            RAA(M,MD) = 0.D+0
          ENDDO
        ENDDO
        RAP(IEQAX) = 1.D+0
      ENDIF
!
!---    Load Jacobian Matrix  ---
!
        CALL JCBL_GT( RAS,RAP,RAA,N,I,J,K,IEQAX )
!
!---  Continue to Next Node  ---
!
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBA_GT group
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBL_GT( RSS,RSP,RSA,N,I,J,K,MEQ )
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
!     Load the Jacobian matrix.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
#ifdef petsc
!----------------------PETSc Modules-----------------------------------!
!
      USE PETSC_STOMP
!
#endif
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#ifdef petsc
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#endif
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RSP(LUK),RSA(LUK,6)
#ifdef petsc
      REAL*8 BUFFER(LUK)
      INTEGER IROW(1), ICOL(LUK)
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBL_GT'
!
!---  Banded solver  ---
!
      IF( ILES.EQ.1 ) THEN
!
!---    Node  ---
!
        NMD = IXP(N)
        MP = IM(MEQ,NMD)
        DO M = 1,ISVC
          MCOL = IM(M,NMD)
          MROW = MP-MCOL+MDC
          ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSP(M)-RSS)/DNR(M,N)
        ENDDO
        BLU(MP) = BLU(MP) - RSS
        RSDL(MEQ,N) = BLU(MP)
!
!---    Bottom ---
!
        IF( K.NE.1 ) THEN
          NB = N-IJFLD
          IF( IXP(NB).EQ.0 .OR. INBS(1,N).GT.0 ) GOTO 210
          NMD = ABS(IXP(NB))
          DO M = 1,ISVC
            DNRX = DNR(M,NB)
            MCOL = IM(M,NMD)
            MROW = MP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSA(M,1)-RSS)/DNRX
          ENDDO
        ENDIF
  210   CONTINUE
!
!---    South ---
!
        IF( J.NE.1 ) THEN
          NS = N-IFLD
          IF( IXP(NS).EQ.0 .OR. INBS(2,N).GT.0 ) GOTO 310
          NMD = ABS(IXP(NS))
          DO M = 1,ISVC
            DNRX = DNR(M,NS)
            MCOL = IM(M,NMD)
            MROW = MP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSA(M,2)-RSS)/DNRX
          ENDDO
        ENDIF
  310   CONTINUE
!
!---    West ---
!
        IF( I.NE.1 ) THEN
          NW = N-1
          IF( IXP(NW).EQ.0 .OR. INBS(3,N).GT.0 ) GOTO 410
          NMD = ABS(IXP(NW))
          DO M = 1,ISVC
            DNRX = DNR(M,NW)
            MCOL = IM(M,NMD)
            MROW = MP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSA(M,3)-RSS)/DNRX
          ENDDO
        ENDIF
  410   CONTINUE
!
!---    East ---
!
        IF( I.NE.IFLD ) THEN
          NE = N+1
          IF( IXP(NE).EQ.0 .OR. INBS(4,N).GT.0 ) GOTO 510
          NMD = ABS(IXP(NE))
          DO M = 1,ISVC
            DNRX = DNR(M,NE)
            MCOL = IM(M,NMD)
            MROW = MP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSA(M,4)-RSS)/DNRX
          ENDDO
        ENDIF
  510   CONTINUE
!
!---    North ---
!
        IF( J.NE.JFLD ) THEN
          NN = N+IFLD
          IF( IXP(NN).EQ.0 .OR. INBS(5,N).GT.0 ) GOTO 610
          NMD = ABS(IXP(NN))
          DO M = 1,ISVC
            DNRX = DNR(M,NN)
            MCOL = IM(M,NMD)
            MROW = MP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSA(M,5)-RSS)/DNRX
          ENDDO
        ENDIF
  610   CONTINUE
!
!---    Top ---
!
        IF( K.NE.KFLD ) THEN
          NT = N+IJFLD
          IF( IXP(NT).EQ.0 .OR. INBS(6,N).GT.0 ) GOTO 710
          NMD = ABS(IXP(NT))
          DO M = 1,ISVC
            DNRX = DNR(M,NT)
            MCOL = IM(M,NMD)
            MROW = MP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RSA(M,6)-RSS)/DNRX
          ENDDO
        ENDIF
  710   CONTINUE
!
!---  SPLib solver  ---
!
      ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---    Node  ---
!
        NMD = IXP(N)
        MP = IM(MEQ,NMD)
        MA = 0
        DO M = 1,ISVC
          MCOL = KLU(MP,M+MA)
          DLU(MCOL) = DLU(MCOL) + (RSP(M)-RSS)/DNR(M,N)
        ENDDO
        BLU(MP) = BLU(MP) - RSS
        RSDL(MEQ,N) = BLU(MP)
        MA = MA + ISVC
!
!---    Bottom ---
!
        IF( K.NE.1 ) THEN
          NB = N-IJFLD
          IF( IXP(NB).EQ.0 .OR. INBS(1,N).GT.0 ) GOTO 2210
          NMD = ABS(IXP(NB))
          DO M = 1,ISVC
            DNRX = DNR(M,NB)
            MCOL = KLU(MP,M+MA)
            DLU(MCOL) = DLU(MCOL) + (RSA(M,1)-RSS)/DNRX
          ENDDO
          MA = MA + ISVC
        ENDIF
 2210   CONTINUE
!
!---    South ---
!
        IF( J.NE.1 ) THEN
          NS = N-IFLD
          IF( IXP(NS).EQ.0 .OR. INBS(2,N).GT.0 ) GOTO 2310
          NMD = ABS(IXP(NS))
          DO M = 1,ISVC
            DNRX = DNR(M,NS)
            MCOL = KLU(MP,M+MA)
            DLU(MCOL) = DLU(MCOL) + (RSA(M,2)-RSS)/DNRX
          ENDDO
          MA = MA + ISVC
        ENDIF
 2310   CONTINUE
!
!---    West ---
!
        IF( I.NE.1 ) THEN
          NW = N-1
          IF( IXP(NW).EQ.0 .OR. INBS(3,N).GT.0 ) GOTO 2410
          NMD = ABS(IXP(NW))
          DO M = 1,ISVC
            DNRX = DNR(M,NW)
            MCOL = KLU(MP,M+MA)
            DLU(MCOL) = DLU(MCOL) + (RSA(M,3)-RSS)/DNRX
          ENDDO
          MA = MA + ISVC
        ENDIF
 2410   CONTINUE
!
!---    East ---
!
        IF( I.NE.IFLD ) THEN
          NE = N+1
          IF( IXP(NE).EQ.0 .OR. INBS(4,N).GT.0 ) GOTO 2510
          NMD = IXP(NE)
          DO M = 1,ISVC
            DNRX = DNR(M,NE)
            MCOL = KLU(MP,M+MA)
            DLU(MCOL) = DLU(MCOL) + (RSA(M,4)-RSS)/DNRX
          ENDDO
          MA = MA + ISVC
        ENDIF
 2510   CONTINUE
!
!---    North ---
!
        IF( J.NE.JFLD ) THEN
          NN = N+IFLD
          IF( IXP(NN).EQ.0 .OR. INBS(5,N).GT.0 ) GOTO 2610
          NMD = ABS(IXP(NN))
          DO M = 1,ISVC
            DNRX = DNR(M,NN)
            MCOL = KLU(MP,M+MA)
            DLU(MCOL) = DLU(MCOL) + (RSA(M,5)-RSS)/DNRX
          ENDDO
          MA = MA + ISVC
        ENDIF
 2610   CONTINUE
!
!---    Top ---
!
        IF( K.NE.KFLD ) THEN
          NT = N+IJFLD
          IF( IXP(NT).EQ.0 .OR. INBS(6,N).GT.0 ) GOTO 2710
          NMD = ABS(IXP(NT))
          DO M = 1,ISVC
            DNRX = DNR(M,NT)
            MCOL = KLU(MP,M+MA)
            DLU(MCOL) = DLU(MCOL) + (RSA(M,6)-RSS)/DNRX
          ENDDO
          MA = MA + ISVC
        ENDIF
 2710   CONTINUE
#ifdef petsc
!
!---  PETSc solver  ---
!
      ELSEIF( ILES.EQ.5 ) THEN
!
!---    Node  ---
!
        NMD = IXP(N)
        MP = IM(MEQ,NMD)
        IROW(1) = MP-1
        DO M = 1,ISVC
          ICOL(M) = IM(M,NMD)-1
          BUFFER(M) = (RSP(M)-RSS)/DNR(M,N)
        ENDDO
        CALL MatSetValues(F_MAT,1,IROW,ISVC,ICOL,BUFFER,
     &                    ADD_VALUES,IERR)
        BLU(MP) = BLU(MP) - RSS
        RSDL(MEQ,N) = BLU(MP)
!
!---    Bottom ---
!
        IF( K.NE.1 ) THEN
          NB = N-IJFLD
          IF( IXP(NB).EQ.0 .OR. INBS(1,N).GT.0 ) GOTO 3210
          NMD = ABS(IXP(NB))
          DO M = 1,ISVC
            DNRX = DNR(M,NB)
            ICOL(M) = IM(M,NMD)-1
            BUFFER(M) = (RSA(M,1)-RSS)/DNRX
          ENDDO
          CALL MatSetValues(F_MAT,1,IROW,ISVC,ICOL,BUFFER,
     &                      ADD_VALUES,IERR)
        ENDIF
 3210   CONTINUE
!
!---    South ---
!
        IF( J.NE.1 ) THEN
          NS = N-IFLD
          IF( IXP(NS).EQ.0 .OR. INBS(2,N).GT.0 ) GOTO 3310
          NMD = ABS(IXP(NS))
          DO M = 1,ISVC
            DNRX = DNR(M,NS)
            ICOL(M) = IM(M,NMD)-1
            BUFFER(M) = (RSA(M,2)-RSS)/DNRX
          ENDDO
          CALL MatSetValues(F_MAT,1,IROW,ISVC,ICOL,BUFFER,
     &                      ADD_VALUES,IERR)
        ENDIF
 3310   CONTINUE
!
!---    West ---
!
        IF( I.NE.1 ) THEN
          NW = N-1
          IF( IXP(NW).EQ.0 .OR. INBS(3,N).GT.0 ) GOTO 3410
          NMD = ABS(IXP(NW))
          DO M = 1,ISVC
            DNRX = DNR(M,NW)
            ICOL(M) = IM(M,NMD)-1
            BUFFER(M) = (RSA(M,3)-RSS)/DNRX
          ENDDO
          CALL MatSetValues(F_MAT,1,IROW,ISVC,ICOL,BUFFER,
     &                      ADD_VALUES,IERR)
        ENDIF
 3410   CONTINUE
!
!---    East ---
!
        IF( I.NE.IFLD ) THEN
          NE = N+1
          IF( IXP(NE).EQ.0 .OR. INBS(4,N).GT.0 ) GOTO 3510
          NMD = IXP(NE)
          DO M = 1,ISVC
            DNRX = DNR(M,NE)
            ICOL(M) = IM(M,NMD)-1
            BUFFER(M) = (RSA(M,4)-RSS)/DNRX
          ENDDO
          CALL MatSetValues(F_MAT,1,IROW,ISVC,ICOL,BUFFER,
     &                      ADD_VALUES,IERR)
        ENDIF
 3510   CONTINUE
!
!---    North ---
!
        IF( J.NE.JFLD ) THEN
          NN = N+IFLD
          IF( IXP(NN).EQ.0 .OR. INBS(5,N).GT.0 ) GOTO 3610
          NMD = ABS(IXP(NN))
          DO M = 1,ISVC
            DNRX = DNR(M,NN)
            ICOL(M) = IM(M,NMD)-1
            BUFFER(M) = (RSA(M,5)-RSS)/DNRX
          ENDDO
          CALL MatSetValues(F_MAT,1,IROW,ISVC,ICOL,BUFFER,
     &                      ADD_VALUES,IERR)
        ENDIF
 3610   CONTINUE
!
!---    Top ---
!
        IF( K.NE.KFLD ) THEN
          NT = N+IJFLD
          IF( IXP(NT).EQ.0 .OR. INBS(6,N).GT.0 ) GOTO 3710
          NMD = ABS(IXP(NT))
          DO M = 1,ISVC
            DNRX = DNR(M,NT)
            ICOL(M) = IM(M,NMD)-1
            BUFFER(M) = (RSA(M,6)-RSS)/DNRX
          ENDDO
          CALL MatSetValues(F_MAT,1,IROW,ISVC,ICOL,BUFFER,
     &                      ADD_VALUES,IERR)
        ENDIF
 3710   CONTINUE
#endif
      ELSE
        INDX = 3
        CHMSG = 'Unknown Linear Equation Solver'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBL_GT group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLB_GT( RS,N,MEQ )
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
!     Modify the Jacobian matrix for boundary conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
#ifdef petsc
!----------------------PETSc Modules-----------------------------------!
!
      USE PETSC_STOMP
!
#endif
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE JACOB
      USE GRID
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#ifdef petsc
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#endif
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
#ifdef petsc
      REAL*8 BUFFER(LUK)
      INTEGER IROW(1), ICOL(LUK)
#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLB_GT'
      NMD = ABS(IXP(N))
      MP = IM(MEQ,NMD)
      MA = 0
!
!---  Banded solver  ---
!
      IF( ILES.EQ.1 ) THEN
        DO M = 1,ISVC
          MCOL = IM(M,NMD)
          MROW = MP-MCOL+MDC
          DNRX = DNR(M,N)
          ALU(MROW,MCOL) = ALU(MROW,MCOL) + (RS(M+1)-RS(1))/DNRX
        ENDDO
!
!---  SPLib solver  ---
!
      ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
        DO M = 1,ISVC
          MCOL = KLU(MP,M+MA)
          DNRX = DNR(M,N)
          DLU(MCOL) = DLU(MCOL) + (RS(M+1)-RS(1))/DNRX
        ENDDO
#ifdef petsc
!
!---  PETSc solver  ---
!
      ELSEIF( ILES.EQ.5 ) THEN
        IROW(1) = IM(MEQ,NMD)-1
        DO M = 1,ISVC
          DNRX = DNR(M,N)
          ICOL(M) = IM(M,NMD)-1
          BUFFER(M) = (RS(M+1)-RS(1))/DNRX
        ENDDO
        CALL MatSetValues(F_MAT,1,IROW,ISVC,ICOL,BUFFER,
     &                    ADD_VALUES,IERR)
#endif
!
!---  Unknown solver  ---
!
      ELSE
        INDX = 3
        CHMSG = 'Unknown Linear Equation Solver'
        CALL WRMSGS( INDX )
      ENDIF
      BLU(MP) = BLU(MP) - RS(1)
      RSDL(MEQ,N) = BLU(MP)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLB_GT group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBS_GT
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
!     Load the Jacobian matrix for the salt equation aqueous-phase
!     contributions (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STSX(LUK+1),RSP(LUK),RSA(LUK,6),FS(LSFV,6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBS_GT'
      IEQSX = IEQS
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO N = 1,NFLD
        I = ID(N)
        J = JD(N)
        K = KD(N)
        IF( IXP(N).EQ.0 ) CYCLE
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
!
!---    First-order, forward-difference, time differential  ---
!
        DO M = 1,ISVC+1
          MP = M + 1
          STSX(M) = (TMS(MP,N)-TMS(1,N))*VOL(N)*DTI
        ENDDO
!
!---  Initialize surface fluxes  ---
!
        DO MD = 1,6
          DO M = 1,ISVF
            FS(M,MD) = 0.D+0
          ENDDO
        ENDDO
!
!---  Bottom surface flux component  ---
!
        IF( K.NE.1 ) THEN
          NB = N-IJFLD
          IF( IXP(NB).EQ.0 .OR. INBS(1,N).GT.0 ) GOTO 110
          DO M = 1,ISVF
            FS(M,1) = -AFZ(NPZ)*WS(M,NPZ)
          ENDDO
        ENDIF
  110   CONTINUE
!
!---  South surface flux component ---
!
        IF( J.NE.1 ) THEN
          NS = N-IFLD
          IF( IXP(NS).EQ.0 .OR. INBS(2,N).GT.0 ) GOTO 210
          DO M = 1,ISVF
            FS(M,2) = -AFY(NPY)*VS(M,NPY)
          ENDDO
        ENDIF
  210   CONTINUE
!
!---  West surface flux component ---
!
        IF( I.NE.1 ) THEN
          NW = N-1
          IF( IXP(NW).EQ.0 .OR. INBS(3,N).GT.0 ) GOTO 310
          DO M = 1,ISVF
            FS(M,3) = -AFX(NPX)*US(M,NPX)
          ENDDO
        ENDIF
  310   CONTINUE
!
!---  East surface flux component  ---
!
        IF( I.NE.IFLD ) THEN
          NE = N+1
          IF( IXP(NE).EQ.0 .OR. INBS(4,N).GT.0 ) GOTO 410
          DO M = 1,ISVF
            MF = MFLX(M)
            FS(M,4) = AFX(NQX)*US(MF,NQX)
          ENDDO
        ENDIF
  410   CONTINUE
!
!---  North surface flux component  ---
!
        IF( J.NE.JFLD ) THEN
          NN = N+IFLD
          IF( IXP(NN).EQ.0 .OR. INBS(5,N).GT.0 ) GOTO 510
          DO M = 1,ISVF
            MF = MFLX(M)
            FS(M,5) = AFY(NQY)*VS(MF,NQY)
          ENDDO
        ENDIF
  510   CONTINUE
!
!---  Top surface flux component  ---
!
        IF( K.NE.KFLD ) THEN
          NT = N+IJFLD
          IF( IXP(NT).EQ.0 .OR. INBS(6,N).GT.0 ) GOTO 610
          DO M = 1,ISVF
            MF = MFLX(M)
            FS(M,6) = AFZ(NQZ)*WS(MF,NQZ)
          ENDDO
        ENDIF
  610   CONTINUE
!
!---  Salt equation residuals  ---
!
        RSS = STSX(1) - SRCS(2,N)
        DO MD = 1,6
          RSS = RSS + FS(1,MD)
        ENDDO
        DO M = 1,ISVC
          RSP(M) = STSX(M+1) - SRCS(M+2,N)
          MM = 2*M
          DO MD = 1,6
            RSP(M) = RSP(M) + FS(MM,MD)
          ENDDO
        ENDDO
        DO M = 1,ISVC
          MM = 2*M + 1
          DO MD = 1,6
            RSA(M,MD) = RSS - FS(1,MD) + FS(MM,MD)
          ENDDO
        ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        RSS = 0.D+0
        DO M = 1,ISVC
          RSP(M) = 0.D+0
          DO MD = 1,6
            RSA(M,MD) = 0.D+0
          ENDDO
        ENDDO
        RSP(IEQSX) = 1.D+0
      ENDIF
!
!---  Jacobian matrix loader  --
!
        CALL JCBL_GT( RSS,RSP,RSA,N,I,J,K,IEQSX )
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBS_GT group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBT_GT
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
!     Load the Jacobian matrix for the energy equation with
!     aqueous-phase and gas-phase contributions
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVS
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STTX(LUK+1),RTP(LUK),RTA(LUK,6),FT(LSFV,6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBT_GT'
      IEQAX = IEQA
      IEQSX = IEQS
      IEQTX = IEQT
      IEQWX = IEQW
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO N = 1,NFLD
        I = ID(N)
        J = JD(N)
        K = KD(N)
        IF( IXP(N).EQ.0 ) CYCLE
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
        IZN = IZ(N)
!
!---    First-order, forward-difference, time differential  ---
!
        STT1 = (1.D+0-POR(2,IZN))*RHOS(IZN)*CPS(IZN)*T(1,N) +
     &    PORD(1,N)*SS(1,N)*RHOSP(1,N)*HSP(1,N) +
     &    (PORT(1,N)-PORD(1,N))*RHOL(1,N)*HL(1,N) +
     &    PORD(1,N)*(SL(1,N)*RHOL(1,N)*HL(1,N) +
     &    SG(1,N)*RHOG(1,N)*UEG(1,N))
        DO M = 1,ISVC+1
          MP = M + 1
          STT0 = (1.D+0-POR(2,IZN))*RHOS(IZN)*CPS(IZN)*T(MP,N) +
     &      PORD(MP,N)*SS(MP,N)*RHOSP(MP,N)*HSP(MP,N) +
     &      (PORT(MP,N)-PORD(MP,N))*RHOL(MP,N)*HL(MP,N) +
     &      PORD(MP,N)*(SL(MP,N)*RHOL(MP,N)*HL(MP,N) +
     &      SG(MP,N)*RHOG(MP,N)*UEG(MP,N))
          STTX(M) = (STT0-STT1)*DTI*VOL(N)
        ENDDO
!
!---  Initialize surface fluxes  ---
!
        DO MD = 1,6
          DO M = 1,ISVF
            FT(M,MD) = 0.D+0
          ENDDO
        ENDDO
!
!---  Compute surface fluxes  ---
!
!---  Bottom ---
!
        IF( K.NE.1 ) THEN
          NB = N-IJFLD
          IF( IXP(NB).EQ.0 .OR. INBS(1,N).GT.0 ) GOTO 310
          DO M = 1,ISVF
            FT(M,1) = -AFZ(NPZ)*WQ(M,NPZ)
          ENDDO
        ENDIF
  310   CONTINUE
!
!---  South ---
!
        IF( J.NE.1 ) THEN
          NS = N-IFLD
          IF( IXP(NS).EQ.0 .OR. INBS(2,N).GT.0 ) GOTO 410
          DO M = 1,ISVF
            FT(M,2) = -AFY(NPY)*VQ(M,NPY)
          ENDDO
        ENDIF
  410   CONTINUE
!
!---  West ---
!
        IF( I.NE.1 ) THEN
          NW = N-1
          IF( IXP(NW).EQ.0 .OR. INBS(3,N).GT.0 ) GOTO 510
          DO M = 1,ISVF
            FT(M,3) = -AFX(NPX)*UQ(M,NPX)
          ENDDO
        ENDIF
  510   CONTINUE
!
!---  East ---
!
        IF( I.NE.IFLD ) THEN
          NE = N+1
          IF( IXP(NE).EQ.0 .OR. INBS(4,N).GT.0 ) GOTO 610
          DO M = 1,ISVF
            MF = MFLX(M)
            FT(M,4) = AFX(NQX)*UQ(MF,NQX)
          ENDDO
        ENDIF
  610   CONTINUE
!
!---  North ---
!
        IF( J.NE.JFLD ) THEN
          NN = N+IFLD
          IF( IXP(NN).EQ.0 .OR. INBS(5,N).GT.0 ) GOTO 710
          DO M = 1,ISVF
            MF = MFLX(M)
            FT(M,5) = AFY(NQY)*VQ(MF,NQY)
          ENDDO
        ENDIF
  710   CONTINUE
!
!---  Top ---
!
        IF( K.NE.KFLD ) THEN
          NT = N+IJFLD
          IF( IXP(NT).EQ.0 .OR. INBS(6,N).GT.0 ) GOTO 810
          DO M = 1,ISVF
            MF = MFLX(M)
            FT(M,6) = AFZ(NQZ)*WQ(MF,NQZ)
          ENDDO
        ENDIF
  810   CONTINUE
!
!---  Compute energy equation residuals  ---
!
        RTS = STTX(1) - SRCT(2,N)
        DO MD = 1,6
          RTS = RTS + FT(1,MD)
        ENDDO
        DO M = 1,ISVC
          RTP(M) = STTX(M+1) - SRCT(M+2,N)
          MM = 2*M
          DO MD = 1,6
            RTP(M) = RTP(M) + FT(MM,MD)
          ENDDO
        ENDDO
        DO M = 1,ISVC
          MM = 2*M + 1
          DO MD = 1,6
            RTA(M,MD) = RTS - FT(1,MD) + FT(MM,MD)
          ENDDO
        ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        IF( IEQWX.NE.0 ) RTP(IEQWX) = RTS
        IF( IEQAX.NE.0 ) RTP(IEQAX) = RTS
        IF( IEQSX.NE.0 ) RTP(IEQSX) = RTS
        DO MD = 1,6
          IF( IEQWX.NE.0 ) RTA(IEQWX,MD) = RTS
          IF( IEQAX.NE.0 ) RTA(IEQAX,MD) = RTS
          IF( IEQSX.NE.0 ) RTA(IEQSX,MD) = RTS
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
        CALL JCBL_GT( RTS,RTP,RTA,N,I,J,K,IEQTX )
!
!---  Continue to Next Node  ---
!
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBT_GT group
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBW_GT
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
!     Load the Jacobian matrix for the water equation with
!     aqueous-phase, gas-phase and salt contributions
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
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
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/JCBW_GT'
      IEQWX = IEQW
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO N = 1,NFLD
        I = ID(N)
        J = JD(N)
        K = KD(N)
        IF( IXP(N).EQ.0 ) CYCLE
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
!
!---    First-order, forward-difference, time differential  ---
!
        STW1 = PORD(1,N)*(XLW(1,N)*RHOL(1,N)*SL(1,N) +
     &      XGW(1,N)*RHOG(1,N)*SG(1,N))
        DO M = 1,ISVC+1
          MP = M + 1
          STW0 = PORD(MP,N)*(XLW(MP,N)*RHOL(MP,N)*SL(MP,N) +
     &      XGW(MP,N)*RHOG(MP,N)*SG(MP,N))
          STWX(M) = (STW0-STW1)*DTI*VOL(N)
        ENDDO
!
!---  Initialize surface fluxes  ---
!
        DO MD = 1,6
          DO M = 1,ISVF
            FW(M,MD) = 0.D+0
          ENDDO
        ENDDO
!
!---  Compute surface fluxes  ---
!
!---  Bottom ---
!
        IF( K.NE.1 ) THEN
          NB = N-IJFLD
          IF( IXP(NB).EQ.0 .OR. INBS(1,N).GT.0 ) GOTO 310
          DO M = 1,ISVF
            MB = MADJ(M)
            MP = MNOD(M)
            FLWB = XLW(MB,NB)*RHOL(MB,NB)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWB,FLWP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
            FGWB = XGW(MB,NB)*RHOG(MB,NB)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            INDX = 3
            FGW = DIFMN( FGWB,FGWP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
            WLW(M,NPZ) = WL(M,NPZ)*FLW - WTMW*WDLA(M,NPZ)
     &        - WTMW*WDS(M,NPZ)/WTMS
            WGW(M,NPZ) = WG(M,NPZ)*FGW + WTMW*WDGW(M,NPZ)
            FW(M,1) = -AFZ(NPZ)*(WL(M,NPZ)*FLW + WG(M,NPZ)*FGW
     &        + WTMW*WDGW(M,NPZ) - WTMW*WDLA(M,NPZ)
     &        - WTMW*WDS(M,NPZ)/WTMS)
          ENDDO
        ENDIF
  310   CONTINUE
!
!---  South ---
!
        IF( J.NE.1 ) THEN
          NS = N-IFLD
          IF( IXP(NS).EQ.0 .OR. INBS(2,N).GT.0 ) GOTO 410
          DO M = 1,ISVF
            MS = MADJ(M)
            MP = MNOD(M)
            FLWS = XLW(MS,NS)*RHOL(MS,NS)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWS,FLWP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
            FGWS = XGW(MS,NS)*RHOG(MS,NS)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            INDX = 3
            FGW = DIFMN( FGWS,FGWP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
            VLW(M,NPY) = VL(M,NPY)*FLW - WTMW*VDLA(M,NPY)
     &        - WTMW*VDS(M,NPY)/WTMS
            VGW(M,NPY) = VG(M,NPY)*FGW + WTMW*VDGW(M,NPY)
            FW(M,2) = -AFY(NPY)*(VL(M,NPY)*FLW + VG(M,NPY)*FGW
     &        + WTMW*VDGW(M,NPY) - WTMW*VDLA(M,NPY)
     &        - WTMW*VDS(M,NPY)/WTMS)
          ENDDO
        ENDIF
  410   CONTINUE
!
!---  West ---
!
        IF( I.NE.1 ) THEN
          NW = N-1
          IF( IXP(NW).EQ.0 .OR. INBS(3,N).GT.0 ) GOTO 510
          DO M = 1,ISVF
            MW = MADJ(M)
            MP = MNOD(M)
            FLWW = XLW(MW,NW)*RHOL(MW,NW)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWW,FLWP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            FGWW = XGW(MW,NW)*RHOG(MW,NW)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            INDX = 3
            FGW = DIFMN( FGWW,FGWP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            ULW(M,NPX) = UL(M,NPX)*FLW + WTMW*UDGW(M,NPX)
     &        - WTMW*UDS(M,NPX)/WTMS
            UGW(M,NPX) = UG(M,NPX)*FGW + WTMW*UDGW(M,NPX)
            FW(M,3) = -AFX(NPX)*(UL(M,NPX)*FLW + UG(M,NPX)*FGW
     &        + WTMW*UDGW(M,NPX) - WTMW*UDLA(M,NPX)
     &        - WTMW*UDS(M,NPX)/WTMS)
          ENDDO
        ENDIF
  510   CONTINUE
!
!---  East ---
!
        IF( I.NE.IFLD ) THEN
          NE = N+1
          IF( IXP(NE).EQ.0 .OR. INBS(4,N).GT.0 ) GOTO 610
          DO M = 1,ISVF
            ME = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FLWE = XLW(ME,NE)*RHOL(ME,NE)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWP,FLWE,DXGF(N),DXGF(NE),UL(1,NQX),INDX )
            FGWE = XGW(ME,NE)*RHOG(ME,NE)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            INDX = 3
            FGW = DIFMN( FGWP,FGWE,DXGF(N),DXGF(NE),UG(1,NQX),INDX )
            FW(M,4) = AFX(NQX)*(UL(MF,NQX)*FLW + UG(MF,NQX)*FGW
     &        + WTMW*UDGW(MF,NQX) - WTMW*UDLA(MF,NQX)
     &        - WTMW*UDS(MF,NQX)/WTMS)
          ENDDO
        ENDIF
  610   CONTINUE
!
!---  North ---
!
        IF( J.NE.JFLD ) THEN
          NN = N+IFLD
          IF( IXP(NN).EQ.0 .OR. INBS(5,N).GT.0 ) GOTO 710
          DO M = 1,ISVF
            MN = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FLWN = XLW(MN,NN)*RHOL(MN,NN)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWP,FLWN,DYGF(N),DYGF(NN),VL(1,NQY),INDX )
            FGWN = XGW(MN,NN)*RHOG(MN,NN)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            INDX = 3
            FGW = DIFMN( FGWP,FGWN,DYGF(N),DYGF(NN),VG(1,NQY),INDX )
            FW(M,5) = AFY(NQY)*(VL(MF,NQY)*FLW + VG(MF,NQY)*FGW
     &        + WTMW*VDGW(MF,NQY) - WTMW*VDLA(MF,NQY)
     &        - WTMW*VDS(MF,NQY)/WTMS)
          ENDDO
        ENDIF
  710   CONTINUE
!
!---  Top ---
!
        IF( K.NE.KFLD ) THEN
          NT = N+IJFLD
          IF( IXP(NT).EQ.0 .OR. INBS(6,N).GT.0 ) GOTO 810
          DO M = 1,ISVF
            MT = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FLWT = XLW(MT,NT)*RHOL(MT,NT)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            INDX = 2
            FLW = DIFMN( FLWP,FLWT,DZGF(N),DZGF(NT),WL(1,NQZ),INDX )
            FGWT = XGW(MT,NT)*RHOG(MT,NT)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            INDX = 3
            FGW = DIFMN( FGWP,FGWT,DZGF(N),DZGF(NT),WG(1,NQZ),INDX )
            FW(M,6) = AFZ(NQZ)*(WL(MF,NQZ)*FLW + WG(MF,NQZ)*FGW
     &        + WTMW*WDGW(MF,NQZ) - WTMW*WDLA(MF,NQZ)
     &        - WTMW*WDS(MF,NQZ)/WTMS)
          ENDDO
        ENDIF
  810   CONTINUE
!
!---  Compute water equation residuals  ---
!
        RWS = STWX(1) - SRCW(2,N)
        DO MD = 1,6
          RWS = RWS + FW(1,MD)
        ENDDO
        DO M = 1,ISVC
          RWP(M) = STWX(M+1) - SRCW(M+2,N)
          MM = 2*M
          DO MD = 1,6
            RWP(M) = RWP(M) + FW(MM,MD)
          ENDDO
        ENDDO
        DO M = 1,ISVC
          MM = 2*M + 1
          DO MD = 1,6
            RWA(M,MD) = RWS - FW(1,MD) + FW(MM,MD)
          ENDDO
        ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        RWS = 0.D+0
        DO M = 1,ISVC
          RWP(M) = 0.D+0
          DO MD = 1,6
            RWA(M,MD) = 0.D+0
          ENDDO
        ENDDO
        RWP(IEQWX) = 1.D+0
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
        CALL JCBL_GT( RWS,RWP,RWA,N,I,J,K,IEQWX )
!
!---  Continue to Next Node  ---
!
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBW_GT group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGAB_GT( N,NB,NPZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, air equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGAB_GT'
      DO M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLAB = XGAB(MB,NB)*RHOGB(MB,NB)
        FLAP = XGA(MB,N)*RHOG(MB,N)
        INDX = 3
        FLA = DIFMN( FLAB,FLAP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
        RS(M) = -AFZ(NPZ)*(WG(MP,NPZ)*FLA - WTMA*WDGW(MP,NPZ))
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQA )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGAB_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGAS_GT( N,NB,NPY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, air equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGAS_GT'
      DO M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLAB = XGAB(MB,NB)*RHOGB(MB,NB)
        FLAP = XGA(MB,N)*RHOG(MB,N)
        INDX = 3
        FLA = DIFMN( FLAB,FLAP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
        RS(M) = -AFY(NPY)*(VG(MP,NPY)*FLA - WTMA*VDGW(MP,NPY))
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQA )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGAS_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGAW_GT( N,NB,NPX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, air equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGAW_GT'
      DO M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLAB = XGAB(MB,NB)*RHOGB(MB,NB)
        FLAP = XGA(MB,N)*RHOG(MB,N)
        INDX = 3
        FLA = DIFMN( FLAB,FLAP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
        RS(M) = -AFX(NPX)*(UG(MP,NPX)*FLA - WTMA*UDGW(MP,NPX))
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQA )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGAW_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGAE_GT( N,NB,NQX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, air equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGAE_GT'
      DO M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLAB = XGAB(MB,NB)*RHOGB(MB,NB)
        FLAP = XGA(MB,N)*RHOG(MB,N)
        INDX = 3
        FLA = DIFMN( FLAP,FLAB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
        RS(M) = AFX(NQX)*(UG(MN,NQX)*FLA - WTMA*UDGW(MN,NQX))
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQA )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGAE_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGAN_GT( N,NB,NQY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, air equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGAN_GT'
      DO M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLAB = XGAB(MB,NB)*RHOGB(MB,NB)
        FLAP = XGA(MB,N)*RHOG(MB,N)
        INDX = 3
        FLA = DIFMN( FLAP,FLAB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
        RS(M) = AFY(NQY)*(VG(MN,NQY)*FLA - WTMA*VDGW(MN,NQY))
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQA )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGAN_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGAT_GT( N,NB,NQZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, air equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGAT_GT'
      DO M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLAB = XGAB(MB,NB)*RHOGB(MB,NB)
        FLAP = XGA(MB,N)*RHOG(MB,N)
        INDX = 3
        FLA = DIFMN( FLAP,FLAB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
        RS(M) = AFZ(NQZ)*(WG(MN,NQZ)*FLA - WTMA*WDGW(MN,NQZ))
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQA )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGAT_GT group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGWB_GT( N,NB,NPZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, water equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGWB_GT'
      DO M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLWB = XGWB(MB,NB)*RHOGB(MB,NB)
        FLWP = XGW(MB,N)*RHOG(MB,N)
        INDX = 3
        FLW = DIFMN( FLWB,FLWP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
        WGW(MP,NPZ) = WG(MP,NPZ)*FLW + WTMW*WDGW(MP,NPZ)
        RS(M) = -AFZ(NPZ)*(WG(MP,NPZ)*FLW + WTMW*WDGW(MP,NPZ))
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQW )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGWB_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGWS_GT( N,NB,NPY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, water equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGWS_GT'
      DO M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLWB = XGWB(MB,NB)*RHOGB(MB,NB)
        FLWP = XGW(MB,N)*RHOG(MB,N)
        INDX = 3
        FLW = DIFMN( FLWB,FLWP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
        VGW(MP,NPY) = VG(MP,NPY)*FLW + WTMW*VDGW(MP,NPY)
        RS(M) = -AFY(NPY)*(VG(MP,NPY)*FLW + WTMW*VDGW(MP,NPY))
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQW )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGWS_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGWW_GT( N,NB,NPX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, water equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGWW_GT'
      DO M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLWB = XGWB(MB,NB)*RHOGB(MB,NB)
        FLWP = XGW(MB,N)*RHOG(MB,N)
        INDX = 3
        FLW = DIFMN( FLWB,FLWP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
        UGW(MP,NPX) = UG(MP,NPX)*FLW + WTMW*UDGW(MP,NPX)
        RS(M) = -AFX(NPX)*(UG(MP,NPX)*FLW + WTMW*UDGW(MP,NPX))
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQW )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGWW_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGWE_GT( N,NB,NQX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, water equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGWE_GT'
      DO M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLWB = XGWB(MB,NB)*RHOGB(MB,NB)
        FLWP = XGW(MB,N)*RHOG(MB,N)
        INDX = 3
        FLW = DIFMN( FLWP,FLWB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
        UGW(MN,NQX) = UG(MN,NQX)*FLW + WTMW*UDGW(MN,NQX)
        RS(M) = AFX(NQX)*(UG(MN,NQX)*FLW + WTMW*UDGW(MN,NQX))
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQW )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGWE_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGWN_GT( N,NB,NQY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, water equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGWN_GT'
      DO M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLWB = XGWB(MB,NB)*RHOGB(MB,NB)
        FLWP = XGW(MB,N)*RHOG(MB,N)
        INDX = 3
        FLW = DIFMN( FLWP,FLWB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
        VGW(MN,NQY) = VG(MN,NQY)*FLW + WTMW*VDGW(MN,NQY)
        RS(M) = AFY(NQY)*(VG(MN,NQY)*FLW + WTMW*VDGW(MN,NQY))
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQW )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGWN_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBGWT_GT( N,NB,NQZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (gas boundary, water equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBGWT_GT'
      DO M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLWB = XGWB(MB,NB)*RHOGB(MB,NB)
        FLWP = XGW(MB,N)*RHOG(MB,N)
        INDX = 3
        FLW = DIFMN( FLWP,FLWB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
        WGW(MN,NQZ) = WG(MN,NQZ)*FLW + WTMW*WDGW(MN,NQZ)
        RS(M) = AFZ(NQZ)*(WG(MN,NQZ)*FLW + WTMW*WDGW(MN,NQZ))
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQW )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBGWT_GT group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLAB_GT( N,NB,NPZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, air equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLAB_GT'
      DO M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLAB = XLAB(MB,NB)*RHOLB(MB,NB)
        FLAP = XLA(MB,N)*RHOL(MB,N)
        INDX = 2
        FLA = DIFMN( FLAB,FLAP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
        RS(M) = -AFZ(NPZ)*(WL(MP,NPZ)*FLA + WDLA(MP,NPZ)*WTMA)
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQA )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLAB_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLAS_GT( N,NB,NPY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, air equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLAS_GT'
      DO M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLAB = XLAB(MB,NB)*RHOLB(MB,NB)
        FLAP = XLA(MB,N)*RHOL(MB,N)
        INDX = 2
        FLA = DIFMN( FLAB,FLAP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
        RS(M) = -AFY(NPY)*(VL(MP,NPY)*FLA + VDLA(MP,NPY)*WTMA)
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQA )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLAS_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLAW_GT( N,NB,NPX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, air equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLAW_GT'
      DO M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLAB = XLAB(MB,NB)*RHOLB(MB,NB)
        FLAP = XLA(MB,N)*RHOL(MB,N)
        INDX = 2
        FLA = DIFMN( FLAB,FLAP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
        RS(M) = -AFX(NPX)*(UL(MP,NPX)*FLA + UDLA(MP,NPX)*WTMA)
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQA )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLAW_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLAE_GT( N,NB,NQX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, air equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLAE_GT'
      DO M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLAB = XLAB(MB,NB)*RHOLB(MB,NB)
        FLAP = XLA(MB,N)*RHOL(MB,N)
        INDX = 2
        FLA = DIFMN( FLAP,FLAB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
        RS(M) = AFX(NQX)*(UL(MN,NQX)*FLA + UDLA(MN,NQX)*WTMA)
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQA )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLAE_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLAN_GT( N,NB,NQY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, air equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLAN_GT'
      DO M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLAB = XLAB(MB,NB)*RHOLB(MB,NB)
        FLAP = XLA(MB,N)*RHOL(MB,N)
        INDX = 2
        FLA = DIFMN( FLAP,FLAB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
        RS(M) = AFY(NQY)*(VL(MN,NQY)*FLA + VDLA(MN,NQY)*WTMA)
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQA )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLAN_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLAT_GT( N,NB,NQZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, air equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLAT_GT'
      DO M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLAB = XLAB(MB,NB)*RHOLB(MB,NB)
        FLAP = XLA(MB,N)*RHOL(MB,N)
        INDX = 2
        FLA = DIFMN( FLAP,FLAB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
        RS(M) = AFZ(NQZ)*(WL(MN,NQZ)*FLA + WDLA(MN,NQZ)*WTMA )
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQA )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLAT_GT group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWB_GT( N,NB,NPZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, water equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLWB_GT'
      DO M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLWB = XLWB(MB,NB)*RHOLB(MB,NB)
        FLWP = XLW(MB,N)*RHOL(MB,N)
        INDX = 2
        FLW = DIFMN( FLWB,FLWP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
        WLW(MP,NPZ) = WL(MP,NPZ)*FLW - WTMW*WDLA(MP,NPZ)
     &    - WTMW*WDS(MP,NPZ)/WTMS
        RS(M) = -AFZ(NPZ)*(WL(MP,NPZ)*FLW - WTMW*WDLA(MP,NPZ)
     &    - WTMW*WDS(MP,NPZ)/WTMS)
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQW )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWB_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWS_GT( N,NB,NPY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, water equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLWS_GT'
      DO M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLWB = XLWB(MB,NB)*RHOLB(MB,NB)
        FLWP = XLW(MB,N)*RHOL(MB,N)
        INDX = 2
        FLW = DIFMN( FLWB,FLWP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
        VLW(MP,NPY) = VL(MP,NPY)*FLW - WTMW*VDLA(MP,NPY)
     &    - WTMW*VDS(MP,NPY)/WTMS
        RS(M) = -AFY(NPY)*(VL(MP,NPY)*FLW - WTMW*VDLA(MP,NPY)
     &    - WTMW*VDS(MP,NPY)/WTMS)
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQW )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWS_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWW_GT( N,NB,NPX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, water equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLWW_GT'
      DO M = 1,ISVC+1
        MP = MPOSB(M)
        MB = M + 1
        FLWB = XLWB(MB,NB)*RHOLB(MB,NB)
        FLWP = XLW(MB,N)*RHOL(MB,N)
        INDX = 2
        FLW = DIFMN( FLWB,FLWP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
        ULW(MP,NPX) = UL(MP,NPX)*FLW - WTMW*UDLA(MP,NPX)
     &    - WTMW*UDS(MP,NPX)/WTMS
        RS(M) = -AFX(NPX)*(UL(MP,NPX)*FLW - WTMW*UDLA(MP,NPX)
     &    - WTMW*UDS(MP,NPX)/WTMS)
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQW )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWW_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWE_GT( N,NB,NQX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, water equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLWE_GT'
      DO M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLWB = XLWB(MB,NB)*RHOLB(MB,NB)
        FLWP = XLW(MB,N)*RHOL(MB,N)
        INDX = 2
        FLW = DIFMN( FLWP,FLWB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
        ULW(MN,NQX) = UL(MN,NQX)*FLW - WTMW*UDLA(MN,NQX)
     &    - WTMW*UDS(MN,NQX)/WTMS
        RS(M) = AFX(NQX)*(UL(MN,NQX)*FLW - WTMW*UDLA(MN,NQX)
     &    - WTMW*UDS(MN,NQX)/WTMS)
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQW )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWE_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWN_GT( N,NB,NQY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, water equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLWN_GT'
      DO M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLWB = XLWB(MB,NB)*RHOLB(MB,NB)
        FLWP = XLW(MB,N)*RHOL(MB,N)
        INDX = 2
        FLW = DIFMN( FLWP,FLWB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
        VLW(MN,NQY) = VL(MN,NQY)*FLW - WTMW*VDLA(MN,NQY)
     &        - WTMW*VDS(MN,NQY)/WTMS
        RS(M) = AFY(NQY)*(VL(MN,NQY)*FLW - WTMW*VDLA(MN,NQY)
     &        - WTMW*VDS(MN,NQY)/WTMS)
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQW )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWN_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLWT_GT( N,NB,NQZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (aqueous boundary, water equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLWT_GT'
      DO M = 1,ISVC+1
        MN = MNEGB(M)
        MB = M + 1
        FLWB = XLWB(MB,NB)*RHOLB(MB,NB)
        FLWP = XLW(MB,N)*RHOL(MB,N)
        INDX = 2
        FLW = DIFMN( FLWP,FLWB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
        WLW(MN,NQZ) = WL(MN,NQZ)*FLW - WTMW*WDLA(MN,NQZ)
     &    - WTMW*WDS(MN,NQZ)/WTMS
        RS(M) = AFZ(NQZ)*(WL(MN,NQZ)*FLW - WTMW*WDLA(MN,NQZ)
     &    - WTMW*WDS(MN,NQZ)/WTMS)
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQW )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLWT_GT group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBSB_GT( N,NB,NPZ )
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
!     Modify the surfactant equations of the Jacobian matrix for
!     boundary conditions on the bottom surface in the aqueous phase.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBSB_GT'
      DO M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFZ(NPZ)*WS(MP,NPZ)
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQS )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBSB_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBSS_GT( N,NB,NPY )
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
!     Modify the surfactant equations of the Jacobian matrix for
!     boundary conditions on the south surface in the aqueous phase.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBSS_GT'
      DO M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFY(NPY)*VS(MP,NPY)
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQS )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBSS_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBSW_GT( N,NB,NPX )
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
!     Modify the surfactant equations of the Jacobian matrix for
!     boundary conditions on the west surface in the aqueous phase.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBSW_GT'
      DO M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFX(NPX)*US(MP,NPX)
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQS )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBSW_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBSE_GT( N,NB,NQX )
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
!     Modify the surfactant equations of the Jacobian matrix for
!     boundary conditions on the east surface in the aqueous phase.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBSE_GT'
      DO M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFX(NQX)*US(MN,NQX)
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQS )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBSE_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBSN_GT( N,NB,NQY )
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
!     Modify the surfactant equations of the Jacobian matrix for
!     boundary conditions on the north surface in the aqueous phase.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBSN_GT'
      DO M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFY(NQY)*VS(MN,NQY)
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQS )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBSN_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBST_GT( N,NB,NQZ )
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
!     Modify the surfactant equations of the Jacobian matrix for
!     boundary conditions on the top surface in the aqueous phase.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBST_GT'
      DO M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFZ(NQZ)*WS(MN,NQZ)
      ENDDO
!
!---  Hot-dry-rock option  ---
!
      IF( ISLC(38).EQ.1 ) THEN
        DO M = 1,ISVC
          RS(M) = 0.D+0
        ENDDO
      ENDIF
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQS )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBST_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTB_GT( N,NB,NPZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (thermal boundary, energy equation, bottom surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBTB_GT'
      DO M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFZ(NPZ)*WQ(MP,NPZ)
      ENDDO
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQT )
!
!---  End of JCBTB_GT group  ---
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTS_GT( N,NB,NPY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (thermal boundary, energy equation, south surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBTS_GT'
      DO M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFY(NPY)*VQ(MP,NPY)
      ENDDO
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTS_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTW_GT( N,NB,NPX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (thermal boundary, energy equation, west surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBTW_GT'
      DO M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFX(NPX)*UQ(MP,NPX)
      ENDDO
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTW_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTE_GT( N,NB,NQX )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (thermal boundary, energy equation, east surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBTE_GT'
      DO M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFX(NQX)*UQ(MN,NQX)
      ENDDO
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTE_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTN_GT( N,NB,NQY )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (thermal boundary, energy equation, north surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBTN_GT'
      DO M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFY(NQY)*VQ(MN,NQY)
      ENDDO
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTN_GT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTT_GT( N,NB,NQZ )
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
!     Modify the Jacobian matrix for boundary conditions.
!     (thermal boundary, energy equation, top surface)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 February 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBTT_GT'
      DO M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFZ(NQZ)*WQ(MN,NQZ)
      ENDDO
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_GT( RS,N,IEQT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTT_GT group  ---
!
      RETURN
      END

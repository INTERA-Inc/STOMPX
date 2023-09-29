!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBL_HYDT_KE( RSS,RSP,RSA,N,I,J,K,MEQ )
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
!     Written by MD White, PNNL, 17 June 2011.
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
      SUB_LOG(ISUB_LOG) = '/JCBL_HYDT_KE'
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
!---  End of JCBL_HYDT_KE group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLB_HYDT_KE( RS,N,MEQ )
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
!     Written by MD White, PNNL, 11 July 2011
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
      SUB_LOG(ISUB_LOG) = '/JCBLB_HYDT_KE'
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
!---  End of JCBLB_HYDT_KE group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBT_HYDT_KE
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
!     aqueous, gas, liquid CO2, and hydrate contributions
!     (zero flux boundary conditions).
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 September 2010
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
      USE FDVI
      USE FDVH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STT(LUK+1),RTP(LUK),RTA(LUK,6),FT(LSFV,6)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBT_HYDT_KE'
      IEQTX = IEQT
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
        NQX = NPX + 1
        NQY = NPY + IFLD
        NQZ = NPZ + IJFLD
        IZN = IZ(N)
!
!---  First-order, forward-difference, time differential  ---
!
        USTOX = RHOS(IZN)*CPS(IZN)*(1.D+0-PORD(1,N))*T(1,N)
        USLOX = PORD(1,N)*SL(1,N)*RHOL(1,N)*UEL(1,N)
        USGOX = PORD(1,N)*SG(1,N)*RHOG(1,N)*UEG(1,N)
        USNOX = PORD(1,N)*SN(1,N)*RHON(1,N)*UEN(1,N)
        USHOX = PORD(1,N)*SH(1,N)*RHOH(1,N)*HH(1,N)
        USIOX = PORD(1,N)*SI(1,N)*RHOI(1,N)*HI(1,N)
        USSOX = PORD(1,N)*SS(1,N)*RHOSP(1,N)*HSP(1,N)
        DO M = 1,ISVC+1
          MP = M + 1
          USTX = RHOS(IZN)*CPS(IZN)*(1.D+0-PORD(MP,N))*T(MP,N)
          USLX = PORD(MP,N)*SL(MP,N)*RHOL(MP,N)*UEL(MP,N)
          USGX = PORD(MP,N)*SG(MP,N)*RHOG(MP,N)*UEG(MP,N)
          USNX = PORD(MP,N)*SN(MP,N)*RHON(MP,N)*UEN(MP,N)
          USHX = PORD(MP,N)*SH(MP,N)*RHOH(MP,N)*HH(MP,N)
          USIX = PORD(MP,N)*SI(MP,N)*RHOI(MP,N)*HI(MP,N)
          USSX = PORD(MP,N)*SS(MP,N)*RHOSP(MP,N)*HSP(MP,N)
          DUSTX = (USTX-USTOX)
          DUSLX = (USLX-USLOX)
          DUSGX = (USGX-USGOX)
          DUSNX = (USNX-USNOX)
          DUSHX = (USHX-USHOX)
          DUSIX = (USIX-USIOX)
          DUSSX = (USSX-USSOX)
          STT(M) = DUSTX + DUSLX + DUSGX + DUSNX + DUSHX + DUSIX + DUSSX
          STT(M) = STT(M)*DTI*VOL(N)
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
!---  Compute air equation residuals  ---
!
        RTS = STT(1) - SRCT(2,N)
        DO MD = 1,6
          RTS = RTS + FT(1,MD)
        ENDDO
        DO M = 1,ISVC
          RTP(M) = STT(M+1) - SRCT(M+2,N)
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
!---    Load Jacobian Matrix  ---
!
        CALL JCBL_HYDT_KE( RTS,RTP,RTA,N,I,J,K,IEQTX )
!
!---  Continue to Next Node  ---
!
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBT_HYDT_KE group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTB( N,NB,NPZ )
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
!     Written by MD White, PNNL, 14 September 2010
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
      SUB_LOG(ISUB_LOG) = '/JCBTB'
      NBX = NB
      DO M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFZ(NPZ)*WQ(MP,NPZ)
      ENDDO
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_HYDT_KE( RS,N,IEQT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTS( N,NB,NPY )
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
!     Written by MD White, PNNL, 14 September 2010
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
      SUB_LOG(ISUB_LOG) = '/JCBTS'
      NBX = NB
      DO M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFY(NPY)*VQ(MP,NPY)
      ENDDO
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_HYDT_KE( RS,N,IEQT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTW( N,NB,NPX )
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
!     Written by MD White, PNNL, 14 September 2010
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
      SUB_LOG(ISUB_LOG) = '/JCBTW'
      NBX = NB
      DO M = 1,ISVC+1
        MP = MPOSB(M)
        RS(M) = -AFX(NPX)*UQ(MP,NPX)
      ENDDO
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_HYDT_KE( RS,N,IEQT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTE( N,NB,NQX )
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
!     Written by MD White, PNNL, 14 September 2010
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
      SUB_LOG(ISUB_LOG) = '/JCBTE'
      NBX = NB
      DO M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFX(NQX)*UQ(MN,NQX)
      ENDDO
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_HYDT_KE( RS,N,IEQT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTN( N,NB,NQY )
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
!     Written by MD White, PNNL, 14 September 2010
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
      SUB_LOG(ISUB_LOG) = '/JCBTN'
      NBX = NB
      DO M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFY(NQY)*VQ(MN,NQY)
      ENDDO
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_HYDT_KE( RS,N,IEQT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBTT( N,NB,NQZ )
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
!     Written by MD White, PNNL, 14 September 2010
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
      SUB_LOG(ISUB_LOG) = '/JCBTT'
      NBX = NB
      DO M = 1,ISVC+1
        MN = MNEGB(M)
        RS(M) = AFZ(NQZ)*WQ(MN,NQZ)
      ENDDO
!
!---  Load Jacobian Matrix  ---
!
      CALL JCBLB_HYDT_KE( RS,N,IEQT )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBTT group  ---
!
      RETURN
      END

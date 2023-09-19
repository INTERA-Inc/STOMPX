!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGW
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
!     Compute water vapor mole diffusion rates.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
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
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGW'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
!
!---  X-direction vapor mole diffusion, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO 200 K = 1,KFLD
        DO 200 J = 1,JFLD
        DO 200 I = 2,IFLD
          N = ND(I,J,K)
          NW = N-1
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
          NPX = NSX(N)
          DXMGW = XMGW(2,NW)-XMGW(2,N)
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &        DFGW(MP,N)*RHOMG(MP,N)
            DFW = TORG(MN,NW)*PORD(MN,NW)*(SG(MN,NW)-SGT(MN,NW))*
     &        DFGW(MN,NW)*RHOMG(MN,NW)
            INDX = 12
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXMGW,INDX )
            UDGW(M,NPX) = DFM*(XMGW(MN,NW)-XMGW(MP,N))/DXGP(NPX)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            FGWW = XGW(MN,NW)*RHOG(MN,NW)
            INDX = 3
            FGW = DIFMN( FGWW,FGWP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            UGW(M,NPX) = UG(M,NPX)*FGW + WTMW*UDGW(M,NPX)
  100     CONTINUE
  200   CONTINUE
      ENDIF
!
!---  Y-direction vapor mole diffusion, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO 400 I = 1,IFLD
        DO 400 K = 1,KFLD
        DO 400 J = 2,JFLD
          N = ND(I,J,K)
          NS = N-IFLD
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
          NPY = NSY(N)
          DXMGW = XMGW(2,NS)-XMGW(2,N)
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &        DFGW(MP,N)*RHOMG(MP,N)
            DFS = TORG(MN,NS)*PORD(MN,NS)*(SG(MN,NS)-SGT(MN,NS))*
     &        DFGW(MN,NS)*RHOMG(MN,NS)
            INDX = 12
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXMGW,INDX )
            VDGW(M,NPY) = DFM*(XMGW(MN,NS)-XMGW(MP,N))/DYGP(NPY)/RP(I)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            FGWS = XGW(MN,NS)*RHOG(MN,NS)
            INDX = 3
            FGW = DIFMN( FGWS,FGWP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
            VGW(M,NPY) = VG(M,NPY)*FGW + WTMW*VDGW(M,NPY)
  300     CONTINUE
  400   CONTINUE
      ENDIF
!
!---  Z-direction vapor mole diffusion, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO 600 J = 1,JFLD
        DO 600 I = 1,IFLD
        DO 600 K = 2,KFLD
          N = ND(I,J,K)
          NB = N-IJFLD
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
          NPZ = NSZ(N)
          DXMGW = XMGW(2,NB)-XMGW(2,N)
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &        DFGW(MP,N)*RHOMG(MP,N)
            DFB = TORG(MN,NB)*PORD(MN,NB)*(SG(MN,NB)-SGT(MN,NB))*
     &        DFGW(MN,NB)*RHOMG(MN,NB)
            INDX = 12
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXMGW,INDX )
            WDGW(M,NPZ) = DFM*(XMGW(MN,NB)-XMGW(MP,N))/DZGP(NPZ)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            FGWB = XGW(MN,NB)*RHOG(MN,NB)
            INDX = 3
            FGW = DIFMN( FGWB,FGWP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
            WGW(M,NPZ) = WG(M,NPZ)*FGW + WTMW*WDGW(M,NPZ)
  500     CONTINUE
  600   CONTINUE
      ENDIF
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
!
!---  X-direction vapor mole diffusion, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO 1200 K = 1,KFLD
        DO 1200 J = 1,JFLD
        DO 1200 I = 2,IFLD
          N = ND(I,J,K)
          NW = N-1
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 1200
          NPX = NSX(N)
          DXMGW = XMGW(2,NW)*RHOMG(2,NW)-XMGW(2,N)*RHOMG(2,N)
          DO 1100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &        DFGW(MP,N)
            DFW = TORG(MN,NW)*PORD(MN,NW)*(SG(MN,NW)-SGT(MN,NW))*
     &        DFGW(MN,NW)
            INDX = 12
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXMGW,INDX )
            UDGW(M,NPX) = DFM*(XMGW(MN,NW)*RHOMG(MN,NW)
     &        - XMGW(MP,N)*RHOMG(MP,N))/DXGP(NPX)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            FGWW = XGW(MN,NW)*RHOG(MN,NW)
            INDX = 3
            FGW = DIFMN( FGWW,FGWP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            UGW(M,NPX) = UG(M,NPX)*FGW + WTMW*UDGW(M,NPX)
 1100     CONTINUE
 1200   CONTINUE
      ENDIF
!
!---  Y-direction vapor mole diffusion, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO 1400 I = 1,IFLD
        DO 1400 K = 1,KFLD
        DO 1400 J = 2,JFLD
          N = ND(I,J,K)
          NS = N-IFLD
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 1400
          NPY = NSY(N)
          DXMGW = XMGW(2,NS)*RHOMG(2,NS)-XMGW(2,N)*RHOMG(2,N)
          DO 1300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &        DFGW(MP,N)
            DFS = TORG(MN,NS)*PORD(MN,NS)*(SG(MN,NS)-SGT(MN,NS))*
     &        DFGW(MN,NS)
            INDX = 12
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXMGW,INDX )
            VDGW(M,NPY) = DFM*(XMGW(MN,NS)*RHOMG(MN,NS)
     &        - XMGW(MP,N)*RHOMG(MP,N))/DYGP(NPY)/RP(I)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            FGWS = XGW(MN,NS)*RHOG(MN,NS)
            INDX = 3
            FGW = DIFMN( FGWS,FGWP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
            VGW(M,NPY) = VG(M,NPY)*FGW + WTMW*VDGW(M,NPY)
 1300     CONTINUE
 1400   CONTINUE
      ENDIF
!
!---  Z-direction vapor mole diffusion, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO 1600 J = 1,JFLD
        DO 1600 I = 1,IFLD
        DO 1600 K = 2,KFLD
          N = ND(I,J,K)
          NB = N-IJFLD
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 1600
          NPZ = NSZ(N)
          DXMGW = XMGW(2,NB)*RHOMG(2,NB)-XMGW(2,N)*RHOMG(2,N)
          DO 1500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &        DFGW(MP,N)
            DFB = TORG(MN,NB)*PORD(MN,NB)*(SG(MN,NB)-SGT(MN,NB))*
     &        DFGW(MN,NB)
            INDX = 12
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXMGW,INDX )
            WDGW(M,NPZ) = DFM*(XMGW(MN,NB)*RHOMG(MN,NB)
     &        - XMGW(MP,N)*RHOMG(MP,N))/DZGP(NPZ)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            FGWB = XGW(MN,NB)*RHOG(MN,NB)
            INDX = 3
            FGW = DIFMN( FGWB,FGWP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
            WGW(M,NPZ) = WG(M,NPZ)*FGW + WTMW*WDGW(M,NPZ)
 1500     CONTINUE
 1600   CONTINUE
      ENDIF
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGW group  ---
!
      RETURN
      END

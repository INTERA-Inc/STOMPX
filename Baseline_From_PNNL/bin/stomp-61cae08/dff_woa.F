!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGA
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
!     Compute molar diffusive air gas fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 May 2008.
!     Last modified by MD White, PNNL, 14 May 2008.
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
      SUB_LOG(ISUB_LOG) = '/DFFGA'
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
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NW).LT.0 ) MN = 2
            FGAP = XGA(MP,N)*RHOG(MP,N)
            FGAW = XGA(MN,NW)*RHOG(MN,NW)
            INDX = 3
            FGA = DIFMN( FGAW,FGAP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            UGA(M,NPX) = UG(M,NPX)*FGA - WTMA*(UDGW(M,NPX)+UDGO(M,NPX))
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
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NS).LT.0 ) MN = 2
            FGAP = XGA(MP,N)*RHOG(MP,N)
            FGAS = XGA(MN,NS)*RHOG(MN,NS)
            INDX = 3
            FGA = DIFMN( FGAS,FGAP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
            VGA(M,NPY) = VG(M,NPY)*FGA - WTMA*(VDGW(M,NPY)+VDGO(M,NPY))
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
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NB).LT.0 ) MN = 2
            FGAP = XGA(MP,N)*RHOG(MP,N)
            FGAB = XGA(MN,NB)*RHOG(MN,NB)
            INDX = 3
            FGA = DIFMN( FGAB,FGAP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
            WGA(M,NPZ) = WG(M,NPZ)*FGA - WTMA*(WDGW(M,NPZ)+WDGO(M,NPZ))
  500     CONTINUE
  600   CONTINUE
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGA group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGO
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
!     Compute molar diffusive oil gas fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 May 2008.
!     Last modified by MD White, PNNL, 14 May 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE FDVG
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
      SUB_LOG(ISUB_LOG) = '/DFFGO'
!
!---  X-direction vapor mole diffusion, excluding boundaries
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
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
          NPX = NSX(N)
          DXMGO = XMGO(2,NW) - XMGO(2,N)
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGO(MP,N)*RHOMG(MP,N)
            DFW = TORG(MN,NW)*PORD(MN,NW)*SG(MN,NW)*
     &        DFGO(MN,NW)*RHOMG(MN,NW)
            INDX = 12
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXMGO,INDX )
            UDGO(M,NPX) = DFM*(XMGO(MN,NW)-XMGO(MP,N))/DXGP(NPX)
            FGOP = XGO(MP,N)*RHOG(MP,N)
            FGOW = XGO(MN,NW)*RHOG(MN,NW)
            INDX = 3
            FGO = DIFMN( FGOW,FGOP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            UGO(M,NPX) = UG(M,NPX)*FGO + WTMO*UDGO(M,NPX)
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
          DXMGO = XMGO(2,NS) - XMGO(2,N)
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGO(MP,N)*RHOMG(MP,N)
            DFS = TORG(MN,NS)*PORD(MN,NS)*SG(MN,NS)*
     &        DFGO(MN,NS)*RHOMG(MN,NS)
            INDX = 12
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXMGO,INDX )
            VDGO(M,NPY) = DFM*(XMGO(MN,NS)-XMGO(MP,N))/
     &        (DYGP(NPY)*RP(I))
            FGOP = XGO(MP,N)*RHOG(MP,N)
            FGOS = XGO(MN,NS)*RHOG(MN,NS)
            INDX = 3
            FGO = DIFMN( FGOS,FGOP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
            VGO(M,NPY) = VG(M,NPY)*FGO + WTMO*VDGO(M,NPY)
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
          DXMGO = XMGO(2,NB) - XMGO(2,N)
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGO(MP,N)*RHOMG(MP,N)
            DFB = TORG(MN,NB)*PORD(MN,NB)*SG(MN,NB)*
     &        DFGO(MN,NB)*RHOMG(MN,NB)
            INDX = 12
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXMGO,INDX )
            WDGO(M,NPZ) = DFM*(XMGO(MN,NB)-XMGO(MP,N))/
     &        DZGP(NPZ)
            FGOP = XGO(MP,N)*RHOG(MP,N)
            FGOB = XGO(MN,NB)*RHOG(MN,NB)
            INDX = 3
            FGO = DIFMN( FGOB,FGOP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
            WGO(M,NPZ) = WG(M,NPZ)*FGO + WTMO*WDGO(M,NPZ)
  500     CONTINUE
  600   CONTINUE
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGO group  ---
!
      RETURN
      END

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
!     Compute molar diffusive water gas fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 May 2008.
!     Last modified by MD White, PNNL, 14 May 2008.
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
          DXMGW = XMGW(2,NW) - XMGW(2,N)
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGW(MP,N)*RHOMG(MP,N)
            DFW = TORG(MN,NW)*PORD(MN,NW)*SG(MN,NW)*
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
          DXMGW = XMGW(2,NS) - XMGW(2,N)
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGW(MP,N)*RHOMG(MP,N)
            DFS = TORG(MN,NS)*PORD(MN,NS)*SG(MN,NS)*
     &        DFGW(MN,NS)*RHOMG(MN,NS)
            INDX = 12
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXMGW,INDX )
            VDGW(M,NPY) = DFM*(XMGW(MN,NS)-XMGW(MP,N))/(DYGP(NPY)*RP(I))
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
          DXMGW = XMGW(2,NB) - XMGW(2,N)
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGW(MP,N)*RHOMG(MP,N)
            DFB = TORG(MN,NB)*PORD(MN,NB)*SG(MN,NB)*
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
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLA
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
!     Compute molar diffusive air aqueous fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 May 2008.
!     Last modified by MD White, PNNL, 14 May 2008.
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
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLA'
!
!---  X-direction molar diffusion, excluding boundaries
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
          DXLA = XMLA(2,NW) - XMLA(2,N)
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            IF( IXP(NW).LT.0 ) MN = 2
            IF( IXP(N).LT.0 ) MP = 2
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*
     &        DFLA(MP,N)*RHOML(MP,N)
            DFW = TORL(MN,NW)*PORD(MN,NW)*SL(MN,NW)*
     &        DFLA(MN,NW)*RHOML(MN,NW)
            INDX = 14
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXLA,INDX)
            UDLA(M,NPX) = DFM*(XMLA(MN,NW)-XMLA(MP,N))/DXGP(NPX)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            FLAW = XLA(MN,NW)*RHOL(MN,NW)
            INDX = 2
            FLA = DIFMN( FLAW,FLAP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            ULA(M,NPX) = UL(M,NPX)*FLA + WTMA*UDLA(M,NPX)
  100     CONTINUE
  200 CONTINUE
      ENDIF
!
!---  Y-direction molar diffusion, excluding boundaries
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
          DXLA = XMLA(2,NS) - XMLA(2,N)
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NS).LT.0 ) MN = 2
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*
     &        DFLA(MP,N)*RHOML(MP,N)
            DFS = TORL(MN,NS)*PORD(MN,NS)*SL(MN,NS)*
     &        DFLA(MN,NS)*RHOML(MN,NS)
            INDX = 14
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXLA,INDX )
            VDLA(M,NPY) = DFM*(XMLA(MN,NS)-XMLA(MP,N))/(DYGP(NPY)*RP(I))
            FLAP = XLA(MP,N)*RHOL(MP,N)
            FLAS = XLA(MN,NS)*RHOL(MN,NS)
            INDX = 2
            FLA = DIFMN( FLAS,FLAP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
            VLA(M,NPY) = VL(M,NPY)*FLA + WTMA*VDLA(M,NPY)
  300   CONTINUE
  400 CONTINUE
      ENDIF
!
!---  Z-direction molar diffusion, excluding boundaries
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
          DXLA = XMLA(2,NB) - XMLA(2,N)
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NB).LT.0 ) MN = 2
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*
     &        DFLA(MP,N)*RHOML(MP,N)
            DFB = TORL(MN,NB)*PORD(MN,NB)*SL(MN,NB)*
     &        DFLA(MN,NB)*RHOML(MN,NB)
            INDX = 14
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXLA,INDX)
            WDLA(M,NPZ) = DFM*(XMLA(MN,NB)-XMLA(MP,N))/DZGP(NPZ)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            FLAB = XLA(MN,NB)*RHOL(MN,NB)
            INDX = 2
            FLA = DIFMN( FLAB,FLAP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
            WLA(M,NPZ) = WL(M,NPZ)*FLA + WTMA*WDLA(M,NPZ)
  500   CONTINUE
  600 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLA group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLO
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
!     Compute molar diffusive oil aqueous fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 May 2008.
!     Last modified by MD White, PNNL, 14 May 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLO'
!
!---  X-direction molar diffusion, excluding boundaries
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
          DXLO = XMLO(2,NW) - XMLO(2,N)
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            IF( IXP(NW).LT.0 ) MN = 2
            IF( IXP(N).LT.0 ) MP = 2
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*
     &        DFLO(MP,N)*RHOML(MP,N)
            DFW = TORL(MN,NW)*PORD(MN,NW)*SL(MN,NW)*
     &        DFLO(MN,NW)*RHOML(MN,NW)
            INDX = 14
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXLO,INDX)
            UDLO(M,NPX) = DFM*(XMLO(MN,NW)-XMLO(MP,N))/DXGP(NPX)
            FLOP = XLO(MP,N)*RHOL(MP,N)
            FLOW = XLO(MN,NW)*RHOL(MN,NW)
            INDX = 2
            FLO = DIFMN( FLOW,FLOP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            ULO(M,NPX) = UL(M,NPX)*FLO + WTMO*UDLO(M,NPX)
  100     CONTINUE
  200 CONTINUE
      ENDIF
!
!---  Y-direction molar diffusion, excluding boundaries
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
          DXLO = XMLO(2,NS) - XMLO(2,N)
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NS).LT.0 ) MN = 2
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*
     &        DFLO(MP,N)*RHOML(MP,N)
            DFS = TORL(MN,NS)*PORD(MN,NS)*SL(MN,NS)*
     &        DFLO(MN,NS)*RHOML(MN,NS)
            INDX = 14
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXLO,INDX )
            VDLO(M,NPY) = DFM*(XMLO(MN,NS)-XMLO(MP,N))/(DYGP(NPY)*RP(I))
            FLOP = XLO(MP,N)*RHOL(MP,N)
            FLOS = XLO(MN,NS)*RHOL(MN,NS)
            INDX = 2
            FLO = DIFMN( FLOS,FLOP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
            VLO(M,NPY) = VL(M,NPY)*FLO + WTMO*VDLO(M,NPY)
  300   CONTINUE
  400 CONTINUE
      ENDIF
!
!---  Z-direction molar diffusion, excluding boundaries
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
          DXLO = XMLO(2,NB) - XMLO(2,N)
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NB).LT.0 ) MN = 2
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*
     &        DFLO(MP,N)*RHOML(MP,N)
            DFB = TORL(MN,NB)*PORD(MN,NB)*SL(MN,NB)*
     &        DFLO(MN,NB)*RHOML(MN,NB)
            INDX = 14
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXLO,INDX)
            WDLO(M,NPZ) = DFM*(XMLO(MN,NB)-XMLO(MP,N))/DZGP(NPZ)
            FLOP = XLO(MP,N)*RHOL(MP,N)
            FLOB = XLO(MN,NB)*RHOL(MN,NB)
            INDX = 2
            FLO = DIFMN( FLOB,FLOP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
            WLO(M,NPZ) = WL(M,NPZ)*FLO + WTMO*WDLO(M,NPZ)
  500   CONTINUE
  600 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLO group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLW
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
!     Compute molar diffusive water aqueous fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 May 2008.
!     Last modified by MD White, PNNL, 14 May 2008.
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
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLW'
!
!---  X-direction molar diffusion, excluding boundaries
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
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NW).LT.0 ) MN = 2
            FLWP = XLW(MP,N)*RHOL(MP,N)
            FLWW = XLW(MN,NW)*RHOL(MN,NW)
            INDX = 2
            FLW = DIFMN( FLWW,FLWP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            ULW(M,NPX) = UL(M,NPX)*FLW - WTMW*(UDLA(M,NPX) 
     &        + UDLO(M,NPX))
  100     CONTINUE
  200 CONTINUE
      ENDIF
!
!---  Y-direction molar diffusion, excluding boundaries
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
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NS).LT.0 ) MN = 2
            FLWP = XLW(MP,N)*RHOL(MP,N)
            FLWS = XLW(MN,NS)*RHOL(MN,NS)
            INDX = 2
            FLW = DIFMN( FLWS,FLWP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
            VLW(M,NPY) = VL(M,NPY)*FLW - WTMW*(VDLA(M,NPY) 
     &        + VDLO(M,NPY))
  300   CONTINUE
  400 CONTINUE
      ENDIF
!
!---  Z-direction molar diffusion, excluding boundaries
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
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NB).LT.0 ) MN = 2
            FLWP = XLW(MP,N)*RHOL(MP,N)
            FLWB = XLW(MN,NB)*RHOL(MN,NB)
            INDX = 2
            FLW = DIFMN( FLWB,FLWP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
            WLW(M,NPZ) = WL(M,NPZ)*FLW - WTMW*(WDLA(M,NPZ) 
     &        + WDLO(M,NPZ))
  500   CONTINUE
  600 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNO
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
!     Compute molar diffusive oil NAPL fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 May 2008.
!     Last modified by MD White, PNNL, 14 May 2008.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVA
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
      SUB_LOG(ISUB_LOG) = '/DFFNO'
!
!---  X-direction, excluding boundaries
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
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NW).LT.0 ) MN = 2
            FNOP = XNO(MP,N)*RHON(MP,N)
            FNOW = XNO(MN,NW)*RHON(MN,NW)
            INDX = 4
            FNO = DIFMN( FNOW,FNOP,DXGF(NW),DXGF(N),UN(1,NPX),INDX )
            UNA(M,NPX) = UN(M,NPX)*FNO
  100     CONTINUE
  200   CONTINUE
      ENDIF
!
!---  Y-direction, excluding boundaries
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
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NS).LT.0 ) MN = 2
            FNOP = XNO(MP,N)*RHON(MP,N)
            FNOS = XNO(MN,NS)*RHON(MN,NS)
            INDX = 4
            FNO = DIFMN( FNOS,FNOP,DYGF(NS),DYGF(N),VN(1,NPY),INDX )
            VNA(M,NPY) = VN(M,NPY)*FNO
  300     CONTINUE
  400   CONTINUE
      ENDIF
!
!---  Z-direction, excluding boundaries
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
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NB).LT.0 ) MN = 2
            FNOP = XNO(MP,N)*RHON(MP,N)
            FNOB = XNO(MN,NB)*RHON(MN,NB)
            INDX = 4
            FNO = DIFMN( FNOB,FNOP,DZGF(NB),DZGF(N),WN(1,NPZ),INDX )
            WNA(M,NPZ) = WN(M,NPZ)*FNO
  500     CONTINUE
  600   CONTINUE
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNO group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAB( N,NB )
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
!     Compute advective and molar diffusive air gas fluxes on 
!     a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFGAB'
      K = KD(N)
      NPZ = NSZ(N)
      DXMGA = XMGAB(2,NB) - XMGA(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FGAP = XGA(MP,N)*RHOG(MP,N)
        FGAB = XGAB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGA = DIFMN( FGAB,FGAP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( WG(1,NPZ).LT.-EPSL ) THEN
            WGA(M,NPZ) = WG(M,NPZ)*FGA
          ELSE
            WGA(M,NPZ) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.15 ) THEN
          IF( WG(1,NPZ).GT.EPSL ) THEN
            WGA(M,NPZ) = WG(M,NPZ)*FGA
          ELSE
            WGA(M,NPZ) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          WGA(M,NPZ) = WG(M,NPZ)*FGA - WTMA*(WDGO(M,NPZ)+WDGW(M,NPZ))
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAE( N,NB )
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
!     Compute advective and molar diffusive air gas fluxes on 
!     an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFGAE'
      I = ID(N)
      NQX = NSX(N)+1
      DXMGA = XMGA(2,N) - XMGAB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FGAP = XGA(MN,N)*RHOG(MN,N)
        FGAB = XGAB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGA = DIFMN( FGAP,FGAB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( UG(1,NQX).GT.EPSL ) THEN
            UGA(M,NQX) = UG(M,NQX)*FGA
          ELSE
            UGA(M,NQX) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.15 ) THEN
          IF( UG(1,NQX).LT.-EPSL ) THEN
            UGA(M,NQX) = UG(M,NQX)*FGA
          ELSE
            UGA(M,NQX) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          UGA(M,NQX) = UG(M,NQX)*FGA - WTMA*(UDGO(M,NQX)+UDGW(M,NQX))
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAN( N,NB )
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
!     Compute advective and molar diffusive air gas fluxes on 
!     a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFGAN'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMGA = XMGA(2,N) - XMGAB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FGAP = XGA(MN,N)*RHOG(MN,N)
        FGAB = XGAB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGA = DIFMN( FGAP,FGAB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( VG(1,NQY).GT.EPSL ) THEN
            VGA(M,NQY) = VG(M,NQY)*FGA
          ELSE
            VGA(M,NQY) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.15 ) THEN
          IF( VG(1,NQY).LT.-EPSL ) THEN
            VGA(M,NQY) = VG(M,NQY)*FGA
          ELSE
            VGA(M,NQY) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          VGA(M,NQY) = VG(M,NQY)*FGA - WTMA*(VDGO(M,NQY)+VDGW(M,NQY))
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAS( N,NB )
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
!     Compute advective and molar diffusive air gas fluxes on 
!     south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFGAS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMGA = XMGAB(2,NB) - XMGA(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FGAP = XGA(MP,N)*RHOG(MP,N)
        FGAB = XGAB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGA = DIFMN( FGAB,FGAP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( VG(1,NPY).LT.-EPSL ) THEN
            VGA(M,NPY) = VG(M,NPY)*FGA
          ELSE
            VGA(M,NPY) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.15 ) THEN
          IF( VG(1,NPY).GT.EPSL ) THEN
            VGA(M,NPY) = VG(M,NPY)*FGA
          ELSE
            VGA(M,NPY) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          VGA(M,NPY) = VG(M,NPY)*FGA - WTMA*(VDGO(M,NPY)-VDGW(M,NPY))
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAT( N,NB )
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
!     Compute advective and molar diffusive air gas fluxes on 
!     a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFGAT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMGA = XMGA(2,N)-XMGAB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FGAP = XGA(MN,N)*RHOG(MN,N)
        FGAB = XGAB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGA = DIFMN( FGAP,FGAB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( WG(1,NQZ).GT.EPSL ) THEN
            WGA(M,NQZ) = WG(M,NQZ)*FGA
          ELSE
            WGA(M,NQZ) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.15 ) THEN
          IF( WG(1,NQZ).LT.-EPSL ) THEN
            WGA(M,NQZ) = WG(M,NQZ)*FGA
          ELSE
            WGA(M,NQZ) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          WGA(M,NQZ) = WG(M,NQZ)*FGA - WTMA*(WDGO(M,NQZ)-WDGW(M,NQZ))
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAW( N,NB )
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
!     Compute advective and molar diffusive air gas fluxes on 
!     a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFGAW'
      I = ID(N)
      NPX = NSX(N)
      DXMGA = XMGAB(2,NB)-XMGA(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FGAP = XGA(MP,N)*RHOG(MP,N)
        FGAB = XGAB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGA = DIFMN( FGAB,FGAP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( UG(1,NPX).LT.-EPSL ) THEN
            UGA(M,NPX) = UG(M,NPX)*FGA
          ELSE
            UGA(M,NPX) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.15 ) THEN
          IF( UG(1,NPX).GT.EPSL ) THEN
            UGA(M,NPX) = UG(M,NPX)*FGA
          ELSE
            UGA(M,NPX) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          UGA(M,NPX) = UG(M,NPX)*FGA - WTMA*(UDGO(M,NPX)+UDGW(M,NPX))
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGOB( N,NB )
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
!     Compute advective and molar diffusive air gas fluxes on 
!     a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE CONST
      USE BCVP
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFGOB'
      K = KD(N)
      NPZ = NSZ(N)
      DXMGO = XMGOB(2,NB) - XMGO(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FGOP = XGO(MP,N)*RHOG(MP,N)
        FGOB = XGOB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGO = DIFMN( FGOB,FGOP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( WG(1,NPZ).LT.-EPSL ) THEN
            WGO(M,NPZ) = WG(M,NPZ)*FGO
          ELSE
            WGO(M,NPZ) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.15 ) THEN
          IF( WG(1,NPZ).GT.EPSL ) THEN
            WGO(M,NPZ) = WG(M,NPZ)*FGO
          ELSE
            WGO(M,NPZ) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          WGO(M,NPZ) = WG(M,NPZ)*FGO + WTMO*WDGO(M,NPZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGOB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGOE( N,NB )
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
!     Compute advective and molar diffusive air gas fluxes on 
!     an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE CONST
      USE BCVP
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFGOE'
      I = ID(N)
      NQX = NSX(N)+1
      DXMGO = XMGO(2,N) - XMGOB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FGOP = XGO(MN,N)*RHOG(MN,N)
        FGOB = XGOB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGO = DIFMN( FGOP,FGOB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( UG(1,NQX).GT.EPSL ) THEN
            UGO(M,NQX) = UG(M,NQX)*FGO
          ELSE
            UGO(M,NQX) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.15 ) THEN
          IF( UG(1,NQX).LT.-EPSL ) THEN
            UGO(M,NQX) = UG(M,NQX)*FGO
          ELSE
            UGO(M,NQX) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          UGO(M,NQX) = UG(M,NQX)*FGO + WTMO*UDGO(M,NQX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGAE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGON( N,NB )
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
!     Compute advective and molar diffusive air gas fluxes on 
!     a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE CONST
      USE BCVP
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFGON'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMGO = XMGO(2,N) - XMGOB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FGOP = XGO(MN,N)*RHOG(MN,N)
        FGOB = XGOB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGO = DIFMN( FGOP,FGOB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( VG(1,NQY).GT.EPSL ) THEN
            VGO(M,NQY) = VG(M,NQY)*FGO
          ELSE
            VGO(M,NQY) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.15 ) THEN
          IF( VG(1,NQY).LT.-EPSL ) THEN
            VGO(M,NQY) = VG(M,NQY)*FGO
          ELSE
            VGO(M,NQY) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          VGO(M,NQY) = VG(M,NQY)*FGO + WTMO*VDGO(M,NQY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGON group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGOS( N,NB )
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
!     Compute advective and molar diffusive air gas fluxes on 
!     south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE CONST
      USE BCVP
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFGOS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMGO = XMGOB(2,NB) - XMGO(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FGOP = XGO(MP,N)*RHOG(MP,N)
        FGOB = XGOB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGO = DIFMN( FGOB,FGOP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( VG(1,NPY).LT.-EPSL ) THEN
            VGO(M,NPY) = VG(M,NPY)*FGO
          ELSE
            VGO(M,NPY) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.15 ) THEN
          IF( VG(1,NPY).GT.EPSL ) THEN
            VGO(M,NPY) = VG(M,NPY)*FGO
          ELSE
            VGO(M,NPY) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          VGO(M,NPY) = VG(M,NPY)*FGO + WTMO*VDGO(M,NPY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGOS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGOT( N,NB )
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
!     Compute advective and molar diffusive air gas fluxes on 
!     a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE CONST
      USE BCVP
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFGOT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMGO = XMGO(2,N)-XMGOB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FGOP = XGO(MN,N)*RHOG(MN,N)
        FGOB = XGOB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGO = DIFMN( FGOP,FGOB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( WG(1,NQZ).GT.EPSL ) THEN
            WGO(M,NQZ) = WG(M,NQZ)*FGO
          ELSE
            WGO(M,NQZ) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.15 ) THEN
          IF( WG(1,NQZ).LT.-EPSL ) THEN
            WGO(M,NQZ) = WG(M,NQZ)*FGO
          ELSE
            WGO(M,NQZ) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          WGO(M,NQZ) = WG(M,NQZ)*FGO + WTMO*WDGO(M,NQZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGOT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGOW( N,NB )
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
!     Compute advective and molar diffusive air gas fluxes on 
!     a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE CONST
      USE BCVP
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFGOW'
      I = ID(N)
      NPX = NSX(N)
      DXMGO = XMGOB(2,NB)-XMGO(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FGOP = XGO(MP,N)*RHOG(MP,N)
        FGOB = XGOB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGO = DIFMN( FGOB,FGOP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( UG(1,NPX).LT.-EPSL ) THEN
            UGO(M,NPX) = UG(M,NPX)*FGO
          ELSE
            UGO(M,NPX) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.15 ) THEN
          IF( UG(1,NPX).GT.EPSL ) THEN
            UGO(M,NPX) = UG(M,NPX)*FGO
          ELSE
            UGO(M,NPX) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          UGO(M,NPX) = UG(M,NPX)*FGO + WTMO*UDGO(M,NPX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGOW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWB( N,NB )
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
!     Compute advective and molar diffusive oil gas fluxes on
!     a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCVP
      USE BCVG
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFGWB'
      K = KD(N)
      NPZ = NSZ(N)
      DXMGW = XMGWB(2,NB) - XMGW(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FGWP = XGW(MP,N)*RHOG(MP,N)
        FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGW = DIFMN( FGWB,FGWP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( WG(1,NPZ).LT.-EPSL ) THEN
            WGW(M,NPZ) = WG(M,NPZ)*FGW
          ELSE
            WGW(M,NPZ) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.15 ) THEN
          IF( WG(1,NPZ).GT.EPSL ) THEN
            WGW(M,NPZ) = WG(M,NPZ)*FGW
          ELSE
            WGW(M,NPZ) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))
     &      *DFGW(MP,N)*RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)
     &      *DFGWB(MP,NB)*RHOMGB(MP,NB)
          INDX = 13
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMGW,INDX )
          WDGW(M,NPZ) = DFM*(XMGWB(MP,NB)
     &      - XMGW(MP,N))/(5.D-1*DZGF(N))
          WGW(M,NPZ) = WG(M,NPZ)*FGW + WTMW*WDGW(M,NPZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWE( N,NB )
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
!     Compute advective and molar diffusive oil gas fluxes on
!     an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
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
      USE BCVP
      USE BCVG
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFGWE'
      I = ID(N)
      NQX = NSX(N)+1
      DXMGW = XMGW(2,N) - XMGWB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FGWP = XGW(MN,N)*RHOG(MN,N)
        FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGW = DIFMN( FGWP,FGWB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( UG(1,NQX).GT.EPSL ) THEN
            UGW(M,NQX) = UG(M,NQX)*FGW
          ELSE
            UGW(M,NQX) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.15 ) THEN
          IF( UG(1,NQX).LT.-EPSL ) THEN
            UGW(M,NQX) = UG(M,NQX)*FGW
          ELSE
            UGW(M,NQX) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))
     &      *DFGW(MN,N)*RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)
     &      *DFGWB(MN,NB)*RHOMGB(MN,NB)
          INDX = 13
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMGW,INDX )
          UDGW(M,NQX) = DFM*(XMGW(MN,N)
     &      - XMGWB(MN,NB))/(5.D-1*DXGF(N))
          UGW(M,NQX) = UG(M,NQX)*FGW + WTMW*UDGW(M,NQX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWN( N,NB )
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
!     Compute advective and molar diffusive oil gas fluxes on
!     a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
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
      USE BCVP
      USE BCVG
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFGWN'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMGW = XMGW(2,N) - XMGWB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FGWP = XGW(MN,N)*RHOG(MN,N)
        FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGW = DIFMN( FGWP,FGWB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( VG(1,NQY).GT.EPSL ) THEN
            VGW(M,NQY) = VG(M,NQY)*FGW
          ELSE
            VGW(M,NQY) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.15 ) THEN
          IF( VG(1,NQY).LT.-EPSL ) THEN
            VGW(M,NQY) = VG(M,NQY)*FGW
          ELSE
            VGW(M,NQY) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))
     &      *DFGW(MN,N)*RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)
     &      *DFGWB(MN,NB)*RHOMGB(MN,NB)
          INDX = 13
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMGW,INDX )
          VDGW(M,NQY) = DFM*(XMGW(MN,N)
     &      - XMGWB(MN,NB))/(5.D-1*DYGF(N)*RP(I))
          VGW(M,NQY) = VG(M,NQY)*FGW + WTMW*VDGW(M,NQY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWS( N,NB )
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
!     Compute advective and molar diffusive oil gas fluxes on
!     a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
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
      USE BCVP
      USE BCVG
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFGWS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMGW = XMGWB(2,NB) - XMGW(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FGWP = XGW(MP,N)*RHOG(MP,N)
        FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGW = DIFMN( FGWB,FGWP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( VG(1,NPY).LT.-EPSL ) THEN
            VGW(M,NPY) = VG(M,NPY)*FGW
          ELSE
            VGW(M,NPY) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.15 ) THEN
          IF( VG(1,NPY).GT.EPSL ) THEN
            VGW(M,NPY) = VG(M,NPY)*FGW
          ELSE
            VGW(M,NPY) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))
     &      *DFGW(MP,N)*RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)
     &      *DFGWB(MP,NB)*RHOMGB(MP,NB)
          INDX = 13
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMGW,INDX )
          VDGW(M,NPY) = DFM*(XMGWB(MP,NB)
     &      - XMGW(MP,N))/(5.D-1*DYGF(N)*RP(I))
          VGW(M,NPY) = VG(M,NPY)*FGW + WTMW*VDGW(M,NPY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWT( N,NB )
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
!     Compute advective and molar diffusive oil gas fluxes on
!     a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
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
      USE BCVP
      USE BCVG
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFGWT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMGW = XMGW(2,N)-XMGWB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FGWP = XGW(MN,N)*RHOG(MN,N)
        FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGW = DIFMN( FGWP,FGWB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( WG(1,NQZ).GT.EPSL ) THEN
            WGW(M,NQZ) = WG(M,NQZ)*FGW
          ELSE
            WGW(M,NQZ) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.15 ) THEN
          IF( WG(1,NQZ).LT.-EPSL ) THEN
            WGW(M,NQZ) = WG(M,NQZ)*FGW
          ELSE
            WGW(M,NQZ) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*(SG(MN,N)-SGT(MN,N))
     &      *DFGW(MN,N)*RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)
     &      *DFGWB(MN,NB)*RHOMGB(MN,NB)
          INDX = 13
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMGW,INDX )
          WDGW(M,NQZ) = DFM*(XMGW(MN,N)
     &     - XMGWB(MN,NB))/(5.D-1*DZGF(N))
          WGW(M,NQZ) = WG(M,NQZ)*FGW + WTMW*WDGW(M,NQZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGWW( N,NB )
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
!     Compute advective and molar diffusive oil gas fluxes on
!     a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
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
      USE BCVP
      USE BCVG
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFGWW'
      I = ID(N)
      NPX = NSX(N)
      DXMGW = XMGWB(2,NB)-XMGW(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FGWP = XGW(MP,N)*RHOG(MP,N)
        FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGW = DIFMN( FGWB,FGWP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( UG(1,NPX).LT.-EPSL ) THEN
            UGW(M,NPX) = UG(M,NPX)*FGW
          ELSE
            UGW(M,NPX) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.15 ) THEN
          IF( UG(1,NPX).GT.EPSL ) THEN
            UGW(M,NPX) = UG(M,NPX)*FGW
          ELSE
            UGW(M,NPX) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))
     &      *DFGW(MP,N)*RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)
     &      *DFGWB(MP,NB)*RHOMGB(MP,NB)
          INDX = 13
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMGW,INDX )
          UDGW(M,NPX) = DFM*(XMGWB(MP,NB)
     &      - XMGW(MP,N))/(5.D-1*DXGF(N))
          UGW(M,NPX) = UG(M,NPX)*FGW + WTMW*UDGW(M,NPX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAB( N,NB )
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
!     Compute advective and molar diffusive air aqueous fluxes on 
!     a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFLAB'
      K = KD(N)
      NPZ = NSZ(N)
      DXMLA = XMLAB(2,NB) - XMLA(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FLAP = XLA(MP,N)*RHOL(MP,N)
        FLAB = XLAB(MP,NB)*RHOLB(MP,NB)
        INDX = 3
        FLA = DIFMN( FLAB,FLAP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( WL(1,NPZ).LT.-EPSL ) THEN
            WLA(M,NPZ) = WL(M,NPZ)*FLA
          ELSE
            WLA(M,NPZ) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.14 ) THEN
          IF( WL(1,NPZ).GT.EPSL ) THEN
            WLA(M,NPZ) = WL(M,NPZ)*FLA
          ELSE
            WLA(M,NPZ) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
     &      *RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMLA,INDX )
          WDLA(M,NPZ) = DFM*(XMLAB(MP,NB)
     &      - XMLA(MP,N))/(5.D-1*DZGF(N))
          WLA(M,NPZ) = WL(M,NPZ)*FLA + WTMA*WDLA(M,NPZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAE( N,NB )
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
!     Compute advective and molar diffusive air aqueous fluxes on 
!     an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFLAE'
      I = ID(N)
      NQX = NSX(N)+1
      DXMLA = XMLA(2,N) - XMLAB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FLAP = XLA(MN,N)*RHOL(MN,N)
        FLAB = XLAB(MN,NB)*RHOLB(MN,NB)
        INDX = 3
        FLA = DIFMN( FLAP,FLAB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( UL(1,NQX).GT.EPSL ) THEN
            ULA(M,NQX) = UL(M,NQX)*FLA
          ELSE
            ULA(M,NQX) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.14 ) THEN
          IF( UL(1,NQX).LT.-EPSL ) THEN
            ULA(M,NQX) = UL(M,NQX)*FLA
          ELSE
            ULA(M,NQX) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
     &      *RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMLA,INDX )
          UDLA(M,NQX) = DFM*(XMLA(MN,N)
     &      - XMLAB(MN,NB))/(5.D-1*DXGF(N))
          ULA(M,NQX) = UL(M,NQX)*FLA + WTMA*UDLA(M,NQX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAN( N,NB )
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
!     Compute advective and molar diffusive air aqueous fluxes on 
!     a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFLAN'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMLA = XMLA(2,N) - XMLAB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FLAP = XLA(MN,N)*RHOL(MN,N)
        FLAB = XLAB(MN,NB)*RHOLB(MN,NB)
        INDX = 3
        FLA = DIFMN( FLAP,FLAB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( VL(1,NQY).GT.EPSL ) THEN
            VLA(M,NQY) = VL(M,NQY)*FLA
          ELSE
            VLA(M,NQY) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.14 ) THEN
          IF( VL(1,NQY).LT.-EPSL ) THEN
            VLA(M,NQY) = VL(M,NQY)*FLA
          ELSE
            VLA(M,NQY) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
     &      *RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMLA,INDX )
          VDLA(M,NQY) = DFM*(XMLA(MN,N)
     &      - XMLAB(MN,NB))/((5.D-1*DYGF(N))*RP(I))
          VLA(M,NQY) = VL(M,NQY)*FLA + WTMA*VDLA(M,NQY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAS( N,NB )
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
!     Compute advective and molar diffusive air aqueous fluxes on 
!     a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFLAS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMLA = XMLAB(2,NB) - XMLA(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FLAP = XLA(MP,N)*RHOL(MP,N)
        FLAB = XLAB(MP,NB)*RHOLB(MP,NB)
        INDX = 3
        FLA = DIFMN( FLAB,FLAP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( VL(1,NPY).LT.-EPSL ) THEN
            VLA(M,NPY) = VL(M,NPY)*FLA
          ELSE
            VLA(M,NPY) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.14 ) THEN
          IF( VL(1,NPY).GT.EPSL ) THEN
            VLA(M,NPY) = VL(M,NPY)*FLA
          ELSE
            VLA(M,NPY) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
     &      *RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMLA,INDX )
          VDLA(M,NPY) = DFM*(XMLAB(MP,NB)
     &      - XMLA(MP,N))/((5.D-1*DYGF(N))*RP(I))
          VLA(M,NPY) = VL(M,NPY)*FLA + WTMA*VDLA(M,NPY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAT( N,NB )
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
!     Compute advective and molar diffusive air aqueous fluxes on 
!     a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFLAT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMLA = XMLA(2,N)-XMLAB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FLAP = XLA(MN,N)*RHOL(MN,N)
        FLAB = XLAB(MN,NB)*RHOLB(MN,NB)
        INDX = 3
        FLA = DIFMN( FLAP,FLAB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( WL(1,NQZ).GT.EPSL ) THEN
            WLA(M,NQZ) = WL(M,NQZ)*FLA
          ELSE
            WLA(M,NQZ) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.14 ) THEN
          IF( WL(1,NQZ).LT.-EPSL ) THEN
            WLA(M,NQZ) = WL(M,NQZ)*FLA
          ELSE
            WLA(M,NQZ) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
     &      *RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMLA,INDX )
          WDLA(M,NQZ) = DFM*(XMLA(MN,N)
     &      - XMLAB(MN,NB))/(5.D-1*DZGF(N))
          WLA(M,NQZ) = WL(M,NQZ)*FLA + WTMA*WDLA(M,NQZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAW( N,NB )
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
!     Compute advective and molar diffusive air aqueous fluxes on 
!     a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFLAW'
      I = ID(N)
      NPX = NSX(N)
      DXMLA = XMLAB(2,NB)-XMLA(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
         FLAP = XLA(MP,N)*RHOL(MP,N)
        FLAB = XLAB(MP,NB)*RHOLB(MP,NB)
        INDX = 3
        FLA = DIFMN( FLAB,FLAP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( UL(1,NPX).LT.-EPSL ) THEN
            ULA(M,NPX) = UL(M,NPX)*FLA
          ELSE
            ULA(M,NPX) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.14 ) THEN
          IF( UL(1,NPX).GT.EPSL ) THEN
            ULA(M,NPX) = UL(M,NPX)*FLA
          ELSE
            ULA(M,NPX) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
     &      *RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMLA,INDX )
          UDLA(M,NPX) = DFM*(XMLAB(MP,NB)
     &      - XMLA(MP,N))/(5.D-1*DXGF(N))
          ULA(M,NPX) = UL(M,NPX)*FLA + WTMA*UDLA(M,NPX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLOB( N,NB )
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
!     Compute advective and molar diffusive air aqueous fluxes on 
!     a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE CONST
      USE BCVP
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFLOB'
      K = KD(N)
      NPZ = NSZ(N)
      DXMLO = XMLOB(2,NB) - XMLO(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FLOP = XLO(MP,N)*RHOL(MP,N)
        FLOB = XLOB(MP,NB)*RHOLB(MP,NB)
        INDX = 3
        FLO = DIFMN( FLOB,FLOP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( WL(1,NPZ).LT.-EPSL ) THEN
            WLO(M,NPZ) = WL(M,NPZ)*FLO
          ELSE
            WLO(M,NPZ) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.14 ) THEN
          IF( WL(1,NPZ).GT.EPSL ) THEN
            WLO(M,NPZ) = WL(M,NPZ)*FLO
          ELSE
            WLO(M,NPZ) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLO(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLOB(MP,NB)
     &      *RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMLO,INDX )
          WDLO(M,NPZ) = DFM*(XMLOB(MP,NB)
     &      - XMLO(MP,N))/(5.D-1*DZGF(N))
          WLO(M,NPZ) = WL(M,NPZ)*FLO + WTMO*WDLO(M,NPZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLOB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLOE( N,NB )
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
!     Compute advective and molar diffusive air aqueous fluxes on 
!     an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE CONST
      USE BCVP
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFLOE'
      I = ID(N)
      NQX = NSX(N)+1
      DXMLO = XMLO(2,N) - XMLOB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FLOP = XLO(MN,N)*RHOL(MN,N)
        FLOB = XLOB(MN,NB)*RHOLB(MN,NB)
        INDX = 3
        FLO = DIFMN( FLOP,FLOB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( UL(1,NQX).GT.EPSL ) THEN
            ULO(M,NQX) = UL(M,NQX)*FLO
          ELSE
            ULO(M,NQX) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.14 ) THEN
          IF( UL(1,NQX).LT.-EPSL ) THEN
            ULO(M,NQX) = UL(M,NQX)*FLO
          ELSE
            ULO(M,NQX) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLO(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLOB(MN,NB)
     &      *RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMLO,INDX )
          UDLO(M,NQX) = DFM*(XMLO(MN,N)
     &      - XMLOB(MN,NB))/(5.D-1*DXGF(N))
          ULO(M,NQX) = UL(M,NQX)*FLO + WTMO*UDLO(M,NQX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLOE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLON( N,NB )
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
!     Compute advective and molar diffusive air aqueous fluxes on 
!     a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE CONST
      USE BCVP
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFLON'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMLO = XMLO(2,N) - XMLOB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FLOP = XLO(MN,N)*RHOL(MN,N)
        FLOB = XLOB(MN,NB)*RHOLB(MN,NB)
        INDX = 3
        FLO = DIFMN( FLOP,FLOB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( VL(1,NQY).GT.EPSL ) THEN
            VLO(M,NQY) = VL(M,NQY)*FLO
          ELSE
            VLO(M,NQY) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.14 ) THEN
          IF( VL(1,NQY).LT.-EPSL ) THEN
            VLO(M,NQY) = VL(M,NQY)*FLO
          ELSE
            VLO(M,NQY) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLO(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLOB(MN,NB)
     &      *RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMLO,INDX )
          VDLO(M,NQY) = DFM*(XMLO(MN,N)
     &      - XMLOB(MN,NB))/((5.D-1*DYGF(N))*RP(I))
          VLO(M,NQY) = VL(M,NQY)*FLO + WTMO*VDLO(M,NQY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLON group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLOS( N,NB )
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
!     Compute advective and molar diffusive air aqueous fluxes on 
!     a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE CONST
      USE BCVP
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFLOS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMLO = XMLOB(2,NB) - XMLO(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FLOP = XLO(MP,N)*RHOL(MP,N)
        FLOB = XLOB(MP,NB)*RHOLB(MP,NB)
        INDX = 3
        FLO = DIFMN( FLOB,FLOP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( VL(1,NPY).LT.-EPSL ) THEN
            VLO(M,NPY) = VL(M,NPY)*FLO
          ELSE
            VLO(M,NPY) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.14 ) THEN
          IF( VL(1,NPY).GT.EPSL ) THEN
            VLO(M,NPY) = VL(M,NPY)*FLO
          ELSE
            VLO(M,NPY) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLO(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLOB(MP,NB)
     &      *RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMLO,INDX )
          VDLO(M,NPY) = DFM*(XMLOB(MP,NB)
     &      - XMLO(MP,N))/((5.D-1*DYGF(N))*RP(I))
          VLO(M,NPY) = VL(M,NPY)*FLO + WTMO*VDLO(M,NPY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLOS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLOT( N,NB )
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
!     Compute advective and molar diffusive air aqueous fluxes on 
!     a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE CONST
      USE BCVP
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFLOT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMLO = XMLO(2,N)-XMLOB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FLOP = XLO(MN,N)*RHOL(MN,N)
        FLOB = XLOB(MN,NB)*RHOLB(MN,NB)
        INDX = 3
        FLO = DIFMN( FLOP,FLOB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( WL(1,NQZ).GT.EPSL ) THEN
            WLO(M,NQZ) = WL(M,NQZ)*FLO
          ELSE
            WLO(M,NQZ) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.14 ) THEN
          IF( WL(1,NQZ).LT.-EPSL ) THEN
            WLO(M,NQZ) = WL(M,NQZ)*FLO
          ELSE
            WLO(M,NQZ) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLO(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLOB(MN,NB)
     &      *RHOMLB(MN,NB)
          INDX = 14
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMLO,INDX )
          WDLO(M,NQZ) = DFM*(XMLO(MN,N)
     &      - XMLOB(MN,NB))/(5.D-1*DZGF(N))
          WLO(M,NQZ) = WL(M,NQZ)*FLO + WTMO*WDLO(M,NQZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLOT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLOW( N,NB )
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
!     Compute advective and molar diffusive air aqueous fluxes on 
!     a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE CONST
      USE BCVP
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFLOW'
      I = ID(N)
      NPX = NSX(N)
      DXMLO = XMLOB(2,NB)-XMLO(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
         FLOP = XLO(MP,N)*RHOL(MP,N)
        FLOB = XLOB(MP,NB)*RHOLB(MP,NB)
        INDX = 3
        FLO = DIFMN( FLOB,FLOP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(2,NB).EQ.7 .OR. IBCT(2,NB).EQ.19 ) THEN
          IF( UL(1,NPX).LT.-EPSL ) THEN
            ULO(M,NPX) = UL(M,NPX)*FLO
          ELSE
            ULO(M,NPX) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.13 .OR. IBCT(2,NB).EQ.14 ) THEN
          IF( UL(1,NPX).GT.EPSL ) THEN
            ULO(M,NPX) = UL(M,NPX)*FLO
          ELSE
            ULO(M,NPX) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLO(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLOB(MP,NB)
     &      *RHOMLB(MP,NB)
          INDX = 14
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMLO,INDX )
          UDLO(M,NPX) = DFM*(XMLOB(MP,NB)
     &      - XMLO(MP,N))/(5.D-1*DXGF(N))
          ULO(M,NPX) = UL(M,NPX)*FLO + WTMO*UDLO(M,NPX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLOW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLWB( N,NB )
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
!     Compute advective and molar diffusive water aqueous fluxes on 
!     on bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLWB'
      K = KD(N)
      NPZ = NSZ(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FLWP = XLW(MP,N)*RHOL(MP,N)
        FLWB = XLWB(MP,NB)*RHOLB(MP,NB)
        INDX = 3
        FLW = DIFMN( FLWB,FLWP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
        WLW(M,NPZ) = WL(M,NPZ)*FLW - WTMW*(WDLA(M,NPZ)+WDLO(M,NPZ))
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLWB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLWE( N,NB )
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
!     Compute advective and molar diffusive water aqueous fluxes on 
!     an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLWE'
      I = ID(N)
      NQX = NSX(N)+1
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FLWP = XLW(MN,N)*RHOL(MN,N)
        FLWB = XLWB(MN,NB)*RHOLB(MN,NB)
        INDX = 3
        FLW = DIFMN( FLWP,FLWB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
        ULW(M,NQX) = UL(M,NQX)*FLW - WTMW*(UDLA(M,NQX)+UDLO(M,NQX))
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLWE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLWN( N,NB )
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
!     Compute advective and molar diffusive water aqueous fluxes on 
!     a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLWN'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FLWP = XLW(MN,N)*RHOL(MN,N)
        FLWB = XLWB(MN,NB)*RHOLB(MN,NB)
        INDX = 3
        FLW = DIFMN( FLWP,FLWB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
        VLW(M,NQY) = VL(M,NQY)*FLW - WTMW*(VDLA(M,NQY)+VDLO(M,NQY))
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLWN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLWS( N,NB )
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
!     Compute advective and molar diffusive water aqueous fluxes on 
!     a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLWS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FLWP = XLW(MP,N)*RHOL(MP,N)
        FLWB = XLWB(MP,NB)*RHOLB(MP,NB)
        INDX = 3
        FLW = DIFMN( FLWB,FLWP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
        VLW(M,NPY) = VL(M,NPY)*FLW - WTMW*(VDLA(M,NPY)+VDLO(M,NPY))
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLWS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLWT( N,NB )
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
!     Compute advective and molar diffusive water aqueous fluxes on 
!     a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLWT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FLWP = XLW(MN,N)*RHOL(MN,N)
        FLWB = XLWB(MN,NB)*RHOLB(MN,NB)
        INDX = 3
        FLW = DIFMN( FLWP,FLWB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
        WLW(M,NQZ) = WL(M,NQZ)*FLW - WTMW*(WDLA(M,NQZ)+WDLO(M,NQZ))
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLWT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLWW( N,NB )
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
!     Compute advective and molar diffusive water aqueous fluxes on 
!     a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLWW'
      I = ID(N)
      NPX = NSX(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FLWP = XLW(MP,N)*RHOL(MP,N)
        FLWB = XLWB(MP,NB)*RHOLB(MP,NB)
        INDX = 3
        FLW = DIFMN( FLWB,FLWP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
        ULW(M,NPX) = UL(M,NPX)*FLW - WTMW*(UDLA(M,NPX)+UDLO(M,NPX))
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLWW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNOB( N,NB )
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
!     Compute advective oil NAPL fluxes on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVA
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFNOB'
      K = KD(N)
      NPZ = NSZ(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FNOP = XNO(MP,N)*RHON(MP,N)
        FNOB = XNOB(MP,NB)*RHONB(MP,NB)
        INDX = 3
        FNO = DIFMN( FNOB,FNOP,DZGF(N),DZGF(N),WN(1,NPZ),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(3,NB).EQ.7 .OR. IBCT(3,NB).EQ.19 ) THEN
          IF( WN(1,NPZ).LT.-EPSL ) THEN
            WNA(M,NPZ) = WN(M,NPZ)*FNO
          ELSE
            WNA(M,NPZ) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.13 .OR. IBCT(3,NB).EQ.16 ) THEN
          IF( WN(1,NPZ).GT.EPSL ) THEN
            WNA(M,NPZ) = WN(M,NPZ)*FNO
          ELSE
            WNA(M,NPZ) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          WNA(M,NPZ) = WN(M,NPZ)*FNO
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNOB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNOE( N,NB )
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
!     Compute advective oil NAPL fluxes on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVA
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFNOE'
      I = ID(N)
      NQX = NSX(N)+1
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FNOP = XNO(MN,N)*RHON(MN,N)
        FNOB = XNOB(MN,NB)*RHONB(MN,NB)
        INDX = 3
        FNO = DIFMN( FNOP,FNOB,DXGF(N),DXGF(N),UN(1,NQX),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(3,NB).EQ.7 .OR. IBCT(3,NB).EQ.19 ) THEN
          IF( UN(1,NQX).GT.EPSL ) THEN
            UNA(M,NQX) = UN(M,NQX)*FNO
          ELSE
            UNA(M,NQX) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.13 .OR. IBCT(3,NB).EQ.16 ) THEN
          IF( UN(1,NQX).LT.-EPSL ) THEN
            UNA(M,NQX) = UN(M,NQX)*FNO
          ELSE
            UNA(M,NQX) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          UNA(M,NQX) = UN(M,NQX)*FNO
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNOE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNON( N,NB )
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
!     Compute advective oil NAPL fluxes on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVA
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFNON'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FNOP = XNO(MN,N)*RHON(MN,N)
        FNOB = XNOB(MN,NB)*RHONB(MN,NB)
        INDX = 3
        FNO = DIFMN( FNOP,FNOB,DYGF(N),DYGF(N),VN(1,NQY),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(3,NB).EQ.7 .OR. IBCT(3,NB).EQ.19 ) THEN
          IF( VN(1,NQY).GT.EPSL ) THEN
            VNA(M,NQY) = VN(M,NQY)*FNO
          ELSE
            VNA(M,NQY) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.13 .OR. IBCT(3,NB).EQ.16 ) THEN
          IF( VN(1,NQY).LT.-EPSL ) THEN
            VNA(M,NQY) = VN(M,NQY)*FNO
          ELSE
            VNA(M,NQY) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          VNA(M,NQY) = VN(M,NQY)*FNO
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNON group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNOS( N,NB )
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
!     Compute advective oil NAPL fluxes on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVA
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFNOS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FNOP = XNO(MP,N)*RHON(MP,N)
        FNOB = XNOB(MP,NB)*RHONB(MP,NB)
        INDX = 3
        FNO = DIFMN( FNOB,FNOP,DYGF(N),DYGF(N),VN(1,NPY),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(3,NB).EQ.7 .OR. IBCT(3,NB).EQ.19 ) THEN
          IF( VN(1,NPY).LT.-EPSL ) THEN
            VNA(M,NPY) = VN(M,NPY)*FNO
          ELSE
            VNA(M,NPY) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.13 .OR. IBCT(3,NB).EQ.16 ) THEN
          IF( VN(1,NPY).GT.EPSL ) THEN
            VNA(M,NPY) = VN(M,NPY)*FNO
          ELSE
            VNA(M,NPY) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          VNA(M,NPY) = VN(M,NPY)*FNO
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNOS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNOT( N,NB )
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
!     Compute advective oil NAPL fluxes on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVA
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFNOT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FNOP = XNO(MN,N)*RHON(MN,N)
        FNOB = XNOB(MN,NB)*RHONB(MN,NB)
        INDX = 3
        FNO = DIFMN( FNOP,FNOB,DZGF(N),DZGF(N),WN(1,NQZ),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(3,NB).EQ.7 .OR. IBCT(3,NB).EQ.19 ) THEN
          IF( WN(1,NQZ).GT.EPSL ) THEN
            WNA(M,NQZ) = WN(M,NQZ)*FNO
          ELSE
            WNA(M,NQZ) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.13 .OR. IBCT(3,NB).EQ.16 ) THEN
          IF( WN(1,NQZ).LT.-EPSL ) THEN
            WNA(M,NQZ) = WN(M,NQZ)*FNO
          ELSE
            WNA(M,NQZ) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          WNA(M,NQZ) = WN(M,NQZ)*FNO
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNOT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNOW( N,NB )
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
!     Compute advective oil NAPL fluxes on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_woa.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVA
      USE CONST
      USE BCVP
      USE BCVA
      USE BCV
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
      SUB_LOG(ISUB_LOG) = '/DFFNOW'
      I = ID(N)
      NPX = NSX(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FNOP = XNO(MP,N)*RHON(MP,N)
        FNOB = XNOB(MP,NB)*RHONB(MP,NB)
        INDX = 3
        FNO = DIFMN( FNOB,FNOP,DXGF(N),DXGF(N),UN(1,NPX),INDX )
!
!---    Outflow-type boundary condition  ---
!
        IF( IBCT(3,NB).EQ.7 .OR. IBCT(3,NB).EQ.19 ) THEN
          IF( UN(1,NPX).LT.-EPSL ) THEN
            UNA(M,NPX) = UN(M,NPX)*FNO
          ELSE
            UNA(M,NPX) = 0.D+0
          ENDIF
!
!---    Inflow-type boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.13 .OR. IBCT(3,NB).EQ.16 ) THEN
          IF( UN(1,NPX).GT.EPSL ) THEN
            UNA(M,NPX) = UN(M,NPX)*FNO
          ELSE
            UNA(M,NPX) = 0.D+0
          ENDIF
!
!---    Not outflow- nor inflow-type boundary conditions  ---
!
        ELSE
          UNA(M,NPX) = UN(M,NPX)*FNO
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNOW group  ---
!
      RETURN
      END


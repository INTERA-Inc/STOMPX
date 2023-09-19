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
!     Compute diffusive CO2 gas fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_hyd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
          DXGA = XMGA(2,NW) - XMGA(2,N)
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            IF( IXP(NW).LT.0 ) MN = 2
            IF( IXP(N).LT.0 ) MP = 2
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGA(MP,N)*RHOMG(MP,N)
            DFW = TORG(MN,NW)*PORD(MN,NW)*SG(MN,NW)*
     &        DFGA(MN,NW)*RHOMG(MN,NW)
            INDX = 20
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXGA,INDX )
            UDGA(M,NPX) = DFM*(XMGA(MN,NW)
     &        - XMGA(MP,N))/DXGP(NPX)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            FGAW = XGA(MN,NW)*RHOG(MN,NW)
            INDX = 3
            FGA = DIFMN( FGAW,FGAP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            UGA(M,NPX) = UG(M,NPX)*FGA + WTMA*UDGA(M,NPX)
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
          DXGA = XMGA(2,NS) - XMGA(2,N)
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NS).LT.0 ) MN = 2
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGA(MP,N)*RHOMG(MP,N)
            DFS = TORG(MN,NS)*PORD(MN,NS)*SG(MN,NS)*
     &        DFGA(MN,NS)*RHOMG(MN,NS)
            INDX = 20
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXGA,INDX )
            VDGA(M,NPY) = DFM*(XMGA(MN,NS)
     &        - XMGA(MP,N))/(DYGP(NPY)*RP(I))
            FGAP = XGA(MP,N)*RHOG(MP,N)
            FGAS = XGA(MN,NS)*RHOG(MN,NS)
            INDX = 3
            FGA = DIFMN( FGAS,FGAP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
            VGA(M,NPY) = VG(M,NPY)*FGA + WTMA*VDGA(M,NPY)
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
          DXGA = XMGA(2,NB) - XMGA(2,N)
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NB).LT.0 ) MN = 2
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGA(MP,N)*RHOMG(MP,N)
            DFB = TORG(MN,NB)*PORD(MN,NB)*SG(MN,NB)*
     &        DFGA(MN,NB)*RHOMG(MN,NB)
            INDX = 20
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXGA,INDX )
            WDGA(M,NPZ) = DFM*(XMGA(MN,NB)
     &        - XMGA(MP,N))/DZGP(NPZ)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            FGAB = XGA(MN,NB)*RHOG(MN,NB)
            INDX = 3
            FGA = DIFMN( FGAB,FGAP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
            WGA(M,NPZ) = WG(M,NPZ)*FGA + WTMA*WDGA(M,NPZ)
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
!     Compute diffusive oil gas fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_hyd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
          NPX = NSX(N)
          DXGO = XMGO(2,NW) - XMGO(2,N)
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NW).LT.0 ) MN = 2
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGO(MP,N)*RHOMG(MP,N)
            DFW = TORG(MN,NW)*PORD(MN,NW)*SG(MN,NW)*
     &        DFGO(MN,NW)*RHOMG(MN,NW)
            INDX = 13
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXGO,INDX )
            UDGO(M,NPX) = DFM*(XMGO(MN,NW)
     &        - XMGO(MP,N))/DXGP(NPX)
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
          DXGO = XMGO(2,NS) - XMGO(2,N)
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NS).LT.0 ) MN = 2
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGO(MP,N)*RHOMG(MP,N)
            DFS = TORG(MN,NS)*PORD(MN,NS)*SG(MN,NS)*
     &        DFGO(MN,NS)*RHOMG(MN,NS)
            INDX = 13
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXGO,INDX )
            VDGO(M,NPY) = DFM*(XMGO(MN,NS)
     &        - XMGO(MP,N))/(DYGP(NPY)*RP(I))
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
          DXGO = XMGO(2,NB) - XMGO(2,N)
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NB).LT.0 ) MN = 2
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGO(MP,N)*RHOMG(MP,N)
            DFB = TORG(MN,NB)*PORD(MN,NB)*SG(MN,NB)*
     &        DFGO(MN,NB)*RHOMG(MN,NB)
            INDX = 13
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXGO,INDX )
            WDGO(M,NPZ) = DFM*(XMGO(MN,NB)
     &        - XMGO(MP,N))/DZGP(NPZ)
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
!     Compute diffusive CO2 gas fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_hyd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NW).LT.0 ) MN = 2
            FGWP = XGW(MP,N)*RHOG(MP,N)
            FGWW = XGW(MN,NW)*RHOG(MN,NW)
            INDX = 3
            FGW = DIFMN( FGWW,FGWP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            UDGW(M,NPX) = -(UDGA(M,NPX)+UDGO(M,NPX))
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
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NS).LT.0 ) MN = 2
            FGWP = XGW(MP,N)*RHOG(MP,N)
            FGWS = XGW(MN,NS)*RHOG(MN,NS)
            INDX = 3
            FGW = DIFMN( FGWS,FGWP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
            VDGW(M,NPY) = -(VDGA(M,NPY)+VDGO(M,NPY))
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
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NB).LT.0 ) MN = 2
            FGWP = XGW(MP,N)*RHOG(MP,N)
            FGWB = XGW(MN,NB)*RHOG(MN,NB)
            INDX = 3
            FGW = DIFMN( FGWB,FGWP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
            WDGW(M,NPZ) = -(WDGA(M,NPZ)+WDGO(M,NPZ))
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
!     Compute dissolved air molar diffusion rates through the
!     aqueous phase.
!     dff_hyd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_hyd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
!     Compute dissolved oil fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_hyd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
          DXMLO = (XMLO(2,NW)-XMLO(2,N))
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NW).LT.0 ) MN = 2
            DFSP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLO(MP,N)
     &        *RHOML(MP,N)
            DFSW = TORL(MN,NW)*SL(MN,NW)*PORD(MN,NW)*DFLO(MN,NW)
     &        *RHOML(MN,NW)
            INDX = 15
            DFM = DIFMN( DFSW,DFSP,DXGF(NW),DXGF(N),DXMLO,INDX )
            UDLO(M,NPX) = DFM*(XMLO(MN,NW)
     &        - XMLO(MP,N))/DXGP(NPX)
            FLOP = XLO(MP,N)*RHOL(MP,N)
            FLOW = XLO(MN,NW)*RHOL(MN,NW)
            INDX = 2
            FLO = DIFMN( FLOW,FLOP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            ULO(M,NPX) = UL(M,NPX)*FLO + WTMO*UDLO(M,NPX)
  100     CONTINUE
  200   CONTINUE
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
          DXMLO = (XMLO(2,NS)-XMLO(2,N))
          DO 300 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NS).LT.0 ) MN = 2
            DFSP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLO(MP,N)
     &        *RHOML(MP,N)
            DFSS = TORL(MN,NS)*SL(MN,NS)*PORD(MN,NS)*DFLO(MN,NS)
     &        *RHOML(MN,NS)
            INDX = 15
            DFM = DIFMN( DFSS,DFSP,DYGF(NS),DYGF(N),DXMLO,INDX )
            VDLO(M,NPY) = DFM*(XMLO(MN,NS)
     &        - XMLO(MP,N))/(DYGP(NPY)*RP(I))
            FLOP = XLO(MP,N)*RHOL(MP,N)
            FLOS = XLO(MN,NS)*RHOL(MN,NS)
            INDX = 2
            FLO = DIFMN( FLOS,FLOP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
            VLO(M,NPY) = VL(M,NPY)*FLO + WTMO*VDLO(M,NPY)
  300     CONTINUE
  400   CONTINUE
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
          DXMLO = (XMLO(2,NB)-XMLO(2,N))
          DO 500 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NB).LT.0 ) MN = 2
            DFSP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLO(MP,N)
     &        *RHOML(MP,N)
            DFSB = TORL(MN,NB)*SL(MN,NB)*PORD(MN,NB)*DFLO(MN,NB)
     &        *RHOML(MN,NB)
            INDX = 15
            DFM = DIFMN( DFSB,DFSP,DZGF(NB),DZGF(N),DXMLO,INDX )
            WDLO(M,NPZ) = DFM*(XMLO(MN,NB)
     &        - XMLO(MP,N))/DZGP(NPZ)
            FLOP = XLO(MP,N)*RHOL(MP,N)
            FLOB = XLO(MN,NB)*RHOL(MN,NB)
            INDX = 2
            FLO = DIFMN( FLOB,FLOP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
            WLO(M,NPZ) = WL(M,NPZ)*FLO + WTMO*WDLO(M,NPZ)
  500     CONTINUE
  600   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLO group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLS
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
!     Compute salt aqueous-phase fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_hyd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
      USE FLUXP
      USE FDVS
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
      SUB_LOG(ISUB_LOG) = '/DFFLS'
!
!---  X-direction Darcy velocities, excluding boundaries
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
!
!---        Diffusion coefficients  ---
!
            IF( IEDLS.EQ.1 ) THEN
              TCOR = (T(MP,N)+TABS)/TSPRF
              SMDLP = DFLS(MP,N)*TCOR*(VISRL/VISL(MP,N))
              DFCLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
              TCOR = (T(MN,NW)+TABS)/TSPRF
              SMDLP = DFLS(MN,NW)*TCOR*(VISRL/VISL(MN,NW))
              DFCLW = TORL(MN,NW)*SL(MN,NW)*PORD(MN,NW)*SMDLP
            ELSEIF( IEDLS.EQ.2 ) THEN
              DFCLP = SDCLS(1,IZ(N))*SDCLS(2,IZ(N))*
     &          EXP(SL(MP,N)*PORD(MP,N)*SDCLS(3,IZ(N)))
              DFCLW = SDCLS(1,IZ(NW))*SDCLS(2,IZ(NW))*
     &          EXP(SL(MN,NW)*PORD(MN,NW)*SDCLS(3,IZ(NW)))
            ELSEIF( IEDLS.EQ.3 ) THEN
              DFCLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
              DFCLW = TORL(MN,NW)*SL(MN,NW)*PORD(MN,NW)*DFLS(MN,NW)
            ENDIF
            INDX = 18
            DFCLW = DIFMN(DFCLW,DFCLP,DXGF(NW),DXGF(N),UL(1,NPX),INDX)
!
!---        Salt aqueous flux by advection and diffusion  ---
!
            DDLW = DFCLW/DXGP(NPX)
            AL = MAX( UL(M,NPX),ZERO ) +
     &        DDLW*MAX((ONE-(TENTH*ABS(UL(M,NPX))/
     &        (DDLW+SMALL)))**5,ZERO)
            ALP = MAX( -UL(M,NPX),ZERO ) +
     &        DDLW*MAX((ONE-(TENTH*ABS(UL(M,NPX))/
     &        (DDLW+SMALL)))**5,ZERO)
            US(M,NPX) = XLS(MN,NW)*RHOL(MN,NW)*AL -
     &        XLS(MP,N)*RHOL(MP,N)*ALP
            UDS(M,NPX) = DDLW*(XLS(MN,NW)*RHOL(MN,NW) -
     &        XLS(MP,N)*RHOL(MP,N))
  100     CONTINUE
  200   CONTINUE
      ENDIF
!
!---  Y-direction Darcy velocities, excluding boundaries
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
            MN = MNEG(M)
            MP = MPOS(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NS).LT.0 ) MN = 2
!
!---        Diffusion coefficients  ---
!
            IF( IEDLS.EQ.1 ) THEN
              TCOR = (T(MP,N)+TABS)/TSPRF
              SMDLP = DFLS(MP,N)*TCOR*(VISRL/VISL(MP,N))
              DFCLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
              TCOR = (T(MN,NS)+TABS)/TSPRF
              SMDLP = DFLS(MN,NS)*TCOR*(VISRL/VISL(MN,NS))
              DFCLS = TORL(MN,NS)*SL(MN,NS)*PORD(MN,NS)*SMDLP
            ELSEIF( IEDLS.EQ.2 ) THEN
              DFCLP = SDCLS(1,IZ(N))*SDCLS(2,IZ(N))*
     &          EXP(SL(MP,N)*PORD(MP,N)*SDCLS(3,IZ(N)))
              DFCLS = SDCLS(1,IZ(NS))*SDCLS(2,IZ(NS))*
     &          EXP(SL(MN,NS)*PORD(MN,NS)*SDCLS(3,IZ(NS)))
            ELSEIF( IEDLS.EQ.3 ) THEN
              DFCLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
              DFCLS = TORL(MN,NS)*SL(MN,NS)*PORD(MN,NS)*DFLS(MN,NS)
            ENDIF
            INDX = 18
            DFCLS = DIFMN(DFCLS,DFCLP,DYGF(NS),DYGF(N),VL(1,NPY),INDX)
!
!---        Salt aqueous flux by advection and diffusion  ---
!
            DDLS = DFCLS /(DYGP(NPY)*RP(I))
            AL = MAX( VL(M,NPY),ZERO ) +
     &        DDLS*MAX((ONE-(TENTH*ABS(VL(M,NPY))/
     &        (DDLS+SMALL)))**5,ZERO)
            ALP = MAX( -VL(M,NPY),ZERO ) +
     &        DDLS*MAX((ONE-(TENTH*ABS(VL(M,NPY))/
     &        (DDLS+SMALL)))**5,ZERO)
            VS(M,NPY) = (XLS(MN,NS)*RHOL(MN,NS)*AL -
     &        XLS(MP,N)*RHOL(MP,N)*ALP)
            VDS(M,NPY) = DDLS*(XLS(MN,NS)*RHOL(MN,NS) -
     &        XLS(MP,N)*RHOL(MP,N))
  300     CONTINUE
  400   CONTINUE
      ENDIF
!
!---  Z-direction Darcy velocities, excluding boundaries
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
            MN = MNEG(M)
            MP = MPOS(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NB).LT.0 ) MN = 2
!
!---        Diffusion coefficients  ---
!
            IF( IEDLS.EQ.1 ) THEN
              TCOR = (T(MP,N)+TABS)/TSPRF
              SMDLP = DFLS(MP,N)*TCOR*(VISRL/VISL(MP,N))
              DFCLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
              TCOR = (T(MN,NB)+TABS)/TSPRF
              SMDLP = DFLS(MN,NB)*TCOR*(VISRL/VISL(MN,NB))
              DFCLB = TORL(MN,NB)*SL(MN,NB)*PORD(MN,NB)*SMDLP
            ELSEIF( IEDLS.EQ.2 ) THEN
              DFCLP = SDCLS(1,IZ(N))*SDCLS(2,IZ(N))*
     &          EXP(SL(MP,N)*PORD(MP,N)*SDCLS(3,IZ(N)))
              DFCLB = SDCLS(1,IZ(NB))*SDCLS(2,IZ(NB))*
     &          EXP(SL(MN,NB)*PORD(MN,NB)*SDCLS(3,IZ(NB)))
            ELSEIF( IEDLS.EQ.3 ) THEN
              DFCLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
              DFCLB = TORL(MN,NB)*SL(MN,NB)*PORD(MN,NB)*DFLS(MN,NB)
            ENDIF
            INDX = 18
            DFCLB = DIFMN(DFCLB,DFCLP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX)
!
!---        Salt aqueous flux by advection and diffusion  ---
!
            DDLB = DFCLB/DZGP(NPZ)
            AL = MAX( WL(M,NPZ),ZERO ) +
     &        DDLB*MAX((ONE-(TENTH*ABS(WL(M,NPZ))/
     &        (DDLB+SMALL)))**5,ZERO)
            ALP = MAX( -WL(M,NPZ),ZERO ) +
     &        DDLB*MAX((ONE-(TENTH*ABS(WL(M,NPZ))/
     &        (DDLB+SMALL)))**5,ZERO)
            WS(M,NPZ) = (XLS(MN,NB)*RHOL(MN,NB)*AL -
     &        XLS(MP,N)*RHOL(MP,N)*ALP)
            WDS(M,NPZ) = DDLB*(XLS(MN,NB)*RHOL(MN,NB) -
     &        XLS(MP,N)*RHOL(MP,N))
  500     CONTINUE
  600   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLS group  ---
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
!     Compute dissolved air molar diffusion rates through the
!     aqueous phase.
!     dff_hyd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_hyd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
            ULW(M,NPX) = UL(M,NPX)*FLW 
     &        - WTMW*(UDLA(M,NPX)+UDLO(M,NPX)+UDS(M,NPX)/WTMS)
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
            VLW(M,NPY) = VL(M,NPY)*FLW 
     &        - WTMW*(VDLA(M,NPY)+VDLO(M,NPY)+VDS(M,NPY)/WTMS)
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
            WLW(M,NPZ) = WL(M,NPZ)*FLW 
     &        - WTMW*(WDLA(M,NPZ)+WDLO(M,NPZ)+WDS(M,NPZ)/WTMS)
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
      SUBROUTINE DFFNA
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
!     Compute diffusive CO2 liquid-CO2 fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 March 2010.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFNA'
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
            FNAP = XNC(1,MP,N)*RHON(MP,N)
            FNAW = XNC(1,MN,NW)*RHON(MN,NW)
            INDX = 3
            FNA = DIFMN( FNAW,FNAP,DXGF(NW),DXGF(N),UN(1,NPX),INDX )
            UDNC(1,M,NPX) = -(UDNC(3,M,NPX)+UDNC(2,M,NPX))
            UNC(1,M,NPX) = UN(M,NPX)*FNA + WTMA*UDNC(1,M,NPX)
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
            FNAP = XNC(1,MP,N)*RHON(MP,N)
            FNAS = XNC(1,MN,NS)*RHON(MN,NS)
            INDX = 3
            FNA = DIFMN( FNAS,FNAP,DYGF(NS),DYGF(N),VN(1,NPY),INDX )
            VDNC(1,M,NPY) = -(VDNC(3,M,NPY)+VDNC(2,M,NPY))
            VNC(1,M,NPY) = VN(M,NPY)*FNA + WTMA*VDNC(1,M,NPY)
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
            FNAP = XNC(1,MP,N)*RHON(MP,N)
            FNAB = XNC(1,MN,NB)*RHON(MN,NB)
            INDX = 3
            FNA = DIFMN( FNAB,FNAP,DZGF(NB),DZGF(N),WN(1,NPZ),INDX )
            WDNC(1,M,NPZ) = -(WDNC(3,M,NPZ)+WDNC(2,M,NPZ))
            WNC(1,M,NPZ) = WN(M,NPZ)*FNA + WTMA*WDNC(1,M,NPZ)
  500     CONTINUE
  600   CONTINUE
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNA group  ---
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
!     Compute diffusive CH4 liquid-CO2 fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 March 2010.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVN
      USE FDVGC
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
          DXNO = XMNC(2,2,NW) - XMNC(2,2,N)
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NW).LT.0 ) MN = 2
            DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &        DFNC(2,MP,N)*RHOMN(MP,N)
            DFW = TORN(MN,NW)*PORD(MN,NW)*SN(MN,NW)*
     &        DFNC(2,MN,NW)*RHOMN(MN,NW)
            INDX = 13
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXNO,INDX )
            UDNC(2,M,NPX) = DFM*(XMNC(2,MN,NW)
     &        - XMNC(2,MP,N))/DXGP(NPX)
            FNOP = XNC(2,MP,N)*RHON(MP,N)
            FNOW = XNC(2,MN,NW)*RHON(MN,NW)
            INDX = 3
            FNO = DIFMN( FNOW,FNOP,DXGF(NW),DXGF(N),UN(1,NPX),INDX )
            UNC(2,M,NPX) = UN(M,NPX)*FNO + WTMO*UDNC(2,M,NPX)
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
          DXNO = XMNC(2,2,NS) - XMNC(2,2,N)
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NS).LT.0 ) MN = 2
            DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &        DFNC(2,MP,N)*RHOMN(MP,N)
            DFS = TORN(MN,NS)*PORD(MN,NS)*SN(MN,NS)*
     &        DFNC(2,MN,NS)*RHOMN(MN,NS)
            INDX = 13
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXNO,INDX )
            VDNC(2,M,NPY) = DFM*(XMNC(2,MN,NS)
     &        - XMNC(2,MP,N))/(DYGP(NPY)*RP(I))
            FNOP = XNC(2,MP,N)*RHON(MP,N)
            FNOS = XNC(2,MN,NS)*RHON(MN,NS)
            INDX = 3
            FNO = DIFMN( FNOS,FNOP,DYGF(NS),DYGF(N),VN(1,NPY),INDX )
            VNC(2,M,NPY) = VN(M,NPY)*FNO + WTMO*VDNC(2,M,NPY)
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
          DXNO = XMNC(2,2,NB) - XMNC(2,2,N)
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NB).LT.0 ) MN = 2
            DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &        DFNC(2,MP,N)*RHOMN(MP,N)
            DFB = TORN(MN,NB)*PORD(MN,NB)*SN(MN,NB)*
     &        DFNC(2,MN,NB)*RHOMN(MN,NB)
            INDX = 13
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXNO,INDX )
            WDNC(2,M,NPZ) = DFM*(XMNC(2,MN,NB)
     &        - XMNC(2,MP,N))/DZGP(NPZ)
            FNOP = XNC(2,MP,N)*RHON(MP,N)
            FNOB = XNC(2,MN,NB)*RHON(MN,NB)
            INDX = 3
            FNO = DIFMN( FNOB,FNOP,DZGF(NB),DZGF(N),WN(1,NPZ),INDX )
            WNC(2,M,NPZ) = WN(M,NPZ)*FNO + WTMO*WDNC(2,M,NPZ)
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
      SUBROUTINE DFFNW
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
!     Compute diffusive H2O liquid-CO2 fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 March 2010.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVN
      USE FDVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFNW'
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
          DXNW = XMNC(3,2,NW) - XMNC(3,2,N)
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            IF( IXP(NW).LT.0 ) MN = 2
            IF( IXP(N).LT.0 ) MP = 2
            DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &        DFNC(3,MP,N)*RHOMN(MP,N)
            DFW = TORN(MN,NW)*PORD(MN,NW)*SN(MN,NW)*
     &        DFNC(3,MN,NW)*RHOMN(MN,NW)
            INDX = 12
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXNW,INDX )
            UDNC(3,M,NPX) = DFM*(XMNC(3,MN,NW)
     &        - XMNC(3,MP,N))/DXGP(NPX)
            FNWP = XNC(3,MP,N)*RHON(MP,N)
            FNWW = XNC(3,MN,NW)*RHON(MN,NW)
            INDX = 3
            FNW = DIFMN( FNWW,FNWP,DXGF(NW),DXGF(N),UN(1,NPX),INDX )
            UNC(3,M,NPX) = UN(M,NPX)*FNW + WTMW*UDNC(3,M,NPX)
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
          DXNW = XMNC(3,2,NS) - XMNC(3,2,N)
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NS).LT.0 ) MN = 2
            DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &        DFNC(3,MP,N)*RHOMN(MP,N)
            DFS = TORN(MN,NS)*PORD(MN,NS)*SN(MN,NS)*
     &        DFNC(3,MN,NS)*RHOMN(MN,NS)
            INDX = 12
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXNW,INDX )
            VDNC(3,M,NPY) = DFM*(XMNC(3,MN,NS)
     &        - XMNC(3,MP,N))/(DYGP(NPY)*RP(I))
            FNWP = XNC(3,MP,N)*RHON(MP,N)
            FNWS = XNC(3,MN,NS)*RHON(MN,NS)
            INDX = 3
            FNW = DIFMN( FNWS,FNWP,DYGF(NS),DYGF(N),VN(1,NPY),INDX )
            VNC(3,M,NPY) = VN(M,NPY)*FNW + WTMW*VDNC(3,M,NPY)
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
          DXNW = XMNC(3,2,NB) - XMNC(3,2,N)
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            IF( IXP(N).LT.0 ) MP = 2
            IF( IXP(NB).LT.0 ) MN = 2
            DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &        DFNC(3,MP,N)*RHOMN(MP,N)
            DFB = TORN(MN,NB)*PORD(MN,NB)*SN(MN,NB)*
     &        DFNC(3,MN,NB)*RHOMN(MN,NB)
            INDX = 12
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXNW,INDX )
            WDNC(3,M,NPZ) = DFM*(XMNC(3,MN,NB)
     &        - XMNC(3,MP,N))/DZGP(NPZ)
            FNWP = XNC(3,MP,N)*RHON(MP,N)
            FNWB = XNC(3,MN,NB)*RHON(MN,NB)
            INDX = 3
            FNW = DIFMN( FNWB,FNWP,DZGF(NB),DZGF(N),WN(1,NPZ),INDX )
            WNC(3,M,NPZ) = WN(M,NPZ)*FNW + WTMW*WDNC(3,M,NPZ)
  500     CONTINUE
  600   CONTINUE
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNW group  ---
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
!     Compute CO2 gas fluxes on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE BCVP
      USE BCVH
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
      SUB_LOG(ISUB_LOG) = '/DFFGAB'
      K = KD(N)
      NPZ = NSZ(N)
      DXMGA = XMGAB(2,NB) - XMGA(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Vapor micro emulsion  ---
!
        IF( ZLAB(MP,NB).LT.-EPSL ) THEN
          FLAP = 0.D+0
          FLAB = XGAB(MP,NB)*RHOGB(MP,NB)*(-ZLAB(MP,NB))
          INDX = 3
          FLA = DIFMN( FLAB,FLAP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLA = 0.D+0
        ENDIF
        FGAP = XGA(MP,N)*RHOG(MP,N)
        FGAB = XGAB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGA = DIFMN( FGAB,FGAP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.26 ) THEN
          IF( WG(1,NPZ).LT.-EPSL ) THEN
            WGA(M,NPZ) = WG(M,NPZ)*FGA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.42 ) THEN
          IF( WG(1,NPZ).GT.EPSL ) THEN
            WGA(M,NPZ) = WG(M,NPZ)*FGA
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGA(MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGAB(MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 20
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMGA,INDX )
          WDGA(M,NPZ) = DFM*(XMGAB(MP,NB)
     &      - XMGA(MP,N))/(5.D-1*DZGF(N))
          WGA(M,NPZ) = WL(M,NPZ)*FLA + WG(M,NPZ)*FGA + WTMA*WDGA(M,NPZ)
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
!     Compute CO2 gas fluxes on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE BCVP
      USE BCVH
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
      SUB_LOG(ISUB_LOG) = '/DFFGAE'
      I = ID(N)
      NQX = NSX(N)+1
      DXMGA = XMGA(2,N) - XMGAB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Vapor micro emulsion  ---
!
        IF( ZLAB(MN,NB).LT.-EPSL ) THEN
          FLAP = 0.D+0
          FLAB = XGAB(MN,NB)*RHOGB(MN,NB)*(-ZLAB(MN,NB))
          INDX = 3
          FLA = DIFMN( FLAP,FLAB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLA = 0.D+0
        ENDIF
        FGAP = XGA(MN,N)*RHOG(MN,N)
        FGAB = XGAB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGA = DIFMN( FGAP,FGAB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.26 ) THEN
          IF( UG(1,NQX).GT.EPSL ) THEN
            UGA(M,NQX) = UG(M,NQX)*FGA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.42 ) THEN
          IF( UG(1,NQX).LT.-EPSL ) THEN
            UGA(M,NQX) = UG(M,NQX)*FGA
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGA(MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGAB(MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 20
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMGA,INDX )
          UDGA(M,NQX) = DFM*(XMGA(MN,N)
     &      - XMGAB(MN,NB))/(5.D-1*DXGF(N))
          UGA(M,NQX) = UL(M,NQX)*FLA + UG(M,NQX)*FGA + WTMA*UDGA(M,NQX)
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
!     Compute CO2 gas fluxes on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE BCVP
      USE BCVH
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
      SUB_LOG(ISUB_LOG) = '/DFFGAN'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMGA = XMGA(2,N) - XMGAB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Vapor micro emulsion  ---
!
        IF( ZLAB(MN,NB).LT.-EPSL ) THEN
          FLAP = 0.D+0
          FLAB = XGAB(MN,NB)*RHOGB(MN,NB)*(-ZLAB(MN,NB))
          INDX = 3
          FLA = DIFMN( FLAP,FLAB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLA = 0.D+0
        ENDIF
        FGAP = XGA(MN,N)*RHOG(MN,N)
        FGAB = XGAB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGA = DIFMN( FGAP,FGAB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.26 ) THEN
          IF( VG(1,NQY).GT.EPSL ) THEN
            VGA(M,NQY) = VG(M,NQY)*FGA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.42 ) THEN
          IF( VG(1,NQY).LT.-EPSL ) THEN
            VGA(M,NQY) = VG(M,NQY)*FGA
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGA(MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGAB(MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 20
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMGA,INDX )
          VDGA(M,NQY) = DFM*(XMGA(MN,N)
     &      - XMGAB(MN,NB))/((5.D-1*DYGF(N))*RP(I))
          VGA(M,NQY) = VL(M,NQY)*FLA + VG(M,NQY)*FGA + WTMA*VDGA(M,NQY)
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
!     Compute CO2 gas fluxes on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE BCVP
      USE BCVH
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
      SUB_LOG(ISUB_LOG) = '/DFFGAS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMGA = XMGAB(2,NB) - XMGA(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Vapor micro emulsion  ---
!
        IF( ZLAB(MP,NB).LT.-EPSL ) THEN
          FLAP = 0.D+0
          FLAB = XGAB(MP,NB)*RHOGB(MP,NB)*(-ZLAB(MP,NB))
          INDX = 3
          FLA = DIFMN( FLAB,FLAP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLA = 0.D+0
        ENDIF
        FGAP = XGA(MP,N)*RHOG(MP,N)
        FGAB = XGAB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGA = DIFMN( FGAB,FGAP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.26 ) THEN
          IF( VG(1,NPY).LT.-EPSL ) THEN
            VGA(M,NPY) = VG(M,NPY)*FGA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.42 ) THEN
          IF( VG(1,NPY).GT.EPSL ) THEN
            VGA(M,NPY) = VG(M,NPY)*FGA
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGA(MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGAB(MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 20
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMGA,INDX )
          VDGA(M,NPY) = DFM*(XMGAB(MP,NB)
     &      - XMGA(MP,N))/((5.D-1*DYGF(N))*RP(I))
          VGA(M,NPY) = VL(M,NPY)*FLA + VG(M,NPY)*FGA + WTMA*VDGA(M,NPY)
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
!     Compute CO2 gas fluxes on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE BCVP
      USE BCVH
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
      SUB_LOG(ISUB_LOG) = '/DFFGAT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMGA = XMGA(2,N)-XMGAB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Vapor micro emulsion  ---
!
        IF( ZLAB(MN,NB).LT.-EPSL ) THEN
          FLAP = 0.D+0
          FLAB = XGAB(MN,NB)*RHOGB(MN,NB)*(-ZLAB(MN,NB))
          INDX = 3
          FLA = DIFMN( FLAP,FLAB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLA = 0.D+0
        ENDIF
        FGAP = XGA(MN,N)*RHOG(MN,N)
        FGAB = XGAB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGA = DIFMN( FGAP,FGAB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.26 ) THEN
          IF( WG(1,NQZ).GT.EPSL ) THEN
            WGA(M,NQZ) = WG(M,NQZ)*FGA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.42 ) THEN
          IF( WG(1,NQZ).LT.-EPSL ) THEN
            WGA(M,NQZ) = WG(M,NQZ)*FGA
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGA(MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGAB(MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 20
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMGA,INDX )
          WDGA(M,NQZ) = DFM*(XMGA(MN,N)
     &      - XMGAB(MN,NB))/(5.D-1*DZGF(N))
          WGA(M,NQZ) = WL(M,NQZ)*FLA + WG(M,NQZ)*FGA + WTMA*WDGA(M,NQZ)
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
!     Compute CO2 gas fluxes on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE BCVP
      USE BCVH
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
      SUB_LOG(ISUB_LOG) = '/DFFGAW'
      I = ID(N)
      NPX = NSX(N)
      DXMGA = XMGAB(2,NB)-XMGA(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Vapor micro emulsion  ---
!
        IF( ZLAB(MP,NB).LT.-EPSL ) THEN
          FLAP = 0.D+0
          FLAB = XGAB(MP,NB)*RHOGB(MP,NB)*(-ZLAB(MP,NB))
          INDX = 3
          FLA = DIFMN( FLAB,FLAP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLA = 0.D+0
        ENDIF
        FGAP = XGA(MP,N)*RHOG(MP,N)
        FGAB = XGAB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGA = DIFMN( FGAB,FGAP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.26 ) THEN
          IF( UG(1,NPX).LT.-EPSL ) THEN
            UGA(M,NPX) = UG(M,NPX)*FGA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.42 ) THEN
          IF( UG(1,NPX).GT.EPSL ) THEN
            UGA(M,NPX) = UG(M,NPX)*FGA
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGA(MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGAB(MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 20
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMGA,INDX )
          UDGA(M,NPX) = DFM*(XMGAB(MP,NB)
     &      - XMGA(MP,N))/(5.D-1*DXGF(N))
          UGA(M,NPX) = UL(M,NPX)*FLA + UG(M,NPX)*FGA + WTMA*UDGA(M,NPX)
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
!     Compute CH4 gas fluxes on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE CONST
      USE BCVP
      USE BCVH
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
      SUB_LOG(ISUB_LOG) = '/DFFGOB'
      K = KD(N)
      NPZ = NSZ(N)
      DXMGO = XMGOB(2,NB) - XMGO(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Vapor micro emulsion  ---
!
        IF( ZLAB(MP,NB).LT.-EPSL ) THEN
          FLOP = 0.D+0
          FLOB = XGOB(MP,NB)*RHOGB(MP,NB)*(-ZLAB(MP,NB))
          INDX = 3
          FLO = DIFMN( FLOB,FLOP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLO = 0.D+0
        ENDIF
        FGOP = XGO(MP,N)*RHOG(MP,N)
        FGOB = XGOB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGO = DIFMN( FGOB,FGOP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.26 ) THEN
          IF( WG(1,NPZ).LT.-EPSL ) THEN
            WGO(M,NPZ) = WG(M,NPZ)*FGO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.42 ) THEN
          IF( WG(1,NPZ).GT.EPSL ) THEN
            WGO(M,NPZ) = WG(M,NPZ)*FGO
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGO(MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGOB(MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 13
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMGO,INDX )
          WDGO(M,NPZ) = DFM*(XMGOB(MP,NB)
     &      - XMGO(MP,N))/(5.D-1*DZGF(N))
          WGO(M,NPZ) = WL(M,NPZ)*FLO + WG(M,NPZ)*FGO + WTMO*WDGO(M,NPZ)
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
!     Compute CH4 gas fluxes on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE CONST
      USE BCVP
      USE BCVH
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
      SUB_LOG(ISUB_LOG) = '/DFFGOE'
      I = ID(N)
      NQX = NSX(N)+1
      DXMGO = XMGO(2,N) - XMGOB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Vapor micro emulsion  ---
!
        IF( ZLAB(MN,NB).LT.-EPSL ) THEN
          FLOP = 0.D+0
          FLOB = XGOB(MN,NB)*RHOGB(MN,NB)*(-ZLAB(MN,NB))
          INDX = 3
          FLO = DIFMN( FLOP,FLOB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLO = 0.D+0
        ENDIF
        FGOP = XGO(MN,N)*RHOG(MN,N)
        FGOB = XGOB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGO = DIFMN( FGOP,FGOB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.26 ) THEN
          IF( UG(1,NQX).GT.EPSL ) THEN
            UGO(M,NQX) = UG(M,NQX)*FGO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.42 ) THEN
          IF( UG(1,NQX).LT.-EPSL ) THEN
            UGO(M,NQX) = UG(M,NQX)*FGO
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGO(MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGOB(MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 13
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMGO,INDX )
          UDGO(M,NQX) = DFM*(XMGO(MN,N)
     &      - XMGOB(MN,NB))/(5.D-1*DXGF(N))
          UGO(M,NQX) = UL(M,NQX)*FLO + UG(M,NQX)*FGO + WTMO*UDGO(M,NQX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGOE group  ---
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
!     Compute CH4 gas fluxes on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE CONST
      USE BCVP
      USE BCVH
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
      SUB_LOG(ISUB_LOG) = '/DFFGON'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMGO = XMGO(2,N) - XMGOB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Vapor micro emulsion  ---
!
        IF( ZLAB(MN,NB).LT.-EPSL ) THEN
          FLOP = 0.D+0
          FLOB = XGOB(MN,NB)*RHOGB(MN,NB)*(-ZLAB(MN,NB))
          INDX = 3
          FLO = DIFMN( FLOP,FLOB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLO = 0.D+0
        ENDIF
        FGOP = XGO(MN,N)*RHOG(MN,N)
        FGOB = XGOB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGO = DIFMN( FGOP,FGOB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.26 ) THEN
          IF( VG(1,NQY).GT.EPSL ) THEN
            VGO(M,NQY) = VG(M,NQY)*FGO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.42 ) THEN
          IF( VG(1,NQY).LT.-EPSL ) THEN
            VGO(M,NQY) = VG(M,NQY)*FGO
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGO(MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGOB(MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 13
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMGO,INDX )
          VDGO(M,NQY) = DFM*(XMGO(MN,N)
     &      - XMGOB(MN,NB))/((5.D-1*DYGF(N))*RP(I))
          VGO(M,NQY) = VL(M,NQY)*FLO + VG(M,NQY)*FGO + WTMO*VDGO(M,NQY)
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
!     Compute CH4 gas fluxes on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE CONST
      USE BCVP
      USE BCVH
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
      SUB_LOG(ISUB_LOG) = '/DFFGOS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMGO = XMGOB(2,NB) - XMGO(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Vapor micro emulsion  ---
!
        IF( ZLAB(MP,NB).LT.-EPSL ) THEN
          FLOP = 0.D+0
          FLOB = XGOB(MP,NB)*RHOGB(MP,NB)*(-ZLAB(MP,NB))
          INDX = 3
          FLO = DIFMN( FLOB,FLOP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLO = 0.D+0
        ENDIF
        FGOP = XGO(MP,N)*RHOG(MP,N)
        FGOB = XGOB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGO = DIFMN( FGOB,FGOP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.26 ) THEN
          IF( VG(1,NPY).LT.-EPSL ) THEN
            VGO(M,NPY) = VG(M,NPY)*FGO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.42 ) THEN
          IF( VG(1,NPY).GT.EPSL ) THEN
            VGO(M,NPY) = VG(M,NPY)*FGO
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGO(MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGOB(MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 13
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMGO,INDX )
          VDGO(M,NPY) = DFM*(XMGOB(MP,NB)
     &      - XMGO(MP,N))/((5.D-1*DYGF(N))*RP(I))
          VGO(M,NPY) = VL(M,NPY)*FLO + VG(M,NPY)*FGO + WTMO*VDGO(M,NPY)
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
!     Compute CH4 gas fluxes on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE CONST
      USE BCVP
      USE BCVH
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
      SUB_LOG(ISUB_LOG) = '/DFFGOT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMGO = XMGO(2,N)-XMGOB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Vapor micro emulsion  ---
!
        IF( ZLAB(MN,NB).LT.-EPSL ) THEN
          FLOP = 0.D+0
          FLOB = XGOB(MN,NB)*RHOGB(MN,NB)*(-ZLAB(MN,NB))
          INDX = 3
          FLO = DIFMN( FLOP,FLOB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLO = 0.D+0
        ENDIF
        FGOP = XGO(MN,N)*RHOG(MN,N)
        FGOB = XGOB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGO = DIFMN( FGOP,FGOB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.26 ) THEN
          IF( WG(1,NQZ).GT.EPSL ) THEN
            WGO(M,NQZ) = WG(M,NQZ)*FGO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.42 ) THEN
          IF( WG(1,NQZ).LT.-EPSL ) THEN
            WGO(M,NQZ) = WG(M,NQZ)*FGO
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGO(MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGOB(MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 13
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMGO,INDX )
          WDGO(M,NQZ) = DFM*(XMGO(MN,N)
     &      - XMGOB(MN,NB))/(5.D-1*DZGF(N))
          WGO(M,NQZ) = WL(M,NQZ)*FLO + WG(M,NQZ)*FGO + WTMO*WDGO(M,NQZ)
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
!     Compute CH4 gas fluxes on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE CONST
      USE BCVP
      USE BCVH
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
      SUB_LOG(ISUB_LOG) = '/DFFGOW'
      I = ID(N)
      NPX = NSX(N)
      DXMGO = XMGOB(2,NB)-XMGO(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Vapor micro emulsion  ---
!
        IF( ZLAB(MP,NB).LT.-EPSL ) THEN
          FLOP = 0.D+0
          FLOB = XGOB(MP,NB)*RHOGB(MP,NB)*(-ZLAB(MP,NB))
          INDX = 3
          FLO = DIFMN( FLOB,FLOP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLO = 0.D+0
        ENDIF
        FGOP = XGO(MP,N)*RHOG(MP,N)
        FGOB = XGOB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGO = DIFMN( FGOB,FGOP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.26 ) THEN
          IF( UG(1,NPX).LT.-EPSL ) THEN
            UGO(M,NPX) = UG(M,NPX)*FGO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.42 ) THEN
          IF( UG(1,NPX).GT.EPSL ) THEN
            UGO(M,NPX) = UG(M,NPX)*FGO
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGO(MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGOB(MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 13
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMGO,INDX )
          UDGO(M,NPX) = DFM*(XMGOB(MP,NB)
     &      - XMGO(MP,N))/(5.D-1*DXGF(N))
          UGO(M,NPX) = UL(M,NPX)*FLO + UG(M,NPX)*FGO + WTMO*UDGO(M,NPX)
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
!     Compute H2O gas fluxes on bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 30 October 2007.
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
      USE BCVH
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
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Vapor micro emulsion  ---
!
        IF( ZLAB(MP,NB).LT.-EPSL ) THEN
          FLWP = 0.D+0
          FLWB = XGWB(MP,NB)*RHOGB(MP,NB)*(-ZLAB(MP,NB))
          INDX = 3
          FLW = DIFMN( FLWB,FLWP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLW = 0.D+0
        ENDIF
        FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
        FGWP = XGW(MP,N)*RHOG(MP,N)
        INDX = 3
        FGW = DIFMN( FGWB,FGWP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.26 ) THEN
          IF( WG(1,NPZ).LT.-EPSL ) THEN
            WGW(M,NPZ) = WG(M,NPZ)*FGW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.42 ) THEN
          IF( WG(1,NPZ).GT.EPSL ) THEN
            WGW(M,NPZ) = WG(M,NPZ)*FGW
          ENDIF
        ELSE
          WDGW(M,NPZ) = -(WDGA(M,NPZ)+WDGO(M,NPZ))
          WGW(M,NPZ) = WL(M,NPZ)*FLW + WG(M,NPZ)*FGW + WTMW*WDGW(M,NPZ)
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
!     Compute H2O gas fluxes on east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 30 October 2007.
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
      USE BCVH
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
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Vapor micro emulsion  ---
!
        IF( ZLAB(MN,NB).LT.-EPSL ) THEN
          FLWP = 0.D+0
          FLWB = XGWB(MN,NB)*RHOGB(MN,NB)*(-ZLAB(MN,NB))
          INDX = 3
          FLW = DIFMN( FLWP,FLWB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLW = 0.D+0
        ENDIF
        FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
        FGWP = XGW(MN,N)*RHOG(MN,N)
        INDX = 3
        FGW = DIFMN( FGWP,FGWB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.26 ) THEN
          IF( UG(1,NQX).GT.EPSL ) THEN
            UGW(M,NQX) = UG(M,NQX)*FGW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.42 ) THEN
          IF( UG(1,NQX).LT.-EPSL ) THEN
            UGW(M,NQX) = UG(M,NQX)*FGW
          ENDIF
        ELSE
          UDGW(M,NQX) = -(UDGA(M,NQX)+UDGO(M,NQX))
          UGW(M,NQX) = UL(M,NQX)*FLW + UG(M,NQX)*FGW + WTMW*UDGW(M,NQX)
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
!     Compute H2O gas fluxes on north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 30 October 2007.
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
      USE BCVH
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
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Vapor micro emulsion  ---
!
        IF( ZLAB(MN,NB).LT.-EPSL ) THEN
          FLWP = 0.D+0
          FLWB = XGWB(MN,NB)*RHOGB(MN,NB)*(-ZLAB(MN,NB))
          INDX = 3
          FLW = DIFMN( FLWP,FLWB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLW = 0.D+0
        ENDIF
        FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
        FGWP = XGW(MN,N)*RHOG(MN,N)
        INDX = 3
        FGW = DIFMN( FGWP,FGWB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.26 ) THEN
          IF( VG(1,NQY).GT.EPSL ) THEN
            VGW(M,NQY) = VG(M,NQY)*FGW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.42 ) THEN
          IF( VG(1,NQY).LT.-EPSL ) THEN
            VGW(M,NQY) = VG(M,NQY)*FGW
          ENDIF
        ELSE
          VDGW(M,NQY) = -(VDGA(M,NQY)+VDGO(M,NQY))
          VGW(M,NQY) = VL(M,NQY)*FLW + VG(M,NQY)*FGW + WTMW*VDGW(M,NQY)
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
!     Compute H2O gas fluxes on south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 30 October 2007.
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
      USE BCVH
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
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Vapor micro emulsion  ---
!
        IF( ZLAB(MP,NB).LT.-EPSL ) THEN
          FLWP = 0.D+0
          FLWB = XGWB(MP,NB)*RHOGB(MP,NB)*(-ZLAB(MP,NB))
          INDX = 3
          FLW = DIFMN( FLWB,FLWP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLW = 0.D+0
        ENDIF
        FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
        FGWP = XGW(MP,N)*RHOG(MP,N)
        INDX = 3
        FGW = DIFMN( FGWB,FGWP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.26 ) THEN
          IF( VG(1,NPY).LT.-EPSL ) THEN
            VGW(M,NPY) = VG(M,NPY)*FGW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.42 ) THEN
          IF( VG(1,NPY).GT.EPSL ) THEN
            VGW(M,NPY) = VG(M,NPY)*FGW
          ENDIF
        ELSE
          VDGW(M,NPY) = -(VDGA(M,NPY)+VDGO(M,NPY))
          VGW(M,NPY) = VL(M,NPY)*FLW + VG(M,NPY)*FGW + WTMW*VDGW(M,NPY)
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
!     Compute H2O gas fluxes on top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 30 October 2007.
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
      USE BCVH
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
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Vapor micro emulsion  ---
!
        IF( ZLAB(MN,NB).LT.-EPSL ) THEN
          FLWP = 0.D+0
          FLWB = XGWB(MN,NB)*RHOGB(MN,NB)*(-ZLAB(MN,NB))
          INDX = 3
          FLW = DIFMN( FLWP,FLWB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLW = 0.D+0
        ENDIF
        FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
        FGWP = XGW(MN,N)*RHOG(MN,N)
        INDX = 3
        FGW = DIFMN( FGWP,FGWB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.26 ) THEN
          IF( WG(1,NQZ).GT.EPSL ) THEN
            WGW(M,NQZ) = WG(M,NQZ)*FGW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.42 ) THEN
          IF( WG(1,NQZ).LT.-EPSL ) THEN
            WGW(M,NQZ) = WG(M,NQZ)*FGW
          ENDIF
        ELSE
          WDGW(M,NQZ) = -(WDGA(M,NQZ)+WDGO(M,NQZ))
          WGW(M,NQZ) = WL(M,NQZ)*FLW + WG(M,NQZ)*FGW + WTMW*WDGW(M,NQZ)
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
!     Compute H2O gas fluxes on west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 30 October 2007.
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
      USE BCVH
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
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Vapor micro emulsion  ---
!
        IF( ZLAB(MP,NB).LT.-EPSL ) THEN
          FLWP = 0.D+0
          FLWB = XLWB(MP,NB)*RHOLB(MP,NB)*(-ZLAB(MP,NB))
          INDX = 3
          FLW = DIFMN( FLWB,FLWP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLW = 0.D+0
        ENDIF
        FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
        FGWP = XGW(MP,N)*RHOG(MP,N)
        INDX = 3
        FGW = DIFMN( FGWB,FGWP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(3,NB).EQ.26 ) THEN
          IF( UG(1,NPX).LT.-EPSL ) THEN
            UGW(M,NPX) = UG(M,NPX)*FGW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(3,NB).EQ.42 ) THEN
          IF( UG(1,NPX).GT.EPSL ) THEN
            UGW(M,NPX) = UG(M,NPX)*FGW
          ENDIF
        ELSE
          UDGW(M,NPX) = -(UDGA(M,NPX)+UDGO(M,NPX))
          UGW(M,NPX) = UL(M,NPX)*FLW + UG(M,NPX)*FGW + WTMW*UDGW(M,NPX)
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
!     Compute CO2 aqeuous fluxes on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE BCVH
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
        ZLABX = (1.D+0-ABS(ZLAB(MP,NB)))
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MP,NB).GT.EPSL ) THEN
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)*ZLABX
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MP,NB).LT.-EPSL ) THEN
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)*ZLABX
!
!---    No micro emulsion  ---
!
        ELSE
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)
        ENDIF
        INDX = 3
        FLA = DIFMN( FLAB,FLAP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.26 ) THEN
          IF( WL(1,NPZ).LT.-EPSL ) THEN
            WLA(M,NPZ) = WL(M,NPZ)*FLA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.42 ) THEN
          IF( WL(1,NPZ).GT.EPSL ) THEN
            WLA(M,NPZ) = WL(M,NPZ)*FLA
          ENDIF
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)*
     &      ZLABX*RHOMLB(MP,NB)
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
!     Compute CO2 aqeuous fluxes on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE BCVH
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
        ZLABX = (1.D+0-ABS(ZLAB(MN,NB)))
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MN,NB).GT.EPSL ) THEN
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)*ZLABX
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MN,NB).LT.-EPSL ) THEN
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)*ZLABX
!
!---    No micro emulsion  ---
!
        ELSE
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)
        ENDIF
        INDX = 3
        FLA = DIFMN( FLAP,FLAB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.26 ) THEN
          IF( UL(1,NQX).GT.EPSL ) THEN
            ULA(M,NQX) = UL(M,NQX)*FLA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.42 ) THEN
          IF( UL(1,NQX).LT.-EPSL ) THEN
            ULA(M,NQX) = UL(M,NQX)*FLA
          ENDIF
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)*
     &      ZLABX*RHOMLB(MN,NB)
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
!     Compute CO2 aqeuous fluxes on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE BCVH
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
        ZLABX = (1.D+0-ABS(ZLAB(MN,NB)))
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MN,NB).GT.EPSL ) THEN
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)*ZLABX
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MN,NB).LT.-EPSL ) THEN
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)*ZLABX
!
!---    No micro emulsion  ---
!
        ELSE
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)
        ENDIF
        INDX = 3
        FLA = DIFMN( FLAP,FLAB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.26 ) THEN
          IF( VL(1,NQY).GT.EPSL ) THEN
            VLA(M,NQY) = VL(M,NQY)*FLA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.42 ) THEN
          IF( VL(1,NQY).LT.-EPSL ) THEN
            VLA(M,NQY) = VL(M,NQY)*FLA
          ENDIF
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)*
     &      ZLABX*RHOMLB(MN,NB)
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
!     Compute CO2 aqeuous fluxes on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE BCVH
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
        ZLABX = (1.D+0-ABS(ZLAB(MP,NB)))
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MP,NB).GT.EPSL ) THEN
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)*ZLABX
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MP,NB).LT.-EPSL ) THEN
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)*ZLABX
!
!---    No micro emulsion  ---
!
        ELSE
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)
        ENDIF
        INDX = 3
        FLA = DIFMN( FLAB,FLAP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.26 ) THEN
          IF( VL(1,NPY).LT.-EPSL ) THEN
            VLA(M,NPY) = VL(M,NPY)*FLA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.42 ) THEN
          IF( VL(1,NPY).GT.EPSL ) THEN
            VLA(M,NPY) = VL(M,NPY)*FLA
          ENDIF
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)*
     &      ZLABX*RHOMLB(MP,NB)
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
!     Compute CO2 aqeuous fluxes on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE BCVH
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
        ZLABX = (1.D+0-ABS(ZLAB(MN,NB)))
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MN,NB).GT.EPSL ) THEN
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)*ZLABX
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MN,NB).LT.-EPSL ) THEN
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)*ZLABX
!
!---    No micro emulsion  ---
!
        ELSE
          FLAP = XLA(MN,N)*RHOL(MN,N)
          FLAB = XLAB(MN,NB)*RHOLB(MN,NB)
        ENDIF
        INDX = 3
        FLA = DIFMN( FLAP,FLAB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.26 ) THEN
          IF( WL(1,NQZ).GT.EPSL ) THEN
            WLA(M,NQZ) = WL(M,NQZ)*FLA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.42 ) THEN
          IF( WL(1,NQZ).LT.-EPSL ) THEN
            WLA(M,NQZ) = WL(M,NQZ)*FLA
          ENDIF
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)*
     &      ZLABX*RHOMLB(MN,NB)
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
!     Compute CO2 aqeuous fluxes on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE BCVH
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
        ZLABX = (1.D+0-ABS(ZLAB(MP,NB)))
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MP,NB).GT.EPSL ) THEN
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)*ZLABX
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MP,NB).LT.-EPSL ) THEN
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)*ZLABX
!
!---    No micro emulsion  ---
!
        ELSE
          FLAP = XLA(MP,N)*RHOL(MP,N)
          FLAB = XLAB(MP,NB)*RHOLB(MP,NB)
        ENDIF
        INDX = 3
        FLA = DIFMN( FLAB,FLAP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.26 ) THEN
          IF( UL(1,NPX).LT.-EPSL ) THEN
            ULA(M,NPX) = UL(M,NPX)*FLA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.42 ) THEN
          IF( UL(1,NPX).GT.EPSL ) THEN
            ULA(M,NPX) = UL(M,NPX)*FLA
          ENDIF
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)*
     &      ZLABX*RHOMLB(MP,NB)
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
!     Compute CH4 aqeuous fluxes on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE BCVH
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
        ZLABX = (1.D+0-ABS(ZLAB(MP,NB)))
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MP,NB).GT.EPSL ) THEN
          FLOP = XLO(MP,N)*RHOL(MP,N)
          FLOB = XLOB(MP,NB)*RHOLB(MP,NB)*ZLABX
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MP,NB).LT.-EPSL ) THEN
          FLOP = XLO(MP,N)*RHOL(MP,N)
          FLOB = XLOB(MP,NB)*RHOLB(MP,NB)*ZLABX
!
!---    No micro emulsion  ---
!
        ELSE
          FLOP = XLO(MP,N)*RHOL(MP,N)
          FLOB = XLOB(MP,NB)*RHOLB(MP,NB)
        ENDIF
        INDX = 3
        FLO = DIFMN( FLOB,FLOP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.26 ) THEN
          IF( WL(1,NPZ).LT.-EPSL ) THEN
            WLO(M,NPZ) = WL(M,NPZ)*FLO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.42 ) THEN
          IF( WL(1,NPZ).GT.EPSL ) THEN
            WLO(M,NPZ) = WL(M,NPZ)*FLO
          ENDIF
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLO(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLOB(MP,NB)*
     &      ZLABX*RHOMLB(MP,NB)
          INDX = 15
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
!     Compute CH4 aqeuous fluxes on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE BCVH
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
        ZLABX = (1.D+0-ABS(ZLAB(MN,NB)))
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MN,NB).GT.EPSL ) THEN
          FLOP = XLO(MN,N)*RHOL(MN,N)
          FLOB = XLOB(MN,NB)*RHOLB(MN,NB)*ZLABX
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MN,NB).LT.-EPSL ) THEN
          FLOP = XLO(MN,N)*RHOL(MN,N)
          FLOB = XLOB(MN,NB)*RHOLB(MN,NB)*ZLABX
!
!---    No micro emulsion  ---
!
        ELSE
          FLOP = XLO(MN,N)*RHOL(MN,N)
          FLOB = XLOB(MN,NB)*RHOLB(MN,NB)
        ENDIF
        INDX = 3
        FLO = DIFMN( FLOP,FLOB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.26 ) THEN
          IF( UL(1,NQX).GT.EPSL ) THEN
            ULO(M,NQX) = UL(M,NQX)*FLO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.42 ) THEN
          IF( UL(1,NQX).LT.-EPSL ) THEN
            ULO(M,NQX) = UL(M,NQX)*FLO
          ENDIF
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLO(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLOB(MN,NB)*
     &      ZLABX*RHOMLB(MN,NB)
          INDX = 15
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
!     Compute CH4 aqeuous fluxes on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE BCVH
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
        ZLABX = (1.D+0-ABS(ZLAB(MN,NB)))
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MN,NB).GT.EPSL ) THEN
          FLOP = XLO(MN,N)*RHOL(MN,N)
          FLOB = XLOB(MN,NB)*RHOLB(MN,NB)*ZLABX
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MN,NB).LT.-EPSL ) THEN
          FLOP = XLO(MN,N)*RHOL(MN,N)
          FLOB = XLOB(MN,NB)*RHOLB(MN,NB)*ZLABX
!
!---    No micro emulsion  ---
!
        ELSE
          FLOP = XLO(MN,N)*RHOL(MN,N)
          FLOB = XLOB(MN,NB)*RHOLB(MN,NB)
        ENDIF
        INDX = 3
        FLO = DIFMN( FLOP,FLOB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.26 ) THEN
          IF( VL(1,NQY).GT.EPSL ) THEN
            VLO(M,NQY) = VL(M,NQY)*FLO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.42 ) THEN
          IF( VL(1,NQY).LT.-EPSL ) THEN
            VLO(M,NQY) = VL(M,NQY)*FLO
          ENDIF
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLO(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLOB(MN,NB)*
     &      ZLABX*RHOMLB(MN,NB)
          INDX = 15
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
!     Compute VOC vapor mole diffusion rates on a south boundary.
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
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE CONST
      USE BCVP
      USE BCVH
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
        ZLABX = (1.D+0-ABS(ZLAB(MP,NB)))
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MP,NB).GT.EPSL ) THEN
          FLOP = XLO(MP,N)*RHOL(MP,N)
          FLOB = XLOB(MP,NB)*RHOLB(MP,NB)*ZLABX
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MP,NB).LT.-EPSL ) THEN
          FLOP = XLO(MP,N)*RHOL(MP,N)
          FLOB = XLOB(MP,NB)*RHOLB(MP,NB)*ZLABX
!
!---    No micro emulsion  ---
!
        ELSE
          FLOP = XLO(MP,N)*RHOL(MP,N)
          FLOB = XLOB(MP,NB)*RHOLB(MP,NB)
        ENDIF
        INDX = 3
        FLO = DIFMN( FLOB,FLOP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.26 ) THEN
          IF( VL(1,NPY).LT.-EPSL ) THEN
            VLO(M,NPY) = VL(M,NPY)*FLO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.42 ) THEN
          IF( VL(1,NPY).GT.EPSL ) THEN
            VLO(M,NPY) = VL(M,NPY)*FLO
          ENDIF
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLO(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLOB(MP,NB)*
     &      ZLABX*RHOMLB(MP,NB)
          INDX = 15
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
!     Compute VOC vapor mole diffusion rates on a top boundary.
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
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE CONST
      USE BCVP
      USE BCVH
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
        ZLABX = (1.D+0-ABS(ZLAB(MN,NB)))
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MN,NB).GT.EPSL ) THEN
          FLOP = XLO(MN,N)*RHOL(MN,N)
          FLOB = XLOB(MN,NB)*RHOLB(MN,NB)*ZLABX
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MN,NB).LT.-EPSL ) THEN
          FLOP = XLO(MN,N)*RHOL(MN,N)
          FLOB = XLOB(MN,NB)*RHOLB(MN,NB)*ZLABX
!
!---    No micro emulsion  ---
!
        ELSE
          FLOP = XLO(MN,N)*RHOL(MN,N)
          FLOB = XLOB(MN,NB)*RHOLB(MN,NB)
        ENDIF
        INDX = 3
        FLO = DIFMN( FLOP,FLOB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.26 ) THEN
          IF( WL(1,NQZ).GT.EPSL ) THEN
            WLO(M,NQZ) = WL(M,NQZ)*FLO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.42 ) THEN
          IF( WL(1,NQZ).LT.-EPSL ) THEN
            WLO(M,NQZ) = WL(M,NQZ)*FLO
          ENDIF
        ELSE
          DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLO(MN,N)
     &      *RHOML(MN,N)
          DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLOB(MN,NB)*
     &      ZLABX*RHOMLB(MN,NB)
          INDX = 15
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMLO,INDX )
          WDLO(M,NQZ) = DFM*(XMLO(MN,N)
     &     - XMLOB(MN,NB))/(5.D-1*DZGF(N))
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
!     Compute VOC vapor mole diffusion rates on a west boundary.
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
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXD
      USE FDVP
      USE CONST
      USE BCVP
      USE BCVH
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
        ZLABX = (1.D+0-ABS(ZLAB(MP,NB)))
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MP,NB).GT.EPSL ) THEN
          FLOP = XLO(MP,N)*RHOL(MP,N)
          FLOB = XLOB(MP,NB)*RHOLB(MP,NB)*ZLABX
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MP,NB).LT.-EPSL ) THEN
          FLOP = XLO(MP,N)*RHOL(MP,N)
          FLOB = XLOB(MP,NB)*RHOLB(MP,NB)*ZLABX
!
!---    No micro emulsion  ---
!
        ELSE
          FLOP = XLO(MP,N)*RHOL(MP,N)
          FLOB = XLOB(MP,NB)*RHOLB(MP,NB)
        ENDIF
        INDX = 3
        FLO = DIFMN( FLOB,FLOP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.26 ) THEN
          IF( UL(1,NPX).LT.-EPSL ) THEN
            ULO(M,NPX) = UL(M,NPX)*FLO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.42 ) THEN
          IF( UL(1,NPX).GT.EPSL ) THEN
            ULO(M,NPX) = UL(M,NPX)*FLO
          ENDIF
        ELSE
          DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLO(MP,N)
     &      *RHOML(MP,N)
          DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLOB(MP,NB)*
     &      ZLABX*RHOMLB(MP,NB)
          INDX = 15
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
      SUBROUTINE DFFLSB( N,NB )
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
!     Compute salt aqueous-phase fluxes on bottom boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_hyd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
      USE FLUXP
      USE FDVS
      USE FDVP
      USE CONST
      USE BCVS
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLSB'
      K = KD(N)
      IZN = IZ(N)
      NPZ = NSZ(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---  Diffusion coefficients  ---
!
        IF( IEDLS.EQ.1 ) THEN
          TCOR = (T(MP,N)+TABS)/TSPRF
          SMDLP = DFLS(MP,N)*TCOR*(VISRL/VISL(MP,N))
          DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
          TCOR = (TB(MP,NB)+TABS)/TSPRF
          SMDLB = DFLSB(MP,NB)*TCOR*(VISRL/VISLB(MP,NB))
          DFFLB = TORLB(MP,NB)*SLB(MP,NB)*PORDB(MP,NB)*SMDLB
        ELSEIF( IEDLS.EQ.2 ) THEN
          DFFLP = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SL(MP,N)*
     &      PORD(MP,N)*SDCLS(3,IZN))
          DFFLB = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SLB(MP,NB)*
     &      PORDB(MP,NB)*SDCLS(3,IZN))
        ELSEIF( IEDLS.EQ.3 ) THEN
          DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
          DFFLB = TORLB(MP,NB)*SLB(MP,NB)*PORDB(MP,NB)*DFLSB(MP,NB)
        ENDIF
        INDX = 18
        DFFLB = DIFMN(DFFLB,DFFLP,DZGF(N),DZGF(N),WL(1,NPZ),INDX)
!
!---   Outflow boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.7 ) THEN
          ALP = MAX( -WL(M,NPZ),ZERO )
          WS(M,NPZ) = -XLS(MP,N)*RHOL(MP,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.35 .OR. IBCT(IEQS,NB).EQ.37 ) THEN
          AL = MAX( WL(M,NPZ),ZERO )
          WS(M,NPZ) = XLSB(MP,NB)*RHOLB(MP,NB)*AL
!
!---   Dirichlet boundary types  ---
!
        ELSE
          DDLB = DFFLB/(5.D-1*DZGF(N))
          AL = MAX( WL(M,NPZ),ZERO ) +
     &     DDLB*MAX((ONE-(TENTH*ABS(WL(M,NPZ))/(DDLB+SMALL)))**5,ZERO)
          ALP = MAX( -WL(M,NPZ),ZERO ) +
     &     DDLB*MAX((ONE-(TENTH*ABS(WL(M,NPZ))/(DDLB+SMALL)))**5,ZERO)
          WS(M,NPZ) = (XLSB(MP,NB)*RHOLB(MP,NB)*AL
     &     - XLS(MP,N)*RHOL(MP,N)*ALP)
          WDS(M,NPZ) = DDLB*(XLSB(MP,NB)*RHOLB(MP,NB) -
     &      XLS(MP,N)*RHOL(MP,N))
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLSB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSS( N,NB )
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
!     Compute salt aqueous-phase fluxes on south boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_hyd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
      USE FLUXP
      USE FDVS
      USE FDVP
      USE CONST
      USE BCVS
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLSS'
      J = JD(N)
      I = ID(N)
      IZN = IZ(N)
      NPY = NSY(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---  Diffusion coefficients  ---
!
        IF( IEDLS.EQ.1 ) THEN
          TCOR = (T(MP,N)+TABS)/TSPRF
          SMDLP = DFLS(MP,N)*TCOR*(VISRL/VISL(MP,N))
          DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
          TCOR = (TB(MP,NB)+TABS)/TSPRF
          SMDLP = DFLSB(MP,NB)*TCOR*(VISRL/VISLB(MP,NB))
          DFFLS = TORLB(MP,NB)*SLB(MP,NB)*PORDB(MP,NB)*SMDLP
        ELSEIF( IEDLS.EQ.2 ) THEN
          DFFLP = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SL(MP,N)*
     &      PORD(MP,N)*SDCLS(3,IZN))
          DFFLS = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SLB(MP,NB)*
     &      PORDB(MP,NB)*SDCLS(3,IZN))
        ELSEIF( IEDLS.EQ.3 ) THEN
          DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
          DFFLS = TORLB(MP,NB)*SLB(MP,NB)*PORDB(MP,NB)*DFLSB(MP,NB)
        ENDIF
        INDX = 18
        DFFLS = DIFMN(DFFLS,DFFLP,DYGF(N),DYGF(N),VL(1,NPY),INDX)
!
!---   Outflow boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.7 ) THEN
          ALP = MAX( -VL(M,NPY),ZERO )
          VS(M,NPY) = -XLS(MP,N)*RHOL(MP,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.35 .OR. IBCT(IEQS,NB).EQ.37 ) THEN
          AL = MAX( VL(M,NPY),ZERO )
          VS(M,NPY) = XLSB(MP,NB)*RHOLB(MP,NB)*AL
!
!---   Dirichlet boundary types  ---
!
        ELSE
          DDLS = DFFLS/RP(I)/(5.D-1*DYGF(N))
          AL = MAX( VL(M,NPY),ZERO ) +
     &     DDLS*MAX((ONE-(TENTH*ABS(VL(M,NPY))/(DDLS+SMALL)))**5,ZERO)
          ALP = MAX( -VL(M,NPY),ZERO ) +
     &     DDLS*MAX((ONE-(TENTH*ABS(VL(M,NPY))/(DDLS+SMALL)))**5,ZERO)
          VS(M,NPY) = (XLSB(MP,NB)*RHOLB(MP,NB)*AL
     &      - XLS(MP,N)*RHOL(MP,N)*ALP)
          VDS(M,NPY) = DDLS*(XLSB(MP,NB)*RHOLB(MP,NB) -
     &      XLS(MP,N)*RHOL(MP,N))
        ENDIF

  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLSS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSW( N,NB )
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
!     Compute salt aqueous-phase fluxes on west boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_hyd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
      USE FLUXP
      USE FDVS
      USE FDVP
      USE CONST
      USE BCVS
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLSW'
      I = ID(N)
      IZN = IZ(N)
      NPX = NSX(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---  Diffusion coefficients  ---
!
        IF( IEDLS.EQ.1 ) THEN
          TCOR = (T(MP,N)+TABS)/TSPRF
          SMDLP = DFLS(MP,N)*TCOR*(VISRL/VISL(MP,N))
          DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*SMDLP
          TCOR = (TB(MP,NB)+TABS)/TSPRF
          SMDLP = DFLSB(MP,NB)*TCOR*(VISRL/VISLB(MP,NB))
          DFFLW = TORLB(MP,NB)*SLB(MP,NB)*PORDB(MP,NB)*SMDLP
        ELSEIF( IEDLS.EQ.2 ) THEN
          DFFLP = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SL(MP,N)*
     &      PORD(MP,N)*SDCLS(3,IZN))
          DFFLW = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SLB(MP,NB)*
     &      PORDB(MP,NB)*SDCLS(3,IZN))
        ELSEIF( IEDLS.EQ.3 ) THEN
          DFFLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
          DFFLW = TORLB(MP,NB)*SLB(MP,NB)*PORDB(MP,NB)*DFLSB(MP,NB)
        ENDIF
        INDX = 18
        DFFLW = DIFMN(DFFLW,DFFLP,DXGF(N),DXGF(N),UL(1,NPX),INDX)
!
!---   Outflow boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.7 ) THEN
          ALP = MAX( -UL(M,NPX),ZERO )
          US(M,NPX) = -XLS(MP,N)*RHOL(MP,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.35 .OR. IBCT(IEQS,NB).EQ.37 ) THEN
          AL = MAX( UL(M,NPX),ZERO )
          US(M,NPX) = XLSB(MP,NB)*RHOLB(MP,NB)*AL
!
!---   Dirichlet boundary types  ---
!
        ELSE
          DDLW = DFFLW/(5.D-1*DXGF(N))
          AL = MAX( UL(M,NPX),ZERO ) +
     &     DDLW*MAX((ONE-(TENTH*ABS(UL(M,NPX))/(DDLW+SMALL)))**5,ZERO)
          ALP = MAX( -UL(M,NPX),ZERO ) +
     &     DDLW*MAX((ONE-(TENTH*ABS(UL(M,NPX))/(DDLW+SMALL)))**5,ZERO)
          US(M,NPX) = (XLSB(MP,NB)*RHOLB(MP,NB)*AL
     &      - XLS(MP,N)*RHOL(MP,N)*ALP)
          UDS(M,NPX) = DDLW*(XLSB(MP,NB)*RHOLB(MP,NB) -
     &      XLS(MP,N)*RHOL(MP,N))
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLSW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSE( N,NB )
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
!     Compute salt aqueous-phase fluxes on west boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_hyd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
      USE FLUXP
      USE FDVS
      USE FDVP
      USE CONST
      USE BCVS
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLSE'
      I = ID(N)
      IZN = IZ(N)
      NQX = NSX(N)+1
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---  Diffusion coefficients  ---
!
        IF( IEDLS.EQ.1 ) THEN
          TCOR = (T(MN,N)+TABS)/TSPRF
          SMDLP = DFLS(MN,N)*TCOR*(VISRL/VISL(MN,N))
          DFFLP = TORL(MN,N)*SL(MN,N)*PORD(MN,N)*SMDLP
          TCOR = (TB(MN,NB)+TABS)/TSPRF
          SMDLP = DFLSB(MN,NB)*TCOR*(VISRL/VISLB(MN,NB))
          DFFLE = TORLB(MN,NB)*SLB(MN,NB)*PORDB(MN,NB)*SMDLP
        ELSEIF( IEDLS.EQ.2 ) THEN
          DFFLP = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SL(MN,N)*
     &      PORD(MN,N)*SDCLS(3,IZN))
          DFFLE = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SLB(MN,NB)*
     &      PORDB(MN,NB)*SDCLS(3,IZN))
        ELSEIF( IEDLS.EQ.3 ) THEN
          DFFLP = TORL(MN,N)*SL(MN,N)*PORD(MN,N)*DFLS(MN,N)
          DFFLE = TORLB(MN,NB)*SLB(MN,NB)*PORDB(MN,NB)*DFLSB(MN,NB)
        ENDIF
        INDX = 18
        DFFLE = DIFMN(DFFLP,DFFLE,DXGF(N),DXGF(N),UL(1,NQX),INDX)
!
!---   Outflow boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.7 ) THEN
          ALP = MAX( UL(M,NQX),ZERO )
          US(M,NQX) = XLS(MN,N)*RHOL(MN,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.35 .OR. IBCT(IEQS,NB).EQ.37 ) THEN
          AL = MAX( -UL(M,NQX),ZERO )
          US(M,NQX) = -XLSB(MN,NB)*RHOLB(MN,NB)*AL
!
!---   Dirichlet boundary types  ---
!
        ELSE
          DDLE = DFFLE/(5.D-1*DXGF(N))
          AL = MAX( -UL(M,NQX),ZERO ) +
     &     DDLE*MAX((ONE-(TENTH*ABS(UL(M,NQX))/(DDLE+SMALL)))**5,ZERO)
          ALP = MAX( UL(M,NQX),ZERO ) +
     &     DDLE*MAX((ONE-(TENTH*ABS(UL(M,NQX))/(DDLE+SMALL)))**5,ZERO)
          US(M,NQX) = (XLS(MN,N)*RHOL(MN,N)*ALP
     &     - XLSB(MN,NB)*RHOLB(MN,NB)*AL)
          UDS(M,NQX) = DDLE*(XLS(MN,N)*RHOL(MN,N) -
     &     XLSB(MN,NB)*RHOLB(MN,NB))
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLSE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLSN( N,NB )
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
!     Compute salt aqueous-phase fluxes on north boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_hyd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
      USE FLUXP
      USE FDVS
      USE FDVP
      USE CONST
      USE BCVS
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLSN'
      J = JD(N)
      I = ID(N)
      IZN = IZ(N)
      NQY = NSY(N)+IFLD
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---  Diffusion coefficients  ---
!
        IF( IEDLS.EQ.1 ) THEN
          TCOR = (T(MN,N)+TABS)/TSPRF
          SMDLP = DFLS(MN,N)*TCOR*(VISRL/VISL(MN,N))
          DFFLP = TORL(MN,N)*SL(MN,N)*PORD(MN,N)*SMDLP
          TCOR = (TB(MN,NB)+TABS)/TSPRF
          SMDLP = DFLSB(MN,NB)*TCOR*(VISRL/VISLB(MN,NB))
          DFFLN = TORLB(MN,NB)*SLB(MN,NB)*PORDB(MN,NB)*SMDLP
        ELSEIF( IEDLS.EQ.2 ) THEN
          DFFLP = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SL(MN,N)*
     &      PORD(MN,N)*SDCLS(3,IZN))
          DFFLN = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SLB(MN,NB)*
     &      PORDB(MN,NB)*SDCLS(3,IZN))
        ELSEIF( IEDLS.EQ.3 ) THEN
          DFFLP = TORL(MN,N)*SL(MN,N)*PORD(MN,N)*DFLS(MN,N)
          DFFLN = TORLB(MN,NB)*SLB(MN,NB)*PORDB(MN,NB)*DFLSB(MN,NB)
        ENDIF
        INDX = 18
        DFFLN = DIFMN(DFFLP,DFFLN,DYGF(N),DYGF(N),VL(1,NQY),INDX)
!
!---   Outflow boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.7 ) THEN
          ALP = MAX( VL(M,NQY),ZERO )
          VS(M,NQY) = XLS(MN,N)*RHOL(MN,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.35 .OR. IBCT(IEQS,NB).EQ.37 ) THEN
          AL = MAX( -VL(M,NQY),ZERO )
          VS(M,NQY) = -XLSB(MN,NB)*RHOLB(MN,NB)*AL
!
!---   Dirichlet boundary types  ---
!
        ELSE
          DDLN = DFFLN/RP(I)/(5.D-1*DYGF(N))
          AL = MAX( -VL(M,NQY),ZERO ) +
     &     DDLN*MAX((ONE-(TENTH*ABS(VL(M,NQY))/(DDLN+SMALL)))**5,ZERO)
          ALP = MAX( VL(M,NQY),ZERO ) +
     &     DDLN*MAX((ONE-(TENTH*ABS(VL(M,NQY))/(DDLN+SMALL)))**5,ZERO)
          VS(M,NQY) = (XLS(MN,N)*RHOL(MN,N)*ALP -
     &      XLSB(MN,NB)*RHOLB(MN,NB)*AL)
          VDS(M,NQY) = DDLN*(XLS(MN,N)*RHOL(MN,N) -
     &      XLSB(MN,NB)*RHOLB(MN,NB))
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLSN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLST( N,NB )
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
!     Compute salt aqueous-phase fluxes on top boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_hyd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXS
      USE FLUXP
      USE FDVS
      USE FDVP
      USE CONST
      USE BCVS
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLST'
      K = KD(N)
      IZN = IZ(N)
      NQZ = NSZ(N)+IJFLD
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---  Diffusion coefficients  ---
!
        IF( IEDLS.EQ.1 ) THEN
          TCOR = (T(MN,N)+TABS)/TSPRF
          SMDLP = DFLS(MN,N)*TCOR*(VISRL/VISL(MN,N))
          DFFLP = TORL(MN,N)*SL(MN,N)*PORD(MN,N)*SMDLP
          TCOR = (TB(MN,NB)+TABS)/TSPRF
          SMDLP = DFLSB(MN,NB)*TCOR*(VISRL/VISLB(MN,NB))
          DFFLT = TORLB(MN,NB)*SLB(MN,NB)*PORDB(MN,NB)*SMDLP
        ELSEIF( IEDLS.EQ.2 ) THEN
          DFFLP = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SL(MN,N)*
     &      PORD(MN,N)*SDCLS(3,IZN))
          DFFLT = SDCLS(1,IZN)*SDCLS(2,IZN)*EXP(SLB(MN,NB)*
     &      PORDB(MN,NB)*SDCLS(3,IZN))
        ELSEIF( IEDLS.EQ.3 ) THEN
          DFFLP = TORL(MN,N)*SL(MN,N)*PORD(MN,N)*DFLS(MN,N)
          DFFLT = TORLB(MN,NB)*SLB(MN,NB)*PORDB(MN,NB)*DFLSB(MN,NB)
        ENDIF
        INDX = 18
        DFFLT = DIFMN(DFFLP,DFFLT,DZGF(N),DZGF(N),WL(1,NQZ),INDX)
!
!---   Outflow boundary types  ---
!
        IF( IBCT(IEQS,NB).EQ.7 ) THEN
          ALP = MAX( WL(M,NQZ),ZERO )
          WS(M,NQZ) = XLS(MN,N)*RHOL(MN,N)*ALP
!
!---   Inflow boundary types  ---
!
        ELSEIF( IBCT(IEQS,NB).EQ.35 .OR. IBCT(IEQS,NB).EQ.37 ) THEN
          AL = MAX( -WL(M,NQZ),ZERO )
          WS(M,NQZ) = -XLSB(MN,NB)*RHOLB(MN,NB)*AL
!
!---   Dirichlet boundary types  ---
!
        ELSE
          DDLT = DFFLT/(5.D-1*DZGF(N))
          AL = MAX( -WL(M,NQZ),ZERO ) +
     &     DDLT*MAX((ONE-(TENTH*ABS(WL(M,NQZ))/(DDLT+SMALL)))**5,ZERO)
          ALP = MAX( WL(M,NQZ),ZERO ) +
     &     DDLT*MAX((ONE-(TENTH*ABS(WL(M,NQZ))/(DDLT+SMALL)))**5,ZERO)
          WS(M,NQZ) = (XLS(MN,N)*RHOL(MN,N)*ALP
     &      - XLSB(MN,NB)*RHOLB(MN,NB)*AL)
          WDS(M,NQZ) = DDLT*(XLS(MN,N)*RHOL(MN,N) -
     &      XLSB(MN,NB)*RHOLB(MN,NB))
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLST group  ---
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
!     Aqueous water flux on bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_hyd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCVH
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
      SUB_LOG(ISUB_LOG) = '/DFFLWB'
      K = KD(N)
      NPZ = NSZ(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        ZLABX = (1.D+0-ABS(ZLAB(MP,NB)))
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MP,NB).GT.EPSL ) THEN
          FLWP = XLW(MP,N)*RHOL(MP,N)
          FLWB = XLWB(MP,NB)*RHOLB(MP,NB)*ZLABX
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MP,NB).LT.-EPSL ) THEN
          FLWP = XLW(MP,N)*RHOL(MP,N)
          FLWB = XLWB(MP,NB)*RHOLB(MP,NB)*ZLABX
!
!---    No micro emulsion  ---
!
        ELSE
          FLWP = XLW(MP,N)*RHOL(MP,N)
          FLWB = XLWB(MP,NB)*RHOLB(MP,NB)
        ENDIF
        INDX = 3
        FLW = DIFMN( FLWB,FLWP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.26 ) THEN
          IF( WL(1,NPZ).LT.-EPSL ) THEN
            WLW(M,NPZ) = WL(M,NPZ)*FLW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.42 ) THEN
          IF( WL(1,NPZ).GT.EPSL ) THEN
            WLW(M,NPZ) = WL(M,NPZ)*FLW
          ENDIF
        ELSE
          WLW(M,NPZ) = WL(M,NPZ)*FLW - WTMW*WDLA(M,NPZ)
     &      - WTMW*WDS(M,NPZ)/WTMS - WTMW*WDLO(M,NPZ)
        ENDIF
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
!     Aqueous water flux on east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_hyd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCVH
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
      SUB_LOG(ISUB_LOG) = '/DFFLWE'
      I = ID(N)
      NQX = NSX(N)+1
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        ZLABX = (1.D+0-ABS(ZLAB(MN,NB)))
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MN,NB).GT.EPSL ) THEN
          FLWP = XLW(MN,N)*RHOL(MN,N)
          FLWB = XLWB(MN,NB)*RHOLB(MN,NB)*ZLABX
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MN,NB).LT.-EPSL ) THEN
          FLWP = XLW(MN,N)*RHOL(MN,N)
          FLWB = XLWB(MN,NB)*RHOLB(MN,NB)*ZLABX
!
!---    No micro emulsion  ---
!
        ELSE
          FLWP = XLW(MN,N)*RHOL(MN,N)
          FLWB = XLWB(MN,NB)*RHOLB(MN,NB)
        ENDIF
        INDX = 3
        FLW = DIFMN( FLWP,FLWB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.26 ) THEN
          IF( UL(1,NQX).GT.EPSL ) THEN
            ULW(M,NQX) = UL(M,NQX)*FLW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.42 ) THEN
          IF( UL(1,NQX).LT.-EPSL ) THEN
            ULW(M,NQX) = UL(M,NQX)*FLW
          ENDIF
        ELSE
          ULW(M,NQX) = UL(M,NQX)*FLW - WTMW*UDLA(M,NQX)
     &      - WTMW*UDS(M,NQX)/WTMS - WTMW*UDLO(M,NQX)
        ENDIF
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
!     Aqueous water flux on north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_hyd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCVH
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
      SUB_LOG(ISUB_LOG) = '/DFFLWN'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        ZLABX = (1.D+0-ABS(ZLAB(MN,NB)))
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MN,NB).GT.EPSL ) THEN
          FLWP = XLW(MN,N)*RHOL(MN,N)
          FLWB = XLWB(MN,NB)*RHOLB(MN,NB)*ZLABX
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MN,NB).LT.-EPSL ) THEN
          FLWP = XLW(MN,N)*RHOL(MN,N)
          FLWB = XLWB(MN,NB)*RHOLB(MN,NB)*ZLABX
!
!---    No micro emulsion  ---
!
        ELSE
          FLWP = XLW(MN,N)*RHOL(MN,N)
          FLWB = XLWB(MN,NB)*RHOLB(MN,NB)
        ENDIF
        INDX = 3
        FLW = DIFMN( FLWP,FLWB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.26 ) THEN
          IF( VL(1,NQY).GT.EPSL ) THEN
            VLW(M,NQY) = VL(M,NQY)*FLW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.42 ) THEN
          IF( VL(1,NQY).LT.-EPSL ) THEN
            VLW(M,NQY) = VL(M,NQY)*FLW
          ENDIF
        ELSE
          VLW(M,NQY) = VL(M,NQY)*FLW - WTMW*VDLA(M,NQY)
     &      - WTMW*VDS(M,NQY)/WTMS - WTMW*VDLO(M,NQY)
        ENDIF
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
!     Aqueous water flux on south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_hyd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCVH
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
      SUB_LOG(ISUB_LOG) = '/DFFLWS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        ZLABX = (1.D+0-ABS(ZLAB(MP,NB)))
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MP,NB).GT.EPSL ) THEN
          FLWP = XLW(MP,N)*RHOL(MP,N)
          FLWB = XLWB(MP,NB)*RHOLB(MP,NB)*ZLABX
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MP,NB).LT.-EPSL ) THEN
          FLWP = XLW(MP,N)*RHOL(MP,N)
          FLWB = XLWB(MP,NB)*RHOLB(MP,NB)*ZLABX
!
!---    No micro emulsion  ---
!
        ELSE
          FLWP = XLW(MP,N)*RHOL(MP,N)
          FLWB = XLWB(MP,NB)*RHOLB(MP,NB)
        ENDIF
        INDX = 3
        FLW = DIFMN( FLWB,FLWP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.26 ) THEN
          IF( VL(1,NPY).LT.-EPSL ) THEN
            VLW(M,NPY) = VL(M,NPY)*FLW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.42 ) THEN
          IF( VL(1,NPY).GT.EPSL ) THEN
            VLW(M,NPY) = VL(M,NPY)*FLW
          ENDIF
        ELSE
          VLW(M,NPY) = VL(M,NPY)*FLW - WTMW*VDLA(M,NPY)
     &      - WTMW*VDS(M,NPY)/WTMS - WTMW*VDLO(M,NPY)
        ENDIF
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
!     Aqueous water flux on top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_hyd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCVH
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
      SUB_LOG(ISUB_LOG) = '/DFFLWT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        ZLABX = (1.D+0-ABS(ZLAB(MN,NB)))
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MN,NB).GT.EPSL ) THEN
          FLWP = XLW(MN,N)*RHOL(MN,N)
          FLWB = XLWB(MN,NB)*RHOLB(MN,NB)*ZLABX
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MN,NB).LT.-EPSL ) THEN
          FLWP = XLW(MN,N)*RHOL(MN,N)
          FLWB = XLWB(MN,NB)*RHOLB(MN,NB)*ZLABX
!
!---    No micro emulsion  ---
!
        ELSE
          FLWP = XLW(MN,N)*RHOL(MN,N)
          FLWB = XLWB(MN,NB)*RHOLB(MN,NB)
        ENDIF
        INDX = 3
        FLW = DIFMN( FLWP,FLWB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.26 ) THEN
          IF( WL(1,NQZ).GT.EPSL ) THEN
            WLW(M,NQZ) = WL(M,NQZ)*FLW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.42 ) THEN
          IF( WL(1,NQZ).LT.-EPSL ) THEN
            WLW(M,NQZ) = WL(M,NQZ)*FLW
          ENDIF
        ELSE
          WLW(M,NQZ) = WL(M,NQZ)*FLW - WTMW*WDLA(M,NQZ)
     &      - WTMW*WDS(M,NQZ)/WTMS - WTMW*WDLO(M,NQZ)
        ENDIF
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
!     Aqueous water flux on west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!     Last modified by MD White, PNNL, 28 March 2005.
!     dff_hyd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCVH
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
      SUB_LOG(ISUB_LOG) = '/DFFLWW'
      I = ID(N)
      NPX = NSX(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        ZLABX = (1.D+0-ABS(ZLAB(MP,NB)))
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MP,NB).GT.EPSL ) THEN
          FLWP = XLW(MP,N)*RHOL(MP,N)
          FLWB = XLWB(MP,NB)*RHOLB(MP,NB)*ZLABX
!
!---    Vapor micro emulsion  ---
!
        ELSEIF( ZLAB(MP,NB).LT.-EPSL ) THEN
          FLWP = XLW(MP,N)*RHOL(MP,N)
          FLWB = XLWB(MP,NB)*RHOLB(MP,NB)*ZLABX
!
!---    No micro emulsion  ---
!
        ELSE
          FLWP = XLW(MP,N)*RHOL(MP,N)
          FLWB = XLWB(MP,NB)*RHOLB(MP,NB)
        ENDIF
        INDX = 3
        FLW = DIFMN( FLWB,FLWP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.26 ) THEN
          IF( UL(1,NPX).LT.-EPSL ) THEN
            ULW(M,NPX) = UL(M,NPX)*FLW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.42 ) THEN
          IF( UL(1,NPX).GT.EPSL ) THEN
            ULW(M,NPX) = UL(M,NPX)*FLW
          ENDIF
        ELSE
          ULW(M,NPX) = UL(M,NPX)*FLW - WTMW*UDLA(M,NPX)
     &      - WTMW*UDS(M,NPX)/WTMS - WTMW*UDLO(M,NPX)
        ENDIF
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
      SUBROUTINE DFFNAB( N,NB )
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
!     Compute CO2 liquid-CO2 fluxes on bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 30 October 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE CONST
      USE BCVP
      USE BCVH
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFNAB'
      K = KD(N)
      NPZ = NSZ(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MP,NB).GT.EPSL ) THEN
          FLAP = 0.D+0
          FLAB = XNCB(1,MP,NB)*RHONB(MP,NB)*ZLAB(MP,NB)
          INDX = 3
          FLA = DIFMN( FLAB,FLAP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLA = 0.D+0
        ENDIF
        FNAB = XNCB(1,MP,NB)*RHONB(MP,NB)
        FNAP = XNC(1,MP,N)*RHON(MP,N)
        INDX = 3
        FNA = DIFMN( FNAB,FNAP,DZGF(N),DZGF(N),WN(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(4,NB).EQ.26 ) THEN
          IF( WN(1,NPZ).LT.-EPSL ) THEN
            WNC(1,M,NPZ) = WN(M,NPZ)*FNA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(4,NB).EQ.42 ) THEN
          IF( WN(1,NPZ).GT.EPSL ) THEN
            WNC(1,M,NPZ) = WN(M,NPZ)*FNA
          ENDIF
        ELSE
          WDNC(1,M,NPZ) = -(WDNC(2,M,NPZ)+WDNC(3,M,NPZ))
          WNC(1,M,NPZ) = WL(M,NPZ)*FLA + WN(M,NPZ)*FNA 
     &      + WTMA*WDNC(1,M,NPZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNAB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNAE( N,NB )
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
!     Compute CO2 liquid-CO2 fluxes on east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 30 October 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE CONST
      USE BCVP
      USE BCVH
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFNAE'
      I = ID(N)
      NQX = NSX(N)+1
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MN,NB).GT.EPSL ) THEN
          FLAP = 0.D+0
          FLAB = XNCB(1,MN,NB)*RHONB(MN,NB)*ZLAB(MN,NB)
          INDX = 3
          FLA = DIFMN( FLAP,FLAB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLA = 0.D+0
        ENDIF
        FNAB = XNCB(1,MN,NB)*RHONB(MN,NB)
        FNAP = XNC(1,MN,N)*RHON(MN,N)
        INDX = 3
        FNA = DIFMN( FNAP,FNAB,DXGF(N),DXGF(N),UN(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(4,NB).EQ.26 ) THEN
          IF( UN(1,NQX).GT.EPSL ) THEN
            UNC(1,M,NQX) = UN(M,NQX)*FNA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(4,NB).EQ.42 ) THEN
          IF( UN(1,NQX).LT.-EPSL ) THEN
            UNC(1,M,NQX) = UN(M,NQX)*FNA
          ENDIF
        ELSE
          UDNC(1,M,NQX) = -(UDNC(2,M,NQX)+UDNC(3,M,NQX))
          UNC(1,M,NQX) = UL(M,NQX)*FLA + UN(M,NQX)*FNA 
     &      + WTMA*UDNC(1,M,NQX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNAE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNAN( N,NB )
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
!     Compute CO2 liquid-CO2 fluxes on north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 30 October 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE CONST
      USE BCVP
      USE BCVH
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFNAN'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MN,NB).GT.EPSL ) THEN
          FLAP = 0.D+0
          FLAB = XNCB(1,MN,NB)*RHONB(MN,NB)*ZLAB(MN,NB)
          INDX = 3
          FLA = DIFMN( FLAP,FLAB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLA = 0.D+0
        ENDIF
        FNAB = XNCB(1,MN,NB)*RHONB(MN,NB)
        FNAP = XNC(1,MN,N)*RHON(MN,N)
        INDX = 3
        FNA = DIFMN( FNAP,FNAB,DYGF(N),DYGF(N),VN(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(4,NB).EQ.26 ) THEN
          IF( VN(1,NQY).GT.EPSL ) THEN
            VNC(1,M,NQY) = VN(M,NQY)*FNA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(4,NB).EQ.42 ) THEN
          IF( VN(1,NQY).LT.-EPSL ) THEN
            VNC(1,M,NQY) = VN(M,NQY)*FNA
          ENDIF
        ELSE
          VDNC(1,M,NQY) = -(VDNC(2,M,NQY)+VDNC(3,M,NQY))
          VNC(1,M,NQY) = VL(M,NQY)*FLA + VN(M,NQY)*FNA 
     &      + WTMA*VDNC(1,M,NQY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNAN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNAS( N,NB )
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
!     Compute CO2 liquid-CO2 fluxes on south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 30 October 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE CONST
      USE BCVP
      USE BCVH
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFNAS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MP,NB).GT.EPSL ) THEN
          FLAP = 0.D+0
          FLAB = XNCB(1,MP,NB)*RHONB(MP,NB)*ZLAB(MP,NB)
          INDX = 3
          FLA = DIFMN( FLAB,FLAP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLA = 0.D+0
        ENDIF
        FNAB = XNCB(1,MP,NB)*RHONB(MP,NB)
        FNAP = XNC(1,MP,N)*RHON(MP,N)
        INDX = 3
        FNA = DIFMN( FNAB,FNAP,DYGF(N),DYGF(N),VN(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(4,NB).EQ.26 ) THEN
          IF( VN(1,NPY).LT.-EPSL ) THEN
            VNC(1,M,NPY) = VN(M,NPY)*FNA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(4,NB).EQ.42 ) THEN
          IF( VN(1,NPY).GT.EPSL ) THEN
            VNC(1,M,NPY) = VN(M,NPY)*FNA
          ENDIF
        ELSE
          VDNC(1,M,NPY) = -(VDNC(2,M,NPY)+VDNC(3,M,NPY))
          VNC(1,M,NPY) = VL(M,NPY)*FLA + VN(M,NPY)*FNA 
     &      + WTMA*VDNC(1,M,NPY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNAS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNAT( N,NB )
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
!     Compute CO2 liquid-CO2 fluxes on top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 30 October 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE CONST
      USE BCVP
      USE BCVH
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFNAT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MN,NB).GT.EPSL ) THEN
          FLAP = 0.D+0
          FLAB = XNCB(1,MN,NB)*RHONB(MN,NB)*ZLAB(MN,NB)
          INDX = 3
          FLA = DIFMN( FLAP,FLAB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLA = 0.D+0
        ENDIF
        FNAB = XNCB(1,MN,NB)*RHONB(MN,NB)
        FNAP = XNC(1,MN,N)*RHON(MN,N)
        INDX = 3
        FNA = DIFMN( FNAP,FNAB,DZGF(N),DZGF(N),WN(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(4,NB).EQ.26 ) THEN
          IF( WN(1,NQZ).GT.EPSL ) THEN
            WNC(1,M,NQZ) = WN(M,NQZ)*FNA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(4,NB).EQ.42 ) THEN
          IF( WN(1,NQZ).LT.-EPSL ) THEN
            WNC(1,M,NQZ) = WN(M,NQZ)*FNW
          ENDIF
        ELSE
          WDNC(1,M,NQZ) = -(WDNC(2,M,NQZ)+WDNC(3,M,NQZ))
          WNC(1,M,NQZ) = WL(M,NQZ)*FLA + WN(M,NQZ)*FNA 
     &      + WTMA*WDNC(1,M,NQZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNAT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNAW( N,NB )
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
!     Compute CO2 liquid-CO2 fluxes on west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 30 October 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE CONST
      USE BCVP
      USE BCVH
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFNAW'
      I = ID(N)
      NPX = NSX(N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MP,NB).GT.EPSL ) THEN
          FLAP = 0.D+0
          FLAB = XLWB(MP,NB)*RHOLB(MP,NB)*ZLAB(MP,NB)
          INDX = 3
          FLA = DIFMN( FLAB,FLAP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLA = 0.D+0
        ENDIF
        FNAB = XNCB(1,MP,NB)*RHONB(MP,NB)
        FNAP = XNC(1,MP,N)*RHON(MP,N)
        INDX = 3
        FNA = DIFMN( FNAB,FNAP,DXGF(N),DXGF(N),UN(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(4,NB).EQ.26 ) THEN
          IF( UN(1,NPX).LT.-EPSL ) THEN
            UNC(1,M,NPX) = UN(M,NPX)*FNA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(4,NB).EQ.42 ) THEN
          IF( UN(1,NPX).GT.EPSL ) THEN
            UNC(1,M,NPX) = UN(M,NPX)*FNA
          ENDIF
        ELSE
          UDNC(1,M,NPX) = -(UDNC(2,M,NPX)+UDNC(3,M,NPX))
          UNC(1,M,NPX) = UL(M,NPX)*FLA + UN(M,NPX)*FNA 
     &      + WTMA*UDNC(1,M,NPX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNAW group  ---
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
!     Compute CH4 liquid-CO2 fluxes on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVN
      USE FDVGC
      USE CONST
      USE BCVP
      USE BCVN
      USE BCVH
      USE BCVGC
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
      DXMNO = XMNCB(2,2,NB) - XMNC(2,2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MP,NB).GT.EPSL ) THEN
          FLOP = 0.D+0
          FLOB = XNCB(2,MP,NB)*RHONB(MP,NB)*ZLAB(MP,NB)
          INDX = 3
          FLO = DIFMN( FLOB,FLOP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLO = 0.D+0
        ENDIF
        FNOP = XNC(2,MP,N)*RHON(MP,N)
        FNOB = XNCB(2,MP,NB)*RHONB(MP,NB)
        INDX = 3
        FNO = DIFMN( FNOB,FNOP,DZGF(N),DZGF(N),WN(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(4,NB).EQ.26 ) THEN
          IF( WN(1,NPZ).LT.-EPSL ) THEN
            WNC(2,M,NPZ) = WN(M,NPZ)*FNO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(4,NB).EQ.42 ) THEN
          IF( WN(1,NPZ).GT.EPSL ) THEN
            WNC(2,M,NPZ) = WN(M,NPZ)*FNO
          ENDIF
        ELSE
          DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*DFNC(2,MP,N)
     &      *RHOMN(MP,N)
          DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNCB(2,MP,NB)
     &      *RHOMNB(MP,NB)
          INDX = 13
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMNO,INDX )
          WDNC(2,M,NPZ) = DFM*(XMNCB(2,MP,NB)
     &      - XMNC(2,MP,N))/(5.D-1*DZGF(N))
          WNC(2,M,NPZ) = WL(M,NPZ)*FLO + WN(M,NPZ)*FNO 
     &      + WTMO*WDNC(2,M,NPZ)
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
!     Compute CH4 liquid-CO2 fluxes on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVN
      USE FDVGC
      USE CONST
      USE BCVP
      USE BCVN
      USE BCVH
      USE BCVGC
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
      DXMNO = XMNC(2,2,N) - XMNCB(2,2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MN,NB).GT.EPSL ) THEN
          FLOP = 0.D+0
          FLOB = XNCB(2,MN,NB)*RHONB(MN,NB)*ZLAB(MN,NB)
          INDX = 3
          FLO = DIFMN( FLOP,FLOB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLO = 0.D+0
        ENDIF
        FNOP = XNC(2,MN,N)*RHON(MN,N)
        FNOB = XNCB(2,MN,NB)*RHONB(MN,NB)
        INDX = 3
        FNO = DIFMN( FNOP,FNOB,DXGF(N),DXGF(N),UN(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(4,NB).EQ.26 ) THEN
          IF( UN(1,NQX).GT.EPSL ) THEN
            UNC(2,M,NQX) = UN(M,NQX)*FNO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(4,NB).EQ.42 ) THEN
          IF( UN(1,NQX).LT.-EPSL ) THEN
            UNC(2,M,NQX) = UN(M,NQX)*FNO
          ENDIF
        ELSE
          DFP = TORN(MN,N)*PORD(MN,N)*SN(MN,N)*DFNC(2,MN,N)
     &      *RHOMN(MN,N)
          DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNCB(2,MN,NB)
     &      *RHOMNB(MN,NB)
          INDX = 13
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMNO,INDX )
          UDNC(2,M,NQX) = DFM*(XMNC(2,MN,N)
     &      - XMNCB(2,MN,NB))/(5.D-1*DXGF(N))
          UNC(2,M,NQX) = UL(M,NQX)*FLO + UN(M,NQX)*FNO 
     &      + WTMO*UDNC(2,M,NQX)
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
!     Compute CH4 liquid-CO2 fluxes on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVN
      USE FDVGC
      USE CONST
      USE BCVP
      USE BCVN
      USE BCVH
      USE BCVGC
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
      DXMNO = XMNC(2,2,N) - XMNCB(2,2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MN,NB).GT.EPSL ) THEN
          FLOP = 0.D+0
          FLOB = XNCB(2,MN,NB)*RHONB(MN,NB)*ZLAB(MN,NB)
          INDX = 3
          FLO = DIFMN( FLOP,FLOB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLO = 0.D+0
        ENDIF
        FNOP = XNC(2,MN,N)*RHON(MN,N)
        FNOB = XNCB(2,MN,NB)*RHONB(MN,NB)
        INDX = 3
        FNO = DIFMN( FNOP,FNOB,DYGF(N),DYGF(N),VN(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(4,NB).EQ.26 ) THEN
          IF( VN(1,NQY).GT.EPSL ) THEN
            VNC(2,M,NQY) = VN(M,NQY)*FNO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(4,NB).EQ.42 ) THEN
          IF( VN(1,NQY).LT.-EPSL ) THEN
            VNC(2,M,NQY) = VN(M,NQY)*FNO
          ENDIF
        ELSE
          DFP = TORN(MN,N)*PORD(MN,N)*SN(MN,N)*DFNC(2,MN,N)
     &      *RHOMN(MN,N)
          DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNCB(2,MN,NB)
     &      *RHOMNB(MN,NB)
          INDX = 13
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMNO,INDX )
          VDNC(2,M,NQY) = DFM*(XMNC(2,MN,N)
     &      - XMNCB(2,MN,NB))/((5.D-1*DYGF(N))*RP(I))
          VNC(2,M,NQY) = VL(M,NQY)*FLO + VN(M,NQY)*FNO 
     &      + WTMO*VDNC(2,M,NQY)
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
!     Compute CH4 liquid-CO2 fluxes on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVN
      USE FDVGC
      USE CONST
      USE BCVP
      USE BCVN
      USE BCVH
      USE BCVGC
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
      DXMNO = XMNCB(2,2,NB) - XMNC(2,2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MP,NB).GT.EPSL ) THEN
          FLOP = 0.D+0
          FLOB = XNCB(2,MP,NB)*RHONB(MP,NB)*ZLAB(MP,NB)
          INDX = 3
          FLO = DIFMN( FLOB,FLOP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLO = 0.D+0
        ENDIF
        FNOP = XNC(2,MP,N)*RHON(MP,N)
        FNOB = XNCB(2,MP,NB)*RHONB(MP,NB)
        INDX = 3
        FNO = DIFMN( FNOB,FNOP,DYGF(N),DYGF(N),VN(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(4,NB).EQ.26 ) THEN
          IF( VN(1,NPY).LT.-EPSL ) THEN
            VNC(2,M,NPY) = VN(M,NPY)*FNO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(4,NB).EQ.42 ) THEN
          IF( VN(1,NPY).GT.EPSL ) THEN
            VNC(2,M,NPY) = VN(M,NPY)*FNO
          ENDIF
        ELSE
          DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*DFNC(2,MP,N)
     &      *RHOMN(MP,N)
          DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNCB(2,MP,NB)
     &      *RHOMNB(MP,NB)
          INDX = 13
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMNO,INDX )
          VDNC(2,M,NPY) = DFM*(XMNCB(2,MP,NB)
     &      - XMNC(2,MP,N))/((5.D-1*DYGF(N))*RP(I))
          VNC(2,M,NPY) = VL(M,NPY)*FLO + VN(M,NPY)*FNO 
     &      + WTMO*VDNC(2,M,NPY)
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
!     Compute CH4 liquid-CO2 fluxes on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVN
      USE FDVGC
      USE CONST
      USE BCVP
      USE BCVN
      USE BCVH
      USE BCVGC
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
      DXMNO = XMNC(2,2,N)-XMNCB(2,2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MN,NB).GT.EPSL ) THEN
          FLOP = 0.D+0
          FLOB = XNCB(2,MN,NB)*RHONB(MN,NB)*ZLAB(MN,NB)
          INDX = 3
          FLO = DIFMN( FLOP,FLOB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLO = 0.D+0
        ENDIF
        FNOP = XNC(2,MN,N)*RHON(MN,N)
        FNOB = XNCB(2,MN,NB)*RHONB(MN,NB)
        INDX = 3
        FNO = DIFMN( FNOP,FNOB,DZGF(N),DZGF(N),WN(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(4,NB).EQ.26 ) THEN
          IF( WN(1,NQZ).GT.EPSL ) THEN
            WNC(2,M,NQZ) = WN(M,NQZ)*FNO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(4,NB).EQ.42 ) THEN
          IF( WN(1,NQZ).LT.-EPSL ) THEN
            WNC(2,M,NQZ) = WN(M,NQZ)*FNO
          ENDIF
        ELSE
          DFP = TORN(MN,N)*PORD(MN,N)*SN(MN,N)*DFNC(2,MN,N)
     &      *RHOMN(MN,N)
          DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNCB(2,MN,NB)
     &      *RHOMNB(MN,NB)
          INDX = 13
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMNO,INDX )
          WDNC(2,M,NQZ) = DFM*(XMNC(2,MN,N)
     &      - XMNCB(2,MN,NB))/(5.D-1*DZGF(N))
          WNC(2,M,NQZ) = WL(M,NQZ)*FLO + WN(M,NQZ)*FNO 
     &      + WTMO*WDNC(2,M,NQZ)
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
!     Compute CH4 liquid-CO2 fluxes on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
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
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVN
      USE FDVGC
      USE CONST
      USE BCVP
      USE BCVN
      USE BCVH
      USE BCVGC
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
      DXMNO = XMNCB(2,2,NB)-XMNC(2,2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MP,NB).GT.EPSL ) THEN
          FLOP = 0.D+0
          FLOB = XNCB(2,MP,NB)*RHONB(MP,NB)*ZLAB(MP,NB)
          INDX = 3
          FLO = DIFMN( FLOB,FLOP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLO = 0.D+0
        ENDIF
        FNOP = XNC(2,MP,N)*RHON(MP,N)
        FNOB = XNCB(2,MP,NB)*RHONB(MP,NB)
        INDX = 3
        FNO = DIFMN( FNOB,FNOP,DXGF(N),DXGF(N),UN(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(4,NB).EQ.26 ) THEN
          IF( UN(1,NPX).LT.-EPSL ) THEN
            UNC(2,M,NPX) = UN(M,NPX)*FNO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(4,NB).EQ.42 ) THEN
          IF( UN(1,NPX).GT.EPSL ) THEN
            UNC(2,M,NPX) = UN(M,NPX)*FNO
          ENDIF
        ELSE
          DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*DFNC(2,MP,N)
     &      *RHOMN(MP,N)
          DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNCB(2,MP,NB)
     &      *RHOMNB(MP,NB)
          INDX = 13
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMNO,INDX )
          UDNC(2,M,NPX) = DFM*(XMNCB(2,MP,NB)
     &      - XMNC(2,MP,N))/(5.D-1*DXGF(N))
          UNC(2,M,NPX) = UL(M,NPX)*FLO + UN(M,NPX)*FNO 
     &      + WTMO*UDNC(2,M,NPX)
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

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNWB( N,NB )
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
!     Compute CH4 liquid-CO2 fluxes on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVN
      USE FDVGC
      USE CONST
      USE BCVP
      USE BCVN
      USE BCVH
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFNWB'
      K = KD(N)
      NPZ = NSZ(N)
      DXMNW = XMNCB(3,2,NB) - XMNC(3,2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MP,NB).GT.EPSL ) THEN
          FLWP = 0.D+0
          FLWB = XNCB(3,MP,NB)*RHONB(MP,NB)*ZLAB(MP,NB)
          INDX = 3
          FLW = DIFMN( FLWB,FLWP,DZGF(N),DZGF(N),WL(1,NPZ),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLW = 0.D+0
        ENDIF
        FNWP = XNC(3,MP,N)*RHON(MP,N)
        FNWB = XNCB(3,MP,NB)*RHONB(MP,NB)
        INDX = 3
        FNW = DIFMN( FNWB,FNWP,DZGF(N),DZGF(N),WN(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(4,NB).EQ.26 ) THEN
          IF( WN(1,NPZ).LT.-EPSL ) THEN
            WNC(3,M,NPZ) = WN(M,NPZ)*FNW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(4,NB).EQ.42 ) THEN
          IF( WN(1,NPZ).GT.EPSL ) THEN
            WNC(3,M,NPZ) = WN(M,NPZ)*FNW
          ENDIF
        ELSE
          DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*DFNC(3,MP,N)
     &      *RHOMN(MP,N)
          DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNCB(3,MP,NB)
     &      *RHOMNB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMNW,INDX )
          WDNC(3,M,NPZ) = DFM*(XMNCB(3,MP,NB)
     &      - XMNC(3,MP,N))/(5.D-1*DZGF(N))
          WNC(3,M,NPZ) = WL(M,NPZ)*FLW + WN(M,NPZ)*FNW 
     &      + WTMW*WDNC(3,M,NPZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNWB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNWE( N,NB )
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
!     Compute CH4 liquid-CO2 fluxes on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVN
      USE FDVGC
      USE CONST
      USE BCVP
      USE BCVN
      USE BCVH
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFNWE'
      I = ID(N)
      NQX = NSX(N)+1
      DXMNW = XMNC(3,2,N) - XMNCB(3,2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MN,NB).GT.EPSL ) THEN
          FLWP = 0.D+0
          FLWB = XNCB(3,MN,NB)*RHONB(MN,NB)*ZLAB(MN,NB)
          INDX = 3
          FLW = DIFMN( FLWP,FLWB,DXGF(N),DXGF(N),UL(1,NQX),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLW = 0.D+0
        ENDIF
        FNWP = XNC(3,MN,N)*RHON(MN,N)
        FNWB = XNCB(3,MN,NB)*RHONB(MN,NB)
        INDX = 3
        FNW = DIFMN( FNWP,FNWB,DXGF(N),DXGF(N),UN(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(4,NB).EQ.26 ) THEN
          IF( UN(1,NQX).GT.EPSL ) THEN
            UNC(3,M,NQX) = UN(M,NQX)*FNW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(4,NB).EQ.42 ) THEN
          IF( UN(1,NQX).LT.-EPSL ) THEN
            UNC(3,M,NQX) = UN(M,NQX)*FNW
          ENDIF
        ELSE
          DFP = TORN(MN,N)*PORD(MN,N)*SN(MN,N)*DFNC(3,MN,N)
     &      *RHOMN(MN,N)
          DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNCB(3,MN,NB)
     &      *RHOMNB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMNW,INDX )
          UDNC(3,M,NQX) = DFM*(XMNC(3,MN,N)
     &      - XMNCB(3,MN,NB))/(5.D-1*DXGF(N))
          UNC(3,M,NQX) = UL(M,NQX)*FLW + UN(M,NQX)*FNW 
     &      + WTMW*UDNC(3,M,NQX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNWE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNWN( N,NB )
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
!     Compute CH4 liquid-CO2 fluxes on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVN
      USE FDVGC
      USE CONST
      USE BCVP
      USE BCVN
      USE BCVH
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFNWN'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMNW = XMNC(3,2,N) - XMNCB(3,2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MN,NB).GT.EPSL ) THEN
          FLWP = 0.D+0
          FLWB = XNCB(3,MN,NB)*RHONB(MN,NB)*ZLAB(MN,NB)
          INDX = 3
          FLW = DIFMN( FLWP,FLWB,DYGF(N),DYGF(N),VL(1,NQY),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLW = 0.D+0
        ENDIF
        FNWP = XNC(3,MN,N)*RHON(MN,N)
        FNWB = XNCB(3,MN,NB)*RHONB(MN,NB)
        INDX = 3
        FNW = DIFMN( FNWP,FNWB,DYGF(N),DYGF(N),VN(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(4,NB).EQ.26 ) THEN
          IF( VN(1,NQY).GT.EPSL ) THEN
            VNC(3,M,NQY) = VN(M,NQY)*FNW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(4,NB).EQ.42 ) THEN
          IF( VN(1,NQY).LT.-EPSL ) THEN
            VNC(3,M,NQY) = VN(M,NQY)*FNW
          ENDIF
        ELSE
          DFP = TORN(MN,N)*PORD(MN,N)*SN(MN,N)*DFNC(3,MN,N)
     &      *RHOMN(MN,N)
          DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNCB(3,MN,NB)
     &      *RHOMNB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMNW,INDX )
          VDNC(3,M,NQY) = DFM*(XMNC(3,MN,N)
     &      - XMNCB(3,MN,NB))/((5.D-1*DYGF(N))*RP(I))
          VNC(3,M,NQY) = VL(M,NQY)*FLW + VN(M,NQY)*FNW 
     &      + WTMW*VDNC(3,M,NQY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNWN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNWS( N,NB )
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
!     Compute CH4 liquid-CO2 fluxes on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVN
      USE FDVGC
      USE CONST
      USE BCVP
      USE BCVN
      USE BCVH
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFNWS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMNW = XMNCB(3,2,NB) - XMNC(3,2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MP,NB).GT.EPSL ) THEN
          FLWP = 0.D+0
          FLWB = XNCB(3,MP,NB)*RHONB(MP,NB)*ZLAB(MP,NB)
          INDX = 3
          FLW = DIFMN( FLWB,FLWP,DYGF(N),DYGF(N),VL(1,NPY),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLW = 0.D+0
        ENDIF
        FNWP = XNC(3,MP,N)*RHON(MP,N)
        FNWB = XNCB(3,MP,NB)*RHONB(MP,NB)
        INDX = 3
        FNW = DIFMN( FNWB,FNWP,DYGF(N),DYGF(N),VN(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(4,NB).EQ.26 ) THEN
          IF( VN(1,NPY).LT.-EPSL ) THEN
            VNC(3,M,NPY) = VN(M,NPY)*FNW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(4,NB).EQ.42 ) THEN
          IF( VN(1,NPY).GT.EPSL ) THEN
            VNC(3,M,NPY) = VN(M,NPY)*FNW
          ENDIF
        ELSE
          DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*DFNC(3,MP,N)
     &      *RHOMN(MP,N)
          DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNCB(3,MP,NB)
     &      *RHOMNB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMNW,INDX )
          VDNC(3,M,NPY) = DFM*(XMNCB(3,MP,NB)
     &      - XMNC(3,MP,N))/((5.D-1*DYGF(N))*RP(I))
          VNC(3,M,NPY) = VL(M,NPY)*FLW + VN(M,NPY)*FNW 
     &      + WTMW*VDNC(3,M,NPY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNWS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNWT( N,NB )
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
!     Compute CH4 liquid-CO2 fluxes on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVN
      USE FDVGC
      USE CONST
      USE BCVP
      USE BCVN
      USE BCVH
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFNWT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMNW = XMNC(3,2,N)-XMNCB(3,2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MN,NB).GT.EPSL ) THEN
          FLWP = 0.D+0
          FLWB = XNCB(3,MN,NB)*RHONB(MN,NB)*ZLAB(MN,NB)
          INDX = 3
          FLW = DIFMN( FLWP,FLWB,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLW = 0.D+0
        ENDIF
        FNWP = XNC(3,MN,N)*RHON(MN,N)
        FNWB = XNCB(3,MN,NB)*RHONB(MN,NB)
        INDX = 3
        FNW = DIFMN( FNWP,FNWB,DZGF(N),DZGF(N),WN(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(4,NB).EQ.26 ) THEN
          IF( WN(1,NQZ).GT.EPSL ) THEN
            WNC(3,M,NQZ) = WN(M,NQZ)*FNW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(4,NB).EQ.42 ) THEN
          IF( WN(1,NQZ).LT.-EPSL ) THEN
            WNC(3,M,NQZ) = WN(M,NQZ)*FNW
          ENDIF
        ELSE
          DFP = TORN(MN,N)*PORD(MN,N)*SN(MN,N)*DFNC(3,MN,N)
     &      *RHOMN(MN,N)
          DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNCB(3,MN,NB)
     &      *RHOMNB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMNW,INDX )
          WDNC(3,M,NQZ) = DFM*(XMNC(3,MN,N)
     &      - XMNCB(3,MN,NB))/(5.D-1*DZGF(N))
          WNC(3,M,NQZ) = WL(M,NQZ)*FLW + WN(M,NQZ)*FNW 
     &      + WTMW*WDNC(3,M,NQZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNWT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNWW( N,NB )
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
!     Compute CH4 liquid-CO2 fluxes on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 March 2005.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVN
      USE FDVGC
      USE CONST
      USE BCVP
      USE BCVN
      USE BCVH
      USE BCVGC
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
      SUB_LOG(ISUB_LOG) = '/DFFNWW'
      I = ID(N)
      NPX = NSX(N)
      DXMNW = XMNCB(3,2,NB)-XMNC(3,2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
!
!---    Liquid micro emulsion  ---
!
        IF( ZLAB(MP,NB).GT.EPSL ) THEN
          FLWP = 0.D+0
          FLWB = XNCB(3,MP,NB)*RHONB(MP,NB)*ZLAB(MP,NB)
          INDX = 3
          FLW = DIFMN( FLWB,FLWP,DXGF(N),DXGF(N),UL(1,NPX),INDX )
!
!---    No micro emulsion  ---
!
        ELSE
          FLW = 0.D+0
        ENDIF
        FNWP = XNC(3,MP,N)*RHON(MP,N)
        FNWB = XNCB(3,MP,NB)*RHONB(MP,NB)
        INDX = 3
        FNW = DIFMN( FNWB,FNWP,DXGF(N),DXGF(N),UN(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(4,NB).EQ.26 ) THEN
          IF( UN(1,NPX).LT.-EPSL ) THEN
            UNC(3,M,NPX) = UN(M,NPX)*FNW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(4,NB).EQ.42 ) THEN
          IF( UN(1,NPX).GT.EPSL ) THEN
            UNC(3,M,NPX) = UN(M,NPX)*FNW
          ENDIF
        ELSE
          DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*DFNC(3,MP,N)
     &      *RHOMN(MP,N)
          DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNCB(3,MP,NB)
     &      *RHOMNB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMNW,INDX )
          UDNC(3,M,NPX) = DFM*(XMNCB(3,MP,NB)
     &      - XMNC(3,MP,N))/(5.D-1*DXGF(N))
          UNC(3,M,NPX) = UL(M,NPX)*FLW + UN(M,NPX)*FNW 
     &      + WTMW*UDNC(3,M,NPX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNWW group  ---
!
      RETURN
      END


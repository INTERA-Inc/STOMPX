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
!     Diffusive CO2 gas fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
        DO N = 1,NFLD
          I = ID(N)
          IF( I.EQ.1 ) CYCLE
          J = JD(N)
          K = KD(N)
          NW = N-1
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) CYCLE
          NPX = NSX(N)
          DXMGA = XMGA(2,NW) - XMGA(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGA(MP,N)*RHOMG(MP,N)
            DFW = TORG(MN,NW)*PORD(MN,NW)*SG(MN,NW)*
     &        DFGA(MN,NW)*RHOMG(MN,NW)
            INDX = 12
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXMGA,INDX )
            UDGA(M,NPX) = DFM*(XMGA(MN,NW)-XMGA(MP,N))/DXGP(NPX)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            FGAW = XGA(MN,NW)*RHOG(MN,NW)
            INDX = 3
            FGA = DIFMN( FGAW,FGAP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            UGA(M,NPX) = UG(M,NPX)*FGA + WTMA*UDGA(M,NPX)
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction vapor mole diffusion, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) CYCLE
          K = KD(N)
          NS = N-IFLD
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) CYCLE
          NPY = NSY(N)
          DXMGA = XMGA(2,NS) - XMGA(2,N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGA(MP,N)*RHOMG(MP,N)
            DFS = TORG(MN,NS)*PORD(MN,NS)*SG(MN,NS)*
     &        DFGA(MN,NS)*RHOMG(MN,NS)
            INDX = 12
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXMGA,INDX )
            VDGA(M,NPY) = DFM*(XMGA(MN,NS)-XMGA(MP,N))/(DYGP(NPY)*RP(I))
            FGAP = XGA(MP,N)*RHOG(MP,N)
            FGAS = XGA(MN,NS)*RHOG(MN,NS)
            INDX = 3
            FGA = DIFMN( FGAS,FGAP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
            VGA(M,NPY) = VG(M,NPY)*FGA + WTMA*VDGA(M,NPY)
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction vapor mole diffusion, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) CYCLE
          NB = N-IJFLD
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) CYCLE
          NPZ = NSZ(N)
          DXMGA = XMGA(2,NB) - XMGA(2,N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGA(MP,N)*RHOMG(MP,N)
            DFB = TORG(MN,NB)*PORD(MN,NB)*SG(MN,NB)*
     &        DFGA(MN,NB)*RHOMG(MN,NB)
            INDX = 12
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXMGA,INDX )
            WDGA(M,NPZ) = DFM*(XMGA(MN,NB)-XMGA(MP,N))/DZGP(NPZ)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            FGAB = XGA(MN,NB)*RHOG(MN,NB)
            INDX = 3
            FGA = DIFMN( FGAB,FGAP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
            WGA(M,NPZ) = WG(M,NPZ)*FGA + WTMA*WDGA(M,NPZ)
          ENDDO
        ENDDO
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
      SUBROUTINE DFFGN
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
!     Diffusive N2 gas fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      USE FDVH
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
      SUB_LOG(ISUB_LOG) = '/DFFGN'
!
!---  X-direction vapor mole diffusion, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          IF( I.EQ.1 ) CYCLE
          J = JD(N)
          K = KD(N)
          NW = N-1
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) CYCLE
          NPX = NSX(N)
          DXMGN = XMGN(2,NW) - XMGN(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGN(MP,N)*RHOMG(MP,N)
            DFW = TORG(MN,NW)*PORD(MN,NW)*SG(MN,NW)*
     &        DFGN(MN,NW)*RHOMG(MN,NW)
            INDX = 12
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXMGN,INDX )
            UDGN(M,NPX) = DFM*(XMGN(MN,NW)-XMGN(MP,N))/DXGP(NPX)
            FGNP = XGN(MP,N)*RHOG(MP,N)
            FGNW = XGN(MN,NW)*RHOG(MN,NW)
            INDX = 3
            FGN = DIFMN( FGNW,FGNP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            UGN(M,NPX) = UG(M,NPX)*FGN + WTMN*UDGN(M,NPX)
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction vapor mole diffusion, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) CYCLE
          K = KD(N)
          NS = N-IFLD
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) CYCLE
          NPY = NSY(N)
          DXMGN = XMGN(2,NS) - XMGN(2,N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGN(MP,N)*RHOMG(MP,N)
            DFS = TORG(MN,NS)*PORD(MN,NS)*SG(MN,NS)*
     &        DFGN(MN,NS)*RHOMG(MN,NS)
            INDX = 12
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXMGN,INDX )
            VDGN(M,NPY) = DFM*(XMGN(MN,NS)-XMGN(MP,N))/(DYGP(NPY)*RP(I))
            FGNP = XGN(MP,N)*RHOG(MP,N)
            FGNS = XGN(MN,NS)*RHOG(MN,NS)
            INDX = 3
            FGN = DIFMN( FGNS,FGNP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
            VGN(M,NPY) = VG(M,NPY)*FGN + WTMN*VDGN(M,NPY)
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction vapor mole diffusion, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) CYCLE
          NB = N-IJFLD
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) CYCLE
          NPZ = NSZ(N)
          DXMGN = XMGN(2,NB) - XMGN(2,N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGN(MP,N)*RHOMG(MP,N)
            DFB = TORG(MN,NB)*PORD(MN,NB)*SG(MN,NB)*
     &        DFGN(MN,NB)*RHOMG(MN,NB)
            INDX = 12
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXMGN,INDX )
            WDGN(M,NPZ) = DFM*(XMGN(MN,NB)-XMGN(MP,N))/DZGP(NPZ)
            FGNP = XGN(MP,N)*RHOG(MP,N)
            FGNB = XGN(MN,NB)*RHOG(MN,NB)
            INDX = 3
            FGN = DIFMN( FGNB,FGNP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
            WGN(M,NPZ) = WG(M,NPZ)*FGN + WTMN*WDGN(M,NPZ)
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGN group  ---
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
!     Diffusive CH4 gas fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
        DO N = 1,NFLD
          I = ID(N)
          IF( I.EQ.1 ) CYCLE
          J = JD(N)
          K = KD(N)
          NW = N-1
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) CYCLE
          NPX = NSX(N)
          DXMGO = XMGO(2,NW) - XMGO(2,N)
          DO M = 1,ISVF
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
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction vapor mole diffusion, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) CYCLE
          K = KD(N)
          NS = N-IFLD
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) CYCLE
          NPY = NSY(N)
          DXMGO = XMGO(2,NS) - XMGO(2,N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGO(MP,N)*RHOMG(MP,N)
            DFS = TORG(MN,NS)*PORD(MN,NS)*SG(MN,NS)*
     &        DFGO(MN,NS)*RHOMG(MN,NS)
            INDX = 12
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXMGO,INDX )
            VDGO(M,NPY) = DFM*(XMGO(MN,NS)-XMGO(MP,N))/(DYGP(NPY)*RP(I))
            FGOP = XGO(MP,N)*RHOG(MP,N)
            FGOS = XGO(MN,NS)*RHOG(MN,NS)
            INDX = 3
            FGO = DIFMN( FGOS,FGOP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
            VGO(M,NPY) = VG(M,NPY)*FGO + WTMO*VDGO(M,NPY)
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction vapor mole diffusion, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) CYCLE
          NB = N-IJFLD
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) CYCLE
          NPZ = NSZ(N)
          DXMGO = XMGO(2,NB) - XMGO(2,N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGO(MP,N)*RHOMG(MP,N)
            DFB = TORG(MN,NB)*PORD(MN,NB)*SG(MN,NB)*
     &        DFGO(MN,NB)*RHOMG(MN,NB)
            INDX = 12
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXMGO,INDX )
            WDGO(M,NPZ) = DFM*(XMGO(MN,NB)-XMGO(MP,N))/DZGP(NPZ)
            FGOP = XGO(MP,N)*RHOG(MP,N)
            FGOB = XGO(MN,NB)*RHOG(MN,NB)
            INDX = 3
            FGO = DIFMN( FGOB,FGOP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
            WGO(M,NPZ) = WG(M,NPZ)*FGO + WTMO*WDGO(M,NPZ)
          ENDDO
        ENDDO
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
!     Diffusive water gas fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
        DO N = 1,NFLD
          I = ID(N)
          IF( I.EQ.1 ) CYCLE
          J = JD(N)
          K = KD(N)
          NW = N-1
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) CYCLE
          NPX = NSX(N)
          DXMGW = XMGW(2,NW) - XMGW(2,N)
          DO M = 1,ISVF
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
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction vapor mole diffusion, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) CYCLE
          K = KD(N)
          NS = N-IFLD
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) CYCLE
          NPY = NSY(N)
          DXMGW = XMGW(2,NS) - XMGW(2,N)
          DO M = 1,ISVF
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
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction vapor mole diffusion, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) CYCLE
          NB = N-IJFLD
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) CYCLE
          NPZ = NSZ(N)
          DXMGW = XMGW(2,NB) - XMGW(2,N)
          DO M = 1,ISVF
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
          ENDDO
        ENDDO
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
!     Diffusive CO2 gas fluxes on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      DO M = 1,ISVF
        MP = MPOS(M)
        FGAP = XGA(MP,N)*RHOG(MP,N)
        FGAB = XGAB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGA = DIFMN( FGAB,FGAP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WG(1,NPZ).LT.-EPSL ) THEN
            WGA(M,NPZ) = WG(M,NPZ)*FGA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WG(1,NPZ).GT.EPSL ) THEN
            WGA(M,NPZ) = WG(M,NPZ)*FGA
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGA(MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGAB(MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMGA,INDX )
          WDGA(M,NPZ) = 2.D+0*DFM*(XMGAB(MP,NB)-XMGA(MP,N))/DZGF(N)
          WGA(M,NPZ) = WG(M,NPZ)*FGA + WTMA*WDGA(M,NPZ)
        ENDIF
      ENDDO
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
!     Diffusive CO2 gas fluxes on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      DO M = 1,ISVF
        MN = MNEG(M)
        FGAP = XGA(MN,N)*RHOG(MN,N)
        FGAB = XGAB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGA = DIFMN( FGAP,FGAB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UG(1,NQX).GT.EPSL ) THEN
            UGA(M,NQX) = UG(M,NQX)*FGA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UG(1,NQX).LT.-EPSL ) THEN
            UGA(M,NQX) = UG(M,NQX)*FGA
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGA(MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGAB(MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMGA,INDX )
          UDGA(M,NQX) = 2.D+0*DFM*(XMGA(MN,N)-XMGAB(MN,NB))/DXGF(N)
          UGA(M,NQX) = UG(M,NQX)*FGA + WTMA*UDGA(M,NQX)
        ENDIF
      ENDDO
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
!     Diffusive CO2 gas fluxes on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      DO M = 1,ISVF
        MN = MNEG(M)
        FGAP = XGA(MN,N)*RHOG(MN,N)
        FGAB = XGAB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGA = DIFMN( FGAP,FGAB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VG(1,NQY).GT.EPSL ) THEN
            VGA(M,NQY) = VG(M,NQY)*FGA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VG(1,NQY).LT.-EPSL ) THEN
            VGA(M,NQY) = VG(M,NQY)*FGA
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGA(MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGAB(MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMGA,INDX )
          VDGA(M,NQY) = 2.D+0*DFM*(XMGA(MN,N)-XMGAB(MN,NB))
     &      /(DYGF(N)*RP(I))
          VGA(M,NQY) = VG(M,NQY)*FGA + WTMA*VDGA(M,NQY)
        ENDIF
      ENDDO
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
!     Diffusive CO2 gas fluxes on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      DO M = 1,ISVF
        MP = MPOS(M)
        FGAP = XGA(MP,N)*RHOG(MP,N)
        FGAB = XGAB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGA = DIFMN( FGAB,FGAP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VG(1,NPY).LT.-EPSL ) THEN
            VGA(M,NPY) = VG(M,NPY)*FGA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VG(1,NPY).GT.EPSL ) THEN
            VGA(M,NPY) = VG(M,NPY)*FGA
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGA(MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGAB(MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMGA,INDX )
          VDGA(M,NPY) = 2.D+0*DFM*(XMGAB(MP,NB)-XMGA(MP,N))
     &      /(DYGF(N)*RP(I))
          VGA(M,NPY) = VG(M,NPY)*FGA + WTMA*VDGA(M,NPY)
        ENDIF
      ENDDO
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
!     Diffusive CO2 gas fluxes on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      DO M = 1,ISVF
        MN = MNEG(M)
        FGAP = XGA(MN,N)*RHOG(MN,N)
        FGAB = XGAB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGA = DIFMN( FGAP,FGAB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WG(1,NQZ).GT.EPSL ) THEN
            WGA(M,NQZ) = WG(M,NQZ)*FGA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WG(1,NQZ).LT.-EPSL ) THEN
            WGA(M,NQZ) = WG(M,NQZ)*FGA
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGA(MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGAB(MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMGA,INDX )
          WDGA(M,NQZ) = 2.D+0*DFM*(XMGA(MN,N)-XMGAB(MN,NB))/DZGF(N)
          WGA(M,NQZ) = WG(M,NQZ)*FGA + WTMA*WDGA(M,NQZ)
        ENDIF
      ENDDO
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
!     Diffusive CO2 gas fluxes on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      DO M = 1,ISVF
        MP = MPOS(M)
        FGAP = XGA(MP,N)*RHOG(MP,N)
        FGAB = XGAB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGA = DIFMN( FGAB,FGAP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UG(1,NPX).LT.-EPSL ) THEN
            UGA(M,NPX) = UG(M,NPX)*FGA
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UG(1,NPX).GT.EPSL ) THEN
            UGA(M,NPX) = UG(M,NPX)*FGA
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGA(MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGAB(MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMGA,INDX )
          UDGA(M,NPX) = 2.D+0*DFM*(XMGAB(MP,NB)-XMGA(MP,N))/DXGF(N)
          UGA(M,NPX) = UG(M,NPX)*FGA + WTMA*UDGA(M,NPX)
        ENDIF
      ENDDO
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
      SUBROUTINE DFFGNB( N,NB )
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
!     Diffusive N2 gas fluxes on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      USE FDVH
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
      SUB_LOG(ISUB_LOG) = '/DFFGNB'
      K = KD(N)
      NPZ = NSZ(N)
      DXMGN = XMGNB(2,NB) - XMGN(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FGNP = XGN(MP,N)*RHOG(MP,N)
        FGNB = XGNB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGN = DIFMN( FGNB,FGNP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WG(1,NPZ).LT.-EPSL ) THEN
            WGN(M,NPZ) = WG(M,NPZ)*FGN
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WG(1,NPZ).GT.EPSL ) THEN
            WGN(M,NPZ) = WG(M,NPZ)*FGN
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGN(MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGNB(MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMGN,INDX )
          WDGN(M,NPZ) = 2.D+0*DFM*(XMGNB(MP,NB)-XMGN(MP,N))/DZGF(N)
          WGN(M,NPZ) = WG(M,NPZ)*FGN + WTMN*WDGN(M,NPZ)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGNB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGNE( N,NB )
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
!     Diffusive N2 gas fluxes on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      USE FDVH
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
      SUB_LOG(ISUB_LOG) = '/DFFGNE'
      I = ID(N)
      NQX = NSX(N)+1
      DXMGN = XMGN(2,N) - XMGNB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FGNP = XGN(MN,N)*RHOG(MN,N)
        FGNB = XGNB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGN = DIFMN( FGNP,FGNB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UG(1,NQX).GT.EPSL ) THEN
            UGN(M,NQX) = UG(M,NQX)*FGN
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UG(1,NQX).LT.-EPSL ) THEN
            UGN(M,NQX) = UG(M,NQX)*FGN
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGN(MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGNB(MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMGN,INDX )
          UDGN(M,NQX) = 2.D+0*DFM*(XMGN(MN,N)-XMGNB(MN,NB))/DXGF(N)
          UGN(M,NQX) = UG(M,NQX)*FGN + WTMN*UDGN(M,NQX)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGNE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGNN( N,NB )
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
!     Diffusive N2 gas fluxes on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      USE FDVH
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
      SUB_LOG(ISUB_LOG) = '/DFFGNN'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMGN = XMGN(2,N) - XMGNB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FGNP = XGN(MN,N)*RHOG(MN,N)
        FGNB = XGNB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGN = DIFMN( FGNP,FGNB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VG(1,NQY).GT.EPSL ) THEN
            VGN(M,NQY) = VG(M,NQY)*FGN
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VG(1,NQY).LT.-EPSL ) THEN
            VGN(M,NQY) = VG(M,NQY)*FGN
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGN(MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGNB(MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMGN,INDX )
          VDGN(M,NQY) = 2.D+0*DFM*(XMGN(MN,N)-XMGNB(MN,NB))
     &      /(DYGF(N)*RP(I))
          VGN(M,NQY) = VG(M,NQY)*FGN + WTMN*VDGN(M,NQY)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGNN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGNS( N,NB )
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
!     Diffusive N2 gas fluxes on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      USE FDVH
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
      SUB_LOG(ISUB_LOG) = '/DFFGNS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMGN = XMGNB(2,NB) - XMGN(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FGNP = XGN(MP,N)*RHOG(MP,N)
        FGNB = XGNB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGN = DIFMN( FGNB,FGNP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VG(1,NPY).LT.-EPSL ) THEN
            VGN(M,NPY) = VG(M,NPY)*FGN
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VG(1,NPY).GT.EPSL ) THEN
            VGN(M,NPY) = VG(M,NPY)*FGN
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGN(MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGNB(MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMGN,INDX )
          VDGN(M,NPY) = 2.D+0*DFM*(XMGNB(MP,NB)-XMGN(MP,N))
     &      /(DYGF(N)*RP(I))
          VGN(M,NPY) = VG(M,NPY)*FGN + WTMN*VDGN(M,NPY)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGNS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGNT( N,NB )
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
!     Diffusive N2 gas fluxes on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      USE FDVH
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
      SUB_LOG(ISUB_LOG) = '/DFFGNT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMGN = XMGN(2,N)-XMGNB(2,NB)
      DO M = 1,ISVF
        MN = MNEG(M)
        FGNP = XGN(MN,N)*RHOG(MN,N)
        FGNB = XGNB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGN = DIFMN( FGNP,FGNB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WG(1,NQZ).GT.EPSL ) THEN
            WGN(M,NQZ) = WG(M,NQZ)*FGN
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WG(1,NQZ).LT.-EPSL ) THEN
            WGN(M,NQZ) = WG(M,NQZ)*FGN
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGN(MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGNB(MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMGN,INDX )
          WDGN(M,NQZ) = 2.D+0*DFM*(XMGN(MN,N)-XMGNB(MN,NB))/DZGF(N)
          WGN(M,NQZ) = WG(M,NQZ)*FGN + WTMN*WDGN(M,NQZ)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGNT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGNW( N,NB )
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
!     Diffusive N2 gas fluxes on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      USE FDVH
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
      SUB_LOG(ISUB_LOG) = '/DFFGNW'
      I = ID(N)
      NPX = NSX(N)
      DXMGN = XMGNB(2,NB)-XMGN(2,N)
      DO M = 1,ISVF
        MP = MPOS(M)
        FGNP = XGN(MP,N)*RHOG(MP,N)
        FGNB = XGNB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGN = DIFMN( FGNB,FGNP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UG(1,NPX).LT.-EPSL ) THEN
            UGN(M,NPX) = UG(M,NPX)*FGN
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UG(1,NPX).GT.EPSL ) THEN
            UGN(M,NPX) = UG(M,NPX)*FGN
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGN(MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGNB(MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMGN,INDX )
          UDGN(M,NPX) = 2.D+0*DFM*(XMGNB(MP,NB)-XMGN(MP,N))/DXGF(N)
          UGN(M,NPX) = UG(M,NPX)*FGN + WTMN*UDGN(M,NPX)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGNW group  ---
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
!     Diffusive CH4 gas fluxes on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      DO M = 1,ISVF
        MP = MPOS(M)
        FGOP = XGO(MP,N)*RHOG(MP,N)
        FGOB = XGOB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGO = DIFMN( FGOB,FGOP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WG(1,NPZ).LT.-EPSL ) THEN
            WGO(M,NPZ) = WG(M,NPZ)*FGO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WG(1,NPZ).GT.EPSL ) THEN
            WGO(M,NPZ) = WG(M,NPZ)*FGO
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGO(MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGOB(MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMGO,INDX )
          WDGO(M,NPZ) = 2.D+0*DFM*(XMGOB(MP,NB)-XMGO(MP,N))/DZGF(N)
          WGO(M,NPZ) = WG(M,NPZ)*FGO + WTMO*WDGO(M,NPZ)
        ENDIF
      ENDDO
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
!     Diffusive CH4 gas fluxes on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      DO M = 1,ISVF
        MN = MNEG(M)
        FGOP = XGO(MN,N)*RHOG(MN,N)
        FGOB = XGOB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGO = DIFMN( FGOP,FGOB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UG(1,NQX).GT.EPSL ) THEN
            UGO(M,NQX) = UG(M,NQX)*FGO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UG(1,NQX).LT.-EPSL ) THEN
            UGO(M,NQX) = UG(M,NQX)*FGO
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGO(MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGOB(MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMGO,INDX )
          UDGO(M,NQX) = 2.D+0*DFM*(XMGO(MN,N)-XMGOB(MN,NB))/DXGF(N)
          UGO(M,NQX) = UG(M,NQX)*FGO + WTMO*UDGO(M,NQX)
        ENDIF
      ENDDO
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
!     Diffusive CH4 gas fluxes on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      DO M = 1,ISVF
        MN = MNEG(M)
        FGOP = XGO(MN,N)*RHOG(MN,N)
        FGOB = XGOB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGO = DIFMN( FGOP,FGOB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VG(1,NQY).GT.EPSL ) THEN
            VGO(M,NQY) = VG(M,NQY)*FGO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VG(1,NQY).LT.-EPSL ) THEN
            VGO(M,NQY) = VG(M,NQY)*FGO
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGO(MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGOB(MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMGO,INDX )
          VDGO(M,NQY) = 2.D+0*DFM*(XMGO(MN,N)-XMGOB(MN,NB))
     &      /(DYGF(N)*RP(I))
          VGO(M,NQY) = VG(M,NQY)*FGO + WTMO*VDGO(M,NQY)
        ENDIF
      ENDDO
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
!     Diffusive CH4 gas fluxes on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      DO M = 1,ISVF
        MP = MPOS(M)
        FGOP = XGO(MP,N)*RHOG(MP,N)
        FGOB = XGOB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGO = DIFMN( FGOB,FGOP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VG(1,NPY).LT.-EPSL ) THEN
            VGO(M,NPY) = VG(M,NPY)*FGO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VG(1,NPY).GT.EPSL ) THEN
            VGO(M,NPY) = VG(M,NPY)*FGO
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGO(MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGOB(MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMGO,INDX )
          VDGO(M,NPY) = 2.D+0*DFM*(XMGOB(MP,NB)-XMGO(MP,N))
     &      /(DYGF(N)*RP(I))
          VGO(M,NPY) = VG(M,NPY)*FGO + WTMO*VDGO(M,NPY)
        ENDIF
      ENDDO
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
!     Diffusive CH4 gas fluxes on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      DO M = 1,ISVF
        MN = MNEG(M)
        FGOP = XGO(MN,N)*RHOG(MN,N)
        FGOB = XGOB(MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGO = DIFMN( FGOP,FGOB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WG(1,NQZ).GT.EPSL ) THEN
            WGO(M,NQZ) = WG(M,NQZ)*FGO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WG(1,NQZ).LT.-EPSL ) THEN
            WGO(M,NQZ) = WG(M,NQZ)*FGO
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGO(MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGOB(MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMGO,INDX )
          WDGO(M,NQZ) = 2.D+0*DFM*(XMGO(MN,N)-XMGOB(MN,NB))/DZGF(N)
          WGO(M,NQZ) = WG(M,NQZ)*FGO + WTMO*WDGO(M,NQZ)
        ENDIF
      ENDDO
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
!     Diffusive CH4 gas fluxes on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      DO M = 1,ISVF
        MP = MPOS(M)
        FGOP = XGO(MP,N)*RHOG(MP,N)
        FGOB = XGOB(MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGO = DIFMN( FGOB,FGOP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UG(1,NPX).LT.-EPSL ) THEN
            UGO(M,NPX) = UG(M,NPX)*FGO
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UG(1,NPX).GT.EPSL ) THEN
            UGO(M,NPX) = UG(M,NPX)*FGO
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGO(MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGOB(MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMGO,INDX )
          UDGO(M,NPX) = 2.D+0*DFM*(XMGOB(MP,NB)-XMGO(MP,N))/DXGF(N)
          UGO(M,NPX) = UG(M,NPX)*FGO + WTMO*UDGO(M,NPX)
        ENDIF
      ENDDO
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
!     Diffusive water gas fluxes on bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      DO M = 1,ISVF
        MP = MPOS(M)
        FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
        FGWP = XGW(MP,N)*RHOG(MP,N)
        INDX = 3
        FGW = DIFMN( FGWB,FGWP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WG(1,NPZ).LT.-EPSL ) THEN
            WGW(M,NPZ) = WG(M,NPZ)*FGW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WG(1,NPZ).GT.EPSL ) THEN
            WGW(M,NPZ) = WG(M,NPZ)*FGW
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGW(MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGWB(MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMGW,INDX )
          WDGW(M,NPZ) = 2.D+0*DFM*(XMGWB(MP,NB)-XMGW(MP,N))/DZGF(N)
          WGW(M,NPZ) = WG(M,NPZ)*FGW + WTMW*WDGW(M,NPZ)
        ENDIF
      ENDDO
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
!     Diffusive water gas fluxes on east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      DO M = 1,ISVF
        MN = MNEG(M)
        FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
        FGWP = XGW(MN,N)*RHOG(MN,N)
        INDX = 3
        FGW = DIFMN( FGWP,FGWB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UG(1,NQX).GT.EPSL ) THEN
            UGW(M,NQX) = UG(M,NQX)*FGW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UG(1,NQX).LT.-EPSL ) THEN
            UGW(M,NQX) = UG(M,NQX)*FGW
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGW(MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGWB(MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMGW,INDX )
          UDGW(M,NQX) = 2.D+0*DFM*(XMGW(MN,N)-XMGWB(MN,NB))/DXGF(N)
          UGW(M,NQX) = UG(M,NQX)*FGW + WTMW*UDGW(M,NQX)
        ENDIF
      ENDDO
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
!     Diffusive water gas fluxes on north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      DO M = 1,ISVF
        MN = MNEG(M)
        FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
        FGWP = XGW(MN,N)*RHOG(MN,N)
        INDX = 3
        FGW = DIFMN( FGWP,FGWB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VG(1,NQY).GT.EPSL ) THEN
            VGW(M,NQY) = VG(M,NQY)*FGW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VG(1,NQY).LT.-EPSL ) THEN
            VGW(M,NQY) = VG(M,NQY)*FGW
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGW(MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGWB(MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMGW,INDX )
          VDGW(M,NQY) = 2.D+0*DFM*(XMGW(MN,N)-XMGWB(MN,NB))
     &      /(DYGF(N)*RP(I))
          VGW(M,NQY) = VG(M,NQY)*FGW + WTMW*VDGW(M,NQY)
        ENDIF
      ENDDO
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
!     Diffusive water gas fluxes on south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      DO M = 1,ISVF
        MP = MPOS(M)
        FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
        FGWP = XGW(MP,N)*RHOG(MP,N)
        INDX = 3
        FGW = DIFMN( FGWB,FGWP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( VG(1,NPY).LT.-EPSL ) THEN
            VGW(M,NPY) = VG(M,NPY)*FGW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( VG(1,NPY).GT.EPSL ) THEN
            VGW(M,NPY) = VG(M,NPY)*FGW
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGW(MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGWB(MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMGW,INDX )
          VDGW(M,NPY) = 2.D+0*DFM*(XMGWB(MP,NB)-XMGW(MP,N))
     &      /(DYGF(N)*RP(I))
          VGW(M,NPY) = VG(M,NPY)*FGW + WTMW*VDGW(M,NPY)
        ENDIF
      ENDDO
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
!     Diffusive water gas fluxes on top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      DO M = 1,ISVF
        MN = MNEG(M)
        FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
        FGWP = XGW(MN,N)*RHOG(MN,N)
        INDX = 3
        FGW = DIFMN( FGWP,FGWB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( WG(1,NQZ).GT.EPSL ) THEN
            WGW(M,NQZ) = WG(M,NQZ)*FGW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( WG(1,NQZ).LT.-EPSL ) THEN
            WGW(M,NQZ) = WG(M,NQZ)*FGW
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGW(MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGWB(MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMGW,INDX )
          WDGW(M,NQZ) = 2.D+0*DFM*(XMGW(MN,N)-XMGWB(MN,NB))/DZGF(N)
          WGW(M,NQZ) = WG(M,NQZ)*FGW + WTMW*WDGW(M,NQZ)
        ENDIF
      ENDDO
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
!     Diffusive water gas fluxes on west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 June 2012.
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
      DO M = 1,ISVF
        MP = MPOS(M)
        FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
        FGWP = XGW(MP,N)*RHOG(MP,N)
        INDX = 3
        FGW = DIFMN( FGWB,FGWP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---    Dirichlet-Outflow boundary condition  ---
!
        IF( IBCT(2,NB).EQ.15 .OR. IBCT(2,NB).EQ.26 .OR. 
     &    IBCT(2,NB).EQ.30 ) THEN
          IF( UG(1,NPX).LT.-EPSL ) THEN
            UGW(M,NPX) = UG(M,NPX)*FGW
          ENDIF
!
!---    Dirichlet-Inflow boundary condition  ---
!
        ELSEIF( IBCT(2,NB).EQ.14 .OR. IBCT(2,NB).EQ.29 .OR. 
     &    MOD(IBCT(2,NB),10).EQ.7 ) THEN
          IF( UG(1,NPX).GT.EPSL ) THEN
            UGW(M,NPX) = UG(M,NPX)*FGW
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGW(MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGWB(MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMGW,INDX )
          UDGW(M,NPX) = 2.D+0*DFM*(XMGWB(MP,NB)-XMGW(MP,N))/DXGF(N)
          UGW(M,NPX) = UG(M,NPX)*FGW + WTMW*UDGW(M,NPX)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGWW group  ---
!
      RETURN
      END


!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGAW
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
!     Diffusive CO2 and water gas fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE HYST
      USE GRID
      USE FLUX
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
      SUB_LOG(ISUB_LOG) = '/DFFGAW'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(28).EQ.1 ) THEN
!
!---  Gas water diffusive flux, excluding boundaries
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
!
!---    West surface  ---
!
        NW = ICM(3,N)
        IF( NW.NE.0 ) THEN
          DXMGW = XMGW(2,NW) - XMGW(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &        DFGW(MP,N)*RHOMG(MP,N)
            DFW = TORG(MN,NW)*PORD(MN,NW)*(SG(MN,NW)-SGT(MN,NW))*
     &        DFGW(MN,NW)*RHOMG(MN,NW)
            INDX = 12
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXMGW,INDX )
            UDGW(M,1,N) = DFM*(XMGW(MN,NW)-XMGW(MP,N))/DXGP(1,N)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            FGWW = XGW(MN,NW)*RHOG(MN,NW)
            INDX = 3
            FGW = DIFMN( FGWW,FGWP,DXGF(NW),DXGF(N),UG(1,1,N),INDX )
            UGW(M,1,N) = UG(M,1,N)*FGW + WTMW*UDGW(M,1,N)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            FGAW = XGA(MN,NW)*RHOG(MN,NW)
            INDX = 3
            FGA = DIFMN( FGAW,FGAP,DXGF(NW),DXGF(N),UG(1,1,N),INDX )
            UGA(M,1,N) = UG(M,1,N)*FGA - WTMA*UDGW(M,1,N)
          ENDDO
          DO M = 1,ISVF
            UDGW(M,2,NW) = UDGW(M,1,N)
            UGW(M,2,NW) = UGW(M,1,N)
            UGA(M,2,NW) = UGA(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            UDGW(M,1,N) = 0.D+0
            UGW(M,1,N) = 0.D+0
            UGA(M,1,N) = 0.D+0
          ENDDO
        ENDIF
!
!---    South surface  ---
!
        NS = ICM(2,N)
        IF( NS.NE.0 ) THEN
          DXMGW = XMGW(2,NS) - XMGW(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &        DFGW(MP,N)*RHOMG(MP,N)
            DFS = TORG(MN,NS)*PORD(MN,NS)*(SG(MN,NS)-SGT(MN,NS))*
     &        DFGW(MN,NS)*RHOMG(MN,NS)
            INDX = 12
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXMGW,INDX )
            VDGW(M,1,N) = DFM*(XMGW(MN,NS)-XMGW(MP,N))/
     &        (DYGP(1,N)*RP(N))
            FGWP = XGW(MP,N)*RHOG(MP,N)
            FGWS = XGW(MN,NS)*RHOG(MN,NS)
            INDX = 3
            FGW = DIFMN( FGWS,FGWP,DYGF(NS),DYGF(N),VG(1,1,N),INDX )
            VGW(M,1,N) = VG(M,1,N)*FGW + WTMW*VDGW(M,1,N)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            FGAS = XGA(MN,NS)*RHOG(MN,NS)
            INDX = 3
            FGA = DIFMN( FGAS,FGAP,DYGF(NS),DYGF(N),VG(1,1,N),INDX )
            VGA(M,1,N) = VG(M,1,N)*FGA - WTMA*VDGW(M,1,N)
          ENDDO
          DO M = 1,ISVF
            VDGW(M,2,NS) = VDGW(M,1,N)
            VGW(M,2,NS) = VGW(M,1,N)
            VGA(M,2,NS) = VGA(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            VDGW(M,1,N) = 0.D+0
            VGW(M,1,N) = 0.D+0
            VGA(M,1,N) = 0.D+0
          ENDDO
        ENDIF
!
!---    Bottom surface  ---
!
        NB = ICM(1,N)
        IF( NB.NE.0 ) THEN
          DXMGW = XMGW(2,NB) - XMGW(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &        DFGW(MP,N)*RHOMG(MP,N)
            DFB = TORG(MN,NB)*PORD(MN,NB)*(SG(MN,NB)-SGT(MN,NB))*
     &        DFGW(MN,NB)*RHOMG(MN,NB)
            INDX = 12
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXMGW,INDX )
            WDGW(M,1,N) = DFM*(XMGW(MN,NB)-XMGW(MP,N))/DZGP(1,N)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            FGWB = XGW(MN,NB)*RHOG(MN,NB)
            INDX = 3
            FGW = DIFMN( FGWB,FGWP,DZGF(NB),DZGF(N),WG(1,1,N),INDX )
            WGW(M,1,N) = WG(M,1,N)*FGW + WTMW*WDGW(M,1,N)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            FGAB = XGA(MN,NB)*RHOG(MN,NB)
            INDX = 3
            FGA = DIFMN( FGAB,FGAP,DZGF(NB),DZGF(N),WG(1,1,N),INDX )
            WGA(M,1,N) = WG(M,1,N)*FGA - WTMA*WDGW(M,1,N)
          ENDDO
          DO M = 1,ISVF
            WDGW(M,2,NB) = WDGW(M,1,N)
            WGW(M,2,NB) = WGW(M,1,N)
            WGA(M,2,NB) = WGA(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            WDGW(M,1,N) = 0.D+0
            WGW(M,1,N) = 0.D+0
            WGA(M,1,N) = 0.D+0
          ENDDO
        ENDIF
      ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
!
!---  Gas water diffusive flux, excluding boundaries
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
!
!---    West surface  ---
!
        NW = ICM(3,N)
        IF( NW.NE.0 ) THEN
          DXMGW = XMGW(2,NW)*RHOMG(2,NW) - XMGW(2,N)*RHOMG(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &        DFGW(MP,N)
            DFW = TORG(MN,NW)*PORD(MN,NW)*(SG(MN,NW)-SGT(MN,NW))*
     &        DFGW(MN,NW)
            INDX = 12
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXMGW,INDX )
            UDGW(M,1,N) = DFM*(XMGW(MN,NW)*RHOMG(MN,NW) -
     &        XMGW(MP,N)*RHOMG(MP,N))/DXGP(1,N)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            FGWW = XGW(MN,NW)*RHOG(MN,NW)
            INDX = 3
            FGW = DIFMN( FGWW,FGWP,DXGF(NW),DXGF(N),UG(1,1,N),INDX )
            UGW(M,1,N) = UG(M,1,N)*FGW + WTMW*UDGW(M,1,N)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            FGAW = XGA(MN,NW)*RHOG(MN,NW)
            INDX = 3
            FGA = DIFMN( FGAW,FGAP,DXGF(NW),DXGF(N),UG(1,1,N),INDX )
            UGA(M,1,N) = UG(M,1,N)*FGA - WTMA*UDGW(M,1,N)
          ENDDO
          DO M = 1,ISVF
            UDGW(M,2,NW) = UDGW(M,1,N)
            UGW(M,2,NW) = UGW(M,1,N)
            UGA(M,2,NW) = UGA(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            UDGW(M,1,N) = 0.D+0
            UGW(M,1,N) = 0.D+0
            UGA(M,1,N) = 0.D+0
          ENDDO
        ENDIF
!
!---    South surface  ---
!
        NS = ICM(2,N)
        IF( NS.NE.0 ) THEN
          DXMGW = XMGW(2,NS)*RHOMG(2,NS) - XMGW(2,N)*RHOMG(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &        DFGW(MP,N)
            DFS = TORG(MN,NS)*PORD(MN,NS)*(SG(MN,NS)-SGT(MN,NS))*
     &        DFGW(MN,NS)
            INDX = 12
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXMGW,INDX )
            VDGW(M,1,N) = DFM*(XMGW(MN,NS)*RHOMG(MN,NS) -
     &        XMGW(MP,N)*RHOMG(MP,N))/(DYGP(1,N)*RP(N))
            FGWP = XGW(MP,N)*RHOG(MP,N)
            FGWS = XGW(MN,NS)*RHOG(MN,NS)
            INDX = 3
            FGW = DIFMN( FGWS,FGWP,DYGF(NS),DYGF(N),VG(1,1,N),INDX )
            VGW(M,1,N) = VG(M,1,N)*FGW + WTMW*VDGW(M,1,N)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            FGAS = XGA(MN,NS)*RHOG(MN,NS)
            INDX = 3
            FGA = DIFMN( FGAS,FGAP,DYGF(NS),DYGF(N),VG(1,1,N),INDX )
            VGA(M,1,N) = VG(M,1,N)*FGA - WTMA*VDGW(M,1,N)
          ENDDO
          DO M = 1,ISVF
            VDGW(M,2,NS) = VDGW(M,1,N)
            VGW(M,2,NS) = VGW(M,1,N)
            VGA(M,2,NS) = VGA(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            VDGW(M,1,N) = 0.D+0
            VGW(M,1,N) = 0.D+0
            VGA(M,1,N) = 0.D+0
          ENDDO
        ENDIF
!
!---    Bottom surface  ---
!
        NB = ICM(1,N)
        IF( NB.NE.0 ) THEN
          DXMGW = XMGW(2,NB)*RHOMG(2,NB) - XMGW(2,N)*RHOMG(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORG(MP,N)*PORD(MP,N)*(SG(MP,N)-SGT(MP,N))*
     &        DFGW(MP,N)
            DFB = TORG(MN,NB)*PORD(MN,NB)*(SG(MN,NB)-SGT(MN,NB))*
     &        DFGW(MN,NB)
            INDX = 12
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXMGW,INDX )
            WDGW(M,1,N) = DFM*(XMGW(MN,NB)*RHOMG(MN,NB) -
     &        XMGW(MP,N)*RHOMG(MP,N))/DZGP(1,N)
            FGWP = XGW(MP,N)*RHOG(MP,N)
            FGWB = XGW(MN,NB)*RHOG(MN,NB)
            INDX = 3
            FGW = DIFMN( FGWB,FGWP,DZGF(NB),DZGF(N),WG(1,1,N),INDX )
            WGW(M,1,N) = WG(M,1,N)*FGW + WTMW*WDGW(M,1,N)
            FGAP = XGA(MP,N)*RHOG(MP,N)
            FGAB = XGA(MN,NB)*RHOG(MN,NB)
            INDX = 3
            FGA = DIFMN( FGAB,FGAP,DZGF(NB),DZGF(N),WG(1,1,N),INDX )
            WGA(M,1,N) = WG(M,1,N)*FGA - WTMA*WDGW(M,1,N)
          ENDDO
          DO M = 1,ISVF
            WDGW(M,2,NB) = WDGW(M,1,N)
            WGW(M,2,NB) = WGW(M,1,N)
            WGA(M,2,NB) = WGA(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            WDGW(M,1,N) = 0.D+0
            WGW(M,1,N) = 0.D+0
            WGA(M,1,N) = 0.D+0
          ENDDO
        ENDIF
      ENDDO
      ENDIF
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
      SUBROUTINE DFFLAW
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
!     Diffusive CO2, water and salt aqueous fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE GRID
      USE FLUX
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
      SUB_LOG(ISUB_LOG) = '/DFFLAW'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
!
!---  Aqueous diffusive flux, excluding boundaries
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
!
!---    West surface  ---
!
        NW = ICM(3,N)
        IF( NW.NE.0 ) THEN
          DXMLA = XMLA(2,NW) - XMLA(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*
     &        DFLA(MP,N)*RHOML(MP,N)
            DFW = TORL(MN,NW)*PORD(MN,NW)*SL(MN,NW)*
     &        DFLA(MN,NW)*RHOML(MN,NW)
            INDX = 14
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXMLA,INDX)
            UDLA(M,1,N) = DFM*(XMLA(MN,NW)-XMLA(MP,N))/DXGP(1,N)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            FLAW = XLA(MN,NW)*RHOL(MN,NW)
            INDX = 2
            FLA = DIFMN( FLAW,FLAP,DXGF(NW),DXGF(N),UL(1,1,N),INDX )
            ULA(M,1,N) = UL(M,1,N)*FLA + WTMA*UDLA(M,1,N)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            FLWW = XLW(MN,NW)*RHOL(MN,NW)
            INDX = 2
            FLW = DIFMN( FLWW,FLWP,DXGF(NW),DXGF(N),UL(1,1,N),INDX )
            ULW(M,1,N) = UL(M,1,N)*FLW - WTMW*(UDLA(M,1,N) +
     &        UDS(M,1,N)/WTMS)
          ENDDO
          DO M = 1,ISVF
            UDLA(M,2,NW) = UDLA(M,1,N)
            ULA(M,2,NW) = ULA(M,1,N)
            ULW(M,2,NW) = ULW(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            UDLA(M,1,N) = 0.D+0
            ULA(M,1,N) = 0.D+0
            ULW(M,1,N) = 0.D+0
          ENDDO
        ENDIF
!
!---    South surface  ---
!
        NS = ICM(2,N)
        IF( NS.NE.0 ) THEN
          DXMLA = XMLA(2,NS) - XMLA(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*
     &        DFLA(MP,N)*RHOML(MP,N)
            DFS = TORL(MN,NS)*PORD(MN,NS)*SL(MN,NS)*
     &        DFLA(MN,NS)*RHOML(MN,NS)
            INDX = 14
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXMLA,INDX )
            VDLA(M,1,N) = DFM*(XMLA(MN,NS)-XMLA(MP,N))/(DYGP(1,N)*RP(N))
            FLAP = XLA(MP,N)*RHOL(MP,N)
            FLAS = XLA(MN,NS)*RHOL(MN,NS)
            INDX = 2
            FLA = DIFMN( FLAS,FLAP,DYGF(NS),DYGF(N),VL(1,1,N),INDX )
            VLA(M,1,N) = VL(M,1,N)*FLA + WTMA*VDLA(M,1,N)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            FLWS = XLW(MN,NS)*RHOL(MN,NS)
            INDX = 2
            FLW = DIFMN( FLWS,FLWP,DYGF(NS),DYGF(N),VL(1,1,N),INDX )
            VLW(M,1,N) = VL(M,1,N)*FLW - WTMW*(VDLA(M,1,N) +
     &        VDS(M,1,N)/WTMS)
          ENDDO
          DO M = 1,ISVF
            VDLA(M,2,NS) = VDLA(M,1,N)
            VLA(M,2,NS) = VLA(M,1,N)
            VLW(M,2,NS) = VLW(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            VDLA(M,1,N) = 0.D+0
            VLA(M,1,N) = 0.D+0
            VLW(M,1,N) = 0.D+0
          ENDDO
        ENDIF
!
!---    Bottom surface  ---
!
        NB = ICM(1,N)
        IF( NB.NE.0 ) THEN
          DXMLA = XMLA(2,NB) - XMLA(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*
     &        DFLA(MP,N)*RHOML(MP,N)
            DFB = TORL(MN,NB)*PORD(MN,NB)*SL(MN,NB)*
     &        DFLA(MN,NB)*RHOML(MN,NB)
            INDX = 14
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXMLA,INDX )
            WDLA(M,1,N) = DFM*(XMLA(MN,NB)-XMLA(MP,N))/DZGP(1,N)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            FLAB = XLA(MN,NB)*RHOL(MN,NB)
            INDX = 2
            FLA = DIFMN( FLAB,FLAP,DZGF(NB),DZGF(N),WL(1,1,N),INDX )
            WLA(M,1,N) = WL(M,1,N)*FLA + WTMA*WDLA(M,1,N)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            FLWB = XLW(MN,NB)*RHOL(MN,NB)
            INDX = 2
            FLW = DIFMN( FLWB,FLWP,DZGF(NB),DZGF(N),WL(1,1,N),INDX )
            WLW(M,1,N) = WL(M,1,N)*FLW - WTMW*(WDLA(M,1,N) +
     &        WDS(M,1,N)/WTMS)
          ENDDO
          DO M = 1,ISVF
            WDLA(M,2,NB) = WDLA(M,1,N)
            WLA(M,2,NB) = WLA(M,1,N)
            WLW(M,2,NB) = WLW(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            WDLA(M,1,N) = 0.D+0
            WLA(M,1,N) = 0.D+0
            WLW(M,1,N) = 0.D+0
          ENDDO
        ENDIF
      ENDDO
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
!
!---  Aqueous diffusive flux, excluding boundaries
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
!
!---    West surface  ---
!
        NW = ICM(3,N)
        IF( NW.NE.0 ) THEN
          DXMLA = XMLA(2,NW)*RHOML(2,NW) - XMLA(2,N)*RHOML(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*
     &        DFLA(MP,N)
            DFW = TORL(MN,NW)*PORD(MN,NW)*SL(MN,NW)*
     &        DFLA(MN,NW)
            INDX = 14
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXMLA,INDX)
            UDLA(M,1,N) = DFM*(XMLA(MN,NW)*RHOML(MN,NW) -
     &        XMLA(MP,N)*RHOML(MP,N))/DXGP(1,N)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            FLAW = XLA(MN,NW)*RHOL(MN,NW)
            INDX = 2
            FLA = DIFMN( FLAW,FLAP,DXGF(NW),DXGF(N),UL(1,1,N),INDX )
            ULA(M,1,N) = UL(M,1,N)*FLA + WTMA*UDLA(M,1,N)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            FLWW = XLW(MN,NW)*RHOL(MN,NW)
            INDX = 2
            FLW = DIFMN( FLWW,FLWP,DXGF(NW),DXGF(N),UL(1,1,N),INDX )
            ULW(M,1,N) = UL(M,1,N)*FLW - WTMW*(UDLA(M,1,N) +
     &        UDS(M,1,N)/WTMS)
          ENDDO
          DO M = 1,ISVF
            UDLA(M,2,NW) = UDLA(M,1,N)
            ULA(M,2,NW) = ULA(M,1,N)
            ULW(M,2,NW) = ULW(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            UDLA(M,1,N) = 0.D+0
            ULA(M,1,N) = 0.D+0
            ULW(M,1,N) = 0.D+0
          ENDDO
        ENDIF
!
!---    South surface  ---
!
        NS = ICM(2,N)
        IF( NS.NE.0 ) THEN
          DXMLA = XMLA(2,NS)*RHOML(2,NS) - XMLA(2,N)*RHOML(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*
     &        DFLA(MP,N)
            DFS = TORL(MN,NS)*PORD(MN,NS)*SL(MN,NS)*
     &        DFLA(MN,NS)
            INDX = 14
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXMLA,INDX )
            VDLA(M,1,N) = DFM*(XMLA(MN,NS)*RHOML(MN,NS) -
     &        XMLA(MP,N)*RHOML(MP,N))/(DYGP(1,N)*RP(N))
            FLAP = XLA(MP,N)*RHOL(MP,N)
            FLAS = XLA(MN,NS)*RHOL(MN,NS)
            INDX = 2
            FLA = DIFMN( FLAS,FLAP,DYGF(NS),DYGF(N),VL(1,1,N),INDX )
            VLA(M,1,N) = VL(M,1,N)*FLA + WTMA*VDLA(M,1,N)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            FLWS = XLW(MN,NS)*RHOL(MN,NS)
            INDX = 2
            FLW = DIFMN( FLWS,FLWP,DYGF(NS),DYGF(N),VL(1,1,N),INDX )
            VLW(M,1,N) = VL(M,1,N)*FLW - WTMW*(VDLA(M,1,N) +
     &        VDS(M,1,N)/WTMS)
          ENDDO
          DO M = 1,ISVF
            VDLA(M,2,NS) = VDLA(M,1,N)
            VLA(M,2,NS) = VLA(M,1,N)
            VLW(M,2,NS) = VLW(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            VDLA(M,1,N) = 0.D+0
            VLA(M,1,N) = 0.D+0
            VLW(M,1,N) = 0.D+0
          ENDDO
        ENDIF
!
!---    Bottom surface  ---
!
        NB = ICM(1,N)
        IF( NB.NE.0 ) THEN
          DXMLA = XMLA(2,NB)*RHOML(2,NB) - XMLA(2,N)*RHOML(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*
     &        DFLA(MP,N)
            DFB = TORL(MN,NB)*PORD(MN,NB)*SL(MN,NB)*
     &        DFLA(MN,NB)
            INDX = 14
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXMLA,INDX )
            WDLA(M,1,N) = DFM*(XMLA(MN,NB)*RHOML(MN,NB) -
     &        XMLA(MP,N)*RHOML(MP,N))/DZGP(1,N)
            FLAP = XLA(MP,N)*RHOL(MP,N)
            FLAB = XLA(MN,NB)*RHOL(MN,NB)
            INDX = 2
            FLA = DIFMN( FLAB,FLAP,DZGF(NB),DZGF(N),WL(1,1,N),INDX )
            WLA(M,1,N) = WL(M,1,N)*FLA + WTMA*WDLA(M,1,N)
            FLWP = XLW(MP,N)*RHOL(MP,N)
            FLWB = XLW(MN,NB)*RHOL(MN,NB)
            INDX = 2
            FLW = DIFMN( FLWB,FLWP,DZGF(NB),DZGF(N),WL(1,1,N),INDX )
            WLW(M,1,N) = WL(M,1,N)*FLW - WTMW*(WDLA(M,1,N) +
     &        WDS(M,1,N)/WTMS)
          ENDDO
          DO M = 1,ISVF
            WDLA(M,2,NB) = WDLA(M,1,N)
            WLA(M,2,NB) = WLA(M,1,N)
            WLW(M,2,NB) = WLW(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            WDLA(M,1,N) = 0.D+0
            WLA(M,1,N) = 0.D+0
            WLW(M,1,N) = 0.D+0
          ENDDO
        ENDIF
      ENDDO
      ENDIF
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
!     Diffusive salt/inhibitor aqueous fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE GRID
      USE FLUX
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
      SUB_LOG(ISUB_LOG) = '/DFFLS'
!
!---  Aqueous salt/inhibitor diffusive flux, excluding boundaries
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
!
!---    West surface  ---
!
        NW = ICM(3,N)
        IF( NW.NE.0 ) THEN
          DXMLS = XMLS(2,NW) - XMLS(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
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
            ELSEIF( IEDLS.EQ.3 ) THEN
              DFCLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
              DFCLW = TORL(MN,NW)*SL(MN,NW)*PORD(MN,NW)*DFLS(MN,NW)
            ENDIF
            INDX = 18
            DFCLW = DIFMN(DFCLW,DFCLP,DXGF(NW),DXGF(N),UL(1,1,N),INDX)
!
!---        Salt aqueous flux by advection and diffusion  ---
!
            DDLW = DFCLW/DXGP(1,N)
            AL = MAX( UL(M,1,N),ZERO ) +
     &        DDLW*MAX((ONE-(TENTH*ABS(UL(M,1,N))/
     &        (DDLW+SMALL)))**5,ZERO)
            ALP = MAX( -UL(M,1,N),ZERO ) +
     &        DDLW*MAX((ONE-(TENTH*ABS(UL(M,1,N))/
     &        (DDLW+SMALL)))**5,ZERO)
            US(M,1,N) = XLS(MN,NW)*RHOL(MN,NW)*AL -
     &        XLS(MP,N)*RHOL(MP,N)*ALP
            UDS(M,1,N) = DDLW*(XLS(MN,NW)*RHOL(MN,NW) -
     &        XLS(MP,N)*RHOL(MP,N))
          ENDDO
          DO M = 1,ISVF
            UDS(M,2,NW) = UDS(M,1,N)
            US(M,2,NW) = US(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            UDS(M,1,N) = 0.D+0
            US(M,1,N) = 0.D+0
          ENDDO
        ENDIF
!
!---    South surface  ---
!
        NS = ICM(2,N)
        IF( NS.NE.0 ) THEN
          DXMLS = XMLS(2,NS) - XMLS(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
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
            ELSEIF( IEDLS.EQ.3 ) THEN
              DFCLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
              DFCLS = TORL(MN,NS)*SL(MN,NS)*PORD(MN,NS)*DFLS(MN,NS)
            ENDIF
            INDX = 18
            DFCLS = DIFMN(DFCLS,DFCLP,DYGF(NS),DYGF(N),VL(1,1,N),INDX)
!
!---        Salt aqueous flux by advection and diffusion  ---
!
            DDLS = DFCLS/(DYGP(1,N)*RP(N))
            AL = MAX( VL(M,1,N),ZERO ) +
     &        DDLS*MAX((ONE-(TENTH*ABS(VL(M,1,N))/
     &        (DDLS+SMALL)))**5,ZERO)
            ALP = MAX( -VL(M,1,N),ZERO ) +
     &        DDLS*MAX((ONE-(TENTH*ABS(VL(M,1,N))/
     &        (DDLS+SMALL)))**5,ZERO)
            VS(M,1,N) = (XLS(MN,NS)*RHOL(MN,NS)*AL -
     &        XLS(MP,N)*RHOL(MP,N)*ALP)
            VDS(M,1,N) = DDLS*(XLS(MN,NS)*RHOL(MN,NS) -
     &        XLS(MP,N)*RHOL(MP,N))
          ENDDO
          DO M = 1,ISVF
            VDS(M,2,NS) = VDS(M,1,N)
            VS(M,2,NS) = VS(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            VDS(M,1,N) = 0.D+0
            VS(M,1,N) = 0.D+0
          ENDDO
        ENDIF
!
!---    Bottom surface  ---
!
        NB = ICM(1,N)
        IF( NB.NE.0 ) THEN
          DXMLS = XMLS(2,NB) - XMLS(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
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
            ELSEIF( IEDLS.EQ.3 ) THEN
              DFCLP = TORL(MP,N)*SL(MP,N)*PORD(MP,N)*DFLS(MP,N)
              DFCLB = TORL(MN,NB)*SL(MN,NB)*PORD(MN,NB)*DFLS(MN,NB)
            ENDIF
            INDX = 18
            DFCLB = DIFMN(DFCLB,DFCLP,DZGF(NB),DZGF(N),WL(1,1,N),INDX)
!
!---        Salt aqueous flux by advection and diffusion  ---
!
            DDLB = DFCLB/DZGP(1,N)
            AL = MAX( WL(M,1,N),ZERO ) +
     &        DDLB*MAX((ONE-(TENTH*ABS(WL(M,1,N))/
     &        (DDLB+SMALL)))**5,ZERO)
            ALP = MAX( -WL(M,1,N),ZERO ) +
     &        DDLB*MAX((ONE-(TENTH*ABS(WL(M,1,N))/
     &        (DDLB+SMALL)))**5,ZERO)
            WS(M,1,N) = (XLS(MN,NB)*RHOL(MN,NB)*AL -
     &        XLS(MP,N)*RHOL(MP,N)*ALP)
            WDS(M,1,N) = DDLB*(XLS(MN,NB)*RHOL(MN,NB) -
     &        XLS(MP,N)*RHOL(MP,N))
          ENDDO
          DO M = 1,ISVF
            WDS(M,2,NB) = WDS(M,1,N)
            WS(M,2,NB) = WS(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            WDS(M,1,N) = 0.D+0
            WS(M,1,N) = 0.D+0
          ENDDO
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVG
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
!     Compute the gas-phase Darcy flux from pressure gradients
!     and gravitational body forces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE GRID
      USE FLUX
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
      REAL*8 KGM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DRCVG'
!
!---  Gas Darcy velocities, excluding boundaries
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
!
!---    West surface  ---
!
        NW = ICM(3,N)
        IF( NW.NE.0 ) THEN
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            HDGX = PSO(MN,NW)-PSO(MP,N)-0.5D+0*GRVX(1,N)
     &        *(RHOG(MN,NW)*DXGF(N)+RHOG(MP,N)*DXGF(NW))
            IF( M.EQ.1 ) HDG = HDGX
!
!---        Permeability reduction factor
!
            PERM_WX = PERMRF(MN,NW)*PERM(1,NW)
            PERM_PX = PERMRF(MP,N)*PERM(1,N)
            INDX = 11
            KGM = DIFMN(PERM_WX,PERM_PX,DXGF(NW),DXGF(N),HDG,INDX)
            IF( PERM_WX/EPSL.LT.EPSL ) KGM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KGM = 0.D+0
            INDX = 9
            RKGM = DIFMN(RKG(MN,NW),RKG(MP,N),DXGF(NW),DXGF(N),
     &        HDG,INDX)
            INDX = 6
            VGM = DIFMN(VISG(MN,NW),VISG(MP,N),DXGF(NW),DXGF(N),
     &        HDG,INDX)
            UG(M,1,N) = KGM*RKGM*HDGX/(DXGP(1,N)*VGM)
          ENDDO
          DO M = 1,ISVF
            UG(M,2,NW) = UG(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            UG(M,1,N) = 0.D+0
          ENDDO
        ENDIF
!
!---    South surface  ---
!
        NS = ICM(2,N)
        IF( NS.NE.0 ) THEN
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            HDGY = PSO(MN,NS)-PSO(MP,N)-0.5D+0*GRVY(1,N)
     &        *(RHOG(MN,NS)*DYGF(N)+RHOG(MP,N)*DYGF(NS))
            IF( M.EQ.1 ) HDG = HDGY
!
!---        Permeability reduction factor
!
            PERM_SX = PERMRF(MN,NS)*PERM(2,NS)
            PERM_PX = PERMRF(MP,N)*PERM(2,N)
            INDX = 11
            KGM = DIFMN(PERM_SX,PERM_PX,DYGF(NS),DYGF(N),HDG,INDX)
            IF( PERM_SX/EPSL.LT.EPSL ) KGM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KGM = 0.D+0
            INDX = 9
            RKGM = DIFMN(RKG(MN,NS),RKG(MP,N),DYGF(NS),DYGF(N),
     &        HDG,INDX)
            INDX = 6
            VGM = DIFMN(VISG(MN,NS),VISG(MP,N),DYGF(NS),DYGF(N),
     &        HDG,INDX)
            VG(M,1,N) = KGM*RKGM*HDGY/(DYGP(1,N)*VGM*RP(N))
          ENDDO
          DO M = 1,ISVF
            VG(M,2,NS) = VG(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            VG(M,1,N) = 0.D+0
          ENDDO
        ENDIF
!
!---    Bottom surface  ---
!
        NB = ICM(1,N)
        IF( NB.NE.0 ) THEN
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            HDGZ = PSO(MN,NB)-PSO(MP,N)-0.5D+0*GRVZ(1,N)
     &        *(RHOG(MN,NB)*DZGF(N)+RHOG(MP,N)*DZGF(NB))
            IF( M.EQ.1 ) HDG = HDGZ
!
!---        Permeability reduction factor
!
            PERM_BX = PERMRF(MN,NB)*PERM(3,NB)
            PERM_PX = PERMRF(MP,N)*PERM(3,N)
            INDX = 11
            KGM = DIFMN(PERM_BX,PERM_PX,DZGF(NB),DZGF(N),HDG,INDX)
            IF( PERM_BX/EPSL.LT.EPSL ) KGM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KGM = 0.D+0
            INDX = 9
            RKGM = DIFMN(RKG(MN,NB),RKG(MP,N),DZGF(NB),DZGF(N),
     &        HDG,INDX)
            INDX = 6
            VGM = DIFMN(VISG(MN,NB),VISG(MP,N),DZGF(NB),DZGF(N),
     &        HDG,INDX)
            WG(M,1,N) = KGM*RKGM*HDGZ/(DZGP(1,N)*VGM)
          ENDDO
          DO M = 1,ISVF
            WG(M,2,NB) = WG(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            WG(M,1,N) = 0.D+0
          ENDDO
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVG group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DRCVL
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
!     Compute the aqueous-phase Darcy flux from pressure gradients
!     and gravitational body forces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE GRID
      USE FLUX
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
      REAL*8 KLM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DRCVL'
!
!---  Aqueous Darcy velocities, excluding boundaries
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive nodes  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
!
!---    West surface  ---
!
        NW = ICM(3,N)
        IF( NW.NE.0 ) THEN
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            HDLX = PL(MN,NW)-PL(MP,N)-0.5D+0*GRVX(1,N)
     &        *(RHOL(MN,NW)*DXGF(N)+RHOL(MP,N)*DXGF(NW))
            IF( M.EQ.1 ) HDL = HDLX
!
!---        Permeability reduction factor
!
            PERM_WX = PERMRF(MN,NW)*PERM(1,NW)
            PERM_PX = PERMRF(MP,N)*PERM(1,N)
            INDX = 11
            KLM = DIFMN(PERM_WX,PERM_PX,DXGF(NW),DXGF(N),HDL,INDX)
            IF( PERM_WX/EPSL.LT.EPSL ) KLM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KLM = 0.D+0
            INDX = 8
            RKLM = DIFMN(RKL(MN,NW),RKL(MP,N),DXGF(NW),DXGF(N),
     &        HDL,INDX)
            INDX = 5
            VLM = DIFMN(VISL(MN,NW),VISL(MP,N),DXGF(NW),DXGF(N),
     &        HDL,INDX)
            UL(M,1,N) = KLM*RKLM*HDLX/(DXGP(1,N)*VLM)
          ENDDO
          DO M = 1,ISVF
            UL(M,2,NW) = UL(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            UL(M,1,N) = 0.D+0
          ENDDO
        ENDIF
!
!---    South surface  ---
!
        NS = ICM(2,N)
        IF( NS.NE.0 ) THEN
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            HDLY = PL(MN,NS)-PL(MP,N)-0.5D+0*GRVY(1,N)
     &        *(RHOL(MN,NS)*DYGF(N)+RHOL(MP,N)*DYGF(NS))
            IF( M.EQ.1 ) HDL = HDLY
!
!---        Permeability reduction factor
!
            PERM_SX = PERMRF(MN,NS)*PERM(2,NS)
            PERM_PX = PERMRF(MP,N)*PERM(2,N)
            INDX = 11
            KLM = DIFMN(PERM_SX,PERM_PX,DYGF(NS),DYGF(N),HDL,INDX)
            IF( PERM_SX/EPSL.LT.EPSL ) KLM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KLM = 0.D+0
            INDX = 8
            RKLM = DIFMN(RKL(MN,NS),RKL(MP,N),DYGF(NS),DYGF(N),HDL,INDX)
            INDX = 5
            VLM = DIFMN(VISL(MN,NS),VISL(MP,N),DYGF(NS),DYGF(N),
     &        HDL,INDX)
            VL(M,1,N) = KLM*RKLM*HDLY/(DYGP(1,N)*VLM*RP(N))
          ENDDO
          DO M = 1,ISVF
            VL(M,2,NS) = VL(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            VL(M,1,N) = 0.D+0
          ENDDO
        ENDIF
!
!---    Bottom surface  ---
!
        NB = ICM(1,N)
        IF( NB.NE.0 ) THEN
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            HDLZ = PL(MN,NB)-PL(MP,N)-0.5D+0*GRVZ(1,N)
     &        *(RHOL(MN,NB)*DZGF(N)+RHOL(MP,N)*DZGF(NB))
            IF( M.EQ.1 ) HDL = HDLZ
!
!---        Permeability reduction factor
!
            PERM_BX = PERMRF(MN,NB)*PERM(3,NB)
            PERM_PX = PERMRF(MP,N)*PERM(3,N)
            INDX = 11
            KLM = DIFMN(PERM_BX,PERM_PX,DZGF(NB),DZGF(N),HDL,INDX)
            IF( PERM_BX/EPSL.LT.EPSL ) KLM = 0.D+0
            IF( PERM_PX/EPSL.LT.EPSL ) KLM = 0.D+0
            INDX = 8
            RKLM = DIFMN(RKL(MN,NB),RKL(MP,N),DZGF(NB),DZGF(N),
     &        HDL,INDX)
            INDX = 5
            VLM = DIFMN(VISL(MN,NB),VISL(MP,N),DZGF(NB),DZGF(N),
     &        HDL,INDX)
            WL(M,1,N) = KLM*RKLM*HDLZ/(DZGP(1,N)*VLM)
          ENDDO
          DO M = 1,ISVF
            WL(M,2,NB) = WL(M,1,N)
          ENDDO
        ELSE
          DO M = 1,ISVF
            WL(M,1,N) = 0.D+0
          ENDDO
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DRCVL group  ---
!
      RETURN
      END


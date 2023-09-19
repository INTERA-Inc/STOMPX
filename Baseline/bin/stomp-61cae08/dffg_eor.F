!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGC( IGC )
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
!     Diffusive petroleum component gas fluxes on interior surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE FDVG
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGC'
!
!---  X-direction vapor mole diffusion, excluding boundaries
!
      DO 200 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 200
        DO 150 NC = 1,4
          NW = ICM(NC,3,N)
          IF( NW.EQ.0 ) EXIT
!
!---      Multiple west connections  ---
!
          IF( INBS(3,N).LT.0 ) THEN
            NPX = NSSX(NW)
          ELSE
            NPX = NSX(N)
          ENDIF
          DXMGC = XMGC(IGC,2,NW) - XMGC(IGC,2,N)
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGC(IGC,MP,N)*RHOMG(MP,N)
            DFW = TORG(MN,NW)*PORD(MN,NW)*SG(MN,NW)*
     &        DFGC(IGC,MN,NW)*RHOMG(MN,NW)
            INDX = 12
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXMGC,INDX )
            UDGC(IGC,M,NPX) = DFM*(XMGC(IGC,MN,NW)-XMGC(IGC,MP,N))/
     &        DXGP(NPX)
            FGCP = XGC(IGC,MP,N)*RHOG(MP,N)
            FGCW = XGC(IGC,MN,NW)*RHOG(MN,NW)
            INDX = 3
            FGC = DIFMN( FGCW,FGCP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            UGC(IGC,M,NPX) = UG(M,NPX)*FGC + GCPP(1,IGC)*UDGC(IGC,M,NPX)
  100     CONTINUE
  150   CONTINUE
  200 CONTINUE
!
!---  Y-direction vapor mole diffusion, excluding boundaries
!
      DO 400 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 400
        DO 350 NC = 1,4
          NS = ICM(NC,2,N)
          IF( NS.EQ.0 ) EXIT
          I = ID(N)
!
!---      Multiple south connections  ---
!
          IF( INBS(2,N).LT.0 ) THEN
            NPY = NSSY(NS)
          ELSE
            NPY = NSY(N)
          ENDIF
          DXMGC = XMGC(IGC,2,NS) - XMGC(IGC,2,N)
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGC(IGC,MP,N)*RHOMG(MP,N)
            DFS = TORG(MN,NS)*PORD(MN,NS)*SG(MN,NS)*
     &        DFGC(IGC,MN,NS)*RHOMG(MN,NS)
            INDX = 12
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXMGC,INDX )
            VDGC(IGC,M,NPY) = DFM*(XMGC(IGC,MN,NS)-XMGC(IGC,MP,N))/
     &        (DYGP(NPY)*RP(I))
            FGCP = XGC(IGC,MP,N)*RHOG(MP,N)
            FGCS = XGC(IGC,MN,NS)*RHOG(MN,NS)
            INDX = 3
            FGC = DIFMN( FGCS,FGCP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
            VGC(IGC,M,NPY) = VG(M,NPY)*FGC + GCPP(1,IGC)*VDGC(IGC,M,NPY)
  300     CONTINUE
  350   CONTINUE
  400 CONTINUE
!
!---  Z-direction vapor mole diffusion, excluding boundaries
!
      DO 600 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 600
        DO 550 NC = 1,4
          NB = ICM(NC,1,N)
          IF( NB.EQ.0 ) EXIT
!
!---      Multiple bottom connections  ---
!
          IF( INBS(1,N).LT.0 ) THEN
            NPZ = NSSZ(NB)
          ELSE
            NPZ = NSZ(N)
          ENDIF
          DXMGC = XMGC(IGC,2,NB) - XMGC(IGC,2,N)
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*
     &        DFGC(IGC,MP,N)*RHOMG(MP,N)
            DFB = TORG(MN,NB)*PORD(MN,NB)*SG(MN,NB)*
     &        DFGC(IGC,MN,NB)*RHOMG(MN,NB)
            INDX = 12
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXMGC,INDX )
            WDGC(IGC,M,NPZ) = DFM*(XMGC(IGC,MN,NB)-XMGC(IGC,MP,N))/
     &        DZGP(NPZ)
            FGCP = XGC(IGC,MP,N)*RHOG(MP,N)
            FGCB = XGC(IGC,MN,NB)*RHOG(MN,NB)
            INDX = 3
            FGC = DIFMN( FGCB,FGCP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
            WGC(IGC,M,NPZ) = WG(M,NPZ)*FGC + GCPP(1,IGC)*WDGC(IGC,M,NPZ)
  500     CONTINUE
  550   CONTINUE
  600 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGC group  ---
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
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
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
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGW'
!
!---  X-direction vapor mole diffusion, excluding boundaries
!
      DO 200 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 200
        DO 150 NC = 1,4
          NW = ICM(NC,3,N)
          IF( NW.EQ.0 ) EXIT
!
!---      Multiple west connections  ---
!
          IF( INBS(3,N).LT.0 ) THEN
            NPX = NSSX(NW)
          ELSE
            NPX = NSX(N)
          ENDIF
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
  150   CONTINUE
  200 CONTINUE
!
!---  Y-direction vapor mole diffusion, excluding boundaries
!
      DO 400 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 400
        DO 350 NC = 1,4
          NS = ICM(NC,2,N)
          IF( NS.EQ.0 ) EXIT
          I = ID(N)
!
!---      Multiple south connections  ---
!
          IF( INBS(2,N).LT.0 ) THEN
            NPY = NSSY(NS)
          ELSE
            NPY = NSY(N)
          ENDIF
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
  350   CONTINUE
  400 CONTINUE
!
!---  Z-direction vapor mole diffusion, excluding boundaries
!
      DO 600 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 600
        DO 550 NC = 1,4
          NB = ICM(NC,1,N)
          IF( NB.EQ.0 ) EXIT
!
!---      Multiple bottom connections  ---
!
          IF( INBS(1,N).LT.0 ) THEN
            NPZ = NSSZ(NB)
          ELSE
            NPZ = NSZ(N)
          ENDIF
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
  550   CONTINUE
  600 CONTINUE
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
      SUBROUTINE DFFGCB( IGC,N,NB )
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
!     Diffusive petroleum component gas fluxes on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE FDVG
      USE CONST
      USE CCP
      USE BCVP
      USE BCVGC
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGCB'
      ITFX = MOD(IBCT(2,NB),10)
      K = KD(N)
      NPZ = NSZ(N)
      DXMGC = XMGCB(IGC,2,NB) - XMGC(IGC,2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FGCP = XGC(IGC,MP,N)*RHOG(MP,N)
        FGCB = XGCB(IGC,MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGC = DIFMN( FGCB,FGCP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( WG(1,NPZ).LT.-EPSL ) THEN
            WGC(IGC,M,NPZ) = WG(M,NPZ)*FGC
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( WG(1,NPZ).GT.EPSL ) THEN
            WGC(IGC,M,NPZ) = WG(M,NPZ)*FGC
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGC(IGC,MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGCB(IGC,MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMGC,INDX )
          WDGC(IGC,M,NPZ) = 2.D+0*DFM*(XMGCB(IGC,MP,NB)-XMGC(IGC,MP,N))/
     &      DZGF(N)
          WGC(IGC,M,NPZ) = WG(M,NPZ)*FGC + GCPP(1,IGC)*WDGC(IGC,M,NPZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGCB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGCE( IGC,N,NB )
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
!     Diffusive petroleum component gas fluxes on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE FDVG
      USE CONST
      USE CCP
      USE BCVP
      USE BCVGC
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGCE'
      ITFX = MOD(IBCT(2,NB),10)
      I = ID(N)
      NQX = NSX(N)+1
      DXMGC = XMGC(IGC,2,N) - XMGCB(IGC,2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FGCP = XGC(IGC,MN,N)*RHOG(MN,N)
        FGCB = XGCB(IGC,MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGC = DIFMN( FGCP,FGCB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( UG(1,NQX).GT.EPSL ) THEN
            UGC(IGC,M,NQX) = UG(M,NQX)*FGC
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( UG(1,NQX).LT.-EPSL ) THEN
            UGC(IGC,M,NQX) = UG(M,NQX)*FGC
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGC(IGC,MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGCB(IGC,MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMGC,INDX )
          UDGC(IGC,M,NQX) = 2.D+0*DFM*(XMGC(IGC,MN,N)-XMGCB(IGC,MN,NB))/
     &      DXGF(N)
          UGC(IGC,M,NQX) = UG(M,NQX)*FGC + GCPP(1,IGC)*UDGC(IGC,M,NQX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGCE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGCN( IGC,N,NB )
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
!     Diffusive petroleum component gas fluxes on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE FDVG
      USE CONST
      USE CCP
      USE BCVP
      USE BCVGC
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGCN'
      ITFX = MOD(IBCT(2,NB),10)
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMGC = XMGC(IGC,2,N) - XMGCB(IGC,2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FGCP = XGC(IGC,MN,N)*RHOG(MN,N)
        FGCB = XGCB(IGC,MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGC = DIFMN( FGCP,FGCB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( VG(1,NQY).GT.EPSL ) THEN
            VGC(IGC,M,NQY) = VG(M,NQY)*FGC
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( VG(1,NQY).LT.-EPSL ) THEN
            VGC(IGC,M,NQY) = VG(M,NQY)*FGC
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGC(IGC,MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGCB(IGC,MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMGC,INDX )
          VDGC(IGC,M,NQY) = 2.D+0*DFM*(XMGC(IGC,MN,N)-XMGCB(IGC,MN,NB))
     &      /(DYGF(N)*RP(I))
          VGC(IGC,M,NQY) = VG(M,NQY)*FGC + GCPP(1,IGC)*VDGC(IGC,M,NQY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGCN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGCS( IGC,N,NB )
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
!     Diffusive petroleum component gas fluxes on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE FDVG
      USE CONST
      USE CCP
      USE BCVP
      USE BCVGC
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGCS'
      ITFX = MOD(IBCT(2,NB),10)
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMGC = XMGCB(IGC,2,NB) - XMGC(IGC,2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FGCP = XGC(IGC,MP,N)*RHOG(MP,N)
        FGCB = XGCB(IGC,MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGC = DIFMN( FGCB,FGCP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( VG(1,NPY).LT.-EPSL ) THEN
            VGC(IGC,M,NPY) = VG(M,NPY)*FGC
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( VG(1,NPY).GT.EPSL ) THEN
            VGC(IGC,M,NPY) = VG(M,NPY)*FGC
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGC(IGC,MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGCB(IGC,MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMGC,INDX )
          VDGC(IGC,M,NPY) = 2.D+0*DFM*(XMGCB(IGC,MP,NB)-XMGC(IGC,MP,N))
     &      /(DYGF(N)*RP(I))
          VGC(IGC,M,NPY) = VG(M,NPY)*FGC + GCPP(1,IGC)*VDGC(IGC,M,NPY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGCS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGCT( IGC,N,NB )
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
!     Diffusive petroleum component gas fluxes on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE FDVG
      USE CONST
      USE CCP
      USE BCVP
      USE BCVGC
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGCT'
      ITFX = MOD(IBCT(2,NB),10)
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMGC = XMGC(IGC,2,N)-XMGCB(IGC,2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FGCP = XGC(IGC,MN,N)*RHOG(MN,N)
        FGCB = XGCB(IGC,MN,NB)*RHOGB(MN,NB)
        INDX = 3
        FGC = DIFMN( FGCP,FGCB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( WG(1,NQZ).GT.EPSL ) THEN
            WGC(IGC,M,NQZ) = WG(M,NQZ)*FGC
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( WG(1,NQZ).LT.-EPSL ) THEN
            WGC(IGC,M,NQZ) = WG(M,NQZ)*FGC
          ENDIF
        ELSE
          DFP = TORG(MN,N)*PORD(MN,N)*SG(MN,N)*DFGC(IGC,MN,N)
     &      *RHOMG(MN,N)
          DFB = TORGB(MN,NB)*PORDB(MN,NB)*SGB(MN,NB)*DFGCB(IGC,MN,NB)
     &      *RHOMGB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMGC,INDX )
          WDGC(IGC,M,NQZ) = 2.D+0*DFM*(XMGC(IGC,MN,N)-XMGCB(IGC,MN,NB))/
     &      DZGF(N)
          WGC(IGC,M,NQZ) = WG(M,NQZ)*FGC + GCPP(1,IGC)*WDGC(IGC,M,NQZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGCT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFGCW( IGC,N,NB )
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
!     Diffusive petroleum component gas fluxes on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXP
      USE FLUXGC
      USE FDVP
      USE FDVGC
      USE FDVG
      USE CONST
      USE CCP
      USE BCVP
      USE BCVGC
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGCW'
      ITFX = MOD(IBCT(2,NB),10)
      I = ID(N)
      NPX = NSX(N)
      DXMGC = XMGCB(IGC,2,NB)-XMGC(IGC,2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FGCP = XGC(IGC,MP,N)*RHOG(MP,N)
        FGCB = XGCB(IGC,MP,NB)*RHOGB(MP,NB)
        INDX = 3
        FGC = DIFMN( FGCB,FGCP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( UG(1,NPX).LT.-EPSL ) THEN
            UGC(IGC,M,NPX) = UG(M,NPX)*FGC
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( UG(1,NPX).GT.EPSL ) THEN
            UGC(IGC,M,NPX) = UG(M,NPX)*FGC
          ENDIF
        ELSE
          DFP = TORG(MP,N)*PORD(MP,N)*SG(MP,N)*DFGC(IGC,MP,N)
     &      *RHOMG(MP,N)
          DFB = TORGB(MP,NB)*PORDB(MP,NB)*SGB(MP,NB)*DFGCB(IGC,MP,NB)
     &      *RHOMGB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMGC,INDX )
          UDGC(IGC,M,NPX) = 2.D+0*DFM*(XMGCB(IGC,MP,NB)-XMGC(IGC,MP,N))/
     &      DXGF(N)
          UGC(IGC,M,NPX) = UG(M,NPX)*FGC + GCPP(1,IGC)*UDGC(IGC,M,NPX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFGCW group  ---
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
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGWB'
      ITFX = MOD(IBCT(2,NB),10)
      K = KD(N)
      NPZ = NSZ(N)
      DXMGW = XMGWB(2,NB) - XMGW(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
        FGWP = XGW(MP,N)*RHOG(MP,N)
        INDX = 3
        FGW = DIFMN( FGWB,FGWP,DZGF(N),DZGF(N),WG(1,NPZ),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( WG(1,NPZ).LT.-EPSL ) THEN
            WGW(M,NPZ) = WG(M,NPZ)*FGW
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
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
!     Diffusive water gas fluxes on east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGWE'
      ITFX = MOD(IBCT(2,NB),10)
      I = ID(N)
      NQX = NSX(N)+1
      DXMGW = XMGW(2,N) - XMGWB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
        FGWP = XGW(MN,N)*RHOG(MN,N)
        INDX = 3
        FGW = DIFMN( FGWP,FGWB,DXGF(N),DXGF(N),UG(1,NQX),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( UG(1,NQX).GT.EPSL ) THEN
            UGW(M,NQX) = UG(M,NQX)*FGW
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
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
!     Diffusive water gas fluxes on north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGWN'
      ITFX = MOD(IBCT(2,NB),10)
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMGW = XMGW(2,N) - XMGWB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
        FGWP = XGW(MN,N)*RHOG(MN,N)
        INDX = 3
        FGW = DIFMN( FGWP,FGWB,DYGF(N),DYGF(N),VG(1,NQY),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( VG(1,NQY).GT.EPSL ) THEN
            VGW(M,NQY) = VG(M,NQY)*FGW
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
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
!     Diffusive water gas fluxes on south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGWS'
      ITFX = MOD(IBCT(2,NB),10)
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMGW = XMGWB(2,NB) - XMGW(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
        FGWP = XGW(MP,N)*RHOG(MP,N)
        INDX = 3
        FGW = DIFMN( FGWB,FGWP,DYGF(N),DYGF(N),VG(1,NPY),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( VG(1,NPY).LT.-EPSL ) THEN
            VGW(M,NPY) = VG(M,NPY)*FGW
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
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
!     Diffusive water gas fluxes on top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGWT'
      ITFX = MOD(IBCT(2,NB),10)
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMGW = XMGW(2,N)-XMGWB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FGWB = XGWB(MN,NB)*RHOGB(MN,NB)
        FGWP = XGW(MN,N)*RHOG(MN,N)
        INDX = 3
        FGW = DIFMN( FGWP,FGWB,DZGF(N),DZGF(N),WG(1,NQZ),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( WG(1,NQZ).GT.EPSL ) THEN
            WGW(M,NQZ) = WG(M,NQZ)*FGW
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
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
!     Diffusive water gas fluxes on west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFGWW'
      ITFX = MOD(IBCT(2,NB),10)
      I = ID(N)
      NPX = NSX(N)
      DXMGW = XMGWB(2,NB)-XMGW(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FGWB = XGWB(MP,NB)*RHOGB(MP,NB)
        FGWP = XGW(MP,N)*RHOG(MP,N)
        INDX = 3
        FGW = DIFMN( FGWB,FGWP,DXGF(N),DXGF(N),UG(1,NPX),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( UG(1,NPX).LT.-EPSL ) THEN
            UGW(M,NPX) = UG(M,NPX)*FGW
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
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


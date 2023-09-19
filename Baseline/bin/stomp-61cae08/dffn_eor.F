!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNC( IGC )
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
!     Diffusive petroleum component nonaqueous-liquid fluxes on 
!     interior surfaces.
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
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVN
      USE FDVGC
      USE CCP
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
      SUB_LOG(ISUB_LOG) = '/DFFNC'
!
!---  X-direction nonaqueous-liquid molar diffusion, excluding boundaries
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
          DXMNC = XMNC(IGC,2,NW) - XMNC(IGC,2,N)
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &        DFNC(IGC,MP,N)*RHOMN(MP,N)
            DFW = TORN(MN,NW)*PORD(MN,NW)*SN(MN,NW)*
     &        DFNC(IGC,MN,NW)*RHOMN(MN,NW)
            INDX = 20
            DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXMNC,INDX )
            UDNC(IGC,M,NPX) = DFM*(XMNC(IGC,MN,NW)-XMNC(IGC,MP,N))/
     &        DXGP(NPX)
            FNCP = XNC(IGC,MP,N)*RHON(MP,N)
            FNCW = XNC(IGC,MN,NW)*RHON(MN,NW)
            INDX = 4
            FNC = DIFMN( FNCW,FNCP,DXGF(NW),DXGF(N),UN(1,NPX),INDX )
            UNC(IGC,M,NPX) = UN(M,NPX)*FNC + GCPP(1,IGC)*UDNC(IGC,M,NPX)
  100     CONTINUE
  150   CONTINUE
  200 CONTINUE
!
!---  Y-direction nonaqueous-liquid molar diffusion, excluding boundaries
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
          DXMNC = XMNC(IGC,2,NS) - XMNC(IGC,2,N)
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &        DFNC(IGC,MP,N)*RHOMN(MP,N)
            DFS = TORN(MN,NS)*PORD(MN,NS)*SN(MN,NS)*
     &        DFNC(IGC,MN,NS)*RHOMN(MN,NS)
            INDX = 20
            DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXMNC,INDX )
            VDNC(IGC,M,NPY) = DFM*(XMNC(IGC,MN,NS)-XMNC(IGC,MP,N))/
     &        (DYGP(NPY)*RP(I))
            FNCP = XNC(IGC,MP,N)*RHON(MP,N)
            FNCS = XNC(IGC,MN,NS)*RHON(MN,NS)
            INDX = 4
            FNC = DIFMN( FNCS,FNCP,DYGF(NS),DYGF(N),VN(1,NPY),INDX )
            VNC(IGC,M,NPY) = VN(M,NPY)*FNC + GCPP(1,IGC)*VDNC(IGC,M,NPY)
  300     CONTINUE
  350   CONTINUE
  400 CONTINUE
!
!---  Z-direction nonaqueous-liquid molar diffusion, excluding boundaries
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
          DXMNC = XMNC(IGC,2,NB) - XMNC(IGC,2,N)
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*
     &        DFNC(IGC,MP,N)*RHOMN(MP,N)
            DFB = TORN(MN,NB)*PORD(MN,NB)*SN(MN,NB)*
     &        DFNC(IGC,MN,NB)*RHOMN(MN,NB)
            INDX = 20
            DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXMNC,INDX )
            WDNC(IGC,M,NPZ) = DFM*(XMNC(IGC,MN,NB)-XMNC(IGC,MP,N))/
     &        DZGP(NPZ)
            FNCP = XNC(IGC,MP,N)*RHON(MP,N)
            FNCB = XNC(IGC,MN,NB)*RHON(MN,NB)
            INDX = 4
            FNC = DIFMN( FNCB,FNCP,DZGF(NB),DZGF(N),WN(1,NPZ),INDX )
            WNC(IGC,M,NPZ) = WN(M,NPZ)*FNC + GCPP(1,IGC)*WDNC(IGC,M,NPZ)
  500     CONTINUE
  550   CONTINUE
  600 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNCB( IGC,N,NB )
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
!     Diffusive petroleum component nonaqueous-liquid fluxes on a 
!     bottom boundary.
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
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVN
      USE FDVGC
      USE CONST
      USE CCP
      USE BCVP
      USE BCVN
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFNCB'
      ITFX = MOD(IBCT(2,NB),10)
      K = KD(N)
      NPZ = NSZ(N)
      DXMNC = XMNCB(IGC,2,NB) - XMNC(IGC,2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FNCP = XNC(IGC,MP,N)*RHON(MP,N)
        FNCB = XNCB(IGC,MP,NB)*RHONB(MP,NB)
        INDX = 3
        FNC = DIFMN( FNCB,FNCP,DZGF(N),DZGF(N),WN(1,NPZ),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( WN(1,NPZ).LT.-EPSL ) THEN
            WNC(IGC,M,NPZ) = WN(M,NPZ)*FNC
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( WN(1,NPZ).GT.EPSL ) THEN
            WNC(IGC,M,NPZ) = WN(M,NPZ)*FNC
          ENDIF
        ELSE
          DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*DFNC(IGC,MP,N)
     &      *RHOMN(MP,N)
          DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNCB(IGC,MP,NB)
     &      *RHOMNB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMNC,INDX )
          WDNC(IGC,M,NPZ) = 2.D+0*DFM*(XMNCB(IGC,MP,NB)-XMNC(IGC,MP,N))/
     &      DZGF(N)
          WNC(IGC,M,NPZ) = WN(M,NPZ)*FNC + GCPP(1,IGC)*WDNC(IGC,M,NPZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNCB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNCE( IGC,N,NB )
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
!     Diffusive petroleum component nonaqueous-liquid fluxes on an 
!     east boundary.
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
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVN
      USE FDVGC
      USE CONST
      USE CCP
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNCE'
      ITFX = MOD(IBCT(2,NB),10)
      I = ID(N)
      NQX = NSX(N)+1
      DXMNC = XMNC(IGC,2,N) - XMNCB(IGC,2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FNCP = XNC(IGC,MN,N)*RHON(MN,N)
        FNCB = XNCB(IGC,MN,NB)*RHONB(MN,NB)
        INDX = 3
        FNC = DIFMN( FNCP,FNCB,DXGF(N),DXGF(N),UN(1,NQX),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( UN(1,NQX).GT.EPSL ) THEN
            UNC(IGC,M,NQX) = UN(M,NQX)*FNC
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( UN(1,NQX).LT.-EPSL ) THEN
            UNC(IGC,M,NQX) = UN(M,NQX)*FNC
          ENDIF
        ELSE
          DFP = TORN(MN,N)*PORD(MN,N)*SN(MN,N)*DFNC(IGC,MN,N)
     &      *RHOMN(MN,N)
          DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNCB(IGC,MN,NB)
     &      *RHOMNB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMNC,INDX )
          UDNC(IGC,M,NQX) = 2.D+0*DFM*(XMNC(IGC,MN,N)-XMNCB(IGC,MN,NB))/
     &      DXGF(N)
          UNC(IGC,M,NQX) = UN(M,NQX)*FNC + GCPP(1,IGC)*UDNC(IGC,M,NQX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNCE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNCN( IGC,N,NB )
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
!     Diffusive petroleum component nonaqueous-liquid fluxes on a 
!     north boundary.
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
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVN
      USE FDVGC
      USE CONST
      USE CCP
      USE BCVP
      USE BCVN
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFNCN'
      ITFX = MOD(IBCT(2,NB),10)
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DXMNC = XMNC(IGC,2,N) - XMNCB(IGC,2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FNCP = XNC(IGC,MN,N)*RHON(MN,N)
        FNCB = XNCB(IGC,MN,NB)*RHONB(MN,NB)
        INDX = 3
        FNC = DIFMN( FNCP,FNCB,DYGF(N),DYGF(N),VN(1,NQY),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( VN(1,NQY).GT.EPSL ) THEN
            VNC(IGC,M,NQY) = VN(M,NQY)*FNC
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( VN(1,NQY).LT.-EPSL ) THEN
            VNC(IGC,M,NQY) = VN(M,NQY)*FNC
          ENDIF
        ELSE
          DFP = TORN(MN,N)*PORD(MN,N)*SN(MN,N)*DFNC(IGC,MN,N)
     &      *RHOMN(MN,N)
          DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNCB(IGC,MN,NB)
     &      *RHOMNB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMNC,INDX )
          VDNC(IGC,M,NQY) = 2.D+0*DFM*(XMNC(IGC,MN,N)-XMNCB(IGC,MN,NB))
     &      /(DYGF(N)*RP(I))
          VNC(IGC,M,NQY) = VN(M,NQY)*FNC + GCPP(1,IGC)*VDNC(IGC,M,NQY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNCN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNCS( IGC,N,NB )
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
!     Diffusive petroleum component nonaqueous-liquid fluxes on a 
!     south boundary.
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
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVN
      USE FDVGC
      USE CONST
      USE CCP
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNCS'
      ITFX = MOD(IBCT(2,NB),10)
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DXMNC = XMNCB(IGC,2,NB) - XMNC(IGC,2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FNCP = XNC(IGC,MP,N)*RHON(MP,N)
        FNCB = XNCB(IGC,MP,NB)*RHONB(MP,NB)
        INDX = 3
        FNC = DIFMN( FNCB,FNCP,DYGF(N),DYGF(N),VN(1,NPY),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( VN(1,NPY).LT.-EPSL ) THEN
            VNC(IGC,M,NPY) = VN(M,NPY)*FNC
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( VN(1,NPY).GT.EPSL ) THEN
            VNC(IGC,M,NPY) = VN(M,NPY)*FNC
          ENDIF
        ELSE
          DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*DFNC(IGC,MP,N)
     &      *RHOMN(MP,N)
          DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNCB(IGC,MP,NB)
     &      *RHOMNB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMNC,INDX )
          VDNC(IGC,M,NPY) = 2.D+0*DFM*(XMNCB(IGC,MP,NB)-XMNC(IGC,MP,N))
     &      /(DYGF(N)*RP(I))
          VNC(IGC,M,NPY) = VN(M,NPY)*FNC + GCPP(1,IGC)*VDNC(IGC,M,NPY)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNCS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNCT( IGC,N,NB )
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
!     Diffusive petroleum component nonaqueous-liquid fluxes on a 
!     top boundary.
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
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVN
      USE FDVGC
      USE CONST
      USE CCP
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNCT'
      ITFX = MOD(IBCT(2,NB),10)
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DXMNC = XMNC(IGC,2,N)-XMNCB(IGC,2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        FNCP = XNC(IGC,MN,N)*RHON(MN,N)
        FNCB = XNCB(IGC,MN,NB)*RHONB(MN,NB)
        INDX = 3
        FNC = DIFMN( FNCP,FNCB,DZGF(N),DZGF(N),WN(1,NQZ),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( WN(1,NQZ).GT.EPSL ) THEN
            WNC(IGC,M,NQZ) = WN(M,NQZ)*FNC
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( WN(1,NQZ).LT.-EPSL ) THEN
            WNC(IGC,M,NQZ) = WN(M,NQZ)*FNC
          ENDIF
        ELSE
          DFP = TORN(MN,N)*PORD(MN,N)*SN(MN,N)*DFNC(IGC,MN,N)
     &      *RHOMN(MN,N)
          DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNCB(IGC,MN,NB)
     &      *RHOMNB(MN,NB)
          INDX = 12
          DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMNC,INDX )
          WDNC(IGC,M,NQZ) = 2.D+0*DFM*(XMNC(IGC,MN,N)-XMNCB(IGC,MN,NB))/
     &      DZGF(N)
          WNC(IGC,M,NQZ) = WN(M,NQZ)*FNC + GCPP(1,IGC)*WDNC(IGC,M,NQZ)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNCT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNCW( IGC,N,NB )
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
!     Diffusive petroleum component nonaqueous-liquid fluxes on a 
!     west boundary.
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
      USE FLUXN
      USE FLUXGC
      USE FDVP
      USE FDVN
      USE FDVGC
      USE CONST
      USE CCP
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/DFFNCW'
      ITFX = MOD(IBCT(2,NB),10)
      I = ID(N)
      NPX = NSX(N)
      DXMNC = XMNCB(IGC,2,NB)-XMNC(IGC,2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        FNCP = XNC(IGC,MP,N)*RHON(MP,N)
        FNCB = XNCB(IGC,MP,NB)*RHONB(MP,NB)
        INDX = 3
        FNC = DIFMN( FNCB,FNCP,DXGF(N),DXGF(N),UN(1,NPX),INDX )
!
!---    Dirichlet-Outflow or Neumann-Outflow boundary condition  ---
!
        IF( ITFX.EQ.6 .OR. ITFX.EQ.9 ) THEN
          IF( UN(1,NPX).LT.-EPSL ) THEN
            UNC(IGC,M,NPX) = UN(M,NPX)*FNC
          ENDIF
!
!---    Dirichlet-Inflow or Neumann-Inflow boundary condition  ---
!
        ELSEIF( ITFX.EQ.5 .OR. ITFX.EQ.7 ) THEN
          IF( UN(1,NPX).GT.EPSL ) THEN
            UNC(IGC,M,NPX) = UN(M,NPX)*FNC
          ENDIF
        ELSE
          DFP = TORN(MP,N)*PORD(MP,N)*SN(MP,N)*DFNC(IGC,MP,N)
     &      *RHOMN(MP,N)
          DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNCB(IGC,MP,NB)
     &      *RHOMNB(MP,NB)
          INDX = 12
          DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMNC,INDX )
          UDNC(IGC,M,NPX) = 2.D+0*DFM*(XMNCB(IGC,MP,NB)-XMNC(IGC,MP,N))/
     &      DXGF(N)
          UNC(IGC,M,NPX) = UN(M,NPX)*FNC + GCPP(1,IGC)*UDNC(IGC,M,NPX)
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNCW group  ---
!
      RETURN
      END

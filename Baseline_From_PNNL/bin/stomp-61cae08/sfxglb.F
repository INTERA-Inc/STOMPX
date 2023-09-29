!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFXGLB( NSL )
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
!     Compute solute gas-phase fluxes on boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1995.
!     Last Modified by MD White, Battelle, PNL, October 13, 1995.
!     sfxglb.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
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
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SFXGLB'
!
!---  Loop over number of specified boundary conditions  ---
!
      DO 200 NB = 1,NBC
        TMZ = TM
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) GOTO 200
        IF( IBCM(NB).GT.1 .AND. TMZ.GT.BC(1,IBCM(NB),MB) ) GOTO 200
        IF( IBCM(NB) .EQ. 1 ) THEN
          BCX = BC(NSL+LBCU,1,MB)
        ELSE
          DO 100 M = 2,IBCM(NB)
            IF( TMZ.LE.BC(1,M,MB) ) THEN
             DTBC = MIN( BC(1,M,MB)-TM,DT )
             TFBC = (TM-5.D-1*DTBC-BC(1,M-1,MB))/
     &         (BC(1,M,MB)-BC(1,M-1,MB))
             BCX = BC(NSL+LBCU,M-1,MB) +
     &         TFBC*(BC(NSL+LBCU,M,MB)-BC(NSL+LBCU,M-1,MB))
             GOTO 110
            ENDIF
  100     CONTINUE
          GOTO 200
        ENDIF
  110   CONTINUE
        N = IBCN(NB)
        IZN = IZ(N)
        MP = IXP(N)
        I = ID(N)
        J = JD(N)
        K = KD(N)
!
!---    Diffusion coefficients at node adjacent to boundary  ---
!
        TCOR = (T(2,N)+TABS)/TSPRF
        PCOR = (PG(2,N)+PATM)/PATM
        SDFGP = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
        DGP = TORG(2,N)*(SG(2,N)-SGT(2,N))*PORD(2,N)*SDFGP
!
!---    Phase fraction factors at node adjacent to boundary  ---
!
        XVLP = SG(2,N)*PORD(2,N)
        FCGP = 0.D+0
        IF( XVLP.GT.EPSL ) FCGP = YG(N,NSL)/XVLP
!
!---    Compute boundary phase fractions  ---
!
        IF( IPCL(NSL).EQ.2 ) THEN
          XVSB = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORTB(2,NB))*SLB(2,NB)
        ELSE
          XVSB = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORTB(2,NB))
        ENDIF
        XVLB = SLB(2,NB)*PORDB(2,NB)
        XVGB = SGB(2,NB)*PORDB(2,NB)
!
!---    Constant gas-aqueous partition coefficient  ---
!
        IF( IPCGL(NSL).EQ.0 ) THEN
          PCGLX = PCGL(1,NSL)
!
!---    Temperature dependent gas-aqueous partition coefficient  ---
!
        ELSEIF( IPCGL(NSL).EQ.1 ) THEN
          TK = TB(2,NB)+TABS
          PCGLX = EXP( PCGL(1,NSL) + PCGL(2,NSL)/TK
     &      + PCGL(3,NSL)*LOG(TK) + PCGL(4,NSL)*TK + PCGL(5,NSL)*TK**2 )
!
!---    Water-vapor equilibrium gas-aqueous partition coefficient  ---
!
        ELSEIF( IPCGL(NSL).EQ.2 ) THEN
          PCGLX = RHOG(2,N)*XGW(2,N)/(RHOL(2,N)*XLW(2,N))
        ENDIF
        PCGLX = MAX( PCGLX,1.D-20 )
        PCGLX = MIN( PCGLX,1.D+20 )
!
!---    Phase-volumetric concentration ratios  ---
!
        FCLB = 1.D+0/(XVSB + XVLB + XVGB*PCGLX)
        FCGB = 1.D+0/((XVSB + XVLB)/PCGLX + XVGB)
!
!---    Phase mole fractions  ---
!
        YLB(NB,NSL) = XVLB*FCLB
        YGB(NB,NSL) = XVGB*FCGB
!
!---    Convert boundary concentrations  ---
!
        IF( IBCT(IEQC,NB).EQ.8 ) THEN
          BCX = BCX/(FCLB+SMALL)
        ELSEIF( IBCT(IEQC,NB).EQ.9 ) THEN
          BCX = BCX/(FCGB+SMALL)
        ENDIF
!
!---    Bottom boundary  ---
!
        IF( IBCD(NB) .EQ. -3 ) THEN
          NPZ = NSZ(N)
!
!---      Hydraulic dispersion  ---
!
          IF( IDISP.EQ.1 ) THEN
            UGBX = (0.5D+0*(UG(1,NSX(N))+UG(1,NSX(N)+1)))**2
            VGBX = (0.5D+0*(VG(1,NSY(N))+VG(1,NSY(N)+IFLD)))**2
            WGBX = (WG(1,NPZ))**2
            ZGBX = SQRT(UGBX + VGBX + WGBX)
            DPGB = (DISPL(IZN)*WGBX + DISPT(IZN)*(UGBX+VGBX))/
     &        (ZGBX+SMALL)
          ELSE
            DPGB = 0.D+0
          ENDIF
!
!---  Dirichlet  ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &      IBCT(IEQC,NB).EQ.9 .OR. IBCT(IEQC,NB).EQ.10 .OR.
     &      IBCT(IEQC,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            DGB = TORGB(2,NB)*XVGB*SMDGB
            INDX = 16
            DGZ = DIFMN(DGB,DGP,DZGF(N),DZGF(N),WG(1,NPZ),INDX)
            DGZ = (DGZ+DPGB)/(5.D-1*DZGF(N))
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 )  THEN
              AB = DGZ
              AP = DGZ
            ELSE
              AB = MAX( WG(1,NPZ),ZERO ) +
     &         DGZ*MAX((ONE-(TENTH*ABS(WG(1,NPZ))/(DGZ+SMALL)))**5,ZERO)
              AP = MAX( -WG(1,NPZ),ZERO ) +
     &         DGZ*MAX((ONE-(TENTH*ABS(WG(1,NPZ))/(DGZ+SMALL)))**5,ZERO)
            ENDIF
            WC(NPZ,NSL) = WC(NPZ,NSL) + (BCX*AB*FCGB - C(N,NSL)*AP*FCGP)
!
!---  Outflow ---
!
          ELSEIF( IBCT(IEQC,NB) .EQ. 7 ) THEN
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 )  THEN
              AP = 0.D+0
            ELSE
              AP = MAX( -WG(1,NPZ),ZERO )
            ENDIF
            WC(NPZ,NSL) = WC(NPZ,NSL) - C(N,NSL)*AP*FCGP
!
!---  Inflow ---
!
          ELSEIF( IBCT(IEQC,NB).GE.13 .AND. IBCT(IEQC,NB).LE.16 ) THEN
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 )  THEN
              AB = 0.D+0
            ELSE
              AB = MAX( WG(1,NPZ),ZERO )
            ENDIF
            WC(NPZ,NSL) = WC(NPZ,NSL) + BCX*AB*FCGB
          ENDIF
!
!---  South boundary  ---
!
        ELSEIF( IBCD(NB) .EQ. -2 ) THEN
          NPY = NSY(N)
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            UGSX = (0.5D+0*(UG(1,NSX(N))+UG(1,NSX(N)+1)))**2
            VGSX = VG(1,NPY)**2
            WGSX = (0.5D+0*(WG(1,NSZ(N))+WG(1,NSZ(N)+IJFLD)))**2
            ZGS = SQRT(UGSX + VGSX + WGSX)
            DPGS = (DISPL(IZ(N))*VGSX + DISPT(IZ(N))*(UGSX+WGSX))/
     &        (ZGS+SMALL)
          ELSE
            DPGS = 0.D+0
          ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &      IBCT(IEQC,NB).EQ.9 .OR. IBCT(IEQC,NB).EQ.10 .OR.
     &      IBCT(IEQC,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            DGB = TORGB(2,NB)*XVGB*SMDGB
            INDX = 16
            DGY = DIFMN(DGB,DGP,DYGF(N),DYGF(N),VG(1,NPY),INDX)
            DGY = (DGY+DPGS)/RP(I)/(5.D-1*DYGF(N))
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 )  THEN
              AS = DGY
              AP = DGY
            ELSE
              AS = MAX( VG(1,NPY),ZERO ) +
     &         DGY*MAX((ONE-(TENTH*ABS(VG(1,NPY))/(DGY+SMALL)))**5,ZERO)
              AP = MAX( -VG(1,NPY),ZERO ) +
     &         DGY*MAX((ONE-(TENTH*ABS(VG(1,NPY))/(DGY+SMALL)))**5,ZERO)
            ENDIF
            VC(NPY,NSL) = VC(NPY,NSL) + (BCX*AS*FCGB - C(N,NSL)*AP*FCGP)
!
!---  Outflow ---
!
          ELSEIF( IBCT(IEQC,NB) .EQ. 7 ) THEN
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 )  THEN
              AP = 0.D+0
            ELSE
              AP = MAX( -VG(1,NPY),ZERO )
            ENDIF
            VC(NPY,NSL) = VC(NPY,NSL) - C(N,NSL)*AP*FCGP
!
!---  Inflow ---
!
          ELSEIF( IBCT(IEQC,NB).GE.13 .AND. IBCT(IEQC,NB).LE.16 ) THEN
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 )  THEN
              AS = 0.D+0
            ELSE
              AS = MAX( VG(1,NPY),ZERO )
            ENDIF
            VC(NPY,NSL) = VC(NPY,NSL) - BCX*AS*FCGB
          ENDIF
!
!---  West boundary  ---
!
        ELSEIF( IBCD(NB) .EQ. -1 ) THEN
          NPX = NSX(N)
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            UGX = UG(1,NPX)**2
            VGX = (0.5D+0*(VG(1,NSY(N))+VG(1,NSY(N)+IFLD)))**2
            WGX = (0.5D+0*(WG(1,NSZ(N))+WG(1,NSZ(N)+IJFLD)))**2
            ZGW = SQRT(UGX + VGX + WGX)
            DPGW = (DISPL(IZ(N))*UGX + DISPT(IZ(N))*(WGX+VGX))/
     &          (ZGW+SMALL)
            ELSE
              DPGW = 0.D+0
            ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &      IBCT(IEQC,NB).EQ.9 .OR. IBCT(IEQC,NB).EQ.10 .OR.
     &      IBCT(IEQC,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            DGB = TORGB(2,NB)*XVGB*SMDGB
            INDX = 16
            DGX = DIFMN(DGB,DGP,DXGF(N),DXGF(N),UG(1,NPX),INDX)
            DGX = (DGX+DPGW)/(5.D-1*DXGF(N))
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 )  THEN
              AW = DGX
              AP = DGX
            ELSE
              AW = MAX( UG(1,NPX),ZERO ) +
     &         DGX*MAX((ONE-(TENTH*ABS(UG(1,NPX))/(DGX+SMALL)))**5,ZERO)
              AP = MAX( -UG(1,NPX),ZERO ) +
     &         DGX*MAX((ONE-(TENTH*ABS(UG(1,NPX))/(DGX+SMALL)))**5,ZERO)
            ENDIF
            UC(NPX,NSL) = UC(NPX,NSL) + (BCX*AW*FCGB - C(N,NSL)*AP*FCGP)
!
!---  Outflow ---
!
          ELSEIF( IBCT(IEQC,NB) .EQ. 7 ) THEN
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 )  THEN
              AP = 0.D+0
            ELSE
              AP = MAX( -UG(1,NPX),ZERO )
            ENDIF
            UC(NPX,NSL) = UC(NPX,NSL) -  C(N,NSL)*AP*FCGP
!
!---  Inflow ---
!
          ELSEIF( IBCT(IEQC,NB).GE.13 .AND. IBCT(IEQC,NB).LE.16 ) THEN
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 )  THEN
              AW = 0.D+0
            ELSE
              AW = MAX( UG(1,NPX),ZERO )
            ENDIF
            UC(NPX,NSL) = UC(NPX,NSL) + BCX*AW*FCGB
          ENDIF
!
!---  East boundary  ---
!
        ELSEIF( IBCD(NB) .EQ. 1 ) THEN
          NQX = NSX(N)+1
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            UGEX = UG(1,NQX)**2
            VGEX = (0.5D+0*(VG(1,NSY(N))+VG(1,NSY(N)+IFLD)))**2
            WGEX = (0.5D+0*(WG(1,NSZ(N))+WG(1,NSZ(N)+IJFLD)))**2
            ZGE = SQRT(UGEX + VGEX + WGEX)
            DPGE = (DISPL(IZ(N))*UGEX + DISPT(IZ(N))*(WGEX+VGEX))/
     &        (ZGE+SMALL)
          ELSE
            DPGE = 0.D+0
          ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &      IBCT(IEQC,NB).EQ.9 .OR. IBCT(IEQC,NB).EQ.10 .OR.
     &      IBCT(IEQC,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            DGB = TORGB(2,NB)*XVGB*SMDGB
            INDX = 16
            DGX = DIFMN(DGP,DGB,DXGF(N),DXGF(N),UG(1,NQX),INDX)
            DGX = (DGX+DPGE)/(5.D-1*DXGF(N))
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 )  THEN
              AE = DGX
              AP = DGX
            ELSE
              AE = MAX( -UG(1,NQX),ZERO ) +
     &         DGX*MAX((ONE-(TENTH*ABS(UG(1,NQX))/(DGX+SMALL)))**5,ZERO)
              AP = MAX( UG(1,NQX),ZERO ) +
     &         DGX*MAX((ONE-(TENTH*ABS(UG(1,NQX))/(DGX+SMALL)))**5,ZERO)
            ENDIF
            UC(NQX,NSL) = UC(NQX,NSL)+(C(N,NSL)*AP*FCGP-BCX*AE*FCGB)
!
!---  Outflow ---
!
          ELSEIF( IBCT(IEQC,NB) .EQ. 7 ) THEN
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 )  THEN
              AP = 0.D+0
            ELSE
              AP = MAX( UG(1,NQX),ZERO )
            ENDIF
            UC(NQX,NSL) = UC(NQX,NSL) + C(N,NSL)*AP*FCGP
!
!---  Inflow ---
!
          ELSEIF( IBCT(IEQC,NB).GE.13 .AND. IBCT(IEQC,NB).LE.16 ) THEN
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 )  THEN
              AE = 0.D+0
            ELSE
              AE = MAX( -UG(1,NQX),ZERO )
            ENDIF
            UC(NQX,NSL) = UC(NQX,NSL) - BCX*AE*FCGB
          ENDIF
!
!---  North boundary  ---
!
        ELSEIF( IBCD(NB) .EQ. 2 ) THEN
          NQY = NSY(N)+IFLD
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            UGNX = (0.5D+0*(UG(1,NSX(N))+UG(1,NSX(N)+1)))**2
            VGNX = VG(1,NQY)**2
            WGNX = (0.5D+0*(WG(1,NSZ(N))+WG(1,NSZ(N)+IJFLD)))**2
            ZGN = SQRT(UGNX + VGNX + WGNX)
            DPGN = (DISPL(IZ(N))*VGNX + DISPT(IZ(N))*(UGNX+WGNX))/
     &        (ZGN+SMALL)
          ELSE
            DPGN = 0.D+0
          ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &      IBCT(IEQC,NB).EQ.9 .OR. IBCT(IEQC,NB).EQ.10 .OR.
     &      IBCT(IEQC,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            DGB = TORGB(2,NB)*XVGB*SMDGB
            INDX = 16
            DGY = DIFMN(DGP,DGB,DYGF(N),DYGF(N),VG(1,NQY),INDX)
            DGY = (DGY+DPGN)/RP(I)/(5.D-1*DYGF(N))
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 )  THEN
              AN = DGY
              AP = DGY
            ELSE
              AN = MAX( -VG(1,NQY),ZERO ) +
     &         DGY*MAX((ONE-(TENTH*ABS(VG(1,NQY))/(DGY+SMALL)))**5,ZERO)
              AP = MAX( VG(1,NQY),ZERO ) +
     &         DGY*MAX((ONE-(TENTH*ABS(VG(1,NQY))/(DGY+SMALL)))**5,ZERO)
            ENDIF
            VC(NQY,NSL) = VC(NQY,NSL) + (C(N,NSL)*AP*FCGP - BCX*AN*FCGB)
!
!---  Outflow ---
!
          ELSEIF( IBCT(IEQC,NB) .EQ. 7 ) THEN
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 )  THEN
              AP = 0.D+0
            ELSE
              AP = MAX( VG(1,NQY),ZERO )
            ENDIF
            VC(NQY,NSL) = VC(NQY,NSL) + C(N,NSL)*AP*FCGP
!
!---  Inflow ---
!
          ELSEIF( IBCT(IEQC,NB).GE.13 .AND. IBCT(IEQC,NB).LE.16 ) THEN
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 )  THEN
              AN = 0.D+0
            ELSE
              AN = MAX( -VG(1,NQY),ZERO )
            ENDIF
            VC(NQY,NSL) = VC(NQY,NSL) - BCX*AN*FCGB
          ENDIF
!
!---  Top boundary
!
        ELSEIF( IBCD(NB) .EQ. 3 ) THEN
          NQZ = NSZ(N)+IJFLD
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            UGTX = (0.5D+0*(UG(1,NSX(N))+UG(1,NSX(N)+1)))**2
            VGTX = (0.5D+0*(VG(1,NSY(N))+VG(1,NSY(N)+IFLD)))**2
            WGTX = (WG(1,NQZ))**2
            ZGT = SQRT(UGTX + VGTX + WGTX)
            DPGT = (DISPL(IZ(N))*WGTX + DISPT(IZ(N))*(UGTX+VGTX))/
     &        (ZGT+SMALL)
          ELSE
            DPGT = 0.D+0
          ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &      IBCT(IEQC,NB).EQ.9 .OR. IBCT(IEQC,NB).EQ.10 .OR.
     &      IBCT(IEQC,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            DGB = TORGB(2,NB)*XVGB*SMDGB
            INDX = 16
            DGZ = DIFMN(DGP,DGB,DZGF(N),DZGF(N),WG(1,NQZ),INDX)
            DGZ = (DGZ+DPGT)/(5.D-1*DZGF(N))
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 )  THEN
              AT = DGZ
              AP = DGZ
            ELSE
              AT = MAX( -WG(1,NQZ),ZERO ) +
     &         DGZ*MAX((ONE-(TENTH*ABS(WG(1,NQZ))/(DGZ+SMALL)))**5,ZERO)
              AP = MAX( WG(1,NQZ),ZERO ) +
     &         DGZ*MAX((ONE-(TENTH*ABS(WG(1,NQZ))/(DGZ+SMALL)))**5,ZERO)
            ENDIF
            WC(NQZ,NSL) = WC(NQZ,NSL) + (C(N,NSL)*AP*FCGP - BCX*AT*FCGB)
!
!---  Outflow ---
!
          ELSEIF( IBCT(IEQC,NB) .EQ. 7 ) THEN
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 )  THEN
              AP = 0.D+0
            ELSE
              AP = MAX( WG(1,NQZ),ZERO )
            ENDIF
            WC(NQZ,NSL) = WC(NQZ,NSL) + C(N,NSL)*AP*FCGP
!
!---  Inflow ---
!
          ELSEIF( IBCT(IEQC,NB).GE.13 .AND. IBCT(IEQC,NB).LE.16 ) THEN
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 )  THEN
              AT = 0.D+0
            ELSE
              AT = MAX( -WG(1,NQZ),ZERO )
            ENDIF
            WC(NQZ,NSL) = WC(NQZ,NSL) - BCX*AT*FCGB
          ENDIF
        ENDIF
  200 CONTINUE
!
!---  End of SFXGLB group  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END


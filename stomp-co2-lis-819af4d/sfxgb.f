!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFXGB( NSL )
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
!     sfxgb.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 BCX(LSPBC+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SFXGB'
      SVGB = 0.D+0
      FCGB = 0.D+0
!
!---  Loop over number of specified boundary conditions  ---
!
      NBCT = MIN( NSL+LUK,NSOLU+LUK+1 )
      DO 200 NB = 1,NBC
        TMZ = TM
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) GOTO 200
        IF( IBCM(NB).GT.1 .AND. TMZ.GT.BC(1,IBCM(NB),MB) ) GOTO 200
        IF( IBCM(NB) .EQ. 1 ) THEN
!
!---      Solute transport  ---
!
          IF( NSL.LE.NSOLU ) THEN
            BCX(1) = BC(NSL+LBCU,1,MB)
!
!---      Reactive species transport  ---
!
          ELSE
            BCX(1) = 0.D+0
            DO 10 NSPX = 1,IBCSP(1,NB)
              MX = NSOLU+LBCU+NSPX
              BCX(NSPX+1) = BC(MX,1,MB)
   10       CONTINUE
          ENDIF
        ELSE
          DO 100 M = 2,IBCM(NB)
            IF( TMZ.LE.BC(1,M,MB) ) THEN
              DTBC = MIN( BC(1,M,MB)-TM,DT )
              TFBC = (TM-5.D-1*DTBC-BC(1,M-1,MB))/
     &          (BC(1,M,MB)-BC(1,M-1,MB))
!
!---          Solute transport  ---
!
              IF( NSL.LE.NSOLU ) THEN
                BCX(1) = BC(NSL+LBCU,M-1,MB) +
     &            TFBC*(BC(NSL+LBCU,M,MB)-BC(NSL+LBCU,M-1,MB))
                IF( IBCT(NBCT,NB).EQ.12 ) BCX(1) = CBO(NB,NSL)
!
!---          Reactive species transport  ---
!
              ELSE
                BCX(1) = 0.D+0
                DO 20 NSPX = 1,IBCSP(1,NB)
                  MX = NSOLU+LBCU+NSPX
                  BCX(NSPX+1) = BC(MX,M-1,MB) +
     &              TFBC*(BC(MX,M,MB)-BC(MX,M-1,MB))
                  IF( IBCT(NBCT,NB).EQ.12 ) BCX(NSPX) = CBO(NB,NSL)
   20           CONTINUE
              ENDIF
              GOTO 110
            ENDIF
  100     CONTINUE
          GOTO 200
        ENDIF
  110   CONTINUE
        N = IBCN(NB)
        MF = 1
        IZN = IZ(N)
        MP = IXP(N)
        I = ID(N)
        J = JD(N)
        K = KD(N)
!
!---  Compute adjacent node phase fractions  ---
!
        SVGP = SG(2,N)*PORD(2,N)
        FCGP = 0.D+0
        IF( SVGP.GT.SMALL ) FCGP = YG(N,NSL)/SVGP
!
!---    Solute transport only, skip calculations for reactive
!       species transport  ---
!
        IF( NSL.LE.NSOLU ) THEN
!
!---    Compute boundary phase fractions  ---
!
          IF( IPCL(NSL).EQ.2 ) THEN
            SVSB = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORTB(2,NB))*
     &        SLB(2,NB)
          ELSE
            SVSB = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORTB(2,NB))
          ENDIF
          SVLB = SLB(2,NB)*PORDB(2,NB)
          SVGB = SGB(2,NB)*PORDB(2,NB)
          SVNB = SNB(2,NB)*PORDB(2,NB)
!
!---      Constant gas-aqueous partition coefficient  ---
!
          IF( IPCGL(NSL).EQ.0 ) THEN
            PCGLX = PCGL(1,NSL)
!
!---      Temperature dependent gas-aqueous partition coefficient  ---
!
          ELSEIF( IPCGL(NSL).EQ.1 ) THEN
            TK = TB(2,NB)+TABS
            PCGLX = EXP( PCGL(1,NSL) + PCGL(2,NSL)/TK
     &        + PCGL(3,NSL)*LOG(TK) + PCGL(4,NSL)*TK
     &        + PCGL(5,NSL)*TK**2 )
!
!---      Water-vapor equilibrium gas-aqueous partition coefficient  ---
!
          ELSEIF( IPCGL(NSL).EQ.2 ) THEN
            PCGLX = RHOG(2,N)*XGW(2,N)/(RHOL(2,N)*XLW(2,N))
          ENDIF
          PCGLX = MAX( PCGLX,1.D-20 )
          PCGLX = MIN( PCGLX,1.D+20 )
!
!---      Phase-volumetric concentration ratios  ---
!
          FCLB = 1.D+0/(SVSB + SVLB + SVNB/PCLN(1,NSL)
     &       + SVGB*PCGLX)
          FCGB = 1.D+0/((SVSB + SVLB + SVNB)/PCGLX + SVGB)
          FCNB = 1.D+0/((SVSB + SVLB + SVGB*PCGLX)*PCLN(1,NSL) + SVNB)
! 
!---      Phase mole fractions  ---
!
          YLB(NB,NSL) = SVLB*FCLB
          YGB(NB,NSL) = SVGB*FCGB
          YNB(NB,NSL) = SVNB*FCNB
!
!---      Convert boundary concentrations  ---
!
          IF( IBCT(NBCT,NB).EQ.8 .OR. IBCT(NBCT,NB).EQ.14
     &      .OR. IBCT(NBCT,NB).EQ.23 ) THEN
            BCX(1) = BCX(1)/(FCLB+SMALL)
          ELSEIF( IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.15 ) THEN
            BCX(1) = BCX(1)/(FCGB+SMALL)
          ELSEIF( IBCT(NBCT,NB).EQ.10 .OR. IBCT(NBCT,NB).EQ.16 ) THEN
            BCX(1) = BCX(1)/(FCNB+SMALL)
          ENDIF
        ELSE
!
!---      Skip for initial condition type boundary condition  ---
!
          IF( IBCT(NBCT,NB).NE.12 ) THEN
!
!---        Convert species concentrations to total-component
!           concentrations  ---
!
            IF( NSL.LE.NSOLU+NEQC ) THEN
              NEQ = NSL-NSOLU
              DO 130 NSP = 1,IEQ_C(1,NEQ)
                DO 120 NSPX = 1,IBCSP(1,NB)
                  IF( IBCSP(NSPX+1,NB).EQ.IEQ_C(NSP+1,NEQ) ) THEN
                    BCX(1) = BCX(1) + EQ_C(NSP,NEQ)*BCX(NSPX+1)
                  ENDIF
  120           CONTINUE
  130         CONTINUE
!
!---        Convert species concentrations to total-kinetic
!           concentrations  ---
!
            ELSEIF( NSL.LE.NSOLU+NEQC+NEQK ) THEN
              NEQ = NSL-NSOLU-NEQC
              DO 150 NSP = 1,IEQ_K(1,NEQ)
                DO 140 NSPX = 1,IBCSP(1,NB)
                  IF( IBCSP(NSPX+1,NB).EQ.IEQ_K(NSP+1,NEQ) ) THEN
                    BCX(1) = BCX(1) + EQ_K(NSP,NEQ)*BCX(NSPX+1)
                  ENDIF
  140           CONTINUE
  150         CONTINUE
            ENDIF
          ENDIF
          SVLB = SLB(2,NB)*PORDB(2,NB)
          YLB(NB,NSL) = 1.D+0
          FCLB = 0.D+0
          IF( SVLB.GT.SMALL ) FCLB = YLB(NB,NSL)/SVLB
!
!---      Convert boundary phase concentrations to
!         volumetric concentrations  ---
!
          IF( IBCT(NBCT,NB).EQ.8 .OR. IBCT(NBCT,NB).EQ.14
     &      .OR. IBCT(NBCT,NB).EQ.23 ) THEN
            BCX(1) = BCX(1)*SVLB
          ENDIF
        ENDIF
!
!---    Diffusion coefficients at nodes adjacent to boundaries  ---
!
        IF( IEDL(NSL).EQ.1 ) THEN
          TCOR = (T(2,N)+TABS)/TSPRF
          PCOR = (PG(2,N)+PATM)/PATM
          SDFGP = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
          DGP = TORG(2,N)*(SG(2,N)-SGT(2,N))*PORD(2,N)*SDFGP
        ENDIF
!
!---    Bottom boundary  ---
!
        IF( IBCD(NB) .EQ. -3 ) THEN
          NPZ = NSZ(N)
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVBB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGBX,VGBX,WGBX,N,MF )
            CALL SHDP( WGBX,UGBX,VGBX,DISPL(IZN),DISPT(IZN),DPGB )
          ELSE
            DPGB = 0.D+0
          ENDIF
!
!---      Dirichlet ---
!
          IF( IBCT(NBCT,NB).EQ.1 .OR. IBCT(NBCT,NB).EQ.8 .OR.
     &      IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.10 .OR.
     &      IBCT(NBCT,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            PCOR = (PGB(2,NB)+PATM)/PATM
            SDFGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            DGB = TORGB(2,NB)*SVGB*SDFGB
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
            WC(NPZ,NSL) = WC(NPZ,NSL) + (BCX(1)*AB*FCGB
     &        - C(N,NSL)*AP*FCGP)
!
!---      Outflow ---
!
          ELSEIF( IBCT(NBCT,NB) .EQ. 7 ) THEN
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 )  THEN
              AP = 0.D+0
            ELSE
              AP = MAX( -WG(1,NPZ),ZERO )
            ENDIF
            WC(NPZ,NSL) = WC(NPZ,NSL) - C(N,NSL)*AP*FCGP
!
!---      Inflow ---
!
          ELSEIF( IBCT(NBCT,NB).GE.13 .AND. IBCT(NBCT,NB).LE.16 ) THEN
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 )  THEN
              AB = 0.D+0
            ELSE
              AB = MAX( WL(1,NPZ),ZERO )
            ENDIF
            WC(NPZ,NSL) = WC(NPZ,NSL) + BCX(1)*AB*FCGB
          ENDIF
!
!---    South boundary  ---
!
        ELSEIF( IBCD(NB) .EQ. -2 ) THEN
          NPY = NSY(N)
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVSB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGSX,VGSX,WGSX,N,MF )
            CALL SHDP( VGSX,WGSX,UGSX,DISPL(IZN),DISPT(IZN),DPGS )
          ELSE
            DPGS = 0.D+0
          ENDIF
!
!---      Dirichlet ---
!
          IF( IBCT(NBCT,NB).EQ.1 .OR. IBCT(NBCT,NB).EQ.8 .OR.
     &      IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.10 .OR.
     &      IBCT(NBCT,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            PCOR = (PGB(2,NB)+PATM)/PATM
            SDFGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            DGB = TORGB(2,NB)*SVGB*SDFGB
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
            VC(NPY,NSL) = VC(NPY,NSL) + (BCX(1)*AS*FCGB
     &        - C(N,NSL)*AP*FCGP)
!
!---      Outflow ---
!
          ELSEIF( IBCT(NBCT,NB) .EQ. 7 ) THEN
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 )  THEN
              AP = 0.D+0
            ELSE
              AP = MAX( -VG(1,NPY),ZERO )
            ENDIF
            VC(NPY,NSL) = VC(NPY,NSL) - C(N,NSL)*AP*FCGP
!
!---      Inflow ---
!
          ELSEIF( IBCT(NBCT,NB).GE.13 .AND. IBCT(NBCT,NB).LE.16 ) THEN
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 )  THEN
              AS = 0.D+0
            ELSE
              AS = MAX( VG(1,NPY),ZERO )
            ENDIF
            VC(NPY,NSL) = VC(NPY,NSL) - BCX(1)*AS*FCGB
          ENDIF
!
!---    West boundary  ---
!
        ELSEIF( IBCD(NB) .EQ. -1 ) THEN
          NPX = NSX(N)
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVWB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGX,VGX,WGX,N,MF )
            CALL SHDP( UGX,VGX,WGX,DISPL(IZN),DISPT(IZN),DPGW )
            ELSE
              DPGW = 0.D+0
            ENDIF
!
!---      Dirichlet ---
!
          IF( IBCT(NBCT,NB).EQ.1 .OR. IBCT(NBCT,NB).EQ.8 .OR.
     &      IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.10 .OR.
     &      IBCT(NBCT,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            PCOR = (PGB(2,NB)+PATM)/PATM
            SDFGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            DGB = TORGB(2,NB)*SVGB*SDFGB
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
            UC(NPX,NSL) = UC(NPX,NSL) + (BCX(1)*AW*FCGB
     &        - C(N,NSL)*AP*FCGP)
!
!---      Outflow ---
!
          ELSEIF( IBCT(NBCT,NB) .EQ. 7 ) THEN
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 )  THEN
              AP = 0.D+0
            ELSE
              AP = MAX( -UG(1,NPX),ZERO )
            ENDIF
            UC(NPX,NSL) = UC(NPX,NSL) -  C(N,NSL)*AP*FCGP
!
!---      Inflow ---
!
          ELSEIF( IBCT(NBCT,NB).GE.13 .AND. IBCT(NBCT,NB).LE.16 ) THEN
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 )  THEN
              AW = 0.D+0
            ELSE
              AW = MAX( UG(1,NPX),ZERO )
            ENDIF
            UC(NPX,NSL) = UC(NPX,NSL) + BCX(1)*AW*FCGB
          ENDIF
!
!---    East boundary  ---
!
        ELSEIF( IBCD(NB) .EQ. 1 ) THEN
          NQX = NSX(N)+1
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVEB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGEX,VGEX,WGEX,N,MF )
            CALL SHDP( UGEX,VGEX,WGEX,DISPL(IZN),DISPT(IZN),DPGE )
          ELSE
            DPGE = 0.D+0
          ENDIF
!
!---      Dirichlet ---
!
          IF( IBCT(NBCT,NB).EQ.1 .OR. IBCT(NBCT,NB).EQ.8 .OR.
     &      IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.10 .OR.
     &      IBCT(NBCT,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            PCOR = (PGB(2,NB)+PATM)/PATM
            SDFGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            DGB = TORGB(2,NB)*SVGB*SDFGB
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
            UC(NQX,NSL) = UC(NQX,NSL)+(C(N,NSL)*AP*FCGP-BCX(1)*AE*FCGB)
!
!---      Outflow ---
!
          ELSEIF( IBCT(NBCT,NB) .EQ. 7 ) THEN
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 )  THEN
              AP = 0.D+0
            ELSE
              AP = MAX( UG(1,NQX),ZERO )
            ENDIF
            UC(NQX,NSL) = UC(NQX,NSL) + C(N,NSL)*AP*FCGP
!
!---      Inflow ---
!
          ELSEIF( IBCT(NBCT,NB).GE.13 .AND. IBCT(NBCT,NB).LE.16 ) THEN
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 )  THEN
              AE = 0.D+0
            ELSE
              AE = MAX( -UG(1,NQX),ZERO )
            ENDIF
            UC(NQX,NSL) = UC(NQX,NSL) - BCX(1)*AE*FCGB
          ENDIF
!
!---    North boundary  ---
!
        ELSEIF( IBCD(NB) .EQ. 2 ) THEN
          NQY = NSY(N)+IFLD
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVNB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGNX,VGNX,WGNX,N,MF )
            CALL SHDP( VGNX,WGNX,UGNX,DISPL(IZN),DISPT(IZN),DPGN )
          ELSE
            DPGN = 0.D+0
          ENDIF
!
!---      Dirichlet ---
!
          IF( IBCT(NBCT,NB).EQ.1 .OR. IBCT(NBCT,NB).EQ.8 .OR.
     &      IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.10 .OR.
     &      IBCT(NBCT,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            PCOR = (PGB(2,NB)+PATM)/PATM
            SDFGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            DGB = TORGB(2,NB)*SVGB*SDFGB
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
            VC(NQY,NSL) = VC(NQY,NSL) + (C(N,NSL)*AP*FCGP
     &        - BCX(1)*AN*FCGB)
!
!---      Outflow ---
!
          ELSEIF( IBCT(NBCT,NB) .EQ. 7 ) THEN
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 )  THEN
              AP = 0.D+0
            ELSE
              AP = MAX( VG(1,NQY),ZERO )
            ENDIF
            VC(NQY,NSL) = VC(NQY,NSL) + C(N,NSL)*AP*FCGP
!
!---      Inflow ---
!
          ELSEIF( IBCT(NBCT,NB).GE.13 .AND. IBCT(NBCT,NB).LE.16 ) THEN
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 )  THEN
              AN = 0.D+0
            ELSE
              AN = MAX( -VG(1,NQY),ZERO )
            ENDIF
            VC(NQY,NSL) = VC(NQY,NSL) - BCX(1)*AN*FCGB
          ENDIF
!
!---    Top boundary
!
        ELSEIF( IBCD(NB) .EQ. 3 ) THEN
          NQZ = NSZ(N)+IJFLD
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVTB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGTX,VGTX,WGTX,N,MF )
            CALL SHDP( WGTX,UGTX,VGTX,DISPL(IZN),DISPT(IZN),DPGT )
          ELSE
            DPGT = 0.D+0
          ENDIF
!
!---      Dirichlet ---
!
          IF( IBCT(NBCT,NB).EQ.1 .OR. IBCT(NBCT,NB).EQ.8 .OR.
     &      IBCT(NBCT,NB).EQ.9 .OR. IBCT(NBCT,NB).EQ.10 .OR.
     &      IBCT(NBCT,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            PCOR = (PGB(2,NB)+PATM)/PATM
            SDFGB = SMDG(NSL)*(TCOR**1.75D+0)/PCOR
            DGB = TORGB(2,NB)*SVGB*SDFGB
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
            WC(NQZ,NSL) = WC(NQZ,NSL) + (C(N,NSL)*AP*FCGP
     &        - BCX(1)*AT*FCGB)
!
!---      Outflow ---
!
          ELSEIF( IBCT(NBCT,NB) .EQ. 7 ) THEN
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 )  THEN
              AP = 0.D+0
            ELSE
              AP = MAX( WG(1,NQZ),ZERO )
            ENDIF
            WC(NQZ,NSL) = WC(NQZ,NSL) + C(N,NSL)*AP*FCGP
!
!---      Inflow ---
!
          ELSEIF( IBCT(NBCT,NB).GE.13 .AND. IBCT(NBCT,NB).LE.16 ) THEN
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 )  THEN
              AT = 0.D+0
            ELSE
              AT = MAX( -WG(1,NQZ),ZERO )
            ENDIF
            WC(NQZ,NSL) = WC(NQZ,NSL) - BCX(1)*AT*FCGB
          ENDIF
        ENDIF
  200 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFXGB group  ---
!
      RETURN
      END



!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFXNB( NSL )
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
!     Compute solute NAPL fluxes on boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1995.
!     Last Modified by MD White, Battelle, PNL, October 13, 1995.
!     sfxnb.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE NAPL
      USE GRID
      USE FLUXP
      USE FLUXN
      USE FDVP
      USE FDVN
      USE CONST
      USE BCVP
      USE BCVN
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
      SUB_LOG(ISUB_LOG) = '/SFXNB'
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
        MF = 1
        IZN = IZ(N)
        MP = IXP(N)
        I = ID(N)
        J = JD(N)
        K = KD(N)
!
!---  Compute adjacent node phase fractions  ---
!
        SVNP = SN(2,N)*PORD(2,N)
        FCNP = 0.D+0
        IF( SVNP.GT.SMALL ) FCNP = YN(N,NSL)/SVNP
!
!---  Compute boundary phase fractions  ---
!
        IF( IPCL(NSL).EQ.2 ) THEN
          SVSB = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORTB(2,NB))*SLB(2,NB)
        ELSE
          SVSB = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORTB(2,NB))
        ENDIF
        SVLB = SLB(2,NB)*PORDB(2,NB)
        SVGB = SGB(2,NB)*PORDB(2,NB)
        SVNB = SNB(2,NB)*PORDB(2,NB)
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
     &      + PCGL(3,NSL)*LOG(TK) + PCGL(4,NSL)*TK
     &      + PCGL(5,NSL)*TK**2 )
!
!---    Water-vapor equilibrium gas-aqueous partition coefficient  ---
!
        ELSEIF( IPCGL(NSL).EQ.2 ) THEN
          PCGLX = RHOG(2,N)*XGW(2,N)/(RHOL(2,N)*XLW(2,N))
        ENDIF
        PCGLX = MAX( PCGLX,1.D-20 )
        PCGLX = MIN( PCGLX,1.D+20 )
!
!---  Phase-volumetric concentration ratios  ---
!
        FCLB = 1.D+0/(SVSB + SVLB + SVNB/PCLN(1,NSL) + SVGB*PCGLX)
        FCGB = 1.D+0/((SVSB + SVLB + SVNB)/PCGLX + SVGB)
        FCNB = 1.D+0/((SVSB + SVLB + SVGB*PCGLX)*PCLN(1,NSL) + SVNB)
!
!---  Phase mole fractions  ---
!
        YLB(NB,NSL) = SVLB*FCLB
        YGB(NB,NSL) = SVGB*FCGB
        YNB(NB,NSL) = SVNB*FCNB
!
!---  Convert boundary concentrations  ---
!
        IF( IBCT(IEQC,NB).EQ.8 .OR. IBCT(IEQC,NB).EQ.14 ) THEN
          BCX = BCX/(FCLB+SMALL)
        ELSEIF( IBCT(IEQC,NB).EQ.9 .OR. IBCT(IEQC,NB).EQ.15 ) THEN
          BCX = BCX/(FCGB+SMALL)
        ELSEIF( IBCT(IEQC,NB).EQ.10 .OR. IBCT(IEQC,NB).EQ.16 ) THEN
          BCX = BCX/(FCNB+SMALL)
        ENDIF
!
!---  Diffusion coefficients at nodes adjacent to boundaries  ---
!
        TCOR = (T(2,N)+TABS)/TSPRF
        SDFNP = SMDN(NSL)*TCOR*(VISRN/VISN(2,N))
        DNP = TORN(2,N)*SN(2,N)*PORD(2,N)*SDFNP
!
!---  Bottom boundary  ---
!
        IF( IBCD(NB) .EQ. -3 ) THEN
          NPZ = NSZ(N)
!
!---  Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVBB( PORD(2,N),PORDB(2,NB),SN(2,N),SNB(2,NB),
     &        UN,VN,WN,UNBX,VNBX,WNBX,N,MF )
            CALL SHDP( WNBX,UNBX,VNBX,DISPL(IZN),DISPT(IZN),DPNB )
!            UNBX = (0.5D+0*(UN(1,NSX(N))+UN(1,NSX(N)+1)))**2
!            VNBX = (0.5D+0*(VN(1,NSY(N))+VN(1,NSY(N)+IFLD)))**2
!            WNBX = (WN(1,NPZ))**2
!            ZNB = SQRT(UNBX + VNBX + WNBX)
!            DPNB = (DISPL(IZN)*WNBX + DISPT(IZN)*(UNBX+VNBX))/
!     &        (ZNB+SMALL)
          ELSE
            DPNB = 0.D+0
          ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &      IBCT(IEQC,NB).EQ.9 .OR. IBCT(IEQC,NB).EQ.10 .OR.
     &      IBCT(IEQC,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SDFNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*SVNB*SDFNB
            INDX = 16
            DNZ = DIFMN(DNB,DNP,DZGF(N),DZGF(N),WN(1,NPZ),INDX)
            DNZ = (DNZ+DPNB)/(5.D-1*DZGF(N))
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 )  THEN
              AB = DNZ
              AP = DNZ
            ELSE
              AB = MAX( WN(1,NPZ),ZERO ) +
     &         DNZ*MAX((ONE-(TENTH*ABS(WN(1,NPZ))/(DNZ+SMALL)))**5,ZERO)
              AP = MAX( -WN(1,NPZ),ZERO ) +
     &         DNZ*MAX((ONE-(TENTH*ABS(WN(1,NPZ))/(DNZ+SMALL)))**5,ZERO)
            ENDIF
            WC(NPZ,NSL) = WC(NPZ,NSL) + (BCX*AB*FCNB - C(N,NSL)*AP*FCNP)
!
!---  Outflow ---
!
          ELSEIF( IBCT(IEQC,NB) .EQ. 7 ) THEN
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 )  THEN
              AP = 0.D+0
            ELSE
              AP = MAX( -WN(1,NPZ),ZERO )
            ENDIF
            WC(NPZ,NSL) = WC(NPZ,NSL) - C(N,NSL)*AP*FCNP
!
!---  Inflow ---
!
          ELSEIF( IBCT(IEQC,NB).GE.13 .AND. IBCT(IEQC,NB).LE.16 ) THEN
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 )  THEN
              AB = 0.D+0
            ELSE
              AB = MAX( WL(1,NPZ),ZERO )
            ENDIF
            WC(NPZ,NSL) = WC(NPZ,NSL) + BCX*AB*FCNP
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
            CALL ADVSB( PORD(2,N),PORDB(2,NB),SN(2,N),SNB(2,NB),
     &        UN,VN,WN,UNSX,VNSX,WNSX,N,MF )
            CALL SHDP( VNSX,WNSX,UNSX,DISPL(IZN),DISPT(IZN),DPNS )
!            UNSX = (0.5D+0*(UN(1,NSX(N))+UN(1,NSX(N)+1)))**2
!            VNSX = VN(1,NPY)**2
!            WNSX = (0.5D+0*(WN(1,NSZ(N))+WN(1,NSZ(N)+IJFLD)))**2
!            ZNS = SQRT(UNSX + VNSX + WNSX)
!            DPNS = (DISPL(IZN)*VNSX + DISPT(IZN)*(UNSX+WNSX))/
!     &        (ZNS+SMALL)
          ELSE
            DPNS = 0.D+0
          ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &      IBCT(IEQC,NB).EQ.9 .OR. IBCT(IEQC,NB).EQ.10 .OR.
     &      IBCT(IEQC,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SDFNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*SVNB*SDFNB
            INDX = 16
            DNY = DIFMN(DNB,DNP,DYGF(N),DYGF(N),VN(1,NPY),INDX)
            DNY = (DNY+DPNS)/RP(I)/(5.D-1*DYGF(N))
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 )  THEN
              AS = DNY
              AP = DNY
            ELSE
              AS = MAX( VN(1,NPY),ZERO ) +
     &         DNY*MAX((ONE-(TENTH*ABS(VN(1,NPY))/(DNY+SMALL)))**5,ZERO)
              AP = MAX( -VN(1,NPY),ZERO ) +
     &         DNY*MAX((ONE-(TENTH*ABS(VN(1,NPY))/(DNY+SMALL)))**5,ZERO)
            ENDIF
            VC(NPY,NSL) = VC(NPY,NSL) + (BCX*AS*FCNP - C(N,NSL)*AP*FCNP)
!
!---  Outflow ---
!
          ELSEIF( IBCT(IEQC,NB) .EQ. 7 ) THEN
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 )  THEN
              AP = 0.D+0
            ELSE
              AP = MAX( -VN(1,NPY),ZERO )
            ENDIF
            VC(NPY,NSL) = VC(NPY,NSL) - C(N,NSL)*AP*FCNP
!
!---  Inflow ---
!
          ELSEIF( IBCT(IEQC,NB).GE.13 .AND. IBCT(IEQC,NB).LE.16 ) THEN
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 )  THEN
              AS = 0.D+0
            ELSE
              AS = MAX( VN(1,NPY),ZERO )
            ENDIF
            VC(NPY,NSL) = VC(NPY,NSL) - BCX*AS*FCNP
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
            CALL ADVWB( PORD(2,N),PORDB(2,NB),SN(2,N),SNB(2,NB),
     &        UN,VN,WN,UNWX,VNWX,WNWX,N,MF )
            CALL SHDP( UNWX,VNWX,WNWX,DISPL(IZN),DISPT(IZN),DPNW )
!            UNWX = UN(1,NPX)**2
!            VNWX = (0.5D+0*(VN(1,NSY(N))+VN(1,NSY(N)+IFLD)))**2
!            WNWX = (0.5D+0*(WN(1,NSZ(N))+WN(1,NSZ(N)+IJFLD)))**2
!            ZNW = SQRT(UNWX + VNWX + WNWX)
!            DPNW = (DISPL(IZN)*UNWX + DISPT(IZN)*(WNWX+VNWX))/
!     &          (ZNW+SMALL)
            ELSE
              DPNW = 0.D+0
            ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &      IBCT(IEQC,NB).EQ.9 .OR. IBCT(IEQC,NB).EQ.10 .OR.
     &      IBCT(IEQC,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SDFNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*SVNB*SDFNB
            INDX = 16
            DNX = DIFMN(DNB,DNP,DXGF(N),DXGF(N),UN(1,NPX),INDX)
            DNX = (DNX+DPNW)/(5.D-1*DXGF(N))
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 )  THEN
              AW = DNX
              AP = DNX
            ELSE
              AW = MAX( UN(1,NPX),ZERO ) +
     &         DNX*MAX((ONE-(TENTH*ABS(UN(1,NPX))/(DNX+SMALL)))**5,ZERO)
              AP = MAX( -UN(1,NPX),ZERO ) +
     &         DNX*MAX((ONE-(TENTH*ABS(UN(1,NPX))/(DNX+SMALL)))**5,ZERO)
            ENDIF
            UC(NPX,NSL) = UC(NPX,NSL) + (BCX*AW*FCNP - C(N,NSL)*AP*FCNP)
!
!---  Outflow ---
!
          ELSEIF( IBCT(IEQC,NB) .EQ. 7 ) THEN
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 )  THEN
              AP = 0.D+0
            ELSE
              AP = MAX( -UN(1,NPX),ZERO )
            ENDIF
            UC(NPX,NSL) = UC(NPX,NSL) -  C(N,NSL)*AP*FCNP
!
!---  Inflow ---
!
          ELSEIF( IBCT(IEQC,NB).GE.13 .AND. IBCT(IEQC,NB).LE.16 ) THEN
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 )  THEN
              AW = 0.D+0
            ELSE
              AW = MAX( UN(1,NPX),ZERO )
            ENDIF
            UC(NPX,NSL) = UC(NPX,NSL) + BCX*AW*FCNP
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
            CALL ADVEB( PORD(2,N),PORDB(2,NB),SN(2,N),SNB(2,NB),
     &        UN,VN,WN,UNEX,VNEX,WNEX,N,MF )
            CALL SHDP( UNEX,VNEX,WNEX,DISPL(IZN),DISPT(IZN),DPNE )
!            UNEX = UN(1,NQX)**2
!            VNEX = (0.5D+0*(VN(1,NSY(N))+VN(1,NSY(N)+IFLD)))**2
!            WNEX = (0.5D+0*(WN(1,NSZ(N))+WN(1,NSZ(N)+IJFLD)))**2
!            ZNE = SQRT(UNEX + VNEX + WNEX)
!            DPNE = (DISPL(IZN)*UNEX + DISPT(IZN)*(WNEX+VNEX))/
!     &        (ZNE+SMALL)
          ELSE
            DPNE = 0.D+0
          ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &      IBCT(IEQC,NB).EQ.9 .OR. IBCT(IEQC,NB).EQ.10 .OR.
     &      IBCT(IEQC,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SDFNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*SVNB*SDFNB
            INDX = 16
            DNX = DIFMN(DNP,DNB,DXGF(N),DXGF(N),UN(1,NQX),INDX)
            DNX = (DNX+DPNE)/(5.D-1*DXGF(N))
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 )  THEN
              AE = DNX
              AP = DNX
            ELSE
              AE = MAX( -UN(1,NQX),ZERO ) +
     &         DNX*MAX((ONE-(TENTH*ABS(UN(1,NQX))/(DNX+SMALL)))**5,ZERO)
              AP = MAX( UN(1,NQX),ZERO ) +
     &         DNX*MAX((ONE-(TENTH*ABS(UN(1,NQX))/(DNX+SMALL)))**5,ZERO)
            ENDIF
            UC(NQX,NSL) = UC(NQX,NSL)+(C(N,NSL)*AP*FCNP-BCX*AE*FCNP)
!
!---  Outflow ---
!
          ELSEIF( IBCT(IEQC,NB) .EQ. 7 ) THEN
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 )  THEN
              AP = 0.D+0
            ELSE
              AP = MAX( UN(1,NQX),ZERO )
            ENDIF
            UC(NQX,NSL) = UC(NQX,NSL) + C(N,NSL)*AP*FCNP
!
!---  Inflow ---
!
          ELSEIF( IBCT(IEQC,NB).GE.13 .AND. IBCT(IEQC,NB).LE.16 ) THEN
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 )  THEN
              AE = 0.D+0
            ELSE
              AE = MAX( -UN(1,NQX),ZERO )
            ENDIF
            UC(NQX,NSL) = UC(NQX,NSL) - BCX*AE*FCNP
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
            CALL ADVNB( PORD(2,N),PORDB(2,NB),SN(2,N),SNB(2,NB),
     &        UN,VN,WN,UNNX,VNNX,WNNX,N,MF )
            CALL SHDP( VNNX,WNNX,UNNX,DISPL(IZN),DISPT(IZN),DPNN )
!            UNNX = (0.5D+0*(UN(1,NSX(N))+UN(1,NSX(N)+1)))**2
!            VNNX = VN(1,NQY)**2
!            WNNX = (0.5D+0*(WN(1,NSZ(N))+WN(1,NSZ(N)+IJFLD)))**2
!            ZNN = SQRT(UNNX + VNNX + WNNX)
!            DPNN = (DISPL(IZN)*VNNX + DISPT(IZN)*(UNNX+WNNX))/
!     &        (ZNN+SMALL)
          ELSE
            DPNN = 0.D+0
          ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &      IBCT(IEQC,NB).EQ.9 .OR. IBCT(IEQC,NB).EQ.10 .OR.
     &      IBCT(IEQC,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SDFNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*SVNB*SDFNB
            INDX = 16
            DNY = DIFMN(DNP,DNB,DYGF(N),DYGF(N),VN(1,NQY),INDX)
            DNY = (DNY+DPNN)/RP(I)/(5.D-1*DYGF(N))
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 )  THEN
              AN = DNY
              AP = DNY
            ELSE
              AN = MAX( -VN(1,NQY),ZERO ) +
     &         DNY*MAX((ONE-(TENTH*ABS(VN(1,NQY))/(DNY+SMALL)))**5,ZERO)
              AP = MAX( VN(1,NQY),ZERO ) +
     &         DNY*MAX((ONE-(TENTH*ABS(VN(1,NQY))/(DNY+SMALL)))**5,ZERO)
            ENDIF
            VC(NQY,NSL) = VC(NQY,NSL) + (C(N,NSL)*AP*FCNP - BCX*AN*FCNP)
!
!---  Outflow ---
!
          ELSEIF( IBCT(IEQC,NB) .EQ. 7 ) THEN
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 )  THEN
              AP = 0.D+0
            ELSE
              AP = MAX( VN(1,NQY),ZERO )
            ENDIF
            VC(NQY,NSL) = VC(NQY,NSL) + C(N,NSL)*AP*FCNP
!
!---  Inflow ---
!
          ELSEIF( IBCT(IEQC,NB).GE.13 .AND. IBCT(IEQC,NB).LE.16 ) THEN
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 )  THEN
              AN = 0.D+0
            ELSE
              AN = MAX( -VN(1,NQY),ZERO )
            ENDIF
            VC(NQY,NSL) = VC(NQY,NSL) - BCX*AN*FCNP
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
            CALL ADVTB( PORD(2,N),PORDB(2,NB),SN(2,N),SNB(2,NB),
     &        UN,VN,WN,UNTX,VNTX,WNTX,N,MF )
            CALL SHDP( WNTX,UNTX,VNTX,DISPL(IZN),DISPT(IZN),DPNT )
!            UNTX = (0.5D+0*(UN(1,NSX(N))+UN(1,NSX(N)+1)))**2
!            VNTX = (0.5D+0*(VN(1,NSY(N))+VN(1,NSY(N)+IFLD)))**2
!            WNTX = (WN(1,NQZ))**2
!            ZNT = SQRT(UNTX + VNTX + WNTX)
!            DPNT = (DISPL(IZN)*WNTX + DISPT(IZN)*(UNTX+VNTX))/
!     &        (ZNT+SMALL)
          ELSE
            DPNT = 0.D+0
          ENDIF
!
!---  Dirichlet ---
!
          IF( IBCT(IEQC,NB).EQ.1 .OR. IBCT(IEQC,NB).EQ.8 .OR.
     &      IBCT(IEQC,NB).EQ.9 .OR. IBCT(IEQC,NB).EQ.10 .OR.
     &      IBCT(IEQC,NB).EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SDFNB = SMDN(NSL)*TCOR*(VISRN/VISNB(2,NB))
            DNB = TORNB(2,NB)*SVNB*SDFNB
            INDX = 16
            DNZ = DIFMN(DNP,DNB,DZGF(N),DZGF(N),WN(1,NQZ),INDX)
            DNZ = (DNZ+DPNT)/(5.D-1*DZGF(N))
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 )  THEN
              AT = DNZ
              AP = DNZ
            ELSE
              AT = MAX( -WN(1,NQZ),ZERO ) +
     &         DNZ*MAX((ONE-(TENTH*ABS(WN(1,NQZ))/(DNZ+SMALL)))**5,ZERO)
              AP = MAX( WN(1,NQZ),ZERO ) +
     &         DNZ*MAX((ONE-(TENTH*ABS(WN(1,NQZ))/(DNZ+SMALL)))**5,ZERO)
            ENDIF
            WC(NQZ,NSL) = WC(NQZ,NSL) + (C(N,NSL)*AP*FCNP - BCX*AT*FCNP)
!
!---  Outflow ---
!
          ELSEIF( IBCT(IEQC,NB) .EQ. 7 ) THEN
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 )  THEN
              AP = 0.D+0
            ELSE
              AP = MAX( WN(1,NQZ),ZERO )
            ENDIF
            WC(NQZ,NSL) = WC(NQZ,NSL) + C(N,NSL)*AP*FCNP
!
!---  Inflow ---
!
          ELSEIF( IBCT(IEQC,NB).GE.13 .AND. IBCT(IEQC,NB).LE.16 ) THEN
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 )  THEN
              AT = 0.D+0
            ELSE
              AT = MAX( -WN(1,NQZ),ZERO )
            ENDIF
            WC(NQZ,NSL) = WC(NQZ,NSL) - BCX*AT*FCNP
          ENDIF
        ENDIF
  200 CONTINUE
!
!---  End of SFXNB group  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END


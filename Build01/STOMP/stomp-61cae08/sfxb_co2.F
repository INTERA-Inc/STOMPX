!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFXB32( NSL )
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
!     Compute solute aqueous-phase fluxes on boundary surfaces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNNL, January, 2006.
!     Last Modified by MD White, Battelle, PNNL, 18 January 2006.
!     sfxb_co2.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/SFXB32'
!
!---  Loop over number of specified boundary conditions  ---
!
      DO 200 NB = 1,NBC
        TMZ = TM
        MB = IBCIN(NB)
        IF( IBCC(NB).EQ.1 ) TMZ = MOD( TM,BC(1,IBCM(NB),MB) )
        IF( TMZ.LE.BC(1,1,MB) ) GOTO 200
        IF( IBCM(NB).GT.1 .AND. TMZ.GT.BC(1,IBCM(NB),MB) ) GOTO 200
!
!---    Solute transport  ---
!
        IF( NSL.LE.NSOLU ) THEN
          IBCTX = IBCT(NSL+LUK,NB)
!
!---    Reactive species transport  ---
!
        ELSE
          IBCTX = IBCT(NSOLU+LUK+1,NB)
        ENDIF
!
!---    Single boundary condition time  ---
!
        IF( IBCM(NB).EQ.1 ) THEN
!
!---      Solute transport  ---
!
          IF( NSL.LE.NSOLU ) THEN
            BCX(1) = BC(NSL+LBCU,1,MB)
            IF( IBCTX.EQ.12 ) BCX(1) = CBO(NB,NSL)
!
!---      Reactive species transport  ---
!
          ELSE
            BCX(1) = 0.D+0
            DO 10 NSPX = 1,IBCSP(1,NB)
              NSP = IBCSP(NSPX+1,NB)
              MX = NSOLU+LBCU+NSPX
              BCX(NSPX+1) = BC(MX,1,MB)
!
!---          Aqueous species ---
!
              IF( NSP.LE.NSPL ) THEN
                IF( IBCT(NSOLU+LUK+1,NB).EQ.12 ) 
     &            BCX(NSPX+1) = SP_CBO(NB,NSP)
!
!---          Gas species ---
!
              ELSE
                IF( IBCT(NSOLU+LUK+2,NB).EQ.12 ) 
     &            BCX(NSPX+1) = SP_CBO(NB,NSP)
              ENDIF
   10       CONTINUE
          ENDIF
!
!---    Multiple boundary condition times  ---
!
        ELSE
          DO 100 M = 2,IBCM(NB)
            IF( TMZ.LE.BC(1,M,MB) ) THEN
              TDBC = (BC(1,M,MB)-BC(1,M-1,MB))
              DTBC = MIN( BC(1,M,MB)-TMZ,DT )
              TFBC = (TMZ-5.D-1*DTBC-BC(1,M-1,MB))/TDBC
!
!---          Solute transport  ---
!
              IF( NSL.LE.NSOLU ) THEN
                BCX(1) = BC(NSL+LBCU,M-1,MB) +
     &            TFBC*(BC(NSL+LBCU,M,MB)-BC(NSL+LBCU,M-1,MB))
                IF( IBCT(NSL+LUK,NB).EQ.12 ) BCX(1) = CBO(NB,NSL)
!
!---          Reactive species transport  ---
!
              ELSE
                BCX(1) = 0.D+0
                DO 20 NSPX = 1,IBCSP(1,NB)
                  NSP = IBCSP(NSPX+1,NB)
                  MX = NSOLU+LBCU+NSPX
                  BCX(NSPX+1) = BC(MX,M-1,MB) +
     &              TFBC*(BC(MX,M,MB)-BC(MX,M-1,MB))
!
!---              Aqueous species ---
!
                  IF( NSP.LE.NSPL ) THEN
                    IF( IBCT(NSOLU+LUK+1,NB).EQ.12 ) 
     &                BCX(NSPX+1) = SP_CBO(NB,NSP)
!
!---              Gas species ---
!
                  ELSE
                    IF( IBCT(NSOLU+LUK+2,NB).EQ.12 ) 
     &                BCX(NSPX+1) = SP_CBO(NB,NSP)
                  ENDIF
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
!---    Diffusion coefficients at node adjacent to boundary  ---
!
        TCOR = (T(2,N)+TABS)/TSPRF
        SMDLP = SMDL(NSL)*TCOR*(VISRL/VISL(2,N))
        DLP = TORL(2,N)*SL(2,N)*PORD(2,N)*SMDLP
        PCOR = (PG(2,N)+PATM)/PATM
        SMDGP = SMDG(NSL)*(TCOR**1.75)/PCOR
        DGP = TORG(2,N)*(SG(2,N)-SGT(2,N))*PORD(2,N)*SMDGP
!
!---    Phase fraction factors at node adjacent to boundary  ---
!
        XVLP = SL(2,N)*PORD(2,N)
        FCLP = 0.D+0
        IF( XVLP.GT.SMALL ) FCLP = YL(N,NSL)/XVLP
        XVGP = SG(2,N)*PORD(2,N)
        FCGP = 0.D+0
        IF( XVGP.GT.SMALL ) FCGP = YG(N,NSL)/XVGP
!
!---    Solute transport only, skip calculations for reactive
!       species transport  ---
!
        XVLB = SLB(2,NB)*PORDB(2,NB)
        XVGB = SGB(2,NB)*PORDB(2,NB)
        IF( NSL.LE.NSOLU ) THEN
!
!---      Compute boundary phase fractions  ---
!
          IF( IPCL(NSL).EQ.2 ) THEN
            SVSB = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORTB(2,NB))*
     &        SLB(2,NB)
          ELSE
            SVSB = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORTB(2,NB))
          ENDIF
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
          FCL = 1.D+0/(SVSB + XVLB + XVGB*PCGLX)
          FCG = 1.D+0/((SVSB + XVLB)/PCGLX + XVGB)
!
!---      Phase mole fractions  ---
!
          YLB(NB,NSL) = XVLB*FCL
          YGB(NB,NSL) = XVGB*FCG
!
!---      Convert boundary phase concentrations to
!         volumetric concentrations  ---
!
          IF( IBCT(NSL+LUK,NB).EQ.8 .OR.
     &      IBCT(NSL+LUK,NB).EQ.14 .OR.
     &      IBCT(NSL+LUK,NB).EQ.23 ) THEN
            BCX(1) = BCX(1)/(FCL+SMALL)
          ELSEIF( IBCT(NSL+LUK,NB).EQ.9 .OR.
     &      IBCT(NSL+LUK,NB).EQ.15 .OR.
     &      IBCT(NSL+LUK,NB).EQ.43 ) THEN
            BCX(1) = BCX(1)/(FCG+SMALL)
          ENDIF
#ifdef ecke
        ELSE
!
!---      Convert species concentrations to total-component
!         concentrations  ---
!
          IF( NSL.LE.NSOLU+NEQC ) THEN
            NEQ = NSL-NSOLU
            YSPLX = 0.D+0
            DO 130 NSP = 1,IEQ_C(1,NEQ)
              DO 120 NSPX = 1,IBCSP(1,NB)
                IF( IBCSP(NSPX+1,NB).EQ.IEQ_C(NSP+1,NEQ) ) THEN
!
!---              Aqueous species ---
!
                  IF( NSP.LE.NSPL ) THEN
                    IF( IBCT(NSOLU+LUK+1,NB).EQ.8 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.14 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.23 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVLB
                    ELSEIF( IBCT(NSOLU+LUK+1,NB).EQ.9 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.15 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.43 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVGB
                    ENDIF
                    YSPLX = YSPLX + EQ_C(NSP,NEQ)*BCX(NSPX+1)
!
!---              Gas species ---
!
                  ELSE
                    IF( IBCT(NSOLU+LUK+2,NB).EQ.8 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.14 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.23 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVLB
                    ELSEIF( IBCT(NSOLU+LUK+2,NB).EQ.9 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.15 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.43 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVGB
                    ENDIF
                  ENDIF                    
                  BCX(1) = BCX(1) + EQ_C(NSP,NEQ)*BCX(NSPX+1)
                ENDIF
  120         CONTINUE
  130       CONTINUE
!
!---        Convert species concentrations to total-kinetic
!           concentrations  ---
!
          ELSEIF( NSL.LE.NSOLU+NEQC+NEQK ) THEN
            NEQ = NSL-NSOLU-NEQC
            YSPLX = 0.D+0
            DO 150 NSP = 1,IEQ_K(1,NEQ)
              DO 140 NSPX = 1,IBCSP(1,NB)
                IF( IBCSP(NSPX+1,NB).EQ.IEQ_K(NSP+1,NEQ) ) THEN
!
!---              Aqueous species ---
!
                  IF( NSP.LE.NSPL ) THEN
                    IF( IBCT(NSOLU+LUK+1,NB).EQ.8 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.14 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.23 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVLB
                    ELSEIF( IBCT(NSOLU+LUK+1,NB).EQ.9 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.15 .OR.
     &                IBCT(NSOLU+LUK+1,NB).EQ.43 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVGB
                    ENDIF
                    YSPLX = YSPLX + EQ_K(NSP,NEQ)*BCX(NSPX+1)
!
!---              Gas species ---
!
                  ELSE
                    IF( IBCT(NSOLU+LUK+2,NB).EQ.8 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.14 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.23 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVLB
                    ELSEIF( IBCT(NSOLU+LUK+2,NB).EQ.9 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.15 .OR.
     &                IBCT(NSOLU+LUK+2,NB).EQ.43 ) THEN
                      BCX(NSPX+1) = BCX(NSPX+1)*XVGB
                    ENDIF
                  ENDIF                    
                  BCX(1) = BCX(1) + EQ_K(NSP,NEQ)*BCX(NSPX+1)
                ENDIF
  140         CONTINUE
  150       CONTINUE
          ENDIF
          IF( ABS(BCX(1))/EPSL.LT.EPSL ) THEN
            YSPLX = 0.D+0
          ELSE
            YSPLX = YSPLX/BCX(1)
          ENDIF
!
!---      Phase-volumetric concentration ratios  ---
!
          YLBX = MAX( MIN( 1.D+0,YSPLX ),0.D+0 )
          YGBX = 1.D+0-YLBX
          FCL = 0.D+0
          IF( XVLB/EPSL.GT.EPSL ) FCL = YLBX/XVLB
          FCG = 0.D+0
          IF( XVGB/EPSL.GT.EPSL ) FCG = YGBX/XVGB
!
!---      Phase mole fractions  ---
!
          YLB(NB,NSL) = YLBX
          YGB(NB,NSL) = YGBX
#endif
        ENDIF
!
!---    Bottom boundary  ---
!
        IF( IBCD(NB).EQ.-3 ) THEN
          NPZ = NSZ(N)
          NQZ = NPZ+IJFLD
          IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVBB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULBX,VLBX,WLBX,N,MF )
            CALL SHDP( WLBX,ULBX,VLBX,DISPL(IZN),DISPT(IZN),DPLB )
            CALL ADVBB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGBX,VGBX,WGBX,N,MF )
            CALL SHDP( WGBX,UGBX,VGBX,DISPL(IZN),DISPT(IZN),DPGB )
          ELSE
            DPLB = 0.D+0
            DPGB = 0.D+0
          ENDIF
          FLB = AFZ(NPZ)*WL(1,NPZ)
          FGB = AFZ(NPZ)*WG(1,NPZ)
          CRLB = ABS( WL(1,NPZ) )*DT/(DZGF(N)*XVLB+SMALL)
          CRGB = ABS( WG(1,NPZ) )*DT/(DZGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            INDX = 16
            DLZ = DIFMN(DLB,DLP,DZGF(N),DZGF(N),WL(1,NPZ),INDX)
            DLZ = AFZ(NPZ)*(DLZ+DPLB)/(5.D-1*DZGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGZ = DIFMN(DGB,DGP,DZGF(N),DZGF(N),WG(1,NPZ),INDX)
            DGZ = AFZ(NPZ)*(DGZ+DPGB)/(5.D-1*DZGF(N))
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 )  THEN
              IF( FLB.GE.ZERO ) THEN
                WCLZ = BCX(1)*FCL*FLB
              ELSEIF( FLB.LT.ZERO .AND. K.LT.KFLD ) THEN
                NBT = N+IJFLD
                FCLT = YL(NBT,NSL)/(SL(2,NBT)*PORD(2,NBT)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBT,NSL)*FCLT)
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
                THETA = FLIMIT( R,CRLB,ISLC(1) )
                WCLZ = BCX(1)*FLB*THETA*FCL
     &            + C(N,NSL)*FLB*(1.D+0-THETA)*FCLP
              ELSEIF( FLB.LT.ZERO .AND. K.EQ.KFLD ) THEN
                WCLZ = C(N,NSL)*FLB*FCLP
              ENDIF
              IF( FGB.GE.ZERO ) THEN
                WCGZ = BCX(1)*FCG*FGB
              ELSEIF( FGB.LT.ZERO .AND. K.LT.KFLD ) THEN
                NBT = N+IJFLD
                FCGT = YG(NBT,NSL)/(SG(2,NBT)*PORD(2,NBT)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBT,NSL)*FCGT)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
                THETA = FLIMIT( R,CRGB,ISLC(1) )
                WCGZ = BCX(1)*FGB*THETA*FCG
     &            + C(N,NSL)*FGB*(1.D+0-THETA)*FCGP
              ELSEIF( FGB.LT.ZERO .AND. K.EQ.KFLD ) THEN
                WCGZ = C(N,NSL)*FGB*FCGP
              ENDIF
              WC(NPZ,NSL) = WC(NPZ,NSL) + (WCLZ+WCGZ)/AFZ(NPZ)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FLT = AFZ(NQZ)*WL(1,NQZ)
              IF( FLT.GE.ZERO ) THEN
                NBT = N+IJFLD
                XVLX = SL(2,NBT)*PORD(2,NBT)
                FCLT = YL(NBT,NSL)/(XVLX+SMALL)
                CRLT = ABS( WL(1,NQZ) )*DT/DZGP(NQZ)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
     &            /(C(NBT,NSL)*FCLT-C(N,NSL)*FCLP+SMALL))
     &            *((DZGF(NBT)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRLT,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBT))
                WCLZ = C(N,NSL)*FLT*(1.D+0-THETA*DZF)*FCLP
     &            + C(NBT,NSL)*FLT*THETA*DZF*FCLT
                WCLZF = CO(N,NSL)*FLT*FCLP
                WC(NQZ,NSL) = WC(NQZ,NSL) + (WCLZ-WCLZF)/AFZ(NQZ)
              ENDIF
              FGT = AFZ(NQZ)*WG(1,NQZ)
              IF( FGT.GE.ZERO ) THEN
                NBT = N+IJFLD
                XVGX = SG(2,NBT)*PORD(2,NBT)
                FCGT = YG(NBT,NSL)/(XVGX+SMALL)
                CRGT = ABS( WG(1,NQZ) )*DT/DZGP(NQZ)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
     &            /(C(NBT,NSL)*FCGT-C(N,NSL)*FCGP+SMALL))
     &            *((DZGF(NBT)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRGT,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBT))
                WCGZ = C(N,NSL)*FGT*(1.D+0-THETA*DZF)*FCGP
     &            + C(NBT,NSL)*FGT*THETA*DZF*FCGT
                WCGZF = CO(N,NSL)*FGT*FCGP
                WC(NQZ,NSL) = WC(NQZ,NSL) + (WCGZ-WCGZF)/AFZ(NQZ)
              ENDIF
            ELSE
              ALB = MAX( FLB,ZERO ) +
     &          DLZ*MAX((ONE-(TENTH*ABS(FLB)/(DLZ+SMALL)))**5,ZERO)
              AGB = MAX( FGB,ZERO ) +
     &          DGZ*MAX((ONE-(TENTH*ABS(FGB)/(DGZ+SMALL)))**5,ZERO)
              AP = (ALB-FLB)*FCLP + (AGB-FGB)*FCGP
              AB = ALB*FCL + AGB*FCG
              WC(NPZ,NSL) = WC(NPZ,NSL) + (BCX(1)*AB - C(N,NSL)*AP)
            ENDIF
!
!---      Outflow aqueous  ---
!
          ELSEIF( IBCTX.EQ.7 .OR.
     &      ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (WL(1,NPZ)/EPSL.LT.-EPSL)) ) THEN
!
!---      TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 )  THEN
              WCLZ = 0.D+0
              IF( FLB.LT.ZERO .AND. K.LT.KFLD ) THEN
                NBT = N+IJFLD
                FCLT = YL(NBT,NSL)/(SL(2,NBT)*PORD(2,NBT)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBT,NSL)*FCLT)
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
                THETA = FLIMIT( R,CRLB,ISLC(1) )
                WCLZ = BCX(1)*FLB*THETA*FCL
     &            + C(N,NSL)*FLB*(1.D+0-THETA)*FCLP
              ELSEIF( FLB.LT.ZERO .AND. K.EQ.KFLD ) THEN
                WCLZ = C(N,NSL)*FLB*FCLP
              ENDIF
              WC(NPZ,NSL) = WC(NPZ,NSL) + WCLZ/AFZ(NPZ)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FLT = AFZ(NQZ)*WL(1,NQZ)
              IF( FLT.GE.ZERO ) THEN
                NBT = N+IJFLD
                XVLX = SL(2,NBT)*PORD(2,NBT)
                FCLT = YL(NBT,NSL)/(XVLX+SMALL)
                CRLT = ABS( WL(1,NQZ) )*DT/DZGP(NQZ)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
     &            /(C(NBT,NSL)*FCLT-C(N,NSL)*FCLP+SMALL))
     &            *((DZGF(NBT)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRLT,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBT))
                WCLZ = C(N,NSL)*FLT*(1.D+0-THETA*DZF)*FCLP
     &            + C(NBT,NSL)*FLT*THETA*DZF*FCLT
                WCLZF = CO(N,NSL)*FLT*FCLP
                WC(NQZ,NSL) = WC(NQZ,NSL) + (WCLZ-WCLZF)/AFZ(NQZ)
              ENDIF
            ELSE
              ALB = MAX( FLB,ZERO )
              AP = (ALB-FLB)*FCLP
              AB = ALB*FCL
              WC(NPZ,NSL) = WC(NPZ,NSL) + (BCX(1)*AB - C(N,NSL)*AP)
            ENDIF
!
!---      Outflow gas ---
!
          ELSEIF( IBCTX.EQ.7 .OR.
     &      ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (WG(1,NPZ)/EPSL.LT.-EPSL)) ) THEN
!
!---      TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 )  THEN
              WCGZ = 0.D+0
              IF( FGB.LT.ZERO .AND. K.LT.KFLD ) THEN
                NBT = N+IJFLD
                FCGT = YG(NBT,NSL)/(SG(2,NBT)*PORD(2,NBT)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBT,NSL)*FCGT)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBT)))
                THETA = FLIMIT( R,CRGB,ISLC(1) )
                WCGZ = BCX(1)*FGB*THETA*FCG
     &            + C(N,NSL)*FGB*(1.D+0-THETA)*FCGP
              ELSEIF( FGB.LT.ZERO .AND. K.EQ.KFLD ) THEN
                WCGZ = C(N,NSL)*FGB*FCGP
              ENDIF
              WC(NPZ,NSL) = WC(NPZ,NSL) + WCGZ/AFZ(NPZ)
!
!---          TVD Transport for interior surface adjacent
!             to boundary  ---
!
              FGT = AFZ(NQZ)*WG(1,NQZ)
              IF( FGT.GE.ZERO ) THEN
                NBT = N+IJFLD
                XVGX = SG(2,NBT)*PORD(2,NBT)
                FCGT = YG(NBT,NSL)/(XVGX+SMALL)
                CRGT = ABS( WG(1,NQZ) )*DT/DZGP(NQZ)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
     &            /(C(NBT,NSL)*FCGT-C(N,NSL)*FCGP+SMALL))
     &            *((DZGF(NBT)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRGT,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBT))
                WCGZ = C(N,NSL)*FGT*(1.D+0-THETA*DZF)*FCGP
     &            + C(NBT,NSL)*FGT*THETA*DZF*FCGT
                WCGZF = CO(N,NSL)*FGT*FCGP
                WC(NQZ,NSL) = WC(NQZ,NSL) + (WCGZ-WCGZF)/AFZ(NQZ)
              ENDIF
            ELSE
              AGB = MAX( FGB,ZERO )
              AP = (AGB-FGB)*FCGP
              AB = AGB*FCG
              WC(NPZ,NSL) = WC(NPZ,NSL) + (BCX(1)*AB - C(N,NSL)*AP)
            ENDIF
!
!---      Inflow aqueous ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (WL(1,NPZ)/EPSL.GT.EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 )  THEN
              WCLZ = 0.D+0
              IF( FLB.GE.ZERO ) WCLZ = BCX(1)*FCL*FLB
              WC(NPZ,NSL) = WC(NPZ,NSL) + WCLZ/AFZ(NPZ)
!
!---          TVD Transport for interior surface 
!             adjacent to boundary  ---
!
              FLT = AFZ(NQZ)*WL(1,NQZ)
              IF( FLT.GE.ZERO ) THEN
                NBT = N+IJFLD
                XVLX = SL(2,NBT)*PORD(2,NBT)
                FCLT = YL(NBT,NSL)/(XVLX+SMALL)
                CRLT = ABS( WL(1,NQZ) )*DT/DZGP(NQZ)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
     &            /(C(NBT,NSL)*FCLT-C(N,NSL)*FCLP+SMALL))
     &            *((DZGF(NBT)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRLT,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBT))
                WCLZ = C(N,NSL)*FLT*(1.D+0-THETA*DZF)*FCLP
     &            + C(NBT,NSL)*FLT*THETA*DZF*FCLT
                WCLZF = CO(N,NSL)*FLT*FCLP
                WC(NQZ,NSL) = WC(NQZ,NSL) + (WCLZ-WCLZF)/AFZ(NQZ)
              ENDIF
            ELSE
              ALB = MAX( FLB,ZERO )
              AP = (ALB-FLB)*FCLP
              AB = ALB*FCL
              WC(NPZ,NSL) = WC(NPZ,NSL) + (BCX(1)*AB - C(N,NSL)*AP)
            ENDIF
!
!---      Inflow gas ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (WG(1,NPZ)/EPSL.GT.EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 )  THEN
              WCGZ = 0.D+0
              IF( FGB.GE.ZERO ) WCGZ = BCX(1)*FCG*FGB
              WC(NPZ,NSL) = WC(NPZ,NSL) + WCGZ/AFZ(NPZ)
!
!---          TVD Transport for interior surface 
!             adjacent to boundary  ---
!
              FGT = AFZ(NQZ)*WG(1,NQZ)
              IF( FGT.GE.ZERO ) THEN
                NBT = N+IJFLD
                XVGX = SG(2,NBT)*PORD(2,NBT)
                FCGT = YG(NBT,NSL)/(XVGX+SMALL)
                CRGT = ABS( WG(1,NQZ) )*DT/DZGP(NQZ)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
     &            /(C(NBT,NSL)*FCGT-C(N,NSL)*FCGP+SMALL))
     &            *((DZGF(NBT)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRGT,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBT))
                WCGZ = C(N,NSL)*FGT*(1.D+0-THETA*DZF)*FCGP
     &            + C(NBT,NSL)*FGT*THETA*DZF*FCGT
                WCGZF = CO(N,NSL)*FGT*FCGP
                WC(NQZ,NSL) = WC(NQZ,NSL) + (WCGZ-WCGZF)/AFZ(NQZ)
              ENDIF
            ELSE
              AGB = MAX( FGB,ZERO )
              AP = (AGB-FGB)*FCGP
              AB = AGB*FCG
              WC(NPZ,NSL) = WC(NPZ,NSL) + (BCX(1)*AB - C(N,NSL)*AP)
            ENDIF
          ENDIF
!
!---    South boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
          NPY = NSY(N)
          NQY = NPY+IFLD
          IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVSB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULSX,VLSX,WLSX,N,MF )
            CALL SHDP( VLSX,WLSX,ULSX,DISPL(IZN),DISPT(IZN),DPLS )
            CALL ADVSB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGSX,VGSX,WGSX,N,MF )
            CALL SHDP( VGSX,WGSX,UGSX,DISPL(IZN),DISPT(IZN),DPGS )
          ELSE
            DPLS = 0.D+0
            DPGS = 0.D+0
          ENDIF
          FLS = AFY(NPY)*VL(1,NPY)
          FGS = AFY(NPY)*VG(1,NPY)
          CRLS = ABS( VL(1,NPY) )*DT/(RP(I)*DYGF(N)*XVLB+SMALL)
          CRGS = ABS( VG(1,NPY) )*DT/(RP(I)*DYGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            INDX = 16
            DLY = DIFMN(DLB,DLP,DYGF(N),DYGF(N),VL(1,NPY),INDX)
            DLY = AFY(NPY)*(DLY+DPLB)/RP(I)/(5.D-1*DYGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGY = DIFMN(DGB,DGP,DYGF(N),DYGF(N),VG(1,NPY),INDX)
            DGY = AFY(NPY)*(DGY+DPGB)/RP(I)/(5.D-1*DYGF(N))
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 )  THEN
              IF( FLS.GE.ZERO ) THEN
                VCLY = BCX(1)*FCL*FLS
              ELSEIF( FLS.LT.ZERO .AND. J.LT.JFLD ) THEN
                NBN = N+IFLD
                FCLN = YL(NBN,NSL)/(SL(2,NBN)*PORD(2,NBN)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBN,NSL)*FCLN)
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
                THETA = FLIMIT( R,CRLS,ISLC(1) )
                VCLY = BCX(1)*FLS*THETA*FCL
     &            + C(N,NSL)*FLS*(1.D+0-THETA)*FCLP
              ELSEIF( FLS.LT.ZERO .AND. J.EQ.JFLD ) THEN
                 VCLY = C(N,NSL)*FLS*FCLP
             ENDIF
              IF( FGS.GE.ZERO ) THEN
                VCGY = BCX(1)*FCG*FGS
              ELSEIF( FGS.LT.ZERO .AND. J.LT.JFLD ) THEN
                NBN = N+IFLD
                FCGN = YG(NBN,NSL)/(SG(2,NBN)*PORD(2,NBN)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBN,NSL)*FCGN)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
                THETA = FLIMIT( R,CRGS,ISLC(1) )
                VCGY = BCX(1)*FGS*THETA*FCG
     &            + C(N,NSL)*FGS*(1.D+0-THETA)*FCGP
              ELSEIF( FGS.LT.ZERO .AND. J.EQ.JFLD ) THEN
                VCGY = C(N,NSL)*FGS*FCGP
              ENDIF
              AS = DLY*FCL + DGY*FCG
              AP = DGY*FCGP + DLY*FCLP
              VC(NPY,NSL) = VC(NPY,NSL) + (VCLY+VCGY)/AFY(NPY)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FLN = AFY(NQY)*VL(1,NQY)
              IF( FLN.GE.ZERO ) THEN
                NBN = N+IFLD
                XVLX = SL(2,NBN)*PORD(2,NBN)
                FCLN = YL(NBN,NSL)/(XVLX+SMALL)
                CRLN = ABS( VL(1,NQY) )*DT/DYGP(NQY)/(RP(I)*XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
     &            /(C(NBN,NSL)*FCLN-C(N,NSL)*FCLP+SMALL))
     &            *((DYGF(NBN)+DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRLN,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBN))
                VCLY = C(N,NSL)*FLN*(1.D+0-THETA*DYF)*FCLP
     &            + C(NBN,NSL)*FLN*THETA*DYF*FCLN
                VCLYF = CO(N,NSL)*FLN*FCLP
                VC(NQY,NSL) = VC(NQY,NSL) + (VCLY-VCLYF)/AFY(NQY)
              ENDIF
              FGN = AFY(NQY)*VG(1,NQY)
              IF( FGN.GE.ZERO ) THEN
                NBN = N+IFLD
                XVGX = SG(2,NBN)*PORD(2,NBN)
                FCGN = YG(NBN,NSL)/(XVGX+SMALL)
                CRGN = ABS( VG(1,NQY) )*DT/DYGP(NQY)/(RP(I)*XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
     &            /(C(NBN,NSL)*FCGN-C(N,NSL)*FCGP+SMALL))
     &            *((DYGF(NBN)+DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRGN,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBN))
                VCGY = C(N,NSL)*FGN*(1.D+0-THETA*DYF)*FCGP
     &            + C(NBN,NSL)*FGN*THETA*DYF*FCGN
                VCGYF = CO(N,NSL)*FGN*FCGP
                VC(NQY,NSL) = VC(NQY,NSL) + (VCGY-VCGYF)/AFY(NQY)
              ENDIF
            ELSE
              ALS = MAX( FLS,ZERO ) +
     &          DLY*MAX((ONE-(TENTH*ABS(FLS)/(DLY+SMALL)))**5,ZERO)
              AGS = MAX( FGS,ZERO ) +
     &          DGY*MAX((ONE-(TENTH*ABS(FGS)/(DGY+SMALL)))**5,ZERO)
              AP = (ALS-FLS)*FCLP + (AGS-FGS)*FCGP
              AS = ALS*FCL + AGS*FCG
              VC(NPY,NSL) = VC(NPY,NSL) + (BCX(1)*AS - C(N,NSL)*AP)
            ENDIF
!
!---      Outflow aqueous  ---
!
          ELSEIF( IBCTX.EQ.7 .OR.
     &      ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (VL(1,NPY)/EPSL.LT.-EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 )  THEN
              VCLY = 0.D+0
              IF( FLS.LT.ZERO .AND. J.LT.JFLD ) THEN
                NBN = N+IFLD
                FCLN = YL(NBN,NSL)/(SL(2,NBN)*PORD(2,NBN)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBN,NSL)*FCLN)
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
                THETA = FLIMIT( R,CRLS,ISLC(1) )
                VCLY = BCX(1)*FLS*THETA*FCL
     &            + C(N,NSL)*FLS*(1.D+0-THETA)*FCLP
              ELSEIF( FLS.LT.ZERO .AND. J.EQ.JFLD ) THEN
                 VCLY = C(N,NSL)*FLS*FCLP
              ENDIF
              VC(NPY,NSL) = VC(NPY,NSL) + VCLY/AFY(NPY)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FLN = AFY(NQY)*VL(1,NQY)
              IF( FLN.GE.ZERO ) THEN
                NBN = N+IFLD
                XVLX = SL(2,NBN)*PORD(2,NBN)
                FCLN = YL(NBN,NSL)/(XVLX+SMALL)
                CRLN = ABS( VL(1,NQY) )*DT/DYGP(NQY)/(RP(I)*XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
     &            /(C(NBN,NSL)*FCLN-C(N,NSL)*FCLP+SMALL))
     &            *((DYGF(NBN)+DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRLN,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBN))
                VCLY = C(N,NSL)*FLN*(1.D+0-THETA*DYF)*FCLP
     &            + C(NBN,NSL)*FLN*THETA*DYF*FCLN
                VCLYF = CO(N,NSL)*FLN*FCLP
                VC(NQY,NSL) = VC(NQY,NSL) + (VCLY-VCLYF)/AFY(NQY)
              ENDIF
            ELSE
              ALS = MAX( FLS,ZERO )
              AP = (ALS-FLS)*FCLP
              AS = ALS*FCL
            ENDIF
            VC(NPY,NSL) = VC(NPY,NSL) + (BCX(1)*AS - C(N,NSL)*AP)
!
!---      Outflow gas  ---
!
          ELSEIF( IBCTX.EQ.7 .OR.
     &      ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (VG(1,NPY)/EPSL.LT.-EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 )  THEN
              VCGY = 0.D+0
              IF( FGS.LT.ZERO .AND. J.LT.JFLD ) THEN
                NBN = N+IFLD
                FCGN = YG(NBN,NSL)/(SG(2,NBN)*PORD(2,NBN)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBN,NSL)*FCGN)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBN)))
                THETA = FLIMIT( R,CRGS,ISLC(1) )
                VCGY = BCX(1)*FGS*THETA*FCG
     &            + C(N,NSL)*FGS*(1.D+0-THETA)*FCGP
              ELSEIF( FGS.LT.ZERO .AND. J.EQ.JFLD ) THEN
                VCGY = C(N,NSL)*FGS*FCGP
              ENDIF
              VC(NPY,NSL) = VC(NPY,NSL) + VCGY/AFY(NPY)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FGN = AFY(NQY)*VG(1,NQY)
              IF( FGN.GE.ZERO ) THEN
                NBN = N+IFLD
                XVGX = SG(2,NBN)*PORD(2,NBN)
                FCGN = YG(NBN,NSL)/(XVGX+SMALL)
                CRGN = ABS( VG(1,NQY) )*DT/DYGP(NQY)/(RP(I)*XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
     &            /(C(NBN,NSL)*FCGN-C(N,NSL)*FCGP+SMALL))
     &            *((DYGF(NBN)+DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRGN,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBN))
                VCGY = C(N,NSL)*FGN*(1.D+0-THETA*DYF)*FCGP
     &            + C(NBN,NSL)*FGN*THETA*DYF*FCGN
                VCGYF = CO(N,NSL)*FGN*FCGP
                VC(NQY,NSL) = VC(NQY,NSL) + (VCGY-VCGYF)/AFY(NQY)
              ENDIF
            ELSE
              AGS = MAX( FGS,ZERO )
              AP = (AGS-FGS)*FCGP
              AS = AGS*FCG
            ENDIF
            VC(NPY,NSL) = VC(NPY,NSL) + (BCX(1)*AS - C(N,NSL)*AP)
!
!---      Inflow aqueous ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (VL(1,NPY)/EPSL.GT.EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 )  THEN
              VCLY = 0.D+0
              IF( FLS.GE.ZERO ) VCLY = BCX(1)*FCL*FLS
              VC(NPY,NSL) = VC(NPY,NSL) + VCLY/AFY(NPY)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FLN = AFY(NQY)*VL(1,NQY)
              IF( FLN.GE.ZERO ) THEN
                NBN = N+IFLD
                XVLX = SL(2,NBN)*PORD(2,NBN)
                FCLN = YL(NBN,NSL)/(XVLX+SMALL)
                CRLN = ABS( VL(1,NQY) )*DT/DYGP(NQY)/(RP(I)*XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
     &            /(C(NBN,NSL)*FCLN-C(N,NSL)*FCLP+SMALL))
     &            *((DYGF(NBN)+DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRLN,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBN))
                VCLY = C(N,NSL)*FLN*(1.D+0-THETA*DYF)*FCLP
     &            + C(NBN,NSL)*FLN*THETA*DYF*FCLN
                VCLYF = CO(N,NSL)*FLN*FCLP
                VC(NQY,NSL) = VC(NQY,NSL) + (VCLY-VCLYF)/AFY(NQY)
              ENDIF
            ELSE
              ALS = MAX( FLS,ZERO )
              AP = (ALS-FLS)*FCLP
              AS = ALS*FCL
            ENDIF
            VC(NPY,NSL) = VC(NPY,NSL) + (BCX(1)*AS - C(N,NSL)*AP)
!
!---      Inflow gas ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (VG(1,NPY)/EPSL.GT.EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 )  THEN
              VCGY = 0.D+0
              IF( FGS.GE.ZERO ) VCGY = BCX(1)*FCG*FGS
              VC(NPY,NSL) = VC(NPY,NSL) + VCGY/AFY(NPY)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FGN = AFY(NQY)*VG(1,NQY)
              IF( FGN.GE.ZERO ) THEN
                NBN = N+IFLD
                XVGX = SG(2,NBN)*PORD(2,NBN)
                FCGN = YG(NBN,NSL)/(XVGX+SMALL)
                CRGN = ABS( VG(1,NQY) )*DT/DYGP(NQY)/(RP(I)*XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
     &            /(C(NBN,NSL)*FCGN-C(N,NSL)*FCGP+SMALL))
     &            *((DYGF(NBN)+DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRGN,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBN))
                VCGY = C(N,NSL)*FGN*(1.D+0-THETA*DYF)*FCGP
     &            + C(NBN,NSL)*FGN*THETA*DYF*FCGN
                VCGYF = CO(N,NSL)*FGN*FCGP
                VC(NQY,NSL) = VC(NQY,NSL) + (VCGY-VCGYF)/AFY(NQY)
              ENDIF
            ELSE
              AGS = MAX( FGS,ZERO )
              AP = (AGS-FGS)*FCGP
              AS = AGS*FCG
            ENDIF
            VC(NPY,NSL) = VC(NPY,NSL) + (BCX(1)*AS - C(N,NSL)*AP)
          ENDIF
!
!---    West boundary  ---
!
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
          NPX = NSX(N)
          NQX = NPX+1
          IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVWB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULX,VLX,WLX,N,MF )
            CALL SHDP( ULX,VLX,WLX,DISPL(IZN),DISPT(IZN),DPLW )
            CALL ADVWB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGX,VGX,WGX,N,MF )
            CALL SHDP( UGX,VGX,WGX,DISPL(IZN),DISPT(IZN),DPGW )
            ELSE
              DPLW = 0.D+0
              DPGW = 0.D+0
            ENDIF
            FLW = AFX(NPX)*UL(1,NPX)
            FGW = AFX(NPX)*UG(1,NPX)
            CRLW = ABS( UL(1,NPX) )*DT/(DXGF(N)*XVLB+SMALL)
            CRGW = ABS( UG(1,NPX) )*DT/(DXGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            INDX = 16
            DLX = DIFMN(DLB,DLP,DXGF(N),DXGF(N),UL(1,NPX),INDX)
            DLX = AFX(NPX)*(DLX+DPLW)/(5.D-1*DXGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGX = DIFMN(DGB,DGP,DXGF(N),DXGF(N),UG(1,NPX),INDX)
            DGX = AFX(NPX)*(DGX+DPGW)/(5.D-1*DXGF(N))
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 )  THEN
              IF( FLW.GE.ZERO ) THEN
                UCLX = BCX(1)*FCL*FLW
              ELSEIF( FLW.LT.ZERO .AND. I.LT.IFLD ) THEN
                NBE = N+1
                FCLE = YL(NBE,NSL)/(SL(2,NBE)*PORD(2,NBE)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBE,NSL)*FCLE)
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
                THETA = FLIMIT( R,CRLW,ISLC(1) )
                UCLX = C(N,NSL)*FLW*(1.D+0-THETA)*FCLP
     &            + BCX(1)*FLW*THETA*FCL
              ELSEIF( FLW.LT.ZERO .AND. I.EQ.IFLD ) THEN
                UCLX = C(N,NSL)*FLW*FCLP
              ENDIF
              IF( FGW.GE.ZERO ) THEN
                UCGX = BCX(1)*FCG*FGW
              ELSEIF( FGW.LT.ZERO .AND. I.LT.IFLD ) THEN
                NBE = N+1
                FCGE = YG(NBE,NSL)/(SG(2,NBE)*PORD(2,NBE)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBE,NSL)*FCGE)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
                THETA = FLIMIT( R,CRGW,ISLC(1) )
                UCGX = C(N,NSL)*FGW*(1.D+0-THETA)*FCGP
     &            + BCX(1)*FGW*THETA*FCG
              ELSEIF( FGW.LT.ZERO .AND. I.EQ.IFLD ) THEN
                UCGX = C(N,NSL)*FGW*FCGP
              ENDIF
              AW = DLX*FCL + DGX*FCG
              AP = DLX*FCLP + DGX*FCGP
              UC(NPX,NSL) = UC(NPX,NSL) + (UCLX+UCGX)/AFX(NPX)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FLE = AFX(NQX)*UL(1,NQX)
              IF( FLE.GE.ZERO ) THEN
                NBE = N+1
                XVLX = SL(2,NBE)*PORD(2,NBE)
                FCLE = YL(NBE,NSL)/(XVLX+SMALL)
                CRLE = ABS( UL(1,NQX) )*DT/DXGP(NQX)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
     &            /(C(NBE,NSL)*FCLE-C(N,NSL)*FCLP+SMALL))
     &            *((DXGF(NBE)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRLE,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBE))
                UCLX = C(N,NSL)*FLE*(1.D+0-THETA*DXF)*FCLP
     &            + C(NBE,NSL)*FLE*THETA*DXF*FCLE
                UCLXF = CO(N,NSL)*FLE*FCLP
                UC(NQX,NSL) = UC(NQX,NSL) + (UCLX-UCLXF)/AFX(NQX)
              ENDIF
              FGE = AFX(NQX)*UG(1,NQX)
              IF( FGE.GE.ZERO ) THEN
                NBE = N+1
                XVGX = SG(2,NBE)*PORD(2,NBE)
                FCGE = YG(NBE,NSL)/(XVGX+SMALL)
                CRGE = ABS( UG(1,NQX) )*DT/DXGP(NQX)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
     &            /(C(NBE,NSL)*FCGE-C(N,NSL)*FCGP+SMALL))
     &            *((DXGF(NBE)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRGE,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBE))
                UCGX = C(N,NSL)*FGE*(1.D+0-THETA*DXF)*FCGP
     &            + C(NBE,NSL)*FGE*THETA*DXF*FCGE
                UCGXF = CO(N,NSL)*FGE*FCGP
                UC(NQX,NSL) = UC(NQX,NSL) + (UCGX-UCGXF)/AFX(NQX)
              ENDIF
            ELSE
              ALW = MAX(FLW,ZERO)
     &          + DLX*MAX((ONE-(TENTH*ABS(FLW)/(DLX+SMALL)))**5,ZERO)
              AGW = MAX(FGW,ZERO)
     &          + DGX*MAX((ONE-(TENTH*ABS(FGW)/(DGX+SMALL)))**5,ZERO)
              AP = (ALW-FLW)*FCLP + (AGW-FGW)*FCGP
              AW = ALW*FCL + AGW*FCG
              UC(NPX,NSL) = UC(NPX,NSL) + (BCX(1)*AW - C(N,NSL)*AP)
            ENDIF
!
!---      Outflow aqueous  ---
!
          ELSEIF( IBCTX.EQ.7 .OR.
     &      ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (UL(1,NPX)/EPSL.LT.-EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 )  THEN
              UCLX = 0.D+0
              IF( FLW.LT.ZERO .AND. I.LT.IFLD ) THEN
                NBE = N+1
                FCLE = YL(NBE,NSL)/(SL(2,NBE)*PORD(2,NBE)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBE,NSL)*FCLE)
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
                THETA = FLIMIT( R,CRLW,ISLC(1) )
                UCLX = C(N,NSL)*FLW*(1.D+0-THETA)*FCLP
     &            + BCX(1)*FLW*THETA*FCL
              ELSEIF( FLW.LT.ZERO .AND. I.EQ.IFLD ) THEN
                UCLX = C(N,NSL)*FLW*FCLP
              ENDIF
              UC(NPX,NSL) = UC(NPX,NSL) + UCLX/AFX(NPX)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FLE = AFX(NQX)*UL(1,NQX)
              IF( FLE.GE.ZERO ) THEN
                NBE = N+1
                XVLX = SL(2,NBE)*PORD(2,NBE)
                FCLE = YL(NBE,NSL)/(XVLX+SMALL)
                CRLE = ABS( UL(1,NQX) )*DT/DXGP(NQX)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
     &            /(C(NBE,NSL)*FCLE-C(N,NSL)*FCLP+SMALL))
     &            *((DXGF(NBE)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRLE,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBE))
                UCLX = C(N,NSL)*FLE*(1.D+0-THETA*DXF)*FCLP
     &            + C(NBE,NSL)*FLE*THETA*DXF*FCLE
                UCLXF = CO(N,NSL)*FLE*FCLP
                UC(NQX,NSL) = UC(NQX,NSL) + (UCLX-UCLXF)/AFX(NQX)
              ENDIF
            ELSE
              ALW = MAX(FLW,ZERO)
              AP = (ALW-FLW)*FCLP
              AW = ALW*FCL
              UC(NPX,NSL) = UC(NPX,NSL) + (BCX(1)*AW - C(N,NSL)*AP)
            ENDIF
!
!---      Outflow gas  ---
!
          ELSEIF( IBCTX.EQ.7 .OR.
     &      ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (UG(1,NPX)/EPSL.LT.-EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 )  THEN
              UCGX = 0.D+0
              IF( FGW.LT.ZERO .AND. I.LT.IFLD ) THEN
                NBE = N+1
                FCGE = YG(NBE,NSL)/(SG(2,NBE)*PORD(2,NBE)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBE,NSL)*FCGE)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBE)))
                THETA = FLIMIT( R,CRGW,ISLC(1) )
                UCGX = C(N,NSL)*FGW*(1.D+0-THETA)*FCGP
     &            + BCX(1)*FGW*THETA*FCG
              ELSEIF( FGW.LT.ZERO .AND. I.EQ.IFLD ) THEN
                UCGX = C(N,NSL)*FGW*FCGP
              ENDIF
              UC(NPX,NSL) = UC(NPX,NSL) + UCGX/AFX(NPX)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FGE = AFX(NQX)*UG(1,NQX)
              IF( FGE.GE.ZERO ) THEN
                NBE = N+1
                XVGX = SG(2,NBE)*PORD(2,NBE)
                FCGE = YG(NBE,NSL)/(XVGX+SMALL)
                CRGE = ABS( UG(1,NQX) )*DT/DXGP(NQX)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
     &            /(C(NBE,NSL)*FCGE-C(N,NSL)*FCGP+SMALL))
     &            *((DXGF(NBE)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRGE,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBE))
                UCGX = C(N,NSL)*FGE*(1.D+0-THETA*DXF)*FCGP
     &            + C(NBE,NSL)*FGE*THETA*DXF*FCGE
                UCGXF = CO(N,NSL)*FGE*FCGP
                UC(NQX,NSL) = UC(NQX,NSL) + (UCGX-UCGXF)/AFX(NQX)
              ENDIF
            ELSE
              AGW = MAX(FGW,ZERO)
              AP = (AGW-FGW)*FCGP
              AW = AGW*FCG
              UC(NPX,NSL) = UC(NPX,NSL) + (BCX(1)*AW - C(N,NSL)*AP)
            ENDIF
!
!---      Inflow aqueous ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (UL(1,NPX)/EPSL.GT.EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 )  THEN
              UCLX = 0.D+0
              IF( FLW.GE.ZERO ) UCLX = BCX(1)*FCL*FLW
              UC(NPX,NSL) = UC(NPX,NSL) + UCLX/AFX(NPX)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FLE = AFX(NQX)*UL(1,NQX)
              IF( FLE.GE.ZERO ) THEN
                NBE = N+1
                XVLX = SL(2,NBE)*PORD(2,NBE)
                FCLE = YL(NBE,NSL)/(XVLX+SMALL)
                CRLE = ABS( UL(1,NQX) )*DT/DXGP(NQX)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
     &            /(C(NBE,NSL)*FCLE-C(N,NSL)*FCLP+SMALL))
     &            *((DXGF(NBE)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRLE,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBE))
                UCLX = C(N,NSL)*FLE*(1.D+0-THETA*DXF)*FCLP
     &            + C(NBE,NSL)*FLE*THETA*DXF*FCLE
                UCLXF = CO(N,NSL)*FLE*FCLP
                UC(NQX,NSL) = UC(NQX,NSL) + (UCLX-UCLXF)/AFX(NQX)
              ENDIF
            ELSE
              ALW = MAX(FLW,ZERO)
              AP = (ALW-FLW)*FCLP
              AW = ALW*FCL
              UC(NPX,NSL) = UC(NPX,NSL) + (BCX(1)*AW - C(N,NSL)*AP)
            ENDIF
!
!---      Inflow gas ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (UG(1,NPX)/EPSL.GT.EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 )  THEN
              UCGX = 0.D+0
              IF( FGW.GE.ZERO ) UCGX = BCX(1)*FCG*FGW
              UC(NPX,NSL) = UC(NPX,NSL) + UCGX/AFX(NPX)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FGE = AFX(NQX)*UG(1,NQX)
              IF( FGE.GE.ZERO ) THEN
                NBE = N+1
                XVGX = SG(2,NBE)*PORD(2,NBE)
                FCGE = YG(NBE,NSL)/(XVGX+SMALL)
                CRGE = ABS( UG(1,NQX) )*DT/DXGP(NQX)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
     &            /(C(NBE,NSL)*FCGE-C(N,NSL)*FCGP+SMALL))
     &            *((DXGF(NBE)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRGE,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBE))
                UCGX = C(N,NSL)*FGE*(1.D+0-THETA*DXF)*FCGP
     &            + C(NBE,NSL)*FGE*THETA*DXF*FCGE
                UCGXF = CO(N,NSL)*FGE*FCGP
                UC(NQX,NSL) = UC(NQX,NSL) + (UCGX-UCGXF)/AFX(NQX)
              ENDIF
            ELSE
              AGW = MAX(FGW,ZERO)
              AP = (AGW-FGW)*FCGP
              AW = AGW*FCG
              UC(NPX,NSL) = UC(NPX,NSL) + (BCX(1)*AW - C(N,NSL)*AP)
            ENDIF
          ENDIF
!
!---    East boundary
!
        ELSEIF( IBCD(NB).EQ.1 ) THEN
          NQX = NSX(N) + 1
          IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
          NPX = NSX(N)
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVEB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULEX,VLEX,WLEX,N,MF )
            CALL SHDP( ULEX,VLEX,WLEX,DISPL(IZN),DISPT(IZN),DPLE )
            CALL ADVEB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGEX,VGEX,WGEX,N,MF )
            CALL SHDP( UGEX,VGEX,WGEX,DISPL(IZN),DISPT(IZN),DPGE )
          ELSE
            DPLE = 0.D+0
            DPGE = 0.D+0
          ENDIF
          FLE = AFX(NQX)*UL(1,NQX)
          FGE = AFX(NQX)*UG(1,NQX)
          CRLE = ABS( UL(1,NQX) )*DT/(DXGF(N)*XVLB+SMALL)
          CRGE = ABS( UG(1,NQX) )*DT/(DXGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            INDX = 16
            DLX = DIFMN(DLP,DLB,DXGF(N),DXGF(N),UL(1,NQX),INDX)
            DLX = AFX(NQX)*(DLX+DPLE)/(5.D-1*DXGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGX = DIFMN(DGP,DGB,DXGF(N),DXGF(N),UG(1,NQX),INDX)
            DGX = AFX(NQX)*(DGX+DPGE)/(5.D-1*DXGF(N))
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 ) THEN
              IF( FLE.LT.ZERO ) THEN
                UCLX = BCX(1)*FCL*FLE
              ELSEIF( FLE.GE.ZERO .AND. I.GT.1 ) THEN
                NBW = N-1
                FCLW = YL(NBW,NSL)/(SL(2,NBW)*PORD(2,NBW)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBW,NSL)*FCLW)
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
                THETA = FLIMIT( R,CRLE,ISLC(1) )
                UCLX =  C(N,NSL)*FLE*(1.D+0-THETA)*FCLP
     &          +  BCX(1)*FLE*THETA*FCL
              ELSEIF( FLE.GE.ZERO .AND. I.EQ.1 ) THEN
                UCLX =  C(N,NSL)*FLE*FCLP
              ENDIF
              IF( FGE.LT.ZERO ) THEN
                UCGX = BCX(1)*FCG*FGE
              ELSEIF( FGE.GE.ZERO .AND. I.GT.1 ) THEN
                NBW = N-1
                FCGW = YG(NBW,NSL)/(SG(2,NBW)*PORD(2,NBW)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBW,NSL)*FCGW)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
                THETA = FLIMIT( R,CRGE,ISLC(1) )
                UCGX =  C(N,NSL)*FGE*(1.D+0-THETA)*FCGP
     &          +  BCX(1)*FGE*THETA*FCG
              ELSEIF( FGE.GE.ZERO .AND. I.EQ.1 ) THEN
                UCGX =  C(N,NSL)*FGE*FCGP
              ENDIF
              AE = DLX*FCL + DGX*FCG
              AP = DLX*FCLP + DGX*FCGP
              UC(NQX,NSL) = UC(NQX,NSL) + (UCLX+UCGX)/AFX(NQX)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FLW = AFX(NPX)*UL(1,NPX)
              IF( FLW.LT.ZERO ) THEN
                NBW = N-1
                XVLX = SL(2,NBW)*PORD(2,NBW)
                CRLW = ABS( UL(1,NPX) )*DT/DXGP(NPX)/(XVLX+SMALL)
                FCLW = YL(NBW,NSL)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
     &            /(C(NBW,NSL)*FCLW-C(N,NSL)*FCLP+SMALL))
     &            *((DXGF(NBW)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRLW,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBW))
                UCLX = C(N,NSL)*FLW*(1.D+0-THETA*DXF)*FCLP
     &            + C(NBW,NSL)*FLW*THETA*DXF*FCLW
                UCLXF = CO(N,NSL)*FLW*FCLP
                UC(NPX,NSL) = UC(NPX,NSL) + (UCLX-UCLXF)/AFX(NPX)
              ENDIF
              IF( FGW.LT.ZERO ) THEN
                NBW = N-1
                XVGX = SG(2,NBW)*PORD(2,NBW)
                CRGW = ABS( UG(1,NPX) )*DT/DXGP(NPX)/(XVGX+SMALL)
                FCGW = YG(NBW,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
     &            /(C(NBW,NSL)*FCGW-C(N,NSL)*FCGP+SMALL))
     &            *((DXGF(NBW)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRGW,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBW))
                UCGX = C(N,NSL)*FGW*(1.D+0-THETA*DXF)*FCGP
     &            + C(NBW,NSL)*FGW*THETA*DXF*FCGW
                UCGXF = CO(N,NSL)*FGW*FCGP
                UC(NPX,NSL) = UC(NPX,NSL) + (UCGX-UCGXF)/AFX(NPX)
              ENDIF
            ELSE
              ALE = MAX( -FLE,ZERO ) +
     &          DLX*MAX((ONE-(TENTH*ABS(FLE)/(DLX+SMALL)))**5,ZERO)
              AGE = MAX( -FGE,ZERO ) +
     &          DGX*MAX((ONE-(TENTH*ABS(FGE)/(DGX+SMALL)))**5,ZERO)
              AP = (ALE+FLE)*FCLP + (AGE+FGE)*FCGP
              AE = ALE*FCL + AGE*FCG
              UC(NQX,NSL) = UC(NQX,NSL)+(C(N,NSL)*AP-BCX(1)*AE)
            ENDIF
!
!---      Outflow aqueous  ---
!
          ELSEIF( IBCTX.EQ.7 .OR.
     &      ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (UL(1,NQX)/EPSL.GT.EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 ) THEN
              UCLX = 0.D+0
              IF( FLE.GE.ZERO .AND. I.GT.1 ) THEN
                NBW = N-1
                FCLW = YL(NBW,NSL)/(SL(2,NBW)*PORD(2,NBW)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBW,NSL)*FCLW)
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
                THETA = FLIMIT( R,CRLE,ISLC(1) )
                UCLX =  C(N,NSL)*FLE*(1.D+0-THETA)*FCLP
     &          +  BCX(1)*FLE*THETA*FCL
              ELSEIF( FLE.GE.ZERO .AND. I.EQ.1 ) THEN
                UCLX =  C(N,NSL)*FLE*FCLP
              ENDIF
              UC(NQX,NSL) = UC(NQX,NSL) + UCLX/AFX(NQX)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FLW = AFX(NPX)*UL(1,NPX)
              IF( FLW.LT.ZERO ) THEN
                NBW = N-1
                XVLX = SL(2,NBW)*PORD(2,NBW)
                CRLW = ABS( UL(1,NPX) )*DT/DXGP(NPX)/(XVLX+SMALL)
                FCLW = YL(NBW,NSL)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
     &            /(C(NBW,NSL)*FCLW-C(N,NSL)*FCLP+SMALL))
     &            *((DXGF(NBW)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRLW,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBW))
                UCLX = C(N,NSL)*FLW*(1.D+0-THETA*DXF)*FCLP
     &            + C(NBW,NSL)*FLW*THETA*DXF*FCLW
                UCLXF = CO(N,NSL)*FLW*FCLP
                UC(NPX,NSL) = UC(NPX,NSL) + (UCLX-UCLXF)/AFX(NPX)
              ENDIF
            ELSE
              ALE = MAX( -FLE,ZERO )
              AP = (ALE+FLE)*FCLP
              AE = ALE*FCL
              UC(NQX,NSL) = UC(NQX,NSL)+(C(N,NSL)*AP-BCX(1)*AE)
            ENDIF
!
!---      Outflow gas  ---
!
          ELSEIF( IBCTX.EQ.7 .OR.
     &      ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (UG(1,NQX)/EPSL.GT.EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 ) THEN
              UCGX = 0.D+0
              IF( FGE.GE.ZERO .AND. I.GT.1 ) THEN
                NBW = N-1
                FCGW = YG(NBW,NSL)/(SG(2,NBW)*PORD(2,NBW)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBW,NSL)*FCGW)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DXGF(N)/(DXGF(N)+DXGF(NBW)))
                THETA = FLIMIT( R,CRGE,ISLC(1) )
                UCGX =  C(N,NSL)*FGE*(1.D+0-THETA)*FCGP
     &          +  BCX(1)*FGE*THETA*FCG
              ELSEIF( FGE.GE.ZERO .AND. I.EQ.1 ) THEN
                UCGX =  C(N,NSL)*FGE*FCGP
              ENDIF
              UC(NQX,NSL) = UC(NQX,NSL) + UCGX/AFX(NQX)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              IF( FGW.LT.ZERO ) THEN
                NBW = N-1
                XVGX = SG(2,NBW)*PORD(2,NBW)
                CRGW = ABS( UG(1,NPX) )*DT/DXGP(NPX)/(XVGX+SMALL)
                FCGW = YG(NBW,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
     &            /(C(NBW,NSL)*FCGW-C(N,NSL)*FCGP+SMALL))
     &            *((DXGF(NBW)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRGW,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBW))
                UCGX = C(N,NSL)*FGW*(1.D+0-THETA*DXF)*FCGP
     &            + C(NBW,NSL)*FGW*THETA*DXF*FCGW
                UCGXF = CO(N,NSL)*FGW*FCGP
                UC(NPX,NSL) = UC(NPX,NSL) + (UCGX-UCGXF)/AFX(NPX)
              ENDIF
            ELSE
              AGE = MAX( -FGE,ZERO )
              AP = (AGE+FGE)*FCGP
              AE = AGE*FCG
              UC(NQX,NSL) = UC(NQX,NSL)+(C(N,NSL)*AP-BCX(1)*AE)
            ENDIF
!
!---      Inflow aqueous ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (UL(1,NQX)/EPSL.LT.-EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 ) THEN
              UCLX = 0.D+0
              IF( FLE.LT.ZERO ) UCLX = BCX(1)*FCL*FLE
              UC(NQX,NSL) = UC(NQX,NSL) + UCLX/AFX(NQX)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              IF( FGW.LT.ZERO ) THEN
                NBW = N-1
                XVGX = SG(2,NBW)*PORD(2,NBW)
                CRGW = ABS( UG(1,NPX) )*DT/DXGP(NPX)/(XVGX+SMALL)
                FCGW = YG(NBW,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
     &            /(C(NBW,NSL)*FCGW-C(N,NSL)*FCGP+SMALL))
     &            *((DXGF(NBW)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRGW,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBW))
                UCGX = C(N,NSL)*FGW*(1.D+0-THETA*DXF)*FCGP
     &            + C(NBW,NSL)*FGW*THETA*DXF*FCGW
                UCGXF = CO(N,NSL)*FGW*FCGP
                UC(NPX,NSL) = UC(NPX,NSL) + (UCGX-UCGXF)/AFX(NPX)
              ENDIF
            ELSE
              ALE = MAX( -FLE,ZERO )
              AP = (ALE+FLE)*FCLP
              AE = ALE*FCL
              UC(NQX,NSL) = UC(NQX,NSL)+(C(N,NSL)*AP-BCX(1)*AE)
            ENDIF
!
!---      Inflow gas ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (UG(1,NQX)/EPSL.LT.-EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. IFLD.GT.1 ) THEN
              UCGX = 0.D+0
              IF( FGE.LT.ZERO ) UCGX = BCX(1)*FCG*FGE
              UC(NQX,NSL) = UC(NQX,NSL) + UCGX/AFX(NQX)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              IF( FGW.LT.ZERO ) THEN
                NBW = N-1
                XVGX = SG(2,NBW)*PORD(2,NBW)
                CRGW = ABS( UG(1,NPX) )*DT/DXGP(NPX)/(XVGX+SMALL)
                FCGW = YG(NBW,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
     &            /(C(NBW,NSL)*FCGW-C(N,NSL)*FCGP+SMALL))
     &            *((DXGF(NBW)+DXGF(N))/DXGF(N))
                THETA = FLIMIT( R,CRGW,ISLC(1) )
                DXF = DXGF(N)/(DXGF(N)+DXGF(NBW))
                UCGX = C(N,NSL)*FGW*(1.D+0-THETA*DXF)*FCGP
     &            + C(NBW,NSL)*FGW*THETA*DXF*FCGW
                UCGXF = CO(N,NSL)*FGW*FCGP
                UC(NPX,NSL) = UC(NPX,NSL) + (UCGX-UCGXF)/AFX(NPX)
              ENDIF
            ELSE
              AGE = MAX( -FGE,ZERO )
              AP = (AGE+FGE)*FCGP
              AE = AGE*FCG
              UC(NQX,NSL) = UC(NQX,NSL)+(C(N,NSL)*AP-BCX(1)*AE)
            ENDIF
          ENDIF
!
!---    North boundary  ---
!
        ELSEIF( IBCD(NB).EQ.2 ) THEN
          NQY = NSY(N) + IFLD
          IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
          NPY = NSY(N)
!
!---      Hydraulic dispersion  ---
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVNB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULNX,VLNX,WLNX,N,MF )
            CALL SHDP( VLNX,WLNX,ULNX,DISPL(IZN),DISPT(IZN),DPLN )
            CALL ADVNB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGNX,VGNX,WGNX,N,MF )
            CALL SHDP( VGNX,WGNX,UGNX,DISPL(IZN),DISPT(IZN),DPGN )
          ELSE
            DPLN = 0.D+0
            DPGN = 0.D+0
          ENDIF
          FLN = AFY(NQY)*VL(1,NQY)
          FGN = AFY(NQY)*VG(1,NQY)
          CRLN = ABS( VL(1,NQY) )*DT/(RP(I)*DYGF(N)*XVLB+SMALL)
          CRGN = ABS( VG(1,NQY) )*DT/(RP(I)*DYGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            INDX = 16
            DLY = DIFMN(DLP,DLB,DYGF(N),DYGF(N),VL(1,NQY),INDX)
            DLY = AFY(NQY)*(DLY+DPLN)/RP(I)/(5.D-1*DYGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGY = DIFMN(DGP,DGB,DYGF(N),DYGF(N),VG(1,NQY),INDX)
            DGY = AFY(NQY)*(DGY+DPGN)/RP(I)/(5.D-1*DYGF(N))
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 ) THEN
              IF( FLN.LT.ZERO ) THEN
                VCLY = BCX(1)*FCL*FLN
              ELSEIF( FLN.GE.ZERO .AND. J.GT.1 ) THEN
                NBS = N-IFLD
                FCLS = YL(NBS,NSL)/(SL(2,NBS)*PORD(2,NBS)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBS,NSL)*FCLS)
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
                THETA = FLIMIT( R,CRLN,ISLC(1) )
                VCLY =  BCX(1)*FLN*THETA*FCL
     &            + C(N,NSL)*FLN*(1.D+0-THETA)*FCLP
              ELSEIF( FLN.GE.ZERO .AND. J.EQ.1 ) THEN
                VCLY =  C(N,NSL)*FLN*FCLP
              ENDIF
              IF( FGN.LT.ZERO ) THEN
                VCGY = BCX(1)*FCG*FGN
              ELSEIF( FGN.GE.ZERO .AND. J.GT.1 ) THEN
                NBS = N-IFLD
                FCGS = YG(NBS,NSL)/(SG(2,NBS)*PORD(2,NBS)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBS,NSL)*FCGS)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
                THETA = FLIMIT( R,CRGN,ISLC(1) )
                VCGY =  BCX(1)*FGN*THETA*FCG
     &            + C(N,NSL)*FGN*(1.D+0-THETA)*FCGP
              ELSEIF( FGN.GE.ZERO .AND. J.EQ.1 ) THEN
                VCGY =  C(N,NSL)*FGN*(1.D+0-THETA)*FCGP
              ENDIF
              AN = DLY*FCL + DGY*FCG
              AP = DLY*FCLP + DGY*FCGP
              VC(NQY,NSL) = VC(NQY,NSL) + (VCLY+VCGY)/AFY(NQY)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FLS = AFY(NPY)*VL(1,NPY)
              IF( FLS.LT.ZERO ) THEN
                NBS = N-IFLD
                XVLX = SL(2,NBS)*PORD(2,NBS)
                CRLS = ABS( VL(1,NPY) )*DT/DYGP(NPY)/(XVLX*RP(I)+SMALL)
                FCLS = YL(NBS,NSL)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
     &            /(C(NBS,NSL)*FCLS-C(N,NSL)*FCLP+SMALL))
     &            *((DYGF(NBS)-DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRLS,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBS))
                VCLY = C(N,NSL)*FLS*(1.D+0-THETA*DYF)*FCLP
     &            + C(NBS,NSL)*FLS*THETA*DYF*FCLS
                VCLYF = CO(N,NSL)*FLS*FCLP
                VC(NPY,NSL) = VC(NPY,NSL) + (VCLY-VCLYF)/AFY(NPY)
              ENDIF
              FGS = AFY(NPY)*VG(1,NPY)
              IF( FGS.LT.ZERO ) THEN
                NBS = N-IFLD
                XVGX = SG(2,NBS)*PORD(2,NBS)
                CRGS = ABS( VG(1,NPY) )*DT/DYGP(NPY)/(XVGX*RP(I)+SMALL)
                FCGS = YG(NBS,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
     &            /(C(NBS,NSL)*FCGS-C(N,NSL)*FCGP+SMALL))
     &            *((DYGF(NBS)-DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRGS,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBS))
                VCGY = C(N,NSL)*FGS*(1.D+0-THETA*DYF)*FCGP
     &            + C(NBS,NSL)*FGS*THETA*DYF*FCGS
                VCGYF = CO(N,NSL)*FGS*FCGP
                VC(NPY,NSL) = VC(NPY,NSL) + (VCGY-VCGYF)/AFY(NPY)
              ENDIF
            ELSE
              ALN = MAX( -FLN,ZERO ) +
     &          DLY*MAX((ONE-(TENTH*ABS(FLN)/(DLY+SMALL)))**5,ZERO)
              AGN = MAX( -FGN,ZERO ) +
     &          DGY*MAX((ONE-(TENTH*ABS(FGN)/(DGY+SMALL)))**5,ZERO)
              AP = (ALN+FLN)*FCLP + (AGN+FGN)*FCGP
              AN = ALN*FCL + AGN*FCG
              VC(NQY,NSL) = VC(NQY,NSL) + (C(N,NSL)*AP - BCX(1)*AN)
            ENDIF
!
!---      Outflow aqueous  ---
!
          ELSEIF( IBCTX.EQ.7 .OR.
     &      ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (VL(1,NQY)/EPSL.GT.EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 ) THEN
              VCLY = 0.D+0
              IF( FLN.GE.ZERO .AND. J.GT.1 ) THEN
                NBS = N-IFLD
                FCLS = YL(NBS,NSL)/(SL(2,NBS)*PORD(2,NBS)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBS,NSL)*FCLS)
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
                THETA = FLIMIT( R,CRLN,ISLC(1) )
                VCLY =  BCX(1)*FLN*THETA*FCL
     &            + C(N,NSL)*FLN*(1.D+0-THETA)*FCLP
              ELSEIF( FLN.GE.ZERO .AND. J.EQ.1 ) THEN
                VCLY =  C(N,NSL)*FLN*FCLP
              ENDIF
              VC(NQY,NSL) = VC(NQY,NSL) + VCLY/AFY(NQY)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FLS = AFY(NPY)*VL(1,NPY)
              IF( FLS.LT.ZERO ) THEN
                NBS = N-IFLD
                XVLX = SL(2,NBS)*PORD(2,NBS)
                CRLS = ABS( VL(1,NPY) )*DT/DYGP(NPY)/(XVLX*RP(I)+SMALL)
                FCLS = YL(NBS,NSL)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
     &            /(C(NBS,NSL)*FCLS-C(N,NSL)*FCLP+SMALL))
     &            *((DYGF(NBS)-DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRLS,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBS))
                VCLY = C(N,NSL)*FLS*(1.D+0-THETA*DYF)*FCLP
     &            + C(NBS,NSL)*FLS*THETA*DYF*FCLS
                VCLYF = CO(N,NSL)*FLS*FCLP
                VC(NPY,NSL) = VC(NPY,NSL) + (VCLY-VCLYF)/AFY(NPY)
              ENDIF
            ELSE
              ALN = MAX( -FLN,ZERO )
              AP = (ALN+FLN)*FCLP
              AN = ALN*FCL
              VC(NQY,NSL) = VC(NQY,NSL) + (C(N,NSL)*AP - BCX(1)*AN)
            ENDIF
!
!---      Outflow gas  ---
!
          ELSEIF( IBCTX.EQ.7 .OR.
     &      ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (VG(1,NQY)/EPSL.GT.EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 ) THEN
              VCGY = 0.D+0
              IF( FGN.GE.ZERO .AND. J.GT.1 ) THEN
                NBS = N-IFLD
                FCGS = YG(NBS,NSL)/(SG(2,NBS)*PORD(2,NBS)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBS,NSL)*FCGS)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DYGF(N)/(DYGF(N)+DYGF(NBS)))
                THETA = FLIMIT( R,CRGN,ISLC(1) )
                VCGY =  BCX(1)*FGN*THETA*FCG
     &            + C(N,NSL)*FGN*(1.D+0-THETA)*FCGP
              ELSEIF( FGN.GE.ZERO .AND. J.EQ.1 ) THEN
                VCGY =  C(N,NSL)*FGN*(1.D+0-THETA)*FCGP
              ENDIF
              VC(NQY,NSL) = VC(NQY,NSL) + VCGY/AFY(NQY)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FGS = AFY(NPY)*VG(1,NPY)
              IF( FGS.LT.ZERO ) THEN
                NBS = N-IFLD
                XVGX = SG(2,NBS)*PORD(2,NBS)
                CRGS = ABS( VG(1,NPY) )*DT/DYGP(NPY)/(XVGX*RP(I)+SMALL)
                FCGS = YG(NBS,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
     &            /(C(NBS,NSL)*FCGS-C(N,NSL)*FCGP+SMALL))
     &            *((DYGF(NBS)-DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRGS,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBS))
                VCGY = C(N,NSL)*FGS*(1.D+0-THETA*DYF)*FCGP
     &            + C(NBS,NSL)*FGS*THETA*DYF*FCGS
                VCGYF = CO(N,NSL)*FGS*FCGP
                VC(NPY,NSL) = VC(NPY,NSL) + (VCGY-VCGYF)/AFY(NPY)
              ENDIF
            ELSE
              AGN = MAX( -FGN,ZERO )
              AP = (AGN+FGN)*FCGP
              AN = AGN*FCG
              VC(NQY,NSL) = VC(NQY,NSL) + (C(N,NSL)*AP - BCX(1)*AN)
            ENDIF
!
!---      Inflow aqueous ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (VL(1,NQY)/EPSL.LT.-EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 ) THEN
              VCLY = 0.D+0
              IF( FLN.LT.ZERO ) VCLY = BCX(1)*FCL*FLN
              VC(NQY,NSL) = VC(NQY,NSL) + VCLY/AFY(NQY)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FLS = AFY(NPY)*VL(1,NPY)
              IF( FLS.LT.ZERO ) THEN
                NBS = N-IFLD
                XVLX = SL(2,NBS)*PORD(2,NBS)
                CRLS = ABS( VL(1,NPY) )*DT/DYGP(NPY)/(XVLX*RP(I)+SMALL)
                FCLS = YL(NBS,NSL)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
     &            /(C(NBS,NSL)*FCLS-C(N,NSL)*FCLP+SMALL))
     &            *((DYGF(NBS)-DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRLS,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBS))
                VCLY = C(N,NSL)*FLS*(1.D+0-THETA*DYF)*FCLP
     &            + C(NBS,NSL)*FLS*THETA*DYF*FCLS
                VCLYF = CO(N,NSL)*FLS*FCLP
                VC(NPY,NSL) = VC(NPY,NSL) + (VCLY-VCLYF)/AFY(NPY)
              ENDIF
            ELSE
              ALN = MAX( -FLN,ZERO )
              AP = (ALN+FLN)*FCLP
              AN = ALN*FCL
              VC(NQY,NSL) = VC(NQY,NSL) + (C(N,NSL)*AP - BCX(1)*AN)
            ENDIF
!
!---      Inflow gas ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (VG(1,NQY)/EPSL.LT.-EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. JFLD.GT.1 ) THEN
              VCGY = 0.D+0
              IF( FGN.LT.ZERO ) VCGY = BCX(1)*FCG*FGN
              VC(NQY,NSL) = VC(NQY,NSL) + VCGY/AFY(NQY)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FGS = AFY(NPY)*VG(1,NPY)
              IF( FGS.LT.ZERO ) THEN
                NBS = N-IFLD
                XVGX = SG(2,NBS)*PORD(2,NBS)
                CRGS = ABS( VG(1,NPY) )*DT/DYGP(NPY)/(XVGX*RP(I)+SMALL)
                FCGS = YG(NBS,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
     &            /(C(NBS,NSL)*FCGS-C(N,NSL)*FCGP+SMALL))
     &            *((DYGF(NBS)-DYGF(N))/DYGF(N))
                THETA = FLIMIT( R,CRGS,ISLC(1) )
                DYF = DYGF(N)/(DYGF(N)+DYGF(NBS))
                VCGY = C(N,NSL)*FGS*(1.D+0-THETA*DYF)*FCGP
     &            + C(NBS,NSL)*FGS*THETA*DYF*FCGS
                VCGYF = CO(N,NSL)*FGS*FCGP
                VC(NPY,NSL) = VC(NPY,NSL) + (VCGY-VCGYF)/AFY(NPY)
              ENDIF
            ELSE
              AGN = MAX( -FGN,ZERO )
              AP = (AGN+FGN)*FCGP
              AN = AGN*FCG
              VC(NQY,NSL) = VC(NQY,NSL) + (C(N,NSL)*AP - BCX(1)*AN)
            ENDIF
          ENDIF
!
!---    Top boundary
!
        ELSEIF( IBCD(NB).EQ.3 ) THEN
          NQZ = NSZ(N) + IJFLD
          IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
          NPZ = NSZ(N)
!
!---      Hydraulic dispersion
!
          IF( IDISP.EQ.1 ) THEN
            CALL ADVTB( PORD(2,N),PORDB(2,NB),SL(2,N),SLB(2,NB),
     &        UL,VL,WL,ULTX,VLTX,WLTX,N,MF )
            CALL SHDP( WLTX,ULTX,VLTX,DISPL(IZN),DISPT(IZN),DPLT )
            CALL ADVTB( PORD(2,N),PORDB(2,NB),SG(2,N),SGB(2,NB),
     &        UG,VG,WG,UGTX,VGTX,WGTX,N,MF )
            CALL SHDP( WGTX,UGTX,VGTX,DISPL(IZN),DISPT(IZN),DPGT )
          ELSE
            DPLT = 0.D+0
            DPGT = 0.D+0
          ENDIF
          FLT = AFZ(NQZ)*WL(1,NQZ)
          FGT = AFZ(NQZ)*WG(1,NQZ)
          CRLT = ABS( WL(1,NQZ) )*DT/(DZGF(N)*XVLB+SMALL)
          CRGT = ABS( WG(1,NQZ) )*DT/(DZGF(N)*XVGB+SMALL)
!
!---      Dirichlet ---
!
          IF( IBCTX.EQ.1 .OR. IBCTX.EQ.8 .OR.
     &      IBCTX.EQ.9 .OR. IBCTX.EQ.12 ) THEN
            TCOR = (TB(2,NB)+TABS)/TSPRF
            SMDLB = SMDL(NSL)*TCOR*(VISRL/VISLB(2,NB))
            DLB = TORLB(2,NB)*SLB(2,NB)*PORDB(2,NB)*SMDLB
            INDX = 16
            DLZ = DIFMN(DLP,DLB,DZGF(N),DZGF(N),WL(1,NQZ),INDX)
            DLZ = AFZ(NQZ)*(DLZ+DPLT)/(5.D-1*DZGF(N))
            PCOR = (PGB(2,NB)+PATM)/PATM
            SMDGB = SMDG(NSL)*(TCOR**1.75)/PCOR
            DGB = TORGB(2,NB)*SGB(2,NB)*PORDB(2,NB)*SMDGB
            INDX = 16
            DGZ = DIFMN(DGP,DGB,DZGF(N),DZGF(N),WG(1,NQZ),INDX)
            DGZ = AFZ(NQZ)*(DGZ+DPGT)/(5.D-1*DZGF(N))
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 ) THEN
              IF( FLT.LT.ZERO ) THEN
                WCLZ = BCX(1)*FCL*FLT
              ELSEIF( FLT.GE.ZERO .AND. K.GT.1 ) THEN
                NBB = N-IJFLD
                FCLB = YL(NBB,NSL)/(SL(2,NBB)*PORD(2,NBB)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBB,NSL)*FCLB)
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
                THETA = FLIMIT( R,CRLT,ISLC(1) )
                WCLZ =  C(N,NSL)*FLT*(1.D+0-THETA)*FCLP
     &            + BCX(1)*FLT*THETA*FCL
              ELSEIF( FLT.GE.ZERO .AND. K.EQ.1 ) THEN
                WCLZ =  C(N,NSL)*FLT*FCLP
              ENDIF
              IF( FGT.LT.ZERO ) THEN
                WCGZ = BCX(1)*FCG*FGT
              ELSEIF( FGT.GE.ZERO .AND. K.GT.1 ) THEN
                NBB = N-IJFLD
                FCGB = YG(NBB,NSL)/(SG(2,NBB)*PORD(2,NBB)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBB,NSL)*FCGB)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
                THETA = FLIMIT( R,CRGT,ISLC(1) )
                WCGZ =  C(N,NSL)*FGT*(1.D+0-THETA)*FCGP
     &            + BCX(1)*FGT*THETA*FCG
              ELSEIF( FGT.GE.ZERO .AND. K.EQ.1 ) THEN
                WCGZ =  C(N,NSL)*FGT*FCGP
              ENDIF
              AT = DLZ*FCL + DGZ*FCG
              AP = DLZ*FCLP + DGZ*FCGP
              WC(NQZ,NSL) = WC(NQZ,NSL) + (WCLZ+WCGZ)/AFZ(NQZ)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FLB = AFZ(NPZ)*WL(1,NPZ)
              IF( FLB.LT.ZERO ) THEN
                NBB = N-IJFLD
                XVLX = SL(2,NBB)*PORD(2,NBB)
                CRLB = ABS( WL(1,NPZ) )*DT/DZGP(NPZ)/(XVLX+SMALL)
                FCLB = YL(NBB,NSL)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
     &            /(C(NBB,NSL)*FCLB-C(N,NSL)*FCLP+SMALL))
     &            *((DZGF(NBB)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRLB,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
                WCLZ = C(N,NSL)*FLB*(1.D+0-THETA*DZF)*FCLP
     &            + C(NBB,NSL)*FLB*THETA*DZF*FCLB
                WCLZF = CO(N,NSL)*FLB*FCLP
                WC(NPZ,NSL) = WC(NPZ,NSL) + (WCLZ-WCLZF)/AFZ(NPZ)
              ENDIF
              FGB = AFZ(NPZ)*WG(1,NPZ)
              IF( FGB.LT.ZERO ) THEN
                NBB = N-IJFLD
                XVGX = SG(2,NBB)*PORD(2,NBB)
                CRGB = ABS( WG(1,NPZ) )*DT/DZGP(NPZ)/(XVGX+SMALL)
                FCGB = YG(NBB,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
     &            /(C(NBB,NSL)*FCGB-C(N,NSL)*FCGP+SMALL))
     &            *((DZGF(NBB)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRGB,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
                WCGZ = C(N,NSL)*FGB*(1.D+0-THETA*DZF)*FCGP
     &            + C(NBB,NSL)*FGB*THETA*DZF*FCGB
                WCGZF = CO(N,NSL)*FGB*FCGP
                WC(NPZ,NSL) = WC(NPZ,NSL) + (WCGZ-WCGZF)/AFZ(NPZ)
              ENDIF
            ELSE
              ALT = MAX( -FLT,ZERO ) +
     &          DLZ*MAX((ONE-(TENTH*ABS(FLT)/(DLZ+SMALL)))**5,ZERO)
              AGT = MAX( -FGT,ZERO ) +
     &          DGZ*MAX((ONE-(TENTH*ABS(FGT)/(DGZ+SMALL)))**5,ZERO)
              AP = (ALT+FLT)*FCLP + (AGT+FGT)*FCGP
              AT = ALT*FCL + AGT*FCG
              WC(NQZ,NSL) = WC(NQZ,NSL) + (C(N,NSL)*AP - BCX(1)*AT)
            ENDIF
!
!---      Outflow aqueous  ---
!
          ELSEIF( IBCTX.EQ.7 .OR.
     &      ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (WL(1,NQZ)/EPSL.GT.EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 ) THEN
              WCLZ = 0.D+0
              IF( FLT.GE.ZERO .AND. K.GT.1 ) THEN
                NBB = N-IJFLD
                FCLB = YL(NBB,NSL)/(SL(2,NBB)*PORD(2,NBB)+SMALL)
                R = ((C(N,NSL)*FCLP-C(NBB,NSL)*FCLB)
     &            /(BCX(1)*FCL-C(N,NSL)*FCLP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
                THETA = FLIMIT( R,CRLT,ISLC(1) )
                WCLZ =  C(N,NSL)*FLT*(1.D+0-THETA)*FCLP
     &            + BCX(1)*FLT*THETA*FCL
              ELSEIF( FLT.GE.ZERO .AND. K.EQ.1 ) THEN
                WCLZ =  C(N,NSL)*FLT*FCLP
              ENDIF
              WC(NQZ,NSL) = WC(NQZ,NSL) + WCLZ/AFZ(NQZ)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FLB = AFZ(NPZ)*WL(1,NPZ)
              IF( FLB.LT.ZERO ) THEN
                NBB = N-IJFLD
                XVLX = SL(2,NBB)*PORD(2,NBB)
                CRLB = ABS( WL(1,NPZ) )*DT/DZGP(NPZ)/(XVLX+SMALL)
                FCLB = YL(NBB,NSL)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
     &            /(C(NBB,NSL)*FCLB-C(N,NSL)*FCLP+SMALL))
     &            *((DZGF(NBB)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRLB,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
                WCLZ = C(N,NSL)*FLB*(1.D+0-THETA*DZF)*FCLP
     &            + C(NBB,NSL)*FLB*THETA*DZF*FCLB
                WCLZF = CO(N,NSL)*FLB*FCLP
                WC(NPZ,NSL) = WC(NPZ,NSL) + (WCLZ-WCLZF)/AFZ(NPZ)
              ENDIF
            ELSE
              ALT = MAX( -FLT,ZERO )
              AP = (ALT+FLT)*FCLP
              AT = ALT*FCL
              WC(NQZ,NSL) = WC(NQZ,NSL) + (C(N,NSL)*AP - BCX(1)*AT)
            ENDIF
!
!---      Outflow gas  ---
!
          ELSEIF( IBCTX.EQ.7 .OR.
     &      ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (WG(1,NQZ)/EPSL.GT.EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 ) THEN
              WCGZ = 0.D+0
              IF( FGT.GE.ZERO .AND. K.GT.1 ) THEN
                NBB = N-IJFLD
                FCGB = YG(NBB,NSL)/(SG(2,NBB)*PORD(2,NBB)+SMALL)
                R = ((C(N,NSL)*FCGP-C(NBB,NSL)*FCGB)
     &            /(BCX(1)*FCG-C(N,NSL)*FCGP+SMALL))
     &            *(DZGF(N)/(DZGF(N)+DZGF(NBB)))
                THETA = FLIMIT( R,CRGT,ISLC(1) )
                WCGZ =  C(N,NSL)*FGT*(1.D+0-THETA)*FCGP
     &            + BCX(1)*FGT*THETA*FCG
              ELSEIF( FGT.GE.ZERO .AND. K.EQ.1 ) THEN
                WCGZ =  C(N,NSL)*FGT*FCGP
              ENDIF
              WC(NQZ,NSL) = WC(NQZ,NSL) + WCLZ/AFZ(NQZ)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FGB = AFZ(NPZ)*WG(1,NPZ)
              IF( FGB.LT.ZERO ) THEN
                NBB = N-IJFLD
                XVGX = SG(2,NBB)*PORD(2,NBB)
                CRGB = ABS( WG(1,NPZ) )*DT/DZGP(NPZ)/(XVGX+SMALL)
                FCGB = YG(NBB,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
     &            /(C(NBB,NSL)*FCGB-C(N,NSL)*FCGP+SMALL))
     &            *((DZGF(NBB)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRGB,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
                WCGZ = C(N,NSL)*FGB*(1.D+0-THETA*DZF)*FCGP
     &            + C(NBB,NSL)*FGB*THETA*DZF*FCGB
                WCGZF = CO(N,NSL)*FGB*FCGP
                WC(NPZ,NSL) = WC(NPZ,NSL) + (WCGZ-WCGZF)/AFZ(NPZ)
              ENDIF
            ELSE
              AGT = MAX( -FGT,ZERO )
              AP = (AGT+FGT)*FCGP
              AT = AGT*FCG
              WC(NQZ,NSL) = WC(NQZ,NSL) + (C(N,NSL)*AP - BCX(1)*AT)
            ENDIF
!
!---      Inflow aqueous ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (WL(1,NQZ)/EPSL.LT.-EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 ) THEN
              WCLZ = 0.D+0
              IF( FLT.LT.ZERO ) WCLZ = BCX(1)*FCL*FLT
              WC(NQZ,NSL) = WC(NQZ,NSL) + WCLZ/AFZ(NQZ)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FLB = AFZ(NPZ)*WL(1,NPZ)
              IF( FLB.LT.ZERO ) THEN
                NBB = N-IJFLD
                XVLX = SL(2,NBB)*PORD(2,NBB)
                CRLB = ABS( WL(1,NPZ) )*DT/DZGP(NPZ)/(XVLX+SMALL)
                FCLB = YL(NBB,NSL)/(XVLX+SMALL)
                R = ((C(N,NSL)*FCLP-BCX(1)*FCL)
     &            /(C(NBB,NSL)*FCLB-C(N,NSL)*FCLP+SMALL))
     &            *((DZGF(NBB)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRLB,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
                WCLZ = C(N,NSL)*FLB*(1.D+0-THETA*DZF)*FCLP
     &            + C(NBB,NSL)*FLB*THETA*DZF*FCLB
                WCLZF = CO(N,NSL)*FLB*FCLP
                WC(NPZ,NSL) = WC(NPZ,NSL) + (WCLZ-WCLZF)/AFZ(NPZ)
              ENDIF
            ELSE
              ALT = MAX( -FLT,ZERO )
              AP = (ALT+FLT)*FCLP
              AT = ALT*FCL
              WC(NQZ,NSL) = WC(NQZ,NSL) + (C(N,NSL)*AP - BCX(1)*AT)
            ENDIF
!
!---      Inflow gas ---
!
          ELSEIF( IBCTX.EQ.13 .OR. IBCTX.EQ.14 .OR. IBCTX.EQ.15
     &       .OR. ((IBCTX.EQ.19 .OR. IBCTX.EQ.23 .OR. IBCTX.EQ.43)
     &      .AND. (WG(1,NQZ)/EPSL.LT.-EPSL)) ) THEN
!
!---        TVD Transport for the boundary surface  ---
!
            IF( ISLC(1).GE.1 .AND. KFLD.GT.1 ) THEN
              WCGZ = 0.D+0
              IF( FGT.LT.ZERO ) WCGZ = BCX(1)*FCG*FGT
              WC(NQZ,NSL) = WC(NQZ,NSL) + WCGZ/AFZ(NQZ)
!
!---          TVD Transport for interior surface adjacent 
!             to boundary  ---
!
              FGB = AFZ(NPZ)*WG(1,NPZ)
              IF( FGB.LT.ZERO ) THEN
                NBB = N-IJFLD
                XVGX = SG(2,NBB)*PORD(2,NBB)
                CRGB = ABS( WG(1,NPZ) )*DT/DZGP(NPZ)/(XVGX+SMALL)
                FCGB = YG(NBB,NSL)/(XVGX+SMALL)
                R = ((C(N,NSL)*FCGP-BCX(1)*FCG)
     &            /(C(NBB,NSL)*FCGB-C(N,NSL)*FCGP+SMALL))
     &            *((DZGF(NBB)+DZGF(N))/DZGF(N))
                THETA = FLIMIT( R,CRGB,ISLC(1) )
                DZF = DZGF(N)/(DZGF(N)+DZGF(NBB))
                WCGZ = C(N,NSL)*FGB*(1.D+0-THETA*DZF)*FCGP
     &            + C(NBB,NSL)*FGB*THETA*DZF*FCGB
                WCGZF = CO(N,NSL)*FGB*FCGP
                WC(NPZ,NSL) = WC(NPZ,NSL) + (WCGZ-WCGZF)/AFZ(NPZ)
              ENDIF
            ELSE
              AGT = MAX( -FGT,ZERO )
              AP = (AGT+FGT)*FCGP
              AT = AGT*FCG
              WC(NQZ,NSL) = WC(NQZ,NSL) + (C(N,NSL)*AP - BCX(1)*AT)
            ENDIF
          ENDIF
        ENDIF
  200 CONTINUE
!
!---  End of SFXB32 group  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

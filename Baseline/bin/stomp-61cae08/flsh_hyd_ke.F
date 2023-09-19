!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_IC1( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for initial condition #1
!
!     SHX = 0.0
!     SNX = 0.0
!     SGX = 0.0
!
!     Declared variables:
!
!     TX - temperature, C
!     PLX - aqueous pressure, Pa
!     YMMOX - mobile CH4 mole fraction of formers (input as PVOX )
!     WLFX - relative saturation of dissolved hydrate 
!             formers (input as PVAX )
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 August 2010.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NCG_PT
      USE NAPL
      USE HYST
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
      REAL*8 RKLX(3)
      SAVE IALLOC,AJ,BJ,GX,IJ
      DATA IALLOC /0/
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: AJ,GX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: BJ
      INTEGER, DIMENSION(:), ALLOCATABLE :: IJ
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_IC1'
!
!---  Dynamic memory allocation  ---
!
      IF( IALLOC.EQ.0 ) THEN
        ALLOCATE( AJ(1:2,1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: AJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( BJ(1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: BJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( GX(1:2,1:3),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: GX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( IJ(1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: IJ'
          CALL WRMSGS( INDX )
        ENDIF
        IALLOC = 1
      ENDIF
!
!---  Assign mobile CH4 mole fraction of formers and
!     dissolved CH4 relative saturation  ---
!
      YMMOX = MIN( MAX( PVOX,0.D+0 ),1.D+0 )
      WLFX = MIN( MAX( PVAX,0.D+0 ),1.D+0 )
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure without hydrate barrier adjustment  ---
!
      ENPR = SCHR(2,IZN)*RHORL*GRAV
!
!---  Gas pressure  ---
!
      PGX = PLX+ENPR/BGL
!
!---  Liquid CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      PNX = PNC
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Guess CO2 and CH4 gas partial pressure  ---
!
      CALL SP_B( TX,XLSX,PSBX )
      PVWX = PSBX
      PVX = PGX+PATM-PVWX
      PVAX = (PVX*WLFX)*(1.D+0-YMMOX)
      PVOX = (PVX*WLFX)*YMMOX
      IF( WLFX/EPSL.LT.EPSL ) GOTO 40
!
!---  Two-variable Newton-Raphson loop (PVAX,PVOX)  ---
!
      NC = 0
   10 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 17
        CHMSG = 'Unconverged FLH_IC1: Mobile CO2 Vapor ' //
     &    'Partial Pressure and Mobile CH4 Vapor Partial ' //
     &    'Pressure'
        CALL WRMSGS( INDX )
      ENDIF
      DO 20 M = 1,3
        PVAY = PVAX
        PVOY = PVOX
        DPVOX = SIGN( 1.D-1,5.D-1*PVX-PVOX )
        DPVAX = SIGN( 1.D-1,5.D-1*PVX-PVAX )
        IF( M.EQ.2 ) PVAY = PVAX + DPVAX
        IF( M.EQ.3 ) PVOY = PVOX + DPVOX
        PVOZ = (PVOY/(PVOY+PVAY))*PVX
        PVAZ = (PVAY/(PVOY+PVAY))*PVX
        INCG = 1
        CALL SOL_LA( TX,PVAY,XLSX,XMLAX )
        INCG = 2
        CALL SOL_LO( TX,PVOY,XLSX,XMLOX )
        GX(1,M) = YMMOX - XMLOX/(XMLOX+XMLAX+SMALL)
        GX(2,M) = WLFX - (PVAY+PVOY)/(PVX+SMALL)
   20 CONTINUE
!
!---  Load solution vector and Jacobian matrix for
!     two-phase conditions  ---
!
      DO 30 M = 1,2
        AJ(M,1) = (GX(M,2)-GX(M,1))/DPVAX
        AJ(M,2) = (GX(M,3)-GX(M,1))/DPVOX
        BJ(M) = -GX(M,1)
  30  CONTINUE
!
!---  Solve linear system  ---
!
      JP = 2
      KP = 2
      CALL LUDCMP( AJ,JP,KP,IJ,DJ )
      CALL LUBKSB( AJ,JP,KP,IJ,BJ )
!
!---  Update primary unknowns  ---
!
      DPVAX = BJ(1)
      DPVOX = BJ(2)
      PVAX = PVAX + DPVAX
      PVOX = PVOX + DPVOX
!
!---  Upper limit on CH4 and CH2 vapor partial pressures  ---
!
      PVAX = MIN( PVAX,PVX )
      PVOX = MIN( PVOX,PVX )
!
!---  Lower limit on CH4 and CH2 vapor partial pressures  ---
!
      IF( PVAX/EPSL.LT.EPSL .AND. DPVAX.LT.0.D+0 ) THEN
        PVAX = 0.D+0
        DPVAX = 0.D+0
      ENDIF
      IF( PVOX/EPSL.LT.EPSL .AND. DPVOX.LT.0.D+0 ) THEN
        PVOX = 0.D+0
        DPVOX = 0.D+0
      ENDIF
!
!---  Liquid-CO2 conditions  ---
!
      CALL SP_A( TX,PSAX )
      TKX = TX + TABS
      IF( PVAX.GT.PSAX .AND. PVAX.LT.PCRA
     &    .AND. TKX.GT.TTPA .AND. TKX.LT.TCRA ) THEN
        INDX = 11
        RLMSG = PVAX - PSAX
        CHMSG = 'FLH_IC1: Liquid-CO2 Conditions: PVAX - PSAX = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Convergence check  ---
!
      IF( ABS(DPVAX).GT.1.D-1 ) GOTO 10
      IF( ABS(DPVOX).GT.1.D-1 ) GOTO 10
!
!---  Check for hydrate conditions  ---
!
   40 CONTINUE
!
!---  Hydrate equilibrium pressure  ---
!
      IF( (PVAX+PVOX).GT.EPSL ) THEN
        PVHOX = PVOX
        PVHAX = PVAX
        YMHGOX = PVHOX/(PVHAX+PVHOX)
        PEQHX = PVHOX+PVHAX
!
!---    Inhibitor correction  ---
!
        IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
          CALL HYD_PO( YMHGOX,POHX )
        ENDIF
        CALL HYD_T( TEQHX,YMHGOX,XLSX,PEQHX,POHX,RHOHX,XHAX,XHOX,XHWX )
      ENDIF
!
!---  Hydrate conditions  ---
!
      IF( TX.LE.TEQHX ) THEN
        INDX = 11
        RLMSG = TX - TEQHX
        CHMSG = 'FLH_IC1: Hydrate Conditions: T - TEQH = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Hydrate properties for no-hydrate conditions  ---
!
      TMHAX = 0.D+0
      TMHOX = 0.D+0
      YMHGOX = 0.D+0
      TEQHX = -TABS
      XHAX = 0.D+0
      XHOX = 0.D+0
      XHWX = 1.D+0
      RHOHX = 0.D+0
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_IC1 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_IC2( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for initial condition #2
!
!     SHX = 0.0
!     SNX = 0.0
!     SGX > 0.0
!
!     Declared variables:
!
!     TX - temperature, C
!     PLX - aqueous pressure, Pa
!     YMMOX - mobile CH4 mole fraction of formers (input as PVOX )
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 August 2010.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NCG_PT
      USE NAPL
      USE HYST
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
      REAL*8 RKLX(3)
      SAVE IALLOC,AJ,BJ,GX,IJ
      DATA IALLOC /0/
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: AJ,GX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: BJ
      INTEGER, DIMENSION(:), ALLOCATABLE :: IJ#endif
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_IC2'
!
!---  Dynamic memory allocation  ---
!
      IF( IALLOC.EQ.0 ) THEN
        ALLOCATE( AJ(1:2,1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: AJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( BJ(1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: BJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( GX(1:2,1:3),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: GX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( IJ(1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: IJ'
          CALL WRMSGS( INDX )
        ENDIF
        IALLOC = 1
      ENDIF
!
!---  Assign mobile CH4 mole fraction of formers  ---
!
      YMMOX = MIN( MAX( PVOX,0.D+0 ),1.D+0 )
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure without hydrate barrier adjustment  ---
!
      ENPR = SCHR(2,IZN)*RHORL*GRAV
!
!---  Check for no gas phase conditions  ---
!
      IF( PGX.LE.(PLX+ENPR/BGL) ) THEN
        INDX = 11
        RLMSG = PGX - (PLX+ENPR/BGL)
        CHMSG = 'FLH_IC2: No Gas: PG - (PL + ENPR/BGN) = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Liquid CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      PNX = PNC
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Guess CO2 and CH4 gas partial pressure  ---
!
      CALL SP_B( TX,XLSX,PSBX )
      PPLX = MAX( PLX+PATM,PSBX )
      CALL DENS_B( TX,PPLX,XLSX,RHOBX )
      PCX = MAX( BGL*(PGX-PLX),0.D+0 )
      IF( ISLC(44).EQ.1 ) THEN
        PVWX = PSBX
      ELSE
        CALL VPL_B( TX,PSBX,PCX,RHOBX,PVWX,XLSX )
      ENDIF
      PVX = PGX+PATM-PVWX
      PVAX = PVX*(1.D+0-YMMOX)
      PVOX = PVX*YMMOX
!
!---  Single-variable Newton-Raphson loop (PVOX)  ---
!
      NC = 0
   10 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 17
        CHMSG = 'Unconverged FLH_IC2: CH4 Vapor Partial Pressure'
        CALL WRMSGS( INDX )
      ENDIF
      DO 20 M = 1,2
        PVOY = PVOX
        DPVOX = SIGN( 1.D-1,5.D-1*PVX-PVOX )
        IF( M.EQ.2 ) PVOY = PVOX + DPVOX
        PVAX = PVX-PVOY
        INCG = 1
        CALL SOL_LA( TX,PVAX,XLSX,XMLAX )
        INCG = 2
        CALL SOL_LO( TX,PVOY,XLSX,XMLOX )
        WTMLX = (XMLOX*WTMO + XMLAX*WTMA + 
     &    WTMW*(1.D+0-XMLOX-XMLAX))/
     &    (1.D+0-XLSX+WTMW*XLSX/WTMS)
        XLAX = XMLAX*WTMA/WTMLX
        XLOX = XMLOX*WTMO/WTMLX
        CALL DENS_B( TX,PX,XLSX,RHOBX )
        CALL DENS_L( TX,RHOBX,XLAX,XLOX,RHOLX )
        XLWX = MAX( 1.D+0-XLAX-XLOX-XLSX,0.D+0 )
        XMLSX = XLSX/WTMS/WTMLX
        XMLWX = XLWX/WTMW/WTMLX
        RHOMLX = RHOLX/WTMLX
        ISRX = 2
        CALL DENS_W( TX,PVWX,RHOX,RHOGWX,ISRX )
        CALL DENS_A( TX,PVAX,RHOGAX,ISRX )
        CALL DENS_O( TX,PVOY,RHOGOX )
        RHOGX = RHOGWX+RHOGAX+RHOGOX
        XGAX = RHOGAX/RHOGX
        XGOX = RHOGOX/RHOGX
        XGWX = RHOGWX/RHOGX
        WTMGX = XGAX/WTMA + XGOX/WTMO + XGWX/WTMW
        XMGAX = XGAX/WTMA/WTMGX
        XMGOX = XGOX/WTMO/WTMGX
        XMGWX = XGWX/WTMW/WTMGX
        RHOMGX = RHOGX*WTMGX
        XMTOX = XMLOX*SLX*RHOMLX + XMGOX*SGX*RHOMGX
        XMTAX = XMLAX*SLX*RHOMLX + XMGAX*SGX*RHOMGX
        GX(1,M) = YMMOX - XMTOX/(XMTOX+XMTAX)
   20 CONTINUE
!
!---  Linearization  ---
!
      DFX = (GX(1,2)-GX(1,1))/DPVOX
      FX = -GX(1,1)
!
!---  Update primary unknowns  ---
!
      DPVOX = FX/DFX
      DPX = MAX( 1.D+3,1.D-1*PVX )
      DPVOX = SIGN( MIN( ABS(DPX),ABS(DPVOX) ),DPVOX )
      PVOX = PVOX + DPVOX
!
!---  Upper limit on CH4 vapor partial pressure  ---
!
      PVOX = MIN( PVOX,PVX )
!
!---  Lower limit on CH4 vapor partial pressure  ---
!
      IF( PVOX/EPSL.LT.EPSL .AND. DPVOX.LT.0.D+0 ) THEN
        PVOX = 0.D+0
        DPVOX = 0.D+0
      ENDIF
!
!---  Convergence check  ---
!
      IF( ABS(DPVOX).GT.1.D-1 ) GOTO 10
!
!---  Gas pressure  ---
!
      PVAX = MAX( PVX-PVOX,0.D+0 )
!
!---  Liquid-CO2 conditions  ---
!
      CALL SP_A( TX,PSAX )
      TKX = TX + TABS
      IF( PVAX.GT.PSAX .AND. PVAX.LT.PCRA
     &    .AND. TKX.GT.TTPA .AND. TKX.LT.TCRA ) THEN
        INDX = 11
        RLMSG = PVAX - PSAX
        CHMSG = 'FLH_IC2: Liquid-CO2 Conditions: PVAX - PSAX = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Hydrate equilibrium pressure  ---
!
      IF( (PVAX+PVOX).GT.EPSL ) THEN
        PVHOX = PVOX
        PVHAX = PVAX
        YMHGOX = PVHOX/(PVHAX+PVHOX)
        PEQHX = PVHOX+PVHAX
!
!---    Inhibitor correction  ---
!
        IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
          CALL HYD_PO( YMHGOX,POHX )
        ENDIF
        CALL HYD_T( TEQHX,YMHGOX,XLSX,PEQHX,POHX,RHOHX,XHAX,XHOX,XHWX )
      ENDIF
!
!---  Hydrate conditions  ---
!
      IF( TX.LE.TEQHX ) THEN
        INDX = 11
        RLMSG = TX - TEQHX
        CHMSG = 'FLH_IC2: Hydrate Conditions: T - TEQH = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Hydrate properties for no-hydrate conditions  ---
!
      TMHAX = 0.D+0
      TMHOX = 0.D+0
      YMHGOX = 0.D+0
      TEQHX = -TABS
      XHAX = 0.D+0
      XHOX = 0.D+0
      XHWX = 1.D+0
      RHOHX = 0.D+0
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_IC2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_IC3( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for initial condition #3
!
!     SHX > 0.0
!     SNX = 0.0
!     SGX = 0.0
!
!     Declared variables:
!
!     TX - temperature, C
!     PLX - aqueous pressure, Pa
!     YMMOX - mobile CH4 mole fraction of formers (input as PVOX)
!     WLFX - relative saturation of dissolved hydrate 
!             formers (input as PVAX )
!     YMHOX - hydrate CH4 mole fraction of formers (input as PVHOX)
!     SHX - hydrate saturation
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 August 2010.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NCG_PT
      USE NAPL
      USE HYST
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLX(3)
      SAVE IALLOC,AJ,BJ,GX,IJ
      DATA IALLOC /0/
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: AJ,GX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: BJ
      INTEGER, DIMENSION(:), ALLOCATABLE :: IJ
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_IC3'
!
!---  Dynamic memory allocation  ---
!
      IF( IALLOC.EQ.0 ) THEN
        ALLOCATE( AJ(1:2,1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: AJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( BJ(1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: BJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( GX(1:2,1:3),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: GX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( IJ(1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: IJ'
          CALL WRMSGS( INDX )
        ENDIF
        IALLOC = 1
      ENDIF
!
!---  Assign mobile CH4 mole fraction of formers,
!     dissolved CH4 relative saturation, and
!     hydrate CH4 mole fraction of formers  ---
!
      YMMOX = MIN( MAX( PVOX,0.D+0 ),1.D+0 )
      WLFX = MIN( MAX( PVAX,0.D+0 ),1.D+0 )
      YMHOX = MIN( MAX( PVHOX,0.D+0 ),1.D+0 )
!
!---  Limit hydrate saturation  ---
!
      SHX = MIN( MAX( SHX,0.D+0 ),1.D+0 )
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure with hydrate barrier adjustment  ---
!
      SHZ = MAX( SHX,1.D-9 )
      HCPBFX = 1.D+0
      IF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.12 ) THEN
        HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &    SIN(SCHR(12,IZN)*SHZ + SCHR(13,IZN)) )
      ENDIF
      ENPR = SCHR(2,IZN)*RHORL*GRAV*HCPBFX
!
!---  Pore-space equilibrium temperature  ---
!
      TEQHX = TX
!
!---  Guess hydrate-gas CH4 mole fraction of formers  ---
!
      YMHGOX = YMHOX
!
!---  Single-variable Newton-Raphson loop (YMHGOX)  ---
!
      NC = 0
   10 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 17
        CHMSG = 'Unconverged FLH_IC3: Hydrate-Gas CH4 ' //
     &     'Mole Fraction of Formers'
        CALL WRMSGS( INDX )
      ENDIF
      DO 20 M = 1,2
        YMHGOY = YMHGOX
        DYMHGOX = SIGN( 1.D-6,(1.D-1-YMHGOX) )
        IF( M.EQ.2 ) YMHGOY = YMHGOX + DYMHGOX
        CALL HYD_P( TEQHX,YMHGOY,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
        XMHOX = XHOX/WTMO
        XMHAX = XHAX/WTMA
        GX(M,1) = YMHOX - XMHOX/(XMHOX+XMHAX)
   20 CONTINUE
!
!---  Linearization  ---
!
      DFX = (GX(2,1)-GX(1,1))/DYMHGOX
      FX = -GX(1,1)
!
!---  Update primary unknowns  ---
!
      DYMHGOX = FX/DFX
      YMHGOX = YMHGOX + DYMHGOX
!
!---  Limit hydrate-gas CH4 mole fraction of formers  ---
!
      YMHGOX = MIN( MAX( YMHGOX,0.D+0 ),1.D+0 )
!
!---  Convergence check  ---
!
      IF( ABS(DYMHGOX).GT.1.D-7 ) GOTO 10
!
!---  Converged solution  ---
!
      CALL HYD_P( TEQHX,YMHGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
      PVHX = PEQHX
      PVHOX = PVHX*YMHGOX
      PVHAX = PVHX-PVHOX
!
!---  Gas pressure  ---
!
      PGX = PLX+ENPR/BGL
!
!---  Liquid CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      PNX = PNC
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Guess CO2 and CH4 gas partial pressure  ---
!
      CALL SP_B( TX,XLSX,PSBX )
      PPLX = MAX( PLX+PATM,PSBX )
      CALL DENS_B( TX,PPLX,XLSX,RHOBX )
      PVWX = PSBX
      PVX = PGX+PATM-PVWX
      PVAX = (PVX*WLFX)*(1.D+0-YMMOX)
      PVOX = (PVX*WLFX)*YMMOX
      IF( WLFX/EPSL.LT.EPSL ) GOTO 60
!
!---  Two-variable Newton-Raphson loop (PVAX,PVOX)  ---
!
      NC = 0
   30 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 17
        CHMSG = 'Unconverged FLH_IC3: Mobile CO2 Vapor ' //
     &    'Partial Pressure and Mobile CH4 Vapor Partial ' //
     &    'Pressure'
        CALL WRMSGS( INDX )
      ENDIF
      DO 40 M = 1,3
        PVAY = PVAX
        PVOY = PVOX
        DPVOX = SIGN( 1.D-1,5.D-1*PVX-PVOX )
        DPVAX = SIGN( 1.D-1,5.D-1*PVX-PVAX )
        IF( M.EQ.2 ) PVAY = PVAX + DPVAX
        IF( M.EQ.3 ) PVOY = PVOX + DPVOX
        PVOZ = (PVOY/(PVOY+PVAY))*PVX
        PVAZ = (PVAY/(PVOY+PVAY))*PVX
        INCG = 1
        CALL SOL_LA( TX,PVAY,XLSX,XMLAX )
!        CALL SOL_LA( TX,PVAZ,XLSX,XMLAZ )
        INCG = 2
        CALL SOL_LO( TX,PVOY,XLSX,XMLOX )
!        CALL SOL_LO( TX,PVOZ,XLSX,XMLOZ )
        GX(1,M) = YMMOX - XMLOX/(XMLOX+XMLAX+SMALL)
        GX(2,M) = WLFX - (PVAY+PVOY)/(PVX+SMALL)
!        IF( PVOX.EQ.0 ) THEN
!          GX(2,M) = WLFX - XMLAX/XMLAZ
!        ELSE
!          GX(2,M) = WLFX - XMLOX/XMLOZ
!        ENDIF
   40 CONTINUE
!
!---  Load solution vector and Jacobian matrix for
!     two-phase conditions  ---
!
      DO 50 M = 1,2
        AJ(M,1) = (GX(M,2)-GX(M,1))/DPVAX
        AJ(M,2) = (GX(M,3)-GX(M,1))/DPVOX
        BJ(M) = -GX(M,1)
  50  CONTINUE
!
!---  Solve linear system  ---
!
      JP = 2
      KP = 2
      CALL LUDCMP( AJ,JP,KP,IJ,DJ )
      CALL LUBKSB( AJ,JP,KP,IJ,BJ )
!
!---  Update primary unknowns  ---
!
      DPVAX = BJ(1)
      DPVOX = BJ(2)
      PVAX = PVAX + DPVAX
      PVOX = PVOX + DPVOX
!
!---  Upper limit on CH4 and CH2 vapor partial pressures  ---
!
      PVAX = MIN( PVAX,PVX )
      PVOX = MIN( PVOX,PVX )
!
!---  Lower limit on CH4 and CH2 vapor partial pressures  ---
!
      IF( PVAX/EPSL.LT.EPSL .AND. DPVAX.LT.0.D+0 ) THEN
        PVAX = 0.D+0
        DPVAX = 0.D+0
      ENDIF
      IF( PVOX/EPSL.LT.EPSL .AND. DPVOX.LT.0.D+0 ) THEN
        PVOX = 0.D+0
        DPVOX = 0.D+0
      ENDIF
!
!---  Convergence check  ---
!
      IF( ABS(DPVAX).GT.1.D-1 ) GOTO 30
      IF( ABS(DPVOX).GT.1.D-1 ) GOTO 30
!
!---  Converged solution  ---
!
   60 CONTINUE
!
!---  Liquid-CO2 conditions  ---
!
      CALL SP_A( TX,PSAX )
      TKX = TX + TABS
      IF( PVAX.GT.PSAX .AND. PVAX.LT.PCRA
     &    .AND. TKX.GT.TTPA .AND. TKX.LT.TCRA ) THEN
        INDX = 11
        RLMSG = PVAX - PSAX
        CHMSG = 'FLH_IC3: Liquid-CO2 Conditions: PVAX - PSAX = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_IC3 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_IC4( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for initial condition #4
!
!     SHX > 0.0
!     SNX = 0.0
!     SGX > 0.0
!
!     Declared variables:
!
!     TX - temperature, C
!     PLX - aqueous pressure, Pa
!     YMMOX - mobile CH4 mole fraction of formers (input as PVOX )
!     YMHOX - hydrate CH4 mole fraction of formers (input as PVHOX)
!     SHX - hydrate saturation
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 August 2010.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NCG_PT
      USE NAPL
      USE HYST
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
      REAL*8 RKLX(3)
      SAVE IALLOC,AJ,BJ,GX,IJ
      DATA IALLOC /0/
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: AJ,GX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: BJ
      INTEGER, DIMENSION(:), ALLOCATABLE :: IJ
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_IC4'
!
!---  Dynamic memory allocation  ---
!
      IF( IALLOC.EQ.0 ) THEN
        ALLOCATE( AJ(1:2,1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: AJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( BJ(1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: BJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( GX(1:2,1:3),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: GX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( IJ(1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: IJ'
          CALL WRMSGS( INDX )
        ENDIF
        IALLOC = 1
      ENDIF
!
!---  Assign mobile CH4 mole fraction of formers and
!     hydrate CH4 mole fraction of formers  ---
!
      YMMOX = MIN( MAX( PVOX,0.D+0 ),1.D+0 )
      YMHOX = MIN( MAX( PVHOX,0.D+0 ),1.D+0 )
!
!---  Limit hydrate saturation  ---
!
      SHX = MIN( MAX( SHX,0.D+0 ),1.D+0 )
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure with hydrate barrier adjustment  ---
!
      SHZ = MAX( SHX,1.D-9 )
      HCPBFX = 1.D+0
      IF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.12 ) THEN
        HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &    SIN(SCHR(12,IZN)*SHZ + SCHR(13,IZN)) )
      ENDIF
      ENPR = SCHR(2,IZN)*RHORL*GRAV*HCPBFX
!
!---  Limit minimum aqueous pressure  ---
!
      PLX = MAX( PLX,(PGX-HDOD*HCPBFX*RHORL*GRAV/BGL) )
!
!---  Check for no gas phase conditions  ---
!
      IF( PGX.LE.(PLX+ENPR/BGL) ) THEN
        INDX = 11
        RLMSG = PGX - (PLX+ENPR/BGL)
        CHMSG = 'FLH_IC4: No Gas: PG - (PL + ENPR/BGN) = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Pore-space equilibrium temperature  ---
!
      TEQHX = TX
!
!---  Guess hydrate-gas CH4 mole fraction of formers  ---
!
      YMHGOX = YMHOX
!
!---  Single-variable Newton-Raphson loop (YMHGOX)  ---
!
      NC = 0
   10 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 17
        CHMSG = 'Unconverged FLH_IC4: Hydrate-Gas CH4 ' //
     &     'Mole Fraction of Formers'
        CALL WRMSGS( INDX )
      ENDIF
      DO 20 M = 1,2
        YMHGOY = YMHGOX
        DYMHGOX = SIGN( 1.D-6,(1.D-1-YMHGOX) )
        IF( M.EQ.2 ) YMHGOY = YMHGOX + DYMHGOX
        CALL HYD_P( TEQHX,YMHGOY,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
        XMHOX = XHOX/WTMO
        XMHAX = XHAX/WTMA
        GX(M,1) = YMHOX - XMHOX/(XMHOX+XMHAX)
   20 CONTINUE
!
!---  Linearization  ---
!
      DFX = (GX(2,1)-GX(1,1))/DYMHGOX
      FX = -GX(1,1)
!
!---  Update primary unknowns  ---
!
      DYMHGOX = FX/DFX
      YMHGOX = YMHGOX + DYMHGOX
!
!---  Limit hydrate-gas CH4 mole fraction of formers  ---
!
      YMHGOX = MIN( MAX( YMHGOX,0.D+0 ),1.D+0 )
!
!---  Convergence check  ---
!
      IF( ABS(DYMHGOX).GT.1.D-7 ) GOTO 10
!
!---  Converged solution  ---
!
      CALL HYD_P( TEQHX,YMHGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
      PVHX = PEQHX
      PVHOX = PVHX*YMHGOX
      PVHAX = PVHX-PVHOX
!
!---  Liquid CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      PNX = PNC
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Guess CO2 and CH4 gas partial pressure  ---
!
      CALL SP_B( TX,XLSX,PSBX )
      PPLX = MAX( PLX+PATM,PSBX )
      CALL DENS_B( TX,PPLX,XLSX,RHOBX )
      PCX = MAX( BGL*(PGX-PLX),0.D+0 )
      IF( ISLC(44).EQ.1 ) THEN
        PVWX = PSBX
      ELSE
        CALL VPL_B( TX,PSBX,PCX,RHOBX,PVWX,XLSX )
      ENDIF
      PVX = PGX+PATM-PVWX
      PVAX = PVX*(1.D+0-YMMOX)
      PVOX = PVX*YMMOX
!
!---  Single-variable Newton-Raphson loop (PVOX)  ---
!
      NC = 0
   30 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 17
        CHMSG = 'Unconverged FLH_IC4: Mobile CH4 Vapor ' //
     &    'Partial Pressure'
        CALL WRMSGS( INDX )
      ENDIF
      DO 40 M = 1,2
        PVOY = PVOX
        DPVOX = SIGN( 1.D-1,5.D-1*PVX-PVOX )
        IF( M.EQ.2 ) PVOY = PVOX + DPVOX
        PVAX = PVX-PVOY
        INCG = 1
        CALL SOL_LA( TX,PVAX,XLSX,XMLAX )
        INCG = 2
        CALL SOL_LO( TX,PVOY,XLSX,XMLOX )
        WTMLX = (XMLOX*WTMO + XMLAX*WTMA + 
     &    WTMW*(1.D+0-XMLOX-XMLAX))/
     &    (1.D+0-XLSX+WTMW*XLSX/WTMS)
        XLAX = XMLAX*WTMA/WTMLX
        XLOX = XMLOX*WTMO/WTMLX
        CALL DENS_B( TX,PX,XLSX,RHOBX )
        CALL DENS_L( TX,RHOBX,XLAX,XLOX,RHOLX )
        XLWX = MAX( 1.D+0-XLAX-XLOX-XLSX,0.D+0 )
        XMLSX = XLSX/WTMS/WTMLX
        XMLWX = XLWX/WTMW/WTMLX
        RHOMLX = RHOLX/WTMLX
        ISRX = 2
        CALL DENS_W( TX,PVWX,RHOX,RHOGWX,ISRX )
        CALL DENS_A( TX,PVAX,RHOGAX,ISRX )
        CALL DENS_O( TX,PVOY,RHOGOX )
        RHOGX = RHOGWX+RHOGAX+RHOGOX
        XGAX = RHOGAX/RHOGX
        XGOX = RHOGOX/RHOGX
        XGWX = RHOGWX/RHOGX
        WTMGX = XGAX/WTMA + XGOX/WTMO + XGWX/WTMW
        XMGAX = XGAX/WTMA/WTMGX
        XMGOX = XGOX/WTMO/WTMGX
        XMGWX = XGWX/WTMW/WTMGX
        RHOMGX = RHOGX*WTMGX
        XMTOX = XMLOX*SLX*RHOMLX + XMGOX*SGX*RHOMGX
        XMTAX = XMLAX*SLX*RHOMLX + XMGAX*SGX*RHOMGX
        GX(1,M) = YMMOX - XMTOX/(XMTOX+XMTAX)
   40 CONTINUE
!
!---  Load solution vector and Jacobian matrix for
!     two-phase conditions  ---
!
      DFX = (GX(1,2)-GX(1,1))/DPVOX
      FX = -GX(1,1)
!
!---  Update primary unknowns  ---
!
      DPVOX = FX/DFX
      PVOX = PVOX + DPVOX
!
!---  Upper limit on CH4 vapor partial pressure  ---
!
      PVOX = MIN( PVOX,PVX )
!
!---  Lower limit on CH4 vapor partial pressure  ---
!
      IF( PVOX/EPSL.LT.EPSL .AND. DPVOX.LT.0.D+0 ) THEN
        PVOX = 0.D+0
        DPVOX = 0.D+0
      ENDIF
!
!---  Convergence check  ---
!
      IF( ABS(DPVOX).GT.1.D-1 ) GOTO 30
!
!---  Gas pressure  ---
!
      PVAX = MAX( PVX-PVOX,0.D+0 )
!
!---  Liquid-CO2 conditions  ---
!
      CALL SP_A( TX,PSAX )
      TKX = TX + TABS
      IF( PVAX.GT.PSAX .AND. PVAX.LT.PCRA
     &    .AND. TKX.GT.TTPA .AND. TKX.LT.TCRA ) THEN
        INDX = 11
        RLMSG = PVAX - PSAX
        CHMSG = 'FLH_IC4: Liquid-CO2 Conditions: PVAX - PSAX = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_IC4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_IC5( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for initial condition #5
!
!     SHX > 0.0
!     SNX = 0.0
!     SGX ³ 0.0
!     Hydrate-mobile phase equilibrium
!
!     Declared variables:
!
!     TX - temperature, C
!     PLX - aqueous pressure, Pa
!     YMHOX - hydrate CH4 mole fraction of formers (input as PVHOX)
!     SHX - hydrate saturation
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 August 2010.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NCG_PT
      USE NAPL
      USE HYST
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLX(3)
      SAVE IALLOC,AJ,BJ,GX,IJ
      DATA IALLOC /0/
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: AJ,GX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: BJ
      INTEGER, DIMENSION(:), ALLOCATABLE :: IJ
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_IC5'
!
!---  Dynamic memory allocation  ---
!
      IF( IALLOC.EQ.0 ) THEN
        ALLOCATE( AJ(1:2,1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: AJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( BJ(1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: BJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( GX(1:2,1:3),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: GX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( IJ(1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: IJ'
          CALL WRMSGS( INDX )
        ENDIF
        IALLOC = 1
      ENDIF
!
!---  Assign hydrate CH4 mole fraction of formers  ---
!
      YMHOX = MIN( MAX( PVHOX,0.D+0 ),1.D+0 )
!
!---  Limit hydrate saturation  ---
!
      SHX = MIN( MAX( SHX,0.D+0 ),1.D+0 )
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure with hydrate barrier adjustment  ---
!
      SHZ = MAX( SHX,1.D-9 )
      HCPBFX = 1.D+0
      IF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.12 ) THEN
        HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &    SIN(SCHR(12,IZN)*SHZ + SCHR(13,IZN)) )
      ENDIF
      ENPR = SCHR(2,IZN)*RHORL*GRAV*HCPBFX
!
!---  Pore-space equilibrium temperature  ---
!
      TEQHX = TX
!
!---  Guess hydrate-gas CH4 mole fraction of formers  ---
!
      YMHGOX = YMHOX
!
!---  Single-variable Newton-Raphson loop (YMHGOX)  ---
!
      NC = 0
   10 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 17
        CHMSG = 'Unconverged FLH_IC5: Hydrate-Gas CH4 ' //
     &     'Mole Fraction of Formers'
        CALL WRMSGS( INDX )
      ENDIF
      DO 20 M = 1,2
        YMHGOY = YMHGOX
        DYMHGOX = SIGN( 1.D-6,(1.D-1-YMHGOX) )
        IF( M.EQ.2 ) YMHGOY = YMHGOX + DYMHGOX
        CALL HYD_P( TEQHX,YMHGOY,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
        XMHOX = XHOX/WTMO
        XMHAX = XHAX/WTMA
        GX(M,1) = YMHOX - XMHOX/(XMHOX+XMHAX)
   20 CONTINUE
!
!---  Linearization  ---
!
      DFX = (GX(2,1)-GX(1,1))/DYMHGOX
      FX = -GX(1,1)
!
!---  Update primary unknowns  ---
!
      DYMHGOX = FX/DFX
      YMHGOX = YMHGOX + DYMHGOX
!
!---  Limit hydrate-gas CH4 mole fraction of formers  ---
!
      YMHGOX = MIN( MAX( YMHGOX,0.D+0 ),1.D+0 )
!
!---  Convergence check  ---
!
      IF( ABS(DYMHGOX).GT.1.D-7 ) GOTO 10
!
!---  Converged solution  ---
!
      CALL HYD_P( TEQHX,YMHGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
      PVHX = PEQHX
      PVHOX = PVHX*YMHGOX
      PVHAX = PVHX-PVHOX
!
!---  Hydrate-mobile phase equilibrium conditions  ---
!
      PVOX = PVHOX
      PVAX = PVHAX
!
!---  Dissolved CH4 and CO2 concentrations  ---
!
      INCG = 1
      CALL SOL_LA( TX,PVAX,XLSX,XMLAX )
      INCG = 2
      CALL SOL_LO( TX,PVOX,XLSX,XMLOX )
      XMLSX = XLSX*(XMLOX*WTMO + XMLAX*WTMA + 
     &  (1.D+0-XMLOX-XMLAX)*WTMW)/((1.D+0-XLSX)*WTMS + XLSX*WTMW)
      XMLWX = MAX( 1.D+0-XMLSX-XMLAX-XMLOX,0.D+0 )
      WTMLX = XMLSX*WTMS + XMLAX*WTMA + XMLOX*WTMO + XMLWX*WTMW
      XLAX = XMLAX*WTMA/WTMLX
      XLOX = XMLOX*WTMO/WTMLX
!
!---  Vapor-pressure lowering  ---
!
      CALL SP_B( TX,XLSX,PSBX )
      PPLX = MAX( PLX+PATM,PSBX )
      CALL DENS_B( TX,PPLX,XLSX,RHOBX )
      IF( ISLC(44).EQ.1 ) THEN
        PVWX = PSBX
      ELSE
        CALL VPL_BX( TX,PLX,PVAX,PVOX,PVWX,PSBX,RHOBX,XLSX )
      ENDIF
      PGX = MAX( PVOX+PVAX+PVWX-PATM,PLX+(ENPR/BGL) )
!
!---  Limit minimum aqueous pressure  ---
!
      PLX = MAX( PLX,(PGX-HDOD*HCPBFX*RHORL*GRAV/BGL) )
!
!---  Liquid CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      PNX = PNC
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_IC5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_IC6( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for initial condition #6
!
!     SHX = 0.0
!     SNX > 0.0
!     SGX = 0.0
!
!     Declared variables:
!
!     TX - temperature, C
!     PLX - aqueous pressure, Pa
!     PNX - liquid-CO2 pressure, Pa
!     WLOX - relative saturation of dissolved CH4 (input as PVOX )
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 August 2010.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NCG_PT
      USE HYST
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
      REAL*8 RKLX(3)
      SAVE IALLOC,AJ,BJ,GX,IJ
      DATA IALLOC /0/
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: AJ,GX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: BJ
      INTEGER, DIMENSION(:), ALLOCATABLE :: IJ
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_IC6'
!
!---  Dynamic memory allocation  ---
!
      IF( IALLOC.EQ.0 ) THEN
        ALLOCATE( AJ(1:2,1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: AJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( BJ(1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: BJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( GX(1:2,1:3),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: GX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( IJ(1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: IJ'
          CALL WRMSGS( INDX )
        ENDIF
        IALLOC = 1
      ENDIF
!!
!!---  Assign mobile CH4 mole fraction of formers  ---
!!
!      YMMOX = MIN( MAX( PVOX,0.D+0 ),1.D+0 )
!
!---  Assign relative saturation of dissolved CH4  ---
!
      WLOX = MIN( MAX( PVOX,0.D+0 ),1.D+0 )
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Saturated CO2 conditions  ---
!      
      CALL SP_A( TX,PSAX )
      PVAX = PSAX
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure without hydrate barrier adjustment  ---
!
      ENPR = SCHR(2,IZN)*RHORL*GRAV
!
!---  Limit minimum aqueous pressure  ---
!
      PLX = MAX( PLX,(PNX-HDOD*RHORL*GRAV/BNL) )
!
!---  Gas pressure  ---
!
      PGX = MAX( PNX+ENPR/BGN,PLX+ENPR/BGL )
!
!---  Guess CO2 and CH4 gas partial pressure  ---
!
      CALL SP_B( TX,XLSX,PSBX )
      PPLX = MAX( PLX+PATM,PSBX )
      CALL DENS_B( TX,PPLX,XLSX,RHOBX )
      PCX = MAX( BNL*(PNX-PLX),0.D+0 )
      IF( ISLC(44).EQ.1 ) THEN
        PVWX = PSBX
      ELSE
        CALL VPL_B( TX,PSBX,PCX,RHOBX,PVWX,XLSX )
      ENDIF
!
!---  Vapor pressure limit  ---
!
      PVX = PGX + PATM - PVWX
      PVOX = WLOX*(PVX-PVAX)
!
!---  Check for gas phase conditions  ---
!
      IF( (PVAX+PVOX+PVWX-PATM).GT.PGX ) THEN
        INDX = 11
        RLMSG = (PVAX+PVOX+PVWX-PATM) - PGX
        CHMSG = 'FLH_IC6: Gas Phase Conditions: PV - PG = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Liquid CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
!
!---  Check for no liquid-CO2 phase conditions  ---
!
      IF( PNX.LE.PNC ) THEN
        INDX = 11
        RLMSG = PNX - PNC
        CHMSG = 'FLH_IC6: No Liquid-CO2: PN - PNC = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Hydrate equilibrium pressure  ---
!
      IF( (PVAX+PVOX).GT.EPSL ) THEN
        PVHOX = PVOX
        PVHAX = PVAX
        YMHGOX = PVHOX/(PVHAX+PVHOX)
        PEQHX = PVHOX+PVHAX
!
!---    Inhibitor correction  ---
!
        IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
          CALL HYD_PO( YMHGOX,POHX )
        ENDIF
        CALL HYD_T( TEQHX,YMHGOX,XLSX,PEQHX,POHX,RHOHX,XHAX,XHOX,XHWX )
      ENDIF
!
!---  Hydrate conditions  ---
!
      IF( TX.LE.TEQHX ) THEN
        INDX = 11
        RLMSG = TX - TEQHX
        CHMSG = 'FLH_IC6: Hydrate Conditions: T - TEQH = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Hydrate properties for no-hydrate conditions  ---
!
      TMHAX = 0.D+0
      TMHOX = 0.D+0
      YMHGOX = 0.D+0
      TEQHX = -TABS
      XHAX = 0.D+0
      XHOX = 0.D+0
      XHWX = 1.D+0
      RHOHX = 0.D+0
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_IC6 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_IC7( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for initial condition #7
!
!     SHX = 0.0
!     SNX > 0.0
!     SGX > 0.0
!
!     Declared variables:
!
!     TX - temperature, C
!     PGX - gas pressure, Pa
!     PLX - aqueous pressure, Pa
!     PNX - liquid-CO2 pressure, Pa
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 August 2010.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NCG_PT
      USE HYST
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
      REAL*8 RKLX(3)
      SAVE IALLOC,AJ,BJ,GX,IJ
      DATA IALLOC /0/
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: AJ,GX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: BJ
      INTEGER, DIMENSION(:), ALLOCATABLE :: IJ
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_IC7'
!
!---  Dynamic memory allocation  ---
!
      IF( IALLOC.EQ.0 ) THEN
        ALLOCATE( AJ(1:2,1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: AJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( BJ(1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: BJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( GX(1:2,1:3),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: GX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( IJ(1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: IJ'
          CALL WRMSGS( INDX )
        ENDIF
        IALLOC = 1
      ENDIF
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Saturated CO2 conditions  ---
!      
      CALL SP_A( TX,PSAX )
      PVAX = PSAX
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure without hydrate barrier adjustment  ---
!
      ENPR = SCHR(2,IZN)*RHORL*GRAV
!
!---  Limit minimum aqueous pressure  ---
!
      PLX = MAX( PLX,(PNX-HDOD*RHORL*GRAV/BNL) )
!
!---  CH4 gas partial pressure  ---
!
      CALL SP_B( TX,XLSX,PSBX )
      PPLX = MAX( PLX+PATM,PSBX )
      CALL DENS_B( TX,PPLX,XLSX,RHOBX )
      PCX = MAX( BNL*(PNX-PLX),0.D+0 )
      IF( ISLC(44).EQ.1 ) THEN
        PVWX = PSBX
      ELSE
        CALL VPL_B( TX,PSBX,PCX,RHOBX,PVWX,XLSX )
      ENDIF
      PVOX = PGX+PATM-PVWX-PVAX
!
!---  Check for gas pressure less than vapor pressure  ---
!
      IF( PVOX.LT.0.D+0 ) THEN
        INDX = 11
        RLMSG = PVOX
        CHMSG = 'FLH_IC7: Negative CH4 Partial Pressure: PVO = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Check for no gas phase conditions  ---
!
      IF( PGX.LE.(PNX+ENPR/BGN) ) THEN
        INDX = 11
        RLMSG = PGX - (PNX+ENPR/BGN)
        CHMSG = 'FLH_IC7: No Gas: PG - (PN + ENPR/BGN) = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Critical liquid-CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
!
!---  Check for no liquid-CO2 phase conditions  ---
!
      IF( PNX.LE.PNC ) THEN
        INDX = 11
        RLMSG = PNX - PNC
        CHMSG = 'FLH_IC7: No Liquid-CO2: PN - PNC = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Hydrate equilibrium pressure  ---
!
      IF( (PVAX+PVOX).GT.EPSL ) THEN
        PVHOX = PVOX
        PVHAX = PVAX
        YMHGOX = PVHOX/(PVHAX+PVHOX)
        PEQHX = PVHOX+PVHAX
!
!---    Inhibitor correction  ---
!
        IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
          CALL HYD_PO( YMHGOX,POHX )
        ENDIF
        CALL HYD_T( TEQHX,YMHGOX,XLSX,PEQHX,POHX,RHOHX,XHAX,XHOX,XHWX )
      ENDIF
!
!---  Hydrate conditions  ---
!
      IF( TX.LE.TEQHX ) THEN
        INDX = 11
        RLMSG = TX - TEQHX
        CHMSG = 'FLH_IC7: Hydrate Conditions: T - TEQH = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Hydrate properties for no-hydrate conditions  ---
!
      TMHAX = 0.D+0
      TMHOX = 0.D+0
      YMHGOX = 0.D+0
      TEQHX = -TABS
      XHAX = 0.D+0
      XHOX = 0.D+0
      XHWX = 1.D+0
      RHOHX = 0.D+0
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_IC7 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_IC8( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for initial condition #8
!
!     SHX > 0.0
!     SNX > 0.0
!     SGX = 0.0
!
!     Declared variables:
!
!     TX - temperature, C
!     PLX - aqueous pressure, Pa
!     PNX - liquid-CO2 pressure, Pa
!     SHX - hydrate saturation
!     WLOX - relative saturation of dissolved CH4 (input as PVOX )
!     YHMOX - hydrate CH4 mole fraction of formers (input as PVHOX )
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 August 2010.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NAPL
      USE HYST
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
      REAL*8 RKLX(3)
      SAVE IALLOC,AJ,BJ,GX,IJ
      DATA IALLOC /0/
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: AJ,GX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: BJ
      INTEGER, DIMENSION(:), ALLOCATABLE :: IJ
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_IC8'
!
!---  Dynamic memory allocation  ---
!
      IF( IALLOC.EQ.0 ) THEN
        ALLOCATE( AJ(1:2,1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: AJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( BJ(1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: BJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( GX(1:2,1:3),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: GX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( IJ(1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: IJ'
          CALL WRMSGS( INDX )
        ENDIF
        IALLOC = 1
      ENDIF
!
!---  Assign relative saturation of dissolved CH4 and
!     hydrate CH4 mole fraction of formers  ---
!
      WLOX = MIN( MAX( PVOX,0.D+0 ),1.D+0 )
      YMHOX = MIN( MAX( PVHOX,0.D+0 ),1.D+0 )
!
!---  Limit hydrate saturation  ---
!
      SHX = MIN( MAX( SHX,0.D+0 ),1.D+0 )
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Saturated CO2 conditions  ---
!      
      CALL SP_A( TX,PSAX )
      PVAX = PSAX
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure with hydrate barrier adjustment  ---
!
      SHZ = MAX( SHX,1.D-9 )
      HCPBFX = 1.D+0
      IF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.12 ) THEN
        HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &    SIN(SCHR(12,IZN)*SHZ + SCHR(13,IZN)) )
      ENDIF
      ENPR = SCHR(2,IZN)*RHORL*GRAV*HCPBFX
!
!---  Limit minimum aqueous pressure  ---
!
      PLX = MAX( PLX,(PNX-HDOD*HCPBFX*RHORL*GRAV/BNL) )
!
!---  Pore-space equilibrium temperature  ---
!
      TEQHX = TX
!
!---  Guess hydrate-gas CH4 mole fraction of formers  ---
!
      YMHGOX = YMHOX
!
!---  Single-variable Newton-Raphson loop (YMHGOX)  ---
!
      NC = 0
   10 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 17
        CHMSG = 'Unconverged FLH_IC8: Hydrate-Gas CH4 ' //
     &     'Mole Fraction of Formers'
        CALL WRMSGS( INDX )
      ENDIF
      DO 20 M = 1,2
        YMHGOY = YMHGOX
        DYMHGOX = SIGN( 1.D-6,(1.D-1-YMHGOX) )
        IF( M.EQ.2 ) YMHGOY = YMHGOX + DYMHGOX
        CALL HYD_P( TEQHX,YMHGOY,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
        XMHOX = XHOX/WTMO
        XMHAX = XHAX/WTMA
        GX(M,1) = YMHOX - XMHOX/(XMHOX+XMHAX)
   20 CONTINUE
!
!---  Linearization  ---
!
      DFX = (GX(2,1)-GX(1,1))/DYMHGOX
      FX = -GX(1,1)
!
!---  Update primary unknowns  ---
!
      DYMHGOX = FX/DFX
      YMHGOX = YMHGOX + DYMHGOX
!
!---  Limit hydrate-gas CH4 mole fraction of formers  ---
!
      YMHGOX = MIN( MAX( YMHGOX,0.D+0 ),1.D+0 )
!
!---  Convergence check  ---
!
      IF( ABS(DYMHGOX).GT.1.D-7 ) GOTO 10
!
!---  Converged solution  ---
!
      CALL HYD_P( TEQHX,YMHGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
      PVHX = PEQHX
      PVHOX = PVHX*YMHGOX
      PVHAX = PVHX-PVHOX
!
!---  Gas pressure  ---
!
      PGX = MAX( PNX+ENPR/BGN,PLX+ENPR/BGL )
!
!---  Guess CO2 and CH4 gas partial pressure  ---
!
      CALL SP_B( TX,XLSX,PSBX )
      PPLX = MAX( PLX+PATM,PSBX )
      CALL DENS_B( TX,PPLX,XLSX,RHOBX )
      PCX = MAX( BNL*(PNX-PLX),0.D+0 )
      IF( ISLC(44).EQ.1 ) THEN
        PVWX = PSBX
      ELSE
        CALL VPL_B( TX,PSBX,PCX,RHOBX,PVWX,XLSX )
      ENDIF
!
!---  Vapor pressure limit  ---
!
      PVX = PGX + PATM - PVWX
      PVOX = WLOX*(PVX-PVAX)
!
!---  Check for gas phase conditions  ---
!
      IF( (PVAX+PVOX+PVWX-PATM).GT.PGX ) THEN
        INDX = 11
        RLMSG = (PVAX+PVOX+PVWX-PATM) - PGX
        CHMSG = 'FLH_IC8: Gas Phase Conditions: PV - PG = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Liquid CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
!
!---  Check for no liquid-CO2 phase conditions  ---
!
      IF( PNX.LE.PNC ) THEN
        INDX = 11
        RLMSG = PNX - PNC
        CHMSG = 'FLH_IC8: No Liquid-CO2: PN - PNC = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_IC8 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_IC9( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for initial condition #7
!
!     SHX > 0.0
!     SNX > 0.0
!     SGX > 0.0
!
!     Declared variables:
!
!     TX - temperature, C
!     PGX - gas pressure, Pa
!     PLX - aqueous pressure, Pa
!     PNX - liquid-CO2 pressure, Pa
!     SHX - hydrate saturation
!     YMHOX - hydrate CH4 mole fraction of formers (input as PVHOX)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 August 2010.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NAPL
      USE HYST
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
      REAL*8 RKLX(3)
      SAVE IALLOC,AJ,BJ,GX,IJ
      DATA IALLOC /0/
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: AJ,GX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: BJ
      INTEGER, DIMENSION(:), ALLOCATABLE :: IJ
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_IC9'
!
!---  Dynamic memory allocation  ---
!
      IF( IALLOC.EQ.0 ) THEN
        ALLOCATE( AJ(1:2,1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: AJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( BJ(1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: BJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( GX(1:2,1:3),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: GX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( IJ(1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: IJ'
          CALL WRMSGS( INDX )
        ENDIF
        IALLOC = 1
      ENDIF
!
!---  Assign hydrate CH4 mole fraction of formers  ---
!
      YMHOX = MIN( MAX( PVHOX,0.D+0 ),1.D+0 )
!
!---  Limit hydrate saturation  ---
!
      SHX = MIN( MAX( SHX,0.D+0 ),1.D+0 )
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Saturated CO2 conditions  ---
!      
      CALL SP_A( TX,PSAX )
      PVAX = PSAX
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure with hydrate barrier adjustment  ---
!
      SHZ = MAX( SHX,1.D-9 )
      HCPBFX = 1.D+0
      IF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.12 ) THEN
        HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &    SIN(SCHR(12,IZN)*SHZ + SCHR(13,IZN)) )
      ENDIF
      ENPR = SCHR(2,IZN)*RHORL*GRAV*HCPBFX
!
!---  Limit minimum aqueous pressure  ---
!
      PLX = MAX( PLX,(PNX-HDOD*HCPBFX*RHORL*GRAV/BNL) )
!
!---  Pore-space equilibrium temperature  ---
!
      TEQHX = TX
!
!---  Guess hydrate-gas CH4 mole fraction of formers  ---
!
      YMHGOX = YMHOX
!
!---  Single-variable Newton-Raphson loop (YMHGOX)  ---
!
      NC = 0
   10 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 17
        CHMSG = 'Unconverged FLH_IC9: Hydrate-Gas CH4 ' //
     &     'Mole Fraction of Formers'
        CALL WRMSGS( INDX )
      ENDIF
      DO 20 M = 1,2
        YMHGOY = YMHGOX
        DYMHGOX = SIGN( 1.D-6,(1.D-1-YMHGOX) )
        IF( M.EQ.2 ) YMHGOY = YMHGOX + DYMHGOX
        CALL HYD_P( TEQHX,YMHGOY,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
        XMHOX = XHOX/WTMO
        XMHAX = XHAX/WTMA
        GX(M,1) = YMHOX - XMHOX/(XMHOX+XMHAX)
   20 CONTINUE
!
!---  Linearization  ---
!
      DFX = (GX(2,1)-GX(1,1))/DYMHGOX
      FX = -GX(1,1)
!
!---  Update primary unknowns  ---
!
      DYMHGOX = FX/DFX
      YMHGOX = YMHGOX + DYMHGOX
!
!---  Limit hydrate-gas CH4 mole fraction of formers  ---
!
      YMHGOX = MIN( MAX( YMHGOX,0.D+0 ),1.D+0 )
!
!---  Convergence check  ---
!
      IF( ABS(DYMHGOX).GT.1.D-7 ) GOTO 10
!
!---  Converged solution  ---
!
      CALL HYD_P( TEQHX,YMHGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
      PVHX = PEQHX
      PVHOX = PVHX*YMHGOX
      PVHAX = PVHX-PVHOX
!
!---  CH4 gas partial pressure  ---
!
      CALL SP_B( TX,XLSX,PSBX )
      PPLX = MAX( PLX+PATM,PSBX )
      CALL DENS_B( TX,PPLX,XLSX,RHOBX )
      PCX = MAX( BNL*(PNX-PLX),0.D+0 )
      IF( ISLC(44).EQ.1 ) THEN
        PVWX = PSBX
      ELSE
        CALL VPL_B( TX,PSBX,PCX,RHOBX,PVWX,XLSX )
      ENDIF
      PVOX = PGX+PATM-PVWX-PVAX
!
!---  Check for gas pressure less than vapor pressure  ---
!
      IF( PVOX.LT.0.D+0 ) THEN
        INDX = 11
        RLMSG = PVOX
        CHMSG = 'FLH_IC9: Negative CH4 Partial Pressure: PVO = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Check for no gas phase conditions  ---
!
      IF( PGX.LE.(PNX+ENPR/BGN) ) THEN
        INDX = 11
        RLMSG = PGX - (PNX+ENPR/BGN)
        CHMSG = 'FLH_IC9: No Gas: PG - (PN + ENPR/BGN) = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Critical liquid-CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
!
!---  Check for no liquid-CO2 phase conditions  ---
!
      IF( PNX.LE.PNC ) THEN
        INDX = 11
        RLMSG = PNX - PNC
        CHMSG = 'FLH_IC9: No Liquid-CO2: PN - PNC = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_IC9 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_IC10( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for initial condition #10
!
!     SHX > 0.0
!     SNX > 0.0
!     SGX ³ 0.0
!     Hydrate-mobile phase equilibrium
!
!     Declared variables:
!
!     TX - temperature, C
!     PLX - aqueous pressure, Pa
!     PNX - liquid-CO2 pressure, Pa
!     YMHOX - hydrate CH4 mole fraction of formers (input as PVHOX)
!     SHX - hydrate saturation
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 August 2010.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NCG_PT
      USE NAPL
      USE HYST
      USE FDVP
      USE FDVH
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLX(3)
      SAVE IALLOC,AJ,BJ,GX,IJ
      DATA IALLOC /0/
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: AJ,GX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: BJ
      INTEGER, DIMENSION(:), ALLOCATABLE :: IJ
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_IC10'
!
!---  Dynamic memory allocation  ---
!
      IF( IALLOC.EQ.0 ) THEN
        ALLOCATE( AJ(1:2,1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: AJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( BJ(1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: BJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( GX(1:2,1:3),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: GX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( IJ(1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: IJ'
          CALL WRMSGS( INDX )
        ENDIF
        IALLOC = 1
      ENDIF
!
!---  Assign hydrate CH4 mole fraction of formers  ---
!
      YMHOX = MIN( MAX( PVHOX,0.D+0 ),1.D+0 )
!
!---  Limit hydrate saturation  ---
!
      SHX = MIN( MAX( SHX,0.D+0 ),1.D+0 )
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Saturated CO2 conditions  ---
!      
      CALL SP_A( TX,PSAX )
      PVAX = PSAX
      PVHAX = PVAX
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure with hydrate barrier adjustment  ---
!
      SHZ = MAX( SHX,1.D-9 )
      HCPBFX = 1.D+0
      IF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.12 ) THEN
        HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &    SIN(SCHR(12,IZN)*SHZ + SCHR(13,IZN)) )
      ENDIF
      ENPR = SCHR(2,IZN)*RHORL*GRAV*HCPBFX
!
!---  Limit minimum aqueous pressure  ---
!
      PLX = MAX( PLX,(PNX-HDOD*HCPBFX*RHORL*GRAV/BNL) )
!
!---  Pore-space equilibrium temperature  ---
!
      TEQHX = TX
!
!---  Equilibrium pressure range  ---
!
      YMHGOX = 0.D+0
      CALL HYD_P( TEQHX,YMHGOX,XLSX,PEQH1X,RHOHX,XHAX,XHOX,XHWX )
!
!---  Equilibrium pressure range  ---
!
      YMHGOX = 1.D+0
      CALL HYD_P( TEQHX,YMHGOX,XLSX,PEQH2X,RHOHX,XHAX,XHOX,XHWX )
!
!---  Check equilibrium condition   ---
!
      IF( PVAX.GT.PEQH1X ) THEN
        INDX = 11
        RLMSG = PVAX - PEQH1X 
        CHMSG = 'FLH_IC10: No Equilibrium Condition: ' // 
     &    'PVA > PEQH (YMHGO=0): PVA - PEQH (YMHGO=0) = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Guess YMHGO  ---
!
      YMHGOX = 5.D-1
!
!---  Single-variable Newton-Raphson loop (YMHGO)  ---
!
      NC = 0
   10 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 17
        CHMSG = 'Unconverged FLH_IC10: Hydrate-Gas CH4 ' //
     &     'Mole Fraction of Formers'
        CALL WRMSGS( INDX )
      ENDIF
      DO 20 M = 1,2
        YMHGOY = YMHGOX
        DYMHGOX = SIGN( 1.D-6,(1.D-1-YMHGOX) )
        IF( M.EQ.2 ) YMHGOY = YMHGOX + DYMHGOX
        CALL HYD_P( TEQHX,YMHGOY,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
        GX(M,1) = PVHAX - (1.D+0-YMHGOY)*PEQHX
   20 CONTINUE
!
!---  Linearization  ---
!
      DFX = (GX(2,1)-GX(1,1))/DYMHGOX
      FX = -GX(1,1)
!
!---  Update primary unknowns  ---
!
      DYMHGOX = FX/DFX
      YMHGOX = YMHGOX + DYMHGOX
!
!---  Limit hydrate-gas CH4 mole fraction of formers  ---
!
      YMHGOX = MIN( MAX( YMHGOX,0.D+0 ),1.D+0 )
!
!---  Convergence check  ---
!
      IF( ABS(DYMHGOX).GT.1.D-7 ) GOTO 10
!
!---  Converged solution  ---
!
      CALL HYD_P( TEQHX,YMHGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
      PVHX = PEQHX
      PVHOX = PEQHX - PVHAX
      YMHOX = PVHOX/(PVHAX+PVHOX)
!
!---  Hydrate-mobile phase equilibrium conditions  ---
!
      PVOX = PVHOX
!
!---  Dissolved CH4 and CO2 concentrations  ---
!
      INCG = 1
      CALL SOL_LA( TX,PVAX,XLSX,XMLAX )
      INCG = 2
      CALL SOL_LO( TX,PVOX,XLSX,XMLOX )
      XMLSX = XLSX*(XMLOX*WTMO + XMLAX*WTMA + 
     &  (1.D+0-XMLOX-XMLAX)*WTMW)/((1.D+0-XLSX)*WTMS + XLSX*WTMW)
      XMLWX = MAX( 1.D+0-XMLSX-XMLAX-XMLOX,0.D+0 )
      WTMLX = XMLSX*WTMS + XMLAX*WTMA + XMLOX*WTMO + XMLWX*WTMW
      XLAX = XMLAX*WTMA/WTMLX
      XLOX = XMLOX*WTMO/WTMLX
!
!---  Vapor-pressure lowering  ---
!
      CALL SP_B( TX,XLSX,PSBX )
      PPLX = MAX( PLX+PATM,PSBX )
      CALL DENS_B( TX,PPLX,XLSX,RHOBX )
      PCX = MAX( BNL*(PNX-PLX),0.D+0 )
      IF( ISLC(44).EQ.1 ) THEN
        PVWX = PSBX
      ELSE
        CALL VPL_B( TX,PSBX,PCX,RHOBX,PVWX,XLSX )
      ENDIF
!
!---  Gas pressure  ---
!
      PGX = MAX( PVOX+PVAX+PVWX-PATM,PLX+(ENPR/BGL) )
!
!---  Liquid CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
!
!---  Check for no liquid-CO2 phase conditions  ---
!
      IF( PNX.LE.PNC ) THEN
        INDX = 11
        RLMSG = PNX - PNC
        CHMSG = 'FLH_IC10: No Liquid-CO2: PN - PNC = '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_IC10 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_PC1( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for phase condition #1
!
!     Aqueous, no hydrate, no liquid CO2
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     mobile CO2 mass - CO2 vapor pressure
!     mobile CH4 mass - CH4 vapor pressure
!     hydrate CO2 mass - hydrate CO2 mass
!     hydrate CH4 mass - hydrate CH4 mass
!     NaCl mass - total NaCl aqueous mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 October 2007.
!     Last Modified by MD White, PNNL, 12 October 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NCG_PT
      USE HYST
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_PC1'
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure without hydrate barrier adjustment  ---
!
      ENPR = SCHR(2,IZN)*RHORL*GRAV
!
!---  Brine vapor pressure  ---
!
      IF( PSBX.GT.0.D+0 ) THEN
        CALL SP_B( TX,XLSX,PSBX )
        PPLX = MAX( PLX+PATM,PSBX )
        CALL DENS_B( TX,PPLX,XLSX,RHOBX )
        IF( ISLC(44).EQ.1 ) THEN
          PVWX = PSBX
        ELSE
          CALL VPL_BX( TX,PLX,PVAX,PVOX,PVWX,PSBX,RHOBX,XLSX )
        ENDIF
      ENDIF
      PSBX = ABS(PSBX)
!
!---  Gas pressure  ---
!
!      PGX = MAX( PVAX+PVOX+PVWX-PATM,PLX+ENPR/BGL )
      PGX = PLX+ENPR/BGL
!
!---  Hydrate equilibrium temperature for no-hydrate conditions  ---
!
      TEQHX = -TABS
      YMHGOX = 0.D+0
      XHAX = 0.D+0
      XHOX = 0.D+0
      XHWX = 1.D+0
      RHOHX = 0.D+0
!
!---  Hydrate equilibrium temperature at the vapor pressure conditions  ---
!
      PVHOX = PVOX
      PVHAX = PVAX
      IF( (PVAX+PVOX).GT.EPSL ) THEN
        YMHGOX = PVHOX/(PVHAX+PVHOX)
        PEQHX = PVHOX+PVHAX
!
!---    Inhibitor correction  ---
!
        IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
          CALL HYD_PO( YMHGOX,POHX )
        ENDIF
        CALL HYD_T( TEQHX,YMHGOX,XLSX,PEQHX,POHX,RHOHX,XHAX,XHOX,XHWX )
!
!---    Hydrate equilibrium pressure at the actual temperature  ---
!
        IF( TX.LE.TEQHX ) THEN
          CALL HYD_P( TX,YMHGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
          PVHOX = YMHGOX*PEQHX
          PVHAX = PEQHX - PVHOX
!
!---    Hydrate equilibrium temperature for no-hydrate conditions  ---
!
        ELSE
          TEQHX = -TABS
          YMHGOX = 0.D+0
          XHAX = 0.D+0
          XHOX = 0.D+0
          XHWX = 1.D+0
          RHOHX = 0.D+0
        ENDIF
      ENDIF
!
!---  Liquid CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      PNX = PNC
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_PC1 group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_PC2( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for phase condition #2
!
!     Aqueous, hydrate, no liquid CO2, PVHO > PVHA
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     Mobile CO2 mass - CO2 vapor pressure
!     Mobile CH4 mass - CH4 vapor pressure
!     Hydrate CO2 mass - hydrate CH4 mole fraction of formers
!     Hydrate CH4 mass - hydrate saturation
!     NaCl mass - total NaCl aqueous mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 October 2007.
!     Last Modified by MD White, PNNL, 12 October 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE HYST
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_PC2'
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure with hydrate barrier adjustment  ---
!
      HCPBFX = 1.D+0
      IF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.12 ) THEN
        HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &    SIN(SCHR(12,IZN)*SHX + SCHR(13,IZN)) )
      ENDIF
      ENPR = SCHR(2,IZN)*RHORL*GRAV*HCPBFX
!
!---  Brine vapor pressure  ---
!
      IF( PSBX.GT.0.D+0 ) THEN
        CALL SP_B( TX,XLSX,PSBX )
        PPLX = MAX( PLX+PATM,PSBX )
        CALL DENS_B( TX,PPLX,XLSX,RHOBX )
        IF( ISLC(44).EQ.1 ) THEN
          PVWX = PSBX
        ELSE
          CALL VPL_BX( TX,PLX,PVAX,PVOX,PVWX,PSBX,RHOBX,XLSX )
        ENDIF
      ENDIF
      PSBX = ABS(PSBX)
!
!---  Gas pressure  ---
!
!      PGX = MAX( PVAX+PVOX+PVWX-PATM,PLX+ENPR/BGL )
      PGX = PLX+ENPR/BGL
!
!---  Pore-space equilibrium temperature  ---
!
      TEQHX = TX
!
!---  Hydrate equilibrium pressure  ---
!
      CALL HYD_P( TEQHX,YMHGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
      PVHOX = YMHGOX*PEQHX
      PVHAX = PEQHX - PVHOX
!
!---  Liquid CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      PNX = PNC
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Hydrate guest-molecule mass  ---
!
      TMHAX = XHAX*RHOHX*SHX
      TMHOX = XHOX*RHOHX*SHX
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_PC2 group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_PC3( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for phase condition #3
!
!     Aqueous, hydrate, no liquid CO2, PVHO < PVHA
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     Mobile CO2 mass - CO2 vapor pressure
!     Mobile CH4 mass - CH4 vapor pressure
!     Hydrate CO2 mass - hydrate saturation
!     Hydrate CH4 mass - hydrate CH4 mole fraction of formers
!     NaCl mass - total NaCl aqueous mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 October 2007.
!     Last Modified by MD White, PNNL, 12 October 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE HYST
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_PC3'
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure with hydrate barrier adjustment  ---
!
      HCPBFX = 1.D+0
      IF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.12 ) THEN
        HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &    SIN(SCHR(12,IZN)*SHX + SCHR(13,IZN)) )
      ENDIF
      ENPR = SCHR(2,IZN)*RHORL*GRAV*HCPBFX
!
!---  Brine vapor pressure  ---
!
      IF( PSBX.GT.0.D+0 ) THEN
        CALL SP_B( TX,XLSX,PSBX )
        PPLX = MAX( PLX+PATM,PSBX )
        CALL DENS_B( TX,PPLX,XLSX,RHOBX )
        IF( ISLC(44).EQ.1 ) THEN
          PVWX = PSBX
        ELSE
          CALL VPL_BX( TX,PLX,PVAX,PVOX,PVWX,PSBX,RHOBX,XLSX )
        ENDIF
      ENDIF
      PSBX = ABS(PSBX)
!
!---  Gas pressure  ---
!
!      PGX = MAX( PVAX+PVOX+PVWX-PATM,PLX+ENPR/BGL )
      PGX = PLX+ENPR/BGL
!
!---  Pore-space equilibrium temperature  ---
!
      TEQHX = TX
!
!---  Hydrate equilibrium pressure  ---
!
      CALL HYD_P( TEQHX,YMHGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
      PVHOX = YMHGOX*PEQHX
      PVHAX = PEQHX - PVHOX
!
!---  Liquid CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      PNX = PNC
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Hydrate guest-molecule mass  ---
!
      TMHAX = XHAX*RHOHX*SHX
      TMHOX = XHOX*RHOHX*SHX
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_PC3 group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_PC4( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for phase condition #4
!
!     Aqueous, no hydrate, liquid CO2
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     mobile CO2 mass - liquid-CO2 pressure
!     mobile CH4 mass - CH4 vapor pressure
!     hydrate CO2 mass - hydrate CO2 mass
!     hydrate CH4 mass - hydrate CH4 mass
!     NaCl mass - total NaCl aqueous mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 October 2007.
!     Last Modified by MD White, PNNL, 12 October 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NCG_PT
      USE HYST
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_PC4'
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure without hydrate barrier adjustment  ---
!
      ENPR = SCHR(2,IZN)*RHORL*GRAV
!
!---  Maximum capillary pressure  ---
!
      PLX = MAX( PNX-(HDOD*RHORL*GRAV/BNL),PLX )
!
!---  Brine vapor pressure  ---
!
      IF( PSBX.GT.0.D+0 ) THEN
        CALL SP_B( TX,XLSX,PSBX )
        PPLX = MAX( PLX+PATM,PSBX )
        CALL DENS_B( TX,PPLX,XLSX,RHOBX )
        PCX = MAX( BNL*(PNX-PLX),0.D+0 )
        IF( ISLC(44).EQ.1 ) THEN
          PVWX = PSBX
        ELSE
          CALL VPL_B( TX,PSBX,PCX,RHOBX,PVWX,XLSX )
        ENDIF
      ENDIF
      PSBX = ABS(PSBX)
!
!---  Gas pressure  ---
!
!      PGX = MAX( PVAX+PVOX+PVWX-PATM,PLX+ENPR/BGL,PNX+ENPR/BGN )
      PGX = MAX( PLX+ENPR/BGL,PNX+ENPR/BGN )
!
!---  Hydrate equilibrium temperature for no-hydrate conditions  ---
!
      TEQHX = -TABS
      YMHGOX = 0.D+0
      XHAX = 0.D+0
      XHOX = 0.D+0
      XHWX = 1.D+0
      RHOHX = 0.D+0
!
!---  Hydrate equilibrium temperature at the vapor pressure conditions  ---
!
      PVHOX = PVOX
      PVHAX = PVAX
      IF( (PVAX+PVOX).GT.EPSL ) THEN
        YMHGOX = PVHOX/(PVHAX+PVHOX)
        PEQHX = PVHOX+PVHAX
!
!---    Inhibitor correction  ---
!
        IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
          CALL HYD_PO( YMHGOX,POHX )
        ENDIF
        CALL HYD_T( TEQHX,YMHGOX,XLSX,PEQHX,POHX,RHOHX,XHAX,XHOX,XHWX )
!
!---    Hydrate equilibrium pressure at the actual temperature  ---
!
        IF( TX.LE.TEQHX ) THEN
          CALL HYD_P( TX,YMHGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
          PVHOX = YMHGOX*PEQHX
          PVHAX = PEQHX - PVHOX
!
!---    Hydrate equilibrium temperature for no-hydrate conditions  ---
!
        ELSE
          TEQHX = -TABS
          YMHGOX = 0.D+0
          XHAX = 0.D+0
          XHOX = 0.D+0
          XHWX = 1.D+0
          RHOHX = 0.D+0
        ENDIF
      ENDIF
!
!---  Liquid CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      PNX = MAX( PNX,PNC )
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_PC4 group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_PC5( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for phase condition #5
!
!     Aqueous, hydrate, liquid CO2, PVHO > PVHA
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     Mobile CO2 mass - liquid-CO2 pressure
!     Mobile CH4 mass - CH4 vapor pressure
!     Hydrate CO2 mass - hydrate CH4 mole fraction of formers
!     Hydrate CH4 mass - hydrate saturation
!     NaCl mass - total NaCl aqueous mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 October 2007.
!     Last Modified by MD White, PNNL, 12 October 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE HYST
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_PC5'
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure with hydrate barrier adjustment  ---
!
      HCPBFX = 1.D+0
      IF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.12 ) THEN
        HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &    SIN(SCHR(12,IZN)*SHX + SCHR(13,IZN)) )
      ENDIF
      ENPR = SCHR(2,IZN)*RHORL*GRAV*HCPBFX
!
!---  Maximum capillary pressure  ---
!
      PLX = MAX( PNX-(HDOD*RHORL*GRAV/BNL),PLX )
!
!---  Brine vapor pressure  ---
!
      IF( PSBX.GT.0.D+0 ) THEN
        CALL SP_B( TX,XLSX,PSBX )
        PPLX = MAX( PLX+PATM,PSBX )
        CALL DENS_B( TX,PPLX,XLSX,RHOBX )
        PCX = MAX( BNL*(PNX-PLX),0.D+0 )
        IF( ISLC(44).EQ.1 ) THEN
          PVWX = PSBX
        ELSE
          CALL VPL_B( TX,PSBX,PCX,RHOBX,PVWX,XLSX )
        ENDIF
      ENDIF
      PSBX = ABS(PSBX)
!
!---  Gas pressure  ---
!
!      PGX = MAX( PVAX+PVOX+PVWX-PATM,PLX+ENPR/BGL,PNX+ENPR/BGN )
      PGX = MAX( PLX+ENPR/BGL,PNX+ENPR/BGN )
!
!---  Pore-space equilibrium temperature  ---
!
      TEQHX = TX
!
!---  Hydrate equilibrium pressure  ---
!
      CALL HYD_P( TEQHX,YMHGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
      PVHOX = YMHGOX*PEQHX
      PVHAX = PEQHX - PVHOX
!
!---  Liquid CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      PNX = MAX( PNX,PNC )
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Hydrate guest-molecule mass  ---
!
      TMHAX = XHAX*RHOHX*SHX
      TMHOX = XHOX*RHOHX*SHX
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_PC5 group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_PC6( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for phase condition #6
!
!     Aqueous, hydrate, liquid CO2, PVHO < PVHA
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     Mobile CO2 mass - liquid-CO2 pressure
!     Mobile CH4 mass - CH4 vapor pressure
!     Hydrate CO2 mass - hydrate saturation
!     Hydrate CH4 mass - hydrate CH4 mole fraction of formers
!     NaCl mass - total NaCl aqueous mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 October 2007.
!     Last Modified by MD White, PNNL, 12 October 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE HYST
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_PC6'
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure with hydrate barrier adjustment  ---
!
      HCPBFX = 1.D+0
      IF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.12 ) THEN
        HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &    SIN(SCHR(12,IZN)*SHX + SCHR(13,IZN)) )
      ENDIF
      ENPR = SCHR(2,IZN)*RHORL*GRAV*HCPBFX
!
!---  Maximum capillary pressure  ---
!
      PLX = MAX( PNX-(HDOD*RHORL*GRAV/BNL),PLX )
!
!---  Brine vapor pressure  ---
!
      IF( PSBX.GT.0.D+0 ) THEN
        CALL SP_B( TX,XLSX,PSBX )
        PPLX = MAX( PLX+PATM,PSBX )
        CALL DENS_B( TX,PPLX,XLSX,RHOBX )
        PCX = MAX( BNL*(PNX-PLX),0.D+0 )
        IF( ISLC(44).EQ.1 ) THEN
          PVWX = PSBX
        ELSE
          CALL VPL_B( TX,PSBX,PCX,RHOBX,PVWX,XLSX )
        ENDIF
      ENDIF
      PSBX = ABS(PSBX)
!
!---  Gas pressure  ---
!
!      PGX = MAX( PVAX+PVOX+PVWX-PATM,PLX+ENPR/BGL,PNX+ENPR/BGN )
      PGX = MAX( PLX+ENPR/BGL,PNX+ENPR/BGN )
!
!---  Pore-space equilibrium temperature  ---
!
      TEQHX = TX
!
!---  Hydrate equilibrium pressure  ---
!
      CALL HYD_P( TEQHX,YMHGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
      PVHOX = YMHGOX*PEQHX
      PVHAX = PEQHX - PVHOX
!
!---  Liquid CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      PNX = MAX( PNX,PNC )
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Hydrate guest-molecule mass  ---
!
      TMHAX = XHAX*RHOHX*SHX
      TMHOX = XHOX*RHOHX*SHX
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_PC6 group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_PC7( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for phase condition #11
!
!     Aqueous-gas, no hydrate, no liquid CO2
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     mobile CO2 mass - gas pressure
!     mobile CH4 mass - CH4 vapor pressure
!     hydrate CO2 mass - hydrate CO2 mass
!     hydrate CH4 mass - hydrate CH4 mass
!     NaCl mass - total NaCl aqueous mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 October 2007.
!     Last Modified by MD White, PNNL, 12 October 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NCG_PT
      USE HYST
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_PC7'
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure without hydrate barrier adjustment  ---
!
      ENPR = SCHR(2,IZN)*RHORL*GRAV
!
!---  Maximum capillary pressure  ---
!
      PLX = MAX( PGX-(HDOD*RHORL*GRAV/BGL),PLX )
!
!---  Brine vapor pressure  ---
!
      IF( PSBX.GT.0.D+0 ) THEN
        CALL SP_B( TX,XLSX,PSBX )
        PPLX = MAX( PLX+PATM,PSBX )
        CALL DENS_B( TX,PPLX,XLSX,RHOBX )
        PCX = MAX( BGL*(PGX-PLX),0.D+0 )
        IF( ISLC(44).EQ.1 ) THEN
          PVWX = PSBX
        ELSE
          CALL VPL_B( TX,PSBX,PCX,RHOBX,PVWX,XLSX )
        ENDIF
      ENDIF
      PSBX = ABS(PSBX)
      PVX = PGX+PATM-PVWX
!
!---  PVA > PVO  ---
!     Mobile CO2 mass - gas pressure
!     Mobile CH4 mass - CH4 vapor pressure
!
      IF( NPHAZ(2,N).LT.0 ) THEN
        PVOX = MIN( PVX,PVOX )
        PVAX = MAX( PVX-PVOX,0.D+0 )
!
!---  PVA ² PVO  ---
!     Mobile CO2 mass - CO2 vapor pressure
!     Mobile CH4 mass - gas pressure
!
      ELSE
        PVAX = MIN( PVX,PVAX )
        PVOX = MAX( PVX-PVAX,0.D+0 )
      ENDIF
!
!---  Hydrate equilibrium temperature for no-hydrate conditions  ---
!
      TEQHX = -TABS
      YMHGOX = 0.D+0
      XHAX = 0.D+0
      XHOX = 0.D+0
      XHWX = 1.D+0
      RHOHX = 0.D+0
!
!---  Hydrate equilibrium temperature at the vapor pressure conditions  ---
!
      PVHOX = PVOX
      PVHAX = PVAX
      IF( (PVAX+PVOX).GT.EPSL ) THEN
        YMHGOX = PVHOX/(PVHAX+PVHOX)
        PEQHX = PVHOX+PVHAX
!
!---    Inhibitor correction  ---
!
        IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
          CALL HYD_PO( YMHGOX,POHX )
        ENDIF
        CALL HYD_T( TEQHX,YMHGOX,XLSX,PEQHX,POHX,RHOHX,XHAX,XHOX,XHWX )
!
!---    Hydrate equilibrium pressure at the actual temperature  ---
!
        IF( TX.LE.TEQHX ) THEN
          CALL HYD_P( TX,YMHGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
          PVHOX = YMHGOX*PEQHX
          PVHAX = PEQHX - PVHOX
!
!---    Hydrate equilibrium temperature for no-hydrate conditions  ---
!
        ELSE
          TEQHX = -TABS
          YMHGOX = 0.D+0
          XHAX = 0.D+0
          XHOX = 0.D+0
          XHWX = 1.D+0
          RHOHX = 0.D+0
        ENDIF
      ENDIF
!
!---  Liquid CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      PNX = PNC
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_PC7 group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_PC8( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for phase condition #12
!
!     Aqueous-gas, hydrate, no liquid CO2, PVHO > PVHA
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     Mobile CO2 mass - gas pressure
!     Mobile CH4 mass - CH4 vapor pressure
!     Hydrate CO2 mass - hydrate CH4 mole fraction of formers
!     Hydrate CH4 mass - hydrate saturation
!     NaCl mass - total NaCl aqueous mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 October 2007.
!     Last Modified by MD White, PNNL, 12 October 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE HYST
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_PC8'
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure with hydrate barrier adjustment  ---
!
      HCPBFX = 1.D+0
      IF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.12 ) THEN
        HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &    SIN(SCHR(12,IZN)*SHX + SCHR(13,IZN)) )
      ENDIF
      ENPR = SCHR(2,IZN)*RHORL*GRAV*HCPBFX
!
!---  Maximum capillary pressure  ---
!
      PLX = MAX( PGX-(HDOD*RHORL*GRAV/BGL),PLX )
!
!---  Pore-space equilibrium temperature  ---
!
      TEQHX = TX
!
!---  Hydrate equilibrium pressure  ---
!
      CALL HYD_P( TEQHX,YMHGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
      PVHOX = YMHGOX*PEQHX
      PVHAX = PEQHX - PVHOX
!
!---  Brine vapor pressure  ---
!
      IF( PSBX.GT.0.D+0 ) THEN
        CALL SP_B( TX,XLSX,PSBX )
        PPLX = MAX( PLX+PATM,PSBX )
        CALL DENS_B( TX,PPLX,XLSX,RHOBX )
        PCX = MAX( BGL*(PGX-PLX),0.D+0 )
        IF( ISLC(44).EQ.1 ) THEN
          PVWX = PSBX
        ELSE
          CALL VPL_B( TX,PSBX,PCX,RHOBX,PVWX,XLSX )
        ENDIF
      ENDIF
      PSBX = ABS(PSBX)
      PVX = PGX+PATM-PVWX
!
!---  PVA > PVO  ---
!     Mobile CO2 mass - gas pressure
!     Mobile CH4 mass - CH4 vapor pressure
!
      IF( NPHAZ(2,N).LT.0 ) THEN
        PVOX = MIN( PVX,PVOX )
        PVAX = MAX( PVX-PVOX,0.D+0 )
!
!---  PVA ² PVO  ---
!     Mobile CO2 mass - CO2 vapor pressure
!     Mobile CH4 mass - gas pressure
!
      ELSE
        PVAX = MIN( PVX,PVAX )
        PVOX = MAX( PVX-PVAX,0.D+0 )
      ENDIF
!
!---  Liquid CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      PNX = PNC
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Hydrate guest-molecule mass  ---
!
      TMHAX = XHAX*RHOHX*SHX
      TMHOX = XHOX*RHOHX*SHX
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_PC8 group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_PC9( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for phase condition #13
!
!     Aqueous-gas, hydrate, no liquid CO2, PVHO < PVHA
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     Mobile CO2 mass - gas pressure
!     Mobile CH4 mass - CH4 vapor pressure
!     Hydrate CO2 mass - hydrate saturation
!     Hydrate CH4 mass - hydrate CH4 mole fraction of formers
!     NaCl mass - total NaCl aqueous mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 October 2007.
!     Last Modified by MD White, PNNL, 12 October 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE HYST
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_PC9'
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure with hydrate barrier adjustment  ---
!
      HCPBFX = 1.D+0
      IF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.12 ) THEN
        HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &    SIN(SCHR(12,IZN)*SHX + SCHR(13,IZN)) )
      ENDIF
      ENPR = SCHR(2,IZN)*RHORL*GRAV*HCPBFX
!
!---  Maximum capillary pressure  ---
!
      PLX = MAX( PGX-(HDOD*RHORL*GRAV/BGL),PLX )
!
!---  Pore-space equilibrium temperature  ---
!
      TEQHX = TX
!
!---  Hydrate equilibrium pressure  ---
!
      CALL HYD_P( TEQHX,YMHGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
      PVHOX = YMHGOX*PEQHX
      PVHAX = PEQHX - PVHOX
!
!---  Brine vapor pressure  ---
!
      IF( PSBX.GT.0.D+0 ) THEN
        CALL SP_B( TX,XLSX,PSBX )
        PPLX = MAX( PLX+PATM,PSBX )
        CALL DENS_B( TX,PPLX,XLSX,RHOBX )
        PCX = MAX( BGL*(PGX-PLX),0.D+0 )
        IF( ISLC(44).EQ.1 ) THEN
          PVWX = PSBX
        ELSE
          CALL VPL_B( TX,PSBX,PCX,RHOBX,PVWX,XLSX )
        ENDIF
      ENDIF
      PSBX = ABS(PSBX)
      PVX = PGX+PATM-PVWX
!
!---  PVA > PVO  ---
!     Mobile CO2 mass - gas pressure
!     Mobile CH4 mass - CH4 vapor pressure
!
      IF( NPHAZ(2,N).LT.0 ) THEN
        PVOX = MIN( PVX,PVOX )
        PVAX = MAX( PVX-PVOX,0.D+0 )
!
!---  PVA ² PVO  ---
!     Mobile CO2 mass - CO2 vapor pressure
!     Mobile CH4 mass - gas pressure
!
      ELSE
        PVAX = MIN( PVX,PVAX )
        PVOX = MAX( PVX-PVAX,0.D+0 )
      ENDIF
!
!---  Liquid CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      PNX = PNC
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Hydrate guest-molecule mass  ---
!
      TMHAX = XHAX*RHOHX*SHX
      TMHOX = XHOX*RHOHX*SHX
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_PC9 group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_PC10( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for phase condition #14
!
!     Aqueous-gas, no hydrate, liquid CO2
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     mobile CO2 mass - liquid-CO2 pressure
!     mobile CH4 mass - gas pressure
!     hydrate CO2 mass - hydrate CO2 mass
!     hydrate CH4 mass - hydrate CH4 mass
!     NaCl mass - total NaCl aqueous mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 October 2007.
!     Last Modified by MD White, PNNL, 12 October 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NCG_PT
      USE HYST
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_PC10'
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure without hydrate barrier adjustment  ---
!
      ENPR = SCHR(2,IZN)*RHORL*GRAV
!
!---  Maximum capillary pressure  ---
!
      PLX = MAX( PNX-(HDOD*RHORL*GRAV/BNL),PLX )
!
!---  Brine vapor pressure  ---
!
      IF( PSBX.GT.0.D+0 ) THEN
        CALL SP_B( TX,XLSX,PSBX )
        PPLX = MAX( PLX+PATM,PSBX )
        CALL DENS_B( TX,PPLX,XLSX,RHOBX )
        PCX = MAX( BNL*(PNX-PLX),0.D+0 )
        IF( ISLC(44).EQ.1 ) THEN
          PVWX = PSBX
        ELSE
          CALL VPL_B( TX,PSBX,PCX,RHOBX,PVWX,XLSX )
        ENDIF
      ENDIF
      PSBX = ABS(PSBX)
      PVOX = MAX( PGX+PATM-PVWX-PVAX,0.D+0 )
!
!---  Hydrate equilibrium temperature for no-hydrate conditions  ---
!
      TEQHX = -TABS
      YMHGOX = 0.D+0
      XHAX = 0.D+0
      XHOX = 0.D+0
      XHWX = 1.D+0
      RHOHX = 0.D+0
!
!---  Hydrate equilibrium temperature at the vapor pressure conditions  ---
!
      PVHOX = PVOX
      PVHAX = PVAX
      IF( (PVAX+PVOX).GT.EPSL ) THEN
        YMHGOX = PVHOX/(PVHAX+PVHOX)
        PEQHX = PVHOX+PVHAX
!
!---    Inhibitor correction  ---
!
        IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
          CALL HYD_PO( YMHGOX,POHX )
        ENDIF
        CALL HYD_T( TEQHX,YMHGOX,XLSX,PEQHX,POHX,RHOHX,XHAX,XHOX,XHWX )
!
!---    Hydrate equilibrium pressure at the actual temperature  ---
!
        IF( TX.LE.TEQHX ) THEN
          CALL HYD_P( TX,YMHGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
          PVHOX = YMHGOX*PEQHX
          PVHAX = PEQHX - PVHOX
!
!---    Hydrate equilibrium temperature for no-hydrate conditions  ---
!
        ELSE
          TEQHX = -TABS
          YMHGOX = 0.D+0
          XHAX = 0.D+0
          XHOX = 0.D+0
          XHWX = 1.D+0
          RHOHX = 0.D+0
        ENDIF
      ENDIF
!
!---  Liquid CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      PNX = MAX( PNX,PNC )
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_PC10 group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_PC11( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for phase condition #15
!
!     Aqueous-gas, hydrate, liquid CO2, PVHO > PVHA
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     Mobile CO2 mass - liquid-CO2 pressure
!     Mobile CH4 mass - gas pressure
!     Hydrate CO2 mass - hydrate CH4 mole fraction of formers
!     Hydrate CH4 mass - hydrate saturation
!     NaCl mass - total NaCl aqueous mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 October 2007.
!     Last Modified by MD White, PNNL, 12 October 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE HYST
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_PC11'
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure with hydrate barrier adjustment  ---
!
      HCPBFX = 1.D+0
      IF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.12 ) THEN
        HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &    SIN(SCHR(12,IZN)*SHX + SCHR(13,IZN)) )
      ENDIF
      ENPR = SCHR(2,IZN)*RHORL*GRAV*HCPBFX
!
!---  Maximum capillary pressure  ---
!
      PLX = MAX( PNX-(HDOD*RHORL*GRAV/BNL),PLX )
!
!---  Brine vapor pressure  ---
!
      IF( PSBX.GT.0.D+0 ) THEN
        CALL SP_B( TX,XLSX,PSBX )
        PPLX = MAX( PLX+PATM,PSBX )
        CALL DENS_B( TX,PPLX,XLSX,RHOBX )
        PCX = MAX( BNL*(PNX-PLX),0.D+0 )
        IF( ISLC(44).EQ.1 ) THEN
          PVWX = PSBX
        ELSE
          CALL VPL_B( TX,PSBX,PCX,RHOBX,PVWX,XLSX )
        ENDIF
      ENDIF
      PSBX = ABS(PSBX)
      PVOX = MAX( PGX+PATM-PVWX-PVAX,0.D+0 )
!
!---  Pore-space equilibrium temperature  ---
!
      TEQHX = TX
!
!---  Hydrate equilibrium pressure  ---
!
      CALL HYD_P( TEQHX,YMHGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
      PVHOX = YMHGOX*PEQHX
      PVHAX = PEQHX - PVHOX
!
!---  Liquid CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      PNX = MAX( PNX,PNC )
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Hydrate guest-molecule mass  ---
!
      TMHAX = XHAX*RHOHX*SHX
      TMHOX = XHOX*RHOHX*SHX
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_PC11 group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_PC12( PEQHX,PGX,PIX,PLX,PNX,PSBX,PVAX,PVHAX,PVOX,
     &  PVHOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,TX,TEQHX,
     &  TMHAX,TMHOX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMHGOX,IZN,N )
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
!     Flash calculation for phase condition #16
!
!     Aqueous-gas, hydrate, liquid CO2, PVHO < PVHA
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     Mobile CO2 mass - liquid-CO2 pressure
!     Mobile CH4 mass - gas pressure
!     Hydrate CO2 mass - hydrate saturation
!     Hydrate CH4 mass - hydrate CH4 mole fraction of formers
!     NaCl mass - total NaCl aqueous mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 October 2007.
!     Last Modified by MD White, PNNL, 12 October 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE HYST
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_PC12'
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
      BNL = SCHR(8,IZN)/(SIG_NL*CA_NL)
      BGN = SCHR(8,IZN)/(SIG_GN*CA_GN)
      BHL = SCHR(8,IZN)/SIG_HL
      BIL = SCHR(8,IZN)/SIG_HL
!
!---  Entry pressure with hydrate barrier adjustment  ---
!
      HCPBFX = 1.D+0
      IF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.12 ) THEN
        HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &    SIN(SCHR(12,IZN)*SHX + SCHR(13,IZN)) )
      ENDIF
      ENPR = SCHR(2,IZN)*RHORL*GRAV*HCPBFX
!
!---  Maximum capillary pressure  ---
!
      PLX = MAX( PNX-(HDOD*RHORL*GRAV/BNL),PLX )
!
!---  Brine vapor pressure  ---
!
      IF( PSBX.GT.0.D+0 ) THEN
        CALL SP_B( TX,XLSX,PSBX )
        PPLX = MAX( PLX+PATM,PSBX )
        CALL DENS_B( TX,PPLX,XLSX,RHOBX )
        PCX = MAX( BNL*(PNX-PLX),0.D+0 )
        IF( ISLC(44).EQ.1 ) THEN
          PVWX = PSBX
        ELSE
          CALL VPL_B( TX,PSBX,PCX,RHOBX,PVWX,XLSX )
        ENDIF
      ENDIF
      PSBX = ABS(PSBX)
      PVOX = MAX( PGX+PATM-PVWX-PVAX,0.D+0 )
!
!---  Pore-space equilibrium temperature  ---
!
      TEQHX = TX
!
!---  Hydrate equilibrium pressure  ---
!
      CALL HYD_P( TEQHX,YMHGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
      PVHOX = YMHGOX*PEQHX
      PVHAX = PEQHX - PVHOX
!
!---  Liquid CO2 pressure  ---
!
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      PNX = MAX( PNX,PNC )
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP_HYD_KE( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Hydrate guest-molecule mass  ---
!
      TMHAX = XHAX*RHOHX*SHX
      TMHOX = XHOX*RHOHX*SHX
!
!---  Mobile mole fraction of hydrate formers  ---
!
      IF( PVOX/EPSL.LT.EPSL ) THEN
        YMGOX = 0.D+0
      ELSE
        YMGOX = PVOX/(PVAX+PVOX)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_PC12 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_RS4( PLX,PNX,PVOX,TX,YLSX,N )
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
!     Phase resolver for phase condition #4
!
!     Aqueous, no hydrate, liquid CO2
!
!     mobile CO2 mass - liquid-CO2 pressure
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 8 February 2011.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE JACOB
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FDVT
      USE FDVS
      USE FDVP
      USE FDVGC
      USE FDVG
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 DNRX(LUK)
      REAL*8 STAX(LUK+1),STOX(LUK+1),STTX(LUK+1),STWX(LUK+1)
      REAL*8 RAP(LUK),ROP(LUK),RTP(LUK),RWP(LUK)
      REAL*8 FA(6),FO(6),FT(6),FW(6)
      REAL*8 UEL(LSV),UEN(LSV)
      REAL*8 AJ(LUK,LUK),BJ(LUK)
      INTEGER IJ(LUK)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_RS4'
!
!---  Store input primary variables  ---
!
      PLX_I = PLX
      PNX_I = PNX
      PVOX_I = PVOX
      TX_I = TX
      YLSX_I = YLSX
      ISVCX = ISVC
      IEQWX = 2
      IEQAX = 1
      IEQDOX = 3
      IF( ISLC(32).EQ.0 ) THEN
        IEQSX = 4
        IEQTX = 5
        ISVC = ISVCX-4
      ELSE
        IEQTX = 4
        ISVC = ISVCX-5
      ENDIF 
!
!---  Grid surfaces  ---
!
      NPX = NSX(N)
      NPY = NSY(N)
      NPZ = NSZ(N)
      NQX = NPX + 1
      NQY = NPY + IFLD
      NQZ = NPZ + IJFLD
!
!---  Mobile CO2 flux  ---
!
      FA(1) = -AFZ(NPZ)*(WLA(1,NPZ)+WGA(1,NPZ)+WNC(1,1,NPZ))
      FA(2) = -AFY(NPY)*(VLA(1,NPY)+VGA(1,NPY)+VNC(1,1,NPY))
      FA(3) = -AFX(NPX)*(ULA(1,NPX)+UGA(1,NPX)+UNC(1,1,NPX))
      FA(4) = AFX(NQX)*(ULA(1,NQX)+UGA(1,NQX)+UNC(1,1,NQX))
      FA(5) = AFY(NQY)*(VLA(1,NQY)+VGA(1,NQY)+VNC(1,1,NQY))
      FA(6) = AFZ(NQZ)*(WLA(1,NQZ)+WGA(1,NQZ)+WNC(1,1,NQZ))
      FAX = 0.D+0
      DO 10 MD = 1,6
        FAX = FAX + FA(MD)
   10 CONTINUE
!      PRINT *,'STAX(1) = ',STAX(1)
!      PRINT *,'SRCA(2,N) = ',SRCA(2,N)
!      PRINT *,'TMBP_A(2,N) = ',TMBP_A(2,N)
!      PRINT *,'FAX = ',FAX
!
!---  Newton Raphson iteration  ---
!
      NC = 0
  100 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        PNX = PNX_I
        GOTO 900
      ENDIF
!
!---    Energy (temperature), isothermal option  ---
!
        DNRX(IEQTX) = 1.D-7
!
!---    Assign dissolved-salt mass fraction increments,
!       isobrine option  ---
!
        IF( ISLC(32).EQ.0 ) THEN
          CALL SOL_LS( TX,XLSMX )
          XLSX = MIN( YLSX,XLSMX )
          DNRX(IEQS) = 1.D-6*XLSMX
        ENDIF
!
!---  Capillary pressure control  ---
!
      PCGNX = PG(2,N)-PNX
      PCNLX = PNX-PLX
      DNRNLX = MAX( 1.D-3,1.D-6*PCNLX )
      DNRNX = DNRNLX
      DNRLX = DNRNLX
!
!---  Water mass (aqueous pressure), decrement for
!     low capillary pressure, increment for high
!     capillary pressure  ---
!
      DNRX(IEQWX) = SIGN( DNRLX,(PCODX-PCNLX) )
!
!---  CO2 mass (liquid-CO2 pressure), decrement
!     for larger PCNLX and increment for larger PCNGX  ---
!
      DNRX(IEQAX) = SIGN( DNRNX,PCGNX-PCNLX )
!
!---  Mobile CH4 mass (CH4 vapor pressure), increment  ---
!
      DNRX(IEQDOX) = 1.D-2
!
!---  Increment the primary variables  ---
!
      DO 200 M = 2,ISVC+2
        T(M,N) = TX
        PL(M,N) = PLX
        PN(M,N) = PNX
        PVO(M,N) = PVOX
        PVA(M,N) = PVA(2,N)
        YLS(M,N) = YLSX
!        IF( M.EQ.IEQTX+2 ) THEN
!          T(M,N) = T(M,N) + DNRX(IEQTX)
!        ELSEIF( M.EQ.IEQWX+2 ) THEN
!          PL(M,N) = PL(M,N) + DNRX(IEQWX)
!        ELSEIF( M.EQ.IEQAX+2 ) THEN
!          PN(M,N) = PN(M,N) + DNRX(IEQAX)
!        ELSEIF( M.EQ.IEQDOX+2 ) THEN
!          PVO(M,N) = PVO(M,N) + DNRX(IEQDOX)
!!
!!---    Isobrine option  ---
!!
!        ELSEIF( M.EQ.IEQS+2 .AND. ISLC(32).EQ.0 ) THEN
!          YLS(M,N) = YLS(M,N) + DNRX(IEQS)
!        ENDIF
        IF( M.EQ.IEQAX+2 ) PN(M,N) = PN(M,N) + DNRX(IEQAX)
  200 CONTINUE
!
!---  Compute secondary variables  ---
!
      NLX = N
      NHX = N
      CALL PROP_HYD_KE( NLX,NHX )
!      PRINT *,'PN(2,N) = ',PN(2,N)
!      PRINT *,'SN(2,N) = ',SN(2,N)
!      PRINT *,'RHON(2,N) = ',RHON(2,N)
!      PRINT *,'SL(2,N) = ',SL(2,N)
!      PRINT *,'RHOL(2,N) = ',RHOL(2,N)
!!
!!---  Energy capacitance  ---
!!
!      UEL(1) = HL(1,N) - (PG(1,N)+PATM)/RHOL(1,N)
!      UEN(1) = HN(1,N) - (PG(1,N)+PATM)/RHON(1,N)
!      USTOX = (1.D+0-PORD(1,N))*T(1,N)
!      USLOX = PORD(1,N)*SL(1,N)*RHOL(1,N)*UEL(1)
!      USGOX = PORD(1,N)*SG(1,N)*RHOG(1,N)*UEG(1,N)
!      USNOX = PORD(1,N)*SN(1,N)*RHON(1,N)*UEN(1)
!      USHOX = PORD(1,N)*SH(1,N)*RHOH(1,N)*HH(1,N)
!      USIOX = PORD(1,N)*SI(1,N)*RHOI(1,N)*HI(1,N)
!      DO 200 M = 1,ISVC+1
!        STTX(M) = 0.D+0
!        MP = M + 1
!        UEL(MP) = HL(MP,N) - (PG(MP,N)+PATM)/RHOL(MP,N)
!        UEN(MP) = HN(MP,N) - (PG(MP,N)+PATM)/RHON(MP,N)
!        USTX = (1.D+0-PORD(MP,N))*T(MP,N)
!        USLX = PORD(MP,N)*SL(MP,N)*RHOL(MP,N)*UEL(MP)
!        USGX = PORD(MP,N)*SG(MP,N)*RHOG(MP,N)*UEG(MP,N)
!        USNX = PORD(MP,N)*SN(MP,N)*RHON(MP,N)*UEN(MP)
!        USHX = PORD(MP,N)*SH(MP,N)*RHOH(MP,N)*HH(MP,N)
!        USIX = PORD(MP,N)*SI(MP,N)*RHOI(MP,N)*HI(MP,N)
!        STTX(M) = STTX(M) + RHOS(IZN)*CPS(IZN)*(USTX-USTOX)
!        STTX(M) = STTX(M) + (USLX-USLOX)
!        STTX(M) = STTX(M) + (USGX-USGOX)
!        STTX(M) = STTX(M) + (USNX-USNOX)
!        STTX(M) = STTX(M) + (USHX-USHOX)
!        STTX(M) = STTX(M) + (USIX-USIOX)
!        STTX(M) = STTX(M)*DTI*VOL(N)
!  200 CONTINUE
!!
!!---  Energy fluxes  ---
!!
!      FT(1) = -AFZ(NPZ)*WQ(1,NPZ)
!      FT(2) = -AFY(NPY)*VQ(1,NPY)
!      FT(3) = -AFX(NPX)*UQ(1,NPX)
!      FT(4) = AFX(NQX)*UQ(1,NQX)
!      FT(5) = AFY(NQY)*VQ(1,NQY)
!      FT(6) = AFZ(NQZ)*WQ(1,NQZ)
!      FTX = 0.D+0
!      DO 220 MD = 1,6
!        FTX = FTX + FT(MD)
!  220 CONTINUE
!!
!!---  Energy equation residuals  ---
!!
!      RTS = STTX(1) - SRCT(2,N) + FTX
!!      PRINT *,'RTS = ',RTS
!      DO 230 M = 1,ISVC
!        RTP(M) = STTX(M+1) - SRCT(M+2,N) + FTX
!!        PRINT *,'RTP(',M,') = ',RTP(M)
!  230 CONTINUE
!!
!!---  Water capacitance  ---
!!
!      USLX = PORD(1,N)*XLW(1,N)*RHOL(1,N)*SL(1,N)
!      USGX = PORD(1,N)*XGW(1,N)*RHOG(1,N)*SG(1,N)
!      USNX = PORD(1,N)*XNC(3,1,N)*RHON(1,N)*SN(1,N)
!      USHX = PORD(1,N)*XHW(1,N)*RHOH(1,N)*SH(1,N)
!      USIX = PORD(1,N)*RHOI(1,N)*SI(1,N)
!      DO 300 M = 1,ISVC+1
!        STWX(M) = 0.D+0
!        MP = M + 1
!        STWX(M) = STWX(M) + 
!     &    (PORD(MP,N)*XLW(MP,N)*RHOL(MP,N)*SL(MP,N)-USLX)
!        STWX(M) = STWX(M) + 
!     &    (PORD(MP,N)*XGW(MP,N)*RHOG(MP,N)*SG(MP,N)-USGX)
!        STWX(M) = STWX(M) + 
!     &    (PORD(MP,N)*XNW(MP,N)*RHON(MP,N)*SN(MP,N)-USNX)
!        STWX(M) = STWX(M) + 
!     &    (PORD(MP,N)*XHW(MP,N)*RHOH(MP,N)*SH(MP,N)-USHX)
!        STWX(M) = STWX(M) + 
!     &    (PORD(MP,N)*RHOI(MP,N)*SI(MP,N)-USIX)
!        STWX(M) = STWX(M)*DTI*VOL(N)
!  300 CONTINUE
!!
!!---  Water fluxes  ---
!!
!      FW(1) = -AFZ(NPZ)*(WLW(1,NPZ)+WGW(1,NPZ)+WNC(3,1,NPZ))
!      FW(2) = -AFY(NPY)*(VLW(1,NPY)+VGW(1,NPY)+VNC(3,1,NPY))
!      FW(3) = -AFX(NPX)*(ULW(1,NPX)+UGW(1,NPX)+UNC(3,1,NPX))
!      FW(4) = AFX(NQX)*(ULW(1,NQX)+UGW(1,NQX)+UNC(3,1,NQX))
!      FW(5) = AFY(NQY)*(VLW(1,NQY)+VGW(1,NQY)+VNC(3,1,NQY))
!      FW(6) = AFZ(NQZ)*(WLW(1,NQZ)+WGW(1,NQZ)+WNC(3,1,NQZ))
!      FWX = 0.D+0
!      DO 320 MD = 1,6
!        FWX = FWX + FW(MD)
!  320 CONTINUE
!!
!!---  Water equation residuals  ---
!!
!      RWS = STWX(1) - SRCW(2,N) + FWX
!      PRINT *,'RWS = ',RWS
!      DO 330 M = 1,ISVC
!        RWP(M) = STWX(M+1) - SRCW(M+2,N) + FWX
!        PRINT *,'RWP(',M,') = ',RWP(M)
!  330 CONTINUE
!!
!!---  Mobile CH4 capacitance  ---
!!
!      USLX = PORD(1,N)*XLO(1,N)*RHOL(1,N)*SL(1,N)
!      USGX = PORD(1,N)*XGO(1,N)*RHOG(1,N)*SG(1,N)
!      USNX = PORD(1,N)*XNC(2,1,N)*RHON(1,N)*SN(1,N)
!      DO 510 M = 1,ISVC+1
!        STOX(M) = 0.D+0
!        MP = M + 1
!        STOX(M) = STOX(M) + 
!     &    (PORD(MP,N)*XLO(MP,N)*RHOL(MP,N)*SL(MP,N)-USLX)
!        STOX(M) = STOX(M) + 
!     &    (PORD(MP,N)*XGO(MP,N)*RHOG(MP,N)*SG(MP,N)-USGX)
!        STOX(M) = STOX(M) + 
!     &    (PORD(MP,N)*XNC(2,MP,N)*RHON(MP,N)*SN(MP,N)-USNX)
!        STOX(M) = STOX(M)*DTI*VOL(N)
!  510 CONTINUE
!!
!!---  Mobile CH4 flux  ---
!!
!      FO(1) = -AFZ(NPZ)*(WLO(1,NPZ)+WGO(1,NPZ)+WNC(2,1,NPZ))
!      FO(2) = -AFY(NPY)*(VLO(1,NPY)+VGO(1,NPY)+VNC(2,1,NPY))
!      FO(3) = -AFX(NPX)*(ULO(1,NPX)+UGO(1,NPX)+UNC(2,1,NPX))
!      FO(4) = AFX(NQX)*(ULO(1,NQX)+UGO(1,NQX)+UNC(2,1,NQX))
!      FO(5) = AFY(NQY)*(VLO(1,NQY)+VGO(1,NQY)+VNC(2,1,NQY))
!      FO(6) = AFZ(NQZ)*(WLO(1,NQZ)+WGO(1,NQZ)+WNC(2,1,NQZ))
!      FOX = 0.D+0
!      DO 520 MD = 1,6
!        FOX = FOX + FO(MD)
!  520 CONTINUE
!!
!!---  Mobile CH4 equation residuals  ---
!!
!      ROS = STOX(1) - SRCO(2,N) + TMBP_O(2,N) + FOX
!      PRINT *,'ROS = ',ROS
!      DO 530 M = 1,ISVC
!        ROP(M) = STOX(M+1) - SRCO(M+2,N) + TMBP_O(M+2,N) + FOX
!        PRINT *,'ROP(',M,') = ',ROP(M)
!  530 CONTINUE
!
!---  Mobile CO2 capacitance  ---
!
      USLX = PORD(1,N)*XLA(1,N)*RHOL(1,N)*SL(1,N)
      USGX = PORD(1,N)*XGA(1,N)*RHOG(1,N)*SG(1,N)
      USNX = PORD(1,N)*XNC(1,1,N)*RHON(1,N)*SN(1,N)
      DO 410 M = 1,ISVC+1
        STAX(M) = 0.D+0
        MP = M + 1
        STAX(M) = STAX(M) + 
     &    (PORD(MP,N)*XLA(MP,N)*RHOL(MP,N)*SL(MP,N)-USLX)
        STAX(M) = STAX(M) + 
     &    (PORD(MP,N)*XGA(MP,N)*RHOG(MP,N)*SG(MP,N)-USGX)
        STAX(M) = STAX(M) + 
     &    (PORD(MP,N)*XNC(1,MP,N)*RHON(MP,N)*SN(MP,N)-USNX)
        STAX(M) = STAX(M)*DTI*VOL(N)
  410 CONTINUE
!
!---  Mobile CO2 equation residuals  ---
!
      RAS = STAX(1) - SRCA(2,N) + TMBP_A(2,N) + FAX
!      PRINT *,'RAS = ',RAS
!      PRINT *,'RSDL(IEQA,N) = ',RSDL(IEQA,N)
      DO 430 M = 1,ISVC
        RAP(M) = STAX(M+1) - SRCA(M+2,N) + TMBP_A(M+2,N) + FAX
!        PRINT *,'RAP(',M,') = ',RAP(M)
  430 CONTINUE
!
!---  Load solution vector and Jacobian  ---
!
!      BJ(IEQTX) = -RTS
!      BJ(IEQWX) = -RWS
      BJ(IEQAX) = -RAS
!      BJ(IEQDOX) = -ROS
      DO 600 M = 1,ISVC
!        AJ(IEQTX,M) = (RTP(M)-RTS)/DNRX(M)
!        AJ(IEQWX,M) = (RWP(M)-RWS)/DNRX(M)
        AJ(IEQAX,M) = (RAP(M)-RAS)/DNRX(M)
!        AJ(IEQDOX,M) = (ROP(M)-ROS)/DNRX(M)
  600 CONTINUE
!      DO 620 I = 1,ISVC
!        DO 610 J = 1,ISVC
!          PRINT *,'AJ(I,J) = ',AJ(I,J)
!  610   CONTINUE
!        PRINT *,'BJ(I) = ',BJ(I)
!  620 CONTINUE
!
!---  Solve linear system  ---
!
      JP = ISVC
      KP = 6
      CALL LUDCMP( AJ,JP,KP,IJ,DJ )
      CALL LUBKSB( AJ,JP,KP,IJ,BJ )
!
!---  Primary updates  ---
!
!      DPT = BJ(IEQTX)
!      DPW = BJ(IEQWX)
      DPX = 1.D+3
      DPA = BJ(IEQAX)
      DPA = SIGN( MIN( ABS(DPX),ABS(DPA) ),DPA )
!      DPO = BJ(IEQDOX)
!      PRINT *,'DPT = ',DPT
!      PRINT *,'DPW = ',DPW
!      PRINT *,'DPA = ',DPA
!      PRINT *,'DPO = ',DPO
!      PRINT *,'NC = ',NC
!      TX = TX + DPT
!      PLX = PLX + DPW
      IF( PNX+DPA.LT.PLX ) DPA = 6.D-1*DPA
      PNX = MAX( PNX+DPA,PLX )
!      IF( PVOX.LT.EPSL .AND. DPO.LT.EPSL ) DPO = 0.D+0
!      PVOX = PVOX + DPO
      IF( (ABS(DPA)/(PNX+PATM)).GT.1.D-6 ) GOTO 100
!
!---  Reset increment looping index  ---
!
  900 CONTINUE
      ISVC = ISVCX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_RS4 group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_RS7( PLX,PGX,PVAX,PVOX,TX,YLSX,N )
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
!     Phase resolver for phase condition #7
!
!     Aqueous-gas, no hydrate, no liquid-CO2
!
!     Mobile CO2 mass - CO2 vapor pressure
!     Mobile CH4 mass - gas pressure
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 8 February 2011.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE PORMED
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
      USE FLUXGC
      USE FLUXD
      USE FDVT
      USE FDVS
      USE FDVP
      USE FDVGC
      USE FDVG
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 DNRX(LUK)
      REAL*8 STAX(LUK+1),STOX(LUK+1),STTX(LUK+1),STWX(LUK+1)
      REAL*8 RAP(LUK),ROP(LUK),RTP(LUK),RWP(LUK)
      REAL*8 FA(6),FO(6),FT(6),FW(6)
      REAL*8 UEL(LSV),UEN(LSV)
      REAL*8 AJ(LUK,LUK),BJ(LUK)
      INTEGER IJ(LUK)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_RS7'
!
!---  Store input primary variables  ---
!
      IZN = IZ(N)
      PLX_I = PLX
      PGX_I = PGX
      PVAX_I = PVAX
      PVOX_I = PVOX
      RVAOX = 0.D+0
      IF( (PVAX+PVOX)/EPSL.GT.EPSL ) RVAOX = PVAX/(PVAX+PVOX)
      TX_I = TX
      YLSX_I = YLSX
!
!---  Brine salt mass fraction  ---
!      
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      IF( ISLC(7).EQ.2 .OR. ISLC(7).EQ.3 )
     &  CALL SFT_L( TX,XLSX,SIG_GL )
      CA_NL = ((SIG_NL**2) - (SIG_GN**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_NL)
      CA_GN = ((SIG_GN**2) - (SIG_NL**2) + (SIG_GL**2))
     &  /(2.D+0*SIG_GL*SIG_GN)
!
!---  Interphase scaling factors  ---
!
      BGL = SCHR(8,IZN)/SIG_GL
!
!---  Brine vapor pressure  ---
!
      CALL SP_B( TX,XLSX,PSBX )
      PPLX = MAX( PLX+PATM,PSBX )
      CALL DENS_B( TX,PPLX,XLSX,RHOBX )
      ISVCX = ISVC
      ISVC = 1
      IEQDOX = 1
      IEQAX = 2
!
!---  Grid surfaces  ---
!
      NPX = NSX(N)
      NPY = NSY(N)
      NPZ = NSZ(N)
      NQX = NPX + 1
      NQY = NPY + IFLD
      NQZ = NPZ + IJFLD
!
!---  Mobile CO2 flux  ---
!
      FA(1) = -AFZ(NPZ)*(WLA(1,NPZ)+WGA(1,NPZ)+WNC(1,1,NPZ))
      FA(2) = -AFY(NPY)*(VLA(1,NPY)+VGA(1,NPY)+VNC(1,1,NPY))
      FA(3) = -AFX(NPX)*(ULA(1,NPX)+UGA(1,NPX)+UNC(1,1,NPX))
      FA(4) = AFX(NQX)*(ULA(1,NQX)+UGA(1,NQX)+UNC(1,1,NQX))
      FA(5) = AFY(NQY)*(VLA(1,NQY)+VGA(1,NQY)+VNC(1,1,NQY))
      FA(6) = AFZ(NQZ)*(WLA(1,NQZ)+WGA(1,NQZ)+WNC(1,1,NQZ))
      FAX = 0.D+0
      DO 10 MD = 1,6
        FAX = FAX + FA(MD)
   10 CONTINUE
!
!---  Mobile CH4 flux  ---
!
      FO(1) = -AFZ(NPZ)*(WLO(1,NPZ)+WGO(1,NPZ)+WNC(2,1,NPZ))
      FO(2) = -AFY(NPY)*(VLO(1,NPY)+VGO(1,NPY)+VNC(2,1,NPY))
      FO(3) = -AFX(NPX)*(ULO(1,NPX)+UGO(1,NPX)+UNC(2,1,NPX))
      FO(4) = AFX(NQX)*(ULO(1,NQX)+UGO(1,NQX)+UNC(2,1,NQX))
      FO(5) = AFY(NQY)*(VLO(1,NQY)+VGO(1,NQY)+VNC(2,1,NQY))
      FO(6) = AFZ(NQZ)*(WLO(1,NQZ)+WGO(1,NQZ)+WNC(2,1,NQZ))
      FOX = 0.D+0
      DO 20 MD = 1,6
        FOX = FOX + FO(MD)
   20 CONTINUE
!
!---  Newton Raphson iteration  ---
!
      NC = 0
  100 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        PGX = PGX_I
        GOTO 900
      ENDIF
!
!---  Capillary pressure control  ---
!
      PCODX = 5.D-1*HDOD*RHORL*GRAV
      PCGLX = PGX-PLX
      DNRGLX = MAX( 1.D-3,1.D-6*PCGLX )
!
!---  Mobile CH4 mass (gas pressure), increment for
!     low capillary pressure, decrement for high
!     capillary pressure  ---
!
      DNRX(IEQDOX) = SIGN( DNRGLX,(PCODX-PCGLX) )
!
!---  Increment the primary variables  ---
!
      DO 200 M = 2,ISVC+2
        T(M,N) = TX
        PL(M,N) = PLX
        PN(M,N) = PNX
        PVA(M,N) = PVAX
        PVO(M,N) = PVOX
        PG(M,N) = PGX
        YLS(M,N) = YLSX
        IF( M.EQ.IEQDOX+2 ) THEN
          PG(M,N) = PG(M,N) + DNRX(IEQDOX)
        ENDIF
        PCX = MAX( BGL*(PG(M,N)-PLX),0.D+0 )
        IF( ISLC(44).EQ.1 ) THEN
          PVWX = PSBX
        ELSE
          CALL VPL_B( TX,PSBX,PCX,RHOBX,PVWX,XLSX )
        ENDIF
        PVX = MAX( PG(M,N)+PATM-PVWX,0.D+0 )
        IF( PVX/EPSL.GT.EPSL ) THEN
          PVA(M,N) = MIN( PVX*RVAOX,PVX )
          PVO(M,N) = MAX( (PVX-PVA(M,N)),0.D+0 )
        ELSE
          PVA(M,N) = 0.D+0
          PVO(M,N) = 0.D+0
        ENDIF
  200 CONTINUE
!
!---  Compute secondary variables  ---
!
      NLX = N
      NHX = N
      CALL PROP_HYD_KE( NLX,NHX )
!
!---  Mobile CO2 capacitance  ---
!
      USLX = PORD(1,N)*XLA(1,N)*RHOL(1,N)*SL(1,N)
      USGX = PORD(1,N)*XGA(1,N)*RHOG(1,N)*SG(1,N)
      USNX = PORD(1,N)*XNC(1,1,N)*RHON(1,N)*SN(1,N)
      DO 410 M = 1,ISVC+1
        STAX(M) = 0.D+0
        MP = M + 1
        STAX(M) = STAX(M) + 
     &    (PORD(MP,N)*XLA(MP,N)*RHOL(MP,N)*SL(MP,N)-USLX)
        STAX(M) = STAX(M) + 
     &    (PORD(MP,N)*XGA(MP,N)*RHOG(MP,N)*SG(MP,N)-USGX)
        STAX(M) = STAX(M) + 
     &    (PORD(MP,N)*XNC(1,MP,N)*RHON(MP,N)*SN(MP,N)-USNX)
        STAX(M) = STAX(M)*DTI*VOL(N)
  410 CONTINUE
!
!---  Mobile CO2 equation residuals  ---
!
      RAS = STAX(1) - SRCA(2,N) + TMBP_A(2,N) + FAX
      DO 430 M = 1,ISVC
        RAP(M) = STAX(M+1) - SRCA(M+2,N) + TMBP_A(M+2,N) + FAX
  430 CONTINUE
!
!---  Mobile CH4 capacitance  ---
!
      USLX = PORD(1,N)*XLO(1,N)*RHOL(1,N)*SL(1,N)
      USGX = PORD(1,N)*XGO(1,N)*RHOG(1,N)*SG(1,N)
      USNX = PORD(1,N)*XNC(2,1,N)*RHON(1,N)*SN(1,N)
      DO 510 M = 1,ISVC+1
        STOX(M) = 0.D+0
        MP = M + 1
        STOX(M) = STOX(M) + 
     &    (PORD(MP,N)*XLO(MP,N)*RHOL(MP,N)*SL(MP,N)-USLX)
        STOX(M) = STOX(M) + 
     &    (PORD(MP,N)*XGO(MP,N)*RHOG(MP,N)*SG(MP,N)-USGX)
        STOX(M) = STOX(M) + 
     &    (PORD(MP,N)*XNC(2,MP,N)*RHON(MP,N)*SN(MP,N)-USNX)
        STOX(M) = STOX(M)*DTI*VOL(N)
  510 CONTINUE
!
!---  Mobile CH4 equation residuals  ---
!
      ROS = STOX(1) - SRCO(2,N) + TMBP_O(2,N) + FOX
      DO 530 M = 1,ISVC
        ROP(M) = STOX(M+1) - SRCO(M+2,N) + TMBP_O(M+2,N) + FOX
  530 CONTINUE
!
!---  Load solution vector and Jacobian  ---
!
!      BJ(IEQAX) = -RAS
      BJ(IEQDOX) = -ROS -RAS
      DO 600 M = 1,ISVC
!        AJ(IEQAX,M) = (RAP(M)-RAS)/DNRX(M)
        AJ(IEQDOX,M) = (ROP(M)-ROS)/DNRX(M) + (RAP(M)-RAS)/DNRX(M)
  600 CONTINUE
!
!---  Solve linear system  ---
!
      JP = ISVC
      KP = 6
      CALL LUDCMP( AJ,JP,KP,IJ,DJ )
      CALL LUBKSB( AJ,JP,KP,IJ,BJ )
!
!---  Primary updates  ---
!
      DPX = 1.D+3
      DPO = BJ(IEQDOX)
      DPO = SIGN( MIN( ABS(DPX),ABS(DPO) ),DPO )
      IF( PGX+DPO.LT.PLX ) DPO = 6.D-1*DPO
      PGX = MAX( PGX+DPO,PLX )
      IF( (ABS(DPO)/(PGX+PATM)).GT.1.D-6 ) GOTO 100
!
!---  Converged gas pressure  ---
!
      PCX = MAX( BGL*(PGX-PLX),0.D+0 )
      IF( ISLC(44).EQ.1 ) THEN
        PVWX = PSBX
      ELSE
        CALL VPL_B( TX,PSBX,PCX,RHOBX,PVWX,XLSX )
      ENDIF
      PVX = MAX( PGX+PATM-PVWX,0.D+0 )
      IF( PVX/EPSL.GT.EPSL ) THEN
        PVAX = MIN( PVX*RVAOX,PVX )
        PVOX = MAX( (PVX-PVAX),0.D+0 )
      ELSE
        PVAX = 0.D+0
        PVOX = 0.D+0
      ENDIF
!
!---  Reset increment looping index  ---
!
  900 CONTINUE
      ISVC = ISVCX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_RS7 group  ---
!
      RETURN
      END
      


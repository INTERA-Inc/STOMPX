!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ATMOS_C( PA,RH,RN,TA,UZ,IERR )
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
!     Atmosphseric conditions
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 20 August 2002.
!     Last Modified by MD White, PNNL, 20 August 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PLT_ATM
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
      SUB_LOG(ISUB_LOG) = '/ATMOS_C'
!
!---  Set local time variable  ---
!
      IERR = 0
      TMZ = TM
      TMZO = TM-DT
      IF( NSTEP-NRST.EQ.0 ) THEN
        TMZ = TMZ*(1.D+0+EPSL)+EPSL
        TMZO = TMZ
      ENDIF
!
!---  Cyclic atmospheric conditions  ---
!
      IF( IATM_C.EQ.1 ) THEN
        TMZ = MOD( TM,ATMOS(NATM_T,1) )
        TMZO =  MOD( TM-DT,ATMOS(NATM_T,1) )
      ENDIF
!
!---  Simulation time less than starting atmospheric
!     conditions time  ---
!
      IF( TMZ.LE.ATMOS(1,1) ) THEN
        IERR = 1
        GOTO 150
      ENDIF
!
!---  Interpolate atmospheric conditions  ---
!
      IF( NATM_T.EQ.1 ) THEN
        TA = ATMOS(1,2)
        PA = ATMOS(1,3)
        RH = ATMOS(1,4)
        RN = ATMOS(1,5)
        UZ = ATMOS(1,6)
      ELSE
        DO 100 NT = 2,NATM_T
          IF( TMZ.LE.ATMOS(NT,1) ) THEN
           TDX = ATMOS(NT,1)-ATMOS(NT-1,1)
           TFX = (TMZ-ATMOS(NT-1,1))/TDX
           TA = ATMOS(NT-1,2) + TFX*(ATMOS(NT,2)-ATMOS(NT-1,2))
           PA = ATMOS(NT-1,3) + TFX*(ATMOS(NT,3)-ATMOS(NT-1,3))
           RH = ATMOS(NT-1,4) + TFX*(ATMOS(NT,4)-ATMOS(NT-1,4))
           RN = ATMOS(NT-1,5) + TFX*(ATMOS(NT,5)-ATMOS(NT-1,5))
           UZ = ATMOS(NT-1,6) + TFX*(ATMOS(NT,6)-ATMOS(NT-1,6))
           GOTO 110
          ENDIF
  100   CONTINUE
  110   CONTINUE
!
!---    Integrate solar radiation  ---
!
        RNSZ = 0.D+0
        TMSZ = 0.D+0
        IERR = 1
        DO 120 NT = 2,NATM_T
          IF( TMZO.LE.ATMOS(NT,1) ) THEN
           TDX = ATMOS(NT,1)-ATMOS(NT-1,1)
           TFX = (TMZO-ATMOS(NT-1,1))/TDX
           RNSX = ATMOS(NT-1,5) + TFX*(ATMOS(NT,5)-ATMOS(NT-1,5))
           TMSX = TMZO
           NTX = NT
           IERR = 0
           GOTO 130
          ENDIF
 120    CONTINUE
 130    CONTINUE
!
!---    Simulation time greater than ending atmospheric
!       conditions time  ---
!
        IF( IERR.EQ.1 ) GOTO 150
        DO 140 NT = NTX,NATM_T
          IF( TMZ.LE.ATMOS(NT,1) ) THEN
           TDX = ATMOS(NT,1)-ATMOS(NT-1,1)
           TFX = (TMZ-ATMOS(NT-1,1))/TDX
           RNEX = ATMOS(NT-1,5) + TFX*(ATMOS(NT,5)-ATMOS(NT-1,5))
           TMEX = TMZ
           RNSZ = RNSZ + 5.D-1*(RNSX+RNEX)*(TMEX-TMSX)
           TMSZ = TMSZ + (TMEX-TMSX)
           IF( TMSZ.GT.EPSL ) THEN
             RN = RNSZ/TMSZ
           ELSE
             RN = 0.D+0
           ENDIF
           GOTO 150
          ELSE
           RNSZ = RNSZ + 5.D-1*(RNSX+ATMOS(NT,5))*(ATMOS(NT,1)-TMSX)
           TMSZ = TMSZ + (ATMOS(NT,1)-TMSX)
           TMSX = ATMOS(NT,1)
           RNSX = ATMOS(NT,5)
          ENDIF
  140   CONTINUE
!
!---    Simulation time greater than ending atmospheric
!       conditions time  ---
!
        IERR = 1
      ENDIF
  150 CONTINUE
!
!---  Minimum wind speed, 1 cm/s  ---
!
      UZ = MAX( UZ,1.D-2 )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ATMOS_C group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BARE_S( BCX,DB,GB,M,N,NB )
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
!     Bare soil boundary condition.
!
!     Atmospheric variables
!
!     TA  Atmospheric temperature, C
!     RH  Atmospheric relative humidity, %
!     PA  Atmospheric pressure, Pa
!     RN  Solar radiation, W/m^2
!     UZ  Wind speed, m/s
!
!     Atmospheric properties
!
!     E_A - water-vapor density, kg/m^3
!     PSW_A - saturated water-vapor pressure, Pa
!     PGW_A - water-vapor partial pressure, Pa
!     HGW_A - water-vapor enthalpy, J/kg
!     UGW_A - water-vapor internal energy, J/kg
!
!     Radiation variables
!
!     RN_SNS - net short-wave radiation into the ground surface, W/m^2
!     RN_LD - downward long-wave radiation, W/m^2
!     RN_LDS - downward long-wave radiation to the ground surface, W/m^2
!     RN_LUS - upward long-wave radiation from the ground surface, W/m^2
!     RN_LNS - net long-wave radiation into the ground surface, W/m^2
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 August 2002.
!     Last Modified by MD White, PNNL, 12 August 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOURC
      USE SOLTN
      USE PORMED
      USE PLT_ATM
      USE HYST
      USE GRID
      USE FILES
      USE FDVT
      USE FDVP
      USE FDVG
      USE CONST
      USE BCVT
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
      REAL*8 RW_S(5),RE_S(5)
      REAL*8 L_SAMX,L_SAX
      REAL*8 DNRX(4)
      REAL*8 BCX(LBCV)
      REAL*8 L_NS,L_SA,K_VK,KVIS_A
      REAL*8 RKL_S(3)
      REAL*8 FX(2)
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: BJ
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: AJ
      INTEGER, DIMENSION(:), ALLOCATABLE ::  IJ
!
!----------------------Data Statements---------------------------------!
!
      SAVE K_VK,ED_DC
      SAVE A_SA,CP_A,CPG_A,D,D_O,D_X,DZ_ANG,E_A,EM_CS,F_CC,GBX,
     &  HG_A,HG_AI,HGA_A,HGA_AI,HGW_A,HGW_AI,HL_A,HLW_A,IBD,IERR,
     &  INDX,J_DAY,K,KVIS_A,L_SA,MYR,NSTEP_O,PA,PGA_A,PGW_A,
     &  PR_A,PSW_A,RA_BSMX,RA_BSX,RH,RHO_A,RHOG_A,RHOG_AI,
     &  RHOGA_A,RHOGA_AI,RHOGW_A,RHOGW_AI
      SAVE RHOL_A,RN,RN_LD,RN_LDS,SIGMA,T_DP,TA,TA_K,TA_Y,TAI,THD_A,
     &  THK_A,TMZ,TMZ_O,UGA_A,UGA_AI,UGW_A,UGW_AI,UZ,UZ_FX,
     &  VIS_A,VLG_1,VLG_2,XGA_A,XGA_AI,XGW_A,XGW_AI,XLA_A,XLW_A,
     &  XMLA_A,XMLW_A,MPX
      SAVE AJ,BJ,IJ
      DATA K_VK,ED_DC / 4.1D-1,2.5D+0 /
      DATA TMZ_O,D_O / -1.D+20,0.D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BARE_S'
!
!---  Dynamic memory allocation  ---
!
      KP = 2+LPLANT*3
      IF( .NOT.ALLOCATED(AJ) ) THEN
        ALLOCATE( AJ(1:KP,1:KP),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: AJ'
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
      IF( .NOT.ALLOCATED(BJ) ) THEN
        ALLOCATE( BJ(1:KP),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: BJ'
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
      IF( .NOT.ALLOCATED(IJ) ) THEN
        ALLOCATE( IJ(1:KP),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: IJ'
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
!
!---  Skip atmospheric calculations  ---
!
      MPOX = 1
      NC = 0
      IF( M.EQ.2 ) MPX = 3
      IF( NSTEP-NRST.EQ.0 ) NSTEP_O = -1
      IF( M.GT.2 .OR. NSTEP.EQ.NSTEP_O ) GOTO 300
      NSTEP_O = NSTEP
!
!---  Initialize variables  ---
!
      GBX = GB/DB
      IBD = ABS(IBCD(NB))
      K = KD(N)
      IF( IBCD(NB).EQ.-3 ) THEN
        A_SA = AFZ(NSZ(N))
      ELSEIF( IBCD(NB).EQ.-2 ) THEN
        A_SA = AFY(NSY(N))
      ELSEIF( IBCD(NB).EQ.-1 ) THEN
        A_SA = AFX(NSX(N))
      ELSEIF( IBCD(NB).EQ.1 ) THEN
        A_SA = AFX(NSX(N)+1)
      ELSEIF( IBCD(NB).EQ.2 ) THEN
        A_SA = AFY(NSY(N)+IFLD)
      ELSEIF( IBCD(NB).EQ.3 ) THEN
        A_SA = AFZ(NSZ(N)+IJFLD)
      ENDIF
!
!---  Time of year  ---
!
      TMZ = TM
      IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
      IF( TMZ.EQ.TMZ_O ) THEN
        D = D_O
      ELSE
        D = TMZ + ATMST
        MYR = 1900
  100   CONTINUE
        MYR = MYR+1
        D_X = D
        D = D - 3.1536D+7
        IF( MOD(MYR,1000).EQ.0 ) THEN
          IF( MOD(MYR,400).EQ.0 ) THEN
            D = D - 8.64D+4
          ENDIF
        ELSEIF( MOD(MYR,4).EQ.0 ) THEN
          D = D - 8.64D+4
        ENDIF
        IF( D.GT.0.D+0 ) GOTO 100
        D = D_X
        D_O = D
        TMZ_O = TMZ
        J_DAY = INT(D/8.64D+4) + 1
      ENDIF
!
!---  Set atmospheric conditions ---
!
      CALL ATMOS_C( PA,RH,RN,TA,UZ,IERR )
!
!---  Simulation time greater than ending atmospheric conditions time
!     skip Shuttleworth-Wallace calculations  ---
!
      IF( IERR.EQ.1 )  GOTO 1000
!
!---  Downward long-wave radiation,
!     W/m^2 ground surface  ---
!
!     RN_LD - downward long-wave radiation, W/m^2
!     EM_CS - clear-sky emissivity
!     F_CC - fractional cloud cover
!     SIGMA - Stefan-Boltzmann constant, W/m^2 K^4
!     TA_K - atmospheric temperature, K
!
      CALL DEW_PT( RH,TA,T_DP )
      CALL CLOUD( RN,F_CC,J_DAY,DZ_ANG )
      EM_CS = MIN( 0.741D+0 + 0.0062*T_DP,1.D+0 )
      SIGMA = 5.67D-8
      TA_K = TA + TABS
      RN_LD = EM_CS*(1.D+0-0.84D+0*F_CC)*SIGMA*(TA_K**4)
     &  + 0.84D+0*F_CC*SIGMA*(TA_K**4)
!
!---  Downward long-wave radiation to the ground surface,
!     W/m^2 ground surface  ---
!
!     RN_LD - downward long-wave radiation, W/m^2
!     RN_LDS - downward long-wave radiation to ground surface, W/m^2
!
      RN_LDS = RN_LD
!
!---  Atmospheric properties  ---
!
!     E_A - water-vapor density, kg/m^3
!     PSW_A - saturated water-vapor pressure, Pa
!     PGW_A - atmospheric water-vapor pressure, Pa
!     PGA_A - atmospheric air partial pressure, Pa
!     CPG_A - atmospheric specific heat, J/kg
!     RHOGW_A - atmospheric water-vapor density, kg/m^3
!     RHOGA_A - atmospheric air density, kg/m^3
!     RHOG_A - atmospheric gas density, kg/m^3
!     HGW_A - atmospheric water-vapor enthalpy, J/kg
!     UGW_A - atmospheric water-vapor internal energy, J/kg
!     HGA_A - atmospheric air enthalpy, J/kg
!     UGA_A - atmospheric air internal energy, J/kg
!     RHO_A - dry-air density, kg/m^3
!     VIS_A - dry-air viscosity, Pa s
!     KVIS_A - dry-air kinematic viscosity, m^2/s
!     CP_A - dry-air specific heat, J/kg
!     THK_A - dry-air thermal conductivity, W/m K
!     THD_A - dry-air thermal diffusivity, m^2/s
!     PR_A - dry-air Prandtl number
!
      INDX = 0
      CALL REGION_4( TA,PSW_A,INDX )
      PGW_A = RH*PSW_A
      INDX = 1
      CALL WATGSD( TA,PGW_A,E_A,INDX )
      CALL WATGSH( TA,PSW_A,E_A,HGW_A,UGW_A )
      PGA_A = MAX( PA-PGW_A,0.D+0 )
      INDX = 1
      CALL WATGSD( TA,PGW_A,RHOGW_A,INDX )
      CALL AIRGSD( TA,PGA_A,RHOGA_A )
      RHOG_A = RHOGW_A+RHOGA_A
      XGW_A = RHOGW_A/RHOG_A
      XGA_A = RHOGA_A/RHOG_A
      CALL WATGSH( TA,PSW_A,RHOGW_A,HGW_A,UGW_A )
      CALL AIRGSH( TA,PGA_A,HGA_A,UGA_A )
      HG_A = XGW_A*HGW_A + XGA_A*HGA_A
      TAI = TA + 1.D-3
      INDX = 1
      CALL WATGSD( TAI,PGW_A,RHOGW_AI,INDX )
      CALL AIRGSD( TAI,PGA_A,RHOGA_AI )
      RHOG_AI = RHOGW_AI+RHOGA_AI
      XGW_AI = RHOGW_AI/RHOG_AI
      XGA_AI = RHOGA_AI/RHOG_AI
      CALL WATGSH( TAI,PSW_A,RHOGW_AI,HGW_AI,UGW_AI )
      CALL AIRGSH( TAI,PGA_A,HGA_AI,UGA_AI )
      HG_AI = XGW_AI*HGW_AI + XGA_AI*HGA_AI
      CPG_A = (HG_AI-HG_A)/1.D-3
      CALL AIRGSD( TA,PA,RHO_A )
      CALL AIRGSV( TA,VIS_A )
      KVIS_A = VIS_A/RHO_A
      CALL AIRGSC( TA,CP_A )
      CALL AIRGSK( TA,THK_A )
      THD_A = THK_A/(RHO_A*CP_A)
      PR_A = KVIS_A/THD_A
!
!---  Precipitation properties  ---
!
      XMLA_A = PGA_A/HCAW
      XMLW_A = MAX( 1.D+0-XMLA_A,0.D+0 )
      XLW_A = XMLW_A*WTMW/(XMLA_A*WTMA + XMLW_A*WTMW)
      XLA_A = MAX( 1.D+0-XLW_A,0.D+0 )
      TA_Y = MAX( TA,1.D-2 )
      CALL WATLQD( TA_Y,PA,RHOL_A )
      CALL WATLQH( TA_Y,PGA_A,HLW_A )
      HL_A = XLW_A*HLW_A + XLA_A*HGA_A
!
!---  Eddy diffusion resistance for bare surface  ---
!
!     ATMC(1) - wind speed reference height, m
!     K_VK - von Karman constant, 0.4
!     UZ  Wind speed, m/s
!     UZ_FX  Friction velocity, m/s
!     RA_BSX - aerodynamic resistance between ground and
!       reference height for bare surface, s/m
!     RA_BSMX - maximum aerodynamic resistance between ground and
!       reference height for bare surface, s/m
!     ATMC(6) - momentum transport surface roughness length, m
!     ATMC(7) - heat and mass transport surface roughness length, m
!
!      RA_BS = ((LOG(ATMC(1)/ATMC(6)))**2)/((K_VK**2)*UZ)
      RA_BSMX = RHOG_A*CPG_A*ATMC(2)/THK_A
      VLG_1 = LOG((ATMC(1)+ATMC(6))/ATMC(6))
      VLG_2 = LOG((ATMC(2)+ATMC(7))/ATMC(7))
      UZ_FX = K_VK*UZ/VLG_1
      RA_BSX = MIN( VLG_2/(K_VK*UZ_FX),RA_BSMX )
!
!---  Precipitation  ---
!
      IF( BCX(2)/EPSL.GT.EPSL ) THEN
        L_SA = -BCX(2)
      ELSE
        L_SA = 0.D+0
      ENDIF
!
!---  Initialize convergence variables  ---
!
  300 CONTINUE
      DP_SO = 0.D+0
      DP_SOO = 0.D+0
      RHOL_SO = RHOLB(1,NB)
      RHOG_SO = RHOGB(1,NB)
      SL_SO = SLB(1,NB)
      HL_SO = HLB(1,NB)
      UG_SO = UEGB(1,NB)
      T_SOX = TB(1,NB)
      XGW_SO = XGWB(1,NB)
      XLW_SO = XLWB(1,NB)
!
!---  Skip iterative calculations  ---
!
      IF( (M.NE.2 .OR. NITER.GT.1) .AND. ISLC(24).EQ.10 ) MPX = 1
      IF( NITER.GT.(NRIMX/2) ) MPX = 1
!
!---  Guess unknowns  ---
!
      PG_S = PA
      IF( M.EQ.2 .AND. (NSTEP-NRST).EQ.0 ) THEN
        T_S = TA
        RHX = SQRT(SQRT(RH))
        CALL WATLQD( T_S,PG_S,RHOL_S )
        IF( RH.LT.9.D-1 ) THEN
          PL_S = LOG(RHX)*(RHOL_S*RCW*(T_S+TABS)) + PG_S
        ELSEIF( BCX(2)/EPSL.GT.EPSL ) THEN
          PL_S = 0.99D+0*PG_S
        ELSE
          PL_S = RH*PG_S
        ENDIF
        PL_S = MIN( PL_S,PG_S+SWBCD*RHOL_A*GRAV )
      ELSE
        T_S = T_SO(2,NB)
!
!---    Guess saturation under rainfall conditions  ---
!
        IF( M.EQ.2 .AND. NITER.EQ.1 .AND. L_SA.LT.-EPSL ) THEN
          SLX = MIN( 1.D+0,ABS(L_SA)*VISRL/(PERM(3,IZ(N))*RHORL*GRAV) )
          CALL CAP_GT( IZ(N),SLX,CPGLX )
          PL_S = MAX( PG_S-CPGLX,PL_SO(2,NB) )
        ELSE
          PL_S = PL_SO(2,NB)
        ENDIF
        PL_S = MIN( PL_S,PG_S+SWBCD*RHOL_A*GRAV )
      ENDIF
      T_SMN = 0.D+0
      PL_SMN = 0.D+0
      SL_SMN = 0.D+0
  310 CONTINUE
      NC = NC+1
!
!---  Primary variable increments  ---
!
      IF( MPX.GT.0 ) THEN
        DNRX(1) = 1.D-4
        DNRX(2) = -MAX( 1.D-1,1.D-6*MAX(PG_S-PL_S,0.D+0) )
      ELSE
        DNRX(1) = -MAX( 1.D-1,1.D-6*MAX(PG_S-PL_S,0.D+0) )
        DNRX(2) = 1.D-4
      ENDIF
!
!---  Loop over number of primary unknowns  ---
!
!     T_S - surface temperature, C
!     PL_S - surface aqueous pressure, Pa (abs)
!     MPX (+2) fix surface aqueous pressure 
!     MPX (-2) fix surface temperature
!     MPX (+3) solve surface aqueous pressure and temperature
!
      DO 400 MP = 1,ABS(MPX)
!
!---    Increment primary variables  ---
!
        T_SX = T_S
        PL_SX = PL_S
        IF( MPX.GT.0 ) THEN
          IF( MP.EQ.2 ) T_SX = T_S + DNRX(1)
          IF( MP.EQ.3 ) PL_SX = PL_S + DNRX(2)
        ELSE
          IF( MP.EQ.2 ) PL_SX = PL_S + DNRX(1)
          IF( MP.EQ.3 ) T_SX = T_S + DNRX(2)
        ENDIF
        D_NS = ABS(DB)
!
!---    Limit property data temperature to above freezing  ---
!
        T_SY = MAX( T_SX,1.D-1 )
!
!---    Water vapor pressure at ground surface
!       w/ vapor pressure lowering  ---
!
        CALL WATLQD( T_SY,PG_S,RHOL_S )
        INDX = 0
        CALL REGION_4( T_SX,PSW_S,INDX )
        PC_S = MAX( PG_S-PL_SX,0.D+0 )
        IF( ISLC(44).EQ.1 ) THEN
          PGW_S = PSW_S
        ELSE
          CALL VPLWR( PC_S,PSW_S,T_SX,RHOL_S,RCW,PGW_S )
        ENDIF
!
!---    Water-vapor density at ground surface, kg/m^3  ---
!
        PGA_S = MAX( PG_S-PGW_S,0.D+0 )
        XMLA_S = PGA_S/HCAW
        XMLW_S = MAX( 1.D+0-XMLA_S,0.D+0 )
        XLW_S = XMLW_S*WTMW/(XMLA_S*WTMA + XMLW_S*WTMW)
        XLA_S = MAX( 1.D+0-XLW_S,0.D+0 )
        INDX = 1
        CALL WATGSD( T_SX,PGW_S,RHOGW_S,INDX )
        CALL AIRGSD( T_SX,PGA_S,RHOGA_S )
        RHOG_S = RHOGW_S+RHOGA_S
        XGW_S = RHOGW_S/RHOG_S
        XGA_S = RHOGA_S/RHOG_S
        E_S = RHOGW_S
!
!---    Gas viscosity at ground surface  ---
!
        CALL WATGSV( T_SX,VISW_S )
        CALL AIRGSV( T_SX,VISA_S )
        XMGW_S = PGW_S/PG_S
        XMGA_S = PGA_S/PG_S
        CALL GASVIS( XMGW_S,ZERO,XMGA_S,VISW_S,SMALL,VISA_S,VISG_S )
!
!---    Aqueous pressure at ground surface, from vapor-pressure
!       lowering function  ---
!
        CALL WATLQD( T_SY,PG_S,RHOL_S )
        INDX = 0
        CALL REGION_4( T_SX,PSW_S,INDX )
!
!---    Saturated water-vapor density at ground surface
!       temperature, kg/m^3  ---
!
        CALL WATGSD( T_SX,PSW_S,PE_S,INDX )
!
!---    Aqueous viscosity at the ground surface, Pa s  ---
!
        CALL WATLQV( T_SY,PG_S,PSW_S,VISL_S )
!
!---    Saturation and relative permeability at ground surface  ---
!
        CALL KSP_GT( N,IZ(N),M,PG_S,PL_SX,SL_S,SG_S,RKL_S,RKG_S )
        IF( MP.EQ.1 ) SL_SX = SL_S
!
!---    Net short-wave radiation intercepted by ground surface,
!       W/m^2 ground surface  ---
!
!       RN_SNS - net short-wave radiation transitted to the
!         ground surface, W/m^2 ground
!       RN - short-wave radiation, W/m^2
!       ALB_S - short-wave albedo of ground surface (0.1 - 0.4)
!
!---    Moisture ground-surface solar angle model  ---
!
        IF( IALB(IZ(N)).GT.0 ) THEN
          ALB_S = ALBEDO(2,IZ(N)) + 
     &      (ALBEDO(1,IZ(N))-ALBEDO(2,IZ(N)))*
     &      EXP(-ALBEDO(3,IZ(N))*SL_S*PORD(2,N))
        ENDIF
!
!---    Plem and Xiu ground-surface solar angle model  ---
!
        IF( IALB(IZ(N)).EQ.1 ) THEN
          ALB_S = ALB_S + 1.D-2*(EXP(3.286D-3*(DZ_ANG**1.5D+0))-1.D+0)
!
!---    Briegleb ground-surface solar angle model  ---
!
        ELSEIF( IALB(IZ(N)).EQ.2 ) THEN
          Z_ANG = DZ_ANG*GPI/1.8D+2
          ALB_S = ALB_S + ALBEDO(4,IZ(N))*
     &      (1.D+0 + ALBEDO(5,IZ(N)))/
     &      (1.D+0 + 2.D+0*ALBEDO(5,IZ(N))*COS(Z_ANG))
!
!---    Wang ground-surface solar angle model  ---
!
        ELSEIF( IALB(IZ(N)).EQ.3 ) THEN
          Z_ANG = DZ_ANG*GPI/1.8D+2
          G_1 = (-7.574D-3) + (-7.0987D-2)*(Z_ANG**2) + 
     &      (3.07588D-1)*(Z_ANG**3)
          G_2 = (-1.284909D+0) + (-1.66314D-1)*(Z_ANG**2) +
     &      (4.184D-2)*(Z_ANG**3)
          G_1R = 2.67808D-1
          G_2R = -1.419244D+0
          ALB_S = ALB_S + ALBEDO(4,IZ(N))*(1.D+0 +
     &      0.346D+0*(G_1-G_1R) + 0.063*(G_2-G_2R))
        ELSE
          ALB_S = 0.25D+0
        ENDIF
        RN_SNS = RN*(1.D+0-ALB_S)
!
!---    Gas tortuosity at the ground surface  ---
!
        CALL TORTU( IZ(N),SL_S,SG_S,ZERO,PORD(2,N),
     &    TORL_S,TORG_S,TORN_S )
!
!---    Water-vapor diffusion coefficient at the ground surface  ---
!
        IF( ISLC(2).EQ.1 ) THEN
          DFGW_S = DFGWC
        ELSEIF( ISLC(2).EQ.2 ) THEN
          CALL BNDFAW( T_SY,PG_S,DFGW_S )
        ELSEIF( ISLC(2).EQ.3 ) THEN
          CALL BNDFAW( T_SY,PG_S,DFGW_S )
          CMFF = 1.D+0 + 2.6D+0/(DFGWC**0.5)
          AMC = PORD(2,N)*SL_S
          ENHF = 9.5D+0 + 6.D+0*(AMC) -
     &      8.5D+0/EXP((CMFF*AMC)**4)
          DFGW_S = ENHF*DFGW_S
        ELSEIF( ISLC(2).EQ.4 ) THEN
          CALL BNDFAW( T_SY,PG_S,DFGW_S )
          ENHF = DFEF(1,IZ(N))+DFEF(2,IZ(N))*SL_S-
     &      (DFEF(1,IZ(N))-DFEF(4,IZ(N)))
     &      *EXP(-((DFEF(3,IZ(N))*SL_S)**DFEF(5,IZ(N))))
          DFGW_S = ENHF*DFGW_S
        ENDIF
!
!---    Effective thermal conductivity coefficient  ---
!
        CALL WATLQK( T_SY,PG_S,PSW_S,THKL_S )
        INDX = IBD
        TK_N = THKL(M,N)
        TK_S = THKE_L( IZ(N),SL_S,THKL_S,PORD(2,N),PORT(2,N),INDX )
        INDX = -1
        TK_NS = DIFMN(TK_N,TK_S,D_NS,D_NS,ZERO,INDX)
!
!---    Sensible heat flux from the field node to
!       the ground surface, W/m^2  ---
!
        H_NS = (T(M,N)-T_SX)*TK_NS/D_NS
!
!---    Successive substitution, following Campbell [1985],
!       to determine aerodynamic resistance  ---
!
        H_SA = (T_SX-TA)*RHOG_A*CPG_A/RA_BSX
        VARX = K_VK*ATMC(2)*GRAV/(RHOG_A*CPG_A*TA_K)
        ASP = -VARX*H_SA/(UZ_FX**3)
        DO 330 NCX = 1,4
!
!---      Surface temperature higher than atmospheric
!         temperature, unstable conditions  ---
!
          IF( (T_SX-TA).GT.EPSL ) THEN
            SPCFH = -2.D+0*LOG(5.D-1*(1.D+0+SQRT(
     &        MAX(1.D+0-(1.6D+1*ASP),0.D+0))))
            SPCFM = 6.D-1*SPCFH
!
!---      Surface temperature lower than atmospheric
!         temperature, stable conditions  ---
!
          ELSEIF( (T_SX-TA).LT.-EPSL ) THEN
            SPCFH = 4.7D+0*ASP
            SPCFM = 4.7D+0*ASP
!
!---      Surface temperature equal to atmospheric
!         temperature, skip calculations  ---
!
          ELSE
            SPCFH = 0.D+0
            SPCFM = 0.D+0
          ENDIF
!
!---      Corrections limited to 20%  ---
!
          SPCFH = SIGN( MIN( ABS(SPCFH),2.D-1*VLG_2 ),SPCFH )
          SPCFM = SIGN( MIN( ABS(SPCFM),2.D-1*VLG_1 ),SPCFM )
          UZ_F = K_VK*UZ/(VLG_1+SPCFM)
          RA_BS = MIN( (VLG_2+SPCFH)/(K_VK*UZ_F),RA_BSMX )
!
!---      Sensible heat flux from the ground surface to the
!         atmosphere, W/m^2  ---
!
          H_SA = (T_SX-TA)*RHOG_A*CPG_A/RA_BS
!
!---      Atmospheric stability parameter  ---
!
          ASP = -VARX*H_SA/(UZ_F**3)
  330   CONTINUE
!
!---    Aqueous volumetric flux from the node to the
!       ground surface, m/s  ---
!
        HDL_S = PL_SX + D_NS*GBX*5.D-1*(RHOL(M,N)+RHOL_S)
        HDL_N = PL(M,N) + PATM
        HDL_NS = HDL_N - HDL_S
        INDX = -1
        VISL_NS = DIFMN(VISL(M,N),VISL_S,D_NS,D_NS,HDL_NS,INDX)
        INDX = -4
        RKL_NS = DIFMN(RKL(IBD,M,N),RKL_S(IBD),D_NS,D_NS,HDL_NS,INDX)
        L_NS = (RKL_NS*PERM(IBD,IZ(N))/VISL_NS)*(HDL_NS/D_NS)
!
!---    Diffusive water-vapor mass flux from the ground
!       surface to the atmosphere, kg/m^2 s  ---
!
        E_SA = (E_S-E_A)/RA_BS
        PE_SA = (PE_S-E_A)/RA_BS
!
!---    Gas volumetric flux from the node to the ground
!       surface, m/s  ---
!
        HDG_S = PG_S + D_NS*GBX*5.D-1*(RHOG(M,N)+RHOG_S)
        IF( 1.D+0-SL(M,N).LT.1.D-9 ) THEN
          HDG_N = HDG_S
        ELSE
          HDG_N = PG(M,N) + PATM
        ENDIF
        HDG_NS = HDG_N - HDG_S
        INDX = -1
        VISG_NS = DIFMN(VISG(M,N),VISG_S,D_NS,D_NS,HDG_NS,INDX)
        INDX = -4
        RKG_NS = DIFMN(RKG(M,N),RKG_S,D_NS,D_NS,HDG_NS,INDX)
!
!---    Gas flow from node to ground surface, m^3/s  ---
!
        G_NS = (RKG_NS*PERM(IBD,IZ(N))/VISG_NS)*((HDG_N-HDG_S)/(D_NS))
!
!---    Diffusive water-vapor mass flux from the node to
!       the surface, kg/m^2 s  ---
!
        DF_S = TORG_S*PORD(2,N)*SG_S*DFGW_S
        DF_N = TORG(M,N)*PORD(M,N)*(SG(M,N)-SGT(M,N))*DFGW(M,N)
        INDX = -1
        DF_NS = DIFMN(DF_N,DF_S,D_NS,D_NS,G_NS,INDX)
        E_N = RHOG(M,N)*XGW(M,N)
        E_NS = DF_NS*(E_N-E_S)/D_NS
!
!---    Enthalpy of water-vapor, air and gas at the surface, J/kg  ---
!
        CALL WATGSH( T_SX,PSW_S,RHOGW_S,HGW_S,UGW_S )
        CALL AIRGSH( T_SX,PGA_S,HGA_S,UGA_S )
        HG_S = XGW_S*HGW_S + XGA_S*HGA_S
        UG_S = XGW_S*UGW_S + XGA_S*UGA_S
!
!---    Enthalpy of liquid-water and aqueous at the
!       ground surface, J/kg  ---
!
        CALL WATLQH( T_SX,PGA_S,HLW_S )
        HL_S = XLW_S*HLW_S + XLA_S*HGA_S
!
!---    Enthalpy of gas at the node, J/kg  ---
!
        HG_N = XGW(M,N)*HGW(M,N) + XGA(M,N)*HGA(M,N)
!
!---    Upward long-wave radiation from the ground surface  ---
!
        TS_K = T_SX + TABS
        EM_S = 9.D-1 + 1.8D-1*(SL_S*PORD(2,N))
        RN_LUS = (1.D+0-EM_S)*RN_LDS + EM_S*SIGMA*(TS_K**4)
!
!---    Net long-wave radiation into the ground surface  ---
!
        RN_LNS = RN_LDS - RN_LUS
!
!---    Net long- and short-wave radiation into the
!       ground surface  ---
!
        RN_NS = RN_LNS + RN_SNS
!
!---    Component and phase density interfacial averages  ---
!
        RHOGA_N = RHOG(M,N)*XGA(M,N)
        RHOGW_N = RHOG(M,N)*XGW(M,N)
        RHOLW_N = RHOL(M,N)*XLW(M,N)
        RHOLW_S = RHOL_S*XLW_S
        RHOLW_A = RHOL_A*XLW_A
        INDX = -4
        RHOGW_SA = DIFMN(RHOGW_S,RHOGW_A,D_NS,D_NS,G_NS,INDX)
        RHOGW_NS = DIFMN(RHOGW_N,RHOGW_S,D_NS,D_NS,G_NS,INDX)
        RHOGA_SA = DIFMN(RHOGA_S,RHOGA_A,D_NS,D_NS,G_NS,INDX)
        RHOGA_NS = DIFMN(RHOGA_N,RHOGA_S,D_NS,D_NS,G_NS,INDX)
        RHOG_SA = DIFMN(RHOG_S,RHOG_A,D_NS,D_NS,G_NS,INDX)
        RHOG_NS = DIFMN(RHOG(M,N),RHOG_S,D_NS,D_NS,G_NS,INDX)
        RHOLW_NS = DIFMN(RHOLW_N,RHOLW_S,D_NS,D_NS,L_NS,INDX)
        RHOLW_SA = DIFMN(RHOLW_S,RHOLW_A,D_NS,D_NS,L_SA,INDX)
        RHOLW_PA = DIFMN(RHOLW_S,RHOLW_A,D_NS,D_NS,L_SA,INDX)
!
!---    Gas flow from ground surface to atmosphere, m^3/s
!       using air mass balance at ground surface  ---
!
        G_SA = G_NS*RHOGA_NS/RHOGA_SA
!
!---    Component and phase interfacial averages  ---
!
        INDX = -4
        HGW_SA = DIFMN(HGW_S,HGW_A,D_NS,D_NS,E_SA,INDX)
        HGW_NS = DIFMN(HGW(M,N),HGW_S,D_NS,D_NS,E_NS,INDX)
        HGA_SA = DIFMN(HGA_S,HGA_A,D_NS,D_NS,E_SA,INDX)
        HGA_NS = DIFMN(HGA(M,N),HGA_S,D_NS,D_NS,E_NS,INDX)
        HG_SA = DIFMN(HG_S,HG_A,D_NS,D_NS,G_SA,INDX)
        HG_NS = DIFMN(HG_N,HG_S,D_NS,D_NS,G_NS,INDX)
        HL_NS = DIFMN(HL(M,N),HL_S,D_NS,D_NS,L_NS,INDX)
        HL_SA = DIFMN(HL_S,HL_A,D_NS,D_NS,L_SA,INDX)
!
!---    Ponding height limit  ---
!
!        L_SAMX = -PERM(3,IZ(N))*(RHOLW_SA**2)*GRAV*
!     &    (1.D+0-SWBCD/DB)/VISL_S
        HDL_SX = PG_S + SWBCD*RHOL_A*GRAV + 
     &    D_NS*GBX*5.D-1*(RHOL(M,N)+RHOL_S)
        HDL_NX = PL(M,N) + PATM
        HDL_NSX = HDL_NX - HDL_SX
        L_SAMX = MIN((PERM(3,IZ(N))/VISL_S)*(HDL_NSX/D_NS),0.D+0 )
        L_SAX = L_SA*RHOLW_SA - RWRO(1)*RHOLW_SA/A_SA
        IF( ABS(L_SAX).GT.ABS(L_SAMX) ) THEN
          L_SAX = L_SAMX
        ENDIF
!
!---    Residual of the water balance at the ground surface  ---
!
        RW_S(MP) = E_NS + G_NS*RHOGW_NS + L_NS*RHOLW_NS
     &    - E_SA - G_SA*RHOGW_SA - L_SAX
!
!---    Residual of the energy balance at the ground surface,
!       ignoring dissolved air mass in aqueous  ---
!
        RE_S(MP) = RN_NS + E_NS*HGW_NS + G_NS*RHOG_NS*HG_NS
     &    + L_NS*RHOLW_NS*HL_NS + H_NS
     &    - E_SA*HGW_SA - G_SA*RHOG_SA*HG_SA
     &    - L_SAX*HL_SA - H_SA
  400 CONTINUE
      IF( ABS(MPX).EQ.1 ) GOTO 430
!
!---  Jacobian coefficient matrix  ---
!
      IF( MPX.GT.0 ) THEN
        BJ(1) = -RE_S(1)
        BJ(2) = -RW_S(1)
        DO 410 JP = 1,ABS(MPX)-1
          AJ(1,JP) = (RE_S(JP+1)-RE_S(1))/DNRX(JP)
          AJ(2,JP) = (RW_S(JP+1)-RW_S(1))/DNRX(JP)
  410   CONTINUE
      ELSE
        BJ(1) = -RW_S(1)
        BJ(2) = -RE_S(1)
        DO 412 JP = 1,ABS(MPX)-1
          AJ(1,JP) = (RW_S(JP+1)-RW_S(1))/DNRX(JP)
          AJ(2,JP) = (RE_S(JP+1)-RE_S(1))/DNRX(JP)
  412   CONTINUE
      ENDIF
!
!---  Solve linear system  ---
!
      JP = ABS(MPX)-1
      CALL LUDCMP( AJ,JP,KP,IJ,DJ )
      CALL LUBKSB( AJ,JP,KP,IJ,BJ )
!
!---  Limit changes in surface saturation  ---
!
      IF( MPX.EQ.3 ) THEN
        DP_S = MAX( 1.D+5,ABS(PG_S-PL_S) )
        DP_S = MIN( DP_S,ABS(BJ(2)) )
        DP_S = SIGN( DP_S,BJ(2) )
      ELSEIF( MPX.EQ.-2 ) THEN
        DP_S = MAX( 1.D+5,ABS(PG_S-PL_S) )
        DP_S = MIN( DP_S,ABS(BJ(1)) )
        DP_S = SIGN( DP_S,BJ(1) )
      ELSE
        DP_S = 0.D+0
      ENDIF
!
!---  Surface saturation above residual  ---
!
      IF( SL_SX.GT.SCHR(4,IZ(N)) ) THEN
!
!---    Limit rising saturation changes to 0.05  ---
!
        IF( DP_S.GT.0.D+0 .AND. SL_SX.LT.9.5D-1 ) THEN
          CALL CAP_GT( IZ(N),SL_SX,CPGLX )
          SLX = SL_SX+5.D-2
          CALL CAP_GT( IZ(N),SLX,CPGL )
          DP_S = MIN( DP_S,(CPGLX-CPGL) )
!
!---    Limit falling saturation changes to 0.05  ---
!
        ELSEIF( DP_S.LT.0.D+0 ) THEN
          CALL CAP_GT( IZ(N),SL_SX,CPGLX )
          SLX = MAX( SL_SX-5.D-2,9.5D-1*SCHR(4,IZ(N)),5.D-1*SL_SX )
          CALL CAP_GT( IZ(N),SLX,CPGL )
!
!---      Transition from saturated conditions  ---
!
          CPGL = MAX( PL_S-PG_S+CPGL,CPGL )
          DP_S = MAX( DP_S,(CPGLX-CPGL) )
        ENDIF
!
!---  Surface saturation below residual  ---
!
      ELSE
!
!---    Limit rising saturation changes to 0.05
!       of the residual saturation  ---
!
        IF( DP_S.GT.0.D+0 ) THEN
          CALL CAP_GT( IZ(N),SL_SX,CPGLX )
          SLX = SL_SX + 5.D-2*SCHR(4,IZ(N))
          CALL CAP_GT( IZ(N),SLX,CPGL )
          DP_S = MIN( DP_S,(CPGLX-CPGL) )
!
!---    Limit falling saturation changes to 0.05
!       of the residual saturation  ---
!
        ELSEIF( DP_S.LT.0.D+0 ) THEN
          CALL CAP_GT( IZ(N),SL_SX,CPGLX )
          SLX = MAX( SL_SX-5.D-2*SCHR(4,IZ(N)),5.D-1*SL_SX )
          CALL CAP_GT( IZ(N),SLX,CPGL )
          DP_S = MAX( DP_S,(CPGLX-CPGL) )
        ENDIF
      ENDIF
!
!---  Limit changes in temperature at the ground surface  ---
!
      IF( MPX.GT.0 ) THEN
        DT_S = MIN( 5.D-1,ABS(BJ(1)) )
        DT_S = SIGN( DT_S,BJ(1) )
      ELSE
        DT_S = 0.D+0
      ENDIF
!
!---  Relax solution for flopping surface pressure conditions
!     or transitions  ---
!
      IFLOP = 0
      IF( (DP_SO*DP_S).LT.0.D+0 .AND.
     &  ABS( 1.D+0-ABS(DP_SO/(DP_S+EPSL)) ).LT.1.D-1 .AND.
     &  (DP_SOO*DP_SO).LT.0.D+0 .AND.
     &  ABS( 1.D+0-ABS(DP_SOO/(DP_SO+EPSL)) ).LT.1.D-1) THEN
        IFLOP = 1
      ENDIF
      IF( IFLOP.EQ.1 ) THEN
        PL_S = PL_S + 6.D-1*DP_S
        T_S = T_S + 6.D-1*DT_S
      ELSE
        PL_S = PL_S + DP_S
        T_S = T_S + DT_S
      ENDIF
      PL_S = MIN( PL_S,PG_S+SWBCD*RHOL_A*GRAV )
!
!---  Store correction to surface pressure and
!     surface capillary pressure  ---
!
      DP_SOO = DP_SO
      DP_SO = DP_S
!
!---  Maximum residual  ---
!
      RSDWX = 1.D-20
      RSDEX = 1.D-20
      IF( MPX.GT.0 ) RSDEX = ABS(BJ(1))/TABS
      IF( MPX.EQ.3 ) RSDWX = ABS(BJ(2))/PATM
      IF( MPX.EQ.-2 ) THEN
        CALL KSP_GT( N,IZ(N),M,PG_S,PL_S,SLX,SG_S,RKL_S,RKG_S )
        RSDWX = MIN( ABS(BJ(1))/PATM,ABS(SL_S-SLX)*1.D-2 )
      ENDIF
      RSDX = MAX( RSDWX,RSDEX )
      IF( NC.GT.24 ) THEN
        T_SMN = T_SMN + T_S
        PL_SMN = PL_SMN + PL_S
        SL_SMN = SL_SMN + SL_S
      ENDIF
!
!---  No convergence on water mass and energy
!     balances ---
!
      IF( NC.EQ.48 ) THEN
        REALX = REAL(NC-24)
        SL_S = SL_SMN/REALX
!
!---    Low residual, select the average temperature
!       and aqueous pressure ---
!
        IF( RSDX.LT.1.D-4 ) THEN
          T_S = T_SMN/REALX
          PL_S = PL_SMN/REALX
          CALL CAP_GT( IZ(N),SL_S,CPGLX )
          PL_S = PG_S - CPGLX
          PL_S = MIN( PL_S,PG_S+SWBCD*RHOL_A*GRAV )
          GOTO 420
!
!---    Fix surface aqueous pressure ---
!
        ELSEIF( SL_S.LT.5.D-1 ) THEN
          IF( M.EQ.2 ) THEN
            PL_S = PL_SMN/REALX
            CALL CAP_GT( IZ(N),SL_S,CPGLX )
            PL_S = PG_S - CPGLX
            PL_S = MIN( PL_S,PG_S+SWBCD*RHOL_A*GRAV )
          ENDIF
          T_S = T_SMN/REALX
          MPX = 2
!
!---    Fix surface temperature ---
!
        ELSE
          IF( M.EQ.2 ) THEN
            T_S = T_SMN/REALX
          ENDIF
          PL_S = PL_SMN/REALX
          CALL CAP_GT( IZ(N),SL_S,CPGLX )
          PL_S = PG_S - CPGLX
          PL_S = MIN( PL_S,PG_S+SWBCD*RHOL_A*GRAV )
          MPX = -2
!!
!!---      Convergence failure cut time step ---
!!
!          ISWBCF = 1
!          ISUB_LOG = ISUB_LOG-1
!          RETURN
        ENDIF
      ENDIF
      IF( NC.EQ.72 ) THEN
        WRITE(ISC,'(10X,A)') '---  Shuttleworth-Wallace Boundary ' //
     &    'Condition Failure: Bare Surface: Used Averaged Values  ---'
        WRITE(IWR,'(10X,A)') '---  Shuttleworth-Wallace Boundary ' //
     &    'Condition Failure: Bare Surface: Used Averaged Values  ---'
        T_S = T_SMN/REALX
        PL_S = PL_SMN/REALX
        CALL CAP_GT( IZ(N),SL_S,CPGLX )
        PL_S = PG_S - CPGLX
        PL_S = MIN( PL_S,PG_S+SWBCD*RHOL_A*GRAV )
        GOTO 420
      ENDIF
!
!---  Check solution convergence  ---
!
      IF( RSDX.GE.1.D-6 ) GOTO 310
  420 CONTINUE
!
!---  Convergence on mass and energy balances at the
!     ground surface  ---
!
      T_SO(2,NB) = T_S
      PL_SO(2,NB) = PL_S
!
!---  Converged solution, calculate fluxes  ---
!
      MPOX = MPX
      MPX = 1
      NC = NC-1
      GOTO 310
  430 CONTINUE
      MPX = MPOX
      IF( M.EQ.2 ) THEN
!
!---    Load variables for output  ---
!
!       RHON(1,N) - Water mass balance at ground surface, kg/s
!       RHON(2,N) - Actual evaporation rate, kg/s
!       RHON(3,N) - Potential evaporation rate, kg/s
!       RHON(4,N) - Actual transpiration rate, kg/s
!       RHON(5,N) - Potential transpiration rate, kg/s
!
        RHON(1,N) = (E_NS + G_NS*RHOGW_NS + L_NS*RHOLW_NS)*A_SA
        RHON(2,N) = (E_SA + G_SA*RHOGW_SA)*A_SA
        RHON(3,N) = (PE_SA + G_SA*RHOGW_SA)*A_SA
        RHON(4,N) = 0.D+0
        RHON(5,N) = 0.D+0
!
!---    Load variables for output
!
!       SN(1,N) - Surface-atmosphere latent heat flux, W/m^2
!       SN(2,N) - Surface-atmosphere sensible heat flux, W/m^2
!       SN(3,N) - Surface-atmosphere net long-wave radiation, W/m^2
!       SN(4,N) - Surface-atmosphere net short-wave radiation, W/m^2
!       SN(5,N) - Ground-surface sensible and latent heat flux, W/m^2
!
        SN(1,N) = E_SA*(HGW_SA-HL_S)
     &    + L_SAX*(HL_SA-HL_S)
     &    + G_SA*RHOGW_SA*(HGW_SA-HL_S)
        SN(2,N) = G_SA*RHOGA_SA*(HGA_SA-HGA_S) + H_SA
        SN(3,N) = RN_LNS
        SN(4,N) = RN_SNS
        SN(5,N) = E_NS*(HGW_NS-HL_S)
     &    + L_NS*RHOLW_NS*(HL_NS-HL_S)
     &    + G_NS*RHOGW_NS*(HGW_NS-HL_S)
     &    + G_NS*RHOGA_NS*(HGA_NS-HGA_S) + H_NS
!
!---    Load variables for output  ---
!
!       XLO(1,N) - Surface temperature, C
!       XLO(2,N) - Surface aqueous pressure, Pa
!       XLO(3,N) - Surface gas pressure, Pa
!       XLO(4,N) - Surface aqueous saturation
!       XLO(5,N) - Surface water-vapor pressure, Pa
!
        XLO(1,N) = T_S
        XLO(2,N) = PL_S
        XLO(3,N) = PG_S
        XLO(4,N) = SL_S
        XLO(5,N) = PGW_S
!
!---    Load variables for output  ---
!
!       PN(1,N) - Atmospheric temperature, C
!       PN(2,N) - Atmospheric pressure, Pa
!       PN(3,N) - Atmospheric relative humidity
!       PN(4,N) - Atmospheric wind speed, m/s
!       PN(5,N) - Surface solar radiation, W/m^2
!
        PN(1,N) = TA
        PN(2,N) = PA
        PN(3,N) = RH
        PN(4,N) = UZ
        PN(5,N) = RN
!
!---    Load variables for output  ---
!
!       XGO(1,N) - Rain-water volumetric runoff rate, m^3/s
!       XGO(2,N) - Aerodynamic resistance between ground and
!                   reference height for bare surface, s/m
!       XGO(3,N) - Surface volumetric precipitation rate, m^3/s
!       XGO(4,N) - Surface mass precipitation rate, kg/s
!       XGO(5,N) - Atmospheric mass precipitation rate, kg/s
!!       XGO(5,N) - Atmospheric stability parameter
!
        RWROX = (L_SAX/RHOLW_SA - L_SA)*A_SA
        RWRO(2) = RWRO(1) + RWROX
        IF( ABS(RWRO(2)).LT.1.D-20 ) RWRO(2) = 0.D+0
        XGO(1,N) = RWROX
        XGO(2,N) = RA_BS
        XGO(3,N) = (L_SAX/RHOLW_SA)*A_SA
        XGO(4,N) = L_SAX*A_SA
        XGO(5,N) = -BCX(2)*RHOLW_SA*A_SA
!        XGO(5,N) = -K_VK*ATMC(2)*GRAV*H_NS/
!     &    (RHOG_A*CPG_A*TA_K*(UZ_F**3))
!
!---    Surface area  ---
!
        IF( IBCD(NB).EQ.-3 ) THEN
          A_SA = AFZ(NSZ(N))
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
          A_SA = AFY(NSY(N))
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
          A_SA = AFX(NSX(N))
        ELSEIF( IBCD(NB).EQ.1 ) THEN
          A_SA = AFX(NSX(N)+1)
        ELSEIF( IBCD(NB).EQ.2 ) THEN
          A_SA = AFY(NSY(N)+IFLD)
        ELSEIF( IBCD(NB).EQ.3 ) THEN
          A_SA = AFZ(NSZ(N)+IJFLD)
        ENDIF
      ENDIF
!
!---  Load boundary condition variables
!
      TB(M,NB) = T_S
      PLB(M,NB) = PL_S - PATM
      PGB(M,NB) = PG_S - PATM
      PVWB(M,NB) = PGW_S
      PVAB(M,NB) = PGA_S
      SLB(M,NB) = SL_S
      RKLB(1,M,NB) = RKL_S(1)
      RKLB(2,M,NB) = RKL_S(2)
      RKLB(3,M,NB) = RKL_S(3)
      SGB(M,NB) = SG_S
      RKGB(M,NB) = RKG_S
      XLWB(M,NB) = XLW_S
      XLAB(M,NB) = XLA_S
      XMLWB(M,NB) = XMLW_S
      XMLAB(M,NB) = XMLA_S
      XGWB(M,NB) = XGW_S
      XGAB(M,NB) = XGA_S
      XMGWB(M,NB) = XMGW_S
      XMGAB(M,NB) = XMGA_S
      RHOLB(M,NB) = RHOL_S
      RHOGB(M,NB) = RHOG_S
      HGWB(M,NB) = HGW_S
      HGAB(M,NB) = HGA_S
      HLB(M,NB) = HL_S
      UEGB(M,NB) = UG_S
!
!---  Impose boundary conditions through sources
!     at the boundary node.
!
!     SRCT(M,N) - thermal energy rate into the node, W
!     SRCA(M,N) - air mass rate into the node, kg/s
!     SRCW(M,N) - water mass rate into the node, kg/s
!
      SRCT(M,N) = SRCT(M,N)-(E_NS*HGW_NS + G_NS*RHOG_NS*HG_NS
     &    + L_NS*RHOLW_NS*HL_NS + H_NS)*A_SA
      SRCA(M,N) = SRCA(M,N)-(G_NS*RHOGA_NS)*A_SA
      SRCW(M,N) = SRCW(M,N)-(E_NS + G_NS*RHOGW_NS + L_NS*RHOLW_NS)*A_SA
!
!---  Archive unknowns  ---
!
      IF( M.EQ.2 ) THEN
        T_SO(2,NB) = T_S
        PL_SO(2,NB) = PL_S
      ENDIF
 1000 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BARE_S group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BC_SW( BCX,DB,GB,M,N,NB )
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
!     Shuttleworth-Wallace energy boundary condition.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 August 2002.
!     Last Modified by MD White, PNNL, 12 August 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PLT_ATM
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 BCX(LBCV)
      REAL*8 PAIX,LAIX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BC_SW'
      IPLANT = 0
      DO 100 IP = 1,NPLANT
        IPX = (IP-1)*2 + 4
        PAIX = BCX(IPX)
        LAIX = BCX(IPX-1)
        IF( PAIX*LAIX.GT.0.D+0 ) IPLANT = 1
  100 CONTINUE
!
!---  Bare soil surface, no plants  ---
!
      IF( NPLANT.EQ.0 .OR. IPLANT.EQ.0 .OR. IBCT(IEQT,NB).EQ.31 ) THEN
        CALL BARE_S( BCX,DB,GB,M,N,NB )
!
!---  Single plant temperature option, solve for surface temperature
!     and surface aqueous pressure, plant temperature, air temperature
!     at the mean canopy flow height, and water-vapor partial pressure
!     at the mean canopy flow height  ---
!
      ELSEIF( MOD(ISLC(24),10).EQ.1 ) THEN
        CALL SNGL_T( BCX,DB,GB,M,N,NB )
!
!---  Multiple plant temperature option, solve for surface temperature
!     and water-vapor density, plant specie temperature, specie air
!     temperature at the mean canopy flow height, and specie water-vapor
!      density at the mean canopy flow height  ---
!
!      ELSEIF( MOD(ISLC(24),10).EQ.2 ) THEN
!        CALL MULT_T( BCX,DB,GB,M,N,NB )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BC_SW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CLOUD( RN,F_CC,J_DAY,DZ_ANG )
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
!     Fractional cloud cover.
!
!     Llasat, M.C., and R.L. Snyder. 1998.  "Data error effects on net
!     radiation and evapotranspiration estimation."  Agricultural
!     and Forest Meterology, 91:209-221.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 31 October 2003.
!     Last Modified by MD White, PNNL, 31 October 2003.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PLT_ATM
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Data Statements---------------------------------!
!
      SAVE G_SC
      DATA G_SC / 1.367D+3 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CLOUD'
!
!---  Incident solar radiation zero, skip calculations  ---
!
      IF( RN/EPSL.LT.EPSL ) GOTO 100
!
!---  Equation of time, hr  ---
!
!      B_EQT = 2.D+0*GPI*REAL(J_DAY-81)/364.D+0
!      T_EQT = 1.645D-1*SIN(2.D+0*B_EQT) - 1.255D-1*COS(B_EQT) -
!     &  2.5D-2*SIN(B_EQT)
      IF( J_DAY.LT.181 ) THEN
        EJ_DAY = J_DAY*1.D-2
        T_EQT = -0.04056D+0 - 0.74503D+0*EJ_DAY + 0.08823D+0*(EJ_DAY**2)
     &    + 2.0516D+0*(EJ_DAY**3) - 1.8111D+0*(EJ_DAY**4)
     &    + 0.42832D+0*(EJ_DAY**5)
      ELSE
        EJ_DAY = (J_DAY-180)*1.D-2
        T_EQT = -0.05039D+0 - 0.33954D+0*EJ_DAY + 0.04084D+0*(EJ_DAY**2)
     &    + 1.8928D+0*(EJ_DAY**3) - 1.7619D+0*(EJ_DAY**4)
     &    + 0.4224D+0*(EJ_DAY**5)
      ENDIF
!
!---  Time of solar noon, hr
!     ATMC(3) longitude, radians
!     ATMC(4) latitude, radians
!     ATMC(5) meridian, radians  ---
!
      T_SN = 1.2D+1 + (ATMC(3)-ATMC(5))*1.2D+1/GPI - T_EQT
!
!---  Time of day, hr  ---
!
      T_DAY = MOD(TM+ATMST,8.64D+4)/3.6D+3
!
!---  Solar time angle, radians  ---
!
      S_TANG = (T_DAY-T_SN)*GPI/1.2D+1
!
!---  Solar declination, radians  ---
!
      REALX = REAL(J_DAY-172)
      S_DEC = 4.101D-1*COS(2.D+0*GPI*REALX/3.65D+2)
      SJ_DAY = J_DAY*1.D-2
      S_DEC = -0.37726D+0 - 0.10564D+0*SJ_DAY + 1.2458D+0*(SJ_DAY**2)
     &  - 0.75478*(SJ_DAY**3) + 0.13627*(SJ_DAY**4)
     &  - 0.00572*(SJ_DAY**5)
      S_DEC = ASIN(S_DEC)
!
!---  Sun altitude, radians  ---
!
      S_ANG = ASIN( SIN(S_DEC)*SIN(ATMC(4)) +
     &  COS(S_DEC)*COS(ATMC(4))*COS(S_TANG) )
!
!---  Sun altitude, deg  ---
!
      DS_ANG = S_ANG*1.8D+2/GPI
!
!---  Solar zenith angle, deg  ---
!
      DZ_ANG = MIN( MAX( 9.D+1-DS_ANG,0.D+0 ),9.D+1 )
!
!---  Sun angle less than 10 deg, skip calculations ---
!
      IF( DS_ANG.LT.1.D+1 ) GOTO 100
!
!---  Incident clear sky irradiance ---
!     RN_EX - extra terrestrial radiation on a horizontal surface, W/m^2
!
      RN_EX = 1367.D+0*SIN(S_ANG)
      RN_A = (0.79D+0 - 3.75D+0/DS_ANG)*RN_EX
      F_CC = 1.088D+0*(MAX( (1.D+0 - (RN/RN_A)),0.D+0 )**(0.294D+0))
      F_CC = MIN( F_CC,1.D+0 )
      F_CC = MAX( F_CC,0.D+0 )
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CLOUD group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DEW_PT( RH,TA,T_DP )
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
!     Dew point
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 31 October 2003.
!     Last Modified by MD White, PNNL, 31 October 2003.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 PSX(2),FX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DEW_PT'
!
!---  Guess dew-point temperature  ---
!
      T_DP = TA
      DT_DP = -1.D-4
      INDX = 0
      CALL REGION_4( T_DP,PSAX,INDX )
!
!---  Newton iteration to compute dew-point temperature  ---
!
      NC = 0
  100 CONTINUE
      NC = NC + 1
      DO 200 M = 1,2
        T_DPX = T_DP
        IF( M.EQ.2 ) T_DPX = T_DP + DT_DP
        INDX = 0
        CALL REGION_4( T_DPX,PSX(M),INDX )
        FX(M) = PSAX*RH - PSX(M)
  200 CONTINUE
      DFX = (FX(2)-FX(1))/DT_DP
      DT_DPX = -FX(1)/DFX
      T_DP = T_DP + DT_DPX
!
!---  No convergence on phase compositions  ---
!
      IF( NC.GT.32 ) THEN
        NC = NC
      ENDIF
      IF( NC.GT.64 ) THEN
        INDX = 14
        CHMSG = 'Unconverged Dew Point Temperature'
        RLMSG = T_DP
        CALL WRMSGS( INDX )
      ENDIF
      IF( ABS(DT_DPX).GT.1.E-6 ) GOTO 100
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DEW_PT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE POLINT( XA,YA,N,X,Y,DY,IERR )
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
!     Given arrays XA and YA, each of length N, and given a value X,
!     this routine returns a value Y, and an error estimate DY.
!     If P(X) is the polynomial of degree N-1 such that
!     P(XA(I)) = YA(I),
!     I = 1,N, then the returned value Y = P(X).
!
!     Press, W.H., B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling.
!     1986.  Numerical Recipes, The Art of Scientific Computing.
!     Cambridge University Press, Cambridge.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, February, 1999.
!     Last Modified by Mark White, Battelle, February 19, 1999.
!     roots.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Parameter Statements----------------------------!
!
      PARAMETER( NMAX=10, EPSL=1.D-14 )
      REAL*8 XA(N),YA(N),C(NMAX),D(NMAX)
!
!----------------------Executable Lines--------------------------------!
!
      IERR = 0
      NS = 1
      DIF = ABS(X-XA(1))
!
!---  Find the index NS of the closest table entry  ---
!
      DO 11 I = 1,N
        DIFT = ABS(X-XA(I))
        IF( DIFT.LT.DIF ) THEN
          NS = I
          DIF = DIFT
        ENDIF
!
!---  Initialize the tableau fo C's and D's  ---
!
        C(I) = YA(I)
        D(I) = YA(I)
   11 CONTINUE
!
!---  Initial approximation to Y  ---
!
      Y = YA(NS)
      NS = NS-1
!
!---  Loop over the columns in the tableau and
!     update the C's and D's  ---
!
      DO 13 M = 1,N-1
        DO 12 I = 1,N-M
          HO = XA(I)-X
          HP = XA(I+M)-X
          W = C(I+1)-D(I)
          DEN = HO-HP
!
!---  Exit subroutine if two input XA's are
!     identical, within roundoff  ---
!
          IF( ABS( DEN )/EPSL.LT.EPSL ) THEN
            IERR = 1
            RETURN
          ENDIF
          DEN = W/DEN
!
!---  Update C's and D's  ---
!
          D(I) = HP*DEN
          C(I) = HO*DEN
   12   CONTINUE
!
!---  After each column in the tableau is completed, decide
!     which direction C or D to add the accumulating value
!     of Y, (i.e., which path to take through the tableau -
!     forking up or down).  Do this in such a way as
!     to take the most "straight line" route through the
!     tableau to its apex, updating NS accordingly to keep track.
!     This route keeps the partial approximations centered
!     (insofar as possible) on the target X.  The
!     last DY added is thus the error indication.  ---
!
        IF( 2*NS.LT.N-M ) THEN
          DY = C(NS+1)
        ELSE
          DY = D(NS)
          NS = NS-1
        ENDIF
        Y = Y+DY
   13 CONTINUE
!
!---  End of POLINT group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE QROMB( FUNC,A,B,SSX,IERR,INDX )
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
!     Returns as S the integral of the function FUNC from A to B.
!     Integration is performed by Romberg's method of order 2K, where
!     e.g., K=2 is Simpson's rule.
!
!     Press, W.H., B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling.
!     1986.  Numerical Recipes, The Art of Scientific Computing.
!     Cambridge University Press, Cambridge.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, November 19, 1999.
!     Last Modified by Mark White, Battelle, November 19, 1999.
!     roots.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Parameter Statements----------------------------!
!
      PARAMETER( EPS=1.D-6, EPSL=1.D-14 )
      PARAMETER( JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1 )
!
!----------------------Type Declarations-------------------------------!
!
      EXTERNAL FUNC
      REAL*8 S(JMAXP),H(JMAXP)
!
!----------------------Executable Lines--------------------------------!
!
      IERR = 0
      ZERO = 0.D+0
!
!---  S and H store the successive trapezodial approximations and their
!     relative step-sizes.  ---
!
      H(1) = 1.D+0
      DO 10 J = 1,JMAX
        CALL TRAPZD( FUNC,A,B,S(J),J,INDX )
        IF( J.GE.K ) THEN
          CALL POLINT( H(J-KM),S(J-KM),K,ZERO,SSX,DSS,IERR )
          IF( (ABS(DSS)-EPS*ABS(SSX))/EPSL.LT.EPSL ) RETURN
        ENDIF
        S(J+1) = S(J)
        H(J+1) = 2.5D-1*H(J)
   10 CONTINUE
      IERR = 1
!
!---  End of QROMB group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ROOTS( UPTAKE,ZSX,ZEX,IP,N )
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
!     Stressed root-water uptake.
!
!----------------------Authors-----------------------------------------!
!
!     Written by EJ Freeman and MD White, PNNL, 20 June 2002.
!     Last Modified by MD White, PNNL, 20 June 2002.
!     Last Modified by MD White, PNNL, 19 August 2002.
!     Last Modified by MD White, PNNL, 20 August 2002.
!     Last Modified by MD White, PNNL, 1 December 2003.
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PLT_ATM
      USE GRID
      USE FDVP
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      EXTERNAL VRUGT
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ROOTS'
!
!---  Set plant variables  ---
!
      ZM = PARMS_P(1,IP)
!      ZSX = PARMS_P(4,IP)
      PZ = PARMS_P(7,IP)
      H1 = PARMS_P(12,IP)
      H2 = PARMS_P(13,IP)
      H3 = PARMS_P(14,IP)
      H4 = PARMS_P(15,IP)
!
!---  Vrugt root-water uptake model  ---
!
      I = ID(N)
      J = JD(N)
      K = KD(N)
      NK = ND(I,J,KFLD)
!      ZSX = MIN( ZM,ZE(5,NK)-ZE(5,N) )
!      ZEX = MIN( ZM,ZE(5,NK)-ZE(5,N) )
      CALL QROMB( VRUGT,ZSX,ZEX,RWUX,IERR,IP )
      IF( IERR.EQ.1 ) THEN
        INDX = 12
        IMSG = ND(I,J,K)
        CHMSG = 'Unconverged Romberg Integration: ' //
     &    'Source Root-Water Uptake Integration: Node: '
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Vrugt water stress function  ---
!
      IF( IRSM_P(IP).EQ.1 ) THEN
        HDGL = (PG(1,N)-PL(1,N))/(RHORL*GRAV)
        H1 = PARMS_P(12,IP)
        H2 = PARMS_P(13,IP)
        H3 = PARMS_P(14,IP)
        H4 = PARMS_P(15,IP)
        IF( HDGL.GE.H4 .OR. HDGL.LE.H1 )THEN
          WSFX = 0.D+0
        ELSEIF( HDGL.LE.H3 .AND. HDGL.GE.H2 )THEN
          WSFX = 1.D+0
        ELSEIF( HDGL.GT.H1 .AND. HDGL.LT.H2 )THEN
          WSFX = (1.D+0/(H2-H1))*(HDGL-H1)
        ELSE
          WSFX = (-1.D+0/(H4-H3))*(HDGL-H3) + 1.D+0
        ENDIF
!
!---  Jarvis water stress function  ---
!
      ELSEIF( IRSM_P(IP).EQ.2 ) THEN
        WC1 = PARMS_P(12,IP)
        WC2 = PARMS_P(13,IP)
        WC3 = PARMS_P(14,IP)
        WC4 = PARMS_P(15,IP)
        WCX = SL(1,N)*PORD(1,N)
        WSFX = (WCX-WC1)/(WC4-WC1+SMALL)
        IF( WSFX.GT.WC2 .AND. (WSFX-EPSL).LE.1.D+0 ) THEN
          WSFX = (1.D+0-WSFX)/(1.D+0-WC2+SMALL)
        ELSEIF( WSFX.GE.WC1 .AND. WSFX.LE.WC2 ) THEN
          WSFX = 1.D+0
        ELSEIF( (WSFX+EPSL).GE.0.D+0 .AND. WSFX.LT.WC1 ) THEN
          WSFX = WSFX/(WC1+SMALL)
        ELSE
          INDX = 12
          IMSG = ND(I,J,K)
          CHMSG = 'Out of Range Stress Index: Node: '
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
      UPTAKE = MAX( WSFX*RWUX,0.D+0 )
 1000 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ROOTS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SNGL_T( BCX,DB,GB,M,N,NB )
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
!     Single-temperature plant boundary condition.
!
!     Atmospheric variables
!
!     TA  Atmospheric temperature, C
!     RH  Atmospheric relative humidity, %
!     PA  Atmospheric pressure, Pa
!     RN  Solar radiation, W/m^2
!     UZ  Wind speed, m/s
!
!     Atmospheric properties
!
!     E_A - water-vapor density, kg/m^3
!     PSW_A - saturated water-vapor pressure, Pa
!     PGW_A - water-vapor partial pressure, Pa
!     HGW_A - water-vapor enthalpy, J/kg
!     UGW_A - water-vapor internal energy, J/kg
!
!     Radiation variables
!
!     RN_SNS - net short-wave radiation into the ground surface, W/m^2
!     RN_SNP - net short-wave radiation into the plant, W/m^2
!     RN_LD - downward long-wave radiation, W/m^2
!     RN_LDS - downward long-wave radiation to the ground surface, W/m^2
!     RN_LUS - upward long-wave radiation from the ground surface, W/m^2
!     RN_LNS - net long-wave radiation into the ground surface, W/m^2
!     RN_LNP - net long-wave radiation into the plant, W/m^2
!
!     Rainfall interception variables
!
!     RFIM_P - rainfall interception or condensation water mass,
!              kg/m^2 ground
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 August 2002.
!     Last Modified by MD White, PNNL, 12 August 2002.
!     Last Modified by MD White, PNNL, 12 December 2003.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOURC
      USE SOLTN
      USE PORMED
      USE PLT_ATM
      USE HYST
      USE GRID
      USE FILES
      USE FDVT
      USE FDVP
      USE FDVG
      USE CONST
      USE BCVT
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
      REAL*8 L_NS,L_SA,K_VK,KVIS_A,NU_MCF,L_SPX,L_PAX
      REAL*8 L_SAPMX,L_SAPX
      REAL*8 RKL_S(3),RFIM_PX(2),HLW_PX(2)
      REAL*8 BCX(LBCV)
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: RN_SNP,RN_LNP
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: RA_PC,RS_PC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: RW_S,RE_S
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: RW_C
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: RE_C,RE_P
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: C_COEF,RWU_S,RWU_I
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: BJ
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DNRX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: L_PA,L_SP
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: PAI,LAI
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: E_PC,PE_PC,T_PC,PT_PC
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: RFIC_P
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: RWF,AJ
      INTEGER, DIMENSION(:), ALLOCATABLE ::  IJ,KMN,IPLANT
      EXTERNAL VRUGT
!
!----------------------Data Statements---------------------------------!
!
      SAVE K_VK,ED_DC
      SAVE TMZ_O,D_O
      SAVE F_CC
      SAVE A_SA,ALB_P,ALB_P1,ALB_P2,ALB_P3,ALB_P4,ALB_P5,BETAX,
     &  C_COEF,C_COEF1,C_COEF2,C_COEF3,C_COEF4,C_COEF5,CP_A,
     &  CPG_A,D,D_O,D_X,D1,D2,D3,D4,D5,DMMX_P,DZ_ANG,E_A,EM_CS,
     &  F_1,F_CC,GBX,HG_A,HG_AI,HGA_A,HGA_AI,HGW_A,HGW_AI,
     &  HL_A,HLW_A,I,IBD,IERR,INDX,IP,IPX,J,J_DAY,K,KMN
      SAVE KVIS_A,L_PA,L_PAX,L_SA,LAI,MYR,NU_MCF,NX,NSTEP_O,
     &  PA,PAI_P,PAI,PAIX,PGA_A,PGW_A,PH_AVG,PR_A,PSW_A,
     &  RA_BSMX,RA_BSX,RA_CA,RA_PC,RA_SC,RA_X,RB_P,RE_MCF,
     &  RFIC_P,RH,RHO_A,RHOG_A,RHOG_AI,RHOGA_A,RHOGA_AI,
     &  RHOGW_A,RHOGW_AI,RHOL_A,RL_P,RN,RN_LD,RN_LDS,RN_SDS
      SAVE RN_SNP,RS_MIN,RS_PC,RWF,RWU_I,
     &  RWU_S,SIGMA,T_DP,TA,TA_K,TA_Y,TAI,THD_A,THK_A,TKAX,
     &  TKCX,TKHX,TKOX,TMZ,TMZ_O,UGA_A,UGA_AI,UGW_A,UGW_AI,
     &  UPTAKE,UZ,UZ_F,UZ_FX,UZ_MCF,VIS_A,VLG_1,VLG_2,
     &  WL_P,XGA_A,XGA_AI,XGW_A,XGW_AI,XLA_A,XLW_A,XMLA_A, 
     &  XMLW_A,ZEX,ZM,ZO,ZP_D
      SAVE RN_LNP,RW_S,RE_S,RW_C,RE_C,RE_P,BJ,
     &  DNRX,L_SP,E_PC,PE_PC,T_PC,PT_PC,AJ,IJ,MPX,IPLANT
      DATA K_VK,ED_DC / 4.1D-1,2.5D+0 /
      DATA TMZ_O,D_O / -1.D+20,0.D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SNGL_T'
!
!---  Dynamic memory allocation  ---
!
      KP = 2+LPLANT*3
      IF( NSTEP-NRST.EQ.0 .AND. M.EQ.2 ) THEN
        ALLOCATE( RN_SNP(1:LPLANT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: RN_SNP'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( RN_LNP(1:LPLANT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: RN_LNP'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( RA_PC(1:LPLANT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: RA_PC'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( RS_PC(1:LPLANT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: RS_PC'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( RW_S(1:3+LPLANT*3),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: RW_S'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( RE_S(1:3+LPLANT*3),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: RE_S'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( RW_C(1:3+LPLANT*3),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: RW_C'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( RE_C(1:3+LPLANT*3),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: RE_C'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( RE_P(1:3+LPLANT*3),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: RE_P'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( C_COEF(1:LPLANT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: C_COEF'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( RWU_S(1:LPLANT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: RWU_S'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( RWU_I(1:LPLANT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: RWU_I'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( RWF(1:LPLANT,1:LFZ),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: RWF'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( AJ(1:KP,1:KP),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: AJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( BJ(1:KP),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: BJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( DNRX(1:KP),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: DNRX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( L_PA(1:LPLANT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: L_PA'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( L_SP(1:LPLANT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: L_SP'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( PAI(1:LPLANT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: PAI'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( LAI(1:LPLANT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: LAI'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( E_PC(1:LPLANT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: E_PC'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( PE_PC(1:LPLANT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: PE_PC'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( T_PC(1:LPLANT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: T_PC'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( PT_PC(1:LPLANT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: PT_PC'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( RFIC_P(1:LPLANT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: RFIC_P'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( IJ(1:KP),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: IJ'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( KMN(1:LPLANT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: KMN'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( IPLANT(1:LPLANT),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: IPLANT'
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
!
!---  Initialize variables  ---
!
      GBX = GB/DB
      IBD = ABS(IBCD(NB))
      K = KD(N)
      IF( IBCD(NB).EQ.-3 ) THEN
        A_SA = AFZ(NSZ(N))
      ELSEIF( IBCD(NB).EQ.-2 ) THEN
        A_SA = AFY(NSY(N))
      ELSEIF( IBCD(NB).EQ.-1 ) THEN
        A_SA = AFX(NSX(N))
      ELSEIF( IBCD(NB).EQ.1 ) THEN
        A_SA = AFX(NSX(N)+1)
      ELSEIF( IBCD(NB).EQ.2 ) THEN
        A_SA = AFY(NSY(N)+IFLD)
      ELSEIF( IBCD(NB).EQ.3 ) THEN
        A_SA = AFZ(NSZ(N)+IJFLD)
      ENDIF
!
!---  Skip atmospheric calculations  ---
!
      MPOX = 1
      NC = 0
      IF( M.EQ.2 ) MPX = 6
      IF( NSTEP-NRST.EQ.0 ) NSTEP_O = -1
      IF( M.GT.2 .OR. NSTEP.EQ.NSTEP_O ) GOTO 300
      NSTEP_O = NSTEP
!
!---  Time of year  ---
!
      TMZ = TM
      IF( NSTEP-NRST.EQ.0 ) TMZ = TMZ*(1.D+0+EPSL)+EPSL
      IF( TMZ.EQ.TMZ_O ) THEN
        D = D_O
      ELSE
        D = TMZ + ATMST
        MYR = 1900
   90   CONTINUE
        MYR = MYR+1
        D_X = D
        D = D - 3.1536D+7
        IF( MOD(MYR,1000).EQ.0 ) THEN
          IF( MOD(MYR,400).EQ.0 ) THEN
            D = D - 8.64D+4
          ENDIF
        ELSEIF( MOD(MYR,4).EQ.0 ) THEN
          D = D - 8.64D+4
        ENDIF
        IF( D.GT.0.D+0 ) GOTO 90
        D = D_X
        D_O = D
        TMZ_O = TMZ
        J_DAY = INT(D/8.64D+4) + 1
      ENDIF
!
!---  Set atmospheric conditions ---
!
      CALL ATMOS_C( PA,RH,RN,TA,UZ,IERR )
!
!---  Simulation time greater than ending atmospheric conditions time
!     skip Shuttleworth-Wallace calculations  ---
!
      IF( IERR.EQ.1 )  GOTO 1000
  100 CONTINUE
!
!---  Net short-wave radiation at ground surface, considering
!     interception by plants, W/m^2 ground  ---
!
!     RN - downward, short-wave radiation, W/m^2 ground
!     RN_SDS - downward short-wave radiation
!       at ground surface, W/m^2 ground
!     RN_SUS - upward short-wave radiation
!       at ground surface, W/m^2 ground
!     RN_SNS - net short-wave radiation at ground surface, W/m^2 ground
!     PAI_P - summation of the plant-area index of all plant species
!     LAI(IP) - leaf-area index of plant specie IP,
!       m^2 leaf/m^2 plant ground
!     PAI(IP) - plant-area index of plant specie IP,
!       m^2 plant ground/m^2 ground
!     ALB_S - short-wave albedo of ground surface
!
      PAI_P = 0.D+0
      RN_SDS = 0.D+0
      DO 110 IP = 1,NPLANT
        IPLANT(IP) = 1
        IPX = (IP-1)*2 + 4
        PAI(IP) = BCX(IPX)
        LAI(IP) = BCX(IPX-1)
        IF( PAI(IP)*LAI(IP)/EPSL.LT.EPSL ) THEN
          IPLANT(IP) = 0
          GOTO 110
        ENDIF
        PAI_P = PAI_P + PAI(IP)
        RN_SDS = RN_SDS + PAI(IP)*EXP(-7.D-1*LAI(IP))
  110 CONTINUE
      RN_SDS = RN*((1.D+0-PAI_P) + RN_SDS)
      IF( PAI_P.GT.1.D+0 ) THEN
        INDX = 14
        N_DB = N
        CHMSG = 'Combined Plant Area Index Greater Than One: '
        RLMSG = PAI_P
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Downward long-wave radiation,
!     W/m^2 ground surface  ---
!
!     RN_LD - downward long-wave radiation, W/m^2
!     EM_CS - clear-sky emissivity
!     F_CC - fractional cloud cover
!     SIGMA - Stefan-Boltzmann constant, W/m^2 K^4
!     TA_K - atmospheric temperature, K
!
      CALL DEW_PT( RH,TA,T_DP )
      CALL CLOUD( RN,F_CC,J_DAY,DZ_ANG )
      EM_CS = MIN( 0.741D+0 + 0.0062*T_DP,1.D+0 )
      SIGMA = 5.67D-8
      TA_K = TA + TABS
      RN_LD = EM_CS*(1.D+0-0.84D+0*F_CC)*SIGMA*(TA_K**4)
     &  + 0.84D+0*F_CC*SIGMA*(TA_K**4)
!
!---  Downward long-wave radiation to the ground surface,
!     excluding emission from plants, W/m^2 ground surface  ---
!
!     RN_LD - downward long-wave radiation, W/m^2
!     RN_LDS - downward long-wave radiation to ground surface, W/m^2
!     PAI(IP) - plant-area index of plant specie IP
!     LAI(IP) - leaf-area index of plant specie IP
!
      RN_LDS = 0.D+0
      DO 120 IP = 1,NPLANT
        IF( IPLANT(IP).EQ.0 ) GOTO 120
        RN_LDS = RN_LDS + PAI(IP)*EXP(-7.D-1*LAI(IP))
  120 CONTINUE
      RN_LDS = RN_LD*((1.D+0-PAI_P) + RN_LDS)
!
!---  Atmospheric properties  ---
!
!     E_A - water-vapor density, kg/m^3
!     PSW_A - saturated water-vapor pressure, Pa
!     PGW_A - atmospheric water-vapor pressure, Pa
!     PGA_A - atmospheric air partial pressure, Pa
!     CPG_A - atmospheric specific heat, J/kg
!     RHOGW_A - atmospheric water-vapor density, kg/m^3
!     RHOGA_A - atmospheric air density, kg/m^3
!     RHOG_A - atmospheric gas density, kg/m^3
!     HGW_A - atmospheric water-vapor enthalpy, J/kg
!     UGW_A - atmospheric water-vapor internal energy, J/kg
!     HGA_A - atmospheric air enthalpy, J/kg
!     UGA_A - atmospheric air internal energy, J/kg
!     RHO_A - dry-air density, kg/m^3
!     VIS_A - dry-air viscosity, Pa s
!     KVIS_A - dry-air kinematic viscosity, m^2/s
!     CP_A - dry-air specific heat, J/kg
!     THK_A - dry-air thermal conductivity, W/m K
!     THD_A - dry-air thermal diffusivity, m^2/s
!     PR_A - dry-air Prandtl number
!
      INDX = 0
      CALL REGION_4( TA,PSW_A,INDX )
      PGW_A = RH*PSW_A
      INDX = 1
      CALL WATGSD( TA,PGW_A,E_A,INDX )
      CALL WATGSH( TA,PSW_A,E_A,HGW_A,UGW_A )
      PGA_A = MAX( PA-PGW_A,0.D+0 )
      INDX = 1
      CALL WATGSD( TA,PGW_A,RHOGW_A,INDX )
      CALL AIRGSD( TA,PGA_A,RHOGA_A )
      RHOG_A = RHOGW_A+RHOGA_A
      XGW_A = RHOGW_A/RHOG_A
      XGA_A = RHOGA_A/RHOG_A
      CALL WATGSH( TA,PSW_A,RHOGW_A,HGW_A,UGW_A )
      CALL AIRGSH( TA,PGA_A,HGA_A,UGA_A )
      HG_A = XGW_A*HGW_A + XGA_A*HGA_A
      TAI = TA + 1.D-3
      INDX = 1
      CALL WATGSD( TAI,PGW_A,RHOGW_AI,INDX )
      CALL AIRGSD( TAI,PGA_A,RHOGA_AI )
      RHOG_AI = RHOGW_AI+RHOGA_AI
      XGW_AI = RHOGW_AI/RHOG_AI
      XGA_AI = RHOGA_AI/RHOG_AI
      CALL WATGSH( TAI,PSW_A,RHOGW_AI,HGW_AI,UGW_AI )
      CALL AIRGSH( TAI,PGA_A,HGA_AI,UGA_AI )
      HG_AI = XGW_AI*HGW_AI + XGA_AI*HGA_AI
      CPG_A = (HG_AI-HG_A)/1.D-3
      CALL AIRGSD( TA,PA,RHO_A )
      CALL AIRGSV( TA,VIS_A )
      KVIS_A = VIS_A/RHO_A
      CALL AIRGSC( TA,CP_A )
      CALL AIRGSK( TA,THK_A )
      THD_A = THK_A/(RHO_A*CP_A)
      PR_A = KVIS_A/THD_A
!
!---  Precipitation properties  ---
!
      XMLA_A = PGA_A/HCAW
      XMLW_A = MAX( 1.D+0-XMLA_A,0.D+0 )
      XLW_A = XMLW_A*WTMW/(XMLA_A*WTMA + XMLW_A*WTMW)
      XLA_A = MAX( 1.D+0-XLW_A,0.D+0 )
      TA_Y = MAX( TA,1.D-2 )
      CALL WATLQD( TA_Y,PA,RHOL_A )
      CALL WATLQH( TA_Y,PGA_A,HLW_A )
      HL_A = XLW_A*HLW_A + XLA_A*HGA_A
!
!---  Eddy diffusion resistance above and below the average
!     plant height  ---
!
!     PH_AVG - average plant height, m
!     PAI(IP) - plant-area index of plant specie IP
!     PARMS_P(11,IP) - height of plant specie IP
!     ZP_D - zero-plane displacement, m
!     RL_P - plant roughness length, m
!     RA_SC - aerodynamic resistance between ground and
!             plant height, s/m
!     RA_CA - aerodynamic resistance between mean canopy flow and
!             reference height, s/m
!     ATMC(1) - wind speed reference height, m
!     K_VK - von Karman constant, 0.4
!     ED_DC - eddy diffusion decay constant
!     UZ  Wind speed, m/s
!     ATMC(6) - momentum transport surface roughness length, m
!     ATMC(7) - heat and mass transport surface roughness length, m
!     RA_BSMX - maximum aerodynamic resistance between ground and
!       reference height for bare surface, s/m
!
      PH_AVG = 0.D+0
      PAIX = 0.D+0
      DO 130 IP = 1,NPLANT
        IF( IPLANT(IP).EQ.0 ) GOTO 130
        PH_AVG = PH_AVG + PAI(IP)*PARMS_P(11,IP)
        PAIX = PAIX + PAI(IP)
  130 CONTINUE
      PH_AVG = PH_AVG/PAIX
      IF( PH_AVG.GT.EPSL ) THEN
        ZP_D = 0.63D+0*PH_AVG
        RL_P = 0.13D+0*PH_AVG
        RA_X = (LOG((ATMC(1)-ZP_D)/RL_P)/((K_VK**2)*UZ))*
     &    (PH_AVG/(ED_DC*(PH_AVG-ZP_D)))
        RA_SC = RA_X*
     &    (EXP(ED_DC) - EXP(ED_DC*(1.D+0-((ZP_D+RL_P)/PH_AVG))))
        RA_CA = RA_X*
     &    ((EXP(ED_DC*(1.D+0-((ZP_D+RL_P)/PH_AVG)))-1.D+0) +
     &    LOG((ATMC(1)-ZP_D)/(PH_AVG-ZP_D)))
        RA_BSMX = RHOG_A*CPG_A*ATMC(2)/THK_A
        VLG_1 = LOG((ATMC(1)+ATMC(6))/ATMC(6))
        VLG_2 = LOG((ATMC(2)+ATMC(7))/ATMC(7))
        UZ_FX = K_VK*UZ/VLG_1
        RA_BSX = MIN( VLG_2/(K_VK*UZ_FX),RA_BSMX )
      ELSE
!
!---  Eddy diffusion resistance for bare surface  ---
!
!     ATMC(1) - wind speed reference height, m
!     K_VK - von Karman constant, 0.4
!     UZ  Wind speed, m/s
!     UZ_FX  Friction velocity, m/s
!     RA_BSX - aerodynamic resistance between ground and
!       reference height for bare surface, s/m
!     RA_BSMX - maximum aerodynamic resistance between ground and
!       reference height for bare surface, s/m
!     ATMC(6) - momentum transport surface roughness length, m
!     ATMC(7) - heat and mass transport surface roughness length, m
!
        RA_BSMX = RHOG_A*CPG_A*ATMC(2)/THK_A
        VLG_1 = LOG((ATMC(1)+ATMC(6))/ATMC(6))
        VLG_2 = LOG((ATMC(2)+ATMC(7))/ATMC(7))
        UZ_FX = K_VK*UZ/VLG_1
        RA_BSX = MIN( VLG_2/(K_VK*UZ_FX),RA_BSMX )
        RA_SC = 5.D-1*RA_BSX
        RA_CA = 5.D-1*RA_BSX
      ENDIF
!
!---  Aerodynamic resistance between the plant leaves and
!     mean canopy flow  ---
!
!     UZ_MCF - wind speed at the mean canopy flow height, m/s
!     RE_MCF - Reynolds number at the mean canopy flow height
!     PR_A - dry-air Prandtl number
!     NU_MCF - Nusselt number at the mean canopy flow height
!     RHO_A - dry-air density, kg/m^3
!     VIS_A - dry-air viscosity, Pa s
!     KVIS_A - dry-air kinematic viscosity, m^2/s
!     CPG_A - dry-air specific heat, J/kg
!     THK_A - dry-air thermal conductivity, W/m K
!     THD_A - dry-air thermal diffusivity, m^2/s
!     RB_P - leaf boundary layer resistance, s/m
!     RA_PC - aerodynamic resistance between plant leaves and
!             mean canopy flow, s/m
!     WL_P - leaf width, m
!
      IF( PH_AVG.GT.EPSL ) THEN
        UZ_F = UZ*K_VK/LOG((ATMC(1)-ZP_D+RL_P)/RL_P)
        UZ_MCF = (UZ_F/K_VK)*LOG((PH_AVG-ZP_D+RL_P)/RL_P)
!
!---  Use dry air properties to compute Reynolds, Prandtl,
!     and Nusselt numbers ---
!
        DO 140 IP = 1,NPLANT
          IF( IPLANT(IP).EQ.0 ) GOTO 140
          WL_P = 1.D-2*LAI(IP)
          RE_MCF = WL_P*UZ_MCF/KVIS_A
          NU_MCF = (2.D+0/3.D+0)*(RE_MCF**5.D-1)*(PR_A**(THIRD))
          RB_P = 7.D-1*WL_P/(THD_A*NU_MCF)
          RA_PC(IP) = RB_P/LAI(IP)
  140   CONTINUE
      ENDIF
!
!---  Stomatal resistance  ---
!
!     RS_PC - stomatal resistance between plant leaf interior
!             and exterior, s/m
!     RS_MIN - minimum stomatal resistance between plant leaf interior
!             and exterior, s/m
!     F_1 - influence of solar radiation
!
!     Hicks stomatal resistance model
!
!     PARMS_P(2,IP) - minimum stomatal resistance, s/m
!     PARMS_P(3,IP) - light response coefficient, W/m^2
!     PARMS_P(27,IP) - minimum temperature for stomatal opening, K
!     PARMS_P(28,IP) - maximum temperature for stomatal opening, K
!     PARMS_P(29,IP) - optimum temperature for stomatal opening, K
!
!---  Hicks stomatal resistance model  ---
!
      DO 150 IP = 1,NPLANT
        IF( ISRM_P(IP).EQ.1 ) THEN
          RS_MIN = PARMS_P(2,IP)
          BETAX = PARMS_P(3,IP)
          TKAX = TA+TABS
          TKCX = PARMS_P(27,IP)
          TKHX = PARMS_P(28,IP)
          TKOX = PARMS_P(29,IP)
          F_1 = ((TKAX-TKCX)/(TKOX-TKCX))*
     &      (((TKHX-TKAX)/(TKHX-TKOX))**((TKHX-TKOX)/(TKOX-TKCX)))
          F_1 = MAX( F_1,EPSL )
          IF( IPLANT(IP).EQ.0 ) THEN
            RS_PC(IP) = RS_MIN
          ELSE
            RS_PC(IP) = RS_MIN*(1.D+0 + BETAX/(RN+SMALL))/F_1/LAI(IP)
          ENDIF
        ELSE
          RS_MIN = 5.D+1
          F_1 = (4.D+2+RN)/(1.4D+0*RN+EPSL)
          IF( IPLANT(IP).EQ.0 ) THEN
            RS_PC(IP) = RS_MIN
          ELSE
            RS_PC(IP) = RS_MIN*F_1/LAI(IP)
          ENDIF
        ENDIF
  150 CONTINUE
!
!---  Crop growth index  ---
!
      DO 160 IP = 1,NPLANT
        IF( IPLANT(IP).EQ.0 ) THEN
          C_COEF(IP) = 0.D+0
        ELSE
          C_COEF1 = PARMS_P(17,IP)
          D1 = PARMS_P(18,IP)
          C_COEF2 = PARMS_P(19,IP)
          D2 = PARMS_P(20,IP)
          C_COEF3 = PARMS_P(21,IP)
          D3 = PARMS_P(22,IP)
          C_COEF4 = PARMS_P(23,IP)
          D4 = PARMS_P(24,IP)
          C_COEF5 = PARMS_P(25,IP)
          D5 = PARMS_P(26,IP)
          IF( D.GE.D1 .AND. D.LT.D2 )THEN
            C_COEF(IP) = ((C_COEF2-C_COEF1)/(D2-D1))*(D-D1)+C_COEF1
          ELSEIF( D.GE.D2 .AND. D.LT.D3 )THEN
            C_COEF(IP) = ((C_COEF3-C_COEF2)/(D3-D2))*(D-D2)+C_COEF2
          ELSEIF( D.GE.D3 .AND. D.LT.D4 )THEN
            C_COEF(IP) = ((C_COEF4-C_COEF3)/(D4-D3))*(D-D3)+C_COEF3
          ELSEIF( D.GE.D4 .AND. D.LT.D5 )THEN
            C_COEF(IP) = ((C_COEF5-C_COEF4)/(D5-D4))*(D-D4)+C_COEF4
          ELSE
            C_COEF(IP) = 0.D+0
          ENDIF
        ENDIF
  160 CONTINUE
!
!---  Rainfall via boundary condition input  ---
!
!     L_SA - rainfall intensity on ground surface, m^3/s m^2 ground
!     L_PA(IP) - rainfall intensity for plant IP, m^3/s m^2 ground
!     L_PAX - rainfall intensity on plants, m^3/s m^2 ground
!
      IF( BCX(2)/EPSL.GT.EPSL ) THEN
        L_SA = -BCX(2)
      ELSE
        L_SA = 0.D+0
      ENDIF
      L_PAX = 0.D+0
!
!---  Rainfall interception by plants  ---
!
      DO 170 IP = 1,NPLANT
        RFIC_P(IP) = 0.D+0
        L_PA(IP) = 0.D+0
        IF( IPLANT(IP).EQ.0 ) GOTO 170
!
!---    Rainfall interception capacity, kg/m^2 ground  ---
!
!       PARMS_P(16,IP) - maximum dew depth for plant specie IP, m
!       DMMX_P -  maximum dew mass, kg/m^2 ground
!       RFIC_P(IP) - rainfall interception capacity, kg/m^2 ground
!       RFI_CF - rainfall interception capacity coefficient, s/m
!
        IF( ISLC(26).EQ.1 ) THEN
          DMMX_P = PARMS_P(16,IP)*LAI(IP)*PAI(IP)*RHOL_A
          RFIC_P(IP) = DMMX_P
        ENDIF
!
!---    Reduce rainfall intensity by plant area index,
!       m^3/s m^2 ground
!
        L_PA(IP) = L_SA*PAI(IP)
!
!---    Sum over plants, incident rainfall m^3/m^3 ground  ---
!
        L_PAX = L_PAX + L_PA(IP)
  170 CONTINUE
      L_SA = L_SA - L_PAX
!
!---  Skip point for atmospheric conditions  ---
!
  300 CONTINUE
!
!---  Skip root-water uptake calculations  ---
!
      IF( M.GT.2 ) GOTO 340
!
!---  Potential root-water uptake based on vertical 
!     root distribution  ---
!
      DO 330 IP = 1,NPLANT
        IF( IPLANT(IP).EQ.0 ) GOTO 330
        ZO = 0.D+0
        ZM = PARMS_P(1,IP)
        CALL QROMB( VRUGT,ZO,ZM,RWU_I(IP),IERR,IP )
        IF( IERR.EQ.1 ) THEN
          INDX = 3
          CHMSG = 'Unconverged Romberg Integration: ' //
     &      'Root-Water Uptake Integration'
          CALL WRMSGS( INDX )
        ENDIF
!
!---  Actual root-water uptake based on vertical 
!     root distribution and root stress model, using
!     previous time-step capillary pressure and saturation  ---
!
        RWU_S(IP) = 0.D+0
        I = ID(N)
        J = JD(N)
        KMN(IP) = KFLD+1
        KTOP = 0
        DO 310 K = KFLD,1,-1
          NX = ND(I,J,K)
          IF( IXP(NX).EQ.0 ) GOTO 310
          IF( KTOP.EQ.0 ) KTOP = K+1
!          NK = ND(I,J,KTOP)
!          ZSX = MIN( ZM,ZE(1,NK)-ZE(5,N) )
!          ZEX = MIN( ZM,ZE(1,NK)-ZE(1,N) )
          ZSX = MIN( ZM,ZE(5,N)-ZE(5,NX) )
          ZEX = MIN( ZM,ZE(5,N)-ZE(1,NX) )
          CALL ROOTS( UPTAKE,ZSX,ZEX,IP,NX )
          RWU_S(IP) = RWU_S(IP) + UPTAKE
          RWF(IP,K) = UPTAKE
!          RWU_S(IP) = MIN( RWU_I(IP),RWU_S(IP) )
!          NK = ND(I,J,KFLD)
!          ZEX = MIN( ZM,ZE(5,NK)-ZE(1,N) )
          KMN(IP) = K
          IF( (ZM-ZEX)/EPSL.LT.EPSL ) GOTO 312
  310   CONTINUE
  312   CONTINUE
!
!---    Root water fraction, (i.e., nodal fraction
!       of total actual root-water uptake)  ---
!
        KTOP = 0
        DO 320 K = KFLD,1,-1
          NX = ND(I,J,K)
          IF( IXP(NX).EQ.0 ) GOTO 320
          IF( KTOP.EQ.0 ) KTOP = K+1
!          NK = ND(I,J,KTOP)
!          ZSX = MIN( ZM,ZE(1,NK)-ZE(5,N) )
!          ZEX = MIN( ZM,ZE(1,NK)-ZE(1,N) )
          ZSX = MIN( ZM,ZE(5,N)-ZE(5,NX) )
          ZEX = MIN( ZM,ZE(5,N)-ZE(1,NX) )
!          CALL ROOTS( UPTAKE,ZSX,ZEX,IP,NX )
          IF( RWU_S(IP)/EPSL.LT.EPSL ) THEN
            RWF(IP,K) = 0.D+0
          ELSE
            RWF(IP,K) = RWF(IP,K)/RWU_S(IP)
          ENDIF
!          NK = ND(I,J,KFLD)
!          ZEX = MIN( ZM,ZE(5,NK)-ZE(1,N) )
          IF( (ZM-ZEX)/EPSL.LT.EPSL ) GOTO 330
  320   CONTINUE
  330 CONTINUE
!
!---  Skip point for root-water uptake calculations  ---
!
  340 CONTINUE
!
!---  Initialize convergence variables  ---
!
      RHOL_SO = RHOLB(1,NB)
      RHOG_SO = RHOGB(1,NB)
      SL_SO = SLB(1,NB)
      HL_SO = HLB(1,NB)
      UG_SO = UEGB(1,NB)
      T_SOX = TB(1,NB)
      XGW_SO = XGWB(1,NB)
      XLW_SO = XLWB(1,NB)
      DP_SO = 0.D+0
      DP_SOO = 0.D+0
!
!---  Skip iterative calculations  ---
!
      IF( (M.NE.2 .OR. NITER.GT.1) .AND. ISLC(24).EQ.11 ) MPX = 1
      IF( NITER.GT.(NRIMX/2) ) MPX = 1
!
!---  Guess unknowns  ---
!
!     T_S - surface temperature, C
!     PL_S - surface aqueous pressure, Pa
!     T_P - plant temperature, C
!     T_C - air temperature at the mean canopy flow height, C
!     PGW_C - water-vapor partial pressure at the mean canopy
!             flow height, Pa
!
      PG_S = PA
      IF( M.EQ.2 .AND. (NSTEP-NRST).EQ.0 ) THEN
        T_S = TA
        RHX = SQRT(SQRT(MAX(RH,0.D+0)))
        CALL WATLQD( T_S,PG_S,RHOL_S )
        IF( RH.LT.9.D-1 ) THEN
          PL_S = LOG(RHX)*(RHOL_S*RCW*(T_S+TABS)) + PG_S
        ELSEIF( BCX(2)/EPSL.GT.EPSL ) THEN
          PL_S = 0.99D+0*PG_S
        ELSE
          PL_S = RH*PG_S
        ENDIF
        PL_S = MIN( PL_S,PG_S+SWBCD*RHOL_A*GRAV )
        T_P = RH*(T_DP-TA) + TA
        T_C = RH*(T_DP-TA) + TA
        PGW_C = PGW_A
      ELSE
        T_S = T_SO(2,NB)
!
!---    Guess saturation under rainfall conditions  ---
!
        IF( M.EQ.2 .AND. NITER.EQ.1 .AND. L_SA.LT.-EPSL ) THEN
          SLX = MIN( 1.D+0,ABS(L_SA)*VISRL/(PERM(3,IZ(N))*RHORL*GRAV) )
          CALL CAP_GT( IZ(N),SLX,CPGLX )
          PL_S = MAX( PG_S-CPGLX,PL_SO(2,NB) )
        ELSE
          PL_S = PL_SO(2,NB)
        ENDIF
        PL_S = MIN( PL_S,PG_S+SWBCD*RHOL_A*GRAV )
        T_P = T_PO(2,1,NB)
        T_C = T_CO(2,1,NB)
        PGW_C = PGW_CO(2,1,NB)
      ENDIF
      T_SMN = 0.D+0
      PL_SMN = 0.D+0
      T_PMN = 0.D+0
      T_CMN = 0.D+0
      PGW_CMN = 0.D+0
      SL_SMN = 0.D+0
  350 CONTINUE
      NC = NC+1
!
!---  Primary variable increments  ---
!
      DNRX(1) = 1.D-4
      DNRX(2) = 1.D-4
      DNRX(3) = 1.D-1
      IF( MPX.GT.0 ) THEN
        DNRX(4) = -MAX( 1.D-1,1.D-6*MAX(PG_S-PL_S,0.D+0) )
        DNRX(5) = 1.D-4
      ELSE
        DNRX(4) = 1.D-4
        DNRX(5) = -MAX( 1.D-1,1.D-6*MAX(PG_S-PL_S,0.D+0) )
      ENDIF
!
!---  Loop over number of primary unknowns  ---
!
!     T_S - surface temperature, C
!     PL_S - surface water-vapor partial pressure, Pa (abs)
!     T_P - plant temperature, C
!     T_C - plant canopy height temperature, C
!     PGW_C - plant canopy water-vapor partial pressure, Pa (abs)
!     MPX (-5) fix surface aqueous pressure 
!     MPX (+5) fix surface temperature
!     MPX (+6) solve surface aqueous pressure and temperature
!
      DO 430 MP = 1,ABS(MPX)
!
!---    Increment primary variables  ---
!
        T_SX = T_S
        PL_SX = PL_S
        T_PX = T_P
        T_CX = T_C
        PGW_CX = PGW_C
        IF( MP.EQ.2 ) T_PX = T_P + DNRX(1)
        IF( MP.EQ.3 ) T_CX = T_C + DNRX(2)
        IF( MP.EQ.4 ) PGW_CX = PGW_C + DNRX(3)
        IF( MPX.GT.0 ) THEN
          IF( MP.EQ.5 ) PL_SX = PL_S + DNRX(4)
          IF( MP.EQ.6 ) T_SX = T_S + DNRX(5)
        ELSE
          IF( MP.EQ.5 ) T_SX = T_S + DNRX(4)
          IF( MP.EQ.6 ) PL_SX = PL_S + DNRX(5)
        ENDIF
        D_NS = ABS(DB)
!
!---    Limit property data temperature to above freezing  ---
!
        T_SY = MAX( T_SX,1.D-1 )
        T_PY = MAX( T_PX,1.D-1 )
        T_CY = MAX( T_CX,1.D-1 )
!
!---    Water vapor pressure at ground surface
!       w/ vapor pressure lowering  ---
!
        CALL WATLQD( T_SY,PG_S,RHOL_S )
        INDX = 0
        CALL REGION_4( T_SX,PSW_S,INDX )
        PC_S = MAX( PG_S-PL_SX,0.D+0 )
        IF( ISLC(44).EQ.1 ) THEN
          PGW_S = PSW_S
        ELSE
          CALL VPLWR( PC_S,PSW_S,T_SX,RHOL_S,RCW,PGW_S )
        ENDIF
!
!---    Water-vapor density at ground surface, kg/m^3  ---
!
        PGA_S = MAX( PG_S-PGW_S,0.D+0 )
        XMLA_S = PGA_S/HCAW
        XMLW_S = MAX( 1.D+0-XMLA_S,0.D+0 )
        XLW_S = XMLW_S*WTMW/(XMLA_S*WTMA + XMLW_S*WTMW)
        XLA_S = MAX( 1.D+0-XLW_S,0.D+0 )
        INDX = 1
        CALL WATGSD( T_SX,PGW_S,RHOGW_S,INDX )
        CALL AIRGSD( T_SX,PGA_S,RHOGA_S )
        RHOG_S = RHOGW_S+RHOGA_S
        XGW_S = RHOGW_S/RHOG_S
        XGA_S = RHOGA_S/RHOG_S
        E_S = RHOGW_S
!
!---    Gas viscosity at ground surface  ---
!
        CALL WATGSV( T_SX,VISW_S )
        CALL AIRGSV( T_SX,VISA_S )
        XMGW_S = PGW_S/PG_S
        XMGA_S = PGA_S/PG_S
        CALL GASVIS( XMGW_S,ZERO,XMGA_S,VISW_S,SMALL,VISA_S,VISG_S )
!
!---    Aqueous pressure at ground surface, from vapor-pressure
!       lowering function  ---
!
        CALL WATLQD( T_SY,PG_S,RHOL_S )
        INDX = 0
        CALL REGION_4( T_SX,PSW_S,INDX )
!
!---    Saturated water-vapor density at ground surface
!       temperature, kg/m^3  ---
!
        CALL WATGSD( T_SX,PSW_S,PE_S,INDX )
!
!---    Aqueous viscosity at the ground surface, Pa s  ---
!
        CALL WATLQV( T_SY,PG_S,PSW_S,VISL_S )
!
!---    Saturation and relative permeability at ground surface  ---
!
        CALL KSP_GT( N,IZ(N),M,PG_S,PL_SX,SL_S,SG_S,RKL_S,RKG_S )
        IF( MP.EQ.1 ) THEN
          SL_SX = SL_S
        ENDIF
!
!---    Net short-wave radiation intercepted by ground surface,
!       W/m^2 ground surface  ---
!
!       RN_SNS - net short-wave radiation transitted to the
!         ground surface, W/m^2 ground
!       RN_SDS - downward short-wave radiation
!         at ground surface, W/m^2 ground
!       RN_SUS - upward short-wave radiation
!         at ground surface, W/m^2 ground
!       ALB_S - short-wave albedo of ground surface (0.1 - 0.4)
!
!---    Moisture ground-surface solar angle model  ---
!
        IF( IALB(IZ(N)).GT.0 ) THEN
          ALB_S = ALBEDO(2,IZ(N)) + 
     &      (ALBEDO(1,IZ(N))-ALBEDO(2,IZ(N)))*
     &      EXP(-ALBEDO(3,IZ(N))*SL_S*PORD(2,N))
        ENDIF
!
!---    Plem and Xiu ground-surface solar angle model  ---
!
        IF( IALB(IZ(N)).EQ.1 ) THEN
          ALB_S = ALB_S + 1.D-2*(EXP(3.286D-3*(DZ_ANG**1.5D+0))-1.D+0)
!
!---    Briegleb ground-surface solar angle model  ---
!
        ELSEIF( IALB(IZ(N)).EQ.2 ) THEN
          Z_ANG = DZ_ANG*GPI/1.8D+2
          ALB_S = ALB_S + ALBEDO(4,IZ(N))*
     &      (1.D+0 + ALBEDO(5,IZ(N)))/
     &      (1.D+0 + 2.D+0*ALBEDO(5,IZ(N))*COS(Z_ANG))
!
!---    Wang ground-surface solar angle model  ---
!
        ELSEIF( IALB(IZ(N)).EQ.3 ) THEN
          Z_ANG = DZ_ANG*GPI/1.8D+2
          G_1 = (-7.574D-3) + (-7.0987D-2)*(Z_ANG**2) + 
     &      (3.07588D-1)*(Z_ANG**3)
          G_2 = (-1.284909D+0) + (-1.66314D-1)*(Z_ANG**2) +
     &      (4.184D-2)*(Z_ANG**3)
          G_1R = 2.67808D-1
          G_2R = -1.419244D+0
          ALB_S = ALB_S + ALBEDO(4,IZ(N))*(1.D+0 +
     &      0.346D+0*(G_1-G_1R) + 0.063*(G_2-G_2R))
        ELSE
          ALB_S = 0.25D+0
        ENDIF
        RN_SUS = RN_SDS*ALB_S
        RN_SNS = (RN_SDS-RN_SUS)
!
!---    Net short-wave radiation at plants, considering
!       single reflection from ground surface, W/m^2 ground  ---
!
!       RN_SNP(IP) - net short-wave radiation into plant
!         specie IP, W/m^2 ground
!       ALB_P - short-wave albedo of plant specie IP
!       LAI(IP) - leaf-area index of plant specie IP,
!         m^2 leaf/m^2 plant ground
!       PAI(IP) - plant-area index of plant specie IP,
!         m^2 plant ground/m^2 ground
!       PAI_P - plant-area index of all plant species
!
        DO 360 IP = 1,NPLANT
          RN_SNP(IP) = 0.D+0
          IF( IPLANT(IP).EQ.0 ) GOTO 360
!
!---      Temporal plant solar albedo  ---
!
          IF( IALB_P(IP).EQ.1 ) THEN
            ALB_P1 = PARMS_P(5,IP)
            D1 = PARMS_P(18,IP)
            ALB_P2 = PARMS_P(6,IP)
            D2 = PARMS_P(20,IP)
            ALB_P3 = PARMS_P(8,IP)
            D3 = PARMS_P(22,IP)
            ALB_P4 = PARMS_P(9,IP)
            D4 = PARMS_P(24,IP)
            ALB_P5 = PARMS_P(10,IP)
            D5 = PARMS_P(26,IP)
            IF( D.GE.D1 .AND. D.LT.D2 )THEN
              ALB_P = ((ALB_P2-ALB_P1)/(D2-D1))*(D-D1)+ALB_P1
            ELSEIF( D.GE.D2 .AND. D.LT.D3 )THEN
              ALB_P = ((ALB_P3-ALB_P2)/(D3-D2))*(D-D2)+ALB_P2
            ELSEIF( D.GE.D3 .AND. D.LT.D4 )THEN
              ALB_P = ((ALB_P4-ALB_P3)/(D4-D3))*(D-D3)+ALB_P3
            ELSEIF( D.GE.D4 .AND. D.LT.D5 )THEN
              ALB_P = ((ALB_P5-ALB_P4)/(D5-D4))*(D-D4)+ALB_P4
            ELSE
              ALB_P = PARMS_P(10,IP)
            ENDIF
!
!---      Constant plant solar albedo  ---
!
          ELSE
            ALB_P = PARMS_P(10,IP)
          ENDIF
          RN_SNP(IP) = (RN+RN_SUS)*PAI(IP)*
     &      (1.D+0-EXP(-7.D-1*LAI(IP)))*(1.D+0-ALB_P)
  360   CONTINUE
!
!---    Gas tortuosity at the ground surface  ---
!
        CALL TORTU( IZ(N),SL_S,SG_S,ZERO,PORD(2,N),
     &    TORL_S,TORG_S,TORN_S )
!
!---    Water-vapor diffusion coefficient at the ground surface  ---
!
        IF( ISLC(2).EQ.1 ) THEN
          DFGW_S = DFGWC
        ELSEIF( ISLC(2).EQ.2 ) THEN
          CALL BNDFAW( T_SY,PG_S,DFGW_S )
        ELSEIF( ISLC(2).EQ.3 ) THEN
          CALL BNDFAW( T_SY,PG_S,DFGW_S )
          CMFF = 1.D+0 + 2.6D+0/(DFGWC**0.5)
          AMC = PORD(2,N)*SL_S
          ENHF = 9.5D+0 + 6.D+0*(AMC) -
     &      8.5D+0/EXP((CMFF*AMC)**4)
          DFGW_S = ENHF*DFGW_S
        ELSEIF( ISLC(2).EQ.4 ) THEN
          CALL BNDFAW( T_SY,PG_S,DFGW_S )
          ENHF = DFEF(1,IZ(N))+DFEF(2,IZ(N))*SL_S-
     &      (DFEF(1,IZ(N))-DFEF(4,IZ(N)))
     &      *EXP(-((DFEF(3,IZ(N))*SL_S)**DFEF(5,IZ(N))))
          DFGW_S = ENHF*DFGW_S
        ENDIF
!
!---    Effective thermal conductivity coefficient  ---
!
        CALL WATLQK( T_SY,PG_S,PSW_S,THKL_S )
        INDX = IBD
        TK_N = THKL(M,N)
        TK_S = THKE_L( IZ(N),SL_S,THKL_S,PORD(2,N),PORT(2,N),INDX )
        INDX = -1
        TK_NS = DIFMN(TK_N,TK_S,D_NS,D_NS,ZERO,INDX)
!
!---    Aqueous volumetric flux from the node to the
!       ground surface, m/s  ---
!
        HDL_S = PL_SX + D_NS*GBX*5.D-1*(RHOL(M,N)+RHOL_S)
        HDL_N = PL(M,N) + PATM
        HDL_NS = HDL_N - HDL_S
        INDX = -1
        VISL_NS = DIFMN(VISL(M,N),VISL_S,D_NS,D_NS,HDL_NS,INDX)
        INDX = -4
        RKL_NS = DIFMN(RKL(IBD,M,N),RKL_S(IBD),D_NS,D_NS,HDL_NS,INDX)
        L_NS = (RKL_NS*PERM(IBD,IZ(N))/VISL_NS)*(HDL_NS/D_NS)
!
!---    Gas volumetric flux from the node to the ground
!       surface, m/s  ---
!
        HDG_S = PG_S + D_NS*GBX*5.D-1*(RHOG(M,N)+RHOG_S)
        IF( 1.D+0-SL(M,N).LT.1.D-9 ) THEN
          HDG_N = HDG_S
        ELSE
          HDG_N = PG(M,N) + PATM
        ENDIF
        HDG_NS = HDG_N - HDG_S
        INDX = -1
        VISG_NS = DIFMN(VISG(M,N),VISG_S,D_NS,D_NS,HDG_NS,INDX)
        INDX = -4
        RKG_NS = DIFMN(RKG(M,N),RKG_S,D_NS,D_NS,HDG_NS,INDX)
!
!---    Gas flow from node to ground surface, m^3/s  ---
!
        G_NS = (RKG_NS*PERM(IBD,IZ(N))/VISG_NS)*((HDG_N-HDG_S)/(D_NS))
!
!---    Diffusive water-vapor mass flux from the node to
!       the surface, kg/m^2 s  ---
!
        DF_S = TORG_S*PORD(2,N)*SG_S*DFGW_S
        DF_N = TORG(M,N)*PORD(M,N)*(SG(M,N)-SGT(M,N))*DFGW(M,N)
        INDX = -1
        DF_NS = DIFMN(DF_N,DF_S,D_NS,D_NS,G_NS,INDX)
        E_N = RHOG(M,N)*XGW(M,N)
        E_NS = DF_NS*(E_N-E_S)/D_NS
!
!---    Successive substitution, following Campbell [1985],
!       to determine aerodynamic resistance  ---
!
        H_SA = (1.D+0-PAI_P)*(T_SX-TA)*RHOG_A*CPG_A/RA_BSX
        VARX = K_VK*ATMC(2)*GRAV/(RHOG_A*CPG_A*TA_K)
        ASP = -VARX*H_SA/(UZ_FX**3)
        DO 370 NCX = 1,4
!
!---      Surface temperature higher than atmospheric
!         temperature, unstable conditions  ---
!
          IF( (T_SX-TA).GT.EPSL ) THEN
            SPCFH = -2.D+0*LOG(5.D-1*(1.D+0+SQRT(
     &        MAX(1.D+0-(1.6D+1*ASP),0.D+0))))
            SPCFM = 6.D-1*SPCFH
!
!---      Surface temperature lower than atmospheric
!         temperature, stable conditions  ---
!
          ELSEIF( (T_SX-TA).LT.-EPSL ) THEN
            SPCFH = 4.7D+0*ASP
            SPCFM = 4.7D+0*ASP
!
!---      Surface temperature equal to atmospheric
!         temperature, skip calculations  ---
!
          ELSE
            SPCFH = 0.D+0
            SPCFM = 0.D+0
          ENDIF
!
!---      Corrections limited to 20%  ---
!
          SPCFH = SIGN( MIN( ABS(SPCFH),2.D-1*VLG_2 ),SPCFH )
          SPCFM = SIGN( MIN( ABS(SPCFM),2.D-1*VLG_1 ),SPCFM )
          UZ_F = K_VK*UZ/(VLG_1+SPCFM)
          RA_BS = MIN( (VLG_2+SPCFH)/(K_VK*UZ_F),RA_BSMX )
!
!---      Sensible heat flux from the ground surface to the
!         atmosphere, W/m^2  ---
!
          H_SA = (1.D+0-PAI_P)*(T_SX-TA)*RHOG_A*CPG_A/RA_BS
!
!---      Atmospheric stability parameter  ---
!
          ASP = -VARX*H_SA/(UZ_F**3)
  370   CONTINUE
!
!---    Diffusive water-vapor mass flux from the ground
!       surface to the atmosphere, kg/m^2 s  ---
!
        E_SA = (1.D+0-PAI_P)*(E_S-E_A)/RA_BS
        PE_SA = (1.D+0-PAI_P)*(PE_S-E_A)/RA_BS
!
!---    Enthalpy of water-vapor, air and gas at the surface, J/kg  ---
!
        CALL WATGSH( T_SX,PSW_S,RHOGW_S,HGW_S,UGW_S )
        CALL AIRGSH( T_SX,PGA_S,HGA_S,UGA_S )
        HG_S = XGW_S*HGW_S + XGA_S*HGA_S
        UG_S = XGW_S*UGW_S + XGA_S*UGA_S
!
!---    Enthalpy of liquid-water and aqueous at the
!       ground surface, J/kg  ---
!
        CALL WATLQH( T_SX,PGA_S,HLW_S )
        HL_S = XLW_S*HLW_S + XLA_S*HGA_S
!
!---    Enthalpy of gas at the node, J/kg  ---
!
        HG_N = XGW(M,N)*HGW(M,N) + XGA(M,N)*HGA(M,N)
!
!---    Upward long-wave radiation from the ground surface,
!       including reflected long-wave radiation from the sky
!       and emission from the ground surface, W/m^2 ground  ---
!
        TS_K = T_SX + TABS
        EM_S = 9.D-1 + 1.8D-1*(SL_S*PORD(M,N))
        RN_LUS = (1.D+0-EM_S)*RN_LDS + EM_S*SIGMA*(TS_K**4)
!
!---    Net long-wave radiation into the ground surface,
!       W/m^2 ground  ---
!
        RN_LNS = RN_LDS - RN_LUS
!
!---    Net long-wave radiation into the plants and
!       net downward long-wave radiation from plants to
!       ground surface, W/m^2 ground  ---
!
!       EM_P - Emissivity of plants [Llasat and Snyder, 1998]
!
        DO 380 IP = 1,NPLANT
          RN_LNP(IP) = 0.D+0
          IF( IPLANT(IP).EQ.0 ) GOTO 380
          TP_K = T_PX + TABS
          EM_P = 9.8D-1
          RN_LNP(IP) = PAI(IP)*(RN_LD + RN_LUS
     &      - 2.D+0*SIGMA*(TP_K**4))*EM_P*(1.D+0-EXP(-7.D-1*LAI(IP)))
          RN_LNS = RN_LNS + PAI(IP)*SIGMA*(TP_K**4)*EM_P
     &      *(1.D+0-EXP(-7.D-1*LAI(IP)))*EM_S
  380   CONTINUE
!
!---    Net long- and short-wave radiation into the
!       ground surface, W/m^2 ground  ---
!
        RN_NS = RN_LNS + RN_SNS
!
!---    Net long- and short-wave radiation into the
!       plants, W/m^2 ground  ---
!
        RN_NP = 0.D+0
        DO 390 IP = 1,NPLANT
          IF( IPLANT(IP).EQ.0 ) GOTO 390
          RN_NP = RN_NP + RN_LNP(IP) + RN_SNP(IP)
  390   CONTINUE
!
!---    Sensible heat flux from the field node to
!       the ground surface, W/m^2  ---
!
        CALL WATLQK( T_SY,PG_S,PSW_S,THKL_S )
        INDX = IBD
        TK_S = THKE_L( IZ(N),SL_S,THKL_S,PORD(M,N),PORT(M,N),INDX )
        TK_N = THKL(M,N)
        INDX = -1
        TK_NS = DIFMN(TK_N,TK_S,D_NS,D_NS,ZERO,INDX)
        H_NS = (T(M,N)-T_SX)*TK_NS/D_NS
!
!---    Plant condensate density  ---
!
        PG_P = PA
        CALL WATLQD( T_PY,PG_P,RHOL_P )
        INDX = 0
        CALL REGION_4( T_PX,PSW_P,INDX )
        PGW_P = PSW_P
        PGA_P = MAX( PG_P-PGW_P,0.D+0 )
        XMLA_P = PGA_P/HCAW
        XMLW_P = MAX( 1.D+0-XMLA_P,0.D+0 )
        XLW_P = XMLW_P*WTMW/(XMLA_P*WTMA + XMLW_P*WTMW)
        XLA_P = MAX( 1.D+0-XLW_P,0.D+0 )
        RHOLW_P = RHOL_P*XLW_P
!
!---    Component and phase density interfacial averages  ---
!
        RHOGA_N = RHOG(M,N)*XGA(M,N)
        RHOGW_N = RHOG(M,N)*XGW(M,N)
        RHOLW_N = RHOL(M,N)*XLW(M,N)
        RHOLW_S = RHOL_S*XLW_S
        RHOLW_A = RHOL_A*XLW_A
        INDX = -4
        RHOGW_SA = DIFMN(RHOGW_S,RHOGW_A,D_NS,D_NS,G_NS,INDX)
        RHOGW_NS = DIFMN(RHOGW_N,RHOGW_S,D_NS,D_NS,G_NS,INDX)
        RHOGA_SA = DIFMN(RHOGA_S,RHOGA_A,D_NS,D_NS,G_NS,INDX)
        RHOGA_NS = DIFMN(RHOGA_N,RHOGA_S,D_NS,D_NS,G_NS,INDX)
        RHOG_SA = DIFMN(RHOG_S,RHOG_A,D_NS,D_NS,G_NS,INDX)
        RHOG_NS = DIFMN(RHOG(M,N),RHOG_S,D_NS,D_NS,G_NS,INDX)
        RHOLW_NS = DIFMN(RHOLW_N,RHOLW_S,D_NS,D_NS,L_NS,INDX)
        RHOLW_SA = DIFMN(RHOLW_S,RHOLW_A,D_NS,D_NS,L_SA,INDX)
        RHOLW_PA = DIFMN(RHOLW_P,RHOLW_A,D_NS,D_NS,L_PAX,INDX)
!
!---    Gas flow from ground surface to atmosphere, m^3/s
!       using air mass balance at ground surface  ---
!
        G_SA = G_NS*RHOGA_NS/RHOGA_SA
!
!---    Component and phase interfacial averages  ---
!
        INDX = -4
        HGW_SA = DIFMN(HGW_S,HGW_A,D_NS,D_NS,E_SA,INDX)
        HGW_NS = DIFMN(HGW(M,N),HGW_S,D_NS,D_NS,E_NS,INDX)
        HGA_SA = DIFMN(HGA_S,HGA_A,D_NS,D_NS,E_SA,INDX)
        HGA_NS = DIFMN(HGA(M,N),HGA_S,D_NS,D_NS,E_NS,INDX)
        HG_SA = DIFMN(HG_S,HG_A,D_NS,D_NS,G_SA,INDX)
        HG_NS = DIFMN(HG_N,HG_S,D_NS,D_NS,G_NS,INDX)
        HL_NS = DIFMN(HL(M,N),HL_S,D_NS,D_NS,L_NS,INDX)
        HL_SA = DIFMN(HL_S,HL_A,D_NS,D_NS,L_SA,INDX)
        HL_PA = DIFMN(HL_P,HL_A,D_NS,D_NS,L_PAX,INDX)
!
!---    Water-vapor density at canopy height, kg/m^3  ---
!
        PGA_C = MAX( PA-PGW_CX,0.D+0 )
        XMLA_C = PGA_C/HCAW
        XMLW_C = MAX( 1.D+0-XMLA_C,0.D+0 )
        XLW_C = XMLW_C*WTMW/(XMLA_C*WTMA + XMLW_C*WTMW)
        XLA_C = MAX( 1.D+0-XLW_C,0.D+0 )
        INDX = 1
        CALL WATGSD( T_CX,PGW_CX,RHOGW_C,INDX )
        CALL AIRGSD( T_CX,PGA_C,RHOGA_C )
        RHOG_C = RHOGW_C+RHOGA_C
        XGW_C = RHOGW_C/RHOG_C
        XGA_C = RHOGA_C/RHOG_C
        E_C = RHOGW_C
!
!---    Water-vapor density at the plant, assuming
!       saturated vapor conditions inside the plant
!       leaf, kg/m^3  ---
!
        INDX = 0
        CALL REGION_4( T_PX,PSW_P,INDX )
        PGW_P = PSW_P
        PGA_P = MAX( PA-PGW_P,0.D+0 )
        XMLA_P = PGA_P/HCAW
        XMLW_P = MAX( 1.D+0-XMLA_P,0.D+0 )
        XLW_P = XMLW_P*WTMW/(XMLA_P*WTMA + XMLW_P*WTMW)
        XLA_P = MAX( 1.D+0-XLW_P,0.D+0 )
        INDX = 1
        CALL WATGSD( T_PX,PGW_P,RHOGW_P,INDX )
        CALL AIRGSD( T_PX,PGA_P,RHOGA_P )
        RHOG_P = RHOGW_P+RHOGA_P
        XGW_P = RHOGW_P/RHOG_P
        XGA_P = RHOGA_P/RHOG_P
        E_P = RHOGW_P
!
!---    Diffusive water-vapor mass flux from the ground
!       surface to the canopy height, kg/m^2 s  ---
!
        PE_SC = PAI_P*(PE_S-E_C)/RA_SC
        E_SC = PAI_P*(E_S-E_C)/RA_SC
!
!---    Sensible heat flux from the plant to the
!       canopy height, W/m^2  ---
!
        H_PC = 0.D+0
        DO 400 IP = 1,NPLANT
          IF( IPLANT(IP).EQ.0 ) GOTO 400
          H_PC = H_PC + PAI(IP)/RA_PC(IP)
  400   CONTINUE
        H_PC = H_PC*(T_PX-T_CX)*RHOG_A*CPG_A
!
!---    Sensible heat flux from the ground surface to the
!       canopy height, W/m^2  ---
!
        H_SC = PAI_P*(T_SX-T_CX)*RHOG_A*CPG_A/RA_SC
!
!---    Sensible heat flux from the canopy height
!       to the atmosphere, W/m^2  ---
!
        H_CA = PAI_P*(T_CX-TA)*RHOG_A*CPG_A/RA_CA
!
!---    Enthalpy of water-vapor, air and gas
!       at the canopy height, J/kg  ---
!
        CALL WATGSH( T_CX,PGW_C,RHOGW_C,HGW_C,UGW_C )
        CALL AIRGSH( T_CX,PGA_C,HGA_C,UGA_C )
        HG_C = XGW_C*HGW_C + XGA_C*HGA_C
!
!---    Enthalpy of water-vapor, air and gas at the plant, J/kg  ---
!
        CALL WATGSH( T_PX,PSW_P,RHOGW_P,HGW_P,UGW_P )
        CALL AIRGSH( T_PX,PGA_P,HGA_P,UGA_P )
        HG_P = XGW_P*HGW_P + XGA_P*HGA_P
!
!---    Enthalpy of liquid-water at the plant, J/kg  ---
!
!       HLW_P(1,IP,NB) - enthalpy of liquid-water at old time step,
!                        for plant IP, boundary condition NB, J/kg
!       HLW_P(2,IP,NB) - enthalpy of liquid-water at current time,
!                        for plant IP, boundary condition NB, J/kg
!       HLW_PX(1) - enthalpy of liquid-water at old time step,
!                   for all plants, boundary condition NB, J/kg
!       HLW_PX(2) - enthalpy of liquid-water at current time,
!                   for all plants, boundary condition NB, J/kg
!
        DO 410 IP = 1,NPLANT
          IF( IPLANT(IP).EQ.0 ) GOTO 410
          HLW_PX(1) = HLW_P(1,IP,NB)
          CALL WATLQH( T_PX,PGA_P,HLW_P(2,IP,NB) )
          HLW_PX(2) = HLW_P(2,IP,NB)
  410   CONTINUE
!
!---    Enthalpy of liquid-water at the node, J/kg  ---
!
        HLW_N = (HL(M,N)-XLA(M,N)*HGA(M,N))/XLW(M,N)
!
!---    Diffusive water-vapor mass flux from the plant
!       to the canopy height, using aerodynamic and
!       stomatal resistances; attenuated by the crop
!       growth index and normalized stressed root uptake,
!       kg/m^2 s  ---
!
        E_PCX = 0.D+0
        PE_PCX = 0.D+0
        T_PCX = 0.D+0
        PT_PCX = 0.D+0
        RFIM_PX(1) = 0.D+0
        RFIM_PX(2) = 0.D+0
        L_SPX = 0.D+0
        DO 420 IP = 1,NPLANT
          IF( IPLANT(IP).EQ.0 ) GOTO 420
          DE_PC = E_P-E_C
!
!---      Rainfall interception mass on plant leaves from
!         rainfall, kg/m^2 ground
!
!         RFIM_P(2,IP,NB) - rainfall interception mass on plant
!                           leaves at current time for plant IP,
!                           for boundary NB, kg/m^2 ground
!         RFIM_P(1,IP,NB) - rainfall interception mass on plant
!                           leaves at previous time step for plant IP,
!                           for boundary NB, kg/m^2 ground
!         DT - time step, s
!         L_PA(IP) - rainfall intensity on plant IP, m^3/s m^2 ground
!         RHOL_A - density of atmospheric water, kg/m^3
!
          RFIM_P(2,IP,NB) = RFIM_P(1,IP,NB) - L_PA(IP)*RHOL_A*DT
!
!---      Evaporation and transpiration from the plant
!         to the canopy  ---
!
          IF( DE_PC.GT.0.D+0 ) THEN
!
!---        Transpiration, including stomatal resistance and
!           root-water uptake stress ---
!
            PT_PC(IP) = PAI(IP)*DE_PC/(RA_PC(IP)+RS_PC(IP))
            PT_PCX = PT_PCX + PT_PC(IP)
            T_PC(IP) = PT_PC(IP)*C_COEF(IP)*RWU_S(IP)/RWU_I(IP)
            T_PCX = T_PCX + T_PC(IP)
!!
!!---        Evaporation of rainfall interception or condensate water ---
!!
!            IF( RFIM_P(2,IP,NB).GT.EPSL ) THEN
!              E_LC = PAI(IP)*DE_PC/RA_PC(IP)
!              IF( E_LC.GT.(RFIM_P(2,IP,NB)/DT) ) THEN
!                E_LC = RFIM_P(2,IP,NB)/DT
!                RFIM_P(2,IP,NB) = 0.D+0
!              ELSE
!                RFIM_P(2,IP,NB) = RFIM_P(2,IP,NB)-E_LC*DT
!              ENDIF
!              E_PC(IP) = E_LC
!              E_PCX = E_PCX + E_PC(IP)
!            ENDIF
!!
!!---      Condensation on the plant surface,
!!         no stomatal resistance, no root water uptake stress,
!!         moderated condensation rate near the dew point  ---
!!
!          ELSEIF( T_PX.LE.T_DP ) THEN
!            PE_PC(IP) = PAI(IP)*DE_PC*(1.D+0-EXP(T_PX-T_DP))
!     &        /RA_PC(IP)
!            PE_PCX = PE_PCX + PE_PC(IP)
!            E_PC(IP) = PE_PC(IP)
!            E_PCX = E_PCX + E_PC(IP)
!            T_PC(IP) = 0.D+0
!!
!!---        Rainfall interception mass on plant leaves from
!!           condensation, kg/m^2 ground
!!
!!           RFIM_P(2,IP,NB) - rainfall interception mass on plant
!!                             leaves at current time for plant IP,
!!                             for boundary NB, kg/m^2 ground
!!           RFIM_P(1,IP,NB) - rainfall interception mass on plant
!!                             leaves at previous time step for plant IP,
!!                             for boundary NB, kg/m^2 ground
!!           E_PC(IP) - condensation rate on plant IP, kg/s m^2 ground
!!
!            RFIM_P(2,IP,NB) = RFIM_P(2,IP,NB) - E_PC(IP)*DT
!!
!!---      Condensation within the plant,
!!         no root water uptake stress  ---
!!
!          ELSE
!            PT_PC(IP) = PAI(IP)*DE_PC/(RA_PC(IP)+RS_PC(IP))
!            PT_PCX = PT_PCX + PT_PC(IP)
!            T_PC(IP) = PT_PC(IP)
!            T_PCX = T_PCX + T_PC(IP)
          ENDIF
!
!---      Rainfall and condenstation greater than the maximum
!         rainfall interception capacity, water is shed
!         to the ground surface, kg/s m^2 ground
!
!         L_SP(IP) - shed rainfall and condensation from plant IP,
!         kg/s m^2 ground
!
          IF( RFIM_P(2,IP,NB).GT.RFIC_P(IP) ) THEN
            L_SP(IP) = -(RFIM_P(2,IP,NB)-RFIC_P(IP))/DT
            RFIM_P(2,IP,NB) = RFIC_P(IP)
          ELSE
            L_SP(IP) = 0.D+0
          ENDIF
!
!---      Sum over plants, shed rainfall and condensate
!         mass kg/m^2 ground  ---
!
          L_SPX = L_SPX + L_SP(IP)
!
!---      Sum rainfall interception mass and shed water mass over
!         plant species
!
          RFIM_PX(1) = RFIM_PX(1) + RFIM_P(1,IP,NB)
          RFIM_PX(2) = RFIM_PX(2) + RFIM_P(2,IP,NB)
  420   CONTINUE
!
!---    Diffusive water-vapor mass flux from the canopy height
!       to the atmosphere, kg/m^2 s  ---
!
        E_CA = PAI_P*(E_C-E_A)/RA_CA
!
!---    Component and phase interfacial averages  ---
!
        INDX = -4
        HGW_SC = DIFMN(HGW_S,HGW_C,ONE,ONE,E_SC,INDX)
        HGWE_PC = DIFMN(HGW_P,HGW_C,ONE,ONE,E_PCX,INDX)
        HGWT_PC = DIFMN(HGW_P,HGW_C,ONE,ONE,T_PCX,INDX)
        HGW_CA = DIFMN(HGW_C,HGW_A,ONE,ONE,E_CA,INDX)
        HLW_SP = DIFMN(HLW_S,HLW_PX(2),ONE,ONE,L_SPX,INDX)
        HLW_NP = DIFMN(HLW_N,HLW_PX(2),ONE,ONE,T_PCX,INDX)
!
!---    Ponding height limit  ---
!
        L_SAPMX = -PERM(3,IZ(N))*(RHOLW_SA**2)*GRAV*
     &    (1.D+0-SWBCD/DB)/VISL_S
        L_SAPX = L_SA*RHOLW_SA + L_SPX - RWRO(1)*RHOLW_SA/A_SA
        IF( ABS(L_SAPX).GT.ABS(L_SAPMX) ) THEN
          L_SAPX = L_SAPMX
        ENDIF
!
!---    Residual of the water balance at the canopy height,
!       kg/s m^2 ground  ---
!
        RW_C(MP) = E_SC + E_PCX + T_PCX - E_CA
!
!---    Residual of the energy balance at the canopy height,
!       W/m^2 ground  ---
!
        RE_C(MP) = H_SC + E_SC*HGW_SC
     &    + H_PC + E_PCX*HGWE_PC + T_PCX*HGWT_PC
     &    - H_CA - E_CA*HGW_CA
!
!---    Residual of the energy balance at the plant,
!       W/m^2 ground  ---
!
!       RN_NP - solar radiation, W/m^2 ground
!       L_PAX*RHOLW_PA*HL_PA - rainfall interception, W/m^2 ground
!       T_PCX*HLW_NP - rootwater uptake, W/m^2 ground
!       (RFIM_PX(2)*HLW_PX(2)-RFIM_PX(1)*HLW_PX(1))/DT - rainfall
!         interception and condensate water storage, , W/m^2 ground
!       H_PC - advective heat transport, W/m^2 ground
!       E_PCX*HGWE_PC - evaporation of rainfall interception and
!         condensate water, , W/m^2 ground
!       T_PCX*HGWT_PC - transpiration, W/m^2 ground
!       L_SPX*HLW_SP - shed rainfall and condensate, W/m^2 ground
!
        RE_P(MP) = RN_NP - L_PAX*RHOLW_PA*HL_PA + T_PCX*HLW_NP
     &    - (RFIM_PX(2)*HLW_PX(2)-RFIM_PX(1)*HLW_PX(1))/DT
     &    - H_PC - E_PCX*HGWE_PC - T_PCX*HGWT_PC + L_SPX*HLW_SP
!
!---    Residual of the water balance at the ground surface,
!       kg/s m^2 ground  ---
!
!        RW_S(MP) = E_NS + G_NS*RHOGW_NS + L_NS*RHOLW_NS
!     &    - E_SA - G_SA*RHOGW_SA - L_SA*RHOLW_SA - E_SC
!     &    - L_SPX
        RW_S(MP) = E_NS + G_NS*RHOGW_NS + L_NS*RHOLW_NS
     &    - E_SA - G_SA*RHOGW_SA - L_SAPX - E_SC
!
!---    Residual of the energy balance at the ground surface,
!       ignoring dissolved air mass in aqueous,
!       W/m^2 ground  ---
!
!        RE_S(MP) = RN_NS + E_NS*HGW_NS + G_NS*RHOG_NS*HG_NS
!     &    + L_NS*RHOLW_NS*HL_NS + H_NS
!     &    - E_SA*HGW_SA - G_SA*RHOG_SA*HG_SA
!     &    - L_SA*RHOLW_SA*HL_SA - H_SA
!     &    - H_SC - E_SC*HGW_SC - L_SPX*HLW_SP
        RE_S(MP) = RN_NS + E_NS*HGW_NS + G_NS*RHOG_NS*HG_NS
     &    + L_NS*RHOLW_NS*HL_NS + H_NS
     &    - E_SA*HGW_SA - G_SA*RHOG_SA*HG_SA
     &    - L_SAPX*HL_SA - H_SA
     &    - H_SC - E_SC*HGW_SC
  430 CONTINUE
      IF( ABS(MPX).EQ.1 ) GOTO 470
!
!---  Jacobian coefficient matrix  ---
!
      IF( MPX.GT.0 ) THEN
        BJ(1) = -RE_P(1)
        BJ(2) = -RE_C(1)
        BJ(3) = -RW_C(1)
        BJ(4) = -RW_S(1)
        BJ(5) = -RE_S(1)
        DO 440 JP = 1,ABS(MPX)-1
          AJ(1,JP) = (RE_P(JP+1)-RE_P(1))/DNRX(JP)
          AJ(2,JP) = (RE_C(JP+1)-RE_C(1))/DNRX(JP)
          AJ(3,JP) = (RW_C(JP+1)-RW_C(1))/DNRX(JP)
          AJ(4,JP) = (RW_S(JP+1)-RW_S(1))/DNRX(JP)
          AJ(5,JP) = (RE_S(JP+1)-RE_S(1))/DNRX(JP)
  440   CONTINUE
      ELSE
        BJ(1) = -RE_P(1)
        BJ(2) = -RE_C(1)
        BJ(3) = -RW_C(1)
        BJ(4) = -RE_S(1)
        BJ(5) = -RW_S(1)
        DO 450 JP = 1,ABS(MPX)-1
          AJ(1,JP) = (RE_P(JP+1)-RE_P(1))/DNRX(JP)
          AJ(2,JP) = (RE_C(JP+1)-RE_C(1))/DNRX(JP)
          AJ(3,JP) = (RW_C(JP+1)-RW_C(1))/DNRX(JP)
          AJ(4,JP) = (RE_S(JP+1)-RE_S(1))/DNRX(JP)
          AJ(5,JP) = (RW_S(JP+1)-RW_S(1))/DNRX(JP)
  450   CONTINUE
      ENDIF
!
!---  Solve linear system  ---
!
      JP = ABS(MPX)-1
      CALL LUDCMP( AJ,JP,KP,IJ,DJ )
      CALL LUBKSB( AJ,JP,KP,IJ,BJ )
!
!---  Limit changes in surface saturation  ---
!
      IF( MPX.GT.0 ) THEN
        DP_S = MAX( 1.D+5,ABS(PG_S-PL_S) )
        DP_S = MIN( DP_S,ABS(BJ(4)) )
        DP_S = SIGN( DP_S,BJ(4) )
      ELSE
        DP_S = 0
      ENDIF
!
!---  Surface saturation above residual  ---
!
      IF( SL_SX.GT.SCHR(4,IZ(N)) ) THEN
!
!---    Limit rising saturation changes to 0.05  ---
!
        IF( DP_S.GT.0.D+0 .AND. SL_SX.LT.9.5D-1 ) THEN
          CALL CAP_GT( IZ(N),SL_SX,CPGLX )
          SLX = SL_SX+5.D-2
          CALL CAP_GT( IZ(N),SLX,CPGL )
          DP_S = MIN( DP_S,(CPGLX-CPGL) )
!
!---    Limit falling saturation changes to 0.05  ---
!
        ELSEIF( DP_S.LT.0.D+0 ) THEN
          CALL CAP_GT( IZ(N),SL_SX,CPGLX )
          SLX = MAX( SL_SX-5.D-2,9.5D-1*SCHR(4,IZ(N)),5.D-1*SL_SX )
          CALL CAP_GT( IZ(N),SLX,CPGL )
!
!---      Transition from saturated conditions  ---
!
          CPGL = MAX( PL_S-PG_S+CPGL,CPGL )
          DP_S = MAX( DP_S,(CPGLX-CPGL) )
        ENDIF
!
!---  Surface saturation below residual  ---
!
      ELSE
!
!---    Limit rising saturation changes to 0.05
!       of the residual saturation  ---
!
        IF( DP_S.GT.0.D+0 ) THEN
          CALL CAP_GT( IZ(N),SL_SX,CPGLX )
          SLX = SL_SX + 5.D-2*SCHR(4,IZ(N))
          CALL CAP_GT( IZ(N),SLX,CPGL )
          DP_S = MIN( DP_S,(CPGLX-CPGL) )
!
!---    Limit falling saturation changes to 0.05
!       of the residual saturation  ---
!
        ELSEIF( DP_S.LT.0.D+0 ) THEN
          CALL CAP_GT( IZ(N),SL_SX,CPGLX )
          SLX = MAX( SL_SX-5.D-2*SCHR(4,IZ(N)),5.D-1*SL_SX )
          CALL CAP_GT( IZ(N),SLX,CPGL )
          DP_S = MAX( DP_S,(CPGLX-CPGL) )
        ENDIF
      ENDIF
!
!---  Limit changes in temperature at the ground surface  ---
!
      IF( MPX.EQ.6 ) THEN
        DT_S = MIN( 5.D-1,ABS(BJ(5)) )
        DT_S = SIGN( DT_S,BJ(5) )
      ELSEIF( MPX.EQ.-5 ) THEN
        DT_S = MIN( 5.D-1,ABS(BJ(4)) )
        DT_S = SIGN( DT_S,BJ(4) )
      ELSE
        DT_S = 0.D+0
      ENDIF
!
!---  Limit changes in temperature at the plant  ---
!
      DT_P = MIN( 5.D-1,ABS(BJ(1)) )
      DT_P = SIGN( DT_P,BJ(1) )
!
!---  Limit changes in temperature at the canopy height  ---
!
      DT_C = MIN( 5.D-1,ABS(BJ(2)) )
      DT_C = SIGN( DT_C,BJ(2) )
!
!---  Update water-vapor partial pressure at the canopy
!     height  ---
!
      DPGW_C = MIN( 1.D-1*PSW_A,ABS(BJ(3)) )
      DPGW_C = SIGN( DPGW_C,BJ(3) )
!
!---  Relax solution for flopping surface pressure conditions
!     or transitions across the freezing point  ---
!
      IFLOP = 0
      IF( (DP_SO*DP_S).LT.0.D+0 .AND.
     &  ABS( 1.D+0-ABS(DP_SO/(DP_S+EPSL)) ).LT.1.D-1 .AND.
     &  (DP_SOO*DP_SO).LT.0.D+0 .AND.
     &  ABS( 1.D+0-ABS(DP_SOO/(DP_SO+EPSL)) ).LT.1.D-1) THEN
        IFLOP = 1
      ENDIF
      IF( IFLOP.EQ.1 ) THEN
        PL_S = PL_S + 6.D-1*DP_S
        T_S = T_S + 6.D-1*DT_S
        T_P = T_P + 6.D-1*DT_P
        T_C = T_C + 6.D-1*DT_C
        PGW_C = PGW_C + 6.D-1*DPGW_C
      ELSE
        PL_S = PL_S + DP_S
        T_S = T_S + DT_S
        T_P = T_P + DT_P
        T_C = T_C + DT_C
        PGW_C = PGW_C + DPGW_C
      ENDIF
      PL_S = MIN( PL_S,PG_S+SWBCD*RHOL_A*GRAV )
      PGW_C = MAX( PGW_C,0.D+0 )
!
!---  Store correction to surface pressure  ---
!
      DP_SOO = DP_SO
      DP_SO = DP_S
!
!---  Maximum residual  ---
!
      RSDEPX = 1.D-20
      RSDECX = 1.D-20
      RSDWCX = 1.D-20
      RSDWSX = 1.D-20
      RSDESX = 1.D-20
      RSDEPX = ABS(BJ(1))/TABS
      RSDECX = ABS(BJ(2))/TABS
      RSDWCX = ABS(BJ(3))/PATM
      IF( MPX.EQ.6 ) THEN
        RSDWSX = ABS(BJ(4))/PATM
        RSDESX = ABS(BJ(5))/TABS
      ENDIF
      IF( MPX.EQ.5 ) THEN
        CALL KSP_GT( N,IZ(N),M,PG_S,PL_S,SLX,SG_S,RKL_S,RKG_S )
        RSDWX = MIN( ABS(BJ(4))/PATM,ABS(SL_S-SLX)*1.D-2 )
      ENDIF
      IF( MPX.EQ.-5 ) RSDESX = ABS(BJ(4))/TABS
      RSDX = MAX( RSDEPX,RSDECX,RSDWCX,RSDWSX,RSDESX )
      IF( NC.GT.24 ) THEN
        T_SMN = T_SMN + T_S
        PL_SMN = PL_SMN + PL_S
        SL_SMN = SL_SMN + SL_S
        T_PMN = T_PMN + T_P
        T_CMN = T_CMN + T_C
        PGW_CMN = PGW_CMN + PGW_C
      ENDIF
!
!---  No convergence on water mass and energy
!     balances ---
!
      IF( NC.EQ.48 ) THEN
        REALX = REAL(NC-24)
        SL_S = SL_SMN/REALX
!
!---    Low residual, select the average temperature
!       and aqueous pressure ---
!
        IF( RSDX.LT.1.D-4 ) THEN
          IF( SL_S.LT.9.5D-1 ) THEN
            SL_S = SL_SMN/REALX
            CALL CAP_GT( IZ(N),SL_S,CPGLX )
            PL_S = PG_S - CPGLX
          ELSE
            PL_S = PL_SMN/REALX
          ENDIF
          PL_S = MIN( PL_S,PG_S+SWBCD*RHOL_A*GRAV )
          T_S = T_SMN/REALX
          T_P = T_PMN/REALX
          T_C = T_CMN/REALX
          PGW_C = PGW_CMN/REALX
          GOTO 460
!
!---    Fix surface aqueous pressure ---
!
        ELSEIF( SL_S.LT.5.D-1 ) THEN
          IF( M.EQ.2 ) THEN
            SL_S = SL_SMN/REALX
            CALL CAP_GT( IZ(N),SL_S,CPGLX )
            PL_S = PG_S - CPGLX
            PL_S = MIN( PL_S,PG_S+SWBCD*RHOL_A*GRAV )
          ENDIF
          T_S = T_SMN/REALX
          T_P = T_PMN/REALX
          T_C = T_CMN/REALX
          PGW_C = PGW_CMN/REALX
          MPX = -5
!!
!!---    Fix surface aqueous pressure ---
!!
!        ELSEIF( SL_S.GT.9.5D-1 ) THEN
!          IF( M.EQ.2 ) THEN
!            PL_S = PL_SMN/REALX
!            PL_S = MIN( PL_S,PG_S+SWBCD*RHOL_A*GRAV )
!          ENDIF
!          T_S = T_SMN/REALX
!          T_P = T_PMN/REALX
!          T_C = T_CMN/REALX
!          PGW_C = PGW_CMN/REALX
!          MPX = -5
!
!---    Fix surface temperature  ---
!
        ELSE
          IF( M.EQ.2 ) THEN
            T_S = T_SMN/REALX
          ENDIF
          IF( SL_S.LT.9.5D-1 ) THEN
            SL_S = SL_SMN/REALX
            CALL CAP_GT( IZ(N),SL_S,CPGLX )
            PL_S = PG_S - CPGLX
          ELSE
            PL_S = PL_SMN/REALX
          ENDIF
          PL_S = MIN( PL_S,PG_S+SWBCD*RHOL_A*GRAV )
          T_P = T_PMN/REALX
          T_C = T_CMN/REALX
          PGW_C = PGW_CMN/REALX
          MPX = 5
        ENDIF
      ENDIF
      IF( NC.EQ.72 ) THEN
        REALX = REAL(NC-24)
        T_P = T_PMN/REALX
        T_C = T_CMN/REALX
        PGW_C = PGW_CMN/REALX
        WRITE(ISC,'(10X,A)') '---  Shuttleworth-Wallace Boundary ' //
     &    'Condition Failure: Single Plant: Used Averaged Values  ---'
        WRITE(IWR,'(10X,A)') '---  Shuttleworth-Wallace Boundary ' //
     &    'Condition Failure: Single Plant: Used Averaged Values  ---'
        IF( MPX.EQ.5 ) THEN
          WRITE(ISC,'(10X,A)') 'Fixed Surface Temperature'
          WRITE(IWR,'(10X,A)') 'Fixed Surface Temperature'
          SL_S = SL_SMN/REALX
          CALL CAP_GT( IZ(N),SL_S,CPGLX )
          PL_S = PG_S - CPGLX
          PL_S = MIN( PL_S,PG_S+SWBCD*RHOL_A*GRAV )
        ELSEIF( MPX.EQ.-5 ) THEN
          WRITE(ISC,'(10X,A)') 'Fixed Surface Pressure'
          WRITE(IWR,'(10X,A)') 'Fixed Surface Pressure'
          T_S = T_SMN/REALX
        ENDIF
        WRITE(ISC,'(10X,A,I6,A,I6)') 'Node = ',N,'Boundary = ',NB
        WRITE(IWR,'(10X,A,I6,A,I6)') 'Node = ',N,'Boundary = ',NB
        WRITE(ISC,'(10X,A,1PE12.5)') 'Surface Pressure, Pa = ',
     &    (PL_S+PATM)
        WRITE(IWR,'(10X,A,1PE12.5)') 'Surface Pressure, Pa = ',
     &    (PL_S+PATM)
        WRITE(ISC,'(10X,A,1PE12.5)') 'Surface Temperature, C = ',T_S
        WRITE(IWR,'(10X,A,1PE12.5)') 'Surface Temperature, C = ',T_S
        WRITE(IWR,'(10X,A,1PE12.5)') 'Surface Saturation = ',SL_S
        WRITE(ISC,'(10X,A,1PE12.5)') 'Surface Saturation = ',SL_S
        WRITE(ISC,'(10X,A,1PE12.5)') 'Plant Temperature, C = ',T_P
        WRITE(IWR,'(10X,A,1PE12.5)') 'Plant Temperature, C = ',T_P
        WRITE(ISC,'(10X,A,1PE12.5)') 'Crop Temperature, C = ',T_C
        WRITE(IWR,'(10X,A,1PE12.5)') 'Crop Temperature, C = ',T_C
        WRITE(ISC,'(10X,A,1PE12.5)') 'Crop Vapor Pressure, Pa = ',PGW_C
        WRITE(IWR,'(10X,A,1PE12.5)') 'Crop Vapor Pressure, Pa = ',PGW_C
        GOTO 460
      ENDIF
!
!---  Check solution convergence  ---
!
      IF( RSDX.GE.1.D-6 ) GOTO 350
  460 CONTINUE
!
!---  Convergence on mass and energy balances at the
!     ground surface  ---
!
      T_SO(2,NB) = T_S
      PL_SO(2,NB) = PL_S
      T_PO(2,1,NB) = T_P
      T_CO(2,1,NB) = T_C
      PGW_CO(2,1,NB) = PGW_C
!
!---  Converged solution, calculate fluxes  ---
!
      MPOX = MPX
      MPX = 1
      NC = NC-1
      GOTO 350
  470 CONTINUE
      MPX = MPOX
      IF( M.EQ.2 ) THEN
!
!---    Load variables for output  ---
!
!       RHON(1,N) - Water mass balance at ground surface, kg/s
!       RHON(2,N) - Actual evaporation rate, kg/s
!       RHON(3,N) - Potential evaporation rate, kg/s
!       RHON(4,N) - Actual transpiration rate, kg/s
!       RHON(5,N) - Potential transpiration rate, kg/s
!
        RHON(1,N) = (E_NS+G_NS*RHOGW_NS+L_NS*RHOLW_NS)*A_SA
        RHON(2,N) = (E_SC + E_SA + G_SA*RHOGW_SA)*A_SA
        RHON(3,N) = (PE_SC + PE_SA + G_SA*RHOGW_SA)*A_SA
        RHON(4,N) = T_PCX*A_SA
        RHON(5,N) = PT_PCX*A_SA
!
!---    Load variables for output
!
!       SN(1,N) - Surface-atmosphere latent heat flux, W/m^2
!       SN(2,N) - Surface-atmosphere sensible heat flux, W/m^2
!       SN(3,N) - Surface-atmosphere net long-wave radiation, W/m^2
!       SN(4,N) - Surface-atmosphere net short-wave radiation, W/m^2
!       SN(5,N) - Ground-surface sensible and latent heat flux, W/m^2
!
        SN(1,N) = E_SA*(HGW_SA-HL_S)
     &    + G_SA*RHOGW_SA*(HGW_SA-HL_S)
        IF( ABS(L_SA).GT.EPSL ) THEN
          SN(1,N) = SN(1,N)
     &      + L_SA*(L_SAPX/(L_SA*RHOLW_SA+L_SPX))*RHOLW_SA*(HL_SA-HL_S)
        ENDIF
        SN(2,N) = G_SA*RHOGA_SA*(HGA_SA-HGA_S) + H_SA
        SN(3,N) = RN_LNS
        SN(4,N) = RN_SNS
        SN(5,N) = E_NS*(HGW_NS-HL_S)
     &    + L_NS*RHOLW_NS*(HL_NS-HL_S)
     &    + G_NS*RHOGW_NS*(HGW_NS-HL_S)
     &    + G_NS*RHOGA_NS*(HGA_NS-HGA_S) + H_NS
!
!---    Load variables for output  ---
!
!       XLO(1,N) - Surface temperature, C
!       XLO(2,N) - Surface aqueous pressure, Pa
!       XLO(3,N) - Surface gas pressure, Pa
!       XLO(4,N) - Surface aqueous saturation
!       XLO(5,N) - Surface water-vapor pressure, Pa
!
        XLO(1,N) = T_S
        XLO(2,N) = PL_S
        XLO(3,N) = PG_S
        XLO(4,N) = SL_S
        XLO(5,N) = PGW_S
!
!---    Load variables for output  ---
!
!       PN(1,N) - Atmospheric temperature, C
!       PN(2,N) - Atmospheric pressure, Pa
!       PN(3,N) - Atmospheric relative humidity
!       PN(4,N) - Atmospheric wind speed, m/s
!       PN(5,N) - Surface solar radiation, W/m^2
!
        PN(1,N) = TA
        PN(2,N) = PA
        PN(3,N) = RH
        PN(4,N) = UZ
        PN(5,N) = RN
!
!---    Load variables for output  ---
!
!       XMLO(1,N) - Plant temperature, C
!
        XMLO(1,N) = T_P
!
!---    Load variables for output  ---
!
!       XGO(1,N) - Rain-water volumetric runoff rate, m^3/s
!       XGO(2,N) - Aerodynamic resistance between ground and
!                   reference height for bare surface, s/m
!       XGO(3,N) - Surface volumetric precipitation rate, m^3/s
!       XGO(4,N) - Surface mass precipitation rate, kg/s
!       XGO(5,N) - Atmospheric mass precipitation rate, kg/s
!!       XGO(5,N) - Atmospheric stability parameter
!
        RWROX = (L_SAPX/RHOLW_SA - L_SA - L_SPX/RHOLW_SA)*A_SA
        RWRO(2) = RWRO(1) + RWROX
        IF( ABS(RWRO(2)).LT.1.D-20 ) RWRO(2) = 0.D+0
        XGO(1,N) = RWROX
        XGO(2,N) = RA_BS
        XGO(3,N) = (L_SA + L_SPX/RHOLW_SA)*A_SA
        XGO(4,N) = (L_SA*RHOLW_SA + L_SPX)*A_SA
        XGO(5,N) = -BCX(2)*RHOLW_SA*A_SA
!        XGO(5,N) = -K_VK*ATMC(2)*GRAV*H_NS/
!     &    (RHOG_A*CPG_A*TA_K*(UZ_F**3))
!
!---    Surface area  ---
!
        IF( IBCD(NB).EQ.-3 ) THEN
          A_SA = AFZ(NSZ(N))
        ELSEIF( IBCD(NB).EQ.-2 ) THEN
          A_SA = AFY(NSY(N))
        ELSEIF( IBCD(NB).EQ.-1 ) THEN
          A_SA = AFX(NSX(N))
        ELSEIF( IBCD(NB).EQ.1 ) THEN
          A_SA = AFX(NSX(N)+1)
        ELSEIF( IBCD(NB).EQ.2 ) THEN
          A_SA = AFY(NSY(N)+IFLD)
        ELSEIF( IBCD(NB).EQ.3 ) THEN
          A_SA = AFZ(NSZ(N)+IJFLD)
        ENDIF
      ENDIF
!
!---  Load boundary condition variables
!
      TB(M,NB) = T_S
      PLB(M,NB) = PL_S - PATM
      PGB(M,NB) = PG_S - PATM
      PVWB(M,NB) = PGW_S
      PVAB(M,NB) = PGA_S
      SLB(M,NB) = SL_S
      RKLB(1,M,NB) = RKL_S(1)
      RKLB(2,M,NB) = RKL_S(2)
      RKLB(3,M,NB) = RKL_S(3)
      SGB(M,NB) = SG_S
      RKGB(M,NB) = RKG_S
      XLWB(M,NB) = XLW_S
      XLAB(M,NB) = XLA_S
      XMLWB(M,NB) = XMLW_S
      XMLAB(M,NB) = XMLA_S
      XGWB(M,NB) = XGW_S
      XGAB(M,NB) = XGA_S
      XMGWB(M,NB) = XMGW_S
      XMGAB(M,NB) = XMGA_S
      RHOLB(M,NB) = RHOL_S
      RHOGB(M,NB) = RHOG_S
      HGWB(M,NB) = HGW_S
      HGAB(M,NB) = HGA_S
      HLB(M,NB) = HL_S
      UEGB(M,NB) = UG_S
!
!---  Impose boundary conditions through sources
!     at the boundary node  ---
!
!     SRCT(M,N) - thermal energy rate into the node, W
!     SRCA(M,N) - air mass rate into the node, kg/s
!     SRCW(M,N) - water mass rate into the node, kg/s
!
      SRCT(M,N) = SRCT(M,N)-(E_NS*HGW_NS + G_NS*RHOG_NS*HG_NS
     &    + L_NS*RHOLW_NS*HL_NS + H_NS)*A_SA
      SRCA(M,N) = SRCA(M,N)-(G_NS*RHOGA_NS)*A_SA
      SRCW(M,N) = SRCW(M,N)-(E_NS + G_NS*RHOGW_NS + L_NS*RHOLW_NS)*A_SA
!
!---  Remove root water and associated energy, adjusting
!     root-water enthalpy at each node  ---
!
      DO 490 IP = 1,NPLANT
        IF( IPLANT(IP).EQ.0 ) GOTO 490
        ZM = PARMS_P(1,IP)
        RWFX = 0.D+0
        DO 480 K = KMN(IP),KFLD
          NPX = ND(I,J,K)
          IF( IXP(NPX).EQ.0 ) GOTO 480
          IF( K.EQ.KMN(IP) ) THEN
            SRCWX = -MAX( RWF(IP,K)*T_PC(IP)*A_SA,0.D+0 )
            SRCW(M,NPX) = SRCW(M,NPX) + SRCWX
            SRCT(M,NPX) = SRCT(M,NPX) + SRCWX*HL(2,NPX)
          ELSE
            NBX = ND(I,J,K-1)
            SRCWX = -MAX( RWF(IP,K)*T_PC(IP)*A_SA,0.D+0 )
            SRCW(M,NPX) = SRCW(M,NPX) + SRCWX
            SRCT(M,NPX) = SRCT(M,NPX) + SRCWX*HL(2,NPX)
     &        + MAX( RWFX*T_PC(IP)*A_SA,0.D+0 )*(HL(2,NBX)-HL(2,NPX))
          ENDIF
          RWFX = MIN( RWFX+RWF(IP,K),1.D+0 )
  480   CONTINUE
  490 CONTINUE
!
!---  Archive unknowns  ---
!
      IF( M.EQ.2 ) THEN
        T_SO(2,NB) = T_S
        PL_SO(2,NB) = PL_S
        T_PO(2,1,NB) = T_P
        T_CO(2,1,NB) = T_C
        PGW_CO(2,1,NB) = PGW_C
      ENDIF
 1000 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SNGL_T group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SOLR_C( RN,RN_DB,RN_ID,S_ZA,J_DAY )
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
!     Split total horizontal solar radiation into direct-beam and
!     indirect components.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 10 September 2002.
!     Last Modified by MD White, PNNL, 10 September 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PLT_ATM
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Data Statements---------------------------------!
!
      SAVE G_SC
      DATA G_SC / 1.367D+3 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SOLR_C'
!
!---  Initialize direct beam and indirect solar radiation  ---
!
      RN_DB = 0.D+0
      RN_ID = 0.D+0
!
!---  Incident solar radiation zero, skip calculations  ---
!
      IF( RN/EPSL.LT.EPSL ) GOTO 100
!
!---  Equation of time, hr  ---
!
      REALX = REAL(J_DAY-81)
      B_EQT = 2.D+0*GPI*REALX/364.D+0
      T_EQT = 1.645D-1*SIN(2.D+0*B_EQT) - 1.255D-1*COS(B_EQT) -
     &  2.5D-2*SIN(B_EQT)
!
!---  Time of solar noon, hr
!     ATMC(3) longitude, radians
!     ATMC(4) latitude, radians
!     ATMC(5) meridian, radians  ---
!
      T_SN = 1.2D+1 + (ATMC(3)-ATMC(5))*1.2D+1/GPI - T_EQT
!
!---  Time of day, hr  ---
!
      T_DAY = MOD(TM+ATMST,8.64D+4)/3.6D+3
!
!---  Solar time angle, radians  ---
!
      S_TANG = (T_DAY-T_SN)*GPI/1.2D+1
!
!---  Solar declination, radians  ---
!
      REALX = REAL(J_DAY-172)
      S_DEC = 4.101D-1*COS(2.D+0*GPI*REALX/3.65D+2)
!
!---  Sun angle, radians  ---
!
      S_ANG = ASIN( SIN(S_DEC)*SIN(ATMC(4)) +
     &  COS(S_DEC)*COS(ATMC(4))*COS(S_TANG) )
!
!---  Sunset time angle, radians  ---
!
      SS_ANG = ACOS( -TAN(ATMC(3))*TAN(S_DEC) )
!
!---  Sun zenith angle, radians  ---
!
      S_ZA = ACOS( SIN(S_DEC)*SIN(ATMC(4)) +
     &  COS(S_DEC)*COS(ATMC(4))*COS(S_TANG) )
!
!---  Sun angle greater than sunrise or sunset time angle,
!     incident radiation is diffuse  ---
!
      IF( ABS(S_ANG).GT.ABS(SS_ANG) ) THEN
        RN_ID = RN
        GOTO 100
      ENDIF
!
!---  Relative earth-sun distance  ---
!
      REALX = REAL(J_DAY)
      D_ESR = 1.D+0 + 3.3D-2*COS(GPI*2.D+0*REALX/3.65D+2)
!
!---  Extraterrestrial irradiance for a one hour period
!     around the current time  ---
!
      S_ANG1 = MAX( S_ANG - (GPI/4.8D+1),-SS_ANG )
      S_ANG2 = MIN( S_ANG + (GPI/4.8D+1),SS_ANG )
      R_SA = (1.2D+1/GPI)*G_SC*D_ESR*(COS(ATMC(3))*COS(S_DEC)*
     &  (SIN(S_ANG1)-SIN(S_ANG2)) +
     &  (S_ANG1-S_ANG2)*SIN(ATMC(3))*SIN(S_DEC))
!
!---  Ratio of diffuse to incident radiation as a function of
!     the ratio of incident to extraterresterial radiation  ---
!
      TAU_T = RN/R_SA
      CRX = 8.47D-1 - 1.61D+0*SIN(S_ANG) + 1.04D+0*(SIN(S_ANG)**2)
      CKX = (1.47D+0-CRX)/1.66D+0
      IF( TAU_T.LE.2.2D-1 ) THEN
        RN_ID = RN
      ELSEIF( TAU_T.GT.2.2D-1 .AND. TAU_T.LE.3.5D-1 ) THEN
        RN_ID = MIN( (1.D+0 - 6.4D+0*((TAU_T-2.2D-1)**2)),1.D+0 )*RN
        RN_DB = MAX( RN-RN_ID,0.D+0 )
      ELSEIF( TAU_T.GT.3.5D-1 .AND. TAU_T.LE.CKX ) THEN
        RN_ID = MIN( (1.47D+0 - 1.66D+0*TAU_T),1.D+0 )*RN
        RN_DB = MAX( RN-RN_ID,0.D+0 )
      ELSEIF( TAU_T.GT.CKX ) THEN
        RN_ID = MIN( CRX,1.D+0 )*RN
        RN_DB = MAX( RN-RN_ID,0.D+0 )
      ENDIF
      TAU_D = RN_ID/R_SA
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SOLR_C group  ---
!
      RETURN
      END


!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TRAPZD( FUNC,A,B,S,N,INDX )
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
!     This routine computes the N'th stage of refinement of an
!     extended trapezoid rule.
!
!     Press, W.H., B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling.
!     1986.  Numerical Recipes, The Art of Scientific Computing.
!     Cambridge University Press, Cambridge.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, February, 1999.
!     Last Modified by Mark White, Battelle, February 19, 1999.
!     roots.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      EXTERNAL FUNC
      SAVE IT
!
!----------------------Executable Lines--------------------------------!
!
      IF( N.EQ.1 ) THEN
        S = 5.D-1*(B-A)*(FUNC(A,INDX)+FUNC(B,INDX))
!
!---  IT is the number of points to be added on the next call  ---
!
        IT = 1
      ELSE
        REALX = REAL(IT)
        TNM = REALX
!
!---  Spacing of the points to be added.  ---
!
        DEL = (B-A)/TNM
        X = A + 5.D-1*DEL
        SUM = 0.D+0
        DO 100 J = 1,IT
          SUM = SUM + FUNC(X,INDX)
          X = X + DEL
  100   CONTINUE
!
!---  Replace S by its refined value  ---
!
        S = 5.D-1*(S+(B-A)*SUM/TNM)
        IT = 2*IT
      ENDIF
!
!---  End of TRAPZD group
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION VRUGT( ZX,IP )
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
!     Vrugt one-dimensional root water uptake model function.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 19 August 2002.
!     Last Modified by Mark White, PNNL, 19 August 2002.
!     roots.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PLT_ATM
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VRUGT'
!
!---  Set maximum rooting depth, null root depth, and root depth
!     fitting parameter  ---
!
      ZM = PARMS_P(1,IP)
      ZSX = PARMS_P(4,IP)
      PZ = PARMS_P(7,IP)
      VRUGT = (1.D+0 - ZX/(ZM+1.D-20))*EXP(-(PZ/ZM)*ABS(ZSX-ZX))
!
!---  End of VRUGT group
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

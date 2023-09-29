!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_IC1( BTGLX,PGX,PLX,PVAX,SGTX,SLX,TX,TMSX,XMLAX,
     &  YLSX,ICBRNX,ICAIRX,ICSX,N )
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
!      Copyright Battelle Memorial Institute, 1996
!          All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Flash calculation for initial conditions for
!     phase conditions #1 or #3
!
!     SL + SGT = 1.0
!     Aqueuous saturated w/ or w/o trapped gas
!
!     Input variables:
!
!     TX - temperature, C
!     PGX - pressure, abs, Pa
!     PLX - pressure, abs, Pa
!     PVAX - aqueous air relative saturation, or
!     aqueous air mass fraction
!     TMSX - aqueous salt relative saturation, or
!     aqueous salt mass fraction
!
!     Input indices
!
!     ICBRNX - initial condition index for brine concentration
!     ICAIRX - initial condition index for air concentration
!
!     Output variables:
!
!     XMLAX - aqueous mass fraction of air
!     YLSX - total salt brine mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 27 March 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 GX(2,3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_IC1'
!
!---  Absolute temperature and pressure  ---
!
      TKX = TX + TABS
!
!---  Gas pressure and aqueous saturation  ---
!
      IF( ICSX.EQ.1 ) THEN
!
!---    Minimum aqueous saturation computed by CAP_GT  ---
!
        CALL SCF_GL( BTGLX,PGX )
        ASLMINX = -1.D+0
        IZN = IZ(N)
        CALL CAP_GT( ASLMINX,BTGLX,PCX,SLX,SGTX,IZN )
        PLX = PGX - PCX
!
!---    Top of single-variable Newton-Raphson loop  ---
!
        NC = 0
   10   CONTINUE
        NC = NC + 1
        IF( NC.GT.32 ) THEN
          INDX = 17
          RLMSG = PGX
          CHMSG = 'FLH_IC1 Initial Conditions w/  : ' //
     &      'Gas Pressure + Aqueous Saturation : ' // 
     &      'Unconverged Aqu. Pressure : ' // 
     &      'Gas Pressure = '
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Set increments  ---
!
        DPLY = -1.D-1
!
!---    Loop over increments  ---
!
        DO 12 M = 1,2
          PLY = PLX
          IF( M.EQ.2 ) PLY = PLX + DPLY
          CALL SCF_GL( BTGLX,PLY )
          INDX = 1
          IZN = IZ(N)
          CALL CAP_GT( ASLMINX,BTGLX,PCX,SLX,SGTX,IZN )
          GX(1,M) = PGX - (PLY+PCX)
   12   CONTINUE
        FX = -GX(1,1)
        DFX = (GX(1,2)-GX(1,1))/DPLY
        DPLX = FX/(DFX+SMALL)
        PLX = MIN( PLX+DPLX,PGX )
        IF( ABS(DPLX).GT.1.D-1 ) GOTO 10
        CALL SCF_GL( BTGLX,PLX )
!
!---  Aqueous pressure and aqueous saturation  ---
!
      ELSEIF( ICSX.EQ.2 ) THEN
!
!---    Minimum aqueous saturation computed by CAP_GT  ---
!
        CALL SCF_GL( BTGLX,PLX )
        ASLMINX = -1.D+0
        IZN = IZ(N)
        CALL CAP_GT( ASLMINX,BTGLX,PCX,SLX,SGTX,IZN )
        PGX = PLX + PCX
!
!---  Gas and aqueous pressure  ---
!
      ELSE
        CALL SCF_GL( BTGLX,PLX )
      ENDIF
!
!---  Total pressure  ---
!
      PX = PLX
!
!---  Saturated water vapor pressure  ---
!
      INDX = 0
      CALL REGION_4( TX,PSWX,INDX )
      IF( PLX.LT.PSWX ) THEN
        INDX = 17
        RLMSG = PLX
        CHMSG = 'FLH_IC1 Initial Conditions: Aqu. Pressure ' //
     &    '< Water Sat. Pressure : Aqu. Pressure, Pa = '
         CALL WRMSGS( INDX )
      ENDIF
      PWX = PLX
!
!---  Aqueous-salt concentration  ---
!
      IF( ICBRNX.EQ.1 ) THEN
!
!---    Aqueous-air concentration  ---
!
        IF( ICAIRX.EQ.1 ) THEN
          RHOLSX = TMSX
          RHOLAX = PVAX
          CALL P_IAPWS( TX,PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
          CALL SOL_BRNS( TX,PWX,XLSMX )
!
!---      Guess aqueous-salt mass fractions  ---
!
          XLSX = RHOLSX/(RHOLWX + RHOLSX)
!
!---      Top of single-variable Newton-Raphson loop  ---
!
          NC = 0
   20     CONTINUE
          NC = NC + 1
          IF( NC.GT.32 ) THEN
            INDX = 17
            RLMSG = RHOLSX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Conc. + Aqu. Salt Conc. : ' //
     &        'Unconverged : Aqu. Salt Conc., kg/m^3 = '
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Set increments  ---
!
          DXLSY = SIGN( 1.D-6,(5.D-1*XLSMX-XLSX) )
!
!---      Loop over increments  ---
!
          DO 22 M = 1,2
            XLSY = XLSX
            IF( M.EQ.2 ) XLSY = XLSX + DXLSY
            CALL DENS_B( XLSY,RHOBX,RHOLWX,TX )
            GX(1,M) = RHOLSX - RHOBX*XLSY
   22     CONTINUE
          FX = -GX(1,1)
          DFX = (GX(1,2)-GX(1,1))/DXLSY
          DXLSX = FX/(DFX+SMALL)
          XLSX = MAX( MIN( XLSX+DXLSX,XLSMX ),0.D+0 )
          IF( ABS(DXLSX).GT.1.D-6 ) GOTO 20
          CALL SP_B( TX,XLSX,PSBX )
          CALL DENS_B( XLSX,RHOBX,RHOLWX,TX )
          XLAX = RHOLAX/RHOBX
          XLWX = MAX( 1.D+0-XLSX-XLAX,0.D+0 )
          WTMLX = XLWX/WTMW + XLAX/WTMA + XLSX/WTMS
          XMLAX = (XLAX/WTMA)/WTMLX       
          IF( PX.LT.PSBX ) THEN
            INDX = 17
            RLMSG = PSBX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Conc. + Aqu. Salt Conc. : Pressure ' //
     &        '< Saturated Pressure : Saturated Pressure, Pa = '
            CALL WRMSGS( INDX )
          ENDIF
!
!---    Aqueous-air relative saturation  ---
!
        ELSEIF( ICAIRX.EQ.2 ) THEN
          RHOLSX = TMSX
          PHILAX = PVAX
          IF( PHILAX.LT.0.D+0 .OR. PHILAX.GT.1.D+0 ) THEN
            INDX = 17
            RLMSG = PHILAX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Rel. Sat. + Aqu. Salt Conc. : ' //
     &        '1 < Aqu. Air Rel. Sat. < 0 : ' //
     *        'Aqu. Air Rel. Sat. = '
            CALL WRMSGS( INDX )
          ENDIF
          CALL P_IAPWS( TX,PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
          CALL SOL_BRNS( TX,PWX,XLSMX )
!
!---      Guess aqueous-salt mass fractions  ---
!
          XLSX = RHOLSX/(RHOLWX + RHOLSX)
!
!---      Top of single-variable Newton-Raphson loop  ---
!
          NC = 0
   30     CONTINUE
          NC = NC + 1
          IF( NC.GT.32 ) THEN
            INDX = 17
            RLMSG = RHOLSX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Rel. Sat. + Aqu. Salt Conc. : ' //
     &        'Unconverged : Aqu. Salt Conc., kg/m^3 = '
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Set increments  ---
!
          DXLSY = SIGN( 1.D-6,(5.D-1*XLSMX-XLSX) )
!
!---      Loop over increments  ---
!
          DO 32 M = 1,2
            XLSY = XLSX
            IF( M.EQ.2 ) XLSY = XLSX + DXLSY
            CALL DENS_B( XLSY,RHOBX,RHOLWX,TX )
            GX(1,M) = RHOLSX - RHOBX*XLSY
   32     CONTINUE
          FX = -GX(1,1)
          DFX = (GX(1,2)-GX(1,1))/DXLSY
          DXLSX = FX/(DFX+SMALL)
          XLSX = MAX( MIN( XLSX+DXLSX,XLSMX ),0.D+0 )
          IF( ABS(DXLSX).GT.1.D-6 ) GOTO 30
          CALL SP_B( TX,XLSX,PSBX )
          CALL DENS_B( XLSX,RHOBX,RHOLWX,TX )
          XMLAX = PHILAX*PWX/HCAW    
          IF( PX.LT.PSBX ) THEN
            INDX = 17
            RLMSG = PSBX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Rel. Sat. + Aqu. Salt Conc. : Pressure ' //
     &        '< Saturated Pressure : Saturated Pressure, Pa = '
            CALL WRMSGS( INDX )
          ENDIF
!
!---    Aqueous-air mass fraction  ---
!
        ELSEIF( ICAIRX.EQ.3 .OR. ICAIRX.EQ.0 ) THEN
          RHOLSX = TMSX
          XLAX = PVAX
          CALL P_IAPWS( TX,PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
          CALL SOL_BRNS( TX,PWX,XLSMX )
!
!---      Guess aqueous-salt mass fractions  ---
!
          XLSX = RHOLSX/(RHOLWX + RHOLSX)
!
!---      Top of single-variable Newton-Raphson loop  ---
!
          NC = 0
   40     CONTINUE
          NC = NC + 1
          IF( NC.GT.32 ) THEN
            INDX = 17
            RLMSG = RHOLSX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Mass Frac. + Aqu. Salt Conc. : ' //
     &        'Unconverged : Aqu. Salt Conc., kg/m^3 = '
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Set increments  ---
!
          DXLSY = SIGN( 1.D-6,(5.D-1*XLSMX-XLSX) )
!
!---      Loop over increments  ---
!
          DO 42 M = 1,2
            XLSY = XLSX
            IF( M.EQ.2 ) XLSY = XLSX + DXLSY
            CALL DENS_B( XLSY,RHOBX,RHOLWX,TX )
            GX(1,M) = RHOLSX - RHOBX*XLSY
   42     CONTINUE
          FX = -GX(1,1)
          DFX = (GX(1,2)-GX(1,1))/DXLSY
          DXLSX = FX/(DFX+SMALL)
          XLSX = MAX( MIN( XLSX+DXLSX,XLSMX ),0.D+0 )
          IF( ABS(DXLSX).GT.1.D-6 ) GOTO 40
          CALL SP_B( TX,XLSX,PSBX )
          CALL DENS_B( XLSX,RHOBX,RHOLWX,TX )
          XLWX = MAX( 1.D+0-XLSX-XLAX,0.D+0 )
          WTMLX = XLWX/WTMW + XLAX/WTMA + XLSX/WTMS
          XMLAX = (XLAX/WTMA)/WTMLX       
          IF( PX.LT.PSBX ) THEN
            INDX = 17
            RLMSG = PSBX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Mass Frac. + Aqu. Salt Conc. : ' //
     &        'Pressure < Saturated Press. : Saturated Pressure, Pa = '
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
!
!---  Aqueous-salt relative saturation  ---
!
      ELSEIF( ICBRNX.EQ.2 ) THEN
!
!---    Aqueous-air concentration  ---
!
        IF( ICAIRX.EQ.1 ) THEN
          PHILSX = TMSX
          RHOLAX = PVAX
          IF( PHILSX.LT.0.D+0 .OR. PHILSX.GT.1.D+0 ) THEN
            INDX = 17
            RLMSG = PHILSX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Conc. + Aqu. Salt Rel. Sat. :' //
     &        '1 < Aqu. Salt Rel. Sat. < 0 : Aqu. Salt Rel. Sat. = '
            CALL WRMSGS( INDX )
          ENDIF
          CALL P_IAPWS( TX,PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
          CALL SOL_BRNS( TX,PWX,XLSMX )
          XLSX = PHILSX*XLSMX
          CALL SP_B( TX,XLSX,PSBX )
          CALL DENS_B( XLSX,RHOBX,RHOLWX,TX )
          XLAX = RHOLAX/RHOBX
          XLWX = MAX( 1.D+0-XLSX-XLAX,0.D+0 )
          WTMLX = XLWX/WTMW + XLAX/WTMA + XLSX/WTMS
          XMLAX = (XLAX/WTMA)/WTMLX       
          IF( PX.LT.PSBX ) THEN
            INDX = 17
            RLMSG = PSBX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Conc. + Aqu. Salt Rel. Sat. : Pressure ' //
     &        '< Saturated Pressure : Saturated Pressure, Pa = '
            CALL WRMSGS( INDX )
          ENDIF
!
!---    Aqueous-air relative saturation  ---
!
        ELSEIF( ICAIRX.EQ.2 ) THEN
          PHILSX = TMSX
          IF( PHILSX.LT.0.D+0 .OR. PHILSX.GT.1.D+0 ) THEN
            INDX = 17
            RLMSG = PHILSX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Rel. Sat. + Aqu. Salt Rel. Sat. :' //
     &        '1 < Aqu. Salt Rel. Sat. < 0 : Aqu. Salt. Rel. Sat. = '
            CALL WRMSGS( INDX )
          ENDIF
          PHILAX = PVAX
          IF( PHILAX.LT.0.D+0 .OR. PHILAX.GT.1.D+0 ) THEN
            INDX = 17
            RLMSG = PHILAX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Rel. Sat. + Aqu. Salt Rel. Sat. :' //
     &        '1 < Aqu. Air Rel. Sat. < 0 : Aqu. Air Rel. Sat. = '
            CALL WRMSGS( INDX )
          ENDIF
          CALL P_IAPWS( TX,PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
          CALL SOL_BRNS( TX,PWX,XLSMX )
          XLSX = PHILSX*XLSMX
          CALL SP_B( TX,XLSX,PSBX )
          CALL DENS_B( XLSX,RHOBX,RHOLWX,TX )
          XMLAX = PHILAX*PWX/HCAW    
          IF( PX.LT.PSBX ) THEN
            INDX = 17
            RLMSG = PSBX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Rel. Sat. + Aqu. Salt Rel. Sat. : Pressure ' //
     &        '< Saturated Pressure : Saturated Pressure, Pa = '
            CALL WRMSGS( INDX )
          ENDIF
!
!---    Aqueous-air mass fraction  ---
!
        ELSEIF( ICAIRX.EQ.3 .OR. ICAIRX.EQ.0 ) THEN
          PHILSX = TMSX
          IF( PHILSX.LT.0.D+0 .OR. PHILSX.GT.1.D+0 ) THEN
            INDX = 17
            RLMSG = PHILSX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Mass Frac. + Aqu. Salt Rel. Sat. :' //
     &        '1 < Aqu. Salt Rel. Sat. < 0 : Aqu. Salt Rel. Sat. = '
            CALL WRMSGS( INDX )
          ENDIF
          XLAX = PVAX
          CALL P_IAPWS( TX,PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
          CALL SOL_BRNS( TX,PWX,XLSMX )
          XLSX = PHILSX*XLSMX
          CALL SP_B( TX,XLSX,PSBX )
          CALL DENS_B( XLSX,RHOBX,RHOLWX,TX )
          XLWX = MAX( 1.D+0-XLSX-XLAX,0.D+0 )
          WTMLX = XLWX/WTMW + XLAX/WTMA + XLSX/WTMS
          XMLAX = (XLAX/WTMA)/WTMLX       
          IF( PX.LT.PSBX ) THEN
            INDX = 17
            RLMSG = PSBX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Mass Frac. + Aqu. Salt Rel. Sat. : Pressure ' //
     &        '< Saturated Pressure : Saturated Pressure, Pa = '
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
!
!---  Aqueous-salt mass fraction  ---
!
      ELSEIF( ICBRNX.EQ.3 .OR. ICBRNX.EQ.0 ) THEN
!
!---    Aqueous-air concentration  ---
!
        IF( ICAIRX.EQ.1 ) THEN
          YLSX = TMSX
          RHOLAX = PVAX
          CALL P_IAPWS( TX,PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
          CALL SOL_BRNS( TX,PWX,XLSMX )
          IF( YLSX.GT.XLSMX ) THEN
            INDX = 17
            RLMSG = XLSMX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Conc. + Aqu. Salt Mass Frac. : ' //
     &        'Aqu. Salt Mass Frac. > Solubility Limit : ' //
     &        'Aqu. Salt Mass Frac. Solubility Limit = '
            CALL WRMSGS( INDX )
          ENDIF
          XLSX = YLSX
          CALL SP_B( TX,XLSX,PSBX )
          CALL DENS_B( XLSX,RHOBX,RHOLWX,TX )
          XLAX = RHOLAX/RHOBX
          XLWX = MAX( 1.D+0-XLSX-XLAX,0.D+0 )
          WTMLX = XLWX/WTMW + XLAX/WTMA + XLSX/WTMS
          XMLAX = (XLAX/WTMA)/WTMLX       
          IF( PX.LT.PSBX ) THEN
            INDX = 17
            RLMSG = PSBX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Conc. + Aqu. Salt Rel. Sat. : Pressure ' //
     &        '< Saturated Pressure : Saturated Pressure, Pa = '
            CALL WRMSGS( INDX )
          ENDIF
!
!---    Aqueous-air relative saturation  ---
!
        ELSEIF( ICAIRX.EQ.2 ) THEN
          YLSX = TMSX
          PHILAX = PVAX
          IF( PHILAX.LT.0.D+0 .OR. PHILAX.GT.1.D+0 ) THEN
            INDX = 17
            RLMSG = PHILAX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Rel. Sat. + Aqu. Salt Mass Frac. : ' //
     &        '1 < Aqu. Air Rel. Sat. < 0 : Aqu. Air Rel. Sat. = '
            CALL WRMSGS( INDX )
          ENDIF
          CALL P_IAPWS( TX,PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
          CALL SOL_BRNS( TX,PWX,XLSMX )
          IF( YLSX.GT.XLSMX ) THEN
            INDX = 17
            RLMSG = XLSMX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Conc. + Aqu. Salt Mass Frac. : ' //
     &        'Aqu. Salt Mass Frac. > Solubility Limit : ' //
     &        'Aqu. Salt Mass Frac. Solubility Limit = '
            CALL WRMSGS( INDX )
          ENDIF
          XLSX = YLSX
          CALL SP_B( TX,XLSX,PSBX )
          CALL DENS_B( XLSX,RHOBX,RHOLWX,TX )
          XMLAX = PHILAX*PWX/HCAW    
          IF( PX.LT.PSBX ) THEN
            INDX = 17
            RLMSG = PSBX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Rel. Sat. + Aqu. Salt Rel. Sat. : Pressure ' //
     &        '< Saturated Pressure : Saturated Pressure, Pa = '
            CALL WRMSGS( INDX )
          ENDIF
!
!---    Aqueous-air mass fraction  ---
!
        ELSEIF( ICAIRX.EQ.3 .OR. ICAIRX.EQ.0 ) THEN
          YLSX = TMSX
          XLAX = PVAX
          CALL P_IAPWS( TX,PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
          CALL SOL_BRNS( TX,PWX,XLSMX )
          IF( YLSX.GT.XLSMX ) THEN
            INDX = 17
            RLMSG = XLSMX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Conc. + Aqu. Salt Mass Frac. : ' //
     &        'Aqu. Salt Mass Frac. > Solubility Limit : ' //
     &        'Aqu. Salt Mass Frac. Solubility Limit = '
            CALL WRMSGS( INDX )
          ENDIF
          XLSX = YLSX
          CALL SP_B( TX,XLSX,PSBX )
          CALL DENS_B( XLSX,RHOBX,RHOLWX,TX )
          XLWX = MAX( 1.D+0-XLSX-XLAX,0.D+0 )
          WTMLX = XLWX/WTMW + XLAX/WTMA + XLSX/WTMS
          XMLAX = (XLAX/WTMA)/WTMLX       
          IF( PX.LT.PSBX ) THEN
            INDX = 17
            RLMSG = PSBX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Mass Frac. + Aqu. Salt Rel. Sat. : ' //
     &        '< Saturated Pressure : Saturated Pressure, Pa = '
            CALL WRMSGS( INDX )
          ENDIF
       ENDIF
!
!---  Aqueous-salt molality (mol salt/kg water)  ---
!
      ELSEIF( ICBRNX.EQ.4 ) THEN
!
!---    Aqueous-air concentration  ---
!
        IF( ICAIRX.EQ.1 ) THEN
          YLSX = TMSX
          RHOLAX = PVAX
          CALL P_IAPWS( TX,PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
          CALL SOL_BRNS( TX,PWX,XLSMX )
          IF( YLSX.GT.EPSL ) THEN
            YLSX = 1.D+0/(1.D+0 + 1.D+0/(YLSX*WTMS))
          ELSE
            YLSX = 0.D+0
          ENDIF
          IF( YLSX.GT.XLSMX ) THEN
            INDX = 17
            RLMSG = XLSMX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Conc. + Aqu. Salt Molality : ' //
     &        'Aqu. Salt Mass Frac. > Solubility Limit : ' //
     &        'Aqu. Salt Mass Frac. Solubility Limit = '
            CALL WRMSGS( INDX )
          ENDIF
          XLSX = YLSX
          CALL SP_B( TX,XLSX,PSBX )
          CALL DENS_B( XLSX,RHOBX,RHOLWX,TX )
          XLAX = RHOLAX/RHOBX
          XLWX = MAX( 1.D+0-XLSX-XLAX,0.D+0 )
          WTMLX = XLWX/WTMW + XLAX/WTMA + XLSX/WTMS
          XMLAX = (XLAX/WTMA)/WTMLX       
          IF( PX.LT.PSBX ) THEN
            INDX = 17
            RLMSG = PSBX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Conc. + Aqu. Salt Molality : Pressure ' //
     &        '< Saturated Pressure : Saturated Pressure, Pa = '
            CALL WRMSGS( INDX )
          ENDIF
!
!---    Aqueous-air relative saturation  ---
!
        ELSEIF( ICAIRX.EQ.2 ) THEN
          YLSX = TMSX
          PHILAX = PVAX
          IF( PHILAX.LT.0.D+0 .OR. PHILAX.GT.1.D+0 ) THEN
            INDX = 17
            RLMSG = PHILAX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Rel. Sat. + Aqu. Salt Molality :' //
     &        '1 < Aqu. Air Rel. Sat. < 0 : Aqu. Air Rel. Sat. = '
            CALL WRMSGS( INDX )
          ENDIF
          CALL P_IAPWS( TX,PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
          CALL SOL_BRNS( TX,PWX,XLSMX )
          IF( YLSX.GT.EPSL ) THEN
            YLSX = 1.D+0/(1.D+0 + 1.D+0/(YLSX*WTMS))
          ELSE
            YLSX = 0.D+0
          ENDIF
          IF( YLSX.GT.XLSMX ) THEN
            INDX = 17
            RLMSG = XLSMX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Conc. + Aqu. Salt Molality : ' //
     &        'Aqu. Salt Mass Frac. > Solubility Limit : ' //
     &        'Aqu. Salt Mass Frac. Solubility Limit = '
            CALL WRMSGS( INDX )
          ENDIF
          XLSX = YLSX
          CALL SP_B( TX,XLSX,PSBX )
          CALL DENS_B( XLSX,RHOBX,RHOLWX,TX )
          XMLAX = PHILAX*PWX/HCAW    
          IF( PX.LT.PSBX ) THEN
            INDX = 17
            RLMSG = PSBX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Rel. Sat. + Aqu. Salt Molality : Pressure ' //
     &        '< Saturated Pressure : Saturated Pressure, Pa = '
            CALL WRMSGS( INDX )
          ENDIF
!
!---    Aqueous-air mass fraction  ---
!
        ELSEIF( ICAIRX.EQ.3 .OR. ICAIRX.EQ.0 ) THEN
          XLAX = PVAX
          YLSX = TMSX
          CALL P_IAPWS( TX,PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
          CALL SOL_BRNS( TX,PWX,XLSMX )
          IF( YLSX.GT.EPSL ) THEN
            YLSX = 1.D+0/(1.D+0 + 1.D+0/(YLSX*WTMS))
          ELSE
            YLSX = 0.D+0
          ENDIF
          IF( YLSX.GT.XLSMX ) THEN
            INDX = 17
            RLMSG = XLSMX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Conc. + Aqu. Salt Molality :' //
     &        'Aqu. Salt Mass Frac. > Solubility Limit : ' //
     &        'Aqu. Salt Mass Frac. Solubility Limit = '
            CALL WRMSGS( INDX )
          ENDIF
          XLSX = YLSX
          CALL SP_B( TX,XLSX,PSBX )
          CALL DENS_B( XLSX,RHOBX,RHOLWX,TX )
          XLWX = MAX( 1.D+0-XLSX-XLAX,0.D+0 )
          WTMLX = XLWX/WTMW + XLAX/WTMA + XLSX/WTMS
          XMLAX = (XLAX/WTMA)/WTMLX       
          IF( PX.LT.PSBX ) THEN
            INDX = 17
            RLMSG = PSBX
            CHMSG = 'FLH_IC1 Initial Conditions w/ ' // 
     &        'Aqu. Air Mass Frac. + Aqu. Salt Molality : Pressure ' //
     &        '< Saturated Pressure : Saturated Pressure, Pa = '
            CALL WRMSGS( INDX )
          ENDIF
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
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
      SUBROUTINE FLH_IC2( BTGLX,PGX,PLX,SGTX,SLX,TX,TMSX,XMLAX,YLSX,
     &  ICBRNX,ICSX,N )
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
!      Copyright Battelle Memorial Institute, 1996
!          All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Flash calculation for initial conditions for
!     phase conditions #2
!
!     0.0 < SL + SGT < 1.0
!     Aqueuous saturated
!
!     Input variables:
!
!     TX - temperature, C
!     PGX - pressure, abs, Pa
!     PLX - pressure, abs, Pa
!     TMSX - aqueous salt relative saturation, or
!     aqueous salt mass fraction
!
!     Input indices
!
!     ICBRNX - initial condition index for brine concentration
!
!     Output variables:
!
!     XMLAX - aqueous mass fraction of air
!     YLSX - total salt brine mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 30 March 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 GX(2,3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_IC2'
      IZN = IZ(N)
!
!---  Absolute temperature and pressure  ---
!
      TKX = TX + TABS
!
!---  Saturated water vapor pressure  ---
!
      INDX = 0
      CALL REGION_4( TX,PSWX,INDX )
      CALL P_IAPWS( TX,PSWX,RHOX,RHOLWX,HX,HLWX,UX,ULWX )
      CALL SCF_GL( BTGLX,PSWX )
!
!---  Gas pressure and aqueous saturation  ---
!
      IF( ICSX.EQ.1 ) THEN
!
!---    Minimum aqueous saturation computed by CAP_GT  ---
!
        ASLMINX = -1.D+0
        IZN = IZ(N)
        CALL CAP_GT( ASLMINX,BTGLX,PCX,SLX,SGTX,IZN )
        PLX = PGX - PCX
!
!---    Top of single-variable Newton-Raphson loop  ---
!
        NC = 0
   10   CONTINUE
        NC = NC + 1
        IF( NC.GT.32 ) THEN
          INDX = 17
          RLMSG = SLX
          CHMSG = 'FLH_IC2 Initial Conditions w/  : ' //
     &      'Gas Pressure + Aqueous Saturation : ' // 
     &      'Unconverged Aqu. Pressure : ' // 
     &      'Aqu. Saturation = '
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Set increments  ---
!
        DPLY = -1.D-1
!
!---    Loop over increments  ---
!
        DO 12 M = 1,2
          PLY = PLX
          IF( M.EQ.2 ) PLY = PLX + DPLY
          PCX = PGX - PLY
          CALL VPL_B( TX,PSWX,PCX,RHOLWX,PVWX )
          CALL SCF_GL( BTGLX,PVWX )
          INDX = 1
          CALL SP_GT( ASLFX,ASLMX,ASLX,ASLMINX,ELSX,ESGTX,ESGTMX,
     &      PGX,PLY,BTGLX,SGX,SGTX,SLY,SLFX,
     &      SLMX,SLRX,INDX,IZN )
          GX(1,M) = SLX - SLY
   12   CONTINUE
        FX = -GX(1,1)
        DFX = (GX(1,2)-GX(1,1))/DPLY
        DPLX = FX/(DFX+SMALL)
        PLX = MIN( PLX+DPLX,PGX )
        IF( ABS(DPLX).GT.1.D-1 ) GOTO 10
        PCX = PGX - PLX
        CALL VPL_B( TX,PSWX,PCX,RHOLWX,PVWX )
        CALL SCF_GL( BTGLX,PVWX )
        INDX = 1
        CALL SP_GT( ASLFX,ASLMX,ASLX,ASLMINX,ELSX,ESGTX,ESGTMX,
     &    PGX,PLX,BTGLX,SGX,SGTX,SLX,SLFX,
     &    SLMX,SLRX,INDX,IZN )
!
!---  Aqueous pressure and aqueous saturation  ---
!
      ELSEIF( ICSX.EQ.2 ) THEN
!
!---    Minimum aqueous saturation computed by CAP_GT  ---
!
        ASLMINX = -1.D+0
        IZN = IZ(N)
        CALL CAP_GT( ASLMINX,BTGLX,PCX,SLX,SGTX,IZN )
        PGX = PLX + PCX
!
!---    Top of single-variable Newton-Raphson loop  ---
!
        NC = 0
   20   CONTINUE
        NC = NC + 1
        IF( NC.GT.32 ) THEN
          INDX = 17
          RLMSG = SLX
          CHMSG = 'FLH_IC2 Initial Conditions w/ : ' //
     &      'Gas Pressure + Aqueous Saturation : ' // 
     &      'Unconverged Aqu. Pressure : ' // 
     &      'Aqu. Saturation = '
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Set increments  ---
!
        DPGY = 1.D-1
!
!---    Loop over increments  ---
!
        DO 22 M = 1,2
          PGY = PGX
          IF( M.EQ.2 ) PGY = PGX + DPGY
          PCX = PGY - PLX
          CALL VPL_B( TX,PSWX,PCX,RHOLWX,PVWX )
          CALL SCF_GL( BTGLX,PVWX )
          INDX = 1
          CALL SP_GT( ASLFX,ASLMX,ASLX,ASLMINX,ELSX,ESGTX,ESGTMX,
     &      PGY,PLX,BTGLX,SGX,SGTX,SLY,SLFX,
     &      SLMX,SLRX,INDX,IZN )
          GX(1,M) = SLX - SLY
   22   CONTINUE
        FX = -GX(1,1)
        DFX = (GX(1,2)-GX(1,1))/DPGY
        DPGX = FX/(DFX+SMALL)
        PGX = MAX( PGX+DPGX,PLX )
        IF( ABS(DPGX).GT.1.D-1 ) GOTO 20
        PCX = PGX - PLX
        CALL VPL_B( TX,PSWX,PCX,RHOLWX,PVWX )
        CALL SCF_GL( BTGLX,PVWX )
        INDX = 1
        CALL SP_GT( ASLFX,ASLMX,ASLX,ASLMINX,ELSX,ESGTX,ESGTMX,
     &    PGX,PLX,BTGLX,SGX,SGTX,SLX,SLFX,
     &    SLMX,SLRX,INDX,IZN )
!
!---  Gas and aqueous pressure  ---
!
      ELSE
        PCX = PGX - PLX
        CALL VPL_B( TX,PSWX,PCX,RHOLWX,PVWX )
        CALL SCF_GL( BTGLX,PVWX )
        INDX = 1
        CALL SP_GT( ASLFX,ASLMX,ASLX,ASLMINX,ELSX,ESGTX,ESGTMX,
     &    PGX,PLX,BTGLX,SGX,SGTX,SLX,SLFX,
     &    SLMX,SLRX,INDX,IZN )
      ENDIF
      PX = PGX
      PWX = PVWX
!
!---  Aqueous-salt concentration  ---
!
      IF( ICBRNX.EQ.1 ) THEN
        RHOLSX = TMSX
        CALL P_IAPWS( TX,PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
        CALL SOL_BRNS( TX,PWX,XLSMX )
!
!---    Guess aqueous-salt mass fractions  ---
!
        XLSX = RHOLSX/(RHOLWX + RHOLSX)
!
!---    Top of single-variable Newton-Raphson loop  ---
!
        NC = 0
   30   CONTINUE
        NC = NC + 1
        IF( NC.GT.32 ) THEN
          INDX = 17
          RLMSG = RHOLSX
          CHMSG = 'FLH_IC2 Initial Conditions w/ : ' // 
     &      'Aqu. Salt Conc. : Unconverged : ' //
     &      'Aqu. Salt Conc., kg/m^3 = '
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Set increments  ---
!
        DXLSY = SIGN( 1.D-6,(5.D-1*XLSMX-XLSX) )
!
!---    Loop over increments  ---
!
        DO 32 M = 1,2
          XLSY = XLSX
          IF( M.EQ.2 ) XLSY = XLSX + DXLSY
          CALL DENS_B( XLSY,RHOBX,RHOLWX,TX )
          GX(1,M) = RHOLSX - RHOBX*XLSY
   32   CONTINUE
        FX = -GX(1,1)
        DFX = (GX(1,2)-GX(1,1))/DXLSY
        DXLSX = FX/(DFX+SMALL)
        XLSX = MAX( MIN( XLSX+DXLSX,XLSMX ),0.D+0 )
        IF( ABS(DXLSX).GT.1.D-6 ) GOTO 30
        CALL SP_B( TX,XLSX,PSBX )
        CALL DENS_B( XLSX,RHOBX,RHOLWX,TX )
        PCX = PGX - PLX
        CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX )
        IF( PGX.LT.PVBX ) THEN
          INDX = 17
          RLMSG = PVBX
          CHMSG = 'FLH_IC2 Initial Conditions w/ ' // 
     &      'Aqu. Salt Conc. : Gas Pressure < Reduced ' //
     &      'Saturated Pressure : Reduced Saturated Pressure, Pa = '
          CALL WRMSGS( INDX )
        ENDIF
        PVAX = MAX( PGX-PVBX,0.D+0 )
        XMLAX = PVAX/HCAW
!
!---  Aqueous-salt relative saturation  ---
!
      ELSEIF( ICBRNX.EQ.2 ) THEN
        PHILSX = TMSX
        IF( PHILSX.LT.0.D+0 .OR. PHILSX.GT.1.D+0 ) THEN
          INDX = 17
          RLMSG = PHILSX
          CHMSG = 'FLH_IC2 Initial Conditions w/ ' // 
     &      '1 < Aqu. Salt Rel. Sat. < 0 :' //
     &      'Aqu. Salt Rel. Sat. = '
          CALL WRMSGS( INDX )
        ENDIF
        CALL P_IAPWS( TX,PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
        CALL SOL_BRNS( TX,PWX,XLSMX )
        XLSX = PHILSX*XLSMX
        CALL SP_B( TX,XLSX,PSBX )
        CALL DENS_B( XLSX,RHOBX,RHOLWX,TX )
        PCX = PGX - PLX
        CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX )
        IF( PGX.LT.PVBX ) THEN
          INDX = 17
          RLMSG = PVBX
          CHMSG = 'FLH_IC2 Initial Conditions w/ ' // 
     &      'Aqu. Salt Rel. Sat. : Gas Pressure < Reduced ' //
     &      'Saturated Pressure : Reduced Saturated Pressure, Pa = '
          CALL WRMSGS( INDX )
        ENDIF
        PVAX = MAX( PGX-PVBX,0.D+0 )
        XMLAX = PVAX/HCAW
!
!---  Aqueous-salt mass fraction  ---
!
      ELSEIF( ICBRNX.EQ.3 .OR. ICBRNX.EQ.0 ) THEN
        YLSX = TMSX
        CALL P_IAPWS( TX,PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
        CALL SOL_BRNS( TX,PWX,XLSMX )
        IF( YLSX.GT.XLSMX ) THEN
          INDX = 17
          RLMSG = XLSMX
          CHMSG = 'FLH_IC2 Initial Conditions w/ ' // 
     &      'Aqu. Salt Mass Frac. : Aqu. Salt Mass Frac. > ' //
     &      'Solubility Limit : Aqu. Salt Mass Frac. Sol. Limit = '
          CALL WRMSGS( INDX )
        ENDIF
        XLSX = YLSX
        CALL SP_B( TX,XLSX,PSBX )
        CALL DENS_B( XLSX,RHOBX,RHOLWX,TX )
        PCX = PGX - PLX
        CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX )
        IF( PGX.LT.PVBX ) THEN
          INDX = 17
          RLMSG = PVBX
          CHMSG = 'FLH_IC2 Initial Conditions w/ ' // 
     &      'Aqu. Salt Rel. Sat. : Gas Pressure < Reduced ' //
     &      'Saturated Pressure : Reduced Saturated Pressure, Pa = '
          CALL WRMSGS( INDX )
        ENDIF
        PVAX = MAX( PGX-PVBX,0.D+0 )
        XMLAX = PVAX/HCAW
!
!---  Aqueous-salt molality (mol salt/kg water)  ---
!
      ELSEIF( ICBRNX.EQ.4 ) THEN
        YLSX = TMSX
        CALL P_IAPWS( TX,PWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
        CALL SOL_BRNS( TX,PWX,XLSMX )
        IF( YLSX.GT.EPSL ) THEN
          YLSX = 1.D+0/(1.D+0 + 1.D+0/(YLSX*WTMS))
        ELSE
          YLSX = 0.D+0
        ENDIF
        IF( YLSX.GT.XLSMX ) THEN
          INDX = 17
          RLMSG = XLSMX
          CHMSG = 'FLH_IC2 Initial Conditions w/ ' // 
     &      'w/ Aqu. Salt Molality : Aqu. Salt Mass Frac. > ' //
     &      'Solubility Limit : Aqu. Salt Mass Frac. Sol. Limit = '
          CALL WRMSGS( INDX )
        ENDIF
        XLSX = YLSX
        CALL SP_B( TX,XLSX,PSBX )
        CALL DENS_B( XLSX,RHOBX,RHOLWX,TX )
        PCX = PGX - PLX
        CALL VPL_B( TX,PSBX,PCX,RHOBX,PVBX )
        IF( PGX.LT.PVBX ) THEN
          INDX = 17
          RLMSG = PVBX
          CHMSG = 'FLH_IC2 Initial Conditions w/ ' // 
     &      'Aqu. Salt Molality : Gas Pressure < Reduced ' //
     &      'Saturated Pressure : Reduced Saturated Pressure, Pa = '
          CALL WRMSGS( INDX )
        ENDIF
        PVAX = MAX( PGX-PVBX,0.D+0 )
        XMLAX = PVAX/HCAW
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
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
      SUBROUTINE FLH_IC4( BTGLX,PGX,PLX,PVAX,PVWX,SGTX,SLX,TX,TMSX,
     &  XMGAX,YLSX,ICAIRX,ICSX,N )
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
!      Copyright Battelle Memorial Institute, 1996
!          All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Flash calculation for initial conditions for
!     phase conditions #4
!
!     SL = 0.0
!     Aqueuous unsaturated
!
!     Input variables:
!
!     TX - temperature, C
!     PGX - pressure, abs, Pa
!     PVWX - gas-water relative saturation, or
!     gas-water mass fraction, or gas-water concentration
!     TMSX - aqueous salt relative saturation
!
!     Input indices
!
!     ICBRNX - initial condition index for brine concentration
!     ICAIRX - initial condition index for water concentration
!
!     Output variables:
!
!     SLX - aqueous saturation 
!     PVAX - air partial pressure
!     PLX - pressure, abs, Pa
!     YLSX - total salt brine mass fraction
!     TMSX - total salt mass
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 27 March 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 GX(2,3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_IC4'
!
!---  Absolute temperature and pressure  ---
!
      TKX = TX + TABS
      PX = PGX
      IF( TMSX.GT.EPSL ) THEN
        INDX = 17
        RLMSG = TMSX
        CHMSG = 'FLH_IC4 Initial Conditions w/ ' // 
     &    'Non-Zero Salt Concentration Declared : ' //
     &    'Salt Concentration = '
        CALL WRMSGS( INDX )
      ENDIF
      YLSX = 0.D+0
      TMSX = 0.D+0
!
!---  Gas-water concentration  ---
!
      IF( ICAIRX.EQ.1 ) THEN
        RHOGWX = PVWX
!
!---    Guess water vapor partial pressure  ---
!
        CALL AIRGSD( TX,PGX,RHOGAX )
        RHOMGWX = (RHOGWX/WTMW)/(RHOGWX/WTMW + RHOGAX/WTMA)
        PVWX = PGX*RHOMGWX
!
!---    Top of single-variable Newton-Raphson loop  ---
!
        NC = 0
   10   CONTINUE
        NC = NC + 1
        IF( NC.GT.32 ) THEN
          INDX = 17
          RLMSG = RHOGWX
          CHMSG = 'FLH_IC4 Initial Conditions w/ ' // 
     &      'Gas-Water Conc. : Unconverged : Gas-Water Conc., kg/m^3 = '
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Set increments  ---
!
        DPVWY = 1.D-6*PGX
!
!---    Loop over increments  ---
!
        DO M = 1,2
          PVWY = PVWX
          IF( M.EQ.2 ) PVWY = PVWX + DPVWY
          CALL P_IAPWS( TX,PVWY,RHOGWY,RHOLWX,HGWX,HLWX,UGWX,ULWX )
          PGAX = PGX - PVWY
          CALL AIRGSD( TX,PGAX,RHOGAX )
          RHOGX = RHOGAX + RHOGWY
          GX(1,M) = RHOGWX - RHOGWY
        ENDDO
        FX = -GX(1,1)
        DFX = (GX(1,2)-GX(1,1))/DPVWY
        DPWX = FX/(DFX+SMALL)
        PVWX = MAX( MIN( PVWX+DPWX,PGX ),0.D+0 )
        IF( ABS(DPWX).GT.1.D-3 ) GOTO 10
        CALL P_IAPWS( TX,PVWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
        PGAX = PGX - PVWX
        CALL AIRGSD( TX,PGAX,RHOGAX )
        RHOGX = RHOGAX + RHOGWX
        PVAX = MAX( PGX-PVWX,0.D+0 )
!
!---  Gas-water relative saturation  ---
!
      ELSEIF( ICAIRX.EQ.2 ) THEN
        PHIGWX = PVWX
        IF( PHIGWX.LT.0.D+0 .OR. PHIGWX.GT.1.D+0 ) THEN
          INDX = 17
          RLMSG = PHIGWX
          CHMSG = 'FLH_IC4 Initial Conditions w/ ' // 
     &      '1 < Gas-Water Rel. Sat. < 0 : Gas-Water Rel. Sat. = '
          CALL WRMSGS( INDX )
        ENDIF
        PVWX = PHIGWX*PGX
        PVAX = PGX - PVWX
        CALL P_IAPWS( TX,PVWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
        CALL AIRGSD( TX,PVAX,RHOGAX )
        RHOGX = RHOGAX + RHOGWX
        XGAX = RHOGAX/RHOGX
        XGWX = RHOGWX/RHOGX
        WTMGX = XGWX/WTMW + XGAX/WTMA
        XMGAX = XGAX/WTMA/WTMGX
!
!---  Gas-water mass fraction  ---
!
      ELSEIF( ICAIRX.EQ.3 .OR. ICAIRX.EQ.0 ) THEN
        XGWX = PVWX
        IF( XGWX.LT.0.D+0 .OR. XGWX.GT.1.D+0 ) THEN
          INDX = 17
          RLMSG = XGWX
          CHMSG = 'FLH_IC4 Initial Conditions w/ ' // 
     &      '1 < Gas-Water Mass Frac. < 0 : Gas-Water Mass Frac. = '
          CALL WRMSGS( INDX )
        ENDIF
        XGAX = 1.D+0-XGWX
        WTMGX = XGWX/WTMW + XGAX/WTMA
        XMGWX = XGWX/WTMW/WTMGX
        XMGAX = XGAX/WTMA/WTMGX
        PVWX = PGX*XMGWX
        PVAX = PGX*XMGAX
        CALL P_IAPWS( TX,PVWX,RHOGWX,RHOLWX,HGWX,HLWX,UGWX,ULWX )
        CALL AIRGSD( TX,PVAX,RHOGAX )
        RHOGX = RHOGAX + RHOGWX
      ENDIF
!
!---  Capillary pressure scaling  ---
!
      CALL SCF_GL( BTGLX,PVWX )
!
!---  Gas pressure and aqueous saturation  ---
!
      IF( ICSX.EQ.1 ) THEN
!
!---    Minimum aqueous saturation computed by CAP_GT  ---
!
        ASLMINX = -1.D+0
        IZN = IZ(N)
        CALL CAP_GT( ASLMINX,BTGLX,PCX,SLX,SGTX,IZN )
        PLX = PGX - PCX
!
!---  Aqueous pressure and aqueous saturation  ---
!
      ELSEIF( ICSX.EQ.2 ) THEN
        CALL SCF_GL( BTGLX,PVWX )
!
!---    Minimum aqueous saturation computed by CAP_GT  ---
!
        ASLMINX = -1.D+0
        IZN = IZ(N)
        CALL CAP_GT( ASLMINX,BTGLX,PCX,SLX,SGTX,IZN )
        PGX = PLX + PCX
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_IC4 group  ---
!
      RETURN
      END



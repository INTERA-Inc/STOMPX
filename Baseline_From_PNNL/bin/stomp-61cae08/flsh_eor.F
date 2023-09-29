!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BO_IC1( PX,RHOLX,TMSX,TX,XLSX,XLWX,YLSX,IC8X )
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
!     SL = 1.0
!     SG+SN = 0.0
!     Aqueuous saturated
!
!     Input variables:
!
!     TX - temperature, C
!     PX - pressure, abs, Pa
!     TMSX - aqueous salt relative saturation, or
!            aqueous salt mass fraction
!
!     Output variables:
!
!     RHOLX - aqueous density, kg/m^3
!     XLSX - aqueous mass fraction of salt
!     XLWX - aqueous mass fraction of water
!     YLSX - total salt brine mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 19 May 2014.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 GX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BO_IC1'
!
!---  Absolute temperature  ---
!
      TKX = TX + TABS
!
!---  Aqueous salt mass fraction specified, solve for XBSX  ---
!
      IF( IC8X.EQ.3 ) THEN
        XLSX = TMSX
!
!---    Guest salt brine mass fraction  ---
!
        XBSX = XLSX
!
!---    Zero aqueous gas mass fraction  ---
!
        XLAX = 0.D+0
        PVAX = 0.D+0
!
!---    Salt brine mass fraction increment  ---
!
        DXY = 1.D-7
!
!---    Top of Newton-Raphson loop for salt brine mass fraction
!       and CO2 partial pressure  ---
!
  100   CONTINUE
        DO 110 M = 1,2
          XBSY = XBSX
          IF( M.EQ.2 ) XBSY = XBSX + DXY
          CALL SP_B( TX,XBSY,PSBX )
          CALL EQUIL( TX,PX,PVAX,XBSY,XLAX,XLSY,XLWY )
          GX(M) = XLSX - XLSY
  110   CONTINUE
!
!---    Load solution vector and Jacobian matrix  ---
!
        FX = GX(1)
        DFX = (GX(2)-GX(1))/DXY
!
!---    Update primary unknowns  ---
!
        DXBSX = -FX/DFX
        XBSX = XBSX + DXBSX
        XBSX = MIN( MAX( XBSX,0.D+0 ),1.D+0 )
!
!---    Convergence check  ---
!
        IF( ABS(DXBSX).GT.1.D-8 ) GOTO 100
!
!---  Aqueous salt relative saturation specified, solve for XBSX  ---
!
      ELSEIF( IC8X.EQ.2 ) THEN
        WLSX = TMSX
        CALL SOL_LS( TX,XBSMX )
!
!---    Guest salt brine mass fraction  ---
!
        XBSX = WLSX*XBSMX
!
!---    Zero aqueous gas mass fraction  ---
!
        XLAX = 0.D+0
        PVAX = 0.D+0
!
!---    Salt brine mass fraction increment  ---
!
        DXY = 1.D-7
!
!---    Top of Newton-Raphson loop for salt brine mass fraction  ---
!
  300   CONTINUE
        DO 310 M = 1,2
          XBSY = XBSX
          IF( M.EQ.2 ) XBSY = XBSX + DXY
          CALL SP_B( TX,XBSY,PSBX )
          CALL EQUIL( TX,PX,PVAX,XBSY,XLAX,XLSY,XLWY )
          CALL EQUIL( TX,PX,PX,XBSMX,XLAX,XLSZ,XLWZ )
          WLSY = XLSY/XLSZ
          GX(M) = WLSX - WLSY
  310   CONTINUE
!
!---    Load solution vector and Jacobian matrix  ---
!
        FX = GX(1)
        DFX = (GX(2)-GX(1))/DXY
!
!---    Update primary unknowns  ---
!
        DXBSX = -FX/DFX
        XBSX = XBSX + DXBSX
        XBSX = MIN( MAX( XBSX,0.D+0 ),1.D+0 )
!
!---    Convergence check  ---
!
        IF( ABS(DXBSX).GT.1.D-8 ) GOTO 300
      ENDIF
!
!---  Aqueous mass fractions  ---
!
      CALL SP_B( TX,XBSX,PSBX )
      CALL EQUIL( TX,PX,PVAX,XBSX,XLAX,XLSX,XLWX )
      WTMLX = XLSX/WTMS + XLWX/WTMW
      XMLSX = XLSX/WTMS/WTMLX
      XMLWX = XLWX/WTMW/WTMLX
      XMLSX = XMLSX/(XMLSX+XMLWX)
      XMLWX = XMLWX/(XMLSX+XMLWX)
      WTMLX = XMLSX*WTMS + XMLWX*WTMW
      XLSX = XMLSX*WTMS/WTMLX
      XLWX = XMLWX*WTMW/WTMLX
!
!---  Brine density based on CO2-free brine  ---
!
      CALL DENS_B( TX,PX,XBSX,RHOBX )
      YLSX = XBSX
!
!---  Aqueous density based on dissoved CO2 and salt  ---
!
      CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BO_IC1 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BO_IC2( POSMX,PX,RHOGX,RHOLX,RHONX,SLX,
     &  TMSX,TX,XGCX,XGWX,XLSX,XLWX,XNCX,YLSX,IC8X,IZN )
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
!---  Initial condition state #2
!
!     SG = 0.0
!     SL > 0.0
!     SN > 0.0
!     Aqueous and nonaqueous liquid saturated
!
!     Declared variables:
!
!     TX - temperature, C
!     PX - pressure, absolute, Pa
!     SLX - aqueous saturation
!     POSMX - gas-oil ratio
!     TMSX - aqueous salt relative saturation, or
!            aqueous salt mass fraction
!
!---  Initial condition state #3
!
!     SG > 0.0
!     SL > 0.0
!     SN > 0.0
!     Three-phase
!
!     Declared variables:
!
!     TX - temperature, C
!     PX - pressure, Pa
!     SLX - aqueous saturation
!     SGX - gas saturation
!     TMSX - aqueous salt relative saturation, or
!            aqueous salt mass fraction
!
!     Output variables:
!
!     PGX - gas pressure, guage Pa
!     PLX - aqueous pressure, guage, Pa
!     PNX - nonaqueous liquid pressure, guage Pa
!     RHOGX - gas density, kg/m^3
!     RHOLX - aqueous density, kg/m^3
!     RHONX - nonaqueous density, kg/m^3
!     XGCX - gas mass fraction of petroleum components
!     XGWX - gas mass fraction of water
!     XLSX - aqueous mass fraction of salt
!     XLWX - aqueous mass fraction of water
!     XNCX - nonaqueous liquid mass fraction of petroleum components
!     YLSX - brine mass fraction of total salt
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 July 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GLB_PAR
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XGCX(LNGC),XNCX(LNGC)
      REAL*8 GX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BO_IC2'
!
!---  Absolute temperature  ---
!
      TKX = TX + TABS
      PVAX = 0.D+0
      XLAX = 0.D+0
!
!---  Aqueous salt mass fraction, solve for XBSX  ---
!
      IF( IC8X.EQ.3 ) THEN
        XLSX = TMSX
!
!---    Guest salt brine mass fraction  ---
!
        XBSX = XLSX
        DXY = 1.D-7
!
!---    Top of Newton-Raphson loop for salt brine mass fraction  ---
!
  100   CONTINUE
        DO 110 M = 1,2
          XBSY = XBSX
          IF( M.EQ.2 ) XBSY = XBSX + DXY
          CALL SP_B( TX,XBSY,PSBX )
!
!---      Vapor-pressure lowering  ---
!
          IF( IZN.NE.0 ) THEN
            CALL VPL( PSBX,PVWX,SLX,TX,IZN )
!
!---      No vapor-pressure lowering, used in source and well model  ---
!
          ELSE
            PVWX = PSBX
          ENDIF
          PVX = PX - PVWX
          CALL SP_B( TX,XBSY,PSBX )
          PTX = PSBX + PVAX
          CALL EQUIL( TX,PTX,PVAX,XBSY,XLAX,XLSY,XLWY )
          GX(M) = XLSX - XLSY
  110   CONTINUE
!
!---    Load solution vector and Jacobian matrix  ---
!
        DFX = (GX(2)-GX(1))/DXY
        FX = -GX(1)
!
!---    Update primary unknowns  ---
!
        DXBSX = FX/DFX
        XBSX = XBSX + DXBSX
        XBSX = MIN( MAX( XBSX,0.D+0 ),1.D+0 )
!
!---    Convergence check  ---
!
        IF( ABS(DXBSX).GT.1.D-8 ) GOTO 100
!
!---  Aqueous salt relative saturation, solve for XBSX  ---
!
      ELSEIF( IC8X.EQ.2 ) THEN
        WLSX = TMSX
        CALL SOL_LS( TX,XBSMX )
!
!---    Guest salt brine mass fraction  ---
!
        XBSX = WLSX*XBSMX
        DXY = 1.D-7
!
!---    Top of Newton-Raphson loop for salt brine mass fraction  ---
!
  200   CONTINUE
        DO 210 M = 1,2
          XBSY = XBSX
          IF( M.EQ.2 ) XBSY = XBSX + DXY
          CALL SP_B( TX,XBSY,PSBX )
!
!---      Vapor-pressure lowering  ---
!
          IF( IZN.NE.0 ) THEN
            CALL VPL( PSBX,PVWX,SLX,TX,IZN )
!
!---      No vapor-pressure lowering, used in source and well model  ---
!
          ELSE
            PVWX = PSBX
          ENDIF
          PVX = PX - PVWX
          PTX = PSBX + PVAX
          CALL EQUIL( TX,PTX,PVAX,XBSY,XLAX,XLSY,XLWY )
          CALL EQUIL( TX,PTX,PVAX,XBSMX,XLAX,XLSZ,XLWZ )
          WLSY = XLSY/XLSZ
          GX(M) = WLSX - WLSY
  210   CONTINUE
!
!---    Load solution vector and Jacobian matrix  ---
!
        DFX = (GX(2)-GX(1))/DXY
        FX = -GX(1)
!
!---    Update primary unknowns  ---
!
        DXBSX = FX/DFX
        XBSX = XBSX + DXBSX
        XBSX = MIN( MAX( XBSX,0.D+0 ),1.D+0 )
!
!---    Convergence check  ---
!
        IF( ABS(DXBSX).GT.1.D-8 ) GOTO 200
      ENDIF
!
!---  Aqueous mass fractions  ---
!
      CALL SP_B( TX,XBSX,PSBX )
!
!---  Vapor-pressure lowering  ---
!
      IF( IZN.NE.0 ) THEN
        CALL VPL( PSBX,PVWX,SLX,TX,IZN )
!
!---  No vapor-pressure lowering, used in source and well model  ---
!
      ELSE
        PVWX = PSBX
      ENDIF
      PVX = PX - PVWX
!
!---  Gas density and composition  ---
!
      CALL BO_EQ_GAS( PX,PVWX,RHOGX,TX,XGCX,XGWX )
!
!---  Nonaqueous-liquid density and composition  ---
!
      CALL BO_EQ_OIL( RHONX,PX,POSMX,TX,XNCX )
      CALL SP_B( TX,XBSX,PSBX )
      PTX = PSBX + PVAX
      CALL EQUIL( TX,PTX,PVAX,XBSX,XLAX,XLSX,XLWX )
      WTMLX = XLSX/WTMS + XLWX/WTMW
      XMLSX = XLSX/WTMS/WTMLX
      XMLWX = XLWX/WTMW/WTMLX
      XMLSX = XMLSX/(XMLSX+XMLWX)
      XMLWX = XMLWX/(XMLSX+XMLWX)
      WTMLX = XMLSX*WTMS + XMLWX*WTMW
      XLSX = XMLSX*WTMS/WTMLX
      XLWX = XMLWX*WTMW/WTMLX
!
!---  Brine density based on CO2-free brine  ---
!
      CALL DENS_B( TX,PX,XBSX,RHOBX )
      YLSX = XBSX
!
!---  Aqueous density based on dissoved CO2 and salt  ---
!
      CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BO_IC2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BO_IC3( POSMX,PVWX,PX,RHOGX,RHONX,TX,
     &  XGCX,XGWX,XNCX,YLSX,IZN )
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
!---  Flash calculation for initial condition state #4
!
!     SG = 0.0
!     SL = 0.0
!     SN = 1.0
!     Nonaqueous liquid saturated
!
!     Declared variables:
!
!     TX - temperature, C
!     PSOX - pressure, Pa
!     POSMX - gas-oil ratio
!     PVWX - relative saturation of water vapor
!
!---  Flash calculation for initial condition state #5
!
!     SG > 0.0
!     SL = 0.0
!     SN > 0.0
!     Nonaqueous liquid and gas saturated
!
!     Declared variables:
!
!     TX - temperature, C
!     PSOX - pressure, Pa
!     SGX - gas saturation
!     PVWX - relative saturation of water vapor
!
!     Output variables:
!
!     RHOGX - gas density, kmol/m^3
!     RHONX - nonaqueous-liquid density, kmol/m^3
!     SGX - gas saturation
!     SNX - nonaqueous-liquid saturation
!     XGCX - gas mass fraction of petroleum components
!     XNCX - nonaqueous-liquid mass fraction of petroleum components
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 July 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GLB_PAR
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XGCX(LNGC),XNCX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BO_IC3'
!
!---  Absolute temperature  ---
!
      TKX = TX + TABS
!
!---  Maximum solubility of salt in pure water  ---
!
      CALL SOL_LS( TX,XBSMX )
!
!---  Salt mass fraction in brine  ---
!
      XBSX = MIN( YLSX,XBSMX )
!
!---  Water vapor pressure in equilibrium with brine at saturated
!     conditions  ---
!
      CALL SP_B( TX,XBSX,PSBX )
!
!---  Water vapor pressure adjusted for aqueous saturation
!     in geologic media  ---
!
      SLX = 0.D+0
!
!---  Water vapor relative saturation  ---
!
      PSBX = PSBX*PVWX
!
!---  Vapor-pressure lowering  ---
!
      IF( IZN.NE.0 ) THEN
        CALL VPL( PSBX,PVWX,SLX,TX,IZN )
!
!---  No vapor-pressure lowering, used in source and well model  ---
!
      ELSE
        PVWX = PSBX
      ENDIF
!
!---  Total nonaqueous vapor pressure  ---
!
      PVX = PX - PVWX
!
!---  Gas density and composition  ---
!
      CALL BO_EQ_GAS( PX,PVWX,RHOGX,TX,XGCX,XGWX )
!
!---  Nonaqueous-liquid density and composition  ---
!
      CALL BO_EQ_OIL( RHONX,PX,POSMX,TX,XNCX )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BO_IC3 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BO_PC1( PORDX,POSMX,PVWX,PX,RHOGX,RHOLX,RHONX,SGX,SLX,
     &  SNX,TMSX,TX,VOLX,XGCX,XGWX,XLSX,XLWX,XNCX,YLSX )
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
!     Flash calculation for phase condition #1 for the black-oil option
!
!     SL = 1.0
!     SG+SN = 0.0
!     Aqueuous saturated
!
!     energy - temperature, C - T
!     water mass - aqueous (system) pressure, Pa gauge - PL,PSO
!     gas mass - total gas mass - TMC(1,2,N)
!     oil mass - total oil mass - TMC(2,2,N)
!     salt mass - total salt brine mass fraction - YLS
!
!     Output variables:
!
!     RHOGX - gas density, kg/m^3
!     RHOLX - aqueous density, kg/m^3
!     RHONX - nonaqueous-liquid density, kg/m^3
!     SGX - gas saturation
!     SNX - nonaqueous-liquid saturation
!     XGCX - gas mass fraction of petroleum components
!     XGWX - gas mass fraction of water
!     XLSX - aqueous mass fraction of salt
!     XLWX - aqueous mass fraction of water
!     XNCX - nonaqueous-liquid mass fraction of petroleum components
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 21 May 2014
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GLB_PAR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XGCX(LNGC),XNCX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BO_PC1'
!
!---  No dissolved gas in aqueous  ---
!
      PVAX = 0.D+0
      XLAX = 0.D+0
!
!---  Maximum solubility of salt in pure water  ---
!
      CALL SOL_LS( TX,XBSMX )
!
!---  Salt mass fraction in brine  ---
!
      XBSX = MIN( YLSX,XBSMX )
!
!---  Water vapor pressure in equilibrium with brine at saturated
!     conditions  ---
!
      CALL SP_B( TX,XBSX,PSBX )
!
!---  Water vapor pressure adjusted for aqueous saturation
!     in geologic media  ---
!
      PVWX = PSBX
!
!---  Brine density  ---
!
      CALL DENS_B( TX,PX,XBSX,RHOBX )
!
!---  Aqueous density  ---
!
      CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
!
!---  Aqueous component mass fractions  ---
!
      XLSX = XBSX
      XLWX = 1.D+0-XLSX
!
!---  Total mass of salt  ---
!
      TMSX = YLSX*RHOBX*PORDX*VOLX
!
!---  Nonaqueous mole fractions, densities, and saturations  ---
!
      CALL BO_GOR( PX,POSMX,TX )
      CALL BO_EQ_GAS( PX,PVBX,RHOGX,TX,XGCX,XGWX )
      CALL BO_EQ_OIL( RHONX,PX,POSMX,TX,XNCX )
      SLX = 1.D+0
      SGX = 0.D+0
      SNX = 0.D+0
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BO_PC1 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BO_PC2( PORDX,POSMX,PVWX,PX,RHOGX,RHOLX,RHONX,SGX,SLX,
     &  SNX,TMCX,TMSX,TX,VOLX,XGCX,XGWX,XLSX,XLWX,XNCX,YLSX,IZN )
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
!     Flash calculation for phase condition #2 for the black-oil option
!
!     SG = 0.0
!     SL > 0.0
!     SN > 0.0
!     Aqueous and nonaqueous liquid saturated
!
!     energy - temperature, C - T
!     water mass - aqueous saturation, SL
!     gas mass - gas to oil ratio, POSM
!     oil mass - nonaqueous (system) pressure, Pa gauge - PN,PSO
!     salt mass - total salt brine mass fraction - YLS
!
!     Output variables:
!
!     RHOGX - gas density, kg/m^3
!     RHOLX - aqueous density, kg/m^3
!     RHONX - nonaqueous-liquid density, kg/m^3
!     SGX - gas saturation
!     SNX - nonaqueous-liquid saturation
!     XGCX - gas mass fraction of petroleum components
!     XGWX - gas mass fraction of water
!     XLSX - aqueous mass fraction of salt
!     XLWX - aqueous mass fraction of water
!     XNCX - nonaqueous-liquid mass fraction of petroleum components
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 20 May 2014
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GLB_PAR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 TMCX(LNGC),XGCX(LNGC),XNCX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BO_PC2'
!
!---  No dissolved gas in aqueous  ---
!
      PVAX = 0.D+0
      XLAX = 0.D+0
!
!---  Maximum solubility of salt in pure water  ---
!
      CALL SOL_LS( TX,XBSMX )
!
!---  Salt mass fraction in brine  ---
!
      XBSX = MIN( YLSX,XBSMX )
!
!---  Water vapor pressure in equilibrium with brine at saturated
!     conditions  ---
!
      CALL SP_B( TX,XBSX,PSBX )
!
!---  Water vapor pressure adjusted for aqueous saturation
!     in geologic media  ---
!
      CALL VPL( PSBX,PVWX,SLX,TX,IZN )
!
!---  Nonaqueous mole fractions, densities, and saturations  ---
!
      CALL BO_EQ_GAS( PX,PVWX,RHOGX,TX,XGCX,XGWX )
      CALL BO_EQ_OIL( RHONX,PX,POSMX,TX,XNCX )
      SGX = 0.D+0
      SNX = 1.D+0-SLX
!
!---  Brine density  ---
!
      CALL DENS_B( TX,PX,XBSX,RHOBX )
!
!---  Aqueous density  ---
!
      CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
!
!---  Aqueous component mass fractions  ---
!
      XLSX = XBSX
      XLWX = 1.D+0-XLSX
!
!---  Total mass of salt  ---
!
      TMSX = YLSX*RHOBX*PORDX*VOLX
!
!---  Total mass of gas  ---
!
      TMCX(1) = PORDX*VOLX*SNX*XNCX(1)*RHONX
!
!---  Total mass oil  ---
!
      TMCX(2) = PORDX*VOLX*SNX*XNCX(2)*RHONX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BO_PC2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BO_PC3( PORDX,POSMX,PVWX,PX,RHOGX,RHOLX,RHONX,SGX,SLX,
     &  SNX,TMCX,TMSX,TX,VOLX,XGCX,XGWX,XLSX,XLWX,XNCX,YLSX,IZN )
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
!     Flash calculation for phase condition #3 for the black-oil option
!
!     SG > 0.0
!     SL > 0.0
!     SN > 0.0
!     Three-phase
!
!     energy - temperature, C - T
!     water mass - aqueous saturation, SL
!     gas mass - gas (system) pressure, Pa gauge - PG,PSO
!     oil mass - nonaqueous-liquid saturation, SN
!     salt mass - total salt brine mass fraction - YLS
!
!     Output variables:
!
!     POSMX - gas to oil ratio, m^3/m^3
!     RHOGX - gas density, kg/m^3
!     RHOLX - aqueous density, kg/m^3
!     RHONX - nonaqueous-liquid density, kg/m^3
!     SGX - gas saturation
!     SNX - nonaqueous-liquid saturation
!     XGCX - gas mass fraction of petroleum components
!     XGWX - gas mass fraction of water
!     XLSX - aqueous mass fraction of salt
!     XLWX - aqueous mass fraction of water
!     XNCX - nonaqueous-liquid mass fraction of petroleum components
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 20 May 2014
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GLB_PAR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 TMCX(LNGC),XGCX(LNGC),XNCX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BO_PC3'
!
!---  No dissolved gas in aqueous  ---
!
      PVAX = 0.D+0
      XLAX = 0.D+0
!
!---  Maximum solubility of salt in pure water  ---
!
      CALL SOL_LS( TX,XBSMX )
!
!---  Salt mass fraction in brine  ---
!
      XBSX = MIN( YLSX,XBSMX )
!
!---  Water vapor pressure in equilibrium with brine at saturated
!     conditions  ---
!
      CALL SP_B( TX,XBSX,PSBX )
!
!---  Water vapor pressure adjusted for aqueous saturation
!     in geologic media  ---
!
      CALL VPL( PSBX,PVWX,SLX,TX,IZN )
!
!---  Nonaqueous mole fractions, densities, and saturations  ---
!
      CALL BO_GOR( PX,POSMX,TX )
      CALL BO_EQ_GAS( PX,PVWX,RHOGX,TX,XGCX,XGWX )
      CALL BO_EQ_OIL( RHONX,PX,POSMX,TX,XNCX )
      SNX = 1.D+0-SGX-SLX
!
!---  Brine density  ---
!
      CALL DENS_B( TX,PX,XBSX,RHOBX )
!
!---  Aqueous density  ---
!
      CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
!
!---  Aqueous component mass fractions  ---
!
      XLSX = XBSX
      XLWX = 1.D+0-XLSX
!
!---  Total mass of salt  ---
!
      TMSX = YLSX*RHOBX*PORDX*VOLX
!
!---  Total mass of gas  ---
!
      TMCX(1) = PORDX*VOLX*(SGX*XGCX(1)*RHOGX + SNX*XNCX(1)*RHONX)
!
!---  Total mass oil  ---
!
      TMCX(2) = PORDX*VOLX*SNX*XNCX(2)*RHONX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BO_PC3 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BO_PC4( PORDX,POSMX,PVWX,PX,RHOGX,RHOLX,RHONX,SGX,SLX,
     &  SNX,TMCX,TMSX,TX,VOLX,XGCX,XGWX,XLSX,XLWX,XNCX,YLSX,IZN )
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
!     Flash calculation for phase condition #4 for the black-oil option
!
!     SG = 0.0
!     SL = 0.0
!     SN = 1.0
!     Nonaqueous liquid saturated
!
!     energy - temperature, C - T
!     water mass - aqueous saturation, SL
!     gas mass - gas to oil ratio, POSM
!     oil mass - nonaqueous (system) pressure, Pa gauge - PN,PSO
!     salt mass - total salt mass - TMS(2,N)
!
!     Output variables:
!
!     RHOGX - gas density, kg/m^3
!     RHOLX - aqueous density, kg/m^3
!     RHONX - nonaqueous-liquid density, kg/m^3
!     SGX - gas saturation
!     SNX - nonaqueous-liquid saturation
!     XGCX - gas mass fraction of petroleum components
!     XGWX - gas mass fraction of water
!     XLSX - aqueous mass fraction of salt
!     XLWX - aqueous mass fraction of water
!     XNCX - nonaqueous-liquid mass fraction of petroleum components
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 20 May 2014
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GLB_PAR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 TMCX(LNGC),XGCX(LNGC),XNCX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BO_PC4'
!
!---  No dissolved gas in aqueous  ---
!
      PVAX = 0.D+0
      XLAX = 0.D+0
!
!---  Maximum solubility of salt in pure water  ---
!
      CALL SOL_LS( TX,XBSMX )
!
!---  Salt mass fraction in brine  ---
!
      XBSX = MIN( YLSX,XBSMX )
!
!---  Water vapor pressure in equilibrium with brine at saturated
!     conditions  ---
!
      CALL SP_B( TX,XBSX,PSBX )
!
!---  Water vapor pressure adjusted for aqueous saturation
!     in geologic media  ---
!
      CALL VPL( PSBX,PVWX,SLX,TX,IZN )
!
!---  Nonaqueous mole fractions, densities, and saturations  ---
!
      CALL BO_EQ_GAS( PX,PVWX,RHOGX,TX,XGCX,XGWX )
      CALL BO_EQ_OIL( RHONX,PX,POSMX,TX,XNCX )
      SGX = 0.D+0
      SNX = 1.D+0-SLX
!
!---  Brine density  ---
!
      CALL DENS_B( TX,PX,XBSX,RHOBX )
!
!---  Aqueous density  ---
!
      CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
!
!---  Aqueous component mass fractions  ---
!
      XLSX = XBSX
      XLWX = 1.D+0-XLSX
!
!---  Total mass of salt  ---
!
      TMSX = YLSX*RHOBX*PORDX*VOLX
!
!---  Total mass of gas  ---
!
      TMCX(1) = PORDX*VOLX*SNX*XNCX(1)*RHONX
!
!---  Total mass oil  ---
!
      TMCX(2) = PORDX*VOLX*SNX*XNCX(2)*RHONX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BO_PC4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BO_PC5( PORDX,POSMX,PVWX,PX,RHOGX,RHOLX,RHONX,SGX,SLX,
     &  SNX,TMCX,TMSX,TX,VOLX,XGCX,XGWX,XLSX,XLWX,XNCX,YLSX,IZN )
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
!     Flash calculation for phase condition #5 for the black-oil option
!
!     SG > 0.0
!     SL = 0.0
!     SN > 0.0
!     Nonaqueous liquid and gas saturated
!
!     energy - temperature, C - T
!     water mass - aqueous saturation, SL
!     gas mass - gas (system) pressure, Pa gauge - PG,PSO
!     oil mass - nonaqueous-liquid saturation, SN
!     salt mass - total salt mass - TMS(2,N)
!
!     Output variables:
!
!     POSMX - gas to oil ratio, m^3/m^3
!     RHOGX - gas density, kg/m^3
!     RHOLX - aqueous density, kg/m^3
!     RHONX - nonaqueous-liquid density, kg/m^3
!     SGX - gas saturation
!     SNX - nonaqueous-liquid saturation
!     XGCX - gas mass fraction of petroleum components
!     XGWX - gas mass fraction of water
!     XLSX - aqueous mass fraction of salt
!     XLWX - aqueous mass fraction of water
!     XNCX - nonaqueous-liquid mass fraction of petroleum components
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 20 May 2014
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GLB_PAR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 TMCX(LNGC),XGCX(LNGC),XNCX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BO_PC5'
!
!---  No dissolved gas in aqueous  ---
!
      PVAX = 0.D+0
      XLAX = 0.D+0
!
!---  Maximum solubility of salt in pure water  ---
!
      CALL SOL_LS( TX,XBSMX )
!
!---  Salt mass fraction in brine  ---
!
      XBSX = MIN( YLSX,XBSMX )
!
!---  Water vapor pressure in equilibrium with brine at saturated
!     conditions  ---
!
      CALL SP_B( TX,XBSX,PSBX )
!
!---  Water vapor pressure adjusted for aqueous saturation
!     in geologic media  ---
!
      CALL VPL( PSBX,PVWX,SLX,TX,IZN )
!
!---  Nonaqueous mole fractions, densities, and saturations  ---
!
      CALL BO_GOR( PX,POSMX,TX )
      CALL BO_EQ_GAS( PX,PVWX,RHOGX,TX,XGCX,XGWX )
      CALL BO_EQ_OIL( RHONX,PX,POSMX,TX,XNCX )
      SNX = 1.D+0-SGX-SLX
!
!---  Brine density  ---
!
      CALL DENS_B( TX,PX,XBSX,RHOBX )
!
!---  Aqueous density  ---
!
      CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
!
!---  Aqueous component mass fractions  ---
!
      XLSX = XBSX
      XLWX = 1.D+0-XLSX
!
!---  Total mass of salt  ---
!
      TMSX = YLSX*RHOBX*PORDX*VOLX
!
!---  Total mass of gas  ---
!
      TMCX(1) = PORDX*VOLX*(SGX*XGCX(1)*RHOGX + SNX*XNCX(1)*RHONX)
!
!---  Total mass oil  ---
!
      TMCX(2) = PORDX*VOLX*SNX*XNCX(2)*RHONX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BO_PC5 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_MMP( PX,TX,ZMCX )
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
!     Prediction of CO2 Minimum Miscibility Pressures
!     T.H. Ahmed, 1994
!     SPE 27032
!
!     Declared variables:
!
!     TX - temperature, C
!     ZMCX - total mole fraction of petroleum components
!            or total mass fraction of petroleum components
!
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 17 October 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GLB_PAR
      USE FILES
      USE EOR
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Parameter Statements----------------------------!
!
      PARAMETER(LMCX=1000)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XMCX(LNGC),YMCX(LNGC),ZMCX(LNGC)
      REAL*8 XMCZ(LNGC),YMCZ(LNGC)
      REAL*8 XMCIX(LNGC,LMCX),YMCIX(LNGC,LMCX),ZMCIX(LNGC,LMCX)
      REAL*8 FKX(LNGC),BETAX,TLDX(LMCX)
      REAL*8 PKX(100),TLKX(LNGC,100)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_MMP'
!
!---  Absolute temperature  ---
!
      TKX = TX + TABS
!
!---  Initialize pressure index and number of active components  ---
!
      IPX = 0
      NGCX = 0
!
!---  Initialize gas and oil  ---
!
      DO 10 IGC = 1,NGC+2
        XMCX(IGC) = ZMCX(IGC)
        YMCX(IGC) = GC_MMP(IGC)
        IF( XMCX(IGC).GT.0.D+0 .OR. YMCX(IGC).GT.0.D+0 ) NGCX = NGCX + 1
   10 CONTINUE
!
!---  Top of pressure loop  ---
!
   20 CONTINUE
      IPX = IPX + 1
      PKX(IPX) = PX
!
!---  Initialize number of contacts  ---
!
      NCX = 0
!
!---  Top of contacts loop  ---
!
   30 CONTINUE
!
!---  Increment number of contacts  ---
!
      NCX = NCX + 1
      IF( NCX.GE.LMCX ) THEN
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Set mixing cell concentrations  ---
!
      DO 100 NX = 1,NCX
        IF( NX.EQ.1 ) THEN
          DO 40 IGC = 1,NGC+2
            YMCIX(IGC,NX) = YMCX(IGC)
   40     CONTINUE
        ENDIF
        IF( NX.EQ.NCX ) THEN
          DO 50 IGC = 1,NGC+2
            XMCIX(IGC,NX) = XMCX(IGC)
   50     CONTINUE
        ENDIF
        DO 60 IGC = 1,NGC+2
          ZMCIX(IGC,NX) = XMCIX(IGC,NX) + 
     &      5.D-1*(YMCIX(IGC,NX)-XMCIX(IGC,NX))
   60   CONTINUE
  100 CONTINUE
!
!---  Next tier of mixing cell concentrations  ---
!
      DO 200 NX = 1,NCX
        CALL STABILITY( FKX,PX,TKX,XMCZ,YMCZ,ZMCIX(1,NX),INDX )
        CALL NEG_FLASH( BETAX,CGX,CNX,FKX,PX,TKX,
     &    XMCIX(1,NX),YMCIX(1,NX+1),ZMCIX(1,NX),ZGX,ZNX )
        TLDX(NX) = 0.D+0
        DO 110 IGC = 1,NGC+2
          TLDX(NX) = TLDX(NX) + (YMCIX(IGC,NX+1)-XMCIX(IGC,NX))**2
  110   CONTINUE
        TLDX(NX) = SQRT(TLDX(NX))
  200 CONTINUE
!
!---  Find key tie-lines  ---
!
      KCX = 0
      KSX = 0
      DO 300 NX = 1,NCX-2
        DTLDX = 5.D-1*(-3.D+0*TLDX(NX) + 4.D+0*TLDX(NX+1) - TLDX(NX+2))
        DTLD2X = TLDX(NX) - 2.D+0*TLDX(NX+1) + TLDX(NX+2)
        CURVX = ABS(DTLD2X)/((1.D+0 + DTLDX**2)**(1.5D+0))
        CURV2X = MAX( ABS(TLDX(NX)-TLDX(NX+1)),ABS(TLDX(NX)-TLDX(NX+2)),
     &    ABS(TLDX(NX+1)-TLDX(NX+2)) )
!
!---    Low curvature, start search for tie-line series  ---
!
        IF( CURVX.LT.1.D-7 .AND. CURV2X.LT.1.D-4 .AND. KSX.EQ.0 ) THEN
          KSX = NX
!
!---    Low curvature ends or mixing cell tier ends,
!       find average tie-line  ---
!
        ELSEIF( ((CURVX.GE.1.D-7 .OR. CURV2X.GE.1.D-4) .AND. KSX.GT.0)
     &    .OR. (KSX.GT.0 .AND. NX.EQ.(NCX-2)) ) THEN
          KCX = KCX + 1
          TLKX(KCX,IPX) = 0.D+0
          DO 210 MX = KSX,NX+1
            TLKX(KCX,IPX) = TLKX(KCX,IPX) + TLDX(MX)
  210     CONTINUE
          TLKX(KCX,IPX) = TLKX(KCX,IPX)/REAL(NX+2-KSX)
          KSX = 0
        ENDIF
  300 CONTINUE
      IF( KCX.LT.(NGCX-1) ) GOTO 30
      IF( IPX.EQ.1 ) THEN
        WRITE(IWR,'(//,A)') ' ---  Minimum Miscibility Pressure ' // 
     &    'Output Record  ---'
        WRITE(ISC,'(//,A)') ' ---  Minimum Miscibility Pressure ' // 
     &    'Output Record  ---'
        WRITE(IWR,'(A,$)') 'Pressure, Pa '
        WRITE(ISC,'(A,$)') 'Pressure, Pa '
        DO 310 KX = 1,KCX
          WRITE(ISC,'(2X,A,I2,$)') 'Tie-Line #',KX
          WRITE(IWR,'(2X,A,I2,$)') 'Tie-Line #',KX
  310   CONTINUE
        WRITE(IWR,'(2X,A)') 'MMP Est., Pa  Number of Cells'
        WRITE(ISC,'(2X,A)') 'MMP Est., Pa  Number of Cells'
      ENDIF
      WRITE(ISC,'(1PE12.5,$)') PX
      WRITE(IWR,'(1PE12.5,$)') PX
      TLKMNX = 1.D+2
      DO 320 KX = 1,KCX
        WRITE(ISC,'(2X,1PE12.5,$)') TLKX(KX,IPX)
        WRITE(IWR,'(2X,1PE12.5,$)') TLKX(KX,IPX)
        IF( TLKX(KX,IPX).LT.TLKMNX ) THEN
          KMX = KX
          TLKMNX = TLKX(KX,IPX)
        ENDIF
  320 CONTINUE
      IF( IPX.GT.1 ) THEN
        PMMX = PX + TLKX(KMX,IPX)*(PX-PKX(IPX-1))/
     &  (TLKX(KMX,IPX-1)-TLKX(KMX,IPX))
        WRITE(ISC,'(2X,1PE12.5,2X,I3)') PMMX,NCX
        WRITE(IWR,'(2X,1PE12.5,2X,I3)') PMMX,NCX
      ELSE
        WRITE(ISC,'(A)') ' '
        WRITE(IWR,'(A)') ' '
      ENDIF
      IF( PX.GT.15.0D+6 ) THEN
        PX = PX + 5.D+4
      ELSE
        PX = PX + 5.D+5
      ENDIF
      IF( PX.LT.15.8D+6 ) GOTO 20 
      WRITE(ISC,'(2X,A)') 'Tie-Lines'
      DO 330 NX = 1,NCX
        WRITE(IWR,'(2X,1PE12.5)') TLDX(NX)
  330 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_MMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_NA( BETAX,FKX,GNIFTX,PX,PVWX,RHOMGX,RHOMNX,SGX,
     &  SNX,TX,XMGCX,XMGWX,XMNCX,ZGX,ZMCX,ZNX,IBETAX,NPHAZX )
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
!     Flash calculation for nonaqueous phases
!
!     Declared variables:
!
!     TX - temperature, C
!     PX - pressure, absolute, Pa
!     PVWX - water vapor partial pressure, Pa
!     ZMCX - total mole fraction of petroleum components
!            or total mass fraction of petroleum components
!
!     Output variables:
!
!     RHOMGX - molar gas density, kmol/m^3
!     RHOMNX - molar nonaqueous-liquid density, kmol/m^3
!     SGX - gas saturation
!     SNX - nonaqueous-liquid saturation
!     YMCX - gas mole fraction of petroleum components
!     XCMX - nonaqueous-liquid mole fraction of petroleum components
!     XMGCX - gas mole fraction of petroleum components
!     XMGWX - gas mole fraction of water
!     XMNCX - nonaqueous-liquid mole fraction of petroleum components
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 10 July 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GLB_PAR
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 AIX(LNGC),FKX(LNGC)
      REAL*8 XMCX(LNGC),YMCX(LNGC),ZMCX(LNGC)
      REAL*8 XMGCX(LNGC),XMNCX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_NA'
!
!---  Absolute temperature and pressure  ---
!
      TKX = TX + TABS
      SIGMAX = GNIFTX
!
!---  Fixed beta  ---
!
      IF( IBETAX.EQ.1 ) THEN
        IF( SIGMAX.LT.EPSL ) THEN
          IPHX = 1
          DO IGC = 1,NGC+2
            XMCX(IGC) = ZMCX(IGC)
            YMCX(IGC) = ZMCX(IGC)
          ENDDO
        ELSE
          IPHX = 2
          DO IGC = 1,NGC+2
            XMCX(IGC) = ZMCX(IGC)/(1.D+0 + BETAX*(FKX(IGC)-1.D+0))
            YMCX(IGC) = XMCX(IGC)*FKX(IGC)
          ENDDO
        ENDIF
!
!---    Soave-Redlich-Kwong temperature dependent component 
!       coefficients  ---
!
        IF( ISLC(68).EQ.2 ) THEN
          CALL SRK_CT( AIX,TKX )
!
!---    Peng-Robinson temperature dependent component 
!       coefficients  ---
!
        ELSE
          CALL PR_CT( AIX,TKX )
        ENDIF
!
!---  The Isothermal Two-Phase Flash  ---
!
      ELSE
        CALL ISOTHFL( AIX,BETAX,FKX,PX,SIGMAX,TKX,XMCX,YMCX,ZMCX,IPHX )
      ENDIF
!
!---  Gas-nonaqueous liquid interfacial tension check  ---
!
      IF( IPHX.EQ.2 ) THEN
        CALL COMP_LIQ( AIX,BNX,CNX,PX,TKX,XMCX,ZNX )
        RHOMNX = PX/(ZNX*TKX*RCU - PX*CNX)
        CALL COMP_VAP( AIX,BGX,CGX,PX,TKX,YMCX,ZGX )
        RHOMGX = PX/(ZGX*TKX*RCU - PX*CGX)
        GNIFTX = 0.D+0
        DO IGC = 1,NGC+2
          GNIFTX = GNIFTX + GCPP(22,IGC)*(XMCX(IGC)*RHOMNX 
     &      - YMCX(IGC)*RHOMGX)
        ENDDO
!
!---    Interfacial tension in dynes/cm  ---
!
        GNIFTX = (1.D-3*GNIFTX)**4
        IF( GNIFTX.LT.EPSL ) GNIFTX = 0.D+0
!
!---    Low gas-nonaqueous liquid interfacial tension with two-phase
!       index, switch to single-phase conditions  ---
!
        IF( GNIFTX.LT.1.D-2 ) THEN
          IPHX = 1
        ENDIF
      ELSE
        GNIFTX = 0.D+0
      ENDIF
!
!---  Single-phase conditions  ---
!
      IF( ABS(IPHX).EQ.1 ) THEN
        NPHAZX = (NPHAZX/100)*100 + 10 + MOD(NPHAZX,10)
        CALL COMP_LIQ( AIX,BNX,CNX,PX,TKX,XMCX,ZNX )
        RHOMNX = PX/(ZNX*TKX*RCU - PX*CNX)
        VARNX = 1.D+0/(RHOMNX*BNX)
        CALL COMP_VAP( AIX,BGX,CGX,PX,TKX,YMCX,ZGX )
        RHOMGX = PX/(ZGX*TKX*RCU - PX*CGX)
        VARGX = 1.D+0/(RHOMGX*BGX)
        VARX = 5.D-1*(VARNX+VARGX)
!
!---    Logistic function to transition between liquid and gas
!       for single-phase conditions  ---
!
        IF( IBETAX.EQ.0 ) THEN
          BETAX = 1.D+0/(1.D+0 + EXP(-1.D+1*(VARX-1.75D+0)))
        ENDIF
!!
!!---  Single-phase liquid conditions  ---
!!
!      IF( IPHX.EQ.-1 ) THEN
!        CALL COMP_LIQ( AIX,BNX,CNX,PX,TKX,XMCX,ZNX )
!        RHOMNX = PX/(ZNX*TKX*RCU - PX*CNX)
!        VARNX = 1.D+0/(RHOMNX*BNX)
!        BETAX = MIN( MAX( (VARNX-1.25D+0),0.D+0 ),1.D+0 )
!        IF( BETAX.GT.EPSL ) THEN
!          CALL COMP_VAP( AIX,BGX,CGX,PX,TKX,YMCX,ZGX )
!          RHOMGX = PX/(ZGX*TKX*RCU - PX*CGX)
!        ELSE
!          ZGX = ZNX
!          RHOMGX = RHOMNX
!        ENDIF
!!
!!---  Single-phase vapor conditions  ---
!!
!      ELSEIF( IPHX.EQ.1 ) THEN
!        CALL COMP_VAP( AIX,BGX,CGX,PX,TKX,YMCX,ZGX )
!        RHOMGX = PX/(ZGX*TKX*RCU - PX*CGX)
!        VARGX = 1.D+0/(RHOMGX*BGX)
!        BETAX = MIN( MAX( (VARGX-1.25D+0),0.D+0 ),1.D+0 )
!        IF( (1.D+0-BETAX).GT.EPSL ) THEN
!          CALL COMP_LIQ( AIX,BNX,CNX,PX,TKX,XMCX,ZNX )
!          RHOMNX = PX/(ZNX*TKX*RCU - PX*CNX)
!        ELSE
!          ZNX = ZGX
!          RHOMNX = RHOMGX
!        ENDIF
!
!---  Two-phase conditions  ---
!
      ELSEIF( IPHX.EQ.2 ) THEN
        NPHAZX = (NPHAZX/100)*100 + 20 + MOD(NPHAZX,10)
        CALL COMP_LIQ( AIX,BNX,CNX,PX,TKX,XMCX,ZNX )
        RHOMNX = PX/(ZNX*TKX*RCU - PX*CNX)
        CALL COMP_VAP( AIX,BGX,CGX,PX,TKX,YMCX,ZGX )
        RHOMGX = PX/(ZGX*TKX*RCU - PX*CGX)
        VARNX = 1.D+0/(RHOMNX*BNX)
        VARGX = 1.D+0/(RHOMGX*BGX)
        VARX = BETAX*VARGX + (1.D+0-BETAX)*VARNX
        IF( IBETAX.EQ.0 ) THEN
          BETA1X = 1.D+0/(1.D+0 + EXP(-1.D+1*(VARX-1.75D+0)))
          BETA2X = BETAX
          FX = 1.D+0 - EXP(-5.D-1*GNIFTX)
          BETAX = (1.D+0-FX)*BETA1X + FX*BETA2X
        ENDIF
!
!---  Unknown phase condition  ---
!
      ELSE
        INDX = 14
        CHMSG = 'Isothermal Flash/Stability Failure: Beta = '
        RLMSG = BETAX
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Interfacial tension in N/m  ---
!
      GNIFTX = 1.D+3*GNIFTX
!
!---  Water vapor density, kg/m^3 ---
!
      ISRX = 2
      CALL DENS_W( TX,PVWX,RHOLWX,RHOGWX,ISRX )
!
!---  Water vapor molar density, kmol/m^3 ---
!
      RHOMGWX = RHOGWX/WTMW
!
!---  Include water in gas mole fractions and molar density  ---
!
      DO 120 IGC = 1,NGC+2
        XMGCX(IGC) = YMCX(IGC)*RHOMGX/(RHOMGX+RHOMGWX)
        XMNCX(IGC) = XMCX(IGC)
  120 CONTINUE
      XMGWX = RHOMGWX/(RHOMGX+RHOMGWX)
      RHOMGX = RHOMGX+RHOMGWX
      SGX = (BETAX*RHOMNX)/(BETAX*RHOMNX + 
     &  (1.D+0-BETAX)*(1.D+0-XMGWX)*RHOMGX)
      TOLX = 1.D-12
      IF( SGX.LT.TOLX ) SGX = 0.D+0
      IF( SGX.GT.(1.D+0-TOLX)) SGX = 1.D+0
      SNX = 1.D+0 - SGX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_NA group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_IC1( PX,PVAX,RHOLX,RHOMLX,TMSX,TX,XLAX,XLSX,
     &  XLWX,XMLAX,XMLSX,XMLWX,YLSX,ZMCX,IC7X,IC8X )
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
!     SL = 1.0
!     SG+SN = 0.0
!     Aqueuous saturated
!
!     Input variables:
!
!     TX - temperature, C
!     PX - pressure, abs, Pa
!     PVAX - aqueous CO2 relative saturation, or
!           aqueous CO2 mass fraction
!     TMSX - aqueous salt relative saturation, or
!         aqueous salt mass fraction
!
!     Output variables:
!
!     RHOLX - aqueous density, kg/m^3
!     RHOMLX - aqueous molar density, kg/m^3
!     TMCX - total mass of petroleum components, kg
!     XLAX - aqueous mass fraction of CO2
!     XLSX - aqueous mass fraction of salt
!     XLWX - aqueous mass fraction of water
!     YLSX - total salt brine mass fraction
!     ZMCX - total mole fraction of petroleum components
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 05 July 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GLB_PAR
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 ZMCX(LNGC)
      REAL*8 AJ(2,2),BJ(2),GX(2,3)
      INTEGER IJ(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_IC1'
!
!---  Absolute temperature  ---
!
      TKX = TX + TABS
!
!---  Aqueous salt mass fraction and aqueous CO2 mass fraction 
!     specified, solve for PVAX and XBSX  ---
!
      IF( IC8X.EQ.3 .AND. IC7X.EQ.3 ) THEN
        XLSX = TMSX
!
!---    Guest salt brine mass fraction  ---
!
        XBSX = XLSX
        XLAX = PVAX
        PVAX = PX
        CALL EQUIL( TX,PX,PVAX,XBSX,XLAZ,XLSZ,XLWZ )
!
!---    Guest CO2 partial pressure  ---
!
        PVAX = PX*(XLAX/XLAZ)
        DPY = 1.D-1
        DXY = 1.D-7
!
!---    Top of Newton-Raphson loop for salt brine mass fraction
!       and CO2 partial pressure  ---
!
  100   CONTINUE
        DO 110 M = 1,3
          PVAY = PVAX
          XBSY = XBSX
          IF( M.EQ.2 ) PVAY = PVAX + DPY
          IF( M.EQ.3 ) XBSY = XBSX + DXY
          CALL SP_B( TX,XBSY,PSBX )
          CALL EQUIL( TX,PX,PVAY,XBSY,XLAY,XLSY,XLWY )
          GX(1,M) = XLAX - XLAY
          GX(2,M) = XLSX - XLSY
  110   CONTINUE
!
!---    Load solution vector and Jacobian matrix  ---
!
        DO 120 M = 1,2
          AJ(M,1) = (GX(M,2)-GX(M,1))/DPY
          AJ(M,2) = (GX(M,3)-GX(M,1))/DXY
          BJ(M) = -GX(M,1)
 120    CONTINUE
        IF( ABS(AJ(1,1))/EPSL.LT.EPSL .OR. 
     &    ABS(AJ(2,2))/EPSL.LT.EPSL ) THEN
          PRINT *,'Aqueous salt mass fraction and aqueous ' // 
     &      'CO2 mass fraction: PVAX = ',PVAX,' XBSX = ',XBSX
          PRINT *,'XLAX = ',XLAX,' XLAY = ',XLAY
          PRINT *,'XLSX = ',XLSX,' XLSY = ',XLSY
          STOP
        ENDIF
!
!---    Solve linear system  ---
!
        JP = 2
        KP = 2
        CALL LUDCMP( AJ,JP,KP,IJ,DJ )
        CALL LUBKSB( AJ,JP,KP,IJ,BJ )
!
!---    Update primary unknowns  ---
!
        DPGAX = BJ(1)
        DXBSX = BJ(2)
        PVAX = PVAX + DPGAX
        XBSX = XBSX + DXBSX
        PVAX = MAX( PVAX,0.D+0 )
        XBSX = MIN( MAX( XBSX,0.D+0 ),1.D+0 )
!
!---    Convergence check  ---
!
        IF( ABS(DPGAX).GT.1.D-2 ) GOTO 100
        IF( ABS(DXBSX).GT.1.D-8 ) GOTO 100
!
!---    Gas pressures exceed aqueous pressure  ---
!
        IF( PSBX+PVAX.GT.PX ) THEN
          INDX = 10
          CHMSG = 'Initial Condition State 1: ' //
     &      'Gas Pressure > Total Pressure'
          RLMSG = PSBX + PVAX
          CALL WRMSGS( INDX )
        ENDIF
!
!---  Aqueous salt mass fraction and aqueous CO2 relative saturation 
!     specified, solve for PVAX and XBSX  ---
!
      ELSEIF( IC8X.EQ.3 .AND. IC7X.EQ.2 ) THEN
        XLSX = TMSX
!
!---    Guest salt brine mass fraction  ---
!
        XBSX = XLSX
        WLAX = PVAX
!
!---    Guest CO2 partial pressure  ---
!
        PVAX = WLAX*PX
        DPY = 1.D-1
        DXY = 1.D-7
!
!---    Top of Newton-Raphson loop for salt brine mass fraction
!       and CO2 partial pressure  ---
!
  200   CONTINUE
        DO 210 M = 1,3
          PVAY = PVAX
          XBSY = XBSX
          IF( M.EQ.2 ) PVAY = PVAX + DPY
          IF( M.EQ.3 ) XBSY = XBSX + DXY
          CALL SP_B( TX,XBSY,PSBX )
          CALL EQUIL( TX,PX,PVAY,XBSY,XLAY,XLSY,XLWY )
          CALL EQUIL( TX,PX,PX,XBSY,XLAZ,XLSZ,XLWZ )
          IF( XLAZ.GT.EPSL ) THEN
            WLAY = XLAY/XLAZ
          ELSE
            WLAY = 0.D+0
          ENDIF
          GX(1,M) = WLAX - WLAY
          GX(2,M) = XLSX - XLSY
  210   CONTINUE
!
!---    Load solution vector and Jacobian matrix  ---
!
        DO 220 M = 1,2
          AJ(M,1) = (GX(M,2)-GX(M,1))/DPY
          AJ(M,2) = (GX(M,3)-GX(M,1))/DXY
          BJ(M) = -GX(M,1)
 220    CONTINUE
        IF( ABS(AJ(1,1))/EPSL.LT.EPSL .OR. 
     &    ABS(AJ(2,2))/EPSL.LT.EPSL ) THEN
          PRINT *,'Aqueous salt mass fraction and aqueous ' // 
     &      'CO2 relative saturation: PVAX = ',PVAX,' XBSX = ',XBSX
          PRINT *,'WLAX = ',WLAX,' WLAY = ',WLAY
          PRINT *,'XLSX = ',XLSX,' XLSY = ',XLSY
          STOP
        ENDIF
!
!---    Solve linear system  ---
!
        JP = 2
        KP = 2
        CALL LUDCMP( AJ,JP,KP,IJ,DJ )
        CALL LUBKSB( AJ,JP,KP,IJ,BJ )
!
!---    Update primary unknowns  ---
!
        DPGAX = BJ(1)
        DXBSX = BJ(2)
        PVAX = PVAX + DPGAX
        XBSX = XBSX + DXBSX
        PVAX = MAX( PVAX,0.D+0 )
        XBSX = MIN( MAX( XBSX,0.D+0 ),1.D+0 )
!
!---    Convergence check  ---
!
        IF( ABS(DPGAX).GT.1.D-2 ) GOTO 200
        IF( ABS(DXBSX).GT.1.D-8 ) GOTO 200
!
!---  Aqueous salt relative saturation and aqueous CO2 mass fraction 
!     specified, solve for PVAX and XBSX  ---
!
      ELSEIF( IC8X.EQ.2 .AND. IC7X.EQ.3 ) THEN
        WLSX = TMSX
        CALL SOL_LS( TX,XBSMX )
!
!---    Guest salt brine mass fraction  ---
!
        XBSX = WLSX*XBSMX
        XLAX = PVAX
        PVAX = PX
        CALL EQUIL( TX,PX,PVAX,XBSX,XLAZ,XLSZ,XLWZ )
!
!---    Guest CO2 partial pressure  ---
!
        PVAX = PX*(XLAX/XLAZ)
        DPY = 1.D-1
        DXY = 1.D-7
!
!---    Top of Newton-Raphson loop for salt brine mass fraction
!       and CO2 partial pressure  ---
!
  300   CONTINUE
        DO 310 M = 1,3
          PVAY = PVAX
          XBSY = XBSX
          IF( M.EQ.2 ) PVAY = PVAX + DPY
          IF( M.EQ.3 ) XBSY = XBSX + DXY
          CALL SP_B( TX,XBSY,PSBX )
          CALL EQUIL( TX,PX,PVAY,XBSY,XLAY,XLSY,XLWY )
          CALL EQUIL( TX,PX,PVAY,XBSMX,XLAZ,XLSZ,XLWZ )
          WLSY = XLSY/XLSZ
          GX(1,M) = XLAX - XLAY
          GX(2,M) = WLSX - WLSY
  310   CONTINUE
!
!---    Load solution vector and Jacobian matrix  ---
!
        DO 320 M = 1,2
          AJ(M,1) = (GX(M,2)-GX(M,1))/DPY
          AJ(M,2) = (GX(M,3)-GX(M,1))/DXY
          BJ(M) = -GX(M,1)
 320    CONTINUE
        IF( ABS(AJ(1,1))/EPSL.LT.EPSL .OR. 
     &    ABS(AJ(2,2))/EPSL.LT.EPSL ) THEN
          PRINT *,'Aqueous salt relative saturation and aqueous ' // 
     &      'CO2 mass fraction: PVAX = ',PVAX,' XBSX = ',XBSX
          PRINT *,'XLAX = ',XLAX,' XLAY = ',XLAY
          PRINT *,'WLSX = ',WLSX,' WLSY = ',WLSY
          STOP
        ENDIF
!
!---    Solve linear system  ---
!
        JP = 2
        KP = 2
        CALL LUDCMP( AJ,JP,KP,IJ,DJ )
        CALL LUBKSB( AJ,JP,KP,IJ,BJ )
!
!---    Update primary unknowns  ---
!
        DPGAX = BJ(1)
        DXBSX = BJ(2)
        PVAX = PVAX + DPGAX
        XBSX = XBSX + DXBSX
        PVAX = MAX( PVAX,0.D+0 )
        XBSX = MIN( MAX( XBSX,0.D+0 ),1.D+0 )
!
!---    Convergence check  ---
!
        IF( ABS(DPGAX).GT.1.D-2 ) GOTO 300
        IF( ABS(DXBSX).GT.1.D-8 ) GOTO 300
!
!---    Gas pressures exceed aqueous pressure  ---
!
        IF( PSBX+PVAX.GT.PX ) THEN
          INDX = 10
          CHMSG = 'Initial Condition State 1: ' //
     &      'Gas Pressure > Total Pressure'
          RLMSG = PSBX + PVAX
          CALL WRMSGS( INDX )
        ENDIF
!
!---  Aqueous salt relative saturation and aqueous CO2 relative  
!     saturation specified, solve for PVAX and XBSX  ---
!
      ELSEIF( IC8X.EQ.2 .AND. IC7X.EQ.2 ) THEN
        WLSX = TMSX
        CALL SOL_LS( TX,XBSMX )
!
!---    Guest salt brine mass fraction  ---
!
        XBSX = WLSX*XBSMX
        WLAX = PVAX
!
!---    Guest CO2 partial pressure  ---
!
        PVAX = WLAX*PX
        DPY = 1.D-1
        DXY = 1.D-7
!
!---    Top of Newton-Raphson loop for salt brine mass fraction
!       and CO2 partial pressure  ---
!
  400   CONTINUE
        DO 410 M = 1,3
          PVAY = PVAX
          XBSY = XBSX
          IF( M.EQ.2 ) PVAY = PVAX + DPY
          IF( M.EQ.3 ) XBSY = XBSX + DXY
          CALL SP_B( TX,XBSY,PSBX )
          CALL EQUIL( TX,PX,PVAY,XBSY,XLAY,XLSY,XLWY )
          CALL EQUIL( TX,PX,PX,XBSY,XLAZ,XLSZ,XLWZ )
          IF( XLAZ.GT.EPSL ) THEN
            WLAY = XLAY/XLAZ
          ELSE
            WLAY = 0.D+0
          ENDIF
          CALL EQUIL( TX,PX,PVAY,XBSMX,XLAZ,XLSZ,XLWZ )
          IF( XLSZ.GT.EPSL ) THEN
            WLSY = XLSY/XLSZ
          ELSE
            WLSY = 0.D+0
          ENDIF
          GX(1,M) = WLAX - WLAY
          GX(2,M) = WLSX - WLSY
  410   CONTINUE
!
!---    Load solution vector and Jacobian matrix  ---
!
        DO 420 M = 1,2
          AJ(M,1) = (GX(M,2)-GX(M,1))/DPY
          AJ(M,2) = (GX(M,3)-GX(M,1))/DXY
          BJ(M) = -GX(M,1)
 420    CONTINUE
        IF( ABS(AJ(1,1))/EPSL.LT.EPSL .OR. 
     &    ABS(AJ(2,2))/EPSL.LT.EPSL ) THEN
          PRINT *,'Aqueous salt relative saturation and aqueous ' // 
     &      'CO2 relative saturation: PVAX = ',PVAX,' XBSX = ',XBSX
          PRINT *,'WLAX = ',WLAX,' WLAY = ',WLAY
          PRINT *,'WLSX = ',WLSX,' WLSY = ',WLSY
          STOP
        ENDIF
!
!---    Solve linear system  ---
!
        JP = 2
        KP = 2
        CALL LUDCMP( AJ,JP,KP,IJ,DJ )
        CALL LUBKSB( AJ,JP,KP,IJ,BJ )
!
!---    Update primary unknowns  ---
!
        DPGAX = BJ(1)
        DXBSX = BJ(2)
        PVAX = PVAX + DPGAX
        XBSX = XBSX + DXBSX
        PVAX = MAX( PVAX,0.D+0 )
        XBSX = MIN( MAX( XBSX,0.D+0 ),1.D+0 )
!
!---    Convergence check  ---
!
        IF( ABS(DPGAX).GT.1.D-2 ) GOTO 400
        IF( ABS(DXBSX).GT.1.D-8 ) GOTO 400
!
!---    Gas pressures exceed aqueous pressure  ---
!
        IF( PSBX+PVAX.GT.PX ) THEN
          INDX = 10
          CHMSG = 'Initial Condition State 1: ' //
     &      'Gas Pressure > Total Pressure'
          RLMSG = PSBX + PVAX
          CALL WRMSGS( INDX )
        ENDIF
      ENDIF
!
!---  Aqueous mass fractions  ---
!
      CALL SP_B( TX,XBSX,PSBX )
      CALL EQUIL( TX,PX,PVAX,XBSX,XLAX,XLSX,XLWX )
      WTMLX = XLAX/WTMA + XLSX/WTMS + XLWX/WTMW
      XMLAX = XLAX/WTMA/WTMLX
      XMLSX = XLSX/WTMS/WTMLX
      XMLWX = XLWX/WTMW/WTMLX
!      XMLAX = XMLAX*PVAX/PX
!      XMLSX = (1.D+0-XMLAX)*XMLSX/(XMLSX+XMLWX)
!      XMLWX = (1.D+0-XMLAX)*XMLWX/(XMLSX+XMLWX)
      WTMLX = XMLAX*WTMA + XMLSX*WTMS + XMLWX*WTMW
      XLAX = XMLAX*WTMA/WTMLX
      XLSX = XMLSX*WTMS/WTMLX
      XLWX = XMLWX*WTMW/WTMLX
!
!---  Brine density based on CO2-free brine  ---
!
      CALL DENS_B( TX,PX,XBSX,RHOBX )
      YLSX = XBSX
!
!---  Aqueous density based on dissoved CO2 and salt  ---
!
      CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
!
!---  Aqueous molar density  ---
!
      RHOMLX = RHOLX*WTMLX
!
!---  Total petroleum component mass fractions  ---
!
      ZMCX(1) = 1.D+0
      DO 500 IGC = 2,NGC+2
        ZMCX(IGC) = 0.D+0
  500 CONTINUE
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
      SUBROUTINE FLH_IC2( PX,RHOGX,RHOMGX,RHOLX,RHOMLX,RHONX,RHOMNX,
     &  SGX,SLX,SNX,TMSX,TX,XGCX,XMGCX,XGWX,XMGWX,XLAX,XLSX,XLWX,
     &  XMLAX,XMLSX,XMLWX,XNCX,XMNCX,YLSX,ZGX,ZMCX,ZNX,IC8X,IZN,NPHAZX )
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
!     SL > 1.0
!     SG+SN > 0.0
!     Three-phase (total petroleum components)
!
!     Declared variables:
!
!     TX - temperature, C
!     PX - pressure, absolute, Pa
!     SLX - aqueous saturation
!     ZMCX - total mole fraction of petroleum components
!     TMSX - aqueous salt relative saturation, or
!         aqueous salt mass fraction
!
!     Output variables:
!
!     PGX - gas pressure, guage Pa
!     PLX - aqueous pressure, guage, Pa
!     PNX - nonaqueous liquid pressure, guage Pa
!     RHOGX - gas density, kg/m^3
!     RHOMGX - gas molar density, kg/m^3
!     RHOLX - aqueous density, kg/m^3
!     RHOMLX - aqueous molar density, kg/m^3
!     RHONX - nonaqueous density, kg/m^3
!     RHOMNX - nonaqueous molar density, kg/m^3
!     SGX - gas saturation
!     SLX - aqueous saturation
!     SNX - nonaqueous liquid saturation
!     XGCX - gas mass fraction of petroleum components
!     XMGCX - gas mole fraction of petroleum components
!     XGWX - gas mass fraction of water
!     XMGWX - gas mole fraction of water
!     XLAX - aqueous mass fraction of CO2
!     XLSX - aqueous mass fraction of salt
!     XLWX - aqueous mass fraction of water
!     XNCX - nonaqueous liquid mass fraction of petroleum components
!     XMNCX - nonaqueous liquid mole fraction of petroleum components
!     YLSX - brine mass fraction of total salt
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 July 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GLB_PAR
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 ZMCX(LNGC),FKX(LNGC)
      REAL*8 XGCX(LNGC),XNCX(LNGC)
      REAL*8 XMGCX(LNGC),XMNCX(LNGC)
      REAL*8 GX(1,2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_IC2'
!
!---  Absolute temperature  ---
!
      TKX = TX + TABS
!
!---  Aqueous salt mass fraction, solve for XBSX  ---
!
      IF( IC8X.EQ.3 ) THEN
        XLSX = TMSX
!
!---    Guest salt brine mass fraction  ---
!
        XBSX = XLSX
        DXY = 1.D-7
!
!---    Top of Newton-Raphson loop for salt brine mass fraction  ---
!
  100   CONTINUE
        DO 110 M = 1,2
          XBSY = XBSX
          IF( M.EQ.2 ) XBSY = XBSX + DXY
          CALL SP_B( TX,XBSY,PSBX )
!
!---      Vapor-pressure lowering  ---
!
          IF( IZN.NE.0 ) THEN
            CALL VPL( PSBX,PVWX,SLX,TX,IZN )
!
!---      No vapor-pressure lowering, used in source and well model  ---
!
          ELSE
            PVWX = PSBX
          ENDIF
          PVX = PX - PVWX
          GNIFTX = 0.D+0
          IBETAX = 0
          CALL FLH_NA( BETAX,FKX,GNIFTX,PVX,PVWX,RHOMGX,RHOMNX,SGX,
     &      SNX,TX,XMGCX,XMGWX,XMNCX,ZGX,ZMCX,ZNX,IBETAX,NPHAZX )
          PVAX = PVX*XMGCX(1)
          CALL SP_B( TX,XBSY,PSBX )
          PTX = PSBX + PVAX
          CALL EQUIL( TX,PTX,PVAX,XBSY,XLAY,XLSY,XLWY )
          GX(1,M) = XLSX - XLSY
  110   CONTINUE
!
!---    Load solution vector and Jacobian matrix  ---
!
        DFX = (GX(1,2)-GX(1,1))/DXY
        FX = -GX(1,1)
!
!---    Update primary unknowns  ---
!
        DXBSX = FX/DFX
        XBSX = XBSX + DXBSX
        XBSX = MIN( MAX( XBSX,0.D+0 ),1.D+0 )
!
!---    Convergence check  ---
!
        IF( ABS(DXBSX).GT.1.D-8 ) GOTO 100
!
!---  Aqueous salt relative saturation, solve for XBSX  ---
!
      ELSEIF( IC8X.EQ.2 ) THEN
        WLSX = TMSX
        CALL SOL_LS( TX,XBSMX )
!
!---    Guest salt brine mass fraction  ---
!
        XBSX = WLSX*XBSMX
        DXY = 1.D-7
!
!---    Top of Newton-Raphson loop for salt brine mass fraction  ---
!
  200   CONTINUE
        DO 210 M = 1,2
          XBSY = XBSX
          IF( M.EQ.2 ) XBSY = XBSX + DXY
          CALL SP_B( TX,XBSY,PSBX )
!
!---      Vapor-pressure lowering  ---
!
          IF( IZN.NE.0 ) THEN
            CALL VPL( PSBX,PVWX,SLX,TX,IZN )
!
!---      No vapor-pressure lowering, used in source and well model  ---
!
          ELSE
            PVWX = PSBX
          ENDIF
          PVX = PX - PVWX
          GNIFTX = 0.D+0
          IBETAX = 0
          CALL FLH_NA( BETAX,FKX,GNIFTX,PVX,PVWX,RHOMGX,RHOMNX,SGX,
     &      SNX,TX,XMGCX,XMGWX,XMNCX,ZGX,ZMCX,ZNX,IBETAX,NPHAZX )
          PVAX = PVX*XMGCX(1)
          PTX = PSBX + PVAX
          CALL EQUIL( TX,PTX,PVAX,XBSY,XLAY,XLSY,XLWY )
          CALL EQUIL( TX,PTX,PVAX,XBSMX,XLAZ,XLSZ,XLWZ )
          WLSY = XLSY/XLSZ
          GX(1,M) = WLSX - WLSY
  210   CONTINUE
!
!---    Load solution vector and Jacobian matrix  ---
!
        DFX = (GX(1,2)-GX(1,1))/DXY
        FX = -GX(1,1)
!
!---    Update primary unknowns  ---
!
        DXBSX = FX/DFX
        XBSX = XBSX + DXBSX
        XBSX = MIN( MAX( XBSX,0.D+0 ),1.D+0 )
!
!---    Convergence check  ---
!
        IF( ABS(DXBSX).GT.1.D-8 ) GOTO 200
      ENDIF
!
!---  Aqueous mass fractions  ---
!
      CALL SP_B( TX,XBSX,PSBX )
!
!---  Vapor-pressure lowering  ---
!
      IF( IZN.NE.0 ) THEN
        CALL VPL( PSBX,PVWX,SLX,TX,IZN )
!
!---  No vapor-pressure lowering, used in source and well model  ---
!
      ELSE
        PVWX = PSBX
      ENDIF
      PVX = PX - PVWX
      GNIFTX = 0.D+0
      IBETAX = 0
      CALL FLH_NA( BETAX,FKX,GNIFTX,PVX,PVWX,RHOMGX,RHOMNX,SGX,SNX,
     &  TX,XMGCX,XMGWX,XMNCX,ZGX,ZMCX,ZNX,IBETAX,NPHAZX )
      SGX = (1.D+0-SLX)*SGX
      SNX = (1.D+0-SLX)*SNX
      PVAX = PVX*XMGCX(1)
      CALL SP_B( TX,XBSX,PSBX )
      PTX = PSBX + PVAX
      CALL EQUIL( TX,PTX,PVAX,XBSX,XLAX,XLSX,XLWX )
      WTMLX = XLAX/WTMA + XLSX/WTMS + XLWX/WTMW
      XMLAX = XLAX/WTMA/WTMLX
      XMLSX = XLSX/WTMS/WTMLX
      XMLWX = XLWX/WTMW/WTMLX
!      XMLAX = XMLAX*PVAX/PTX
!      XMLSX = (1.D+0-XMLAX)*XMLSX/(XMLSX+XMLWX)
!      XMLWX = (1.D+0-XMLAX)*XMLWX/(XMLSX+XMLWX)
      WTMLX = XMLAX*WTMA + XMLSX*WTMS + XMLWX*WTMW
      XLAX = XMLAX*WTMA/WTMLX
      XLSX = XMLSX*WTMS/WTMLX
      XLWX = XMLWX*WTMW/WTMLX
!
!---  Brine density based on CO2-free brine  ---
!
      CALL DENS_B( TX,PX,XBSX,RHOBX )
      YLSX = XBSX
!
!---  Aqueous density based on dissoved CO2 and salt  ---
!
      CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
!
!---  Aqueous molar density  ---
!
      RHOMLX = RHOLX*WTMLX
!
!---  Mass fractions of the nonaqueous phases  ---
!
      WTMGX = XMGWX*WTMW
      WTMNX = 0.D+0
      DO 310 IGC = 1,NGC+2
        WTMGX = WTMGX + XMGCX(IGC)*GCPP(1,IGC)
        WTMNX = WTMNX + XMNCX(IGC)*GCPP(1,IGC)
  310 CONTINUE
      XGWX = XMGWX*WTMW/WTMGX
      DO 320 IGC = 1,NGC+2
        XGCX(IGC) = XMGCX(IGC)*GCPP(1,IGC)/WTMGX
        XNCX(IGC) = XMNCX(IGC)*GCPP(1,IGC)/WTMNX
  320 CONTINUE
!
!---  Density of the nonaqueous phases  ---
!
      RHOGX = RHOMGX*WTMGX
      RHONX = RHOMNX*WTMNX
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
      SUBROUTINE FLH_IC3( PVWX,PX,RHOGX,RHOMGX,RHONX,RHOMNX,SGX,SNX,TX,
     &  XGCX,XGWX,XMGCX,XMGWX,XNCX,XMNCX,YLSX,ZGX,ZMCX,ZNX,IZN,NPHAZX )
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
!     SL = 0.0
!     SG+SN = 1.0
!     Gas and nonaqueous-liquid saturated (total petroleum components)
!
!     Declared variables:
!
!     TX - temperature, C
!     PX - pressure, absolute, Pa
!     PVWX - water vapor relative saturation
!     ZMCX - total mole fraction of petroleum components
!            or total mass fraction of petroleum components
!     YLSX - 
!
!     Output variables:
!
!     RHOGX - gas density, kmol/m^3
!     RHOMGX - molar gas density, kmol/m^3
!     RHONX - nonaqueous-liquid density, kmol/m^3
!     RHOMNX - molar nonaqueous-liquid density, kmol/m^3
!     SGX - gas saturation
!     SNX - nonaqueous-liquid saturation
!     XGCX - gas mass fraction of petroleum components
!     XMGCX - gas mole fraction of petroleum components
!     XNCX - nonaqueous-liquid mass fraction of petroleum components
!     XMNCX - nonaqueous-liquid mole fraction of petroleum components
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 July 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GLB_PAR
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 ZMCX(LNGC),FKX(LNGC)
      REAL*8 XGCX(LNGC),XMGCX(LNGC),XNCX(LNGC),XMNCX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_IC3'
!
!---  Absolute temperature  ---
!
      TKX = TX + TABS
!
!---  Maximum solubility of salt in pure water  ---
!
      CALL SOL_LS( TX,XBSMX )
!
!---  Salt mass fraction in brine  ---
!
      XBSX = MIN( YLSX,XBSMX )
!
!---  Water vapor pressure in equilibrium with brine at saturated
!     conditions  ---
!
      CALL SP_B( TX,XBSX,PSBX )
!
!---  Water vapor pressure adjusted for aqueous saturation
!     in geologic media  ---
!
      SLX = 0.D+0
!
!---  Water vapor relative saturation  ---
!
      PSBX = PSBX*PVWX
!
!---  Vapor-pressure lowering  ---
!
      IF( IZN.NE.0 ) THEN
        CALL VPL( PSBX,PVWX,SLX,TX,IZN )
!
!---  No vapor-pressure lowering, used in source and well model  ---
!
      ELSE
        PVWX = PSBX
      ENDIF
!
!---  Total nonaqueous vapor pressure  ---
!
      PVX = PX - PVWX
!
!---  Nonaqueous mole fractions, densities, and saturations  ---
!
      GNIFTX = 0.D+0
      IBETAX = 0
      CALL FLH_NA( BETAX,FKX,GNIFTX,PVX,PVWX,RHOMGX,RHOMNX,SGX,SNX,
     &  TX,XMGCX,XMGWX,XMNCX,ZGX,ZMCX,ZNX,IBETAX,NPHAZX )
!
!---  Mass fractions of the nonaqueous phases  ---
!
      WTMGX = XMGWX*WTMW
      WTMNX = 0.D+0
      DO 310 IGC = 1,NGC+2
        WTMGX = WTMGX + XMGCX(IGC)*GCPP(1,IGC)
        WTMNX = WTMNX + XMNCX(IGC)*GCPP(1,IGC)
  310 CONTINUE
      XGWX = XMGWX*WTMW/WTMGX
      DO 320 IGC = 1,NGC+2
        XGCX(IGC) = XMGCX(IGC)*GCPP(1,IGC)/WTMGX
        XNCX(IGC) = XMNCX(IGC)*GCPP(1,IGC)/WTMNX
  320 CONTINUE
!
!---  Density of the nonaqueous phases  ---
!
      RHOGX = RHOMGX*WTMGX
      RHONX = RHOMNX*WTMNX
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
      SUBROUTINE FLH_PC1( BETAX,FKX,GNIFTX,PORDX,PVAX,PVWX,PX,
     &  RHOGX,RHOLX,RHOMGX,RHOMLX,RHOMNX,RHONX,
     &  SGX,SLX,SNX,TMCX,TMSX,TX,VOLX,XGCX,
     &  XGWX,XLAX,XLSX,XLWX,XNCX,XMGCX,XMGWX,XMLAX,XMLSX,XMLWX,XMNCX,
     &  YLSX,ZGX,ZMCX,ZNX,IBETAX,IZN,NPHAZX )
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
!     SL = 1.0, SG+SN = 0.0
!
!     energy - temperature, C - T
!     water mass - aqueous (system) pressure, Pa gauge - PL,PSO
!     CO2 mass - CO2 vapor pressure - PVA
!     CH4 mass - total CH4 mass - TMC
!     petroleum component mass - total petroleum comp. mass - TMC
!     salt mass - total salt brine mass fraction - YLS
!
!     Output variables:
!
!     RHOGX - gas density, kg/m^3
!     RHOLX - aqueous density, kg/m^3
!     RHONX - nonaqueous-liquid density, kg/m^3
!     RHOMGX - gas molar density, kmol/m^3
!     RHOMLX - aqueous dmolar density, kmol/m^3
!     RHOMNX - nonaqueous-liquid molar density, kmol/m^3
!     SGX - gas saturation
!     SNX - nonaqueous-liquid saturation
!     XGCX - gas mass fraction of petroleum components
!     XGWX - gas mass fraction of water
!     XLAX - aqueous mass fraction of CO2
!     XLSX - aqueous mass fraction of salt
!     XLWX - aqueous mass fraction of water
!     XNCX - nonaqueous-liquid mass fraction of petroleum components
!     ZMCX - total mole fraction of petroleum components
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GLB_PAR
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 TMCX(LNGC),XGCX(LNGC),XMGCX(LNGC),XMNCX(LNGC)
      REAL*8 XNCX(LNGC),ZMCX(LNGC),FKX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_PC1'
!
!---  Maximum solubility of salt in pure water  ---
!
      CALL SOL_LS( TX,XBSMX )
!
!---  Salt mass fraction in brine  ---
!
      XBSX = MIN( YLSX,XBSMX )
!
!---  Solubility of CO2 in aqueous phase  ---
!
      CALL SP_B( TX,XBSX,PSBX )
      CALL EQUIL( TX,PX,PVAX,XBSX,XLAX,XLSX,XLWX )
!
!---  Brine density  ---
!
      CALL DENS_B( TX,PX,XBSX,RHOBX )
!
!---  Aqueous density  ---
!
      CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
!
!---  Aqueous component mole fractions  ---
!
      WTMLX = XLAX/WTMA + XLSX/WTMS + XLWX/WTMW
      XMLAX = (XLAX/WTMA)/WTMLX
      XMLSX = (XLSX/WTMS)/WTMLX
      XMLWX = (XLWX/WTMW)/WTMLX
!
!---  Aqueous molar density  ---
!
      RHOMLX = RHOLX*WTMLX
!
!---  Total mass of CO2  ---
!
      TMCX(1) = XLAX*RHOLX*PORDX*VOLX
!
!---  Total mass of salt  ---
!
      TMSX = YLSX*RHOBX*PORDX*VOLX
!
!---  Total petroleum molar mass  ---
!
      STMMCX = 0.D+0
      DO 10 IGC = 1,NGC+2
        STMMCX = STMMCX + TMCX(IGC)/GCPP(1,IGC)
   10 CONTINUE
      IF( STMMCX.GT.EPSL ) THEN
        ZMCMX = 0.D+0
        DO 20 IGC = 1,NGC+2
          ZMCX(IGC) = TMCX(IGC)/GCPP(1,IGC)/STMMCX
   20   CONTINUE
      ELSE
        ZMCX(1) = 1.D+0
        DO 30 IGC = 2,NGC+2
          ZMCX(IGC) = 0.D+0
   30   CONTINUE
      ENDIF
!
!---  Nonaqueous mole fractions, densities, and saturations  ---
!
      CALL SP_B( TX,XBSX,PSBX )
      CALL VPL( PSBX,PVWX,SLX,TX,IZN )
      PVX = PX - PVWX
      CALL FLH_NA( BETAX,FKX,GNIFTX,PVX,PVWX,RHOMGX,RHOMNX,SGX,SNX,
     &  TX,XMGCX,XMGWX,XMNCX,ZGX,ZMCX,ZNX,IBETAX,NPHAZX )
      SGX = (1.D+0-SLX)*SGX
      SNX = (1.D+0-SLX)*SNX
!
!---  Nonaqueous mass fractions  ---
!
      WTMGX = XMGWX*WTMW
      WTMNX = 0.D+0
      DO 50 IGC = 1,NGC+2
        WTMGX = WTMGX + XMGCX(IGC)*GCPP(1,IGC)
        WTMNX = WTMNX + XMNCX(IGC)*GCPP(1,IGC)
   50 CONTINUE
      XGWX = XMGWX*WTMW/WTMGX
      DO 60 IGC = 1,NGC+2
        XGCX(IGC) = XMGCX(IGC)*GCPP(1,IGC)/WTMGX
        XNCX(IGC) = XMNCX(IGC)*GCPP(1,IGC)/WTMNX
   60 CONTINUE
!
!---  Noaqueous densities  ---
!
      RHOGX = RHOMGX*WTMGX
      RHONX = RHOMNX*WTMNX
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
      SUBROUTINE FLH_PC2( BETAX,FKX,GNIFTX,PORDX,PVAX,PVWX,PX,
     &  RHOGX,RHOLX,RHOMGX,RHOMLX,RHOMNX,RHONX,
     &  SGX,SLX,SNX,TMCX,TMSX,TX,VOLX,XGCX,
     &  XGWX,XLAX,XLSX,XLWX,XNCX,XMGCX,XMGWX,XMLAX,XMLSX,XMLWX,XMNCX,
     &  YLSX,ZGX,ZMCX,ZNX,IBETAX,IZM,IZN,NPHAZX )
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
!     SL < 1.0, ST > 0.0, SG = 0.0
!
!     energy - temperature, C - T
!     water mass - aqueous saturation, SL
!     CO2 mass
!       if( ZMC(1,2,N) > max( ZMC(J,2,N) ) J /= 1
!         - nonaqueous (system) pressure, Pa gauge - PN,PSO
!       else
!         - nonaqueous mole fraction of CO2 - ZMC(1,2,N)
!       endif
!     CH4 mass
!       if( ZMC(2,2,N) > max( ZMC(J,2,N) ) J /= 2
!         - nonaqueous (system) pressure, Pa gauge - PN,PSO
!       else
!         - nonaqueous mole fraction of CO2 - ZMC(2,2,N)
!       endif
!     petroleum component mass
!       if( ZMC(I,2,N) > max( ZMC(J,2,N) ) J /= I
!         - nonaqueous (system) pressure, Pa gauge - PN,PSO
!       else
!         - nonaqueous mole fraction of CO2 - ZMC(I,2,N)
!       endif
!     salt mass - total salt brine mass fraction - YLS
!
!     Output variables:
!
!     RHOGX - gas density, kg/m^3
!     RHOLX - aqueous density, kg/m^3
!     RHONX - nonaqueous-liquid density, kg/m^3
!     RHOMGX - gas molar density, kmol/m^3
!     RHOMLX - aqueous dmolar density, kmol/m^3
!     RHOMNX - nonaqueous-liquid molar density, kmol/m^3
!     SGX - gas saturation
!     SNX - nonaqueous-liquid saturation
!     XGCX - gas mass fraction of petroleum components
!     XGWX - gas mass fraction of water
!     XLAX - aqueous mass fraction of CO2
!     XLSX - aqueous mass fraction of salt
!     XLWX - aqueous mass fraction of water
!     XNCX - nonaqueous-liquid mass fraction of petroleum components
!     ZMCX - total mole fraction of petroleum components
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GLB_PAR
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 TMCX(LNGC),XGCX(LNGC),XMGCX(LNGC),XMNCX(LNGC)
      REAL*8 XNCX(LNGC),ZMCX(LNGC),FKX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_PC2'
!
!---  Total petroleum molar mass  ---
!
      IF( IZM.NE.0 ) THEN
        ZMCX(IZM) = 1.D+0
        DO 10 IGC = 1,NGC+2
          IF( IGC.NE.IZM ) ZMCX(IZM) = ZMCX(IZM) - ZMCX(IGC)
   10   CONTINUE
        ZMCX(IZM) = MAX( MIN( ZMCX(IZM),1.D+0 ),0.D+0 )
      ENDIF
!
!---  Maximum solubility of salt in pure water  ---
!
      CALL SOL_LS( TX,XBSMX )
!
!---  Salt mass fraction in brine  ---
!
      XBSX = MIN( YLSX,XBSMX )
!
!---  Water vapor pressure in equilibrium with brine at saturated
!     conditions  ---
!
      CALL SP_B( TX,XBSX,PSBX )
!
!---  Water vapor pressure adjusted for aqueous saturation
!     in geologic media  ---
!
      CALL VPL( PSBX,PVWX,SLX,TX,IZN )
!
!---  Total nonaqueous vapor pressure  ---
!
      PVX = PX - PVWX
!
!---  Nonaqueous mole fractions, densities, and saturations  ---
!
      CALL FLH_NA( BETAX,FKX,GNIFTX,PVX,PVWX,RHOMGX,RHOMNX,SGX,SNX,
     &  TX,XMGCX,XMGWX,XMNCX,ZGX,ZMCX,ZNX,IBETAX,NPHAZX )
      SGX = (1.D+0-SLX)*SGX
      SNX = (1.D+0-SLX)*SNX
!
!---  CO2 vapor pressure  ---
!
      PVAX = PVX*XMGCX(1)
!
!---  Solubility of CO2 in aqueous phase  ---
!
      CALL SP_B( TX,XBSX,PSBX )
      CALL EQUIL( TX,PX,PVAX,XBSX,XLAX,XLSX,XLWX )
!
!---  Brine density  ---
!
      CALL DENS_B( TX,PX,XBSX,RHOBX )
!
!---  Aqueous density  ---
!
      CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
!
!---  Aqueous component mole fractions  ---
!
      WTMLX = XLAX/WTMA + XLSX/WTMS + XLWX/WTMW
      XMLAX = (XLAX/WTMA)/WTMLX
      XMLSX = (XLSX/WTMS)/WTMLX
      XMLWX = (XLWX/WTMW)/WTMLX
!
!---  Aqueous molar density  ---
!
      RHOMLX = RHOLX*WTMLX
!
!---  Nonaqueous mass fractions  ---
!
      WTMGX = XMGWX*WTMW
      WTMNX = 0.D+0
      DO 30 IGC = 1,NGC+2
        WTMGX = WTMGX + XMGCX(IGC)*GCPP(1,IGC)
        WTMNX = WTMNX + XMNCX(IGC)*GCPP(1,IGC)
   30 CONTINUE
      XGWX = XMGWX*WTMW/WTMGX
      DO 40 IGC = 1,NGC+2
        XGCX(IGC) = XMGCX(IGC)*GCPP(1,IGC)/WTMGX
        XNCX(IGC) = XMNCX(IGC)*GCPP(1,IGC)/WTMNX
   40 CONTINUE
!
!---  Noaqueous densities  ---
!
      RHOGX = RHOMGX*WTMGX
      RHONX = RHOMNX*WTMNX
!
!---  Total mass of salt  ---
!
      TMSX = YLSX*RHOBX*PORDX*VOLX
!
!---  Total mass of CO2  ---
!
      TMCX(1) = PORDX*VOLX*(SLX*XLAX*RHOLX + SGX*XGCX(1)*RHOGX + 
     &    SNX*XNCX(1)*RHONX)
!
!---  Total mass of non-CO2 petroleum components  ---
!
      DO 50 IGC = 2,NGC+2
        TMCX(IGC) = PORDX*VOLX*(SGX*XGCX(IGC)*RHOGX + 
     &    SNX*XNCX(IGC)*RHONX)
   50 CONTINUE
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
      SUBROUTINE FLH_PC3( BETAX,FKX,GNIFTX,PORDX,PVAX,PVWX,PX,
     &  RHOGX,RHOLX,RHOMGX,RHOMLX,RHOMNX,RHONX,
     &  SGX,SLX,SNX,TMCX,TMSX,TX,VOLX,XGCX,
     &  XGWX,XLAX,XLSX,XLWX,XNCX,XMGCX,XMGWX,XMLAX,XMLSX,XMLWX,XMNCX,
     &  YLSX,ZGX,ZMCX,ZNX,IBETAX,IZM,IZN,NPHAZX )
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
!     SL < 1.0, ST > 0.0, SG > 0.0
!
!     energy - temperature, C - T
!     water mass - aqueous saturation, SL
!     CO2 mass
!       if( ZMC(1,2,N) > max( ZMC(J,2,N) ) J /= 1
!         - gas (system) pressure, Pa gauge - PG,PSO
!       else
!         - nonaqueous mole fraction of CO2 - ZMC(1,2,N)
!       endif
!     CH4 mass
!       if( ZMC(2,2,N) > max( ZMC(J,2,N) ) J /= 2
!         - gas (system) pressure, Pa gauge - PG,PSO
!       else
!         - nonaqueous mole fraction of CO2 - ZMC(2,2,N)
!       endif
!     petroleum component mass
!       if( ZMC(I,2,N) > max( ZMC(J,2,N) ) J /= I
!         - gas (system) pressure, Pa gauge - PG,PSO
!       else
!         - nonaqueous mole fraction of CO2 - ZMC(I,2,N)
!       endif
!     salt mass - total salt brine mass fraction - YLS
!
!     Output variables:
!
!     RHOGX - gas density, kg/m^3
!     RHOLX - aqueous density, kg/m^3
!     RHONX - nonaqueous-liquid density, kg/m^3
!     RHOMGX - gas molar density, kmol/m^3
!     RHOMLX - aqueous dmolar density, kmol/m^3
!     RHOMNX - nonaqueous-liquid molar density, kmol/m^3
!     SGX - gas saturation
!     SNX - nonaqueous-liquid saturation
!     XGCX - gas mass fraction of petroleum components
!     XGWX - gas mass fraction of water
!     XLAX - aqueous mass fraction of CO2
!     XLSX - aqueous mass fraction of salt
!     XLWX - aqueous mass fraction of water
!     XNCX - nonaqueous-liquid mass fraction of petroleum components
!     ZMCX - total mole fraction of petroleum components
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GLB_PAR
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 TMCX(LNGC),XGCX(LNGC),XMGCX(LNGC),XMNCX(LNGC)
      REAL*8 XNCX(LNGC),ZMCX(LNGC),FKX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_PC3'
!
!---  Total petroleum molar mass  ---
!
      IF( IZM.NE.0 ) THEN
        ZMCX(IZM) = 1.D+0
        DO 10 IGC = 1,NGC+2
          IF( IGC.NE.IZM ) ZMCX(IZM) = ZMCX(IZM) - ZMCX(IGC)
   10   CONTINUE
        ZMCX(IZM) = MAX( MIN( ZMCX(IZM),1.D+0 ),0.D+0 )
      ENDIF
!
!---  Maximum solubility of salt in pure water  ---
!
      CALL SOL_LS( TX,XBSMX )
!
!---  Salt mass fraction in brine  ---
!
      XBSX = MIN( YLSX,XBSMX )
!
!---  Water vapor pressure in equilibrium with brine at saturated
!     conditions  ---
!
      CALL SP_B( TX,XBSX,PSBX )
!
!---  Water vapor pressure adjusted for aqueous saturation
!     in geologic media  ---
!
      CALL VPL( PSBX,PVWX,SLX,TX,IZN )
!
!---  Total nonaqueous vapor pressure  ---
!
      PVX = PX - PVWX
!
!---  Nonaqueous mole fractions, densities, and saturations  ---
!
      CALL FLH_NA( BETAX,FKX,GNIFTX,PVX,PVWX,RHOMGX,RHOMNX,SGX,SNX,
     &  TX,XMGCX,XMGWX,XMNCX,ZGX,ZMCX,ZNX,IBETAX,NPHAZX )
      SGX = (1.D+0-SLX)*SGX
      SNX = (1.D+0-SLX)*SNX
!
!---  CO2 vapor pressure  ---
!
      PVAX = PVX*XMGCX(1)
!
!---  Solubility of CO2 in aqueous phase  ---
!
      CALL SP_B( TX,XBSX,PSBX )
      PTX = PSBX + PVAX
      CALL EQUIL( TX,PTX,PVAX,XBSX,XLAX,XLSX,XLWX )
!
!---  Brine density  ---
!
      CALL DENS_B( TX,PX,XBSX,RHOBX )
!
!---  Aqueous density  ---
!
      CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
!
!---  Aqueous component mole fractions  ---
!
      WTMLX = XLAX/WTMA + XLSX/WTMS + XLWX/WTMW
      XMLAX = (XLAX/WTMA)/WTMLX
      XMLSX = (XLSX/WTMS)/WTMLX
      XMLWX = (XLWX/WTMW)/WTMLX
!
!---  Aqueous molar density  ---
!
      RHOMLX = RHOLX*WTMLX
!
!---  Nonaqueous mass fractions  ---
!
      WTMGX = XMGWX*WTMW
      WTMNX = 0.D+0
      DO 30 IGC = 1,NGC+2
        WTMGX = WTMGX + XMGCX(IGC)*GCPP(1,IGC)
        WTMNX = WTMNX + XMNCX(IGC)*GCPP(1,IGC)
   30 CONTINUE
      XGWX = XMGWX*WTMW/WTMGX
      DO 40 IGC = 1,NGC+2
        XGCX(IGC) = XMGCX(IGC)*GCPP(1,IGC)/WTMGX
        XNCX(IGC) = XMNCX(IGC)*GCPP(1,IGC)/WTMNX
   40 CONTINUE
!
!---  Noaqueous densities  ---
!
      RHOGX = RHOMGX*WTMGX
      RHONX = RHOMNX*WTMNX
!
!---  Total mass of salt  ---
!
      TMSX = YLSX*RHOBX*PORDX*VOLX
!
!---  Total mass of CO2  ---
!
      TMCX(1) = PORDX*VOLX*(SLX*XLAX*RHOLX + SGX*XGCX(1)*RHOGX + 
     &    SNX*XNCX(1)*RHONX)
!
!---  Total mass of non-CO2 petroleum components  ---
!
      DO 50 IGC = 2,NGC+2
        TMCX(IGC) = PORDX*VOLX*(SGX*XGCX(IGC)*RHOGX + 
     &    SNX*XNCX(IGC)*RHONX)
   50 CONTINUE
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
      SUBROUTINE FLH_PC4( BETAX,FKX,GNIFTX,PORDX,PVAX,PVWX,PX,
     &  RHOGX,RHOLX,RHOMGX,RHOMLX,RHOMNX,RHONX,
     &  SGX,SLX,SNX,TMCX,TMSX,TX,VOLX,XGCX,
     &  XGWX,XLAX,XLSX,XLWX,XNCX,XMGCX,XMGWX,XMLAX,XMLSX,XMLWX,XMNCX,
     &  YLSX,ZGX,ZMCX,ZNX,IBETAX,IZM,IZN,NPHAZX )
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
!     SL = 0.0, ST = 1.0, SG = 0.0
!
!     energy - temperature, C - T
!     water mass - aqueous saturation, SL
!     CO2 mass
!       if( ZMC(1,2,N) > max( ZMC(J,2,N) ) J /= 1
!         - nonaqueous (system) pressure, Pa gauge - PN,PSO
!       else
!         - nonaqueous mole fraction of CO2 - ZMC(1,2,N)
!       endif
!     CH4 mass
!       if( ZMC(2,2,N) > max( ZMC(J,2,N) ) J /= 2
!         - nonaqueous (system) pressure, Pa gauge - PN,PSO
!       else
!         - nonaqueous mole fraction of CO2 - ZMC(2,2,N)
!       endif
!     petroleum component mass
!       if( ZMC(I,2,N) > max( ZMC(J,2,N) ) J /= I
!         - nonaqueous (system) pressure, Pa gauge - PN,PSO
!       else
!         - nonaqueous mole fraction of CO2 - ZMC(I,2,N)
!       endif
!     salt mass - total salt mass - TMS
!
!     Output variables:
!
!     RHOGX - gas density, kg/m^3
!     RHOLX - aqueous density, kg/m^3
!     RHONX - nonaqueous-liquid density, kg/m^3
!     RHOMGX - gas molar density, kmol/m^3
!     RHOMLX - aqueous dmolar density, kmol/m^3
!     RHOMNX - nonaqueous-liquid molar density, kmol/m^3
!     SGX - gas saturation
!     SNX - nonaqueous-liquid saturation
!     XGCX - gas mass fraction of petroleum components
!     XGWX - gas mass fraction of water
!     XLAX - aqueous mass fraction of CO2
!     XLSX - aqueous mass fraction of salt
!     XLWX - aqueous mass fraction of water
!     XNCX - nonaqueous-liquid mass fraction of petroleum components
!     ZMCX - total mole fraction of petroleum components
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GLB_PAR
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 TMCX(LNGC),XGCX(LNGC),XMGCX(LNGC),XMNCX(LNGC)
      REAL*8 XNCX(LNGC),ZMCX(LNGC),FKX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_PC4'
!
!---  Total petroleum molar mass  ---
!
      IF( IZM.NE.0 ) THEN
        ZMCX(IZM) = 1.D+0
        DO 10 IGC = 1,NGC+2
          IF( IGC.NE.IZM ) ZMCX(IZM) = ZMCX(IZM) - ZMCX(IGC)
   10   CONTINUE
        ZMCX(IZM) = MAX( MIN( ZMCX(IZM),1.D+0 ),0.D+0 )
      ENDIF
!
!---  Maximum solubility of salt in pure water  ---
!
      CALL SOL_LS( TX,XBSMX )
      IF( TMSX.GT.0.D+0 ) THEN
        YLSX = XBSMX
      ELSE
        YSLX = 0.D+0
      ENDIF
!
!---  Salt mass fraction in brine  ---
!
      XBSX = MIN( YLSX,XBSMX )
!
!---  Water vapor pressure in equilibrium with brine at saturated
!     conditions  ---
!
      CALL SP_B( TX,XBSX,PSBX )
!
!---  Water vapor pressure adjusted for aqueous saturation
!     in geologic media  ---
!
      CALL VPL( PSBX,PVWX,SLX,TX,IZN )
!
!---  Total nonaqueous vapor pressure  ---
!
      PVX = PX - PVWX
!
!---  Nonaqueous mole fractions, densities, and saturations  ---
!
      CALL FLH_NA( BETAX,FKX,GNIFTX,PVX,PVWX,RHOMGX,RHOMNX,SGX,SNX,
     &  TX,XMGCX,XMGWX,XMNCX,ZGX,ZMCX,ZNX,IBETAX,NPHAZX )
      SGX = (1.D+0-SLX)*SGX
      SNX = (1.D+0-SLX)*SNX
!
!---  CO2 vapor pressure  ---
!
      PVAX = PVX*XMGCX(1)
!
!---  Solubility of CO2 in aqueous phase  ---
!
      CALL SP_B( TX,XBSX,PSBX )
      PTX = PSBX + PVAX
      CALL EQUIL( TX,PTX,PVAX,XBSX,XLAX,XLSX,XLWX )
!
!---  Brine density  ---
!
      CALL DENS_B( TX,PX,XBSX,RHOBX )
!
!---  Aqueous density  ---
!
      CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
!
!---  Aqueous component mole fractions  ---
!
      WTMLX = XLAX/WTMA + XLSX/WTMS + XLWX/WTMW
      XMLAX = (XLAX/WTMA)/WTMLX
      XMLSX = (XLSX/WTMS)/WTMLX
      XMLWX = (XLWX/WTMW)/WTMLX
!
!---  Aqueous molar density  ---
!
      RHOMLX = RHOLX*WTMLX
!
!---  Nonaqueous mass fractions  ---
!
      WTMGX = XMGWX*WTMW
      WTMNX = 0.D+0
      DO 30 IGC = 1,NGC+2
        WTMGX = WTMGX + XMGCX(IGC)*GCPP(1,IGC)
        WTMNX = WTMNX + XMNCX(IGC)*GCPP(1,IGC)
   30 CONTINUE
      XGWX = XMGWX*WTMW/WTMGX
      DO 40 IGC = 1,NGC+2
        XGCX(IGC) = XMGCX(IGC)*GCPP(1,IGC)/WTMGX
        XNCX(IGC) = XMNCX(IGC)*GCPP(1,IGC)/WTMNX
   40 CONTINUE
!
!---  Noaqueous densities  ---
!
      RHOGX = RHOMGX*WTMGX
      RHONX = RHOMNX*WTMNX
!
!---  Total mass of CO2  ---
!
      TMCX(1) = PORDX*VOLX*(SLX*XLAX*RHOLX + SGX*XGCX(1)*RHOGX + 
     &    SNX*XNCX(1)*RHONX)
!
!---  Total mass of non-CO2 petroleum components  ---
!
      DO 50 IGC = 2,NGC+2
        TMCX(IGC) = PORDX*VOLX*(SGX*XGCX(IGC)*RHOGX + 
     &    SNX*XNCX(IGC)*RHONX)
   50 CONTINUE
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
      SUBROUTINE FLH_PC5( BETAX,FKX,GNIFTX,PORDX,PVAX,PVWX,PX,
     &  RHOGX,RHOLX,RHOMGX,RHOMLX,RHOMNX,RHONX,
     &  SGX,SLX,SNX,TMCX,TMSX,TX,VOLX,XGCX,
     &  XGWX,XLAX,XLSX,XLWX,XNCX,XMGCX,XMGWX,XMLAX,XMLSX,XMLWX,XMNCX,
     &  YLSX,ZGX,ZMCX,ZNX,IBETAX,IZM,IZN,NPHAZX )
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
!     SL = 0.0, ST = 1.0, SG > 0.0
!
!     energy - temperature, C - T
!     water mass - aqueous saturation, SL
!     CO2 mass
!       if( ZMC(1,2,N) > max( ZMC(J,2,N) ) J /= 1
!         - gas (system) pressure, Pa gauge - PG,PSO
!       else
!         - nonaqueous mole fraction of CO2 - ZMC(1,2,N)
!       endif
!     CH4 mass
!       if( ZMC(2,2,N) > max( ZMC(J,2,N) ) J /= 2
!         - gas (system) pressure, Pa gauge - PG,PSO
!       else
!         - nonaqueous mole fraction of CO2 - ZMC(2,2,N)
!       endif
!     petroleum component mass
!       if( ZMC(I,2,N) > max( ZMC(J,2,N) ) J /= I
!         - gas (system) pressure, Pa gauge - PG,PSO
!       else
!         - nonaqueous mole fraction of CO2 - ZMC(I,2,N)
!       endif
!     salt mass - total salt mass - TMS
!
!     Output variables:
!
!     RHOGX - gas density, kg/m^3
!     RHOLX - aqueous density, kg/m^3
!     RHONX - nonaqueous-liquid density, kg/m^3
!     RHOMGX - gas molar density, kmol/m^3
!     RHOMLX - aqueous dmolar density, kmol/m^3
!     RHOMNX - nonaqueous-liquid molar density, kmol/m^3
!     SGX - gas saturation
!     SNX - nonaqueous-liquid saturation
!     XGCX - gas mass fraction of petroleum components
!     XGWX - gas mass fraction of water
!     XLAX - aqueous mass fraction of CO2
!     XLSX - aqueous mass fraction of salt
!     XLWX - aqueous mass fraction of water
!     XNCX - nonaqueous-liquid mass fraction of petroleum components
!     ZMCX - total mole fraction of petroleum components
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GLB_PAR
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 TMCX(LNGC),XGCX(LNGC),XMGCX(LNGC),XMNCX(LNGC)
      REAL*8 XNCX(LNGC),ZMCX(LNGC),FKX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_PC5'
!
!---  Total petroleum molar mass  ---
!
      IF( IZM.NE.0 ) THEN
        ZMCX(IZM) = 1.D+0
        DO 10 IGC = 1,NGC+2
          IF( IGC.NE.IZM ) ZMCX(IZM) = ZMCX(IZM) - ZMCX(IGC)
   10   CONTINUE
        ZMCX(IZM) = MAX( MIN( ZMCX(IZM),1.D+0 ),0.D+0 )
      ENDIF
!
!---  Maximum solubility of salt in pure water  ---
!
      CALL SOL_LS( TX,XBSMX )
      IF( TMSX.GT.0.D+0 ) THEN
        YLSX = XBSMX
      ELSE
        YSLX = 0.D+0
      ENDIF
!
!---  Salt mass fraction in brine  ---
!
      XBSX = MIN( YLSX,XBSMX )
!
!---  Water vapor pressure in equilibrium with brine at saturated
!     conditions  ---
!
      CALL SP_B( TX,XBSX,PSBX )
!
!---  Water vapor pressure adjusted for aqueous saturation
!     in geologic media  ---
!
      CALL VPL( PSBX,PVWX,SLX,TX,IZN )
!
!---  Total nonaqueous vapor pressure  ---
!
      PVX = PX - PVWX
!
!---  Nonaqueous mole fractions, densities, and saturations  ---
!
      CALL FLH_NA( BETAX,FKX,GNIFTX,PVX,PVWX,RHOMGX,RHOMNX,SGX,SNX,
     &  TX,XMGCX,XMGWX,XMNCX,ZGX,ZMCX,ZNX,IBETAX,NPHAZX )
      SGX = (1.D+0-SLX)*SGX
      SNX = (1.D+0-SLX)*SNX
!
!---  CO2 vapor pressure  ---
!
      PVAX = PVX*XMGCX(1)
!
!---  Solubility of CO2 in aqueous phase  ---
!
      CALL SP_B( TX,XBSX,PSBX )
      PTX = PSBX + PVAX
      CALL EQUIL( TX,PTX,PVAX,XBSX,XLAX,XLSX,XLWX )
!
!---  Brine density  ---
!
      CALL DENS_B( TX,PX,XBSX,RHOBX )
!
!---  Aqueous density  ---
!
      CALL DENS_L( TX,RHOBX,XLAX,RHOLX )
!
!---  Aqueous component mole fractions  ---
!
      WTMLX = XLAX/WTMA + XLSX/WTMS + XLWX/WTMW
      XMLAX = (XLAX/WTMA)/WTMLX
      XMLSX = (XLSX/WTMS)/WTMLX
      XMLWX = (XLWX/WTMW)/WTMLX
!
!---  Aqueous molar density  ---
!
      RHOMLX = RHOLX*WTMLX
!
!---  Nonaqueous mass fractions  ---
!
      WTMGX = XMGWX*WTMW
      WTMNX = 0.D+0
      DO 30 IGC = 1,NGC+2
        WTMGX = WTMGX + XMGCX(IGC)*GCPP(1,IGC)
        WTMNX = WTMNX + XMNCX(IGC)*GCPP(1,IGC)
   30 CONTINUE
      XGWX = XMGWX*WTMW/WTMGX
      DO 40 IGC = 1,NGC+2
        XGCX(IGC) = XMGCX(IGC)*GCPP(1,IGC)/WTMGX
        XNCX(IGC) = XMNCX(IGC)*GCPP(1,IGC)/WTMNX
   40 CONTINUE
!
!---  Noaqueous densities  ---
!
      RHOGX = RHOMGX*WTMGX
      RHONX = RHOMNX*WTMNX
!
!---  Total mass of CO2  ---
!
      TMCX(1) = PORDX*VOLX*(SLX*XLAX*RHOLX + SGX*XGCX(1)*RHOGX + 
     &    SNX*XNCX(1)*RHONX)
!
!---  Total mass of non-CO2 petroleum components  ---
!
      DO 50 IGC = 2,NGC+2
        TMCX(IGC) = PORDX*VOLX*(SGX*XGCX(IGC)*RHOGX + 
     &    SNX*XNCX(IGC)*RHONX)
   50 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_PC5 group  ---
!
      RETURN
      END



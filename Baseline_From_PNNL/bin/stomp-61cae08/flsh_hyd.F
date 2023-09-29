!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_305( PGX,PHX,PLX,PNX,PVAX,PVOX,TX,
     &  RHOHX,RHOLX,SHX,SIX,SLX,YHOX,YLSX,YMGOX,XLAX,XLOX,XHAX,
     &  XHOX,XHWX,XMLAX,XMLOX,INDX,IZN,N )
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
!     Givens:
!
!     Hydrate saturation
!     Hydrate CH4 mass fraction of formers
!     Aqueous salt mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 October 2007.
!     Last Modified by MD White, PNNL, 15 October 2007.
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
      PARAMETER(KP=6)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 GX(KP,KP-1),AJ(KP,KP),BJ(KP),RKLX(3)
      INTEGER IJ(KP)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_305'
      INDX = 0
!
!---  Dissolved-salt aqueous mass fraction  ---
!
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
      CALL SP_B( TX,XLSX,PSBX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      CALL SFT_L( TX,XLSX,SIG_GL )
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
!---  Guess gas CH4 mole fraction of formers  ---
!
      YHAX = 1.D+0-YHOX
      YMGOX = (YHOX/WTMO)/((YHOX/WTMO)+(YHAX/WTMA))
      SHZ = MAX( SHX,1.D-9 )
      NC = 0
!
!---  Entry pressure with hydrate barrier adjustment  ---
!
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
!---  Single-variable Newton-Raphson loop (YMGO)  ---
!
   10 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 17
        CHMSG = 'Unconverged FLH_305: Dissolved-Salt  ' //
     &    'Aqueous Mass Fraction and Total CH4 Mass ' //
     &    'Fraction of Hydrate Formers: '
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      DO 20 M = 1,2
        YMGOY = YMGOX
        DYMGOX = SIGN( 1.D-6,(1.D-1-YMGOX) )
        IF( M.EQ.2 ) YMGOY = YMGOX + DYMGOX
        CALL HYD_P( TEQHX,YMGOY,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
        GX(M,1) = YHOX - XHOX/(XHOX+XHAX)
   20 CONTINUE
!
!---  Load linear system  ---
!
      AJ(1,1) = (GX(2,1)-GX(1,1))/DYMGOX
      BJ(1) = -GX(1,1)
!
!---  Solve linear system  ---
!
      JP = 1
      CALL LUDCMP( AJ,JP,KP,IJ,DJ )
      CALL LUBKSB( AJ,JP,KP,IJ,BJ )
      YMGOX = MIN( MAX( YMGOX + BJ(1),0.D+0 ),1.D+0 )
      IF( ABS(BJ(1)).GT.1.D-7 ) GOTO 10
!
!---  Converged solution  ---
!
      CALL HYD_P( TEQHX,YMGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
      PVX = PEQHX
      PVOX = PVX*YMGOX
      PVAX = PVX-PVOX
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
      PPLX = MAX( PLX+PATM,PSBX )
      CALL DENS_B( TX,PPLX,XLSX,RHOBX )
      CALL DENS_L( TX,RHOBX,XLAX,XLOX,RHOLX )
      IF( ISLC(44).EQ.1 ) THEN
        PVWX = PSBX
      ELSE
        CALL VPL_BX( TX,PLX,PVAX,PVOX,PVWX,PSBX,RHOLX,XLSX )
      ENDIF
      PGX = MAX( PVOX+PVAX+PVWX-PATM,PLX+(ENPR/BGL) )
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
      INDX = 0
      CALL KSP37( PGX,PIX,PLX,PNC,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
!
!---  Check for non-zero liquid-CO2 saturation   ---
!
      IF( PNX.GT.PNC ) THEN
        INDX = 17
        RLMSG = PNX - PNC
        CHMSG = 'FLH_305: Non-zero Liquid-CO2 Saturation: PN - PNC = '
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_305 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_350( PGX,PHX,PLX,PNX,PVAX,PVOX,TX,
     &  RHOHX,RHOLX,SHX,SIX,SLX,YHAX,YLSX,YMGOX,XLAX,XLOX,XHAX,
     &  XHOX,XHWX,XMLAX,XMLOX,INDX,IZN,N )
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
!     Givens:
!
!     Hydrate saturation
!     Hydrate CO2 mass fraction of formers
!     Aqueous salt mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 October 2007.
!     Last Modified by MD White, PNNL, 15 October 2007.
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
      PARAMETER(KP=6)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 GX(KP,KP-1),AJ(KP,KP),BJ(KP),RKLX(3)
      INTEGER IJ(KP)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_350'
      INDX = 0
!
!---  Dissolved-salt aqueous mass fraction  ---
!
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
      CALL SP_B( TX,XLSX,PSBX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      CALL SFT_L( TX,XLSX,SIG_GL )
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
!---  Guess gas CH4 mole fraction of formers  ---
!
      YHOX = 1.D+0-YHAX
      YMGOX = (YHOX/WTMO)/((YHOX/WTMO)+(YHAX/WTMA))
      NC = 0
!
!---  Pore-space equilibrium temperature  ---
!
      TEQHX = TX
!
!---  Single-variable Newton-Raphson loop (YMGO)  ---
!
   10 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 17
        CHMSG = 'Unconverged FLH_350: Dissolved-Salt  ' //
     &    'Aqueous Mass Fraction and Total CH4 Mass ' //
     &    'Fraction of Hydrate Formers: '
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      DO 20 M = 1,2
        YMGOY = YMGOX
        DYMGOX = SIGN( 1.D-6,(1.D-1-YMGOX) )
        IF( M.EQ.2 ) YMGOY = YMGOX + DYMGOX
        CALL HYD_P( TEQHX,YMGOY,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
        PVX = PEQHX
        PVOX = PVX*YMGOY
        PVAX = PVX-PVOX
        INCG = 1
        CALL SOL_LA( TX,PVAX,XLSX,XMLAX )
        INCG = 2
        CALL SOL_LO( TX,PVOX,XLSX,XMLOX )
        XMLSX = XLSX*(XMLOX*WTMO + XMLAX*WTMA + 
     &    (1.D+0-XMLOX-XMLAX)*WTMW)/((1.D+0-XLSX)*WTMS + XLSX*WTMW)
        XMLWX = MAX( 1.D+0-XMLSX-XMLAX-XMLOX,0.D+0 )
        WTMLX = XMLSX*WTMS + XMLAX*WTMA + XMLOX*WTMO + XMLWX*WTMW
        XLAX = XMLAX*WTMA/WTMLX
        XLOX = XMLOX*WTMO/WTMLX
        PPLX = MAX( PLX+PATM,PSBX )
        CALL DENS_B( TX,PPLX,XLSX,RHOBX )
        CALL DENS_L( TX,RHOBX,XLAX,XLOX,RHOLX )
        IF( ISLC(44).EQ.1 ) THEN
          PVWX = PSBX
        ELSE
          CALL VPL_BX( TX,PLX,PVAX,PVOX,PVWX,PSBX,RHOLX,XLSX )
        ENDIF
        PGX = MAX( PVOX+PVAX+PSBX-PATM,PLX+(ENPR/BGL) )
        PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
        PCX = MAX( PGX-PLX,0.D+0 )
        CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
        INDX = 0
        CALL KSP37( PGX,PIX,PLX,PNC,SGX,SHX,SIX,SLX,SNX,
     &    RKGX,RKLX,RKNX,TX,INDX,IZN,N )
        IF( ISLC(44).EQ.1 ) THEN
          PVWX = PSBX
        ELSE
          CALL VPL_B( TX,PSBX,PCX,RHOLX,PVWX,XLSX )
        ENDIF
        ISRX = 2
        CALL DENS_W( TX,PVWX,RHOX,RHOGWX,ISRX )
        CALL DENS_A( TX,PVAX,RHOGAX,ISRX )
        CALL DENS_O( TX,PVOX,RHOGOX )
        RHOGX = RHOGWX+RHOGAX+RHOGOX
        XGAX = RHOGAX/RHOGX
        XGOX = RHOGOX/RHOGX
        XGWX = RHOGWX/RHOGX
        TMAX = XHAX*RHOHX*SHX
        TMOX = XHOX*RHOHX*SHX
        GX(M,1) = YHAX - TMAX/(TMOX+TMAX)
   20 CONTINUE
!
!---  Load linear system  ---
!
      AJ(1,1) = (GX(2,1)-GX(1,1))/DYMGOX
      BJ(1) = -GX(1,1)
!
!---  Solve linear system  ---
!
      JP = 1
      CALL LUDCMP( AJ,JP,KP,IJ,DJ )
      CALL LUBKSB( AJ,JP,KP,IJ,BJ )
      YMGOX = MIN( MAX( YMGOX + BJ(1),0.D+0 ),1.D+0 )
      IF( ABS(BJ(1)).GT.1.D-7 ) GOTO 10
!
!---  Converged solution  ---
!
      CALL HYD_P( TEQHX,YMGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
      PVX = PEQHX
      PVOX = PVX*YMGOX
      PVAX = PVX-PVOX
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
      PPLX = MAX( PLX+PATM,PSBX )
      CALL DENS_B( TX,PPLX,XLSX,RHOBX )
      CALL DENS_L( TX,RHOBX,XLAX,XLOX,RHOLX )
      IF( ISLC(44).EQ.1 ) THEN
        PVWX = PSBX
      ELSE
        CALL VPL_BX( TX,PLX,PVAX,PVOX,PVWX,PSBX,RHOLX,XLSX )
      ENDIF
      PGX = MAX( PVOX+PVAX+PSBX-PATM,PLX+(ENPR/BGL) )
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
      INDX = 0
      CALL KSP37( PGX,PIX,PLX,PNC,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
!
!---  Check for non-zero liquid-CO2 saturation   ---
!
      IF( PNX.GT.PNC ) THEN
        INDX = 17
        RLMSG = PNX - PNC
        CHMSG = 'FLH_350: Non-zero Liquid-CO2 Saturation: PN - PNC = '
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_350 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_378( PGX,PHX,PLX,PNX,PVAX,PVOX,PVWX,TX,RHOHX,
     &  RHOLX,SHX,SIX,SLX,YLSX,YMGAX,YMGOX,XLAX,XLOX,XHAX,
     &  XHOX,XHWX,XMLAX,XMLOX,INDX,IZN,N )
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
!     Givens:
!
!     Gas pressure
!     Gas CO2 mass fraction of formers
!     Aqueous salt mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 October 2007.
!     Last Modified by MD White, PNNL, 15 October 2007.
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
      PARAMETER(KP=6)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_378'
      INDX = 0
!
!---  Dissolved-salt aqueous mass fraction  ---
!
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
      CALL SP_B( TX,XLSX,PSBX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      CALL SFT_L( TX,XLSX,SIG_GL )
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
      SHX = 0.D+0
      HCPBFX = 1.D+0
      IF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.12 ) THEN
        HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &    SIN(SCHR(12,IZN)*SHX + SCHR(13,IZN)) )
      ENDIF
      ENPR = SCHR(2,IZN)*RHORL*GRAV*HCPBFX
!
!---  Water vapor pressure ---
!
      PX = MAX( PGX,PLX )+PATM
      PGX = MAX( PGX,PLX+(ENPR/BGL) )
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      ISRX = 1
      CALL DENS_W( TX,PX,RHOLWX,RHOGWX,ISRX )
      CALL DENS_B( TX,PX,XLSX,RHOBX )
      PCX = MAX( PGX-PLX,0.D+0 )
      IF( ISLC(44).EQ.1 ) THEN
        PVWX = PSBX
      ELSE
        CALL VPL_B( TX,PSBX,PCX,RHOBX,PVWX,XLSX )
      ENDIF
      PVX = MAX( PGX+PATM-PVWX,0.D+0 )
      YMGOX = 1.D+0-YMGAX
      PVOX = MAX( MIN( YMGOX,1.D+0 ),0.D+0)*PVX
      PVAX = PVX-PVOX
      CALL SOL_LO( TX,PVOX,XLSX,XMLOX )
      CALL SOL_LA( TX,PVAX,XLSX,XMLAX )
      XMLSX = XLSX*(XMLAX*WTMA + XMLOX*WTMO + 
     &  (1.D+0-XMLAX-XMLOX)*WTMW)/((1.D+0-XLSX)*WTMS + XLSX*WTMW)
      XMLWX = MAX( 1.D+0-XMLSX-XMLOX-XMLAX,0.D+0 )
      WTMLX = XMLSX*WTMS + XMLOX*WTMO + XMLAX*WTMA + XMLWX*WTMW
      XLAX = XMLAX*WTMA/WTMLX
      XLOX = XMLOX*WTMO/WTMLX
      CALL DENS_B( TX,PX,XLSX,RHOBX )
      CALL DENS_L( TX,RHOBX,XLAX,XLOX,RHOLX )
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
      INDX = 0
      CALL KSP37( PGX,PIX,PLX,PNC,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
!
!---  Check for hydrate conditions   ---
!
      IF( (PVAX+PVOX).GT.EPSL ) THEN
        YMGOX = PVOX/(PVAX+PVOX)
        PEQHX = PVOX+PVAX
!
!---    Inhibitor correction  ---
!
        IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
          CALL HYD_PO( YMGOX,POHX )
        ENDIF
        CALL HYD_T( TEQHX,YMGOX,XLSX,PEQHX,POHX,RHOHX,XHAX,XHOX,XHWX )
        IF( TX.LE.TEQHX ) THEN
          INDX = 17
          RLMSG = TEQHX
          CHMSG = 'FLH_378: Non-Zero Hydrate Saturation: ' //
     &       'Temperature < Hydrate Equilibrium Temperature = '
          ISUB_LOG = ISUB_LOG-1
         RETURN
       ENDIF
      ENDIF
!
!---  Check for non-zero liquid-CO2 saturation   ---
!
      IF( PNX.GT.PNC ) THEN
        INDX = 17
        RLMSG = PNX - PNC
        CHMSG = 'FLH_378: Non-zero Liquid-CO2 Saturation: PN - PNC = '
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_378 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_387( PGX,PHX,PLX,PNX,PVAX,PVOX,PVWX,TX,RHOHX,
     &  RHOLX,SHX,SIX,SLX,YLSX,YMGOX,XLAX,XLOX,XHAX,XHOX,XHWX,XMLAX,
     &  XMLOX,INDX,IZN,N )
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
!     Givens:
!
!     Gas pressure
!     Gas CH4 mass fraction of formers
!     Aqueous salt mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 October 2007.
!     Last Modified by MD White, PNNL, 15 October 2007.
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
      PARAMETER(KP=6)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RKLX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_387'
      INDX = 0
!
!---  Dissolved-salt aqueous mass fraction  ---
!
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
      CALL SP_B( TX,XLSX,PSBX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      CALL SFT_L( TX,XLSX,SIG_GL )
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
      SHX = 0.D+0
      HCPBFX = 1.D+0
      IF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.12 ) THEN
        HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &    SIN(SCHR(12,IZN)*SHX + SCHR(13,IZN)) )
      ENDIF
      ENPR = SCHR(2,IZN)*RHORL*GRAV*HCPBFX
!
!---  Water vapor pressure  ---
!
      PGX = MAX( PGX,PLX+(ENPR/BGL) )
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      PX = MAX( PGX,PLX )+PATM
      CALL DENS_B( TX,PX,XLSX,RHOBX )
      PCX = MAX( PGX-PLX,0.D+0 )
      IF( ISLC(44).EQ.1 ) THEN
        PVWX = PSBX
      ELSE
        CALL VPL_B( TX,PSBX,PCX,RHOBX,PVWX,XLSX )
      ENDIF
      PVX = MAX( PGX+PATM-PVWX,0.D+0 )
      PVOX = MAX( MIN( YMGOX,1.D+0 ),0.D+0)*PVX
      PVAX = PVX-PVOX
      CALL SOL_LO( TX,PVOX,XLSX,XMLOX )
      CALL SOL_LA( TX,PVAX,XLSX,XMLAX )
      XMLSX = XLSX*(XMLAX*WTMA + XMLOX*WTMO + 
     &  (1.D+0-XMLAX-XMLOX)*WTMW)/((1.D+0-XLSX)*WTMS + XLSX*WTMW)
      XMLWX = MAX( 1.D+0-XMLSX-XMLOX-XMLAX,0.D+0 )
      WTMLX = XMLSX*WTMS + XMLOX*WTMO + XMLAX*WTMA + XMLWX*WTMW
      XLAX = XMLAX*WTMA/WTMLX
      XLOX = XMLOX*WTMO/WTMLX
      CALL DENS_B( TX,PX,XLSX,RHOBX )
      CALL DENS_L( TX,RHOBX,XLAX,XLOX,RHOLX )
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
      INDX = 0
      CALL KSP37( PGX,PIX,PLX,PNC,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
!
!---  Check for hydrate conditions   ---
!
      IF( (PVAX+PVOX).GT.EPSL ) THEN
        YMGOX = PVOX/(PVAX+PVOX)
        PEQHX = PVOX+PVAX
!
!---    Inhibitor correction  ---
!
        IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
          CALL HYD_PO( YMGOX,POHX )
        ENDIF
        CALL HYD_T( TEQHX,YMGOX,XLSX,PEQHX,POHX,RHOHX,XHAX,XHOX,XHWX )
        IF( TX.LE.TEQHX ) THEN
          INDX = 17
          RLMSG = TEQHX
          CHMSG = 'FLH_387: Non-Zero Hydrate Saturation: ' //
     &       'Temperature < Hydrate Equilibrium Temperature = '
          ISUB_LOG = ISUB_LOG-1
         RETURN
       ENDIF
      ENDIF
!
!---  Check for non-zero liquid-CO2 saturation   ---
!
      IF( PNX.GT.PNC ) THEN
        INDX = 17
        RLMSG = PNX - PNC
        CHMSG = 'FLH_387: Non-zero Liquid-CO2 Saturation: PN - PNC = '
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_387 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_390( PGX,PHX,PLX,PNX,PVAX,PVOX,TX,
     &  RHOHX,RHOLX,SHX,SIX,SLX,SNX,YHOX,YLSX,YMGOX,XLAX,XLOX,XHAX,
     &  XHOX,XHWX,XMLAX,XMLOX,INDX,IZN,N )
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
!     Givens:
!
!     Hydrate saturation
!     Hydrate CH4 mass fraction of formers
!     Aqueous salt mass fraction
!     Liquid CO2 pressure > aqueous pressure
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 October 2007.
!     Last Modified by MD White, PNNL, 15 October 2007.
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
      PARAMETER(KP=6)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 GX(KP,KP-1),AJ(KP,KP),BJ(KP),RKLX(3)
      INTEGER IJ(KP)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_390'
      INDX = 0
!
!---  Dissolved-salt aqueous mass fraction  ---
!
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
      CALL SP_B( TX,XLSX,PSBX )
      CALL SP_A( TX,PSAX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      CALL SFT_L( TX,XLSX,SIG_GL )
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
!---  Guess gas CH4 mole fraction of formers  ---
!
      YMGOX = 5.D-1
      SHZ = MAX( SHX,1.D-9 )
      NC = 0
!
!---  Entry pressure with hydrate barrier adjustment  ---
!
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
!---  Assume liquid-CO2 conditions  ---
!
   10 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 10
        CHMSG = 'Unconverged FLH_390: Dissolved-Salt ' //
     &    'Aqueous Mass Fraction and Liquid-CO2 Pressure: '
        RLMSG = YMGOX
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      DO 20 M = 1,2
        YMGOY = YMGOX
        DYMGOY = SIGN( 1.D-6,(1.D-1-YMGOX) )
        IF( M.EQ.2 ) YMGOY = YMGOX + DYMGOY
        CALL HYD_P( TEQHX,YMGOY,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
        PVX = PEQHX
        PVOX = PVX*YMGOY
        PVAX = PVX-PVOX
        GX(M,1) = PSAX - PVAX
   20 CONTINUE
!
!---  Load linear system  ---
!
      AJ(1,1) = (GX(2,1)-GX(1,1))/DYMGOY
      BJ(1) = -GX(1,1)
!
!---  Solve linear system  ---
!
      JP = 1
      CALL LUDCMP( AJ,JP,KP,IJ,DJ )
      CALL LUBKSB( AJ,JP,KP,IJ,BJ )
      DYMGOX = SIGN( MIN(ABS(BJ(1)),1.D-1),BJ(1) )
      YMGOX = MIN( MAX( YMGOX+DYMGOX,0.D+0 ),1.D+0 )
      IF( ABS(BJ(1)).GT.1.D-7 ) GOTO 10
!
!---  Converged solution  ---
!
      CALL HYD_P( TEQHX,YMGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
      PVX = PEQHX
      PVOX = PVX*YMGOX
      PVAX = PVX-PVOX
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
      PPLX = MAX( PLX+PATM,PSBX )
      CALL DENS_B( TX,PPLX,XLSX,RHOBX )
      CALL DENS_L( TX,RHOBX,XLAX,XLOX,RHOLX )
      IF( ISLC(44).EQ.1 ) THEN
        PVWX = PSBX
      ELSE
        CALL VPL_BX( TX,PLX,PVAX,PVOX,PVWX,PSBX,RHOLX,XLSX )
      ENDIF
      PGX = MAX( PVOX+PVAX+PVWX-PATM,PLX+(ENPR/BGL) )
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
      INDX = 0
      CALL KSP37( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
!
!---  Check for no liquid-CO2 conditions  ---
!
      IF( PNX.LE.PNC ) THEN
        INDX = 10
        RLMSG = PNX - PNC
        CHMSG = 'FLH_390: No Liquid-CO2: PN - PNC = '
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_390 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_395( PGX,PHX,PLX,PNX,PVAX,PVOX,TX,
     &  RHOHX,RHOLX,SHX,SIX,SLX,SNX,YHOX,YLSX,YMGOX,XLAX,XLOX,XHAX,
     &  XHOX,XHWX,XMLAX,XMLOX,INDX,IZN,N )
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
!     Givens:
!
!     Hydrate saturation
!     Hydrate CH4 mass fraction of formers
!     Aqueous salt mass fraction
!     Liquid CO2 pressure > aqueous pressure
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 October 2007.
!     Last Modified by MD White, PNNL, 15 October 2007.
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
      PARAMETER(KP=6)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 GX(KP,KP-1),AJ(KP,KP),BJ(KP),RKLX(3)
      INTEGER IJ(KP)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_395'
      INDX = 0
!
!---  Dissolved-salt aqueous mass fraction  ---
!
      CALL SOL_LS( TX,XLSMX )
      XLSX = MIN( YLSX,XLSMX )
      CALL SP_B( TX,XLSX,PSBX )
      CALL SP_A( TX,PSAX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      CALL SFT_L( TX,XLSX,SIG_GL )
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
!---  Guess gas CH4 mole fraction of formers  ---
!
      YHAX = 1.D+0-YHOX
      YMGOX = (YHOX/WTMO)/((YHOX/WTMO)+(YHAX/WTMA))
      SHZ = MAX( SHX,1.D-9 )
      NC = 0
!
!---  Entry pressure with hydrate barrier adjustment  ---
!
      HCPBFX = 1.D+0
      IF( ISCHR(IZN).EQ.11 .OR. ISCHR(IZN).EQ.12 ) THEN
        HCPBFX = MAX( 1.D+0,SCHR(10,IZN) + SCHR(11,IZN)*
     &    SIN(SCHR(12,IZN)*SHZ + SCHR(13,IZN)) )
      ENDIF
      ENPR = SCHR(2,IZN)*RHORL*GRAV*HCPBFX
!
!---  Assume liquid-CO2 conditions  ---
!
   10 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 10
        CHMSG = 'Unconverged FLH_395: Dissolved-Salt ' //
     &    'Aqueous Mass Fraction + Liquid-CO2 Pressure ' //
     &    '+ CH4 Hydrate Mass Fraction of Formers: '
        RLMSG = YMGOX
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      DO 20 M = 1,2
        YMGOY = YMGOX
        DYMGOY = SIGN( 1.D-6,(1.D-1-YMGOX) )
        IF( M.EQ.2 ) YMGOY = YMGOX + DYMGOY
        PEQHX = PSAX/(1.D+0-YMGOY)
!
!---    Inhibitor correction  ---
!
        IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
          CALL HYD_PO( YMGOY,POHX )
        ENDIF
        CALL HYD_T( TEQHX,YMGOY,XLSX,PEQHX,POHX,RHOHX,XHAX,XHOX,XHWX )
        GX(M,1) = YHOX - XHOX/(XHOX+XHAX)
   20 CONTINUE
!
!---  Load linear system  ---
!
      AJ(1,1) = (GX(2,1)-GX(1,1))/DYMGOY
      BJ(1) = -GX(1,1)
!
!---  Solve linear system  ---
!
      JP = 1
      CALL LUDCMP( AJ,JP,KP,IJ,DJ )
      CALL LUBKSB( AJ,JP,KP,IJ,BJ )
      DYMGOX = SIGN( MIN(ABS(BJ(1)),1.D-1),BJ(1) )
      YMGOX = MIN( MAX( YMGOX+DYMGOX,0.D+0 ),1.D+0 )
      IF( ABS(BJ(1)).GT.1.D-7 ) GOTO 10
!
!---  Converged solution  ---
!
      PEQHX = PSAX/(1.D+0-YMGOX)
      PVX = PEQHX
!
!---  Inhibitor correction  ---
!
      IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
        CALL HYD_PO( YMGOX,POHX )
      ENDIF
      CALL HYD_T( TEQHX,YMGOX,XLSX,PEQHX,POHX,RHOHX,XHAX,XHOX,XHWX )
      PVAX = PSAX
      PVOX = PVX-PVAX
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
      PPLX = MAX( PLX+PATM,PSBX )
      CALL DENS_B( TX,PPLX,XLSX,RHOBX )
      CALL DENS_L( TX,RHOBX,XLAX,XLOX,RHOLX )
      IF( ISLC(44).EQ.1 ) THEN
        PVWX = PSBX
      ELSE
        CALL VPL_BX( TX,PLX,PVAX,PVOX,PVWX,PSBX,RHOLX,XLSX )
      ENDIF
      PGX = MAX( PVOX+PVAX+PVWX-PATM,PLX+(ENPR/BGL) )
      PNC = MAX( (BNL*PLX+BGN*PGX)/(BNL+BGN),(PLX*BNL+ENPR)/BNL )
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
      INDX = 0
      CALL KSP37( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
!
!---  Check for no liquid-CO2 conditions  ---
!
      IF( PNX.LE.PNC ) THEN
        INDX = 10
        RLMSG = PNX - PNC
        CHMSG = 'FLH_395: No Liquid-CO2: PN - PNC = '
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_395 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_PC1( PEQHX,PGX,PHX,PIX,PLX,PNX,PSBX,PVAX,
     &  PVOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,
     &  SLX,SNX,TX,TEQHX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,IZN,N )
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
!     Aqueous or aqueous-gas, no hydrate, no liquid CO2
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     CO2 mass - CO2 vapor pressure
!     CH4 mass - CH4 vapor pressure
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
      CALL SFT_L( TX,XLSX,SIG_GL )
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
!---  Hydrate equilibrium temperature  ---
!
      IF( (PVAX+PVOX).GT.EPSL ) THEN
        YMGOX = PVOX/(PVAX+PVOX)
        PEQHX = PVOX+PVAX
!
!---    Inhibitor correction  ---
!
        IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
          CALL HYD_PO( YMGOX,POHX )
        ENDIF
        CALL HYD_T( TEQHX,YMGOX,XLSX,PEQHX,POHX,RHOHX,XHAX,XHOX,XHWX )
      ELSE
        YMGOX = 0.D+0
        TEQHX = -TABS
      ENDIF
      XHAX = 0.D+0
      XHOX = 0.D+0
      XHWX = 1.D+0
      RHOHX = 0.D+0
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
      PGX = MAX( PVAX+PVOX+PVWX-PATM,PLX+ENPR/BGL )
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
      CALL KSP37( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_PC1 group  ---
!
      RETURN
      END
      
#ifdef eq_hyd
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_PC2( PEQHX,PGX,PHX,PIX,PLX,PNX,PSBX,PVAX,
     &  PVOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,
     &  TX,TEQHX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,IZN,N )
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
!     Aqueous or aqueous-gas, hydrate, no liquid CO2
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     CO2 mass - CO2 vapor pressure
!     CH4 mass - CH4 vapor pressure
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
      USE GRID
      USE FDVH
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER(KP=6)
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
      CALL SFT_L( TX,XLSX,SIG_GL )
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
!---  Pore-space equilibrium temperature  ---
!
      TEQHX = TX
!
!---  Hydrate equilibrium pressure  ---
!
      CALL HYD_P( TEQHX,YMGOX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
      PVOX = YMGOX*PEQHX
      PVAX = PEQHX - PVOX
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
      PGX = MAX( PVAX+PVOX+PVWX-PATM,PLX+ENPR/BGL )
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
      CALL KSP37( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_PC2 group  ---
!
      RETURN
      END
      
#endif
#ifdef k1_hyd
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_PC2( PEQHX,PGX,PHX,PIX,PLX,PNX,PSBX,PVAX,
     &  PVOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,
     &  TX,TEQHX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMGOHX,YMGOOX,IZN,N )
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
!     Aqueous or aqueous-gas, hydrate, no liquid CO2
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     CO2 mass - CO2 vapor pressure
!     CH4 mass - CH4 vapor pressure
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
      USE GRID
      USE FDVH
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER(KP=6)
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
      CALL SFT_L( TX,XLSX,SIG_GL )
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
!---  Pore-space equilibrium temperature  ---
!
      TEQHX = TX
!
!---  Hydrate equilibrium pressure  ---
!
      YMGOHX = (YMGOOX-YMGOX)*EXP(-DT*CHKN(4)) + YMGOX
      CALL HYD_P( TEQHX,YMGOHX,XLSX,PEQHX,RHOHX,XHAX,XHOX,XHWX )
      PVOX = YMGOX*PEQHX
      PVAX = PEQHX - PVOX
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
      PGX = MAX( PVAX+PVOX+PVWX-PATM,PLX+ENPR/BGL )
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
      CALL KSP37( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_PC2 group  ---
!
      RETURN
      END
      
#endif
#ifdef k2_hyd
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_PC2( PEQHX,PGX,PHX,PIX,PLX,PNX,PSBX,
     &  PVAX,PVAHX,PVAHOX,PVOX,PVOHX,PVOHOX,PVWX,
     &  RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,
     &  TX,TEQHX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,IZN,N )
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
!     Aqueous or aqueous-gas, hydrate, no liquid CO2
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     CO2 mass - CO2 vapor pressure
!     CH4 mass - CH4 vapor pressure
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
      USE GRID
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
      SAVE IALLOC
      DATA IALLOC /0/
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: AJ,GX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: BJ
      INTEGER, DIMENSION(:), ALLOCATABLE :: IJ
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_PC2'
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
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      CALL SFT_L( TX,XLSX,SIG_GL )
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
!---  Guess gas CH4 and CH2 vapor pressures  ---
!
      NC = 0
      DPVAY = MAX( 1.D-1,1.D-6*PVAX )
      DPVOY = MAX( 1.D-1,1.D-6*PVOX )
  200 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 12
        IMSG = N_DB
        CHMSG = 'Flash Phase Condition #2 Convergence Failure: ' // 
     &    'Newton-Raphson Loop @ Node'
        CALL WRMSGS( INDX )
      ENDIF
      DO 100 M = 1,3
        PVAY = PVAX
        PVOY = PVOX
        IF( M.EQ.2 ) PVAY = PVAX + DPVAY
        IF( M.EQ.3 ) PVOY = PVOX + DPVOY
        PVAHY = (PVAHOX-PVAY)*EXP(-DT*CHKN(4)) + PVAY
        PVOHY = (PVOHOX-PVOY)*EXP(-DT*CHKN(4)) + PVOY
!        PVAHX = (PVAHY-PVAHOX*EXP(-DT*CHKN(4)))/(1.D+0-EXP(-DT*CHKN(4)))
!        PVOHX = (PVOHY-PVOHOX*EXP(-DT*CHKN(4)))/(1.D+0-EXP(-DT*CHKN(4)))
        YMGOHX = PVOHY/(PVOHY+PVAHY)
        PEQHX = PVOY+PVAY
!
!---  Inhibitor correction  ---
!
        IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
          CALL HYD_PO( YMGOHX,POHX )
        ENDIF
        CALL HYD_T( TEQHX,YMGOHX,XLSX,PEQHX,POHX,RHOHX,XHAX,XHOX,XHWX )
        GX(1,M) = YMGOX - PVOY/(PVOY+PVAY)
        GX(2,M) = TX - TEQHX
  100 CONTINUE
!
!---  Load solution vector and Jacobian matrix for
!     two-phase conditions  ---
!
      DO 110 M = 1,2
        AJ(M,1) = (GX(M,2)-GX(M,1))/DPVAY
        AJ(M,2) = (GX(M,3)-GX(M,1))/DPVOY
        BJ(M) = -GX(M,1)
  110 CONTINUE
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
!---  Limit gas CH4 and CH2 vapor pressures  ---
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
      IF( ABS(DPVAX).GT.1.D-1 ) GOTO 200
      IF( ABS(DPVOX).GT.1.D-1 ) GOTO 200
!
!---  Converged solution  ---
!
      PVAHX = (PVAHOX-PVAX)*EXP(-DT*CHKN(4)) + PVAX
      PVOHX = (PVOHOX-PVOX)*EXP(-DT*CHKN(4)) + PVOX
!      PVAX = (PVAHX-PVAHOX*EXP(-DT*CHKN(4)))/(1.D+0-EXP(-DT*CHKN(4)))
!      PVOX = (PVOHX-PVOHOX*EXP(-DT*CHKN(4)))/(1.D+0-EXP(-DT*CHKN(4)))
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
      PGX = MAX( PVAX+PVOX+PVWX-PATM,PLX+ENPR/BGL )
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
      CALL KSP37( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_PC2 group  ---
!
      RETURN
      END
      
#endif
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_PC3( PEQHX,PGX,PHX,PIX,PLX,PNX,PSBX,PVAX,
     &  PVOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,
     &  SLX,SNX,TX,TEQHX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,IZN,N )
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
!     Aqueous or aqueous-gas, no hydrate, low liquid CO2
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     CO2 mass - liquid CO2 saturation
!     CH4 mass - CH4 vapor pressure
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
      CALL SFT_L( TX,XLSX,SIG_GL )
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
!---  Hydrate equilibrium temperature  ---
!
      IF( (PVAX+PVOX).GT.EPSL ) THEN
        YMGOX = PVOX/(PVAX+PVOX)
        PEQHX = PVOX+PVAX
!
!---    Inhibitor correction  ---
!
        IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
          CALL HYD_PO( YMGOX,POHX )
        ENDIF
        CALL HYD_T( TEQHX,YMGOX,XLSX,PEQHX,POHX,RHOHX,XHAX,XHOX,XHWX )
      ELSE
        YMGOX = 0.D+0
        TEQHX = -TABS
      ENDIF
      XHAX = 0.D+0
      XHOX = 0.D+0
      XHWX = 1.D+0
      RHOHX = 0.D+0
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
      PGX = MAX( PVAX+PVOX+PVWX-PATM,PLX+ENPR/BGL,PNX+ENPR/BGN )
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 1
      CALL KSP37( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
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
      SUBROUTINE FLH_PC4( PEQHX,PGX,PHX,PIX,PLX,PNX,PSBX,PVAX,
     &  PVOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,
     &  SLX,SNX,TX,TEQHX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,IZN,N )
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
!     Aqueous or aqueous-gas, hydrate, low liquid CO2
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     CO2 mass - liquid CO2 saturation
!     CH4 mass - hydrate saturation
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
!---  Gas CH4 mole fraction of formers, 
!     using Newton-Raphson iteration  ---
!
      CALL SP_A( TX,PVAX )
      CALL SOL_LS( TX,XLSMX )
      XLSY = MIN( YLSX,XLSMX )
!
!---  CO2 vapor pressure exceeds hydrate equilibrium pressure  ---
!
      YMGOZ = 0.D+0
      CALL HYD_P( TX,YMGOZ,XLSX,PEQHZ,RHOHX,XHAX,XHOX,XHWX )
      IF( PVAX.GT.PEQHZ ) GOTO 110
      NC = 0
  100 CONTINUE
      NC = NC + 1
      YMGOY = YMGOX
      PVOY = YMGOY*PVAX/(1.D+0-YMGOY)
      PEQHY = PVOY+PVAX
!
!---  Inhibitor correction  ---
!
      IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
        CALL HYD_PO( YMGOY,POHX )
      ENDIF
      CALL HYD_T( TEQHY,YMGOY,XLSX,PEQHY,POHX,RHOHX,XHAX,XHOX,XHWX )
      DYMGOZ = SIGN( 1.D-6,5.D-1-YMGOX )
      YMGOZ = YMGOX + DYMGOZ
      PVOZ = YMGOZ*PVAX/(1.D+0-YMGOZ)
      PEQHZ = PVOZ+PVAX
!
!---  Inhibitor correction  ---
!
      IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
        CALL HYD_PO( YMGOZ,POHX )
      ENDIF
      CALL HYD_T( TEQHZ,YMGOZ,XLSX,PEQHZ,POHX,RHOHX,XHAX,XHOX,XHWX )
      FX = TX - TEQHY
      FZ = TX - TEQHZ
      DFX = (FZ-FX)/DYMGOZ
      DYMGOX = SIGN( MIN( 1.D-1,ABS(FX/DFX) ),-FX/DFX )
      YMGOX = YMGOX + DYMGOX
      YMGOX = MAX( MIN( YMGOX,1.D+0 ),0.D+0 )
      IF( NC.GT.32 ) GOTO 110
      IF( ABS(DYMGOX).GT.1.D-8 ) GOTO 100
  110 CONTINUE
      PVOX = YMGOX*PVAX/(1.D+0-YMGOX)
      PEQHX = PVOX+PVAX
!
!---  Inhibitor correction  ---
!
      IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
        CALL HYD_PO( YMGOX,POHX )
      ENDIF
      CALL HYD_T( TEQHX,YMGOX,XLSX,PEQHX,POHX,RHOHX,XHAX,XHOX,XHWX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      CALL SFT_L( TX,XLSX,SIG_GL )
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
      PGX = MAX( PVAX+PVOX+PVWX-PATM,PLX+ENPR/BGL,PNX+ENPR/BGN )
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 1
      CALL KSP37( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
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
      SUBROUTINE FLH_PC5( PEQHX,PGX,PHX,PIX,PLX,PNX,PSBX,PVAX,
     &  PVOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,
     &  SLX,SNX,TX,TEQHX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,IZN,N )
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
!
!     Aqueous or aqueous-gas, no hydrate, liquid CO2
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     CO2 mass - liquid CO2 pressure
!     CH4 mass - CH4 vapor pressure
!     NaCl mass - total NaCl aqueous mass fraction
!
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
      CALL SFT_L( TX,XLSX,SIG_GL )
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
!---  Hydrate equilibrium temperature  ---
!
      IF( (PVAX+PVOX).GT.EPSL ) THEN
        YMGOX = PVOX/(PVAX+PVOX)
        PEQHX = PVOX+PVAX
!
!---    Inhibitor correction  ---
!
        IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
          CALL HYD_PO( YMGOX,POHX )
        ENDIF
        CALL HYD_T( TEQHX,YMGOX,XLSX,PEQHX,POHX,RHOHX,XHAX,XHOX,XHWX )
      ELSE
        YMGOX = 0.D+0
        TEQHX = -TABS
      ENDIF
      XHAX = 0.D+0
      XHOX = 0.D+0
      XHWX = 1.D+0
      RHOHX = 0.D+0
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
      PGX = MAX( PVAX+PVOX+PVWX-PATM,PLX+ENPR/BGL,PNX+ENPR/BGN )
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP37( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_PC5 group  ---
!
      RETURN
      END

#ifdef eq_hyd
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_PC6( PEQHX,PGX,PHX,PIX,PLX,PNX,PSBX,PVAX,
     &  PVOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,
     &  TX,TEQHX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,IZN,N )
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
!     Aqueous or aqueous-gas, hydrate, liquid CO2
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     CO2 mass - liquid CO2 pressure
!     CH4 mass - hydrate saturation
!     NaCl mass - total NaCl aqueous mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 4 January 2008.
!     Last Modified by MD White, PNNL, 4 January 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NCG_PT
      USE HYST
      USE GRID
      USE FDVH
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER(KP=6)
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
!---  Gas CH4 mole fraction of formers, 
!     using Newton-Raphson iteration  ---
!
      CALL SP_A( TX,PVAX )
      CALL SOL_LS( TX,XLSMX )
      XLSY = MIN( YLSX,XLSMX )
!!
!!---  CO2 vapor pressure exceeds hydrate equilibrium pressure  ---
!!
!      YMGOZ = 0.D+0
!      CALL HYD_P( TX,YMGOZ,XLSX,PEQHZ,RHOHX,XHAX,XHOX,XHWX )
!      IF( PVAX.GT.PEQHZ ) GOTO 110
      NC = 0
  100 CONTINUE
      NC = NC + 1
      YMGOY = YMGOX
      PVOY = YMGOY*PVAX/(1.D+0-YMGOY)
      PEQHY = PVOY+PVAX
!
!---  Inhibitor correction  ---
!
      IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
        CALL HYD_PO( YMGOY,POHX )
      ENDIF
      CALL HYD_T( TEQHY,YMGOY,XLSX,PEQHY,POHX,RHOHX,XHAX,XHOX,XHWX )
      DYMGOZ = SIGN( 1.D-6,5.D-1-YMGOX )
      YMGOZ = YMGOX + DYMGOZ
      PVOZ = YMGOZ*PVAX/(1.D+0-YMGOZ)
      PEQHZ = PVOZ+PVAX
!
!---  Inhibitor correction  ---
!
      IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
        CALL HYD_PO( YMGOZ,POHX )
      ENDIF
      CALL HYD_T( TEQHZ,YMGOZ,XLSX,PEQHZ,POHX,RHOHX,XHAX,XHOX,XHWX )
      FX = TX - TEQHY
      FZ = TX - TEQHZ
      DFX = (FZ-FX)/DYMGOZ
      DYMGOX = SIGN( MIN( 1.D-1,ABS(FX/DFX) ),-FX/DFX )
      YMGOX = YMGOX + DYMGOX
      YMGOX = MAX( MIN( YMGOX,1.D+0 ),0.D+0 )
      IF( NC.GT.32 ) GOTO 110
      IF( ABS(DYMGOX).GT.1.D-8 ) GOTO 100
  110 CONTINUE
      PVOX = YMGOX*PVAX/(1.D+0-YMGOX)
      PEQHX = PVOX+PVAX
!
!---  Inhibitor correction  ---
!
      IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
        CALL HYD_PO( YMGOX,POHX )
      ENDIF
      CALL HYD_T( TEQHX,YMGOX,XLSX,PEQHX,POHX,RHOHX,XHAX,XHOX,XHWX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      CALL SFT_L( TX,XLSX,SIG_GL )
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
      PGX = MAX( PVAX+PVOX+PVWX-PATM,PNX+ENPR/BGN )
!      PGX = MAX( PVAX+PVOX+PVWX-PATM,PLX+ENPR/BGL,PNX+ENPR/BGN )
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP37( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_PC6 group  ---
!
      RETURN
      END
      
#endif
#ifdef k1_hyd
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_PC6( PEQHX,PGX,PHX,PIX,PLX,PNX,PSBX,PVAX,
     &  PVOX,PVWX,RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,
     &  TX,TEQHX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,YMGOHX,YMGOOX,IZN,N )
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
!     Aqueous or aqueous-gas, hydrate, liquid CO2
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     CO2 mass - liquid CO2 pressure
!     CH4 mass - hydrate saturation
!     NaCl mass - total NaCl aqueous mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 4 January 2008.
!     Last Modified by MD White, PNNL, 4 January 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NCG_PT
      USE HYST
      USE GRID
      USE FDVH
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER(KP=6)
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
!---  Gas CH4 mole fraction of formers, 
!     using Newton-Raphson iteration  ---
!
      CALL SP_A( TX,PVAX )
      CALL SOL_LS( TX,XLSMX )
      XLSY = MIN( YLSX,XLSMX )
!!
!!---  CO2 vapor pressure exceeds hydrate equilibrium pressure  ---
!!
!      YMGOZ = 0.D+0
!      CALL HYD_P( TX,YMGOZ,XLSX,PEQHZ,RHOHX,XHAX,XHOX,XHWX )
!      IF( PVAX.GT.PEQHZ ) GOTO 110
      NC = 0
  100 CONTINUE
      NC = NC + 1
      YMGOY = YMGOX
      YMGOHX = (YMGOOX-YMGOY)*EXP(-DT*CHKN(4)) + YMGOY
      PVOY = YMGOY*PVAX/(1.D+0-YMGOY)
      PEQHY = PVOY+PVAX
!
!---  Inhibitor correction  ---
!
      IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
        CALL HYD_PO( YMGOHX,POHX )
      ENDIF
      CALL HYD_T( TEQHY,YMGOHX,XLSX,PEQHY,POHX,RHOHX,XHAX,XHOX,XHWX )
      DYMGOZ = SIGN( 1.D-6,5.D-1-YMGOX )
      YMGOZ = YMGOX + DYMGOZ
      YMGOHX = (YMGOOX-YMGOZ)*EXP(-DT*CHKN(4)) + YMGOZ
      PVOZ = YMGOZ*PVAX/(1.D+0-YMGOZ)
      PEQHZ = PVOZ+PVAX
!
!---  Inhibitor correction  ---
!
      IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
        CALL HYD_PO( YMGOHX,POHX )
      ENDIF
      CALL HYD_T( TEQHZ,YMGOHX,XLSX,PEQHZ,POHX,RHOHX,XHAX,XHOX,XHWX )
      FX = TX - TEQHY
      FZ = TX - TEQHZ
      DFX = (FZ-FX)/DYMGOZ
      DYMGOX = SIGN( MIN( 1.D-1,ABS(FX/DFX) ),-FX/DFX )
      YMGOX = YMGOX + DYMGOX
      YMGOX = MAX( MIN( YMGOX,1.D+0 ),0.D+0 )
      IF( NC.GT.32 ) GOTO 110
      IF( ABS(DYMGOX).GT.1.D-8 ) GOTO 100
  110 CONTINUE
      YMGOHX = (YMGOOX-YMGOX)*EXP(-DT*CHKN(4)) + YMGOX
      PVOX = YMGOX*PVAX/(1.D+0-YMGOX)
      PEQHX = PVOX+PVAX
!
!---  Inhibitor correction  ---
!
      IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
        CALL HYD_PO( YMGOHX,POHX )
      ENDIF
      CALL HYD_T( TEQHX,YMGOHX,XLSX,PEQHX,POHX,RHOHX,XHAX,XHOX,XHWX )
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      CALL SFT_L( TX,XLSX,SIG_GL )
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
      PGX = MAX( PVAX+PVOX+PVWX-PATM,PNX+ENPR/BGN )
!      PGX = MAX( PVAX+PVOX+PVWX-PATM,PLX+ENPR/BGL,PNX+ENPR/BGN )
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP37( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_PC6 group  ---
!
      RETURN
      END

#endif
#ifdef k2_hyd
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLH_PC6( PEQHX,PGX,PHX,PIX,PLX,PNX,PSBX,
     &  PVAX,PVAHX,PVAHOX,PVOX,PVOHX,PVOHOX,PVWX,
     &  RHOHX,RKGX,RKLX,RKNX,SHX,SIX,SGX,SLX,SNX,
     &  TX,TEQHX,XHAX,XHOX,XHWX,XLSX,YLSX,YMGOX,IZN,N )
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
!     Aqueous or aqueous-gas, hydrate, liquid CO2
!
!     Energy - temperature
!     Water mass - aqueous pressure
!     CO2 mass - liquid CO2 pressure
!     CH4 mass - hydrate saturation
!     NaCl mass - total NaCl aqueous mass fraction
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 4 January 2008.
!     Last Modified by MD White, PNNL, 4 January 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE NCG_PT
      USE HYST
      USE GRID
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
      SAVE IALLOC
      DATA IALLOC /0/
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: AJ,GX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: BJ
      INTEGER, DIMENSION(:), ALLOCATABLE :: IJ
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLH_PC6'
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
!---  Gas CH4 mole fraction of formers, 
!     using Newton-Raphson iteration  ---
!
      CALL SP_A( TX,PVAX )
      CALL SOL_LS( TX,XLSMX )
      XLSY = MIN( YLSX,XLSMX )
      PVAHX = (PVAHOX-PVAX)*EXP(-DT*CHKN(4)) + PVAX
!
!---  Guess gas CH4 vapor pressures  ---
!
      NC = 0
      DPVOY = MAX( 1.D-1,1.D-6*PVOX )
  200 CONTINUE
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 12
        IMSG = N_DB
        CHMSG = 'Flash Phase Condition #6 Convergence Failure: ' // 
     &    'Newton-Raphson Loop @ Node'
        CALL WRMSGS( INDX )
      ENDIF
      DO 100 M = 1,2
        PVOY = PVOX
        IF( M.EQ.2 ) PVOY = PVOX + DPVOY
        PVOHY = (PVOHOX-PVOY)*EXP(-DT*CHKN(4)) + PVOY
        YMGOHX = PVOHY/(PVOHY+PVAHX)
        PEQHX = PVOY+PVAX
!
!---  Inhibitor correction  ---
!
        IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
          CALL HYD_PO( YMGOHX,POHX )
        ENDIF
        CALL HYD_T( TEQHX,YMGOHX,XLSX,PEQHX,POHX,RHOHX,XHAX,XHOX,XHWX )
        GX(1,M) = TX - TEQHX
  100 CONTINUE
!
!---  Load solution vector and Jacobian matrix for
!     two-phase conditions  ---
!
      AJ(1,1) = (GX(1,2)-GX(1,1))/DPVOY
      BJ(1) = -GX(1,1)
!
!---  Update primary unknowns  ---
!
      DPVOX = BJ(1)/AJ(1,1)
      PVOX = PVOX + DPVOX
!
!---  Limit gas CH4 vapor pressures  ---
!
      IF( PVOX/EPSL.LT.EPSL .AND. DPVOX.LT.0.D+0 ) THEN
        PVOX = 0.D+0
        DPVOX = 0.D+0
      ENDIF
!
!---  Convergence check  ---
!
      IF( ABS(DPVOX).GT.1.D-1 ) GOTO 200
!
!---  Converged solution  ---
!
      YMGOX = PVOX/(PVOX+PVAX)
      PVOHX = (PVOHOX-PVOX)*EXP(-DT*CHKN(4)) + PVOX
!
!---  Interfacial surface tensions  ---
!
      SIG_GL = SCHR(5,IZN)
      SIG_NL = SCHR(6,IZN)
      SIG_GN = SCHR(7,IZN)
      SIG_HL = SCHR(9,IZN)
      SIG_IL = SCHR(9,IZN)
      CALL SFT_L( TX,XLSX,SIG_GL )
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
      PGX = MAX( PVAX+PVOX+PVWX-PATM,PNX+ENPR/BGN )
!
!---  Ice pressure  ---
!
      CALL ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
!
!---  Phase saturations  ---
!
      INDX = 0
      CALL KSP37( PGX,PIX,PLX,PNX,SGX,SHX,SIX,SLX,SNX,
     &  RKGX,RKLX,RKNX,TX,INDX,IZN,N )
      INDX = 0
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLH_PC6 group  ---
!
      RETURN
      END
#endif

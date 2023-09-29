!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCIH_PI( PX,VAR_HT,VAR2P_HT,VARX,ZMCX )
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
!     Bicubic spline interpolation. Cubic splines are created at the
!     tabulated values of CO2 and at the CO2-CH4 diagonal.  A cubic
!     spline is then created using the CO2 splines at the CH4 mole
!     fraction.  The returned value is computed from the cubic spline
!     created at the CH4 mole fraction, evaluated at the CO2 mole
!     fraction.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 December 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYDT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VAR_HT(LHE_HT,LCH_HT),VAR2P_HT(LHE_HT,LCH_HT)
      REAL*8 ZMCX(LHF_HT)
      REAL*8 YCH4X(LCN_HT),XCH4X(LCN_HT),Y2CH4X(LCN_HT)
      REAL*8 YCO2X(LCN_HT+1),XCO2X(LCN_HT+1),Y2CO2X(LCN_HT+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BCIH_PI'
      TOLX = 1.D-12
      ZB1X = 0.D+0
      ZB2X = 0.D+0
!
!---  Locate the CO2 index  ---
!
      IX = NINT(ZMIH_HT(1)/DZHC_HT) + 1
      ZX = REAL(IX-1)*DZHC_HT + ZB1X
      IMX = 2
      IPX = 2
      IF( (IX-IMX).LT.1 ) THEN
        I1X = 1
        I2X = 5
      ELSEIF( (IX+IPX).GT.NZH_HT ) THEN
        I1X = NZH_HT - 4
        I2X = NZH_HT
      ELSE
        I1X = IX - IMX
        I2X = IX + IPX
      ENDIF
!
!---  Locate the CH4 index  ---
!
      JX = NINT(ZMIH_HT(2)/DZHC_HT) + 1
      ZX = REAL(JX-1)*DZHC_HT + ZB2X
      JMX = 2
      JPX = 2
      IF( (JX-JMX).LT.1 ) THEN
        J1X = 1
        J2X = 5
      ELSEIF( (JX+JPX).GT.NZH_HT ) THEN
        J1X = NZH_HT - 4
        J2X = NZH_HT
      ELSE
        J1X = JX - JMX
        J2X = JX + JPX
      ENDIF
!
!---  CO2 concentration at the CH4 diagonal  ---
!
      XCO2DX = 1.D+0-ZMCX(2)
!
!---  Create cubic splines at the tabulated CO2 concentrations, and 
!     along the diagonal of the tabulated CH4-CO2 concentrations  ---
!
      NCI = 0
      DO 40 I = I1X,I2X
!
!---    Splines at the tabulated CO2 concentrations  ---
!
        NCJ = 0
        DO 20 J = J1X,J2X
          NCJ = NCJ + 1
          XCH4X(NCJ) = REAL(J-1)*DZHC_HT
          NC = IZH_HT(I,J)
          NP = NHEP_HT(NC)
          IF( PX.LE.PHE_HT(1,NC) ) THEN
            YCH4X(NCJ) = VAR_HT(1,NC)
          ELSEIF( PX.GE.PHE_HT(NP,NC) ) THEN
            YCH4X(NCJ) = VAR_HT(NP,NC)
          ELSE
            CALL SPLINT( PHE_HT(1,NC),VAR_HT(1,NC),VAR2P_HT(1,NC),
     &        NP,PX,YCH4X(NCJ))
          ENDIF
   20   CONTINUE
        CALL SPLINE( XCH4X,YCH4X,NCJ,Y2CH4X )
        NCI = NCI + 1
        XCO2X(NCI) = REAL(I-1)*DZHC_HT + ZB1X
        CALL SPLINT( XCH4X,YCH4X,Y2CH4X,NCJ,ZMCX(2),YCO2X(NCI) )
!
!---    Diagonal spline  ---
!
        XCO2PX = REAL(I)*DZHC_HT + ZB1X
        IF( XCO2DX.GT.XCO2X(NCI) .AND. XCO2DX.LT.XCO2PX ) THEN
!
!---      Splines along the diagonal of the tabulated CH4-CO2 
!         concentrations  ---
!
          NCJ = 0
          DO 30 J = 1,NZH_HT
            NCJ = NCJ + 1
            XCH4X(NCJ) = REAL(J-1)*DZHC_HT
            IX = NZH_HT - J + 1
            NC = IZH_HT(IX,J)
            NP = NHEP_HT(NC)
            IF( PX.LE.PHE_HT(1,NC) ) THEN
              YCH4X(NCJ) = VAR_HT(1,NC)
            ELSEIF( PX.GE.PHE_HT(NP,NC) ) THEN
              YCH4X(NCJ) = VAR_HT(NP,NC)
            ELSE
              CALL SPLINT( PHE_HT(1,NC),VAR_HT(1,NC),VAR2P_HT(1,NC),
     &          NP,PX,YCH4X(NCJ))
            ENDIF
   30     CONTINUE
          CALL SPLINE( XCH4X,YCH4X,NCJ,Y2CH4X )
          NCI = NCI + 1
          XCO2X(NCI) = XCO2DX
          CALL SPLINT( XCH4X,YCH4X,Y2CH4X,NCJ,ZMCX(2),YCO2X(NCI) )
        ENDIF
   40 CONTINUE
!
!---  Create a cubic spline across the tabulated CO2 concentrations
!     and CO2 concentration at the CH4 diagonal 
!     at the CH4 concentration  ---
!
      CALL SPLINE( XCO2X,YCO2X,NCI,Y2CO2X )
!
!---  Compute the hydrate equilibrium value at the CO2 concentration,
!     using the cubic spline at the CH4 concentration  ---
!
      CALL SPLINT( XCO2X,YCO2X,Y2CO2X,NCI,ZMCX(1),VARX )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCIH_PI group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCIH_TI( TKX,VAR_HT,VAR2T_HT,VARX,ZMCX )
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
!     Bicubic spline interpolation. Cubic splines are created at the
!     tabulated values of CO2 and at the CO2-CH4 diagonal.  A cubic
!     spline is then created using the CO2 splines at the CH4 mole
!     fraction.  The returned value is computed from the cubic spline
!     created at the CH4 mole fraction, evaluated at the CO2 mole
!     fraction.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 December 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYDT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VAR_HT(LHE_HT,LCH_HT),VAR2T_HT(LHE_HT,LCH_HT)
      REAL*8 ZMCX(LHF_HT)
      REAL*8 YCH4X(LCN_HT),XCH4X(LCN_HT),Y2CH4X(LCN_HT)
      REAL*8 YCO2X(LCN_HT+1),XCO2X(LCN_HT+1),Y2CO2X(LCN_HT+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BCIH_TI'
      TOLX = 1.D-12
      ZB1X = 0.D+0
      ZB2X = 0.D+0
!
!---  Locate the CO2 index  ---
!
      IX = NINT(ZMIH_HT(1)/DZHC_HT) + 1
      ZX = REAL(IX-1)*DZHC_HT + ZB1X
      IMX = 2
      IPX = 2
      IF( (IX-IMX).LT.1 ) THEN
        I1X = 1
        I2X = 5
      ELSEIF( (IX+IPX).GT.NZH_HT ) THEN
        I1X = NZH_HT - 4
        I2X = NZH_HT
      ELSE
        I1X = IX - IMX
        I2X = IX + IPX
      ENDIF
!
!---  Locate the CH4 index  ---
!
      JX = NINT(ZMIH_HT(2)/DZHC_HT) + 1
      ZX = REAL(JX-1)*DZHC_HT + ZB2X
      JMX = 2
      JPX = 2
      IF( (JX-JMX).LT.1 ) THEN
        J1X = 1
        J2X = 5
      ELSEIF( (JX+JPX).GT.NZH_HT ) THEN
        J1X = NZH_HT - 4
        J2X = NZH_HT
      ELSE
        J1X = JX - JMX
        J2X = JX + JPX
      ENDIF
!
!---  CO2 concentration at the CH4 diagonal  ---
!
      XCO2DX = 1.D+0-ZMCX(2)
!
!---  Create cubic splines at the tabulated CO2 concentrations, and 
!     along the diagonal of the tabulated CH4-CO2 concentrations  ---
!
      NCI = 0
      DO 40 I = I1X,I2X
!
!---    Splines at the tabulated CO2 concentrations  ---
!
        NCJ = 0
        DO 20 J = J1X,J2X
          NCJ = NCJ + 1
          XCH4X(NCJ) = REAL(J-1)*DZHC_HT
          NC = IZH_HT(I,J)
          NP = NHEP_HT(NC)
          IF( TKX.LE.THE_HT(1,NC) ) THEN
            YCH4X(NCJ) = VAR_HT(1,NC)
          ELSEIF( TKX.GE.THE_HT(NP,NC) ) THEN
            YCH4X(NCJ) = VAR_HT(NP,NC)
          ELSE
            CALL SPLINT( THE_HT(1,NC),VAR_HT(1,NC),VAR2T_HT(1,NC),
     &        NP,TKX,YCH4X(NCJ) )
          ENDIF
   20   CONTINUE
        CALL SPLINE( XCH4X,YCH4X,NCJ,Y2CH4X )
        NCI = NCI + 1
        XCO2X(NCI) = REAL(I-1)*DZHC_HT + ZB1X
        CALL SPLINT( XCH4X,YCH4X,Y2CH4X,NCJ,ZMCX(2),YCO2X(NCI) )
!
!---    Diagonal spline  ---
!
        XCO2PX = REAL(I)*DZHC_HT + ZB1X
        IF( XCO2DX.GT.XCO2X(NCI) .AND. XCO2DX.LT.XCO2PX ) THEN
!
!---      Splines along the diagonal of the tabulated CH4-CO2 
!         concentrations  ---
!
          NCJ = 0
          DO 30 J = 1,NZH_HT
            NCJ = NCJ + 1
            XCH4X(NCJ) = REAL(J-1)*DZHC_HT
            IX = NZH_HT - J + 1
            NC = IZH_HT(IX,J)
            NP = NHEP_HT(NC)
            IF( TKX.LE.THE_HT(1,NC) ) THEN
              YCH4X(NCJ) = VAR_HT(1,NC)
            ELSEIF( TKX.GE.THE_HT(NP,NC) ) THEN
              YCH4X(NCJ) = VAR_HT(NP,NC)
            ELSE
              CALL SPLINT( THE_HT(1,NC),VAR_HT(1,NC),VAR2T_HT(1,NC),
     &          NP,TKX,YCH4X(NCJ) )
            ENDIF
   30     CONTINUE
          CALL SPLINE( XCH4X,YCH4X,NCJ,Y2CH4X )
          NCI = NCI + 1
          XCO2X(NCI) = XCO2DX
          CALL SPLINT( XCH4X,YCH4X,Y2CH4X,NCJ,ZMCX(2),YCO2X(NCI) )
        ENDIF
   40 CONTINUE
!
!---  Create a cubic spline across the tabulated CO2 concentrations
!     and CO2 concentration at the CH4 diagonal 
!     at the CH4 concentration  ---
!
      CALL SPLINE( XCO2X,YCO2X,NCI,Y2CO2X )
!
!---  Compute the hydrate equilibrium value at the CO2 concentration,
!     using the cubic spline at the CH4 concentration  ---
!
      CALL SPLINT( XCO2X,YCO2X,Y2CO2X,NCI,ZMCX(1),VARX )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCIH_TI group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCIP_SC( VAR_HT,VARX,ZMCX )
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
!     Bicubic spline interpolation. Cubic splines are created at the
!     tabulated values of CO2 and at the CO2-CH4 diagonal.  A cubic
!     spline is then created using the CO2 splines at the CH4 mole
!     fraction.  The returned value is computed from the cubic spline
!     created at the CH4 mole fraction, evaluated at the CO2 mole
!     fraction.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 December 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYDT
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VAR_HT(LCP_HT)
      REAL*8 ZMCX(LHF_HT)
      REAL*8 YCH4X(LCN_HT),XCH4X(LCN_HT),Y2CH4X(LCN_HT)
      REAL*8 YCO2X(LCN_HT+1),XCO2X(LCN_HT+1),Y2CO2X(LCN_HT+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BCIP_SC'
      TOLX = 1.D-12
      ZB1X = 2.D-1
      ZB2X = 0.D+0
!
!---  Locate the CO2 index  ---
!
      IX = INT( (ZMCX(1)-ZB1X+TOLX)/DZPC_HT ) + 1
      ZX = (REAL(IX)+5.D-1)*DZPC_HT
      IF( (ZMCX(1)-ZX).GT.EPSL ) IX = IX + 1
      IF( (IX-2).LT.1 ) THEN
        I1X = 1
        I2X = 5
      ELSEIF( (IX+2).GT.NZP_HT ) THEN
        I1X = NZP_HT - 4
        I2X = NZP_HT
      ELSE
        I1X = IX - 2
        I2X = IX + 2
      ENDIF
!
!---  Locate the CH4 index  ---
!
      JX = INT( (ZMCX(2)-ZB2X+TOLX)/DZPC_HT ) + 1
      ZX = (REAL(JX)+5.D-1)*DZPC_HT
      IF( (ZMCX(2)-ZX).GT.EPSL ) JX = JX + 1
      IF( (JX-2).LT.1 ) THEN
        J1X = 1
        J2X = 5
      ELSEIF( (JX+2).GT.NZP_HT ) THEN
        J1X = NZP_HT - 4
        J2X = NZP_HT
      ELSE
        J1X = JX - 2
        J2X = JX + 2
      ENDIF
!
!---  CO2 concentration at the CH4 diagonal  ---
!
      XCO2DX = 1.D+0-ZMCX(2)
!
!---  Create cubic splines at the tabulated CO2 concentrations
!     and potentially the CO2 concentration at the CH4 diagonal  ---
!
      NCI = 0
      DO 30 I = I1X,I2X
!
!---    CO2 concentration at the CH4 diagonal is not a tabulated value, 
!       include diagonal in the CH4 cubic spline  ---
!
        IF( I.GT.I1X ) THEN
          XCO2MX = REAL(I-2)*DZPC_HT + ZB1X
          XCO2PX = XCO2MX + DZPC_HT
          IF( ((XCO2DX-XCO2MX).GT.EPSL) .AND. 
     &      ((XCO2PX-XCO2DX).GT.EPSL) ) THEN
            NCJ = 0
            DO 10 J = J1X,J2X
              NCJ = NCJ + 1
              XCH4X(NCJ) = REAL(J-1)*DZPC_HT
              K = NZP_HT+1-J
              YCH4X(NCJ) = VAR_HT(IZP_HT(K,J))
   10       CONTINUE
            CALL SPLINE( XCH4X,YCH4X,NCJ,Y2CH4X )
            NCI = NCI + 1
            XCO2X(NCI) = XCO2DX
            CALL SPLINT( XCH4X,YCH4X,Y2CH4X,NCJ,ZMCX(2),YCO2X(NCI) )
          ENDIF
        ENDIF
!
!---    CO2 concentration at tabulated values  ---
!
        NCJ = 0
        DO 20 J = J1X,J2X
          NCJ = NCJ + 1
          XCH4X(NCJ) = REAL(J-1)*DZPC_HT
          YCH4X(NCJ) = VAR_HT(IZP_HT(I,J))
   20   CONTINUE
        CALL SPLINE( XCH4X,YCH4X,NCJ,Y2CH4X )
        NCI = NCI + 1
        XCO2X(NCI) = REAL(I-1)*DZPC_HT + ZB1X
        CALL SPLINT( XCH4X,YCH4X,Y2CH4X,NCJ,ZMCX(2),YCO2X(NCI) )
   30 CONTINUE
!
!---  Create a cubic spline across the tabulated CO2 concentrations
!     and CO2 concentration at the CH4 diagonal 
!     at the CH4 concentration  ---
!
      CALL SPLINE( XCO2X,YCO2X,NCI,Y2CO2X )
!
!---  Compute the hydrate equilibrium value at the CO2 concentration,
!     using the cubic spline at the CH4 concentration  ---
!
      CALL SPLINT( XCO2X,YCO2X,Y2CO2X,NCI,ZMCX(1),VARX )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCIP_SC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCIP_TS( D2VAR_HT,VAR_HT,VARX,TCTX,TKX,T_HT,ZMCX,
     &  N_HT )
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
!     Bicubic spline interpolation. Cubic splines are created at the
!     tabulated values of CO2 and at the CO2-CH4 diagonal.  A cubic
!     spline is then created using the CO2 splines at the CH4 mole
!     fraction.  The returned value is computed from the cubic spline
!     created at the CH4 mole fraction, evaluated at the CO2 mole
!     fraction.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 December 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYDT
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 D2VAR_HT(LPE_HT,LCP_HT),VAR_HT(LPE_HT,LCP_HT)
      REAL*8 T_HT(LPE_HT,LCP_HT)
      REAL*8 ZMCX(LHF_HT)
      REAL*8 YCH4X(LCN_HT),XCH4X(LCN_HT),Y2CH4X(LCN_HT)
      REAL*8 YCO2X(LCN_HT+1),XCO2X(LCN_HT+1),Y2CO2X(LCN_HT+1)
      INTEGER N_HT(LCP_HT)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BCIP_TS'
      TOLX = 1.D-12
      ZB1X = 2.D-1
      ZB2X = 0.D+0
!
!---  Locate the CO2 index  ---
!
      IX = INT( (ZMCX(1)-ZB1X+TOLX)/DZPC_HT ) + 1
      ZX = (REAL(IX)+5.D-1)*DZPC_HT
      IF( (ZMCX(1)-ZX).GT.EPSL ) IX = IX + 1
      IF( (IX-2).LT.1 ) THEN
        I1X = 1
        I2X = 5
      ELSEIF( (IX+2).GT.NZP_HT ) THEN
        I1X = NZP_HT - 4
        I2X = NZP_HT
      ELSE
        I1X = IX - 2
        I2X = IX + 2
      ENDIF
!
!---  Locate the CH4 index  ---
!
      JX = INT( (ZMCX(2)-ZB2X+TOLX)/DZPC_HT ) + 1
      ZX = (REAL(JX)+5.D-1)*DZPC_HT
      IF( (ZMCX(2)-ZX).GT.EPSL ) JX = JX + 1
      IF( (JX-2).LT.1 ) THEN
        J1X = 1
        J2X = 5
      ELSEIF( (JX+2).GT.NZP_HT ) THEN
        J1X = NZP_HT - 4
        J2X = NZP_HT
      ELSE
        J1X = JX - 2
        J2X = JX + 2
      ENDIF
!
!---  CO2 concentration at the CH4 diagonal  ---
!
      XCO2DX = 1.D+0-ZMCX(2)
!
!---  Create cubic splines at the tabulated CO2 concentrations
!     and potentially the CO2 concentration at the CH4 diagonal  ---
!
      NCI = 0
      DO 30 I = I1X,I2X
!
!---    CO2 concentration at the CH4 diagonal is not a tabulated value, 
!       include diagonal in the CH4 cubic spline  ---
!
        IF( I.GT.I1X ) THEN
          XCO2MX = REAL(I-2)*DZPC_HT + ZB1X
          XCO2PX = XCO2MX + DZPC_HT
          IF( ((XCO2DX-XCO2MX).GT.EPSL) .AND. 
     &      ((XCO2PX-XCO2DX).GT.EPSL) ) THEN
            NCJ = 0
            DO 10 J = J1X,J2X
              NCJ = NCJ + 1
              XCH4X(NCJ) = REAL(J-1)*DZPC_HT
              K = NZP_HT+1-J
              TCTY = TCT_HT(IZP_HT(K,J))
              TKZ = TCTY*TKX/TCTX
              IZX = IZP_HT(K,J)
              CALL SPLINT( T_HT(1,IZX),VAR_HT(1,IZX),
     &          D2VAR_HT(1,IZX),N_HT(IZX),TKZ,YCH4X(NCJ) )
   10       CONTINUE
            CALL SPLINE( XCH4X,YCH4X,NCJ,Y2CH4X )
            NCI = NCI + 1
            XCO2X(NCI) = XCO2DX
            CALL SPLINT( XCH4X,YCH4X,Y2CH4X,NCJ,ZMCX(2),YCO2X(NCI) )
          ENDIF
        ENDIF
!
!---    CO2 concentration at tabulated values  ---
!
        NCJ = 0
        DO 20 J = J1X,J2X
          NCJ = NCJ + 1
          XCH4X(NCJ) = REAL(J-1)*DZPC_HT
          TCTY = TCT_HT(IZP_HT(I,J))
          TKZ = TCTY*TKX/TCTX
          IZX = IZP_HT(I,J)
          CALL SPLINT( T_HT(1,IZX),VAR_HT(1,IZX),
     &      D2VAR_HT(1,IZX),N_HT(IZX),TKZ,YCH4X(NCJ) )
   20   CONTINUE
        CALL SPLINE( XCH4X,YCH4X,NCJ,Y2CH4X )
        NCI = NCI + 1
        XCO2X(NCI) = REAL(I-1)*DZPC_HT + ZB1X
        CALL SPLINT( XCH4X,YCH4X,Y2CH4X,NCJ,ZMCX(2),YCO2X(NCI) )
   30 CONTINUE
!
!---  Create a cubic spline across the tabulated CO2 concentrations
!     and CO2 concentration at the CH4 diagonal 
!     at the CH4 concentration  ---
!
      CALL SPLINE( XCO2X,YCO2X,NCI,Y2CO2X )
!
!---  Compute the hydrate equilibrium value at the CO2 concentration,
!     using the cubic spline at the CH4 concentration  ---
!
      CALL SPLINT( XCO2X,YCO2X,Y2CO2X,NCI,ZMCX(1),VARX )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCIP_TS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCIP_TPS( D2V2P_HT,V2P_HT,VARX,PPELX,PPEUX,PX,
     &  TCTX,TKX,ZMCX,N_HT )
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
!     Bicubic spline interpolation. Cubic splines are created at the
!     tabulated values of CO2 and at the CO2-CH4 diagonal.  A cubic
!     spline is then created using the CO2 splines at the CH4 mole
!     fraction.  The returned value is computed from the cubic spline
!     created at the CH4 mole fraction, evaluated at the CO2 mole
!     fraction.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, September 3 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYDT
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 D2V2P_HT(LTP_HT,LPP_HT,LCP_HT),V2P_HT(LTP_HT,LPP_HT,LCP_HT)
      REAL*8 ZMCX(LHF_HT)
      REAL*8 YCH4X(LCN_HT),XCH4X(LCN_HT),Y2CH4X(LCN_HT)
      REAL*8 YCO2X(LCN_HT+1),XCO2X(LCN_HT+1),Y2CO2X(LCN_HT+1)
      INTEGER N_HT(LCP_HT)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BCIP_TPS'
      TOLX = 1.D-12
      ZB1X = 2.D-1
      ZB2X = 0.D+0
!
!---  Locate the CO2 index  ---
!
      IX = INT( (ZMCX(1)+TOLX-ZB1X)/DZPC_HT ) + 1
      ZX = (REAL(IX)+5.D-1)*DZPC_HT
      IF( (ZMCX(1)-ZX).GT.EPSL ) IX = IX + 1
      IF( (IX-2).LT.1 ) THEN
        I1X = 1
        I2X = 5
      ELSEIF( (IX+2).GT.NZP_HT ) THEN
        I1X = NZP_HT - 4
        I2X = NZP_HT
      ELSE
        I1X = IX - 2
        I2X = IX + 2
      ENDIF
!
!---  Locate the CH4 index  ---
!
      JX = INT( (ZMCX(2)+TOLX-ZB2X)/DZPC_HT ) + 1
      ZX = (REAL(JX)+5.D-1)*DZPC_HT
      IF( (ZMCX(2)-ZX).GT.EPSL ) JX = JX + 1
      IF( (JX-2).LT.1 ) THEN
        J1X = 1
        J2X = 5
      ELSEIF( (JX+2).GT.NZP_HT ) THEN
        J1X = NZP_HT - 4
        J2X = NZP_HT
      ELSE
        J1X = JX - 2
        J2X = JX + 2
      ENDIF
!
!---  CO2 concentration at the CH4 diagonal  ---
!
      XCO2DX = 1.D+0-ZMCX(2)
!
!---  Create cubic splines at the tabulated CO2 concentrations
!     and potentially the CO2 concentration at the CH4 diagonal  ---
!
      NCI = 0
      DO 30 I = I1X,I2X
!
!---    CO2 concentration at the CH4 diagonal is not a tabulated value, 
!       include diagonal in the CH4 cubic spline  ---
!
        IF( I.GT.I1X ) THEN
          XCO2MX = REAL(I-2)*DZPC_HT + ZB1X
          XCO2PX = XCO2MX + DZPC_HT
          IF( ((XCO2DX-XCO2MX).GT.EPSL) .AND. 
     &      ((XCO2PX-XCO2DX).GT.EPSL) ) THEN
            NCJ = 0
            DO 10 J = J1X,J2X
              NCJ = NCJ + 1
              XCH4X(NCJ) = REAL(J-1)*DZPC_HT
              K = NZP_HT+1-J
              IZX = IZP_HT(K,J)
              CALL BCTPS( D2V2P_HT,V2P_HT,PPELX,PPEUX,PX,TCTX,TKX,
     &          YCH4X(NCJ),IZX,N_HT )
   10       CONTINUE
            CALL SPLINE( XCH4X,YCH4X,NCJ,Y2CH4X )
            NCI = NCI + 1
            XCO2X(NCI) = XCO2DX
            CALL SPLINT( XCH4X,YCH4X,Y2CH4X,NCJ,ZMCX(2),YCO2X(NCI) )
          ENDIF
        ENDIF
!
!---    CO2 concentration at tabulated values  ---
!
        NCJ = 0
        DO 20 J = J1X,J2X
          NCJ = NCJ + 1
          XCH4X(NCJ) = REAL(J-1)*DZPC_HT
          IZX = IZP_HT(I,J)
          CALL BCTPS( D2V2P_HT,V2P_HT,PPELX,PPEUX,PX,TCTX,TKX,
     &        YCH4X(NCJ),IZX,N_HT )
   20   CONTINUE
        CALL SPLINE( XCH4X,YCH4X,NCJ,Y2CH4X )
        NCI = NCI + 1
        XCO2X(NCI) = REAL(I-1)*DZPC_HT + ZB1X
        CALL SPLINT( XCH4X,YCH4X,Y2CH4X,NCJ,ZMCX(2),YCO2X(NCI) )
   30 CONTINUE
!
!---  Create a cubic spline across the tabulated CO2 concentrations
!     and CO2 concentration at the CH4 diagonal 
!     at the CH4 concentration  ---
!
      CALL SPLINE( XCO2X,YCO2X,NCI,Y2CO2X )
!
!---  Compute the hydrate equilibrium value at the CO2 concentration,
!     using the cubic spline at the CH4 concentration  ---
!
      CALL SPLINT( XCO2X,YCO2X,Y2CO2X,NCI,ZMCX(1),VARX )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCIP_TPS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCTPS( D2V2P_HT,V2P_HT,PPELX,PPEUX,PX,TCTX,TKX,VARX,
     &  IZX,N_HT )
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
!     Bicubic temperature-pressure splines, within the two-phase
!     envelope.
!
!     VARX - interpolated value
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, September 17 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYDT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 D2V2P_HT(LTP_HT,LPP_HT,LCP_HT),V2P_HT(LTP_HT,LPP_HT,LCP_HT)
      INTEGER N_HT(LCP_HT)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BCTPS'
!
!---  Temperature scaled by the cricondenterm temperature  ---
!
      TCTY = TCT_HT(IZX)
      TKZ = TCTY*TKX/TCTX
!
!---  Pressure scaled between lower and upper phase envelope
!     pressures  ---
!
      CALL SPLINT( TLE_HT(1,IZX),PLE_HT(1,IZX),
     &  D2PLE_HT(1,IZX),N_HT(IZX),TKZ,PPELY )
      CALL SPLINT( TUE_HT(1,IZX),PUE_HT(1,IZX),
     &  D2PUE_HT(1,IZX),N_HT(IZX),TKZ,PPEUY )
      PZ = (PPEUY-PPELY)*(PX-PPELX)/(PPEUX-PPELX) + PPELY
      CALL SPLIN2( T2P_HT(1,IZX),P2P_HT(1,IZX),
     &  V2P_HT(1,1,IZX),D2V2P_HT(1,1,IZX),TKZ,PZ,
     &  VARX,LTP_HT,LPP_HT,NTP_HT(IZX),NPP_HT(IZX) )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCTPS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BNDF_G( TX,PX,DFCX,I,J )
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
!     Calculates the component-pair binary diffusion coefficient from
!     the Wilke and Lee theory. pp. 587.
!
!     Liquid molar volume is computed from critical molar volume with
!     the Tyn and Calus method. pp. 53-54.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
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
      REAL*8 A(8)
!
!----------------------Data Statements---------------------------------!
!
      DATA A /1.06036D+0,1.5610D-1,1.9300D-1,4.7635D-1,1.03587D+0,
     &        1.52996D+0,1.76474D+0,3.89411D+0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BNDF_G'
      TKX = TX+TABS
      PBX = PX*1.D-5
      TPIX = 1.15D+0*HCPP(2,I)
      TPJX = 1.15D+0*HCPP(2,J)
      TPX = TKX/SQRT(TPIX*TPJX)
      SIGIX = 1.18D+0*(1.D+3*HCPP(9,I))**3.333D-1
      SIGJX = 1.18D+0*(1.D+3*HCPP(9,J))**3.333D-1
      SIGX = 5.D-1*(SIGIX+SIGJX)
      OMGX = (A(1)/TPX**A(2)) + (A(3)/EXP(A(4)*TPX))
     &  + (A(5)/EXP(A(6)*TPX))  + (A(7)/EXP(A(8)*TPX))
      WTMX = 2.D+0/((1.D+0/HCPP(1,I)) + (1.D+0/HCPP(1,J)))
      DFCX = (3.03D+0 - (9.8D-1/SQRT(WTMX)))*1.D-3*(TKX**1.5D+0)/
     &  (PBX*SQRT(WTMX)*(SIGX**2)*(OMGX))*1.D-4
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BNDF_G group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DENS_A( TX,PX,RHOAX,INDX )
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
!     Density of pure CO2.
!
!     Peng-Robinson compressibility.
!
!     INDX = 1 Gas
!     INDX = 2 Liquid
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 21 September 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FDVH
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
      SUB_LOG(ISUB_LOG) = '/DENS_A'
      TKX = TX + TABS
      TRX = TKX/HCPP(3,1)
      AX = HCPP(13,1)*(1.D+0 + HCPP(15,1)*(1.D+0 - SQRT(TRX)))**2
      BX = HCPP(14,1)
      CAX = AX*PX/((RCU*TKX)**2)
      CBX = BX*PX/(RCU*TKX)
!
!---  Peng-Robinson cubic function for Z  ---
!
      C3X = 1.D+0
      C2X = -(1.D+0-CBX)
      C1X = (CAX - 3.D+0*(CBX**2) - 2.D+0*CBX)
      C0X = -(CAX*CBX - (CBX**2) - (CBX**3))
!
!---  Direct cubic equation solve  ---
!
      CALL NICKALLS( C3X,C2X,C1X,C0X,R1X,R2X,R3X )
!
!---  Gas  ---
!
      IF( INDX.EQ.1 ) THEN
        ZX = MAX(R1X,R2X,R3X)
!
!---  Liquid  ---
!
      ELSEIF( INDX.EQ.2 ) THEN
        ZX = MIN(R1X,R2X,R3X)
!
!---  Ideal gas  ---
!
      ELSE
        ZX = 1.D+0
      ENDIF
!
!---  Gas density, kg/m^3  ---
!
      RHOAX = HCPP(1,1)*PX/(ZX*TKX*RCU)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DENS_A group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DENS_B( TX,PX,XLSX,RHOBX )
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
!     This subroutine calculates the density of NaCl brine as a
!     function of temperature (C), pressure (Pa), and NaCl
!     concentration (molality, mol NaCl/kg H2O).
!
!     Haas Jr., J.L.  1976.  Physical Properties of the Coexisting
!     Phases and Thermochemical Properties of the H2O Component
!     in Boiling NaCl Solutions.  Preliminary Steam Tables for
!     NaCl Solutions.  Geological Survey Bulletin 1421-A.
!
!     Phillips, S.L., H. Ozbek, and L.F. Silvester.  1983.
!     Density of Sodium Chloride Solutions at High Temperatures and
!     Pressures, LBL-16275, Lawrence Berkeley Laboratory, University
!     of California, Berkeley, California.
!
!     Temperature Range: 0 - 350 C
!     Pressure Range:  0.1 - 100 MPa
!     NaCl Concentration Range:  0 - 5 Molal (mol NaCl/kg H2O)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 4 April 2002.
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
      REAL*8 CCX(4),CHX(10),SCX(3),SAX(3)
!
!----------------------Data Statements---------------------------------!
!
      DATA CCX / -3.033405D+0, 10.128163D+0, -8.750567D+0, 2.663107D+0 /
      DATA CHX / -167.219D+0, 448.55D+0, -261.07D+0, -13.644D+0,
     &  13.97D+0, -0.315154D+0, -1.203374D-3, 7.48908D-13,
     &  0.1342489D+0, -3.946963D-3 /
      DATA SCX / -9.9559D+0, 7.0845D+0, 3.9093D+0 /
      DATA SAX / -4.539D-3, -1.638D-4, 2.551D-5 /
      DATA VCX / 3.1975D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DENS_B'
!
!---  Convert pressure to bar and mass fraction to molality  ---
!
      PBX = 1.D-5*PX
      GLSX = 1.D+3*XLSX/(WTMS*(1.D+0-XLSX))
!
!---  Restrict temperature to above 0.01 C  ---
!
      TY = MAX( TX,1.D-2 )
!
!---  Compressed or vapor-saturated density of pure water using the
!     ASME formulations  ---
!
      CALL SP_W( TY,PSWX )
      ISRX = 1
      PWX = MAX( PX,PSWX )
      CALL DENS_W( TY,PWX,RHOLWX,RHOX,ISRX )
!
!---  Convert density units to gm/cm^3  ---
!
      RHOLWX = 1.D-3*RHOLWX
      VOX = 1.D+0/RHOLWX
!
!---  Limiting apparent molal volume (cm^3/mol) of NaCl in solution as
!     the concentration goes to zero  ---
!
      PHIPX = CHX(1) + CHX(2)*VOX + CHX(3)*(VOX**2)
!
!---  Apparent molal volume (cm^3/mol) of NaCl in solution  ---
!
      PHIX = PHIPX + (CHX(4)+CHX(5)*VOX)*((VOX/(VCX-VOX))**2)*SQRT(GLSX)
!
!---  Brine density (gm/cm^3)  ---
!
      RHOBX = (1.D+3 + GLSX*WTMS)/(1.D+3*VOX + GLSX*PHIX)
!
!---  Compressed or vapor-saturated density of pure water (gm/cm^3)
!     using the Phillips et al. formulations  ---
!
!      PBWX = 1.D-5*PWX
!      SXX = SCX(1) + SCX(2)*EXP(SAX(2)*TY) +
!     &  SCX(3)*EXP(SAX(3)*PBWX)
!      RHOLWPX = CCX(1) + CCX(2)*SXX + CCX(3)*(SXX**2) + CCX(4)*(SXX**3)
!
!---  Density of NaCl brine (gm/cm^3) using the Phillips et al.
!     formulations  ---
!
!      SXX = SCX(1)*EXP(SAX(1)*GLSX) + SCX(2)*EXP(SAX(2)*TY) +
!     &  SCX(3)*EXP(SAX(3)*PBX)
!      RHOBPX = CCX(1) + CCX(2)*SXX + CCX(3)*(SXX**2) +  CCX(4)*(SXX**3)
!
!---  Apparent molal volume of NaCl in solution  ---
!
!      SPHIX = 0.D+0
!      IF( GLSX/EPSL.GT.EPSL ) SPHIX = (RHOLWPX*(1.D+3 + GLSX*WTMS) -
!     &  RHOBPX*1.D+3)/(GLSX*RHOLWPX*RHOBPX)
!
!---  Normalize brine density (gm/cm^3) to ASME formulation  ---
!
!      RHOBX = (1.D+3 + GLSX*WTMS)/((1.D+3/RHOLWX) + GLSX*SPHIX)
!      RHOBX = RHOBPX + (RHOLWX - RHOLWPX)*EXP(-GLSX)
!
!---  Convert density units to kg/m^3  ---
!
      RHOBX = 1.D+3*RHOBX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DENS_B group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DENS_I( TX,RHOIX )
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
!     Calculate the density of ice as a function of temperature using a
!     using polynomial fit of data from ASHRAE (1977).
!
!     ASHRAE, 1977; Psychrometric Tables given by ASHRAE Handbook 1977
!     Fundamentals, American Society of Heating, Refrigerating and
!     Air-Conditioning Engineering, Inc., 345 East 47th Street, New
!     York, N.Y.
!
!     The temperature is limited in this subroutine to the following
!     values:  -100.0 C < T <= 0.0 C
!
!----------------------Authors-----------------------------------------!
!
!     Written by WE Nichols, PNNL, September 1994.
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
      REAL*8 A(3)
!
!----------------------Data Statements---------------------------------!
!
      DATA A / 1.06472D-3,1.23173D-8,3.0203D-10 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DENS_I'
      TKX = TX+TABS
!
!---  Polynomial fit for specific volume  ---
!
      SPVI = A(1) + A(2)*TKX + A(3)*(TKX**2)
      RHOIX = 1.D+0/SPVI
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DENS_I group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DENS_L( TX,RHOBX,XLAX,XLOX,XLNX,RHOLX )
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
!     This subroutine calculates the density of CO2 gas dissolved
!     in NaCl brine as a function of the brine density and dissolved
!     gas mass fraction.
!
!     Alendal, G., and H. Drange.  2001.  "Two-phase, near-field,
!     modeling of purposefully released CO2 in the ocean."  Journal
!     of Geophysical Research, 106(C1):1085-1096.
!
!     Anderson, G.M., and D.A. Crerar.  1992.  Thermodynamics in
!     Geochemistry: The Equilibrium Model, Oxford University Press.
!
!     Variable definitions
!
!     PMV_A - molar volume of CO2 (m^3/kmol)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 April 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NAPL
      USE HYDT
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------TYPE DECLARATIONS-------------------------------!
!
      REAL*8 PMV_A
      REAL*8 PMV_C(5)
      REAL*8 CX(LHF_HT)
!
!----------------------DATA STATEMENTS---------------------------------!
!
!      DATA PMV_A / 34.D-3 /
!
!----------------------Data Statements---------------------------------!
!
      DATA PMV_C / 37.36D-3,-7.109D-5,-3.812D-8,3.296D-9,-3.702D-12 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DENS_L'
!
!---  Restrict temperature to above 0.01 C  ---
!
      TY = MAX( TX,1.D-2 )
!
!---  Partial molar volume of CO2 (m^3/kmol), formulation of
!     Anderson et al. (1992)  ---
!
      PMV_A = 0.D+0
      DO 10 M = 1,5
        PMV_A = PMV_A + PMV_C(M)*(TY**(M-1))
   10 CONTINUE
      CX(1) = PMV_A*RHOBX*XLAX/WTMA
!
!---  Partial molar volume of CH4 (m^3/kmol), at 25 C  ---
!
      PMV_O = 37.D-3
      CX(2) = PMV_O*RHOBX*XLOX/WTMO
!
!---  Partial molar volume of N2 (m^3/kmol), at 25 C  ---
!
      PMV_N = 40.D-3
      CX(3) = PMV_N*RHOBX*XLNX/WTMN
!
!---  Poynting correction  ---
!
      PCX = 1.D+0
      DO 20 NF = 1,NHF_HT
        PCX = PCX + CX(NF)
   20 CONTINUE
      PCX = PCX - XLAX - XLOX - XLNX
      RHOLX = RHOBX/PCX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DENS_L group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DENS_O( TX,PX,RHOOX,INDX )
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
!     Density of pure CH4.
!
!     Peng-Robinson compressibility.
!
!     INDX = 1 Gas
!     INDX = 2 Liquid
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 21 September 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FDVH
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
      SUB_LOG(ISUB_LOG) = '/DENS_O'
      TKX = TX + TABS
      TRX = TKX/HCPP(3,2)
      AX = HCPP(13,2)*(1.D+0 + HCPP(15,2)*(1.D+0 - SQRT(TRX)))**2
      BX = HCPP(14,2)
      CAX = AX*PX/((RCU*TKX)**2)
      CBX = BX*PX/(RCU*TKX)
!
!---  Peng-Robinson cubic function for Z  ---
!
      C3X = 1.D+0
      C2X = -(1.D+0-CBX)
      C1X = (CAX - 3.D+0*(CBX**2) - 2.D+0*CBX)
      C0X = -(CAX*CBX - (CBX**2) - (CBX**3))
!
!---  Direct cubic equation solve  ---
!
      CALL NICKALLS( C3X,C2X,C1X,C0X,R1X,R2X,R3X )
!
!---  Gas  ---
!
      IF( INDX.EQ.1 ) THEN
        ZX = MAX(R1X,R2X,R3X)
!
!---  Liquid  ---
!
      ELSEIF( INDX.EQ.2 ) THEN
        ZX = MIN(R1X,R2X,R3X)
!
!---  Ideal gas  ---
!
      ELSE
        ZX = 1.D+0
      ENDIF
!
!---  Gas density, kg/m^3  ---
!
      RHOOX = HCPP(1,2)*PX/(ZX*TKX*RCU)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DENS_O group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DENS_N( TX,PX,RHONX,INDX )
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
!     Density of pure N2.
!
!     Peng-Robinson compressibility.
!
!     INDX = 1 Gas
!     INDX = 2 Liquid
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 21 September 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FDVH
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
      SUB_LOG(ISUB_LOG) = '/DENS_N'
      TKX = TX + TABS
      TRX = TKX/HCPP(3,3)
      AX = HCPP(13,3)*(1.D+0 + HCPP(15,3)*(1.D+0 - SQRT(TRX)))**2
      BX = HCPP(14,3)
      CAX = AX*PX/((RCU*TKX)**2)
      CBX = BX*PX/(RCU*TKX)
!
!---  Peng-Robinson cubic function for Z  ---
!
      C3X = 1.D+0
      C2X = -(1.D+0-CBX)
      C1X = (CAX - 3.D+0*(CBX**2) - 2.D+0*CBX)
      C0X = -(CAX*CBX - (CBX**2) - (CBX**3))
!
!---  Direct cubic equation solve  ---
!
      CALL NICKALLS( C3X,C2X,C1X,C0X,R1X,R2X,R3X )
!
!---  Gas  ---
!
      IF( INDX.EQ.1 ) THEN
        ZX = MAX(R1X,R2X,R3X)
!
!---  Liquid  ---
!
      ELSEIF( INDX.EQ.2 ) THEN
        ZX = MIN(R1X,R2X,R3X)
!
!---  Ideal gas  ---
!
      ELSE
        ZX = 1.D+0
      ENDIF
!
!---  Gas density, kg/m^3  ---
!
      RHONX = HCPP(1,3)*PX/(ZX*TKX*RCU)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DENS_N group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DENS_S( TX,PX,RHOSX )
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
!     Density of precipitated NaCl.
!
!     Battistelli, A., C. Claudio, and K. Pruess.  1997.  "The simulator
!     TOUGH2/EWASG for modelling geothermal reservoirs with brines and
!     noncondensible gas."  Geothermics, 26(4): 437-464.
!
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 1 May 2002.
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DENS_S'
!
!---  Formulation of Battistelli et al.  ---
!
      RHOSX = 2.165D+3*EXP(-1.2D-4*TX + 4.D-11*PX)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DENS_S group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DENS_W( TX,PX,RHOLX,RHOVX,ISRX )
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
!     Density (kg/m^3) of pure water as a function of temperature and
!     pressure.
!
!     Pressure Range: 0 - 100 MPa (1000 Bar)
!     Temperature Range:  273.16 K (0.01 C) to 1073.15K (800 C)
!
!     Meyer, C.A., R.B. McClintock, G.J. Silvestri, and R.C. Spencer
!     1993.  ASME Steam Tables, The American Society of Mechanical
!     Engineers, New York.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 March 2002.
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
      REAL*8 CAX(23),SAX(12),CBX(31),SBX(5),LFCX(3),SLX
      REAL*8 TX,PX,RHOLX,RHOVX
      INTEGER ISRX,INX(8),IZX(8,3),ITX(8),IXX(8,2)
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Data Statements---------------------------------!
!
      DATA CAX / 6.824687741D+3, -5.422063673D+2, -2.096666205D+4,
     &  3.941286787D+4, -6.733277739D+4, 9.902381028D+4,
     &  -1.093911774D+5, 8.590841667D+4, -4.511168742D+4,
     &  1.418138926D+4, -2.017271113D+3, 7.982692717D+0,
     &  -2.616571843D-2, 1.522411790D-3, 2.284279054D-2,
     &  2.421647003D+2, 1.269716088D-10, 2.074838328D-7,
     &  2.174020350D-8, 1.105710498D-9, 1.293441934D+1,
     &  1.308119072D-5, 6.047626338D-14 /
      DATA SAX / 8.438375405D-1, 5.362162162D-4, 1.720000000D+0,
     &  7.342278489D-2, 4.975858870D-2, 6.537154300D-1,
     &  1.150000000D-6, 1.150800000D-5, 1.418800000D-1,
     &  7.002753165D+0, 2.995284926D-4, 2.040000000D-1 /
      DATA CBX / 1.683599274D+1, 2.856067796D+1, -5.438923329D+1,
     &  4.330662834D-1, -6.547711697D-1, 8.565182058D-2,
     &  6.670375918D-2, 1.388983801D+0, 8.390104328D-2,
     &  2.614670893D-2, -3.373439453D-2, 4.520918904D-1,
     &  1.069036614D-1, -5.975336707D-1, -8.847535804D-2,
     &  5.958051609D-1, -5.159303373D-1, 2.075021122D-1,
     &  1.190610271D-1, -9.867174132D-2, 1.683998803D-1,
     &  -5.809438001D-2, 6.552390126D-3, 5.710218649D-4,
     &  1.936587558D+2, -1.388522425D+3, 4.126607219D+3,
     &  -6.508211677D+3, 5.745984054D+3, -2.693088365D+3,
     &  5.235718623D+2/
      DATA SBX / 7.633333333D-1, 4.006073948D-1, 8.636081627D-2,
     &  -8.532322921D-1, 3.460208861D-1 /
      DATA LFCX / 1.574373327D+1, -3.417061978D+1, 1.931380707D+1 /
      DATA SLX / 4.260321148D+0 /
      DATA INX / 2, 3, 2, 2, 3, 2, 2, 2 /
      DATA IZX / 13, 18, 18, 25, 32, 12, 24, 24,
     &  3, 2, 10, 14, 28, 11, 18, 14,
     &  0, 1, 0, 0, 24, 0, 0, 0 /
      DATA ITX / 0, 0, 0, 0, 0, 1, 1, 2 /
      DATA IXX / 0, 0, 0, 0, 0, 14, 19, 54,
     &  0, 0, 0, 0, 0, 0, 0, 27 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DENS_W'
!
!---  Restrict temperature to above 0.01 C and below 800 C  ---
!
      TY = MIN( MAX( TX,1.D-2 ),8.D+2 )
!
!---  Restrict pressure to above 0.0 MPa and below 100 MPa  ---
!
      PY = MIN( MAX( PX,0.D+0 ),1.D+8 )
!
!---  Zero pressure condition  ---
!
      IF( PY/EPSL.LT.EPSL ) THEN
        RHOLX = 0.D+0
        RHOVX = 0.D+0
        GOTO 100
      ENDIF
!
!---  Reduced temperature and pressure  ---
!
      RHOLX = 0.D+0
      RHOVX = 0.D+0
      THETAX = (TY+TABS)/TCRW
      BETAX = PY/PCRW
!
!---  L-function  ---
!
      BETALX = LFCX(1) + LFCX(2)*THETAX + LFCX(3)*(THETAX**2)
!
!---  Subregions 1 or 6  ---
!
      IF( ISRX.EQ.1 .OR. ISRX.EQ.6 ) THEN
        CYX = 1.D+0 - SAX(1)*(THETAX**2) - SAX(2)/(THETAX**6)
        CZX = CYX + SQRT((SAX(3)*(CYX**2)) - (2.D+0*SAX(4)*THETAX) +
     &    (2.D+0*SAX(5)*BETAX))
        RVX = CAX(12)*SAX(5)*(CZX**(-5.D+0/17.D+0))
        RVX = RVX +  (CAX(13) + CAX(14)*THETAX + CAX(15)*(THETAX**2) +
     &    CAX(16)*((SAX(6)-THETAX)**10) + CAX(17)/(SAX(7)+(THETAX**19)))
        RVX = RVX - (CAX(18) + 2.D+0*CAX(19)*BETAX +
     &    3.D+0*CAX(20)*(BETAX**2))/(SAX(8)+(THETAX**11))
        RVX = RVX - CAX(21)*(THETAX**18)*(SAX(9)+(THETAX**2))*
     &    (-3.D+0/((SAX(10)+BETAX)**4))
        RVX = RVX + 3.D+0*CAX(22)*(SAX(12)-THETAX)*(BETAX**2)
        RVX = RVX + 4.D+0*CAX(23)*(BETAX**3)/(THETAX**20)
        RHOLX = (1.D+3*WTMW)/(RVX*VCRW)
      ENDIF
!
!---  Subregions 2 or 6  ---
!
      IF( ISRX.EQ.2 .OR. ISRX.EQ.6 ) THEN
        CXX = EXP(SBX(1)*(1.D+0-THETAX))
        RVX = SLX*THETAX/BETAX
        ICBX = 6
        DO 30 I = 1,5
          RVAX = 0.D+0
          DO 20 J = 1,INX(I)
            ICBX = ICBX + 1
            RVAX = RVAX + CBX(ICBX)*(CXX**IZX(I,J))
   20     CONTINUE
          REALX = REAL(I)
          RVX = RVX - REALX*(BETAX**(I-1))*RVAX
   30   CONTINUE
        ICBX = 18
        ISBX = 1
        DO 60 I = 6,8
          RVAX = 0.D+0
          DO 40 J = 1,INX(I)
            ICBX = ICBX + 1
            RVAX = RVAX + CBX(ICBX)*(CXX**IZX(I,J))
   40     CONTINUE
          RVBX = 0.D+0
          DO 50 J = 1,ITX(I)
            ISBX = ISBX + 1
            RVBX = RVBX + SBX(ISBX)*(CXX**IXX(I,J))
   50     CONTINUE
          REALX = REAL(I)
          RVX = RVX - ((REALX-2.D+0)*(BETAX**(1-I))*RVAX)/
     &      (((BETAX**(2-I))+RVBX)**2)
   60   CONTINUE
        ICBX = 24
        RVAX = 0.D+0
        DO 70 I = 0,6
          ICBX = ICBX + 1
          RVAX = RVAX + CBX(ICBX)*(CXX**I)
   70   CONTINUE
        RVX = RVX + 1.1D+1*((BETAX/BETALX)**10)*RVAX
        RHOVX = (1.D+3*WTMW)/(RVX*VCRW)
      ENDIF
!
!---  Subregions 3 or 5  ---
!
      IF( ISRX.EQ.3 .OR. ISRX.EQ.5 ) THEN
        INDX = 4
        IMSGX = 0
        NMSGX = 0
        SUBLOGX = 'DENS_W'
        CHMSGX = 'Steam Table Subregion 3/5: Exceeds Temperature ' // 
     &    'Range: '
        RLMSGX = TX
        CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
      ENDIF
!
!---  Subregions 4 or 5  ---
!
      IF( ISRX.EQ.4 .OR. ISRX.EQ.5 ) THEN
        INDX = 4
        IMSGX = 0
        NMSGX = 0
        SUBLOGX = 'DENS_W'
        CHMSGX = 'Steam Table Subregion 4/5: Exceeds Temperature ' // 
     &    'Range: '
        RLMSGX = TX
        CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
      ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DENS_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DIFC_G( TX,PX,DFGAX,DFGOX,DFGNX,DFGWX,
     &  XMGAX,XMGOX,XMGNX,XMGWX )
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
!     Component gas diffusion coefficient, using Blanc's law
!     and binary diffusion coefficients, pp. 597.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 21 September 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYDT
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 DFCX(LHF_HT+1,LHF_HT+1),DFGX(LHF_HT+1),XMGX(LHF_HT+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DIFC_G'
!
!---  Translate gas mole fractions  ---
!
      XMGX(1) = XMGAX
      XMGX(2) = XMGOX
      XMGX(3) = XMGNX
      XMGX(4) = XMGWX
!
!---  Binary diffusion coefficients  ---
!
      DO 20 J = 1,NHF_HT+1
        DO 10 I = 1,NHF_HT+1
          IF( I.EQ.J ) THEN
            DFCX(I,J) = 1.D+0
          ELSEIF( I.GT.J ) THEN
            CALL BNDF_G( TX,PX,DFCX(I,J),I,J )
          ELSE
            DFCX(I,J) = DFCX(J,I)
          ENDIF
   10   CONTINUE
   20 CONTINUE
!
!---  Gas component diffusion coefficients  ---
!
      DO 40 I = 1,NHF_HT+1
        IF( (1.D+0-XMGX(I)).LT.EPSL ) THEN
          DFGX(I) = 0.D+0
          DO 30 J = 1,NHF_HT+1
            IF( I.NE.J ) DFGX(I) = DFGX(I) + 1.D+0/DFCX(I,J)
   30     CONTINUE
          DFGX(I) = REAL(NHF_HT)/(DFGX(I)+SMALL)
        ELSE
          DFGX(I) = 0.D+0
          DO 32 J = 1,NHF_HT+1
            IF( I.NE.J ) DFGX(I) = DFGX(I) + XMGX(J)/DFCX(I,J)
   32     CONTINUE
          DFGX(I) = (1.D+0-XMGX(I))/(DFGX(I)+SMALL)
        ENDIF
   40 CONTINUE
!
!---  Translate gas diffusion coefficients  ---
!
      DFGAX = DFGX(1)
      DFGOX = DFGX(2)
      DFGNX = DFGX(3)
      DFGWX = DFGX(4)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DIFC_G group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DIFC_LA( TX,VISLX,VISGAX,DFLAX )
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
!     Calculates the diffusion coefficient (m^2/s) for CO2 through NaCl
!     aqueous solutions.
!
!     Renner, T.A. 1988. "Measurement and Correlation of Diffusion
!     Coefficients for CO2 and Rich-Gas Applications."  SPE Reservoir
!     Engineering, pp:517-523.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 1 May 2002.
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DIFC_LA'
!
!---  Restrict temperature to above 0.01 C  ---
!
      TY = MAX( TX,1.D-2 )
!
!---  Convert viscosities to centipoise  ---
!
      DFLAX = 6.391D+3*(((VISLX*1.D+3)**(-1.584D-1))*
     &  ((VISGAX*1.D+3)**(6.911D+0)))
!
!---  Limit diffusion coefficient to the upper limit
!       reported in the reference  ---
!
      DFLAX = MIN( DFLAX,10.D-9 )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DIFC_LA group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DIFC_LN( TX,VISLX,DFLNX )
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
!     Calculates binary liquid diffusion coefficient for N2 in liquid
!     water, using the Wilke and Chang estimation method; pp. 598.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 24 September 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FDVH
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
      SUB_LOG(ISUB_LOG) = '/DIFC_LN'
!
!---  Restrict temperature to above 0.01 C  ---
!
      TY = MAX( TX,1.D-2 )
      VBO = (1.D+3*HCPP(9,3))
      DFLNX = 7.4D-15*((2.6D+0*WTMW)**5.D-1)*(TY+TABS)/
     &  (VISLX*(VBO**6.D-1))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DIFC_LN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DIFC_LO( TX,VISLX,DFLOX )
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
!     Calculates binary liquid diffusion coefficient for CH4 in liquid
!     water, using the Wilke and Chang estimation method; pp. 598.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 1 May 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FDVH
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
      SUB_LOG(ISUB_LOG) = '/DIFC_LO'
!
!---  Restrict temperature to above 0.01 C  ---
!
      TY = MAX( TX,1.D-2 )
      VBO = (1.D+3*HCPP(9,2))
      DFLOX = 7.4D-15*((2.6D+0*WTMW)**5.D-1)*(TY+TABS)/
     &  (VISLX*(VBO**6.D-1))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DIFC_LO group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DIFC_LS( TX,XLSX,VISLX,DFLSX )
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
!     Calculates the diffusion coefficient (m^2/s) for NaCl in
!     aqueous solutions, following the method of Nernst-Haskell for
!     dilute solutions and the method of Gordon for concentrated
!     solutions; where the mean ionic activity is computed according
!     to the method of Bromley.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids. pp. 620-621.
!     McGraw-Hill, New York, New York
!
!     Bromley, L.A.  1973.  "Thermodynamic properties of strong
!     electrolytes in aqueous solutions."  AIChE Journal, 19(2):313-320
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 2 May 2002.
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
!----------------------Data Statements---------------------------------!
!
      DATA CBX / 0.0547D+0 /
      DATA TRX / 25.D+0 /
      DATA VISWRX / 0.8904339807D-3 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DIFC_LS'
!
!---  Restrict temperature to above 0.01 C  ---
!
      TY = MAX( TX,1.D-2 )
!
!---  Convert halite mass fraction to molaity   ---
!
      GLSX = 1.D+3*XLSX/(WTMS*(1.D+0-XLSX))
      TKX = TY+TABS
!
!---  Partial derivative of the natural logarithm  of the mean ionic
!     activity with respect to the molaity at 298 K   ---
!
      IF( GLSX.GT.EPSL ) THEN
        DLNGX = (-0.2555D+0/(SQRT(GLSX)*(1.D+0+SQRT(GLSX))) +
     &    0.2555D+0/((1.D+0+SQRT(GLSX))**2) +
     &    (6.D-2 + 6.D-1*CBX)/((1.D+0+1.5D+0*GLSX)**2) -
     &    3.D+0*(6.D-2 + 6.D-1*CBX)*GLSX/((1.D+0+1.5D+0*GLSX)**3) +
     &    CBX)*2.302585D+0
      ELSE
        DLNGX = -5.883D-6
      ENDIF
!
!---  Diffusion coefficient (m^2/s) for dilute NaCl aqueous solutions
!     at 298 K   ---
!
      DFLSX = 1.61D-9
!
!---  Viscosity of brine as a function of NaCl mass fraction
!     at 298.15 K   ---
!
      CALL VISC_B( TRX,XLSX,VISWRX,VISBRX )
!
!---  Diffusion coefficient for concentrated NaCl aqueous solutions
!     at 298.15 K   ---
!
      DFLSX = DFLSX*(VISWRX/VISBRX)*(1.D+0 + GLSX*(DLNGX))
!
!---  Viscosity of brine as a function of NaCl mass fraction  ---
!
      CALL VISC_B( TKX,XLSX,VISLX,VISBX )
!
!---  Correct diffusion coefficient for temperature  ---
!
      DFLSX = DFLSX*(TKX/2.9815D+2)*(VISBRX/VISBX)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DIFC_LS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DIFC_N( DFNCX,TX,VISNX,WTMNX )
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
!     Calculates the diffusion coefficient for a gas component
!     in NAPL, using the Wilke and Chang estimation method; pp. 598.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 August 2009
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
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
      REAL*8 DFNCX(LNHC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DIFC_N'
!
!---  NAPL diffusion coefficients  ---
!
      DO 100 IGC = 1,NGC
        VBCX = 2.85D-1*HCPP(5,IGC)**1.048D+0
        DFNCX(IGC) = 7.4D-15*((2.6D+0*WTMNX)**5.D-1)*(TX+TABS)/
     &    (VISNX*(VBCX**6.D-1))
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DIFC_N group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ENTH_CH4( TX,PX,HGAX,HNAX,UGAX,UNAX )
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
!     Enthalpy and internal energy of CH4.
!
!     Setzmann, U. and W. Wagner.  1991.  A New Equation of State and
!     Tables of Thermodynamic Properties for Methane Covering the
!     Range from the Melting Line to 625 K at Pressures up to 1000 MPa.
!     J. Phys. Chem. Ref. Data, 20(6):1061-1155.
! 
!     I_VX = 1 subcritical gas
!     I_VX = 2 subcritical liquid
!     I_VX = 3 supercritical gas
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 April 2002.
!     Last Modified by MD White, PNNL, 25 April 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER I_PX(2),I_TX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ENTH_CH4'
!
!---  Table lookup and bilinear interpolation  ---
!
      NF = 2
      PTX = MAX( PX,P_TA(1,NF)*1.D+6 )
      CALL ITL_A( PTX,TX,I_PX,I_TX,I_VX,NF)
      CALL PTL_A( PTX,TX,H_TA,HAX,I_PX,I_TX,NF )
      CALL PTL_A( PTX,TX,U_TA,UAX,I_PX,I_TX,NF )
!
!---  Subcritical gas conditions, set liquid enthalpy and internal
!     energy to saturated liquid conditions at the same temperature ---
!
      IF( I_VX.EQ.1 ) THEN
        HGAX = HAX
        UGAX = UAX
!
!---    Find saturation indices  ---
!
        TKX = TX + TABS
        IF( TKX.LT.T_LV(1,NF) ) TKX = T_LV(1,NF)
        IF( TKX.GT.T_LV(I_LV(NF),NF) ) TKX = T_LV(I_LV(NF),NF)
        ITLX = 1
        ITUX = I_LV(NF)
   10   IF( ITUX-ITLX.GT.1 ) THEN
          ITM = (ITLX+ITUX)/2
          IF( (T_LV(I_LV(NF),NF).GT.T_LV(1,NF)).EQV.
     &      (TKX.GT.T_LV(ITM,NF)) ) THEN
            ITLX = ITM
          ELSE
            ITUX = ITM
          ENDIF
          GOTO 10
        ENDIF
        TKLX = T_LV(ITLX,NF)
        TKUX = T_LV(ITLX+1,NF)
        HLX = HL_LV(ITLX,NF)
        HUX = HL_LV(ITLX+1,NF)
        HNAX = ((TKX-TKLX)/(TKUX-TKLX))*(HUX-HLX) + HLX
        ULX = UL_LV(ITLX,NF)
        UUX = UL_LV(ITLX+1,NF)
        UNAX = ((TKX-TKLX)/(TKUX-TKLX))*(UUX-ULX) + ULX
!
!---  Subcritical liquid conditions, set gas enthalpy and internal
!     energy to saturated gas conditions at the same temperature ---
!
      ELSEIF( I_VX.EQ.2 ) THEN
        HNAX = HAX
        UNAX = UAX
!
!---    Find saturation indices  ---
!
        TKX = TX + TABS
        IF( TKX.LT.T_LV(1,NF) ) TKX = T_LV(1,NF)
        IF( TKX.GT.T_LV(I_LV(NF),NF) ) TKX = T_LV(I_LV(NF),NF)
        ITLX = 1
        ITUX = I_LV(NF)
   20   IF( ITUX-ITLX.GT.1 ) THEN
          ITM = (ITLX+ITUX)/2
          IF( (T_LV(I_LV(NF),NF).GT.T_LV(1,NF)).EQV.
     &      (TKX.GT.T_LV(ITM,NF)) ) THEN
            ITLX = ITM
          ELSE
            ITUX = ITM
          ENDIF
          GOTO 20
        ENDIF
        TKLX = T_LV(ITLX,NF)
        TKUX = T_LV(ITLX+1,NF)
        HLX = HV_LV(ITLX,NF)
        HUX = HV_LV(ITLX+1,NF)
        HGAX = ((TKX-TKLX)/(TKUX-TKLX))*(HUX-HLX) + HLX
        ULX = UV_LV(ITLX,NF)
        UUX = UV_LV(ITLX+1,NF)
        UGAX = ((TKX-TKLX)/(TKUX-TKLX))*(UUX-ULX) + ULX
!
!---  Supercritical conditions, set gas and liquid enthalpy and
!     and internal energy to the gas value  ---
!
      ELSEIF( I_VX.EQ.3 ) THEN
        HGAX = HAX
        UGAX = UAX
        HNAX = HAX
        UNAX = UAX
      ENDIF    
!
!---  Convert enthalpy and internal energy units to J/kg  ---
!
      HGAX = 1.D+3*HGAX
      HNAX = 1.D+3*HNAX
      UGAX = 1.D+3*UGAX
      UNAX = 1.D+3*UNAX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTH_CH4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ENTH_CO2( TX,PX,HGAX,HNAX,UGAX,UNAX )
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
!     Enthalpy and internal energy of CO2.
!
!     Span, R., and W. Wagner.  1996.  A New Equation of State for
!     Carbon Dioxide Covering the Fluid Region from the Triple-Point
!     to 1100 K at Pressures up to 800 MPa.  J. Phys. Chem. Ref. Data
!     25(6):1509-1588.
! 
!     I_VX = 1 subcritical gas
!     I_VX = 2 subcritical liquid
!     I_VX = 3 supercritical gas
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 April 2002.
!     Last Modified by MD White, PNNL, 25 April 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
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
      INTEGER I_PX(2),I_TX(2)
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ENTH_CO2'
!
!---  Table lookup and bilinear interpolation  ---
!
      NF = 1
      PTX = MAX( PX,P_TA(1,NF)*1.D+6 )
      CALL ITL_A( PTX,TX,I_PX,I_TX,I_VX,NF)
      CALL PTL_A( PTX,TX,H_TA,HAX,I_PX,I_TX,NF )
      CALL PTL_A( PTX,TX,U_TA,UAX,I_PX,I_TX,NF )
!
!---  Subcritical gas conditions, set liquid enthalpy and internal
!     energy to saturated liquid conditions at the same temperature ---
!
      IF( I_VX.EQ.1 ) THEN
        HGAX = HAX
        UGAX = UAX
!
!---    Find saturation indices  ---
!
        TKX = TX + TABS
        IF( TKX.LT.T_LV(1,NF) ) TKX = T_LV(1,NF)
        IF( TKX.GT.T_LV(I_LV(NF),NF) ) TKX = T_LV(I_LV(NF),NF)
        ITLX = 1
        ITUX = I_LV(NF)
   10   IF( ITUX-ITLX.GT.1 ) THEN
          ITM = (ITLX+ITUX)/2
          IF( (T_LV(I_LV(NF),NF).GT.T_LV(1,NF)).EQV.
     &      (TKX.GT.T_LV(ITM,NF)) ) THEN
            ITLX = ITM
          ELSE
            ITUX = ITM
          ENDIF
          GOTO 10
        ENDIF
        TKLX = T_LV(ITLX,NF)
        TKUX = T_LV(ITLX+1,NF)
        HLX = HL_LV(ITLX,NF)
        HUX = HL_LV(ITLX+1,NF)
        HNAX = ((TKX-TKLX)/(TKUX-TKLX))*(HUX-HLX) + HLX
        ULX = UL_LV(ITLX,NF)
        UUX = UL_LV(ITLX+1,NF)
        UNAX = ((TKX-TKLX)/(TKUX-TKLX))*(UUX-ULX) + ULX
!
!---  Subcritical liquid conditions, set gas enthalpy and internal
!     energy to saturated gas conditions at the same temperature ---
!
      ELSEIF( I_VX.EQ.2 ) THEN
        HNAX = HAX
        UNAX = UAX
!
!---    Find saturation indices  ---
!
        TKX = TX + TABS
        IF( TKX.LT.T_LV(1,NF) ) TKX = T_LV(1,NF)
        IF( TKX.GT.T_LV(I_LV(NF),NF) ) TKX = T_LV(I_LV(NF),NF)
        ITLX = 1
        ITUX = I_LV(NF)
   20   IF( ITUX-ITLX.GT.1 ) THEN
          ITM = (ITLX+ITUX)/2
          IF( (T_LV(I_LV(NF),NF).GT.T_LV(1,NF)).EQV.
     &      (TKX.GT.T_LV(ITM,NF)) ) THEN
            ITLX = ITM
          ELSE
            ITUX = ITM
          ENDIF
          GOTO 20
        ENDIF
        TKLX = T_LV(ITLX,NF)
        TKUX = T_LV(ITLX+1,NF)
        HLX = HV_LV(ITLX,NF)
        HUX = HV_LV(ITLX+1,NF)
        HGAX = ((TKX-TKLX)/(TKUX-TKLX))*(HUX-HLX) + HLX
        ULX = UV_LV(ITLX,NF)
        UUX = UV_LV(ITLX+1,NF)
        UGAX = ((TKX-TKLX)/(TKUX-TKLX))*(UUX-ULX) + ULX
!
!---  Supercritical conditions, set gas and liquid enthalpy and
!     and internal energy to the gas value  ---
!
      ELSEIF( I_VX.EQ.3 ) THEN
        HGAX = HAX
        UGAX = UAX
        HNAX = HAX
        UNAX = UAX
!
!---  Subcritical gas conditions below the triple-point  ---
!
      ELSE
        INDX = 4
        RLMSGX = 0.D+0
        IMSGX = 0
        NMSGX = N
        SUBLOGX = 'ENTH_CO2'
        CHMSGX = 'Solid Phase for CO2: T = '
        RLMSGX = TX
        CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
      ENDIF    
!
!---  Convert enthalpy and internal energy units to J/kg  ---
!
      HGAX = 1.D+3*HGAX
      HNAX = 1.D+3*HNAX
      UGAX = 1.D+3*UGAX
      UNAX = 1.D+3*UNAX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTH_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ENTH_B( TX,XLSX,HLWX,HBX )
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
!     This subroutine calculates the enthalpy of NaCl solutions
!     as a function of temperature and NaCl concentration.
!
!     tx - temperature, C
!     xlsx - mass fraction of aqueous NaCl
!     tkbx - thermal conductivity of NaCl brine, W/m K
!
!     Michaelides, E.E.  1981. "Thermodynamic properties of geothermal
!     fluids."  Geothermal Resources Council, Transactions 5:361-364.
!
!     Gudmundsson, J.S., and H. Thrainsson.  1989.  "Power potential of
!     two-phase geothermal wells."  Geothermics 18(3):357-366.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 8 April 2002.
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
      REAL*8 SAX(12),SBX(3),SCX(4),CAX(6),CBX(15)
!
!----------------------Data Statements---------------------------------!
!
      DATA SAX / 9633.6D+0, -4080.0D+0, 286.49D+0, 166.58D+0,
     &  68.577D+0, -4.6856D+0, -0.90963D+0, -0.36524D+0,
     &  0.249667D-1, 0.17965D-2, 0.71924D-3, -0.4900D-4 /
      DATA SBX / -0.83624D-3, 0.16792D+0, -25.9293D+0 /
      DATA SCX / 0.12453D-4, -0.45137D-2, 4.81155D+0, -29.578D+0 /
      DATA CAX / 25.19D+0, 0.1973D+0, -6.0114D-04, 8.81505D-7,
     &  -4.76500D-10, -1.923188214D+5 /
      DATA CBX / -104.51D+0, 81.086D+0, -308.22D+0, -1.6952D+0,
     &  -16.65D+0, -8.6385D+0, 0.010618D+0, 0.029634D+0, 0.61366D+0,
     &  -2.4977D+0, -1.9876D-05, -0.00032121D+0, 0.0022773D+0,
     &  -0.015262D+0, 0.081865D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ENTH_B'
!
!---  Convert mass fraction to molality and weight percent  ---
!
      TKX = TX + TABS
      GLSX = 1.D+3*XLSX/(WTMS*(1.D+0-XLSX))
      YLSX = XLSX*1.D+2
!!
!!---  Enthalpy of mixing  ---
!!
!      HMX = 0.D+0
!      NC = 0
!      DO 20 I = 0,4
!        DO 10 J = 0,I
!          NC = NC + 1
!          HMX = HMX + CBX(NC)*(TX**(I-J))*(YLSX**J)
!   10   CONTINUE
!   20 CONTINUE
!!
!!---  Enthalpy of pure sodium chloride (halite)  ---
!!
!      HSX = 0.D+0
!      DO 30 I = 1,5
!        REALX = REAL(I)
!        HSX = HSX + CAX(I)*(TKX**I)/REALX
!   30 CONTINUE
!      HSX = HSX*1.D+3/WTMS + CAX(6)
!!
!!---  Enthalpy of pure water at vapor-saturated conditions  ---
!!
!      ISRX = 1
!      CALL ENTH_W( TX,PSWX,HLWX,HVWX,ISRX )
!!
!!---  Enthalpy of brine  ---
!!
!      HBX = (1.D+0-XLSX)*HLWX + XLSX*HSX + HMX
!
!---  Enthalpy of mixing  ---
!
      DHMX = 0.D+0
      NC = 0
      DO 20 I = 0,3
        DO 10 J = 0,2
          NC = NC + 1
          DHMX = DHMX + SAX(NC)*(TX**I)*(GLSX**J)
   10   CONTINUE
   20 CONTINUE
      DHMX = (4.184D+3/(1.D+3+WTMS*GLSX))*DHMX
!
!---  Enthalpy of pure sodium chloride (halite)  ---
!
      HSX = 4.184D+3*(SBX(1)*(TX**3) + SBX(2)*(TX**2) +
     &  SBX(3)*TX)/WTMS
!
!---  Enthalpy of pure water at vapor-saturated conditions  ---
!
!      HLWX = 1.D+3*(SCX(1)*(TX**3) + SCX(2)*(TX**2) + SCX(3)*TX +
!     &   SCX(4))
!
!---  Enthalpy of brine  ---
!
      HBX = (1.D+0-XLSX)*HLWX + XLSX*HSX + GLSX*DHMX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTH_B group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ENTH_H( TX,XHAX,XHOX,XHNX,XHWX,HGAX,HGOX,HGNX,
     &  HHX,HLWX )
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
!     Enthalpy of CH4-CO2 hydrate based on liquid water, gaseous CH4,
!     and gaseous CO2 at 273.15 K.
!
!     Yoon, J, Y. Yamamoto, T. Komai, and H. Haneda.  2003.  "Rigorous
!     approach to the predicition of the heat of dissociation of gas
!     hydrates."  Ind. Eng. Chem. Res.  42:1111-1114.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 November 2004.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NAPL
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
      SUB_LOG(ISUB_LOG) = '/ENTH_H'
      TKX = TX+TABS
!
!---  Hydrate mole fractions  ---
!
      WTMX = 1.D+0/(XHAX/WTMA + XHOX/WTMO + XHNX/WTMN + XHWX/WTMW)
      XMHAX = (XHAX/WTMA)*WTMX
      XMHOX = (XHOX/WTMO)*WTMX
      XMHNX = (XHNX/WTMN)*WTMX
      XMHWX = (XHWX/WTMW)*WTMX
!
!---  Hydrate mole fractions of formers  ---
!
      YMHAX = XMHAX/(XMHAX+XMHOX+XMHNX+SMALL)
      YMHOX = XMHOX/(XMHAX+XMHOX+XMHNX+SMALL)
      YMHNX = XMHNX/(XMHAX+XMHOX+XMHNX+SMALL)
!
!---  Hydrate molecular weight, where 1 mole of hydrate
!     is M-guest + n*M-H2O, where M-guest is the guest molecule 
!     molecular weight, n is the hydration number, and M-H20 is the
!     water molecular weight.
!
      IF( (1.D+0-XMHWX).GT.EPSL ) THEN
        WTMHX = YMHOX*WTMO + YMHAX*WTMA + YMHNX*WTMN + 
     &    (XMHWX/(1.D+0-XMHWX))*WTMW
      ELSE
        WTMHX = WTMW
      ENDIF
!
!---  Enthalpy contribution of liquid water, J/kg hydrate  ---
!
      HHX = XHWX*HLWX
!
!---  Enthalpy contribution of gaseous CO2, J/kg hydrate  ---
!
      HHX = HHX + XHAX*HGAX
!
!---  Enthalpy contribution of gaseous CH4, J/kg hydrate  ---
!
      HHX = HHX + XHOX*HGOX
!
!---  Enthalpy contribution of gaseous N2, J/kg hydrate  ---
!
      HHX = HHX + XHNX*HGNX
!
!---  Heat of dissociation for pure CO2 hydrate to liquid water
!     and gaseous CO2 as a function of temperature
!     from Kamath, J/mol hydrate, corrected to Yoon value
!     at 273.15 of 57.66 KJ/mol, n = 6.21 ---
!
      HDAX = 4.184D+0*(17.86466D+3 - 14.95D+0*TKX)
!
!---  Heat of dissociation for pure CH4 hydrate to liquid water
!     and gaseous CH4 as a function of temperature
!     from Kamath J/mol hydrate, corrected to Yoon value
!     at 273.15 of 53.81 KJ/mol, n = 6.07 ---
!
      HDOX = 4.184D+0*(13.95896D+3 - 4.02D+0*TKX)
!
!---  Heat of dissociation (J/mol) for pure N2 hydrate to liquid water
!     and gaseous N2 at 273.15 of 49.54 KJ/mol, n = 6.28  ---
!
      HDNX = 49.54D+3
!
!---  Heat of dissociation from hydrate to liquid water and gaseous
!     CO2, CH4, and N2 J/kg hydrate  ---
!
      HHX = HHX + (-1.D+3)*(YMHAX*HDAX + YMHOX*HDOX + YMHNX*HDNX)/WTMHX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTH_H group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ENTH_I( TX,HIX )
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
!     Calculate the enthalpy of ice as a function of temperature using a
!     using polynomial fit of data from ASHRAE (1977) with a 
!     reference point of liquid water at 273.15 K.
!
!     ASHRAE, 1977; Psychrometric Tables given by ASHRAE Handbook 1977
!     Fundamentals, American Society of Heating, Refrigerating and
!     Air-Conditioning Engineering, Inc., 345 East 47th Street, New
!     York, N.Y.
!
!     The temperature is limited in this subroutine to the following
!     values:  -100.0 C < T <= 0.0 C
!
!----------------------Authors-----------------------------------------!
!
!     Written by WE Nichols, PNNL, September 1994.
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
      REAL*8 A(3)
!
!----------------------Data Statements---------------------------------!
!
      DATA A / -6.36443D+5,1.14468D+2,3.64104D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ENTH_I'
      TKX = MIN( TX+TABS,TABS )
!
!---  Polynomial fit for reference point of liquid water
!     at 273.15 K (includes the heat of dissociation 333.7 kJ/kg)  ---
!
      HIX = A(1) + A(2)*TKX + A(3)*(TKX**2)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTH_I group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ENTH_ID( TX,HIDX )
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
!     Molar ideal gas enthalpy of a mixture, using a reference
!     temperature of 273.15 K.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 January 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYDT
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
      REAL*8 HIDX(LHF_HT)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ENTH_ID'
      TKX = TX+TABS
!
!---  Loop over components  ---
!
      DO 10 NF = 1,NHF_HT
        HIDX(NF) = (HCPP(18,NF)*(TKX-TABS) +
     &    HCPP(19,NF)*((TKX**2)-(TABS**2))/2.D+0 +
     &    HCPP(20,NF)*((TKX**3)-(TABS**3))/3.D+0 +
     &    HCPP(21,NF)*((TKX**4)-(TABS**4))/4.D+0)
   10 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTH_ID group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ENTH_L( TX,XLSX,XLAX,XLOX,XLNX,HBX,HGAX,HGOX,HGNX,HLX )
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
!     This subroutine calculates the enthalpy of H2O-NaCl-CO2 solutions
!     as a function of brine enthalpy, gaseous C02 enthalpy,
!     temperature, mass fraction of dissolved NaCl, mass fraction
!     of dissolved CO2 and mass fraction of dissolved CH4.
!
!     Battistelli, A., C. Claudio, and K. Pruess.  1997.  "The simulator
!     TOUGH2/EWASG for modelling geothermal reservoirs with brines and
!     noncondensible gas."  Geothermics, 26(4): 437-464.
!
!     Himmelblau, D. M.  1959.  "Partial molal heats and entropies of
!     solution for gases dissolved in water from the freezing point
!     to near the critical point."  Journal of Physical Chemistry,
!     63:1803-1808.
!
!     TX - temperature, C
!     XLSX - dissolved NaCl mass fraction
!     XLAX - dissolved CO2 mass fraction
!     XLOX - dissolved CH4 mass fraction
!     XLOX - dissolved N2 mass fraction
!     HBX - enthalpy of brine (NaCl-H2O), J/kg
!     HGAX - enthalpy of gaseous and nonaqueous CO2, J/kg
!     HGOX - enthalpy of gaseous and nonaqueous CH4, J/kg
!     HGNX - enthalpy of gaseous and nonaqueous N2, J/kg
!     HLX - enthalpy of aqueous phase, J/kg
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 16 April 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NAPL
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
      SUB_LOG(ISUB_LOG) = '/ENTH_L'
!
!---  Partial differential of Henry's constant for CO2 at constant 
!     pressure with respect to temperature  ---
!
      DTX = 1.D-6
      CALL HC_LA( TX,XLSX,HCX,DHCX )
      DHCX = DHCX/HCX
!
!---  Heat of solution for CO2 from Himmelblau (1959), J/kg  ---
!
      TKX = TX + TABS
      HSAX = -RCU*(TKX**2)*DHCX/WTMA
!
!---  Partial differential of Henry's constant for CH4 at constant 
!     pressure with respect to temperature  ---
!
      DTX = 1.D-6
      CALL HC_LO( TX,XLSX,HCX,DHCX )
      DHCX = DHCX/HCX
!
!---  Heat of solution for CH4 from Himmelblau (1959), J/kg  ---
!
      TKX = TX + TABS
      HSOX = -RCU*(TKX**2)*DHCX/WTMO
!
!---  Partial differential of Henry's constant for N2 at constant 
!     pressure with respect to temperature  ---
!
      DTX = 1.D-6
      CALL HC_LN( TX,XLSX,HCX,DHCX )
      DHCX = DHCX/HCX
!
!---  Heat of solution for N2 from Himmelblau (1959), J/kg  ---
!
      TKX = TX + TABS
      HSNX = -RCU*(TKX**2)*DHCX/WTMN
!
!---  Aqueous enthalpy  ---
!
!      HLX = (1.D+0-XLAX-XLOX-XLNX)*HBX
!     &  + XLAX*(HGAX + HSAX) + XLOX*(HGOX + HSOX)  + XLNX*(HGNX + HSNX)
      HLX = (1.D+0-XLAX-XLOX-XLNX)*HBX
     &  + XLAX*HGAX + XLOX*HGOX  + XLNX*HGNX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTH_L group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ENTH_N2( TX,PX,HGAX,HNAX,UGAX,UNAX )
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
!     Enthalpy and internal energy of CO2.
!
!     Span, R., E. W. Lemmon, R. T. Jacobsen, W. Wagner, and 
!     A. Yokozeki. 2000.  A reference equation of state for the 
!     thermodynamic properties of nitrogen for temperatures from 
!     63.151 to 1000 K and pressures to 2200 MPa.  
!     J. Phys. Chem. Ref. Data, 29(6):1361-1433.
! 
!     I_VX = 1 subcritical gas
!     I_VX = 2 subcritical liquid
!     I_VX = 3 supercritical gas
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 April 2002.
!     Last Modified by MD White, PNNL, 25 April 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
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
      INTEGER I_PX(2),I_TX(2)
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ENTH_N2'
!
!---  Table lookup and bilinear interpolation  ---
!
      NF = 3
      PTX = MAX( PX,P_TA(1,NF)*1.D+6 )
      CALL ITL_A( PTX,TX,I_PX,I_TX,I_VX,NF)
      CALL PTL_A( PTX,TX,H_TA,HAX,I_PX,I_TX,NF )
      CALL PTL_A( PTX,TX,U_TA,UAX,I_PX,I_TX,NF )
!
!---  Subcritical gas conditions, set liquid enthalpy and internal
!     energy to saturated liquid conditions at the same temperature ---
!
      IF( I_VX.EQ.1 ) THEN
        HGAX = HAX
        UGAX = UAX
!
!---    Find saturation indices  ---
!
        TKX = TX + TABS
        IF( TKX.LT.T_LV(1,NF) ) TKX = T_LV(1,NF)
        IF( TKX.GT.T_LV(I_LV(NF),NF) ) TKX = T_LV(I_LV(NF),NF)
        ITLX = 1
        ITUX = I_LV(NF)
   10   IF( ITUX-ITLX.GT.1 ) THEN
          ITM = (ITLX+ITUX)/2
          IF( (T_LV(I_LV(NF),NF).GT.T_LV(1,NF)).EQV.
     &      (TKX.GT.T_LV(ITM,NF)) ) THEN
            ITLX = ITM
          ELSE
            ITUX = ITM
          ENDIF
          GOTO 10
        ENDIF
        TKLX = T_LV(ITLX,NF)
        TKUX = T_LV(ITLX+1,NF)
        HLX = HL_LV(ITLX,NF)
        HUX = HL_LV(ITLX+1,NF)
        HNAX = ((TKX-TKLX)/(TKUX-TKLX))*(HUX-HLX) + HLX
        ULX = UL_LV(ITLX,NF)
        UUX = UL_LV(ITLX+1,NF)
        UNAX = ((TKX-TKLX)/(TKUX-TKLX))*(UUX-ULX) + ULX
!
!---  Subcritical liquid conditions, set gas enthalpy and internal
!     energy to saturated gas conditions at the same temperature ---
!
      ELSEIF( I_VX.EQ.2 ) THEN
        HNAX = HAX
        UNAX = UAX
!
!---    Find saturation indices  ---
!
        TKX = TX + TABS
        IF( TKX.LT.T_LV(1,NF) ) TKX = T_LV(1,NF)
        IF( TKX.GT.T_LV(I_LV(NF),NF) ) TKX = T_LV(I_LV(NF),NF)
        ITLX = 1
        ITUX = I_LV(NF)
   20   IF( ITUX-ITLX.GT.1 ) THEN
          ITM = (ITLX+ITUX)/2
          IF( (T_LV(I_LV(NF),NF).GT.T_LV(1,NF)).EQV.
     &      (TKX.GT.T_LV(ITM,NF)) ) THEN
            ITLX = ITM
          ELSE
            ITUX = ITM
          ENDIF
          GOTO 20
        ENDIF
        TKLX = T_LV(ITLX,NF)
        TKUX = T_LV(ITLX+1,NF)
        HLX = HV_LV(ITLX,NF)
        HUX = HV_LV(ITLX+1,NF)
        HGAX = ((TKX-TKLX)/(TKUX-TKLX))*(HUX-HLX) + HLX
        ULX = UV_LV(ITLX,NF)
        UUX = UV_LV(ITLX+1,NF)
        UGAX = ((TKX-TKLX)/(TKUX-TKLX))*(UUX-ULX) + ULX
!
!---  Supercritical conditions, set gas and liquid enthalpy and
!     and internal energy to the gas value  ---
!
      ELSEIF( I_VX.EQ.3 ) THEN
        HGAX = HAX
        UGAX = UAX
        HNAX = HAX
        UNAX = UAX
!
!---  Subcritical gas conditions below the triple-point  ---
!
      ELSE
        INDX = 4
        CHMSGX = 'Solid Phase for N2: T = '
        IMSGX = 0
        NMSGX = 0
        SUBLOGX = 'ENTH_N2'
        RLMSG = TX
        CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
      ENDIF    
!
!---  Convert enthalpy and internal energy units to J/kg  ---
!
      HGAX = 1.D+3*HGAX
      HNAX = 1.D+3*HNAX
      UGAX = 1.D+3*UGAX
      UNAX = 1.D+3*UNAX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTH_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ENTH_RES( PSOX,TX,HRESX,ZMCX,INDX )
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
!     Molar residual gas enthalpy of a mixture, using a reference
!     temperature of 273.15 K.
!
!     INDX = 1 Gas
!     INDX = 2 Liquid
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 January 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYDT
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
      REAL*8 ZMCX(LHF_HT),HRESX(LHF_HT),AIX(LHF_HT)
      REAL*8 DAIXDT(LHF_HT),DAXDT(LHF_HT)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ENTH_RES'
      TKX = TX+TABS
!
!---  Peng-Robinson temperature dependent component 
!     coefficients  ---
!
      CALL PRC_T( AIX,TKX )
!
!---  Peng-Robinson nonaqueous mixture coefficients  ---
!
      CALL PRC_X( AX,AIX,BX,ZMCX )
!
!---  Peng-Robinson mixture compressibility ---
!
      CALL PRC_Z( AX,BX,CAX,CBX,PSOX,TKX,ZX,INDX )
      V7X = (ZX + (SQRT(2.D+0) + 1.D+0)*CBX)
      V8X = (ZX - (SQRT(2.D+0) - 1.D+0)*CBX)
!
!---  Partial derivative of component a with respect to temperature  ---
!
      DO 100 NF = 1,NHF_HT
        TRX = TKX/HCPP(3,NF)
        DAIXDT(NF) = -HCPP(13,NF)*(1.D+0 + HCPP(15,NF)*(1.D+0 -
     &    SQRT(TRX)))*(HCPP(15,NF)/(SQRT(TRX)*HCPP(3,NF)))
  100 CONTINUE
!
!---  Partial derivative of mixture a with respect to temperature  ---
!
      DO 120 I = 1,NHF_HT
        DAXDT(I) = 0.D+0
        DO 110 J = 1,NHF_HT
          DAXDT(I) = DAXDT(I) + (5.D-1*ZMCX(J)/SQRT(AIX(I)*AIX(J)))*
     &      (DAIXDT(I)*AIX(J) + DAIXDT(J)*AIX(I))
  110   CONTINUE
  120 CONTINUE
!
!---  Component residual enthalpy  ---
!
      DO 130 NF = 1,NHF_HT
        HRESX(NF) = (ZX - 1.D+0)*RCU*TKX + 
     &    ((TKX*DAXDT(NF) - AX)/(SQRT(8.D+0)*BX))*LOG(V7X/V8X)
  130 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTH_RES group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ENTH_S( TX,HSX )
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
!     Enthalpy of precipitated NaCl.
!
!     Lide, D.R. and H.V. Kehiaian.  1994.  CRC Handbook of
!     Thermophysical and Thermochemical Data, CRC Press, Inc.,
!     Boca Raton, Florida, pp. 97-98.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 1 May 2002.
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
      REAL*8 CAX(5)
!
!----------------------Data Statements---------------------------------!
!
      DATA CAX / 25.19D+0, 0.1973D+0, -6.0114D-4, 8.81505D-7,
     &  -4.765D-10 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ENTH_S'
!
!---  0 C Reference  ---
!
      TKX = TX+TABS
      HSX = -1.24858D-4
      DO 10 I = 1,5
        REALX = REAL(I)
        HSX = HSX + CAX(I)*(TKX**I)/REALX
   10 CONTINUE
!
!---  Convert from J/mol to J/kg  ---
!
      HSX = 1.D+3*HSX/WTMS
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTH_S group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ENTH_W( TX,PX,HLX,HVX,ISRX )
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
!     Enthalpy and internal energy (J/kg) of pure water as a function
!     of temperature and pressure.
!
!     Pressure Range: 0 - 100 MPa (1000 Bar)
!     Temperature Range:  273.16 K (0.01 C) to 1073.15K (800 C)
!
!     Meyer, C.A., R.B. McClintock, G.J. Silvestri, and R.C. Spencer
!     1993.  ASME Steam Tables, The American Society of Mechanical
!     Engineers, New York.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 March 2002.
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
      REAL*8 CAX(23),SAX(12),CBX(31),SBX(5),LFCX(3),SLX
      REAL*8 TX,PX
      INTEGER ISRX,INX(8),IZX(8,3),ITX(8),IXX(8,2)
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Data Statements---------------------------------!
!
      DATA CAX / 6.824687741D+3, -5.422063673D+2, -2.096666205D+4,
     &  3.941286787D+4, -6.733277739D+4, 9.902381028D+4,
     &  -1.093911774D+5, 8.590841667D+4, -4.511168742D+4,
     &  1.418138926D+4, -2.017271113D+3, 7.982692717D+0,
     &  -2.616571843D-2, 1.522411790D-3, 2.284279054D-2,
     &  2.421647003D+2, 1.269716088D-10, 2.074838328D-7,
     &  2.174020350D-8, 1.105710498D-9, 1.293441934D+1,
     &  1.308119072D-5, 6.047626338D-14 /
      DATA SAX / 8.438375405D-1, 5.362162162D-4, 1.720000000D+0,
     &  7.342278489D-2, 4.975858870D-2, 6.537154300D-1,
     &  1.150000000D-6, 1.150800000D-5, 1.418800000D-1,
     &  7.002753165D+0, 2.995284926D-4, 2.040000000D-1 /
      DATA CBX / 1.683599274D+1, 2.856067796D+1, -5.438923329D+1,
     &  4.330662834D-1, -6.547711697D-1, 8.565182058D-2,
     &  6.670375918D-2, 1.388983801D+0, 8.390104328D-2,
     &  2.614670893D-2, -3.373439453D-2, 4.520918904D-1,
     &  1.069036614D-1, -5.975336707D-1, -8.847535804D-2,
     &  5.958051609D-1, -5.159303373D-1, 2.075021122D-1,
     &  1.190610271D-1, -9.867174132D-2, 1.683998803D-1,
     &  -5.809438001D-2, 6.552390126D-3, 5.710218649D-4,
     &  1.936587558D+2, -1.388522425D+3, 4.126607219D+3,
     &  -6.508211677D+3, 5.745984054D+3, -2.693088365D+3,
     &  5.235718623D+2/
      DATA SBX / 7.633333333D-1, 4.006073948D-1, 8.636081627D-2,
     &  -8.532322921D-1, 3.460208861D-1 /
      DATA LFCX / 1.574373327D+1, -3.417061978D+1, 1.931380707D+1 /
      DATA SLX / 4.260321148D+0 /
      DATA INX / 2, 3, 2, 2, 3, 2, 2, 2 /
      DATA IZX / 13, 18, 18, 25, 32, 12, 24, 24,
     &  3, 2, 10, 14, 28, 11, 18, 14,
     &  0, 1, 0, 0, 24, 0, 0, 0 /
      DATA ITX / 0, 0, 0, 0, 0, 1, 1, 2 /
      DATA IXX / 0, 0, 0, 0, 0, 14, 19, 54,
     &  0, 0, 0, 0, 0, 0, 0, 27 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ENTH_W'
!
!---  Restrict temperature to above 0.01 C and below 800 C  ---
!
      TY = MIN( MAX( TX,1.D-2 ),8.D+2 )
!
!---  Restrict pressure to above 0.0 MPa and below 100 MPa  ---
!
      PY = MIN( MAX( PX,0.D+0 ),1.D+8 )
!
!---  Reduced temperature and pressure  ---
!
      THETAX = (TY+TABS)/TCRW
      BETAX = PY/PCRW
!
!---  L-function  ---
!
      BETALX = LFCX(1) + LFCX(2)*THETAX + LFCX(3)*(THETAX**2)
      BETALPX = LFCX(2) + 2.D+0*LFCX(3)*THETAX
!
!---  Subregions 1 or 6  ---
!
      IF( ISRX.EQ.1 .OR. ISRX.EQ.6 ) THEN
        CYX = 1.D+0 - SAX(1)*(THETAX**2) - SAX(2)/(THETAX**6)
        CYPX = -2.D+0*SAX(1)*THETAX + 6.D+0*SAX(2)/(THETAX**7)
        CZX = SAX(3)*(CYX**2)
        CZX = CZX - 2.D+0*SAX(4)*THETAX
        CZX = CZX + 2.D+0*SAX(5)*BETAX
        IF( (CZX/EPSL).GT.EPSL ) THEN
          CZX = CYX + SQRT(CZX)
        ELSE
          CZX = CYX
        ENDIF
        RHX = CAX(1)*THETAX
        DO 10 I = 1,10
          REALX = REAL(I)
          RHX = RHX - (REALX-2.D+0)*CAX(I+1)*(THETAX**(I-1))
   10   CONTINUE
        RHX = RHX + CAX(12)*(CZX*(1.7D+1*((CZX/2.9D+1)-(CYX/1.2D+1)) +
     &    5.D+0*THETAX*(CYPX/1.2D+1)) + SAX(4)*THETAX -
     &    (SAX(3)-1.D+0)*THETAX*CYX*CYPX)/(CZX**(5.D+0/17.D+0))
        RHX = RHX + (CAX(13) - CAX(15)*(THETAX**2) +
     &    CAX(16)*(9.D+0*THETAX + SAX(6))*((SAX(6)-THETAX)**9) +
     &    CAX(17)*(2.D+1*(THETAX**19) +
     &    SAX(7))/((SAX(7) + (THETAX**19))**2))*BETAX
        RHX = RHX - (1.2D+1*(THETAX**11) + SAX(8))*
     &    (CAX(18)*BETAX + CAX(19)*(BETAX**2) + CAX(20)*(BETAX**3))
     &    /((SAX(8) + (THETAX**11))**2)
        RHX = RHX + CAX(20)*(THETAX**18)*
     &    (1.7D+1*SAX(9) + 1.9D+1*(THETAX**2))*
     &    (((SAX(10) + BETAX)**3) + SAX(11)*BETAX)
        RHX = RHX + CAX(22)*SAX(12)*(BETAX**3)
        RHX = RHX + 2.1D+1*CAX(23)*(BETAX**4)/(THETAX**20)
        HLX = 1.D-3*RHX*PCRW*VCRW/WTMW
      ENDIF
!
!---  Below 0.01 C  ---
!
      IF( TX.LT.1.D-2 ) HLX = HLX - 4.202405*(1.D-2-TX)
!
!---  Subregions 2 or 6  ---
!
      IF( ISRX.EQ.2 .OR. ISRX.EQ.6 ) THEN
        CXX = EXP(SBX(1)*(1.D+0-THETAX))
        RHX = CBX(1)*THETAX
        ICBX = 1
        DO 20 I = 1,5
          ICBX = ICBX + 1
          REALX = REAL(I)
          RHX = RHX - CBX(ICBX)*(REALX-2.D+0)*(THETAX**(I-1))
   20   CONTINUE
        ICBX = 6
        DO 40 I = 1,5
          RHAX = 0.D+0
          DO 30 J = 1,INX(I)
            ICBX = ICBX + 1
            RHAX = RHAX + CBX(ICBX)*(1.D+0 + IZX(I,J)*SBX(1)*THETAX)*
     &        (CXX**IZX(I,J))
   30     CONTINUE
          RHX = RHX - (BETAX**I)*RHAX
   40   CONTINUE
        ICBX = 18
        DO 70 I = 6,8
          RHBX = 0.D+0
          RHCX = 0.D+0
          ISBX = 1
          DO 50 K = 1,ITX(I)
            ISBX = ISBX + 1
            RHBX = RHBX + IXX(I,K)*SBX(ISBX)*(CXX**IXX(I,K))
            RHCX = RHCX + SBX(ISBX)*(CXX**IXX(I,K))
   50     CONTINUE
          RHBX = RHBX*SBX(1)*THETAX
          RHCX = RHCX + (BETAX**(2-I))
          DO 60 J = 1,INX(I)
            ICBX = ICBX + 1
            RHAX = CBX(ICBX)*(CXX**IZX(I,J))
            RHAX = RHAX*((1.D+0 + IZX(I,J)*SBX(1)*THETAX) - RHBX/RHCX )
   60     CONTINUE
          RHX = RHX - SBX(1)*RHAX/RHCX
   70   CONTINUE
        RHAX = BETAX*((BETAX/BETALX)**10)
        RHBX = 1.D+1*BETALPX/BETALX
        ICBX = 24
        DO 80 I = 0,6
          ICBX = ICBX + 1
          REALX = REAL(I)
          RHX = RHX + RHAX*(1.D+0 + THETAX*(RHBX + REALX*SBX(1))*
     &      CBX(ICBX)*(CXX**I))
   80   CONTINUE
        HVX = 1.D-3*RHX*PCRW*VCRW/WTMW
      ENDIF
!
!---  Subregions 3 or 5  ---
!
      IF( ISRX.EQ.3 .OR. ISRX.EQ.5 ) THEN
        INDX = 4
        CHMSGX = 'Steam Table Subregion 3/5: Exceeds Temperature ' // 
     &    'Range: '
        IMSGX = 0
        NMSGX = 0
        SUBLOGX = 'ENTH_W'
        RLMSGX = TX
        CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
      ENDIF
!
!---  Subregions 4 or 5  ---
!
      IF( ISRX.EQ.4 .OR. ISRX.EQ.5 ) THEN
        INDX = 34
        CHMSGX = 'Steam Table Subregion 4/5: Exceeds Temperature ' // 
     &    'Range: '
        IMSGX = 0
        NMSGX = 0
        SUBLOGX = 'ENTH_W'
        RLMSGX = TX
        CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTH_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FUG( AX,AIX,BX,CAX,CBX,FUGX,FUGLX,XMCX,ZX )
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
!     Fugacity coefficient.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 March 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYDT
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
      REAL*8 AIX(LHF_HT),FUGX(LHF_HT),FUGLX(LHF_HT),XMCX(LHF_HT)
      REAL*8 FPRKIJ(3,3)
!
!----------------------Data Statements---------------------------------!
!
      DATA FPRKIJ / 0.D+0, 0.109267D+0, 0.0195133D+0, 0.109267D+0, 
     &  0.D+0, 0.0262942D+0, 0.0195133D+0, 0.0262942D+0, 0.D+0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FUG'
      IF( CAX.GT.EPSL ) THEN
        DO 100 IHF = 1,NHF_HT
          V1X = -LOG(ZX-CBX)
          V2X = (ZX-1.D+0)*HCPP(14,IHF)/BX
          V3X = CAX/((2.D+0**1.5D+0)*CBX)
          V4X = 2.D+0*SQRT(AIX(IHF))/AX
          V5X = 0.D+0
          DO 10 JHF = 1,NHF_HT
            V5X = V5X + XMCX(JHF)*SQRT(AIX(JHF))*(1.D+0-FPRKIJ(IHF,JHF))
   10     CONTINUE
          V6X = HCPP(14,IHF)/BX
          V7X = (ZX + (SQRT(2.D+0) + 1.D+0)*CBX)
          V8X = (ZX - (SQRT(2.D+0) - 1.D+0)*CBX)
          V3X = V3X*(V4X*V5X-V6X)*LOG( V7X/V8X )
          FUGLX(IHF) = V1X + V2X - V3X
          FUGX(IHF) = EXP( FUGLX(IHF) )
  100   CONTINUE
      ELSE
        DO 110 IHF = 1,NHF_HT
          FUGLX(IHF) = 0.D+0
          FUGX(IHF) = 1.D+0
  110   CONTINUE
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FUG group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE HC_LA( TX,XLSX,HCX,DHCX )
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
!     Henry's coefficient as a function of temperature for the
!     dissolution of CO2 in NaCl aqueous solutions.
!
!     Battistelli, A., C. Claudio, and K. Pruess.  1997.  "The simulator
!     TOUGH2/EWASG for modelling geothermal reservoirs with brines and
!     noncondensible gas."  Geothermics, 26(4): 437-464.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 March 2002.
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
      REAL*8 COBX(6),COCX(5)
      REAL*8 CHBX(4),CHCX(4)
!
!----------------------Data Statements---------------------------------!
!
      DATA COBX / 7.83666D+7, 1.96025D+6, 8.20574D+4, -7.40674D+2,
     &  2.18380D+0, -2.20999D-3 /
      DATA COCX / 1.19784D-1, -7.17823D-4, 4.93854D-6, -1.03826D-8,
     &  1.08233D-11 /
      DATA CHBX / -3.65183D+2, 1.81067D+4, 4.97554D+1, -2.85D-4 /
      DATA CHCX / 1.668D-1, -1.09236D-3, 7.81043D-6, -1.52343D-8 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/HC_LA'
!
!---  Empirical formulation by Battistelli for Henry's coefficient
!     for pure water as a function of temperature  ---
!
      HCX = 0.D+0
      DO 10 I = 0,5
        HCX = HCX + COBX(I+1)*(TX**I)
   10 CONTINUE
      DHCX = 0.D+0
      DO 20 I = 1,5
        DHCX = DHCX + REAL(I)*COBX(I+1)*(TX**(I-1))
   20 CONTINUE
!
!---  Empirical formulation by Battistelli for salting-out
!     coefficient for NaCl aqueous solutions as a function
!     of temperature  ---
!
      SKBX = 0.D+0
      DO 30 I = 0,4
        SKBX = SKBX + COCX(I+1)*(TX**I)
   30 CONTINUE
      DSKBX = 0.D+0
      DO 40 I = 1,4
        DSKBX = DSKBX + REAL(I)*COCX(I+1)*(TX**(I-1))
   40 CONTINUE
!
!---  Empirical formulation by Battistelli for Henry's coefficient
!     coefficient for NaCl aqueous solutions as a function
!     of temperature and salt molality  ---
!
      GLSX = 1.D+3*XLSX/(WTMS*(1.D+0-XLSX))
      HCX = HCX*(1.D+1**(GLSX*SKBX))
      DHCX = DHCX*(1.D+1**(GLSX*SKBX)) + HCX*GLSX*DSKBX*LOG(1.D+1)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of HC_LA group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE HC_LN( TX,XLSX,HCX,DHCX )
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
!     Henry's coefficient as a function of temperature for the
!     dissolution of N2 in NaCl aqueous solutions.
!
!     Sloan, E. D., Jr.  1997.  Clathrate Hydrates of Natural Gases.
!     Marcel Dekker, Inc., New York.
!
!     Battistelli, A., C. Claudio, and K. Pruess.  1997.  "The simulator
!     TOUGH2/EWASG for modelling geothermal reservoirs with brines and
!     noncondensible gas."  Geothermics, 26(4): 437-464.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 8 October 2004.
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
      REAL*8 CHBX(4),CHCX(4)
!
!----------------------Data Statements---------------------------------!
!
      DATA CHBX / -3.27850D+2, 1.67576D+4, 4.28400D+1, 1.6765D-2 /
      DATA CHCX / 1.668D-1, -1.09236D-3, 7.81043D-6, -1.52343D-8 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/HC_LN'
!
!---  The formulation by Sloan for Henry's coefficient of nitrogen 
!     for pure water as a function of temperature  ---
!
      TKX = TX + TABS
      RCX = 1.9872
      HCX = (CHBX(1)/RCX) + (CHBX(2)/(RCX*TKX)) + 
     &  ((CHBX(3)/RCX)*LOG(TKX)) + (CHBX(4)*TKX/RCX)
      DHCX = (CHBX(2)/(RCX*(TKX**2))) - 
     &  (CHBX(3)/(RCX*TKX)) - (CHBX(4)/RCX)
      HCX = EXP(-HCX)
      HCX = HCX*PATM
      DHCX = DHCX*HCX
!
!---  Empirical formulation for salting-out coefficient for 
!     NaCl aqueous solutions as a function of temperature  
!
      SKBX = 0.D+0
      DO 10 I = 0,3
        SKBX = SKBX + CHCX(I+1)*(TX**I)
   10 CONTINUE
      DSKBX = 0.D+0
      DO 20 I = 1,3
        DSKBX = DSKBX + REAL(I)*CHCX(I+1)*(TX**(I-1))
   20 CONTINUE
!
!---  Empirical formulation by Battistelli for Henry's coefficient
!     coefficient for NaCl aqueous solutions as a function
!     of temperature and salt molality  ---
!
      GLSX = 1.D+3*XLSX/(WTMS*(1.D+0-XLSX))
      HCX = HCX*(1.D+1**(GLSX*SKBX))
      DHCX = DHCX*(1.D+1**(GLSX*SKBX)) + HCX*GLSX*DSKBX*LOG(1.D+1)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of HC_LN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE HC_LO( TX,XLSX,HCX,DHCX )
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
!     Henry's coefficient as a function of temperature for the
!     dissolution of CH4 in NaCl aqueous solutions.
!
!     Sloan, E. D., Jr.  1997.  Clathrate Hydrates of Natural Gases.
!     Marcel Dekker, Inc., New York.
!
!     Battistelli, A., C. Claudio, and K. Pruess.  1997.  "The simulator
!     TOUGH2/EWASG for modelling geothermal reservoirs with brines and
!     noncondensible gas."  Geothermics, 26(4): 437-464.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 8 October 2004.
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
      REAL*8 CHBX(4),CHCX(4)
!
!----------------------Data Statements---------------------------------!
!
      DATA CHBX / -3.65183D+2, 1.81067D+4, 4.97554D+1, -2.85D-4 /
      DATA CHCX / 1.668D-1, -1.09236D-3, 7.81043D-6, -1.52343D-8 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/HC_LO'
!
!---  The formulation by Sloan for Henry's coefficient of methane 
!     for pure water as a function of temperature  ---
!
      TKX = TX + TABS
      RCX = 1.9872
      HCX = (CHBX(1)/RCX) + (CHBX(2)/(RCX*TKX)) + 
     &  ((CHBX(3)/RCX)*LOG(TKX)) + (CHBX(4)*TKX/RCX)
      DHCX = (CHBX(2)/(RCX*(TKX**2))) - 
     &  (CHBX(3)/(RCX*TKX)) - (CHBX(4)/RCX)
      HCX = EXP(-HCX)
      HCX = HCX*PATM
      DHCX = DHCX*HCX
!
!---  Empirical formulation for salting-out coefficient for 
!     NaCl aqueous solutions as a function of temperature  
!
      SKBX = 0.D+0
      DO 10 I = 0,3
        SKBX = SKBX + CHCX(I+1)*(TX**I)
   10 CONTINUE
      DSKBX = 0.D+0
      DO 20 I = 1,3
        DSKBX = DSKBX + REAL(I)*CHCX(I+1)*(TX**(I-1))
   20 CONTINUE
!
!---  Empirical formulation by Battistelli for Henry's coefficient
!     coefficient for NaCl aqueous solutions as a function
!     of temperature and salt molality  ---
!
      GLSX = 1.D+3*XLSX/(WTMS*(1.D+0-XLSX))
      HCX = HCX*(1.D+1**(GLSX*SKBX))
      DHCX = DHCX*(1.D+1**(GLSX*SKBX)) + HCX*GLSX*DSKBX*LOG(1.D+1)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of HC_LO group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE HYD_X( TKX,ZMCX,XLSX,RHOHX,XHAX,
     &  XHOX,XHNX,XHWX,XMHAX,XMHOX,XMHNX,XMHWX )
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
!     Compute hydrate composition from temperature and hydrate gas
!     composition.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 August 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
      USE NAPL
      USE HYDT
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
      REAL*8 FOCX(2,LHF_HT)
      REAL*8 TCPWX(2)
      REAL*8 ZMCX(LHF_HT)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/HYD_X'
!
!---  Inhibitor correction  ---
!
      DTHX = 0.D+0
      IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
        CALL BCIH_TI( THX,PHE_HT,PHE2T_HT,PHX,ZMCX )
        XLSFX = 1.D+2*XLSX
        POX = 1.D-3*PATM
        DTHX = (CINH(1,I_INH)*XLSFX + CINH(2,I_INH)*(XLSFX**2)
     &    + CINH(3,I_INH)*(XLSFX**3))*
     &    (CINH(4,I_INH)*LOG(PHX*1.D-3) + CINH(5,I_INH))*
     &    (CINH(6,I_INH)*(POX-1.D+3)+1.D+0)
      ENDIF
      THX = TKX + DTHX
!
!---  Small-cage occupancies using bilinear interpolation on
!     concentration of hydrate formers  ---
!
      IF( ZMCX(1)/EPSL.LT.EPSL ) THEN
        XSCAX = 0.D+0
      ELSE
        CALL BCIH_TI( THX,XSCA_HT,XSCA2T_HT,XSCAX,ZMCX )
      ENDIF
      IF( ZMCX(2)/EPSL.LT.EPSL ) THEN
        XSCOX = 0.D+0
      ELSE
        CALL BCIH_TI( THX,XSCO_HT,XSCO2T_HT,XSCOX,ZMCX )
      ENDIF
      IF( ZMCX(3)/EPSL.LT.EPSL ) THEN
        XSCNX = 0.D+0
      ELSE
        CALL BCIH_TI( THX,XSCN_HT,XSCN2T_HT,XSCNX,ZMCX )
      ENDIF
!
!---  Large-cage occupancies using bilinear interpolation on
!     concentration of hydrate formers  ---
!
      IF( ZMCX(1)/EPSL.LT.EPSL ) THEN
        XLCAX = 0.D+0
      ELSE
        CALL BCIH_TI( THX,XLCA_HT,XLCA2T_HT,XLCAX,ZMCX )
      ENDIF
      IF( ZMCX(2)/EPSL.LT.EPSL ) THEN
        XLCOX = 0.D+0
      ELSE
        CALL BCIH_TI( THX,XLCO_HT,XLCO2T_HT,XLCOX,ZMCX )
      ENDIF
      IF( ZMCX(3)/EPSL.LT.EPSL ) THEN
        XLCNX = 0.D+0
      ELSE
        CALL BCIH_TI( THX,XLCN_HT,XLCN2T_HT,XLCNX,ZMCX )
      ENDIF
!
!---  Hydrate mole fraction of formers using bilinear interpolation on
!     concentration of hydrate formers  ---
!
       SXSCX = XSCAX+XSCOX+XSCNX
       SXLCX = XLCAX+XLCOX+XLCNX
       YMHAX = (XSCAX + 3.D+0*XLCAX)/(SXSCX + 3.D+0*SXLCX)
       YMHOX = (XSCOX + 3.D+0*XLCOX)/(SXSCX + 3.D+0*SXLCX)
       YMHNX = (XSCNX + 3.D+0*XLCNX)/(SXSCX + 3.D+0*SXLCX)
!
!---  H2O hydrate mole fraction  ---
!
      IF( SXSCX+SXLCX.GT.EPSL ) THEN
        XMHWX = 23.D+0/((3.D+0*(SXLCX))+SXSCX)
        XMHWX = XMHWX/(1.D+0+XMHWX)
        XMHAX = (1.D+0-XMHWX)*YMHAX
        XMHOX = (1.D+0-XMHWX)*YMHOX
        XMHNX = (1.D+0-XMHWX)*YMHNX
      ELSE
        XMHWX = 1.D+0
        XMHAX = 0.D+0
        XMHOX = 0.D+0
        XMHNX = 0.D+0
      ENDIF
!
!---  Hydrate mass fractions  ---
!
      WTMHX = XMHAX*WTMA + XMHOX*WTMO  + XMHNX*WTMN + XMHWX*WTMW
      XHAX = XMHAX*WTMA/WTMHX
      XHOX = XMHOX*WTMO/WTMHX
      XHNX = XMHNX*WTMN/WTMHX
      XHWX = XMHWX*WTMW/WTMHX
!
!---  Hydrate density, assuming sI hydrate crystal cell structure,
!     
!     WMPUC: number of water molecules per unit cell
!     AVGNX: Avagadro's number, molecules/mol
!     HCPP(1,J): molecular weight of component J, g/mol
!     FOCX(I,J): fractional occupancy of cavity I by component J
!     TCPWX(I): number of type I cavities per water molecule 
!       per unit cell
!     VUCX: volume of unit cell, cm^3
!     NCTX: number of cavity types
!     NCHX: number of components in hydrate phase
!     RHOHX: hydrate density, kg/m^3
!
      WMPUC = 46.D+0
      AVGNX = 6.023D+23
      FOCX(1,1) = XSCAX
      FOCX(2,1) = XLCAX
      FOCX(1,2) = XSCOX
      FOCX(2,2) = XLCOX
      FOCX(1,3) = XSCNX
      FOCX(2,3) = XLCNX
      TCPWX(1) = 2.D+0
      TCPWX(2) = 6.D+0
      VUCX = 12.D-8**3
      NCTX = 2
      NCHX = 3
      RHOHX = WMPUC*HCPP(1,4)
      DO 210 J = 1,NCHX
      DO 200 I = 1,NCTX
        RHOHX = RHOHX + FOCX(I,J)*TCPWX(I)*HCPP(1,J)
  200 CONTINUE
  210 CONTINUE
      RHOHX = 1.D+3*RHOHX/(AVGNX*VUCX)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of HYD_X group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE HYD_P( TKX,ZMCX,XLSX,PHX )
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
!     Compute hydrate equilibrium pressure as a function of temperature
!     vapor mole fractions, and inhibitor concentration.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 4 November 2014.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
      USE HYDT
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 ZMCX(LHF_HT)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/HYD_P'
!
!---  Inhibitor correction  ---
!
      DTHX = 0.D+0
      IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
        CALL BCIH_TI( TKX,PHE_HT,PHE2T_HT,PHX,ZMCX )
        XLSFX = 1.D+2*XLSX
        POX = 1.D-3*PATM
        DTHX = (CINH(1,I_INH)*XLSFX + CINH(2,I_INH)*(XLSFX**2)
     &    + CINH(3,I_INH)*(XLSFX**3))*
     &    (CINH(4,I_INH)*LOG(PHX*1.D-3) + CINH(5,I_INH))*
     &    (CINH(6,I_INH)*(POX-1.D+3)+1.D+0)
      ENDIF
      THX = TKX + DTHX
!
!---  Equilibrium pressure using bilinear interpolation on
!     concentration of hydrate formers  ---
!
      CALL BCIH_TI( THX,PHE_HT,PHE2T_HT,PHX,ZMCX )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of HYD_P group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE HYD_C( TKX,ZMCX,XLSX,RHOHX,XHAX,XHOX,
     &  XHNX,XHWX,XMHAX,XMHOX,XMHNX,XMHWX )
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
!     Compute hydrate equilibrium pressure, hydrate composition,
!     and hydrate density
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 10 August 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
      USE NAPL
      USE HYDT
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
      REAL*8 FOCX(2,LHF_HT)
      REAL*8 TCPWX(2)
      REAL*8 ZMCX(LHF_HT)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/HYD_PC'
!
!---  Inhibitor correction  ---
!
      DTHX = 0.D+0
      IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
        CALL BCIH_TI( TKX,PHE_HT,PHE2T_HT,PHX,ZMCX )
        XLSFX = 1.D+2*XLSX
        POX = 1.D-3*PATM
        DTHX = (CINH(1,I_INH)*XLSFX + CINH(2,I_INH)*(XLSFX**2)
     &    + CINH(3,I_INH)*(XLSFX**3))*
     &    (CINH(4,I_INH)*LOG(PHX*1.D-3) + CINH(5,I_INH))*
     &    (CINH(6,I_INH)*(POX-1.D+3)+1.D+0)
      ENDIF
      THX = TKX + DTHX
!
!---  Small-cage occupancies using bilinear interpolation on
!     concentration of hydrate formers  ---
!
      IF( ZMCX(1)/EPSL.LT.EPSL ) THEN
        XSCAX = 0.D+0
      ELSE
        CALL BCIH_TI( THX,XSCA_HT,XSCA2T_HT,XSCAX,ZMCX )
        XSCAX = MIN( MAX( 0.D+0,XSCAX ),1.D+0 )
      ENDIF
      IF( ZMCX(2)/EPSL.LT.EPSL ) THEN
        XSCOX = 0.D+0
      ELSE
        CALL BCIH_TI( THX,XSCO_HT,XSCO2T_HT,XSCOX,ZMCX )
        XSCOX = MIN( MAX( 0.D+0,XSCOX ),1.D+0 )
      ENDIF
      IF( ZMCX(3)/EPSL.LT.EPSL ) THEN
        XSCNX = 0.D+0
      ELSE
        CALL BCIH_TI( THX,XSCN_HT,XSCN2T_HT,XSCNX,ZMCX )
        XSCNX = MIN( MAX( 0.D+0,XSCNX ),1.D+0 )
      ENDIF
!
!---  Large-cage occupancies using bilinear interpolation on
!     concentration of hydrate formers  ---
!
      IF( ZMCX(1)/EPSL.LT.EPSL ) THEN
        XLCAX = 0.D+0
      ELSE
        CALL BCIH_TI( THX,XLCA_HT,XLCA2T_HT,XLCAX,ZMCX )
        XLCAX = MIN( MAX( 0.D+0,XLCAX ),1.D+0 )
      ENDIF
      IF( ZMCX(2)/EPSL.LT.EPSL ) THEN
        XLCOX = 0.D+0
      ELSE
        CALL BCIH_TI( THX,XLCO_HT,XLCO2T_HT,XLCOX,ZMCX )
        XLCOX = MIN( MAX( 0.D+0,XLCOX ),1.D+0 )
      ENDIF
      IF( ZMCX(3)/EPSL.LT.EPSL ) THEN
        XLCNX = 0.D+0
      ELSE
        CALL BCIH_TI( THX,XLCN_HT,XLCN2T_HT,XLCNX,ZMCX )
        XLCNX = MIN( MAX( 0.D+0,XLCNX ),1.D+0 )
      ENDIF
!
!---  Hydrate mole fraction of formers  ---
!
       SXSCX = XSCAX+XSCOX+XSCNX
       SXLCX = XLCAX+XLCOX+XLCNX
       YMHAX = (XSCAX + 3.D+0*XLCAX)/(SXSCX + 3.D+0*SXLCX)
       YMHOX = (XSCOX + 3.D+0*XLCOX)/(SXSCX + 3.D+0*SXLCX)
       YMHNX = (XSCNX + 3.D+0*XLCNX)/(SXSCX + 3.D+0*SXLCX)
!
!---  H2O hydrate mole fraction  ---
!
      IF( SXSCX+SXLCX.GT.EPSL ) THEN
        XMHWX = 23.D+0/((3.D+0*(SXLCX))+SXSCX)
        XMHWX = XMHWX/(1.D+0+XMHWX)
        XMHAX = (1.D+0-XMHWX)*YMHAX
        XMHOX = (1.D+0-XMHWX)*YMHOX
        XMHNX = (1.D+0-XMHWX)*YMHNX
      ELSE
        XMHWX = 1.D+0
        XMHAX = 0.D+0
        XMHOX = 0.D+0
        XMHNX = 0.D+0
      ENDIF
!
!---  Hydrate mass fractions  ---
!
      WTMHX = XMHAX*WTMA + XMHOX*WTMO  + XMHNX*WTMN + XMHWX*WTMW
      XHAX = XMHAX*WTMA/WTMHX
      XHOX = XMHOX*WTMO/WTMHX
      XHNX = XMHNX*WTMN/WTMHX
      XHWX = XMHWX*WTMW/WTMHX
!
!---  Hydrate density, assuming sI hydrate crystal cell structure,
!     
!     WMPUC: number of water molecules per unit cell
!     AVGNX: Avagadro's number, molecules/mol
!     HCPP(1,J): molecular weight of component J, g/mol
!     FOCX(I,J): fractional occupancy of cavity I by component J
!     TCPWX(I): number of type I cavities per water molecule 
!       per unit cell
!     VUCX: volume of unit cell, cm^3
!     NCTX: number of cavity types
!     NCHX: number of components in hydrate phase
!     RHOHX: hydrate density, kg/m^3
!
      WMPUC = 46.D+0
      AVGNX = 6.023D+23
      FOCX(1,1) = XSCAX
      FOCX(2,1) = XLCAX
      FOCX(1,2) = XSCOX
      FOCX(2,2) = XLCOX
      FOCX(1,3) = XSCNX
      FOCX(2,3) = XLCNX
      TCPWX(1) = 2.D+0
      TCPWX(2) = 6.D+0
      VUCX = 12.D-8**3
      NCTX = 2
      NCHX = 3
      RHOHX = WMPUC*HCPP(1,4)
      DO 210 J = 1,NCHX
      DO 200 I = 1,NCTX
        RHOHX = RHOHX + FOCX(I,J)*TCPWX(I)*HCPP(1,J)
  200 CONTINUE
  210 CONTINUE
      RHOHX = 1.D+3*RHOHX/(AVGNX*VUCX)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of HYD_PC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ICE( TX,PLX,XLSX,TFPX,PIX )
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
!     Compute ice freezing point and ice pressure as a function of
!     temperature, pressure, and dissolved aqueous components.
!
!     Pure water freezing point as a function of pressure from
!     Wagner W., A. Saul, and A. Pruss.  1994.  "International equations
!     for the pressure along the melting and along the sublimation
!     curve of ordinary water substance."  Journal of Physical and 
!     Chemical Reference Data, 23:515-527.
!
!     Ice-aqueous radius of curvature from Brun, M., A. Lallemand,
!     J.F. Quinson, and C. Eyraud.  A new method for the simultaneous
!     determination of the size and shape of pores:
!     The thermoporometry."  Thermochimica Acta, 21:59-88.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 2 August 2005.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYST
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
      SUB_LOG(ISUB_LOG) = '/ICE'
!
!---  Ice I freezing point (273.16 - 251.165 K) from Wagner et al., 
!     1994, using Newton iteration and a starting guess of 0 C ---
!
      TNX = 273.16D+0
      PNX = 6.11657D+2
      TFPX = TABS
      THETAX = TFPX/TNX
   10 CONTINUE
      PIX = 1.D+0 - 0.626D+6*(1.D+0-(1.D+0/(THETAX**3))) + 
     &  0.197135D+6*(1.D+0-(THETAX**(21.2D+0)))
      DPIX = -1.878D+6/(THETAX**4) - 4.179262D+6*(THETAX**20.2D+0)
      FX = PIX*PNX-(PLX+PATM)
      DFX = DPIX*PNX
      DTHETAX = -FX/DFX
      THETAX = THETAX + DTHETAX
      IF( ABS(DTHETAX).GT.1.D-12 ) GOTO 10
      TFPX = THETAX*TNX - TABS
      CMSX = 1.D+3*XLSX/WTMS
      TFPX = TFPX - 1.858D+0*CMSX
!
!---  Ice-aqueous interfacial radius from Brun et al., 1977  ---
!
      IF( TX.LT.TFPX ) THEN
!        RADX = 1.D-9*((64.67D+0/(TFPX-TX)) + 0.57D+0)
!        RADX = 1.D-10*((6467.D+0/(TFPX-TX)) - 2.3D+0)
        RADX = 3.6D-9*(1.D+0 - (6.6D-1/LOG10((TX+TABS)/(TFPX+TABS))))
!        RADX = 1.0000D-7/(TFPX-TX)
        PIX = PLX + 2.D+0*SIG_IL/RADX
      ELSE
        PIX = PLX
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ICE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ITL_A( PX,TX,I_PX,I_TX,I_VX,NF )
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
!     This subroutine determines the lookup table indices for
!     computing pure CO2, CH4, or N2 properties.
!
!     I_VX = 0 solid
!     I_VX = 1 subcritical gas
!     I_VX = 2 subcritical liquid
!     I_VX = 3 supercritical gas
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 22 April 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
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
      INTEGER I_PX(2),I_TX(2),I_VX
      CHARACTER*25 CHMSGX(3)
      CHARACTER(64) :: SUBLOGX
!
!----------------------Data Statements---------------------------------!
!
      DATA CHMSGX / 'Solid Phase for CO2: T = ',
     &  'Solid Phase for CH4: T = ',
     &  'Solid Phase for N2: T =  '/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ITL_A'
!
!---  Convert pressure to MPa and temperature to degrees Kelvin  ---
!
      PMX = 1.D-6*PX
      TKX = TX + TABS
      PM_TPX = P_TA(IPTP(NF),NF)
      PM_CRX = P_TA(IPCR(NF),NF)
      TK_CRX = T_TA(IV_TA(IPCR(NF),NF),IPCR(NF),NF)
      TK_TPX = T_TA(IV_TA(IPTP(NF),NF),IPTP(NF),NF)
!
!---  Find pressure indices  ---
!
      IF( PMX.LT.P_TA(1,NF) ) PMX = P_TA(1,NF)
      IF( PMX.GT.P_TA(IP_TA(NF),NF) ) PMX = P_TA(IP_TA(NF),NF)
      IPLX = 1
      IPUX = IP_TA(NF)
   10 IF( IPUX-IPLX.GT.1 ) THEN
        IPM = (IPLX+IPUX)/2
        IF( (P_TA(IP_TA(NF),NF).GT.P_TA(1,NF)).EQV.
     &    (PMX.GT.P_TA(IPM,NF)) ) THEN
          IPLX = IPM
        ELSE
          IPUX = IPM
        ENDIF
        GOTO 10
      ENDIF
!
!---  Temperature above the critical point ---
!
      IF( TKX.GT.TK_CRX ) THEN
        ISLX = 1
        IELX = IT_TA(IPLX,NF)
        ISUX = 1
        IEUX = IT_TA(IPUX,NF)
        I_VX = 3
!
!---  Temperature below the critical point  ---
!
      ELSE
!
!---    Determine saturation-point temperature at the 
!       specified pressure, via linear interpolation  ---
!
        PLX = P_TA(IPLX,NF)
        PUX = P_TA(IPUX,NF)
        TLX = T_TA(IV_TA(IPLX,NF),IPLX,NF)
        TUX = T_TA(IV_TA(IPUX,NF),IPUX,NF)
        TK_SPX = ((PMX-PLX)/(PUX-PLX))*(TUX-TLX) + TLX
        IF( PMX.GE.PM_CRX ) TK_SPX = TK_CRX
!
!---    Determine melting-point temperature at the 
!       specified pressure, via linear interpolation  ---
!
        TLX = T_TA(1,IPLX,NF)
        TUX = T_TA(1,IPUX,NF)
        TK_MPX = ((PMX-PLX)/(PUX-PLX))*(TUX-TLX) + TLX
!
!---    Interpolation on gas side, searching from
!       the saturation line  ---
!
        IF( TKX.GE.TK_SPX ) THEN
          ISLX = IV_TA(IPLX,NF)
          IELX = IT_TA(IPLX,NF)
          ISUX = IV_TA(IPUX,NF)
          IEUX = IT_TA(IPUX,NF)
          I_VX = 1
!
!---    Pressure above critical pressure, search over all
!       temperatures  ---
!
        ELSEIF( PMX.GT.PM_CRX .AND. TKX.GE.TK_MPX ) THEN
          ISLX = 1
          IELX = IT_TA(IPLX,NF)
          ISUX = 1
          IEUX = IT_TA(IPUX,NF)
          I_VX = 2
!
!---    Interpolation on liquid side, searching to
!       the saturation line  ---
!
        ELSEIF( PMX.GT.PM_TPX .AND. TKX.GE.TK_MPX ) THEN
          ISLX = 1
          IELX = IV_TA(IPLX,NF)-1
          ISUX = 1
          IEUX = IV_TA(IPUX,NF)-1
          I_VX = 2
!
!---    Temperature less than melting point, solid  ---
!
        ELSE
          INDX = 4
          IMSGX = 0
          NMSGX = 0
          SUBLOGX = 'ITL_A'
          RLMSGX = TX
          CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX(NF),IMSGX,NMSGX,INDX )
        ENDIF
      ENDIF      
!
!---  Lower and upper pressure indices  ---
!
      I_PX(1) = IPLX
      I_PX(2) = IPUX
!
!---  Find lower-pressure temperature indices  ---
!
      IF( TKX.LT.T_TA(ISLX,IPLX,NF) ) THEN
        I_TX(1) = ISLX
      ELSEIF( TKX.GT.T_TA(IELX,IPLX,NF) ) THEN
        I_TX(1) = IELX-1
      ELSE
        ITLX = 1
        ITUX = IT_TA(IPLX,NF)
   20   IF( ITUX-ITLX.GT.1 ) THEN
          ITM = (ITLX+ITUX)/2
          IF( (T_TA(IELX,IPLX,NF).GT.T_TA(ISLX,IPLX,NF))
     &      .EQV.(TKX.GT.T_TA(ITM,IPLX,NF)) ) THEN
            ITLX = ITM
          ELSE
            ITUX = ITM
          ENDIF
          GOTO 20
        ENDIF
      I_TX(1) = ITLX
      ENDIF
!
!---  Find upper-pressure temperature indices  ---
!
      IF( TKX.LT.T_TA(ISUX,IPUX,NF) ) THEN
        I_TX(2) = ISUX
      ELSEIF( TKX.GT.T_TA(IEUX,IPUX,NF) ) THEN
        I_TX(2) = IEUX-1
      ELSE
        ITLX = 1
        ITUX = IT_TA(IPUX,NF)
   30   IF( ITUX-ITLX.GT.1 ) THEN
          ITM = (ITLX+ITUX)/2
          IF( (T_TA(IEUX,IPUX,NF).GT.T_TA(ISUX,IPUX,NF))
     &      .EQV.(TKX.GT.T_TA(ITM,IPUX,NF)) ) THEN
            ITLX = ITM
          ELSE
            ITUX = ITM
          ENDIF
          GOTO 30
        ENDIF
        I_TX(2) = ITLX
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ITL_A group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE KF_2PE( BTX,PPELX,PPEUX,PX,TCTX,TKX,ZMCX )
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
!     K-factors and beta within the two-phase pressure envelope.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 3 September 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYDT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
!      REAL*8 FKX(LHF_HT)
      REAL*8 ZMCX(LHF_HT)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/KF_2PE'
!
!---  Beta for the specified composition, using bicubic interpolation
!     over triangles across the CO2 (x-direction) and 
!     CH4 (y-direction) table of phase envelope data, and using 
!     bi-cubic spline interpolation on temperature and pressure, 
!     defined relative to the cricondenterm temperature, lower phase
!     envelope pressure, and upper phase envelope pressure
!     for the specified concentration  ---
!
      CALL BCIP_TPS( D2B2P_HT,B2P_HT,BTX,PPELX,PPEUX,PX,
     &  TCTX,TKX,ZMCX,NPEP_HT )
      BTX = MAX( MIN( BTX,1.D+0 ),0.D+0 )
!
!---  K-factors for the specified composition, using bicubic 
!     interpolation over triangles across the CO2 (x-direction) and 
!     CH4 (y-direction) table of phase envelope data, and using 
!     bi-cubic spline interpolation on temperature and pressure, 
!     defined relative to the cricondenterm temperature, lower phase
!     envelope pressure, and upper phase envelope pressure
!     for the specified concentration  ---
!
!      DO 10 NF = 1,NHF_HT
!        CALL BCIP_TPS( D2FK2P_HT(1,1,1,NF),FK2P_HT(1,1,1,NF),FKX(NF),
!     &    PPELX,PPEUX,PX,TCTX,TKX,ZMCX,NPEP_HT )
!   10 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of KF_2PE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE KF_LPE( TCTX,TKX,ZMCX,FKX )
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
!     K-factors at the lower pressure envelope.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 17 August 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYDT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 FKX(LHF_HT),ZMCX(LHF_HT)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/KF_LPE'
      DO 10 NF = 1,NHF_HT
        CALL BCIP_TS( D2FKLE_HT(1,1,NF),FKLE_HT(1,1,NF),FKX(NF),
     &    TCTX,TKX,TLE_HT,ZMCX,NPEP_HT )
   10 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of KF_LPE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE KF_UPE( TCTX,TKX,ZMCX,FKX )
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
!     K-factors at the upper pressure envelope.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 17 August 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYDT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 FKX(LHF_HT),ZMCX(LHF_HT)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/KF_UPE'
      DO 10 NF = 1,NHF_HT
        CALL BCIP_TS( D2FKUE_HT(1,1,NF),FKUE_HT(1,1,NF),FKX(NF),
     &    TCTX,TKX,TUE_HT,ZMCX,NPEP_HT )
   10 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of KF_UPE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE NICKALLS( CAX,CBX,CCX,CDX,R1X,R2X,R3X )
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
!     This subroutine calculates the roots of a cubic equation using
!     the Nickalls formulation.
!
!     Nickalls, R.W.D. 1993. "A new approach to solving the cubic:
!     Cardan's solution revealed."  Math. Gaz. 77:354-359.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 June 2010.
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/NICKALLS'
!
!---  Nickalls cubic polynomial solver  ---
!
      XNX = -CBX/(3.D+0*CAX)
      YNX = CAX*(XNX**3) + CBX*(XNX**2) + CCX*XNX + CDX
      YN2X = YNX**2
      D2X = ((CBX**2)-(3.D+0*CAX*CCX))/((3.D+0*CAX)**2)
      IF( D2X.LE.0.D+0 ) THEN
        HX = 0.D+0
      ELSE
        DX = SQRT(D2X)
        HX = -2.D+0*(DX**3)
      ENDIF
      H2X = 4.D+0*(CAX**2)*(D2X**3)
!
!---  YN2X > H2X, 1 real root  ---
!
      IF( YN2X-H2X.GT.EPSL ) THEN
        R1X = (5.D-1/CAX)*(-YNX+SQRT(YN2X-H2X))
        R2X = (5.D-1/CAX)*(-YNX-SQRT(YN2X-H2X))
        IF( R1X.LT.0.D+0 .AND. R2X.LT.0.D+0 ) THEN
          VX = XNX + (ABS(R1X)**(1.D+0/3.D+0))
     &      + (ABS(R2X)**(1.D+0/3.D+0))
          CHKX = CAX*(VX**3) + CBX*(VX**2) + CCX*VX + CDX
        ENDIF
        R3X = XNX + SIGN((ABS(R1X)**(1.D+0/3.D+0)),R1X)
     &    + SIGN((ABS(R2X)**(1.D+0/3.D+0)),R2X)
        R1X = R3X
        R2X = R3X
!
!---  YN2X < H2X, 3 distinct real roots  ---
!
      ELSEIF( YN2X-H2X.LT.-EPSL ) THEN
        THETAX = (ACOS(YNX/HX))/3.D+0
        R1X = XNX + 2.D+0*DX*COS(THETAX)
        R2X = XNX + 2.D+0*DX*COS(((2.D+0*GPI)/3.D+0)+THETAX)
        R3X = XNX + 2.D+0*DX*COS(((4.D+0*GPI)/3.D+0)+THETAX)
!
!---  YN2X = H2X, 3 real roots (two or three equal roots)  ---
!
      ELSE
!
!---    HX /= 0 (two equal roots)  ---
!
        IF( ABS(HX)/EPSL.GT.EPSL ) THEN
          DX = YNX/(2.D+0*CAX)
          DX = SIGN((ABS(DX)**(1.D+0/3.D+0)),DX)
          R1X = XNX + DX
          R2X = XNX + DX
          R3X = XNX - 2.D+0*DX
!
!---    HX = 0 (three equal roots)  ---
!
        ELSE
          R1X = XNX
          R2X = XNX
          R3X = XNX
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PRC_Z group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      FUNCTION POLY( TX,WX,LX )
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
!     Polynomial function.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 24 November 2004.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 WX(LX)
!
!----------------------Executable Lines--------------------------------!
!
      POLY = WX(1)
      DO 100 I = 2,LX
        POLY = POLY + WX(I)*(TX**(I-1))
  100 CONTINUE
!
!---  End of POLY group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PTL_A( PX,TX,VAR_TA,VARX,I_PX,I_TX,NF )
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
!     This subroutine computes CO2 properties as a function of
!     pressure and temperature, using bilinear interpolation of
!     tabular data.
!
!     I_PX(1) - lower pressure index
!     I_PX(2) - upper pressure index
!     I_TX(1) - lower temperature index at lower pressure
!     I_TX(2) - lower temperature index at upper pressure
!     VAR_TA(LT_TA,LP_TA,LNHC) - pressure, temperature ordered 
!       array of property values
!     VARX - interpolated property value
!     DVARX - partial derivative of the property value with respect
!             to pressure
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 April 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VAR_TA(LT_TA,LP_TA,LNHC)
      INTEGER I_PX(2),I_TX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PTL_A'
!
!---  Convert pressure to MPa and temperature to degrees Kelvin  ---
!
      PMX = 1.D-6*PX
      TKX = TX + TABS
!
!---  Interpolation parameters  ---
!
      VPLTLX = VAR_TA(I_TX(1),I_PX(1),NF)
      VPLTUX = VAR_TA(I_TX(1)+1,I_PX(1),NF)
      VPUTLX = VAR_TA(I_TX(2),I_PX(2),NF)
      VPUTUX = VAR_TA(I_TX(2)+1,I_PX(2),NF)
      PLX = P_TA(I_PX(1),NF)
      PUX = P_TA(I_PX(2),NF)
      TLPLX = T_TA(I_TX(1),I_PX(1),NF)
      TUPLX = T_TA(I_TX(1)+1,I_PX(1),NF)
      TLPUX = T_TA(I_TX(2),I_PX(2),NF)
      TUPUX = T_TA(I_TX(2)+1,I_PX(2),NF)
!
!---  Bilinear interpolation scheme  ---
!
      VPLX = ((VPLTUX-VPLTLX)*(TKX-TLPLX)/(TUPLX-TLPLX)) + VPLTLX
      VPUX = ((VPUTUX-VPUTLX)*(TKX-TLPUX)/(TUPUX-TLPUX)) + VPUTLX
      IF( ABS(PUX-PLX).LT.EPSL ) THEN
        DVARX = 0.D+0
        VARX = VPLX
      ELSE
        DVARX = (VPUX-VPLX)/(PUX-PLX)
        VARX = (DVARX*(PMX-PLX)) + VPLX
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PTL_A group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PRC_T( AIX,TKX )
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
!     Peng-Robinson temperature-dependent constants for mixtures.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 24 August 2009
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYDT
      USE FDVH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 AIX(LHF_HT)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PRC_T'
      DO 100 NF = 1,NHF_HT
        TRX = TKX/HCPP(3,NF)
        AIX(NF) = HCPP(13,NF)*
     &    (1.D+0 + HCPP(15,NF)*(1.D+0 - SQRT(TRX)))**2
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PRC_T group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PRC_X( AX,AIX,BX,XMCX )
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
!     Peng-Robinson composition-dependent constants for mixtures.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 March 2012
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYDT
      USE FDVH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 AIX(LHF_HT),XMCX(LHF_HT)
      REAL*8 FPRKIJ(3,3)
!
!----------------------Data Statements---------------------------------!
!
      DATA FPRKIJ / 0.D+0, 0.109267D+0, 0.0195133D+0, 0.109267D+0, 
     &  0.D+0, 0.0262942D+0, 0.0195133D+0, 0.0262942D+0, 0.D+0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PRC_X'
      AX = 0.D+0
      BX = 0.D+0
      DO 120 IHF = 1,NHF_HT
        DO 110 JHF = 1,NHF_HT
          AIJX = SQRT(AIX(IHF)*AIX(JHF))*(1.D+0-FPRKIJ(IHF,JHF))
          AX = AX + XMCX(IHF)*XMCX(JHF)*AIJX
  110   CONTINUE
        BX = BX + XMCX(IHF)*HCPP(14,IHF)
  120 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PRC_X group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PRC_Z( AX,BX,CAX,CBX,PX,TKX,ZX,INDX )
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
!     Peng-Robinson compressibility.
!
!     INDX = 1 Gas
!     INDX = 2 Liquid
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 March 2012
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PRC_Z'
      CAX = AX*PX/((RCU*TKX)**2)
      CBX = BX*PX/(RCU*TKX)
!
!---  Peng-Robinson cubic function for Z  ---
!
      C3X = 1.D+0
      C2X = -(1.D+0-CBX)
      C1X = (CAX - 3.D+0*(CBX**2) - 2.D+0*CBX)
      C0X = -(CAX*CBX - (CBX**2) - (CBX**3))
!
!---  Direct cubic equation solve  ---
!
      CALL NICKALLS( C3X,C2X,C1X,C0X,R1X,R2X,R3X )
!
!---  Gas  ---
!
      IF( INDX.EQ.1 ) THEN
        ZX = MAX(R1X,R2X,R3X)
!
!---  Liquid  ---
!
      ELSEIF( INDX.EQ.2 ) THEN
        R1X = MAX( R1X,1.000001D+0*CBX )
        R2X = MAX( R2X,1.000001D+0*CBX )
        R3X = MAX( R3X,1.000001D+0*CBX )
        ZX = MIN( R1X,R2X,R3X )
!
!---  Ideal gas  ---
!
      ELSE
        ZX = 1.D+0
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PRC_Z group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE P_STATE( TCRX,TCTX,TKX,PPELX,PPEUX,PX,ZMCX,NP_HT )
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
!     Determine the phase state and concentration data sets
!
!     NP_HT = 0 - gas, temperature > cricondentherm temperature
!     NP_HT = 1 - liquid, temperature < critical temperature,
!                         pressure > upper phase envelope pressure
!     NP_HT = 2 - two-phase, temperature < cricondentherm temperature,
!                            pressure ≥ lower phase envelope pressure,
!                            pressure ≤ upper phase envelope pressure
!     NP_HT = 3 - gas, temperature < cricondentherm temperature
!                      temperature > critical temperature,
!                      pressure > upper phase envelope pressure                  
!     NP_HT = 4 - gas, temperature < cricondentherm temperature
!                      pressure < lower phase envelope pressure
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 03 August 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYDT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 ZMCX(LHF_HT)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/P_STATE'
      NP_HT = -1
      TOLX = 1.D-12
!
!---  If CO2 concentration is less than table minimum, set
!     phase condition to gas  ---
!
      IF( ZMCX(1).LE.ZPC_HT(1,1) ) THEN
        NP_HT = 0
        TCTX = 0.D+0
        TCRX = 0.D+0
        PPELX = 0.D+0
        PPEUX = 1.D+20
        GOTO 100
      ENDIF
!
!---  Cricondenterm temperature for the specified composition  ---
!
      CALL BCIP_SC( TCT_HT,TCTX,ZMCX )
      TCTX = MAX( TCTX,0.D+0 )
!
!---  Critical temperature for the specified composition  ---
!
      CALL BCIP_SC( TCR_HT,TCRX,ZMCX )
      TCRX = MAX( TCRX,0.D+0 )
!
!---  Gas, temperature > cricondentherm temperature  ---
!
      IF( TKX.GT.TCTX ) THEN
        NP_HT = 0
        PPELX = 0.D+0
        PPEUX = 1.D+20
        GOTO 100
!
!---  Temperature less than cricondenterm temperature  ---
!
      ELSE
!
!---    Lower phase envelope pressure for the specified composition
!       and temperature  ---
!
        CALL BCIP_TS( D2PLE_HT,PLE_HT,PPELX,TCTX,TKX,TLE_HT,ZMCX,
     &    NPEP_HT )
!
!---    Upper phase envelope pressure for the specified composition
!       and temperature  ---
!
        CALL BCIP_TS( D2PUE_HT,PUE_HT,PPEUX,TCTX,TKX,TUE_HT,ZMCX,
     &    NPEP_HT )
      ENDIF
!
!---  Gas, temperature < cricondentherm temperature
!     pressure < lower phase envelope pressure  ---
!
      IF( PX.LE.PPELX ) THEN
        NP_HT = 4
!
!---  Gas, temperature < cricondentherm temperature
!     temperature > critical temperature,
!     pressure > upper phase envelope pressure  ---
!
      ELSEIF( PX.GE.PPEUX .AND. TKX.GT.TCRX ) THEN
        NP_HT = 3
!
!---  Liquid, temperature < critical temperature,
!     pressure > upper phase envelope pressure  ---
!
      ELSEIF( PX.GE.PPEUX .AND. TKX.LE.TCRX ) THEN
        NP_HT = 1
!
!---  Two-phase, temperature < cricondentherm temperature,
!     pressure ≥ lower phase envelope pressure,
!     pressure ≤ upper phase envelope pressure  ---
!
      ELSE
        NP_HT = 2
      ENDIF
  100 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of P_STATE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDPF_HYDT_KE
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
!     H2O-NaCl-CO2 Equation of state.  This subroutine reads tables of
!     density, enthalpy, internal energy, and fugacity coefficient of
!     pure noncondensible gas as a function of temperature and
!     pressure.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 27 March 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
      USE HYDT
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      LOGICAL FCHK
      CHARACTER*512 CHDUM
      CHARACTER*64 FDUM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDPF_HYDT_KE'
      TOLX = 1.D-6
!
!---  Number of hydrate formers  ---
!
      NHF_HT = 3
!
!---  Check for phase-envelope file  ---
!
      INQUIRE( FILE='co2_ch4_n2_pe.dat', FORM=FDUM, EXIST=FCHK )
      IF( .NOT.FCHK ) THEN
        INDX = 4
        CHMSG = 'Missing Phase-Envelope File: co2_ch4_n2_pe.dat'
        CALL WRMSGS( INDX )
      ELSEIF( FDUM.EQ.'unformatted' ) THEN
        INDX = 4
        CHMSG = 'Unformatted Phase-Envelope File: co2_ch4_n2_pe.dat'
        CALL WRMSGS( INDX )
      ENDIF
      OPEN(UNIT=26, FILE='co2_ch4_n2_pe.dat', STATUS='OLD', 
     &    FORM='FORMATTED')
!
!---  Read Id tag line ---
!
      READ(26,'(A)') CHDUM
!
!---  Loop over the number of phase equilibria table compositions  ---
!
      DO 200 NC = 1,LCP_HT
!
!---    Read data header line  ---
!
        READ(26,'(A)') CHDUM
!
!---    Read hydrate former mole fractions  ---
!
        READ(26,'(3(1PE16.9,1X))') (ZPC_HT(NF,NC),NF=1,NHF_HT)
!
!---    Table pointer  ---
!
        IX = INT( (ZPC_HT(1,NC)/DZPC_HT)+EPSL ) - 1
        JX = INT( (ZPC_HT(2,NC)/DZPC_HT)+EPSL ) + 1
        IZP_HT(IX,JX) = NC
!
!---    Read data header line  ---
!
        READ(26,'(A)') CHDUM
!
!---    Read composition critical temperature, critical pressure,
!       cricondenbar temperature, cricondenbar pressure,
!       cricondenterm temperature, and cricondenterm pressure  ---
!
        READ(26,'(6(1PE16.9,1X))') TCR_HT(NC),PCR_HT(NC),
     &    TCB_HT(NC),PCB_HT(NC),TCT_HT(NC),PCT_HT(NC)
!
!---    Read data header line  ---
!
        READ(26,'(A)') CHDUM
!
!---    Read number of phase envelope table points  ---
!
        READ(26,'(I6)') NPEP_HT(NC)
!
!---    Read data header line  ---
!
        READ(26,'(A)') CHDUM
!
!---    Loop over the number of phase envelope table points  ---
!
        DO 110 NE = 1,NPEP_HT(NC)
!
!---      Read temperature for the lower phase envelope,
!         pressure for the lower phase envelope,
!         second derivative of pressure for lower phase envelope
!         K-factors for the lower phase envelope,
!         second derivative of K-factors for lower phase envelope  ---
!
          READ(26,'(9(1PE16.9,1X))') TLE_HT(NE,NC),PLE_HT(NE,NC),
     &      D2PLE_HT(NE,NC),(FKLE_HT(NE,NC,NF),NF=1,NHF_HT),
     &      (D2FKLE_HT(NE,NC,NF),NF=1,NHF_HT)
  110   CONTINUE
!
!---    Read data header line  ---
!
        READ(26,'(A)') CHDUM
!
!---    Read number of phase envelope table points  ---
!
        READ(26,'(I6)') NPEP_HT(NC)
!
!---    Read data header line  ---
!
        READ(26,'(A)') CHDUM
!
!---    Loop over the number of phase envelope table points  ---
!
        DO 120 NE = 1,NPEP_HT(NC)
!
!---      Read temperature for the upper phase envelope,
!         pressure for the upper phase envelope,
!         second derivative of pressure for upper phase envelope
!         K-factors for the upper phase envelope,
!         second derivative of K-factors for upper phase envelope  ---
!
          READ(26,'(9(1PE16.9,1X))') TUE_HT(NE,NC),PUE_HT(NE,NC),
     &      D2PUE_HT(NE,NC),(FKUE_HT(NE,NC,NF),NF=1,NHF_HT),
     &      (D2FKUE_HT(NE,NC,NF),NF=1,NHF_HT)
  120   CONTINUE
!
!---    Read data header line  ---
!
        READ(26,'(A)') CHDUM
!
!---    Read number of two-phase temperature points  ---
!
        READ(26,'(I6)') NTP_HT(NC)
!
!---    Read data header line  ---
!
        READ(26,'(A)') CHDUM
!
!---    Read two-phase temperature points  ---
!
        READ(26,'(10(1PE16.9,1X))') (T2P_HT(NT,NC),NT=1,NTP_HT(NC))
!
!---    Read data header line  ---
!
        READ(26,'(A)') CHDUM
!
!---    Read number of two-phase pressure points  ---
!
        READ(26,'(I6)') NPP_HT(NC)
!
!---    Read data header line  ---
!
        READ(26,'(A)') CHDUM
!
!---    Read two-phase pressure points  ---
!
        READ(26,'(10(1PE16.9,1X))') (P2P_HT(NP,NC),NP=1,NPP_HT(NC))
!
!---    Read data header line  ---
!
        READ(26,'(A)') CHDUM
!
!---    Read two-phase beta values at temperature and 
!       pressure points  ---
!
        READ(26,'(10(1PE16.9,1X))') ((B2P_HT(NT,NP,NC),
     &    NT=1,NTP_HT(NC)),NP=1,NPP_HT(NC))
!
!---    Read data header line  ---
!
        READ(26,'(A)') CHDUM
!
!---    Read two-phase second derivative beta values at temperature and 
!       pressure points  ---
!
        READ(26,'(10(1PE16.9,1X))') ((D2B2P_HT(NT,NP,NC),
     &    NT=1,NTP_HT(NC)),NP=1,NPP_HT(NC))
!
!---    Read data header line  ---
!
        READ(26,'(A)') CHDUM
!
!---    Read two-phase K-factor values at temperature and 
!       pressure points for each hydrate former  ---
!
        READ(26,'(10(1PE16.9,1X))') (((FK2P_HT(NT,NP,NC,NF),
     &    NT=1,NTP_HT(NC)),NP=1,NPP_HT(NC)),NF=1,NHF_HT)
!
!---    Read data header line  ---
!
        READ(26,'(A)') CHDUM
!
!---    Read two-phase second derivative K-factor values at 
!       temperature and pressure points  ---
!
        READ(26,'(10(1PE16.9,1X))') (((D2FK2P_HT(NT,NP,NC,NF),
     &    NT=1,NTP_HT(NC)),NP=1,NPP_HT(NC)),NF=1,NHF_HT)
  200 CONTINUE
!
!---  Close phase-envelope file  ---
!
      CLOSE(26)
!
!---  Create mirror image of phase envelope table  ---
!
      NZP_HT = 9
      DO 204 I = 2,NZP_HT
        DO 202 J = NZP_HT+2-I,NZP_HT
          IZP_HT(I,J) = IZP_HT(NZP_HT+1-J,NZP_HT+1-I)
  202   CONTINUE
  204 CONTINUE
!
!---  Check for CO2-CH4-N2 hydrate-equilibrium file  ---
!
      INQUIRE( FILE='hyd_co2_ch4_n2_eq.dat', FORM=FDUM, EXIST=FCHK )
      IF( .NOT.FCHK ) THEN
        INDX = 4
        CHMSG = 'Missing Phase-Envelope File: hyd_co2_ch4_n2_eq.dat'
         CALL WRMSGS( INDX )
      ELSEIF( FDUM.EQ.'unformatted' ) THEN
        INDX = 4
        CHMSG = 'Unformatted Phase-Envelope File: ' // 
     &    'hyd_co2_ch4_n2_eq.dat'
        CALL WRMSGS( INDX )
      ENDIF
      OPEN(UNIT=26, FILE='hyd_co2_ch4_n2_eq.dat', STATUS='OLD', 
     &    FORM='FORMATTED')
!
!---  Read Id tag line ---
!
      READ(26,'(A)') CHDUM
!
!---  Loop over the number of hydrate former compositions  ---
!
      DO 300 NC = 1,LCH_HT
!
!---    Read data header line  ---
!
        READ(26,'(A)') CHDUM
!
!---    Read hydrate former mole fractions  ---
!
        READ(26,'(3(1PE12.5,1X))') (ZHC_HT(NF,NC),NF=1,NHF_HT)
!
!---    Table pointer  ---
!
        IX = INT( (ZHC_HT(1,NC)/DZHC_HT)+EPSL ) + 1
        JX = INT( (ZHC_HT(2,NC)/DZHC_HT)+EPSL ) + 1
        IZH_HT(IX,JX) = NC
!
!---    Read number of hydrate equilibrium table points  ---
!
        READ(26,'(I6)') NHEP_HT(NC)
!
!---    Loop over the number of hydrate equilibrium table points  ---
!
        DO 210 NE = 1,NHEP_HT(NC)
!
!---      Read temperature (K), pressure (bar), small cage occupancy
!         for CO2, CH4 and N2, and large cage occupancy for CO2, CH4,
!         and N2 ---
!
          READ(26,'(8(1PE12.5,1X))') THE_HT(NE,NC),PHE_HT(NE,NC),
     &      XSCA_HT(NE,NC),XSCO_HT(NE,NC),XSCN_HT(NE,NC),
     &      XLCA_HT(NE,NC),XLCO_HT(NE,NC),XLCN_HT(NE,NC)
          PHE_HT(NE,NC) = 1.D+5*PHE_HT(NE,NC)
  210   CONTINUE
!
!---    Smooth equilibrium curves at the freezing point of water  ---
!
        PHE_HT(61,NC) = 5.D-1*(PHE_HT(60,NC)+PHE_HT(62,NC))
        XSCA_HT(61,NC) = 5.D-1*(XSCA_HT(60,NC)+XSCA_HT(62,NC))
        XSCO_HT(61,NC) = 5.D-1*(XSCO_HT(60,NC)+XSCO_HT(62,NC))
        XSCN_HT(61,NC) = 5.D-1*(XSCN_HT(60,NC)+XSCN_HT(62,NC))
        XLCA_HT(61,NC) = 5.D-1*(XLCA_HT(60,NC)+XLCA_HT(62,NC))
        XLCO_HT(61,NC) = 5.D-1*(XLCO_HT(60,NC)+XLCO_HT(62,NC))
        XLCN_HT(61,NC) = 5.D-1*(XLCN_HT(60,NC)+XLCN_HT(62,NC))
!
!---    Cubic spline of temperature versus pressure  ---
!
        CALL SPLINE( PHE_HT(1,NC),THE_HT(1,NC),NHEP_HT(NC),
     &    THE2P_HT(1,NC) )
!
!---    Cubic spline of CO2 small cage occupancy versus temperature  ---
!
        CALL SPLINE( THE_HT(1,NC),XSCA_HT(1,NC),NHEP_HT(NC),
     &    XSCA2T_HT(1,NC) )
!
!---    Cubic spline of CH4 small cage occupancy versus temperature  ---
!
        CALL SPLINE( THE_HT(1,NC),XSCO_HT(1,NC),NHEP_HT(NC),
     &    XSCO2T_HT(1,NC) )
!
!---    Cubic spline of N2 small cage occupancy versus temperature  ---
!
        CALL SPLINE( THE_HT(1,NC),XSCN_HT(1,NC),NHEP_HT(NC),
     &    XSCA2T_HT(1,NC) )
!
!---    Cubic spline of CO2 large cage occupancy versus temperature  ---
!
        CALL SPLINE( THE_HT(1,NC),XLCA_HT(1,NC),NHEP_HT(NC),
     &    XLCA2T_HT(1,NC) )
!
!---    Cubic spline of CH4 large cage occupancy versus temperature  ---
!
        CALL SPLINE( THE_HT(1,NC),XLCO_HT(1,NC),NHEP_HT(NC),
     &    XLCO2T_HT(1,NC) )
!
!---    Cubic spline of N2 large cage occupancy versus temperature  ---
!
        CALL SPLINE( THE_HT(1,NC),XLCN_HT(1,NC),NHEP_HT(NC),
     &    XLCN2T_HT(1,NC) )
!
!---    Cubic spline of pressure versus temperature  ---
!
        CALL SPLINE( THE_HT(1,NC),PHE_HT(1,NC),NHEP_HT(NC),
     &    PHE2T_HT(1,NC) )
  300 CONTINUE
!
!---  Close hydrate equilibrium file  ---
!
      CLOSE(26)
!
!---  Create mirror image of hydrate equilibrium table  ---
!
      NZH_HT = 11
      DO 304 I = 2,NZH_HT
        DO 302 J = NZH_HT+2-I,NZH_HT
          IZH_HT(I,J) = IZH_HT(NZH_HT+1-J,NZH_HT+1-I)
  302   CONTINUE
  304 CONTINUE
!
!---  Check for inhibitor correlation file  ---
!
      INQUIRE( FILE='inhibitor.dat', FORM=FDUM, EXIST=FCHK )
      IF( .NOT.FCHK ) THEN
        INDX = 4
        CHMSG = 'Missing Inhibitor Correlation File: inhibitor.dat'
        CALL WRMSGS( INDX )
      ELSEIF( FDUM.EQ.'unformatted' ) THEN
        INDX = 4
        CHMSG = 'Unformatted Inhibitor Correlation File: inhibitor.dat'
        CALL WRMSGS( INDX )
      ENDIF
      OPEN(UNIT=26, FILE='inhibitor.dat', STATUS='OLD', 
     &  FORM='FORMATTED')
!
!---  Read inhibitor correlation file  ---
!
      DO INHX = 1,LINH
        READ(26,'(A10,7(E10.3))') INHNM(INHX),(CINH(M,INHX),M=1,6),
     &    XLIMX(INHX)
      ENDDO
      CLOSE(26)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDPF_HYDT_KE group  ---
!
      RETURN
      END


!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFT_L( TX,XLSX,SFTLX )
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
!     This subroutine calculates the surface tension of aqueous
!     solutions of sodium chloride as a function of temperature
!     and salt concentration.
!
!     Abramzon, A.A., and R.D. Gaukhberg.  1993.  "Surface tension of
!     salt solutions."  Russian Journal of Applied Chemistry,
!     66(6):1139-1146.
!
!     Lide, D.R., and H.V. Kehiaian.  1994.  CRC Handbook of
!     Thermophysical and Thermochemical Data.  CRC Press, Boca Raton,
!     Florida.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 May 2002.
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SFT_L'
!
!---  Convert mass fraction to molality  ---
!
      GLSX = 1.D+3*XLSX/(WTMS*(1.D+0-XLSX))
!
!---  Pure water vapor surface tension as a function of temperature
!     by D. Lide and H. Kehiaian
!
      SFTWX = 1.D-3*(75.6592D+0 - 1.40959D-1*TX - 2.66317D-4*(TX**2))
!
!---  Function by A. Abramzon and R. Gaukhberg  ---
!
      SFTLX = SFTWX + 1.57D-3*GLSX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFT_L group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SOL_LA( TX,FUGAX,XLSX,XMLAX )
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
!     Solubility of CO2 in NaCl aqueous solutions as a function
!     of temperature, salt mass fraction, and CO2 fugacity.
!
!     Battistelli, A., C. Claudio, and K. Pruess.  1997.  "The simulator
!     TOUGH2/EWASG for modelling geothermal reservoirs with brines and
!     noncondensible gas."  Geothermics, 26(4): 437-464.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 6 August 2012.
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SOL_LA'
!
!---  Henry's coefficient for NaCl aqueous solutions as a function
!     of temperature and salt molality  ---
!
      CALL HC_LA( TX,XLSX,HCX,DHCX )
      XMLAX = MIN( MAX( (FUGAX/HCX),0.D+0 ),1.D+0 )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SOL_LA group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SOL_LS( TX,XLSMX )
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
!     This subroutine calculates the mass fraction of NaCl salt in
!     saturated aqueous solutions.
!
!     McKibbin, R., and A. McNabb.  1993.  "Modeling the phase
!     boundaries and fluid properties of the system H2O-NaCl at high
!     temperatures and pressures."  Proceedings 15th NZ Geothermal
!     Workshop, University of Auckland, New Zealand.
!
!     Temperature range 0-800 C.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 2 April 2002.
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
      REAL*8 SFX(3)
!
!----------------------Data Statements---------------------------------!
!
      DATA SFX / 2.6218D-1, 7.2D-5, 1.06D-6 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SOL_LS'
!
!---  Maximum NaCl saturation  ---
!
      XLSMX = SFX(1) + SFX(2)*TX + SFX(3)*(TX**2)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SOL_LS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SOL_LN( TX,FUGNX,XLSX,XMLNX )
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
!     Solubility of N2 in NaCl aqueous solutions as a function
!     of temperature, salt mass fraction, and CO2 fugacity.
!
!     Battistelli, A., C. Claudio, and K. Pruess.  1997.  "The simulator
!     TOUGH2/EWASG for modelling geothermal reservoirs with brines and
!     noncondensible gas."  Geothermics, 26(4): 437-464.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 6 August 2012.
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SOL_LN'
!
!---  Henry's coefficient for NaCl aqueous solutions as a function
!     of temperature and salt molality  ---
!
      CALL HC_LN( TX,XLSX,HCX,DHCX )
      XMLNX = MIN( MAX( (FUGNX/HCX),0.D+0 ),1.D+0 )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SOL_LN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SOL_LO( TX,FUGOX,XLSX,XMLOX )
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
!     Solubility of CH4 in NaCl aqueous solutions as a function
!     of temperature, salt mass fraction, and CH4 fugacity.
!
!     Sloan, E. D., Jr.  1997.  Clathrate Hydrates of Natural Gases.
!     Marcel Dekker, Inc., New York.
!
!     Battistelli, A., C. Claudio, and K. Pruess.  1997.  "The simulator
!     TOUGH2/EWASG for modelling geothermal reservoirs with brines and
!     noncondensible gas."  Geothermics, 26(4): 437-464.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 6 August 2012.
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SOL_LO'
!
!---  Henry's coefficient for NaCl aqueous solutions as a function
!     of temperature and salt molality  ---
!
      CALL HC_LO( TX,XLSX,HCX,DHCX )
      XMLOX = MIN( MAX( (FUGOX/HCX),0.D+0 ),1.D+0 )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SOL_LO group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SP_A( TX,PVAX )
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
!     Sublimation-vapor pressure of CO2.
!
!     Span, R., and W. Wagner.  1996.  A New Equation of State for
!     Carbon Dioxide Covering the Fluid Region from the Triple-Point
!     to 1100 K at Pressures up to 800 MPa.  J. Phys. Chem. Ref. Data
!     25(6):1509-1588.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 2 September 2004.
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
      REAL*8 A_C(3),B_C(4)
!
!----------------------Data Statements---------------------------------!
!
      DATA A_C / -14.740846D+0,2.4327015D+0,-5.3061778D+0 /
      DATA B_C / -7.0602087D+0,1.9391218D+0,-1.6463597D+0,
     &  -3.2995634D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SP_A'
!
!---  Absolute temperature  ---
!
      TKX = TX + TABS
!
!---  Triple-point and critical temperature  ---
!
      TTX = 216.592D+0
      TCX = 304.1282D+0
!
!---  Sublimation pressure  ---
!
      IF( TKX.LT.TTX ) THEN
        PTX = 0.51795D+6
        TRX = TKX/TTX
        PVAX = PTX*EXP((TTX/TKX)*(A_C(1)*(1.D+0-TRX) + 
     &    A_C(2)*((1.D+0-TRX)**1.9D+0) + A_C(3)*((1.D+0-TRX)**2.9D+0)))
!
!---  Vapor pressure  ---
!
      ELSE
        PCX = 7.3773D+6
        TRX = TKX/TCX
        IF( TKX.LT.TCX ) THEN
          PVAX = PCX*EXP((TCX/TKX)*(B_C(1)*(1.D+0-TRX) + 
     &      B_C(2)*((1.D+0-TRX)**1.5D+0) + 
     &      B_C(3)*((1.D+0-TRX)**2.0D+0) +
     &      B_C(4)*((1.D+0-TRX)**4.0D+0)))
        ELSE
          PVAX = PCX
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SP_A group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SP_B( TX,XLSX,PSBX )
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
!     This subroutine calculates the saturation pressure of NaCl brine
!     as a function of temperature and salt concentration.
!
!     tx - temperature, C
!     xlsx - mass fraction of NaCl salt
!     psbx - saturation pressure of brine, Pa
!
!     Haas, J.L., Jr.  1976.  Physical Properties of the Coexisting
!     Phases and Thermochemical Properties of the H2O Component in
!     Boiling NaCl Solutions, U.S. Geological Survey Bulletin, 1421-A,
!     United States Government Printing Office, Washington.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 1 April 2002.
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
      REAL*8 SAX(3),SBX(5)
!
!----------------------Data Statements---------------------------------!
!
      DATA SAX / 5.93582D-6, -5.19386D-5, 1.23156D-5 /
      DATA SBX / 1.15420D-6, 1.41254D-7, -1.92476D-8, -1.70717D-9,
     &  1.05390D-10 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SP_B'
!
!---  Convert temperature to Kelvin and mass fraction to molality  ---
!
      TKX = TX + TABS
      GLSX = 1.D+3*XLSX/(WTMS*(1.D+0-XLSX))
!
!---  Concentration dependent coefficients  ---
!
      AX = 1.D+0
      DO 10 I = 1,3
        AX = AX + SAX(I)*(GLSX**I)
   10 CONTINUE
      BX = 0.D+0
      DO 20 I = 1,5
        BX = BX + SBX(I)*(GLSX**I)
   20 CONTINUE
!
!---  Temperature dependent coefficient  ---
!
      CX = 1.D+0/(AX + BX*TKX)
!
!---  Equivalent pure water temperature  ---
!
      TWX = EXP( CX*LOG(TKX) ) - TABS
      CALL SP_W( TWX,PSBX )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SP_B group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SP_O( TX,PVOX )
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
!     Vapor pressure of CH4.
!
!     Setzmann, U., and W. Wagner.  1991.  A New Equation of State and
!     Tables of Thermodynamic Properties for Methane Covering the
!     Range from the Melting Line to 625 K at Pressures up to 1000 MPa.
!     J. Phys. Chem. Ref. Data 20(6):1061-1155.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 8 October 2004.
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
      REAL*8 A_C(4)
!
!----------------------Data Statements---------------------------------!
!
      DATA TCROX,PCROX / 190.564D+0,4.5992D+6 /
      DATA A_C / -6.036219D+0,1.409353D+0,-0.4945199D+0,-1.443048D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SP_O'
!
!---  Absolute temperature  ---
!
      TKX = TX + TABS
      TRX = TKX/TCROX
      THX = 1.D+0 - TRX
!
!---  Vapor pressure  ---
!
      PVOX = PCROX*EXP( (TCROX/TKX)*(A_C(1)*THX + A_C(2)*(THX**1.5D+0)
     &  + A_C(3)*(THX**2) + A_C(4)*(THX**4.5D+0)) )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SP_O group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SP_W( TX,PSWX )
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
!     Saturation pressure (Pa) of pure water as a function of
!     temperature.
!
!     Meyer, C.A., R.B. McClintock, G.J. Silvestri, and R.C. Spencer
!     1993.  ASME Steam Tables, The American Society of Mechanical
!     Engineers, New York.
!
!     Buck, A. L.  1981.  "New equations for computing vapor pressure
!     and enhancement factor."  Journal of Applied Meterology, 
!     20:1527-1532.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 March 2002.
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
      REAL*8 K(9),PSWX,THETAX,THETAY,TX
      INTEGER I
!
!----------------------Data Statements---------------------------------!
!
      DATA K / -7.691234564D+0, -2.608023696D+1, -1.681706546D+2,
     &  6.423285504D+1, -1.189646225D+2, 4.167117320D+0,
     &  2.097506760D+1, 1.D+9, 6.D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SP_W'
!
!---  Merge functions above and below freezing point  ---
!
      VARX = MAX( MIN( 5.D-1*(1.D+0+TX),1.D+0 ), 0.D+0 )
!
!---  The K-function (saturation line)  ---
!
      THETAX = (TX+TABS)/TCRW
      THETAY = (1.D+0-THETAX)
      PSWX = 0.D+0
      DO 10 I = 1,5
        PSWX = PSWX + K(I)*(THETAY**I)
   10 CONTINUE
      PSWX = PSWX/((1.D+0 + K(6)*THETAY + K(7)*(THETAY**2))*THETAX)
      PSWX = PSWX - THETAY/(K(8)*(THETAY**2)+K(9))
      PSWX = VARX*EXP(PSWX)*PCRW
!
!---  Buck formulation  ---
!
      PSWX = PSWX + (1.D+0-VARX)*6.1121D+2*EXP(17.502*TX/(240.97+TX))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SP_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SPLIN2( X1A,X2A,YA,Y2A,X1,X2,Y,LM,LN,M,N )
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
!     Given an M by N tabulated function YA, and tabulated independent
!     variables X1A(M values) and X2A (N values), and the second 
!     derivatives in the array, Y2A of the one dimensional natural
!     cubic splines of the rows of YA; and given a desired interpolating
!     point X1, X2; this routine returns an interpolated function value
!     Y by bicubic spline interpolation.
!
!     Press, W.H., B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling.
!     1986.  Numerical Recipes, The Art of Scientific Computing,
!     Cambridge University Press, Cambridge.  pp. 101.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 14 May 2012.
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Parameter Statements----------------------------!
!
      PARAMETER(LSPLINE=256)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 X1A(M),X2A(N)
      REAL*8 YA(LM,LN),Y2A(LM,LN)
      REAL*8 YTMP(LSPLINE),Y2TMP(LSPLINE),YYTMP(LSPLINE)
!
!----------------------Executable Lines--------------------------------!
!
      DO 12 J = 1,M
        DO 11 K = 1,N
          YTMP(K) = YA(J,K)
          Y2TMP(K) = Y2A(J,K)
   11   CONTINUE
        CALL SPLINT( X2A,YTMP,Y2TMP,N,X2,YYTMP(J) )
   12 CONTINUE
      CALL SPLINE( X1A,YYTMP,M,Y2TMP )
      CALL SPLINT( X1A,YYTMP,Y2TMP,M,X1,Y )
!
!---  End of SPLIN2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SPLINE( X,Y,N,Y2 )
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
!     Cubic spline second derivative.
!
!     Press, W.H., B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling.
!     1986.  Numerical Recipes, The Art of Scientific Computing,
!     Cambridge University Press, Cambridge.  pp. 86-89.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 27 July 2002.
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Parameter Statements----------------------------!
!
      PARAMETER(LSPLINE=256)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 X(N),Y(N),Y2(N)
      REAL*8 U(LSPLINE)
!
!----------------------Executable Lines--------------------------------!
!
      IF( N.EQ.1 ) THEN
        Y2(1) = 0.D+0
        RETURN
      ENDIF
      Y2(1) = 0.D+0
      U(1) = 0.D+0
      DO 100 I = 2,N-1
        SIG = (X(I)-X(I-1))/(X(I+1)-X(I-1))
        P = SIG*Y2(I-1)+2.D+0
        Y2(I) = (SIG-1.D+0)/P
        U(I) = (6.D+0*((Y(I+1)-Y(I))/(X(I+1)-X(I)) -
     &    (Y(I)-Y(I-1))/(X(I)-X(I-1)))/
     &    (X(I+1)-X(I-1)) - SIG*U(I-1))/P
  100 CONTINUE
      QN = 0.D+0
      UN = 0.D+0
      Y2(N) = (UN-QN*U(N-1))/(QN*Y2(N-1)+1.D+0)
      DO 200 K = N-1,1,-1
        Y2(K) = Y2(K)*Y2(K+1)+U(K)
  200 CONTINUE
!
!---  End of SPLINE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SPLINT( XA,YA,Y2A,N,X,Y )
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
!     Cubic spline interpolation.
!
!     Press, W.H., B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling.
!     1986.  Numerical Recipes, The Art of Scientific Computing,
!     Cambridge University Press, Cambridge.  pp. 86-89.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 27 July 2002.
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XA(N),YA(N),Y2A(N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( N.EQ.1 ) THEN
        Y = YA(1)
        RETURN
      ELSEIF( X.LE.XA(1) ) THEN
        Y = YA(1)
        RETURN
      ELSEIF( X.GE.XA(N) ) THEN
        Y = YA(N)
        RETURN
      ENDIF
      KLO = 1
      KHI = N
   10 CONTINUE
      IF( KHI-KLO.GT.1 ) THEN
        K = (KHI+KLO)/2
        IF( XA(K).GT.X ) THEN
          KHI = K
        ELSE
          KLO = K
        ENDIF
        GOTO 10
      ENDIF
      H = XA(KHI)-XA(KLO)
      A = (XA(KHI)-X)/H
      B = (X-XA(KLO))/H
      Y = A*YA(KLO)+B*YA(KHI)+
     &    ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.D+0
!
!---  End of SPLINT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SR_W( TX,PX,ISRX )
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
!     Formulation subregion as a function of temperature and pressure.
!
!     Meyer, C.A., R.B. McClintock, G.J. Silvestri, and R.C. Spencer
!     1993.  ASME Steam Tables, The American Society of Mechanical
!     Engineers, New York.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 March 2002.
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
      REAL*8 LFCX(3),PSWX,PSRBX,THETAX,TOLX
      INTEGER ISRX
!
!----------------------Data Statements---------------------------------!
!
      DATA LFCX / 1.574373327D+1, -3.417061978D+1, 1.931380707D+1 /
      DATA TOLX / 1.D-2 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SR_W'
!
!---  Subcritical-critical temperature  ---
!
      TKX = TX+TABS
      IF( TKX.LE.TCRW ) THEN
!
!---    The K-function (saturation line)  ---
!
        CALL SP_W( TX,PSWX )
!
!---    Subregions 5/6 and 1/4 boundaries  ---
!
        IF( TX.LE.350.D+0 ) THEN
          IF( (PX-PSWX).GE.TOLX ) THEN
            ISRX = 1
          ELSEIF( (PX-PSWX).LE.-TOLX ) THEN
            ISRX = 2
          ELSE
            ISRX = 6
          ENDIF
        ELSE
          IF( (PX-PSWX).GE.TOLX ) THEN
            ISRX = 4
          ELSEIF( (PX-PSWX).LE.-TOLX ) THEN
            ISRX = 2
          ELSE
            ISRX = 5
          ENDIF
        ENDIF
!
!---  Supercritical temperature  ---
!
      ELSE
!
!---    The L-function (subregions 2/3 boundary)  ---
!
        THETAX = TKX/TCRW
        PSRBX = PCRW*(LFCX(1) + LFCX(2)*THETAX + LFCX(3)*(THETAX**2))
        ISRX = 2
        IF( PX.GT.PSRBX ) ISRX = 3
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SR_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THK_CO2( TX,PX,RHOAX,THKAX )
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
!     Calculation of CO2 thermal conductivity.
!
!     Vesovic, V., W.A. Wakeman, G.A. Olchowy, J.V. Sengers,
!     J.T.R. Watson, and J. Millat.  1990.  The Transport Properties
!     of Carbon Dioxide.  J. Phys. Chem. Ref. Data, 19(3):763-808.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 April 2002.
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
      PARAMETER( LTX=7,LPX=29 )
      REAL*8 SCX(5),SBX(8),SDX(4)
      REAL*8 TX_TA(LTX),PX_TA(LPX),THKX_TA(LTX,LPX),
     &  PSX_TA(LTX),THKLSX_TA(LTX),THKGSX_TA(LTX)
      REAL*8 TIX(2),PIX(2,2),THKIX(2,2)
!
!----------------------Data Statements---------------------------------!
!
      DATA SBX /  0.4226159D+0, 0.6280115D+0, -0.5387661D+0,
     &  0.6735941D+0, 0.0D+0, 0.0D+0, -0.4362677D+0, 0.2255388D+0  /
      DATA SCX / 2.387869D-2, 4.350794D+0, -10.33404D+0, 7.981590D+0,
     &  -1.940558D+0 /
      DATA SDX /  2.447164D-2, 8.705605D-5, -6.547950D-8,
     &  6.594919D-11  /
      DATA ESPX / 251.196D+0 /
      DATA RHOCRA / 467.69D+0 /
      DATA TX_TA / 298.D+0,   300.D+0,   302.D+0,   304.D+0,   306.D+0,
     &  308.D+0,   310.D+0 /
      DATA PX_TA / 0.1D+6,  0.5D+6,  1.D+6,   1.5D+6,  2.D+6,   2.5D+6,
     &  3.D+6,   3.5D+6,  4.D+6,   4.5D+6,  5.D+6,   5.5D+6,  6.D+6,
     &  6.5D+6,  7.D+6,   7.5D+6,  8.D+6,   8.5D+6,  9.D+6,   9.5D+6,
     &  10.D+6,   10.5D+6,  11.D+6,   11.5D+6,  12.D+6,   12.5D+6,
     &  13.D+6,   13.5D+6,  14.D+6 /
      DATA THKX_TA /
     &  16.61D+0,16.77D+0,16.93D+0,17.10D+0,17.26D+0,17.42D+0,17.59D+0,
     &  16.81D+0,16.97D+0,17.13D+0,17.29D+0,17.45D+0,17.61D+0,17.77D+0,
     &  17.10D+0,17.25D+0,17.41D+0,17.56D+0,17.72D+0,17.88D+0,18.04D+0,
     &  17.45D+0,17.60D+0,17.75D+0,17.90D+0,18.05D+0,18.20D+0,18.35D+0,
     &  17.88D+0,18.02D+0,18.16D+0,18.30D+0,18.44D+0,18.59D+0,18.73D+0,
     &  18.41D+0,18.54D+0,18.66D+0,18.79D+0,18.92D+0,19.05D+0,19.18D+0,
     &  19.07D+0,19.17D+0,19.27D+0,19.38D+0,19.49D+0,19.60D+0,19.72D+0,
     &  19.91D+0,19.97D+0,20.03D+0,20.11D+0,20.19D+0,20.27D+0,20.37D+0,
     &  20.97D+0,20.97D+0,20.99D+0,21.01D+0,21.05D+0,21.10D+0,21.15D+0,
     &  22.38D+0,22.28D+0,22.21D+0,22.16D+0,22.13D+0,22.11D+0,22.12D+0,
     &  24.31D+0,24.03D+0,23.81D+0,23.63D+0,23.50D+0,23.39D+0,23.31D+0,
     &  27.18D+0,26.52D+0,26.01D+0,25.61D+0,25.29D+0,25.03D+0,24.82D+0,
     &  32.26D+0,30.50D+0,29.29D+0,28.41D+0,27.74D+0,27.21D+0,26.79D+0,
     &  82.87D+0,39.22D+0,35.15D+0,32.87D+0,31.36D+0,30.27D+0,29.44D+0,
     &  83.94D+0,81.65D+0,59.63D+0,42.52D+0,37.70D+0,35.04D+0,33.28D+0,
     &  85.42D+0,82.70D+0,80.43D+0,81.57D+0,56.67D+0,44.52D+0,39.66D+0,
     &  86.92D+0,84.16D+0,81.52D+0,79.25D+0,78.64D+0,83.25D+0,53.66D+0,
     &  88.36D+0,85.65D+0,82.98D+0,80.39D+0,78.07D+0,76.57D+0,76.76D+0,
     &  89.74D+0,87.10D+0,84.46D+0,81.85D+0,79.30D+0,76.94D+0,74.95D+0,
     &  91.05D+0,88.48D+0,85.91D+0,83.33D+0,80.76D+0,78.26D+0,75.88D+0,
     &  92.29D+0,89.80D+0,87.29D+0,84.77D+0,82.24D+0,79.73D+0,77.27D+0,
     &  93.48D+0,91.05D+0,88.61D+0,86.16D+0,83.68D+0,81.20D+0,78.75D+0,
     &  94.62D+0,92.25D+0,89.87D+0,87.48D+0,85.06D+0,82.63D+0,80.22D+0,
     &  95.70D+0,93.39D+0,91.07D+0,88.73D+0,86.38D+0,84.01D+0,81.65D+0,
     &  96.75D+0,94.48D+0,92.22D+0,89.94D+0,87.64D+0,85.33D+0,83.02D+0,
     &  97.76D+0,95.54D+0,93.32D+0,91.09D+0,88.84D+0,86.59D+0,84.33D+0,
     &  98.73D+0,96.55D+0,94.38D+0,92.19D+0,89.99D+0,87.79D+0,85.58D+0,
     &  99.67D+0,97.53D+0,95.40D+0,93.26D+0,91.10D+0,88.94D+0,86.78D+0,
     &  100.58D+0,98.47D+0,96.38D+0,94.28D+0,92.17D+0,90.05D+0,87.94D+0
     &  /
      DATA PSX_TA / 6.4121D+6, 6.7131D+6, 7.0268D+6, 7.3555D+6,
     &  0.D+0, 0.D+0, 0.D+0 /
      DATA THKLSX_TA / 83.46D+0, 82.3D+0, 84.9D+0, 187.1D+0, 0.D+0,
     &  0.D+0, 0.D+0 /
      DATA THKGSX_TA / 45.34D+0, 53.13D+0, 47.61D+0, 39.02D+0, 0.D+0,
     &  0.D+0, 0.D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THK_CO2'
!
!---  Absolute and reduced temperature  ---
!
      TKX = TX+TABS
      TRX = TKX/ESPX
!
!---  Near-critical-point conditions  ---
!
      IF( TKX.GE.TX_TA(1) .AND. TKX.LE.TX_TA(LTX) .AND.
     &  PX.GE.PX_TA(1) .AND. PX.LE.PX_TA(LPX) ) THEN
!
!---    Temperature table index  ---
!
        ITLX = 1
        ITUX = LTX
   10   IF( ITUX-ITLX.GT.1 ) THEN
          ITM = (ITLX+ITUX)/2
          IF( TKX.GT.TX_TA(ITM) ) THEN
            ITLX = ITM
          ELSE
            ITUX = ITM
          ENDIF
          GOTO 10
        ENDIF
!
!---    Pressure table index  ---
!
        IPLX = 1
        IPUX = LPX
   20   IF( IPUX-IPLX.GT.1 ) THEN
          IPM = (IPLX+IPUX)/2
          IF( PX.GT.PX_TA(IPM) ) THEN
            IPLX = IPM
          ELSE
            IPUX = IPM
          ENDIF
          GOTO 20
        ENDIF
!
!---    Check for saturation boundary  ---
!
        TIX(1) = TX_TA(ITLX)
        TIX(2) = TX_TA(ITUX)
        PIX(1,1) = PX_TA(IPLX)
        PIX(2,1) = PX_TA(IPLX)
        PIX(1,2) = PX_TA(IPUX)
        PIX(2,2) = PX_TA(IPUX)
        THKIX(1,1) = THKX_TA(ITLX,IPLX)
        THKIX(2,1) = THKX_TA(ITUX,IPLX)
        THKIX(1,2) = THKX_TA(ITLX,IPUX)
        THKIX(2,2) = THKX_TA(ITUX,IPUX)
        IF( (PSX_TA(ITLX).GT.PX_TA(IPLX)) .AND.
     &    (PSX_TA(ITLX).LT.PX_TA(IPUX)) ) THEN
          IF( RHOAX.LE.RHOCRA ) THEN
            PIX(1,2) = PSX_TA(ITLX)
            THKIX(1,2) = THKGSX_TA(ITLX)
          ELSE
            PIX(1,1) = PSX_TA(ITLX)
            THKIX(1,1) = THKLSX_TA(ITLX)
          ENDIF
        ENDIF
        IF( (PSX_TA(ITUX).GT.PX_TA(IPLX)) .AND.
     &    (PSX_TA(ITUX).LT.PX_TA(IPUX)) ) THEN
          IF( RHOAX.LE.RHOCRA ) THEN
            PIX(2,2) = PSX_TA(ITUX)
            THKIX(2,2) = THKGSX_TA(ITUX)
          ELSE
            PIX(2,1) = PSX_TA(ITUX)
            THKIX(2,1) = THKLSX_TA(ITUX)
          ENDIF
        ENDIF
!
!---    Bilinear interpolation  ---
!
        VTLX = (THKIX(1,2)-THKIX(1,1))*(PX-PIX(1,1))/(PIX(1,2)-PIX(1,1))
     &    + THKIX(1,1)
        VTUX = (THKIX(2,2)-THKIX(2,1))*(PX-PIX(2,1))/(PIX(2,2)-PIX(2,1))
     &    + THKIX(2,1)
        THKAX = ((VTUX-VTLX)*(TKX-TIX(1))/(TIX(2)-TIX(1))) + VTLX
!
!---  Non-near-critical-point conditions  ---
!
      ELSE
!
!---    Zero-density-limit component  ---
!
        COKX = 0.D+0
        DO 40 I = 1,5
          COKX = COKX + SCX(I)*((TKX/1.D+2)**(2-I))
   40   CONTINUE
        COKX = 1.D+0 + EXP(-183.5D+0/TKX)*COKX
        SRX = SQRT( 2.D+0*COKX/5.D+0 )
        ZETAX = 0.D+0
        DO 50 I = 0,7
          ZETAX = ZETAX + SBX(I+1)/(TRX**I)
   50   CONTINUE
        THK_ZD = 0.475598D+0*SQRT(TKX)*(1.D+0+(SRX**2))/ZETAX
!
!---    Density-dependent component  ---
!
        THK_DD = 0.D+0
        DO 60 I = 1,4
          THK_DD = THK_DD + SDX(I)*(RHOAX**I)
   60   CONTINUE
!
!---    Sum components and convert to W/m K  ---
!
        THKAX = (THK_ZD + THK_DD)*1.D-3
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THK_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THK_B( TX,XLSX,THKLWX,THKBX )
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
!     This subroutine calculates the thermal conductivity of pure
!     water as a function of temperature.
!
!     tx - temperature, C
!     xlsx - mass fraction of aqueous NaCl
!     tkbx - thermal conductivity of NaCl brine, W/m K
!
!     Ozbek, H. and S.L. Phillips.  1980.  "Thermal conductivity of
!     aqueous sodium chloride solutions from 20 to 330 C."
!     J. Chem. Engr. Data, 25:263-267.
!
!     Temperature Range:  20 - 330 C
!     NaCl Concentration Range:  0 - saturation
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 4 April 2002.
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
      REAL*8 SCX(6)
!
!----------------------Data Statements---------------------------------!
!
      DATA SCX / 2.3434D-3, -7.924D-6, 3.924D-8, 1.06D-5, -2.D-8,
     &  1.2D-10 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THK_B'
!
!---  Yusufova correlation  ---
!
      YLSX = 1.D+2*XLSX
      THKBX = 1.D+0 - (SCX(1) + SCX(2)*TX + SCX(3)*(TX**2))*YLSX +
     &  (SCX(4) + SCX(5)*TX + SCX(6)*(TX**2))*(YLSX**2)
      THKBX = THKBX*THKLWX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THK_B group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THK_G( PHIX,THKAX,THKOX,THKNX,THKWX,
     &  XMAX,XMOX,XMNX,XMWX,THKGX )
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
!     Calculate the gas thermal conductivity using the modification of
!     Mason and Saxena to the Waasiljewa equation.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York. pp: 530-531.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 April 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYDT
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 PHIX(LHF_HT+1,LHF_HT+1),THKX(LHF_HT+1),XMGX(LHF_HT+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THK_G'
      THKX(1) = THKAX
      THKX(2) = THKOX
      THKX(3) = THKNX
      THKX(4) = THKWX
      XMGX(1) = XMAX
      XMGX(2) = XMOX
      XMGX(3) = XMNX
      XMGX(4) = XMWX
!
!---  Loop over component pairs  ---
!
      THKGX = 0.D+0
      DO 40 I = 1,NHF_HT+1
        CHIX = 0.D+0
        DO 30 J = 1,NHF_HT+1
          CHIX = CHIX + XMGX(J)*PHIX(I,J)
   30   CONTINUE
        THKGX = THKGX + XMGX(I)*THKX(I)/(CHIX+SMALL)
   40 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THK_G group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THK_H( TX,THKHX )
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
!     Calculate the thermal conductivity of methane hydrate.
!
!     American Geophysical Union Conference in San Francisco, 
!     December 8-12, 2003
!  
!     http://woodshole.er.usgs.gov/operations/hi_fi/pure_hydrate
!       /index.html
!
!     Waite, W.F., L.A. Stern, S.H. Kirby, W.J., Winters, and D.H. Mason
!     2007. Simultaneous determination of thermal conductivity,
!     thermal diffusivity, and specific heat in sI methane hydrate,
!     Geophysical Journal International, v. 169, p. 767–774, 
!     doi:10.1111/j.1365-246X.2007.03382.x.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 November 2004.
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THK_H'
!
!---  Weak temperature dependence  ---
!
      THKHX = 5.D-1 - 3.75D-4*TX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THK_H group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THK_I( TX,THKIX )
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
!     Calculate the thermal conductivity of ice.
!
!     Fourth International Conference on Gas Hydrates as part of a 
!     paper on thermal conductivity in methane hydrate: 
!     Waite, W., Pinkston, J., Kirby, S.H., 2002, Preliminary 
!     laboratory thermal conductivity measurements in pure methane 
!     hydrate and hydrate-sediment mixtures; A progress report: 
!     Proceedings of the Fourth International Conference on Gas 
!     Hydrates, Yokohama, Japan, May 19-23, p. 728-733.
!
!     http://woodshole.er.usgs.gov/operations/hi_fi/ice_results
!       /index.html
!
!     Calculate the thermal conductivity of ice as a function of
!     temperature using polynomial fit of data from Dickerson (1969).
!
!     Dickerson, R. W., Jr., Thermal properties of food, in The Freezing
!     Preservation of Foods, 4th ed., Vol. 2, D. K. Tressler, W. B. Van
!     Arnsdel, and M. J. Copley (Editors), AVI Publishing Co., Westport,
!     Conn.
!
!     The temperature is limited in this subroutine to the following
!     values:  -100.0 C < T <= 0.0 C
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 November 2004.
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
      REAL*8 A(3)
!
!----------------------Data Statements---------------------------------!
!
      DATA A / 7.39519D+0,-2.86936D-2,3.54452D-5 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THK_I'
!      THKIX = 2.13D+0 - 1.32D-2*TX
      TKX = TX+TABS
      THKIX = A(1) +A(2)*TKX +A(3)*(TKX**2)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THK_I group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THK_N( THKNAX,THKNOX,THKNNX,XMNAX,XMNOX,XMNNX,THKNX )
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
!     Calculate the nonaqueous liquid thermal conductivity with the 
!     method of Li, using critical molar volumes for the Vi terms.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York. pp: 565.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 24 September 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYDT
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
      REAL*8 XMNX(LHF_HT),OMGX(LHF_HT),THKX(LHF_HT)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THK_N'
      XMNX(1) = XMNAX
      XMNX(2) = XMNOX
      XMNX(3) = XMNNX
      THKX(1) = THKNAX
      THKX(2) = THKNOX
      THKX(3) = THKNNX
      SOMGX = 0.D+0
      DO 10 I = 1,NHF_HT
        SOMGX = SOMGX + XMNX(I)*HCPP(5,I)
   10 CONTINUE
      DO 20 I = 1,NHF_HT
        OMGX(I) = XMNX(I)*HCPP(5,I)/(SOMGX+SMALL)
   20 CONTINUE
      THKNX = 0.D+0
      DO 40 I = 1,NHF_HT
        DO 30 J = 1,NHF_HT
          THKNX = THKNX + 2.D+0*OMGX(I)*OMGX(J)/
     &      ((1.D+0/THKX(I)) + (1.D+0/THKX(J)))
   30   CONTINUE
   40 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THK_N group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THK_N2( TX,RHONX,THKNX )
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
!     Calculation of N2 thermal conductivity.
!
!     Lemmon, E.W., and R.T. Jacobsen. 2004. "Viscosity and thermal
!     conductivity equations for nitrogen, oxygen, argon, and air,"
!     International Journal of Thermophysics 25(1):21-69.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 20 September 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
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
      REAL*8 BX(5)
      REAL*8 CGX(9),CNX(9),CTX(9)
      INTEGER IDX(9),ILX(9)
!
!----------------------Data Statements---------------------------------!
!
      DATA BX / 0.431D+0,-0.4623D+0,0.08406D+0,0.005341D+0,-0.00331D+0 /
      DATA CGX / 0.D+0,0.D+0,0.D+0,0.D+0,0.D+0,1.D+0,1.D+0,1.D+0,1.D+0 /
      DATA CNX / 1.511D+0, 2.117D+0, -3.332D+0, 8.862D+0, 31.11D+0,
     &  -73.13D+0, 20.03D+0, -0.7096D+0, 0.2672D+0 /
      DATA CTX / 0.D+0, -1.D+0, -0.7D+0, 0.D+0, 0.03D+0, 0.2D+0,
     &  0.8D+0, 0.6D+0, 1.9D+0 /
      DATA IDX / 0,0,0,1,2,3,4,8,10 /
      DATA ILX / 0,0,0,0,0,1,2,2,2 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THK_N2'
!
!---  Absolute temperature, K  ---
!
      TKX = TX + TABS
!
!---  Lennard-Jones energy parameter, K  ---
!
      EPLJ = 98.94D+0
!
!---  Collision integral  ---
!
      TRX = TKX/EPLJ
      OMEGA = 0.D+0
      DO 10 M = 1,5
        OMEGA = OMEGA + BX(M)*((LOG(TRX))**(M-1))
   10 CONTINUE
      OMEGA = EXP( OMEGA )
!
!---  Lennard-Jones size parameter, nm  ---
!
      SIGMA = 0.3656D+0
!
!---  Dilute gas viscosity, uPa s  ---
!
      VISNOX = 0.0266958D+0*SQRT(WTMN*TKX)/((SIGMA**2)*OMEGA)
!
!---  HCPP(3,3) critical temperature, K
!     HCPP(5,3) critical molar volume, m^3/kmol  ---
!
      TAU = HCPP(3,3)/TKX
      RHOMNX = RHONX/WTMN
      DELTA = RHOMNX*HCPP(5,3)
!
!---  Dilute gas thermal conductivity, mW/m K  ---
!
      THKNOX = CNX(1)*VISNOX + CNX(2)*(TAU**CTX(2)) 
     &  + CNX(3)*(TAU**CTX(3))
!
!---  Residual gas thermal conductivity, mW/m K  ---
!
      THKNRX = 0.D+0
      DO 20 M = 4,9
        THKNRX = THKNRX + CNX(M)*(TAU**CTX(M))*(DELTA**IDX(M))*
     &    EXP(-CGX(M)*(DELTA**ILX(M))) 
   20 CONTINUE
!
!---  Critical enhancement set to zero as minimum temperaure of 
!     223.15 K is well above N2 critical temperature of 126.192K,
!     mW/m K  ---
!
      THKNCX = 0.D+0
!
!---  Sum of components, W/m K  ---
!
      THKNX = 1.D-3*(THKNOX + THKNRX + THKNCX)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THK_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THK_CH4( TX,RHOOX,THKOX )
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
!     CH4 thermal conductivity.
!
!     Hanley, H.J.M., W.M. Haynes, and R.D. McCarty. 1997. "The
!     viscosity and thermal conductivity coefficients for dense
!     gaseous and liquid methane." J. Phys. Chem. Ref. Data,
!     6(2):597-609.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 24 September 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NAPL
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
      REAL*8 DX(4),FX(7),GTX(9)
!
!----------------------Data Statements---------------------------------!
!
      DATA DX / -0.25276292D+0, 0.33432859D+0, 1.12D+0, 168.D+0 /
      DATA FX / -7.0403639907D+0, 12.319512908D+0, -8.8525979933D+2,
     &  72.835897919D+0, 0.74421462902D+0, -2.9706914540D+0,
     &  2.2209758501D+3 /
      DATA GTX / -2.147621D+5, 2.190461D+5, -8.618097D+4, 1.496099D+4,
     &  -4.730660D+2, -2.331178D+2, 3.778439D+1, -2.320481D+0,
     &  5.311764D-2 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THK_CH4'
!
!---  Absolute temperature, K ---
!
      TKX = TX + TABS
!
!---  Critical density and density, gm/cm^3 ---
!
      RHOCRX = 1.D-3*WTMO/HCPP(5,2)
      RHOX = 1.D-3*RHOOX
!
!---  Zero-density-limit component, mW/m K  ---
!
      THKO0X = GTX(1)*(TKX**(-1.D+0)) + GTX(2)*(TKX**(-2.D+0/3.D+0))
     &   + GTX(3)*(TKX**(-1.D+0/3.D+0)) + GTX(4)
     &   + GTX(5)*(TKX**(1.D+0/3.D+0)) + GTX(6)*(TKX**(2.D+0/3.D+0))
     &   + GTX(7)*TKX + GTX(8)*(TKX**(4.D+0/3.D+0))
     &   + GTX(9)*(TKX**(5.D+0/3.D+0))
!
!---  First density correction, mW/m K  ---
!
      THKO1X = RHOX*(DX(1) + DX(2)*((DX(3) - LOG(TKX/DX(4)))**2))
!
!---  Remainder density correction, mW/m K  ---
!
      THETAX = (RHOX - RHOCRX)/RHOCRX
      THKO2X = EXP( FX(1) + FX(4)/TKX )*EXP( (RHOX**1.D-1)*(FX(2) + 
     &  FX(3)/(TKX**1.5D+0)) + THETAX*SQRT(RHOX)*(FX(5) + (FX(6)/TKX)
     &  + FX(7)/(TKX**2)) ) - 1.D+0
!
!---  Sum of components, W/m K  ---
!
      THKOX = 1.D-3*(THKO0X + THKO1X + THKO2X)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THK_CH4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THK_W( TX,PX,RHOX,THKX )
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
!     Thermal conductivity (W/m K) of pure water as a function of
!     temperature and density.
!
!     Meyer, C.A., R.B. McClintock, G.J. Silvestri, and R.C. Spencer
!     1993.  ASME Steam Tables, The American Society of Mechanical
!     Engineers, New York.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 27 March 2002.
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
      REAL*8 SAX(4),SBX(3),CBX(2),SDX(4),CCX(6)
!
!----------------------Data Statements---------------------------------!
!
      DATA TREF / 6.4727D+2 /
      DATA RHOREF / 3.17763D+2 /
      DATA PREF / 2.2115D+7 /
      DATA THKREF / 1.D+0 /
      DATA SAX / 0.0102811D+0, 0.0299621D+0, 0.0156146D+0,
     &  -0.00422464 /
      DATA SBX / -0.397070D+0, 0.400302D+0, 1.060000D+0 /
      DATA CBX / -0.171587D+0, 2.392190D+0 /
      DATA SDX / 0.0701309D+0, 0.0118520D+0, 0.00169937D+0, -1.0200D+0 /
      DATA CCX / 0.642857D+0, -4.11717D+0, -6.17937D+0, 0.00308976D+0,
     &  0.0822994D+0, 10.0932D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THK_W'
!
!---  Restrict temperature to above 0.01 C  ---
!
      TY = MAX( TX,1.D-2 )
!
!---  Zero density  ---
!
      IF( RHOX.LT.EPSL ) THEN
        THKX = 0.D+0
!
!---    Reset subroutine string sequence  ---
!
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Reduced temperature, density, and pressure  ---
!
      THETAX = (TY+TABS)/TREF
      RHOBX = RHOX/RHOREF
      BETAX = PX/PREF
!
!---  Zero term  ---
!
      THKX = 0.D+0
      ISAX = 0
      DO 10 I = 0,3
        ISAX = ISAX + 1
        THKX = THKX + SAX(ISAX)*(THETAX**I)
   10 CONTINUE
      THKX = SQRT(THETAX)*THKX
!
!---  First term  ---
!
      THKX = THKX + SBX(1) + SBX(2)*RHOBX +
     &  SBX(3)*EXP(CBX(1)*((RHOBX+CBX(2))**2))
!
!---  Second term  ---
!
      DTHETAX = ABS(THETAX-1.D+0) + CCX(4)
      CQX = 2.D+0 + CCX(5)/(DTHETAX**(3.D+0/5.D+0))
      IF( DTHETAX.GE.1.D+0 ) THEN
        CSX = 1.D+0/DTHETAX
      ELSE
        CSX = CCX(6)/(DTHETAX**(3.D+0/5.D+0))
      ENDIF
      THKX = THKX + ((SDX(1)/(THETAX**10))+SDX(2))*
     &  (RHOBX**(9.D+0/5.D+0))*
     &  EXP(CCX(1)*(1.D+0-(RHOBX**(14.D+0/5.D+0))))
      THKX = THKX + SDX(3)*CSX*(RHOBX**CQX)*
     &  EXP((CQX/(1.D+0+CQX))*(1.D+0-(RHOBX**(1.D+0+CQX))))
      VARX = CCX(2)*(THETAX**(3.D+0/2.D+0)) + CCX(3)/(RHOBX**5)
      VARX = MAX( VARX,-50.D+0 )
      TKX = TKX + SDX(4)*EXP(VARX)
!
!---  Dimensionalize thermal conductivity  ---
!
      THKX = THKX*THKREF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THK_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VISC_CO2( TX,RHOAX,VISAX )
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
!     Calculation of CO2 viscosity using the formulation of
!     Fenghour et al. in the temperature range 200K < T < 1500K
!     and densities up to 1400 kg/m^3.
!
!     Fenghour, A., W. A. Wakeham, V. Vesovic.  1998.  "The viscosity
!     of carbon dioxide."  J. Phys. Chem. Ref. Data, 27(1):31-41.
!
!     Calculation of CH4 vapor viscosity with the corresponding 
!     states method.  pp 397.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 April 2002.
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
      REAL*8 SAX(5),SBX(5)
!
!----------------------Data Statements---------------------------------!
!
      DATA SAX / 0.235156D+0, -0.491266D+0,
     &  5.211155D-2, 5.347906D-2, -1.537102D-2 /
      DATA SBX / 0.4071119D-2, 0.7198037D-4,
     &  0.2411697D-16, 0.2971072D-22, -0.1627888D-22 /
      DATA ESPX / 251.196D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VISC_CO2'
!
!---  Zero-density-limit component  ---
!
      TKX = TX+TABS
      TRX = TKX/ESPX
      ECSX = 0.D+0
      DO I = 0,4
        ECSX = ECSX + SAX(I+1)*(LOG(TRX)**I)
      END DO
      RECSX = EXP(ECSX)
      VIS_ZD = 1.00697D+0*SQRT(TKX)/RECSX
!
!---  Excess-viscosity component  ---
!
      VIS_EX = SBX(1)*RHOAX + SBX(2)*(RHOAX**2) +
     &  SBX(3)*(RHOAX**6)/(TRX**3) +
     &  SBX(4)*(RHOAX**8) + SBX(5)*(RHOAX**8)/TRX
!
!---  Sum components and convert to Pa s  ---
!
      VISAX = (VIS_ZD + VIS_EX)*1.D-6
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VISC_CO2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VISC_B( TX,XLSX,VISWX,VISBX )
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
!     Viscosity (Pa s) of brine as a function of temperature (C),
!     NaCl brine mass fraction and pure water viscosity (Pa s)
!
!     Phillips, S.L., A. Igbene, J.A. Fair, H. Ozbek, and M. Tavana.
!     1981.  A Technical Databook for Geothermal Energy Utilization
!     LBL-12810, UC-66a, Lawrence Berkeley Laboratory, University of
!     California, Berkeley, California.
!
!     Temperature Range 10-350 C
!     Pressure Range 0.1-50 MPa
!     NaCl Concentration Range 0-5 Molal (mol NaCl/kg H2O)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 4 April 2002.
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
      REAL*8 SAX(5)
!
!----------------------Data Statements---------------------------------!
!
      DATA SAX / 0.0816D+0, 0.0122D+0, 0.000128D+0, 0.000629D+0,
     &  -0.7D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VISC_B'
!
!---  Convert mass fraction to molality  ---
!
      GLSX = 1.D+3*XLSX/(WTMS*(1.D+0-XLSX))
!
!---  Formulation of Phillips et al.  ---
!
      VISBX = VISWX*(1.D+0 + SAX(1)*GLSX + SAX(2)*(GLSX**2) +
     &  SAX(3)*(GLSX**3) + SAX(4)*TX*(1.D+0-EXP(SAX(5)*GLSX)))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VISC_B group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VISC_G( PHIX,VISAX,VISOX,VISNX,VISWX,
     &  XMAX,XMOX,XMNX,XMWX,VISGX )
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
!     Calculate the gas viscosity using the method of Wilke, p.407.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York. pp: 332-337.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 April 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYDT
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
      REAL*8 PHIX(LHF_HT+1,LHF_HT+1),VISX(LHF_HT+1),XMGX(LHF_HT+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VISC_G'
      VISX(1) = VISAX
      VISX(2) = VISOX
      VISX(3) = VISNX
      VISX(4) = VISWX
      XMGX(1) = XMAX
      XMGX(2) = XMOX
      XMGX(3) = XMNX
      XMGX(4) = XMWX
!
!---  Loop over component pairs  ---
!
      DO 20 J = 1,NHF_HT+1
        DO 10 I = 1,NHF_HT+1
          PHIX(I,J) = ((1.D+0+SQRT(VISX(I)/(VISX(J)+SMALL))*
     &      SQRT(SQRT(HCPP(1,J)/HCPP(1,I))))**2)/
     &      SQRT(8.D+0*(1.D+0 + HCPP(1,I)/HCPP(1,J)))
   10   CONTINUE
   20 CONTINUE
      VISGX = 0.D+0
      DO 40 I = 1,NHF_HT+1
        CHIX = 0.D+0
        DO 30 J = 1,NHF_HT+1
          CHIX = CHIX + XMGX(J)*PHIX(I,J)
   30   CONTINUE
        VISGX = VISGX + XMGX(I)*VISX(I)/(CHIX+SMALL)
   40 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VISC_G group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VISC_L( XMLAX,XMLOX,XMLNX,
     &  VISBX,VISAX,VISOX,VISNX,VISLX )
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
!     Viscosity (Pa s) of aqueous NaCl solutions with dissolved CO2,
!     following the method of Grunberg and Nissan, calibrated against
!     the data of Kumagai and Yokoyama.
!
!     Kumagai, A., and C. Yokoyama.  1999.  "Viscosities of aqueous
!     NaCl solutions containing CO2 at high pressures."  J. Chem. Eng.
!     Data, 44:227-229.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids. pp. 474-475.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 30 April 2002.
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VISC_L'
!
!---  Formulation of Grunberg and Nissan  ---
!
      VISLX = (1.D+0-XMLAX-XMLOX-XMLNX)*LOG(VISBX)
      IF( VISAX.GT.EPSL ) VISLX = VISLX + XMLAX*LOG(VISAX)
      IF( VISOX.GT.EPSL ) VISLX = VISLX + XMLOX*LOG(VISOX)
      IF( VISNX.GT.EPSL ) VISLX = VISLX + XMLNX*LOG(VISNX)
      VISLX = EXP( VISLX )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VISC_L group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VISC_N2( TX,RHONX,VISNX )
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
!     Calculation of N2 viscosity.
!
!     Lemmon, E.W., and R.T. Jacobsen. 2004. "Viscosity and thermal
!     conductivity equations for nitrogen, oxygen, argon, and air,"
!     International Journal of Thermophysics 25(1):21-69.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 20 September 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
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
      REAL*8 BX(5)
      REAL*8 CGX(5),CNX(5),CTX(5)
      INTEGER IDX(5),ILX(5)
!
!----------------------Data Statements---------------------------------!
!
      DATA BX / 0.431D+0,-0.4623D+0,0.08406D+0,0.005341D+0,-0.00331D+0 /
      DATA CGX / 0.D+0,1.D+0,1.D+0,1.D+0,1.D+0 /
      DATA CNX / 10.72D+0,0.03989D+0,0.001208D+0,-7.402D+0,4.620D+0 /
      DATA CTX / 0.1D+0,0.25D+0,3.2D+0,0.9D+0,0.3D+0 /
      DATA IDX / 2,10,12,2,1 /
      DATA ILX / 0,1,1,2,3 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VISC_N2'
!
!---  Absolute temperature, K  ---
!
      TKX = TX + TABS
!
!---  Lennard-Jones energy parameter, K  ---
!
      EPLJ = 98.94D+0
!
!---  Collision integral  ---
!
      TRX = TKX/EPLJ
      OMEGA = 0.D+0
      DO 10 M = 1,5
        OMEGA = OMEGA + BX(M)*((LOG(TRX))**(M-1))
   10 CONTINUE
      OMEGA = EXP( OMEGA )
!
!---  Lennard-Jones size parameter, nm  ---
!
      SIGMA = 0.3656D+0
!
!---  Dilute gas viscosity, Pa s  ---
!
      VISNOX = 0.0266958D-6*SQRT(WTMN*TKX)/((SIGMA**2)*OMEGA)
!
!---  Residual fluid viscosity, Pa s  
!     HCPP(3,3) critical temperature, K
!     HCPP(5,3) critical molar volume, m^3/kmol ---
!
      TAU = HCPP(3,3)/TKX
      RHOMNX = RHONX/WTMN
      DELTA = RHOMNX*HCPP(5,3)
      VISNRX = 0.D+0
      DO 20 M = 1,5
        VISNRX = VISNRX + CNX(M)*(TAU**CTX(M))*(DELTA**IDX(M))*
     &    EXP(-CGX(M)*(DELTA**ILX(M)))
   20 CONTINUE
      VISNRX = 1.D-6*VISNRX
!
!---  Sum of components  ---
!
      VISNX = VISNOX + VISNRX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VISC_N2 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VISC_CH4( TX,RHOOX,VISOX )
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
!     Hanley, H.J.M., W.M. Haynes, and R.D. McCarty. 1997. "The
!     viscosity and thermal conductivity coefficients for dense
!     gaseous and liquid methane." J. Phys. Chem. Ref. Data,
!     6(2):597-609.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 April 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NAPL
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
      REAL*8 DX(4),FX(7),GVX(9)
!
!----------------------Data Statements---------------------------------!
!
      DATA DX / 1.696985927D+0, -0.133372346D+0, 1.4D+0, 168.D+0 /
      DATA FX / -1.035060586D+1, 1.7571599671D+1, -3.0193918656D+3,
     &  1.8873011594D+2, 4.2903609488D-2, 1.4529023444D+2,
     &  6.1276818706D+3 /
      DATA GVX / -2.090975D+5, 2.647269D+5, -1.472818D+5, 4.716740D+4,
     &  -9.491872D+3, 1.219979D+3, -9.627993D+1, 4.274152D+0,
     &  -8.141531D-2 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VISC_CH4'
!
!---  Absolute temperature, K  ---
!
      TKX = TX + TABS
!
!---  Critical density and density, gm/cm^3 ---
!
      RHOCRX = 1.D-3*WTMO/HCPP(5,2)
      RHOX = 1.D-3*RHOOX
!
!---  Zero-density-limit component, ug/(cm s)  ---
!
      VISO0X = GVX(1)*(TKX**(-1.D+0)) + GVX(2)*(TKX**(-2.D+0/3.D+0))
     &   + GVX(3)*(TKX**(-1.D+0/3.D+0)) + GVX(4)
     &   + GVX(5)*(TKX**(1.D+0/3.D+0)) + GVX(6)*(TKX**(2.D+0/3.D+0))
     &   + GVX(7)*TKX + GVX(8)*(TKX**(4.D+0/3.D+0))
     &   + GVX(9)*(TKX**(5.D+0/3.D+0))
!
!---  First density correction, ug/(cm s)  ---
!
      VISO1X = RHOX*(DX(1) + DX(2)*((DX(3) - LOG(TKX/DX(4)))**2))
!
!---  Remainder density correction, ug/(cm s)  ---
!
      THETAX = (RHOX - RHOCRX)/RHOCRX
      VISO2X = EXP( FX(1) + FX(4)/TKX )*EXP( (RHOX**1.D-1)*(FX(2) + 
     &  FX(3)/(TKX**1.5D+0)) + THETAX*SQRT(RHOX)*(FX(5) + (FX(6)/TKX)
     &  + FX(7)/(TKX**2)) ) - 1.D+0
!
!---  Sum of components, Pa s  ---
!
      VISOX = 1.D-7*(VISO0X + VISO1X + VISO2X)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VISC_CH4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VISC_W( TX,PX,RHOWX,VISWX )
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
!     Viscosity (Pa s) of pure water as a function of temperature and
!     density.
!
!     Meyer, C.A., R.B. McClintock, G.J. Silvestri, and R.C. Spencer
!     1993.  ASME Steam Tables, The American Society of Mechanical
!     Engineers, New York.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 27 March 2002.
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
      REAL*8 CHX(46)
!
!----------------------Data Statements---------------------------------!
!
      DATA TREF / 6.4727D+2 /
      DATA RHOREF / 3.17763D+2 /
      DATA PREF / 2.2115D+7 /
      DATA VISREF / 5.5071D+1 /
      DATA CHX / 1.D+0, 9.78197D-1, 5.79829D-1, -2.02354D-1,
     &  5.132047D-1, 3.205656D-1, 0.D+0, 0.D+0, -7.782567D-1,
     &  1.885447D-1, 2.151778D-1, 7.317883D-1, 1.241044D+0,
     &  1.476783D+0, 0.D+0, 0.D+0, -2.818107D-1, -1.070786D+0,
     & -1.263184D+0, 0.D+0, 0.D+0, 0.D+0, 1.778064D-1,
     &  4.605040D-1, 2.340379D-1, -4.924179D-1, 0.D+0, 0.D+0,
     &  -4.176610D-2, 0.D+0, 0.D+0, 1.600435D-1, 0.D+0, 0.D+0,
     &  0.D+0, -1.578386D-2, 0.D+0, 0.D+0, 0.D+0, 0.D+0, 0.D+0,
     &  0.D+0, 0.D+0, -3.629481D-3,  0.D+0, 0.D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VISC_W'
!
!---  Restrict temperature to above 0.01 C  ---
!
      TY = MAX( TX,1.D-2 )
!
!---  Reduced temperature, density, and pressure  ---
!
      THETAX = (TY+TABS)/TREF
      RHOBX = RHOWX/RHOREF
      BETAX = PX/PREF
!
!---  Zero term  ---
!
      VISWX = 0.D+0
      ICH = 0
      DO 10 I = 0,3
        ICH = ICH + 1
        VISWX = VISWX + CHX(ICH)/(THETAX**I)
   10 CONTINUE
      VISWX = SQRT(THETAX)/VISWX
!
!---  First term  ---
!
      VISAX = 0.D+0
      DO 30 I = 0,5
        DO 20 J = 0,6
          ICH = (J*6) + I + 5
          VISAX = VISAX + CHX(ICH)*(((1.D+0/THETAX)-1.D+0)**I)*
     &      ((RHOBX-1.D+0)**J)
   20   CONTINUE
   30 CONTINUE
      VISWX = VISWX*EXP(RHOBX*VISAX)
!
!---  Second term  ---
!
      IF( THETAX.GE.0.9970 .AND. THETAX.LE.1.0082 .AND.
     &  RHOBX.GE.0.755 .AND. RHOBX.LE.1.290 ) THEN
        DPX = 1.D-1
        DBETAX = DPX/PREF
        PIX = PX+DPX
        CALL SR_W( TY,PIX,ISRX )
        CALL DENS_W( TY,PIX,RHOLX,RHOVX,ISRX )
        IF( (1.D+0-ABS(RHOLX/RHOWX)).LT.(1.D+0-ABS(RHOVX/RHOWX)) ) THEN
          RHOBIX = RHOLX/RHOREF
        ELSE
          RHOBIX = RHOVX/RHOREF
        ENDIF
        CHIX = RHOBX*(RHOBIX-RHOBX)/DBETAX
        IF( CHIX.GE.21.93 ) VISWX = VISWX*0.922D+0*(CHIX**0.0263D+0)
      ENDIF
!
!---  Dimensionalize viscosity  ---
!
      VISWX = 1.D-6*VISWX*VISREF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VISC_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VOCDFG( XMW,XMO,XMA,DFGAO,DFGOW,DFGOX )
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
!     Calculates VOC diffusion coefficient by the Wilke method. pp. 34.
!
!     Falta, R.W., K. Pruess, I. Javandel, and P.A. Witherspoon. 1990.
!     Numerical Modeling of Steam Injection for the Removal of
!     Nonaqueous Phase Liquids from the Subsurface: 1 Numerical
!     Formulation. LBL-29615, Lawrence Berkeley Laboratory.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VOCDFG'
      IF( (1.D+0-XMO) .LE. 1.D-20 ) THEN
        DFGOX = 2.D+0*DFGAO*DFGOW/(DFGAO+DFGOW)
      ELSE
        DFGOX = (1.D+0-XMO)/((XMA/DFGAO) + (XMW/DFGOW))
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VOCDFG group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VOLM_LW( TX,PX,VMLWX )
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
!     Molar volume of water in aqueous.
!
!     Klauda, J. B., and S. I. Sandler.  2000.  "A fugacity model
!     for gas hydrate phase equilibria.  Ind. Eng. Chem. Res.
!     39:3377-3386.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 16 November 2004.
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VOLM_LW'
!
!---  Molar volume, m^3/kmol  ---
!
      PMX = (PX-PATM)*1.D-6
      VMLWX = 1.D+3*EXP( -10.9241D+0 + 2.5D-4*TX - 3.532D-4*PMX
     &    + 1.559D-7*(PMX**2))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VOLM_LW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VPL( PSBX,PVBX,SLX,TX,IZN )
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
!     Vapor pressure lowering versus aqueous saturation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 1 August 2013.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
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
      SUB_LOG(ISUB_LOG) = '/VPL'
!
!---  Absolute temperature  ---
!
      TKX = TX + TABS
!
!---  van Genuchten saturation function
!
      IF( ISCHR(IZN).EQ.1 .OR. ISCHR(IZN).EQ.11 .OR.
     &  ISCHR(IZN).EQ.13 ) THEN
        CAX = SCHR(1,IZN)
        CNX = MAX( SCHR(3,IZN),SMALL )
        IF( SCHR(14,IZN).LE.ZERO ) THEN
          IF( IRPL(IZN).EQ.2 ) THEN
            CMX = 1.D+0 - 2.D+0/CNX
          ELSE
            CMX = 1.D+0 - 1.D+0/CNX
         ENDIF
        ELSE
          CMX = SCHR(14,IZN)
        ENDIF
        SMPX = SCHR(16,IZN)
        HMPX = SCHR(17,IZN)
!
!---    Aqueous saturation below the aqueous matching point,
!       use Webb extension  ---
!
        IF( SLX.LT.SMPX ) THEN
          HCX = EXP( (SLX*(LOG(HMPX)-LOG(HDOD)) + SMPX*LOG(HDOD))/SMPX )
        ELSE
          SLRX = SCHR(4,IZN)
          ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
          HCX = (((1.D+0/ESLX)**(1.D+0/CMX)-1.D+0)**(1.D+0/CNX))/CAX
        ENDIF
!
!---  Brooks and Corey saturation function  ---
!
      ELSEIF( ISCHR(IZN).EQ.2 .OR. ISCHR(IZN).EQ.12 .OR.
     &  ISCHR(IZN).EQ.14 ) THEN
        CPX = SCHR(2,IZN)
        CLX = MAX( SCHR(3,IZN),SMALL )
        SMPX = SCHR(16,IZN)
        HMPX = SCHR(17,IZN)
!
!---    Aqueous saturation below the aqueous matching point,
!       use Webb extension  ---
!
        IF( SLX.LT.SMPX ) THEN
          HCX = EXP( (SLX*(LOG(HMPX)-LOG(HDOD)) + SMPX*LOG(HDOD))/SMPX )
        ELSE
          SLRX = SCHR(4,IZN)
          ESLX = MIN( MAX( (SLX-SLRX)/(1.D+0-SLRX),0.D+0 ),1.D+0 )
          HCX = CPX*(1.D+0/ESLX)**(1.D+0/CLX)
        ENDIF
      ENDIF
!
!---  Kelvin's equation  ---
!
      PVBX = PSBX*EXP( -WTMW*HCX*GRAV/(RCU*TKX) )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VPL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VPL_B( TX,PSBX,PCX,RHOBX,PVBX,XLSX )
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
!     Vapor pressure lowering of brine.
!
!     Battistelli, A., C. Claudio, and K. Pruess.  1997.  "The simulator
!     TOUGH2/EWASG for modelling geothermal reservoirs with brines and
!     noncondensible gas."  Geothermics, 26(4): 437-464.
!
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 April 2002.
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VPL_B'
!
!---  Kelvin's equation  ---
!
      TKX = TX + TABS
      WTMX = XLSX*WTMS + (1.D+0-XLSX)*WTMW
      HDCX = MIN( HDOD,(PCX/RHORL/GRAV) )
      HDRX = (HDOD-HDCX)/HDOD
      PVBX = PSBX*EXP( -WTMX*PCX/(RHOBX*RCU*TKX) )*HDRX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VPL_B group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VPL_BX( TX,PLX,PSOX,PVWX,PSWX,RHOBX,XLSX )
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
!     Vapor pressure lowering of brine.
!
!     Battistelli, A., C. Claudio, and K. Pruess.  1997.  "The simulator
!     TOUGH2/EWASG for modelling geothermal reservoirs with brines and
!     noncondensible gas."  Geothermics, 26(4): 437-464.
!
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 29 August 2006
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYST
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
      SUB_LOG(ISUB_LOG) = '/VPL_BX'
!
!---  Newton scheme on Kevin's equation  ---
!
      TKX = TX + TABS
      WTMX = XLSX*WTMS + (1.D+0-XLSX)*WTMW
      PVWX = PSWX
  100 CONTINUE
      PCX = MAX( BGL*(PSOX+PVWX-PLX-PATM),0.D+0 )
      HDCX = MIN( HDOD,(PCX/RHORL/GRAV) )
      HDRX = (HDOD-HDCX)/HDOD
      FX = PVWX - PSWX*EXP( -WTMX*PCX/(RHOBX*RCU*TKX) )*HDRX
      DFX = 1.D+0 + (PSWX*WTMX/(RHOBX*RCU*TKX))*
     &  EXP( -WTMX*PCX/(RHOBX*RCU*TKX) )*HDRX
     &  + PSWX*EXP( -WTMX*PCX/(RHOBX*RCU*TKX) )/(RHORL*GRAV*HDOD)
      DPVWX = -FX/DFX
      PVWX = PVWX+DPVWX
      IF( ABS(DPVWX).GT.1.D-6 ) GOTO 100
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VPL_BX group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WATDFG( XMW,XMO,XMA,DFGAW,DFGOW,DFGWX )
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
!     Calculates water diffusion coefficient by the Wilke
!     method. pp. 34.
!
!     Falta, R.W., K. Pruess, I. Javandel, and P.A. Witherspoon. 1990.
!     Numerical Modeling of Steam Injection for the Removal of
!     Nonaqueous Phase Liquids from the Subsurface: 1 Numerical
!     Formulation. LBL-29615, Lawrence Berkeley Laboratory.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WATDFG'
      IF( (1.D+0-XMW) .LE. 1.D-20 ) THEN
        DFGWX = 2.D+0*DFGAW*DFGOW/(DFGAW+DFGOW)
      ELSE
        DFGWX = (1.D+0-XMW)/((XMA/DFGAW) + (XMO/DFGOW))
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WATDFG group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WATLQH( TX,PX,H )
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
!     Calculate the subcooled or saturated enthalpy, as a function of
!     temperature and pressure per the steam table equations
!     as given by the 1967 International Formulation Committee:
!     Formulation for Industrial Use.
!
!     Thermodynamic and Transport Properties of Steam.
!     1967. ASME Steam Tables.
!     The American Society of Mechanical Engineers.
!     United Engineering Center, 345 East 47th Street, New York, N.Y.
!
!     The temperature is limited in this subroutine to the following
!     values:  0.01 C < T > 364.0 !
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
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
      REAL*8 A(23),SA(12),B(4)
!
!----------------------Data Statements---------------------------------!
!
      DATA A /6.824687741D+3,-5.422063673D+2,-2.096666205D+4,
     &3.941286787D+4,-6.733277739D+4,9.902381028D+4,-1.093911774D+5,
     &8.590841667D+4,-4.511168742D+4,1.418138926D+4,-2.017271113D+3,
     &7.982692717D+0,-2.616571843D-2,1.522411790D-3,2.284279054D-2,
     &2.421647003D+2,1.269716088D-10,2.074838328D-7,2.174020350D-8,
     &1.105710498D-9,1.293441934D+1,1.308119072D-5,6.047626338D-14/
      DATA SA /8.438375405D-1,5.362162162D-4,1.720,7.342278489D-2,
     &4.975858870D-2,6.537154300D-1,1.15D-6,1.1508D-5,1.4188D-1,
     &7.002753165D+0,2.995284926D-4,2.040D-1/
      DATA B / -6.36443D+5,1.14468D+2,3.64104D+0,-3.3561464D+5 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WATLQH'
!
!---  Liquid water enthalpy extended below freezing point  ---
!
      TR = (MAX(1.D-2,TX)+TABS)/TCRW
      PR = PX/PCRW
      CY = 1.0D+0 - SA(1)*TR*TR-SA(2)/(TR**6)
      CDY = -2.D+0*SA(1)*TR + 6.D+0*SA(2)/(TR**7)
      CZ = CY +
     &  SQRT(MAX( ZERO,(SA(3)*CY*CY-2.D+0*SA(4)*TR+2.D+0*SA(5)*PR)))
      HX = TR*A(1)
      DO 100 I = 1,10
        HX = HX - (I-2)*A(I+1)*(TR**(I-1))
  100 CONTINUE
      HX = HX + A(12)*(CZ*(1.7D+1*(CZ/2.9D+1 - CY/1.2D+1) +
     &  4.1667D-1*TR*CDY) + SA(4)*TR - (SA(3) - 1.D+0)*TR*CY*CDY)/
     &  (CZ**2.9412D-1)
      HX = HX + (A(13) - A(15)*TR*TR
     &  + A(16)*(9.D+0*TR + SA(6))*((SA(6) - TR)**9)
     &  + A(17)*(2.D+1*(TR**19) + SA(7))/((SA(7) + (TR**19))**2))*PR
      HX = HX - ((1.2D+1*(TR**11) + SA(8))/((SA(8) + (TR**11))**2))*
     &  (A(18)*PR + A(19)*PR*PR + A(20)*PR*PR*PR)
      HX = HX + A(21)*(TR**18)*(1.7D+1*SA(9) + 1.9D+1*TR*TR)*
     &  ((1.D+0/((SA(10) + PR)**3) + SA(11)*PR))
      HX = HX + A(22)*SA(12)*PR*PR*PR + 2.1D+1*A(23)*(PR**4)/(TR**20)
      H = HX*PCRW*VCRW*1.D-3/WTMW
!
!---  Extend the enthalpy function for freezing conditions
!     by incorporating the latent heat of fusion over the
!     temperature range -0.01 C to -1.0 C and using a polynomial
!     fit for the enthalpy of ice, referenced to -1.0 C  ---
!
      IF( TX.LT.1.D-2 ) THEN
        H = H + (MAX(-1.D+0,TX)-1.D-2)*333.7D+3/1.01D+0
      ENDIF
      IF( TX.LT.-1.D+0 ) THEN
        TK = TX + TABS
        H = H + B(1) + B(2)*TK + B(3)*TK**2 - B(4)
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WATLQH group  ---
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WEBB_BC( N )
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
!     Webb saturation and capillary pressure matching points for
!     the Brooks-Corey capillary pressure-saturation
!     function.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 16 July 2010
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
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
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WEBB_BC'
      IZN = IZ(N)
!
!---  Find the saturation matching point  ---
!
      SRX = SCHR(4,IZN)
      PSIX = SCHR(2,IZN)
      CLX = MAX( SCHR(3,IZN),SMALL )
!
!---  Initial guess of the matching-point head  ---
!
      SMPX = 1.D-2*(1.D+0-SRX) + SRX
      ESMPX = (SMPX-SRX)/(1.D+0-SRX)
      HMPX = PSIX/(ESMPX**(1.D+0/CLX))
      HMPX = 1.D+1**(4.D-1*LOG10(HDOD))
!
!---  Newton-Raphson iteration for the matching-point head  ---
!
      NC = 0
  100 CONTINUE
      NC = NC + 1
      DHX = MAX( 1.D-4,1.D-6*HMPX )
      DHX = SIGN( DHX,5.D-1*HDOD-HMPX )
      AX = -((PSIX/HMPX)**CLX)*(1.D+0-SRX)
      BX = HMPX*(LOG(HDOD)-LOG(HMPX))
      F1X = (AX-SRX)/BX - AX*CLX/HMPX
      HMPY = HMPX + DHX
      AX = -((PSIX/HMPY)**CLX)*(1.D+0-SRX)
      BX = HMPY*(LOG(HDOD)-LOG(HMPY))
      F2X = (AX-SRX)/BX - AX*CLX/HMPY
      DFX = (F2X-F1X)/DHX
      DHMPX = -F1X/DFX
      HMPX = MIN( MAX( HMPX+DHMPX,PSIX ),HDOD-1.D+0 )
!
!---  No convergence on saturation matching point  ---
!
      IF( NC.GT.64 ) THEN
        INDX = 7
        RLMSGX = 0.D+0
        IMSGX = N
        NMSGX = N
        SUBLOGX = 'WEBB_BC'
        CHMSGX = 'No Convergence on Brooks and Corey'
     &    // ' Matching Point Head @ Node: '
        CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
      ENDIF
      IF( ABS(DHMPX).GT.1.D-4 ) GOTO 100
!
!---  Find the capillary matching point saturation  ---
!
      ESMPX = (PSIX/HMPX)**CLX
      SMPX = ESMPX*(1.D+0-SRX) + SRX
!
!---  Load the matching point saturation and head  ---
!
      SCHR(16,IZN) = SMPX
      SCHR(17,IZN) = HMPX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WEBB_BC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WEBB_VG( N )
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
!     Webb saturation and capillary pressure matching points for
!     the van Genuchten capillary pressure-saturation
!     function.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 16 July 2010
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
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
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WEBB_VG'
      IZN = IZ(N)
!
!---  Find the saturation matching point  ---
!
      SRX = SCHR(4,IZN)
      ALPHAX = SCHR(1,IZN)
      CNX = SCHR(3,IZN)
      CMX = SCHR(14,IZN)
!
!---  Initial guess of the matching-point head  ---
!
      SMPX = 1.D-2*(1.D+0-SRX) + SRX
      ESMPX = (SMPX-SRX)/(1.D+0-SRX)
      HMPX = (1.D+0/ALPHAX)*
     &  (((1.D+0/(ESMPX**(1.D+0/CMX)))-1.D+0)**(1.D+0/CNX))
      HMPX = 1.D+1**(4.D-1*LOG10(HDOD))
!
!---  Newton-Raphson iteration for the matching-point head  ---
!
      NC = 0
  100 CONTINUE
      NC = NC + 1
      DHX = MAX( 1.D-4,1.D-6*HMPX )
      DHX = SIGN( DHX,5.D-1*HDOD-HMPX )
      AX = ((1.D+0 + (ALPHAX*HMPX)**CNX)**(-CMX))*(1.D+0-SRX) + SRX
      BX = HMPX*(LOG(HDOD)-LOG(HMPX))
      CX = (1.D+0 + (ALPHAX*HMPX)**CNX)**(-CMX-1.D+0)
      DX = CNX*CMX*((ALPHAX*HMPX)**CNX)*(1.D+0-SRX)/HMPX
      F1X = AX/BX - CX*DX
      HMPY = HMPX + DHX
      AX = ((1.D+0 + (ALPHAX*HMPY)**CNX)**(-CMX))*(1.D+0-SRX) + SRX
      BX = HMPY*(LOG(HDOD)-LOG(HMPY))
      CX = (1.D+0 + (ALPHAX*HMPY)**CNX)**(-CMX-1.D+0)
      DX = CNX*CMX*((ALPHAX*HMPY)**CNX)*(1.D+0-SRX)/HMPY
      F2X = AX/BX - CX*DX
      DFX = (F2X-F1X)/DHX
      DHMPX = -F1X/DFX
      HMPX = MIN( MAX( HMPX+DHMPX,1.D-6 ),HDOD-1.D+0 )
!
!---  No convergence on saturation matching point  ---
!
      IF( NC.GT.64 ) THEN
        INDX = 7
        RLMSGX = 0.D+0
        IMSGX = N
        NMSGX = N
        SUBLOGX = 'WEBB_VG'
        CHMSGX = 'No Convergence on van Genuchten'
     &    // ' Matching Point Head @ Node: '
        CALL WRMSGX( RLMSGX,SUBLOGX,CHMSGX,IMSGX,NMSGX,INDX )
      ENDIF
      IF( ABS(DHMPX).GT.1.D-4 ) GOTO 100
!
!---  Find the capillary matching point saturation  ---
!
      ESMPX = (1.D+0 + (ALPHAX*HMPX)**CNX)**(-CMX)
      SMPX = ESMPX*(1.D+0-SRX) + SRX
!
!---  Load the matching point saturation and head  ---
!
      SCHR(16,IZN) = SMPX
      SCHR(17,IZN) = HMPX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WEBB_VG group  ---
!
      RETURN
      END

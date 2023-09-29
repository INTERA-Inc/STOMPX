!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRNDBH( TX,PX,APHI )
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
!     Compute the Debye-Huckel parameter as a function of temperature
!     and pressure. Below 373.15 K (100 oC), the parameters
!     apply to a standard pressure of 0.1 MPa (1 bar). Above 373.15 K,
!     they apply to water under its own vapor presssure.
!
!     Bradley, D.J., and K.S.Pitzer. 1983. Thermodynamics of 
!     electrolytes XII. Dielectric proprties of water and Debye-Huckel 
!     parameters to 350 oC and 1kbar.
!     J. Phys. Chem 83:1599,1979;87:3798.
!
!     ELC = electronic charge,
!     ANO = Avagadro's number, molecules/mole
!     BZM = Boltzmann's constant,
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, November 1995.
!     Last Modified by MD White, Battelle, November 16, 1995.
!     Last Modified by AL Ward, Battelle, July 23, 1999.
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
      SAVE ELC,ANO,BZM,EPO
      DATA ELC,ANO,BZM / 1.60210D-19,6.02252D+23,1.38054D-23 /
      DATA EPO / 8.854D-12 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BRNDBH'
!
!---  Convert input temperature to absolute temperature ---
!
        TK = TX+TABS
!
!---  Compute the liquid water density as a function of P and T ---
!
        CALL WATLQD( TX,PX,RHO )
!
!---  Compute relative permittivity of pore water at temperature T---
!
        ER = EPR(TX)
!
!---  Compute Debye-Huckel parameter for the osmotic coefficient ---
!
        VAR1 = SQRT(2.D+0*GPI*ANO*RHO)
        VAR2 = ((ELC**2)/(4.D+0*GPI*EPO*ER*BZM*TK))**1.5D+0
        APHI = VAR1*VAR2/3.D+0
!
!---  End of BRNDBH group  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!
!----------------------Subroutine--------------------------------------!
!
      FUNCTION EPR( TX )
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
!     Compute the the relative dielectric constant as a function of
!     temperature assuming a pressure of (</=100 MPa, or 1 kbar). At
!     T< 373.15 K (100 oC), the parameters apply to a standard pressure
!     of 0.1 MPa (1 bar). At T >373.15 K, they apply to water under its
!     own vapor presssure.
!
!     Bradley, D.J., and K.S.Pitzer. 1983. Thermodynamics of 
!     electrolytes XII. Dielectric proprties of water and Debye-Huckel 
!     parameters to 350 oC and 1kbar. 
!     J. Phys. Chem 83:1599,1979;87:3798.
!
!----------------------Authors-----------------------------------------!
!
!     Written by AL Ward, Battelle, PNNL, July 23, 1999.
!     Last Modified by AL Ward, Battelle, July 24, 1999.
!
!----------------------Fortran 90 Modules------------------------------!
!
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
      REAL*8 U(9)
!
!----------------------Data Statements---------------------------------!
!
      SAVE U
      DATA U / 3.4279D+2, -5.0866D-3, 9.4690D-7, -2.0525D+0, 3.1159D+3,
     &  -1.8289D+2, -8.0325D+3, 4.2142D+6, 2.1417D+0 /
!
!----------------------Executable Lines--------------------------------!
!
!
      TK = TX+TABS
      EPK = U(1)* EXP(U(2)*TK+U(3)*TK*TK)
      CX = U(4)+U(5)/(U(6)+TK)
      BX = U(7)+U(8)/TK+U(9)*TK
      VPX = PATM
      IF( TK.GT.3.73149D+2 ) CALL WATSP( TX,VPX )
      EPR = EPK+CX*LOG((BX+VPX*1.D-5)/(BX+1.D+3))
!
!---  End of EPR group  ---
!
      RETURN
      END


!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRNOSM ( TX,XLSX,APHI,BOC,SOC,PVWX,POSX )
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
!     Compute the the osmotic coefficient as a function of salt mass
!     fraction and temperature using the Pitzer equation
!
!     Pitzer, K J.C. 1991. Ion interaction approach: Theory and data
!     correlation. In Activity coefficeints in electrolyte solutions,
!     2nd Ed.., CRC Press: Boston.
!
!     Li, Z-B,Y-G Li, and J-F Lu. 1999. Surface Tension Model for
!     Concentrated Electrolyte Aqueous Solutions by the Pitzer Equation.
!     Ind. Eng. Chem. Res. 38, 1133-1139
!
!     PTPS(1) = B0
!     PTPS(2) = B1
!     PTPS(3) = B2
!     PTPS(4) = CPHI
!     PTPS(5) = GOC
!     IPTPS(1) = IVPOS
!     IPTPS(2) = IVNEG
!     IPTPS(3) = IZPOS
!     IPTPS(4) = IZNEG
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, November 1995.
!     Last Modified by MD White, Battelle, November 16, 1995.
!     Last Modified by AL Ward, Battelle, July 23, 1999.
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
      SAVE SOMV
      DATA SOMV / 1.805D-5 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BRNOSM'
      TK = TX+TABS
!
!---  Compute stoichiometric and charge balance terms ---
!
      ZPOS = REAL(IPTPS(3))
      ZNEG = REAL(IPTPS(4))
      VPOS = REAL(IPTPS(1))
      VNEG = REAL(IPTPS(2))
      ZMX = ABS(ZPOS*ZNEG)
      VPN = VPOS+VNEG
      VMX = VPOS*VNEG
      OCB = 1.2D+0
!
!---  Compute mole fraction and molality in bulk and surface phases ---
!
      XMLSX = XLSX/WTMS/(XLSX/WTMS + (1.D+0-XLSX)/WTMW)
      XGLSX = 1.D+3*XMLSX/(WTMW - XMLSX*WTMW)
      SXGLSX = PTPS(5)*XGLSX
!
!---  Compute the bulk phase ionic strength and its square root ---
!
      BSI = 5.D-1*(XGLSX*ZPOS*ZPOS + XGLSX*ZNEG*ZNEG)
      DBSI = SQRT(BSI)
      SSI = 0.5D+0*(SXGLSX*ZPOS*ZPOS + SXGLSX*ZNEG*ZNEG)
      DSSI = SQRT(SSI)
!
!---  Compute bulk phase osmotic coefficent ---
!
      IF( XGLSX.GT.EPSL ) THEN
        IF( IPTPS(1).EQ.1 .AND. IPTPS(2).GE.1 .AND. IPTPS(2).LE.3 ) THEN
          A1 = 2.D+0
          BOC = 1.D+0 + ZMX*APHI*DBSI/(1.D+0+OCB*DBSI) +
     &      2.D+0*XGLSX*VMX*(PTPS(1)+PTPS(2)*EXP(-A1*DBSI))/VPN +
     &      2.D+0*(XGLSX**2)*(VMX**1.5D+0)*PTPS(4)/VPN
        ELSEIF( IPTPS(1).GE.2 .AND. IPTPS(2).GE.2 ) THEN
          A1 = 1.4D+0
          A2 = 1.2D+1
          BOC = 1.D+0 + ZMX*APHI*DBSI/(1.D+0+OCB*DBSI) +
     &      2.D+0*XGLSX*VMX*(PTPS(1)+PTPS(2)*EXP(-A1*DBSI))/VPN +
     &      2.D+0*(XGLSX**2)*(VMX**1.5D+0)*PTPS(4)/VPN +
     &      PTPS(3)*EXP(-A2*DBSI)
        END IF
      ELSE
        BOC = 0.D+0
      END IF
!
!---  Compute surface phase osmotic coefficent ---
!
      IF( XGLSX.GT.EPSL ) THEN
        IF( IPTPS(1).EQ.1 .AND. IPTPS(2).GE.1 .AND. IPTPS(2).LE.3 ) THEN
          A1 = 2.D+0
          SOC = 1.D+0 + ZMX*APHI*DSSI/(1.D+0+OCB*DSSI) +
     &      2.D+0*SXGLSX*VMX*(PTPS(1)+PTPS(2)*EXP(-A1*DSSI))/VPN +
     &      2.D+0*(SXGLSX**2)*(VMX**1.5D+0)*PTPS(4)/VPN
        ELSEIF( IPTPS(1).GE.2 .AND. IPTPS(2).GE.2 )THEN
          A1 = 1.4D+0
          A2 = 1.2D+1
          SOC = 1.D+0 + ZMX*APHI*DSSI/(1.D+0+OCB*DSSI) +
     &      2.D+0*SXGLSX*VMX*(PTPS(1)+PTPS(2)*EXP(-A1*DSSI))/VPN +
     &      2.D+0*(SXGLSX**2)*(VMX**1.5D+0)*PTPS(4)/VPN +
     &      PTPS(3)*EXP(-A2*DSSI)
        END IF
      ELSE
        SOC = 0.D+0
      END IF
!
!---  Compute overall osmotic coefficient
!
      OC = BOC - PTPS(5)*SOC
!
!---  Compute the osmotic pressure of the pore water ---
!
      POSX = 0.D+0
      IF( ISLC(7).EQ.1 .OR. ISLC(7).EQ.2 ) THEN
        POSX = -OC*(RCU/1.D+3)*TK*VPN*XGLSX*WTMW/(SOMV*1.D+3)
      END IF
!
!---  Compute the vapor pressure of the pore water ---
!
      AAC = EXP(-3.D+0*OC*WTMW*XGLSX/1.D+3 )
      PVWX =PVWX*AAC
!
!---  End of BRNOSM group  ---
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRNSTN( TX,XLSX,SOC,BOC,ST )
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
!     Compute the pore water surface tension as a function of salt mass
!     fraction and temperature using the Pitzer equation
!
!     Adamson, A.W. 1976. The physical Chemistry of Surfaces. 3rd. Ed
!     Wiley--Interscience, NY.
!
!     SOMV = water molar volume, m^3/mol
!     STW = water-air surface tension, N/m
!     PTPS(5) = GOC
!
!----------------------Authors-----------------------------------------!
!
!     Written by AL Ward, Battelle, PNL, July 17, 1999.
!     Last Modified by AL Ward, Battelle, PNL, July 19, 1999.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_FD
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
      SAVE SOMO,SOMV,STW,ANO
      DATA SOMO,SOMV,STW,ANO / 55.50837D+1,1.805D-5,7.2D-2,6.02252D+23 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BRNSTN'
!
!---  Compute the molar surface area of pure water ---
!
      SAMW = (ANO**(THIRD))*(SOMV**(2.D+0/3.D+0))
!
!---  Compute mole fraction and molality from wt % in the bulk phase ---
!
      XMLSX = (XLSX/WTMS)/(XLSX/WTMS + (1.D+0-XLSX)/WTMW)
      XGLSX = 1.D+3*XMLSX/(WTMW - XMLSX*WTMW)
!
!---  Compute the van't Hoff coefficient ---
!
      TK = TX + TABS
      VPOS = REAL(IPTPS(1))
      VNEG = REAL(IPTPS(2))
      VHOF =((VPOS+VNEG)*(RCU/1.D+3)*TK*XGLSX)/(SOMO*SAMW)
!
!---  Compute the contribution of the bulk liquid phase ---
!
      STB = VHOF*BOC
!
!---  Compute the contribution of the surface phase ---
!
      STS = VHOF*PTPS(5)*SOC
!
!---  Compute the surface tension of solution ---
!
      ST = STW + (STB - STS)
!
!---  End of BRNSTN group
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END


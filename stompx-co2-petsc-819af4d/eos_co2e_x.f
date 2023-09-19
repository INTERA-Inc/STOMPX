!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DENS_A( TX,PX,RHOAX,I_VX )
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
!     Density of liquid or vapor CO2.
!
!     ISRX liquid or vapor index: 1 - liquid 2 - vapor or supercritical
!
!     Span, R., and W. Wagner.  1996.  A New Equation of State for
!     Carbon Dioxide Covering the Fluid Region from the Triple-Point
!     to 1100 K at Pressures up to 800 MPa.  J. Phys. Chem. Ref. Data
!     25(6):1509-1588.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 April 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER I_PX(2),I_SX,I_TX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DENS_A'
      INCG = 1
      PTX = MAX( PX,P_TA(1,INCG)*1.D+6 )
      CALL ITL_A( PTX,PTPA,PCRA,TX,I_PX,I_SX,I_TX,I_VX )
      CALL PTL_A( PTX,TX,RHO_TA,RHOAX,DRHOAX,I_PX,I_SX,I_TX,I_VX )
!
!---  Ideal gas law  ---
!
      IF( PX.LT.PTX ) THEN
        RHOAX = PX*RHOAX/PTX
      ENDIF
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
!     Written by M.D. White, PNNL, 4 April 2002.
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
!---  Compressed or vapor-saturated density of pure water using the
!     ASME formulations  ---
!
      CALL SP_W( TX,PSWX )
      ISRX = 1
      PWX = MAX( PX,PSWX )
      CALL DENS_W( TX,PWX,RHOLWX,RHOVWX,ISRX )
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
!      SXX = SCX(1) + SCX(2)*EXP(SAX(2)*TX) +
!     &  SCX(3)*EXP(SAX(3)*PBWX)
!      RHOLWPX = CCX(1) + CCX(2)*SXX + CCX(3)*(SXX**2) + CCX(4)*(SXX**3)
!
!---  Density of NaCl brine (gm/cm^3) using the Phillips et al.
!     formulations  ---
!
!      SXX = SCX(1)*EXP(SAX(1)*GLSX) + SCX(2)*EXP(SAX(2)*TX) +
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
      SUBROUTINE DENS_L( TX,RHOBX,XLAX,RHOLX )
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
!     Alendal, G., and H. Drange.  2001.  Two-phase, near-field,
!     modeling of purposefully released CO2 in the ocean.  Journal
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
!     Written by M.D. White, PNNL, 12 April 2002.
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
!----------------------TYPE DECLARATIONS-------------------------------!
!
      REAL*8 PMV_A
      REAL*8 PMV_C(5)
!
!----------------------DATA STATEMENTS---------------------------------!
!
!      SAVE PMV_A
!      DATA PMV_A / 34.D-3 /
!
!----------------------Data Statements---------------------------------!
!
!      SAVE PMV_C
      DATA PMV_C / 37.36D-3,-7.109D-5,-3.812D-8,3.296D-9,-3.702D-12 /
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(9).EQ.1 ) THEN
        RHOLX = RHOLI
!
!---    Reset subroutine string sequence  ---
!
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DENS_L'
!
!---  Partial molar volume of CO2 (m^3/kmol), formulation of
!     Anderson et al. (1992)  ---
!
      PMV_A = 0.D+0
      DO 10 M = 1,5
        PMV_A = PMV_A + PMV_C(M)*(TX**(M-1))
   10 CONTINUE
      CXA = PMV_A*RHOBX*XLAX/WTMA
!
!---  Poynting correction  ---
!
      RHOLX = RHOBX/(1.D+0 + CXA - XLAX)
!      RHOLX = RHOBX + WTMA*CX - CX*RHOBX*PMV_A
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
!     Battistelli, A., C. Claudio, and K. Pruess.  1997.  The simulator
!     TOUGH2/EWASG for modelling geothermal reservoirs with brines and
!     noncondensible gas.  Geothermics, 26(4): 437-464.
!
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 1 May 2002.
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
!     Written by M.D. White, PNNL, 26 March 2002.
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
      REAL*8 CAX(23),SAX(12),CBX(31),SBX(5),LFCX(3),SLX
      REAL*8 TX,PX,RHOLX,RHOVX
      INTEGER ISRX,INX(8),IZX(8,3),ITX(8),IXX(8,2)
!
!----------------------Data Statements---------------------------------!
!
!      SAVE CAX,SAX,CBX,SBX
!      SAVE LFCX,SLX,INX,IZX,ITX,IXX
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
!---  Check temperature and pressure  ---
!
      IF( TX.LT.1.D-2 .OR. TX.GT.8.D+2 ) THEN
        TX = MAX( TX,1.D-2 )
        TX = MIN( TX,8.D+2 )
      ENDIF
      IF( PX.LT.0.D+0 .OR. PX.GT.1.D+8 ) THEN
        PX = MAX( PX,0.D+0 )
        PX = MIN( PX,1.D+8 )
      ENDIF
      IF( PX.LT.EPSL ) THEN
        RHOLX = 0.D+0
        RHOVX = 0.D+0
        GOTO 100
      ENDIF
!
!---  Reduced temperature and pressure  ---
!
      RHOLX = 0.D+0
      RHOVX = 0.D+0
      THETAX = (TX+TABS)/TCRW
      BETAX = PX/PCRW
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
        M_ERR(1) = 'Steam Table Subregion 3/5: Exceeds Temperature ' // 
     &    'Range: '
        M_ERR(2) = ' Temperature = '
        IF( N_DB.LT.0 ) THEN
          M_ERR(3) = ' at Boundary Surface: '
        ELSE
          M_ERR(3) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = TX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 2
        I_ERR(3) = 3
        I_ERR(4) = ID
      ENDIF
!
!---  Subregions 4 or 5  ---
!
      IF( ISRX.EQ.4 .OR. ISRX.EQ.5 ) THEN
        M_ERR(1) = 'Steam Table Subregion 4/5: Exceeds Temperature ' // 
     &    'Range: '
        M_ERR(2) = ' Temperature = '
        IF( N_DB.LT.0 ) THEN
          M_ERR(3) = ' at Boundary Surface: '
        ELSE
          M_ERR(3) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = TX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 2
        I_ERR(3) = 3
        I_ERR(4) = ID
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
      SUBROUTINE DIFC_GW( TX,PX,DFGWX )
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
!     Calculates the CO2-H2O binary diffusion coefficient from
!     the Wilke and Lee theory. pp. 587
!
!     Liquid molar volume is computed from critical molar volume with
!     the Tyn and Calus method. pp. 53-54
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids. pp. 587.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 29 April 2002.
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
      REAL*8 CAX(8),CLJ_W(2),CLJ_A(2)
!
!----------------------Data Statements---------------------------------!
!
!      SAVE CAX,CLJ_W,CLJ_A
      DATA CAX /1.06036D+0,1.5610D-1,1.9300D-1,4.7635D-1,1.03587D+0,
     &  1.52996D+0,1.76474D+0,3.89411D+0/
      DATA CLJ_W / 3.190008977D+0, 429.18D+0 /
      DATA CLJ_A / 3.795165630D+0, 95.85245D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DIFC_GW'
!
!---  Absolute temperature and pressure in bars  ---
!
      TKX = TX+TABS
      PBX = PX*1.D-5
!
!---  Mean characteristic Lennard-Jones energy and length (angstrom) ---
!
      EPSX = SQRT(CLJ_W(2)*CLJ_A(2))
      SIGX = 5.D-1*(CLJ_W(1)+CLJ_A(1))
      TPX = TKX/EPSX
      OMGX = (CAX(1)/TPX**CAX(2)) + (CAX(3)/EXP(CAX(4)*TPX))
     &  + (CAX(5)/EXP(CAX(6)*TPX))  + (CAX(7)/EXP(CAX(8)*TPX))
      WTMX = 2.D+0/((1.D+0/WTMA) + (1.D+0/WTMW))
      DFGWX = (3.03D+0 - (9.8D-1/SQRT(WTMX)))*1.D-3*(TKX**1.5D+0)/
     &  (PBX*SQRT(WTMX)*(SIGX**2)*(OMGX))*1.D-4
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DIFC_GW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DIFC_LA( TX,XLSX,DFLAX )
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
!     Cadogan, S.P., G.C. Maitland, and J.P. Martin Trulser. 2014.
!     Diffusion Coefficients of CO2 and N2 in Water at Temperatures
!     between 298.15 K and 423.15 K at Pressures up to 45 MPa.
!     dx.doi.org/10.1021/je401008s
!     J. Chem. Eng. Data 2014, 59, 519−525
!
!     Belgodere, C. et al. 2014. Experimental determination of CO2 
!     diffusion coefficient in aqueous solutions under pressure via 
!     Raman spectroscopy at room temperature: impact of salinity 
!     (NaCl) on dissolved CO2 diffusivity. 
!     Abstract for 11th GeoRaman International Conference, 
!     June 15-19, 2014, St. Louis, Missourri, USA
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 1 May 2002.
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
      SUB_LOG(ISUB_LOG) = '/DIFC_LA'
!
!---  Pure water diffusion coefficient as a function of temperature
!     Cadogan et al. (2014), m^2/s  ---
!
      TKX = TX + TABS
      DFLAX = 3.5984D+0 - 6.5113D-2*TKX + 2.0282D-4*(TKX**2)
      DFLAX = DFLAX*1.D-9
!
!---  Correction for salinity (mol NaCL/kg H2O)  ---
!
      SLNTY = 1.D+3*(XLSX/WTMS)/(1.D+0-XLSX)
      DFLAX = DFLAX*(1.6678D+0 - 1.2531D-1*SLNTY)/1.6678D+0
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
!     Bromley, L.A.  1973.  Thermodynamic properties of strong
!     electrolytes in aqueous solutions.  AIChE Journal, 19(2):313-320
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 May 2002.
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
!      SAVE CBX,TRX,VISWRX
      DATA CBX / 0.0547D+0 /
      DATA TRX / 25.D+0 /
      DATA VISWRX / 0.8904339807D-3 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DIFC_LS'
!
!---  Convert halite mass fraction to molality   ---
!
      GLSX = 1.D+3*XLSX/(WTMS*(1.D+0-XLSX))
      TKX = TX+TABS
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
        DLNGX = 0.D+0
      ENDIF
!
!---  Diffusion coefficient (m^2/s) for dilute NaCl aqueous solutions
!     at 298 K   ---
!
      DFLSX = 2.254D-9
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
      SUBROUTINE ENTH_A( TX,PX,HAX,UAX,I_VX )
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
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 April 2002.
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
      INTEGER I_PX(2),I_SX,I_TX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ENTH_A'
      INCG = 1
      PTX = MAX( PX,P_TA(1,INCG)*1.D+6 )
      CALL ITL_A( PTX,PTPA,PCRA,TX,I_PX,I_SX,I_TX,I_VX )
      CALL PTL_A( PTX,TX,H_TA,HAX,DHAX,I_PX,I_SX,I_TX,I_VX )
      CALL PTL_A( PTX,TX,U_TA,UAX,DUAX,I_PX,I_SX,I_TX,I_VX )
!
!---  Convert units to J/kg  ---
!
      HAX = 1.D+3*HAX
      UAX = 1.D+3*UAX
!
!---  Ideal gas law  ---
!
      IF( PX.LT.PTX ) THEN
        TKX = TX + TABS
        HAX = UAX + RCU*TKX/WTMA
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTH_A group  ---
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
!     Michaelides, E.E.  1981. Thermodynamic properties of geothermal
!     fluids.  Geothermal Resources Council, Transactions 5:361-364.
!
!     Gudmundsson, J.S., and H. Thrainsson.  1989.  Power potential of
!     two-phase geothermal wells.  Geothermics 18(3):357-366.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 8 April 2002.
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
!      SAVE SAX,SBX,SCX,CAX,CBX
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
      SUBROUTINE ENTH_L( TX,XLSX,XLAX,HBX,HGAX,HLX )
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
!     temperature, mass fraction of aqueous NaCl, mass fraction
!     of aqueous CO2 and mass fraction of aqueous CH4.
!
!     Battistelli, A., C. Claudio, and K. Pruess.  1997.  The simulator
!     TOUGH2/EWASG for modelling geothermal reservoirs with brines and
!     noncondensible gas.  Geothermics, 26(4): 437-464.
!
!     Himmelblau, D. M.  1959.  Partial molal heats and entropies of
!     solution for gases dissolved in water from the freezing point
!     to near the critical point.  Journal of Physical Chemsitry,
!     63:1803-1808.
!
!     TX - temperature, C
!     XLSX - aqueous NaCl mass fraction
!     XLAX - aqueous CO2 mass fraction
!     HBX - enthalpy of brine (NaCl-H2O), J/kg
!     HGAX - enthalpy of gaseous CO2, J/kg
!     HLX - enthalpy of aqueous phase, J/kg
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 16 April 2002.
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
      SUB_LOG(ISUB_LOG) = '/ENTH_L'
!
!---  Partial differential of Henry's constant for CO2 at constant 
!     pressure with respect to temperature  ---
!
      DTX = 1.D-6
      CALL HC_LA( TX,XLSX,HCX )
      TY = TX + DTX
      CALL HC_LA( TY,XLSX,HCY )
      DHCX = LOG(HCY/HCX)/DTX
!
!---  Heat of solution for CO2 from Himmelblau (1959)  ---
!
      TKX = TX + TABS
      HSAX = -RCU*(TKX**2)*DHCX/WTMA
!
!---  Aqueous enthalpy  ---
!
      HLX = MAX(1.D+0-XLAX,0.D+0)*HBX + XLAX*(HGAX + HSAX)
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
!     Written by M.D. White, PNNL, 1 May 2002.
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
!      SAVE CAX
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
!     Written by M.D. White, PNNL, 26 March 2002.
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
      REAL*8 CAX(23),SAX(12),CBX(31),SBX(5),LFCX(3),SLX
      REAL*8 TX,PX
      INTEGER ISRX,INX(8),IZX(8,3),ITX(8),IXX(8,2)
!
!----------------------Data Statements---------------------------------!
!
!      SAVE CAX,SAX,CBX,SBX
!      SAVE LFCX,SLX,INX,IZX,ITX,IXX
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
!---  Restrict pressure to above 0.01 Pa and below 100 MPa  ---
!
      PY = MIN( MAX( PX,1.D-2 ),1.D+8 )
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
        M_ERR(1) = 'Steam Table Subregion 3/5: Exceeds Temperature ' // 
     &    'Range: '
        M_ERR(2) = ' Temperature = '
        IF( N_DB.LT.0 ) THEN
          M_ERR(3) = ' at Boundary Surface: '
        ELSE
          M_ERR(3) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = TX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 2
        I_ERR(3) = 3
        I_ERR(4) = ID
      ENDIF
!
!---  Subregions 4 or 5  ---
!
      IF( ISRX.EQ.4 .OR. ISRX.EQ.5 ) THEN
        M_ERR(1) = 'Steam Table Subregion 4/5: Exceeds Temperature ' // 
     &    'Range: '
        M_ERR(2) = ' Temperature = '
        IF( N_DB.LT.0 ) THEN
          M_ERR(3) = ' at Boundary Surface: '
        ELSE
          M_ERR(3) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = TX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 2
        I_ERR(3) = 3
        I_ERR(4) = ID
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
      SUBROUTINE ENTHG_I( HGX,HRGOX,PVAX,PVWX,QMX,QTHX,RHOGX,
     &  THKGX,TPX,TX,VOLW,XGAX,XGWX )
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
!     H2O-NaCl-CO2 Equation of state.  This subroutine finds the
!     temperature at a given gas enthalpy, CO2 partial pressure,
!     and water partial pressure.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 27 July 2007.
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ENTHG_I'
      ISRX = 2
      TIX = TX
!
!---  Newton scheme  ---
!
      NC = 0
  100 CONTINUE
      NC = NC + 1
      CALL DENS_W( TX,PVWX,RHOLWX,RHOGWX,ISRX )
      CALL DENS_A( TX,PVAX,RHOGAX,ISRX )
      RHOGX = RHOGWX+RHOGAX
      XGAX = RHOGAX/RHOGX
      XGWX = RHOGWX/RHOGX
      CALL ENTH_A( TX,PVAX,HGAX,UEGAX,ISRX )
      CALL ENTH_W( TX,PVWX,HLWX,HGWX,ISRX )
      HGY = XGAX*HGAX + XGWX*HGWX
      RY = QMX*(HGY-HGX) + (HGY*RHOGX-HRGOX)*VOLW*DTI + (TX-TPX)*THKGX
      TZ = TX + 1.D-7
      CALL DENS_W( TZ,PVWX,RHOLWZ,RHOGWZ,ISRX )
      CALL DENS_A( TZ,PVAX,RHOGAZ,ISRX )
      RHOGZ = RHOGWZ+RHOGAZ
      XGAZ = RHOGAZ/RHOGZ
      XGWZ = RHOGWZ/RHOGZ
      CALL ENTH_A( TZ,PVAX,HGAZ,UEGAZ,ISRX )
      CALL ENTH_W( TZ,PVWX,HLWZ,HGWZ,ISRX )
      HGZ = XGAZ*HGAZ + XGWZ*HGWZ
      RZ = QMX*(HGZ-HGX) + (HGZ*RHOGZ-HRGOX)*VOLW*DTI + (TZ-TPX)*THKGX
      FX = RY
      DFX = 1.D+7*(RZ-RY)
      DTX = -FX/DFX
      IF( ABS(TX+TABS-TCRA).LT.2.D+0 ) THEN
        DTX = SIGN( MIN(ABS(DTX),(ABS(TX+TABS-TCRA)/8.D-1)),DTX )
      ELSE
        DTX = 8.D-1*DTX
      ENDIF
      TX = TX + DTX
      IF( NC.GT.32 ) THEN
        TX = TIX
        M_ERR(1) = 'Unconverged Well Temperature: Temperature = '
        IF( N_DB.LT.0 ) THEN
          M_ERR(2) = ' at Boundary Surface: '
        ELSE
          M_ERR(2) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = TX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 1
        I_ERR(3) = 2
        I_ERR(4) = ID
      ENDIF
      IF( ABS(DTX).GT.1.D-6 ) GOTO 100
!
!---  Use the inlet enthalpy and outlet temperature to calculate 
!     energy transport into the formation  ---
!
      QTHX = QMX*HGX + (TX-TPX)*THKGX
      CALL DENS_W( TX,PVWX,RHOLWX,RHOGWX,ISRX )
      CALL DENS_A( TX,PVAX,RHOGAX,ISRX )
      RHOGX = RHOGWX+RHOGAX
      CALL ENTH_A( TX,PVAX,HGAX,UEGAX,ISRX )
      CALL ENTH_W( TX,PVWX,HLWX,HGWX,ISRX )
      HGX = XGAX*HGAX + XGWX*HGWX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTHG_I group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE EQUIL( TX,PX,PGAX,PGWX,PSBX,PVBX,
     &  XGAX,XGWX,XLAX,XLSX,XLWX,XMGAX,XMGWX,XMLAX,XMLSX,XMLWX )
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
!     This subroutine calculates the 1 state between
!     CO2 and water, following the phase-partitioning model
!     developed by Spycher and Pruess.
!
!     1 - H2O
!     2 - CO2
!
!     Spycher, N., and K. Pruess.  2010.  A phase-partitioning model
!     for CO2-brine mixtures at elevated temperatures and pressures:
!     Application to CO2-enchanced geothermal systems. Transport
!     in Porous Media, 82:173-196.
!
!     Salts
!
!     1 - Na
!     2 - K
!     3 - Ca
!     4 - Mg
!     5 - Cl
!     6 - SO4
!
!     Arguments
!
!     TX - temperature, C
!     PX - pressure, Pa abs
!     PGAX - partial pressure of CO2, Pa abs
!     PGWX - partial pressure of water, Pa abs
!     PSBX - saturated vapor pressure of brine, Pa abs
!     PVBX - reduced vapor pressure of brine, Pa abs
!     XGAX - mass fraction of CO2 in gas phase
!     XGWX - mass fraction of water in gas phase
!     XLAX - mass fraction of CO2 in aqueous phase
!     XLSX - (in) mass fraction of salt in brine (i.e., salt + water)
!     XLSX - (out) mass fraction of salt in aqueous phase
!     XLWX - mass fraction of water in aqueous phase
!     XMGAX - mole fraction of CO2 in gas phase
!     XMGWX - mole fraction of water in gas phase
!     XMLAX - mole fraction of CO2 in aqueous phase
!     XMLSX - mole fraction of salt in aqueous phase
!     XMLWX - mole fraction of water in aqueous phase
!
!     ISLC(32):  0    Nonisobrine
!                1    Isobrine
!                2    Isobrine with simplified Spycher 1
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 11 June 2010.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
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
      REAL*8 CACX(2),CAWX(2),CACWX(2)
      REAL*8 CLKWX(4),CLKCNX(4),CLKCGX(4)
      REAL*8 CAMX,CPRX
      REAL*8 CLMBX(3),CXIX(3)
      REAL*8 DACX(2),DAWX(2)
      REAL*8 DKWCX(2),DKCWX(2)
      REAL*8 DLKWX(5),DLKCX(5)
      REAL*8 DVMCX(2),DVMWX(2)
      REAL*8 DAMX(2),DPRX(5)
      REAL*8 AIX(2,2),BIX(2),FUGX(2),YX(2),XX(2),CKX(2,2)
      REAL*8 SKX(2,2),SAX(2),ACTVX(2)
      REAL*8 FUGLX(2),FUGHX(2)
      REAL*8 YZ(2),XZ(2)
      REAL*8 AJ(2,2),BJ(2),GX(2,3)
      REAL*8 CXMGW(17),CXMLA(8)
      INTEGER IJ(2)
!
!----------------------Data Statements---------------------------------!
!
      DATA CACX / 7.54D+7,-4.13D+4 /
      DATA CAWX / 0.D+0,0.D+0 /
      DATA CACWX / 7.89D+7,0.D+0 /
      DATA CBCX / 27.8D+0 /
      DATA CBWX / 18.18D+0 /
      DATA CLKWX / -2.209D+0,3.097D-2,-1.098D-4,2.048D-7 /
      DATA CLKCNX / 1.169D+0,1.368D-2,-5.380D-5,0.D+0 /
      DATA CLKCGX / 1.189D+0,1.304D-2,-5.446D-5,0.D+0 /
      DATA CVCX / 32.6D+0 /
      DATA CVWX / 18.1D+0 /
      DATA CAMX / 0.D+0 /
      DATA CPRX / 1.D+0 /
      DATA CLMBX / 2.217D-4,1.074D+0,2.648D+3 /
      DATA CXIX / 1.3D-5,-2.012D+1,5.259D+3 /
      DATA WTMNAX / 22.9898D+0 /
      DATA WTMCLX / 35.453D+0 /
      DATA DACX / 8.008D+7,-4.984D+4 /
      DATA DAWX / 1.337D+8,-1.4D+4 /
      DATA DBCX / 28.25D+0 /
      DATA DBWX / 15.70D+0 /
      DATA DKWCX / 1.427D-2,-4.037D-4 /
      DATA DKCWX / 0.4228D+0,-7.422D-4 /
!      DATA DLKWX / -2.209D+0,3.097D-2,-1.098D-4,2.048D-7,0.D+0 /
      DATA DLKWX / -2.1077D+0,2.8127D-2,-8.4298D-5,1.4969D-7,
     &             -1.1812D-10 /
      DATA DLKCX / 1.668D+0,3.992D-3,-1.156D-5,1.593D-9,0.D+0 /
      DATA DVMCX / 32.6D+0,3.413D-2 /
      DATA DVMWX / 18.1D+0,3.137D-2 /
      DATA DAMX / -3.084D-2,1.927D-5 /
      DATA DPRX / -1.9906D-1,2.0471D-3,1.0152D-4,-1.4234D-6,1.4168D-8 /
      DATA CXMGW / 3.2217D-3,1.225D-8,3.067D+0,-9.7198D-3,5.1621D+0,
     &             3.485D+2,7.7053D+1,1.0928D-2,3.663D+2,-1.9472D+0,
     &             1.3937D+0,2.4992D+1,2.5343D+2,1.4677D+1,3.7952D-2,
     &             2.2122D+3,-1.8936D+0 /
      DATA CXMLA / 4.1265D-2,1.0715D+1,2.6542D+1,2.8708D+2,
     &             2.5478D-2,-3.0218D-4,1.3776D-6,-2.2457D-09 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/EQUIL'
      INCG = 1
!
!---  Temperature in Kelvin, pressure in bar, 
!     gas constant in cm^3 bar/K mol   ---
!
      TKX = TX+TABS
      PBX = MAX(PX,PATM)*1.D-5
      RCUBX = RCU*1.D-2
      TLX = 99.D+0
      THX = 101.D+0
      TKLX = TLX+TABS
      TKHX = THX+TABS
!
!---  Salt effects   ---
!
      XIX = CXIX(1)*TKX + CXIX(2)/TKX + CXIX(3)/(TKX**2)
      CLX = CLMBX(1)*TKX + CLMBX(2)/TKX + CLMBX(3)/(TKX**2)
!
!---  Convert NaCl concentration to molality   ---
!
      XMOLSX = 1.D+3*(XLSX/(1.D+0-XLSX))/WTMS
      XMOLWX = 1.D+3/WTMW
      XMLSX = XMOLSX/(XMOLSX+XMOLWX)
      XMOLNAX = XMOLSX
      XMOLCLX = XMOLSX
      APCX = (1.D+0+(XMOLNAX+XMOLCLX)/XMOLWX)*
     &  EXP(2.D+0*CLX*XMOLNAX + XIX*XMOLCLX*XMOLNAX)
!
!---  Simplified form of the Spycher equation, applicable to 275˚C,
!     but liquid CO2 equilibria are only approximate  ---
!
      IF( ISLC(32).EQ.2 ) THEN
!
!---    Water mole fraction in CO2-rich phase, using a double
!       exponential equation form  ---
!
        Y0X = CXMGW(1) + CXMGW(2)*(TX**CXMGW(3))
        A1X = CXMGW(4) + CXMGW(5)/(1.D+0 + EXP(-(TX-CXMGW(6))/
     &    CXMGW(7)))
        TAU1X = CXMGW(8) + CXMGW(9)*(TX**CXMGW(10))
        A2X = CXMGW(11) + CXMGW(12)/(1.D+0 + EXP(-(TX-CXMGW(13))/
     &    CXMGW(14)))
        TAU2X = CXMGW(15) + CXMGW(16)*(TX**CXMGW(17))
        XMGWX = Y0X + A1X*EXP(-TAU1X*PBX) + A2X*EXP(-TAU2X*PBX)
        XMGWX = MAX( MIN( XMGWX,1.D+0 ),0.D+0 )
!
!
!---    CO2 mole fraction in aqueous phase, temperature dependent
!       Henry's constant  ---
!
        HCAX = 6.305D-4*EXP(2.4D+3*((1.D+0/TKX)-(1.D+0/298.15)))
        PVAX = MAX( PX-PVBX,0.D+0 )/1.D+5
        XMLAX = PVAX*HCAX
!        Y0X = CXMLA(1) + (CXMLA(2)-CXMLA(1))/
!     &    (1.D+0 + (CXMLA(4)/TX)**CXMLA(3))
!        TAU1X = CXMLA(5) + CXMLA(6)*TX + CXMLA(7)*(TX**2) + 
!     &    CXMLA(8)*(TX**3)
!        XMLAX = Y0X*(1.D+0 - EXP(-TAU1X*PVAX))
        XMLAX = MAX( MIN( XMLAX,1.D+0 ),0.D+0 )
        XMLSX = (1.D+0-XMLAX)*XMLSX
        XMLWX = 1.D+0 - XMLAX - XMLSX
!
!---    Water vapor pressure lowering, reduce water mole fraction in 
!       the CO2-rich phase by the ratio of the reduced to 
!       saturated water vapor pressure  ---
!
        XMGWX = XMGWX*PVBX/PSBX
        XMGAX = 1.D+0-XMGWX
!
!---  Low temperature regime, direct solution  ---
!
      ELSEIF( TX.LT.TLX .OR. ISLC(32).EQ.3 ) THEN
        AX = CACX(1) + CACX(2)*TKX
        BX = CBCX
!
!---    Polynomial coefficients for Redlich-Kwong equation of state  ---
!
        CAX = 1.D+0
        CBX = -RCUBX*TKX/PBX
        CCX = -((RCUBX*TKX*BX/PBX) - (AX/(PBX*SQRT(TKX))) + (BX**2))
        CDX = -(AX*BX/(PBX*SQRT(TKX)))
        CALL NICKALLS( CAX,CBX,CCX,CDX,V1X,V2X,V3X )
!
!---    Molar volume, cm^3/mol  ---
!
        VGX = MAX( V1X,V2X,V3X )
        VNX = MIN( V1X,V2X,V3X )
!
!---    Gas conditions  ---
!
        IF( TKX.GE.TCRA ) THEN
          VX = VGX
        ELSE
          W1X = PBX*(VGX-VNX)
          W2X = RCUBX*TKX*LOG((VGX-BX)/(VNX-BX)) + (AX/(SQRT(TKX)*BX))*
     &      LOG(((VGX+BX)*VNX)/(VNX+BX)*VGX)
!
!---      Gas phase condition  ---
!
          IF( (W2X-W1X)/EPSL.GE.EPSL ) THEN
            VX = VGX
!
!---      Liquid phase condition  ---
!
          ELSEIF( (W1X-W2X)/EPSL.GE.EPSL ) THEN
            VX = VNX
!
!---      Gas-liquid phase condition, use gas molar volume  ---
!
          ELSE
            VX = VGX
          ENDIF
        ENDIF
!
!---    Fugacity coefficients  ---
!
        AIX(1,1) = CAWX(1)
        AIX(2,1) = CACWX(1)
        AIX(1,2) = CACWX(1)
        AIX(2,2) = CACX(1) + CACX(2)*TKX
        BIX(1) = CBWX
        BIX(2) = CBCX
        YX(1) = 0.D+0
        YX(2) = 1.D+0
        DO 120 K = 1,2
          SUMX = 0.D+0
          DO 110 I = 1,2
            SUMX = SUMX + YX(I)*AIX(I,K)
  110     CONTINUE
          FUGX(K) = LOG(VX/(VX-BX)) + (BIX(K)/(VX-BX))
          FUGX(K) = FUGX(K) - (2.D+0*SUMX/(RCUBX*(TKX**1.5D+0)*BX))*
     &      LOG((VX+BX)/VX) + (AX*BIX(K)/(RCUBX*(TKX**1.5D+0)*(BX**2)))*
     &      (LOG((VX+BX)/VX)-(BX/(VX+BX))) - LOG(PBX*VX/(RCUBX*TKX))
          FUGX(K) = EXP(FUGX(K))
  120   CONTINUE
!
!---    A parameter  ---
!
        EQKWOX = 0.D+0
        PREFX = CPRX
        DO 130 I = 1,5
          EQKWOX = EQKWOX + DLKWX(I)*(TX**(I-1))
 130    CONTINUE
        EQKWOX = 1.D+1**EQKWOX
        EQKWX = EQKWOX*EXP((PBX-PREFX)*CVWX/(RCUBX*TKX))
        CAX = (EQKWX/(FUGX(1)*PBX))
!
!---    B parameter  ---
!
        EQKCOX = 0.D+0
        PREFX = CPRX
        IF( TX.LT.31.D+0 .AND. VX.LT.94.D+0 ) THEN
          DO 140 I = 1,4
            EQKCOX = EQKCOX + CLKCNX(I)*(TX**(I-1))
 140      CONTINUE
        ELSE
          DO 150 I = 1,4
            EQKCOX = EQKCOX + CLKCGX(I)*(TX**(I-1))
 150      CONTINUE
        ENDIF
        EQKCOX = 1.D+1**EQKCOX
        EQKCX = EQKCOX*EXP((PBX-PREFX)*CVCX/(RCUBX*TKX))
        CBX = (FUGX(2)*PBX)/(XMOLWX*APCX*EQKCX)
!
!---    Equilibrium mole fractions  ---
!
        XMGWX = (1.D+0-CBX)*XMOLWX/(((1.D+0/CAX)-CBX)*
     &    (XMOLNAX + XMOLCLX + XMOLWX) + (XMOLNAX + XMOLCLX)*CBX)
        XMLAX = CBX*(1.D+0-XMGWX)*FSFLA
        XMLSX = (1.D+0-XMLAX)*XMLSX
        XMLWX = 1.D+0 - XMLAX - XMLSX
!
!---    Water vapor pressure lowering, reduce water mole fraction in 
!       the CO2-rich phase by the ratio of the reduced to 
!       saturated water vapor pressure  ---
!
        XMGWX = XMGWX*PVBX/PSBX
        XMGAX = 1.D+0-XMGWX
!
!---  High temperature regime, iterative solution  ---
!
      ELSEIF( TX.GT.THX ) THEN
!
!---    Initial guess of water gas mole fraction,
!       and CO2 aqueous mole fraction  ---
!
        CALL SP_W( TX,PSWX )
        YX(1) = PSWX/PX
        YX(2) = MAX( MIN( 1.D+0-YX(1),1.D+0 ),0.D+0 )
        XX(2) = 9.D-3
        XX(1) = 1.D+0-XX(2)
        Y1X = YX(1)
        Y2X = YX(2)
        X1X = XX(1)
        X2X = XX(2)
        CKX(1,1) = 0.D+0
        CKX(2,1) = DKCWX(1) + DKCWX(2)*TKX
        CKX(1,2) = DKWCX(1) + DKWCX(2)*TKX
        CKX(2,2) = 0.D+0
        SAX(1) = DAWX(1) + DAWX(2)*TKX
        SAX(2) = DACX(1) + DACX(2)*TKX
        BIX(1) = DBWX
        BIX(2) = DBCX
        EQKWOX = 0.D+0
        EQKCOX = 0.D+0
        PREFX = 0.D+0
        DTKX = MAX( (TKX-373.15D+0),0.D+0 )
        AMX = DAMX(1)*DTKX + DAMX(2)*(DTKX**2)
        VMWX = DVMWX(1) + DVMWX(2)*DTKX
        VMCX = DVMCX(1) + DVMCX(2)*DTKX      
!
!---    Equilibrium constants  ---
!
        DO 200 I = 1,5
          EQKWOX = EQKWOX + DLKWX(I)*(TX**(I-1))
          EQKCOX = EQKCOX + DLKCX(I)*(TX**(I-1))
          PREFX = PREFX + DPRX(I)*(TX**(I-1))          
  200   CONTINUE
        EQKWOX = 1.D+1**EQKWOX
        EQKCOX = 1.D+1**EQKCOX
        EQKWX = EQKWOX*EXP((PBX-PREFX)*VMWX/(RCUBX*TKX))
        EQKCX = EQKCOX*EXP((PBX-PREFX)*VMCX/(RCUBX*TKX))
!
!---    Newton-Raphson iteration  ---
!
        NC = 0
  210   CONTINUE
        NC = NC + 1
        IF( NC.GT.32 ) THEN
          YX(1) = Y1X
          YX(2) = Y2X
          XX(1) = X1X
          XX(2) = X2X
          M_ERR(1) = 'Mutual Solubility Convergence Failure: '
          M_ERR(2) = ' High Temperature Regime: Temperature = '
          IF( N_DB.LT.0 ) THEN
            M_ERR(3) = ' at Boundary Surface: '
          ELSE
            M_ERR(3) = ' at Node: '
          ENDIF
          CALL PATH
          R_ERR = TX
          I_ERR(1) = ABS(N_DB)
          I_ERR(2) = 2
          I_ERR(3) = 3
          I_ERR(4) = ID
        ENDIF
!
!---    Loop over increments  ---
!
        DO 280 M = 1,3
          YZ(1) = YX(1)
          XZ(2) = XX(2)
          IF( M.EQ.2 ) YZ(1) = YX(1) + 1.D-9
          IF( M.EQ.3 ) XZ(2) = XX(2) + 1.D-9
          YZ(2) = MAX( MIN( 1.D+0-YZ(1),1.D+0 ),0.D+0 )
          XZ(1) = MAX( MIN( 1.D+0-XZ(2),1.D+0 ),0.D+0 )
!
!---      Compute interaction parameters  ---
!
          DO I = 1,2
          DO J = 1,2
            IF( I.EQ.J ) THEN
              SKX(I,J) = 0.D+0
            ELSE
              SKX(I,J) = CKX(I,J)*YZ(I) + CKX(J,I)*YZ(J)
            ENDIF
          ENDDO
          ENDDO
!
!---      a-ij parameters  ---
!
          DO I = 1,2
          DO J = 1,2
            AIX(I,J) = SQRT(SAX(I)*SAX(J))*(1.D+0-SKX(I,J))
          ENDDO
          ENDDO
!
!---      a-mix and b-mix parameters  ---
!
          AX = 0.D+0
          BX = 0.D+0
          DO 250 I = 1,2
            DO 240 J = 1,2
              AX = AX + YZ(I)*YZ(J)*AIX(I,J)
  240       CONTINUE
            BX = BX + YZ(I)*BIX(I) 
  250     CONTINUE   
!
!---      Polynomial coefficients for Redlich-Kwong equation of state  ---
!
          CAX = 1.D+0
          CBX = -RCUBX*TKX/PBX
          CCX = -((RCUBX*TKX*BX/PBX) - (AX/(PBX*SQRT(TKX))) + (BX**2))
          CDX = -(AX*BX/(PBX*SQRT(TKX)))
          CALL NICKALLS( CAX,CBX,CCX,CDX,V1X,V2X,V3X )
!
!---      Gas (mixture) molar volume, cm^3/mol  ---
!
          VX = MAX( V1X,V2X,V3X )
!
!---      Fugacity coefficients  ---
!
          DO 270 K = 1,2
            FUGX(K) = (BIX(K)/BX)*((PBX*VX/(RCUBX*TKX))-1.D+0) -
     &        LOG(PBX*(VX-BX)/(RCUBX*TKX))
            VAR1X = 0.D+0
            DO 260 I = 1,2
              VAR1X = VAR1X + YZ(I)*(AIX(I,K)+AIX(K,I))
  260       CONTINUE
            FUGX(K) = FUGX(K) + (VAR1X/AX - BIX(K)/BX)*
     &        (AX/(BX*RCUBX*(TKX**1.5D+0)))*LOG(VX/(VX+BX))
            FUGX(K) = EXP(FUGX(K))
  270     CONTINUE
!
!---      Activities  ---
!
          ACTVX(1) = EXP( (AMX-2.D+0*AMX*XZ(1))*(XZ(2)**2) )
          ACTVX(2) = EXP( 2.D+0*AMX*XZ(2)*(XZ(1)**2) )
!
!---      A parameter  ---
!
          CAX = (EQKWX*ACTVX(1)/(FUGX(1)*PBX))
!
!---      B parameter  ---
!
          CBX = (FUGX(2)*PBX)/(XMOLWX*ACTVX(2)*APCX*EQKCX)
!
!---      Update mole fractions ---
!
          GX(1,M) = YZ(1) - (1.D+0-CBX)/((1.D+0/CAX)-CBX)
          GX(2,M) = XZ(2) - CBX*(1.D+0-YZ(1))
  280   CONTINUE
!
!---    Load solution vector and Jacobian matrix for
!       two-phase conditions  ---
!
        DO 290 M = 1,2
          AJ(M,1) = (GX(M,2)-GX(M,1))/1.D-9
          AJ(M,2) = (GX(M,3)-GX(M,1))/1.D-9
          BJ(M) = -GX(M,1)
  290   CONTINUE
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
        DYX = SIGN(MIN(ABS(BJ(1)),5.D-2),BJ(1))
        DXX = SIGN(MIN(ABS(BJ(2)),5.D-2),BJ(2))
        YX(1) = MAX( MIN( YX(1)+DYX,1.D+0 ),0.D+0 )
        XX(2) = MAX( MIN( XX(2)+DXX,1.D+0 ),0.D+0 )
!
!---    Check for convergence ---
!
        ERR1X = ABS((DYX))
        ERR2X = ABS((DXX))
        IF( ERR1X.GT.1.D-6 .OR. ERR2X.GT.1.D-6 ) GOTO 210  
        YX(2) = MAX( MIN( 1.D+0-YX(1),1.D+0 ),0.D+0 )
        XX(1) = MAX( MIN( 1.D+0-XX(2),1.D+0 ),0.D+0 )
!
!---    Compute interaction parameters  ---
!
        DO I = 1,2
        DO J = 1,2
          IF( I.EQ.J ) THEN
            SKX(I,J) = 0.D+0
          ELSE
            SKX(I,J) = CKX(I,J)*YX(I) + CKX(J,I)*YX(J)
          ENDIF
        ENDDO
        ENDDO
!
!---    a-ij parameters  ---
!
        DO I = 1,2
        DO J = 1,2
          AIX(I,J) = SQRT(SAX(I)*SAX(J))*(1.D+0-SKX(I,J))
        ENDDO
        ENDDO
!
!---    a-mix and b-mix parameters  ---
!
        AX = 0.D+0
        BX = 0.D+0
        DO 295 I = 1,2
          DO 294 J = 1,2
            AX = AX + YX(I)*YX(J)*AIX(I,J)
  294     CONTINUE
          BX = BX + YX(I)*BIX(I) 
  295   CONTINUE   
!
!---    Polynomial coefficients for Redlich-Kwong equation of state  ---
!
        CAX = 1.D+0
        CBX = -RCUBX*TKX/PBX
        CCX = -((RCUBX*TKX*BX/PBX) - (AX/(PBX*SQRT(TKX))) + (BX**2))
        CDX = -(AX*BX/(PBX*SQRT(TKX)))
        CALL NICKALLS( CAX,CBX,CCX,CDX,V1X,V2X,V3X )
!
!---    Gas (mixture) molar volume, cm^3/mol  ---
!
        VX = MAX( V1X,V2X,V3X )
!
!---    Fugacity coefficients  ---
!
        DO 297 K = 1,2
          FUGX(K) = (BIX(K)/BX)*((PBX*VX/(RCUBX*TKX))-1.D+0) -
     &      LOG(PBX*(VX-BX)/(RCUBX*TKX))
          VAR1X = 0.D+0
          DO 296 I = 1,2
            VAR1X = VAR1X + YX(I)*(AIX(I,K)+AIX(K,I))
  296     CONTINUE
          FUGX(K) = FUGX(K) + (VAR1X/AX - BIX(K)/BX)*
     &      (AX/(BX*RCUBX*(TKX**1.5D+0)))*LOG(VX/(VX+BX))
          FUGX(K) = EXP(FUGX(K))
  297   CONTINUE
!
!---    Activities  ---
!
        ACTVX(1) = EXP( (AMX-2.D+0*AMX*XX(1))*(XX(2)**2) )
        ACTVX(2) = EXP( 2.D+0*AMX*XX(2)*(XX(1)**2) )
!
!---    A parameter  ---
!
        CAX = (EQKWX*ACTVX(1)/(FUGX(1)*PBX))
!
!---    B parameter  ---
!
        CBX = (FUGX(2)*PBX)/(XMOLWX*ACTVX(2)*APCX*EQKCX)
!
!---    Equilibrium mole fractions  ---
!
        XMGWX = (1.D+0-CBX)*XMOLWX/(((1.D+0/CAX)-CBX)*
     &    (XMOLNAX + XMOLCLX + XMOLWX) + (XMOLNAX + XMOLCLX)*CBX)
        XMLAX = CBX*(1.D+0-XMGWX)*FSFLA
        XMLSX = (1.D+0-XMLAX)*XMLSX
        XMLWX = 1.D+0 - XMLAX - XMLSX
!
!---    Water vapor pressure lowering, reduce water mole fraction in 
!       the CO2-rich phase by the ratio of the reduced to 
!       saturated water vapor pressure  ---
!
        XMGWX = XMGWX*PVBX/PSBX
        XMGAX = 1.D+0-XMGWX
!
!---  Intermediate temperature regime ( 99.0 C < TX < 109.0 C ),
!     linear interpolation of low- and high-temperature regimes  ---
!
      ELSE
        AX = CACX(1) + CACX(2)*TKLX
        BX = CBCX
!
!---    Polynomial coefficients for Redlich-Kwong equation of state  ---
!
        CAX = 1.D+0
        CBX = -RCUBX*TKLX/PBX
        CCX = -((RCUBX*TKLX*BX/PBX) - (AX/(PBX*SQRT(TKLX))) + (BX**2))
        CDX = -(AX*BX/(PBX*SQRT(TKLX)))
        CALL NICKALLS( CAX,CBX,CCX,CDX,V1X,V2X,V3X )
!
!---    Molar volume, cm^3/mol  ---
!
        VX = MAX( V1X,V2X,V3X )
!
!---    Fugacity coefficients  ---
!
        AIX(1,1) = CAWX(1)
        AIX(2,1) = CACWX(1)
        AIX(1,2) = CACWX(1)
        AIX(2,2) = CACX(1) + CACX(2)*TKLX
        BIX(1) = CBWX
        BIX(2) = CBCX
        YX(1) = 0.D+0
        YX(2) = 1.D+0
        DO 320 K = 1,2
          SUMX = 0.D+0
          DO 310 I = 1,2
            SUMX = SUMX + YX(I)*AIX(I,K)
  310     CONTINUE
          FUGLX(K) = LOG(VX/(VX-BX)) + (BIX(K)/(VX-BX))
          FUGLX(K) = FUGLX(K) - (2.D+0*SUMX/(RCUBX*(TKLX**1.5D+0)*BX))*
     &      LOG((VX+BX)/VX)+(AX*BIX(K)/(RCUBX*(TKLX**1.5D+0)*(BX**2)))*
     &      (LOG((VX+BX)/VX)-(BX/(VX+BX))) - LOG(PBX*VX/(RCUBX*TKLX))
          FUGLX(K) = EXP(FUGLX(K))
  320   CONTINUE
!
!---    A parameter  ---
!
        EQKWOX = 0.D+0
        PREFX = CPRX
        DO 330 I = 1,5
          EQKWOX = EQKWOX + DLKWX(I)*(TLX**(I-1))
 330    CONTINUE
        EQKWOX = 1.D+1**EQKWOX
        EQKWLX = EQKWOX*EXP((PBX-PREFX)*CVWX/(RCUBX*TKLX))
!
!---    B parameter  ---
!
        EQKCOX = 0.D+0
        PREFX = CPRX
        DO 340 I = 1,4
          EQKCOX = EQKCOX + CLKCGX(I)*(TLX**(I-1))
 340    CONTINUE
        EQKCOX = 1.D+1**EQKCOX
        EQKCLX = EQKCOX*EXP((PBX-PREFX)*CVCX/(RCUBX*TKLX))
!
!---    Initial guess of water gas mole fraction,
!       and CO2 aqueous mole fraction  ---
!
        CALL SP_W( THX,PSWX )
        YX(1) = PSWX/PX
        YX(2) = MAX( MIN( 1.D+0-YX(1),1.D+0 ),0.D+0 )
        XX(2) = 9.D-3
        XX(1) = 1.D+0-XX(2)
        Y1X = YX(1)
        Y2X = YX(2)
        X1X = XX(1)
        X2X = XX(2)
        CKX(1,1) = 0.D+0
        CKX(2,1) = DKCWX(1) + DKCWX(2)*TKHX
        CKX(1,2) = DKWCX(1) + DKWCX(2)*TKHX
        CKX(2,2) = 0.D+0
        SAX(1) = DAWX(1) + DAWX(2)*TKHX
        SAX(2) = DACX(1) + DACX(2)*TKHX
        BIX(1) = DBWX
        BIX(2) = DBCX
        EQKWOX = 0.D+0
        EQKCOX = 0.D+0
        PREFX = 0.D+0
        DTKX = MAX( (TKHX-373.15D+0),0.D+0 )
        AMX = DAMX(1)*DTKX + DAMX(2)*(DTKX**2)
        VMWX = DVMWX(1) + DVMWX(2)*DTKX
        VMCX = DVMCX(1) + DVMCX(2)*DTKX      
!
!---    Equilibrium constants  ---
!
        DO 400 I = 1,5
          EQKWOX = EQKWOX + DLKWX(I)*(THX**(I-1))
          EQKCOX = EQKCOX + DLKCX(I)*(THX**(I-1))
          PREFX = PREFX + DPRX(I)*(THX**(I-1))          
  400   CONTINUE
        EQKWOX = 1.D+1**EQKWOX
        EQKCOX = 1.D+1**EQKCOX
        EQKWHX = EQKWOX*EXP((PBX-PREFX)*VMWX/(RCUBX*TKHX))
        EQKCHX = EQKCOX*EXP((PBX-PREFX)*VMCX/(RCUBX*TKHX))
!
!---    Sucessive-substitution iteration  ---
!
        NC = 0
  410   CONTINUE
        NC = NC + 1
        IF( NC.GT.32 ) THEN
          YX(1) = Y1X
          YX(2) = Y2X
          XX(1) = X1X
          XX(2) = X2X
          M_ERR(1) = 'Mutual Solubility Convergence Failure: '
          M_ERR(2) = ' Intermediate Temperature Regime: Temperature = '
          IF( N_DB.LT.0 ) THEN
            M_ERR(3) = ' at Boundary Surface: '
          ELSE
            M_ERR(3) = ' at Node: '
          ENDIF
          CALL PATH
          R_ERR = TX
          I_ERR(1) = ABS(N_DB)
          I_ERR(2) = 2
          I_ERR(3) = 3
          I_ERR(4) = ID
        ENDIF
        YXOX = YX(1)
!
!---    Compute interaction parameters  ---
!
        DO I = 1,2
        DO J = 1,2
          IF( I.EQ.J ) THEN
            SKX(I,J) = 0.D+0
          ELSE
            SKX(I,J) = CKX(I,J)*YX(I) + CKX(J,I)*YX(J)
          ENDIF
        ENDDO
        ENDDO
!
!---    a-ij parameters  ---
!
        DO I = 1,2
        DO J = 1,2
          AIX(I,J) = SQRT(SAX(I)*SAX(J))*(1.D+0-SKX(I,J))
        ENDDO
        ENDDO
!
!---    a-mix and b-mix parameters  ---
!
        AX = 0.D+0
        BX = 0.D+0
        DO 450 I = 1,2
          DO 440 J = 1,2
            AX = AX + YX(I)*YX(J)*AIX(I,J)
  440     CONTINUE
          BX = BX + YX(I)*BIX(I) 
  450   CONTINUE   
!
!---    Polynomial coefficients for Redlich-Kwong equation of state  ---
!
        CAX = 1.D+0
        CBX = -RCUBX*TKHX/PBX
        CCX = -((RCUBX*TKHX*BX/PBX) - (AX/(PBX*SQRT(TKHX))) + (BX**2))
        CDX = -(AX*BX/(PBX*SQRT(TKHX)))
        CALL NICKALLS( CAX,CBX,CCX,CDX,V1X,V2X,V3X )
!
!---    Gas molar volume, cm^3/mol  ---
!
        VX = MAX( V1X,V2X,V3X )
!
!---    Fugacity coefficients  ---
!
        DO 480 K = 1,2
          FUGHX(K) = (BIX(K)/BX)*((PBX*VX/(RCUBX*TKHX))-1.D+0) -
     &      LOG(PBX*(VX-BX)/(RCUBX*TKHX))
          VAR1X = 0.D+0
          DO 470 I = 1,2
            VAR1X = VAR1X + YX(I)*(AIX(I,K)+AIX(K,I))
  470     CONTINUE
          FUGHX(K) = FUGHX(K) + (VAR1X/AX - BIX(K)/BX)*
     &      (AX/(BX*RCUBX*(TKHX**1.5D+0)))*LOG(VX/(VX+BX))
          FUGHX(K) = EXP(FUGHX(K))
  480   CONTINUE
!
!---    Activities  ---
!
        ACTVX(1) = EXP( (AMX-2.D+0*AMX*XX(1))*(XX(2)**2) )
        ACTVX(2) = EXP( 2.D+0*AMX*XX(2)*(XX(1)**2) )
!
!---    Linear interpolation  ---
!
        DO 490 K = 1,2
          FUGX(K) = (TX-TLX)*(FUGHX(K)-FUGLX(K))/(THX-TLX) + FUGLX(K)
 490    CONTINUE
        EQKWX = (TX-TLX)*(EQKWHX-EQKWLX)/(THX-TLX) + EQKWLX
        EQKCX = (TX-TLX)*(EQKCHX-EQKCLX)/(THX-TLX) + EQKCLX
        ACTVX(1) = (TX-TLX)*(ACTVX(1)-1.D+0)/(THX-TLX) + 1.D+0
        ACTVX(2) = (TX-TLX)*(ACTVX(2)-1.D+0)/(THX-TLX) + 1.D+0
!
!---    A parameter  ---
!
        CAX = (EQKWX*ACTVX(1)/(FUGX(1)*PBX))
!
!---    B parameter  ---
!
        CBX = (FUGX(2)*PBX)/(XMOLWX*ACTVX(2)*APCX*EQKCX)
!
!---    Update mole fractions ---
!
        YX(1) = (1.D+0-CBX)/((1.D+0/CAX)-CBX)
        XX(2) = CBX*(1.D+0-YX(1))
        YX(2) = MAX( MIN( 1.D+0-YX(1),1.D+0 ),0.D+0 )
        XX(1) = MAX( MIN( 1.D+0-XX(2),1.D+0 ),0.D+0 )
!
!---    Check for convergence ---
!
        ERRX = ABS((YX(1)-YXOX)/(PSWX/PX))
        IF( ERRX.GT.1.D-6 ) GOTO 410  
!
!---    Equilibrium mole fractions  ---
!
        XMGWX = (1.D+0-CBX)*XMOLWX/(((1.D+0/CAX)-CBX)*
     &    (XMOLNAX + XMOLCLX + XMOLWX) + (XMOLNAX + XMOLCLX)*CBX)
        XMLAX = CBX*(1.D+0-XMGWX)*FSFLA
        XMLSX = (1.D+0-XMLAX)*XMLSX
        XMLWX = 1.D+0 - XMLAX - XMLSX
!
!---    Water vapor pressure lowering, reduce water mole fraction in 
!       the CO2-rich phase by the ratio of the reduced to 
!       saturated water vapor pressure  ---
!
        XMGWX = XMGWX*PVBX/PSBX
        XMGAX = 1.D+0-XMGWX
      ENDIF
!
!---  Limit mole fractions  ---
!
      IF( XMLWX.LT.1.D-16 ) XMLWX = 0.D+0
      IF( XMLAX.LT.1.D-16 ) XMLAX = 0.D+0
      IF( XMLSX.LT.1.D-16 ) XMLSX = 0.D+0
      IF( XMGWX.LT.1.D-16 ) XMGWX = 0.D+0
      IF( XMGAX.LT.1.D-16 ) XMGAX = 0.D+0
!
!---  Gas partial pressures, Pa  ---
!
      PGAX = XMGAX*PX
      PGWX = XMGWX*PX
!
!---  Equilibrium mass fractions  ---
!
      WTMGX = XMGWX*WTMW + XMGAX*WTMA
      XGWX = XMGWX*WTMW/WTMGX
      XGAX = XMGAX*WTMA/WTMGX
      WTMLX = XMLWX*WTMW + XMLAX*WTMA + XMLSX*WTMS
      XLWX = XMLWX*WTMW/WTMLX
      XLAX = XMLAX*WTMA/WTMLX
      XLSX = XMLSX*WTMS/WTMLX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of EQUIL group  ---
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION FSPLNY( FX,ITS,ITE )
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
!     Given arrays TBLX and TBLY containing a tabulated function
!     (i.e., y(n) = f(x(n))), with x(1) < x(2) < ... < x(n), and given
!     the array TBLDDY, which is output from subroutine SPLINE, and
!     given a value of FX this function returns a cubic spline
!     interpolated value.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, Battelle, PNL, September, 1994.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TABL
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
      CHARACTER(64) :: SUBLOGX
      CHARACTER(256) :: CHMSGX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FSPLNY'
!
!---  Ascending table order  ---
!
      IF( TBLX(ITE).GT.TBLX(ITS) ) THEN
        IF( FX.LE.TBLX(ITS) ) THEN
          FSPLNY = TBLY(ITS)
!
!---      Reset subroutine string sequence  ---
!
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ELSEIF( FX.GE.TBLX(ITE) ) THEN
          FSPLNY = TBLY(ITE)
!
!---      Reset subroutine string sequence  ---
!
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ENDIF
!
!---  Find the right place in the table by means of bisection  ---
!
        KLO = ITS
        KHI = ITE
   10   CONTINUE
        IF( KHI-KLO.GT.1 ) THEN
          K = (KHI+KLO)/2
          IF( TBLX(K).GT.FX ) THEN
            KHI = K
          ELSE
            KLO = K
          ENDIF
          GOTO 10
        ENDIF
        H = TBLX(KHI)-TBLX(KLO)
        IF( ABS(H)/EPSL.LT.EPSL ) THEN
          H = 1.D+0
          M_ERR(1) = 'Invalid Table in FSPLNY '
          CALL PATH
          I_ERR(2) = 0
          I_ERR(3) = 0
          I_ERR(4) = ID
        ENDIF
!
!---  Evaluate cubic spline  ---
!
        A = (TBLX(KHI)-FX)/H
        B = (FX-TBLX(KLO))/H
        FSPLNY = A*TBLY(KLO)+B*TBLY(KHI)+
     &    ((A**3-A)*TBLDDY(KLO)+(B**3-B)*TBLDDY(KHI))*(H**2)/6.D+0
!
!---  Descending table order  ---
!
      ELSEIF( TBLX(ITE).LT.TBLX(ITS) ) THEN
        IF( FX.GE.TBLX(ITS) ) THEN
          FSPLNY = TBLY(ITS)
!
!---      Reset subroutine string sequence  ---
!
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ELSEIF( FX.LE.TBLX(ITE) ) THEN
          FSPLNY = TBLY(ITE)
!
!---      Reset subroutine string sequence  ---
!
          ISUB_LOG = ISUB_LOG-1
          RETURN
        ENDIF
!
!---  Find the right place in the table by means of bisection  ---
!
        KLO = ITS
        KHI = ITE
   20   CONTINUE
        IF( KHI-KLO.GT.1 ) THEN
          K = (KHI+KLO)/2
          IF( TBLX(K).LT.FX ) THEN
            KHI = K
          ELSE
            KLO = K
          ENDIF
          GOTO 20
        ENDIF
        H = TBLX(KLO)-TBLX(KHI)
        IF( ABS(H)/EPSL.LT.EPSL ) THEN
          H = 1.D+0
          M_ERR(1) = 'Invalid Table in FSPLNY '
          CALL PATH
          I_ERR(2) = 0
          I_ERR(3) = 0
          I_ERR(4) = ID
        ENDIF
!
!---  Evaluate cubic spline  ---
!
        A = (TBLX(KLO)-FX)/H
        B = (FX-TBLX(KHI))/H
        FSPLNY = A*TBLY(KHI)+B*TBLY(KLO)+
     &    ((A**3-A)*TBLDDY(KHI)+(B**3-B)*TBLDDY(KLO))*(H**2)/6.D+0
      ELSE
        M_ERR(1) = 'Invalid Table in FSPLNY '
        CALL PATH
        I_ERR(2) = 0
        I_ERR(3) = 0
        I_ERR(4) = ID
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FSPLNY group
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION FNTBLX( FY,N1,N2,INDX )
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
!     Linearly interpolate from a table of values.
!
!     INDX = 0 : Table truncation beyond table limits.
!     INDX = 1 : Table extrapolation beyond table limits.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, Battelle, PNL, March, 1993.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TABL
      USE GRID
      USE SOLTN
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
      SUB_LOG(ISUB_LOG) = '/FNTBLX'
!
!---  Ascending table order  ---
!
      IF( TBLY(N2).GT.TBLY(N1) ) THEN
        IF( FY.LE.TBLY(N1) ) THEN
          IF( INDX.EQ.0 ) THEN
            FNTBLX = TBLX(N1)
          ELSE
            FNTBLX = (FY-TBLY(N1))*(TBLX(N1+1)-TBLX(N1))/
     &        (TBLY(N1+1)-TBLY(N1)) + TBLX(N1)
          ENDIF
        ELSEIF( FY.GE.TBLY(N2) ) THEN
          IF( INDX.EQ.0 ) THEN
            FNTBLX = TBLX(N2)
          ELSE
            FNTBLX = (FY-TBLY(N2))*(TBLX(N2)-TBLX(N2-1))/
     &        (TBLY(N2)-TBLY(N2-1)) + TBLX(N2)
          ENDIF
        ELSE
          DO N = N1+1,N2
            IF( FY.LE.TBLY(N) ) THEN
              FNTBLX = (FY-TBLY(N-1))*(TBLX(N)-TBLX(N-1))/
     &          (TBLY(N)-TBLY(N-1)) + TBLX(N-1)
              EXIT
            ENDIF
          ENDDO
        ENDIF
!
!---  Descending table order  ---
!
      ELSEIF( TBLY(N2) .LT. TBLY(N1) ) THEN
        IF( FY.LE.TBLY(N2) ) THEN
          IF( INDX.EQ.0 ) THEN
            FNTBLX = TBLX(N2)
          ELSE
            FNTBLX = (FY-TBLY(N2))*(TBLX(N2-1)-TBLX(N2))/
     &        (TBLY(N2-1)-TBLY(N2)) + TBLX(N2)
          ENDIF
        ELSEIF( FY.GE.TBLY(N1) ) THEN
          IF( INDX.EQ.0 ) THEN
            FNTBLX = TBLX(N1)
          ELSE
            FNTBLX = (FY-TBLY(N1))*(TBLX(N1)-TBLX(N1+1))/
     &        (TBLY(N1)-TBLY(N1+1)) + TBLX(N1)
          ENDIF
        ELSE
          DO N = N1+1,N2
            IF( FY.GE.TBLY(N) ) THEN
              FNTBLX = (FY-TBLY(N))*(TBLX(N-1)-TBLX(N))/
     &          (TBLY(N-1)-TBLY(N)) + TBLX(N)
              EXIT
            ENDIF
          ENDDO
        ENDIF
      ELSE
        A = 5.D-1
        B = 5.D-1
        FNTBLX = A*TBLX(N1) + B*TBLX(N2)
        M_ERR(1) = 'Invalid Table in FNTBLX '
        CALL PATH
        I_ERR(2) = 0
        I_ERR(3) = 0
        I_ERR(4) = ID
        ICNV = 4
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FNTBLX group
!
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION FNTBLY( FX,N1,N2,INDX )
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
!
!----------------------Description-------------------------------------!
!
!     Linearly interpolate from a table of values.
!
!     INDX = 0 : Table truncation beyond table limits.
!     INDX = 1 : Table extrapolation beyond table limits.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, Battelle, PNL, March, 1993.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE TABL
      USE GRID
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
      SUB_LOG(ISUB_LOG) = '/FNTBLY'
!
!---  Ascending table order  ---
!
      IF( TBLX(N2).GT.TBLX(N1) ) THEN
        IF( FX .LE. TBLX(N1) ) THEN
          IF( INDX.EQ.0 ) THEN
            FNTBLY = TBLY(N1)
          ELSE
            FNTBLY = (FX-TBLX(N1))*(TBLY(N1+1)-TBLY(N1))/
     &        (TBLX(N1+1)-TBLX(N1)) + TBLY(N1)
          ENDIF
        ELSEIF( FX .GE. TBLX(N2) ) THEN
          IF( INDX.EQ.0 ) THEN
            FNTBLY = TBLY(N2)
          ELSE
            FNTBLY = (FX-TBLX(N2))*(TBLY(N2)-TBLY(N2-1))/
     &        (TBLX(N2)-TBLX(N2-1)) + TBLY(N2)
          ENDIF
        ELSE
          DO N = N1+1,N2
            IF( FX .LE. TBLX(N) ) THEN
              FNTBLY = (FX-TBLX(N-1))*(TBLY(N)-TBLY(N-1))/
     &          (TBLX(N)-TBLX(N-1)) + TBLY(N-1)
              EXIT
            ENDIF
          ENDDO
        ENDIF
!
!---  Descending table order  ---
!
      ELSEIF( TBLX(N2).LT.TBLX(N1) ) THEN
        IF( FX .LE. TBLX(N2) ) THEN
          IF( INDX.EQ.0 ) THEN
            FNTBLY = TBLY(N2)
          ELSE
            FNTBLY = (FX-TBLX(N2))*(TBLY(N2-1)-TBLY(N2))/
     &        (TBLX(N2-1)-TBLX(N2)) + TBLY(N2)
          ENDIF
        ELSEIF( FX .GE. TBLX(N1) ) THEN
          IF( INDX.EQ.0 ) THEN
            FNTBLY = TBLY(N1)
          ELSE
            FNTBLY = (FX-TBLX(N1))*(TBLY(N1)-TBLY(N1+1))/
     &        (TBLX(N1)-TBLX(N1+1)) + TBLY(N1)
          ENDIF
        ELSE
          DO  N = N1+1,N2
            IF( FX .GE. TBLX(N) ) THEN
              FNTBLY = (FX-TBLX(N))*(TBLY(N-1)-TBLY(N))/
     &          (TBLX(N-1)-TBLX(N)) + TBLY(N)
              EXIT
            ENDIF
          ENDDO
        ENDIF
      ELSE
        A = 5.D-1
        B = 5.D-1
        FNTBLY = A*TBLY(N1) + B*TBLY(N2)
        M_ERR(1) = 'Invalid Table in FNTBLY '
        CALL PATH
        I_ERR(2) = 0
        I_ERR(3) = 0
        I_ERR(4) = ID
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FNTBLY group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE HC_LA( TX,XLSX,HCX )
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
!     Battistelli, A., C. Claudio, and K. Pruess.  1997.  The simulator
!     TOUGH2/EWASG for modelling geothermal reservoirs with brines and
!     noncondensible gas.  Geothermics, 26(4): 437-464.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 15 March 2002.
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
      REAL*8 CBX(6),CCX(5)
!
!----------------------Data Statements---------------------------------!
!
!      SAVE CBX,CCX
      DATA CBX / 7.83666D+7, 1.96025D+6, 8.20574D+4, -7.40674D+2,
     &  2.18380D+0, -2.20999D-3 /
      DATA CCX / 1.19784D-1, -7.17823D-4, 4.93854D-6, -1.03826D-8,
     &  1.08233D-11 /
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
        HCX = HCX + CBX(I+1)*(TX**I)
   10 CONTINUE
!
!---  Empirical formulation by Battistelli for salting-out
!     coefficient for NaCl aqueous solutions as a function
!     of temperature  ---
!
      SKBX = 0.D+0
      DO 20 I = 0,4
        SKBX = SKBX + CCX(I+1)*(TX**I)
   20 CONTINUE
!
!---  Empirical formulation by Battistelli for Henry's coefficient
!     coefficient for NaCl aqueous solutions as a function
!     of temperature and salt molality  ---
!
      GLSX = 1.D+3*XLSX/(WTMS*(1.D+0-XLSX))
      HCX = HCX*(1.D+1**(GLSX*SKBX))
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
      SUBROUTINE ITL_A( PX,P_TPX,P_CRX,TX,I_PX,I_SX,I_TX,I_VX )
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
!     computing CO2 properties.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 22 April 2002.
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
      INTEGER I_PX(2),I_TX(2),I_VX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ITL_A'
!
!---  Convert pressure to MPa and temperature to degrees Kelvin  ---
!
      PMX = 1.D-6*PX
      PM_TPX = 1.D-6*P_TPX
      PM_CRX = 1.D-6*P_CRX
      TKX = TX + TABS
      IPTPX = MIN( MAX( IPTP(INCG)-1,1 ),IP_TA(INCG) )
      IPCRX = MIN( MAX( IPCR(INCG)+1,1 ),IP_TA(INCG) )
!
!---  Find pressure indices  ---
!
      IF( PMX.LT.P_TA(1,INCG) ) PMX = P_TA(1,INCG)
      IF( PMX.GT.P_TA(IP_TA(INCG),INCG) ) PMX = P_TA(IP_TA(INCG),INCG)
      IPLX = 1
      IPUX = IP_TA(INCG)
   10 IF( IPUX-IPLX.GT.1 ) THEN
        IPM = (IPLX+IPUX)/2
        IF( (P_TA(IP_TA(INCG),INCG).GT.P_TA(1,INCG)).EQV.
     &    (PMX.GT.P_TA(IPM,INCG)) ) THEN
          IPLX = IPM
        ELSE
          IPUX = IPM
        ENDIF
        GOTO 10
      ENDIF
!
!---  Pressure below the triple-point pressure,
!     search the across the table  ---
!
      IF( PMX.LT.P_TA(IPTPX,INCG) ) THEN
        ISLX = 1
        IELX = IT_TA(IPLX,INCG)
        ISUX = 1
        IEUX = IT_TA(IPUX,INCG)
        I_VX = 0
!
!---  Pressure above the critical-point pressure,
!     search the across the table  ---
!
      ELSEIF( PMX.GT.P_TA(IPCRX,INCG) ) THEN
        ISLX = 1
        IELX = IT_TA(IPLX,INCG)
        ISUX = 1
        IEUX = IT_TA(IPUX,INCG)
        I_VX = 2
!
!---  Pressure greater than or equal to the triple-point pressure
!     and pressure less than or equal to the 
!     critical-point pressure  ---
!
      ELSE
!
!---    Determine saturation temperature at the specified pressure,
!       via linear interpolation  ---
!
        PLX = P_TA(IPLX,INCG)
        PUX = P_TA(IPLX+1,INCG)
        TLX = T_TA(IV_TA(IPLX,INCG),IPLX,INCG)
        TUX = T_TA(IV_TA(IPUX,INCG),IPUX,INCG)
        TSX = ((PMX-PLX)/(PUX-PLX))*(TUX-TLX) + TLX
        IF( PMX.GE.PM_CRX ) TSX = TCRA
!
!---    Interpolation on gas side, searching from
!       the saturation line  ---
!
        IF( TKX.GE.TSX ) THEN
          ISLX = IV_TA(IPLX,INCG)
          IELX = IT_TA(IPLX,INCG)
          ISUX = IV_TA(IPUX,INCG)
          IEUX = IT_TA(IPUX,INCG)
          I_VX = 0
!
!---    Interpolation on liquid side, searching to
!       the saturation line  ---
!
        ELSE
          IPLX = IPTPX
          IPUX = IPCRX
          ISLX = 1
          IELX = IT_TA(IPLX,INCG)
          ISUX = 1
          IEUX = IT_TA(IPUX,INCG)
          I_VX = 1
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
      IF( TKX.LT.T_TA(ISLX,IPLX,INCG) ) THEN
        I_TX(1) = ISLX
      ELSEIF( TKX.GT.T_TA(IELX,IPLX,INCG) ) THEN
        I_TX(1) = IELX-1
      ELSE
        ITLX = 1
        ITUX = IT_TA(IPLX,INCG)
   20   IF( ITUX-ITLX.GT.1 ) THEN
          ITM = (ITLX+ITUX)/2
          IF( (T_TA(IELX,IPLX,INCG).GT.T_TA(ISLX,IPLX,INCG))
     &      .EQV.(TKX.GT.T_TA(ITM,IPLX,INCG)) ) THEN
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
      IF( TKX.LT.T_TA(ISUX,IPUX,INCG) ) THEN
        I_TX(2) = ISUX
      ELSEIF( TKX.GT.T_TA(IEUX,IPUX,INCG) ) THEN
        I_TX(2) = IEUX-1
      ELSE
        ITLX = 1
        ITUX = IT_TA(IPUX,INCG)
   30   IF( ITUX-ITLX.GT.1 ) THEN
          ITM = (ITLX+ITUX)/2
          IF( (T_TA(IEUX,IPUX,INCG).GT.T_TA(ISUX,IPUX,INCG))
     &      .EQV.(TKX.GT.T_TA(ITM,IPUX,INCG)) ) THEN
            ITLX = ITM
          ELSE
            ITUX = ITM
          ENDIF
          GOTO 30
        ENDIF
        I_TX(2) = ITLX
      ENDIF
!
!---  Find saturation line index  ---
!
      IPLX = MIN( MAX( IPTP(INCG),1 ),IP_TA(INCG) )
      IPUX = MIN( MAX( IPCR(INCG),1 ),IP_TA(INCG) )
      TLX = T_TA(IV_TA(IPLX,INCG),IPLX,INCG)
      TUX = T_TA(IV_TA(IPUX,INCG),IPUX,INCG)
      IF( TKX.LT.TLX ) THEN
        I_SX = 0
      ELSEIF( TKX.GT.TUX ) THEN
        I_SX = 0
      ELSE
   40   IF( IPUX-IPLX.GT.1 ) THEN
          IPMX = (IPLX+IPUX)/2
          TMX = T_TA(IV_TA(IPMX,INCG),IPMX,INCG)
          IF( (TUX.GT.TLX).EQV.(TKX.GT.TMX) ) THEN
            IPLX = IPMX
          ELSE
            IPUX = IPMX
          ENDIF
          GOTO 40
        ENDIF
        I_SX = IPLX
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
      SUBROUTINE LUDCMP( AX,NX,NPX,IX,DX )
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
!            Copyright Battelle Memorial Institute, 1996.
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Numerical Recipes, The Art of Scientific Computing
!     W.H. Press, B.P. Flannery, Saul A. Teukolsky, and W.T. Vetterling
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, August 1, 2000.
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
      REAL*8 AX(NPX,NPX),VV(NPX)
      INTEGER IX(NPX)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/LUDCMP'
      DX = 1.D+0
      DO I = 1,NX
        AAMAX = 0.D+0
        DO J = 1,NX
          IF( ABS(AX(I,J)).GT.AAMAX ) AAMAX = ABS(AX(I,J))
        ENDDO
        IF( ABS(AAMAX)/EPSL.LT.EPSL ) THEN
          M_ERR(1) = 'LUDCMP: Singular Matrix'
          CALL PATH
          I_ERR(1) = 0
          I_ERR(2) = 0
          I_ERR(3) = 0
          I_ERR(4) = ID
          DO J = 1,NX
            AX(I,J) = 0.D+0
            IF( I.EQ.J ) AX(I,J) = 1.D+0
          ENDDO
        ENDIF
        VV(I) = 1.D+0/AAMAX
      ENDDO
      IMAX = 0
      DO J = 1,NX
        DO I = 1,J-1
          SUM = AX(I,J)
          DO K = 1,I-1
            SUM = SUM - AX(I,K)*AX(K,J)
          ENDDO
          AX(I,J) = SUM
        ENDDO
        AAMAX = 0.D+0
        DO I = J,NX
          SUM = AX(I,J)
          DO K = 1,J-1
            SUM = SUM - AX(I,K)*AX(K,J)
          ENDDO
          AX(I,J) = SUM
          DUM = VV(I)*ABS(SUM)
          IF( DUM.GE.AAMAX ) THEN
            IMAX = I
            AAMAX = DUM
          ENDIF
        ENDDO
        IF( J.NE.IMAX ) THEN
          DO K = 1,NX
            DUM = AX(IMAX,K)
            AX(IMAX,K) = AX(J,K)
            AX(J,K) = DUM
          ENDDO
          DX = -DX
          VV(IMAX) = VV(J)
        ENDIF
        IX(J) = IMAX
        IF( ABS(AX(J,J))/EPSL.LT.EPSL ) AX(J,J) = 1.D-30
        IF( J.NE.NX ) THEN
          DUM = 1.D+0/AX(J,J)
          DO I = J+1,NX
            AX(I,J) = AX(I,J)*DUM
          ENDDO
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LUDCMP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LUBKSB( AX,NX,NPX,IX,BX )
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
!            Copyright Battelle Memorial Institute, 1996.
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Numerical Recipes, The Art of Scientific Computing
!     W.H. Press, B.P. Flannery, Saul A. Teukolsky, and W.T. Vetterling
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, August 1, 2000.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
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
      REAL*8 AX(NPX,NPX),BX(NPX)
      INTEGER IX(NPX)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/LUBKSB'
      II = 0
      DO 12 I = 1,NX
        IL = IX(I)
        SUM = BX(IL)
        BX(IL) = BX(I)
        IF( II.NE.0 ) THEN
          DO 11 J = II,I-1
            SUM = SUM - AX(I,J)*BX(J)
   11     CONTINUE
        ELSEIF( ABS(SUM)/EPSL.GT.EPSL ) THEN
          II = I
        ENDIF
        BX(I) = SUM
   12 CONTINUE
      DO 14 I = NX,1,-1
        SUM = BX(I)
        IF( I.LT.NX ) THEN
          DO 13 J = I+1,NX
            SUM = SUM - AX(I,J)*BX(J)
   13     CONTINUE
        ENDIF
        BX(I) = SUM/AX(I,I)
   14 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LUBKSB group  ---
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
!     Nickalls, R.W.D. 1993. A new approach to solving the cubic:
!     Cardan's solution revealed.  Math. Gaz. 77:354-359.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 14 June 2010.
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/NICKALLS'
      INCG = 1
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
!---  End of NICKALLS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PERM_LG( PERMX,N )
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
!     Calculation of average intrinsic permeability around well,
!     using logarithmic averaging.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 26 July 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/PERM_LG'
!
!---  Average intrinsic permeability  ---
!
      PERMX = 0.D+0
      NPERMX = 0
      IF( (PERM(1,N)/EPSL).GT.EPSL ) THEN
        PERMX = PERMX + LOG(PERM(1,N))
        NPERMX = NPERMX+1
      ENDIF
      IF( (PERM(2,N)/EPSL).GT.EPSL ) THEN
        PERMX = PERMX + LOG(PERM(2,N))
        NPERMX = NPERMX+1
      ENDIF
      IF( (PERM(3,N)/EPSL).GT.EPSL ) THEN
        PERMX = PERMX + LOG(PERM(3,N))
        NPERMX = NPERMX+1
      ENDIF
      IF( NPERMX.GT.0 ) THEN
        REALX = REAL(NPERMX)
        PERMX = EXP( PERMX/REALX )
      ELSE
        PERMX = 0.D+0
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PERM_LG group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PERM_R( SSX,PERMRFX,PORDX,N )
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
!     Calculation of permeability reduction factor.
!
!     Verma, A., and K. Pruess.  1988.  Thermohydrological Conditions
!     and Silica Redistribution Near High-Level Nuclear Wastes
!     Emplaced in Saturated Geological Formations.  Journal of
!     Geophysical Research, 93(B2):1159-1173.
!
!     Pruess, K., and J. Garcia.  2002.  Multiphase flow dynamics
!     during CO2 disposal into saline aquifers.  Environmental Geology
!     http://link.springer.de/link/service/journals/00254/contents
!     /01/00498/paper/s00254-001-0498-3ch110.html
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 29 May 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PERM_R'
!
!---  Reduced porosity with salt precipitation  ---
!
      PORDX = MAX( PORDX*(1.D+0-SSX),PORDX*PERM(5,N),1.D-12 )
!
!---  Normalized porosity  ---
!
      PORD_NX = MAX( (1.D+0-SSX-PERM(5,N))/(1.D+0-PERM(5,N)),0.D+0 )
!
!---  Tube area ratio  ---
!
      OMEGAX = 1.D+0 + (1.D+0/PERM(4,N))/((1.D+0/PERM(5,N))-1.D+0)
!
!---  Permeability reduction factor  ---
!
      PERMRFX = (PORD_NX**2)*(1.D+0-PERM(4,N)+PERM(4,N)/(OMEGAX**2))
     &  /(1.D+0-PERM(4,N)+PERM(4,N)*
     &  ((PORD_NX/(PORD_NX+OMEGAX-1.D+0))**2))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PERM_R group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!      
      SUBROUTINE PERM_V( SSX,PERMRFX,PORDX,PORM1,PORM0,N )
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
!     Calculation of permeability reduction factor.
!
!     Verma, A., and K. Pruess.  1988.  Thermohydrological Conditions
!     and Silica Redistribution Near High-Level Nuclear Wastes
!     Emplaced in Saturated Geological Formations.  Journal of
!     Geophysical Research, 93(B2):1159-1173.
!
!     Pruess, K., and J. Garcia.  2002.  Multiphase flow dynamics
!     during CO2 disposal into saline aquifers.  Environmental Geology
!     http://link.springer.de/link/service/journals/00254/contents
!     /01/00498/paper/s00254-001-0498-3ch110.html
!
!----------------------Authors-----------------------------------------!
!
!     Written by DH Bacon, PNNL, 31 October 2012!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PERM_V'
!
!---  Reduced porosity with salt precipitation  ---
!
      PORDX = MAX( PORDX*(1.D+0-SSX),PORDX*PERM(5,N),1.D-12 )
!
!--- Precipitated mineral saturation
!
      SPM = 0.D+0
      IF( ISLC(43).EQ.1 ) THEN
        SPM = 1.D+0 - PORM1/PORM0
      ENDIF
!
!---  Normalized porosity  ---
!
      PORD_NX = MAX( (1.D+0-SSX-SPM-PERM(5,N))/(1.D+0-PERM(5,N)),
     &  0.D+0 )
!
!---  Tube area ratio  ---
!
      OMEGAX = 1.D+0 + (1.D+0/PERM(4,N))/((1.D+0/PERM(5,N))-1.D+0)
!
!---  Permeability reduction factor  ---
!
      PERMRFX = (PORD_NX**2)*(1.D+0-PERM(4,N)+PERM(4,N)/(OMEGAX**2))
     &  /(1.D+0-PERM(4,N)+PERM(4,N)*
     &  ((PORD_NX/(PORD_NX+OMEGAX-1.D+0))**2))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PERM_V group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PORSTY_CO2E( PX,PREFX,PORDX,PORTX,N )
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
!     Compute diffusive and total porosities.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, Battelle, PNNL, 28 October 2021
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE GRID
      USE GEO_MECH
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PORSTY_CO2E'
      DPX = PX-PREFX
      ISLC50X = ABS(ISLC(50))
!
!---  Geomechanics simulations  --
!
      IF( ISLC(50).NE.0 .AND. MOD(ISLC50X,10).NE.2 .AND. 
     &  MOD(ISLC50X,10).NE.4 ) THEN
!
!---    Davis-Davis porosity versus mean effective stress model  ---
!
        IF( IPROP_GM(1,N).EQ.1 ) THEN
!
!---      Diffusive and total mechanical porosity, at iterate level k  ---
!
          PORDMCHKX = (POR(1,N)-PROP_GM(5,N))*
     &      EXP(-PROP_GM(6,N)*SIGV_GM(2,N)) +  PROP_GM(5,N)
          PORTMCHKX = (POR(1,N)-PROP_GM(5,N))*
     &      EXP(-PROP_GM(6,N)*SIGV_GM(2,N)) +  PROP_GM(5,N)
          IDX = MIN( ID+1,2 )
!          IF( N.EQ.IDX ) PRINT *,'ND(N) = ',ND(N),
!     &      'POR(1,N) = ',POR(1,N),'PROP_GM(5,N) = ',PROP_GM(5,N),
!     &      'PROP_GM(6,N) = ',PROP_GM(6,N),
!     &      'SIGV_GM(2,N) = ',SIGV_GM(2,N),'ID = ',ID
!
!---      Pressure differential, at iterate level k+1  ---
!
          DPK1X = PX-P_GM(2,N)
!
!---      Diffusive and total flow porosity, at iterate level k+1  ---
!
          DPORX = (POR(1,N)-PROP_GM(5,N))*PROP_GM(6,N)*
     &      PROP_GM(3,N)*EXP(-PROP_GM(6,N)*SIGV_GM(2,N))
          PORDX = PORDMCHKX + DPORX*DPK1X
          PORTX = PORTMCHKX + DPORX*DPK1X
!
!---    Classical model  ---
!
        ELSE
!
!---      Drained bulk modulus  --
!
          BLKDRNX = PROP_GM(1,N)/(3.0D+0*(1.0D+0-2.0D+0*PROP_GM(2,N)))
!
!---      1/N  ---
!
          OONMODX = (PROP_GM(3,N)-POR(1,N))*(1.0D+0-PROP_GM(3,N))
     &      /BLKDRNX
!
!---      Volumetric strain differential, at iterate level k  ---
!
          DEVKX = EPSV_GM(2,N)-EPSV_CMP(N)
!
!---      Pressure differential, at iterate level k  ---
!
          DPKX = P_GM(2,N)-PCMP(N)
!
!---      Diffusive and total mechanical porosity, 
!         at iterate level k  ---
!
          PORDMCHKX = POR(1,N) - PROP_GM(3,N)*DEVKX + OONMODX*DPKX
          PORTMCHKX = POR(2,N) - PROP_GM(3,N)*DEVKX + OONMODX*DPKX
!
!---      Pressure differential, at iterate level k+1  ---
!
          DPK1X = PX-P_GM(2,N)
!
!---      Diffusive and total flow porosity, at iterate level k+1  ---
!
          PORDX = PORDMCHKX + ((PROP_GM(3,N)**2)/BLKDRNX+OONMODX)*DPK1X
          PORTX = PORTMCHKX + ((PROP_GM(3,N)**2)/BLKDRNX+OONMODX)*DPK1X
         ENDIF
!
!---  No geomechanics  --
!
      ELSE        
!
!---    Pore compressibility w/ fixed bulk volume  ---
!
        IF( ISLC(15).EQ.1 ) THEN
!
!---      Reactive transport porosity alteration  ---
!
          IF( ISLC(43).EQ.1 ) THEN
            PORTX = POR0(1,N)*EXP(DPX*CMP(1,N))
            PORDX = POR0(2,N)*EXP(DPX*CMP(1,N))
          ELSE
            PORTX = POR(1,N)*EXP(DPX*CMP(1,N))
            PORDX = POR(2,N)*EXP(DPX*CMP(1,N))
          ENDIF
!
!---    Bulk compressibility w/ variable bulk volume  ---
!
        ELSEIF( ISLC(15).EQ.10 ) THEN
!
!---      Reactive transport porosity alteration  ---
!
          IF( ISLC(43).EQ.1 ) THEN
            PORTX = POR0(1,N)*EXP(DPX*CMP(1,N)*(1.D+0-POR0(1,N)/
     &        POR0(1,N)))
            PORDX = POR0(2,N)*EXP(DPX*CMP(1,N)*(1.D+0-POR0(2,N)/
     &        POR0(2,N)))
          ELSE
            PORTX = POR(1,N)*EXP(DPX*CMP(1,N)*(1.D+0-POR(1,N)/
     &        POR(1,N)))
            PORDX = POR(2,N)*EXP(DPX*CMP(1,N)*(1.D+0-POR(2,N)/
     &        POR(2,N)))
          ENDIF
!
!---    Pore compressibility w/ variable bulk volume  ---
!
        ELSEIF( ISLC(15).EQ.11 ) THEN
!
!---      Reactive transport porosity alteration  ---
!
          IF( ISLC(43).EQ.1 ) THEN
            PORTX = POR0(1,N)*EXP(DPX*CMP(1,N)*(1.D+0-POR0(1,N)))
            PORDX = POR0(2,N)*EXP(DPX*CMP(1,N)*(1.D+0-POR0(2,N)))
          ELSE
            PORTX = POR(1,N)*EXP(DPX*CMP(1,N)*(1.D+0-POR(1,N)))
            PORDX = POR(2,N)*EXP(DPX*CMP(1,N)*(1.D+0-POR(2,N)))
          ENDIF
!
!---    Bulk compressibility w/ fixed bulk volume  ---
!
        ELSE
!
!---      Reactive transport porosity alteration  ---
!
          IF( ISLC(43).EQ.1 ) THEN
            PORTX = POR0(1,N) + DPX*CMP(1,N)
            PORDX = POR0(2,N) + DPX*CMP(1,N)
          ELSE
            PORTX = POR(1,N) + DPX*CMP(1,N)
            PORDX = POR(2,N) + DPX*CMP(1,N)
          ENDIF
        ENDIF
      ENDIF
      PORTX = MAX( MIN( PORTX,1.D+0 ),1.D-6 )
      PORDX = MAX( MIN( PORDX,1.D+0 ),1.D-6 )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PORSTY_CO2E group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PORSTY_GM
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
!     Update porosity.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 28 February 2017 (The erroneous word
!     "dord" is discovered in the Webster's New International 
!     Dictionary, Second Edition, prompting an investigation in 1939).
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE JACOB
      USE GRID
      USE FDVP
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
      SUB_LOG(ISUB_LOG) = '/PORSTY_GM'
!
!---  Loop over the finite elements (i.e., grid cells) ---
!
      DO N = 1,NFCGC(ID+1)
!
!---    Skip for inactive  ---
!
        IF( IXP(N).EQ.0 ) CYCLE
!
!---    Loop over increment indices  ---
!
        DO M = 2,ISVC+2
          PX = MAX( PG(M,N),PL(M,N) ) + PATM
          CALL PORSTY_CO2E( PX,PCMP(N),PORD(M,N),PORT(M,N),N )
          IDX = MIN( ID+1,2 )
          PORD(M,N) = MAX( PORD(M,N),EPSL )
          PORT(M,N) = MAX( PORT(M,N),PORD(M,N) )
        ENDDO
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PORSTY_GM group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PTL_A( PX,TX,VAR_TA,VARX,DVARX,I_PX,I_SX,I_TX,I_VX )
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
!     I_VX = 0 - subcritical pressure gas, bilinear interpolation
!     I_VX = 1 - subcritical pressure liquid, liquid-side interpolation
!     I_VX = 2 - supercritical pressure gas/liq, bilinear interpolation
!     VAR_TA(LT_TA,LP_TA,LNNGC) - pressure, temperature ordered 
!       array of property values
!     VARX - interpolated property value
!     DVARX - partial derivative of the property value with respect
!             to pressure
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 23 April 2002.
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
      REAL*8 VAR_TA(LT_TA,LP_TA,LNNGC)
      INTEGER I_PX(2),I_SX,I_TX(2)
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
      VPLTLX = VAR_TA(I_TX(1),I_PX(1),INCG)
      VPLTUX = VAR_TA(I_TX(1)+1,I_PX(1),INCG)
      VPUTLX = VAR_TA(I_TX(2),I_PX(2),INCG)
      VPUTUX = VAR_TA(I_TX(2)+1,I_PX(2),INCG)
      PLX = P_TA(I_PX(1),INCG)
      PUX = P_TA(I_PX(2),INCG)
      TLPLX = T_TA(I_TX(1),I_PX(1),INCG)
      TUPLX = T_TA(I_TX(1)+1,I_PX(1),INCG)
      TLPUX = T_TA(I_TX(2),I_PX(2),INCG)
      TUPUX = T_TA(I_TX(2)+1,I_PX(2),INCG)
!
!---  Liquid-side interpolation scheme  ---
!
      IF( I_VX.EQ.1 ) THEN
        ILSX = IV_TA(I_SX,INCG)-1
        IUSX = IV_TA(I_SX+1,INCG)-1
        VLSX = VAR_TA(ILSX,I_SX,INCG)
        VUSX = VAR_TA(IUSX,I_SX+1,INCG)
        TLSX = T_TA(ILSX,I_SX,INCG)
        TUSX = T_TA(IUSX,I_SX+1,INCG)
        PLSX = P_TA(I_SX,INCG)
        PUSX = P_TA(I_SX+1,INCG)
        VSX = ((TKX-TLSX)/(TUSX-TLSX))*(VUSX-VLSX) + VLSX
        PSX = ((TKX-TLSX)/(TUSX-TLSX))*(PUSX-PLSX) + PLSX
        VPUX = ((VPUTUX-VPUTLX)*(TKX-TLPUX)/(TUPUX-TLPUX)) + VPUTLX
        DVARX = (VPUX-VSX)/(PUX-PSX)
        VARX = (DVARX*(PMX-PSX)) + VSX
!
!---  Bilinear interpolation scheme  ---
!
      ELSE
        VPLX = ((VPLTUX-VPLTLX)*(TKX-TLPLX)/(TUPLX-TLPLX)) + VPLTLX
        VPUX = ((VPUTUX-VPUTLX)*(TKX-TLPUX)/(TUPUX-TLPUX)) + VPUTLX
        IF( ABS(PUX-PLX).LT.EPSL ) THEN
          DVARX = 0.D+0
          VARX = VPLX
        ELSE
          DVARX = (VPUX-VPLX)/(PUX-PLX)
          VARX = (DVARX*(PMX-PLX)) + VPLX
        ENDIF
      ENDIF
!
!---  Convert partial derivative with respect to pressure from
!     1/MPa to 1/Pa  ---
!
      DVARX = DVARX*1.D-6
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
      SUBROUTINE RHOG_I( TX,PX,XGWX,RHOGX )
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
!     H2O-NaCl-CO2 Equation of state.  This subroutine finds the gas
!     density at fixed gas water/CO2 mass fraction, temperature,
!     and pressure.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 26 July 2007.
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
      SUB_LOG(ISUB_LOG) = '/RHOG_I'
      ISRX = 2
!
!---  Estimate partial pressures  ---
!
      XGAX = MAX( 1.D+0-XGWX,0.D+0 )
      XMGWX = (XGWX/WTMW)/(XGWX/WTMW + XGAX/WTMA)
      PVWX = XMGWX*PX
!
!---  Newton scheme  ---
!
  100 CONTINUE
      PVAX = MAX( PX-PVWX,0.D+0 )
      CALL DENS_W( TX,PVWX,RHOLWX,RHOGWX,ISRX )
      CALL DENS_A( TX,PVAX,RHOGAX,ISRX )
      RHOGX = RHOGWX+RHOGAX
      DPVWY = SIGN( 1.D-1,(5.D-1*PX-PVWX) )
      PVWY = PVWX + DPVWY
      PVAY = MAX( PX-PVWY,0.D+0 )
      CALL DENS_W( TX,PVWY,RHOLWY,RHOGWY,ISRX )
      CALL DENS_A( TX,PVAY,RHOGAY,ISRX )
      RHOGY = RHOGWY+RHOGAY
      FX = XGWX - RHOGWX/RHOGX
      DFX = (RHOGWX/RHOGX - RHOGWY/RHOGY)/DPVWY
      DPVWX = -FX/DFX
      PVWX = PVWX + DPVWX
      IF( ABS(DPVWX).GT.1.D-4 ) GOTO 100
      PVAX = MAX( PX-PVWX,0.D+0 )
      CALL DENS_W( TX,PVWX,RHOLWX,RHOGWX,ISRX )
      CALL DENS_A( TX,PVAX,RHOGAX,ISRX )
      RHOGX = RHOGWX+RHOGAX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RHOG_I group  ---
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
!     Abramzon, A.A., and R.D. Gaukhberg.  1993.  Surface tension of
!     salt solutions.  Russian Journal of Applied Chemistry,
!     66(6):1139-1146.
!
!     Lide, D.R., and H.V. Kehiaian.  1994.  CRC Handbook of
!     Thermophysical and Thermochemical Data.  CRC Press, Boca Raton,
!     Florida.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 14 May 2002.
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
!     McKibbin, R., and A. McNabb.  1993.  Modeling the phase
!     boundaries and fluid properties of the system H2O-NaCl at high
!     temperatures and pressures.  Proceedings 15th NZ Geothermal
!     Workshop, University of Auckland, New Zealand.
!
!     Temperature range 0-800 C.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 April 2002.
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
!     Written by M.D. White, PNNL, 1 April 2002.
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
!      SAVE SAX,SBX
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
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 18 March 2002.
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
!      SAVE K
      DATA K / -7.691234564D+0, -2.608023696D+1, -1.681706546D+2,
     &  6.423285504D+1, -1.189646225D+2, 4.167117320D+0,
     &  2.097506760D+1, 1.D+9, 6.D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SP_W'
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
      PSWX = EXP(PSWX)*PCRW
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
!     Written by M.D. White, PNNL, 27 July 2002.
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
      REAL*8 X(N),Y(N),Y2(N),UX(N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SPLINE'
      Y2(1) = 0.D+0
      UX(1) = 0.D+0
      DO 100 I = 2,N-1
        SIGX = (X(I)-X(I-1))/(X(I+1)-X(I-1))
        PX = SIGX*Y2(I-1)+2.D+0
        Y2(I) = (SIGX-1.D+0)/PX
        UX(I) = (6.D+0*((Y(I+1)-Y(I))/(X(I+1)-X(I)) -
     &    (Y(I)-Y(I-1))/(X(I)-X(I-1)))/
     &    (X(I+1)-X(I-1)) - SIGX*UX(I-1))/PX
  100 CONTINUE
      QNX = 0.D+0
      UNX = 0.D+0
      Y2(N) = (UNX-QNX*UX(N-1))/(QNX*Y2(N-1)+1.D+0)
      DO 200 K = N-1,1,-1
        Y2(K) = Y2(K)*Y2(K+1)+UX(K)
  200 CONTINUE
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
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
!     Written by M.D. White, PNNL, 27 July 2002.
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
      REAL*8 XA(N),YA(N),Y2A(N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SPLINT'
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
!      Y1A = (YA(KHI)-YA(KLO))/H -
!     &  ((3.D+0*(A**2)-1.D+0)/6.D+0)*H*Y2A(KLO) +
!     &  ((3.D+0*(B**2)-1.D+0)/6.D+0)*H*Y2A(KHI)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
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
!     Written by M.D. White, PNNL, 18 March 2002.
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
!      SAVE LFCX,TOLX
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
      SUBROUTINE THK_A( TX,PX,RHOAX,TKAX )
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
!     Written by M.D. White, PNNL, 26 April 2002.
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
      REAL*8 TX_TA(LTX),PX_TA(LPX),TKX_TA(LTX,LPX),
     &  PSX_TA(LTX),TKLSX_TA(LTX),TKGSX_TA(LTX)
      REAL*8 TIX(2),PIX(2,2),TKIX(2,2)
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
      DATA TKX_TA /
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
      DATA TKLSX_TA / 83.46D+0, 82.3D+0, 84.9D+0, 187.1D+0, 0.D+0,
     &  0.D+0, 0.D+0 /
      DATA TKGSX_TA / 45.34D+0, 53.13D+0, 47.61D+0, 39.02D+0, 0.D+0,
     &  0.D+0, 0.D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THK_A'
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
        TKIX(1,1) = TKX_TA(ITLX,IPLX)
        TKIX(2,1) = TKX_TA(ITUX,IPLX)
        TKIX(1,2) = TKX_TA(ITLX,IPUX)
        TKIX(2,2) = TKX_TA(ITUX,IPUX)
        IF( (PSX_TA(ITLX).GT.PX_TA(IPLX)) .AND.
     &    (PSX_TA(ITLX).LT.PX_TA(IPUX)) ) THEN
          IF( RHOAX.LE.RHOCRA ) THEN
            PIX(1,2) = PSX_TA(ITLX)
            TKIX(1,2) = TKGSX_TA(ITLX)
          ELSE
            PIX(1,1) = PSX_TA(ITLX)
            TKIX(1,1) = TKLSX_TA(ITLX)
          ENDIF
        ENDIF
        IF( (PSX_TA(ITUX).GT.PX_TA(IPLX)) .AND.
     &    (PSX_TA(ITUX).LT.PX_TA(IPUX)) ) THEN
          IF( RHOAX.LE.RHOCRA ) THEN
            PIX(2,2) = PSX_TA(ITUX)
            TKIX(2,2) = TKGSX_TA(ITUX)
          ELSE
            PIX(2,1) = PSX_TA(ITUX)
            TKIX(2,1) = TKLSX_TA(ITUX)
          ENDIF
        ENDIF
!
!---    Bilinear interpolation  ---
!
        VTLX = (TKIX(1,2)-TKIX(1,1))*(PX-PIX(1,1))/(PIX(1,2)-PIX(1,1)) +
     &    TKIX(1,1)
        VTUX = (TKIX(2,2)-TKIX(2,1))*(PX-PIX(2,1))/(PIX(2,2)-PIX(2,1)) +
     &    TKIX(2,1)
        TKAX = ((VTUX-VTLX)*(TKX-TIX(1))/(TIX(2)-TIX(1))) + VTLX
        TKAX = 1.D-3*TKAX
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
        TK_ZD = 0.475598D+0*SQRT(TKX)*(1.D+0+(SRX**2))/ZETAX
!
!---    Density-dependent component  ---
!
        TK_DD = 0.D+0
        DO 60 I = 1,4
          TK_DD = TK_DD + SDX(I)*(RHOAX**I)
   60   CONTINUE
!
!---    Sum components and convert to W/m K  ---
!
        TKAX = (TK_ZD + TK_DD)*1.D-3
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THK_A group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THK_B( TX,XLSX,TKLWX,TKBX )
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
!     Ozbek, H. and S.L. Phillips.  1980.  Thermal conductivity of
!     aqueous sodium chloride solutions from 20 to 330 C.
!     J. Chem. Engr. Data, 25:263-267.
!
!     Temperature Range:  20 - 330 C
!     NaCl Concentration Range:  0 - saturation
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 4 April 2002.
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
!      SAVE SCX
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
      TKBX = 1.D+0 - (SCX(1) + SCX(2)*TX + SCX(3)*(TX**2))*YLSX +
     &  (SCX(4) + SCX(5)*TX + SCX(6)*(TX**2))*(YLSX**2)
      TKBX = TKBX*TKLWX
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
      SUBROUTINE THK_G( TX,TKAX,TKWX,XMAX,XMWX,TKGX )
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
!     Written by M.D. White, PNNL, 29 April 2002.
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
      SUB_LOG(ISUB_LOG) = '/THK_G'
!
!---  Reduced, inverse thermal conductivity  ---
!
      TKX = TX+TABS
      GAMMA_W = 210.D+0*(TCRW*(WTMW**3)/((PCRW*1.D-5)**4))
     &  **(1.D+0/6.D+0)
      GAMMA_A = 210.D+0*(TCRA*(WTMA**3)/((PCRA*1.D-5)**4))
     &  **(1.D+0/6.D+0)
      TKR_A = TKX/TCRA
      TKR_W = TKX/TCRW
      TKR_A = GAMMA_A*(EXP(4.64D-2*TKR_A)-EXP(-2.412D-1*TKR_A))
      TKR_W = GAMMA_W*(EXP(4.64D-2*TKR_W)-EXP(-2.412D-1*TKR_W))
!
!---  Zero-density-limit component  ---
!
      IF( TKR_W.GT.EPSL ) THEN
        PHIAWX = ((1+SQRT(TKR_A/TKR_W)*((WTMW/WTMA)**2.5D-1))**2)
     &      /SQRT(8.D+0*(1.D+0 + WTMA/WTMW))
      ELSE
        PHIAWX = 0.D+0
      ENDIF
      IF( TKR_A.GT.EPSL ) THEN
        PHIWAX = ((1+SQRT(TKR_W/TKR_A)*((WTMA/WTMW)**2.5D-1))**2)
     &      /SQRT(8.D+0*(1.D+0 + WTMW/WTMA))
      ELSE
        PHIWAX = 0.D+0
      ENDIF
      CHIWX = XMWX + XMAX*PHIWAX
      CHIAX = XMWX*PHIAWX + XMAX
      TKGX = XMWX*TKWX/CHIWX + XMAX*TKAX/CHIAX
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
      SUBROUTINE THK_W( TX,PX,RHOX,TKX )
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
!     Written by M.D. White, PNNL, 27 March 2002.
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
!      SAVE SAX,SBX,CBX,SDX,CCX
!      SAVE TREF,RHOREF,PREF,TKREF
      DATA TREF / 6.4727D+2 /
      DATA RHOREF / 3.17763D+2 /
      DATA PREF / 2.2115D+7 /
      DATA TKREF / 1.D+0 /
      DATA SAX / 0.0102811D+0, 0.0299621D+0, 0.0156146D+0,
     &  -0.00422464D+0 /
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
!---  Reduced temperature, density, and pressure  ---
!
      THETAX = (TX+TABS)/TREF
      RHOBX = RHOX/RHOREF
      BETAX = PX/PREF
!
!---  Zero term  ---
!
      TKX = 0.D+0
      ISAX = 0
      DO 10 I = 0,3
        ISAX = ISAX + 1
        TKX = TKX + SAX(ISAX)*(THETAX**I)
   10 CONTINUE
      TKX = SQRT(THETAX)*TKX
!
!---  First term  ---
!
      TKX = TKX + SBX(1) + SBX(2)*RHOBX +
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
      TKX = TKX + ((SDX(1)/(THETAX**10))+SDX(2))*
     &  (RHOBX**(9.D+0/5.D+0))*
     &  EXP(CCX(1)*(1.D+0-(RHOBX**(14.D+0/5.D+0))))
      TKX = TKX + SDX(3)*CSX*(RHOBX**CQX)*
     &  EXP((CQX/(1.D+0+CQX))*(1.D+0-(RHOBX**(1.D+0+CQX))))
      VARX = CCX(2)*(THETAX**(3.D+0/2.D+0)) + CCX(3)/(RHOBX**5)
      VARX = MAX( VARX,-50.D+0 )
      TKX = TKX + SDX(4)*EXP(VARX)
!
!---  Dimensionalize thermal conductivity  ---
!
      TKX = TKX*TKREF
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
      SUBROUTINE THKS_LG( THKSX,N )
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
!     Calculation of average intrinsic permeability around well,
!     using logarithmic averaging.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 26 July 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
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
      SUB_LOG(ISUB_LOG) = '/THKS_LG'
!
!---  Average intrinsic permeability  ---
!
      THKSX = 0.D+0
      NTHKSX = 0
      IF( (THKS(1,N)/EPSL).GT.EPSL ) THEN
        THKSX = THKSX + LOG(THKS(1,N))
        NTHKSX = NTHKSX+1
      ENDIF
      IF( (THKS(2,N)/EPSL).GT.EPSL ) THEN
        THKSX = THKSX + LOG(THKS(2,N))
        NTHKSX = NTHKSX+1
      ENDIF
      IF( (THKS(3,N)/EPSL).GT.EPSL ) THEN
        THKSX = THKSX + LOG(THKS(3,N))
        NTHKSX = NTHKSX+1
      ENDIF
      IF( NTHKSX.GT.0 ) THEN
        REALX = REAL(NTHKSX)
        THKSX = EXP( THKSX/REALX )
      ELSE
        THKSX = 0.D+0
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THKS_LG group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TORTU( SLX,SGX,PORDX,TORLX,TORGX,N )
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
!     Compute phase tortuosity.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, Battelle, PNNL, 28 October 2021.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE MPI
      USE GLB_PAR
      USE GRID
      USE SOLTN
      USE PROP
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
      SUB_LOG(ISUB_LOG) = '/TORTU'
!
!---  Constant phase tortuosity  ---
!
      IF( ITOR(N).EQ.1 ) THEN
        TORLX = TOR(1,N)
        TORGX = TOR(2,N)
!
!---  Millington and Quirk tortuosity model  ---
!
      ELSEIF( ITOR(N).EQ.2 ) THEN
        IF( SLX*PORDX.LT.EPSL ) THEN
          TORLX = 0.D+0
        ELSE
          TORLX = (PORDX*(SLX**7))**(THIRD)
        ENDIF
        IF( SGX*PORDX.LT.EPSL ) THEN
          TORGX = 0.D+0
        ELSE
          TORGX = (PORDX*(SGX**7))**(THIRD)
        ENDIF
!
!---  Millington and Quirk (free gas) tortuosity model  ---
!
      ELSEIF( ITOR(N).EQ.3 ) THEN
        IF( SLX*PORDX.LT.EPSL ) THEN
          TORLX = 0.D+0
        ELSE
          TORLX = (PORDX*(SLX**7))**(THIRD)
        ENDIF
        SGFX = MAX( SGX,0.D+0 )
        IF( SGFX*PORDX.LT.EPSL ) THEN
          TORGX = 0.D+0
        ELSE
          TORGX = (PORDX*(SGFX**7))**(THIRD)
        ENDIF
!
!---  Marshal tortuosity model  ---
!
      ELSEIF( ITOR(N).EQ.4 ) THEN
        IF( SLX*PORDX.LT.EPSL ) THEN
          TORLX = 0.D+0
        ELSE
          TORLX = SQRT(PORDX*SLX)
        ENDIF
        IF( SGX*PORDX.LT.EPSL ) THEN
          TORGX = 0.D+0
        ELSE
          TORGX = SQRT(PORDX*SGX)
        ENDIF
!
!---  Constant aqueous/Millington gas tortuosity model  ---
!
      ELSEIF( ITOR(N).EQ.5 ) THEN
        TORLX = TOR(1,N)
        IF( SGX*PORDX.LT.EPSL ) THEN
          TORGX = 0.D+0
        ELSE
          TORGX = (PORDX*(SGX**7))**(THIRD)
        ENDIF
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TORTU group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE VISC_A( TX,RHOAX,VISAX )
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
!     Fenghour, A., W. A. Wakeham, V. Vesovic.  1998.  The viscosity
!     of carbon dioxide.  J. Phys. Chem. Ref. Data, 27(1):31-41.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 April 2002.
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
!      SAVE SAX,SBX,ESPX
      DATA SAX / 0.235156D+0, -0.491266D+0,
     &  5.211155D-2, 5.347906D-2, -1.537102D-2 /
      DATA SBX / 0.4071119D-2, 0.7198037D-4,
     &  0.2411697D-16, 0.2971072D-22, -0.1627888D-22 /
      DATA ESPX / 251.196D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VISC_A'
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
!---  End of VISC_A group  ---
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
!     Written by M.D. White, PNNL, 4 April 2002.
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
      SUBROUTINE VISC_G( VISAX,VISWX,XMAX,XMWX,VISGX )
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
!     Written by M.D. White, PNNL, 25 April 2002.
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
      IF( ISLC(9).EQ.1 ) THEN
        VISGX = VISGI
        RETURN
      ENDIF
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VISC_G'
!
!---  Zero-density-limit component  ---
!
      PHIAWX = ((1+SQRT(VISAX/VISWX)*((WTMW/WTMA)**2.5D-1))**2)
     &      /SQRT(8.D+0*(1.D+0 + WTMA/WTMW))
      PHIWAX = ((1+SQRT(VISWX/VISAX)*((WTMA/WTMW)**2.5D-1))**2)
     &      /SQRT(8.D+0*(1.D+0 + WTMW/WTMA))
      CHIWX = XMWX + XMAX*PHIWAX
      CHIAX = XMWX*PHIAWX + XMAX
      VISGX = XMWX*VISWX/CHIWX + XMAX*VISAX/CHIAX
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
      SUBROUTINE VISC_L( XMLAX,VISBX,VISAX,VISLX )
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
!     Kumagai, A., and C. Yokoyama.  1999.  Viscosities of aqueous
!     NaCl solutions containing CO2 at high pressures.  J. Chem. Eng.
!     Data, 44:227-229.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids. pp. 474-475.
!     McGraw-Hill, New York, New York
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 30 April 2002.
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
      IF( ISLC(9).EQ.1 ) THEN
        VISLX = VISLI
        RETURN
      ENDIF
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VISC_L'
!
!---  Formulation of Grunberg and Nissan  ---
!
      VISLX = (1.D+0-XMLAX)*LOG(VISBX)
      IF( VISAX.GT.EPSL ) VISLX = VISLX + XMLAX*LOG(VISAX)
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
!     Written by M.D. White, PNNL, 27 March 2002.
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
!---  Reduced temperature, density, and pressure  ---
!
      THETAX = (TX+TABS)/TREF
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
        CALL SR_W( TX,PIX,ISRX )
        CALL DENS_W( TX,PIX,RHOLX,RHOVX,ISRX )
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
      SUBROUTINE VPL_B( TX,PSBX,PCX,RHOBX,PVBX,XLSX,NX )
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
!     Battistelli, A., C. Claudio, and K. Pruess.  1997.  The simulator
!     TOUGH2/EWASG for modelling geothermal reservoirs with brines and
!     noncondensible gas.  Geothermics, 26(4): 437-464.
!
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 25 April 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PROP
      USE GRID
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
      IF( PCX.GT.EPSL ) THEN
        IF( ISM(NX).EQ.0 ) THEN
          PVBX = PSBX*EXP( -WTMW*(PCX**1.25D+0)/(RHOBX*RCU*TKX) )
        ELSE
          PVBX = PSBX*EXP( -WTMW*PCX/(RHOBX*RCU*TKX) )
        ENDIF
      ELSE
        PVBX = PSBX
      ENDIF
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
      SUBROUTINE WEBB_BC( CLX,HMPX,PSIX,SMPX,SRX )
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
!     Written by M.D. White, PNNL, 16 July 2010
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WEBB_BC'
!
!---  Initial guess of the matching-point head  ---
!
      SMPX = 1.D-2*(1.D+0-SRX) + SRX
      SMPIX = SMPX
      ESMPX = (SMPX-SRX)/(1.D+0-SRX)
      HMPX = PSIX/(ESMPX**(1.D+0/CLX))
      HMPX = 1.D+1**(4.D-1*LOG10(HDOD))
      HMPIX = HMPX
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
        SMPX = SMPIX
        M_ERR(1) = 'No Convergence on Saturation Matching Point: '
        M_ERR(2) = ' Brooks and Corey: Initial Saturation Guess = '
        IF( N_DB.LT.0 ) THEN
          M_ERR(3) = ' at Boundary Surface: '
        ELSE
          M_ERR(3) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = SMPIX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 2
        I_ERR(3) = 3
        I_ERR(4) = ID
        DSMPX = 0.D+0
      ENDIF
      IF( ABS(DHMPX).GT.1.D-4 ) GOTO 100
!
!---  Find the capillary matching point saturation  ---
!
      ESMPX = (PSIX/HMPX)**CLX
      SMPX = ESMPX*(1.D+0-SRX) + SRX
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
      SUBROUTINE WEBB_VG( ALPHAX,CMX,CNX,HMPX,SMPX,SRX )
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
!     Written by M.D. White, PNNL, 16 July 2010
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WEBB_VG'
!
!---  Use the matrix saturation at 0.4 LOG10(HDOD) as
!     the initial guess  ---
!
      HDX = 1.D+1**(4.D-1*LOG10(HDOD))
      HDIX = HDX
      SMPX = (1.D+0/(1.D+0 + (ALPHAX*HDX)**CNX))**CMX
      SMPX = SMPX*(1.D+0-SRX) + SRX
      SMPX = 1.D-1*(1.D+0-SRX) + SRX
      SMPIX = SMPX
!
!---  Newton-Raphson iteration for the matrix saturation
!     matching point  ---
!
      NC = 0
  100 CONTINUE
      NC = NC + 1
      SEMPX = (SMPX-SRX)/(1.D+0-SRX)
      ESEMPX = (1.D+0/SEMPX)**(1.D+0/CMX)
      FX1 = LOG10(HDOD)
      FX1 = FX1 - LOG10(((ESEMPX-1.D+0)**(1.D+0/CNX))/ALPHAX)
      FX1 = FX1 - (SMPX/(SMPX-SRX))/(LOG(1.D+1)*CNX*CMX*
     &  (1.D+0-(SEMPX**(1.D+0/CMX))))
      SMPY = SMPX + 1.D-8
      SEMPX = (SMPY-SRX)/(1.D+0-SRX)
      ESEMPX = (1.D+0/SEMPX)**(1.D+0/CMX)
      FX2 = LOG10(HDOD)
      FX2 = FX2 - LOG10(((ESEMPX-1.D+0)**(1.D+0/CNX))/ALPHAX)
      FX2 = FX2 - (SMPY/(SMPY-SRX))/(LOG(1.D+1)*CNX*CMX*
     &  (1.D+0-(SEMPX**(1.D+0/CMX))))
      DFX = (FX2-FX1)/1.D-8
      DSMPX = -FX1/DFX
      SMPX = MAX( SMPX+DSMPX,SRX+1.D-12 )
!
!---  No convergence on matrix saturation matching point  ---
!
      IF( NC.GT.32 ) THEN
        SMPX = SMPIX
        M_ERR(1) = 'No Convergence on Saturation Matching Point: '
        M_ERR(2) = ' van Genuchten: Initial Saturation Guess = '
        IF( N_DB.LT.0 ) THEN
          M_ERR(3) = ' at Boundary Surface: '
        ELSE
          M_ERR(3) = ' at Node: '
        ENDIF
        CALL PATH
        R_ERR = SMPIX
        I_ERR(1) = ABS(N_DB)
        I_ERR(2) = 2
        I_ERR(3) = 3
        I_ERR(4) = ID
        DSMPX = 0.D+0
      ENDIF
      IF( ABS(DSMPX).GT.1.D-9 ) GOTO 100
!
!---  Find the matrix capillary head matching point  ---
!
      SEMPX = (SMPX-SRX)/(1.D+0-SRX)
      HMPX = (1.D+0/ALPHAX)*
     &  (((1.D+0/(SEMPX**(1.D+0/CMX)))-1.D+0)**(1.D+0/CNX))
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WEBB_VG group  ---
!
      RETURN
      END


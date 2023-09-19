!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ACTV_W( ACTVX,TX,XMLAX,XMLOX,XMLSX,XMLWX )
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
!     Activity coefficient of water in aqueous.
!
!     Vanderbeken, I., S. Ye, B. Bouyssiere, H. Carrier, and P. Xans.
!     1999.  "Ability of the MHV2 mixing rule to describe the effect
!     of salt on gas solubility in brines at high temperature and
!     high pressure.  High Temperatures - High Pressures, 31:653-663.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 October 2004.
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
      REAL*8 CAX(3),COX(3)
      REAL*8 AAX(3,3),AOX(3,3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ACTV_W'
!
!---  Modified Wilson parameters (1-gas, 2-water, 3-salt)  ---
!
      TKX = TX+TABS
      CAX(1) = 290.4261D+0 - 0.233109D+0*TKX
      CAX(2) = 0.663874D+0 + 348.599D+0/TKX
      CAX(3) = -434.8795D+0 + 147416.58D+0/TKX
      COX(1) = -5.837922D+0 + 2.9336D-2*TKX - 2.83D-5*TKX*TKX
      COX(2) = -50.092813D+0 - 4.00308D-1*TKX + 8.45D-4*TKX*TKX
      COX(3) = -479.1109D+0 + 2.55665D+0*TKX - 4.987D-3*TKX*TKX
     &  + 3.51D-6*(TKX**3)
      AAX(1,1) = 1.D+0
      AAX(1,2) = 0.987774D+0 - 5.357714D+0/TKX
      AAX(2,1) = 1.303747D+0 +379.84D+0/TKX
      AAX(2,2) = 1.D+0
      AAX(1,3) = -1.096508D+0 + 1393.4636D+0/TKX 
     &  - 246332.86D+0/(TKX**2)
      AAX(3,1) = -1.064705D+0 + 770.11D+0/TKX
      AAX(2,3) = -7.362276D+0 + 1757.82D+0/TKX
      AAX(3,2) = -0.358477D+0 + 478.32D+0/TKX
      AAX(3,3) = 1.D+0
      AOX(1,1) = 1.D+0
      AOX(1,2) = 0.789660D+0 + 3.567D-3*TKX - 3.29D-6*TKX*TKX
      AOX(2,1) = 1.122929D+0 - 3.98D-4*TKX + 3.84D-7*TKX*TKX
      AOX(2,2) = 1.D+0
      AOX(1,3) = 1.350959D+0 - 8.424D-3*TKX + 1.87D-5*TKX*TKX
      AOX(3,1) = 32.48964D+0 - 0.171947D+0*TKX + 2.4D-4*TKX*TKX
      AOX(2,3) = 10.380108D+0 - 7.4383D-2*TKX + 1.85D-4*TKX*TKX
     &  - 1.54D-7*(TKX**3)
      AOX(3,2) = 87.848001D+0 - 0.509548D+0*TKX + 7.65D-4*TKX*TKX
      AOX(3,3) = 1.D+0
      XMLNX = XMLAX+XMLOX
      ACVAX = XMLNX*CAX(1) + XMLWX*CAX(2) + XMLSX*CAX(3)
      ACVAX = ACVAX - CAX(2)*LOG( XMLNX*AAX(2,1) + XMLWX*AAX(2,2)
     &  + XMLSX*AAX(2,3) )
      ACVAX = ACVAX - CAX(1)*XMLNX*AAX(1,2)/
     &  (XMLNX*AAX(1,1) + XMLWX*AAX(1,2) + XMLSX*AAX(1,3) + SMALL)
      ACVAX = ACVAX - CAX(2)*XMLWX*AAX(2,2)/
     &  (XMLNX*AAX(2,1) + XMLWX*AAX(2,2) + XMLSX*AAX(2,3) + SMALL)
      ACVAX = ACVAX - CAX(3)*XMLSX*AAX(3,2)/
     &  (XMLNX*AAX(3,1) + XMLWX*AAX(3,2) + XMLSX*AAX(3,3) + SMALL)
      ACVAX = EXP( ACVAX )      
      ACVOX = XMLNX*COX(1) + XMLWX*COX(2) + XMLSX*COX(3)
      ACVOX = ACVOX - COX(2)*LOG( XMLNX*AOX(2,1) + XMLWX*AOX(2,2)
     &  + XMLSX*AOX(2,3) )
      ACVOX = ACVOX - COX(1)*XMLNX*AOX(1,2)/
     &  (XMLNX*AOX(1,1) + XMLWX*AOX(1,2) + XMLSX*AOX(1,3) + SMALL)
      ACVOX = ACVOX - COX(2)*XMLWX*AOX(2,2)/
     &  (XMLNX*AOX(2,1) + XMLWX*AOX(2,2) + XMLSX*AOX(2,3) + SMALL)
      ACVOX = ACVOX - COX(3)*XMLSX*AOX(3,2)/
     &  (XMLNX*AOX(3,1) + XMLWX*AOX(3,2) + XMLSX*AOX(3,3) + SMALL)
      ACVOX = EXP( ACVOX )   
      ACTVX = (XMLAX*ACVAX + XMLOX*ACVOX)/(XMLNX+SMALL)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ACTV_W group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BNDFAO( TX,PX,DIFX )
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
!     Calculates the air-oil binary diffusion coefficient from
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
      USE NAPL
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
      SAVE A
      DATA A /1.06036D+0,1.5610D-1,1.9300D-1,4.7635D-1,1.03587D+0,
     &        1.52996D+0,1.76474D+0,3.89411D+0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BNDFAO'
      TKX = TX+TABS
      PB = PX*1.D-5
      TPAX = 1.15D+0*TBA
      TPOX = 1.15D+0*TBO
      TP = TKX/SQRT(TPAX*TPOX)
      SIGAX = 1.18D+0*(2.85D-1*(VCRA**1.048D+0))**3.333D-1
      SIGOX = 1.18D+0*(2.85D-1*(VCRO**1.048D+0))**3.333D-1
      SIG = 5.D-1*(SIGAX+SIGOX)
      OMG = (A(1)/TP**A(2)) + (A(3)/EXP(A(4)*TP))
     &  + (A(5)/EXP(A(6)*TP))  + (A(7)/EXP(A(8)*TP))
      WTMX = 2.D+0/((1.D+0/WTMA) + (1.D+0/WTMO))
      DIFX = (3.03D+0 - (9.8D-1/SQRT(WTMX)))*1.D-3*(TKX**1.5D+0)/
     &  (PB*SQRT(WTMX)*(SIG**2)*(OMG))*1.D-4
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BNDFAO group  ---
!
      RETURN
      END


!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BNDFAW( TX,PX,DIFX )
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
!     Calculates the air-water binary diffusion coefficient from
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
      REAL*8 A(8)
!
!----------------------Data Statements---------------------------------!
!
      SAVE A
      DATA A /1.06036D+0,1.5610D-1,1.9300D-1,4.7635D-1,1.03587D+0,
     &        1.52996D+0,1.76474D+0,3.89411D+0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BNDFAW'
      TKX = TX+TABS
      PB = PX*1.D-5
      TPAX = 1.15D+0*TBA
      TPWX = 1.15D+0*TBW
      TP = TKX/SQRT(TPAX*TPWX)
      SIGAX = 1.18D+0*(2.85D-1*(VCRA**1.048D+0))**3.333D-1
      SIGWX = 1.18D+0*(2.85D-1*(VCRW**1.048D+0))**3.333D-1
      SIG = 5.D-1*(SIGAX+SIGWX)
      OMG = (A(1)/TP**A(2)) + (A(3)/EXP(A(4)*TP))
     &  + (A(5)/EXP(A(6)*TP))  + (A(7)/EXP(A(8)*TP))
      WTMX = 2.D+0/((1.D+0/WTMA) + (1.D+0/WTMW))
      DIFX = (3.03D+0 - (9.8D-1/SQRT(WTMX)))*1.D-3*(TKX**1.5D+0)/
     &  (PB*SQRT(WTMX)*(SIG**2)*(OMG))*1.D-4
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BNDFAW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BNDFOW( TX,PX,DIFX )
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
!     Calculates the oil-water binary diffusion coefficient from
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
      USE NAPL
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
      SAVE A
      DATA A /1.06036D+0,1.5610D-1,1.9300D-1,4.7635D-1,1.03587D+0,
     &        1.52996D+0,1.76474D+0,3.89411D+0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BNDFOW'
      TKX = TX+TABS
      PB = PX*1.D-5
      TPOX = 1.15D+0*TBO
      TPWX = 1.15D+0*TBW
      TP = TKX/SQRT(TPOX*TPWX)
      SIGOX = 1.18D+0*(2.85D-1*(VCRO**1.048D+0))**3.333D-1
      SIGWX = 1.18D+0*(2.85D-1*(VCRW**1.048D+0))**3.333D-1
      SIG = 5.D-1*(SIGOX+SIGWX)
      OMG = (A(1)/TP**A(2)) + (A(3)/EXP(A(4)*TP))
     &  + (A(5)/EXP(A(6)*TP))  + (A(7)/EXP(A(8)*TP))
      WTMX = 2.D+0/((1.D+0/WTMO) + (1.D+0/WTMW))
      DIFX = (3.03D+0 - (9.8D-1/SQRT(WTMX)))*1.D-3*(TKX**1.5D+0)/
     &  (PB*SQRT(WTMX)*(SIG**2)*(OMG))*1.D-4
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BNDFOW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DENS_A( TX,PX,RHOAX,ISRX )
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
!     Written by MD White, PNNL, 25 April 2002.
!     Last Modified by MD White, PNNL, 25 April 2002.
!     Last Modified by MD White, PNNL, 7 September 2004.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
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
      INTEGER I_PX(2),I_SX,I_TX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DENS_A'
      INCG = 1
      P_TPX = GCPP(17,INCG)
      P_CRX = GCPP(4,INCG)
!
!---  Liquid density  ---
!
      IF( ISRX.EQ.1 ) THEN
!
!---    Table lookup and bilinear interpolation  ---
!
        PTX = MAX( PX,P_TA(1,INCG)*1.D+6 )
        I_VX = 0
        CALL ITL_A( PTX,P_TPX,P_CRX,TX,I_PX,I_SX,I_TX,I_VX )
        CALL PTL_A( PTX,TX,RHO_TA,RHOAX,DRHOAX,I_PX,I_SX,I_TX,I_VX )
!
!---  Vapor or supercritical density  ---
!
      ELSE
!
!---    Table lookup and bilinear interpolation  ---
!
        IF( (PX*1.D-6).GT.P_TA(1,INCG) ) THEN
          I_VX = 1
          CALL ITL_A( PX,P_TPX,P_CRX,TX,I_PX,I_SX,I_TX,I_VX )
          CALL PTL_A( PX,TX,RHO_TA,RHOAX,DRHOAX,I_PX,I_SX,I_TX,I_VX )
!
!---    Ideal gas law  ---
!
        ELSE
          TKX = TX+TABS
          RHOAX = (PX*WTMA)/(RCU*TKX)
        ENDIF
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
      SAVE CCX,CHX,SCX,SAX,VCX
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
      SUBROUTINE DENS_G( PVAX,PVOX,PVWX,RHOGX,RHOMGX,RHOGAX,RHOGOX,
     &  RHOGWX,TX,WTMGX,XGAX,XMGAX,XGOX,XMGOX,XGWX,XMGWX )
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
!     Gas density and component fractions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 17 December 2010
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VVX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DENS_G'
      ISRX = 2
      CALL DENS_W( TX,PVWX,RHOX,RHOGWX,ISRX )
      CALL DENS_A( TX,PVAX,RHOGAX,ISRX )
      CALL DENS_O( TX,PVOX,RHOGOX )
      RHOGX = RHOGWX+RHOGAX+RHOGOX
      XGAX = RHOGAX/RHOGX
      XGOX = RHOGOX/RHOGX
      XGWX = RHOGWX/RHOGX
      WTMGX = XGAX/WTMA + XGOX/WTMO + XGWX/WTMW
      XMGAX = XGAX/WTMA/WTMGX
      XMGOX = XGOX/WTMO/WTMGX
      XMGWX = XGWX/WTMW/WTMGX
      RHOMGX = RHOGX*WTMGX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DENS_G group  ---
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
      SAVE A
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
      SUBROUTINE DENS_L( TX,RHOBX,XLAX,XLOX,RHOLX )
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
      SAVE PMV_C
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
      CXA = PMV_A*RHOBX*XLAX/WTMA
!      CXA = RHOBX*XLAX/WTMA
!
!---  Partial molar volume of CH4 (m^3/kmol), at 25 C  ---
!
      PMV_O = 37.D-3
      CXO = PMV_O*RHOBX*XLOX/WTMO
!      CXO = RHOBX*XLOX/WTMO
!
!---  Poynting correction  ---
!
      RHOLX = RHOBX/(1.D+0 + CXA + CXO - XLAX - XLOX)
!      RHOLX = RHOBX+WTMA*CXA-CXA*RHOBX*PMV_A+WTMO*CXO-CXO*RHOBX*PMV_O
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
      SUBROUTINE DENS_N( AX,BX,PX,RHONX,RHOMNX,TX,WTMX,INDX )
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
!     Density (kg/m^3) of NAPL as a function of temperature, pressure
!     and NAPL molar concentration.  Nonlinear solution of the
!     Peng-Robinson equation of state, using zero binary
!     interaction parameters.
!
!     Use 0.3 kmol/m^3 as the initial guess for the NAPL molar volume.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 17 August 2009
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
      REAL*8 VVX(3)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DENS_N'
!
!---  Direct cubic equation solve for Peng-Robinson 
!     equation of state  ---
!
      TKX = TX + TABS
!
!---  Ideal gas law reduction  ---
!
      IF( AX/EPSL.LT.EPSL ) THEN
        VX = RCU*TKX/PX
        RHOMNX = 1.D+0/VX
        RHONX = WTMX/VX
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Polynomial coefficients for Peng-Robinson equation of state  ---
!
      CAX = 1.D+0
      CBX = (PX*BX - RCU*TKX)/PX
      CCX = (AX - 2.D+0*RCU*TKX*BX - 3.D+0*PX*(BX**2))/PX
      CDX = (PX*(BX**3) + RCU*TKX*(BX**2) - AX*BX)/PX
!
!---  Nickalls cubic polynomial solver  ---
!
      XNX = -CBX/(3.D+0*CAX)
      YNX = CAX*(XNX**3) + CBX*(XNX**2) + CCX*XNX + CDX
      YN2X = YNX**2
      D2X = ((CBX**2)-(3.D+0*CAX*CCX))/((3.D+0*CAX)**2)
      IF( D2X.LE.0.D+0 ) THEN
        DX = 0.D+0
        HX = 0.D+0
      ELSE
        DX = SQRT(D2X)
        HX = -2.D+0*(DX**3)
      ENDIF
      H2X = 4.D+0*(CAX**2)*(D2X**3)
!
!---  YN2X > H2X, 1 real root  ---
!
      IF( YN2X.GT.H2X ) THEN
        V1X = (5.D-1/CAX)*(-YNX+SQRT(YN2X-H2X))
        V2X = (5.D-1/CAX)*(-YNX-SQRT(YN2X-H2X))
        VX = XNX + SIGN((ABS(V1X)**(1.D+0/3.D+0)),V1X)
     &    + SIGN((ABS(V2X)**(1.D+0/3.D+0)),V2X)
!
!---  YN2X < H2X, 3 distinct real roots  ---
!
      ELSEIF( YN2X.LT.H2X ) THEN
!        DX = SQRT(D2X)
!        HX = -2.D+0*(DX**3)
        THETAX = (ACOS(YNX/HX))/3.D+0
        V1X = XNX + 2.D+0*DX*COS(THETAX)
        V2X = XNX + 2.D+0*DX*COS(((2.D+0*GPI)/3.D+0)+THETAX)
        V3X = XNX + 2.D+0*DX*COS(((4.D+0*GPI)/3.D+0)+THETAX)
        VX = MIN(V1X,V2X,V3X)
!
!---  YN2X = H2X, 3 real roots (two or three equal roots)  ---
!
      ELSE
!
!---    HX /= 0 (two equal roots)  ---
!
!        HX = -2.D+0*(DX**3)
        IF( ABS(HX)/EPSL.GT.EPSL ) THEN
          DX = YNX/(2.D+0*CAX)
          DX = SIGN((ABS(DX)**(1.D+0/3.D+0)),DX)
          V1X = XNX + DX
          V2X = XNX + DX
          V3X = XNX - 2.D+0*DX
          VX = MIN(V1X,V2X,V3X)
!
!---    HX = 0 (three equal roots)  ---
!
        ELSE
          V1X = XNX
          V2X = XNX
          V3X = XNX
          VX = MIN(V1X,V2X,V3X)
        ENDIF
      ENDIF
      RHOMNX = 1.D+0/VX
      RHONX = WTMX*RHOMNX
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
      SUBROUTINE DENS_O( TX,PX,RHOOX )
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
!     Density of CH4 gas.
!
!     Setzmann, U., and W. Wagner.  1991.  A New Equation of State and
!     Tables of Thermodynamic Properties for Methane Covering the
!     Range from the Melting Line to 625 K at Pressures up to 1000 MPa.
!     J. Phys. Chem. Ref. Data 20(6):1061-1155.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 25 April 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
      USE NAPL
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
      INTEGER I_PX(2),I_SX,I_TX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DENS_O'
      INCG = 2
      P_TPX = GCPP(17,INCG)
      P_CRX = GCPP(4,INCG)
!
!---  Table lookup and bilinear interpolation  ---
!
      IF( (PX*1.D-6).GT.P_TA(1,INCG) ) THEN
        I_VX = 1
        CALL ITL_A( PX,P_TPX,P_CRX,TX,I_PX,I_SX,I_TX,I_VX )
        CALL PTL_A( PX,TX,RHO_TA,RHOOX,DRHOOX,I_PX,I_SX,I_TX,I_VX )
!
!---  Ideal gas law  ---
!
      ELSE
        TKX = TX+TABS
        RHOOX = (PX*WTMO)/(RCU*TKX)
      ENDIF
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
!
!----------------------Data Statements---------------------------------!
!
      SAVE CAX,SAX,CBX,SBX
      SAVE LFCX,SLX,INX,IZX,ITX,IXX
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
        INDX = 14
        CHMSG = 'Steam Table Subregion 3/5: Exceeds Temperature Range: '
        RLMSG = TX
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Subregions 4 or 5  ---
!
      IF( ISRX.EQ.4 .OR. ISRX.EQ.5 ) THEN
        INDX = 14
        CHMSG = 'Steam Table Subregion 4/5: Exceeds Temperature Range: '
        RLMSG = TX
        CALL WRMSGS( INDX )
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
      USE NCG_PT
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
      SUB_LOG(ISUB_LOG) = '/DIFC_LA'
!
!---  Restrict temperature to above 0.01 C  ---
!
      TY = MAX( TX,1.D-2 )
!
!---  CO2 properties  ---
!
      IF( INCG.EQ.1 ) THEN
!
!---    Convert viscosities to centipoise  ---
!
        DFLAX = 6.391D+3*(((VISLX*1.D+3)**(-1.584D-1))*
     &    ((VISGAX*1.D+3)**(6.911D+0)))
!
!---    Limit diffusion coefficient to the upper limit
!       reported in the reference  ---
!
        DFLAX = MIN( DFLAX,10.D-9 )
!
!---  CH4 properties  ---
!
      ELSEIF( INCG.EQ.2 ) THEN
        VBO = 2.85D-1*VCRO**1.048D+0
        DFLAX = 7.4D-15*((2.6D+0*WTMW)**5.D-1)*(TY+TABS)/
     &    (VISLX*(VBO**6.D-1))
      ENDIF
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
      SAVE CBX,TRX,VISWRX
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
!     in NAPL, using the Wilke and Chang estimation
!      method; pp. 598.
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
      REAL*8 DFNCX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DIFC_N'
!
!---  NAPL diffusion coefficients  ---
!
      DO 100 IGC = 1,NGC
        VBCX = 2.85D-1*GCPP(5,IGC)**1.048D+0
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
      SUBROUTINE ENTH_A( TX,PX,HAX,UAX,ISRX )
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
!     Written by MD White, PNNL, 25 April 2002.
!     Last Modified by MD White, PNNL, 25 April 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
      USE CCP
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
      RHOCRA = 467.6D+0
      P_TPX = GCPP(17,INCG)
      P_CRX = GCPP(4,INCG)
!
!---  Liquid enthalpy  ---
!
      IF( ISRX.EQ.1 ) THEN
!
!---    Table lookup and bilinear interpolation  ---
!
        I_VX = 0
        PTX = MAX( PX,P_TA(1,INCG)*1.D+6 )
        CALL ITL_A( PTX,P_TPX,P_CRX,TX,I_PX,I_SX,I_TX,I_VX)
        CALL PTL_A( PTX,TX,H_TA,HAX,DHAX,I_PX,I_SX,I_TX,I_VX )
        UAX = HAX
        DUAX = DHAX
!
!---  Vapor and supercritical enthalpy  ---
!
      ELSE  
!
!---    Table lookup and bilinear interpolation  ---
!
        I_VX = 1
        PTX = MAX( PX,P_TA(1,INCG)*1.D+6 )
        CALL ITL_A( PTX,P_TPX,P_CRX,TX,I_PX,I_SX,I_TX,I_VX )
        CALL PTL_A( PTX,TX,H_TA,HAX,DHAX,I_PX,I_SX,I_TX,I_VX )
        CALL PTL_A( PTX,TX,U_TA,UAX,DUAX,I_PX,I_SX,I_TX,I_VX )
      ENDIF
!
!---  Convert units to J/kg  ---
!
      HAX = 1.D+3*HAX
      UAX = 1.D+3*UAX
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
      SAVE SAX,SBX,SCX,CAX,CBX
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
      SUBROUTINE ENTH_H( TX,XHAX,XHOX,XHWX,HGAX,HGOX,HHX,HLWX )
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
!     Handa, Y.P., 1986.  "Compositions, enthalpies of dissociation,
!     and heat capacities in the range 85 to 270 K for clathrate
!     hydrates of methane, ethane, and propane, and enthalpy of 
!     dissociation of isobutane hydrate as determined by a heat-flow
!     calorimeter.  J. Chem. Thermodynamics.  18:915-921.
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 CPX(4)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CPX
      DATA CPX /6.6D+0,1.4538D+0,-3.64D-3,6.312D-6/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ENTH_H'
      TKX = TX+TABS
!
!---  Hydrate mole fractions  ---
!
      WTMX = 1.D+0/(XHAX/WTMA + XHOX/WTMO + XHWX/WTMW)
      XMHAX = (XHAX/WTMA)*WTMX
      XMHOX = (XHOX/WTMO)*WTMX
      XMHWX = (XHWX/WTMW)*WTMX
!
!---  Hydrate mole fractions of formers  ---
!
      YMHOX = XMHOX/(XMHOX+XMHAX+SMALL) 
      YMHAX = XMHAX/(XMHOX+XMHAX+SMALL)
!
!---  Hydrate molecular weight, where 1 mole of hydrate
!     is MnH20, where M is the guest molecule and
!     n is the hydration number  ---
!
      WTMHX = YMHOX*WTMO + YMHAX*WTMA + (XMHWX/(1.D+0-XMHWX))*WTMW
!!
!!---  Enthalpy of liquid water at the reference dissociation
!!     state (i.e., 273.15 K, 101325 Pa)  ---
!!
!      HHX = XHWX*1.D+1
!!
!!---  Enthalpy of gaseous CO2 at the reference dissociation
!!     state (i.e., 273.15 K, 101325 Pa)  ---
!!
!      HHX = HHX + XHAX*(-21.898910D+3)
!!
!!---  Enthalpy of gaseous CH4 at the reference dissociation
!!     state (i.e., 273.15 K, 101325 Pa)  ---
!!
!      HHX = HHX + XHOX*(-56.103140D+3)
!!
!!---  Heat of dissociation from hydrate to liquid water and
!!     gas at 273.15 K, from Yoon (J/kg)  ---
!!
!      HHX = HHX + (-1.D+6)*(YMHAX*57.66D+0 + YMHOX*53.81D+0)/WTMHX
!!
!!---  Specific heat contribution for hydrate from 273.15 K,
!!     from Handa (J/kg)  ---
!!
!      HHX = HHX + 1.D+3*(CPX(1)*TKX + CPX(2)*(TKX**2)/2.D+0 +
!     &  CPX(3)*(TKX**3)/3.D+0 + CPX(4)*(TKX**4)/4.D+0 -
!     &  CPX(1)*TABS - CPX(2)*(TABS**2)/2.D+0 -
!     &  CPX(3)*(TABS**3)/3.D+0 - CPX(4)*(TABS**4)/4.D+0)/WTMHX
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
!---  Heat of dissociation for pure CO2 hydrate to liquid water
!     and gaseous CO2 as a function of temperature
!     from Kamath, J/mol hydrate, corrected to Yoon value
!     at 273.15 of 57.66 KJ/mol, n = 6.21 ---
!
!      HDAX = 4.184D+0*(19.199D+3 - 14.95D+0*TKX)
      HDAX = 4.184D+0*(17.86466D+3 - 14.95D+0*TKX)
!
!---  Heat of dissociation for pure CH4 hydrate to liquid water
!     and gaseous CH4 as a function of temperature
!     from Kamath J/mol hydrate, corrected to Yoon value
!     at 273.15 of 53.81 KJ/mol, n = 6.07 ---
!
!      HDOX = 4.184D+0*(13.521D+3 - 4.02D+0*TKX)
      HDOX = 4.184D+0*(13.95896D+3 - 4.02D+0*TKX)
!
!---  Heat of dissociation from hydrate to liquid water and gaseous
!     CO2 and CH4, J/kg hydrate  ---
!
      HHX = HHX + (-1.D+3)*(YMHAX*HDAX + YMHOX*HDOX)/WTMHX
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
      SAVE A
      DATA A / -6.36443D+5,1.14468D+2,3.64104D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ENTH_I'
      TKX = TX+TABS
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
      SUBROUTINE ENTH_L( TX,XLSX,XLAX,XLOX,HBX,HGAX,HGOX,HLX )
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
!     HBX - enthalpy of brine (NaCl-H2O), J/kg
!     HGAX - enthalpy of gaseous CO2, J/kg
!     HGOX - enthalpy of gaseous CO2, J/kg
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
      TY = TX + DTX
      CALL HC_LA( TY,XLSX,HCY,DHCY )
      DHCY = LOG(HCY/HCX)/DTX
!
!---  Heat of solution for CO2 from Himmelblau (1959)  ---
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
      TY = TX + DTX
      CALL HC_LO( TY,XLSX,HCY,DHCY )
      DHCY = LOG(HCY/HCX)/DTX
!
!---  Heat of solution for CH4 from Himmelblau (1959)  ---
!
      TKX = TX + TABS
      HSOX = -RCU*(TKX**2)*DHCX/WTMO
!
!---  Aqueous enthalpy  ---
!
      XBAX = XLAX/(1.D+0-XLSX)
      XBOX = XLOX/(1.D+0-XLSX)
      HLX = (1.D+0-XBAX-XBOX)*HBX
     &  + XBAX*(HGAX + HSAX) + XBOX*(HGOX + HSOX)
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
      SUBROUTINE ENTH_O( TX,PX,HOX,UOX )
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
!     Enthalpy and internal energy of CH4 gas.
!
!     Setzmann, U., and W. Wagner.  1991.  A New Equation of State and
!     Tables of Thermodynamic Properties for Methane Covering the
!     Range from the Melting Line to 625 K at Pressures up to 1000 MPa.
!     J. Phys. Chem. Ref. Data 20(6):1061-1155.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 2 December 2004.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
      USE CCP
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
      SUB_LOG(ISUB_LOG) = '/ENTH_O'
      INCG = 2
      P_TPX = GCPP(17,INCG)
      P_CRX = GCPP(4,INCG)
!
!---  Table lookup and bilinear interpolation  ---
!
      I_VX = 1
      PTX = MAX( PX,P_TA(1,INCG)*1.D+6 )
      CALL ITL_A( PTX,P_TPX,P_CRX,TX,I_PX,I_SX,I_TX,I_VX )
      CALL PTL_A( PTX,TX,H_TA,HOX,DHOX,I_PX,I_SX,I_TX,I_VX )
      CALL PTL_A( PTX,TX,U_TA,UOX,DUOX,I_PX,I_SX,I_TX,I_VX )
!
!---  Convert units to J/kg  ---
!
      HOX = 1.D+3*HOX
      UOX = 1.D+3*UOX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ENTH_O group  ---
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
      SAVE CAX
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
!
!----------------------Data Statements---------------------------------!
!
      SAVE CAX,SAX,CBX,SBX
      SAVE LFCX,SLX,INX,IZX,ITX,IXX
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
        INDX = 14
        CHMSG = 'Steam Table Subregion 3/5: Exceeds Temperature Range: '
        RLMSG = TX
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Subregions 4 or 5  ---
!
      IF( ISRX.EQ.4 .OR. ISRX.EQ.5 ) THEN
        INDX = 14
        CHMSG = 'Steam Table Subregion 4/5: Exceeds Temperature Range: '
        RLMSG = TX
        CALL WRMSGS( INDX )
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
      SUBROUTINE FUG( AX,AIX,BX,FUGX,PX,RHOMX,TX )
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
!     Written by MD White, PNNL, 24 August 2009
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
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
      REAL*8 AIX(LNGC),FUGX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FUG'
      IF( AX.GT.EPSL ) THEN
        D1X = SQRT(2.D+0) + 1.D+0
        D2X = 1.D+0 - SQRT(2.D+0)
        ALPHAX = SQRT(AX)
        TKX = TX + TABS
        VX = 1.D+0/RHOMX
        Q0X = -LOG((PX/(RCU*TKX))*(VX-BX))
        Q1X = -2.D+0*(ALPHAX/(BX*RCU*TKX*(D1X-D2X)))
     &    *LOG((VX+D1X*BX)/(VX+D2X*BX))
        Q2X = (((PX*VX/(RCU*TKX))-1.D+0)/BX) - (ALPHAX*Q1X/(2.D+0*BX))
        DO 100 IGC = 1,NGC
          FUGX(IGC) = EXP( Q0X + Q1X*SQRT(AIX(IGC)) + Q2X*GCPP(14,IGC) )
  100   CONTINUE
      ELSE
        DO 110 IGC = 1,NGC
          FUGX(IGC) = 0.D+0
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
      REAL*8 COBX(6),COCX(5)
      REAL*8 CHBX(4),CHCX(4)
!
!----------------------Data Statements---------------------------------!
!
      SAVE COBX,COCX
      DATA COBX / 7.83666D+7, 1.96025D+6, 8.20574D+4, -7.40674D+2,
     &  2.18380D+0, -2.20999D-3 /
      DATA COCX / 1.19784D-1, -7.17823D-4, 4.93854D-6, -1.03826D-8,
     &  1.08233D-11 /
      SAVE CHBX,CHCX
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
      HCX = HCX*CHKN(1)
      DHCX = DHCX*CHKN(1)
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
      REAL*8 CHBX(4),CHCX(4)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CHBX,CHCX
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
      HCX = HCX*CHKN(1)
      DHCX = DHCX*CHKN(1)
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
      SUBROUTINE HYD_P( TX,YMGOX,XLSX,PHX,RHOHX,XHAX,XHOX,XHWX )
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
!     Compute hydrate equilibrium pressure and hydrate-H2O, -CO2, and
!     -CH4 mass fractions
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 August 2004.
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
      REAL*8 UVWX(8),WTMX(LNNGC),FOCX(2,LNNGC)
      REAL*8 TCPWX(2)
      INTEGER I_AX(2),I_TX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/HYD_P'
!
!---  Limit temperature and mass fraction of hydrate formers  ---
!
      TY = MAX( TX,T_PH(1) )
      TY = MIN( TY,T_PH(LT_PH) )
      YMGOX = MAX( YMGOX,YMGO_PH(1) )
      YMGOX = MIN( YMGOX,YMGO_PH(LO_PH) )
!
!---  Find Lower Temperature Index
!
      I_TL = 1
      I_TU = LT_PH
  110 IF( I_TU-I_TL.GT.1 ) THEN
        I_TM = (I_TL+I_TU)/2
        IF( (T_PH(LT_PH).GT.T_PH(1)).EQV.
     &    (TY.GT.T_PH(I_TM)) ) THEN
          I_TL = I_TM
        ELSE
          I_TU = I_TM
        ENDIF
        GOTO 110
      ENDIF
      I_TX(1) = I_TL
      I_TX(2) = I_TU
!
!---  Find Lower CH4 Gas Mole Fraction Index
!
      I_AL = 1
      I_AU = LO_PH
  120 IF( I_AU-I_AL.GT.1 ) THEN
        I_AM = (I_AL+I_AU)/2
        IF( (YMGO_PH(LO_PH).GT.YMGO_PH(1)).EQV.
     &    (YMGOX.GT.YMGO_PH(I_AM)) ) THEN
          I_AL = I_AM
        ELSE
          I_AU = I_AM
        ENDIF
        GOTO 120
      ENDIF
      I_AX(1) = I_AL
      I_AX(2) = I_AU
!
!---  Bilinear interpolation scheme  ---
!
      UX = ( TY-T_PH(I_TX(1)) )/( T_PH(I_TX(2))-T_PH(I_TX(1)) )
      WX = ( YMGOX-YMGO_PH(I_AX(1)) )/
     &  ( YMGO_PH(I_AX(2))-YMGO_PH(I_AX(1)) )
      UVWX(1) = (1.D+0-UX)*(1.D+0-WX)
      UVWX(2) = UX*(1.D+0-WX)
      UVWX(3) = (1.D+0-UX)*WX
      UVWX(4) = UX*WX
!
!---  Pressure  ---
!
      IC = 0
      PHX = 0.D+0
      DO 134 IAX = I_AX(1),I_AX(2)
      DO 130 ITX = I_TX(1),I_TX(2)
        IC = IC + 1
        PHX = PHX + UVWX(IC)*P_PH(ITX,IAX)
  130 CONTINUE
  134 CONTINUE
!
!---  Inhibitor correction  ---
!
      IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
        CALL HYD_PO( YMGOX,POHX )
      ENDIF
!
!---  Converge on pressure using Newton-Raphson iteration  ---
!
      TEQHY = TY
      TEQHZ = TY
      NC = 0
  140 CONTINUE
      NC = NC + 1
      PEQHY = PHX
      CALL HYD_T( TEQHY,YMGOX,XLSX,PEQHY,POHX,RHOHX,XHAX,XHOX,XHWX )
      PEQHZ = PHX + 1.D-4
      CALL HYD_T( TEQHZ,YMGOX,XLSX,PEQHZ,POHX,RHOHX,XHAX,XHOX,XHWX )
      FX = TY - TEQHY
      IF( ABS(FX).GT.EPSL ) THEN
        FZ = TY - TEQHZ
        DFX = (FZ-FX)*1.D+4
        DPHX = -FX/DFX
        PHX = PHX + DPHX
        IF( NC.GT.32 ) GOTO 144
        IF( ABS(DPHX).GT.1.D-4 ) GOTO 140
  144   CONTINUE
      ENDIF
      CALL HYD_T( TEQHX,YMGOX,XLSX,PHX,POHX,RHOHX,XHAX,XHOX,XHWX )
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
      SUBROUTINE HYD_PO( YMGOX,PHX )
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
!     Hydrate equilibrium pressure at 0 C and zero inhibitor
!     concentration as a function of the CH4 gas mole fraction
!     of formers.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 22 April 2010.
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
      REAL*8 UVWX(8),WTMX(LNNGC),FOCX(2,LNNGC)
      REAL*8 TCPWX(2)
      INTEGER I_AX(2),I_TX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/HYD_PO'
!
!---  Limit temperature and mass fraction of hydrate formers  ---
!
      TX = 0.D+0
      XLSX = 0.D+0
      TY = MAX( TX,T_PH(1) )
      TY = MIN( TY,T_PH(LT_PH) )
      YMGOX = MAX( YMGOX,YMGO_PH(1) )
      YMGOX = MIN( YMGOX,YMGO_PH(LO_PH) )
!
!---  Find Lower Temperature Index
!
      I_TL = 1
      I_TU = LT_PH
  110 IF( I_TU-I_TL.GT.1 ) THEN
        I_TM = (I_TL+I_TU)/2
        IF( (T_PH(LT_PH).GT.T_PH(1)).EQV.
     &    (TY.GT.T_PH(I_TM)) ) THEN
          I_TL = I_TM
        ELSE
          I_TU = I_TM
        ENDIF
        GOTO 110
      ENDIF
      I_TX(1) = I_TL
      I_TX(2) = I_TU
!
!---  Find Lower CH4 Gas Mole Fraction Index
!
      I_AL = 1
      I_AU = LO_PH
  120 IF( I_AU-I_AL.GT.1 ) THEN
        I_AM = (I_AL+I_AU)/2
        IF( (YMGO_PH(LO_PH).GT.YMGO_PH(1)).EQV.
     &    (YMGOX.GT.YMGO_PH(I_AM)) ) THEN
          I_AL = I_AM
        ELSE
          I_AU = I_AM
        ENDIF
        GOTO 120
      ENDIF
      I_AX(1) = I_AL
      I_AX(2) = I_AU
!
!---  Bilinear interpolation scheme  ---
!
      UX = ( TY-T_PH(I_TX(1)) )/( T_PH(I_TX(2))-T_PH(I_TX(1)) )
      WX = ( YMGOX-YMGO_PH(I_AX(1)) )/
     &  ( YMGO_PH(I_AX(2))-YMGO_PH(I_AX(1)) )
      UVWX(1) = (1.D+0-UX)*(1.D+0-WX)
      UVWX(2) = UX*(1.D+0-WX)
      UVWX(3) = (1.D+0-UX)*WX
      UVWX(4) = UX*WX
!
!---  Pressure  ---
!
      IC = 0
      PHX = 0.D+0
      DO 134 IAX = I_AX(1),I_AX(2)
      DO 130 ITX = I_TX(1),I_TX(2)
        IC = IC + 1
        PHX = PHX + UVWX(IC)*P_PH(ITX,IAX)
  130 CONTINUE
  134 CONTINUE
!
!---  Converge on pressure using Newton-Raphson iteration  ---
!
      TEQHY = TY
      TEQHZ = TY
      POHX = 0.D+0
      NC = 0
  140 CONTINUE
      NC = NC + 1
      PEQHY = PHX
      CALL HYD_T( TEQHY,YMGOX,XLSX,PEQHY,POHX,RHOHX,XHAX,XHOX,XHWX )
      PEQHZ = PHX + 1.D-4
      CALL HYD_T( TEQHZ,YMGOX,XLSX,PEQHZ,POHX,RHOHX,XHAX,XHOX,XHWX )
      FX = TY - TEQHY
      IF( ABS(FX).GT.EPSL ) THEN
        FZ = TY - TEQHZ
        DFX = (FZ-FX)*1.D+4
        DPHX = -FX/DFX
        PHX = PHX + DPHX
        IF( NC.GT.32 ) GOTO 144
        IF( ABS(DPHX).GT.1.D-4 ) GOTO 140
  144   CONTINUE
      ENDIF
      CALL HYD_T( TEQHX,YMGOX,XLSX,PHX,POHX,RHOHX,XHAX,XHOX,XHWX )
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of HYD_PO group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE HYD_T( TX,YMGOX,XLSX,PHX,POHX,RHOHX,XHAX,XHOX,XHWX )
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
!     Compute hydrate equilibrium temperature and hydrate-H2O, -CO2, 
!     and -CH4 mass fractions
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 13 August 2004.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
      USE NAPL
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 UVWX(8),WTMX(LNNGC),FOCX(2,LNNGC)
      REAL*8 TCPWX(2)
      INTEGER I_PX(2),I_AX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/HYD_T'
!
!---  Check input pressure, CO2 gas mole fraction,
!     and dissolved-NaCl aqueous mass fraction  ---
!
      IF( PHX.LT.P_TH(1) ) THEN
        TX = -1.D+20
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      IF( PHX.GT.P_TH(LT_TH) ) THEN
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      IF( YMGOX.LT.YMGO_TH(1) .OR. YMGOX.GT.YMGO_TH(LO_TH) ) THEN
        INDX = 17
        CHMSG = 'CH4 Gas Mole Fraction (Excluding H2O) Out of Range: '
        RLMSG = YMGOX
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Find Lower Pressure Index
!
      I_PL = 1
      I_PU = LT_TH
  110 IF( I_PU-I_PL.GT.1 ) THEN
        I_PM = (I_PL+I_PU)/2
        IF( (P_TH(LT_TH).GT.P_TH(1)).EQV.
     &    (PHX.GT.P_TH(I_PM)) ) THEN
          I_PL = I_PM
        ELSE
          I_PU = I_PM
        ENDIF
        GOTO 110
      ENDIF
      I_PX(1) = I_PL
      I_PX(2) = I_PU
!
!---  Find Lower CH4 Gas Mole Fraction Index
!
      I_AL = 1
      I_AU = LO_TH
  130 IF( I_AU-I_AL.GT.1 ) THEN
        I_AM = (I_AL+I_AU)/2
        IF( (YMGO_TH(LO_TH).GT.YMGO_TH(1)).EQV.
     &    (YMGOX.GT.YMGO_TH(I_AM)) ) THEN
          I_AL = I_AM
        ELSE
          I_AU = I_AM
        ENDIF
        GOTO 130
      ENDIF
      I_AX(1) = I_AL
      I_AX(2) = I_AU
!
!---  Bilinear interpolation scheme  ---
!
      UX = ( PHX-P_TH(I_PX(1)) )/( P_TH(I_PX(2))-P_TH(I_PX(1)) )
      WX = ( YMGOX-YMGO_TH(I_AX(1)) )/
     &  ( YMGO_TH(I_AX(2))-YMGO_TH(I_AX(1)) )
      UVWX(1) = (1.D+0-UX)*(1.D+0-WX)
      UVWX(2) = UX*(1.D+0-WX)
      UVWX(3) = (1.D+0-UX)*WX
      UVWX(4) = UX*WX
!
!---  Temperature  ---
!
      IC = 0
      TX = 0.D+0
      DO 144 IAX = I_AX(1),I_AX(2)
      DO 140 IPX = I_PX(1),I_PX(2)
        IC = IC + 1
        TX = TX + UVWX(IC)*T_TH(IPX,IAX)
  140 CONTINUE
  144 CONTINUE
!
!---  Inhibitor correction  ---
!
      DTX = 0.D+0
      IF( I_INH.NE.0 .AND. ISLC(41).EQ.1 .AND. XLSX.GT.EPSL ) THEN
        XLSFX = 1.D+2*XLSX
        POX = 1.D-3*(POHX + PATM)
        DTX = (CINH(1,I_INH)*XLSFX + CINH(2,I_INH)*(XLSFX**2)
     &    + CINH(3,I_INH)*(XLSFX**3))*
     &    (CINH(4,I_INH)*LOG(PHX*1.D-3) + CINH(5,I_INH))*
     &    (CINH(6,I_INH)*(POX-1.D+3)+1.D+0)
      ENDIF
      TX = TX - DTX
!
!---  CH4 Hydrate mole fraction (excluding water)  ---
!
      IC = 0
      YMHOX = 0.D+0
      DO 154 IAX = I_AX(1),I_AX(2)
      DO 150 IPX = I_PX(1),I_PX(2)
        IC = IC + 1
        YMHOX = YMHOX + UVWX(IC)*YMHO_TH(IPX,IAX)
  150 CONTINUE
  154 CONTINUE
!
!---  CO2 Hydrate mole fraction (excluding water)  ---
!
      YMHAX = MAX( 1.D+0-YMHOX,0.D+0 )
!
!---  CO2 small-cage occupancy number  ---
!
      IC = 0
      XSCAX = 0.D+0
      DO 164 IAX = I_AX(1),I_AX(2)
      DO 160 IPX = I_PX(1),I_PX(2)
        IC = IC + 1
        XSCAX = XSCAX + UVWX(IC)*XSCA_TH(IPX,IAX)
  160 CONTINUE
  164 CONTINUE
!
!---  CO2 large-cage occupancy number  ---
!
      IC = 0
      XLCAX = 0.D+0
      DO 174 IAX = I_AX(1),I_AX(2)
      DO 170 IPX = I_PX(1),I_PX(2)
        IC = IC + 1
        XLCAX = XLCAX + UVWX(IC)*XLCA_TH(IPX,IAX)
  170 CONTINUE
  174 CONTINUE
!
!---  CH4 small-cage occupancy number  ---
!
      IC = 0
      XSCOX = 0.D+0
      DO 184 IAX = I_AX(1),I_AX(2)
      DO 180 IPX = I_PX(1),I_PX(2)
        IC = IC + 1
        XSCOX = XSCOX + UVWX(IC)*XSCO_TH(IPX,IAX)
  180 CONTINUE
  184 CONTINUE
!
!---  CH4 large-cage occupancy number  ---
!
      IC = 0
      XLCOX = 0.D+0
      DO 194 IAX = I_AX(1),I_AX(2)
      DO 190 IPX = I_PX(1),I_PX(2)
        IC = IC + 1
        XLCOX = XLCOX + UVWX(IC)*XLCO_TH(IPX,IAX)
  190 CONTINUE
  194 CONTINUE
!
!---  H2O hydrate mole fraction  ---
!
      IF( XLCAX+XLCOX+XSCOX+XSCAX.GT.EPSL ) THEN
        XMHWX = 23.D+0/((3.D+0*(XLCAX+XLCOX))+XSCOX+XSCAX)
        XMHWX = XMHWX/(1.D+0+XMHWX)
        XMHOX = (1.D+0-XMHWX)*YMHOX
        XMHAX = MAX(1.D+0-XMHWX-XMHOX,0.D+0)
      ELSE
        XMHWX = 1.D+0
        XMHOX = 0.D+0
        XMHAX = 0.D+0
      ENDIF
!
!---  Hydrate mass fractions  ---
!
      WTMHX = XMHWX*WTMW + XMHAX*WTMA + XMHOX*WTMO
      XHWX = XMHWX*WTMW/WTMHX
      XHAX = XMHAX*WTMA/WTMHX
      XHOX = XMHOX*WTMO/WTMHX
!
!---  Hydrate density, assuming sI hydrate crystal cell structure,
!     
!     WMPUC: number of water molecules per unit cell
!     AVGNX: Avagadro's number, molecules/mol
!     WMTMX: molecular weight of component, g/mol
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
      WTMX(1) = WTMA
      WTMX(2) = WTMO
      FOCX(1,1) = XSCAX
      FOCX(2,1) = XLCAX
      FOCX(1,2) = XSCOX
      FOCX(2,2) = XLCOX
      TCPWX(1) = 2.D+0
      TCPWX(2) = 6.D+0
      VUCX = 12.D-8**3
      NCTX = 2
      NCHX = 2
      RHOHX = WMPUC*WTMW
      DO 210 J = 1,NCHX
      DO 200 I = 1,NCTX
        RHOHX = RHOHX + FOCX(I,J)*TCPWX(I)*WTMX(J)
  200 CONTINUE
  210 CONTINUE
      RHOHX = 1.D+3*RHOHX/(AVGNX*VUCX)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of HYD_T group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ICE( TX,PLX,XLAX,XLOX,XLSX,TFPX,PIX )
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
!        RADX = 3.6D-9*(1.D+0 - (6.6D-1/LOG10((TX+TABS)/(TFPX+TABS))))
        RADX = 1.0000D-7/(TFPX-TX)
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

!!----------------------Subroutine--------------------------------------!
!!
!      SUBROUTINE ISO_FUG
!!
!!-------------------------Disclaimer-----------------------------------!
!!
!!     This material was prepared as an account of work sponsored by
!!     an agency of the United States Government. Neither the
!!     United States Government nor the United States Department of
!!     Energy, nor Battelle, nor any of their employees, makes any
!!     warranty, express or implied, or assumes any legal liability or
!!     responsibility for the accuracy, completeness, or usefulness
!!     of any information, apparatus, product, software or process
!!     disclosed, or represents that its use would not infringe
!!     privately owned rights.
!!
!!----------------------Acknowledgement---------------------------------!
!!
!!     This software and its documentation were produced with Government
!!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!!     United Department of Energy. The Government retains a paid-up
!!     non-exclusive, irrevocable worldwide license to reproduce,
!!     prepare derivative works, perform publicly and display publicly
!!     by or for the Government, including the right to distribute to
!!     other Government contractors.
!!
!!---------------------Copyright Notices--------------------------------!
!!
!!            Copyright Battelle Memorial Institute, 1996
!!                    All Rights Reserved.
!!
!!----------------------Description-------------------------------------!
!!
!!     Generate data for the fugacity coefficient as a function of
!!     pressure at constant temperature, dissolved NaCl, and
!!     dissolved CO2.
!!
!!----------------------Authors-----------------------------------------!
!!
!!     Written by MD White, PNNL, 27 July 2002.
!!
!!----------------------Fortran 90 Modules------------------------------!
!!
!      USE GLB_PAR
!      USE SOLTN
!      USE CONST
!      USE CCP
!
!!----------------------Implicit Double Precision-----------------------!
!!
!      IMPLICIT REAL*8 (A-H,O-Z)
!      IMPLICIT INTEGER (I-N)
!!
!!----------------------Type Declarations-------------------------------!
!!
!      REAL*8 PVAX(LP_TA),FUGAX(LP_TA),FUGAX2(LP_TA)
!      INTEGER I_PX(2),I_SX,I_TX(2)
!!
!!----------------------Executable Lines--------------------------------!
!!
!      ISUB_LOG = ISUB_LOG+1
!      SUB_LOG(ISUB_LOG) = '/ISO_FUG'
!!
!!---  Assign values and open file  ---
!!
!      TX = 25.D+0
!      P_TPX = GCPP(17,INCG)
!      P_CRX = GCPP(4,INCG)
!      OPEN(UNIT=26, FILE='iso_fug', STATUS='UNKNOWN', FORM='FORMATTED')
!      CLOSE(UNIT=26,STATUS='DELETE')
!      OPEN(UNIT=26, FILE='iso_fug', STATUS='NEW', FORM='FORMATTED')
!      DO 100 M = 1,IP_TA(INCG)
!        I_VX = 1
!        PVAX(M) = P_TA(M,INCG)*1.D+6
!        CALL ITL_A( PVAX(M),P_TPX,P_CRX,TX,I_PX,I_SX,I_TX,I_VX )
!        CALL PTL_A( PVAX(M),TX,FUG_TA,FUGAX(M),DFUGAX,I_PX,I_SX,
!     &    I_TX,I_VX )
!  100 CONTINUE
!      CALL SPLINE( PVAX,FUGAX,IP_TA(INCG),FUGAX2 )
!      DO 200 M = 2,IP_TA(INCG)
!      DO 200 N = 1,4
!        REALX = REAL(N-1)
!        PVX = (REALX/4.D+0)*(PVAX(M)-PVAX(M-1)) + PVAX(M-1)
!        CALL SPLINT( PVAX,FUGAX,FUGAX1,FUGAX2,IP_TA(INCG),PVX,FUGX )
!        WRITE(26,'(3(1PE11.4,1X))') PVX,FUGX,FUGAX1
!  200 CONTINUE
!      CLOSE(UNIT=26)
!!
!!---  Reset subroutine string sequence  ---
!!
!      ISUB_LOG = ISUB_LOG-1
!!
!!---  End of ISO_FUG group  ---
!!
!      RETURN
!      END

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
!     Written by MD White, PNNL, 22 April 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
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
!     search the across the table for gas
!     and liquid ---
!
      IF( PMX.LT.P_TA(IPTPX,INCG) ) THEN
        ISLX = 1
        IELX = IT_TA(IPLX,INCG)
        ISUX = 1
        IEUX = IT_TA(IPUX,INCG)
        I_VX = 0
!
!---  Pressure above the critical-point pressure,
!     search the across the table for gas
!     and liquid ---
!
      ELSEIF( PMX.GT.P_TA(IPCRX,INCG) ) THEN
        ISLX = 1
        IELX = IT_TA(IPLX,INCG)
        ISUX = 1
        IEUX = IT_TA(IPUX,INCG)
        I_VX = 0
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
        IF( PMX.GE.PM_CRX ) TSX = GCPP(3,INCG)
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
          IF( INCG.EQ.1 ) THEN
            I_VX = 1
          ELSEIF( INCG.EQ.2 ) THEN
            I_VX = 1
          ENDIF
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
      SUBROUTINE PERM_R( SSX,PERMRFX,N )
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
!     Verma, A., and K. Pruess.  1988.  "Thermohydrological Conditions
!     and Silica Redistribution Near High-Level Nuclear Wastes
!     Emplaced in Saturated Geological Formations."  Journal of
!     Geophysical Research, 93(B2):1159-1173.
!
!     Pruess, K., and J. Garcia.  2002.  "Multiphase flow dynamics
!     during CO2 disposal into saline aquifers."  Environmental Geology
!     http://link.springer.de/link/service/journals/00254/contents
!     /01/00498/paper/s00254-001-0498-3ch110.html
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 May 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
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
      IZN = IZ(N)
!
!---  Normalized porosity  ---
!
      PORD_NX = (1.D+0-SSX-PERM(5,IZN))/(1.D+0-PERM(5,IZN))
!
!---  Tube area ratio  ---
!
      OMEGAX = 1.D+0 + (1.D+0/PERM(4,IZN))/((1.D+0/PERM(5,IZN))-1.D+0)
!
!---  Permeability reduction factor  ---
!
      PERMRFX = (PORD_NX**2)*(1.D+0-PERM(4,IZN)+PERM(4,IZN)/(OMEGAX**2))
     &  /(1.D+0-PERM(4,IZN)+PERM(4,IZN)*
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
!     I_VX = 0 - bilinear interpolation
!     I_VX = 1 - liquid-side interpolation
!     VAR_TA(LT_TA,LP_TA,LNNGC) - pressure, temperature ordered 
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
      SUBROUTINE PRC_T( AIX,TX )
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
      REAL*8 AIX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PRC_T'
      TKX = TX + TABS
      DO 100 IGC = 1,NGC
        TRX = TKX/GCPP(3,IGC)
        AIX(IGC) = GCPP(13,IGC)*
     &    (1.D+0 + GCPP(15,IGC)*(1.D+0 - SQRT(TRX)))**2
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
      SUBROUTINE PRC_X( AX,ALPHAX,AIX,BX,WTMX,XMCX )
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
!     Written by MD White, PNNL, 24 August 2009
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 AIX(LNGC),XMCX(LNGC)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PRC_X'
      AX = 0.D+0
      BX = 0.D+0
      WTMX = 0.D+0
      DO 120 IGC = 1,NGC
        DO 110 JGC = 1,NGC
          AIJX = SQRT(AIX(IGC)*AIX(JGC))
          AX = AX + XMCX(IGC)*XMCX(JGC)*AIJX
  110   CONTINUE
        ALPHAX = SQRT(AX)
        BX = BX + XMCX(IGC)*GCPP(14,IGC)
        WTMX = WTMX + XMCX(IGC)*GCPP(1,IGC)
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
      SUBROUTINE RDPF_HYD_KE
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
      LOGICAL FCHK
      CHARACTER*64 FDUM
      CHARACTER*64 FNM(LNNGC)
      CHARACTER*512 CHDUM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDPF_HYD_KE'
!
!---  Load file names  ---
!
      FNM(1) = 'co2_prop.dat'
      FNM(2) = 'ch4_prop.dat'
!
!---  Loop over noncondensible gases  ---
!
      DO 300 INCG = 1,2
!
!---    Check for property file  ---
!
        NCH = INDEX( FNM(INCG)(1:),'  ' )-1
        INQUIRE( FILE=FNM(INCG)(1:NCH), FORM=FDUM, EXIST=FCHK )
        IF( .NOT.FCHK ) THEN
          INDX = 4
          CHMSG = 'Missing Property File: ' // FNM(INCG)(1:NCH)
          CALL WRMSGS( INDX )
        ELSEIF( FDUM.EQ.'unformatted' ) THEN
          INDX = 4
          CHMSG = 'Unformatted Property File: ' // FNM(INCG)(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
        OPEN(UNIT=26, FILE=FNM(INCG)(1:NCH), STATUS='OLD', 
     &    FORM='FORMATTED')
!
!---    Read Id tag line ---
!
        READ(26,'(A)') CHDUM
!
!---    Read property file  ---
!
        READ(26,*) IP_TA(INCG)
        READ(26,*) (P_TA(I,INCG),I=1,IP_TA(INCG))
        READ(26,*) (IT_TA(I,INCG),I=1,IP_TA(INCG))
        DO 10 IPX = 1,IP_TA(INCG)
          READ(26,*) (T_TA(ITX,IPX,INCG),ITX=1,IT_TA(IPX,INCG))
   10   CONTINUE
        DO 20 IPX = 1,IP_TA(INCG)
          READ(26,*) (RHO_TA(ITX,IPX,INCG),ITX=1,IT_TA(IPX,INCG))
   20   CONTINUE
        DO 30 IPX = 1,IP_TA(INCG)
          READ(26,*) (H_TA(ITX,IPX,INCG),ITX=1,IT_TA(IPX,INCG))
   30  CONTINUE
        DO 40 IPX = 1,IP_TA(INCG)
          READ(26,*) (U_TA(ITX,IPX,INCG),ITX=1,IT_TA(IPX,INCG))
   40   CONTINUE
        DO 50 IPX = 1,IP_TA(INCG)
          READ(26,*) (FUG_TA(ITX,IPX,INCG),ITX=1,IT_TA(IPX,INCG))
   50   CONTINUE
        DO 60 IPX = 1,IP_TA(INCG)
          READ(26,*) (S_TA(ITX,IPX,INCG),ITX=1,IT_TA(IPX,INCG))
   60   CONTINUE
        READ(26,*) (IV_TA(I,INCG),I=1,IP_TA(INCG))
        READ(26,*) I_LV(INCG)
        DO 70 IPX = 1,I_LV(INCG)
          READ(26,*) T_LV(IPX,INCG),P_LV(IPX,INCG),RHOL_LV(IPX,INCG),
     &      HL_LV(IPX,INCG),UL_LV(IPX,INCG),SL_LV(IPX,INCG),
     &      RHOV_LV(IPX,INCG),HV_LV(IPX,INCG),UV_LV(IPX,INCG),
     &      SV_LV(IPX,INCG),FUG_LV(IPX,INCG)
   70   CONTINUE
        CLOSE(26)
!
!---  Load pressure triple-point and critical-point indices  ---
!
      DO 80 I=1,IP_TA(INCG)
        PTPAX = 1.D-6*GCPP(17,INCG)
        PCRAX = 1.D-6*GCPP(4,INCG)
        IF( ABS(P_TA(I,INCG)-PTPAX).LT.EPSL ) IPTP(INCG) = I
        IF( ABS(P_TA(I,INCG)-PCRAX).LT.EPSL ) IPCR(INCG) = I
   80 CONTINUE
!!
!!---  Vapor temperature splines  ---
!!
!      DO 110 IPX = 1,IP_TA(INCG)
!        N = IT_TA(IPX,INCG)-IV_TA(IPX,INCG)+1
!        CALL SPLINE( T_TA(IV_TA(IPX,INCG),IPX,INCG),
!     &    RHO_TA(IV_TA(IPX,INCG),IPX,INCG),N,
!     &    RHO_ST(IV_TA(IPX,INCG),IPX,INCG) )
!        CALL SPLINE( T_TA(IV_TA(IPX,INCG),IPX,INCG),
!     &    H_TA(IV_TA(IPX,INCG),IPX,INCG),N,
!     &    H_ST(IV_TA(IPX,INCG),IPX,INCG) )
!        CALL SPLINE( T_TA(IV_TA(IPX,INCG),IPX,INCG),
!     &    U_TA(IV_TA(IPX,INCG),IPX,INCG),N,
!     &    U_ST(IV_TA(IPX,INCG),IPX,INCG) )
!        CALL SPLINE( T_TA(IV_TA(IPX,INCG),IPX,INCG),
!     &    FUG_TA(IV_TA(IPX,INCG),IPX,INCG),N,
!     &    FUG_ST(IV_TA(IPX,INCG),IPX,INCG) )
!        CALL SPLINE( T_TA(IV_TA(IPX,INCG),IPX,INCG),
!     &    S_TA(IV_TA(IPX,INCG),IPX,INCG),N,
!     &    S_ST(IV_TA(IPX,INCG),IPX,INCG) )
!  110 CONTINUE
!!
!!---  Liquid temperature splines  ---
!!
!      DO 210 IPX = 1,IP_TA(INCG)
!        N = IV_TA(IPX,INCG)-1
!        IF( N.GT.0 ) THEN
!          CALL SPLINE( T_TA(1,IPX,INCG),
!     &      RHO_TA(1,IPX,INCG),N,
!     &      RHO_ST(1,IPX,INCG) )
!          CALL SPLINE( T_TA(1,IPX,INCG),
!     &      H_TA(1,IPX,INCG),N,
!     &      H_ST(1,IPX,INCG) )
!          CALL SPLINE( T_TA(1,IPX,INCG),
!     &      U_TA(1,IPX,INCG),N,
!     &      U_ST(1,IPX,INCG) )
!          CALL SPLINE( T_TA(1,IPX,INCG),
!     &      FUG_TA(1,IPX,INCG),N,
!     &      FUG_ST(1,IPX,INCG) )
!          CALL SPLINE( T_TA(1,IPX,INCG),
!     &      S_TA(1,IPX,INCG),N,
!     &      S_ST(1,IPX,INCG) )
!        ENDIF
!  210 CONTINUE
  300 CONTINUE
!
!---  Check for hydrate equilibrium-pressure file  ---
!
      INQUIRE( FILE='hydrate_peq.dat', FORM=FDUM, EXIST=FCHK )
      IF( .NOT.FCHK ) THEN
        INDX = 4
        CHMSG = 'Missing Hydrate File: hydrate_peq.dat'
        CALL WRMSGS( INDX )
      ELSEIF( FDUM.EQ.'unformatted' ) THEN
        INDX = 4
        CHMSG = 'Unformatted Hydrate File: hydrate_peq.dat'
        CALL WRMSGS( INDX )
      ENDIF
      OPEN(UNIT=26, FILE='hydrate_peq.dat', STATUS='OLD', 
     &  FORM='FORMATTED')
!
!---  Read Id tag line ---
!
      READ(26,'(A)') CHDUM
!
!---  Read hydrate equilibrium-pressure file  ---
!
      READ(26,*) IT_PH
      READ(26,*) IO_PH
      READ(26,*) (T_PH(ITX),ITX=1,IT_PH)
      READ(26,*) (YMGO_PH(IOX),IOX=1,IO_PH)
      DO 140 IOX = 1,IO_PH
      DO 120 ITX = 1,IT_PH
        READ(26,*) P_PH(ITX,IOX),YMHO_PH(ITX,IOX),
     &    XSCO_PH(ITX,IOX),XLCO_PH(ITX,IOX),
     &    XSCA_PH(ITX,IOX),XLCA_PH(ITX,IOX)
  120 CONTINUE
  140 CONTINUE
      CLOSE(26)
!
!---  Check for hydrate equilibrium-temperature file  ---
!
      INQUIRE( FILE='hydrate_teq.dat', FORM=FDUM, EXIST=FCHK )
      IF( .NOT.FCHK ) THEN
        INDX = 4
        CHMSG = 'Missing Hydrate File: hydrate_teq.dat'
        CALL WRMSGS( INDX )
      ELSEIF( FDUM.EQ.'unformatted' ) THEN
        INDX = 4
        CHMSG = 'Unformatted Hydrate File: hydrate_teq.dat'
        CALL WRMSGS( INDX )
      ENDIF
      OPEN(UNIT=26, FILE='hydrate_teq.dat', STATUS='OLD', 
     &  FORM='FORMATTED')
!
!---  Read Id tag line ---
!
      READ(26,'(A)') CHDUM
!
!---  Read hydrate equilibrium-temperature file  ---
!
      READ(26,*) IT_TH
      READ(26,*) IO_TH
      READ(26,*) (P_TH(ITX),ITX=1,IT_TH)
      READ(26,*) (YMGO_TH(IOX),IOX=1,IO_TH)
      DO 170 IOX = 1,IO_TH
      DO 150 ITX = 1,IT_TH
        READ(26,*) T_TH(ITX,IOX),YMHO_TH(ITX,IOX),
     &    XSCO_TH(ITX,IOX),XLCO_TH(ITX,IOX),
     &    XSCA_TH(ITX,IOX),XLCA_TH(ITX,IOX)
  150 CONTINUE
  160 CONTINUE
  170 CONTINUE
      CLOSE(26)
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
!---  Read Id tag line ---
!
      READ(26,'(A)') CHDUM
!
!---  Read inhibitor correlation file  ---
!
      DO 180 INHX = 1,LINH
        READ(26,*) INHNM(INHX),CINH(1,INHX),CINH(2,INHX),CINH(3,INHX),
     &    CINH(4,INHX),CINH(5,INHX),CINH(6,INHX),XLIMX(INHX)
  180 CONTINUE
      CLOSE(26)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDPF_HYD_KE group  ---
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
      SUBROUTINE SOL_LA( TX,PVAX,XLSX,XMLAX )
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
!     H2O-NaCl-CO2 Equation of state.  This subroutine computes the
!     equilibrium mole fraction of dissolved noncondensible gas in the
!     aqueous phase as a function of temperature, pressure, water-vapor
!     pressure and fugacity.  The temperature and salt concentration
!     dependence on Henry's coefficient is computed external to this
!     routine.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York. pp: 332-337.
!
!     Anderson, G.M., and D.A. Crerar.  1992.  Thermodynamics in
!     Geochemistry: The Equilibrium Model, Oxford University Press.
!
!     PMV_NG - partial molar volume (m^3/kmol) of CO2 at infinite
!       dilution in the aqueous phase (pure water value)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 1 April 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
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
      REAL*8 PMV_C(5)
      INTEGER I_PX(2),I_SX,I_TX(2)
!
!----------------------Data Statements---------------------------------!
!
      SAVE PMV_C
      DATA PMV_C / 37.36D-3,-7.109D-5,-3.812D-8,3.296D-9,-3.702D-12 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SOL_LA'
      INCG = 1
      P_TPX = GCPP(17,INCG)
      P_CRX = GCPP(4,INCG)
!
!---  Partial molar volume of CO2 (m^3/kmol), formulation of
!     Anderson et al. (1992)  ---
!
      PMV_A = 0.D+0
      DO 10 M = 1,5
        PMV_A = PMV_A + PMV_C(M)*(TX**(M-1))
   10 CONTINUE
!
!---  Henry's constant for CO2 dissolved in NaCl aqueous solutions  ---
!
      CALL HC_LA( TX,XLSX,HCX,DHCX )
!
!---  Fugacity coefficient for CO2  ---
!
      FUGAX = 1.D+0
      IF( (PVAX*1.D-6).GT.P_TA(1,INCG) ) THEN
        I_VX = 1
        CALL ITL_A( PVAX,P_TPX,P_CRX,TX,I_PX,I_SX,I_TX,I_VX )
        CALL PTL_A( PVAX,TX,FUG_TA,FUGAX,DFUGAX,I_PX,I_SX,I_TX,I_VX )
      ENDIF
      TKX = TX+TABS
!
!---  Poynting effect for CO2 solubility  ---
!
      IF( ISLC(25).EQ.1 ) THEN
        HCX = EXP( LOG(HCX) + (PMV_A*PVAX/(RCU*TKX)) )
      ENDIF
      XMLAX = MIN( MAX( (PVAX*FUGAX/HCX),0.D+0 ),1.D+0 )
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
      SAVE SFX
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
      SUBROUTINE SOL_LO( TX,PVOX,XLSX,XMLOX )
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
!     H2O-NaCl-CO2 Equation of state.  This subroutine computes the
!     equilibrium mole fraction of dissolved noncondensible gas in the
!     aqueous phase as a function of temperature, pressure, water-vapor
!     pressure and fugacity.  The temperature and salt concentration
!     dependence on Henry's coefficient is computed external to this
!     routine.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York. pp: 332-337.
!
!     Anderson, G.M., and D.A. Crerar.  1992.  Thermodynamics in
!     Geochemistry: The Equilibrium Model, Oxford University Press.
!
!     PMV_NG - partial molar volume (m^3/kmol) of CO2 at infinite
!       dilution in the aqueous phase (pure water value)
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 1 April 2002.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NCG_PT
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
      INTEGER I_PX(2),I_SX,I_TX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SOL_LO'
      INCG = 2
      P_TPX = GCPP(17,INCG)
      P_CRX = GCPP(4,INCG)
      PMV_O = 37.D-3
!
!---  Henry's constant for CH4 dissolved in NaCl aqueous solutions  ---
!
      CALL HC_LO( TX,XLSX,HCX,DHCX )
!
!---  Fugacity coefficient for CH4  ---
!
      FUGOX = 1.D+0
      IF( (PVOX*1.D-6).GT.P_TA(1,INCG) ) THEN
        I_VX = 1
        CALL ITL_A( PVOX,P_TPX,P_CRX,TX,I_PX,I_SX,I_TX,I_VX )
        CALL PTL_A( PVOX,TX,FUG_TA,FUGOX,DFUGOX,I_PX,I_SX,I_TX,I_VX )
      ENDIF
      TKX = TX+TABS
!
!---  Poynting effect for CH4 solubility  ---
!
      IF( ISLC(25).EQ.1 ) THEN
        HCX = EXP( LOG(HCX) + (PMV_O*PVOX/(RCU*TKX)) )
      ENDIF
      XMLOX = MIN( MAX( (PVOX*FUGOX/HCX),0.D+0 ),1.D+0 )
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
      SUBROUTINE SOL_N( PX,RHOGX,RHOMGX,RHONX,RHOMNX,TX,WTMNX,
     &  XMGAX,XMGOX,XMGWX,XNAX,XNOX,XNWX )
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
!     Phase equilibria.
!
!     Calculate the gas and NAPL molar volume fractions; component
!     gas and NAPL mole fractions.
!
!     Component index
!       1 - CO2 (carbon dioxide)
!       2 - CH4 (methane)
!       3 - H2O (water)
!     Inputs
!       PX - pressure, Pa
!       TX - temperature, C
!       RHOGX - gas density, kg/m^3 gas
!       RHOMGX - gas molar density, kmol/m^3 gas
!       XMGCX(IGCX) - gas mole fraction of component IGCX
!     Outputs
!       RHONX - NAPL density, kg/m^3 NAPL
!       RHOMNX - NAPL molar density, kmol/m^3 NAPL
!       XMNCX(IGCX) - NAPL mole fraction of component IGCX
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 28 December 2009
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
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
      SAVE IALLOC,AJ,BJ,GX,DXMNCY,XMGCX,XMNCX,XMNCY,XNCX
      SAVE AIX,FUGGX,FUGNX,IJ
      DATA IALLOC /0/
      REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: AJ,GX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: BJ,DXMNCY
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: XMGCX,XMNCX,XMNCY
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: XNCX
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: AIX,FUGGX,FUGNX
      INTEGER, DIMENSION(:), ALLOCATABLE :: IJ
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SOL_N'
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
        ALLOCATE( DXMNCY(1:2),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: DXMNCY'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( FUGGX(1:LNGC),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: FUGGX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( FUGNX(1:LNGC),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: FUGNX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( AIX(1:LNGC),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: AIX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( XNCX(1:LNGC),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: XNCX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( XMGCX(1:LNGC),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: XMGCX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( XMNCX(1:LNGC),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: XMNCX'
          CALL WRMSGS( INDX )
        ENDIF
        ALLOCATE( XMNCY(1:LNGC),STAT=ISTAT )
        IF( ISTAT.NE.0 ) THEN
          INDX = 3
          CHMSG = 'Allocation Error: XMNCY'
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
!---  Peng-Robinson temperature-dependent constants  ---
!
!     (i) At the given TX and PX, evaluate AIX and BIX
!
      CALL PRC_T( AIX,TX )
!
!---  Calculate AFNX and BNX for pure CO2  ---
!
      XMNCY(1) = 1.D+0
      XMNCY(2) = 0.D+0
      XMNCY(3) = 0.D+0
      CALL PRC_X( ANX,AFNX,AIX,BNX,WTMNX,XMNCY )
      CALL DENS_N( ANX,BNX,PX,RHONPX,RHOMNPX,TX,WTMNX,INDX )
!
!---  Skip calculations for pure CO2  ---
!
      IF( (XMGOX+XMGWX)/EPSL.LT.EPSL ) THEN
        XMNCX(1) = 1.D+0
        XMNCX(2) = 0.D+0
        XMNCX(3) = 0.D+0
        GOTO 252
      ENDIF
!
!---  Peng-Robinson composition-dependent constants for gas phase
!
!     (i) and estimate values of AFGX, BGX, and XMGX  ---
!
      XMGCX(1) = XMGAX
      XMGCX(2) = XMGOX
      XMGCX(3) = XMGWX
      CALL PRC_X( AGX,AFGX,AIX,BGX,WTMGX,XMGCX )
!
!---  Newton-Raphson iteration  ---
!
      NC = 0
      XMNCX(2) = 0.D+0
      XMNCX(3) = 0.D+0
  200 CONTINUE
      NC = NC + 1
      IF( NC.GT.64 ) THEN
        INDX = 12
        IMSG = N_DB
        CHMSG = 'Liquid-CO2 Solubility Convergence Failure: ' // 
     &    'Newton-Raphson Loop @ Node'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Increments to primary unknowns
!     mole fraction of water, CO2, and CH4
!     in liquid CO2  ---
!
      DXMNCY(1) = SIGN(1.D-7,5.D-1-XMNCX(2))
      DXMNCY(2) = SIGN(1.D-7,5.D-1-XMNCX(3))
      DO 240 M = 1,3
        XMNCY(2) = XMNCX(2)
        XMNCY(3) = XMNCX(3)
        IF( M.EQ.2 ) XMNCY(2) = XMNCX(2) + DXMNCY(1)
        IF( M.EQ.3 ) XMNCY(3) = XMNCX(3) + DXMNCY(2)
        XMNCY(1) = 1.D+0 - XMNCY(2) - XMNCY(3)
!
!---    Calculate AFNX and BNX  ---
!
        CALL PRC_X( ANX,AFNX,AIX,BNX,WTMNX,XMNCY )
!
!---    Liquid-CO2 molar density  ---
!
        WTMNX = 0.D+0
        DO 120 IGC = 1,NGC
          WTMNX = WTMNX + XMNCY(IGC)*GCPP(1,IGC)
  120   CONTINUE
        RHOMNX = RHONX/WTMNX
!
!---    Component fugacity
!
        CALL FUG( AGX,AIX,BGX,FUGGX,PX,RHOMGX,TX )
        CALL FUG( ANX,AIX,BNX,FUGNX,PX,RHOMNX,TX )
!
!---    Check functions  ---
!
        GX(1,M) = XMNCY(2) - XMGCX(2)*FUGGX(2)/FUGNX(2)
        GX(2,M) = XMNCY(3) - XMGCX(3)*FUGGX(3)/FUGNX(3)
  240 CONTINUE
!
!---  Load solution vector and Jacobian matrix for
!     two-phase conditions  ---
!
      DO 250 M = 1,2
        AJ(M,1) = (GX(M,2)-GX(M,1))/DXMNCY(1)
        AJ(M,2) = (GX(M,3)-GX(M,1))/DXMNCY(2)
        BJ(M) = -GX(M,1)
  250 CONTINUE
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
      DXMNOX = SIGN( MIN(ABS(BJ(1)),1.D-1),BJ(1) )
      DXMNWX = SIGN( MIN(ABS(BJ(2)),1.D-1),BJ(2) )
      XMNCX(2) = XMNCX(2) + DXMNOX
      XMNCX(3) = XMNCX(3) + DXMNWX
      XMNCX(1) = 1.D+0 - XMNCX(2) - XMNCX(3)
      XMNCX(1) = MIN( MAX( XMNCX(1),0.D+0 ),1.D+0 )
      XMNCX(2) = MIN( MAX( XMNCX(2),0.D+0 ),1.D+0 )
      XMNCX(3) = MIN( MAX( XMNCX(3),0.D+0 ),1.D+0 )
!
!---  Limit CH4 and H2O liquid-CO2 mole fractions  ---
!
      IF( XMNCX(2)/EPSL.LT.EPSL .AND. DXMNOX.LT.0.D+0 ) THEN
        DXMNOX = 0.D+0
        XMNCX(2) = 0.D+0
      ENDIF
      IF( XMNCX(3)/EPSL.LT.EPSL .AND. DXMNWX.LT.0.D+0 ) THEN
        DXMNWX = 0.D+0
        XMNCX(3) = 0.D+0
      ENDIF
!
!---  Convergence check  ---
!
      IF( ABS(DXMNOX).GT.1.D-6 ) GOTO 200
      IF( ABS(DXMNWX).GT.1.D-6 ) GOTO 200
  252 CONTINUE
!
!---  Normalize liquid-CO2 mole fractions  ---
!
      SXMNCX = 0.D+0
      DO 260 IGC = 1,NGC
        SXMNCX = SXMNCX + XMNCX(IGC)
  260 CONTINUE
      DO 262 IGC = 1,NGC
        XMNCX(IGC) = XMNCX(IGC)/SXMNCX
  262 CONTINUE
!
!---  Liquid-CO2 molecular weight  ---
!
      WTMNX = 0.D+0
      DO 264 IGC = 1,NGC
        WTMNX = WTMNX + XMNCX(IGC)*GCPP(1,IGC)
  264 CONTINUE
      CALL DENS_N( ANX,BNX,PX,RHONMX,RHOMNMX,TX,WTMNX,INDX )
      RHONX = RHONX*(RHONMX/RHONPX)
!
!---  Liquid-CO2 mass fraction and mass density  ---
!
      DO 266 IGC = 1,NGC
        XNCX(IGC) = XMNCX(IGC)*GCPP(1,IGC)/WTMNX
  266 CONTINUE
      XNAX = XNCX(1)
      XNOX = XNCX(2)
      XNWX = XNCX(3)
      RHOMNX = RHONX/WTMNX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SOL_N group  ---
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
      SAVE A_C,B_C
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
      SAVE SAX,SBX
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
      SAVE TCROX,PCROX,A_C
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
      SAVE K
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

!!----------------------Subroutine--------------------------------------!
!!
!      SUBROUTINE SPLINE( X,Y,N,Y2 )
!!
!!-------------------------Disclaimer-----------------------------------!
!!
!!     This material was prepared as an account of work sponsored by
!!     an agency of the United States Government. Neither the
!!     United States Government nor the United States Department of
!!     Energy, nor Battelle, nor any of their employees, makes any
!!     warranty, express or implied, or assumes any legal liability or
!!     responsibility for the accuracy, completeness, or usefulness
!!     of any information, apparatus, product, software or process
!!     disclosed, or represents that its use would not infringe
!!     privately owned rights.
!!
!!----------------------Acknowledgement---------------------------------!
!!
!!     This software and its documentation were produced with Government
!!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!!     United Department of Energy. The Government retains a paid-up
!!     non-exclusive, irrevocable worldwide license to reproduce,
!!     prepare derivative works, perform publicly and display publicly
!!     by or for the Government, including the right to distribute to
!!     other Government contractors.
!!
!!---------------------Copyright Notices--------------------------------!
!!
!!            Copyright Battelle Memorial Institute, 1996
!!                    All Rights Reserved.
!!
!!----------------------Description-------------------------------------!
!!
!!     Cubic spline second derivative.
!!
!!     Press, W.H., B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling.
!!     1986.  Numerical Recipes, The Art of Scientific Computing,
!!     Cambridge University Press, Cambridge.  pp. 86-89.
!!
!!----------------------Authors-----------------------------------------!
!!
!!     Written by M.D. White, PNNL, 29 August 2006
!!
!!----------------------Fortran 90 Modules------------------------------!
!!
!      USE GLB_PAR
!
!!----------------------Implicit Double Precision-----------------------!
!!
!      IMPLICIT REAL*8 (A-H,O-Z)
!      IMPLICIT INTEGER (I-N)
!!
!!----------------------Type Declarations-------------------------------!
!!
!      REAL*8 X(N),Y(N),Y2(N),U(N)
!!
!!----------------------Executable Lines--------------------------------!
!!
!      Y2(1) = 0.D+0
!      U(1) = 0.D+0
!      DO 100 I = 2,N-1
!        SIG = (X(I)-X(I-1))/(X(I+1)-X(I-1))
!        P = SIG*Y2(I-1)+2.D+0
!        Y2(I) = (SIG-1.D+0)/P
!        U(I) = (6.D+0*((Y(I+1)-Y(I))/(X(I+1)-X(I)) -
!     &    (Y(I)-Y(I-1))/(X(I)-X(I-1)))/
!     &    (X(I+1)-X(I-1)) - SIG*U(I-1))/P
!  100 CONTINUE
!      QN = 0.D+0
!      UN = 0.D+0
!      Y2(N) = (UN-QN*U(N-1))/(QN*Y2(N-1)+1.D+0)
!      DO 200 K = N-1,1,-1
!        Y2(K) = Y2(K)*Y2(K+1)+U(K)
!  200 CONTINUE
!!
!!---  End of SPLINE group  ---
!!
!      RETURN
!      END

!!----------------------Subroutine--------------------------------------!
!!
!      SUBROUTINE SPLINT( XA,YA,Y1A,Y2A,N,X,Y )
!!
!!-------------------------Disclaimer-----------------------------------!
!!
!!     This material was prepared as an account of work sponsored by
!!     an agency of the United States Government. Neither the
!!     United States Government nor the United States Department of
!!     Energy, nor Battelle, nor any of their employees, makes any
!!     warranty, express or implied, or assumes any legal liability or
!!     responsibility for the accuracy, completeness, or usefulness
!!     of any information, apparatus, product, software or process
!!     disclosed, or represents that its use would not infringe
!!     privately owned rights.
!!
!!----------------------Acknowledgement---------------------------------!
!!
!!     This software and its documentation were produced with Government
!!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!!     United Department of Energy. The Government retains a paid-up
!!     non-exclusive, irrevocable worldwide license to reproduce,
!!     prepare derivative works, perform publicly and display publicly
!!     by or for the Government, including the right to distribute to
!!     other Government contractors.
!!
!!---------------------Copyright Notices--------------------------------!
!!
!!            Copyright Battelle Memorial Institute, 1996
!!                    All Rights Reserved.
!!
!!----------------------Description-------------------------------------!
!!
!!     Cubic spline interpolation.
!!
!!     Press, W.H., B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling.
!!     1986.  Numerical Recipes, The Art of Scientific Computing,
!!     Cambridge University Press, Cambridge.  pp. 86-89.
!!
!!----------------------Authors-----------------------------------------!
!!
!!     Written by M.D. White, PNNL, 29 August 2006
!!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
!
!!----------------------Implicit Double Precision-----------------------!
!!
!      IMPLICIT REAL*8 (A-H,O-Z)
!      IMPLICIT INTEGER (I-N)
!!
!!----------------------Type Declarations-------------------------------!
!!
!      REAL*8 XA(N),YA(N),Y2A(N)
!!
!!----------------------Executable Lines--------------------------------!
!!
!      KLO = 1
!      KHI = N
!   10 CONTINUE
!      IF( KHI-KLO.GT.1 ) THEN
!        K = (KHI+KLO)/2
!        IF( XA(K).GT.X ) THEN
!          KHI = K
!        ELSE
!          KLO = K
!        ENDIF
!        GOTO 10
!      ENDIF
!      H = XA(KHI)-XA(KLO)
!      A = (XA(KHI)-X)/H
!      B = (X-XA(KLO))/H
!      Y = A*YA(KLO)+B*YA(KHI)+
!     &    ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.D+0
!      Y1A = (YA(KHI)-YA(KLO))/H -
!     &  ((3.D+0*(A**2)-1.D+0)/6.D+0)*H*Y2A(KLO) +
!     &  ((3.D+0*(B**2)-1.D+0)/6.D+0)*H*Y2A(KHI)
!!
!!---  End of SPLINT group  ---
!!
!      RETURN
!      END

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
      SAVE LFCX,TOLX
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
      REAL*8 TX_TA(LTX),PX_TA(LPX),TKX_TA(LTX,LPX),
     &  PSX_TA(LTX),TKLSX_TA(LTX),TKGSX_TA(LTX)
      REAL*8 TIX(2),PIX(2,2),TKIX(2,2)
!
!----------------------Data Statements---------------------------------!
!
      SAVE SBX,SCX,SDX,ESPX,RHOCRA,TX_TA,PX_TA,PSX_TA,TKLSX_TA,TKGSX_TA
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
      SAVE SCX
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
      SUBROUTINE THK_G( TX,TKAX,TKOX,TKWX,XMAX,XMOX,XMWX,TKGX )
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
      SUB_LOG(ISUB_LOG) = '/THK_G'
!
!---  Reduced, inverse thermal conductivity  ---
!
      TKX = TX+TABS
      GAMMA_W = 210.D+0*(TCRW*(WTMW**3)/((PCRW*1.D-5)**4))
     &  **(1.D+0/6.D+0)
      GAMMA_A = 210.D+0*(TCRA*(WTMA**3)/((PCRA*1.D-5)**4))
     &  **(1.D+0/6.D+0)
      GAMMA_O = 210.D+0*(TCRO*(WTMO**3)/((PCRO*1.D-5)**4))
     &  **(1.D+0/6.D+0)
      TKR_A = TKX/TCRA
      TKR_O = TKX/TCRO
      TKR_W = TKX/TCRW
      TKR_A = GAMMA_A*(EXP(4.64D-2*TKR_A)-EXP(-2.412D-1*TKR_A))
      TKR_O = GAMMA_O*(EXP(4.64D-2*TKR_O)-EXP(-2.412D-1*TKR_O))
      TKR_W = GAMMA_W*(EXP(4.64D-2*TKR_W)-EXP(-2.412D-1*TKR_W))
!
!---  Zero-density-limit component  ---
!
      IF( TKR_W.GT.EPSL ) THEN
        PHIAWX = ((1+SQRT(TKR_A/TKR_W)*((WTMW/WTMA)**2.5D-1))**2)
     &      /SQRT(8.D+0*(1.D+0 + WTMA/WTMW))
        PHIOWX = ((1+SQRT(TKR_O/TKR_W)*((WTMW/WTMO)**2.5D-1))**2)
     &      /SQRT(8.D+0*(1.D+0 + WTMO/WTMW))
      ELSE
        PHIAWX = 0.D+0
        PHIOWX = 0.D+0
      ENDIF
      IF( TKR_A.GT.EPSL ) THEN
        PHIWAX = ((1+SQRT(TKR_W/TKR_A)*((WTMA/WTMW)**2.5D-1))**2)
     &      /SQRT(8.D+0*(1.D+0 + WTMW/WTMA))
        PHIOAX = ((1+SQRT(TKR_O/TKR_A)*((WTMA/WTMO)**2.5D-1))**2)
     &      /SQRT(8.D+0*(1.D+0 + WTMO/WTMA))
      ELSE
        PHIWAX = 0.D+0
        PHIOAX = 0.D+0
      ENDIF
      IF( TKR_O.GT.EPSL ) THEN
        PHIAOX = ((1+SQRT(TKR_A/TKR_O)*((WTMO/WTMA)**2.5D-1))**2)
     &      /SQRT(8.D+0*(1.D+0 + WTMA/WTMO))
        PHIWOX = ((1+SQRT(TKR_W/TKR_O)*((WTMO/WTMW)**2.5D-1))**2)
     &      /SQRT(8.D+0*(1.D+0 + WTMW/WTMO))
      ELSE
        PHIAOX = 0.D+0
        PHIWOX = 0.D+0
      ENDIF
      CHIWX = XMWX + XMAX*PHIWAX + XMOX*PHIWOX
      CHIAX = XMWX*PHIAWX + XMAX + XMOX*PHIAOX
      CHIOX = XMWX*PHIOWX + XMAX*PHIOAX + XMOX
      TKGX = XMWX*TKWX/CHIWX + XMAX*TKAX/CHIAX + XMOX*TKOX/CHIOX
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
      SAVE A
      DATA A / 7.39519D+0,-2.86936D-2,3.54452D-5 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THK_I'
!      THKIX = 2.13D+0 - 1.32D-2*TX
      TKX = TX+TABS
      THKIX = A(1) +A(2)*TKX +A(3)*(TKX**2)
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THK_I group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THK_O( TX,RHOOX,THKOX,VISOOX )
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
!     Kiselev, S.B., and M. L. Huber.  1998.  "Transport properties
!     of carbon dioxide + ethane and methane + ethane mixtures in
!     the extended critical region."  Fluid Phase Equilibria,
!     142:253-280.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 1 December 2004.
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
!----------------------Type Declarations-------------------------------!
!
      REAL*8 AX(4),BX(6),GX(7)
!
!----------------------Data Statements---------------------------------!
!
      SAVE AX,BX,GX
      DATA AX / 1.925D+1, 5.213D-2, 1.197D-5, -1.132D-8 /
      DATA BX / 21.30511D+0, 11.07070D+0, -5.744289D+0,
     &  0.140431D+0, 1.670839D+0, -0.095669D+0 /
      DATA GX / 4.D+0, 3.00573D+0, 1.87D+3, 1.82529D+0, 2.18D+3,
     &  4.07954D+0, 4.17D+3 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THK_O'
!
!---  Ideal gas heat capacity, J/mol K  ---
!
      TKX = TX+TABS
      RCUX = 1.D-3*RCU
!      CPX = GX(1)
!      DO 10 J = 1,3
!        CPY = GX(2*J+1)/TKX
!        CPX = CPX + GX(2*J)*EXP(CPY)*((CPY/(EXP(CPY)-1.D+0))**2)
!   10 CONTINUE
!      CPX = CPX*RCUX
      CPX = AX(1) + AX(2)*TKX + AX(3)*(TKX**2) + AX(4)*(TKX**3)
      F1X = 1.35558587D+0 - 0.11306676*TCRO/TKX
!
!---  Zero-density-limit component, mW/m K  ---
!
      THKOOX = (F1X*(CPX-(5.D+0*RCUX/2.D+0)) + 1.5D+1*RCUX/4.D+0)
     &  *(VISOOX*1.D+6/WTMO)
!
!---  Molar density, mol/L or kmol/m^3 ---
!
      RMOX = RHOOX/WTMO
      RMCROX = 1.D+3/VCRO
!
!---  Excess-thermal conductivity component, mW/m K  ---
!
      RMRX = RMOX/RMCROX
      THKOEX = RMRX*(BX(1) + BX(2)*(RMRX**2) + 
     &  (BX(3) + (BX(4)*TCRO/TKX))*(RMRX**3) + 
     &  (BX(5) + (BX(6)*TCRO/TKX))*(RMRX**4))
!
!---  Sum of components  ---
!
      THKOX = 1.D-3*(THKOOX + THKOEX)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THK_O group  ---
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
      SAVE SAX,SBX,CBX,SDX,CCX
      SAVE TREF,RHOREF,PREF,TKREF
      DATA TREF / 6.4727D+2 /
      DATA RHOREF / 3.17763D+2 /
      DATA PREF / 2.2115D+7 /
      DATA TKREF / 1.D+0 /
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
        TKX = 0.D+0
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
      SAVE SAX,SBX,ESPX
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
      SAVE SAX
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
      SUBROUTINE VISC_G( VISAX,VISOX,VISWX,XMAX,XMOX,XMWX,VISGX )
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
      SUB_LOG(ISUB_LOG) = '/VISC_G'
!
!---  Zero-density-limit component  ---
!
      PHIAW = ((1.D+0+SQRT(VISAX/(VISWX+SMALL))*
     &      SQRT(SQRT(WTMW/WTMA)))**2)/SQRT(8.D+0*(1.D+0 + WTMA/WTMW))
      PHIWA = ((1.D+0+SQRT(VISWX/(VISAX+SMALL))*
     &      SQRT(SQRT(WTMA/WTMW)))**2)/SQRT(8.D+0*(1.D+0 + WTMW/WTMA))
      PHIOW = ((1.D+0+SQRT(VISOX/(VISWX+SMALL))*
     &      SQRT(SQRT(WTMW/WTMO)))**2)/SQRT(8.D+0*(1.D+0 + WTMO/WTMW))
      PHIWO = ((1.D+0+SQRT(VISWX/(VISOX+SMALL))*
     &      SQRT(SQRT(WTMO/WTMW)))**2)/SQRT(8.D+0*(1.D+0 + WTMW/WTMO))
      PHIAO = ((1.D+0+SQRT(VISAX/(VISOX+SMALL))*
     &      SQRT(SQRT(WTMO/WTMA)))**2)/SQRT(8.D+0*(1.D+0 + WTMA/WTMO))
      PHIOA = ((1.D+0+SQRT(VISOX/(VISAX+SMALL))*
     &      SQRT(SQRT(WTMA/WTMO)))**2)/SQRT(8.D+0*(1.D+0 + WTMO/WTMA))
      CHIW = XMWX + XMAX*PHIWA + XMOX*PHIWO
      CHIO = XMWX*PHIOW + XMAX*PHIOA + XMOX
      CHIA = XMWX*PHIAW + XMAX + XMOX*PHIAO
      VISGX = XMWX*VISWX/CHIW + XMOX*VISOX/CHIO + XMAX*VISAX/CHIA
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
      SUBROUTINE VISC_L( XMLAX,XMLOX,VISBX,VISAX,VISOX,VISLX )
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
      VISLX = (1.D+0-XMLAX-XMLOX)*LOG(VISBX)
      IF( VISOX.GT.EPSL ) VISLX = VISLX + XMLOX*LOG(VISOX)
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
      SUBROUTINE VISC_O( TX,RHOOX,VISOOX,VISOEX,VISOX )
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
!     Calculation of CH4 vapor viscosity with the corresponding 
!     states method.  pp 397.
!
!     Reid, R.C., J.M. Prausnitz, and B.E. Poling. 1987.
!     The Properties of Gases and Liquids.
!     McGraw-Hill, New York, New York
!
!     Kiselev, S.B., and M. L. Huber.  1998.  "Transport properties
!     of carbon dioxide + ethane and methane + ethane mixtures in
!     the extended critical region."  Fluid Phase Equilibria,
!     142:253-280.
!
!     Sakonidou, E.P., H.R. van den Berg, C.A. ten Seldam, and
!     J.V. Sengers.  1996.  "The thermal conductivity of methane
!     in the critical region."  J. Chem. Phys. 105(23):10535-10555.
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
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 AX(7),BX(5),CX(6)
!
!----------------------Data Statements---------------------------------!
!
      SAVE AX,BX,CX
      DATA AX / -11.1460D+0, 442.232D+0, 11.9729D+0, -40000.8D+0,
     &  0.01285D+0, 13.2309D+0, 1873.14D+0 /
      DATA BX / 0.226946D+0, -0.46460D+0, -0.63653D-2,
     &  0.10925D+0, -0.32954D-1 /
      DATA CX / 1.16145D+0, 0.14874D+0, 0.52487D+0, 0.77320D+0,
     &  2.16178D+0, 2.43787D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/VISC_O'
!!
!!---  Corresponding states method  ---
!!
!      TR = (TX+TABS)/TCRO
!      DPMR = 5.246D+1*(DPMO**2)*(PCRO/1.D+5)/TCRO
!      IF( DPMR .LT. 2.2D-2 ) THEN
!        FPO = 1.D+0
!      ELSEIF( DPMR .LT. 7.5D-2 ) THEN
!        FPO = 1.D+0 + 3.055D+1*((2.92D-1-ZCRO)**1.72D+0)
!      ELSE
!        FPO = 1.D+0 + 3.055D+1*((2.92D-1-ZCRO)**1.72D+0)*
!     &    ABS(9.6D-1 + 1.D-1*(TR-7.D-1))
!      ENDIF
!      EPSILON = (1.76D-1*((TCRO/((WTMO**3)*((PCRO/1.D+5)**4)))
!     &  **(1.6667D-1)))
!      VISOX = (0.807D+0*(TR**0.618D+0) - 0.357D+0*EXP(-0.449D+0*TR)
!     &  + 0.340D+0*EXP(-4.058D+0*TR) + 0.018D+0)*FPO/EPSILON
!      VISOX = VISOX*1.D-7
!
!---  Zero-density-limit component  ---
!
      TKX = TX+TABS
      TRX = TKX/148.6D+0
      OMEGA = CX(1)/(TRX**CX(2)) + CX(3)*EXP(-CX(4)*TRX)
     &  + CX(5)*EXP(-CX(6)*TRX)
      SIGMA = 0.3758D+0
      VISOOX = 26.69167D-9*SQRT(WTMO*TKX)/OMEGA/(SIGMA**2)
!
!---  Molar density, mol/L or kmol/m^3 ---
!
      RMOX = RHOOX/WTMO
      RMCROX = 1.D+3/VCRO
!
!---  Excess-viscosity component  ---
!
      VISOEX = 1.D-7*EXP(AX(1) + AX(2)/TKX)*
     &  (EXP((AX(3) + AX(4)/SQRT(TKX**3))*(RMOX**1.D-1) +
     &  ((RMOX/RMCROX)-1.D+0)*SQRT(RMOX)*
     &  (AX(5) + AX(6)/TKX + AX(7)/(TKX**2))) - 1.D+0)
!
!---  Sum of components  ---
!
      VISOX = VISOOX + VISOEX
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of VISC_O group  ---
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
      SAVE CHX,TREF,RHOREF,PREF,VISREF
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
      SUBROUTINE VPL_BX( TX,PLX,PVAX,PVOX,PVWX,PSWX,RHOBX,XLSX )
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
      PCX = MAX( BGL*(PVAX+PVOX+PVWX-PLX-PATM),0.D+0 )
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
      SAVE A,SA,B
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WEBB_BC'
      IZN = IZ(N)
!
!---  Find the saturation matching point  ---
!
      SRX = SCHR(4,IZN)
      PSIX = SCHR(1,IZN)
      CLX = MAX( SCHR(3,IZN),SMALL )
!
!---  Use the saturation at 0.4 LOG10(HDOD) as
!     the initial guess  ---
!
      HDX = 1.D+1**(4.D-1*LOG10(HDOD))
      SMPX = (PSIX/HDX)**CLX
      SMPX = SMPX*(1.D+0-SRX) + SRX
!
!---  Newton-Raphson iteration for the saturation
!     matching point  ---
!
      NC = 0
  100 CONTINUE
      NC = NC + 1
      SEMPX = (SMPX-SRX)/(1.D+0-SRX)
      FX = LOG10(HDOD) - LOG10((PSIX/(SEMPX**(1.D+0/CLX))))
     &  - 1.D+0/(LOG(1.D+1)*CLX*(SMPX-SRX))
      DFX = 1.D+0/(LOG(1.D+1)*CLX*(SMPX-SRX))
     &  + 1.D+0/(LOG(1.D+1)*CLX*((SMPX-SRX)**2))
      DSMPX = -FX/DFX
      SMPX = MAX( SMPX+DSMPX,SRX+1.D-12 )
!
!---  No convergence on saturation matching point  ---
!
      IF( NC.GT.32 ) THEN
        INDX = 7
        IMSG = N
        CHMSG = 'No Convergence on Saturation '
     &    // 'Matching Point @ Node: '
        CALL WRMSGS( INDX )
      ENDIF
      IF( ABS(DSMPX).GT.1.D-9 ) GOTO 100
      SCHR(16,IZN) = SMPX
!
!---  Find the capillary head matching point  ---
!
      SEMPX = (SMPX-SRX)/(1.D+0-SRX)
      SCHR(17,IZN) = PSIX/(SEMPX**(1.D+0/CLX))
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
!---  Use the saturation at 0.4 LOG10(HDOD) as
!     the initial guess  ---
!
      HDX = 1.D+1**(4.D-1*LOG10(HDOD))
      SMPX = (1.D+0/(1.D+0 + (ALPHAX*HDX)**CNX))**CMX
      SMPX = SMPX*(1.D+0-SRX) + SRX
!
!---  Newton-Raphson iteration for the saturation
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
!---  No convergence on saturation matching point  ---
!
      IF( NC.GT.32 ) THEN
        INDX = 7
        IMSG = N
        CHMSG = 'No Convergence on Saturation '
     &    // 'Matching Point @ Node: '
        CALL WRMSGS( INDX )
      ENDIF
      IF( ABS(DSMPX).GT.1.D-9 ) GOTO 100
      SCHR(16,IZN) = SMPX
!
!---  Find the capillary head matching point  ---
!
      SEMPX = (SMPX-SRX)/(1.D+0-SRX)
      SCHR(17,IZN) = (1.D+0/ALPHAX)*
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
